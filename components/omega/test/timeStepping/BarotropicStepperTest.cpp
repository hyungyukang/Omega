//===-- Test driver for the OMEGA barotropic steppers ----------*- C++ -*-===/
//
/// \file
/// \brief Test driver for the barotropic steppers of split time stepping
///
/// The barotropic steppers are tested on their own, on the sphere test mesh,
/// with a column flux thickness and a depth-mean specific volume prescribed
/// in the auxiliary state and the equation of state. The driver checks
///  - the options read for each BtrTimeStepper
///  - that a state at rest (uniform pressure anomaly, no velocity, no
///    forcing) stays at rest
///  - that every stepper conserves the area integral of the pressure anomaly
///    (barotropic mass), and that for the SUNDIALS steppers the returned
///    time-averaged flux reproduces the pressure-anomaly change exactly,
///    B'(W) - B'(0) = W (-div FluxAvg), the consistency the transport
///    velocity correction relies on
///  - the convergence of the pressure anomaly at the end of a window with a
///    gravity-wave bump to an RK4 reference: second order for RK2 and for
///    the implicit Crank-Nicolson stepper (two iterations per step, as in
///    SplitExplicitRK2), third order for RK3, and at least first order for
///    the forward-backward predictor-corrector
///  - that the implicit stepper converges with PCG and BiCGStab, with and
///    without Jacobi preconditioning, to the same solution
/// The Coriolis parameter of the test mesh is replaced by the Earth's.
//
//===-----------------------------------------------------------------------===/

#include "BarotropicStepper.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "SundialsBarotropicImplicitStepper.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <string>
#include <vector>

using namespace OMEGA;

namespace {

constexpr int NVertLayers    = 4;
constexpr Real ColumnDepth   = 1000.;  // column pseudo thickness [m]
constexpr Real BumpAmplitude = 0.01;   // sea-surface bump [m]
constexpr Real BumpRadius    = 1.0e6;  // e-folding radius of the bump [m]
constexpr Real ForcingValue  = 1.0e-6; // barotropic forcing [m/s^2]
constexpr Real TestWindow    = 1920.;  // window of the tests [s]
constexpr Real ReferenceStep = 7.5;    // step of the RK4 reference [s]
constexpr Real RhoGravity    = RhoSw * Gravity;

// Inputs and outputs of one barotropic advance
struct BtrArrays {
   Array1DReal VelStart, PressStart, PressRef, VelRef, Forcing;
   Array1DReal Vel, Press, Flux;
};

BtrArrays makeArrays(const std::string &Tag) {
   const auto *Mesh = HorzMesh::getDefault();
   BtrArrays A;
   A.VelStart   = Array1DReal("VelStart" + Tag, Mesh->NEdgesSize);
   A.PressStart = Array1DReal("PressStart" + Tag, Mesh->NCellsSize);
   A.PressRef   = Array1DReal("PressRef" + Tag, Mesh->NCellsSize);
   A.VelRef     = Array1DReal("VelRef" + Tag, Mesh->NEdgesSize);
   A.Forcing    = Array1DReal("Forcing" + Tag, Mesh->NEdgesSize);
   A.Vel        = Array1DReal("Vel" + Tag, Mesh->NEdgesSize);
   A.Press      = Array1DReal("Press" + Tag, Mesh->NCellsSize);
   A.Flux       = Array1DReal("Flux" + Tag, Mesh->NEdgesSize);
   deepCopy(A.Vel, 0._Real);
   deepCopy(A.Press, 0._Real);
   deepCopy(A.Flux, 0._Real);
   return A;
}

// A gravity-wave bump at rest with uniform forcing, or a state at rest with a
// uniform pressure anomaly and no forcing
void setInitialState(BtrArrays &A, bool Rest) {
   const auto *Mesh = HorzMesh::getDefault();

   Array1DReal PressStart = A.PressStart;
   if (Rest) {
      deepCopy(PressStart, 123._Real);
      deepCopy(A.Forcing, 0._Real);
   } else {
      // Bump centered at 30N, 60E
      const Real Radius = std::sqrt(Mesh->XCellH(0) * Mesh->XCellH(0) +
                                    Mesh->YCellH(0) * Mesh->YCellH(0) +
                                    Mesh->ZCellH(0) * Mesh->ZCellH(0));
      const Real Lat0 = M_PI / 6., Lon0 = M_PI / 3.;
      const Real X0 = Radius * std::cos(Lat0) * std::cos(Lon0);
      const Real Y0 = Radius * std::cos(Lat0) * std::sin(Lon0);
      const Real Z0 = Radius * std::sin(Lat0);
      OMEGA_SCOPE(XCell, Mesh->XCell);
      OMEGA_SCOPE(YCell, Mesh->YCell);
      OMEGA_SCOPE(ZCell, Mesh->ZCell);
      parallelFor(
          {Mesh->NCellsSize}, KOKKOS_LAMBDA(int ICell) {
             const Real DX    = XCell(ICell) - X0;
             const Real DY    = YCell(ICell) - Y0;
             const Real DZ    = ZCell(ICell) - Z0;
             const Real Dist2 = DX * DX + DY * DY + DZ * DZ;
             PressStart(ICell) =
                 RhoGravity * BumpAmplitude *
                 Kokkos::exp(-Dist2 / (BumpRadius * BumpRadius));
          });
      deepCopy(A.Forcing, ForcingValue);
   }
   deepCopy(A.PressRef, A.PressStart);
   deepCopy(A.VelStart, 0._Real);
   deepCopy(A.VelRef, 0._Real);
   Halo::getDefault()->exchangeFullArrayHalo(A.PressStart, OnCell);
   Halo::getDefault()->exchangeFullArrayHalo(A.PressRef, OnCell);
}

void advanceOnce(BarotropicStepper &Stepper, const BtrArrays &A, Real Window,
                 I4 NSubsteps) {
   BarotropicInput In;
   In.VelStart   = A.VelStart;
   In.PressStart = A.PressStart;
   In.PressRef   = A.PressRef;
   In.VelRef     = A.VelRef;
   In.Forcing    = A.Forcing;
   In.Window     = Window;
   In.NSubsteps  = NSubsteps;
   BarotropicOutput Out;
   Out.Vel     = A.Vel;
   Out.Press   = A.Press;
   Out.FluxAvg = A.Flux;
   Stepper.advance(In, Out);
}

// Chains NSteps implicit steps over the window, with two iterations per step
// as in SplitExplicitRK2: the second iteration centers the Coriolis term with
// the first iteration's velocity. The flux-thickness reference is kept at the
// initial pressure anomaly so that the equations match the reference run.
void advanceImplicitChain(BarotropicStepper &Stepper, BtrArrays &A, Real Window,
                          I4 NSteps) {
   const auto *Mesh = HorzMesh::getDefault();
   Array1DReal VelCur("ChainVel", Mesh->NEdgesSize);
   Array1DReal PressCur("ChainPress", Mesh->NCellsSize);
   Array1DReal VelIter("ChainVelIter", Mesh->NEdgesSize);
   Array1DReal PressIter("ChainPressIter", Mesh->NCellsSize);
   Array1DReal FluxIter("ChainFluxIter", Mesh->NEdgesSize);
   deepCopy(VelCur, A.VelStart);
   deepCopy(PressCur, A.PressStart);

   const Real Step = Window / NSteps;
   for (I4 N = 0; N < NSteps; ++N) {
      BarotropicInput In;
      In.VelStart   = VelCur;
      In.PressStart = PressCur;
      In.PressRef   = A.PressRef;
      In.VelRef     = VelCur;
      In.Forcing    = A.Forcing;
      In.Window     = Step;
      In.NSubsteps  = 1;
      BarotropicOutput Out;
      Out.Vel     = VelIter;
      Out.Press   = PressIter;
      Out.FluxAvg = FluxIter;
      Stepper.advance(In, Out);

      In.VelRef   = VelIter;
      Out.Vel     = A.Vel;
      Out.Press   = A.Press;
      Out.FluxAvg = A.Flux;
      Stepper.advance(In, Out);

      deepCopy(VelCur, A.Vel);
      deepCopy(PressCur, A.Press);
   }
}

MPI_Comm testComm() { return MachEnv::getDefault()->getComm(); }

// Area integral of a cell field over the owned cells of all tasks
Real areaIntegral(const Array1DReal &Field) {
   const auto *Mesh = HorzMesh::getDefault();
   auto FieldH      = createHostMirrorCopy(Field);
   Real Local       = 0;
   for (int ICell = 0; ICell < Mesh->NCellsOwned; ++ICell)
      Local += Mesh->AreaCellH(ICell) * FieldH(ICell);
   Real Global = 0;
   MPI_Allreduce(&Local, &Global, 1, MPI_RealKind, MPI_SUM, testComm());
   return Global;
}

// Maximum over owned elements of |A - B|
Real maxAbsDiff(const Array1DReal &A, const Array1DReal &B, int NOwned) {
   auto AH    = createHostMirrorCopy(A);
   auto BH    = createHostMirrorCopy(B);
   Real Local = 0;
   for (int I = 0; I < NOwned; ++I)
      Local = std::max(Local, std::abs(AH(I) - BH(I)));
   Real Global = 0;
   MPI_Allreduce(&Local, &Global, 1, MPI_RealKind, MPI_MAX, testComm());
   return Global;
}

Real maxAbs(const Array1DReal &A, int NOwned) {
   auto AH    = createHostMirrorCopy(A);
   Real Local = 0;
   for (int I = 0; I < NOwned; ++I)
      Local = std::max(Local, std::abs(AH(I)));
   Real Global = 0;
   MPI_Allreduce(&Local, &Global, 1, MPI_RealKind, MPI_MAX, testComm());
   return Global;
}

// B'(0) + W (-div Flux) on owned cells
Array1DReal pressFromFlux(const BtrArrays &A, Real Window) {
   const auto *Mesh   = HorzMesh::getDefault();
   const auto *VCoord = VertCoord::getDefault();
   Array1DReal Result("PressFromFlux", Mesh->NCellsSize);
   Array1DReal PressStart = A.PressStart;
   Array1DReal Flux       = A.Flux;
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(EdgeSignOnCell, Mesh->EdgeSignOnCell);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   parallelFor(
       {Mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          Real Div = 0;
          for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
             const int JEdge = EdgesOnCell(ICell, J);
             if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                continue;
             Div += EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * Flux(JEdge);
          }
          Result(ICell) = PressStart(ICell) + Window * Div / AreaCell(ICell);
       });
   return Result;
}

std::unique_ptr<BarotropicStepper>
makeStepper(const BarotropicStepperConfig &Cfg, const std::string &Name) {
   return BarotropicStepper::create(Cfg, Name, AuxiliaryState::getDefault(),
                                    HorzMesh::getDefault(), Halo::getDefault(),
                                    VertCoord::getDefault());
}

BarotropicStepperConfig explicitConfig(const std::string &Method) {
   BarotropicStepperConfig Cfg;
   Cfg.ExplicitMethod = Method;
   Cfg.Method         = Method == "ForwardBackward"
                            ? BarotropicMethod::ForwardBackward
                            : BarotropicMethod::SundialsExplicit;
   return Cfg;
}

BarotropicStepperConfig implicitConfig(BarotropicLinearSolver Solver,
                                       BarotropicPreconditioner Prec) {
   BarotropicStepperConfig Cfg;
   Cfg.Method         = BarotropicMethod::SundialsImplicit;
   Cfg.LinearSolver   = Solver;
   Cfg.Preconditioner = Prec;
   Cfg.RelTolerance   = 1.e-13;
   Cfg.MaxIterations  = 2000;
   return Cfg;
}

//------------------------------------------------------------------------------
int initBarotropicTest(const std::string &MeshFile) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();
   initLogging(DefEnv);
   LOG_INFO("------ Barotropic Stepper Unit Test ------");

   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();
   Clock *ModelClock = TimeStepper::getDefault()->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);
   Halo::init();
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   HorzMesh::init(ModelClock);
   VertCoord::init(false, NVertLayers);
   Tracers::init();
   VertAdv::init();
   AuxiliaryState::init();
   Eos::init();

   // The test mesh has no Coriolis parameter; use the Earth's so that the
   // Coriolis terms of the steppers are exercised
   auto *Mesh = HorzMesh::getDefault();
   OMEGA_SCOPE(FEdge, Mesh->FEdge);
   OMEGA_SCOPE(LatEdge, Mesh->LatEdge);
   parallelFor(
       {Mesh->NEdgesSize}, KOKKOS_LAMBDA(int IEdge) {
          FEdge(IEdge) = 2. * 7.292115e-5 * Kokkos::sin(LatEdge(IEdge));
       });

   // Column flux thickness and depth-mean specific volume seen by the
   // barotropic steppers
   auto *AuxState = AuxiliaryState::getDefault();
   deepCopy(AuxState->PseudoThicknessAux.FluxPseudoThickEdge,
            ColumnDepth / NVertLayers);
   deepCopy(Eos::getInstance()->DepthMeanSpecificVolume, 1._Real / RhoSw);
   return 0;
}

//------------------------------------------------------------------------------
int testReadConfig() {
   int Err = 0;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   OmegaConfig->get(TimeIntConfig);
   Config ShareConfig("ModeSplitShare");
   TimeIntConfig.get(ShareConfig);
   if (!ShareConfig.existsGroup("Barotropic")) {
      Config NewGroup("Barotropic");
      NewGroup.add("ExplicitMethod", std::string("ForwardBackward"));
      ShareConfig.add(NewGroup);
   }
   Config BtrConfig("Barotropic");
   ShareConfig.get(BtrConfig);

   auto setOrAdd = [&](const std::string &Key, auto Value) {
      if (BtrConfig.existsVar(Key))
         BtrConfig.set(Key, Value);
      else
         BtrConfig.add(Key, Value);
   };

   auto Cfg = BarotropicStepper::readConfig("Predictor-Corrector");
   if (Cfg.Method != BarotropicMethod::ForwardBackward) {
      ++Err;
      LOG_ERROR("Predictor-Corrector did not select ForwardBackward");
   }

   setOrAdd("ExplicitMethod", std::string("RK2"));
   Cfg = BarotropicStepper::readConfig("Explicit");
   if (Cfg.Method != BarotropicMethod::SundialsExplicit ||
       Cfg.ExplicitMethod != "RK2") {
      ++Err;
      LOG_ERROR("Explicit with RK2 read wrong options");
   }

   setOrAdd("Theta", R8(0.6));
   setOrAdd("LinearSolver", std::string("BiCGStab"));
   setOrAdd("Preconditioner", std::string("None"));
   setOrAdd("RelTolerance", R8(1.e-8));
   setOrAdd("AbsTolerance", R8(1.e-3));
   setOrAdd("MaxIterations", I4(50));
   Cfg = BarotropicStepper::readConfig("Implicit");
   if (Cfg.Method != BarotropicMethod::SundialsImplicit ||
       std::abs(Cfg.Theta - 0.6) > 1.e-12 ||
       Cfg.LinearSolver != BarotropicLinearSolver::BiCGStab ||
       Cfg.Preconditioner != BarotropicPreconditioner::None ||
       std::abs(Cfg.RelTolerance - 1.e-8) > 1.e-20 ||
       std::abs(Cfg.AbsTolerance - 1.e-3) > 1.e-15 || Cfg.MaxIterations != 50) {
      ++Err;
      LOG_ERROR("Implicit read wrong options");
   }

   // Restore the defaults for later steppers
   setOrAdd("ExplicitMethod", std::string("ForwardBackward"));
   setOrAdd("Theta", R8(0.5));
   setOrAdd("LinearSolver", std::string("PCG"));
   setOrAdd("Preconditioner", std::string("Jacobi"));
   setOrAdd("RelTolerance", R8(1.e-10));
   setOrAdd("AbsTolerance", R8(0.));
   setOrAdd("MaxIterations", I4(1000));

   if (Err == 0)
      LOG_INFO("BarotropicStepperTest: reading options PASS");
   return Err;
}

//------------------------------------------------------------------------------
// A state at rest stays at rest
int testRest(BarotropicStepper &Stepper, const std::string &Name,
             bool Implicit) {
   const auto *Mesh = HorzMesh::getDefault();
   BtrArrays A      = makeArrays("Rest" + Name);
   setInitialState(A, true);

   if (Implicit)
      advanceOnce(Stepper, A, TestWindow, 1);
   else
      advanceOnce(Stepper, A, TestWindow, 32);

   const Real VelMax   = maxAbs(A.Vel, Mesh->NEdgesOwned);
   const Real FluxMax  = maxAbs(A.Flux, Mesh->NEdgesOwned);
   const Real PressErr = maxAbsDiff(A.Press, A.PressStart, Mesh->NCellsOwned);
   if (VelMax > 0 || FluxMax > 0 || PressErr > 1.e-12 * 123.) {
      LOG_ERROR("BarotropicStepperTest: {} does not keep a state at rest: "
                "max |ubar| {}, max |flux| {}, max |dB'| {}",
                Name, VelMax, FluxMax, PressErr);
      return 1;
   }
   LOG_INFO("BarotropicStepperTest: {} state at rest PASS", Name);
   return 0;
}

//------------------------------------------------------------------------------
// Mass conservation, and for the SUNDIALS steppers the flux consistency
int testConservation(BarotropicStepper &Stepper, const std::string &Name,
                     bool Implicit, bool CheckFlux) {
   const auto *Mesh = HorzMesh::getDefault();
   BtrArrays A      = makeArrays("Cons" + Name);
   setInitialState(A, false);

   advanceOnce(Stepper, A, TestWindow, Implicit ? 1 : 32);

   int Err            = 0;
   const Real Mass0   = areaIntegral(A.PressStart);
   const Real Mass1   = areaIntegral(A.Press);
   const Real MassErr = std::abs(Mass1 - Mass0) / std::abs(Mass0);
   const Real Change  = maxAbsDiff(A.Press, A.PressStart, Mesh->NCellsOwned);
   if (MassErr > 1.e-12) {
      ++Err;
      LOG_ERROR("BarotropicStepperTest: {} mass error {}", Name, MassErr);
   }
   if (Change <= 0) {
      ++Err;
      LOG_ERROR("BarotropicStepperTest: {} did not move the bump", Name);
   }

   Real FluxErr = 0;
   if (CheckFlux) {
      Array1DReal FromFlux = pressFromFlux(A, TestWindow);
      FluxErr = maxAbsDiff(FromFlux, A.Press, Mesh->NCellsOwned) / Change;
      if (FluxErr > 1.e-10) {
         ++Err;
         LOG_ERROR("BarotropicStepperTest: {} flux does not reproduce the "
                   "pressure change, relative error {}",
                   Name, FluxErr);
      }
   }
   LOG_INFO("BarotropicStepperTest: {} mass error {:.3e}, flux consistency "
            "error {:.3e}, max |dB'| {:.3e} Pa {}",
            Name, MassErr, FluxErr, Change, Err == 0 ? "PASS" : "FAIL");
   return Err;
}

//------------------------------------------------------------------------------
// Convergence of B' at the end of the window to the RK4 reference
int testConvergence(BarotropicStepper &Stepper, const std::string &Name,
                    bool Implicit, const std::vector<Real> &Steps,
                    const Array1DReal &Reference, Real ExpectedOrder, Real Tol,
                    bool AtLeast) {
   const auto *Mesh = HorzMesh::getDefault();
   std::vector<Real> Errors;
   for (Real Step : Steps) {
      BtrArrays A = makeArrays("Conv" + Name);
      setInitialState(A, false);
      const I4 N = static_cast<I4>(std::lround(TestWindow / Step));
      if (Implicit)
         advanceImplicitChain(Stepper, A, TestWindow, N);
      else
         advanceOnce(Stepper, A, TestWindow, N);
      Errors.push_back(maxAbsDiff(A.Press, Reference, Mesh->NCellsOwned));
   }
   const Real Rate = std::log2(Errors[Errors.size() - 2] / Errors.back());
   const bool Pass = AtLeast ? Rate > ExpectedOrder - Tol
                             : std::abs(Rate - ExpectedOrder) < Tol;
   std::string ErrStr;
   for (size_t I = 0; I < Steps.size(); ++I)
      ErrStr += fmt::format(" dt={} s: {:.3e}", Steps[I], Errors[I]);
   LOG_INFO("BarotropicStepperTest: {} errors{}; rate {:.3f} (expected {}{}) "
            "{}",
            Name, ErrStr, Rate, AtLeast ? ">= " : "", ExpectedOrder,
            Pass ? "PASS" : "FAIL");
   if (!Pass)
      LOG_ERROR("BarotropicStepperTest: {} converges at rate {:.3f}", Name,
                Rate);
   return Pass ? 0 : 1;
}

//------------------------------------------------------------------------------
// The implicit stepper with PCG and BiCGStab, with and without Jacobi
// preconditioning
int testSolvers() {
   const auto *Mesh = HorzMesh::getDefault();
   int Err          = 0;

   struct Variant {
      std::string Name;
      BarotropicLinearSolver Solver;
      BarotropicPreconditioner Prec;
   };
   const std::vector<Variant> Variants = {
       {"PCG+Jacobi", BarotropicLinearSolver::PCG,
        BarotropicPreconditioner::Jacobi},
       {"PCG", BarotropicLinearSolver::PCG, BarotropicPreconditioner::None},
       {"BiCGStab+Jacobi", BarotropicLinearSolver::BiCGStab,
        BarotropicPreconditioner::Jacobi},
       {"BiCGStab", BarotropicLinearSolver::BiCGStab,
        BarotropicPreconditioner::None}};

   Array1DReal First("SolverFirst", Mesh->NCellsSize);
   for (size_t I = 0; I < Variants.size(); ++I) {
      const auto &V = Variants[I];
      auto Cfg      = implicitConfig(V.Solver, V.Prec);
      auto Stepper  = makeStepper(Cfg, "BtrTest" + V.Name);
      auto *Imp =
          dynamic_cast<SundialsBarotropicImplicitStepper *>(Stepper.get());

      BtrArrays A = makeArrays("Solver" + V.Name);
      setInitialState(A, false);
      advanceOnce(*Stepper, A, TestWindow, 1);

      const int Iters = Imp->getLastNumIterations();
      const Real Res  = Imp->getLastResidualNorm();
      const Real Tol  = Imp->getLastTolerance();
      Real Diff       = 0;
      if (I == 0) {
         deepCopy(First, A.Press);
      } else {
         const Real Change =
             maxAbsDiff(A.Press, A.PressStart, Mesh->NCellsOwned);
         Diff = maxAbsDiff(A.Press, First, Mesh->NCellsOwned) / Change;
      }
      const bool Pass = Iters < Cfg.MaxIterations && Res <= Tol && Diff < 1.e-8;
      LOG_INFO("BarotropicStepperTest: implicit {} solve: {} iterations, "
               "residual {:.3e} (tolerance {:.3e}), difference from {} "
               "{:.3e} {}",
               V.Name, Iters, Res, Tol, Variants[0].Name, Diff,
               Pass ? "PASS" : "FAIL");
      if (!Pass) {
         ++Err;
         LOG_ERROR("BarotropicStepperTest: implicit {} solve failed", V.Name);
      }
   }
   return Err;
}

//------------------------------------------------------------------------------
void finalizeBarotropicTest() {
   Tracers::clear();
   TimeStepper::clear();
   Eos::destroyInstance();
   AuxiliaryState::clear();
   VertAdv::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Dimension::clear();
   Field::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//------------------------------------------------------------------------------
int runBarotropicTests() {

   int Err = testReadConfig();

   // RK4 reference for the convergence tests
   const auto *Mesh = HorzMesh::getDefault();
   Array1DReal Reference("BtrReference", Mesh->NCellsSize);
   {
      auto Ref    = makeStepper(explicitConfig("RK4"), "BtrTestRK4Ref");
      BtrArrays A = makeArrays("Reference");
      setInitialState(A, false);
      advanceOnce(*Ref, A, TestWindow,
                  static_cast<I4>(std::lround(TestWindow / ReferenceStep)));
      deepCopy(Reference, A.Press);
   }

   struct ExplicitCase {
      std::string Method;
      Real Order;
      Real Tol;
      bool AtLeast;
      bool CheckFlux;
   };
   const std::vector<ExplicitCase> Cases = {
       {"ForwardBackward", 1., 0.2, true, false},
       {"RK2", 2., 0.3, false, true},
       {"RK3", 3., 0.3, false, true}};

   for (const auto &C : Cases) {
      auto Stepper =
          makeStepper(explicitConfig(C.Method), "BtrTest" + C.Method);
      Err += testRest(*Stepper, C.Method, false);
      Err += testConservation(*Stepper, C.Method, false, C.CheckFlux);
      Err += testConvergence(*Stepper, C.Method, false, {120., 60.}, Reference,
                             C.Order, C.Tol, C.AtLeast);
   }

   {
      auto Cfg     = implicitConfig(BarotropicLinearSolver::PCG,
                                    BarotropicPreconditioner::Jacobi);
      auto Stepper = makeStepper(Cfg, "BtrTestCN");
      Err += testRest(*Stepper, "CrankNicolson", true);
      Err += testConservation(*Stepper, "CrankNicolson", true, true);
      Err += testConvergence(*Stepper, "CrankNicolson", true, {240., 120.},
                             Reference, 2., 0.3, false);
   }

   Err += testSolvers();

   return Err;
}

//------------------------------------------------------------------------------
int barotropicStepperTest(const std::string &MeshFile) {

   int Err = initBarotropicTest(MeshFile);
   // The steppers and arrays of the tests are released before finalizing
   Err += runBarotropicTests();

   if (Err == 0) {
      LOG_INFO("------ Barotropic Stepper Unit Test Successful ------");
   } else {
      LOG_ERROR("------ Barotropic Stepper Unit Test Failed: {} errors ------",
                Err);
   }

   finalizeBarotropicTest();
   return Err;
}

} // namespace

int main(int argc, char *argv[]) {
   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += barotropicStepperTest("OmegaSphereMesh.nc");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   return RetVal;
}
//===-----------------------------------------------------------------------===/

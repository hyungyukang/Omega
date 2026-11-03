//===-- Test driver for OMEGA time steppers -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA time steppers
///
/// This driver tests the OMEGA time stepping module. It checks that every
/// implemented time stepping scheme converges to an exact solution at
/// the expected theoretical rate, including several methods selected by name
/// for the Sundials stepper, and that the Sundials stepper gives the same
/// result as RungeKutta4 and RungeKutta2 with the same methods. To only
/// measure time convergence, the initial conditions and tendency terms are set
/// up such that for every mesh element an independent ODE is solved.
///
//
//===-----------------------------------------------------------------------===/

#include "TimeStepper.h"
#include "../ocn/OceanTestCommon.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "VertMix.h"
#include "mpi.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>
#include <vector>

using namespace OMEGA;

// Only one vertical layer is needed
constexpr int NVertLayers = 1;

// Baroclinic time step at the coarsest refinement level
constexpr Real BaseTimeStepSeconds = 0.2;
constexpr int BtrSubcyclesPerStep  = 10;

// Custom tendency for normal velocity
// du/dt = -coeff * u
struct DecayVelocityTendency {
   Real Coeff = 0.5;

   // exact solution assumes that this is the only tendency active
   // the solution is exponential decay
   Real exactSolution(Real Time) { return std::exp(-Coeff * Time); }

   void operator()(Array2DReal NormalVelTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) const {

      auto *Mesh                = HorzMesh::getDefault();
      auto NVertLayers          = NormalVelTend.extent_int(1);
      Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);

      OMEGA_SCOPE(LocCoeff, Coeff);

      parallelFor(
          {Mesh->NEdgesAll, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int K) {
             NormalVelTend(IEdge, K) -= LocCoeff * NormalVelEdge(IEdge, K);
          });
   }
};

int initState() {
   int Err = 0;

   auto *Mesh              = HorzMesh::getDefault();
   auto *State             = OceanState::get("TestState");
   auto *VCoord            = VertCoord::getDefault();
   Array3DReal TracerArray = Tracers::getAll(0);

   Array2DReal PseudoThickCell = State->getPseudoThickness(0);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(0);

   // Initially set thickness and velocity and tracers to 1
   deepCopy(PseudoThickCell, 1);
   deepCopy(NormalVelEdge, 1);
   deepCopy(TracerArray, 1);

   // Split-explicit time steppers need to have valid values for these
   deepCopy(VCoord->SurfacePressure, 0._Real);
   deepCopy(VCoord->BottomGeomDepth, 500._Real);

   return Err;
}

int createExactSolution(Real TimeEnd) {
   int Err = 0;

   auto *DefHalo           = Halo::getDefault();
   auto *DefMesh           = HorzMesh::getDefault();
   Array3DReal TracerArray = Tracers::getAll(0);

   auto *ExactState =
       OceanState::create("Exact", DefMesh, DefHalo, NVertLayers, 1);

   Array2DReal PseudoThickCell = ExactState->getPseudoThickness(0);
   Array2DReal NormalVelEdge   = ExactState->getNormalVelocity(0);

   // There are no thickness tendencies in this test, so exact thickness ==
   // initial thickness
   deepCopy(PseudoThickCell, 1);
   // Normal velocity decays exponentially
   deepCopy(NormalVelEdge, DecayVelocityTendency{}.exactSolution(TimeEnd));
   // No tracer tendenciesk, final tracers == initial tracers
   deepCopy(TracerArray, 1);

   return Err;
}

ErrorMeasures computeErrors() {
   const auto *DefMesh = HorzMesh::getDefault();

   const auto *State      = OceanState::get("TestState");
   const auto *ExactState = OceanState::get("Exact");

   Array2DReal NormalVelEdge      = State->getNormalVelocity(0);
   Array2DReal ExactNormalVelEdge = ExactState->getNormalVelocity(0);

   // Only velocity errors matters, because thickness remains constant
   ErrorMeasures VelErrors;
   computeErrors(VelErrors, NormalVelEdge, ExactNormalVelEdge, DefMesh, OnEdge);

   return VelErrors;
}

// Returns the Sundials subgroup of the in-memory config, where changes made
// by the tests are seen by Sundials time steppers created afterwards
Config getSundialsConfig() {
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   OmegaConfig->get(TimeIntConfig);
   Config SundialsConfig("Sundials");
   TimeIntConfig.get(SundialsConfig);
   return SundialsConfig;
}

// Selects the method of Sundials time steppers created afterwards
void setSundialsMethod(const std::string &Method) {
   Config SundialsConfig = getSundialsConfig();
   SundialsConfig.set("Method", Method);
}

// Adds Kutta's third-order method as the table used by the Custom method
void addSundialsCustomTable() {
   Config TableConfig("CustomTable");
   TableConfig.add("Stages", 3);
   TableConfig.add("Order", 3);
   TableConfig.add("A", std::vector<R8>{0, 0, 0, 0.5, 0, 0, -1, 2, 0});
   TableConfig.add("B", std::vector<R8>{1. / 6, 2. / 3, 1. / 6});
   TableConfig.add("C", std::vector<R8>{0, 0.5, 1});

   Config SundialsConfig = getSundialsConfig();
   SundialsConfig.add(TableConfig);
}

//------------------------------------------------------------------------------
// The initialization routine for time stepper testing
int initTimeStepperTest(const std::string &mesh) {
   int Err = 0;
   Error Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Default init

   initLogging(DefEnv);
   LOG_INFO("------ Time Stepper Unit Test ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Override BtrTimeStep in the in-memory config so the split-explicit
   // stepper runs BtrSubcyclesPerStep barotropic subcycles per baroclinic
   // step, This is test-only.
   Config *TestOmegaConfig = Config::getOmegaConfig();
   Config TimeIntTestConfig("TimeIntegration");
   TestOmegaConfig->get(TimeIntTestConfig);
   Config ModeSplitShareTestConfig("ModeSplitShare");
   TimeIntTestConfig.get(ModeSplitShareTestConfig);
   ModeSplitShareTestConfig.set(
       "BtrTimeStep",
       std::to_string(BaseTimeStepSeconds / BtrSubcyclesPerStep));

   // A table for the Sundials stepper's Custom method
   addSundialsCustomTable();

   // Note that the default time stepper is not used in subsequent tests
   // but is initialized here because the number of time levels is needed
   // to initialize the Tracers and a model clock is needed for stream IO.
   // If a later timestepper test uses more time levels than the default,
   // this unit test may fail.
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(mesh);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default halo");
   }

   // Initialize IO streams
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   // Read and initialize horizontal mesh
   HorzMesh::init(ModelClock);

   // initialize vertical coordinate, do not read stream and use local
   // NVertLayers value
   VertCoord::init(false, NVertLayers);
   auto *DefVertCoord = VertCoord::getDefault();

   Tracers::init();

   VertAdv::init();
   auto *DefVAdv = VertAdv::getDefault();

   AuxiliaryState::init();
   Eos::init();
   PressureGrad::init();
   VertMix::init();
   Tendencies::init();

   // finish initializing default time stepper
   TimeStepper::init2();

   // Default time stepper never used and time stepper tests require
   // the no calendar option, so we reset calendar here
   Calendar::reset();
   Calendar::init("No Calendar");

   // Non-default init
   // Creating non-default state and auxiliary state to use only one vertical
   // layer

   auto *DefMesh  = HorzMesh::getDefault();
   auto *DefHalo  = Halo::getDefault();
   auto *DefEos   = Eos::getInstance();
   auto *DefPGrad = PressureGrad::getDefault();
   auto *DefVMix  = VertMix::getInstance();

   int NTracers            = Tracers::getNumTracers();
   const int NTimeLevels   = 2;
   const bool UseModeSplit = true;
   auto *TestOceanState    = OceanState::create(
       "TestState", DefMesh, DefHalo, NVertLayers, NTimeLevels, UseModeSplit);

   if (!TestOceanState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test state");
   }

   TimeInterval PosTimeStep(
       1, TimeUnits::Seconds); // Arbitrary positive time step (unused)
   auto *TestAuxState =
       AuxiliaryState::create("TestAuxState", DefMesh, DefHalo, DefVertCoord,
                              DefVAdv, NTracers, PosTimeStep);

   Config *OmegaConfig = Config::getOmegaConfig();
   TestAuxState->readConfigOptions(OmegaConfig);

   if (!TestAuxState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test auxiliary state");
   }

   Config Options;

   // Creating non-default tendencies with custom velocity tendencies
   auto *TestTendencies = Tendencies::create(
       "TestTendencies", DefMesh, DefVertCoord, DefVAdv, DefPGrad, DefEos,
       DefVMix, NTracers, PosTimeStep, &Options,
       Tendencies::CustomTendencyType{}, DecayVelocityTendency{});
   if (!TestTendencies) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test tendencies");
   }

   // Disable all other tendencies
   TestTendencies->PseudoThicknessFluxDiv.Enabled = false;
   TestTendencies->PotentialVortHAdv.Enabled      = false;
   TestTendencies->KEGrad.Enabled                 = false;
   TestTendencies->SSHGrad.Enabled                = false;
   TestTendencies->VelocityDiffusion.Enabled      = false;
   TestTendencies->VelocityHyperDiff.Enabled      = false;
   TestTendencies->TracerHorzAdv.Enabled          = false;
   TestTendencies->TracerDiffusion.Enabled        = false;
   TestTendencies->TracerHyperDiff.Enabled        = false;
   TestTendencies->SfcStressForcing.Enabled       = false;
   TestTendencies->SfcTracerForcing.Enabled       = false;
   TestTendencies->SfcThicknessForcing.Enabled    = false;
   TestTendencies->SurfaceTracerRestoring.Enabled = false;
   TestTendencies->ExplicitBottomDrag.Enabled     = false;
   DefVAdv->ThickVertAdvEnabled                   = false;
   DefVAdv->VelVertAdvEnabled                     = false;
   DefVAdv->TracerVertAdvEnabled                  = false;
   DefVMix->VelVertMixSetup.Enabled               = false;
   DefVMix->TracerVertMixSetup.Enabled            = false;
   TestTendencies->VelVertMixTendencyEnable       = false;
   TestTendencies->TracerVertMixTendencyEnable    = false;

   return Err;
}

// Slightly adjust time step so that it evenly divides TimeEnd and return number
// of steps
int adjustTimeStep(TimeStepper *Stepper, Real TimeEnd) {
   TimeInterval TimeStep = Stepper->getTimeStep();
   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   const int NSteps = std::ceil(TimeEnd / TimeStepSeconds);

   TimeStepSeconds = TimeEnd / NSteps;
   TimeStep.set(TimeStepSeconds, TimeUnits::Seconds);
   Stepper->changeTimeStep(TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

   return NSteps;
}

void timeLoop(TimeInstant TimeStart, Real TimeEnd) {
   auto *Stepper = TimeStepper::get("TestTimeStepper");
   auto *State   = OceanState::get("TestState");

   const int NSteps            = adjustTimeStep(Stepper, TimeEnd);
   const TimeInterval TimeStep = Stepper->getTimeStep();
   TimeInstant CurTime         = Stepper->getStartTime();

   // Time loop
   for (int Step = 0; Step < NSteps; ++Step) {
      Stepper->doStep(State, CurTime);
   }
}

void finalizeTimeStepperTest() {
   Forcing::clear();
   Tracers::clear();
   TimeStepper::clear();
   PressureGrad::clear();
   VertMix::destroyInstance();
   Eos::destroyInstance();
   Tendencies::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   VertAdv::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Dimension::clear();
   Field::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int testTimeStepper(const std::string &Name, TimeStepperType Type,
                    Real ExpectedOrder, Real ATol) {
   int Err = 0;

   // Set pointers to data
   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefVCoord      = VertCoord::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");

   // Set time information
   const TimeInstant TimeStart(0, 0, 0, 0, 0, 0);

   TimeStepperStartType StartType = TimeStepperStartType::StartUp;
   TimeStepperStopType StopType   = TimeStepperStopType::AtTime;
   const Real TimeEnd             = 1;
   TimeInstant TimeEndTI(0, 0, 0, 0, 0, 1);

   TimeInterval TimeStepTI(BaseTimeStepSeconds, TimeUnits::Seconds);

   auto *TestTimeStepper = TimeStepper::create(
       "TestTimeStepper", Type, TestTendencies, TestAuxState, DefMesh,
       DefVCoord, DefHalo, TimeStepTI, StartType, TimeStart, StopType,
       TimeEndTI);

   if (!TestTimeStepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test time stepper {}", Name);
   }

   int NRefinements = 2;
   std::vector<ErrorMeasures> Errors(NRefinements);

   // This creates global exact solution and needs to be done only once
   const static bool CallOnlyOnce = [=]() {
      createExactSolution(TimeEnd);
      return true;
   }();

   R8 TimeStepSeconds = BaseTimeStepSeconds;

   // Convergence loop
   for (int RefLevel = 0; RefLevel < NRefinements; ++RefLevel) {
      TestTimeStepper->changeTimeStep(
          TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

      Err += initState();

      TestTimeStepper->initializeStateFromInput(OceanState::get("TestState"),
                                                false);

      timeLoop(TimeStart, TimeEnd);

      Errors[RefLevel] = computeErrors();

      TimeStepSeconds /= 2;
   }

   std::vector<Real> ConvRates(NRefinements - 1);
   for (int RefLevel = 0; RefLevel < NRefinements - 1; ++RefLevel) {
      ConvRates[RefLevel] =
          std::log2(Errors[RefLevel].LInf / Errors[RefLevel + 1].LInf);
   }

   if (std::abs(ConvRates.back() - ExpectedOrder) > ATol) {
      Err++;
      LOG_ERROR(
          "Wrong convergence rate for time stepper {}, got {:.3f}, expected {}",
          Name, ConvRates.back(), ExpectedOrder);
   }

   TimeStepper::erase("TestTimeStepper");

   return Err;
}

// Final prognostic fields of the test state after a run
struct FinalFields {
   HostArray2DReal Thickness;
   HostArray2DReal Velocity;
   HostArray3DReal TracerArray;
};

// Copies an array to new host memory. createHostMirrorCopy alone returns the
// array itself when it is already on the host, and the state arrays are
// overwritten by later runs.
template <class ArrayType> auto copyToHost(const ArrayType &Arr) {
   auto Mirror = createHostMirrorCopy(Arr);
   auto Copy   = Kokkos::create_mirror(Mirror);
   deepCopy(Copy, Mirror);
   return Copy;
}

// Runs the test problem to the end time at the base time step with a new
// stepper. RungeKutta4Stepper registers a provisional state under its name
// that outlives the stepper, so each stepper created here needs its own name.
FinalFields runToEnd(const std::string &StepperName, TimeStepperType Type) {
   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefVCoord      = VertCoord::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");
   auto *State          = OceanState::get("TestState");

   const TimeInstant TimeStart(0, 0, 0, 0, 0, 0);
   const Real TimeEnd = 1;
   TimeInstant TimeEndTI(0, 0, 0, 0, 0, 1);
   TimeInterval TimeStepTI(BaseTimeStepSeconds, TimeUnits::Seconds);

   auto *Stepper = TimeStepper::create(
       StepperName, Type, TestTendencies, TestAuxState, DefMesh, DefVCoord,
       DefHalo, TimeStepTI, TimeStepperStartType::StartUp, TimeStart,
       TimeStepperStopType::AtTime, TimeEndTI);

   initState();
   Stepper->initializeStateFromInput(State, false);
   const int NSteps    = adjustTimeStep(Stepper, TimeEnd);
   TimeInstant CurTime = Stepper->getStartTime();
   for (int Step = 0; Step < NSteps; ++Step) {
      Stepper->doStep(State, CurTime);
   }

   FinalFields Final{copyToHost(State->getPseudoThickness(0)),
                     copyToHost(State->getNormalVelocity(0)),
                     copyToHost(Tracers::getAll(0))};

   TimeStepper::erase(StepperName);
   return Final;
}

// Checks that the Sundials stepper with Method gives the same result as an
// Omega stepper using the same Runge-Kutta method, to round-off
int testSundialsMatchesNative(const std::string &Method,
                              const std::string &NativeName,
                              TimeStepperType NativeType) {
   int Err = 0;

   const FinalFields Native = runToEnd("Native" + NativeName, NativeType);
   setSundialsMethod(Method);
   const FinalFields Sundials =
       runToEnd("Sundials" + Method, TimeStepperType::Sundials);

   // Compare owned values, which the steps update and which have halo copies
   const auto *Mesh  = HorzMesh::getDefault();
   const int NLayers = Native.Thickness.extent_int(1);
   Real MaxDiff      = 0;
   for (int K = 0; K < NLayers; ++K) {
      for (int ICell = 0; ICell < Mesh->NCellsOwned; ++ICell) {
         MaxDiff = std::max(MaxDiff, std::abs(Native.Thickness(ICell, K) -
                                              Sundials.Thickness(ICell, K)));
         for (int L = 0; L < Native.TracerArray.extent_int(0); ++L) {
            MaxDiff =
                std::max(MaxDiff, std::abs(Native.TracerArray(L, ICell, K) -
                                           Sundials.TracerArray(L, ICell, K)));
         }
      }
      for (int IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
         MaxDiff = std::max(MaxDiff, std::abs(Native.Velocity(IEdge, K) -
                                              Sundials.Velocity(IEdge, K)));
      }
   }

   // The fields are of order one
   const Real Tol = 100 * std::numeric_limits<Real>::epsilon();
   LOG_INFO("TimeStepperTest: Sundials {} vs {}: max difference {:.3e}", Method,
            NativeName, MaxDiff);
   if (MaxDiff > Tol) {
      Err++;
      LOG_ERROR("TimeStepperTest: Sundials {} differs from {} by {:.3e} > "
                "{:.3e}",
                Method, NativeName, MaxDiff, Tol);
   }

   return Err;
}

int testOptionalStopTime(const std::string &Name, TimeStepperType Type) {
   int Err = 0;

   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefVCoord      = VertCoord::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");

   TimeInstant TimeStart(0, 0, 0, 0, 0, 0);
   TimeInterval TimeStep(0.2, TimeUnits::Seconds);

   // 2-phase create without StopTime — used by the coupled driver
   TimeStepperStartType StartType = TimeStepperStartType::StartUp;
   TimeStepperStopType StopType   = TimeStepperStopType::Coupled;
   auto *Stepper = TimeStepper::create("CoupledTestStepper", Type, TimeStep,
                                       StartType, TimeStart, StopType);

   if (!Stepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test stepper StopTime {}",
                Name);
      return Err;
   }

   if (Stepper->getStepCount() != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: {}: getStepCount() should be 0 before any "
                "doStep calls",
                Name);
   }

   // Verify doStep works without an EndAlarm (coupled run loop never uses it)
   Stepper->attachData(TestTendencies, TestAuxState, DefMesh, DefVCoord,
                       DefHalo);
   auto *State = OceanState::get("TestState");
   initState();
   Stepper->initializeStateFromInput(State, false);
   TimeInstant CurTime = TimeStart;
   Stepper->doStep(State, CurTime);

   if (Stepper->getStepCount() != 1) {
      Err++;
      LOG_ERROR("TimeStepperTest: {}: getStepCount() should be 1 after one "
                "doStep call",
                Name);
   }

   // Coupled driver calls doStep repeatedly across separate ocnRun
   // invocations, so StepCount must persist and keep incrementing rather
   // than reset
   Stepper->doStep(State, CurTime);

   if (Stepper->getStepCount() != 2) {
      Err++;
      LOG_ERROR("TimeStepperTest: {}: getStepCount() should be 2 after two "
                "doStep calls",
                Name);
   }

   TimeStepper::erase("CoupledTestStepper");

   return Err;
}

int timeStepperTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initTimeStepperTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Error initializing");
   }

   // Test convergence rate of different time steppers

   Real ExpectedOrder = 4;
   Real ATol          = 0.1;
   Err += testTimeStepper("RungeKutta4", TimeStepperType::RungeKutta4,
                          ExpectedOrder, ATol);

   ExpectedOrder = 1;
   ATol          = 0.1;
   Err += testTimeStepper("ForwardBackward", TimeStepperType::ForwardBackward,
                          ExpectedOrder, ATol);

   ExpectedOrder = 2;
   ATol          = 0.1;
   Err += testTimeStepper("RungeKutta2", TimeStepperType::RungeKutta2,
                          ExpectedOrder, ATol);

   ExpectedOrder = 1;
   ATol          = 0.15;
   Err += testTimeStepper("SE-RK2", TimeStepperType::SplitExplicitRK2,
                          ExpectedOrder, ATol);

   ExpectedOrder = 2;
   ATol          = 0.1;
   Err += testTimeStepper("UnsplitRK2", TimeStepperType::UnsplitRK2,
                          ExpectedOrder, ATol);

   // Methods selected by name for the Sundials stepper: the shorthands,
   // SUNDIALS tables (Bogacki-Shampine reuses its last stage, Zonneveld has
   // five stages) and a custom table
   const std::vector<std::pair<std::string, Real>> SundialsMethods = {
       {"ForwardEuler", 1},
       {"RK2", 2},
       {"RK3", 3},
       {"RK4", 4},
       {"ARKODE_HEUN_EULER_2_1_2", 2},
       {"ARKODE_BOGACKI_SHAMPINE_4_2_3", 3},
       {"ARKODE_ZONNEVELD_5_3_4", 4},
       {"Custom", 3}};
   ATol = 0.1;
   for (const auto &[Method, Order] : SundialsMethods) {
      setSundialsMethod(Method);
      Err += testTimeStepper("Sundials " + Method, TimeStepperType::Sundials,
                             Order, ATol);
   }

   // The same Runge-Kutta methods through SUNDIALS and natively
   Err += testSundialsMatchesNative("RK4", "RungeKutta4",
                                    TimeStepperType::RungeKutta4);
   Err += testSundialsMatchesNative("RK2", "RungeKutta2",
                                    TimeStepperType::RungeKutta2);

   // Verify stepper operates correctly without StopTime/EndAlarm

   Err += testOptionalStopTime("ForwardBackward",
                               TimeStepperType::ForwardBackward);
   Err += testOptionalStopTime("SE-RK2", TimeStepperType::SplitExplicitRK2);
   Err += testOptionalStopTime("UnsplitRK2", TimeStepperType::UnsplitRK2);
   setSundialsMethod("RK4");
   Err += testOptionalStopTime("Sundials", TimeStepperType::Sundials);

   if (Err == 0) {
      LOG_INFO("------ Time Stepper Unit Test Successful ------");
   }

   finalizeTimeStepperTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += timeStepperTest("OmegaSphereMesh.nc");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/

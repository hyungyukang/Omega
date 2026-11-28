//===-- SundialsBarotropicImplicitStepper.cpp - theta barotropic -*- C++
//-*-===//
//
// Implicit (theta / Crank-Nicolson) barotropic stepper: the pseudo thickness
// (barotropic pressure anomaly) is solved implicitly with a SUNDIALS Krylov
// solver and the barotropic velocity is then computed from it explicitly.
//
//===----------------------------------------------------------------------===//

#include "SundialsBarotropicImplicitStepper.h"
#include "AuxiliaryState.h"
#include "Eos.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "Pacer.h"

#include <cmath>
#include <nvector/nvector_mpimanyvector.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>

// SUNDIALS is built once, in the precision of Omega's Real; see
// SundialsStepper.cpp
#if defined(OMEGA_SINGLE_PRECISION) == defined(SUNDIALS_SINGLE_PRECISION)
#define OMEGA_SUNDIALS_MATCHES_REAL
#endif

namespace OMEGA {

namespace {

// Operator and preconditioner callbacks of the SUNDIALS Krylov solver
int helmholtzATimes(void *AData, N_Vector V, N_Vector Z) {
   static_cast<SundialsBarotropicImplicitStepper *>(AData)->applyOperator(V, Z);
   return 0;
}

int jacobiPSolve(void *PData, N_Vector R, N_Vector Z, sunrealtype Tol, int LR) {
   static_cast<SundialsBarotropicImplicitStepper *>(PData)->applyJacobi(R, Z);
   return 0;
}

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
// Device data of an owned-cell vector, an MPI ManyVector of one piece
Real *ownedData(N_Vector Vec) {
   return N_VGetDeviceArrayPointer(N_VGetSubvector_MPIManyVector(Vec, 0));
}
#endif

const char *linearSolverName(BarotropicLinearSolver Solver) {
   return Solver == BarotropicLinearSolver::PCG ? "PCG" : "BiCGStab";
}

} // namespace

//------------------------------------------------------------------------------
SundialsBarotropicImplicitStepper::SundialsBarotropicImplicitStepper(
    const std::string &InName, const BarotropicStepperConfig &Cfg,
    const AuxiliaryState *InAuxState, const HorzMesh *InMesh, Halo *InMeshHalo,
    const VertCoord *InVCoord)
    : BarotropicStepper(InName, "Theta", InAuxState, InMesh, InMeshHalo,
                        InVCoord),
      Config(Cfg) {

#ifndef OMEGA_SUNDIALS_MATCHES_REAL
   ABORT_ERROR("{}: the SUNDIALS barotropic stepper needs SUNDIALS built in "
               "the precision of Omega's Real",
               Name);
#endif

   MethodName = Config.Theta == 0.5_Real ? "CrankNicolson" : "Theta";

   FluxPressEdge   = Array1DReal("BtrImpFluxPress" + Name, Mesh->NEdgesSize);
   OperatorWeight  = Array1DReal("BtrImpOpWeight" + Name, Mesh->NEdgesSize);
   VelTilde        = Array1DReal("BtrImpVelTilde" + Name, Mesh->NEdgesSize);
   MeanSpecVolEdge = Array1DReal("BtrImpSpecVol" + Name, Mesh->NEdgesSize);
   OperatorArg     = Array1DReal("BtrImpOpArg" + Name, Mesh->NCellsSize);
   PressNew        = Array1DReal("BtrImpPressNew" + Name, Mesh->NCellsSize);
   JacobiDiag      = Array1DReal("BtrImpJacobiDiag" + Name, Mesh->NCellsSize);
   deepCopy(FluxPressEdge, 0._Real);
   deepCopy(OperatorWeight, 0._Real);
   deepCopy(VelTilde, 0._Real);
   deepCopy(MeanSpecVolEdge, 0._Real);
   deepCopy(OperatorArg, 0._Real);
   deepCopy(PressNew, 0._Real);
   deepCopy(JacobiDiag, 1._Real);

   LOG_INFO("{}: implicit barotropic stepper, Theta={} ({}), {} with {} "
            "preconditioner, RelTolerance={}, AbsTolerance={}, "
            "MaxIterations={}",
            Name, Config.Theta, MethodName,
            linearSolverName(Config.LinearSolver),
            Config.Preconditioner == BarotropicPreconditioner::Jacobi ? "Jacobi"
                                                                      : "no",
            Config.RelTolerance, Config.AbsTolerance, Config.MaxIterations);
}

//------------------------------------------------------------------------------
SundialsBarotropicImplicitStepper::~SundialsBarotropicImplicitStepper() {
   if (LS)
      SUNLinSolFree(LS);
   if (ResVec)
      N_VDestroy(ResVec);
   if (CorrVec)
      N_VDestroy(CorrVec);
   if (SolVec)
      N_VDestroy(SolVec);
   if (RhsVec)
      N_VDestroy(RhsVec);
   SolPiece.reset();
   RhsPiece.reset();
   if (Ctx)
      SUNContext_Free(&Ctx);
}

//------------------------------------------------------------------------------
void SundialsBarotropicImplicitStepper::createSolver() {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   const MPI_Comm Comm = MeshHalo->getComm();
   if (SUNContext_Create(Comm, &Ctx) != SUN_SUCCESS)
      ABORT_ERROR("{}: SUNContext_Create failed", Name);

   // Owned-cell vectors with their own memory. Omega orders the owned cells
   // first, so these are the first NCellsOwned values of a cell array.
   auto ownedVector = [&](const std::string &Label) {
      SundialsVector::view_type View(Label + Name, Mesh->NCellsOwned);
      auto Vec = std::make_unique<SundialsVector>(View, hostView(View), Ctx);
      enableFusedOps(*Vec);
      return Vec;
   };
   SolPiece = ownedVector("BtrImpSol");
   RhsPiece = ownedVector("BtrImpRhs");

   N_Vector SolPieces[1] = {SolPiece->get()};
   N_Vector RhsPieces[1] = {RhsPiece->get()};
   SolVec                = N_VMake_MPIManyVector(Comm, 1, SolPieces, Ctx);
   RhsVec                = N_VMake_MPIManyVector(Comm, 1, RhsPieces, Ctx);
   OMEGA_REQUIRE(SolVec && RhsVec, "{}: N_VMake_MPIManyVector failed", Name);
   ResVec  = N_VClone(SolVec);
   CorrVec = N_VClone(SolVec);
   OMEGA_REQUIRE(ResVec && CorrVec, "{}: N_VClone failed", Name);

   const bool UsePrec =
       Config.Preconditioner == BarotropicPreconditioner::Jacobi;
   if (Config.LinearSolver == BarotropicLinearSolver::PCG) {
      // Any preconditioning type selects symmetric preconditioning for PCG
      LS = SUNLinSol_PCG(SolVec, UsePrec ? SUN_PREC_LEFT : SUN_PREC_NONE,
                         Config.MaxIterations, Ctx);
   } else {
      // Right preconditioning keeps the convergence test on the true residual
      LS = SUNLinSol_SPBCGS(SolVec, UsePrec ? SUN_PREC_RIGHT : SUN_PREC_NONE,
                            Config.MaxIterations, Ctx);
   }
   OMEGA_REQUIRE(LS, "{}: creating the {} solver failed", Name,
                 linearSolverName(Config.LinearSolver));

   if (SUNLinSolSetATimes(LS, this, helmholtzATimes) != SUN_SUCCESS)
      ABORT_ERROR("{}: SUNLinSolSetATimes failed", Name);
   if (UsePrec && SUNLinSolSetPreconditioner(LS, this, nullptr, jacobiPSolve) !=
                      SUN_SUCCESS)
      ABORT_ERROR("{}: SUNLinSolSetPreconditioner failed", Name);
   if (SUNLinSolInitialize(LS) != SUN_SUCCESS)
      ABORT_ERROR("{}: SUNLinSolInitialize failed", Name);
#endif
}

//------------------------------------------------------------------------------
// Z = A V with
//    (A V)_i = A_i V_i + Theta^2 W^2 sum_e s_ie w_e (V_c1(e) - V_c0(e))
void SundialsBarotropicImplicitStepper::applyOperator(N_Vector V, N_Vector Z) {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   const I4 NCellsOwned = Mesh->NCellsOwned;
   Array1DReal VOwned(ownedData(V), NCellsOwned);
   Array1DReal ZOwned(ownedData(Z), NCellsOwned);
   Array1DReal Arg = OperatorArg;

   parallelFor(
       "btrImpCopyArg", {NCellsOwned},
       KOKKOS_LAMBDA(I4 ICell) { Arg(ICell) = VOwned(ICell); });
   MeshHalo->exchangeFullArrayHalo(Arg, OnCell);

   const Real Coef    = OperatorCoef;
   Array1DReal Weight = OperatorWeight;
   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(EdgeSignOnCell, Mesh->EdgeSignOnCell);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);

   parallelFor(
       "btrImpOperator", {NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          Real Sum = 0._Real;
          for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
             const I4 JEdge = EdgesOnCell(ICell, J);
             const I4 Cell0 = CellsOnEdge(JEdge, 0);
             const I4 Cell1 = CellsOnEdge(JEdge, 1);
             Sum += EdgeSignOnCell(ICell, J) * Weight(JEdge) *
                    (Arg(Cell1) - Arg(Cell0));
          }
          ZOwned(ICell) = AreaCell(ICell) * Arg(ICell) + Coef * Sum;
       });
#endif
}

//------------------------------------------------------------------------------
void SundialsBarotropicImplicitStepper::applyJacobi(N_Vector R, N_Vector Z) {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   const I4 NCellsOwned = Mesh->NCellsOwned;
   Array1DReal ROwned(ownedData(R), NCellsOwned);
   Array1DReal ZOwned(ownedData(Z), NCellsOwned);
   Array1DReal Diag = JacobiDiag;
   parallelFor(
       "btrImpJacobi", {NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          ZOwned(ICell) = ROwned(ICell) / Diag(ICell);
       });
#endif
}

//------------------------------------------------------------------------------
void SundialsBarotropicImplicitStepper::advance(const BarotropicInput &In,
                                                const BarotropicOutput &Out) {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   if (In.Window <= 0)
      LOG_CRITICAL("{}: invalid barotropic window {} s", Name, In.Window);

   Pacer::start("BtrImp:advance", 2);

   if (!LS)
      createSolver();

   const Real W              = In.Window;
   const Real Theta          = Config.Theta;
   constexpr Real RhoGravity = RhoSw * Gravity;
   const bool Upwind         = upwindFluxThickness();
   OperatorCoef              = Theta * Theta * W * W;

   computeColumnFluxThickness();

   Eos *EqState = Eos::getInstance();

   Array1DReal VelStart       = In.VelStart;
   Array1DReal PressStart     = In.PressStart;
   Array1DReal PressRef       = In.PressRef;
   Array1DReal VelRef         = In.VelRef;
   Array1DReal Forcing        = In.Forcing;
   Array1DReal LocColumnThick = ColumnFluxThickEdge;
   Array1DReal LocFluxPress   = FluxPressEdge;
   Array1DReal LocWeight      = OperatorWeight;
   Array1DReal LocVelTilde    = VelTilde;
   Array1DReal LocSpecVolEdge = MeanSpecVolEdge;

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(EdgeSignOnCell, Mesh->EdgeSignOnCell);
   OMEGA_SCOPE(NEdgesOnEdge, Mesh->NEdgesOnEdge);
   OMEGA_SCOPE(EdgesOnEdge, Mesh->EdgesOnEdge);
   OMEGA_SCOPE(WeightsOnEdge, Mesh->WeightsOnEdge);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(FEdge, Mesh->FEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(DepthMeanSpecVol, EqState->DepthMeanSpecificVolume);

   // 1. Edge quantities at the start of the window: flux thickness times
   // rho0 g, operator weights and the explicit part of the new velocity.
   // These read VelRef and PressRef, which the outputs may alias.
   Pacer::start("BtrImp:setup", 3);
   parallelFor(
       "btrImpEdgeSetup", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(I4 IEdge) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          if (MaxLayerEdgeTop(IEdge) < KMin) {
             LocFluxPress(IEdge)   = 0._Real;
             LocWeight(IEdge)      = 0._Real;
             LocVelTilde(IEdge)    = 0._Real;
             LocSpecVolEdge(IEdge) = 0._Real;
             return;
          }
          const I4 Cell0  = CellsOnEdge(IEdge, 0);
          const I4 Cell1  = CellsOnEdge(IEdge, 1);
          const Real Mask = EdgeMask(IEdge, KMin);

          const Real VelOld = VelStart(IEdge);
          const Real SpecVolEdge =
              0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));
          const Real PressEdge =
              RhoGravity * LocColumnThick(IEdge) +
              btrDeltaPressureEdge(PressStart(Cell0) - PressRef(Cell0),
                                   PressStart(Cell1) - PressRef(Cell1), VelOld,
                                   Upwind);

          // Explicit Coriolis term, centered between the start and the
          // provisional velocity
          Real CoriolisTend = 0._Real;
          for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
             const I4 JEdge = EdgesOnEdge(IEdge, J);
             CoriolisTend += WeightsOnEdge(IEdge, J) * 0.5_Real *
                             (VelStart(JEdge) + VelRef(JEdge)) * FEdge(JEdge);
          }
          const Real PressGradOld =
              (PressStart(Cell1) - PressStart(Cell0)) / DcEdge(IEdge);

          LocFluxPress(IEdge)   = PressEdge;
          LocSpecVolEdge(IEdge) = SpecVolEdge;
          LocWeight(IEdge) =
              Mask * DvEdge(IEdge) * PressEdge * SpecVolEdge / DcEdge(IEdge);
          LocVelTilde(IEdge) =
              Mask *
              (VelOld + W * (CoriolisTend -
                             (1._Real - Theta) * SpecVolEdge * PressGradOld +
                             Forcing(IEdge)));
       });

   // 2. Right-hand side, operator diagonal and initial guess B'^n on owned
   // cells
   const I4 NCellsOwned = Mesh->NCellsOwned;
   Array1DReal RhsOwned(ownedData(RhsVec), NCellsOwned);
   Array1DReal SolOwned(ownedData(SolVec), NCellsOwned);
   Array1DReal Diag = JacobiDiag;
   const Real Coef  = OperatorCoef;
   parallelFor(
       "btrImpRhs", {NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          Real FluxDiv   = 0._Real;
          Real WeightSum = 0._Real;
          for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
             const I4 JEdge = EdgesOnCell(ICell, J);
             if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                continue;
             const Real Flux =
                 LocFluxPress(JEdge) * (Theta * LocVelTilde(JEdge) +
                                        (1._Real - Theta) * VelStart(JEdge));
             FluxDiv += EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * Flux;
             WeightSum += LocWeight(JEdge);
          }
          RhsOwned(ICell) = AreaCell(ICell) * PressStart(ICell) + W * FluxDiv;
          Diag(ICell)     = AreaCell(ICell) + Coef * WeightSum;
          SolOwned(ICell) = PressStart(ICell);
       });
   Pacer::stop("BtrImp:setup", 3);

   // 3. Krylov solve, to RelTolerance relative to the right-hand side norm
   // plus AbsTolerance, on the residual 2-norm
   Pacer::start("BtrImp:solve", 3);
   const Real RhsNorm = std::sqrt(N_VDotProd(RhsVec, RhsVec));
   LastTolerance      = Config.RelTolerance * RhsNorm + Config.AbsTolerance;

   // Solve A dx = b - A x0 for the correction to the initial guess x0 = B'^n
   // from a zero guess; the residual of dx is the residual of x0 + dx
   applyOperator(SolVec, ResVec);
   N_VLinearSum(1._Real, RhsVec, -1._Real, ResVec, ResVec);
   N_VConst(0._Real, CorrVec);

   SUNLinSolSetZeroGuess(LS, SUNTRUE);
   if (SUNLinSolSetup(LS, nullptr) != SUN_SUCCESS)
      ABORT_ERROR("{}: SUNLinSolSetup failed", Name);
   const int Flag = SUNLinSolSolve(LS, nullptr, CorrVec, ResVec, LastTolerance);
   LastNumIterations = SUNLinSolNumIters(LS);
   LastResidualNorm  = SUNLinSolResNorm(LS);
   N_VLinearSum(1._Real, SolVec, 1._Real, CorrVec, SolVec);
   Pacer::stop("BtrImp:solve", 3);

   if (Flag < 0)
      ABORT_ERROR("{}: the {} solve of the barotropic pseudo-thickness "
                  "equation failed with flag {}",
                  Name, linearSolverName(Config.LinearSolver), Flag);
   if (Flag != SUN_SUCCESS)
      LOG_WARN("{}: the {} solve did not converge in {} iterations "
               "(residual {}, tolerance {}, flag {})",
               Name, linearSolverName(Config.LinearSolver), LastNumIterations,
               LastResidualNorm, LastTolerance, Flag);
   LOG_DEBUG("{}: {} solve in {} iterations, residual {}, tolerance {}", Name,
             linearSolverName(Config.LinearSolver), LastNumIterations,
             LastResidualNorm, LastTolerance);

   // 4. New pseudo thickness with halo, then the new velocity explicitly from
   // it, and the flux formed from the final velocity
   Array1DReal LocPressNew = PressNew;
   parallelFor(
       "btrImpPressNew", {NCellsOwned},
       KOKKOS_LAMBDA(I4 ICell) { LocPressNew(ICell) = SolOwned(ICell); });
   MeshHalo->exchangeFullArrayHalo(LocPressNew, OnCell);

   Array1DReal OutVel   = Out.Vel;
   Array1DReal OutFlux  = Out.FluxAvg;
   Array1DReal OutPress = Out.Press;
   parallelFor(
       "btrImpVelFlux", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(I4 IEdge) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          if (MaxLayerEdgeTop(IEdge) < KMin) {
             OutVel(IEdge)  = 0._Real;
             OutFlux(IEdge) = 0._Real;
             return;
          }
          const I4 Cell0  = CellsOnEdge(IEdge, 0);
          const I4 Cell1  = CellsOnEdge(IEdge, 1);
          const Real Mask = EdgeMask(IEdge, KMin);
          const Real PressGradNew =
              (LocPressNew(Cell1) - LocPressNew(Cell0)) / DcEdge(IEdge);
          const Real VelNew = LocVelTilde(IEdge) - Mask * Theta * W *
                                                       LocSpecVolEdge(IEdge) *
                                                       PressGradNew;
          OutVel(IEdge) = VelNew;
          OutFlux(IEdge) =
              LocFluxPress(IEdge) *
              (Theta * VelNew + (1._Real - Theta) * VelStart(IEdge));
       });

   // 5. B'^{n+1} from the divergence of that flux, which conserves mass and
   // matches the flux returned for the transport correction exactly
   parallelFor(
       "btrImpPressUpdate", {NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          Real FluxDiv = 0._Real;
          for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
             const I4 JEdge = EdgesOnCell(ICell, J);
             if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                continue;
             FluxDiv +=
                 EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * OutFlux(JEdge);
          }
          OutPress(ICell) = PressStart(ICell) + W * FluxDiv / AreaCell(ICell);
       });

   MeshHalo->exchangeFullArrayHalo(OutVel, OnEdge);
   MeshHalo->exchangeFullArrayHalo(OutFlux, OnEdge);
   MeshHalo->exchangeFullArrayHalo(OutPress, OnCell);

   Pacer::stop("BtrImp:advance", 2);
#else
   ABORT_ERROR("{}: SUNDIALS precision does not match Real", Name);
#endif
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//

//===-- SundialsBarotropicExplicitStepper.cpp - ERK barotropic -*- C++ -*-===//
//
// Explicit Runge-Kutta barotropic stepper built on the SUNDIALS ARKODE
// ERKStep integrator.
//
//===----------------------------------------------------------------------===//

#include "SundialsBarotropicExplicitStepper.h"
#include "AuxiliaryState.h"
#include "Eos.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "Pacer.h"

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep.h>
#include <map>
#include <nvector/nvector_mpimanyvector.h>

// SUNDIALS is built once, in the precision of Omega's Real; see
// SundialsStepper.cpp
#if defined(OMEGA_SINGLE_PRECISION) == defined(SUNDIALS_SINGLE_PRECISION)
#define OMEGA_SUNDIALS_MATCHES_REAL
#endif

namespace OMEGA {

namespace {

// Shorthand method names and the SUNDIALS tables they select, as for the
// Sundials time stepper
const std::map<std::string, std::string> ErkShorthands = {
    {"RK2", "ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2"},
    {"RK3", "ARKODE_SHU_OSHER_3_2_3"},
    {"RK4", "ARKODE_KUTTA_RK4a_4_4"}};

void checkArkode(int Flag, const std::string &Call) {
   if (Flag < 0)
      ABORT_ERROR("SundialsBarotropicExplicitStepper: {} failed with {}", Call,
                  ARKodeGetReturnFlagName(Flag));
}

// Right-hand side called by ARKODE
int barotropicRHS(sunrealtype T, N_Vector Y, N_Vector YDot, void *UserData) {
   static_cast<SundialsBarotropicExplicitStepper *>(UserData)->computeRHS(Y,
                                                                          YDot);
   return 0;
}

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
// Device data of one piece of a state vector
Real *pieceData(N_Vector Vec, int Piece) {
   return N_VGetDeviceArrayPointer(N_VGetSubvector_MPIManyVector(Vec, Piece));
}
#endif

} // namespace

//------------------------------------------------------------------------------
SundialsBarotropicExplicitStepper::SundialsBarotropicExplicitStepper(
    const std::string &InName, const BarotropicStepperConfig &Cfg,
    const AuxiliaryState *InAuxState, const HorzMesh *InMesh, Halo *InMeshHalo,
    const VertCoord *InVCoord)
    : BarotropicStepper(InName, Cfg.ExplicitMethod, InAuxState, InMesh,
                        InMeshHalo, InVCoord) {

#ifndef OMEGA_SUNDIALS_MATCHES_REAL
   ABORT_ERROR("{}: the SUNDIALS barotropic stepper needs SUNDIALS built in "
               "the precision of Omega's Real",
               Name);
#endif

   auto Shorthand = ErkShorthands.find(Cfg.ExplicitMethod);
   MethodName     = Shorthand != ErkShorthands.end() ? Shorthand->second
                                                     : Cfg.ExplicitMethod;
   Table          = ARKodeButcherTable_LoadERKByName(MethodName.c_str());
   if (!Table)
      ABORT_ERROR("{}: unknown barotropic ExplicitMethod {}; use "
                  "ForwardBackward, RK2, RK3, RK4 or the name of a SUNDIALS "
                  "explicit Butcher table",
                  Name, Cfg.ExplicitMethod);

   if (Table->stages <= 2)
      LOG_WARN("{}: barotropic method {} has {} stages; one- and two-stage "
               "explicit Runge-Kutta methods do not include the imaginary "
               "axis in their stability region, so undamped gravity waves "
               "grow slowly. RK3 or RK4 are stable for oscillatory modes.",
               Name, MethodName, Table->stages);

   Vel          = Array1DReal("BtrErkVel" + Name, Mesh->NEdgesSize);
   Press        = Array1DReal("BtrErkPress" + Name, Mesh->NCellsSize);
   FluxIntegral = Array1DReal("BtrErkFluxIntegral" + Name, Mesh->NEdgesSize);
   VelIntegral  = Array1DReal("BtrErkVelIntegral" + Name, Mesh->NEdgesSize);
   PressRef     = Array1DReal("BtrErkPressRef" + Name, Mesh->NCellsSize);
   deepCopy(Vel, 0._Real);
   deepCopy(Press, 0._Real);
   deepCopy(FluxIntegral, 0._Real);
   deepCopy(VelIntegral, 0._Real);
   deepCopy(PressRef, 0._Real);

   LOG_INFO("{}: SUNDIALS explicit barotropic stepper {} with {} stages, "
            "order {}",
            Name, MethodName, Table->stages, Table->q);
}

//------------------------------------------------------------------------------
SundialsBarotropicExplicitStepper::~SundialsBarotropicExplicitStepper() {
   if (ArkMem)
      ARKodeFree(&ArkMem);
   if (StateVec)
      N_VDestroy(StateVec);
   VelVec.reset();
   PressVec.reset();
   FluxVec.reset();
   VelIntVec.reset();
   if (Table)
      ARKodeButcherTable_Free(Table);
   if (Ctx)
      SUNContext_Free(&Ctx);
}

int SundialsBarotropicExplicitStepper::getNumStages() const {
   return Table->stages;
}

int SundialsBarotropicExplicitStepper::getOrder() const { return Table->q; }

//------------------------------------------------------------------------------
void SundialsBarotropicExplicitStepper::createIntegrator() {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   const MPI_Comm Comm = MeshHalo->getComm();
   if (SUNContext_Create(Comm, &Ctx) != SUN_SUCCESS)
      ABORT_ERROR("{}: SUNContext_Create failed", Name);

   VelVec    = wrapArray(Vel, Ctx);
   PressVec  = wrapArray(Press, Ctx);
   FluxVec   = wrapArray(FluxIntegral, Ctx);
   VelIntVec = wrapArray(VelIntegral, Ctx);

   N_Vector Pieces[4] = {VelVec->get(), PressVec->get(), FluxVec->get(),
                         VelIntVec->get()};
   StateVec           = N_VMake_MPIManyVector(Comm, 4, Pieces, Ctx);
   OMEGA_REQUIRE(StateVec, "{}: N_VMake_MPIManyVector failed", Name);

   ArkMem = ERKStepCreate(barotropicRHS, 0.0, StateVec, Ctx);
   OMEGA_REQUIRE(ArkMem, "{}: ERKStepCreate failed", Name);
   checkArkode(ERKStepSetTable(ArkMem, Table), "ERKStepSetTable");
   checkArkode(ARKodeSetUserData(ArkMem, this), "ARKodeSetUserData");
   // One step per call needs no interpolation
   checkArkode(ARKodeSetInterpolantType(ArkMem, ARK_INTERP_NONE),
               "ARKodeSetInterpolantType");
#endif
}

//------------------------------------------------------------------------------
// Right-hand side of the barotropic system and of the flux and velocity
// integrals. Halos of ubar and B' are exchanged in place, which only
// refreshes values the stage already holds on other tasks.
void SundialsBarotropicExplicitStepper::computeRHS(N_Vector Y, N_Vector YDot) {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   Pacer::start("BtrErk:rhs", 3);

   const I4 NEdgesSize  = Mesh->NEdgesSize;
   const I4 NCellsSize  = Mesh->NCellsSize;
   const I4 NEdgesOwned = Mesh->NEdgesOwned;
   const I4 NCellsOwned = Mesh->NCellsOwned;

   Array1DReal YVel(pieceData(Y, 0), NEdgesSize);
   Array1DReal YPress(pieceData(Y, 1), NCellsSize);
   Array1DReal DVel(pieceData(YDot, 0), NEdgesSize);
   Array1DReal DPress(pieceData(YDot, 1), NCellsSize);
   Array1DReal DFlux(pieceData(YDot, 2), NEdgesSize);
   Array1DReal DVelInt(pieceData(YDot, 3), NEdgesSize);

   Pacer::start("BtrErk:halo", 3);
   MeshHalo->exchangeFullArrayHalo(YVel, OnEdge);
   MeshHalo->exchangeFullArrayHalo(YPress, OnCell);
   Pacer::stop("BtrErk:halo", 3);

   Eos *EqState                   = Eos::getInstance();
   const bool UpwindFluxThickness = upwindFluxThickness();
   constexpr Real RhoGravity      = RhoSw * Gravity;

   Array1DReal LocForcing     = Forcing;
   Array1DReal LocPressRef    = PressRef;
   Array1DReal LocColumnThick = ColumnFluxThickEdge;

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

   // Velocity tendency and the integrands of the flux and velocity integrals
   // on owned edges; zero elsewhere so ARKODE leaves those values unchanged
   parallelFor(
       "btrErkEdgeTend", {NEdgesSize}, KOKKOS_LAMBDA(I4 IEdge) {
          DVel(IEdge)    = 0._Real;
          DFlux(IEdge)   = 0._Real;
          DVelInt(IEdge) = 0._Real;
          if (IEdge >= NEdgesOwned)
             return;
          const I4 KMin = MinLayerEdgeBot(IEdge);
          if (MaxLayerEdgeTop(IEdge) < KMin)
             return;

          const I4 Cell0  = CellsOnEdge(IEdge, 0);
          const I4 Cell1  = CellsOnEdge(IEdge, 1);
          const Real Mask = EdgeMask(IEdge, KMin);

          Real CoriolisTend = 0._Real;
          for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
             const I4 JEdge = EdgesOnEdge(IEdge, J);
             CoriolisTend +=
                 WeightsOnEdge(IEdge, J) * YVel(JEdge) * FEdge(JEdge);
          }
          const Real MeanSpecVolEdge =
              0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));
          const Real PressGrad =
              (YPress(Cell1) - YPress(Cell0)) / DcEdge(IEdge);

          DVel(IEdge)    = Mask * (CoriolisTend - MeanSpecVolEdge * PressGrad +
                                LocForcing(IEdge));
          DFlux(IEdge)   = btrEdgeFlux(RhoGravity, LocColumnThick(IEdge),
                                       YPress(Cell0) - LocPressRef(Cell0),
                                       YPress(Cell1) - LocPressRef(Cell1),
                                       YVel(IEdge), UpwindFluxThickness);
          DVelInt(IEdge) = YVel(IEdge);
       });

   // Pressure-anomaly tendency, -div of the flux, on owned cells
   parallelFor(
       "btrErkCellTend", {NCellsSize}, KOKKOS_LAMBDA(I4 ICell) {
          DPress(ICell) = 0._Real;
          if (ICell >= NCellsOwned)
             return;
          Real FluxDivTend = 0._Real;
          for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
             const I4 JEdge = EdgesOnCell(ICell, J);
             if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                continue;
             const I4 Cell0  = CellsOnEdge(JEdge, 0);
             const I4 Cell1  = CellsOnEdge(JEdge, 1);
             const Real Flux = btrEdgeFlux(RhoGravity, LocColumnThick(JEdge),
                                           YPress(Cell0) - LocPressRef(Cell0),
                                           YPress(Cell1) - LocPressRef(Cell1),
                                           YVel(JEdge), UpwindFluxThickness);
             FluxDivTend += EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * Flux;
          }
          DPress(ICell) = FluxDivTend / AreaCell(ICell);
       });

   Pacer::stop("BtrErk:rhs", 3);
#else
   ABORT_ERROR("{}: SUNDIALS precision does not match Real", Name);
#endif
}

//------------------------------------------------------------------------------
void SundialsBarotropicExplicitStepper::advance(const BarotropicInput &In,
                                                const BarotropicOutput &Out) {

   if (In.NSubsteps < 1 || In.Window <= 0)
      LOG_CRITICAL("{}: invalid barotropic window {} s with {} substeps", Name,
                   In.Window, In.NSubsteps);

   Pacer::start("BtrErk:advance", 2);

   if (!ArkMem)
      createIntegrator();

   // Inputs held fixed over the window. The reference is copied because the
   // pressure output may alias it.
   computeColumnFluxThickness();
   deepCopy(PressRef, In.PressRef);
   Forcing = In.Forcing;

   // Start of the window, with zero flux and velocity integrals
   deepCopy(Vel, In.VelStart);
   deepCopy(Press, In.PressStart);
   deepCopy(FluxIntegral, 0._Real);
   deepCopy(VelIntegral, 0._Real);

   const Real StepSize = In.Window / In.NSubsteps;
   checkArkode(ARKodeReset(ArkMem, 0.0, StateVec), "ARKodeReset");
   checkArkode(ARKodeSetFixedStep(ArkMem, StepSize), "ARKodeSetFixedStep");

   sunrealtype TReached = 0.0;
   for (I4 Step = 0; Step < In.NSubsteps; ++Step) {
      checkArkode(ARKodeEvolve(ArkMem, (Step + 1) * StepSize, StateVec,
                               &TReached, ARK_ONE_STEP),
                  "ARKodeEvolve");
   }

   // Time averages over the window and the pressure anomaly at its end
   const Real InvWindow        = 1._Real / In.Window;
   Array1DReal OutVel          = Out.Vel;
   Array1DReal OutFlux         = Out.FluxAvg;
   Array1DReal OutPress        = Out.Press;
   Array1DReal LocFluxIntegral = FluxIntegral;
   Array1DReal LocVelIntegral  = VelIntegral;
   Array1DReal LocPress        = Press;

   parallelFor(
       "btrErkAverage", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
          OutVel(IEdge)  = LocVelIntegral(IEdge) * InvWindow;
          OutFlux(IEdge) = LocFluxIntegral(IEdge) * InvWindow;
       });
   parallelFor(
       "btrErkPressEnd", {Mesh->NCellsOwned},
       KOKKOS_LAMBDA(I4 ICell) { OutPress(ICell) = LocPress(ICell); });

   MeshHalo->exchangeFullArrayHalo(OutVel, OnEdge);
   MeshHalo->exchangeFullArrayHalo(OutFlux, OnEdge);
   MeshHalo->exchangeFullArrayHalo(OutPress, OnCell);

   Pacer::stop("BtrErk:advance", 2);
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//

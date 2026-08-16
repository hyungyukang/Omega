//===-- SplitExplicitRK2Stepper.cpp - split-explicit RK2 ------*- C++ -*-===//
//
// Framework for split-explicit RK2 time stepping. Stage 1 and Stage 3 live
// here, with Stage 2 delegated to the configured barotropic substepper.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitRK2Stepper.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "Pacer.h"
#include "SplitExplicitInit.h"
#include "VertAdv.h"
#include "VertMix.h"

namespace OMEGA {

//------------------------------------------------------------------------------
SplitExplicitRK2Stepper::SplitExplicitRK2Stepper(
    const std::string &InName, TimeStepperType InType,
    const TimeInterval &InTimeStep, const TimeInstant &InStartTime,
    std::optional<TimeInstant> InStopTime)
    : TimeStepper(InName, InType, 2, InTimeStep, InStartTime, InStopTime),
      SEConfig(SplitExplicitInit::readConfigOptions(
          InTimeStep, InType == TimeStepperType::UnsplitRK2)) {}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::finalizeInit() {

   if (!Tend)
      LOG_CRITICAL("Tendency not initialized");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");
   if (!MeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");

   // The velocity tendency is shared with the unsplit time steppers. Tell it
   // to leave the linear Coriolis acceleration out of the vorticity flux, so
   // that doBaroclinicCoriolisIteration can iterate it, and which barotropic
   // weight to use for the surface pressure gradient.
   Tend->setModeSplit(CoriolisTendMode::Separate, SEConfig.SplitFactor);

   SplitExplicitInit::allocateScratch(SEScratch, Mesh, VCoord, Name);
   initBarotropicStepper();
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::initializeStateFromInput(OceanState *State,
                                                       bool ReadRestart) const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   constexpr I4 CurLevel  = 0;
   constexpr I4 NextLevel = 1;

   Array3DReal CurTracerArray = Tracers::getAll(CurLevel);
   AuxState->computeMomVertAux(State, CurTracerArray, CurLevel);

   if (SEConfig.SplitFactor == 0._Real) {
      SplitExplicitInit::computeUnsplitVelocitySplit(State, Mesh, VCoord,
                                                     CurLevel, NextLevel);
   } else if (!ReadRestart) {
      SplitExplicitInit::initializeBarotropicPressure(SEScratch, State, Mesh,
                                                      VCoord, CurLevel);
      SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, CurLevel);
   }

   initializeNextState(State, CurLevel, NextLevel, SEConfig.SplitFactor, false);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::initBarotropicStepper() {

   if (SEConfig.BtrTimeStepper ==
       SplitExplicitBarotropicStepperType::PredictorCorrector) {
      BarotropicUpdate = [this](OceanState *State, I4 CurLevel, I4 NextLevel,
                                const TimeInstant &StageTime,
                                const TimeInterval &StageTimeStep) {
         BarotropicPCStepper.doBarotropicVelocityUpdate(
             State, AuxState, SEScratch, SEConfig, Mesh, MeshHalo, VCoord,
             CurLevel, NextLevel, StageTime, StageTimeStep);
      };
      return;
   }

   LOG_CRITICAL("Invalid split-explicit barotropic time stepper");
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doBaroclinicVelocityUpdate(
    OceanState *State, const Array3DReal &TendencyTracerArray, I4 CurLevel,
    I4 NextLevel, const TimeInstant &StageTime,
    const TimeInterval &StageTimeStep) const {

   Pacer::start("SE-RK2:stage1Bcl", 2);

   // TODO: A placeholder at this moment.
   // prescribeState(State, CurLevel, State, CurLevel, StageTime);

   // Compute baroclinic velocity tendencies and update baroclinic velocity for
   // the first half of the stage time step. This is the same entry point the
   // unsplit time steppers use; the split-specific terms are selected by the
   // mode set in finalizeInit.
   Tend->computeVelocityTendencies(State, AuxState, TendencyTracerArray,
                                   NextLevel, NextLevel, NextLevel, StageTime,
                                   StageTimeStep);

   // Save the baroclinic velocity tendencies before the baroclinic velocity
   // update
   // TODO: this can be optimized in the future.
   deepCopy(SEScratch.BaseVelocityTend, Tend->NormalVelocityTend);

   // Perform baroclinic velocity iteration with Coriolis acceleration.
   doBaroclinicCoriolisIteration(State, SEScratch.BaseVelocityTend, CurLevel,
                                 NextLevel, StageTimeStep);

   Pacer::stop("SE-RK2:stage1Bcl", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doThicknessTracerUpdate(
    OceanState *State, const Array3DReal &CurTracerArray,
    const Array3DReal &NextTracerArray, I4 CurLevel, I4 NextLevel,
    const TimeInstant &StageTime, const TimeInterval &StageTimeStep,
    bool FinalIteration) const {

   Pacer::start("SE-RK2:stage3TrThick", 2);

   // TODO: A placeholder at this moment
   // prescribeState(State, NextLevel, State, CurLevel,
   //               StageTime + 0.5 * StageTimeStep);

   // NormalTransportVelocity is the corrected transport velocity provided by
   // computeTransportVelocity before this stage.
   const Array2DReal NormalTransportVelocity =
       SEScratch.NormalTransportVelocity;
   // Compute thickness auxiliary variables at the new time level
   AuxState->computePseudoThicknessTracerAux(State, NextTracerArray, NextLevel,
                                             NormalTransportVelocity,
                                             StageTimeStep);

   computeVerticalPseudoVelocity(State, NextLevel, NormalTransportVelocity,
                                 StageTimeStep);

   // Compute thickness and tracer tendencies at the new time level
   Tend->computePseudoThicknessTendenciesOnly(
       State, AuxState, NextLevel, NextLevel, NormalTransportVelocity,
       StageTime + 0.5 * StageTimeStep);

   Tend->computeTracerTendenciesOnly(
       State, AuxState, NextTracerArray, NextLevel, NormalTransportVelocity,
       StageTime + 0.5 * StageTimeStep, StageTimeStep);

   if (FinalIteration) {
      // Retain the full-step conservative update on the final iteration.
      updateThicknessByTend(State, NextLevel, State, CurLevel, StageTimeStep);
      updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel,
                          State, CurLevel, StageTimeStep);
   } else {
      // Construct the RK2 midpoint state. Tracer concentration is the average
      // of its old and provisional full-step values, rather than a conservative
      // update over half a time step.
      updateThicknessByTend(State, NextLevel, State, CurLevel,
                            0.5 * StageTimeStep);
      updateTracersToMidpoint(NextTracerArray, CurTracerArray, State, CurLevel,
                              StageTimeStep);
   }

   finalizeTimeStepIterationState(State, CurLevel, NextLevel, FinalIteration);

   Pacer::stop("SE-RK2:stage3TrThick", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doBarotropicVelocityUpdate(
    OceanState *State, I4 CurLevel, I4 NextLevel, const TimeInstant &StageTime,
    const TimeInterval &StageTimeStep) const {

   if (!BarotropicUpdate)
      LOG_CRITICAL("Split-explicit barotropic time stepper not initialized");

   BarotropicUpdate(State, CurLevel, NextLevel, StageTime, StageTimeStep);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::computeTransportVelocity(OceanState *State,
                                                       I4 TimeLevel) const {

   Pacer::start("SE-RK2:computeTransportVelocity", 2);

   Array2DReal NormalVel          = State->getNormalVelocity(TimeLevel);
   Array2DReal NormalTransportVel = SEScratch.NormalTransportVelocity;
   Array2DReal NormalBclVel    = State->getNormalBaroclinicVelocity(TimeLevel);
   Array1DReal NormalBtrVel    = State->getNormalBarotropicVelocity(TimeLevel);
   Array2DReal PseudoThickCell = State->getPseudoThickness(TimeLevel);
   Array1DReal BtrFlux         = SEScratch.BarotropicFlux;
   const Real LocSplitFactor   = SEConfig.SplitFactor;
   constexpr Real RhoGravity   = RhoSw * Gravity;

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);

   parallelForOuter(
       "computeSplitExplicitTransportVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin      = MinLayerEdgeBot(IEdge);
          const I4 KMax      = MaxLayerEdgeTop(IEdge);
          Real VelCorrection = 0._Real;

          if (LocSplitFactor != 0._Real) {

             if (KMax >= KMin) {
                const I4 Cell0 = CellsOnEdge(IEdge, 0);
                const I4 Cell1 = CellsOnEdge(IEdge, 1);

                Real ThickSum     = 0._Real;
                Real TransportSum = 0._Real;

                parallelReduceInner(
                    Team, Range{KMin, KMax},
                    INNER_LAMBDA(I4 K, Real & ThickAccum, Real & FluxAccum) {
                       const Real TransportVel =
                           NormalBtrVel(IEdge) + NormalBclVel(IEdge, K);
                       const Real ThickEdge =
                           0.5_Real * (PseudoThickCell(Cell0, K) +
                                       PseudoThickCell(Cell1, K));
                       ThickAccum += ThickEdge;
                       FluxAccum += ThickEdge * TransportVel;
                    },
                    ThickSum, TransportSum);

                VelCorrection =
                    (BtrFlux(IEdge) / RhoGravity - TransportSum) / ThickSum;
             }

          } // if SplitFactor

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(I4 K) {
                 const Real TotalVel =
                     NormalBtrVel(IEdge) + NormalBclVel(IEdge, K);

                 NormalVel(IEdge, K)          = TotalVel;
                 NormalTransportVel(IEdge, K) = TotalVel + VelCorrection;
              });
       });

   Pacer::stop("SE-RK2:computeTransportVelocity", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doBaroclinicCoriolisIteration(
    OceanState *State, const Array2DReal &BaseVelocityTend, I4 CurLevel,
    I4 NextLevel, const TimeInterval &StageTimeStep) const {

   Pacer::start("SE-RK2:bclCoriolisIter", 2);

   const Array1DReal &FEdge = Mesh->FEdge;
   for (I4 Iter = 0; Iter < SEConfig.NBclCoriolisIteration; ++Iter) {

      const bool FinalCoriolisIter = Iter + 1 == SEConfig.NBclCoriolisIteration;

      if (Iter > 0) {
         deepCopy(Tend->NormalVelocityTend, BaseVelocityTend);
      }

      Array2DReal NormalBclVelEdge =
          State->getNormalBaroclinicVelocity(NextLevel);

      // Compute the baroclinic part of the Coriolis acceleration
      Tend->computeCoriolisAccelerationOnEdge(Tend->NormalVelocityTend,
                                              NormalBclVelEdge, FEdge);

      // Compute the barotropic forcing
      computeBarotropicForcing(State, CurLevel, NextLevel, StageTimeStep);

      // Update the baroclinic velocity by tendency at (n+1/2)
      updateBaroclinicVelocityByTend(State, NextLevel, State, CurLevel,
                                     0.5 * StageTimeStep);

      if (!FinalCoriolisIter) {
         Pacer::start("SE-RK2:haloBclVelIter", 3);
         Array2DReal NormalBclVelIter =
             State->getNormalBaroclinicVelocity(NextLevel);
         MeshHalo->exchangeFullArrayHalo(NormalBclVelIter, OnEdge);
         Pacer::stop("SE-RK2:haloBclVelIter", 3);
      }
   }

   if (SEConfig.SplitFactor != 0._Real) {
      Pacer::start("SE-RK2:haloBtrForcing", 3);
      Array1DReal BtrForcing = SEScratch.BarotropicForcing;
      MeshHalo->exchangeFullArrayHalo(BtrForcing, OnEdge);
      Pacer::stop("SE-RK2:haloBtrForcing", 3);
   }

   Pacer::stop("SE-RK2:bclCoriolisIter", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::computeBarotropicForcing(
    OceanState *State, I4 CurLevel, I4 NextLevel,
    const TimeInterval &StageTimeStep) const {

   const Real LocSplitFactor = SEConfig.SplitFactor;

   // Return if the unsplit time stepper (i.e., BarotropicForcing = 0)
   if (LocSplitFactor == 0._Real)
      return;

   Pacer::start("SE-RK2:barotropicForcing", 2);

   Array2DReal NormalBclVelCur = State->getNormalBaroclinicVelocity(CurLevel);
   Array2DReal PseudoThickCell = State->getPseudoThickness(NextLevel);
   Array2DReal NormalVelTend   = Tend->NormalVelocityTend;
   Array1DReal BtrForcing      = SEScratch.BarotropicForcing;

   R8 DtSecondsR8;
   StageTimeStep.get(DtSecondsR8, TimeUnits::Seconds);
   const Real DtSeconds    = DtSecondsR8;
   const Real InvDtSeconds = 1._Real / DtSeconds;

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);

   parallelForOuter(
       "computeBarotropicForcing", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          if (KMax >= KMin) {

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             Real ThicknessSum = 0._Real;
             parallelReduceInner(
                 Team, Range{KMin, KMax},
                 INNER_LAMBDA(I4 K, Real & ThickAccum) {
                    const Real ThickEdge =
                        0.5_Real *
                        (PseudoThickCell(Cell0, K) + PseudoThickCell(Cell1, K));

                    ThickAccum += ThickEdge;
                 },
                 ThicknessSum);

             Real NormalThicknessFluxSum = 0._Real;
             parallelReduceInner(
                 Team, Range{KMin, KMax},
                 INNER_LAMBDA(I4 K, Real & FluxAccum) {
                    const Real ProvisionalBclVel =
                        NormalBclVelCur(IEdge, K) +
                        DtSeconds * NormalVelTend(IEdge, K);

                    const Real ThickEdge =
                        0.5_Real *
                        (PseudoThickCell(Cell0, K) + PseudoThickCell(Cell1, K));

                    FluxAccum += (ThickEdge / ThicknessSum) * ProvisionalBclVel;
                 },
                 NormalThicknessFluxSum);

             const Real Forcing =
                 LocSplitFactor * NormalThicknessFluxSum * InvDtSeconds;

             Kokkos::single(
                 PerTeam(Team),
                 INNER_LAMBDA() { BtrForcing(IEdge) = Forcing; });

             parallelForInner(
                 Team, Range{KMin, KMax},
                 INNER_LAMBDA(I4 K) { NormalVelTend(IEdge, K) -= Forcing; });

          } else {

             Kokkos::single(
                 PerTeam(Team),
                 INNER_LAMBDA() { BtrForcing(IEdge) = 0._Real; });
          }
       });

   Pacer::stop("SE-RK2:barotropicForcing", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::updateBaroclinicVelocityByTend(
    OceanState *State1, I4 TimeLevel1, OceanState *State2, I4 TimeLevel2,
    TimeInterval Coeff) const {

   Array2DReal NormalBclVel1 = State1->getNormalBaroclinicVelocity(TimeLevel1);
   Array2DReal NormalBclVel2 = State2->getNormalBaroclinicVelocity(TimeLevel2);

   R8 CoeffSeconds;
   Coeff.get(CoeffSeconds, TimeUnits::Seconds);

   OMEGA_SCOPE(NormalVelTend, Tend->NormalVelocityTend);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "updateBclVelByTend", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 NormalBclVel1(IEdge, K) =
                     NormalBclVel2(IEdge, K) +
                     CoeffSeconds * NormalVelTend(IEdge, K);
              });
       });
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::updateTracersToMidpoint(
    const Array3DReal &MidpointTracers, const Array3DReal &CurTracers,
    OceanState *State, I4 CurLevel, const TimeInterval &StageTimeStep) const {

   const Array2DReal CurPseudoThickness = State->getPseudoThickness(CurLevel);

   R8 DtSeconds;
   StageTimeStep.get(DtSeconds, TimeUnits::Seconds);

   OMEGA_SCOPE(PseudoThicknessTend, Tend->PseudoThicknessTend);
   OMEGA_SCOPE(TracerTend, Tend->TracerTend);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   const I4 NTracers = TracerTend.extent_int(0);
   parallelForOuter(
       "updateTracersToMidpoint", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 L, I4 ICell, const TeamMember &Team) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(I4 K) {
                 const Real CurThickness = CurPseudoThickness(ICell, K);
                 const Real EndThickness =
                     CurThickness + DtSeconds * PseudoThicknessTend(ICell, K);
                 const Real EndTracer =
                     (CurTracers(L, ICell, K) * CurThickness +
                      DtSeconds * TracerTend(L, ICell, K)) /
                     EndThickness;

                 MidpointTracers(L, ICell, K) =
                     0.5_Real * (CurTracers(L, ICell, K) + EndTracer);
              });
       });
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::initializeNextState(
    OceanState *State, I4 CurLevel, I4 NextLevel, Real SplitFactor,
    bool ReinitSplitVelocity) const {

   const bool RecomputeSplit = ReinitSplitVelocity && SplitFactor != 0._Real;
   if (RecomputeSplit) {
      SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, CurLevel);
   }

   Array2DReal PseudoThickCur   = State->getPseudoThickness(CurLevel);
   Array2DReal PseudoThickNext  = State->getPseudoThickness(NextLevel);
   Array2DReal NormalVelCur     = State->getNormalVelocity(CurLevel);
   Array2DReal NormalVelNext    = State->getNormalVelocity(NextLevel);
   Array2DReal NormalBclVelCur  = State->getNormalBaroclinicVelocity(CurLevel);
   Array2DReal NormalBclVelNext = State->getNormalBaroclinicVelocity(NextLevel);
   Array1DReal NormalBtrVelCur  = State->getNormalBarotropicVelocity(CurLevel);
   Array1DReal NormalBtrVelNext = State->getNormalBarotropicVelocity(NextLevel);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   if (!RecomputeSplit) {
      parallelForOuter(
          "initializeNormalBaroclinicVelocity", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin = MinLayerEdgeBot(IEdge);
             const int KMax = MaxLayerEdgeTop(IEdge);

             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    NormalBclVelCur(IEdge, K) =
                        NormalVelCur(IEdge, K) - NormalBtrVelCur(IEdge);
                 });
          });
   }

   deepCopy(NormalBclVelNext, NormalBclVelCur);
   deepCopy(PseudoThickNext, PseudoThickCur);
   deepCopy(NormalVelNext, NormalVelCur);

   if (SplitFactor != 0._Real) {
      Array1DReal BtrPressAnomalyCur =
          State->getBarotropicPressureAnomaly(CurLevel);
      Array1DReal BtrPressAnomalyNext =
          State->getBarotropicPressureAnomaly(NextLevel);
      deepCopy(NormalBtrVelNext, NormalBtrVelCur);
      deepCopy(BtrPressAnomalyNext, BtrPressAnomalyCur);
   }
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::reconstructNormalVelocity(
    OceanState *State, I4 CurLevel, I4 NextLevel, bool FinalIteration) const {

   Array2DReal NormalVelNext    = State->getNormalVelocity(NextLevel);
   Array2DReal NormalBclVelCur  = State->getNormalBaroclinicVelocity(CurLevel);
   Array2DReal NormalBclVelNext = State->getNormalBaroclinicVelocity(NextLevel);
   Array1DReal NormalBtrVelNext = State->getNormalBarotropicVelocity(NextLevel);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   if (FinalIteration) {

      parallelForOuter(
          "reconstructFinalNormalVelocity", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin = MinLayerEdgeBot(IEdge);
             const int KMax = MaxLayerEdgeTop(IEdge);

             // Reconstruct NormalBclVel at n+1 if FinalIteration
             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    NormalVelNext(IEdge, K) =
                        2._Real * NormalBclVelNext(IEdge, K) -
                        NormalBclVelCur(IEdge, K) + NormalBtrVelNext(IEdge);
                 });
          });

   } else {

      parallelForOuter(
          "reconstructFinalNormalVelocity", {Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin = MinLayerEdgeBot(IEdge);
             const int KMax = MaxLayerEdgeTop(IEdge);

             // Reconstruct NormalVel at n+0.5
             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    NormalVelNext(IEdge, K) =
                        NormalBclVelNext(IEdge, K) + NormalBtrVelNext(IEdge);
                 });
          });
   } // FinalIteration
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::finalizeTimeStepIterationState(
    OceanState *State, I4 CurLevel, I4 NextLevel, bool FinalIteration) const {

   Pacer::start("SE-RK2:finalizeTimeStepIterationState", 2);

   // Reconstruction of NormalVelocity Next
   reconstructNormalVelocity(State, CurLevel, NextLevel, FinalIteration);

   if (SEConfig.SplitFactor == 0._Real) {

      Pacer::stop("SE-RK2:finalizeTimeStepIterationState", 2);
      return;

   } else {

      // Keep the barotropic pressure anomaly consistent with the most recently
      // updated column pseudo-thickness.
      Array2DReal PseudoThickNext = State->getPseudoThickness(NextLevel);
      VCoord->computeTotalPseudoThickness(PseudoThickNext);

      Array1DReal BtrPressAnomalyNext =
          State->getBarotropicPressureAnomaly(NextLevel);
      constexpr Real RhoGravity = RhoSw * Gravity;

      OMEGA_SCOPE(TotalPseudoThickness, VCoord->TotalPseudoThickness);
      OMEGA_SCOPE(BottomGeomDepth, VCoord->BottomGeomDepth);

      parallelFor(
          "resetBarotropicPressureAnomaly", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             BtrPressAnomalyNext(ICell) =
                 RhoGravity *
                 (TotalPseudoThickness(ICell) - BottomGeomDepth(ICell));
          });

      Pacer::stop("SE-RK2:finalizeTimeStepIterationState", 2);
   }
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::computeVerticalPseudoVelocity(
    OceanState *State, I4 ThickTimeLevel, I4 VelTimeLevel,
    TimeInterval StageTimeStep) const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);
   computeVerticalPseudoVelocity(State, ThickTimeLevel, NormalVelEdge,
                                 StageTimeStep);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::computeVerticalPseudoVelocity(
    OceanState *State, I4 ThickTimeLevel, const Array2DReal &NormalVelEdge,
    TimeInterval StageTimeStep) const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);

   R8 DtSeconds;
   StageTimeStep.get(DtSeconds, TimeUnits::Seconds);

   VertAdv *VertAdvection = VertAdv::getDefault();
   if (!VertAdvection)
      LOG_CRITICAL("Invalid vertical advection");

   const auto &FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   VertAdvection->computeVerticalPseudoVelocity(
       NormalVelEdge, FluxPseudoThickEdge, PseudoThickCell, DtSeconds);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doStep(OceanState *State,
                                     TimeInstant &SimTime) const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   const MPI_Comm Comm = MeshHalo->getComm();

   const int CurLevel  = 0;
   const int NextLevel = 1;
   int NTracers        = Tracers::getNumTracers();

   Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
   Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

   VertMix *VMix = VertMix::getInstance();

   // Initialize NextLevel from CurLevel
   // TODO: This can be optimized in the future.
   initializeNextState(State, CurLevel, NextLevel, SEConfig.SplitFactor,
                       SEConfig.ReinitSplitVelocity);
   deepCopy(NextTracerArray, CurTracerArray);

   const TimeInstant StageTime = SimTime;
   for (I4 TimeStepIteration = 0;
        TimeStepIteration < SEConfig.NTimeStepIteration; ++TimeStepIteration) {

      const bool FinalIteration =
          TimeStepIteration + 1 == SEConfig.NTimeStepIteration;

      // The first iteration evaluates the momentum right-hand side at the
      // n state copied into NextLevel; every later iteration sees the
      // midpoint state left by its predecessor, so time-dependent terms must
      // be sampled at n+1/2 to keep the predictor-corrector second order.
      const TimeInstant VelStageTime =
          TimeStepIteration == 0 ? StageTime : StageTime + 0.5 * TimeStep;

      // Stage 1: Baroclinic velocity advance, with long time step
      doBaroclinicVelocityUpdate(State, NextTracerArray, CurLevel, NextLevel,
                                 VelStageTime, TimeStep);

      Pacer::timingBarrier("SE-RK2:haloStage1Barrier", 3, Comm);
      Pacer::start("SE-RK2:haloStage1", 3);
      Array2DReal NormalBclVelNext =
          State->getNormalBaroclinicVelocity(NextLevel);
      MeshHalo->exchangeFullArrayHalo(NormalBclVelNext, OnEdge);
      Pacer::stop("SE-RK2:haloStage1", 3);

      if (SEConfig.SplitFactor != 0._Real) {
         // Stage 2: Barotropic velocity advance, explicitly subcycled
         doBarotropicVelocityUpdate(State, CurLevel, NextLevel,
                                    StageTime + 0.5 * TimeStep, TimeStep);
      }

      // Compute physical total velocity and the corrected transport velocity
      // used in Stage 3.
      computeTransportVelocity(State, NextLevel);

      // Stage 3: Update thickness, tracers, other diagnostics
      doThicknessTracerUpdate(State, CurTracerArray, NextTracerArray, CurLevel,
                              NextLevel, StageTime, TimeStep, FinalIteration);

      if (TimeStepIteration + 1 < SEConfig.NTimeStepIteration) {
         Pacer::timingBarrier("SE-RK2:haloTimeStepIterationBarrier", 3, Comm);
         Pacer::start("SE-RK2:haloTimeStepIteration", 3);
         State->exchangeHalo(NextLevel);
         MeshHalo->exchangeFullArrayHalo(NextTracerArray, OnCell);
         Pacer::stop("SE-RK2:haloTimeStepIteration", 3);
      }
   }

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   Pacer::timingBarrier("SE-RK2:haloExchBarrier", 3, Comm);
   Pacer::start("SE-RK2:haloExch", 3);
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
   Pacer::stop("SE-RK2:haloExch", 3);

   // Refresh kinetic diagnostics from the completed n+dt velocity before
   // validation and history output.
   const Array2DReal NormalVelCur = State->getNormalVelocity(CurLevel);
   OMEGA_SCOPE(LocKineticAux, AuxState->KineticAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   parallelForOuter(
       "refreshFinalKineticAux", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocKineticAux.computeVarsOnCell(ICell, KChunk, NormalVelCur);
              });
       });

   // Apply implicit vertical mixing
   CurTracerArray = Tracers::getAll(CurLevel);
   if (VMix->VelVertMixSetup.Enabled or VMix->TracerVertMixSetup.Enabled) {
      VMix->VertMixImplicit(State, AuxState, CurTracerArray, NTracers,
                            CurLevel);

      // Re-exchange halos after vertical mixing
      Pacer::timingBarrier("SE-RK2:vMixHaloExchBarrier", 3, Comm);
      Pacer::start("SE-RK2:vMixHaloExch", 3);
      State->exchangeHalo(CurLevel);
      Tracers::exchangeHalo(CurLevel);
      Pacer::stop("SE-RK2:vMixHaloExch", 3);
   }

   validateOceanState(State, AuxState, VertCoord::getDefault(), CurLevel);

   // Advance the clock and update the simulation time
   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
   ++StepCount;
}

} // namespace OMEGA

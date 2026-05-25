//===-- SplitExplicitRK2Stepper.cpp - split-explicit RK2 ------*- C++ -*-===//
//
// Framework for split-explicit RK2 time stepping. Stage 1 and Stage 3 live
// here, with Stage 2 delegated to the configured barotropic substepper.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitRK2Stepper.h"
#include "Logging.h"
#include "Pacer.h"
#include "SplitExplicitInit.h"
#include "VertAdv.h"

namespace OMEGA {

//------------------------------------------------------------------------------
SplitExplicitRK2Stepper::SplitExplicitRK2Stepper(
    const std::string &InName, const TimeInterval &InTimeStep,
    const TimeInstant &InStartTime, std::optional<TimeInstant> InStopTime)
    : TimeStepper(InName, TimeStepperType::SplitExplicitRK2, 2, InTimeStep,
                  InStartTime, InStopTime),
      SEConfig(SplitExplicitInit::readConfigOptions(InTimeStep)) {}

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

   SplitExplicitInit::allocateScratch(SEScratch, Mesh, VCoord, Name);
   initBarotropicStepper();
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::initBarotropicStepper() {

   if (SEConfig.BtrTimeStepper ==
       SplitExplicitBarotropicStepperType::PredictorCorrector) {
      BarotropicStage2 =
          [this](OceanState *State, I4 TimeLevel,
                 const TimeInstant &StageTime,
                 const TimeInterval &StageTimeStep) {
             BarotropicPCStepper.doSplitStage2(
                 State, SEScratch, SEConfig, Mesh, VCoord, TimeLevel,
                 StageTime, StageTimeStep);
          };
      return;
   }

   LOG_CRITICAL("Invalid split-explicit barotropic time stepper");
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doSplitStage1(
    OceanState *State, const Array3DReal &CurTracerArray,
    const Array3DReal &NextTracerArray, I4 CurLevel, I4 NextLevel,
    const TimeInstant &StageTime, const TimeInterval &StageTimeStep) const {

   Pacer::start("SE-RK2:stage1Bcl", 2);

   SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, CurLevel);

   prescribeState(State, CurLevel, State, CurLevel, StageTime);

   Tend->computeBaroclinicVelocityTendencies(
       State, AuxState, CurTracerArray, CurLevel, CurLevel, CurLevel, CurLevel,
       SEConfig.SplitFactor, 0.5 * StageTimeStep);

   deepCopy(SEScratch.BaseVelocityTend, Tend->NormalVelocityTend);

   updateBaroclinicVelocityByTend(State, NextLevel, State, CurLevel,
                                  0.5 * StageTimeStep);

   doBaroclinicCoriolisIteration(State, SEScratch.BaseVelocityTend, CurLevel,
                                 NextLevel, StageTimeStep);

   Pacer::stop("SE-RK2:stage1Bcl", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doSplitStage3(
    OceanState *State, const Array3DReal &CurTracerArray,
    const Array3DReal &NextTracerArray, I4 CurLevel, I4 NextLevel,
    const TimeInstant &StageTime, const TimeInterval &StageTimeStep) const {

   Pacer::start("SE-RK2:stage3TrThick", 2);

   prescribeState(State, NextLevel, State, CurLevel,
                  StageTime + 0.5 * StageTimeStep);

   reconstructNormalVelocity(State, NextLevel);

   AuxState->computeThicknessTracerAux(State, NextTracerArray, NextLevel,
                                       NextLevel);
   computeVerticalVelocity(State, NextLevel, NextLevel, 0.5 * StageTimeStep);
   Tend->computeThicknessTendenciesOnly(
       State, AuxState, NextLevel, NextLevel,
       StageTime + 0.5 * StageTimeStep);
   Tend->computeTracerTendenciesOnly(State, AuxState, NextTracerArray,
                                     NextLevel, NextLevel,
                                     StageTime + 0.5 * StageTimeStep);

   updateThicknessByTend(State, NextLevel, State, CurLevel, StageTimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel,
                       State, CurLevel, StageTimeStep);

   SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, NextLevel);

   Pacer::stop("SE-RK2:stage3TrThick", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doSplitStage2(
    OceanState *State, I4 TimeLevel, const TimeInstant &StageTime,
    const TimeInterval &StageTimeStep) const {

   if (!BarotropicStage2)
      LOG_CRITICAL("Split-explicit barotropic time stepper not initialized");

   BarotropicStage2(State, TimeLevel, StageTime, StageTimeStep);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doBaroclinicCoriolisIteration(
    OceanState *State, const Array2DReal &BaseVelocityTend, I4 CurLevel,
    I4 NextLevel, const TimeInterval &StageTimeStep) const {

   Pacer::start("SE-RK2:bclCoriolisIter", 2);

   const Array1DReal &FEdge = Mesh->FEdge;
   for (I4 Iter = 0; Iter < SEConfig.NBclCoriolisIteration; ++Iter) {
      if (Iter > 0) {
         deepCopy(Tend->NormalVelocityTend, BaseVelocityTend);
      }

      Array2DReal NormalBclVelEdge =
          State->getNormalBaroclinicVelocity(NextLevel);
      Tend->computeCoriolisAccelerationOnEdge(Tend->NormalVelocityTend,
                                              NormalBclVelEdge, FEdge);

      updateBaroclinicVelocityByTend(State, NextLevel, State, CurLevel,
                                     0.5 * StageTimeStep);
   }

   Pacer::stop("SE-RK2:bclCoriolisIter", 2);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::updateBaroclinicVelocityByTend(
    OceanState *State1, I4 TimeLevel1, OceanState *State2, I4 TimeLevel2,
    TimeInterval Coeff) const {

   Array2DReal NormalBclVel1 =
       State1->getNormalBaroclinicVelocity(TimeLevel1);
   Array2DReal NormalBclVel2 =
       State2->getNormalBaroclinicVelocity(TimeLevel2);

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
void SplitExplicitRK2Stepper::reconstructNormalVelocity(OceanState *State,
                                                        I4 TimeLevel) const {

   SplitExplicitInit::combineVelocitySplit(State, Mesh, VCoord, TimeLevel);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::computeVerticalVelocity(
    OceanState *State, I4 ThickTimeLevel, I4 VelTimeLevel,
    TimeInterval StageTimeStep) const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal LayerThickCell = State->getLayerThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge  = State->getNormalVelocity(VelTimeLevel);

   R8 DtSeconds;
   StageTimeStep.get(DtSeconds, TimeUnits::Seconds);

   VertAdv *VertAdvection = VertAdv::getDefault();
   if (!VertAdvection)
      LOG_CRITICAL("Invalid vertical advection");

   const auto &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   VertAdvection->computeVerticalVelocity(
       NormalVelEdge, FluxLayerThickEdge, LayerThickCell, DtSeconds);
}

//------------------------------------------------------------------------------
void SplitExplicitRK2Stepper::doStep(OceanState *State, TimeInstant &SimTime)
    const {

   if (!State)
      LOG_CRITICAL("Invalid State");

   const int CurLevel  = 0;
   const int NextLevel = 1;

   const MPI_Comm Comm = MeshHalo->getComm();
   const TimeInterval TimeStepIterationTimeStep =
       (1._Real / static_cast<Real>(SEConfig.NTimeStepIteration)) * TimeStep;

   TimeInstant StageTime = SimTime;
   for (I4 TimeStepIteration = 0;
        TimeStepIteration < SEConfig.NTimeStepIteration; ++TimeStepIteration) {
      Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
      Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

      doSplitStage1(State, CurTracerArray, NextTracerArray, CurLevel,
                    NextLevel, StageTime, TimeStepIterationTimeStep);

      Pacer::timingBarrier("SE-RK2:haloStage1Barrier", 3, Comm);
      Pacer::start("SE-RK2:haloStage1", 3);
      MeshHalo->exchangeFullArrayHalo(
          State->getNormalBaroclinicVelocity(NextLevel), OnEdge);
      Pacer::stop("SE-RK2:haloStage1", 3);

      if (SEConfig.SplitFactor != 0._Real) {
         doSplitStage2(State, NextLevel,
                       StageTime + 0.5 * TimeStepIterationTimeStep,
                       TimeStepIterationTimeStep);
      }

      doSplitStage3(State, CurTracerArray, NextTracerArray, CurLevel,
                    NextLevel, StageTime, TimeStepIterationTimeStep);

      if (TimeStepIteration + 1 < SEConfig.NTimeStepIteration) {
         Pacer::timingBarrier("SE-RK2:haloTimeStepIterationBarrier", 3, Comm);
         Pacer::start("SE-RK2:haloTimeStepIteration", 3);
         State->exchangeHalo(NextLevel);
         MeshHalo->exchangeFullArrayHalo(NextTracerArray, OnCell);
         Pacer::stop("SE-RK2:haloTimeStepIteration", 3);

         State->updateTimeLevels();
         Tracers::updateTimeLevels();
         StageTime = StageTime + TimeStepIterationTimeStep;
      }
   }

   Pacer::timingBarrier("SE-RK2:haloExchBarrier", 3, Comm);
   Pacer::start("SE-RK2:haloExch", 3);
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
   Pacer::stop("SE-RK2:haloExch", 3);

   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
}

} // namespace OMEGA

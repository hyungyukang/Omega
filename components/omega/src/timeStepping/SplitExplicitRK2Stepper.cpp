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

   SplitExplicitInit::allocateBuffers(SEBuffers, Mesh, Name);
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
       State, AuxState, CurTracerArray, CurLevel, CurLevel,
       0.5 * StageTimeStep);

   const TimeInterval ZeroTimeStep = 0._Real * StageTimeStep;
   updateThicknessByTend(State, NextLevel, State, CurLevel, ZeroTimeStep);
   updateVelocityByTend(State, NextLevel, State, CurLevel,
                        0.5 * StageTimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel,
                       State, CurLevel, ZeroTimeStep);

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

   AuxState->computeAll(State, NextTracerArray, NextLevel, NextLevel,
                        0.5 * StageTimeStep);
   Tend->computeThicknessTendenciesOnly(
       State, AuxState, NextLevel, NextLevel,
       StageTime + 0.5 * StageTimeStep);
   Tend->computeVelocityTendenciesOnly(State, AuxState, NextTracerArray,
                                       NextLevel, NextLevel, NextLevel,
                                       StageTime + 0.5 * StageTimeStep);
   Tend->computeTracerTendenciesOnly(State, AuxState, NextTracerArray,
                                     NextLevel, NextLevel,
                                     StageTime + 0.5 * StageTimeStep);

   updateStateByTend(State, NextLevel, State, CurLevel, StageTimeStep);
   updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel,
                       State, CurLevel, StageTimeStep);

   SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, NextLevel);

   Pacer::stop("SE-RK2:stage3TrThick", 2);
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
      State->exchangeHalo(NextLevel);
      MeshHalo->exchangeFullArrayHalo(NextTracerArray, OnCell);
      Pacer::stop("SE-RK2:haloStage1", 3);

      SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, NextLevel);
      switch (SEConfig.BtrTimeStepper) {
      case SplitExplicitBarotropicStepperType::PredictorCorrector:
         BarotropicPCStepper.doSplitStage2(
             State, SEBuffers, SEConfig, Mesh, VCoord, NextLevel,
             StageTime + 0.5 * TimeStepIterationTimeStep,
             TimeStepIterationTimeStep);
         break;
      case SplitExplicitBarotropicStepperType::Invalid:
      default:
         LOG_CRITICAL("Invalid split-explicit barotropic time stepper");
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

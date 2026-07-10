#ifndef OMEGA_TS_SPLIT_EXPLICIT_RK2_H
#define OMEGA_TS_SPLIT_EXPLICIT_RK2_H
//===-- SplitExplicitRK2Stepper.h - split-explicit RK2 step ---*- C++ -*-===//
//
/// \file
/// \brief Contains the framework for split-explicit RK2 time stepping.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "SplitExplicitTypes.h"
#include "TimeStepper.h"

#include <functional>

namespace OMEGA {

class SplitExplicitRK2Stepper : public TimeStepper {
 public:
   SplitExplicitRK2Stepper(
       const std::string &InName,      ///< [in] name of time stepper
       const TimeInterval &InTimeStep, ///< [in] time step
       const TimeInstant &InStartTime, ///< [in] start time for time stepping
       ///< [in] stop time for time stepping, missing in coupled mode
       std::optional<TimeInstant> InStopTime = std::nullopt);

   /// Advance the state by one split-explicit RK2 step.
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;

   /// Initialize split-explicit state after initial/restart input is read.
   void initializeStateFromInput(
       OceanState *State, ///< [inout] model state after input has been read
       bool ReadRestart   ///< [in] true if restart input initialized the state
   ) const override;

 protected:
   /// Performs additional initialization for split-explicit scratch fields.
   void finalizeInit() override;

 private:
   void doSplitStage1(
       OceanState *State,                  ///< [inout] model state
       const Array3DReal &TendencyTracerArray, ///< [in] tracers for tendencies
       I4 CurLevel,                        ///< [in] current time level
       I4 NextLevel,                       ///< [in] next time level
       const TimeInstant &StageTime,       ///< [in] current stage time
       const TimeInterval &StageTimeStep   ///< [in] current stage time step
   ) const;

   void doSplitStage3(
       OceanState *State,                  ///< [inout] model state
       const Array3DReal &CurTracerArray,  ///< [in] current tracers
       const Array3DReal &NextTracerArray, ///< [out] next tracers
       I4 CurLevel,                        ///< [in] current time level
       I4 NextLevel,                       ///< [in] next time level
       const TimeInstant &StageTime,       ///< [in] current stage time
       const TimeInterval &StageTimeStep,  ///< [in] current stage time step
       bool FinalIteration ///< [in] true on the final time-step iteration
   ) const;

   using BarotropicStage2Function = std::function<void(
       OceanState *, I4, I4, const TimeInstant &, const TimeInterval &)>;

   void initBarotropicStepper();

   void doSplitStage2(OceanState *State, I4 CurLevel, I4 NextLevel,
                      const TimeInstant &StageTime,
                      const TimeInterval &StageTimeStep) const;

   void computeTransportVelocity(
       OceanState *State, ///< [inout] model state
       I4 TimeLevel       ///< [in] time level for split velocity
   ) const;

   void doBaroclinicCoriolisIteration(
       OceanState *State,                   ///< [inout] model state
       const Array2DReal &BaseVelocityTend, ///< [in] non-Coriolis tendency
       I4 CurLevel,                         ///< [in] current time level
       I4 NextLevel,                        ///< [in] next time level
       const TimeInterval &StageTimeStep    ///< [in] current stage time step
   ) const;

   void computeBarotropicForcing(
       OceanState *State,                ///< [inout] model state
       I4 CurLevel,                      ///< [in] current time level
       I4 NewLevel,                      ///< [in] next time level
       const TimeInterval &StageTimeStep ///< [in] current stage time step
   ) const;

   void updateBaroclinicVelocityByTend(
       OceanState *State1, ///< [out] updated state
       I4 TimeLevel1,      ///< [in] time level index for new time
       OceanState *State2, ///< [in] state for current time
       I4 TimeLevel2,      ///< [in] time level index for current time
       TimeInterval Coeff  ///< [in] time-related coeff for tendency
   ) const;

   void
   initializeNextState(OceanState *State, ///< [inout] model state
                       I4 CurLevel,       ///< [in] current time level
                       I4 NextLevel,      ///< [in] next time level
                       I4 SplitFactor     ///< [in] split factor
   ) const;

   void reconstructNormalVelocity(
       OceanState *State, ///< [inout] model state
       I4 CurLevel,       ///< [in] current time level
       I4 NextLevel,      ///< [in] next time level
       bool FinalIteration ///< [in] true on the final time-step iteration
   ) const;

   void finalizeTimeStepIterationState(
       OceanState *State, ///< [inout] model state
       I4 CurLevel,       ///< [in] current time level
       I4 NextLevel,      ///< [in] next time level
       bool FinalIteration ///< [in] true on the final time-step iteration
   ) const;

   void computeVerticalPseudoVelocity(
       OceanState *State,         ///< [inout] model state
       I4 ThickTimeLevel,         ///< [in] thickness time level
       I4 VelTimeLevel,           ///< [in] reconstructed velocity time level
       TimeInterval StageTimeStep ///< [in] current stage time step
   ) const;
   void computeVerticalPseudoVelocity(
       OceanState *State,                 ///< [inout] model state
       I4 ThickTimeLevel,                 ///< [in] thickness time level
       const Array2DReal &NormalVelEdge,  ///< [in] velocity on edges
       TimeInterval StageTimeStep         ///< [in] current stage time step
   ) const;

   SplitExplicitConfig SEConfig;
   mutable SplitExplicitScratch SEScratch;
   BarotropicStage2Function BarotropicStage2;
   SplitExplicitBarotropicPCStepper BarotropicPCStepper;
};

} // namespace OMEGA

#endif

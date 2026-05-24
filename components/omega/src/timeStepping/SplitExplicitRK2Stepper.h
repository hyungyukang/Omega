#ifndef OMEGA_TS_SPLIT_EXPLICIT_RK2_H
#define OMEGA_TS_SPLIT_EXPLICIT_RK2_H
//===-- SplitExplicitRK2Stepper.h - split-explicit RK2 step ---*- C++ -*-===//
//
/// \file
/// \brief Contains the framework for split-explicit RK2 time stepping.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "SplitExplicitData.h"
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

 protected:
   /// Performs additional initialization for split-explicit scratch fields.
   void finalizeInit() override;

 private:
   void doSplitStage1(
       OceanState *State,                   ///< [inout] model state
       const Array3DReal &CurTracerArray,   ///< [in] current tracers
       const Array3DReal &NextTracerArray,  ///< [out] provisional tracers
       I4 CurLevel,                         ///< [in] current time level
       I4 NextLevel,                        ///< [in] next time level
       const TimeInstant &StageTime,        ///< [in] current stage time
       const TimeInterval &StageTimeStep    ///< [in] current stage time step
   ) const;

   void doSplitStage3(
       OceanState *State,                   ///< [inout] model state
       const Array3DReal &CurTracerArray,   ///< [in] current tracers
       const Array3DReal &NextTracerArray,  ///< [out] next tracers
       I4 CurLevel,                         ///< [in] current time level
       I4 NextLevel,                        ///< [in] next time level
       const TimeInstant &StageTime,        ///< [in] current stage time
       const TimeInterval &StageTimeStep    ///< [in] current stage time step
   ) const;

   using BarotropicStage2Function = std::function<void(
       OceanState *, I4, const TimeInstant &, const TimeInterval &)>;

   void initBarotropicStepper();
   void doSplitStage2(OceanState *State, I4 TimeLevel,
                      const TimeInstant &StageTime,
                      const TimeInterval &StageTimeStep) const;

   SplitExplicitConfig SEConfig;
   mutable SplitExplicitScratch SEScratch;
   BarotropicStage2Function BarotropicStage2;
   SplitExplicitBarotropicPCStepper BarotropicPCStepper;
};

} // namespace OMEGA

#endif

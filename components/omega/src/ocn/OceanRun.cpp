//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "OceanDriver.h"
#include "OceanState.h"
#include "TimeStepper.h"

namespace OMEGA {

int ocnRun(TimeInstant &CurrTime, ///< [inout] current sim time
           Alarm &EndAlarm        ///< [inout] alarm to end simulation
) {

   // error code
   I4 Err = 0;

   // fetch default OceanState and TimeStepper
   OceanState *DefOceanState   = OceanState::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();

   TimeInterval TimeStep, ZeroInterval;

   // set simulation clock and attach EndAlarm
   TimeStep = DefTimeStepper->getTimeStep();
   Clock OmegaClock(CurrTime, TimeStep);
   Err = OmegaClock.attachAlarm(&EndAlarm);

   if (TimeStep == ZeroInterval) {
      LOG_ERROR("ocnRun: TimeStep must be initialized");
      ++Err;
   }

   I8 IStep = 0;

   ///////////////////////////////
   TimeInstant StartTime = OmegaClock.getStartTime();

   TimeFrac ElapsedTime = StartTime.getElapsedTime();

   R8 ElapsedTimeSec;
   ElapsedTimeSec = ElapsedTime.getSeconds();

   TimeInterval RefInterval;

   Err = RefInterval.set(ElapsedTimeSec,TimeUnits::Seconds);
   ///////////////////////////////

   ///////////////////////////////////////
   Err = IO_init(StartTime);
   ///////////////////////////////////////

   // time loop, integrate until EndAlarm or error encountered
   while (Err == 0 && !(EndAlarm.isRinging())) {

      // advance clock
      OmegaClock.advance();
      ++IStep;

      // call forcing routines, anything needed pre-timestep

      // do forward time step
      TimeInstant SimTime = OmegaClock.getPreviousTime();
      ///////////////////////////////
      //TimeInstant SimTime = OmegaClock.getPreviousTime() - RefInterval;
      //TimeInstant SimTime = OmegaClock.getPreviousTime() - StartTime;

      TimeInstant CurrRunDuration = SimTime - RefInterval;
      ///////////////////////////////
      //
      DefTimeStepper->doStep(DefOceanState, CurrRunDuration);

      // write restart file/output, anything needed post-timestep

      CurrTime = OmegaClock.getCurrentTime();
      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               CurrTime.getString(4, 4, "-"));
   }


   // write restart file/output
   ///////////////////////////////////////
   Err = IO_write(CurrTime);
   ///////////////////////////////////////

   return Err;

} // end ocnRun

} // end namespace OMEGA

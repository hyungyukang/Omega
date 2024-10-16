//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "OceanDriver.h"
#include "OceanState.h"
#include "TimeStepper.h"
#include "IOStream.h"

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

   // Ocean IO info init
   Err = ocnIOInit(OmegaClock);

   if (TimeStep == ZeroInterval) {
      LOG_ERROR("ocnRun: TimeStep must be initialized");
      ++Err;
   }

   I8 IStep = 0;

   // For the manufactured solution ---------------------
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   I4 NCellsAll = DefHorzMesh->NCellsAll;
   I4 NCellsSize = DefHorzMesh->NCellsSize;
   Array1DR8 ssh("ssh",NCellsSize);
   // Check if use the manufactured solution test
   bool UseManufacturedSolution = false;
   Config *OmegaConfig = Config::getOmegaConfig();
   Config MSConfig("ManufacturedSolution");
   Err = OmegaConfig->get(MSConfig);
   if ( Err != 0 ) {
      LOG_ERROR("Tendencies:: ManufacturedSolution group not found in Config");
      return Err;
   }

   Err = MSConfig.get("UseManufacturedSolution", UseManufacturedSolution);
   if ( Err != 0 ) {
      LOG_ERROR("Tendencies:: ManufacturedSolution: UseManufacturedSolution "
                "not found in Config");
      return Err;
   }
   R8 RestingThickness = 0.0;
   if ( Err != 0 ) {
      LOG_ERROR("Tendencies:: ManufacturedSolution: UseManufacturedSolution "
                "not found in Config");
      return Err;
   } else {
      Err = MSConfig.get("RestingThickness", RestingThickness);
   }
   //----------------------------------------------------

   // time loop, integrate until EndAlarm or error encountered
   while (Err == 0 && !(EndAlarm.isRinging())) {

      // advance clock
      OmegaClock.advance();
      ++IStep;

      // call forcing routines, anything needed pre-timestep

      // do forward time step
      TimeInstant SimTime = OmegaClock.getPreviousTime();
      DefTimeStepper->doStep(DefOceanState, SimTime);


      //SSH field update
      if (UseManufacturedSolution) {
         for (int ICell = 0 ; ICell < NCellsAll; ++ICell) {
            ssh(ICell) = DefOceanState->LayerThickness[0](ICell,0) - RestingThickness;
         }
         Err = Field::attachFieldData<Array1DR8>("ssh", ssh);
      }

      // write restart file/output, anything needed post-timestep
      Err = IOStream::writeAll(OmegaClock);

      CurrTime = OmegaClock.getCurrentTime();
      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               CurrTime.getString(4, 4, "-"));
   }

   // Finalize IOStream
   IOStream::finalize(OmegaClock);


   return Err;

} // end ocnRun

} // end namespace OMEGA

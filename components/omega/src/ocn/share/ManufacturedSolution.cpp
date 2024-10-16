//===-- ocn/share/ManufacturedSolution.cpp - Manufactured Solution -*- C++ -*-===//
//
// \file
// \brief Contains a function for initializing the manufactured tendencies
//
//===-------------------------------------------------------------------------===//

#include "ManufacturedSolution.h"

//===--------------------------------------------------------------------------===/

//using namespace OMEGA;
namespace OMEGA {
//===--------------------------------------------------------------------------===/

int ManufacturedSolution::init() {

   int Err = 0;

   // Get ManufacturedSolution Config group & some constants from Config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config MSConfig("ManufacturedSolution");

   Err = OmegaConfig->get(MSConfig);
   if ( Err != 0 ) {
      LOG_ERROR("Error in ManufacturedSolution group Config");
   }

   // Get parameters from Config
   Err = MSConfig.get("Amplitude",        ManufacturedThickTend.Amplitude);
   Err = MSConfig.get("Amplitude",        ManufacturedVelTend.Amplitude);
   Err = MSConfig.get("RestingThickness", ManufacturedThickTend.RestingThickness);
   Err = MSConfig.get("RestingThickness", ManufacturedVelTend.RestingThickness);

   if ( Err != 0 ) {
      LOG_ERROR("ManufacturedSolution::Error reading ManufacturedSolution "
                " parameters in Config");
   }


   LOG_INFO("ElapsedInMSinit");

   // Compute RefInterval since reference time
   // by retrieving start time from config

   // Get TimeManagement Config group
   Config TimeConfig("TimeManagement");
   Err = OmegaConfig->get(TimeConfig);
   std::string ConfigCalStr;
   Err = TimeConfig.get("CalendarType", ConfigCalStr);

   CalendarKind ConfigCalKind = CalendarUnknown;
   I4 ICalType                = CalendarUnknown;
   for (I4 I = 0; I < NUM_SUPPORTED_CALENDARS; ++I) {
      if (ConfigCalStr == CalendarKindName[I]) {
         ICalType      = I;
         ConfigCalKind = (CalendarKind)(ICalType + 1);
         break;
      }
   }
   Calendar MSCal;
   MSCal.~Calendar();
   MSCal = Calendar(ConfigCalStr, ConfigCalKind);

   // retrieve start time from config
   std::string StartTimeStr;
   Err = TimeConfig.get("StartTime", StartTimeStr);
   if (Err != 0) {
      LOG_CRITICAL("ocnInit: StartTime not found in TimeConfig");
      return Err;
   }
   TimeInstant MSStartTime = TimeInstant(&MSCal, StartTimeStr);


   // Get reference elapsed time
   ManufacturedThickTend.RefElapsedTime = MSStartTime.getElapsedTime();
   ManufacturedVelTend.RefElapsedTime = MSStartTime.getElapsedTime();

   MSCal.~Calendar();
} // end init

//===--------------------------------------------------------------------------===/
} // namespace OMEGA

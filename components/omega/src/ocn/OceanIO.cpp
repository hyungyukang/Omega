//===-- ocn/OceanInit.cpp - Ocean Initialization ----------------*- C++ -*-===//
//
// This file contians ocnInit and associated methods which initialize Omega.
// The ocnInit process reads the config file and uses the config options to
// initialize time management and call all the individual initialization
// routines for each module in Omega.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "IOStream.h"
#include "Dimension.h"
#include "Tracers.h"

#include "mpi.h"

namespace OMEGA {


//------------------------------------------------------------------------------
// A simple test evaluation function
template <typename T>
void TestEval(const std::string &TestName, T TestVal, T ExpectVal, int &Error) {

   if (TestVal == ExpectVal) {
      LOG_INFO("{}: PASS", TestName);
   } else {
      LOG_ERROR("{}: FAIL", TestName);
      ++Error;
   }
}

int ocnIOInit(Clock &OmegaClock) {

   I4 Err = 0;

   LOG_INFO("Entering ocnIOInit");

   ///////////////////////////////////////////////////////////////////
   Err = IOStream::init(OmegaClock);
   ///////////////////////////////////////////////////////////////////

   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   auto *DefState = OceanState::getDefault();

   I4 NCellsOwned  = DefHorzMesh->NCellsOwned;
   I4 NCellsAll   = DefHorzMesh->NCellsAll;
   I4 NCellsSize  = DefHorzMesh->NCellsSize;
   I4 NEdgesSize  = DefHorzMesh->NEdgesSize;
   I4 NVertLevels = DefHorzMesh->NVertLevels;

   LOG_INFO("NVertLevels {}",NVertLevels);

   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);

   std::shared_ptr<Dimension> TimeDim =
       Dimension::create("Time", 1);

   std::shared_ptr<Dimension> StrDim =
       Dimension::create("StrLen", 64);



   I4 Err1 = 0;
   I4 ErrRef = 0;

   //-----------------------------------

   std::shared_ptr<Field> SimField  = Field::get(SimMeta);

   TimeInstant SimStartTime = OmegaClock.getStartTime();

   std::string StartTimeStr = SimStartTime.getString(4, 2, "_");
   Err1 = SimField->addMetadata("SimStartTime", StartTimeStr);
   TestEval("Add SimStartTime metadata", Err1, ErrRef, Err);

   Real FillValue = -1.2345e-30;
   std::string TimeFillValue;
   TimeFillValue = "9999-99-99_99:99:99";

   // Define Xtime ---------------------------------
   //std::vector<std::string> DimNames1(1);
   //DimNames1[0] = "Time";
   //auto TimeField = Field::create(
   //    "xtime", "model time, with format \'YYYY-MM-DD_HH:MM:SS\'", "none",
   //    "model_time",  "0000-00-00_00:00:00", "9999-99-99_99:99:99", TimeFillValue, 1, DimNames1);

   // Define state fields --------------------------
   std::vector<std::string> DimNames3(3);
   // NormalVelocity
   DimNames3[0] = "Time";
   DimNames3[1] = "NEdges";
   DimNames3[2] = "NVertLevels";
   // 2D Fields on device
   DimNames3[0]    = "Time";
   DimNames3[1]    = "NEdges";
   DimNames3[2]    = "NVertLevels";
   auto VelField = Field::create(
       "normalVelocity", "Normal velocity at edges", "m/s",
       "normal_velocity",  -100.0, 100.0, FillValue, 3, DimNames3);

   // LayerThickness
   DimNames3[0] = "Time";
   DimNames3[1] = "NCells";
   DimNames3[2] = "NVertLevels";
   // 2D Fields on device
   DimNames3[0] = "Time";
   DimNames3[1]    = "NCells";
   DimNames3[2]    = "NVertLevels";
   auto ThickField =
       Field::create("layerThickness", "Layer thickness at cell centers", "",
                     "layer_thickness", 0.0, 100.0, FillValue, 3, DimNames3);


   // Sea surface height
   std::vector<std::string> DimNames2(2);
   // NormalVelocity
   DimNames2[0] = "Time";
   DimNames2[1] = "NCells";
   // 2D Fields on device
   DimNames2[0] = "Time";
   DimNames2[1] = "NCells";
   auto SshField = Field::create(
       "ssh", "Sea surface height at cells", "m",
       "sea_surface_height",  0.0, 100.0, FillValue, 2, DimNames2);

//   // Define tracer fields --------------------------
//   DimNames[0] = "Time";
//   DimNames[1] = "NCells";
//   DimNames[2] = "NVertLevels";
//
//   // Temperature
//   auto TempField = Field::create(
//       "Temperature", "Potential temperature at cell centers", "deg C",
//       "sea_water_pot_tem", -3.0, 100.0, FillValue, 3, DimNames);
//
//   // Salinity
//   auto SaltField =
//       Field::create("Salinity", "Salinity at cell centers", "",
//                     "sea_water_salinity", 0.0, 100.0, FillValue, 3, DimNames);
//
//   // tracer1
//   auto tracer1Field =
//       Field::create("tracer1", "tracer1 at cell centers", "",
//                     "none", 0.0, 100.0, FillValue, 3, DimNames);
//

   // Create Time group --------------------------
   //auto TimeGroup = FieldGroup::create("Time");
   //Err1 = TimeGroup->addField("xtime");
   //TestEval("Add xtime to time group", Err1, ErrRef, Err);

   // Create State group --------------------------
   auto StateGroup = FieldGroup::create("States");

   // Add fields to State group
   Err1 = StateGroup->addField("normalVelocity");
   TestEval("Add normalVelocity to state group", Err1, ErrRef, Err);
   Err1 = StateGroup->addField("layerThickness");
   TestEval("Add layerThickness to state group", Err1, ErrRef, Err);
   Err1 = StateGroup->addField("ssh");
   TestEval("Add ssh to state group", Err1, ErrRef, Err);

//   // Create Tracer group --------------------------
//   auto TracerGroup = FieldGroup::create("Tracers");
//
//   // Add fields to tracer group
//   Err1 = TracerGroup->addField("Temperature");
//   TestEval("Add Temperature to tracer group", Err1, ErrRef, Err);
//   Err1 = TracerGroup->addField("Salinity");
//   TestEval("Add Salinity to tracer group", Err1, ErrRef, Err);
//   Err1 = TracerGroup->addField("tracer1");
//   TestEval("Add tracer1 to tracer group", Err1, ErrRef, Err);

//--------------------------------------------

   // Also create a Restart group
   auto RestartGroup = FieldGroup::create("Restart");

   // Add fields to restart group
   Err1 = RestartGroup->addField("normalVelocity");
   TestEval("Add normalVelocity to restart group", Err1, ErrRef, Err);
   Err1 = RestartGroup->addField("layerThickness");
   TestEval("Add layerThickness to restart group", Err1, ErrRef, Err);

//   Err1 = RestartGroup->addField("Temperature");
//   TestEval("Add Temperature to restart group", Err1, ErrRef, Err);
//   Err1 = RestartGroup->addField("Salinity");
//   TestEval("Add Salinity to restart group", Err1, ErrRef, Err);

//-----------------------------------

   // Retrieve dimension lengths
   //I4 NCellsSize  = Dimension::getDimLengthLocal("NCells");
   //I4 NVertLevels = Dimension::getDimLengthLocal("NVertLevels");

   // Create data arrays

   //Array2DR8 Temp("Temp", NCellsSize, NVertLevels);
   //Array2DR8 Salt("Salt", NCellsSize, NVertLevels);


   // Retrieve RestingThickness for the manufactured solution
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

   // Attach data arrays to fields ------------------------------------
   //std::vector<std::string> Xtime(1);
   //std::string Xtime;
   //Err1 = Field::attachFieldData<std::string>("xtime",Xtime);

   Err1 = Field::attachFieldData<Array2DR8>("normalVelocity", DefState->NormalVelocity[0]);
   TestEval("Attach normalVelocity data to field", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData<Array2DR8>("layerThickness", DefState->LayerThickness[0]);
   TestEval("Attach layerThickness data to field", Err1, ErrRef, Err);

////////////////////////////////////
   Array1DR8 ssh("ssh",NCellsSize);
   //SSH field update
   if (UseManufacturedSolution) {
      for (int ICell = 0 ; ICell < NCellsAll; ++ICell) {
         ssh(ICell) = DefState->LayerThickness[0](ICell,0) - RestingThickness;
      }
      Err = Field::attachFieldData<Array1DR8>("ssh", ssh);
   }
   TestEval("Attach ssh data to field", Err1, ErrRef, Err);
////////////////////////////////////

   //Err1 = Field::attachFieldData<Array2DR8>("Temperature", Tracers::getFieldByIndex(0)->getDataArray<Array2DReal>());
   //TestEval("Attach temperature data to field", Err1, ErrRef, Err);
   //Err1 = Field::attachFieldData<Array2DR8>("Salinity", Tracers::getFieldByIndex(1)->getDataArray<Array2DReal>());
   //TestEval("Attach salinity data to field", Err1, ErrRef, Err);
   //Err1 = Field::attachFieldData<Array2DR8>("tracer1", Tracers::getFieldByIndex(2)->getDataArray<Array2DReal>());
   //TestEval("Attach tracer1 data to field", Err1, ErrRef, Err);

   //for (I4 Cell = 0; Cell < NCellsOwned; Cell++) {
   //   for (I4 Vert = 0; Vert < NVertLevels; Vert++) {
   //   LOG_INFO("Temp {}", TestFieldData(Cell,Vert));
   //   }
   //}

   // Validate all streams (Mesh stream already validated in HorzMesh?)
   bool AllValidated = IOStream::validateAll();
   TestEval("IOStream Validation", AllValidated, true, Err);

   // Try to writeAll for OnStartup
   Err = IOStream::writeAll(OmegaClock);

}


} // end namespace OMEGA
//===----------------------------------------------------------------------===//

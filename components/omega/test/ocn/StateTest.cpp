//===-- Test driver for OMEGA State -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA state class
///
/// This driver tests that the OMEGA state class member variables are read in
/// correctly from a sample shperical mesh file. Also tests that the time level
/// update works as expected.
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for State testing. It calls various
// init routines, including the creation of the default decomposition.

int initStateTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the IO system
   Err = OMEGA::IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("State: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default decomposition");

   // Initialize the default halo
   Err = OMEGA::Halo::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default halo");

   // Initialize the default mesh
   Err = OMEGA::HorzMesh::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default mesh");

   // Initialize the default state
   Err = OMEGA::OceanState::init();
   if (Err != 0)
      LOG_ERROR("State: error initializing default State");

   return Err;
}

//------------------------------------------------------------------------------
// The test driver for State -> This tests the time level update of state
// variables and verifies the state is read in correctly.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // Call initialization routine to create the default decomposition
      int Err = initStateTest();
      if (Err != 0)
         LOG_CRITICAL("State: Error initializing");

      // Get MPI vars if needed
      OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
      MPI_Comm Comm          = DefEnv->getComm();
      OMEGA::I4 MyTask       = DefEnv->getMyTask();
      OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster          = DefEnv->isMasterTask();

      OMEGA::HorzMesh *DefHorzMesh = OMEGA::HorzMesh::getDefault();
      OMEGA::Decomp *DefDecomp     = OMEGA::Decomp::getDefault();
      OMEGA::Halo *DefHalo         = OMEGA::Halo::getDefault();

      // These hard-wired variables need to be upated
      // with retrivals/config options
      int NVertLevels = 60;
      int NTimeLevels = 2;

      // Create "test" state
      OMEGA::OceanState DefOceanState("Test", DefHorzMesh, DefDecomp, DefHalo,
                                      NVertLevels, NTimeLevels);

      // Test retrieval of the default state
      OMEGA::OceanState *DefState = OMEGA::OceanState::getDefault();
      if (DefState) { // true if non-null ptr
         LOG_INFO("State: Default state retrieval PASS");
      } else {
         LOG_INFO("State: Default state retrieval FAIL");
      }

      OMEGA::OceanState *TestState = OMEGA::OceanState::get("Test");
      if (TestState) { // true if non-null ptr
         LOG_INFO("State: Test state retrieval PASS");
      } else {
         LOG_INFO("State: Test state retrieval FAIL");
      }

      // auto LayerThicknessH_0 = DefState->LayerThicknessH[0];
      // auto LayerThicknessH_1 = DefState->LayerThicknessH[1];
      // for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
      //    for (int Level = 0; Level < DefState->NVertLevels; Level++) {
      //       LOG_INFO(LayerThicknessH_0(Cell, Level));
      //       LOG_INFO(LayerThicknessH_1(Cell, Level));
      //    }
      // }

      // Test that reasonable values have been read in for LayerThickness
      int count = 0;
      for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
         int colCount = 0;
         for (int Level = 0; Level < DefState->NVertLevels; Level++) {
            OMEGA::R8 val = DefState->LayerThicknessH[0](Cell, Level);
            if (val > 0.0 && val < 300.0) {
               colCount++;
            }
         }
         if (colCount < 2) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("State: State read PASS");
      } else {
         LOG_INFO("State: State read FAIL");
      }

      // Test that initally the 0 time levels of the
      // Def and Test state arrays match
      count                     = 0;
      auto LayerThicknessH_def  = DefState->LayerThicknessH[0];
      auto LayerThicknessH_test = TestState->LayerThicknessH[0];
      for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
         for (int Level = 0; Level < DefState->NVertLevels; Level++) {
            if (LayerThicknessH_def(Cell, Level) !=
                LayerThicknessH_test(Cell, Level)) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("State: Default test state comparison PASS");
      } else {
         LOG_INFO("State: Default test state comparison FAIL");
      }

      // Test that the time level update is correct.
      DefState->updateTimeLevels();

      count                = 0;
      LayerThicknessH_def  = DefState->LayerThicknessH[1];
      LayerThicknessH_test = TestState->LayerThicknessH[0];
      for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
         for (int Level = 0; Level < DefState->NVertLevels; Level++) {
            if (LayerThicknessH_def(Cell, Level) !=
                LayerThicknessH_test(Cell, Level)) {
               count++;
            }
         }
      }

      LayerThicknessH_def  = DefState->LayerThicknessH[0];
      LayerThicknessH_test = TestState->LayerThicknessH[1];
      for (int Cell = 0; Cell < DefState->NCellsAll; Cell++) {
         for (int Level = 0; Level < DefState->NVertLevels; Level++) {
            if (LayerThicknessH_def(Cell, Level) !=
                LayerThicknessH_test(Cell, Level)) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("State: time level update PASS");
      } else {
         LOG_INFO("State: time level update FAIL");
      }

      // Test time level update on device
      int count1;
      auto LayerThickness_def  = DefState->LayerThickness[1];
      auto LayerThickness_test = TestState->LayerThickness[0];
      OMEGA::parallelReduce(
          "reduce", {DefState->NCellsAll, DefState->NVertLevels},
          KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
             if (LayerThickness_def(Cell, Level) !=
                 LayerThickness_test(Cell, Level)) {
                Accum++;
             }
          },
          count1);

      int count2;
      LayerThickness_def  = DefState->LayerThickness[0];
      LayerThickness_test = TestState->LayerThickness[1];
      OMEGA::parallelReduce(
          "reduce", {DefState->NCellsAll, DefState->NVertLevels},
          KOKKOS_LAMBDA(int Cell, int Level, int &Accum) {
             if (LayerThickness_def(Cell, Level) !=
                 LayerThickness_test(Cell, Level)) {
                Accum++;
             }
          },
          count2);

      if (count1 == 0 && count2 == 0) {
         LOG_INFO("State: time level update (GPU) PASS");
      } else {
         LOG_INFO("State: time level update (GPU) FAIL");
      }

      // Finalize Omega objects
      OMEGA::HorzMesh::clear();
      OMEGA::Decomp::clear();
      OMEGA::MachEnv::removeAll();
      OMEGA::OceanState::clear();
      OMEGA::IOField::clear();

      // MPI_Status status;
      if (Err == 0)
         LOG_INFO("State: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/

//===--   -----------------------------*- C++ -*-===/
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "mpi.h"
#include "sw_common.h"
#include "sw_constants.h"

#include <iostream>
#include <cmath>

using namespace OMEGA;

//===-----------------------------------------------------------------------===/

//constexpr Geometry Geom = Geometry::Spherical;
constexpr char DefaultMeshFile[] = "OmegaSphereMesh.nc";

//===-----------------------------------------------------------------------===/

int sw_initCore(const std::string &MeshFile) {
   int Err = 0;
  
   //--------------------------------------------------------------------------/

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();
 
   int IOErr = IO::init(DefComm);
   if (IOErr != 0 ) {
      Err++;
      LOG_ERROR("SW: error initializing IO");
   }
   int DecompErr = Decomp::init(MeshFile);
   if (DecompErr != 0 ) {
      Err++;
      LOG_ERROR("SW: error initializing Decomp");
   }
   int HaloErr = Halo::init();
   if (HaloErr != 0 ) {
      Err++;
      LOG_ERROR("SW: error initializing Halo");
   }
   int MeshErr = HorzMesh::init(); 
   if (MeshErr != 0 ) {
      Err++;
      LOG_ERROR("SW: error initializing HorzMesh");
   }
   int StateErr = OceanState::init();
   if (StateErr != 0 ) {
      Err++;
      LOG_ERROR("SW: error initializing State");
   }

   //--------------------------------------------------------------------------/

   OceanState *DefState = OMEGA::OceanState::getDefault();
   NCellsAll = DefState->NCellsAll;

   //LOG_INFO("NCellsAll {}", NCellsAll);

   //--------------------------------------------------------------------------/
    
   return Err;
}

//===-----------------------------------------------------------------------===/

//int sw_initState() {

//   int Err = 0;

   //--------------------------------------------------------------------------/

   //const auto &Mesh = HorzMesh::getDefault();

   
   
   //Mesh = OceanState:
   
   /*
   const I4 NCellsOwned = Mesh->NCellsOwned; 
   const I4 NEdgesOwned = Mesh->NEdgesOwned; 
   const I4 NVerticesOwned = Mesh->NVerticesOwned; 

   const I4 NCellsAll = Mesh->NCellsAll; 
   const I4 NEdgesAll = Mesh->NEdgesAll; 
   const I4 NVerticesAll = Mesh->NVerticesAll; 
   */

   //Array2DReal normalVelocity("normalVelocity", NEdgesAll, NVertLevels);
   //Err += setVectorCell(
   //    KOKKOS_LAMBDA(Real(&

   //Array2DReal layerThickness("layerThickness", NCellsAll, NVertLevels);

   /*
   LOG_INFO("NCellsOwned {}", NCellsOwned);
   LOG_INFO("NEdgesOwned {}", NEdgesOwned);
   LOG_INFO("NVerticesOwned {}", NVerticesOwned);
   LOG_INFO("NCellsAll {}", NCellsAll);
   LOG_INFO("NEdgesAll {}", NEdgesAll);
   LOG_INFO("NVerticesAll {}", NVerticesAll);
   */

//   return Err;
//}

//===-----------------------------------------------------------------------===/

void sw_finalize() {
   OceanState::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//===-----------------------------------------------------------------------===/

void sw_init(const std::string &MeshFile = DefaultMeshFile) {

   int initCoreErr = sw_initCore(MeshFile);

   LOG_INFO("NCellsAll in init {}", NCellsAll);

   if (initCoreErr != 0 ) {
      LOG_CRITICAL("Error sw_initCore");
   }

   if (initCoreErr == 0 ) {
      LOG_INFO("Successful sw_initCore");
   }

}

//===-----------------------------------------------------------------------===/

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   
   sw_init();
   sw_finalize();

   Kokkos::finalize();
   MPI_Finalize();
}  // end of main

//===-----------------------------------------------------------------------===/

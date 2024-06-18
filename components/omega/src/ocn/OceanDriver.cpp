//===--   -----------------------------*- C++ -*-===/
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "OceanDriver.h"
#include "ShallowWaterCore.h"
#include "TendencyTerms.h"
#include "mpi.h"
//#include "sw_common.h"
//#include "sw_constants.h"

#include <iostream>
#include <cmath>

//===-----------------------------------------------------------------------===/

namespace OMEGA{

//-----------------------------------------------------------------------------/

void OceanDriver::
get_comm(int argc, char *argv[]) {
   LOG_INFO("OMEGA::initialize");

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc,argv);
} // get_comm

//-----------------------------------------------------------------------------/

void OceanDriver::
init_IO(){
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();
   
   int IOErr = IO::init(DefComm);
   if (IOErr != 0 ) {
      LOG_ERROR("SW: error initializing IO");
   }

} // init_IO 

//-----------------------------------------------------------------------------/

void OceanDriver::
init_Decomp(const std::string &MeshFile){
   int DecompErr = Decomp::init(MeshFile);
   if (DecompErr != 0 ) {
      LOG_ERROR("SW: error initializing Decomp");
   }
} // init_Decomp 

//-----------------------------------------------------------------------------/

void OceanDriver::
init_Halo(){
   int HaloErr = Halo::init();
   if (HaloErr != 0 ) {
      LOG_ERROR("SW: error initializing Halo");
   }
} // init_Halo

//-----------------------------------------------------------------------------/

void OceanDriver::
init_HorzMesh(){
   int MeshErr = HorzMesh::init();
   if (MeshErr != 0 ) {
      LOG_ERROR("SW: error initializing HorzMesh");
   }
} // init_HorzMesh

//-----------------------------------------------------------------------------/

void OceanDriver::
init_OceanState(){
   int StateErr = OceanState::init();
   if (StateErr != 0 ) {
      LOG_ERROR("SW: error initializing State");
   }
} // init_OceanState

//-----------------------------------------------------------------------------/

// A wrapper of the above
void OceanDriver::
initialize (int argc,
            char *argv[],
            const std::string &MeshFile){

   get_comm(argc,argv);
   init_IO();
   init_Decomp(MeshFile);
   init_Halo();
   init_HorzMesh();
   init_OceanState();
   set_params();

}

//-----------------------------------------------------------------------------/

void OceanDriver::
set_params(){

   HorzMesh *DefMesh = HorzMesh::getDefault();
   OceanState *DefState = OceanState::getDefault();

   auto NCellsOwned = DefMesh->NCellsOwned;
   auto NCellsAll   = DefMesh->NCellsAll;
   auto NCellsSize  = DefMesh->NCellsSize; 

   auto NEdgesOwned = DefMesh->NEdgesOwned;
   auto NEdgesAll   = DefMesh->NEdgesAll;
   auto NEdgesSize  = DefMesh->NEdgesSize;

   auto MaxCellsOnEdge = DefMesh->MaxCellsOnEdge;
   auto MaxEdges       = DefMesh->MaxEdges;
   auto MaxEdges2      = DefMesh->MaxEdges2;

   auto NVerticesOwned = DefMesh->NVerticesOwned;
   auto NVerticesAll   = DefMesh->NVerticesAll;
   auto NVerticesSize  = DefMesh->NVerticesSize;
   auto VertexDegree   = DefMesh->VertexDegree;

   auto NVertLevels    = DefState->NVertLevels;

   //LOG_INFO(NEdgesAll);

}

//-----------------------------------------------------------------------------/

void OceanDriver::
run() {
   //LOG_INFO("OceanDriver::run");

   ShallowWater SW;
   SW.sw_run();
   
   //LOG_INFO(layerThickness[0](1,1));
}

//-----------------------------------------------------------------------------/

void OceanDriver::
finalize() {
   LOG_INFO("OMEGA::finalize");

   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
   OceanState::clear();
   IOField::clear();
   Kokkos::finalize();
   MPI_Finalize();
} // finalize

//-----------------------------------------------------------------------------/

} // namespace OMEGA

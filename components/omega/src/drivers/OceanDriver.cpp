//===-----------------------------------------------------------------------===/

#include "OceanDriver.h"

//===-----------------------------------------------------------------------===/

namespace OMEGA{

//-----------------------------------------------------------------------------/

void OceanDriver::
get_comm(int argc, char *argv[]) {

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
} // initialize

//-----------------------------------------------------------------------------/

void OceanDriver::
set_params(){

   HorzMesh *DefMesh = HorzMesh::getDefault();
   OceanState *DefState = OceanState::getDefault();

} // seet_params

//-----------------------------------------------------------------------------/

void OceanDriver::
run() {
   //LOG_INFO("OceanDriver::run");

   MachEnv *Env = MachEnv::getDefaultEnv();
   MPI_Comm Comm = Env->getComm();
   Halo *DefHalo = Halo::getDefault();
   HorzMesh *Mesh = HorzMesh::getDefault();
   Decomp *DefDecomp = Decomp::getDefault();
   OceanState *State = OceanState::getDefault();

   //if ( core == "ShallowWater" ) {;  // <- from Config
   ShallowWaterCore SWCore;
   SWCore.sw_model(Comm, Env, DefDecomp, DefHalo, Mesh, State);
   // }
   
} // run

//-----------------------------------------------------------------------------/

void OceanDriver::
finalize() {
   //LOG_INFO("OMEGA::finalize");

   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
   OceanState::clear();
   IOField::clear();
   Kokkos::finalize();
   MPI_Finalize();
} // finalize

//===-----------------------------------------------------------------------===/
} // namespace OMEGA


//===-----------------------------------------------------------------------===/
//=== OMEGA MAIN -----------------------------------------------------------===/
//===-----------------------------------------------------------------------===/

using namespace OMEGA;

int main(int argc, char *argv[]) {

    // Mesh file name <- From Config/YAML later and can be placed anywhere
    constexpr char DefMeshFile[] = "OmegaSphereMesh.nc";
   
    OceanDriver OD;

    // Model initialize
    OD.initialize(argc,argv,DefMeshFile);

    // Model run
    OD.run();

    // Model finalize
    OD.finalize();

    return 0;

} // end of main

//-----------------------------------------------------------------------------/


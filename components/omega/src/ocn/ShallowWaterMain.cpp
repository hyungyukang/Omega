//===--   -----------------------------*- C++ -*-===/
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
//#include "DataTypes.h"
//#include "Decomp.h"
//#include "Halo.h"
//#include "HorzMesh.h"
//#include "HorzOperators.h"
//#include "IO.h"
//#include "Logging.h"
//#include "MachEnv.h"
//#include "OceanState.h"
//#include "OmegaKokkos.h"
#include "OceanDriver.h"
#include "ShallowWaterCore.h"
//#include "TendencyTerms.h"
#include "mpi.h"
//#include "sw_common.h"
//#include "sw_constants.h"

#include <iostream>
#include <cmath>


//===-----------------------------------------------------------------------===/

using namespace OMEGA;

constexpr char DefMeshFile[] = "OmegaSphereMesh.nc";

//===-----------------------------------------------------------------------===/

int main(int argc, char *argv[]) {
    
    OceanDriver OD;

    OD.initialize(argc,argv,DefMeshFile);

    OD.run();

    OD.finalize();
 
    return 0;
//===-----------------------------------------------------------------------===/

} // end of main

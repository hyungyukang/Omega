#ifndef OMEGA_OCEANDRIVER_H
#define OMEGA_OCEANDRIVER_H

//===-----------------------------------------------------------------------===/

#include "Config.h"
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
#include "TendencyTerms.h"
#include "model/ShallowWaterCore.h"
#include "mpi.h"

//===-----------------------------------------------------------------------===/

namespace OMEGA{

class OceanDriver
{

public:
   OceanDriver () = default;
    
   // Get default communicator
   void get_comm (int argc, char *argv[]);

   // Initialize IO
   void init_IO ();
  
   // Initialize Decomp with MeshFile
   void init_Decomp(const std::string &MeshFile);

   // Initialize Halo
   void init_Halo();

   // Initialize Horizontal mesh
   void init_HorzMesh();

   // Initialize ocean state
   void init_OceanState();

   // A wrapper of the above
   void initialize(int argc,
                   char *argv[],
                   const std::string &MeshFile);

   // Set base parameters (now defunct)
   void set_params();

   // Run Omega
   void run();

   // Finalize Omega
   void finalize ();

//-----------------------------------------------------------------------------/

}; // class OmegaDriver

} // namespace OMEGA
#endif

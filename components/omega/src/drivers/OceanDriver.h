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
    
   void get_comm (int argc, char *argv[]);

   void init_IO ();
  
   void init_Decomp(const std::string &MeshFile);

   void init_Halo();

   void init_HorzMesh();

   void init_OceanState();

   void init_Create();

   // A wrapper of the above
   void initialize(int argc,
                   char *argv[],
                   const std::string &MeshFile);

   void set_params();

   void run();

   void finalize ();

//-----------------------------------------------------------------------------/

}; // class OmegaDriver

} // namespace OMEGA

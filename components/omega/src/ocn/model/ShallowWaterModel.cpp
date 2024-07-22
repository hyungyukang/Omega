//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_model(int Comm, const MachEnv *Env, const Decomp *DefDecomp,
         Halo *DefHalo, HorzMesh *Mesh, OceanState *State) {

   // SW var initialization ---------------------------------------------------/
   sw_init_var(Mesh, State);

   // SW IO initialization ----------------------------------------------------/
   sw_init_io(DefDecomp, Mesh, State);

   // Initial field create ----------------------------------------------------/
   sw_init_field(TestCase, Env, DefHalo, Mesh, State );

   // Time stepping -----------------------------------------------------------/
   for (int step = 0; step < nsteps; ++step) {
       // Compute current time
       CurrentTime = (step) * dt;
  
       // Do step
       sw_time_stepper(time_integrator, CurrentTime, dt, Env, DefHalo, Mesh, State);

       // Check current status
       sw_check_status(printInterval, CurrentTime, dt, Comm, Env, Mesh, State);
   }
   //--------------------------------------------------------------------------/

   // Writing final state fields ----------------------------------------------/
   CurrentTime = endTime;//-0.05*dt;
   sw_io_write(DefDecomp, Mesh, State);

} // sw_model

//===-----------------------------------------------------------------------===/

} // namespace OMEGA

//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_check_status(const int printInterval, R8 t, const R8 dt, int Comm, 
                const MachEnv *Env, const HorzMesh *Mesh, OceanState *State) {

   //-----------------------------------------------------------------------------------/
   // Check current status
   //-----------------------------------------------------------------------------------/

   // If this test is using time-dependent solutions
   if ( TimeDependentSolution ) {
      sw_time_dependent_solution(TestCase, t+dt, Mesh, State);
   }

   int ErrSWVelError = 0;
   ErrorMeasures SWVelError;
   ErrSWVelError += computeErrors(SWVelError, State->NormalVelocity[0],
                                  NormalVelocitySolution, Mesh, OnEdge, NVertLevels);

   int ErrSWThickError = 0;
   ErrorMeasures SWThickError;
   ErrSWThickError += computeErrors(SWThickError, State->LayerThickness[0],
                                    LayerThicknessSolution, Mesh, OnCell, NVertLevels);

   if ( Env->getMyTask() == 0 ) {
      R8 time_print = (int) (t+dt) % printInterval;
      if ( time_print == 0 ) {
         //LOG_INFO("{} {} {}", (t+dt)/86400.0, SWVelError.L2, SWThickError.L2);

         cout<<fixed<<setprecision(15);
         cout << (t+dt)/86400.0 << " " << SWVelError.L2 << " " << SWThickError.L2 << '\n';
      }
   }

//   R8 maxNormVelLocal = maxVal(State->NormalVelocity[0]);
//   R8 maxLayThickLocal= maxVal(State->LayerThickness[0]);
//
//   Real maxNormVel, maxLayThick;
//
//   int ErrMaxVal = MPI_Allreduce(&maxNormVelLocal,&maxNormVel,1,MPI_RealKind,MPI_MAX, Comm);
//   int ErrMaxThick = MPI_Allreduce(&maxLayThickLocal,&maxLayThick,1,MPI_RealKind,MPI_MAX, Comm);
//
//   if ( Env->getMyTask() == 0 ) {
//      R8 time_print = (int)(t+dt) % 3600;
//      if ( time_print == 0 ) {
//         LOG_INFO("sw_timeStepper:: current Time = {}", t+dt);
//         LOG_INFO("sw_timeStepper:: after Min. NormalVelocity = {}", maxNormVel);
//         LOG_INFO("sw_timeStepper:: after Min. LayerThickness = {}", maxLayThick);
//      }
//   }

   //--------------------------------------------------------------------------/

} // sw_check_status

} // namespace OMEGA

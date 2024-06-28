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
   R8 time_print = (int) (t+dt) % printInterval;
   if ( time_print == 0 ) {

      // Print current time first
      if ( Env->getMyTask() == 0 ) {
         cout<<fixed<<setprecision(15);
         cout << (t+dt)/86400.0 ;
      }

      // If this test is using time-dependent solutions --------------------------------/
      if ( TimeDependentSolution ) {
         sw_time_dependent_solution(TestCase, t+dt, Mesh, State);
      }

      // Compute L2 and Linf error norms -----------------------------------------------/
      if ( ComputeNormError ) {
         int ErrSWVelError = 0;
         ErrorMeasures SWVelError;
         ErrSWVelError += computeErrors(SWVelError, State->NormalVelocity[0],
                                        NormalVelocitySolution, Mesh, OnEdge, NVertLevels);
   
         int ErrSWThickError = 0;
         ErrorMeasures SWThickError;
         ErrSWThickError += computeErrors(SWThickError, State->LayerThickness[0],
                                          LayerThicknessSolution, Mesh, OnCell, NVertLevels);
   
         if ( Env->getMyTask() == 0 ) {
               cout<<fixed<<setprecision(15);
               cout << " " << SWVelError.L2 << " " << SWThickError.L2;
            }
      } // ComputeNormError


      // Compute total energy ----------------------------------------------------------/
      if ( ComputeTotalEnergy ) {
         Real TotalPE = 0.0_Real;
         Real TotalKE = 0.0_Real;
   
         // KE part
         for (int IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
            const Real AreaEdge = 0.5_Real * Mesh->DvEdge(IEdge) * Mesh->DcEdge(IEdge);
            TotalKE += AreaEdge * LayerThicknessEdge(IEdge,0) *
                       0.5_Real * pow(State->NormalVelocity[0](IEdge,0),2);
         }

         // PE part 
         for (int ICell = 0; ICell < Mesh->NCellsOwned; ++ICell) {
            const Real AreaCell = Mesh->AreaCell(ICell);
            const Real LayerThick = State->LayerThickness[0](ICell,0);
            TotalPE += AreaCell * (  gravity * LayerThick *
                                   (0.5_Real * LayerThick + BottomTopography(ICell)) );
         }
   
         Real TotalEnergy = TotalKE + TotalPE;
         Real TotalEnergyGlobal = 0.0_Real;

         int Err1 = MPI_Allreduce(&TotalEnergy,&TotalEnergyGlobal,1,MPI_RealKind,MPI_SUM, Comm);
   
         if ( Env->getMyTask() == 0 ) {
            R8 time_print = (int) (t+dt) % printInterval;
            if ( time_print == 0 ) {
               cout << " " << TotalEnergyGlobal;
            }
         }
      } // ComputeTotalEnergy

      // Line breaking
      if ( Env->getMyTask() == 0 ) {
         cout << " " << '\n';
      }

   // Print Max val
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

   } // if time_print

   //--------------------------------------------------------------------------/

} // sw_check_status

} // namespace OMEGA

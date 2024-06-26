//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

//===-----------------------------------------------------------------------===/
//=== Shallow Water Time Stepper -------------------------------------------===/
//===-----------------------------------------------------------------------===/
//
// Choices: forward-euler, forward-backward, heuns, ssp-rk3
//

void ShallowWaterCore::
sw_time_stepper(const std::string &time_integrator, R8 t, const R8 dt, const MachEnv *Env,
                const Halo *Halo, const HorzMesh *Mesh, OceanState *State) {

   //===--------------------------------------------------------------------===/
   if ( time_integrator == "forward-euler" ) {
   //===--------------------------------------------------------------------===/

      // Thickness tendency at (n) --------------------------------------------/
      sw_tend_thick(0, 0, Env, Halo, Mesh, State);

      // Normal Velocity tendency at (n) --------------------------------------/
      sw_tend_vel(0, 0, Env, Halo, Mesh, State);

      // State advance to (n+0.5) ---------------------------------------------/
      // LayerThickness
      parallelFor(
         {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = (State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel));
         });

      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = (State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();

   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "forward-backward" ) {
   //===--------------------------------------------------------------------===/

      // Thickness tendency at (n) --------------------------------------------/
         // Thick time level = 0 ; Vel time level = 0
      sw_tend_thick(0, 0, Env, Halo, Mesh, State);

         // State advance to (n+1)
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = (State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel));
         });

      // Normal Velocity tendency at (n) with h_n+1 ---------------------------/
         // Thick time level = 1 ; Vel time level = 0
      sw_tend_vel(1, 0, Env, Halo, Mesh, State);

         // State advance to (n+1)
      parallelFor(
         {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = (State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();


   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "heuns" ) {
   //===--------------------------------------------------------------------===/

      // Stage 1 --------------------------------------------------------------/
      // Thickness tendency
      sw_tend_thick(0, 0, Env, Halo,Mesh,State);

      // State adavance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel);
         });

      // Normal Velocity tendency
      sw_tend_vel(1, 0, Env, Halo,Mesh,State);

      // State adavance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 2 --------------------------------------------------------------/
      // Thickness tendency
      sw_tend_thick(1, 1, Env, Halo,Mesh,State);

      // State advance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
               + 0.5*(State->LayerThickness[1](ICell,KLevel) - State->LayerThickness[0](ICell,KLevel)
                      + dt * LayerThicknessTend(ICell,KLevel));
         });

      // Normal Velocity tendency
      sw_tend_vel(1, 1, Env, Halo,Mesh,State);

      // State advance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
               + 0.5*(State->NormalVelocity[1](IEdge,KLevel) - State->NormalVelocity[0](IEdge,KLevel)
                      + dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();



   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "ssp-rk3" ) {
   //===--------------------------------------------------------------------===/

      // Stage 1 --------------------------------------------------------------/
      // Tendencies
      sw_tend_thick(0, 0, Env, Halo, Mesh, State);
      sw_tend_vel(  0, 0, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel);
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 2 --------------------------------------------------------------/
      // Tendencies
      sw_tend_thick(1, 1, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) =
               (3/4.)*   State->LayerThickness[0](ICell,KLevel)
             + (1/4.)* ( State->LayerThickness[1](ICell,KLevel)
                        +dt * LayerThicknessTend(ICell,KLevel));
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) =
               (3/4.)*   State->NormalVelocity[0](IEdge,KLevel)
             + (1/4.)* ( State->NormalVelocity[1](IEdge,KLevel)
                        +dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Stage 3 --------------------------------------------------------------/
      // Tendencies
      sw_tend_thick(1, 1, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsOwned,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) =
               (1/3.)*   State->LayerThickness[0](ICell,KLevel)
             + (2/3.)* ( State->LayerThickness[1](ICell,KLevel)
                        +dt * LayerThicknessTend(ICell,KLevel));
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) =
               (1/3.)*   State->NormalVelocity[0](IEdge,KLevel)
             + (2/3.)* ( State->NormalVelocity[1](IEdge,KLevel)
                        +dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();

   //--------------------------------------------------------------------------/
   } else {
      LOG_ERROR("Invalid choice of time_integrator \
                (Choices: forward-euler, forward-backward, heuns, ssp-rk3)");
   } // if time_integrator


   //--------------------------------------------------------------------------/
   // Check current status
   //--------------------------------------------------------------------------/

//   int ErrSWVelError = 0;
//   int ErrSWThickError = 0;
//
//   // parallelFor(
//   //    {Mesh->NCellsOwned,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
//   //       //LOG_INFO(SWVar->LayerThicknessInit(ICell,KLevel));
//   //       LOG_INFO(LayerThicknessInit(ICell,KLevel));
//   //    });
//   ErrorMeasures SWVelError;
//   ErrSWVelError += computeErrors(SWVelError, State->NormalVelocity[0],
//                        NormalVelocityInit, Mesh, OnEdge, NVertLevels);
//
//   ErrorMeasures SWThickError;
//   ErrSWThickError += computeErrors(SWThickError, State->LayerThickness[0],
//                        LayerThicknessInit, Mesh, OnCell, NVertLevels);
//
//   if ( Env->getMyTask() == 0 ) {
//      R8 time_print = (int) (t+dt) % 3600;
//      if ( time_print == 0 ) {
//         //LOG_INFO("{} {} {}", (t+dt)/86400.0, SWVelError.L2, SWThickError.L2);
//
//         cout<<fixed<<setprecision(15);
//         cout << (t+dt)/86400.0 << " " << SWVelError.L2 << " " << SWThickError.L2 << '\n';
//      }
//   }

} // sw_time_stepper

} // namespace OMEGA

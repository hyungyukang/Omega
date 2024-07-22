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
                Halo *Halo, HorzMesh *Mesh, OceanState *State) {

   //===--------------------------------------------------------------------===/
   if ( time_integrator == "forward-euler" ) {
   //===--------------------------------------------------------------------===/

      R8 nowTime = t; 

      // Thickness tendency at (n) --------------------------------------------/
      sw_tend_thick(0, 0, nowTime, Env, Halo, Mesh, State);

      // Normal Velocity tendency at (n) --------------------------------------/
      sw_tend_vel(0, 0, nowTime, Env, Halo, Mesh, State);

      // State advance to (n+1) -----------------------------------------------/
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

      R8 nowTime = t; 

      // Thickness tendency at (n) --------------------------------------------/
         // Thick time level = 0 ; Vel time level = 0
      sw_tend_thick(0, 0, nowTime, Env, Halo, Mesh, State);

         // State advance to (n+1)
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = (State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel));
         });

      // Normal Velocity tendency at (n) with h_n+1 ---------------------------/
         // Thick time level = 1 ; Vel time level = 0
      sw_tend_vel(1, 0, nowTime, Env, Halo, Mesh, State);

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
      R8 nowTime = t; 

      // Thickness tendency
      sw_tend_thick(0, 0, nowTime, Env, Halo,Mesh,State);

      // State adavance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel);
                                                   //+ 0.5* dt * LayerThicknessTend(ICell,KLevel);
         });

      // Normal Velocity tendency
      sw_tend_vel(1, 0, nowTime, Env, Halo,Mesh,State);

      // State adavance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel);
                                                   //+ 0.5* dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 2 --------------------------------------------------------------/
      //nowTime = t+0.5*dt; 
      nowTime = t+dt; 

      // Thickness tendency
      sw_tend_thick(1, 1, nowTime, Env, Halo,Mesh,State);

      // State advance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            //State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
            //                                       + dt * LayerThicknessTend(ICell,KLevel);
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
               + 0.5*(State->LayerThickness[1](ICell,KLevel) - State->LayerThickness[0](ICell,KLevel)
                      + dt * LayerThicknessTend(ICell,KLevel));
         });

      // Normal Velocity tendency
      sw_tend_vel(1, 1, nowTime, Env, Halo,Mesh,State);

      // State advance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            //State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
            //          + dt * NormalVelocityTend(IEdge,KLevel);
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
               + 0.5*(State->NormalVelocity[1](IEdge,KLevel) - State->NormalVelocity[0](IEdge,KLevel)
                      + dt * NormalVelocityTend(IEdge,KLevel));
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();

   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "midpoint" ) {
   //===--------------------------------------------------------------------===/

      // Stage 1 --------------------------------------------------------------/
      R8 nowTime = t; 

      // Thickness tendency
      sw_tend_thick(0, 0, nowTime, Env, Halo,Mesh,State);

      // Normal Velocity tendency
      sw_tend_vel(0, 0, nowTime, Env, Halo,Mesh,State);

      // State adavance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + 0.5* dt * LayerThicknessTend(ICell,KLevel);
         });

      // State adavance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + 0.5* dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 2 --------------------------------------------------------------/
      nowTime = t+0.5*dt; 

      // Thickness tendency
      sw_tend_thick(1, 1, nowTime, Env, Halo,Mesh,State);

      // Normal Velocity tendency
      sw_tend_vel(1, 1, nowTime, Env, Halo,Mesh,State);

      // State advance: LayerThickness
      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel);
         });

      // State advance: NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                      + dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();


   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "ssp-rk3" ) {
   //===--------------------------------------------------------------------===/

      // Stage 1 --------------------------------------------------------------/
      R8 nowTime = t; 

      // Tendencies
      sw_tend_thick(0, 0, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  0, 0, nowTime, Env, Halo, Mesh, State);

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
      nowTime = t + dt;

      // Tendencies
      sw_tend_thick(1, 1, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, nowTime, Env, Halo, Mesh, State);

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

      // Update halo ----------------------------------------------------------/
      State->copyToHost(1);
      Halo->exchangeFullArrayHalo(State->NormalVelocityH[1], OMEGA::OnEdge);
      Halo->exchangeFullArrayHalo(State->LayerThicknessH[1], OMEGA::OnCell);
      State->copyToDevice(1);


      // Stage 3 --------------------------------------------------------------/
      nowTime = t + 0.5*dt;

      // Tendencies
      sw_tend_thick(1, 1, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, nowTime, Env, Halo, Mesh, State);

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

   //===--------------------------------------------------------------------===/
   } else if ( time_integrator == "rk4" ) {
   //===--------------------------------------------------------------------===/

      // Stage 1 --------------------------------------------------------------/
      R8 nowTime = t; 

      // Tendencies
      sw_tend_thick(0, 0, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  0, 0, nowTime, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            LayerThicknessRKTemp(ICell,KLevel) = LayerThicknessTend(ICell,KLevel) / 6.0;
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + 0.5 * dt * LayerThicknessTend(ICell,KLevel);
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            NormalVelocityRKTemp(IEdge,KLevel) = NormalVelocityTend(IEdge,KLevel) / 6.0;
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + 0.5 * dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 2 --------------------------------------------------------------/
      nowTime = t + 0.5*dt;

      // Tendencies
      sw_tend_thick(1, 1, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, nowTime, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            LayerThicknessRKTemp(ICell,KLevel) += LayerThicknessTend(ICell,KLevel) / 3.0;
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + 0.5 * dt * LayerThicknessTend(ICell,KLevel);
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            NormalVelocityRKTemp(IEdge,KLevel) += NormalVelocityTend(IEdge,KLevel) / 3.0;
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                  + 0.5 * dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Update halo ----------------------------------------------------------/
      State->copyToHost(1);
      Halo->exchangeFullArrayHalo(State->NormalVelocityH[1], OMEGA::OnEdge);
      Halo->exchangeFullArrayHalo(State->LayerThicknessH[1], OMEGA::OnCell);
      State->copyToDevice(1);

      // Stage 3 --------------------------------------------------------------/
      nowTime = t + 0.5*dt;

      // Tendencies
      sw_tend_thick(1, 1, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, nowTime, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsAll,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            LayerThicknessRKTemp(ICell,KLevel) += LayerThicknessTend(ICell,KLevel) / 3.0;
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessTend(ICell,KLevel);
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            NormalVelocityRKTemp(IEdge,KLevel) += NormalVelocityTend(IEdge,KLevel) / 3.0;
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityTend(IEdge,KLevel);
         });

      // Stage 4 --------------------------------------------------------------/
      nowTime = t + dt;

      // Tendencies
      sw_tend_thick(1, 1, nowTime, Env, Halo, Mesh, State);
      sw_tend_vel(  1, 1, nowTime, Env, Halo, Mesh, State);

      // LayerThickness
      parallelFor(
         {Mesh->NCellsOwned,   NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            LayerThicknessRKTemp(ICell,KLevel) += LayerThicknessTend(ICell,KLevel) / 6.0;
            State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel)
                                                   + dt * LayerThicknessRKTemp(ICell,KLevel);
         });
      // NormalVelocity
      parallelFor(
         {Mesh->NEdgesOwned,   NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            NormalVelocityRKTemp(IEdge,KLevel) += NormalVelocityTend(IEdge,KLevel) / 6.0;
            State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                                   + dt * NormalVelocityRKTemp(IEdge,KLevel);
         });

      // Update time level ----------------------------------------------------/
      State->updateTimeLevels();

   //--------------------------------------------------------------------------/
   } else {
      LOG_ERROR("Invalid choice of time_integrator \
                (Choices: forward-euler, forward-backward, heuns, ssp-rk3)");
   } // if time_integrator

} // sw_time_stepper

} // namespace OMEGA

//===-----------------------------------------------------------------------===/
#include "ShallowWaterCore.h"
//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/
//=== Tendency for layer thickness -----------------------------------------===/
//===-----------------------------------------------------------------------===/
void ShallowWaterCore::
sw_tend_thick(int ThickCurTimeLevel, int VelCurTimeLevel, const MachEnv *Env,
                   const Halo *Halo, const HorzMesh *Mesh, OceanState *State) {

   Config *TendConfig;

   // Initialize Tendency
   parallelFor(
      {Mesh->NCellsSize, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         LayerThicknessTend(ICell,KLevel) = 0.0_Real;
      });

   // LayerThickEdge
   LayerThicknessAuxVars LayerThicknessAux(Mesh, NVertLevels);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         LayerThicknessAux.computeVarsOnEdge(IEdge, KLevel,
            State->LayerThickness[ThickCurTimeLevel], State->NormalVelocity[VelCurTimeLevel]);
      });
   LayerThicknessEdge = LayerThicknessAux.MeanLayerThickEdge;

   //if ( !ThickTend )
   // return;

   // NormalVelocity * LayerThickEdge
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         ThicknessFlux(IEdge,KLevel) = State->NormalVelocity[VelCurTimeLevel](IEdge,KLevel)
                                     * LayerThicknessAux.MeanLayerThickEdge(IEdge,KLevel);
                                     //* LayerThicknessAux.FluxLayerThickEdge(IEdge,KLevel);
      });

   ThicknessFluxDivOnCell ThickFluxDivOnC(Mesh, TendConfig);
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         ThickFluxDivOnC(LayerThicknessTend, ICell, KLevel, ThicknessFlux);
      });
} // sw_tend_thick


//===-----------------------------------------------------------------------===/
//=== Tendency for normal velocity -----------------------------------------===/
//===-----------------------------------------------------------------------===/
void ShallowWaterCore::
sw_tend_vel(int ThickCurTimeLevel, int VelCurTimeLevel, const MachEnv *Env,
            const Halo *Halo, const HorzMesh *Mesh, OceanState *State) {

   Config *TendConfig;

   // Initialize Tendency
   parallelFor(
      {Mesh->NEdgesSize, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         NormalVelocityTend(IEdge, KLevel) = 0.0_Real;
      });

   //if ( !VelTend )
   //return ;

   // AurxVars ----------------------------------------------------------------/

   // VorticityAuxVars - RelativeVorticityVertex, and Normalized by layerThickness
   VorticityAuxVars VorticityAux(Mesh, NVertLevels);
   parallelFor(
      {Mesh->NVerticesAll, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int KLevel) {
         VorticityAux.computeVarsOnVertex(IVertex, KLevel,
            State->LayerThickness[ThickCurTimeLevel], State->NormalVelocity[VelCurTimeLevel]);
      });

   // Vorticities : Vertex to edge
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         VorticityAux.computeVarsOnEdge(IEdge, KLevel);
      });

   // KineticAuxVars
   KineticAuxVars KineticAux(Mesh, NVertLevels);
   parallelFor(
       {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          KineticAux.computeVarsOnCell(ICell, KLevel, State->NormalVelocity[VelCurTimeLevel]);
       });

   // Update total depth
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         TotalDepthKECell(ICell,KLevel) =
              (  gravity * State->LayerThickness[ThickCurTimeLevel](ICell,KLevel)
               + BottomTopography(ICell) ) + KineticAux.KineticEnergyCell(ICell,KLevel);
      });

   // TendencyTerms -----------------------------------------------------------/

   // Potential Vorticity * h * uPerp
   PotentialVortHAdvOnEdge PotVortHAdvOnE(Mesh,TendConfig);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         PotVortHAdvOnE(NormalVelocityTend, IEdge, KLevel,
                        VorticityAux.NormRelVortEdge,VorticityAux.NormPlanetVortEdge,
                        LayerThicknessEdge,State->NormalVelocity[VelCurTimeLevel]);
      });

   // SSHGradOnEdge  --> Computation merged in KEGradOnEdge (TotalDepthKECell)
   //SSHGradOnEdge SSHGradOnE(Mesh, TendConfig);
   //parallelFor(
   //   {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
   //      SSHGradOnE(NormalVelocityTend, IEdge, KLevel, State->LayerThickness[ThickCurTimeLevel]);
   //      //SSHGradOnE(NormalVelocityTend, IEdge, KLevel, TotalDepth);
   //   });

   // KEGradOnEdge
   KEGradOnEdge KEGradOnE(Mesh, TendConfig);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         KEGradOnE(NormalVelocityTend, IEdge, KLevel, TotalDepthKECell);
      });
   // Del2 velocity
   //VelocityDiffusionOnEdge VelDiffOnE(Mesh, TendConfig);
   //parallelFor(
   //   {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
   //      VelDiffOnE(NormalVelocityTend, IEdge, KLevel,
   //                 KineticAux.VelocityDivCell, VorticityAux.RelVortVertex);
   //   });

} // sw_tend_vel

} // namespace OMEGA

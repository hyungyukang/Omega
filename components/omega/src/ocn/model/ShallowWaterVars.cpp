//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_init_var(HorzMesh *Mesh, const OceanState *State) {

   //--------------------------------------------------------------------------/

   // SW constants

   NVertLevels    = State->NVertLevels;

   NEdgesSize     = Mesh->NEdgesSize;
   NCellsSize     = Mesh->NCellsSize;
   NVerticesSize  = Mesh->NVerticesSize;
   NEdgesAll      = Mesh->NEdgesAll;
   NCellsAll      = Mesh->NCellsAll;
   NVerticesAll   = Mesh->NVerticesAll;
   NEdgesOwned    = Mesh->NEdgesOwned;
   NCellsOwned    = Mesh->NCellsOwned;
   NVerticesOwned = Mesh->NVerticesOwned;

   //--------------------------------------------------------------------------/

   // SW arrays

   Array2DReal normalVelocitySolution("normaVelocitySolution",NEdgesSize,NVertLevels);
   NormalVelocitySolution = normalVelocitySolution;

   Array2DReal layerThicknessSolution("layerThicknessSolution", NCellsSize, NVertLevels);
   LayerThicknessSolution = layerThicknessSolution;

   Array2DReal totalDepthKECell("totalDepthKECell", NCellsSize, NVertLevels);
   TotalDepthKECell = totalDepthKECell;

   Array2DReal layerThicknessEdge("layerThicknessEdge",NEdgesSize,NVertLevels);
   LayerThicknessEdge = layerThicknessEdge;

   Array2DReal thicknessFlux("thicknessFlux",NEdgesSize,NVertLevels);
   ThicknessFlux = thicknessFlux;

   Array2DReal tangentialVelocity("tangentialVelocity",NEdgesSize,NVertLevels);
   TangentialVelocity = tangentialVelocity;

   Array2DReal normalVelocityTend("normalVelocityTend",NEdgesSize,NVertLevels);
   NormalVelocityTend = normalVelocityTend;

   Array2DReal layerThicknessTend("layerThicknessTend", NCellsSize, NVertLevels);
   LayerThicknessTend = layerThicknessTend;

   Array2DReal normalVelocityRKTemp("normalVelocityRKTemp",NEdgesSize,NVertLevels);
   NormalVelocityRKTemp = normalVelocityRKTemp;

   Array2DReal layerThicknessRKTemp("layerThicknessRKTemp", NCellsSize, NVertLevels);
   LayerThicknessRKTemp = layerThicknessRKTemp;

   //-----------------------------/ 

   Array1DReal bottomTopography("bottomTopography", NCellsSize);
   BottomTopography = bottomTopography;

   //--------------------------------------------------------------------------/

} // sw_init_var

} // namespace OMEGA

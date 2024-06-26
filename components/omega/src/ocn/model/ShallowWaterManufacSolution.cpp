//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
#include <complex>
#define J dcomp(0.0,1.0)
typedef complex<double> dcomp;

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_manufactured_solution(const int TestCase, R8 t, 
                         const HorzMesh *Mesh, OceanState *State) {

   //--------------------------------------------------------------------------/
   // Select test case

   //const int TestCase = 2;
   //--------------------------------------------------------------------------/

   //--------------------------------------------------------------------------/
   // Test case 21: Solid body rotation
   //--------------------------------------------------------------------------/

   if ( TestCase == 21 ) {

      const R8 u0    = 2.0 * pii * Rearth / (12.0 * 86400);
      const R8 alpha = 0.0;
      const R8 k1 = 133681;
      const R8 omegaT = omega * t;

      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            const R8 sinAsinLat = sin(alpha) * sin(Mesh->LatEdgeH(IEdge));
            const R8 cosLonCosOT = cos(Mesh->LonEdgeH(IEdge)) * cos(omegaT);
            const R8 sinLonSinOT = sin(Mesh->LonEdgeH(IEdge)) * sin(omegaT);
            const R8 cosACosLat  = cos(alpha) * cos(Mesh->LatEdgeH(IEdge));

            uZonal(IEdge,KLevel) = u0 * ( sinAsinLat * (cosLonCosOT - sinLonSinOT)
                                         +cosACosLat) ;

            const R8 sinLonCosOT = sin(Mesh->LonEdgeH(IEdge)) * cos(omegaT);
            const R8 cosLonSinOT = cos(Mesh->LonEdgeH(IEdge)) * sin(omegaT);

            vMerid(IEdge,KLevel) = -u0 * sin(alpha) 
                                       * (sinLonCosOT + cosLonSinOT);
         });

      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            const R8 sinACosLat = sin(alpha) * cos(Mesh->LatCellH(ICell));
            const R8 cosLonCosOT = cos(Mesh->LonCellH(ICell)) * cos(omegaT);
            const R8 sinLonSinOT = sin(Mesh->LonCellH(ICell)) * sin(omegaT);
            const R8 cosASinLat  = cos(alpha) * sin(Mesh->LatCellH(ICell));
            const R8 ROSinLat    = Rearth * omega * sin(Mesh->LatCellH(ICell)); 

            const R8 htemp1 = u0 * ( sinACosLat * (-cosLonCosOT+sinLonSinOT)
                                    +cosASinLat ) + ROSinLat;

            const R8 htemp1s = pow(htemp1,2);
            const R8 ROSinLats    = pow(ROSinLat,2);
            
            const R8 gh = -0.5 * htemp1s + 0.5 * ROSinLats + k1;
            const R8 height = gh / gravity;

            if ( t == 0.0 ) {
               State->LayerThicknessH[0](ICell,KLevel) = height;
               State->LayerThicknessH[1](ICell,KLevel) = height;
            }

            LayerThicknessSolution(ICell,KLevel) = height;
          });

      // uZonal & vMerid -> NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            const R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge))*uZonal(IEdge,KLevel);
            const R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge))*vMerid(IEdge,KLevel);
            const R8 normalComp      = zonalNormalComp + meridNormalComp;

            if ( t == 0.0 ) {
               State->NormalVelocityH[0](IEdge,KLevel) = normalComp;
               State->NormalVelocityH[1](IEdge,KLevel) = normalComp;
            }

            NormalVelocitySolution(IEdge,KLevel) = normalComp;
         });

      // Bottom topography
      if ( t == 0 ) {
         parallelFor(
            {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
            const R8 sinLat2 = pow(sin(Mesh->LatCellH(ICell)),2);
            BottomTopography(ICell) = Rearth * omega * sinLat2 / (2.0*gravity);
         });
   
         State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
         State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
         State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
         State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
         BottomTopography         = createDeviceMirrorCopy(BottomTopography);
      }

   //--------------------------------------------------------------------------/
   // Use initial condition from the input file
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 0) {

   } else {
      LOG_ERROR("Invalid choice of TestCase with manufactured solutions \
                 (Choices: 21, 0)"); 
   } // if TestCase

} // sw_manufacturedSolution

} // namespace OMEGA

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
sw_init_field(const int TestCase, const MachEnv *Env, const Halo *Halo, 
              HorzMesh *Mesh, OceanState *State) {

   //--------------------------------------------------------------------------/
   // Test case 2: Global steady state nonlinear geostrophic flow
   //--------------------------------------------------------------------------/

   if ( TestCase == 2 ) {

      const R8 u0    = 2.0 * pii * Rearth / (12.0 * 86400);
      const R8 gh0   = 2.94e4;

      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      // Bottom topography
      parallelFor(
         {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
         BottomTopography(ICell) = 0.0;
      });

      // Zonal & Meridional velocities
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            uZonal(IEdge,KLevel) = u0 * cos(Mesh->LatEdgeH(IEdge));
            vMerid(IEdge,KLevel) = 0.0;
         });

      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            const R8 gh = gh0 - (Rearth*omega*u0 + u0*u0/2.0) * pow(sin(Mesh->LatCellH(ICell)),2);
            const R8 height = gh / gravity;

            State->LayerThicknessH[0](ICell,KLevel) = height;
            State->LayerThicknessH[1](ICell,KLevel) = height;

            LayerThicknessSolution(ICell,KLevel) = height;
          });

      // uZonal & vMerid -> NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            const R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge))*uZonal(IEdge,KLevel);
            const R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge))*vMerid(IEdge,KLevel);
            const R8 normalComp      = zonalNormalComp + meridNormalComp;

            // For some tests...
            //const R8 zonalTangentialComp = -sin(Mesh->AngleEdge(IEdge))*uZonal(IEdge,KLevel);
            //const R8 meridTangentialComp = cos(Mesh->AngleEdge(IEdge))*vMerid(IEdge,KLevel);
            //const R8 tangentialComp = zonalTangentialComp + meridTangentialComp;
            //const R8 zonalComp = normalComp * cos(Mesh->AngleEdge(IEdge))
            //                   - tangentialComp * sin(Mesh->AngleEdge(IEdge));
            //const R8 meridComp = 0.0;
            //TangentialVelocity(IEdge,KLevel) = tangentialComp;

            State->NormalVelocityH[0](IEdge,KLevel) = normalComp;
            State->NormalVelocityH[1](IEdge,KLevel) = normalComp;

            NormalVelocitySolution(IEdge,KLevel) = normalComp;
         });

      State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
      State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
      State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
      State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
      BottomTopography         = createDeviceMirrorCopy(BottomTopography);


   //--------------------------------------------------------------------------/
   // Test case 5: Zonal flow over an isolated mountain
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 5 ) {

      const R8 u0   = 20.0;
      const R8 gh0  = 5960.0 * gravity;
      const R8 hs0  = 2000.0;
      const R8 r0   =  pii / 9.0;
      const R8 clon = -pii / 2.0 + 2.0*pii;
      const R8 clat =  pii / 6.0;

      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      // Bottom topography
      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            const R8 sinClatSinLat        = sin(clat) * sin(Mesh->LatCellH(ICell));
            const R8 cosClatCosLatCosDlon = cos(clat) * cos(Mesh->LatCellH(ICell))
                                                      * cos(Mesh->LonCellH(ICell) - clon);
            const R8 dist = pow((Mesh->LonCellH(ICell)-clon),2)
                          + pow((Mesh->LatCellH(ICell)-clat),2);
            const R8 r0s = pow(r0,2);
            const R8 rsquare = min(r0s,dist);
            const R8 r = sqrt(rsquare);

            BottomTopography(ICell) = hs0 * (1.0-r/r0);
         });

      // Zonal & Meridional velocities
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            uZonal(IEdge,KLevel) = u0 * cos(Mesh->LatEdgeH(IEdge));
            vMerid(IEdge,KLevel) = 0.0;
         });

      parallelFor(
         {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
            const R8 gh = gh0 - (Rearth*omega*u0 + u0*u0/2.0) 
                        * pow(sin(Mesh->LatCellH(ICell)),2);
            const R8 height = gh / gravity;

            // Fluid thickness = Total depth - BottomTopography
            State->LayerThicknessH[0](ICell,KLevel) = height-BottomTopography(ICell);
            State->LayerThicknessH[1](ICell,KLevel) = height-BottomTopography(ICell);

            LayerThicknessSolution(ICell,KLevel)    = height-BottomTopography(ICell);
          });

      // uZonal & vMerid -> NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            const R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge))*uZonal(IEdge,KLevel);
            const R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge))*vMerid(IEdge,KLevel);
            const R8 normalComp      = zonalNormalComp + meridNormalComp;

            State->NormalVelocityH[0](IEdge,KLevel) = normalComp;
            State->NormalVelocityH[1](IEdge,KLevel) = normalComp;

            NormalVelocitySolution(IEdge,KLevel)    = normalComp;
         });

      State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
      State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
      State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
      State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
      BottomTopography         = createDeviceMirrorCopy(BottomTopography);

   //--------------------------------------------------------------------------/
   // Solid body rotation with time-dependent solution
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 21) {

       sw_time_dependent_solution(TestCase, "init", false, false, 0.0, Mesh, State);

   //--------------------------------------------------------------------------/
   // Manfactured solution with time-dependent solution
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 22) {

       sw_time_dependent_solution(TestCase, "init", true, true,  0.0, Mesh, State);

   //--------------------------------------------------------------------------/
   // Use initial condition from the input file
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 0) {

   //--------------------------------------------------------------------------/
   // Spherical harmonics: Y(m=3,n=4) - test purpose
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 99) {

      Array2DR8 streamf("streamf", Mesh->NCellsSize,NVertLevels);
      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      const int m = 3;

      // Spherical harmonics Y(m=4,n=3)
      parallelFor(
         {Mesh->NCellsAll,NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
             const R8 z = sin(Mesh->LatCellH(ICell));
             const R8 powZ = pow((1.0-pow(z,2)),3.0/2.0);
             const R8 p43  = -105.0 * z * powZ;

             auto c = exp(m*Mesh->LonCellH(ICell)*J);
             auto cc = c * p43;

             R8 height = real(cc);

             State->LayerThickness[0](ICell,KLevel) = height;
             State->LayerThickness[1](ICell,KLevel) = height;
             State->LayerThicknessH[0](ICell,KLevel) = height;
             State->LayerThicknessH[1](ICell,KLevel) = height;

          });
   } else {
      LOG_ERROR("Invalid choice of TestCase \
                 (Choices: 0, 2)"); 
   } // if TestCase

} // sw_init_field

} // namespace OMEGA

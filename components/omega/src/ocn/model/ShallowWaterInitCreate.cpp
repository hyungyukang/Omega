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
              const HorzMesh *Mesh, OceanState *State) {

   //--------------------------------------------------------------------------/
   // Select test case

   //const int TestCase = 2;
   //--------------------------------------------------------------------------/

   //--------------------------------------------------------------------------/
   // Test case 2: Global steady state nonlinear geostrophic flow
   //--------------------------------------------------------------------------/

   if ( TestCase == 2 ) {

      const R8 u0    = 2.0 * pii * Rearth / (12.0 * 86400);
      const R8 gh0   = 2.94e4;
      //const R8 om    = 7.292e-5;
      //const R8 alpha = 0.0;

      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

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

            //SWVar->LayerThicknessInit(ICell,KLevel) = height;
            LayerThicknessInit(ICell,KLevel) = height;
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

            //SWVar->NormalVelocityInit(IEdge,KLevel) = normalComp;
            NormalVelocityInit(IEdge,KLevel) = normalComp;
         });

      // Bottom topography
      parallelFor(
         {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
         BottomTopography(ICell) = 0.0;
      });

      State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
      State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
      State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
      State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
      BottomTopography         = createDeviceMirrorCopy(BottomTopography);


   //--------------------------------------------------------------------------/
   // Use initial condition from the input file
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 0) {

   //--------------------------------------------------------------------------/

   // Spherical harmonics: Y(m=3,n=4) - test purpose
   } else if ( TestCase == 99) {

      Array2DR8 streamf("streamf", Mesh->NCellsSize,NVertLevels);
      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      const int m = 3;
      //const int n = 4;

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

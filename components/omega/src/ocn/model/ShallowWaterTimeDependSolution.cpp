//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <unordered_map>

using namespace std;
#include <complex>
#define J dcomp(0.0,1.0)
typedef complex<double> dcomp;

void parseFile(const std::string& filename, std::unordered_map<std::string, std::string>& data) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    std::string label;
    std::string value;
    while (file >> label >> value) {
        data[label] = value;
    }
}
//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_time_dependent_solution(const int TestCase, const std::string &solutionOpt,
                           bool tendThickSourceTerm, bool tendVelSourceTerm,
                           const R8 nowTime, HorzMesh *Mesh, OceanState *State) {

   //--------------------------------------------------------------------------/
   // Test case 21: Solid body rotation with time-dependent solutions
   //--------------------------------------------------------------------------/

   if ( TestCase == 21 ) {

      const R8 u0    = 2.0 * pii * Rearth / (12.0 * 86400);
      const R8 alpha = 0.0;
      const R8 k1 = 133681;
      const R8 omegaT = omega * nowTime;

      Array2DR8 uZonal("uZonal", Mesh->NEdgesSize, NVertLevels);
      Array2DR8 vMerid("vMerid", Mesh->NEdgesSize, NVertLevels);

      // Bottom topography
      if ( nowTime == 0 ) {
         parallelFor(
            {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
            const R8 sinLat2 = pow(sin(Mesh->LatCellH(ICell)),2);
            BottomTopography(ICell) = Rearth * omega * sinLat2 / (2.0*gravity);
         });
      }

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

            // Thickness = Total depth - BottomTopo
            if ( nowTime == 0.0 ) {
               State->LayerThicknessH[0](ICell,KLevel) = height-BottomTopography(ICell);
               State->LayerThicknessH[1](ICell,KLevel) = height-BottomTopography(ICell);
            }

            LayerThicknessSolution(ICell,KLevel) = height-BottomTopography(ICell);
          });

      // uZonal & vMerid -> NormalVelocity
      parallelFor(
         {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
            const R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge))*uZonal(IEdge,KLevel);
            const R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge))*vMerid(IEdge,KLevel);
            const R8 normalComp      = zonalNormalComp + meridNormalComp;

            if ( nowTime == 0.0 ) {
               State->NormalVelocityH[0](IEdge,KLevel) = normalComp;
               State->NormalVelocityH[1](IEdge,KLevel) = normalComp;
            }

            NormalVelocitySolution(IEdge,KLevel) = normalComp;
         });

      if ( nowTime == 0 ) {
         State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
         State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
         State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
         State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
         BottomTopography         = createDeviceMirrorCopy(BottomTopography);
      }

   //--------------------------------------------------------------------------/
   // Test case 22: Nonlinear manufactured solution
   //--------------------------------------------------------------------------/

   } else if ( TestCase == 22 ) {

     // Set default values

     //const R8 H0 = 1000.0;
     H0 = 1000.0 ; // not const as it is defined as a global variable.
     const R8 etaHat1 = 1.00;
     const R8 f0 = 1.0e-4;
     const I4 mx = 2;
     const I4 my = 2;

     const R8 g = gravity;
     const R8 lX = 10000.0 * 1000.0;   // maybe from Mesh
     const R8 lY = lX*sqrt(3.0)/2.0; // maybe from Mesh
     const R8 kX1 = mx* 2.0*pii/lX;
     const R8 kY1 = my* 2.0*pii/lY;
     const R8 omega1 = sqrt(g*H0*(kX1*kX1 + kY1*kY1));
     const R8 omegaT = omega1 * nowTime;

     if ( solutionOpt == "init" ) {

        //R8 H0 = 0.0;
        //R8 etaHat1 = 0.0;
        //R8 f0 = 0.0;
        //I4 mx = 0;
        //I4 my = 0;

        //std::unordered_map<std::string, std::string> data;
        //parseFile("manufactured_solution.cfg", data);
    
        //// Convert and store values in appropriate types
        //const R8 H0 = std::stod(data["H0"]);
        //const R8 etaHat1 = std::stod(data["etaHat1"]);
        //const R8 f0 = std::stod(data["f0"]);
        //const I4 mx = std::stoi(data["mx"]);
        //const I4 my = std::stoi(data["my"]);
    
        //// Print the values to verify
        //std::cout << "H0: " << H0 << std::endl;
        //std::cout << "etaHat1: " << etaHat1 << std::endl;
        //std::cout << "f0: " << f0 << std::endl;
        //std::cout << "mx: " << mx << std::endl;
        //std::cout << "my: " << my << std::endl;

       // Bottom topography
       parallelFor(
          {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          BottomTopography(ICell) = 0.0;
       });

       // FVertex - Coriolis parameter
       parallelFor(
          {Mesh->NVerticesAll}, KOKKOS_LAMBDA(int IVertex) {
           Mesh->FVertexH(IVertex) = f0 ;
       });
     }


     if ( tendThickSourceTerm ) {
       // surface elevation
       parallelFor(
          {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
             const R8 x = Mesh->XCellH(ICell);
             const R8 y = Mesh->YCellH(ICell);

             const R8 phase = kX1*x + kY1*y - omegaT;
             const R8 eta = etaHat1 * sin(phase); 
             const R8 etaSourceTerm = ( etaHat1*(-H0*(kX1+kY1) * sin(phase)
                                       -omega1*cos(phase) + etaHat1*(kX1+kY1)*cos(2.0*phase)));

             if ( solutionOpt == "init" && !initFieldFromFile ) {
                State->LayerThicknessH[0](ICell,KLevel) = eta + H0;
                State->LayerThicknessH[1](ICell,KLevel) = eta + H0;
             } else if ( solutionOpt == "tend" ) {
              LayerThicknessTend(ICell,KLevel) += etaSourceTerm;
             } else if ( solutionOpt == "solution") { 
                LayerThicknessSolution(ICell,KLevel)  = eta + H0;
             }
          });
     }

     // Zonal & merid veolcity -> normal velocity
     if ( tendVelSourceTerm ) {
       parallelFor(
          {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
             const R8 x = Mesh->XEdgeH(IEdge);
             const R8 y = Mesh->YEdgeH(IEdge);

             const R8 phase = kX1*x + kY1*y - omegaT;
             const R8 u = etaHat1 * cos(phase);
             const R8 v = etaHat1 * cos(phase);
             const R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge)) * u;
             const R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge)) * v;
             const R8 normalComp      = zonalNormalComp + meridNormalComp;

             const R8 SourceTerm0 = omega1*sin(phase)
                                  - 0.5*etaHat1*(kX1 + kY1) * sin(2.0*phase);
             const R8 uSourceTerm = etaHat1*((-f0 + g*kX1) * cos(phase) + SourceTerm0);
             const R8 vSourceTerm = etaHat1*(( f0 + g*kY1) * cos(phase) + SourceTerm0);

             const R8 zonalNormalCompSourceTerm = cos(Mesh->AngleEdge(IEdge)) * uSourceTerm;
             const R8 meridNormalCompSourceTerm = sin(Mesh->AngleEdge(IEdge)) * vSourceTerm;
             const R8 normalCompSourceTerm      = zonalNormalCompSourceTerm + meridNormalCompSourceTerm;

             if ( solutionOpt == "init" && !initFieldFromFile ) {
                State->NormalVelocityH[0](IEdge,KLevel) = normalComp;
                State->NormalVelocityH[1](IEdge,KLevel) = normalComp;
             } else if ( solutionOpt == "tend" ) {
                NormalVelocityTend(IEdge,KLevel) += normalCompSourceTerm;
             } else if ( solutionOpt == "solution") { 
                NormalVelocitySolution(IEdge,KLevel)  = normalComp;
             }

          });
     }

     if ( solutionOpt == "init" ) {
        State->NormalVelocity[0] = createDeviceMirrorCopy(State->NormalVelocityH[0]);
        State->NormalVelocity[1] = createDeviceMirrorCopy(State->NormalVelocityH[1]);
        State->LayerThickness[0] = createDeviceMirrorCopy(State->LayerThicknessH[0]);
        State->LayerThickness[1] = createDeviceMirrorCopy(State->LayerThicknessH[1]);
        BottomTopography         = createDeviceMirrorCopy(BottomTopography);
        Mesh->FVertex            = createDeviceMirrorCopy(Mesh->FVertexH);
     }

   //--------------------------------------------------------------------------/
   // Use initial condition from the input file
   //--------------------------------------------------------------------------/
   } else if ( TestCase == 0) {

   } else {
      LOG_ERROR("Invalid choice of TestCase with time-dependent solutions \
                 (Choices: 21, 22, 0)"); 
   } // if TestCase

} // sw_time_dependent_solution

} // namespace OMEGA



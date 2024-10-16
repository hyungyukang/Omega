#ifndef OMEGA_MANUFACTUREDSOLUTION_H
#define OMEGA_MANUFACTUREDSOLUTION_H
//===-- ocn/share/ManufacturedSolution.h - Manufactured Solution ---*- C++ -*-===//
//
// Customized tendencies for computing the thickness and normal velocity tendencies
// for the manufactured solution test. See Bishnu et al. 2024
// (https://doi.org/10.1029/2022MS003545) and the manufactured solution test case
// in Polaris.
//
//===-------------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "TimeMgr.h"

namespace OMEGA {

//===--------------------------------------------------------------------------===/
// Manufactured solution for thickness tendency
//===--------------------------------------------------------------------------===/
struct ManufacturedThicknessTendency {

   R8 Amplitude;
   R8 RestingThickness;
   TimeFrac RefElapsedTime;

   void operator()(Array2DReal ThicknessTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) {

     // Get ElapsedTime
     TimeFrac ElapsedTime = Time.getElapsedTime() - RefElapsedTime;
     R8 ElapsedTimeSec = ElapsedTime.getSeconds();

     LOG_INFO("InMS ElapsedTime {}", ElapsedTimeSec);

     R8 etaHat1 = Amplitude; // Amplitude
     R8 H0 = RestingThickness; // Depth

     R8 f0 = 1.0e-4; // Coriolis parameter
     I4 mx = 2;      // X-dir wavenumber
     I4 my = 2;      // Y-dir wavenumber

     R8 pii = 3.141592653589793;  // PI
     R8 g = 9.80665_Real;         // Gravity acceleration
     R8 lX = 10000.0 * 1000.0;    // X-dir domain extent
     R8 lY = lX*sqrt(3.0)/2.0;    // Y-dir domain extent
     R8 kX1 = mx * 2.0*pii/lX;    // Wave in X-dir
     R8 kY1 = my * 2.0*pii/lY;    // Wave in Y-dir

     R8 omega1 = sqrt(g*H0*(kX1*kX1 + kY1*kY1));
     R8 omegaT = omega1 * ElapsedTimeSec;

     auto *Mesh                = HorzMesh::getDefault();
     auto NVertLevels          = ThicknessTend.extent_int(1);

     parallelFor(
        {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
           R8 x = Mesh->XCellH(ICell);
           R8 y = Mesh->YCellH(ICell);

           R8 phase = kX1*x + kY1*y - omegaT;
           R8 eta = etaHat1 * sin(phase);
           R8 etaSourceTerm = ( etaHat1*(-H0*(kX1+kY1) * sin(phase)
                                     -omega1*cos(phase) + etaHat1*(kX1+kY1)*cos(2.0*phase)));

           ThicknessTend(ICell,KLevel) += etaSourceTerm;
        });
   }
}; // end struct ManufacturedThicknessTendency

//===--------------------------------------------------------------------------===/
// Manufactured solution for normal velocity tendency
//===--------------------------------------------------------------------------===/
struct ManufacturedVelocityTendency {

   R8 Amplitude;
   R8 RestingThickness;
   TimeFrac RefElapsedTime;

   void operator()(Array2DReal NormalVelTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) {

     // Get ElapsedTime
     TimeFrac ElapsedTime = Time.getElapsedTime() - RefElapsedTime;
     R8 ElapsedTimeSec = ElapsedTime.getSeconds();
     LOG_INFO("InMS ElapsedTime");

     R8 etaHat1 = Amplitude; // Amplitude
     R8 H0 = RestingThickness;    // Depth
     R8 f0 = 1.0e-4; // Coriolis parameter
     I4 mx = 2;      // X-dir wavenumber
     I4 my = 2;      // Y-dir wavenumber

     R8 pii = 3.141592653589793;  // PI
     R8 g = 9.80665_Real;         // Gravity acceleration
     R8 lX = 10000.0 * 1000.0;    // X-dir domain extent
     R8 lY = lX*sqrt(3.0)/2.0;    // Y-dir domain extent
     R8 kX1 = mx * 2.0*pii/lX;    // Wave in X-dir
     R8 kY1 = my * 2.0*pii/lY;    // Wave in Y-dir

     R8 omega1 = sqrt(g*H0*(kX1*kX1 + kY1*kY1));
     R8 omegaT = omega1 * ElapsedTimeSec;

     auto *Mesh                = HorzMesh::getDefault();
     auto NVertLevels          = NormalVelTend.extent_int(1);

     parallelFor(
        {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
           R8 x = Mesh->XEdgeH(IEdge);
           R8 y = Mesh->YEdgeH(IEdge);

           R8 phase = kX1*x + kY1*y - omegaT;
           R8 u = etaHat1 * cos(phase);
           R8 v = etaHat1 * cos(phase);
           R8 zonalNormalComp = cos(Mesh->AngleEdge(IEdge)) * u;
           R8 meridNormalComp = sin(Mesh->AngleEdge(IEdge)) * v;
           R8 normalComp      = zonalNormalComp + meridNormalComp;

           R8 SourceTerm0 = omega1*sin(phase)
                          - 0.5*etaHat1*(kX1 + kY1) * sin(2.0*phase);
           R8 uSourceTerm = etaHat1*((-f0 + g*kX1) * cos(phase) + SourceTerm0);
           R8 vSourceTerm = etaHat1*(( f0 + g*kY1) * cos(phase) + SourceTerm0);

           R8 zonalNormalCompSourceTerm = cos(Mesh->AngleEdge(IEdge)) * uSourceTerm;
           R8 meridNormalCompSourceTerm = sin(Mesh->AngleEdge(IEdge)) * vSourceTerm;
           R8 normalCompSourceTerm      = zonalNormalCompSourceTerm + meridNormalCompSourceTerm;

           NormalVelTend(IEdge,KLevel) += normalCompSourceTerm;

        });
   }
}; // end struct ManufacturedVelocityTendency

//===--------------------------------------------------------------------------===/
// A class for the manufactured tendencies
//===--------------------------------------------------------------------------===/
class ManufacturedSolution {

   public:
      // Instances of manufactured tendencies
      ManufacturedThicknessTendency ManufacturedThickTend;
      ManufacturedVelocityTendency ManufacturedVelTend;

      int init();
      //Calendar MSCal;

};
//===--------------------------------------------------------------------------===/

} // namespace OMEGA
#endif

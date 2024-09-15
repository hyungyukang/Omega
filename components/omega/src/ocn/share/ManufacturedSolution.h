#ifndef OMEGA_MANUFACTUREDSOLUTION_H
#define OMEGA_MANUFACTUREDSOLUTION_H

#include "OceanConstans.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"

//===-----------------------------------------------------------------------===/

using namespace OMEGA;
//===-----------------------------------------------------------------------===/

// Manufactured solution for thickness tendency
struct ManufacturedThicknessTendency {
   void operator()(Array2DReal ThicknessTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) const {

     // Get ElapsedTime
     R8 ElapsedTimeSec;
     I4 Err = StartTime.get(ElapsedTimeSec, TimeUnits::Seconds);


     R8 H0 = 1000.0;
     R8 etaHat1 = 1.00;
     I4 mx = 2;
     I4 my = 2;

     R8 g = Grav;
     R8 lX = 10000.0 * 1000.0;   // maybe from Mesh
     R8 lY = lX*sqrt(3.0)/2.0; // maybe from Mesh
     R8 kX1 = mx* 2.0*pii/lX;
     R8 kY1 = my* 2.0*pii/lY;
     R8 omega1 = sqrt(g*H0*(kX1*kX1 + kY1*kY1));
     R8 omegaT = omega1 * ElapsedTimeSec;

     auto *Mesh                = HorzMesh::getDefault();
     auto NVertLevels          = ThicknessTend.extent_int(1);
     //const auto &NormalVelEdge = State->NormalVelocity[VelTimeLevel];

     parallelFor(
        {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
           const R8 x = Mesh->XCellH(ICell);
           const R8 y = Mesh->YCellH(ICell);

           const R8 phase = kX1*x + kY1*y - omegaT;
           const R8 eta = etaHat1 * sin(phase);
           const R8 etaSourceTerm = ( etaHat1*(-H0*(kX1+kY1) * sin(phase)
                                     -omega1*cos(phase) + etaHat1*(kX1+kY1)*cos(2.0*phase)));

           ThicknessTend(ICell,KLevel) += etaSourceTerm;
        });
   }
};

// Manufactured solution for normal velocity tendency
struct ManufacturedVelocityTendency {
   void operator()(Array2DReal NormalVelTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) const {

     // Get ElapsedTime
     R8 ElapsedTimeSec;
     I4 Err = StartTime.get(ElapsedTimeSec, TimeUnits::Seconds);


     R8 H0 = 1000.0;
     R8 etaHat1 = 1.00;
     R8 f0 = 1.0e-4;
     I4 mx = 2;
     I4 my = 2;

     R8 g = Grav;
     R8 lX = 10000.0 * 1000.0;   // maybe from Mesh
     R8 lY = lX*sqrt(3.0)/2.0; // maybe from Mesh
     R8 kX1 = mx* 2.0*pii/lX;
     R8 kY1 = my* 2.0*pii/lY;
     R8 omega1 = sqrt(g*H0*(kX1*kX1 + kY1*kY1));
     R8 omegaT = omega1 * ElapsedTimeSec;

     auto *Mesh                = HorzMesh::getDefault();
     auto NVertLevels          = NormalVelTend.extent_int(1);

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

           NormalVelTend(IEdge,KLevel) += normalCompSourceTerm;

        });
   }
};

//===-----------------------------------------------------------------------===/
#endif

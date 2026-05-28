//===-- SplitExplicitBarotropicPCStepper.cpp - SE stage 2 -----*- C++ -*-===//
//
// Framework for the explicitly subcycled forward-backward
// predictor-corrector barotropic velocity update.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "SplitExplicitInit.h"

#include <utility>

namespace OMEGA {

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::doSplitStage2(
    OceanState *State, SplitExplicitScratch &Scratch,
    const SplitExplicitConfig &SEConfig, const HorzMesh *Mesh,
    const VertCoord *VCoord, I4 TimeLevel, const TimeInstant &StageTime,
    const TimeInterval &StageTimeStep) const {

   if (SEConfig.NBtrSubcycles < 1)
      LOG_CRITICAL("Invalid split-explicit barotropic subcycle count");

   // Keep the arguments live in this framework implementation. The barotropic
   // tendency terms will use them once the full split-explicit equations land.
   (void)StageTime;

   Eos *EqState = Eos::getInstance();

   R8 StageTimeStepSeconds;
   StageTimeStep.get(StageTimeStepSeconds, TimeUnits::Seconds);

   const I4 NBtrSubcycles       = 2 * SEConfig.NBtrSubcycles;
   const Real BtrDt             = StageTimeStepSeconds / SEConfig.NBtrSubcycles;
   const Real InvBtrVelAvgCount = 1._Real / (NBtrSubcycles + 1);
   const Real InvBtrFluxAvgCount = 1._Real / NBtrSubcycles;
   constexpr Real Gamma1         = 0.5333_Real;
   constexpr Real Gamma2         = 0.5333_Real;
   constexpr Real Gamma3         = 1.0_Real;
   constexpr Real RhoGravity     = RhoSw * Gravity;

   Pacer::start("SE-RK2:stage2BtrPC", 2);

   // Get array required
   Array1DReal NormalBtrVel    = State->getNormalBarotropicVelocity(TimeLevel);
   Array1DReal BtrPressAnomaly = State->getBarotropicPressureAnomaly(TimeLevel);
   Array1DReal NormalBtrVelSubcycleCur =
       Scratch.NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBtrVelSubcycleNew =
       Scratch.NormalBarotropicVelocitySubcycleNew;
   Array1DReal BtrPressAnomalySubcycleCur =
       Scratch.BarotropicPressureAnomalySubcycleCur;
   Array1DReal BtrPressAnomalySubcycleNew =
       Scratch.BarotropicPressureAnomalySubcycleNew;
   Array1DReal BtrPressure = Scratch.BarotropicPressure;
   Array1DReal BtrForcing  = Scratch.BarotropicForcing;
   Array1DReal BtrFlux     = Scratch.BarotropicFlux;

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(EdgeSignOnCell, Mesh->EdgeSignOnCell);
   OMEGA_SCOPE(NEdgesOnEdge, Mesh->NEdgesOnEdge);
   OMEGA_SCOPE(EdgesOnEdge, Mesh->EdgesOnEdge);
   OMEGA_SCOPE(WeightsOnEdge, Mesh->WeightsOnEdge);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(FEdge, Mesh->FEdge);
   OMEGA_SCOPE(BottomDepth, VCoord->BottomDepth);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgetop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(DepthMeanSpecVol, EqState->DepthMeanSpecificVolume);

   // Initilize barotropic vars
   deepCopy(NormalBtrVelSubcycleCur, NormalBtrVel);
   deepCopy(BtrPressAnomalySubcycleCur, BtrPressAnomaly);
   deepCopy(BtrFlux, 0._Real);

   for (I4 Subcycle = 0; Subcycle < NBtrSubcycles; ++Subcycle) {

      // Barotropic velocity predictor
      parallelFor(
          "btrVelocityPredictor", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real Mask = EdgeMask(IEdge, 0);

             Real CoriolisTend = 0._Real;
             for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
                const I4 JEdge = EdgesOnEdge(IEdge, J);
                CoriolisTend += WeightsOnEdge(IEdge, J) *
                                NormalBtrVelSubcycleCur(JEdge) * FEdge(JEdge);
             }

             const Real MeanSpecVolEdge =
                 0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));

             const Real BtrPressAnomalyGrad =
                 (BtrPressAnomalySubcycleCur(Cell1) -
                  BtrPressAnomalySubcycleCur(Cell0)) /
                 DcEdge(IEdge);

             NormalBtrVelSubcycleNew(IEdge) =
                 Mask * (NormalBtrVelSubcycleCur(IEdge) +
                         BtrDt * (CoriolisTend -
                                  MeanSpecVolEdge * BtrPressAnomalyGrad +
                                  BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly predictor
      parallelFor(
          "btrPressurePredictor", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);
                const Real BtrPressureEdge =
                    0.5_Real * (BtrPressAnomalySubcycleCur(Cell0) +
                                RhoGravity * BottomDepth(Cell0) +
                                BtrPressAnomalySubcycleCur(Cell1) +
                                RhoGravity * BottomDepth(Cell1));
                const Real NormalBtrVelEdge =
                    (1._Real - Gamma1) * NormalBtrVelSubcycleCur(JEdge) +
                    Gamma1 * NormalBtrVelSubcycleNew(JEdge);
                const Real BtrFlux = BtrPressureEdge * NormalBtrVelEdge;

                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * BtrFlux;
             }

             // TODO: Surface freshwater flux (rho0*g*Q) is assumed to be 0 for
             // now.
             BtrPressAnomalySubcycleNew(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });

      // Barotropic velocity corrector
      parallelFor(
          "btrVelocityCorrector", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(I4 IEdge) {
             const Real Mask = EdgeMask(IEdge, 0);
             const I4 Cell0  = CellsOnEdge(IEdge, 0);
             const I4 Cell1  = CellsOnEdge(IEdge, 1);

             Real CoriolisTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
                const I4 JEdge = EdgesOnEdge(IEdge, J);
                CoriolisTend += WeightsOnEdge(IEdge, J) *
                                NormalBtrVelSubcycleNew(JEdge) * FEdge(JEdge);
             }

             const Real BtrPressure0 =
                 (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell0) +
                 Gamma2 * BtrPressAnomalySubcycleNew(Cell0);

             const Real BtrPressure1 =
                 (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell1) +
                 Gamma2 * BtrPressAnomalySubcycleNew(Cell1);

             const Real MeanSpecVolEdge =
                 0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));

             const Real BtrPressureGrad =
                 (BtrPressure1 - BtrPressure0) / DcEdge(IEdge);

             NormalBtrVelSubcycleNew(IEdge) =
                 Mask *
                 (NormalBtrVelSubcycleCur(IEdge) +
                  BtrDt * (CoriolisTend - MeanSpecVolEdge * BtrPressureGrad +
                           BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly corrector
      parallelFor(
          "btrPressureCorrector", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);

                const Real BtrPressure0 =
                    (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell0) +
                    Gamma2 * BtrPressAnomalySubcycleNew(Cell0);

                const Real BtrPressure1 =
                    (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell1) +
                    Gamma2 * BtrPressAnomalySubcycleNew(Cell1);

                const Real BtrPressureEdge =
                    0.5_Real * (BtrPressure0 + RhoGravity * BottomDepth(Cell0) +
                                BtrPressure1 + RhoGravity * BottomDepth(Cell1));

                const Real NormalBtrVelEdge =
                    (1._Real - Gamma3) * NormalBtrVelSubcycleCur(JEdge) +
                    Gamma3 * NormalBtrVelSubcycleNew(JEdge);

                const Real BtrFlux = BtrPressureEdge * NormalBtrVelEdge;
                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * BtrFlux;
             }

             BtrPressAnomalySubcycleNew(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });
   } // Subcycle

   deepCopy(NormalBtrVel, Scratch.NormalBarotropicVelocitySubcycleCur);
   SplitExplicitInit::combineVelocitySplit(State, Mesh, VCoord, TimeLevel);

   Pacer::stop("SE-RK2:stage2BtrPC", 2);
}

} // namespace OMEGA

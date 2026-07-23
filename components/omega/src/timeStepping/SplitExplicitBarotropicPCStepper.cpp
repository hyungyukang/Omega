//===-- SplitExplicitBarotropicPCStepper.cpp - SE stage 2 -----*- C++ -*-===//
//
// Framework for the explicitly subcycled forward-backward
// predictor-corrector barotropic velocity update.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "AuxiliaryState.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"

#include <utility>

namespace OMEGA {

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::doSplitStage2(
    OceanState *State, const AuxiliaryState *AuxState,
    SplitExplicitScratch &Scratch,
    const SplitExplicitConfig &SEConfig, const HorzMesh *Mesh,
    Halo *MeshHalo, const VertCoord *VCoord, I4 CurLevel, I4 NextLevel,
    const TimeInstant &StageTime, const TimeInterval &StageTimeStep) const {

   if (SEConfig.NBtrSubcycles < 1)
      LOG_CRITICAL("Invalid split-explicit barotropic subcycle count");
   if (!MeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");
   if (!AuxState)
      LOG_CRITICAL("Invalid auxiliary state");

   // Keep the arguments live in this framework implementation. The barotropic
   // tendency terms will use them once the full split-explicit equations land.
   (void)StageTime;

   Eos *EqState = Eos::getInstance();

   R8 StageTimeStepSeconds;
   StageTimeStep.get(StageTimeStepSeconds, TimeUnits::Seconds);

   const I4 NBtrSubcycles        = 2 * SEConfig.NBtrSubcycles;
   const Real BtrDt              = StageTimeStepSeconds / SEConfig.NBtrSubcycles;
   const Real InvBtrVelAvgCount  = 1._Real / (NBtrSubcycles + 1);
   const Real InvBtrFluxAvgCount = 1._Real / NBtrSubcycles;
   constexpr Real Gamma1         = 0.5333_Real;
   constexpr Real Gamma2         = 0.5333_Real;
   constexpr Real Gamma3         = 1.0_Real;
   constexpr Real RhoGravity     = RhoSw * Gravity;

   Pacer::start("SE-RK2:stage2BtrPC", 2);

   // Get array required
   Array1DReal NormalBtrVelCur     = State->getNormalBarotropicVelocity(CurLevel);
   Array1DReal BtrPressAnomalyCur  = State->getBarotropicPressureAnomaly(CurLevel);
   Array1DReal NormalBtrVelNext    = State->getNormalBarotropicVelocity(NextLevel);
   Array1DReal BtrPressAnomalyNext = State->getBarotropicPressureAnomaly(NextLevel);

   Array1DReal NormalBtrVelSubcycleCur =
       Scratch.NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBtrVelSubcycleNew =
       Scratch.NormalBarotropicVelocitySubcycleNew;
   Array1DReal BtrPressAnomalySubcycleCur =
       Scratch.BarotropicPressureAnomalySubcycleCur;
   Array1DReal BtrPressAnomalySubcycleNew =
       Scratch.BarotropicPressureAnomalySubcycleNew;
   Array1DReal BtrForcing = Scratch.BarotropicForcing;
   Array1DReal BtrFlux    = Scratch.BarotropicFlux;
   Array1DReal BaroclinicPseudoThicknessEdge =
       Scratch.BaroclinicPseudoThicknessEdge;
   const Array2DReal FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   const bool UpwindFluxThickness =
       AuxState->PseudoThicknessAux.FluxThickEdgeChoice ==
       FluxThickEdgeOption::Upwind;

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
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(DepthMeanSpecVol, EqState->DepthMeanSpecificVolume);

   // Initialize barotropic vars
   deepCopy(NormalBtrVelSubcycleCur, NormalBtrVelCur);
   deepCopy(BtrPressAnomalySubcycleCur, BtrPressAnomalyCur);
   deepCopy(NormalBtrVelNext, NormalBtrVelSubcycleCur);
   deepCopy(BtrFlux, 0._Real);

   // MPAS-O's baroclinicThickness is the vertical sum of the provisional
   // layer thickness used for fluxes. Reuse Omega's corresponding edge
   // pseudo-thickness so that partial bottom cells contribute only over the
   // vertical range shared by both cells on an edge.
   parallelForOuter(
       "computeBtrBaroclinicPseudoThicknessEdge", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          Real PseudoThicknessSum = 0._Real;
          parallelReduceInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(I4 K, Real &Sum) {
                 Sum += FluxPseudoThickEdge(IEdge, K);
              },
              PseudoThicknessSum);

          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 BaroclinicPseudoThicknessEdge(IEdge) = PseudoThicknessSum;
              });
       });

   for (I4 Subcycle = 0; Subcycle < NBtrSubcycles; ++Subcycle) {

      // Barotropic velocity predictor
      parallelFor(
          "btrVelocityPredictor", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcycleNew(IEdge) = 0._Real;
                return;
             }

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real Mask = EdgeMask(IEdge, KMin);

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
      MeshHalo->exchangeFullArrayHalo(NormalBtrVelSubcycleNew, OnEdge);

      // Barotropic pressure anomaly predictor
      parallelFor(
          "btrPressurePredictor", {Mesh->NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);
                const Real NormalBtrVelEdge =
                    (1._Real - Gamma1) * NormalBtrVelSubcycleCur(JEdge) +
                    Gamma1 * NormalBtrVelSubcycleNew(JEdge);

                const Real DeltaBtrPressure0 =
                    BtrPressAnomalySubcycleCur(Cell0) -
                    BtrPressAnomalyNext(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressAnomalySubcycleCur(Cell1) -
                    BtrPressAnomalyNext(Cell1);
                Real DeltaBtrPressureEdge =
                    0.5_Real * (DeltaBtrPressure0 + DeltaBtrPressure1);
                if (UpwindFluxThickness) {
                   if (NormalBtrVelEdge > 0._Real)
                      DeltaBtrPressureEdge = DeltaBtrPressure0;
                   else if (NormalBtrVelEdge < 0._Real)
                      DeltaBtrPressureEdge = DeltaBtrPressure1;
                   else
                      DeltaBtrPressureEdge =
                          Kokkos::max(DeltaBtrPressure0, DeltaBtrPressure1);
                }

                // This is rho0*g*wctEdge. The pressure-anomaly difference is
                // rho0*g times MPAS-O's deltaSSH.
                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;
                const Real PredictorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;

                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) *
                    PredictorBtrFlux;
             }

             // TODO: Surface freshwater flux (rho0*g*Q) is assumed to be 0 for
             // now.
             BtrPressAnomalySubcycleNew(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });
      MeshHalo->exchangeFullArrayHalo(BtrPressAnomalySubcycleNew, OnCell);

      // Barotropic velocity corrector
      parallelFor(
          "btrVelocityCorrector", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcycleNew(IEdge) = 0._Real;
                return;
             }

             const Real Mask = EdgeMask(IEdge, KMin);
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
      MeshHalo->exchangeFullArrayHalo(NormalBtrVelSubcycleNew, OnEdge);


      // Barotropic pressure anomaly corrector
      parallelFor(
          "btrPressureCorrector", {Mesh->NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);

                const Real BtrPressure0 =
                    (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell0) +
                    Gamma2 * BtrPressAnomalySubcycleNew(Cell0);

                const Real BtrPressure1 =
                    (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell1) +
                    Gamma2 * BtrPressAnomalySubcycleNew(Cell1);

                const Real NormalBtrVelEdge =
                    (1._Real - Gamma3) * NormalBtrVelSubcycleCur(JEdge) +
                    Gamma3 * NormalBtrVelSubcycleNew(JEdge);

                const Real DeltaBtrPressure0 =
                    BtrPressure0 - BtrPressAnomalyNext(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressure1 - BtrPressAnomalyNext(Cell1);
                Real DeltaBtrPressureEdge =
                    0.5_Real * (DeltaBtrPressure0 + DeltaBtrPressure1);
                if (UpwindFluxThickness) {
                   if (NormalBtrVelEdge > 0._Real)
                      DeltaBtrPressureEdge = DeltaBtrPressure0;
                   else if (NormalBtrVelEdge < 0._Real)
                      DeltaBtrPressureEdge = DeltaBtrPressure1;
                   else
                      DeltaBtrPressureEdge =
                          Kokkos::max(DeltaBtrPressure0, DeltaBtrPressure1);
                }

                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;

                const Real CorrectorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;
                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) *
                    CorrectorBtrFlux;
             }

             BtrPressAnomalySubcycleNew(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });
      MeshHalo->exchangeFullArrayHalo(BtrPressAnomalySubcycleNew, OnCell);


      // Accumulate the corrector barotropic flux for the time-averaged
      // transport used by the barotropic-baroclinic velocity correction.
      parallelFor(
          "btrFluxCorrector", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
             if (MaxLayerEdgeTop(IEdge) < MinLayerEdgeBot(IEdge))
                return;

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real BtrPressure0 =
                 (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell0) +
                 Gamma2 * BtrPressAnomalySubcycleNew(Cell0);

             const Real BtrPressure1 =
                 (1._Real - Gamma2) * BtrPressAnomalySubcycleCur(Cell1) +
                 Gamma2 * BtrPressAnomalySubcycleNew(Cell1);

             const Real NormalBtrVelEdge =
                 (1._Real - Gamma3) * NormalBtrVelSubcycleCur(IEdge) +
                 Gamma3 * NormalBtrVelSubcycleNew(IEdge);

             const Real DeltaBtrPressure0 =
                 BtrPressure0 - BtrPressAnomalyNext(Cell0);
             const Real DeltaBtrPressure1 =
                 BtrPressure1 - BtrPressAnomalyNext(Cell1);
             Real DeltaBtrPressureEdge =
                 0.5_Real * (DeltaBtrPressure0 + DeltaBtrPressure1);
             if (UpwindFluxThickness) {
                if (NormalBtrVelEdge > 0._Real)
                   DeltaBtrPressureEdge = DeltaBtrPressure0;
                else if (NormalBtrVelEdge < 0._Real)
                   DeltaBtrPressureEdge = DeltaBtrPressure1;
                else
                   DeltaBtrPressureEdge =
                       Kokkos::max(DeltaBtrPressure0, DeltaBtrPressure1);
             }

             const Real BtrPressureEdge =
                 RhoGravity * BaroclinicPseudoThicknessEdge(IEdge) +
                 DeltaBtrPressureEdge;

             BtrFlux(IEdge) += BtrPressureEdge * NormalBtrVelEdge;
          });

      parallelFor("btrVelocityAccumulate", {Mesh->NEdgesOwned},
                  KOKKOS_LAMBDA(I4 IEdge) {
                     NormalBtrVelNext(IEdge) +=
                         NormalBtrVelSubcycleNew(IEdge);
                  });

      deepCopy(NormalBtrVelSubcycleCur, NormalBtrVelSubcycleNew);
      deepCopy(BtrPressAnomalySubcycleCur, BtrPressAnomalySubcycleNew);
   } // Subcycle

   parallelFor(
       "btrSubcycleAverage", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
          NormalBtrVelNext(IEdge) *= InvBtrVelAvgCount;
          BtrFlux(IEdge) *= InvBtrFluxAvgCount;
       });

   deepCopy(BtrPressAnomalyNext, BtrPressAnomalySubcycleCur);
   MeshHalo->exchangeFullArrayHalo(NormalBtrVelNext, OnEdge);
   MeshHalo->exchangeFullArrayHalo(BtrFlux, OnEdge);
   MeshHalo->exchangeFullArrayHalo(BtrPressAnomalyNext, OnCell);

   Pacer::stop("SE-RK2:stage2BtrPC", 2);
}

} // namespace OMEGA

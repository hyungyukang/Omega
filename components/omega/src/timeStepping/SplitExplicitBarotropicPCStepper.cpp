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
SplitExplicitBarotropicPCStepper::SplitExplicitBarotropicPCStepper(
    const std::string &InName, const AuxiliaryState *InAuxState,
    const HorzMesh *InMesh, Halo *InMeshHalo, const VertCoord *InVCoord)
    : BarotropicStepper(InName, "ForwardBackward", InAuxState, InMesh,
                        InMeshHalo, InVCoord) {

   VelSubcycleCur   = Array1DReal("NormalBarotropicVelocitySubcycleCur" + Name,
                                  Mesh->NEdgesSize);
   VelSubcyclePre   = Array1DReal("NormalBarotropicVelocitySubcyclePre" + Name,
                                  Mesh->NEdgesSize);
   VelSubcycleCor   = Array1DReal("NormalBarotropicVelocitySubcycleCor" + Name,
                                  Mesh->NEdgesSize);
   PressSubcycleCur = Array1DReal("BarotropicPressureAnomalySubcycleCur" + Name,
                                  Mesh->NCellsSize);
   PressSubcyclePre = Array1DReal("BarotropicPressureAnomalySubcyclePre" + Name,
                                  Mesh->NCellsSize);
   PressSubcycleCor = Array1DReal("BarotropicPressureAnomalySubcycleCor" + Name,
                                  Mesh->NCellsSize);
   deepCopy(VelSubcycleCur, 0._Real);
   deepCopy(VelSubcyclePre, 0._Real);
   deepCopy(VelSubcycleCor, 0._Real);
   deepCopy(PressSubcycleCur, 0._Real);
   deepCopy(PressSubcyclePre, 0._Real);
   deepCopy(PressSubcycleCor, 0._Real);
}

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::advance(const BarotropicInput &In,
                                               const BarotropicOutput &Out) {

   if (In.NSubsteps < 1)
      LOG_CRITICAL("Invalid barotropic subcycle count {}", In.NSubsteps);

   Eos *EqState = Eos::getInstance();

   // NSubsteps predictor-corrector steps over the window, which for the
   // split-explicit stepper is twice the baroclinic step
   const I4 NBtrSubcycles        = In.NSubsteps;
   const Real BtrDt              = In.Window / NBtrSubcycles;
   const Real InvBtrVelAvgCount  = 1._Real / (NBtrSubcycles + 1);
   const Real InvBtrFluxAvgCount = 1._Real / NBtrSubcycles;
   // Forward-backward feedback weights for the predictor-corrector subcycle.
   constexpr Real Gamma1     = 0.5333_Real;
   constexpr Real Gamma2     = 0.5333_Real;
   constexpr Real Gamma3     = 1.0_Real;
   constexpr Real RhoGravity = RhoSw * Gravity;

   Pacer::start("SE-RK2:stage2BtrPC", 2);

   // Start of the window, reference pressure anomaly of the flux thickness,
   // and the outputs. The reference is read throughout the subcycle and the
   // pressure output is only written after it, so the two may alias.
   Array1DReal NormalBtrVelCur     = In.VelStart;
   Array1DReal BtrPressAnomalyCur  = In.PressStart;
   Array1DReal BtrPressAnomalyRef  = In.PressRef;
   Array1DReal NormalBtrVelNext    = Out.Vel;
   Array1DReal BtrPressAnomalyNext = Out.Press;

   // Cur: the state at the start of the subcycle
   // Pre: the predictor output
   // Cor: the corrector output
   Array1DReal NormalBtrVelSubcycleCur    = VelSubcycleCur;
   Array1DReal NormalBtrVelSubcyclePre    = VelSubcyclePre;
   Array1DReal NormalBtrVelSubcycleCor    = VelSubcycleCor;
   Array1DReal BtrPressAnomalySubcycleCur = PressSubcycleCur;
   Array1DReal BtrPressAnomalySubcyclePre = PressSubcyclePre;
   Array1DReal BtrPressAnomalySubcycleCor = PressSubcycleCor;

   Array1DReal BtrForcing = In.Forcing;
   Array1DReal BtrFlux    = Out.FluxAvg;

   Array1DReal BaroclinicPseudoThicknessEdge = ColumnFluxThickEdge;

   const bool UpwindFluxThickness = upwindFluxThickness();

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

   // Compute the baroclinic pseudo thickness on edges
   computeColumnFluxThickness();

   //-----------------------------------------------------------------------
   // Communication-avoiding barotropic subcycling:
   //   - Each kernel can fill one halo layer less than its inputs cover,
   //     so the ranges below shrink in turn, and the last of them sets the
   //     floor on HaloWidth. Halo exchange is conducted once at the top
   //     of each subcycle.
   const I4 HaloWidth = static_cast<I4>(Mesh->NCellsHaloH.extent(0));

   // Every kernel downstream loses the layers its stencil reaches past
   // what its inputs cover.
   constexpr I4 VelPredShrink   = 0;
   constexpr I4 PressPredShrink = 1;
   constexpr I4 VelCorrShrink   = 1;
   constexpr I4 PressCorrShrink = 2;

   constexpr I4 MinHaloWidth =
       1 + Kokkos::max(Kokkos::max(VelPredShrink, PressPredShrink),
                       Kokkos::max(VelCorrShrink, PressCorrShrink));

   if (HaloWidth < MinHaloWidth)
      ABORT_ERROR("Split-explicit barotropic subcycling needs a Decomp "
                  "HaloWidth of at least {}, got {}",
                  MinHaloWidth, HaloWidth);

   // A halo deeper than the minimum just leaves every range deeper.
   const I4 HaloAll            = HaloWidth - 1;
   const I4 VelPredEdgeRange   = Mesh->NEdgesHaloH(HaloAll - VelPredShrink);
   const I4 PressPredCellRange = Mesh->NCellsHaloH(HaloAll - PressPredShrink);
   const I4 VelCorrEdgeRange   = Mesh->NEdgesHaloH(HaloAll - VelCorrShrink);
   const I4 PressCorrCellRange = Mesh->NCellsHaloH(HaloAll - PressCorrShrink);

   // Start the barotrolic subcycling
   for (I4 Subcycle = 0; Subcycle < NBtrSubcycles; ++Subcycle) {

      // Restore buffers. This is the only communication inside the subcycle.
      MeshHalo->exchangeFullArrayHalo(NormalBtrVelSubcycleCur, OnEdge);
      MeshHalo->exchangeFullArrayHalo(BtrPressAnomalySubcycleCur, OnCell);

      // Barotropic velocity predictor
      parallelFor(
          "btrVelocityPredictor", {VelPredEdgeRange}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcyclePre(IEdge) = 0._Real;
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

             NormalBtrVelSubcyclePre(IEdge) =
                 Mask * (NormalBtrVelSubcycleCur(IEdge) +
                         BtrDt * (CoriolisTend -
                                  MeanSpecVolEdge * BtrPressAnomalyGrad +
                                  BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly predictor
      parallelFor(
          "btrPressurePredictor", {PressPredCellRange},
          KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);
                const Real NormalBtrVelEdge =
                    btrFeedbackBlend(Gamma1, NormalBtrVelSubcycleCur(JEdge),
                                     NormalBtrVelSubcyclePre(JEdge));

                const Real DeltaBtrPressure0 =
                    BtrPressAnomalySubcycleCur(Cell0) -
                    BtrPressAnomalyRef(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressAnomalySubcycleCur(Cell1) -
                    BtrPressAnomalyRef(Cell1);
                const Real DeltaBtrPressureEdge =
                    btrDeltaPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                         NormalBtrVelEdge, UpwindFluxThickness);

                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;
                const Real PredictorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;

                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * PredictorBtrFlux;
             }

             constexpr Real SurfaceFreshwaterFlux = 0._Real;
             BtrPressAnomalySubcyclePre(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * (BtrFluxDivTend / AreaCell(ICell) +
                          RhoGravity * SurfaceFreshwaterFlux);
          });

      // Barotropic velocity corrector
      parallelFor(
          "btrVelocityCorrector", {VelCorrEdgeRange}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcycleCor(IEdge) = 0._Real;
                return;
             }

             const Real Mask = EdgeMask(IEdge, KMin);
             const I4 Cell0  = CellsOnEdge(IEdge, 0);
             const I4 Cell1  = CellsOnEdge(IEdge, 1);

             Real CoriolisTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
                const I4 JEdge = EdgesOnEdge(IEdge, J);
                CoriolisTend += WeightsOnEdge(IEdge, J) *
                                NormalBtrVelSubcyclePre(JEdge) * FEdge(JEdge);
             }

             const Real BtrPressure0 =
                 btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                                  BtrPressAnomalySubcyclePre(Cell0));

             const Real BtrPressure1 =
                 btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                                  BtrPressAnomalySubcyclePre(Cell1));

             const Real MeanSpecVolEdge =
                 0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));

             const Real BtrPressureGrad =
                 (BtrPressure1 - BtrPressure0) / DcEdge(IEdge);

             NormalBtrVelSubcycleCor(IEdge) =
                 Mask *
                 (NormalBtrVelSubcycleCur(IEdge) +
                  BtrDt * (CoriolisTend - MeanSpecVolEdge * BtrPressureGrad +
                           BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly corrector
      parallelFor(
          "btrPressureCorrector", {PressCorrCellRange},
          KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);

                const Real BtrPressure0 =
                    btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                                     BtrPressAnomalySubcyclePre(Cell0));

                const Real BtrPressure1 =
                    btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                                     BtrPressAnomalySubcyclePre(Cell1));

                const Real NormalBtrVelEdge =
                    btrFeedbackBlend(Gamma3, NormalBtrVelSubcycleCur(JEdge),
                                     NormalBtrVelSubcycleCor(JEdge));

                const Real DeltaBtrPressure0 =
                    BtrPressure0 - BtrPressAnomalyRef(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressure1 - BtrPressAnomalyRef(Cell1);
                const Real DeltaBtrPressureEdge =
                    btrDeltaPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                         NormalBtrVelEdge, UpwindFluxThickness);

                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;

                const Real CorrectorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;
                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * CorrectorBtrFlux;
             }

             BtrPressAnomalySubcycleCor(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });

      // Accumulate the corrector barotropic velocity and flux for the
      // time-averaged transport used by the barotropic-baroclinic velocity
      // correction.
      parallelFor(
          "btrCorrectorAccumulate", {Mesh->NEdgesOwned},
          KOKKOS_LAMBDA(I4 IEdge) {
             NormalBtrVelNext(IEdge) += NormalBtrVelSubcycleCor(IEdge);

             if (MaxLayerEdgeTop(IEdge) < MinLayerEdgeBot(IEdge))
                return;

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real BtrPressure0 =
                 btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                                  BtrPressAnomalySubcycleCor(Cell0));

             const Real BtrPressure1 =
                 btrFeedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                                  BtrPressAnomalySubcycleCor(Cell1));

             const Real NormalBtrVelEdge =
                 btrFeedbackBlend(Gamma3, NormalBtrVelSubcycleCur(IEdge),
                                  NormalBtrVelSubcycleCor(IEdge));

             const Real DeltaBtrPressure0 =
                 BtrPressure0 - BtrPressAnomalyRef(Cell0);
             const Real DeltaBtrPressure1 =
                 BtrPressure1 - BtrPressAnomalyRef(Cell1);
             const Real DeltaBtrPressureEdge =
                 btrDeltaPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                      NormalBtrVelEdge, UpwindFluxThickness);

             const Real BtrPressureEdge =
                 RhoGravity * BaroclinicPseudoThicknessEdge(IEdge) +
                 DeltaBtrPressureEdge;

             BtrFlux(IEdge) += BtrPressureEdge * NormalBtrVelEdge;
          });

      // Rotate the subcycle buffers: the corrector output becomes the state
      // at the start of the next subcycle.
      std::swap(NormalBtrVelSubcycleCur, NormalBtrVelSubcycleCor);
      std::swap(BtrPressAnomalySubcycleCur, BtrPressAnomalySubcycleCor);
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

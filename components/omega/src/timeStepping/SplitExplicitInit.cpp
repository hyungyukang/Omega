//===-- SplitExplicitInit.cpp - split-explicit initialization --*- C++ -*-===//
//
// Utilities for initializing split-explicit time stepping state.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitInit.h"
#include "Config.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"

#include <algorithm>
#include <cmath>

namespace OMEGA {

//------------------------------------------------------------------------------
SplitExplicitConfig
SplitExplicitInit::readConfigOptions(const TimeInterval &TimeStep) {

   SplitExplicitConfig SEConfig;
   SEConfig.BtrTimeStep = TimeStep;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   Error Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   bool IsUnsplit = false;
   std::string TimeStepperStr;
   if (TimeIntConfig.get("TimeStepper", TimeStepperStr).isSuccess()) {
      if (TimeStepperStr.rfind("Unsplit", 0) == 0) {
         IsUnsplit            = true;
         SEConfig.SplitFactor = 0._Real;
      }
   }

   if (!IsUnsplit) {
      std::string BtrTimeStepStr;
      if (TimeIntConfig.get("BtrTimeStep", BtrTimeStepStr).isSuccess()) {
         SEConfig.BtrTimeStep = TimeInterval(BtrTimeStepStr);
      }
   }

   I4 NTimeStepIteration = 1;
   if (TimeIntConfig.get("NTimeStepIteration", NTimeStepIteration)
           .isSuccess()) {
      if (NTimeStepIteration < 1) {
         ABORT_ERROR("NTimeStepIteration must be greater than zero");
      }
      SEConfig.NTimeStepIteration = NTimeStepIteration;
   }

   I4 NBclCoriolisIteration = 2;
   if (TimeIntConfig.get("NBclCoriolisIteration", NBclCoriolisIteration)
           .isSuccess()) {
      if (NBclCoriolisIteration < 1) {
         ABORT_ERROR("NBclCoriolisIteration must be greater than zero");
      }
      SEConfig.NBclCoriolisIteration = NBclCoriolisIteration;
   }
   if (IsUnsplit) {
      SEConfig.NBclCoriolisIteration = 1;
   }

   if (!IsUnsplit) {
      SEConfig.NBtrSubcycles =
          computeSubcycleCount(TimeStep, SEConfig.BtrTimeStep);
   }

   return SEConfig;
}

//------------------------------------------------------------------------------
SplitExplicitBarotropicStepperType
SplitExplicitInit::getBtrTimeStepperFromStr(const std::string &InString) {

   if (InString == "Predictor-Corrector") {
      return SplitExplicitBarotropicStepperType::PredictorCorrector;
   }

   ABORT_ERROR("BtrTimeStepper should be 'Predictor-Corrector' but got {}",
               InString);
   return SplitExplicitBarotropicStepperType::Invalid;
}

//------------------------------------------------------------------------------
I4 SplitExplicitInit::computeSubcycleCount(const TimeInterval &TimeStep,
                                           const TimeInterval &BtrTimeStep) {

   R8 TimeStepSeconds;
   R8 BtrTimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   BtrTimeStep.get(BtrTimeStepSeconds, TimeUnits::Seconds);

   if (BtrTimeStepSeconds <= 0.) {
      ABORT_ERROR("BtrTimeStep must be greater than zero");
   }

   return std::max<I4>(
       1, static_cast<I4>(std::ceil(TimeStepSeconds / BtrTimeStepSeconds)));
}

//------------------------------------------------------------------------------
void SplitExplicitInit::allocateScratch(SplitExplicitScratch &Scratch,
                                        const HorzMesh *Mesh,
                                        const VertCoord *VCoord,
                                        const std::string &Name) {

   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");

   Scratch.NormalBarotropicVelocitySubcycleCur = Array1DReal(
       "NormalBarotropicVelocitySubcycleCur" + Name, Mesh->NEdgesSize);
   Scratch.NormalBarotropicVelocitySubcycleNew = Array1DReal(
       "NormalBarotropicVelocitySubcycleNew" + Name, Mesh->NEdgesSize);
   Scratch.NormalBarotropicVelocityNew =
       Array1DReal("NormalBarotropicVelocityNew" + Name, Mesh->NEdgesSize);
   Scratch.BarotropicPressureAnomalySubcycleCur = Array1DReal(
       "BarotropicPressureAnomalySubcycleCur" + Name, Mesh->NCellsSize);
   Scratch.BarotropicPressureAnomalySubcycleNew = Array1DReal(
       "BarotropicPressureAnomalySubcycleNew" + Name, Mesh->NCellsSize);
   Scratch.BarotropicPressure =
       Array1DReal("BarotropicPressure" + Name, Mesh->NCellsSize);
   Scratch.BarotropicForcing =
       Array1DReal("BarotropicForcing" + Name, Mesh->NEdgesSize);
   Scratch.BarotropicFlux =
       Array1DReal("BarotropicFlux" + Name, Mesh->NEdgesSize);
   Scratch.BaroclinicPseudoThicknessEdge = Array1DReal(
       "BaroclinicPseudoThicknessEdge" + Name, Mesh->NEdgesSize);
   Scratch.BaseVelocityTend = Array2DReal(
       "BaseVelocityTend" + Name, Mesh->NEdgesSize, VCoord->NVertLayers);
   Scratch.NormalTransportVelocity =
       Array2DReal("NormalTransportVelocity" + Name, Mesh->NEdgesSize,
                   VCoord->NVertLayers);

   parallelFor(
       "initializeCell1D", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(I4 ICell) {
          Scratch.BarotropicPressureAnomalySubcycleCur(ICell) = 0._Real;
          Scratch.BarotropicPressureAnomalySubcycleNew(ICell) = 0._Real;
          Scratch.BarotropicPressure(ICell) = 0._Real;
       });

   parallelFor(
       "initializeCell1D", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge) {
          Scratch.NormalBarotropicVelocitySubcycleCur(IEdge) = 0._Real;
          Scratch.NormalBarotropicVelocitySubcycleNew(IEdge) = 0._Real;
          Scratch.NormalBarotropicVelocityNew(IEdge) = 0._Real;
          Scratch.BarotropicForcing(IEdge) = 0._Real;
          Scratch.BarotropicFlux(IEdge) = 0._Real;
          Scratch.BaroclinicPseudoThicknessEdge(IEdge) = 0._Real;
       });

   deepCopy(Scratch.BaseVelocityTend, 0.);
   deepCopy(Scratch.NormalTransportVelocity, 0.);
}

//------------------------------------------------------------------------------
void SplitExplicitInit::computeVelocitySplit(OceanState *State,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocity = State->getNormalVelocity(TimeLevel);
   Array2DReal PseudoThickness = State->getPseudoThickness(TimeLevel);
   Array2DReal NormalBaroclinicVelocity =
       State->getNormalBaroclinicVelocity(TimeLevel);
   Array1DReal NormalBarotropicVelocity =
       State->getNormalBarotropicVelocity(TimeLevel);

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);

   deepCopy(NormalBaroclinicVelocity, 0.);
   deepCopy(NormalBarotropicVelocity, 0.);

   parallelForOuter(
       "SplitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {

          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          Real BarotropicVelocity = 0._Real;

          if ( KMax >= KMin ) {
             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             Real ThicknessSum = 0._Real;
             Real FluxSum      = 0._Real;

             parallelReduceInner(
                Team, Range{KMin, KMax},
                INNER_LAMBDA(const int K, Real &ThickAccum, Real &FluxAccum) {
                   const Real ThickEdge = 0.5_Real * (PseudoThickness(Cell0, K) +
                                                      PseudoThickness(Cell1, K));

                   ThickAccum += ThickEdge;
                   FluxAccum += ThickEdge * NormalVelocity(IEdge, K);
                },
                ThicknessSum, FluxSum);

             BarotropicVelocity = FluxSum / ThicknessSum;

             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    NormalBarotropicVelocity(IEdge) =
                        BarotropicVelocity * EdgeMask(IEdge, KMin);
                 });

          } else {

             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    NormalBarotropicVelocity(IEdge) = 0._Real;
                 });

          }

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(I4 K) {
                  NormalBaroclinicVelocity(IEdge, K) =
                      NormalVelocity(IEdge, K) - BarotropicVelocity;
              });
       });

}

//------------------------------------------------------------------------------
void SplitExplicitInit::computeUnsplitVelocitySplit(OceanState *State,
                                                    const HorzMesh *Mesh,
                                                    const VertCoord *VCoord,
                                                    I4 CurLevel, I4 NextLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocityCur = State->getNormalVelocity(CurLevel);
   Array2DReal NormalVelocityNew = State->getNormalVelocity(NextLevel);
   Array2DReal NormalBaroclinicVelocityCur =
       State->getNormalBaroclinicVelocity(CurLevel);
   Array2DReal NormalBaroclinicVelocityNew =
       State->getNormalBaroclinicVelocity(NextLevel);
   Array1DReal NormalBarotropicVelocityCur =
       State->getNormalBarotropicVelocity(CurLevel);
   Array1DReal NormalBarotropicVelocityNew =
       State->getNormalBarotropicVelocity(NextLevel);

   deepCopy(NormalBaroclinicVelocityCur, 0.);
   deepCopy(NormalBaroclinicVelocityNew, 0.);

   deepCopy(NormalBarotropicVelocityCur, 0.);
   deepCopy(NormalBarotropicVelocityNew, 0.);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "UnsplitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {

          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(I4 K) {
                  NormalBaroclinicVelocityCur(IEdge, K) =
                  NormalVelocityCur(IEdge, K);
                  NormalBaroclinicVelocityNew(IEdge, K) =
                  NormalVelocityNew(IEdge, K);
              });
       });

}

//------------------------------------------------------------------------------
void SplitExplicitInit::initializeBarotropicPressure(
    SplitExplicitScratch &Scratch, OceanState *State, const HorzMesh *Mesh,
    const VertCoord *VCoord, I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");

   Array1DReal BtrPressure        = Scratch.BarotropicPressure;
   Array1DReal BtrPressAnomaly    = State->getBarotropicPressureAnomaly(TimeLevel);
   Array1DReal SurfacePressure    = VCoord->SurfacePressure;
   Array2DReal PressureInterface  = VCoord->PressureInterface;
   Array1DReal BottomGeomDepth    = VCoord->BottomGeomDepth;
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelFor(
       "initializeBarotropicPressure", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell) {
          const I4 KMax = MaxLayerCell(ICell);

          const Real Pressure =
              PressureInterface(ICell, KMax + 1) - SurfacePressure(ICell);
          BtrPressure(ICell)     = Pressure;
          BtrPressAnomaly(ICell) = Pressure - RhoSw * Gravity * BottomGeomDepth(ICell);
          Scratch.BarotropicPressureAnomalySubcycleCur(ICell) = BtrPressAnomaly(ICell);
          Scratch.BarotropicPressureAnomalySubcycleNew(ICell) = BtrPressAnomaly(ICell);
       });
}

//------------------------------------------------------------------------------
void SplitExplicitInit::combineVelocitySplit(OceanState *State,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocity = State->getNormalVelocity(TimeLevel);
   Array2DReal NormalBaroclinicVelocity =
       State->getNormalBaroclinicVelocity(TimeLevel);
   Array1DReal NormalBarotropicVelocity =
       State->getNormalBarotropicVelocity(TimeLevel);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelFor(
       "combineSplitExplicitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          const Real BarotropicVelocity = NormalBarotropicVelocity(IEdge);
          for (I4 K = KMin; K <= KMax; ++K) {
             NormalVelocity(IEdge, K) =
                 NormalBaroclinicVelocity(IEdge, K) + BarotropicVelocity;
          }
       });
}

} // namespace OMEGA

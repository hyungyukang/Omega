//===-- SplitExplicitInit.cpp - split-explicit initialization --*- C++ -*-===//
//
// Utilities for initializing split-explicit time stepping state.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitInit.h"
#include "Config.h"
#include "Error.h"
#include "Logging.h"
#include "OmegaKokkos.h"

#include <algorithm>
#include <cmath>

namespace OMEGA {

//------------------------------------------------------------------------------
SplitExplicitConfig
SplitExplicitInit::readConfigOptions(const TimeInterval &TimeStep) {

   SplitExplicitConfig Options;
   Options.BtrTimeStep = TimeStep;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   Error Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   std::string TimeStepperStr;
   if (TimeIntConfig.get("TimeStepper", TimeStepperStr).isSuccess()) {
      if (TimeStepperStr == "Unsplit-RK2" || TimeStepperStr == "UnsplitRK2") {
         Options.SplitFactor = 0._Real;
      }
   }

   std::string BtrTimeStepStr;
   if (TimeIntConfig.get("BtrTimeStep", BtrTimeStepStr).isSuccess()) {
      Options.BtrTimeStep = TimeInterval(BtrTimeStepStr);
   }

   std::string BtrTimeStepperStr = "Predictor-Corrector";
   if (TimeIntConfig.get("BtrTimeStepper", BtrTimeStepperStr).isSuccess()) {
      Options.BtrTimeStepper =
          getBtrTimeStepperFromStr(BtrTimeStepperStr);
   }

   I4 NTimeStepIteration = 1;
   if (TimeIntConfig.get("NTimeStepIteration", NTimeStepIteration).isSuccess()) {
      if (NTimeStepIteration < 1) {
         ABORT_ERROR("NTimeStepIteration must be greater than zero");
      }
      Options.NTimeStepIteration = NTimeStepIteration;
   }

   I4 NBclCoriolisIteration = 2;
   if (TimeIntConfig.get("NBclCoriolisIteration",
                         NBclCoriolisIteration).isSuccess()) {
      if (NBclCoriolisIteration < 1) {
         ABORT_ERROR("NBclCoriolisIteration must be greater than zero");
      }
      Options.NBclCoriolisIteration = NBclCoriolisIteration;
   }

   Options.NBtrSubcycles =
       computeSubcycleCount(TimeStep, Options.BtrTimeStep);

   return Options;
}

//------------------------------------------------------------------------------
SplitExplicitBarotropicStepperType
SplitExplicitInit::getBtrTimeStepperFromStr(const std::string &InString) {

   if (InString == "Predictor-Corrector" ||
       InString == "PredictorCorrector" || InString == "PC" ||
       InString == "FBPC" || InString == "Forward-Backward-Predictor-Corrector") {
      return SplitExplicitBarotropicStepperType::PredictorCorrector;
   }

   ABORT_ERROR("BtrTimeStepper should be one of 'Predictor-Corrector' or "
               "'PC' but got {}",
               InString);
   return SplitExplicitBarotropicStepperType::Invalid;
}

//------------------------------------------------------------------------------
I4 SplitExplicitInit::computeSubcycleCount(
    const TimeInterval &TimeStep, const TimeInterval &BtrTimeStep) {

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
   Scratch.BaseVelocityTend =
       Array2DReal("BaseVelocityTend" + Name, Mesh->NEdgesSize,
                   VCoord->NVertLayers);

   deepCopy(Scratch.NormalBarotropicVelocitySubcycleCur, 0.);
   deepCopy(Scratch.NormalBarotropicVelocitySubcycleNew, 0.);
   deepCopy(Scratch.BarotropicPressureAnomalySubcycleCur, 0.);
   deepCopy(Scratch.BarotropicPressureAnomalySubcycleNew, 0.);
   deepCopy(Scratch.BarotropicPressure, 0.);
   deepCopy(Scratch.BarotropicForcing, 0.);
   deepCopy(Scratch.BarotropicFlux, 0.);
   deepCopy(Scratch.BaseVelocityTend, 0.);
}

//------------------------------------------------------------------------------
void SplitExplicitInit::computeVelocitySplit(OceanState *State,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocity = State->getNormalVelocity(TimeLevel);
   Array2DReal LayerThickness = State->getLayerThickness(TimeLevel);
   Array2DReal NormalBaroclinicVelocity =
       State->getNormalBaroclinicVelocity(TimeLevel);
   Array1DReal NormalBarotropicVelocity =
       State->getNormalBarotropicVelocity(TimeLevel);

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   deepCopy(NormalBaroclinicVelocity, 0.);
   deepCopy(NormalBarotropicVelocity, 0.);

   parallelFor(
       "computeSplitExplicitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          if (KMax < KMin) {
             NormalBarotropicVelocity(IEdge) = 0._Real;
             return;
          }

          const I4 Cell1 = CellsOnEdge(IEdge, 0);
          const I4 Cell2 = CellsOnEdge(IEdge, 1);

          Real ThicknessSum = 0._Real;
          Real FluxSum      = 0._Real;
          for (I4 K = KMin; K <= KMax; ++K) {
             const Real EdgeThickness =
                 0.5_Real *
                 (LayerThickness(Cell1, K) + LayerThickness(Cell2, K));
             ThicknessSum += EdgeThickness;
             FluxSum += EdgeThickness * NormalVelocity(IEdge, K);
          }

          const Real BarotropicVelocity =
              ThicknessSum > 0._Real ? FluxSum / ThicknessSum : 0._Real;
          NormalBarotropicVelocity(IEdge) = BarotropicVelocity;

          for (I4 K = KMin; K <= KMax; ++K) {
             NormalBaroclinicVelocity(IEdge, K) =
                 NormalVelocity(IEdge, K) - BarotropicVelocity;
          }
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

          if (KMax < KMin) {
             return;
          }

          const Real BarotropicVelocity = NormalBarotropicVelocity(IEdge);
          for (I4 K = KMin; K <= KMax; ++K) {
             NormalVelocity(IEdge, K) =
                 NormalBaroclinicVelocity(IEdge, K) + BarotropicVelocity;
          }
       });
}

} // namespace OMEGA

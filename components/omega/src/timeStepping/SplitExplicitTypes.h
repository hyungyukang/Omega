#ifndef OMEGA_SPLIT_EXPLICIT_TYPES_H
#define OMEGA_SPLIT_EXPLICIT_TYPES_H
//===-- SplitExplicitTypes.h - split-explicit shared types -----*- C++ -*-===//
//
/// \file
/// \brief Contains shared configuration and scratch arrays for SE time
/// stepping.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "TimeMgr.h"

namespace OMEGA {

enum class SplitExplicitBarotropicStepperType { PredictorCorrector, Invalid };

struct SplitExplicitConfig {
   TimeInterval BtrTimeStep;
   SplitExplicitBarotropicStepperType BtrTimeStepper =
       SplitExplicitBarotropicStepperType::PredictorCorrector;
   I4 NBtrSubcycles         = 1;
   I4 NTimeStepIteration    = 1;
   I4 NBclCoriolisIteration = 2;
   bool ReinitializeVelocitySplitEachTimeStep = false;
   Real SplitFactor         = 1._Real;
};

struct SplitExplicitScratch {
   Array1DReal NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBarotropicVelocitySubcycleNew;
   Array1DReal BarotropicPressureAnomalySubcycleCur;
   Array1DReal BarotropicPressureAnomalySubcycleNew;
   Array1DReal BarotropicPressure;
   Array1DReal BarotropicForcing;
   Array1DReal BarotropicFlux;
   Array1DReal BaroclinicPseudoThicknessEdge;
   Array2DReal BaseVelocityTend;
   Array2DReal NormalTransportVelocity;
};

} // namespace OMEGA

#endif

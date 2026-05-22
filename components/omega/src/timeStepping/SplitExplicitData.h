#ifndef OMEGA_SPLIT_EXPLICIT_DATA_H
#define OMEGA_SPLIT_EXPLICIT_DATA_H
//===-- SplitExplicitData.h - split-explicit shared data -------*- C++ -*-===//
//
/// \file
/// \brief Contains shared configuration and scratch arrays for SE time stepping.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "TimeMgr.h"

namespace OMEGA {

enum class SplitExplicitBarotropicStepperType {
   PredictorCorrector,
   Invalid
};

struct SplitExplicitConfig {
   TimeInterval BtrTimeStep;
   SplitExplicitBarotropicStepperType BtrTimeStepper =
       SplitExplicitBarotropicStepperType::PredictorCorrector;
   I4 NBtrSubcycles = 1;
   I4 NSuperCycle   = 1;
};

struct SplitExplicitBuffers {
   Array1DReal NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBarotropicVelocitySubcycleNew;
   Array1DReal BarotropicPressure;
   Array1DReal BarotropicForcing;
   Array1DReal BarotropicFlux;
};

} // namespace OMEGA

#endif

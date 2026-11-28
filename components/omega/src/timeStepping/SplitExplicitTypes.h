#ifndef OMEGA_SPLIT_EXPLICIT_TYPES_H
#define OMEGA_SPLIT_EXPLICIT_TYPES_H
//===-- SplitExplicitTypes.h - split-explicit shared types -----*- C++ -*-===//
//
/// \file
/// \brief Contains shared configuration and scratch arrays for SE time
/// stepping.
//
//===----------------------------------------------------------------------===//

#include "BarotropicStepper.h"
#include "DataTypes.h"
#include "TimeMgr.h"

namespace OMEGA {

struct SplitExplicitConfig {
   TimeInterval BtrTimeStep;
   // Barotropic stepper options (BtrTimeStepper and the Barotropic subgroup)
   BarotropicStepperConfig BtrConfig;
   I4 NBtrSubcycles = 1;
   // Two iterations are required for the RK2 predictor-corrector structure;
   // a single iteration leaves only the predictor and degrades the scheme to
   // forward Euler.
   I4 NTimeStepIteration    = 2;
   I4 NBclCoriolisIteration = 2;
   bool ReinitSplitVelocity = false;
   Real SplitFactor         = 1._Real;
};

struct SplitExplicitScratch {
   // The barotropic subcycle buffers are owned by the barotropic stepper.
   Array1DReal BarotropicForcing;
   Array1DReal BarotropicFlux;
   Array2DReal IterVelocityTend;
   Array2DReal TransportVelocityAdd;
};

} // namespace OMEGA

#endif

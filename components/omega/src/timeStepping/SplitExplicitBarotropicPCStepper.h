#ifndef OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
#define OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
//===-- SplitExplicitBarotropicPCStepper.h - SE stage 2 ------*- C++ -*-===//
//
/// \file
/// \brief Forward-backward predictor-corrector barotropic subcycle framework.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "OceanState.h"
#include "SplitExplicitData.h"
#include "TimeMgr.h"
#include "VertCoord.h"

namespace OMEGA {

class SplitExplicitBarotropicPCStepper {
 public:
   void doSplitStage2(
       OceanState *State,             ///< [inout] model state
       SplitExplicitScratch &Scratch, ///< [inout] split-explicit scratch data
       const SplitExplicitConfig &Config, ///< [in] split-explicit options
       const HorzMesh *Mesh,              ///< [in] horizontal mesh
       const VertCoord *VCoord,           ///< [in] vertical coordinate
       I4 TimeLevel,                      ///< [in] state time level to update
       const TimeInstant &StageTime,      ///< [in] current stage time
       const TimeInterval &StageTimeStep  ///< [in] current stage time step
   ) const;
};

} // namespace OMEGA

#endif

#ifndef OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
#define OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
//===-- SplitExplicitBarotropicPCStepper.h - SE stage 2 ------*- C++ -*-===//
//
/// \file
/// \brief Forward-backward predictor-corrector barotropic subcycle framework.
//
//===----------------------------------------------------------------------===//

#include "Halo.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "SplitExplicitTypes.h"
#include "TimeMgr.h"
#include "VertCoord.h"

namespace OMEGA {

class AuxiliaryState;

class SplitExplicitBarotropicPCStepper {
 public:
   void doBarotropicVelocityUpdate(
       OceanState *State,              ///< [inout] model state
       const AuxiliaryState *AuxState, ///< [in] provisional auxiliary state
       SplitExplicitScratch &Scratch,  ///< [inout] split-explicit scratch data
       const SplitExplicitConfig &Config, ///< [in] split-explicit options
       const HorzMesh *Mesh,              ///< [in] horizontal mesh
       Halo *MeshHalo,                    ///< [in] mesh halo exchange
       const VertCoord *VCoord,           ///< [in] vertical coordinate
       I4 CurLevel,                       ///< [in] state time level to update
       I4 NextLevel,                      ///< [in] state time level to update
       const TimeInstant &StageTime,      ///< [in] current stage time
       const TimeInterval &StageTimeStep  ///< [in] current stage time step
   ) const;
};

} // namespace OMEGA

#endif

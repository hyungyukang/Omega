#ifndef OMEGA_SPLIT_EXPLICIT_INIT_H
#define OMEGA_SPLIT_EXPLICIT_INIT_H
//===-- SplitExplicitInit.h - split-explicit initialization ----*- C++ -*-===//
//
/// \file
/// \brief Initialization helpers for split-explicit time stepping.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "OceanState.h"
#include "SplitExplicitData.h"
#include "VertCoord.h"

namespace OMEGA {

class SplitExplicitInit {
 public:
   static SplitExplicitConfig readConfigOptions(
       const TimeInterval &TimeStep ///< [in] full baroclinic time step
   );

   static void allocateBuffers(
       SplitExplicitBuffers &Buffers, ///< [inout] split-explicit scratch data
       const HorzMesh *Mesh,          ///< [in] horizontal mesh
       const std::string &Name        ///< [in] owning time stepper name
   );

   static void computeVelocitySplit(
       OceanState *State,       ///< [inout] ocean state
       const HorzMesh *Mesh,    ///< [in] horizontal mesh
       const VertCoord *VCoord, ///< [in] vertical coordinate
       I4 TimeLevel             ///< [in] state time level to split
   );

   static void combineVelocitySplit(
       OceanState *State,       ///< [inout] ocean state
       const HorzMesh *Mesh,    ///< [in] horizontal mesh
       const VertCoord *VCoord, ///< [in] vertical coordinate
       I4 TimeLevel             ///< [in] state time level to update
   );

 private:
   static SplitExplicitBarotropicStepperType getBtrTimeStepperFromStr(
       const std::string &InString ///< [in] barotropic stepping method
   );

   static I4 computeSubcycleCount(
       const TimeInterval &TimeStep,   ///< [in] full baroclinic time step
       const TimeInterval &BtrTimeStep ///< [in] requested barotropic time step
   );
};

} // namespace OMEGA

#endif

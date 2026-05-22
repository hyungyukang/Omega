//===-- SplitExplicitBarotropicPCStepper.cpp - SE stage 2 -----*- C++ -*-===//
//
// Framework for the explicitly subcycled forward-backward
// predictor-corrector barotropic velocity update.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "SplitExplicitInit.h"

namespace OMEGA {

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::doSplitStage2(
    OceanState *State, SplitExplicitBuffers &Buffers,
    const SplitExplicitConfig &Config, const HorzMesh *Mesh,
    const VertCoord *VCoord, I4 TimeLevel, const TimeInstant &StageTime,
    const TimeInterval &StageTimeStep) const {

   if (!State)
      LOG_CRITICAL("Invalid State");
   if (Config.NBtrSubcycles < 1)
      LOG_CRITICAL("Invalid split-explicit barotropic subcycle count");

   // Keep the arguments live in this framework implementation. The barotropic
   // tendency terms will use them once the full split-explicit equations land.
   (void)StageTime;
   (void)StageTimeStep;

   Pacer::start("SE-RK2:stage2BtrPC", 2);

   Array1DReal NormalBarotropicVelocity =
       State->getNormalBarotropicVelocity(TimeLevel);

   deepCopy(Buffers.NormalBarotropicVelocitySubcycleCur,
            NormalBarotropicVelocity);

   for (I4 Subcycle = 0; Subcycle < Config.NBtrSubcycles; ++Subcycle) {
      // Placeholder for the barotropic velocity/SSH predictor-corrector. For
      // now the framework preserves the initialized barotropic mode exactly.
      deepCopy(Buffers.NormalBarotropicVelocitySubcycleNew,
               Buffers.NormalBarotropicVelocitySubcycleCur);
      deepCopy(Buffers.NormalBarotropicVelocitySubcycleCur,
               Buffers.NormalBarotropicVelocitySubcycleNew);
   }

   deepCopy(NormalBarotropicVelocity,
            Buffers.NormalBarotropicVelocitySubcycleCur);
   SplitExplicitInit::combineVelocitySplit(State, Mesh, VCoord, TimeLevel);

   Pacer::stop("SE-RK2:stage2BtrPC", 2);
}

} // namespace OMEGA

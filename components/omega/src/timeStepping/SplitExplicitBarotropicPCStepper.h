#ifndef OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
#define OMEGA_SPLIT_EXPLICIT_BAROTROPIC_PC_STEPPER_H
//===-- SplitExplicitBarotropicPCStepper.h - SE stage 2 ------*- C++ -*-===//
//
/// \file
/// \brief Forward-backward predictor-corrector barotropic subcycle framework.
///
/// This is the Omega forward-backward predictor-corrector barotropic stepper
/// (BtrTimeStepper Predictor-Corrector, or Explicit with ExplicitMethod
/// ForwardBackward). It subcycles the barotropic system with NSubsteps
/// predictor-corrector steps over the window and returns the average of the
/// corrected barotropic velocity and pressure flux over the subcycles.
//
//===----------------------------------------------------------------------===//

#include "BarotropicStepper.h"

namespace OMEGA {

class SplitExplicitBarotropicPCStepper : public BarotropicStepper {
 public:
   SplitExplicitBarotropicPCStepper(
       const std::string &InName,        ///< [in] name of the owner
       const AuxiliaryState *InAuxState, ///< [in] provisional auxiliary state
       const HorzMesh *InMesh,           ///< [in] horizontal mesh
       Halo *InMeshHalo,                 ///< [in] mesh halo exchange
       const VertCoord *InVCoord         ///< [in] vertical coordinate
   );

   /// Subcycles the barotropic system over the window of In
   void advance(const BarotropicInput &In,  ///< [in] inputs
                const BarotropicOutput &Out ///< [out] outputs
                ) override;

 private:
   // The subcycle keeps three buffers per field: the state at the start of
   // the subcycle (Cur), the predictor output (Pre) and the corrector output
   // (Cor).
   Array1DReal VelSubcycleCur;
   Array1DReal VelSubcyclePre;
   Array1DReal VelSubcycleCor;
   Array1DReal PressSubcycleCur;
   Array1DReal PressSubcyclePre;
   Array1DReal PressSubcycleCor;
};

} // namespace OMEGA
#endif

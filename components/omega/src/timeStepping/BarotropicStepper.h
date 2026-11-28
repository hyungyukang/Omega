#ifndef OMEGA_BAROTROPIC_STEPPER_H
#define OMEGA_BAROTROPIC_STEPPER_H
//===-- BarotropicStepper.h - barotropic split steppers -------*- C++ -*-===//
//
/// \file
/// \brief Interface of the barotropic (fast mode) steppers of split time
/// stepping
///
/// A barotropic stepper advances the barotropic normal velocity ubar (edges)
/// and the barotropic pressure anomaly B' (cells) over a window [0, Window]
/// with the barotropic forcing Gbar held fixed,
///    d ubar / dt = -f ubar_perp - alphabar grad B' + Gbar
///    d B' / dt   = -div( P(B') ubar )
/// where P(B') is the column pseudo-thickness on edges times rho0 g, corrected
/// by the departure of B' from a reference B'*. It returns the barotropic
/// velocity, the pressure anomaly at the end of the window and the
/// time-averaged barotropic flux used for the transport-velocity correction.
///
/// Every input and output is passed as an argument, so the same stepper can
/// be called from any baroclinic stepper: the Omega SplitExplicitRK2 stepper
/// calls it once per outer iteration, and a SUNDIALS baroclinic stepper can
/// call it from its stage right-hand side with a stage-dependent window.
///
/// The steppers are selected in the ModeSplitShare subgroup of the
/// TimeIntegration configuration:
/// \ConfigInput
/// TimeIntegration:
///   ModeSplitShare:
///     # Predictor-Corrector (the Omega forward-backward predictor-corrector,
///     # kept for compatibility), Explicit or Implicit
///     BtrTimeStepper: Explicit
///     # Requested barotropic time step of the explicit steppers. Ignored by
///     # the implicit stepper, which takes one step per baroclinic step.
///     BtrTimeStep: 0000_00:00:20
///     Barotropic:
///       # Explicit steppers: ForwardBackward (the Omega predictor-corrector),
///       # RK2 (explicit midpoint), RK3 (SSP RK3), RK4, or the name of any
///       # SUNDIALS explicit Butcher table, integrated with ARKODE ERKStep
///       ExplicitMethod: ForwardBackward
///       # Implicit stepper: weight of the new time level (0.5 is
///       # Crank-Nicolson), Krylov solver (PCG or BiCGStab) for the
///       # pseudo-thickness (B') equation, its preconditioner (None or
///       # Jacobi), relative and absolute tolerances on the residual 2-norm
///       # and maximum number of iterations
///       Theta: 0.5
///       LinearSolver: PCG
///       Preconditioner: Jacobi
///       RelTolerance: 1.0e-10
///       AbsTolerance: 0.0
///       MaxIterations: 1000
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <memory>
#include <string>

namespace OMEGA {

class AuxiliaryState;

/// Barotropic stepping algorithms
enum class BarotropicMethod {
   ForwardBackward,  ///< Omega forward-backward predictor-corrector
   SundialsExplicit, ///< explicit Runge-Kutta with SUNDIALS ERKStep
   SundialsImplicit, ///< theta method with a SUNDIALS Krylov solver
   Invalid
};

/// Krylov solvers of the implicit barotropic stepper
enum class BarotropicLinearSolver { PCG, BiCGStab };

/// Preconditioners of the implicit barotropic stepper
enum class BarotropicPreconditioner { None, Jacobi };

/// Options of the barotropic steppers
struct BarotropicStepperConfig {
   BarotropicMethod Method = BarotropicMethod::ForwardBackward;
   /// Name of the explicit method (ForwardBackward or an ERK table)
   std::string ExplicitMethod = "ForwardBackward";
   /// Weight of the new time level in the implicit stepper
   Real Theta                              = 0.5_Real;
   BarotropicLinearSolver LinearSolver     = BarotropicLinearSolver::PCG;
   BarotropicPreconditioner Preconditioner = BarotropicPreconditioner::Jacobi;
   Real RelTolerance                       = 1.e-10_Real;
   Real AbsTolerance                       = 0._Real;
   I4 MaxIterations                        = 1000;
};

/// Inputs of one barotropic advance. All arrays must have valid halos.
struct BarotropicInput {
   Array1DReal VelStart;   ///< ubar at the start of the window (edges)
   Array1DReal PressStart; ///< B' at the start of the window (cells)
   Array1DReal PressRef;   ///< reference B'* of the flux thickness (cells)
   Array1DReal VelRef;     ///< provisional ubar at the window end (edges),
                           ///< used to center the implicit Coriolis term
   Array1DReal Forcing;    ///< barotropic forcing Gbar (edges)
   Real Window  = 0;       ///< window length [s]
   I4 NSubsteps = 1;       ///< explicit substeps over the window
};

/// Outputs of one barotropic advance. Vel may alias VelRef and Press may
/// alias PressRef, but no output may alias VelStart or PressStart. Halos are
/// valid on return.
struct BarotropicOutput {
   Array1DReal Vel;     ///< barotropic velocity: the window average for the
                        ///< explicit steppers, the end value for the
                        ///< implicit stepper (edges)
   Array1DReal Press;   ///< B' at the end of the window (cells)
   Array1DReal FluxAvg; ///< time-averaged barotropic flux (edges)
};

/// Base class of the barotropic steppers
class BarotropicStepper {
 public:
   virtual ~BarotropicStepper() = default;

   /// Reads the barotropic stepper options from the ModeSplitShare subgroup
   /// of the TimeIntegration configuration for the given BtrTimeStepper
   static BarotropicStepperConfig
   readConfig(const std::string &BtrTimeStepperStr ///< [in] BtrTimeStepper
   );

   /// Creates the barotropic stepper selected by Cfg
   static std::unique_ptr<BarotropicStepper>
   create(const BarotropicStepperConfig &Cfg, ///< [in] options
          const std::string &Name,            ///< [in] name of the owner
          const AuxiliaryState *AuxState,     ///< [in] provisional aux state
          const HorzMesh *Mesh,               ///< [in] horizontal mesh
          Halo *MeshHalo,                     ///< [in] halo exchange
          const VertCoord *VCoord             ///< [in] vertical coordinate
   );

   /// Advances the barotropic system over the window of In
   virtual void advance(const BarotropicInput &In,  ///< [in] inputs
                        const BarotropicOutput &Out ///< [out] outputs
                        ) = 0;

   /// True for steppers that take one implicit step over the window
   virtual bool isImplicit() const { return false; }

   /// Window and number of substeps for one baroclinic step of length
   /// TimeStepSeconds: [0, 2 dt] in 2 NBtrSubcycles substeps for the explicit
   /// steppers (the extended averaging window of the split-explicit scheme)
   /// and [0, dt] in one step for the implicit stepper
   void getDefaultWindow(Real TimeStepSeconds, ///< [in] baroclinic step [s]
                         I4 NBtrSubcycles,     ///< [in] subcycles per step
                         Real &Window,         ///< [out] window length [s]
                         I4 &NSubsteps         ///< [out] substeps
   ) const;

   /// Name of the method
   const std::string &getMethodName() const { return MethodName; }

 protected:
   BarotropicStepper(const std::string &InName, const std::string &InMethod,
                     const AuxiliaryState *InAuxState, const HorzMesh *InMesh,
                     Halo *InMeshHalo, const VertCoord *InVCoord);

   /// Sums the flux pseudo thickness of the provisional auxiliary state over
   /// the active layers of every edge into ColumnFluxThickEdge
   void computeColumnFluxThickness() const;

   /// True if the flux thickness is upwinded rather than centered
   bool upwindFluxThickness() const;

   std::string Name;
   std::string MethodName;
   const AuxiliaryState *AuxState = nullptr;
   const HorzMesh *Mesh           = nullptr;
   Halo *MeshHalo                 = nullptr;
   const VertCoord *VCoord        = nullptr;

   /// Column sum of the flux pseudo thickness on edges
   Array1DReal ColumnFluxThickEdge;
};

/// Blends an old and a new value with the forward-backward feedback weight
/// Gamma, i.e. (1-Gamma)*Old + Gamma*New.
KOKKOS_INLINE_FUNCTION Real btrFeedbackBlend(Real Gamma, Real Old, Real New) {
   return (1._Real - Gamma) * Old + Gamma * New;
}

/// Edge value of the barotropic pressure deviation from the reference
/// column, centered or upwinded to match the flux-thickness choice. Where the
/// velocity is exactly zero the upwind direction is undefined and the larger
/// value is taken.
KOKKOS_INLINE_FUNCTION Real btrDeltaPressureEdge(Real Delta0, Real Delta1,
                                                 Real NormalBtrVelEdge,
                                                 bool Upwind) {
   if (!Upwind)
      return 0.5_Real * (Delta0 + Delta1);
   if (NormalBtrVelEdge > 0._Real)
      return Delta0;
   if (NormalBtrVelEdge < 0._Real)
      return Delta1;
   return Kokkos::max(Delta0, Delta1);
}

/// Barotropic flux P(B') ubar on an edge, with P the column flux thickness
/// times rho0 g plus the edge value of the deviation of B' from its reference
KOKKOS_INLINE_FUNCTION Real btrEdgeFlux(Real RhoGravity, Real ColumnThickEdge,
                                        Real Delta0, Real Delta1,
                                        Real NormalBtrVelEdge, bool Upwind) {
   const Real PressEdge =
       RhoGravity * ColumnThickEdge +
       btrDeltaPressureEdge(Delta0, Delta1, NormalBtrVelEdge, Upwind);
   return PressEdge * NormalBtrVelEdge;
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif

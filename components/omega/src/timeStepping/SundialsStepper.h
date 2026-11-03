#ifndef OMEGA_SUNDIALSSTEPPER_H
#define OMEGA_SUNDIALSSTEPPER_H
//===-- SundialsStepper.h - explicit Runge-Kutta from SUNDIALS --*- C++ -*-===//
//
/// \file
/// \brief Contains the class for explicit Runge-Kutta methods from SUNDIALS
///
/// The SundialsStepper advances the state with the explicit Runge-Kutta
/// stepper (ERKStep) of the SUNDIALS ARKODE package, so the method is chosen
/// at run time from the configuration: any explicit Butcher table SUNDIALS
/// provides, or one given in the configuration. The method is set in the
/// Sundials subgroup of the TimeIntegration configuration:
/// \ConfigInput
/// TimeIntegration:
///   TimeStepper: Sundials
///   Sundials:
///     # Explicit Runge-Kutta method. One of the shorthands ForwardEuler,
///     # RK2 (explicit midpoint), RK3 (SSP RK3 of Shu and Osher) or RK4
///     # (classic fourth order), the name of any SUNDIALS explicit Butcher
///     # table (eg ARKODE_HEUN_EULER_2_1_2), or Custom for the table given
///     # in CustomTable
///     Method: RK4
///     # Only read when Method is Custom. A is the Butcher matrix by rows
///     # and must be strictly lower triangular. The coefficients are checked
///     # against the order conditions for the given order.
///     CustomTable:
///       Stages: 3
///       Order: 3
///       A: [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, -1.0, 2.0, 0.0]
///       B: [0.16666666666666667, 0.66666666666666667, 0.16666666666666667]
///       C: [0.0, 0.5, 1.0]
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "SundialsVector.h"
#include "TimeStepper.h"

#include <arkode/arkode_butcher.h>
#include <array>
#include <map>
#include <memory>
#include <sundials/sundials_context.h>
#include <utility>
#include <vector>

namespace OMEGA {

class SundialsStepper : public TimeStepper {
 public:
   /// Constructor creates an instance of the SUNDIALS stepper, fills in some
   /// time information and loads the Butcher table chosen in the
   /// configuration. Data pointers are added later.
   SundialsStepper(
       const std::string &InName,              ///< [in] name of time stepper
       const TimeInterval &InTimeStep,         ///< [in] time step
       const TimeStepperStartType InStartType, ///< [in] option for starting
       const TimeInstant &InStartTime,         ///< [in] start time for full sim
       const TimeStepperStopType InStopType,   ///< [in] option for stopping
       std::optional<TimeInstant> InStopTime,  ///< [in] stop time if AtTime
       std::optional<TimeInterval> InDuration  ///< [in] duration if AfterDur
   );

   /// Destructor frees the SUNDIALS integrator, vectors and table
   ~SundialsStepper() override;

   /// Advance the state by one step of the chosen Runge-Kutta method
   void doStep(OceanState *State,   ///< [inout] model state
               TimeInstant &SimTime ///< [inout] current simulation time
   ) const override;

   /// Name of the method, the SUNDIALS table name or Custom
   const std::string &getMethodName() const;

   /// Number of stages of the method
   int getNumStages() const;

   /// Order of accuracy of the method
   int getOrder() const;

   // These are called by ARKODE during a step. They should be private, they
   // are public only because of CUDA limitations.

   /// Prepares an intermediate stage, which ARKODE has formed in the next
   /// time level: prescribes the state, converts the thickness-weighted
   /// tracers to concentrations and exchanges halos
   void prepareStage(int Stage ///< [in] stage index
   ) const;

   /// Computes the tendencies at a stage and copies them to YDot
   void computeStageTendencies(int Stage,    ///< [in] stage index
                               N_Vector YDot ///< [out] SUNDIALS tendencies
   ) const;

   /// Copies the tendencies into the SUNDIALS vector YDot, with zeros in
   /// inactive layers and padding so that ARKODE leaves those unchanged
   void packTendencies(N_Vector YDot ///< [out] SUNDIALS tendencies
   ) const;

   /// Index of the stage ARKODE is working on
   int getStageIndex() const;

   /// The vector ARKODE forms intermediate stages and the result in
   N_Vector getStepEnd() const;

 protected:
   /// Allocates the stage tracer array
   void finalizeInit() override;

 private:
   /// Name of the method
   std::string MethodName;

   /// Butcher table of the method
   ARKodeButcherTable Table = nullptr;

   /// Fraction of the time step at which each stage is evaluated (the c
   /// coefficients of the Butcher table)
   std::vector<Real> StageFrac;

   /// Fraction of the time step over which each stage projects the
   /// pseudo-thickness: the next stage's StageFrac, and 1 for the last stage
   std::vector<Real> ProjFrac;

   /// Tracer concentrations at an intermediate stage
   Array3DReal StageTracers;

   /// SUNDIALS context, created with the integrator on the first step
   mutable SUNContext Ctx = nullptr;

   /// ARKODE integrator memory
   mutable void *ArkMem = nullptr;

   /// State and start time of the step in progress, and the vector ARKODE
   /// forms the stages and result in
   mutable OceanState *StepState = nullptr;
   mutable TimeInstant StepStart;
   mutable N_Vector StepEnd = nullptr;

   /// SUNDIALS vectors wrapping Omega arrays, by data pointer and size. Omega
   /// rotates time levels without moving arrays, so these are created on the
   /// first steps and reused.
   mutable std::map<std::pair<const Real *, size_t>,
                    std::unique_ptr<SundialsVector>>
       ArrayVectors;

   /// State vectors {thickness, velocity, thickness-weighted tracers}, by
   /// the data pointers of their three arrays
   mutable std::map<std::array<const Real *, 3>, N_Vector> StateVectors;

   /// Loads the Butcher table chosen in the configuration
   void loadTable();

   /// Returns the state vector made of three Omega arrays
   N_Vector getStateVector(const Array2DReal &Thickness,
                           const Array2DReal &Velocity,
                           const Array3DReal &WeightedTracers) const;

   /// Creates the SUNDIALS context and integrator on the first step
   void createIntegrator(N_Vector Template) const;
};

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif

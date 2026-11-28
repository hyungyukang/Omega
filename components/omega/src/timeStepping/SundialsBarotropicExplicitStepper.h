#ifndef OMEGA_SUNDIALS_BAROTROPIC_EXPLICIT_STEPPER_H
#define OMEGA_SUNDIALS_BAROTROPIC_EXPLICIT_STEPPER_H
//===-- SundialsBarotropicExplicitStepper.h - ERK barotropic --*- C++ -*-===//
//
/// \file
/// \brief Explicit Runge-Kutta barotropic stepper built on SUNDIALS ERKStep
///
/// The barotropic system is integrated over the window with NSubsteps fixed
/// steps of an explicit Runge-Kutta method of the SUNDIALS ARKODE ERKStep
/// integrator (ExplicitMethod RK2, RK3, RK4 or any SUNDIALS explicit Butcher
/// table). The ODE state is
///    { ubar, B', PhiF, PhiU }
/// where PhiF and PhiU integrate the barotropic flux P(B') ubar and the
/// velocity ubar in time (d PhiF/dt = P ubar, d PhiU/dt = ubar). At the end
/// of the window PhiF / Window and PhiU / Window are the time-averaged flux
/// and velocity. Because the flux accumulator is integrated with the same
/// Runge-Kutta weights as B', B'(Window) = B'(0) + Window * (-div FluxAvg)
/// holds exactly, which is the consistency the transport-velocity correction
/// relies on.
///
/// Each piece is an Omega array wrapped as a SUNDIALS Kokkos vector and the
/// pieces are joined in an MPI ManyVector, as in the Sundials time stepper.
/// The right-hand side exchanges the ubar and B' halos once per stage and
/// computes the tendencies on owned edges and cells.
//
//===----------------------------------------------------------------------===//

#include "BarotropicStepper.h"
#include "SundialsVector.h"

#include <arkode/arkode_butcher.h>
#include <sundials/sundials_context.h>

namespace OMEGA {

class SundialsBarotropicExplicitStepper : public BarotropicStepper {
 public:
   SundialsBarotropicExplicitStepper(
       const std::string &InName,          ///< [in] name of the owner
       const BarotropicStepperConfig &Cfg, ///< [in] options
       const AuxiliaryState *InAuxState,   ///< [in] provisional aux state
       const HorzMesh *InMesh,             ///< [in] horizontal mesh
       Halo *InMeshHalo,                   ///< [in] mesh halo exchange
       const VertCoord *InVCoord           ///< [in] vertical coordinate
   );

   ~SundialsBarotropicExplicitStepper() override;

   /// Integrates the barotropic system over the window of In
   void advance(const BarotropicInput &In,  ///< [in] inputs
                const BarotropicOutput &Out ///< [out] outputs
                ) override;

   /// Number of stages and order of the Runge-Kutta method
   int getNumStages() const;
   int getOrder() const;

   /// Right-hand side, called by ARKODE through a static wrapper. Public only
   /// because of CUDA limitations on lambdas in private member functions.
   void computeRHS(N_Vector Y,   ///< [in] state {ubar, B', PhiF, PhiU}
                   N_Vector YDot ///< [out] tendencies
   );

 private:
   /// Butcher table of the method
   ARKodeButcherTable Table = nullptr;

   /// SUNDIALS context and ERKStep memory, created on the first advance
   SUNContext Ctx = nullptr;
   void *ArkMem   = nullptr;

   /// ODE state arrays and the wrapped state vector
   Array1DReal Vel;
   Array1DReal Press;
   Array1DReal FluxIntegral;
   Array1DReal VelIntegral;
   std::unique_ptr<SundialsVector> VelVec;
   std::unique_ptr<SundialsVector> PressVec;
   std::unique_ptr<SundialsVector> FluxVec;
   std::unique_ptr<SundialsVector> VelIntVec;
   N_Vector StateVec = nullptr;

   /// Inputs held fixed over the window, captured by the right-hand side
   Array1DReal PressRef; ///< copy of In.PressRef
   Array1DReal Forcing;  ///< In.Forcing

   /// Creates the context, the wrapped state vector and the integrator
   void createIntegrator();
};

} // namespace OMEGA
#endif

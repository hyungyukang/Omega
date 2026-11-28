#ifndef OMEGA_SUNDIALS_BAROTROPIC_IMPLICIT_STEPPER_H
#define OMEGA_SUNDIALS_BAROTROPIC_IMPLICIT_STEPPER_H
//===-- SundialsBarotropicImplicitStepper.h - theta barotropic *- C++ -*-===//
//
/// \file
/// \brief Implicit (theta / Crank-Nicolson) barotropic stepper solved with a
/// SUNDIALS Krylov solver
///
/// The barotropic system is advanced over the window W in one step of the
/// theta method (Theta = 0.5 is Crank-Nicolson),
///    ubar^{n+1} = utilde - Theta W alphabar grad B'^{n+1}
///    utilde     = ubar^n + W [ Cor - (1-Theta) alphabar grad B'^n + Gbar ]
///    B'^{n+1}   = B'^n - W div( P [Theta ubar^{n+1} + (1-Theta) ubar^n] )
/// Substituting the first equation into the last gives a symmetric positive
/// definite Helmholtz equation for the pseudo thickness, written for B',
///    A_i B'^{n+1}_i + Theta^2 W^2 sum_e s_ie w_e (B'^{n+1}_c1 - B'^{n+1}_c0)
///       = A_i [ B'^n - W div( P (Theta utilde + (1-Theta) ubar^n) ) ]_i
/// with w_e = dv_e P_e alphabar_e / dc_e and A_i the cell area, which is
/// solved with the SUNDIALS PCG or BiCGStab solver, used without an ARKODE
/// integrator. The barotropic velocity is then computed explicitly from the
/// new pseudo thickness, the flux P (Theta ubar^{n+1} + (1-Theta) ubar^n) is
/// formed from the final velocity and B'^{n+1} is recomputed from its
/// divergence, so that mass is conserved and the flux returned for the
/// transport-velocity correction matches the B' update exactly, whatever the
/// solver tolerance.
///
/// The flux thickness P is held at its value at the start of the window and
/// the Coriolis term is explicit, evaluated with the average of the start
/// velocity and the provisional velocity VelRef. In the split-explicit
/// stepper VelRef is the start velocity in the first outer iteration and the
/// first iteration's result in the second, which makes the Coriolis term
/// trapezoidal over the two iterations.
///
/// The Krylov vectors hold the owned cells only, so their inner products are
/// global without double counting halo cells; the operator exchanges halos
/// of its argument in a work array. The solver is given the residual of the
/// initial guess B'^n and solves for the correction from a zero guess, since
/// right-preconditioned BiCGStab in SUNDIALS needs a zero initial guess.
//
//===----------------------------------------------------------------------===//

#include "BarotropicStepper.h"
#include "SundialsVector.h"

#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>

namespace OMEGA {

class SundialsBarotropicImplicitStepper : public BarotropicStepper {
 public:
   SundialsBarotropicImplicitStepper(
       const std::string &InName,          ///< [in] name of the owner
       const BarotropicStepperConfig &Cfg, ///< [in] options
       const AuxiliaryState *InAuxState,   ///< [in] provisional aux state
       const HorzMesh *InMesh,             ///< [in] horizontal mesh
       Halo *InMeshHalo,                   ///< [in] mesh halo exchange
       const VertCoord *InVCoord           ///< [in] vertical coordinate
   );

   ~SundialsBarotropicImplicitStepper() override;

   /// Advances the barotropic system over the window of In in one step
   void advance(const BarotropicInput &In,  ///< [in] inputs
                const BarotropicOutput &Out ///< [out] outputs
                ) override;

   bool isImplicit() const override { return true; }

   /// Diagnostics of the last solve
   int getLastNumIterations() const { return LastNumIterations; }
   Real getLastResidualNorm() const { return LastResidualNorm; }
   Real getLastTolerance() const { return LastTolerance; }

   // Called by SUNDIALS through static wrappers. Public only because of CUDA
   // limitations on lambdas in private member functions.

   /// Z = A V for the Helmholtz operator
   void applyOperator(N_Vector V, ///< [in] owned-cell vector
                      N_Vector Z  ///< [out] owned-cell vector
   );

   /// Z = D^{-1} R with D the diagonal of the operator
   void applyJacobi(N_Vector R, ///< [in] owned-cell vector
                    N_Vector Z  ///< [out] owned-cell vector
   );

 private:
   BarotropicStepperConfig Config;

   /// SUNDIALS context, Krylov solver and owned-cell vectors
   SUNContext Ctx     = nullptr;
   SUNLinearSolver LS = nullptr;
   std::unique_ptr<SundialsVector> SolPiece;
   std::unique_ptr<SundialsVector> RhsPiece;
   N_Vector SolVec  = nullptr;
   N_Vector RhsVec  = nullptr;
   N_Vector ResVec  = nullptr; ///< residual of the initial guess
   N_Vector CorrVec = nullptr; ///< correction to the initial guess

   /// Edge work arrays: flux thickness times rho0 g, operator weights,
   /// explicit part of the new velocity and edge mean specific volume
   Array1DReal FluxPressEdge;
   Array1DReal OperatorWeight;
   Array1DReal VelTilde;
   Array1DReal MeanSpecVolEdge;

   /// Cell work arrays: operator argument with halo, new pseudo thickness
   /// with halo and the operator diagonal
   Array1DReal OperatorArg;
   Array1DReal PressNew;
   Array1DReal JacobiDiag;

   /// Theta^2 W^2 of the current solve
   Real OperatorCoef = 0;

   int LastNumIterations = 0;
   Real LastResidualNorm = 0;
   Real LastTolerance    = 0;

   /// Creates the context, the vectors and the Krylov solver
   void createSolver();
};

} // namespace OMEGA
#endif

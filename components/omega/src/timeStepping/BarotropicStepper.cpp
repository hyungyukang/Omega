//===-- BarotropicStepper.cpp - barotropic split steppers -----*- C++ -*-===//
//
// Options, creation and shared utilities of the barotropic steppers used by
// split time stepping.
//
//===----------------------------------------------------------------------===//

#include "BarotropicStepper.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "Error.h"
#include "Logging.h"
#include "SplitExplicitBarotropicPCStepper.h"
#include "SundialsBarotropicExplicitStepper.h"
#include "SundialsBarotropicImplicitStepper.h"

namespace OMEGA {

namespace {

BarotropicLinearSolver getLinearSolverFromStr(const std::string &Str) {
   if (Str == "PCG")
      return BarotropicLinearSolver::PCG;
   if (Str == "BiCGStab")
      return BarotropicLinearSolver::BiCGStab;
   ABORT_ERROR("ModeSplitShare::Barotropic::LinearSolver should be 'PCG' or "
               "'BiCGStab' but got {}",
               Str);
   return BarotropicLinearSolver::PCG;
}

BarotropicPreconditioner getPreconditionerFromStr(const std::string &Str) {
   if (Str == "None")
      return BarotropicPreconditioner::None;
   if (Str == "Jacobi")
      return BarotropicPreconditioner::Jacobi;
   ABORT_ERROR("ModeSplitShare::Barotropic::Preconditioner should be 'None' "
               "or 'Jacobi' but got {}",
               Str);
   return BarotropicPreconditioner::None;
}

} // namespace

//------------------------------------------------------------------------------
// Reads the barotropic stepper options for the given BtrTimeStepper. The
// Barotropic subgroup is optional; missing keys keep their defaults.
BarotropicStepperConfig
BarotropicStepper::readConfig(const std::string &BtrTimeStepperStr) {

   BarotropicStepperConfig Cfg;

   Config BtrConfig("Barotropic");
   bool HasBtrConfig   = false;
   Config *OmegaConfig = Config::getOmegaConfig();
   if (OmegaConfig) {
      Config TimeIntConfig("TimeIntegration");
      Config ModeSplitShareConfig("ModeSplitShare");
      if (OmegaConfig->get(TimeIntConfig).isSuccess() &&
          TimeIntConfig.get(ModeSplitShareConfig).isSuccess() &&
          ModeSplitShareConfig.get(BtrConfig).isSuccess())
         HasBtrConfig = true;
   }

   if (BtrTimeStepperStr == "Predictor-Corrector") {
      // The original name of the Omega forward-backward predictor-corrector
      Cfg.Method         = BarotropicMethod::ForwardBackward;
      Cfg.ExplicitMethod = "ForwardBackward";
      std::string Method;
      if (HasBtrConfig && BtrConfig.get("ExplicitMethod", Method).isSuccess() &&
          Method != "ForwardBackward")
         LOG_WARN("BtrTimeStepper Predictor-Corrector ignores "
                  "Barotropic::ExplicitMethod {}; use BtrTimeStepper "
                  "Explicit to select it",
                  Method);
   } else if (BtrTimeStepperStr == "Explicit") {
      if (HasBtrConfig)
         BtrConfig.get("ExplicitMethod", Cfg.ExplicitMethod);
      Cfg.Method = Cfg.ExplicitMethod == "ForwardBackward"
                       ? BarotropicMethod::ForwardBackward
                       : BarotropicMethod::SundialsExplicit;
   } else if (BtrTimeStepperStr == "Implicit") {
      Cfg.Method = BarotropicMethod::SundialsImplicit;
   } else {
      ABORT_ERROR("BtrTimeStepper should be 'Predictor-Corrector', 'Explicit' "
                  "or 'Implicit' but got {}",
                  BtrTimeStepperStr);
   }

   if (Cfg.Method == BarotropicMethod::SundialsImplicit && HasBtrConfig) {
      R8 Value = 0;
      if (BtrConfig.get("Theta", Value).isSuccess())
         Cfg.Theta = Value;
      std::string Str;
      if (BtrConfig.get("LinearSolver", Str).isSuccess())
         Cfg.LinearSolver = getLinearSolverFromStr(Str);
      if (BtrConfig.get("Preconditioner", Str).isSuccess())
         Cfg.Preconditioner = getPreconditionerFromStr(Str);
      if (BtrConfig.get("RelTolerance", Value).isSuccess())
         Cfg.RelTolerance = Value;
      if (BtrConfig.get("AbsTolerance", Value).isSuccess())
         Cfg.AbsTolerance = Value;
      BtrConfig.get("MaxIterations", Cfg.MaxIterations);
   }

   if (Cfg.Theta < 0.5_Real || Cfg.Theta > 1._Real)
      ABORT_ERROR("ModeSplitShare::Barotropic::Theta must be in [0.5, 1] for "
                  "an unconditionally stable implicit step, got {}",
                  Cfg.Theta);
   if (Cfg.RelTolerance < 0 || Cfg.AbsTolerance < 0 ||
       Cfg.RelTolerance + Cfg.AbsTolerance <= 0)
      ABORT_ERROR("ModeSplitShare::Barotropic tolerances must be "
                  "non-negative and not both zero, got RelTolerance={} "
                  "AbsTolerance={}",
                  Cfg.RelTolerance, Cfg.AbsTolerance);
   if (Cfg.MaxIterations < 1)
      ABORT_ERROR("ModeSplitShare::Barotropic::MaxIterations must be "
                  "positive, got {}",
                  Cfg.MaxIterations);

   return Cfg;
}

//------------------------------------------------------------------------------
std::unique_ptr<BarotropicStepper>
BarotropicStepper::create(const BarotropicStepperConfig &Cfg,
                          const std::string &Name,
                          const AuxiliaryState *AuxState, const HorzMesh *Mesh,
                          Halo *MeshHalo, const VertCoord *VCoord) {

   if (!AuxState)
      LOG_CRITICAL("Invalid auxiliary state");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!MeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");

   switch (Cfg.Method) {
   case BarotropicMethod::ForwardBackward:
      return std::make_unique<SplitExplicitBarotropicPCStepper>(
          Name, AuxState, Mesh, MeshHalo, VCoord);
   case BarotropicMethod::SundialsExplicit:
      return std::make_unique<SundialsBarotropicExplicitStepper>(
          Name, Cfg, AuxState, Mesh, MeshHalo, VCoord);
   case BarotropicMethod::SundialsImplicit:
      return std::make_unique<SundialsBarotropicImplicitStepper>(
          Name, Cfg, AuxState, Mesh, MeshHalo, VCoord);
   default:
      ABORT_ERROR("Invalid barotropic stepper for {}", Name);
   }
   return nullptr;
}

//------------------------------------------------------------------------------
BarotropicStepper::BarotropicStepper(const std::string &InName,
                                     const std::string &InMethod,
                                     const AuxiliaryState *InAuxState,
                                     const HorzMesh *InMesh, Halo *InMeshHalo,
                                     const VertCoord *InVCoord)
    : Name(InName), MethodName(InMethod), AuxState(InAuxState), Mesh(InMesh),
      MeshHalo(InMeshHalo), VCoord(InVCoord) {

   ColumnFluxThickEdge =
       Array1DReal("BtrColumnFluxThickEdge" + Name, Mesh->NEdgesSize);
   deepCopy(ColumnFluxThickEdge, 0._Real);
}

//------------------------------------------------------------------------------
void BarotropicStepper::getDefaultWindow(Real TimeStepSeconds, I4 NBtrSubcycles,
                                         Real &Window, I4 &NSubsteps) const {
   if (isImplicit()) {
      Window    = TimeStepSeconds;
      NSubsteps = 1;
   } else {
      Window    = 2._Real * TimeStepSeconds;
      NSubsteps = 2 * NBtrSubcycles;
   }
}

//------------------------------------------------------------------------------
bool BarotropicStepper::upwindFluxThickness() const {
   return AuxState->PseudoThicknessAux.FluxThickEdgeChoice ==
          FluxThickEdgeOption::Upwind;
}

//------------------------------------------------------------------------------
void BarotropicStepper::computeColumnFluxThickness() const {

   const Array2DReal FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   Array1DReal LocColumnThick = ColumnFluxThickEdge;

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "computeBtrColumnFluxThickness", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          Real PseudoThicknessSum = 0._Real;
          parallelReduceInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(I4 K, Real & Sum) {
                 Sum += FluxPseudoThickEdge(IEdge, K);
              },
              PseudoThicknessSum);

          Kokkos::single(
              PerTeam(Team),
              INNER_LAMBDA() { LocColumnThick(IEdge) = PseudoThicknessSum; });
       });
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//

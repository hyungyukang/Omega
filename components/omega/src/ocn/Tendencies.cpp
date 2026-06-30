//===-- ocn/Tendencies.cpp - Tendencies ------------------*- C++ -*-===//
//
// The Tendencies class is responsible for managing tendencies of state
// variables and tracers. It contains arrays that store the tendency data
// and provides methods for computing different tendency groups for use
// within the timestepping algorithm. At initialization, it determines which
// tendency terms are enabled.
//
//===----------------------------------------------------------------------===//

#include "Tendencies.h"
#include "CustomTendencyTerms.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
#include "HorzOperators.h"
#include "KPPMix.h"
#include "OceanState.h"
#include "PGrad.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertMix.h"
#include <string>

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh, VertCoord, VertAdv, and
// TimeStepper  has already been initialized.
void Tendencies::init() {
   Error Err; // error code

   HorzMesh *DefHorzMesh       = HorzMesh::getDefault();
   VertCoord *DefVertCoord     = VertCoord::getDefault();
   VertAdv *DefVertAdv         = VertAdv::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   Eos *DefEos                 = Eos::getInstance();
   PressureGrad *DefPGrad      = PressureGrad::getDefault();
   VertMix *DefVertMix         = VertMix::getInstance();

   I4 NTracers = Tracers::getNumTracers();

   // Get TendConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   // Check if use the customized tendencies. If it is not found in the
   // config, we assume it is not used (false)
   bool UseCustomTendency = false;
   Err += TendConfig.get("UseCustomTendency", UseCustomTendency);

   /// Instances of custom tendencies - empty by default
   CustomTendencyType CustomThickTend;
   CustomTendencyType CustomVelTend;

   if (UseCustomTendency) {
      // Check if use manufactured tendency terms if it is not found in
      // the config file, we will assume it is not used (false)
      bool ManufacturedTend = false;
      Error ManufacturedTendErr =
          TendConfig.get("ManufacturedSolutionTendency", ManufacturedTend);

      if (ManufacturedTend) {
         ManufacturedSolution ManufacturedSol;
         ManufacturedSol.init();

         CustomThickTend = ManufacturedSol.ManufacturedThickTend;
         CustomVelTend   = ManufacturedSol.ManufacturedVelTend;

      } // if ManufacturedTend

   } // end if UseCustomTendency

   TimeInterval TimeStep = DefTimeStepper->getTimeStep();

   // Ceate default tendencies
   Tendencies::DefaultTendencies =
       create("Default", DefHorzMesh, DefVertCoord, DefVertAdv, DefPGrad,
              DefEos, DefVertMix, NTracers, TimeStep, &TendConfig,
              CustomThickTend, CustomVelTend);

   DefaultTendencies->readConfig(OmegaConfig);

   if (DefaultTendencies->TracerHorzAdv.Enabled)
      DefaultTendencies->TracerHorzAdv.init();
} // end init

//------------------------------------------------------------------------------
// Destroys the tendencies
Tendencies::~Tendencies() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Removes all tendencies instances before exit
void Tendencies::clear() {
   AllTendencies.clear();
   DefaultTendencies = nullptr; // prevent dangling pointer
} // end clear

//------------------------------------------------------------------------------
// Removes tendencies from list by name
void Tendencies::erase(const std::string &Name) {

   AllTendencies.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get default tendencies
Tendencies *Tendencies::getDefault() {

   return Tendencies::DefaultTendencies;

} // end get default

//------------------------------------------------------------------------------
// Get tendencies by name
Tendencies *Tendencies::get(const std::string &Name ///< [in] Name of tendencies
) {

   auto it = AllTendencies.find(Name);

   if (it != AllTendencies.end()) {
      return it->second.get();
   } else {
      LOG_ERROR(
          "Tendencies::get: Attempt to retrieve non-existent tendencies:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get tendencies

//------------------------------------------------------------------------------
// Check if tendencies exists by name
bool Tendencies::exists(const std::string &Name ///< [in] Name of tendencies
) {

   return AllTendencies.find(Name) != AllTendencies.end();

} // end exists

//------------------------------------------------------------------------------
// read and set config options
void Tendencies::readConfig(Config *OmegaConfig ///< [in] Omega config
) {
   Error Err; // error code

   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   Err += TendConfig.get("ThicknessFluxTendencyEnable",
                         this->PseudoThicknessFluxDiv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: ThicknessFluxTendencyEnable not found in TendConfig");

   Err += TendConfig.get("PVTendencyEnable", this->PotentialVortHAdv.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: PVTendencyEnable not found in TendConfig");

   Err += TendConfig.get("KETendencyEnable", this->KEGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: KETendencyEnable not found in TendConfig");

   Err += TendConfig.get("SSHTendencyEnable", this->SSHGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: SSHTendencyEnable not found in TendConfig");

   Err +=
       TendConfig.get("VelDiffTendencyEnable", this->VelocityDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelDiffTendencyEnable not found in TendConfig");

   Err += TendConfig.get("VelHyperDiffTendencyEnable",
                         this->VelocityHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelHyperDiffTendencyEnable not found in TendConfig");

   if (this->VelocityDiffusion.Enabled) {
      Err += TendConfig.get("ViscDel2", this->VelocityDiffusion.ViscDel2);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel2 not found in TendConfig");
   }

   if (this->VelocityHyperDiff.Enabled) {
      Err += TendConfig.get("ViscDel4", this->VelocityHyperDiff.ViscDel4);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel4 not found in TendConfig");
      Err += TendConfig.get("DivFactor", this->VelocityHyperDiff.DivFactor);
      CHECK_ERROR_ABORT(Err, "Tendencies: DivFactor not found in TendConfig");
   }
   Err += TendConfig.get("TracerHorzAdvTendencyEnable",
                         this->TracerHorzAdv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerHorzAdvTendencyEnable not found in TendConfig");
   if (this->TracerHorzAdv.Enabled) {
      Config AdvectConfig("Advection");
      Err += OmegaConfig->get(AdvectConfig);
      CHECK_ERROR_ABORT(Err, "Tendencies: Advection group not in Config");

      I4 Order = 0;
      Err += AdvectConfig.get("HorzTracerFluxOrder", Order);
      CHECK_ERROR_ABORT(
          Err, "Tendencies: HorzTracerFluxOrder not found in AdvectConfig");
      OMEGA_REQUIRE(Order >= 2 && Order <= 4,
                    "HorzTracerFluxOrder: Only values are 2, 3, 4, found {}",
                    Order);

      if (Order == 2) {
         this->TracerHorzAdv.ForceLowOrder = true;
         this->TracerHorzAdv.Coef3rdOrder  = 0;
      }
      if (Order == 3) {
         Err += AdvectConfig.get("Coef3rdOrder", TracerHorzAdv.Coef3rdOrder);
         CHECK_ERROR_ABORT(
             Err, "Tendencies: Coef3rdOrder not found in AdvectConfig");
      }
      if (Order == 4) {
         this->TracerHorzAdv.Coef3rdOrder = 0;
      }
   }
   Err += TendConfig.get("TracerDiffTendencyEnable",
                         this->TracerDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerDiffTendencyEnable not found in TendConfig");

   Err += TendConfig.get("SfcStressForcingTendencyEnable",
                         this->SfcStressForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: SfcStressForcingTendencyEnable not found in TendConfig");

   Err += TendConfig.get("BottomDragTendencyEnable", this->BottomDrag.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendencyEnable not found in TendConfig");

   Err += TendConfig.get("BottomDragCoeff", this->BottomDrag.Coeff);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: BottomDragCoeff not found in TendConfig");

   if (this->TracerDiffusion.Enabled) {
      Err += TendConfig.get("EddyDiff2", this->TracerDiffusion.EddyDiff2);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff2 not found in TendConfig");
   }

   Err += TendConfig.get("TracerHyperDiffTendencyEnable",
                         this->TracerHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: TracerHyperDiffTendencyEnable not found in TendConfig");

   if (this->TracerHyperDiff.Enabled) {
      Err += TendConfig.get("EddyDiff4", this->TracerHyperDiff.EddyDiff4);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff4 not found in TendConfig");
   }

   Err += TendConfig.get("PressureGradTendencyEnable", this->PGrad->Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: PressureGradTendencyEnable not found in TendConfig");

   Err += TendConfig.get("SurfaceTracerRestoringEnable",
                         this->SurfaceTracerRestoring.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: SurfaceTracerRestoringEnable not found in TendConfig");
   if (this->SurfaceTracerRestoring.Enabled) {
      Config SurfRestConfig("SurfaceRestoring");
      Err += OmegaConfig->get(SurfRestConfig);
      Err += SurfRestConfig.get("PistonVelocity",
                                this->SurfaceTracerRestoring.PistonVelocity);
      CHECK_ERROR_ABORT(
          Err,
          "Tendencies: PistonVelocity not found in SurfaceRestoringConfig");

      std::vector<std::string> TracersToRestore;
      SurfRestConfig.get("TracersToRestore", TracersToRestore);

      // Enable restoring for specified individual tracers
      I4 NumInvalidTracers = 0;
      std::vector<I4> TracerIdsToRestoreVec;
      for (const auto &TracerName : TracersToRestore) {
         I4 TracerIndex = -1;
         Tracers::getIndex(TracerIndex, TracerName);
         if (TracerIndex == -1) {
            LOG_ERROR("Tendencies: Tracer {} in TracersToRestore is not "
                      "defined",
                      TracerName);
            NumInvalidTracers++;
         } else {
            bool IsDuplicate = false;
            for (const auto ExistingTracerIndex : TracerIdsToRestoreVec) {
               if (ExistingTracerIndex == TracerIndex) {
                  IsDuplicate = true;
                  break;
               }
            }
            if (!IsDuplicate) {
               TracerIdsToRestoreVec.push_back(TracerIndex);
            }
         }
      }
      if (NumInvalidTracers > 0) {
         ABORT_ERROR("Tendencies: {} invalid tracer(s) in TracersToRestore "
                     "configuration",
                     NumInvalidTracers);
      }

      this->SurfaceTracerRestoring.NTracersToRestore =
          static_cast<I4>(TracerIdsToRestoreVec.size());
      this->SurfaceTracerRestoring.TracerIdsToRestore = Array1DI4(
          "TracerIdsToRestore", this->SurfaceTracerRestoring.NTracersToRestore);
      deepCopy(this->SurfaceTracerRestoring.TracerIdsToRestore,
               HostArray1DI4(TracerIdsToRestoreVec.data(),
                             TracerIdsToRestoreVec.size()));
   }

   // Validate VertMix tendency
   Err += TendConfig.get("VelVertMixTendencyEnable",
                         this->VMix->VelVertMixSetup.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelVertMixTendencyEnable not found in TendConfig");

   Err += TendConfig.get("TracerVertMixTendencyEnable",
                         this->VMix->TracerVertMixSetup.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerVertMixTendencyEnable not found in TendConfig");

   if (this->VMix->VelVertMixSetup.Enabled ||
       this->VMix->TracerVertMixSetup.Enabled) {

      if (!this->EqState) {
         ABORT_ERROR("Tendencies: Eos must be initialized when"
                     "vertical mixing tendencies are enabled");
      }
      if (!this->VMix) {
         ABORT_ERROR("Tendencies: VertMix must be initialized when"
                     "vertical mixing tendencies are enabled");
      }
      // Optional KPP non-local tracer tendency: no abort if missing
      Error TracerNonLocalErr = TendConfig.get(
          "TracerNonLocalFluxTendencyEnable", this->TracerNonLocalFluxEnabled);
      if (!TracerNonLocalErr.isSuccess()) {
         TracerNonLocalErr.reset();
         this->TracerNonLocalFluxEnabled = false;
      }

      Error TempFluxBridgeErr =
          TendConfig.get("UseTempSurfaceTracerFluxBridge",
                         this->UseTempSurfaceTracerFluxBridge);
      if (!TempFluxBridgeErr.isSuccess()) {
         TempFluxBridgeErr.reset();
         this->UseTempSurfaceTracerFluxBridge = true;
      }

      Error TempTopLayerBridgeErr =
          TendConfig.get("UseTempTopLayerFluxTendencyBridge",
                         this->UseTempTopLayerFluxTendencyBridge);
      if (!TempTopLayerBridgeErr.isSuccess()) {
         TempTopLayerBridgeErr.reset();
         this->UseTempTopLayerFluxTendencyBridge = true;
      }

      // Conversion factors: heat flux to buoyancy, freshwater flux to buoyancy.
      // Fall back to physically-derived defaults if not in config.
      const Real HFluxFac  = 1._Real / (RhoSw * CpSw);
      const Real FwFluxFac = 1._Real / RhoFw;
      Error HeatFacErr     = TendConfig.get("KPPHeatFluxToBuoyancyFactor",
                                            this->KPPHeatFluxToBuoyancyFactor);
      if (!HeatFacErr.isSuccess()) {
         HeatFacErr.reset();
         this->KPPHeatFluxToBuoyancyFactor = Gravity * 2.0e-4_Real * HFluxFac;
      }
      Error ThickFacErr =
          TendConfig.get("KPPThicknessFluxToBuoyancyFactor",
                         this->KPPThicknessFluxToBuoyancyFactor);
      // FIXME - This has a magic number that needs to go away
      if (!ThickFacErr.isSuccess()) {
         ThickFacErr.reset();
         this->KPPThicknessFluxToBuoyancyFactor =
             -Gravity * 8.0e-4_Real * OcnRefSal * FwFluxFac;
      }
   }
}

//------------------------------------------------------------------------------
// Define fields associated with tendencies
void Tendencies::defineFields() {
   std::string PseudoThicknessTendFieldName = "PseudoThicknessTend";
   std::string NormalVelocityTendFieldName  = "NormalVelocityTend";
   std::string TracerTendFieldName          = "TracerTend";
   std::string SurfaceTracerFluxFieldName   = "SurfaceTracerFlux";
   if (Name != "Default") {
      PseudoThicknessTendFieldName.append(Name);
      NormalVelocityTendFieldName.append(Name);
      TracerTendFieldName.append(Name);
      SurfaceTracerFluxFieldName.append(Name);
   }

   int NDims = 2;
   std::vector<std::string> DimNamesThickness(NDims);
   DimNamesThickness[0] = "NCells";
   DimNamesThickness[1] = "NVertLayers";
   auto PseudoThicknessTendField =
       Field::create(PseudoThicknessTendFieldName, "Pseudo-thickness tendency",
                     "m/s", "cell_thickness_tendency", -9.99E+10, 9.99E+10,
                     -9.99E+30, NDims, DimNamesThickness);
   NDims = 3;
   std::vector<std::string> DimNamesTracer(NDims);
   DimNamesTracer[0]    = "NTracers";
   DimNamesTracer[1]    = "NCells";
   DimNamesTracer[2]    = "NVertLayers";
   auto TracerTendField = Field::create(
       TracerTendFieldName, "Tracer tendency", "kg/m^3/s", "tracer_tendency",
       -9.99E+10, 9.99E+10, -9.99E+30, NDims, DimNamesTracer);
   NDims = 2;
   std::vector<std::string> DimNamesVelocity(NDims);
   DimNamesVelocity[0] = "NEdges";
   DimNamesVelocity[1] = "NVertLayers";
   auto NormalVelocityTendField =
       Field::create(NormalVelocityTendFieldName, "Normal velocity tendency",
                     "m/s^2", "sea_water_velocity_tendency", -9.99E+10,
                     9.99E+10, -9.99E+30, NDims, DimNamesVelocity);

   NDims = 2;
   std::vector<std::string> DimNamesSurfaceFlux(NDims);
   DimNamesSurfaceFlux[0]      = "NTracers";
   DimNamesSurfaceFlux[1]      = "NCells";
   auto SurfaceTracerFluxField = Field::create(
       SurfaceTracerFluxFieldName, "Surface tracer flux", "1", "", -9.99E+10,
       9.99E+10, -9.99E+30, NDims, DimNamesSurfaceFlux);

   std::string TendGroupName = "Tendencies";
   if (Name != "Default") {
      TendGroupName.append(Name);
   }
   auto TendGroup = FieldGroup::create(TendGroupName);

   TendGroup->addField(PseudoThicknessTendFieldName);
   TendGroup->addField(NormalVelocityTendFieldName);
   TendGroup->addField(TracerTendFieldName);
   TendGroup->addField(SurfaceTracerFluxFieldName);

   PseudoThicknessTendField->attachData<Array2DReal>(PseudoThicknessTend);
   NormalVelocityTendField->attachData<Array2DReal>(NormalVelocityTend);
   TracerTendField->attachData<Array3DReal>(TracerTend);
   SurfaceTracerFluxField->attachData<Array2DReal>(SurfaceTracerFlux);

} // end defineFields

//------------------------------------------------------------------------------
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name_, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,     ///< [in] Horizontal mesh
                       VertCoord *VCoord,        ///< [in] Vertical coordinate
                       VertAdv *VAdv,            ///< [in] Vertical advection
                       PressureGrad *PGrad,      ///< [in] Pressure gradient
                       Eos *EqState,             ///< [in] Equation of state
                       VertMix *VMix,            ///< [in] Vertical mixing
                       int NTracersIn,           ///< [in] Number of tracers
                       TimeInterval TimeStepIn,  ///< [in] Time step
                       Config *Options,          ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : Mesh(Mesh), VCoord(VCoord), VAdv(VAdv),
      PseudoThicknessFluxDiv(Mesh, VCoord), PotentialVortHAdv(Mesh, VCoord),
      KEGrad(Mesh, VCoord), SSHGrad(Mesh, VCoord),
      VelocityDiffusion(Mesh, VCoord), VelocityHyperDiff(Mesh, VCoord),
      SfcStressForcing(Mesh, VCoord), BottomDrag(Mesh, VCoord),
      TracerDiffusion(Mesh, VCoord), TracerHyperDiff(Mesh, VCoord),
      TracerHorzAdv(Mesh, VCoord), SurfaceTracerRestoring(Mesh),
      PresGradZ(Mesh, VCoord), PresGradForce(Mesh, VCoord),
      GeoptGrad(Mesh, VCoord), CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend), EqState(EqState), PGrad(PGrad),
      VMix(VMix) {

   // Tendency arrays
   PseudoThicknessTend = Array2DReal("PseudoThicknessTend", Mesh->NCellsSize,
                                     VCoord->NVertLayers);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, VCoord->NVertLayers);
   TracerTend = Array3DReal("TracerTend", NTracersIn, Mesh->NCellsSize,
                            VCoord->NVertLayers);
   SurfaceTracerFlux =
       Array2DReal("SurfaceTracerFlux", NTracersIn, Mesh->NCellsAll);
   deepCopy(SurfaceTracerFlux, 0.0_Real);

   Name = Name_;

   NTracers = NTracersIn;
   TimeStep = TimeStepIn;

   defineFields();

} // end constructor

Tendencies::Tendencies(const std::string &Name_, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,     ///< [in] Horizontal mesh
                       VertCoord *VCoord,        ///< [in] Vertical coordinate
                       VertAdv *VAdv,            ///< [in] Vertical advection
                       PressureGrad *PGrad,      ///< [in] Pressure gradient
                       Eos *EqState,             ///< [in] Equation of state
                       VertMix *VMix,            ///< [in] Vertical mixing
                       int NTracersIn,           ///< [in] Number of tracers
                       TimeInterval TimeStepIn,  ///< [in] Time step
                       Config *Options)          ///< [in] Configuration options
    : Tendencies(Name_, Mesh, VCoord, VAdv, PGrad, EqState, VMix, NTracersIn,
                 TimeStepIn, Options, CustomTendencyType{},
                 CustomTendencyType{}) {}

//------------------------------------------------------------------------------
// Compute tendencies for the pseudo-thickness equation
void Tendencies::computePseudoThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocPseudoThicknessTend, PseudoThicknessTend);
   OMEGA_SCOPE(LocThicknessFluxDiv, PseudoThicknessFluxDiv);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);

   Pacer::start("Tend:computePseudoThicknessTendenciesOnly", 1);

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(int K) { LocPseudoThicknessTend(ICell, K) = 0; });
       });

   // Compute pseudo-thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;

   if (LocThicknessFluxDiv.Enabled) {
      Pacer::start("Tend:thicknessFluxDiv", 2);
      parallelForOuter(
          {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocThicknessFluxDiv(LocPseudoThicknessTend, ICell, KChunk,
                                        ThickFluxEdge, NormalVelEdge);
                 });
          });
      Pacer::stop("Tend:thicknessFluxDiv", 2);
   }

   Pacer::start("Tend:computePseudoThicknessVAdvTend", 2);
   // Compute thickness tendency from vertical advection
   VAdv->computePseudoThicknessVAdvTend(PseudoThicknessTend);
   Pacer::stop("Tend:computePseudoThicknessVAdvTend", 2);

   if (CustomThicknessTend) {
      Pacer::start("Tend:customThicknessTend", 2);
      CustomThicknessTend(LocPseudoThicknessTend, State, AuxState,
                          ThickTimeLevel, VelTimeLevel, Time);
      Pacer::stop("Tend:customThicknessTend", 2);
   }

   Pacer::stop("Tend:computePseudoThicknessTendenciesOnly", 1);

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocPotentialVortHAdv, PotentialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);
   OMEGA_SCOPE(LocSfcStressForcing, SfcStressForcing);
   OMEGA_SCOPE(LocBottomDrag, BottomDrag);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocSshCell, VCoord->SshCell);

   Pacer::start("Tend:computeVelocityTendenciesOnly", 1);

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(int K) { LocNormalVelocityTend(IEdge, K) = 0; });
       });

   // Compute potential vorticity horizontal advection
   const Array2DReal &FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   Array2DReal NormVelEdge          = State->getNormalVelocity(VelTimeLevel);
   if (LocPotentialVortHAdv.Enabled) {
      Pacer::start("Tend:PotentialVortHAdv", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocPotentialVortHAdv(LocNormalVelocityTend, IEdge, KChunk,
                                         NormRVortEdge, NormFEdge,
                                         FluxPseudoThickEdge, NormVelEdge);
                 });
          });
      Pacer::stop("Tend:PotentialVortHAdv", 2);
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      Pacer::start("Tend:KEGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocKEGrad(LocNormalVelocityTend, IEdge, KChunk, KECell);
                 });
          });
      Pacer::stop("Tend:KEGrad", 2);
   }

   // Compute sea surface height gradient
   const Array1DReal &SSHCell = LocSshCell;
   if (LocSSHGrad.Enabled) {
      Pacer::start("Tend:SSHGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocSSHGrad(LocNormalVelocityTend, IEdge, KChunk, SSHCell);
                 });
          });
      Pacer::stop("Tend:SSHGrad", 2);
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      Pacer::start("Tend:velocityDiffusion", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocVelocityDiffusion(LocNormalVelocityTend, IEdge, KChunk,
                                         DivCell, RVortVertex);
                 });
          });
      Pacer::stop("Tend:velocityDiffusion", 2);
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      Pacer::start("Tend:velocityHyperDiff", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocVelocityHyperDiff(LocNormalVelocityTend, IEdge, KChunk,
                                         Del2DivCell, Del2RVortVertex);
                 });
          });
      Pacer::stop("Tend:velocityHyperDiff", 2);
   }

   Pacer::start("Tend:computeVelocityVAdvTend", 2);
   // Compute velocity tendency from vertical advection
   VAdv->computeVelocityVAdvTend(NormalVelocityTend, NormVelEdge,
                                 FluxPseudoThickEdge);
   Pacer::stop("Tend:computeVelocityVAdvTend", 2);

   const auto &MeanPseudoThickEdge =
       AuxState->PseudoThicknessAux.MeanPseudoThickEdge;

   // Compute surface stress forcing
   if (LocSfcStressForcing.Enabled) {
      const auto *ForcingState = Forcing::getDefault();
      if (!ForcingState) {
         ABORT_ERROR("Tendencies::computeVelocityTendenciesOnly: surface "
                     "stress forcing is enabled but Forcing has not been "
                     "initialized");
      }
      const auto &NormalStressEdge =
          ForcingState->SfcStressForcing.NormalStressEdge;
      Pacer::start("Tend:sfcStressForcing", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocSfcStressForcing(LocNormalVelocityTend, IEdge, KChunk,
                                        NormalStressEdge, MeanPseudoThickEdge);
                 });
          });
      Pacer::stop("Tend:sfcStressForcing", 2);
   }

   // Compute bottom drag
   if (LocBottomDrag.Enabled) {
      Pacer::start("Tend:bottomDrag", 2);
      parallelFor(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
             LocBottomDrag(LocNormalVelocityTend, IEdge, NormVelEdge, KECell,
                           MeanPseudoThickEdge);
          });
      Pacer::stop("Tend:bottomDrag", 2);
   }

   if (CustomVelocityTend) {
      Pacer::start("Tend:customVelocityTend", 2);
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, ThickTimeLevel,
                         VelTimeLevel, Time);
      Pacer::stop("Tend:customVelocityTend", 2);
   }

   // Compute pressure gradient
   if (PGrad->Enabled) {

      Pacer::start("Tend:pressureGradTerm", 2);
      Array2DReal PseudoThick       = State->getPseudoThickness(ThickTimeLevel);
      const auto &PressureMid       = VCoord->PressureMid;
      const auto &PressureInterface = VCoord->PressureInterface;
      const auto &SpecVol           = EqState->SpecVol;
      const auto &GeomZInterface    = VCoord->GeomZInterface;
      PGrad->computePressureGrad(LocNormalVelocityTend, PressureMid,
                                 PressureInterface, SpecVol, GeomZInterface,
                                 PseudoThick);
      Pacer::stop("Tend:pressureGradTerm", 2);
   }

   // Compute pressure gradient from sloping surfaces (PresGradZ)
   OMEGA_SCOPE(LocPresGradZ, PresGradZ);
   if (LocPresGradZ.Enabled) {
      Pacer::start("Tend:presGradZ", 2);
      const auto &SpecVol           = EqState->SpecVol;
      const auto &LayerThickCell    = State->getPseudoThickness(ThickTimeLevel);
      const auto &PressureInterface = VCoord->PressureInterface;
      const auto &MeanPseudoThickEdge =
          AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocPresGradZ(LocNormalVelocityTend, IEdge, KChunk, SpecVol,
                                 MeanPseudoThickEdge, LayerThickCell,
                                 PressureInterface);
                 });
          });
      Pacer::stop("Tend:presGradZ", 2);
   }

   // Compute pressure gradient force (PresGradForce)
   OMEGA_SCOPE(LocPresGradForce, PresGradForce);
   if (LocPresGradForce.Enabled) {
      Pacer::start("Tend:presGradForce", 2);
      const auto &SpecVol        = EqState->SpecVol;
      const auto &LayerThickCell = State->getPseudoThickness(ThickTimeLevel);
      const auto &PressureMid    = VCoord->PressureMid;
      const auto &MeanPseudoThickEdge =
          AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocPresGradForce(LocNormalVelocityTend, IEdge, KChunk,
                                     SpecVol, MeanPseudoThickEdge,
                                     LayerThickCell, PressureMid);
                 });
          });
      Pacer::stop("Tend:presGradForce", 2);
   }

   // Compute geopotential gradient (GeoptGrad)
   OMEGA_SCOPE(LocGeoptGrad, GeoptGrad);
   if (LocGeoptGrad.Enabled) {
      Pacer::start("Tend:geoptGrad", 2);
      const auto &GeoptMid = VCoord->GeopotentialMid;
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocGeoptGrad(LocNormalVelocityTend, IEdge, KChunk,
                                 GeoptMid);
                 });
          });
      Pacer::stop("Tend:geoptGrad", 2);
   }

   Pacer::stop("Tend:computeVelocityTendenciesOnly", 1);

} // end velocity tendency compute

void Tendencies::computeTracerTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocTracerHorzAdv, TracerHorzAdv);
   OMEGA_SCOPE(LocTracerDiffusion, TracerDiffusion);
   OMEGA_SCOPE(LocTracerHyperDiff, TracerHyperDiff);
   OMEGA_SCOPE(LocSurfaceTracerRestoring, SurfaceTracerRestoring);
   OMEGA_SCOPE(LocSurfaceTracerFlux, SurfaceTracerFlux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computeTracerTendenciesOnly", 1);

   parallelForOuter(
       {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(int K) { LocTracerTend(L, ICell, K) = 0; });
       });

   // compute tracer horizotal advection
   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);
   const Array2DReal &FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   if (LocTracerHorzAdv.Enabled) {
      Pacer::start("Tend:tracerHorzAdv", 2);
      parallelForOuter(
          {NTracers, Mesh->NEdgesAll},
          KOKKOS_LAMBDA(int L, int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerHorzAdv(L, IEdge, KChunk, TracerArray,
                                     FluxPseudoThickEdge, NormalVelEdge);
                 });
          });
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerHorzAdv(LocTracerTend, L, ICell, KChunk);
                 });
          });
      Pacer::stop("Tend:tracerHorzAdv", 2);
   }

   // compute tracer diffusion
   const Array2DReal &MeanPseudoThickEdge =
       AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
   if (LocTracerDiffusion.Enabled) {
      Pacer::start("Tend:tracerDiffusion", 2);
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerDiffusion(LocTracerTend, L, ICell, KChunk,
                                       TracerArray, MeanPseudoThickEdge);
                 });
          });
      Pacer::stop("Tend:tracerDiffusion", 2);
   }

   // compute tracer hyperdiffusion
   const Array3DReal &Del2TracersCell = AuxState->TracerAux.Del2TracersCell;
   if (LocTracerHyperDiff.Enabled) {
      Pacer::start("Tend:tracerHyperDiff", 2);
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerHyperDiff(LocTracerTend, L, ICell, KChunk,
                                       Del2TracersCell);
                 });
          });
      Pacer::stop("Tend:tracerHyperDiff", 2);
   }

   // compute tracer tendency from vertical advection
   Pacer::start("Tend:computeTracerVAdvTend", 2);
   Array2DReal ThicknessForVAdv;
   if (VAdv->VertAdvChoice == VertAdvOption::Standard) {
      ThicknessForVAdv = State->getPseudoThickness(ThickTimeLevel);
   } else if (VAdv->VertAdvChoice == VertAdvOption::FCT) {
      ThicknessForVAdv = AuxState->PseudoThicknessAux.ProvPseudoThickness;
   }
   VAdv->computeTracerVAdvTend(LocTracerTend, TracerArray, ThicknessForVAdv,
                               TimeStep);
   Pacer::stop("Tend:computeTracerVAdvTend", 2);

   // compute tracer surface restoring
   const Array2DReal &TracersMonthlySurfClimo =
       AuxState->SurfTracerRestAux.TracersMonthlySurfClimoCell;
   const I4 NTracersToRestore = LocSurfaceTracerRestoring.NTracersToRestore;
   const auto &TracerIdsToRestore =
       LocSurfaceTracerRestoring.TracerIdsToRestore;
   if (LocSurfaceTracerRestoring.Enabled && NTracersToRestore > 0) {
      Pacer::start("Tend:surfaceTracerRestoring", 2);
      parallelFor(
          {NTracersToRestore, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int R, int ICell) {
             const int KMin = MinLayerCell(ICell);
             const int L    = TracerIdsToRestore(R);
             LocSurfaceTracerRestoring(LocTracerTend, L, ICell, KMin,
                                       TracersMonthlySurfClimo, TracerArray);
          });
      Pacer::stop("Tend:surfaceTracerRestoring", 2);
   }

   // Compute KPP non-local tracer tendency
   if (TracerNonLocalFluxEnabled) {
      KPPMix *KPPInstance = KPPMix::getInstance();
      if (KPPInstance && KPPInstance->Enabled) {
         Pacer::start("Tend:tracerNonLocalFlux", 2);
         OMEGA_SCOPE(LocNonLocalFlux, KPPInstance->VertNonLocalFlux);
         parallelForOuter(
             {NTracers, Mesh->NCellsAll},
             KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
                const int KMin   = MinLayerCell(ICell);
                const int KMax   = MaxLayerCell(ICell);
                const int KRange = vertRangeChunked(KMin, KMax);
                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KChunk) {
                       const I4 KStart = chunkStart(KChunk, KMin);
                       const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                       for (int KVec = 0; KVec < KLen; ++KVec) {
                          const I4 K = KStart + KVec;
                          LocTracerTend(L, ICell, K) +=
                              LocSurfaceTracerFlux(L, ICell) *
                                  LocNonLocalFlux(ICell, K) -
                              LocSurfaceTracerFlux(L, ICell) *
                                  LocNonLocalFlux(ICell, K + 1);
                       }
                    });
             });
         Pacer::stop("Tend:tracerNonLocalFlux", 2);
      }
   }

   if (UseTempSurfaceTracerFluxBridge && UseTempTopLayerFluxTendencyBridge) {
      I4 TempIdx = -1;
      if (Tracers::getIndex(TempIdx, "Temperature") == 0) {
         const I4 TempTracerIndex = TempIdx;
         parallelFor(
             "TempTopLayerFluxTendencyBridge", {Mesh->NCellsAll},
             KOKKOS_LAMBDA(I4 ICell) {
                const I4 KMin = MinLayerCell(ICell);
                LocTracerTend(TempTracerIndex, ICell, KMin) +=
                    LocSurfaceTracerFlux(TempTracerIndex, ICell);
             });
      }
   }

   Pacer::stop("Tend:computeTracerTendenciesOnly", 1);
} // end tracer tendency compute

void Tendencies::computePseudoThicknessTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   // only need PseudoThicknessAux on edge
   Array2DReal PseudoThick = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormVel     = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(PseudoThicknessAux, AuxState->PseudoThicknessAux);
   OMEGA_SCOPE(PseudoThickCell, PseudoThick);
   OMEGA_SCOPE(NormalVelEdge, NormVel);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computePseudoThicknessTendencies", 1);

   Pacer::start("Tend:computePseudoThickAux", 2);
   parallelForOuter(
       "computePseudoThickAux", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 PseudoThicknessAux.computeVarsOnEdge(
                     IEdge, KChunk, PseudoThickCell, NormalVelEdge);
              });
       });
   Pacer::stop("Tend:computePseudoThickAux", 2);

   computePseudoThicknessTendenciesOnly(State, AuxState, ThickTimeLevel,
                                        VelTimeLevel, Time);

   Pacer::stop("Tend:computePseudoThicknessTendencies", 1);
}

void Tendencies::computeVelocityTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt ///< [in] Time interval for projection over the
                        ///< current time stepper stage
) {
   Pacer::start("Tend:computeVelocityTendencies", 1);

   if (StageVerticalMixingEnabled) {
      computeStageVerticalMixing(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel);
   }
   AuxState->computeMomAux(State, TracerArray, ThickTimeLevel, VelTimeLevel,
                           ProjDt);
   computeVelocityTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel, TracerTimeLevel, Time);

   Pacer::stop("Tend:computeVelocityTendencies", 1);
}

void Tendencies::computeTracerTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(TracerAux, AuxState->TracerAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computeTracerTendencies", 1);

   if (StageVerticalMixingEnabled) {
      computeStageVerticalMixing(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel);
   }
   const auto &MeanPseudoThickEdge =
       AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
   Pacer::start("Tend:computeTracerAuxCell", 2);
   parallelForOuter(
       "computeTracerAuxCell", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 TracerAux.computeVarsOnCells(LTracer, ICell, KChunk,
                                              MeanPseudoThickEdge, TracerArray);
              });
       });
   Pacer::stop("Tend:computeTracerAuxCell", 2);

   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);

   Pacer::stop("Tend:computeTracerTendencies", 1);
}

//------------------------------------------------------------------------------
// Compute both pseudo-thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt ///< [in] Time interval for projection over the
                        ///< current time stepper stage
) {
   AuxState->computeAll(State, TracerArray, ThickTimeLevel, VelTimeLevel,
                        ProjDt);

   if (StageVerticalMixingEnabled) {
      computeStageVerticalMixing(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel);
   }
   computePseudoThicknessTendenciesOnly(State, AuxState, ThickTimeLevel,
                                        VelTimeLevel, Time);
   computeVelocityTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel, TracerTimeLevel, Time);
   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);
} // end all tendency compute

//------------------------------------------------------------------------------
// Set surface tracer flux for use by KPP non-local tracer tendency
void Tendencies::setSurfaceTracerFlux(const Array2DReal &Flux) {
   OMEGA_REQUIRE(Flux.extent(0) == SurfaceTracerFlux.extent(0),
                 "Tendencies::setSurfaceTracerFlux: tracer dimension mismatch");
   OMEGA_REQUIRE(Flux.extent(1) == SurfaceTracerFlux.extent(1),
                 "Tendencies::setSurfaceTracerFlux: cell dimension mismatch");
   Kokkos::deep_copy(SurfaceTracerFlux, Flux);
}

//------------------------------------------------------------------------------
// Run KPP mixing for the current stage. VertMixImplicit merges the KPP
// coefficients into VertMix after the base coefficients are recomputed.
void Tendencies::computeStageVerticalMixing(const OceanState *State,
                                            const AuxiliaryState *AuxState,
                                            const Array3DReal &TracerArray,
                                            int ThickTimeLevel,
                                            int VelTimeLevel) {
   (void)AuxState;
   KPPMix *KPPInstance = KPPMix::getInstance();

   if (!EqState || !KPPInstance || !KPPInstance->Enabled)
      return;

   I4 TempIdx = -1;
   I4 SaltIdx = -1;
   if (Tracers::getIndex(TempIdx, "Temperature") != 0 ||
       Tracers::getIndex(SaltIdx, "Salinity") != 0) {
      LOG_WARN("Tendencies::computeStageVerticalMixing: Temperature/Salinity "
               "tracers not found, skipping KPP stage update");
      return;
   }

   const I4 NCellsAll   = Mesh->NCellsAll;
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal ConservTemp("KPP-ConservTemp", NCellsAll, NVertLayers);
   Array2DReal AbsSalinity("KPP-AbsSalinity", NCellsAll, NVertLayers);
   parallelFor(
       "KPP-ExtractTS", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          ConservTemp(ICell, K) = TracerArray(TempIdx, ICell, K);
          AbsSalinity(ICell, K) = TracerArray(SaltIdx, ICell, K);
       });

   Array2DReal LayerThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge  = State->getNormalVelocity(VelTimeLevel);

   Array1DReal SurfacePressure("KPP-SurfacePressure", NCellsAll);
   deepCopy(SurfacePressure, 1.0e5_Real);
   const_cast<VertCoord *>(VCoord)->computePressure(LayerThickCell,
                                                    SurfacePressure);

   OMEGA_SCOPE(PressureMid, VCoord->PressureMid);
   Array2DReal PressureMidDbar("KPP-PressureMidDbar", NCellsAll, NVertLayers);
   parallelFor(
       "KPP-PressureToDbar", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          PressureMidDbar(ICell, K) = PressureMid(ICell, K) * 1.0e-4_Real;
       });

   EqState->computeSpecVol(ConservTemp, AbsSalinity, PressureMidDbar);
   EqState->computeBruntVaisalaFreqSq(ConservTemp, AbsSalinity, PressureMidDbar,
                                      EqState->SpecVol);

   Array2DReal PotentialDensity("KPP-PotentialDensity", NCellsAll, NVertLayers);
   OMEGA_SCOPE(SpecVol, EqState->SpecVol);
   parallelFor(
       "KPP-PotentialDensity", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          PotentialDensity(ICell, K) =
              1.0_Real / Kokkos::max(1.0e-12_Real, SpecVol(ICell, K));
       });

   // Compute tangential velocity on edges (same pattern as
   // VertMix::VertMixImplicit)
   Array2DReal TangentialVelEdge("KPP-TangentialVelEdge", Mesh->NEdgesSize,
                                 NVertLayers);
   {
      TangentialReconOnEdge TanReconEdge(Mesh);
      OMEGA_SCOPE(LocTangentialVelEdge, TangentialVelEdge);
      OMEGA_SCOPE(MinLayerEdgeTop, VCoord->MinLayerEdgeTop);
      OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeTop(IEdge);
             const int KMax   = MaxLayerEdgeBot(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    TanReconEdge(LocTangentialVelEdge, IEdge, KChunk,
                                 NormalVelEdge);
                 });
          });
   }

   Array1DReal IceFraction("KPP-IceFraction", NCellsAll);

   // Use KPPInstance's member arrays so these values are retained for output
   OMEGA_SCOPE(LocSurfaceFrictionVelocity,
               KPPInstance->SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfaceBuoyancyFlux, KPPInstance->SurfaceBuoyancyFlux);

   const Real LocKPPHeatFluxToBuoyancyFactor = KPPHeatFluxToBuoyancyFactor;
   const Real LocKPPThicknessFluxToBuoyancyFactor =
       KPPThicknessFluxToBuoyancyFactor;

   const auto *ForcingState = Forcing::getDefault();
   if (!ForcingState) {
      LOG_WARN("Tendencies::computeStageVerticalMixing: Forcing has not "
               "been initialized, skipping KPP stage update");
      return;
   }

   const auto &SfcStressForcing = ForcingState->SfcStressForcing;
   OMEGA_SCOPE(ZonalStressCell, SfcStressForcing.ZonalStressCell);
   OMEGA_SCOPE(MeridStressCell, SfcStressForcing.MeridStressCell);
   OMEGA_SCOPE(LocLatentHeatFlux, SfcStressForcing.LatentHeatFlux);
   OMEGA_SCOPE(LocSensibleHeatFlux, SfcStressForcing.SensibleHeatFlux);
   OMEGA_SCOPE(LocShortWaveHeatFlux, SfcStressForcing.ShortWaveHeatFlux);
   OMEGA_SCOPE(LocEvaporationFlux, SfcStressForcing.EvaporationFlux);
   OMEGA_SCOPE(LocRainFlux, SfcStressForcing.RainFlux);
   OMEGA_SCOPE(LocRiverRunoffFlux, SfcStressForcing.RiverRunoffFlux);
   OMEGA_SCOPE(LocIceRunoffFlux, SfcStressForcing.IceRunoffFlux);
   OMEGA_SCOPE(LocSubglacialRunoffFlux, SfcStressForcing.SubglacialRunoffFlux);
   OMEGA_SCOPE(LocIcebergFreshWaterFlux,
               SfcStressForcing.IcebergFreshWaterFlux);
   OMEGA_SCOPE(LocSurfaceTracerFlux, SurfaceTracerFlux);

   const bool LocUseTempSurfaceTracerFluxBridge =
       TracerNonLocalFluxEnabled && UseTempSurfaceTracerFluxBridge;
   const I4 TempTracerIndex              = TempIdx;
   const Real HeatFluxToTracerFluxFactor = 1._Real / (RhoSw * CpSw);

   if (LocUseTempSurfaceTracerFluxBridge) {
      deepCopy(SurfaceTracerFlux, 0.0_Real);
   }

   parallelFor(
       "KPP-SurfaceForcing", {NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const Real tau_x   = ZonalStressCell(ICell);
          const Real tau_y   = MeridStressCell(ICell);
          const Real tau_mag = Kokkos::sqrt(tau_x * tau_x + tau_y * tau_y);
          LocSurfaceFrictionVelocity(ICell) =
              Kokkos::sqrt(Kokkos::max(0.0_Real, tau_mag / RhoSw));
          const Real heat_flux = LocLatentHeatFlux(ICell) +
                                 LocSensibleHeatFlux(ICell) +
                                 LocShortWaveHeatFlux(ICell);
          const Real freshwater_flux =
              LocRainFlux(ICell) + LocRiverRunoffFlux(ICell) +
              LocIceRunoffFlux(ICell) + LocSubglacialRunoffFlux(ICell) +
              LocIcebergFreshWaterFlux(ICell) - LocEvaporationFlux(ICell);
          LocSurfaceBuoyancyFlux(ICell) =
              heat_flux * LocKPPHeatFluxToBuoyancyFactor +
              freshwater_flux * LocKPPThicknessFluxToBuoyancyFactor;
          if (LocUseTempSurfaceTracerFluxBridge) {
             LocSurfaceTracerFlux(TempTracerIndex, ICell) =
                 heat_flux * HeatFluxToTracerFluxFactor;
          }
          IceFraction(ICell) = 0.0_Real;
       });

   Array1DReal WindSpeed10m;
   KPPInstance->computeKPPMix(
       PotentialDensity, NormalVelEdge, TangentialVelEdge,
       KPPInstance->SurfaceFrictionVelocity, KPPInstance->SurfaceBuoyancyFlux,
       EqState->BruntVaisalaFreqSq, IceFraction, WindSpeed10m);
   // KPP VertDiff/VertVisc are merged into VertMix inside VertMixImplicit(),
   // after computeVertMix() resets them, so no merge needed here.
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//

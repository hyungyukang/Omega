//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "share/ManufacturedSolution.h"

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh as alread been initialized.
int Tendencies::init() {

   int Err = 0;

   HorzMesh *DefHorzMesh = HorzMesh::getDefault();

   // Get TendConfig group if available
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   if (OmegaConfig->existsGroup("Tendencies")) {
      Err = OmegaConfig->get(TendConfig);
   }

   int NVertLevels = 60;

   // Retrieve NVertLevels from Config if available
   Config DimConfig("Dimension");
   if (OmegaConfig->existsGroup("Dimension")) {
      Err = OmegaConfig->get(DimConfig);
      if (DimConfig.existsVar("NVertLevels")) {
         Err = DimConfig.get("NVertLevels", NVertLevels);
      }
   }

   // TODO: move setMasks to HorzMesh constructor
   DefHorzMesh->setMasks(NVertLevels);

   Tendencies::DefaultTendencies =
       create("Default", DefHorzMesh, NVertLevels, &TendConfig);

   // Check if use the manufactured solution test
   bool UseManufacturedSolution = false;
   Config ManufacturedSolutionConfig("ManufacturedSolution");
   if (OmegaConfig->existsGroup("ManufacturedSolution")) {
      Err = OmegaConfig->get(ManufacturedSolutionConfig);
      if (ManufacturedSolutionConfig.existsVar("UseManufacturedSolution")) {
         Err = ManufacturedSolutionConfig.get("UseManufacturedSolution", UseManufacturedSolution);
      }
   }

   if (UseManufacturedSolution) {
      // Clear tendencies
      Tendencies::clear();

      //Tendencies::CustomTendencyType{} = ManufacturedThicknessTendency{};
      //Tendencies::CustomVelocityType{} = ManufacturedVelocityTendency{};

      // Re-create tendencies with the manufactured solution
      Tendencies::DefaultTendencies =
          create("Default", DefHorzMesh, NVertLevels, &TendConfig,
                 ManufacturedThicknessTendency{}, ManufacturedVelocityTendency{});
                // Tendencies::CustomTendencyType{},Tendencies::CustomVelocityType{});
                // CustomTendencyType{}, ManufacturedVelocityTendency{});

      //auto *TestTendencies = Tendencies::create(
      //   "TestTendencies", DefHorzMesh, NVertLevels, &TendConfig,
      //   ManufacturedThicknessTendency{}, ManufacturedVelocityTendency{});
   }

   return Err;

} // end init

//------------------------------------------------------------------------------
// Destroys the tendencies
Tendencies::~Tendencies() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Removes all tendencies instances before exit
void Tendencies::clear() { AllTendencies.clear(); } // end clear

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
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       Config *Options, ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : ThicknessFluxDiv(Mesh, Options), PotientialVortHAdv(Mesh, Options),
      KEGrad(Mesh, Options), SSHGrad(Mesh, Options),
      VelocityDiffusion(Mesh, Options), VelocityHyperDiff(Mesh, Options),
      CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend) {

   // Tendency arrays
   LayerThicknessTend =
       Array2DReal("LayerThicknessTend", Mesh->NCellsSize, NVertLevels);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, NVertLevels);

   // Array dimension lengths
   NCellsAll = Mesh->NCellsAll;
   NEdgesAll = Mesh->NEdgesAll;
   NChunks   = NVertLevels / VecLength;

} // end constructor

Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       Config *Options) ///< [in] Configuration options
    : Tendencies(Name, Mesh, NVertLevels, Options, CustomTendencyType{},
                 CustomTendencyType{}) {}

//------------------------------------------------------------------------------
// Compute tendencies for layer thickness equation
void Tendencies::computeThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocLayerThicknessTend, LayerThicknessTend);
   OMEGA_SCOPE(LocThicknessFluxDiv, ThicknessFluxDiv);
   const Array2DReal &NormalVelEdge = State->NormalVelocity[VelTimeLevel];

   deepCopy(LocLayerThicknessTend, 0);

   // Compute thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;

   if (LocThicknessFluxDiv.Enabled) {
      parallelFor(
          {NCellsAll, NChunks}, KOKKOS_LAMBDA(int ICell, int KChunk) {
             LocThicknessFluxDiv(LocLayerThicknessTend, ICell, KChunk,
                                 ThickFluxEdge, NormalVelEdge);
          });
   }

   if (CustomThicknessTend) {
      CustomThicknessTend(LocLayerThicknessTend, State, AuxState,
                          ThickTimeLevel, VelTimeLevel, Time);
   }

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocPotientialVortHAdv, PotientialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);

   deepCopy(LocNormalVelocityTend, 0);

   // Compute potential vorticity horizontal advection
   const Array2DReal &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   const Array2DReal &NormVelEdge   = State->NormalVelocity[VelTimeLevel];
   if (LocPotientialVortHAdv.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocPotientialVortHAdv(LocNormalVelocityTend, IEdge, KChunk,
                                   NormRVortEdge, NormFEdge, FluxLayerThickEdge,
                                   NormVelEdge);
          });
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocKEGrad(LocNormalVelocityTend, IEdge, KChunk, KECell);
          });
   }

   // Compute sea surface height gradient
   const Array2DReal &SSHCell = AuxState->LayerThicknessAux.SshCell;
   if (LocSSHGrad.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocSSHGrad(LocNormalVelocityTend, IEdge, KChunk, SSHCell);
          });
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityDiffusion(LocNormalVelocityTend, IEdge, KChunk, DivCell,
                                  RVortVertex);
          });
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      parallelFor(
          {NEdgesAll, NChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityHyperDiff(LocNormalVelocityTend, IEdge, KChunk,
                                  Del2DivCell, Del2RVortVertex);
          });
   }

   if (CustomVelocityTend) {
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, ThickTimeLevel,
                         VelTimeLevel, Time);
   }

} // end velocity tendency compute

void Tendencies::computeThicknessTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   // only need LayerThicknessAux on edge
   OMEGA_SCOPE(LayerThicknessAux, AuxState->LayerThicknessAux);
   OMEGA_SCOPE(LayerThickCell, State->LayerThickness[ThickTimeLevel]);
   OMEGA_SCOPE(NormalVelEdge, State->NormalVelocity[VelTimeLevel]);

   parallelFor(
       "computeLayerThickAux", {NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });

   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);
}

void Tendencies::computeVelocityTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   AuxState->computeAll(State, ThickTimeLevel, VelTimeLevel);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);
}

//------------------------------------------------------------------------------
// Compute both layer thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   AuxState->computeAll(State, ThickTimeLevel, VelTimeLevel);
   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);

} // end all tendency compute

// TODO: Implement Config options for all constructors
ThicknessFluxDivOnCell::ThicknessFluxDivOnCell(const HorzMesh *Mesh,
                                               Config *TendConfig)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) {

   if (TendConfig->existsVar("ThicknessFluxTendencyEnable")) {
      TendConfig->get("ThicknessFluxTendencyEnable", Enabled);
   }
}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 Config *TendConfig)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge) {

   if (TendConfig->existsVar("PVTendencyEnable")) {
      TendConfig->get("PVTendencyEnable", Enabled);
   }
}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, Config *TendConfig)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge) {

   if (TendConfig->existsVar("KETendencyEnable")) {
      TendConfig->get("KETendencyEnable", Enabled);
   }
}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, Config *TendConfig)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge) {

   if (TendConfig->existsVar("SSHTendencyEnable")) {
      TendConfig->get("SSHTendencyEnable", Enabled);
   }
}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 Config *TendConfig)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(Mesh->EdgeMask) {

   if (TendConfig->existsVar("VelDiffTendencyEnable")) {
      TendConfig->get("VelDiffTendencyEnable", Enabled);
   }
   if (TendConfig->existsVar("ViscDel2")) {
      TendConfig->get("ViscDel2", ViscDel2);
   }
}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 Config *TendConfig)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(Mesh->EdgeMask) {

   if (TendConfig->existsVar("VelHyperDiffTendencyEnable")) {
      TendConfig->get("VelHyperDiffTendencyEnable", Enabled);
   }
   if (TendConfig->existsVar("ViscDel4")) {
      TendConfig->get("ViscDel4", ViscDel4);
   }
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//

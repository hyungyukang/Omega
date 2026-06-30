//===-- ocn/KPPMix.cpp - KPP Boundary Layer Mixing Implementation --*- C++
//-*-===//
//
/// \file
/// \brief Implementation of KPP boundary layer mixing scheme
///
/// This file implements the KPPMix class for computing ocean boundary layer
/// mixing following Large et al. (1994) with optional Langmuir enhancement.
//
//===----------------------------------------------------------------------===//

#include "KPPMix.h"
#include "DataTypes.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "KPPConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <limits>

namespace OMEGA {

// Singleton instance
KPPMix *KPPMix::Instance = nullptr;

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &Name_in, const HorzMesh *Mesh_in,
               const VertCoord *VCoord_in)
    : Name(Name_in), Mesh(Mesh_in), VCoord(VCoord_in) {

   // Allocate output arrays
   VertDiff = Array2DReal("VertDiff", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   VertVisc = Array2DReal("VertVisc", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   BoundaryLayerDepth = Array1DReal("BoundaryLayerDepth", Mesh->NCellsAll);
   IndexBoundaryLayerDepth =
       Array1DI4("IndexBoundaryLayerDepth", Mesh->NCellsAll);
   VertNonLocalFlux       = Array2DReal("VertNonLocalFlux", Mesh->NCellsAll,
                                        VCoord->NVertLayers + 1);
   BulkRichardsonNumber   = Array2DReal("BulkRichardsonNumber", Mesh->NCellsAll,
                                        VCoord->NVertLayers + 1);
   TurbulentVelocityScale = Array2DReal(
       "TurbulentVelocityScale", Mesh->NCellsAll, VCoord->NVertLayers + 1);

   // Set field names
   VertDiffFldName          = "VertDiff";
   VertViscFldName          = "VertVisc";
   OBLDepthFldName          = "BoundaryLayerDepth";
   NonLocalFluxFldName      = "VertNonLocalFlux";
   BulkRichardsonFldName    = "BulkRichardsonNumber";
   TurbulentVelScaleFldName = "TurbulentVelocityScale";

   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      OBLDepthFldName.append(Name);
      NonLocalFluxFldName.append(Name);
      BulkRichardsonFldName.append(Name);
      TurbulentVelScaleFldName.append(Name);
   }

   defineFields();
}

/// Destructor for KPPMix
KPPMix::~KPPMix() {}

/// Get singleton instance
KPPMix *KPPMix::getInstance() { return Instance; }

/// Destroy singleton instance
void KPPMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initialize KPPMix from configuration
void KPPMix::init() {
   if (!Instance) {
      Instance = new KPPMix("Default", HorzMesh::getDefault(),
                            VertCoord::getDefault());
   }

   Error Err;
   KPPMix *DefKPPMix   = KPPMix::getInstance();
   Config *OmegaConfig = Config::getOmegaConfig();

   // Get VertMix config group
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "KPPMix::init: VertMix group not found in Config");

   // Get KPP config subgroup
   Config KPPConfig("KPP");
   Err += VertMixConfig.get(KPPConfig);
   if (Err.isFail()) {
      LOG_INFO("KPPMix::init: KPP subgroup not found, using defaults");
      return; // Continue with defaults
   }

   // Read KPP parameters
   bool enable = true;
   Err += KPPConfig.get("Enable", enable);
   DefKPPMix->Enabled = enable;

   Err += KPPConfig.get("CriticalBulkRichardsonNumber",
                        DefKPPMix->CriticalRichardson);
   Err += KPPConfig.get("StopOBLSearch", DefKPPMix->StopOBLSearchMult);
   Err += KPPConfig.get("SurfaceLayerExtent", DefKPPMix->SurfaceLayerExtent);

   // KPP matching/profile semantics. Keep legacy key for backward
   // compatibility.
   Error MatchErr =
       KPPConfig.get("MatchTechnique", DefKPPMix->MatchTechniqueStr);
   if (!MatchErr.isSuccess()) {
      MatchErr.reset();
      std::string legacy_scheme = DefKPPMix->MatchTechniqueStr;
      Error LegacyMatchErr =
          KPPConfig.get("BoundaryLayerDepthScheme", legacy_scheme);
      if (LegacyMatchErr.isSuccess()) {
         DefKPPMix->MatchTechniqueStr = legacy_scheme;
      } else {
         LegacyMatchErr.reset();
      }
   }
   Error InterpErr = KPPConfig.get("InterpType2", DefKPPMix->InterpType2Str);
   if (!InterpErr.isSuccess()) {
      InterpErr.reset();
   }
   Error EnhancedErr =
       KPPConfig.get("UseEnhancedDiffusion", DefKPPMix->UseEnhancedDiffusion);
   if (!EnhancedErr.isSuccess()) {
      EnhancedErr.reset();
   }

   // Keep active options focused on what is used in OMEGA.
   if (DefKPPMix->MatchTechniqueStr == "MatchGradient") {
      LOG_INFO("KPPMix::init: MatchGradient is deprecated/unused in OMEGA; "
               "mapping to SimpleShapes");
      DefKPPMix->MatchTechniqueStr = "SimpleShapes";
   }
   if (DefKPPMix->MatchTechniqueStr != "SimpleShapes" &&
       DefKPPMix->MatchTechniqueStr != "MatchBoth" &&
       DefKPPMix->MatchTechniqueStr != "ParabolicNonLocal") {
      LOG_INFO(
          "KPPMix::init: Unsupported MatchTechnique='{}', using SimpleShapes",
          DefKPPMix->MatchTechniqueStr);
      DefKPPMix->MatchTechniqueStr = "SimpleShapes";
   }

   // Wave and flux options
   Err += KPPConfig.get("UseLangmuirCirculation",
                        DefKPPMix->UseLangmuirCirculation);
   Err += KPPConfig.get("UseNonLocalFlux", DefKPPMix->UseNonLocalFlux);
   Err += KPPConfig.get("IceFractionThresholdForLangmuir",
                        DefKPPMix->IceFractionThresholdForLangmuir);
   Err += KPPConfig.get("IceFractionThresholdForMinimumOBL",
                        DefKPPMix->IceFractionThresholdForMinimumOBL);
   Err +=
       KPPConfig.get("MinimumOBLUnderSeaIce", DefKPPMix->MinimumOBLUnderSeaIce);
   Error DebugErr =
       KPPConfig.get("DebugDiagnostics", DefKPPMix->DebugDiagnostics);
   if (!DebugErr.isSuccess()) {
      DebugErr.reset();
      DefKPPMix->DebugDiagnostics = false;
   }

   // Background mixing
   Err += KPPConfig.get("BackgroundViscosity", DefKPPMix->BackgroundVisc);
   Err += KPPConfig.get("BackgroundDiffusivity", DefKPPMix->BackgroundDiff);

   LOG_INFO("KPPMix::init: KPP initialized");
}

/// Main computation routine
void KPPMix::computeKPPMix(const Array2DReal &PotentialDensity,
                           const Array2DReal &NormalVelocity,
                           const Array2DReal &TangentialVelocity,
                           const Array1DReal &SurfaceFrictionVelocity,
                           const Array1DReal &SurfaceBuoyancyFlux,
                           const Array2DReal &BruntVaisalaFreqSq,
                           const Array1DReal &IceFraction,
                           const Array1DReal &WindSpeed10m) {

   if (!Enabled) {
      return; // Skip if disabled
   }

   // =======================================================================
   // Stage 1: Compute OBL Depth
   // =======================================================================
   computeOBLDepth(PotentialDensity, NormalVelocity, TangentialVelocity,
                   SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                   BruntVaisalaFreqSq, IceFraction, WindSpeed10m);

   // =======================================================================
   // Stage 2: Compute Mixing Coefficients
   // =======================================================================
   computeMixingCoefficients(PotentialDensity, NormalVelocity,
                             TangentialVelocity, SurfaceFrictionVelocity,
                             SurfaceBuoyancyFlux);

   if (DebugDiagnostics) {
      logDiagnostics(PotentialDensity, NormalVelocity, TangentialVelocity,
                     SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                     WindSpeed10m);
   }
}

void KPPMix::logDiagnostics(const Array2DReal &PotentialDensity,
                            const Array2DReal &NormalVelocity,
                            const Array2DReal &TangentialVelocity,
                            const Array1DReal &SurfaceFrictionVelocity,
                            const Array1DReal &SurfaceBuoyancyFlux,
                            const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   const auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   const auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);
   const auto ZInterfaceH   = createHostMirrorCopy(VCoord->GeomZInterface);
   const auto DensityH      = createHostMirrorCopy(PotentialDensity);
   const auto UVelH         = createHostMirrorCopy(NormalVelocity);
   const auto VVelH         = createHostMirrorCopy(TangentialVelocity);
   const auto UStarH        = createHostMirrorCopy(SurfaceFrictionVelocity);
   const auto B0H           = createHostMirrorCopy(SurfaceBuoyancyFlux);
   const auto OBLDepthH     = createHostMirrorCopy(BoundaryLayerDepth);
   const auto OBLIndexH     = createHostMirrorCopy(IndexBoundaryLayerDepth);

   const int NCellsAll = Mesh->NCellsAll;
   if (NCellsAll <= 0) {
      return;
   }

   const int ICell       = 0;
   const int KMin        = MinLayerCellH(ICell);
   const int KMax        = MaxLayerCellH(ICell);
   const int NVertLayers = VCoord->NVertLayers;

   if (KMin > KMax) {
      return;
   }

   const Real rho_ref    = 1025.0_Real;
   const int KSurf       = Kokkos::min(KMin, NVertLayers - 1);
   const Real rho_surf   = DensityH(ICell, KSurf);
   const Real u_surf     = UVelH(ICell, KSurf);
   const Real v_surf     = VVelH(ICell, KSurf);
   const Real u_star     = UStarH(ICell);
   const Real u_star_eff = Kokkos::fmax(KPP::MIN_USTAR, u_star);
   const Real b0         = B0H(ICell);
   Real u10              = 0.0_Real;
   if (WindSpeed10m.extent(0) > 0) {
      const auto Wind10mH = createHostMirrorCopy(WindSpeed10m);
      u10                 = Wind10mH(ICell);
   }
   const Real langmuir_factor =
       UseLangmuirCirculation ? ComputeEnhancementFactor(u10, u_star_eff, 50.0)
                              : 1.0_Real;
   const Real b0_eff = b0 * langmuir_factor;

   LOG_INFO("KPP debug: cell={} h_obl={} m k_obl={} u*={} b0={} b0_eff={} "
            "langmuir={}",
            ICell, OBLDepthH(ICell), OBLIndexH(ICell), u_star, b0, b0_eff,
            langmuir_factor);

   const int KTop   = Kokkos::min(KMax, KMin + 3);
   const int k_obl  = OBLIndexH(ICell);
   const Real h_obl = OBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int kCell    = Kokkos::min(K, NVertLayers - 1);
      const int kInt     = Kokkos::min(K + 1, NVertLayers);
      const Real z_depth = Kokkos::abs(ZInterfaceH(ICell, kInt));

      const Real rho_k = DensityH(ICell, kCell);
      const Real u_k   = UVelH(ICell, kCell);
      const Real v_k   = VVelH(ICell, kCell);

      const Real delta_rho = rho_k - rho_surf;
      const Real delta_b   = Gravity * delta_rho / rho_ref;
      const Real du        = u_k - u_surf;
      const Real dv        = v_k - v_surf;
      const Real shear2    = du * du + dv * dv;
      const Real w_turb =
          ComputeTurbulentVelocityScale(u_star_eff, b0_eff, z_depth);
      const Real vel_scale2 = shear2 + w_turb * w_turb + 1.0e-12_Real;
      const Real ri_b       = delta_b * z_depth / vel_scale2;

      Real sigma = 0.0_Real;
      if (K <= k_obl) {
         sigma = -1.0_Real * static_cast<Real>(K - KMin) /
                 static_cast<Real>(k_obl - KMin + 1);
         sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma));
      }

      const Real z_local = -sigma * h_obl;
      Real zeta          = 0.0_Real;
      const Real denom   = VonKar * b0;
      if (Kokkos::abs(denom) > 1.0e-16_Real) {
         const Real l_mo = (u_star_eff * u_star_eff * u_star_eff) / denom;
         if (Kokkos::abs(l_mo) > 1.0e-16_Real) {
            zeta = z_local / l_mo;
         }
      }

      const Real phi_m = KPP::KPPProfileM2(zeta);
      const Real phi_s = KPP::KPPProfileS2(zeta);

      LOG_INFO(
          "KPP debug top: cell={} k={} z={} ri_b={} zeta={} phi_m={} phi_s={}",
          ICell, K, z_depth, ri_b, zeta, phi_m, phi_s);
   }
}

/// Stage 1: Compute OBL depth using bulk Richardson search
void KPPMix::computeOBLDepth(const Array2DReal &PotentialDensity,
                             const Array2DReal &NormalVelocity,
                             const Array2DReal &TangentialVelocity,
                             const Array1DReal &SurfaceFrictionVelocity,
                             const Array1DReal &SurfaceBuoyancyFlux,
                             const Array2DReal &BruntVaisalaFreqSq,
                             const Array1DReal &IceFraction,
                             const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsAll);
   const Real LocIceFracThresholdForLangmuir = IceFractionThresholdForLangmuir;
   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          const Real iceFrac = IceFraction(ICell);
          if (UseLangmuirCirculation &&
              iceFrac < LocIceFracThresholdForLangmuir) {
             const Real uStar = SurfaceFrictionVelocity(ICell);
             const Real u10 =
                 (WindSpeed10m.extent(0) > 0) ? WindSpeed10m(ICell) : 0.0_Real;
             LangmuirFactor(ICell) = ComputeEnhancementFactor(u10, uStar, 50.0);
          } else {
             LangmuirFactor(ICell) = 1.0;
          }
       });

   // =======================================================================
   // Stage 1: Compute OBL depth - simplified implementation
   // =======================================================================

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(ZInterface, VCoord->GeomZInterface);
   OMEGA_SCOPE(ZMid, VCoord->GeomZMid);
   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   OMEGA_SCOPE(LocNormalVelocity, NormalVelocity);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
   OMEGA_SCOPE(LocBruntVaisalaFreqSq, BruntVaisalaFreqSq);
   OMEGA_SCOPE(LocIceFraction, IceFraction);
   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocBulkRichardson, BulkRichardsonNumber);
   const Real LocIceFracThresholdForMinOBL = IceFractionThresholdForMinimumOBL;
   const Real LocMinimumOBLUnderSeaIce     = MinimumOBLUnderSeaIce;

   deepCopy(BulkRichardsonNumber, 0.0_Real);

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          const Real u_star =
              Kokkos::fmax(0.0_Real, SurfaceFrictionVelocity(ICell));
          const Real b0 = SurfaceBuoyancyFlux(ICell);

          const I4 KMin     = MinLayerCell(ICell);
          const I4 KMax     = MaxLayerCell(ICell);
          const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
          const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

          // Keep OBL search active under ice to match MPAS behavior.
          // Only apply a minimum-depth constraint later via ConstrainOBLDepth.
          const Real iceFrac = LocIceFraction(ICell);

          // Full bulk-Richardson profile using MPAS-style surface-layer
          // averages for buoyancy and velocity.
          const Real rho_ref    = 1025.0_Real;
          Real obl_depth        = Kokkos::abs(ZInterface(ICell, KIntDeep));
          I4 k_cross            = -1;
          const Real ri_crit    = CriticalRichardson;
          const Real ri_scaling = 1.0_Real - 0.5_Real * SurfaceLayerExtent;
          const Real b0_eff     = b0 * LocLangmuirFactor(ICell);

          // CVMix default unresolved-shear constants.
          const Real c_s_unres = 24.0_Real * Kokkos::sqrt(17.0_Real);
          const Real vtc =
              Kokkos::sqrt(
                  0.2_Real /
                  Kokkos::max(1.0e-12_Real, c_s_unres * SurfaceLayerExtent)) /
              (VonKar * VonKar);

          I4 k_surface_avg     = KMin;
          const Real thick_top = Kokkos::abs(ZInterface(ICell, KMin + 1) -
                                             ZInterface(ICell, KMin));
          Real sum_thickness   = Kokkos::max(thick_top, 1.0e-12_Real);
          Real sum_rho = LocPotentialDensity(ICell, KMin) * sum_thickness;
          Real sum_u   = LocNormalVelocity(ICell, KMin) * sum_thickness;
          Real sum_v   = LocTangentialVelocity(ICell, KMin) * sum_thickness;

          for (I4 k = KMin; k <= KMax; ++k) {
             const I4 kCell      = Kokkos::min(k, NVertLayers - 1);
             const I4 kInt       = Kokkos::min(k + 1, NVertLayers);
             const Real z_depth  = Kokkos::abs(ZInterface(ICell, kInt));
             const Real z_center = Kokkos::abs(ZMid(ICell, kCell));
             if (z_depth < 1.0e-12)
                continue;

             const Real surf_layer_depth = SurfaceLayerExtent * z_depth;
             while (k_surface_avg < k &&
                    Kokkos::abs(ZInterface(ICell, k_surface_avg + 1)) <
                        surf_layer_depth) {
                ++k_surface_avg;
                const Real dk =
                    Kokkos::abs(ZInterface(ICell, k_surface_avg + 1) -
                                ZInterface(ICell, k_surface_avg));
                const Real thick_k = Kokkos::max(dk, 1.0e-12_Real);
                sum_thickness += thick_k;
                sum_rho += LocPotentialDensity(ICell, k_surface_avg) * thick_k;
                sum_u += LocNormalVelocity(ICell, k_surface_avg) * thick_k;
                sum_v += LocTangentialVelocity(ICell, k_surface_avg) * thick_k;
             }

             const Real inv_sum_thickness =
                 1.0_Real / Kokkos::max(sum_thickness, 1.0e-12_Real);
             const Real rho_avg_surf = sum_rho * inv_sum_thickness;
             const Real u_avg_surf   = sum_u * inv_sum_thickness;
             const Real v_avg_surf   = sum_v * inv_sum_thickness;

             const Real rho_k = LocPotentialDensity(ICell, kCell);
             const Real u_k   = LocNormalVelocity(ICell, kCell);
             const Real v_k   = LocTangentialVelocity(ICell, kCell);

             const Real delta_rho = rho_k - rho_avg_surf;
             const Real delta_b   = Gravity * delta_rho / rho_ref;

             const Real du     = u_k - u_avg_surf;
             const Real dv     = v_k - v_avg_surf;
             const Real shear2 = du * du + dv * dv;

             const Real sigma_loc = Kokkos::fmin(
                 1.0_Real, Kokkos::fmax(0.0_Real, SurfaceLayerExtent));

             Real w_turb = 0.0_Real;
             if (u_star > 1.0e-12_Real) {
                const Real u3        = u_star * u_star * u_star;
                const Real zeta      = sigma_loc * z_depth * VonKar * b0_eff /
                                       Kokkos::max(u3, 1.0e-20_Real);
                const Real phi_inv_s = KPP::KPPProfileS2(zeta);
                w_turb = VonKar * u_star * Kokkos::max(phi_inv_s, 0.0_Real);
             } else if (b0_eff < 0.0_Real) {
                const Real c_s = 16.0_Real;
                const Real ws3 = -c_s * sigma_loc * z_depth * VonKar * b0_eff;
                w_turb = VonKar * Kokkos::pow(Kokkos::max(ws3, 0.0_Real),
                                              1.0_Real / 3.0_Real);
             }
             const Real n_cntr = Kokkos::sqrt(
                 Kokkos::max(0.0_Real, LocBruntVaisalaFreqSq(ICell, kInt)));
             const Real cv  = (n_cntr < 0.002_Real)
                                  ? (2.1_Real - 200.0_Real * n_cntr)
                                  : 1.7_Real;
             const Real vt2 = Kokkos::max(
                 1.0e-10_Real, cv * vtc * z_center * n_cntr * w_turb /
                                   Kokkos::max(ri_crit, 1.0e-12_Real));

             const Real vel_scale2 = shear2 + vt2;

             const Real ri_b = ri_scaling * delta_b * z_center /
                               Kokkos::max(vel_scale2, 1.0e-12_Real);
             LocBulkRichardson(ICell, kInt) = ri_b;

             if (k_cross < 0 && ri_b >= ri_crit) {
                k_cross = k;
             }
          }

          if (k_cross >= KMin) {
             const I4 kIntCross = Kokkos::min(k_cross + 1, NVertLayers);
             if (k_cross > KMin) {
                const Real z_above  = Kokkos::abs(ZInterface(ICell, k_cross));
                const Real z_below  = Kokkos::abs(ZInterface(ICell, kIntCross));
                const Real ri_above = LocBulkRichardson(ICell, k_cross);
                const Real ri_below = LocBulkRichardson(ICell, kIntCross);
                const Real d_ri     = ri_below - ri_above;
                if (Kokkos::abs(d_ri) > 1.0e-12_Real) {
                   const Real frac = Kokkos::fmax(
                       0.0_Real,
                       Kokkos::fmin(1.0_Real, (ri_crit - ri_above) / d_ri));
                   obl_depth = z_above + frac * (z_below - z_above);
                } else {
                   obl_depth = z_below;
                }
             } else {
                obl_depth = Kokkos::abs(ZInterface(ICell, kIntCross));
             }
          } else {
             obl_depth = Kokkos::abs(ZInterface(ICell, KIntDeep));
          }

          // NOTE: MatchBoth in CVMix controls boundary-layer/interior shape
          // matching for coefficients, not OBL-depth blending. Keep OBL depth
          // set solely by the bulk-Richardson crossing and physical bounds.

          const Real top_layer_thickness =
              Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
          const Real min_obl_depth = 0.5_Real * top_layer_thickness;
          const Real water_depth   = Kokkos::abs(ZInterface(ICell, KIntDeep));
          // Apply physical bounds using configurable thresholds.
          obl_depth = Kokkos::fmax(obl_depth, min_obl_depth);
          if (iceFrac > LocIceFracThresholdForMinOBL) {
             obl_depth = Kokkos::fmax(obl_depth, LocMinimumOBLUnderSeaIce);
          }
          obl_depth = Kokkos::fmin(obl_depth, 0.95_Real * water_depth);

          I4 k_final = KMax;
          for (I4 k = KMin; k < KMax; ++k) {
             const Real z_above = Kokkos::abs(ZInterface(ICell, k));
             const Real z_below = Kokkos::abs(ZInterface(ICell, k + 1));
             if (obl_depth >= z_above && obl_depth <= z_below) {
                k_final = k;
                break;
             }
          }

          LocBoundaryLayerDepth(ICell)      = obl_depth;
          LocIndexBoundaryLayerDepth(ICell) = k_final;
       });

   LOG_INFO("KPPMix::computeOBLDepth: OBL depth computed");
}

/// Stage 2: Compute mixing coefficients within and below OBL
void KPPMix::computeMixingCoefficients(
    const Array2DReal &PotentialDensity, const Array2DReal &NormalVelocity,
    const Array2DReal &TangentialVelocity,
    const Array1DReal &SurfaceFrictionVelocity,
    const Array1DReal &SurfaceBuoyancyFlux) {

   using namespace KPP;

   (void)PotentialDensity;
   (void)NormalVelocity;
   (void)TangentialVelocity;

   I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Capture data for use in lambda
   // =======================================================================
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocTurbulentVelocityScale, TurbulentVelocityScale);
   OMEGA_SCOPE(LocSurfaceFrictionVelocity, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfaceBuoyancyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(ZInterface, VCoord->GeomZInterface);

   // Capture member variables for use in lambda
   Real LocBackgroundDiff  = BackgroundDiff;
   Real LocBackgroundVisc  = BackgroundVisc;
   bool LocUseNonLocalFlux = UseNonLocalFlux;
   I4 LocMatchTechnique = 0; // 0=SimpleShapes, 1=MatchBoth, 2=ParabolicNonLocal
   if (MatchTechniqueStr == "MatchBoth") {
      LocMatchTechnique = 1;
   } else if (MatchTechniqueStr == "ParabolicNonLocal") {
      LocMatchTechnique = 2;
   }
   bool LocUseEnhancedDiffusion = UseEnhancedDiffusion;
   I4 LocInterpType2            = 3; // 0=Linear, 1=Quadratic, 2=Cubic, 3=LMD94
   if (InterpType2Str == "Linear" || InterpType2Str == "linear") {
      LocInterpType2 = 0;
   } else if (InterpType2Str == "Quadratic" || InterpType2Str == "quadratic") {
      LocInterpType2 = 1;
   } else if (InterpType2Str == "Cubic" || InterpType2Str == "cubic") {
      LocInterpType2 = 2;
   } else {
      LocInterpType2 = 3;
   }
   const Real LocKappa = VonKar;

   // =======================================================================
   // Initialize with background mixing
   // =======================================================================
   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          LocVertDiff(ICell, K)               = LocBackgroundDiff;
          LocVertVisc(ICell, K)               = LocBackgroundVisc;
          LocVertNonLocalFlux(ICell, K)       = 0.0;
          LocTurbulentVelocityScale(ICell, K) = 0.0;
       });

   // =======================================================================
   // Stage 2: Compute KPP profile-based mixing coefficients
   // =======================================================================

   parallelFor(
       "KPP-MixingCoeffs", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          // Get OBL information for this cell
          Real h_obl = LocBoundaryLayerDepth(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          // =============================================================
          // Compute turbulent velocity scale w_s
          // =============================================================
          Real u_star = LocSurfaceFrictionVelocity(ICell);
          Real b0     = LocSurfaceBuoyancyFlux(ICell);

          // =============================================================
          // Compute mixing coefficients at each interface
          // =============================================================
          for (I4 k = KMin; k <= KMax + 1; ++k) {
             const I4 k_iface   = Kokkos::min(Kokkos::max(k, 0), NVertLayers);
             const Real z_depth = Kokkos::abs(ZInterface(ICell, k_iface));

             // Check if within OBL using geometric depth.
             if (z_depth <= h_obl && h_obl > 0.0_Real) {
                // Normalized depth in Omega sign convention: sigma in [-1,0].
                Real sigma = -z_depth / h_obl;
                sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma));

                // CVMix-style turbulent scales: w = kappa*u*/phi in general,
                // with explicit free-convection limits when u*=0.
                const Real sigma_coord = -sigma; // [0,1]
                const Real sigma_loc   = Kokkos::fmin(
                    SurfaceLayerExtent, Kokkos::fmax(0.0_Real, sigma_coord));

                Real zeta     = 0.0_Real;
                Real w_m_turb = 0.0_Real;
                Real w_s_turb = 0.0_Real;

                if (u_star > 0.0_Real) {
                   const Real u3 = u_star * u_star * u_star;
                   zeta          = sigma_loc * h_obl * b0 * LocKappa /
                                   Kokkos::max(u3, 1.0e-20_Real);

                   // KPPProfileM2/S2 return phi^{-1}; do not invert again.
                   const Real phi_inv_m = KPP::KPPProfileM2(zeta);
                   const Real phi_inv_s = KPP::KPPProfileS2(zeta);

                   w_m_turb =
                       LocKappa * u_star * Kokkos::max(phi_inv_m, 0.0_Real);
                   w_s_turb =
                       LocKappa * u_star * Kokkos::max(phi_inv_s, 0.0_Real);
                } else if (b0 < 0.0_Real) {
                   // Free-convection edge case (u*=0, unstable forcing).
                   const Real c_m = 16.0_Real;
                   const Real c_s = 16.0_Real;
                   const Real wm3 = -c_m * sigma_loc * h_obl * LocKappa * b0;
                   const Real ws3 = -c_s * sigma_loc * h_obl * LocKappa * b0;
                   w_m_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, wm3),
                                                     1.0_Real / 3.0_Real);
                   w_s_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, ws3),
                                                     1.0_Real / 3.0_Real);
                }

                // ========================================================
                // Momentum mixing: K_m = u_star * w_s * M1(σ) * M2(zeta)
                // ========================================================
                Real m1 = KPP::KPPProfileM1(sigma);
                LocVertVisc(ICell, k) =
                    Kokkos::fmax(KPP::MIN_COEFFICIENT, h_obl * w_m_turb * m1);

                // ========================================================
                // Tracer mixing: K_s = u_star * w_s * S1(σ) * S2(zeta)
                // ========================================================
                Real s1 = KPP::KPPProfileS1(sigma);
                LocVertDiff(ICell, k) =
                    Kokkos::fmax(KPP::MIN_COEFFICIENT, h_obl * w_s_turb * s1);
                LocTurbulentVelocityScale(ICell, k) = w_s_turb;

                // ========================================================
                // Non-local flux: G(σ)
                // ========================================================
                if (LocUseNonLocalFlux) {
                   Real g_sigma = 0.0_Real;
                   if (LocMatchTechnique == 2) {
                      g_sigma = KPP::KPPProfileGParabolicNonLocal(sigma);
                   } else if (LocMatchTechnique == 1) {
                      g_sigma = KPP::KPPProfileGMatchBoth(sigma);
                   } else {
                      g_sigma = KPP::KPPProfileG(sigma);
                   }
                   LocVertNonLocalFlux(ICell, k) = g_sigma;
                } else {
                   LocVertNonLocalFlux(ICell, k) = 0.0;
                }

             } else {
                // Below OBL: use background values
                LocVertDiff(ICell, k)               = LocBackgroundDiff;
                LocVertVisc(ICell, k)               = LocBackgroundVisc;
                LocVertNonLocalFlux(ICell, k)       = 0.0;
                LocTurbulentVelocityScale(ICell, k) = 0.0;
             }
          }

          // Optional enhanced diffusion/viscosity treatment at OBL base.
          // This keeps OMEGA C++-native while approximating the CVMix-style
          // increased base mixing with configurable interpolation behavior.
          if (LocUseEnhancedDiffusion && h_obl > 0.0_Real) {
             const I4 k_obl = Kokkos::max(
                 KMin, Kokkos::min(LocIndexBoundaryLayerDepth(ICell), KMax));
             const I4 k_base_iface  = Kokkos::min(k_obl + 1, KMax + 1);
             const I4 k_above_iface = Kokkos::max(KMin, k_base_iface - 1);
             const I4 k_below_iface = Kokkos::min(k_base_iface + 1, KMax + 1);

             const Real z_above = Kokkos::abs(ZInterface(ICell, k_above_iface));
             const Real z_base  = Kokkos::abs(ZInterface(ICell, k_base_iface));
             const Real dz_base = Kokkos::max(z_base - z_above, 1.0e-12_Real);
             const Real frac    = Kokkos::fmax(
                 0.0_Real, Kokkos::fmin(1.0_Real, (h_obl - z_above) / dz_base));

             const Real diff_above = LocVertDiff(ICell, k_above_iface);
             const Real visc_above = LocVertVisc(ICell, k_above_iface);
             const Real diff_base  = LocVertDiff(ICell, k_base_iface);
             const Real visc_base  = LocVertVisc(ICell, k_base_iface);

             Real diff_interp =
                 (1.0_Real - frac) * diff_above + frac * diff_base;
             Real visc_interp =
                 (1.0_Real - frac) * visc_above + frac * visc_base;

             // Quadratic/cubic path uses a 3-point Lagrange fit when available.
             if ((LocInterpType2 == 1 || LocInterpType2 == 2) &&
                 k_below_iface > k_base_iface) {
                const Real x          = frac;
                const Real diff_below = LocVertDiff(ICell, k_below_iface);
                const Real visc_below = LocVertVisc(ICell, k_below_iface);

                const Real l0 = ((x - 1.0_Real) * (x - 2.0_Real)) / 2.0_Real;
                const Real l1 = -x * (x - 2.0_Real);
                const Real l2 = (x * (x - 1.0_Real)) / 2.0_Real;

                diff_interp =
                    l0 * diff_above + l1 * diff_base + l2 * diff_below;
                visc_interp =
                    l0 * visc_above + l1 * visc_base + l2 * visc_below;
             }

             LocVertDiff(ICell, k_base_iface) =
                 Kokkos::max(LocVertDiff(ICell, k_base_iface),
                             Kokkos::max(diff_interp, KPP::MIN_COEFFICIENT));
             LocVertVisc(ICell, k_base_iface) =
                 Kokkos::max(LocVertVisc(ICell, k_base_iface),
                             Kokkos::max(visc_interp, KPP::MIN_COEFFICIENT));
          }
       });

   LOG_INFO("KPPMix::computeMixingCoefficients: Phase 2 mixing coefficients "
            "computed");
}

/// Register fields with I/O system
void KPPMix::defineFields() {
   const Real FillValue = -9.99e30;

   // BoundaryLayerDepth on cells
   std::vector<std::string> CellDims(1);
   CellDims[0] = "NCells";
   auto OBLDepthField =
       Field::create(OBLDepthFldName,                  // field name
                     "ocean boundary layer depth",     // long name
                     "m",                              // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue,                        // fill value
                     1,                                // number of dims
                     CellDims);

   // KPP non-local tracer flux profile on cell-layer interfaces
   std::vector<std::string> FluxDims(2);
   FluxDims[0] = "NCells";
   FluxDims[1] = "NVertLayersP1";
   auto NonLocalFluxField =
       Field::create(NonLocalFluxFldName,                 // field name
                     "KPP non-local tracer flux profile", // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     FillValue,                           // fill value
                     2,                                   // number of dims
                     FluxDims);

   auto BulkRichardsonField =
       Field::create(BulkRichardsonFldName,               // field name
                     "bulk Richardson number",            // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     FillValue,                           // fill value
                     2,                                   // number of dims
                     FluxDims);

   auto TurbulentVelScaleField =
       Field::create(TurbulentVelScaleFldName,         // field name
                     "KPP turbulent velocity scale",   // long name
                     "m s-1",                          // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue,                        // fill value
                     2,                                // number of dims
                     FluxDims);

   // Group KPP-specific outputs for convenient stream selection.
   auto KPPGroup = FieldGroup::create("KPPMix");
   KPPGroup->addField(OBLDepthFldName);
   KPPGroup->addField(NonLocalFluxFldName);
   KPPGroup->addField(BulkRichardsonFldName);
   KPPGroup->addField(TurbulentVelScaleFldName);

   OBLDepthField->attachData<Array1DReal>(BoundaryLayerDepth);
   NonLocalFluxField->attachData<Array2DReal>(VertNonLocalFlux);
   BulkRichardsonField->attachData<Array2DReal>(BulkRichardsonNumber);
   TurbulentVelScaleField->attachData<Array2DReal>(TurbulentVelocityScale);

   LOG_INFO("KPPMix::defineFields: registered {}, {}, {}, {}", OBLDepthFldName,
            NonLocalFluxFldName, BulkRichardsonFldName,
            TurbulentVelScaleFldName);
}

} // namespace OMEGA

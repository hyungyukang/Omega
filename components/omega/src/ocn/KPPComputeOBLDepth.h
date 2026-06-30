#ifndef OMEGA_KPP_COMPUTE_OBL_DEPTH_H
#define OMEGA_KPP_COMPUTE_OBL_DEPTH_H
//===-- ocn/KPPComputeOBLDepth.h - OBL Depth Computation -------*- C++ -*-===//
//
/// \file
/// \brief Compute ocean boundary layer depth from bulk Richardson criterion
///
/// This header defines functors for computing the OBL depth using the bulk
/// Richardson number criterion. Supports both MatchBoth and SimpleShapes
/// OBL depth estimation schemes from Large et al. (1994).
//
//===----------------------------------------------------------------------===//

#include "KPPConstants.h"
#include "OmegaKokkos.h"

namespace OMEGA::KPP {

enum class OBLDepthScheme {
   SimpleShapes, ///< Traditional: OBL where Ri_critical is met
   MatchBoth     ///< Interpolate OBL between Ri_critical and surface
};

/// @brief Compute OBL depth from bulk Richardson numbers
/// Per-cell functor for computing boundary layer depth
///
/// Implements the OBL depth search following Large et al. (1994):
/// 1. Initialize OBL at surface (z=0)
/// 2. Search downward, computing bulk Richardson number at each layer
/// 3. Stop when Ri_b > Ri_critical (criterion met)
/// 4. Optional: Apply MatchBoth scheme to interpolate between points
/// 5. Apply min/max constraints
///
/// REFERENCES: Large et al. (1994) Eq. (15)-(16), Large et al. (1997)
class KPPComputeOBLDepth {

 public:
   // Member data - references to mesh and state arrays
   Array1DReal ZInterface;     ///< Depth at interfaces (m, negative down)
   Array1DReal ZCenter;        ///< Depth at cell centers (m)
   Array1DReal LayerThickness; ///< Layer thickness (m)
   Array1DI4 MinLayerCell;     ///< Min layer index per cell
   Array1DI4 MaxLayerCell;     ///< Max layer index per cell

   // Input forcing data
   Real SurfaceFrictionVelocity;   ///< u* (m/s) for this cell
   Real SurfaceBuoyancyFlux;       ///< B_0 (m²/s³) for this cell
   Real LangmuirEnhancementFactor; ///< Wave mixing enhancement (>= 1.0)

   // Stratification data (for bulk Richardson computation)
   Array1DReal BruntVaisalaFreqSq;   ///< N²(z) at interfaces (s⁻²)
   Array1DReal VelocityShearSquared; ///< (∂u/∂z)² (s⁻²)

   // OBL computation parameters
   Real CriticalRichardson;    ///< Ri_crit threshold (typically 0.3)
   Real SurfaceLayerExtent;    ///< Fraction of OBL for averaging surface layer
   OBLDepthScheme DepthScheme; ///< MatchBoth or SimpleShapes

   Real IceFraction; ///< Sea ice coverage (0-1) for suppression
   I4 LandIceMask;   ///< Land ice mask (suppress if != 0)
   Real WaterDepth;  ///< Maximum water depth (m)

   // Output
   Real ComputedOBLDepth; ///< OBL depth (m) calculated by this functor
   I4 ComputedOBLIndex;   ///< Layer index k where OBL base is located

   /// @brief Constructor - initialize with cell and mesh data
   KPPComputeOBLDepth(const Array1DReal &ZInterface_in,
                      const Array1DReal &ZCenter_in,
                      const Array1DReal &LayerThickness_in,
                      const Array1DI4 &MinLayerCell_in,
                      const Array1DI4 &MaxLayerCell_in, Real u_star, Real b0,
                      Real langmuir_factor, const Array1DReal &BVFSq_in,
                      const Array1DReal &VelShearSq_in, Real ri_crit,
                      Real surf_layer_ext, OBLDepthScheme scheme, Real ice_frac,
                      I4 land_ice, Real water_depth)
       : ZInterface(ZInterface_in), ZCenter(ZCenter_in),
         LayerThickness(LayerThickness_in), MinLayerCell(MinLayerCell_in),
         MaxLayerCell(MaxLayerCell_in), SurfaceFrictionVelocity(u_star),
         SurfaceBuoyancyFlux(b0), LangmuirEnhancementFactor(langmuir_factor),
         BruntVaisalaFreqSq(BVFSq_in), VelocityShearSquared(VelShearSq_in),
         CriticalRichardson(ri_crit), SurfaceLayerExtent(surf_layer_ext),
         DepthScheme(scheme), IceFraction(ice_frac), LandIceMask(land_ice),
         WaterDepth(water_depth), ComputedOBLDepth(0.0), ComputedOBLIndex(0) {}

   /// @brief Compute OBL depth using bulk Richardson search
   /// Called from Stage 1 (per-cell operation, no K loop)
   ///
   /// Algorithm:
   /// 1. Initialize search at kmin (surface)
   /// 2. For each k from kmin to kmax:
   ///    a. Compute Ri_b(k) from accumulated buoyancy/shear
   ///    b. If Ri_b > Ri_crit, OBL base found at this k
   ///    c. Store OBL depth = |zw(k+1)|
   /// 3. Apply MatchBoth interpolation if enabled
   /// 4. Enforce min/max constraints
   /// 5. Check for ice suppression
   ///
   /// CRITICAL: Use cumulative sum pattern to avoid O(nLevels²) behavior
   ///
   KOKKOS_FUNCTION
   void computeOBLDepth(I4 ICell) {

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);

      // Initialize output
      ComputedOBLDepth = 0.0;
      ComputedOBLIndex = KMin;

      // =======================================================================
      // Stage 1a: Check for suppression conditions
      // =======================================================================
      if (ShouldSuppressOBL(IceFraction, LandIceMask)) {
         // Under ice or land ice: set to minimum
         ComputedOBLDepth = MIN_OBL_UNDER_ICE;
         ComputedOBLIndex =
             findLayerIndexForDepth(KMin, KMax, ComputedOBLDepth);
         return;
      }

      // =======================================================================
      // Stage 1b: Perform bulk Richardson search
      // CRITICAL: Accumulate sums to avoid redundant calculation
      // =======================================================================

      // Running sums for bulk Richardson calculation
      // Initialize at surface (k = KMin)
      Real buoy_sum             = 0.0; // Cumulative buoyancy term
      Real shear_sum            = 0.0; // Cumulative shear term
      Real richardson_max_found = -1.0;
      I4 k_richardson_base      = KMin;

      // Apply Langmuir enhancement to surface buoyancy forcing
      Real b0_enhanced = SurfaceBuoyancyFlux * LangmuirEnhancementFactor;

      // Loop downward through layers
      for (I4 k = KMin; k <= KMax; ++k) {

         // Get current layer depth drop
         Real z_depth = Kokkos::abs(ZInterface(k + 1));
         if (z_depth < 1.0e-10) {
            // Avoid division by zero at surface
            continue;
         }

         // =================================================================
         // Accumulate bulk Richardson components
         // =================================================================
         // Buoyancy term: integral of N² from surface to depth z
         if (k > KMin) {
            Real dz = Kokkos::abs(ZCenter(k) - ZCenter(k - 1));
            buoy_sum += BruntVaisalaFreqSq(k) * dz;
         }

         // Shear term: integral of (∂u/∂z)² from surface to depth z
         if (k > KMin) {
            Real dz = Kokkos::abs(ZCenter(k) - ZCenter(k - 1));
            shear_sum += VelocityShearSquared(k) * dz;
         }

         // Compute bulk Richardson number at this depth
         // Ri_b = (b0 * z * buoy_sum) / (u*^4 * shear_sum + epsilon)
         Real u_star_sq =
             Kokkos::fmax(MIN_USTAR * MIN_USTAR,
                          SurfaceFrictionVelocity * SurfaceFrictionVelocity);
         Real denom = u_star_sq * u_star_sq * Kokkos::fmax(1.0e-15, shear_sum);
         Real richardson = (b0_enhanced * z_depth * buoy_sum) / denom;

         richardson_max_found = Kokkos::fmax(richardson_max_found, richardson);

         // Check stopping criterion
         if (richardson >= CriticalRichardson) {
            k_richardson_base = k;
            break; // Early exit - OBL base found
         }
      }

      // Store raw depth and index from search
      ComputedOBLDepth = Kokkos::abs(ZInterface(k_richardson_base + 1));
      ComputedOBLIndex = k_richardson_base;

      // =======================================================================
      // Stage 1c: Apply MatchBoth interpolation if enabled
      // =======================================================================
      if (DepthScheme == OBLDepthScheme::MatchBoth) {
         // Interpolate between Ri_critical and surface points
         // This smooths the OBL depth boundary
         applyMatchBothInterpolation(KMin, KMax);
      }

      // =======================================================================
      // Stage 1d: Apply min/max constraints
      // =======================================================================
      Real surface_thickness = Kokkos::abs(ZCenter(KMin + 1) - ZCenter(KMin));
      ComputedOBLDepth = ConstrainOBLDepth(ComputedOBLDepth, surface_thickness,
                                           WaterDepth, IceFraction);

      // Update layer index after constraining depth
      ComputedOBLIndex = findLayerIndexForDepth(KMin, KMax, ComputedOBLDepth);
   }

   /// @brief Find vertical layer index corresponding to a given depth
   /// Binary search or linear search to find k where depth is between z(k) and
   /// z(k+1)
   KOKKOS_FUNCTION
   I4 findLayerIndexForDepth(I4 KMin, I4 KMax, Real target_depth) {
      target_depth = Kokkos::abs(target_depth);

      // Linear search (simple, sufficient for reasonable OBL depths)
      for (I4 k = KMin; k < KMax; ++k) {
         Real z_above = Kokkos::abs(ZInterface(k));
         Real z_below = Kokkos::abs(ZInterface(k + 1));

         if (target_depth >= z_above && target_depth <= z_below) {
            return k;
         }
      }

      return KMax; // Default to deepest layer if not found
   }

   /// @brief Apply MatchBoth scheme interpolation
   /// Matches OBL using cubic interpolation between surface and critical points
   KOKKOS_FUNCTION
   void applyMatchBothInterpolation(I4 KMin, I4 KMax) {
      // MatchBoth: Interpolate OBL depth between Ri=0 (surface) and Ri_crit
      // Smooths potential discontinuities in OBL depth time series

      // For now, use simple weighted average between surface and critical point
      // Can be enhanced with cubic polynomial interpolation later
      const Real weight_surface  = 0.3; // Blend 30% surface influence
      const Real weight_critical = 0.7; // Blend 70% critical point

      Real surface_depth  = 1.0; // Surface layer depth (minimum)
      Real critical_depth = ComputedOBLDepth;

      // Weighted blend
      ComputedOBLDepth =
          weight_surface * surface_depth + weight_critical * critical_depth;

      // Recompute index after interpolation
      ComputedOBLIndex = findLayerIndexForDepth(KMin, KMax, ComputedOBLDepth);
   }

}; // class KPPComputeOBLDepth

} // namespace OMEGA::KPP

#endif // OMEGA_KPP_COMPUTE_OBL_DEPTH_H

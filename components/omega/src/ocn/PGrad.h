#ifndef OMEGA_PGRAD_H
#define OMEGA_PGRAD_H
//===-- ocn/PGrad.h - Pressure Gradient -----------------*- C++ -*-===//
///
/// Implements the PressureGrad class which provides a centered,
/// TEOS-10 pressure-integrated, high-order pressure gradient option
/// and dispatches computations to functor objects. This follows the
/// patterns used in Eos.h/Eos.cpp.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <map>
#include <memory>
#include <string>

namespace OMEGA {

enum class PressureGradType {
   Centered,
   Integrated,
   HighOrder1,
   HighOrder2
};

//------------------------------------------------------------------------------
// Centered pressure gradient functor
class PressureGradCentered {
 public:
   bool Enabled;

   // constructor declaration
   PressureGradCentered(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                        const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   // Compute centered pressure gradient contribution for given edge and
   // vertical chunk. This appends results into the Tend array (in-place).
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &PressureMid,
                                   const Array2DReal &PressureInterface,
                                   const Array2DReal &GeomZInterface,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1.0_Real / DcEdge(IEdge);

      Real GradGeoPot =
          (TidalPotential(ICell1) - TidalPotential(ICell0)) * InvDcEdge +
          (SelfAttractionLoading(ICell1) - SelfAttractionLoading(ICell0)) *
              InvDcEdge;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Real MontPotCell0K =
             PressureInterface(ICell0, K) * SpecVol(ICell0, K) +
             Gravity * GeomZInterface(ICell0, K);
         Real MontPotCell1K =
             PressureInterface(ICell1, K) * SpecVol(ICell1, K) +
             Gravity * GeomZInterface(ICell1, K);
         Real GradMontPotK = (MontPotCell1K - MontPotCell0K) * InvDcEdge;

         Real MontPotCell0Kp1 =
             PressureInterface(ICell0, K + 1) * SpecVol(ICell0, K) +
             Gravity * GeomZInterface(ICell0, K + 1);
         Real MontPotCell1Kp1 =
             PressureInterface(ICell1, K + 1) * SpecVol(ICell1, K) +
             Gravity * GeomZInterface(ICell1, K + 1);
         Real GradMontPotKp1 = (MontPotCell1Kp1 - MontPotCell0Kp1) * InvDcEdge;
         Real GradMontPot    = 0.5_Real * (GradMontPotK + GradMontPotKp1);

         Real PGradAlpha =
             0.5_Real * (PressureMid(ICell1, K) + PressureMid(ICell0, K)) *
             (SpecVol(ICell1, K) - SpecVol(ICell0, K)) * InvDcEdge;
         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * (-GradMontPot + PGradAlpha - GradGeoPot);
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};


//------------------------------------------------------------------------------
// TEOS-10 pressure-integrated pressure-gradient functor
//
// At an edge, the pressure/geopotential increment at each layer interface is
//
//   Integral_{p0}^{p1} alpha(SA_e, CT_e, p) dp + g (z1 - z0).
//
// This avoids forming alpha*grad(p) and grad(Phi) as two independently
// discretized large terms.  SA_e and CT_e are centered edge values.
class PressureGradIntegrated {
 public:
   bool Enabled = false;

   PressureGradIntegrated(
       const HorzMesh *Mesh,   ///< [in] Horizontal mesh
       const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   KOKKOS_FUNCTION void
   operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
              const Array2DReal &PressureInterface,
              const Array2DReal &GeomZInterface,
              const Array1DReal &TidalPotential,
              const Array1DReal &SelfAttractionLoading,
              const Array2DReal &ConservTemp,
              const Array2DReal &AbsSalinity) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1.0_Real / DcEdge(IEdge);

      const Real GradGeoPot =
          (TidalPotential(ICell1) - TidalPotential(ICell0)) * InvDcEdge +
          (SelfAttractionLoading(ICell1) -
           SelfAttractionLoading(ICell0)) *
              InvDcEdge;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         // First centered implementation: reconstruct one thermodynamic state
         // at the edge and integrate only through pressure space.
         const Real CtEdge =
             0.5_Real * (ConservTemp(ICell0, K) + ConservTemp(ICell1, K));
         const Real SaEdge =
             0.5_Real * (AbsSalinity(ICell0, K) + AbsSalinity(ICell1, K));

         const Real IntAlphaDpTop = calcSpecVolPressureIntegral(
             SaEdge, CtEdge, PressureInterface(ICell0, K),
             PressureInterface(ICell1, K));

         const Real IntAlphaDpBot = calcSpecVolPressureIntegral(
             SaEdge, CtEdge, PressureInterface(ICell0, K + 1),
             PressureInterface(ICell1, K + 1));

         const Real DeltaPhiTop =
             Gravity * (GeomZInterface(ICell1, K) -
                        GeomZInterface(ICell0, K));
         const Real DeltaPhiBot =
             Gravity * (GeomZInterface(ICell1, K + 1) -
                        GeomZInterface(ICell0, K + 1));

         const Real ResidualTop = IntAlphaDpTop + DeltaPhiTop;
         const Real ResidualBot = IntAlphaDpBot + DeltaPhiBot;

         const Real PGrad =
             -0.5_Real * (ResidualTop + ResidualBot) * InvDcEdge;

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * (PGrad - GradGeoPot);
      }
   }

 private:
   /// Analytic integral of the Omega TEOS-10 75-term specific-volume
   /// polynomial with respect to gauge pressure in Pa.
   ///
   /// Return units: J kg^-1 = m^2 s^-2.
   KOKKOS_FUNCTION Real
   calcSpecVolPressureIntegral(const Real Sa, const Real Ct,
                               const Real P0Pa, const Real P1Pa) const {

      Real SpecVolPCoeffs[6 * VecLength];
      Teos10Evaluator.calcPCoeffs(SpecVolPCoeffs, 0, Ct, Sa);

      // Reference-pressure-profile coefficients used in Teos10Eos::
      // calcRefProfile().  The normalized pressure is
      // X = pressure(Pa) * 1.e-8.
      constexpr Real V00 = -4.4015007269e-05;
      constexpr Real V01 =  6.9232335784e-06;
      constexpr Real V02 = -7.5004675975e-07;
      constexpr Real V03 =  1.7009109288e-08;
      constexpr Real V04 = -1.6884162004e-08;
      constexpr Real V05 =  1.9613503930e-09;

      // alpha(X) = A0 + A1*X + ... + A6*X^6.
      Real A[7];
      A[0] = SpecVolPCoeffs[0];
      A[1] = SpecVolPCoeffs[1] + V00;
      A[2] = SpecVolPCoeffs[2] + V01;
      A[3] = SpecVolPCoeffs[3] + V02;
      A[4] = SpecVolPCoeffs[4] + V03;
      A[5] = SpecVolPCoeffs[5] + V04;
      A[6] = V05;

      constexpr Real PaToNormalizedPressure = 1.0e-8;
      constexpr Real NormalizedPressureToPa = 1.0e8;

      const Real X0 = P0Pa * PaToNormalizedPressure;
      const Real X1 = P1Pa * PaToNormalizedPressure;
      const Real DX = X1 - X0;

      // Stable evaluation of
      // (X1^(N+1) - X0^(N+1)) / (X1 - X0).
      Real PowerQuotient = 1.0_Real;
      Real X0Power       = 1.0_Real;
      Real IntegralMean  = A[0];

      for (int N = 1; N <= 6; ++N) {
         X0Power *= X0;
         PowerQuotient = X1 * PowerQuotient + X0Power;
         IntegralMean +=
             A[N] * PowerQuotient / static_cast<Real>(N + 1);
      }

      return NormalizedPressureToPa * DX * IntegralMean;
   }

   Teos10Eos Teos10Evaluator;
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

//------------------------------------------------------------------------------
// High-order pressure gradient functor (placeholder)
class PressureGradHighOrder {
 public:
   bool Enabled;

   // constructor declaration
   PressureGradHighOrder(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                         const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &PressureMid,
                                   const Array2DReal &PressureInterface,
                                   const Array2DReal &GeomZInterface,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading,
                                   const Array2DReal &SpecVol) const {

      // Placeholder: for now, no-op (future high-order implementation)
      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += 0.0_Real;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

//------------------------------------------------------------------------------
// Pressure gradient manager class
class PressureGrad {
 public:
   // Flag to indicate if pressure gradient term is enabled
   bool Enabled;

   // Initialize the default instance
   static void init();

   // Create a new pressure gradient object and add to map
   static PressureGrad *create(const std::string &Name, const HorzMesh *Mesh,
                               const VertCoord *VCoord, Config *Options);

   // Get the default instance
   static PressureGrad *getDefault();

   // Get instance by name
   static PressureGrad *
   get(const std::string &Name ///< [in] Name of PressureGrad
   );

   // Deallocates arrays and deletes instance
   static void clear();

   // Remove pressure gradient object by name
   static void erase(const std::string &Name ///< [in] Name of PressureGrad
   );

   // Destructor
   ~PressureGrad();

   // Compute pressure gradient tendencies and add into Tend array
   void computePressureGrad(Array2DReal &Tend, const Array2DReal &PressureMid,
                            const Array2DReal &PressureInterface,
                            const Array2DReal &SpecVol,
                            const Array2DReal &GeomZInterface,
                            const Array2DReal &PseudoThick,
                            const Array2DReal &ConservTemp,
                            const Array2DReal &AbsSalinity) const;

 private:
   // Construct a new pressure gradient object
   PressureGrad(const HorzMesh *Mesh, const VertCoord *VCoord, Config *Options);

   // forbid copy and move construction
   PressureGrad(const PressureGrad &) = delete;
   PressureGrad(PressureGrad &&)      = delete;

   // Pointer to default pressure gradient object
   static PressureGrad *DefaultPGrad;

   // Mesh-related sizes
   I4 NEdgesAll     = 0;
   I4 NEdgesOwned   = 0;
   I4 NVertLayers   = 0;
   I4 NVertLayersP1 = 0;

   // Data required for computation (stored copies of VCoord arrays)
   Array1DI4 MinLayerEdgeBot; ///< min vertical layer on each edge
   Array1DI4 MaxLayerEdgeTop; ///< max vertical layer on each edge

   // Temporary: to be moveed to tidal forcing module in future
   Array1DReal TidalPotential; ///< Tidal potential for tidal forcing
   Array1DReal
       SelfAttractionLoading; ///< Self attraction and loading for tidal forcing

   // Instances of functors
   PressureGradCentered CenteredPGrad;
   PressureGradIntegrated IntegratedPGrad;
   PressureGradHighOrder HighOrderPGrad;

   // Choice from config
   PressureGradType PressureGradChoice = PressureGradType::Centered;

   // Map of all pressure gradient objects by name
   static std::map<std::string, std::unique_ptr<PressureGrad>> AllPGrads;

}; // end class PressureGrad

} // namespace OMEGA
#endif

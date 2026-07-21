#ifndef OMEGA_PGRAD_H
#define OMEGA_PGRAD_H
//===-- ocn/PGrad.h - Pressure Gradient -----------------*- C++ -*-===//
///
/// Implements the PressureGrad class which provides a centered and
/// high-order pressure gradient option and dispatches computations to
/// functor objects. This follows the patterns used in Eos.h/Eos.cpp.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <memory>

namespace OMEGA {

enum class PressureGradType { Centered, HighOrder1, HighOrder2 };

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
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &GeomZMid,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading) const {

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
         const Real SpecVolEdge =
             0.5_Real * (SpecVol(ICell0, K) + SpecVol(ICell1, K));
         const Real GradPressure =
             (PressureMid(ICell1, K) - PressureMid(ICell0, K)) * InvDcEdge;
         const Real GradGeomZ =
             (GeomZMid(ICell1, K) - GeomZMid(ICell0, K)) * InvDcEdge;

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) *
             (-SpecVolEdge * GradPressure - Gravity * GradGeomZ - GradGeoPot);
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// High-order pressure gradient functor
class PressureGradHighOrder {
 public:
   bool Enabled;

   // constructor declaration
   PressureGradHighOrder(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                         const VertCoord *VCoord ///< [in] Vertical coordinate
   );

   /// Reconstruct hydrostatically consistent geometric height at the
   /// quadrature points used by the pressure-gradient integration.
   void computeGeometry(const Array2DReal &PressureInterface,
                        const Array2DReal &ConservTemp,
                        const Array2DReal &AbsSalinity,
                        const Array2DReal &SpecVol,
                        EosType EosChoice) const;

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &PressureInterface,
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading,
                                   EosType EosChoice) const {
      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1.0_Real / DcEdge(IEdge);

      const Real GradGeoPot =
          ((TidalPotential(ICell1) - TidalPotential(ICell0)) +
           (SelfAttractionLoading(ICell1) -
            SelfAttractionLoading(ICell0))) *
          InvDcEdge;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Real PGrad = 0._Real;

         for (int Q = 0; Q < NQuad; ++Q) {
            const Real Sigma = quadraturePoint(Q);
            const Real Weight = quadratureWeight(Q);

            const Real Pressure0 =
                (1._Real - Sigma) * PressureInterface(ICell0, K) +
                Sigma * PressureInterface(ICell0, K + 1);
            const Real Pressure1 =
                (1._Real - Sigma) * PressureInterface(ICell1, K) +
                Sigma * PressureInterface(ICell1, K + 1);
            const Real ConservTemp0 =
                reconstruct(ConservTemp, ICell0, K, Sigma);
            const Real ConservTemp1 =
                reconstruct(ConservTemp, ICell1, K, Sigma);
            const Real AbsSalinity0 =
                reconstruct(AbsSalinity, ICell0, K, Sigma);
            const Real AbsSalinity1 =
                reconstruct(AbsSalinity, ICell1, K, Sigma);
            const Real SpecVol0 = reconstruct(SpecVol, ICell0, K, Sigma);
            const Real SpecVol1 = reconstruct(SpecVol, ICell1, K, Sigma);

            Real SpecVolGradPressure = 0._Real;
            for (int R = 0; R < NQuad; ++R) {
               const Real Lambda = quadraturePoint(R);
               const Real EdgeWeight = quadratureWeight(R);
               const Real Pressure =
                   (1._Real - Lambda) * Pressure0 + Lambda * Pressure1;

               Real SpecVolAtQuad;
               if (EosChoice == EosType::Teos10Eos) {
                  const Real ConservTempAtQuad =
                      (1._Real - Lambda) * ConservTemp0 +
                      Lambda * ConservTemp1;
                  const Real AbsSalinityAtQuad =
                      (1._Real - Lambda) * AbsSalinity0 +
                      Lambda * AbsSalinity1;
                  SpecVolAtQuad = Teos10.evaluate(
                      ConservTempAtQuad, AbsSalinityAtQuad, Pressure);
               } else {
                  SpecVolAtQuad =
                      (1._Real - Lambda) * SpecVol0 + Lambda * SpecVol1;
               }

               SpecVolGradPressure +=
                   EdgeWeight * SpecVolAtQuad *
                   (Pressure1 - Pressure0) * InvDcEdge;
            }

            const Real GradGeomZ =
                (GeomZQuadrature(ICell1, K, Q) -
                 GeomZQuadrature(ICell0, K, Q)) *
                InvDcEdge;

            PGrad -=
                Weight * (SpecVolGradPressure + Gravity * GradGeomZ);
         }

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * (PGrad - GradGeoPot);
      }
   }

 private:
   static constexpr I4 NQuad = 3;

   KOKKOS_FUNCTION static Real quadraturePoint(const I4 Q) {
      constexpr Real Offset = 0.38729833462074168852_Real;
      return Q == 0 ? 0.5_Real - Offset
                    : (Q == 1 ? 0.5_Real : 0.5_Real + Offset);
   }

   KOKKOS_FUNCTION static Real quadratureWeight(const I4 Q) {
      return Q == 1 ? 4._Real / 9._Real : 5._Real / 18._Real;
   }

   KOKKOS_FUNCTION Real limitedSlope(const Array2DReal &Field, const I4 ICell,
                                     const I4 K) const {
      if (K <= MinLayerCell(ICell) || K >= MaxLayerCell(ICell)) {
         return 0._Real;
      }

      const Real DeltaTop = Field(ICell, K) - Field(ICell, K - 1);
      const Real DeltaBot = Field(ICell, K + 1) - Field(ICell, K);
      if (DeltaTop * DeltaBot <= 0._Real) {
         return 0._Real;
      }

      const Real Centered = 0.5_Real * (DeltaTop + DeltaBot);
      const Real Sign = Centered < 0._Real ? -1._Real : 1._Real;
      return Sign *
             Kokkos::min(Kokkos::abs(Centered),
                         Kokkos::min(2._Real * Kokkos::abs(DeltaTop),
                                     2._Real * Kokkos::abs(DeltaBot)));
   }

   KOKKOS_FUNCTION Real reconstruct(const Array2DReal &Field, const I4 ICell,
                                    const I4 K, const Real Sigma) const {
      return Field(ICell, K) +
             (Sigma - 0.5_Real) * limitedSlope(Field, ICell, K);
   }

   KOKKOS_FUNCTION static Real lagrangeAntiderivative(
       const Real Sigma, const Real XOther0, const Real XOther1) {
      return Sigma * Sigma * Sigma / 3._Real -
             0.5_Real * (XOther0 + XOther1) * Sigma * Sigma +
             XOther0 * XOther1 * Sigma;
   }

   KOKKOS_FUNCTION static Real integratedBasis(const I4 I,
                                               const Real Sigma) {
      const I4 J = I == 0 ? 1 : 0;
      const I4 L = I == 2 ? 1 : 2;
      const Real XI = quadraturePoint(I);
      const Real XJ = quadraturePoint(J);
      const Real XL = quadraturePoint(L);
      const Real Denom = (XI - XJ) * (XI - XL);
      return (lagrangeAntiderivative(1._Real, XJ, XL) -
              lagrangeAntiderivative(Sigma, XJ, XL)) /
             Denom;
   }

   KOKKOS_FUNCTION Real evaluateSpecVol(
       const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
       const Array2DReal &SpecVol, const I4 ICell, const I4 K,
       const Real Sigma, const Real Pressure, EosType EosChoice) const {
      if (EosChoice == EosType::Teos10Eos) {
         return Teos10.evaluate(reconstruct(ConservTemp, ICell, K, Sigma),
                                reconstruct(AbsSalinity, ICell, K, Sigma),
                                Pressure);
      }
      return reconstruct(SpecVol, ICell, K, Sigma);
   }

   I4 NCellsAll;
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DReal BottomGeomDepth;
   Array3DReal GeomZQuadrature;
   Teos10Eos Teos10;
};

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
                            const Array2DReal &GeomZMid,
                            const Array2DReal &ConservTemp,
                            const Array2DReal &AbsSalinity,
                            EosType EosChoice) const;

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
   PressureGradHighOrder HighOrderPGrad;

   // Choice from config
   PressureGradType PressureGradChoice = PressureGradType::Centered;

   // Map of all pressure gradient objects by name
   static std::map<std::string, std::unique_ptr<PressureGrad>> AllPGrads;

}; // end class PressureGrad

} // namespace OMEGA
#endif

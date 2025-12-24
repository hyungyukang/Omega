#ifndef OMEGA_VERTMIX_H
#define OMEGA_VERTMIX_H
//===-- ocn/VertMix.cpp - Vertical Mixing Coefficients -----------*- C++
//-*-===//
//
/// \file
/// \brief Contains functors for calculating vertical diffusivity and viscosity
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the vertical diffusivity and viscosity
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

class ConvectiveMix {
 public:
   bool Enabled = true; ///< Enable convective mixing flag

   // Convective mixing parameters
   Real ConvDiff =
       1.0; ///< Convective vertical viscosity and diffusivity (m^2 s^-1)
   Real ConvTriggerBVF = 0.0; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// Constructor for ConvectiveMix
   ConvectiveMix(const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(Array2DReal VertDiff, Array2DReal VertVisc,
                                   I4 ICell, I4 KChunk,
                                   const Array2DReal &BruntVaisalaFreq) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         if (K == MinLayerCell(ICell)) {
            VertVisc(ICell, MinLayerCell(ICell)) = 0.0_Real;
            VertDiff(ICell, MinLayerCell(ICell)) = 0.0_Real;
         } else {
            if (BruntVaisalaFreq(ICell, K) < ConvTriggerBVF) {
               VertDiff(ICell, K) += ConvDiff;
               VertVisc(ICell, K) += ConvDiff;
            }
         }
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

class ShearMix {
 public:
   bool Enabled = true; ///< Enable shear mixing flag

   // Shear mixing parameters
   Real ShearNuZero =
       0.005; ///< Numerator of Pacanowski and Philander (1981) Eq (1).
   Real ShearAlpha = 5.0; ///< Alpha value used in Pacanowski and Philander
                          ///< (1981) Eqs (1) and (2).
   Real ShearExponent =
       2.0; /// Exponent value used in Pacanowski and Philander (1981) Eqs (1).

   /// Constructor for ShearMix
   ShearMix(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(Array2DReal VertDiff, Array2DReal VertVisc,
                                   I4 ICell, I4 KChunk,
                                   const Array2DReal &NormalVelocity,
                                   const Array2DReal &TangentialVelocity,
                                   const Array2DReal &BruntVaisalaFreq) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;

      Real ShearSquared[VecLength] = {0};
      const Real InvAreaCell       = 1.0_Real / AreaCell(ICell);

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         const Real Factor =
             0.5_Real * DcEdge(JEdge) * DvEdge(JEdge) * InvAreaCell;

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;

            const Real DelNormVel =
                NormalVelocity(JEdge, K - 1) - NormalVelocity(JEdge, K);
            const Real DelTangVel =
                TangentialVelocity(JEdge, K - 1) - TangentialVelocity(JEdge, K);
            ShearSquared[KVec] +=
                Factor * (DelNormVel * DelNormVel + DelTangVel * DelTangVel);
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;

         if (K == MinLayerCell(ICell)) {
            VertVisc(ICell, MinLayerCell(ICell)) = 0.0_Real;
            VertDiff(ICell, MinLayerCell(ICell)) = 0.0_Real;
         } else {
            const Real DelZMid = ZMid(ICell, K - 1) - ZMid(ICell, K);
            ShearSquared[KVec] /= (DelZMid * DelZMid);
            const Real RichardsonNum =
                BruntVaisalaFreq(ICell, K) /
                Kokkos::max(1.0e-12_Real, ShearSquared[KVec]);
            const Real ShearViscVal =
                ShearNuZero / Kokkos::pow(1.0_Real + ShearAlpha * RichardsonNum,
                                          ShearExponent);
            VertVisc(ICell, K) += ShearViscVal;
            VertDiff(ICell, K) +=
                ShearViscVal / (1.0_Real + ShearAlpha * RichardsonNum);
         }
      }
   }

 private:
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal ZMid;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Class for Vertical Mixing Coefficient (VertMix) calculations
class VertMix {
 public:
   /// Get instance of VertMix
   static VertMix *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   Array2DReal VertDiff; ///< Vertical diffusivity field (m^2 s^-1)
   Array2DReal VertVisc; ///< Vertical viscosity field (m^2 s^-1)

   std::string VertDiffFldName;  ///< Field name for vertical diffusivity
   std::string VertViscFldName;  ///< Field name for vertical viscosity
   std::string VertMixGroupName; ///< VertMix group name (for config)
   std::string Name;             ///< Name of this VertMix instance

   // TODO: Temporary array allocations for tridiagonal solver
   //       This array allocation makes debugging easier, but requires more
   //       memory and possibly degrades performance.
   //       Once debugging is done, tridiagonal solver with Kokkos Scratch
   //       should be implemented.
   Array2DReal GWorkEdge;
   Array2DReal HWorkEdge;
   Array2DReal XWorkEdge;
   Array2DReal GWorkCell;
   Array2DReal HWorkCell;
   Array2DReal XWorkCell;

   // Background mixing parameters
   Real BackDiff = 1.0e-5; ///< Background vertical diffusivity (m^2 s^-1)
   Real BackVisc = 1.0e-4; ///< Background vertical viscosity (m^2 s^-1)

   ConvectiveMix
       ComputeVertMixConv;       ///< Functor for Convective VertMix calculation
   ShearMix ComputeVertMixShear; ///< Functor for Shear VertMix calculation

   /// Compute vertical diffusivity and viscosity for all cells/layers
   void computeVertMix(const Array2DReal &NormalVelocity,
                       const Array2DReal &TangentialVelocity,
                       const Array2DReal &BruntVaisalaFreq);

   /// Initialize VertMix from config and mesh
   static void init();

 private:
   /// Private constructor
   VertMix(const std::string &Name, const HorzMesh *Mesh,
           const VertCoord *VCoord);

   /// Private destructor
   ~VertMix();

   static VertMix *Instance; ///< Instance pointer
   const HorzMesh *Mesh;
   const VertCoord *VCoord;

   // Delete copy and move constructors and assignment operators
   VertMix(const VertMix &)            = delete;
   VertMix &operator=(const VertMix &) = delete;
   VertMix(VertMix &&)                 = delete;
   VertMix &operator=(VertMix &&)      = delete;

   I4 NCellsAll;   ///< Number of horizontal cells
   I4 NChunks;     ///< Number of vertical chunks (for vectorization)
   I4 NVertLayers; ///< Number of vertical layers

   // Define fields and metadata
   void defineFields();

}; // End class VertMix

} // namespace OMEGA
#endif

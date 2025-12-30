#ifndef OMEGA_TENDENCYTERMS_H
#define OMEGA_TENDENCYTERMS_H
//===-- ocn/TendencyTerms.h - Tendency Terms --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating tendency terms
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate tendencies used to update state variables.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "VertCoord.h"

#include <cmath> // for std::copysign

namespace OMEGA {

/// Divergence of thickness flux at cell centers, for updating layer thickness
/// arrays
class ThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   /// constructor declaration
   ThicknessFluxDivOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes cell index, vertical chunk index, and thickness flux
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell, I4 KChunk,
                                   const Array2DReal &ThicknessFlux,
                                   const Array2DReal &NormalVelEdge) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            DivTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                            ThicknessFlux(JEdge, K) * NormalVelEdge(JEdge, K) *
                            InvAreaCell;
         }
      }

      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(ICell, K) -= DivTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Horizontal advection of potential vorticity defined on edges, for
/// momentum equation
class PotentialVortHAdvOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   PotentialVortHAdvOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// normalized relative vorticity, normalized planetary vorticity, layer
   /// thickness on edges, and normal velocity on edges as inputs,
   /// outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &NormFEdge,
                                   const Array2DReal &FluxLayerThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      Real VortTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const I4 K    = KStart + KVec;
            Real NormVort = (NormRVortEdge(IEdge, K) + NormFEdge(IEdge, K) +
                             NormRVortEdge(JEdge, K) + NormFEdge(JEdge, K)) *
                            0.5_Real;

            VortTmp[KVec] += WeightsOnEdge(IEdge, J) *
                             FluxLayerThickEdge(JEdge, K) *
                             NormVelEdge(JEdge, K) * NormVort;
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += EdgeMask(IEdge, K) * VortTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of kinetic energy defined on edges, for momentum equation
class KEGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and kinetic energy
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &KECell) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -= EdgeMask(IEdge, K) *
                           (KECell(JCell1, K) - KECell(JCell0, K)) * InvDcEdge;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of sea surface height defined on edges multipled by gravitational
/// acceleration, for momentum equation
class SSHGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and array of
   /// layer thickness/SSH, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &SshCell) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -= EdgeMask(IEdge, K) * Gravity *
                           (SshCell(ICell1, K) - SshCell(ICell0, K)) *
                           InvDcEdge;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Laplacian horizontal mixing, for momentum equation
class VelocityDiffusionOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel2;

   /// constructor declaration
   VelocityDiffusionOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// divergence of horizontal velocity (defined at cell centers) and relative
   /// vorticity (defined at vertices), outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &DivCell,
                                   const Array2DReal &RVortVertex) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             ((DivCell(ICell1, K) - DivCell(ICell0, K)) * DcEdgeInv -
              (RVortVertex(IVertex1, K) - RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * ViscDel2 * MeshScalingDel2(IEdge) * Del2U;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel2;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Biharmonic horizontal mixing, for momentum equation
class VelocityHyperDiffOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel4;
   Real DivFactor;

   /// Constructor declaration
   VelocityHyperDiffOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// the laplacian of divergence of horizontal velocity and the laplacian of
   /// the relative vorticity, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &Del2DivCell,
                                   const Array2DReal &Del2RVortVertex) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             (DivFactor * (Del2DivCell(ICell1, K) - Del2DivCell(ICell0, K)) *
                  DcEdgeInv -
              (Del2RVortVertex(IVertex1, K) - Del2RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) -=
             EdgeMask(IEdge, K) * ViscDel4 * MeshScalingDel4(IEdge) * Del2U;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel4;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Velocity vertical mixing
class VelVertMixSetupOnEdge {
 public:
   bool Enabled;
   Real LocRhoSw;

   /// constructor declaration
   VelVertMixSetupOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// layer specific volume, layer thickness on edge,
   /// interface pressure, and outputs tendency array
   KOKKOS_FUNCTION void
   operator()(I4 IEdge, I4 KChunk, Real DT, const Array2DReal &SpecVol,
              const Array2DReal &LayerThickEdge, const Array2DReal &VertVisc,
              const Array2DReal &NormalVelEdge, const Array2DReal &GWorkEdge,
              const Array2DReal &HWorkEdge,
              const Array2DReal &XWorkEdge) const {

      const I4 NVertLayers1 = NVertLayers - 1;
      const I4 KStart       = chunkStart(KChunk, 0);
      const I4 KLen         = chunkLength(KChunk, KStart, NVertLayers1);

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         if (K < KMin) {
            // For K <= KMin, set unit diagonal
            GWorkEdge(IEdge, K) = 0.0_Real;
            HWorkEdge(IEdge, K) = 1.0_Real;
            XWorkEdge(IEdge, K) = 1.0_Real;
         } else if (K == KMin) {
            // For K <= KMin, set unit diagonal
            GWorkEdge(IEdge, K) = 0.0_Real;
            HWorkEdge(IEdge, K) = 1.0_Real;
            XWorkEdge(IEdge, K) = NormalVelEdge(IEdge, K);
         } else if (K == KMax) {
            const Real LayerThickEdgeTop =
                0.5 * (LayerThickEdge(IEdge, K - 1) + LayerThickEdge(IEdge, K));
            const Real SpecVolEdgeTop =
                0.5 * (0.5 * (SpecVol(ICell0, K - 1) + SpecVol(ICell1, K - 1)) +
                       0.5 * (SpecVol(ICell0, K) + SpecVol(ICell1, K)));
            const Real ViscAlphaEdgeTop =
                0.5 * (VertVisc(ICell0, K) + VertVisc(ICell1, K)) /
                (LocRhoSw * SpecVolEdgeTop);
            const Real LayerThickEdgeTopDT = LayerThickEdgeTop / DT;

            GWorkEdge(IEdge, K - 1) =
                ViscAlphaEdgeTop / (LayerThickEdgeTopDT * LayerThickEdgeTop);

            HWorkEdge(IEdge, K) = 1.0_Real;

            XWorkEdge(IEdge, K) = NormalVelEdge(IEdge, K);

         } else if (K > KMax) {
            // For K >= KMax, set unit diagonal
            GWorkEdge(IEdge, K) = 0.0_Real;
            HWorkEdge(IEdge, K) = 1.0_Real;
            XWorkEdge(IEdge, K) = 1.0_Real;
         } else {

            const Real LayerThickEdgeTop =
                0.5 * (LayerThickEdge(IEdge, K - 1) + LayerThickEdge(IEdge, K));
            const Real SpecVolEdgeTop =
                0.5 * (0.5 * (SpecVol(ICell0, K - 1) + SpecVol(ICell1, K - 1)) +
                       0.5 * (SpecVol(ICell0, K) + SpecVol(ICell1, K)));
            const Real ViscAlphaEdgeTop =
                0.5 * (VertVisc(ICell0, K) + VertVisc(ICell1, K)) /
                (LocRhoSw * SpecVolEdgeTop);

            const Real LayerThickEdgeTopDT = LayerThickEdgeTop / DT;

            GWorkEdge(IEdge, K - 1) =
                ViscAlphaEdgeTop / (LayerThickEdgeTopDT * LayerThickEdgeTop);

            HWorkEdge(IEdge, K) = 1.0_Real;

            XWorkEdge(IEdge, K) = NormalVelEdge(IEdge, K);
         }
      }
   }

 private:
   I4 NVertLayers;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Pressure gradient contribution from sloping surfaces
class PresGradZOnEdge {
 public:
   bool Enabled;
   Real LocRhoSw;

   /// constructor declaration
   PresGradZOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// layer specific volume, layer thickness on edge,
   /// interface pressure, and outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &LayerThickEdge,
                                   const Array2DReal &LayerThickCell,
                                   const Array2DReal &PressureInterface) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      const Real InvRhoSwGravity = 1._Real / (LocRhoSw * Gravity);

      // This tendency term has 0 value at Floor and Surface by the boundary
      // condition.

      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         const I4 KM1Cell0 = Kokkos::max(K - 1, MinLayerCell(ICell0));
         const I4 KM1Cell1 = Kokkos::max(K - 1, MinLayerCell(ICell1));
         const I4 KP1Cell0 = Kokkos::min(K + 1, MaxLayerCell(ICell0));
         const I4 KP1Cell1 = Kokkos::min(K + 1, MaxLayerCell(ICell1));

         // LayerThickCell0 at Top (K)
         const Real LayerThickCell0TopK =
             0.5_Real *
             (LayerThickCell(ICell0, KM1Cell0) + LayerThickCell(ICell0, K));

         // LayerThickCell1 at Top (K)
         const Real LayerThickCell1TopK =
             0.5_Real *
             (LayerThickCell(ICell1, KM1Cell1) + LayerThickCell(ICell1, K));

         // LayerThickCell0 at Top (K+1)
         const Real LayerThickCell0TopKP1 =
             0.5_Real *
             (LayerThickCell(ICell0, K) + LayerThickCell(ICell0, KP1Cell0));

         // LayerThickCell1 at Top (K+1)
         const Real LayerThickCell1TopKP1 =
             0.5_Real *
             (LayerThickCell(ICell1, K) + LayerThickCell(ICell1, KP1Cell1));

         // SpecVolCell0 at Top (K)
         const Real SpecVolCell0TopK =
             0.5_Real *
             (SpecVol(ICell0, KM1Cell0) * LayerThickCell(ICell0, KM1Cell0) +
              SpecVol(ICell0, K) * LayerThickCell(ICell0, K)) /
             LayerThickCell0TopK;

         // SpecVolCell1 at Top (K)
         const Real SpecVolCell1TopK =
             0.5_Real *
             (SpecVol(ICell1, KM1Cell1) * LayerThickCell(ICell1, KM1Cell1) +
              SpecVol(ICell1, K) * LayerThickCell(ICell1, K)) /
             LayerThickCell1TopK;

         // SpecVolCell0 at Top (K+1)
         const Real SpecVolCell0TopKP1 =
             0.5_Real *
             (SpecVol(ICell0, K) * LayerThickCell(ICell0, K) +
              SpecVol(ICell0, KP1Cell0) * LayerThickCell(ICell0, KP1Cell0)) /
             LayerThickCell0TopK;

         // SpecVolCell1 at Top (K+1)
         const Real SpecVolCell1TopKP1 =
             0.5_Real *
             (SpecVol(ICell1, K) * LayerThickCell(ICell1, K) +
              SpecVol(ICell1, KP1Cell1) * LayerThickCell(ICell1, KP1Cell1)) /
             LayerThickCell1TopK;

         // Pressue * SpecVol at Top (K)
         const Real PSpecVolEdgeTopK =
             0.5_Real * (PressureInterface(ICell0, K) * SpecVolCell0TopK +
                         PressureInterface(ICell1, K) * SpecVolCell1TopK);

         // Pressue * SpecVol at Top (K+1)
         const Real PSpecVolEdgeTopKP1 =
             0.5_Real *
             (PressureInterface(ICell0, KP1Cell0) * SpecVolCell0TopKP1 +
              PressureInterface(ICell1, KP1Cell1) * SpecVolCell1TopKP1);

         // Compute grad(\tilde{z}) = grad(-p) / (Rho0 * Gravity)
         const Real GradZTildeTopK =
             (-PressureInterface(ICell1, K) + PressureInterface(ICell0, K)) *
             InvRhoSwGravity;
         const Real GradZTildeTopKP1 = (-PressureInterface(ICell1, KP1Cell1) +
                                        PressureInterface(ICell0, KP1Cell0)) *
                                       InvRhoSwGravity;

         const Real InvLayerThickEdge = 1._Real / LayerThickEdge(IEdge, K);

         Real PresSpecVolGradZK   = PSpecVolEdgeTopK * GradZTildeTopK;
         Real PresSpecVolGradZKP1 = PSpecVolEdgeTopKP1 * GradZTildeTopKP1;

         // 0 at surface
         if (K == KMin)
            PresSpecVolGradZK = 0._Real;
         // 0 at floor
         if (K == KMax)
            PresSpecVolGradZKP1 = 0._Real;

         const Real ZGradTerm = -InvDcEdge * InvLayerThickEdge *
                                (PresSpecVolGradZK - PresSpecVolGradZKP1);

         Tend(IEdge, K) += EdgeMask(IEdge, K) * ZGradTerm;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Pressure gradient force
class PresGradForceOnEdge {
 public:
   bool Enabled;

   /// constructor declaration
   PresGradForceOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// layer specific volume, layer thickness on edge, layer thickness on cell,
   /// layer-mean pressure, and outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &LayerThickEdge,
                                   const Array2DReal &LayerThickCell,
                                   const Array2DReal &PressureMid) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         // LayerThick * SpecVol * Pressure : (A)
         const Real LayerSpecVolPresCell0 = LayerThickCell(ICell0, K) *
                                            SpecVol(ICell0, K) *
                                            PressureMid(ICell0, K);

         const Real LayerSpecVolPresCell1 = LayerThickCell(ICell1, K) *
                                            SpecVol(ICell1, K) *
                                            PressureMid(ICell1, K);

         // -grad(A) / LayerThick
         const Real PGFTerm = -InvDcEdge *
                              (LayerSpecVolPresCell1 - LayerSpecVolPresCell0) /
                              LayerThickEdge(IEdge, K);

         Tend(IEdge, K) += EdgeMask(IEdge, K) * PGFTerm;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of geopotential
class GeoptGradOnEdge {
 public:
   bool Enabled;

   /// constructor declaration
   GeoptGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// geopotential, and outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &GeoptMid) const {

      const I4 KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         // -grad(Geopotential)
         const Real GeoptGradTerm =
             -InvDcEdge * (GeoptMid(ICell1, K) - GeoptMid(ICell0, K));

         Tend(IEdge, K) += EdgeMask(IEdge, K) * GeoptGradTerm;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Wind forcing
class WindForcingOnEdge {
 public:
   bool Enabled = false;
   Real LocRhoSw;

   /// constructor declaration
   WindForcingOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// normal wind stress and edge layer thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array1DReal &NormalStressEdge,
                                   const Array2DReal &LayerThickEdge) const {
      if (KChunk == 0) {
         const I4 K = MinLayerEdgeBot(IEdge);

         const Real InvThickEdge = 1._Real / LayerThickEdge(IEdge, K);
         Tend(IEdge, K) += EdgeMask(IEdge, K) * InvThickEdge *
                           NormalStressEdge(IEdge) / LocRhoSw;
      }
   }

 private:
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
};

/// Bottom drag
class BottomDragOnEdge {
 public:
   bool Enabled = false;
   Real Coeff;

   /// constructor declaration
   BottomDragOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index and arrays for
   /// horizontal velocity, kinetic energy,
   /// and edge layer thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormalVelEdge,
                                   const Array2DReal &KECell,
                                   const Array2DReal &LayerThickEdge) const {
      const I4 KBot = MaxLayerEdgeTop(IEdge);

      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);

      const Real VelNormEdge =
          Kokkos::sqrt(KECell(JCell0, KBot) + KECell(JCell1, KBot));

      const Real InvThickEdge = 1._Real / LayerThickEdge(IEdge, KBot);
      Tend(IEdge, KBot) -= EdgeMask(IEdge, KBot) * Coeff * VelNormEdge *
                           InvThickEdge * NormalVelEdge(IEdge, KBot);
   }

 private:
   I4 NVertLayers;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeMask;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer horizontal advection term
class TracerHorzAdvOnCell {
 public:
   bool Enabled = false;

   TracerHorzAdvOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, const I4 L,
                                   const I4 ICell, const I4 KChunk,
                                   const Array2DReal &NormVelEdge,
                                   const Array3DReal &HTracersOnEdge) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell  = 1._Real / AreaCell(ICell);
      Real HAdvTmp[VecLength] = {0};
      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            HAdvTmp[KVec] -= EdgeMask(JEdge, K) * DvEdge(JEdge) *
                             EdgeSignOnCell(ICell, J) *
                             HTracersOnEdge(L, JEdge, K) *
                             NormVelEdge(JEdge, K) * InvAreaCell;
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(L, ICell, K) -= HAdvTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer high order horizontal advection term
class TracerHighOrderHorzAdvOnCell {
 public:
   bool Enabled = false;
   // coefficient for blending high-order terms
   Real coef3rdOrder = 0.25;

   TracerHighOrderHorzAdvOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);
   void init();

   KOKKOS_FUNCTION void operator()(const I4 L, const I4 IEdge, const I4 KChunk,
                                   const Array3DReal &TracerCell,
                                   const Array2DReal &FluxLayerThickEdge,
                                   const Array2DReal &NormVelEdge) const {
      const I4 KStart = KChunk * VecLength;
      const I4 KEnd   = KStart + VecLength;
      for (int K = KStart; K < KEnd; ++K)
         HighOrderFlxHorz(L, IEdge, K) = 0;
      if (AdvMaskHighOrder(IEdge)) {
         for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
            const I4 ICell = AdvCellsForEdge(IEdge, I);
            for (int K = KStart; K < KEnd; ++K) {
               const Real normalThicknessFlux =
                   FluxLayerThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
               const Real tracerWgt =
                   (AdvCoefs(I, IEdge) +
                    coef3rdOrder * std::copysign(1._Real, normalThicknessFlux) *
                        AdvCoefs3rd(I, IEdge)) *
                   normalThicknessFlux;
               HighOrderFlxHorz(L, IEdge, K) +=
                   tracerWgt * TracerCell(L, ICell, K);
            }
         }
      } else {
         for (int K = KStart; K < KEnd; ++K) {
            const I4 JCell0 = CellsOnEdge(IEdge, 0);
            const I4 JCell1 = CellsOnEdge(IEdge, 1);
            const Real normalThicknessFlux =
                FluxLayerThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
            const Real tracerWgt =
                DvEdge(IEdge) * 0.5_Real * normalThicknessFlux;
            HighOrderFlxHorz(L, IEdge, K) +=
                tracerWgt *
                (TracerCell(L, JCell1, K) + TracerCell(L, JCell0, K));
         }
      }
   }

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, const I4 L,
                                   const I4 ICell, const I4 KChunk) const {
      const I4 KStart        = KChunk * VecLength;
      const I4 KEnd          = KStart + VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      for (int I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);
         for (int K = KStart; K < KEnd; ++K) {
            Tend(L, ICell, K) += EdgeSignOnCell(ICell, I) *
                                 HighOrderFlxHorz(L, IEdge, K) * InvAreaCell;
         }
      }
   }

 private:
   const HorzMesh *HorzontalMesh;
   Array1DI4 NAdvCellsForEdge;
   Array2DI4 AdvCellsForEdge;
   Array1DI4 AdvMaskHighOrder;
   Array2DReal AdvCoefs;
   Array2DReal AdvCoefs3rd;
   Array3DReal HighOrderFlxHorz;

   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
};

// Tracer horizontal diffusion term
class TracerDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff2;

   TracerDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 L, I4 ICell, I4 KChunk,
              const Array3DReal &TracerCell,
              const Array2DReal &MeanLayerThickEdge) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DiffTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel2(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            const Real TracerGrad =
                (TracerCell(L, JCell1, K) - TracerCell(L, JCell0, K));

            DiffTmp[KVec] -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                             RTemp * MeanLayerThickEdge(JEdge, K) * TracerGrad;
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(L, ICell, K) += EddyDiff2 * DiffTmp[KVec] * InvAreaCell;
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel2;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer biharmonic horizontal mixing term
class TracerHyperDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff4;

   TracerHyperDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk,
                                   const Array3DReal &TrDel2Cell) const {

      const I4 KStartCell = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLenCell = chunkLength(KChunk, KStartCell, MaxLayerCell(ICell));
      const I4 KEndCell = KStartCell + KLenCell - 1;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HypTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge      = EdgesOnCell(ICell, J);
         const I4 KStartEdge = Kokkos::max(KStartCell, MinLayerEdgeBot(JEdge));
         const I4 KEndEdge   = Kokkos::min(KEndCell, MaxLayerEdgeTop(JEdge));

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel4(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const I4 KVec = K - KStartCell;
            const Real Del2TrGrad =
                (TrDel2Cell(L, JCell1, K) - TrDel2Cell(L, JCell0, K));

            HypTmp[KVec] -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                            RTemp * Del2TrGrad;
         }
      }
      for (int KVec = 0; KVec < KLenCell; ++KVec) {
         const I4 K = KStartCell + KVec;
         Tend(L, ICell, K) -= EddyDiff4 * HypTmp[KVec] * InvAreaCell;
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel4;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

// Tracer vertical mixing term
class TracerVertMixSetupOnCell {
 public:
   bool Enabled;
   Real LocRhoSw;

   TracerVertMixSetupOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(I4 L, I4 ICell, I4 KChunk, R8 DT, const Array2DReal &SpecVol,
              const Array2DReal &LayerThickCell, const Array2DReal &VertDiff,
              const Array3DReal &TracersOnCell, const Array2DReal &GWorkCell,
              const Array2DReal &HWorkCell,
              const Array2DReal &XWorkCell) const {

      const I4 NVertLayers1 = NVertLayers - 1;
      const I4 KStart       = chunkStart(KChunk, 0);
      const I4 KLen         = chunkLength(KChunk, KStart, NVertLayers1);

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         if (K < KMin) {
            // For K < KMin, set unit diagonal
            GWorkCell(ICell, K) = 0.0_Real;
            HWorkCell(ICell, K) = 1.0_Real;
            XWorkCell(ICell, K) = 1.0_Real;
         } else if (K == KMin) {
            GWorkCell(ICell, K) = 0.0_Real;
            HWorkCell(ICell, K) = 1.0_Real;
            XWorkCell(ICell, K) = TracersOnCell(L, ICell, K);
         } else if (K == KMax) {
            const Real LayerThickCellTop =
                0.5 * (LayerThickCell(ICell, K - 1) + LayerThickCell(ICell, K));
            const Real SpecVolCellTop =
                0.5 * (SpecVol(ICell, K - 1) + SpecVol(ICell, K));
            const Real DiffAlphaCellTop =
                VertDiff(ICell, K) / (LocRhoSw * SpecVolCellTop);

            GWorkCell(ICell, K-1) =
                DT * DiffAlphaCellTop /
                (LayerThickCellTop * LayerThickCell(ICell,K));

            HWorkCell(ICell, K) = 1.0_Real;

            XWorkCell(ICell, K) = TracersOnCell(L, ICell, K);

         } else if (K > KMax) {
            // For K > KMax, set unit diagonal
            GWorkCell(ICell, K) = 0.0_Real;
            HWorkCell(ICell, K) = 1.0_Real;
            XWorkCell(ICell, K) = 1.0_Real;
         } else {

            const Real LayerThickCellTop =
                0.5 * (LayerThickCell(ICell, K-1 ) + LayerThickCell(ICell, K));
            const Real SpecVolCellTop =
                0.5 * (SpecVol(ICell, K-1 ) + SpecVol(ICell, K));
            const Real DiffAlphaCellTop =
                VertDiff(ICell, K) / (LocRhoSw * SpecVolCellTop);

            GWorkCell(ICell, K-1) =
                DT * DiffAlphaCellTop /
                (LayerThickCellTop * LayerThickCell(ICell,K));

            HWorkCell(ICell, K) = 1.0_Real;

            XWorkCell(ICell, K) = TracersOnCell(L, ICell, K);
         }
      }
   }

 private:
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif

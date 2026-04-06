#ifndef OMEGA_AUX_KINETIC_H
#define OMEGA_AUX_KINETIC_H

#include "DataTypes.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class KineticAuxVars {
 public:
   Real LocRhoSw;

   Array2DReal KineticEnergyCell;
   Array2DReal VelocityDivCell;
   Array2DReal ProjVelDivCell;

   KineticAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                  const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int KChunk,
                     const Array2DReal &NormalVelEdge,
                     const Array2DReal &LayerThickCell,
                     const Array2DReal &PressureInterface) const {

      const int KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real KineticEnergyCellTmp[VecLength] = {0};
      Real VelocityDivCellTmp[VecLength]   = {0};
      Real ProjVelDivCellTmp[VecLength]    = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge     = EdgesOnCell(ICell, J);
         const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            KineticEnergyCellTmp[KVec] += AreaEdge * 0.5_Real * InvAreaCell *
                                          NormalVelEdge(JEdge, K) *
                                          NormalVelEdge(JEdge, K);
            VelocityDivCellTmp[KVec] -= DvEdge(JEdge) * InvAreaCell *
                                        EdgeSignOnCell(ICell, J) *
                                        NormalVelEdge(JEdge, K);

            const Real GradZTildeEdgeTop = (PressureInterface(JCell1,K) -
                                           PressureInterface(JCell0,K) ) /
                                          (Gravity * LocRhoSw * DcEdge(JEdge));
            const Real LayerThickEdgeKM1
                = 0.5 * (LayerThickCell(JCell0,K-1) +
                         LayerThickCell(JCell1,K-1));
            const Real LayerThickEdge
                = 0.5 * (LayerThickCell(JCell0,K) +
                         LayerThickCell(JCell1,K));
            const Real NormalVelEdgeTop
                = 0.5 * (NormalVelEdge(JEdge,K-1) * LayerThickEdgeKM1 +
                         NormalVelEdge(JEdge,K) * LayerThickEdge) /
                  (0.5 * (LayerThickEdgeKM1 + LayerThickEdge));
            ProjVelDivCellTmp[KVec] -= DvEdge(JEdge) * InvAreaCell *
                                       EdgeSignOnCell(ICell, J) *
                                       NormalVelEdgeTop *
                                       GradZTildeEdgeTop;
         }
      }
      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K                 = KStart + KVec;
         KineticEnergyCell(ICell, K) = KineticEnergyCellTmp[KVec];
         VelocityDivCell(ICell, K)   = VelocityDivCellTmp[KVec];

         if (K == MinLayerCell(ICell)) {
            ProjVelDivCell(ICell, K)   = 0._Real;
         } else {
            ProjVelDivCell(ICell, K)   = ProjVelDivCellTmp[KVec];
         }
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif

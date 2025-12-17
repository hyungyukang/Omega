#ifndef OMEGA_AUX_TANGENT_H
#define OMEGA_AUX_TANGENT_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class TangentAuxVars {
 public:
   Array2DReal TangentVelEdge;

   TangentAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                  const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = chunkStart(KChunk, MinLayerEdgeTop(IEdge));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeBot(IEdge));

      // Compute tangential velocity
      Real ReconEdgeTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const int JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            ReconEdgeTmp[KVec] +=
                WeightsOnEdge(IEdge, J) * NormalVelEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K              = KStart + KVec;
         TangentVelEdge(IEdge, K) = ReconEdgeTmp[KVec];
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeBot;
   Array2DReal WeightsOnEdge;
};

} // namespace OMEGA
#endif

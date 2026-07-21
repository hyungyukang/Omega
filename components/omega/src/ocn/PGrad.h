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
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &GeomZMid,
                                   const Array1DReal &TidalPotential,
                                   const Array1DReal &SelfAttractionLoading) const {

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
                            const Array2DReal &SpecVol,
                            const Array2DReal &GeomZMid) const;

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

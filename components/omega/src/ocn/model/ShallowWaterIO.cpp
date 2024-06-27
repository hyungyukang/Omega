//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"

#include <iostream>
#include <cmath>
#include <iomanip>

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_init_io(const Decomp *DefDecomp, const HorzMesh *Mesh, const OceanState *State) {

   //--------------------------------------------------------------------------/

   // IO initialize -----------------------------------------------------------/
   //int OutFileID;
   OutFileID = 0;
   int Err0 = IO::openFile(
         OutFileID, "SWoutput.nc", IO::ModeWrite, IO::FmtDefault,
         IO::IfExists::Replace, IO::Precision::Double);
   if ( Err0 != 0)
      LOG_ERROR("IO: error opening file for output FAILE");

   //--------------------------------------------------------------------------/

//   // Define IO field
//   auto EdgeDim          = MetaDim::create("NEdges", Mesh->NEdgesSize);
//   auto CellDim          = MetaDim::create("NCells", Mesh->NCellsSize);
//   auto VertDim          = MetaDim::create("NVertLevels", NVertLevels);
//   std::vector<std::shared_ptr<OMEGA::MetaDim>> CellDim2DR8{CellDim, VertDim};
//   std::vector<std::shared_ptr<OMEGA::MetaDim>> EdgeDim2DR8{EdgeDim, VertDim};
//   auto FieldNormalVelocity = ArrayMetaData::create(
//       "NormalVelocity",
//       "Normal Velocity on Edge", /// long Name or description
//       "m/s",                         /// units
//       "R8StdName",                 /// CF standard Name
//       -12345.0,                    /// min valid value
//       12345.0,                     /// max valid value
//       -9.99E+30,                   /// scalar used for undefined entries
//       2,                           /// number of dimensions
//       EdgeDim2DR8                  /// dim pointers
//   );
//   int Err1 = IOField::define("NormalVelocity");
//   int Err2 = IOField::attachData<OMEGA::Array2DR8>("NormalVelocity",State->NormalVelocity[0]);
//
//
//   auto FieldLayerThickness = ArrayMetaData::create(
//       "LayerThickness",
//       "Layer thickness on Cell", /// long Name or description
//       "m",                           /// units
//       "R8StdName",                 /// CF standard Name
//       0,                           /// min valid value
//       12345.0,                     /// max valid value
//       -9.99E+30,                   /// scalar used for undefined entries
//       2,                           /// number of dimensions
//       CellDim2DR8                  /// dim pointers
//   );
//   int Err3 = IOField::define("LayerThickness");
//   int Err4 = IOField::attachData<OMEGA::Array2DR8>("LayerThickness",State->LayerThickness[0]);
//

   //--------------------------------------------------------------------------/

   //if ( Err1 != 0 || Err2 != 0 || Err3 != 0 || Err4 != 0) {
   //   LOG_CRITICAL("IOField: Initializing State IO field FAIL");
   //}
   //
   int Err;

   //--------------------------------------------------------------------------/
   std::vector<int> OffsetEdge(Mesh->NEdgesSize * NVertLevels, -1);
   for (int Edge = 0; Edge < Mesh->NEdgesOwned; ++Edge) {
      int GlobalEdgeAdd = DefDecomp->EdgeIDH(Edge) - 1;
      for (int k = 0; k < NVertLevels; ++k) {
         int VectorAdd         = Edge * NVertLevels + k;
         OffsetEdge[VectorAdd] = GlobalEdgeAdd * NVertLevels + k;
      }
   }

   std::vector<int> EdgeDims{DefDecomp->NEdgesGlobal, NVertLevels};
   int ErrDecomp = IO::createDecomp(DecompEdgeR8, IO::IOTypeR8, 2,
                                 EdgeDims, Mesh->NEdgesSize * NVertLevels, OffsetEdge,
                                 IO::DefaultRearr);

   int DimVertID;
   int ErrDimVert = OMEGA::IO::defineDim(OutFileID, "NVertLevels", NVertLevels,
                                     DimVertID);

   int DimEdgeID;
   int ErrDim = OMEGA::IO::defineDim(OutFileID, "NEdges", DefDecomp->NEdgesGlobal, DimEdgeID);
   int EdgeDimIDs[2] = {DimEdgeID, DimVertID};

   //int NormalVelocityID;
   int ErrVar = IO::defineVar(OutFileID, "NormalVelocity", IO::IOTypeR8, 2, EdgeDimIDs, NormalVelocityID);

   //R8 VarMetaR8Ref      = 2.23456789;
   //int ErrMeta = IO::writeMeta("NormalVelocity", VarMetaR8Ref, OutFileID,
   //                           NormalVelocityID);

   // Write array
   //R8 FillR8 = -1.23456789e30;
   //int ErrWrite = IO::writeArray(State->NormalVelocity[0].data(), Mesh->NEdgesSize * NVertLevels,
   //                             &FillR8, OutFileID, DecompEdgeR8, NormalVelocityID);
   //--------------------------------------------------------------------------/

   std::vector<int> OffsetCell(Mesh->NCellsSize * NVertLevels, -1);
   for (int Cell = 0; Cell < Mesh->NCellsOwned; ++Cell) {
      int GlobalCellAdd = DefDecomp->CellIDH(Cell) - 1;
      for (int k = 0; k < NVertLevels; ++k) {
         int VectorAdd         = Cell * NVertLevels + k;
         OffsetCell[VectorAdd] = GlobalCellAdd * NVertLevels + k;
      }
   }

   std::vector<int> CellDims{DefDecomp->NCellsGlobal, NVertLevels};
   //int DecompCellR8;
   Err = IO::createDecomp(DecompCellR8, IO::IOTypeR8, 2,
                                 CellDims, Mesh->NCellsSize * NVertLevels, OffsetCell,
                                 IO::DefaultRearr);

   

   int DimCellID;
   Err = OMEGA::IO::defineDim(OutFileID, "NCells", DefDecomp->NCellsGlobal, DimCellID);
   int CellDimIDs[2] = {DimCellID, DimVertID};

   //int LayerThicknessID;
   Err = IO::defineVar(OutFileID, "LayerThickness", IO::IOTypeR8, 2, CellDimIDs, LayerThicknessID);

   //R8 VarMetaR8Ref      = 2.23456789;
   //int ErrMeta = IO::writeMeta("LayerThickness", VarMetaR8Ref, OutFileID,
   //                             LayerThicknessID);

   // Write array
   //Err = IO::writeArray(State->LayerThickness[0].data(), Mesh->NCellsSize * NVertLevels,
   //                            &FillR8, OutFileID, DecompCellR8, LayerThicknessID);

   //--------------------------------------------------------------------------/

} // sw_init_io

void ShallowWaterCore::
sw_io_write(const Decomp *DefDecomp, const HorzMesh *Mesh, const OceanState *State) {

//   TangentialReconOnEdge2D TanReconEdge2D(Mesh);
//   parallelFor(
//      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
//         const R8 tangentialComp = TanReconEdge2D(IEdge,KLevel, State->NormalVelocity[0]);
//         TangentialVelocity(IEdge,KLevel) = tangentialComp;
//      });
//   parallelFor(
//      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
//         const R8 normalComp = State->NormalVelocity[0](IEdge,KLevel);
//         const R8 tangentialComp = TangentialVelocity(IEdge,KLevel);
//
//         const R8 zonalComp = normalComp * cos(Mesh->AngleEdge(IEdge))
//                            - tangentialComp * sin(Mesh->AngleEdge(IEdge));
//
//         //const R8 meridComp = normalComp * cos(Mesh->AngleEdge(IEdge))
//         //                   - tangentialComp * sin(Mesh->AngleEdge(IEdge));
//
//         //const R8 zonalComp = normalComp * cos(Mesh->AngleEdge(IEdge))
//         //                   - tangentialComp * sin(Mesh->AngleEdge(IEdge));
//         //const R8 meridComp = 0.0;
//
//         State->NormalVelocity[0](IEdge,KLevel) = zonalComp;
//         State->NormalVelocity[1](IEdge,KLevel) = 0.0;
//
//      });

   // Save total depth ( H = h + b )
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         State->LayerThicknessH[0](ICell,KLevel) = State->LayerThicknessH[0](ICell,KLevel)+BottomTopography(ICell);
       });   

   int ErrWriteVel   = IO::writeArray(State->NormalVelocity[0].data(), Mesh->NEdgesSize * NVertLevels,
                                      &FillR8, OutFileID, DecompEdgeR8, NormalVelocityID);
   int ErrWriteThick = IO::writeArray(State->LayerThickness[0].data(), Mesh->NCellsSize * NVertLevels,
                                      &FillR8, OutFileID, DecompCellR8, LayerThicknessID);

   // Finished writing, close file
   int ErrClose = IO::closeFile(OutFileID);
   if (ErrClose != 0)
      LOG_ERROR("IO: error closing output file FAIL");

} // sw_io_write

} // namespace OMEGA

//===-----------------------------------------------------------------------===/

#include "ShallowWaterCore.h"
#include "pio.h"

#include <iostream>
#include <cmath>
#include <iomanip>

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWaterCore::
sw_init_io(const Decomp *DefDecomp, HorzMesh *Mesh, OceanState *State) {

   //--------------------------------------------------------------------------/
   int Err = 0;

   // IO initialize -----------------------------------------------------------/
   //int OutFileID;
   OutFileID = 0;
   int Err0 = IO::openFile(
         OutFileID, "output.nc", IO::ModeWrite, IO::FmtDefault,
         IO::IfExists::Replace, IO::Precision::Double);
   if ( Err0 != 0)
      LOG_ERROR("IO: error opening file for output FAILE");

   //--------------------------------------------------------------------------/

   // Define IO field
   auto TimeDim          = MetaDim::create("Time", NTimeLevels);
   auto EdgeDim          = MetaDim::create("NEdges", Mesh->NEdgesSize);
   auto CellDim          = MetaDim::create("NCells", Mesh->NCellsSize);
   auto VertDim          = MetaDim::create("NVertLevels", NVertLevels);
   std::vector<std::shared_ptr<OMEGA::MetaDim>> CellDim2DR8{TimeDim,CellDim};
   std::vector<std::shared_ptr<OMEGA::MetaDim>> CellDim3DR8{TimeDim,CellDim,VertDim};
   std::vector<std::shared_ptr<OMEGA::MetaDim>> EdgeDim3DR8{TimeDim,EdgeDim,VertDim};

   auto FieldNormalVelocity = ArrayMetaData::create(
       "NormalVelocity",
       "Normal Velocity on Edge", /// long Name or description
       "m/s",                         /// units
       "R8StdName",                 /// CF standard Name
       -12345.0,                    /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       3,                           /// number of dimensions
       EdgeDim3DR8                  /// dim pointers
   );
   int Err1 = IOField::define("NormalVelocity");
   int Err2 = IOField::attachData<OMEGA::Array3DR8>("NormalVelocity",NormalVelocityOut);


   auto FieldLayerThickness = ArrayMetaData::create(
       "LayerThickness",
       "Layer thickness on Cell", /// long Name or description
       "m",                           /// units
       "R8StdName",                 /// CF standard Name
       0,                           /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       3,                           /// number of dimensions
       CellDim3DR8                  /// dim pointers
   );
   int Err3 = IOField::define("LayerThickness");
   int Err4 = IOField::attachData<OMEGA::Array3DR8>("LayerThickness",LayerThicknessOut);


   auto FieldSsh = ArrayMetaData::create(
       "ssh",
       "ssh on Cell", /// long Name or description
       "m",                           /// units
       "R8StdName",                 /// CF standard Name
       0,                           /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       3,                           /// number of dimensions
       CellDim3DR8                  /// dim pointers
   );
   int Err5 = IOField::define("ssh");
   int Err6 = IOField::attachData<OMEGA::Array2DR8>("ssh",SshOut);

   //--------------------------------------------------------------------------/

   std::vector<int> OffsetStr(NStrLen, -1);

   int DimStrID;
   int ErrDimStr = OMEGA::IO::defineDim(OutFileID, "StrLen", NStrLen, DimStrID);

   int DimTimeID;
   int ErrDimTime = OMEGA::IO::defineDim(OutFileID, "Time", NTimeLevels, DimTimeID);

   int StrDimIDs[2] = {DimTimeID, DimStrID};
   int ErrStr = IO::defineVar(OutFileID, "xtime", IO::IOTypeChar, 2, StrDimIDs, XtimeID);

   int DimVertID;
   int ErrDimVert = OMEGA::IO::defineDim(OutFileID, "nVertLevels", NVertLevels,DimVertID);

   int DimEdgeID;
   int ErrDim = OMEGA::IO::defineDim(OutFileID, "nEdges", DefDecomp->NEdgesGlobal, DimEdgeID);

   int DimCellID;
   Err = OMEGA::IO::defineDim(OutFileID, "nCells", DefDecomp->NCellsGlobal, DimCellID);

   //--------------------------------------------------------------------------/

   std::vector<int> OffsetEdge(Mesh->NEdgesSize * NVertLevels, -1);
   for (int Edge = 0; Edge < Mesh->NEdgesOwned; ++Edge) {
      int GlobalEdgeAdd = DefDecomp->EdgeIDH(Edge) - 1;
      for (int k = 0; k < NVertLevels; ++k) {
         int VectorAdd         = Edge * NVertLevels + k;
         OffsetEdge[VectorAdd] = GlobalEdgeAdd * NVertLevels + k;
      }
   }

   std::vector<int> EdgeDims{NTimeLevels, DefDecomp->NEdgesGlobal, NVertLevels};
   int ErrDecomp = IO::createDecomp(DecompEdgeR8, IO::IOTypeR8, 3,
                                 EdgeDims, Mesh->NEdgesSize * NVertLevels, OffsetEdge,
                                 IO::DefaultRearr);

   int EdgeDimIDs[3] = {DimTimeID, DimEdgeID, DimVertID};
   int ErrVar = IO::defineVar(OutFileID, "normalVelocity", IO::IOTypeR8, 3, EdgeDimIDs, NormalVelocityID);

   //--------------------------------------------------------------------------/

   std::vector<int> OffsetCell2D(Mesh->NCellsSize * NVertLevels, -1);
   for (int Cell = 0; Cell < Mesh->NCellsOwned; ++Cell) {
      int GlobalCellAdd = DefDecomp->CellIDH(Cell) - 1;
      for (int k = 0; k < NVertLevels; ++k) {
         int VectorAdd         = Cell * NVertLevels + k;
         OffsetCell2D[VectorAdd] = GlobalCellAdd * NVertLevels + k;
      }
   }

   std::vector<int> Cell3Dims{NTimeLevels, DefDecomp->NCellsGlobal, NVertLevels};
   Err = IO::createDecomp(DecompCell3DR8, IO::IOTypeR8, 3,
                                 Cell3Dims, Mesh->NCellsSize * NVertLevels, OffsetCell2D,
                                 IO::DefaultRearr);

   int Cell3DimIDs[3] = {DimTimeID, DimCellID, DimVertID};
   Err = IO::defineVar(OutFileID, "layerThickness", IO::IOTypeR8, 3, Cell3DimIDs, LayerThicknessID);

   //--------------------------------------------------------------------------/

   std::vector<int> OffsetCell1D(Mesh->NCellsSize, -1);
   for (int Cell = 0; Cell < Mesh->NCellsOwned; ++Cell) {
      int GlobalCellAdd = DefDecomp->CellIDH(Cell) - 1;
      int VectorAdd         = Cell ;
      OffsetCell1D[VectorAdd] = GlobalCellAdd ;
   }

   std::vector<int> Cell2Dims{NTimeLevels, DefDecomp->NCellsGlobal};
   Err = IO::createDecomp(DecompCell2DR8, IO::IOTypeR8, 2,
                                 Cell2Dims, Mesh->NCellsSize, OffsetCell1D,
                                 IO::DefaultRearr);
   int Cell2DimIDs[2] = {DimTimeID, DimCellID};

   Err = IO::defineVar(OutFileID, "ssh", IO::IOTypeR8, 2, Cell2DimIDs, SshID);

   //--------------------------------------------------------------------------/

} // sw_init_io

void ShallowWaterCore::
sw_io_write(const Decomp *DefDecomp, HorzMesh *Mesh, OceanState *State) {

   //sw_time_dependent_solution(TestCase, "solution", true, true, CurrentTime, Mesh, State);

//\/*
   TangentialReconOnEdge2D TanReconEdge2D(Mesh);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         const R8 tangentialComp = TanReconEdge2D(IEdge,KLevel, State->NormalVelocity[0]);
         State->NormalVelocity[1](IEdge,KLevel) = tangentialComp;
      });
   //parallelFor(
   //   {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
   //      const R8 tangentialComp = TanReconEdge2D(IEdge,KLevel, NormalVelocitySolution);
   //      TangentialVelocity(IEdge,KLevel) = tangentialComp;
   //   });

   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         const R8 normalComp = State->NormalVelocity[0](IEdge,KLevel);
         const R8 tangentialComp = State->NormalVelocity[1](IEdge,KLevel);
         const R8 zonalComp = normalComp * cos(Mesh->AngleEdge(IEdge))
                            - tangentialComp * sin(Mesh->AngleEdge(IEdge));
         const R8 meridComp = normalComp * sin(Mesh->AngleEdge(IEdge))
                            + tangentialComp * cos(Mesh->AngleEdge(IEdge));
      
         //const R8 normalCompSol = NormalVelocitySolution(IEdge,KLevel);
         //const R8 tangentialCompSol = TangentialVelocity(IEdge,KLevel);
         //const R8 zonalCompSol = normalCompSol * cos(Mesh->AngleEdge(IEdge))
         //                      - tangentialCompSol * sin(Mesh->AngleEdge(IEdge));
         //const R8 meridCompSol = normalCompSol * sin(Mesh->AngleEdge(IEdge))
         //                      + tangentialCompSol * cos(Mesh->AngleEdge(IEdge));


         //const R8 zonalComp = normalComp * cos(Mesh->AngleEdge(IEdge))
         //                   - tangentialComp * sin(Mesh->AngleEdge(IEdge));
         //const R8 meridComp = 0.0;

         //State->NormalVelocity[0](IEdge,KLevel) = zonalComp - zonalCompSol;
         State->NormalVelocity[0](IEdge,KLevel) = zonalComp;
         State->NormalVelocity[1](IEdge,KLevel) = zonalComp;
         NormalVelocityOut(0,IEdge,KLevel) = zonalComp;
         //State->NormalVelocity[0](IEdge,KLevel) = zonalCompSol;
         //State->NormalVelocity[0](IEdge,KLevel) = meridComp;
         //State->NormalVelocity[1](IEdge,KLevel) = 0.0;

      });

   // Save total depth ( H = h + b )
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         State->LayerThicknessH[0](ICell,KLevel) = (State->LayerThicknessH[0](ICell,KLevel)+BottomTopography(ICell));
         LayerThicknessOut(0,ICell,KLevel) = (State->LayerThicknessH[1](ICell,KLevel)+BottomTopography(ICell));
       });   

   if ( TestCase == 22 ) {
   parallelFor(
      {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
         SshOut(0,ICell) = (State->LayerThicknessH[1](ICell,0)-H0);
       });   
   } else {
   parallelFor(
      {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
         SshOut(0,ICell) = (State->LayerThicknessH[1](ICell,0));
       });   
   }

   // Time info
   char xtime[State->NTimeLevels][64] = {{' '}};

   strncpy (xtime[0],"0001-01-01_10:00:00", strlen("0001-01-01_10:00:00")+1);
   //strncpy (xtime[1],"0001-01-01_10:00:00", strlen("0001-01-01_10:00:00")+1);

   // Write Vars
   
   int ErrWriteStr = PIOc_put_var_text(OutFileID, XtimeID, &xtime[0][0]);

   int ErrWriteVel   = IO::writeArray(NormalVelocityOut.data(), Mesh->NEdgesSize * NVertLevels,
                                      &FillR8, OutFileID, DecompEdgeR8, NormalVelocityID);
   int ErrWriteThick = IO::writeArray(LayerThicknessOut.data(), Mesh->NCellsSize * NVertLevels,
                                      &FillR8, OutFileID, DecompCell3DR8, LayerThicknessID);
   int ErrWriteSsh   = IO::writeArray(SshOut.data(), Mesh->NCellsSize,
                                      &FillR8, OutFileID, DecompCell2DR8, SshID);

   // Finished writing, close file
   int ErrClose = IO::closeFile(OutFileID);

   if (ErrClose != 0)
      LOG_ERROR("IO: error closing output file FAIL");

} // sw_io_write

} // namespace OMEGA

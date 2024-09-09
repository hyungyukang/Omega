//===-----------------------------------------------------------------------===/

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "TimeStepper.h"

#include "pio.h"

#include <iostream>
#include <cmath>
#include <iomanip>

//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

   int OutFileID;

   // Define IO field
   int TimeDimID;
   int VertDimID;
   int EdgeDimID;
   int CellDimID;
   int StrlDimID;
   int XtimeID;

   int NormalVelocityID;
   int LayerThicknessID;
   int SshID;

   int DecompCell2DR8;
   int DecompCell3DR8;
   int DecompEdgeR8;

   int NTimeLevelsIO = 1;

   int NStrLen = 64;

   R8 FillR8 = -1.23456789e30;


//===-----------------------------------------------------------------------===/

int IO_init(TimeInstant &Time) {

   int Err;

   MachEnv *DefEnv = MachEnv::getDefault();
   Decomp *DefDecomp     = Decomp::getDefault();
   HorzMesh *Mesh = HorzMesh::getDefault();
   OceanState *State  = OceanState::getDefault();

   //--------------------------------------------------------------------------/

   Err = Field::init();

   //--------------------------------------------------------------------------/

   // IO initialize -----------------------------------------------------------/
   OutFileID = 0;
   int Err0 = IO::openFile(
         OutFileID, "output.nc", IO::ModeWrite, IO::FmtDefault,
         IO::IfExists::Replace); //, IO::Precision::Double);
   if ( Err0 != 0)
      LOG_ERROR("IO: error opening file for output FAILE");

   //--------------------------------------------------------------------------/

   //Err = IO::defineDim(OutFileID,"Time",State->NTimeLevels,TimeDimID);
   Err = IO::defineDim(OutFileID,"StrLen",NStrLen,StrlDimID);
   Err = IO::defineDim(OutFileID,"Time",NTimeLevelsIO,TimeDimID);
   Err = IO::defineDim(OutFileID,"nVertLevels",State->NVertLevels,VertDimID);
   Err = IO::defineDim(OutFileID,"nCells",DefDecomp->NCellsGlobal,CellDimID);
   Err = IO::defineDim(OutFileID,"nEdges",DefDecomp->NEdgesGlobal,EdgeDimID);

   int TimeStrDim[2] = {TimeDimID,StrlDimID};
   Err = IO::defineVar(OutFileID, "xtime", IO::IOTypeChar, 2, TimeStrDim,XtimeID );

   int Edge3DimIDs[3] = {TimeDimID, EdgeDimID, VertDimID};
   Err = IO::defineVar(OutFileID, "normalVelocity", IO::IOTypeR8, 3, Edge3DimIDs, NormalVelocityID);

   int Cell3DimIDs[3] = {TimeDimID, CellDimID, VertDimID};
   Err = IO::defineVar(OutFileID, "layerThickness", IO::IOTypeR8, 3, Cell3DimIDs, LayerThicknessID);
   LOG_INFO("LayerID1 {}",LayerThicknessID);

   int Cell2DimIDs[2] = {TimeDimID, CellDimID};
   Err = IO::defineVar(OutFileID, "ssh", IO::IOTypeR8, 2, Cell2DimIDs, SshID);

     //--------------------------------------------------------------------------/

     std::vector<int> OffsetEdge(Mesh->NEdgesSize * State->NVertLevels, -1);
     for (int Edge = 0; Edge < Mesh->NEdgesOwned; ++Edge) {
        int GlobalEdgeAdd = DefDecomp->EdgeIDH(Edge) - 1;
        for (int k = 0; k < State->NVertLevels; ++k) {
           int VectorAdd         = Edge * State->NVertLevels + k;
           OffsetEdge[VectorAdd] = GlobalEdgeAdd * State->NVertLevels + k;
        }
     }

     std::vector<int> EdgeDims{NTimeLevelsIO, DefDecomp->NEdgesGlobal, State->NVertLevels};
     int ErrDecomp = IO::createDecomp(DecompEdgeR8, IO::IOTypeR8, 3,
                                   EdgeDims, Mesh->NEdgesSize * State->NVertLevels, OffsetEdge,
                                   IO::DefaultRearr);

     //--------------------------------------------------------------------------/

     std::vector<int> OffsetCell2D(Mesh->NCellsSize * State->NVertLevels, -1);
     for (int Cell = 0; Cell < Mesh->NCellsOwned; ++Cell) {
        int GlobalCellAdd = DefDecomp->CellIDH(Cell) - 1;
        for (int k = 0; k < State->NVertLevels; ++k) {
           int VectorAdd         = Cell * State->NVertLevels + k;
           OffsetCell2D[VectorAdd] = GlobalCellAdd * State->NVertLevels + k;
        }
     }

     std::vector<int> Cell3Dims{NTimeLevelsIO, DefDecomp->NCellsGlobal, State->NVertLevels};
     Err = IO::createDecomp(DecompCell3DR8, IO::IOTypeR8, 3,
                                   Cell3Dims, Mesh->NCellsSize * State->NVertLevels, OffsetCell2D,
                                   IO::DefaultRearr);

     //--------------------------------------------------------------------------/

     std::vector<int> OffsetCell1D(Mesh->NCellsSize, -1);
     for (int Cell = 0; Cell < Mesh->NCellsOwned; ++Cell) {
        int GlobalCellAdd = DefDecomp->CellIDH(Cell) - 1;
        int VectorAdd         = Cell ;
        OffsetCell1D[VectorAdd] = GlobalCellAdd ;
     }

     std::vector<int> Cell2Dims{NTimeLevelsIO, DefDecomp->NCellsGlobal};
     Err = IO::createDecomp(DecompCell2DR8, IO::IOTypeR8, 2,
                                   Cell2Dims, Mesh->NCellsSize, OffsetCell1D,
                                   IO::DefaultRearr);

//   //--------------------------------------------------------------------------/

} // IO_init

int IO_write(TimeInstant &Time) {

   int Err;

   MachEnv *DefEnv = MachEnv::getDefault();
   Decomp *DefDecomp     = Decomp::getDefault();
   HorzMesh *Mesh = HorzMesh::getDefault();
   OceanState *State  = OceanState::getDefault();

   // Define Out arrays
   Array3DR8 NormalVelocityOut("NormalVelocityOut", NTimeLevelsIO, Mesh->NEdgesSize, State->NVertLevels);
   Array3DR8 LayerThicknessOut("LayerThicknessOut", NTimeLevelsIO, Mesh->NCellsSize, State->NVertLevels);
   Array2DR8 SshOut("SshOut", NTimeLevelsIO, Mesh->NCellsSize);

   // NormalVelocityOut save
   parallelFor(
      {Mesh->NEdgesAll, State->NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         NormalVelocityOut(0,IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel);
      });
   // LayerThicknessOut save
   parallelFor(
      {Mesh->NCellsAll, State->NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         LayerThicknessOut(0,ICell,KLevel) = State->LayerThickness[0](ICell,KLevel);
         LOG_INFO("LayerThickness {}",LayerThicknessOut(0,ICell,KLevel));
       });
   // SshOut save
   parallelFor(
      {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
         SshOut(0,ICell) = State->LayerThickness[0](ICell,0);
       });


   // If Manufactured Solution, SshOut = SshOut - 1000.0
   bool UseManufacturedSolution = false;
   Config *OmegaConfig = Config::getOmegaConfig();
   Config ManufacturedSolutionConfig("ManufacturedSolution");
   if (OmegaConfig->existsGroup("ManufacturedSolution")) {
      Err = OmegaConfig->get(ManufacturedSolutionConfig);
      if (ManufacturedSolutionConfig.existsVar("UseManufacturedSolution")) {
         Err = ManufacturedSolutionConfig.get("UseManufacturedSolution", UseManufacturedSolution);

         parallelFor(
            {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
               SshOut(0,ICell) = SshOut(0,ICell) - 1000.0;
            });
      }
   }

   // Time info
   char xtime[NTimeLevelsIO][64] = {{' '}};
   std::string CurrTimeString = Time.getString(4,0,"_");
   char *cstr = CurrTimeString.data();
   strncpy (xtime[0],cstr,strlen(cstr)+1);

   // Write Vars
   Err = PIOc_put_var_text(OutFileID, XtimeID, &xtime[0][0]);

   Err = IO::writeArray(NormalVelocityOut.data(), NTimeLevelsIO * Mesh->NEdgesSize * State->NVertLevels,
                        &FillR8, OutFileID, DecompEdgeR8, NormalVelocityID);
   Err = IO::writeArray(LayerThicknessOut.data(), Mesh->NCellsSize * State->NVertLevels,
                        &FillR8, OutFileID, DecompCell3DR8, LayerThicknessID);
   Err = IO::writeArray(SshOut.data(), Mesh->NCellsSize,
                        &FillR8, OutFileID, DecompCell2DR8, SshID);
   // Finished writing, close file
   Err = IO::closeFile(OutFileID);

   if (Err != 0)
      LOG_ERROR("IO: error closing output file FAIL");

} // IO_write

} // namespace OMEGA

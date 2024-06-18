//===--   -----------------------------*- C++ -*-===/
//
//===-----------------------------------------------------------------------===/

#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "OceanDriver.h"
#include "ShallowWaterCore.h"
#include "TendencyTerms.h"
#include "auxiliaryVars/KineticAuxVars.h"
#include "auxiliaryVars/LayerThicknessAuxVars.h"
#include "auxiliaryVars/VorticityAuxVars.h"
#include "mpi.h"
//#include "sw_common.h"
//#include "sw_constants.h"

#include <iostream>
#include <cmath>


//===-----------------------------------------------------------------------===/

namespace OMEGA {

//===-----------------------------------------------------------------------===/

void ShallowWater::
sw_run() {
    
   LOG_INFO("sw_run");

   //HorzMesh *Mesh = HorzMesh::getDefault();
   const auto &Env = MachEnv::getDefaultEnv();
   const auto &MPI_Comm = Env->getComm();
   const auto &Halo = Halo::getDefault();
   const auto &Mesh = HorzMesh::getDefault();
   //State = OceanState::getDefault();
   Decomp *DefDecomp = Decomp::getDefault();
   OceanState *State = OceanState::getDefault();

   NVertLevels = State->NVertLevels;
  
   Array2DR8 RVortVert("RVortVert", Mesh->NVerticesAll, NVertLevels);
   Array2DR8 RVortEdge("RVortEdge", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 AVortVert("AVortVert", Mesh->NVerticesAll, NVertLevels);
   Array2DR8 PVortVert("PVortVert", Mesh->NVerticesAll, NVertLevels);
   Array2DR8 hVertex("hVertex", Mesh->NVerticesAll, NVertLevels);
   Array2DR8 hEdge("hEdge", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 totalHEdge("totalHEdge", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 PVortEdge("PVorteEdge", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 KECell("KECell", Mesh->NCellsAll, NVertLevels);
   Array2DR8 hTend("hTend", Mesh->NCellsAll, NVertLevels);
   Array2DR8 uTend("uTend", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 tanVel("tanVel", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 tanVelLayerThickEdge("tanVelLayerThickEdge", Mesh->NEdgesAll, NVertLevels);
   Array2DR8 normalVelThick("normalVelThick", Mesh->NEdgesAll, NVertLevels);

   Array1DR8 botDepthEdge("botDepthEdge", Mesh->NEdgesAll);

   Array2DR8 totalDepth("totalDepth", Mesh->NCellsAll, NVertLevels);

   RelativeVorticityVertex = RVortVert;
   RelativeVorticityEdge = RVortEdge;
   AbsoluteVorticityVertex = AVortVert;
   PotentialVorticityVertex = PVortVert;
   KineticEnergyCell = KECell;

   TotalDepth     = totalDepth;
   LayerThicknessEdge = hEdge;
   LayerThicknessTend = hTend;
   NormalVelocityTend = uTend;
   TangentialVelocity = tanVel;
   TangentialVelThick = tanVelLayerThickEdge;
   ThicknessFlux = normalVelThick;

   // IO initialize -----------------------------------------------------------/

   //int OutFileID;
   OutFileID = 0;
   int Err0 = IO::openFile(
         OutFileID, "SWoutput.nc", IO::ModeWrite, IO::FmtDefault,
         IO::IfExists::Replace, IO::Precision::Double);
   if ( Err0 != 0)
      LOG_ERROR("IO: error opening file for output FAILE"); 

   //--------------------------------------------------------------------------/

   // Define IO field
   auto EdgeDim          = MetaDim::create("NEdges", Mesh->NEdgesSize);
   auto CellDim          = MetaDim::create("NCells", Mesh->NCellsSize);
   auto VertDim          = MetaDim::create("NVertLevels", NVertLevels);
   std::vector<std::shared_ptr<OMEGA::MetaDim>> CellDim2DR8{CellDim, VertDim};
   std::vector<std::shared_ptr<OMEGA::MetaDim>> EdgeDim2DR8{EdgeDim, VertDim};
   auto FieldNormalVelocity = ArrayMetaData::create(
       "NormalVelocity",
       "Normal Velocity on Edge", /// long Name or description
       "m/s",                         /// units
       "R8StdName",                 /// CF standard Name
       -12345.0,                    /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       2,                           /// number of dimensions
       EdgeDim2DR8                  /// dim pointers
   );
   int Err1 = IOField::define("NormalVelocity");
   int Err2 = IOField::attachData<OMEGA::Array2DR8>("NormalVelocity",State->NormalVelocity[0]);


   auto FieldLayerThickness = ArrayMetaData::create(
       "LayerThickness",
       "Layer thickness on Cell", /// long Name or description
       "m",                           /// units
       "R8StdName",                 /// CF standard Name
       0,                           /// min valid value
       12345.0,                     /// max valid value
       -9.99E+30,                   /// scalar used for undefined entries
       2,                           /// number of dimensions
       CellDim2DR8                  /// dim pointers
   );
   int Err3 = IOField::define("LayerThickness");
   int Err4 = IOField::attachData<OMEGA::Array2DR8>("LayerThickness",State->LayerThickness[0]);

   //--------------------------------------------------------------------------/

   if ( Err1 != 0 || Err2 != 0 || Err3 != 0 || Err4 != 0) {
      LOG_CRITICAL("IOField: Initializing State IO field FAIL");
   } 
   
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
   int DecompEdgeR8;
   int ErrDecomp = IO::createDecomp(DecompEdgeR8, IO::IOTypeR8, 2,
                                 EdgeDims, Mesh->NEdgesSize * NVertLevels, OffsetEdge,
                                 IO::DefaultRearr);

   int DimVertID;
   int ErrDimVert = OMEGA::IO::defineDim(OutFileID, "NVertLevels", NVertLevels,
                                     DimVertID);

   int DimEdgeID;
   int ErrDim = OMEGA::IO::defineDim(OutFileID, "NEdges", DefDecomp->NEdgesGlobal, DimEdgeID);
   int EdgeDimIDs[2] = {DimEdgeID, DimVertID};

   int NormalVelocityID;
   int ErrVar = IO::defineVar(OutFileID, "NormalVelocity", IO::IOTypeR8, 2, EdgeDimIDs, NormalVelocityID);

   //R8 VarMetaR8Ref      = 2.23456789;
   //int ErrMeta = IO::writeMeta("NormalVelocity", VarMetaR8Ref, OutFileID,
   //                           NormalVelocityID);

   // Write array
   R8 FillR8 = -1.23456789e30;
   int ErrWrite = IO::writeArray(State->NormalVelocity[0].data(), Mesh->NEdgesSize * NVertLevels,
                                &FillR8, OutFileID, DecompEdgeR8, NormalVelocityID);
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
   int DecompCellR8;
   Err = IO::createDecomp(DecompCellR8, IO::IOTypeR8, 2,
                                 CellDims, Mesh->NCellsSize * NVertLevels, OffsetCell,
                                 IO::DefaultRearr);

   int DimCellID;
   Err = OMEGA::IO::defineDim(OutFileID, "NCells", DefDecomp->NCellsGlobal, DimCellID);
   int CellDimIDs[2] = {DimCellID, DimVertID};

   int LayerThicknessID;
   Err = IO::defineVar(OutFileID, "LayerThickness", IO::IOTypeR8, 2, CellDimIDs, LayerThicknessID);

   //R8 VarMetaR8Ref      = 2.23456789;
   //int ErrMeta = IO::writeMeta("LayerThickness", VarMetaR8Ref, OutFileID,
   //                             LayerThicknessID);

   // Write array
   Err = IO::writeArray(State->LayerThickness[0].data(), Mesh->NCellsSize * NVertLevels,
                                &FillR8, OutFileID, DecompCellR8, LayerThicknessID);

   //--------------------------------------------------------------------------/


   //LOG_ERROR("IOTest: error creating cell decomp R4 FAIL");
   //Err = OMEGA::IO::createDecomp(DecompCellR8, OMEGA::IO::IOTypeR8, 2,
    //                             CellDims, CellArraySize, OffsetCell,
     //                            OMEGA::IO::DefaultRearr);

   // Time stepping
   for (int step = 0; step < nsteps; ++step){
       const R8 currentTime = step * dt;
       sw_timeStepper(currentTime,dt,Env, Halo,Mesh,State);
   }

   // Finished writing, close file
   int Err5 = IO::closeFile(OutFileID);
   if (Err5 != 0)
      LOG_ERROR("IO: error closing output file FAIL");

   /*
   Array2DI4 &CellsOnCell;      ///< Indx of cells that neighbor each cell
   Array2DI4 &EdgesOnCell;      ///< Indx of edges that border each cell
   Array1DI4 &NEdgesOnCell;      ///< Num of active edges around each cell
   Array2DI4 &VerticesOnCell;      ///< Indx of vertices bordering each cell
   Array2DI4 &CellsOnEdge;      ///< Indx of cells straddling each edge
   Array2DI4 &EdgesOnEdge;      ///< Indx of edges around cells across each edge
   Array1DI4 &NEdgesOnEdge;      ///< Num of edges around the cells across edge
   Array2DI4 &VerticesOnEdge;      ///< Indx of vertices straddling each edge
   Array2DI4 &CellsOnVertex;      ///< Indx of cells that share a vertex
   Array2DI4 &EdgesOnVertex;      ///< Indx of edges sharing vertex as endpoint
   */

   //--- Time stepper ---------------------------------------------------------/

   //--------------------------------------------------------------------------/

} // sw_run

void ShallowWater::
sw_timeStepper(R8 t, R8 dt, const MachEnv *Env, const Halo *Halo, const HorzMesh *Mesh, OceanState *State){

   LOG_INFO("sw_timeStepper");

   // Initialize Tendency
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         NormalVelocityTend(IEdge, KLevel) = 0.0_Real;
      });
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         LayerThicknessTend(IEdge,KLevel) = 0.0_Real;
      });

   // NormalVelocity Tend -----------------------------------------------------/

   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         const I4 ICell0 = Mesh->CellsOnEdge(IEdge,0);
         const I4 ICell1 = Mesh->CellsOnEdge(IEdge,1);
         LayerThicknessEdge(IEdge,KLevel) = 0.5 * (State->LayerThickness[0](ICell0,KLevel)
                                                  +State->LayerThickness[0](ICell1,KLevel));
      });

   // LayerThickEdge
   LayerThicknessAuxVars LayerThicknessAux(Mesh, NVertLevels);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         LayerThicknessAux.computeVarsOnEdge(IEdge, KLevel, 
                                             State->LayerThickness[0], State->NormalVelocity[0]);
      });
   const auto &LayerThicknessEdge = LayerThicknessAux.MeanLayerThickEdge;
   
   // Total depth
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         TotalDepth(ICell,KLevel) = State->LayerThickness[0](ICell,KLevel) + Mesh->BottomDepth(ICell);
      });

   // VorticityAuxVars - RelativeVorticityVertex, and Normalized by layerThickness
   VorticityAuxVars VorticityAux(Mesh, NVertLevels);
   parallelFor(
      {Mesh->NVerticesAll, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int KLevel) {
         VorticityAux.computeVarsOnVertex(IVertex, KLevel,
                                          State->LayerThickness[0], State->NormalVelocity[0]);
      });

   const auto &RelVortVertex        = VorticityAux.RelVortVertex;
   const auto &NormRelVortVertex    = VorticityAux.NormRelVortVertex;
   const auto &NormPlanetVortVertex = VorticityAux.NormPlanetVortVertex;

   // Vorticities : Vertex to edge
   parallelFor(
      {Mesh->NVerticesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         VorticityAux.computeVarsOnEdge(IEdge, KLevel);
      });

   const auto &NormRelVortEdge = VorticityAux.NormRelVortEdge;
   const auto &NormPlanetVortEdge = VorticityAux.NormPlanetVortEdge;

   // KineticAuxVars
   KineticAuxVars KineticAux(Mesh, NVertLevels);
   parallelFor(
       {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          KineticAux.computeVarsOnCell(ICell, KLevel, State->NormalVelocity[0]);
       });
   const auto &KineticEnergyCell = KineticAux.KineticEnergyCell;
   const auto &VelocityDivCell   = KineticAux.VelocityDivCell;

   // Potential Vorticity * h * uPerp
   Config *TendConfig;
   PotentialVortHAdvOnEdge PotVortHAdvOnE(Mesh,TendConfig);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         PotVortHAdvOnE(NormalVelocityTend, IEdge, KLevel, NormRelVortEdge,
                        NormPlanetVortEdge, LayerThicknessEdge, State->NormalVelocity[0]);
      });

   // KEGradOnEdge
   KEGradOnEdge KEGradOnE(Mesh, TendConfig);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         KEGradOnE(NormalVelocityTend, IEdge, KLevel, KineticEnergyCell);
      });

   // SSHGradOnEdge
   SSHGradOnEdge SSHGradOnE(Mesh, TendConfig);
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         SSHGradOnE(NormalVelocityTend, IEdge, KLevel, TotalDepth);
      });


   // LayerThickness Tend -----------------------------------------------------/

   // NormalVelocity * LayerThickEdge
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         ThicknessFlux(IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                      + LayerThicknessEdge(IEdge,KLevel);
      });

   ThicknessFluxDivOnCell ThickFluxDivOnC(Mesh, TendConfig);
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         ThickFluxDivOnC(LayerThicknessTend, ICell, KLevel, ThicknessFlux);
      });


   
   // State advance -----------------------------------------------------------/

   // LayerThickness
   parallelFor(
      {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
         State->LayerThickness[1](ICell,KLevel) = State->LayerThickness[0](ICell,KLevel) 
                                                + LayerThicknessTend(ICell,KLevel);
      });

   // NormalVelocity
   parallelFor(
      {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
         State->NormalVelocity[1](IEdge,KLevel) = State->NormalVelocity[0](IEdge,KLevel)
                                         + NormalVelocityTend(IEdge,KLevel);
      });
   
   // Update time level -------------------------------------------------------/
   LOG_INFO("Before update");
   State->updateTimeLevels();
   
   //ErrW = IO::writeArray(State->NormalVelocity[0], NCellsSize * NVertLevels,
   //                     -1.23456789e30, OutFileID

   //LOG_INFO(RelVortVertex(1,1));

   
//   // RVort
//   CurlOnVertex2D CurlVertex2D(Mesh);
//   for (int KLevel = 0; KLevel < NVertLevels-1 ; KLevel++){
//      parallelFor(
//         {Mesh->NVerticesAll}, KOKKOS_LAMBDA(int IVertex) {
//         RelativeVorticityVertex(IVertex,KLevel) = CurlVertex2D(IVertex, KLevel, State->NormalVelocity[0]);
//      });
//   } 
//
//   parallelFor(
//      {Mesh->NEdgesAll,  NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
//         RelativeVorticityEdge(IEdge,KLevel) = 1.0; 
//      });
//
//

   //parallelFor(
   //   {Mesh->NVerticessAll, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int KLevel) {
   //      PotentialVorticityVertex(IVertex,KLevel) = Relative 
  

   //LOG_INFO(RelativeVorticityVertex(1,1));
   
//   Array1DReal L2Vertex("L2Vertex", Mesh->NVerticesOwned);
//   Array1DReal L2ScaleVertex("L2ScaleVertex", Mesh->NVerticesOwned);
//   CurlOnVertex CurlVertex(Mesh);
//   parallelFor(
//       {Mesh->NVerticesOwned}, KOKKOS_LAMBDA(int IVertex) {
//          // Numerical result
//          const Real CurlNum = CurlVertex(IVertex, VecEdge);
//
//          // Exact result
//          const Real X         = XVertex(IVertex);
//          const Real Y         = YVertex(IVertex);
//          const Real CurlExact = Setup.exactCurlVec(X, Y);
//
//          // Errors
//          LInfVertex(IVertex)      = std::abs(CurlNum - CurlExact);
//          LInfScaleVertex(IVertex) = std::abs(CurlExact);
//          L2Vertex(IVertex) =
//              AreaTriangle(IVertex) * LInfVertex(IVertex) * LInfVertex(IVertex);
//          L2ScaleVertex(IVertex) = AreaTriangle(IVertex) *
//                                   LInfScaleVertex(IVertex) *
//                                   LInfScaleVertex(IVertex);
//       });
//


   // get current state : normalVelocity, layerThickness
    
   // compute vorticity

   // compute tendency

   // stepper + update state
   
   

} // sw_timeStepper

//===-----------------------------------------------------------------------===/

} // namespace OMEGA

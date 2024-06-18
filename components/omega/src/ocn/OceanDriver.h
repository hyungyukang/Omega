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
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OceanConstants.h"
#include "OmegaKokkos.h"
#include "TendencyTerms.h"
#include "mpi.h"
//#include "sw_common.h"
//#include "sw_constants.h"

#include <iostream>
#include <cmath>


namespace OMEGA{

class OceanDriver
{
public:
   OceanDriver () = default;
   //OceanDriver (int argc, char *argv[]);

   //~OceanDriver ();
    
   void get_comm (int argc, char *argv[]);

   void init_IO ();
  
   void init_Decomp(const std::string &MeshFile);

   void init_Halo();

   void init_HorzMesh();

   void init_OceanState();

   void init_Create();

   // A wrapper of the above
   void initialize(int argc,
                   char *argv[],
                   const std::string &MeshFile);

   void set_params();

   void run();

   void finalize ();

//-----------------------------------------------------------------------------/


   I4 NCellsOwned; ///< Number of cells owned by this task
   I4 NCellsAll;   ///< Total number of local cells (owned + all halo)
   I4 NCellsSize;  ///< Array size (incl padding, bndy cell) for cell arrays

   I4 NEdgesOwned;    ///< Number of edges owned by this task
   I4 NEdgesAll;      ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;     ///< Array length (incl padding, bndy) for edge dim

   I4 MaxCellsOnEdge; ///< Max number of cells sharing an edge
   I4 MaxEdges;       ///< Max number of edges around a cell
   I4 MaxEdges2;      ///< Max number of edges around a cell x2

   I4 NVerticesOwned; ///< Number of vertices owned by this task
   I4 NVerticesAll;   ///< Total number (owned+halo) of local vertices
   I4 NVerticesSize;  ///< Array length (incl padding, bndy) for vrtx dim
   I4 VertexDegree;   ///< Number of cells that meet at each vertex

   I4 NVertLevels;

   // Mesh connectivity

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

// int NVertLevels;

//-----------------------------------------------------------------------------/

//-----------------------------------------------------------------------------/



}; // class OmegaDriver



} // namespace OMEGA

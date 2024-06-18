//===--   -----------------------------*- C++ -*-===/
//
//===-----------------------------------------------------------------------===/

//#include "Config.h"
//#include "DataTypes.h"
//#include "Decomp.h"
//#include "Halo.h"
#include "HorzMesh.h"
//#include "HorzOperators.h"
//#include "IO.h"
//#include "Logging.h"
//#include "MachEnv.h"
#include "OceanState.h"
//#include "OceanDriver.h"
//#include "OceanConstants.h"
//#include "OmegaKokkos.h"
//#include "TendencyTerms.h"
//#include "mpi.h"
//#include "sw_common.h"
//#include "sw_constants.h"

#include <iostream>
#include <cmath>


namespace OMEGA{

class ShallowWater
{
public:
   ShallowWater () = default;
    
   //void sw_initialize ();

   void sw_run ();

   void sw_timeStepper(R8 t, R8 dt, const MachEnv *Env, const Halo *Halo, const HorzMesh *Mesh, OceanState *State);

   const R8 bottomDepth = 1000.0;

   const R8 dt = 100.0;
   const R8 initTime = 0.0;
   const R8 endTime = 300.0;

   const I4 nsteps = std::ceil(endTime / dt);

   I4 NVertLevels;

   //Array2DReal normalVelocity;
   //Array2DReal RelativeVorticityVertex

   //Array2DReal RelativeVorticityVertex("RelativeVorticityVertex", Mesh->NVerticesAll, NVertLevels);

   Array2DR8 RelativeVorticityVertex;
   Array2DR8 RelativeVorticityEdge;
   Array2DR8 AbsoluteVorticityVertex;
   Array2DR8 PotentialVorticityVertex;
   Array2DR8 PotentialVorticityEdge;
   Array2DR8 KineticEnergyCell;
   Array2DR8 LayerThicknessVertex; 
   Array2DR8 LayerThicknessEdge; 
   Array2DR8 NormalVelocityTend; 
   Array2DR8 LayerThicknessTend; 

   Array2DR8 TangentialVelocity;
   Array2DR8 TangentialVelThick;

   Array2DR8 TotalDepth;
   Array2DR8 ThicknessFlux;

   Array1DR8 BottomDepthEdge;

   int OutFileID;

   //Array2DReal *RelativeVorticityEdge;
   //Array2DReal *PlanetrayVorticityEdge;
   //Array2DReal *PotentialVorticityEdge;

//-----------------------------------------------------------------------------/

}; // class ShallowWater



} // namespace OMEGA

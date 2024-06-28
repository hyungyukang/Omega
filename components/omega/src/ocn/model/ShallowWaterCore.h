#ifndef OMEGA_SHALLOWWATERCORE_H
#define OMEGA_SHALLOWWATERCORE_H

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
#include "OceanConstants.h"
#include "OceanTestCommon.h"
#include "TendencyTerms.h"
#include "auxiliaryVars/KineticAuxVars.h"
#include "auxiliaryVars/LayerThicknessAuxVars.h"
#include "auxiliaryVars/VorticityAuxVars.h"
#include "mpi.h"

//===-----------------------------------------------------------------------===/

namespace OMEGA{

class ShallowWaterCore {

public:
   ShallowWaterCore () = default;

   //------------------------------------------------------------------------//
   // All from Config later...
   const R8 dt       = 100.0;     // Time step size (sec)
   const R8 initTime = 0.0;       // Model initial time (sec)
   const R8 endTime  = 5*86400;   // Model eEnd time (sec)

   const I4 printInterval = 3600; // Time interval of status check
   const I4 nsteps = std::ceil((endTime) / dt); // Number of time steps

   // Time stepper choices  -------------------------------------------------/
   const char *time_integrator = "heuns"; // (RK2)
   //const char *time_integrator = "forward-euler";
   //const char *time_integrator = "forward-backward";
   //const char *time_integrator = "ssp-rk3";

    // Test cases -----------------------------------------------------------/
      // Stationary time solutions -----------------
   const bool TimeDependentSolution = false; // true If this test is using time-dendent solutions
   //const int TestCase = 0; // Use initial conditions in input file
   const int TestCase = 2; // Global steady-state nonlinear flow
   //const int TestCase = 5; // Zonal flow over an isolated moutain

      // Time-dependent solution -------------------
   //const bool TimeDependentSolution = true; // true If this test is using time-dendent solutions
   //const int TestCase = 21; // Solid body rotation; time-dependent solution = true

   // Check status ----------------------------------------------------------/
   const bool ComputeNormError = true;   // compute L2 norm error
   const bool ComputeTotalEnergy = true; // compute total energy

   //------------------------------------------------------------------------//

   // SW initialize vars
   void sw_init_var(const HorzMesh *Mesh, const OceanState *State);

   // SW initialize IO
   void sw_init_io(const Decomp *DefDecomp, const HorzMesh *Mesh, const OceanState *State);

   // SW initial field create
   void sw_init_field (const int TestCase, const MachEnv *Env, const Halo *Halo, 
                       const HorzMesh *Mesh, OceanState *State);

   // SW time-dependent solution
   void sw_time_dependent_solution(const int TestCase, R8 t,
                                 const HorzMesh *Mesh, OceanState *State);

   // SW model core
   void sw_model(int Comm, const MachEnv *Env,const Decomp *DefDecomp, const Halo *DefHalo, 
                 const HorzMesh *Mesh, OceanState *State);

   // SW time stepper
   void sw_time_stepper(const std::string &time_integrator, R8 t, const R8 dt, const MachEnv *Env, 
                        const Halo *Halo, const HorzMesh *Mesh, OceanState *State);

   // SW layerThickness tendency
   void sw_tend_thick(int ThickCurTimeLevel, int VelCurTimeLevel, const MachEnv *Env, 
                      const Halo *Halo, const HorzMesh *Mesh, OceanState *State);

   // SW normalVelocity tendency
   void sw_tend_vel(int ThickCurTimeLevel, int VelCurTimeLevel, const MachEnv *Env, 
                    const Halo *Halo, const HorzMesh *Mesh, OceanState *State);

   // SW IO write
   void sw_io_write(const Decomp *DefDecomp, const HorzMesh *Mesh, const OceanState *State);

   // Status check
   void sw_check_status(const int printInterval, R8 t, const R8 dt, int Comm,
                        const MachEnv *Env, const HorzMesh *Mesh, OceanState *State);
  
   //------------------------------------------------------------------------//
   // SW constants

   I4 NVertLevels;

   I4 NEdgesSize;
   I4 NCellsSize;
   I4 NVerticesSize;
   I4 NEdgesAll;
   I4 NCellsAll;
   I4 NVerticesAll;
   I4 NEdgesOwned;
   I4 NCellsOwned;
   I4 NVerticesOwned;

   //------------------------------------------------------------------------//
   // SW vars
   Array2DReal LayerThicknessSolution;
   Array2DReal NormalVelocitySolution;
   Array2DReal TotalDepthKECell;
   Array2DReal LayerThicknessEdge;
   Array2DReal ThicknessFlux;
   Array2DReal TangentialVelocity;
   Array2DReal NormalVelocityTend;
   Array2DReal LayerThicknessTend;

   Array1DReal BottomTopography;

   //------------------------------------------------------------------------//
   // IO-related vars

   int OutFileID;
   int DecompCellR8;
   int DecompEdgeR8;
   int LayerThicknessID;
   int NormalVelocityID;
   R8 FillR8 = -1.23456789e30;

}; // class ShallowWaterCore

} // namespace OMEGA
#endif

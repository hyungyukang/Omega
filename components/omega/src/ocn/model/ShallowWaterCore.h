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
#include "TimeMgr.h"
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
   const R8 dt       =  5.0;    // Time step size (sec)
   const R8 initTime = 0.0;       // Model initial time (sec)
   const R8 endTime  = 10.0*3600.0;   // Model eEnd time (sec)
   //const R8 endTime  = 15*86400;   // Model eEnd time (sec)
   //const R8 endTime  = 10.0*3600.0;   // Model eEnd time (sec)
   //const R8 endTime  = dt;   // Model eEnd time (sec)
   //const I4 printInterval = endTime; // Time interval of status check
   const I4 printInterval = 3600.0; // Time interval of status check
   const I4 nsteps = std::ceil((endTime) / dt); // Number of time steps
   R8 CurrentTime = 0.0;

   // Time manager  ---------------------------------------------------------/
   //Calendar CalGreg("Gregorian", OMEGA::CalendarGregorian);
   //TimeInstant Time0(&CalGreg, 1, 1, 1, 0, 0, 0.0);
   //TimeInstant StartTime = Time0;
   //TimeInstant EndTime = (&CalGreg, 1, 1, 1, 0, 0, endTime);

   //const I4 nsteps = 100; // Number of time steps
   //const R8 endTime  = nsteps*dt;   // Model eEnd time (sec)
   //const I4 nsteps = std::ceil((endTime) / dt); // Number of time steps

   // Time stepper choices  -------------------------------------------------/
   const char *time_integrator = "midpoint"; // (~RK2)
   //const char *time_integrator = "heuns"; // (~RK2)
   //const char *time_integrator = "forward-euler";
   //const char *time_integrator = "forward-backward";
   //const char *time_integrator = "ssp-rk3";
   //const char *time_integrator = "rk4";

    // Test cases -----------------------------------------------------------/

    const bool initFieldFromFile = true; // If true, use initial conditions in input file

      // Stationary time solutions -----------------
   //const bool TimeDependentSolution = false; // true If this test is using time-dendent solutions
   //const int TestCase = 0; // Use initial conditions in input file
   //const int TestCase = 2; // Global steady-state nonlinear flow
   //const int TestCase = 5; // Zonal flow over an isolated moutain

      // Time-dependent solution -------------------
   const bool TimeDependentSolution = true; // true If this test is using time-dendent solutions
   //const int TestCase = 21; // Solid body rotation; time-dependent solution = true
   const int TestCase = 22; // Nonlinear manufactured solution; time-dependent solution = true

   // Check status ----------------------------------------------------------/
   const bool ComputeNormError = true;   // compute L2 norm error
   const bool ComputeTotalEnergy = true; // compute total energy

   //------------------------------------------------------------------------//

   // SW initialize vars
   void sw_init_var(HorzMesh *Mesh, const OceanState *State);

   // SW initialize IO
   void sw_init_io(const Decomp *DefDecomp, HorzMesh *Mesh, OceanState *State);

   // SW initial field create
   void sw_init_field (const int TestCase, const MachEnv *Env, const Halo *Halo, 
                       HorzMesh *Mesh, OceanState *State);

   // SW time-dependent solution
   void sw_time_dependent_solution(const int TestCase, const std::string &solutionOpt,
                                   bool tendThickSourceTerm, bool tendVelSourceTerm,
                                   const R8 nowTime, HorzMesh *Mesh, OceanState *State);

   // SW model core
   void sw_model(int Comm, const MachEnv *Env,const Decomp *DefDecomp, Halo *DefHalo, 
                 HorzMesh *Mesh, OceanState *State);

   // SW time stepper
   void sw_time_stepper(const std::string &time_integrator, R8 t, const R8 dt, const MachEnv *Env, 
                        Halo *Halo, HorzMesh *Mesh, OceanState *State);

   // SW layerThickness tendency
   void sw_tend_thick(int ThickCurTimeLevel, int VelCurTimeLevel, const R8 nowTime, const MachEnv *Env, 
                      const Halo *Halo, HorzMesh *Mesh, OceanState *State);

   // SW normalVelocity tendency
   void sw_tend_vel(int ThickCurTimeLevel, int VelCurTimeLevel, const R8 nowTime, const MachEnv *Env, 
                    const Halo *Halo, HorzMesh *Mesh, OceanState *State);

   // SW IO write
   void sw_io_write(const Decomp *DefDecomp, HorzMesh *Mesh, OceanState *State);

   // Status check
   void sw_check_status(const int printInterval, R8 t, const R8 dt, int Comm,
                        const MachEnv *Env, HorzMesh *Mesh, OceanState *State);
  
   //------------------------------------------------------------------------//
   // SW constants

   I4 NStrLen = 64;
   I4 NTimeLevels = 1;
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

   R8 H0;

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
   Array2DReal NormalVelocityRKTemp;
   Array2DReal LayerThicknessRKTemp;

   Array2DReal SshOut;
   Array3DReal LayerThicknessOut;
   Array3DReal NormalVelocityOut;

   Array1DReal BottomTopography;

   //------------------------------------------------------------------------//
   // IO-related vars

   int OutFileID;
   int DecompStr;
   int DecompCell2DR8;
   int DecompCell3DR8;
   int DecompEdgeR8;
   int XtimeID;
   int LayerThicknessID;
   int NormalVelocityID;
   int SshID;
   R8 FillR8 = -1.23456789e30;

}; // class ShallowWaterCore

} // namespace OMEGA
#endif

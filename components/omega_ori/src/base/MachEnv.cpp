//===-- base/MachEnv.cpp - machine environment methods ----------*- C++ -*-===//
//
// The machine environment defines a number of parameters associated with
// the message-passing and threading environments. It also can describe
// the node-level hardware environment, including any useful accelerator
// or processor-level constants. Multiple machine environments can be defined,
// mostly associated with subsets of the processor decomposition, but
// a default machine env must be created very early in model
// initialization. All environments can be retrieved by name, though a
// specific retrieval for the most common default environment is provided
// for efficiency.
//
//===----------------------------------------------------------------------===//

#include "MachEnv.h"
#include "mpi.h"

#include <map>
#include <string>
// Note that we should replace iostream and std::cerr with the logging
// capability once that is enabled.
#include <iostream>

namespace OMEGA {

// create the static class members
MachEnv *MachEnv::DefaultEnv = nullptr;
std::map<std::string, MachEnv> MachEnv::AllEnvs;

// Constructors
//------------------------------------------------------------------------------
// Constructor that creates an environment using an input communicator.
// This is typically used to set the default OMEGA communicator

MachEnv::MachEnv(const std::string Name, // [in] name of environment
                 const MPI_Comm InComm   // [in] parent MPI communicator
) {

   // Check to see if an environment of the same name already exists and
   // if so, exit with an error
   if (AllEnvs.find(Name) != AllEnvs.end()) {
      std::cerr << "Attempted to create a MachEnv with name " << Name
                << " but an Env of that name already exists ";
      return;
   }

   // Set the communicator to the input communicator by duplicating it
   MPI_Comm_dup(InComm, &Comm);

   // get task ID and set local MPI task
   MPI_Comm_rank(Comm, &MyTask);

   // get total number of MPI tasks
   MPI_Comm_size(Comm, &NumTasks);

   // Set task 0 as master
   MasterTask = 0;

   // determine if this task is the master task
   if (MyTask == MasterTask) {
      MasterTaskFlag = true;
   } else {
      MasterTaskFlag = false;
   }

   // All tasks are members of this communicator's group
   MemberFlag = true;

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NumThreads = omp_get_num_threads();
#else
   NumThreads = 1;
#endif

   // Add this environment to the list of environments
   AllEnvs.emplace(Name, *this);

} // end constructor with MPI communicator

//------------------------------------------------------------------------------
// Construct a new environment from a contiguous subset of an
// existing environment starting from same root task. The master task
// can be changed from the default (0) as the optional argument InMasterTask.

MachEnv::MachEnv(const std::string Name, // [in] name of environment
                 const MachEnv *InEnv,   // [in] parent MachEnv
                 const int NewSize,      // [in] num tasks in new env
                 int InMasterTask        // [in] optionally set Master Task
) {

   // Check to see if an environment of the same name already exists and
   // if so, exit with an error
   if (AllEnvs.find(Name) != AllEnvs.end()) {
      std::cerr << "Attempted to create a MachEnv with name " << Name
                << " but an Env of that name already exists ";
      return;
   }

   // Check for valid master task input
   if (InMasterTask < 0 || InMasterTask >= NewSize) {
      std::cerr << "Invalid MasterTask " << InMasterTask
                << " sent to constructor: must be [0,NewSize-1]" << std::endl;
      return;
   }

   // Get communicator from old environment
   MPI_Comm InComm = InEnv->getComm();

   // Error checks on new size
   int OldSize;
   MPI_Comm_size(InComm, &OldSize);
   if (NewSize > OldSize) {
      std::cerr << "Invalid NewSize in MachEnv constructor "
                << "NewSize = " << NewSize << " is larger than old size"
                << OldSize << std::endl;
      return;
   }

   // First retrieve the group associated with the old environment
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Define the range of tasks in new group (0, NewSize-1)
   int NRanges     = 1;
   int LastTask    = NewSize - 1;
   int Range[1][3] = {0, LastTask, 1};

   // Create a new group with the new range
   MPI_Group NewGroup;
   MPI_Group_range_incl(InGroup, NRanges, Range, &NewGroup);

   // Create the communicator for the new group
   int Err = MPI_Comm_create(InComm, NewGroup, &Comm);
   if (Err != MPI_SUCCESS) {
      std::cerr << "Error creating new communicator in MachEnv constructor"
                << std::endl;
      return;
   }

   // if this task is part of the new group/communicator,
   // initialize all values
   if (Comm != MPI_COMM_NULL) {

      MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(Comm, &MyTask);
      // get total number of MPI tasks
      MPI_Comm_size(Comm, &NumTasks);
      // Set master task to either the optional input or 0 as default
      MasterTask = InMasterTask;
      // determine if this task is the master task
      if (MyTask == MasterTask) {
         MasterTaskFlag = true;
      } else {
         MasterTaskFlag = false;
      }

      // otherwise initialize all values to bogus values
   } else {

      MemberFlag     = false;
      MyTask         = -999;
      NumTasks       = -999;
      MasterTask     = -999;
      MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NumThreads = omp_get_num_threads();
#else
   NumThreads = 1;
#endif

   // Add this environment to the list of environments
   AllEnvs.emplace(Name, *this);

} // end subset constructor with contiguous range

//------------------------------------------------------------------------------
// Construct a new environment from a strided subset of an
// existing environment. The master task can be changed from the default (0)
// with an optional input argument InMasterTask.

MachEnv::MachEnv(const std::string Name, // [in] name of environment
                 const MachEnv *InEnv,   // [in] parent MachEnv
                 const int NewSize,      // [in] num tasks in new env
                 const int Begin,        // [in] starting parent task
                 const int Stride,       // [in] stride for tasks to incl
                 const int InMasterTask  // [in] optionally set Master Task
) {

   // Check to see if an environment of the same name already exists and
   // if so, exit with an error
   if (AllEnvs.find(Name) != AllEnvs.end()) {
      std::cerr << "Attempted to create a MachEnv with name " << Name
                << " but an Env of that name already exists ";
      return;
   }

   // Check for valid master task input
   if (InMasterTask < 0 || InMasterTask >= NewSize) {
      std::cerr << "Invalid MasterTask " << InMasterTask
                << " sent to constructor: must be [0,NewSize-1]" << std::endl;
      return;
   }

   // First retrieve the group associated with the parent environment
   MPI_Comm InComm = InEnv->getComm();
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Define the range of tasks in new group based on the input
   // start task and stride
   int NRanges     = 1;
   int LastTask    = Begin + (NewSize - 1) * Stride;
   int Range[1][3] = {Begin, LastTask, Stride};

   // Error checks for valid range
   int OldSize;
   MPI_Comm_size(InComm, &OldSize);
   if (LastTask >= OldSize) {
      std::cerr << "Invalid range in strided MachEnv constructor "
                << "LastTask = " << LastTask << " is larger than old size"
                << OldSize << std::endl;
      return;
   }

   // Create a new group with the new range
   MPI_Group NewGroup;
   MPI_Group_range_incl(InGroup, NRanges, Range, &NewGroup);

   // Create the communicator for the new group
   int Err = MPI_Comm_create(InComm, NewGroup, &Comm);
   if (Err != MPI_SUCCESS) {
      std::cerr << "Error creating new communicator in MachEnv constructor"
                << std::endl;
      return;
   }

   // if this task is part of the new group/communicator
   // initialize all members
   if (Comm != MPI_COMM_NULL) { // this task is part of the group

      MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(Comm, &MyTask);
      // get total number of MPI tasks
      MPI_Comm_size(Comm, &NumTasks);
      // Set master task to either the input value or default to 0
      MasterTask = InMasterTask;
      // determine if this task is the master task
      if (MyTask == MasterTask) {
         MasterTaskFlag = true;
      } else {
         MasterTaskFlag = false;
      }

      // otherwise, set all members to bogus values
   } else {

      MemberFlag     = false;
      MyTask         = -999;
      NumTasks       = -999;
      MasterTask     = -999;
      MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NumThreads = omp_get_num_threads();
#else
   NumThreads = 1;
#endif

   // Add this environment to the list of environments
   AllEnvs.emplace(Name, *this);

} // end constructor using strided range

//------------------------------------------------------------------------------
// Construct a new environment from a custom subset of an
// existing environment, supplying list of parent tasks to include
// The master task can be changed from the default (0) with an optional
// argument InMasterTask.

MachEnv::MachEnv(const std::string Name, // [in] name of environment
                 const MachEnv *InEnv,   // [in] parent MPI communicator
                 const int NewSize,      // [in] num tasks in new env
                 const int Tasks[],      // [in] vector of parent tasks to incl
                 const int InMasterTask  // [in] optionally set Master Task
) {

   // Check to see if an environment of the same name already exists and
   // if so, exit with an error
   if (AllEnvs.find(Name) != AllEnvs.end()) {
      std::cerr << "Attempted to create a MachEnv with name " << Name
                << " but an Env of that name already exists ";
      return;
   }

   // Check for valid master task input
   if (InMasterTask < 0 || InMasterTask >= NewSize) {
      std::cerr << "Invalid MasterTask " << InMasterTask
                << " sent to constructor: must be [0,NewSize-1]" << std::endl;
      return;
   }

   // Check the input task list to make sure all are valid
   MPI_Comm InComm = InEnv->getComm();
   int OldSize;
   MPI_Comm_size(InComm, &OldSize);
   for (int i = 0; i < NewSize; ++i) {
      int ThisTask = Tasks[i];
      if (ThisTask < 0 || ThisTask >= OldSize) {
         std::cerr << "Invalid task in MachEnv constructor "
                   << "task = " << ThisTask << " is < 0 or larger than old size"
                   << OldSize << std::endl;
         return;
      }
   }

   // First retrieve the group associated with the input communicator
   MPI_Group InGroup;
   MPI_Comm_group(InComm, &InGroup);

   // Create a new group with the selected tasks
   MPI_Group NewGroup;
   MPI_Group_incl(InGroup, NewSize, Tasks, &NewGroup);

   // Create the communicator for the new group
   int Err = MPI_Comm_create(InComm, NewGroup, &Comm);
   if (Err != MPI_SUCCESS) {
      std::cerr << "Error creating new communicator in MachEnv constructor"
                << std::endl;
      return;
   }

   // If this task is part of the new group/communicator
   // initialize all class members
   if (Comm != MPI_COMM_NULL) { // member of new group

      MemberFlag = true;
      // get task ID for local MPI task
      MPI_Comm_rank(Comm, &MyTask);
      // get total number of MPI tasks
      MPI_Comm_size(Comm, &NumTasks);
      // set master task to either input value or default to 0
      MasterTask = InMasterTask;
      // determine if this task is the master task
      if (MyTask == MasterTask) {
         MasterTaskFlag = true;
      } else {
         MasterTaskFlag = false;
      }

      // otherwise, set all members to bogus values
   } else {

      MemberFlag     = false;
      MyTask         = -999;
      NumTasks       = -999;
      MasterTask     = -999;
      MasterTaskFlag = false;
   }

#ifdef OMEGA_THREADED
   // total number of OpenMP threads
   NumThreads = omp_get_num_threads();
#else
   NumThreads = 1;
#endif

   // Add this environment to the list of environments
   AllEnvs.emplace(Name, *this);

} // end constructor with selected tasks

//------------------------------------------------------------------------------
// Initializes the Machine Environment by creating the DefaultEnv for Omega

void MachEnv::init(const MPI_Comm InComm // [in] communicator to use
) {

   // Create the default environment based on the input communicator
   MachEnv DefEnv("Default", InComm);

   // Retrieve this environment and set pointer to DefaultEnv
   MachEnv::DefaultEnv = getEnv("Default");

} // end init Mach Env

// Remove/delete functions
//------------------------------------------------------------------------------
// Remove environment

void MachEnv::removeEnv(const std::string Name // [in] name of env to remove
) {

   AllEnvs.erase(Name);

} // end removeEnv

//------------------------------------------------------------------------------
// Remove all environments

void MachEnv::removeAll() {

   AllEnvs.clear(); // removes all environments by removing them from map

} // end removeAll

// Retrieval functions
//------------------------------------------------------------------------------
// Get default environment
MachEnv *MachEnv::getDefaultEnv() { return MachEnv::DefaultEnv; }

//------------------------------------------------------------------------------
// Get environment by name
MachEnv *MachEnv::getEnv(const std::string Name ///< [in] Name of environment
) {

   // look for an instance of this name
   auto it = AllEnvs.find(Name);

   // if found, return the environment pointer
   if (it != AllEnvs.end()) {
      return &(it->second);

      // otherwise print an error and return a null pointer
   } else {
      std::cerr << "Attempt to retrieve a non-existent MachEnv " << Name
                << " has not been defined or has been removed";
      return nullptr;
   }

} // end getEnv

//------------------------------------------------------------------------------
// Get communicator for an environment
MPI_Comm MachEnv::getComm() const { return Comm; }

//------------------------------------------------------------------------------
// Get local task/rank ID
int MachEnv::getMyTask() const { return MyTask; }

//------------------------------------------------------------------------------
// Get total number of MPI tasks/ranks
int MachEnv::getNumTasks() const { return NumTasks; }

//------------------------------------------------------------------------------
// Get task ID for the master task (typically 0)
int MachEnv::getMasterTask() const { return MasterTask; }

//------------------------------------------------------------------------------
// Determine whether local task is the master

bool MachEnv::isMasterTask() const { return MasterTaskFlag; }

//------------------------------------------------------------------------------
// Determine whether local task is in this communicator's group

bool MachEnv::isMember() const { return MemberFlag; }

//------------------------------------------------------------------------------
// Set task ID for the master task (if not 0)

int MachEnv::setMasterTask(const int TaskID) {

   int Err = 0;

   // If called from outside the group, don't do anything
   if (!MemberFlag)
      return Err;

   // Check for valid input and reset the master if valid
   if (TaskID >= 0 && TaskID < NumTasks) {
      MasterTask = TaskID;
      if (MyTask == MasterTask) {
         MasterTaskFlag = true;
      } else {
         MasterTaskFlag = false;
      }
   } else {
      std::cerr << "Error: invalid TaskID sent to MachEnv.setMasterTask "
                << TaskID << std::endl;
      Err = -1;
   }
   return Err;

} // end setMasterTask

//------------------------------------------------------------------------------
// Print all members of a MachEnv

void MachEnv::print() const {

   std::cout << "  Comm           = " << Comm << std::endl;
   std::cout << "  MyTask         = " << MyTask << std::endl;
   std::cout << "  NumTasks       = " << NumTasks << std::endl;
   std::cout << "  MasterTask     = " << MasterTask << std::endl;
   std::cout << "  MasterTaskFlag = " << MasterTaskFlag << std::endl;
   std::cout << "  MemberFlag     = " << MemberFlag << std::endl;
   std::cout << "  NumThreads     = " << NumThreads << std::endl;
   std::cout << "  VecLength      = " << VecLength << std::endl;

} // end print

} // end namespace OMEGA

//===----------------------------------------------------------------------===//

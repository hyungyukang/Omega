//===-- SundialsStepper.cpp - Runge-Kutta steps from SUNDIALS ---*- C++ -*-===//
//
// Methods for the time stepper that takes explicit Runge-Kutta steps with the
// ERKStep integrator of the SUNDIALS ARKODE package.
//
// Each step mirrors RungeKutta4Stepper, with ARKODE forming the stages:
//  - the start state {h, u, h*T} is the current time level, with the
//    thickness-weighted tracers held in the next tracer time level, which is
//    free until the result is written there
//  - ARKODE is reset to that state at t = 0 and takes one fixed step of
//    length TimeStep into a vector wrapping the next time level, so every
//    intermediate stage is formed directly in Omega's arrays
//  - before each intermediate stage is used, prepareStage prescribes the
//    state, computes tracer concentrations and exchanges halos, and the
//    right-hand side computes Omega's tendencies at the stage
//  - the work done once per step (KPP fields, implicit vertical mixing,
//    validation) stays outside ARKODE
//
//===----------------------------------------------------------------------===//

#include "SundialsStepper.h"
#include "Config.h"
#include "Logging.h"
#include "Pacer.h"

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_mpimanyvector.h>

// SUNDIALS is built once, in the precision of Omega's Real. The
// single-precision copy of OmegaLib that some tests build cannot share
// vectors with it, so there this stepper only reports the mismatch.
#if defined(OMEGA_SINGLE_PRECISION) == defined(SUNDIALS_SINGLE_PRECISION)
#define OMEGA_SUNDIALS_MATCHES_REAL
#endif

namespace OMEGA {

namespace {

// Shorthand method names and the SUNDIALS tables they select
const std::map<std::string, std::string> MethodShorthands = {
    {"ForwardEuler", "ARKODE_FORWARD_EULER_1_1"},
    {"RK2", "ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2"},
    {"RK3", "ARKODE_SHU_OSHER_3_2_3"},
    {"RK4", "ARKODE_KUTTA_RK4a_4_4"}};

// Aborts with ARKODE's name for a failed call's return flag
void checkArkode(int Flag, const std::string &Call) {
   if (Flag < 0)
      ABORT_ERROR("SundialsStepper: {} failed with {}", Call,
                  ARKodeGetReturnFlagName(Flag));
}

// Builds the Butcher table given in the CustomTable configuration group
ARKodeButcherTable createCustomTable(Config &SundialsConfig) {
   Config TableConfig("CustomTable");
   Error Err = SundialsConfig.get(TableConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration::Sundials::CustomTable group not "
                          "found in Config; it is required when Method is "
                          "Custom");

   I4 NStages = 0;
   I4 Order   = 0;
   std::vector<R8> A;
   std::vector<R8> B;
   std::vector<R8> C;
   Err += TableConfig.get("Stages", NStages);
   Err += TableConfig.get("Order", Order);
   Err += TableConfig.get("A", A);
   Err += TableConfig.get("B", B);
   Err += TableConfig.get("C", C);
   CHECK_ERROR_ABORT(Err, "Sundials::CustomTable requires Stages, Order, A, "
                          "B and C");

   if (NStages < 1 || A.size() != static_cast<size_t>(NStages * NStages) ||
       B.size() != static_cast<size_t>(NStages) ||
       C.size() != static_cast<size_t>(NStages))
      ABORT_ERROR("Sundials::CustomTable with {} stages needs {} values in A "
                  "and {} in B and C, but has {}, {} and {}",
                  NStages, NStages * NStages, NStages, A.size(), B.size(),
                  C.size());

   for (int I = 0; I < NStages; ++I) {
      for (int J = I; J < NStages; ++J) {
         if (A[I * NStages + J] != 0)
            ABORT_ERROR("Sundials::CustomTable is not explicit: A({},{}) = {} "
                        "must be zero",
                        I, J, A[I * NStages + J]);
      }
   }

   std::vector<sunrealtype> ATable(A.begin(), A.end());
   std::vector<sunrealtype> BTable(B.begin(), B.end());
   std::vector<sunrealtype> CTable(C.begin(), C.end());
   ARKodeButcherTable Table = ARKodeButcherTable_Create(
       NStages, Order, 0, CTable.data(), ATable.data(), BTable.data(), nullptr);
   OMEGA_REQUIRE(Table, "ARKodeButcherTable_Create failed for CustomTable");

   // A negative flag means the coefficients do not reach the stated order
   int MeasuredOrder = 0;
   int EmbedOrder    = 0;
   if (ARKodeButcherTable_CheckOrder(Table, &MeasuredOrder, &EmbedOrder,
                                     nullptr) < 0)
      ABORT_ERROR("Sundials::CustomTable satisfies the order conditions only "
                  "to order {}, not the stated order {}",
                  MeasuredOrder, Order);

   return Table;
}

// Right-hand side called by ARKODE: Omega's tendencies at the current stage.
// ARKODE evaluates the first stage on its own copy of the start state and
// later stages on the vector passed to ARKodeEvolve, which wraps the next
// time level.
int stepperRHS(sunrealtype T, N_Vector Y, N_Vector YDot, void *UserData) {
   const auto *Stepper = static_cast<const SundialsStepper *>(UserData);
   const int Stage     = Stepper->getStageIndex();
   OMEGA_REQUIRE((Stage == 0) == (Y != Stepper->getStepEnd()),
                 "SundialsStepper: ARKODE evaluated stage {} on an "
                 "unexpected vector",
                 Stage);
   Stepper->computeStageTendencies(Stage, YDot);
   return 0;
}

// Called by ARKODE once it has formed an intermediate stage
int stepperStage(sunrealtype T, N_Vector Y, void *UserData) {
   const auto *Stepper = static_cast<const SundialsStepper *>(UserData);
   OMEGA_REQUIRE(Y == Stepper->getStepEnd(),
                 "SundialsStepper: ARKODE formed a stage outside the next "
                 "time level");
   Stepper->prepareStage(Stepper->getStageIndex());
   return 0;
}

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
// Device data of one piece of a state vector
Real *pieceData(N_Vector Vec, int Piece) {
   return N_VGetDeviceArrayPointer(N_VGetSubvector_MPIManyVector(Vec, Piece));
}
#endif

} // namespace

//------------------------------------------------------------------------------
// Constructor creates an instance of the SUNDIALS stepper, fills in some time
// information and loads the Butcher table chosen in the configuration. Data
// pointers are added later.
SundialsStepper::SundialsStepper(
    const std::string &InName,              // [in] name of time stepper
    const TimeInterval &InTimeStep,         // [in] time step
    const TimeStepperStartType InStartType, // [in] option for starting
    const TimeInstant &InStartTime,         // [in] start time for full sim
    const TimeStepperStopType InStopType,   // [in] option for stopping
    std::optional<TimeInstant> InStopTime,  // [in] stop time if option AtTime
    std::optional<TimeInterval> InDuration  // [in] duration if AfterDuration
    )
    : TimeStepper(InName, TimeStepperType::Sundials, 2, InTimeStep, InStartType,
                  InStartTime, InStopType, InStopTime, InDuration),
      StepStart(InStartTime) {

#ifndef OMEGA_SUNDIALS_MATCHES_REAL
   ABORT_ERROR("TimeStepper {}: the Sundials time stepper needs SUNDIALS "
               "built in the precision of Omega's Real",
               Name);
#endif

   loadTable();
}

//------------------------------------------------------------------------------
// Destructor frees the SUNDIALS integrator, vectors and table. The state
// vectors only reference the wrapped arrays, which are freed last.
SundialsStepper::~SundialsStepper() {
   if (ArkMem)
      ARKodeFree(&ArkMem);
   for (auto &Entry : StateVectors) {
      N_VDestroy(Entry.second);
   }
   ArrayVectors.clear();
   if (Table)
      ARKodeButcherTable_Free(Table);
   if (Ctx)
      SUNContext_Free(&Ctx);
}

//------------------------------------------------------------------------------
// Loads the Butcher table chosen in the configuration and derives the stage
// and projection fractions of the time step from it
void SundialsStepper::loadTable() {

   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in SundialsStepper");
   Config TimeIntConfig("TimeIntegration");
   Error Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   Config SundialsConfig("Sundials");
   Err = TimeIntConfig.get(SundialsConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration::Sundials group not found in "
                          "Config; it is required by the Sundials time "
                          "stepper");

   std::string Method;
   Err = SundialsConfig.get("Method", Method);
   CHECK_ERROR_ABORT(Err, "Method not found in TimeIntegration::Sundials "
                          "Config");

   if (Method == "Custom") {
      MethodName = Method;
      Table      = createCustomTable(SundialsConfig);
   } else {
      auto Shorthand = MethodShorthands.find(Method);
      MethodName =
          Shorthand != MethodShorthands.end() ? Shorthand->second : Method;
      Table = ARKodeButcherTable_LoadERKByName(MethodName.c_str());
      if (!Table)
         ABORT_ERROR("TimeStepper {}: unknown Sundials Method {}; use "
                     "ForwardEuler, RK2, RK3, RK4, Custom or the name of a "
                     "SUNDIALS explicit Butcher table, eg "
                     "ARKODE_HEUN_EULER_2_1_2",
                     Name, Method);
   }

   // Omega projects the pseudo-thickness from the start of the step to the
   // next stage, or to the end of the step after the last stage, so each
   // stage after the first must lie after the start of the step
   const int NStages = Table->stages;
   StageFrac.resize(NStages);
   ProjFrac.resize(NStages);
   for (int Stage = 0; Stage < NStages; ++Stage) {
      StageFrac[Stage] = Table->c[Stage];
      ProjFrac[Stage]  = Stage + 1 < NStages ? Table->c[Stage + 1] : 1;
      if (ProjFrac[Stage] <= 0)
         ABORT_ERROR("TimeStepper {}: method {} has c[{}] = {}, but Omega "
                     "needs a positive time fraction for every stage after "
                     "the first",
                     Name, MethodName, Stage + 1, ProjFrac[Stage]);
   }

   LOG_INFO("TimeStepper {}: SUNDIALS explicit Runge-Kutta method {} with {} "
            "stages, order {}",
            Name, MethodName, NStages, Table->q);
}

//------------------------------------------------------------------------------
// Accessors
const std::string &SundialsStepper::getMethodName() const { return MethodName; }

int SundialsStepper::getNumStages() const { return Table->stages; }

int SundialsStepper::getOrder() const { return Table->q; }

N_Vector SundialsStepper::getStepEnd() const { return StepEnd; }

int SundialsStepper::getStageIndex() const {
   int Stage     = 0;
   int MaxStages = 0;
   checkArkode(ARKodeGetStageIndex(ArkMem, &Stage, &MaxStages),
               "ARKodeGetStageIndex");
   return Stage;
}

//------------------------------------------------------------------------------
// Allocates the stage tracer array once the tendencies are known
void SundialsStepper::finalizeInit() {

   if (!Tend)
      LOG_CRITICAL("Tendency not initialized");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!MeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");

   StageTracers = Array3DReal(
       "SundialsStageTracers", Tend->TracerTend.extent_int(0),
       Tend->TracerTend.extent_int(1), Tend->TracerTend.extent_int(2));
}

//------------------------------------------------------------------------------
// Returns the state vector made of three Omega arrays, creating it the first
// time those arrays are seen
N_Vector
SundialsStepper::getStateVector(const Array2DReal &Thickness,
                                const Array2DReal &Velocity,
                                const Array3DReal &WeightedTracers) const {

   const std::array<const Real *, 3> Key = {Thickness.data(), Velocity.data(),
                                            WeightedTracers.data()};
   auto Found                            = StateVectors.find(Key);
   if (Found != StateVectors.end())
      return Found->second;

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   // packTendencies relies on the pieces matching the tendency arrays
   OMEGA_REQUIRE(Thickness.size() == Tend->PseudoThicknessTend.size() &&
                     Velocity.size() == Tend->NormalVelocityTend.size() &&
                     WeightedTracers.size() == Tend->TracerTend.size(),
                 "TimeStepper {}: state arrays do not match the tendency "
                 "arrays",
                 Name);

   auto arrayVector = [&](const auto &Arr) {
      auto &Vec = ArrayVectors[{Arr.data(), Arr.size()}];
      if (!Vec)
         Vec = wrapArray(Arr, Ctx);
      return Vec->get();
   };
   N_Vector Pieces[3] = {arrayVector(Thickness), arrayVector(Velocity),
                         arrayVector(WeightedTracers)};

   N_Vector Vec = N_VMake_MPIManyVector(MeshHalo->getComm(), 3, Pieces, Ctx);
   OMEGA_REQUIRE(Vec, "TimeStepper {}: N_VMake_MPIManyVector failed", Name);
   StateVectors.emplace(Key, Vec);
   return Vec;
#else
   ABORT_ERROR("TimeStepper {}: SUNDIALS precision does not match Real", Name);
   return nullptr;
#endif
}

//------------------------------------------------------------------------------
// Creates the ERKStep integrator on the first step, with Template giving the
// layout of the state vector
void SundialsStepper::createIntegrator(N_Vector Template) const {

   ArkMem = ERKStepCreate(stepperRHS, 0.0, Template, Ctx);
   OMEGA_REQUIRE(ArkMem, "TimeStepper {}: ERKStepCreate failed", Name);

   checkArkode(ERKStepSetTable(ArkMem, Table), "ERKStepSetTable");
   checkArkode(ARKodeSetUserData(ArkMem, const_cast<SundialsStepper *>(this)),
               "ARKodeSetUserData");
   checkArkode(ARKodeSetPostprocessStageFn(ArkMem, stepperStage),
               "ARKodeSetPostprocessStageFn");
   // Taking one step per call needs no interpolation, which would add work
   // vectors
   checkArkode(ARKodeSetInterpolantType(ArkMem, ARK_INTERP_NONE),
               "ARKodeSetInterpolantType");
}

//------------------------------------------------------------------------------
// Prepares an intermediate stage, which ARKODE has formed in the next time
// level, as RungeKutta4Stepper prepares its provisional state
void SundialsStepper::prepareStage(int Stage) const {

   const int NextLevel         = 1;
   const TimeInstant StageTime = StepStart + StageFrac[Stage] * TimeStep;

   // As after Omega's own provisional update: no velocity outside the active
   // layers, and the prescribed state at the stage time
   VCoord->applyEdgeLayerMask(StepState->getNormalVelocity(NextLevel),
                              Mesh->NEdgesAll);
   prescribeState(StepState, NextLevel, StepState, NextLevel, StageTime);

   // Tracer concentrations from the thickness-weighted tracers
   Array2DReal Thickness       = StepState->getPseudoThickness(NextLevel);
   Array3DReal WeightedTracers = Tracers::getAll(NextLevel);
   OMEGA_SCOPE(LocStageTracers, StageTracers);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   const int NTracers = StageTracers.extent_int(0);

   parallelForOuter(
       "sundialsStageTracers", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 LocStageTracers(L, ICell, K) =
                     WeightedTracers(L, ICell, K) / Thickness(ICell, K);
              });
       });

   const MPI_Comm Comm = MeshHalo->getComm();
   Pacer::timingBarrier("Sundials:haloExchStageBarrier", 3, Comm);
   Pacer::start("Sundials:haloExchStage", 3);
   StepState->exchangeHalo(NextLevel);
   MeshHalo->exchangeFullArrayHalo(StageTracers, OnCell);
   Pacer::stop("Sundials:haloExchStage", 3);
}

//------------------------------------------------------------------------------
// Computes the tendencies at a stage and copies them to YDot
void SundialsStepper::computeStageTendencies(int Stage, N_Vector YDot) const {

   const int CurLevel          = 0;
   const int NextLevel         = 1;
   const TimeInstant StageTime = StepStart + StageFrac[Stage] * TimeStep;
   const TimeInterval ProjDt   = ProjFrac[Stage] * TimeStep;

   if (Stage == 0) {
      // ARKODE's copy of the start state equals the current time level
      Tend->computeAllTendencies(StepState, AuxState, Tracers::getAll(CurLevel),
                                 CurLevel, CurLevel, CurLevel, StageTime,
                                 ProjDt);
   } else {
      Tend->computeAllTendencies(StepState, AuxState, StageTracers, NextLevel,
                                 NextLevel, NextLevel, StageTime, ProjDt);
   }

   packTendencies(YDot);
}

//------------------------------------------------------------------------------
// Copies the tendencies into YDot. Inactive layers and padding get zero, so
// ARKODE leaves them as they were at the start of the step, while the active
// range is the one Omega's own updates write.
void SundialsStepper::packTendencies(N_Vector YDot) const {

#ifdef OMEGA_SUNDIALS_MATCHES_REAL
   Pacer::start("Sundials:packTendencies", 3);

   const Array2DReal ThickTend  = Tend->PseudoThicknessTend;
   const Array2DReal VelTend    = Tend->NormalVelocityTend;
   const Array3DReal TracerTend = Tend->TracerTend;

   // The pieces of YDot have the sizes and layouts of the tendency arrays
   const Array2DReal ThickDot(pieceData(YDot, 0), ThickTend.extent(0),
                              ThickTend.extent(1));
   const Array2DReal VelDot(pieceData(YDot, 1), VelTend.extent(0),
                            VelTend.extent(1));
   const Array3DReal TracerDot(pieceData(YDot, 2), TracerTend.extent(0),
                               TracerTend.extent(1), TracerTend.extent(2));

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   const int NCellsAll = Mesh->NCellsAll;
   const int NEdgesAll = Mesh->NEdgesAll;

   parallelFor(
       "sundialsPackThickness",
       {ThickTend.extent_int(0), ThickTend.extent_int(1)},
       KOKKOS_LAMBDA(int ICell, int K) {
          const bool Active = ICell < NCellsAll && K >= MinLayerCell(ICell) &&
                              K <= MaxLayerCell(ICell);
          ThickDot(ICell, K) = Active ? ThickTend(ICell, K) : 0;
       });

   parallelFor(
       "sundialsPackVelocity", {VelTend.extent_int(0), VelTend.extent_int(1)},
       KOKKOS_LAMBDA(int IEdge, int K) {
          const bool Active = IEdge < NEdgesAll &&
                              K >= MinLayerEdgeBot(IEdge) &&
                              K <= MaxLayerEdgeTop(IEdge);
          VelDot(IEdge, K) = Active ? VelTend(IEdge, K) : 0;
       });

   parallelFor(
       "sundialsPackTracers",
       {TracerTend.extent_int(0), TracerTend.extent_int(1),
        TracerTend.extent_int(2)},
       KOKKOS_LAMBDA(int L, int ICell, int K) {
          const bool Active = ICell < NCellsAll && K >= MinLayerCell(ICell) &&
                              K <= MaxLayerCell(ICell);
          TracerDot(L, ICell, K) = Active ? TracerTend(L, ICell, K) : 0;
       });

   Pacer::stop("Sundials:packTendencies", 3);
#else
   ABORT_ERROR("TimeStepper {}: SUNDIALS precision does not match Real", Name);
#endif
}

//------------------------------------------------------------------------------
// Advance the state by one step of the chosen Runge-Kutta method
void SundialsStepper::doStep(OceanState *State,   // model state
                             TimeInstant &SimTime // current simulation time
) const {

   const MPI_Comm Comm = MeshHalo->getComm();

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
   Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

   // Work done before the first stage, as in RungeKutta4Stepper. ARKODE
   // advances thickness-weighted tracers, which are kept in the next tracer
   // time level.
   weightTracers(NextTracerArray, CurTracerArray, State, CurLevel);
   prescribeState(State, CurLevel, State, CurLevel, SimTime);
   updateKPPFields(State, CurLevel, CurLevel, CurLevel);

   if (!Ctx) {
      if (SUNContext_Create(Comm, &Ctx) != SUN_SUCCESS)
         ABORT_ERROR("TimeStepper {}: SUNContext_Create failed", Name);
   }
   N_Vector StartVec =
       getStateVector(State->getPseudoThickness(CurLevel),
                      State->getNormalVelocity(CurLevel), NextTracerArray);
   N_Vector EndVec =
       getStateVector(State->getPseudoThickness(NextLevel),
                      State->getNormalVelocity(NextLevel), NextTracerArray);
   if (!ArkMem)
      createIntegrator(StartVec);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   // ARKODE copies the start state, so the stages can overwrite the
   // weighted tracers it shares with the end state
   StepState = State;
   StepStart = SimTime;
   StepEnd   = EndVec;
   checkArkode(ARKodeReset(ArkMem, 0.0, StartVec), "ARKodeReset");
   checkArkode(ARKodeSetFixedStep(ArkMem, TimeStepSeconds),
               "ARKodeSetFixedStep");

   Pacer::start("Sundials:arkodeEvolve", 2);
   sunrealtype TReached = 0.0;
   checkArkode(
       ARKodeEvolve(ArkMem, TimeStepSeconds, EndVec, &TReached, ARK_ONE_STEP),
       "ARKodeEvolve");
   Pacer::stop("Sundials:arkodeEvolve", 2);

   StepState = nullptr;
   StepEnd   = nullptr;

   // As after Omega's own final update
   VCoord->applyEdgeLayerMask(State->getNormalVelocity(NextLevel),
                              Mesh->NEdgesAll);
   finalizeTracersUpdate(NextTracerArray, State, NextLevel);

   // Update time levels (New -> Old) of prognostic variables with halo
   // exchanges
   Pacer::timingBarrier("Sundials:haloExchBarrier", 3, Comm);
   Pacer::start("Sundials:haloExch", 3);
   State->updateTimeLevels();
   Tracers::updateTimeLevels();
   Pacer::stop("Sundials:haloExch", 3);

   applyImplicitVerticalMixing(State, CurLevel, CurLevel, CurLevel, "Sundials");

   validateOceanState(State, AuxState, VertCoord::getDefault(), CurLevel);

   // Advance the clock and update the simulation time
   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
   ++StepCount;
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//

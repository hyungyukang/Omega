//===-- Test driver for the SUNDIALS library in OMEGA -----------*- C++ -*-===/
//
/// \file
/// \brief Test driver for calling SUNDIALS from OMEGA
///
/// This driver checks that the SUNDIALS library built under external/sundials
/// is compiled, linked and callable from OMEGA, and tests the pieces the
/// Sundials time stepper is built from. It
///  - loads the SUNDIALS table ARKODE_KUTTA_RK4a_4_4 by name and checks that
///    it has exactly the coefficients of OMEGA's RungeKutta4Stepper
///  - integrates the linear system dy_i/dt = Lambda_i y_i with that table on
///    an OMEGA host array wrapped, without copying, in a serial N_Vector
///  - integrates the same system on two OMEGA arrays in OMEGA's memory space,
///    wrapped as SUNDIALS Kokkos vectors and joined in an MPI ManyVector, with
///    the right-hand side computed by an OMEGA kernel. As in the Sundials
///    time stepper, a second pair of arrays plays the next time level, and
///    each step resets ARKODE to one level and writes the result to the other.
///  - checks the fused linear combination added to the Kokkos vectors against
///    SUNDIALS' generic definition
/// After NSteps steps each component must equal y0 * R(Lambda_i * h)^NSteps,
/// where
///    R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24
/// is the amplification factor of classical RK4. That holds only if ARKODE
/// takes exactly the four-stage RK4 step, so the check is to round-off.
///
//
//===-----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "SundialsVector.h"
#include "mpi.h"

#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_mpimanyvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_core.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

using namespace OMEGA;

// Problem and step parameters shared by the host and device integrations.
// On the device the system is split into two pieces of NVals1 and NVals2
// values.
constexpr int NVals      = 8;
constexpr int NVals1     = 5;
constexpr int NVals2     = NVals - NVals1;
constexpr int NSteps     = 20;
constexpr int NRK4Stages = 4;
constexpr sunrealtype T0 = 0.0;
constexpr sunrealtype H  = 0.1;
constexpr sunrealtype Y0 = 1.0;

// Right-hand side dy_i/dt = Lambda_i * y_i on the host
static int linearRHSHost(sunrealtype T, N_Vector Y, N_Vector YDot,
                         void *UserData) {
   const HostArray1DReal &Lambda = *static_cast<HostArray1DReal *>(UserData);
   const sunrealtype *YData      = N_VGetArrayPointer(Y);
   sunrealtype *YDotData         = N_VGetArrayPointer(YDot);
   for (int I = 0; I < NVals; ++I) {
      YDotData[I] = Lambda(I) * YData[I];
   }
   return 0;
}

// Right-hand side dy_i/dt = Lambda_i * y_i on each piece of a ManyVector,
// computed by an Omega kernel on unmanaged arrays over the SUNDIALS data
int linearRHSDevice(sunrealtype T, N_Vector Y, N_Vector YDot, void *UserData) {
   const Array1DReal *Lambda = static_cast<Array1DReal *>(UserData);
   for (int Piece = 0; Piece < 2; ++Piece) {
      N_Vector YPiece    = N_VGetSubvector_MPIManyVector(Y, Piece);
      N_Vector YDotPiece = N_VGetSubvector_MPIManyVector(YDot, Piece);
      const int Length   = N_VGetLength(YPiece);
      const Array1DReal LambdaPiece = Lambda[Piece];
      const Array1DReal YArr(N_VGetDeviceArrayPointer(YPiece), Length);
      const Array1DReal YDotArr(N_VGetDeviceArrayPointer(YDotPiece), Length);
      parallelFor(
          {Length},
          KOKKOS_LAMBDA(int I) { YDotArr(I) = LambdaPiece(I) * YArr(I); });
   }
   return 0;
}

// Creates an ERKStep integrator for the built-in classical RK4 with a fixed
// step H, starting at T0 from the initial condition in Y
static void *createRK4(ARKRhsFn RHS, N_Vector Y, void *UserData,
                       SUNContext Ctx) {
   void *ArkMem = ERKStepCreate(RHS, T0, Y, Ctx);
   if (!ArkMem)
      ABORT_ERROR("SundialsTest: FAIL ERKStepCreate");
   if (ERKStepSetTableName(ArkMem, "ARKODE_KUTTA_RK4a_4_4") != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ERKStepSetTableName");
   if (ARKodeSetUserData(ArkMem, UserData) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ARKodeSetUserData");
   if (ARKodeSetFixedStep(ArkMem, H) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ARKodeSetFixedStep");
   return ArkMem;
}

// Aborts unless ARKODE took NSteps steps with four RHS evaluations each
static void checkCounters(const std::string &Label, void *ArkMem) {
   long int NStepsTaken = 0;
   long int NRhsEvals   = 0;
   ARKodeGetNumSteps(ArkMem, &NStepsTaken);
   ARKodeGetNumRhsEvals(ArkMem, 0, &NRhsEvals);
   LOG_INFO("SundialsTest: {} steps = {}, RHS evaluations = {}", Label,
            NStepsTaken, NRhsEvals);

   if (NStepsTaken != NSteps || NRhsEvals != NRK4Stages * NSteps)
      ABORT_ERROR("SundialsTest: FAIL {} took {} steps with {} RHS "
                  "evaluations, expected {} with {}",
                  Label, NStepsTaken, NRhsEvals, NSteps, NRK4Stages * NSteps);
}

// Aborts unless Result is y0 * R(Lambda_i * H)^NSteps, the exact result of
// NSteps classical RK4 steps, to within round-off
template <class HostArray>
static void checkRK4Result(const std::string &Label, const HostArray &Result,
                           const HostArray1DReal &Lambda) {
   const sunrealtype Tol = 1000 * std::numeric_limits<sunrealtype>::epsilon();
   sunrealtype MaxDiff   = 0.0;
   for (int I = 0; I < Result.extent_int(0); ++I) {
      const sunrealtype Z = Lambda(I) * H;
      const sunrealtype Amp =
          1 + Z + Z * Z / 2 + Z * Z * Z / 6 + Z * Z * Z * Z / 24;
      const sunrealtype YRK4 = Y0 * std::pow(Amp, NSteps);
      MaxDiff = std::max(MaxDiff, std::abs(Result(I) - YRK4) / YRK4);
   }
   LOG_INFO("SundialsTest: {} max rel diff from RK4 amplification = {:.3e}",
            Label, MaxDiff);
   if (MaxDiff > Tol)
      ABORT_ERROR("SundialsTest: FAIL {} result does not match RK4 "
                  "({:.3e} > {:.3e})",
                  Label, MaxDiff, Tol);
}

// Checks the built-in classical RK4 table, loaded by name, against the
// coefficients of RungeKutta4Stepper written as a full Butcher tableau
static void checkRK4Table() {
   const sunrealtype C[NRK4Stages]             = {0.0, 0.5, 0.5, 1.0};
   const sunrealtype A[NRK4Stages][NRK4Stages] = {{0.0, 0.0, 0.0, 0.0},
                                                  {0.5, 0.0, 0.0, 0.0},
                                                  {0.0, 0.5, 0.0, 0.0},
                                                  {0.0, 0.0, 1.0, 0.0}};
   const sunrealtype B[NRK4Stages] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0,
                                      1.0 / 6.0};

   ARKodeButcherTable RK4Table =
       ARKodeButcherTable_LoadERKByName("ARKODE_KUTTA_RK4a_4_4");
   if (!RK4Table)
      ABORT_ERROR("SundialsTest: FAIL ARKodeButcherTable_LoadERKByName");

   int MethodOrder   = 0;
   int EmbeddedOrder = 0;
   // A positive return only says there is no embedding, which RK4 lacks
   if (ARKodeButcherTable_CheckOrder(RK4Table, &MethodOrder, &EmbeddedOrder,
                                     nullptr) < 0 ||
       MethodOrder != 4)
      ABORT_ERROR("SundialsTest: FAIL RK4 table order is {}, not 4",
                  MethodOrder);
   if (RK4Table->stages != NRK4Stages)
      ABORT_ERROR("SundialsTest: FAIL RK4 table has {} stages, not {}",
                  RK4Table->stages, NRK4Stages);

   int NDiffs = 0;
   for (int I = 0; I < NRK4Stages; ++I) {
      if (RK4Table->b[I] != B[I] || RK4Table->c[I] != C[I])
         ++NDiffs;
      for (int J = 0; J < NRK4Stages; ++J) {
         if (RK4Table->A[I][J] != A[I][J])
            ++NDiffs;
      }
   }
   if (NDiffs != 0)
      ABORT_ERROR("SundialsTest: FAIL {} coefficients of "
                  "ARKODE_KUTTA_RK4a_4_4 differ from RungeKutta4Stepper",
                  NDiffs);
   LOG_INFO("SundialsTest: ARKODE_KUTTA_RK4a_4_4 matches the "
            "RungeKutta4Stepper coefficients");
   ARKodeButcherTable_Free(RK4Table);
}

// Host: the state is an Omega host array that ARKODE updates in place
static void testHost(const HostArray1DReal &Lambda, SUNContext Ctx) {
   HostArray1DReal State("State", NVals);
   Kokkos::deep_copy(State, Y0);

   N_Vector Y = N_VMake_Serial(NVals, State.data(), Ctx);
   if (!Y)
      ABORT_ERROR("SundialsTest: FAIL N_VMake_Serial");

   void *ArkMem =
       createRK4(linearRHSHost, Y, const_cast<HostArray1DReal *>(&Lambda), Ctx);
   sunrealtype TCur = T0;
   for (int Step = 0; Step < NSteps; ++Step) {
      if (ARKodeEvolve(ArkMem, T0 + NSteps * H, Y, &TCur, ARK_ONE_STEP) !=
          ARK_SUCCESS)
         ABORT_ERROR("SundialsTest: FAIL host ARKodeEvolve at step {}", Step);
   }
   checkCounters("host", ArkMem);
   checkRK4Result("host", State, Lambda);

   ARKodeFree(&ArkMem);
   N_VDestroy(Y);
}

// Omega's memory space: two time levels, each two arrays wrapped as Kokkos
// vectors and joined in an MPI ManyVector
static void testKokkosVectors(const HostArray1DReal &Lambda, SUNContext Ctx,
                              MPI_Comm Comm) {

   // Lambda split into the two pieces
   HostArray1DReal Lambda1("Lambda1", NVals1);
   HostArray1DReal Lambda2("Lambda2", NVals2);
   for (int I = 0; I < NVals; ++I) {
      if (I < NVals1)
         Lambda1(I) = Lambda(I);
      else
         Lambda2(I - NVals1) = Lambda(I);
   }
   Array1DReal LambdaDev[2] = {createDeviceMirrorCopy(Lambda1),
                               createDeviceMirrorCopy(Lambda2)};

   // Level 0 starts at Y0; level 1 is the other time level
   Array1DReal Level[2][2] = {{Array1DReal("Level0Piece1", NVals1),
                               Array1DReal("Level0Piece2", NVals2)},
                              {Array1DReal("Level1Piece1", NVals1),
                               Array1DReal("Level1Piece2", NVals2)}};
   deepCopy(Level[0][0], Y0);
   deepCopy(Level[0][1], Y0);

   std::unique_ptr<SundialsVector> Wrapped[2][2];
   N_Vector Y[2];
   for (int L = 0; L < 2; ++L) {
      Wrapped[L][0] = wrapArray(Level[L][0], Ctx);
      Wrapped[L][1] = wrapArray(Level[L][1], Ctx);
      if (N_VGetDeviceArrayPointer(*Wrapped[L][0]) != Level[L][0].data() ||
          N_VGetDeviceArrayPointer(*Wrapped[L][1]) != Level[L][1].data())
         ABORT_ERROR("SundialsTest: FAIL wrapArray does not share memory");
      N_Vector Pieces[2] = {*Wrapped[L][0], *Wrapped[L][1]};
      Y[L]               = N_VMake_MPIManyVector(Comm, 2, Pieces, Ctx);
      if (!Y[L])
         ABORT_ERROR("SundialsTest: FAIL N_VMake_MPIManyVector");
   }

   void *ArkMem = createRK4(linearRHSDevice, Y[0], LambdaDev, Ctx);
   if (ARKodeSetInterpolantType(ArkMem, ARK_INTERP_NONE) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ARKodeSetInterpolantType");

   // Each step starts from one level and writes the other, as the time levels
   // of an Omega state alternate
   for (int Step = 0; Step < NSteps; ++Step) {
      const int Cur    = Step % 2;
      const int Next   = 1 - Cur;
      sunrealtype TCur = 0;
      if (ARKodeReset(ArkMem, 0.0, Y[Cur]) != ARK_SUCCESS ||
          ARKodeEvolve(ArkMem, H, Y[Next], &TCur, ARK_ONE_STEP) != ARK_SUCCESS)
         ABORT_ERROR("SundialsTest: FAIL Kokkos ARKodeEvolve at step {}", Step);
   }
   Kokkos::fence();
   checkCounters("Kokkos", ArkMem);

   const int Last = NSteps % 2;
   HostArray1DReal Result("Result", NVals);
   auto Piece1 = createHostMirrorCopy(Level[Last][0]);
   auto Piece2 = createHostMirrorCopy(Level[Last][1]);
   for (int I = 0; I < NVals; ++I) {
      Result(I) = I < NVals1 ? Piece1(I) : Piece2(I - NVals1);
   }
   checkRK4Result("Kokkos", Result, Lambda);

   ARKodeFree(&ArkMem);
   N_VDestroy(Y[0]);
   N_VDestroy(Y[1]);
}

// Checks the fused linear combination against SUNDIALS' generic definition,
// on vectors cloned from a wrapped array, which must keep the fused operation
void testFusedLinearCombination(SUNContext Ctx) {
   constexpr int Length = 1000;
   constexpr int NTerms = 5;

   Array1DReal Base("Base", Length);
   std::unique_ptr<SundialsVector> Wrapped = wrapArray(Base, Ctx);

   N_Vector X[NTerms];
   for (int J = 0; J < NTerms; ++J) {
      X[J] = N_VClone(*Wrapped);
      if (X[J]->ops->nvlinearcombination != fusedLinearCombination)
         ABORT_ERROR("SundialsTest: FAIL a clone lost the fused operation");
      const Array1DReal XArr(N_VGetDeviceArrayPointer(X[J]), Length);
      parallelFor(
          {Length}, KOKKOS_LAMBDA(int I) {
             XArr(I) = Kokkos::sin(0.01 * (I + 1) * (J + 1)) + J;
          });
   }
   N_Vector ZFused   = N_VClone(*Wrapped);
   N_Vector ZGeneric = N_VClone(*Wrapped);

   // Includes a zero coefficient, which the fused version skips
   sunrealtype C[NTerms] = {0.5, 0.0, -1.25, 2.0, 1.0};
   N_VLinearCombination(NTerms, C, X, ZFused);
   N_VScale(C[0], X[0], ZGeneric);
   for (int J = 1; J < NTerms; ++J) {
      N_VLinearSum(C[J], X[J], 1.0, ZGeneric, ZGeneric);
   }

   // The result may also overwrite one of the terms
   N_Vector Aliased = N_VClone(*Wrapped);
   N_VScale(1.0, X[NTerms - 1], Aliased);
   N_Vector XAliased[NTerms];
   std::copy(X, X + NTerms, XAliased);
   XAliased[NTerms - 1] = Aliased;
   N_VLinearCombination(NTerms, C, XAliased, Aliased);

   Kokkos::fence();
   const Array1DReal FusedArr(N_VGetDeviceArrayPointer(ZFused), Length);
   const Array1DReal GenericArr(N_VGetDeviceArrayPointer(ZGeneric), Length);
   const Array1DReal AliasedArr(N_VGetDeviceArrayPointer(Aliased), Length);
   auto FusedHost   = createHostMirrorCopy(FusedArr);
   auto GenericHost = createHostMirrorCopy(GenericArr);
   auto AliasedHost = createHostMirrorCopy(AliasedArr);

   // The two sum the terms in different orders
   const sunrealtype Tol = 100 * std::numeric_limits<sunrealtype>::epsilon();
   sunrealtype MaxDiff   = 0;
   for (int I = 0; I < Length; ++I) {
      MaxDiff = std::max(MaxDiff, std::abs(FusedHost(I) - GenericHost(I)));
      MaxDiff = std::max(MaxDiff, std::abs(AliasedHost(I) - FusedHost(I)));
   }
   LOG_INFO("SundialsTest: fused linear combination max diff from generic = "
            "{:.3e}",
            MaxDiff);
   if (MaxDiff > Tol)
      ABORT_ERROR("SundialsTest: FAIL fused linear combination differs from "
                  "the generic one ({:.3e} > {:.3e})",
                  MaxDiff, Tol);

   for (int J = 0; J < NTerms; ++J) {
      N_VDestroy(X[J]);
   }
   N_VDestroy(ZFused);
   N_VDestroy(ZGeneric);
   N_VDestroy(Aliased);
}

int main(int argc, char *argv[]) {

   // initialize environments
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);
   LOG_INFO("------ SUNDIALS Unit Tests ------");

   {
      // SUNDIALS was built with the same precision as Omega's Real, so an
      // Omega array can be handed to SUNDIALS without conversion
      static_assert(std::is_same_v<sunrealtype, Real>,
                    "SUNDIALS precision does not match Omega Real");

      char Version[32];
      if (SUNDIALSGetVersion(Version, static_cast<int>(sizeof(Version))) !=
          SUN_SUCCESS)
         ABORT_ERROR("SundialsTest: FAIL SUNDIALSGetVersion");
      LOG_INFO("SundialsTest: SUNDIALS version {}", Version);

      // The context carries the Omega communicator
      SUNContext Ctx = nullptr;
      if (SUNContext_Create(DefEnv->getComm(), &Ctx) != SUN_SUCCESS)
         ABORT_ERROR("SundialsTest: FAIL SUNContext_Create");

      checkRK4Table();

      HostArray1DReal Lambda("Lambda", NVals);
      for (int I = 0; I < NVals; ++I) {
         Lambda(I) = -0.25 * (I + 1);
      }

      testHost(Lambda, Ctx);
      testKokkosVectors(Lambda, Ctx, DefEnv->getComm());
      testFusedLinearCombination(Ctx);

      SUNContext_Free(&Ctx);
   }

   LOG_INFO("------ SUNDIALS Unit Tests Successful ------");
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, return successfully

} // end of main
//===-----------------------------------------------------------------------===/

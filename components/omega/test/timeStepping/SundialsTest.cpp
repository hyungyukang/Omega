//===-- Test driver for the SUNDIALS library in OMEGA -----------*- C++ -*-===/
//
/// \file
/// \brief Test driver for calling SUNDIALS from OMEGA
///
/// This driver checks that the SUNDIALS library built under external/sundials
/// is compiled, linked and callable from OMEGA, using ARKODE's ERKStep with
/// the classical fourth-order Runge-Kutta method. It
///  - checks that the built-in SUNDIALS table ARKODE_KUTTA_RK4a_4_4 has
///    exactly the coefficients of OMEGA's RungeKutta4Stepper
///  - integrates the linear system dy_i/dt = Lambda_i y_i with that table on
///    an OMEGA host array wrapped, without copying, in a serial N_Vector
///  - on HIP builds, integrates the same system on OMEGA device arrays wrapped
///    in a HIP N_Vector, with the right-hand side computed by an OMEGA kernel
///    and the SUNDIALS kernels run on the Kokkos stream. After ERKStep has
///    read the initial condition, N_VSetDeviceArrayPointer points the vector
///    at a second array, as for the next time level, so ARKODE must write its
///    result there and leave the first array alone.
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
#include "mpi.h"

#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_core.h>
#ifdef OMEGA_ENABLE_HIP
#include <nvector/nvector_hip.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

using namespace OMEGA;

// Problem and step parameters shared by the host and device integrations
constexpr int NVals      = 8;
constexpr int NSteps     = 20;
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

#ifdef OMEGA_ENABLE_HIP
// Right-hand side dy_i/dt = Lambda_i * y_i computed on the device by an Omega
// kernel. The N_Vector data are wrapped in unmanaged Omega arrays, so nothing
// is copied.
static int linearRHSDevice(sunrealtype T, N_Vector Y, N_Vector YDot,
                           void *UserData) {
   const Array1DReal Lambda = *static_cast<Array1DReal *>(UserData);
   const Array1DReal YArr(N_VGetDeviceArrayPointer(Y), NVals);
   const Array1DReal YDotArr(N_VGetDeviceArrayPointer(YDot), NVals);
   parallelFor(
       {NVals}, KOKKOS_LAMBDA(int I) { YDotArr(I) = Lambda(I) * YArr(I); });
   return 0;
}
#endif

// Creates an ERKStep integrator for the built-in classical RK4 with a fixed
// step H, starting at T0 from the initial condition in Y
static void *createRK4(ARKRhsFn RHS, N_Vector Y, void *UserData,
                       SUNContext Ctx) {
   void *ArkMem = ERKStepCreate(RHS, T0, Y, Ctx);
   if (!ArkMem)
      ABORT_ERROR("SundialsTest: FAIL ERKStepCreate");
   if (ERKStepSetTableNum(ArkMem, ARKODE_KUTTA_RK4a_4_4) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ERKStepSetTableNum");
   if (ARKodeSetUserData(ArkMem, UserData) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ARKodeSetUserData");
   if (ARKodeSetFixedStep(ArkMem, H) != ARK_SUCCESS)
      ABORT_ERROR("SundialsTest: FAIL ARKodeSetFixedStep");
   return ArkMem;
}

// Takes NSteps steps, leaving the result in Y. Each ARKodeEvolve call takes
// one step, as an Omega doStep would, so ARKODE never interpolates to an
// output time.
static void takeRK4Steps(const std::string &Label, void *ArkMem, N_Vector Y) {
   sunrealtype TCur = T0;
   for (int Step = 0; Step < NSteps; ++Step) {
      if (ARKodeEvolve(ArkMem, T0 + NSteps * H, Y, &TCur, ARK_ONE_STEP) !=
          ARK_SUCCESS)
         ABORT_ERROR("SundialsTest: FAIL {} ARKodeEvolve at step {}", Label,
                     Step);
   }

   long int NStepsTaken = 0;
   long int NRhsEvals   = 0;
   ARKodeGetNumSteps(ArkMem, &NStepsTaken);
   ARKodeGetNumRhsEvals(ArkMem, 0, &NRhsEvals);
   LOG_INFO("SundialsTest: {} t = {}, steps = {}, RHS evaluations = {}", Label,
            TCur, NStepsTaken, NRhsEvals);

   if (NStepsTaken != NSteps)
      ABORT_ERROR("SundialsTest: FAIL {} took {} steps, expected {}", Label,
                  NStepsTaken, NSteps);
}

// Aborts unless Result is y0 * R(Lambda_i * H)^NSteps, the exact result of
// NSteps classical RK4 steps, to within round-off
template <class HostArray>
static void checkRK4Result(const std::string &Label, const HostArray &Result,
                           const HostArray1DReal &Lambda) {
   const sunrealtype Tol = 1000 * std::numeric_limits<sunrealtype>::epsilon();
   sunrealtype MaxDiff   = 0.0;
   for (int I = 0; I < NVals; ++I) {
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

      // Compare the built-in classical RK4 table with the coefficients of
      // RungeKutta4Stepper, written here as a full Butcher tableau
      constexpr int NStages                 = 4;
      const sunrealtype C[NStages]          = {0.0, 0.5, 0.5, 1.0};
      const sunrealtype A[NStages][NStages] = {{0.0, 0.0, 0.0, 0.0},
                                               {0.5, 0.0, 0.0, 0.0},
                                               {0.0, 0.5, 0.0, 0.0},
                                               {0.0, 0.0, 1.0, 0.0}};
      const sunrealtype B[NStages]          = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0,
                                               1.0 / 6.0};



// string version by name
      ARKodeButcherTable RK4Table =
          ARKodeButcherTable_LoadERK(ARKODE_KUTTA_RK4a_4_4);
      if (!RK4Table)
         ABORT_ERROR("SundialsTest: FAIL ARKodeButcherTable_LoadERK");

      int MethodOrder   = 0;
      int EmbeddedOrder = 0;
      // A positive return only says there is no embedding, which RK4 lacks
      if (ARKodeButcherTable_CheckOrder(RK4Table, &MethodOrder, &EmbeddedOrder,
                                        nullptr) < 0 ||
          MethodOrder != 4)
         ABORT_ERROR("SundialsTest: FAIL RK4 table order is {}, not 4",
                     MethodOrder);
      if (RK4Table->stages != NStages)
         ABORT_ERROR("SundialsTest: FAIL RK4 table has {} stages, not {}",
                     RK4Table->stages, NStages);

      int NDiffs = 0;
      for (int I = 0; I < NStages; ++I) {
         if (RK4Table->b[I] != B[I] || RK4Table->c[I] != C[I])
            ++NDiffs;
         for (int J = 0; J < NStages; ++J) {
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

      HostArray1DReal Lambda("Lambda", NVals);
      for (int I = 0; I < NVals; ++I) {
         Lambda(I) = -0.25 * (I + 1);
      }

      // Host: the state is an Omega host array that ARKODE updates in place
      {
         HostArray1DReal State("State", NVals);
         Kokkos::deep_copy(State, Y0);

         N_Vector Y = N_VMake_Serial(NVals, State.data(), Ctx);
         if (!Y)
            ABORT_ERROR("SundialsTest: FAIL N_VMake_Serial");

         void *ArkMem = createRK4(linearRHSHost, Y, &Lambda, Ctx);
         takeRK4Steps("host", ArkMem, Y);
         checkRK4Result("host", State, Lambda);

         ARKodeFree(&ArkMem);
         N_VDestroy(Y);
      }

#ifdef OMEGA_ENABLE_HIP

/// kokkos_vector
/// The NVECTOR_MPIMANYVECTOR Module


      // Device: two Omega device arrays play the current and next time levels
      {
         Array1DReal LambdaDev = createDeviceMirrorCopy(Lambda);
         Array1DReal StateCur("StateCur", NVals);
         Array1DReal StateNext("StateNext", NVals);
         Kokkos::deep_copy(StateCur, Y0);

         // N_VMake_Hip wraps both pointers without taking ownership, so the
         // Omega arrays keep their memory. Its host array is used only by
         // N_VCopyFromDevice_Hip and N_VCopyToDevice_Hip.
         HostArray1DReal StateHost("StateHost", NVals);
         N_Vector Y =
             N_VMake_Hip(NVals, StateHost.data(), StateCur.data(), Ctx);
         if (!Y)
            ABORT_ERROR("SundialsTest: FAIL N_VMake_Hip");

         // Run the SUNDIALS kernels on the stream Omega's kernels use, so each
         // waits for the other without explicit synchronization. Vectors that
         // ARKODE clones from Y inherit these policies.
         const hipStream_t Stream = ExecSpace().hip_stream();
         SUNHipThreadDirectExecPolicy StreamPolicy(256, Stream);
         SUNHipBlockReduceExecPolicy ReducePolicy(256, 0, Stream);
         if (N_VSetKernelExecPolicy_Hip(Y, &StreamPolicy, &ReducePolicy) !=
             SUN_SUCCESS)
            ABORT_ERROR("SundialsTest: FAIL N_VSetKernelExecPolicy_Hip");

         void *ArkMem = createRK4(linearRHSDevice, Y, &LambdaDev, Ctx);

         // ERKStepCreate has copied the initial condition out of StateCur;
         // the steps write their result to whichever array Y points at
         N_VSetDeviceArrayPointer(StateNext.data(), Y);
         if (N_VGetDeviceArrayPointer(Y) != StateNext.data())
            ABORT_ERROR("SundialsTest: FAIL N_VSetDeviceArrayPointer");

         takeRK4Steps("device", ArkMem, Y);
         Kokkos::fence();

         auto NextHost = createHostMirrorCopy(StateNext);
         auto CurHost  = createHostMirrorCopy(StateCur);
         checkRK4Result("device", NextHost, Lambda);

         int NChanged = 0;
         for (int I = 0; I < NVals; ++I) {
            if (CurHost(I) != Y0)
               ++NChanged;
         }
         if (NChanged != 0)
            ABORT_ERROR("SundialsTest: FAIL {} values of the array Y no "
                        "longer points at were overwritten",
                        NChanged);

         ARKodeFree(&ArkMem);
         N_VDestroy(Y);
      }
#endif

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

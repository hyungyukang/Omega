#ifndef OMEGA_SUNDIALSVECTOR_H
#define OMEGA_SUNDIALSVECTOR_H
//===-- SundialsVector.h - Omega arrays as SUNDIALS vectors -----*- C++ -*-===//
//
/// \file
/// \brief Wraps Omega arrays as SUNDIALS vectors without copying them
///
/// SUNDIALS integrators work on N_Vector objects. Here an Omega array is
/// wrapped as a SUNDIALS Kokkos vector, which shares the array's memory and
/// runs its kernels in Omega's execution space, so SUNDIALS and Omega kernels
/// are ordered without extra synchronization. The Kokkos vector has no fused
/// linear combination, which ARKODE uses to form each stage and the step
/// result, so one is added here.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Error.h"
#include "OmegaKokkos.h"

#include <memory>
#include <nvector/nvector_kokkos.hpp>
#include <type_traits>

namespace OMEGA {

/// SUNDIALS Kokkos vector in Omega's memory space. Vectors that SUNDIALS
/// clones from it allocate their own memory in that space.
using SundialsVector = sundials::kokkos::Vector<ExecSpace, MemSpace>;

/// Fused linear combination Z = sum_j C[j] * X[j] in one kernel, for
/// SundialsVector arguments. Terms with zero coefficients are skipped, and
/// the last term is added first: ARKODE's explicit stepper lists the step's
/// start state last, so the sum is then formed in the same order as Omega's
/// own Runge-Kutta updates.
SUNErrCode fusedLinearCombination(int NVec, sunrealtype *C, N_Vector *X,
                                  N_Vector Z);

/// Gives a SundialsVector the fused linear combination. Vectors cloned from
/// it, as ARKODE does for its work vectors, get it too.
void enableFusedOps(N_Vector Vec);

/// Wraps a contiguous Omega array of any rank as a SundialsVector that
/// shares the array's memory. The array must outlive the vector.
template <class ArrayType>
std::unique_ptr<SundialsVector> wrapArray(const ArrayType &Arr,
                                          SUNContext Ctx) {
   static_assert(
       std::is_same_v<typename ArrayType::non_const_value_type, sunrealtype>,
       "SUNDIALS precision does not match the Omega array");
   static_assert(std::is_same_v<typename ArrayType::memory_space, MemSpace>,
                 "Only arrays in Omega's memory space can be wrapped");
   OMEGA_REQUIRE(Arr.span_is_contiguous() && Arr.span() == Arr.size(),
                 "wrapArray: array {} is not contiguous", Arr.label());

   // A view built from a pointer does not own the memory
   SundialsVector::view_type View(Arr.data(), Arr.size());

   // The Kokkos vector would otherwise allocate a host copy of a device
   // array, which nothing here uses. Host memory is shared as it is.
   SundialsVector::host_view_type HostView;
   if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                            MemSpace>::accessible) {
      HostView = Kokkos::create_mirror_view(View);
   }

   auto Vec = std::make_unique<SundialsVector>(View, HostView, Ctx);
   enableFusedOps(*Vec);
   return Vec;
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif

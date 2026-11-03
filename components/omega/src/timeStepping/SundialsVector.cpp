//===-- SundialsVector.cpp - Omega arrays as SUNDIALS vectors ---*- C++ -*-===//
//
// Fused vector operations for SUNDIALS Kokkos vectors in Omega's memory space
//
//===----------------------------------------------------------------------===//

#include "SundialsVector.h"

namespace OMEGA {

namespace {

// Most terms one kernel adds. ERKStep passes one per stage plus the start
// state, and the largest SUNDIALS explicit table has 16 stages.
constexpr int MaxFusedTerms = 32;

// Coefficients and data of the nonzero terms, passed to the kernel by value
struct FusedTerms {
   int NTerms = 0;
   sunrealtype Coeff[MaxFusedTerms];
   const sunrealtype *Data[MaxFusedTerms];
};

sunrealtype *deviceData(N_Vector Vec) {
   return sundials::kokkos::GetVec<SundialsVector>(Vec)->View().data();
}

// Clone that keeps the fused operations, so they carry over to every
// vector ARKODE derives from a wrapped Omega array
N_Vector cloneWithFusedOps(N_Vector Template) {
   SundialsVector *Clone = SundialsVector::Clone(
       *sundials::kokkos::GetVec<SundialsVector>(Template), Template->sunctx);
   N_Vector Vec = Clone->get();
   enableFusedOps(Vec);
   return Vec;
}

} // namespace

//------------------------------------------------------------------------------
// Fused linear combination Z = sum_j C[j] * X[j] in one kernel
SUNErrCode fusedLinearCombination(int NVec, sunrealtype *C, N_Vector *X,
                                  N_Vector Z) {

   // Beyond the kernel's capacity, use the generic SUNDIALS definition
   if (NVec > MaxFusedTerms) {
      N_VScale(C[0], X[0], Z);
      for (int J = 1; J < NVec; ++J) {
         N_VLinearSum(C[J], X[J], 1.0, Z, Z);
      }
      return SUN_SUCCESS;
   }

   FusedTerms Terms;
   auto addTerm = [&](int J) {
      if (C[J] != 0) {
         Terms.Coeff[Terms.NTerms] = C[J];
         Terms.Data[Terms.NTerms]  = deviceData(X[J]);
         ++Terms.NTerms;
      }
   };
   addTerm(NVec - 1);
   for (int J = 0; J < NVec - 1; ++J) {
      addTerm(J);
   }

   sunrealtype *ZData = deviceData(Z);
   const int NVals    = N_VGetLength(Z);

   // Each value is read before it is written, so Z may also be one of X
   parallelFor(
       "SundialsLinearCombination", {NVals}, KOKKOS_LAMBDA(int I) {
          sunrealtype Sum = 0;
          for (int J = 0; J < Terms.NTerms; ++J) {
             Sum += Terms.Coeff[J] * Terms.Data[J][I];
          }
          ZData[I] = Sum;
       });

   return SUN_SUCCESS;
}

//------------------------------------------------------------------------------
// Gives a SundialsVector the fused linear combination
void enableFusedOps(N_Vector Vec) {
   Vec->ops->nvlinearcombination = fusedLinearCombination;
   Vec->ops->nvclone             = cloneWithFusedOps;
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//

//===-- Test driver for OMEGA GSW-C library -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA GSW-C external library
///
/// This driver tests that the GSW-C library can be called
/// and returns expected value (as published in Roquet et al 2015)
//
//===-----------------------------------------------------------------------===/

#include "Eos.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Halo.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "mpi.h"

// added for debug
#include "AuxiliaryState.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "VertCoord.h"

#include <gswteos-10.h>

using namespace OMEGA;

/// Test constants and expected values
constexpr int NVertLayers = 60;

/// Published values (TEOS-10 and linear) to test against
const Real TeosSVExpValue =
    0.0009732819628; // Expected value for TEOS-10 specific volume
const Real LinearExpValue =
    0.0009784735812133072; // Expected value for Linear specific volume
const Real TeosBVFExpValue =
    0.020911744281268286; // Expected value for TEOS-10 Brunt-Vaisala frequency
const Real LinearBVFExpValue =
    0.017833905406889013; // Expected value for Linear Brunt-Vaisala frequency
const Real GswBVFExpValue =
    0.02081197958166906; // Expected value from GSW-C library

/// Test input values
const Real Sa = 30.0;   // Absolute Salinity in g/kg
const Real Ct = 10.0;   // Conservative Temperature in degC
const Real P  = 1000.0; // Pressure in dbar

const I4 KDisp  = 1;     // Displate parcel to K=1 for TEOS-10 eos
const Real RTol = 1e-10; // Relative tolerance for isApprox checks

/// The initialization routine for Eos testing. It calls various
/// init routines, including the creation of the default decomposition.
void initEosTest(const std::string &mesh) {

   int Err = 0;

   /// Initialize the Machine Environment class - this also creates
   /// the default MachEnv. Then retrieve the default environment and
   /// some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   /// Initialize logging
   initLogging(DefEnv);
   LOG_INFO("------ EOS Unit Tests ------");

   /// Open and read config file
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   /// Initialize parallel IO
   IO::init(DefComm);

   /// Initialize decomposition
   Decomp::init(mesh);

   // Initialize streams
   IOStream::init();

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      LOG_ERROR("VertCoordTest: error initializing default halo");

   /// Initialize vertical coordinate (phase 1)
   VertCoord::init1();

   /// Initialize mesh
   HorzMesh::init();

   // Initialize the vertical coordinate (phase 2)
   VertCoord::init2();
   auto VCoord         = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;

   /// Initialize Eos
   Eos::init();

   /// Retrieve Eos
   Eos *DefEos = Eos::getInstance();
   if (!DefEos)
      ABORT_ERROR("EosTest: Eos retrieval FAIL");
}

/// Test Linear EOS calculation for all cells/layers
void testEosLinear() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Compute specific volume
   TestEos->computeSpecVol(TArray, SArray, PArray);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches   = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduceOuter(
       "CheckSpecVolMatrix-linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), LinearExpValue, RTol)) {
                    InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   auto SpecVolH = createHostMirrorCopy(SpecVol);

   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (!isApprox(SpecVolH(ICell, K), LinearExpValue, RTol)) {
               LOG_ERROR("EosTest: SpecVol Linear Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         ICell, K, SpecVolH(ICell, K), LinearExpValue);
            }
         }
      }
      ABORT_ERROR("EosTest: SpecVol Linear FAIL with {} bad values",
                  NumMismatches);
   }
   return;
}

/// Test Linear EOS calculation with vertical displacement
void testEosLinearDisplaced() {
   /// Get mesh and coord info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVolDisplaced, 0.0);

   /// Compute displaced specific volume
   TestEos->computeSpecVolDisp(TArray, SArray, PArray, KDisp);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches            = 0;
   Array2DReal SpecVolDisplaced = TestEos->SpecVolDisplaced;
   parallelReduceOuter(
       "CheckSpecVolDispMatrix-linear", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVolDisplaced(ICell, K), LinearExpValue,
                               RTol)) {
                    InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);

   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (!isApprox(SpecVolDisplacedH(ICell, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Linear Displaced Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         ICell, K, SpecVolDisplacedH(ICell, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: Linear SpecVolDisp FAIL with {} bad values ",
                  NumMismatches);
   }
   return;
}

/// Test linear Brunt-Vaisala frequency calculation for all cells/layers
void testBruntVaisalaFreqLinear() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::LinearEos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use deep copy to initialize results to zero
   deepCopy(TestEos->SpecVol, 0.0);
   deepCopy(TestEos->BruntVaisalaFreq, 0.0);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   // fill remaining entries with sample values that should lead to ref result
   // for K = 1.
   OMEGA_SCOPE(ZMid, VCoord->ZMid);

   parallelForOuter(
       "populateArrays", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const int K = KMin + KChunk;

                 if (K == 0) {
                    ZMid(ICell, 0)   = -992.1173890198451_Real;
                    SArray(ICell, 0) = Sa - 1.0_Real;
                    TArray(ICell, 0) = Ct + 15.0_Real;
                    PArray(ICell, 0) = P;
                 } else if (K == 1) {
                    ZMid(ICell, 1)   = -993.1071379053125_Real;
                    SArray(ICell, 1) = Sa;
                    TArray(ICell, 1) = Ct + 10.0_Real;
                    PArray(ICell, 1) = P + 1.0_Real;
                 } else if (K == 2) {
                    ZMid(ICell, 2)   = -994.0968821072275_Real;
                    SArray(ICell, 2) = Sa + 1.0_Real;
                    TArray(ICell, 2) = Ct + 5.0_Real;
                    PArray(ICell, K) = P + 2.0_Real;
                 } else { // fill rest to valid junk to avoid NaNs or Inf
                    ZMid(ICell, K)   = -994.0968821072275_Real - 0.1_Real * K;
                    SArray(ICell, K) = Sa + 1.0_Real + 0.1_Real * K;
                    TArray(ICell, K) = Ct + 5.0_Real - 0.01_Real * K;
                    PArray(ICell, K) = P + 2.0_Real + 0.1_Real * K;
                 }
              });
       });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreq(TArray, SArray, PArray, SpecVol);

   /// Check all array values against expected value
   int NumMismatches = 0;
   OMEGA_SCOPE(BruntVaisalaFreq, TestEos->BruntVaisalaFreq);

   parallelReduceOuter(
       "CheckBruntVaisala-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;

                 if (K == 0) { // should be zero
                    if (BruntVaisalaFreq(ICell, K) != 0.0)
                       InnerCount++;
                 } else if (K == 1) { // should be ref value
                    if (!isApprox(BruntVaisalaFreq(ICell, K), LinearBVFExpValue,
                                  RTol))
                       InnerCount++;
                 } else { // just check for unreasonable values
                    if (BruntVaisalaFreq(ICell, K) == 0.0 or
                        Kokkos::isnan(BruntVaisalaFreq(ICell, K)) or
                        Kokkos::isinf(BruntVaisalaFreq(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   // If test fails, print bad values and abort
   auto BruntVaisalaFreqH = createHostMirrorCopy(BruntVaisalaFreq);
   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (K == 0) {
               // top layer should be zero
               if (BruntVaisalaFreqH(ICell, 0) != 0.0)
                  LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                            "BruntVaisala({},{}) = {}; Expected {}",
                            ICell, 0, BruntVaisalaFreqH(ICell, 0), 0.0);
            } else if (K == 1) {
               // K = 1 should be ref value
               if (!isApprox(BruntVaisalaFreqH(ICell, 1), LinearBVFExpValue,
                             RTol))
                  LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                            "BruntVaisala({},{}) = {}; Expected {}",
                            ICell, 1, BruntVaisalaFreqH(ICell, 1),
                            LinearBVFExpValue);
            } else {
               // remaining values just check for other conditions
               if (BruntVaisalaFreqH(ICell, K) == 0.0 or
                   Kokkos::isnan(BruntVaisalaFreqH(ICell, K)) or
                   Kokkos::isinf(BruntVaisalaFreqH(ICell, K)))
                  LOG_ERROR("EosTest: Brunt-Vaisala Linear Bad Value: "
                            "BruntVaisala({},{}) = {}",
                            ICell, K, BruntVaisalaFreqH(ICell, K));
            }
         }
      }
      ABORT_ERROR("EosTest: BruntVaisala Linear FAIL with {} bad values",
                  NumMismatches);
   }
   return;
}

/// Test TEOS-10 EOS calculation for all cells/layers
void testEosTeos10() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Compute specific volume
   TestEos->computeSpecVol(TArray, SArray, PArray);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches   = 0;
   Array2DReal SpecVol = TestEos->SpecVol;
   parallelReduceOuter(
       "CheckSpecVolMatrix-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVol(ICell, K), TeosSVExpValue, RTol)) {
                    InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   auto SpecVolH = createHostMirrorCopy(SpecVol);

   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (!isApprox(SpecVolH(ICell, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol TEOS Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         ICell, K, SpecVolH(ICell, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol TEOS FAIL with {} bad values",
                  NumMismatches);
   }
   return;
}

/// Test TEOS-10 EOS calculation with vertical displacement
void testEosTeos10Displaced() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use Kokkos::deep_copy to fill the entire view with the ref value
   deepCopy(SArray, Sa);
   deepCopy(TArray, Ct);
   deepCopy(PArray, P);
   deepCopy(TestEos->SpecVolDisplaced, 0.0);

   /// Compute displaced specific volume
   TestEos->computeSpecVolDisp(TArray, SArray, PArray, KDisp);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches            = 0;
   Array2DReal SpecVolDisplaced = TestEos->SpecVolDisplaced;
   parallelReduceOuter(
       "CheckSpecVolDispMatrix-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;
                 if (!isApprox(SpecVolDisplaced(ICell, K), TeosSVExpValue,
                               RTol)) {
                    InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   auto SpecVolDisplacedH = createHostMirrorCopy(SpecVolDisplaced);

   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (!isApprox(SpecVolDisplacedH(ICell, K), LinearExpValue, RTol))
               LOG_ERROR("EosTest: SpecVol Displaced TEOS Bad Value: "
                         "SpecVol({},{}) = {}; Expected {}",
                         ICell, K, SpecVolDisplacedH(ICell, K), LinearExpValue);
         }
      }
      ABORT_ERROR("EosTest: SpecVol Displaced TEOS FAIL with {} bad values",
                  NumMismatches);
   }
   return;
}

/// Test TEOS-10 Brunt-Vaisala frequency calculation for all cells/layer
void testBruntVaisalaFreqTeos10() {
   /// Get mesh and coordinate info
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   VCoord->NVertLayers = NVertLayers;
   I4 NCellsAll        = Mesh->NCellsAll;
   /// Get Eos instance to test
   Eos *TestEos       = Eos::getInstance();
   TestEos->EosChoice = EosType::Teos10Eos;

   /// Create and fill ocean state arrays
   Array2DReal SArray = Array2DReal("SArray", NCellsAll, NVertLayers);
   Array2DReal TArray = Array2DReal("TArray", NCellsAll, NVertLayers);
   Array2DReal PArray = Array2DReal("PArray", NCellsAll, NVertLayers);
   /// Use deep copy to initialize results to zero
   deepCopy(TestEos->BruntVaisalaFreq, 0.0);
   deepCopy(TestEos->SpecVol, 0.0);

   /// Fill inputs with values that should lead to ref result for K=1
   OMEGA_SCOPE(ZMid, VCoord->ZMid);
   parallelFor(
       "populateArrays", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          if (K == 0) {
             ZMid(ICell, 0)   = -992.1173890198451_Real;
             SArray(ICell, 0) = Sa - 1.0_Real;
             TArray(ICell, 0) = Ct + 15.0_Real;
             PArray(ICell, 0) = P;
          } else if (K == 1) {
             ZMid(ICell, 1)   = -993.1071379053125_Real;
             SArray(ICell, 1) = Sa;
             TArray(ICell, 1) = Ct + 10.0_Real;
             PArray(ICell, 1) = P + 1.0_Real;
          } else if (K == 2) {
             ZMid(ICell, 2)   = -994.0968821072275_Real;
             SArray(ICell, 2) = Sa + 1.0_Real;
             TArray(ICell, 2) = Ct + 5.0_Real;
             PArray(ICell, K) = P + 2.0_Real;
          } else { // fill rest with valid junk to avoid Nans and Inf
             ZMid(ICell, K)   = -994.0968821072275_Real - 0.1_Real * K;
             SArray(ICell, K) = Sa + 1.0_Real + 0.1_Real * K;
             TArray(ICell, K) = Ct + 5.0_Real - 0.01_Real * K;
             PArray(ICell, K) = P + 2.0_Real + 0.1_Real * K;
          }
       });

   /// Compute specific volume first
   TestEos->computeSpecVol(TArray, SArray, PArray);
   Array2DReal SpecVol = TestEos->SpecVol;

   /// Compute Brunt-Vaisala frequency
   TestEos->computeBruntVaisalaFreq(TArray, SArray, PArray, SpecVol);

   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   /// Check all array values against expected value
   int NumMismatches = 0;
   OMEGA_SCOPE(BruntVaisalaFreq, TestEos->BruntVaisalaFreq);

   parallelReduceOuter(
       "CheckBruntVaisala-Teos", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team, int &OuterCount) {
          int NumMismatchesCol;
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, int &InnerCount) {
                 const int K = KMin + KOff;

                 if (K == 0) { // should be zero at top
                    if (BruntVaisalaFreq(ICell, K) != 0.0)
                       InnerCount++;
                 } else if (K == 1) { // should be ref value
                    if (!isApprox(BruntVaisalaFreq(ICell, K), TeosBVFExpValue,
                                  RTol))
                       InnerCount++;
                 } else { // just check for unreasonable values
                    if (BruntVaisalaFreq(ICell, K) == 0.0 or
                        Kokkos::isnan(BruntVaisalaFreq(ICell, K)) or
                        Kokkos::isinf(BruntVaisalaFreq(ICell, K)))
                       InnerCount++;
                 }
              },
              NumMismatchesCol);

          Kokkos::single(PerTeam(Team),
                         [&]() { OuterCount += NumMismatchesCol; });
       },
       NumMismatches);

   // If test fails, print bad values and abort
   auto BruntVaisalaFreqH = createHostMirrorCopy(BruntVaisalaFreq);

   if (NumMismatches != 0) {

      for (int ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const int KMin = MinLayerCell(ICell);
         const int KMax = MaxLayerCell(ICell);

         for (int K = KMin; K <= KMax; ++K) {
            if (K == 0) {
               if (BruntVaisalaFreqH(ICell, K) != 0.0)
                  LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                            "BruntVaisala({},{}) = {}; Expected {}",
                            ICell, 0, BruntVaisalaFreqH(ICell, 0), 0.0);
            } else if (K == 1) {
               if (!isApprox(BruntVaisalaFreqH(ICell, 1), TeosBVFExpValue,
                             RTol))
                  LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                            "BruntVaisala({},{}) = {}; Expected {}",
                            ICell, 1, BruntVaisalaFreqH(ICell, 1),
                            TeosBVFExpValue);
            } else {
               if (BruntVaisalaFreqH(ICell, K) == 0.0 or
                   Kokkos::isnan(BruntVaisalaFreqH(ICell, K)) or
                   Kokkos::isinf(BruntVaisalaFreqH(ICell, K)))
                  LOG_ERROR("EosTest: Brunt-Vaisala TEOS Bad Value: "
                            "BruntVaisala({},{}) = {}",
                            ICell, K, BruntVaisalaFreqH(ICell, 1));
            }
         }
      }
      ABORT_ERROR("EosTest: BruntVaisala TEOS FAIL with {} bad values",
                  NumMismatches);
   }
   return;
}

/// Finalize and clean up all test infrastructure
void finalizeEosTest() {
   IOStream::finalize();
   TimeStepper::clear();
   Eos::destroyInstance();
   HorzMesh::clear();
   VertCoord::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   MachEnv::removeAll();
}

/// Test that the external GSW-C library returns the expected specific volume
void checkValueGswcSpecVol() {
   const Real RTol = 1e-10;

   /// Get specific volume from GSW-C library
   double SpecVol = gsw_specvol(Sa, Ct, P);
   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(SpecVol, TeosSVExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcSpecVol: SpecVol FAIL, expected {}, got {}",
                  TeosSVExpValue, SpecVol);
   }
   return;
}

/// Test that the external GSW-C library returns the expected N2
void checkValueGswcN2() {
   const Real RTol = 1e-10;

   // Number of intervals (nz)
   int Nz = 2;

   // Input arrays: length nz+1
   double Salt[4]  = {Sa - 1.0, Sa, Sa + 1.0}; // Absolute Salinity (g/kg)
   double Temp[4]  = {Ct + 15.0, Ct + 10.0,
                      Ct + 5.0};            // Conservative Temperature (deg C)
   double Press[4] = {P, P + 1.0, P + 2.0}; // Pressure (dbar)

   // Latitude (degrees north)
   double Latitude[4] = {0.0, 0.0, 0.0};

   // Output arrays: length nz
   double N2[Nz];   // Brunt–Väisälä frequency squared
   double PMid[Nz]; // Midpoint pressure

   /// Get specific volume from GSW-C library
   gsw_nsquared(Salt, Temp, Press, Latitude, Nz, N2, PMid);

   /// Check the value against the expected TEOS-10 value
   bool Check = isApprox(N2[0], GswBVFExpValue, RTol);
   if (!Check) {
      ABORT_ERROR("checkValueGswcN2: N2 FAIL, expected {}, got {}",
                  GswBVFExpValue, N2[0]);
   }
   return;
}

// the main tests (all in one to have the same log):
// Single value test:
// --> test calls the external GSW-C library
// and compares the specific volume to the published value
// Full array tests:
// --> one tests the value on a Eos with linear option
// --> next checks the value on a Eos with linear displaced option
// --> next checks the value of the linear Brunt Vaisala Freq. calculation
// --> next checks the value on a Eos with TEOS-10 option
// --> next checks the value on a Eos with TEOS-10 displaced option
// --> last checks the value of the TOES-10 Brunt Vaisala Freq. calculation
void eosTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initEosTest(MeshFile);
   const auto &Mesh = HorzMesh::getDefault();

   checkValueGswcSpecVol();
   checkValueGswcN2();

   testEosLinear();
   testEosLinearDisplaced();
   testBruntVaisalaFreqLinear();
   testEosTeos10();
   testEosTeos10Displaced();
   testBruntVaisalaFreqTeos10();

   finalizeEosTest();

   return;
}

// The test driver for Eos testing
int main(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   eosTest();

   LOG_INFO("------ EOS Unit Tests Successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   // If we made it here, test is successful
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/

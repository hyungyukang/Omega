# Omega Top-level CMakeLists.txt
#
# - supports two build modes:
#    1) standalone build
#    2) e3sm build


# Included scripts do automatic cmake_policy() PUSH and POP
cmake_policy(SET CMP0011 NEW)

# list command no longer ignores empty elements
cmake_policy(SET CMP0007 NEW)

# Only interpret if() arguments as variables or keywords when unquoted
cmake_policy(SET CMP0054 NEW)

#find_package() uses <PackageName>_ROOT variables.
cmake_policy(SET CMP0074 NEW)

###########################################################
# STEP 1: Setup                                           #
#                                                         #
# Assign values to the Omega build-controlling variables  #
###########################################################

# define variables and macros for Omega build
include(${CMAKE_CURRENT_SOURCE_DIR}/OmegaBuild.cmake)

# common pre-processing for standalone and e3sm mode
common()

# handle build modes
if (NOT DEFINED PROJECT_NAME)
  # enter standalone build

  cmake_minimum_required(VERSION 3.21) # the minimum version for HIP support

  # Collect machine and compiler info from CIME
  init_standalone_build()

  set(CMAKE_CXX_STANDARD 20) # used in E3SM
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  project(${OMEGA_PROJECT_NAME}
          LANGUAGES C CXX
         )

  set(LINKER_LANGUAGE CXX)   # needed to support Cray C compiler wrapper

  # update variables for standalone build
  setup_standalone_build()

else()
  # enter embedded build

  if(DEFINED E3SM_SOURCE_DIR) # E3SM_SOURCE_DIR = E3SM_ROOT/components

    # update variables for e3sm build
    setup_e3sm_build()

  else()
    message(FATAL_ERROR "UNKNOWN OMEGA BUILD MODE")

  endif()

endif()

###########################################################
# STEP 2: Update                                          #
#                                                         #
# update cmake & kokkos variables, and adjust               #
# build-controlling variables                             #
###########################################################

# update CMake and Kokkos variables from the build-controlling variables
update_variables()

###########################################
# place holder for further build-controls #
###########################################

# check if all-required variables have proper values
check_setup()


###########################################################
# STEP 3: Build                                           #
#                                                         #
# build Omega, and generate tests                          #
###########################################################

add_subdirectory(external)
add_subdirectory(src)

if(OMEGA_BUILD_TEST)

  enable_testing()

  if (OMEGA_TEST_CDASH)
    include(CTest)
  endif()

  target_link_libraries(
      ${OMEGA_LIB_NAME}_SP
      PRIVATE
      OmegaLibFlags
  )

  target_compile_definitions(
      ${OMEGA_LIB_NAME}_SP
      PUBLIC
      OMEGA_SINGLE_PRECISION=1
  )
endif()


# build Omega as an E3SM component
if("${OMEGA_BUILD_MODE}" STREQUAL "E3SM")

    set(OCN_SRC
      drivers/coupled/ocn_comp_mct.F90
      drivers/coupled/omega_cxx2f_interface.cpp
      drivers/coupled/omega_f2cxx_interface.F90
      drivers/coupled/omega_cpl_indices.F90
    )

    # MOAB is optional: only present when COMP_INTERFACE=moab (buildlib_cmake
    # passes OMEGA_MOAB_ROOT in that case). The non-MOAB build must be
    # unaffected when it isn't found.
    if(OMEGA_MOAB_ROOT)
      find_package(MOAB CONFIG PATHS ${OMEGA_MOAB_ROOT} NO_DEFAULT_PATH)
    endif()
    if(MOAB_FOUND)
      list(APPEND OCN_SRC drivers/coupled/MoabInterface.cpp)
    endif()

    # create ocn lib
    add_library(ocn ${OCN_SRC})

    # Set f90 compiled modules dir
    set_target_properties(
      ocn
      PROPERTIES
      Fortran_MODULE_DIRECTORY
      ${CMAKE_BINARY_DIR}/cmake/ocn
    )
    target_include_directories(ocn PUBLIC ${CMAKE_BINARY_DIR}/cmake/ocn)

    # Link libraries
    target_link_libraries(
      ocn
      PRIVATE
      ${OMEGA_LIB_NAME}
      OmegaLibFlags
      csm_share
    )

  add_subdirectory(test)

    if(MOAB_FOUND)
      target_compile_definitions(ocn PRIVATE HAVE_MOAB)
      target_include_directories(ocn PRIVATE ${MOAB_INCLUDE_DIRS})
      target_link_libraries(ocn PRIVATE ${MOAB_LIBRARIES})
    endif()

endif()

###########################################################
# STEP 4: Output                                          #
#                                                         #
# prepare outputs according to the output-controlling     #
# variables                                               #
###########################################################

wrap_outputs()

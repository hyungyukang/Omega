#!/bin/sh

export PARMETIS_ROOT=/global/cfs/cdirs/e3sm/software/polaris/pm-cpu/spack/dev_polaris_0_4_0_nvidia_mpich/var/spack/environments/dev_polaris_0_4_0_nvidia_mpich/.spack-env/view

cmake \
   -DOMEGA_CIME_COMPILER=nvidia \
   -DOMEGA_BUILD_TYPE=Release \
   -DOMEGA_CIME_MACHINE=pm-cpu \
   -DOMEGA_PARMETIS_ROOT=${PARMETIS_ROOT}\
   -Wno-dev \
   -S ../ -B .

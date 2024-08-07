#!/bin/sh


# Perlmutter pm-cpu

module load cmake

PARMETIS_ROOT=/global/cfs/cdirs/e3sm/software/polaris/pm-cpu/spack/dev_polaris_0_3_0_gnu_mpich/var/spack/environments/dev_polaris_0_3_0_gnu_mpich/.spack-env/view
cmake \
   -DOMEGA_CIME_COMPILER=gnu  \
   -DOMEGA_CIME_MACHINE=pm-cpu \
   -DOMEGA_PARMETIS_ROOT=${PARMETIS_ROOT}\
   -Wno-dev \
   ../

./omega_build.sh

cp -rf ./default_inputs ./src/

   #-DOMEGA_BUILD_TEST=ON \
   #-S /global/homes/m/mpeterse/repos/omega/pr_tendterms_boneill/components/omega -B .
#./omega_build.sh

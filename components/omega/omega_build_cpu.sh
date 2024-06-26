#!/bin/sh

cmake=/lustre/cyclone/nwp501/proj-shared/hgkang/programs/cmake-3.22.6/bin/cmake

$cmake  -DOMEGA_CIME_COMPILER=nvidiagpu \
        -DOMEGA_CIME_MACHINE=miller-gpu \
        -DOMEGA_PARMETIS_ROOT=/lustre/storm/nwp501/proj-shared/hgkang/programs/ParMETIS \
        ../

#        -DOMEGA_ENABLE_OPENMP=ON \
#        -DOMEGA_ARCH=SERIAL \

#       -DOMEGA_ARCH=SERIAL \
#        -DOMEGA_ARCH=CUDA \
#        -DOMEGA_BUILD_TEST=ON \
#        -DOMEGA_METIS_ROOT=/lustre/storm/nwp501/proj-shared/hgkang/programs/METIS \
#        -DOMEGA_GKLIB_ROOT=/lustre/storm/nwp501/proj-shared/hgkang/programs/GKlib \
#        -DPIO_ENABLE_DOC:BOOL=OFF \
#        -DOMEGA_ARCH=CUDA \

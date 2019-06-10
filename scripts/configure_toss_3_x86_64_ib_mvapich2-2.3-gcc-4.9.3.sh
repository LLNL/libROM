#/bin/bash

mvapich="mvapich2-2.3"
gcc="gcc-4.9.3"

./configure --with-CXX=/usr/tce/packages/mvapich2/${mvapich}-${gcc}/bin/mpicxx --with-hdf5=/usr/tce/packages/hdf5/hdf5-serial-1.8.18-${gcc}/bin/h5cc --with-lapack='-L/usr/gapps/bdiv/toss_3_x86_64_ib/gnu-4.9-${mvapich}/lapack/3.5.0 -llapack'

#/bin/bash

# choose versions

mvapich="mvapich2-2.3"
gcc="gcc-4.9.3"
hdf5=hdf5-serial-1.8.18

# derived vars below here

gnu_major=$( echo $gcc | sed 's/gcc/gnu/' | cut -d. -f1-2 )
gapps=/usr/gapps/bdiv/toss_3_x86_64_ib/${gnu_major}-${mvapich}

./configure \
  --with-CXX=/usr/tce/packages/mvapich2/${mvapich}-${gcc}/bin/mpicxx \
  --with-hdf5=/usr/tce/packages/hdf5/hdf5-serial-1.8.18-${gcc}/bin/h5cc \
  --with-lapack="-L${gapps}/{gnu_major}-${mvapich}/lapack/3.5.0 -llapack" \
  --with-gtest=${gapps}/gmock/1.7.0/gtest

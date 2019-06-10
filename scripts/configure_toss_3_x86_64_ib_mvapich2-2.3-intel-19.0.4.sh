#/bin/bash

# choose versions

mvapich="mvapich2-2.3"
intel="intel-18.0.1"
hdf5="hdf5-serial-1.8.18"

# derived vars below here

intel_major=$( echo $intel | cut -d. -f1 )
gapps=/usr/gapps/bdiv/toss_3_x86_64_ib/${intel_major}-${mvapich}

./configure \
  --with-CXX=/usr/tce/packages/mvapich2/${mvapich}-${intel}/bin/mpicxx \
  --with-hdf5=/usr/tce/packages/hdf5/${hdf5}-${intel}/bin/h5cc \
  --with-lapack='-L${gapps}/lapack/3.5.0 -llapack' \
  --with-gtest=${gapps}/gmock/1.7.0/gtest

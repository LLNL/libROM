#!/usr/bin/env bash

./configure --with-CXX="/usr/tce/packages/mvapich2/mvapich2-2.2-intel-17.0.2/bin/mpicxx -gxx-name=/usr/tce/packages/gcc/gcc-7.1.0/bin/g++" --with-lapack=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-17-mvapich2-2.2/lapack/3.5.0/lib --with-hdf5=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-17-mvapich2-2.2/hdf5/1.8.10p1 --with-zlib=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-17-mvapich2-2.2/zlib/1.2.3/lib --with-gtest=/usr/gapps/bdiv/toss_3_x86_64_ib/intel-17-mvapich2-2.2/gtest/1.8.0 --with-elemental=no --enable-opt=yes --enable-debug=no "$@"

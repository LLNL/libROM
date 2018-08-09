#!/usr/bin/env bash

FC="/usr/tce/packages/mvapich2/mvapich2-2.2-clang-3.9.1/bin/mpif90" ./configure --with-CXX="/usr/tce/packages/mvapich2/mvapich2-2.2-clang-3.9.1/bin/mpicxx" --with-lapack=/usr/gapps/bdiv/toss_3_x86_64_ib/clang-3.9-mvapich2-2.2/lapack/3.5.0/lib --with-hdf5=/usr/gapps/bdiv/toss_3_x86_64_ib/clang-3.9-mvapich2-2.2/hdf5/1.8.10p1 --with-zlib=/usr/gapps/bdiv/toss_3_x86_64_ib/clang-3.9-mvapich2-2.2/zlib/1.2.3/lib --with-gtest=/usr/gapps/bdiv/toss_3_x86_64_ib/clang-3.9-mvapich2-2.2/gtest/1.8.0 --with-elemental=no --enable-opt=yes --enable-debug=no "$@"

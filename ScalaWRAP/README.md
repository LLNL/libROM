ScalaWRAP
=========
This project aims to implement a modern C++ interface to the functionality
contained in the ScaLAPACK library. ScaLAPACK is a collection of routines for
distributed linear algebra; it contains the BLACS (Basic Linear Algebra
Communication Subprograms) and a prototype PBLAS (parallel BLAS), along with
high-level blocked routines for solving linear systems and computing singular
value and eigenvalue decompositions. The original ScaLAPACK is implemented in
a mixture of Fortran 77 and C and is available from
[netlib](http://www.netlib.org/scalapack/); Intel has implemented the ScaLAPACK
standard in the MKL as well. If the MKL implementation is not available to you,
it is suggested to install ScaLAPACK with the
[spack](https://spack.io/) package manager (the package name is
`netlib-scalapack`).

Objectives
----------
The first objective is performance; the abstractions introduced by the library
should not add any noticeable overhead over directly calling the Fortran
subroutines from ScaLAPACK. The second is clarity: working with a ScaLAPACK
matrix should feel for the most part just like working with your favorite
serial linear algebra library's matrix object. Where the two objectives collide,
performance is preferred by default; the main place you may run into this is
that copy construction and assignment are disabled for ScaLAPACK matrices.
The move constructor is not disabled and does a simple pointer copy from one
matrix to another; the moved-from matrix is still valid but does not own its
data any more. Changes to one matrix are reflected in another.

Overview
--------
The main objects used in interacting with the library are the `LocalMatrix`
class, which is a very simple wrapper around a raw pointer on a single process
along with the matrix size and information about memory order, and the
`ScalaMAT` class, which contains all of the functionality needed to interact
with ScaLAPACK and provides an abstracted interface. Data distribution and
collect (un-distribution?) are key parts of working with the distributed
library, and every effort has been made to ensure that the semantics of these
operations are obvious. As an example, the following constructs a 1000 x 1000
ScaLAPACK matrix and then assigns a 100 x 1000 block stored in C (row-major)
order from each process to it:

    #include "ScalaMat.hpp"
    using namespace ScalaWRAP;

    ScalaMat A(1000, 1000);

    // Assuming that local_data is a pointer to the proper array on each rank
    // and there are ten ranks
    for (int rank = 0; rank < 10; ++rank) {
        LocalMatrix local(local_data, 100, 1000, rank);

        // Indexing is 1 based because of the Fortran roots of ScaLAPACK
        A(
            rowrange(rank*100+1, rank*100+100),
            colrange(1, 1000)
         ) = local;
    }

Fine grained control of data distribution is available as well, but the default
parameters should be suitable for a range of problems and internally are
overridden where it is advantageous.

Common linear algebra operations as implemented in the PBLAS are (*will be*)
implemented both in the form of optimized calls for common operations (such as
`pdgemm`) and as overloaded operators, where appropriate. Some operations are
represented lazily either using private class data or a simple wrapper class.
These include the operations of transposition and scaling of a `ScalaMat`, which
do not actually do any extra computation until it is required by another
operation.

High level driver routines are organized in modules; this list is a WIP as
functionality is being added as needed by client projects.

- Eigensolvers (`Eigensolvers.hpp`): Export a high level interface to the
eigensolvers in ScaLAPACK. These include solving the basic problem
`A * x = q * x` for symmetric `A` (interface for complex numbers is not
implemented) as well as the generalized problem `A * x = q * B * x`.

Building
--------
The build system for ScalaWRAP uses CMake; a basic build where all required
libraries (BLAS, LAPACK, MPI, ScaLAPACK) may be found using environment
variables is as simple as

    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make

after which you will find `libscalawrap.a` in the build directory, and the
`all-tests` executable in the `test` subdirectory.

The CMake variables you are likely to want to set are:

* `CMAKE_BUILD_TYPE`: One of "Debug", "Release"
* `CMAKE_CXX_COMPILER`: To override the default compiler detected, which is
probably Gnu if it's in the path. If you override this variable, you should also
set `CMAKE_C_COMPILER` and `CMAKE_Fortran_COMPILER` to the corresponding
compilers so that the ABI is correctly determined.
* `OPTIM_FLAGS`: By default the build will give `-O3` to the compiler, but if
you have specific knowledge of the target architecture you can probably do
better. Specify your optimization flags using this variable.
* `LAPACK_VENDOR`: Set a vendor string as detailed in the CMake docs for the
`FindLAPACK` module; if the requested LAPACK implementation is not found the
configuration fails. Note the default order is to look for MKL, then OpenBLAS,
then all others, and sequential MKL is serached for first.
* `LAPACK_EXECUTION`: One of "threaded", "sequential". Specifying "threaded"
prefers a multithreaded LAPACK/BLAS implementation; by default I assume that
process-level parallelism is preferred and link a sequential LAPACK.
* `SCALAPACK_LIBNAMES`: CMake will search for the MKL ScaLAPACK implementation
if `MKLROOT` is set in the environment, and otherwise search for a library
`libscalapack`. To override this behavior, specify a space-separated list of
library names (e.g. "mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64") that need to
be linked for ScaLAPACK functionality; if they are not found the configuration
fails.
* `SCALAPACK_LIB_DIRS`: If ScaLAPACK libraries are not located in default search
paths, specify a space-separated list of directories where they are located in
this variable.

Copyright
---------
This work Copyright (c) 2019, Lawrence Livermore National Security, LLC. See the
top level COPYRIGHT file for detail.

SPDX-License-Identifier: MIT

The unit tests rely on the [Catch2](https://github.com/catchorg/Catch2)
C++ unit testing library; the single header version included in our repository
is licensed under the
[Boost Software License](https://www.boost.org/LICENSE_1_0.txt).
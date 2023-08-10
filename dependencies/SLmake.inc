############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 15, 2000
#
#  Modified:        October 13, 2011
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS         = -DAdd_

#
#  The fortran and C compilers, loaders, and their flags
#

# LIBROM instruction: copy this file to scalapack-2.2.0 directory and type make.

# -fallow-argument-mismatch is needed to compile some legacy fortran code that does not conform to modern gnu compiler standard.
# -fPIC is needed to build the static scalapack library.

FC            = mpif90
CC            = mpicc 
NOOPT         = -O0 -fallow-argument-mismatch -fPIC
FCFLAGS       = -O3 -fallow-argument-mismatch -fPIC
CCFLAGS       = -O3 -fPIC
FCLOADER      = $(FC)
CCLOADER      = $(CC)
FCLOADFLAGS   = $(FCFLAGS)
CCLOADFLAGS   = $(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#

ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = ranlib

#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

BLASLIB       = -lblas
LAPACKLIB     = -llapack
LIBS          = $(LAPACKLIB) $(BLASLIB)

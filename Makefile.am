###############################################################################
#
#  Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
#  and other libROM project developers. See the top-level COPYRIGHT
#  file for details.
#
#  SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
###############################################################################
ACLOCAL_AMFLAGS = -I config
AM_CPPFLAGS = -I. $(INC_PATH)
AM_LDFLAGS = $(LIB_PATH) $(LIB_NAME)

lib_LIBRARIES = librom.a
librom_a_SOURCES =
librom_a_SOURCES += Matrix.C
librom_a_SOURCES += Vector.C
librom_a_SOURCES += BasisReader.C
librom_a_SOURCES += BasisWriter.C
librom_a_SOURCES += SVDBasisGenerator.C
librom_a_SOURCES += SVDSampler.C
librom_a_SOURCES += SVD.C
librom_a_SOURCES += IncrementalSVD.C
librom_a_SOURCES += IncrementalSVDFastUpdate.C
librom_a_SOURCES += IncrementalSVDStandard.C
librom_a_SOURCES += IncrementalSVDBasisGenerator.C
librom_a_SOURCES += IncrementalSVDSampler.C
librom_a_SOURCES += StaticSVD.C
librom_a_SOURCES += StaticSVDBasisGenerator.C
librom_a_SOURCES += StaticSVDSampler.C
librom_a_SOURCES += Database.C
librom_a_SOURCES += HDFDatabase.C
librom_a_SOURCES += Utilities.C
librom_a_SOURCES += ParallelBuffer.C
librom_a_SOURCES += DEIM.C
librom_a_SOURCES += QDEIM.C
librom_a_SOURCES += GNAT.C
librom_a_SOURCES += scalapack_c_wrapper.c
librom_a_SOURCES += scalapack_f_wrapper.f90

include_HEADERS =
include_HEADERS += BasisReader.h
include_HEADERS += BasisWriter.h
include_HEADERS += CAROM_config.h
include_HEADERS += DEIM.h
include_HEADERS += Database.h
include_HEADERS += HDFDatabase.h
include_HEADERS += IncrementalSVD.h
include_HEADERS += IncrementalSVDBasisGenerator.h
include_HEADERS += IncrementalSVDSampler.h
include_HEADERS += IncrementalSVDStandard.h
include_HEADERS += Matrix.h
include_HEADERS += ParallelBuffer.h
include_HEADERS += QDEIM.h
include_HEADERS += SVD.h
include_HEADERS += SVDBasisGenerator.h
include_HEADERS += SVDSampler.h
include_HEADERS += StaticSVD.h
include_HEADERS += StaticSVDBasisGenerator.h
include_HEADERS += StaticSVDSampler.h
include_HEADERS += Utilities.h
include_HEADERS += Vector.h
include_HEADERS += scalapack_wrapper.h

noinst_PROGRAMS =
noinst_PROGRAMS += uneven_dist
noinst_PROGRAMS += random_test
noinst_PROGRAMS += weak_scaling
noinst_PROGRAMS += smoke_test
noinst_PROGRAMS += smoke_static
noinst_PROGRAMS += load_samples
noinst_PROGRAMS += tests/test_Matrix
noinst_PROGRAMS += tests/test_Vector
noinst_PROGRAMS += tests/test_SVD
noinst_PROGRAMS += tests/test_IncrementalSVD
noinst_PROGRAMS += tests/test_StaticSVD

# Static library linking requires different flag ordering than
# dynamic/shared library linking. Static library linking
# tends to require that:
#
# 1) Libraries with symbols referenced in an object file must come
#    *after* the object file in link commands.
#
# 2) If library B references symbols in library A, then library A must
#    come *after* library B.
#
# 3) If libraries B and A reference each other, then flags for at
#    least one of them may need to be repeated (e.g., "-lB -lA -lB").
#
# Automake uses exactly one of the following two variables to
# determine linker flags that come before the executable and object
# file in the linker line: AM_LDFLAGS *or* target_LDFLAGS; the latter
# takes precedence. Flags in target_LDADD are appended to the linker
# line. To accommodate rule (1), set target_LDFLAGS to the empty
# string, which sets linker flags before the executable and object
# file to the empty string, and move all flags to target_LDADD, which
# puts them after the executable and object file.
#
# To accommodate rule (2), put librom after $(LIB_PATH) (which
# contains all of the -L flags to modify the linker search path), but
# before $(LIB_NAME) (which contains all of the -l flags that specify
# which libraries to link).
#
uneven_dist_SOURCES = uneven_dist.C
uneven_dist_LDFLAGS =
uneven_dist_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

tests_test_StaticSVD_SOURCES = tests/test_StaticSVD.C
tests_test_StaticSVD_LDFLAGS =
tests_test_StaticSVD_LDADD = $(LIB_PATH) librom.a $(LIB_NAME)

random_test_SOURCES = random_test.C
random_test_LDFLAGS =
random_test_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

weak_scaling_SOURCES = weak_scaling.C
weak_scaling_LDFLAGS =
weak_scaling_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

smoke_test_SOURCES = smoke_test.C
smoke_test_LDFLAGS =
smoke_test_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

tests_test_Matrix_SOURCES = tests/test_Matrix.C
tests_test_Matrix_LDFLAGS =
tests_test_Matrix_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

tests_test_Vector_SOURCES = tests/test_Vector.C
tests_test_Vector_LDFLAGS =
tests_test_Vector_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

tests_test_SVD_SOURCES = tests/test_SVD.C
tests_test_SVD_LDFLAGS =
tests_test_SVD_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

tests_test_IncrementalSVD_SOURCES = tests/test_IncrementalSVD.C
tests_test_IncrementalSVD_LDFLAGS =
tests_test_IncrementalSVD_LDADD   = $(LIB_PATH) librom.a $(LIB_NAME)

dox:
	(cd docs && $(MAKE) $@)

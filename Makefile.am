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

noinst_PROGRAMS =
noinst_PROGRAMS += uneven_dist
noinst_PROGRAMS += random_test
noinst_PROGRAMS += weak_scaling
noinst_PROGRAMS += smoke_test
noinst_PROGRAMS += tests/test_Matrix
noinst_PROGRAMS += tests/test_Vector
noinst_PROGRAMS += tests/test_SVD

uneven_dist_SOURCES = uneven_dist.C
uneven_dist_LDADD = librom.a

random_test_SOURCES = random_test.C
random_test_LDADD = librom.a

weak_scaling_SOURCES = weak_scaling.C
weak_scaling_LDADD = librom.a

smoke_test_SOURCES = smoke_test.C
smoke_test_LDADD = librom.a

tests_test_Matrix_SOURCES = tests/test_Matrix.C
tests_test_Matrix_LDADD = -lgtest -lpthread librom.a

tests_test_Vector_SOURCES = tests/test_Vector.C
tests_test_Vector_LDADD = -lgtest -lpthread librom.a

tests_test_SVD_SOURCES = tests/test_SVD.C
tests_test_SVD_LDADD = -lgtest -lpthread librom.a

dox:
	(cd docs && $(MAKE) $@)

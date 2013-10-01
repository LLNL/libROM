.SUFFIXES: .o .C

CXX = mpiicpc
CXXFLAGS = -g -Wall
LD = mpiicpc
LDFLAGS = -g

PETSC_INCDIR = /usr/local/tools/petsc-3.2/include
PETSC_LIBDIR = /usr/local/tools/petsc-3.2/lib
PETSC_LIB = petsc

SLEPC_INCDIR = /usr/gapps/carom/slepc-3.2/include
SLEPC_LIBDIR = /usr/gapps/carom/slepc-3.2/lib
SLEPC_LIB = slepc

LIBROM = libROM.a

INCS = -I$(SLEPC_INCDIR) -I$(PETSC_INCDIR)

LIBS = $(LIBROM) -L$(SLEPC_LIBDIR) -Wl,-R$(SLEPC_LIBDIR) -l$(SLEPC_LIB) -L$(PETSC_LIBDIR) -Wl,-R$(PETSC_LIBDIR) -l$(PETSC_LIB) -llapack

OBJS = SLEPcManager.o \
       incremental_svd.o \
       incremental_svd_rom.o \
       incremental_svd_time_stepper.o \
       static_svd.o \
       static_svd_rom.o \
       static_svd_time_stepper.o

default: lib

all: lib random_test smoke_test

lib: $(OBJS)
	ar ru $(LIBROM) $(OBJS)

random_test: lib random_test.o
	$(LD) $(LDFLAGS) -o random_test random_test.o $(LIBS)

smoke_test: lib smoke_test.o
	$(LD) $(LDFLAGS) -o smoke_test smoke_test.o $(LIBS)

.C.o :; $(CXX) $(CXXFLAGS) $(INCS) -c $*.C

clean:
	rm -f *.o $(LIBROM) random_test smoke_test

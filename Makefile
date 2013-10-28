.SUFFIXES: .o .C

CXX = mpiicpc
CXXFLAGS = -g -Wall
LD = mpiicpc
LDFLAGS = -g

LIBROM = libROM.a

LIBS = $(LIBROM) -llapack

OBJS = matrix.o \
       vector.o \
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

.C.o :; $(CXX) $(CXXFLAGS) -c $*.C

clean:
	rm -f *.o $(LIBROM) random_test smoke_test

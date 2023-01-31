FROM ubuntu:18.04

ENV USERNAME=env

# install sudo
RUN apt-get -yq update && apt-get -yq install sudo

# create and switch to a user
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN useradd --no-log-init --home-dir /home/$USERNAME --create-home --shell /bin/bash $USERNAME
RUN adduser $USERNAME sudo
USER $USERNAME
WORKDIR /home/$USERNAME

# install packages
RUN sudo apt-get install -yq git
RUN sudo apt-get install --no-install-recommends -yq make gcc gfortran libssl-dev
RUN sudo apt-get install -yq libopenblas-dev libmpich-dev libblas-dev liblapack-dev libscalapack-mpi-dev libhdf5-serial-dev
RUN sudo apt-get install -yq vim
RUN sudo apt-get install -yq git-lfs
RUN sudo apt-get install -yq valgrind
RUN sudo apt-get install -yq wget

RUN sudo apt-get clean -q

# download dependencies
ENV LIB_DIR=/home/$USERNAME/dependencies
WORKDIR $LIB_DIR

RUN wget https://github.com/Kitware/CMake/releases/download/v3.13.0/cmake-3.13.0.tar.gz
RUN wget -O mfem-4.5.tar.gz https://github.com/mfem/mfem/archive/refs/tags/v4.5.tar.gz
RUN wget -O hypre-2.20.0.tar.gz https://github.com/hypre-space/hypre/archive/refs/tags/v2.20.0.tar.gz
RUN wget -O parmetis-4.0.3.tar.gz http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
RUN wget -O gslib-1.0.7.tar.gz https://github.com/gslib/gslib/archive/v1.0.7.tar.gz
#RUN wget -O glvis-4.2.tar.gz https://github.com/GLVis/glvis/archive/refs/tags/v4.2.tar.gz
#RUN wget -O metis-4.0.3.tar.gz https://github.com/mfem/tpls/raw/gh-pages/metis-4.0.3.tar.gz

# install cmake-3.13.0
RUN tar -zxvf cmake-3.13.0.tar.gz
WORKDIR ./cmake-3.13.0
RUN ./bootstrap
RUN make
RUN sudo make install
WORKDIR $LIB_DIR

ENV CFLAGS="-fPIC"
ENV CPPFLAGS="-fPIC"
ENV CXXFLAGS="-fPIC"

# install hypre
RUN tar -zxvf hypre-2.20.0.tar.gz
RUN mv hypre-2.20.0 hypre
WORKDIR ./hypre/src/
RUN ./configure --disable-fortran
RUN make -j
WORKDIR $LIB_DIR

# install gslib
ENV MG=YES
RUN tar -zxvf gslib-1.0.7.tar.gz
RUN mv gslib-1.0.7 gslib
WORKDIR ./gslib
RUN make CC=mpicc -j
WORKDIR $LIB_DIR

## install metis
#RUN tar -zxvf metis-4.0.3.tar.gz
#WORKDIR ./metis-4.0.3
#RUN make OPTFLAGS=-Wno-error=implicit-function-declaration
#WORKDIR $LIB_DIR
#RUN ln -s metis-4.0.3 metis-4.0

# install parmetis
RUN tar -zxvf parmetis-4.0.3.tar.gz
WORKDIR ./parmetis-4.0.3
RUN make config
RUN make

# These environment variables are for mfem build.
# Need different values for libROM build.
ENV METIS_DIR=${LIB_DIR}/parmetis-4.0.3
ENV METIS_OPT=-I${METIS_DIR}/metis/include
WORKDIR ${METIS_DIR}/build
# Currently docker cannot save command results to environment variables.
# These variables are for libROM build.
ENV MACHINE_ARCH=Linux-x86_64
RUN ln -s $MACHINE_ARCH lib
WORKDIR ${METIS_DIR}/build/lib/libparmetis
RUN ln -s ./ lib
RUN ln -s ${METIS_DIR}/metis metis
WORKDIR ${METIS_DIR}/build/lib/libmetis
RUN ln -s ./ lib
ENV METIS_LIB="-L${METIS_DIR}/build/lib/libparmetis -lparmetis -L${METIS_DIR}/build/lib/libmetis -lmetis"

ENV CFLAGS=
ENV CPPFLAGS=
ENV CXXFLAGS=

WORKDIR $LIB_DIR

# install mfem
#RUN git clone https://github.com/mfem/mfem.git mfem_parallel
RUN tar -zxvf mfem-4.5.tar.gz
RUN mv mfem-4.5 mfem_parallel
WORKDIR ./mfem_parallel
#RUN git pull
#RUN make serial -j 4
RUN make parallel -j 4 STATIC=NO SHARED=YES MFEM_USE_MPI=YES MFEM_USE_GSLIB=${MG} MFEM_USE_METIS=YES MFEM_USE_METIS_5=YES METIS_DIR="$METIS_DIR" METIS_OPT="$METIS_OPT" METIS_LIB="$METIS_LIB"
RUN ln -s ./ lib
RUN ln -s ./ include
WORKDIR $LIB_DIR
RUN ln -s mfem_parallel mfem

# clean up
WORKDIR $LIB_DIR
RUN rm *.tar.gz

# cmake toolchain file for librom
RUN echo 'set(CMAKE_C_COMPILER mpicc)\n\
set(CMAKE_CXX_COMPILER mpicxx)\n\
set(CMAKE_Fortran_COMPILER mpif90)\n\
set(LIB_DIR /home/env/dependencies)\n\
set(MFEM_DIR ${LIB_DIR}/mfem)\n\
set(HYPRE_DIR ${LIB_DIR}/hypre/src/hypre)\n\
set(PARMETIS_DIR ${LIB_DIR}/parmetis-4.0.3/build/lib/libparmetis)\n\
set(METIS_DIR ${LIB_DIR}/parmetis-4.0.3/build/lib/libmetis)' > ./librom_env.cmake

# flags for libROM cmake
ENV TOOLCHAIN_FILE=$LIB_DIR/librom_env.cmake
ENV BUILD_TYPE=Optimized
ENV USE_MFEM=On
ENV MFEM_USE_GSLIB=On
ENV MFEM_DIR=$LIB_DIR/mfem
ENV HYPRE_DIR=$LIB_DIR/hypre/src/hypre
ENV PARMETIS_DIR=$LIB_DIR/parmetis-4.0.3/build/lib/libparmetis
ENV METIS_DIR=$LIB_DIR/parmetis-4.0.3/build/lib/libmetis

WORKDIR /home

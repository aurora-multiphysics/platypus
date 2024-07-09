# Get base image
FROM ubuntu:22.04

# This is needed or it mpiexec complains because docker runs as root
# Discussion on this issue https://github.com/open-mpi/ompi/issues/4451
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# By default four cores are used to compile
ARG compile_cores=4

# By default we install everything under /opt
ARG WORKDIR=opt

# Install apt dependencies
# hadolint ignore=DL3008
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    build-essential \
    bison \
    flex \
    clang-format \
    clang-tidy \
    curl \
    doxygen \
    git \
    ninja-build \
    pre-commit \
    python3 \
    python3-dev \
    python3-distutils \
    python-is-python3 \
    python3-pip \
    libboost-all-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglew-dev \
    libglm-dev \
    libhdf5-mpich-dev \
    libnetcdf-dev \
    libomp-14-dev \
    libpng-dev \
    libsdl2-dev \
    libssl-dev \
    rsync \
    xxd && \
    rm -rf /var/lib/apt/lists/*

# Install pip packages
# hadolint ignore=DL3042,DL3013
RUN pip install packaging && \
    pip install pyyaml && \
    pip install setuptools && \
    pip install jinja2

# Install CMake
WORKDIR /$WORKDIR
RUN git clone https://github.com/Kitware/CMake.git
WORKDIR $WORKDIR/CMake/
RUN git checkout 78ae759 && \
    ./bootstrap --parallel="$compile_cores" && \
    make -j"$compile_cores" && \
    make install && \
    alias cmake='/$WORKDIR/CMake/bin/cmake'

# Download SuperLU
WORKDIR /$WORKDIR
RUN curl -LJO https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz

# Download PETSc
WORKDIR /$WORKDIR
# hadolint ignore=SC2006
RUN mkdir petsc
WORKDIR /$WORKDIR/petsc
RUN curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.3.tar.gz && \
    tar -xf petsc-3.19.3.tar.gz -C .

# Build PETSc
WORKDIR /$WORKDIR/petsc/petsc-3.19.3
# hadolint ignore=SC2006
RUN ./configure \
    --prefix=/$WORKDIR/petsc \
    CXXOPTFLAGS='-O3 -funroll-loops' COPTFLAGS='-O3 -funroll-loops' FOPTFLAGS='-O3' \
    --with-debugging=0 \
    --with-ssl=0 \
    --with-pic=1 \
    --with-openmp=1 \
    --with-mpi=1 \
    --with-shared-libraries=1 \
    --with-fortran-bindings=0 \
    --with-sowing=0 \
    --download-hypre=1 \
    --download-fblaslapack=1 \
    --download-metis=1 \
    --download-ptscotch=1 \
    --download-parmetis=1 \
    --download-superlu_dist=/$WORKDIR/superlu_dist-8.1.0.tar.gz \
    --download-scalapack=1 \
    --download-mumps=1 \
    --download-slepc=1 \
    --with-64-bit-indices=1 \
    PETSC_DIR=`pwd` PETSC_ARCH=arch-linux-c-opt && \
    make

# Install PETSc
WORKDIR /$WORKDIR/petsc/petsc-3.19.3
RUN make PETSC_DIR=/$WORKDIR/petsc/petsc-3.19.3 PETSC_ARCH=arch-linux-c-opt install && \
    make PETSC_DIR=/$WORKDIR/petsc/petsc-3.19.3 PETSC_ARCH=arch-linux-c-opt check
WORKDIR /$WORKDIR
RUN export PETSC_DIR=/$WORKDIR/petsc

# Build MFEM and common miniapp
WORKDIR /$WORKDIR
RUN git clone https://github.com/Heinrich-BR/mfem.git
WORKDIR /$WORKDIR/mfem 
RUN git checkout SubmeshBoundary && mkdir build
WORKDIR /$WORKDIR/mfem/build
RUN cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POSITION_INDEPENDENT_CODE=YES \
    -DMFEM_USE_OPENMP=NO \
    -DMFEM_THREAD_SAFE=NO \
    -DHYPRE_DIR=/$WORKDIR/petsc/ \
    -DBLAS_LIBRARIES=/$WORKDIR/petsc/lib/libfblas.a \
    -DLAPACK_LIBRARIES=/$WORKDIR/petsc/lib/libflapack.a \
    -DMFEM_USE_MPI=YES \
    -DMFEM_USE_METIS_5=YES \
    -DMETIS_DIR=/$WORKDIR/petsc/ \
    -DParMETIS_DIR=/$WORKDIR/petsc/ \
    -DMFEM_USE_SUPERLU=YES \
    -DSuperLUDist_DIR=/$WORKDIR/petsc/ \
    -DMFEM_USE_NETCDF=YES \
    -DHDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/mpich/ \
    make -j"$compile_cores"
WORKDIR /$WORDIR/mfem/build/miniapps/common
RUN make -j"$compile_cores"

# Build GLVis
WORKDIR /$WORKDIR
RUN git clone https://github.com/GLVis/glvis.git
WORKDIR /$WORKDIR/glvis
RUN mkdir build
WORKDIR /$WORKDIR/glvis/build/
RUN cmake -G Ninja -DMFEM_DIR=/$WORKDIR/mfem/build .. && \
    ninja

# Configure MOOSE
WORKDIR /$WORKDIR
RUN export MOOSE_JOBS=$compile_cores && \
    git clone https://github.com/idaholab/moose
WORKDIR /$WORKDIR/moose
RUN git checkout master && \
    export PETSC_DIR=/$WORKDIR/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    export CC=mpicc && \
    export CXX=mpicxx && \
    export F90=mpif90 && \
    export F77=mpif77 && \
    export FC=mpif90  && \
    METHODS=opt ./scripts/update_and_rebuild_libmesh.sh --with-mpi && \
    ./configure --with-derivative-size=200

# Build WASP
WORKDIR /$WORKDIR/moose/scripts
RUN METHODS=opt ./update_and_rebuild_wasp.sh

# Build MOOSE
WORKDIR /$WORKDIR/moose/framework
RUN METHOD=opt make -j$compile_cores

WORKDIR /$WORKDIR/moose/test
RUN METHOD=opt make -j$compile_cores
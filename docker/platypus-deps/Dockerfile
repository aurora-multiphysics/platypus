# Ignore warnings for not pinning versions in pip, apt-get install.
# hadolint global ignore=DL3008,DL3013

# Get base image
FROM ubuntu:24.04

# This is needed or it mpiexec complains because docker runs as root
# Discussion on this issue https://github.com/open-mpi/ompi/issues/4451
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# By default four cores are used to compile
ARG compile_cores=4

# By default we install everything under /opt
ARG WORKDIR=opt

# Install apt dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends \
    build-essential \
    bison \
    flex \
    clang-format \
    clang-tidy \
    cmake \
    curl \
    doxygen \
    gfortran \
    git \
    graphviz \
    ninja-build \
    pre-commit \
    python3-deepdiff \
    python3-jinja2 \
    python3-livereload \
    python3-pybtex \
    python3-pylatexenc \
    python3-xmltodict \
    python-is-python3 \
    libboost-all-dev \
    libfontconfig1-dev \
    libglew-dev \
    libglm-dev \
    libhdf5-openmpi-dev \
    libnetcdf-dev \
    libsdl2-dev \
    libtirpc-dev \
    rsync \
    xxd && \
    rm -rf /var/lib/apt/lists/*

# Download SuperLU
WORKDIR /$WORKDIR
RUN curl -LJO https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz

# Download PETSc
WORKDIR /$WORKDIR
RUN mkdir petsc
WORKDIR /$WORKDIR/petsc
RUN curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.3.tar.gz && \
    tar -xf petsc-3.19.3.tar.gz -C .

# Build and install PETSc
WORKDIR /$WORKDIR/petsc/petsc-3.19.3
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
    --with-64-bit-indices=1 && \
    make && make install && make check

# Download and build Conduit
WORKDIR /$WORKDIR
RUN git clone --recursive https://github.com/llnl/conduit.git
WORKDIR /$WORKDIR/conduit
RUN mkdir build
WORKDIR /$WORKDIR/conduit/build
RUN cmake  ../src \
    -DCMAKE_INSTALL_PREFIX=/$WORKDIR/conduit/installed \
    -DCMAKE_BUILD_TYPE="Release" \
    -DENABLE_MPI=YES
RUN make install -j$compile_cores

# Download and build CEED
WORKDIR /$WORKDIR
RUN git clone https://github.com/CEED/libCEED.git
WORKDIR /$WORKDIR/libCEED
# This is the v0.12 commit
# The most recent main pushes won't build with MFEM
RUN git checkout 4018a20
RUN make install prefix=/$WORKDIR/libCEED/build

# Build MFEM and common miniapp
WORKDIR /$WORKDIR
RUN git clone https://github.com/mfem/mfem.git
WORKDIR /$WORKDIR/mfem 
RUN git checkout master && mkdir build
WORKDIR /$WORKDIR/mfem/build
RUN cmake .. \
    -DCMAKE_INSTALL_PREFIX=/$WORKDIR/mfem/installed \
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
    -DMFEM_USE_CONDUIT=YES \
    -DMFEM_USE_CEED=YES \
    -DCEED_DIR=/$WORKDIR/libCEED/build \
    -DCONDUIT_DIR=/$WORKDIR/conduit/installed \
    -DHDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ && \
    make install -j$compile_cores
WORKDIR /$WORKDIR/mfem/build/miniapps/common
RUN make install -j$compile_cores

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
RUN git clone https://github.com/idaholab/moose
WORKDIR /$WORKDIR/moose
RUN git checkout master && \
    export MOOSE_JOBS=$compile_cores && \
    export PETSC_DIR=/$WORKDIR/petsc && \
    export CC=mpicc && \
    export CXX=mpicxx && \
    export FC=mpif90  && \
    METHODS=opt ./scripts/update_and_rebuild_libmesh.sh && \
    METHODS=opt ./scripts/update_and_rebuild_wasp.sh && \
    ./configure --with-derivative-size=200

# Build MOOSE
WORKDIR /$WORKDIR/moose/framework
RUN METHOD=opt make -j$compile_cores

WORKDIR /$WORKDIR/moose/test
RUN METHOD=opt make -j$compile_cores

# Build HEAT_TRANSFER module 
WORKDIR /$WORKDIR/moose/modules/heat_transfer
RUN METHOD=opt make -j$compile_cores

# Build SOLID_MECHANICS module 
WORKDIR /$WORKDIR/moose/modules/solid_mechanics
RUN METHOD=opt make -j$compile_cores

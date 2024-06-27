# Build image containing Hephaestus dependencies
# Get base image
FROM ubuntu:22.04

# This is needed or it mpiexec complains because docker runs as root
# Discussion on this issue https://github.com/open-mpi/ompi/issues/4451
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# By default one core is used to compile
ARG compile_cores=4

# By default checkout mfem-moose branch
ARG WORKDIR="opt"

# Install apt dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y \
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
    xxd

# Install pip packages
RUN pip install packaging && \
    pip install pyyaml && \
    pip install setuptools && \
    pip install jinja2

# Install GSLib
RUN cd /$WORKDIR && \
    git clone https://github.com/Nek5000/gslib.git && \
    cd gslib/ && \
    make CFLAGS='-O2 -fPIC' DEBUG=1

# Install CMake
RUN cd /$WORKDIR && \
    git clone https://github.com/Kitware/CMake.git && \
    cd CMake/ && \
    git checkout 78ae759 && \
    ./bootstrap --parallel="$compile_cores" && \
    make -j"$compile_cores" && \
    make install && \
    alias cmake='/$WORKDIR/CMake/bin/cmake'

# Download SuperLU
RUN cd /$WORKDIR && \
    curl -LJO https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz

# Download PETSc
RUN cd /$WORKDIR && \
    mkdir petsc && \
    cd petsc && \
    curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.3.tar.gz && \
    tar -xf petsc-3.19.3.tar.gz -C .

# Build PETSc
RUN cd /$WORKDIR/petsc/petsc-3.19.3 && \
    ./configure \
    --prefix=/$WORKDIR/petsc \
    CXXOPTFLAGS='-O0' COPTFLAGS='-O0' FOPTFLAGS='-O0' \
    --with-debugging=1 \
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
RUN cd /$WORKDIR/petsc/petsc-3.19.3 && \
    make PETSC_DIR=/$WORKDIR/petsc/petsc-3.19.3 PETSC_ARCH=arch-linux-c-opt install && \
    make PETSC_DIR=/$WORKDIR/petsc/petsc-3.19.3 PETSC_ARCH=arch-linux-c-opt check && \
    cd .. && \
    cd .. && \
    export PETSC_DIR=/$WORKDIR/petsc

# Build MFEM and common miniapp
RUN cd /$WORKDIR && \
    git clone https://github.com/Heinrich-BR/mfem.git && \
    cd mfem && \
    git checkout SubmeshBoundary && \
    mkdir build && \
    cd build && \
    cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
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
    -DMFEM_USE_GSLIB=YES \
    -DGSLIB_DIR=/$WORKDIR/gslib/build && \
    make -j"$compile_cores" && \
    cd miniapps/common && \
    make -j"$compile_cores"

# Build GLVis
RUN cd /$WORKDIR && \
    git clone https://github.com/GLVis/glvis.git && \
    cd glvis/ && \
    mkdir build && \
    cd build/ && \
    cmake -G Ninja -DMFEM_DIR=/$WORKDIR/mfem/build .. && \
    ninja

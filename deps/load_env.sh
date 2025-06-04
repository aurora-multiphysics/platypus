#!/bin/bash

printf 'Make sure that spack has been loaded and activated from the home directory\n'

spack load gcc

export SLU_DIR=$(spack location -i superlu-dist)
export HDF5_DIR=$(spack location -i hdf5)
export SLEPC_DIR=$(spack location -i slepc)
export PETSC_DIR=$(spack location -i petsc)
export CONDUIT_DIR=$(spack location -i conduit)
export TIRPC_DIR=$(spack location -i libtirpc)
export CEED_DIR=$(spack location -i libceed)
export GSLIB_DIR=$(spack location -i gslib)
export TIRPC_DIR=$(spack location -i libtirpc)

spack load py-deepdiff
spack load py-jinja2
spack load py-packaging
spack load py-pyyaml
spack load py-setuptools
spack load py-xmltodict
spack load py-pip

export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CACHE_PATH=$PWD

export OMPI_CC=$(which gcc)
export OMPI_CXX=$(which g++)
export CXX=mpic++
export CC=mpicc
export F90=mpif90
export F77=mpif77
export FC=mpif90

export MOOSE_DIR=$PWD/moose
export MFEM_DIR=$PWD/mfem/installed
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MFEM_DIR}/lib

export CPPFLAGS="${CPPFLAGS} -I${TIRPC_DIR}/include/tirpc -gz=none"
export CXXFLAGS="${CXXFLAGS} -D_GLIBCXX_USE_CXX11_ABI=1 -I${TIRPC_DIR}/include/tirpc -gz=none"
export LDFLAGS="${LDFLAGS} -L${TIRPC_DIR}/lib"

export compile_cores=64
export MOOSE_JOBS=$compile_cores
export LIBMESH_JOBS=$compile_cores
export METHOD="dbg"
export CUDA_MFEM="NO"

spack load petsc
spack load cmake
spack load bzip2
spack load petsc

export BUILD_PATH=$PWD

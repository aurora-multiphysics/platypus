#!/bin/bash

printf 'This is for csd3 only!!\n'
printf 'Make sure this script is run from within the deps folder!!\n'

module load gcc/9.4.0/gcc-11.2.0-72sgv5z
source $PWD/spack/share/spack/setup-env.sh
spacktivate platypus

export SLU_DIR=$(spack location -i superlu-dist)
export HDF5_DIR=$(spack location -i hdf5)
export SLEPC_DIR=$(spack location -i slepc)
export PETSC_DIR=$(spack location -i petsc)
export CONDUIT_DIR=$(spack location -i conduit)
export TIRPC_DIR=$(spack location -i libtirpc)
export CEED_DIR=$(spack location -i libceed)
export GSLIB_DIR=$(spack location -i gslib)
export TIRPC_DIR=$(spack location -i libtirpc)

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

export CPPFLAGS="${CPPFLAGS} -I${TIRPC_DIR}/include/tirpc -fpermissive"
export CXXFLAGS="${CXXFLAGS} -fpermissive"
export LDFLAGS="${LDFLAGS} -L${TIRPC_DIR}/lib"

export compile_cores=64
export MOOSE_JOBS=$compile_cores
export LIBMESH_JOBS=$compile_cores
export METHOD="dbg"
export CUDA_MFEM="YES"

spack load petsc

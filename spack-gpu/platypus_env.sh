export SLU_DIR=$(spack find --format "{prefix}" superlu-dist)
export HDF5_DIR=$(spack find --format "{prefix}" hdf5)
export SLEPC_DIR=$(spack find --format "{prefix}" slepc)
export PETSC_DIR=$(spack find --format "{prefix}" petsc)

export compile_cores=16
export OMPI_CC=clang
export OMPI_CXX=clang++
export CXX=mpic++
export CC=mpicc
export F90=mpif90
export F77=mpif77
export FC=mpif90
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BUILD_PATH}/mfem/build:${BUILD_PATH}/mfem/build/miniapps/common
export MOOSE_JOBS=$compile_cores
export LIBMESH_JOBS=$compile_cores
export METHOD="opt"

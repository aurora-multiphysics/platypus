#!/bin/bash
set -e

# By default checkout main branch from aurora-multiphysics/platypus
build_git_sha=main
build_git_repo=aurora-multiphysics/platypus

git clone https://github.com/spack/spack.git -b v0.22.1

# Disable shellcheck lint as file not in repository
# shellcheck disable=SC1091
. spack/share/spack/setup-env.sh
spack compiler find
spack env create platypus spack.yaml
spack env activate platypus
spack install

compile_cores=$(nproc)
export SPACK_VIEW=$SPACK_ENV/.spack-env/view
export METHOD=dbg
export LIBMESH_JOBS=$compile_cores
export MOOSE_JOBS=$compile_cores
export LDFLAGS="-L$LD_LIBRARY_PATH"
export ADDITIONAL_CPPFLAGS="-DLIBMESH_ENABLE_UNIQUE_ID"
export FC=mpif90
export F77=mpif77
export F90=mpif90
export CC=mpicc
export CXX=mpicxx

# build mfem
git clone https://github.com/mfem/mfem.git
cd mfem
git checkout master
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_POSITION_INDEPENDENT_CODE=YES \
    -DMFEM_THREAD_SAFE=NO \
    -DMFEM_USE_OPENMP=NO \
    -DMFEM_USE_MPI=YES \
    -DMFEM_USE_METIS_5=YES \
    -DMFEM_USE_SUPERLU=YES \
    -DMFEM_USE_NETCDF=YES \
    -DMFEM_USE_CONDUIT=YES \
    -DHYPRE_DIR="$SPACK_VIEW" \
    -DSuperLUDist_DIR="$SPACK_VIEW" \
    -DHDF5_DIR="$SPACK_VIEW" \
    -DCONDUIT_DIR="$SPACK_VIEW" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=YES \
    -DCMAKE_INSTALL_PREFIX=./installed
cmake --build build -j "$compile_cores" --verbose
cmake --build build/miniapps/common/ -j "$compile_cores" --verbose
cmake --install build
cd ..

# build libmesh, wasp and moose
git clone https://github.com/idaholab/moose
git checkout master
cd moose
./configure --with-derivative-size=200
./scripts/update_and_rebuild_libmesh.sh \
    --with-cxx=mpicxx \
    --with-cc=mpicc \
    --with-fc=mpif90 \
    --with-f77=mpif77 \
    --with-f90=mpif90 \
    --with-mpi="$SPACK_VIEW" \
    --with-mpi-include="$SPACK_VIEW" \
    --with-mpi-lib="$SPACK_VIEW" \
    --with-netcdf="$SPACK_VIEW"
./scripts/update_and_rebuild_wasp.sh
make -C framework -j "$compile_cores" -B
cd ..

# build platypus
git clone https://github.com/"$build_git_repo"
cd platypus
git checkout "$build_git_sha"
make -j "$compile_cores" -B
cd ..

# test platypus
cd platypus
make test -j "$compile_cores"
make -C unit test -j "$compile_cores"
cd ..

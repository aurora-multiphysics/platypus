set -e

# By default checkout main branch from aurora-multiphysics/platypus
# ARG build_git_sha=main
# ARG build_git_repo=aurora-multiphysics/platypus

git clone https://github.com/spack/spack.git -b v0.22.1
. spack/share/spack/setup-env.sh

spack compiler find
spack env create platypus spack.yaml
spack env activate platypus
spack install

export SPACK_VIEW=$SPACK_ENV/.spack-env/view
export METHOD=dbg
export LIBMESH_JOBS=$(nproc)
export MOOSE_JOBS=$(nproc)
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
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_POSITION_INDEPENDENT_CODE=YES \
-DMFEM_THREAD_SAFE=NO \
-DMFEM_USE_OPENMP=NO \
-DMFEM_USE_MPI=YES \
-DMFEM_USE_METIS_5=YES \
-DMFEM_USE_SUPERLU=YES \
-DMFEM_USE_NETCDF=YES\
-DHYPRE_DIR=$SPACK_VIEW \
-DSuperLUDist_DIR=$SPACK_VIEW \
-DHDF5_DIR=$SPACK_VIEW
cmake --build build -j$(nproc) --verbose
cmake --build build/miniapps/common/ -j$(nproc) --verbose
cd ..

# build libmesh, wasp and moose
git clone https://github.com/idaholab/moose
cd moose
./configure --with-derivative-size=200
./scripts/update_and_rebuild_libmesh.sh \
    --with-mpi=$SPACK_ENV/.spack-env/view/ \
    --with-mpi-include=$SPACK_ENV/.spack-env/view/include/ \
    --with-mpi-lib=$SPACK_ENV/.spack-env/view/lib/
./scripts/update_and_rebuild_wasp.sh
make -C framework -j $(nproc) -B
cd ..

# build platypus
git clone https://github.com/$build_git_repo
cd platypus
git checkout $build_git_sha
make -j $(nproc) -B

# test platypus
make test -j $(nproc)
make -C unit test -j$(nproc)
cd ..

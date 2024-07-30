#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=01:30:00
#SBATCH --mail-type=none
#SBATCH -p sapphire
#SBATCH -A UKAEA-AP001-CPU
#SBATCH --cpus-per-task=1
#SBATCH -o out_%j_%A_%a
#SBATCH --exclusive

# shellcheck source=/dev/null
. /etc/profile.d/modules.sh

# TODO:
# - Update the platypus-commit to master once the Platypus PR which makes
#   the Hephaestus code no longer a dependency is merged-in.

function load_modules() {
    module purge
    module load rhel8/default-icl cmake
    module load intel-oneapi-mkl/2022.1.0/intel/mngj3ad6
    module load gcc/11
    module load curl
    module load python/3.8
    module load ninja

    module load py-pyaml
    module load py-jinja2
    module load py-packaging
    module load py-setuptools
    module load py-deepdiff
    module load py-xmltodict

    STACK_SRC=$(mktemp -d /tmp/moose_stack_src.XXXXXX)
    export STACK_SRC

    WORKDIR=$(pwd)
    export WORKDIR

    export NUM_COMPILE_CORES=8

    export OMPI_MCA_mca_base_component_show_load_errors=0

    USER=$(whoami)
    BUILD_DIR_NAME=${WORKDIR}

    ROOT_PATH=/home/${USER}/rds/rds-ukaea-ap001/${USER}
    BUILD_PATH=${ROOT_PATH}/${BUILD_DIR_NAME}

    HDF5_MAJ_VER=1.10
    HDF5_MIN_VER=10
    HDF5_DIR_NAME=hdf5
    HDF5_INSTALL_PATH=${WORKDIR}/${HDF5_DIR_NAME}

    PETSC_DIR_NAME=petsc

    MOOSE_COMMIT=4e99faf9804480e7be302895ff9b8ded5b9944ea

    PLATYPUS_COMMIT=7cfad0fe4488da867c2c949447a2cbf037032e33

    export PATH=${BUILD_PATH}:${PATH}

    cd "{$WORKDIR}" || exit

    #Need to set some compiler flags via config file"
    echo "-std=c++17" >> icpx.cfg
    echo "-Wno-tautological-constant-compare" >> icpx.cfg
    export ICPXCFG=${WORKDIR}/icpx.cfg

    export CC=mpiicc
    export CXX=mpiicpc
    export FC=mpiifort
    export F77=mpiifort

    export I_MPI_CC=icx
    export I_MPI_CXX=icpx
    export I_MPI_F90=ifx
    export I_MPI_FC=ifx
    export I_MPI_F77=ifx
}

function build_hdf5() {
    cd "${WORKDIR}" || exit
    mkdir -p ${HDF5_DIR_NAME} ||
        {
            echo "Failed to create ${HDF5_DIR_NAME}"
            exit
        }

    HDF5_MAJ_VER=1.10
    HDF5_MIN_VER=10
    HDF5_VER=${HDF5_MAJ_VER}.${HDF5_MIN_VER}

    echo "Downloading HDF5"
    curl -kLJO \
        https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_MAJ_VER}/hdf5-${HDF5_VER}/src/hdf5-${HDF5_VER}.tar.gz ||
        {
            echo "Failed to download hdf5"
            exit
        }
    tar -xf hdf5-${HDF5_VER}.tar.gz

    cd hdf5-${HDF5_VER} || exit
    make clean
    ./configure --prefix="${HDF5_INSTALL_PATH}" --enable-cxx --enable-fortran --enable-build-mode=production
    make install -j ${NUM_COMPILE_CORES}
    if [ $? -eq 2 ]; then
        echo "HDF5 Build failed"
        exit
    fi
    echo "HDF5 built"
}

function download_superlu_dist() {
    cd "$WORKDIR" || exit
    echo "Downloading SuperLU_dist"
    curl -kLJO https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v8.1.0.tar.gz
}

function build_petsc() {
    cd "$WORKDIR" || exit
    mkdir -p petsc
    cd petsc || exit
    curl -kL -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.3.tar.gz
    tar -xf petsc-3.19.3.tar.gz -C .
    cd petsc-3.19.3 || exit
    ./configure \
        --with-cc=$CC --with-cxx=$CXX --with-fc=$FC -CXXPP=cpp \
        --prefix="${WORKDIR}"/${PETSC_DIR_NAME} \
        --download-hypre=1 \
        --with-shared-libraries \
        --with-debugging=no \
        --with-hdf5-dir="${WORKDIR}"/${HDF5_DIR_NAME} \
        --with-blaslapack-dir="${MKLROOT}" \
        --download-metis=1 \
        --download-parmetis=1 \
        --download-ptscotch=1 \
        --download-mumps=1 \
        --download-superlu_dist="${WORKDIR}"/superlu_dist-8.1.0.tar.gz \
        --download-scalapack=1 \
        --download-slepc=1 \
        --with-mpi=1 \
        --with-cxx-dialect=C++17 \
        --with-fortran-bindings=0 \
        --with-sowing=0 \
        --with-64-bit-indices \
        --with-make-np="${SLURM_NTASKS}" \
        COPTFLAGS='-O3 -fno-slp-vectorize' \
        CXXOPTFLAGS='-O3 -fno-slp-vectorize' \
        FOPTFLAGS='-O3 -fno-slp-vectorize' \
        PETSC_DIR="$(pwd)" PETSC_ARCH=arch-linux-c-opt
    make
    make PETSC_DIR="${WORKDIR}"/${PETSC_DIR_NAME}/petsc-3.19.3 PETSC_ARCH=arch-linux-c-opt install ||
        {
            echo "Failed to build petsc"
            exit
        }
    cd ../..
    export PETSC_DIR=$WORKDIR/petsc
}

function build_libtirpc() {
    # Build libtirpc (libmesh dependency)
    wget https://sourceforge.net/projects/libtirpc/files/libtirpc/1.3.3/libtirpc-1.3.3.tar.bz2 -O- | tar -xj
    mkdir libtirpc
    mv libtirpc-* libtirpc/build
    cd libtirpc/build || exit
    ./configure --prefix="$WORKDIR"/libtirpc --disable-gssapi --disable-static && make && make install
}

function build_moose() {
    export MOOSE_JOBS=32
    cd "$WORKDIR" || exit
    git clone https://github.com/idaholab/moose
    cd moose || exit
    git checkout ${MOOSE_COMMIT} ||
        {
            echo "Checkout failed"
            exit
        }

    if [ ! -f "$WORKDIR/petsc/lib/libpetsc.so" ]; then
        echo "PETSc Install Unsuccessful"
        return
    fi

    export PETSC_DIR=$WORKDIR/petsc
    export PETSC_ARCH=arch-linux-c-opt

    #libmesh configure fails with the llvm compilers so we set these back here
    export I_MPI_CXX=icpc
    export I_MPI_FC=ifort
    export I_MPI_F90=ifort
    export I_MPI_F77=ifort
    export I_MPI_C=icc

    export CPPFLAGS=-I$WORKDIR/libtirpc/include/tirpc
    export LDFLAGS=-L$WORKDIR/libtirpc/lib

    ./scripts/update_and_rebuild_libmesh.sh \
        --with-cxx-std=2017 \
        --with-cc=$CC \
        --with-cxx=$CXX \
        --with-fc=$FC

    export I_MPI_CXX=icpx
    export I_MPI_FC=ifx
    export I_MPI_F90=ifx
    export I_MPI_F77=ifx
    export I_MPI_C=icx

    ./configure --with-derivative-size=200 --with-ad-indexing-type=global
    METHODS='opt' ./scripts/update_and_rebuild_wasp.sh
    cd framework || exit
    METHOD=opt make -j$NUM_COMPILE_CORES
    cd ../modules || exit
    METHOD=opt make -j$NUM_COMPILE_CORES
    cd ../test || exit
    METHOD=opt make -j$NUM_COMPILE_CORES
    ./run_tests -j$NUM_COMPILE_CORES
    cd ../..
}

function build_gslib() {
    cd "$WORKDIR" || exit
    if [ -d "$WORKDIR/gslb" ]; then
        return
    fi
    git clone https://github.com/Nek5000/gslib.git
    cd gslib || exit
    make CFLAGS='-O2 -fPIC'
}

function build_mfem() {
    cd "$WORKDIR" || exit
    git clone https://github.com/Heinrich-BR/mfem.git
    cd mfem || exit
    git checkout SubmeshBoundary
    sed -i "s | list | # list|g" "$WORKDIR"/mfem/config/cmake/modules/FindNetCDF.cmake
    mkdir build
    cd build || exit
    echo "Building MFEM"
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=YES \
        -DMFEM_USE_OPENMP=NO \
        -DMFEM_THREAD_SAFE=NO \
        -DHYPRE_DIR="$WORKDIR"/petsc/ \
        -DMFEM_USE_LAPACK=YES \
        -DMFEM_USE_MPI=YES \
        -DMFEM_USE_METIS_5=YES \
        -DMETIS_DIR="$WORKDIR"/petsc/ \
        -DParMETIS_DIR="$WORKDIR"/petsc/ \
        -DMFEM_USE_SUPERLU=YES \
        -DSuperLUDist_DIR="$WORKDIR"/petsc/ \
        -DSuperLUDist_VERSION_OK=YES \
        -DMFEM_USE_NETCDF=YES \
        -DNETCDF_LIBRARIES="$WORKDIR"/moose/libmesh/installed/lib/libnetcdf.so \
        -DNETCDF_INCLUDE_DIRS="$WORKDIR"/moose/libmesh/contrib/netcdf/netcdf-c-4.6.2/include \
        -DHDF5_DIR="$WORKDIR"/${HDF5_DIR_NAME}/ \
        -DMFEM_USE_GSLIB=YES \
        -DGSLIB_DIR="$WORKDIR"/gslib/build
    make -j$NUM_COMPILE_CORES
    cd miniapps/common || exit
    make -j$NUM_COMPILE_CORES
}

function build_platypus() {
    cd "$WORKDIR" || exit
    git clone https://github.com/aurora-multiphysics/platypus.git
    cd platypus || exit
    git checkout ${PLATYPUS_COMMIT}
    make -j$NUM_COMPILE_CORES
}

load_modules
build_hdf5
download_superlu_dist
build_petsc
build_libtirpc
build_moose
build_gslib
build_mfem
build_platypus

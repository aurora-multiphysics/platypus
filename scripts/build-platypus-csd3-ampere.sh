#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mail-type=none
#SBATCH -p ampere
#SBATCH -A ukaea-ap001-GPU
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --output=platypus_gpu_build.%j.out
#SBATCH --error=platypus_gpu_build.%j.err

## WARNING: THIS SCRIPT WILL UNINSTALL ALL SPACK MODULES ASSOCIATED WITH
## THE ARCHITECTURE DEFINED IN THE ARCH VARIABLE. IF YOU DO NOT WISH TO DO
## THAT, COMMENT OUT THE SPACK UNINSTALL LINE BEFORE SUBMITTING THE SCRIPT
## -> UNINSTALL LINE IN THE install_spack_deps() FUNCTION
ARCH="linux-rocky8-zen"

export compile_cores=32

load_modules() {

    # Load modules

    # shellcheck source=/dev/null
    . /etc/profile.d/modules.sh                # Leave this line (enables the module command)
    module purge
    module load rhel8/slurm
    module use /usr/local/software/spack/spack-modules/rocky8-a100-20230831/linux-rocky8-zen3
    module load cuda/11.7.1

}

set_paths() {

    USER=$(whoami)
    BUILD_PREFIX=platypus_gpu
    BUILD_DIR_NAME=${BUILD_PREFIX}_build

    ROOT_PATH=/home/${USER}/rds/rds-ukaea-ap001/${USER}
    BUILD_PATH=${ROOT_PATH}/${BUILD_DIR_NAME}

    echo "Building in ${BUILD_PATH}"
    mkdir -p "${BUILD_PATH}" || {
                                  echo "Failed to create ${BUILD_PATH}"
                                                                          exit 1
    }

    cd "${BUILD_PATH}" || exit 1

}

check_spack() {

    cd "${ROOT_PATH}" || exit 1

    if [ "$(command -v spack)" ]; then
        echo "Spack command detected. Using pre-loaded spack."
    elif [ -f "${ROOT_PATH}"/spack/share/spack/setup-env.sh ]; then
        echo "Spack detected in root directory. Loading."
        # shellcheck source=/dev/null
        . spack/share/spack/setup-env.sh
    else
        echo "No spack detected. Building from source."
        git clone --depth=100 https://github.com/spack/spack.git
        # shellcheck source=/dev/null
        . spack/share/spack/setup-env.sh
    fi

}

install_spack_deps() {

    # Cleaning up everything to start with a new environment
    spack uninstall -ay arch=${ARCH}

    spack external find cuda@11.7.1

    echo "Installing Petsc..."
    # Spack's petsc doesn't like openmpi, but it works with mpich
    spack install petsc +cuda cuda_arch=80 +fortran +hdf5 +hypre +metis +mpi \
        ^mpich +cuda cuda_arch=80 \
        ^hdf5 +cxx +fortran +hl +mpi +shared \
        ^hypre +mpi +shared +cuda cuda_arch=80 +superlu-dist +cublas +gpu-aware-mpi \
        ^superlu-dist +cuda cuda_arch=80 +parmetis +shared
    spack load petsc arch=${ARCH}

    echo "Installing SLEPc..."
    spack install slepc +cuda cuda_arch=80
    spack load slepc arch=${ARCH}

    echo "Installing netcdf..."
    spack install netcdf-c +parallel-netcdf
    spack load netcdf-c arch=${ARCH}

    echo "Installing ninja..."
    spack install ninja
    spack load ninja arch=${ARCH}

    echo "Adding python modules..."

    spack install py-pyaml
    spack load py-pyaml arch=${ARCH}

    spack install py-jinja2
    spack load py-jinja2 arch=${ARCH}

    spack install py-packaging
    spack load py-packaging arch=${ARCH}

    spack install py-setuptools
    spack load py-setuptools arch=${ARCH}

}

install_gslib() {

    echo "Installing gslib..."
    cd "${BUILD_PATH}" || exit 1
    git clone https://github.com/Nek5000/gslib.git
    cd gslib || exit 1
    make CC=mpicc CFLAGS="-O2 -fPIC" -j"$compile_cores"
}

install_mfem() {

    export CXX=mpic++
    export CC=mpicc
    export F90=mpif90
    export F77=mpif77
    export FC=mpif90

    # Build MFEM
    cd "${BUILD_PATH}" || exit 1
    git clone https://github.com/mfem/mfem.git
    cd mfem || exit 1
    # This is just until MFEM merges Edward's changes. Without this, GPU build crashes!
    git checkout EdwardPalmer99/add-missing-header-to-exodus-writer-fix
    mkdir build
    cd build || exit 1
    echo "Building MFEM"
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=YES \
        -DMFEM_USE_OPENMP=NO \
        -DMFEM_THREAD_SAFE=YES \
        -DMFEM_ENABLE_EXAMPLES=YES \
        -DMFEM_ENABLE_MINIAPPS=YES \
        -DMFEM_USE_MPI=YES \
        -DMFEM_USE_CUDA=YES \
        -DCUDA_ARCH=sm_80 \
        -DMFEM_USE_METIS_5=YES \
        -DMFEM_USE_SUPERLU=YES \
        -DMFEM_USE_NETCDF=YES \
        -DMFEM_USE_GSLIB=YES \
        -DGSLIB_DIR="${BUILD_PATH}/gslib/build"

    if [ $? -eq 2 ]; then
        echo "MFEM config failed"
        exit 1
    fi

    make -j"$compile_cores"

    if [ $? -eq 2 ]; then
        echo "MFEM build failed"
        exit 1
    fi

    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BUILD_PATH}/mfem/build:${BUILD_PATH}/mfem/build/miniapps/common
}

install_moose() {

    # Some of the variables needed
    export MOOSE_JOBS=$compile_cores
    export LIBMESH_JOBS=$compile_cores
    export METHOD="opt"
    SLEPC_DIR=$(spack find --format "{prefix}" slepc arch=${ARCH})
    export SLEPC_DIR

    cd "${BUILD_PATH}" || exit 1
    git clone https://github.com/idaholab/moose
    cd moose || exit 1

    echo "Building libmesh..."
    ./scripts/update_and_rebuild_libmesh.sh --with-mpi
    if [ $? -eq 2 ]; then
        echo "libmesh build failed"
        exit 1
    fi

    echo "Building WASP..."
    ./scripts/update_and_rebuild_wasp.sh
    if [ $? -eq 2 ]; then
        echo "WASP build failed"
        exit 1
    fi

    ./configure --with-derivative-size=200
    if [ $? -eq 2 ]; then
        echo "MOOSE configure failed"
        exit 1
    fi

    cd framework || exit 1
    make -j"$compile_cores"
    if [ $? -eq 2 ]; then
        echo "MOOSE framework build failed"
        exit 1
    fi

    cd ../modules || exit 1
    make -j"$compile_cores"
    if [ $? -eq 2 ]; then
        echo "MOOSE modules build failed"
        exit 1
    fi

    # This takes very long! Only run the tests if you really need to!
    #cd ../test || exit 1
    #make -j"$compile_cores"
    #if [ $? -eq 2 ]; then
    #    echo "MOOSE test build failed"
    #    exit 1
    #fi

    #./run_tests -j"$compile_cores"
}

install_platypus() {

    cd "${BUILD_PATH}" || exit 1

    echo "Building platypus..."
    git clone https://github.com/aurora-multiphysics/platypus.git
    cd platypus || exit 1
    git submodule update --init --recursive
    cd contrib/hephaestus/ || exit 1
    mkdir build
    cd build || exit 1
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR="${BUILD_PATH}/mfem/build" ..
    ninja
    cd "${BUILD_PATH}"/platypus || exit 1
    make -j"$compile_cores"

}

load_modules
set_paths
check_spack
install_spack_deps
install_gslib
install_mfem
install_moose
install_platypus
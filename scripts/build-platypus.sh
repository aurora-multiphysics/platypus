#!/bin/bash

help() {
    printf "\n This bash script utilises spack to build and install Platypus and its dependencies, including MOOSE and MFEM.\n"
    printf "\nUsage: ./build-platypus [options]\n"
    printf "\n Available options:\n"
    printf "\n -h, --help\n Displays this help message.\n"
    printf "\n -cl, --clean\n Deletes the entire build directory along with all dependencies and installed files.\n"
    printf "\n -g, --gpu\n Defines a GPU build. If this option is not added, a CPU build is assumed.\n"
    printf "\n -b=[...], --gpu-backend=[...]\n Defines the GPU backend to be used. Options are cuda and rocm for NVIDIA and AMD GPUs, respectively.\n"
    printf "\n -a=[...], --gpu-arch=[...]\n Defines the target GPU architecture. For CUDA backends, use only the number. For instance, to target a GPU whose arch code is sm_80, you would add -a=80.\n"
    printf "\n -mpicxx=[<path>], --ompi-cxx=[<path>]\n Path to a C++ compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds.\n"
    printf "\n -mpicc=[<path>], --ompi-cc=[<path>]\n Path to a C compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds.\n"
    printf "\n -p=[<name> <version> <path>], --package=[<name> <version> <path>]\n Adds an external package to the spack environment so that it is not built by spack. It is possible to add any number of packages.\n"
    printf "\n -c=<path>, --compiler=<path>\n Adds an external compiler to the spack environment to be used to build the spack packages. \n"
    printf "\n Example usage:\n"

    printf "\n ./scripts/build-platypus -g \\ \n"
    printf "                          -b=rocm \\ \n"
    printf "                          -a=gfx942 \\ \n"
    printf "                          -mpicxx=/opt/rocm/bin/amdclang++ \\ \n"
    printf "                          -mpicc=/opt/rocm-6.2.4/bin/amdclang \\ \n"
    printf "                          -p=\"hip 6.2.4 /opt/rocm-6.2.4/\" \\ \n"
    printf "                          -p=\"rocrand 6.2.4 /opt/rocm-6.2.4/\" \\ \n"
    printf "                          -c=\"/opt/llvm\" \n\n"
    exit 0
}

clean() {
    rm -rf "${BUILD_PATH}"
    exit 0
}

replace_in_file() {
    # First argument is the file
    # Second argument is the word to be replaced
    # Third argument is what it is to be replaced by
    sed -i "s/@$2@/$3/" "$1"
}

parse_options() {
    for arg in "$@"; do
        case $arg in
            -h | --help)
            help
            ;;
            -cl | --clean)
            clean
            ;;
            -g | --gpu)
            GPU_BUILD=1
            ;;
            -b=* | --gpu-backend=*)
            GPU_BACKEND="${arg#*=}"
            ;;
            -a=* | --gpu-arch=*)
            GPU_ARCH="${arg#*=}"
            ;;
            -ompi-v=* | --openmpi-version=*)
            OPENMPI_VER="${arg#*=}"
            ;;
            -llvm-v=* | --llvm-version=*)
            LLVM_VER="${arg#*=}"
            ;;
            -amdllvm-v=* | --llvm-amdgpu-version=*)
            AMDLLVM_VER="${arg#*=}"
            ;;
            -p=* | --package=*)
            PACKAGES+=("${arg#*=}")
            ;;
            -c=* | --compiler=*)
            SPACK_COMPILER_PATH="${arg#*=}"
            ;;
            -mpicxx=* | --ompi-cxx=*)
            OMPICXX="${arg#*=}"
            ;;
            -mpicc=* | --ompi-cc=*)
            OMPICC="${arg#*=}"
            ;;
            *)
            OTHER_ARGUMENTS+=("${arg}")
            ;;
        esac
    done

    export_config_file "$@"
}

export_config_file() {

    INV_COMMAND="$*"
    {
        printf 'Invocation command:\n\n ./%s %s\n\n' "$(basename "$0")" "${INV_COMMAND}"
        printf 'Unused arguments:\n\n'
        for u in "${OTHER_ARGUMENTS[@]}"; do
            printf '%s\n' "${u}"
        done

        printf '\nOptions:\n\n'
        printf 'GPU_BUILD = %s\n' "${GPU_BUILD}"
        printf 'GPU_BACKEND = %s\n' "${GPU_BACKEND}"
        printf 'GPU_ARCH = %s\n' "${GPU_ARCH}"
        printf 'OMPI_CXX = %s\n' "${OMPICXX}"
        printf 'OMPI_CC = %s\n\n' "${OMPICC}"
    } >> "${BUILD_PATH}"/${CONFIG_FILE}

}

load_spack() {
    git clone --depth=100 https://github.com/spack/spack.git
    # shellcheck source=/dev/null
    . spack/share/spack/setup-env.sh
}

make_spack_env() {

    if [ "${GPU_BUILD}" -eq 1 ]; then
        printf "GPU build detected\n"
        replace_in_file ${SPACK_MOD} "gpu_aware_mpi" "+gpu-aware-mpi"

        if [[ -z ${GPU_BACKEND} || (${GPU_BACKEND} != "cuda" && ${GPU_BACKEND} != "rocm")         ]]; then
            printf "Please set the GPU backend with --gpu-backend=[...]. Options are cuda and rocm.\n"
            exit 1
        else
            printf 'GPU backend %s detected\n' "${GPU_BACKEND}"
            if [ "${GPU_BACKEND}" = "cuda" ]; then
                export LLVM_TYPE="llvm"
                replace_in_file ${SPACK_MOD} "llvm_version" "@${LLVM_VER}"
                replace_in_file ${SPACK_MOD} "openmpi" "openmpi@openmpi_version@ @gpu@ @gpu_arch@"
                replace_in_file ${SPACK_MOD} "ucx" "ucx @gpu@ @gpu_arch@ +gdrcopy"
            else
                export LLVM_TYPE="llvm-amdgpu"
                replace_in_file ${SPACK_MOD} "llvm_version" "@${AMDLLVM_VER}"
                replace_in_file ${SPACK_MOD} "openmpi" "openmpi@openmpi_version@"
                replace_in_file ${SPACK_MOD} "ucx" "ucx @gpu@"
            fi
            replace_in_file ${SPACK_MOD} "gpu" "+${GPU_BACKEND}"
        fi

        if [ -z "${GPU_ARCH}" ]; then
            printf "GPU arch target not detected. Not set\n"
            replace_in_file ${SPACK_MOD} "gpu_arch" ""
        else
            printf 'GPU arch target %s detected\n' "${GPU_ARCH}"
            if [ "${GPU_BACKEND}" = "cuda" ]; then
                replace_in_file ${SPACK_MOD} "gpu_arch" "cuda_arch=${GPU_ARCH}"
            else
                replace_in_file ${SPACK_MOD} "gpu_arch" "amdgpu_target=${GPU_ARCH}"
            fi
        fi
    else
        printf "CPU build detected\n"

        export LLVM_TYPE="llvm"

        # Clean up all GPU options
        replace_in_file ${SPACK_MOD} "openmpi" "openmpi@openmpi_version@"
        replace_in_file ${SPACK_MOD} "ucx" "ucx"
        replace_in_file ${SPACK_MOD} "gpu_aware_mpi" ""
        replace_in_file ${SPACK_MOD} "gpu" ""
        replace_in_file ${SPACK_MOD} "gpu_arch" ""
        replace_in_file ${SPACK_MOD} "blas" ""
        replace_in_file ${SPACK_MOD} "llvm_version" "@${LLVM_VER}"
    fi

    replace_in_file ${SPACK_MOD} "openmpi_version" "@${OPENMPI_VER}"
    replace_in_file ${SPACK_MOD} "llvm" ${LLVM_TYPE}

}

add_package() {

    # First argument is the package name
    # Second argument is the version
    # Third argument is the path

    if grep -q "$1@$2" "${SPACK_MOD}"; then
        echo "$1 module found in spack environment"
    else
        printf '    %s:\n      externals:\n      - spec: %s\n        prefix: %s\n      buildable: False\n' "$1" "$1@$2" "$3" >> "${SPACK_MOD}"
    fi
}

add_external_packages() {
    if [ -z "${PACKAGES[*]}" ]; then
        printf "No external packages added\n"
    else
        printf "  packages:\n"          >> ${SPACK_MOD}
        printf 'External packages:\n' >> ${CONFIG_FILE}
        for p in "${PACKAGES[@]}"; do
            read -ra STR_ARRAY <<< "$p"
            printf "\nExternal package added"
            printf '\nName: %s\n' "${STR_ARRAY[0]}"  | tee -a ${CONFIG_FILE}
            printf 'Version: %s\n' "${STR_ARRAY[1]}" | tee -a ${CONFIG_FILE}
            printf 'Path: %s\n' "${STR_ARRAY[2]}"    | tee -a ${CONFIG_FILE}
            add_package "${STR_ARRAY[0]}" "${STR_ARRAY[1]}" "${STR_ARRAY[2]}"

            case ${STR_ARRAY[0]} in
                openmpi)
                export OPENMPI_VER=${STR_ARRAY[1]}
                ;;
                llvm)
                export LLVM_VER=${STR_ARRAY[1]}
                ;;
                llvm-amdgpu)
                export AMDLLVM_VER=${STR_ARRAY[1]}
                ;;
            esac

        done
    fi

}

set_environment_vars() {
    SLU_DIR=$(spack location -i superlu-dist)
    HDF5_DIR=$(spack location -i hdf5)
    SLEPC_DIR=$(spack location -i slepc)
    PETSC_DIR=$(spack location -i petsc)
    CONDUIT_DIR=$(spack location -i conduit)
    TIRPC_DIR=$(spack location -i libtirpc)
    CEED_DIR=$(spack location -i libceed)

    export HDF5_DIR
    export SLEPC_DIR
    export PETSC_DIR
    export CONDUIT_DIR
    export CEED_DIR

    export CPPFLAGS="${CPPFLAGS} -I${TIRPC_DIR}/include/tirpc"
    export LDFLAGS="${LDFLAGS} -L${TIRPC_DIR}/lib"
    export CXXFLAGS="${CXXFLAGS} -D_GLIBCXX_USE_CXX11_ABI=1"

    if [ -z "${OMPICXX}" ]; then
        OMPI_CXX=$(spack location -i ${LLVM_TYPE})/bin/clang++
    else
        OMPI_CXX=${OMPICXX}
    fi

    if [ -z "${OMPICC}" ]; then
        OMPI_CC=$(spack location -i ${LLVM_TYPE})/bin/clang
    else
        OMPI_CC=${OMPICC}
    fi

    export OMPI_CXX
    export OMPI_CC

    export compile_cores=16
    export CXX=mpic++
    export CC=mpicc
    export F90=mpif90
    export F77=mpif77
    export FC=mpif90
    export MOOSE_DIR=${BUILD_PATH}/moose
    export MFEM_DIR=${BUILD_PATH}/mfem/installed
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MFEM_DIR}/lib
    export MOOSE_JOBS=$compile_cores
    export LIBMESH_JOBS=$compile_cores
    export METHOD="dbg"

    if [ "${GPU_BACKEND}" = "cuda" ]; then
        export CUDA_MFEM="YES"
    else
        export CUDA_MFEM="NO"
    fi

    if [ "${GPU_BACKEND}" = "rocm" ]; then
        export HIP_MFEM="YES"
    else
        export HIP_MFEM="NO"
    fi

}

install_mfem() {
    cd "${BUILD_PATH}" || exit 1
    git clone https://github.com/mfem/mfem.git
    cd mfem || exit 1
    git checkout master
    spack load cmake
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_STANDARD=14 \
        -DCMAKE_INSTALL_PREFIX="${BUILD_PATH}"/mfem/installed \
        -DBUILD_SHARED_LIBS=YES \
        -DMFEM_USE_OPENMP=NO \
        -DMFEM_THREAD_SAFE=NO \
        -DMFEM_ENABLE_EXAMPLES=YES \
        -DMFEM_ENABLE_MINIAPPS=YES \
        -DMFEM_USE_MPI=YES \
        -DMFEM_USE_CUDA="${CUDA_MFEM}" \
        -DMFEM_USE_HIP="${HIP_MFEM}" \
        -DCUDA_ARCH=sm_"${GPU_ARCH}" \
        -DHIP_ARCH="${GPU_ARCH}" \
        -DMFEM_USE_METIS_5=YES \
        -DMFEM_USE_SUPERLU=YES \
        -DMFEM_USE_NETCDF=YES \
        -DMFEM_USE_GSLIB=YES \
        -DMFEM_USE_CONDUIT=YES \
        -DMFEM_USE_CEED=YES \
        -DGSLIB_DIR="${GSLIB_DIR}" \
        -DCONDUIT_DIR="${CONDUIT_DIR}" \
        -DHDF5_DIR="${HDF5_DIR}" \
        -DCEED_DIR="${CEED_DIR}" \
        -DSuperLUDist_DIR="${SLU_DIR}" \
        -DSuperLUDist_VERSION_OK=YES \
        -DHYPRE_VERSION=23200
    cd build || exit
    make install -j $compile_cores
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX="${BUILD_PATH}"/mfem/installed -DBUILD_SHARED_LIBS=YES -DMFEM_USE_OPENMP=NO -DMFEM_THREAD_SAFE=NO -DMFEM_ENABLE_EXAMPLES=YES -DMFEM_ENABLE_MINIAPPS=YES -DMFEM_USE_MPI=YES -DMFEM_USE_CUDA=NO -DMFEM_USE_METIS_5=YES -DMFEM_USE_SUPERLU=YES -DMFEM_USE_NETCDF=YES -DMFEM_USE_GSLIB=YES -DMFEM_USE_CONDUIT=YES -DGSLIB_DIR="${GSLIB_DIR}" -DCONDUIT_DIR="${CONDUIT_DIR}" -DHDF5_DIR="${HDF5_DIR}" -DCEED_DIR="${CEED_DIR}" -DSuperLUDist_DIR="${SLU_DIR}" -DSuperLUDist_VERSION_OK=YES -DHYPRE_VERSION=23200
}

install_axom() {
    cd ${BUILD_PATH} || exit 1
    git clone --recursive https://github.com/LLNL/axom.git axom-repo
    cd axom-repo || exit 1

    # the mfem repo has a cmake file that we can use
    python3 ./config-build.py -hc ${BUILD_PATH}/mfem/miniapps/tribol/axom-gcc-notpl.cmake -bt Debug -DCMAKE_INSTALL_PREFIX=${BUILD_PATH}/axom -DAXOM_ENABLE_MIR=NO -DAXOM_ENABLE_SINA=NO || exit 1
    cd build-axom-gcc-notpl-debug || exit 1
    make -j $compile_cores install || exit 1
}

install_tribol() {
    # must be done strictly after configuring mfem so that we have something
    # to link to. Then we go back and make sure that mfem is configured with
    # tribol.
    cd ${BUILD_PATH} || exit 1
    # use my fork to evade the axom namespace issue
    # cf https://github.com/LLNL/Tribol/issues/133
    git clone https://github.com/sean-baccas/Tribol --recursive tribol-repo
    cd tribol-repo || exit 1
    python3 ./config-build.py -hc ${BUILD_PATH}/mfem/miniapps/tribol/tribol-gcc-basictpl.cmake -bt Debug -DCMAKE_INSTALL_PREFIX=${BUILD_PATH}/tribol -DMFEM_DIR=${BUILD_PATH}/mfem/build/ || exit 1
    cd build-tribol-gcc-basictpl-debug || exit 1
    make -j $compile_cores install || exit 1

    # finally reconfigure mfem
    cd ${BUILD_PATH}/mfem
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_CXX_STANDARD=14 \
        -DCMAKE_INSTALL_PREFIX="${BUILD_PATH}"/mfem/installed \
        -DBUILD_SHARED_LIBS=YES \
        -DMFEM_USE_OPENMP=NO \
        -DMFEM_THREAD_SAFE=YES \
        -DMFEM_ENABLE_EXAMPLES=YES \
        -DMFEM_ENABLE_MINIAPPS=YES \
        -DMFEM_USE_MPI=YES \
        -DMFEM_USE_CUDA="${CUDA_MFEM}" \
        -DMFEM_USE_HIP="${HIP_MFEM}" \
        -DCUDA_ARCH=sm_"${GPU_ARCH}" \
        -DHIP_ARCH="${GPU_ARCH}" \
        -DMFEM_USE_METIS_5=YES \
        -DMFEM_USE_SUPERLU=YES \
        -DMFEM_USE_NETCDF=YES \
        -DMFEM_USE_GSLIB=YES \
        -DMFEM_USE_CONDUIT=YES \
        -DMFEM_USE_CEED=YES \
        -DMFEM_USE_TRIBOL=YES \
        -DGSLIB_DIR="${GSLIB_DIR}" \
        -DCONDUIT_DIR="${CONDUIT_DIR}" \
        -DHDF5_DIR="${HDF5_DIR}" \
        -DCEED_DIR="${CEED_DIR}" \
        -DSuperLUDist_DIR="${SLU_DIR}" \
        -DSuperLUDist_VERSION_OK=YES \
        -DHYPRE_VERSION=23200
    cd build
    make install -j $compile_cores
}


install_moose() {

    spack load py-deepdiff
    spack load py-jinja2
    spack load py-packaging
    spack load py-pyyaml
    spack load py-setuptools
    spack load py-xmltodict
    spack load py-pip

    cd "${BUILD_PATH}" || exit 1
    #git clone https://github.com/idaholab/moose
    git clone https://github.com/Heinrich-BR/moose.git
    cd moose || exit 1
    #git checkout master
    ./scripts/update_and_rebuild_libmesh.sh --with-mpi
    ./scripts/update_and_rebuild_wasp.sh

    ./configure --with-derivative-size=200
    cd framework || exit 1
    make -j"$compile_cores"
    cd ../modules || exit 1
    make -j"$compile_cores"
}

install_platypus() {
    cd "${ROOT_PATH}" || exit 1

    echo "Building platypus..."
    make -j"$compile_cores"
    ./run_tests -j"$compile_cores"
}

# Template file for the spack environment
SPACK_FILE="scripts/spack-env.txt"

# Name of the file to be used for spack environment
SPACK_MOD=".spack_env_platypus.yaml"

# Name of the config file where we print the invocation options
CONFIG_FILE="build_platypus_config.txt"

# Location for the spack hidden directory
export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CACHE_PATH="deps"

GPU_BUILD=0
GPU_BACKEND=""
GPU_ARCH=""
LLVM_VER="18.1.8"
AMDLLVM_VER="6.2.4"
OPENMPI_VER="5.0.6"
OMPICXX=""
OMPICC=""
SPACK_COMPILER_PATH=""
PACKAGES=()
OTHER_ARGUMENTS=()

export BUILD_DIR_NAME="deps"
ROOT_PATH=$(pwd)
export ROOT_PATH
export BUILD_PATH=${ROOT_PATH}/${BUILD_DIR_NAME}

parse_options "$@"

mkdir -p "${BUILD_PATH}"

# Create modifiable spack environment file
cp ${SPACK_FILE} "${BUILD_PATH}"/${SPACK_MOD}

cd "${BUILD_PATH}" || exit 1

add_external_packages
load_spack
make_spack_env

# If an external compiler is not added from the command line, spack will try to find existing compilers on the machine
if [ -z "${SPACK_COMPILER_PATH}" ]; then
    printf "Spack compiler not set. Automatically finding compilers...\n"
    spack compiler find
else
    spack compiler add "${SPACK_COMPILER_PATH}"
fi

spack install bzip2
spack load bzip2

printf "Creating spack environment\n"
spack env create platypus ${SPACK_MOD}
spack env activate platypus

spack concretize -f

# Using vim instead of emacs to avoid a parsing bug in the autoconf installation
EMACS=vim spack install
# spack load petsc

# # Cleaning intermediary files to use less space
# # spack clean -a

# set_environment_vars
# install_mfem
# install_axom
# install_tribol
# install_moose
# install_platypus

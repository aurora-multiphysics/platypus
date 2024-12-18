#!/bin/bash

help() {
    printf "\n This bash script utilises spack to build and install Platypus and its dependencies, including MOOSE and MFEM.\n"
    printf "\nUsage: ./build-platypus [options]\n"
    printf "\n Available options:\n"
    printf "\n -h, --help\n Displays this help message.\n"
    printf "\n -g, --gpu\n Defines a GPU build. If this option is not added, a CPU build is assumed.\n"
    printf "\n -b=[...], --gpu-backend=[...]\n Defines the GPU backend to be used. Options are cuda and rocm for NVIDIA and AMD GPUs, respectively.\n"
    printf "\n -a=[...], --gpu-arch=[...]\n Defines the target GPU architecture. For CUDA backends, use only the number. For instance, to target a GPU whose arch code is sm_80, you would add -a=80\n"
    printf "\n -mpicxx=[<path>], --ompi-cxx=[<path>]\n Path to a C++ compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds.\n"
    printf "\n -mpicc=[<path>], --ompi-cc=[<path>]\n Path to a C compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds.\n"
    printf "\n -p=[<name> <version> <path>], --package=[<name> <version> <path>]\n Adds an external package to the spack environment so that it is not built by spack. It is possible to add any number of packages.\n"
    printf "\n -c=[<name> <version> <options>], --compiler=[<name> <version> <options>]\n Adds an external compiler to the spack environment. It is possible to add any number of compilers. In <options>, one would include the path to CC, CXX, F77, and FC compilers. It is not necessary to fill them all. See example below.\n"
    printf "\n Example usage:\n"

    printf "\n ./build-platypus -g \\ \n"
    printf "                  -b=rocm \\ \n"
    printf "                  -a=gfx942 \\ \n"
    printf "                  -mpicxx=/opt/rocm/bin/amdclang++ \\ \n"
    printf "                  -mpicc=/opt/rocm-6.2.4/bin/amdclang \\ \n"
    printf "                  -p=\"hip 6.2.4 /opt/rocm-6.2.4/\" \\ \n"
    printf "                  -p=\"rocrand 6.2.4 /opt/rocm-6.2.4/\" \\ \n"
    printf "                  -c=\"clang 16.0.0 CXX=/opt/llvm/bin/clang++ CC=/opt/llvm/bin/clang F77=/opt/llvm/bin/flang\" \n\n"
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
            -g | --gpu)
            GPU_BUILD=1
            ;;
            -b=* | --gpu-backend=*)
            GPU_BACKEND="${arg#*=}"
            ;;
            -a=* | --gpu-arch=*)
            GPU_ARCH="${arg#*=}"
            ;;
            -p=* | --package=*)
            PACKAGES+=("${arg#*=}")
            ;;
            -c=* | --compiler=*)
            COMPILERS+=("${arg#*=}")
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
    } >> ${CONFIG_FILE}

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
            replace_in_file ${SPACK_MOD} "gpu" "+${GPU_BACKEND}"
            if [ "${GPU_BACKEND}" = "cuda" ]; then
                replace_in_file ${SPACK_MOD} "blas" "+cublas"
                replace_in_file ${SPACK_MOD} "amdgpu" ""
                replace_in_file ${SPACK_MOD} "llvm_version" "18.1.8"
            else
                replace_in_file ${SPACK_MOD} "blas" "+rocblas"
                replace_in_file ${SPACK_MOD} "amdgpu" "-amdgpu"
                replace_in_file ${SPACK_MOD} "llvm_version" "6.2.4"
            fi
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

        # Clean up all GPU options
        replace_in_file ${SPACK_MOD} "gpu_aware_mpi" ""
        replace_in_file ${SPACK_MOD} "gpu" ""
        replace_in_file ${SPACK_MOD} "gpu_arch" ""
        replace_in_file ${SPACK_MOD} "blas" ""
        replace_in_file ${SPACK_MOD} "amdgpu" ""
        replace_in_file ${SPACK_MOD} "llvm_version" "18.1.8"
    fi

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
        done
    fi

}

parse_compiler_options() {
    CC_PATH=""
    CXX_PATH=""
    F77_PATH=""
    FC_PATH=""

    for arg in "$@"; do
        case $arg in
            cc=* | CC=*)
            CC_PATH="${arg#*=}"
            ;;
            cxx=* | CXX=*)
            CXX_PATH="${arg#*=}"
            ;;
            f77=* | F77=*)
            F77_PATH="${arg#*=}"
            ;;
            fc=* | FC=*)
            FC_PATH=("${arg#*=}")
            ;;
            *)
            OTHER_ARGUMENTS+=("${arg}")
            ;;
        esac
    done

    COMP_ARRAY=("${CC_PATH[@]}" "${CXX_PATH[@]}" "${F77_PATH[@]}" "${FC_PATH[@]}")
    for arg in "${COMP_ARRAY[@]}"; do
        if [ -z "$arg" ]; then
            arg="None"
        fi
    done
}

add_compiler() {
    # First argument is the package name
    # Second argument is the version
    # Third argument is the CC path
    # Fourth argument is the CXX path
    # Fifth argument is the F77 path
    # Sixth argument is the FC path

    if grep -q "$1@$2" "${SPACK_MOD}"; then
        printf '%s compiler found in spack environment' "$1"
    else
        printf '  - compiler:\n      spec: %s\n      operating_system: %s\n      modules: []\n      paths:\n        cc: %s\n        cxx: %s\n        f77: %s\n        fc: %s\n' "$1@=$2" "$(spack arch -o)" "$3" "$4" "$5" "$6" >> "${SPACK_MOD}"
    fi
}

add_external_compilers() {
    if [ -z "${COMPILERS[*]}" ]; then
        printf "No external compilers added\n"
    else
        printf "  compilers:\n" >> ${SPACK_MOD}
        printf '\nExternal compilers:\n' >> ${CONFIG_FILE}
        for c in "${COMPILERS[@]}"; do
            read -ra STR_ARRAY <<< "$c"
            parse_compiler_options "${STR_ARRAY[@]}"
            printf "\nExternal compiler added"
            printf '\nName: %s\n' "${STR_ARRAY[0]}"  | tee -a ${CONFIG_FILE}
            printf 'Version: %s\n' "${STR_ARRAY[1]}" | tee -a ${CONFIG_FILE}
            printf 'CC_PATH: %s\n' "${CC_PATH[*]}"   | tee -a ${CONFIG_FILE}
            printf 'CXX_PATH: %s\n' "${CXX_PATH[*]}" | tee -a ${CONFIG_FILE}
            printf 'F77_PATH: %s\n' "${F77_PATH[*]}" | tee -a ${CONFIG_FILE}
            printf 'FC_PATH: %s\n' "${FC_PATH[*]}"   | tee -a ${CONFIG_FILE}
            add_compiler "${STR_ARRAY[0]}" "${STR_ARRAY[1]}" "${CC_PATH[*]}" "${CXX_PATH[*]}" "${F77_PATH[*]}" "${FC_PATH[*]}"
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

    export HDF5_DIR
    export SLEPC_DIR
    export PETSC_DIR
    export CONDUIT_DIR

    export CPPFLAGS="${CPPFLAGS} -I${TIRPC_DIR}/include/tirpc"
    export LDFLAGS="${LDFLAGS} -L${TIRPC_DIR}/lib"

    if [ -z "${OMPICXX}" ]; then
        OMPI_CXX=clang++
    else
        OMPI_CXX=${OMPICXX}
    fi

    if [ -z "${OMPICC}" ]; then
        OMPI_CC=clang
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
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BUILD_PATH}/mfem/build:${BUILD_PATH}/mfem/build/miniapps/common
    export MOOSE_JOBS=$compile_cores
    export LIBMESH_JOBS=$compile_cores
    export METHOD="opt"

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
        -DCMAKE_BUILD_TYPE=Release \
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
        -DGSLIB_DIR="${GSLIB_DIR}" \
        -DCONDUIT_DIR="${CONDUIT_DIR}" \
        -DHDF5_DIR="${HDF5_DIR}" \
        -DSuperLUDist_DIR="${SLU_DIR}" \
        -DSuperLUDist_VERSION_OK=YES \
        -DHYPRE_VERSION=23200
    cmake --build build -j"$compile_cores"
}

install_moose() {
    cd "${BUILD_PATH}" || exit 1
    git clone https://github.com/idaholab/moose
    cd moose || exit 1
    ./scripts/update_and_rebuild_libmesh.sh --with-mpi
    ./scripts/update_and_rebuild_wasp.sh

    ./configure --with-derivative-size=200
    cd framework || exit 1
    make -j"$compile_cores"
    cd ../modules || exit 1
    make -j"$compile_cores"
}

install_platypus() {
    cd "${BUILD_PATH}" || exit 1

    echo "Building platypus..."
    git clone https://github.com/aurora-multiphysics/platypus.git
    cd platypus || exit 1
    make -j"$compile_cores"
    ./run_tests -j"$compile_cores"
}

# Template file for the spack environment
SPACK_FILE="spack-env.txt"

# Name of the file to be used for spack environment
SPACK_MOD=".spack_env_platypus.yaml"

# Name of the config file where we print the invocation options
CONFIG_FILE="build_platypus_config.txt"

GPU_BUILD=0
GPU_BACKEND=""
GPU_ARCH=""
OMPICXX=""
OMPICC=""
PACKAGES=()
COMPILERS=()
OTHER_ARGUMENTS=()

parse_options "$@"

export BUILD_DIR_NAME="platypus_build"
ROOT_PATH=$(pwd)
export ROOT_PATH
export BUILD_PATH=${ROOT_PATH}/${BUILD_DIR_NAME}
mkdir -p "${BUILD_PATH}"

# Create modifiable spack environment file
cp ${SPACK_FILE} "${BUILD_PATH}"/${SPACK_MOD}

cd "${BUILD_PATH}" || exit 1

add_external_packages
add_external_compilers
load_spack
make_spack_env

# Will try to find a pre-installed compiler. If no compilers are found, you might need to add them manually
spack compiler find

spack install bzip2
spack load bzip2

printf "Creating spack environment\n"
spack env create platypus ${SPACK_MOD}
spack env activate platypus

spack concretize -f

# Using vim instead of emacs to avoid a parsing bug in the autoconf installation
EMACS=vim spack install
spack load petsc

# Cleaning intermediary files to use less space
spack clean -a

set_environment_vars
install_mfem
install_moose
install_platypus

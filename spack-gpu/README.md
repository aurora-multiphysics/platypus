# Build platypus depedencies with CUDA for GPU

First install spack.

```bash
git clone --depth=100 https://github.com/spack/spack.git
```

Activate spack

```bash
. spack/share/spack/setup-env.sh
```

Load gcc

```bash
module load gcc/some-version
```

Install llvm-18 with gcc

```
spack env create llvm
spack compiler find
spack add llvm@18 %gcc
spack concretize -f
spack install
```

Unload gcc

```bash
module unload gcc
```

Add the compiler block specifying llvm compiler to the environment's yaml.

```bash
spack compiler find
```

Copy that compiler block in `$SPACK_ENV/spack.yaml` for adding to other
environment `spack.yaml` later.

Load the modules (to get cuda-aware openmpi, we will use a different llvm later)

```bash
module load llvm/18.1.8
module load cuda/12.6.2
module load openmpi/4.1.6
```

Create and build platypus environment:

```bash
spack env create platypus spack-platypus-gpu.yaml
spack env activate platypus
spack concretize -f
spack install
```

Source env setup script:

```bash
. platypus_env.sh
```

Make the build directory:
```bash
mkdir -p "${BUILD_PATH}"
```

clone and build mfem:

```bash
cd "${BUILD_PATH}"
git clone https://github.com/mfem/mfem.git
cd mfem
git checkout master
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=YES \
    -DMFEM_USE_OPENMP=NO \
    -DMFEM_THREAD_SAFE=YES \
    -DMFEM_ENABLE_EXAMPLES=YES \
    -DMFEM_ENABLE_MINIAPPS=YES \
    -DMFEM_USE_MPI=YES \
    -DMFEM_USE_CUDA=YES \
    -DCUDA_ARCH=sm_90 \
    -DMFEM_USE_METIS_5=YES \
    -DMFEM_USE_SUPERLU=YES \
    -DMFEM_USE_NETCDF=YES \
    -DMFEM_USE_GSLIB=YES \
    -DGSLIB_DIR="${BUILD_PATH}/gslib/build" \
    -DHDF5_DIR=${HDF5_DIR} \
    -DSuperLUDist_DIR="${SLU_DIR}" \
    -DSuperLUDist_VERSION_OK=YES
cmake --build build -j"$compile_cores"
```

build moose:

```bash
cd "${BUILD_PATH}"
git clone https://github.com/idaholab/moose
cd moose
./scripts/update_and_rebuild_libmesh.sh --with-mpi
./scripts/update_and_rebuild_wasp.sh

./configure --with-derivative-size=200
cd framework
make -j"$compile_cores"
cd ../modules
make -j"$compile_cores"
```

build platypus:

```bash
cd "${BUILD_PATH}"
git clone https://github.com/aurora-multiphysics/platypus.git
cd platypus
make -j"$compile_cores"
```

# Set up state again:

```bash
module load llvm/18.1.8
module load cuda/12.6.2
module load openmpi/4.1.6
. spack/share/spack/setup-env.sh
spack env activate platypus
. platypus_env.sh
```

Useful if you login again after installing spack and need to build
MFEM/MOOSE/platypus.

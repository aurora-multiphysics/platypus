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

clone and build mfem:

```bash

```

build moose:

```bash

```

build platypus:

```bash

```

# Set up state again:

```bash
module load llvm/18.1.8
module load cuda/12.6.2
module load openmpi/4.1.6
spack env activate platypus
. spack/share/spack/setup-env.sh
. platypus_env.sh
```

Useful if you login again after installing spack and need to build
MFEM/MOOSE/platypus.

## Platypus on CSD3

First, navigate to the location you wish Platypus to be stored, and clone the repository using

``` {.sh}
git clone https://github.com/aurora-multiphysics/platypus.git
```

Platypus can then be built using the build scripts in `platypus/scripts`.

### GPU-Ready Platypus

For use of Platypus on the Ampere nodes of CSD3, first copy the build script from
`https://github.com/aurora-multiphysics/platypus/blob/main/scripts/build-platypus-csd3-ampere.sh`
and build Platypus using

``` {.sh}
sbatch build-platypus-csd3-ampere.sh
```

!alert! warning title=Existing Spack Installation

By default, this script will uninstall all spack modules currently associated with the architecture
defined in the `ARCH` variable. If this is not desired, comment out the
`spack uninstall -ay arch=${ARCH}` line in the `install_spack_deps()` function before use.

!alert-end!

### CPU-Only Platypus

You may wish to run Platypus using only the Sapphire Rapids nodes of CSD3,
and do not wish to install dependencies associated with use of available GPUs.
In this case, copy the build script from
`https://github.com/aurora-multiphysics/platypus/blob/main/scripts/build-platypus-csd3-sapphire.sh`
and build Platypus using

``` {.sh}
sbatch build-platypus-csd3-sapphire.sh
```

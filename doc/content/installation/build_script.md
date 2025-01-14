## Using the generic build script

Platypus contains a build script which can be used to install all of its dependencies and Platypus itself on any machine, according to user-defined preferences. To make use of it, first make sure you are building on the machine where you intend to run, as cross-compiling is currently not supported. Then, navigate to Platypus' root directory and run

``` {.sh}
./scripts/build-platypus.sh
```
followed by flags defining the specifics of the build. The `-h` or `--help` flag can be used to show all options. They are also listed here:

| Flag | Description |
| ----------- | ----------- |
| `-h`, `--help` | Shows the script usage manual. |
|  `-g`, `--gpu` | Defines a GPU build. If this option is not added, a CPU build is assumed. |
|  `-b=[...]`, `--gpu-backend=[...]` | Defines the GPU backend to be used. Options are `cuda` and `rocm` for NVIDIA and AMD GPUs, respectively. |
|  `-a=[...]`, `--gpu-arch=[...]` | Defines the target GPU architecture. For CUDA backends, use only the number. For instance, to target a GPU whose arch code is sm_80, you would add `-a=80` |
|  `-mpicxx=[<path>]`, `--ompi-cxx=[<path>]` | Path to a C++ compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds. |
|  `-mpicc=[<path>]`, `--ompi-cc=[<path>]` | Path to a C compiler binary in case you wish to wrap the MPI compiler with one that is different to the one it was built with for the MFEM, MOOSE and Platypus builds. |
| `-p=[<name> <version> <path>]`,` --package=[<name> <version> <path>]` | Adds an external package to the spack environment so that it is not built by spack. It is possible to add any number of packages. |
| `-c=[<name> <version> <options>]`,` --compiler=[<name> <version> <options>]` | Adds an external compiler to the spack environment. It is possible to add any number of compilers. In `<options>`, one would include the path to CC, CXX, F77, and FC compilers. It is not necessary to fill them all. See example below. |

As an example, should we wish to build Platypus on a system with AMD GPUs and with some pre-defined compilers and packages already present on the machine, one possible invocation command would be:

``` {.sh}
./scripts/build-platypus -g \
                         -b=rocm \
                         -a=gfx942 \
                         -mpicxx=/opt/rocm/bin/amdclang++ \
                         -mpicc=/opt/rocm-6.2.4/bin/amdclang \
                         -p="hip 6.2.4 /opt/rocm-6.2.4/" \
                         -p="rocrand 6.2.4 /opt/rocm-6.2.4" \
                         -c="clang 16.0.0 CXX=/opt/llvm/bin/clang++ CC=/opt/llvm/bin/clang F77=/opt/llvm/bin/flang"
```


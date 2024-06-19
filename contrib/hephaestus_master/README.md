# Hephaestus
Hephaestus is an MFEM-based library, created to provide access to a range of finite element electromagnetics formulations and scalable solvers for large-scale engineering analysis on complex geometries.

Hephaestus depends on the [MFEM](https://mfem.org/) finite element library, and is used in the [MOOSE](https://github.com/idaholab/moose) application [Apollo](https://github.com/aurora-multiphysics/apollo) to enable coupling between MOOSE-based multiphysics simulations and MFEM-based electromagnetic solves, with flexibility in the choice of electromagnetic formulation.

Hephaestus is still under active development and is being updated frequently.
# Getting started
Docker images of Hephaestus for Ubuntu with all dependencies are built weekly and uploaded to DockerHub, and
can be downloaded via
```
docker pull alexanderianblair/hephaestus:master
```
Once downloaded, the image can be run in interactive mode with the command
```
docker run -it alexanderianblair/hephaestus:master
```
Additional information and options for using Docker can be found at this [tutorial](https://docs.docker.com/get-started/) on the Docker website.

Alternatively, up-to-date images of only the current dependencies for Hephaestus can be downloaded from
```
docker pull alexanderianblair/hephaestus-deps:master
```
for those who with to build Hephaestus themselves.

Dockerfiles used to build these images can be found in the `hephaestus/docker` directory.
# Building software and tests
Hephaestus can be built with the following commands from the top level `hephaestus` directory in either of the above containers:

    mkdir build
    cd build
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/opt/mfem/build ..
    ninja
    ninja test

Running `ninja test` after `hephaestus` is built will run the entire set of tests found in `hephaestus/test`. Running a specific test is possible by providing the test name when running the corresponding test executable; for example:
```
/opt/hephaestus/bin/integration_tests TestAFormSource
```
# Examples
In addition to the set of tests in `hephaestus/test`, longer running examples of Hephaestus being used can be found in `hephaestus/examples`, with associated executables found in `hephaestus/examples/bin`.
```
 mpirun -n 4 -genv OMP_NUM_THREADS=1 /opt/hephaestus/examples/bin/team7
```
# Related Projects
Hephaestus is used extensively by the MOOSE application [Apollo](https://github.com/aurora-multiphysics/apollo). More information can be found on the Apollo GitHub pages.

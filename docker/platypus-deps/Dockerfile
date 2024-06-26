# Build Platypus dependencies
# Get base image
FROM alexanderianblair/hephaestus-deps:master

# By default four cores are used to compile
ARG compile_cores=4

# By default we install everything under /opt
ARG WORKDIR=opt

# Configure MOOSE
WORKDIR /$WORKDIR
# hadolint ignore=DL3003
RUN export MOOSE_JOBS=$compile_cores && \
    git clone https://github.com/idaholab/moose && \
    cd moose && \
    git checkout master && \
    export PETSC_DIR=/$WORKDIR/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    export CC=mpicc && \
    export CXX=mpicxx && \
    export F90=mpif90 && \
    export F77=mpif77 && \
    export FC=mpif90  && \
    METHODS=opt ./scripts/update_and_rebuild_libmesh.sh --with-mpi && \
    ./configure --with-derivative-size=200

# Build WASP
WORKDIR /$WORKDIR/moose/scripts
RUN METHODS=opt ./update_and_rebuild_wasp.sh

# Build MOOSE
WORKDIR /$WORKDIR/moose/framework
RUN METHOD=opt make -j$compile_cores

WORKDIR /$WORKDIR/moose/test
RUN METHOD=opt make -j$compile_cores

# This is needed or it mpiexec complains because docker runs as root
# Discussion on this issue https://github.com/open-mpi/ompi/issues/4451
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# WORKDIR /$WORKDIR/moose/test
# RUN ./run_tests -j"$compile_cores" --re='(?!initial_condition.newton_with_exact_initialization)'

# ENV OMPI_ALLOW_RUN_AS_ROOT=
# ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=

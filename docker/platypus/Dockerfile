# Build and test Platypus image
# Get base image
FROM alexanderianblair/platypus-deps:main

# By default four cores are used to compile
ARG compile_cores=4

# By default checkout main branch from aurora-multiphysics/platypus
ARG build_git_sha=main
ARG build_git_repo=aurora-multiphysics/platypus

# By default we install everything under /opt
ARG WORKDIR=opt

# By default, build without code coverage flags
ARG coverage=false

# Get Platypus with Hephaestus
WORKDIR /$WORKDIR 
# hadolint ignore=DL3003
RUN git clone https://github.com/$build_git_repo && \
    cd platypus && \
    git checkout $build_git_sha && \
    git submodule update --init --recursive

# Build Hephaestus
WORKDIR /$WORKDIR/platypus/contrib/hephaestus/
# hadolint ignore=DL3003
RUN mkdir build && \
    cd build && \
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/$WORKDIR/mfem/build .. && \
    ninja

# Test Hephaestus
WORKDIR /$WORKDIR/platypus/contrib/hephaestus/build
RUN ninja test

# Build Platypus
WORKDIR /$WORKDIR/platypus
RUN make -j$compile_cores coverage=$coverage

# Test Platypus
WORKDIR /$WORKDIR/platypus
RUN make test linkcoverage=$coverage

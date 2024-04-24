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

# Get Platypus with Hephaestus
RUN cd /$WORKDIR && \
    git clone https://github.com/$build_git_repo && \
    cd platypus && \
    git checkout $build_git_sha && \
    git submodule update --init --recursive

# Build Hephaestus
RUN cd /$WORKDIR/platypus/contrib/hephaestus/ && \
    mkdir build && \
    cd build && \
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/$WORKDIR/mfem/build .. && \
    ninja

# Test Hephaestus
RUN cd /$WORKDIR/platypus/contrib/hephaestus/build && \
    ninja test

# Build Platypus
RUN cd /$WORKDIR/platypus && \
    make -j$compile_cores

# Test Platypus
RUN cd /$WORKDIR/platypus && \
    make test

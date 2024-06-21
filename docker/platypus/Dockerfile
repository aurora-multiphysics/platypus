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

# Get Platypus
RUN cd /$WORKDIR && \
    git clone https://github.com/$build_git_repo && \
    cd platypus && \
    git checkout $build_git_sha

# Build Platypus
RUN cd /$WORKDIR/platypus && \
    make -j$compile_cores

# Test Platypus
RUN cd /$WORKDIR/platypus && \
    make test

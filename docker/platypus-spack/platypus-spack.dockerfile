FROM ubuntu:22.04

ARG UID=1000
ARG GID=1000
ARG UNAME=pspack

RUN apt-get update && \
	# install:
	# build-essential and gfortran to build stuff
	# git to clone stuff
	# python for spack
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
		build-essential \
		gfortran \
		git \
		python3 python3-pip python-is-python3 \
		libtirpc-dev \
		unzip && \
	# clean apt-cache to reduce size
	rm -rf /var/lib/apt/lists/*

RUN userdel -r ubuntu

RUN groupadd -g $GID -o $UNAME
RUN useradd -lm -u $UID -g $GID -s /bin/bash $UNAME

USER $UNAME
ENV HOME=/home/$UNAME
ENV USER=$UNAME
WORKDIR $HOME

SHELL ["/bin/bash", "-c"]

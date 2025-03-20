## Platypus Docker Container

Docker images of Platypus for Ubuntu with all dependencies are built weekly
and uploaded to DockerHub, and can be downloaded via

``` {.sh}
docker pull alexanderianblair/platypus:main
```

Once downloaded, the image can be run in interactive mode with the command

``` {.sh}
docker run -it alexanderianblair/platypus:main
```

The Platypus executable can then be found at `/opt/platypus/platypus-opt` inside the container.

Additional information and options for using Docker can be found at
this [tutorial](https://docs.docker.com/get-started/) on the Docker site.

!alert! tip title=Platypus Dependencies
Alternatively, a container containing up-to-date images of only
the current dependencies for Platypus can be downloaded from

``` {.sh}
docker pull alexanderianblair/platypus-deps:main
```

for those who wish to build Platypus themselves.

!alert-end!


### Running Tests

Platypus can then be built with the following commands from the top level `platypus`
directory in either of the above containers.

To build and run the entire set of regression tests for platypus found in `/opt/platypus/test`, run

``` {.sh}
METHOD=opt make test
```

from the `/opt/platypus` directory. Running a specific test (or tests) is possible
by using the `/opt/platypus/run_tests` python script and using the `--re`
command line argument to pass in a regular expression satisfied by the names
of the tests you wish to run; for example:

``` {.sh}
METHOD=opt /opt/platypus/run_tests --re=MFEMDiffusion
```

Similarly, unit tests can be built and run in a similar way, by running `make test` from `/opt/platypus/unit`.
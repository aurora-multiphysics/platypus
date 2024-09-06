Test the spack install script works

Build the docker image
```bash
docker build --build-arg "UID=$(id -u)" --build-arg "GID=$(id -g)" -t platypus-spack -f platypus-spack.dockerfile --no-cache .
```

Start a new interactive container:
```bash
docker run -v $PWD:/opt/platypus -it ubuntu:24.04
```

Run test install test script:
```bash
./spack-install.sh
```

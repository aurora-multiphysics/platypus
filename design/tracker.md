# MOOSE/libMesh FE Limitations Issue Tracker

## Lack of separation between mesh elements and finite element space

### Issue Description

- Mesh element order determines maximum polynomial order simulation can be run with;
- Assumption that higher order finite elements contain a superset of the low order FE nodes when computing the necessary constraints at non-conforming interfaces, e.g. AMR interfaces;
- Limited number of kernels using `LAGRANGE_VEC`; most MOOSE modules use three scalar `LAGRANGE` variables added to the problem, affecting problem ordering.

### Affected Problems

- All non-trivial problems since, in general, problems benefit quicker from p-refinement than from h-refinement, i.e. for the same error it is usually more efficient to use a higher order set of basis functions than to refine the underlying mesh.

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

## Limited support for H(curl) and H(div) FE types

### Issue Description

- Currently support for only lowest order ND1 and RT FEs;
- These must be run on second order `TET14`, `HEX27`, `TRI6` and `QUAD9` elements, with the exception of first order ND1 FEs, which can be run on `TET10` and `HEX20` meshes, and that, in 2d, the user can also use `TRI7` and `QUAD8` elements since the implementation does not use the node at the centroid.

### Affected Problems

- Any non-trivial problem in H(curl) or H(div) discretised with ND1 or RT FEs, respectively, since these are likely to benefit from p-refinement for the reasons above.

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

## Orientation Information for hypre preconditioners

### Issue Description

- Orientation information needed when setting up hypre's AMS and ADS preconditioners not provided, since orientation information is currently calculated on-the-fly in libMesh. These preconditioners are thus currently unavailable to the user in both libMesh and MOOSE.

### Affected Problems

- Any non-trivial problem in H(curl) or H(div) discretised with ND1 or RT FEs, respectively.

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

## Dirichlet BCs for ND1 and RT FE types

### Issue Description

- ND1 and RT FEs can currently only use integrated BCs. Dirichlet BCs are approximated using penalty methods, which can be unstable as they can lead to poor conditioning;
- The suggested solution of overriding the local element matrix and residual vector currently does not work with the AD kernels in MOOSE as these directly write to global data structures, completely bypassing the local data structures currently being intercepted.

### Affected Problems

- Any kind of problem which uses ND1 or RT FEs, such as electromagnetic formulations.

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

## Poor support for systems of complex variables

### Issue Description

- Systems of complex variables are difficult to set up in MOOSE, and separate kernels for real and imaginary components must be created - highly error prone.

### Affected Problems

- Any kind of problem we might like to solve in the frequency domain.

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

## No support for General Field Vector Transfers

### Issue Description

- Vector variables can only be transferred using `CopyTransfer`s in MOOSE; particularly relevant for transfers of volumetric force densities (`MONOMIAL_VEC`) and nodal forces (`LAGRANGE_VEC`);
- A possible workaround would be to decompose vector field variables in each of their components and then to use MOOSE's more robust/flexible support for transferring the resulting scalar field variables instead.

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues/PRs

- Issue #
- Issue #
- Issue #

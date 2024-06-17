# MOOSE/libMesh FE Limitations Issue Tracker

## Lack of separation between mesh elements and finite element space

### Issue Description

- Mesh element order determines maximum polynomial order simulation can be run with
- Assumption that higher order finite elements contain a superset of the low order FE nodes
- Limited number of kernels using LAGRANGE_VEC; most MOOSE modules use three scalar LAGRANGE variables added to the problem, affecting problem ordering.

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

## Limited support for H(curl) and H(Div) FE types

### Issue Description

- Currently support for only lowest order ND1 and RT FEs, which must be run on second order TET/HEX/TRI/QUAD elements

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

## Orientation Information for Hypre preconditioners

### Issue Description

- Orientation information needed when setting up HypreAMS and HypreADS preconditioners not provided, since orientation information is currently calculated on-the-fly in libMesh. These preconditioners are poor as a result, as currently implemented

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

## Dirichlet BCs for ND and RT FE types

### Issue Description

- ND and RT FEs can currently only use integrated BCs. Dirichlet BCs are approximated using penalty methods, which can be unstable as they can lead to poor conditioning.
- Seems to be related to AD implementation in libMesh

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

## Poor support for systems of complex variables

### Issue Description

- Systems of complex variables are difficult to set up in MOOSE, and separate kernels for real and imaginary components must be created - highly error prone

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

## No support for General Field Vector Transfers

### Issue Description

- Vector variables can only be transferred using CopyTransfers in MOOSE; particularly relevant for transfers of volumetric force densities (MONOMIAL_VEC) and nodal forces (LAGRANGE_VEC)

### Affected Problems

### Relevant MOOSE/libMesh GitHub issues

- Issue #
- Issue #
- Issue #

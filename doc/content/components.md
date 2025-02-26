# Key Components

## Description

Platypus components are intended to be used together in MOOSE simulations with the
[`MFEMProblem`](source/problem/MFEMProblem.md) `Problem` type. In general, regular MOOSE objects and
types should not be used in a sub-app using [`MFEMProblem`](source/problem/MFEMProblem.md), with the
notable exception of MOOSE `Functions`, which may be used and included in Platypus problems.
Documentation on the main components comprising Platypus can be found in the following files:

## Problem Assembly

- [PlatypusApp](source/base/PlatypusApp.md): Defines relative dependencies between Actions when
  parsing user input files during problem set-up, and thus the order in which they are executed.

- [MFEMExecutioner](source/executioners/MFEMExecutioner.md): Controls the assembly and execution of
  an MFEM FE problem. Choice of device (CPU/GPU) and assembly level is controlled here.

- [MFEMProblem](source/problem/MFEMProblem.md): Builder responsible for constructing and adding user
  requested objects to the FE problem.

- [MFEMProblemData](source/problem/MFEMProblemData.md): Struct containing data associated with the
  state of the MFEMProblem.

- [ProblemOperator](source/problem_operators/problem_operator.md): Responsible for updating the
  state of the system. For transient problems, the `ProblemOperator` is passed to the timestepper.

## Problem Data

- [MFEMMesh](source/mesh/MFEMMesh.md): Builds an `mfem::ParMesh` object from the provided mesh input
  file.

- [MFEMFESpace](source/fespaces/MFEMFESpace.md): Base class for defining the finite element
  space with respect to which an  `MFEMVariable` is defined.

- [MFEMVariable](source/variables/MFEMVariable.md): Wrapper for set-up of `mfem::ParGridFunction`
  objects representing trial variables; stores a set of (true) degrees of freedom for a variable
  with respect to an `MFEMFESpace`.

- [MFEMKernel](source/kernels/MFEMKernel.md): Base class defining integrators contributing to the
  weak form built inside `EquationSystem`. Contains methods returning
  `mfem::NonlinearFormIntegrators` and `mfem::LinearFormIntegrators` to add to the weak form, along
  with marker arrays labelling the volumes the (domain) integrators are to be applied to.

- [MFEMBoundaryCondition](source/bcs/MFEMBoundaryCondition.md): Base class defining essential
  (Dirichlet) and integrated boundary conditions to apply to the equation system.

- [EquationSystem](source/equation_systems/equation_system.md): Class responsible for defining and
  assembling the weak form into an `mfem::Operator`, that can be used to solve for the set of trial
  variables at a given time or timestep, which can be passed to an `mfem::NewtonSolver`.

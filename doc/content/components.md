# Key Components

## Description

Documentation on the main components comprising Platypus can be found in:

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

- [MFEMFESpace](source/fespaces/MFEMFESpace.md): Class responsible for defining the finite element
  space that an `MFEMVariable` is defined with respect to.

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
  assembling the weak form into an `mfem::Operator` that solves an iteration of the FE problem and
  can be passed to an `mfem::NewtonSolver`.
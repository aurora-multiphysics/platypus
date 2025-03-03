# Examples

## Description

As part of the Platypus test suite, we solve a number of sample problems that may serve as a useful
starting point for users to adapt:

## Thermal Problems

- [Steady State Diffusion](examples/Diffusion.md): Steady state diffusion problem for the
  temperature profile across a mug with fixed temperatures on two boundaries.

- [Transient Heat Transfer](examples/HeatTransfer.md): Transient heat conduction problem with a
  boundary parameterized by a heat transfer coefficient that exchanges heat with a thermal
  reservoir.

## Mechanical Problems

- [Linear Elasticity](examples/LinearElasticity.md): Solves a 3D linear elasticity problem for
  the deformation of a multi-material cantiliver beam. This example
  is based on [MFEM Example 2](https://mfem.org/examples/).

## Electromagnetic Problems

- [Definite Maxwell](examples/DefiniteMaxwell.md): Solves a 3D electromagnetic diffusion problem for
  the electric field on a cube missing an octant, discretized using $H(\mathrm{curl})$ conforming
  Nédélec elements. This example is based on [MFEM Example 3](https://mfem.org/examples/).

- [GradDiv](examples/GradDiv.md): Solves a diffusion problem for a vector field
  on a cuboid domain, discretized using $H(\mathrm{div})$ conforming
  Raviart-Thomas elements. This example is based on [MFEM Example 4](https://mfem.org/examples/)
  and is relevant for solving Maxwell's equations using potentials without the Coulomb gauge.

## Fluids

- [Darcy Problem](examples/Darcy.md): Solves a simple 2D/3D mixed Darcy problem corresponding to the saddle point system.
We discretize with Raviart-Thomas finite elements (velocity \( \mathbf{u} \)) and piecewise discontinuous
polynomials (pressure \( p \)).
This example demonstrates the use of transpose in the input file used in mixed problems with different trial and
test variables when setting up block matrices.
This example is based on [MFEM Example 5](https://mfem.org/examples/).

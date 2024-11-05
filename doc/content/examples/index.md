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

## Electromagnetic Problems

- [DefiniteMaxwell](examples/DefiniteMaxwell.md): Solves a 3D electromagnetic diffusion problem for
  the electric field on a cube missing an octant, discretized using $H(\mathrm{curl})$ conforming
  Nédélec elements. This example is based on [MFEM Example 3](https://mfem.org/examples/).

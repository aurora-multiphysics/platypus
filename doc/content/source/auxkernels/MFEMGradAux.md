# MFEMGradAux

## Summary

!syntax description /AuxKernels/MFEMGradAux

## Overview

AuxKernel for calculating the gradient of a scalar $H^1$ conforming source variable and storing it in
a scalar elemental result variable defined on a $H(\mathrm{curl})$ conforming ND FE space.

The result may be scaled by an optional (global) scalar factor.

!equation
v = \lambda \vec\nabla u

where $u \in H^1$, $\vec v \in H(\mathrm{curl})$ and $\lambda$ is a scalar constant.

## Example Input File Syntax

!listing kernels/diffusion.i

!listing kernels/irrotational.i

!syntax parameters /AuxKernels/MFEMGradAux

!syntax inputs /AuxKernels/MFEMGradAux

!syntax children /AuxKernels/MFEMGradAux

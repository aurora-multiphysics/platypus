# MFEMDiffusionKernel

## Summary

!syntax description /Kernels/MFEMDiffusionKernel

## Overview

Adds the domain integrator for integrating the bilinear form

!equation
(k\vec\nabla u, \vec\nabla v)_\Omega \,\,\, \forall v \in V

where $u, v \in H^1$.

This term arises from the weak form of the Laplacian operator

!equation
-k \nabla^2 u

## Example Input File Syntax

!listing kernels/diffusion.i

!syntax parameters /Kernels/MFEMDiffusionKernel

!syntax inputs /Kernels/MFEMDiffusionKernel

!syntax children /Kernels/MFEMDiffusionKernel
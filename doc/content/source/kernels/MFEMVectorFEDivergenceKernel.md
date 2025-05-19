# MFEMVectorFEDivergenceKernel

## Summary

!syntax description /Kernels/MFEMVectorFEDivergenceKernel

## Overview

Adds the domain integrator for integrating the mixed bilinear form

!equation
(\sigma \nabla V, u')_\Omega

where $\sigma$ is a coefficient, $V$ is a test function, and $u'$ is a trial function.

This term arises from the weak form of the divergence operator.

## Example Input File Syntax

!syntax parameters /Kernels/MFEMVectorFEDivergenceKernel

!syntax inputs /Kernels/MFEMVectorFEDivergenceKernel

!syntax children /Kernels/MFEMVectorFEDivergenceKernel

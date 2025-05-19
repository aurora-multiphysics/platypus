# MFEMDomainLFKernel

## Summary

!syntax description /Kernels/MFEMDomainLFKernel

## Overview

Adds the domain integrator for integrating the linear form

!equation
(f, v)_\Omega \,\,\, \forall v \in V

where $v \in H^1$ or $L^2$ is the test variable and $f$ is a
forcing coefficient.

This term arises from the weak form of the forcing term

!equation
f

## Example Input File Syntax

!listing kernels/darcy.i block=Kernels

!syntax parameters /Kernels/MFEMDomainLFKernel

!syntax inputs /Kernels/MFEMDomainLFKernel

!syntax children /Kernels/MFEMDomainLFKernel

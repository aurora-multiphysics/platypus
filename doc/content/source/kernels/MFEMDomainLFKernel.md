# MFEMDomainLFKernel

## Summary

!syntax description /Kernels/MFEMDomainLFKernel

## Overview

Adds the domain integrator for the **linear form**:

!equation
( \mathbf{f}, \mathbf{u'} )_{\Omega}

where:
- \( \mathbf{f} \) is a given **vector-valued function**
- \( \mathbf{u'} \) is a **test function** from the vector finite element space.

where $f$ is a given function, and $u'$ is a test function. The term arises from the weak form for integrating a given function against the test function.

## Example Input File Syntax

!syntax parameters /Kernels/MFEMDomainLFKernel

!syntax inputs /Kernels/MFEMDomainLFKernel

!syntax children /Kernels/MFEMDomainLFKernel

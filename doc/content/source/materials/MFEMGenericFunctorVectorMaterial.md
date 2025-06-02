# MFEMGenericFunctorVectorMaterial

## Summary

!syntax description /Materials/MFEMGenericFunctorVectorMaterial

## Overview

`MFEMGenericFunctorVectorMaterial` defines one or more vector material properties with values
obtained from coefficients on one or more subdomains of the mesh, given by the `blocks` parameter
if provided, or applied to the entire mesh if missing. The vector material properties are named
according to members in the `prop_names` parameter, with respective coefficients used to get property
values given by the members of `prop_values`. The coefficients in `prop_names` must be vector-valued.

!syntax parameters /Materials/MFEMGenericFunctorVectorMaterial

!syntax inputs /Materials/MFEMGenericFunctorVectorMaterial

!syntax children /Materials/MFEMGenericFunctorVectorMaterial

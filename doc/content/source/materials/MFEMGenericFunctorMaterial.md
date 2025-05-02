# MFEMGenericFunctorMaterial

## Summary

!syntax description /Materials/MFEMGenericFunctorMaterial

## Overview

`MFEMGenericFunctorMaterial` defines one or more scalar material properties with values obtained from a coefficient on
one or more subdomains of the mesh, given by the `blocks` parameter if provided, or applied to the
entire mesh if missing. The scalar material properties are named according to members in the
`prop_names` parameter, with respective coefficients used to get property values given by the members of `prop_values`.

!syntax parameters /Materials/MFEMGenericFunctorMaterial

!syntax inputs /Materials/MFEMGenericFunctorMaterial

!syntax children /Materials/MFEMGenericFunctorMaterial

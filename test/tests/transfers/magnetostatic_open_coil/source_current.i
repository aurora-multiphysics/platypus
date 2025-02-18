[Mesh]
  type = MFEMMesh
  file = coil.gen
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HCurlFESpace]
    type = MFEMFESpace
    fec_type = ND
    fec_order = FIRST
  []
[]

[Variables]
  [electric_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []  
[]

[AuxVariables]
  [neg_grad_electric_potential]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
[]

[AuxKernels]
  [Gradient]
    type = MFEMGradAux
    source = electric_potential
    variable = neg_grad_electric_potential
    scale_factor = -1.0
  []  
[]

[BCs]
  [high_terminal]
    type = MFEMScalarDirichletBC
    variable = electric_potential
    boundary = 1
    value = 1.0
  []
  [low_terminal]
    type = MFEMScalarDirichletBC
    variable = electric_potential
    boundary = 2
    value = -1.0
  []
[]

[Materials]
  [Coil]
    type = MFEMGenericConstantMaterial
    prop_names = electric_conductivity
    prop_values = 62.83185
    block = 1
  []
  [Vacuum]
    type = MFEMGenericConstantMaterial
    prop_names = electric_conductivity
    prop_values = 62.83185e-6
    block = 2
  []
  [Core]
    type = MFEMGenericConstantMaterial
    prop_names = electric_conductivity
    prop_values = 62.83185e-6
    block = 3
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = electric_potential
    coefficient = electric_conductivity
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = boomeramg
  l_tol = 1e-16
  l_max_its = 1000  
[]


[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/DiffusionSub
    vtk_format = ASCII
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

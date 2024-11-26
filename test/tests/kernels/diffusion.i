[Mesh]
  type = MFEMMesh
  file = gold/mug.e
  dim = 3
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
[]

[Variables]
  [diffused]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[ICs]
  [diffused_ic]
    type = MFEMScalarIC
    coefficient = one
    variable = diffused
  []
[]

[Functions]
  [one]
    type = ParsedFunction
    expression = 1.0
  []
  [zero]
    type = ParsedFunction
    expression = 0.0
  []
[]

[BCs]
  [bottom]
    type = MFEMScalarDirichletBC
    variable = diffused
    boundary = '1'
    coefficient = one
  []
  [low_terminal]
    type = MFEMScalarDirichletBC
    variable = diffused
    boundary = '2'
    coefficient = zero
  []
[]

[Materials]
  [Substance]
    type = MFEMGenericConstantMaterial
    prop_names = diffusivity
    prop_values = 1.0
  []
[]

[Coefficients]
  [zero]
    type = MFEMFunctionCoefficient
    function = zero
  []
  [one]
    type = MFEMFunctionCoefficient
    function = one
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = diffused
    coefficient = diffusivity
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

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Diffusion
    vtk_format = ASCII
  []
[]

[Mesh]
  type = MFEMMesh
  file = gold/mug.e
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

[Materials]
  [DiffusiveMaterial]
    type = MFEMGenericConstantMaterial
    prop_names = diffusivity
    prop_values = 1.0
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = diffused
    coefficient = diffusivity
  []
[]

[Coefficients]
  [TopValue]
    type = MFEMConstantCoefficient
    value = 0.0
  []
  [BottomValue]
    type = MFEMConstantCoefficient
    value = 1.0
  []
[]

[BCs]
  [MugBottom]
    type = MFEMScalarDirichletBC
    variable = diffused
    boundary = 1
    coefficient = BottomValue
  []
  [MugTop]
    type = MFEMScalarDirichletBC
    variable = diffused
    boundary = 2
    coefficient = TopValue
  []
[]

[Preconditioner]
  [BoomerAMG]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = BoomerAMG
  l_tol = 1e-16
  l_max_its = 1000  
[]

[Executioner]
  type = MFEMSteady
  device = cpu
  assembly_level = legacy
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/MFEM/Diffusion
    vtk_format = ASCII
    execute_on = 'final'
  []
[]

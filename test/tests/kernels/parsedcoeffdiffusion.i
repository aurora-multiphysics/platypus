[Mesh]
  type = MFEMMesh
  file = gold/star.mesh
  dim = 2
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []  
[]

[Variables]
  [concentration]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [concentration_gradient]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  [grad]
    type = MFEMGradAux
    variable = concentration_gradient
    source = concentration
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  [exterior_boundary]
    type = MFEMScalarDirichletBC
    variable = concentration
    boundary = '1'
    value = 1.0
  []
[]

[Materials]
  [Substance_a]
    type = MFEMParsedCoeffMaterial
    expression = 'concentration'
    var_names = 'concentration'
    prop_name = source_coeff
  []
  [Substance_b]
    type = MFEMGenericConstantMaterial
    prop_names = diffusivity
    prop_values = 1.0
  []
[]

[Kernels]
  [mass]
   type = MFEMMassKernel
   variable = concentration
   coefficient = diffusivity
  []
  [diff]
    type = MFEMDiffusionKernel
    variable = concentration
    coefficient = diffusivity
  []
  [source]
    type = MFEMScalarDomainLFKernel
    variable = concentration
    coefficient = source_coeff
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
    file_base = OutputData/ParsedCoeffNLDiffusion
    vtk_format = ASCII
  []
[]

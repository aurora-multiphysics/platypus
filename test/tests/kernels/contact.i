[Mesh]
  type = MFEMMesh
  file = gold/two-hex.mesh
  dim = 3
  parallel_refine = 2
  displacement = 'displacement'
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Kernels]
  [diff]
    type = MFEMLinearElasticityKernel
    variable = displacement
    lambda = lambda
    mu = mu
  []
[]

[Materials]
  [Rigidium]
    type = MFEMGenericConstantMaterial
    prop_names = 'lambda mu'
    prop_values = '50.0 50.0'
    block = '1 2'
  []
[]

[BCs]
  [interface]
    type = MFEMContactBC
    variable = displacement
  []
[]

[Executioner]
  type = MFEMContact
  device = cpu
[]

[Solver]
  type = MFEMMinResSolver
  preconditioner = hyprediag
  rel_tol   = 1e-12
  l_max_its = 1500  
[]

[Preconditioner]
  [hyprediag]
    type = MFEMHypreDiagScale
  []
[]



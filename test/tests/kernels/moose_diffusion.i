[Mesh]
  type = FileMesh
  file = gold/mug.e
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [diffused]
    order = FIRST
    family = LAGRANGE    
  []
[]

[Materials]
  [DiffusiveMaterial]
    type = GenericConstantMaterial
    prop_names = diffusivity
    prop_values = 1.0
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = diffused
    diffusivity = diffusivity
  []
[]

[BCs]
  [MugBottom]
    type = DirichletBC
    variable = diffused
    boundary = 1
    value = 1.0
  []
  [MugTop]
    type = DirichletBC
    variable = diffused
    boundary = 2
    value = 0.0
  []
[]

[Executioner]
  type = Steady
  solve_type = 'LINEAR'
  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type'
  petsc_options_value = 'gmres hypre boomeramg'
[]

[Outputs]
  [ExodusOutput]
    type = Exodus
    file_base = OutputData/MOOSE/Diffusion
    execute_on = 'final'
  []
[]

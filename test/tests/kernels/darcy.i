[Mesh]
  type = MFEMMesh
  file = gold/star.mesh
  uniform_refine = 4
[]

[Problem]
  type = MFEMProblem
[]

[Functions]
  [exact_velocity]
    type = ParsedVectorFunction
    expression_x = '-exp(x) * sin(y)'
    expression_y = '-exp(x) * cos(y)'
  []
  [exact_pressure]
    type = ParsedFunction
      expression = 'exp(x) * sin(y)'
  []
  [f_bdr]
    type = ParsedFunction
    expression = '-exp(x) * sin(y)'
  []    
  [g_force]
    type = ParsedFunction
    expression = '0.0'
  []
[]

[FESpaces]
  [L2FESpace]
    type = MFEMFESpace
    fec_type = L2
    fec_order = SECOND
  []
  [HDivFESpace]
    type = MFEMFESpace
    fec_type = RT
    fec_order = SECOND
  []  
[]

[Variables]
  [velocity]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [pressure]
    type = MFEMVariable
    fespace = L2FESpace
  []  
[]

[BCs]
  [flux_boundaries]
    type =   MFEMVectorFEFunctionBoundaryFluxIntegratedBC
    variable = velocity
    function = f_bdr
    boundary = 1
  []
[]

[Kernels]
  [VelocityMass]
    type = MFEMVectorFEMassKernel
    variable = velocity
    coefficient = one
  []
  [PressureGrad]
    type = MFEMVectorFEDivergenceKernel
    trial_variable = pressure
    variable = velocity
    coefficient = neg_one
    transpose = true
  []
  [VelocityDiv]
    type = MFEMVectorFEDivergenceKernel
    trial_variable = velocity
    coefficient = neg_one
    variable = pressure
  []
  [PressureForcing]
    type = MFEMDomainLFKernel
    coefficient = g_force
    variable = pressure
  []    
[]

[Materials]
  [Substance]
    type = MFEMGenericConstantMaterial
    prop_names = 'one neg_one'
    prop_values = '1.0 -1.0'
  []
  [Forcing]
    type = MFEMGenericFunctionMaterial
    prop_names = 'g_force'
    prop_values = 'g_force'
  []  
[]

[Solver]
  type = MFEMSuperLU
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [pressure_error]
    type = MFEML2Error
    variable = pressure
    function = exact_pressure
    execution_order_group = 1
  []
  [velocity_error]
    type = MFEMVectorL2Error
    variable = velocity
    function = exact_velocity
    execution_order_group = 1
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Darcy
    vtk_format = ASCII
  []
  [DarcyErrorCSV]
    type = CSV
    file_base = OutputData/Darcy
  []  
[]

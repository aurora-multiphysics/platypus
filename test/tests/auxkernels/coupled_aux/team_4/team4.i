[Mesh]
  type = MFEMMesh
  file = team4_symmetrized.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMFESpace
    fec_type = RT
    fec_order = CONSTANT
  []  
[]

[Variables]
  [magnetic_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
[]

[AuxVariables]
  [grad_magnetic_potential]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [current_density]
    type = MFEMVariable
    fespace = HDivFESpace
  []  
[]

[AuxKernels]
  [current_density]
    type = MFEMCurlAux
    variable = current_density
    source = magnetic_field
    execute_on = TIMESTEP_END
  []
[]

[Functions]
  # Here, externally applied B field = grad psi
  [source_field_magnitude]
    type = ParsedFunction
    expression = -(B0/tau)*exp(-t/tau)
    symbol_names = 'B0 tau'
    symbol_values = '0.1 0.0119'
  []
[]

[Kernels]
  [curl_rho_curl_H]
    type = MFEMCurlCurlKernel
    variable = magnetic_field
    coefficient = electric_resistivity
  []
  [mu_dH_dt]
    type = MFEMVectorFETimeDerivativeMassKernel
    variable = magnetic_field
    coefficient = magnetic_permeability
  []
  [source_field]
    type = MFEMVectorFEDomainLFCoupledAuxKernel
    trial_variable = grad_magnetic_potential
    variable = magnetic_field
    coefficient = source_field_magnitude
  []
[]

[BCs]
  [tangential_H_bdr]
    type = MFEMVectorTangentialDirichletBC
    variable = magnetic_field
    boundary = '1 2 5 6'
    values = '0.0 0.0 0.0'
  []
[]

[Materials]
  [Brick]
    type = MFEMGenericConstantMaterial
    prop_names = 'electric_resistivity magnetic_permeability'
    prop_values = '3.940e-8 1.25663706e-6'
    block = 1
  []
  [Vacuum]
    type = MFEMGenericConstantMaterial
    prop_names = 'electric_resistivity magnetic_permeability'
    prop_values = '1.0 1.25663706e-6'
    block = 2
  []
  [SourceFieldMagnitude]
    type = MFEMGenericFunctionMaterial
    prop_names = 'source_field_magnitude'
    prop_values = 'source_field_magnitude'
  []  
[]

[MultiApps]
  [./subapp]
    type = FullSolveMultiApp
    input_files = div_free_source.i
    execute_on = INITIAL
  [../]
[]

[Preconditioner]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = ams
  l_tol = 1e-16
  l_max_its = 300  
[]

[Executioner]
  type = MFEMTransient
  dt = 0.001
  start_time = 0.0
  end_time = 0.02
  device = cpu
[]

[Transfers]
  [./from_sub]
      type = MultiAppMFEMCopyTransfer
      source_variable = grad_magnetic_potential
      variable = grad_magnetic_potential
      from_multi_app = subapp
      execute_on = INITIAL
  [../]    
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TEAM4
    vtk_format = ASCII
  []
[]


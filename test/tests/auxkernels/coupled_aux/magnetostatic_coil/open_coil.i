[Mesh]
  type = MFEMMesh
  file = coil.gen
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
  [magnetic_vector_potential]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
[]

[AuxVariables]
  [neg_grad_electric_potential]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [magnetic_flux_density]
    type = MFEMVariable
    fespace = HDivFESpace
  []  
[]

[AuxKernels]
  [curlA]
    type = MFEMCurlAux
    variable = magnetic_flux_density
    source = magnetic_vector_potential
    execute_on = TIMESTEP_END
  []
[]

[Kernels]
  [curlcurlA]
    type = MFEMCurlCurlKernel
    variable = magnetic_vector_potential
    coefficient = permeability
  []
  [source_current]
    type = MFEMVectorFEDomainLFCoupledAuxKernel
    trial_variable = neg_grad_electric_potential
    variable = magnetic_vector_potential
    coefficient = electric_conductivity
  []
[]

[BCs]
  [tangential_A_bdr]
    type = MFEMVectorTangentialDirichletBC
    variable = magnetic_vector_potential
    boundary = '1 2 4'
    values = '0.0 0.0 0.0'
  []
[]

[Materials]
  [Permability]
    type = MFEMGenericConstantMaterial
    prop_names = permeability
    prop_values = 1.0
  []
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

[MultiApps]
  [./subapp]
    type = FullSolveMultiApp
    input_files = source_current.i
    execute_on = INITIAL
  [../]
[]

[Preconditioner]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    singular = true
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = ams
  l_tol = 1e-6
  l_max_its = 100  
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Transfers]
  [./from_sub]
      type = MultiAppMFEMCopyTransfer
      source_variable = neg_grad_electric_potential
      variable = neg_grad_electric_potential
      from_multi_app = subapp
      execute_on = INITIAL
  [../]
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/Diffusion
    vtk_format = ASCII
  []
[]


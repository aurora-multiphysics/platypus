[Mesh]
  type = MFEMMesh
  file = team4_symmetrized.e
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
  [magnetic_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []  
[]

[AuxVariables]
  [grad_magnetic_potential]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
[]

[AuxKernels]
  [Gradient]
    type = MFEMGradAux
    source = magnetic_potential
    variable = grad_magnetic_potential
  []  
[]

[Functions]
  # Here, externally applied B field = grad psi
  [psi_bdr_func]
    type = ParsedFunction
    expression = z
  []
[]

[BCs]
  [normal_B_bc]
    type = MFEMScalarFunctionDirichletBC
    variable = magnetic_potential
    boundary = '1 2'
    function  = psi_bdr_func
  []
[]

[Materials]
  [Brick]
    type = MFEMGenericConstantMaterial
    prop_names = 'electric_resistivity rel_magnetic_permeability'
    prop_values = '3.940e-8 1.0'
    block = 1
  []
  [Vacuum]
    type = MFEMGenericConstantMaterial
    prop_names = 'electric_resistivity rel_magnetic_permeability'
    prop_values = '1.0 1.0'
    block = 2
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = magnetic_potential
    coefficient = rel_magnetic_permeability
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
  l_max_its = 100  
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

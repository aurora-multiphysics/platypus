# Darcy problem based on MFEM Example 5.

[Mesh]
  type = MFEMMesh
  file = gold/beam-tet.mesh
  dim = 3
  uniform_refine = 1
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [VelocityFESpace]
    type = MFEMFESpace
    fec_type = RT
    fec_order = CONSTANT
    vdim = 1
    ordering = "vdim"
  [PressureFESpace]
    type = MFEMFESpace
    fec_type = L2
    fec_order = 1
    vdim = 1
    ordering = "vdim"
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = VelocityFESpace
  [p]
    type = MFEMVariable
    fespace = PressureFESpace
[]

[Functions]
  [uFun_ex]
    type = ParsedVectorFunction
    expression_x = '-exp(x) * sin(y) * cos(z)'
    expression_y = '-exp(x) * cos(y) * cos(z)'
    expression_z = 'exp(x) * sin(y) * sin(z)'
    symbol_names = x, y, z
    symbol_values = 0.0, 0.0, 0.0
  []
  [pFun_ex]
    type = ParsedScalarFunction
    expression = 'exp(x) * sin(y) * cos(z)'
    symbol_names = x, y, z
    symbol_values = 0.0, 0.0, 0.0
  []
  [fFun]
    type = ParsedVectorFunction
    expression_x = '0.0'
    expression_y = '0.0'
    expression_z = '0.0'
    symbol_names = x, y, z
    symbol_values = 0.0, 0.0, 0.0
  []
  [gFun]
    type = ParsedScalarFunction
    expression = '-exp(x) * sin(y) * cos(z)'

    symbol_names = x, y, z
    symbol_values = 0.0, 0.0, 0.0
  []
  [f_natural]
    type = ParsedScalarFunction
    expression = '-exp(x) * sin(y) * cos(z)'

    symbol_names = x, y, z
    symbol_values = 0.0, 0.0, 0.0
  []
[]

[BCs]
  [dirichlet]
    type = MFEMScalarDirichletBC
    variable = p
    boundary = '1 2 3'
    function = f_natural
  []
[]

// [Materials]
//   [Beamium]
//     type = MFEMGenericConstantMaterial
//     prop_names = 'alpha beta'
//     prop_values = '1.0 1.0'
//     block = '1 2'
//   []
// []

[Kernels]
  [mass]
    type = MFEMVectorFEMassKernel
    variable = u
    coefficient = kappa
  []
  [mixedScalarGrad]
    type = MFEMMixedScalarGradientKernel
    variable = u
    trial_variable = p
    coefficient = -1.0
  []
  [vectorDiv]
    type = MFEMVectorFEDivergenceKernel
    variable = p
    trial_variable = u
    coefficient = -1.0
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = F
    function = f
  []
[]

// [Preconditioner]
//   [ADS]
//     type = MFEMHypreADS
//     fespace = HDivFESpace
//   []
// []

// [Solver]
//   type = MFEMCGSolver
//   preconditioner = ADS
//   l_tol = 1e-16
//   l_max_its = 1000
//   print_level = 2
// []

// [Executioner]
//   type = MFEMSteady
//   device = "cpu"
// []

// [Outputs]
//   [ParaViewDataCollection]
//     type = MFEMParaViewDataCollection
//     file_base = OutputData/GradDiv
//     vtk_format = ASCII
//   []
// []

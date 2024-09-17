[Mesh]
  type = MFEMMesh
  file = gold/beam-tet.mesh
  dim = 3
  uniform_refine = 2
[]

[Problem]
  type = MFEMProblem
  device = "cpu"
[]

[FESpaces]
  [H1FESpace]
    type = MFEMFESpace
    fec_type = H1
    fec_order = FIRST
    vdim = 3
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [dirichlet]
    type = MFEMVectorDirichletBC
    variable = displacement
    boundary = '1'
    vector_coefficient = FixedValue
  []
  [pull_down]
    type = MFEMVectorBoundaryIntegratedBC
    variable = displacement
    boundary = '2'
    vector_coefficient = PullDownValue
  []
[]

[Coefficients]
  [lambda]
    type = MFEMConstantCoefficient
    value = 1.0
  []
  [mu]
    type = MFEMConstantCoefficient
    value = 1.0
  []
[]

[VectorCoefficients]
  [FixedValue]
    type = MFEMVectorConstantCoefficient
    value_x = 0.0
    value_y = 0.0
    value_z = 0.0
  []
  [PullDownValue]
    type = MFEMVectorConstantCoefficient
    value_x = 0.0
    value_y = -0.01
    value_z = 0.0
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

[Solver]
  type = MFEMHyprePCG
  l_tol = 1e-6
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/LinearElasticity
    vtk_format = ASCII
  []
[]

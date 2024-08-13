#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"

namespace platypus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms)
*/
struct EquationSystemData
{

  friend class EquationSystemAssembler;
  friend class DiagonalEquationSystemAssembler;
  friend class DenseEquationSystemAssembler;

private:

  using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
  using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
  using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
  using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

  // Test variables are associated with LinearForms,
  // whereas trial variables are associated with gridfunctions.

  // Names of all variables corresponding to gridfunctions. This may differ
  // from test_var_names when time derivatives are present.
  std::vector<std::string> _trial_var_names;
  // Names of all test variables corresponding to linear forms in this equation
  // system
  std::vector<std::string> _test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> _test_pfespaces;

  // Components of weak form. // Named according to test variable
  platypus::NamedFieldsMap<mfem::ParBilinearForm> _blfs;
  platypus::NamedFieldsMap<mfem::ParLinearForm> _lfs;
  //platypus::NamedFieldsMap<mfem::ParNonlinearForm> _nlfs;
  platypus::NamedFieldsMap<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _mblfs; // named according to trial variable

  // Assembly level for the equation system. May be LEGACY, FULL, ELEMENT, or PARTIAL
  mfem::AssemblyLevel _assembly_level;

  std::vector<mfem::Array<int>> _ess_tdof_lists;


  // gridfunctions for setting Dirichlet BCs
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _bc_gridfunc;

  // Array for mixed systems. Can only be non-block-diagonal when @_assembly_level is LEGACY
  mfem::Array2D<const mfem::HypreParMatrix *> _h_blocks;

  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>> _blf_kernels_map;

  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMLinearFormKernel>>> _lf_kernels_map;

  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMNonlinearFormKernel>>> _nlf_kernels_map;

  platypus::NamedFieldsMap<
      platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>>
      _mblf_kernels_map_map;

  //mutable mfem::OperatorHandle _jacobian;

  // Variables for time-dependent equation systems
  mfem::ConstantCoefficient _dt_coef; // Coefficient for timestep scaling
  std::vector<std::string> _trial_var_time_derivative_names;

};


} // namespace platypus

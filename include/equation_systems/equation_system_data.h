#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"

namespace platypus
{

using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

/*
Struct to hold all of the data related to the equation system needed for the solver to use in static
systems
*/
struct EquationSystemData
{

  friend class EquationSystemOperatorBase;
  friend class EquationSystemOperator;
  friend class TimeDependentEquationSystemOperator;
  friend class EquationSystemProblemOperator;
  friend class TimeDomainEquationSystemProblemOperator;

protected:
  enum MatrixType
  {
    DIAGONAL,
    DENSE
  };

  MatrixType _matrix_type;
  mfem::AssemblyLevel _assembly_level;

  // Test variables are associated with LinearForms,
  // whereas trial variables are associated with gridfunctions.

  // Names of all variables corresponding to gridfunctions. This may differ
  // from test_var_names when time derivatives are present.
  std::vector<std::string> _trial_var_names;
  // Pointers to trial variables.
  platypus::GridFunctions _trial_variables;
  // Names of all test variables corresponding to linear forms in this equation
  // system
  std::vector<std::string> _test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> _test_pfespaces;

  // Components of weak form. // Named according to test variable
  platypus::NamedFieldsMap<mfem::ParBilinearForm> _blfs;
  platypus::NamedFieldsMap<mfem::ParLinearForm> _lfs;
  platypus::NamedFieldsMap<mfem::ParNonlinearForm> _nlfs;
  platypus::NamedFieldsMap<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _mblfs; // named according to trial variable

  std::vector<mfem::Array<int>> _ess_tdof_lists;

  // gridfunctions for setting Dirichlet BCs
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _xs;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _dxdts;

  mfem::Array2D<mfem::HypreParMatrix *> _h_blocks;

  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>> _blf_kernels_map;

  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMLinearFormKernel>>> _lf_kernels_map;

  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMNonlinearFormKernel>>> _nlf_kernels_map;

  platypus::NamedFieldsMap<
      platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>>
      _mblf_kernels_map_map;

  mutable mfem::OperatorHandle _jacobian;
};

/*
Struct to hold all of the data related to the equation system needed for the solver to use in time
dependent systems
*/
struct TimeDependentEquationSystemData : public EquationSystemData
{

  friend class EquationSystemOperatorBase;
  friend class EquationSystemOperator;
  friend class TimeDependentEquationSystemOperator;
  friend class EquationSystemProblemOperator;
  friend class TimeDomainEquationSystemProblemOperator;

protected:
  mfem::ConstantCoefficient _dt_coef; // Coefficient for timestep scaling
  std::vector<std::string> _trial_var_time_derivative_names;

  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>>
      _td_blf_kernels_map;
  // Container to store contributions to weak form of the form (F du/dt, v)
  platypus::NamedFieldsMap<mfem::ParBilinearForm> _td_blfs;
};

} // namespace platypus
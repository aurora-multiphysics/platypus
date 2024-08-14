#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"
#include "equation_system_modifier.h"

namespace platypus
{
/*
Operator containing the necessary methods to supply to the solver
*/

class EquationSystemOperator : public mfem::Operator
{

public:
  EquationSystemOperator() = default;
  ~EquationSystemOperator() = default;

  void Init(std::shared_ptr<EquationSystemData> data);
  void Mult(const mfem::Vector & u, mfem::Vector & residual);
  mfem::Operator & GetGradient(const mfem::Vector & u);
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions);

  [[nodiscard]] EquationSystemData * getData()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }

    return _equation_system_data.get();
  }

protected:
  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};
};

} // namespace platypus

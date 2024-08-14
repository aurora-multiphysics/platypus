#include "equation_system_operator.h"

namespace platypus
{

void
EquationSystemOperator::Init(std::shared_ptr<EquationSystemData> data)
{
  _equation_system_data = data;
  height = getData()->_jacobian_height;
  width = getData()->_jacobian_width;
}

void
EquationSystemOperator::Mult(const mfem::Vector & x, mfem::Vector & residual)
{
  getData()->_jacobian->Mult(x, residual);
}

mfem::Operator &
EquationSystemOperator::GetGradient(const mfem::Vector & u)
{
  return *(getData()->_jacobian);
}

void
EquationSystemOperator::RecoverFEMSolution(mfem::BlockVector & trueX,
                                           platypus::GridFunctions & gridfunctions)
{
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto & test_var_name = getData()->_test_var_names.at(i);
    trueX.GetBlock(i).SyncAliasMemory(trueX);
    gridfunctions.Get(test_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
}

} // namespace platypus
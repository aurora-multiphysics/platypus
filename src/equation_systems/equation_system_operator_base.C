#include "equation_system_operator_base.h"

namespace platypus
{

// Base methods
bool
EquationSystemOperatorBase::VectorContainsName(const std::vector<std::string> & the_vector,
                                           const std::string & name) const
{

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

} // namespace platypus
#pragma once
#include "MFEMMaterial.h"
#include "MFEMContainers.h"

/**
 * Declares material properties based on names and functions prescribed by input parameters.
 *
 * This is identical in function to the GenericFunctionMaterial in Moose.
 */
class MFEMGenericFunctorMaterial : public MFEMMaterial
{
public:
  static InputParameters validParams();

  MFEMGenericFunctorMaterial(const InputParameters & parameters);
  virtual ~MFEMGenericFunctorMaterial();

protected:
  const std::vector<std::string> & _prop_names;
  const std::vector<MFEMScalarCoefficientName> & _prop_values;
  unsigned int _num_props;
};

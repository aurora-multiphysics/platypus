#include "MFEMGenericFunctorMaterial.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMGenericFunctorMaterial);

InputParameters
MFEMGenericFunctorMaterial::validParams()
{
  InputParameters params = MFEMMaterial::validParams();
  params.addClassDescription("Declares material scalar properties based on names and functors "
                             "prescribed by input parameters.");
  params.addRequiredParam<std::vector<std::string>>(
      "prop_names", "The names of the properties this material will have");
  params.addRequiredParam<std::vector<MFEMScalarCoefficientName>>(
      "prop_values",
      "The corresponding names of functors associated with the named properties. A functor is any "
      "of the following: a variable, an MFEM material property, a function, or a post-processor.");

  return params;
}

MFEMGenericFunctorMaterial::MFEMGenericFunctorMaterial(const InputParameters & parameters)
  : MFEMMaterial(parameters),
    _prop_names(getParam<std::vector<std::string>>("prop_names")),
    _prop_values(getParam<std::vector<MFEMScalarCoefficientName>>("prop_values"))
{
  unsigned int num_names = _prop_names.size();
  unsigned int num_values = _prop_values.size();

  if (num_names != num_values)
    mooseError(
        "Number of prop_names must match the number of prop_values for a GenericConstantMaterial!");

  _num_props = num_names;
  for (unsigned int i = 0; i < _num_props; i++)
  {
    _properties.declareScalarProperty(
        _prop_names[i], subdomainsToStrings(_block_ids), _prop_values[i]);
  }
}

MFEMGenericFunctorMaterial::~MFEMGenericFunctorMaterial() {}

#pragma once
#include "MFEMMaterial.h"
#include "FunctionParserUtils.h"

/**
 * Declares material properties based on names and values prescribed by input parameters.
 *
 * This is identical in function to the GenericConstantMaterial in Moose.
 */
class MFEMParsedCoeffMaterial : public MFEMMaterial, public FunctionParserUtils<false>
{
public:
  static InputParameters validParams();

  MFEMParsedCoeffMaterial(const InputParameters & parameters);
  virtual ~MFEMParsedCoeffMaterial();

protected:
  /// function expression
  std::string _function;
  const std::vector<std::string> & _var_names;
  const std::vector<std::string> & _prop_names;
  const std::vector<Real> & _prop_values;
  unsigned int _num_props;
  /// import coordinates and time
  const bool _use_xyzt;
  /// coordinate and time variable names
  const std::vector<std::string> _xyzt;
  const MFEMProblemData & _problem_data;
  /// function parser object for the resudual and on-diagonal Jacobian
  SymFunctionPtr _func_F;

};

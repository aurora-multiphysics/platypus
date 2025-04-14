#include "MFEMParsedCoeffMaterial.h"
#include "MFEMScalarParsedCoeff.h"
#include <vector>
#include <memory>

registerMooseObject("PlatypusApp", MFEMParsedCoeffMaterial);

InputParameters
MFEMParsedCoeffMaterial::validParams()
{
  InputParameters params = MFEMMaterial::validParams();
  params += FunctionParserUtils<false>::validParams();
  params.addClassDescription("Declares constant material properties based on names and values "
                             "prescribed by input parameters.");
  params.addRequiredCustomTypeParam<std::string>(
                              "function", "FunctionExpression", "Parsed function expression to compute");
  params.deprecateParam("function", "expression", "02/07/2024");
  params.addRequiredParam<std::vector<std::string>>(
    "var_names", "The names of the function variable names");
  params.addParam<bool>(
    "use_xyzt",
    false,
    "Make coordinate (x,y,z) and time (t) variables available in the function expression.");
  return params;
}

MFEMParsedCoeffMaterial::MFEMParsedCoeffMaterial(const InputParameters & parameters)
  : MFEMMaterial(parameters),
    FunctionParserUtils(parameters),
    _function(getParam<std::string>("expression")),    
    _var_names(getParam<std::vector<std::string>>("var_names")),
    _use_xyzt(getParam<bool>("use_xyzt")),
    _xyzt({"x", "y", "z", "t"}),
    _problem_data(getMFEMProblem().getProblemData())
{

    // build variables argument
    std::string variables;

    // coupled field variables
    for (const auto i : index_range(_var_names))
      variables += (i == 0 ? "" : ",") + _var_names[i];

    // positions and time
    if (_use_xyzt)
      for (auto & v : _xyzt)
        variables += (variables.empty() ? "" : ",") + v;

    // base function object
    _func_F = std::make_shared<SymFunction>();

    // set FParser internal feature flags
    setParserFeatureFlags(_func_F);

    // parse function
    if (_func_F->Parse(_function, variables) >= 0)
        mooseError(
            "Invalid function\n", _function, "\nin ParsedAux ", name(), ".\n", _func_F->ErrorMsg());

    // optimize
    if (!_disable_fpoptimizer)
      _func_F->Optimize();

    // just-in-time compile
    if (_enable_jit)
      {
        // let rank 0 do the JIT compilation first
        if (_communicator.rank() != 0)
          _communicator.barrier();

        _func_F->JITCompile();

        // wait for ranks > 0 to catch up
        if (_communicator.rank() == 0)
          _communicator.barrier();
      }

    // reserve storage for parameter passing buffer
    _func_params.resize(_var_names.size() + (_use_xyzt ? 4 : 0));

  //std::shared_ptr<mfem::ParGridFunction> _shared_grid_function = _problem_data._gridfunctions.GetShared("hello");
  //mfem::ParGridFunction & _var_grid_function = _problem_data._gridfunctions.GetRef("hello");

  //using GFMapType = std::vector<std::shared_ptr<mfem::ParGridFunction>>;
  // GFMapType grid_functions;

  // for (unsigned int i = 0; i < _var_names.size(); i++)
  // grid_functions.push_back(_problem_data._gridfunctions.GetShared(_var_names[i]));

  // MFEMScalarParsedCoeff( _problem_data._gridfunctions, _var_names
  //  , std::function<double(std::vector<double>)> func_);

    _properties.declareScalar<MFEMScalarParsedCoeff>(
        "pared_material", subdomainsToStrings(_block_ids), _problem_data._gridfunctions, _var_names, _func_F);

}

MFEMParsedCoeffMaterial::~MFEMParsedCoeffMaterial() {}

#include "MFEMEstimator.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMEstimator);

// static method
InputParameters
MFEMEstimator::validParams()
{
  InputParameters params = MFEMGeneralUserObject::validParams();
  params.registerBase("Estimator");
  return params;
}


MFEMEstimator::MFEMEstimator(const InputParameters & params)
  : MFEMGeneralUserObject(params), _test_var_name(getParam<std::string>("test_variable")),
    _kernel_name(getParam<std::string>("kernel_name"))
{}

// static method
std::shared_ptr<MFEMEstimator>
MFEMEstimator::setUp(std::string estimator_type, std::string estimator_name, InputParameters estimator_params)
{
  getMFEMProblem().addUserObject( estimator_type, estimator_name, estimator_params );
  const UserObject * est_uo = &(getMFEMProblem().getUserObjectBase(estimator_name));

  if ( dynamic_cast<const MFEMEstimator *>(est_uo) != nullptr )
  {
    std::shared_ptr<MooseObject> object_ptr  = getMFEMProblem().getUserObject<MFEMEstimator>(estimator_name).getSharedPtr();
    std::shared_ptr<MFEMEstimator> estimator = std::dynamic_pointer_cast<MFEMEstimator>(object_ptr);
    return std::move(estimator);
  }

  else
  {
    mooseError("Cannot add ", estimator_type);
    return {};
  }
}

#include "MFEMVectorBoundaryIntegratedBC.h"

registerMooseObject("PlatypusApp", MFEMVectorBoundaryIntegratedBC);

InputParameters
MFEMVectorBoundaryIntegratedBC::validParams()
{
  InputParameters params = MFEMBoundaryCondition::validParams();
  params.addRequiredParam<UserObjectName>(
      "vector_coefficient", "The vector MFEM coefficient which will be used in the integrated BC");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorBoundaryIntegratedBC::MFEMVectorBoundaryIntegratedBC(const InputParameters & parameters)
  : MFEMBoundaryCondition(parameters),
    _vec_coef(const_cast<MFEMVectorCoefficient *>(
        &getUserObject<MFEMVectorCoefficient>("vector_coefficient")))
{

  _boundary_condition = std::make_shared<platypus::IntegratedBC>(
      getParam<std::string>("variable"),
      bdr_attr,
      std::make_unique<mfem::VectorBoundaryLFIntegrator>(*_vec_coef->getVectorCoefficient()));
}

#pragma once

#include "MFEMGeneralUserObject.h"
#include "MFEMContainers.h"

/*
Class to construct an MFEM integrator to apply to the equation system.
*/
class MFEMKernel : public MFEMGeneralUserObject
{
public:
  static InputParameters validParams();

  MFEMKernel(const InputParameters & parameters);

  virtual ~MFEMKernel() = default;

  // Create a new MFEM integrator to apply to the weak form. Ownership managed by the caller.
  virtual mfem::LinearFormIntegrator * createResidualContribution() { return nullptr; }
  virtual mfem::BilinearFormIntegrator * createJacobianContribution() { return nullptr; }

  // Get name of the test variable labelling the weak form this kernel is added to
  const std::string & getTestVariableName() const { return _test_var_name; }

  // Get name of the trial variable (gridfunction) the kernel acts on.
  // Defaults to the name of the test variable labelling the weak form.
  virtual const std::string & getTrialVariableName() const { return _test_var_name; }

  const bool isSubdomainRestricted() { return _subdomain_names.size(); }

  mfem::Array<int> & getSubdomains() { return _subdomain_attributes; }

protected:
  // Name of (the test variable associated with) the weak form that the kernel is applied to.
  std::string _test_var_name;
  std::vector<SubdomainName> _subdomain_names;
  mfem::Array<int> _subdomain_attributes;
};

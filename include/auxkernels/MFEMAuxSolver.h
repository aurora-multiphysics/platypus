#pragma once

#include "MFEMGeneralUserObject.h"
#include "coefficients.h"
#include "auxsolver_base.h"

class MFEMAuxSolver : public MFEMGeneralUserObject
{
public:
  static InputParameters validParams();

  MFEMAuxSolver(const InputParameters & parameters);
  virtual ~MFEMAuxSolver();

  inline virtual std::shared_ptr<platypus::AuxSolver> getAuxSolver() const { return _auxsolver; }

  virtual void storeCoefficients(platypus::Coefficients & coefficients) {}

protected:
  std::shared_ptr<platypus::AuxSolver> _auxsolver{nullptr};
};

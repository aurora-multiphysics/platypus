#pragma once

#include "GeneralUserObject.h"
#include "problem_builder_base.h"

class MFEMFormulation : public GeneralUserObject
{
public:
  static InputParameters validParams();

  MFEMFormulation(const InputParameters & parameters);
  virtual ~MFEMFormulation();

  virtual void execute() override {}
  virtual void initialize() override {}
  virtual void finalize() override {}

  virtual std::shared_ptr<platypus::ProblemBuilder> getProblemBuilder()
  {
    mooseError(
        "Base class MFEMFormulation cannot return a valid ProblemBuilder. Use a child class.");
  }
};

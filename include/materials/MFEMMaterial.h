#pragma once

#include <mfem.hpp>
#include "MFEMGeneralUserObject.h"
#include "CoefficientManager.h"

class MFEMMaterial : public MFEMGeneralUserObject
{
public:
  static InputParameters validParams();
  static std::vector<std::string> subdomainsToStrings(std::vector<SubdomainName> blocks);
  static libMesh::Point pointFromMFEMVector(const mfem::Vector & vec);

  MFEMMaterial(const InputParameters & parameters);
  virtual ~MFEMMaterial();

  const std::vector<SubdomainName> & getBlocks() const { return _block_ids; }

protected:
  std::vector<SubdomainName> _block_ids;
  platypus::CoefficientManager & _properties;
};

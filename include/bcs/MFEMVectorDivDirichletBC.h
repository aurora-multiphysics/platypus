#pragma once

#include "MFEMVectorDirichletBCBase.h"

class MFEMVectorDivDirichletBC : public MFEMVectorDirichletBCBase
{
public:
  MFEMVectorDivDirichletBC(const InputParameters & parameters);
  ~MFEMVectorDivDirichletBC() override = default;
};

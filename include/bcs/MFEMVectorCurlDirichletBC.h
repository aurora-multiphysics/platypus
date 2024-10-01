#pragma once

#include "MFEMVectorDirichletBCBase.h"

class MFEMVectorCurlDirichletBC : public MFEMVectorDirichletBCBase
{
public:
  MFEMVectorCurlDirichletBC(const InputParameters & parameters);
  ~MFEMVectorCurlDirichletBC() override = default;
};

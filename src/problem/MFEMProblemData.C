#include "MFEMProblemData.h"

void
MFEMProblemData::updateFESpaces()
{
  for (const auto & fe_space_pair : _fespaces)
  {
    fe_space_pair.second->Update();
  }
  for (const auto & gridfunction_pair : _gridfunctions)
  {
    gridfunction_pair.second->Update();
  }
  _eqn_system->UpdateEquationSystem();
}

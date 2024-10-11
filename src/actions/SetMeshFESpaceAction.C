#include "SetMeshFESpaceAction.h"

registerMooseAction("PlatypusApp", SetMeshFESpaceAction, "set_mesh_fe_space");

InputParameters
SetMeshFESpaceAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set the mesh nodal finite element space to the same as the mesh "
                             "displacement variable, if one is specified.");
  return params;
}

SetMeshFESpaceAction::SetMeshFESpaceAction(const InputParameters & parameters) : Action(parameters)
{
}

void
SetMeshFESpaceAction::act()
{
  auto & mfem_problem = static_cast<MFEMProblem &>(*_problem.get());
  mfem::ParMesh & mesh = mfem_problem.mesh().getMFEMParMesh();
  mfem::ParGridFunction const * displacement = mfem_problem.getMeshDisplacementGridFunction();
  if (displacement)
  {
    mesh.SetNodalFESpace(displacement->ParFESpace());
  }
}

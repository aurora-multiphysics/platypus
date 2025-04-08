#pragma once

#include "MFEMFESpace.h"
#include "MFEMGeneralUserObject.h"
#include "MooseFunctor.h"

/**
 * Constructs and stores an mfem::ParGridFunction object.
 */
class MFEMVariable : public MFEMGeneralUserObject, private Moose::FunctorBase<libMesh::Real>
{
public:
  static InputParameters validParams();

  MFEMVariable(const InputParameters & parameters);

  /// Returns a shared pointer to the constructed gridfunction.
  inline std::shared_ptr<mfem::ParGridFunction> getGridFunction() const { return _gridfunction; }

  libMesh::Real evaluate(const Moose::ElemArg & elem, const Moose::StateArg & state) const
  {
    return 0.0;
  }
  libMesh::Real evaluate(const Moose::FaceArg & face, const Moose::StateArg & state) const
  {
    return 0.0;
  }
  libMesh::Real evaluate(const Moose::ElemQpArg & qp, const Moose::StateArg & state) const
  {
    return 0.0;
  }
  libMesh::Real evaluate(const Moose::ElemSideQpArg & side_qp, const Moose::StateArg & state) const
  {
    return 0.0;
  }
  libMesh::Real evaluate(const Moose::ElemPointArg & elem_point, const Moose::StateArg & state) const
  {
    return 0.0;
  }
  libMesh::Real evaluate(const Moose::NodeArg & node, const Moose::StateArg & state) const
  {
    return 0.0;
  }

  bool supportsFaceArg() const
  {
    return false;
  }
  bool supportsElemSideQpArg() const
  {
    return false;
  }

protected:
  const MFEMFESpace & _fespace;

private:
  /// Constructs the gridfunction.
  const std::shared_ptr<mfem::ParGridFunction> buildGridFunction();

  /// Stores the constructed gridfunction.
  const std::shared_ptr<mfem::ParGridFunction> _gridfunction{nullptr};
    
};

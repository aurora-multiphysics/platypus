#pragma once
#include "MFEMGeneralUserObject.h"
#include "mfem/miniapps/common/pfem_extras.hpp"

namespace platypus
{

/// Integrator which scales its results by a constant value
class MFEMCustomIntegrator : public mfem::BilinearFormIntegrator
{
public:
  static InputParameters validParams();

  MFEMCustomIntegrator(mfem::Coefficient & q);

  // Real
  // MFEMCustomIntegrator::computeQpResidual();

  mfem::real_t computeQpJacobian();

  // void AssembleElementVector(const FiniteElement &el,
  //                         ElementTransformation &Tr,
  //                         const Vector &elfun, Vector &elvect) override;

  void AssembleElementMatrix(const mfem::FiniteElement & el,
                             mfem::ElementTransformation & Trans,
                             mfem::DenseMatrix & _local_ke) override;

  void AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                              const mfem::FiniteElement & test_fe,
                              mfem::ElementTransformation & Trans,
                              mfem::DenseMatrix & _local_ke) override;

  const mfem::IntegrationRule & GetRule(const mfem::FiniteElement & trial_fe,
                                        const mfem::FiniteElement & test_fe,
                                        mfem::ElementTransformation & Trans);

protected:
  mfem::Coefficient * Q;
  mfem::VectorCoefficient * VQ;
  mfem::MatrixCoefficient * MQ;

private:
  mfem::Vector vec, vecdxt, pointflux;
  mfem::Vector _test, _phi;
  mfem::real_t _q, _JxW;
  unsigned int _i, _j;
#ifndef MFEM_THREAD_SAFE
  mfem::DenseMatrix dshape, dshapedxt, invdfdx, M, dshapedxt_m;
  mfem::DenseMatrix te_dshape, te_dshapedxt;
  mfem::Vector D;
#endif
};

// ComputeQPResidual
// _test[_i][_qp] * (_velocity * _grad_u[_qp]);
// ComputeJac
// _test[_i][_qp] * (_velocity * _grad_phi[_j][_qp]);

// BdFidxT = Q \cdot \nabla u = _velocity * _grad_phi[_j][_qp]) * weight
// dshape[j] = _grad_test[_j][_qp])
// shape[i] = v = _test[_i][_qp]
// elfun = grad_u[qp] = sum(c[j] _grad_phi[j])

// DenseMatrix _local_ke;
// AssembleElementMatrix(el, Tr, _local_ke);
// elvect.SetSize(_local_ke.Height());
// _local_ke.Mult(elfun, elvect); -

// AssembleElementMatrix2(
//    const FiniteElement &trial_fe, const FiniteElement &test_fe,
//    ElementTransformation &Trans, DenseMatrix &_local_ke)

} // namespace platypus

#pragma once
#include "MFEMGeneralUserObject.h"
#include "MFEMIntegratorInterface.h"
#include "mfem/miniapps/common/pfem_extras.hpp"

namespace platypus
{

/// Integrator which scales its results by a constant value
class MFEMCustomIntegrator : public mfem::BilinearFormIntegrator
{
public:
  static InputParameters validParams();

  MFEMCustomIntegrator(platypus::MFEMIntegratorInterface * mfem_ii);

  // Real
  // MFEMCustomIntegrator::computeQpResidual();

  // void AssembleElementVector(const FiniteElement &el,
  //                         ElementTransformation &Tr,
  //                         const Vector &elfun, Vector &elvect) override;

  void AssembleElementMatrix(const mfem::FiniteElement & el,
                             mfem::ElementTransformation & Trans,
                             mfem::DenseMatrix & local_ke) override;

  void AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                              const mfem::FiniteElement & test_fe,
                              mfem::ElementTransformation & Trans,
                              mfem::DenseMatrix & local_ke) override;

  const mfem::IntegrationRule & GetRule(const mfem::FiniteElement & trial_fe,
                                        const mfem::FiniteElement & test_fe,
                                        mfem::ElementTransformation & Trans);

protected:
  platypus::MFEMIntegratorInterface * _mfem_ii;
  mfem::Coefficient * Q;
  mfem::VectorCoefficient * VQ;
  mfem::MatrixCoefficient * MQ;

private:
  mfem::Vector &_test, &_phi;
  unsigned int &_i, &_j;
  mfem::real_t & _JxW;
};

// ComputeQPResidual
// _test[_i][_qp] * (_velocity * _grad_u[_qp]);
// ComputeJac
// _test[_i][_qp] * (_velocity * _grad_phi[_j][_qp]);

// BdFidxT = Q \cdot \nabla u = _velocity * _grad_phi[_j][_qp]) * weight
// dshape[j] = _grad_test[_j][_qp])
// shape[i] = v = _test[_i][_qp]
// elfun = grad_u[qp] = sum(c[j] _grad_phi[j])

// DenseMatrix local_ke;
// AssembleElementMatrix(el, Tr, local_ke);
// elvect.SetSize(local_ke.Height());
// local_ke.Mult(elfun, elvect); -

// AssembleElementMatrix2(
//    const FiniteElement &trial_fe, const FiniteElement &test_fe,
//    ElementTransformation &Trans, DenseMatrix &local_ke)

} // namespace platypus

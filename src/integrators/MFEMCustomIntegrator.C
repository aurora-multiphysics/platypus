#include "MFEMCustomIntegrator.h"

namespace platypus
{

MFEMCustomIntegrator::MFEMCustomIntegrator(mfem::Coefficient & q)
  : mfem::BilinearFormIntegrator(), Q(&q)
{
}

// Real
// MFEMCustomIntegrator::computeQpResidual()
// {
//   // velocity * _grad_u[_qp] is actually doing a dot product
//   // return _test[_i][_qp] * (_velocity * _grad_u[_qp]);
//   return _test[_i][_qp] * (_velocity * _grad_u[_qp]);
// }

mfem::real_t
MFEMCustomIntegrator::computeQpJacobian(int & i, int & j)
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
  // return _test[_i][_qp] * (_velocity * _grad_phi[_j][_qp]);

  return _test(i) * _q * _phi(j);
}

void
MFEMCustomIntegrator::AssembleElementMatrix(const mfem::FiniteElement & el,
                                            mfem::ElementTransformation & Trans,
                                            mfem::DenseMatrix & elmat)
{
  AssembleElementMatrix2(el, el, Trans, elmat);
}

void
MFEMCustomIntegrator::AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                                             const mfem::FiniteElement & test_fe,
                                             mfem::ElementTransformation & Trans,
                                             mfem::DenseMatrix & elmat)
{
  int tr_nd = trial_fe.GetDof();
  int te_nd = test_fe.GetDof();
  mfem::real_t w;
  elmat.SetSize(te_nd, tr_nd);
  _phi.SetSize(tr_nd);
  _test.SetSize(te_nd);

  const mfem::IntegrationRule * ir = IntRule ? IntRule : &GetRule(trial_fe, test_fe, Trans);

  elmat = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++)
  {
    const mfem::IntegrationPoint & ip = ir->IntPoint(i);
    Trans.SetIntPoint(&ip);

    trial_fe.CalcPhysShape(Trans, _phi);
    test_fe.CalcPhysShape(Trans, _test);
    // trial_fe.CalcDShape(ip, _grad_phi); // _grad_phi
    // test_fe.CalcDShape(ip, _grad_test); //_grad_test

    // CalcAdjugate(Trans.Jacobian(), invdfdx);
    // w = Trans.Weight();
    // w = ip.weight / (square ? w : w*w*w);
    // Mult(_grad_phi, invdfdx, _grad_phidxt);
    // Mult(_grad_test, invdfdx, _grad_testdxt);

    w = Trans.Weight() * ip.weight;
    if (Q)
    {
      _q = Q->Eval(Trans, ip);
    }
    else
    {
      _q = 1.0;
    }
    _test *= w;
    const int test_size = _test.Size();
    const int phi_size = _phi.Size();

    const int m = _test.Size(), n = _phi.Size();
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        elmat(i, j) += computeQpJacobian(i, j);
      }
    }
  }
}

//  void AssembleElementVector(const FiniteElement &el,
//                             ElementTransformation &Tr,
//                             const Vector &elfun, Vector &elvect) override;

const mfem::IntegrationRule &
MFEMCustomIntegrator::GetRule(const mfem::FiniteElement & trial_fe,
                              const mfem::FiniteElement & test_fe,
                              mfem::ElementTransformation & Trans)
{
  // int order = trial_fe.GetOrder() + test_fe.GetOrder();
  const int order = trial_fe.GetOrder() + test_fe.GetOrder() + Trans.OrderW();

  if (trial_fe.Space() == mfem::FunctionSpace::rQk)
  {
    return mfem::RefinedIntRules.Get(trial_fe.GetGeomType(), order);
  }
  return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
}

} // namespace platypus
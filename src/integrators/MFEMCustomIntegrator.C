#include "MFEMCustomIntegrator.h"

namespace platypus
{

MFEMCustomIntegrator::MFEMCustomIntegrator(platypus::MFEMIntegratorInterface * mfem_ii)
  : mfem::BilinearFormIntegrator(),
    _mfem_ii(mfem_ii),
    _test(mfem_ii->_test),
    _phi(mfem_ii->_phi),
    _i(mfem_ii->_i),
    _j(mfem_ii->_j),
    _JxW(mfem_ii->_JxW)
{
}

// Real
// MFEMCustomIntegrator::computeQpResidual()
// {
//   // velocity * _grad_u[_qp] is actually doing a dot product
//   // return _test[_i][_qp] * (_velocity * _grad_u[_qp]);
//   return _test[_i][_qp] * (_velocity * _grad_u[_qp]);
// }

// mfem::real_t
// MFEMCustomIntegrator::computeQpJacobian()
// {
//   // the partial derivative of _grad_u is just _grad_phi[_j]
//   // return _test[_i][_qp] * (_velocity * _grad_phi[_j][_qp]);

//   return _test[_i] * _phi[_j];
// }

void
MFEMCustomIntegrator::AssembleElementMatrix(const mfem::FiniteElement & el,
                                            mfem::ElementTransformation & Trans,
                                            mfem::DenseMatrix & local_ke)
{
  AssembleElementMatrix2(el, el, Trans, local_ke);
}

void
MFEMCustomIntegrator::AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                                             const mfem::FiniteElement & test_fe,
                                             mfem::ElementTransformation & Trans,
                                             mfem::DenseMatrix & local_ke)
{
  int tr_nd = trial_fe.GetDof();
  int te_nd = test_fe.GetDof();
  local_ke.SetSize(te_nd, tr_nd);
  _phi.SetSize(tr_nd);
  _test.SetSize(te_nd);

  const mfem::IntegrationRule * ir = IntRule ? IntRule : &GetRule(trial_fe, test_fe, Trans);

  local_ke = 0.0;
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
    // if (Q)
    // {
    //   _q = Q->Eval(Trans, ip);
    // }
    // else
    // {
    //   _q = 1.0;
    // }

    _JxW = Trans.Weight() * ip.weight;
    const unsigned int test_size = _test.Size();
    const unsigned int phi_size = _phi.Size();
    for (_i = 0; _i < test_size; _i++)
    {
      for (_j = 0; _j < phi_size; _j++)
      {
        local_ke(_i, _j) += _JxW * _mfem_ii->computeQpJacobian();
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
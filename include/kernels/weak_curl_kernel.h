#pragma once
#include "kernel_base.h"

namespace platypus
{

/*
(αv_{n}, ∇×u')
*/
class WeakCurlKernel : public Kernel<mfem::ParLinearForm>
{
public:
  WeakCurlKernel(const platypus::InputParameters & params);

  ~WeakCurlKernel() override = default;

  void Init(platypus::GridFunctions & gridfunctions,
            const platypus::FESpaces & fespaces,
            platypus::BCMap & bc_map,
            platypus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * lf) override;

  std::string _hcurl_gf_name, _hdiv_gf_name;
  std::string _coef_name;
  mfem::ParGridFunction * _u{nullptr};
  mfem::ParGridFunction * _v{nullptr};
  mfem::Coefficient * _coef{nullptr};

  std::unique_ptr<mfem::ParMixedBilinearForm> _weak_curl;
};

} // namespace platypus

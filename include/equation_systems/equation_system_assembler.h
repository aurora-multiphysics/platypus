#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"

namespace platypus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystemAssembler
{
public:
  using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
  using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
  using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
  using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

  EquationSystemAssembler() = default;
  ~EquationSystemAssembler() = default;


  [[nodiscard]] EquationSystemData * getData()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }

    return _equation_system_data.get();

  };

private:

  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};

};


} // namespace platypus

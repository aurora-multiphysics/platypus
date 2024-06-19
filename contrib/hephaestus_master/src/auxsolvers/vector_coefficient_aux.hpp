#pragma once
#include "coefficient_aux.hpp"

namespace hephaestus
{

// Project a stored vector Coefficient onto a (vector) GridFunction
class VectorCoefficientAux : public AuxSolver
{
public:
  VectorCoefficientAux(std::string gf_name,
                       std::string vec_coef_name,
                       hephaestus::InputParameters solver_options = hephaestus::InputParameters());

  ~VectorCoefficientAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  virtual void BuildBilinearForm();
  virtual void BuildLinearForm();
  void Solve(double t = 0.0) override;

protected:
  const std::string _gf_name;       // name of the variable
  const std::string _vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction * _gf{nullptr};
  mfem::VectorCoefficient * _vec_coef{nullptr};

  // Pointer to store test FE space. Assumed to be same as trial FE space.
  mfem::ParFiniteElementSpace * _test_fes{nullptr};

  // Bilinear and linear forms
  std::unique_ptr<mfem::ParBilinearForm> _a{nullptr};
  std::unique_ptr<mfem::ParLinearForm> _b{nullptr};

private:
  const hephaestus::InputParameters _solver_options;

  // Operator matrices
  std::unique_ptr<mfem::HypreParMatrix> _a_mat{nullptr};

  // Solver
  std::unique_ptr<hephaestus::DefaultJacobiPCGSolver> _solver{nullptr};
};

} // namespace hephaestus

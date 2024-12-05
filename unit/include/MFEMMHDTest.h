#include "time_domain_problem_operator.h"
#include "MFEMProblemData.h"

class LiMHDPreconditioner : public platypus::TimeDomainProblemOperator
{
public:
  LiMHDPreconditioner(MFEMProblemData & problem) : platypus::TimeDomainProblemOperator(problem) {}
  void Assemble();
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}
  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
};
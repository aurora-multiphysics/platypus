/***********************************************\
!        Construct a non-linear parsed          !
!        grid function coefficient              !
!                                               !
\***********************************************/

// Define a coefficient that, given a grid function u,
// function func, returns func(u)

#include "mfem.hpp"
#include "MFEMContainers.h"
#include <vector>
#include <functional>

class MFEMScalarParsedCoeff : public mfem::Coefficient
{
  private:
    platypus::GridFunctions  & gFuncs;
    std::vector<std::string> & inputs;
    std::function<double(std::vector<double>)> func;

  public:
    MFEMScalarParsedCoeff( platypus::GridFunctions  & gFuncs_
                 , std::vector<std::string> & inputs_
                 , std::function<double(std::vector<double>)> func_);

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

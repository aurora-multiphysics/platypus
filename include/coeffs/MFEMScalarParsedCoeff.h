/***********************************************\
!        Construct a non-linear parsed          !
!        grid function coefficient              !
!                                               !
\***********************************************/

// Define a coefficient that, given a grid function u,
// function func, returns func(u)

#include "mfem.hpp"
#include "MFEMContainers.h"
#include "FunctionParserUtils.h"
#include <vector>
#include <functional>



class MFEMScalarParsedCoeff : public mfem::Coefficient
{
  private:
    const platypus::GridFunctions  & gFuncs;
    const std::vector<std::string> & inputs;
    const FunctionParserUtils<false>::SymFunctionPtr & func;

  public:
    MFEMScalarParsedCoeff(const platypus::GridFunctions  & gFuncs_
                 , const std::vector<std::string> & inputs_
                 , const FunctionParserUtils<false>::SymFunctionPtr & func_);

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

 #include "MFEMScalarParsedCoeff.h"


MFEMScalarParsedCoeff::MFEMScalarParsedCoeff( platypus::GridFunctions  & gFuncs_
    , std::vector<std::string> & inputs_
    , FunctionParserUtils<false>::SymFunctionPtr & func_)
    : gFuncs(gFuncs_), inputs(inputs_), func(func_){}

double  MFEMScalarParsedCoeff::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
    std::vector<libMesh::Real> InpVals(inputs.size() + 4);
    
    mfem::Vector transip(3);
    T.Transform(ip, transip);

    for(unsigned i=0; i<inputs.size(); i++)
      InpVals[i] = gFuncs.GetRef(inputs[i]).GetValue(T, ip);

    for(unsigned i=0; i<3; i++)
      InpVals[inputs.size()+i] = transip(i);
    
    InpVals[inputs.size()+3] = GetTime();
    
    return func->Eval(InpVals.data());
}
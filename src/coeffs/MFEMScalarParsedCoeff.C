 #include "MFEMScalarParsedCoeff.h"


MFEMScalarParsedCoeff::MFEMScalarParsedCoeff( platypus::GridFunctions  & gFuncs_
    , std::vector<std::string> & inputs_
    , std::function<double(std::vector<double>)> func_)
    : gFuncs(gFuncs_), inputs(inputs_), func(func_){}

double  MFEMScalarParsedCoeff::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
    std::vector<double> InpVals(inputs.size() + 4);   
    for(unsigned I=0; I<inputs.size(); I++)
    {
        //if(){//Checks if its is a scalar GridFunction
        //  InpVals[I] = gFuncs[inputs[I]].GetValue(T, ip);
    }
    
    return func(InpVals);
}
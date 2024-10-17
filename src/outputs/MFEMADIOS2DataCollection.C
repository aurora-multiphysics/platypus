#include "MFEMADIOS2DataCollection.h"
#include "FEProblemBase.h"
#include "MFEMMesh.h"

registerMooseObject("PlatypusApp", MFEMADIOS2DataCollection);

InputParameters
MFEMADIOS2DataCollection::validParams()
{
  InputParameters params = MFEMDataCollection::validParams();
  params.addClassDescription("Output for controlling MFEMADIOS2DataCollection inherited data.");
  params.addParam<unsigned int>("refinements",
                                0,
                                "Number of uniform refinements for oversampling "
                                "(refinement levels beyond any uniform "
                                "refinements)");
  MooseEnum engine_types("BP3 BP4 BP5", "BP3", false);
  params.addParam<MooseEnum>(
      "engine_type", engine_types, "Choice of ADIOS2 engine to be used for output.");
  return params;
}

MFEMADIOS2DataCollection::MFEMADIOS2DataCollection(const InputParameters & parameters)
  : MFEMDataCollection(parameters),
    _refinements(getParam<unsigned int>("refinements")),
    _engine_type(parameters.get<MooseEnum>("engine_type")),
    _adios_dc(_problem_data._comm,
              (_file_base + std::string("/Run") + std::to_string(getFileNumber())).c_str(),
              _problem_data._pmesh.get(),
              _engine_type)
{
  _adios_dc.SetPrecision(9);
  _adios_dc.SetLevelsOfDetail(_refinements + 1);
  registerFields();
}

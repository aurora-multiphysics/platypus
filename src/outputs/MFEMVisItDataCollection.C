#include "MFEMVisItDataCollection.h"

registerMooseObject("PlatypusApp", MFEMVisItDataCollection);

InputParameters
MFEMVisItDataCollection::validParams()
{
  NVTX3_FUNC_RANGE();
  InputParameters params = MFEMDataCollection::validParams();
  params.addClassDescription("Output for controlling MFEMVisItDataCollection inherited data.");
  params.addParam<unsigned int>("refinements",
                                0,
                                "Number of uniform refinements for oversampling "
                                "(refinement levels beyond any uniform "
                                "refinements)");
  return params;
}

MFEMVisItDataCollection::MFEMVisItDataCollection(const InputParameters & parameters)
  : MFEMDataCollection(parameters), _refinements(getParam<unsigned int>("refinements"))
{
  NVTX3_FUNC_RANGE();
}

std::shared_ptr<mfem::DataCollection>
MFEMVisItDataCollection::createDataCollection(const std::string & collection_name) const
{
  NVTX3_FUNC_RANGE();
  auto visit_dc = std::make_shared<mfem::VisItDataCollection>(_file_base.c_str() + collection_name);
  visit_dc->SetLevelsOfDetail(_refinements + 1);

  return visit_dc;
}
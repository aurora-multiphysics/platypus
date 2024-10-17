#pragma once
#include "MFEMDataCollection.h"
#include "mfem.hpp"

/**
 * Class for output information saved in MFEM ADIOS2DataCollections
 */
class MFEMADIOS2DataCollection : public MFEMDataCollection
{
public:
  static InputParameters validParams();

  MFEMADIOS2DataCollection(const InputParameters & parameters);

  virtual mfem::DataCollection & getDataCollection() override { return _adios_dc; }

protected:
  unsigned int _refinements;
  const std::string _engine_type;
  mfem::ADIOS2DataCollection _adios_dc;
};

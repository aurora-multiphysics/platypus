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

  std::shared_ptr<mfem::DataCollection>
  createDataCollection(const std::string & collection_name) const override;

protected:
  unsigned int _refinements;
  const std::string _engine_type;
};

#pragma once

// MOOSE includes
#include "MFEMTransfer.h"
#include "MultiMooseEnum.h"
#include "MultiApp.h"

// Forwards declaration.
class MFEMProblem;

/*
 * This class adds a getMFEMProblem method.
 */
class MFEMToMOOSECopyTransfer : public MFEMTransfer
{
public:
  static InputParameters validParams();

  MFEMToMOOSECopyTransfer(const InputParameters & parameters);

  void initialSetup() override;

  void execute() override;

  /**
   * Use this getter to obtain the MultiApp for transfers with a single direction
   */
  const std::shared_ptr<MultiApp> getMultiApp() const
  {
    if (_from_multi_app && _to_multi_app && _from_multi_app != _to_multi_app)
      mooseError("Unclear which app you want to retrieve from Transfer ", name());
    else if (_from_multi_app)
      return _from_multi_app;
    else if (_to_multi_app)
      return _to_multi_app;
    else
      mooseError("Should not get here, there should be a multiapp");
  }

  /// Get the MultiApp to transfer data from
  const std::shared_ptr<MultiApp> getFromMultiApp() const
  {
    if (!_from_multi_app)
      mooseError(
          "A from_multiapp was requested but is unavailable. Check the from_multi_app parameter");
    else
      return _from_multi_app;
  }

  /// Get the MultiApp to transfer data to
  const std::shared_ptr<MultiApp> getToMultiApp() const
  {
    if (!_to_multi_app)
      mooseError(
          "A to_multiapp was requested but is unavailable. Check the to_multi_app parameter");
    else
      return _to_multi_app;
  }

  /**
   * Get the name of thing being transferred from
   * @return the name of the multiapp or "Parent"
   */
  std::string getFromName() const
  {
    if (_from_multi_app)
      return _from_multi_app->name();
    else
      return "Parent";
  }

  /**
   * Get the name of thing being transferred to
   * @return the name of the multiapp or "Parent"
   */
  std::string getToName() const
  {
    if (_to_multi_app)
      return _to_multi_app->name();
    else
      return "Parent";
  }

protected:
  /**
   * Performs the transfer of a variable between two problems if they have the same mesh.
   */
  void transfer(FEProblemBase & to_problem, FEProblemBase & from_problem);

private:
  /// The MultiApps this Transfer is transferring data to or from
  std::shared_ptr<MultiApp> _from_multi_app;
  std::shared_ptr<MultiApp> _to_multi_app;
};
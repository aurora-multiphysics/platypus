// largely copying the mfemmeshtest here!

/*

  Could also try to inherit from the MFEMSolverTest here, but make
  sure to bypass its constructor
*/

#include "MFEMMesh.h"

#include "MFEMDiffusionKernel.h"
#include "MFEMScalarDirichletBC.h"
#include "MFEMHypreBoomerAMG.h"

#include "MMSTest.h"

class HypreBoomerAMGMMSTest : public MMSTestBase
{
public:
  virtual void RunConvergenceTest()
  {
    for (int i=0; i<3; i++)
    {
      // solver params!
      InputParameters solver_params          = _factory.getValidParams("MFEMHypreBoomerAMG");
      solver_params.set<double>("l_tol")     = _tol; // HypreBoomerAMG cannot set absolute tolerance
      solver_params.set<int>("l_max_its")    = _max_iters;
      MFEMHypreBoomerAMG & solver =
          addObject<MFEMHypreBoomerAMG>("MFEMHypreBoomerAMG", "solver" + i, solver_params);
      
      // Test MFEMSolver returns an solver of the expected type
      auto solver_downcast = std::dynamic_pointer_cast<mfem::HypreBoomerAMG>(solver.getSolver());
      solver_downcast->SetErrorMode(mfem::HypreSolver::ErrorMode::IGNORE_HYPRE_ERRORS);
      TestDiffusionSolve(*solver_downcast.get());
    }
  }
};

TEST_F( HypreBoomerAMGMMSTest, MMSTest )
{
  RunConvergenceTest();
  ASSERT_GE( EstimateConvergenceRate(), 0.0 );
}

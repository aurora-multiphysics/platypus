// largely copying the mfemmeshtest here!

/*

  Could also try to inherit from the MFEMSolverTest here, but make
  sure to bypass its constructor
*/

#include "MFEMMesh.h"

#include "MFEMDiffusionKernel.h"
#include "MFEMScalarDirichletBC.h"
#include "MFEMHypreBoomerAMG.h"
#include "MFEMGMRESSolver.h"
#include "MFEMHypreGMRES.h"
#include "MFEMHyprePCG.h"

#include "MMSTest.h"


class MFEMGMRESMMSTest : public MMSTestBase
{
public:
  virtual mfem::Solver & SetUpSolver() override
  {
    // solver params!
    InputParameters solver_params          = _factory.getValidParams("MFEMGMRESSolver");
    solver_params.set<double>("l_tol")     = 0.0;
    solver_params.set<double>("l_abs_tol") = _tol;
    solver_params.set<int>("l_max_its")    = _max_iters;
    MFEMGMRESSolver & solver =
      // We need a unique name every time we make a new solver object. So we bodge it by
      // appending the finite element order and the current refinement level to the name here

      addObject<MFEMGMRESSolver>("MFEMGMRESSolver", "solver" + std::to_string(_refinement_level) + std::to_string(_mesh_order), solver_params);
    
    // Test MFEMSolver returns an solver of the expected type
    auto solver_downcast = std::dynamic_pointer_cast<mfem::GMRESSolver>(solver.getSolver());

    // ASSERT_NE(solver_downcast.get(), nullptr);
    return *solver_downcast.get();
  }
};


class MFEMHypreGMRESMMSTest : public MMSTestBase
{
public:
  virtual mfem::Solver & SetUpSolver() override
  {
    InputParameters solver_params          = _factory.getValidParams("MFEMHypreGMRES");
    solver_params.set<double>("l_tol")     = 0.0;
    solver_params.set<double>("l_abs_tol") = _tol;
    solver_params.set<int>("l_max_its")    = _max_iters;

    // Construct kernel
    MFEMHypreGMRES & solver = 
      addObject<MFEMHypreGMRES>("MFEMHypreGMRES", "solver1" + std::to_string(_refinement_level) + std::to_string(_mesh_order), solver_params);

    // Test MFEMKernel returns an integrator of the expected type
    auto solver_downcast = std::dynamic_pointer_cast<mfem::HypreGMRES>(solver.getSolver());

    // ASSERT_NE(solver_downcast.get(), nullptr);
    return *solver_downcast.get();
  }
};


class MFEMHyprePCGMMSTest : public MMSTestBase
{
public:
  virtual mfem::Solver & SetUpSolver() override
  {
    InputParameters solver_params          = _factory.getValidParams("MFEMHyprePCG");
    solver_params.set<double>("l_tol")     = 0.0;
    solver_params.set<double>("l_abs_tol") = _tol;
    solver_params.set<int>("l_max_its")    = _max_iters;

    // Construct kernel
    MFEMHyprePCG & solver = 
      addObject<MFEMHyprePCG>("MFEMHyprePCG", "solver1" + std::to_string(_refinement_level) + std::to_string(_mesh_order), solver_params);

    // Test MFEMKernel returns an integrator of the expected type
    auto solver_downcast = std::dynamic_pointer_cast<mfem::HyprePCG>(solver.getSolver());

    // ASSERT_NE(solver_downcast.get(), nullptr);
    return *solver_downcast.get();
  }
};


TEST_F( MFEMGMRESMMSTest, MMSTest )
{
  RunConvergenceTest();

  auto iter = _log_log_gradients.begin();

  for ( int i=_lowest_mesh_order; i<=_max_mesh_order; i++, iter++ )
    ASSERT_NEAR( *iter, 1.0 +(double)i , 0.02 );
}

TEST_F( MFEMHypreGMRESMMSTest, MMSTest )
{
  RunConvergenceTest();

  auto iter = _log_log_gradients.begin();

  for ( int i=_lowest_mesh_order; i<=_max_mesh_order; i++, iter++ )
    ASSERT_NEAR( *iter, 1.0 +(double)i , 0.02 );
}


TEST_F( MFEMHyprePCGMMSTest, MMSTest )
{
  RunConvergenceTest();

  auto iter = _log_log_gradients.begin();

  for ( int i=_lowest_mesh_order; i<=_max_mesh_order; i++, iter++ )
    ASSERT_NEAR( *iter, 1.0 +(double)i , 0.02 );
}


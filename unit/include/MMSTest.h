#pragma once

#include "gtest/gtest.h"

#include "MFEMMesh.h"
#include "MFEMProblem.h"
#include "AppFactory.h"
#include "MooseMain.h"

/**
 *
 * Base class for all MMS problems. Contains member variables to be filled later
 * and some basic functions. 
 *
 *
 */
class MMSTestBase : public ::testing::Test
{
public:
  //! Provide arguments
  MMSTestBase(const std::string & app_name="PlatypusApp", const std::string& mesh_file_name = "data/beam-tet.mesh")
    : _mfem_mesh_ptr(nullptr), _app(Moose::createMooseApp(app_name, 0, nullptr)), _factory(_app->getFactory()),
      _mesh_file_name(mesh_file_name), _dimension(3), _num_refinements(5),
      _refinement_level(0)
  {}

  static constexpr double _tol               = 1e-16;
  static constexpr int    _max_iters         = 25;


  //! Pure virtual function to estimate convergence rate
  virtual double EstimateConvergenceRate() = 0;

  //! Pure virtual method to set up solver
  virtual mfem::Solver & SetUpSolver() = 0;

  virtual void RunConvergenceTest() = 0;

  virtual std::function<double(const mfem::Vector& x, mfem::real_t)> GetUExact() = 0;
  virtual std::function<double(const mfem::Vector& x, mfem::real_t)> GetFExact() = 0;

  void BuildObjects();

  template <typename T>
  T & addObject(const std::string & type, const std::string & name, InputParameters & params);

protected:
  std::unique_ptr<MFEMMesh>    _mfem_mesh_ptr;
  std::shared_ptr<MooseApp>    _app;
  Factory &                    _factory;
  std::shared_ptr<MFEMProblem> _mfem_problem;
  std::string                  _mesh_file_name;

  //! Stores average volume of a mesh element
  std::list<double>            _mesh_element_sizes;
  //! Stores l2 errors for each level of refinement
  std::list<double>            _l2_errors;
  //! Stores estimated gradients of the log-log plots, one for each finite element order, starting
  //! with the lowest order at the front()
  std::list<double>            _log_log_gradients;

  //! Dimensions of the problem; 3 by default  
  const int _dimension;

  //! parameter to describe how many times to refine the mesh for any given finite element
  //! order. Defaults to 5 in the constructor
  int _num_refinements;

  //! current refinement level
  int _refinement_level;


  //! Current finite element order ( 1 <= _fe_order <= _max_fe_order )
  int _fe_order;
};


//! Copied from the object unit test class
void
MMSTestBase::BuildObjects()
{
  InputParameters mesh_params           = _factory.getValidParams("MFEMMesh");
  mesh_params.set<MeshFileName>("file") = _mesh_file_name;
  _mfem_mesh_ptr = _factory.createUnique<MFEMMesh>("MFEMMesh", "moose_mesh" + std::to_string(_fe_order), mesh_params);
  _mfem_mesh_ptr->setMeshBase(_mfem_mesh_ptr->buildMeshBaseObject());
  _mfem_mesh_ptr->buildMesh();

  InputParameters problem_params                  = _factory.getValidParams("MFEMProblem");
  problem_params.set<MooseMesh *>("mesh")         = _mfem_mesh_ptr.get();
  problem_params.set<std::string>("_object_name") = "name2";
 
  _mfem_problem = _factory.create<MFEMProblem>("MFEMProblem", "problem" + std::to_string(_fe_order), problem_params);

  _app->actionWarehouse().problemBase() = _mfem_problem;
}


template <typename T>
T &
MMSTestBase::addObject(const std::string & type,
                              const std::string & name,
                              InputParameters & params)
{
  auto objects = _mfem_problem->addObject<T>(type, name, params);
  mooseAssert(objects.size() == 1, "Doesn't work with threading");
  return *objects[0];
}


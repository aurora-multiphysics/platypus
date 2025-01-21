#pragma once

#include "gtest/gtest.h"

#include "MFEMMesh.h"
#include "MFEMProblem.h"
#include "AppFactory.h"
#include "MooseMain.h"


/**
 *
 * Base class for performing MMS tests on diffusion-like problems
 *
 * Currently only does this in space.
 * 
 * For now - must be able to construct with no arguments!
 *
 */

class MMSTestBase : public ::testing::Test
{
public:
  //! Provide arguments
  MMSTestBase(const std::string & app_name="PlatypusApp", const std::string& mesh_file_name = "data/mug.e")
    : _app(Moose::createMooseApp(app_name, 0, nullptr)), _factory(_app->getFactory()),
      _mesh_file_name(mesh_file_name)
  {
    BuildObjects();
  }

  //! Function for uExact; override this in child class if needs be
  // virtual void 

  // maybe should be static
  double EstimateConvergenceRate();

  //! Runs one loop of the problem and refines the mesh at the end
  //! We don't want this to be overridden!
  virtual void TestDiffusionSolve(mfem::Solver & solver) final;

  //! Pure virtual setup method
  virtual void RunConvergenceTest() = 0;

  static constexpr double     _tol       = 1e-16;
  static constexpr int        _max_iters = 100;

  //! Overwrite these for different forcing/exact functions
  virtual std::function<double(const mfem::Vector& x)> GetUExact();
  virtual std::function<double(const mfem::Vector& x)> GetFExact();

protected:
  std::unique_ptr<MFEMMesh>    _mfem_mesh_ptr;
  std::shared_ptr<MooseApp>    _app;
  Factory &                    _factory;
  std::shared_ptr<MFEMProblem> _mfem_problem;
  std::string                  _mesh_file_name;

  std::list<double>   _mesh_element_sizes;
  std::list<double>   _l2_errors;

  void BuildObjects();

  template <typename T>
  T & addObject(const std::string & type, const std::string & name, InputParameters & params);

};

std::function<double(const mfem::Vector& x)>
MMSTestBase::GetUExact()
{
  return [](const mfem::Vector& x)
    {return sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );};
}

std::function<double(const mfem::Vector& x)>
MMSTestBase::GetFExact()
{
  return [](const mfem::Vector& x)
    {return 8 * M_PI * M_PI * sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );};
}


//! Copied from the object unit test class
void
MMSTestBase::BuildObjects()
{
  InputParameters mesh_params = _factory.getValidParams("MFEMMesh");
  mesh_params.set<MeshFileName>("file") = _mesh_file_name;
  _mfem_mesh_ptr = _factory.createUnique<MFEMMesh>("MFEMMesh", "moose_mesh", mesh_params);
  _mfem_mesh_ptr->setMeshBase(_mfem_mesh_ptr->buildMeshBaseObject());
  _mfem_mesh_ptr->buildMesh();

  InputParameters problem_params = _factory.getValidParams("MFEMProblem");
  problem_params.set<MooseMesh *>("mesh") = _mfem_mesh_ptr.get();
  problem_params.set<std::string>("_object_name") = "name2";
  _mfem_problem = _factory.create<MFEMProblem>("MFEMProblem", "problem", problem_params);

  _app->actionWarehouse().problemBase() = _mfem_problem;
}

void MMSTestBase::TestDiffusionSolve(mfem::Solver & solver)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int ne = 4;
  int dim      = 3;

  // create parmesh from regular mesh
  mfem::ParMesh pmesh( _mfem_mesh_ptr->getMFEMParMesh() );

  _mesh_element_sizes.push_front( pmesh.GetElementSize(0) );

  int order = 2;

  // this is "re" - defining the FE space
  mfem::H1_FECollection fec(order, dim);
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);

  // FROM HERE ONWARDS IS SETTING DIRICHLET BOUNDARY CONDITIONS
  mfem::Array<int> ess_tdof_list, ess_bdr;

  if (pmesh.bdr_attributes.Size())
  {
    ess_bdr.SetSize(pmesh.bdr_attributes.Max());
    ess_bdr = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  }

  // forcing term
  mfem::FunctionCoefficient f( GetFExact() );
  mfem::ParLinearForm       b( &fespace );
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(f));
  b.Assemble();

  mfem::ParBilinearForm a(&fespace);
  mfem::ConstantCoefficient one(1.0);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.Assemble();

  mfem::ParGridFunction x(&fespace);

  // create coeff called "u exact"
  mfem::FunctionCoefficient uex( GetUExact() );
  x = 0.0;

  // x is 0 everywhere, except on the boundary, where we have it matching the exact solution
  x.ProjectBdrCoefficient(uex, ess_bdr);

  mfem::OperatorPtr A;
  mfem::Vector B, X;
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);


  solver.SetOperator(*A);
  solver.Mult(B, X); // this is the actual solve!

  mfem::Vector Y(X.Size());
  A->Mult(X, Y);
  Y -= B; // check AX = B -> store AX in Y

  mfem::ParGridFunction uexact(&fespace);
  uexact.ProjectCoefficient( uex ); // every quadrature point / dof / node(!) will have the uexact value
  x = X;
  x.ComputeL2Error( uex );

  // store the error we calculated
  _l2_errors.push_front( x.ComputeL2Error( uex ) );


  // Finally - refine the mesh before the next time we call this function again!
  _mfem_mesh_ptr->getMFEMParMesh().UniformRefinement();

}

double
MMSTestBase::EstimateConvergenceRate()
{
  return
    _mesh_element_sizes.size() <= 1 or _l2_errors.size() <= 1 ?
    -1 :
    std::log( _l2_errors.back() / _l2_errors.front() ) / std::log( _mesh_element_sizes.back() / _mesh_element_sizes.front() );
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


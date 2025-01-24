#pragma once

#include "gtest/gtest.h"

#include "MFEMMesh.h"
#include "MFEMProblem.h"
#include "AppFactory.h"
#include "MooseMain.h"

// GetNE() on parmesh returns number of elements

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
  MMSTestBase(const std::string & app_name="PlatypusApp", const std::string& mesh_file_name = "data/simple-cube-multi-element-order2.e")
    : _app(Moose::createMooseApp(app_name, 0, nullptr)), _factory(_app->getFactory()),
      _mesh_file_name(mesh_file_name), _element_order(2), _dimension(3)
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

  //! Override these for different forcing/exact functions
  virtual std::function<double(const mfem::Vector& x)> GetUExact();
  virtual std::function<double(const mfem::Vector& x)> GetFExact();

  //! Calculate the average volume of a mesh element by calculating the
  //! entire mesh volume and dividing by the total number of elements
  virtual double CaculateElementVolume( mfem::ParMesh* ) final;

protected:
  std::unique_ptr<MFEMMesh>    _mfem_mesh_ptr;
  std::shared_ptr<MooseApp>    _app;
  Factory &                    _factory;
  std::shared_ptr<MFEMProblem> _mfem_problem;
  std::string                  _mesh_file_name;

  std::list<double>            _mesh_element_sizes;
  std::list<double>            _l2_errors;

  int _element_order;
  int _dimension;

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

  // create parmesh from regular mesh
  mfem::ParMesh pmesh( _mfem_mesh_ptr->getMFEMParMesh() );

  _mesh_element_sizes.push_front(
    CaculateElementVolume( &pmesh )
  );

  // define the FE space
  mfem::H1_FECollection       fec(_element_order, _dimension);
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);

  // Set the dirichlet boundary conditions
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

  // exact solution
  mfem::FunctionCoefficient uex( GetUExact() );
  x = 0.0;

  // x is 0 everywhere, except on the boundary, where we have it matching the exact solution
  x.ProjectBdrCoefficient(uex, ess_bdr);

  mfem::OperatorPtr A;
  mfem::Vector      B, X;
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

// Crudely estimate the slope of a log-log plot by working out the gradient of the line
// that joins the finest datapoint with the coarsest
// Finest goes at the front of the list since we push_front()
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

double
MMSTestBase::CaculateElementVolume( mfem::ParMesh* pmesh )
{
  // Set the space for the high-order mesh nodes.
  mfem::FiniteElementCollection *fec = new mfem::H1_FECollection(2, pmesh->Dimension());
  mfem::FiniteElementSpace nodal_fes(pmesh, fec, pmesh->SpaceDimension());

  mfem::FiniteElementSpace *fespace = new mfem::FiniteElementSpace(pmesh, fec);

  mfem::FunctionCoefficient func([](const mfem::Vector &x) { return 1.0; });

  mfem::GridFunction x(fespace);
  x.ProjectCoefficient(func);

  mfem::BilinearForm *vol = new mfem::BilinearForm(fespace);
  vol->AddDomainIntegrator(new mfem::MassIntegrator());
  vol->Assemble();

  mfem::Array<int> bdr_markers(1);
  bdr_markers[0] = 1;

  mfem::BilinearForm *area = new mfem::BilinearForm(fespace);
  area->AddBoundaryIntegrator(new mfem::MassIntegrator(), bdr_markers);
  area->Assemble();

  double averageMeshElementSize = vol->InnerProduct(x, x) / pmesh->GetNE();

  delete fec;
  delete fespace;
  delete vol;
  delete area;

  // return total vol divided by number of elements
  return averageMeshElementSize;
}

#pragma once

#include "MMSTest.h"


/**
 *
 * Class for performing MMS tests on diffusion-like problems
 * 
 * Following the fenics book, we expect that if we run the refined
 * mesh test for various fe element orders, we should see that the
 * slope of the log-log plot in each case is proportional to 
 * the mesh element order
 * 
 * To produce a new test case, all that needs to be done is overwrite
 * the SetUpSolver() method. This should set up the desired solver using the addObject
 * method and then return a reference to it.
 * 
 * 
 * Everything is handled within the RunConvergenceTest() method. It has two main loops.
 * The first loops over each finite element order that we want to test. The inner loop
 * then builds the mesh from the given mesh file, and refines the mesh to the second-finest
 * level. From there, we pass over to TestDiffusionSolve(), which we call with the pure
 * virtual SetUpSolver() method. 
 * 
 * 
 * The actual test is in the method TestDiffusionSolve(). This defines the finite element space from
 * the mesh we set up during BuildObjects, sets up the linear system, and then solves using our
 * exact functions for f and u. Finally, it compares the L2 error of the final estimation with the 
 * exact solution.
 * 
 */
class MMSDiffusionTestBase : public MMSTestBase
{
public:

  MMSDiffusionTestBase() : MMSTestBase(), _max_fe_order(3) {}

  static constexpr int    _lowest_fe_order = 1;

  //! Calculate gradient of log-log plot
  virtual double EstimateConvergenceRate() override;

  //! Runs one loop of the problem and refines the mesh at the end
  void TestDiffusionSolve(mfem::Solver & solver);

  //! Virtual setup method - we don't want this to be overridden!
  virtual void RunConvergenceTest() override;

  //! Override these for different forcing/exact functions
  virtual std::function<double(const mfem::Vector& x, mfem::real_t)> GetUExact() override;
  virtual std::function<double(const mfem::Vector& x, mfem::real_t)> GetFExact() override;

  //! Calculate the average volume of a mesh element by calculating the
  //! entire mesh volume and dividing by the total number of elements
  virtual double CaculateElementVolume( mfem::ParMesh* ) final;

  //! How many finite element orders to test? We start the finite element order at 1
  //! and increase it this many times...
  int _max_fe_order;


};

//! Standard U exact function for diffusion-like problem that has sinusoidal spatial dependence
//! and no temporal dependence.
std::function<double(const mfem::Vector& x, mfem::real_t t)>
MMSDiffusionTestBase::GetUExact()
{
  return [](const mfem::Vector& x, mfem::real_t t)
    {return sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );};
}

//! Standard forcing function based on above definition of U
std::function<double(const mfem::Vector& x, mfem::real_t t)>
MMSDiffusionTestBase::GetFExact()
{
  return [](const mfem::Vector& x, mfem::real_t t)
    {return 12 * M_PI * M_PI * sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );};
}


void MMSDiffusionTestBase::TestDiffusionSolve(mfem::Solver & solver)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create parmesh from regular mesh
  mfem::ParMesh pmesh( _mfem_mesh_ptr->getMFEMParMesh() );

  _mesh_element_sizes.push_front(
    std::pow( CaculateElementVolume( &pmesh ), 1./_dimension )
  );

  // define the FE space
  mfem::H1_FECollection       fec(_fe_order, _dimension);
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

  x.SetFromTrueDofs(X);
  x.ComputeL2Error( uex );

  // store the error we calculated. push to the front so that we have errors
  // for the finest mesh level at the front.
  _l2_errors.push_front( x.ComputeL2Error( uex ) );

}

//! Crudely estimate the slope of a log-log plot by working out the gradient of the line
//! that joins the finest datapoint with the coarsest
//! Finest goes at the front of the list since we push_front()
double
MMSDiffusionTestBase::EstimateConvergenceRate()
{
  if ( _mesh_element_sizes.size() <= 1 or _l2_errors.size() <= 1 )
    return -1;

  // get the errors for the finest mesh and second-finest
  double l2_error_0 = _l2_errors.front(); _l2_errors.pop_front();
  double l2_error_1 = _l2_errors.front(); _l2_errors.pop_front();

  // ditto for volumes
  double h_0 = _mesh_element_sizes.front(); _mesh_element_sizes.pop_front();
  double h_1 = _mesh_element_sizes.front(); _mesh_element_sizes.pop_front();

  double output = std::log( l2_error_1 / l2_error_0 ) / std::log( h_1 / h_0 );

  // clear these lists; we will reuse them next time we pick a new polynomial order
  _mesh_element_sizes.clear();
  _l2_errors.clear();

  return output;
}

double
MMSDiffusionTestBase::CaculateElementVolume( mfem::ParMesh* pmesh )
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

void
MMSDiffusionTestBase::RunConvergenceTest()
{

  // for mesh element order
  for (int fe_order=_lowest_fe_order; fe_order <= _max_fe_order; fe_order++)
  {
    // set the element order member variable
    _fe_order = fe_order;

    BuildObjects();

    // we only want to test the finest mesh, and the second finest.
    // So let's start by refining a few times

    // refines 3 times
    for (int r=0; r<_num_refinements-2; r++) _mfem_mesh_ptr->getMFEMParMesh().UniformRefinement();
    
    // another loop; each time we need to further refine the mesh
    for (int r=_num_refinements-2; r<_num_refinements; r++)
    {
      // set the member variable so that the current solver object gets a unique name
      _refinement_level = r;

      TestDiffusionSolve(
        SetUpSolver()
      );

      // Finally - refine the mesh before the next time we call this function again!
      _mfem_mesh_ptr->getMFEMParMesh().UniformRefinement();
    }

    _log_log_gradients.push_back( EstimateConvergenceRate() );
  }

}


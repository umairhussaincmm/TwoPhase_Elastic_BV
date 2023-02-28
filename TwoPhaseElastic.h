//
// Created by umair on 23/2/23.
//

#ifndef TWO_PHASE_ELASTICITY_TWOPHASEELASTIC_H
#define TWO_PHASE_ELASTICITY_TWOPHASEELASTIC_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_in.h>

//For Parallel Computation
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/process_grid.h>

#include <deal.II/distributed/solution_transfer.h>

//For Adaptive Mesh Refinement
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <fstream>
#include <iostream>

using namespace dealii;

//Structure to store stress-strain point data
struct PointHistory
{
    std::vector<std::vector<double>> stress; // No. of GP * No of tensor components
    std::vector<std::vector<double>> strain;
};

class TwoPhaseElastic {
public:
    TwoPhaseElastic();
    void run();
    const double M, kappa, c_max, R, T, omega; //Phase field model parameters
    const double lambda, G, beta, ref_conc; //Elastic parameters
    double err;
    double Farad, lam, A, B, cm, cL, kre, n, Vmax, Vmin; //B-V parameters
private:
    void          make_grid_and_dofs(unsigned int timestep_number, unsigned int it);
    void          assemble_system();
    void          solve();
    void          output_results(const unsigned int timestep_number) const;
    void          applying_bc();
    // Stress-strain data
    void    setup_quadrature_point_history();
    void    stress_strain_calc();
    std::vector<double>        fchem(double c_hat);
    std::vector<long double>   BVflux(double c);

    MPI_Comm mpi_communicator;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;
    ConditionalOStream pcout;

    parallel::distributed::Triangulation<2> triangulation; //For AMR
    FESystem<2>          fe;
    DoFHandler<2>        dof_handler;
    GridIn<2>            gridin;
    AffineConstraints<double> hanging_node_constraints;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    std::vector<PointHistory> quadrature_point_history;

    PETScWrappers::MPI::SparseMatrix jacobian_matrix;

    double       time;
    const double final_time, time_step;
    const double theta;

    PETScWrappers::MPI::Vector conv_solution, solution_update, old_solution, solution_update_old, difference;
    PETScWrappers::MPI::Vector system_rhs, residual;
    Vector<double> conv_solution_np, old_solution_np; //To store data in non parallel vector
};

class InitialValues : public Function<2>
{
public:
    InitialValues(): Function<2>(4)
    {}
    virtual void vector_value(const Point<2> &p,
                              Vector<double> & value) const override;
};


#endif //TWO_PHASE_ELASTICITY_TWOPHASEELASTIC_H

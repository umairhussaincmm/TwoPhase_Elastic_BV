//
// Created by umair on 24/2/23.
//
#include "TwoPhaseElastic.h"

void TwoPhaseElastic::solve(){
    SolverControl cn;
    PETScWrappers::SparseDirectMUMPS A_direct(cn, mpi_communicator);
    A_direct.solve(jacobian_matrix, solution_update, system_rhs);

    //Updating the solution
    PETScWrappers::MPI::Vector completely_distributed_solution(locally_owned_dofs,mpi_communicator);
    completely_distributed_solution = conv_solution;
    completely_distributed_solution.add(1, solution_update);
    completely_distributed_solution.compress(VectorOperation::add);
    conv_solution = completely_distributed_solution;
}

void TwoPhaseElastic::output_results(const unsigned int timestep_number) const {
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names;
    solution_names.emplace_back ("u1");
    solution_names.emplace_back ("u2");
    solution_names.emplace_back ("mu");
    solution_names.emplace_back ("c");

    data_out.add_data_vector(old_solution, solution_names);
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    // Adding stress-strain
    Vector<double> stress_x(triangulation.n_active_cells());
    Vector<double> stress_y(triangulation.n_active_cells());
    Vector<double> stress_z(triangulation.n_active_cells());
    Vector<double> stress_xy(triangulation.n_active_cells());
    Vector<double> stress_von_mises(triangulation.n_active_cells());
    Vector<double> strain_x(triangulation.n_active_cells());
    Vector<double> strain_y(triangulation.n_active_cells());
    Vector<double> strain_xy(triangulation.n_active_cells());
    QGauss<2> quadrature_formula(fe.degree + 1);
    for (auto &cell : triangulation.active_cell_iterators()){
        double accumulated_stress_x=0, accumulated_stress_y=0, accumulated_stress_z=0, accumulated_stress_xy=0, accumulated_strain_x=0, accumulated_strain_y=0, accumulated_strain_xy=0;
        for (unsigned int q = 0; q < quadrature_formula.size(); ++q) {
            accumulated_stress_x += quadrature_point_history[cell->active_cell_index()].stress[q][0];
            accumulated_stress_y += quadrature_point_history[cell->active_cell_index()].stress[q][1];
            accumulated_stress_z += quadrature_point_history[cell->active_cell_index()].stress[q][2];
            accumulated_stress_xy += quadrature_point_history[cell->active_cell_index()].stress[q][3];
            accumulated_strain_x += quadrature_point_history[cell->active_cell_index()].strain[q][0];
            accumulated_strain_y += quadrature_point_history[cell->active_cell_index()].strain[q][1];
            accumulated_strain_xy += quadrature_point_history[cell->active_cell_index()].strain[q][3];
        }
        stress_x(cell->active_cell_index()) = (accumulated_stress_x / quadrature_formula.size());
        stress_y(cell->active_cell_index()) = (accumulated_stress_y / quadrature_formula.size());
        stress_z(cell->active_cell_index()) = (accumulated_stress_z / quadrature_formula.size());
        stress_xy(cell->active_cell_index()) = (accumulated_stress_xy / quadrature_formula.size());
        strain_x(cell->active_cell_index()) = (accumulated_strain_x / quadrature_formula.size());
        strain_y(cell->active_cell_index()) = (accumulated_strain_y / quadrature_formula.size());
        strain_xy(cell->active_cell_index()) = (accumulated_strain_xy / quadrature_formula.size());
    }
    data_out.add_data_vector(stress_x, "Stress_x");
    data_out.add_data_vector(stress_y, "Stress_y");
    data_out.add_data_vector(stress_z, "Stress_z");
    data_out.add_data_vector(stress_xy, "Stress_xy");
    data_out.add_data_vector(strain_x, "Total_strain_x");
    data_out.add_data_vector(strain_y, "Total_strain_y");
    data_out.add_data_vector(strain_xy, "Total_strain_xy");

    data_out.build_patches();
    //data_out.write_vtk(output);
    data_out.write_vtu_with_pvtu_record("results/trial7-nomuflux/", "solution",timestep_number,mpi_communicator,2,16);
}
//
// Created by umair on 23/2/23.
//
#include "TwoPhaseElastic.h"
void TwoPhaseElastic::make_grid_and_dofs(unsigned int timestep_number, unsigned int it) {
    if (timestep_number ==1 && it ==1){
        //Reading mesh
        gridin.attach_triangulation(triangulation);
        std::ifstream f("mesh/lithiation3.msh");
        gridin.read_msh(f);
        std::ofstream out ("grid.eps");
        GridOut grid_out;
        grid_out.write_eps (triangulation, out);
        triangulation.reset_manifold(0);
        //triangulation.refine_global(1);
        //GridTools::partition_triangulation(n_mpi_processes, triangulation);
        dof_handler.distribute_dofs(fe);
        //DoFRenumbering::subdomain_wise(dof_handler);
        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
        old_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        conv_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        system_rhs.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        solution_update.reinit(locally_owned_dofs, mpi_communicator);
        hanging_node_constraints.clear();
        hanging_node_constraints.reinit(locally_relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        /*VectorTools::interpolate_boundary_values(dof_handler,
                                                 8,
                                                 ConstantFunction<2>(40000, 6),
                                                 hanging_node_constraints);*/
        hanging_node_constraints.close();
        DynamicSparsityPattern dsp(locally_relevant_dofs);
        DoFTools::make_sparsity_pattern(dof_handler, dsp, hanging_node_constraints, false);

        SparsityTools::distribute_sparsity_pattern(dsp,
                                                   dof_handler.locally_owned_dofs(),
                                                   mpi_communicator,
                                                   locally_relevant_dofs);
        jacobian_matrix.reinit(locally_owned_dofs,
                               locally_owned_dofs,
                               dsp,
                               mpi_communicator);
        conv_solution_np.reinit(dof_handler.n_dofs());
        old_solution_np.reinit(dof_handler.n_dofs());
    }
    else{
        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
        system_rhs.clear();
        solution_update.clear();
        jacobian_matrix.clear();
        pcout<<"Reinitializing system rhs"<<std::endl;
        system_rhs.reinit(locally_owned_dofs, mpi_communicator);
        solution_update.reinit(locally_owned_dofs, mpi_communicator);
        DynamicSparsityPattern dsp(locally_relevant_dofs);
        DoFTools::make_sparsity_pattern(dof_handler, dsp, hanging_node_constraints, false);
        SparsityTools::distribute_sparsity_pattern(dsp,
                                                   dof_handler.locally_owned_dofs(),
                                                   mpi_communicator,
                                                   locally_relevant_dofs);

        jacobian_matrix.reinit(locally_owned_dofs,
                               locally_owned_dofs,
                               dsp,
                               mpi_communicator);
    }
}
//
// Created by umair on 24/2/23.
//
#include "TwoPhaseElastic.h"

void TwoPhaseElastic::run(){
    unsigned int timestep_number = 1;
    make_grid_and_dofs(timestep_number, 1);
    setup_quadrature_point_history();
    //Initialise the solution
    InitialValues initial_value;
    PETScWrappers::MPI::Vector local_vector(locally_owned_dofs,
                                            mpi_communicator);
    VectorTools::interpolate(dof_handler,
                             initial_value,
                             local_vector);
    old_solution = local_vector;
    conv_solution = old_solution;
    //pcout << "No error till here" << std::endl;
    //Applying Boundary Conditions at t=0
    applying_bc();
    //Plotting initial solution
    output_results(0);


    //Starting clock right before time loop
    auto start = std::chrono::steady_clock::now();

    //Time steps begin here:
    for (; time <= final_time; time += time_step, ++timestep_number) {
        pcout << "Time step " << timestep_number << " at t=" << time+time_step
              << std::endl;
        conv_solution.operator=(old_solution); // initialising the newton solution
        //Newton's iterations begin here:
        for (unsigned int it = 1; it <= 100; ++it) {
            pcout << "Newton iteration number:" << it << std::endl;
            if (it == 100) {
                pcout << "Convergence Failure!!!!!!!!!!!!!!!" << std::endl;
                std::exit(0);
            }
            //Saving parallel vectors as non-parallel ones
            conv_solution_np = conv_solution;
            old_solution_np = old_solution;
            //Initialise the delta solution as zero
            VectorTools::interpolate(dof_handler,
                                     ZeroFunction<2>(4),
                                     solution_update);
            solution_update.compress(VectorOperation::insert);
            assemble_system();
            solve();
            //Checking for convergence
            double residual_norm = system_rhs.l2_norm();
            //pcout << "Nothing wrong till here!!!!!!" << std::endl;
            pcout << "the residual is:" << residual_norm << std::endl;
            if (residual_norm <= (1e-4)) {
                pcout << "Solution Converged!" << std::endl;
                break;
            }

        }
        old_solution.operator=(conv_solution);
        old_solution.compress(VectorOperation::add);
        //Plot the solution at the end of Newton iterations
        if (timestep_number%10 == 0){
            // Refining
            //adaptively_refine(timestep_number+1,1);
            stress_strain_calc();
            output_results(timestep_number);
        }
    }

    //Stopping clock right at end of time loop
    auto end = std::chrono::steady_clock::now();

    // Store the time difference between start and end
    auto diff = end - start;
    pcout << std::chrono::duration <double, std::kilo> (diff).count() << " x10^3 s" << std::endl;
}
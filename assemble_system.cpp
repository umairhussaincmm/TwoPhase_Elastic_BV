//
// Created by umair on 23/2/23.
//
#include "TwoPhaseElastic.h"

void TwoPhaseElastic::assemble_system() {
    FEValuesExtractors::Scalar u1(0);
    FEValuesExtractors::Scalar u2(1);
    FEValuesExtractors::Scalar mu(2);
    FEValuesExtractors::Scalar c(3);

    QGauss<2> quadrature_formula(fe.degree + 1);//Increase it for more accuracy
    FEValues<2> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);
    QGauss<1> face_quadrature_formula(fe.degree + 1);
    FEFaceValues<2> fe_face_values(fe, face_quadrature_formula, update_values | update_quadrature_points |
                                                                update_normal_vectors |  update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    //To copy values and gradients from previous iteration
    //Old Newton iteration
    std::vector<Tensor<1, 2>> old_newton_solution_gradients_u1(n_q_points);
    std::vector<double> old_newton_solution_values_u1(n_q_points);
    std::vector<Tensor<1, 2>> old_newton_solution_gradients_u2(n_q_points);
    std::vector<double> old_newton_solution_values_u2(n_q_points);
    std::vector<Tensor<1, 2>> old_newton_solution_gradients_c(n_q_points);
    std::vector<double> old_newton_solution_values_c(n_q_points);
    std::vector<Tensor<1, 2>> old_newton_solution_gradients_mu(n_q_points);
    std::vector<double> old_newton_solution_values_mu(n_q_points);

    //Old time step iteration
    std::vector<Tensor<1, 2>> old_time_solution_gradients_u1(n_q_points);
    std::vector<double> old_time_solution_values_u1(n_q_points);
    std::vector<Tensor<1, 2>> old_time_solution_gradients_u2(n_q_points);
    std::vector<double> old_time_solution_values_u2(n_q_points);
    std::vector<Tensor<1, 2>> old_time_solution_gradients_c(n_q_points);
    std::vector<double> old_time_solution_values_c(n_q_points);
    std::vector<Tensor<1, 2>> old_time_solution_gradients_mu(n_q_points);
    std::vector<double> old_time_solution_values_mu(n_q_points);

    // Face values
    std::vector<Tensor<1, 2>> on_grad_face_c(n_q_points);
    std::vector<double> on_value_face_c(n_q_points);
    std::vector<Tensor<1, 2>> ot_grad_face_c(n_q_points);
    std::vector<double> ot_value_face_c(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    jacobian_matrix.operator=(0.);
    system_rhs.operator=(0.);
    //Creating temporary Non-Ghosted Vector
    PETScWrappers::MPI::Vector local_vector_rhs(locally_owned_dofs, mpi_communicator);
    for (const auto &cell : dof_handler.active_cell_iterators()) {
        //pcout << "Nothing wrong till here!!!!!!" << std::endl;
        if (cell->subdomain_id() == this_mpi_process) {
            cell_matrix = 0;
            cell_rhs = 0;
            fe_values.reinit(cell);

            fe_values[u1].get_function_values(conv_solution_np, old_newton_solution_values_u1);
            fe_values[u1].get_function_gradients(conv_solution_np, old_newton_solution_gradients_u1);
            fe_values[u2].get_function_values(conv_solution_np, old_newton_solution_values_u2);
            fe_values[u2].get_function_gradients(conv_solution_np, old_newton_solution_gradients_u2);
            fe_values[c].get_function_values(conv_solution_np, old_newton_solution_values_c);
            fe_values[c].get_function_gradients(conv_solution_np, old_newton_solution_gradients_c);
            fe_values[mu].get_function_values(conv_solution_np, old_newton_solution_values_mu);
            fe_values[mu].get_function_gradients(conv_solution_np, old_newton_solution_gradients_mu);

            fe_values[u1].get_function_values(old_solution_np, old_time_solution_values_u1);
            fe_values[u1].get_function_gradients(old_solution_np, old_time_solution_gradients_u1);
            fe_values[u2].get_function_values(old_solution_np, old_time_solution_values_u2);
            fe_values[u2].get_function_gradients(old_solution_np, old_time_solution_gradients_u2);
            fe_values[c].get_function_values(old_solution_np, old_time_solution_values_c);
            fe_values[c].get_function_gradients(old_solution_np, old_time_solution_gradients_c);
            fe_values[mu].get_function_values(old_solution_np, old_time_solution_values_mu);
            fe_values[mu].get_function_gradients(old_solution_np, old_time_solution_gradients_mu);

            for (unsigned int q = 0; q < n_q_points; ++q) {
                double u1_on = old_newton_solution_values_u1[q];
                double u1_x_on = old_newton_solution_gradients_u1[q][0];
                double u1_y_on = old_newton_solution_gradients_u1[q][1];
                double u2_on = old_newton_solution_values_u2[q];
                double u2_x_on = old_newton_solution_gradients_u2[q][0];
                double u2_y_on = old_newton_solution_gradients_u2[q][1];
                double c_on = old_newton_solution_values_c[q];
                Tensor<1, 2> grad_c_on = old_newton_solution_gradients_c[q];
                double mu_on = old_newton_solution_values_mu[q];
                Tensor<1, 2> grad_mu_on = old_newton_solution_gradients_mu[q];

                double u1_ot = old_time_solution_values_u1[q];
                double u1_x_ot = old_time_solution_gradients_u1[q][0];
                double u1_y_ot = old_time_solution_gradients_u1[q][1];
                double u2_ot = old_time_solution_values_u2[q];
                double u2_x_ot = old_time_solution_gradients_u2[q][0];
                double u2_y_ot = old_time_solution_gradients_u2[q][1];
                double c_ot = old_time_solution_values_c[q];
                Tensor<1, 2> grad_c_ot = old_time_solution_gradients_c[q];
                double mu_ot = old_time_solution_values_mu[q];
                Tensor<1, 2> grad_mu_ot = old_time_solution_gradients_mu[q];
                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    double eta1_i = fe_values[u1].value(i, q);
                    double eta1_x_i = fe_values[u1].gradient(i, q)[0];
                    double eta1_y_i = fe_values[u1].gradient(i, q)[1];
                    double eta2_i = fe_values[u2].value(i, q);
                    double eta2_x_i = fe_values[u2].gradient(i, q)[0];
                    double eta2_y_i = fe_values[u2].gradient(i, q)[1];
                    double psi_i = fe_values[mu].value(i, q);
                    Tensor<1, 2> grad_psi_i = fe_values[mu].gradient(i, q);
                    double phi_i = fe_values[c].value(i, q);
                    Tensor<1, 2> grad_phi_i = fe_values[c].gradient(i, q);
                    for (unsigned int j = 0; j < dofs_per_cell; ++j){
                        double eta1_j = fe_values[u1].value(j, q);
                        double eta1_x_j = fe_values[u1].gradient(j, q)[0];
                        double eta1_y_j = fe_values[u1].gradient(j, q)[1];
                        double eta2_j = fe_values[u2].value(j, q);
                        double eta2_x_j = fe_values[u2].gradient(j, q)[0];
                        double eta2_y_j = fe_values[u2].gradient(j, q)[1];
                        double psi_j = fe_values[mu].value(j, q);
                        Tensor<1, 2> grad_psi_j = fe_values[mu].gradient(j, q);
                        double phi_j = fe_values[c].value(j, q);
                        Tensor<1, 2> grad_phi_j = fe_values[c].gradient(j, q);

                        double f1_u1 = 0;
                        double f1_u2 = 0;
                        double L_mu = grad_psi_i*(M*c_on*(1-c_on))*grad_psi_j; //do L by do mu
                        double f1_mu = time_step*theta*L_mu;
                        double L_c = grad_psi_j*(M*(1-2*c_on)*grad_mu_on)*phi_j; //do L by do c
                        double H = psi_i*(phi_j);
                        double f1_c = time_step*theta*L_c + H;

                        double muel_u1 = -beta*(3*lambda+2*G)*eta1_x_j;
                        double F_u1 = phi_i*muel_u1;
                        double f2_u1 = F_u1;
                        double muel_u2 = -beta*(3*lambda+2*G)*eta2_y_j;
                        double F_u2 = phi_i*muel_u2;
                        double f2_u2 = F_u2;
                        double N = phi_i*psi_j;
                        double f2_mu = -N;
                        double K = grad_phi_i*(kappa*grad_phi_j);
                        double muel_c = 3*beta*beta*(3*lambda+2*G)*phi_j;
                        double F_c = phi_i*(fchem(c_on)[2]*phi_j) + muel_c;
                        double f2_c = K + F_c;

                        double f3_u1 = eta1_x_i*(lambda+2*G)*eta1_x_j + eta1_y_i*G*eta1_y_j;
                        double f3_u2 = eta1_x_i*lambda*eta2_y_j + eta1_y_i*G*eta2_x_j;
                        double f3_mu = 0;
                        double f3_c = eta1_x_i*(-beta*(3*lambda+2*G))*phi_j;

                        double f4_u1 = eta2_y_i*lambda*eta1_x_j + eta2_x_i*G*eta1_y_j;
                        double f4_u2 = eta2_y_i*(lambda+2*G)*eta2_y_j + eta2_x_i*G*eta2_x_j;
                        double f4_mu = 0;
                        double f4_c = eta2_y_i*(-beta*(3*lambda+2*G))*phi_j;
                        //Assembling Jacobian matrix
                        cell_matrix(i,j) += (f1_u1+f1_u2+f1_c+f1_mu + f2_u1+f2_u2+f2_c+f2_mu + f3_u1+f3_u2+f3_c+f3_mu + f4_u1+f4_u2+f4_c+f4_mu)*fe_values.JxW(q);
                    }
                    //Finding f1, f2, f3 and f4 at previous newton iteration for rhs vector
                    double H_k = psi_i*c_on; //Mass matrix at previous newton iteration
                    double L_k = grad_psi_i*(M*c_on*(1-c_on)*grad_mu_on);
                    double H_t = psi_i*c_ot; // Mass matrix at previous time step
                    double L_t = grad_psi_i*(M*c_ot*(1-c_ot)*grad_mu_ot);
                    double f1n = H_k + theta*time_step*L_k - H_t + (1-theta)*time_step*L_t;

                    double N_k = phi_i*(mu_on); //N matrix at previous newton iteration
                    double K_k = grad_phi_i*(kappa*grad_c_on); // K matrix at previous newton iteration
                    double muel = -beta*(3*lambda+2*G)*(u1_x_on + u2_y_on - 3*beta*(c_on-ref_conc));
                    double F_k = phi_i* (fchem(c_on)[1] + muel);
                    double f2n = (K_k) - N_k + F_k;

                    double f3t1 = eta1_x_i*((lambda+2*G)*u1_x_on + lambda*u2_y_on - beta*(3*lambda+2*G)*(c_on-ref_conc));
                    double f3t2 = eta1_y_i*G*(u1_y_on+u2_x_on);
                    double f3n = f3t1 + f3t2;

                    double f4t1 = eta2_y_i*((lambda+2*G)*u2_y_on + lambda*u1_x_on - beta*(3*lambda+2*G)*(c_on-ref_conc));
                    double f4t2 = eta2_x_i*G*(u1_y_on+u2_x_on);
                    double f4n = f4t1 + f4t2;
                    //Assembling RHS vector
                    cell_rhs(i) -= (f1n + f2n + f3n + f4n)*fe_values.JxW(q);
                }
            }
            for (unsigned int face_number = 0; face_number < GeometryInfo<2>::faces_per_cell; ++face_number){
                if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 3)){
                    fe_face_values.reinit(cell, face_number);
                    for (unsigned int q_pt = 0; q_pt < n_face_q_points; ++q_pt){
                        double nbcmu = .0;
                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            double Qmu = nbcmu*fe_face_values[mu].value(i, q_pt);
                            double fmu_bv = - Qmu;
                            cell_rhs(i) -= (fmu_bv) * fe_values.JxW(q_pt);
                        }
                    }
                }
            }
            cell->get_dof_indices(local_dof_indices);
            hanging_node_constraints.distribute_local_to_global(cell_matrix,
                                                                cell_rhs,
                                                                local_dof_indices,
                                                                jacobian_matrix,
                                                                local_vector_rhs);
        }
    }
    jacobian_matrix.compress(VectorOperation::add);
    local_vector_rhs.compress(VectorOperation::add);
    //Applying zero BC
    ComponentMask u1_mask = fe.component_mask(u1);
    ComponentMask u2_mask = fe.component_mask(u2);
    ComponentMask mu_mask = fe.component_mask(mu);
    ComponentMask c_mask = fe.component_mask(c);
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              3,
                                              ZeroFunction<2>(4),
                                              boundary_values,c_mask);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              4,
                                              ZeroFunction<2>(4),
                                              boundary_values,mu_mask);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              4,
                                              ZeroFunction<2>(4),
                                              boundary_values,u1_mask);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              4,
                                              ZeroFunction<2>(4),
                                              boundary_values,u2_mask);
    MatrixTools::apply_boundary_values(boundary_values,
                                       jacobian_matrix,
                                       solution_update,
                                       local_vector_rhs, false);
    jacobian_matrix.compress(VectorOperation::add);
    local_vector_rhs.compress(VectorOperation::add);
    system_rhs = local_vector_rhs;
}

std::vector<double> TwoPhaseElastic::fchem(double c_hat){
    std::vector<double> f(3);
    f[0] = c_max*R*T*(c_hat*std::log(c_hat+err)+(1-c_hat)*std::log(1-c_hat+err) + omega*c_hat*(1-c_hat)); //value of fchem
    f[1] = c_max*R*T*(std::log((c_hat+err)/(1-c_hat+err)) + omega*(1-2*c_hat)); //first derivative
    f[2] = c_max*R*T*(1/((c_hat+err)*(1-c_hat+err)) - 2*omega); //second derivative
    //pcout << "Value of f_c: " << f[1] << std::endl;
    return f;
}

std::vector<long double> TwoPhaseElastic::BVflux(double c){
    std::vector<long double> J(2);
    double i0 = n*std::pow(cm-c,1-lam)*std::pow(c,lam);
    double i0_c = n*std::pow(c,lam-1)*(lam*cm-c)/(std::pow(cm-c,lam));
    double V = ((Vmax-Vmin)/2)*std::sin(2*3.14*time) + (Vmax+Vmin)/2;
    //std::cout << "Current value of V: " << V << std::endl;
    double U = 4.198 + 0.056*std::tanh(-14.55*(c)+8.6) - 0.027*(1/std::pow(.998-(c),.49)-1.9) - .157*std::exp(-.047*std::pow((c),8)) + .81*std::exp(-40*((c)-.134));
    double U_c = -0.815*std::pow((1/std::cosh(14.55*(c)-8.6)),2) - .013*(std::pow(.998-(c),-1.49)) + .059*std::pow((c),7)*std::exp(-.047*std::pow((c),8)) -32.4*std::exp(-40*((c)-.134));
    //std::cout << "Value of U: " << U << std::endl;
    double eta = V - U;
    double eta_c = -U_c;
    J[0] = (i0/Farad)*(std::exp(A*eta) - std::exp(B*eta)); //value of flux
    J[1] = i0_c*(std::exp(A*eta) - std::exp(B*eta))/Farad + eta_c*(i0/Farad)*(A*std::exp(A*eta) - B*std::exp(B*eta)); //first derivative
    //pcout << "Value of J_c: " << J[1] << std::endl;
    return J;
}

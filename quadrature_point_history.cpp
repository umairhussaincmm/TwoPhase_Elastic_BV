//
// Created by umair on 7/12/22.
//
#include "TwoPhaseElastic.h"

void TwoPhaseElastic::setup_quadrature_point_history(){
    pcout<<"Number of elements is: "<< triangulation.n_active_cells() <<std::endl;
    QGauss<2> quadrature_formula(fe.degree + 1);
    triangulation.clear_user_data();
    {
        std::vector<PointHistory> tmp;
        quadrature_point_history.swap(tmp);
    }
    quadrature_point_history.resize(triangulation.n_active_cells()); //Initialised qph object with size of number of local elements (Changed to total elements due to error with parallel code)
    for (auto &cell : dof_handler.active_cell_iterators()){
        quadrature_point_history[cell->active_cell_index()].stress.resize(quadrature_formula.size(),std::vector<double>(4));
        quadrature_point_history[cell->active_cell_index()].strain.resize(quadrature_formula.size(),std::vector<double>(4));
    }
}

std::vector<double> strain(std::vector<Tensor<1, 2>> disp_grad_u1, std::vector<Tensor<1, 2>> disp_grad_u2, int q_point){
    std::vector<double> strain_ele(4);
    strain_ele[0] = disp_grad_u1[q_point][0];
    strain_ele[1] = disp_grad_u2[q_point][1];
    strain_ele[2] = 0;
    strain_ele[3] = (disp_grad_u1[q_point][1]+disp_grad_u2[q_point][0])/2;
    return strain_ele;
}
std::vector<double> stress(std::vector<double> strain){
    std::vector<double> stress_ele(4);
    TwoPhaseElastic tpe;
    double trace_strain = strain[0]+strain[1];
    std::vector<double> identity{1,1,1,0};
    for(int i=0; i<4; i++){
        stress_ele[i] = tpe.lambda*trace_strain*identity[i] + 2*tpe.G*(strain[i]);
    }
    return stress_ele;
}

void TwoPhaseElastic::stress_strain_calc(){
    FEValuesExtractors::Scalar u1(0);
    FEValuesExtractors::Scalar u2(1);
    QGauss<2> quadrature_formula(fe.degree + 1);
    FEValues<2> fe_values(fe, quadrature_formula, update_values | update_gradients);
    std::vector<Tensor<1, 2>> disp_grad_u1(quadrature_formula.size());
    std::vector<Tensor<1, 2>> disp_grad_u2(quadrature_formula.size());
    std::vector<double> strain_ele(4);
    std::vector<double> stress_ele(4);
    //old_solution_np = old_solution;
    for (auto &cell : dof_handler.active_cell_iterators()){
        if (cell->subdomain_id() == this_mpi_process) {
            fe_values.reinit(cell);
            fe_values[u1].get_function_gradients(old_solution, disp_grad_u1);
            fe_values[u2].get_function_gradients(old_solution, disp_grad_u2);
            for (unsigned int q_point = 0; q_point < quadrature_formula.size(); ++q_point) {
                //Calculating stress-strain
                strain_ele = strain(disp_grad_u1, disp_grad_u2, q_point);
                stress_ele = stress(strain_ele);
                // Storing stress-strain values to history variables
                for (int i = 0; i < 4; i++) {
                    quadrature_point_history[cell->active_cell_index()].strain[q_point][i] += strain_ele[i];
                    quadrature_point_history[cell->active_cell_index()].stress[q_point][i] += stress_ele[i];
                }
            }
        }
    }
}
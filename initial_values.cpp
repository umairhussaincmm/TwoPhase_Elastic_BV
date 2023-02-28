//
// Created by umair on 24/2/23.
//
#include "TwoPhaseElastic.h"
void InitialValues::vector_value(const Point<2> &p,
                                 Vector<double> & values) const {
    double u1_initial = 0.0;
    double u2_initial = 0.0;
    double mu_initial = .0001;
    double c_initial = 0.128;

    values(0)= u1_initial; //Initial u1 of domain
    values(1)= u2_initial; //Initial u2 of domain
    values(2)= mu_initial; //Initial mu of domain
    values(3)= c_initial; //Initial c of domain
}

void TwoPhaseElastic::applying_bc() {
    FEValuesExtractors::Scalar u1(0);
    FEValuesExtractors::Scalar u2(1);
    FEValuesExtractors::Scalar mu(2);
    FEValuesExtractors::Scalar c(3);

    ComponentMask u1_mask = fe.component_mask(u1);
    ComponentMask u2_mask = fe.component_mask(u2);
    ComponentMask c_mask = fe.component_mask(c);
    ComponentMask mu_mask = fe.component_mask(mu);

    std::map<types::global_dof_index,double> boundary_values;
    double c_boundary_value = 0.9;
    double mu_boundary_value = 0.0001;
    double u1_bv = 0;
    double u2_bv = 0;

    /*VectorTools::interpolate_boundary_values (dof_handler,
                                              2,
                                              ConstantFunction<2>(mu_boundary_value, 4),
                                              boundary_values,mu_mask);
    */
    VectorTools::interpolate_boundary_values (dof_handler,
                                              3,
                                              ConstantFunction<2>(c_boundary_value, 4),
                                              boundary_values,c_mask);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              4,
                                              ConstantFunction<2>(u1_bv, 4),
                                              boundary_values,u1_mask);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              4,
                                              ConstantFunction<2>(u2_bv, 4),
                                              boundary_values,u2_mask);

    //Creating temporary Non-Ghosted Vector
    PETScWrappers::MPI::Vector local_vector(locally_owned_dofs, mpi_communicator);
    local_vector=old_solution;
    for (auto &boundary_value : boundary_values)
        local_vector(boundary_value.first) = boundary_value.second;

    local_vector.compress(VectorOperation::insert);
    old_solution=local_vector;
}
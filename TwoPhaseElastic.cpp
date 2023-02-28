//
// Created by umair on 23/2/23.
//

#include "TwoPhaseElastic.h"

TwoPhaseElastic::TwoPhaseElastic()
        : M(10)
        , kappa(.0005)
        , c_max(1)
        , R(1)
        , T(1)
        , omega(2.6)
        , lambda(0)
        , G(100)
        , beta(0.01)
        , ref_conc(0.)
        , err(1e-9)
        , Farad(96500)
        , lam(.5)
        , A((1-lam))
        , B((-lam))
        , cm(2.29e4)
        , cL(1e3)
        , kre(1.9e-5)
        , n(Farad*kre*std::pow(cL,1-lam))
        , Vmax(18.5)
        , Vmin(-10.5)
        , mpi_communicator(MPI_COMM_WORLD)
        , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
        , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
        , pcout(std::cout, (this_mpi_process == 0))
        , triangulation(mpi_communicator,
                        typename Triangulation<2>::MeshSmoothing(
                                Triangulation<2>::smoothing_on_refinement | Triangulation<2>::smoothing_on_coarsening))
        , fe(FE_Q<2>(1), 4)//, FE_Q<2>(2), 3) //<2> -> dimension, (1) -> linear shape function, 4 -> Number of variables in problem
        , dof_handler(triangulation)
        , time(0.0)
        , final_time(1.)
        , time_step(.001)
        , theta(0.5)
        {}
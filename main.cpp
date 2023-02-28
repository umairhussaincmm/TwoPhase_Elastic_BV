#include <iostream>
#include "TwoPhaseElastic.h"

int main(int argc, char **argv) {
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    TwoPhaseElastic twophaseelastic;
    twophaseelastic.run();
    return 0;
}

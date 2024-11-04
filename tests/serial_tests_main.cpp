#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "better_mpi.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    auto result = Catch::Session().run(argc, argv);
    MPI_Finalize();

    return result;
}

#pragma once

#include <array>
#include <vector>

#include "better_mpi.h"

/**
 * \brief Wrapper class for a parallel sample sort
 *
 * \sa https://en.wikipedia.org/wiki/Samplesort
 */
class DistributedMemorySorter {
  public:
    using Coordinate = std::array<float, 3>;
    using Coordinates = std::vector<Coordinate>;

    /**
     * \brief Sort and evenly re-distribute coordinates
     *
     * \param coordinates Coordinates to be sorted
     * \param dimension Index of the coordinates to sort by
     * \param comm MPI communicator to use
     */
    static void sort_and_balance(Coordinates& coordinates, size_t dimension, MPI_Comm comm);

  private:
    /**
     *
     */
    DistributedMemorySorter();
    ~DistributedMemorySorter();

    MPI_Datatype mpi_coordinate_;

    /*
     * Defines the size of each group used for sub sampling in the Sample Sort
     * algorithm.\n\n
     * Example: for 16384 CPUs (full machine), if SUB_SAMPLING_GROUP_SIZE
     * is set as 2048, then we do sampling for [0..2043, 2044..4095, 4096.. 8172, etc].
     * After that, all 8 masters (1 for each group) will do the same procedure with
     * a unique master (rank 0);
     */
    constexpr static int SUB_SAMPLING_GROUP_SIZE = 128;
    constexpr static int SUB_SAMPLING_MINIMUM_CPUS = 4096;

    /**
     * \brief Performs a parallel Sample Sort.
     *
     * Checks consistency of element count and sum across all cores pre and
     * post sorting.
     *
     * \param coordinates Data to sort
     * \param sortDimension Index of the data to sort in
     * \param mpiComm Communicator to use
     */
    void sort(Coordinates& coordinates, size_t sortDimension, MPI_Comm mpiComm);


    /**
     * \brief Balances the amount of data across all the cores in the system
     *
     * Redistributes data such that all ranks have an equal amount of data.
     * Checks the sum and count of data elements before and after balancing
     * and ensures that the count on all cores is equal with a maximum
     * discrepancy of one.
     *
     * \param coordinates Data to balance
     * \param mpiComm Communicator to use
     */
    void balance(Coordinates& coordinates, MPI_Comm mpiComm);
};

#pragma once

#include "better_mpi.h"

class MemoryTools {
  public:
    /**
     * Configures the system to use a different method for handling memory errors when using the new
     * operator. This version will print the offending node's current memory usage
     */
    static void OverrideNewHandler();

    static size_t LocalMemoryUsageByte();

    /**
     * \return memory usage of the local node in Kb
     */
    static double LocalMemoryUsage();


    static size_t TotalAvailableMemoryByte();
    /**
     * returns memory available of the local node
     * \return
     */
    static double TotalAvailableMemoryMB();

    /**
     * Prints a 1 line report showing the maximum and minimum available
     * memory on the system nodes and the standard deviation and mean of
     * memory in all nodes, in the MPI_COMM_WORLD
     * \return memory usage (MB) of the CPU with max memory usage
     */
    static double AvailableMemory();

    /**
     * Prints a 1 line report showing the maximum and minimum available
     * memory on the system nodes and the standard deviation and mean of
     * memory in all nodes, in the parameter MPI_Comm
     * \param mpiComm MPI_Comm to report about
     * \return memory usage (MB) of the CPU with max memory usage
     */
    static double AvailableMemory(MPI_Comm mpiComm);
};

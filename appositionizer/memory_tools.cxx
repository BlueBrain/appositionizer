#include "memory_tools.h"  // Must be on the first substantial line.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>

#if defined(__APPLE__)
#include <mach/mach_init.h>
#include <mach/task.h>
#include <sys/sysctl.h>
#endif


/// global func to catch out of memory events when using the new operator.  Prints message with
/// current usage
void outOfMem() {
    // since this node, has run out of memory, have it report its usage info without MPI
    fprintf(stderr, "Out of Memory error");
    double localMB = MemoryTools::LocalMemoryUsage();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    mpi::show("Out of Memory Error on rank {}, current usage: {} MB.  Terminating", rank, localMB);

    abort();
}


//-----------------------------------------------------------------------------

void MemoryTools::OverrideNewHandler() {
    std::set_new_handler(outOfMem);
}

//-----------------------------------------------------------------------------

size_t MemoryTools::LocalMemoryUsageByte() {
// OSX has neither procfs nor CNK API
// we use the mach micro kernel API
#if defined(__APPLE__)
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, (task_info_t) &t_info, &t_info_count);
    return t_info.resident_size;

// for Linux only cluster we can use procfs
// Impossible to use on BG/Q due to CNK I/O mapping
#elif defined(__linux__)

    unsigned long rss;
    std::ifstream file;
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        file.open("/proc/self/statm");
        file >> rss;  // first is vmm
        file >> rss;  // second is rss
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("Impossible to open and parse /proc/self/statm : ") +
                                 e.what());
    }
    return rss * sysconf(_SC_PAGESIZE);

#else
#error "No valid Memory measurement API found"

#endif
}

double MemoryTools::LocalMemoryUsage() {
    return (static_cast<double>(LocalMemoryUsageByte()) / 1024.0 / 1024.0);
}


size_t MemoryTools::TotalAvailableMemoryByte() {
// Linux and FreeBSD supports sysconf
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
    const long pages = sysconf(_SC_PHYS_PAGES);
    const long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;

// MacOSX and other BSD derivative use classical sysctl interface
#elif defined(__APPLE__)
    int mib[2];
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    int64_t size = 0;
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, nullptr, 0) == 0)
        return static_cast<size_t>(size);
    throw std::runtime_error(
        std::string("sysctl() failure during MemoryTools::TotalAvailableMemoryByte()"));
#else
#error "Unable to find an API to check total available memory"
#endif
}

double MemoryTools::TotalAvailableMemoryMB() {
    return TotalAvailableMemoryByte() / (1024.0 * 1024.0);
}

//-----------------------------------------------------------------------------

// Returns the amount of free memory, in Mbytes
double MemoryTools::AvailableMemory() {
    return AvailableMemory(MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------------

double MemoryTools::AvailableMemory(MPI_Comm mpiComm) {
    int rank, size;
    MPI_Comm_rank(mpiComm, &rank);
    MPI_Comm_size(mpiComm, &size);

    double minUsageMB, maxUsageMB, avgUsageMB, stdevUsageMB;
    double usageMB = LocalMemoryUsage();

    int maxUsageRank, minUsageRank;

    struct {
        double usageMB;
        int rank;
    } myData, recvDataMin, recvDataMax;

    myData.rank = rank;
    myData.usageMB = usageMB;

    // Calculates minimum and maximum values for MB and its rank
    MPI_Reduce(&myData, &recvDataMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, mpiComm);
    MPI_Reduce(&myData, &recvDataMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, mpiComm);

    maxUsageMB = recvDataMax.usageMB;
    maxUsageRank = recvDataMax.rank;

    minUsageMB = recvDataMin.usageMB;
    minUsageRank = recvDataMin.rank;

    // Calculates average MB
    MPI_Reduce(&usageMB, &avgUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, mpiComm);
    avgUsageMB /= size;

    // Calculates std deviation
    MPI_Bcast(&avgUsageMB, 1, MPI_DOUBLE, 0, mpiComm);
    double diffSquared = (usageMB - avgUsageMB) * (usageMB - avgUsageMB);
    MPI_Reduce(&diffSquared, &stdevUsageMB, 1, MPI_DOUBLE, MPI_SUM, 0, mpiComm);
    stdevUsageMB = sqrt(stdevUsageMB / size);

    mpi::rank0::show(
        "Memory usage: Max {:.1f}MB [rank {:d}], Min {:.1f}MB [rank {:d}], Mean {:.1f}MB, stdDev "
        "{:.1f}MB",
        maxUsageMB,
        maxUsageRank,
        minUsageMB,
        minUsageRank,
        avgUsageMB,
        stdevUsageMB);

    double MBytes;
    unsigned long val, st[2];

    st[0] = 123456;
    val = (unsigned long) st;
    MBytes = (double) (val - (unsigned long) sbrk(0)) * 0.00000095367431640625;
    return (MBytes);
}

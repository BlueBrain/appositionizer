#include "distributed_memory_sorter.h"  // Must be on the first substantial line.

#include <numeric>

#include <oneapi/tbb/parallel_sort.h>

#include "instrumentor.h"

using namespace oneapi;

/**
 * \brief Auxiliary methods to cross-check results
 *
 * Used throughout the sorting and balancing process.
 */
namespace {

constexpr static long double EPSILON = 0.000001;

using Coordinate = std::array<float, 3>;
using Coordinates = std::vector<Coordinate>;
using data_t = Coordinate::value_type;

/**
 * \brief Helper class to store global count of values and their sum
 *
 * Used to ensure that the distributed values are consisent while
 * sorting.
 */
struct State {
    long long count = 0;
    long double sum = 0.0;

    State(const Coordinates& coordinates, MPI_Comm comm) {
        long long local_count = coordinates.size();
        long double local_sum = std::accumulate(
            coordinates.begin(), coordinates.end(), 0.0l, [](const long double s, const auto& c) {
                return s + std::accumulate(c.begin(), c.end(), 0.0l);
            });

        MPI_Allreduce(&local_sum, &sum, 1, MPI_LONG_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&local_count, &count, 1, MPI_LONG_LONG, MPI_SUM, comm);
    }
};

/**
 * \brief Checks that elements are distributed evenly across ranks
 *
 * Prints a warning if inconsistencies are detected.
 *
 * \param count Element count to verify
 * \param comm Communicator to use
 */
void checkBalancing(long long count, MPI_Comm comm) {
    long long minCount;
    long long maxCount;

    // Root gets the minimum and maximum number of elements
    MPI_Reduce(&count, &minCount, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
    MPI_Reduce(&count, &maxCount, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);

    if (maxCount - minCount > 1 && mpi::first(comm)) {
        mpi::warn(
            "Load balancing not accurate: the gap between "
            "the highest and lowest number of elements is "
            "{} (expected 0 or 1). (FLAG10)",
            maxCount - minCount);
    }
}


void checkConsistency(State pre_values, State post_values, MPI_Comm comm) {
    // basic checksum
    if (std::abs(pre_values.sum - post_values.sum) > EPSILON && mpi::first(comm)) {
        mpi::warn(
            "Load checksum balancing failed: "
            "the sum of all elements in the system before ({:.15f}) and after the sorting/balancing"
            "({:.15f}) differ - possibly due to numeric overflow or system dependent bug? (FLAG11)",
            pre_values.sum,
            post_values.sum);
    }

    if (pre_values.count != post_values.count) {
        mpi::abort(
            "Load balancing failed: "
            "the number of elements in the system before ({}) and after the sorting/balancing"
            "({}) differ - possibly due to numeric overflow or system dependent bug? (FLAG12)",
            pre_values.count,
            post_values.count);
    }
}

void checkSorting(const Coordinates& coordinates, size_t sortDimension, MPI_Comm comm) {
    const int mpiSize = mpi::size(comm);

    std::vector<data_t> globalMins;
    std::vector<data_t> globalMaxs;

    if (mpi::first(comm)) {
        globalMins.resize(mpiSize);
        globalMaxs.resize(mpiSize);
    }

    // All processes send to processor zero its min and max value (for checking)
    MPI_Gather(&(coordinates.front()[sortDimension]),
               1,
               mpi::datatype<data_t>(),
               globalMins.data(),
               1,
               mpi::datatype<data_t>(),
               0,
               comm);
    MPI_Gather(&(coordinates.back()[sortDimension]),
               1,
               mpi::datatype<data_t>(),
               globalMaxs.data(),
               1,
               mpi::datatype<data_t>(),
               0,
               comm);

    if (!mpi::first(comm)) {
        return;
    }

    for (size_t i = 0; i < globalMins.size(); i++) {
        if (globalMins[i] > globalMaxs[i]) {
            mpi::abort(
                "Sorting failed for dimension {}: rank {:d} has min value {:f} and max {:f}. "
                "(FLAG15)",
                sortDimension,
                i,
                globalMins[i],
                globalMaxs[i]);
        }

        if (i < static_cast<size_t>(mpiSize - 1) && globalMaxs[i] > globalMins[i + 1]) {
            mpi::abort(
                "Sorting failed for dimension {}: rank {:d} has max value {:f} and rank {:d} has "
                "min value "
                "{:f}. (FLAG16)",
                sortDimension,
                i,
                globalMaxs[i],
                i + 1,
                globalMins[i]);
        }
    }
}

}  // namespace

DistributedMemorySorter::DistributedMemorySorter() {
    MPI_Type_contiguous(sizeof(Coordinate), MPI_BYTE, &mpi_coordinate_);
    MPI_Type_commit(&mpi_coordinate_);
}

DistributedMemorySorter::~DistributedMemorySorter() {
    MPI_Type_free(&mpi_coordinate_);
}

void DistributedMemorySorter::sort_and_balance(Coordinates& coordinates,
                                               size_t dimension,
                                               MPI_Comm comm) {
    DistributedMemorySorter dms;
    dms.sort(coordinates, dimension, comm);
    dms.balance(coordinates, comm);
}

/******************************************************************/
/*               Work load balancing method                       */
/******************************************************************/

void DistributedMemorySorter::balance(Coordinates& coordinates, MPI_Comm comm) {
    auto phase_all = Instrumentor::phase("balance");
    // Initialization of MPI variabls
    int mpiSize = -1, mpiRank = -1;
    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    long long count = coordinates.size();
    const State pre_values(coordinates, comm);

    // each CPU gets the current number of elements of all the others
    std::vector<long long> coordinatesPerRank(mpiSize);
    MPI_Allgather(&count, 1, MPI_LONG_LONG, coordinatesPerRank.data(), 1, MPI_LONG_LONG, comm);

    // global index of beginning & end of my data (inclusive)
    long long myStartPos = 0;
    for (int i = 0; i < mpiRank; i++)
        myStartPos += coordinatesPerRank[i];
    long long myEndPos = myStartPos + coordinatesPerRank[mpiRank] - 1;

    // number of final elements per cpu
    std::vector<long long> elementsCountPerCpu(mpiSize);
    for (int i = 0; i < mpiSize; i++) {
        elementsCountPerCpu[i] = floor((long double) pre_values.count / mpiSize) +
                                 (i < (pre_values.count % mpiSize) ? 1 : 0);
    }

    // calculates the number of segments to be sent to each processor
    std::vector<int> sentElemsSize(mpiSize, 0);

    // global index of beginning & end of future slices (inclusive)
    long long indexStart = 0, indexEnd = -1;
    for (int i = 0; i < mpiSize; i++) {
        indexEnd += elementsCountPerCpu[i];

        // if the two intervals have regions in commmon (ie intersect)
        if (indexStart <= myEndPos && indexEnd >= myStartPos) {
            // if they intersect, the data to be sent is the size of the intersection
            //(+1 because beginning and end positions are inclusive
            sentElemsSize[i] = (int) (std::min(indexEnd, myEndPos) -
                                      std::max(indexStart, myStartPos) + 1);
        }

        indexStart += elementsCountPerCpu[i];
    }

    // calculates the offset on his data, for the received data from the other ranks
    std::vector<int> sentElemsOffset(mpiSize + 1, 0);
    std::partial_sum(sentElemsSize.begin(), sentElemsSize.end(), sentElemsOffset.begin() + 1);

    // calculate the amount of data received from each rank
    std::vector<int> recvElemsSize(mpiSize);
    MPI_Alltoall(sentElemsSize.data(), 1, MPI_INT, recvElemsSize.data(), 1, MPI_INT, comm);

    // calculate the offset of the data received
    std::vector<int> recvElemsOffset(mpiSize + 1, 0);
    std::partial_sum(recvElemsSize.begin(), recvElemsSize.end(), recvElemsOffset.begin() + 1);

    // calculate final size for data received and allocates memory
    int myRecvElemsSize = recvElemsOffset[mpiSize - 1] + recvElemsSize[mpiSize - 1];
    Coordinates received(myRecvElemsSize);

    MPI_Alltoallv(coordinates.data(),
                  sentElemsSize.data(),
                  sentElemsOffset.data(),
                  mpi_coordinate_,
                  received.data(),
                  recvElemsSize.data(),
                  recvElemsOffset.data(),
                  mpi_coordinate_,
                  comm);

    coordinates = std::move(received);

    const State post_values(coordinates, comm);
    checkConsistency(pre_values, post_values, comm);
    checkBalancing(coordinates.size(), comm);
}

/******************************************************************/
/*                 Sample Sort - MPI implementation               */
/******************************************************************/

void gatherSamples(std::vector<data_t>& samples, MPI_Comm comm) {
    const int numSamples = samples.size() - 1;  // Depends on global comm size (integer)

    int mpiRank, mpiSize;
    MPI_Comm_rank(comm, &mpiRank);
    MPI_Comm_size(comm, &mpiSize);

    std::vector<data_t> recvSamples;
    if (mpiRank == 0) {
        recvSamples.resize(mpiSize * numSamples);
    }

    /* gather the samples into the master */
    MPI_Gather(samples.data(),
               numSamples,
               mpi::datatype<data_t>(),
               recvSamples.data(),
               numSamples,
               mpi::datatype<data_t>(),
               0,
               comm);

    /* master selects the samples among all received samples, in order to broadcast them */
    if (mpiRank == 0) {
        tbb::parallel_sort(recvSamples.begin(), recvSamples.end());

        for (int i = 1; i <= numSamples; i++)
            samples[i - 1] = recvSamples[i * mpiSize - 1];
        samples[numSamples] = std::numeric_limits<data_t>::max();
    }
}

void DistributedMemorySorter::sort(Coordinates& coordinates, size_t dimension, MPI_Comm comm) {
    auto phase_all = Instrumentor::phase("sort");
    long long count = coordinates.size();

    int mpiSize, mpiRank;

    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    if (coordinates.size() == 0) {
        mpi::abort("Cannot continue sorting operation: rank {:d} has no segments at all", mpiRank);
    }

    const State pre_values(coordinates, comm);

    auto cmp = [=](const Coordinate& a, const Coordinate& b) {
        for (size_t n = dimension; n < a.size(); ++n) {
            if (a[n] < b[n]) {
                return true;
            } else if (a[n] > b[n]) {
                return false;
            }
        }
        return false;
    };

    /* sort local array */
    tbb::parallel_sort(coordinates.begin(), coordinates.end(), cmp);

    /* allocate memory for the arrays that will store the samples */
    std::vector<data_t> samples(mpiSize);

    /* select local p-1 equally spaced elements (samples only for the dimension we consider)*/
    for (int i = 1; i <= mpiSize; i++) {
        size_t pos = (count - 1) * i / mpiSize;
        samples[i - 1] = coordinates[pos][dimension];
    }

    // bug fix: as we couldnt allocate memory for N*N samples in only 1 CPU,
    // this method prevents that by using groups for sampling, i.e. a group of
    // CPU does sampling to a group master, and then group masters do sampling,
    // and send it to root 0 (master of all masters).
    if (mpiSize > SUB_SAMPLING_MINIMUM_CPUS) {
        if (mpi::first(comm)) {
            mpi::show("Using sampling by sub-groups, as MPI size is greater than {:d}...",
                      SUB_SAMPLING_MINIMUM_CPUS);
        }

        // split the original communicator into subgroups to gather samples locally
        MPI_Comm newComm;
        int mpiRankLocal = mpiRank % SUB_SAMPLING_GROUP_SIZE;
        int mpiGroupLocal = mpiRank / SUB_SAMPLING_GROUP_SIZE;
        MPI_Comm_split(comm, mpiGroupLocal, mpiRankLocal, &newComm);
        gatherSamples(samples, newComm);
        MPI_Comm_free(&newComm);

        // sub sampling is now complete, so masters of sub groups will do the same
        // step with master of the original comm group (i.e., here, we filter each
        // master rank of each group, and we use the group number as rank number)
        MPI_Comm_split(comm, mpiRankLocal, mpiGroupLocal, &newComm);
        if (mpiRankLocal == 0) {  // if it was a master before
            gatherSamples(samples, newComm);
        }
        MPI_Comm_free(&newComm);
    } else {
        gatherSamples(samples, comm);
    }

    /* now the samples array contains the global samples */
    MPI_Bcast(samples.data(), mpiSize, mpi::datatype<data_t>(), 0, comm);

    /* compute the number of elements that belong to each bucket */
    std::vector<int> bucketElemsSize(mpiSize);

    // ie gets the index of the bucket that contains our 1st coordinate
    size_t bucketIndex = 0;
    for (size_t i = 0; i < coordinates.size(); ++i) {
        for (; coordinates[i][dimension] >= samples[bucketIndex]; ++bucketIndex) {
        }
        ++bucketElemsSize[bucketIndex];
    }

    /* determine the starting location of each bucket's elements in the data array */
    std::vector<int> bucketOffsets(mpiSize);
    if (mpiSize > 0)
        bucketOffsets[0] = 0;
    for (int i = 1; i < mpiSize; i++)
        bucketOffsets[i] = bucketOffsets[i - 1] + bucketElemsSize[i - 1];
    // *3 in order to include the XYZ gap

    /* Perform an all2all communication to inform the corresponding processes */
    /* of the number of elements they are going to receive. */
    /* This information is stored in bucketElemsCount array */
    std::vector<int> recvElemsSize(mpiSize);
    MPI_Alltoall(bucketElemsSize.data(), 1, MPI_INT, recvElemsSize.data(), 1, MPI_INT, comm);

    /* Based on recvElemsCount determines where in the local array the data from each processor */
    /* will be stored. This array will store the received elements as well as the final */
    /* sorted sequence.*/
    std::vector<int> recvElemsOffset(mpiSize);
    if (mpiSize > 0)
        recvElemsOffset[0] = 0;
    for (int i = 1; i < mpiSize; i++)
        recvElemsOffset[i] = recvElemsOffset[i - 1] + recvElemsSize[i - 1];

    /* how many elements I will get */
    long long myRecvElemsSize = recvElemsOffset[mpiSize - 1] + recvElemsSize[mpiSize - 1];

    Coordinates sorted(myRecvElemsSize);
    /* Each process sends and receives the corresponding elements, using the MPI__Alltoallv */
    /* operation. The arrays bucketElemsCount and bucketOffsets are used to specify the number of
     * elements */
    /* to be sent and where these elements are stored, respectively. The arrays recvElemsCount */
    /* and recvElemsOffset are used to specify the number of elements to be received, and where
     * these */
    /* elements will be stored, respectively. */
    MPI_Alltoallv(coordinates.data(),
                  bucketElemsSize.data(),
                  bucketOffsets.data(),
                  mpi_coordinate_,
                  sorted.data(),
                  recvElemsSize.data(),
                  recvElemsOffset.data(),
                  mpi_coordinate_,
                  comm);

    /* perform the final local sort */
    tbb::parallel_sort(sorted.begin(), sorted.end(), cmp);

    coordinates.swap(sorted);

    /* Validates sorting*/
    const State post_values(coordinates, comm);
    checkConsistency(pre_values, post_values, comm);
    checkSorting(coordinates, dimension, comm);
}

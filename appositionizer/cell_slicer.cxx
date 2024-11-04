#include "cell_slicer.h"  // Must be on the first substantial line.

#include <assert.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <range/v3/numeric/accumulate.hpp>

#include "distributed_memory_sorter.h"
#include "instrumentor.h"
#include "better_mpi.h"
#include "geometry.h"

#include "defines.h"

namespace fs = std::filesystem;

using cell::Section;
using cell::SectionPoint;
using cell::Cell;

auto CellSlicer::bestDecomposition(size_t numberOfCpus) -> std::array<size_t, 3> {
    std::array<size_t, 3> slices;
    size_t tempSum = std::numeric_limits<size_t>::max();

    // TODO: shall the y be the last dimension to be sliced?
    // TODO: the y dimension should be changed from MAX_SLICES_PER_DIM to 1 once we slice the
    // cells based on number of segments instead of number of cells we want to minimize number
    // of slices in Y direction
    for (size_t y = 1; y <= 1; y++) {
        for (size_t z = y; z <= numberOfCpus; z++) {
            for (size_t x = y; x * y * z <= numberOfCpus; x++) {
                // if xyz decompose N, and are the best combination
                if (x * y * z == numberOfCpus && x + y + z < tempSum) {
                    slices[0] = x;
                    slices[1] = y;
                    slices[2] = z;
                    tempSum = x + y + z;
                }
            }
        }
    }

    if (tempSum == std::numeric_limits<size_t>::max())
        throw std::runtime_error("can't distribute load onto ranks, add/remove cores");

    mpi::rank0::show("- Volume will be sub-sliced in {}*{}*{}={} slices.",
                     slices[0],
                     slices[1],
                     slices[2],
                     numberOfCpus);

    return slices;
}

// In this function coordinates should be sorted
std::vector<SliceDetails> CellSlicer::getSubSlicesDetails(int dim,
                                                            SliceDetails universe,
                                                            int slicesCount,
                                                            const Coordinates& coordinates,
                                                            MPI_Comm mpiComm) {
    // return value: the array with all subslices of the universe (parameter)
    std::vector<SliceDetails> slicesDetails(slicesCount);

    // Initialization of variables
    int sliceRank = -1;
    MPI_Comm_rank(mpiComm, &sliceRank);

    int cpusPerSlice = universe.numberOfCPUs / slicesCount;

    // Verify that numberOfCPUs is a multiple of slicesCount
    // It should be by the way we compute slicesCount, but tests it anyway
    assert(slicesCount * cpusPerSlice == universe.numberOfCPUs);

    // populates return data...
    for (int s = 0; s < slicesCount; s++) {
        // number of CPUs:
        slicesDetails[s].numberOfCPUs = cpusPerSlice;

        // first and last rank of each slice:
        slicesDetails[s].firstRank = cpusPerSlice * s;
        slicesDetails[s].lastRank = cpusPerSlice * (s + 1) - 1;

        // coordinates of all sub slices are same as universe (for other dimensions)
        for (int d = 0; d < 3; d++) {
            if (d == dim)
                continue;
            slicesDetails[s].min[d] = universe.min[d];
            slicesDetails[s].max[d] = universe.max[d];
        }

        // coordinates of all subSlices for this dimension will be populated
        if (slicesDetails[s].firstRank == sliceRank)
            slicesDetails[s].min[dim] = coordinates.front()[dim];
        if (slicesDetails[s].lastRank == sliceRank)
            slicesDetails[s].max[dim] = coordinates.back()[dim];

        MPI_Bcast(&(slicesDetails[s].min[dim]), 1, MPI_DOUBLE, slicesDetails[s].firstRank, mpiComm);
        MPI_Bcast(&(slicesDetails[s].max[dim]), 1, MPI_DOUBLE, slicesDetails[s].lastRank, mpiComm);
    }

    // A bug fix will be necessary if we use cylinders (segment) instead of points (center of
    // cell). To avoid gaps between sub-cubes, neighbor ranks will agree on the same slice
    // location

    // Validates final data... in case sth went wrong before
    if (sliceRank == 0)
        for (int s = 0; s < slicesCount; s++) {
            // checks ranking order within same slice
            if (slicesDetails[s].firstRank > slicesDetails[s].lastRank) {
                mpi::warn(
                    "Validation failed for slicing: slice {} has first rank {} and last rank {}",
                    s,
                    slicesDetails[s].firstRank,
                    slicesDetails[s].lastRank);
            }

            // checks spatial ordering
            for (int d = 0; d < 3; d++) {
                if (slicesDetails[s].min[d] > slicesDetails[s].max[d]) {
                    mpi::warn(
                        "Validation failed for slicing: slice {} - ranks {} to {} - was calculated "
                        "between min {} and max {} for dim {}",
                        s,
                        slicesDetails[s].firstRank,
                        slicesDetails[s].lastRank,
                        slicesDetails[s].min[d],
                        slicesDetails[s].max[d],
                        d);
                }
            }

            // checks same conditions between neighbooring slices
            if (s == slicesCount - 1)
                continue;

            if (slicesDetails[s].max[dim] > slicesDetails[s + 1].min[dim]) {
                mpi::warn(
                    "Validation failed for slicing: for dim {}, slice {} has max {} but next slice "
                    "{} has min {}",
                    dim,
                    s,
                    slicesDetails[s].max[dim],
                    s + 1,
                    slicesDetails[s + 1].min[dim]);
            }

            if (slicesDetails[s].lastRank != slicesDetails[s + 1].firstRank - 1) {
                mpi::warn(
                    "Validation failed for slicing: slice {} has last rank {} and slice {} has "
                    "first rank {} - should be sequential",
                    s,
                    slicesDetails[s].lastRank,
                    s + 1,
                    slicesDetails[s + 1].firstRank);
            }
        }

    return slicesDetails;
}

auto CellSlicer::initialize(const cell::Group& cells)
    -> std::pair<SliceDetails, Coordinates> {
    std::pair<SliceDetails, Coordinates> result;

    auto& universe = result.first;
    auto& coordinates = result.second;

    // initializes the min and max values for each dimension
    double myMin[3], myMax[3];

    for (int d = 0; d < 3; d++) {
        myMin[d] = std::numeric_limits<double>::max();
        myMax[d] = std::numeric_limits<double>::min();
    }

    coordinates.reserve(cells.size());
    for (const auto& cell: cells) {
        const auto p = cell[0].extends().center();
        coordinates.push_back({p[0], p[1], p[2]});
        for (size_t i = 0; i < 3; ++i) {
            myMin[i] = std::min<double>(myMin[i], p[i]);
            myMax[i] = std::max<double>(myMax[i], p[i]);
        }
    }

    // Initial and final rank and positions for this universe
    universe.numberOfCPUs = comm_size_;
    universe.firstRank = 0;
    universe.lastRank = comm_size_ - 1;

    // gets the system's (ie global circuit's) minimum and maximum values
    for (int d = 0; d < 3; d++) {
        MPI_Allreduce(&(myMin[d]), &(universe.min[d]), 1, MPI_DOUBLE, MPI_MIN, comm_);
        MPI_Allreduce(&(myMax[d]), &(universe.max[d]), 1, MPI_DOUBLE, MPI_MAX, comm_);
    }
    return result;
}

void CellSlicer::writeSlicesToTextFile(const std::string& outputFolder,
                                         const SlicedCells& sliced_cells) {
    if (not mpi::first())
        return;  // only 1 rank will do the writing.

    fs::path outputFile{outputFolder};
    outputFile /= fmt::format("Slices.{}.txt", comm_size_);

    mpi::show("- Writing slices file {}...", outputFile.native());

    // opens file
    std::ofstream fSlices(outputFile);
    if (!fSlices.is_open()) {
        mpi::abort("Couldn't open file {} in order to write the slicing.", outputFile.native());
    }

    // writes header (info about the volume surrounding the cells)
    fSlices << "#universe bounding box:\n";
    fSlices << fmt::format("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n",
                           universe_.min[0],
                           universe_.max[0],
                           universe_.min[1],
                           universe_.max[1],
                           universe_.min[2],
                           universe_.max[2]);

    fSlices << "#slices: ";
    fSlices << comm_size_ << " (" << slicesCount_[0] << "*" << slicesCount_[1] << "*"
            << slicesCount_[2] << ")\n";
    //
    // writes slices
    fSlices << "#slice dimensions:\n";
    for (size_t r = 0; r < slicesCount_[0] * slicesCount_[1] * slicesCount_[2]; r++) {
        fSlices << fmt::format("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n",
                               cellSlices_[r].min[0],
                               cellSlices_[r].max[0],
                               cellSlices_[r].min[1],
                               cellSlices_[r].max[1],
                               cellSlices_[r].min[2],
                               cellSlices_[r].max[2]);
    }

    // writes slices
    fSlices << "#cpus bounding boxes:\n";
    for (int r = 0; r < comm_size_; r++) {
        fSlices << fmt::format("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n",
                               sliced_cells.getSlices()[r].min[0],
                               sliced_cells.getSlices()[r].max[0],
                               sliced_cells.getSlices()[r].min[1],
                               sliced_cells.getSlices()[r].max[1],
                               sliced_cells.getSlices()[r].min[2],
                               sliced_cells.getSlices()[r].max[2]);
    }

    // closes the file
    fSlices.close();
}

void CellSlicer::sortTileRecursive(int dim,
                                     SliceDetails universe,
                                     Coordinates& coordinates,
                                     MPI_Comm sliceComm) {
    // There's three communicators:
    //   `comm_` which encompasses all MPI ranks that partake in Appositionizer.
    //           Essentially, this is MPI_COMM_WORLD, possibly with a few ranks
    //           excluded.
    //
    //   `sliceComm` which has the MPI ranks of the current (dimension = `dim`)
    //               slice.
    //
    //   `subSliceComm` which has the MPI ranks of the next (dimension = `dim+1`)
    //                  slice.


    // recursivity stop condition
    // All cpus are now on their only slice, so they will broadcast this info
    // and get the info from all others
    if (dim == 3) {
        MPI_Allgather(&universe,
                      sizeof(SliceDetails),
                      MPI_BYTE,
                      cellSlices_.data(),
                      sizeof(SliceDetails),
                      MPI_BYTE,
                      comm_);
        return;
    }

    // Sorts cells in 1 dimension for this subvolume (coordinates and segmentsCount will be
    // updated)
    // Also balances workload so that data will be copied to the CPUs where the respective slice is.
    // Eg: for 8 CPUs, we have 8 slices, therefore:
    // For dim 0: data for slice 0 is in ranks [0..3], for slice 1 in ranks [4..7]
    // For dim 1: slice 0 in [0,1], slice 1 in [2,3], ..., slice 3 in [6.7]
    // For dim 2: slice 0 in [0], slice 1 in [1], ..., slice 7 in [7]
    DistributedMemorySorter::sort_and_balance(coordinates, dim, sliceComm);

    // populates sliceDetails array which holds the details of every sub-slice (and recursive call)
    auto slicesDetails =
        getSubSlicesDetails(dim, universe, slicesCount_[dim], coordinates, sliceComm);

    // ranks that are not part of this slice, continue...
    int sliceRank;
    MPI_Comm_rank(sliceComm, &sliceRank);

    // calls this method recursively for every single sub slice in the next dimension..
    for (size_t s = 0; s < slicesCount_[dim]; s++) {
        if (sliceRank < slicesDetails[s].firstRank || sliceRank > slicesDetails[s].lastRank)
            continue;

        // splits the actual MPI_Comm in sub MPI_Comms
        MPI_Comm subSliceComm;
        int newRank = sliceRank - slicesDetails[s].firstRank;  //[0..N];
        MPI_Comm_split(sliceComm, s /*new group code*/, newRank, &subSliceComm);

        // calls this functions for next dimension, for new sub group, and sub slice
        sortTileRecursive(dim + 1, slicesDetails[s], coordinates, subSliceComm);

        MPI_Comm_free(&subSliceComm);
    }
}

SlicedCells CellSlicer::redistribute(const cell::Group& cells) {
    auto phase_all = Instrumentor::phase("redistribute_cells");

    SlicedCells sliced_cells(comm_);

    {
        // array that will hold coordinates (one dimension at a time)
        // coordinates is on the format [x0,y0,z0, x1,y1,z1, ..., zN,yN,zN]
        Coordinates coordinates;
        std::tie(universe_, coordinates) = initialize(cells);

        mpi::rank0::show(
            "- Somas universe: x=[{:.2f} to {:.2f}], y=[{:.2f} to {:.2f}], z=[{:.2f} to {:.2f}].",
            universe_.min[0],
            universe_.max[0],
            universe_.min[1],
            universe_.max[1],
            universe_.min[2],
            universe_.max[2]);
        sortTileRecursive(0, universe_, coordinates, comm_);
    }

    mpi::rank0::info("Gathering local segments...");
    auto segments = localizeSegments(cells);
    mpi::rank0::info("Populating spatial index...");

    SectionPoint::BoundingBox mainBoundingBox = segments[0].boundingBox();
    for (const auto& seg: segments)
        if (seg.valid(segments))
            mainBoundingBox.merge(seg.boundingBox());

    sliced_cells.createIndex(std::move(segments));

    MPI_Allgather(&mainBoundingBox,
                  sizeof(mainBoundingBox),
                  MPI_BYTE,
                  sliced_cells.getSlices().data(),
                  sizeof(sliced_cells.getSlices().front()),
                  MPI_BYTE,
                  comm_);

    return sliced_cells;
}

std::vector<SectionPoint> CellSlicer::localizeSegments(const cell::Group& cells) {
    std::vector<SectionPoint> segments;
    using section_ref = std::reference_wrapper<const Section>;

    mpi::rank0::show("- Sending cells meta-data to the appropriate CPU, based on STR slicing...");

    std::vector<std::vector<section_ref>> sections(comm_size_);
    for (const auto& cell: cells) {
        bool assigned = false;
        auto center = cell[0].extends().center();

        for (int cpu = 0; cpu < comm_size_; cpu++) {
            SliceDetails& slice = cellSlices_[cpu];
            SectionPoint::BoundingBox box(slice.min, slice.max);

            if (box.contains(center)) {
                for (const auto& section: cell) {
                    sections[cpu].push_back(std::ref(section));
                }
                assigned = true;
                break;
            }
        }
        if (not assigned) {
            throw std::runtime_error("DA FUQ");
        }
    }

    // Count the number of section points to be sent to each rank
    std::vector<int> counts_send;
    counts_send.reserve(comm_size_);
    for (const auto& v_s: sections) {
        // TODO dangerous `int sum`.
        int sum = 0;
        for (const Section& s: v_s) {
            sum += s.size() + 1;  // b.size() is the number of segments
        }
        counts_send.push_back(sum);
    }

    std::vector<int> counts_recv(comm_size_);
    MPI_Alltoall(counts_send.data(), 1, MPI_INT, counts_recv.data(), 1, MPI_INT, comm_);

    size_t total_send = std::accumulate(counts_send.begin(), counts_send.end(), 0ul);
    size_t total_recv = std::accumulate(counts_recv.begin(), counts_recv.end(), 0ul);

    // TODO dangerous `std::vector<int> offsets_send`.
    std::vector<int> offsets_send(comm_size_ + 1, 0);
    std::partial_sum(counts_send.begin(), counts_send.end(), offsets_send.begin() + 1);

    std::vector<int> offsets_recv(comm_size_ + 1, 0);
    std::partial_sum(counts_recv.begin(), counts_recv.end(), offsets_recv.begin() + 1);

    MPI_Datatype section_point_t;
    mpi::ensure(MPI_Type_contiguous, sizeof(SectionPoint), MPI_BYTE, &section_point_t);
    mpi::ensure(MPI_Type_commit, &section_point_t);

    std::vector<SectionPoint> segments_send;

    segments_send.reserve(total_send);
    for (const auto& cpu_sections: sections) {
        for (const Section& section: cpu_sections) {
            std::vector<SectionPoint>& segs = const_cast<std::vector<SectionPoint>&>(section.raw());
            std::move(segs.begin(), segs.end(), std::back_inserter(segments_send));
        }
    }

    if (segments_send.size() != total_send) {
        mpi::abort("Expected to send {} segments, but found only {}.",
                   total_send,
                   segments_send.size());
    }

    segments.resize(total_recv);
    MPI_Alltoallv(segments_send.data(),
                  counts_send.data(),
                  offsets_send.data(),
                  section_point_t,
                  segments.data(),
                  counts_recv.data(),
                  offsets_recv.data(),
                  section_point_t,
                  comm_);

    {  // Rank 0 will print statistics
        size_t my_count = segments.size();
        std::vector<size_t> counts(comm_size_);
        MPI_Gather(&my_count,
                   1,
                   mpi::datatype<size_t>(),
                   counts.data(),
                   1,
                   mpi::datatype<size_t>(),
                   0,
                   comm_);

        if (mpi::first()) {
            std::vector<double> mem;
            std::transform(counts.begin(), counts.end(), std::back_inserter(mem), [](auto n) {
                return n / (1024.0 * 1024.0) * sizeof(SectionPoint);
            });
            double sum = std::accumulate(mem.begin(), mem.end(), 0.);
            double mean = sum / mem.size();
            double squares = std::accumulate(mem.begin(), mem.end(), 0., [mean](double s, int c) {
                return s + (c - mean) * (c - mean);
            });
            double sd = std::sqrt(squares / (counts.size() - 1));
            for (int i = 0; i < comm_size_; ++i) {
                mpi::info("Rank {:5d}: spatial index contains {:6.2f} MB of segments", i, mem[i]);
            }
            mpi::info("Overall:   spatial index contains {:6.2f} ± {:.2f} MB of segments",
                      mean,
                      sd);
        }
    }

    return segments;
}

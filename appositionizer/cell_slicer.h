#pragma once

#include <vector>

#include "better_mpi.h"
#include "geometry.h"

#include "data/cell_group.h"
#include "cell_slice.h"
#include "distributed_memory_sorter.h"
#include "data/segment.h"
#include "defines.h"

/**
 * Represents all boundaries around a slice, in terms of dimensions,
 * ranks and index inside of those ranks, that are handling it.
 * To be used by slicing detection, to specify the details of each slice and subslice
 */
typedef struct sliceDetails {
    double min[3];  // bottom left front corner of the subvolume
    double max[3];  // top right back corner of the subvolume

    int firstRank;  // rank where this slice starts (inclusive)
    int lastRank;   // rank where this slice ends (inclusive)

    int numberOfCPUs;  // number of CPUs handling this slice
} SliceDetails;

/**
 * Class responsible for slicing the initial spatial data into N subvolumes,
 * to be assigned to each CPU.Input can be either the txt file containing the
 * Slicing, or no input file which means that slicing will be computed on the fly.
 */
class CellSlicer {
  public:
    using Coordinates = DistributedMemorySorter::Coordinates;

    /**
     * Create a new slicer
     */
    explicit CellSlicer(MPI_Comm comm)
        : cellSlices_(mpi::size(comm))
        , slicesCount_(bestDecomposition(mpi::size(comm)))
        , comm_(comm)
        , comm_size_(mpi::size(comm)) {}

    /**
     * Reshuffle the sections of `cells` amongst the ranks available to
     * achieve a balanced distribution
     * \param cells The fully loaded (including morphologies) cells
     */
    SlicedCells redistribute(const cell::Group& cells);

    void writeSlicesToTextFile(const std::string& outputFolder,
                               const SlicedCells& sliced_cells);

  private:
    std::vector<SliceDetails> cellSlices_;
    SliceDetails universe_;

    const std::array<size_t, 3> slicesCount_; /** Number of slices for each dimension (x y z)*/

    MPI_Comm comm_;
    int comm_size_;

    /**
     * Fills the local segment storage with sections and segments the
     * current rank is responsible for
     * \param cells cell data to distribute for localization
     * \return segments localized in our rank
     */
    std::vector<cell::SectionPoint> localizeSegments(const cell::Group& cells);

    /// Returns the best slices combination in X, Y and Z dimensions, such that
    /// X*Y*Z=N, where N is the number of cpus.
    static auto bestDecomposition(size_t numberOfCpus) -> std::array<size_t, 3>;

    /**
     * Applies to Sort Title Recursive (STR) algorithm to the input volume.
     * The STR splits one dimension in N slices, then for each slices splits
     * it again in N sub-slices in the next dimension, and finally for the
     * 3rd dimensionm splits it again is sub-sub-slices\n
     * \n
     * The result will be written to Slices.AAA.txt text file, where AAA is the
     * number of CPUs.\n
     * \n
     * MORE INFO: https://bbpteam.epfl.ch/confluence/display/BBPDIAS/Load+Balancing+Summary \n
     *\n
     * CREDIT: Original serial algorithm by Farhan Tauheed - parallel version by Bruno Magalhaes
     *
     * \param dim Var used for recursive purposes (initialized as 0)
     * \param universe struct that represents the boundaries (in terms of ranks, indices and
     *dimension) of this volume
     * \param mySegmentsCount number of segments of this current CPU
     * \param nAppositionizers number of cpus (Apposition detectors) allocated to this slice
     * \param coordinates coordinates for this slice
     */
    void sortTileRecursive(int dim,
                           SliceDetails universe,
                           Coordinates& coordinates,
                           MPI_Comm mpiComm);

    /**
     * Collects all coordinates to be considered and their maximum extends.
     * \param cells fully loaded cells to extract the coordinates from
     * \return a pair containing the extends and coordinates to consider
     */
    auto initialize(const cell::Group& cells) -> std::pair<SliceDetails, Coordinates>;

    /**
     * Returns the details for all subslices of the parameter universe.
     * \param dim dimension [0..2]
     * \param universe universe data (ie, previous dimension slice)
     * \param slicesCount total number of sub slices
     * \param globalSegmentsCount total number of segments in universe
     * \param mySegmentsCount total number of segments in this rank
     * \param coordinates coordinates array
     * \param mpiComm MPI comm
     * \return
     */
    std::vector<SliceDetails> getSubSlicesDetails(int dim,
                                                  SliceDetails universe,
                                                  int slicesCount,
                                                  const Coordinates& coordinates,
                                                  MPI_Comm mpiComm);

    /**
     * MPI collective operation and must be called by all processes in `comm_`.
     * Each process send his slice data to rank 0, who will write it in the
     * slice format on the file 'filename'\n
     * \n
     * File content (text file):\n
     * 0: min max (for dimension 0)\n
     * 1: min max (for dimension 1)\n
     * 2: min max (for dimension 2)\n
     * 3: N (number of slices)\n
     * 4: minX minY minZ maxX maxY maxZ (corners defining slice 1)\n
     * 5: minX minY minZ maxX maxY maxZ (corners defining slice 2)\n
     * 6: minX minY minZ maxX maxY maxZ (corners defining slice 3)\n
     * ...\n
     *
     * \param slice current slice
     * \param filename file name - usually Slices.ABC.txt (where ABC is the number of processes)
     */
    void writeSlicingDataToFile(SliceDetails slice, char* filename);
};

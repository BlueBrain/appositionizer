#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "../better_mpi.h"
#include "../data/cell_group.h"
#include "../cell_slice.h"

namespace fs = std::filesystem;

namespace process {

/**
 * \brief Static configuration settings for the cell processing
 */
struct Configuration {
    const fs::path tmp_dir;
    const fs::path output_dir;
    const std::string output_population;
    const std::string source_population;
    const std::string target_population;
};

/**
 * \brief Performs cell overlap detection with independent MPI ranks
 *
 * This mode of the Appositionizer will avoid most MPI interactions to
 * process cell overlap detection.
 *
 * With an initial MPI exchange to determine which cells overlap the
 * current rank, every rank will individually load cells, detect where
 * appositions are formed in overlaps with the rank-local spatial index
 * content, and write the output_dir to a rank-local file.
 *
 * Arguments:
 *  \param cell_group Fair-share group of pre-synaptic cells, with
 *    morphologies loaded.
 *
 *  \param tmp_dir Path of the directory containing all temporary output_dir.
 *    Sometimes a scratch pad is available, e.g. on BB5 NVME drives. The
 *    temporary files will only be read back from the same MPI rank they
 *    were written with.
 *
 *  \param output_dir Path of the directory containing all final output_dir. This
 *    path must refer to the same directory on all MPI ranks.
 *
 *  \param job_size The inner loop that computes the appositions for each
 *  pre-synpatic cell is parallelized by dividing the range of pre-GIDs into
 *  subranges or "jobs" of size `jobsize`.
 *
 */
void disentangled(const cell::Group& cell_group,
                  const SlicedCells& sliced_cells,
                  const Configuration& cfg,
                  size_t job_size,
                  const InterBoutonInterval& interval,
                  MPI_Comm comm);

}  // namespace process

#pragma once

#include <map>
#include <vector>

#include "better_mpi.h"
#include "io/structures.h"

namespace io {

/**
 * \brief Read files with Apposition data in binary format
 */
class FileReader {
  public:
    /**
     * \brief Construct a reader, reading the index file only
     */
    FileReader(const std::string& dirname, const std::string& extension);

    /** \brief Returns count and offset. */
    std::pair<size_t, size_t> count_and_offset(size_t cell_id) const;

    /**
     * \brief Returns the Apposition count with \a cell_id as the pre-synaptic cell
     */
    size_t count(size_t cell_id) const;

    /**
     * \brief Returns the file offset for the pre-synaptic cell \a cell_id.
     *
     * If \a cell_id is not in the index, return the offset of the next larger
     * GID. If the file is sorted, the appositions for pre GIDs `n, ..., m` (exclusive)
     * occupy the range `offset(n), ..., offset(m)` (exclusive) within the file.
     * This is true even if there are no appositions for `cell_id`.
     */
    size_t offset(size_t cell_id) const;

    /**
     * \brief Returns the total Apposition count for all cells in the index.
     */
    uint64_t count_total() const;

    /**
     * \brief For sorted files, returns all appositions within the range of pre-GIDs.
     *
     * This method reads all appositions with pre-GIDs \a gid_begin, ..., \a gid_end
     * (exclusive) without issueing individual read operations for each GID.
     * This requires that the file is sorted by pre-GID in ascending order.
     */
    std::vector<DataEntry> get_range(size_t gid_begin, size_t gid_end) const;

  private:
    std::map<size_t, IndexEntry> index_;
    mpi::File datafile_;
};

/** \brief Read an index of a appositions file. */
std::map<size_t, io::IndexEntry> read_index(const std::string& filename);

}  // namespace io

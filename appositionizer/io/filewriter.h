#pragma once

#include <tuple>
#include <vector>

#include "data/apposition.h"
#include "better_mpi.h"

#include "io/iowriter.h"
#include "io/structures.h"

namespace cell {

class Section;
class Group;
class Segment;

}  // namespace cell

namespace io {

/**
 * \brief Generate file extensions based on the MPI communicator \a comm
 *
 * If the current rank is the first rank of \a comm, will return the rank in the global
 * communicator, otherwise \c .0 .
 */
std::string generate_extension(MPI_Comm comm);

class FileWriter: public IOWriter {
  public:
    FileWriter(const std::string& outputFolder,
               size_t nextIndexOffset,
               size_t nextDataOffset,
               MPI_Comm comm = MPI_COMM_WORLD);
    ~FileWriter();

    /** \brief Give a group of cells and its appositions regions, write it to the appositions file
     *  \param appositions appositions indexed by pre & post gid
     *  \return A pair containing the number of index and data entries created */
    auto write(const std::map<size_t, std::map<size_t, std::vector<io::DataEntry>>>& appositions)
        -> std::tuple<size_t, size_t, size_t, size_t>;

    /** \brief Give a group of cells and its appositions regions, write it to the appositions file
     *  \param appositions appositions sorted by pre & post gid
     *  \param meta Index of the appositions, with 0-based connectionsBinaryAddress.
     *  \return A tuple containing the number of index and data entries created */
    auto write(const std::vector<io::DataEntry>& appositions, std::vector<IndexEntry> meta)
        -> std::tuple<size_t, size_t, size_t, size_t>;

  private:
    void print_summary();
    static mpi::File open_file(const std::string& filename, MPI_Comm comm);
    static mpi::File open_index_file(const std::string& dirname, MPI_Comm comm);
    static mpi::File open_payload_file(const std::string& dirname, MPI_Comm comm);

    MPI_Comm comm_;

    mpi::File mpiEdgeFile_;
    mpi::File mpiEdgeDataFile_;

    mpi::SyncCounter<2, size_t> offsets_;
    const std::string dirname_;

    /**
     * \brief Write the header information of the index.
     */
    void write_index_header();
};

/**
 * \brief Write the header information of the index.
 */
void write_index_header(mpi::File& file, size_t n_entries);

/** \brief Write \a index to disk. */
void write_index(const std::string& filename, const std::vector<io::IndexEntry>& index);


}  // namespace io

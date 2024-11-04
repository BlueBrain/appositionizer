#include <cmath>
#include <filesystem>
#include <numeric>
#include <sstream>

#include "instrumentor.h"
#include "data/cell_group.h"

#include "defines.h"

#include "filewriter.h"
#include "structures.h"

namespace fs = std::filesystem;

using apposition::Region;

const int META_OFFSET = 0;
const int DATA_OFFSET = 1;

namespace io {

FileWriter::FileWriter(const std::string& output_dir,
                       size_t nextIndexOffset,
                       size_t nextDataOffset,
                       MPI_Comm comm)
    : comm_(comm)
    , mpiEdgeFile_(open_index_file(output_dir, comm))
    , mpiEdgeDataFile_(open_payload_file(output_dir, comm))
    , offsets_({nextIndexOffset, nextDataOffset}, comm)
    , dirname_(output_dir) {}

mpi::File FileWriter::open_file(const std::string& filename, MPI_Comm comm) {
    return mpi::File::open(filename, mpi::File::Mode::WRITE, comm);
}

mpi::File FileWriter::open_index_file(const std::string& dirname, MPI_Comm comm) {
    return open_file(io::index_filename(dirname, comm), comm);
}

mpi::File FileWriter::open_payload_file(const std::string& dirname, MPI_Comm comm) {
    return open_file(io::payload_filename(dirname, comm), comm);
}

FileWriter::~FileWriter() {
    write_index_header();
    print_summary();
}

void FileWriter::print_summary() {
    size_t cells = offsets_.add(0, META_OFFSET);
    size_t appositions = offsets_.add(0, DATA_OFFSET);

    const auto ext = generate_extension(comm_);
    mpi::rank0::info(
        "Wrote output to {}/{{{},{}}}{}\n"
        "          Number of cells written: {}\n"
        "          Number of edges written: {}",
        dirname_,
        INDEX_STEM,
        PAYLOAD_STEM,
        ext,
        cells,
        appositions);
}

void FileWriter::write_index_header() {
    const auto total_count = offsets_.add(0, META_OFFSET);
    // only rank 0 must write the header
    if (mpi::first(comm_)) {
        ::io::write_index_header(mpiEdgeFile_, total_count);
    }
}

auto FileWriter::write(const std::map<size_t, std::map<size_t, std::vector<DataEntry>>>& appositions)
    -> std::tuple<size_t, size_t, size_t, size_t> {
    std::vector<IndexEntry> meta;
    std::vector<DataEntry> data;

    Instrumentor::phase_begin("flatten map<map<..>>");
    meta.reserve(appositions.size());

    size_t binaryAddress = 0ul;
    for (const auto& [gid, mapped_appositions]: appositions) {
        size_t count = 0;
        for (const auto& [post_gid, actual_appositions]: mapped_appositions) {
            count += actual_appositions.size();
            std::copy(actual_appositions.begin(), actual_appositions.end(), std::back_inserter(data));
        }

        meta.emplace_back(gid, count, binaryAddress);
        binaryAddress += count * sizeof(DataEntry);
    }
    Instrumentor::phase_end("flatten map<map<..>>");

    return write(data, meta);
}

auto FileWriter::write(const std::vector<DataEntry>& data, std::vector<IndexEntry> meta)
    -> std::tuple<size_t, size_t, size_t, size_t> {
    size_t meta_offset = offsets_.add(meta.size(), META_OFFSET);
    size_t data_offset = offsets_.add(data.size(), DATA_OFFSET);

    // we set offsets as 0 before, therefore now we will add it to the correct one
    for (auto& m: meta) {
        m.connectionsBinaryAddress += data_offset * sizeof(DataEntry);
    }

    Instrumentor::phase_begin("MPI_File_write_at");

    if (not meta.empty()) {
        size_t offset = sizeof(IndexHeader) + sizeof(IndexEntry) * meta_offset;
        mpiEdgeFile_.write_array_at(meta, offset);
    }

    if (not data.empty()) {
        size_t offset = sizeof(DataEntry) * data_offset;
        mpiEdgeDataFile_.write_array_at(data, offset);
    }

    Instrumentor::phase_end("MPI_File_write_at");

    return {meta.size(), meta_offset, data.size(), data_offset};
}

void write_index_header(mpi::File& file, size_t n_entries) {
    io::IndexHeader index_header;
    strncpy(index_header.version, APPOSITIONIZER_VERSION_STRING, 16);
    index_header.numberOfEntries = n_entries;

    file.write_scalar_at(index_header, 0ul);
}

void write_index(const std::string& filename, const std::vector<io::IndexEntry>& index) {
    auto index_file = mpi::File::open(filename, mpi::File::Mode::WRITE, MPI_COMM_SELF);

    write_index_header(index_file, index.size());

    index_file.write_array_at(index, sizeof(io::IndexHeader));
}


}  // namespace io

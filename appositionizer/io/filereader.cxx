#include "io/filereader.h"

#include <filesystem>
#include <fstream>

#include "instrumentor.h"
#include "defines.h"

namespace fs = std::filesystem;

namespace {

mpi::File open_datafile(fs::path filename, const fs::path& extension) {
    filename /= io::PAYLOAD_STEM;
    filename += extension;

    return mpi::File::open(filename, mpi::File::Mode::READ, MPI_COMM_SELF);
}

}  // namespace


namespace io {

// TODO note, that the return type uses the GID as a key. GIDs are commonly `int` not `size_t`.
std::map<size_t, io::IndexEntry> read_index(const std::string& filename) {
    std::map<size_t, io::IndexEntry> result;

    auto f = mpi::File::open(filename, mpi::File::Mode::READ, MPI_COMM_SELF);
    MPI_Offset size;
    mpi::ensure(MPI_File_get_size, *f, &size);
    io::IndexHeader header;
    f.read_scalar_at(header, 0ul);

    size_t have = static_cast<size_t>(size);
    size_t expected = sizeof(header) + header.numberOfEntries * sizeof(io::IndexEntry);
    if (have != expected) {
        mpi::abort("Mismatch of size information for {}: got {}, expected {}",
                   filename.c_str(),
                   have,
                   expected);
    }

    std::vector<io::IndexEntry> entries(header.numberOfEntries);
    f.read_array_at(std::span{entries.data(), size_t(header.numberOfEntries)}, sizeof(header));

    for (const auto& entry: entries) {
        result[entry.cellId] = entry;
    }

    return result;
}

std::map<size_t, io::IndexEntry> read_index(fs::path filename, const fs::path& extension) {
    filename /= io::INDEX_STEM;
    filename += extension;

    return read_index(filename);
}

FileReader::FileReader(const std::string& dirname, const std::string& extension)
    : index_(read_index(dirname, extension))
    , datafile_(open_datafile(dirname, extension)) {}


std::pair<size_t, size_t> FileReader::count_and_offset(size_t gid) const {
    if (index_.empty()) {
        return {0ul, 0ul};
    }

    auto it = index_.lower_bound(gid);

    if (it == index_.end()) {
        const auto& back = index_.rbegin()->second;
        size_t addr = back.connectionsBinaryAddress + back.synapseCount * sizeof(DataEntry);
        return {0ul, addr};
    } else {
        const auto& entry = it->second;
        size_t addr = entry.connectionsBinaryAddress;

        // TODO improve int correctness for `int(gid)`.
        size_t count = entry.cellId == int(gid) ? entry.synapseCount : 0ul;

        return {count, addr};
    }
}

size_t FileReader::count(size_t cell_id) const {
    return count_and_offset(cell_id).first;
}

size_t FileReader::offset(size_t cell_id) const {
    return count_and_offset(cell_id).second;
}

uint64_t FileReader::count_total() const {
    uint64_t total = 0;
    for (auto kv: index_) {
        total += this->count(kv.first);
    }
    return total;
}

std::vector<DataEntry> FileReader::get_range(size_t pre_gid_begin, size_t pre_gid_end) const {
    if (pre_gid_begin >= pre_gid_end) {
        return {};
    }

    auto file_begin = offset(pre_gid_begin);
    auto file_end = offset(pre_gid_end);

    auto count = (file_end - file_begin) / sizeof(DataEntry);
    auto data = std::vector<DataEntry>(count);
    datafile_.read_array_at(data, file_begin);

    return data;
}


}  // namespace io

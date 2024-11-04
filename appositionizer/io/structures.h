#pragma once

#include <filesystem>

#include "../better_mpi.h"

namespace cell {

class Section;
class Group;
class Segment;

}  // namespace cell

namespace apposition {

class Region;
}

/**
 * Data structures to be written to Apposition files.
 *
 * \warning Any modification to these structures will render the binary files incompatible
 * with previous versions.
 */
namespace io {

const static std::string INDEX_STEM = "edges";
const static std::string PAYLOAD_STEM = "edgesData";

inline std::filesystem::path join_path(const std::filesystem::path& dirname,
                                       const std::filesystem::path& basename) {
    return dirname / basename;
}

/** \brief Internal method for computing the 'extension' of data/index files.
 */
inline std::string generate_extension(MPI_Comm comm) {
    // The assumption seems to be that there are two modes:
    //  a) the one file per MPI rank
    //  b) the single global file
    // When a) is required `comm == MPI_COMM_SELF` otherwise it's
    // just the communicator.
    std::stringstream oss;
    oss << "." << (mpi::first(comm) ? mpi::rank(MPI_COMM_WORLD) : 0);
    return oss.str();
}


/** \brief The name of the file containing the index.
 *
 * See `payload_filename` for the file containing the heavy data
 * that was indexed.
 */
inline std::string index_filename(const std::string& dirname, MPI_Comm comm) {
    return join_path(dirname, INDEX_STEM + generate_extension(comm));
}

/** \brief The name of the file containing the heavy data/payload.
 *
 * See `index_filename` for the file containing the index of this data.
 */
inline std::string payload_filename(const std::string& dirname, MPI_Comm comm) {
    return join_path(dirname, PAYLOAD_STEM + generate_extension(comm));
}


/**
 * \brief "Global" metadata information forming the header of appositions.0
 */
struct IndexHeader {
    double architectureIdentifier = 1.001;
    long long numberOfEntries = 0;
    char version[16] = "";
};

/**
 * \brief Metadata information per cell in appositions.0
 *
 * To be written after the header block
 */
struct IndexEntry {
    int cellId = -1;
    unsigned int synapseCount = 0;
    long long connectionsBinaryAddress = 0;

    IndexEntry() = default;
    IndexEntry(int cell_id, unsigned int synapse_count, long long connections_addr)
        : cellId(cell_id)
        , synapseCount(synapse_count)
        , connectionsBinaryAddress(connections_addr){};
    IndexEntry(const IndexEntry&) = default;

    void add(const apposition::Region& r);

    bool compatible(const apposition::Region& r);

    bool empty() const {
        return synapseCount == 0 and cellId < 0;
    }
};

/**
 * \brief Dataformat of a serialized Apposition, written to AppositionsData.0
 */
struct DataEntry {
    // Use one unsigned char to store two section types (2 bits wide)
    static constexpr size_t SECTION_MASK = 0xF;
    static constexpr size_t SECTION_SHIFT = 4;

    int preCell;
    int preSection;
    int preSegment;
    int postCell;
    int postSection;
    int postSegment;
    int sectionOrder;
    float distance_to_soma;
    float efferent_segment_offset;
    float afferent_segment_offset;
    float preSectionFraction;
    float postSectionFraction;
    float preSurfacePosition[3];
    float postCenterPosition[3];
    float spineLength;
    unsigned char sectionType;
    float preCenterPosition[3];
    float postSurfacePosition[3];

    DataEntry() = default;
    DataEntry(const cell::Segment& pre,
              const cell::Segment& post,
              const cell::Section& pre_section,
              unsigned char pre_fraction,
              unsigned char post_fraction);
};

inline bool operator<(const DataEntry& a, const DataEntry& b) {
    auto a_tuple = std::tuple{a.preCell,
                              a.postCell,
                              a.preSection,
                              a.preSegment,
                              a.preSectionFraction,
                              a.postSection,
                              a.postSegment,
                              a.postSectionFraction};

    auto b_tuple = std::tuple{b.preCell,
                              b.postCell,
                              b.preSection,
                              b.preSegment,
                              b.preSectionFraction,
                              b.postSection,
                              b.postSegment,
                              b.postSectionFraction};

    return a_tuple < b_tuple;
}

}  // namespace io

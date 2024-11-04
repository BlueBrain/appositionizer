#pragma once

#include "better_mpi.h"
#include "io/iowriter.h"
#include "io/structures.h"
#include <hdf5.h>
#include <highfive/H5File.hpp>

namespace io {

static const std::size_t SECTION_OFFSET = 1;

struct DataVectors {
    std::vector<uint32_t> preCell;
    std::vector<uint32_t> preSection;
    std::vector<uint32_t> preSegment;
    std::vector<uint32_t> postCell;
    std::vector<uint32_t> postSection;
    std::vector<uint32_t> postSegment;
    std::vector<uint32_t> sectionOrder;
    std::vector<float> distance_to_soma;
    std::vector<float> efferent_segment_offset;
    std::vector<float> afferent_segment_offset;
    std::vector<float> preSectionFraction;
    std::vector<float> postSectionFraction;
    std::vector<float> preSurfacePosition_1;
    std::vector<float> preSurfacePosition_2;
    std::vector<float> preSurfacePosition_3;
    std::vector<float> postCenterPosition_1;
    std::vector<float> postCenterPosition_2;
    std::vector<float> postCenterPosition_3;
    std::vector<float> spineLength;
    std::vector<uint32_t> preSectionType;
    std::vector<uint32_t> postSectionType;
    std::vector<float> preCenterPosition_1;
    std::vector<float> preCenterPosition_2;
    std::vector<float> preCenterPosition_3;
    std::vector<float> postSurfacePosition_1;
    std::vector<float> postSurfacePosition_2;
    std::vector<float> postSurfacePosition_3;
    explicit DataVectors(size_t size);
};


class SonataWriter: public IOWriter {
  public:
    SonataWriter(const std::string& fileName,
                 const std::string& populationName,
                 const std::string& source_population,
                 const std::string& target_population,
                 uint64_t totalRecords,
                 const MPI_Comm& comm);

    auto write(const std::vector<io::DataEntry>& appositions, std::vector<IndexEntry> meta)
        -> std::tuple<size_t, size_t, size_t, size_t>;
    //    static HighFive::DataTransferProps TRANSFER_PROPS;


  protected:
    const MPI_Comm& mpi_comm_;
    HighFive::File file_;
    const std::string population_name_;
    HighFive::Group population_group_;
    HighFive::Group properties_group_;
    const uint64_t n_records_;

    mpi::SyncCounter<1, size_t> offsets_;


  private:
    HighFive::File open_file(const std::string& fileName);
    HighFive::Group get_create_population_group();
    HighFive::Group get_create_properties_group();
    void create_datasets(const std::string& source_population,
                         const std::string& target_population);
    static DataVectors transform_to_vectors(const std::vector<DataEntry>& appositions);
};

}  // namespace io

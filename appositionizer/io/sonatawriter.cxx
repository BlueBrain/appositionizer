//
// Created by temerev on 4/4/24.
//

#include "sonatawriter.h"

namespace io {

SonataWriter::SonataWriter(const std::string& fileName,
                           const std::string& populationName,
                           const std::string& source_population,
                           const std::string& target_population,
                           const uint64_t totalRecords,
                           const MPI_Comm& comm)
    : mpi_comm_(comm)
    , file_(open_file(fileName))
    , population_name_(populationName)
    , population_group_(get_create_population_group())
    , properties_group_(get_create_properties_group())
    , n_records_(totalRecords)
    , offsets_({0}) {
    create_datasets(source_population, target_population);
}

HighFive::Group SonataWriter::get_create_population_group() {
    return file_.createGroup("edges").createGroup(population_name_);
}

HighFive::Group SonataWriter::get_create_properties_group() {
    return population_group_.createGroup("0");
}

void SonataWriter::create_datasets(const std::string& source_population,
                                   const std::string& target_population) {
    population_group_.createDataSet<uint32_t>("edge_type_id", HighFive::DataSpace({n_records_}));
    auto src = population_group_.createDataSet<uint32_t>("source_node_id",
                                                         HighFive::DataSpace({n_records_}));
    src.createAttribute("node_population", source_population);
    auto dst = population_group_.createDataSet<uint32_t>("target_node_id",
                                                         HighFive::DataSpace({n_records_}));
    dst.createAttribute("node_population", target_population);
    properties_group_.createDataSet<uint32_t>("efferent_section_id",
                                              HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<uint32_t>("efferent_segment_id",
                                              HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<uint32_t>("afferent_section_id",
                                              HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<uint32_t>("afferent_segment_id",
                                              HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<uint32_t>("section_order", HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("distance_soma", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("efferent_segment_offset",
                                           HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_segment_offset",
                                           HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("efferent_section_pos",
                                           HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_section_pos",
                                           HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("efferent_surface_x", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("efferent_surface_y", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("efferent_surface_z", HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("afferent_center_x", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_center_y", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_center_z", HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("spine_length", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<uint32_t>("efferent_section_type",
                                              HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<uint32_t>("afferent_section_type",
                                              HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("efferent_center_x", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("efferent_center_y", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("efferent_center_z", HighFive::DataSpace({n_records_}));

    properties_group_.createDataSet<float>("afferent_surface_x", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_surface_y", HighFive::DataSpace({n_records_}));
    properties_group_.createDataSet<float>("afferent_surface_z", HighFive::DataSpace({n_records_}));
}


auto SonataWriter::write(const std::vector<io::DataEntry>& appositions, std::vector<IndexEntry> meta)
    -> std::tuple<size_t, size_t, size_t, size_t> {
    auto transfer_props = HighFive::DataTransferProps{};
    transfer_props.add(HighFive::UseCollectiveIO{});

    DataVectors data = transform_to_vectors(appositions);  // in-memory
    size_t offset = offsets_.add(appositions.size(), 0);
    size_t count = appositions.size();
    std::vector<uint32_t> zeros(count, 0);

    population_group_.getDataSet("edge_type_id")
        .select({offset}, {count})
        .write(zeros, transfer_props);
    population_group_.getDataSet("source_node_id")
        .select({offset}, {count})
        .write(data.preCell, transfer_props);
    properties_group_.getDataSet("efferent_section_id")
        .select({offset}, {count})
        .write(data.preSection, transfer_props);
    properties_group_.getDataSet("efferent_segment_id")
        .select({offset}, {count})
        .write(data.preSegment, transfer_props);
    population_group_.getDataSet("target_node_id")
        .select({offset}, {count})
        .write(data.postCell, transfer_props);
    properties_group_.getDataSet("afferent_section_id")
        .select({offset}, {count})
        .write(data.postSection, transfer_props);
    properties_group_.getDataSet("afferent_segment_id")
        .select({offset}, {count})
        .write(data.postSegment, transfer_props);

    properties_group_.getDataSet("section_order")
        .select({offset}, {count})
        .write(data.sectionOrder, transfer_props);

    properties_group_.getDataSet("distance_soma")
        .select({offset}, {count})
        .write(data.distance_to_soma, transfer_props);
    properties_group_.getDataSet("efferent_segment_offset")
        .select({offset}, {count})
        .write(data.efferent_segment_offset, transfer_props);
    properties_group_.getDataSet("afferent_segment_offset")
        .select({offset}, {count})
        .write(data.afferent_segment_offset, transfer_props);

    properties_group_.getDataSet("efferent_section_pos")
        .select({offset}, {count})
        .write(data.preSectionFraction, transfer_props);
    properties_group_.getDataSet("afferent_section_pos")
        .select({offset}, {count})
        .write(data.postSectionFraction, transfer_props);

    properties_group_.getDataSet("efferent_surface_x")
        .select({offset}, {count})
        .write(data.preSurfacePosition_1, transfer_props);
    properties_group_.getDataSet("efferent_surface_y")
        .select({offset}, {count})
        .write(data.preSurfacePosition_2, transfer_props);
    properties_group_.getDataSet("efferent_surface_z")
        .select({offset}, {count})
        .write(data.preSurfacePosition_3, transfer_props);

    properties_group_.getDataSet("afferent_center_x")
        .select({offset}, {count})
        .write(data.postCenterPosition_1, transfer_props);
    properties_group_.getDataSet("afferent_center_y")
        .select({offset}, {count})
        .write(data.postCenterPosition_2, transfer_props);
    properties_group_.getDataSet("afferent_center_z")
        .select({offset}, {count})
        .write(data.postCenterPosition_3, transfer_props);

    properties_group_.getDataSet("spine_length")
        .select({offset}, {count})
        .write(data.spineLength, transfer_props);

    properties_group_.getDataSet("efferent_section_type")
        .select({offset}, {count})
        .write(data.preSectionType, transfer_props);
    properties_group_.getDataSet("afferent_section_type")
        .select({offset}, {count})
        .write(data.postSectionType, transfer_props);

    properties_group_.getDataSet("efferent_center_x")
        .select({offset}, {count})
        .write(data.preCenterPosition_1, transfer_props);
    properties_group_.getDataSet("efferent_center_y")
        .select({offset}, {count})
        .write(data.preCenterPosition_2, transfer_props);
    properties_group_.getDataSet("efferent_center_z")
        .select({offset}, {count})
        .write(data.preCenterPosition_3, transfer_props);

    properties_group_.getDataSet("afferent_surface_x")
        .select({offset}, {count})
        .write(data.postSurfacePosition_1, transfer_props);
    properties_group_.getDataSet("afferent_surface_y")
        .select({offset}, {count})
        .write(data.postSurfacePosition_2, transfer_props);
    properties_group_.getDataSet("afferent_surface_z")
        .select({offset}, {count})
        .write(data.postSurfacePosition_3, transfer_props);

    return {count, offset, 0, 0};
}

HighFive::File SonataWriter::open_file(const std::string& fileName) {
    HighFive::FileAccessProps fapl;
    fapl.add(HighFive::MPIOFileAccess{mpi_comm_, MPI_INFO_NULL});
    fapl.add(HighFive::MPIOCollectiveMetadata{});
    return HighFive::File(fileName,
                          HighFive::File::Truncate | HighFive::File::Create |
                              HighFive::File::ReadWrite,
                          fapl);
}

DataVectors SonataWriter::transform_to_vectors(const std::vector<DataEntry>& appositions) {
    DataVectors result(appositions.size());
    for (const auto& entry: appositions) {
        result.preCell.push_back(entry.preCell);
        result.preSection.push_back(entry.preSection);
        result.preSegment.push_back(entry.preSegment);
        result.postCell.push_back(entry.postCell);
        result.postSection.push_back(entry.postSection);
        result.postSegment.push_back(entry.postSegment);
        result.sectionOrder.push_back(entry.sectionOrder);
        result.distance_to_soma.push_back(entry.distance_to_soma);
        result.efferent_segment_offset.push_back(entry.efferent_segment_offset);
        result.afferent_segment_offset.push_back(entry.afferent_segment_offset);
        result.preSectionFraction.push_back(entry.preSectionFraction);
        result.postSectionFraction.push_back(entry.postSectionFraction);
        result.preSurfacePosition_1.push_back(entry.preSurfacePosition[0]);
        result.preSurfacePosition_2.push_back(entry.preSurfacePosition[1]);
        result.preSurfacePosition_3.push_back(entry.preSurfacePosition[2]);
        result.postCenterPosition_1.push_back(entry.postCenterPosition[0]);
        result.postCenterPosition_2.push_back(entry.postCenterPosition[1]);
        result.postCenterPosition_3.push_back(entry.postCenterPosition[2]);
        result.spineLength.push_back(entry.spineLength);
        uint32_t pre_type = ((entry.sectionType >> DataEntry::SECTION_SHIFT) & DataEntry::SECTION_MASK) + SECTION_OFFSET;
        uint32_t post_type = (entry.sectionType & DataEntry::SECTION_MASK) + SECTION_OFFSET;
        result.preSectionType.push_back(pre_type);
        result.postSectionType.push_back(post_type);
        result.preCenterPosition_1.push_back(entry.preCenterPosition[0]);
        result.preCenterPosition_2.push_back(entry.preCenterPosition[1]);
        result.preCenterPosition_3.push_back(entry.preCenterPosition[2]);
        result.postSurfacePosition_1.push_back(entry.postSurfacePosition[0]);
        result.postSurfacePosition_2.push_back(entry.postSurfacePosition[1]);
        result.postSurfacePosition_3.push_back(entry.postSurfacePosition[2]);
    }
    return result;
}

DataVectors::DataVectors(size_t size) {
    preCell.reserve(size);
    preSection.reserve(size);
    preSegment.reserve(size);
    postCell.reserve(size);
    postSection.reserve(size);
    postSegment.reserve(size);
    sectionOrder.reserve(size);
    distance_to_soma.reserve(size);
    efferent_segment_offset.reserve(size);
    afferent_segment_offset.reserve(size);
    preSectionFraction.reserve(size);
    postSectionFraction.reserve(size);
    preSurfacePosition_1.reserve(size);
    preSurfacePosition_2.reserve(size);
    preSurfacePosition_3.reserve(size);
    postCenterPosition_1.reserve(size);
    postCenterPosition_2.reserve(size);
    postCenterPosition_3.reserve(size);
    spineLength.reserve(size);
    preSectionType.reserve(size);
    postSectionType.reserve(size);
    preCenterPosition_1.reserve(size);
    preCenterPosition_2.reserve(size);
    preCenterPosition_3.reserve(size);
    postSurfacePosition_1.reserve(size);
    postSurfacePosition_2.reserve(size);
    postSurfacePosition_3.reserve(size);
}

}  // namespace io

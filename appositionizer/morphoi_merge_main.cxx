// # Design Principles
// The idea is to take morphologies stored in individual files and merge them
// into a single file.
//
// One restriction is that HDF5 requires that metadata must be created and
// modified collectively. Alternatively, one can have one MPI rank create the
// metadata of the file, i.e. it's structure. Then all ranks (re-)open the file.
//
// Hence, we need to be able to send the metadata between MPI ranks. This is
// much easier if the metadata is trivially copiable, e.g. why `dims` is
// statically sized and not a `std::vector`. Moreover, this is why there's one
// `MetaData` struct that aggregates the metadata for all possible datasets.
// Whether a dataset is actually present to be checked via `is_present`.
//
// Finally, there's multiple morphology formats. They seem to (fortunately)
// simply add more datasets to each morphology.
//
// If one wants to apply some function to the metadata of each dataset, then one
// should use `apply_all` or `apply_active` to avoid writing it out manually.
//
// The application has three phases:
//   1. Each MPI rank reads its share of the metadata.
//   2. Rank 0 creates the file structure and preallocates all datasets.
//   3. Each MPI rank fills the file with its share of the morphologies.
//
// There's two container formats:
//   - "merged" which simply copies the contents of each individual HDF5
//   morphology file into a group of the container. The name of the group within
//   the container is the morphology name.
//
//   - "unified" which concatenates the datasets across groups. Then it creates
//   the same group structure as in the merged format. However, the datasets are
//   virtual datasets into subsets of the unified datasets. Furthermore, it
//   stores the offset and size into the unified datasets as attributes. This
//   enables bypassing the VDS when loading directly from the unified datasets.

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <bbp/sonata/nodes.h>
#include <hdf5.h>
#include <highfive/H5File.hpp>

namespace mpi {
int rank(MPI_Comm comm) {
    int mpi_rank = -1;
    MPI_Comm_rank(comm, &mpi_rank);

    return mpi_rank;
}

int size(MPI_Comm comm) {
    int mpi_size = -1;
    MPI_Comm_size(comm, &mpi_size);

    return mpi_size;
}
}  // namespace mpi

// The floating point type for anything morphology related.
using float_t = float;

template <size_t n_dims>
void allocate_dataset(HighFive::Group& group,
                      const std::string& tag,
                      const std::array<size_t, n_dims>& dims,
                      const HighFive::DataType& data_type) {
    auto data_space = HighFive::DataSpace(dims);

    // HDF5 can delay growing the file. However, for our purpose we need
    // it to request the disk space now.
    auto create_props = HighFive::DataSetCreateProps{};
    create_props.add(HighFive::AllocationTime(H5D_ALLOC_TIME_EARLY));

    group.createDataSet(tag, data_space, data_type, create_props);
}

template <class Metadata>
void allocate_merged_dataset(HighFive::Group& group,
                             const std::string& tag,
                             const Metadata& metadata,
                             const std::vector<size_t>& offsets) {
    auto dims = metadata.dims;
    dims[0] = offsets.back();

    auto dataset_group = group.createGroup(tag);
    dataset_group.createDataSet("offsets", offsets);
    allocate_dataset(dataset_group, "values", dims, metadata.data_type());
}

template <class T>
void fill_dataset(HighFive::Group& out,
                  const HighFive::File& h5_in,
                  const std::string& dataset_name,
                  std::vector<std::byte>& buffer) {
    auto datatype = HighFive::create_datatype<T>();
    h5_in.getDataSet(dataset_name).read((T*) buffer.data(), datatype);
    auto dataset = out.getDataSet(dataset_name);
    dataset.write_raw((T*) buffer.data(), datatype);
}

template <class T>
void fill_dataset(HighFive::Group& out,
                  const HighFive::File& h5_in,
                  const std::string& dataset_name,
                  std::vector<std::byte>& buffer,
                  const HighFive::RegularHyperSlab& hyperslab) {
    auto datatype = HighFive::create_datatype<T>();

    h5_in.getDataSet(dataset_name).read((T*) buffer.data(), datatype);

    auto dataset = out.getDataSet(".raw/" + dataset_name + "/values");
    auto selection = dataset.select(HighFive::HyperSlab(hyperslab));

    selection.write_raw((T*) buffer.data(), datatype);
}

template <size_t n_dims>
struct DatasetMetadata {
    bool is_present = false;
    std::array<size_t, n_dims> dims;
    size_t size = 0ul;  ///< total size in bytes.

    DatasetMetadata() = default;

    template <class H5>
    DatasetMetadata(const H5& h5, const std::string& tag, size_t scalar_size) {
        is_present = h5.exist(tag);

        if (is_present) {
            auto dataset = h5.getDataSet(tag);
            auto dataspace = dataset.getSpace();

            auto h5_dims = dataspace.getDimensions();
            std::copy(h5_dims.begin(), h5_dims.end(), dims.begin());

            auto n_elements = std::reduce(h5_dims.begin(),
                                          h5_dims.end(),
                                          size_t(1),
                                          [](auto a, auto b) { return a * b; });
            size = n_elements * scalar_size;
        }
    }
};

template <class Derived, size_t n_dims>
struct DatasetMetadataCRTP: public DatasetMetadata<n_dims> {
    DatasetMetadataCRTP() = default;

    template <class H5>
    explicit DatasetMetadataCRTP(const H5& h5)
        : DatasetMetadata<n_dims>(h5, Derived::path(), sizeof(typename Derived::scalar_type)) {}

    template <class H5>
    void allocate_dataset_merged(H5& h5) const {
        if (this->is_present) {
            ::allocate_dataset(h5, Derived::path(), this->dims, data_type());
        }
    }

    template <class H5>
    void fill_dataset_merged(HighFive::Group& out,
                             const H5& in,
                             std::vector<std::byte>& buffer) const {
        if (this->is_present) {
            if (buffer.size() < this->size) {
                buffer.resize(this->size);
            }

            ::fill_dataset<typename Derived::scalar_type>(out, in, Derived::path(), buffer);
        }
    }

    template <class H5>
    void fill_dataset_unified(HighFive::Group& out,
                              const H5& in,
                              std::vector<std::byte>& buffer,
                              size_t global_offset) const {
        if (this->is_present) {
            if (buffer.size() < this->size) {
                buffer.resize(this->size);
            }

            size_t rank = this->dims.size();
            std::vector<size_t> offsets(rank, 0ul);
            offsets[0] = global_offset;

            std::vector<size_t> stride(rank, 1ul);
            std::vector<size_t> count(rank, 1ul);
            std::vector<size_t> block(this->dims.begin(), this->dims.end());

            auto hyperslab = HighFive::RegularHyperSlab(offsets, count, stride, block);
            ::fill_dataset<typename Derived::scalar_type>(
                out, in, Derived::path(), buffer, hyperslab);
        }
    }

    static HighFive::DataType data_type() {
        return HighFive::create_datatype<typename Derived::scalar_type>();
    }
};

#define MORPHOKIT_REGISTER_ARRAY(ClassName, NDims, ScalarType, Name) \
    struct ClassName: public DatasetMetadataCRTP<ClassName, NDims> { \
      private:                                                       \
        using super = DatasetMetadataCRTP<ClassName, NDims>;         \
                                                                     \
      public:                                                        \
        using super::super;                                          \
                                                                     \
        using scalar_type = ScalarType;                              \
                                                                     \
        static std::string path() {                                  \
            return group_name + Name;                                \
        }                                                            \
    }

namespace base {
const std::string group_name = "";

MORPHOKIT_REGISTER_ARRAY(PointsMetadata, 2, float_t, "points");
MORPHOKIT_REGISTER_ARRAY(StructureMetadata, 2, std::int32_t, "structure");
MORPHOKIT_REGISTER_ARRAY(PerimetersMetadata, 1, float_t, "perimeters");
}  // namespace base

namespace mitochondria {
const std::string group_name = "organelles/mitochondria/";

MORPHOKIT_REGISTER_ARRAY(PointsMetadata, 2, float_t, "points");
MORPHOKIT_REGISTER_ARRAY(StructureMetadata, 2, std::int32_t, "structure");
}  // namespace mitochondria

namespace endoplasmatic_reticulum {
const std::string group_name = "organelles/endoplasmic_reticulum/";

MORPHOKIT_REGISTER_ARRAY(SectionIndexMetadata, 1, std::uint32_t, "section_index");
MORPHOKIT_REGISTER_ARRAY(VolumeMetadata, 1, float_t, "volume");
MORPHOKIT_REGISTER_ARRAY(SurfaceAreaMetadata, 1, float_t, "surface_area");
MORPHOKIT_REGISTER_ARRAY(FilamentCountMetadata, 1, std::uint32_t, "filament_count");
}  // namespace endoplasmatic_reticulum

namespace dendritic_spine {
const std::string group_name = "organelles/postsynaptic_density/";

MORPHOKIT_REGISTER_ARRAY(SectionIdsMetadata, 1, std::int32_t, "section_id");
MORPHOKIT_REGISTER_ARRAY(SegmentIdsMetadata, 1, std::int32_t, "segment_id");
MORPHOKIT_REGISTER_ARRAY(OffsetMetadata, 1, float_t, "offset");
}  // namespace dendritic_spine

#undef MORPHOKIT_REGISTER_ARRAY

using cellfamily_t = std::array<unsigned int, 1>;
using morphology_version_t = std::array<unsigned int, 2>;

// The group is called metadata.
struct MetadataMetadata {
    bool is_present;
    cellfamily_t cell_family;
    morphology_version_t version;

    MetadataMetadata() = default;

    template <class H5>
    explicit MetadataMetadata(const H5& h5) {
        is_present = h5.exist("metadata");

        if (is_present) {
            auto group = h5.getGroup("metadata");

            group.getAttribute("cell_family").read(cell_family);
            group.getAttribute("version").read(version);

            if (version[0] != 1) {
                throw std::runtime_error("Only morphology versions 1.x are supported.");
            }

            if (version[1] > 3) {
                throw std::runtime_error("Only morphology versions 1.x with x <= 3 are supported.");
            }
        }
    }

    template <class H5>
    void write(H5& h5) const {
        if (is_present) {
            auto group = h5.createGroup("metadata");
            group.createAttribute("cell_family", cell_family);
            group.createAttribute("version", version);
        }
    }
};

#define MORPHOKIT_REGISTER_APPLIED_METADATA(md) \
    if (md.path() == path) {                    \
        f(md);                                  \
        return;                                 \
    }

struct MorphologyMetadata {
    MetadataMetadata metadata;

    base::PointsMetadata points;
    base::StructureMetadata structure;
    base::PerimetersMetadata perimeters;

    mitochondria::PointsMetadata mito_points;
    mitochondria::StructureMetadata mito_structure;

    endoplasmatic_reticulum::SectionIndexMetadata er_section_index;
    endoplasmatic_reticulum::SurfaceAreaMetadata er_surface_area;
    endoplasmatic_reticulum::VolumeMetadata er_volume;
    endoplasmatic_reticulum::FilamentCountMetadata er_filament_counts;

    dendritic_spine::SectionIdsMetadata ds_section_ids;
    dendritic_spine::SegmentIdsMetadata ds_segment_ids;
    dendritic_spine::OffsetMetadata ds_offset;

    MorphologyMetadata() = default;

    template <class H5>
    explicit MorphologyMetadata(const H5& in)
        : metadata(in)
        , points(in)
        , structure(in)
        , perimeters(in)
        , mito_points(in)
        , mito_structure(in)
        , er_section_index(in)
        , er_surface_area(in)
        , er_volume(in)
        , er_filament_counts(in)
        , ds_section_ids(in)
        , ds_segment_ids(in)
        , ds_offset(in) {}

    void allocate_merged(HighFive::Group& group) const {
        apply_all([&group](const auto& m) { m.allocate_dataset_merged(group); });
    }

    /** \brief Apply `f` to the metadata for the dataset `path`.
     *
     *  The function `f` should have the following signature:
     *
     *     template<class MC>
     *     void f(const MD& md)
     */
    template <class F>
    void apply(const std::string& path, const F& f) const {
        MORPHOKIT_REGISTER_APPLIED_METADATA(points);
        MORPHOKIT_REGISTER_APPLIED_METADATA(structure);
        MORPHOKIT_REGISTER_APPLIED_METADATA(perimeters);
        MORPHOKIT_REGISTER_APPLIED_METADATA(mito_points);
        MORPHOKIT_REGISTER_APPLIED_METADATA(mito_structure);
        MORPHOKIT_REGISTER_APPLIED_METADATA(er_section_index);
        MORPHOKIT_REGISTER_APPLIED_METADATA(er_surface_area);
        MORPHOKIT_REGISTER_APPLIED_METADATA(er_volume);
        MORPHOKIT_REGISTER_APPLIED_METADATA(er_filament_counts);
        MORPHOKIT_REGISTER_APPLIED_METADATA(ds_section_ids);
        MORPHOKIT_REGISTER_APPLIED_METADATA(ds_segment_ids);
        MORPHOKIT_REGISTER_APPLIED_METADATA(ds_offset);

        throw std::invalid_argument("No metadata for: " + path);
    }

    /** \brief Apply `f` to metadata of each dataset, present or not. */
    template <class F>
    void apply_all(const F& f) const {
        for (const auto& path: paths()) {
            apply(path, f);
        }
    }

    /** \brief Apply `f` to metadata of all present datasets. */
    template <class F>
    void apply_present(const F& f) const {
        for (const auto& path: paths()) {
            apply(path, [&f](const auto& md) {
                if (md.is_present) {
                    f(md);
                }
            });
        }
    }

    std::vector<std::string> paths() const {
        return {points.path(),
                structure.path(),
                perimeters.path(),
                mito_points.path(),
                mito_structure.path(),
                er_section_index.path(),
                er_surface_area.path(),
                er_volume.path(),
                er_filament_counts.path(),
                ds_section_ids.path(),
                ds_segment_ids.path(),
                ds_offset.path()};
    }

    template <class H5>
    void fill_merged(HighFive::Group& out, const H5& in, std::vector<std::byte>& buffer) const {
        apply_present([&](const auto& m) { m.fill_dataset_merged(out, in, buffer); });
    }

    template <class H5>
    void fill_unified(HighFive::Group& out,
                      const H5& in,
                      std::vector<std::byte>& buffer,
                      const std::map<std::string, size_t>& offsets) const {
        apply_present(
            [&](const auto& m) { m.fill_dataset_unified(out, in, buffer, offsets.at(m.path())); });
    }
};

using MorphologiesMetadata = std::vector<MorphologyMetadata>;

// VERSION 1: merged: without minifying names
// VERSION 2: merged: with minifying names
// VERSION 3: merged: page allocated, no minified names.
// VERSION 4: unified: page allocated, merged datasets across groups.
enum class ContainerFormat { merged = 3, unified = 4 };

struct PagedFileParams {
    size_t blocksize;
    size_t pagesize;
    ContainerFormat version;
};

HighFive::FileAccessProps make_default_fapl() {
    HighFive::FileAccessProps file_access_props;
    file_access_props.add(HighFive::FileVersionBounds(H5F_LIBVER_V110, H5F_LIBVER_V110));

    return file_access_props;
}

HighFive::FileAccessProps make_fapl(PagedFileParams params) {
    auto file_access_props = make_default_fapl();

    if (params.blocksize > 0) {
        file_access_props.add(HighFive::MetadataBlockSize(params.blocksize));
    }

    return file_access_props;
}

HighFive::FileCreateProps make_fcpl(PagedFileParams params) {
    HighFive::FileCreateProps file_create_props;
    if (params.pagesize > 0) {
        file_create_props.add(HighFive::FileSpaceStrategy(H5F_FSPACE_STRATEGY_PAGE, false, 0));
        file_create_props.add(HighFive::FileSpacePageSize(params.pagesize));
    }

    return file_create_props;
}

HighFive::File create_file(const std::string& filename, PagedFileParams file_params) {
    auto file_access_props = make_fapl(file_params);
    auto file_create_props = make_fcpl(file_params);
    auto file =
        HighFive::File(filename, HighFive::File::Truncate, file_create_props, file_access_props);
    file.createAttribute("version", int(file_params.version));

    return file;
}

HighFive::File open_file(const std::string& filename, PagedFileParams file_params) {
    auto file_access_props = make_fapl(file_params);
    file_access_props.add(HighFive::MPIOFileAccess(MPI_COMM_WORLD, MPI_INFO_NULL));

    return HighFive::File(filename, HighFive::File::ReadWrite, file_access_props);
}

std::string individual_filename(const std::string& dirname, const std::string& morphology) {
    return std::filesystem::path(dirname) / (morphology + ".h5");
}

std::array<size_t, 2> balanced_chunks(size_t length, size_t n_chunks, size_t k_chunk) {
    size_t chunk_size = length / n_chunks;
    size_t n_large_chunks = length % n_chunks;

    size_t low = k_chunk * chunk_size + std::min(k_chunk, n_large_chunks);
    size_t high = low + chunk_size + (k_chunk < n_large_chunks);

    return {std::min(low, length), std::min(high, length)};
};

template <class TimePoint>
double elapsed_seconds_since(const TimePoint& start) {
    auto duration = std::chrono::steady_clock::now() - start;
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    return elapsed.count() * 1e-6;
}

class ProgressReport {
  public:
    ProgressReport() {
        comm_rank = mpi::rank(MPI_COMM_WORLD);
        reset();
    }

    void reset() {
        time_point = std::chrono::steady_clock::now();
    }

    void write(const std::string& phase, size_t n_processed, size_t n_morphologies) {
        if (comm_rank == 0 && n_processed >= next_output) {
            auto elapsed = elapsed_seconds_since(time_point);

            auto n_since_last = n_processed - last_output;

            auto rate = n_since_last / elapsed;
            auto remainder = (n_morphologies - n_processed) / rate;
            std::cout << phase << " " << n_processed << " morphologies, at " << rate
                      << " morphologies/s, " << remainder << "s remaining\n";

            time_point = std::chrono::steady_clock::now();

            next_output = ((n_processed / output_interval) + 1) * output_interval;
            last_output = n_processed;
        }
    }

  private:
    size_t output_interval = 10000ul;
    size_t last_output = 0ul;
    size_t next_output = output_interval;

    int comm_rank;

    std::chrono::steady_clock::time_point time_point;
};

template <class TimePoint>
void write_stage_progress(const std::string& phase, size_t n_morphologies, const TimePoint& start) {
    auto mpi_rank = mpi::rank(MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        auto elapsed = elapsed_seconds_since(start);

        auto rate = n_morphologies / elapsed;
        std::cout << phase << n_morphologies << " morphologies, at " << rate
                  << " morphologies/s.\n";
    }
}

std::vector<std::string> get_morphology_names_from_nodesfile(const std::string& nodefile) {
    bbp::sonata::NodeStorage store(nodefile);
    const auto& population_names = store.populationNames();

    auto morphology_names = std::vector<std::string>{};
    for (const auto& pop_name: population_names) {
        const auto& pop = store.openPopulation(pop_name);
        const auto& names = pop->getAttribute<std::string>("morphology", pop->selectAll());
        morphology_names.insert(morphology_names.end(), names.begin(), names.end());
    }

    std::sort(morphology_names.begin(), morphology_names.end());
    auto new_end = std::unique(morphology_names.begin(), morphology_names.end());
    morphology_names.erase(new_end, morphology_names.end());

    return morphology_names;
}

MorphologiesMetadata load_metadata(const std::vector<std::string>& morphologies,
                                   const std::array<size_t, 2>& range,
                                   const std::string& dirname) {
    auto comm_size = mpi::size(MPI_COMM_WORLD);
    auto progress_report = ProgressReport();

    MorphologiesMetadata metadata;
    metadata.reserve(range[1] - range[0]);

    for (size_t i = range[0]; i < range[1]; ++i) {
        HighFive::File in(individual_filename(dirname, morphologies[i]));

        metadata.push_back(MorphologyMetadata(in));

        progress_report.write("Gathered metadata for", i * comm_size, morphologies.size());
    }

    return metadata;
}

MorphologiesMetadata gather_metadata(const MorphologiesMetadata& local_metadata,
                                     size_t n_morphologies,
                                     MPI_Comm comm) {
    int comm_rank = mpi::rank(comm);
    int comm_size = mpi::size(comm);

    MPI_Datatype mpi_metadata_type;
    MPI_Type_contiguous(sizeof(MorphologyMetadata), MPI_BYTE, &mpi_metadata_type);
    MPI_Type_commit(&mpi_metadata_type);

    if (comm_rank == 0) {
        auto offsets = std::vector<int>(comm_size);
        auto counts = std::vector<int>(comm_size);

        for (size_t i = 0; i < comm_size; ++i) {
            auto local_range = balanced_chunks(n_morphologies, comm_size, i);
            offsets[i] = local_range[0];
            counts[i] = local_range[1] - local_range[0];
        }

        auto global_metadata = MorphologiesMetadata(n_morphologies);

        MPI_Gatherv(local_metadata.data(),
                    local_metadata.size(),
                    mpi_metadata_type,
                    global_metadata.data(),
                    counts.data(),
                    offsets.data(),
                    mpi_metadata_type,
                    0,
                    comm);

        return global_metadata;

    } else {
        MPI_Gatherv(local_metadata.data(),
                    local_metadata.size(),
                    mpi_metadata_type,
                    nullptr,
                    nullptr,
                    nullptr,
                    mpi_metadata_type,
                    0,
                    comm);

        return MorphologiesMetadata{};
    }

    MPI_Type_free(&mpi_metadata_type);
}

template <class Metadata>
void create_virtual_dataset(const std::string& filename,
                            HighFive::Group& group,
                            const std::string& tag,
                            size_t global_offset,
                            size_t global_size,
                            const Metadata& md) {
    int rank = static_cast<int>(md.dims.size());
    std::vector<hsize_t> dims(md.dims.begin(), md.dims.end());

    std::vector<hsize_t> global_dims = dims;
    global_dims[0] = global_size;

    std::vector<hsize_t> h5_offset(rank);
    h5_offset[0] = global_offset;

    std::vector<hsize_t> h5_stride(rank, 1ul);
    std::vector<hsize_t> h5_count(rank, 1ul);
    std::vector<hsize_t> h5_block = dims;

    // Open VDS
    hid_t dst_space_id = H5Screate_simple(rank, dims.data(), dims.data());
    hid_t src_space_id = H5Screate_simple(rank, global_dims.data(), global_dims.data());
    H5Sselect_hyperslab(src_space_id,
                        H5S_SELECT_SET,
                        h5_offset.data(),
                        h5_stride.data(),
                        h5_count.data(),
                        h5_block.data());

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    std::string src_tag = ".raw/" + tag + "/values";
    H5Pset_virtual(dcpl_id, dst_space_id, filename.c_str(), src_tag.c_str(), src_space_id);

    auto lcpl = HighFive::LinkCreateProps();
    lcpl.add(HighFive::CreateIntermediateGroup(true));

    auto datatype = md.data_type();
    hid_t dset = H5Dcreate(group.getId(),
                           tag.c_str(),
                           datatype.getId(),
                           dst_space_id,
                           lcpl.getId(),
                           dcpl_id,
                           H5P_DEFAULT);

    // Close VDS
    H5Dclose(dset);
    H5Pclose(dcpl_id);
    H5Sclose(src_space_id);
    H5Sclose(dst_space_id);
}


template <class FileParams>
void allocate_merged_file(const std::string& filename,
                          const std::string& dirname,
                          const MorphologiesMetadata& metadata,
                          const std::vector<std::string> morphologies,
                          const FileParams& file_params) {
    auto n_morphologies = morphologies.size();

    auto file = create_file(filename, file_params);

    auto progress_report = ProgressReport();
    for (size_t i = 0; i < morphologies.size(); ++i) {
        const auto& morph = morphologies[i];
        auto group = file.createGroup(morph, true);

        const auto& md = metadata[i];
        md.allocate_merged(group);

        md.metadata.write(group);
        progress_report.write("Allocated", i, morphologies.size());
    }
}

template <class FileParams>
void allocate_unified_file(const std::string& filename,
                           const std::string& dirname,
                           const MorphologiesMetadata& metadata,
                           const std::vector<std::string> morphologies,
                           const FileParams& file_params) {
    auto n_morphologies = morphologies.size();

    auto file = create_file(filename, file_params);

    std::set<std::string> present_dataset_names;
    auto merge_dataset_names = [&present_dataset_names](const auto& m) {
        present_dataset_names.insert(m.path());
    };

    for (size_t i = 0; i < morphologies.size(); ++i) {
        const auto& md = metadata[i];
        md.apply_present(merge_dataset_names);
    }

    auto offsets = std::map<std::string, std::vector<size_t>>{};
    for (const auto& name: present_dataset_names) {
        offsets[name] = std::vector<size_t>(n_morphologies + 1);
    }

    auto progress_report = ProgressReport();
    for (size_t i = 0; i < morphologies.size(); ++i) {
        const auto& md = metadata[i];
        for (const auto& name: present_dataset_names) {
            md.apply(name, [i, &name, &offsets](const auto& m) {
                size_t size = m.is_present ? m.dims[0] : size_t(0);
                offsets[name][i + 1] = offsets[name][i] + size;
            });
        }
    }

    for (size_t i = 0; i < morphologies.size(); ++i) {
        const auto& morph_name = morphologies[i];
        auto group = file.createGroup(morph_name);
        group.createAttribute("morphology_id", i);

        const auto& md = metadata[i];

        for (const auto& name: present_dataset_names) {
            md.apply(name, [i, &name, &filename, &group, &offsets](const auto& m) {
                if (m.is_present) {
                    group.createAttribute(name + "_begin", offsets[name][i]);
                    group.createAttribute(name + "_end", offsets[name][i + 1]);
                    create_virtual_dataset(
                        filename, group, name, offsets[name][i], offsets[name].back(), m);
                }
            });
        }

        md.metadata.write(group);

        progress_report.write("Allocated", i, morphologies.size());
    }

    auto raw_group = file.createGroup(".raw");
    for (const auto& name: present_dataset_names) {
        for (const auto& md: metadata) {
            bool is_present = false;

            md.apply(name, [&is_present, &raw_group, &name, &offsets](const auto& m) {
                if (m.is_present) {
                    is_present = true;
                    allocate_merged_dataset(raw_group, name, m, offsets[name]);
                }
            });

            if (is_present) {
                break;
            }
        }
    }
}

template <class FileParams>
void allocate_file(const std::string& filename,
                   const std::string& dirname,
                   const MorphologiesMetadata& metadata,
                   const std::vector<std::string> morphologies,
                   const FileParams& file_params) {
    if (file_params.version == ContainerFormat::merged) {
        allocate_merged_file(filename, dirname, metadata, morphologies, file_params);
    } else if (file_params.version == ContainerFormat::unified) {
        allocate_unified_file(filename, dirname, metadata, morphologies, file_params);
    } else {
        throw std::runtime_error("Unsupported ContainerFormat: " +
                                 std::to_string(int(file_params.version)));
    }
}

template <class FileParams>
void fill_merged_file(const std::string& filename,
                      const FileParams& file_params,
                      const MorphologiesMetadata& metadata,
                      const std::vector<std::string>& morphologies,
                      const std::array<size_t, 2>& morphologies_range,
                      const std::string& dirname) {
    int comm_size = mpi::size(MPI_COMM_WORLD);

    auto h5_out = open_file(filename, file_params);
    auto buffer = std::vector<std::byte>(1ul << 20);

    auto progress_report = ProgressReport();
    auto n_local_morphologies = morphologies_range[1] - morphologies_range[0];
    for (size_t k = 0; k < n_local_morphologies; ++k) {
        const auto& morph_name = morphologies[morphologies_range[0] + k];
        auto in_name = individual_filename(dirname, morph_name);
        auto h5_in = HighFive::File(in_name, HighFive::File::ReadOnly);

        auto h5_group = h5_out.getGroup(morph_name);
        metadata[k].fill_merged(h5_group, h5_in, buffer);

        progress_report.write("Filled ", k * comm_size, morphologies.size());
    }
}

template <class FileParams>
void fill_unified_file(const std::string& filename,
                       const FileParams& file_params,
                       const MorphologiesMetadata& metadata,
                       const std::vector<std::string>& morphologies,
                       const std::array<size_t, 2>& morphologies_range,
                       const std::string& dirname) {
    int comm_size = mpi::size(MPI_COMM_WORLD);

    auto h5_out = open_file(filename, file_params);
    auto buffer = std::vector<std::byte>(1ul << 20);

    auto offsets = std::map<std::string, std::vector<size_t>>{};
    for (const auto& path: metadata[0].paths()) {
        auto raw_group = h5_out.getGroup(".raw");
        if (raw_group.exist(path)) {
            offsets[path] =
                raw_group.getDataSet(path + "/offsets").template read<std::vector<size_t>>();
        }
    }

    auto progress_report = ProgressReport();
    auto n_local_morphologies = morphologies_range[1] - morphologies_range[0];
    for (size_t k = 0; k < n_local_morphologies; ++k) {
        const auto& morph_name = morphologies[morphologies_range[0] + k];
        auto in_name = individual_filename(dirname, morph_name);
        auto h5_in = HighFive::File(in_name, HighFive::File::ReadOnly);

        size_t morph_id =
            h5_out.getGroup(morph_name).getAttribute("morphology_id").template read<size_t>();

        std::map<std::string, size_t> morph_offsets;
        for (const auto& [key, value]: offsets) {
            morph_offsets[key] = value[morph_id];
        }

        auto h5_grp = h5_out.getGroup("/");
        metadata[k].fill_unified(h5_grp, h5_in, buffer, morph_offsets);
        progress_report.write("Filled ", k * comm_size, morphologies.size());
    }
}

template <class FileParams>
void fill_file(const std::string& filename,
               const FileParams& file_params,
               const MorphologiesMetadata& metadata,
               const std::vector<std::string>& morphologies,
               const std::array<size_t, 2>& morphologies_range,
               const std::string& dirname) {
    if (file_params.version == ContainerFormat::merged) {
        fill_merged_file(
            filename, file_params, metadata, morphologies, morphologies_range, dirname);
    } else if (file_params.version == ContainerFormat::unified) {
        fill_unified_file(
            filename, file_params, metadata, morphologies, morphologies_range, dirname);
    } else {
        throw std::runtime_error("Unsupported ContainerFormat: " +
                                 std::to_string(int(file_params.version)));
    }
}

template <class FileParams>
void merge_morphologies(const std::string& outname,
                        const std::vector<std::string>& morphologies,
                        const std::string& dirname,
                        const FileParams& file_params) {
    if (morphologies.empty()) {
        throw std::runtime_error("no morphologies to process");
    }

    auto comm = MPI_COMM_WORLD;
    auto comm_rank = mpi::rank(comm);
    auto comm_size = mpi::size(comm);

    size_t n_morphologies = morphologies.size();
    auto morphologies_range = balanced_chunks(n_morphologies, comm_size, comm_rank);

    auto start = std::chrono::steady_clock::now();

    auto local_metadata = load_metadata(morphologies, morphologies_range, dirname);
    auto global_metadata = gather_metadata(local_metadata, n_morphologies, comm);

    write_stage_progress("Gathered metadata for ", n_morphologies, start);

    if (comm_rank == 0) {
        auto allocate_start = std::chrono::steady_clock::now();
        allocate_file(outname, dirname, global_metadata, morphologies, file_params);

        write_stage_progress("Allocated ", n_morphologies, allocate_start);
    }

    MPI_Barrier(comm);

    auto fill_start = std::chrono::steady_clock::now();
    fill_file(outname, file_params, local_metadata, morphologies, morphologies_range, dirname);

    write_stage_progress("Filled ", n_morphologies, fill_start);
    write_stage_progress("Merged ", n_morphologies, start);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    CLI::App app{
        "Merge morphologies stored as individual files into a single "
        "HDF5 file."};

    std::string outname = "morphologies.h5";
    std::string nodefile;
    std::vector<std::string> morphology_names;
    std::string morphologies_file;
    std::string dirname;
    std::string container_format_name = "merged";

    size_t blocksize_kb = 0;
    size_t pagesize_kb = 4096;

    app.add_option("--format", container_format_name, "Container format: 'unified' or 'merged'.");
    app.add_option("--blocksize", blocksize_kb, "Metadata blocksize in kB");
    app.add_option("--pagesize",
                   pagesize_kb,
                   "Use HDF5 paged allocation with the specified page size, in kB. Default: 4096");
    app.add_option("-o,--output", outname, "The output path of the morphology database")
        ->check(CLI::NonexistentPath);
    app.add_option("-n,--nodes", nodefile, "SONATA node file to read morpholgoy names from");
    app.add_option("--morphologies",
                   morphology_names,
                   "One or more path to morphology names to be merged.");
    app.add_option("--morphologies-file",
                   morphologies_file,
                   "Path to a file of morphology names, one per line.");
    app.add_option("input", dirname, "The input directory containing the morphologies to merge")
        ->check(CLI::ExistingDirectory)
        ->required();

    app.footer(R"(Examples:

For larger circuits (O1, S1, SSCX & NCX) use N=1 or N=2 nodes on
BB5 and run the following command from the directory containing
the circuit:

    srun -n$((N * 72)) --account=ACCOUNT \
        morphokit_merge morphologies/h5/ \
            --morphologies-file MORPHOLOGIES_FILE \
            -o h5_merged/paged.h5

Note: scaling is limited and requesting more MPI tasks might not improve
time-to-solution. You should expect individual phases to run at about
2k - 10k morphologies per second, overall the rate is about 1k morphologies
per second.)");
    CLI11_PARSE(app, argc, argv);

    ContainerFormat container_format;
    if (container_format_name == "unified") {
        container_format = ContainerFormat::unified;
    } else if (container_format_name == "merged") {
        container_format = ContainerFormat::merged;
    } else {
        throw std::runtime_error("Unknown container format: " + container_format_name);
    }

    auto file_params = PagedFileParams{blocksize_kb * 1024, pagesize_kb * 1024, container_format};


    auto deduce_morphology_names = [&]() {
        if (nodefile != "") {
            return get_morphology_names_from_nodesfile(nodefile);
        } else if (!morphology_names.empty()) {
            return morphology_names;
        } else if (morphologies_file != "") {
            auto f = std::ifstream(morphologies_file);
            auto lines = std::vector<std::string>();

            std::string line;
            while (std::getline(f, line)) {
                lines.push_back(line);
            }

            return lines;
        }
        throw std::runtime_error("cannot deduce morphology names");
    };

    merge_morphologies(outname, deduce_morphology_names(), dirname, file_params);

    MPI_Finalize();
    return 0;
}

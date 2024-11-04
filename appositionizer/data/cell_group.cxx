#include "cell_group.h"  // Must be on the first substantial line.

#include <bbp/sonata/config.h>
#include <bbp/sonata/node_sets.h>

#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <unordered_set>

#include "spatial_grid_hash.h"
#include "defines.h"

// Auxiliary types for Libsonata file reading
using Positions = std::vector<std::array<double, 4>>;  // TODO: Move to {x,y,z} only
using Rotations = std::vector<std::array<double, 4>>;
using NodePopulation = bbp::sonata::NodePopulation;
using Selection = bbp::sonata::Selection;
using Range = Selection::Range;
using NodeSets = bbp::sonata::NodeSets;

// Naming convention for Sonata files
constexpr char default_population_name[] = "All";
constexpr char did_morpho[] = "morphology";
constexpr char did_mtypes[] = "mtype";

// Labels for attribute retrieval
#define const_char const char* const
constexpr const_char POSITIONS_LABELS[] = {"x", "y", "z"};
constexpr const_char ROTATIONS_LABELS[] = {"orientation_x",
                                           "orientation_y",
                                           "orientation_z",
                                           "orientation_w"};

namespace {

/// \brief convert a slice ids[start:end] to a Selection container
auto toSelection(const std::vector<int>& ids, size_t start, size_t end) -> Selection {
    Selection::Ranges ranges;
    while (start < end) {
        int offset = ids[start];
        int i = 1;
        while (start + i < end and ids[start + i] == offset + i)
            ++i;
        ranges.push_back({offset, offset + i});
        start += i;
    }
    return Selection(ranges);
}

decltype(auto) openPopulation(const std::string& filename, std::string pop_name) {
    // If no population is provided, attempt to load the default or a single one
    if (pop_name.empty()) {
        const auto& storage = std::make_unique<bbp::sonata::NodeStorage>(filename);
        const auto& populations = storage->populationNames();
        if (populations.count(default_population_name)) {
            pop_name = default_population_name;
        } else if (populations.size() == 1) {
            pop_name = *populations.begin();
        } else if (populations.empty()) {
            mpi::abort("File '{}' does not contain any population", filename);
        } else {
            mpi::abort(
                "Multiple populations found in file '{}'. "
                "Please, select one explicitly.",
                filename);
        }
    }

    return std::make_shared<NodePopulation>(filename, "", pop_name);
}

decltype(auto) getNodesetSelection(const std::shared_ptr<NodePopulation> pop,
                                   const std::string& filename,
                                   std::string nodeset) {
    const auto& nodesets = NodeSets::fromFile(filename);
    return nodesets.materialize(nodeset, *pop);
}

template <class T, size_t N>
T getCoordinates(const std::shared_ptr<NodePopulation> pop,
                 const_char labels[N],
                 const Selection& selection) {
    T coordinates;
    std::vector<std::vector<double>> values(N);

    for (size_t i = 0; i < N; ++i) {
        values[i] = pop->getAttribute<double>(labels[i], selection);
    }

    const auto num_coordinates = values[0].size();
    coordinates.resize(num_coordinates, {0});

    for (size_t j = 0; j < num_coordinates; ++j) {
        for (size_t i = 0; i < N; ++i) {
            coordinates[j][i] = values[i][j];
        }
    }

    return coordinates;
}

}  // namespace

namespace cell {

Group::GlobalMetadata Group::GlobalMetadata::fromConfig(const std::string& config_path,
                                                        std::string population_name,
                                                        const std::string& nodeset_name) {
    const auto cfg = bbp::sonata::CircuitConfig::fromFile(config_path);
    if (population_name.empty()) {
        const auto pops = cfg.listNodePopulations();
        if (pops.size() == 1) {
            population_name = *(pops.begin());
            mpi::rank0::info("Auto-determination of population name: {}", population_name);
        } else {
            mpi::abort("Cannot determine node population");
        }
    }

    const auto& props = cfg.getNodePopulationProperties(population_name);
    // As per SONATA extension documentation, this should default to a directory using SWC
    // files
    std::string morpho_path = props.morphologiesDir;
    auto it = props.alternateMorphologyFormats.find("h5v1");
    // If H5 morphologies are available, prefer those.
    if (it != props.alternateMorphologyFormats.end()) {
        morpho_path = it->second;
    }

    return {
        props.elementsPath,
        population_name,
        cfg.getNodeSetsPath(),
        nodeset_name,
        morpho_path,
    };
}

Group::GlobalMetadata::GlobalMetadata(const std::string& population_path_,
                                      const std::string& population_name_,
                                      const std::string& nodesets_path_,
                                      const std::string& nodeset_name_,
                                      const std::string& morphology_path_)
    : population(openPopulation(population_path_, population_name_))
    , population_path(population_path_)
    , population_name(population_name_)
    , population_id(indexPopulation(population_path_, population_name_))
    , nodesets_path(nodesets_path_)
    , nodeset_name(nodeset_name_)
    , morphology_path(morphology_path_)
    , selection(createSelection()) {}


bbp::sonata::Selection Group::GlobalMetadata::createSelection() const {
    if (!nodesets_path.empty() && !nodeset_name.empty()) {
        return getNodesetSelection(population, nodesets_path, nodeset_name);
    } else {
        return population->selectAll();
    }
}


Group::Group(const Group& o, const std::vector<int>& ids)
    : group_metadata_(o.group_metadata_)
    , morphologies_(o.morphologies_)
    , spine_lengths_(o.spine_lengths_)
    , apposition_space_(o.apposition_space_)
    , comm_(o.comm_)
    , comm_size_(o.comm_size_)
    , comm_rank_(o.comm_rank_) {
    loadMetadata(openPopulation(group_metadata_.population_path,
                                group_metadata_.population->name()),
                 toSelection(ids, 0, ids.size()));
}

Group::Group(const Group& o, const std::vector<Metadata>& meta)
    : group_metadata_(o.group_metadata_)
    , morphologies_(o.morphologies_)
    , spine_lengths_(o.spine_lengths_)
    , apposition_space_(o.apposition_space_)
    , metadata_(meta)
    , comm_(o.comm_)
    , comm_size_(o.comm_size_)
    , comm_rank_(o.comm_rank_) {
    for (size_t n = 0; n < metadata_.size(); ++n) {
        index_[metadata_[n].cellId] = n;
    }
}

decltype(Group::known_populations_) Group::known_populations_{};

size_t Group::indexPopulation(const std::string& filename, const std::string& population) {
    for (size_t i = 0; i < known_populations_.size(); ++i) {
        if (known_populations_[i].first == filename and
            known_populations_[i].second == population) {
            return i;
        }
    }
    known_populations_.emplace_back(filename, population);
    return known_populations_.size() - 1;
}

auto Group::indexSpineLengths(const std::string& filename,
                              const std::string& population,
                              const SpineLengths& lengths) -> std::vector<float> {
    std::vector<float> spines;

    try {
        const auto pop = openPopulation(filename, population);
        const auto mtypes = pop->enumerationValues(did_mtypes);
        for (const auto& mtype: mtypes) {
            const auto match = lengths.find(mtype);
            if (match == lengths.end()) {
                throw std::runtime_error("No spine length present for type " + mtype);
            }
            spines.push_back(match->second);
        }
    } catch (const std::exception& e) {
        mpi::abort("Invalid circuit '{}': {}", filename, e.what());
    }

    return spines;
}

auto Group::createAppositionSpace() -> std::shared_ptr<const AppositionSpace> {
    return createAppositionSpace("nonautaptic", {"axon"}, {"soma", "dendrite"});
}

auto Group::createAppositionSpace(const std::string& space,
                             const std::vector<std::string>& presynaptic,
                             const std::vector<std::string>& postsynaptic)
    -> std::shared_ptr<const AppositionSpace> {
    // Starts the appropriate Apposition spaces
    std::shared_ptr<const AppositionSpace> result;
    if (space == "nonautaptic") {
        result.reset(new NonautapticAppositionSpace(presynaptic, postsynaptic));
    } else if (space == "autaptic") {
        result.reset(new AutapticAppositionSpace(presynaptic, postsynaptic));
    } else {
        throw std::runtime_error("unsupported Apposition space " + space);
    }
    return result;
}


std::vector<int> Group::loadMetadata(const bool check) {
    // Calculate range of each node and load appropriate metadata
    // TODO: Replace ids for uint64_t to avoid conversion after flattening
    const auto& selection_values = group_metadata_.selection.flatten();
    auto ids = std::vector<int>(selection_values.begin(), selection_values.end());

    // Appositionizer does not support node_count > cell_count
    if (check && ids.size() < static_cast<size_t>(comm_size_)) {
        mpi::abort("Fewer cells ({}) than CPUs ({}) in {} / {}. Decrease CPUs",
                   ids.size(),
                   comm_size_,
                   group_metadata_.population_path,
                   group_metadata_.population->name());
    }

    const auto& indices = computeCellIndexRanges(ids.size());
    const auto& range = indices[comm_rank_];
    loadMetadata(group_metadata_.population, toSelection(ids, range.first, range.second));

    return ids;
}

void Group::loadMetadata(const std::shared_ptr<NodePopulation> pop, const Selection& selection) {
    index_.clear();
    metadata_.clear();
    cells_.clear();

    try {
        const auto attrs = pop->attributeNames();
        const bool hasRotations = (attrs.count("orientation_x") + attrs.count("orientation_y") +
                                   attrs.count("orientation_z") + attrs.count("orientation_w")) ==
                                  4;

        if (!hasRotations) {
            mpi::rank0::warn("File '{}' does not contain rotation data",
                             group_metadata_.population_path);
        }


        // read all necessary meta-data using libsonata
        Positions pos = getCoordinates<Positions, 3>(pop, POSITIONS_LABELS, selection);
        Rotations rot(selection.flatSize());
        if (hasRotations) {
            rot = getCoordinates<Rotations, 4>(pop, ROTATIONS_LABELS, selection);
        }

        const auto morphos = pop->getAttribute<std::string>(did_morpho, selection);
        const auto mtypes = pop->getAttribute<size_t>(did_mtypes, selection);
        assert(morphos.size() == selection.flatSize());

        metadata_.reserve(metadata_.size() + selection.flatSize());
        size_t idx = 0;
        for (const int cell_id: selection.flatten()) {
            index_[cell_id] = metadata_.size();
            metadata_.emplace_back(morphos[idx],
                                   pos[idx],
                                   rot[idx],
                                   cell_id,
                                   group_metadata_.population_id,
                                   mtypes[idx],
                                   0,
                                   0,
                                   hasRotations);
            ++idx;
        }
    } catch (const std::exception& e) {
        mpi::abort("Impossible to read population '{}' from '{}': {}",
                   group_metadata_.population->name(),
                   group_metadata_.population_path,
                   e.what());
    }
}

void Group::loadMorphologies(std::function<bool(const SectionType&)> accept) {
    assert(metadata_.size() > 0);
    assert(cells_.size() == 0);

    cells_.reserve(metadata_.size());
    for (const auto& meta: metadata_) {
        unsigned int morphTypeId = meta.morphology;

        // It can assert if we forge metadata ourselves with reset()
        assert(morphTypeId < spine_lengths_.size());

        const auto spine = spine_lengths_[morphTypeId];
        cells_.emplace_back(meta, morphologies_, spine, accept);
    }
}

auto Group::computeCellIndexRanges(size_t count) const -> std::vector<std::pair<size_t, size_t>> {
    std::vector<std::pair<size_t, size_t>> result(comm_size_);
    for (int i = 0; i < comm_size_; i++) {
        if (i == 0) {
            result[i].first = 0;
        } else {
            result[i].first = result[i - 1].second;
        }

        if (i == comm_size_ - 1) {
            result[i].second = count;
        } else {
            result[i].second = round(count / (double) comm_size_ * (i + 1));
        }
    }
    return result;
}

auto gatherStats(MPI_Comm comm, long long count)
    -> std::tuple<long long, long long, long long, long double> {
    long long max;
    long long min;
    long long sum;
    MPI_Allreduce(&count, &max, 1, MPI_LONG_LONG, MPI_MAX, comm);
    MPI_Allreduce(&count, &min, 1, MPI_LONG_LONG, MPI_MIN, comm);
    MPI_Allreduce(&count, &sum, 1, MPI_LONG_LONG, MPI_SUM, comm);
    return std::make_tuple(max, min, sum, ((long double) sum) / mpi::size(comm));
}

void Group::showDistributionInfo() {
    long long sections = 0;
    long long segments = 0;
    for (const auto& cell: cells_) {
        sections += cell.size();
        for (const auto& section: cell)
            segments += section.size();
    }

    const auto cell_stats = gatherStats(comm_, size());
    const auto section_stats = gatherStats(comm_, sections);
    const auto segment_stats = gatherStats(comm_, segments);

    // Output segment distribution information
    mpi::rank0::show("- Cell  Distribution: Max {}, Min {}, Sum {}, Mean {:.2f}",
                     std::get<0>(cell_stats),
                     std::get<1>(cell_stats),
                     std::get<2>(cell_stats),
                     std::get<3>(cell_stats));
    mpi::rank0::show("- Section  Distribution: Max {}, Min {}, Sum {}, Mean {:.2f}",
                     std::get<0>(section_stats),
                     std::get<1>(section_stats),
                     std::get<2>(section_stats),
                     std::get<3>(section_stats));
    mpi::rank0::show("- Segment Distribution: Max {}, Min {}, Sum {}, Mean {:.2f}",
                     std::get<0>(segment_stats),
                     std::get<1>(segment_stats),
                     std::get<2>(segment_stats),
                     std::get<3>(segment_stats));
}

}  // namespace cell

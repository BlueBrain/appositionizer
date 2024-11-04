#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <bbp/sonata/nodes.h>
#include <fmt/format.h>
#include <morphio/collection.h>

#include "data/cell.h"
#include "better_mpi.h"
#include "defines.h"
#include "recipe.h"
#include "data/apposition_space.h"

namespace cell {

/**
 * \brief Container to handle cell metadata and morphologies
 *
 * Two modus operandi are handled: a collective loading of all cell
 * metadata and morphologies, as well as further selective loading of a
 * smaller set of cells.
 *
 * In the former scenario, the container is first created with general
 * settings, and then loadMetadata is used to collectively load the
 * metadata of the brain circuit. Afterwards, loadMorphologies can be
 * employed to fully load the cells.
 *
 * The second usage scenario follows, where a new cell group is created
 * based on the previous "global" instance, and a list of cell indices to
 * load immediately.
 *
 * After either instantiation, the Group container can be iterated over and
 * provides indexed access to the cells loaded.
 */
class Group {
  public:
    using const_iterator = std::vector<Cell>::const_iterator;

    class GlobalMetadata {
      public:
        std::shared_ptr<bbp::sonata::NodePopulation> population;

        std::string population_path;
        std::string population_name;
        const size_t population_id;

        std::string nodesets_path;
        std::string nodeset_name;

        std::string morphology_path;

        bbp::sonata::Selection selection;

        /**
         * \brief Constructor using legacy arguments
         *
         * \param nodesets_path Path to the NodeSets JSON file
         * \param nodeset Nodeset to filter the population
         */
        GlobalMetadata(const std::string& population_path_,
                       const std::string& population_name_,
                       const std::string& nodesets_path_,
                       const std::string& nodeset_name_,
                       const std::string& morphology_path_);

        /**
         * \brief Primary way of obtaining a group of cells
         *
         * \param nodesets_path Path to the NodeSets JSON file
         * \param nodeset Nodeset to filter the population
         */
        static GlobalMetadata fromConfig(const std::string& config_path,
                                         std::string population_name,
                                         const std::string& nodeset_name);

        inline size_t cell_count() const {
            return selection.flatSize();
        }

      private:
        bbp::sonata::Selection createSelection() const;
    };

    /**
     * \brief Instantiate a new cell group from `GlobalMetadata`.
     *
     * The purpose of this constructor is to enable loading some global
     * metadata first, then make a choice, such as the precise communicator
     * to use and only then continue with the construction of the group.
     *
     * \param global_metadata The already loaded metadata.
     * \param morphos The path to the morphology storage
     * \param lengths A mapping of morphology types to spine lengths
     * \param space The Appositionspace to use
     * \param comm The MPI communicator with all participating MPI ranks.
     */
    Group(const GlobalMetadata& metadata,
          const SpineLengths& lengths,
          const std::shared_ptr<const AppositionSpace>& space,
          MPI_Comm comm)
        : group_metadata_(metadata)
        , morphologies_(morphio::Collection(metadata.morphology_path))
        , spine_lengths_(
              indexSpineLengths(metadata.population_path, metadata.population->name(), lengths))
        , apposition_space_(space)
        , comm_(comm)
        , comm_size_(mpi::size(comm))
        , comm_rank_(mpi::rank(comm)){};

    /**
     * \brief Instantiate a new cell group
     *
     * \param circuit The full path to the file containing the cell metadata
     * \param population The population name used when loading the metadata
     * \param morphos The path to the morphology storage
     * \param lengths A mapping of morphology types to spine lengths
     * \param space The Appositionspace to use
     * \param comm The MPI communicator with all participating MPI ranks.
     */
    Group(const std::string& circuit,
          const std::string& population,
          const std::string& morphos,
          const SpineLengths& lengths,
          const std::shared_ptr<const AppositionSpace>& space,
          MPI_Comm comm)
        : group_metadata_(circuit, population, "", "", morphos)
        , morphologies_(morphio::Collection(morphos))
        , spine_lengths_(indexSpineLengths(circuit, population, lengths))
        , apposition_space_(space)
        , comm_(comm)
        , comm_size_(mpi::size(comm))
        , comm_rank_(mpi::rank(comm)){};

    /**
     * \brief Instantiate a new cell group
     *
     * \param circuit The full path to the file containing the cell metadata
     * \param population The population name used when loading the metadata
     * \param nodesets_path The path to the nodeset definition
     * \param nodeset_name The name for the nodeset to use
     * \param morphos The path to the morphology storage
     * \param lengths A mapping of morphology types to spine lengths
     * \param space The Appositionspace to use
     * \param comm The MPI communicator with all participating MPI ranks.
     */
    Group(const std::string& circuit,
          const std::string& population,
          const std::string& nodesets_path,
          const std::string& nodeset_name,
          const std::string& morphos,
          const SpineLengths& lengths,
          const std::shared_ptr<const AppositionSpace>& space,
          MPI_Comm comm)
        : group_metadata_(circuit, population, nodesets_path, nodeset_name, morphos)
        , morphologies_(morphio::Collection(morphos))
        , spine_lengths_(indexSpineLengths(circuit, population, lengths))
        , apposition_space_(space)
        , comm_(comm)
        , comm_size_(mpi::size(comm))
        , comm_rank_(mpi::rank(comm)){};

    /**
     * \brief Instantiate a new cell group
     *
     * \param circuit The full path to the file containing the cell metadata
     * \param population The population name used when loading the metadata
     * \param morphos The path to the morphology storage
     * \param lengths A mapping of morphology types to spine lengths
     * \param comm The MPI communicator with all participating MPI ranks.
     */
    Group(const std::string& circuit,
          const std::string& population,
          const std::string& morphos,
          const SpineLengths& lengths,
          MPI_Comm comm)
        : Group(circuit, population, morphos, lengths, createAppositionSpace(), comm) {}


    /**
     * \brief Copy the essential data of Group and load some metadata
     *
     * Morphologies are not loaded automatically; This is the
     * responsibility of the user.
     *
     * \param o The group from which to copy essential data
     * \param ids The cell indices to load metadata for
     */
    Group(const Group& o, const std::vector<int>& ids);

    /**
     * \brief Copy the essential data of Group and load metadata
     *
     * Morphologies are not loaded automatically; This is the
     * responsibility of the user.
     *
     * \param o The group from which to copy essential data
     * \param meta The cell metadata to use
     */
    Group(const Group& o, const std::vector<Metadata>& meta);

    Group() = delete;
    Group(const Group&) = delete;
    Group(Group&&) = default;
    ~Group() = default;
    Group& operator=(const Group&) = delete;
    Group& operator=(Group&&) = delete;
    const_iterator begin() const {
        return cells_.cbegin();
    };
    const_iterator end() const {
        return cells_.cend();
    };

    size_t size() const {
        return metadata_.size();
    };

    /**
     * \brief Distribute the selection and load the metadata.
     *
     * This overload can be used when the node population and global node
     * selection have been loaded in an earlier step.
     *
     * \param pop the population defined by `circuit` and `population` strings
     *            passed to the constructor.
     *
     * \param global_selection the global selection of the cells in the
     *            population. These will be distributed.
     *
     * \param check If true the routine checks that there are at least as many
     *              cells in the selection as MPI ranks in the communicator.
     *
     * \return The flat list of global cell IDs.
     */
    std::vector<int> loadMetadata(const bool check = false);

    //! reads all morphologies data from the given folder path
    void loadMorphologies(std::function<bool(const SectionType&)> accept = [](const SectionType&) {
        return true;
    });

    /** Create a AppositionSpace according to the strings passed
     * \param space Can be either
     * \param presynaptic The section type of the Appositioning side, can be
     *                    filled with one or more of the following values:
     *                    * all
     *                    * `*`
     *                    * soma
     *                    * axon
     *                    * dendrite
     *                    * basal
     *                    * apical
     * \param postsynaptic The section type of the Appositioned side, see
     *                     `presynaptic` for possible values
     */
    static auto createAppositionSpace(const std::string& space,
                                 const std::vector<std::string>& presynaptic,
                                 const std::vector<std::string>& postsynaptic)
        -> std::shared_ptr<const AppositionSpace>;

    /**
     * \brief Creates a default Apposition space: nonautaptic, axodendritic
     */
    static auto createAppositionSpace() -> std::shared_ptr<const AppositionSpace>;

    //! Outputs information relative to the segments distribution
    void showDistributionInfo();

    //! returns the apposition space
    const AppositionSpace* appositionSpace() const {
        return apposition_space_.get();
    }

    const Cell& operator[](size_t id) const {
        auto idx = index_.find(id);
        if (idx == index_.end()) {
            throw std::out_of_range(fmt::format("Cell with id {} not in memory", id));
        }
        return cells_[idx->second];
    }

    //! returns the meta data in index N from the cell file
    const std::vector<Metadata>& metadata() const {
        return metadata_;
    }

    const std::string& morphology_path() const {
        return group_metadata_.morphology_path;
    }

    const std::string& population_name() const {
        return group_metadata_.population_name;
    }

  protected:
    //! Obtain the current population ID
    static size_t indexPopulation(const std::string& filename, const std::string& population);

    static std::vector<std::pair<std::string, std::string>> known_populations_;

  private:
    static constexpr size_t page_buffer_size = 512 * (1ul << 20);


    //! Load all metadata information from the given population and selection
    void loadMetadata(const std::shared_ptr<bbp::sonata::NodePopulation>,
                      const bbp::sonata::Selection&);

    auto computeCellIndexRanges(size_t count) const -> std::vector<std::pair<size_t, size_t>>;

    //! Translate from a mtype string-based lookup to the circuit based mtype index
    static auto indexSpineLengths(const std::string& filename,
                                  const std::string& population,
                                  const SpineLengths& lengths) -> std::vector<float>;

    GlobalMetadata group_metadata_;
    morphio::Collection morphologies_;

    //! Represents the spine length added to each segments radius, given its morphology
    const std::vector<float> spine_lengths_;
    const std::shared_ptr<const AppositionSpace> apposition_space_;

    //! Array holding all cell meta data read from file (de-serialized version)
    std::vector<Metadata> metadata_;
    std::vector<Cell> cells_;

    // provides conversion between cellId and position in the array of cells
    std::unordered_map<unsigned int, size_t> index_;

    MPI_Comm comm_;
    int comm_size_;
    int comm_rank_;
};

}  // namespace cell

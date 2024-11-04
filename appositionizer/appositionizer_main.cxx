#include <algorithm>
#include <cstdlib>  // getenv
#include <filesystem>
#include <fstream>
#include <regex>

#include <fmt/ranges.h>

#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/info.h>

#include "CLI/CLI.hpp"
#include "CLI/Validators.hpp"

#include "better_mpi.h"

#include "instrumentor.h"
#include "memory_tools.h"
#include "cell_slicer.h"
#include "data/apposition.h"
#include "defines.h"
#include "io/filewriter.h"
#include "process/disentangled.h"

namespace fs = std::filesystem;

constexpr size_t MIN_TBB_THREADS = 4;

/**
 * \brief Set variables to make sure that we can run and make good use of node-level parallelism.
 *
 * Needs to be called before MPI is initialized.  Not setting these environment variables
 * may lead to degraded performance, as Intel TBB may only see a single thread.
 */
void set_tbb_parallelism_configuration() {
    if (const char* dplace = getenv("__DPLACE_")) {
        std::clog << "Using dplace behavior " << dplace << "\n";
        std::clog << "LD_PRELOAD=" << getenv("LD_PRELOAD") << "\n";
    }
    if (!getenv("OMP_NUM_THREADS")) {
        if (const char* cpus_per_task = getenv("SLURM_CPUS_PER_TASK")) {
            setenv("OMP_NUM_THREADS", cpus_per_task, 0);
        } else {
            // MPI_OPENMP_INTEROP requires OMP_NUM_THREADS to be set, we can't set it
            // ourselves, so return.
            return;
        }
    }
    setenv("MPI_OPENMP_INTEROP", "1", 1);
}

/**
 * \brief Check that we have more than two threads per rank.
 */
void assert_tbb_parallelism_configuration() {
    auto have = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    if (have < MIN_TBB_THREADS) {
        mpi::abort(
            "TBB detected {} threads. Please use a minimum of {} threads with Appositionizer!",
            have,
            MIN_TBB_THREADS);
    }
}

int main(int argc, char* argv[]) {
    using io::FileWriter;

    set_tbb_parallelism_configuration();

    mpi::MPI mpi(&argc, &argv);
    Instrumentor::init_profile();

    MPI_Barrier(MPI_COMM_WORLD);
    Instrumentor::phase_begin("main");

    // Starts the Command Line parser and quits (and dumps help) if no arguments
    std::string output_dir = "./output";
    std::string tmp_dir;
    std::string configfile;
    std::string recipefile;
    std::string appositionspace = "nonautaptic";
    std::vector<std::string> presynaptic = {"all"};
    std::vector<std::string> postsynaptic = {"all"};
    size_t jobsize = 4;

    std::string from;
    std::string from_nodeset;
    std::string to;
    std::string to_nodeset;
    std::string output_population;

    CLI::App app{"Detect appositions between cells of a brain circuit"};
    app.set_version_flag("--version", APPOSITIONIZER_VERSION_STRING);
    app.add_option("--circuit-config", configfile, "Circuit configuration file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("--recipe", recipefile, "Circuit recipe file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("--tmp-dir",
                   tmp_dir,
                   "Path for temporary output. Can point to node local storage.");
    app.add_option("--from", from, "Source population name");
    app.add_option("--from-nodeset", from_nodeset, "Source population nodeset name");
    app.add_option("--to", to, "Target population name");
    app.add_option("--to-nodeset", to_nodeset, "Target population nodeset name");
    app.add_option("--output",
                   output_dir,
                   fmt::format("Output folder for edges (default: {})",
                               output_dir));
    app.add_option("--output-population",
                   output_population,
                   fmt::format("Sonata output population name (default: '{}')", output_population));
    app.add_option("--appositionspace",
                   appositionspace,
                   fmt::format("Override the Apposition space (default: {})", appositionspace))
        ->check(CLI::IsMember({"nonautaptic", "autaptic", "axodendritic", "dendrosomatic"}));
    app.add_option("--presynaptic",
                   presynaptic,
                   "Define which section types are possible presynaptic, comma separated")
        ->delimiter(',')
        ->check(CLI::IsMember({"*", "all", "soma", "axon", "dendrite", "basal", "apical"}));
    app.add_option("--postsynaptic",
                   postsynaptic,
                   "Define which section types are possible postsynaptic, comma separated")
        ->delimiter(',')
        ->check(CLI::IsMember({"*", "all", "soma", "axon", "dendrite", "basal", "apical"}));
    app.add_option("--jobsize",
                   jobsize,
                   fmt::format("Number of cells to process in a batch (default: {})", jobsize));

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        if (mpi::first()) {
            app.exit(e);
        }
        exit(1);
    }

    assert_tbb_parallelism_configuration();

    if (appositionspace == "axodendritic") {
        appositionspace = "nonautaptic";
        presynaptic = {"axon"};
        postsynaptic = {"soma", "dendrite"};
    } else if (appositionspace == "dendrosomatic") {
        appositionspace = "nonautaptic";
        presynaptic = {"dendrite"};
        postsynaptic = {"soma", "dendrite"};
    }

    auto threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

    mpi::rank0::info(
        "Apposition Detector {} started with MPI size {}, {} thread(s) (coordinates data type: {})",
        APPOSITIONIZER_VERSION_STRING,
        mpi::size(),
        threads,
        sizeof(t_coordinate) == sizeof(float) ? "float" : "double");

    MemoryTools::AvailableMemory();

    mpi::rank0::show(
        "- Data Sizes: Point on Section:{}, Segment:{}, Section:{}, MetaData:{}, Apposition:{}, "
        "Region:{}",
        sizeof(cell::SectionPoint),
        sizeof(cell::Segment),
        sizeof(cell::Section),
        sizeof(cell::Metadata),
        sizeof(apposition::Apposition),
        sizeof(apposition::Region));

    mpi::rank0::show("- Appositionspace: {}, presynaptic: {} - postsynaptic: {}",
                     appositionspace,
                     fmt::join(presynaptic, ","),
                     fmt::join(postsynaptic, ","));

    // creates the output directory

    if (tmp_dir == "") {
        tmp_dir = output_dir;
    }

    if (mpi::first()) {
        fs::path p(output_dir);
        if (fs::exists(p)) {
            if (not fs::is_directory(p)) {
                mpi::abort("Can't create output directory {}", p.string());
            } else {
                std::regex edgefile("^(.*/)?edge(Data)?.[0-9]+$");
                for (auto& f: fs::directory_iterator(output_dir)) {
                    if (std::regex_match(f.path().string(), edgefile)) {
                        mpi::abort("Found output file from previous run: {}", f.path().string());
                    }
                }
            }
        } else if (not fs::create_directories(p)) {
            mpi::abort("Can't create output directory {}", p.string());
        }
        mpi::show("- Output folder: {}", output_dir);
        mpi::show("- Temporary output folder: {}", tmp_dir);
    }

    const auto post_metadata =
        cell::Group::GlobalMetadata::fromConfig(configfile, to, to_nodeset);

    int max_mpi_ranks = post_metadata.cell_count();

    auto needs_custom_comm = max_mpi_ranks < mpi::size();

    MPI_Comm comm_world;
    if (needs_custom_comm) {
        auto world_comm_color = int(mpi::rank() < max_mpi_ranks);
        MPI_Comm_split(MPI_COMM_WORLD, world_comm_color, mpi::rank(), &comm_world);

        if (world_comm_color != 1) {
            MPI_Comm_free(&comm_world);
            Instrumentor::phase_end("main");
            Instrumentor::finalize_profile();
            return 0;
        }
    } else {
        comm_world = MPI_COMM_WORLD;
    }

    if (needs_custom_comm) {
        mpi::rank0::info("More MPI ranks than cells detected.");
        mpi::rank0::info("Apposition Detector {} runs with effective MPI size {} and {} threads.",
                         APPOSITIONIZER_VERSION_STRING,
                         mpi::size(comm_world),
                         threads);

        mpi::rank0::info("The effectively available memory is:");
        MemoryTools::AvailableMemory(comm_world);
    }

    const auto recipe = read_recipe(recipefile);
    const auto space = cell::Group::createAppositionSpace(appositionspace, presynaptic, postsynaptic);

    SlicedCells sliced_cells(comm_world);
    {
        // SlicedCells takes ownership (i.e. copies via MPI) of
        // re-distributed segments loaded here.
        cell::Group post_cells(post_metadata, recipe.spines, space, comm_world);
        post_cells.loadMetadata(true);

        mpi::rank0::info("Reading target morphologies from {}...", post_metadata.morphology_path);
        post_cells.loadMorphologies([&](const cell::SectionType& type) {
            return post_cells.appositionSpace()->isPostSynaptic(type);
        });

        mpi::rank0::info("Starting Sort Tile Recursive algorithm for slicing...");
        CellSlicer cellSlicer(comm_world);
        sliced_cells = cellSlicer.redistribute(post_cells);
        cellSlicer.writeSlicesToTextFile(output_dir, sliced_cells);

        MemoryTools::AvailableMemory(comm_world);
    }

    const auto pre_metadata =
        cell::Group::GlobalMetadata::fromConfig(configfile, from, from_nodeset);
    cell::Group pre_cells(pre_metadata, recipe.spines, space, comm_world);
    pre_cells.loadMetadata(false);

    // We're loading a fair-share of the morphologies to later compute
    // (collectively) which MPI-ranks even need to look at which GIDs.
    // TODO move this into `process::disentangled` to get this line closer
    // to its sole reason to exist.
    mpi::rank0::info("Reading source morphologies from {}...", pre_cells.morphology_path());
    pre_cells.loadMorphologies([&](const cell::SectionType& type) {
        return pre_cells.appositionSpace()->isPreSynaptic(type);
    });

    const process::Configuration cfg{tmp_dir,
                                     output_dir,
                                     output_population,
                                     pre_metadata.population_name,
                                     post_metadata.population_name};

    process::disentangled(pre_cells, sliced_cells, cfg, jobsize, recipe.interval, comm_world);

    if (needs_custom_comm) {
        MPI_Comm_free(&comm_world);
    }

    MPI_Barrier(comm_world);
    Instrumentor::phase_end("main");
    Instrumentor::finalize_profile();
    return 0;
}

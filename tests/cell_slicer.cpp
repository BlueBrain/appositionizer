#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "../appositionizer/data/cell_group.h"
#include "../appositionizer/cell_slicer.h"
#include "../appositionizer/process/apposition_filter.h"

namespace fs = std::filesystem;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using nlohmann::json;

using cell::SectionType;
using cell::Group;

static auto datadir = fs::path(__FILE__).parent_path() / "data/";
static auto recipe_path = (datadir / "recipe.json").string();
static auto circuit_path = (datadir / "nodes.h5").string();
static auto morpho_path = (datadir / "regions" / "h5").string();
static auto nodesets_path = (datadir / "nodesets.json").string();

TEST_CASE("check that redistribute works") {
    auto recipe = read_recipe(recipe_path);
    Group g(circuit_path,
            "All",
            nodesets_path,
            "load_test",
            morpho_path,
            recipe.spines,
            Group::createAppositionSpace(),
            MPI_COMM_WORLD);
    g.loadMetadata();
    g.loadMorphologies(
        [&](const SectionType& type) { return g.appositionSpace()->isPostSynaptic(type); });
    CellSlicer slicer(MPI_COMM_WORLD);
    slicer.redistribute(g);
}

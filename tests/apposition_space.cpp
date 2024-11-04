#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "../appositionizer/data/cell_group.h"

namespace fs = std::filesystem;

using cell::SectionType;
using cell::Group;

static auto datadir = fs::path(__FILE__).parent_path() / "data/";
static auto recipe_path = (datadir / "recipe.json").string();
static auto circuit_path1 = (datadir / "nodes_ab.h5").string();
static auto circuit_path2 = (datadir / "nodes_b.h5").string();
static auto morpho_path = (datadir / "regions" / "h5").string();

class TestGroup: public Group {
  public:
    using Group::Group;

    static void reset_population_index() {
        known_populations_.clear();
    }
};

auto load_group(const fs::path& path, const std::string& pop) {
    static const auto recipe = read_recipe(recipe_path);
    static const auto space = TestGroup::createAppositionSpace("nonautaptic", {}, {});
    TestGroup::GlobalMetadata m(path, pop, "", "", morpho_path);
    TestGroup g(m, recipe.spines, space, MPI_COMM_WORLD);
    g.loadMetadata(true);
    g.loadMorphologies([&](const SectionType&) { return true; });
    return g;
}

auto soma(const TestGroup& g, size_t id) {
    return g[id][0][0].id();
}

auto section(const TestGroup& g, size_t id = 0, size_t section = 1, size_t segment = 0) {
    return g[id][section][segment].id();
}

TEST_CASE("check that populations are accounted for when determining overlaps") {
    TestGroup::reset_population_index();

    const auto other_space = TestGroup::createAppositionSpace("nonautaptic", {"all"}, {"all"});
    const auto self_space = TestGroup::createAppositionSpace("autaptic", {"all"}, {"all"});

    auto test_cases = std::vector<std::pair<fs::path, std::string>>{{circuit_path1, "A"},
                                                                    {circuit_path2, "B"}};

    for (const auto& [path, pop]: test_cases) {
        const auto g1 = load_group(circuit_path1, "B");
        const auto g2 = load_group(path, pop);

        // Two soma can never be compared
        CHECK_FALSE(other_space->canSegmentsBeCompared(soma(g1, 0), soma(g2, 0)));
        CHECK_FALSE(other_space->canSegmentsBeCompared(soma(g1, 0), soma(g2, 1)));
        CHECK_FALSE(self_space->canSegmentsBeCompared(soma(g1, 0), soma(g2, 0)));
        CHECK_FALSE(self_space->canSegmentsBeCompared(soma(g1, 0), soma(g2, 1)));

        // Should not create appositions on same cell
        CHECK_FALSE(other_space->canSegmentsBeCompared(soma(g1, 0), section(g1, 0)));
        CHECK_FALSE(other_space->canSegmentsBeCompared(section(g1, 0), soma(g1, 0)));
        CHECK_FALSE(other_space->canSegmentsBeCompared(section(g1, 1, 1), section(g1, 1, 5)));

        // Same population, different cell
        CHECK(other_space->canSegmentsBeCompared(soma(g1, 0), section(g1, 1)));
        CHECK(other_space->canSegmentsBeCompared(section(g1, 0), section(g1, 1)));

        // Different populations
        CHECK(other_space->canSegmentsBeCompared(soma(g1, 0), section(g2, 0)));
        CHECK(other_space->canSegmentsBeCompared(section(g1, 0), section(g2, 0)));
        CHECK(other_space->canSegmentsBeCompared(section(g1, 0), section(g2, 1)));
        CHECK(other_space->canSegmentsBeCompared(section(g1, 0), soma(g2, 0)));

        // Autaptic appositions
        CHECK(self_space->canSegmentsBeCompared(soma(g1, 0), section(g1, 0, 2)));
        CHECK(self_space->canSegmentsBeCompared(section(g1, 0, 3), soma(g1, 0)));
        CHECK(self_space->canSegmentsBeCompared(section(g1, 1, 1), section(g1, 1, 5)));

        // Special case: autaptic appositions need to be separated more
        CHECK_FALSE(self_space->canSegmentsBeCompared(section(g1, 0, 1), section(g1, 0, 2)));

        // Different cells in the same population
        CHECK_FALSE(self_space->canSegmentsBeCompared(section(g1, 1, 0), section(g1, 0, 2)));

        // Different populations
        CHECK_FALSE(self_space->canSegmentsBeCompared(soma(g1, 0), section(g2, 0)));
        CHECK_FALSE(self_space->canSegmentsBeCompared(section(g1, 0), section(g2, 0)));
        CHECK_FALSE(self_space->canSegmentsBeCompared(section(g1, 0), section(g2, 0, 3)));
        CHECK_FALSE(self_space->canSegmentsBeCompared(section(g1, 0), soma(g2, 0)));

        TestGroup::reset_population_index();
    }
}

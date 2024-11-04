#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <filesystem>

#include "../appositionizer/recipe.h"

namespace fs = std::filesystem;

using Catch::Matchers::WithinULP;

TEST_CASE("Parsing the recipe") {
    auto filename = GENERATE("recipe.json", "recipe.yaml");
    auto pth = fs::path(__FILE__).parent_path() / "data" / filename;

    auto res = read_recipe(pth.string());

    SECTION("Inter-bouton Intervals") {
        REQUIRE_THAT(res.interval.minDistance, WithinULP(5., 1));
        REQUIRE_THAT(res.interval.maxDistance, WithinULP(7., 1));
        REQUIRE_THAT(res.interval.regionGap, WithinULP(5., 1));
    }

    SECTION("Spine lengths") {
        const auto& lengths = res.spines;
        REQUIRE(lengths.size() == 9);

        auto val = lengths.find("L23_MC");

        REQUIRE(val != lengths.end());

        REQUIRE_THAT(val->second, WithinULP(2.5, 1));
    }
}

TEST_CASE("Parsing an incomplete recipe") {
    auto filename = fs::path(__FILE__).parent_path() / "data" / "recipe_incomplete.json";

    SECTION("Read attempt") {
        REQUIRE_THROWS(read_recipe(filename.string()));
    }
}

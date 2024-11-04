#include <fstream>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "../appositionizer/data/cell.h"
#include "../appositionizer/data/apposition.h"
#include "../appositionizer/process/apposition_filter.h"

namespace fs = std::filesystem;

using nlohmann::json;

using cell::Metadata;
using apposition::Region;
using apposition::Apposition;

TEST_CASE("Round-robin serialization") {
    auto datadir = fs::path(__FILE__).parent_path() / "data" / "regions" / "case_01";

    SECTION("Cell metadata") {
        auto file = datadir / "metadata.json";
        REQUIRE(fs::exists(file));
        REQUIRE(fs::is_regular_file(file));

        std::ifstream i(file.string());
        json jin;
        i >> jin;

        Metadata m = jin[0];
        json jout = m;

        REQUIRE(jin[0].dump(2) == jout.dump(2));
    }

    SECTION("Apposition") {
        auto file = datadir / "appositions.json";
        REQUIRE(fs::exists(file));
        REQUIRE(fs::is_regular_file(file));

        std::ifstream i(file.string());
        json jin;
        i >> jin;

        Apposition t = jin[0];
        json jout = t;

        REQUIRE(jin[0].dump(2) == jout.dump(2));
    }

    SECTION("Region") {
        auto file = datadir / "regions.json";
        REQUIRE(fs::exists(file));
        REQUIRE(fs::is_regular_file(file));

        std::ifstream i(file.string());
        json jin;
        i >> jin;

        Region r = jin[0];
        json jout = r;

        REQUIRE(jin[0].dump(2) == jout.dump(2));
    }
}

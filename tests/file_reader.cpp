#include <catch2/catch_test_macros.hpp>

namespace fs = std::filesystem;

#include "../appositionizer/io/filereader.h"

static auto datadir = fs::path(__FILE__).parent_path() / "data";
static auto edgedir = datadir / "empty_edges";

TEST_CASE("io::FileReader; empty") {
    auto file_reader = io::FileReader(edgedir, ".0");

    CHECK(file_reader.count(12) == 0ul);
    CHECK(file_reader.offset(12) == 0ul);
}

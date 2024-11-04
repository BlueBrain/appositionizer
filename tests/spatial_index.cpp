#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../appositionizer/data/segment.h"
#include "../appositionizer/spatial_grid_hash.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;

using cell::Segment;
using cell::SegmentId;

TEST_CASE("Query the spatial index") {
    static const size_t range = 20;
    std::vector<Segment> data;
    for (size_t i = 0; i < range; ++i) {
        const auto pos = -0.5f * range + i;
        data.emplace_back(SegmentId(0, 0, 2 * i, 1), Point(pos, 0., 0.), Point(-pos, 0., 0.), 4, 6);
        data.emplace_back(
            SegmentId(0, 0, 2 * i + 1, 1), Point(0., pos, 0.), Point(0., -pos, 0.), 4, 6);
    }

    SpatialGridHash<Segment, t_coordinate> index(5.);
    auto dataCopy = data;
    index.createIndex(std::move(dataCopy));

    SECTION("Query the index for all elements") {
        std::array<t_coordinate, 3> p1{-6., -6., -6.};
        std::array<t_coordinate, 3> p2{6., 6., 6.};
        Segment::BoundingBox box(p1, p2);
        auto result = index.query(box);
        REQUIRE(result.size() == data.size());
    }
}

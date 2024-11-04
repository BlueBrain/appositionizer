#include <catch2/catch_test_macros.hpp>

#include "../appositionizer/data/section.h"
#include "../appositionizer/data/segment.h"

using cell::Section;
using cell::SectionPoint;
using cell::SectionPointStorage;
using cell::SegmentId;

TEST_CASE("Iterating over a Section") {
    std::vector<SectionPoint> bps{
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 1, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 2, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 3, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 4, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 5, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 6, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 7, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 8, 1, 0)},
    };
    Section store{0, 0.0, std::move(bps)};

    SECTION("Consecutive iteration") {
        size_t i = 0;
        for (const auto s: store) {
            REQUIRE(s.segment() == i++);
        }
        REQUIRE(i == 8);
    }
}

TEST_CASE("Iterating over a SectionPointStorage") {
    std::vector<SectionPoint> bps{
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 1, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 2, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(1, 2, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(2, 2, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(2, 2, 0, 1, 0)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(2, 2, 0, 1, 1)},
        {Point{0., 0., 0.}, 0., 0., 0., 0., SegmentId(2, 2, 0, 1, 1)},
    };
    SectionPointStorage store{std::move(bps)};

    auto i = store.begin();

    SECTION("First element") {
        auto s = *i;
        REQUIRE(s.cell() == 1);
        REQUIRE(s.section() == 1);
        REQUIRE(s.population() == 0);
    }

    ++i;

    SECTION("Second element") {
        auto s = *i;
        REQUIRE(s.cell() == 1);
        REQUIRE(s.section() == 1);
        REQUIRE(s.population() == 0);
    }

    ++i;

    SECTION("New section") {
        auto s = *i;
        REQUIRE(s.cell() == 1);
        REQUIRE(s.section() == 2);
        REQUIRE(s.population() == 0);
    }

    ++i;

    SECTION("New cell") {
        auto s = *i;
        REQUIRE(s.cell() == 2);
        REQUIRE(s.section() == 2);
        REQUIRE(s.population() == 0);
    }

    ++i;

    SECTION("New population") {
        auto s = *i;
        REQUIRE(s.cell() == 2);
        REQUIRE(s.section() == 2);
        REQUIRE(s.population() == 1);
    }

    ++i;

    SECTION("Out of range") {
        REQUIRE(i == store.end());
        REQUIRE_THROWS(++i);
    }
}

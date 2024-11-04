#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../appositionizer/geometry.h"

using Catch::Matchers::WithinULP;

TEST_CASE("Bounding box behavior") {
    Box<Point::value_type> b(Point(0, 1, 0), Point(1, 0, 1));

    SECTION("Construction") {
        REQUIRE(b.min == Point(0, 0, 0));
        REQUIRE(b.max == Point(1, 1, 1));
    }

    SECTION("Point in volume") {
        REQUIRE(b.contains(Point(.5, .5, .5)));
    }

    SECTION("Point outside volume") {
        REQUIRE(not b.contains(Point(1.5, .5, .5)));
    }

    SECTION("Overlap with box entirely within") {
        Box<Point::value_type> c(Point(.25, .25, .25), Point(.75, .75, .75));

        REQUIRE(b.overlaps(c));
    }

    SECTION("Overlap with box somewhat within") {
        Box<Point::value_type> c(Point(.25, .25, .25), Point(1.75, 1.75, 1.75));

        REQUIRE(b.overlaps(c));
    }

    SECTION("Overlap with box encompassing self") {
        Box<Point::value_type> c(Point(-.25, -.25, -.25), Point(1.75, 1.75, 1.75));

        REQUIRE(b.overlaps(c));
    }

    SECTION("No overlap") {
        Box<Point::value_type> c(Point(1.25, 1.25, 1.25), Point(1.75, 1.75, 1.75));

        REQUIRE(not b.overlaps(c));
    }

    SECTION("Merging and center") {
        Box<Point::value_type> c(Point(1.25, 1.25, -1.25), Point(1.75, 2.75, 1.75));

        c.merge(b);

        REQUIRE(c.min[0] == 0);
        REQUIRE(c.min[1] == 0);
        REQUIRE(c.min[2] == -1.25);
        REQUIRE(c.max[0] == 1.75);
        REQUIRE(c.max[1] == 2.75);
        REQUIRE(c.max[2] == 1.75);

        auto center = c.center();
        REQUIRE(center[0] == 0.875);
        REQUIRE(center[1] == 1.375);
        REQUIRE(center[2] == 0.25);
    }
}

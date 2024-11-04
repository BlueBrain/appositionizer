#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../appositionizer/data/segment.h"
#include "../appositionizer/data/apposition.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;

using cell::Segment;
using cell::SegmentId;

TEST_CASE("Spine properties") {
    Segment s1(SegmentId(0, 0, 0, 1), Point(0., 0., 0.), Point(0., 0., 1.), 0.5, 0.5);

    SECTION("parallel, no overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., 1., 0.), Point(0., 1., 1.), 0.2, 0.2);

        auto props = Segment::spineProperties(s1, 127, s2, 127);
        auto surface = std::get<1>(props);  // pre-synaptic
        auto center = std::get<2>(props);   // post-synaptic
        auto spine = std::get<4>(props);

        REQUIRE_THAT(spine, WithinAbs(.3, .01));
        REQUIRE_THAT(center.x, WithinAbs(.0, .01));
        REQUIRE_THAT(center.y, WithinAbs(1., .01));
        REQUIRE_THAT(center.z, WithinAbs(.5, .01));
        REQUIRE_THAT(surface.x, WithinAbs(.0, .01));
        REQUIRE_THAT(surface.y, WithinAbs(.5, .01));
        REQUIRE_THAT(surface.z, WithinAbs(.5, .01));
    }

    SECTION("parallel, shifted, no overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., 1., 0.), Point(0., 1., 1.), 0.2, 0.2);

        auto props = Segment::spineProperties(s1, 0, s2, 255);
        auto surface = std::get<1>(props);  // pre-synaptic
        auto center = std::get<2>(props);   // post-synaptic
        auto spine = std::get<4>(props);

        // surface points are straight from the center of segments radially
        // outwards.
        auto dist = std::sqrt(.3 * .3 + 1);
        REQUIRE_THAT(spine, WithinAbs(dist, .01));
        REQUIRE_THAT(center.x, WithinAbs(.0, .01));
        REQUIRE_THAT(center.y, WithinAbs(1., .01));
        REQUIRE_THAT(center.z, WithinAbs(1., .01));
        REQUIRE_THAT(surface.x, WithinAbs(.0, .01));
        REQUIRE_THAT(surface.y, WithinAbs(.5, .01));
        REQUIRE_THAT(surface.z, WithinAbs(.0, .01));
    }

    SECTION("orthogonal, overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(.5, .2, .5), Point(-.5, .2, .5), 0.1, 0.1);

        auto props = Segment::spineProperties(s1, 127, s2, 127);
        auto surface = std::get<1>(props);  // pre-synaptic
        auto center = std::get<2>(props);   // post-synaptic
        auto spine = std::get<4>(props);

        REQUIRE_THAT(spine, WithinAbs(0., .01));
        REQUIRE_THAT(center.x, WithinAbs(.0, .01));
        REQUIRE_THAT(center.y, WithinAbs(.2, .01));
        REQUIRE_THAT(center.z, WithinAbs(.5, .01));
        REQUIRE_THAT(surface.x, WithinAbs(.0, .01));
        REQUIRE_THAT(surface.y, WithinAbs(.5, .01));
        REQUIRE_THAT(surface.z, WithinAbs(.5, .01));
    }
}

TEST_CASE("Calculate Apposition positions") {
    Segment s1(SegmentId(0, 0, 0, 1), Point(0., 0., 0.), Point(20., 0., 0.), 4, 6);

    SECTION("Surface position on segment") {
        auto p = s1.surfacePosition(128, Point(200., 123., 0.));

        REQUIRE_THAT(p.x, WithinAbs(10., 0.1));
        REQUIRE_THAT(p.y, WithinAbs(5., 0.05));
        REQUIRE_THAT(p.z, WithinAbs(0., 0.01));
    }

    SECTION("Spine length for non overlapping segments") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., 30., 0.), Point(20., 30., 0.), 9, 3);

        auto result = Segment::spineProperties(s1, 128, s2, 128);

        REQUIRE_THAT(std::get<4>(result), WithinAbs(19, 0.1));
    }

    SECTION("Spine length for overlapping segments") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., 10., 0.), Point(20., 10., 0.), 9, 3);

        auto result = Segment::spineProperties(s1, 128, s2, 128);

        REQUIRE_THAT(std::get<4>(result), WithinULP(0.f, 1));
    }
}

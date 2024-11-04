#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../appositionizer/data/segment.h"
#include "../appositionizer/data/apposition.h"

using cell::SectionPoint;
using cell::Segment;
using cell::SegmentId;
using apposition::Apposition;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;

TEST_CASE("Coordinate calculations") {
    SECTION("norm") {
        Point p(0., 1., 0.);
        REQUIRE(p.norm() == 1.0);

        p = Point(4., 3., 0.);
        REQUIRE(p.norm() == 5.0);
    }

    SECTION("angle") {
        Point p(0., 1., 0.);
        Point q(1., 1., 0.);
        REQUIRE_THAT(p.angle(q), WithinAbs(M_PI / 4., 0.001));
    }

    SECTION("subtraction") {
        Point p(0., 1., 0.);
        Point q(1., 1., 0.);
        auto s = p - q;
        REQUIRE(s.x == -1.);
        REQUIRE(s.y == 0);
        REQUIRE(s.z == 0.);

        Point t(0., 0., 0.);
        t += s;
        REQUIRE(t.x == -1.);
        REQUIRE(t.y == 0);
        REQUIRE(t.z == 0.);
    }
}

TEST_CASE("Segment and Apposition properties") {
    SectionPoint s1(Point(0., 0., 0.), 0., 30.9, 0.4, 0.8, SegmentId(1, 2, 3, 1));
    SectionPoint s2(Point(1.1, 0., 0.), 1.1, 30.9, 0.6, 0.8, SegmentId(1, 2, 3, 1));

    SECTION("Segment properties") {
        REQUIRE_THAT(s1.radius(), WithinULP(1.2f, 1));
    }

    SECTION("Serialized segment properties") {
        Segment s(s1, s2);

        unsigned char pos = 255;

        REQUIRE_THAT(s.average_radius(), WithinULP(1.3f, 1));
        REQUIRE_THAT(s.radius(0), WithinULP(1.2f, 1));
        REQUIRE_THAT(s.radius(pos), WithinULP(1.4f, 1));

        REQUIRE_THAT(s.axis().norm(), WithinULP(1.1f, 1));
    }

    SECTION("Apposition properties set right") {
        Segment s(s1, s2);

        Apposition t(s, s, {1, 2, 3, 4});
        REQUIRE(t.preCell() == 1);
    }
}

TEST_CASE("Detect overlapping segments") {
    Segment s1(SegmentId(0, 0, 0, 1), Point(0., 0., 0.), Point(0., 0., 1.), 0.5, 0.5);

    bool overlap;
    std::array<unsigned char, 4> vectors;

    SECTION("parallel, no overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., 1., 0.), Point(0., 1., 1.), 0.2, 0.2);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(!overlap);
    }

    SECTION("parallel, half overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., .1, 0.), Point(0., .1, .5), 0.2, 0.2);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);
    }

    SECTION("orthogonal, no overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(.5, 10., .5), Point(-.5, 10., .5), 0.1, 0.1);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(!overlap);
    }

    SECTION("orthogonal, overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(.5, 0., .5), Point(-.5, 0., .5), 0.1, 0.1);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);
    }
}

TEST_CASE("Test overlap setting") {
    Segment s1(SegmentId(0, 0, 0, 1), Point(0., 0., -1.), Point(0., 0., 1.), 0.5, 0.5);

    bool overlap;
    std::array<unsigned char, 4> vectors;

    auto precision = 0.1;

    SECTION("parallel, half overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., .1, -0.5), Point(0., .1, 0.5), 0.2, 0.2);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(std::abs(vectors[0] / 255.0 - 0.25) < precision);
        REQUIRE(std::abs(vectors[1] / 255.0 - 0.75) < precision);
        REQUIRE(std::abs(vectors[2] / 255.0 - 0) < precision);
        REQUIRE(std::abs(vectors[3] / 255.0 - 1) < precision);
    }

    SECTION("parallel, half overlap with larger segment") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., .1, 0), Point(0., .1, 1.5), 1.1, 1.1);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(std::abs(vectors[0] / 255.0 - 0.5) < precision);
        REQUIRE(std::abs(vectors[1] / 255.0 - 1) < precision);
        REQUIRE(std::abs(vectors[2] / 255.0 - 0) < precision);
        REQUIRE(std::abs(vectors[3] / 255.0 - 0.666) < precision);
    }

    SECTION("angled, some overlap with larger segment") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(0., .0, 0.5), Point(4., .0, 0.5), 0.5, 0.5);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(std::abs(vectors[0] / 255.0 - 0.5) < precision);
        REQUIRE(std::abs(vectors[1] / 255.0 - 1) < precision);
        REQUIRE(std::abs(vectors[2] / 255.0) < precision);
        REQUIRE(std::abs(vectors[3] / 255.0 - 0.125) < precision);
    }

    SECTION("orthogonal, overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(1., 0., 0.), Point(-1., 0., 0.), 0.2, 0.2);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(std::abs(vectors[0] / 255.0 - 0.4) < precision);
        REQUIRE(std::abs(vectors[1] / 255.0 - 0.6) < precision);
        REQUIRE(std::abs(vectors[2] / 255.0 - 0.25) < precision);
        REQUIRE(std::abs(vectors[3] / 255.0 - 0.75) < precision);
    }

    SECTION("angled, overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(1., 0., -.5), Point(-1., 0., .5), 0.2, 0.2);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(vectors[0] / 255.0 > 0.1);
        REQUIRE(vectors[1] / 255.0 < 0.9);
        REQUIRE(vectors[2] / 255.0 > 0.2);
        REQUIRE(vectors[3] / 255.0 < 0.8);
    }

    SECTION("both angled, overlap") {
        Segment s2(SegmentId(0, 0, 0, 1), Point(-1., 0., .5), Point(1., 0., 2.5), 0.5, 0.5);
        Segment s3(SegmentId(0, 0, 0, 1), Point(3., 0., .5), Point(1., 0., 2.5), 0.5, 0.5);

        std::tie(overlap, vectors) = s2.overlaps(s3);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);

        REQUIRE(vectors[0] / 255.0 > 0.8);
        REQUIRE(vectors[1] / 255.0 - 1 < precision);
        REQUIRE(vectors[2] / 255.0 > 0.8);
        REQUIRE(vectors[3] / 255.0 - 1 < precision);
    }

    SECTION("with a soma") {
        Segment s2(SegmentId(0, 0, 0, 0), Point(.25, 0., .5), Point(.25, 0., .5), 0.5, 0.5);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);

        REQUIRE(vectors[0] / 255.0 - 0.5 < precision);
        REQUIRE(vectors[1] / 255.0 - 1 < precision);
        REQUIRE(vectors[2] == 0);
        REQUIRE(vectors[3] == 0);
    }
}

TEST_CASE("Detect overlapping large segments") {
    Segment s1(SegmentId(0, 0, 0, 1), Point(0., 0., -50.), Point(0., 0., 50.), 10., 10.);

    bool overlap;
    std::array<unsigned char, 4> vectors;

    auto imax = 100;
    for (int i = 0; i <= imax; ++i) {
        auto y = 49. * sin(.5 * M_PI * i / imax);
        auto z = 49. * cos(.5 * M_PI * i / imax);

        Segment s2(SegmentId(0, 0, 0, 1), Point(18.9, -y, -z), Point(18.5, y, z), 9., 9.);

        std::tie(overlap, vectors) = s1.overlaps(s2);
        REQUIRE(overlap);
        REQUIRE(vectors[0] < vectors[1]);
        REQUIRE(vectors[2] < vectors[3]);
    }
}

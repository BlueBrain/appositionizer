#include "data/segment.h"  // Must be on the first substantial line.

#include <cassert>
#include <iostream>
#include <limits>

#include "better_mpi.h"
#include "geometry.h"

#include "data/section.h"
#include "instrumentor.h"
#include "data/segment_id.h"

#include "partial_collision.h"

#define TINY_NUM 0.00000001

using nlohmann::json;

namespace cell {

Segment::Mode Segment::defaultMode_ = Segment::Mode::sampled;

SectionPoint::BoundingBox SectionPoint::boundingBox() const {
    // basically gets the bounding box of the full segment serialized
    return Segment(*this, *(this + 1)).boundingBox();
}

Segment::BoundingBox Segment::boundingBox() const {
    BoundingBox box(p_, p_ + axis_);
    return box.expand(average_radius());
}

Point Segment::surfacePosition(unsigned char where, const Point& vector) const {
    auto along = axis();
    auto center = p1() + along * Segment::norm * where;
    auto center_radius = bareRadius(where);

    auto vector_norm = vector.norm();

    if (vector_norm <= std::numeric_limits<float>::epsilon() * 64) {
        json j = bitKey_;
        mpi::warn("Connecting norm too small for {}: ({}, {}, {})",
                  j.dump(),
                  vector.x,
                  vector.y,
                  vector.z);
        return center;
    }

    if (sectionType() == static_cast<int>(SectionType::soma)) {
        return center + vector / vector_norm * center_radius;
    }

    auto radial = Point::cross(along, Point::cross(vector, along));
    auto radial_norm = radial.norm();

    if (radial_norm <= std::numeric_limits<float>::epsilon() * 64) {
        json j = bitKey_;
        mpi::warn("Radial norm too small for {} with axis ({}, {}, {}) and vector ({}, {}, {})",
                  j.dump(),
                  along.x,
                  along.y,
                  along.z,
                  vector.x,
                  vector.y,
                  vector.z);
        radial_norm = std::numeric_limits<float>::epsilon() * 64;
    }
    return center + radial * (center_radius / radial_norm);
}

auto Segment::spineProperties(const Segment& pre,
                              unsigned char pre_percentage,
                              const Segment& post,
                              unsigned char post_percentage)
    -> std::tuple<Point, Point, Point, Point, float> {
    auto center_pre = pre.p1() + pre.axis() * pre_percentage * SectionPoint::norm;
    auto center_post = post.p1() + post.axis() * post_percentage * SectionPoint::norm;

    auto connection = center_post - center_pre;

    auto surface_pre = pre.surfacePosition(pre_percentage, connection);
    auto surface_post = post.surfacePosition(post_percentage, connection * -1.);
    auto spine = surface_post - surface_pre;

    float spine_length = (connection * spine >= 0) ? spine.norm() : 0.0;

    return std::make_tuple(center_pre, surface_pre, center_post, surface_post, spine_length);
}

auto Segment::axialDistance(const Segment& o) const
    -> std::tuple<t_coordinate, t_coordinate, t_coordinate> {
    t_coordinate sN = -1, tN = -1;

    auto u = axis();
    auto v = o.axis();
    auto w = p1() - o.p1();

    t_coordinate a = u * u;  // always >= 0
    t_coordinate b = u * v;
    t_coordinate c = v * v;  // always >= 0
    t_coordinate d = u * w;
    t_coordinate e = v * w;
    t_coordinate D = a * c - b * b;  // always >= 0
    t_coordinate sD = D;             // sc = sN / sD, default sD = D >= 0
    t_coordinate tD = D;             // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < TINY_NUM) {  // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else {  // get the closest points on the infinite lines
        sN = (b * e - c * d);
        tN = (a * e - b * d);
        if (sN <= 0.0) {  // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        } else if (sN >= sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN <= 0.0) {  // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d <= 0.0)
            sN = 0.0;
        else if (-d >= a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    } else if (tN >= tD) {  // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) <= 0.0)
            sN = 0;
        else if ((-d + b) >= a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    t_coordinate sc = (std::abs(sN) < TINY_NUM ? 0.0 : sN / sD);
    t_coordinate tc = (std::abs(tN) < TINY_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    auto dP = w + (u * sc) - (v * tc);  // = S1(sc) - S2(tc)

    return std::make_tuple(dP * dP, sc, tc);  // square of the closest distance
}

t_coordinate Segment::distance(const Segment& o) const {
    t_coordinate distance2;
    t_coordinate u;
    t_coordinate v;

    std::tie(distance2, u, v) = axialDistance(o);

    if (0.0 < u and u < 1.0 and 0.0 < v and v < 1.0) {
        return std::max(std::sqrt(distance2) - bareRadius() - o.bareRadius(),
                        static_cast<t_coordinate>(0.));
    } else {
        return std::max(std::sqrt(distance2) - bareRadius() - o.bareRadius(),
                        static_cast<t_coordinate>(0.));
    }
}

void Segment::defaultOverlapMode(Mode m) {
    if (m == Mode::standard) {
        throw std::runtime_error("cannot set default mode to 'standard'");
    }
}

auto Segment::overlaps(const Segment& o, Mode mode) const -> std::pair<bool, Extends> {
    const auto [distance2, u, v] = axialDistance(o);
    const auto radii = average_radius() + o.average_radius();
    // checks if the bounding boxes of the volumes Apposition
    if (!boundingBox().overlaps(o.boundingBox())) {
        return {false, Extends{}};
    } else if (distance2 >= radii * radii) {
        return {false, Extends{}};
    }

    if (mode == Mode::standard) {
        mode = defaultMode_;
    }
    if (mode == Mode::sampled) {
        return pointOverlaps(o);
    } else if (mode == Mode::iterative) {
        return iterativeOverlaps(o);
    } else {
        throw std::runtime_error("cannot use default mode 'standard'");
    }
}

auto Segment::pointOverlaps(const Segment& o) const -> std::pair<bool, Extends> {
    return PartialCollision::overlap(*this, o);
}

auto Segment::iterativeOverlaps(const Segment& o) const -> std::pair<bool, Extends> {
    t_coordinate spine = spineLength() + o.spineLength();
    t_coordinate radii = bareRadius() + o.bareRadius();

    t_coordinate distance2;
    t_coordinate u;
    t_coordinate v;

    std::tie(distance2, u, v) = axialDistance(o);

    if (distance2 > (spine + radii) * (spine + radii)) {
        return std::make_pair(false, Extends());
    }

    auto distance = std::sqrt(distance2) - radii;

    Extends extends;
    std::tie(extends[0], extends[1]) = boundaries(o, u, distance, spine);
    std::tie(extends[2], extends[3]) = o.boundaries(*this, u, distance, spine);

    return std::make_pair(true, extends);
}

auto Segment::boundaries(const Segment& o,
                         t_coordinate midway,
                         t_coordinate distance,
                         t_coordinate spine) const -> std::pair<unsigned char, unsigned char> {
    if (sectionType() == static_cast<int>(SectionType::soma)) {
        return std::make_pair(0, 0);
    }

    auto fraction = static_cast<unsigned char>(
        std::min<t_coordinate>(std::max(midway, static_cast<t_coordinate>(0)) * (DISCRETE_MAX + 1),
                               DISCRETE_MAX));
    auto d_start = fractional(0).distance(o);
    auto d_end = fractional(std::numeric_limits<unsigned char>::max()).distance(o);

    return std::make_pair(boundary(o, 0, fraction, d_start, distance, spine),
                          boundary(o, DISCRETE_MAX, fraction, d_end, distance, spine));
}

auto Segment::boundary(const Segment& o,
                       unsigned char start,
                       unsigned char end,
                       t_coordinate d_start,
                       t_coordinate d_end,
                       t_coordinate spine) const -> unsigned char {
    while (d_start > spine) {
        int diff = start - end;
        if (std::abs(diff) == 1) {
            if (d_start > spine and d_end <= spine)
                return end;
            throw std::runtime_error("violation of input specification");
        } else if (diff == 0) {
            if (d_end <= spine)
                return end;
            throw std::runtime_error("violation of input specification: too close");
        }

        unsigned int pivot = (start + end) / 2;
        t_coordinate d_pivot = fractional(static_cast<unsigned char>(pivot)).distance(o);
        if (d_pivot > spine) {
            start = pivot;
            d_start = d_pivot;
        } else {
            end = pivot;
            d_end = d_pivot;
        }
    }
    return start;
}

auto Segment::fractional(unsigned char where) const -> Segment {
    constexpr t_coordinate scale = 1. / (DISCRETE_MAX + 1);
    return Segment(id(),
                   p_ + axis_ * (scale * where),
                   p_ + axis_ * (scale * where + scale),
                   radius_,
                   radius_ + radius_taper_,
                   spine_length_);
}

void to_json(json& j, const SegmentId& f) {
    j = json{{"cell", f.cell()},
             {"section", f.section()},
             {"segment", f.segment()},
             {"type", f.sectionType()}};
}

void from_json(const json& j, SegmentId& f) {
    f.cell(j["cell"]);
    f.section(j["section"]);
    f.sectionType(j["type"]);
    f.segment(j["segment"]);
}

void to_json(json& j, const Segment& s) {
    j = json{{"key", s.id()},
             {"p1", s.p1()},
             {"p2", s.p2()},
             {"r1", s.bareRadius()},
             {"r2", s.bareRadius(255)},
             {"spine", s.spineLength()}};
}

void from_json(const json& j, Segment& s) {
    s = Segment(j["key"], j["p1"].get<Point>(), j["p2"].get<Point>(), j["r1"], j["r2"], j["spine"]);
}

}  // namespace cell

void to_json(json& j, const Point& p) {
    j = json{{"x", p.x}, {"y", p.y}, {"z", p.z}};
}

void from_json(const json& j, Point& p) {
    p.x = j["x"];
    p.y = j["y"];
    p.z = j["z"];
}

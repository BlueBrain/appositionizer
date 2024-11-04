#include "partial_collision.h"  // Must be on the first substantial line.

#include <cmath>
#include <utility>

#include <Eigen/Dense>

#include "data/segment.h"
#include "assert.h"

#include "data/section.h"

using Eigen::Affine;
using Eigen::AngleAxisf;
using Eigen::Matrix3Xf;
using Eigen::Scaling;
using Eigen::Transform;
using Eigen::Translation3f;
using Eigen::Vector3f;
using Eigen::VectorXf;

using cell::Segment;

using Extends = PartialCollision::Extends;
using ContinuousExtends = std::array<float, 4>;

static const auto PI = 4 * std::atan(1.0);

Transform<float, 3, Affine> build_transform(const Segment& s, const Segment& o) {
    const auto diff = s.p1() - o.p1();

    const Vector3f s_axis{s.axis().x, s.axis().y, s.axis().z};
    const Vector3f o_axis{o.axis().x, o.axis().y, o.axis().z};

    // These lead to rotating the z-axis onto the axis of `s`.
    static const auto unit_z = Vector3f::UnitZ();
    const auto source_angle = std::acos(unit_z.dot(s_axis) / s_axis.norm());
    const auto source_axis = unit_z.cross(s_axis).normalized();

    // The transformation is applied to the sample points which are given on
    // the reference cylinder, i.e. a cylinder with radius of one and axis
    // aligned with the z-axis.
    //
    //   1. Scale the reference cylinder to have the same proportions as `s`.
    //   2. Rotate the scaled reference cylinder such that its new axis is
    //      aligned with the axis of `s`.
    //   3. Transform the points into a coordinate system centered on `o.p1()`.
    //      This means that if `o` is a soma, we can check if points intersect
    //      with `o` by checking if they lie within one radius of the sphere
    //      from the origin.
    //
    // If `o` is a cylinder we need one more step:
    //   4. Transform the point cloud into a coordinate system that
    //      additionally has its z-axis aligned with the axis of `o`.
    //      Therefore, one can check if points in the cloud lie in the cylinder
    //      `o` by checking if they are inside the z-aligned cylinder with
    //      radius `o.radius`.
    //
    Transform<float, 3, Affine> t;
    t = Translation3f(diff.x, diff.y, diff.z) * AngleAxisf(source_angle, source_axis) *
        Scaling(s.average_radius(), s.average_radius(), s.length());
    if (not o.soma()) {
        const auto target_angle = std::acos(unit_z.dot(o_axis) / o_axis.norm());
        const auto target_axis = o_axis.cross(unit_z).normalized();
        t = AngleAxisf(target_angle, target_axis) * t;
    }

    return t;
}

Matrix3Xf transform(const Matrix3Xf& points, const Segment& s, const Segment& o) {
    auto t = build_transform(s, o);
    return t * points;
}


void resize(unsigned char& begin, unsigned char& end) {
    if (begin != end)
        return;

    if (begin < 2) {
        end = 4;
        begin = 0;
    } else if (end > 253) {
        begin -= (255 - end);
        begin = 251;
        end = 255;
    } else {
        begin -= 2;
        end += 2;
    }
}


/** \brief Compute the overlap of two segments `s` and `o`.
 *
 *  The assumption is that `s` has been discretized into points `ps`. The
 *  fractional (range [0, 1]) height of point `ps(i)` measured along the axis of
 *  `s` is given by `s_height(i)`.
 *
 *  The point cloud `ps` doesn't need to cover all of `s`, results can be
 *  combined later.
 *
 *  Note: the naming of `s` and `o` is consistent with `build_transform(s, o)`.
 */
template <class F>
auto overlapSegment(const Matrix3Xf& ps, const F& s_height, const Segment& o)
    -> std::pair<bool, ContinuousExtends> {
    std::array<t_coordinate, 4> extends{1.0, 0.0, 1.0, 0.0};

    const auto xs = ps.row(0).array();
    const auto ys = ps.row(1).array();
    const auto zs = ps.row(2).array() / o.length();
    const auto rs = xs.square() + ys.square();

    const auto r2 = o.average_radius() * o.average_radius();
    bool match = false;

    for (long i = 0; i < rs.size(); ++i) {
        if (rs(i) <= r2 and zs(i) >= 0.0f and zs(i) <= 1.0f) {
            match = true;
            extends[0] = std::min(zs(i), extends[0]);
            extends[1] = std::max(zs(i), extends[1]);
            extends[2] = std::min(s_height(i), extends[2]);
            extends[3] = std::max(s_height(i), extends[3]);
        }
    }

    return std::make_pair(match, extends);
}


/** \brief Compute the overlap of a segment `s` with a soma `o`.
 *
 *  The assumption is that `s` has been discretized into points `ps`. The
 *  fractional (range [0, 1]) height of point `ps(i)` measured along the axis of
 *  `s` is given by `s_height(i)`.
 *
 *  The point cloud `ps` doesn't need to cover all of `s`, results can be
 *  combined later.
 *
 *  Note: the naming of `s` and `o` is consistent with `build_transform(s, o)`.
 *
 *  Note: only computes the overlap on the cylinder.
 */
template <class F>
auto overlapSoma(const Matrix3Xf& ps, const F& s_height, const Segment& o)
    -> std::pair<bool, ContinuousExtends> {
    std::array<t_coordinate, 4> extends{0.0, 0.0, 1.0, 0.0};

    const auto xs = ps.row(0).array();
    const auto ys = ps.row(1).array();
    const auto zs = ps.row(2).array();
    const auto rs = xs.square() + ys.square() + zs.square();

    const auto r2 = o.average_radius() * o.average_radius();
    bool match = false;

    for (long i = 0; i < rs.size(); ++i) {
        if (rs(i) <= r2) {
            match = true;
            extends[2] = std::min(s_height(i), extends[2]);
            extends[3] = std::max(s_height(i), extends[3]);
        }
    }

    return std::make_pair(match, extends);
}

namespace detail {

/** \brief Compute the overlap of two segments.
 *
 *  This CRTP base solves the general problem of computing the discrete overlap
 *  of any two segments, one of which could be a soma, by breaking the problem
 *  down into more standarized tasks.
 *
 *  This CRTP base class expects:
 *
 *    * `Derived::overlapSegment(const Segment& s, const Segment& o)` which
 *      computes the continuous overlap of the segment `s` with the segment `o`.
 *      The segment `s` must be smaller than the segment `o`.
 *
 *    * `Derived::overlapSoma(const Segment& s, const Segment& o)` which
 *      computes the continuous overlap of the segment `s` with the soma `o`.
 */
template <class Derived>
class PartialCollisionCRTP {
  public:
    /**
     * \brief Four elements with the \a start, \a end of each segment
     *        overlap area
     */
    using Extends = cell::Segment::Extends;

    /** \brief Compute the overlap of two segments `s1` and `s2`.
     *
     *  This method prepares the input and delegates to
     *  `Derived::overlap{Segment,Soma}` and then converts the continuous
     *  overlaps into discrete overlaps.
     */
    static auto overlap(const cell::Segment& s1, const cell::Segment& s2)
        -> std::pair<bool, Extends>;
};


/**
 * \brief Computes continuous overlap of segments.
 *
 * This CRTP base class requires:
 *
 *   * `Derived::referencePoints()` which returns a reference to the
 *     reference points, i.e. the sample points on the reference cylinder
 *     which has radius 1.0 and axis aligned with the z-axis.
 *
 *   * `Derived::referenceHeight()` which returns a callable object that return
 *     the z-value of the reference points.
 */

template <class Derived>
class OverlapCRTP {
  public:
    /** \brief Overlap of the smaller segment `s` with the larger segment `o`. */
    static auto overlapSegment(const cell::Segment& s, const cell::Segment& o)
        -> std::pair<bool, ContinuousExtends>;

    /** \brief Overlap of the segment `s` with the soma `o`. */
    static auto overlapSoma(const cell::Segment& s, const cell::Segment& o)
        -> std::pair<bool, ContinuousExtends>;
};

class PartialCollision: public OverlapCRTP<PartialCollision>,
                        public PartialCollisionCRTP<PartialCollision> {
  public:
    static const Eigen::Matrix3Xf& referencePoints();
    static const Eigen::VectorXf& referenceHeights();

    /** \brief Create a point cloud covering the surface of the cylinder.
     *
     *  The number of points on each cap is \a disk_points, while the number of
     *  points on the mantle of the cylinder is \a mantle_points.
     */
    static Eigen::Matrix3Xf create_point_cloud(size_t disk_points, size_t mantle_points);
};


template <class Derived>
auto PartialCollisionCRTP<Derived>::overlap(const Segment& s1, const Segment& s2)
    -> std::pair<bool, Extends> {
    auto s1volume = s1.average_radius() * s1.average_radius() * PI * s1.length();
    auto s2volume = s2.average_radius() * s2.average_radius() * PI * s2.length();

    if (s1.soma() or s2.soma()) {
        if (s1volume == 0 and s2volume == 0)
            return std::make_pair(false, Extends());
    } else if (s1.length() == 0 or s2.length() == 0) {
        return std::make_pair(false, Extends());
    }

    bool intersect;
    ContinuousExtends outcome;

    if (s2.soma()) {
        std::tie(intersect, outcome) = Derived::overlapSoma(s1, s2);
        std::swap(outcome[0], outcome[2]);
        std::swap(outcome[1], outcome[3]);
    } else if (s1.soma()) {
        std::tie(intersect, outcome) = Derived::overlapSoma(s2, s1);
    } else if (s1volume < s2volume) {
        std::tie(intersect, outcome) = Derived::overlapSegment(s1, s2);
        std::swap(outcome[0], outcome[2]);
        std::swap(outcome[1], outcome[3]);
    } else {
        std::tie(intersect, outcome) = Derived::overlapSegment(s2, s1);
    }

    Extends result;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = outcome[i] * 255.999f;
    }

    if (not s1.soma()) {
        resize(result[0], result[1]);
    }
    if (not s2.soma()) {
        resize(result[2], result[3]);
    }

    return std::make_pair(intersect, result);
}

template <class Derived>
auto OverlapCRTP<Derived>::overlapSoma(const Segment& s, const Segment& o)
    -> std::pair<bool, ContinuousExtends> {
    return ::overlapSoma(transform(Derived::referencePoints(), s, o),
                         Derived::referenceHeights(),
                         o);
}

template <class Derived>
auto OverlapCRTP<Derived>::overlapSegment(const Segment& s, const Segment& o)
    -> std::pair<bool, ContinuousExtends> {
    return ::overlapSegment(transform(Derived::referencePoints(), s, o),
                            Derived::referenceHeights(),
                            o);
}

const Matrix3Xf& PartialCollision::referencePoints() {
    static const auto MANTLE_POINTS = 250;
    static const auto DISK_POINTS = 50;
    static auto points = create_point_cloud(DISK_POINTS, MANTLE_POINTS);

    return points;
}

const VectorXf& PartialCollision::referenceHeights() {
    static VectorXf heights = referencePoints().row(2);

    return heights;
}

Matrix3Xf PartialCollision::create_point_cloud(size_t disk_points, size_t mantle_points) {
    // use the golden angle to distribute points angularly [sic].
    static const float alpha = static_cast<float>(M_PI) * (3.0f - std::sqrt(5.0f));

    Matrix3Xf points(3, 2 * disk_points + mantle_points);
    for (size_t i = 0; i < disk_points; ++i) {
        auto radius = std::sqrt(i / static_cast<float>(disk_points));
        auto angle = i * alpha;
        points(0, i) = radius * std::cos(angle);
        points(1, i) = radius * std::sin(angle);
        points(2, i) = 0.;
        points(0, i + disk_points) = -points(1, i);
        points(1, i + disk_points) = -points(0, i);
        points(2, i + disk_points) = 1.;
    }

    for (size_t i = 2 * disk_points; i < 2 * disk_points + mantle_points; ++i) {
        size_t angle = i * alpha;
        points(0, i) = std::cos(angle);
        points(1, i) = std::sin(angle);
        points(2, i) = (i - 2 * disk_points) / static_cast<float>(mantle_points);
    }

    return points;
}

}  // namespace detail

auto PartialCollision::overlap(const cell::Segment& s1, const cell::Segment& s2)
    -> std::pair<bool, Extends> {
    return detail::PartialCollision::overlap(s1, s2);
}

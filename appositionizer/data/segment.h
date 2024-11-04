#pragma once

#include <array>
#include <limits>

#include <nlohmann/json.hpp>

#include "segment_id.h"
#include "geometry.h"

namespace cell {

class Section;

/**
 * \brief Compact representation of a segment in a section.
 *
 * In a section, the end point of one segment is the starting point of the
 * next segment. Therefore, when storing an entire section of segments, one
 * can skip storing the end point of each segment.
 *
 * Recall, Appositioner uses the word *section* for an unsectioning section of a cell.
 *
 * This representation can reduce the size of a section by roughly 33 percent, or
 * less if the section is only a few segments long.
 */
class SectionPoint {
  public:
    static constexpr t_coordinate norm = 1. / std::numeric_limits<unsigned char>::max();

    using BoundingBox = Box<t_coordinate>;

    SectionPoint() = default;

    SectionPoint(const Point& p,
                t_coordinate position,
                t_coordinate length,
                t_coordinate radius,
                t_coordinate spine = 0,
                SegmentId key = SegmentId(0, 0, 0, 1))
        : point_(p)
        , position_(position)
        , sectionlength_(length)
        , radius_(radius)
        , spine_length_(spine)
        , bitKey_(key) {}

    template <typename T>
    SectionPoint(const T& p,
                t_coordinate position,
                t_coordinate length,
                t_coordinate radius,
                t_coordinate spine,
                SegmentId key)
        : point_(p[0], p[1], p[2])
        , position_(position)
        , sectionlength_(length)
        , radius_(radius)
        , spine_length_(spine)
        , bitKey_(key) {}

    SectionPoint(const SectionPoint&) = default;
    SectionPoint(SectionPoint&&) = default;
    SectionPoint& operator=(const SectionPoint&) = default;
    SectionPoint& operator=(SectionPoint&&) = default;

    unsigned int population() const {
        return bitKey_.population();
    }

    unsigned int cell() const {
        return bitKey_.cell();
    }

    unsigned int section() const {
        return bitKey_.section();
    }

    unsigned int sectionType() const {
        return bitKey_.sectionType();
    }

    t_coordinate sectionLength() const {
        return sectionlength_;
    }

    unsigned int segment() const {
        return bitKey_.segment();
    }

    /** \brief Gets the coordinated for this Segment */
    Point xyz() const {
        return point_;
    }

    /**
     * \brief Returns the position of the segment along the section (in μm)
     */
    t_coordinate position() const {
        return position_;
    }

    /** \brief Gets the radius for this segment, with spine length applied */
    t_coordinate radius() const {
        return radius_ + spine_length_;
    }

    /** \brief Gets the radius independent of spinelength for this segment */
    t_coordinate bareRadius() const {
        return radius_;
    }

    /** \brief Gets the spine length */
    t_coordinate spineLength() const {
        return spine_length_;
    }

    /** \brief Gets the bitkey (64bit field) of this segment */
    SegmentId id() const {
        return bitKey_;
    }

    template <typename T>
    bool valid(const T& data) const {
        // A section point is a valid segment only if it has
        // a follower and the follower is from the same segment
        if (this == &(*(std::end(data) - 1)))
            return false;
        if (section() != (this + 1)->section() or cell() != (this + 1)->cell())
            return false;
        return true;
    }

    /** \brief Calculate a bounding box of the segment.
     *
     * This takes the radius into account and since it's faster to compute
     * a slightly larger bounding box, this routine doesn't guarantee to compute
     * the smallest bounding box. See `Segment::boundingBox`.
     *
     * \warning This requires that the `SectionPoint`s are stored contiguously
     * and that there is another `SectionPoint` following this one, e.g., the
     * expression
     *
     *     Segment(*this, *(this + 1));
     *
     * must be valid.
     */
    BoundingBox boundingBox() const;

    Point point_;

    t_coordinate position_;
    t_coordinate sectionlength_;
    t_coordinate radius_;
    t_coordinate spine_length_;
    // Necessary as we can't retrieve cellId, sectionId and segID when we
    // spatially index
    SegmentId bitKey_;
};

static_assert(std::is_trivial<SectionPoint>::value, "SectionPoint should be trivial");
static_assert(std::is_trivially_copyable<SectionPoint>::value,
              "SectionPoint should be trivially copyable");

/** \brief Representation of a segment.
 *
 * This class provides the common functionality performed with segments,
 * such as computing overlap of two segments and distances.
 *
 * Unlike `SectionPoint`, this class stores both the starting point and axis
 * of the cylinder. Along with additional information such as the radius and
 * spine length. It therefore contains a complete description of the cylinder.
 *
 * This data structure is commonly used when querying non-trivial properties
 * of a segment that require it's full shape.
 *
 * \note This can also represent a soma.
 */
class Segment {
  public:
    enum class Mode { standard, sampled, iterative };

    using BoundingBox = typename SectionPoint::BoundingBox;

    /** \brief The overlapping region along the center line.
     *
     *  Uses an `unsigned char` to represent the fractional distance along the
     *  segment, both on the pre- and post-synaptic side. The overlap on a
     *  some is always empty.
     */
    using Extends = std::array<unsigned char, 4>;

    static constexpr t_coordinate norm = SectionPoint::norm;

    Segment() = default;
    Segment(SegmentId key)
        : bitKey_(key){};
    Segment(const SectionPoint& s1, const SectionPoint& s2)
        : bitKey_(s1.id())
        , p_(s1.xyz())
        , axis_(s2.xyz() - s1.xyz())
        , length_(s2.position() - s1.position())
        , position_(s1.position())
        , radius_(s1.bareRadius())
        , radius_taper_(s2.bareRadius() - s1.bareRadius())
        , section_length_(s1.sectionLength())
        , spine_length_(s1.spineLength()){};

    template <typename T>
    Segment(SegmentId key,
            const T& p1,
            const T& p2,
            t_coordinate r1 = 0.,
            t_coordinate r2 = 0.,
            t_coordinate spine = 0.)
        : bitKey_(key)
        , p_{p1[0], p1[1], p1[2]}
        , axis_{p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]}
        , length_(axis_.norm())
        , position_(0.0)
        , radius_(r1)
        , radius_taper_(r2 - r1)
        , spine_length_(spine) {}

    Segment(const Segment&) = default;
    Segment(Segment&&) = default;
    Segment& operator=(const Segment&) = default;
    Segment& operator=(Segment&&) = default;

    unsigned int population() const {
        return bitKey_.population();
    }

    unsigned int cell() const {
        return bitKey_.cell();
    }

    unsigned int section() const {
        return bitKey_.section();
    }

    unsigned int sectionType() const {
        return bitKey_.sectionType();
    }

    unsigned int segment() const {
        return bitKey_.segment();
    }

    /** \brief returns the XYZ of the first coordinate of this Segment */
    inline const Point& p1() const {
        return p_;
    }

    /** \brief returns the XYZ of the second coordinate of this Segment */
    inline const Point p2() const {
        return p_ + axis_;
    }

    /** \brief returns the first radius this Segment, including spine length */
    inline t_coordinate average_radius() const {
        return radius_ + spine_length_ + 0.5 * radius_taper_;
    }

    inline t_coordinate radius(unsigned char where) const {
        return radius_ + spine_length_ + Segment::norm * where * radius_taper_;
    }

    /** \brief returns the radius this Segment */
    inline t_coordinate bareRadius() const {
        return radius_;
    }

    inline t_coordinate bareRadius(unsigned char where) const {
        return radius_ + Segment::norm * where * radius_taper_;
    }

    /** \brief Length of the section that the segment is located on */
    t_coordinate sectionLength() const {
        return section_length_;
    }

    t_coordinate spineLength() const {
        return spine_length_;
    }

    /**
     * \brief Calculates the connecting spine length and positions of two segments
     *
     * Calculations are based on points along the center of the segments
     * ("c" below). In a first step, points on the surface are calculated,
     * radially outward from the center points ("s" below). The distance
     * between the two surface points is then the spine length.
     *
     *  ,-----.     ,-----.
     *  |  c  s     |     |
     *  |     |\    |     |
     *  |     | \   |     |
     *  |     |  \  |     |
     *  |     |   \ |     |
     *  |     |    \s  c  |
     *  `-----'     `-----'
     *
     * \param pre segment to connect from
     * \param pre_percentage position in form of (0, Segment::norm) where
     * the connection should start
     * \param post segment to connect to
     * \param post_percentage position in form of (0, Segment::norm) where
     * the connection should terminate
     * \returns a tuple with the point on the center, surface of both pre
     * and post, and the connection length measured between surface points
     */
    static auto spineProperties(const Segment& pre,
                                unsigned char pre_percentage,
                                const Segment& post,
                                unsigned char post_percentage)
        -> std::tuple<Point, Point, Point, Point, t_coordinate>;

    /** \brief returns the point on the surface in the polar direction vector points at. */
    Point surfacePosition(unsigned char where, const Point& vector) const;

    /** \brief returns the 64 bitkey (held as a SegmentId) of this segment */
    SegmentId id() const {
        return bitKey_;
    }

    /**
     * \brief return if we are a soma segment
     */
    bool soma() const {
        return sectionType() == static_cast<int>(SectionType::soma);
    };

    /** \brief gets Bounding Box representing this cell with its radius */
    BoundingBox boundingBox() const;

    template <typename T>
    bool valid(const T& /* data */) const {
        // A segment is always valid
        return true;
    }

    /** \brief return the cylinder axis vector */
    inline Point axis() const {
        return axis_;
    }

    /** \brief return segment length */
    inline t_coordinate length() const {
        return length_;
    }

    /** \brief return the position along the section (the segment extends from position to position +
     * length) */
    inline t_coordinate position() const {
        return position_;
    }

    auto axialDistance(const Segment&) const
        -> std::tuple<t_coordinate, t_coordinate, t_coordinate>;

    t_coordinate distance(const Segment& o) const;

    static void defaultOverlapMode(Mode m);

    auto overlaps(const Segment&, Mode = Mode::standard) const -> std::pair<bool, Extends>;
    auto rngOverlaps(const Segment&) const -> std::pair<bool, Extends>;
    auto pointOverlaps(const Segment&) const -> std::pair<bool, Extends>;
    auto iterativeOverlaps(const Segment&) const -> std::pair<bool, Extends>;

  private:
    auto fractional(unsigned char) const -> Segment;

    auto boundaries(const Segment& o,
                    t_coordinate midway,
                    t_coordinate distance,
                    t_coordinate spine) const -> std::pair<unsigned char, unsigned char>;

    auto boundary(const Segment& o,
                  unsigned char start,
                  unsigned char end,
                  t_coordinate d_start,
                  t_coordinate d_end,
                  t_coordinate spine) const -> unsigned char;

    constexpr static unsigned int DISCRETE_MAX = std::numeric_limits<unsigned char>::max();

    static Mode defaultMode_;

    SegmentId bitKey_;
    Point p_;
    Point axis_;
    t_coordinate length_ = 0.0;
    t_coordinate position_ = 0.0;
    t_coordinate radius_ = 0.0;
    t_coordinate radius_taper_ = 0.0;
    t_coordinate section_length_ = 0.0;
    t_coordinate spine_length_ = 0.0;
};

void to_json(nlohmann::json& j, const cell::SegmentId& f);
void from_json(const nlohmann::json& j, cell::SegmentId& f);

void to_json(nlohmann::json& j, const cell::Segment& s);
void from_json(const nlohmann::json& j, cell::Segment& s);

}  // namespace cell

void to_json(nlohmann::json& j, const Point& p);
void from_json(const nlohmann::json& j, Point& p);

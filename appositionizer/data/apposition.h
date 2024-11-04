#pragma once

#include <nlohmann/json.hpp>

#include "data/segment.h"
#include "defines.h"

namespace apposition {

using cell::Segment;
using cell::SegmentId;

/** \brief Represents the overlapping region on the pre- and post-synaptic side.
 */
class Apposition {
  public:
    Apposition() = default;

    Apposition(const Segment& pre, const Segment& post, const std::array<unsigned char, 4>& extends);
    Apposition(const Apposition& o) = default;

    Apposition& operator=(const Apposition& o) = default;

    SegmentId preKey() const {
        return _pre_key;
    };
    SegmentId postKey() const {
        return _post_seg.id();
    };

    unsigned int preCell() const {
        return _pre_key.cell();
    };
    unsigned int preSection() const {
        return _pre_key.section();
    };
    unsigned int preSegment() const {
        return _pre_key.segment();
    };

    unsigned int postCell() const {
        return _post_seg.cell();
    };
    unsigned int postSection() const {
        return _post_seg.section();
    };
    unsigned int postSegment() const {
        return _post_seg.segment();
    };

    /** \brief Start of the segment overlap on the efferent cell */
    unsigned char preStart() const {
        return _endApposition[0];
    };
    /** \brief End of the segment overlap on the efferent cell */
    unsigned char preEnd() const {
        return _endApposition[1];
    };
    /** \brief Start of the segment overlap on the afferent cell */
    unsigned char postStart() const {
        return _endApposition[2];
    };
    /** \brief End of the segment overlap on the afferent cell */
    unsigned char postEnd() const {
        return _endApposition[3];
    };

    /** \brief The serialized segment of the afferent cell */
    Segment post() const {
        return _post_seg;
    };

    int postSectionType() const {
        return _post_seg.sectionType();
    };

    bool comparableRegion(const Apposition& other) const {
        return preKey().sameSection(other.preKey()) && postCell() == other.postCell();
    };

    /** \brief the length of the second segment of this collision */
    t_coordinate getSeg2Length() const {
        return _post_seg.length();
    };

    /** \brief the sign of the dot product of the two segments */
    bool getSegAngleDotSign() const {
        return _angleDotSign;
    }

    /** \brief This order is suitable for merging apposition regions.
     *
     *  Creating a `apposition::Region` the pre-synaptic side from individual
     *  appositions can be done efficiently if all appositions in a section are
     *  sorted to enable a simple look ahead algorithm.
     *
     *  \sa `apposition::Filter`
     */
    bool operator<(const Apposition& t) const;

  private:
    /** \brief segments bit keys */
    SegmentId _pre_key;
    Segment _post_seg;

    /** \brief segment 1 and 2 pre and post portions */
    std::array<unsigned char, 4> _endApposition;

    /** \brief The sign of the dot-product between the axes of the two Appositioning segments. */
    bool _angleDotSign;
};

void to_json(nlohmann::json& j, const Apposition& t);
void from_json(const nlohmann::json& j, Apposition& t);

using Gaps = std::vector<std::pair<const Apposition*, const Apposition*>>;

/** \brief Stores a Apposition region and the connections from that region.
 *
 *  On the pre-synaptic side multiple appo on the same section are joined
 *  into a so called Apposition region, if they are less than a user defined
 *  threshold apart.
 *
 *  Within a Apposition region `apposition::Region::Connection`s between the pre- and
 *  post-synaptic side are formed.
 */
class Region {
  public:
    using AppositionIt = std::vector<Apposition>::const_iterator;

    /** \brief A connection between a pre- and post-synaptic segment.
     *
     *  The member `pre` and `post` are the boutons on the pre- and
     *  post-synaptic segment respectively. These use the same integer
     *  base fraction to measure distance within a segment as `Extent`.
     */
    struct Connection {
        struct {
            unsigned int segment;
            unsigned char percentage;
        } pre;
        struct {
            Segment segment;
            unsigned char percentage;
        } post;

        Connection()
            : pre{0, 0}
            , post{Segment(), 0} {};
        Connection(unsigned int s, unsigned char p)
            : pre{s, p}
            , post{Segment(), 0} {};

        bool assignPostFrom(const AppositionIt& begin, const AppositionIt& end);
        bool assignPostFrom(const Gaps& gaps);
    };

    Region() = default;
    Region(const Apposition& t1, const Apposition& t2, size_t idx, size_t segs)
        : _index(idx)
        , _size(segs)
        , _start(t1)
        , _end(t2){};

    void resetIndex(size_t i) {
        _index = i;
    };
    unsigned int getStartAppositionIndex() const {
        return _index;
    };

    int getSegmentsCount() const {
        return (int) _size;
    };

    int getStartPrePortion() const {
        return (int) _start.preStart();
    };
    int getEndPrePortion() const {
        return (int) _end.preEnd();
    };

    int getCell1() const {
        return _start.preCell();
    };  // pre cell
    int getCell2() const {
        return _start.postCell();
    };  // post cell
    int getSection1() const {
        return _start.preSection();
    };

    short int getSeg1Begin() const {
        return _start.preSegment();
    };
    short int getSeg1End() const {
        return _end.preSegment();
    };

    std::vector<Connection> connections;

    const Apposition& start() const {
        return _start;
    };

  private:
    size_t _index;
    size_t _size;
    Apposition _start;
    Apposition _end;
};

void to_json(nlohmann::json& j, const apposition::Region& t);
void from_json(const nlohmann::json& j, apposition::Region& t);

void to_json(nlohmann::json& j, const apposition::Region::Connection& t);
void from_json(const nlohmann::json& j, apposition::Region::Connection& t);

}  // namespace apposition

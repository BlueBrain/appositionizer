#include "data/section.h"
#include "data/segment.h"
#include "data/apposition.h"

#include "structures.h"

using cell::Section;
using cell::Segment;

namespace io {

template <typename T, typename U>
void copy(const T& from, U& to) {
    to[0] = from[0];
    to[1] = from[1];
    to[2] = from[2];
}

float fraction_on_section(const Segment& s, unsigned char fraction) {
    const auto f = (fraction * Segment::norm * s.length() + s.position()) / s.sectionLength();
    if (std::isnan(f)) {
        return 0.0f;
    }
    return f;
}

void IndexEntry::add(const apposition::Region& r) {
    synapseCount += r.connections.size();
}

bool IndexEntry::compatible(const apposition::Region& r) {
    return cellId == r.getCell1();
}

DataEntry::DataEntry(const Segment& pre,
                     const Segment& post,
                     const Section& pre_section,
                     unsigned char pre_fraction,
                     unsigned char post_fraction)
    : preCell(pre.cell())
    , preSection(pre.section())
    , preSegment(pre.segment())
    , postCell(post.cell())
    , postSection(post.section())
    , postSegment(post.segment())
    , sectionOrder(pre_section.order())
    , preSectionFraction(fraction_on_section(pre, pre_fraction))
    , postSectionFraction(fraction_on_section(post, post_fraction))
    , sectionType((post.sectionType() & SECTION_MASK) |
                 ((pre.sectionType() & SECTION_MASK) << SECTION_SHIFT)) {
    efferent_segment_offset = pre_fraction * Segment::norm * pre.length();
    afferent_segment_offset = post_fraction * Segment::norm * post.length();

    distance_to_soma = pre_section.distanceToSoma() + pre.position() + efferent_segment_offset;

    auto props = Segment::spineProperties(pre, pre_fraction, post, post_fraction);
    copy(std::get<0>(props), preCenterPosition);
    copy(std::get<1>(props), preSurfacePosition);
    copy(std::get<2>(props), postCenterPosition);
    copy(std::get<3>(props), postSurfacePosition);
    spineLength = std::get<4>(props);
}

}  // namespace io

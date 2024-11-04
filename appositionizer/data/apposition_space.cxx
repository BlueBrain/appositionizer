#include "data/apposition_space.h"  // Must be on the first substantial line.

#include "data/section.h"

#include "defines.h"

using cell::Section;
using cell::SectionType;
using cell::SegmentId;

size_t AppositionSpace::convert(const std::vector<std::string>& l) {
    size_t res = 0;
    for (const auto& e: l) {
        if (e == "all" or e == "*") {
            res |= 0b1111;
        } else if (e == "soma") {
            res |= 0b0001;
        } else if (e == "axon") {
            res |= 0b0010;
        } else if (e == "dendrite") {
            res |= 0b1100;
        } else if (e == "basal") {
            res |= 0b0100;
        } else if (e == "apical") {
            res |= 0b1000;
        } else {
            throw std::runtime_error("unknown section type " + e);
        }
    }
    return res;
}

bool AppositionSpace::isPreSynaptic(const SectionType& type) const {
    return presynaptic_[static_cast<size_t>(type)];
}

bool AppositionSpace::isPostSynaptic(const SectionType& type) const {
    return postsynaptic_[static_cast<size_t>(type)];
}

bool AutapticAppositionSpace::canSegmentsBeCompared(SegmentId s1, SegmentId s2) const {
    // confirm that the segments come from the same cell
    if (s1.cell() != s2.cell() or s1.population() != s2.population())
        return false;

    // Cannot check soma-soma
    if (s1.sectionType() == s2.sectionType() and
        s1.sectionType() == static_cast<unsigned int>(SectionType::soma))
        return false;

    // Also, do *NOT* compare the same segment, or the adjacent neighbors
    const unsigned int pre = s1.section();
    const unsigned int post = s2.section();

    // since bluebuilder will screen for parent/child relationships, should I get rid of this check?
    if (post == pre || post == pre + 1 || post == pre - 1) {
        return false;
    }

    return true;
}

bool NonautapticAppositionSpace::canSegmentsBeCompared(SegmentId s1, SegmentId s2) const {
    // Cannot check soma-soma
    if (s1.sectionType() == s2.sectionType() and
        s1.sectionType() == static_cast<unsigned int>(SectionType::soma))
        return false;

    return (s1.cell() != s2.cell() or s1.population() != s2.population());
}

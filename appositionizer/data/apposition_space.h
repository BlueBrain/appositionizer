#pragma once

#include <bitset>
#include <vector>

#include "data/segment_id.h"

namespace cell {
class Section;
}

class AppositionSpace {
  public:
    AppositionSpace(size_t presynaptic, size_t postsynaptic)
        : presynaptic_(presynaptic)
        , postsynaptic_(postsynaptic) {}

    virtual ~AppositionSpace() {}

    /**
     * Validate that the given segments are applicable for Apposition detection in this AppositionSpace.
     * (eg: if they come from same cell)
     *
     * \param segKey1 Identifier for segment from presynaptic cell
     * \param segKey2 Identifier for segment from (potential) postsynapstic cell
     * \return true if both segments are valid for further Apposition detection; otherwise, false
     */
    virtual bool canSegmentsBeCompared(cell::SegmentId segKey1,
                                       cell::SegmentId segKey2) const = 0;

    /**
     * Tells whether this sections is of the pre-synaptic type used for this AppositionSpace
     * \param section
     * \return
     */
    virtual bool isPreSynaptic(const cell::SectionType& type) const;

    /**
     * Tells whether this sections is of the post-synaptic type used for this AppositionSpace
     * \param section
     * \return
     */
    virtual bool isPostSynaptic(const cell::SectionType& type) const;

  protected:
    static size_t convert(const std::vector<std::string>& l);

    const std::bitset<4> presynaptic_;
    const std::bitset<4> postsynaptic_;
};

/**
 * An inherited version of AppositionSpace where only segments from the same cell are accepted.  This
 * is to be used for identifying cells for potential self intersections.
 */
class AutapticAppositionSpace: public AppositionSpace {
  public:
    AutapticAppositionSpace(const std::vector<std::string>& presynaptic,
                       const std::vector<std::string>& postsynaptic)
        : AppositionSpace(convert(presynaptic), convert(postsynaptic)) {}

    /**
     * Validate that the given segments are applicable for Apposition detection in this AppositionSpace.
     * (eg: if they come from same cell)
     *
     * \param segKey1 Identifier for segment from presynaptic cell
     * \param segKey2 Identifier for segment from (potential) postsynapstic cell
     * \return true if both segments are valid for further Apposition detection; otherwise, false
     */
    bool canSegmentsBeCompared(cell::SegmentId segKey1, cell::SegmentId segKey2) const override;
};


class NonautapticAppositionSpace: public AppositionSpace {
  public:
    NonautapticAppositionSpace(const std::vector<std::string>& presynaptic,
                          const std::vector<std::string>& postsynaptic)
        : AppositionSpace(convert(presynaptic), convert(postsynaptic)) {}

    /**
     * Validate that the given segments are applicable for Apposition detection in this AppositionSpace.
     * (eg: if they come from same cell)
     *
     * \param segKey1 Identifier for segment from presynaptic cell
     * \param segKey2 Identifier for segment from (potential) postsynapstic cell
     * \return true if both segments are valid for further Apposition detection; otherwise, false
     */
    bool canSegmentsBeCompared(cell::SegmentId segKey1, cell::SegmentId segKey2) const override;
};

#pragma once

#include <stdlib.h>

#include <iostream>

#include <math.h>
#include <vector>

#include <nlohmann/json.hpp>

#include "data/cell.h"

#include "data/apposition.h"
#include "defines.h"
#include "recipe.h"

namespace cell {

class Group;
}

namespace apposition {

class Region;
class Apposition;

/** \brief Computes connections between cell from overlap data.
 *
 *  Given the purely geometric information of which segments overlap, this
 *  class computes the connections between the two sections.
 *
 *  It does this in three phases. First, appositions are joined into Apposition regions.
 *  Next, the position of the boutons on the pre-synaptic section are determined
 *  randomly within a inter bouton interval. The random number generator used
 *  for placing the connections is `RNG`. Finally, the corresponding bouton on
 *  the post-synaptic side is determined.
 *
 *  Please also refer to the User Guide for a more detailed high-level
 *  description.
 */
class Filter {
  public:
    Filter(const cell::Section&, const InterBoutonInterval&);

    /// Combines all `Apposition`es of a section into Apposition regions (with connections).
    auto filter(const std::vector<Apposition>& appositions) -> std::vector<Region>;

    /// Distance between two appositions on the pre-synaptic section.
    double getDistance(int cell,
                       int section,
                       int segment_begin,
                       unsigned char fraction_begin,
                       int segment_end,
                       unsigned char fraction_end) const;

    /** \brief Distance between two appositions on the pre-synaptic section.
     *
     *  Computes the length of the section connecting two appositions. More
     *  precisely, it computes the length from the end of the first Apposition `t1`
     *  to the beginning of the second Apposition `t2`, following the section.
     *
     *  Note, the two appositions must lie on the same section, and be sorted such
     *  that
     *
     *      t1.preSegment() <= t2.preSegment()
     *
     *  for simplified traversal of the pre-synaptic section.
     */
    double getDistance(const Apposition&, const Apposition&) const;

  private:
    t_coordinate getSegmentLength(size_t, size_t, size_t) const;

    /** \brief Joins individual appositions into Apposition regions.
     *
     *  This joins appositions that are on the pre-synaptic section within
     *  `_interval.regionGap` of each other into one region. (It does not
     *  compute the connections between the two sections).
     *
     *  \sa `Region` for how Apposition region are represented.
     *  \sa `Apposition` for details on how appositions are stored.
     */
    auto getRegions(const std::vector<Apposition>&) const -> std::vector<Region>;

    /** \brief add the pre appositions for each pair region */
    void addAppositionsForPreCell(Region& r);

    /** \brief Add the post location for all pre appositions already defined for each region */
    void addAppositionsForPostCell(Region& r, const std::vector<Apposition>& cellAppositions);

    void printRegions(void);

    /** \brief Cumulative sum of the segment length inside the Apposition region.
     *
     *  A `Region` defines a Apposition region in terms of a start segment and
     *  an end segment. There might be several segments. This function
     *  compute, for each segment in the Apposition region the distance from the
     *  beginning of the region to either the end of the segment or the end of
     *  the Apposition region, whichever comes first. Distance is measured long the
     *  section.
     */
    auto getRegionBins(const Region& r) const -> std::vector<float>;
    /** \brief creates the list of the gaps exist between pre segment location for a specific u
     * region */
    auto getRegionGaps(const Region& r, const std::vector<Apposition>& cellAppositions) const -> Gaps;

    /** \brief Section to process */
    const cell::Section& _section;

    /** \brief For Inter-Bouton interval, holds min and max spacing values and the gap length that
     * determines when regions merge together */
    const InterBoutonInterval _interval;
};

}  // namespace apposition

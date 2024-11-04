#pragma once

#include "data/segment.h"


/**
 * \brief Provides methods to calculate the overlap between two segments
 *
 * The overlap between segments is given by the region along the cylinder
 * axis that a cylinder-cylinder intersection can be projected unto.
 *
 * For practical purposes (i.e. historical reasons), the cylinder axis is
 * divided into 256 segments, and the overlap region is returned as (\a start,
 * \a end) values for each of the two segments. The method uses the following
 * methodology:
 *
 * \image html overlap_sampling.png "Using a point cloud to sample cylinder intersection"
 *
 * Given two segments \a A and \a B (a), a number of points on the surface
 * of cylinder \a A is generated, preferably in the reference frame of \a
 * A (b). The point cloud is then rotated and translated into the reference
 * frame of \a B, where the overlap is determined (c). The results are then
 * returned along each cylinder axis (d).
 *
 * If one of the segments is the soma, both \a start and \a end will be set
 * to 0 (zero).
 *
 * \todo Overlap detection between two somas is currently not supported. As
 *       use-cases evolve, this should be added.
 */
class PartialCollision {
  public:
    /**
     * \brief Four elements with the \a start, \a end of each segment
     *        overlap area
     */
    using Extends = cell::Segment::Extends;

    /**
     * \brief Detect segment overlap using a fixed point cloud
     *
     * Prefered overlap detection method. Uses a fixed point cloud,
     * consisting of two spirals for the end-caps, and a helix around the
     * mantel of the segment cylinder. The point cloud is scaled to the
     * actual segment shape, and then translated for overlap detection.
     *
     * \param s1 A cell segment
     * \param s2 A cell segment
     * \return A pair indicating if the segments overlap, and an array
     *         giving the extend of the overlap
     */
    static auto overlap(const cell::Segment& s1, const cell::Segment& s2)
        -> std::pair<bool, Extends>;
};

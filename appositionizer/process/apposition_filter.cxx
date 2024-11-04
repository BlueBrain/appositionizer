#include "apposition_filter.h"  // Must be on the first substantial line.

#include <Random123/threefry.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fmt/format.h>
#include <limits>
#include <random>

#include "data/section.h"
#include "data/cell.h"
#include "data/cell_group.h"
#include "data/segment.h"
#include "data/segment_id.h"
#include "data/apposition.h"


/** \brief Convenience wrapper around Random123.
 *
 *  Random123 is a counter-based random number generator. Essentially, this
 *  means that given a number `i` it evaluates a deterministic `f(i)` such that
 *  the sequence
 *
 *      [f(i) for i in range(n)]
 *
 *  has excellent statistical properties. As a consequence Random123 is both a
 *  parallel pseudo random number generator and computes an easily reproducibly
 *  sequence of numbers. They've published a very nicely written
 *  [research paper](https://doi.org/10.1145/2063384.2063405).
 *
 *  \sa https://github.com/DEShawResearch/random123
 */
class RNG {
  public:
    using result_type = size_t;

    RNG(size_t seed1, size_t seed2) {
        r123::Threefry4x64::ukey_type uk = {seed1, seed2, seed1, seed2};
        key_ = uk;
    };

    inline size_t operator()() {
        if (elem_ == 0) {
            for (size_t i = 0; i < counter_.size(); ++i) {
                if (counter_[i] == max()) {
                    counter_[i] = 0;
                } else {
                    counter_[i] += 1;
                    break;
                }
            }
            elem_ = counter_.size();
            values_ = rng_(counter_, key_);
        }
        return values_[--elem_];
    };

    static constexpr size_t max() {
        return std::numeric_limits<size_t>::max();
    }

    static constexpr size_t min() {
        return std::numeric_limits<size_t>::min();
    }

  private:
    r123::Threefry4x64 rng_ = {};
    r123::Threefry4x64::ctr_type counter_ = {{}};
    r123::Threefry4x64::key_type key_ = {{}};
    r123::Threefry4x64::ctr_type values_ = {{}};
    size_t elem_ = 0;
};


namespace apposition {

Filter::Filter(const cell::Section& section, const InterBoutonInterval& interval)
    : _section(section)
    , _interval(interval) {}

auto Filter::filter(const std::vector<Apposition>& appositions) -> std::vector<Region> {
    auto regions = getRegions(appositions);
    for (auto& r: regions) {
        addAppositionsForPreCell(r);
        addAppositionsForPostCell(r, appositions);
    }
    return regions;
}

// Defines the pre-post regions
auto Filter::getRegions(const std::vector<Apposition>& appositions) const -> std::vector<Region> {
    std::vector<Region> regions;
    size_t i = 0;
    while (i < appositions.size()) {
        size_t segments = 1;
        // look ahead until we find a Apposition we can't add to the current
        // region
        Apposition start = appositions[i];
        Apposition end = appositions[i];
        for (; i + segments < appositions.size(); ++segments) {
            // compare last Apposition in region (region starts at i, spans
            // (segments - 1) appositions) to the next candidate
            Apposition candidate = appositions[i + segments];
            if (not(end.comparableRegion(candidate) and
                    getDistance(end, candidate) < _interval.regionGap)) {
                break;
            }
            // we only want to replace "end" if "candidate" is on a new
            // segment or stretches further than the current "end"
            if (candidate.preSegment() > end.preSegment() or candidate.preEnd() > end.preEnd()) {
                end = candidate;
            }
        }
        regions.emplace_back(start, end, i, segments);
        i += segments;
    }
    return regions;
}

double Filter::getDistance(const Apposition& t1, const Apposition& t2) const {
    // I.e. segments have to be sorted by segment Id
    assert(t1.preSegment() <= t2.preSegment());
    return getDistance(t1.preCell(),
                       t1.preSection(),
                       t1.preSegment(),
                       t1.preEnd(),
                       t2.preSegment(),
                       t2.preStart());
}

double Filter::getDistance(int cell,
                           int section,
                           int segment_begin,
                           unsigned char fraction_begin,
                           int segment_end,
                           unsigned char fraction_end) const {
    if (segment_begin == segment_end) {
        return getSegmentLength(cell, section, segment_begin) * Segment::norm *
               std::max(fraction_end - fraction_begin, 0);
    }
    double dis = (getSegmentLength(cell, section, segment_begin) * (255 - fraction_begin) +
                  getSegmentLength(cell, section, segment_end) * fraction_end) *
                 Segment::norm;
    for (auto i = segment_begin + 1; i < segment_end; ++i)
        dis += getSegmentLength(cell, section, i);
    return dis;
}

t_coordinate Filter::getSegmentLength(size_t nId, size_t brId, size_t sId) const {
    if (_section.empty() || _section[0].cell() != nId || _section[0].section() != brId) {
        throw std::runtime_error("Wrong section processed in apposition::Filter!");
    }
    return _section[sId].length();
}

// Define an array of bins for the segment lengths compose each region
// The array is stored inside the region structure and is (partial sums = addition of the previous
// length)
auto Filter::getRegionBins(const Region& r) const -> std::vector<float> {
    std::vector<float> regionBins;

    const float f = 1. / 255.;

    auto n = r.getCell1();
    auto b = r.getSection1();
    auto s1 = r.getSeg1Begin();
    auto s2 = r.getSeg1End();
    auto a1 = (float) r.getStartPrePortion();
    auto a2 = (float) r.getEndPrePortion();

    // if(s2-s1+1>200 || s2-s1+1<=0) printf("rank=%d bins size=%d\n",rank,s2-s1+1);
    assert(s1 <= s2);

    if (s1 == s2) {
        auto l = getSegmentLength(n, b, s1);
        regionBins.push_back((a2 - a1) * l * f);
    } else {
        float dis = f * getSegmentLength(n, b, s1) * (255 - a1);
        regionBins.push_back(dis);
        for (auto i = s1 + 1; i < s2; i++) {
            dis += getSegmentLength(n, b, i);
            regionBins.push_back(dis);
        }
        dis += f * getSegmentLength(n, b, s2) * a2;
        regionBins.push_back(dis);
    }
    return regionBins;
}

// add the pre appositions for each pair region
void Filter::addAppositionsForPreCell(Region& r) {
    auto preS1 = r.getSeg1Begin();
    int p1 = r.getStartPrePortion();
    int p2 = r.getEndPrePortion();

    // reseed the RNG using the first appositions' pre- and post-segment ID. Since there
    // should be only one Apposition/overlap between two cyclinders, this shoulb be a unique
    // seed.
    const auto& starting_Apposition = r.start();
    RNG rng(starting_Apposition.preKey().compact(),
            starting_Apposition.postKey().compact());
    // distance parameters are taken from the recipe
    std::uniform_real_distribution<> distance(_interval.minDistance, _interval.maxDistance);

    // regionBins includes the length of the region, from 1st segment to segment N (for position N
    // in array);
    auto regionBins = getRegionBins(r);  // sets a bin with length of the segment compose the region
                                         // (partial sums = addition of the previous length)
    auto region_length = regionBins.back();

    double progress = 0;
    std::vector<double> steps;

    while (progress <= region_length) {
        progress += distance(rng);

        if (steps.empty()) {
            // if first random position is less than total length available in region, try to place
            // it closer to the end
            if (progress < region_length) {
                progress *= 0.5;
            } else {
                // if 1st random position is greater than total length available on region, place in
                // center of region and finish
                progress = region_length * 0.5;
                steps.push_back(progress);
                break;
            }
        }

        // if we are continuing to work with this region...
        if (progress <= region_length) {
            steps.push_back(progress);
        }
    }

    if (steps.empty()) {
        throw std::runtime_error(fmt::format(
            "Did not generate any connection between {} and {} for a region with length {}",
            r.getCell1(),
            r.getCell2(),
            region_length));
    }

    // for the pre cells
    for (size_t i = 0; i < steps.size(); i++) {  // for pre
        unsigned char dd;
        auto hit = std::find_if(regionBins.begin(), regionBins.end(), [&](auto b) {
            return b >= steps[i];
        });
        if (hit == regionBins.end()) {
            throw std::runtime_error(
                fmt::format("Algorithmic failure when generating connections between {} and {}",
                            r.getCell1(),
                            r.getCell2()));
        }
        auto k = std::distance(regionBins.begin(), hit);

        if (i >= 1) {                      // not the first Apposition in the region
            if (regionBins.size() != 1) {  // first point
                if (k == 0) {
                    const auto c = r.connections[i - 1];
                    dd = (int) c.pre.percentage +
                         steps[i] * (255.0 - (int) c.pre.percentage) / regionBins[k];
                } else if (k == static_cast<int>(regionBins.size()) - 1) {
                    dd = ((steps[i] - regionBins[k - 1]) / (regionBins[k] - regionBins[k - 1])) *
                         r.getEndPrePortion();  // last point
                } else {
                    dd = (steps[i] - regionBins[k - 1]) * 255.0 /
                         (regionBins[k] - regionBins[k - 1]);
                }
            } else {
                const auto c = r.connections[i - 1];
                dd = (int) c.pre.percentage +
                     steps[i] * (r.getEndPrePortion() - (int) c.pre.percentage) / regionBins[k];
            }
        } else {  // first Apposition
            // dd = regions[u]->getStartPrePortion();

            if (regionBins[k] < 0.0001)
                dd = r.getStartPrePortion();
            else {
                if (k == 0) {
                    if (regionBins.size() == 1)
                        // calculates the exact distance to the Apposition point (by using a percentage
                        // of length) (this is because regionBins only hold the full length of
                        // previous segments while steps hold the percentage of the segment length
                        // where the Apposition is;
                        dd = p1 + steps[i] * (p2 - p1) / regionBins[k];
                    else
                        dd = p1 + steps[i] * (255 - p1) / regionBins[k];
                }
                // else dd=0;
                else if (k == static_cast<int>(regionBins.size()) - 1)
                    dd = ((steps[i] - regionBins[k - 1]) / (regionBins[k] - regionBins[k - 1])) *
                         r.getEndPrePortion();  // last point
                else
                    dd = (steps[i] - regionBins[k - 1]) * 255.0 /
                         (regionBins[k] - regionBins[k - 1]);
            }
        }

        r.connections.emplace_back(preS1 + k, dd);
    }
}

// creates the list of the gaps exist between pre segment location for a specific u region

auto Filter::getRegionGaps(const Region& r, const std::vector<Apposition>& cellAppositions) const -> Gaps {
    Gaps gs;
    for (unsigned int t = r.getStartAppositionIndex();
         t < r.getStartAppositionIndex() + r.getSegmentsCount() - 1;
         ++t) {
        auto AppositionPtr1 = &(cellAppositions.at(t));
        auto AppositionPtr2 = &(cellAppositions.at(t + 1));
        auto s1 = AppositionPtr1->preSegment();
        auto s2 = AppositionPtr2->preSegment();
        auto la1 = (int) AppositionPtr1->preStart();
        auto ha1 = (int) AppositionPtr1->preEnd();
        auto la2 = (int) AppositionPtr2->preStart();
        auto ha2 = (int) AppositionPtr2->preEnd();

        if ((s1 == s2 && la1 <= la2 && la2 > ha1 + 1) ||
            (s1 == s2 && la1 >= la2 && la1 >= ha2 + 1)  // BM: these >= and +1 dont make sense!! see
                                                        // paper!
            || (s1 != s2)) {
            gs.emplace_back(AppositionPtr1, AppositionPtr2);
        }
    }
    return gs;
}

// Add the post location for all pre appositions already defined for each region

void Filter::addAppositionsForPostCell(Region& r, const std::vector<Apposition>& appositions) {
    auto gaps = getRegionGaps(r, appositions);

    auto begin = appositions.begin() + r.getStartAppositionIndex();
    auto end = begin + r.getSegmentsCount();

    for (auto& c: r.connections) {
        // BM 22/11/2013 region gap means that there are consecutive appositions without consecutive
        // segments ie the inter-bouton interval placed a Apposition that corresponds to a segment that
        // initially was included in the region but didn't have any Apposition
        //(one explanation may be neighbor segments have appositions which are close enough to form a
        // region)
        if (not c.assignPostFrom(begin, end)) {
            if (not c.assignPostFrom(gaps)) {
                throw std::runtime_error(
                    fmt::format("Could not match a connection between {} and {}",
                                r.getCell1(),
                                r.getCell2()));
            }
        }
    }
}

}  // namespace apposition

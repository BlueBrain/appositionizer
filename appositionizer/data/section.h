#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "data/segment.h"
#include "defines.h"
#include "geometry.h"

namespace cell {

/**
 * \brief Generalized segment storage
 *
 * Convenience storage structure for several sections. Iterating over the
 * storage will yield valid Segments.
 *
 * As the internal storage object is immutable, iterators into this
 * container should remain valid for the storage lifetime.
 *
 * \warning {
 *   It is assumed that the SectionPoint used when constructing a Section in
 *   consecutive, contiguous order for each section stored.  No checks are
 *   performed to ensure this.
 * }
 */
class SectionPointStorage {
  public:
    class Iterator {
      public:
        using T = std::vector<SectionPoint>::const_iterator;

        explicit Iterator(const T& begin, const T& end)
            : it_(begin)
            , end_(end){};

        Segment operator*() const {
            return {*it_, *(it_ + 1)};
        }

        bool operator==(const Iterator& o) const {
            return it_ == o.it_;
        }

        bool operator!=(const Iterator& o) const {
            return it_ != o.it_;
        }

        Iterator& operator++() {
            if (it_ >= end_) {
                // This would be a gross oversight, anyways
                throw std::out_of_range("Exceeded section boundaries");
            }
            ++it_;
            // end_ should always point to the last SectionPoint in the
            // storage, and be valid.
            while (it_ < end_ && !it_->id().sameSection((it_ + 1)->id())) {
                ++it_;
            }
            return *this;
        }

      private:
        T it_;
        const T end_;
    };

    SectionPointStorage(std::vector<SectionPoint>&& s)
        : points_(std::move(s)){};

    /** \brief Return `true` if the section does not contain enough points to form
     * segments. */
    bool empty() const {
        return points_.size() <= 1;
    };

    Iterator begin() const {
        return Iterator(points_.begin(), points_.end() - 1);
    }

    Iterator end() const {
        return Iterator(points_.end() - 1, points_.end() - 1);
    }

    /**
     * \brief Returns a SectionPoint at the given index
     */
    const SectionPoint& raw(size_t idx) const {
        return points_[idx];
    }

    /**
     * \brief Returns all SectionPoints stored
     */
    const std::vector<SectionPoint>& raw() const {
        return points_;
    }

  protected:
    const std::vector<SectionPoint> points_;
};

/**
 * \brief Representation of a section of a cell
 *
 * Augments the underlying SectionPointStorage by section-specific properties
 * such as the "section order" and the distance of the first segment of the
 * section from the soma.
 *
 * \warning {
 *   It is assumed that the SectionPoint used when constructing a Section
 *   belong to the same section and are in consecutive, contiguous order.
 *   No checks are performed to ensure.
 * }
 */
class Section: public SectionPointStorage {
  public:
    using BoundingBox = typename SectionPoint::BoundingBox;

    Section(int order, t_coordinate distance, std::vector<SectionPoint>&& s)
        : SectionPointStorage(std::move(s))
        , order_(order)
        , distance_(distance) {}

    /** \brief The size of the Section in actual segments. */
    size_t size() const {
        return points_.size() - 1;
    };

    Segment operator[](size_t idx) const {
        return {points_[idx], points_[idx + 1]};
    }

    /** \brief return the type of this section (Axon, Dendrite, Soma, Epic) */
    SectionType type() const {
        return static_cast<SectionType>(points_.front().sectionType());
    }

    const BoundingBox extends() const {
        auto extends = points_.front().boundingBox();
        // we don't take the last one cause the boundingBox of the last one use the next one (off by
        // one)
        for (size_t i = 1; i < points_.size() - 1; ++i) {
            extends.merge(points_[i].boundingBox());
        }
        return extends;
    }

    /** \brief gets the section length */
    t_coordinate length() const {
        return points_.back().position();
    }

    t_coordinate distance(size_t segment) const {
        return points_[segment].position();
    }

    t_coordinate distanceToSoma() const {
        return distance_;
    }

    /** \brief returns the section Id */
    int id() const {
        return points_.front().section();
    }

    /** \brief returns -1 or the section Id of order section */
    int order() const {
        return order_;
    }

    bool overlaps(const BoundingBox& box) const;

  private:
    const int order_;
    const t_coordinate distance_;
};

}  // namespace cell

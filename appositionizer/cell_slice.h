#pragma once

#include "data/cell_group.h"
#include "data/segment.h"
#include "spatial_grid_hash.h"

/**
 * This class contains the output of the cell_slicer
 */
class SlicedCells {
    using SectionPoint = cell::SectionPoint;
    using Segment = cell::Segment;

  public:
    template <typename Container>
    class SegmentStore {
      public:
        class Iterator {
          public:
            using T = typename Container::const_iterator;

            Iterator(const T& it)
                : it_(it) {}
            Segment operator*() const {
                return Segment(**it_, *(*it_ + 1));
            }
            bool operator==(const Iterator& o) const {
                return it_ == o.it_;
            }
            bool operator!=(const Iterator& o) const {
                return it_ != o.it_;
            }
            Iterator& operator++() {
                ++it_;
                return *this;
            }

          private:
            T it_;
        };

        SegmentStore(Container&& c)
            : store_(c){};

        Iterator begin() const {
            return Iterator(store_.begin());
        }

        Iterator end() const {
            return Iterator(store_.end());
        }

        size_t size() const {
            return store_.size();
        }

        Segment operator[](size_t idx) const {
            return Segment(*store_[idx], *(store_[idx] + 1));
        }

      private:
        const Container store_;
    };

    using BoundingBox = typename SectionPoint::BoundingBox;
    using SpatialIndex = SpatialGridHash<SectionPoint, t_coordinate>;
    using Segments = SegmentStore<SpatialIndex::QueryResult>;

    SlicedCells(MPI_Comm comm, t_coordinate gridSpacing = 8.0)
        : slices_(mpi::size(comm))
        , index_(gridSpacing)
        , comm_(comm)
        , comm_size_(mpi::size(comm)) {}

    SlicedCells(const SlicedCells&) = delete;
    SlicedCells& operator=(const SlicedCells&) = delete;

    SlicedCells(SlicedCells&&) = default;
    SlicedCells& operator=(SlicedCells&&) = default;

    void createIndex(std::vector<SectionPoint>&& data) {
        index_.createIndex(std::move(data));
    }

    std::vector<BoundingBox>& getSlices() {
        return slices_;
    }

    const std::vector<BoundingBox>& getSlices() const {
        return slices_;
    }

    /**
     * Get the bounding box of a rank
     * \param rank the MPI rank
     */
    const BoundingBox& slice(size_t rank) const {
        return slices_[rank];
    }

    /**
     * Query local segment storage for potential overlaps with another
     * segment
     * \param seg a segment to overlap
     */
    Segments segmentsOverlapping(const Segment& seg) const {
        return Segments(index_.query(seg.boundingBox()));
    }

    /**
     * Get cells to overlapping the current rank's local spatial index
     */
    std::vector<int> cellsOverlapping(const cell::Group&) const;

  private:
    std::vector<BoundingBox> slices_;
    SpatialIndex index_;

    MPI_Comm comm_;
    int comm_size_;
};

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <unordered_map>

#include <morphio/collection.h>

#include <nlohmann/json.hpp>

#include "defines.h"
#include "section.h"

namespace cell {

// class Section;
class SectionPoint;

/** \brief that defines a data related to cell (read from blueMicroCircuitry.bin file) */
struct Metadata {
    std::array<char, 256> cellName;    // TODO: to save memory this should be a key to a uniquely
                                         // identifier name (like a DB secondary key)
    std::array<double, 4> positionData;  // xyz + 4th coordinate for historical reasons
    std::array<double, 4> rotationData;  // center: xyz y-rotation
    int cellId;                        // unique ID starting from 0
    int populationId;
    int morphology;              // The group type of the Cell
    unsigned int segmentsCount;  // total Segments Count
    unsigned int sectionsCount;  // sections count
    bool hasRotations;           // Whether to contain rotation data

    Metadata() = default;
    ~Metadata() = default;

    // Full constructor, this should be the one used!
    template <typename T>
    Metadata(const std::string& name,
             const T& pos,
             const T& rot,
             int id,
             int p,
             int m,
             unsigned int sCount = 0,
             unsigned int bCount = 0,
             bool r = true)
        : positionData()
        , cellId(id)
        , populationId(p)
        , morphology(m)
        , segmentsCount(sCount)
        , sectionsCount(bCount)
        , hasRotations(r) {
        assert(name.size() < 256);
        std::copy(name.begin(), name.end(), cellName.begin());
        cellName[name.size()] = 0;
        std::copy(pos.begin(), pos.end(), positionData.begin());
        std::copy(rot.begin(), rot.end(), rotationData.begin());
    }

    std::string name() const {
        return std::string(cellName.begin());
    }
};

void to_json(nlohmann::json& j, const Metadata& m);
void from_json(const nlohmann::json& j, Metadata& m);

/** \brief Creates the cell properties such as sections*/
class Cell {
  public:
    class const_iterator {
      public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = Section;
        using difference_type = size_t;
        using pointer = Section*;
        using reference = Section&;

        explicit const_iterator(const std::vector<Section>::const_iterator& it,
                                const std::vector<Section>::const_iterator& end)
            : it_(it)
            , end_(end) {}

        const_iterator operator++() {
            for (++it_; it_ != end_; ++it_) {
                if (!it_->empty()) {
                    break;
                }
            }
            return *this;
        }
        const Section& operator*() const {
            return *it_;
        }
        const Section* operator->() {
            return it_.operator->();
        }
        bool operator==(const const_iterator& other) const {
            return it_ == other.it_;
        }
        bool operator!=(const const_iterator& other) const {
            return it_ != other.it_;
        }

      private:
        std::vector<Section>::const_iterator it_;
        std::vector<Section>::const_iterator end_;
    };

    Cell(const Metadata& meta,
           morphio::Collection& storage,
           t_coordinate spine,
           std::function<bool(const SectionType&)> accept)
        : metadata_(meta)
        , sections_(loadMorphology(meta, storage, spine, accept)){};

    const_iterator begin() const {
        for (auto it = sections_.cbegin(); it != sections_.cend(); ++it) {
            if (!(it->empty()))
                return const_iterator(it, sections_.cend());
        }
        return const_iterator(sections_.cend(), sections_.cend());
    };

    const_iterator end() const {
        return const_iterator(sections_.cend(), sections_.cend());
    };

    /**
     * Gets cell morphological type
     */
    int mtype() const {
        return metadata_.morphology;
    }

    /**
     * gets cell GId - global cell index
     * \return
     */
    unsigned int id() const {
        return metadata_.cellId;
    }

    //! Access section
    const Section& operator[](size_t id) const {
        return sections_[id];
    }

    /**
     * returns the number of sections of this cell
     * \return
     */
    size_t size() const {
        return sections_.size();
    }

  private:
    //! Load morphology from a file
    static auto loadMorphology(const Metadata& meta,
                               morphio::Collection& storage,
                               t_coordinate spine,
                               std::function<bool(const SectionType&)> accept)
        -> std::vector<Section>;
    const Metadata metadata_;
    const std::vector<Section> sections_;
};

}  // namespace cell

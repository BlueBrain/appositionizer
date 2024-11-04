#pragma once

#include <nlohmann/json.hpp>

#include "better_mpi.h"

namespace cell {

/**
 * \brief Helper enum to provide constants for the available section types
 *
 * Note that in the hdf5 files, the numbers are 1 higher, so those values
 * are reduced by 1 before being stored here.
 */
enum class SectionType { soma = 0, axon = 1, dend = 2, apic = 3 };

/**
 * \brief Identifier for a cell segment
 *
 * Will store the following information:
 *   - cell index
 *   - section index within the cell
 *   - segment index along the section
 *   - the type of the section that the segment is on
 *   - population index of the cell
 */
class SegmentId {
  public:
    using SectionId = size_t;

    SegmentId() = default;
    SegmentId(size_t n, size_t b, size_t s, size_t t, size_t p = 0)
        : cell_(n)
        , section_(b)
        , segment_(s)
        , section_type_(t)
        , population_(p) {
        if (n > maxCellIndex or b > maxSectionIndex or s > maxSegmentIndex or
            t > maxSectionType or p > maxPopulationIndex) {
            mpi::abort(
                "Invalid segment ID values (maximum):\n"
                "\tcell      {} ({})\n"
                "\tsection      {} ({})\n"
                "\tsegment     {} ({})\n"
                "\tsection type {} ({})\n",
                "\tpopulation  {} ({})\n",
                n,
                maxCellIndex,
                b,
                maxSectionIndex,
                s,
                maxSegmentIndex,
                t,
                maxSectionType,
                p,
                maxPopulationIndex);
        }
    };
    SegmentId(const SegmentId&) = default;
    SegmentId(SegmentId&&) = default;
    SegmentId& operator=(const SegmentId&) = default;
    SegmentId& operator=(SegmentId&&) = default;

    inline bool operator==(const SegmentId& segmentId) {
        return (cell_ == segmentId.cell()) & (section_ == segmentId.section()) &
               (segment_ == segmentId.segment());
    }

    inline bool operator<(const SegmentId& sid) {
        auto lhs = std::tuple{cell(), section(), segment()};
        auto rhs = std::tuple{sid.cell(), sid.section(), sid.segment()};

        return lhs < rhs;
    }

    inline unsigned int population() const {
        return population_;
    }

    inline unsigned int cell() const {
        return cell_;
    }

    inline void cell(size_t n) {
        if (n > maxCellIndex) {
            mpi::abort(
                "Invalid segment ID value (maximum):\n"
                "\tcell      {} ({})\n",
                n,
                maxCellIndex);
        }
        cell_ = n;
    }

    inline unsigned int section() const {
        return section_;
    }

    inline void section(size_t b) {
        if (b > maxSectionIndex) {
            mpi::abort(
                "Invalid segment ID value (maximum):\n"
                "\tsection      {} ({})\n",
                b,
                maxSectionIndex);
        }
        section_ = b;
    }

    inline unsigned int segment() const {
        return segment_;
    }

    inline void segment(size_t s) {
        if (s > maxSegmentIndex) {
            mpi::abort(
                "Invalid segment ID value (maximum):\n"
                "\tsegment     {} ({})\n",
                s,
                maxSegmentIndex);
        }
        segment_ = s;
    }

    inline unsigned int sectionType() const {
        return section_type_;
    }

    inline void sectionType(size_t t) {
        if (t > maxSectionType) {
            mpi::abort(
                "Invalid segment ID value (maximum):\n"
                "\tsection type {} ({})\n",
                t,
                maxSectionType);
        }
        section_type_ = t;
    }

    inline bool sameSection(const SegmentId& other) {
        return population_ == other.population_ && cell_ == other.cell_ &&
               section_ == other.section_;
    }

    /**
     * \brief Return a compact representation of the Id with 14 bit parts for section/segment index
     * 
     * \todo Change unit tests after dropping population and section type to allow for 16 bit section/segment
     */
    inline std::uint64_t compact() const {
        // return (static_cast<std::uint64_t>(cell_) << 32) | (static_cast<std::uint64_t>(section_) << 18) | (static_cast<std::uint64_t>(segment_) << 14)
        //     | (static_cast<std::uint64_t>(section_type_) << 2) | (static_cast<std::uint64_t>(population_) << 1);
        return (static_cast<std::uint64_t>(cell_) << 0) | (static_cast<std::uint64_t>(section_) << 32) | (static_cast<std::uint64_t>(segment_) << 46)
            | (static_cast<std::uint64_t>(section_type_) << 60) | (static_cast<std::uint64_t>(population_) << 62);

    }

  private:
    std::uint32_t cell_;  // 2^32 = 4B (human brain has 100B)
    std::uint16_t section_;
    std::uint16_t segment_;
    std::uint8_t section_type_;
    std::uint8_t population_;

    static constexpr size_t maxCellIndex = std::numeric_limits<decltype(cell_)>::max();
    static constexpr size_t maxSectionIndex = std::numeric_limits<decltype(section_)>::max();
    static constexpr size_t maxSegmentIndex = std::numeric_limits<decltype(segment_)>::max();
    static constexpr size_t maxPopulationIndex = 1;  // we only handle two cell populations
    static constexpr size_t maxSectionType = 3; // see enum above
};

static_assert(std::is_trivial<SegmentId>::value, "SegmentId should be trivial");
static_assert(std::is_trivially_copyable<SegmentId>::value,
              "SegmentId should be trivially copyable");

}  // namespace cell

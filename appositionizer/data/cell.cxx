#include "cell.h"  // Must be on the first substantial line.

#include <cassert>
#include <cstring>
#include <functional>

#include <morphio/errorMessages.h>
#include <morphio/morphology.h>
#include <morphio/section.h>
#include <morphio/soma.h>

#include <Eigen/Geometry>

#include "../data/segment.h"

using nlohmann::json;

namespace {

template <typename T, typename U>
void rotate(Eigen::Matrix3Xd& segments, const T& cell, const U& quaternion) {
    Eigen::Quaterniond q{quaternion[3], quaternion[0], quaternion[1], quaternion[2]};
    q.normalize();
    Eigen::Translation3d t{cell[0], cell[1], cell[2]};

    segments = (t * q) * segments;
}

template <typename T>
t_coordinate distance(const T& a, const T& b) {
    t_coordinate result = 0.;
    for (size_t i = 0; i < 3; ++i) {
        result += (b[i] - a[i]) * (b[i] - a[i]);
    }
    return std::sqrt(result);
}

/// Returns the distance along sections, indexed by section id excluding the soma.
template <typename T>
std::vector<std::vector<T>> pathlengths(const morphio::Morphology& m) {
    std::vector<std::vector<T>> result;
    result.resize(m.sections().size());

    for (const auto& section: m.sections()) {
        const auto idx = section.id();
        const auto& points = section.points();
        if (points.size() > 0) {
            result[idx].resize(points.size());
            result[idx][0] = 0.0;
            for (size_t i = 1; i < section.points().size(); ++i) {
                result[idx][i] = distance(points[i - 1], points[i]) + result[idx][i - 1];
            }
        }
    }

    return result;
}

/// Returns the distances of sections to the soma, given a morphology and corresponding
/// pathlengths.
template <typename T>
std::vector<T> distances_to_soma(const morphio::Morphology& m,
                                 const std::vector<std::vector<T>>& pathlengths) {
    std::vector<T> result;
    result.resize(m.sections().size());

    for (const auto& section: m.sections()) {
        if (!section.isRoot()) {
            auto pid = section.parent().id();
            result[section.id()] = result[pid] + pathlengths[pid].back();
        }
    }

    return result;
}

}  // namespace

namespace cell {

auto Cell::loadMorphology(const Metadata& meta,
                            morphio::Collection& storage,
                            t_coordinate spine,
                            std::function<bool(const SectionType&)> accept) -> std::vector<Section> {
    morphio::set_ignored_warning(morphio::Warning::ONLY_CHILD, true);

    std::vector<Section> sections;

    auto key =
        SegmentId(meta.cellId, 0, 0, static_cast<int>(SectionType::soma), meta.populationId);

    const auto& morpho = storage.load<morphio::Morphology>(meta.name());

    sections.reserve(morpho.sections().size() + 1);

    ssize_t last_section_type = -1;

    // Add the soma: two points for historical reasons.
    if (accept(SectionType::soma)) {
        auto radius = morpho.soma().maxDistance();

        std::vector<SectionPoint> segs;
        segs.emplace_back(meta.positionData, 0, 0, radius, 0, key);
        key.segment(1);
        segs.emplace_back(meta.positionData, 0, 0, radius, 0, key);
        sections.emplace_back(-1, 0., std::move(segs));
        last_section_type = 0;
    } else {
        // Put empty sections in the vector, because we keep index of position
        sections.emplace_back(0, 0., std::vector<SectionPoint>{});
    }

    const auto& section_pathlengths = pathlengths<t_coordinate>(morpho);
    const auto& section_distances_to_soma = distances_to_soma<t_coordinate>(morpho,
                                                                            section_pathlengths);

    for (const auto& section: morpho.sections()) {
        auto section_type = section.type() - 1;
        if (section_type < last_section_type) {
            throw std::runtime_error("Sections are out of order in " + meta.name());
        }
        last_section_type = section_type;

        if (!accept(static_cast<SectionType>(section_type))) {
            // Put empty sections in the vector, because we keep index of position
            sections.emplace_back(0, 0., std::vector<SectionPoint>{});
            continue;
        }

        const size_t sectionSize = section.points().size();
        if (sectionSize < 2) {
            mpi::abort("Section {} of morphology {} has a section with only {} segment(s)",
                       section.id() + 1,
                       meta.name(),
                       sectionSize);
        }

        key.section(section.id() + 1);
        key.sectionType(section_type);

        bool use_spine = key.sectionType() == static_cast<int>(SectionType::apic) or
                         key.sectionType() == static_cast<int>(SectionType::dend);

        // Put all segments in one matrix and rotate all the segment in-place in one matrix
        // operation
        const auto& points = section.points();
        Eigen::Matrix3Xd _points(3, sectionSize);
        if (meta.hasRotations) {
            for (size_t i = 0; i < sectionSize; ++i)
                _points.col(i) = Eigen::Vector3d(points[i][0], points[i][1], points[i][2]);
            rotate(_points, meta.positionData, meta.rotationData);
        } else {
            auto shift =
                Eigen::Vector3d(meta.positionData[0], meta.positionData[1], meta.positionData[2]);
            for (size_t i = 0; i < sectionSize; ++i)
                _points.col(i) = Eigen::Vector3d(points[i][0], points[i][1], points[i][2]) + shift;
        }

        const auto& diameters = section.diameters();
        const auto& positions = section_pathlengths[section.id()];
        const auto& section_length = positions.empty() ? 0.0 : positions.back();

        std::vector<SectionPoint> segments;
        segments.reserve(sectionSize);
        for (size_t i = 0; i < sectionSize; ++i) {
            key.segment(i);
            segments.emplace_back(_points.col(i),
                                  positions[i],
                                  section_length,
                                  diameters[i] * 0.5,
                                  use_spine ? spine : 0.0,
                                  key);
        }

        int order = 1;
        for (auto i = section.upstream_begin(); i != section.upstream_end(); ++i) {
            ++order;
        }

        sections.emplace_back(order, section_distances_to_soma[section.id()], std::move(segments));
    }
    return sections;
}

void to_json(nlohmann::json& j, const cell::Metadata& m) {
    j = json{{"name", m.cellName},
             {"position", m.positionData},
             {"rotation", m.rotationData},
             {"id", m.cellId},
             {"morphology", m.morphology},
             {"segments", m.segmentsCount},
             {"sections", m.sectionsCount}};
}

void from_json(const nlohmann::json& j, cell::Metadata& m) {
    j.at("name").get_to(m.cellName);
    j.at("position").get_to(m.positionData);
    j.at("rotation").get_to(m.rotationData);
    j.at("id").get_to(m.cellId);
    j.at("morphology").get_to(m.morphology);
    j.at("segments").get_to(m.segmentsCount);
    j.at("sections").get_to(m.sectionsCount);

    // Put default value to avoid unintialized values
    m.populationId = 0;
    m.hasRotations = false;
}

}  // namespace cell

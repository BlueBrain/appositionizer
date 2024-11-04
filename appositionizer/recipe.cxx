#include "recipe.h"  // Must be on the first substantial line.

#include <yaml-cpp/yaml.h>


auto read_recipe(const std::string& filename) -> Recipe {
    Recipe result;

    const auto doc = YAML::LoadFile(filename);

    if (!doc["bouton_interval"].IsDefined()) {
        throw std::runtime_error("No 'bouton_interval' key in " + filename);
    }

    try {
        result.interval.minDistance = doc["bouton_interval"]["min_distance"].as<double>();
        result.interval.maxDistance = doc["bouton_interval"]["max_distance"].as<double>();
        result.interval.regionGap = doc["bouton_interval"]["region_gap"].as<double>();
    } catch (const std::exception& e) {
        throw std::runtime_error("Unparseable 'bouton_interval' in " + filename + ": " + e.what());
    }

    const auto nodes = doc["structural_spine_lengths"];
    for (const auto& res: nodes) {
        float length = res["spine_length"].as<float>();
        std::string mtype = res["mtype"].as<std::string>();

        auto match = result.spines.find(mtype);
        if (match == result.spines.end()) {
            result.spines.emplace(mtype, length);
        } else if (match->second != length) {
            throw std::runtime_error("Conflicting spine lengths for " + mtype + " in " + filename);
        }
    }

    return result;
}

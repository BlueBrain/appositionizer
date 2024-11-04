#pragma once

#include <string>
#include <unordered_map>

using SpineLengths = std::unordered_map<std::string, float>;

struct InterBoutonInterval {
    double minDistance;
    double maxDistance;
    double regionGap;
};

struct Recipe {
    InterBoutonInterval interval;
    SpineLengths spines;
};

Recipe read_recipe(const std::string& filename);

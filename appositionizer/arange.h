#pragma once

#include <vector>

template <class Int>
std::vector<Int> arange(Int a, Int b) {
    if (b < a) {
        throw std::runtime_error("Detected: b < a.");
    }

    size_t n = b - a;
    auto x = std::vector<Int>(n);
    for (size_t i = 0; i < n; ++i) {
        x[i] = a + Int(i);
    }

    return x;
}

template <class Int>
std::vector<Int> arange(size_t n) {
    return arange(Int(0), Int(n));
}

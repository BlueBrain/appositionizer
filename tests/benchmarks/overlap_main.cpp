#include <array>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "../../appositionizer/partial_collision.h"
#include "../../appositionizer/Segment.h"

using cell::Segment;
using cell::SegmentId;


inline static std::vector<Segment> create_segments(size_t num, float r) {
    std::vector<Segment> result(num);

    std::mt19937 rng(123);
    std::uniform_real_distribution<float> xyz(0.0f, 100.0f);
    std::uniform_real_distribution<float> radius(0.0f, r);
    std::uniform_real_distribution<float> spine(0.0f, 1.0f);

    for (size_t i = 0; i < num; ++i) {
        Coord3D<t_coordinate> p1{xyz(rng), xyz(rng), xyz(rng)};
        Coord3D<t_coordinate> p2{xyz(rng), xyz(rng), xyz(rng)};
        result[i] = Segment(SegmentId(0, 0, 0, 2), p1, p2, radius(rng), 0, spine(rng));
    }

    return result;
}


static const auto SEGMENT_COUNT = 200;


static void overlap(Segment::Mode mode, float r) {
    auto SEGMENTS = create_segments(SEGMENT_COUNT, r);
    for (const auto& s1: SEGMENTS) {
        for (const auto& s2: SEGMENTS) {
            s1.overlaps(s2, mode);
        }
    }
}


static void detection__naive(benchmark::State& state) {
    for (auto _: state) {
        overlap(Segment::Mode::iterative, state.range(0) * 0.1f);
    }
}
BENCHMARK(detection__naive)->RangeMultiplier(2)->Range(1, 128);


static void detection__fixed(benchmark::State& state) {
    for (auto _: state) {
        overlap(Segment::Mode::sampled, state.range(0) * 0.1f);
    }
}
BENCHMARK(detection__fixed)->RangeMultiplier(2)->Range(1, 128);


static void detection__random(benchmark::State& state) {
    for (auto _: state) {
        overlap(Segment::Mode::random, state.range(0) * 0.1f);
    }
}
BENCHMARK(detection__random)->RangeMultiplier(2)->Range(1, 128);


BENCHMARK_MAIN();

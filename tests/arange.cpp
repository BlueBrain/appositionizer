#include <catch2/catch_test_macros.hpp>

#include "../appositionizer/arange.h"

TEST_CASE("arange(0)") {
    size_t n = 0;
    auto x = arange<size_t>(n);

    CHECK(x.size() == 0);
}

TEST_CASE("arange(n)") {
    size_t n = 5;
    auto x = arange<size_t>(n);

    CHECK(x.size() == n);
    for (size_t i = 0; i < n; ++i) {
        REQUIRE(x[i] == i);
    }
}

TEST_CASE("arange(a, n)") {
    int a = -3;
    int b = 2;
    size_t n = 5;
    auto x = arange<int>(a, b);

    CHECK(x.size() == n);
    for (size_t i = 0; i < n; ++i) {
        REQUIRE(x[i] == a + int(i));
    }
}

#include <catch2/catch_test_macros.hpp>

#include "../appositionizer/index_magic.h"

TEST_CASE("div_up") {
    CHECK(index_magic::div_up(0, 1) == 0);
    CHECK(index_magic::div_up(1, 1) == 1);
    CHECK(index_magic::div_up(1, 4) == 1);
    CHECK(index_magic::div_up(1, 2) == 1);
    CHECK(index_magic::div_up(1, 2) == 1);
    CHECK(index_magic::div_up(3, 2) == 2);
    CHECK(index_magic::div_up(7, 3) == 3);
}

TEST_CASE("fair_share::element2bucket") {
    SECTION("regular") {
        size_t n_elements = 7;
        size_t n_buckets = 3;

        // | 0 , 1 , 2 | 3 , 4 | 5 , 6 |
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 0) == 0);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 1) == 0);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 2) == 0);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 3) == 1);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 4) == 1);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 5) == 2);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 6) == 2);
    }

    SECTION("partially empty") {
        size_t n_elements = 3;
        size_t n_buckets = 7;

        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 0) == 0);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 1) == 1);
        CHECK(index_magic::fair_share::element2bucket(n_elements, n_buckets, 2) == 2);
    }
}

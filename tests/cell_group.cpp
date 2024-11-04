#include "../appositionizer/data/cell_group.h"
#include <catch2/catch_test_macros.hpp>

class TestGroup: public cell::Group {
  public:
    using cell::Group::indexPopulation;
};

TEST_CASE("check population ID assignment") {
    const auto id1 = TestGroup::indexPopulation("foo.h5", "bar");
    const auto id2 = TestGroup::indexPopulation("foo.h5", "spam");
    REQUIRE(id2 != id1);
    const auto id3 = TestGroup::indexPopulation("pub.h5", "bar");
    REQUIRE((id3 != id1 && id3 != id2));
    const auto id4 = TestGroup::indexPopulation("foo.h5", "bar");
    REQUIRE(id4 == id1);
}

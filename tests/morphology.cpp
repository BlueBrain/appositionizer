#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "../appositionizer/data/cell_group.h"
#include "../appositionizer/process/apposition_filter.h"

namespace fs = std::filesystem;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using nlohmann::json;

using cell::Metadata;
using apposition::Region;
using apposition::Apposition;

static auto datadir = fs::path(__FILE__).parent_path() / "data/";
static auto circuit_path = (datadir / "nodes.h5").string();
static auto recipe = read_recipe((datadir / "recipe.yaml").string());

static MPI_Comm dummy_mpi_comm() {
    return MPI_COMM_WORLD;
}

const cell::Group& cells() {
    static auto g =
        cell::Group(circuit_path, "All", datadir.string(), recipe.spines, dummy_mpi_comm());
    return g;
}

TEST_CASE("Verify detection of out-of-order sections") {
    Metadata m("sections_out_of_order",
               std::vector<double>{1, 1, 1},
               std::vector<double>{1, 0, 0, 0},
               1,
               0,
               0,
               0);

    SECTION("Catch exception upon loading") {
        cell::Group g(cells(), {m});
        REQUIRE_THROWS(g.loadMorphologies());
    }

    SECTION("Catch exception upon loading, skipping the axon") {
        cell::Group g(cells(), {m});
        REQUIRE_THROWS(g.loadMorphologies(
            [](const cell::SectionType& type) { return type != cell::SectionType::axon; }));
    }
}

TEST_CASE("Verify morphology loading") {
    Metadata m1("soma_spherical",
                std::vector<double>{1, 1, 1},
                std::vector<double>{1, 0, 0, 0},
                1,
                0,
                0,
                0);

    Metadata m2(
        "soma_extended",
        std::vector<double>{1, .5, .25},
        std::vector<double>{0.8728715609439696, 0.4364357804719848, 0.2182178902359924, 0.0},
        2,
        0,
        0,
        0);

    cell::Group g(cells(), {m1, m2});
    g.loadMorphologies();

    SECTION("Cell translated correctly") {
        auto n1 = g[1];

        auto soma = n1[0].raw(0);
        REQUIRE_THAT(soma.xyz()[0], WithinAbs(1, 0.00001));
        REQUIRE_THAT(soma.xyz()[1], WithinAbs(1, 0.00001));
        REQUIRE_THAT(soma.xyz()[2], WithinAbs(1, 0.00001));
        REQUIRE_THAT(soma.radius(), WithinAbs(0.18028f, 0.00001));

        auto random_segment = n1[3].raw(12);
        REQUIRE_THAT(random_segment.xyz()[0], WithinAbs(7.98622f, 0.00001));
        REQUIRE_THAT(random_segment.xyz()[1], WithinAbs(13.17931f, 0.00001));
        REQUIRE_THAT(random_segment.xyz()[2], WithinAbs(18.53813f, 0.00001));
        REQUIRE_THAT(random_segment.radius(), WithinAbs(0.3f, 0.00001));
    }

    SECTION("Cell translated and rotated correctly") {
        auto n2 = g[2];

        auto soma = n2[0].raw(0);
        REQUIRE_THAT(soma.xyz()[0], WithinAbs(1, 0.00001));
        REQUIRE_THAT(soma.xyz()[1], WithinAbs(.5, 0.00001));
        REQUIRE_THAT(soma.xyz()[2], WithinAbs(.25, 0.00001));
        REQUIRE_THAT(soma.radius(), WithinAbs(4.41688, 0.00001));

        auto random_segment = n2[3].raw(12);
        REQUIRE_THAT(random_segment.xyz()[0], WithinAbs(25.1764f, 0.00001));
        REQUIRE_THAT(random_segment.xyz()[1], WithinAbs(3.07876f, 0.00001));
        REQUIRE_THAT(random_segment.xyz()[2], WithinAbs(34.16223f, 0.00001));
        REQUIRE_THAT(random_segment.radius(), WithinAbs(2.8f, 0.00001));
    }
}

TEST_CASE("Verify distances within morphologies") {
    Metadata m1("soma_extended",
                std::vector<double>{0, 0, 0},
                std::vector<double>{1, 0, 0, 0},
                1,
                0,
                8,
                1751);
    Metadata m2("soma_spherical",
                std::vector<double>{0, 0, 0},
                std::vector<double>{1, 0, 0, 0},
                2,
                0,
                8,
                1751);

    cell::Group g(cells(), {m2, m1});
    g.loadMorphologies();
    apposition::Filter f(g[1][11], {0.1, 0.2, 0.5});

    SECTION("Segment distance") {
        REQUIRE_THAT(f.getDistance(1, 11, 2, 123, 2, 234), WithinAbs(0.3553589, 1e-5));
        REQUIRE_THAT(f.getDistance(1, 11, 2, 123, 5, 234), WithinAbs(2.3617870766, 1e-5));
        REQUIRE(f.getDistance(1, 11, 2, 235, 2, 234) <= 0.0);
    }

    SECTION("Raw cell properties") {
        auto s1 = g[1][11][2];
        REQUIRE(s1.cell() > 0);
    }

    auto seg1 = g[1][11][2];
    auto seg2 = g[1][11][3];
    auto seg3 = g[1][11][5];

    auto seg4 = g[2][2][5];

    SECTION("Segment properties") {
        REQUIRE(seg1.cell() > 0);
        REQUIRE(seg2.cell() > 0);
        REQUIRE(seg3.cell() > 0);
        REQUIRE(seg4.cell() > 0);
    }

    Apposition t1(seg1, seg4, {123, 254, 0, 0});
    Apposition t2(seg2, seg4, {001, 123, 0, 0});
    Apposition t3(seg3, seg4, {0, 0, 0, 0});

    SECTION("Apposition properties") {
        REQUIRE(t1.preCell() > 0);
        REQUIRE(t2.preCell() > 0);
        REQUIRE(t3.preCell() > 0);
    }

    SECTION("Apposition region length 1") {
        std::vector<Apposition> appositions{t1, t2};
        auto regions = f.filter(appositions);
        REQUIRE(regions.size() == 1);
    }

    SECTION("Apposition region length 2") {
        std::vector<Apposition> appositions{t1, t2, t3};
        auto regions = f.filter(appositions);
        REQUIRE(regions.size() == 2);

        auto& r = regions[0];

        std::vector<short int> cmp;
        for (auto& c: r.connections) {
            cmp.push_back(c.pre.segment);
        }

        std::vector<short int> segs{2, 2, 3, 3};
        REQUIRE(cmp == segs);
    }
}

TEST_CASE("Real-world Apposition to region") {
    auto basedir = fs::path(__FILE__).parent_path() / "data" / "regions";
    auto morpho_path = (basedir / "h5/").string();
    auto group = cell::Group(circuit_path, "All", morpho_path, recipe.spines, dummy_mpi_comm());

    std::vector<std::string> subdirs{"case_01", "case_02", "case_03"};
    for (auto dir: subdirs) {
        auto data_subdir = basedir / dir;

        std::vector<Metadata> meta;
        {
            std::ifstream i((data_subdir / "metadata.json").string());
            json j;
            i >> j;
            for (auto& e: j)
                meta.push_back(e);
        }

        std::vector<Apposition> appositions;
        {
            std::ifstream i((data_subdir / "appositions.json").string());
            json j;
            i >> j;
            for (auto& e: j)
                appositions.push_back(e);
        }

        std::vector<Region> regions;
        {
            std::ifstream i((data_subdir / "regions.json").string());
            json j;
            i >> j;
            for (auto& e: j)
                regions.push_back(e);
        }

        cell::Group g(group, meta);
        g.loadMorphologies();

        DYNAMIC_SECTION("Compare outcome for " + dir) {
            std::vector<apposition::Region> result;

            // Backwards compatible shim:
            //
            // appositions used to be filtered for a batch of cells. appositions
            // for this test are loaded as such.  Break them up into
            // appositions per section, and concatenate the result.
            auto it = appositions.begin();
            while (it != appositions.end()) {
                auto next = it + 1;
                for (; next < appositions.end() && next->preCell() == it->preCell() &&
                       next->preSection() == it->preSection();
                     ++next) {
                }
                std::vector<Apposition> slice(it, next);
                apposition::Filter f(g[it->preCell()][it->preSection()], recipe.interval);
                auto rs = f.filter(slice);
                std::copy(rs.begin(), rs.end(), std::back_inserter(result));
                it = next;
            }

            for (size_t i = 0; i < std::min(regions.size(), result.size()); ++i) {
                // exclude Apposition indices from copmarison
                regions[i].resetIndex(result[i].getStartAppositionIndex());

                const auto& rref = regions[i];
                const auto& rcmp = result[i];

                for (size_t j = 0; j < std::min(rref.connections.size(), rcmp.connections.size());
                     ++j) {
                    json cref = rref.connections[j];
                    json ccmp = rcmp.connections[j];

                    REQUIRE(cref.dump(2) == ccmp.dump(2));
                }

                REQUIRE(rref.connections.size() == rcmp.connections.size());

                json ref = rref;
                json cmp = rcmp;

                ref.erase("length");
                cmp.erase("length");

                ref.erase("connections");
                cmp.erase("connections");

                REQUIRE(ref.dump(2) == cmp.dump(2));
            }

            REQUIRE(result.size() == regions.size());
        }
    }
}

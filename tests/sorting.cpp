#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <vector>

#include "../appositionizer/data/apposition.h"
#include "../appositionizer/io/filewriter.h"
#include "../appositionizer/io/sorting.h"

namespace fs = std::filesystem;


TEST_CASE("sort_into_buckets; example") {
    auto input_dirname = "tmp.input";
    auto input_filename = io::join_path(input_dirname, "edgesData.0");
    fs::create_directories(input_dirname);

    auto output_dirname = "tmp.output";
    auto output_filename = io::join_path(output_dirname, "edgesData.0");
    fs::create_directories(output_dirname);

    auto bucket_sort_dirname = "tmp.bucket_sort";
    auto bucket_sort_filename = io::join_path(bucket_sort_dirname, "edgesData.0");
    fs::create_directories(bucket_sort_dirname);

    auto make_entry = [](size_t pre_gid, size_t post_gid, size_t unique) {
        io::DataEntry entry;
        entry.preCell = pre_gid;
        entry.postCell = post_gid;
        entry.preSegment = unique;

        return entry;
    };

    auto are_equal = [](const io::DataEntry& a, const io::DataEntry& b) {
        return a.preCell == b.preCell && a.postCell == b.postCell &&
               a.preSegment == b.preSegment;
    };

    std::vector<io::DataEntry> input_entries;
    // GID 5: 3 input_entries
    input_entries.push_back(make_entry(5, 0, 500));
    input_entries.push_back(make_entry(5, 1, 510));
    input_entries.push_back(make_entry(5, 1, 511));
    // GID 1: 5 input_entries
    input_entries.push_back(make_entry(1, 0, 100));
    input_entries.push_back(make_entry(1, 1, 111));
    input_entries.push_back(make_entry(1, 2, 120));
    input_entries.push_back(make_entry(1, 1, 110));
    input_entries.push_back(make_entry(1, 3, 110));
    // GID 6: 4 input_entries
    input_entries.push_back(make_entry(6, 0, 600));
    input_entries.push_back(make_entry(6, 1, 611));
    input_entries.push_back(make_entry(6, 1, 610));
    input_entries.push_back(make_entry(6, 3, 610));
    // GID 4: 4 input_entries
    input_entries.push_back(make_entry(4, 0, 400));
    input_entries.push_back(make_entry(4, 1, 411));
    input_entries.push_back(make_entry(4, 1, 410));
    input_entries.push_back(make_entry(4, 3, 410));

    auto gid2bid = [](size_t gid) { return gid / 2; };
    auto bucket_counts = std::vector<size_t>{
        5,  // GIDs 0, 1
        0,  // GIDs 2, 3
        7,  // GIDs 4, 5
        4   // GIDs 6, 7
    };
    auto bucket_boundaries = std::vector<size_t>{0, 5, 5, 12, 16};
    auto gid_boundaries = std::vector<size_t>{0, 0, 5, 5, 5, 9, 12, 16};
    auto input_index =
        std::vector<io::IndexEntry>{io::IndexEntry(5, 3, 0 * sizeof(io::DataEntry)),
                                    io::IndexEntry(1, 5, 3 * sizeof(io::DataEntry)),
                                    io::IndexEntry(6, 4, 8 * sizeof(io::DataEntry)),
                                    io::IndexEntry(4, 4, 12 * sizeof(io::DataEntry))};

    {
        auto writer = io::FileWriter(input_dirname, 0ul, 0ul, MPI_COMM_SELF);
        writer.write(input_entries, input_index);
    }

    io::sort_into_buckets(output_filename, input_filename, bucket_counts, gid2bid);

    auto load_entries = [&input_entries](const std::string& filename) {
        auto file = mpi::File::open(filename, mpi::File::Mode::READ, MPI_COMM_SELF);
        auto entries = std::vector<io::DataEntry>(input_entries.size());
        file.read_array_at(entries, 0ul);

        return entries;
    };

    auto bucketed_entries = load_entries(output_filename);

    for (const auto& ie: input_entries) {
        bool found = false;
        for (size_t i = 0; i < bucketed_entries.size(); ++i) {
            const auto& be = bucketed_entries[i];
            if (are_equal(ie, be)) {
                auto gid = be.preCell;
                auto bid = gid2bid(gid);

                REQUIRE(bucket_boundaries[bid] <= i);
                REQUIRE(i < bucket_boundaries[bid + 1]);

                found = true;
            }
        }

        REQUIRE(found);
    }

    std::sort(input_entries.begin(), input_entries.end());

    io::sort_buckets(output_filename, bucket_counts);
    auto sorted_entries = load_entries(output_filename);

    for (size_t i = 0; i < input_entries.size(); ++i) {
        const auto& ie = input_entries[i];
        const auto& se = sorted_entries[i];

        REQUIRE(are_equal(ie, se));
    }

    auto index = create_index_for_sorted(sorted_entries);
    REQUIRE(index.size() == 4);
    for (const auto& ie: index) {
        auto gid = ie.cellId;
        auto count = gid_boundaries[gid + 1] - gid_boundaries[gid];
        auto offset = gid_boundaries[gid] * sizeof(io::DataEntry);

        REQUIRE(count == ie.synapseCount);
        // TODO improve int correctness.
        REQUIRE((long long int) (offset) == ie.connectionsBinaryAddress);
    }

    io::bucket_sort(bucket_sort_dirname, input_dirname, bucket_counts.size(), gid2bid);

    sorted_entries = load_entries(bucket_sort_filename);
    for (size_t i = 0; i < input_entries.size(); ++i) {
        const auto& ie = input_entries[i];
        const auto& se = sorted_entries[i];

        if (!are_equal(ie, se)) {
            std::cout << "i = " << i << ": " << ie.preSegment << " != " << se.preSegment << "\n";
        }

        REQUIRE(are_equal(ie, se));
    }
}

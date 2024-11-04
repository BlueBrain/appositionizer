#include "sorting.h"

#include <numeric>
#include <oneapi/tbb/parallel_sort.h>

#include "bucketed_writer.h"
#include "filereader.h"
#include "filewriter.h"

#include "../better_mpi.h"

namespace io {

void bucket_sort(const std::string& output_dir,
                 const std::string& input_dir,
                 size_t n_buckets,
                 const std::function<size_t(size_t)>& gid2bid) {
    // This method will in a first pass sort the data into buckets. The
    // bucketed files are stored in the `output_dir`. Then the file is
    // sorted bucket by bucket (inplace).

    // TODO avoid bucket sort if the entire file fits into RAM.
    // TODO perform two rounds of bucket sort if the buckets are too big.

    auto phase = Instrumentor::phase("bucket_sort");

    std::map<size_t, IndexEntry> index = read_index(index_filename(input_dir, MPI_COMM_SELF));

    std::vector<size_t> bucket_counts(n_buckets, 0ul);
    for (auto [gid, index_entry]: index) {
        auto bid = gid2bid(gid);
        bucket_counts[bid] += index_entry.synapseCount;
    }

    auto input_filename = io::payload_filename(input_dir, MPI_COMM_SELF);
    auto output_filename = io::payload_filename(output_dir, MPI_COMM_SELF);

    sort_into_buckets(output_filename, input_filename, bucket_counts, gid2bid);
    sort_buckets(output_filename, bucket_counts);

    // The datafiles need to be accompanied by an index file. Since we know the file
    // is sorted by pre-GID we can compute the new index easily from the old index.
    auto sorted_index = std::vector<IndexEntry>{};
    sorted_index.reserve(index.size());

    long long file_offset = 0;
    for (auto [gid, entry]: index) {
        if (entry.synapseCount > 0) {
            sorted_index.emplace_back(gid, entry.synapseCount, file_offset);
            file_offset += entry.synapseCount * sizeof(DataEntry);
        }
    }

    auto index_filename = io::index_filename(output_dir, MPI_COMM_SELF);
    write_index(index_filename, sorted_index);
}

void sort_into_buckets(const std::string& output_filename,
                       const std::string& input_filename,
                       const std::vector<size_t>& bucket_counts,
                       const std::function<size_t(size_t)>& gid2bid) {
    // Description: The BucketedWriter places the entries into their bucket and
    // writes buckets to disk once enough entries have accumulated in each
    // bucket. This allows us to simply load a large chunk of the input file;
    // pass the entries to the bucketed writer; and continue with the next
    // chunk.
    //
    // The point is to read and write in 'large' blocks to avoid cripling
    // performance on parallel filesystems (such as GPFS).

    size_t file_count = std::accumulate(bucket_counts.begin(), bucket_counts.end(), 0ul);
    size_t file_size = file_count * sizeof(DataEntry);

    auto input_file = mpi::File::open(input_filename, mpi::File::Mode::READ, MPI_COMM_SELF);

    constexpr size_t batch_count = 1ul << 20;
    constexpr size_t batch_size = batch_count * sizeof(DataEntry);
    auto n_batches = (file_size + batch_size - 1) / batch_size;

    std::vector<DataEntry> data_entry_buffer(batch_count);
    auto bucketed_writer = BucketedWriter(output_filename, gid2bid, bucket_counts);

    for (size_t k = 0; k < n_batches; ++k) {
        size_t file_begin = k * batch_size;
        size_t file_end = std::min((k + 1) * batch_size, file_size);
        size_t count = (file_end - file_begin) / sizeof(DataEntry);

        Instrumentor::phase_begin("reading");
        input_file.read_array_at(std::span{data_entry_buffer.data(), count}, file_begin);
        Instrumentor::phase_end("reading");

        Instrumentor::phase_begin("writing");
        for (size_t i = 0; i < count; ++i) {
            bucketed_writer.write(data_entry_buffer[i]);
        }
        Instrumentor::phase_end("writing");
    }
}

void sort_buckets(const std::string& filename, const std::vector<size_t> bucket_counts) {
    auto phase = Instrumentor::phase("inplace_sort");

    auto sorted_file = mpi::File::open(filename, mpi::File::Mode::READ_WRITE, MPI_COMM_SELF);

    auto n_buckets = bucket_counts.size();
    auto max_bucket_count = std::reduce(bucket_counts.begin(),
                                        bucket_counts.end(),
                                        0ul,
                                        [](size_t a, size_t b) { return std::max(a, b); });

    auto data_bucket = std::vector<DataEntry>(max_bucket_count);

    size_t file_begin = 0ul;
    for (size_t k = 0; k < n_buckets; ++k) {
        size_t bucket_count = bucket_counts[k];

        sorted_file.read_array_at(std::span{data_bucket.data(), bucket_count}, file_begin);

        Instrumentor::phase_begin("tbb::parallel_sort");
        tbb::parallel_sort(data_bucket.begin(), data_bucket.begin() + bucket_count);
        Instrumentor::phase_end("tbb::parallel_sort");

        sorted_file.write_array_at(std::span{data_bucket.data(), bucket_count}, file_begin);

        file_begin += bucket_count * sizeof(DataEntry);
    }
}

std::vector<IndexEntry> create_index_for_sorted(const std::vector<DataEntry>& data_entries) {
    auto comp = [](size_t value, const auto& element) {  // TODO improve int correctness.
        return int(value) <= element.preCell;
    };

    auto data_begin = data_entries.begin();
    auto it = data_entries.begin();
    auto data_end = data_entries.cend();

    auto index_entries = std::vector<IndexEntry>{};
    while (it != data_end) {
        auto pre_gid = it->preCell;
        auto next_it = std::upper_bound(it, data_end, pre_gid + 1, comp);

        size_t count = next_it - it;
        size_t offset = (it - data_begin) * sizeof(DataEntry);
        index_entries.emplace_back(pre_gid, count, offset);

        it = next_it;
    }

    return index_entries;
}

}  // namespace io

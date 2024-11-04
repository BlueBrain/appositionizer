#pragma once

#include <functional>
#include <string>
#include <vector>

#include "../instrumentor.h"
#include "../better_mpi.h"
#include "structures.h"

namespace io {
/** \brief Utility class for creating a file sorted by buckets.
 *
 *  The twist is that:
 *    a) the data to be sorted into buckets doesn't fit into RAM,
 *    b) the filesystem doesn't allow writing in small blocks.
 *
 *  This is mostly a utility class for `sort_into_buckets`.
 */
class BucketedWriter {
  public:
    /** \brief Create a bucketed writer.
     *
     * Here `filename` is the path of the output file, `gid2bid_callback` maps
     * pre-GID to bucket ID, and `bucket_counts` is the number of appositions per
     * bucket.
     */
    BucketedWriter(const std::string& filename,
                   std::function<size_t(size_t)> gid2bid_callback,
                   const std::vector<size_t>& bucket_counts)
        : buckets(bucket_counts.size())
        , bucket_file_offset(bucket_counts.size())
        , gid2bid_callback_(gid2bid_callback)
        , file_(mpi::File::open(filename, mpi::File::Mode::WRITE, MPI_COMM_SELF)) {
        auto n_buckets = bucket_counts.size();
        if (n_buckets == 0) {
            throw std::runtime_error("Need at least one bucket.");
        }

        bucket_file_offset[0] = 0ul;
        for (size_t i = 1; i < n_buckets; ++i) {
            size_t bucket_size = bucket_counts[i - 1] * sizeof(DataEntry);
            bucket_file_offset[i] = bucket_file_offset[i - 1] + bucket_size;
        }
    }

    ~BucketedWriter() {
        auto n_buckets = buckets.size();
        for (size_t bid = 0; bid < n_buckets; ++bid) {
            write_bucket(bid);
        }
    }

    /** \brief Write `data_entry` to its bucket.
     *
     *  The implied disk access might be delayed and the `data_entry` is kept in
     *  a buffer until enough data has accumulated to enable efficient file
     *  access.
     */
    void write(const DataEntry& data_entry) {
        auto bid = gid2bid(data_entry.preCell);
        auto& bucket = buckets[bid];
        bucket.push_back(data_entry);

        if (bucket.size() * sizeof(DataEntry) >= max_bucket_size) {
            flush_bucket(bid);
        }
    }

  protected:
    /** \brief Flush bucket with ID `bid` to disk.
     *
     *  Note, after flushing the bucket is ready to receive new entries.
     */
    void flush_bucket(size_t bid) {
        write_bucket(bid);
        clear_bucket(bid);
    }

    /** \brief Only write the bucket with ID to disk. */
    void write_bucket(size_t bid) {
        auto phase = Instrumentor::phase("write_bucket");
        const auto& bucket = buckets[bid];

        size_t count = bucket.size();
        if (count > 0) {
            file_.write_array_at(std::span{bucket.data(), count}, bucket_file_offset[bid]);
            bucket_file_offset[bid] += count * sizeof(DataEntry);
        }
    }

    /** \brief Ensure that the bucket is empty. */
    void clear_bucket(size_t bid) {
        auto& bucket = buckets[bid];
        if (bucket.capacity() * sizeof(DataEntry) >= 3 * max_bucket_size) {
            // A common pattern to force deallocation is swapping with an
            // empty vector and forcing the destructor to run on the swapped
            // into vector.
            auto tmp = std::vector<DataEntry>{};
            bucket.swap(tmp);
        }

        // Note, this causes everything above to be zeroed out. Therefore,
        // this operation is linear in the previous size of the bucket.
        bucket.resize(0ul);
    }

    size_t gid2bid(int gid) {
        return gid2bid_callback_(gid);
    }

    // Threshold at which a bucket should be flushed.
    constexpr static size_t max_bucket_size = 64ul << 20;

  private:
    std::vector<std::vector<DataEntry>> buckets;
    std::vector<size_t> bucket_file_offset;
    std::function<size_t(size_t)> gid2bid_callback_;

    mpi::File file_;
};

}  // namespace io

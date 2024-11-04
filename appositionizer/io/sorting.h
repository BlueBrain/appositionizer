#pragma once

#include <functional>
#include <string>
#include <vector>

#include "structures.h"

namespace io {

/** \brief Bucket sort the temporary output files.
 *
 *  Sorts, by pre- and post-GID, each of the temporary files in `input_dirname`
 *  into `output_dir`. This method will perform the sorting out-of-core if needed.
 *
 *  The out-of-core sorting is done via buckets. The parameter `n_buckets` is
 *  the number of buckets and `gid2bid` computes the bucket given a particular
 *  pre-GID.
 */
void bucket_sort(const std::string& output_dir,
                 const std::string& input_dir,
                 size_t n_buckets,
                 const std::function<size_t(size_t)>& gid2bid);

/** \brief Sort a file of appositions into buckets.
 *
 *  For each pre-GID the callback `gid2bid` defines the bucket
 *  into which all appositions for that pre-GID are put. The number of
 *  appositions must be provided in `bucket_counts`.
 */
void sort_into_buckets(const std::string& output_filename,
                       const std::string& input_filename,
                       const std::vector<size_t>& bucket_counts,
                       const std::function<size_t(size_t)>& gid2bid);

/** \brief Sort the buckets of `filename`.
 *
 * Assume the file of appositions has been pre-sorted into buckets; and buckets are
 * in ascending order. The size of each bucket is given by `bucket_counts`.
 *
 * This function will sort the file, by sorting each of the buckets. Note that
 * a) each of the buckets must fit into RAM; b) the file is sorted inplace.
 */
void sort_buckets(const std::string& filename, const std::vector<size_t> bucket_counts);

/** \brief Create the index for sorted appositions.
 *
 * Note that appositions in `data_entries` must be sorted in ascending order by
 * pre-GID.
 */
std::vector<IndexEntry> create_index_for_sorted(const std::vector<io::DataEntry>& data_entries);


}  // namespace io

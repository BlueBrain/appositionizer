#pragma once

namespace index_magic {

/** \brief Given `n_elements`, how many buckets does one need if each bucket can hold at most
 * `bucket_size` elements?
 *
 * Mathematically this is `ceil(n_elements / bucket_size)`.
 */
inline size_t div_up(size_t n_elements, size_t bucket_size) {
    return (n_elements + bucket_size - 1) / bucket_size;
}

/** \brief Keep buckets as balanced as possible.
 *
 *  IDs are partitioned sequentially into buckets. The size between any two
 *  buckets differs by at most one element. Large buckets are have lower bucket
 *  id than small buckets.
 */
namespace fair_share {

/** \brief Which bucket does element `k_element` belong into? */
inline size_t element2bucket(size_t n_elements, size_t n_buckets, size_t k_element) {
    size_t n_large_buckets = n_elements % n_buckets;
    size_t small_bucket_size = n_elements / n_buckets;
    size_t large_bucket_size = small_bucket_size + 1;

    if (k_element < n_large_buckets * large_bucket_size) {
        return k_element / large_bucket_size;
    }

    return n_large_buckets + (k_element - n_large_buckets * large_bucket_size) / small_bucket_size;
}

}  // namespace fair_share
}  // namespace index_magic

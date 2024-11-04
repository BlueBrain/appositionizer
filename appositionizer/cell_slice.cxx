#include "cell_slice.h"  // Must be on the first substantial line.

std::vector<int> SlicedCells::cellsOverlapping(const cell::Group& cells) const {
    std::vector<int> send_ids;
    std::vector<int> send_counts(comm_size_);

    for (int r = 0; r < comm_size_; ++r) {
        const auto& rank_box = slice(r);
        for (const auto& cell: cells) {
            [&] {
                for (const auto& section: cell) {
                    for (const auto& segment: section) {
                        if (segment.boundingBox().overlaps(rank_box)) {
                            send_ids.push_back(cell.id());
                            ++send_counts[r];
                            return;
                        }
                    }
                }
            }();
        }
    }

    std::vector<int> recv_counts(comm_size_);
    mpi::ensure(
        MPI_Alltoall, send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm_);

    std::vector<int> recv_offsets(comm_size_ + 1);
    std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_offsets.begin() + 1);

    std::vector<int> send_offsets(comm_size_ + 1);
    std::partial_sum(send_counts.begin(), send_counts.end(), send_offsets.begin() + 1);

    std::vector<int> recv_ids(recv_offsets.back());
    mpi::ensure(MPI_Alltoallv,
                send_ids.data(),
                send_counts.data(),
                send_offsets.data(),
                MPI_INT,
                recv_ids.data(),
                recv_counts.data(),
                recv_offsets.data(),
                MPI_INT,
                comm_);

    return recv_ids;
}

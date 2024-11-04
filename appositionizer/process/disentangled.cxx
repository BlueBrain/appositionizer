#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <mutex>
#include <thread>

#include <oneapi/tbb/parallel_sort.h>

#include "disentangled.h"

#include "instrumentor.h"
#include "data/apposition.h"
#include "process/apposition_filter.h"
#include "arange.h"
#include "better_mpi.h"
#include "index_magic.h"
#include "io/filereader.h"
#include "io/filewriter.h"
#include "io/sonatawriter.h"
#include "io/sorting.h"

using namespace oneapi;
using namespace std::literals::chrono_literals;

namespace fs = std::filesystem;

constexpr auto LOG_INTERVAL = 300s;

// namespace {

static const fs::path tmp_unsorted_dir = "tmp";
static const fs::path tmp_sorted_dir = "tmp.sorted";

std::vector<size_t> get_all_ids(const cell::Group& cell_group, MPI_Comm comm) {
    const size_t comm_size = mpi::size(comm);

    const int my_size = cell_group.size();
    std::vector<size_t> my_ids;
    my_ids.reserve(my_size);
    for (const auto& m: cell_group.metadata()) {
        my_ids.push_back(m.cellId);
    }

    std::vector<int> recv_counts(comm_size);
    mpi::ensure(MPI_Allgather, &my_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

    size_t total_count = std::accumulate(std::begin(recv_counts), std::end(recv_counts), 0);
    std::vector<size_t> all_ids(total_count);

    std::vector<int> recv_offsets(comm_size + 1);
    std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_offsets.begin() + 1);

    mpi::ensure(MPI_Allgatherv,
                my_ids.data(),
                my_ids.size(),
                MPI_UNSIGNED_LONG,
                all_ids.data(),
                recv_counts.data(),
                recv_offsets.data(),
                MPI_UNSIGNED_LONG,
                comm);

    return all_ids;
}


void merge_output_files(const cell::Group& cell_group,
                        const process::Configuration& cfg,
                        MPI_Comm comm) {
    const auto phase = Instrumentor::phase("file_merging");
    auto all_ids = get_all_ids(cell_group, comm);

    std::sort(all_ids.begin(), all_ids.end());
    all_ids.erase(std::unique(all_ids.begin(), all_ids.end()), all_ids.end());

    std::string unsorted_dir = cfg.tmp_dir / tmp_unsorted_dir;
    std::string sorted_dir = cfg.tmp_dir / tmp_sorted_dir;

    fs::create_directories(sorted_dir);
    size_t n_cells = all_ids.back() + 1;

    // pre-sort step:
    {
        size_t n_buckets = 128ul;
        auto gid2bid = [n_cells, n_buckets](size_t gid) {
            return index_magic::fair_share::element2bucket(n_cells, n_buckets, gid);
        };
        io::bucket_sort(sorted_dir, unsorted_dir, n_buckets, gid2bid);

        std::filesystem::remove(io::index_filename(unsorted_dir, MPI_COMM_SELF));
        std::filesystem::remove(io::payload_filename(unsorted_dir, MPI_COMM_SELF));
    }

    const size_t comm_size = mpi::size(comm);

    // Needs to match io::FileWriter in `disentangled`
    io::FileReader in(sorted_dir, io::generate_extension(MPI_COMM_SELF));
    uint64_t local_count = in.count_total();
    uint64_t global_count;

    // aggregating the counts

    mpi::ensure(MPI_Allreduce, &local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM, comm);

    std::unique_ptr<io::IOWriter> output_writer;
    if (cfg.output_population.empty()) {
        output_writer = std::make_unique<io::FileWriter>(cfg.output_dir, 0, 0, comm);
    } else {
        std::string file_name = cfg.output_dir / "edges.h5";
        output_writer = std::make_unique<io::SonataWriter>(file_name,
                                                           cfg.output_population,
                                                           cfg.source_population,
                                                           cfg.target_population,
                                                           global_count,
                                                           comm);
    }

    const size_t batchsize = 100;
    const size_t batches = (all_ids.size() + comm_size * batchsize - 1) / (comm_size * batchsize);

    auto exchange_type = mpi::datatype<io::DataEntry>();

    for (size_t batch = 0; batch < batches; ++batch) {
        std::vector<int> send_counts(comm_size);

        size_t batch_idx_begin = (batch * comm_size) * batchsize;
        size_t batch_idx_end = std::min(((batch + 1) * comm_size) * batchsize, all_ids.size());

        // This is just a flat loop from `batch_idx_begin` until `batch_idx_end`.
        for (size_t rank = 0; rank < comm_size; ++rank) {
            size_t idx_begin = batch_idx_begin + rank * batchsize;
            size_t idx_end = std::min(idx_begin + batchsize, all_ids.size());
            for (size_t idx = idx_begin; idx < idx_end; ++idx) {
                send_counts[rank] += in.count(all_ids.at(idx));
            }
        }

        Instrumentor::phase_begin("read data");

        auto send_entries = in.get_range(all_ids.at(batch_idx_begin),
                                         all_ids.at(batch_idx_end - 1) + 1);

        auto send_counts_total = std::accumulate(send_counts.begin(), send_counts.end(), 0ul);
        if (send_counts_total != send_entries.size()) {
            throw std::runtime_error(
                "The number of entries to be sent and the number of entries don't match.");
        }

        Instrumentor::phase_end("read data");
        Instrumentor::phase_begin("exchange data");

        std::vector<int> recv_counts(comm_size);
        mpi::ensure(
            MPI_Alltoall, send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

        const auto recv_counts_total = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);

        std::vector<io::DataEntry> recv_entries(recv_counts_total);

        std::vector<int> recv_offsets(comm_size + 1);
        std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_offsets.begin() + 1);

        std::vector<int> send_offsets(comm_size + 1);
        std::partial_sum(send_counts.begin(), send_counts.end(), send_offsets.begin() + 1);

        mpi::ensure(MPI_Alltoallv,
                    send_entries.data(),
                    send_counts.data(),
                    send_offsets.data(),
                    exchange_type,
                    recv_entries.data(),
                    recv_counts.data(),
                    recv_offsets.data(),
                    exchange_type,
                    comm);

        Instrumentor::phase_end("exchange data");
        Instrumentor::phase_begin("write data");

        Instrumentor::phase_begin("tbb::parallel_sort");
        tbb::parallel_sort(recv_entries.begin(), recv_entries.end());
        Instrumentor::phase_end("tbb::parallel_sort");

        auto recv_meta = create_index_for_sorted(recv_entries);

        output_writer->write(recv_entries, std::move(recv_meta));

        Instrumentor::phase_end("write data");
    }

    std::filesystem::remove(io::index_filename(sorted_dir, MPI_COMM_SELF));
    std::filesystem::remove(io::payload_filename(sorted_dir, MPI_COMM_SELF));
}


void process_cells(const cell::Group& cell_group,
                     const std::string& output,
                     size_t job_size,
                     const SlicedCells& sliced_cells,
                     const InterBoutonInterval& interval,
                     MPI_Comm comm) {
    auto global_phase = Instrumentor::phase("computational_core");

    std::atomic<size_t> job_id{0};
    std::atomic<size_t> processed_cells{0};

    std::mutex m_write;
    std::mutex m_logging;

    auto last_log = std::chrono::steady_clock::now() - LOG_INTERVAL;
    size_t last_processed = 0;

    const auto all_overlapping_cell_ids = sliced_cells.cellsOverlapping(cell_group);

    fs::create_directories(fs::path(output) / tmp_unsorted_dir);
    io::FileWriter writer(fs::path(output) / tmp_unsorted_dir, 0, 0, MPI_COMM_SELF);
    cell::Group all_cells(cell_group, all_overlapping_cell_ids);
    const auto& all_metadata = all_cells.metadata();
    auto all_metadata_ids = [&]() {
        std::vector<std::string> morphology_names;
        morphology_names.reserve(all_metadata.size());
        for (const auto& md: all_metadata) {
            morphology_names.push_back(md.cellName.data());
        }

        auto storage = morphio::Collection(all_cells.morphology_path());
        return storage.argsort(morphology_names);
    }();

    auto process_cells = [&](const tbb::blocked_range<size_t>& range) {
        auto region_phase_name = fmt::format("multi_threaded(id=cnu, tid={})",
                                             tbb::this_task_arena::current_thread_index());
        Instrumentor::phase_begin(region_phase_name.c_str());

        Instrumentor::phase_begin("reading");
        std::vector<cell::Metadata> to_use;
        for (size_t idx = range.begin(); idx != range.end(); ++idx) {
            to_use.push_back(all_metadata[all_metadata_ids[idx]]);
        }
        cell::Group cells(cell_group, to_use);
        cells.loadMorphologies([&](const cell::SectionType& type) {
            return cells.appositionSpace()->isPreSynaptic(type);
        });
        Instrumentor::phase_end("reading");

        Instrumentor::phase_begin("processing");
        std::map<size_t, std::map<size_t, std::vector<io::DataEntry>>> result;
        for (const auto& cell: cells) {
            for (const auto& section: cell) {
                Instrumentor::phase_begin("compute_overlap");
                std::vector<apposition::Apposition> appositions;
                for (const auto& s1: section) {
                    for (const auto& s2: sliced_cells.segmentsOverlapping(s1)) {
                        // check if segments can be compared
                        if (!(cells.appositionSpace()->canSegmentsBeCompared(s1.id(), s2.id()))) {
                            continue;
                        }

                        const auto [overlap, extends] = s1.overlaps(s2);
                        if (overlap) {
                            appositions.emplace_back(s1, s2, extends);
                        }
                    }
                }
                Instrumentor::phase_end("compute_overlap");

                Instrumentor::phase_begin("filter");
                std::sort(appositions.begin(), appositions.end());
                apposition::Filter filter(section, interval);
                auto regions = filter.filter(appositions);
                for (const auto& r: regions) {
                    for (const auto& c: r.connections) {
                        result[r.getCell1()][r.getCell2()].emplace_back(section[c.pre.segment],
                                                                            c.post.segment,
                                                                            section,
                                                                            c.pre.percentage,
                                                                            c.post.percentage);
                    }
                }
                Instrumentor::phase_end("filter");
            }
        }
        Instrumentor::phase_end("processing");

        {
            auto phase_writing = Instrumentor::phase("writing");
            std::lock_guard<std::mutex> lock(m_write);
            writer.write(result);
        }
        auto processed = processed_cells.fetch_add(to_use.size()) + to_use.size();

        {
            std::lock_guard<std::mutex> lock(m_logging);
            auto now = std::chrono::steady_clock::now();
            if (now - last_log > LOG_INTERVAL ||
                (processed - last_processed) / all_metadata.size() > 0.1 ||
                processed == all_metadata.size()) {
                mpi::info("{:5.1f} % :: Rank {:>5} processed job {:>7} with {:>7} cells",
                          100.0 * processed / all_metadata.size(),
                          mpi::rank(comm),
                          job_id++,
                          to_use.size());
                last_log = now;
                last_processed = processed;
            }
        }

        Instrumentor::phase_end(region_phase_name.c_str());
    };

    {
        auto phase = Instrumentor::phase("process_cells");
        tbb::parallel_for(tbb::blocked_range<size_t>(0, all_metadata.size(), 10 * job_size),
                          process_cells,
                          tbb::simple_partitioner());
    }

    MPI_Barrier(comm);
}


namespace process {

void disentangled(const cell::Group& cell_group,
                  const SlicedCells& sliced_cells,
                  const process::Configuration& cfg,
                  size_t job_size,
                  const InterBoutonInterval& interval,
                  MPI_Comm comm) {
    process_cells(cell_group, cfg.tmp_dir, job_size, sliced_cells, interval, comm);

    mpi::rank0::info("Finished processing, beginning file merging");
    merge_output_files(cell_group, cfg, comm);
    mpi::rank0::info("Finished file merging");
}

}  // namespace process

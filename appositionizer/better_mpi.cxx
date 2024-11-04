#include "better_mpi.h"  // Must be on the first substantial line.

namespace mpi {

File File::open(const std::string& filename, Mode mode, MPI_Comm comm) {
    MPI_File file;
    mpi::ensure(MPI_File_open, comm, filename.c_str(), mpi_mode(mode), MPI_INFO_NULL, &file);

    return File(file);
}

void File::free(MPI_File file) {
    mpi::ensure(MPI_File_close, &file);
}

MPI_File File::invalid_handle() noexcept {
    return MPI_FILE_NULL;
}

int File::mpi_mode(Mode mode) {
    if ((mode & Mode::READ) && (mode & Mode::WRITE)) {
        return MPI_MODE_CREATE | MPI_MODE_RDWR;
    } else if (mode & Mode::WRITE) {
        return MPI_MODE_CREATE | MPI_MODE_WRONLY;
    } else if (mode & Mode::READ) {
        return MPI_MODE_RDONLY;
    } else {
        throw std::runtime_error("Invalid mode.");
    }
}

}  // namespace mpi

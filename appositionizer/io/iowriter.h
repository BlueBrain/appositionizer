#pragma once

#include "structures.h"

namespace io {
class IOWriter {
  public:
    IOWriter() = default;
    virtual ~IOWriter() = default;
    virtual auto write(const std::vector<io::DataEntry>& appositions, std::vector<IndexEntry> meta)
        -> std::tuple<size_t, size_t, size_t, size_t> = 0;
};
}  // namespace io

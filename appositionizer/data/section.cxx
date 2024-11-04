#include "section.h"  // Must be on the first substantial line.

#include <type_traits>

namespace cell {

bool Section::overlaps(const BoundingBox& box) const {
    for (const auto& seg: *this) {
        if (seg.boundingBox().overlaps(box)) {
            return true;
        }
    }
    return false;
}

}  // namespace cell

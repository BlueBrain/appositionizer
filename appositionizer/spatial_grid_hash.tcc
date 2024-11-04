#include "spatial_grid_hash.h"  // must be on the first substantial line.

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <list>

#include "better_mpi.h"
#include "geometry.h"

using spatial::t_key;

constexpr t_key filter16bits = 0x000000000000ffff;
constexpr t_key filter32bits = 0x00000000ffffffff;

/**
 * Represents the abstract class for the types of data that can be indexed: ie a point or a segment;
 */

namespace spatial {
inline t_key toKey(long x, long y, long z) {
    // TODO: shall we do this 24 20 20 bits instead?
    t_key key = 0x0000000000000000;
    key |= static_cast<int16_t>(x) & filter16bits;
    key <<= 32;
    key |= static_cast<int32_t>(y) & filter32bits;
    key <<= 16;
    key |= static_cast<int16_t>(z) & filter16bits;
    return key;
}

template <typename T>
auto getKeys(T x1, T y1, T z1, T x2, T y2, T z2, T gridSpacing) -> std::vector<t_key> {
    std::vector<t_key> keys;

    long x_min = x1 / gridSpacing;
    long y_min = y1 / gridSpacing;
    long z_min = z1 / gridSpacing;

    long x_max = x2 / gridSpacing;
    long y_max = y2 / gridSpacing;
    long z_max = z2 / gridSpacing;

    keys.reserve((x_max - x_min + 1) * (y_max - y_min + 1) * (z_max - z_min + 1));

    // for every cell inside the volume mentioned before, test if segment crosses it
    for (auto x = x_min; x <= x_max; ++x) {
        for (auto y = y_min; y <= y_max; ++y) {
            for (auto z = z_min; z <= z_max; ++z) {
                keys.emplace_back(toKey(x, y, z));
            }
        }
    }

    // if this fails, then element was not indexed... bug in the code
    assert(keys.size() > 0);

    return keys;
}

template <typename T>
auto getKeys(Box<T> box, T gridSpacing) -> std::vector<t_key> {
    return getKeys(box.min.x, box.min.y, box.min.z, box.max.x, box.max.y, box.max.z, gridSpacing);
}
}  // namespace spatial

template <typename DataType, typename Coordinate>
SpatialGridHash<DataType, Coordinate>::SpatialGridHash(Coordinate gridSpacing)
    : _gridsSpacing(gridSpacing) {}

template <typename DataType, typename Coordinate>
void SpatialGridHash<DataType, Coordinate>::createIndex(std::vector<DataType>&& data) {
    assert(_spatialIndex.empty());  // Data may only be indexed once

    _data = std::move(data);

    for (const auto& d: _data) {
        if (!d.valid(_data))
            continue;

        auto box = d.boundingBox();
        testKeysRange(box.min, box.max);

        for (const auto& key: spatial::getKeys(box, _gridsSpacing)) {
            _spatialIndex[key].push_back(&d);
        }
    }

    // sorts inserted elements (in fact, sorts by memory position)
    for (auto& cell: _spatialIndex) {
        std::sort(cell.second.begin(), cell.second.end());
    }
}

template <typename DataType, typename Coordinate>
void SpatialGridHash<DataType, Coordinate>::testKeysRange(Coord3D<Coordinate>& coords1,
                                                          Coord3D<Coordinate>& coords2) {
    for (int i = 0; i < 3; i++) {
        int keyCode1 = floor(coords1[i] / _gridsSpacing);
        int keyCode2 = floor(coords2[i] / _gridsSpacing);

        if (keyCode1 >= std::numeric_limits<int16_t>::max() ||
            keyCode1 <= std::numeric_limits<int16_t>::min() ||
            keyCode2 >= std::numeric_limits<int16_t>::max() ||
            keyCode2 <= std::numeric_limits<int16_t>::min())
            mpi::warn(
                "The range of hash keys that your data requires, exceeds the capacity of the 16bit "
                "field in the spatial indexer key (dimension={})",
                i);
    }
}

template <typename DataType, typename Coordinate>
auto SpatialGridHash<DataType, Coordinate>::query(const Box<Coordinate>& box) const
    -> SpatialGridHash::QueryResult {
    QueryResult result;

    std::vector<t_key> keys = spatial::getKeys<Coordinate>(box, _gridsSpacing);
    assert(keys.size() > 0);

    // cells for the final results
    using CellIter = typename SpatialCell::const_iterator;
    using CellList = std::list<std::pair<CellIter, CellIter>>;
    CellList cells;
    // cells.reserve(keys.size());

    size_t max = 0;

    // get cells that belong to my indexer
    for (const auto& key: keys) {
        auto cand = _spatialIndex.find(key);
        // in case there's no data on that cell
        if (cand == _spatialIndex.end())
            continue;

        // append those cells to the results cells
        assert(cand->second.size() > 0);
        cells.emplace_back(cand->second.begin(), cand->second.end());
        max = std::max(max, cand->second.size());
    }

    result.reserve(max);

    using Iter = typename CellList::iterator;
    // Merge *sorted* lists of spatial index cell contents while also
    // dropping duplicated elements
    while (cells.size() > 1) {
        Iter cell = cells.begin();  // Need this to be able to increment++
        Iter min = cell++;
        // Loop over cells and find smallest element
        for (; cell != cells.end(); ++cell) {
            if (*cell->first < *min->first) {
                min = cell;
            } else if (*min->first == *cell->first) {
                // Value has already be seen, skip it next time around
                cell->first++;
                if (cell->first == cell->second) {
                    // Processed all values, don't process this cell any
                    // longer. The decrement will cause the loop to skip to
                    // the "right" next cell
                    cells.erase(cell--);
                }
            }
        }
        result.emplace_back(*(min->first++));
        if (min->first == min->second) {
            cells.erase(min);
        }
    }
    if (cells.size() == 1) {
        // May have removed two cell list entries in one while-loop
        // iteration
        result.insert(result.end(), cells.front().first, cells.front().second);
    }
    return result;
}

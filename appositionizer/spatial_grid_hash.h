#pragma once
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "geometry.h"

/**
 * Represents the implementation of a 3D spatial indexer for STATIC data,
 * using UNIFORM SPATIAL GRID as algorithm.
 * Ie the input volume is split into M*N*O sub-volumes of pre-defined sizes, and
 * data is indexed when constructor is called (therefore called static).
 * Once the data was inserted, the indexer can't be changed
 */
namespace spatial {
using t_key = uint64_t;
}

template <typename DataType, typename Coordinate>
class SpatialGridHash {
  public:
    using t_key = spatial::t_key;
    using SpatialCell = std::vector<const DataType*>;
    using QueryResult = std::vector<const DataType*>;

    /**
     * Constructor of class: constructs a hash grid with a single level
     */
    SpatialGridHash(Coordinate gridSpacing);

    // This class keep pointers internally on itself and is heavy.
    // It should not be copied.
    SpatialGridHash(const SpatialGridHash&) = delete;
    SpatialGridHash& operator=(const SpatialGridHash&) = delete;
    SpatialGridHash(SpatialGridHash&&) = default;
    SpatialGridHash& operator=(SpatialGridHash&&) = default;

    ~SpatialGridHash() = default;

    /**
     * Creates and populates the indexed map of elements given as input
     */
    void createIndex(std::vector<DataType>&& data);

    /**
     * Queries the spatial index for the bounding box provided
     * \param box contents from all cells overlapping this box are returned
     */
    auto query(const Box<Coordinate>& box) const -> QueryResult;

  private:
    std::vector<DataType> _data;
    /**
     * Variables that hold the spatial index itself:
     * ie, the map of elements, and the map saying how many elements exist
     * in each array indexed in the spatial map result
     */
    std::unordered_map<t_key, SpatialCell> _spatialIndex;
    Coordinate _gridsSpacing;

    /**
     * Checks if range of minimum and maximum range coordinates to be inserted
     * can fit into the data type defined as key for the spatial map.\n
     */
    void testKeysRange(Coord3D<Coordinate>& coords1, Coord3D<Coordinate>& coords2);
};

#include "spatial_grid_hash.tcc"

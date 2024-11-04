# Overview of Appositionizer {#overview}
The Appositionizer, in essence, takes a list of unconnected cells and computes
where to place appositions, i.e. places where synapses could form, between the
cells. This can, conceptually, be split into three tasks:

  1. Processing/loading the input data,

  2. Computing which parts of the cells are close enough to nearly touch. In
     Appositionizer we say they *overlap*.

  3. Deciding where to place appositions (synapses) within the
  detected overlapping regions.

Note that these these tasks do not happen sequentially. Rather, they are
interspersed as dictated by optimization and parallelization requirements. We
proceed to provide a high-level overview over each task.


## Input Data
One important thing to understand is that the decision which segments can be
pre- or post-synpatic is made before Appositionizer is run. The mechanisms of specifying
segments are pre- and which are post-synaptic is described in the User Guide in
the %CLI documentation.

In addition to the User Guide: please refer to `AppositionSpace`. Inside of Appositionizer
we do not write any code directly related to loading cells from disk. Instead, any
morphology loading is delegated to [MorphoKit]. Therefore, it can be reused in other part
of the circuit-building pipeline.

[MorphoKit]: https://bbpgitlab.epfl.ch/hpc/morpho-kit

## Computing Overlap {#overview_overlap}
Cells in the brain circuit consist of the soma, modelled as a sphere,
and pre- as well as post-synaptic sections.
Each section can be simplified as a sequence of cylinders, which are
stored compressed as a sequence of `SectionPoint` with a position and diameter
associated, where every two subsequent points form one cylinder or
`Segment`.

Naively finding all pairs of cylinders that have a non-empty intersection has
quadratic, in the number of cylinders, complexity. However in a pre-processing
step, one can place all cylinders into a tree, called spatial index. In doing
so, given any cylinder, one can quickly compute a small list of candidate
cylinders that might overlap with the cylinder. Thereby, the complexity becomes
log-linear in the number of cylinders. Please be aware that in practice the
void avoidance measures are more multi-layered.

The structure of this section is as follows. First we explain the distributed
spatial index; then we cover how the intersection of cylinder and spheres is
computed; and finally we explain how apposition detection is parallelized.


### Distributed Spatial Index {#dist_index}
Every rank of the Appositionizer loads a subset of post-synaptic cells via
`cell::Group`, which first reads the `cell::MetaData` from the circuit
file, and the corresponding spatial data from the morphology path via the
`cell::Cell` class.
The accumulated data is passed to the `cell_slicer`, which calculates the
centers of the bounding box for the first section (the soma) of each cell,
and divides the resulting coordinates up in three dimensions to determine the
distribution of the spatial index on the MPI ranks.

Coordinates are binned in the first dimension using the
`DistributedMemorySorter`, and the coordinates of each bin subsequently binned
in the next dimension, resulting in a tiled distribution with different areas
per bin, but similar number of items in each bin.
The extend of each bin is given by the minimum and maximum of each coordinate
there within.
This algorithm is called Sort Tile Recursion.

Sections are assigned to the rank associated with each final bin if the center
of their bounding box falls into the bin.
On each MPI rank, the individual `Segement`s are placed in a spatial index in
`SlicedCells`.


### Overlap Detection {#overlap_detection}
\note We'd like to draw you attention to the fact that two segments are
considered to overlap if they're within a certain distance of each other. This
is dictated by neuroscience; and therefore is documented in the User Guide.

The aim of overlap detection is two-fold: a) decide if the intersection of two
cylinders or a sphere with a cylinder is non-empty; b) to compute the
projection of the intersection onto the axes of the cylinders. While the
projection of a soma onto the axis of a segment can be non-empty, the converse,
i.e. the projection of the segment on the soma is assumed to be always empty.

Given two shapes that could plausibly overlap, there are two work-avoidance
measures:
1. If the bounding boxes of the two shapes don't overlap, no further work
   is needed.

2. Check if the shortest distance between the axis of the two cylinders is less
   than the sum of the radii of the two cylinders. This is an intersection test
   for two capsules.

At this point we're run out of work avoidance measures and a full overlap
detection must be performed. There are several methods for computing the
overlap of two cylinders. However only one is relevant:

* There's a sampling based algorithm, which is documented in
  `PartialCollision`. This is the default method (and neither
  CLI nor the build system allow configuring a different choice).


### Distributed Apposition Detection {#dist_detection}

#### Static Job Scheduling {#dist_detection_modern}
(This method is run when using the flag `--modern`.)

The static job scheduling approach is much simpler. As before each MPI rank is
responsible for computing overlap of a fixed number of post-synaptic sections,
with any pre-synaptic section.

Each MPI rank must read in all required sections and for each section check if
any of its cylinders overlap with any of the post-synaptic cylinders assigned
to this MPI rank. In order to parallelize the loop over the pre-synaptic
sections, each MPI rank is assigned multiple threads.

Detailed documentation can be found in the API documentation for
`process::disentangled`. Please note that this includes an important work
avoidance mechanism to reduce the number of pre-synaptic cells that need to
be considered on each MPI rank.


## Placing Synapses
The process of deciding where to place synapses given all pairs of intersecting
cylinders is of scientific importance and therefore documented in the User
Guide.

The relevant API documentation includes:
 - `apposition::Filter`
 - `apposition::Apposition`
 - `apposition::Region`

## Output {#outputs}

Following the apposition filtering and redistribution, output is written to both an index
file `edges.0` containing information about the global state and the apposition count and
offset per cell, as well as to a data file `edgesData.0` containing the connections
themselves, stored in contiguous fashion.

Upon successful processing of each batch of cells, the offsets of both
index and data files are stored in a `SaveState` file, which can be used
to resume processing after reconstruction of the distributed spatial
index in case of an unforeseen termination of the program. This is only
available when running [Classic Appositionizer](@ref dist_detection_old).

## See Also
The [Recipe documentation][recipe_docs] for details on the recipe input, and
the [Circuit Building documentation][circuit_building] for a holistic overview
of the brain circuit building.

[recipe_docs]: https://bbpteam.epfl.ch/documentation/projects/circuit-documentation/latest/recipe.html#consumers-and-invocation-order

[circuit_building]: https://bbpteam.epfl.ch/documentation/projects/circuit-build/latest/index.html

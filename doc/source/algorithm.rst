Apposition Detection Algorithm
=========================

A naive sequential version of Appositioner can be described by the following pseudo-code:

::

    cells = load_cells()
    pre_cells = filter(cells, is_presynaptic_cell)
    pre_sections = filter(pre_cells, is_presynaptic_section)

    post_cells = filter(cells, is_postsynaptic_cell)
    post_sections = filter(pre_cells, is_postsynaptic_section)

    overlaps = []
    for pre_segment in post_sections:
        for post_segment in pre_sections:
            # Compute the projection of the intersection of the pre- and
            # post-synaptic segment onto the axis of each segment. (`False` if
            # the two segments don't overlap.)
            overlap = compute_overlap(pre_segment, post_segment)

            if overlap:
                overlaps.append(overlap)

    # This is the neuroscientifically interesting part:
    connections = compute_connections(overlaps)


Scientifically, the most relevant parts of the algorithm are the `Apposition
Criteria`_ and how the purely geometrical information about the projection of
the intersection is converted into connections or more precisely appositions,
see `Filtering appositions`_.

However, let us briefly discuss a few minor points first:

1. Part of the input to Appositioner is choosing which cells can be pre- and which
   post-synaptic. Often all cells can be both pre- and post-synaptic. Then
   the sections of each cell are again filtered into sections that are pre-
   or post-synaptic. The conditions that collectively specify which which
   segments can be pre- and post-synaptic we call a *Apposition space*.  This
   enables computing the connection between axons of one region with dendrite
   of another region. In addition one must specify if segments from the same
   cell can intersect with each other. See :ref:`Command Line Interface` for
   further information on how one can specify each of the notions.

2. Appositioner isn't implemented as a double loop over all segments, since this would
   have complexity that's square in the number of segments. Which would yield
   prohibitive run times.  Instead a data structure, called a spatial index, is
   used to efficiently return all segments within a certain query box.
   Additionally, Appositioner makes use of domain decomposition. Since this additional
   complexity doesn't affect the output, it is of no neuroscientific
   importance. Nevertheless, it's an interesting computer science problem, you
   may find the details in the Developer Guide.

3. The overlap is not an exact projection of the exact intersection of two
   cylinders, instead it's an approximation, that works reliably for overlaps
   of a relevant size, but it may fail to detect cases where the two cylinders
   barely overlap, e.g. when a "corner" of a cylinder barely appositions the other
   cylinder. The precise method is documented in the Developer Guide.


Digital Representation of Cells
---------------------------------
Cells in the brain circuit consist of the soma, modelled as a sphere, and
pre- as well as post-synaptic sections. Each section is modelled as a sequence
of cylinders.
The overlap between two cells is considered to be directional, where a
pre-synaptic side of one cell may Apposition the post-synaptic side of another
cell.

Apposition Criteria
--------------
Segments are approximated as cylinders. A Apposition between two segments occurs
when the round surfaces of the two cylinders are less than one *spine length*
apart. Similarly, a soma, which is approximated as a sphere, appositions a segment
when the distance of their surfaces is less than a spine length. The spine
length is a user configurable parameter, see :ref:`Recipe <Recipe file>`.
Mathematically, these conditions can be translated into inflating the radius of
the cylinder or sphere appropriately. Which is why in Appositioner the notion of Appositioning
and intersecting are used synonymously.


Filtering appositions
-----------------
This section describes the process of placing connections given the purely
geometrical information of how the segments overlap. This corresponds to the
``compute_connections`` in the preceding pseudo-code.

At this point we've computed the projection of every (permissible) intersection
of two segments onto the axis of the pre- and post-synaptic segment. We refer to
these intervals as *appositions*. In a first step we join *appositions* into
*Apposition regions*. On the pre-synaptic side, appositions are merged into a Apposition
region if the two appositions lie on the same section and are less than a
user-specified threshold, called *region gap*, apart. The distance between two
appositions is measured along the center of the section connecting the two appositions.

The algorithm now advanced to the task of placing boutons on the pre-synaptic
section. A *connection* will refer to a pair of boutons, one on the pre- and
post-synaptic section each. A *bouton* refers to the point on the axis of a
segment where a connection begins or ends. On the pre-synaptic side for every
Apposition region, the boutons are distributed randomly, obeying the constraint that
the distance between two boutons must be within the inter-bouton interval. The
underlying random number generator is Random123_.

The next step in the algorithm is to find the bouton on the post-synaptic section
that forms the connection. In the easy case, the pre-synaptic bouton lies
inside a Apposition (not just a Apposition region). In this case the post-synaptic bouton
can be places in the corresponding post-synaptic side of the Apposition, at the same
relative position within the Apposition as on the pre-synaptic side. On the other
hand if the pre-synaptic bouton falls into a gap between two appositions, then the
post-synaptic bouton is either the start or end point of the post-synaptic side
of two appositions surrounding the pre-synaptic bouton.


.. _Random123: https://doi.org/10.1145/2063384.2063405

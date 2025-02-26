.. warning::
   The Blue Brain Project concluded in December 2024, so development has ceased under the BlueBrain GitHub organization.
   Future development will take place at: https://github.com/openbraininstitute/appositionizer

.. image:: doc/source/_static/banner.jpg
   :alt: A nice banner for appositionizer

appositionizer
=============

appositionizer looks at the physical location of cells morphologies,
and detects autaptic appositions between sections (see attachment for a
sample of 2 cells). An additional filtering step is added to
filter/remove appositions that lies within a given inter-bouton interval (in
practice a certain fixed distance between two autaptic appositions
throughout the axonal section).

In addition to this user guide we also provide a detailed
`Developer Guide <_static/doxygen/index.html>`_.

Installation
------------

A complete example of how to build appositionizer can be found in the
[Dockerfile](Dockerfile) in this repository.

For a fully manual build on Ubuntu, install the following essential libraries:

.. code::

   apt-get install build-essential catch2 cmake git libeigen3-dev libhdf5-openmpi-dev librandom123-dev librange-v3-dev libtbb-dev libyaml-cpp-dev ninja-build

Then make sure the following dependencies are installed and accessible via the environment
variable ``$CMAKE_PREFIX_PATH``:

* `libsonata <https://github.com/BlueBrain/libsonata>`_
* `MorphIO <https://github.com/BlueBrain/MorphIO>`_

Afterwards configure and build the project:

.. code::

   cmake -B build -S .
   cmake --build build
   cmake --install build

Notes for running
-----------------

Recipe components required
~~~~~~~~~~~~~~~~~~~~~~~~~~

appositionizer uses parts of a circuit building recipe stored either in JSON or YAML.
See the `recipe documentation
<https://sonata-extension.readthedocs.io/en/latest/recipe.html>`__ for further details.

The necessary components are the bouton interval specification, which controls when
cylinder-cylinder overlap region on the same section are merged (up to a distance of `region_gap` µm between adjacent overlap regions) and how appositions are then distributed along the merged region (spaced randomly using the specified minimum and maximum distances between appositions in µm).  An example configuration:

::

   "bouton_interval": {
     "min_distance": 5.0,
     "max_distance": 7.0,
     "region_gap": 5.0
   }

The second appositionizer specific setting is used to inflate the cylinders representing
apical or basal dendrites.
This allows to account for spines that would form between the actual cylindrical surfaces
of dendrites and axons.
Please note that spine lengths have be specified for all cell mtype values present in
the circuit as follows:

::

    "structural_spine_lengths": [
      {
        "mtype": "L1_SLAC",
        "spine_length": 2.5
      },
      {
        "mtype": "L23_PC",
        "spine_length": 2.5
      }
    ]

Acknowledgment
--------

The development of this software was supported by funding to the Blue Brain Project, a
research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss
government's ETH Board of the Swiss Federal Institutes of Technology.

Copyright (c) 2024 Blue Brain Project/EPFL

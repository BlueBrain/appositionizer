.. _`Using Appositionizer`:

Using Appositionizer
-------------------
The Appositionizer expects a list of unconnected but correctly positioned
cells as input; and (only) outputs where connections between the cell
could be placed. In particular, it doesn't output a fully functional circuit.
Therefore, it is one step of several steps in the circuit building pipeline.
As such Appositioner is not commonly used by itself. Please refer to the documentation of
the `circuit-building pipeline`_.

.. _`circuit-building pipeline`: https://bbpteam.epfl.ch/documentation/ecosystem/circuit-build/introduction.html

If you really do want to use Appositioner on its own, read on.

.. _`Command Line Interface`:

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~
The user interface of Appositioner is a command-line interface. Which takes the form

.. code-block:: shell

    appositionizer [--from POPULATION_NAME]
                  [--to   POPULATION_NAME]
                  [--from-nodeset NODESET_NAME]
                  [--to-nodeset   NODESET_NAME]
                  [--presynaptic SECTION_TYPES,...]
                  [--postsynaptic SECTION_TYPES,...]
                  [--Appositionspace Apposition_SPACE]
                  [--output OUTPUT_DIR]
                  [--modern]
                  --circuit-config CONFIG_FILE
                  --recipe RECIPE_FILE


Let's go through each of these in turn.

The flags ``--from`` and ``--to`` specify two population names, for the source population
``--from`` and the target population ``--to``. Synapses will be formed such that the
pre-synaptic cell is always a member of the source population and the post-synaptic
cell is always in the target population.
If no population name is passed, there must be exactly one population specified in
the circuit configuration file.

The next two flags ``--from-nodeset`` and ``--to-nodeset`` are also optional.
The flags enable setting the nodeset name. A nodeset is a means of selecting cells by
criteria such as the brain region or the M-type, contained in a file specified by the
circuit configuration. Please check the `SONATA Node Sets documentation`_ for details about the nodeset
file.

.. _`SONATA Node Sets documentation`: https://sonata-extension.readthedocs.io/en/latest/sonata_nodeset.html

The flags ``--presynaptic`` and ``--postsynaptic`` enable selecting which
sections of the pre- and post-synaptic cells can be pre- and post-synaptic,
respectively. Both accept a comma-separated list of section types. The
allowed values are ``all``, ``soma``, ``axon``, ``dendrite``, ``basal`` and
``apical``.

The flag ``--Appositionspace`` has two functions, on the one hand it can be set to
either:

- ``nonautaptic`` which will prevent forming synapses within the same cell.
- ``autaptic`` which will allow forming synapses between pre- and post-synaptic
  sections even if the pre- and post-synaptic segments belong to the same
  cell.

Since some combinations of arguments for ``--presynaptic``, ``--postsynaptic``
and ``--Appositionspace`` are exceedingly common, the following shorthand exists.
Setting ``--Appositionspace`` to

* ``axodendritic`` is equivalent to

  ::

    --presynaptic=axon --postsynaptic=soma,dendrite --Appositionspace=nonautaptic

* ``dendrosomatic`` is equivalent to

  ::

    --presynaptic=dendrite --postsynaptic=soma,dendrite --Appositionspace=nonautaptic

Then, there's ``--output`` which defines the output directory. If the output
directory is missing it is created.

The flag ``--modern`` selects the modern implementation of Appositioner. We
recommend using this. However, for reasons of least surprise, Appositioner continues to
use the original implementation if this flag isn't set.

There are two required arguments. One of them is the recipe filename, see
:ref:`Recipe file`. The second is the location of the circuit configuration that contains
all details about the circuit to be processed, as described by them
`SONATA Configuration documentation`_.

.. _`SONATA Configuration documentation`: https://sonata-extension.readthedocs.io/en/latest/sonata_config.html

Since Appositioner uses MorphoKit and MorphoIO to read morphologies defined in the circuit
configuration we refer to upstream documentation (`MorphoKit`_, `MorphIO`_) for further
details regarding morphology formats and storage.

.. _`MorphoKit`: https://bbpgitlab.epfl.ch/hpc/morpho-kit
.. _`MorphIO`: https://morphio.readthedocs.io/en/latest


.. _`Recipe file`:

Recipe
~~~~~~
The recipe is an XML file containing all parameters needed to generate the
connectome of a circuit. We will describe the parts of the recipe required by
Appositioner. Please refer to the upstream `documentation <recipe documentation>`_ for
additional information.

.. _`recipe documentation`: https://bbpteam.epfl.ch/documentation/projects/circuit-documentation/latest/recipe.html

Required Recipe Components
__________________________

Since version 4.4, Appositionizer will search the recipes for information about
the inter bouton interval and spine lengths. The former is used when
consolidating and redistributing appositions that are closely located to each other
on the same cell section. The XML component required follows the form:

::

    <InterBoutonInterval minDistance="5.0" maxDistance="7.0" regionGap="3.13"/>

Spine lengths, used to inflate the reach of cell sections, are found by
attributes only, with no distinct node name required:

::

    <SomeNode id="MY_MTYPE" spineLength="123.456"/>

All nodes having the attributes ``id`` and ``spineLength`` will be processed.
The ``id`` attribute has to contain an mtype identifier, and all spine length
specifications for the same mtype have to be identical.


Building Appositionizer
~~~~~~~~~~~~~~~~~~~~~~
As a user of Appositioner it is unlikely that you need or want to build Appositionizer
from source. Instead you can use the module on BB5. It can be loaded as follows

::

    module load unstable # or a particular archive.
    module load appositionizer

For details on compiling Appositioner from source, please refer to the `Build Instructions`_.

.. include:: developer_guide_links.rst

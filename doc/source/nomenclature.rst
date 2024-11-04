Nomenclature
------------
.. note::

   If you find any of these definitions unhelpful or incomplete, please let us
   know. We'd like this list to be useful.

Computer programs often borrow names from the specific application domain to
describe a digital notion of, hopefully, the same thing. Appositionizer is no
exception. This page tries to assign precise meaning to each term, as it's used
within Appositioner.

* **soma** is a sphere which approximates the shape of the biological soma.

* **segment** is essentially a cylinder. Jointly they form a section, which
  approximates part of a cell.

* **section** is a linear sequence of segments without any forking/sectioning.
  Often this concept is referred to as a section. The union of sections
  approximates the shape of axons or dendrites.

* **cell** refers to collection of sections and a soma, which jointly
  approximate the morphology of a cell. Additionally, cells are considered
  to have a concrete location in space. For the purposes of cells have no
  other properties, e.g. ion channels, etc. They are just a collection of
  simple geometric shapes.

* **Apposition** refers to the projection of the intersection of two
  segments onto their axis; or of a segment with a soma.

* **Apposition region** appositions on the same section are joined into a single Apposition
  region if they're separated by less than a certain threshold.

* **bouton** is the coordinates of the either the starting or end point of
  a connection on the axis of the segment.

* **connection** is a pair of two boutons one on the pre- and post-synaptic
  side each. Connections are where in later stages of the circuit building
  pipeline synapses might be placed. In neuroscience these are referred to
  as *appositions*.

* **synapses**, depending on the point of view, are not actually handled by Appositioner.
  Instead, Appositioner computes the location where synapses might occur (and leaves that
  decision to a later tool).  However, Appositioner does use the term as a synonym for
  connection.

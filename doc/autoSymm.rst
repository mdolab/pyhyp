.. _pyhyp_autoSymm:

Automatic symmetry BCs
======================

Using the :py:data:`unattachedEdgesAreSymmetry` option will automatically apply symmetry boundary conditions to any edges that do not interface with another block.
This page describes how this option works and a rare case where this option can fail.

Automatically applying symmetry boundary conditions involves determining the correct symmetry direction for each unattached edge.
For a standard case, the magnitudes of the coordinates for all nodes of an unattached edge are summed, and the coordinate direction with the minimum sum is set as the symmetry direction.
This works because we expect the geometry to lie flat on the symmetry plane.
However, when an unattached edge is coincident with one of the coordinate axes, this method cannot be used because there will be two directions with a zero sum.

In this case, a secondary check is performed.
The vectors between nodes on the boundary edge and their adjacent nodes on the first interior edge are computed, and the magnitudes of the coordinates are summed over all vectors.
The coordinate direction with the maximum sum is taken to be the symmetry direction.
This assumes that the aforementioned vectors point primarily in the direction normal to the symmetry plane.
This assumption does not always hold, such as for a highly swept wing.
If the vectors deviate enough from the normal direction, the symmetry direction will be incorrectly assigned, and the mesh extrusion will fail.

.. warning::
    
    If you encounter negative volumes near the symmetry plane, set the symmetry boundary conditions manually using :ref:`the BC option<pyhyp_BC>`.

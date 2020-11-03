.. _pyhyp_BC:

Using the BC option
===================

This page describes how to use the ``BC`` option to specify boundary conditions at boundary edges of the surface mesh.

Here is an example of a dictionary that can be used with ``BC``:

.. code-block::

   "BC":{1:{"iLow":"ySymm"}, 2:{"jHigh":"splay"}}

Each entry in the dictionary has a key, a nested key, and a value.
For the first entry, these are ``1``, ``iLow``, and ``ySymm``, respectively.
The key is the 1-based block number, the nested key is the boundary edge specification, and the value is the boundary condition.

The 1-based block number and boundary edge specification for a boundary edge can be determined using Tecplot:

#. Load the surface mesh file into Tecplot.
#. Use the Probe tool to select a point on the boundary edge of interest.
   Use Ctrl+click to snap on to a boundary edge.
#. Select Zone/Cell Info in the toolbar on the right.
#. The number shown after Zone is the 1-based block number.
#. The edge specification depends on the values of I and J.
   The edge is ``iLow`` if I=1, ``iHigh`` if I=I-Max, ``jLow`` if J=1, or ``jHigh`` if J=J-Max.
   Only one of these is true for any boundary edge.

The supported boundary conditions are:

* ``splay`` for free edges
* ``xSymm``, ``ySymm``, ``zSymm`` for symmetry planes
* ``xConst``, ``yConst``, ``zConst``, ``xyConst``, ``yzConst``, ``xzConst`` for constant planes

.. pyHyp documentation master file, created by
   sphinx-quickstart on Sun Oct 13 13:46:01 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _pyhyp:

=====
pyHyp
=====

pyHyp is a hyperbolic mesh generator that automatically generates two or three dimensional meshes around simple geometric configurations.
The basic idea is to start with an initial surface (or curve) corresponding to the geometry of interest and then *grow* or *extrude* the mesh in successive layers until it reaches a sufficient distance from the original surface.
In the process, the entire space surrounding the geometry is meshed.

.. _pyhyp_theory:

An overview of the hyperbolic mesh marching method implemented in pyHyp can be found in Section II.A of `Secco et al <https://arc.aiaa.org/doi/pdf/10.2514/1.J059491>`__.
Most of the theory for pyHyp was taken from `Chan and Steger <https://www.sciencedirect.com/science/article/pii/009630039290073A>`__.

Contents:

.. toctree::
   :maxdepth: 2

   install
   tutorial
   autoSymm
   icem
   BC
   options
   API
   adapt

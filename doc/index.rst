.. pyHyp documentation master file, created by
   sphinx-quickstart on Sun Oct 13 13:46:01 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _pyhyp:

=====
pyHyp
=====

Introduction
============

`pyHyp` is hyperbolic mesh generator that is capable of automatically
generating two or three dimensional 3D meshes around simple geometric
configurations. The basic idea is to start with an initial surface (or
line) corresponding to the geometry of interest and then *grow* or
*extrude* the mesh in successive layers until it reaches a sufficient
distance from the original surface. In the process, the entire space
surrounding the geometry is meshed.

.. _pyhyp_theory:

Most of the theory for `pyHyp` was taken from `Chan and Steger <https://www.sciencedirect.com/science/article/pii/009630039290073A>`_.



Contents:

.. toctree::
   :maxdepth: 2

   install
   tutorial
   icem
   BC
   options
   API

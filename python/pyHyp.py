#!/usr/bin/python
from __future__ import division
"""
pyHyp

The pyHyp module is used for generating hyperbolic grids given a set
of structured quad surface patches. 

Copyright (c) 2013 by G. Kenway
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 02/04/2013$

Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2010)

"""
# =============================================================================
# Standard Python modules
# =============================================================================
import sys, os, time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import import_modules, MPI, mpiPrint, MExt
exec(import_modules('geo_utils','pyGeo'))

# =============================================================================
# MultiBlockMesh class
# =============================================================================

class pyHyp(object):
    """
    This is the main class pyHyp. It is used as the user interface to pyHyp. 
    """

    def __init__(self, fileName, comm=None, Options=None, **kwargs):
        """
        Create the pyHyp object. 

        Input Arguments:
            fileName, str: 2D surface plot3D file name surface to use as 
                the base of the mesh
        Optional Arguments:
            comm, MPI_INTRA_COMM: MPI communication (as obtained 
                from mpi4py) on which to create the pyHyp object. This 
                allows faster computations in parallel. 
                If not provided, MPI_COMM_WORLD is used. 
            Options, dictionary: A dictionary containing the 
                 the options for hyperbolic mesh generator.

         Returns:
             None

             """

        # Set the possible MPI Intracomm
        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm
        # end if

        if Options is not None:
            self.options = Options
        else:
            self.options = {}
        # end if

        # Take the file name and try loading it using the pyGeo
        # architecture. Only do this on one proc; we will MPI-bcast
        # the information back to everyone else when it is necessary.
        if self.comm.rank == 0:
            geo_obj = pyGeo.pyGeo('plot3d', file_name=fileName,
                                  file_type='ascii', order='f')

            # This isn't prety or efficient but does produce the
            # globally reduced list we need and it is greedily
            # reordered. For ~100,000 surface points this takes on the
            # order of a couple of seconds so we're not going to worry
            # too much about efficiency here.
            nodeTol = kwargs.pop('nodeTol',1e-6)
            edgeTol = kwargs.pop('edgeTol',1e-6)
            geo_obj._calcConnectivity(nodeTol, edgeTol)
            sizes = []
            for isurf in xrange(geo_obj.nSurf):
                sizes.append([geo_obj.surfs[isurf].Nu, geo_obj.surfs[isurf].Nv])
            geo_obj.topo.calcGlobalNumbering(sizes)

            # Now we need to check that the surface is closed, In the
            # future it will support symmetry boundary conditions, but not yet

            # Next we need to produced an unstructured-like data
            # structure that points to the neighbours for each
            # node. We will use node-halo information such that
            # regular corners (has 2 neighbours in each of the i and j
            # directions) will be treated exactly the same as an
            # interior node. Extraordinary node, those that have 3, 5
            # or more verticies will be flagged and will be treated
            # specially in the hyperbolic marchign algorithm using a
            # lapalacian operator. 
            

        # Defalut options for mesh warping
        self.options_default = {
            # Number of layers:
            'N': 10, 

            # Offwall-spacing 
            'spacing':'constant',

            # Volume Smooth Factor
            'volSmooth':0.25,

            }

        self._checkOptions(self.options)

        return

      
    def _checkOptions(self, solver_options):
        """
        Check the solver options against the default ones
        and add opt
        ion iff it is NOT in solver_options
        """

        return solver_options

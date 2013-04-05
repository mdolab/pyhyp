#!/usr/bin/python
from __future__ import division
"""
pyHyp2D

The pyHyp2D module is used to generate a hyperbolic grid around a
closed curve

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
from mdo_import_helper import import_modules, MPI, mpiPrint
exec(import_modules('geo_utils','pyGeo'))

# =============================================================================
# MultiBlockMesh class
# =============================================================================

class pyHyp2D(object):
    """
    This is the main class pyHyp. It is used as the user interface to pyHyp2D. 
    """

    def __init__(self, fileName, comm=None, options=None, **kwargs):
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
            # Load the curve from a file using the the following format:
            # <number of points>
            # x1 y1
            # x2 y2 
            # ....
            # xn yn

            f = open(filename, 'r')
            N = int(f.readLine())
            x = []
            y = []
            for i in xrange(N):
                aux = f.readline().split()
                x.append(float(aux[0]))
                y.append(float(aux[1]))
            # end for

            # Done with the file so close
            f.close()

            # Convert the lists to an array
            x = numpy.array(x)
            y = numpy.array(y)

            # Now check if the first and the last point are within
            # tolerance of each other:
            if geo_utils.e_dist([x[0],y[0]],x[-1],y[-1]) < 1e-6:
                # Curve is closed, we're ok
                pass
            else:
                # Curve is not close, add first point at end and let
                # user know we are closing the curve
                x = numpy.append(x,x[0])
                y = numpy.append(y,y[0])
                print '*'*80
                print 'Warning: Curve was not closed! Closing curve with a linear\
 segment. This may or not be what is desired!'
                print '*'*80
            # end if
        else:
            N = None
            x = None
            y = None
        # end if

        # Now broadcase the required info to the other procs:
        self.N = self.comm.bcast(N, root=0)
        self.x = self.comm.bcast(x, root=0)
        self.y = self.comm.bcast(y, root=0)

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

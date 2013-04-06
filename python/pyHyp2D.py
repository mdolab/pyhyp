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
	v. 1.0 - Initial Class Creation (GKffffK, 2010)

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

    def __init__(self, fileName=None, X=None, comm=None, options=None, flip=False, **kwargs):
        """
        Create the pyHyp object. 

        Input Arguments: ONE of:
            fileName, str: file containing xy points 
            X, numpy array, size (N, 2): Numpy array containing N x-y points

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

        # Use supplied options
        if options is not None:
            self.options = options
        else:
            self.options = {}
        # end if

        # Import and set the hyp module
        import hyp
        self.hyp = hyp

        # Check to see how the user has passed in the data:
        if fileName is None and X is None:
            mpiPrint('Error: Either fileName or X must be passed on initialization!')
            return

        if fileName is not None and X is not None:
            mpiPrint('Error: BOTH fileName or X have been passed on initialization!')
            mpiPrint('Error: Only ONE of these are required.')
            return            

        # Take the file name and try loading it using the pyGeo
        # architecture. Only do this on one proc; we will MPI-bcast
        # the information back to everyone else when it is necessary.
        if self.comm.rank == 0:

            if fileName is not None:
                # Load the curve from a file using the the following format:
                # <number of points>
                # x1 y1
                # x2 y2 
                # ....
                # xn yn

                f = open(fileName, 'r')
                N = int(f.readline())
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
                X = numpy.array([x,y]).T
                
            # end if

            # Now check if the first and the last point are within
            # tolerance of each other:
            if geo_utils.e_dist(X[0],X[-1]) < 1e-6:
                # Curve is closed, we're ok. Internally, we work
                # with the unclosed curve and explictly deal with
                # the  0th and Nth points 
                X = X[0:-1,:]
            else:
                # Curve is not closed, print warning
                print '*'*80
                print 'Warning: Curve was not closed! Closing curve with a linear\
 segment. This may or not be what is desired!'
                print '*'*80
            # end if

            if flip: # Reverse direction
                X[:,0] = X[:,0][::-1]
                X[:,1] = X[:,1][::-1]
        # end if

        # Now broadcase the required info to the other procs:
        self.X = self.comm.bcast(X, root=0)
        self.N = len(self.X)



        # Defalut options for mesh warping
        self.options_default = {
            # Number of layers:
            'N': 10, 

            # Initial off-wall spacing
            's0':0.01,
            
            # Grid Spacing Ratio
            'gridRatio':1.15,
            # epsE: The explict smoothing coefficient
            'epsE': 0.5,

            # epsI: The implicit smoothing coefficient
            'epsI': 1.0,

            # theta: The barth implicit smoothing coefficient
            'theta': 1.0,

            # volCoef: The volume smoothing coefficinet
            'volCoef': 0.16,

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 10,
            }

        # Setup the options
        self._checkOptions()
        self._setOptions()
        self.gridGenerated = False
        mpiPrint('Problem sucessfully setup!')

        return

    def run(self):
        """
        Run given using the options given
        """

        self.hyp.run2d(self.X.T)
        self.gridGenerated = True

        return

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""
        
        if self.gridGenerated:
            self.hyp.writeplot3d_2d(fileName)
        else:
            mpiPrint('Error! No grid has been generated! Run the run() \
command before trying to write the grid!')
        # end if

        return

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly \
formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if self.gridGenerated:
            hyp.writecgns_2d(fileName)
        else:
            mpiPrint('Error! No grid has been generated! Run the run() \
command before trying to write the grid!')
        # end if

        return
      
    def _setOptions(self):
        """Internal function to set the options in pyHyp"""
        self.hyp.hypinput.n         = self.options['N']
        self.hyp.hypinput.s0        = self.options['s0']
        self.hyp.hypinput.gridratio = self.options['gridRatio']
        self.hyp.hypinput.epse      = self.options['epsE']
        self.hyp.hypinput.epsi      = self.options['epsI']
        self.hyp.hypinput.theta     = self.options['theta']
        self.hyp.hypinput.volcoef   = self.options['volCoef']
        self.hyp.hypinput.volsmoothiter = self.options['volSmoothIter']

        return

    def _checkOptions(self):
        """
        Check the solver options against the default ones
        and add opt
        ion iff it is NOT in solver_options
        """
        for key in self.options_default.keys():
            if not(key in self.options.keys()):
                self.options[key] = self.options_default[key]
            # end if
        # end for

        return

    def __del__(self):
        """
        Clean up fortran allocated values if necessary
        """
        self.hyp.releasememory()

        return

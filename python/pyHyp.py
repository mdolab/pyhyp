from __future__ import print_function
from __future__ import division
"""
pyHyp

The pyHyp module is used to generate a hyperbolic grid around a 3D
surface

Copyright (c) 2013 by G. Kenway
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 02/04/2013$

Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2013)
"""
# =============================================================================
# Imports
# =============================================================================
import sys
import os
import numpy
import tempfile
import time
from mpi4py import MPI
from . import MExt

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a explicitly raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyHypEllip Error: '
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)
        Exception.__init__(self)

# =============================================================================
# pyHyp class
# =============================================================================
class pyHyp(object):
    """
    This is the main class pyHyp. It is used as the user interface to pyHyp. 
    """

    def __init__(self, comm=None, options=None, debug=False, **kwargs):
        """
        Create the pyHyp object. 

        Parameters
        ----------
        comm : MPI_INTRACOMM
            Comm to use. This is used when running in parallel. If not
            provided, MPI.COMM_WORLD is used by default.

        options : dict
            A dictionary containing the the options for pyHyp.

        debug : bool
            Flag used to specify if debugging. This only needs to be
            set to true when using a symbolic debugger. 
            """
            
        # Default options for hyperbolic generation
        self.options_default = {
            # ---------------------------
            #        Input Information
            # ---------------------------
            # One of '3d' or '2d'
            'dimension':'3d',

            # Input file
            'inputFile':'default.fmt',

            # Type of extrusion: hyperbolic or elliptic
            'mode':'hyperbolic',

            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            # Number of layers:
            'N': 100, 

            # Initial off-wall spacing
            's0':0.01, 
            
            # Rmin: Distance to march in multiples of (automatically
            # computed) initial radius surrounding the body. 
            'rMin': 50, 

            # nodeTol: Tolerance for nodes to be treated as identical.
            'nodeTol':1e-8,

            # mirror: Possible mirroring of mesh. Use 'x', 'y', 'z'
            # if grid needs to be mirrored or None otherwise.
            'mirror':None,

            # panelEps: Distance source panels are "below" nodes. This
            # parameter usually doesn't need to be changed.
            'panelEps':1e-6,
            
            # ------------------------------------------
            #   Pseudo Grid Parameters (Hyperbolic only)
            # ------------------------------------------

            # Maximum permissible ratio of marching direction length
            # to smallest in-plane edge.
            'cMax':6.0, 

            #nonLinear: True/False. Use nonlinear scheme. Not
            #currently working. 
            'nonLinear': False, 

            #slExp: Exponent for Sl calc. Don't change this value
            #unless you know what you are doing! 
            'slExp': .15, 

            # Initial off-wall spacing.
            'ps0':0.01, 
            
            # Pseudo grid Spacing Ratio.
            'pGridRatio':1.15, 

            # ----------------------------------------
            #   Smoothing parameters (Hyperbolic only)
            # ----------------------------------------

            # epsE: The explicit smoothing coefficient
            'epsE': 1.0, 

            # epsI: The implicit smoothing coefficient
            'epsI': 2.0, 

            # theta: The barth implicit smoothing coefficient
            'theta': 3.0, 

            # volCoef: The volume smoothing coefficient for
            # pointJacobi iterations
            'volCoef': .16, 

            # volBlend: The volume blending coefficient to force
            # uniform sizes in farfield
            'volBlend': 0.0001, 

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 10, 

            # -------------------------------
            #   Solution Parameters (Common)
            # ------------------------------
            # kspRelTol: Solution tolerance for linear system
            'kspRelTol': 1e-8, 
            
            # Maximum number of iterations to run for linear system
            'kspMaxIts': 500, 

            # Preconditioner Lag (Hyperbolic only)
            'preConLag': 10, 

            'kspSubspaceSize':50, 

            # ---------------------------
            #   Output Parameters
            # ---------------------------
            # Debugging option to write grid metrics. Hyperbolic only
            'writeMetrics': False, 
            }
   
        # Import and set the hyp module
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.hyp = MExt.MExt('hyp', [curDir], debug=debug)._module
   
        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Initialize PETSc and MPI if not already done so:
        self.hyp.initpetsc(self.comm.py2f())

        # Use supplied options
        if options is None:
            raise Error("The options= keyword argument is *NOT* optional. "
                        "It must always be provided")
        self.options = options

        # Setup the options
        self._checkOptions()
        self._setOptions()
        self.gridGenerated = False

        # Initialize the problem based on dimension
        if self.options['dimension'] == '2d':
            self._init2d()
        else:
            self._init3d2()
            
        if self.comm.rank == 0:
            print('Problem successfully setup!')

    def _init3d2(self, **kwargs):
        """Main initialization routine"""
        if 'xMirror' in kwargs or 'yMirror' in kwargs or 'zMirror' in kwargs:
            raise Error("Mirror specification must be now done in the "
                        "options dictionary using: 'mirror':'z' for "
                        "example.")
        mirror = self.hyp.hypdata.nomirror
        if self.options['mirror'] == 'x':
            mirror = self.hyp.hypdata.xmirror
        elif self.options['mirror'] == 'y':
            mirror = self.hyp.hypdata.ymirror
        elif self.options['mirror'] == 'z':
            mirror = self.hyp.hypdata.zmirror
            
        self.hyp.setup3d(self.options['inputFile'], mirror)

    def run(self):
        """
        Run given using the options given
        """
        if self.options['dimension'] == '2d':
            self.hyp.run2d(self.X.T)
        else:
            self.hyp.run3d()
        self.gridGenerated = True

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""
        
        if self.gridGenerated:
            if self.options['dimension'] == '2d':
                self.hyp.writeplot3d_2d(fileName)
            else:
                self.hyp.writeplot3d_3d(fileName)
        else:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly 
        formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if self.gridGenerated:
            if self.options['dimension'] == '2d':
                self.hyp.writecgns_2d(fileName)
            else:
                self.hyp.writecgns_3d(fileName)
                # Possibly zero mirror plane for mirrored geometries
                # if self.symNodes is not None:
                #     mirrorDims = {'x':1, 'y':2, 'z':3}
                #     mirrorDim = mirrorDims[self.options['mirror']]
                #     self.hyp.zeromirrorplane(fileName, self.symNodes, mirrorDim)
        else:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

    def _setOptions(self):
        """Internal function to set the options in pyHyp"""
        self.hyp.hypinput.n = self.options['N']
        self.hyp.hypinput.s0 = self.options['s0']
        self.hyp.hypinput.rmin = self.options['rMin']
        self.hyp.hypinput.ps0 = self.options['ps0']
        self.hyp.hypinput.pgridratio = self.options['pGridRatio']
        self.hyp.hypinput.slexp = self.options['slExp']
        self.hyp.hypinput.epse = self.options['epsE']
        self.hyp.hypinput.epsi = self.options['epsI']
        self.hyp.hypinput.theta = self.options['theta']
        self.hyp.hypinput.volcoef = self.options['volCoef']
        self.hyp.hypinput.volblend = self.options['volBlend']
        self.hyp.hypinput.cmax = self.options['cMax']
        self.hyp.hypinput.volsmoothiter = self.options['volSmoothIter']
        self.hyp.hypinput.kspreltol = self.options['kspRelTol']
        self.hyp.hypinput.kspmaxits = self.options['kspMaxIts']
        self.hyp.hypinput.preconlag = self.options['preConLag']
        self.hyp.hypinput.nonlinear = self.options['nonLinear']
        self.hyp.hypinput.kspsubspacesize = self.options['kspSubspaceSize']
        self.hyp.hypinput.writemetrics = self.options['writeMetrics']
        self.hyp.hypinput.eps = 1e-8
        self.hyp.hypinput.farfieldtol = 10.0
        self.hyp.hypinput.nodetol = self.options['nodeTol']
    def _checkOptions(self):
        """
        Check the solver options against the default ones
        and add option iff it is NOT in solver_options
        """
        for key in self.options_default.keys():
            if not key in self.options.keys():
                self.options[key] = self.options_default[key]

    def __del__(self):
        """
        Clean up fortran allocated values if necessary
        """
        self.hyp.releasememory()

    def writeLayer(self, fileName, layer=1, meshType='plot3d', partitions=True):
        """
        Write a single mesh layer out to a file for visualization or for other purposes. 
        
        Parameters
        ----------
        fileName : str
            Filename to use. Should have .fmt extension for plot3d or .dat for tecplot
        layer : int
            Index of layer to print. Values greater than 1 are only valid if the mesh
            has already been extruded. 
        meshType : str
            Type of mesh to write. The two valid arguments are 'plot3d' and 'fe'. The plot3d
            will write the mesh in the original plot3d format while the FE mesh is the 
            unstructured internal representation. 
        partitions : bool
            This flag which is only used for the 'fe' mesh type option will write a separate
            zone for each partition on each processor. This is useful for visualizing the 
            parallel mesh decomposition. 
        """
        if meshType.lower() == 'plot3d':
            self.hyp.writelayerplot3d(fileName, layer)
        else:
            self.hyp.writelayerfe(fileName, layer, partitions)


    # def writeFEMesh(self, fileName):
    #     """ Ouput a tecplot FE mesh of the surface. Useful for
    #     debugging numberings.

    #     Parameters
    #     ----------
    #     fileName : str
    #         Filename of tecplot file. Should have .dat extension.
    #     """
    #     X = self.hyp.hypdata.xsurf.T
    #     conn = self.hyp.hypdata.conn.T
    #     nNodes = len(X)
    #     nElem = len(conn)
    #     f = open(fileName, 'w')
    #     f.write("FE Data\n")
    #     f.write('VARIABLES = \"X\", \"Y\", \"z\" \n')
    #     f.write("ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, "
    #             "ZONETYPE=FEQUADRILATERAL\n"%(nNodes, nElem))
    #     for i in xrange(nNodes):
    #         f.write('%f %f %f\n'%(X[i, 0], X[i, 1], X[i, 2]))
    #     for i in xrange(len(conn)):
    #         f.write('%d %d %d %d\n'%(conn[i][0], conn[i][1], 
    #                                  conn[i][2], conn[i][3]))
    #     f.close()
      

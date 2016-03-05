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
import os
import numpy
from mpi4py import MPI
from . import MExt

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a explicitly raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyHyp Error: '
        i = 14
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
        self.defaultOptions = {
            # ---------------------------
            #        Input Information
            # ---------------------------

            # Input surface file. 
            'inputFile':'default.fmt',

            # Input file type
            'fileType':'plot3d',

            # Type of extrusion: hyperbolic or elliptic
            'mode':'hyperbolic',

            # Unattached edges are symmetry plane
            'unattachedEdgesAreSymmetry':True,
    
            # Outerface boundary condition: farfield or overset
            'outerFaceBC':'farfield',

            # Boundary condition information for speicific block edges
            'BC':{},

            # Optional names for wall families
            'families':{},

            # AutoConnect: Run cgnsutilities connect function to add
            # any necessary block to block connectivity
            'autoConnect':True,

            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            # Number of layers:
            'N': 65, 

            # Initial off-wall spacing
            's0':0.01, 
            
            # Number of constant off-wall layers before beginning
            # stretch
            'nConstant':1,

            # sMax': Distance to march. 
            'marchDist':50,

            # nodeTol: Tolerance for nodes to be treated as identical.
            'nodeTol':1e-8,

            # splay: Splay BC spreading factor
            'splay': 0.25,

            # splayEdgeOrthogonality. How hard to try to force
            # orthogonality. Range of 0 to 0.5
            'splayEdgeOrthogonality':0.1,

            # splayCornerOrthogonality
            'splayCornerOrthogonality':0.2,

            # ---------------------------
            #   Elliptic Parameters
            # ---------------------------
            # panelEps: Distance source panels are "below" nodes. This
            # parameter usually doesn't need to be changed.
            'panelEps':1e-8,

            # farFieldTolerance: The multiple of the panel length
            # cutoff to use the approximation formula
            'farFieldTolerance': 4.0,

            # useMatrixFree: Use matrix-free solution
            # technique. (Always matrix free when evalMode is 'fast'.)
            'useMatrixFree': True,

            # evalMode: Type of panel evaluation routine: One of
            # exact, slow or fast. 
            'evalMode': 'fast',

            # sourceStrengthFile: File to use to load/save the source
            # strengths on the surface. 
            'sourceStrengthFile':'panelStrength.source',

            # ------------------------------------------
            #   Pseudo Grid Parameters (Hyperbolic only)
            # ------------------------------------------

            # Maximum permissible ratio of marching direction length
            # to smallest in-plane edge.
            'cMax':1.0, 

            #nonLinear: True/False. Use nonlinear scheme. Not
            #currently working. 
            'nonLinear': False, 

            #slExp: Exponent for Sl calc. Don't change this value
            #unless you know what you are doing! 
            'slExp': .15, 

            # Initial off-wall spacing. Negative values will be
            # calculated automatically.
            'ps0': -1,
            
            # Pseudo grid Spacing Ratio. Negative values will be
            # caculated automatically.
            'pGridRatio': -1,

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
            'volCoef': .25, 

            # volBlend: The volume blending coefficient to force
            # uniform sizes in farfield
            'volBlend': 0.0001, 

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 100, 
            
            # -------------------------------
            #   Solution Parameters (Common)
            # -------------------------------
            # kspRelTol: Solution tolerance for linear system
            'kspRelTol': 1e-8, 
            
            # Maximum number of iterations to run for linear system
            'kspMaxIts': 500, 

            # Subspace size for GMRES.
            'kspSubspaceSize':50, 

            # ---------------------------
            #   Output Parameters
            # ---------------------------
            # Debugging option to write grid metrics. Hyperbolic only
            'writeMetrics': False, 
            }
        self.defaultOptionKeys = set(k.lower() for k in self.defaultOptions)

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

        # Setup the options
        self.options = {}
        self._checkOptions(options)

        # Set the fortan options
        self._setOptions()
        self.gridGenerated = False

        # Convert file type to integer
        fileType = {'cgns':self.hyp.hypinput.cgnsfiletype, 
                    'plot3d':self.hyp.hypinput.plot3dfiletype}
        intFileType = fileType[self._go('fileType').lower()]

        # Check if input file exists.
        if not os.path.isfile(self._go('inputFile')):
            raise Error('Input file \'%s\' not found.'%self._go('inputFile'))

        # Determine the number of blocks we have so we can initialize
        # the BC array:
        nBlocks = self.hyp.getnblocks(self._go('inputFile'), intFileType)

        # The fortran BC information
        fBCs = numpy.zeros((4, nBlocks), order='f')
        fBCs[:, :] = self.hyp.hypinput.bcdefault

        # The python BC information
        BCs = self._go('BC')
        BCMap = {'splay':self.hyp.hypinput.bcsplay,
                 'xsymm':self.hyp.hypinput.bcxsymm,
                 'ysymm':self.hyp.hypinput.bcysymm,
                 'zsymm':self.hyp.hypinput.bczsymm,
                 'xconst':self.hyp.hypinput.bcxconst,
                 'yconst':self.hyp.hypinput.bcyconst,
                 'zconst':self.hyp.hypinput.bczconst,
                 'xyconst':self.hyp.hypinput.bcxyconst,
                 'yzconst':self.hyp.hypinput.bcyzconst,
                 'xzconst':self.hyp.hypinput.bcxzconst}

        edgeMap = {'ilow':self.hyp.hypinput.ilow-1,
                   'ihigh':self.hyp.hypinput.ihigh-1,
                   'jlow':self.hyp.hypinput.jlow-1,
                   'jhigh':self.hyp.hypinput.jhigh-1}
                 
        helpStr = "An example of a boundary specification is:  'BC':{1:{'iLow':'ySymm'}, 2:{'jHigh':'splay'}}"
        for blkBC in BCs:
            if blkBC < 1 or blkBC > nBlocks or not isinstance(blkBC, int):
                raise Error('Keys in BC array must be 1-based integers and less '
                            'than or equal to the total number of blocks. %s'%helpStr)
            for edgeKey in BCs[blkBC]:
                lKey = edgeKey.lower()
                if lKey not in edgeMap.keys():
                    raise Error("Boundary edge specification must be one of: "
                               "'iLow', 'iHigh', 'jLow', or 'jHigh'. %s"%helpStr)
                BCToSet = BCs[blkBC][edgeKey].lower()
                if BCToSet.lower() not in BCMap.keys():
                    raise Error("Boundary condition specification unknown. Must be one of: "
                                "'splay', 'BCXSymm', 'BCYSymm', 'BCZSymm', "
                                "'BCXConst', 'BCYConst', 'BCZConst, 'BCXYConst, "
                                "'BCYZConst or BCXZConst'. %s"%helpStr)
                    
                fBCs[edgeMap[lKey], blkBC-1] = BCMap[BCToSet]

        # Set the boundary condition information into fortran
        self.hyp.hypinput.bcs = fBCs

        # Now process the family information if we have any:
        families = self._go('families')

        fFamilies = []
        # Set default a default name of "wall". 
        for i in range(nBlocks):
            fFamilies.append("Wall")

        if isinstance(families, str):
            for i in range(nBlocks):
                fFamilies[i] = families

        elif  isinstance(families, dict):
            for blkBC in families:
                if blkBC < 1 or blkBC > nBlocks or not isinstance(blkBC, int):
                    raise Error('Keys in families dictionary must be 1-based integers and less '
                                'than or equal to the total number of blocks')
                fFamilies[blkBC-1] = families[blkBC]
        else:
            raise Error("'families' option must be a string or a dictionary. A string will "
                        "set all wall families to the single string. A dictionary with "
                        "one-based keys for the blocks may be used to specify individual "
                        "families for each block.\n Examples: 'families':'fuselage' or"
                        "'families':{'1:'fuselage', 2:'wing'}.")
                
        # Set our family information in fortran
        for i in range(nBlocks):
            self.hyp.setfamily(i+1, fFamilies[i])

        # Now run the fortran setup.
        self.hyp.setup(self._go('inputFile'), intFileType)

    def run(self):
        """
        Run given using the options given
        """
        if self._go('mode').lower() == 'hyperbolic':
            self.hyp.runhyperbolic()
        elif self._go('mode').lower() == 'elliptic':
            self.hyp.setuppanels()
            self.hyp.runelliptic()
        else:
            raise Error("Invalid parameter value for mode. Must be one "
                        "of 'elliptic' or 'hyperbolic'")
        self.gridGenerated = True

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""
        
        if self.gridGenerated:
            self.hyp.writeplot3d(fileName)
        else:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly 
        formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if not self.gridGenerated:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

        self.hyp.writecgns(fileName)

        # Possibly perform autoconnect using cgns_utils
        if self.comm.rank == 0 and self._go('autoConnect'):
            os.system('cgns_utils connect %s'%fileName)
        self.comm.barrier()

    def getSurfaceCoordinates(self):
        """Return the surface coordinates on this processor"""
        coords = numpy.zeros((self.hyp.hypdata.nx, 3))
        self.hyp.getsurfacecoordinates(numpy.ravel(coords))

        return coords

    def setSurfaceCoordinates(self, coords):
        """Set the surface coordinates on this processor"""
        self.hyp.setsurfacecoordinates(numpy.ravel(coords))

    def _setOptions(self):
        """Internal function to set the options in pyHyp"""
        self.hyp.hypinput.n = self._go('N')
        self.hyp.hypinput.nconstant = self._go('nConstant')
        self.hyp.hypinput.s0 = self._go('s0')
        self.hyp.hypinput.marchdist = self._go('marchdist')
        self.hyp.hypinput.ps0 = self._go('ps0')
        self.hyp.hypinput.pgridratio = self._go('pGridRatio')
        self.hyp.hypinput.slexp = self._go('slExp')
        self.hyp.hypinput.epse = self._go('epsE')
        self.hyp.hypinput.epsi = self._go('epsI')
        self.hyp.hypinput.theta = self._go('theta')
        self.hyp.hypinput.volcoef = self._go('volCoef')
        self.hyp.hypinput.volblend = self._go('volBlend')
        self.hyp.hypinput.cmax = self._go('cMax')
        self.hyp.hypinput.volsmoothiter = self._go('volSmoothIter')
        self.hyp.hypinput.splay = self._go('splay')
        self.hyp.hypinput.splayedgeorthogonality = self._go('splayEdgeOrthogonality')
        self.hyp.hypinput.splaycornerorthogonality = self._go('splayCornerOrthogonality')
        self.hyp.hypinput.kspreltol = self._go('kspRelTol')
        self.hyp.hypinput.kspmaxits = self._go('kspMaxIts')
        self.hyp.hypinput.nonlinear = self._go('nonLinear')
        self.hyp.hypinput.kspsubspacesize = self._go('kspSubspaceSize')
        self.hyp.hypinput.writemetrics = self._go('writeMetrics')
        self.hyp.hypinput.nodetol = self._go('nodeTol')
        self.hyp.hypinput.farfieldtol = self._go('farFieldTolerance')
        self.hyp.hypinput.usematrixfree = self._go('useMatrixFree')
        self.hyp.hypinput.unattachededgesaresymmetry = self._go('unattachEdedgesAreSymmetry')
        modes = {'exact':self.hyp.hypinput.eval_exact, 
                 'slow':self.hyp.hypinput.eval_slow,
                 'fast':self.hyp.hypinput.eval_fast}
        self.hyp.hypinput.evalmode = modes[self._go('evalMode').lower()]
        ffType = {'farfield':self.hyp.hypinput.outerfacefarfield,
                  'overset':self.hyp.hypinput.outerfaceoverset}
        self.hyp.hypinput.outerfacetype = ffType[self._go('outerFaceBC').lower()]
        self.hyp.hypinput.sourcestrengthfile[:] = ''
        f = self._go('sourceStrengthFile')
        self.hyp.hypinput.sourcestrengthfile[0:len(f)] = f

    def _checkOptions(self, options):
        """
        Check the solver options against the default ones
        and add option iff it is NOT in options
        """
        # Set existing ones
        for key in options:
            self.setOption(key, options[key])

        # Check for the missing ones
        optionKeys = set(k.lower() for k in options)
        for key in self.defaultOptions:
            if not key.lower() in optionKeys:
                self.setOption(key, self.defaultOptions[key])

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

    def freezeEdge(self, blockID, edge, dstar):
        """
        Specifiy an edge that will be frozen. 
        
        Parameters
        ----------
        blockID : integer
            Index of block IN ONE BASED ORDERING.
        edge : str
            String specified for edge. One of 'ilow', 'ihigh', 'jlow', 'jhigh'
        dstart : float
            How much these nodes will influence points around it. 
        """
        assert edge.lower() in ['ilow', 'ihigh', 'jlow', 'jhigh']
        self.hyp.freezeedge(blockID, edge, dstar)

    def freezeFaces(self, blockIDs, dstar):
        """
        Specifiy one or more faces (blocks) that will be frozen
        
        Parameters
        ----------
        blockIDs : integer or list
            Index of block(s) IN ONE BASED ORDERING.
        dstart : float
            How much these nodes will influence points around it. 
        """
        blockIDs = numpy.array(blockIDs).astype('intc')
        self.hyp.freezefaces(blockIDs, dstar)

    
    def surfaceSmooth(self, nIter, stepSize, surfFile=None):
        """
        Run smoothing iterations on the body surface

        Parameters
        ----------
        nIter : int
            Number of iterations to run
        stepSize : float
            Size of step. Must be < 1. Usually less than 0.1 for stability
            reasons. 
        """
        
        if surfFile is not None:
            try:
                from pygeo import pyGeo
            except:
                raise Error("pyGeo must be available to use the surface "
                            "reprojection object. Try again without specifying "
                            "the surfFile option.")
        
            geoSurf = pyGeo('iges', fileName=surfFile)
            self.hyp.allocatesurfaces(geoSurf.nSurf)
            for iSurf in range(geoSurf.nSurf):
                surf = geoSurf.surfs[iSurf]
                self.hyp.setsurface(iSurf+1, surf.ku, surf.kv, surf.tu, surf.tv, surf.coef.T)

        self.hyp.smoothwrap(nIter, stepSize)
    
    def getOption(self, name):
        """
        Return the value of the requested option. 

        Parameters
        ----------
        name : str
           Name of option to get. Not case sensitive

        Returns
        -------
        value : varries
           Return the curent value of the option.         
        """

        if name.lower() in self.defaultOptionKeys:
            return self.options[name.lower()]
        else:    
            raise Error('getOption: %s is not a volid pyHyp option.'%name)

    def setOption(self, name, value):
        """
        Set the value of the requested option. 

        Parameters
        ----------
        name : str
           Name of option to get. Not case sensitive

        value : varries
           Value to set
        """

        if name.lower() in self.defaultOptionKeys:
            self.options[name.lower()] = value
        else:    
            raise Error('setOption: %s is not a volid pyHyp option.'%name)
            
    def _go(self, name):
        """Internal short-cut function to make text a litle shorter"""
        return self.getOption(name)

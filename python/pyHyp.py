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
from . import hyp, MExt
from collections import OrderedDict
from copy import deepcopy

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
# pyHypMulti class
# =============================================================================
class pyHypMulti(object):

    """
    This is class can be used to run multiple pyHyp cases at once.
    """

    def __init__(self, comm=None, options=None, commonOptions=None, debug=False, skipList = []):

        '''
        The inititalization method will setup, run, and write all the results.
        INPUTS:
        options: ORDERED dictionary or list of dictionaries.
                 This contains options for the extrusion of all several grids. An example of
                 option dictionary is given below:

                 options = {'epsE':4.0,
                            'epsI':8.0,
                            'outputFile':'corner_hyp.cgns',
                            'skip':False}

                 We can set a list of dictionaries as input:

                 options1 = {'epsE':4.0,
                             'epsI':8.0,
                             'outputFile':'corner1_hyp.cgns',
                             'skip':False}
                 options2 = 'cartesian.cgns'
                 options3 = {'epsE':2.0,
                             'epsI':4.0,
                             'outputFile':'corner2_hyp.cgns',
                             'skip':False}
                 options = [options1, options2, options3]

                 OOOORRRR ... We can also set an ORDERED dictionary of dictionaries as input:

                 from collections import OrderedDict
                 options = OrderedDict()
                 options{'case1'} = {'epsE':4.0,
                                     'epsI':8.0,
                                     'outputFile':'corner1_hyp.cgns',
                                     'skip':False}
                 options{'block'} = 'cartesian.cgns'
                 options{'case2'} = {'epsE':2.0,
                                     'epsI':4.0,
                                     'outputFile':'corner2_hyp.cgns',
                                     'skip':False}

                 Each element of the list/dictionary will be considered as a different case.
                 One of the elements can be a string specifying a CGNS file that should be combined
                 with the other grids in the end. pyHyp will not do anything with this file except
                 combine it with the generated grids in the corresponding order.
                 These options will overwrite the default options (defined in the pyHyp class)
                 and the common options (another argument of this method).
                 If the user gives a list, this will be converted to a dictionary with integers
                 as keys. Remember this when setting the skip list for unnamed cases.

        commomOptions: dictionary with options that should be applied to all cases in the
                       options dictionary. See the 'defaulOptions' dictionary defined in
                       the pyHyp class to see the available options.

        skip_list: list containing names of cases that should be skipped.
        '''

        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Get processor ID
        myid = self.comm.Get_rank()

        # Convert input to dictionary even if user gave a single element
        if type(options) is dict:
            raise Error("pyHypMulti only accepts Ordered Dictionaries or Lists as inputs." + \
                        " Declare your options using options=OrderedDict()")
        elif type(options) is list:
            # Set ordered dict
            optionsDict = OrderedDict()
            # Convert list to dictionary using integers as keys
            optionsDict = {k:v for (k,v) in zip(range(len(options)),options[:])}
        else:
            # User gave an ordered dictionary
            optionsDict = deepcopy(options)

        # Add unused common options to each set
        for name in optionsDict:
            # Avoid options that are just strings indicating volume grids
            if type(optionsDict[name]) is not str:
                for key in commonOptions:
                    if key not in optionsDict[name].keys():
                        optionsDict[name][key] = commonOptions[key]

        # Initilize counter
        index = 0

        # Get the number of grids
        self.numGrids = len(optionsDict)

        # Initialize dictionary with results
        self.results = {'name':optionsDict.keys(),
                        'outputFile':[0]*self.numGrids,
                        'gridRatio':[0]*self.numGrids,
                        'minQualityOverall':[0]*self.numGrids,
                        'minVolumeOverall':[0]*self.numGrids,
                        }

        # Loop over all elements in the options list
        for optionName in optionsDict:

            options = optionsDict[optionName]

            if type(options) is str:

                # Just set up relationships to combine it with the other grids
                # later on
                self.results['name'][index] = optionName
                self.results['outputFile'][index] = options
                self.results['gridRatio'][index] = 'N/A'
                self.results['minQualityOverall'][index] = 'N/A'
                self.results['minVolumeOverall'][index] = 'N/A'

                # Increment counter
                index = index + 1

            elif optionName in skipList:
                if myid == 0:
                    print('Skipping case: ',optionName)

                # Get input file name
                try:
                    inputFile = options['inputfile']
                except:
                    inputFile = options['inputFile']

                # Check if output name exists or if we should get the
                # the automatically generated one
                try:
                    outputFile = options['outputfile']
                except:
                    try:
                        outputFile = options['outputFile']
                    except: # User probably did not set neither in options or common options
                        outputFile = None

                if outputFile is None:
                    # Get an automatic one
                    outputFile = generateOutputName(inputFile,outputType='cgns')

                # Save results
                self.results['name'][index] = optionName
                self.results['outputFile'][index] = outputFile
                self.results['gridRatio'][index] = 'skip'
                self.results['minQualityOverall'][index] = 'skip'
                self.results['minVolumeOverall'][index] = 'skip'


                # Increment counter
                index = index + 1

            elif type(options) is dict:

                # Only the root processor will print
                if myid == 0:
                    print('')
                    print('')
                    print('Running case %d : %s'%(index,options['inputFile']))
                    print('')

                # Create pyHyp object using the corresponding options
                hypGrid = pyHyp(comm, options, debug)

                # Run it
                hypGrid.run()

                # Write outputs
                hypGrid.writeOutput()

                # Save results
                self.results['name'][index] = optionName
                self.results['outputFile'][index] = hypGrid.options['outputfile']
                self.results['gridRatio'][index] = float(hypGrid.hyp.hypdata.gridratio)
                self.results['minQualityOverall'][index] = float(hypGrid.hyp.hypdata.minqualityoverall)
                self.results['minVolumeOverall'][index] = float(hypGrid.hyp.hypdata.minvolumeoverall)

                # Delete object to free memory
                del hypGrid

                # Increment counter
                index = index + 1

        # Print the log
        self.writeLog()

    def writeLog(self):
        """
        This will print a log with important information regarding all grids
        """

        # Get processor ID
        myid = self.comm.Get_rank()

        # Only the root processor will print
        if myid == 0:

            print('')
            print('')
            print('')
            print('='*40)
            print(' '*12 + 'pyHyp Summary' + ' '*12)
            print('='*40)

            print('')
            print('+----+-----------------------+-----------+-------------------+------------------+')
            print('| ID |       grid name       | gridRatio | minQualityOverall | minVolumeOverall |')
            print('+----+-----------------------+-----------+-------------------+------------------+')

            for index in range(self.numGrids):

                # Crop filename
                try:
                    filename = self.results['name'][index][:21]
                except:
                    filename = self.results['name'][index]

                # Get remaining data
                gridRatio = self.results['gridRatio'][index]
                minQualityOverall = self.results['minQualityOverall'][index]
                minVolumeOverall = self.results['minVolumeOverall'][index]

                # Format string that will be printed
                log_string1 = '| %02d '%index + '|' + \
                              ' {0: <21} '.format(filename) + '|'
                if type(gridRatio) is str:
                    log_string2 = ' {0: <9} '.format(gridRatio) + '|' + \
                                  ' '*4 + '{0: <10}'.format(minQualityOverall) + ' '*5 + '|' + \
                                  ' '*4 + '{0: <10}'.format(minVolumeOverall) + ' '*4 + '|'
                else:
                    log_string2 = ' {:9.7f} '.format(gridRatio) + '|' + \
                                  ' '*4 + '{:9.8f}'.format(minQualityOverall) + ' '*5 + '|' + \
                                  ' '*4 + '{:8.4e}'.format(minVolumeOverall) + ' '*4 + '|'

                print(log_string1+log_string2)

            print('+----+-----------------------+-----------+-------------------+------------------+')
            print('')
            print('')

    def combineCGNS(self,combinedFile='combined.cgns',additionalGrids=[''],skipList=[]):
        """
        This will gather all newly generated grids and combine them in
        a single cgns file. This only works for CGNS files
        """

        # Start cgns_utils command line
        command = 'cgns_utils combine '

        # Add names of generated grids
        selectedGrids = [self.results['outputFile'][ii] for ii in range(self.numGrids) \
                         if self.results['name'][ii] not in skipList]
        command = command + ' '.join(selectedGrids) + ' '

        # Add names of additional grids
        command = command + ' '.join(additionalGrids) + ' '

        # Add combined file
        command = command + combinedFile

        # Run cgns_utils only on first proc
        # Get processor ID
        myid = self.comm.Get_rank()
        # Only the root processor will print
        if myid == 0:
            print('COMBINE COMMAND:')
            print(command)
            print('')
            os.system(command)

            # Print log
            print('Combined CGNS files into: %s'%combinedFile)

# =============================================================================
# pyHyp class
# =============================================================================

class pyHyp(object):

    def __init__(self, comm=None, options=None, debug=False):
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

        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Default options for hyperbolic generation
        defaultOptions = {
            # ---------------------------
            #        Input Information
            # ---------------------------

            # Input surface file.
            'inputFile':'',

            # Explict numpy arrays of patches
            'patches':[],

            # Input file type
            'fileType':'plot3d',

            # Flag to entirely skip the grid generation of this block
            'skip':False,

            # Type of extrusion: hyperbolic or elliptic
            'mode':'hyperbolic',

            # Unattached edges are symmetry plane
            'unattachedEdgesAreSymmetry':True,

            # Outerface boundary condition: farfield or overset
            'outerFaceBC':'farfield',

            # Boundary condition information for specific block edges
            'BC':{},

            # Optional names for wall families
            'families':{},

            # AutoConnect: Run cgnsutilities connect function to add
            # any necessary block to block connectivity
            'autoConnect':True,

            # noPointReduce. Do not find duplicate nodes along
            # edges. This must ONLY be used with single surface input
            # files.
            'noPointReduce':False,

            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            # Number of layers:
            'N': 65,

            # Initial off-wall spacing
            's0':0.01,

            # Number of constant off-wall layers before beginning
            # stretch
            'nConstantStart':1,

            # Number of constant layers at the end of the march.
            'nConstantEnd':1,

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

            # Maximum convex corner angle necessary to trigger the node averaging.
            # It should be in degrees.
            'cornerAngle':60.0,

            # Coarsen: Automatically coarsen a mesh before starting
            # extrusion. Coarsen=1 gives the same mesh
            # (default). Coarsen=2 gives half the nodes in each direction.
            # Coarsen=3 gives *quarter* the nodes in each direction etc.
            'coarsen':1,

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

            # volSmoothSchedule: If given, use a user-supplied volume
            # smoothing schedule.
            'volSmoothSchedule':None,

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

            # Output format
            'outputType':'cgns',

            # Output filename (if None, an automatic one will be generated by
            # appending "_hyp" to the input filename).
            'outputFile':None,
            }

        # Get keys for every option
        self.defaultOptionKeys = set(k.lower() for k in defaultOptions)

        # Import and set the hyp module
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.hyp = MExt.MExt('hyp', [curDir], debug=debug)._module

        # Initialize PETSc and MPI if not already done so:
        self.hyp.initpetsc(self.comm.py2f())

        # Use supplied options
        if options is None:
            raise Error("The options = keyword argument is *NOT* optional. "
                        "It must always be provided")

        # Setup the options
        self.options = {}
        self._checkOptions(options, defaultOptions)

        # Set the fortan options
        self._setOptions()
        self.gridGenerated = False

        # Convert file type to integer
        fileType = {'cgns':self.hyp.hypinput.cgnsfiletype,
                    'plot3d':self.hyp.hypinput.plot3dfiletype}

        intFileType = fileType[self._go('fileType').lower()]

        # Determine how we are getting data: by Input file or
        # explictly by patches.
        patchInput = False
        patches = self._go('patches')

        if len(patches) > 0:
            patchInput = True
            nBlocks = len(patches)
        if not patchInput:
            if not os.path.isfile(self._go('inputFile')):
                raise Error('Input file \'%s\' not found.'%self._go('inputFile'))

            # Determine the number of blocks we have so we can initialize
            # the BC array:
            nBlocks = self.hyp.getnblocks(self._go('inputFile'), intFileType)

        self.hyp.allocatefamilies(nBlocks)
        if self.getOption('noPointReduce') and nBlocks > 1:
            raise Error('The noPointReduce option may only be true when '
                        'a single surface grid is provided.')

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

        # If we were given a CGNS file we might have families
        # there. So load them and overwrite the default.
        if intFileType == self.hyp.hypinput.cgnsfiletype:
            if self.comm.rank == 0:
                for i in range(nBlocks):
                    family, foundFam  = self.hyp.readfamily(self._go('inputFile'), i+1)
                    if foundFam and len(family.strip()) > 0:
                        fFamilies[i] = family.strip()
            fFamilies = self.comm.bcast(fFamilies)

        # If we have explictly other families given, these will
        # overwrite anything we already have.
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

        # Explictly set patches if necessary
        if patchInput:
            self.hyp.setnumberpatches(len(patches))
            for i in range(len(patches)):
                self.hyp.setpatch(i+1, patches[i])
            intFileType = self.hyp.hypinput.patchinput

        # Now run the fortran setup.
        self.hyp.setup(self._go('inputFile'), intFileType)

    def run(self):
        """
        Run given using the options given
        """
        if not self.options['skip']:
            if self._go('mode').lower() == 'hyperbolic':
                self.hyp.runhyperbolic()
            elif self._go('mode').lower() == 'elliptic':
                self.hyp.setuppanels()
                self.hyp.runelliptic()
            else:
                raise Error("Invalid parameter value for mode. Must be one "
                            "of 'elliptic' or 'hyperbolic'")
            self.gridGenerated = True
        else:
            print('Skipped generation of this grid')
            self.gridGenerated = False

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

    def writeOutput(self, fileName=None, fileType=None):
        """
        This selects the output type based on what is specified
        in the options
        """

        # Get desired output type from options, unless provided
        if fileType is None:
            outputType = self.options['outputtype']
        else:
            outputType = fileType

        # Check if the user specified name in the options
        if fileName is None:
            fileName = self.options['outputfile']

        # If no name is specified even in the options, then
        # we generate one
        if fileName is None:
            fileName = generateOutputName(self.options['inputfile'],outputType=outputType)

        # Update the name stored in the options for future uses
        self.options['outputfile'] = fileName

        if outputType == 'cgns':
            self.writeCGNS(fileName)
        elif outputType == 'plot3d':
            self.writePlot3d(fileName)

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
        self.hyp.hypinput.nconstantstart = self._go('nConstantStart')
        self.hyp.hypinput.nconstantend = self._go('nConstantEnd')
        self.hyp.hypinput.nopointreduce = self._go('noPointReduce')
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
        self.hyp.hypinput.cornerangle = self._go('cornerangle')*numpy.pi/180
        self.hyp.hypinput.coarsen = self._go('coarsen')
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
        self.hyp.hypinput.sourcestrengthfile = (
            self._expandString(self._go('sourceStrengthFile')))

        sch = self._go('volSmoothSchedule')
        if sch is not None:
            sch = numpy.array(sch, 'd')
            # Make sure its normalized
            low = sch[0, 0]
            high = sch[-1, 0]
            sch[:, 0] = (sch[:, 0]-low)/(high-low)
            self.hyp.hypinput.volsmoothschedule = sch

    def _checkOptions(self, options, defaultOptions):
        """
        Check the solver options against the default ones
        and add option iff it is NOT in options
        """
        # Set existing ones
        for key in options:
            self.setOption(key, options[key])

        # Check for the missing ones
        optionKeys = set(k.lower() for k in options)
        for key in defaultOptions:
            if not key.lower() in optionKeys:
                self.setOption(key, defaultOptions[key])

    def _expandString(self, s):
        """Expand a supplied string 's' to be of the constants.maxstring
        length so we can set them in fortran"""
        return s + ' '*(512-len(s))

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
            raise Error('getOption: %s is not a valid pyHyp option.'%name)

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
            raise Error('setOption: %s is not a valid pyHyp option.'%name)

    def _go(self, name):
        """Internal short-cut function to make text a litle shorter"""
        return self.getOption(name)


#=====================================================#

def generateOutputName(inputFile, outputType=None):
    """
    This function will automatically create an output filename by
    appending "_hyp" to the original filename.
    """

    # Get fileFormat from options if none is provided
    if outputType is None:
        outputType = self.options['outputtype']
    else:
        outputType = outputType.lower() # Use just lower case strings

    # Get rid of the extension and add an 'hyp' to the end of the filename
    outputFile = os.path.splitext(os.path.basename(inputFile))[0] + '_hyp'

    # Add extension acording to output file type
    if outputType == 'cgns':
        outputFile = outputFile + '.cgns'
    elif outputType == 'plot3d':
        outputFile = outputFile + '.fmt'
    else:
        raise Error("Output file type not recognized.")

    # Return the generated name
    return outputFile

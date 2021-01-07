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
from collections import OrderedDict
from copy import deepcopy
from baseclasses import BaseSolver
from baseclasses.utils import Error


# =============================================================================
# pyHypMulti class
# =============================================================================
class pyHypMulti(object):

    """
    This is class can be used to run multiple pyHyp cases at once.
    """

    def __init__(self, comm=None, options=None, commonOptions=None, debug=False, skipList=[]):

        """
        The inititalization method will setup, run, and write all the results.

        Parameters
        ----------
        options : object
            ORDERED dictionary or list of dictionaries.
            This contains options for the extrusion of all several grids. An example of
            option dictionary is given below:

            .. code-block:: python

               options = {'epsE':4.0,
                          'epsI':8.0,
                          'outputFile':'corner_hyp.cgns',
                          'skip':False}

            We can set a list of dictionaries as input:

            .. code-block:: python

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

            Alternatively, we can set an ORDERED dictionary of dictionaries as input:

            .. code-block:: python

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

        commomOptions : dict
            Dictionary with options that should be applied to all cases in the
            options dictionary. See the 'defOpts' dictionary defined in
            the pyHyp class to see the available options.

        skip_list : list
            List containing names of cases that should be skipped.
        """

        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Get processor ID
        myid = self.comm.Get_rank()

        # Convert input to dictionary even if user gave a single element
        if type(options) is dict:
            raise Error(
                "pyHypMulti only accepts Ordered Dictionaries or Lists as inputs."
                + " Declare your options using options=OrderedDict()"
            )
        elif type(options) is list:
            # Set ordered dict
            optionsDict = OrderedDict()
            # Convert list to dictionary using integers as keys
            optionsDict = {k: v for (k, v) in zip(range(len(options)), options[:])}
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
        self.results = {
            "name": list(optionsDict.keys()),
            "outputFile": [0] * self.numGrids,
            "gridRatio": [0] * self.numGrids,
            "minQualityOverall": [0] * self.numGrids,
            "minVolumeOverall": [0] * self.numGrids,
        }

        # Loop over all elements in the options list
        for optionName in optionsDict:

            options = optionsDict[optionName]

            if type(options) is str:

                # Just set up relationships to combine it with the other grids
                # later on
                self.results["name"][index] = optionName
                self.results["outputFile"][index] = options
                self.results["gridRatio"][index] = "N/A"
                self.results["minQualityOverall"][index] = "N/A"
                self.results["minVolumeOverall"][index] = "N/A"

                # Increment counter
                index = index + 1

            elif optionName in skipList:
                if myid == 0:
                    print("Skipping case: ", optionName)

                # Get input file name
                try:
                    inputFile = options["inputfile"]
                except KeyError:
                    inputFile = options["inputFile"]

                # Check if output name exists or if we should get the
                # the automatically generated one
                try:
                    outputFile = options["outputfile"]
                except KeyError:
                    try:
                        outputFile = options["outputFile"]
                    except KeyError:  # User probably did not set neither in options or common options
                        outputFile = generateOutputName(inputFile, outputType="cgns")

                # Save results
                self.results["name"][index] = optionName
                self.results["outputFile"][index] = outputFile
                self.results["gridRatio"][index] = "skip"
                self.results["minQualityOverall"][index] = "skip"
                self.results["minVolumeOverall"][index] = "skip"

                # Increment counter
                index = index + 1

            elif type(options) is dict:

                # Only the root processor will print
                if myid == 0:
                    print("")
                    print("")
                    print("Running case %d : %s" % (index, options["inputFile"]))
                    print("")

                # Create pyHyp object using the corresponding options
                hypGrid = pyHyp(comm, options, debug)

                # Run it
                hypGrid.run()

                # Write outputs
                hypGrid.writeOutput()

                # Save results
                self.results["name"][index] = optionName
                self.results["outputFile"][index] = hypGrid.options["outputfile"]
                self.results["gridRatio"][index] = float(hypGrid.hyp.hypdata.gridratio)
                self.results["minQualityOverall"][index] = float(hypGrid.hyp.hypdata.minqualityoverall)
                self.results["minVolumeOverall"][index] = float(hypGrid.hyp.hypdata.minvolumeoverall)

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

            print("")
            print("")
            print("")
            print("=" * 40)
            print(" " * 12 + "pyHyp Summary" + " " * 12)
            print("=" * 40)

            print("")
            print("+----+-----------------------+-----------+-------------------+------------------+")
            print("| ID |       grid name       | gridRatio | minQualityOverall | minVolumeOverall |")
            print("+----+-----------------------+-----------+-------------------+------------------+")

            for index in range(self.numGrids):

                # Crop filename
                try:
                    filename = self.results["name"][index][:21]
                except:  # noqa: E722
                    filename = self.results["name"][index]

                # Get remaining data
                gridRatio = self.results["gridRatio"][index]
                minQualityOverall = self.results["minQualityOverall"][index]
                minVolumeOverall = self.results["minVolumeOverall"][index]

                # Format string that will be printed
                log_string1 = "| %02d " % index + "|" + " {0: <21} ".format(filename) + "|"
                if type(gridRatio) is str:
                    log_string2 = (
                        " {0: <9} ".format(gridRatio)
                        + "|"
                        + " " * 4
                        + "{0: <10}".format(minQualityOverall)
                        + " " * 5
                        + "|"
                        + " " * 4
                        + "{0: <10}".format(minVolumeOverall)
                        + " " * 4
                        + "|"
                    )
                else:
                    log_string2 = (
                        " {:9.7f} ".format(gridRatio)
                        + "|"
                        + " " * 4
                        + "{:9.8f}".format(minQualityOverall)
                        + " " * 5
                        + "|"
                        + " " * 4
                        + "{:8.4e}".format(minVolumeOverall)
                        + " " * 4
                        + "|"
                    )

                print(log_string1 + log_string2)

            print("+----+-----------------------+-----------+-------------------+------------------+")
            print("")
            print("")

    def combineCGNS(self, combinedFile="combined.cgns", additionalGrids=[""], skipList=[], eraseFiles=True):
        """
        This will gather all newly generated grids and combine them in
        a single cgns file. This only works for CGNS output files.
        If eraseFiles=True, we erase the individual files that are combined.
        """

        # Start cgns_utils command line
        command = "cgns_utils combine "

        # Add names of generated grids
        selectedGrids = [
            self.results["outputFile"][ii] for ii in range(self.numGrids) if self.results["name"][ii] not in skipList
        ]
        command = command + " ".join(selectedGrids) + " "

        # Add names of additional grids
        command = command + " ".join(additionalGrids) + " "

        # Add combined file
        command = command + combinedFile

        # Run cgns_utils only on first proc
        # Get processor ID
        myid = self.comm.Get_rank()
        # Only the root processor will print
        if myid == 0:
            print("COMBINE COMMAND:")
            print(command)
            print("")
            os.system(command)

            # Print log
            print("Combined CGNS files into: %s" % combinedFile)

        # Erase input files
        if eraseFiles and (myid == 0):

            selectedGrids = [
                self.results["outputFile"][ii]
                for ii in range(self.numGrids)
                if self.results["name"][ii] not in skipList
            ]

            for filename in selectedGrids:

                os.remove(filename)


# =============================================================================
# pyHyp class
# =============================================================================


class pyHyp(BaseSolver):
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

        name = "pyHyp"
        category = "Hyperbolic mesh generator"

        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD

        # Default options for hyperbolic generation
        defOpts = self._getDefaultOptions()

        # Use supplied options
        if options is None:
            raise Error("The options = keyword argument is *NOT* optional. " "It must always be provided")

        # Initialize the inherited BaseSolver
        super().__init__(name, category, defOpts, options, comm=comm)

        # Import and set the hyp module
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.hyp = MExt.MExt("hyp", [curDir], debug=debug)._module

        # Initialize PETSc and MPI if not already done so:
        self.hyp.initpetsc(self.comm.py2f())

        # Set the fortan options
        self._setOptions()
        self.gridGenerated = False

        # Convert file type to integer
        fileType = {"cgns": self.hyp.hypinput.cgnsfiletype, "plot3d": self.hyp.hypinput.plot3dfiletype}

        intFileType = fileType[self.getOption("fileType")]

        # Determine how we are getting data: by Input file or
        # explictly by patches.
        patchInput = False
        patches = self.getOption("patches")

        if len(patches) > 0:
            patchInput = True
            nBlocks = len(patches)
        if not patchInput:
            if not os.path.isfile(self.getOption("inputFile")):
                raise Error("Input file '%s' not found." % self.getOption("inputFile"))

            # Determine the number of blocks we have so we can initialize
            # the BC array:
            nBlocks = self.hyp.getnblocks(self.getOption("inputFile"), intFileType)

        self.hyp.allocatefamilies(nBlocks)
        if self.getOption("noPointReduce") and nBlocks > 1:
            raise Error("The noPointReduce option may only be true when " "a single surface grid is provided.")

        # The fortran BC information
        fBCs = numpy.zeros((4, nBlocks), order="f")
        fBCs[:, :] = self.hyp.hypinput.bcdefault

        # The python BC information
        BCs = self.getOption("BC")
        BCMap = {
            "splay": self.hyp.hypinput.bcsplay,
            "xsymm": self.hyp.hypinput.bcxsymm,
            "ysymm": self.hyp.hypinput.bcysymm,
            "zsymm": self.hyp.hypinput.bczsymm,
            "xconst": self.hyp.hypinput.bcxconst,
            "yconst": self.hyp.hypinput.bcyconst,
            "zconst": self.hyp.hypinput.bczconst,
            "xyconst": self.hyp.hypinput.bcxyconst,
            "yzconst": self.hyp.hypinput.bcyzconst,
            "xzconst": self.hyp.hypinput.bcxzconst,
        }

        edgeMap = {
            "ilow": self.hyp.hypinput.ilow - 1,
            "ihigh": self.hyp.hypinput.ihigh - 1,
            "jlow": self.hyp.hypinput.jlow - 1,
            "jhigh": self.hyp.hypinput.jhigh - 1,
        }

        helpStr = "An example of a boundary specification is:  'BC':{1:{'iLow':'ySymm'}, 2:{'jHigh':'splay'}}"
        for blkBC in BCs:
            if blkBC < 1 or blkBC > nBlocks or not isinstance(blkBC, int):
                raise Error(
                    "Keys in BC array must be 1-based integers and less "
                    "than or equal to the total number of blocks. %s" % helpStr
                )
            for edgeKey in BCs[blkBC]:
                lKey = edgeKey.lower()
                if lKey not in edgeMap.keys():
                    raise Error(
                        "Boundary edge specification must be one of: "
                        "'iLow', 'iHigh', 'jLow', or 'jHigh'. %s" % helpStr
                    )
                BCToSet = BCs[blkBC][edgeKey].lower()
                if BCToSet.lower() not in BCMap.keys():

                    raise Error(
                        "Boundary condition specification unknown. Must be one of: "
                        "'splay', 'xSymm', 'ySymm', 'zSymm', "
                        "'xConst', 'yConst', 'zConst, 'xyConst, "
                        "'yzConst or xzConst'. %s" % helpStr
                    )

                fBCs[edgeMap[lKey], blkBC - 1] = BCMap[BCToSet]

        # Set the boundary condition information into fortran
        self.hyp.hypinput.bcs = fBCs

        # Now process the family information if we have any:
        families = self.getOption("families")

        fFamilies = []
        # Set default a default name of "wall".
        for i in range(nBlocks):
            fFamilies.append("Wall")

        # If we were given a CGNS file we might have families
        # there. So load them and overwrite the default.
        if intFileType == self.hyp.hypinput.cgnsfiletype:
            if self.comm.rank == 0:
                for i in range(nBlocks):
                    family, foundFam = self.hyp.readfamily(self.getOption("inputFile"), i + 1)
                    if foundFam and len(family.strip()) > 0:
                        fFamilies[i] = family.strip()
            fFamilies = self.comm.bcast(fFamilies)

        # If we have explictly other families given, these will
        # overwrite anything we already have.
        if isinstance(families, str):
            for i in range(nBlocks):
                fFamilies[i] = families

        elif isinstance(families, dict):
            for blkBC in families:
                if blkBC < 1 or blkBC > nBlocks or not isinstance(blkBC, int):
                    raise Error(
                        "Keys in families dictionary must be 1-based integers and less "
                        "than or equal to the total number of blocks"
                    )
                fFamilies[blkBC - 1] = families[blkBC]
        else:
            raise Error(
                "'families' option must be a string or a dictionary. A string will "
                "set all wall families to the single string. A dictionary with "
                "one-based keys for the blocks may be used to specify individual "
                "families for each block.\n Examples: 'families':'fuselage' or"
                "'families':{'1:'fuselage', 2:'wing'}."
            )

        # Set our family information in fortran
        for i in range(nBlocks):
            self.hyp.setfamily(i + 1, fFamilies[i])

        # Explictly set patches if necessary
        if patchInput:
            self.hyp.setnumberpatches(len(patches))
            for i in range(len(patches)):
                self.hyp.setpatch(i + 1, patches[i])
            intFileType = self.hyp.hypinput.patchinput

        # Now run the fortran setup.
        self.hyp.setup(self.getOption("inputFile"), intFileType)

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            # ---------------------------
            #        Input Information
            # ---------------------------
            "inputFile": [str, ""],
            "patches": [list, []],
            "fileType": [str, ["plot3d", "cgns"]],
            "skip": [bool, False],
            "mode": [str, ["hyperbolic", "elliptic"]],
            "unattachedEdgesAreSymmetry": [bool, True],
            "outerFaceBC": [str, ["farfield", "overset"]],
            "BC": [dict, {}],
            "families": [(str, dict), {}],
            "autoConnect": [bool, True],
            "noPointReduce": [bool, False],
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            "N": [int, 65],
            "s0": [float, 0.01],
            "nConstantStart": [int, 1],
            "nConstantEnd": [int, 1],
            "marchDist": [float, 50.0],
            "nodeTol": [float, 1e-8],
            "splay": [float, 0.25],
            "splayEdgeOrthogonality": [float, 0.1],
            "splayCornerOrthogonality": [float, 0.2],
            "cornerAngle": [float, 60.0],
            "coarsen": [int, 1],
            # ---------------------------
            #   Elliptic Parameters
            # ---------------------------
            "panelEps": [float, 1e-8],
            "farFieldTolerance": [float, 4.0],
            "useMatrixFree": [bool, True],
            "evalMode": [str, ["fast", "exact", "slow"]],
            "sourceStrengthFile": [str, "panelStrength.source"],
            # ------------------------------------------
            #   Pseudo Grid Parameters (Hyperbolic only)
            # ------------------------------------------
            "cMax": [float, 1.0],
            "nonLinear": [bool, False],
            "slExp": [float, 0.15],
            "ps0": [float, -1.0],
            "pGridRatio": [float, -1.0],
            # ----------------------------------------
            #   Smoothing parameters (Hyperbolic only)
            # ----------------------------------------
            "epsE": [float, 1.0],
            "epsI": [float, 2.0],
            "theta": [float, 3.0],
            "volCoef": [float, 0.25],
            "volBlend": [float, 0.0001],
            "volSmoothIter": [int, 100],
            "volSmoothSchedule": [(list, type(None)), None],
            # -------------------------------
            #   Solution Parameters (Common)
            # -------------------------------
            "kspRelTol": [float, 1e-8],
            "kspMaxIts": [int, 500],
            "kspSubspaceSize": [int, 50],
            # ---------------------------
            #   Output Parameters
            # ---------------------------
            "writeMetrics": [bool, False],
            "outputType": [str, ["cgns", "plot3d"]],
            "outputFile": [(str, type(None)), None],
        }
        return defOpts

    def run(self):
        """
        Run given using the options given
        """
        if not self.getOption("skip"):
            if self.getOption("mode") == "hyperbolic":
                self.hyp.runhyperbolic()
            elif self.getOption("mode") == "elliptic":
                self.hyp.setuppanels()
                self.hyp.runelliptic()
            self.gridGenerated = True
        else:
            print("Skipped generation of this grid")
            self.gridGenerated = False

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""

        if self.gridGenerated:
            self.hyp.writeplot3d(fileName)
        else:
            raise Error("No grid has been generated! Run the run() " "command before trying to write the grid!")

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly
        formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if not self.gridGenerated:
            raise Error("No grid has been generated! Run the run() " "command before trying to write the grid!")

        self.hyp.writecgns(fileName)

        # Possibly perform autoconnect using cgns_utils
        if self.comm.rank == 0 and self.getOption("autoConnect"):
            error = os.system("cgns_utils connect %s" % fileName)
            if error:
                raise Error(
                    "system command 'cgns_utils connect' failed, \
                            autoconnect was NOT performed"
                )
        self.comm.barrier()

    def writeOutput(self, fileName=None, fileType=None):
        """
        This selects the output type based on what is specified
        in the options
        """

        # Get desired output type from options, unless provided
        if fileType is None:
            outputType = self.getOption("outputtype")
        else:
            outputType = fileType

        # Check if the user specified name in the options
        if fileName is None:
            fileName = self.getOption("outputfile")

        # If no name is specified even in the options, then
        # we generate one
        if fileName is None:
            fileName = generateOutputName(self.getOption("inputfile"), outputType=outputType)

        # Update the name stored in the options for future uses
        self.setOption("outputfile", fileName)

        if outputType == "cgns":
            self.writeCGNS(fileName)
        elif outputType == "plot3d":
            self.writePlot3D(fileName)

    def getSurfaceCoordinates(self):
        """
        Return the surface coordinates on this processor
        """
        coords = numpy.zeros((self.hyp.hypdata.nx, 3))
        self.hyp.getsurfacecoordinates(numpy.ravel(coords))

        return coords

    def setSurfaceCoordinates(self, coords):
        """
        Set the surface coordinates on this processor
        """
        self.hyp.setsurfacecoordinates(numpy.ravel(coords))

    def _setOptions(self):
        """
        Internal function to set the options in pyHyp
        """
        self.hyp.hypinput.n = self.getOption("N")
        self.hyp.hypinput.nconstantstart = self.getOption("nConstantStart")
        self.hyp.hypinput.nconstantend = self.getOption("nConstantEnd")
        self.hyp.hypinput.nopointreduce = self.getOption("noPointReduce")
        self.hyp.hypinput.s0 = self.getOption("s0")
        self.hyp.hypinput.marchdist = self.getOption("marchdist")
        self.hyp.hypinput.ps0 = self.getOption("ps0")
        self.hyp.hypinput.pgridratio = self.getOption("pGridRatio")
        self.hyp.hypinput.slexp = self.getOption("slExp")
        self.hyp.hypinput.epse = self.getOption("epsE")
        self.hyp.hypinput.epsi = self.getOption("epsI")
        self.hyp.hypinput.theta = self.getOption("theta")
        self.hyp.hypinput.volcoef = self.getOption("volCoef")
        self.hyp.hypinput.volblend = self.getOption("volBlend")
        self.hyp.hypinput.cmax = self.getOption("cMax")
        self.hyp.hypinput.volsmoothiter = self.getOption("volSmoothIter")
        self.hyp.hypinput.splay = self.getOption("splay")
        self.hyp.hypinput.splayedgeorthogonality = self.getOption("splayEdgeOrthogonality")
        self.hyp.hypinput.splaycornerorthogonality = self.getOption("splayCornerOrthogonality")
        self.hyp.hypinput.cornerangle = self.getOption("cornerangle") * numpy.pi / 180
        self.hyp.hypinput.coarsen = self.getOption("coarsen")
        self.hyp.hypinput.kspreltol = self.getOption("kspRelTol")
        self.hyp.hypinput.kspmaxits = self.getOption("kspMaxIts")
        self.hyp.hypinput.nonlinear = self.getOption("nonLinear")
        self.hyp.hypinput.kspsubspacesize = self.getOption("kspSubspaceSize")
        self.hyp.hypinput.writemetrics = self.getOption("writeMetrics")
        self.hyp.hypinput.nodetol = self.getOption("nodeTol")
        self.hyp.hypinput.farfieldtol = self.getOption("farFieldTolerance")
        self.hyp.hypinput.usematrixfree = self.getOption("useMatrixFree")
        self.hyp.hypinput.unattachededgesaresymmetry = self.getOption("unattachEdedgesAreSymmetry")
        modes = {
            "exact": self.hyp.hypinput.eval_exact,
            "slow": self.hyp.hypinput.eval_slow,
            "fast": self.hyp.hypinput.eval_fast,
        }
        self.hyp.hypinput.evalmode = modes[self.getOption("evalMode")]
        ffType = {"farfield": self.hyp.hypinput.outerfacefarfield, "overset": self.hyp.hypinput.outerfaceoverset}
        self.hyp.hypinput.outerfacetype = ffType[self.getOption("outerFaceBC")]
        self.hyp.hypinput.sourcestrengthfile = self._expandString("")
        f = self.getOption("sourceStrengthFile")
        self.hyp.hypinput.sourcestrengthfile = self._expandString(f)

        sch = self.getOption("volSmoothSchedule")
        if sch is not None:
            sch = numpy.array(sch, "d")
            # Make sure its normalized
            low = sch[0, 0]
            high = sch[-1, 0]
            sch[:, 0] = (sch[:, 0] - low) / (high - low)
            self.hyp.hypinput.volsmoothschedule = sch

    def _expandString(self, s):
        """Expand a supplied string 's' to be of the constants.maxstring
        length so we can set them in fortran"""
        return s + " " * (512 - len(s))

    def __del__(self):
        """
        Clean up fortran allocated values if necessary
        """
        self.hyp.releasememory()

    def writeLayer(self, fileName, layer=1, meshType="plot3d", partitions=True):
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
        if meshType.lower() == "plot3d":
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
        assert edge.lower() in ["ilow", "ihigh", "jlow", "jhigh"]
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
        blockIDs = numpy.array(blockIDs).astype("intc")
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
            except ImportError:
                raise Error(
                    "pyGeo must be available to use the surface "
                    "reprojection object. Try again without specifying "
                    "the surfFile option."
                )

            geoSurf = pyGeo("iges", fileName=surfFile)
            self.hyp.allocatesurfaces(geoSurf.nSurf)
            for iSurf in range(geoSurf.nSurf):
                surf = geoSurf.surfs[iSurf]
                self.hyp.setsurface(iSurf + 1, surf.ku, surf.kv, surf.tu, surf.tv, surf.coef.T)

        self.hyp.smoothwrap(nIter, stepSize)


# =====================================================#


def generateOutputName(inputFile, outputType):
    """
    This function will automatically create an output filename by
    appending "_hyp" to the original filename.
    """

    outputType = outputType.lower()  # Use just lower case strings

    # Get rid of the extension and add an 'hyp' to the end of the filename
    outputFile = os.path.splitext(os.path.basename(inputFile))[0] + "_hyp"

    # Add extension acording to output file type
    if outputType == "cgns":
        outputFile = outputFile + ".cgns"
    elif outputType == "plot3d":
        outputFile = outputFile + ".fmt"

    # Return the generated name
    return outputFile

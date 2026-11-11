import os
import shutil
import sys
import tempfile
import unittest
from mpi4py import MPI
import numpy

baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "ref")

# Options from the sphere example (examples/sphere/runSphere.py), which can't be imported without running it
SPHERE_OPTIONS = {
    "inputFile": os.path.join(baseDir, "../examples/sphere/even_sphere.fmt"),
    "unattachedEdgesAreSymmetry": True,
    "outerFaceBC": "farfield",
    "autoConnect": False,
    "BC": {},
    "families": "wall",
    "N": 73,
    "s0": 1e-4,
    "marchDist": 10.0,
    "ps0": -1.0,
    "pGridRatio": -1.0,
    "cMax": 5.0,
    "epsE": 1.0,
    "epsI": 2.0,
    "theta": 3.0,
    "volCoef": 0.25,
    "volBlend": 0.0001,
    "volSmoothIter": 100,
}

# Relative tolerance for comparing runs on different numbers of procs. Results don't match to round-off, probably
# because the ASM preconditioner's subdomains depend on the partition and the linear solves stop at KSPRelTol.
PROC_COUNT_REL_TOL = 1e-6


class TestExamples(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        # Insert the repo root in path so that testflo can import examples
        sys.path.insert(0, os.path.join(baseDir, "../"))

    def commonTest(self, testFile, marchDist, relTol=1e-14):
        volumeName = os.path.split(testFile)[1]
        refFile = os.path.join(refDir, volumeName)

        # Set the absolute tolerance for float comparisons based on mesh dimensions
        absTol = marchDist * relTol

        self.assertCGNSFilesMatch(testFile, refFile, absTol)

    def assertCGNSFilesMatch(self, testFile, refFile, absTol):
        """Assert that cgnsdiff finds no data differences larger than absTol between two CGNS files."""
        # Run cgnsdiff and store the terminal output
        # HACK: For some reason subprocess.run doesn't work here on the Intel image
        # so we have to call os.system and then pipe the output to a temporary file
        # we use self.id() to make these unique so there are no race conditions
        # when tests are run in parallel
        TEMPFILE = f"TEMP_{self.id()}.txt"
        if MPI.COMM_WORLD.rank == 0:
            cmd = f"cgnsdiff -d -t {absTol} {testFile} {refFile}"
            exitCode = os.system(f"{cmd} > {TEMPFILE}")
            with open(TEMPFILE) as f:
                output = f.read()
            os.remove(TEMPFILE)
        else:
            exitCode = None
            output = None

        # Broadcast the output to all procs so we can assert on all procs
        exitCode = MPI.COMM_WORLD.bcast(exitCode, root=0)
        output = MPI.COMM_WORLD.bcast(output, root=0)

        # Check that cgnsdiff ran properly
        self.assertEqual(exitCode, 0)

        # Assert that there is no diff
        try:
            self.assertEqual(output, "")
        # Or that only the CGNS version differs
        except AssertionError:
            self.assertEqual(output, "/CGNSLibraryVersion <> /CGNSLibraryVersion : data values differ\n")

    def test2DEuler(self):
        from examples.naca0012.naca0012_euler import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def test2DRans(self):
        from examples.naca0012.naca0012_rans import extrudeDefaultCase

        hyp, volumeFile = extrudeDefaultCase()

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def test2DRansConstantLayers(self):
        from examples.naca0012.naca0012_rans import extrudeConstantLayersCase

        hyp, volumeFile = extrudeConstantLayersCase()

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def test2DRansSchedule(self):
        from examples.naca0012.naca0012_rans import extrudeScheduleCase

        hyp, volumeFile = extrudeScheduleCase()

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def test2DRansExplicit(self):
        from examples.naca0012.naca0012_rans import extrudeExplicitCase

        hyp, volumeFile = extrudeExplicitCase()

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def test717(self):
        from examples.wing_717.run717 import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def testBWB(self):
        from examples.BWB.runBWB import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def testCorner(self):
        from examples.corner.runCorner import volumeFile, commonOptions

        marchDist = commonOptions["marchDist"]
        # The relative tolerance is lowered because the absolute tolerance
        # depends on the march distance (absTol = relTol * marchDist). Since
        # the marchDist is very low here (2.5), the absolute tolerance is so
        # low that this test fails on some docker images.
        self.commonTest(volumeFile, marchDist, relTol=1e-9)

    def testM6(self):
        from examples.m6.runM6 import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def testPlate(self):
        from examples.plate.runPlate import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def testSphere(self):
        from examples.sphere.runSphere import volumeFile, commonOptions

        marchDist = commonOptions["marchDist"]
        self.commonTest(volumeFile, marchDist)

    def extrudeSphereInParallelAndSerial(self, surfaceOffset=None):
        """
        Extrude the even sphere example surface twice: partitioned across all procs, then on rank 0 alone.

        If surfaceOffset is given, the surface is translated by it using setSurfaceCoordinates before each extrusion.
        Returns the parallel and serial pyHyp objects and the CGNS file each one wrote.
        The serial pyHyp object is None on every rank except 0.
        """
        from pyhyp import pyHyp

        def extrude(hyp, outputFile):
            if surfaceOffset is not None:
                hyp.setSurfaceCoordinates(hyp.getSurfaceCoordinates() + surfaceOffset)
            hyp.run()
            hyp.writeCGNS(outputFile)

        comm = MPI.COMM_WORLD
        outputDir = tempfile.mkdtemp() if comm.rank == 0 else None
        outputDir = comm.bcast(outputDir, root=0)
        if comm.rank == 0:
            self.addCleanup(shutil.rmtree, outputDir)
        parallelFile = os.path.join(outputDir, "sphere_parallel.cgns")
        serialFile = os.path.join(outputDir, "sphere_serial.cgns")

        # Run the parallel case first so that PETSc is initialized collectively on all procs before rank 0
        # runs the serial case on its own
        parallelHyp = pyHyp(comm=comm, options=SPHERE_OPTIONS)
        extrude(parallelHyp, parallelFile)

        serialHyp = None
        if comm.rank == 0:
            serialHyp = pyHyp(comm=MPI.COMM_SELF, options=SPHERE_OPTIONS)
            extrude(serialHyp, serialFile)
        comm.barrier()

        return parallelHyp, serialHyp, parallelFile, serialFile

    def test3DMeshIndependentOfProcCount(self):
        """
        Given the even sphere example surface,
        when it is extruded once partitioned across all procs and once on a single proc,
        then the two volume meshes match to within the linear solver's tolerance.
        """
        _, _, parallelFile, serialFile = self.extrudeSphereInParallelAndSerial()

        self.assertCGNSFilesMatch(serialFile, parallelFile, SPHERE_OPTIONS["marchDist"] * PROC_COUNT_REL_TOL)

    def test3DQualityStatsIndependentOfProcCount(self):
        """
        Given the even sphere example surface moved with setSurfaceCoordinates,
        when it is extruded once partitioned across all procs and once on a single proc,
        then both runs report the same minimum cell quality and minimum cell volume.
        """
        # Any change to the surface coordinates will do, a translation leaves the cell quality and volume unchanged
        surfaceOffset = numpy.array([1.0, 2.0, 3.0])
        parallelHyp, serialHyp, _, _ = self.extrudeSphereInParallelAndSerial(surfaceOffset)

        # The min quality and volume are only reduced to rank 0
        if serialHyp is not None:
            parallelData = parallelHyp.hyp.hypdata
            serialData = serialHyp.hyp.hypdata
            self.assertAlmostEqual(
                parallelData.minqualityoverall, serialData.minqualityoverall, delta=PROC_COUNT_REL_TOL
            )
            self.assertAlmostEqual(
                parallelData.minvolumeoverall,
                serialData.minvolumeoverall,
                delta=abs(serialData.minvolumeoverall) * PROC_COUNT_REL_TOL,
            )

    def testSimpleOCart(self):
        from examples.simpleOCart.runSimpleOCart import extrudeDefaultCase

        outFile, hExtra = extrudeDefaultCase()

        self.commonTest(outFile, hExtra)

    def testSimpleOCartNoSurfaceMesh(self):
        from examples.simpleOCart.runSimpleOCart import extrudeNoSurfaceMeshCase

        outFile, hExtra = extrudeNoSurfaceMeshCase()

        self.commonTest(outFile, hExtra)

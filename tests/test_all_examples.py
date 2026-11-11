import os
import sys
import unittest
from mpi4py import MPI

baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "ref")

# The number of procs the reference meshes are generated with
REF_N_PROCS = 2

# The min cell quality pyHyp reports for the 717 example
REF_717_MIN_QUALITY = 0.2898473


class ExampleTests(unittest.TestCase):
    """
    Tests that run each example and compare its mesh to the reference mesh.

    These are run on different numbers of procs by the TestCase subclasses at the bottom of this file, which all compare
    to the same reference meshes.
    """

    def setUp(self):
        # Insert the repo root in path so that testflo can import examples
        sys.path.insert(0, os.path.join(baseDir, "../"))

    def commonTest(self, testFile, marchDist, relTol=1e-14, procCountRelTol=1e-6):
        """
        Compare testFile to the reference mesh with the same name, with absolute tolerance relTol * marchDist.

        When running on a different number of procs than the reference meshes were generated with, the tolerance is
        loosened to at least procCountRelTol. Results on different numbers of procs don't match to round-off, probably
        because the ASM preconditioner's subdomains depend on the partition and the linear solves stop at KSPRelTol.
        """
        if self.N_PROCS != REF_N_PROCS:
            relTol = max(relTol, procCountRelTol)

        volumeName = os.path.split(testFile)[1]
        refFile = os.path.join(refDir, volumeName)

        # Set the absolute tolerance for float comparisons based on mesh dimensions
        absTol = marchDist * relTol

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

        # This example deforms the surface with setSurfaceCoordinates, so also check that the reported min cell quality
        # is computed from up to date coordinates. It is only reduced to the root proc.
        if MPI.COMM_WORLD.rank == 0:
            self.assertAlmostEqual(hyp.hyp.hypdata.minqualityoverall, REF_717_MIN_QUALITY, delta=1e-6)

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
        # The corner amplifies differences in the linear solves, so serial and parallel results are further apart
        self.commonTest(volumeFile, marchDist, relTol=1e-9, procCountRelTol=1e-3)

    def testM6(self):
        from examples.m6.runM6 import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        # This example uses a loose KSPRelTol, so serial and parallel results are further apart
        self.commonTest(volumeFile, marchDist, procCountRelTol=1e-3)

    def testPlate(self):
        from examples.plate.runPlate import volumeFile, hyp

        marchDist = hyp.getUsedMarchDistance()
        self.commonTest(volumeFile, marchDist)

    def testSphere(self):
        from examples.sphere.runSphere import volumeFile, commonOptions

        marchDist = commonOptions["marchDist"]
        self.commonTest(volumeFile, marchDist)

    def testSimpleOCart(self):
        from examples.simpleOCart.runSimpleOCart import extrudeDefaultCase

        outFile, hExtra = extrudeDefaultCase()

        self.commonTest(outFile, hExtra)

    def testSimpleOCartNoSurfaceMesh(self):
        from examples.simpleOCart.runSimpleOCart import extrudeNoSurfaceMeshCase

        outFile, hExtra = extrudeNoSurfaceMeshCase()

        self.commonTest(outFile, hExtra)


class TestExamplesSerial(ExampleTests):
    N_PROCS = 1


class TestExamplesParallel(ExampleTests):
    N_PROCS = REF_N_PROCS


# Remove the base class from the module so that test discovery only runs the subclasses
del ExampleTests

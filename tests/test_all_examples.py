import os
import sys
import unittest
from mpi4py import MPI

baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "ref")


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
        self.commonTest(volumeFile, marchDist)

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

    def testSimpleOCart(self):
        from examples.simpleOCart.runSimpleOCart import extrudeDefaultcase

        outFile, hExtra = extrudeDefaultcase()

        self.commonTest(outFile, hExtra)

    def testSimpleOCartNoSurfaceMesh(self):
        from examples.simpleOCart.runSimpleOCart import extrudeNoSurfaceMeshCase

        outFile, hExtra = extrudeNoSurfaceMeshCase()

        self.commonTest(outFile, hExtra)

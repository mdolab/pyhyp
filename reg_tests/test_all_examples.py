import os
import subprocess
import sys
import unittest

baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "ref")


class TestExamples(unittest.TestCase):

    N_PROCS = 2

    def setUp(self):
        # Insert the repo root in path so that testflo can import examples
        sys.path.insert(0, os.path.join(baseDir, "../"))

    def common_test(self, testFile, marchDist, relTol=1e-14):
        volumeName = os.path.split(testFile)[1]
        refFile = os.path.join(refDir, volumeName)

        # Set the absolute tolerance for float comparisons based on mesh dimensions
        absTol = marchDist * relTol

        # Run cgnsdiff and store the terminal output
        output = subprocess.getoutput(f"cgnsdiff -d -t {absTol} {testFile} {refFile}")

        # Assert that there is no diff
        try:
            self.assertEqual(output, "")
        # Or that only the CGNS version differs
        except AssertionError:
            self.assertEqual(output, "/CGNSLibraryVersion <> /CGNSLibraryVersion : data values differ")

    def test_2D_euler(self):
        from examples.naca0012.naca0012_euler import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_2D_rans(self):
        from examples.naca0012.naca0012_rans import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_717(self):
        from examples.wing_717.run717 import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_BWB(self):
        from examples.BWB.runBWB import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_corner(self):
        from examples.corner.runCorner import volumeFile, commonOptions

        marchDist = commonOptions["marchDist"]
        self.common_test(volumeFile, marchDist)

    def test_M6(self):
        from examples.m6.runM6 import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_plate(self):
        from examples.plate.runPlate import volumeFile, hyp

        marchDist = hyp.getOption("marchDist")
        self.common_test(volumeFile, marchDist)

    def test_sphere(self):
        from examples.sphere.runSphere import volumeFile, commonOptions

        marchDist = commonOptions["marchDist"]
        self.common_test(volumeFile, marchDist)

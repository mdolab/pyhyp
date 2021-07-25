import os
import unittest
from pyhyp import pyHypMulti
from reg_test_utils import run_common_test

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExtrusion(unittest.TestCase):

    N_PROCS = 2

    def setUp(self):

        # Define file paths
        testDir = "sphere"
        surfaceFile1 = os.path.join(baseDir, testDir, "even_sphere.fmt")
        surfaceFile2 = os.path.join(baseDir, testDir, "uneven_sphere.fmt")
        surfaceFile3 = os.path.join(baseDir, testDir, "uneven_sphere_large.fmt")
        self.volumeFile = os.path.join(baseDir, testDir, "combined.cgns")
        self.refFile = os.path.join(baseDir, "ref", "sphere.ref")

        # Define common options
        commonOptions = {
            # ---------------------------
            #        Input Parameters
            # ---------------------------
            "inputFile": surfaceFile1,
            "unattachedEdgesAreSymmetry": True,
            "outerFaceBC": "farfield",
            "autoConnect": True,
            "BC": {},
            "families": "wall",
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            "N": 73,
            "s0": 1e-4,
            "marchDist": 10.0,
            # ---------------------------
            #   Pseudo Grid Parameters
            # ---------------------------
            "ps0": -1.0,
            "pGridRatio": -1.0,
            "cMax": 5.0,
            # ---------------------------
            #   Smoothing parameters
            # ---------------------------
            "epsE": 1.0,
            "epsI": 2.0,
            "theta": 3.0,
            "volCoef": 0.25,
            "volBlend": 0.0001,
            "volSmoothIter": 100,
        }

        # Define specific options that will overwrite the common options
        options1 = {"inputFile": surfaceFile1}
        options2 = {"inputFile": surfaceFile2}
        options3 = {"inputFile": surfaceFile3}

        # Gather specific options in a list
        options = [options1, options2, options3]

        # Run multiple extrusions and combine them into one file
        hyp = pyHypMulti(options=options, commonOptions=commonOptions)
        hyp.combineCGNS(combinedFile=self.volumeFile)

    def test_sphere(self, train=False):
        run_common_test(self.volumeFile, self.refFile, train)

    def train_sphere(self):
        self.test_sphere(train=True)

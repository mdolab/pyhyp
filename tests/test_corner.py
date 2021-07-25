"""
This test shows how to set up dictionaries to run multiple extrusions with pyHypMulti.
We use pyHypMulti to extrude a 90 deg corner twice with different smoothing settings.
"""


import os
import unittest
from pyhyp import pyHypMulti
from reg_test_utils import run_common_test

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExtrusion(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):

        # Define file paths
        testDir = "corner"
        surfaceFile = os.path.join(baseDir, testDir, "corner.cgns")
        self.volumeFile = os.path.join(baseDir, testDir, "combined.cgns")
        self.refFile = os.path.join(baseDir, "ref", "corner.ref")

        # Define common options
        commonOptions = {
            # ---------------------------
            #        Input Parameters
            # ---------------------------
            "inputFile": surfaceFile,
            "fileType": "CGNS",
            "unattachedEdgesAreSymmetry": False,
            "outerFaceBC": "farfield",
            "autoConnect": True,
            "BC": {},
            "families": "wall",
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            "N": 65,
            "s0": 1e-6,
            "marchDist": 2.5,
            "splay": 0.5,
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

        # Define specific options
        options1 = {"outputFile": "corner1_hyp.cgns"}
        options2 = {"epsE": 4.0, "epsI": 8.0, "outputFile": "corner2_hyp.cgns"}

        # Gather specific options in a list
        options = [options1, options2]

        # Run multiple extrusions and combine them into one file
        hyp = pyHypMulti(options=options, commonOptions=commonOptions)
        hyp.combineCGNS(combinedFile=self.volumeFile)

    def test_corner(self, train=False):
        run_common_test(self.volumeFile, self.refFile, train)

    def train_corner(self):
        self.test_corner(train=True)

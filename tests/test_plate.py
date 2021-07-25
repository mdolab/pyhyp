import os
import unittest
from pyhyp import pyHyp
from reg_test_utils import run_common_test

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExtrusion(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):

        # Define file paths
        testDir = "plate"
        surfaceFile = os.path.join(baseDir, testDir, "plate_surf.cgns")
        self.volumeFile = os.path.join(baseDir, testDir, "plate_vol.cgns")
        self.refFile = os.path.join(baseDir, "ref", "plate.ref")

        options = {
            # ---------------------------
            #        Input Parameters
            # ---------------------------
            "inputFile": surfaceFile,
            "fileType": "CGNS",
            "unattachedEdgesAreSymmetry": False,
            "outerFaceBC": "overset",
            "autoConnect": True,
            "BC": {1: {"jLow": "XYConst", "iLow": "XConst", "iHigh": "XConst"}},
            "families": "wall",
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            "N": 67,
            "s0": 1e-6,
            "nTruncate": 65,
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

        # Extrude the volume mesh with pyHyp
        hyp = pyHyp(options=options)
        hyp.run()
        hyp.writeCGNS(self.volumeFile)

    def test_plate(self, train=False):
        run_common_test(self.volumeFile, self.refFile, train)

    def train_plate(self):
        self.test_plate(train=True)

import os
import unittest
from pyhyp import pyHyp
from reg_test_utils import run_common_test

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExtrusion(unittest.TestCase):

    N_PROCS = 2

    def setUp(self):

        # Define file paths
        testDir = "BWB"
        name = "bwb"
        surfaceFile = os.path.join(baseDir, testDir, f"{name}.fmt")
        self.volumeFile = os.path.join(baseDir, testDir, f"{name}.cgns")
        self.refFile = os.path.join(baseDir, "ref", f"{name}.ref")

        options = {
            # ---------------------------
            #        Input Parameters
            # ---------------------------
            "inputFile": surfaceFile,
            "unattachedEdgesAreSymmetry": True,
            "outerFaceBC": "farfield",
            "autoConnect": True,
            "BC": {},
            "families": "wall",
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            "N": 81,
            "s0": 4e-6,
            "marchDist": 1100.0,
            # ---------------------------
            #   Pseudo Grid Parameters
            # ---------------------------
            "ps0": -1.0,
            "pGridRatio": -1.0,
            "cMax": 2.5,
            # ---------------------------
            #   Smoothing parameters
            # ---------------------------
            "epsE": 1.0,
            "epsI": 2.0,
            "theta": 3.0,
            "volCoef": 0.25,
            "volBlend": 0.0002,
            "volSmoothIter": 150,
        }

        # Extrude the volume mesh with pyHyp
        hyp = pyHyp(options=options)
        hyp.run()
        hyp.writeCGNS(self.volumeFile)

    def test_bwb(self, train=False):
        run_common_test(self.volumeFile, self.refFile, train)

    def train_bwb(self):
        self.test_bwb(train=True)

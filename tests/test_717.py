import os
import unittest
import numpy
from pyhyp import pyHyp
from pygeo import DVGeometry
from pyspline import Curve
from reg_test_utils import run_common_test

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExtrusion(unittest.TestCase):

    N_PROCS = 2

    def setUp(self):

        # Define file paths
        testDir = "717"
        surfaceFile = os.path.join(baseDir, testDir, "717_small.fmt")
        ffdFile = os.path.join(baseDir, testDir, "mdo_tutorial_ffd.fmt")
        self.volumeFile = os.path.join(baseDir, testDir, "717.cgns")
        self.refFile = os.path.join(baseDir, "ref", "717.ref")

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
            "s0": 5e-6,
            "marchDist": 325.0,
            # ---------------------------
            #   Pseudo Grid Parameters
            # ---------------------------
            "ps0": -1.0,
            "pGridRatio": -1.0,
            "cMax": 4.0,
            # ---------------------------
            #   Smoothing parameters
            # ---------------------------
            "epsE": 3.0,
            "epsI": 6.0,
            "theta": 3.0,
            "volCoef": 0.25,
            "volBlend": 0.0005,
            "volSmoothIter": 250,
        }

        # Create the pyHyp object
        hyp = pyHyp(options=options)

        # Define winglet deformation
        def winglet(val, geo):
            # Move the last point on ref axis
            C = geo.extractCoef("wing")

            C[-1, 0] += val[0]
            C[-1, 1] += val[1]
            C[-1, 2] += val[2]

            # Also need to get "dihedral" angle for this section
            theta = numpy.arctan(val[1] / (C[-1, 2] - C[-2, 2]))
            geo.rot_x["wing"].coef[-1] = -theta * 180 / numpy.pi
            geo.restoreCoef(C, "wing")

            # Set the chord as well
            geo.scale["wing"].coef[-1] = val[3]

        # Set up geometry object
        coords = hyp.getSurfaceCoordinates()
        DVGeo = DVGeometry(ffdFile)
        coef = DVGeo.FFD.vols[0].coef.copy()

        # First determine the reference chord lengths
        nSpan = coef.shape[2]
        ref = numpy.zeros((nSpan, 3))

        for k in range(nSpan):
            max_x = numpy.max(coef[:, :, k, 0])
            min_x = numpy.min(coef[:, :, k, 0])

            ref[k, 0] = min_x + 0.25 * (max_x - min_x)
            ref[k, 1] = numpy.average(coef[:, :, k, 1])
            ref[k, 2] = numpy.average(coef[:, :, k, 2])

        # Apply winglet deformations
        c0 = Curve(X=ref, k=2)
        DVGeo.addRefAxis("wing", c0)
        DVGeo.addGlobalDV("winglet", [0, 0, 0, 1], winglet, lower=-5, upper=5)
        DVGeo.addPointSet(coords, "coords")
        DVGeo.setDesignVars({"winglet": [1.5, 2.5, -2.0, 0.60]})
        hyp.setSurfaceCoordinates(DVGeo.update("coords"))

        # Extrude the volume mesh with pyHyp
        hyp.run()
        hyp.writeCGNS(self.volumeFile)

    def test_717(self, train=False):
        run_common_test(self.volumeFile, self.refFile, train)

    def train_717(self):
        self.test_717(train=True)

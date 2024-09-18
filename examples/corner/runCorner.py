"""
This example shows how to set up dictionaries to run multiple extrusions with pyHypMulti.
We use pyHypMulti to extrude a 90 deg corner twice with different smoothing settings.
"""

import os
from pyhyp import pyHypMulti

baseDir = os.path.dirname(os.path.abspath(__file__))
surfaceFile = os.path.join(baseDir, "corner.cgns")
volumeFile = os.path.join(baseDir, "corner_vol.cgns")

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

# Now set up specific options
options1 = {"outputFile": "corner1_hyp.cgns"}
options2 = {"epsE": 4.0, "epsI": 8.0, "outputFile": "corner2_hyp.cgns"}
options3 = {"splay": [[0.0, 0.5], [0.5, 0.0], [1.0, 0.5]], "outputFile": "corner3_hyp.cgns"}

# Gather options in a list
options = [options1, options2, options3]

hyp = pyHypMulti(options=options, commonOptions=commonOptions)
hyp.combineCGNS(combinedFile=volumeFile)

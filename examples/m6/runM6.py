import os
from pyhyp import pyHyp

baseDir = os.path.dirname(os.path.abspath(__file__))
surfaceFile = os.path.join(baseDir, "m6_small.fmt")
volumeFile = os.path.join(baseDir, "m6.cgns")

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": surfaceFile,
    "fileType": "PLOT3D",
    "unattachedEdgesAreSymmetry": True,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 81,
    "s0": 1.5e-5,
    "marchDist": 30.0,
    "nConstantStart": 1,
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
    "volBlend": 0.0005,
    "volSmoothIter": 100,
    "kspreltol": 1e-4,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS(volumeFile)

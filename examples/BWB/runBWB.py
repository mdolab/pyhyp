from pyhyp import pyHyp

fileName = "bwb.fmt"

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": fileName,
    "unattachedEdgesAreSymmetry": True,
    "outerFaceBC": "farField",
    "autoConnect": True,
    "BC": {},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 81,
    "s0": 4e-6,
    "marchDist": 1100,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1,
    "pGridRatio": -1,
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

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS("bwb.cgns")

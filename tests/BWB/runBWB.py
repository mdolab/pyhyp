# rst import (start)
from pyhyp import pyHyp

# rst import (end)

fileName = "bwb.fmt"

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": fileName,
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

# rst object
hyp = pyHyp(options=options)
# rst run
hyp.run()
hyp.writeCGNS("bwb.cgns")

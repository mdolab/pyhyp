import os
from pyhyp import pyHypMulti

baseDir = os.path.dirname(os.path.abspath(__file__))
surfaceFile1 = os.path.join(baseDir, "even_sphere.fmt")
surfaceFile2 = os.path.join(baseDir, "uneven_sphere.fmt")
surfaceFile3 = os.path.join(baseDir, "uneven_sphere_large.fmt")
volumeFile = os.path.join(baseDir, "sphere.cgns")

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

# Set up specific options that will overwrite the commonOptions

options1 = {"inputFile": surfaceFile1, "skip": False}
options2 = {"inputFile": surfaceFile2, "skip": False}
options3 = {"inputFile": surfaceFile3, "skip": False}

# Create list with different options
options = [options1, options2, options3]

# set skip list, if any
skipList = []

hyp = pyHypMulti(options=options, commonOptions=commonOptions, skipList=skipList)
hyp.combineCGNS(combinedFile=volumeFile)

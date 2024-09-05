"""
This script uses the NACA 0012 airfoil equation to generate a 2D RANS mesh.
This mesh has a blunt trailing edge.
"""

import os
import numpy
from pyhyp import pyHyp

baseDir = os.path.dirname(os.path.abspath(__file__))
surfaceFile = os.path.join(baseDir, "naca0012_rans.fmt")
volumeFile = os.path.join(baseDir, "naca0012_rans.cgns")

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": surfaceFile,
    "unattachedEdgesAreSymmetry": False,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 129,
    "s0": 1e-6,
    "marchDist": 100.0,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1.0,
    "pGridRatio": -1.0,
    "cMax": 3.0,
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


def generate_surface_file(surface_file):
    alpha = numpy.linspace(0, 2 * numpy.pi, 273)
    x = numpy.cos(alpha) * 0.5 + 0.5
    y = numpy.zeros_like(x)

    for i in range(len(x)):
        if i < len(x) / 2:
            y[i] = 0.6 * (
                0.2969 * numpy.sqrt(x[i]) - 0.1260 * x[i] - 0.3516 * x[i] ** 2 + 0.2843 * x[i] ** 3 - 0.1015 * x[i] ** 4
            )
        else:
            y[i] = -0.6 * (
                0.2969 * numpy.sqrt(x[i]) - 0.1260 * x[i] - 0.3516 * x[i] ** 2 + 0.2843 * x[i] ** 3 - 0.1015 * x[i] ** 4
            )

    # Since the TE is open we need to close it. Close it multiple linear segments.
    delta_y = numpy.linspace(y[-1], y[0], 32, endpoint=True)
    delta_y = delta_y[1:]

    x = numpy.append(x, numpy.ones_like(delta_y))
    y = numpy.append(y, delta_y)

    # Write the plot3d input file:
    f = open(surface_file, "w")
    f.write("1\n")
    f.write("%d %d %d\n" % (len(x), 2, 1))
    for iDim in range(3):
        for j in range(2):
            for i in range(len(x)):
                if iDim == 0:
                    f.write("%g\n" % x[i])
                elif iDim == 1:
                    f.write("%g\n" % y[i])
                else:
                    f.write("%g\n" % (float(j)))
    f.close()


def extrude_volume_mesh(options, volumeFile):
    hyp = pyHyp(options=options)
    hyp.run()
    hyp.writeCGNS(volumeFile)

    return hyp


if __name__ == "__main__":
    generate_surface_file(surfaceFile)
    extrude_volume_mesh(options, volumeFile)

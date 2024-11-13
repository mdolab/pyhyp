"""
This script uses the NACA 0012 airfoil equation to generate a 2D RANS mesh.
This mesh has a blunt trailing edge.
"""

import os
import numpy
from pyhyp import pyHyp
import numpy as np
import argparse

baseDir = os.path.dirname(os.path.abspath(__file__))
surfaceFile = os.path.join(baseDir, "naca0012_rans.fmt")

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


def extrude_default_case():
    """
    This is the default where most values are scalars
    """

    volumeFile = os.path.join(baseDir, "naca0012_rans.cgns")

    generate_surface_file(surfaceFile)
    hyp = extrude_volume_mesh(options, volumeFile)

    return hyp, volumeFile


def extrude_constant_layers_case():
    """
    Here, the first and last layers are kept constant (growth-ratio == 1.0)
    """

    volumeFile = os.path.join(baseDir, "naca0012_rans_constant_layers.cgns")

    options.update(
        {
            "nConstantStart": 5,
            "nConstantEnd": 5,
        }
    )

    generate_surface_file(surfaceFile)
    hyp = extrude_volume_mesh(options, volumeFile)

    return hyp, volumeFile


def extrude_schedule_case():
    """
    Some variables are 'scheduled' which means their value changes depending on
    the extrusion layer. The values are specified on intervals which are
    linearly interpolated by pyHyp.
    """
    volumeFile = os.path.join(baseDir, "naca0012_rans_schedule.cgns")

    options.update(
        {
            "epsE": [[0.0, 1.0], [0.2, 2.0], [1.0, 5.0]],
            "epsI": [[0.0, 2.0], [0.2, 2.0], [1.0, 10.0]],
            "theta": [[0.0, 3.0], [0.2, 2.5], [1.0, 0.0]],
            "volBlend": [[0.0, 0.0001], [1.0, 0.1]],
            "volSmoothIter": [[0.0, 100], [1.0, 500]],
            "volCoef": [[0.0, 0.25], [1.0, 0.5]],
            "growthRatios": [[0.0, 1.05], [1.0, 1.1]],
            "cornerAngle": [[0.0, 110.0], [1.0, 120.0]],
        }
    )

    generate_surface_file(surfaceFile)
    hyp = extrude_volume_mesh(options, volumeFile)

    return hyp, volumeFile


def extrude_explicit_case():
    """
    Some variables are set 'explicitly'. This means, a list of values that
    correspond to each layer is provided.
    """

    def ls(start, stop, dtype=np.float64):
        return np.linspace(start, stop, options["N"] - 1, dtype=dtype).tolist()

    options.update(
        {
            "epsE": ls(1.0, 5.0),
            "epsI": ls(2.0, 10.0),
            "theta": ls(3.0, 0.0),
            "volBlend": ls(0.0001, 0.1),
            "volSmoothIter": ls(100, 500, dtype=np.int32),
            "volCoef": ls(0.25, 0.5),
            "growthRatios": ls(1.05, 1.3),
            "cornerAngle": ls(110.0, 120.0),
        }
    )

    volumeFile = os.path.join(baseDir, "naca0012_rans_explicit.cgns")

    generate_surface_file(surfaceFile)
    hyp = extrude_volume_mesh(options, volumeFile)
    return hyp, volumeFile


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
    parser = argparse.ArgumentParser(description="Process some integers.")
    choices = ["default", "schedule", "growth_ratios"]
    parser.add_argument("--case", choices=choices, default=choices[0])
    args = parser.parse_args()

    if args.case == "default":
        extrude_default_case()
    elif args.case == "schedule":
        extrude_schedule_case()
    else:
        extrude_explicit_case()

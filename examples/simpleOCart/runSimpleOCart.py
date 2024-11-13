import os
from pyhyp.utils import simpleOCart
import argparse

baseDir = os.path.dirname(os.path.abspath(__file__))


# Set the other inputs
dh = 0.05
hExtra = 10.0
nExtra = 49
sym = "y"
mgcycle = 3
userOptions = {"cMax": 5.0}


def extrudeDefaultcase():
    # Usually this would be a nearfield volume file
    # We use the corner surface file for convenience
    nearFile = os.path.join(baseDir, "../corner/corner.cgns")

    outFile = os.path.join(baseDir, "simpleOCart.cgns")

    # Run simpleOCart
    simpleOCart(nearFile, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=userOptions)

    return outFile, hExtra


def extrudeNoSurfaceMeshCase():
    outFile = os.path.join(baseDir, "simpleOCart_no_surface_mesh.cgns")

    xBounds = [[0, 0, 0], [2, 2, 2]]

    # Run simpleOCart
    simpleOCart(None, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=userOptions, xBounds=xBounds)

    return outFile, hExtra


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some integers.")
    choices = ["default", "noSurfaceMesh"]
    parser.add_argument("--case", choices=choices, default=choices[0])
    args = parser.parse_args()

    if args.case == "noSurfaceMesh":
        extrudeNoSurfaceMeshCase()
    else:
        extrudeDefaultcase()

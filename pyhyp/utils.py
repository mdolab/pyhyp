import os
import tempfile
import numpy as np
from mpi4py import MPI
from cgnsutilities.cgnsutilities import readGrid, Block
from .pyHyp import pyHyp


def simpleOCart(nearFile, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=None):
    """
    Generates a Cartesian mesh around the provided grid, surrounded by an O-mesh.

    Parameters
    ----------
    nearFile : str
        The name of the nearfield CGNS file to mesh around.

    dh : float or list of float
        The target edge length of each cell in the Cartesian part of the mesh.
        A list of (x, y, z) lengths can be provided to make non-cubic cells.
        The actual edge lengths will depend on mgcycle.

    hExtra : float
        The distance from the Cartesian mesh boundary to the farfield.

    nExtra : int
        The number of layers to extrude the hyperbolic O-mesh.

    sym : str or list of str
        Axis or plane of symmetry.
        One or more of ('x', 'y', 'z', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax').

    mgcycle : int or list of int
        Number of times mesh should be able to be coarsened for multigrid cycles.
        A list can be provided for nonuniform (x, y, z) coarsening.

    outFile : str
        Output file name.

    userOptions : dict, optional
        The name of the nearfield CGNS file to mesh around.

    """

    # Read the nearfield file
    grid = readGrid(nearFile)

    # First run simpleCart with no extension
    X, dx = grid.simpleCart(dh, 0.0, 0, sym, mgcycle, outFile=None)

    # Pull out the patches from the Cartesian mesh.
    # We have to pay attention to the symmetry and the ordering of the patches
    # to make sure that all the normals are pointing out.
    patches = []

    # First take patches that are opposite from the origin planes
    if "xmax" not in sym:
        patches.append(X[-1, :, :, :])
    if "ymax" not in sym:
        patches.append(X[:, -1, :, :][::-1, :, :])
    if "zmax" not in sym:
        patches.append(X[:, :, -1, :])

    # Then take patches from the origin planes
    if "x" not in sym and "xmin" not in sym:
        patches.append(X[0, :, :, :][::-1, :, :])
    if "y" not in sym and "ymin" not in sym:
        patches.append(X[:, 0, :, :])
    if "z" not in sym and "zmin" not in sym:
        patches.append(X[:, :, 0, :][::-1, :, :])

    # Set up the generic input for pyHyp
    hypOptions = {
        "patches": patches,
        "unattachedEdgesAreSymmetry": True,
        "outerFaceBC": "farfield",
        "autoConnect": True,
        "BC": {},
        "N": nExtra,
        "s0": np.average(dx),
        "marchDist": hExtra,
        "cmax": 3.0,
    }

    # Use user-defined options if provided
    if userOptions is not None:
        hypOptions.update(userOptions)

    # Run pyHyp
    hyp = pyHyp(options=hypOptions)
    hyp.run()

    fName = None
    if MPI.COMM_WORLD.rank == 0:
        dirpath = tempfile.mkdtemp()
        fName = os.path.join(dirpath, "tmp.cgns")

    hyp.writeCGNS(MPI.COMM_WORLD.bcast(fName))

    # Reset symmetry to single axis
    if "x" in sym or "xmin" in sym or "xmax" in sym:
        sym = "x"
    elif "y" in sym or "ymin" in sym or "ymax" in sym:
        sym = "y"
    elif "z" in sym or "zmin" in sym or "zmax" in sym:
        sym = "z"

    if MPI.COMM_WORLD.rank == 0:
        # Read the pyhyp mesh back in and add our additional "X" from above.
        grid = readGrid(fName)
        dims = X.shape[0:3]
        grid.addBlock(Block("interiorBlock", dims, X))
        grid.renameBlocks()
        grid.connect()
        grid.BCs = []
        grid.autoFarfieldBC(sym)
        grid.writeToCGNS(outFile)

        # Delete the temp file
        os.remove(fName)

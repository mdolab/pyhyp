import os
import tempfile
import numpy as np
from mpi4py import MPI
from baseclasses.utils import Error
from cgnsutilities.cgnsutilities import readGrid, Block, simpleCart, Grid
from .pyHyp import pyHyp


def simpleOCart(inputGrid, dh, hExtra, nExtra, sym, mgcycle, outFile, userOptions=None, xBounds=None, useFarfield=True):
    """
    Generates a Cartesian mesh around the provided grid, surrounded by an O-mesh.

    Parameters
    ----------
    inputGrid : cgnsutils Grid object or str
        If a cgnsutils Grid object is provided, we use it as is.
        Alternatively if a string is provided, we treat it as
        the name of the nearfield CGNS file to mesh around.

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
        Custom pyhyp options to be used with this extrusion. If overset BCs are desired
        on the outer face, do not set it in this dictionary because after extrusion
        we overwrite all BCs on the combined grid. See the option useFarfield
        below. The default value (True) results in farfield BCs on the outer face,
        and setting it to false results in overset for the far face. Other pyhyp extrusion
        parameters can be set here.

    xBounds : array (2 x 3), optional
        Optional bounding box coordinates desired for the center cartesian grid.
        The default value can be obtained by: ``xMin, xMax = grid.getBoundingBox()``,
        and then the xBounds array can be set as ``xBounds = [xMin, xMax]``. This option
        allows users to modify the bounding box coordinates rather than simply defaulting
        to the bounding box of the nearfield grid.

    useFarfield : bool, optional
        Optional flag to control the outermost layer's BC. Default, ``True``, will result
        in a farfield outer layer, setting this to ``False`` does overset BCs on the
        outermost layer

    """

    # check if we have a grid object as input or the filename
    if type(inputGrid) == str:
        # Read the nearfield file
        input_filename = inputGrid
        inputGrid = readGrid(input_filename)
    # the input Grid can be None, in this acse, we must have xBounds
    elif inputGrid is None:
        # make sure we have xbounds
        if xBounds is None:
            raise Error("If the inputGrid is None, xBounds must be provided")
    # if the input grid is not provided as a filename, it must be a Grid instance
    elif type(inputGrid) != Grid:
        # if not, raise an error
        raise Error(
            "The inputGrid to simpleOCart must either be the filename of the nearfield grid or a Grid type object from cgnsutilities."
        )

    if xBounds is None:
        # we will automatically determine the bounding box
        X, dx = inputGrid.simpleCart(dh, 0.0, 0, sym, mgcycle, outFile=None)
    else:
        # we are provided the bounding box, skip to the generic simple cart routine
        X, dx = simpleCart(xBounds[0], xBounds[1], dh, 0, 0, sym, mgcycle, outFile=None)

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
        simple_ocart_grid = readGrid(fName)
        dims = X.shape[0:3]
        simple_ocart_grid.addBlock(Block("interiorBlock", dims, X))
        simple_ocart_grid.renameBlocks()
        simple_ocart_grid.connect()
        simple_ocart_grid.BCs = []
        if useFarfield:
            simple_ocart_grid.autoFarfieldBC(sym)
        else:
            simple_ocart_grid.autoNearfieldBC(sym)
        simple_ocart_grid.writeToCGNS(outFile)

        # Delete the temp file
        os.remove(fName)

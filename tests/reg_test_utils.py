import os
from baseclasses import BaseRegTest
from cgnsutilities.cgnsutilities import readGrid

baseDir = os.path.dirname(os.path.abspath(__file__))


def run_common_test(volumeFile, refFile, train):

    # Check the volume mesh with cgnsUtilities
    grid = readGrid(volumeFile)
    nWallCells, nWallNodes = grid.getWallCellsNodes()
    blockInfo = grid.getBlockInfo()

    with BaseRegTest(refFile, train=train) as handler:
        handler.root_add_val("Wall cells", nWallCells, tol=0)
        handler.root_add_val("Wall nodes", nWallNodes, tol=0)
        handler.root_add_dict("Block info", blockInfo, tol=0)

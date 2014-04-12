#!/usr/bin/python
from __future__ import print_function
from __future__ import division
"""
pyHyp

The pyHyp module is used to generate a hyperbolic grid around a 3D
surface

Copyright (c) 2013 by G. Kenway
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 02/04/2013$

Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2013)
"""
# =============================================================================
# Imports
# =============================================================================
import sys
import os
import copy
import numpy
from mpi4py import MPI
from . import MExt

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a explicitly raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyHypEllip Error: '
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)
        Exception.__init__(self)

# =============================================================================
# pyHyp class
# =============================================================================
class pyHyp(object):
    """
    This is the main class pyHyp. It is used as the user interface to pyHyp. 
    """

    def __init__(self, comm=None, options=None, debug=False, **kwargs):
        """
        Create the pyHyp object. 

        Parameters
        ----------
        comm : MPI_INTRACOMM
            Comm to use. This is used when running in parallel. If not
            provided, MPI.COMM_WORLD is used by default.

        options : dict
            A dictionary containing the the options for pyHyp.

        debug : bool
            Flag used to specify if debugging. This only needs to be
            set to true when using a symbolic debugger. 
           
            """
            
        # Default options for hyperbolic generation
        self.options_default = {
            # ---------------------------
            #        Input Information
            # ---------------------------
            # One of '3d' or '2d'
            'dimension':'3d',

            # Input file
            'inputFile':'default',

            # Type of extrusion: hyperbolic or elliptic
            'mode':'hyperbolic',

            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            # Number of layers:
            'N': 100, 

            # Initial off-wall spacing
            's0':0.01, 
            
            # Rmin: Distance to march in multiples of (automatically
            # computed) initial radius surrounding the body. 
            'rMin': 50, 

            # nodeTol: Tolerance for nodes to be treated as identical.
            'nodeTol':1e-8,

            # mirror: Possible mirroring of mesh. Use 'x', 'y', 'z'
            # if grid needs to be mirrored. 
            'mirror':None,

            # panelEps: Distance source panels are "below" nodes. This
            # parameter usually doesn't need to be changed.
            'panelEps':1e-8,
            
            # ------------------------------------------
            #   Pseudo Grid Parameters (Hyperbolic only)
            # ------------------------------------------

            # Maximum permissible ratio of marching direction length
            # to smallest in-plane edge
            'cMax':6.0, 

            #nonLinear: True/False. Use nonlinear scheme. Not
            #currently working. Hyperbolic only
            'nonLinear': False, 

            #slExp: Exponent for Sl calc. Don't change this value
            #unless you know what you are doing! Hyperbolic only.
            'slExp': .15, 

            # Initial off-wall spacing. Hyperbolic only
            'ps0':0.01, 
            
            # Grid Spacing Ratio. 
            'pGridRatio':1.15, 

            # ----------------------------------------
            #   Smoothing parameters (Hyperbolic only)
            # ----------------------------------------

            # epsE: The explicit smoothing coefficient
            'epsE': 1.0, 

            # epsI: The implicit smoothing coefficient
            'epsI': 2.0, 

            # theta: The barth implicit smoothing coefficient
            'theta': 3.0, 

            # volCoef: The volume smoothing coefficient for
            # pointJacobi iterations
            'volCoef': .16, 

            # volBlend: The volume blending coefficient to force
            # uniform sizes in farfield
            'volBlend': 0.0001, 

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 10, 

            # -------------------------------
            #   Solution Parameters (Common)
            # ------------------------------
            # kspRelTol: Solution tolerance for linear system
            'kspRelTol': 1e-8, 
            
            # Maximum number of iterations to run for linear system
            'kspMaxIts': 500, 

            # Preconditioner Lag (Hyperbolic only)
            'preConLag': 10, 

            'kspSubspaceSize':50, 

            # ---------------------------
            #   Output Parameters
            # ---------------------------
            # Debugging option to write grid metrics. 
            'writeMetrics': False, 
            }
   
        # Import and set the hyp module
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.hyp = MExt.MExt('hyp', [curDir], debug=debug)._module
   
        # Set the possible MPI Intracomm
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm

        # Initialize PETSc and MPI if not already done so:
        self.hyp.initpetsc(self.comm.py2f())

        # Use supplied options
        if options is None:
            raise Error("The options= keyword argument is *NOT* optional. "
                        "It must always be provided")
        self.options = options

        # Setup the options
        self._checkOptions()
        self._setOptions()

        # Initialize the problem based on dimension
        if self.options['dimension'] == '2d':
            self._init2d()
        else:
            self._init3d()

        self.gridGenerated = False
        self.symNodes = None
        if self.comm.rank == 0:
            print('Problem successfully setup!')

    def _init2d(self, fileName, X, flip):

        # Check to see how the user has passed in the data:
        if fileName is None and X is None:
            raise Error("Either fileName or X must be passed on "
                        " initialization!")

        if fileName is not None and X is not None:
            raise Error("BOTH fileName or X have been passed on initialization!"
                        " Only one is required.")

        # Take the file name and try loading it.  Only do this on one
        # proc; we will MPI-bcast the information back to everyone
        # else when it is necessary.
        if self.comm.rank == 0:
            if fileName is not None:
                # Load the curve from a file using the the following format:
                # <number of points>
                # x1 y1
                # x2 y2 
                # ....
                # xn yn

                f = open(fileName, 'r')
                N = int(f.readline())
                x = []
                y = []
                for i in xrange(N):
                    aux = f.readline().split()
                    x.append(float(aux[0]))
                    y.append(float(aux[1]))

                # Done with the file so close
                f.close()

                # Convert the lists to an array
                X = numpy.array([x, y]).T

            # Now check if the first and the last point are within
            # tolerance of each other:
            if self._e_dist(X[0], X[-1]) < 1e-6:
                # Curve is closed, we're ok. Internally, we work
                # with the unclosed curve and explictly deal with
                # the  0th and Nth points 
                X = X[0:-1, :]
            else:
                # Curve is not closed, print warning
                print('*'*105)
                print(" Warning: Curve was not closed! Closing curve with a "
                      "linear segment. This may or not be what is desired!")
                print('*'*105)

            if flip: # Reverse direction
                X[:, 0] = X[:, 0][::-1]
                X[:, 1] = X[:, 1][::-1]

        # Now broadcast the required info to the other procs:
        self.X = self.comm.bcast(X, root=0)
        self.N = len(self.X)

    def _init3d(self, **kwargs):

        if 'xMirror' in kwargs or 'yMirror' in kwargs or 'zMirror' in kwargs:
            raise Error("Mirror specification must be now done in the "
                        "options dictionary using: 'mirror':'z' for "
                        "example.")

        # Only need to file i/o and processing on root proc:
        if self.comm.rank == 0:
            delFile = False
            fileName = self.options['inputFile']
            if self.options['mirror'] is not None:

                # We need to first read the file, mirror it, and write it
                # back to a tmp file such that the nominal read below
                # works.

                # Load in the file and determine the free edges, these
                # will be symmetry edges (which turn into symmetry planes)
                surfs, sizes, __, __, l_index = self._readPlot3D(fileName)
                nSurf = len(surfs)

                # To determine the nodes on the symmetry plane, we need to
                # figure out which edges have a face only on one
                # side. Finally, the end goal is to prodce an
                # (unstructured) list that looks like: 
                #
                # [[surfID_1, i_1, j_1], [surfID_2, i_2, j_2] ...]
                # 
                # where the i, j gives the index on surf that is
                # initialy on the mirror plane. This information can
                # then be used to zero this "line" to be precisely on
                # the mirror plane as a post processing step. This
                # method works even if a complete block edge is not on
                # a mirror plane.

                edges = {}
                def getKey(n1, n2):
                    if n1 < n2:
                        return (n1, n2)
                    else:
                        return (n2, n1)

                def incrementEdge(edges, n1, n2):
                    key = getKey(n1, n2)
                    if key in edges:
                        edges[key] += 1
                    else:
                        edges[key] = 1

                    return edges

                # Loop over faces and add edges with keys of (n1, n2)
                # where n1 < n2
                for ii in xrange(len(surfs)):
                    for i in xrange(sizes[ii][0]-1):
                        for j in xrange(sizes[ii][1]-1):
                            edges = incrementEdge(edges, l_index[ii][i  , j  ], l_index[ii][i+1, j  ])
                            edges = incrementEdge(edges, l_index[ii][i+1, j  ], l_index[ii][i+1, j+1])
                            edges = incrementEdge(edges, l_index[ii][i+1, j+1], l_index[ii][i  , j+1])
                            edges = incrementEdge(edges, l_index[ii][i  , j+1], l_index[ii][i  , j  ])

                # Now generate a final list of nodes that are on the
                # mirror plane:
                mirrorNodes = []
                for key in edges:
                    if edges[key] == 1:
                        mirrorNodes.append(key[0])
                        mirrorNodes.append(key[1])

                # Uniquify them
                uniqueMirrorNodes = numpy.unique(mirrorNodes)

                # Now loop back over the boundries of the
                # surfaces. For each node on the boundary check if it
                # is in uniqueMirror nodes.

                self.symNodes = []
                for ii in xrange(len(surfs)):
                    imax = sizes[ii][0]
                    jmax = sizes[ii][1]
                    for i in xrange(imax):
                        if l_index[ii][i, 0] in uniqueMirrorNodes:
                            self.symNodes.append([ii, i  , 0])
                        if l_index[ii][i, jmax-1] in uniqueMirrorNodes:
                            self.symNodes.append([ii, i,  jmax-1])

                    for j in xrange(jmax):
                        if l_index[ii][0, j] in uniqueMirrorNodes:
                            self.symNodes.append([ii, 0, j  ])
                        if l_index[ii][imax-1, j] in uniqueMirrorNodes:
                            self.symNodes.append([ii, imax-1, j  ])

                # Now that we know the symnodes, we can hard zero the
                # nodes on the mirror plane. Note that if your surface
                # isn't closed, some really really stuff will result
                # since those nodes will be zeroed!
                for iSym in xrange(len(self.symNodes)):
                    iSurf = self.symNodes[iSym][0]
                    i     = self.symNodes[iSym][1]
                    j     = self.symNodes[iSym][2]
                    if self.options['mirror'] == 'x':
                        surfs[iSurf][i, j, 0] = 0.0
                    if self.options['mirror'] == 'y':
                        surfs[iSurf][i, j, 1] = 0.0
                    if self.options['mirror'] == 'z':
                        surfs[iSurf][i, j, 2] = 0.0

                # Conver symnodes to array and index everything my 1 
                self.symNodes = numpy.array(self.symNodes) + 1

                # Note that we are guaranteed that symNodes is sorted by
                # patchID. This is important since it means we can loop
                # over the blocked in sequence and all nodes need to 

                # Generate new list of sizes (double the length)
                newSizes = numpy.zeros((nSurf*2, 2), 'intc')
                for i in xrange(nSurf):
                    newSizes[i] = sizes[i]
                    newSizes[i+nSurf] = sizes[i]

                # Now mirror each zone, while flipping the i and j index
                for i in xrange(nSurf):
                    surfs.append(numpy.zeros([newSizes[i+nSurf, 0],
                                              newSizes[i+nSurf, 1], 3]))
                    if self.options['mirror'] == 'x':
                        surfs[i+nSurf][:, :, 0] = -self._reverseRows(surfs[i][:, :, 0])
                    else:
                        surfs[i+nSurf][:, :, 0] = self._reverseRows(surfs[i][:, :, 0])

                    if self.options['mirror'] == 'y':
                        surfs[i+nSurf][:, :, 1] = -self._reverseRows(surfs[i][:, :, 1])
                    else:
                        surfs[i+nSurf][:, :, 1] = self._reverseRows(surfs[i][:, :, 1])

                    if self.options['mirror'] == 'z':
                        surfs[i+nSurf][:, :, 2] = -self._reverseRows(surfs[i][:, :, 2])
                    else:
                        surfs[i+nSurf][:, :, 2] = self._reverseRows(surfs[i][:, :, 2])

                # Dump back out
                f = open('tmp.fmt', 'w')
                f.write('%d\n'%(nSurf*2))
                for i in xrange(nSurf*2):
                    f.write('%d %d %d\n'%(newSizes[i][0], newSizes[i][1], 1))
                for ii in xrange(nSurf*2):
                    for idim in xrange(3):
                        for j in xrange(newSizes[ii][1]):
                            for i in xrange(newSizes[ii][0]):
                                f.write('%20.13g\n'%(surfs[ii][i, j, idim]))
                f.close()
                fileName = 'tmp.fmt'
                delFile = True

            # Now we read (the possibly) mirrored plot3d file
            surfs, sizes, nodes, conn, l_index = self._readPlot3D(fileName)
            f = open('node_test.dat','w')
            f.write('VARIABLES = \"X\", \"Y\", \"Z\" \n')
            f.write("Zone\n")
            for i in range(len(nodes)):
                f.write('%g %g %g\n'% (nodes[i][0], nodes[i][1], nodes[i][2]))
            f.close()
            
            nSurf = len(surfs)
            nGlobal = len(nodes)
      
            if delFile:
                os.remove(fileName)
                
            # Determine the maximum number of "multigrid" levels
            # possible.
            nLevels = self._determineMGLevels(sizes)

            # We really should check if the surface is closed. This check
            # is not implemented yet -- so user beware!

            # Next we need to produced an unstructured-like data
            # structure

            # Now we can invert and get the elements surrounding the nodes:
            nodeToElem = [[] for i in xrange(nGlobal)]
            for iElem in xrange(len(conn)):
                for ii in xrange(4): # 4 nodes on each face
                    # Append the elem Number and the node number of
                    # the node on theface
                    nodeToElem[conn[iElem][ii]].append([iElem, ii]) 

            # Finally we can get the final node pointer structure we need:
            nPtr = [[] for i in range(nGlobal)]
            cPtr = [[] for i in range(nGlobal)]
            for iNode in xrange(nGlobal):
                nFaceNeighbours = len(nodeToElem[iNode])

                # Regular nodes get 4 neightbours, others get double
                if nFaceNeighbours == 2:
                    print('Can\'t do geometries with nodes of valance 2!')
                    sys.exit(0)
                    nPtr[iNode].append(nFaceNeighbours*2)
                elif nFaceNeighbours == 3:
                    nPtr[iNode].append(nFaceNeighbours*2)
                else:
                    nPtr[iNode].append(nFaceNeighbours)
                cPtr[iNode].append(nFaceNeighbours)
                
                # Get the face where this node is node 0
                firstID  = nodeToElem[iNode][0][0]
                nodeID  = numpy.mod(nodeToElem[iNode][0][1]+1, 4)
                nextElemID = None

                while nextElemID != firstID:
                    if nextElemID is None:
                        nextElemID = firstID

                    # Append the next node along that edge:
                    nodeToAdd = conn[nextElemID][numpy.mod(nodeID, 4)]
                    nPtr[iNode].append(nodeToAdd)
                    cPtr[iNode].append(nextElemID)
                    # Add the diagonal if not regular node
                    if nFaceNeighbours == 3:
                        nPtr[iNode].append(conn[nextElemID][
                            numpy.mod(nodeID+1, 4)])

                    # And we need to find the face that contains the
                    # following node:
                    nodeToFind = conn[nextElemID][numpy.mod(nodeID + 2, 4)]

                    found = False
                    jj = -1
                    while not found:
                        jj = jj + 1
                        # Find the next face
                        faceToCheck = nodeToElem[iNode][jj][0]
                        if faceToCheck != nextElemID:
                            # Check the 4 nodes on this face 
                            for ii in xrange(4):
                                if conn[faceToCheck][ii] == nodeToFind:
                                    nextElemID = faceToCheck
                                    nodeID = ii
                                    found = True
                                # end if
                            # end for
                        # end if
                    # end while
                # end if
            # end for

            # Figure out how big we need to make the array for fortran:
            lMax = 0
            for i in xrange(len(nPtr)):
                lMax = max(lMax, len(nPtr[i]))

            tmp = numpy.zeros((len(nPtr), lMax), 'intc')
            for i in xrange(len(nPtr)):
                l = len(nPtr[i])
                tmp[i, 0] = nPtr[i][0]
                tmp[i, 1:l] = numpy.array(nPtr[i][1:]) 
            nPtr = tmp
                           
            # Generate the dual panel mesh for the ellipic method:
            dualPanels = [[] for i in range(nGlobal)]
            for i in range(nGlobal):
                for j in range(cPtr[i][0]):
                    iPanel = cPtr[i][j+1]
                    # Now get the center of this panel:
                    xCen = numpy.zeros(3)
                    for ii in range(4):
                        xCen += nodes[conn[iPanel][ii]]
                    dualPanels[i].append(xCen/4)
            
            # Now compute the centers of DUAL panels
            dualCenters = numpy.zeros((nGlobal, 3))
            xCen = numpy.zeros(3)
            for i in range(nGlobal):
                xCen[:] = 0.0
                nn = cPtr[i][0]
                for j in range(cPtr[i][0]):
                    dualCenters[i] += dualPanels[i][j]/nn
        
            # Now compute the normals of the dual panels:
            normals = numpy.zeros((nGlobal, 3))
            areas = numpy.zeros(nGlobal)

            for i in range(nGlobal):
                if cPtr[i][0] == 3:
                    v1 = dualPanels[i][1] - dualPanels[i][0]
                    v2 = dualPanels[i][2] - dualPanels[i][0]
                    A = 0.5*numpy.linalg.norm(numpy.cross(v1, v2))
                elif cPtr[i][0] == 4:
                    v1 = dualPanels[i][2] - dualPanels[i][0]
                    v2 = dualPanels[i][3] - dualPanels[i][1]
                    A = numpy.linalg.norm(numpy.cross(v1, v2))

                n = numpy.cross(v1, v2)
                n /= numpy.linalg.norm(n)
                normals[i] = n
                areas[i] = A
                # Now "Correct" the dual panels such that they lie
                # "eps" BELOW the actual nodes, along the normal
                eps = self.options['panelEps']
                r = nodes[i] - dualCenters[i] # Displacement vector
                for ii in range(cPtr[i][0]):
                    dualPanels[i][ii] += r - eps*n
            dualPanels = None
            normals = None
            areas = None
        else:
            nPtr = None
            conn = None
            nodes = None
            dualPanels = None
            normals = None
            areas = None
            sizes = None
            l_index = None
            
        # Broadcast the result to everyone:
        (sizes, l_index, nPtr, conn, self.X, dualPanels, normals, areas) = (
            self.comm.bcast((sizes, l_index, nPtr, numpy.array(conn), nodes, 
                             dualPanels, normals, areas)))
        nGlobal = len(nPtr)
        nSurf = len(sizes)
        # Conn needs to be saved for potential calling of writeFEMesh
        self.conn = numpy.array(conn, 'intc')

        # Now we can divy up the nodes -- as evenly as possible.
        evenNodes = nGlobal//self.comm.size
        istart = self.comm.rank*evenNodes
        iend =  istart + evenNodes
        if self.comm.rank == self.comm.size-1:
            iend = nGlobal
        isize = iend - istart
        
        # Partition here:
        self.X = self.X[istart:iend]
        nPtr = nPtr[istart:iend]

        # Now get the ghost information we need:
        ghost = []
        ii = 0
        lnPtr = nPtr.copy()
        for i in range(len(lnPtr)):
            for j in range(lnPtr[i][0]):
                if lnPtr[i][j+1] < istart or lnPtr[i][j+1] >= iend:
                    ghost.append(lnPtr[i][j+1])
                    lnPtr[i, j+1] = isize + ii
                    ii += 1
                else:
                    lnPtr[i, j+1] -= istart
      
        # Now nPtr is the indexing in the global ording. lnPtr uses
        # the local ghosted form.  Convert both to fortran.
        nPtr[:, 1:] += 1
        lnPtr[:, 1:] += 1

        # Now we have all the data we need so we can go ahead and
        # initialze the 3D generation
        self.hyp.hypdata.nx = len(self.X)
        self.hyp.hypdata.nxglobal = nGlobal
        self.hyp.hypdata.xsurf = self.X.T
        self.hyp.hypdata.gnptr = nPtr.T
        self.hyp.hypdata.lnptr = lnPtr.T
        if len(ghost) > 0:
            self.hyp.hypdata.ghost = ghost
        self.hyp.hypdata.nghost = len(ghost)

        self.hyp.init3d(sizes)
        for i in xrange(nSurf):
            self.hyp.setlindex(l_index[i], i+1) # +1 for 0->1 ordering

        # # Now set the panel data individually...a little slow we know
        # for i in range(nGlobal):
        #     self.hyp.setpaneldata(i+1, numpy.array(dualPanels[i]).T,
        #                           normals[i], areas[i])

    def run(self):
        """
        Run given using the options given
        """
        if self.options['dimension'] == '2d':
            self.hyp.run2d(self.X.T)
        else:
            self.hyp.run3d()
        self.gridGenerated = True

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""
        
        if self.gridGenerated:
            if self.options['dimension'] == '2d':
                self.hyp.writeplot3d_2d(fileName)
            else:
                self.hyp.writeplot3d_3d(fileName)
        else:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly 
        formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if self.gridGenerated:
            if self.options['dimension'] == '2d':
                self.hyp.writecgns_2d(fileName)
            else:
                self.hyp.writecgns_3d(fileName)
                # Possibly zero mirror plane for mirrored geometries
                if self.symNodes is not None:
                    mirrorDims = {'x':1,'y':2,'z':3}
                    mirrorDim = mirrorDims[self.options['mirror']]
                    self.hyp.zeromirrorplane(fileName, self.symNodes, mirrorDim)
        else:
            raise Error("No grid has been generated! Run the run() "
                        "command before trying to write the grid!")

    def writeFEMesh(self, fileName):
        """ Ouput a tecplot FE mesh of the surface. Useful for
        debugging numberings.

        Parameters
        ----------
        fileName : str
            Filename of tecplot file. Should have .dat extension.
        """

        nNodes = len(self.X)
        nElem = len(self.conn)
        f = open(fileName, 'w')
        f.write("FE Data\n")
        f.write('VARIABLES = \"X\", \"Y\", \"z\" \n')
        f.write("ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, "
                "ZONETYPE=FEQUADRILATERAL\n"%(nNodes, nElem))
        for i in xrange(nNodes):
            f.write('%f %f %f\n'%(self.X[i, 0], self.X[i, 1], self.X[i, 2]))
        for i in xrange(len(self.conn)):
            f.write('%d %d %d %d\n'%(self.conn[i][0]+1, self.conn[i][1]+1, 
                                     self.conn[i][2]+1, self.conn[i][3]+1))
        f.close()
      
    def _setOptions(self):
        """Internal function to set the options in pyHyp"""
        self.hyp.hypinput.n         = self.options['N']
        self.hyp.hypinput.s0        = self.options['s0']
        self.hyp.hypinput.rmin      = self.options['rMin']
        self.hyp.hypinput.ps0        = self.options['ps0']
        self.hyp.hypinput.pgridratio = self.options['pGridRatio']
        self.hyp.hypinput.slexp     = self.options['slExp']
        self.hyp.hypinput.epse      = self.options['epsE']
        self.hyp.hypinput.epsi      = self.options['epsI']
        self.hyp.hypinput.theta     = self.options['theta']
        self.hyp.hypinput.volcoef   = self.options['volCoef']
        self.hyp.hypinput.volblend  = self.options['volBlend']
        self.hyp.hypinput.cmax      = self.options['cMax']
        self.hyp.hypinput.volsmoothiter = self.options['volSmoothIter']
        self.hyp.hypinput.kspreltol = self.options['kspRelTol']
        self.hyp.hypinput.kspmaxits = self.options['kspMaxIts']
        self.hyp.hypinput.preconlag = self.options['preConLag']
        self.hyp.hypinput.nonlinear = self.options['nonLinear']
        self.hyp.hypinput.kspsubspacesize = self.options['kspSubspaceSize']
        self.hyp.hypinput.writemetrics = self.options['writeMetrics']

    def _checkOptions(self):
        """
        Check the solver options against the default ones
        and add option iff it is NOT in solver_options
        """
        for key in self.options_default.keys():
            if not(key in self.options.keys()):
                self.options[key] = self.options_default[key]

    def __del__(self):
        """
        Clean up fortran allocated values if necessary
        """
        self.hyp.releasememory()

    def _readPlot3D(self, fileName, fileType='ascii', order='f'):
        """Load a plot3D file and create the required unstructured
        data.

        Parameters
        ----------
        fileName : str
            The plot3d filename
    
        fileType : str
            One of 'ascii' or 'binary'. The binary loader is flaky and
            may not work correctly. It is recommended to use ascii format
            if possible.
        order : str
            The internal ordering of the plot3d file. Must be one of
            'f' or 'c'. Fortran ordering ('f') is the usual case. 
        """
        
        if fileType == 'ascii':
            if self.comm.rank == 0:
                print('Loading ascii plot3D file: %s ...'%(fileName))
            binary = False
            f = open(fileName, 'r')
        else:
            if self.comm.rank == 0:
                print('Loading binary plot3D file: %s ...'%(fileName))
            binary = True
            f = open(fileName, 'rb')
       
        if binary:
            itype = self._readNValues(f, 1, 'int', binary)[0]
            nSurf = self._readNValues(f, 1, 'int', binary)[0]
            itype = self._readNValues(f, 1, 'int', binary)[0] # Need these
            itype = self._readNValues(f, 1, 'int', binary)[0] # Need these
            sizes   = self._readNValues(
                f, nSurf*3, 'int', binary).reshape((nSurf, 3))
        else:
            nSurf = self._readNValues(f, 1, 'int', binary)[0]
            sizes   = self._readNValues(
                f, nSurf*3, 'int', binary).reshape((nSurf, 3))

        # ONE of Patch Sizes index must be one
        nPts = 0
        for i in xrange(nSurf):
            if sizes[i, 0] == 1: # Compress back to indices 0 and 1
                sizes[i, 0] = sizes[i, 1]
                sizes[i, 1] = sizes[i, 2] 
            elif sizes[i, 1] == 1:
                sizes[i, 1] = sizes[i, 2]
            elif sizes[i, 2] == 1:
                pass
            else:
                raise Error("Surface %d does not have one block index "
                            "dimension of 1" % i)
            nPts += sizes[i, 0]*sizes[i, 1]

        if self.comm.rank == 0:
            print(' -> nSurf = %d'%(nSurf))
            print(' -> Surface Points: %d'%(nPts))

        surfs = []
        for i in xrange(nSurf):
            cur_size = sizes[i, 0]*sizes[i, 1]
            surfs.append(numpy.zeros([sizes[i, 0], sizes[i, 1], 3]))
            for idim in xrange(3):
                surfs[-1][:, :, idim] = self._readNValues(
                    f, cur_size, 'float', binary).reshape(
                    (sizes[i, 0], sizes[i, 1]), order=order)
        f.close()

        # Final list of patch sizes:
        sizes = []
        for iSurf in xrange(nSurf):
            sizes.append([surfs[iSurf].shape[0], surfs[iSurf].shape[1]])

        # Concatenate the points into a flat array:
        pts = []
        for ii in xrange(nSurf):
            for j in xrange(sizes[ii][1]):
                for i in xrange(sizes[ii][0]):
                    pts.append(surfs[ii][i, j])

        # Find the unique points -- (efficient) spatial search. 
        uniquePts, link, nUnique = self.hyp.pointreducewrap(
            numpy.array(pts).T, self.options['nodeTol'])

        # Just take the actual number of unique points
        uniquePts = uniquePts[:, 0:nUnique].T

        # Convert to 0-based ordering
        link = link - 1

        # We will do a greedy-reordering with substructuring of each
        # of the patch blocks. Basically what we want to do is reorder
        # the nodes such that nodes that are logically next to each
        # other are close in the vector. This will minimize
        # communication during the hyperbolic marching. This is
        # effectively a poor-man's partitioning algorithm.
        MAX_NODE_SIZE = 17
                
        ptInd = -1*numpy.ones(nUnique, dtype='intc')
        newNodes = numpy.zeros((nUnique, 3))
        conn = []
        l_index = []
        pt_offset = 0
        iNode = 0
        for iSurf in range(nSurf):
            l_index.append(numpy.zeros(sizes[iSurf],'intc'))

            iBlocks = sizes[iSurf][0] // MAX_NODE_SIZE + 1
            jBlocks = sizes[iSurf][1] // MAX_NODE_SIZE + 1
       
            # Loop over nodes
            for ii in range(iBlocks):
                for jj in range(jBlocks):
                    iStart = ii*MAX_NODE_SIZE
                    iEnd = min(sizes[iSurf][0], (ii+1)*MAX_NODE_SIZE)

                    jStart = jj*MAX_NODE_SIZE
                    jEnd = min(sizes[iSurf][1], (jj+1)*MAX_NODE_SIZE)

                    for j in range(jStart, jEnd): 
                        for i in range(iStart, iEnd):
                          
                            # This is the original index of this point
                            ind = pt_offset + j*sizes[iSurf][0] + i

                            # Check if we have used this node yet:
                            if ptInd[link[ind]] == -1:
                                ptInd[link[ind]] = iNode
                                newNodes[iNode] = uniquePts[link[ind]]
                                iNode += 1

                            l_index[iSurf][i,j] = ptInd[link[ind]]
                        
            # Loop over faces for the conn array:
            for i in range(sizes[iSurf][0]-1):
                for j in range(sizes[iSurf][1]-1):
                    conn.append([l_index[iSurf][i,   j  ],
                                l_index[iSurf][i+1, j  ],
                                l_index[iSurf][i+1, j+1],
                                l_index[iSurf][i  , j+1]])
                    
            pt_offset += sizes[iSurf][0] * sizes[iSurf][1]
            
        return surfs, sizes, newNodes, conn, l_index

    #-----------------------
    # Some helper routines
    #-----------------------
    def _e_dist(self, x1, x2):
        """Get the eculidean distance between two points"""
        if len(x1) == 3:
            return numpy.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
        elif len(x1) == 2:
            return numpy.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
        elif len(x1) == 1:
            return numpy.abs(x1[0]-x2[0])

    def _reverseRows(self, in_array):
        """Flip Rows (horizontally)"""
        rows = in_array.shape[0]
        cols = in_array.shape[1]
        output = numpy.empty([rows, cols], in_array.dtype)
        for row in xrange(rows):
            output[row] = in_array[row][::-1].copy()

        return output

    def _readNValues(self, handle, N, dtype, binary=False):
        """Read 'N' values of type 'float' or 'int' from file handle
        'handle'"""
        if binary:
            sep = ""
        else:
            sep = ' '

        if dtype == 'int':
            values = numpy.fromfile(handle, dtype='int32', count=N, sep=sep)
        else:
            values = numpy.fromfile(handle, dtype='float', count=N, sep=sep)
        return values

    def _determineMGLevels(self, sizes):
        """Determine the maximum number of 'multi grid' levels
        possible for the given set of sizes. We print some warnings if
        we can't do *any* MG based on the size

        Parameters
        ----------
        sizes : list
            List of lists. Each entry contains the u-v size of the patch. This is the
            NODAL size. 

        Returns
        -------
        levelMax : int
            The maximum number of multigrid levels we can do
        """

        # First get the total size for reference
        nPt = 0
        for size in sizes:
            nPt += size[0]*size[1]

        levelMax = 1000 # Arbitrary large number 
        tmp = copy.deepcopy(sizes) # Copy since we will overwrite
                                   # sizes
        for size in tmp:
            for j in range(levelMax):
                if not numpy.mod(size[0]-1, 2) == 0:
                    levelMax = j + 1
                    break
                if not numpy.mod(size[1]-1, 2) == 0:
                    levelMax = j + 1
                    break
                size[0] = (size[0]-1)/2 + 1
                size[1] = (size[1]-1)/2 + 1

        if nPt > 2000 and levelMax == 1:
            print("Warning: There are more than 2000 panels but the "
                  "node counts do not permit multigrid. This will "
                  "substantially slow down elliptic generation. "
                  "Regridding with a multigrid friendly number of "
                  "nodes is *highly* recommended")

        return levelMax


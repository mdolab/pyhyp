#!/usr/bin/python
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
# Standard Python modules
# =============================================================================
import sys, os, time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================

# =============================================================================
# Global MPI Print Function
# =============================================================================
try:
    from mpi4py import MPI
    USE_MPI = True
except:
    USE_MPI = False

def mpiPrint(print_string, NO_PRINT=False, comm=None):
    if not NO_PRINT:
        if USE_MPI:
            if comm == None:
                comm = MPI.COMM_WORLD
            # end if
            if comm.rank == 0:
                print print_string
            # end if
        # end if
        else:
            print print_string
        # end if
    # end if

    return

# =============================================================================
# pyHyp class
# =============================================================================

class pyHyp(object):
    """
    This is the main class pyHyp. It is used as the user interface to pyHyp. 
    """

    def __init__(self, dimension, fileName=None, X=None, comm=None, 
                 options=None, flip=False, **kwargs):
        """
        Create the pyHyp object. 

        Required Input Argument:
            dimension, str: One of '2d' or '3d'. Specify the dimension 
                            of the problem

        Input Arguments: ONE of:
            fileName, str: file containing xy points 
            X, numpy array, size (N, 2): Numpy array containing N x-y points

        Optional Arguments:
            comm, MPI_INTRA_COMM: MPI communication (as obtained 
                from mpi4py) on which to create the pyHyp object. This 
                allows faster computations in parallel. 
                If not provided, MPI_COMM_WORLD is used. 
            Options, dictionary: A dictionary containing the 
                 the options for hyperbolic mesh generator.

         Returns:
             None

             """

        # Defalut options for hyperbolc generation
        self.options_default = {
            # ---------------------------
            #        Grid Parameters
            # ---------------------------
            # Number of layers:
            'N': 100, 

            # Initial off-wall spacing
            's0':0.01,
            
            # Rmin: Distance to march in multiples of initial radius
            'rMin': 50,
            
            # Maximum permissible ratio of marching direction length
            # to smallest in-plane edge
            'cMax':6.0,

            #nonLinear: True/False. Use nonlinear scheme
            'nonLinear': False,

            #slExp: Exponent for Sl calc. Don't change this value
            #unless you know what you are doing!
            'slExp': .15,

            # ---------------------------
            #   Pseudo Grid Parameters
            # ---------------------------
            # Initial off-wall spacing
            'ps0':0.01,
            
            # Grid Spacing Ratio
            'pGridRatio':1.15,

            # ---------------------------
            #   Smoothing parameters
            # ---------------------------

            # epsE: The explict smoothing coefficient
            'epsE': 0.5,

            # epsI: The implicit smoothing coefficient
            'epsI': 1.5,

            # theta: The barth implicit smoothing coefficient
            'theta': 1.0,

            # volCoef: The volume smoothing coefficinet for
            # pointJacobi iterations
            'volCoef': .16,

            # volBlend: The volume blending coefficient to force
            # uniform sizs in farfield
            'volBlend': 0.0001,

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 10,

            # ---------------------------
            #   Solution Parameters
            # ---------------------------
            # kspRelTol: Solution tolerance for linear system
            'kspRelTol': 1e-8,
            
            # Maximum number of iterations to run for linear system
            'kspMaxIts': 500,

            # Preconditioner Lag
            'preConLag': 10,

            'kspSubspaceSize':50,

            # ---------------------------
            #   Output Parameters
            # ---------------------------
            'writeMetrics': False,
            }

        # Import and set the hyp module
        import hyp
        self.hyp = hyp

        # Set the possible MPI Intracomm
        if comm is None:
            try:
                self.comm = MPI.COMM_WORLD
            except:
                self.comm = None
        else:
            self.comm = comm
        # end if

        # Initialize PETSc if not already done so:
        self.hyp.initpetsc()

        # Set the dimension:
        if not dimension.lower() in ['2d', '3d']:
            mpiPrint('Error: \'dimension\' must be one of \'2d\' or \'3d\'', 
                     comm=self.comm)
        else:
            self.twoD = False
            if dimension.lower() == '2d':
                self.twoD = True
            # end if
        # end if

        # Use supplied options
        if options is not None:
            self.options = options
        else:
            self.options = {}
        # end if

        # No mirroring by default
        self.xMirror = False
        self.yMirror = False
        self.zMirror = False

        # Depending on the dimensionality, we initialize the problem
        # separately. 
        if self.twoD:
            self._init2d(fileName, X, flip)
        else:
            self._init3d(fileName, **kwargs)
        # end if

        # Setup the options
        self._checkOptions()
        self._setOptions()
        self.gridGenerated = False
        mpiPrint('Problem sucessfully setup!')

        return

    def _init2d(self, fileName, X, flip):

        # Check to see how the user has passed in the data:
        if fileName is None and X is None:
            mpiPrint('Error: Either fileName or X must be passed on \
initialization!')
            return

        if fileName is not None and X is not None:
            mpiPrint('Error: BOTH fileName or X have been passed on \
initialization!')
            mpiPrint('Error: Only ONE of these are required.')
            return            

        # Take the file name and try loading it using the pyGeo
        # architecture. Only do this on one proc; we will MPI-bcast
        # the information back to everyone else when it is necessary.
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
                # end for

                # Done with the file so close
                f.close()

                # Convert the lists to an array
                X = numpy.array([x,y]).T
                
            # end if

            # Now check if the first and the last point are within
            # tolerance of each other:
            if self.e_dist(X[0],X[-1]) < 1e-6:
                # Curve is closed, we're ok. Internally, we work
                # with the unclosed curve and explictly deal with
                # the  0th and Nth points 
                X = X[0:-1,:]
            else:
                # Curve is not closed, print warning
                print '*'*105
                print ' Warning: Curve was not closed! Closing curve with a \
linear segment. This may or not be what is desired!'
                print '*'*105
            # end if

            if flip: # Reverse direction
                X[:,0] = X[:,0][::-1]
                X[:,1] = X[:,1][::-1]
        # end if

        # Now broadcast the required info to the other procs:
        self.X = self.comm.bcast(X, root=0)
        self.N = len(self.X)

        return

    def _init3d(self, fileName, **kwargs):

        if 'xMirror' in kwargs and kwargs['xMirror']:
            self.xMirror=True
        if 'yMirror' in kwargs and kwargs['yMirror']:
            self.yMirror=True
        if 'zMirror' in kwargs and kwargs['zMirror']:
            self.zMirror=True

        delFile = False
        if self.xMirror or self.yMirror or self.zMirror:

            # We need to first read the file, mirror it, and write it
            # back to a tmp file such that the nominal read below
            # works.

            # Load in the file and determine the free edges, these
            # will be symmetry edges (which turn into symmetry planes)
            surfs, sizes, nodes, conn, l_index = \
                self._readPlot3D(fileName, **kwargs)
            nSurf = len(surfs)

            # To determine the nodes on the symmetry plane, we need to
            # figure out which edges have a face only on one
            # side. Finally, the end goal is to prodce an
            # (unstructured) list that looks like: 
            #
            # [[surfID_1, i_1, j_1], [surfID_2, i_2, j_2] ...]
            # 
            # where the i,j gives the index on surf that is initialy
            # on the mirror plane. This information can then be used
            # to zero this "line" to be precisely on the mirror plane
            # as a post processing step. This method does not depend
            # on the 
            
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
                # end if

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
                    # end for
                # end for
            # end for

            # Now generate a final list of nodes that are on the mirror plane:
            mirrorNodes = []
            for key in edges:
                if edges[key] == 1:
                    mirrorNodes.append(key[0])
                    mirrorNodes.append(key[1])
                # end if
            # end for

            # Uniquify them
            uniqueMirrorNodes = numpy.unique(mirrorNodes)

            # Now loop back over the boundries of the surfaces. For
            # each node on the boundary check if it is in uniqueMirror nodes. 

            self.symNodes = []
            for ii in xrange(len(surfs)):
                imax = sizes[ii][0]
                jmax = sizes[ii][1]
                for i in xrange(imax):
                    if l_index[ii][i, 0] in uniqueMirrorNodes:
                        self.symNodes.append([ii, i  , 0])
                    if l_index[ii][i, jmax-1] in uniqueMirrorNodes:
                        self.symNodes.append([ii, i,  jmax-1])
                    # end if

                for j in xrange(jmax):
                    if l_index[ii][0, j] in uniqueMirrorNodes:
                        self.symNodes.append([ii, 0, j  ])
                    if l_index[ii][imax-1, j] in uniqueMirrorNodes:
                        self.symNodes.append([ii, imax-1, j  ])
                    # end if
                # end for
            # end for

            # Now that we know the symnodes, we can hard zero the
            # nodes on the mirror plane. Note that if your surface
            # isn't closed, some really really stuff will result since
            # those nodes will be zeroed!
            for iSym in xrange(len(self.symNodes)):
                iSurf = self.symNodes[iSym][0]
                i     = self.symNodes[iSym][1]
                j     = self.symNodes[iSym][2]
                if self.xMirror:
                    surfs[iSurf][i,j,0] = 0.0
                if self.yMirror:
                    surfs[iSurf][i,j,1] = 0.0
                if self.zMirror:
                    surfs[iSurf][i,j,2] = 0.0
                # end if
            # end for
            
            # Conver symnodes to array and index everything my 1 
            self.symNodes = numpy.array(self.symNodes) + 1
            
            # Note that we are guaranteed that symNodes is sorted by
            # patchID. This is important since it means we can loop
            # over the blocked in sequence and all nodes need to 

            # Generate new list of sizes (double the length)
            newSizes = numpy.zeros((nSurf*2, 2),'intc')
            for i in xrange(nSurf):
                newSizes[i] = sizes[i]
                newSizes[i+nSurf] = sizes[i]
            
            # Now mirror each zone, while flipping the i and j index
            for i in xrange(nSurf):
                surfs.append(numpy.zeros([newSizes[i+nSurf,0], newSizes[i+nSurf,1],3]))
                if self.xMirror:
                    surfs[i+nSurf][:,:,0] = -self.reverseRows(surfs[i][:,:,0])
                else:
                    surfs[i+nSurf][:,:,0] = self.reverseRows(surfs[i][:,:,0])
                # end if

                if self.yMirror:
                    surfs[i+nSurf][:,:,1] = -self.reverseRows(surfs[i][:,:,1])
                else:
                    surfs[i+nSurf][:,:,1] = self.reverseRows(surfs[i][:,:,1])
                # end if
                if self.zMirror:
                    surfs[i+nSurf][:,:,2] = -self.reverseRows(surfs[i][:,:,2])
                else:
                    surfs[i+nSurf][:,:,2] = self.reverseRows(surfs[i][:,:,2])
                # end if
            # end for

            # Dump back out
            f = open('tmp.fmt','w')
            f.write('%d\n'%(nSurf*2))
            for i in xrange(nSurf*2):
                f.write('%d %d %d\n'%(newSizes[i][0], newSizes[i][1], 1))
            for ii in xrange(nSurf*2):
                for idim in xrange(3):
                    for j in xrange(newSizes[ii][1]):
                        for i in xrange(newSizes[ii][0]):
                            f.write('%20.13g\n'%(surfs[ii][i,j,idim]))
                        # end for
                    # end for
                # end for
            # end for
            f.close()
            fileName = 'tmp.fmt'
            delFile=True
        else:
            self.symNodes = None
        # end if

        surfs, sizes, nodes, conn, l_index = self._readPlot3D(fileName, **kwargs)
        nSurf = len(surfs)
        self.conn = conn
        self.X = nodes
        self.sizes = sizes
        self.l_index = l_index
        nGlobal = len(nodes)

        if delFile:
            os.remove(fileName)

        # We really should check if the surface is closed. This check
        # is not implemented yet -- so user beware!
    
        # Next we need to produced an unstructured-like data
        # structure

        # Now we can invert and get the elements surrounding the nodes:
        nodeToElem = [[] for i in xrange(nGlobal)]
        for iElem in xrange(len(self.conn)):
            for ii in xrange(4): # 4 nodes on each face
                # Append the elem Number and the node number of the node on theface
                nodeToElem[self.conn[iElem][ii]].append([iElem, ii]) 
            # end for
        # end for

        # Finally we can get the final node pointer structure we need:
        nPtr = [[] for i in xrange(nGlobal)]
        for iNode in xrange(nGlobal):
            nFaceNeighbours = len(nodeToElem[iNode])

            # Regular nodes get 4 neightbours, others get double
            if nFaceNeighbours == 2:
                print 'Can\'t do geometries with nodes of valance 2!'
                sys.exit(0)
                nPtr[iNode].append(nFaceNeighbours*2)
            elif nFaceNeighbours == 3:
                nPtr[iNode].append(nFaceNeighbours*2)
            else:
                nPtr[iNode].append(nFaceNeighbours)
            # end if

            # Get the face where this node is node 0
            firstID  = nodeToElem[iNode][0][0]
            nodeID  = numpy.mod(nodeToElem[iNode][0][1]+1,4)
            nextElemID = None

            while nextElemID <> firstID:
                if nextElemID is None:
                    nextElemID = firstID
                # end if

                # Append the next node along that edge:
                nodeToAdd = self.conn[nextElemID][numpy.mod(nodeID,4)]
                nPtr[iNode].append(nodeToAdd)

                # Add the diagonal if not regular node
                if nFaceNeighbours == 3:
                    nPtr[iNode].append(self.conn[nextElemID][numpy.mod(nodeID+1,4)])
                # end if
                
                # And we need to find the face that contains the following node:
                nodeToFind = self.conn[nextElemID][numpy.mod(nodeID + 2, 4)]

                found = False
                jj = -1
                while not found:
                    jj = jj + 1
                    # Find the next face
                    faceToCheck = nodeToElem[iNode][jj][0]
                    if faceToCheck != nextElemID:
                        # Check the 4 nodes on this face 
                        for ii in xrange(4):
                            if self.conn[faceToCheck][ii] == nodeToFind:
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
            if len(nPtr[i]) > lMax:
                lMax = len(nPtr[i])
            # end if
        # end for

        nPtrArray = numpy.zeros((len(nPtr), lMax), 'intc')
        for i in xrange(len(nPtr)):
            l = len(nPtr[i])
            nPtrArray[i, 0] = nPtr[i][0]
            nPtrArray[i, 1:l] = numpy.array(nPtr[i][1:]) + 1
        # # end for
      
        # Now we have all the data we need so we can go ahead and
        # initialze the 3D generation

        self.hyp.init3d(sizes, nPtrArray.T)
        for i in xrange(nSurf):
            self.hyp.setlindex(l_index[i], i)
        # end for
        return

    def run(self):
        """
        Run given using the options given
        """
        if self.twoD:
            self.hyp.run2d(self.X.T)
        else:
            self.hyp.run3d(self.X.T)
        # end if
        self.gridGenerated = True

        return 

    def writePlot3D(self, fileName):
        """After we have generated a grid, write it out to a plot3d
        file for the user to look at"""
        
        if self.gridGenerated:
            if self.twoD:
                self.hyp.writeplot3d_2d(fileName)
            else:
                self.hyp.writeplot3d_3d(fileName)
            # end if
        else:
            mpiPrint('Error! No grid has been generated! Run the run() \
command before trying to write the grid!')
        # end if

        return

    def writeCGNS(self, fileName):
        """After we have generated a grid, write it out in a properly \
formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if self.gridGenerated:
            if self.twoD:
                self.hyp.writecgns_2d(fileName)
            else:
                self.hyp.writecgns_3d(fileName)
                # Possibly zero mirror plane for mirrored geometries
                if self.symNodes is not None:
                    if self.xMirror:
                        mirrorDim = 1
                    elif self.yMirror:
                        mirrorDim = 2
                    else:
                        mirrorDim = 3
                    # end if
                    self.hyp.zeromirrorplane(fileName, self.symNodes, mirrorDim)
                # end if
            # end if
        else:
            mpiPrint('Error! No grid has been generated! Run the run() \
command before trying to write the grid!')
        # end if

        return

    def writeCGNSOrig(self, fileName):
        """After we have generated a grid, write it out in a properly \
formatted 1-Cell wide CGNS file suitable for running in SUmb."""

        if self.gridGenerated and not self.twoD:
            self.hyp.writecgns_3dorig(fileName)
        else:
            mpiPrint('Error! No grid has been generated or 2D! Run the run() \
command before trying to write the grid!')
        # end if

        return
    def writeFEMesh(self, fileName):
        """ Ouput a tecplot FE mesh of the surface. Useful for debugging numberings"""

        nNodes = len(self.X)
        nElem = len(self.conn)
        f = open(fileName, 'w')
        f.write("FE Data\n")
        f.write('VARIABLES = \"X\", \"Y\", \"z\" \n')
        f.write("ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n"%(nNodes, nElem))
        for i in xrange(nNodes):
            # Compute the normal
            f.write('%f %f %f\n'%(self.X[i,0], self.X[i,1], self.X[i,2]))
        # end for
        for i in xrange(len(self.conn)):
            f.write('%d %d %d %d\n'%(self.conn[i][0]+1, self.conn[i][1]+1, 
                                     self.conn[i][2]+1, self.conn[i][3]+1))
        # end for
        
        # Close file
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
        # Mirror is a special case. 
        if self.xMirror or self.yMirror or self.zMirror:
            self.hyp.hypinput.writemirror = False
        else:
            self.hyp.hypinput.writemirror = True
        # end if

        return

    def _checkOptions(self):
        """
        Check the solver options against the default ones
        and add option iff it is NOT in solver_options
        """
        for key in self.options_default.keys():
            if not(key in self.options.keys()):
                self.options[key] = self.options_default[key]
            # end if
        # end for

        return

    def __del__(self):
        """
        Clean up fortran allocated values if necessary
        """
        self.hyp.releasememory()

        return

    def _readPlot3D(self, file_name, file_type='ascii', order='f', **kwargs):
        '''Load a plot3D file and create the required unstructured data.

        file_name: Required string for file
        file_type: 'ascii' or 'binary'
        order: 'f' for fortran ordering (usual), 'c' for c ordering
        '''
        if file_type == 'ascii':
            mpiPrint('Loading ascii plot3D file: %s ...'%(file_name))
            binary = False
            f = open(file_name, 'r')
        else:
            mpiPrint('Loading binary plot3D file: %s ...'%(file_name))
            binary = True
            f = open(file_name, 'rb')
        # end if
        if binary:
            itype = self.readNValues(f, 1, 'int', binary)[0]
            nSurf = self.readNValues(f, 1, 'int', binary)[0]
            itype = self.readNValues(f, 1, 'int', binary)[0] # Need these
            itype = self.readNValues(f, 1, 'int', binary)[0] # Need these
            sizes   = self.readNValues(
                f, nSurf*3, 'int', binary).reshape((nSurf, 3))
        else:
            nSurf = self.readNValues(f, 1, 'int', binary)[0]
            sizes   = self.readNValues(
                f, nSurf*3, 'int', binary).reshape((nSurf, 3))
        # end if

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
                mpiPrint('Error: One of the plot3d indices must be 1')
            # end if
            nPts += sizes[i, 0]*sizes[i, 1]
        # end for
        mpiPrint(' -> nSurf = %d'%(nSurf))
        mpiPrint(' -> Surface Points: %d'%(nPts))

        surfs = []
        for i in xrange(nSurf):
            cur_size = sizes[i, 0]*sizes[i, 1]
            surfs.append(numpy.zeros([sizes[i, 0], sizes[i, 1], 3]))
            for idim in xrange(3):
                surfs[-1][:, :, idim] = self.readNValues(
                    f, cur_size, 'float', binary).reshape(
                    (sizes[i, 0], sizes[i, 1]), order=order)
            # end for
        # end for

        f.close()

        # Final list of patch sizes:
        sizes = []
        for iSurf in xrange(nSurf):
            sizes.append([surfs[iSurf].shape[0], surfs[iSurf].shape[1]])
        # end for

        # Concatenate the points into a flat array:
        pts = []
        
        for ii in xrange(nSurf):
            for j in xrange(sizes[ii][1]):
                for i in xrange(sizes[ii][0]):
                    pts.append(surfs[ii][i,j])
        # end for

        # Find the unique points
        uniquePts, link, nUnique = self.hyp.pointreducewrap(numpy.array(pts).T, 1e-8)

        # Just take the actual number of unique points
        uniquePts = uniquePts[:,0:nUnique].T

        # Convert to 0-based ordering
        link = link - 1

        # Now we can easily compute the connectivity and the l_index
        conn = []
        l_index = []
        pt_offset = 0
        for ii in xrange(nSurf):
            l_index.append(numpy.zeros(sizes[ii],'intc'))

            # Loop over FACES, take nodes around face and index into link to get
            # the index of the original node (from surfs) on the
            # compacted node list
            for j in xrange(sizes[ii][1]-1):
                for i in xrange(sizes[ii][0]-1):
                    conn.append([link[pt_offset + (j  )*sizes[ii][0] + i  ],
                                 link[pt_offset + (j  )*sizes[ii][0] + i+1],
                                 link[pt_offset + (j+1)*sizes[ii][0] + i+1],
                                 link[pt_offset + (j+1)*sizes[ii][0] + i  ]])
                # end for
            # end for

            # Loop over the NODES, and index directly into link
            for i in xrange(sizes[ii][0]):
                for j in xrange(sizes[ii][1]):
                    l_index[ii][i, j] = link[pt_offset + j*sizes[ii][0] + i]
                    # end for
            # end for

            # Increment pt_offset to account for the number of nodes
            # we've gone through on this patch
            pt_offset += sizes[ii][0] * sizes[ii][1]
        # end if (surf loop)

        return surfs, sizes, uniquePts, conn, l_index

    # Some helper routines
    def e_dist(self, x1, x2):
        '''Get the eculidean distance between two points'''
        if len(x1) == 3:
            return numpy.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
        elif len(x1) == 2:
            return numpy.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
        elif len(x1) == 1:
            return numpy.abs(x1[0]-x2[0])

    def reverseRows(self, in_array):
        '''Flip Rows (horizontally)'''
        rows = in_array.shape[0]
        cols = in_array.shape[1]
        output = numpy.empty([rows, cols], in_array.dtype)
        for row in xrange(rows):
            output[row] = in_array[row][::-1].copy()
        # end for

        return output

    def readNValues(self, handle, N, dtype, binary=False, sep=' '):
        '''Read N values of dtype 'float' or 'int' from file handle'''
        if binary:
            sep = ""
        # end if

        if dtype == 'int':
            values = numpy.fromfile(handle, dtype='int32', count=N, sep=sep)
        else:
            values = numpy.fromfile(handle, dtype='float', count=N, sep=sep)
        return values

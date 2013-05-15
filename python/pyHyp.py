#!/usr/bin/python
from __future__ import division
"""
pyHyp2D

The pyHyp2D module is used to generate a hyperbolic grid around a
closed curve

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
from mdo_import_helper import import_modules, MPI, mpiPrint
exec(import_modules('geo_utils','pyGeo'))

# =============================================================================
# MultiBlockMesh class
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
            'pN': 1000,

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
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm
        # end if

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
            if geo_utils.e_dist(X[0],X[-1]) < 1e-6:
                # Curve is closed, we're ok. Internally, we work
                # with the unclosed curve and explictly deal with
                # the  0th and Nth points 
                X = X[0:-1,:]
            else:
                # Curve is not closed, print warning
                print '*'*80
                print 'Warning: Curve was not closed! Closing curve with a \
linear segment. This may or not be what is desired!'
                print '*'*80
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

        self.xMirror = False
        self.yMirror = False
        self.zMirror = False
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
            geo_obj = pyGeo.pyGeo('plot3d', file_name=fileName,
                                  file_type='ascii', order='f')

            nodeTol = kwargs.pop('nodeTol',1e-6)
            edgeTol = kwargs.pop('edgeTol',1e-6)
            geo_obj._calcConnectivity(nodeTol, edgeTol)

            sizes = []
            nSurf = geo_obj.nSurf
            for isurf in xrange(nSurf):
                sizes.append([geo_obj.surfs[isurf].Nu, geo_obj.surfs[isurf].Nv])
            geo_obj.topo.calcGlobalNumbering(sizes)
            topo = geo_obj.topo

            # Determine the edges that do not have a neighbour
            edgeCount = numpy.zeros(topo.nEdge, 'intc')
            for iFace in xrange(topo.nFace):
                for iEdge in xrange(4):
                    edgeCount[topo.edge_link[iFace][iEdge]] += 1
                # end for
            # end for
                    
            # Now make a list of the patchID, and edge number for
            # later use
            self.symEdges = [] 
            for iFace in xrange(topo.nFace):
                for iEdge in xrange(4):
                    if edgeCount[topo.edge_link[iFace][iEdge]] == 1:
                        self.symEdges.append([iFace,iEdge])
                    # end if
                # end for
            # end for

            surfs = []
            for i in xrange(nSurf):
                surfs.append(geo_obj.surfs[i].X)

            # Generate new list of sizes (double the length)
            newSizes = numpy.zeros((nSurf*2, 2),'intc')
            for i in xrange(nSurf):
                newSizes[i] = sizes[i]
                newSizes[i+nSurf] = sizes[i]
            
            # Now mirror each zone, while flipping the i and j index
            for i in xrange(nSurf):
                surfs.append(numpy.zeros([newSizes[i+nSurf,0], newSizes[i+nSurf,1],3]))
                if self.xMirror:
                    surfs[i+nSurf][:,:,0] = -geo_utils.reverseRows(surfs[i][:,:,0])
                else:
                    surfs[i+nSurf][:,:,0] = geo_utils.reverseRows(surfs[i][:,:,0])
                # end if

                if self.yMirror:
                    surfs[i+nSurf][:,:,1] = -geo_utils.reverseRows(surfs[i][:,:,1])
                else:
                    surfs[i+nSurf][:,:,1] = geo_utils.reverseRows(surfs[i][:,:,1])
                # end if
                if self.zMirror:
                    surfs[i+nSurf][:,:,2] = -geo_utils.reverseRows(surfs[i][:,:,2])
                else:
                    surfs[i+nSurf][:,:,2] = geo_utils.reverseRows(surfs[i][:,:,2])
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
            self.symEdges = None
        # end if

        geo_obj = pyGeo.pyGeo('plot3d', file_name=fileName,
                              file_type='ascii', order='f')
        if delFile:
            os.remove(fileName)
        # This isn't prety or efficient but does produce the
        # globally reduced list we need and it is greedily
        # reordered. For ~100,000 surface points this takes on the
        # order of a couple of seconds so we're not going to worry
        # too much about efficiency here.
        nodeTol = kwargs.pop('nodeTol',1e-6)
        edgeTol = kwargs.pop('edgeTol',1e-6)
        geo_obj._calcConnectivity(nodeTol, edgeTol)
        sizes = []
        for isurf in xrange(geo_obj.nSurf):
            sizes.append([geo_obj.surfs[isurf].Nu, geo_obj.surfs[isurf].Nv])
        geo_obj.topo.calcGlobalNumbering(sizes)
        topo = geo_obj.topo
        # Now we need to check that the surface is closed, In the
        # future it will support symmetry boundary conditions, but not
        # yet. For a closed topology, every edge must be exactly
        # shared:
    
        edgeCount = numpy.zeros(topo.nEdge, 'intc')
        for iFace in xrange(topo.nFace):
            for iEdge in xrange(4):
                edgeCount[topo.edge_link[iFace][iEdge]] += 1
            # end for
        # end for

        if numpy.all(edgeCount == 2): 
            print 'We have a closed surface so we\'re good to go!'
        else:
            print 'Surface is not closed. Can\'t do this yet'
            sys.exit(0)
        # end if

        # Next we need to produced an unstructured-like data
        # structure

        nGlobal = topo.nGlobal
        # Connectivey of elements:  
        self.conn = []
        for iFace in xrange(topo.nFace):
            for j in xrange(topo.l_index[iFace].shape[1]-1):
                for i in xrange(topo.l_index[iFace].shape[0]-1):
                    self.conn.append([topo.l_index[iFace][i  ,j  ],
                                      topo.l_index[iFace][i+1,j  ],
                                      topo.l_index[iFace][i+1,j+1],
                                      topo.l_index[iFace][i  ,j+1]])
                # end for
            # end for
        # end for

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

        # # Special treatment for all nodes ATTACHED to extraordinary
        # # nodes. All these nodes themselves shoudl already have 4
        # # neighbours, (one of which of course is the extraordinary
        # # node). So we have to flag them with nneighbour = -4 to allow
        # # the fortran code to detect these nodes
     
        # for i in xrange(len(nPtr)):
        #     # Check if this node is extraordinary
        #     if nPtr[i][0] > 6 and nPtr[i][0] > 0:
        #         # Flag each of these nodes
        #         for ii in xrange(len(nPtr[i][1:])):
        #             neighbourNode = nPtr[i][ii+1]
        #             # Check if this node has already been flagged:
        #             if nPtr[neighbourNode][0] < 0:
        #                 print 'Somehow this node was already flagged. This is probably  bad news'
        #             nPtr[neighbourNode][0] *= -1
                    
        #             # Also flag the extraordinary node attached to
        #             # neighbourNode so that we know which one to do
        #             # something special with in the fortran code
             

        #             for jj in xrange(-nPtr[neighbourNode][0]):
        #                 if nPtr[neighbourNode][jj+1] == i:
        #                     nPtr[neighbourNode][jj+1] *= -1
        #                 # end if
        #             # end for

        #         # end for
        #     # end if
        # # end for

        # Next, assemble the global X vector:
        self.X = numpy.zeros((topo.nGlobal, 3))
        for iFace in xrange(topo.nFace):
            for j in xrange(topo.l_index[iFace].shape[1]):
                for i in xrange(topo.l_index[iFace].shape[0]):
                    self.X[topo.l_index[iFace][i,j]] = geo_obj.surfs[iFace].X[i,j]
                # end for
            # end for
        # end for

        # Save topology for future reference (writing out FE mesh)
        self.topo = topo

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
        for i in xrange(topo.nFace):
            self.hyp.setlindex(topo.l_index[i], i)
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
                if self.symEdges:
                    if self.xMirror:
                        mirrorDim = 1
                    elif self.yMirror:
                        mirrorDim = 2
                    else:
                        mirrorDim = 3
                    # end if
                    self.hyp.zeromirrorplane(
                        fileName, numpy.array(self.symEdges), mirrorDim)
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

        nNodes = self.topo.nGlobal
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

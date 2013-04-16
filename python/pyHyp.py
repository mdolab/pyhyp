
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
            # Number of layers:
            'N': 10, 

            # Initial off-wall spacing
            's0':0.01,
            
            # Grid Spacing Ratio
            'gridRatio':1.15,

            # epsE: The explict smoothing coefficient
            'epsE': 0.5,

            # epsI: The implicit smoothing coefficient
            'epsI': 1.0,

            # theta: The barth implicit smoothing coefficient
            'theta': 1.0,

            # volCoef: The volume smoothing coefficinet for
            # pointJacobi iterations
            'volCoef': 1.0,

            # volBlend: The volume blending coefficient to force
            # uniform sizs in farfield
            'volBlend': 0.2,

            # volSmoothIter: The number of point-jacobi volume
            # smoothing iterations
            'volSmoothIter': 10,
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
            self._init3d(fileName, flip, **kwargs)
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

    def _init3d(self, fileName, flip, **kwargs):

        geo_obj = pyGeo.pyGeo('plot3d', file_name=fileName,
                              file_type='ascii', order='f')
        
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
    

        # Next we need to produced an unstructured-like data
        # structure that points to the neighbours for each
        # node. We will use node-halo information such that
        # regular corners (has 2 neighbours in each of the i and j
        # directions) will be treated exactly the same as an
        # interior node. Extraordinary node, those that have 3, 5
        # or more verticies will be flagged and will be treated
        # specially in the hyperbolic marchign algorithm using a
        # lapalacian operator. 
        
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

        nGlobal = topo.nGlobal
      
        # Set all nodes as averaging and we just won't add corners
        nPtr = [[] for i in xrange(nGlobal)]
        # Loop over each face and do all points except the
        # corners. For the ones on edges we will use the connectivity
        # to to get the halos. 
  
        edge_mapping = [[] for i in xrange(topo.nEdge)]
        for iFace in xrange(topo.nFace):
            for iEdge in xrange(4):
                uEdge = topo.edge_link[iFace][iEdge]
                if iEdge == 0:
                    vals = topo.l_index[iFace][:, 1].copy()
                elif iEdge == 1:
                    vals = topo.l_index[iFace][:, -2].copy()
                elif iEdge == 2:
                    vals = topo.l_index[iFace][1, :].copy()
                elif iEdge == 3:
                    vals = topo.l_index[iFace][-2, :].copy()
                # end if

                # Flip Direction if necessary
                if topo.edge_dir[iFace][iEdge] == -1: 
                    vals = vals[::-1].copy()
                # end if

                edge_mapping[uEdge].append([iFace, iEdge, vals])
            # end for
        # end for

        for iFace in xrange(topo.nFace):
            
            # Do the interior of each of the 4 edges. We *know* these
            # MUST have halos since the surface is closed
            
            # Edge 0 and 1
            for iEdge in [0,1]:
                j = 0
                if iEdge == 1: 
                    j = sizes[iFace][1]-1
                # end if
                    
                uEdge = topo.edge_link[iFace][iEdge]
                for jj in xrange(2):
                    if edge_mapping[uEdge][jj][0:2] != [iFace, iEdge]:
                        # This is the "other" edge
                        vals = edge_mapping[uEdge][jj][2].copy()
                        
                        # Flip vals if this edge is flipped
                        if topo.edge_dir[iFace][iEdge] == -1:
                            vals = vals[::-1].copy()
                        #end if
                    # end if
                # end for
                for i in xrange(1, sizes[iFace][0] -1):
                    iGlobal = topo.l_index[iFace][i,j]
                    if nPtr[iGlobal] == []: # Not set yet
                        nPtr[iGlobal].extend([0,4])
                        if iEdge == 0: # Right-> Up -> Left -> Down
                            nPtr[iGlobal].append(topo.l_index[iFace][i+1, j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j+1] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i-1, j] + 1)
                            nPtr[iGlobal].append(vals[i] + 1)
                        else:
                            nPtr[iGlobal].append(topo.l_index[iFace][i+1, j] + 1)
                            nPtr[iGlobal].append(vals[i] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i-1, j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j-1] +1)
                            # end if
                    # end if
                # end for
            # end for

            # Edge 2 and 3
            for iEdge in [2,3]:
                i = 0
                if iEdge == 3: 
                    i = sizes[iFace][0]-1
                # end if
                    
                uEdge = topo.edge_link[iFace][iEdge]
                for jj in xrange(2):
                    if edge_mapping[uEdge][jj][0:2] != [iFace, iEdge]:
                        # This is the "other" edge
                        vals = edge_mapping[uEdge][jj][2]
                        
                        # Flip vals if this edge is flipped
                        if topo.edge_dir[iFace][iEdge] == -1:
                            vals = vals[::-1]
                        # end if
                    # end if
                # end for

                for j in xrange(1, sizes[iFace][1] -1):
                    iGlobal = topo.l_index[iFace][i,j]
                    if nPtr[iGlobal]  == []: # Not set yet
                        nPtr[iGlobal].extend([0,4])
                        if iEdge == 2:
                            nPtr[iGlobal].append(topo.l_index[iFace][i+1, j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j+1] + 1)
                            nPtr[iGlobal].append(vals[j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j-1] + 1)
                        else:
                            nPtr[iGlobal].append(vals[j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j+1] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i-1, j] + 1)
                            nPtr[iGlobal].append(topo.l_index[iFace][i, j-1] + 1)
                        # end if
                    # end if
                # end for
            # end for

            # Now do all the interiors. The *always* have the proper neighbours
            for i in xrange(1, sizes[iFace][0]-1):
                for j in xrange(1, sizes[iFace][1]-1):
                    # No need to check if node is added, cannot be
                    # since interior face nodes are unique
                    iGlobal = topo.l_index[iFace][i,j]
                    nPtr[iGlobal].extend([0,4])
                    nPtr[iGlobal].append(topo.l_index[iFace][i+1, j  ] + 1)
                    nPtr[iGlobal].append(topo.l_index[iFace][i  , j+1] + 1)
                    nPtr[iGlobal].append(topo.l_index[iFace][i-1, j  ] + 1)
                    nPtr[iGlobal].append(topo.l_index[iFace][i  , j-1] + 1)
                # end for (Interior J loop)
            # end for (Interior I loop)
        # end for (Face Loop)

        # Everything but the corners are now done.

        # Next, assemble the global X vector:
        self.X = numpy.zeros((topo.nGlobal, 3))
        for iFace in xrange(topo.nFace):
            for j in xrange(topo.l_index[iFace].shape[1]):
                for i in xrange(topo.l_index[iFace].shape[0]):
                    self.X[topo.l_index[iFace][i,j]] = geo_obj.surfs[iFace].X[i,j]
                # end for
            # end for
        # end for

        # Finally we have to figure out the halos for the unstructured
        # nodes. This is a little tricky since we don't know how many
        # neighbours each node will have

        # Create a data structure that stores neighbours of global nodes:
        gnn = [[] for i in xrange(topo.nGlobal)] # GlobalNodeNearest
        gne = [[] for i in xrange(topo.nGlobal)] # GlobalNodeExtra
        gnorm = numpy.zeros((topo.nGlobal, 3))
        for iFace in xrange(topo.nFace):
            for j in xrange(topo.l_index[iFace].shape[1]-1):
                for i in xrange(topo.l_index[iFace].shape[0]-1):
                    # Each element has 4 edges, add the connections to
                    # each global node. Don't worry about duplicates,
                    # we'll take care of that later

                    # Lower edge
                    gnn[topo.l_index[iFace][i, j]].append(topo.l_index[iFace][i+1, j])
                    gnn[topo.l_index[iFace][i+1, j]].append(topo.l_index[iFace][i, j])

                    # upper edge
                    gnn[topo.l_index[iFace][i, j+1]].append(topo.l_index[iFace][i+1, j+1])
                    gnn[topo.l_index[iFace][i+1, j+1]].append(topo.l_index[iFace][i, j+1])

                    # left edge
                    gnn[topo.l_index[iFace][i, j]].append(topo.l_index[iFace][i, j+1])
                    gnn[topo.l_index[iFace][i, j+1]].append(topo.l_index[iFace][i, j])

                    # right edge
                    gnn[topo.l_index[iFace][i+1, j]].append(topo.l_index[iFace][i+1, j+1])
                    gnn[topo.l_index[iFace][i+1, j+1]].append(topo.l_index[iFace][i+1, j])

                    # And also do the diagonals - lower left <-> Upper right
                    # gne[topo.l_index[iFace][i, j]].append(topo.l_index[iFace][i+1, j+1])
                    # gne[topo.l_index[iFace][i+1, j+1]].append(topo.l_index[iFace][i, j])

                    # # And also do the diagonals - lower right <-> upper left
                    # gne[topo.l_index[iFace][i+1, j]].append(topo.l_index[iFace][i, j+1])
                    # gne[topo.l_index[iFace][i, j+1]].append(topo.l_index[iFace][i+1, j])
                    
                    # Also get the norm and scatter to nodes
                    ll = self.X[topo.l_index[iFace][i,j]]
                    lr = self.X[topo.l_index[iFace][i+1,j]]
                    ul = self.X[topo.l_index[iFace][i,j+1]]
                    ur = self.X[topo.l_index[iFace][i+1, j+1]]
                    
                    # Area normal
                    nrm = 0.5*numpy.cross((ur-ll),(ul-lr))
                    
                    # Scatter back to the nodes
                    gnorm[topo.l_index[iFace][i  ,j  ]] = 0.25*nrm
                    gnorm[topo.l_index[iFace][i+1,j  ]] = 0.25*nrm
                    gnorm[topo.l_index[iFace][i  ,j+1]] = 0.25*nrm
                    gnorm[topo.l_index[iFace][i+1,j+1]] = 0.25*nrm                                

                # end for
            # end for
        # end for

        # Now loop back through nPtr and we have the information we
        # need for the corners that have 4 nodes on them
        for i in xrange(len(nPtr)):
            if nPtr[i] == []:
                nearestNeighbours = geo_utils.unique(gnn[i]) 
                extraNeighbours = geo_utils.unique(gne[i]) 
                if len(nearestNeighbours) == 4:
                    nPtr[i].append(0)
                    nPtr[i].append(4)
                else:
                    nPtr[i].append(0)
                    nPtr[i].append(len(nearestNeighbours) + len(extraNeighbours))
                    nearestNeighbours.extend(extraNeighbours)
                # end if

                # Before we can add the neighbours we need to
                # geometriclly sort them by using the gnorm we computed earlier
                nrm = gnorm[i]/numpy.linalg.norm(gnorm[i])
                
                # We can arbitrarily take the first local direction to
                # be from the point to the first first node. 
                
                x1 = self.X[nearestNeighbours[0]] - self.X[i]
                x1 /= numpy.linalg.norm(x1)

                # Get x2 by crossing with nrm which is already normalized
                x2 = numpy.cross(nrm, x1)
                
                # Project each vector onto x1 and x2 to get local
                # coordinates and then determine angle
                theta = []
                for ii in xrange(len(nearestNeighbours)):
                    dr = self.X[nearestNeighbours[ii]] - self.X[i]
                    dx = numpy.dot(dr, x1)
                    dy = numpy.dot(dr, x2)
                    l = numpy.sqrt(dx*dx + dy*dy)
                    if dy > 0:
                        theta.append(numpy.arccos(dx/l))
                    else:
                        theta.append(2*numpy.pi - numpy.arccos(dx/l))
                    # end if
                # end for

                # Argsort the thetas to get the permutation
                tmp = numpy.argsort(numpy.array(theta))
                for ii in xrange(len(tmp)):
                    nPtr[i].append(nearestNeighbours[tmp[ii]] + 1)
                # end for
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
            nPtrArray[i, 0:l] = nPtr[i][:]
        # # end for
      
        # Now we have all the data we need so we can go ahead and
        # initialze the 3D generation

        self.hyp.init3d(sizes, nPtrArray.T)
        for i in xrange(topo.nFace):
            self.hyp.setlindex(topo.l_index[i], i)
        # end for

        #Ouput a tecplot FE mesh
        nNodes = topo.nGlobal
        nElem = 0
        for iFace in xrange(topo.nFace):
            nElem += (topo.l_index[iFace].shape[0]-1)*(topo.l_index[iFace].shape[1]-1)
        f = open('fe_mesh.dat','w')
        f.write("FE Data\n")
        f.write('VARIABLES = \"X\", \"Y\", \"z\" \"nx\" \"ny\" \"nz\" \n')
        f.write("ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n"%(nNodes, nElem))
        for i in xrange(nNodes):
            # Compute the normal
            if nPtrArray[i,0] <> 0:
                jp1 = nPtrArray[i,0]-1
                jm1 = nPtrArray[i,1]-1
                kp1 = nPtrArray[i,2]-1
                km1 = nPtrArray[i,3]-1
                dx = self.X[jp1] - self.X[jm1]
                dy = self.X[kp1] - self.X[km1]
                n = numpy.cross(dx,dy)
            else:
                n = [0,0,0]

            f.write('%f %f %f %f %f %f\n'%(self.X[i,0], self.X[i,1], self.X[i,2],
                                           n[0], n[1], n[2]))
        # end for
        for iFace in xrange(topo.nFace):
            for j in xrange(topo.l_index[iFace].shape[1]-1):
                for i in xrange(topo.l_index[iFace].shape[0]-1):
                    f.write('%d %d %d %d\n'%(topo.l_index[iFace][i, j]+1,
                                             topo.l_index[iFace][i+1, j]+1,
                                             topo.l_index[iFace][i+1, j+1]+1,
                                             topo.l_index[iFace][i, j+1]+1))
                # end for
            # end for
        # end for
        
        # Close file
        f.close()

        # for iFace in xrange(topo.nFace):
        #     print '-------------'
        #     print 'iFace:',iFace
        #     print '-------------'
        #     print topo.l_index[iFace]

        # for i in xrange(len(nPtr)):
        #     print i, nPtr[i]

        

        return

    def run(self):
        """
        Run given using the options given
        """
        if self.twoD:
            self.hyp.run2d(self.X.T)
            Xnew = None
        else:
            Xnew = self.hyp.run3d(self.X.T).T
        # end if
        self.gridGenerated = True

        return Xnew

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
            # end if
        else:
            mpiPrint('Error! No grid has been generated! Run the run() \
command before trying to write the grid!')
        # end if

        return
      
    def _setOptions(self):
        """Internal function to set the options in pyHyp"""
        self.hyp.hypinput.n         = self.options['N']
        self.hyp.hypinput.ndebug    = self.options['Ndebug']
        self.hyp.hypinput.s0        = self.options['s0']
        self.hyp.hypinput.gridratio = self.options['gridRatio']
        self.hyp.hypinput.epse      = self.options['epsE']
        self.hyp.hypinput.epsi      = self.options['epsI']
        self.hyp.hypinput.theta     = self.options['theta']
        self.hyp.hypinput.volcoef   = self.options['volCoef']
        self.hyp.hypinput.volblend  = self.options['volBlend']
        self.hyp.hypinput.volsmoothiter = self.options['volSmoothIter']

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

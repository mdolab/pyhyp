# =============================================================================
# Utility Functions for use in pyHyp
# =============================================================================

import numpy as np

def e_dist(x1, x2):
    '''Get the eculidean distance between two points'''
    if len(x1) == 3:
        return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
    elif len(x1) == 2:
        return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
    elif len(x1) == 1:
        return np.abs(x1[0]-x2[0])

def reverseRows(input):
    '''Flip Rows (horizontally)'''
    rows = input.shape[0]
    cols = input.shape[1]
    output = np.empty([rows, cols], input.dtype)
    for row in xrange(rows):
        output[row] = input[row][::-1].copy()
    # end for

    return output

def readNValues(handle, N, dtype, binary=False, sep=' '):
    '''Read N values of dtype 'float' or 'int' from file handle'''
    if binary == True:
        sep = ""
    # end if

    if dtype == 'int':
        values = np.fromfile(handle, dtype='int32', count=N, sep=sep)
    else:
        values = np.fromfile(handle, dtype='float', count=N, sep=sep)
    return values

def unique(s):
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't np.cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.

    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.

    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

def unique_index(s, s_hash=None):
    '''
    This function is based on unique

    The idea is to take a list s, and reduce it as per unique.

    However, it additionally calculates a linking index arrary that is
    the same size as the original s, and points to where it ends up in
    the the reduced list

    if s_hash is not specified for sorting, s is used

    '''
    if s_hash != None:
        ind = np.argsort(np.argsort(s_hash))
    else:
        ind = np.argsort(np.argsort(s))
    # end if
    n = len(s)
    t = list(s)
    t.sort()
    
    diff = np.zeros(n, 'bool')

    last = t[0]
    lasti = i = 1
    while i < n:
        if t[i] != last:
            t[lasti] = last = t[i]
            lasti += 1
        else:
            diff[i] = True
        # end if
        i += 1
    # end while
    b = np.where(diff == True)[0]
    for i in xrange(n):
        ind[i] -= b.searchsorted(ind[i], side='right')
    # end for

    return t[:lasti], ind

def pointReduce(points, node_tol=1e-4):
    '''Given a list of N points in ndim space, with possible
    duplicates, return a list of the unique points AND a pointer list
    for the original points to the reduced set'''

    # First 
    points = np.array(points)
    N = len(points)
    dists = []
    for ipt in xrange(N): 
        dists.append(np.sqrt(np.dot(points[ipt], points[ipt])))
    # end for
    temp = np.array(dists)
    temp.sort()
    ind = np.argsort(dists)
    i = 0
    cont = True
    new_points = []
    link = np.zeros(N, 'intc')
    link_counter = 0
   
    while cont:
        cont2 = True
        temp_ind = []
        j = i
        while cont2:
            if abs(dists[ind[i]]-dists[ind[j]])<node_tol:
                temp_ind.append(ind[j])
                j = j + 1
                if j == N: # Overrun check
                    cont2 = False
                # end if
            else:
                cont2 = False
            #end if
        # end while
        sub_points = [] # Copy of the list of sub points with the dists
        for ii in xrange(len(temp_ind)):
            sub_points.append(points[temp_ind[ii]])

        # Brute Force Search them 
        sub_unique_pts, sub_link = pointReduceBruteForce(sub_points, node_tol)
        new_points.extend(sub_unique_pts)

        for ii in xrange(len(temp_ind)):
            link[temp_ind[ii]] = sub_link[ii] + link_counter
        # end if
        link_counter += max(sub_link) + 1

        i = j - 1 + 1
        if i == N:
            cont = False
        # end if
    # end while
    return np.array(new_points), np.array(link)

def pointReduceBruteForce(points,  node_tol=1e-4):
    '''Given a list of N points in ndim space, with possible
    duplicates, return a list of the unique points AND a pointer list
    for the original points to the reduced set

    BRUTE FORCE VERSION

    '''
    N = len(points)
    unique_points = [points[0]]
    link = [0]
    for i in xrange(1, N):
        found_it = False
        for j in xrange(len(unique_points)):
            if e_dist(points[i], unique_points[j]) < node_tol:
                link.append(j)
                found_it = True
                break
            # end if
        # end for
        if not found_it:
            unique_points.append(points[i])
            link.append(j+1)
        # end if
    # end for
    return np.array(unique_points), np.array(link)

# --------------------------------------------------------------
#                     Node/Edge Functions
# --------------------------------------------------------------

def nodesFromEdge(edge):
    '''Return the nodes on either edge of a standard edge'''
    if edge == 0:
        return 0, 1
    elif edge == 1:
        return 2, 3
    elif edge == 2:
        return 0, 2
    elif edge == 3:
        return 1, 3
    elif edge == 4:
        return 4, 5
    elif edge == 5:
        return 6, 7
    elif edge == 6:
        return 4, 6
    elif edge == 7:
        return 5, 7
    elif edge == 8:
        return 0, 4
    elif edge == 9:
        return 1, 5
    elif edge == 10:
        return 2, 6
    elif edge == 11:
        return 3, 7

class topology(object):
    '''
    The base topology class from which the BlockTopology,
    SurfaceTology and CuveTopology classes inherit from
    
    The topology object contains all the info required for the block
    topology (most complex) however, simpiler topologies are handled
    accordingly.

    Class Attributes:
        nVol : The number of volumes in the topology (may be 0)
        nFace: The number of unique faces on the topology (may be 0)
        nEdge: The number of uniuqe edges on the topology 
        nNode: The number of unique nodes on the topology

        nEnt: The number of "entities" in the topology class. This may
        be curves, faces or volumes

        mNodeEnt: The number of NODES per entity. For curves it's 2, for
        surfaces 4 and for volumes 8.

        mEdgeEnt: The number of EDGES per entity. For curves it's 1,
        for surfaces, 4 and for volumes, 12

        mFaceEnt: The number of faces per entity. For curves its's 0,
        for surfaces, 1 and for volumes,6

        mVolEnt: The number of volumes per entity. For curves it's 0,
        for surfaces, 0 and for volumnes, 1

        node_link: The array of size nEnt x mNodesEnt which points
                   to the node for each entity
        edge_link: The array of size nEnt x mEdgeEnt which points
                   to the edge for each edge of entity
        face_link: The array of size nEnt x mFaceEnt which points to 
                   the face of each face on an entity

        edge_dir:  The array of size nEnt x mEdgeEnt which detrmines
                   if the intrinsic direction of this edge is
                   opposite of the direction as recorded in the
                   edge list. edge_dir[entity#][#] = 1 means same direction;
                   -1 is opposite direction.
                  
        face_dir:  The array of size nFace x 6 which determines the 
                   intrinsic direction of this face. It is one of 0->7
                   
        l_index:   The local->global list of arrays for each volue
        g_index:   The global->local list points for the entire topology
        edges:     The list of edge objects defining the topology
        simple    : A flag to determine of this is a "simple" topology 
                   which means there are NO degernate Edges, 
                   NO multiple edges sharing the same nodes and NO 
                   edges which loop back and have the same nodes
                   MUST BE SIMPLE
    '''

    def __init__(self):
        # Not sure what should go here...
        return
    def _calcDGs(self, edges, edge_link, edge_link_sorted, edge_link_ind):

        dg_counter = -1
        for i in xrange(self.nEdge):
            if edges[i][2] == -1: # Not set yet
                dg_counter += 1
                edges[i][2] = dg_counter
                self._addDGEdge(i, edges, edge_link, 
                                edge_link_sorted, edge_link_ind)
            # end if
        # end for
        self.nDG = dg_counter + 1
   
    def _addDGEdge(self, i, edges, edge_link, edge_link_sorted, edge_link_ind):
        left  = edge_link_sorted.searchsorted(i, side='left')
        right = edge_link_sorted.searchsorted(i, side='right')
        res   = edge_link_ind[slice(left, right)]

        for j in xrange(len(res)):
            ient = res[j]/self.mEdgeEnt #Integer Division
            iedge = np.mod(res[j], self.mEdgeEnt)

            pEdges = self._getParallelEdges(iedge)
            oppositeEdges = []
            for iii in xrange(len(pEdges)):
                oppositeEdges.append(
                    edge_link[self.mEdgeEnt*ient + pEdges[iii]])
            
            for ii in xrange(len(pEdges)):
                if edges[oppositeEdges[ii]][2] == -1:
                    edges[oppositeEdges[ii]][2] = edges[i][2]
                    if not edges[oppositeEdges[ii]][0] == \
                            edges[oppositeEdges[ii]][1]:
                        self._addDGEdge(oppositeEdges[ii], edges, 
                                        edge_link, edge_link_sorted, 
                                        edge_link_ind)
                # end if
            # end if
        # end for

    def _getParallelEdges(self, iedge):
        '''Return parallel edges for surfaces and volumes'''

        if self.topo_type == 'surface':
            if iedge == 0: return [1]
            if iedge == 1: return [0]
            if iedge == 2: return [3]
            if iedge == 3: return [2]

        if self.topo_type == 'volume':
            if iedge == 0: 
                return [1, 4, 5]
            if iedge == 1: 
                return [0, 4, 5]
            if iedge == 2: 
                return [3, 6, 7]
            if iedge == 3: 
                return [2, 6, 7]
            if iedge == 4: 
                return [0, 1, 5]
            if iedge == 5: 
                return [0, 1, 4]
            if iedge == 6: 
                return [2, 3, 7]
            if iedge == 7: 
                return [2, 3, 6]
            if iedge == 8: 
                return [9, 10, 11]
            if iedge == 9: 
                return [8, 10, 11]
            if iedge == 10: 
                return [8, 9, 11]
            if iedge == 11: 
                return [8, 9, 10]
        if self.topo_type == 'curve':
            return None

    def _getDGList(self):
        '''After calcGlobalNumbering is called with the size
        parameters, we can now produce a list of length ndg with the
        each entry coorsponing to the number N associated with that DG'''

        # This can be run in linear time...just loop over each edge
        # and add to dg list
        N_list = np.zeros(self.nDG, 'intc')
        for iedge in xrange(self.nEdge):
            N_list[self.edges[iedge].dg] = self.edges[iedge].N
        # end for
            
        return N_list

class SurfaceTopology(topology):
    '''
    See topology class for more information
    '''
    def __init__(self, coords=None, face_con=None, file=None, node_tol=1e-4,
                 edge_tol=1e-4, **kwargs):
        '''Initialize the class with data required to compute the topology'''
        topology.__init__(self)
        self.mNodeEnt = 4
        self.mEdgeEnt = 4
        self.mfaceEnt = 1
        self.mVolEnt  = 0
        self.nVol = 0
        self.topo_type = 'surface'
        self.g_index = None
        self.l_index = None
        self.nGlobal = None
        self.edges = None
        self.face_index = None
        self.simple = False

        if not face_con == None: 
            face_con = np.array(face_con)
            midpoints = None
            self.nFace = len(face_con)
            self.nEnt = self.nFace
            self.simple = True
            # Check to make sure nodes are sequential
            self.nNode = len(unique(face_con.flatten()))
            if self.nNode != max(face_con.flatten())+1:
                # We don't have sequential nodes
                mpiPrint("Error: Nodes are not sequential")
                sys.exit(1)
            # end if
            
            edges = []
            edge_hash = []
            for iface in xrange(self.nFace):
                #             n1                ,n2               ,dg,n,degen
                edges.append([face_con[iface][0], face_con[iface][1], -1, 0, 0])
                edges.append([face_con[iface][2], face_con[iface][3], -1, 0, 0])
                edges.append([face_con[iface][0], face_con[iface][2], -1, 0, 0])
                edges.append([face_con[iface][1], face_con[iface][3], -1, 0, 0])
            # end for
            edge_dir = np.ones(len(edges), 'intc')
            for iedge in xrange(self.nFace*4):
                if edges[iedge][0] > edges[iedge][1]:
                    temp = edges[iedge][0]
                    edges[iedge][0] = edges[iedge][1]
                    edges[iedge][1] = temp
                    edge_dir[iedge] = -1
                # end if
                edge_hash.append(
                    edges[iedge][0]*4*self.nFace + edges[iedge][1])
            # end for

            edges, edge_link = unique_index(edges, edge_hash)

            self.nEdge = len(edges)
            self.edge_link = np.array(edge_link).reshape((self.nFace, 4))
            self.node_link = np.array(face_con)
            self.edge_dir  = np.array(edge_dir).reshape((self.nFace, 4))

            edge_link_sorted = np.sort(edge_link)
            edge_link_ind    = np.argsort(edge_link)

        elif not coords == None:
            self.nFace = len(coords)
            self.nEnt  = self.nFace
            # We can use the pointReduce algorithim on the nodes
            node_list, node_link = pointReduce(
                coords[:, 0:4, :].reshape((self.nFace*4, 3)), node_tol)
            node_link = node_link.reshape((self.nFace, 4))
          
            # Next Calculate the EDGE connectivity. -- This is Still
            # Brute Force

            edges = []
            midpoints = []
            edge_link = -1*np.ones(self.nFace*4, 'intc')
            edge_dir  = np.zeros((self.nFace, 4), 'intc')

            for iface in xrange(self.nFace):
                for iedge in xrange(4):
                    n1, n2 = nodesFromEdge(iedge)
                    n1 = node_link[iface][n1]
                    n2 = node_link[iface][n2] 
                    midpoint = coords[iface][iedge + 4]
                    if len(edges) == 0:
                        edges.append([n1, n2, -1, 0, 0])
                        midpoints.append(midpoint)
                        edge_link[4*iface + iedge] = 0
                        edge_dir [iface][iedge] = 1
                    else:
                        found_it = False
                        for i in xrange(len(edges)):
                            if [n1, n2] == edges[i][0:2] and n1 != n2:
                                if e_dist(midpoint, midpoints[i]) < edge_tol:
                                    edge_link[4*iface + iedge] = i
                                    edge_dir [iface][iedge] = 1
                                    found_it = True
                                # end if
                            elif [n2, n1] == edges[i][0:2] and n1 != n2:
                                if e_dist(midpoint, midpoints[i]) < edge_tol:
                                    edge_link[4*iface + iedge] = i
                                    edge_dir[iface][iedge] = -1
                                    found_it = True
                                # end if
                            # end if
                        # end for

                        # We went all the way though the list so add
                        # it at end and return index
                        if not found_it:
                            edges.append([n1, n2, -1, 0, 0])
                            midpoints.append(midpoint)
                            edge_link[4*iface + iedge] = i+1
                            edge_dir [iface][iedge] = 1
                    # end if
                # end for
            # end for

            self.nEdge = len(edges)
            self.edge_link = np.array(edge_link).reshape((self.nFace, 4))
            self.node_link = np.array(node_link)
            self.nNode = len(unique(self.node_link.flatten()))
            self.edge_dir = edge_dir

            edge_link_sorted = np.sort(edge_link.flatten())
            edge_link_ind    = np.argsort(edge_link.flatten())

        # end if
            
        # Next Calculate the Design Group Information
        self._calcDGs(edges, edge_link, edge_link_sorted, edge_link_ind)

        # Set the edge ojects
        self.edges = []
        for i in xrange(self.nEdge): # Create the edge objects
            if midpoints: # If they exist
                if edges[i][0] == edges[i][1] and \
                        e_dist(midpoints[i], node_list[edges[i][0]]) < node_tol:
                    self.edges.append(edge(edges[i][0], edges[i][1], 
                                           0, 1, 0, edges[i][2], edges[i][3]))
                else:
                    self.edges.append(edge(edges[i][0], edges[i][1], 
                                           0, 0, 0, edges[i][2], edges[i][3]))
                # end if
            else:
                self.edges.append(edge(edges[i][0], edges[i][1], 
                                       0, 0, 0, edges[i][2], edges[i][3]))
            # end if
        # end for

        return

    def calcGlobalNumbering(self, sizes, surface_list=None):
        '''Internal function to calculate the global/local numbering
        for each surface'''
        for i in xrange(len(sizes)):
            self.edges[self.edge_link[i][0]].N = sizes[i][0]
            self.edges[self.edge_link[i][1]].N = sizes[i][0]
            self.edges[self.edge_link[i][2]].N = sizes[i][1]
            self.edges[self.edge_link[i][3]].N = sizes[i][1]

        if surface_list == None:
            surface_list = range(0, self.nFace)
        # end if
        
        # ----------------- Start of Edge Computation ---------------------
        counter = 0
        g_index = []
        l_index = []

        assert len(sizes) == len(surface_list), 'Error: The list of sizes and \
the list of surfaces must be the same length'

        # Assign unique numbers to the corners -> Corners are indexed
        # sequentially
        node_index = np.arange(self.nNode)
        counter = len(node_index)
        edge_index = [ [] for i in xrange(len(self.edges))]
     
        # Assign unique numbers to the edges

        for ii in xrange(len(surface_list)):
            cur_size = [sizes[ii][0], sizes[ii][0], sizes[ii][1], sizes[ii][1]]
            isurf = surface_list[ii]
            for iedge in xrange(4):
                edge = self.edge_link[ii][iedge]
                    
                if edge_index[edge] == []:# Not added yet
                    if self.edges[edge].degen == 1:
                        # Get the counter value for this "node"
                        index = node_index[self.edges[edge].n1]
                        for jj in xrange(cur_size[iedge]-2):
                            edge_index[edge].append(index)
                        # end for
                    else:
                        for jj in xrange(cur_size[iedge]-2):
                            edge_index[edge].append(counter)
                            counter += 1
                        # end for
                    # end if
                # end if
            # end for
        # end for
     
        g_index = [ [] for i in xrange(counter)] # We must add [] for
                                                 # each global node
        l_index = []

        # Now actually fill everything up

        for ii in xrange(len(surface_list)):
            isurf = surface_list[ii]
            N = sizes[ii][0]
            M = sizes[ii][1]
            l_index.append(-1*np.ones((N, M), 'intc'))

            for i in xrange(N):
                for j in xrange(M):
                    
                    type, edge, node, index = indexPosition2D(i, j, N, M)

                    if type == 0:           # Interior
                        l_index[ii][i, j] = counter
                        g_index.append([[isurf, i, j]])
                        counter += 1
                    elif type == 1:         # Edge
                       
                        if edge in [0, 1]:
                            # Its a reverse dir
                            if self.edge_dir[ii][edge] == -1:
                                cur_index = edge_index[
                                    self.edge_link[ii][edge]][N-i-2]
                            else:  
                                cur_index = edge_index[
                                    self.edge_link[ii][edge]][i-1]
                            # end if
                        else: # edge in [2, 3]
                            # Its a reverse dir
                            if self.edge_dir[ii][edge] == -1: 
                                cur_index = edge_index[
                                    self.edge_link[ii][edge]][M-j-2]
                            else:  
                                cur_index = edge_index[
                                    self.edge_link[ii][edge]][j-1]
                            # end if
                        # end if
                        l_index[ii][i, j] = cur_index
                        g_index[cur_index].append([isurf, i, j])
                            
                    else:                  # Node
                        cur_node = self.node_link[ii][node]
                        l_index[ii][i, j] = node_index[cur_node]
                        g_index[node_index[cur_node]].append([isurf, i, j])
                    # end for
                # end for (j)
            # end for (i)
        # end for (ii)

        # Reorder the indices with a greedy scheme

        new_indices = np.zeros(len(g_index), 'intc')
        new_indices[:] = -1
        new_g_index = [[] for i in xrange(len(g_index))]
        counter = 0

        # Re-order the l_index
        for ii in xrange(len(surface_list)):
            isurf = surface_list[ii]
            N = sizes[ii][0]
            M = sizes[ii][1]
            for i in xrange(N):
                for j in xrange(M):
                    if new_indices[l_index[ii][i, j]] == -1:
                        new_indices[l_index[ii][i, j]] = counter
                        l_index[ii][i, j] = counter 
                        counter += 1
                    else:
                        l_index[ii][i, j] = new_indices[l_index[ii][i, j]]
                    # end if
                # end for
            # end for
        # end for
       
        # Re-order the g_index
        for ii in xrange(len(g_index)):
            isurf = g_index[ii][0][0]
            i     = g_index[ii][0][1]
            j     = g_index[ii][0][2]
            pt = l_index[isurf][i, j]
            new_g_index[pt] = g_index[ii]
            # end for
        # end for
            
        self.nGlobal = len(g_index)
        self.g_index = new_g_index
        self.l_index = l_index
        
        return 


class edge(object):
    '''A class for edge objects'''

    def __init__(self, n1, n2, cont, degen, intersect, dg, N):
        self.n1        = n1        # Integer for node 1
        self.n2        = n2        # Integer for node 2
        self.cont      = cont      # Integer: 0 for c0 continuity, 1
                                   # for c1 continuity
        self.degen     = degen     # Integer: 1 for degenerate, 0 otherwise
        self.intersect = intersect # Integer: 1 for an intersected
                                   # edge, 0 otherwise
        self.dg        = dg        # Design Group index
        self.N         = N         # Number of control points for this edge
        
    def write_info(self, i, handle):
        handle.write('  %5d        | %5d | %5d | %5d | %5d | %5d |\
  %5d |  %5d |\n'\
                     %(i, self.n1, self.n2, self.cont, self.degen, 
                       self.intersect, self.dg, self.N))

class edge_cmp_object(object):
    '''A temporary class for sorting edge objects'''

    def __init__(self, n1, n2, n1o, n2o, mid_pt, tol):
        self.n1 = n1
        self.n2 = n2
        self.nodes = [n1o, n2o]
        self.mid_pt = mid_pt
        self.tol = tol

    def __repr__(self):
        return 'Node1: %d Node2: %d Mid_pt: %f %f %f'% (
            self.n1, self.n2, self.mid_pt[0], self.mid_pt[1], self.mid_pt[2])

    def __cmp__(self, other):
        # This function should return :
        # -1 if self < other
        #  0 if self == other
        #  1 if self > other

        # Basically we want to make three comparisons: n1, n2 and the
        # mid_pt Its (substantially) faster if we break before all 3
        # comparisons are done
        
        n1_cmp = cmp(self.n1, other.n1)
        if n1_cmp: # n1_cmp is non-zero so return with the result
            return n1_cmp

        n2_cmp = cmp(self.n2, other.n2)

        if n2_cmp: # n2_cmp is non-zero so return 
            return n2_cmp

        x_cmp = cmp(self.mid_pt[0], other .mid_pt[0])
        y_cmp = cmp(self.mid_pt[1], other .mid_pt[1])
        z_cmp = cmp(self.mid_pt[2], other .mid_pt[2])
        
        if e_dist(self.mid_pt, other.mid_pt) < self.tol:
            mid_cmp = 0
        else:
            mid_cmp = x_cmp or y_cmp or z_cmp
        # end if

        return mid_cmp

def indexPosition2D(i, j, N, M):
    '''This function is a generic function which determines if for a grid
    of data NxM with index i going 0->N-1 and j going 0->M-1, it
    determines if i,j is on the interior, on an edge or on a corner

    The funtion return four values: 
    type: this is 0 for interior, 1 for on an edge and 2 for on a corner
    edge: this is the edge number if type==1
    node: this is the node number if type==2 
    index: this is the value index along the edge of interest -- 
    only defined for edges'''

    if i > 0 and i < N - 1 and j > 0 and j < M-1: # Interior
        return 0, None, None, None
    elif i > 0 and i < N - 1 and j == 0:     # Edge 0
        return 1, 0, None, i
    elif i > 0 and i < N - 1 and j == M - 1: # Edge 1
        return 1, 1, None, i
    elif i == 0 and j > 0 and j < M - 1:     # Edge 2
        return 1, 2, None, j
    elif i == N - 1 and j > 0 and j < M - 1: # Edge 3
        return 1, 3, None, j
    elif i == 0 and j == 0:                  # Node 0
        return 2, None, 0, None
    elif i == N - 1 and j == 0:              # Node 1
        return 2, None, 1, None
    elif i == 0 and j == M -1 :              # Node 2
        return 2, None, 2, None
    elif i == N - 1 and j == M - 1:          # Node 3
        return 2, None, 3, None

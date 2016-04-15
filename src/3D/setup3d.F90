subroutine setup(fileName, fileType)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: setup3d is the main routine for setting up the 3D
  !     elliptic or hyperbolic mesh extrusion problem. The only input
  !     is the filename of an ascii-formatted plot3d file and
  !     information about a possible mirror. 
  !    
  !     Ney Secco (2015-2016): Added boundary conditions and CGNS interface
  !
  !     Description of Arguments
  !     Input:
  !     fileName : char* array : Filename of plot3d or cgns file.
  !     fileType : integer : 1-CGNS 2-Plot3d

  use communication
  use hypData
  use hypInput

  implicit none

  ! Input Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: fileType

  ! Working parameters
  real(kind=realType), dimension(:, :), allocatable :: uniquePts, allEdgeNodes
  integer(kind=intType), dimension(:), allocatable :: link, ptInd, nFace, localFace
  integer(kind=intType), dimension(:), allocatable :: nte, ntePtr
  integer(kind=intType), dimension(:, :), allocatable :: fullnPtr
  integer(kind=intType), dimension(:, :), allocatable :: nodeConn
  integer(kind=intType), dimension(:, :), allocatable ::  directedNodeConn
  integer(kind=intType), dimension(:), allocatable ::  nEdge, nDirectedEdge
  
  integer(kind=intType) :: nBlocks, nUnique, corners(4), iCorner
  integer(kind=intType) :: nodeTotal, node1, node2, node3, node4
  integer(kind=intType) :: iNode, evenNodes
  integer(kind=intType) :: NODE_MAX, CELL_MAX
  integer(kind=intType) :: nFaceNeighbours, firstID, nodeID, nextElemID
  integer(kind=intType) :: nodeToAdd, faceToCheck, nodeToFind
  integer(kind=intType) :: i, j, k, ii, jj, iii, jjj, ierr, iStart, iEnd
  integer(kind=intType) :: isize, ind, icell, idim, BCToSet, il, jl, iEdge, iBCToSet
  logical :: found
  character(5) :: faceStr
  integer(kind=intType), dimension(:), pointer :: lPtr1, lPtr2
  real(kind=realType), dimension(:, :), pointer :: xPtr
  real(kind=realType) :: xSum(3), vals(2, 3)
  integer(kind=intType) :: mask(2, 3), nn(2), mm(2), f(3), iFace, jp2

  ! Only the root processor reads the mesh and does the initial processing:
  rootProc: if (myid == 0) then 
     
     fType: if (fileType .eq. cgnsFileType) then ! We have a CGNS file
        call readCGNS(fileName)
     else if (fileType .eq. plot3dFileType) then ! We have a Plot3d file
        call readPlot3d(fileName)
     end if fType

     ! Now we can create the final connectivity
     nodeTotal = 0
     faceTotal = 0
     do ii=1,nPatch
        nodeTotal = nodeTotal + patches(ii)%il*patches(ii)%jl
        faceTotal = faceTotal + (patches(ii)%il-1)*(patches(ii)%jl-1)
     end do

     allocate(allEdgeNodes(3, nodeTotal), uniquePts(3, nodeTotal), link(nodeTotal))

     ! Assign edge nodes to allEdgeNodes
     iNode = 0
     do ii=1, nPatch
        do j=1, patches(ii)%jl, patches(ii)%jl-1
           do i=1, patches(ii)%il
              iNode = iNode + 1
              allEdgeNodes(:, iNode) = patches(ii)%X(:, i, j)
           end do
        end do
        do i=1, patches(ii)%il, patches(ii)%il-1
           do j=2, patches(ii)%jl-1
              iNode = iNode + 1
              allEdgeNodes(:, iNode) = patches(ii)%X(:, i, j)
           end do
        end do
     end do

     ! Now we will call pointReduce to remove any repeated edge nodes.
     ! link maps from allEdgeNodes to uniquePts
     if (noPointReduce) then 
        uniquePts = allEdgeNodes
        nUnique = iNode
        do i=1, iNode
           link(i) = i
        end do
     else
        call pointReduce(allEdgeNodes, iNode, nodeTol, uniquePts, link, nUnique)
     end if

     ! Dump edge/corner into l_index
     iNode = 0
     do ii=1,nPatch
        do j=1, patches(ii)%jl, patches(ii)%jl-1
           do i=1, patches(ii)%il
              iNode = iNode + 1
              patches(ii)%l_index(i, j) = link(iNode)
           end do
        end do
        
        do i=1, patches(ii)%il, patches(ii)%il-1
           do j=2, patches(ii)%jl-1
              iNode = iNode + 1
              patches(ii)%l_index(i, j) = link(iNode)
           end do
        end do
     end do

     ! Now fill in the interior points into uniquePts/l_index
     iNode = nUnique !Up to here nUnique only has the number of unique edge nodes
     do ii=1,nPatch
        do j=2, patches(ii)%jl-1
           do i=2, patches(ii)%il-1
              iNode = iNode + 1
              uniquePts(:, iNode) = patches(ii)%X(:, i, j)
              patches(ii)%l_index(i, j) = iNode
           end do
        end do
     end do
     nUnique = iNode ! Now nUnique has the number of unique edge and interior nodes

     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") " Total Nodes: "
     write(*,"(I7,1x)",advance="no") nodeTotal
     print "(1x)"  
     write(*,"(a)", advance="no") " Unique Nodes:"
     write(*,"(I7,1x)",advance="no") nUnique
     print "(1x)"  
     write(*,"(a)", advance="no") " Total Faces: " 
     write(*,"(I7,1x)",advance="no") faceTotal
     print "(1x)" 
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"   

     ! Free up all allocated memory up to here:
     deallocate(link, allEdgeNodes)

     ! Create the conn array
     allocate(fullConn(4, faceTotal))
     jj = 0
     do ii=1,nPatch
        do j=1, patches(ii)%jl-1
           do i=1, patches(ii)%il-1
              jj = jj + 1
              fullConn(:, jj) = (/&
                   patches(ii)%l_index(i,   j), &
                   patches(ii)%l_index(i+1, j), &
                   patches(ii)%l_index(i+1, j+1), &
                   patches(ii)%l_index(i,   j+1)/)
           end do
        end do
     end do

     ! Determine if the normals are consistent:
     print *,'Normal orientation check ...'

     ! Allocate and initialize matrix that stores nodal connectivity directions
     ! I'm assuming that the maximum number of connections of a single node
     ! is 10
     allocate(directedNodeConn(10, nUnique))
     allocate(nodeConn(10, nUnique))
     allocate(nEdge(nUnique), nDirectedEdge(nUnique))
     nodeConn = 0
     directedNodeConn = 0
     nEdge = 0
     nDirectedEdge = 0
     ! Loop over each face and increment elementes from the connectivity
     ! matrix acording to the connectivity direction.
     do jj=1,faceTotal

        ! Get nodes indices for this face
        node1 = fullConn(1, jj)
        node2 = fullConn(2, jj)
        node3 = fullConn(3, jj)
        node4 = fullConn(4, jj)

        ! Assign connections acording to their orientations
        ! The two routines are in 3D_utilities.F90
        call assignN2NDirected(directedNodeConn, nDirectedEdge, node1, node2, 10, nUnique)
        call assignN2NDirected(directednodeConn, nDirectedEdge, node2, node3, 10, nUnique)
        call assignN2NDirected(directedNodeConn, nDirectedEdge, node3, node4, 10, nUnique)
        call assignN2NDirected(directedNodeConn, nDirectedEdge, node4, node1, 10, nUnique)

        !Here we must add two edges for each 
        call assignN2N(nodeConn, nEdge, node1, node2, 10, nUnique)
        call assignN2N(nodeConn, nEdge, node1, node4, 10, nUnique)
        
        call assignN2N(nodeConn, nEdge, node2, node3, 10, nUnique)
        call assignN2N(nodeConn, nEdge, node2, node1, 10, nUnique)
        
        call assignN2N(nodeConn, nEdge, node3, node4, 10, nUnique)
        call assignN2N(nodeConn, nEdge, node3, node2, 10, nUnique)
        
        call assignN2N(nodeConn, nEdge, node4, node1, 10, nUnique)
        call assignN2N(nodeConn, nEdge, node4, node3, 10, nUnique)
     end do

     ! Now check if we flagged any node during the assignments
     if (minval(directedNodeConn) .eq. -1) then

        ! Say that the normals are wrong
        print *,' ERROR: Normal directions may be wrong'
        print *,'        Check the following coordinates:'

        ! Print coordinates of nodes were the problem occurs
        do i=1, nUnique
           if (directedNodeConn(1, i) .eq. -1) then
              print *,uniquePts(:, i)
           end if
        end do
        stop
     end if

     ! We just determined the number of nodes conected to each
     ! node. Now we determine the number of faces connected to each
     ! node. These two pieces of information will determine the
     ! topology of the node.
     allocate(nFace(nUnique))
     nFace(:) = 0
     do i=1, faceTotal
        do ii=1, 4
           nFace(fullConn(ii, i)) = nFace(fullConn(ii, i)) + 1
        end do
     end do
     
     ! Here we can determine the topology of every node in the
     ! mesh. The topology refers to the number of faces and edges each
     ! node is connected to. 

     allocate(fullTopoType(nUnique))
     allocate(fullBCType(2, nUnique))
     allocate(fullBCVal(2, 2, nUnique)) ! We can have 2 BC, and each
                                        ! can have 2 values. 

     ! Initialze these to internal node, with defaultBC
     fullTopoType = topoInternal
     fullBCType = BCDefault
     fullBCVal = huge(one) ! If we accidently use this and shouldn't,
                           ! this should cause problems.

     do i=1, nUnique
        if (nEdge(i) == 2 .and. nFace(i) == 1) then 
           ! This is a corner
           fullTopoType(i) = topoCorner

        else if (nEdge(i) == 2 .and. nFace(i) == 2) then 
           ! This is a degenerate corner. We can treat this as an
           ! internal node by using diagonals
           fullTopoType(i) = topoInternal

        else if (nEdge(i) == 3 .and. nface(i) == 2) then
           ! This is a regular edge
           fullTopoType(i) = topoEdge

        else if (nEdge(i) == 4 .and. nFace(i) == 2) then 
           ! This is a bow-tie. We don't do that either.
           print *, "Come on. A bow-tie? Really? No, we can't do that'"
           stop

        else if (nEdge(i) == 4 .and. nFace(i) == 3) then 
           ! This is an L corner. In theory, we should be able to do
           ! this...but not yet
           fullTopoType(i) = topoLCorner
           
        else if (nEdge(i) == nFace(i)) then 
           ! This is an internal node so we're ok.
           fullTopoType(i) = topoInternal

        else
           ! This is some really wonky combination which we assume we
           ! can't deal with. 
           print *, "An unknown topology was encountered. I'm assuming we can't do that."
           stop
        end if
     end do

     ! The nte stucture is like this:
     !
     !  ntePtr(i)
     !     |
     !     V
     ! | face 1 | position of node i on face 1 | face 2 | position of node i on face 2 | ...
     !
     ! This is why ntePtr is 2*sum(nFace) long
  
     ! Allocate the nodeToElem (nte) data array and the ptr array:
     allocate(nte(2*sum(nFace)), ntePtr(nUnique+1))
     ntePtr(1) = 1
     do i=1, nUnique
        ntePtr(i+1) = ntePtr(i) + nFace(i)*2
     end do

     ! Now fill-up the array
     nFace(:) = 0
     do i=1, faceTotal
        do ii=1,4
           iNode = fullConn(ii, i)
           jj = nFace(iNode)
           ! We store the index of the face, and index that this node
           ! is on the face
           nte(ntePtr(iNode) + jj*2) = i
           nte(ntePtr(iNode) + jj*2+1) = ii

           ! We need to keep track of the number we've already added. 
           nFace(iNode) = nFace(iNode) + 1
        end do
     end do

     ! nFace is now no longer needed...all the necessary info is in ntePtr
     deallocate(nFace)
     ! Determine the maximum number of nodal neighbours. Node
     ! max must be at least 6 for the 3-corner nodes with a
     ! total of 6 cross-diagonal neighbours
     NODE_MAX = 6
     CELL_MAX = 0
     do i=1, nUnique
        NODE_MAX = max(NODE_MAX, (ntePtr(i+1) - ntePtr(i))/2)
        CELL_MAX = max(CELL_MAX, (ntePtr(i+1) - ntePtr(i))/2)
     end do

     ! Finally we can get the final node pointer structure we need:
     allocate(fullnPtr(NODE_MAX+1, nUnique), fullcPtr(CELL_MAX+1, nUnique))

     nodeLoop: do iNode=1, nUnique
        internalNode: if (fullTopoType(iNode) == topoInternal) then 
           nFaceNeighbours = (ntePtr(iNode+1)-ntePtr(iNode))/2

           ! Low-neighbour count nodes get doubles, others get all neighbours
           if (nFaceNeighbours == 2 .or. nFaceNeighbours ==3) then
              fullnPtr(1, iNode) = nFaceNeighbours*2
           else
              fullnPtr(1, iNode) = nFaceNeighbours
           end if
           fullcPtr(1, iNode) = nFaceNeighbours
        
           firstID = nte(ntePtr(iNode))
           nodeID = mod(nte(ntePtr(iNode)+1), 4) + 1
           nextElemID = -1
           iii = 2
           jjj = 2

           nodeLoopCycle: do while (nextElemID /= firstID)
              if (nextElemID == -1) &
                   nextElemID = firstID

              ! Append the next node along that edge:
              nodeToAdd = fullConn(mod(nodeID-1, 4)+1, nextElemID)  
              fullnPtr(iii, iNode) = nodeToAdd 
              fullcPtr(jjj, iNode) = nextElemID
              iii = iii + 1
              jjj = jjj + 1

              ! Add the diagonal if not regular node
              if (nFaceNeighbours == 3 .or. nFaceNeighbours == 2) then
                 fullnPtr(iii, iNode) =  fullConn(mod(nodeID, 4)+1, nextElemID)
                 iii = iii + 1
              end if
              
              ! And we need to find the face that contains the following node:
              nodeToFind = fullConn(mod(nodeID+1, 4)+1, nextElemID)
              
              found = .False.
              jj = -1
              findNextElement: do while (.not. found)
                 jj = jj + 1
                 
                 ! Find the next face
                 faceToCheck = nte(ntePtr(iNode) + jj*2)
                 if (faceToCheck /= nextElemID) then
                    !  Check the 4 nodes on this face 
                    do ii=1, 4
                       if (fullConn(ii, faceToCheck) == nodeToFind) then
                          nextElemID = faceToCheck
                          nodeID = ii
                          found = .True.
                       end if
                    end do
                 end if
              end do findNextElement
           end do nodeLoopCycle
        else
           ! Not an internal node. Set all values to zero so we know
           ! they have not been processed yet:
           fullNPtr(:, iNode) = 0
           fullCPtr(:, iNode) = 0
        end if internalNode
     end do nodeLoop

     ! Check to make sure that there are only edge and interneral
     ! topology nodes if unattachedEdgesAreSymmetry is true: This can
     ! only be used for mirrored cases where all non interior nodes
     ! MUST be edges
     if (unattachedEdgesAreSymmetry) then 
        do i=1, nUnique
           if (fullTopoType(i) /= topoInternal .and. fullTopoType(i) /= topoEdge) then
              print *,"ERROR: A free corner or other topology that is not an "
              print *,"edge was detected with the unattachedEdgesAreSymmetry option."
              print *,"This option can only be used for configurations that become "
              print *,"closed when mirrored."
              stop
           end if
        end do
     end if

     ! We have now set the fullNodePointer (fullnPtr) and the
     ! fullCellPointer (fullCPtr) for all internal nodes. These were
     ! straightforward since by definition they are fully surrounded
     ! by elements. We now have to set the nodes on boundaries.  Since
     ! we will need to know the local structured topology it is easier
     ! to loop through the edges of the patches themselves and push
     ! the information to global node.

     ! First figure out which of the BCDefault edges are actually
     ! physical edges and assign them to be BCSplay OR Symmetry if the
     ! 'unattachedEdgesAreSymmetry" is true.
     patchLoop0: do ii=1, nPatch
        il = patches(ii)%il
        jl = patches(ii)%jl

        edgeLoop0: do iEdge=1,4
           select case(iEdge)
           case (iLow) 
              iSize = jl
              lPtr1 => patches(ii)%l_index(1, jl:1:-1)
              xPtr => patches(ii)%x(:, 1, jl:1:-1)
           case (iHigh)
              iSize = jl
              lPtr1 => patches(ii)%l_index(il, :)
              xPtr => patches(ii)%x(:, il, :)
           case (jLow) 
              iSize = il
              lPtr1 => patches(ii)%l_index(:, 1)
              xPtr => patches(ii)%x(:, :, 1)
           case (jHigh)
              iSize = il
              lPtr1 => patches(ii)%l_index(il:1:-1, jl)
              xPtr => patches(ii)%x(:, il:1:-1, jl)
           end select
           
           ! Just check the n=2 node:
           if (iSize > 2) then 
              i = 2
              iNode = lPtr1(i)
              if (fullTopoType(iNode) == topoEdge .and. BCs(iEdge, ii) == BCDefault) then 
                 if (unattachedEdgesAreSymmetry) then 
                    ! Determine which symm type we need:
                    xSum(1) = abs(sum(xptr(1, :)))
                    xSum(2) = abs(sum(xptr(2, :)))
                    xSum(3) = abs(sum(xptr(3, :)))
                    ! take dimension with lowest val
                    i = minloc(xSum, 1)
                    if (i== 1) then 
                       BCs(iEdge, ii) = BCXsymm
                    else if(i==2) then 
                       BCs(iEdge, ii) = BCYSymm
                    else if(i==3) then
                       BCs(iEdge, ii) = BCZSymm
                    end if
                    
                 else
                    ! Otherwise it must be a splay
                    BCs(iEdge, ii) = BCSplay
                 end if
              end if
           end if
        end do edgeLoop0
     end do patchLoop0
 
     nAverage = 0
     patchLoop: do ii=1, nPatch
        il = patches(ii)%il
        jl = patches(ii)%jl

        edgeLoop: do iEdge=1,4
           select case(iEdge)
           case (iLow) 
              lPtr1 => patches(ii)%l_index(1, jl:1:-1)
              lPtr2 => patches(ii)%l_index(2, jl:1:-1)
              xPtr => patches(ii)%x(:, 1, jl:1:-1)
              iSize = jl
              faceStr = 'iLow '
           case (iHigh)
              lPtr1 => patches(ii)%l_index(il, :)
              lPtr2 => patches(ii)%l_index(il-1, :)
              xPtr => patches(ii)%x(:, il, :)
              iSize = jl
              faceStr = 'iHigh '
           case (jLow) 
              lPtr1 => patches(ii)%l_index(:, 1)
              lPtr2 => patches(ii)%l_index(:, 2)
              xPtr => patches(ii)%x(:, :, 1)
              iSize = il
              faceStr = 'jLow '
           case (jHigh)
              lPtr1 => patches(ii)%l_index(il:1:-1, jl)
              lPtr2 => patches(ii)%l_index(il:1:-1, jl-1)
              xPtr => patches(ii)%x(:, il:1:-1, jl)
              iSize = il
              faceStr = 'jHigh'
           end select
           
           edgeNodeLoop: do i=1, iSize
              iNode = lPtr1(i)
              checkTopoType: if (fullTopoType(iNode) == topoInternal) then 
                 ! We don't have anything to do for this, BUT if the
                 ! boundary condition is anything but BCDefault than
                 ! it means someone explictly specifid it which is an
                 ! error since this isn't *actually* a boundary condition
                 if (BCs(iEdge, ii) /= BCDefault) then 
101 format (a,a,a,I3,a)
                    print 101, 'ERROR: A boundary condition was specifed for ', &
                         trim(faceStr), ' patch on Block', ii, ' but it is not a physical boundary.'
                    stop
                 end if

              else if (fullTopoType(iNode) == topoEdge) then 
                 ! By definition , edge topology nodes must have 3 neighbours
                 fullnPtr(1, iNode) = 3
                    
                 if (i> 1 .and. i < iSize) then 
                    
                    ! Edge topology, and internal to an edge, we can
                    ! just read off the neighbours.
                    fullnPtr(2, iNode) = lPtr1(i+1)
                    fullnPtr(3, iNode) = lPtr2(i  )
                    fullnPtr(4, iNode) = lPtr1(i-1)
                    fullBCType(1, iNode) = BCs(iEdge, ii)
     
                    ! Set the possible boundary condition value
                    call setBCVal(fullBCType(1, iNode), xPtr(:, i), &
                         fullBCVal(1, :, iNode))
                 else 
                    ! This is an edge topology, but at the logical
                    ! corner of a block. If some other edge has
                    ! already set the BC for us, we can just skip,
                    ! since they did the work.
                    if (BCs(iedge, ii) /= bcDefault) then 
                       if (i == 1) then 
                          fullnPTr(2, iNode) = lPtr1(i+1)
                          fullnPTr(3, iNode) = lPtr2(i  )
                          call addMissing(nodeConn(1:3, iNode), fullnPtr(2, iNode), fullNPtr(3, iNode), &
                               fullNPtr(4, iNode))
                       else ! other end
                          fullnPTr(3, iNode) = lPtr2(i  )
                          fullnPTr(4, iNode) = lPtr1(i-1)
                          call addMissing(nodeConn(1:3, iNode), fullnPtr(4, iNode), fullNPtr(3, iNode), &
                               fullNPtr(2, iNode))
                       end if
                       fullBCType(1, iNode) = BCs(iEdge, ii)

                       ! Set the possible boundary condition value
                       call setBCVal(fullBCType(1, iNode), xPtr(:, i), &
                            fullBCVal(1, :, iNode))
                    end if
                 end if
            
                 ! Setting the cell pointer is common
                 fullcPtr(1, iNode) = 2
                 fullcPtr(2, iNode) = nte(ntePtr(iNode))
                 fullcPtr(3, iNode) = nte(ntePtr(iNode)+2)
                                  
              else if (fullTopoType(iNode) == topoCorner) then 
                 ! Do a sanity check: This should only occur if i==1 or i==iSize
                 if (i/= 1 .and. i /= iSize) then 
                    print *, 'ERROR: Unknown corner topology error detected.', i, isize
                    stop
                 end if
                 
                 ! By definition , corner topology nodes must have 2 neighbours
                 fullnPtr(1, iNode) = 2
                 
                 if (i == 1) then 
                    fullnPtr(2, iNode) = lPtr1(i+1)
                    fullnPtr(3, iNode) = lPtr2(i)
                    iBCToSet = 2
                 else ! Must be iSize
                
                    fullnPtr(2, iNode) = lPtr2(i)
                    fullnPtr(3, iNode) = lPtr1(i-1)
                    iBCToSet = 1
                 end if

                 ! Set BC
                 fullBCType(iBCToSet, iNode) = BCs(iEdge, ii)

                 ! Set the possible boundary condition value
                 call setBCVal(fullBCType(iBCToSet, iNode), xPtr(:, i), &
                      fullBCVal(iBCToSet, :, iNode))

                 ! Set the cell pointer for this node
                 fullcPtr(1, iNode) = 1
                 fullcPtr(2, iNode) = nte(ntePtr(iNode))
        
              end if checkTopoType
           end do edgeNodeLoop
        end do edgeLoop

        ! We now do a special check to see if corners have boundary
        ! conditions that would require the averaging procedure. This
        ! will happen if the two neighbour nodes are constrained to
        ! the same plane. For example, having two ySymm conditions.

        corners(1) = patches(ii)%l_index(1, 1)
        corners(2) = patches(ii)%l_index(il, 1)
        corners(3) = patches(ii)%l_index(il, jl)
        corners(4) = patches(ii)%l_index(1, jl)

        do iCorner=1, 4
           i = corners(iCorner) ! Node index
           nn(1) = fullnPtr(2, i) ! Neighbour 1 index
           nn(2) = fullnPtr(3, i) ! Neighbour 2 index
           
           if (fullTopoType(i) == topoCorner) then 
              ! It is an actual corner topology
              mask = 0
              vals = zero
              do jj=1, 2 ! Loop over the two neighbours and check
                         ! their BC Type
                 select case (fullBCType(1, nn(jj)))

                 case (BCXSymm, BCXConst)
                    ! This means the following: the 'jj' neighbour
                    ! constrains the 'X' direction with value given in
                    ! fullBCVal. The mask(jj, 1) says the 'jj'
                    ! neighbour is constrained in 'x', with the value
                    ! given by vals(jj, 1)
                    mask(jj, 1) = 1
                    vals(jj, 1) = fullBCVal(1, 1, nn(jj))

                 case (BCYSymm, BCYConst)
                    mask(jj, 2) = 1
                    vals(jj, 2) = fullBCVal(1, 1, nn(jj))

                 case (BCZSymm, BCZConst)
                    mask(jj, 3) = 1
                    vals(jj, 3) = fullBCVal(1, 1, nn(jj))

                 case (BCXYConst)
                    mask(jj, 1) = 1
                    mask(jj, 2) = 1
                    vals(jj, 1) = fullBCVal(1, 1, nn(jj))
                    vals(jj, 2) = fullBCVal(1, 2, nn(jj))

                 case(BCYZConst)
                    mask(jj, 2) = 1
                    mask(jj, 3) = 1
                    vals(jj, 2) = fullBCVal(1, 1, nn(jj))
                    vals(jj, 3) = fullBCVal(1, 2, nn(jj))

                 case (BCXZConst)
                    mask(jj, 1) = 1
                    mask(jj, 3) = 1
                    vals(jj, 1) = fullBCVal(1, 1, nn(jj))
                    vals(jj, 3) = fullBCVal(1, 2, nn(jj))
                    
                 end select
              end do
              
              ! Loop over each coordinate direction. 
              do jj=1,3
                 ! Is the same coordinate constrained in each neighbour
                 if (mask(1, jj) == mask(2, jj) .and. mask(1, jj) == 1) then 

                    ! The same coordinate is constrained. Now check if
                    ! the value is the essentially the same
                    if (abs(vals(1, jj) - vals(2, jj)) < nodeTol) then 
                       fullBCType(:, i) = BCAverage
                       nAverage = nAverage + 1
                    end if
                 end if
              end do
           end if
        end do
     end do patchLoop

     do i=1,nUnique
        ! Special loop for the LCorners
        if (fullTopoType(i) == topoLCorner) then 
           
           ! This is quite complicated. It has 4 neighbours. Two
           ! of the neighbours have edge topology, two of the
           ! neighbours have internal topology. 
           !
           !     +---------p2----------+
           !     |         |           |
           !     |   x     |           |
           !     |         |           |
           !     p3--------p0----------p1----->
           !     |         |            
           !     |         |            
           !     |         |            
           !     +---------p4
           !               |
           !               |
           !              \|/
           
           
           ! Determine the two that have internal
           ! topology. These must be part of facex
           ii = 0
           jj = 0
           do j=1,4
              if (fullTopoType(nodeConn(j, i)) == topoInternal) then 
                 ii =ii + 1
                 nn(ii) = nodeConn(j, i) ! Nodes on diagonal face x
              else 
                 jj = jj + 1
                 mm(jj) = nodeConn(j, i) ! Nodes on edge
              end if
              
           end do
         
           ! And the three faces 
           f(1) = nte(ntePtr(i)  )
           f(2) = nte(ntePtr(i)+2)
           f(3) = nte(ntePtr(i)+4)
           
           do iFace=1,3
              ! Check if nn(1) is followed by nn(2)
              do j=1,4
                 jp2 = mod(j+1, 4)+1

                 if (nn(1) == fullConn(j, f(iFace)) .and. nn(2) == fullConn(jp2, f(iFace))) then 
                    fullnPtr(3, i) = nn(1)
                    fullnPtr(4, i) = nn(2)
                 else if (nn(2) == fullConn(j, f(iFace)) .and. nn(1) == fullConn(jp2, f(iFace))) then 
                    fullnPtr(3, i) = nn(2)
                    fullnPtr(4, i) = nn(1)
                 end if
              end do
           end do
           
           ! We have now set the 'p2' and 'p3' nodes. We have to
           ! identify p1 and p4 which should be fairly easy now:
           
           ! This takes the missing one and adds it to the first spot
           call addMissing(directedNodeConn(1:3, i), fullnPtr(3, i), fullnPtr(4, i), fullnPtr(2, i))
           
           ! Add the only other node in mm
           if (mm(1) == fullnPtr(2, i)) then
              fullnPtr(5, i) = mm(2)
           else
              fullnPtr(5, i) = mm(1)
           end if
           
           ! By definition , Lcorner topology nodes must have 4 neighbours
           fullnPtr(1, i) = 4

           ! Set the cell pointer for this node
           fullcPtr(1, i) = 3
           fullcPtr(2, i) = f(1)
           fullcPtr(3, i) = f(2)
           fullcPtr(4, i) = f(3)
           
        end if
     end do

     if (nAverage > 0) then 
        print *, 'Coner averaging activated for ', nAverage, ' nodes.'
     end if

     ! Do a sanity check to make sure that all nodes have their 
     ! fullnPtr populated. 
     do i=1,nUnique
        if (fullNPtr(1, i) == 0) then 
           print *, 'There was a general error with topology computation for node:', uniquePts(:, i)
           stop
        end if
     end do

     ! Another sanity check to make sure all nodes that are physically
     ! at edges have the required number of boundary conditions
     ! their fullnPtr populated.
     do i=1, nUnique
        if (fullTopoType(i) == topoEdge .and. fullBCType(1, i) == BCDefault) then 
           print *, 'There was a boundary condition for edge node:',  uniquePts(:, i)
           stop

        else if (fullTopoType(i) == topoCorner .and. (fullBCType(1, i) == BCDefault .or. &
             fullBCType(2, i) == BCDefault)) then 
           print *, 'There was a boundary condition for corner node:',  uniquePts(:, i)
           stop
        end if

     end do

     ! Free up some more memory
     deallocate(nte, ntePtr, directedNodeConn, nodeConn, nEdge, nDirectedEdge)
  
  end if rootProc

  ! Now we can broadcast the required infomation to everyone. In fact
  ! all we need is the fullnPtr, fullcPtr and the full connectivity list
  call MPI_Bcast(nUnique, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(NODE_MAX, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(CELL_MAX, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(faceTotal, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(nPatch, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Allocate space on all procs EXCEPT the root:
  if (myid /= 0) then
     allocate(fullConn(4, faceTotal))
     allocate(fullnPtr(NODE_MAX+1, nUnique))
     allocate(fullcPtr(CELL_MAX+1, nUnique))
     allocate(fullTopoType(nUnique))
     allocate(fullBCType(2, nUnique))
     allocate(fullBCVal(2, 2, nUnique))
     allocate(uniquePts(3, nUnique))
  end if

  ! Actual broadcasts
  call MPI_Bcast(fullConn, 4*faceTotal, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(fullnPtr, (NODE_MAX+1)*nUnique, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(fullcPtr, (CELL_MAX+1)*nUnique, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(uniquePts, 3*nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(fullTopoType, nUnique, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(fullBCType, 2*nUnique, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(fullBCVal, 2*2*nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! We now do "partitioning". It is essentially tries to divide up the
  ! nodes as evenly as possible (which is always possible up to the
  ! number of procs) without any regard to the number of cuts or
  ! communication costs. 

  evenNodes = nUnique / nProc
  istart = myid*evenNodes + 1
  iend = istart + evenNodes-1
  if (myid+1 == nProc) then
     iend = nUnique
  end if
  isize = iend - istart + 1

  ! Allocate the final space for gnPtr (global node ptr) and lnPtr
  ! (local node ptr) since we now know the range each proc owns

  allocate(lnPtr(NODE_MAX+1, isize), gnPtr(NODE_MAX+1, isize), &
       cPtr(CELL_MAX+1, iSize), topoType(iSize), BCType(2, iSize), &
       BCVal(2, 2, iSize))

  ! Copy in just the required parts:
  do i=istart, iend
     gnPtr(:, i-istart+1) = fullnPtr(:, i)
     lnPtr(:, i-istart+1) = fullnPtr(:, i) ! lnPtr will be modified below
     cPtr(:, i-istart+1) = fullcPtr(:, i) ! cPtr will be modified below
     topoType(i-iStart+1) = fullTopoType(i)
     BCType(:, i-iStart+1) = fullBCType(:, i)
     BCVal(:, :, i-iStart+1) = fullBCVal(:, :, i)
  end do

  ! Compute Ghost Information:

  ! Step 1: Count up the number of ghost nodes: We want all the nodes
  ! on all the elements surrounding the owned nodes. This will be
  ! sufficient to get the faces for each nodes. Note that not *all* of
  ! these ghost nodes may be required for lnPr (diagonals are not
  ! required for lnPtr, but are for the faces)
  nGhost = 0
  allocate(link(nUnique)) ! This will store the local ghost index for
  ! a non-owned node
  link(:) = 0

  do i=istart, iend
     do j=1, fullcPtr(1, i)
        iCell = fullcPtr(j+1, i)
        do k=1, 4
           ind = fullConn(k, iCell)
           
           ghostNode: if (ind < istart .or. ind > iend) then

              notAlreadyUsed: if (link(ind) == 0) then
                 nGhost = nGhost + 1
                 link(ind) = 1
              end if notAlreadyUsed
           end if ghostNode
        end do
     end do
  end do

  allocate(ghost(nGhost))
  ii = 0
  do i=1, nUnique
     if (link(i) /= 0) then
        ii = ii + 1
        ghost(ii) = i-1 !ghost is 0 based! (hence the -1)
        link(i) = ii
     end if
  end do
  
  do i=1, isize
     do j=1, lnPtr(1, i)
        ind = lnPtr(j+1, i)
        if (ind < istart .or. ind > iend) then
           lnPtr(j+1, i) = link(ind) + iSize ! Use ghost value from link
        else
           lnPtr(j+1, i) = lnPtr(j+1, i) - istart + 1 ! Local node
        end if
     end do
  end do

  ! Determine now many local faces we will have, (faces surrounding owned nodes)
  allocate(localFace(faceTotal))
  localFace(:) = 0
  do i=istart, iend
     do j=1, fullcPtr(1, i)
        localFace(fullcPtr(j+1, i)) = 1
     end do
  end do

  ! Now reorder the local faces
  nLocalFace = 0
  do i=1, faceTotal
     if (localFace(i) /= 0) then 
        nLocalFace = nLocalFace + 1
        localFace(i) = nLocalFace
     end if
  end do

  ! Fix cPtr to use the (compacted) local list of faces
  do i=1, isize
     do j=1, cPtr(1, i)
        cPtr(j+1, i) = localFace(cPtr(j+1, i))
     end do
  end do

  ! Finally we need to adjust the conn
  allocate(conn(4, nLocalFace))
  conn = -1
  ii = 0
  do i=1, faceTotal
     if (localFace(i) /= 0) then
        ii = ii + 1
        do j=1, 4
           ind = fullConn(j, i)
           if (ind < istart .or. ind > iend) then
              conn(j, ii) = link(ind) + iSize
           else
              conn(j, ii) = ind - istart + 1
           end if
        end do
     end if
  end do

  ! Set remainder of required information
  nx = isize
  nxglobal = nUnique

  call create3DPetscVars

  ! Copy X-surf into the first grid slot
  call VecGetArrayF90(X(1), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do i=1, nx
     do idim=1, 3
        xx((i-1)*3 + idim) = uniquePts(idim, i+istart-1)
     end do
  end do

  call VecRestoreArrayF90(X(1), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Update ghost values for the quality calc (this will be the -1 level)
  call VecGhostUpdateBegin(X(1), INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecGhostUpdateEnd(X(1), INSERT_VALUES,SCATTER_FORWARD, ierr)

  ! Finally finished with the full set of information so deallocate
  deallocate(fullnPtr, link, localFace, uniquePts, fullTopoType, fullBCType, &
       fullBCVal, fullconn, fullcptr)

end subroutine setup

subroutine getNBlocks(fileName, fileType, nBlocks)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: getNBlocks opens the input file and just reads the
  !     number of blocks and returns them. 
  !
  !     Description of Arguments
  !     Input:
  !     fileName : char* array : Filename of plot3d file. Usually has .fmt extension. 
  !     fileType : integer : 1-CGNS 2-Plot3d
  !     Output:
  !     nBlocks  : number of blocks

  use communication
  use hypData
  use hypInput

  implicit none
  include 'cgnslib_f.h'

  ! Input/Output Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: fileType
  integer(kind=intType), intent(out) :: nBlocks

  ! Working parameters
  integer(kind=intType) :: cg, ierr, base

  fType: if (fileType ==  cgnsFileType) then 

    ! Open and get the number of zones:
     call cg_open_f(trim(fileName), CG_MODE_READ, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
     base = 1_intType
     
     call cg_nzones_f(cg, base, nBlocks, ierr);
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

  else if (fileType .eq. plot3dFileType) then 
     
     ! Open file and read number of patches
     open(unit=7, form='formatted', file=fileName)
     read(7, *) nBlocks
     close(7)
  end if fType

  ! Allocate the families array since we know the number of patches
  allocate(families(nBlocks))

end subroutine getNBlocks

subroutine setFamily(i, fam)

  ! Helper routine to set a family string from python
  use hypInput
  implicit none
  
  character*(*) , intent(in) :: fam
  integer(kind=intType), intent(in) :: i

  families(i) = trim(fam)

end subroutine setFamily

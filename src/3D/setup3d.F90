subroutine setup(fileName)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: setup3d is the main routine for setting up the 3D
  !     elliptic or hyperbolic mesh extrusion problem. The only input
  !     is the filename of an ascii-formatted plot3d file and
  !     information about a possible mirror. 
  !    
  !     Description of Arguments
  !     Input:
  !     fileName : char* array : Filename of plot3d file. Usually has .fmt extension. 
  !     mirror : integer : Specifies if body must be mirrored before extrusion. 

  use communication
  use hypData
  use hypInput
  implicit none

  ! Input Arguments
  character*(*), intent(in) :: fileName

  ! Working parameters
  real(kind=realType), dimension(:, :), allocatable :: uniquePts, allNodes
  real(kind=realType), dimension(3) :: mask, tmp
  real(kind=realType) :: tol2

  integer(kind=intType), dimension(:), allocatable :: link, ptInd, nFace, localFace
  integer(kind=intType), dimension(:), allocatable :: nte, ntePtr
  integer(kind=intType), dimension(:, :), allocatable :: fullnPtr
  integer(kind=intType), dimension(:, :), allocatable :: patchSizes
  integer(kind=intType), dimension(:, :), allocatable :: nodeConn

  integer(kind=intType) :: nBlocks, nUnique
  integer(kind=intType) :: nodeTotal, node1, node2, node3, node4
  integer(kind=intType) :: nodeCounter, iNode, evenNodes
  integer(kind=intType) :: NODE_MAX, CELL_MAX
  integer(kind=intType) :: nFaceNeighbours, firstID, nodeID, nextElemID
  integer(kind=intType) :: nodeToAdd, faceToCheck, nodeToFind
  integer(kind=intType) :: i, j, k, ii, jj, iii, jjj, ierr, iStart, iEnd
  integer(kind=intType) :: isize, ind, icell
  logical :: found
  integer(kind=intType) :: nNode, nElem, shp2(2), shp1(1), idim
  interface 
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none
       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduce

     subroutine assignNode2NodeConn(nodeConn, node1, node2)
       use precision
       implicit none
       integer(kind=intType), intent(in) :: node1, node2
       integer(kind=intType), dimension(:,:), intent(inout) :: nodeConn
     end subroutine assignNode2NodeConn

  end interface
  ! Only the root processor reads the mesh and does the initial processing:
  rootProc: if (myid == 0) then 

     ! Open file and read number of patches
     open(unit=7, form='formatted', file=fileName)
     read(7, *) nBlocks

     ! Allocate space for the patch sizes and patches
     if (mirrorType == noMirror) then 
        nPatch = nBlocks
     else
        nPatch = nBlocks*2
     end if

     allocate(patchSizes(3, nPatch), patches(nPatch))

     read(7,*) (patchSizes(1, i), patchSizes(2, i), patchSizes(3, i), i=1, nBlocks)

     ! Make sure that ALL k-index patch sizes are 1
     do ii=1, nBlocks
        if (patchSizes(3, ii) /= 1) then
           print *,'Plot3d Read Error: k-dimension of all blocks must be precisely 1'
           stop
        end if
     end do

     ! Now allocate and read all the blocks from the plot3d file
     do  ii=1, nBlocks
        patches(ii)%il = patchSizes(1, ii)
        patches(ii)%jl = patchSizes(2, ii)

        ! Allocate space for the grid coordinates on the patch and read
        allocate(patches(ii)%X(3, patches(ii)%il, patches(ii)%jl))
        allocate(patches(ii)%l_index(patches(ii)%il, patches(ii)%jl))
        allocate(patches(ii)%weights(patches(ii)%il, patches(ii)%jl))
        patches(ii)%weights(:, :) = zero

        read(7, *) (( patches(ii)%X(1, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
        read(7, *) (( patches(ii)%X(2, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
        read(7, *) (( patches(ii)%X(3, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
     enddo
     deallocate(patchSizes)

     ! All the reading is done
     close(7)

     ! Need to do more work if we are mirroring
     mirrorCheck: if (mirrorType /= noMirror) then

        ! Just check each node is within tolerance of the sym plane
        mask(:) = one
        if (mirrorType == xmirror) then
           mask(1) = zero
        else if (mirrorType == ymirror) then
           mask(2) = zero
        else if (mirrorType == zmirror) then
           mask(3) = zero
        end if
        if (symTol < zero ) then
           tol2 = nodeTol ** 2
        else
           tol2 = symTol ** 2
        end if
        do ii=1,nBlocks
           patches(ii)%nSym = 0
           ! First do a pass to count:
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 tmp = patches(ii)%X(:, i, j) - patches(ii)%X(:, i, j)*mask
                 if (dot_product(tmp, tmp) < tol2) then
                    patches(ii)%nSym = patches(ii)%nSym + 1
                    patches(ii)%X(:, i, j) = patches(ii)%X(:, i, j)*mask
                 end if
              end do
           end do

           ! Check if the symmetry plane applies along the entire edge of the patch
           if ((patches(ii)%nSym .ne. 0) .and. &
                (patches(ii)%nSym .ne. patches(ii)%il) .and. & 
                (patches(ii)%nSym .ne. patches(ii)%jl)) then
              print *,'WARNING: Possible symmetry nodes extrapolate symmetry plane tolerance'
              print *,'         Increase nodeTol option or check if nodes are too far from'
              print *,'         the symmetry plane. Check the block with following coordinates:'
              print *,patches(ii)%X(:, 1, 1)
              print *,patches(ii)%X(:, 1, patches(ii)%jl)
              print *,patches(ii)%X(:, patches(ii)%il, patches(ii)%jl)
              print *,patches(ii)%X(:, patches(ii)%il, 1)
           end if

           ! Now allocate and store
           allocate(patches(ii)%symNodes(2, patches(ii)%nSym))
           iNode = 0
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 tmp = patches(ii)%X(:, i, j) - patches(ii)%X(:, i, j)*mask
                 if (dot_product(tmp, tmp) < tol2) then
                    iNode = iNode + 1
                    patches(ii)%symNodes(:, iNode) = (/i, j/)
                 end if
              end do
           end do
        end do

        ! Now we can finish mirroring the mesh:
        mask(:) = one
        if (mirrorType == xmirror) then
           mask(1) = -one
        else if (mirrorType == ymirror) then
           mask(2) = -one
        else if (mirrorType == zmirror) then
           mask(3) = -one
        end if

        do ii=1,nBlocks
           allocate(patches(ii+nBlocks)%X(3, patches(ii)%il, patches(ii)%jl))
           allocate(patches(ii+nBlocks)%l_index(patches(ii)%il, patches(ii)%jl))
           allocate(patches(ii+nBlocks)%weights(patches(ii)%il, patches(ii)%jl))
           patches(ii+nBlocks)%weights(:, :) = zero

           patches(ii+nBlocks)%il = patches(ii)%il
           patches(ii+nBlocks)%jl = patches(ii)%jl
           do j=1,patches(ii)%jl
              do i=1,patches(ii)%il
                 patches(ii+nBlocks)%X(:, i, j) = &
                      patches(ii)%X(:, i, j)*mask
              end do
           end do
        end do
     end if mirrorCheck

     ! Now we can create the final connectivity
     nodeTotal = 0
     faceTotal = 0
     do ii=1,nPatch
        nodeTotal = nodeTotal + patches(ii)%il*patches(ii)%jl
        faceTotal = faceTotal + (patches(ii)%il-1)*(patches(ii)%jl-1)
     end do

     allocate(allNodes(3, nodeTotal), uniquePts(3, nodeTotal), link(nodeTotal))
     iNode = 0
     do ii=1, nPatch
        do j=1, patches(ii)%jl, patches(ii)%jl-1
           do i=1, patches(ii)%il
              iNode = iNode + 1
              allNodes(:, iNode) = patches(ii)%X(:, i, j)
           end do
        end do
        do i=1, patches(ii)%il, patches(ii)%il-1
           do j=2, patches(ii)%jl-1
              iNode = iNode + 1
              allNodes(:, iNode) = patches(ii)%X(:, i, j)
           end do
        end do
     end do

     call pointReduce(allNodes, iNode, nodeTol, uniquePts, link, nUnique)

     ! Dump edge/cornder into l_index
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
     iNode = nUnique
     do ii=1,nPatch
        do j=2, patches(ii)%jl-1
           do i=2, patches(ii)%il-1
              iNode = iNode + 1
              uniquePts(:, iNode) = patches(ii)%X(:, i, j)
              patches(ii)%l_index(i, j) = iNode
           end do
        end do
     end do
     nUnique = iNode

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
     deallocate(link, allNodes)

     ! Create the conn array
     allocate(fullConn(4, faceTotal))
     jj = 0
     do ii=1,nPatch
        if (mirrorType==noMirror .or. ii <= nPatch/2) then
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
        else
           ! Patches are flipped orientation
           do j=1, patches(ii)%jl-1
              do i=1, patches(ii)%il-1
                 jj = jj + 1
                 fullConn(:, jj) = (/&
                      patches(ii)%l_index(i,   j), &
                      patches(ii)%l_index(i, j+1), &
                      patches(ii)%l_index(i+1, j+1), &
                      patches(ii)%l_index(i+1, j)/)
              end do
           end do
        end if
     end do

     ! NORMAL ORIENTATION CHECK
     
     ! Print what we are doing
     print *,'Normal orientation check ...'

     ! Allocate and initialize matrix that stores nodal connectivity directions
     ! I'm assuming that the maximum number of connections of a single node
     ! is 10
     allocate(nodeConn(10, nodeTotal))
     nodeConn(:,:) = 0

     ! Loop over each face and increment elementes from the connectivity
     ! matrix acording to the connectivity direction.
     do jj=1,faceTotal

        ! Get nodes indices for this face
        node1 = fullConn(1, jj)
        node2 = fullConn(2, jj)
        node3 = fullConn(3, jj)
        node4 = fullConn(4, jj)

        ! Assign connections acording to their orientations
        ! The assignNode2NodeConn subroutine is defined in
        ! 3D_utilities.F90
        call assignNode2NodeConn(nodeConn, node1, node2)
        call assignNode2NodeConn(nodeConn, node2, node3)
        call assignNode2NodeConn(nodeConn, node3, node4)
        call assignNode2NodeConn(nodeConn, node4, node1)

     end do

     ! Now check if we flagged any node during the assignments
     if (minval(nodeConn) .eq. -1) then

        ! Say that the normals are wrong
        print *,' ERROR: Normal directions may be wrong'
        print *,'        Check the following coordinates:'

        ! Print coordinates of nodes were the problem occurs
        do i=1,nodeTotal
           if (nodeConn(1,i) .eq. -1) then
              print *,uniquePts(:, i)
           end if
        end do

        ! Halt program execution
        stop

     end if

     print *,'Normals seem correct!'

     ! END OF NORMAL ORIENTATION CHECK

     ! Determine the number of number of faces around each node:
     allocate(nFace(nUnique))
     nFace(:) = 0
     do i=1, faceTotal
        do ii=1, 4
           nFace(fullConn(ii, i)) = nFace(fullConn(ii, i)) + 1
        end do
     end do

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
        nFaceNeighbours = (ntePtr(iNode+1)-ntePtr(iNode))/2

        ! Regular nodes get 4 neighbours, others get double
        if (nFaceNeighbours == 2) then
           fullnPtr(1, iNode) = nFaceNeighbours*2
        else if (nFaceNeighbours == 3) then
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
     end do nodeLoop
     deallocate(nte, ntePtr)
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

  ! We now do "partitioning". It is essentially tries to diviy up the
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

  allocate(lnPtr(NODE_MAX+1, isize), gnPtr(NODE_MAX+1, isize), cPtr(CELL_MAX+1, iSize))

  ! Copy in just the required parts of the fullnPtr
  do i=istart, iend
     gnPtr(:, i-istart+1) = fullnPtr(:, i)
     lnPtr(:, i-istart+1) = fullnPtr(:, i) ! lnPtr will be modified below
     cPtr(:, i-istart+1) = fullcPtr(:, i) ! cPtr will be modified below
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
        iCell = fullcPtr(j+1, i) ! 
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

  do i=1,nx
     do idim=1,3
        xx((i-1)*3 + idim) = uniquePts(idim, i+istart-1)
     end do
  end do

  call VecRestoreArrayF90(X(1), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Update ghost values for the quality calc (this will be the -1 level)
  call VecGhostUpdateBegin(X(1), INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecGhostUpdateEnd(X(1), INSERT_VALUES,SCATTER_FORWARD, ierr)

  ! Finally finished with the full set of information so deallocate
  deallocate(fullnPtr, link, localFace, uniquePts)

end subroutine setup




subroutine surfaceSmooth(xVec, nSteps, stepSize)
  use hypData, only : patches, nPatch, cPtr, xx, nx, conn, nghost, nlocalFace, XL_local, X
  use hypInput, only : nodeTol
  use kd_tree
  use communication
  use petsc
  implicit none
#include "include/petscversion.h"
#include "petsc/finclude/petsc.h"

  ! Input Parameters
  integer (kind=intType), intent(in) :: nSteps
  real(kind=realType), intent(in) :: stepSize
  Vec xVec
  ! Working Parameters
  integer(kind=intType) :: i,j, ii, iCell, iter, ierr, nUnique, nFrozen, surfID
  real(kind=realType), dimension(3, nx) :: massCen
  real(kind=realType), dimension(nx) :: distFact
  real(kind=realType), dimension(3, nx+nghost) :: normals
  real(kind=realType), dimension(3) :: v1, v2, s, xdot, dx
  real(kind=realType), dimension(3, nLocalFace) :: cellCen
  real(kind=realType), dimension(nLocalFace) :: areas
  real(kind=realType), dimension(:, :), allocatable :: allNodes, uniqueNodes, tmpNodes
  real(kind=realType), dimension(:), allocatable :: weights
  integer(kind=intType), dimension(:), allocatable :: link
  real(kind=realType) :: aTot, sNorm, fact, dstar
  real(kind=realType), target :: distances(1)
  integer(kind=intType), target :: indexes(1)

  type(tree_master_record), pointer :: mytree

  ! First we have compute the distance factor. We need to get a local
  ! copy of all the nodes onto every proc so that everyone can create
  ! a KDTree with all the nodes


     ! First determine the number of frozen nodes we have...only on the root proc:
     if (myid == 0) then
        nFrozen = 0
        do ii=1, nPatch
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 if (patches(ii)%weights(i, j) /= zero) then
                    nFrozen = nFrozen + 1
                 end if
              end do
           end do
        end do

        ! Allocate space and loop back over to fill up
        allocate(allNodes(3, nFrozen), weights(nFrozen), tmpNodes(3, nFrozen), &
             link(nFrozen))
        nFrozen = 0
        do ii=1, nPatch
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 if (patches(ii)%weights(i, j) /= zero) then
                    nFrozen = nFrozen + 1
                    allNodes(:, nFrozen) = patches(ii)%X(:, i, j)
                 end if
              end do
           end do
        end do
        if (nFrozen > 0) then
           ! Now we will point reduce them to make sure everything is unique
           call pointReduce(allNodes, nFrozen, nodeTol, tmpNodes, link, nUnique)
           allocate(uniqueNodes(3, nUnique))
           do i=1, nUnique
              uniqueNodes(:, i) = tmpNodes(:, i)
           end do
           ! Lastly we use the link array to fill up (nUnique) points on weights
           nFrozen = 0
           do ii=1, nPatch
              do j=1, patches(ii)%jl
                 do i=1, patches(ii)%il
                    if (patches(ii)%weights(i, j) /= zero) then
                       nFrozen = nFrozen + 1
                       weights(link(nFrozen)) = patches(ii)%weights(i, j)
                    end if
                 end do
              end do
           end do
        else
           nUnique = 0
        end if
        deallocate(link, allNodes, tmpNodes)
     end if

     ! Get the unique size on all procs.
     call MPI_Bcast(nUnique, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if (nUnique > 0) then
        ! Note that uniqueNodes is length nFrozen on the root proc by
        ! nUnique on all the other procs.  There is no need to redo stuff
        ! on the root proc.
        if (myid /= 0) then
           allocate(uniqueNodes(3, nUnique), weights(nUnique))
        end if

        call MPI_Bcast(uniqueNodes, 3*nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call MPI_Bcast(weights, nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Now create the global kdTree which will be the same on all procs
        mytree => create_tree(uniqueNodes)

        ! Now we can just search our own locally owned nodes and compute
        ! the required factor.
        call VecGetArrayF90(XVec, xx, ierr)
        call EChk(ierr,__FILE__,__LINE__)
        do i=1, nx
           call n_nearest_to(mytree, xx(3*i-2:3*i), 1, indexes, distances)
           dstar = weights(indexes(1))
           distFact(i) = min((sqrt(distances(1))/dstar)**1.5, one)
        end do

        call VecRestoreArrayF90(XVec, xx, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call destroy_tree(mytree)
        deallocate(uniqueNodes, weights)
     else
        distFact(:) = one
     end if

  masterIteration: do iter=1,nSteps

     ! Update the ghost entries before we start
     call VecGhostUpdateBegin(Xvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGhostUpdateEnd(XVec, INSERT_VALUES,SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get local form for doing stuff
     call VecGhostGetLocalForm(XVec, XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     cellCen(:, :) = zero
     massCen(:, :) = zero
     normals(:, :) = zero

     do i=1, nLocalFace
        ! Compute face center:
        do ii=1,4
           cellCen(:, i) = cellCen(:, i) + fourth*xx(3*conn(ii, i)-2:3*conn(ii, i))
        end do

        ! Compute the normal of this face
        v1 = xx(3*conn(3, i)-2:3*conn(3, i)) - xx(3*conn(1, i)-2:3*conn(1, i))
        v2 = xx(3*conn(4, i)-2:3*conn(4, i)) - xx(3*conn(2, i)-2:3*conn(2, i))

        ! Cross product
        s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
        s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
        s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

        ! Compute the area
        snorm = sqrt(s(1)**2 + s(2)**2 + s(3)**2)

        ! Area:
        areas(i) = half * snorm

        ! Now normalize:
        s = s / snorm

        ! Scatter to each node
        do ii=1,4
           normals(:, conn(ii, i)) = normals(:, conn(ii, i)) + s
        end do
     end do

     ! Now make the second pass over nodes and re-normalize. Also compute
     ! the mass center for each node
     do i=1, nx
        s = normals(:, i)
        normals(:, i) = s / sqrt(s(1)**2 + s(2)**2 + s(3)**2)

        aTot = zero
        ! Loop over each cell around the node
        do ii=1, cPtr(1, i)
           iCell = cPtr(ii+1, i)
           massCen(:, i) = massCen(:, i) + cellCen(:, iCell)*areas(iCell)
           aTot = aTot + areas(iCell)
        end do
        massCen(:, i) = massCen(:, i) / aTot
     end do

     ! Now perform the explict euler udpate
     do i=1, nx
        dx = massCen(:, i) - xx(3*i-2:3*i)
        xDot = dx - dot_product(normals(:, i), dx)*normals(:, i)
        xx(3*i-2:3*i) = xx(3*i-2:3*i) + xDot*stepSize*distFact(i)
     end do

     ! Restore everything to prepare for the next iteration
     call VecRestoreArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGhostRestoreLocalForm(XVec, XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do masterIteration
end subroutine surfaceSmooth

subroutine freezeEdge(blockID, edge, dstar)
  use hypInput
  use hypData
  use communication
  implicit none

  ! Input Parameters
  integer (kind=intType), intent(in) :: blockID
  character*(*) :: edge
  real(kind=realType), intent(in) :: dstar

  ! Working
  integer(kind=intType) :: il, jl, offset

  ! Only the root proc keeps track of the specified weights
  if (myid == 0) then
     il = patches(blockID)%il
     jl = patches(blockID)%jl

     if (trim(edge) == 'ilow') then
        patches(blockID)%weights(1, :) = dstar
     else if (trim(edge) == 'ihigh') then
        patches(blockID)%weights(il, :) = dstar
     else if (trim(edge) == 'jlow') then
        patches(blockID)%weights(:, 1) = dstar
     else if (trim(edge) == 'jhigh') then
        patches(blockID)%weights(:, jl) = dstar
     end if

     ! if (mirrorType /= noMirror) then
     !    offset = nPatch/2
     !    if (trim(edge) == 'ilow') then
     !       patches(blockID+offset)%weights(1, :) = dstar
     !    else if (trim(edge) == 'ihigh') then
     !       patches(blockID+offset)%weights(il, :) = dstar
     !    else if (trim(edge) == 'jlow') then
     !       patches(blockID+offset)%weights(:, 1) = dstar
     !    else if (trim(edge) == 'jhigh') then
     !       patches(blockID+offset)%weights(:, jl) = dstar
     ! end if
     !  end if


  end if

end subroutine freezeEdge

subroutine freezeFaces(blockIDs, nBlockIDs, dstar)
  use hypInput
  use hypData
  use communication
  implicit none

  ! Input Parameters
  integer (kind=intType), intent(in), dimension(nBlockIDs) :: blockIDs
  integer(kind=intType), intent(in) :: nBLockIDs
  real(kind=realType), intent(in) :: dstar

  ! Working
  integer(kind=intType) :: i

  ! Only the root proc keeps track of the specified weights
  if (myid == 0) then
     do i=1, nBlockIDs
        patches(blockIDs(i))%weights(:, :) = dstar
     end do

     ! if (mirrorType /= noMirror) then
     !    do i=1, nBlockIDs
     !       patches(blockIDs(i)+nPatch/2)%weights(:, :) = dstar
     !    end do
     ! end if
  end if
end subroutine freezeFaces


subroutine smoothWrap(nSteps, stepSize)
  use precision
  use hypData
  implicit none

  ! Input Parameters
  integer (kind=intType), intent(in) :: nSteps
  real(kind=realType), intent(in) :: stepSize

  call surfaceSmooth(X(1), nSteps, stepSize)

end subroutine smoothWrap

subroutine surfaceSmooth2(xVec, nSteps, stepSize)
  use hypData, only : patches, nPatch, cPtr, xx, nx, conn, nghost, nlocalFace, XL_local, X
  use hypInput, only : nodeTol
  use kd_tree
  use communication
  use petsc
  implicit none
#include "include/petscversion.h"
#include "petsc/finclude/petsc.h"

  ! Input Parameters
  integer (kind=intType), intent(in) :: nSteps
  real(kind=realType), intent(in) :: stepSize
  Vec xVec
  ! Working Parameters
  integer(kind=intType) :: i,j, ii, iCell, iter, ierr, nUnique, nFrozen, surfID
  real(kind=realType), dimension(3, nx) :: massCen
  real(kind=realType), dimension(nx) :: distFact
  real(kind=realType), dimension(3, nx+nghost) :: normals
  real(kind=realType), dimension(3) :: v1, v2, s, xdot, dx
  real(kind=realType), dimension(3, nLocalFace) :: cellCen
  real(kind=realType), dimension(nLocalFace) :: areas
  real(kind=realType), dimension(:, :), allocatable :: allNodes, uniqueNodes, tmpNodes
  real(kind=realType), dimension(:), allocatable :: weights
  integer(kind=intType), dimension(:), allocatable :: link
  real(kind=realType) :: aTot, sNorm, fact, dstar
  real(kind=realType), target :: distances(1)
  integer(kind=intType), target :: indexes(1)

  type(tree_master_record), pointer :: mytree
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
  end interface

  ! First we have compute the distance factor. We need to get a local
  ! copy of all the nodes onto every proc so that everyone can create
  ! a KDTree with all the nodes


     ! First determine the number of frozen nodes we have...only on the root proc:
     if (myid == 0) then
        nFrozen = 0
        do ii=1, nPatch
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 if (patches(ii)%weights(i, j) /= zero) then
                    nFrozen = nFrozen + 1
                 end if
              end do
           end do
        end do

        ! Allocate space and loop back over to fill up
        allocate(allNodes(3, nFrozen), weights(nFrozen), tmpNodes(3, nFrozen), &
             link(nFrozen))
        nFrozen = 0
        do ii=1, nPatch
           do j=1, patches(ii)%jl
              do i=1, patches(ii)%il
                 if (patches(ii)%weights(i, j) /= zero) then
                    nFrozen = nFrozen + 1
                    allNodes(:, nFrozen) = patches(ii)%X(:, i, j)
                 end if
              end do
           end do
        end do
        if (nFrozen > 0) then
           ! Now we will point reduce them to make sure everything is unique
           call pointReduce(allNodes, nFrozen, nodeTol, tmpNodes, link, nUnique)
           allocate(uniqueNodes(3, nUnique))
           do i=1, nUnique
              uniqueNodes(:, i) = tmpNodes(:, i)
           end do
           ! Lastly we use the link array to fill up (nUnique) points on weights
           nFrozen = 0
           do ii=1, nPatch
              do j=1, patches(ii)%jl
                 do i=1, patches(ii)%il
                    if (patches(ii)%weights(i, j) /= zero) then
                       nFrozen = nFrozen + 1
                       weights(link(nFrozen)) = patches(ii)%weights(i, j)
                    end if
                 end do
              end do
           end do
        else
           nUnique = 0
        end if
        deallocate(link, allNodes, tmpNodes)
     end if

     ! Get the unique size on all procs.
     call MPI_Bcast(nUnique, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if (nUnique > 0) then
        ! Note that uniqueNodes is length nFrozen on the root proc by
        ! nUnique on all the other procs.  There is no need to redo stuff
        ! on the root proc.
        if (myid /= 0) then
           allocate(uniqueNodes(3, nUnique), weights(nUnique))
        end if

        call MPI_Bcast(uniqueNodes, 3*nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call MPI_Bcast(weights, nUnique, MPI_DOUBLE, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Now create the global kdTree which will be the same on all procs
        mytree => create_tree(uniqueNodes)

        ! Now we can just search our own locally owned nodes and compute
        ! the required factor.
        call VecGetArrayF90(XVec, xx, ierr)
        call EChk(ierr,__FILE__,__LINE__)
        do i=1, nx
           call n_nearest_to(mytree, xx(3*i-2:3*i), 1, indexes, distances)
           dstar = weights(indexes(1))
           distFact(i) = min((sqrt(distances(1))/dstar)**1.5, one)
        end do

        call VecRestoreArrayF90(XVec, xx, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call destroy_tree(mytree)
        deallocate(uniqueNodes, weights)
     else
        distFact(:) = one
     end if

  masterIteration: do iter=1,nSteps

     ! Update the ghost entries before we start
     call VecGhostUpdateBegin(Xvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGhostUpdateEnd(XVec, INSERT_VALUES,SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get local form for doing stuff
     call VecGhostGetLocalForm(XVec, XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     cellCen(:, :) = zero
     massCen(:, :) = zero
     normals(:, :) = zero

     do i=1, nLocalFace
        ! Compute face center:
        do ii=1,4
           cellCen(:, i) = cellCen(:, i) + fourth*xx(3*conn(ii, i)-2:3*conn(ii, i))
        end do

        ! Compute the normal of this face
        v1 = xx(3*conn(3, i)-2:3*conn(3, i)) - xx(3*conn(1, i)-2:3*conn(1, i))
        v2 = xx(3*conn(4, i)-2:3*conn(4, i)) - xx(3*conn(2, i)-2:3*conn(2, i))

        ! Cross product
        s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
        s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
        s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

        ! Compute the area
        snorm = sqrt(s(1)**2 + s(2)**2 + s(3)**2)

        ! Area:
        areas(i) = half * snorm

        ! Now normalize:
        s = s / snorm

        ! Scatter to each node
        do ii=1,4
           normals(:, conn(ii, i)) = normals(:, conn(ii, i)) + s
        end do
     end do

     ! Now make the second pass over nodes and re-normalize. Also compute
     ! the mass center for each node
     do i=1, nx
        s = normals(:, i)
        normals(:, i) = s / sqrt(s(1)**2 + s(2)**2 + s(3)**2)

        aTot = zero
        ! Loop over each cell around the node
        do ii=1, cPtr(1, i)
           iCell = cPtr(ii+1, i)
           massCen(:, i) = massCen(:, i) + cellCen(:, iCell)*areas(iCell)
           aTot = aTot + areas(iCell)
        end do
        massCen(:, i) = massCen(:, i) / aTot
     end do

     ! Now perform the explict euler udpate
     do i=1, nx
        dx = massCen(:, i) - xx(3*i-2:3*i)
        xDot = dx !- dot_product(normals(:, i), dx)*normals(:, i)
        xx(3*i-2:3*i) = xx(3*i-2:3*i) + xDot*stepSize*distFact(i)
     end do

     ! Restore everything to prepare for the next iteration
     call VecRestoreArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGhostRestoreLocalForm(XVec, XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do masterIteration
end subroutine surfaceSmooth2

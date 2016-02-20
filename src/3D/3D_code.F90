subroutine runHyperbolic
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: run3D is the main python interface to generate the
  !     3D hyperbolic mesh.
  !    
  use communication
  use hypData
  use hypInput
  implicit none

  ! Working parameters
  integer(kind=intType) :: i, ierr, idim, l, reason
  real(kind=realType) :: tmp
  ! Compute the inital 'Radius' value, radius0
  call computeR(X(1), radius0)

  ! Now the user has specified how far out to march in terms of
  ! multiplies of radius0. This gives us a nominal distance. Since we know
  ! the length of the off-wall direction and the initial spacing, we
  ! find the gridRatio that will satifiy these two conditions.
  call calcGridRatio(N, s0, radius0*rMin, gridRatio)

  ! Write header (defined in 3D_utilities.f90)
  if (myid == 0) then
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") "Grid Ratio:"
     write(*,"(f8.4,1x)",advance="no") gridRatio
     print "(1x)" 
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"
     
     call writeHeader
  end if

  ! Set initial pseudo and real spacings
  deltaS = ps0

  ! Zero out the cumulative distance counter
  scaleDist = zero
  nSubIterPrev = 0
  desiredS = zero
  ! Determine starting CPU Time
  timeStart = dble(mpi_wtime())

  ! This is the master marching direction loop
  marchLoop: do marchIter=2, N

     ! Use 'L' as the working variable
     L = marchIter
     ! Run the "initial guess" function. If the user is running in
     ! 'linear' mode this is all that is done. This function computes
     ! the next step, X(L). It may take multiple sub-steps to get there
     ! which is fine. 
     call initialGuess(X(L))
     
     ! Compute the max radius of the current level, X(L)
     call computeMinR(X(L), Radius)

     ! Check the quality of this layer
     call computeQualityLayer

     ! Compute the stretching
     if (L .gt. 2) then
        call findKStretch(X(L), X(L-1), X(L-2))
     else
        maxKStretch = 0.0
     end if

     ! Possibly write header and write iteration info
     if (myid == 0) then
        if (mod(marchIter, 50) == 0) then
           call writeHeader
        end if
     end if

     ! Write info for this layer
     call writeIteration

  end do marchLoop

end subroutine runHyperbolic

subroutine computeVolumes
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: computeVolumes() determines the nodal volumes for
  !     the owned nodes. An "exact" approach is employed; The area of
  !     each quad element is found and scattered to each of the four
  !     surrounding noes. Also computes the average nodal volume, Vbar.
  !
  !     Ney Secco (2015-2016): Added extrapolation factors
  !
  !     Description of Arguments
  !     Input:
  !     xx (in hypData), array size(3*nx) - Vector of nodes to use
  !
  !     Ouput:
  !     volume (in hypData), array size(3*nx) - Vector of volumes, stored in the third DOF
  !     vBar (in hypData), real : Average nodal volume

  use communication
  use precision
  use hypData
  implicit none

  ! Working Variables
  integer(kind=intType) :: i, j, ipatch, ierr, n1, n2, n3, n4 
  integer(kind=intType) :: nodeArray(4), nodeID, iBC, nPointer, side, nodeType(4)
  real(kind=realType) :: ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS, factor
  real(kind=realType) :: c1,c2,c3,c4, dist
  real(kind=realType) :: Vbar_local, cratio_local
  logical :: search

  call VecZeroEntries(volume, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGhostGetLocalForm(Volume, VolumeLocal, ierr)
  call VecGhostGetLocalForm(XL, XL_local, ierr)
  call VecGetArrayF90(VolumeLocal, Vptr, ierr)
  call VecGetArrayF90(XL_local, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  cRatio_local = zero
  vBar_local = zero

  ! Loop over the owned faces:
  do i=1, nLocalFace

     ! Extract the coordinates of the 4 corners
     n1 = conn(1, i)
     ll = xx(3*n1-2:3*n1) ! lower-left

     n2 = conn(4, i)
     ul = xx(3*n2-2:3*n2) ! upper-left

     n3 = conn(2, i)
     lr = xx(3*n3-2:3*n3) ! lower-right

     n4 = conn(3, i)
     ur = xx(3*n4-2:3*n4) ! upper-right

     ! Compute the area
     v1(:) = ur - ll
     v2(:) = ul - lr

     ! Cross Product
     s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

     areadS = fourth*deltaS*half*sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))

     ! Now we need to check if the nodes are corners or edges
     nodeArray = (/n1, n2, n3, n4/)

     ! Initialize node types
     nodeType = (/0, 0, 0, 0/)

     ! Classify nodes
     ! 0 = interior
     ! 1 = edge
     ! 2 = corner

     do nodeID = 1,4
        ! Get node index
        nPointer = nodeArray(nodeID)
        ! Check if we got a ghost node (index beyond size of BCneighbors)
        ! nX is the number of local nodes
        if (nPointer .le. nX) then
           ! Check if we have BCs
           if (BCneighborsLocal(nPointer,1) .gt. 0) then
              ! Check the number of neighbor faces
              if (cPtr(1,nPointer) .eq. 1) then ! We have a corner
                 ! Assign corner type
                 nodeType(nodeID) = 2
                 ! Stop search loop
                 search = .false.
              else ! Point is edge
                 ! Assign edge type
                 nodeType(nodeID) = 1
                 ! Stop search loop
                 search = .false.
              end if
           end if
        else
           ! We have no BC info on the ghost node so we will take it as interior node
           nodeType(nodeID) = 0
        end if
     end do

     ! Now we can compute area factors
     do nodeID = 1,4
        ! Get node index
        nPointer = nodeArray(nodeID)
        ! For interior nodes, no correction is needed
        factor = one
        ! Identify if we have an edge or corner node
        if (nodeType(nodeID) .eq. 2) then ! We have a corner
           ! Find side that need extrapolation
           if (nodeID .eq. 1) then
              side = 8
           else if (nodeID .eq. 2) then
              side = 5
           else if (nodeID .eq. 3) then
              side = 7
           else if (nodeID .eq. 4) then
              side = 6
           end if
           ! We will consider the extra 3 halo faces to set factor
           call extrapolateArea(ul,ur,lr,ll,side,factor)
        else if (nodeType(nodeID) .eq. 1) then ! Point is edge
           ! Find side that need extrapolation
           if (nodeID .eq. 1) then
              if (nodeType(2) .ne. 0) then !Edge between ll and ul
                 side = 1
              else if (nodeType(3) .ne. 0) then  !Edge between ll and lr
                 side = 4
              end if
           else if (nodeID .eq. 2) then
              if (nodeType(1) .ne. 0) then  !Edge between ul and ll
                 side = 1
              else if (nodeType(4) .ne. 0) then  !Edge between ul and ur
                 side = 2
              end if
           else if (nodeID .eq. 3) then
              if (nodeType(1) .ne. 0) then  !Edge between lr and ll
                 side = 4
              else if (nodeType(4) .ne. 0) then  !Edge between lr and ur
                 side = 3
              end if
           else if (nodeID .eq. 4) then
              if (nodeType(2) .ne. 0) then  !Edge between ur and ul
                 side = 2
              else if (nodeType(3) .ne. 0) then  !Edge between ur and lr
                 side = 3
              end if
           end if
           ! We will consider the extra 3 halo faces to set factor
           call extrapolateArea(ul,ur,lr,ll,side,factor)
        end if
        vPtr(nPointer) = vPtr(nPointer) + factor*areadS
        vBar_local = vBar_local + factor*areadS
     end do

     ! Also compute edge lengths
     c1 = deltaS/dist(ll,lr)
     c2 = deltaS/dist(ul,ur)
     c3 = deltaS/dist(ll,ul)
     c4 = deltaS/dist(lr,ur)

     cRatio_local = max(cRatio_local, c1, c2, c2, c4)
  end do

  ! Restore everything
  call VecRestoreArrayF90(VolumeLocal, Vptr, ierr)
  call VecRestoreArrayF90(XL_local, xx, ierr)
  call VecGhostRestoreLocalForm(Volume, VolumeLocal, ierr)
  call VecGhostRestoreLocalForm(XL, XL_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MPI_allreduce(vBar_local, vBar, 1, MPI_DOUBLE, MPI_SUM, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_allreduce(cRatio_local, cRatio, 1, MPI_DOUBLE, MPI_MAX, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  VBar = VBar / nXGlobal

  call VecGhostUpdateBegin(volume, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecGhostUpdateEnd(volume, INSERT_VALUES,SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
end subroutine computeVolumes

subroutine volumeSmooth
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: volumeSmooth() performs two typs of volume smoothing
  !     on the list of nodal volumes, volume. First it performs a fix
  !     number of point jacobi iterations. Then it performs a
  !     global volume blending based on vBar.
  !
  !     Description of Arguments
  !     Input:
  !     volume (in hypData) - pointer to nodal volumes. New volumes are
  !            in this vector on output.
  !     vBar (in hypData) : Average nodal volume
  !
  !     Ouput:
  !     volume - See above. 

  use communication
  use hypInput
  use hypData
  implicit none

  ! Working Parameters
  real(kind=realType) :: factor, nNeighbors, vSum
  integer(kind=intType) :: i, iSize, iter, ierr, ii,jj
  real(kind=realType), allocatable, dimension(:) :: Vtmp
  real(kind=realType) :: timea, timeb
  logical :: search

  ! Do a Jacobi volume smooth

  call VecGhostGetLocalForm(Volume, VolumeLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetSize(VolumeLocal, iSize, ierr)
  allocate(vTmp(iSize))

  call VecGetArrayF90(VolumeLocal, Vptr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do iter=1,volSmoothIter

     ! Copy vptr to vtmp
     do i=1,isize
        vtmp(i) = vptr(i)
     end do

     do i=1,nx

        vSum = zero
        nNeighbors = lnPtr(1,i)

        do ii=1,lnPtr(1,i)
           ! Sum physical  neighbors
           vSum = vSum + Vtmp(lnPtr(1+ii,i))
        end do

        vptr(i) = (one - volCoef)*Vtmp(i) + volCoef/nNeighbors*vSum

     end do

     call VecGhostUpdateBegin(volume, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecGhostUpdateEnd(volume, INSERT_VALUES,SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do

  deallocate(vtmp)
  call VecRestoreArrayF90(VolumeLocal, Vptr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGhostRestoreLocalForm(Volume, VolumeLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now we do a global volume smooth
  factor = half * (one + ( one - volBlend)**(marchIter-2))

  call VecGetArrayF90(Volume, Vptr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  Vptr = factor*VPtr + (one - factor)*Vbar

  call VecRestoreArrayF90(Volume, Vptr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine volumeSmooth

subroutine initialGuess(Xnew)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: initialGuess determines a good initial guess for the
  !     next grid level. In fact this function implements the a pseudo
  !     stepping version of the Chean and Steger algorithm. If the
  !     nonLinear option is False, only this routine is called. Note
  !     that sensitivty information cannot be performed with the
  !     linear method because the intermediate grid levels are not
  !     stored. 
  !
  !     Ney Secco (2015-2016): Added boundary conditions
  !
  !     Description of Arguments
  !
  !     Ouput:
  !     xNew petsc-vector : Result of linear step.
  !


  use communication
  use hypData, only : nx, deltas, X, XLm1, XL, marchIter, volume
  use hypData, only : cRatio, hypKSP, hypMat, hypRHS, hypDelta, scaleDist
  use hypData, only : kspIts, nsubiter, nsubiterprev, gridRatio, desiredS, radius0

  use hypInput
  implicit none
  real(kind=realType) :: sl
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Output Parameters
  Vec :: Xnew

  ! Working Variables
  integer(kind=intType) :: i, ierr, jp1, jm1, kp1, km1, m, mm, ii, idim
  real(kind=realType) :: fact, onemfact, timeA, timeB
  logical :: keepGoing

  ! ===============
  ! BEGIN EXECUTION
  ! ===============

  ! Desired s:
  desiredS = desiredS + s0*gridRatio**(marchIter-2)

  ! Get values for the very first step; otherwise we will always have
  ! old values
  if (marchIter == 2) then 
     call VecCopy(X(1), XL, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecCopy(X(1), XLm1, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGhostUpdateBegin(XL, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call VecGhostUpdateEnd(XL, INSERT_VALUES,SCATTER_FORWARD, ierr)
     call VecGhostUpdateBegin(XLm1, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call VecGhostUpdateEnd(XLm1, INSERT_VALUES,SCATTER_FORWARD, ierr)

  end if
  keepGoing = .True.
  nsubIter = 0
  do while (keepGoing) 
     nSubIter = nSubIter + 1

     ! Compute volumes based on deltaS
     call computeVolumes
     ! Adjust deltaS if necessary
     if (cratio > cmax) then 
        deltaS = deltaS*(cmax/cratio)
        call computeVolumes
     end if

     ! And then smooth as normal
     call volumeSmooth

     ! Increment our running counter on approx distance from wall
     scaleDist = scaleDist + deltaS
     
     ! Assemble the Jacobaian, second order is true, and also gets the RHS
     call calcResidual()

     call KSPSetOperators(hypKSP, hypMat, hypMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! Now solve the system
     timeA = mpi_wtime()

     call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     timeB = mpi_wtime()

     call KSPGetIterationNumber(hypKsp, kspIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  
     ! Copy to old value
     call VecCopy(XL, XLm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! Get new XX level
     call VecAXPY(XL, one, hypDelta, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now we need to check if we have gone far enough for the next
     ! desired grid level. If so, we compute the factor between xxm1
     ! and xx and set. We make a linear interpolation between pseudo-layers
     ! to find the actual layer that we want.
  
     if (scaleDist >= desiredS) then
     
        fact = (desiredS - (scaleDist - deltaS))/(deltaS)
        onemfact = one-fact

        call VecAXPY(Xnew, fact, XL, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call VecAXPY(Xnew, onemfact, XLm1, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Break for the while loop
        keepGoing = .False.
     end if

     ! Update the ghost values 
     call VecGhostUpdateBegin(XL, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call VecGhostUpdateEnd(XL, INSERT_VALUES,SCATTER_FORWARD, ierr)
     call VecGhostUpdateBegin(XLm1, INSERT_VALUES, SCATTER_FORWARD, ierr)
     call VecGhostUpdateEnd(XLm1, INSERT_VALUES,SCATTER_FORWARD, ierr)

     ! Change deltaS
     deltaS = deltaS * pGridRatio
  end do
end subroutine initialGuess

subroutine calcResidual
  use communication
  use hypData
  use hypInput
  implicit none

  ! Working variables
  integer(kind=intType) :: i, idim, j, jp1, jp2, jm1, kp1, km1, ii, jj, m, MM, iter, ierr
  integer(kind=intType) :: iLow, iHigh, bcType
  !real(kind=realType) :: XwallBC, YwallBC, ZwallBC
  real(kind=realType), dimension(3) :: r_zeta, r_ksi, r_eta, sigma, t, tau
  real(kind=realType), dimension(3) :: r_ksi_ksi, r_eta_eta, rA, rB, rC, f_rhs
  real(kind=realType), dimension(3) :: rMinusjHat, rMinuskHat, r0, rmmr0
  real(kind=realType), dimension(3) :: rPlusjHat, rPluskHat, nhat, G, PinvG
  real(kind=realType), dimension(3,3) :: Q1, Q2, P, eye, PinvQ1, PinvQ2, Pinv
  real(kind=realType) :: a_ksi, a_eta, De(3), factor, alpha, beta, ovrSum1, ovrSum2
  real(kind=realType) :: d_ksi, d_eta, tmp1, tmp2, tmp, thetam, Rksi, Reta, N_ksi, N_eta
  real(kind=realType) :: dist, coef,  deltaVovrDetC, detC, largeWeight, smallWeight, eps_x
  ! This routine must do everything required to evaluate the full
  ! non-linear residual

  real(kind=realType), dimension(3) :: delta_x
  real(kind=realType), dimension(3,3) :: Abc

  ! ===============
  ! BEGIN EXECUTION
  ! ===============

  sl = (scaleDist/(radius0*rMin))**slExp

  ! Record the maximum and minimum grid sensors
  gridSensorMax = zero 
  gridSensorMin = huge(gridSensorMin)

  call VecGetOwnershipRange(XL, iLow, iHigh, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  eye = zero
  eye(1,1) = one
  eye(2,2) = one
  eye(3,3) = one
  call MatZeroEntries(hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecZeroEntries(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Get the local forms of the XL and XLm1
  call VecGhostGetLocalForm(XL, XL_local, ierr)
  call VecGhostGetLocalForm(XLm1, XLm1_local, ierr)

  call VecGetArrayF90(XL_local, xx, ierr)
  call VecGetArrayF90(XLm1_local, xxm1, ierr)
  call VecGetArrayF90(Volume, Vptr, ierr)

  ! Nodal loop
  do i=1, nx
     if (BCneighborsLocal(i,1) .gt. 0) then ! We have boundary conditions on this node
        ! Identify index of the current node on the global block matrix
        ii = i + iLow/3

        ! Read BC type
        bcType = BCneighborsLocal(i,1)

        ! -----------------
        ! Left Hand Side:
        ! -----------------

        ! Remember that the BC indices are defined in hypData.F90

        findBCtype: if (bcType .eq. SplayBCindex) then

           ! Now we will read the neighbors according to the orientation
           ! used to define patch%BCnodes (check readCGNS.F90). This
           ! orientation will depend on the BC type

           ! For the Splay BC:
           ! The i and j indices used in this part of the code are taken
           ! with respect to the marching direction. Therefore, the jp1
           ! node will be towards the interior of the mesh always
           !
           !           jp1
           !            |
           !            |
           !  km1----refNode----kp1
           !     ---------> marching direction
           !
           ! Things are slightly different for corners because we
           ! include the diagonal node:
           !
           !            X
           !            |
           !            |
           !           km1------jp1
           !            |        |
           !            |        |
           !         refNode----kp1----X
           !     ---------> marching direction

           ! Read neighbors in local form
           kp1 = BCneighborsLocal(i, 2)
           jp1 = BCneighborsLocal(i, 3)
           km1 = BCneighborsLocal(i, 4)

           ! Find reference vectors
           rA = xx(3*jp1-2:3*jp1) - xx(3*i-2:3*i)
           rB = (xx(3*kp1-2:3*kp1) - xx(3*km1-2:3*km1))/2
           call cross_prod(rB, rA, rC) !Subroutine defined in 3D_utilities

           ! Assemble matrix
           Abc(1,:) = rA
           Abc(2,:) = rB
           Abc(3,:) = rC

           ! The splay factors come from the user-supplied options

           ! Get global indices of neighbors to assemble the matrix
           kp1 = BCneighborsGlobal(i, 2)
           jp1 = BCneighborsGlobal(i, 3)
           km1 = BCneighborsGlobal(i, 4)

           ! Assemble matrices
           ! Center Point !
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           ! km1 Point !
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, km1-1, &
                -nuSplay*0.5*Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           ! kp1 Point !
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, kp1-1, &
                -nuSplay*0.5*Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
                (/zero, zero, (1.0-sigmaSplay)*(1.0-nuSplay)*deltaS*sqrt(rC(1)*rC(1) + &
                rC(2)*rC(2) + rC(3)*rC(3))/), INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

        else if ((bcType .eq. SymmetryXBCindex) .or.(bcType .eq. SymmetryYBCindex) .or. &
             (bcType .eq. SymmetryZBCindex)) then ! Symmetry BC

           ! Now we will read the neighbors according to the orientation
           ! used to define patch%BCnodes (check readCGNS.F90). This
           ! orientation will depend on the BC type

           ! For the Symmetry BC we need two neighbors:
           ! The i and j indices used in this part of the code are taken
           ! with respect to the marching direction. Therefore, the jp1
           ! node will be towards the interior of the mesh always.
           !
           !           jp2
           !            |
           !            |
           !           jp1
           !            |
           !            |
           !   X----refNode----X
           !     ---------> marching direction

           ! Read neighbors in local form
           jp1 = BCneighborsLocal(i, 2)
           jp2 = BCneighborsLocal(i, 3)

           ! Read coordinates
           rA = xx(3*ii-2:3*ii)
           rB = xx(3*jp1-2:3*jp1)
           rC = xx(3*jp2-2:3*jp2)

           ! Define auxiliary matrix that will blank the symmetric coordinate
           Abc = eye

           ! Compute the right-hand side expression
           f_rhs = 3.0*rA - 4.0*rB + rC

           ! Blank the elements that corresponds to the symmetry variables
           if (bcType .eq. SymmetryXBCindex) then !X-sym
              Abc(1,1) = 0
              f_rhs(1) = 0
           else if (bcType .eq. SymmetryYBCindex) then !Y-sym
              Abc(2,2) = 0
              f_rhs(2) = 0
           else if (bcType .eq. SymmetryZBCindex) then !Z-sym
              Abc(3,3) = 0
              f_rhs(3) = 0
           end if
           
            ! Read neighbors in global form to assemble matrices
           jp1 = BCneighborsGlobal(i, 2)
           jp2 = BCneighborsGlobal(i, 3)

           ! Set the matrix of the reference node ! 
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                -3.0*Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Assign matrix to jp1 node
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, jp1-1, &
                4.0*Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Assign matrix to jp2 node
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, jp2-1, &
                -1.0*Abc, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
                f_rhs, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

        else if ((bcType .eq. ConstXBCindex) .or.(bcType .eq. ConstYBCindex) .or. &
             (bcType .eq. ConstZBCindex)) then ! Wall BC

           ! Now we will read the neighbors according to the orientation
           ! used to define patch%BCnodes (check readCGNS.F90). This
           ! orientation will depend on the BC type

           ! For the Wall BC we need only one neighbor:
           ! The i and j indices used in this part of the code are taken
           ! with respect to the marching direction. Therefore, the jp1
           ! node will be towards the interior of the mesh always
           !
           !           jp1
           !            |
           !            |
           !   X----refNode----X
           !     ---------> marching direction

           ! Read neighbors in local form
           jp1 = BCneighborsLocal(i, 2)

           ! Set the matrix of the reference node ! 
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                eye, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! First neighbor
           Abc = -eye
           ! Set to zero the corresponding coordinate
           if (bcType .eq. ConstXBCindex) then !X-wall
              Abc(1,1) = 0
           else if (bcType .eq. ConstYBCindex) then !Y-Wall
              Abc(2,2) = 0
           else if (bcType .eq. ConstZBCindex) then !Z-Wall
              Abc(3,3) = 0
           end if

           ! Read neighbors in global form to assemble matrices
           jp1 = BCneighborsGlobal(i, 2)

           ! Assign matrix
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, jp1-1, &
                (Abc), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
                (/zero, zero, zero/), INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

        end if findBCtype

     else ! We have an interior point

        if (lnptr(1, i) == 4) then
           ! Extract pointers to make code easier to read
           jp1 = lnPtr(2, i)
           kp1 = lnPtr(3, i)
           jm1 = lnPtr(4, i)
           km1 = lnPtr(5, i)

           ! Compute centered difference with respect to ksi and eta. This is the
           ! nominal calculation that we will use in the interior of the domain
           r_ksi = half*(xx(3*jp1-2:3*jp1) - xx(3*jm1-2:3*jm1))
           r_eta = half*(xx(3*kp1-2:3*kp1) - xx(3*km1-2:3*km1))
           r_ksi_ksi = half*(xx(3*jp1-2:3*jp1) - two*xx(3*i-2:3*i) + xx(3*jm1-2:3*jm1))
           r_eta_eta = half*(xx(3*kp1-2:3*kp1) - two*xx(3*i-2:3*i) + xx(3*km1-2:3*km1))

           ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
           tmp = (dist(xxm1(3*jp1-2:3*jp1), xxm1(3*i-2:3*i)) + &
                dist(xxm1(3*jm1-2:3*jm1), xxm1(3*i-2:3*i))) / &
                (dist(xx(3*jp1-2:3*jp1), xx(3*i-2:3*i)) + &
                dist(xx(3*jm1-2:3*jm1), xx(3*i-2:3*i)))
           d_ksi = max(tmp**(two/sl), 0.1)


           tmp = (dist(xxm1(3*kp1-2:3*kp1), xxm1(3*i-2:3*i)) + &
                dist(xxm1(3*km1-2:3*km1), xxm1(3*i-2:3*i))) / &
                (dist(xx(3*kp1-2:3*kp1), xx(3*i-2:3*i)) + &
                dist(xx(3*km1-2:3*km1), xx(3*i-2:3*i)))

           d_eta = max(tmp**(two/sl), 0.1)
        else 
           ! Lets do this brute force:
           r_ksi = zero
           r_eta = zero
           r_ksi_ksi = zero
           r_eta_eta = zero

           r0 = xx(3*i-2:3*i)
           MM = lnPtr(1, i) 

           tmp1 = zero
           tmp2 = zero
           ovrSum1 = zero
           ovrSum2 = zero
           do m = 0, MM - 1
              thetam = 2*pi*m/MM
              jj = lnptr(2 + m, i)

              rmmr0 = xx(3*jj-2:3*jj)- r0
              r_ksi = r_ksi + (two/MM) * (rmmr0) * cos(thetam)
              r_eta = r_eta + (two/MM) * (rmmr0) * sin(thetam)
              r_ksi_ksi = r_ksi_ksi + (two/MM)*(rmmr0)*(four*cos(thetam)**2 - one)
              r_eta_eta = r_eta_eta + (two/MM)*(rmmr0)*(four*sin(thetam)**2 - one)

              tmp1 = tmp1 + &
                   dist(xxm1(3*jj-2:3*jj), xxm1(3*i-2:3*i)) / &
                   dist(  xx(3*jj-2:3*jj),   xx(3*i-2:3*i))*abs(cos(thetam))

              tmp2 = tmp2 + &
                   dist(xxm1(3*jj-2:3*jj), xxm1(3*i-2:3*i)) / &
                   dist(  xx(3*jj-2:3*jj),   xx(3*i-2:3*i))*abs(sin(thetam))

              ovrSum1 = ovrSum1 + abs(cos(thetam))
              ovrSum2 = ovrSum2 + abs(sin(thetam))
           end do
           d_ksi = max((tmp1/ovrSum1)**(two/sl), 0.1)
           d_eta = max((tmp2/ovrSum2)**(two/sl), 0.1)
        end if

        gridSensorMax = max(gridSensorMax, d_ksi, d_eta)
        gridSensorMin = min(gridSensorMin, d_ksi, d_eta)

        ! Compute the auxiluary variables sigma, t and tau which are
        ! the three cross products of the vectors
        sigma(1) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
        sigma(2) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
        sigma(3) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)

        if (.not. nonLinear) then
           detC = &
                (r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))**2 + &
                (r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))**2 + &
                (r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))**2

           deltaVovrDetC = vPtr(i)/detC

           r_zeta(1) = deltaVovrDetC*(r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))
           r_zeta(2) = deltaVovrDetC*(r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))
           r_zeta(3) = deltaVovrDetC*(r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))
        end if

        t(1) = r_eta(2)*r_zeta(3) - r_zeta(2)*r_eta(3)
        t(2) = r_zeta(1)*r_eta(3) - r_eta(1)*r_zeta(3)
        t(3) = r_eta(1)*r_zeta(2) - r_zeta(1)*r_eta(2)

        tau(1) = r_zeta(2)*r_ksi(3) - r_ksi(2)*r_zeta(3)
        tau(2) = r_ksi(1)*r_zeta(3) - r_zeta(1)*r_ksi(3)
        tau(3) = r_zeta(1)*r_ksi(2) - r_ksi(1)*r_zeta(2)

        Q1(1, :) = r_zeta
        Q1(2, :) = zero
        Q1(3, :) = t

        Q2(1, :) = zero
        Q2(2, :) = r_zeta
        Q2(3, :) = tau

        P(1, :) = r_ksi
        P(2, :) = r_eta
        P(3, :) = sigma

        ! Compute the inverse of P and its multiplcation with Q1 and Q2
        call three_by_three_inverse(P, Pinv)
        PinvQ1 = matmul(Pinv, Q1)
        PinvQ2 = matmul(Pinv, Q2)

        ! Compute the 'g' vector
        g(1) = r_ksi(1)*r_zeta(1) + r_ksi(2)*r_zeta(2) + r_ksi(3)*r_zeta(3)
        g(2) = r_eta(1)*r_zeta(1) + r_eta(2)*r_zeta(2) + r_eta(3)*r_zeta(3)
        g(3) = vPtr(i) + two*(&
             sigma(1)*r_zeta(1) + sigma(2)*r_zeta(2) + sigma(3)*r_zeta(3))

        ! And the residual contribution
        Pinvg = matmul(Pinv, g)

        ! ============================================
        !             Explicit Smoothing 
        ! ============================================

        if (abs(lnptr(1, i)) == 4) then !We have the usual 4 neighbors

           ! Equation 6.9
           rPlusjHat = (xx(3*jp1-2:3*jp1) - xx(3*i-2:3*i))/dist(xx(3*jp1-2:3*jp1), xx(3*i-2:3*i))
           rMinusjHat = (xx(3*jm1-2:3*jm1) - xx(3*i-2:3*i))/dist(xx(3*jm1-2:3*jm1), xx(3*i-2:3*i))

           rPluskHat = (xx(3*kp1-2:3*kp1) - xx(3*i-2:3*i))/dist(xx(3*kp1-2:3*kp1), xx(3*i-2:3*i))
           rMinuskHat = (xx(3*km1-2:3*km1) - xx(3*i-2:3*i))/dist(xx(3*km1-2:3*km1), xx(3*i-2:3*i))

           ! Equation 6.10 (cross_prod defined in 3D_utilities.F90)
           call cross_prod(rPlusjHat-rMinusjHat,rPluskHat-rMinuskHat, nhat)
           nhat = nhat/norm2(nhat)

           ! Equation 6.11
           alpha = acos(dot_product(nhat,rPlusjHat))
           beta = acos(dot_product(nhat,rPluskHat))

           ! Equation 6.12
           if (alpha < 3.1415/2.0) then
              a_ksi = one/(one - cos(alpha)*cos(alpha))
           else
              a_ksi = one
           end if

           if (beta < 3.1415/2.0) then
              a_eta = one/(one - cos(beta)*cos(beta))
           else
              a_eta = one
           end if

        else ! We don't have 4 neighbors
           a_ksi = one
           a_eta = one
        end if

        ! ------------------------------------
        ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
        ! ------------------------------------
        Rksi = Sl * d_ksi * a_ksi
        Reta = Sl * d_eta * a_eta

        ! ------------------------------------
        ! Matrix Norm approximation  (Equation 6.3)
        ! ------------------------------------

        N_ksi = sqrt((r_zeta(1)**2 + r_zeta(2)**2 + r_zeta(3)**2)/ &
             (r_ksi(1)**2 + r_ksi(2)**2 + r_ksi(3)**2))
        N_eta = sqrt((r_zeta(1)**2 + r_zeta(2)**2 + r_zeta(3)**2)/ &
             (r_eta(1)**2 + r_eta(2)**2 + r_eta(3)**2))

        ! ------------------------------------
        ! Final explict dissipation (Equation 6.1 and 6.2)
        ! ------------------------------------
        De = epsE*(Rksi*N_ksi*r_ksi_ksi + Reta*N_eta*r_eta_eta)

        if (marchIter == 2) then
           De = zero
        end if

        ! Assemble the LHS, once again we need to split between regular and extraordinary nodes.
        ii = i + iLow/3
        if (gnptr(1, i) == 4) then
           ! Extract pointers to make code easier to read

           jp1 = gnPtr(2, i)
           kp1 = gnPtr(3, i)
           jm1 = gnPtr(4, i)
           km1 = gnPtr(5, i)

           ! -----------------
           ! Left Hand Side:
           ! -----------------

           ! Point to the left
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, jm1-1, &
                (-(one+theta)*half*PinvQ1 - epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Center Point ! 
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                (eye + four*epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Point to the right
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, jp1-1, &
                ((one+theta)*half*PinvQ1 - epsI*eye),  ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! ----------------------------------------------------------

           ! Point to the bottom
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, km1-1, &
                (-(one+theta)*half*PinvQ2 - epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Point to the right
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, kp1-1, &
                ((one+theta)*half*PinvQ2 - epsI*eye),  ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
                matmul(Pinv,(/zero, zero, vPtr(i)/)) + De, INSERT_VALUES, ierr) !!!
           call EChk(ierr, __FILE__, __LINE__)

        else 
           MM = gnPtr(1, i) 
           do m = 0, MM - 1
              thetam = 2*pi*m/MM
              ! For the 'fm' point in ksi
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, gnPtr(2+m, i)-1, &
                   (one + theta)*PinvQ1*(two/MM)*cos(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in ksi
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                   -(one + theta)*PInvQ1*(two/MM)*cos(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'fm' point in eta
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, gnPtr(2+m, i)-1, &
                   (one + theta)*PinvQ2*(two/MM)*sin(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in ksi
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                   -(one + theta)*PInvQ2*(two/MM)*sin(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! Now we have the second derivative values to do in ksi-ksi

              ! For the 'fm' point in ksi-ksi 
              coef = (two/MM)*(four*cos(thetam)**2 - one)

              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, gnPtr(2+m, i)-1, &
                   -eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in ksi-ksi
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                   eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'fm' point in eta-eta
              coef = (two/MM)*(four*sin(thetam)**2 - one)

              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, gnPtr(2+m, i)-1, &
                   -eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in eta-eta
              call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                   eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do

           ! Finally we need an eye on the diagonal
           call MatSetValuesBlocked(hypMat, 1, ii-1, 1, ii-1, &
                eye, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__) 

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
                matmul(Pinv,(/zero, zero, vPtr(i)/)) + De, INSERT_VALUES, ierr) !!!
           call EChk(ierr, __FILE__, __LINE__)
        end if

        if (writeMetrics .and. metricsAllocated) then
           call VecSetValuesBlocked(metrics(marchIter-1, iX_ksi), 1, (/ii-1/), r_ksi, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iX_eta), 1, (/ii-1/), r_eta, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iX_zeta), 1, (/ii-1/), r_zeta, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iX_ksi_ksi), 1, (/ii-1/), r_ksi_ksi, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iX_eta_eta), 1, (/ii-1/), r_eta_eta, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iX_diss), 1, (/ii-1/), De, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetValuesBlocked(metrics(marchIter-1, iVHist), 1, (/ii-1/), (/zero,zero, vPtr(i)/), INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if
     end if ! BC check
  end do ! NODE LOOP

  call MatAssemblyBegin(hypMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(hypMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyBegin(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (writeMetrics .and. metricsAllocated) then
     do j=1,nMetric
        call vecAssemblyBegin(metrics(marchIter-1, j), ierr)
        call vecAssemblyEnd(metrics(marchIter-1, j), ierr)
     end do
  end if

  ! Finally restore everything
  call VecRestoreArrayF90(XL_local, xx, ierr)
  call VecRestoreArrayF90(XLm1_local, xxm1, ierr)
  call VecRestoreArrayF90(Volume, vPtr, ierr)

  call VecGhostRestoreLocalForm(XL, XL_local, ierr)
  call VecGhostRestoreLocalForm(XLm1, XLm1_local, ierr)

end subroutine calcResidual

subroutine create3DPetscVars
  use communication
  use hypInput
  use hypData

  implicit none
#include "include/petscversion.h"  
  ! Working Variables
  integer(kind=intType) :: ierr, i, j, idim, bs, dummy(1)
  integer(kind=intType), dimension(:), allocatable :: onProc, offProc

  ! ----------------------------------------------------------
  !          Linearized Hyperbolic System Variables
  ! ----------------------------------------------------------
  if (.not. three_d_vars_allocated) then
     
     ! Lets to things in the proper way, create the Mat first
     allocate(onProc(nx), offProc(nx), stat=ierr)
     call EChk(ierr, __FILE__, __LINE__)
     do i=1,nx
        onProc(i) = min(50, nx)
        offProc(i) = min(nxglobal-nx, 50)
     end do

     ! Create a blocked matrix
     bs = 3
     call MatCreateBAIJ(hyp_comm_world, bs, &
          nx*bs, nx*bs, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, onProc, 0, offProc, hypMat, ierr)

     call EChk(ierr, __FILE__, __LINE__)
     deallocate(onProc, offProc, stat=ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! This must be set to allow passing in blocks in native fortran
     ! ordering
     call MatSetOption(hypMat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Then use getVecs to get the vectors we want
#if PETSC_VERSION_MINOR > 5
     call MatCreateVecs(hypMat, hypDelta, hypRHS, ierr)
#else
   call MatGetVecs(hypMat, hypDelta, hypRHS, ierr)
#endif
     call EChk(ierr, __FILE__, __LINE__)

     ! Create the extra state-sized vectors
     call VecDuplicate(hypDelta, hypRes, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  
     ! Create the full list of grid vectors:
     allocate(X(N), stat=ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Lots of additional vectors for metrics if necessary
     if (writeMetrics) then
        allocate(metrics(N, nMetric), stat=ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if

     ! Chreate ghosted vectors for the working pseudo grid levels
     if (nGhost == 0) then
        call VecCreateGhostBlock(hyp_comm_world, 3, 3*nx, PETSC_DETERMINE, &
             nGhost, dummy, XL, ierr)
     else
        call VecCreateGhostBlock(hyp_comm_world, 3, 3*nx, PETSC_DETERMINE, &
             nGhost, ghost, XL, ierr)
     end if
     call EChk(ierr, __FILE__, __LINE__)
     call VecDuplicate(XL, XLm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! The volume array also needs to be ghosted.
     if (nGhost == 0) then
        call VecCreateGhost(hyp_comm_world, nx, PETSC_DETERMINE, &
             nGhost, dummy, Volume, ierr)
     else
        call VecCreateGhost(hyp_comm_world, nx, PETSC_DETERMINE, &
             nGhost, ghost, Volume, ierr)
     end if
     call EChk(ierr, __FILE__, __LINE__)
   
     ! Loop over each "vec" in the X/metric arrays:
     do i=1,N
        call VecDuplicate(XL, X(i), ierr)

        if (writeMetrics) then
           do j=1,nMetric
              call VecDuplicate(hypDelta, metrics(i, j), ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           metricsAllocated = .True.
        end if
     end do
     
     ! Create a scatter context to root:
     call VecScatterCreateToZero(X(1), rootScatter, XLocal, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call vecDuplicate(xLocal, Xlocalm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! This should come as a option, but we'll put it here for now.
     call PetscOptionsSetValue("-sub_pc_factor_levels", "1", ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call PetscOptionsSetValue("-pc_type", "asm", ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! call PetscOptionsSetValue("-ksp_monitor", "", ierr)
     ! call EChk(ierr, __FILE__, __LINE__)

     !call PetscOptionsSetValue("-ksp_monitor_true_residual", "", ierr)
     !call EChk(ierr, __FILE__, __LINE__)
     
     call PetscOptionsSetValue("-sub_pc_asm_overlap", "1", ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Create the KSP object
     call KSPCreate(petsc_comm_world, hypKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetFromOptions(hypKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call KSPGMRESSetRestart(hypKSP, kspSubspacesize, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetOperators(hypKSP, hypMat, hypMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetTolerances(hypKSP, kspRelTol, 1e-30, 1e5, kspMaxIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     three_d_vars_allocated = .True.
  end if
end subroutine create3DPetscVars


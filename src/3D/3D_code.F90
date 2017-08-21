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

  call calcGridRatio(N, nConstantStart, nConstantEnd, s0, marchDist, gridRatio)

  ! Defaults/Error checking
  if (pGridRatio <=0) then 
     pGridRatio = gridRatio
  end if
  
  if (pGridRatio > gridRatio) then 
     if (myid == 0) then 
        print *,'Erorr: The supplied pseudo grid ratio is too large. It must be less than', gridRatio
     end if
     call mpi_barrier(hyp_comm_world, ierr)
     stop
  end if

  ! Set initial pseudo and real spacings
  if (ps0 <= zero) then 
     deltaS = s0/two
  else
     deltaS = ps0
  end if

  if (deltaS > s0) then
     if (myid == 0) then 
        print *,'Erorr: The supplied pseudo grid s0 is too large. It must be less than', s0
     end if
     call mpi_barrier(hyp_comm_world, ierr)
     stop
  end if

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

     ! Compute the how far this layer should go:
     call computeStretch(L)

     ! Run the "initial guess" function. If the user is running in
     ! 'linear' mode this is all that is done. This function computes
     ! the next step, X(L). It may take multiple sub-steps to get there
     ! which is fine. 
     call initialGuess(X(L))
     
     ! Check the quality of this layer
     call computeQualityLayer

     ! Compute the stretching
     if (L .gt. 2) then
        call findKStretch(X(L), X(L-1), X(L-2))
     else
        maxKStretch = zero
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
  real(kind=realType) :: ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS
  real(kind=realType) :: c1,c2,c3,c4, dist
  real(kind=realType) :: Vbar_local, cratio_local
  integer(kind=intType), dimension(:), allocatable :: nodeCount
  
  allocate(nodeCount(nx+nGhost))
  nodeCount = 0

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
     
     areadS = deltaS*half*sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))

     ! Scatter full area back to the nodes. 
     vPtr(n1) = vPtr(n1) + areadS
     vPtr(n2) = vPtr(n2) + areadS
     vPtr(n3) = vPtr(n3) + areadS
     vPtr(n4) = vPtr(n4) + areadS
     nodeCount(n1) = nodeCount(n1) + 1
     nodeCount(n2) = nodeCount(n2) + 1
     nodeCount(n3) = nodeCount(n3) + 1
     nodeCount(n4) = nodeCount(n4) + 1

     vBar_local = vBar_local + areadS

     ! Also compute edge lengths.
     c1 = deltaS/dist(ll,lr)
     c2 = deltaS/dist(ul,ur)
     c3 = deltaS/dist(ll,ul)
     c4 = deltaS/dist(lr,ur)

     cRatio_local = max(cRatio_local, c1, c2, c2, c4)
  end do
  
  ! Finally we have to divide by the count
  do i=1, nx
     Vptr(i) = VPtr(i) / nodeCount(i)
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
  deallocate(nodeCount)
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
  real(kind=realType) :: factor, nNeighbors, vSum, frac, low, high, frac2
  integer(kind=intType) :: i, iSize, iter, ierr, ii,jj, nIter
  real(kind=realType), allocatable, dimension(:) :: Vtmp
  logical :: search

  if (allocated(volSmoothSchedule)) then 
     ! We need to determine how many smoothing iterations to do by
     ! interpolating the volumeSmoothing Schedule
     frac = (marchIter-one)/(N-1)
     
     ! Just do a linear search for the bin:
     do i=1,size(volSmoothSchedule, 1)-1
        if (frac >= volSmoothSchedule(i, 1) .and. frac <= volSmoothSchedule(i+1, 1)) then 
           frac2 = (frac - volSmoothSchedule(i, 1))/ &
                (volSmoothSchedule(i+1, 1) - volSmoothSchedule(i, 1))
           low = volSmoothSchedule(i, 2)
           high = volSmoothSchedule(i+1, 2)
           nIter = int(low + frac2*(high-low))
        end if
     end do
  else
     nIter = volSmoothIter
  end if
  
  ! Do a Jacobi volume smooth
  call VecGhostGetLocalForm(Volume, VolumeLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetSize(VolumeLocal, iSize, ierr)
  allocate(vTmp(iSize))

  call VecGetArrayF90(VolumeLocal, Vptr, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  do iter=1, nIter

     ! Copy vptr to vtmp
     do i=1, isize
        vtmp(i) = vptr(i)
     end do

     do i=1,nx

        vSum = zero
        nNeighbors = lnPtr(1,i)

        do ii=1, lnPtr(1,i)
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
  use hypData, only : kspIts, nsubiter, nsubiterprev, gridRatio, desiredS
  use hypData, only : desiredDeltaS
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
  real(kind=realType) :: fact, onemfact, timeA, timeB, vals(1)
  logical :: keepGoing

  ! Desired s:
  desiredS = desiredS + desiredDeltaS

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

     ! Increment our running counter on the acutal approx distance
     ! from wall
     scaleDist = scaleDist + deltaS
     
     ! Assemble the Jacobaian, second order is true, and also gets the RHS
     call calcResidual()
     call setupPETScKsp()
     
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

     ! Update BCs for this this mesh if necessary. 
     call updateBCs

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

     ! Update the deltaS. We always multiply by but never let it go
     ! higher than desired deltas
     deltaS = deltaS * pGridRatio
     deltaS = min(deltaS, desiredDeltaS)

  end do
end subroutine initialGuess

subroutine calcResidual
  use communication
  use hypData
  use hypInput
  implicit none

  ! Working variables
  integer(kind=intType) :: i, idim, j, jp1, jp2, jm1, kp1, km1, ii, jj, m, MM, iter, ierr
  integer(kind=intType) :: iStart, iEnd
  real(kind=realType), dimension(3) :: xx0, xjp1, xjm1, xkp1, xkm1
  real(kind=realType), dimension(3) :: xx0_m1, xjp1_m1, xjm1_m1, xkp1_m1, xkm1_m1
  real(kind=realType), dimension(3) :: r_zeta, r_ksi, r_eta, sigma, t, tau
  real(kind=realType), dimension(3) :: r_ksi_ksi, r_eta_eta
  real(kind=realType), dimension(3) :: rMinusjHat, rMinuskHat, r0, rmmr0
  real(kind=realType), dimension(3) :: rPlusjHat, rPluskHat, nhat, G, PinvG
  real(kind=realType), dimension(3,3) :: Q1, Q2, P, eye, PinvQ1, PinvQ2, Pinv
  real(kind=realType), dimension(3,3) :: blk0, blk1, blk2, blk3, blk4
  real(kind=realType) :: a_ksi, a_eta, De(3), factor, alpha, beta, ovrSum1, ovrSum2, maxAlpha
  real(kind=realType) :: d_ksi, d_eta, tmp1, tmp2, tmp, thetam, Rksi, Reta, N_ksi, N_eta
  real(kind=realType) :: dist, coef,  deltaVovrDetC, detC, largeWeight, smallWeight, eps_x
  logical :: averageNode

  ! METRIC CORRECTION
  real(kind=realType) :: r_zeta_p(3), r_eta_p(3), r_ksi_p(3), nu_mc
  real(kind=realType) :: rPlusj(3), rMinusj(3), rPlusk(3), rMinusk(3), r_cross(3)

  ! CORNER CHECK
  real(kind=realType) :: normFace1(3), normFace2(3), normFace3(3), normFace4(3)
  real(kind=realType) :: normFace12(3), normFace13(3), normFace24(3), normFace34(3)
  real(kind=realType) :: nAux(3)

  ! This routine must do everything required to evaluate the full
  ! non-linear residual
  sl = (scaleDist/(marchDist))**slExp

  ! Record the maximum and minimum grid sensors
  gridSensorMax = zero 
  gridSensorMin = huge(gridSensorMin)
  nAverage = 0

  call VecGetOwnershipRange(XL, iStart, iEnd, ierr)
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
  masterNodeLoop: do i=1, nx
     averageNode = .False.

     ! Identify index of the current node on the global block matrix
     ii = i + iStart/3

     extraOrdinary: if (topoType(i) == topoInternal .and. lnptr(1, i) /= 4) then 
        ! This is an internal, extraonrdinary node. We do not have to
        ! worry about BCs for thsi one, so use the circular code
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
     else
        ! These nodes have 4 "logical" neighbours, even if some of
        ! them are actually computed halos. 

        ! Extract pointers for each of our neighbours. Some of these
        ! may not be meaningful. 
        jp1 = lnPtr(2, i)
        kp1 = lnPtr(3, i)
        jm1 = lnPtr(4, i)
        km1 = lnPtr(5, i)

        ! Extract our own node, the first two neighbours and the
        ! previous nodes
        xx0    =   xx(3*i-2:3*i)
        xx0_m1 = xxm1(3*i-2:3*i)

        xjp1    =   xx(3*jp1-2:3*jp1)
        xjp1_m1 = xxm1(3*jp1-2:3*jp1)

        xkp1    =   xx(3*kp1-2:3*kp1)
        xkp1_m1 = xxm1(3*kp1-2:3*kp1)

        if (topoType(i) == topoInternal .or. topoType(i) == topoLCorner) then 
           ! Since the extra case was done above, this must be a
           ! regular 4-neighbour node. 

           ! Two more neighbours are significant. 
           xjm1    = xx(3*jm1-2:3*jm1)
           xjm1_m1 = xxm1(3*jm1-2:3*jm1)

           xkm1    = xx(3*km1-2:3*km1)
           xkm1_m1 = xxm1(3*km1-2:3*km1)

        else if (topoType(i) == topoEdge) then

           ! Third neighbours are significant. 
           xjm1 = xx(3*jm1-2:3*jm1)
           xjm1_m1 = xxm1(3*jm1-2:3*jm1)

           ! Get the 4th point based on the BC Type:
           call getBC(BCType(1, i), bcVal(1, :, i), .True., splayEdgeOrthogonality,&
                xx0   , xjp1   , xkp1   , xjm1   , xkm1  )
           call getBC(BCType(1, i), bcVal(1, :, i), .True.,  splayEdgeOrthogonality,&
                xx0_m1, xjp1_m1, xkp1_m1, xjm1_m1, xkm1_m1)

        else if (topoType(i) == topoCorner) then

           ! Get the 3rd and 4th points based on the BCTypes. 

           call getBC(BCType(1, i), bcVal(1, :, i), .False., &
                splayCornerOrthogonality, xx0, xkp1, xjp1, xkm1, xjm1)

           call getBC(BCType(2, i), bcVal(2, :, i), .False., &
                splayCornerOrthogonality, xx0, xjp1, xkp1, xjm1, xkm1)

           call getBC(BCType(1, i), bcVal(1, :, i), .False., &
                splayCornerOrthogonality, xx0_m1, xkp1_m1, xjp1_m1, xkm1_m1, xjm1_m1)

           call getBC(BCType(2, i), bcVal(2, :, i), .False., &
                splayCornerOrthogonality, xx0_m1, xjp1_m1, xkp1_m1, xjm1_m1, xkm1_m1)
        end if
        
        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        r_ksi = half*(xjp1 - xjm1)
        r_eta = half*(xkp1 - xkm1)
        r_ksi_ksi = half*(xjp1 - two*xx0 + xjm1)
        r_eta_eta = half*(xkp1 - two*xx0 + xkm1)
        
        ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
        tmp = (dist(xjp1_m1, xx0_m1) + dist(xjm1_m1, xx0_m1)) /&
             (dist(xjp1, xx0) + dist(xjm1, xx0))
        d_ksi = max(tmp**(two/sl), 0.1)

        tmp = (dist(xkp1_m1, xx0_m1) + dist(xkm1_m1, xx0_m1)) /&
             (dist(xkp1, xx0) + dist(xkm1, xx0))
        d_eta = max(tmp**(two/sl), 0.1)

     end if extraOrdinary
     
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


     if (topoType(i) == topoInternal .and. lnptr(1, i) /= 4) then 
        a_ksi = one
        a_eta = one
     else
        ! Equation 6.9
        rPlusj = xjp1 - xx0
        rMinusj = xjm1 - xx0
        rPlusjHat = rPlusj/norm2(rPlusj)
        rMinusjHat = rMinusj/norm2(rMinusj)
        
        rPlusk = xkp1 - xx0
        rMinusk = xkm1 - xx0
        rPluskHat = rPlusk/norm2(rPlusk)
        rMinuskHat = rMinusk/norm2(rMinusk)

        ! Equation 6.10 (cross_prod defined in 3D_utilities.F90)
        call cross_prod(rPlusjHat-rMinusjHat, rPluskHat-rMinuskHat, nhat)
        nhat = nhat/norm2(nhat)
        
        ! Equation 6.11
        alpha = acos(dot_product(nhat, rPlusjHat))
        beta = acos(dot_product(nhat, rPluskHat))

        ! Equation 6.12
        if (alpha < pi/two) then
           a_ksi = one/(one - cos(alpha)**2)
        else
           a_ksi = one
        end if
        
        if (beta < pi/two) then
           a_eta = one/(one - cos(beta)**2)
        else
           a_eta = one
        end if

        ! Here we detect if we have a *very* sharp corner and modify
        ! what we do at this point. 

        ! Each node has four neighbors, As shown below:
        !
        !        +------(j,k+1)------+
        !        |  face 1 | face 2  |
        !        |         |         |
        !     (j-1,k)----(j,k)----(j+1,k)
        !        |         |         |
        !        |  face 3 | face 4  |
        !        +------(j,k-1)------+
        !
        !

        ! First we will compute approximate face normals for each face.
        ! For instance, for face 1 we compute the cross product of the
        ! vectors connecting (i,j) to (i,j+1) and (i-1,j), and so on

        call cross_prod(rPluskHat, rMinusjHat, normFace1)
        call cross_prod(rPlusjHat, rPluskHat, normFace2)
        call cross_prod(rMinuskHat, rPlusjHat, normFace4)
        call cross_prod(rMinusjHat, rMinuskHat, normFace3)

        ! Initialize maximum corner angle found so far
        maxAlpha = -1.57

        ! Check the corner angles between each pair of faces
        ! Faces 1 and 2
        call computeCornerAngle(normFace1, normFace2, rPluskHat, alpha)
        maxAlpha = max(maxAlpha, alpha)
        ! Faces 2 and 4
        call computeCornerAngle(normFace2, normFace4, rPlusjHat, alpha)
        maxAlpha = max(maxAlpha, alpha)
        ! Faces 4 and 3
        call computeCornerAngle(normFace4, normFace3, rMinuskHat, alpha)
        maxAlpha = max(maxAlpha, alpha)
        ! Faces 3 and 1
        call computeCornerAngle(normFace3, normFace1, rMinusjHat, alpha)
        maxAlpha = max(maxAlpha, alpha)

        ! Check maximum corner angle
        if (maxAlpha > pi - cornerAngle) then 
           ! Flag this node as having to be averaged:
           averageNode = .True. 
           !print *,xx0
        end if

        ! Metric correction

        ! Compute adjusted derivatives (Eq. 7.2)
        r_ksi_p = fourth*(norm2(rPlusj) + norm2(rMinusj))*(rPlusjHat - rMinusjHat)
        r_eta_p = fourth*(norm2(rPlusk) + norm2(rMinusk))*(rPluskHat - rMinuskHat)
           
        ! Compute the vector shown in RHS of Eq. 7.1
        call cross_prod(r_ksi_p, r_eta_p, r_cross)
        
        ! Now use Eq. 7.1
        r_zeta_p = vPtr(i)*r_cross/(r_cross(1)**2 + r_cross(2)**2 + r_cross(3)**2)

        ! Compute blending factor
        nu_mc = one-sl

        ! Now smooth the ksi derivative
        r_zeta = (1-nu_mc)*r_zeta + nu_mc*r_zeta_p
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

     ! Row index is ii
     ii = i + iStart/3

     extraOrdinary2: if (topoType(i) == topoInternal .and. gnptr(1, i) /= 4) then 
        MM = gnPtr(1, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           ! For the 'fm' point in ksi
           call addBlock(ii, gnPtr(2+m, i), (one + theta)*PinvQ1*(two/MM)*cos(thetam))

           ! For the 'f0' point in ksi
           call addBlock(ii, ii, -(one + theta)*PInvQ1*(two/MM)*cos(thetam))

           ! For the 'fm' point in eta
           call addBlock(ii, gnPtr(2+m, i), (one + theta)*PinvQ2*(two/MM)*sin(thetam))
           
           ! For the 'f0' point in ksi
           call addBlock(ii, ii, -(one + theta)*PInvQ2*(two/MM)*sin(thetam))

           ! Now we have the second derivative values to do in ksi-ksi
           
           ! For the 'fm' point in ksi-ksi 
           coef = (two/MM)*(four*cos(thetam)**2 - one)
           call addBlock(ii, gnPtr(2+m, i), -eye*epsI*coef)
           
           ! For the 'f0' point in ksi-ksi
           call addBlock(ii, ii, eye*epsI*coef)
           
           ! For the 'fm' point in eta-eta
           coef = (two/MM)*(four*sin(thetam)**2 - one)
           call addBlock(ii, gnPtr(2+m, i), -eye*epsI*coef)
           
           ! For the 'f0' point in eta-eta
           call addBlock(ii, ii, eye*epsI*coef)
        end do
        
        ! Finally we need an eye on the diagonal
        call addBlock(ii, ii, eye)

     else ! Nominally 'regular' nodes

        if (bcType(1, i) == BCAverage) then 
           averageNode = .True. 
        end if
        
        if (.not. averageNode) then 
           ! First generate the 5 blocks that we will eventually have to set:
           blk0 = (eye + four*epsI*eye) ! Center
           blk1 = ((one+theta)*half*PinvQ1 - epsI*eye)  ! Right (jp1)
           blk2 = ((one+theta)*half*PinvQ2 - epsI*eye)  ! Top (kp1)
           blk3 = (-(one+theta)*half*PinvQ1 - epsI*eye) ! Left (jm1)
           blk4 = (-(one+theta)*half*PinvQ2 - epsI*eye) ! Bottom (km1)
        else
           blk0 = eye
           blk1 = -fourth*eye
           blk2 = -fourth*eye
           blk3 = -fourth*eye
           blk4 = -fourth*eye
           nAverage = nAverage + 1
        end if
           
        ! Extract indices. Not all may be significant.
        jp1 = gnPtr(2, i)
        kp1 = gnPtr(3, i)
        jm1 = gnPtr(4, i)
        km1 = gnPtr(5, i)

        if (topoType(i) == topoInternal .or. topoType(i) == topoLCorner) then 
           ! This must be an internal 4 neighbour node since the
           ! extraordinary nodes with bcInternal are already taken
           ! care of

           call addBlock(ii, ii, blk0)
           call addBlock(ii, jp1, blk1)
           call addBlock(ii, kp1, blk2)
           call addBlock(ii, jm1, blk3)
           call addBlock(ii, km1, blk4)

        else if (topoType(i) == topoEdge) then

           ! Accumulate blk4 into blk0 and blk2
           call getBCBlocks(BCType(1, i), blk4, blk0, blk2)
           
           ! Center and first THREE blocks are real nodes
           call addBlock(ii, ii, blk0)
           call addBlock(ii, jp1, blk1)
           call addBlock(ii, kp1, blk2)
           call addBlock(ii, jm1, blk3)

        else if (topoType(i) == topoCorner) then
           
           ! Accumulate blk3 into blk0 and blk1
           call getBCBlocks(BCType(1, i), blk3, blk0, blk1)

           ! Accumulate blk4 into blk0 and blk2
           call getBCBlocks(BCType(2, i), blk4, blk0, blk2)

           ! Center and first TWO blocks are real nodes
           call addBlock(ii, ii, blk0)
           call addBlock(ii, jp1, blk1)
           call addBlock(ii, kp1, blk2)
        end if
     end if extraOrdinary2

     ! The RHS must tbe set independent of the type of node
     if (.not. averageNode) then 
        call VecSetValuesBlocked(hypRHS, 1, (/ii-1/), &
             matmul(Pinv,(/zero, zero, vPtr(i)/)) + De, INSERT_VALUES, ierr)
     end if
     call EChk(ierr, __FILE__, __LINE__)
     
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
  end do masterNodeLoop

  ! Assemble the matrix and vectors
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

contains
  subroutine addBlock(i, j, blk)
    implicit none
    integer(kind=intType), intent(in) :: i, j
    real(kind=realType), dimension(3,3), intent(in) :: blk
    ! Helper routine for adding values to hypMat. Given indices are
    ! assumed to be 1-based and the zero-based conversion is done
    ! here.
    call MatSetValuesBlocked(hypMat, 1, i-1, 1, j-1, blk, ADD_VALUES, ierr)
    call EChk(ierr, __FILE__, __LINE__)
  end subroutine addBlock
end subroutine calcResidual

subroutine updateBCs

  ! Update XL with correct boundary conditions

  use communication
  use hypData
  use hypInput
  implicit none

  ! Working variables
  integer(kind=intType) :: i, idim, j, jp1, jp2, jm1, kp1, km1, ierr
  integer(kind=intType) :: iStart, iEnd
  real(kind=realType), dimension(3) :: xx0, xjp1, xjm1, xkp1, xkm1
  
  call VecGetOwnershipRange(XL, iStart, iEnd, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGhostGetLocalForm(XL, XL_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(XL_local, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  masterNodeLoop: do i=1, nx

     jp1 = lnPtr(2, i)
     kp1 = lnPtr(3, i)
     jm1 = lnPtr(4, i)
     km1 = lnPtr(5, i)
     
     if (topoType(i) == topoEdge) then 

        ! Extract our own node
        xx0 = xx(3*i-2:3*i)
        
        ! First THREE neighbours are significant. 
        xjp1 = xx(3*jp1-2:3*jp1)
        xkp1 = xx(3*kp1-2:3*kp1)
        xjm1 = xx(3*jm1-2:3*jm1)
        
        ! Get the 4th point based on the BC Type:
        call getBC(BCType(1, i), bcVal(1, :, i), .True., splayEdgeOrthogonality,&
             xx0   , xjp1   , xkp1   , xjm1   , xkm1  )
        xx(3*i-2:3*i) = xx0
        
     else if (topoType(i) == topoCorner) then

        ! Extract our own node
        xx0 = xx(3*i-2:3*i)

        ! First TWO neighbours are significant. 
        xjp1 = xx(3*jp1-2:3*jp1)
        xkp1 = xx(3*kp1-2:3*kp1)
        
        ! Get the 3rd and 4th points based on the BCTypes
        call getBC(BCType(1, i), bcVal(1, :, i), .False., &
             splayCornerOrthogonality, xx0, xkp1, xjp1, xkm1, xjm1)
        
        call getBC(BCType(2, i), bcVal(2, :, i), .False., &
             splayCornerOrthogonality, xx0, xjp1, xkp1, xjm1, xkm1)
        
        ! BC may update our node
        xx(3*i-2:3*i) = xx0

     else if (topoType(i) == topoLCorner) then

        ! Extract our own node
        xx0 = xx(3*i-2:3*i)

        ! Other nodes are not significant. 
        call getBC(BCType(1, i), bcVal(1, :, i), .False., &
             splayCornerOrthogonality, xx0, xkp1, xjp1, xkm1, xjm1)
        
        ! BC may update our node
        xx(3*i-2:3*i) = xx0

     end if
  end do masterNodeLoop

  ! Finally restore everything
  call VecRestoreArrayF90(XL_local, xx, ierr)
  call VecGhostRestoreLocalForm(XL, XL_local, ierr)
end subroutine updateBCs

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
  if (.not. varsAllocated) then
     
     ! Lets to things in the proper way, create the Mat first
     allocate(onProc(nx), offProc(nx), stat=ierr)
     call EChk(ierr, __FILE__, __LINE__)
     do i=1, nx
        onProc(i) = min(10, nx)
        offProc(i) = min(nxglobal-nx, 10)
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
     do i=1, N
        call VecDuplicate(XL, X(i), ierr)

        if (writeMetrics) then
           do j=1, nMetric
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


     ! Create the KSP object
     call KSPCreate(hyp_comm_world, hypKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetFromOptions(hypKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetType(hypKSP, 'gmres', ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPGMRESSetRestart(hypKSP, kspSubspacesize, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetOperators(hypKSP, hypMat, hypMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPSetTolerances(hypKSP, kspRelTol, 1e-30, 1e5, kspMaxIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)


     
     varsAllocated = .True.
     !call MatSetOption(hypMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  end if
end subroutine create3DPetscVars


subroutine setupPETScKSP


  use hypData
  implicit none

  ! Special routine for setting up PETSc hypKSP object. It is getting
  ! progressively more difficult to manually setup KSP objects with
  ! PETSc, hence this routine is now needed for PETSc 3.7.
  integer(kind=intType) :: nlocal, first, ierr


  ! Setup the KSP object. 
  call KSPGetPC(hypKSP, globalPC, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call PCSetType(globalPC, "asm", ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call KSPSetup(hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Extract the ksp objects for each subdomain
  call PCASMGetSubKSP(globalPC, nlocal, first, subksp, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call KSPSetType(subksp, 'preonly', ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Extract the preconditioner for subksp object.
  call KSPGetPC(subksp, subpc, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! The subpc type will almost always be ILU
  call PCSetType(subpc, "ilu", ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Setup the matrix ordering for the subpc object:
  call PCFactorSetMatOrderingtype(subpc, "nd", ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Set the ILU parameters
  call PCFactorSetLevels(subpc, 1, ierr)
  call EChk(ierr, __FILE__, __LINE__) 

end subroutine setupPETScKSP

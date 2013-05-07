subroutine run3D(Xin, nNodes)

  use hypInput
  use hypData
  
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: Xin(3, nNodes)
  integer(kind=intType), intent(in) :: nNodes

  ! Working parameters
  integer(kind=intType) :: i, ierr, idim, l
  real(kind=realType), pointer :: xx(:)
  real(kind=realType) :: vBar

  ! Set the number of nodes on this proc, nx in hypData as nNodes
  nx = nNodes

  ! Write the header
  call writeHeader

  ! Create the PETSc Variables
  call create3DPetscVars

  ! Copy out Xin() into the first grid slot by setting a pointer and
  ! copying Xin

  call VecGetArrayF90(X(1), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do i=1,nx
     do idim=1,3
        xx((i-1)*3 + idim) = Xin(idim, i)
     end do
  end do
 
  call VecRestoreArrayF90(X(1), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute the inital 'Radius' value, R0
  call computeR(X(1), radius0)
  
  factorNext = .True.
  scaleDist = zero

  ! Determine starting CPU Time
  timeStart = mpi_wtime()

  ! This is the master marching direction loop

  marchLoop: do marchIter=2, Nmax

     ! Use 'l' as the working variable
     L = marchIter

     ! Compute the length increment in the marching direction. This is
     ! put in a separate function to allow for potential callbacks to
     ! user supplied functions if necessary

     call computeStretch(L)

     ! Increment our running counter on approx distance from wall
     scaleDist = scaleDist + deltaS

     ! Compute the max radius of the current level, X1
     call computeMinR(X(L), Radius)

     ! Possibly write header and write iteration info
     if (mod(marchIter, 50) == 0) then
        call writeHeader
     end if

     ! Compute volumes for current level
     call computeVolumes(X(L), Volume, VBar)

     ! Smooth the volume
     call volumeSmooth(Volume, vBar, L)

     ! Compute the initial guess for X(L+1) from the normal on X(L)
     call initialGuess(X(L-1), X(L))

     ! Solve for the next grid level
     call solveLevel(X(L-1), X(L), Volume)

     ! Shuffle the volumes backwards
     call VecCopy(volume, volumeOld, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Check the quality of this layer
     !call computeQualityLayer

     call writeIteration

  end do marchLoop


  ! Destroy the PETSc variables
  !call destroyPetscVars
     
 end subroutine run3D

subroutine computeStretch(l)

  use hypInput
  use hypData
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: l

  ! Since we don't have a complex stretching function yet, we will
  ! just use geometric progression which is usually close to what we
  ! actually want in the marching direction anyway

  if (l == 2) then
     ! First step, set deltaS to the desired initial grid spacign
     ! given in HypInput
     deltaS = pS0
  else
     ! Otherwise, just multiply the current deltaS by the pGridRatio
     ! parameter, also in HypInput

     deltaS = deltaS*pGridRatio
  end if

end subroutine computeStretch

subroutine computeVolumes(xVec, VVec, VBar)

  use precision
  use hypData, only : cRatio, nPatch, patches, nx, deltaS
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! Input parameters
  Vec, intent(in) :: xVec

  ! Output Parameters
  Vec, intent(out) :: vVec
  real(kind=realType), intent(out) :: VBar

  ! Working Variables
  integer(kind=intType) :: i, j, jp1, jm1, kp1, km1, nn, ipatch, nNeighbour, ii, ierr
  real(kind=realType) :: deltaV, ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS
  real(kind=realType) :: c1,c2,c3,c4, lenV1, lenV2
  real(kind=realType), dimension(:), pointer :: xx, vv

  ! Zero Volume vector and the normal vector
  call VecSet(vVec, zero, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the pointers  
  call VecGetArrayF90(xVec, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(vVec, vv, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  VBar = zero
  cratio = zero

  ! Loop over the patch faces:
  do iPatch=1,nPatch
     ! Loop over faces on patch
     do j=1,patches(iPatch)%jl-1
        do i=1,patches(iPatch)%il-1

           ! Extract the coordinates of the 4 corners
           ll = xx(3*patches(iPatch)%l_index(i  , j  )-2:patches(iPatch)%l_index(i  , j  ))
           ul = xx(3*patches(iPatch)%l_index(i  , j+1)-2:patches(iPatch)%l_index(i  , j+1))
           lr = xx(3*patches(iPatch)%l_index(i+1, j  )-2:patches(iPatch)%l_index(i+1, j  ))
           ur = xx(3*patches(iPatch)%l_index(i+1, j+1)-2:patches(iPatch)%l_index(i+1, j+1))

           ! Compute the area
           v1(:) = ur - ll
           v2(:) = ul - lr

           ! Cross Product
           s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
           s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
           s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

           ! Fourth of Area times deltaS 
           areadS = fourth*deltaS*half*sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))

           ! Scatter back to the nodes
           vv(patches(iPatch)%l_index(i, j)) = &
                vv(patches(iPatch)%l_index(i, j)) + areadS 
           vv(patches(iPatch)%l_index(i+1, j)) = &
                vv(patches(iPatch)%l_index(i+1, j)) + areadS 
           vv(patches(iPatch)%l_index(i, j+1)) = &
                vv(patches(iPatch)%l_index(i, j+1)) + areadS 
           vv(patches(iPatch)%l_index(i+1, j+1)) = &
                vv(patches(iPatch)%l_index(i+1, j+1)) + areadS 

           VBar = vBar + areadS*four


        end do
     end do
  end do

  VBar = VBar / nx

  ! Reset arrays
  call VecRestoreArrayF90(xVec, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(vVec, vv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine computeVolumes

subroutine volumeSmooth(vVec, vBar, l)

  use hypInput
  use hypData, only : nx, nPtr
  implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Perform averaging on the nodal volumes

  ! Input/Output Parameters
  Vec, intent(inout) :: vVec
  real(kind=realType), intent(in) :: VBar
  integer(kind=intType), intent(in) :: l

  ! Working Parameters
  real(kind=realType) :: factor, Vtmp(nx), oneOvrNNeighbour
  integer(kind=intType) :: i, ii, iter, ierr
  real(kind=realType), dimension(:), pointer :: volume

  ! Get array
  call VecGetArrayF90(vVec, volume, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Do a Jacobai volume smooth
  do i=1,nx
     vtmp(i) = volume(i)
  end do

  do iter=1,volSmoothIter
     do i=1,nx

        Volume(i) = zero
        oneOvrNNeighbour = one/nPtr(1,i)
        do ii=1,nPtr(1,i)
           ! Average around the neighbours
           Volume(i) = Volume(i) + oneOvrNNeighbour*& 
                ((one - volCoef) * Vtmp(i) + &
                volCoef*Vtmp(nPtr(1+ii,i)))
        end do
     end do

     ! Copy Volume back to Vtmp
     do i=1,nx
        vtmp(i) = volume(i)
     end do
  end do

  ! Now we do a global volume smooth
  factor = half * (one + ( one -volBlend)**(l-2))
  volume = factor*Volume + (one - factor)*Vbar

  ! Finally reset array
  call VecRestoreArrayF90(vVec, volume, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine volumeSmooth

subroutine initialGuess(X0, X1)

  use precision
  use hypData, only : nx, deltas, nPtr
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Take the vector in X0, evaluate the two in plane derivatives to
  ! get the normal derivative, then using deltaS estimate an algebraic
  ! estimate of the next layer of cell X1

  ! Input parameters
  Vec, intent(in) :: X0

  ! Output Parameter
  Vec, intent(out) :: X1

  ! Working Variables
  integer(kind=intType) :: i, ierr, jp1, jm1, kp1, km1, m, mm
  real(kind=realType), dimension(:), pointer :: xx0, xx1
  real(kind=realType), dimension(3) :: v1, v2, s, r0, rm
  real(kind=realType) :: thetam
  
  ! Set the pointers  
  call VecGetArrayF90(x0, xx0, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(X1, xx1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Nodal loop
  do i=1,nx
     if (nptr(1, i) == 4) then
        ! Extract pointers to make code easier to read

        jp1 = nPtr(2, i)
        kp1 = nPtr(3, i)
        jm1 = nPtr(4, i)
        km1 = nPtr(5, i)

        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        v1 = (xx0(3*jp1-2:3*jp1) - xx0(3*jm1-2:3*jm1))
        v2 = (xx0(3*kp1-2:3*kp1) - xx0(3*km1-2:3*km1))
      
     else
        v1 = zero
        v2 = zero

        r0 = xx0(3*i-2:3*i)
        MM = nPtr(1, i) 

        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = xx0(3*nptr(2 + m, i)-2:3*nptr(2 + m,i))
           v1 = v1 + (two/MM) * (rm - r0) * cos(thetam)
           v2 = v2 + (two/MM) * (rm - r0) * sin(thetam)
        end do
     end if

     ! Compute the cross product of the in-place derivative vectors
     s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     s(3) = (v1(1)*v2(2) - v1(2)*v2(1))
     
     ! Compute the next layer position by "extruding" a distance of
     ! 'deltaS' along 's'
     xx1(3*i-2:3*i) = xx0(3*i-2:3*i) + deltaS*s
     
  end do
  
  ! reset the pointers  
  call VecRestoreArrayF90(x0, xx0, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(X1, xx1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine initialGuess

! subroutine assembleAndSolve3D(Volume1, nx, l)

!   use hypInput
!   use hypData

!   implicit none

!   ! Input Parameters
!   real(kind=realType), intent(in) :: volume1(nx)
!   integer(kind=intType), intent(in) :: nx, l

!   ! Working parameters
!   real(kind=realType) :: A0(3, 3), B0(3, 3), C0(3, 3)
!   real(Kind=realType) :: A00(3,3)
!   integer(kind=intType) :: j, idim, ierr
!   integer(kind=intType) :: i, jp1, jm1, jm2, jm3, kp1, km1, km2, km3, nAverage, ii, m, MM,  iter
!   real(kind=realType) :: eye(3, 3), CInv(3,3), CInvA(3,3), CInvB(3, 3)
!   real(kind=realType) :: De(3), tmp, tmp1, tmp2,a_ksi, a_eta
!   real(kind=realType) :: exp, Rksi, Reta, N_ksi, N_eta, deltaR(3)
!   real(kind=realType) :: r0(3), rm(3), thetam, coef, d_ksi, d_eta, dmax
!   real(kind=realType) :: r_ksi(3), r_eta(3), r_zeta(3), v1(3), v2(3), s(3), detC
!   real(kind=realType) :: deltaVovrDetC, oneOvrNNeighbour, factor, fact
!   real(kind=realType) :: rplusj(3), rminusj(3), rplusk(3), rminusk(3)
!   real(Kind=realType) :: rplusjhat(3), rpluskhat(3), nhat(3), alpha, beta
!   real(kind=realType) :: amax, ovrSum1, ovrSum2
!   ! Define a 3x3 'I' matrix
!   eye = zero
!   eye(1,1) = one
!   eye(2,2) = one
!   eye(3,3) = one
!   call MatZeroEntries(hypMat, ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   call VecZeroEntries(hypRHS, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   ! ------------------------------------
!   ! Scaling Function (Equation 6.5) -- Modified
!   ! ------------------------------------

!   ! Since we don't acutally know how many iterations it will take to
!   ! get to the farfield, we really don't want to use an N scaling. Sl
!   ! is really quite important, but only near the boundary layer such
!   ! that grid orthogonality is maintained. It is very important that
!   ! Sl continues

!   sl = (scaleDist/(radius0*rMin))**slExp

!   ! Record the maximum and minimum grid sensors
!   gridSensorMax = zero 
!   gridSensorMin = huge(gridSensorMin)

!   amax = zero
!   ! Next we assemble hypMat which is quite straight forward
!   do i=1, nx

!      if (nptr(1, i) == 4) then
!         ! Extract pointers to make code easier to read

!         jp1 = nPtr(2, i)
!         kp1 = nPtr(3, i)
!         jm1 = nPtr(4, i)
!         km1 = nPtr(5, i)

!         ! Compute centered difference with respect to ksi and eta. This is the
!         ! nominal calculation that we will use in the interior of the domain
!         r_ksi = half*(X0(:, jp1) - X0(:, jm1))
!         r_eta  = half*(X0(:, kp1) - X0(:, km1))

!         ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
!         tmp = (dist(Xm1(:, jp1), Xm1(:, i)) + dist(Xm1(:, jm1), Xm1(:, i))) / &
!              (dist(X0(:, jp1), X0(:, i)) + dist(X0(:, jm1), X0(:, i)))
!         d_ksi = max(tmp**(two/sl), 0.1)

!         tmp = (dist(Xm1(:, kp1), Xm1(:, i)) + dist(Xm1(:, km1), Xm1(:, i))) / &
!              (dist(X0(:, kp1), X0(:, i)) + dist(X0(:, km1), X0(:, i)))
!         d_eta = max(tmp**(two/sl), 0.1)

!         gridSensorMax = max(gridSensorMax, d_ksi, d_eta)
!         gridSensorMin = min(gridSensorMin, d_ksi, d_eta)
!         ! Equation 3.5 (Chen and Steger)
!         detC = &
!              (r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))**2 + &
!              (r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))**2 + &
!              (r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))**2

!         deltaVovrDetC = Volume1(i)/detC

!         r_zeta(1) = deltaVovrDetC*(r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))
!         r_zeta(2) = deltaVovrDetC*(r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))
!         r_zeta(3) = deltaVovrDetC*(r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))

!         ! Compute the two coefficient matrices
!         A0(1, 1) = r_zeta(1)
!         A0(1, 2) = r_zeta(2)
!         A0(1, 3) = r_zeta(3)
!         A0(2, 1) = zero
!         A0(2, 2) = zero
!         A0(2, 3) = zero
!         A0(3, 1) = r_eta(2)*r_zeta(3) - r_zeta(2)*r_eta(3)
!         A0(3, 2) = r_zeta(1)*r_eta(3) - r_eta(1)*r_zeta(3)
!         A0(3, 3) = r_eta(1)*r_zeta(2) - r_zeta(1)*r_eta(2)

!         B0(1, :) = zero
!         B0(2, :) = r_zeta
!         B0(3, 1) = r_zeta(2)*r_ksi(3) - r_ksi(2)*r_zeta(3)
!         B0(3, 2) = r_ksi(1)*r_zeta(3) - r_zeta(1)*r_ksi(3)
!         B0(3, 3) = r_zeta(1)*r_ksi(2) - r_ksi(1)*r_zeta(2)

!         C0(1, :) = r_ksi
!         C0(2, :) = r_eta
!         C0(3, 1) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
!         C0(3, 2) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
!         C0(3, 3) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)

!         ! Compute the inverse of C and its multiplcation with A and B
!         call three_by_three_inverse(C0, Cinv)
!         CinvA = matmul(Cinv, A0)
!         CinvB = matmul(Cinv, B0)

!         ! -----------------
!         ! Left Hand Side:
!         ! -----------------

!         ! Point to the left
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm1-1, &
!              (-(one+theta)*half*CinvA - epsI*eye), ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Center Point ! 
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!              (eye + four*epsI*eye), ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Point to the right
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, jp1-1, &
!              ((one+theta)*half*CinvA - epsI*eye),  ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! ----------------------------------------------------------

!         ! Point to the bottom
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, km1-1, &
!              (-(one+theta)*half*CinvB - epsI*eye), ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Point to the right
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, kp1-1, &
!              ((one+theta)*half*CinvB - epsI*eye),  ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! -----------------
!         ! Right Hand Side:
!         ! -----------------

!         ! Assemble the the Cinv comonent of the RHS
!         call VecSetValuesBlocked(hypRHS, 1, i-1, &
!              matmul(Cinv,(/zero, zero, Volume1(i)/)), ADD_VALUES, ierr)

!         call EChk(ierr, __FILE__, __LINE__)

!         ! ============================================
!         !             Explicit Smoothing 
!         ! ============================================

!         ! ------------------------------------
!         ! Compute the grid angle functions (Equation 6.9)
!         ! ------------------------------------
!         rplusj = X0(:, jp1) - X0(:, i)
!         rminusj = X0(:, jm1) -X0(:, i)

!         rplusk = X0(:, kp1) - X0(:, i)
!         rminusk = X0(:, km1) -X0(:, i)

!         v1 = rplusj - rminusj
!         v2 = rplusk - rminusk

!         ! Take cross product
!         nhat(1) = (v1(2)*v2(3) - v1(3)*v2(2))
!         nhat(2) = (v1(3)*v2(1) - v1(1)*v2(3))
!         nhat(3) = (v1(1)*v2(2) - v1(2)*v2(1))

!         ! And normalize
!         nhat = nhat /sqrt(nhat(1)**2 + nhat(2)**2 + nhat(3)**2)

!         ! Normalize rplusj 
!         rplusjhat = rplusj /sqrt(rplusj(1)**2 + rplusj(2)**2 + rplusj(3)**2)

!         alpha = acos(dot_product(nhat,rplusjhat))
!         if (alpha >= 0 .and. alpha < half*pi) then
!            a_ksi = one/(one-cos(alpha)**2)
!         else
!            a_ksi = one
!         end if

!         ! Normalize rplusk 
!         rpluskhat = rplusk /sqrt(rplusk(1)**2 + rplusk(2)**2 + rplusk(3)**2)

!         beta = acos(dot_product(nhat,rplusjhat))
!         if (beta >= 0 .and. beta < half*pi) then
!            a_eta = one/(one-cos(beta)**2)
!         else
!            a_eta = one
!         end if

!         ! ------------------------------------
!         ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
!         ! ------------------------------------
!         Rksi = Sl * d_ksi * a_ksi
!         Reta = Sl * d_eta * a_eta

!         ! ------------------------------------
!         ! Matrix Norm approximation  (Equation 6.3)
!         ! ------------------------------------
!         N_ksi = sqrt((A0(1, 1)**2 + A0(1, 2)**2 + A0(1, 3)**2)/ &
!              (C0(1 ,1)**2 + C0(1, 2)**2 + C0(1, 3)**2))
!         N_eta = sqrt((A0(1, 1)**2 + A0(1, 2)**2 + A0(1, 3)**2)/ &
!              (C0(2 ,1)**2 + C0(2, 2)**2 + C0(2, 3)**2))

!         ! ------------------------------------
!         ! Final explict dissipation (Equation 6.1 and 6.2)
!         ! ------------------------------------

!         De = &
!              epsE*Rksi*N_ksi * (X0(:, jm1) - two*X0(:, i) + X0(:, jp1)) + &
!              epsE*Reta*N_eta * (X0(:, km1) - two*X0(:, i) + X0(:, kp1))

!         if (l == 2) then
!            De = zero
!         end if

!         call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)
!         ! -----------------------------------------------------------------------------------
!      else if (nptr(1, i) .ne. 4) then
!         ! -----------------------------------------------------------------------------------
!         ! We have a node with 6 neighbours.
!         ! Lets do this brute force:
!         r_ksi = zero
!         r_eta = zero
!         r0 = X0(:, i)
!         MM = nPtr(1, i) 

!         tmp1 = zero
!         tmp2 = zero
!         ovrSum1 = zero
!         ovrSum2 = zero
!         do m = 0, MM - 1
!            thetam = 2*pi*m/MM
!            rm = X0(:, nptr(2 + m, i))
!            r_ksi = r_ksi + (two/MM) * (rm - r0) * cos(thetam)
!            r_eta = r_eta + (two/MM) * (rm - r0) * sin(thetam)

!            tmp1 = tmp1 + &
!                 dist(Xm1(:, nptr(2+m, i)), Xm1(:, i)) / &
!                 dist( X0(:, nptr(2+m, i)),  X0(:, i))*abs(cos(thetam))

!            tmp2 = tmp2 + & 
!                 dist(Xm1(:, nptr(2+m, i)), Xm1(:, i)) / &
!                 dist( X0(:, nptr(2+m, i)),  X0(:, i))*abs(sin(thetam))
!            ovrSum1 = ovrSum1 + abs(cos(thetam))
!            ovrSum2 = ovrSum2 + abs(sin(thetam))
!         end do

!         d_ksi = max((tmp1/ovrSum1)**(two/sl), 0.1)
!         d_eta = max((tmp2/ovrSum2)**(two/sl), 0.1)

!         gridSensorMax = max(gridSensorMax, d_ksi, d_eta)
!         gridSensorMin = min(gridSensorMin, d_ksi, d_eta)

!         ! Equation 3.5 (Chen and Steger)
!         detC = &
!              (r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))**2 + &
!              (r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))**2 + &
!              (r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))**2

!         deltaVovrDetC = Volume1(i)/detC

!         r_zeta(1) = deltaVovrDetC*(r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))
!         r_zeta(2) = deltaVovrDetC*(r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))
!         r_zeta(3) = deltaVovrDetC*(r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))

!         ! Compute the two coefficient matrices
!         A0(1, 1) = r_zeta(1)
!         A0(1, 2) = r_zeta(2)
!         A0(1, 3) = r_zeta(3)
!         A0(2, 1) = zero
!         A0(2, 2) = zero
!         A0(2, 3) = zero
!         A0(3, 1) = r_eta(2)*r_zeta(3) - r_zeta(2)*r_eta(3)
!         A0(3, 2) = r_zeta(1)*r_eta(3) - r_eta(1)*r_zeta(3)
!         A0(3, 3) = r_eta(1)*r_zeta(2) - r_zeta(1)*r_eta(2)

!         B0(1, :) = zero
!         B0(2, :) = r_zeta
!         B0(3, 1) = r_zeta(2)*r_ksi(3) - r_ksi(2)*r_zeta(3)
!         B0(3, 2) = r_ksi(1)*r_zeta(3) - r_zeta(1)*r_ksi(3)
!         B0(3, 3) = r_zeta(1)*r_ksi(2) - r_ksi(1)*r_zeta(2)

!         C0(1, :) = r_ksi
!         C0(2, :) = r_eta
!         C0(3, 1) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
!         C0(3, 2) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
!         C0(3, 3) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)
!         !Setting the matrix values are a little tricker here now. 

!         ! (This is the same)
!         call three_by_three_inverse(C0, Cinv)
!         CinvA = matmul(Cinv, A0)
!         CinvB = matmul(Cinv, B0)
!         MM = nPtr(1, i) 
!         do m = 0, MM - 1
!            thetam = 2*pi*m/MM

!            ! For the 'fm' point in ksi
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
!                 (one + theta)*CinvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'f0' point in ksi
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!                 -(one + theta)*cInvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'fm' point in eta
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
!                 (one + theta)*CinvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'f0' point in ksi
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!                 -(one + theta)*cInvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! Now we have the second derivative values to do in ksi-ksi

!            ! ! For the 'fm' point in ksi-ksi 
!            coef = (two/MM)*(four*cos(thetam)**2 - one)

!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
!                 -eye*epsI*coef, ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'f0' point in ksi-ksi
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!                 eye*epsI*coef, ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'fm' point in eta-eta
!            coef = (two/MM)*(four*sin(thetam)**2 - one)

!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
!                 -eye*epsI*coef, ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)

!            ! For the 'f0' point in eta-eta
!            call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!                 eye*epsI*coef, ADD_VALUES, ierr)
!            call EChk(ierr, __FILE__, __LINE__)
!         end do

!         ! Finally we need an eye on the diagonal
!         call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!              eye, ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! -----------------
!         ! Right Hand Side:
!         ! -----------------

!         ! Assemble the the Cinv comonent of the RHS
!         call VecSetValuesBlocked(hypRHS, 1, i-1, &
!              matmul(Cinv,(/zero, zero, Volume1(i)/)), ADD_VALUES, ierr)

!         call EChk(ierr, __FILE__, __LINE__)

!         ! ============================================
!         !             Explicit Smoothing 
!         ! ============================================

!         ! ------------------------------------
!         ! Compute the grid angle functions (Equation 6.9)
!         ! ------------------------------------
!         a_ksi = one
!         a_eta = one

!         ! ------------------------------------
!         ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
!         ! ------------------------------------
!         Rksi = Sl * d_ksi * a_ksi
!         Reta = Sl * d_eta * a_eta

!         ! ------------------------------------
!         ! Matrix Norm approximation  (Equation 6.3)
!         ! ------------------------------------
!         N_ksi = sqrt((A0(1, 1)**2 + A0(1, 2)**2 + A0(1, 3)**2)/ &
!              (C0(1 ,1)**2 + C0(1, 2)**2 + C0(1, 3)**2))
!         N_eta = sqrt((A0(1, 1)**2 + A0(1, 2)**2 + A0(1, 3)**2)/ &
!              (C0(2 ,1)**2 + C0(2, 2)**2 + C0(2, 3)**2))

!         ! ------------------------------------
!         ! Final explict dissipation (Equation 6.1 and 6.2)
!         ! ------------------------------------
!         de = zero
!         r0 = X0(:, i)
!         MM = nPtr(1, i) 
!         do m = 0, MM - 1
!            thetam = 2*pi*m/MM
!            rm = X0(:, nptr(2 + m, i))

!            De = De + &
!                 epsE*Rksi*N_ksi*(two/MM)*(rm - r0)*(four*cos(thetam)**2 - one) + &
!                 epsE*Reta*N_eta*(two/MM)*(rm - r0)*(four*sin(thetam)**2 - one) 
!         end do

!         if (l == 2) then
!            De = zero
!         end if

!         call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!      else
!         ! Two!
!         stop

!         !  !  !We have to do averaging at the corners

!         !  ! ! ! Get the number of nodes to average:
!         !  nAverage = nPtr(1, i) 

!         !  ! ! Loop over neighbours
!         !  do ii = 1, nAverage
!         !     call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(1+ii, i)-1, &
!         !          -one/nAverage*eye, ADD_VALUES, ierr)
!         !  end do

!         ! ! Center Point
!         !  call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!         !       eye, ADD_VALUES, ierr)
!         !  call EChk(ierr, __FILE__, __LINE__)

!         ! !  !Don't set a RHS value for this, it is zero. This is already
!         ! !  !taken care of with the vecSet above
!      end if
!   end do

!   ! Assemble the matrix and vector
!   call MatAssemblyBegin(hypMat, MAT_FINAL_ASSEMBLY, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call MatAssemblyEnd(hypMat, MAT_FINAL_ASSEMBLY, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call VecAssemblyBegin(hypRHS, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call VecAssemblyEnd(hypRHS, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   ! Set the operator so PETSc knows to refactor
!   if (mod(l,preConLag) == 0 .or. factorNext) then
!      call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
!      call EChk(ierr, __FILE__, __LINE__)
!      factorNext = .False.
!   end if

!   ! Now solve the system
!   call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call KSPGetIterationNumber(hypKsp, kspIts, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   ! If we did more iterations that the size of the subspace, flag the
!   ! next iteration to refactor the PC.
!   if (kspIts > kspsubspacesize) then
!      factorNext = .True.
!   end if

!   ! Take the solution in deltaR and add to X0 to get X1
!   do i=1,nx
!      call VecGetValues(hypDelta, 3, (/3*(i-1), 3*(i-1)+1, 3*(i-1)+2/), &
!           deltaR, ierr)
!      call EChk(ierr, __FILE__, __LINE__)
!      X1(:, i) = X0(:, i) + deltaR
!   end do

! contains 

!   function dist(p1, p2)

!     real(kind=realType) :: p1(3), p2(3)
!     real(kind=realType) :: dist
!     dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)

!   end function dist
! end subroutine assembleAndSolve3D

subroutine calcResidual(snes, xVec, rVec, ctx, ierr)

  use precision

  ! PETSc Variables
  SNES snes
  Vec xVec, rVec
  PetscFortranAddr ctx(*)
  integer(kind=intType), intent(out) :: ierr

  ! This routine must do everything required to evaluate the full
  ! non-linear residual

  ierr = 0

end subroutine calcResidual

subroutine calcJacobian(snes, xVec, Jac, B, flag, ctx, ierr)

  use hypInput
  implicit none
#include "include/finclude/petsc.h"

  ! Input (Petsc) Variables
  SNES snes
  Vec xVec
  Mat Jac, B
  PetscInt flag
  PetscFortranAddr ctx(*)
  integer(kind=intType) :: ierr

  ierr = 0

end subroutine calcJacobian

subroutine create3DPetscVars

  use hypInput
  use hypData

  implicit none
  
  ! Working Variables
  integer(kind=intType) :: onProc(nx*3), offProc(nx*3), ierr, i, idim
  external calcResidual, calcJacobian
  real(kind=realType),pointer :: xx(:)

  ! ----------------------------------------------------------
  !          Linearized Hyperbolic System Variables
  ! ----------------------------------------------------------

  ! Lets to things in the proper way, create the Mat first
  onProc = 35
  offProc = 1

  ! Create a blocked matrix
  call MatCreateMPIBAIJ(PETSC_COMM_WORLD, 3, &
       nx*3, nx*3, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! This must be set to allow passing in blocks in native fortran
  ! ordering
  call MatSetOption(hypMat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Then use getVecs to get the vectors we want
  call MatGetVecs(hypMat, hypDelta, hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the extra state-sized vectors
  call VecDuplicate(hypDelta, hypRes, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(hypDelta, Volume, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(hypDelta, VolumeOld, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the full list of grid vectors:
  allocate(X(N), stat=ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Loop over each "vec" in the X array:
  do i=1,N
     call VecDuplicate(hypDelta, X(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

  ! Create the SNES object
  call SNESCreate(PETSC_COMM_WORLD, hypSNES, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Most important thing is to set the "function" that computes the residual "R"
  call SNESSetFunction(hypSNES, hypRes, calcResidual, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Next we must set the jacobian and the function to evaluate the jacobian
  call SNESSetJacobian(hypSNES, hypMat, hypMat, calcJacobian, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the sens tolerances...some of these should be options eventually
  call SNESSetTolerances(hypSNES, 1e-16, kspRelTol, 1e-12, 100, 1000, ierr)
  call EChk(ierr, __FILE__, __LINE__)

    


end subroutine create3DPetscVars

 !  factorNext = .True.
 !  scaleDist = zero
 !  ! Determine starting CPU Time
 !  timeStart = mpi_wtime()
 !  l_0 = -1
 !  smoothIter = 0
 !  ! This is the master marching direction loop

 !  marchLoop: do marchIter=2, Nmax

 !     ! Use 'l' as the working variable
 !     l = marchIter

 !     ! Set pointers to the two levels we are currently working with
 !     X0 => pGrid3D(:, :, l-1)
 !     X1 => pGrid3D(:, :, l  )
 !     if (l > 2) then
 !        Xm1 => pGrid3D(:, :, l-2)
 !     end if

 !     ! Compute the length increment in the marching direction. This is
 !     ! put in a separate function to allow for potential callbacks to
 !     ! user supplied functions if necessary

 !     call computeStretch(l)

 !     ! Compute the nodal volumes on the X0 layer as well as the
 !     ! average volume VVar
 !     call computeVolumes(Volume1, nx, Vbar)

 !     if (cratio > cmax) then 
 !        deltaS = deltas*.5
 !        call computeVolumes(Volume1, nx, Vbar)
 !     end if

 !     ! Increment our running counter on distance from wall
 !     scaleDist = scaleDist + deltaS

 !     ! Now perform the "Volume" smoothing operation.
 !     call VolumeSmooth(Volume1, nx, VBar, l)

 !     ! Assemble and solve
 !     call assembleAndSolve3d(Volume1, nx, l)

 !     ! Compute the max radius of the current level, X1
 !     call computeMinR(nx,Radius)

 !     ! Shuffle the volumes's backwards
 !     do i=1,nx
 !        Volume0(i) = Volume1(i)
 !     end do

 !     ! Possibly write header and write iteration info
 !     if (mod(marchIter, 50) == 0) then
 !        call writeHeader
 !     end if

 !     ! Check the quality of this layer
 !     call computeQualityLayer

 !     call writeIteration

 !     ! Now check if we have gone far enough:
 !     if (scaleDist/radius0 > rmin) then
 !        nLayers = l
 !        exit marchLoop
 !     end if

 !     if (l == Nmax) then
 !        Nlayers = l
 !     end if

 !  end do marchLoop
















! ! We have to do these special nodes:
! ! Loop over the number of nodes
! fact = one/dble(nPtr(1,i))
! do ii =1, nPtr(1, i)
!    jm1 = nPtr((ii-1)*6+2, i)
!    jm2 = nPtr((ii-1)*6+3, i)
!    jm3 = nPtr((ii-1)*6+4, i)
!    km1 = nPtr((ii-1)*6+5, i)
!    km2 = nPtr((ii-1)*6+6, i)
!    km3 = nPtr((ii-1)*6+7, i)

!    ! Compute one sided difference wrt to ksi and eta
!    r_ksi = -X0(:, i) + X0(:, jm1) 
!    r_eta = -X0(:, i) + X0(:, km1) 

!    ! r_ksi = (-three/two)*X0(:, i) + two*X0(:, jm1) - half*X0(:, jm2)
!    ! r_eta = (-three/two)*X0(:, i) + two*X0(:, km1) - half*X0(:, km2)


!   ! 3.5 (Chen and Steger)
!    detC = &
!         (r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))**2 + &
!         (r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))**2 + &
!         (r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))**2

!    deltaVovrDetC = Volume1(i)/detC

!    r_zeta(1) = deltaVovrDetC*(r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))
!    r_zeta(2) = deltaVovrDetC*(r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))
!    r_zeta(3) = deltaVovrDetC*(r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))

!    ! Compute the two coefficient matrices
!    A0(1, :, i) = r_zeta
!    A0(2, :, i) = zero
!    A0(3, 1, i) = r_eta(2)*r_zeta(3) - r_zeta(2)*r_eta(3)
!    A0(3, 2, i) = r_zeta(1)*r_eta(3) - r_eta(1)*r_zeta(3)
!    A0(3, 3, i) = r_eta(1)*r_zeta(2) - r_zeta(1)*r_eta(2)

!    B0(1, :, i) = zero
!    B0(2, :, i) = r_zeta
!    B0(3, 1, i) = r_zeta(2)*r_ksi(3) - r_ksi(2)*r_zeta(3)
!    B0(3, 2, i) = r_ksi(1)*r_zeta(3) - r_zeta(1)*r_ksi(3)
!    B0(3, 3, i) = r_zeta(1)*r_ksi(2) - r_ksi(1)*r_zeta(2)

!    C0(1, :, i) = r_ksi
!    C0(2, :, i) = r_eta
!    C0(3, 1, i) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
!    C0(3, 2, i) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
!    C0(3, 3, i) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)

!    call three_by_three_inverse(C0(:, :, i), Cinv)
!    CinvA = matmul(Cinv, A0(:, :, i))
!    CinvB = matmul(Cinv, B0(:, :, i))

!    ! Now set the matrix
!    ! -----------------
!    ! Left Hand Side:
!    ! -----------------

!    ! ! Center Cell
!    ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!    !      fact*(eye + (one+theta)*cInvA*(-three/two) + &
!    !                  (one+theta)*cInvB*(-three/two) + & 
!    !                  two*epsI*eye), ADD_VALUES, ierr)
!    ! call EChk(ierr, __FILE__, __LINE__)

!    ! ! jm1
!    ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm1-1, &
!    !      fact*((one+theta)*cInvA*two - & 
!    !      two*epsI*eye), ADD_VALUES, ierr)
!    ! call EChk(ierr, __FILE__, __LINE__)

!    ! ! jm2
!    ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm2-1, &
!    !      fact*((one+theta)*cInvA*(-half) + &
!    !      epsI*eye), ADD_VALUES, ierr)
!    ! call EChk(ierr, __FILE__, __LINE__)

!    ! ! km1
!    ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, km1-1, &
!    !      fact*((one+theta)*cInvB*two - &
!    !      two*epsI*eye), ADD_VALUES, ierr)
!    ! call EChk(ierr, __FILE__, __LINE__)

!    ! ! jm2
!    ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, km2-1, &
!    !      fact*((one+theta)*cInvB*(-half) + &
!    !      epsI*eye), ADD_VALUES, ierr)
!    ! call EChk(ierr, __FILE__, __LINE__)

!    ! Center Cell
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!         fact*(eye + (one+theta)*cInvA*(-one) + &
!                     (one+theta)*cInvB*(-one) - & 
!                     two*epsI*eye), ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)

!    ! jm1
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm1-1, &
!         fact*((one+theta)*cInvA*one + & 
!         two*epsI*eye), ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)

!    ! jm2
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm2-1, &
!         -fact*(epsI*eye), ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)

!    ! km1
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, km1-1, &
!         fact*((one+theta)*cInvB*one + &
!         two*epsI*eye), ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)

!    ! jm2
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, km2-1, &
!         -fact*(epsI*eye), ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)



!    ! -----------------
!    ! Right Hand Side:
!    ! -----------------

!    ! Assemble the the Cinv comonent of the RHS
!    call VecSetValuesBlocked(hypRHS, 1, i-1, &
!         matmul(Cinv,(/zero, zero, Volume1(i)/))*fact, ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)

!    ! ============================================
!    !             Explicit Smoothing 
!    ! ============================================
!    a_ksi = one
!    a_eta = one

!    ! ------------------------------------
!    ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
!    ! ------------------------------------
!    ! Compute centered difference with respect to ksi and eta. This is the
!    ! nominal calculation that we will use in the interior of the domain

!    ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
!    ! tmp = dist(Xm1(:, i), Xm1(:, jm1))/ dist(X0(:, i), X0(:, jm1))
!    ! d_ksi = max(tmp**(two/sl), 0.1)

!    ! tmp = dist(Xm1(:, i), Xm1(:, km1))/ dist(X0(:, i), X0(:, km1))
!    ! d_eta = max(tmp**(two/sl), 0.1)

!    Rksi = Sl * d_ksi * a_ksi
!    Reta = Sl * d_eta * a_eta

!    ! ------------------------------------
!    ! Matrix Norm approximation  (Equation 6.3)
!    ! ------------------------------------
!    N_ksi = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
!         (C0(1 ,1 ,i)**2 + C0(1, 2, i)**2 + C0(1, 3, i)**2))
!    N_eta = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
!         (C0(2 ,1 ,i)**2 + C0(2, 2, i)**2 + C0(2, 3, i)**2))

!    ! ------------------------------------
!    ! Final explict dissipation (Equation 6.1 and 6.2)
!    ! ------------------------------------

!    De = &
!         epsE*Rksi*N_ksi * (X0(:, i) - two*X0(:, jm1) + one*X0(:, jm2)) + &
!         epsE*Reta*N_eta * (X0(:, i) - two*X0(:, km1) + one*X0(:, km2))

!     ! De = epsE*Rksi*N_ksi * (two*X0(:, i) - five*X0(:, jm1) + four*X0(:, jm2) - X0(:, jm3)) + &
!     !      epsE*Reta*N_eta * (two*X0(:, i) - five*X0(:, km1) + four*X0(:, km2) - X0(:, km3))

!    if (l == 2) then
!       De = zero
!    end if

!    call VecSetValuesBlocked(hypRHS, 1, i-1, fact*De, ADD_VALUES, ierr)
!    call EChk(ierr, __FILE__, __LINE__)


! ! Laplacian blending
! nAverage = nPtr(1, i) 

! ! ! Loop over neighbours
! do ii = 1, nAverage
!    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(1+ii, i)-1, &
!         -one/nAverage*eye*(one-fact), ADD_VALUES, ierr)
! end do

! ! Center Point
! call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
!      eye*(one-fact), ADD_VALUES, ierr)
! call EChk(ierr, __FILE__, __LINE__)




! if (scaleDist > 4.0) then
!    print *,'laplacian smoothing'
!    ! Do laplacian smoothing
!    Xtmp = X1
!    if (smoothIter < 5) then
!       smoothIter = smoothIter + 1
!    end if
!    ! ALso ramp the factor
!    factor = (smoothIter-1)*0.25/4.0
!    do iter=1,smoothIter
!       do i=1,nx
!          X1(:,i) = zero
!          oneOvrNNeighbour = one/nPtr(1,i)
!          do ii=1,nPtr(1,i)
!             ! Average around the neighbours
!             X1(:,i) = X1(:,i) + oneOvrNNeighbour*& 
!                  ((one - factor) * Xtmp(:,i) + &
!                  factor*Xtmp(:,nPtr(1+ii,i)))
!          end do
!       end do
!       ! Copy X1 backto Xtmp
!       Xtmp = X1
!    end do
! end if


  ! ! Loop over the nodes
  ! volume = zero
  ! do i=1,nx
  !    if (nPtr(1,i) == 4) then
  !       jp1 = nPtr(2, i)
  !       kp1 = nPtr(3, i)
  !       jm1 = nPtr(4, i)
  !       km1 = nPtr(5, i)


  !       v1 = X(:,jp1)-X(:,jm1)
  !       v2 = X(:,kp1)-X(:,km1)

  !       s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
  !       s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
  !       s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

  !       ! Get Length
  !       lenV1 = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
  !       lenV2 = sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2)

  !       volume(i) = fourth*deltaS*half*sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
  !       vBar = vBar + volume(i)

  !       c1 = two*deltaS/lenV1
  !       c2 = two*deltaS/lenV2

  !       cmax = max(cmax,c1,c2)
  !    end if
  ! end do

  ! do i=1,nx
  !    if (nPtr(1,i) .ne.4) then
  !       nNeighbour = nPtr(1, i) 
  !       volume(i) = zero
  !       do ii = 1,nNeighbour
  !          volume(i) = (one/nNeighbour)*volume(nPtr(1+ii,i))
  !       end do
  !       vBar = vBar + volume(i)
  !    end if
  ! end do
  ! nodalVolume = volume


subroutine run3D(Xin, nNodes)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: run3D is the main python interface to generate the
  !     3D hyperbolic mesh.
  !    
  !     Description of Arguments
  !     Input:
  !     Xin - array size (nNodes, 3): The (unstructured) list of nodes that 
  !           this processor owns. 
  !     nNodes - integer : The number of nodes
  !
  !     Ouput: None
  !
  use precision
  use hypData
  use hypInput

  
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: Xin(3, nNodes)
  integer(kind=intType), intent(in) :: nNodes

  ! Working parameters
  integer(kind=intType) :: i, ierr, idim, l
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

  marchLoop: do marchIter=2, N

     ! Use 'l' as the working variable
     L = marchIter

     ! Compute the length increment in the marching direction. This is
     ! put in a separate function to allow for potential callbacks to
     ! user supplied functions if necessary

     call computeStretch(L)

     ! Increment our running counter on approx distance from wall
     scaleDist = scaleDist + deltaS

     ! Possibly write header and write iteration info
     if (mod(marchIter, 50) == 0) then
        call writeHeader
     end if

     ! Compute volumes for current level
     call computeVolumes(X(L-1), Volume, VBar)

     ! Smooth the volume
     call volumeSmooth(Volume, vBar, L)

     ! Compute the initial guess for X(L) from the normals on X(L-1)
     call initialGuess(X(L-1), X(L))

     ! Solve for the next grid level
     if (L==2) then
         ! We don't have two history points yet...
        call solveLevel(X(L-1), X(L-1), X(L), Volume)
     else
        call solveLevel(X(L-2), X(L-1), X(L), Volume)
     end if

     ! Compute the max radius of the current level, X(L)
     call computeMinR(X(L), Radius)

     ! Check the quality of this layer
     !call computeQualityLayer

     call writeIteration

  end do marchLoop

  ! Destroy the PETSc variables
  !call destroyPetscVars
     
 end subroutine run3D

subroutine computeStretch(L)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: computeStretch determines the global stretch value,
  !     deltaS, depending on the grid level L. In the future
  !     this function may be more complex or call a user-supplied function.
  !    
  !     Description of Arguments
  !     Input:
  !     L - integer: Marching 
  !
  !     Ouput:
  !     deltaS: Marching increment set in hypData
  !
  use precision
  use hypInput
  use hypData
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: L

  ! Since we don't have a complex stretching function yet, we will
  ! just use geometric progression which is usually close to what we
  ! actually want in the marching direction anyway

  if (L == 2) then
     ! First step, set deltaS to the desired initial grid spacign
     ! given in HypInput
     deltaS = s0
  else
     ! Otherwise, just multiply the current deltaS by the pGridRatio
     ! parameter, also in HypInput
     deltaS = deltaS*gridRatio
  end if

end subroutine computeStretch

subroutine computeVolumes(xVec, VVec, VBar)
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
  !     xVec - Petsc Vector: Vector of nodes to use
  !
  !     Ouput:
  !     vVec - Petsc Vector: Vector of volumes, stored in the third DOF of vVec
  !     VVar - real : Average nodal volume

  use precision
  use hypData, only : cRatio, nPatch, patches, nx, deltaS, xx, nptr
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input parameters
  Vec, intent(in) :: xVec

  ! Output Parameters
  Vec, intent(out) :: vVec
  real(kind=realType), intent(out) :: VBar

  ! Working Variables
  integer(kind=intType) :: i, j, ipatch, ierr
  real(kind=realType) :: ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS
  real(kind=realType), dimension(:), pointer :: vv

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

           ! Scatter back to the nodes. Put in the third 3dof slot of
           ! vv. The other two slots are unsued. 
           vv(patches(iPatch)%l_index(i, j)*3) = &
                vv(patches(iPatch)%l_index(i, j)*3) + areadS 
           vv(patches(iPatch)%l_index(i+1, j)*3) = &
                vv(patches(iPatch)%l_index(i+1, j)*3) + areadS 
           vv(patches(iPatch)%l_index(i, j+1)*3) = &
                vv(patches(iPatch)%l_index(i, j+1)*3) + areadS 
           vv(patches(iPatch)%l_index(i+1, j+1)*3) = &
                vv(patches(iPatch)%l_index(i+1, j+1)*3) + areadS 

           VBar = vBar + areadS*four
        end do
     end do
  end do

  ! do i=1,nx
  !    if (nPtr(1, i) .ne. 4) then
  !       vv(i*3) = vv(i*3) * five/three
  !    end if
  ! end do

  VBar = VBar / nx

  ! Reset arrays
  call VecRestoreArrayF90(xVec, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  call VecRestoreArrayF90(vVec, vv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine computeVolumes

subroutine volumeSmooth(vVec, vBar, L)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: volumeSmooth() performs two typs of volume smoothing
  !     on the list of nodal volumes, vVec. First it performs a fix
  !     number of point jacobi iterations. Then it performs a
  !     global volume blending based on vBar.
  !
  !     Description of Arguments
  !     Input:
  !     vVec - Petsc Vector: Vector of nodal volumes to smooth. New volumes are
  !            in this vector on output.
  !     vBar - real : Average nodal volume
  !     L    - integer : Marching direction
  !
  !     Ouput:
  !     vVec - See above. 

  use hypInput
  use hypData, only : nx, nPtr
  implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

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

  ! Copy current volumes to vtmp 
  do i=1,nx
     vtmp(i) = volume(3*i)
  end do

  ! Do a Jacobai volume smooth
  do iter=1,volSmoothIter
     do i=1,nx

        Volume(3*i) = zero
        oneOvrNNeighbour = one/nPtr(1,i)
        do ii=1,nPtr(1,i)
           ! Average around the neighbours
           Volume(3*i) = Volume(3*i) + oneOvrNNeighbour*& 
                ((one - volCoef) * Vtmp(i) + &
                volCoef*Vtmp(nPtr(1+ii,i)))
        end do
     end do

     ! Copy Volume back to Vtmp
     do i=1,nx
        vtmp(i) = volume(3*i)
     end do
  end do

  ! Now we do a global volume smooth
  factor = half * (one + ( one - volBlend)**(l-2))
  volume = factor*Volume + (one - factor)*Vbar

  ! Finally reset array
  call VecRestoreArrayF90(vVec, volume, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine volumeSmooth

subroutine initialGuess(X0, X1)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: initialGuess takes the current grid level X0 and
  !     determines a reasonablly good guess for the next level,
  !     X1 by performing a alebraic extrapolation. 
  !
  !     Description of Arguments
  !     Input:
  !     X0 - Petsc Vector: Vector of nodes on current grid level
  !
  !     Ouput:
  !     X1 - Petsc Vector: Initial guess of X for next grid level

  use precision
  use hypData, only : nx, deltas, nPtr, xx, xxp1
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input parameters
  Vec, intent(in) :: X0

  ! Output Parameter
  Vec, intent(out) :: X1

  ! Working Variables
  integer(kind=intType) :: i, ierr, jp1, jm1, kp1, km1, m, mm
  real(kind=realType), dimension(3) :: v1, v2, s, r0, rm
  real(kind=realType) :: thetam
  
  ! Set the pointers  
  call VecGetArrayF90(x0, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(X1, xxp1, ierr)
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
        v1 = (xx(3*jp1-2:3*jp1) - xx(3*jm1-2:3*jm1))
        v2 = (xx(3*kp1-2:3*kp1) - xx(3*km1-2:3*km1))
      
     else
        v1 = zero
        v2 = zero

        r0 = xx(3*i-2:3*i)
        MM = nPtr(1, i) 

        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = xx(3*nptr(2 + m, i)-2:3*nptr(2 + m, i))
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
     xxp1(3*i-2:3*i) = xx(3*i-2:3*i) + deltaS*s
     
  end do
  
  ! reset the pointers  
  call VecRestoreArrayF90(x0, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(X1, xxp1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine initialGuess

subroutine solveLevel(Xm1, X0, Xp1, volume)
 !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: solveLevel is the main solving routine that is
  !     called for each unknown grid level. This routine is
  !     responsible for solving the non-linear system of
  !     equations at each level. 
  !
  !     Description of Arguments
  !     Input:
  !     Xm1 - Petsc Vector: Vector of nodes on last grid level
  !     X0  - Petsc Vector: Vector of nodes on current grid level
  !     volum - Petsc Vector: Vector of nodal volumes on current grid level
  !     Ouput:
  !     X1 - Petsc Vector: Vector of nodes on next grid level. This
  !          is what we are solving for. 

  use hypInput
  use hypData, only : nx, xxm1, xx, xxp1, hypSnes, hypRHS, kspIts

  implicit none
#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input Parameters
  Vec, intent(in) :: Xm1, X0, volume

  ! Output Parameters
  Vec :: Xp1

  ! Working parameters
  integer(kind=intType) :: i, ierr, idim, l, reason
 
  ! First set the pointers for the current and previous grid
  ! level. Pointer to the next grid level Xp1 will be set in the
  ! reidual/jacobian call
  call VecGetArrayF90(Xm1, xxm1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(X0, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSet(hypRHS, zero, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESSolve(hypSnes, hypRHS, Xp1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESGetNumberFunctionEvals(hypsnes, kspits, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESGetIterationNumber(hypSnes, i, ierr)
  print *, 'nonliner its:',i

  ! Finallly restore the pointers for the three grid levels
  call VecRestoreArrayF90(Xm1, xxm1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(X0, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine solveLevel

subroutine calcResidual(snes, xVec, rVec, ctx, ierr)

  use precision
  use hypData, only : xxm1, xx, xxp1, volume, sl, scaleDist, radius0, nx, nptr
  use hypData, only : marchIter, gridSensorMax, gridSensorMin
  use hypInput
  implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! PETSc Variables
  SNES snes
  Vec xVec, rVec
  PetscFortranAddr ctx(*)
  integer(kind=intType), intent(out) :: ierr

  ! Working variables
  integer(kind=intType) :: i, idim, j, jp1, jm1, kp1, km1, ii, m, MM, iter
  real(kind=realType), dimension(3) :: r_zeta, r_ksi, r_eta, sigma, t, tau
  real(kind=realType), dimension(3) :: rplusj, rminusj, rplusk, rminusk, r0, rm
  real(kind=realType), dimension(3) :: rplusjhat, rpluskhat, nhat, G, PinvG
  real(kind=realType), dimension(3,3) :: Q1, Q2, P, eye, PinvQ1, PinvQ2, Pinv
  real(kind=realType) :: a_ksi, a_eta, De(3), factor, alpha, beta, ovrSum1, ovrSum2
  real(kind=realType) :: d_ksi, d_eta, tmp1, tmp2, tmp, thetam, Rksi, Reta, N_ksi, N_eta
  real(kind=realType), pointer :: rr(:), vv(:)
  real(kind=realType) :: dist
  ! This routine must do everything required to evaluate the full
  ! non-linear residual

  call VecGetArrayF90(xVec, xxp1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(rVec, rr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(volume, vv, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  sl = (scaleDist/(radius0*rMin))**slExp
  !sl = sqrt((marchIter-one)/(N-one))

  ! Record the maximum and minimum grid sensors
  gridSensorMax = zero 
  gridSensorMin = huge(gridSensorMin)

  ! Nodal loop
  do i=1, nx
     
     ! First compute the zeta or marching direction derivative. This
     ! will be first order for grid level = 2, (since we don't have
     ! the Xm1 layer) but will be second order for grid level > 2
     
     if (marchIter == 2)  then
        r_zeta = xxp1(3*i-2:3*i) - xx(3*i-2:3*i)
     else
        r_zeta = (three/two)*xxp1(3*i-2:3*i) - two*xx(3*i-2:3*i) + half*xxm1(3*i-2:3*i)
     end if

     if (nptr(1, i) == 4) then

        ! Extract pointers to make code easier to read
        jp1 = nPtr(2, i)
        kp1 = nPtr(3, i)
        jm1 = nPtr(4, i)
        km1 = nPtr(5, i)

        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        r_ksi = half*(xxp1(3*jp1-2:3*jp1) - xxp1(3*jm1-2:3*jm1))
        r_eta = half*(xxp1(3*kp1-2:3*kp1) - xxp1(3*km1-2:3*km1))

   
        ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
        tmp = (dist(xxm1(3*jp1-2:3*jp1), xxm1(3*i-2:3*i)) + &
               dist(xxm1(3*jm1-2:3*jm1), xxm1(3*i-2:3*i))) / &
              (dist(xx  (3*jp1-2:3*jp1), xx  (3*i-2:3*i)) + &
               dist(xx  (3*jm1-2:3*jm1), xx  (3*i-2:3*i)))
        d_ksi = max(tmp**(two/sl), 0.1)

        tmp = (dist(xxm1(3*kp1-2:3*kp1), xxm1(3*i-2:3*i)) + &
               dist(xxm1(3*km1-2:3*km1), xxm1(3*i-2:3*i))) / &
              (dist(xx  (3*kp1-2:3*kp1), xx  (3*i-2:3*i)) + &
               dist(xx  (3*km1-2:3*km1), xx  (3*i-2:3*i)))

        d_ksi = max(tmp**(two/sl), 0.1)
        d_eta = max(tmp**(two/sl), 0.1)

        gridSensorMax = max(gridSensorMax, d_ksi, d_eta)
        gridSensorMin = min(gridSensorMin, d_ksi, d_eta)
     else

        ! We have a node with 6 neighbours.
        ! Lets do this brute force:
        r_ksi = zero
        r_eta = zero

        r0 = xxp1(3*i-2:3*i)
        MM = nPtr(1, i) 

        tmp1 = zero
        tmp2 = zero
        ovrSum1 = zero
        ovrSum2 = zero
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = xxp1(3*nptr(2 + m, i)-2:3*nptr(2 + m, i))
           r_ksi = r_ksi + (two/MM) * (rm - r0) * cos(thetam)
           r_eta = r_eta + (two/MM) * (rm - r0) * sin(thetam)

           tmp1 = tmp1 + &
                dist(xxm1(nptr(2+m, i)*3-2:nptr(2+m,i)*3), xxm1(3*i-2:3*i)) / &
                dist(xx  (nptr(2+m, i)*3-2:nptr(2+m,i)*3), xx  (3*i-2:3*i))*abs(cos(thetam))

           tmp2 = tmp2 + &
                dist(xxm1(nptr(2+m, i)*3-2:nptr(2+m,i)*3), xxm1(3*i-2:3*i)) / &
                dist(xx  (nptr(2+m, i)*3-2:nptr(2+m,i)*3), xx  (3*i-2:3*i))*abs(sin(thetam))

           ovrSum1 = ovrSum1 + abs(cos(thetam))
           ovrSum2 = ovrSum2 + abs(sin(thetam))
        end do

        d_ksi = max((tmp1/ovrSum1)**(two/sl), 0.1)
        d_eta = max((tmp2/ovrSum2)**(two/sl), 0.1)

        gridSensorMax = max(gridSensorMax, d_ksi, d_eta)
        gridSensorMin = min(gridSensorMin, d_ksi, d_eta)
      end if

     ! Compute the auxiluary variables sigma, t and tau which are
     ! the three cross products of the vectors
     sigma(1) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
     sigma(2) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
     sigma(3) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)
        
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
     g(3) = vv(3*i) + two*(&
          sigma(1)*r_zeta(1) + sigma(2)*r_zeta(2) + sigma(3)*r_zeta(3))
     
     ! And the residual contribution
     Pinvg = matmul(Pinv, g)

     ! ============================================
     !             Explicit Smoothing 
     ! ============================================

     a_ksi = one
     a_eta = one
     
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

     if (nptr(1, i) == 4) then
        De = &
             epsE*Rksi*N_ksi*(xxp1(3*jm1-2:3*jm1) - two*xxp1(3*i-2:3*i) + xxp1(3*jp1-2:3*jp1)) + &
             epsE*Reta*N_eta*(xxp1(3*km1-2:3*km1) - two*xxp1(3*i-2:3*i) + xxp1(3*kp1-2:3*kp1))
     else
        de = zero
        r0 = xxp1(3*i-2:3*i)
        MM = nPtr(1, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = xxp1(3*nptr(2 + m, i)-2:3*nptr(2 + m, i))
           De = De + &
                epsE*Rksi*N_ksi*(two/MM)*(rm - r0)*(four*cos(thetam)**2 - one) + &
                epsE*Reta*N_eta*(two/MM)*(rm - r0)*(four*sin(thetam)**2 - one)
        end do

     end if

     if (marchIter == 2) then
        De = zero
     end if
 
     ! ------------------------------------
     ! Final Residual Value
     ! ------------------------------------
      rr(3*i-2:3*i) = PinvG - r_zeta - matmul(PinvQ1, r_ksi) - matmul(PinvQ2, r_eta) + DE

  end do

  

  ! Always remember to restore arrays
  call VecRestoreArrayF90(xVec, xxp1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(rVec, rr, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(volume, vv, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  ierr = 0

  call vecNorm(rVec, NORM_2, tmp, ierr)

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

  ! Just call assembly betgin/end on Jac
  call MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ierr = 0
end subroutine calcJacobian

subroutine create3DPetscVars

  use hypInput
  use hypData

  implicit none
  
  ! Working Variables
  integer(kind=intType) :: onProc(nx*3), offProc(nx*3), ierr, i, idim
  external calcResidual, calcJacobian

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

  ! Note that volume doesn't really have 3 dof per node. But we will
  ! put volume only in the third dof slot since that is where it will
  ! be used for the RHS. Having volume the same size makes things a
  ! little easier
  call VecDuplicate(hypDelta, Volume, ierr)
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
  
  call SNESSetFromOptions(hypSNES, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Most important thing is to set the "function" that computes the residual "R"
  call SNESSetFunction(hypSNES, hypRes, calcResidual, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Next we must set the jacobian and the function to evaluate the jacobian
  call MatCreateSNESMF(hypSnes, hypMatFD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESSetJacobian(hypSNES, hypMatFD, hypMat, calcJacobian, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESGetKSP(hypSnes, hypKSP, ierr)
  call KSPGetPC(hypKSP, hypPC, ierr)
  call PCSetType(hypPC, "none", ierr)


  ! Set the sens tolerances...some of these should be options eventually
  call SNESSetTolerances(hypSNES, 1e-16, kspRelTol, 1e-12, 20, 1000, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine create3DPetscVars





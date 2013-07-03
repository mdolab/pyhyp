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
  integer(kind=intType) :: i, ierr, idim, l, reason

  ! Set the number of nodes on this proc, nx in hypData as nNodes
  nx = nNodes


  ! Create the PETSc Variables
  call create3DPetscVars

  ! Copy out Xin() into the first grid slot by setting a pointer and
  ! copying Xin

  call VecGetArrayF90(X(1), xxm1tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do i=1,nx
     do idim=1,3
        xxm1tmp((i-1)*3 + idim) = Xin(idim, i)
     end do
  end do
 
  call VecRestoreArrayF90(X(1), xxm1tmp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute the inital 'Radius' value, radius0
  call computeR(X(1), radius0)
  
  ! Now the user has specified how far out to march in terms of
  ! multiplies of radius0. This gives us a nominal distance. Since we know
  ! the length of the off-wall direction and the initial spacing, we
  ! find the gridRatio that will satifiy these two conditions.
  call calcGridRatio()
  write(*,"(a)", advance="no") '#--------------------#'
  print "(1x)"  
  write(*,"(a)", advance="no") "Grid Ratio:"
  write(*,"(f8.4,1x)",advance="no") gridRatio
  print "(1x)" 
  write(*,"(a)", advance="no") '#--------------------#'
  print "(1x)"   

  ! Write the header
  call writeHeader

  ! Set initial pseudo and real spacings
  deltaS = ps0

  ! Zero out the cumulative distance counter
  scaleDist = zero
  nSubIterPrev = 0
  desiredS = zero
  ! Determine starting CPU Time
  timeStart = mpi_wtime()

  ! This is the master marching direction loop
  marchLoop: do marchIter=2, N

     ! Use 'L' as the working variable
     L = marchIter

     ! Run the "initial guess" function. If the user is running in
     ! 'linear' mode this is all that is done. This function computes
     ! the next step, X(L). It may take multiple steps to get there
     ! which is fine. 
     call initialGuess(X(L)) 

     ! ! Solve for the next grid level if nonLinear
     ! if (nonLinear) then
     !    call SNESSolve(hypSnes, PETSC_NULL_OBJECT, X(L), ierr)
     !    call EChk(ierr, __FILE__, __LINE__)

     !    call SNESGetNumberFunctionEvals(hypsnes, kspits, ierr)
     !    call EChk(ierr, __FILE__, __LINE__)
     ! end if

     ! Compute the max radius of the current level, X(L)
     call computeMinR(X(L), Radius)

     ! Check the quality of this layer
     call computeQualityLayer

     ! Possibly write header and write iteration info
     if (mod(marchIter, 50) == 0) then
        call writeHeader
     end if

     ! Write info for this layer
     call writeIteration
  end do marchLoop

  ! Destroy the PETSc variables
  !call destroyPetscVars
     
 end subroutine run3D

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

  use precision
  use hypData, only : nPatch, patches, nx, deltaS, xx, volume, nptr, vBar, cratio
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Working Variables
  integer(kind=intType) :: i, j, ipatch, ierr
  real(kind=realType) :: ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS
  real(kind=realType) :: c1,c2,c3,c4, dist
  volume = zero
  VBar = zero
  cratio = zero

  ! Loop over the patch faces:
  do iPatch=1,nPatch
     ! Loop over faces on patch
     do j=1,patches(iPatch)%jl-1
        do i=1,patches(iPatch)%il-1

           ! Extract the coordinates of the 4 corners
           ll = xx(:, patches(iPatch)%l_index(i  , j  ))
           ul = xx(:, patches(iPatch)%l_index(i  , j+1))
           lr = xx(:, patches(iPatch)%l_index(i+1, j  ))
           ur = xx(:, patches(iPatch)%l_index(i+1, j+1))

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
           ! volume. The other two slots are unsued. 
           volume(patches(iPatch)%l_index(i, j)) = &
                volume(patches(iPatch)%l_index(i, j)) + areadS 
           volume(patches(iPatch)%l_index(i+1, j)) = &
                volume(patches(iPatch)%l_index(i+1, j)) + areadS 
           volume(patches(iPatch)%l_index(i, j+1)) = &
                volume(patches(iPatch)%l_index(i, j+1)) + areadS 
           volume(patches(iPatch)%l_index(i+1, j+1)) = &
                volume(patches(iPatch)%l_index(i+1, j+1)) + areadS 

           VBar = vBar + areadS*four

           ! Also compute edge lengths
           c1 = deltaS/dist(ll,lr)
           c2 = deltaS/dist(ul,ur)
           c3 = deltaS/dist(ll,ul)
           c4 = deltaS/dist(lr,ur)

           cRatio = max(cRatio,c1,c2,c2,c4)

        end do
     end do
  end do

  VBar = VBar / nx

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

  use hypInput
  use hypData, only : nx, nPtr, volume, marchIter, vBar
  implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Working Parameters
  real(kind=realType) :: factor, Vtmp(nx), oneOvrNNeighbour
  integer(kind=intType) :: i, ii, iter, ierr

  ! Copy current volumes to vtmp 
  do i=1,nx
     vtmp(i) = volume(i)
  end do

  ! Do a Jacobai volume smooth
  do iter=1,volSmoothIter
     do i=1,nx

        volume(i) = zero
        oneOvrNNeighbour = one/abs(nPtr(1,i))
        do ii=1,abs(nPtr(1,i))
           ! Average around the neighbours
           volume(i) = volume(i) + oneOvrNNeighbour*& 
                ((one - volCoef) * Vtmp(i) + &
                volCoef*Vtmp(abs(nPtr(1+ii,i))))
        end do
     end do

     ! Copy Volume back to Vtmp
     do i=1,nx
        vtmp(i) = volume(i)
     end do
  end do

  ! Now we do a global volume smooth
  factor = half * (one + ( one - volBlend)**(marchIter-2))
  !factor = max(one -( (marchIter -2)/dble(N))**3, 0.75)
  volume = factor*volume + (one - factor)*Vbar

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
  !     Description of Arguments
  !
  !     Ouput:
  !     xNew petsc-vector : Result of linear step.

  use precision
  use hypData, only : nx, deltas, nPtr, X, xxm1, xx, xxm2, xxtmp, marchIter, inds
  use hypData, only : cRatio, hypKSP, hypMat, deltaTmp, hypRHS, hypDelta, scaleDist
  use hypData, only : kspIts, subiter, nsubiter, nsubiterprev, gridRatio, desiredS
  use hypData, only : xxInterp
  use hypInput
  implicit none

#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Output Parameters
  Vec :: Xnew

  ! Working Variables
  integer(kind=intType) :: i, ierr, jp1, jm1, kp1, km1, m, mm, ii, idim
  real(kind=realType) :: fact, onemfact
  logical :: saveMetrics, keepGoing

  ! Desired s:
  desiredS = desiredS + s0*gridRatio**(marchIter-2)
  
  ! Get values for the very first step; otherwise we will always have
  ! old values
  if (marchIter == 2) then 
     call VecGetValues(X(marchIter-1), nx*3, inds, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     xxm1 = xx
  end if

  keepGoing = .True.
  nsubIter = 0
  do while (keepGoing) 
     nSubIter = nSubIter + 1
     saveMetrics = .False.

     ! Compute volumes pased on deltaS
     call computeVolumes

     ! Adjsut deltaS if necessary
     if (cratio > cmax) then 
        deltaS = deltas*(cmax/cratio)
        call computeVolumes
     end if

      ! And then smooth as normal
     call volumeSmooth

     ! Increment our running counter on approx distance from wall
     scaleDist = scaleDist + deltaS
     
     ! Assemble the Jacobaian, second order is true, and also gets the RHS
     call calcResidual(.True., .False., saveMetrics)
     
     call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! Now solve the system
     call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPGetIterationNumber(hypKsp, kspIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  
     call VecGetArrayF90(hypDelta, deltaTmp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Copy xx to xxm1
     xxm1 = xx

     ! We have a delta use this to get updated value of xx
     do i=1,nx
        do idim=1,3
           xx(idim,i) = xx(idim,i) + deltaTmp(3*(i-1)+idim)
        end do
     end do

     call VecRestoreArrayF90(hypDelta, deltaTmp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now we need to check if we have gone far enough for the next
     ! desired grid level. If so, we compute the factor between xxm1
     ! and xx and set. 
     if (scaleDist >= desiredS) then
        fact = (desiredS - (scaleDist - deltas))/(deltaS)
        onemfact = one-fact
        do i=1,nx
           do idim=1,3
              xxInterp(idim,i) = fact*xx(idim,i) + onemfact*xxm1(idim,i)
           end do
        end do

        ! Finally take our last solution in xx and dump into Xnew
        call VecSetValues(Xnew, nx*3, inds, xxInterp, INSERT_VALUES, ierr)
        call EChk(ierr,__FILE__,__LINE__)
        keepGoing = .False.
     end if

     ! Change deltaS
     deltaS = deltaS * pGridRatio
  end do

end subroutine initialGuess

subroutine formFunction(snes, xVec, rVec, ctx, ierr)

  use precision
  use hypData, only : xxtmp, rrtmp, xx, rr, nx, marchIter
  implicit none
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! PETSc Variables
  SNES snes
  Vec xVec, rVec
  PetscFortranAddr ctx(*)
  integer(kind=intType), intent(out) :: ierr
  
  ! Working Variables
  integer(Kind=intType) :: i, idim

  ! Set pointer to xVec and the residual
  call VecGetArrayF90(xVec, xxTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(rVec, rrTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Copy xxTmp into xx
  do i=1,nx
     do idim=1,3
        xx(idim,i) = xxTmp(3*(i-1)+idim)
     end do
  end do

  if (marchIter == 2) then
     call calcResidual(.False., .False., .False.)
  else
     call calcResidual(.False., .True., .False.)
  end if
  
  ! Copy rr into rrTmp
  do i=1,nx
     do idim=1,3
        rrtmp(3*(i-1)+idim) = rr(idim,i)
     end do
  end do

  ! Always remember to restore arrays
  call VecRestoreArrayF90(xVec, xxTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(rVec, rrTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ierr = 0

end subroutine formFunction

subroutine formJacobian(snes, xVec, Jac, B, flag, ctx, ierr)
  use hypData, only : xx, xxtmp, nx, marchIter
  use hypInput 
  implicit none
#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input (Petsc) Variables
  SNES snes
  Vec xVec
  Mat Jac, B
  PetscInt flag
  PetscFortranAddr ctx(*)
  integer(kind=intType) :: ierr

  ! Working variables
  integer(kind=intType) :: I, idim

  ! Just call assembly begin/end on the FD jacobian
  call MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set pointer to xVec
  call VecGetArrayF90(xVec, xxTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Copy xxTmp into xx
  do i=1,nx
     do idim=1,3
        xx(idim,i) = xxTmp(3*(i-1)+idim)
     end do
  end do

  ! Call the residual function and assemble the jacobian
  if (marchIter == 2) then
     call calcResidual(.True., .False., .False.)
  else
     call calcResidual(.True., .True., .False.)
  end if

  ! Always remember to restore arrays
  call VecRestoreArrayF90(xVec, xxTmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ierr = 0
end subroutine formJacobian

subroutine calcResidual(assembleJac, secondOrder, saveMetrics)

  use precision
  use hypData, only : xxm2, xxm1, xx, volume, sl, scaleDist, radius0, nx, nptr
  use hypData, only : marchIter, gridSensorMax, gridSensorMin, rr, hypMat, hypRHS
  use hypData, only : ADD_VALUES, MAT_FINAL_ASSEMBLY, INSERT_VALUES
  use hypData, only : X_ksi, X_eta, X_zeta, X_ksi_ksi, X_eta_eta, X_diss, Vhist
  use hypData, only : subiter, nsubiter
  use hypInput
  implicit none

  ! Input Variables
  logical, intent(in) :: assembleJac, secondOrder, saveMetrics

  ! Working variables
  integer(kind=intType) :: i, idim, j, jp1, jm1, kp1, km1, ii, m, MM, iter, ierr
  integer(kind=intType) :: nAverage
  real(kind=realType), dimension(3) :: r_zeta, r_ksi, r_eta, sigma, t, tau
  real(kind=realType), dimension(3) :: r_ksi_ksi, r_eta_eta
  real(kind=realType), dimension(3) :: rplusj, rminusj, rplusk, rminusk, r0, rmmr0
  real(kind=realType), dimension(3) :: rplusjhat, rpluskhat, nhat, G, PinvG
  real(kind=realType), dimension(3,3) :: Q1, Q2, P, eye, PinvQ1, PinvQ2, Pinv
  real(kind=realType) :: a_ksi, a_eta, De(3), factor, alpha, beta, ovrSum1, ovrSum2
  real(kind=realType) :: d_ksi, d_eta, tmp1, tmp2, tmp, thetam, Rksi, Reta, N_ksi, N_eta
  real(kind=realType) :: dist, coef,  deltaVovrDetC, detC, largeWeight, smallWeight
  ! This routine must do everything required to evaluate the full
  ! non-linear residual

  sl = (scaleDist/(radius0*rMin))**slExp
  !sl = sqrt((marchIter-one)/(N-one))
  !sl = sqrt((marchIter + (dble(subIter)-1)/nSubIter)/(73-one))
  !sl = (scaleDist/(radius0*rMin))**.16
  !print *,'sl:',sl
  ! Record the maximum and minimum grid sensors
  gridSensorMax = zero 
  gridSensorMin = huge(gridSensorMin)
  if (assembleJac) then
     eye = zero
     eye(1,1) = one
     eye(2,2) = one
     eye(3,3) = one
     call MatZeroEntries(hypMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecZeroEntries(hypRHS, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Nodal loop
  do i=1, nx
     
     ! First compute the zeta or marching direction derivative. This
     ! will be first order for grid level = 2, (since we don't have
     ! the Xm1 layer) but will be second order for grid level > 2
     
     if (nonLinear) then
        if (.not. secondOrder) then
           r_zeta = xx(:,i) - xxm1(:,i)
        else
           r_zeta = (three/two)*xx(:,i) - two*xxm1(:,i) + half*xxm2(:, i)
        end if
     end if

     if (abs(nptr(1, i)) == 4) then

        ! Extract pointers to make code easier to read
        jp1 = abs(nPtr(2, i))
        kp1 = abs(nPtr(3, i))
        jm1 = abs(nPtr(4, i))
        km1 = abs(nPtr(5, i))

        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        r_ksi = half*(xx(:,jp1) - xx(:, jm1))
        r_eta = half*(xx(:,kp1) - xx(:, km1))
        r_ksi_ksi = half*(xx(:, jp1) - two*xx(:, i) + xx(:, jm1))
        r_eta_eta = half*(xx(:, kp1) - two*xx(:, i) + xx(:, km1))

        ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
        tmp = (dist(xxm1(:,jp1), xxm1(:,i)) + &
               dist(xxm1(:,jm1), xxm1(:,i))) / &
              (dist(xx  (:,jp1), xx  (:,i)) + &
              dist(xx  (:,jm1), xx  (:,i)))
        d_ksi = max(tmp**(two/sl), 0.1)

        tmp = (dist(xxm1(:,kp1), xxm1(:,i)) + &
               dist(xxm1(:,km1), xxm1(:,i))) / &
              (dist(xx  (:,kp1), xx  (:,i)) + &
               dist(xx  (:,km1), xx  (:,i)))
        d_eta = max(tmp**(two/sl), 0.1)
     else !if (nPtr(1, i) > 0) then

        ! We have a node with 6 neighbours.
        ! Lets do this brute force:
        r_ksi = zero
        r_eta = zero
        r_ksi_ksi = zero
        r_eta_eta = zero

        r0 = xx(:,i)
        MM = nPtr(1, i) 

        tmp1 = zero
        tmp2 = zero
        ovrSum1 = zero
        ovrSum2 = zero
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rmmr0 = xx(:,nptr(2 + m, i)) - r0
           r_ksi = r_ksi + (two/MM) * (rmmr0) * cos(thetam)
           r_eta = r_eta + (two/MM) * (rmmr0) * sin(thetam)
           r_ksi_ksi = r_ksi_ksi + (two/MM)*(rmmr0)*(four*cos(thetam)**2 - one)
           r_eta_eta = r_eta_eta + (two/MM)*(rmmr0)*(four*sin(thetam)**2 - one)

           tmp1 = tmp1 + &
                dist(xxm1(:,nptr(2+m, i)), xxm1(:,i)) / &
                dist(xx  (:,nptr(2+m, i)), xx  (:,i))*abs(cos(thetam))

           tmp2 = tmp2 + &
                dist(xxm1(:,nptr(2+m, i)), xxm1(:,i)) / &
                dist(xx  (:,nptr(2+m, i)), xx  (:,i))*abs(sin(thetam))

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
     
        deltaVovrDetC = Volume(i)/detC
     
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
     g(3) = volume(i) + two*(&
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

     De = epsE*(Rksi*N_ksi*r_ksi_ksi + Reta*N_eta*r_eta_eta)

     if (marchIter == 2) then
        De = zero
     end if
     if ( nPtr(1,i)< 0) then
         de = zero
      end if
     ! ------------------------------------
     ! Final Residual Value
     ! ------------------------------------
     if (.not. assembleJac) then
        rr(:,i) = PinvG - r_zeta - matmul(PinvQ1, r_ksi) - matmul(PinvQ2, r_eta) + DE
     else
        ! Assemble the LHS, once again we need to split between regular and extraordinary nodes.

        if (abs(nptr(1, i)) == 4) then
           ! jp1, jm1 etc are still set from above
           
           ! -----------------
           ! Left Hand Side:
           ! -----------------

           ! Point to the left
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm1-1, &
                (-(one+theta)*half*PinvQ1 - epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Center Point ! 
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                (eye + four*epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           
           ! Point to the right
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, jp1-1, &
                ((one+theta)*half*PinvQ1 - epsI*eye),  ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! ----------------------------------------------------------
           
           ! Point to the bottom
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, km1-1, &
                (-(one+theta)*half*PinvQ2 - epsI*eye), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           
           ! Point to the right
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, kp1-1, &
                ((one+theta)*half*PinvQ2 - epsI*eye),  ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, i-1, &
                matmul(Pinv,(/zero, zero, Volume(i)/)) + De, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

        else if (nPtr(1, i)  > 0) then
           MM = nPtr(1, i) 
           do m = 0, MM - 1
              thetam = 2*pi*m/MM
              ! For the 'fm' point in ksi
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                   (one + theta)*PinvQ1*(two/MM)*cos(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in ksi
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                   -(one + theta)*PInvQ1*(two/MM)*cos(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
              
              ! For the 'fm' point in eta
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                   (one + theta)*PinvQ2*(two/MM)*sin(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
              
              ! For the 'f0' point in ksi
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                   -(one + theta)*PInvQ2*(two/MM)*sin(thetam), ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! Now we have the second derivative values to do in ksi-ksi
              
              ! For the 'fm' point in ksi-ksi 
              coef = (two/MM)*(four*cos(thetam)**2 - one)
              
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                   -eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              ! For the 'f0' point in ksi-ksi
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                   eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
              
              ! For the 'fm' point in eta-eta
              coef = (two/MM)*(four*sin(thetam)**2 - one)
              
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                   -eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
              
              ! For the 'f0' point in eta-eta
              call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                   eye*epsI*coef, ADD_VALUES, ierr)
              call EChk(ierr, __FILE__, __LINE__)
              
           end do

           ! Finally we need an eye on the diagonal
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                eye, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__) 

           ! We also need to set the RHS
           call VecSetValuesBlocked(hypRHS, 1, i-1, &
                matmul(Pinv,(/zero, zero, Volume(i)/)) + De, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           
        else
           print *,'should not be here'
           stop
           ! We have a negative number. That means we have to do the
           ! super funny averaging.
           
           !We have to do averaging at the corners
        
           !  Get the number of nodes to average: REMEMBER THIS IS CURRENTLY NEGATIVE!
           if (nPtr(1,i) < 0) then
              nAverage = -nPtr(1, i) 
              !largeWeight = quarter
              !smallWeight = quarter
              largeWeight = one/nAverage
              smallWeight = one/nAverage

           else
              nAverage = nPtr(1, i)
              largeWeight = one/nAverage
              smallWeight = one/nAverage
           end if
           
           ! Loop over neighbours
           do ii = 1, nAverage
              if (nPtr(i+ii, i) < 0) then
                 call MatSetValuesBlocked(hypMat, 1, i-1, 1, -nPtr(1+ii, i)-1, &
                      -largeWeight*eye, ADD_VALUES, ierr)
              else
                 call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(1+ii, i)-1, &
                      -smallWeight*eye, ADD_VALUES, ierr)
              end if
           end do

           ! Center Point
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                eye, ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if ! Type of node select
     end if

     if (saveMetrics .and. writeMetrics) then
        call VecSetValuesBlocked(X_ksi(marchIter-1), 1, i-1, r_ksi, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(X_eta(marchIter-1), 1, i-1, r_eta, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(X_zeta(marchIter-1), 1, i-1, r_zeta, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(X_ksi_ksi(marchIter-1), 1, i-1, r_ksi_ksi, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(X_eta_eta(marchIter-1), 1, i-1, r_eta_eta, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(X_diss(marchIter-1), 1, i-1, De, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecSetValuesBlocked(Vhist(marchIter-1), 1, i-1, (/zero,zero, volume(i)/), INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

     end if

  end do ! NODE LOOP

  if (assembleJac) then
     call MatAssemblyBegin(hypMat, MAT_FINAL_ASSEMBLY, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  
     call MatAssemblyEnd(hypMat, MAT_FINAL_ASSEMBLY, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecAssemblyBegin(hypRHS, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecAssemblyEnd(hypRHS, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if
  
  if (saveMetrics .and. writeMetrics) then
     call vecAssemblyBegin(X_ksi(marchIter-1), ierr)
     call vecAssemblyEnd(X_ksi(marchIter-1), ierr)

     call vecAssemblyBegin(X_eta(marchIter-1), ierr)
     call vecAssemblyEnd(X_eta(marchIter-1), ierr)

     call vecAssemblyBegin(X_zeta(marchIter-1), ierr)
     call vecAssemblyEnd(X_zeta(marchIter-1), ierr)

     call vecAssemblyBegin(X_ksi_ksi(marchIter-1), ierr)
     call vecAssemblyEnd(X_ksi_ksi(marchIter-1), ierr)

     call vecAssemblyBegin(X_eta_eta(marchIter-1), ierr)
     call vecAssemblyEnd(X_eta_eta(marchIter-1), ierr)

     call vecAssemblyBegin(X_diss(marchIter-1), ierr)
     call vecAssemblyEnd(X_diss(marchIter-1), ierr)
  end if

end subroutine calcResidual

subroutine create3DPetscVars
 use hypInput
  use hypData

  implicit none
  
  ! Working Variables
  integer(kind=intType) :: ierr, i, idim, bs
  integer(kind=intType), dimension(:), allocatable :: onProc, offProc
  external formFunction, formJacobian

  ! ----------------------------------------------------------
  !          Linearized Hyperbolic System Variables
  ! ----------------------------------------------------------

  ! Lets to things in the proper way, create the Mat first
  allocate(onProc(nx), offProc(nx), stat=ierr)
  call EChk(ierr, __FILE__, __LINE__)
  do i=1,nx
     onProc(i) = 35
     offProc(i) = 1
  end do

  ! Create a blocked matrix
  bs = 3
#if PETSC_VERSION_MINOR < 3
  call MatCreateMPIBAIJ(PETSC_COMM_WORLD, bs, &
       nx*bs, nx*bs, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, hypMat, ierr)
#else
  call MatCreateBAIJ(PETSC_COMM_WORLD, bs, &
       nx*bs, nx*bs, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, hypMat, ierr)
#endif
  call EChk(ierr, __FILE__, __LINE__)
  deallocate(onProc, offProc, stat=ierr)
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

  ! Create the full list of grid vectors:
  allocate(X(N), stat=ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (writeMetrics) then
     allocate(X_ksi(N), X_eta(N), X_zeta(N), &
          X_ksi_ksi(N), X_eta_eta(N), X_diss(N), Vhist(N), stat=ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Loop over each "vec" in the X arrays:
  do i=1,N
     call VecDuplicate(hypDelta, X(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)

     if (writeMetrics) then
        call VecDuplicate(hypDelta, X_ksi(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDuplicate(hypDelta, X_eta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDuplicate(hypDelta, X_zeta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDuplicate(hypDelta, X_ksi_ksi(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDuplicate(hypDelta, X_eta_eta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDuplicate(hypDelta, X_diss(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDuplicate(hypDelta, Vhist(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if
  end do

  ! This should come as a option, but we'll put it here for now.
  call PetscOptionsSetValue("-pc_factor_levels","1", ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the KSP object
  call KSPCreate(petsc_comm_world, hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetFromOptions(hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGMRESSetRestart(hypKSP, kspSubspacesize, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetTolerances(hypKSP, kspRelTol, 1e-16, 1e5, kspMaxIts, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the SNES object
  call SNESCreate(PETSC_COMM_WORLD, hypSNES, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call SNESSetFromOptions(hypSNES, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Most important thing is to set the "function" that computes the residual "R"
  call SNESSetFunction(hypSNES, hypRes, formFunction, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Next we must set the jacobian and the function to evaluate the jacobian
  call MatCreateSNESMF(hypSnes, hypMatFD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call SNESSetJacobian(hypSNES, hypMatFD, hypMat, formJacobian, ctx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the sens tolerances...some of these should be options eventually
  call SNESSetTolerances(hypSNES, 1e-16, 1e-8, 1e-12, 50, 5000, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Finally allocate local fortran data arrays
  allocate(xx(3,nx),rr(3,nx), xxm1(3,nx), xxm2(3,nx), inds(3*nx), volume(nx))
  allocate(xxInterp(3,nx))
  do i=1,3*nx
     inds(i) = i-1
  end do

end subroutine create3DPetscVars

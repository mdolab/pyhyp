subroutine run2D(Xin, nNodes)
  
  use hypInput
  use hypData

  implicit none
#include "include/petscversion.h"

  ! Input Parameters
  real(kind=realType) :: Xin(2, nNodes)
  integer(kind=intType) :: nNodes

  ! Working parameters
  integer(kind=intType) :: i, j, l, idim
  real(kind=realType) :: Area1(nNodes), ABar, sMax
  real(kind=realType) :: xMin, xMax, yMin, yMax, fact, onemfact
  logical keepGoing
  ! Set the number of nodes on this proc, nx in hypData as nNodes
  nx = nNodes
  
  ! Deallocate any existing data
  if (two_d_vars_allocated) then
     call releaseMemory
  end if

  allocate(grid2D(2, nx, N))
  allocate(X0(2, nx), X1(2, nx), Xm1(2, nx))
  grid2D = zero
  ! Copy Xin into the first slot of grid2D
  do i=1,nx
     do idim=1,2
        grid2D(idim, i, 1) = Xin(idim, i)
     end do
  end do

  ! Create the petsc variables
  call create2DPetscVars

  ! Compute maximum dimension in x or y and take this is radius0
  xMin = minval(Xin(1, :))
  xMax = maxval(Xin(1, :))
  yMin = minval(Xin(2, :))
  yMax = maxval(Xin(2, :))

  radius0 = max((xMax-xMin), (yMax-yMin))
 
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
 
  ! Compute an approximate estimate of the distance to farfield by
  ! taking s0 and computing
  sMax = s0*gridRatio**(N-1)

  ! Set X0 and Xm1 to the first grid level
  X0  = grid2D(:, :, 1)
  Xm1 = X0

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
     L = marchIter
     desiredS = desiredS + s0*gridRatio**(marchIter-2)
     keepGoing = .True.
     nsubIter = 0
     do while (keepGoing) 
        nSubIter = nSubIter + 1

        ! Compute the nodal areas on the X0 layer as well as the average
        ! area ABar

        ! Compute areas based on deltaS
        call computeAreas(X0, Area1, Abar)
        
        ! Adjsut deltaS if necessary
        if (cratio > cmax) then 
           deltaS = deltas*(cmax/cratio)
           call computeAreas(X0, Area1, Abar)
        end if

        ! Now perform the "Volume" smoothing operation.
        call AreaSmooth(Area1, ABar, l)
        
        ! Increment our running counter on approx distance from wall
        scaleDist = scaleDist + deltaS
     
        call assembleAndSolve(X0, X1, Xm1, Area1)

        ! Now we need to check if we have gone far enough for the next
        ! desired grid level. If so, we compute the factor between xxm1
        ! and xx and set. 
        if (scaleDist >= desiredS) then
           fact = (desiredS - (scaleDist - deltas))/(deltaS)
           onemfact = one-fact

           do i=1,nx
              do idim=1,2
                 grid2d(idim, i, L) = fact*X1(idim, i) + onemfact*X0(idim, i)
              end do
           end do
           keepGoing = .False.
        end if
        
        ! Shuffle the grid and areas backwards
        Xm1 = X0
        X0 = X1

        ! Change deltaS
        deltaS = deltaS * pGridRatio
     end do

     ! Possibly write header and write iteration info
     if (mod(marchIter, 50) == 0) then
        call writeHeader
     end if
     
     ! Write info for this layer
     call writeIteration
  end do marchLoop

end subroutine run2D

subroutine computeAreas(X, nodalArea, ABar)
  use precision
  use hypData, only: nx, deltas, cRatio
  implicit none

  ! Compute the nodal areas of a 2D periodic curve X of length xn,
  ! using a scattering approach. Also compute the average area ABar
  ! at the same time for efficiency

  ! Input Parameters
  real(kind=realType), intent(in) :: X(2, nx)
  
  ! Output Parameters
  real(kind=realType), intent(out) :: nodalArea(nx), ABar
  
  ! Working Variables
  integer(kind=intType) :: i
  real(kind=realType) :: halfA, len, c1

  ! Zero nodalArea and ABar
  nodalArea = zero
  ABar = zero
  cRatio = zero
  ! Loop over the interior segments
  do i=1,nx-1
     len = sqrt( (X(1,i+1)-X(1,i))**2 + (X(2,i+1)-X(2,i))**2 )
     halfA = deltaS*half*len
     nodalArea(i)   = nodalArea(i)   + halfA
     nodalArea(i+1) = nodalArea(i+1) + halfA
     ABar = ABar + halfA
     c1 = deltaS/len
     cRatio = max(cRatio, c1)
  end do
  
  ! Special case for the last segment
  halfA = deltaS*half*sqrt( (X(1,nx)-X(1,1))**2 + (X(2,nx)-X(2,1))**2 )
  nodalArea(nx) = nodalArea(nx) + halfA
  nodalArea(1)  = nodalArea(1)  + halfA
  ABar = ABar + halfA

  ! Finally multiply ABar by two to account for the fact we've been
  ! adding half areas and divide by nx
  ABar = (ABar * 2) / nx

end subroutine computeAreas
  
subroutine areaSmooth(Area, ABar, l)

  use hypInput
  use hypData

  implicit none

  ! Perform averaging on the nodal areas
  
  ! Input/Output Parameters
  real(kind=realType), intent(inout) :: Area(nx)
  real(kind=realType), intent(in) :: ABar
  integer(kind=intType), intent(in) :: l

  ! Working Parameters
  real(kind=realType) :: factor, Atmp(nx)
  integer(kind=intType) :: i, idp1, idm1, ierr, iter
  
  ! Do a local point jacobi smoothing iterations
  Atmp = Area
  do iter=1,volSmoothIter
     do i=1,nx
        if (i==1) then
           idp1 = 2
           idm1 = nx
        else if(i==nx) then
           idp1 = 1
           idm1 = nx-1
        else
           idp1 = i + 1
           idm1 = i - 1
        end if
         
        Area(i) = (one - volCoef) * Atmp(i) + &
             half*volCoef*(Atmp(idp1) + Atmp(idm1))
     end do
     
     ! Copy Area back to Atmp
     Atmp = Area
  end do
  
  ! Next we will do a global area smoothing. This uses Abar which is
  ! the avarage area: (Vector notation below)

  factor = half * (one + ( one - volBlend)**(l-2))
  Area = factor*Area + (one - factor)*Abar

end subroutine areaSmooth

subroutine assembleAndSolve(X0, X1, Xm1, Area1)
  
  use hypInput
  use hypData, only : nx, hypRHS, hypKSP, hypDelta, ADD_VALUES, INSERT_VALUES
  use hypData, only : SAME_NONZERO_PATTERN, hypMat, MAT_FINAL_ASSEMBLY
  use hypData, only : gridSensorMin, gridSensorMax, marchIter, kspIts, sl
  use hypData, only : scaleDist, radius0
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(2, nx), Area1(nx), Xm1(2, nx)
  real(kind=realType), intent(out) :: X1(2, nx)

  ! Working parameters
  real(kind=realType) :: A0(2,2), B0(2,2), r_eta(2), r_ksi(2), r_ksi_ksi(2)
  integer(kind=intType) :: jm1, j, jp1, ierr
  real(kind=realType) :: eye(2, 2), BInv(2,2), BInvA(2,2), deltaR(2)
  real(kind=realType) :: De(2), a_ksi, tmp, d_ksi, N_ksi
  real(kind=realType) :: exp, gamma, ovrGamma, T, D, lambda1, lambda2
  real(kind=realType) :: nVec(2), nHat(2), v1(2), rplusj(2), rminusj(2)

  ! Define a 2x2 'I' matrix
  eye = zero
  eye(1,1) = one
  eye(2,2) = one
  call MatZeroEntries(hypMat, ierr)
  call VecZeroEntries(hypRHS, ierr)

  ! ------------------------------------
  ! Scaling Function (Equation 6.5)
  ! ------------------------------------

  sl = (scaleDist/(radius0*rMin))**slExp
  gridSensorMax = zero 
  gridSensorMin = huge(gridSensorMin)

  ! Next we assemble hypMat which is quite straight forward
  do j=1, nx

     ! Get the indices (fortran ordering) depending on whcih point we
     ! have
     if (j==1) then
        jp1 = 2
        jm1 = nx
     else if(j==nx) then
        jp1 = 1
        jm1 = nx-1
     else
        jp1 = j + 1
        jm1 = j - 1
     end if
     
     ! Compute centered difference with respect to zeta. This is the
     ! nominal calculation that we will use in the interior of the domain
     r_ksi = half*(X0(:, jp1)- X0(:, jm1))
     r_ksi_ksi = half*(X0(:,jp1) - two*X0(:,j) + X0(:,jm1))
  
     ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
     
     tmp = (dist(Xm1(:,jp1), Xm1(:,j) + &
          dist(Xm1(:,jm1), Xm1(:,j))) / &
          (dist(X0(:,jp1), X0(:,j)) + &
          dist(X0 (:,jm1), X0(:,j))))
     d_ksi = max(tmp**(two/sl), 0.1)

     gridSensorMax = max(gridSensorMax, d_ksi)
     gridSensorMin = min(gridSensorMin, d_ksi)

     gamma = r_ksi(1)**2 + r_ksi(2)**2
     ovrGamma = one/gamma

     r_eta(1) = -r_ksi(2)*Area1(j)*ovrGamma
     r_eta(2) =  r_ksi(1)*Area1(j)*ovrGamma
     
     ! Compute the two coefficient matrices
     A0(1, 1) = r_eta(1)
     A0(1, 2) = r_eta(2)
     A0(2, 1) = r_eta(2)
     A0(2, 2) = -r_eta(1)

     B0(1, 1) = r_ksi(1)
     B0(1, 2) = r_ksi(2)
     B0(2, 1) = -r_ksi(2)
     B0(2, 2) = r_ksi(1)

     ! Compute the inverse of B and its multiplcation with A
     call two_by_two_inverse(B0, Binv)
     BinvA = matmul(Binv, A0)
  
     ! -----------------
     ! Left Hand Side:
     ! -----------------

     ! Point to the left
     call MatSetValuesBlocked(hypMat, 1, j-1, 1, jm1-1,&
          -(one+theta)*half*BinvA - epsI*eye, INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Center Point
     call MatSetValuesBlocked(hypMat, 1, j-1, 1, j-1, &
          eye + two*epsI*eye, INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Point to the right
     call MatSetValuesBlocked(hypMat, 1, j-1, 1, jp1-1, &
          (one+theta)*half*BinvA - epsI*eye,  INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! -----------------
     ! Right Hand Side:
     ! -----------------

     ! Assemble the the Binv comonent of the RHS
     call VecSetValuesBlocked(hypRHS, 1, j-1, &
          matmul(Binv,(/zero, Area1(j)/)), INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! ============================================
     !             Explicit Smoothing 
     ! ============================================

     ! ------------------------------------
     ! Compute the grid angle functions (Equation 6.9)
     ! ------------------------------------
     a_ksi = one
     
     ! ------------------------------------
     ! Combine into R_ksi (Equation 6.4)
     ! ------------------------------------
     R_ksi = Sl * d_ksi * a_ksi

     ! ------------------------------------
     ! Matrix Norm approximation 
     ! ------------------------------------
     ! Matrix norm approximation
     ! Compute the eigenvalues of Binv
     N_ksi = sqrt(BinvA(1,1)**2 + Binva(1,2)**2)
  
     ! ------------------------------------
     ! Final explict dissipation (Equation 6.1 and 6.2)
     ! ------------------------------------
     De = epsE*R_ksi*N_ksi *r_ksi_ksi

     if (marchIter == 2) then
        De = zero
     end if

     call VecSetValuesBlocked(hypRHS, 1, j-1, De, ADD_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

  ! Assemble the matrix and vector
  call MatAssemblyBegin(hypMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(hypMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyBegin(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the operator so PETSc knows to refactor
  !if (mod(marchIter,10) == 0) then
#if PETSC_VERSION_MINOR > 4
  call KSPSetOperators(hypKSP, hypMat, hypMat, ierr)
#else
  call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
#endif
  call EChk(ierr, __FILE__, __LINE__)
  !end if
  ! Now solve the system
  call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGetIterationNumber(hypKsp, kspIts, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Take the solution in deltaR and add to X0 to get X1
  do j=1,nx
     call VecGetValues(hypDelta, 2, (/2*(j-1), 2*(j-1)+1/), deltaR, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     X1(:, j) = X0(:, j) + deltaR
  end do
 
  contains 

    function dist(p1, p2)

      real(kind=realType) :: p1(2), p2(2)
      real(kind=realType) :: dist
      dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2)

    end function dist

end subroutine assembleAndSolve

subroutine two_by_two_inverse(A, Ainv)

  use precision

  implicit none

  ! Input/Ouput
  real(kind=realType), intent(in) :: A(2,2)
  real(kind=realType), intent(out) :: Ainv(2,2)

  ! Working
  real(kind=realType) :: idet
  
  idet = one/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  Ainv(1,1) = idet*A(2,2)
  Ainv(1,2) = -idet*A(1,2)
  Ainv(2,1) = -idet*A(2,1)
  Ainv(2,2) = idet*A(1,1)

end subroutine two_by_two_inverse

subroutine create2DPetscVars

  use hypInput
  use hypData
  
  implicit none
    
  ! Working Variables
  integer(kind=intType) :: onProc(nx*2), offProc(nx*2), ierr, bs

  ! ----------------------------------------------------------
  !          Linearized Hyperbolic System Variables
  ! ----------------------------------------------------------

  ! Lets to things in the proper way, create the Mat first
  onProc = 3
  offProc = 0

  bs = 2
  call MatCreateBAIJ(PETSC_COMM_WORLD, bs, &
       nx*bs, nx*bs, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! This must be set to allow passing in blocks in native fortran
  ! ordering
  call MatSetOption(hypMat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Then use getVecs to get the vectors we want
  call MatGetVecs(hypMat, hypDelta, hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the KSP Object
  call KSPCreate(petsc_comm_world, hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetFromOptions(hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#if PETSC_VERSION_MINOR > 4
  call KSPSetOperators(hypKSP, hypMat, hypMat, ierr)
#else
  call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
#endif
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetTolerances(hypKSP, 1e-12, 1e-16, 1e5, 50, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  two_d_vars_allocated = .True.

end subroutine create2DPetscVars


 ! ! We are going to do this smoothing properly by solving the linear
 !  ! system. The equation we will be solving is:
 !  !
 !  ! AA = (1 - volCoef)*Area(i) + half*volCoef*(Area(i-1) + Area(i+1))
 !  !
 !  ! Where AA is the smoothed area. 

 !  if (.not. smoothingAssembled) then
 !     call MatZeroEntries(smoothMat, ierr)
 !     call VecZeroEntries(smoothRHS, ierr)

 !     ! Next we assemble smoothMat which is quite straight forward
 !     do i=1, nx
        
 !        ! Get the indices (fortran ordering) depending on whcih point we
 !        ! have
 !        if (i==1) then
 !           idp1 = 2
 !           idm1 = nx
 !        else if(i==nx) then
 !           idp1 = 1
 !           idm1 = nx-1
 !        else
 !           idp1 = i + 1
 !           idm1 = i - 1
 !        end if

 !        ! -----------------
 !        ! Left Hand Side:
 !        ! -----------------

 !        ! Point to the left
 !        call MatSetValues(smoothMat, 1, i-1, 1, idm1-1, one, INSERT_VALUES, ierr)
 !        call EChk(ierr, __FILE__, __LINE__)

 !        ! Center Point
 !        call MatSetValues(smoothMat, 1, i-1, 1, i-1, -two, INSERT_VALUES, ierr)
 !        call EChk(ierr, __FILE__, __LINE__)

 !        ! Point to the right
 !        call MatSetValues(smoothMat, 1, i-1, 1, idp1-1, one, INSERT_VALUES, ierr)
 !        call EChk(ierr, __FILE__, __LINE__)
 !     end do

 !     ! Assemble the matrix and vector
 !     call MatAssemblyBegin(smoothMat, MAT_FINAL_ASSEMBLY, ierr)
 !     call EChk(ierr, __FILE__, __LINE__)

 !     call MatAssemblyEnd(smoothMat, MAT_FINAL_ASSEMBLY, ierr)
 !     call EChk(ierr, __FILE__, __LINE__)

 !     call MatView(smoothMat, PETSC_VIEWER_STDOUT_WORLD, ierr)
 !     ! Set the operators
 !     call KSPSetOperators(smoothKSP, smoothMat, smoothMat, SAME_NONZERO_PATTERN, ierr)
 !     call EChk(ierr, __FILE__, __LINE__)
     
 !     ! Flag as assembled
 !     smoothingAssembled = .True.
 !  end if

 !  ! Now assemle the RHS:
 !  do i=1,nx
 !     print *,'before:',i, Area(i)
 !     call vecSetValues(smoothRHS, 1, i-1, Area(i), INSERT_VALUES, ierr)
 !     call EChk(ierr, __FILE__, __LINE__)
 !  end do

 !  call VecAssemblyBegin(smoothRHS, ierr)
 !  call EChk(ierr, __FILE__, __LINE__)
  
 !  call VecAssemblyEnd(smoothRHS, ierr)
 !  call EChk(ierr, __FILE__, __LINE__)

 !  ! Now solve the system
 !  call KSPSolve(smoothKSP, smoothRHS, smoothDelta, ierr)
 !  call EChk(ierr, __FILE__, __LINE__)
  
 !  ! The solution is in smoothDelta which isn't actually a delta, but
 !  ! the smoothed areas:

 !  ! Extract the solution
 !  do i=1,nx
 !     call vecGetValues(smoothDelta, 1, i-1, Area(i), ierr)
 !     call EChk(ierr, __FILE__, __LINE__)
 !     print *,'after:',i,Area(i)
 !  end do

subroutine run2D(Xin, nx)
  
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  real(kind=realType) :: Xin(2, nx)
  integer(kind=intType) :: nx

  ! Working parameters
  integer(kind=intType) :: i, j, l, idim
  real(kind=realType), pointer, dimension(:, :) :: X0, X1
  real(kind=realType) :: Area0(nx), Area1(nx), ABar, deltaS
  real(kind=realType) :: A0(2, 2, nx), B0(2, 2, nx) 

  ! First thing we will do is allocate the final grid array:
  if (allocated(grid2D)) then
     deallocate(grid2D)
  end if

  allocate(grid2D(2, nx, N))
  grid2D = zero
  ! Copy Xin into the first slot of grid2D
  do i=1,nx
     do idim=1,2
        grid2D(idim, i, 1) = Xin(idim, i)
     end do
  end do

  ! Compute the original Areas, based on initial s0 value
  X0 => grid2D(:, :, 1)
  call computeAreas(X0, Area0, nx, ABar, s0)

  ! Create the petsc variables
  call create2DPetscVars(nx)

  ! This is the master marching direction loop
  marchLoop: do l=2, N
     print *,'----------- Step: ', l

     ! Set pointers to the two levels we are currently working with
     X0 => grid2D(:, :, l-1)
     X1 => grid2D(:, :, l  )

     ! Compute the length increment in the marching direction. This is
     ! put in a separate function to allow for potential callbacks to
     ! user supplied functions if necessary
     
     call computeStretch(l, deltaS)

     ! Compute the nodal areas on the X0 layer as well as the average
     ! area ABar
     call computeAreas(X0, Area1, nx, Abar, deltaS)

     ! Now perform the "Volume" smoothing operation.
     call AreaSmooth(Area1, nx, ABar, l)

     ! Compute the metrics
     call computeMetrics2D(X0, Area0, Area1, A0, B0, nx)

     ! Assemble and solve
     call assembleAndSolve(X0, X1, Area1, A0, B0, nx)

     ! Shuffle the area's backwards
     Area0 = Area1

  end do marchLoop

  ! Destroy the PETSc variables
  call destroy2DPetscVars

end subroutine run2D

subroutine computeAreas(X, nodalArea, nx, ABar, deltaS)

  use precision

  implicit none

  ! Compute the nodal areas of a 2D periodic curve X of length xn,
  ! using a scattering approach. Also compute the average area ABar
  ! at the same time for efficiency

  ! Input Parameters
  real(kind=realType), intent(in) :: X(2, nx), deltaS
  integer(kind=intType), intent(in) :: nx
  
  ! Output Parameters
  real(kind=realType), intent(out) :: nodalArea(nx), ABar
  
  ! Working Variables
  integer(kind=intType) :: i
  real(kind=realType) :: halfA

  ! Zero nodalArea and ABar
  nodalArea = zero
  ABar = zero

  ! Loop over the interior segments
  do i=1,nx-1
     halfA = deltaS*half*sqrt( (X(1,i+1)-X(1,i))**2 + (X(2,i+1)-X(2,i))**2 )
     nodalArea(i)   = nodalArea(i)   + halfA
     nodalArea(i+1) = nodalArea(i+1) + halfA
     ABar = ABar + halfA
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
  
subroutine areaSmooth(Area0, nx, ABar, l)

  use hypInput

  implicit none

  ! Perform averaging on the nodal areas
  
  ! Input/Output Parameters
  real(kind=realType), intent(inout) :: Area0(nx)
  real(kind=realType), intent(in) :: ABar
  integer(kind=intType), intent(in) :: nx
  integer(kind=intType), intent(in) :: l

  ! Working Parameters
  real(kind=realType) :: factor, Atmp(nx)
  integer(kind=intType) :: i, iter, idp1, idm1

  ! Do point jacobi smoothing iterations
  Atmp = Area0
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
        
        Area0(i) = (one - volCoef) * Atmp(i) + &
             half*volCoef*(Atmp(idp1) + Atmp(idm1))
     end do
     
     ! Copy Area0 back to Atmp
     Atmp = Area0
  end do
end subroutine areaSmooth

subroutine computeMetrics2D(X0, Area0, Area1, A0, B0, nx)

  use precision

  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(2, nx), Area0(nx), Area1(nx)
  real(kind=realType), intent(out) :: A0(2, 2, nx), B0(2, 2, nx)
  integer(kind=intType), intent(in) :: nx

  ! Working Variables
  real(kind=realType) :: alpha, beta, gamma, lambda, ovrGamma, rip1(2), rim1(2)
  real(kind=realType) :: rzeta(2), reta(2)
  integer(kind=intType) :: i

  do i=1, nx
     if (i==1) then
        rip1 = X0(:, 2)
        rim1 = X0(:, nx)
     else if(i==nx) then
        rip1 = X0(:, 1)
        rim1 = X0(:, nx-1)
     else
        rip1 = X0(:, i+1)
        rim1 = X0(:, i-1)
     end if

     ! Compute centered difference with respect to zeta (eq 14)
     rzeta = half*(rip1 - rim1)

     gamma = rzeta(1)**2 + rzeta(2)**2
     ovrGamma = one/gamma

     ! Compute derivative wrt to eta (normal direction) (eq 18)
     reta(1) = -rzeta(2)*Area0(i)/gamma
     reta(2) =  rzeta(1)*Area0(i)/gamma

     ! Compute the two coefficient matrices
     A0(1, 1, i) = reta(1)
     A0(1, 2, i) = reta(2)
     A0(2, 1, i) = reta(2)
     A0(2, 2, i) = -reta(1)

     B0(1, 1, i) = rzeta(1)
     B0(1, 2, i) = rzeta(2)
     B0(2, 1, i) = -rzeta(2)
     B0(2, 2, i) = rzeta(1)

  end do

end subroutine computeMetrics2D

subroutine computeStretch(l, deltaS)

  use hypInput

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: l
  
  ! Output Parameters
  real(kind=realType), intent(inout) :: deltaS

  ! Since we don't have a complex stretching function yet, we will
  ! just use geometric progression which is usually close to what we
  ! actually want in the marching direction anyway

  if (l == 2) then
     ! First step, set deltaS to the desired initial grid spacign
     ! given in HypInput
     deltaS = s0
  else
     ! Otherwise, just multiply the current deltaS by the gridRatio
     ! parameter, also in HypInput
     deltaS = deltaS*gridRatio
  end if

end subroutine computeStretch

subroutine assembleAndSolve(X0, X1, Area1, A0, B0, nx)
  
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(2, nx), Area1(nx)
  real(kind=realType), intent(in) :: A0(2, 2, nx), B0(2, 2, nx)
  integer(kind=intType), intent(in) :: nx
  real(kind=realType), intent(out) :: X1(2, nx)

  ! Working parameters
  integer(kind=intType) :: i, j, l, idim, ierr, idp1, idm1
  real(kind=realType) :: eye(2, 2), BInv(2,2), BInvA(2,2), deltaR(2)

  ! Define a 2x2 'I' matrix
  eye = zero
  eye(1,1) = one
  eye(2,2) = one

  ! Next we assemble A which is quite straight forward
  do i=1, nx

     ! Get the indices (fortran ordering) depending on whcih point we
     ! have
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
     
     ! Compute the inverse of B
     call two_by_two_inverse(B0(:, :, i), Binv)
     
     BinvA = matmul(Binv, A0(:, :, i))

     ! Point to the left
     call MatSetValuesBlocked(A, 1, i-1, 1, idm1-1, -(1+theta)*half*BinvA - epsI*eye, INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Center Point
     call MatSetValuesBlocked(A, 1, i-1, 1, i-1, eye + two*epsI*eye, INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)


     ! Point to the right
     call MatSetValuesBlocked(A, 1, i-1, 1, idp1-1, (1+theta)*half*BinvA - epsI*eye,  INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Finally assemble the RHS
     call VecSetValuesBlocked(RHS, 1, i-1, matmul(Binv,(/zero, Area1(i)/)), INSERT_VALUES, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do

  ! Assemble the matrix and vector
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyBegin(RHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecAssemblyEnd(RHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the KSP Object
  call KSPCreate(petsc_comm_world, kspObj, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetFromOptions(kspObj, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetOperators(kspObj, A, A, SAME_NONZERO_PATTERN, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetTolerances(kspObj, 1e-14, 1e-40, 1e100, 100, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now solve the system
  call KSPSolve(kspObj, RHS, rDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Take the solution and copy into X1, the next layer
  do i=1,nx
     call VecGetValues(rDelta, 2, (/2*(i-1), 2*(i-1)+1/), deltaR, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     X1(:, i) = X0(:, i) + deltaR
  end do
 
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

subroutine create2DPetscVars(nx)

  use hypData
  
  implicit none
  
  ! Input Variables
  integer(kind=intType) :: nx
  
  ! Working Variables
  integer(kind=intType) :: onProc(nx*2), offProc(nx*2), ierr
  
  ! Lets to things in the proper way, create the Mat first
  onProc = 3
  offProc = 1

  call MatCreateMPIAIJ(PETSC_COMM_WORLD,  &
       nx*2, nx*2, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, A, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetBlockSize(A, 2, ierr)
  
  ! This must be set to allow passing in blocks in native fortran
  ! ordering
  call MatSetOption(A, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Then use getVecs to get the vectors we want
  call MatGetVecs(A, rDelta, RHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
end subroutine create2DPetscVars

subroutine destroy2DPetscVars

  use hypData

  implicit none
  
  integer(kind=intType) :: ierr

  ! Now delete the objects we had create
  call KSPDestroy(kspObj, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatDestroy(A, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(rhs, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(rDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine destroy2DPetscVars

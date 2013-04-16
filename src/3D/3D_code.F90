subroutine init3d(sizes, nn, nPtr_in, lmax, nglobal)

  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: sizes(nn, 2), nn, lMax
  integer(kind=intType), intent(in) :: nPtr_in(lMax, nGlobal), nGlobal
  ! Working Parameters
  integer(kind=intType) :: i

  ! Set the number of patches based on the size
  nPatch = nn
 
  if (not(allocated(patches))) then
     allocate(patches(nPatch))
  end if
  
  do i=1,nPatch
     patches(i)%il = sizes(i, 1)
     patches(i)%jl = sizes(i, 2)
  end do

  if (not(allocated(nPtr))) then
     allocate(nptr(lmax, nGlobal))
  end if

  do i=1, nGlobal
     nPtr(:, i) = nPtr_in(:, i) 
  end do

end subroutine init3d

subroutine setlindex(ind_in, nx, ny, iPatch)

  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: ind_in(nx, ny), nx, ny, iPatch

  ! Working Parameters
  integer(kind=intType) :: i, j

  ! Set the ind_in into the proper patch l_index. Note we will convert
  ! to fortran ordering here by adding a 1

  if (not(allocated(patches(iPatch + 1)%l_index))) then
     allocate(patches(iPatch + 1)%l_index( &
     patches(iPatch + 1)%il, patches(iPatch + 1)%jl))
  end if

  do j=1, ny
     do i=1, nx
        patches(iPatch + 1)%l_index(i, j) = ind_in(i, j) + 1
     end do
  end do

end subroutine setlindex

subroutine run3D(Xin, nx, Xout)
  
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: Xin(3, nx)
  integer(kind=intType), intent(in) :: nx
  real(kind=realType), intent(out) :: Xout(3, nx)

  ! Working parameters
  integer(kind=intType) :: i, j, l, idim
  real(kind=realType), pointer, dimension(:, :) :: X0, X1, Xm1
  real(kind=realType) :: Volume0(nx), Volume1(nx), VBar, deltaS
  real(kind=realType) :: A0(3, 3, nx), B0(3, 3, nx), C0(3, 3, nx)
  real(kind=realType) :: sMax

  ! First thing we will do is allocate the final grid array:
  if (allocated(grid3D)) then
     deallocate(grid3D)
  end if

  allocate(grid3D(3, nx, NDebug))
  
  grid3D = zero

  ! Copy Xin into the first slot of grid3D
  do i=1,nx
     do idim=1,3
        grid3D(idim, i, 1) = Xin(idim, i)
     end do
  end do

  ! Compute an approximate estimate of the distance to farfield by
  ! taking s0 and computing
  sMax = s0*gridRatio**(N-1)

  ! Compute the original Volumes, based on initial s0 value
  X0 => grid3D(:, :, 1)

  ! This initializes the volumes for the corners that do not get a
  ! volume
  Volume0 = zero
  Volume1 = zero
  call computeVolumes(X0, Volume0, nx, VBar, s0)

  ! Create the petsc variables
  call create3DPetscVars(nx)

  ! Set the Xm1 pointer, to first layer. It isn't actually used, just
  ! good practice to point it somewhere as it is passed to a function.
  Xm1 => grid3D(:, :, 1) 

  ! This is the master marching direction loop
  marchLoop: do l=2, NDebug
     print *,'----------- Step: ', l

     ! Compute the scaled distance of this layer
     scaleDist = s0*gridRatio**(l-2)/sMax

     ! Set pointers to the two levels we are currently working with
     X0 => grid3D(:, :, l-1)
     X1 => grid3D(:, :, l  )
     if (l > 2) then
        Xm1 => grid3D(:, :, l-2)
     end if

     ! Compute the length increment in the marching direction. This is
     ! put in a separate function to allow for potential callbacks to
     ! user supplied functions if necessary
     call computeStretch(l, deltaS)

     ! Compute the nodal volumes on the X0 layer as well as the
     ! average volume VVar
     call computeVolumes(X0, Volume1, nx, Vbar, deltaS)

     ! Now perform the "Volume" smoothing operation.
     call VolumeSmooth(Volume1, nx, VBar, l)

     ! Compute the metrics
     call computeMetrics3D(X0, Volume0, Volume1, A0, B0, C0, nx, l)

     ! Assemble and solve
     call assembleAndSolve3d(X0, X1, Xm1, Volume1, A0, B0, C0, nx, l)

     ! Shuffle the volumes's backwards
     Volume0 = Volume1

  end do marchLoop

  ! Destroy the PETSc variables
  call destroyPetscVars

  Xout = X1

end subroutine run3D

subroutine computeVolumes(X, nodalVolume, nx, VBar, deltaS)

  use precision
  use hypData
  implicit none

  ! Compute the nodal volumes 

  ! Input Parameters
  real(kind=realType), intent(in) :: X(3, nx), deltaS
  integer(kind=intType), intent(in) :: nx
  
  ! Output Parameters
  real(kind=realType), intent(out) :: nodalVolume(nx), VBar
  
  ! Working Variables
  integer(kind=intType) :: i, j, jp1, jm1, kp1, km1, nn, ipatch
  real(kind=realType) :: deltaV, ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS

  ! Zero nodalVolume and VBar
  nodalVolume = zero
  VBar = zero

  ! Loop over the patfaces:
  do iPatch=1,nPatch
     ! Loop over faces on patch
     do j=1,patches(iPatch)%jl-1
        do i=1,patches(iPatch)%il-1

           ! Extract the coordinates of the 4 corners
           ll = X(:, patches(iPatch)%l_index(i  , j  ))
           ul = X(:, patches(iPatch)%l_index(i  , j+1))
           lr = X(:, patches(iPatch)%l_index(i+1, j  ))
           ur = X(:, patches(iPatch)%l_index(i+1, j+1))

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
           nodalVolume(patches(iPatch)%l_index(i, j)) = &
                nodalVolume(patches(iPatch)%l_index(i, j)) + areadS 
           nodalVolume(patches(iPatch)%l_index(i+1, j)) = &
                nodalVolume(patches(iPatch)%l_index(i+1, j)) + areadS 
           nodalVolume(patches(iPatch)%l_index(i, j+1)) = &
                nodalVolume(patches(iPatch)%l_index(i, j+1)) + areadS 
           nodalVolume(patches(iPatch)%l_index(i+1, j+1)) = &
                nodalVolume(patches(iPatch)%l_index(i+1, j+1)) + areadS 

           VBar = vBar + areadS*four
        end do
     end do
  end do
  VBar = VBar / nx


  ! ! Loop over the global nodes
  ! do i=1,nx

  !    if (.not. nptr(1,i) == 0) then
  !       ! Extract pointers to make code easier to read
  !       jp1 = nPtr(1, i)
  !       jm1 = nPtr(2, i)
  !       kp1 = nPtr(3, i)
  !       km1 = nPtr(4, i)

  !       ! Equation 25 in Gridgen's implementation of Hyperbolic Generation
  !       ! deltaV = deltaS*fourth*(dist(X(:, i), X(:, jm1)) + dist(X(:, i), X(:, jp1)))* &
  !       !      (dist(X(:, i), X(:, km1)) + dist(X(:, i), X(:, kp1)))

  !       !deltaV = deltaS*fourth*(dist(X(:, jp1), X(:, jm1))*dist(X(:, kp1), X(:, km1)))

  !       nodalVolume(i) = deltaV
  !       VBar = VBar + deltaV
    
  !       nn = nn + 1 ! Keep track of the actual number of volumes to
  !                   ! average
  !    end if
  ! end do
  
  ! Finally compute the average volume
           !VBar = VBar / nn

  contains 

    function dist(p1, p2)

      real(kind=realType) :: p1(3), p2(3)
      real(kind=realType) :: dist
      dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)

    end function dist
end subroutine computeVolumes

subroutine volumeSmooth(Volume, nx, VBar, l)

  use hypInput
  use hypData

  implicit none

  ! Perform averaging on the nodal volumes
  
  ! Input/Output Parameters
  real(kind=realType), intent(inout) :: Volume(nx)
  real(kind=realType), intent(in) :: VBar
  integer(kind=intType), intent(in) :: nx
  integer(kind=intType), intent(in) :: l

  ! Working Parameters
  real(kind=realType) :: factor
  integer(kind=intType) :: i, jp1, jm1, kp1, km1
  
  ! Just do a global volume smooth

  factor = half * (one + ( one - volBlend)**(l-2))
  Volume = factor*Volume + (one - factor)*Vbar

end subroutine volumeSmooth

subroutine computeMetrics3D(X0, Volume0, Volume1, A0, B0, C0, nx, l)

  use precision
  use hypData
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(3, nx), Volume0(nx), Volume1(nx)
  real(kind=realType), intent(out) :: A0(3, 3, nx), B0(3, 3, nx), C0(3, 3, nx)
  integer(kind=intType), intent(in) :: nx, l

  ! Working Variables
  real(kind=realType) :: r_ksi(3), r_eta(3), r_zeta(3)
  real(kind=realType) :: detC, deltaVovrDetC, r0(3), rm(3), thetam, s(3), v1(3), v2(3)
  integer(kind=intType) :: i, jp1, jm1, kp1, km1, m, MM
  real(kind=realType) :: d_avg

  do i=1, nx
     if (nptr(2, i) == 4) then
        ! Extract pointers to make code easier to read
        jp1 = nPtr(3, i)
        kp1 = nPtr(4, i)
        jm1 = nPtr(5, i)
        km1 = nPtr(6, i)

        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        r_ksi = half*(X0(:, jp1) - X0(:, jm1))
        r_eta  = half*(X0(:, kp1) - X0(:, km1))
     else

        ! We have a node with 6 neighbours.
        ! Lets do this brute force:
        r_ksi = zero
        r_eta = zero
        r0 = X0(:, i)
        MM = nPtr(2, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = X0(:, nptr(3 + m, i))
           r_ksi = r_ksi + (two/MM) * (rm - r0) * cos(thetam)
           r_eta = r_eta + (two/MM) * (rm - r0) * sin(thetam)
        end do
     end if
     
     v1 = r_ksi
     v2 = r_eta
     s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

     ! Equation 3.5 (Chen and Steger)
     detC = &
          (r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))**2 + &
          (r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))**2 + &
          (r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))**2
        
     deltaVovrDetC = Volume1(i)/detC
     
     r_zeta(1) = deltaVovrDetC*(r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3))
     r_zeta(2) = deltaVovrDetC*(r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3))
     r_zeta(3) = deltaVovrDetC*(r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2))

     ! Compute the two coefficient matrices
     A0(1, :, i) = r_zeta
     A0(2, :, i) = zero
     A0(3, 1, i) = r_eta(2)*r_zeta(3) - r_zeta(2)*r_eta(3)
     A0(3, 2, i) = r_zeta(1)*r_eta(3) - r_eta(1)*r_zeta(3)
     A0(3, 3, i) = r_eta(1)*r_zeta(2) - r_zeta(1)*r_eta(2)
     
     B0(1, :, i) = zero
     B0(2, :, i) = r_zeta
     B0(3, 1, i) = r_zeta(2)*r_ksi(3) - r_ksi(2)*r_zeta(3)
     B0(3, 2, i) = r_ksi(1)*r_zeta(3) - r_zeta(1)*r_ksi(3)
     B0(3, 3, i) = r_zeta(1)*r_ksi(2) - r_ksi(1)*r_zeta(2)
     
     C0(1, :, i) = r_ksi
     C0(2, :, i) = r_eta
     C0(3, 1, i) = r_ksi(2)*r_eta(3) - r_eta(2)*r_ksi(3)
     C0(3, 2, i) = r_eta(1)*r_ksi(3) - r_ksi(1)*r_eta(3)
     C0(3, 3, i) = r_ksi(1)*r_eta(2) - r_eta(1)*r_ksi(2)
     
  end do

end subroutine computeMetrics3D

subroutine assembleAndSolve3D(X0, X1, Xm1, Volume1, A0, B0, C0, nx, l)
  
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(3, nx), Volume1(nx), Xm1(3, nx)
  real(kind=realType), intent(in) :: A0(3, 3, nx), B0(3, 3, nx), C0(3, 3, nx)
  integer(kind=intType), intent(in) :: nx
  real(kind=realType), intent(out) :: X1(3, nx)

  ! Working parameters
  integer(kind=intType) :: j, l, idim, ierr
  integer(kind=intType) :: i, jp1, jm1, kp1, km1, nAverage, ii, m, MM
  real(kind=realType) :: eye(3, 3), CInv(3,3), CInvA(3,3), CInvB(3, 3)
  real(kind=realType) :: De(3), Sl, tmp, a_ksi, a_eta
  real(kind=realType) :: exp, R_ksi, R_eta, N_ksi, N_eta, deltaR(3)
  real(kind=realType) :: r0(3), rm(3), thetam, coef, d_ksi(nx), d_eta(nx),dmax
  ! Define a 3x3 'I' matrix
  eye = zero
  eye(1,1) = one
  eye(2,2) = one
  eye(3,3) = one
  call MatZeroEntries(hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecZeroEntries(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! ------------------------------------
  ! Scaling Function (Equation 6.5)
  ! ------------------------------------
  exp = 1.0
  Sl = ((dble(l)-2)/(N-1))**exp
  

  ! Precompute and store grid convergence sensors:
  d_ksi = zero
  d_eta = zero
  do i=1,nx
     if (nptr(2, i) == 4) then
        ! Extract pointers to make code easier to read
        jp1 = nPtr(3, i)
        kp1 = nPtr(4, i)
        jm1 = nPtr(5, i)
        km1 = nPtr(6, i)
     
        ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
        tmp = (dist(Xm1(:, jp1), Xm1(:, i)) + dist(Xm1(:, jm1), Xm1(:, i))) / &
             (dist(X0(:, jp1), X0(:, i)) + dist(X0(:, jm1), X0(:, i)))
        d_ksi(i) = max(tmp**(2/Sl), 0.1)
        
        tmp = (dist(Xm1(:, kp1), Xm1(:, i)) + dist(Xm1(:, km1), Xm1(:, i))) / &
             (dist(X0(:, kp1), X0(:, i)) + dist(X0(:, km1), X0(:, i)))
        d_eta(i) = max(tmp**(2/Sl), 0.1)
     end if
  end do

  print *,'max ksi, eta:',maxval(d_ksi), maxval(d_eta)

  ! Get grid divergence sensor from neighbours for extraordinary nodes
  do i=1,nx
     if (nptr(2, i) .ne. 4) then

        nAverage = nPtr(2, i) 
        ! Loop over neighbours
        d_ksi(i) = zero
        do ii = 1, nAverage
           d_ksi(i) = d_ksi(i) + d_ksi(nPtr(2+ii, i)) + d_eta(nPtr(2+ii, i))
        end do

        d_ksi(i) = max(0.1, d_ksi(i) / (two * nAverage))
        d_eta(i) = d_ksi(i)

     end if
  end do

  ! Next we assemble hypMat which is quite straight forward
  do i=1, nx

     if (nptr(2, i) == 4) then
        ! Extract pointers to make code easier to read

        jp1 = nPtr(3, i)
        jm1 = nPtr(5, i)
        kp1 = nPtr(4, i)
        km1 = nPtr(6, i)

        ! Compute the inverse of C and its multiplcation with A and B
        call three_by_three_inverse(C0(:, :, i), Cinv)
        CinvA = matmul(Cinv, A0(:, :, i))
        CinvB = matmul(Cinv, B0(:, :, i))

        ! -----------------
        ! Left Hand Side:
        ! -----------------

        ! Point to the left
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, jm1-1, &
             -(one+theta)*half*CinvA - epsI*eye, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Center Point ! Changed two to 4
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
             eye + four*epsI*eye, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Point to the right
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, jp1-1, &
             (one+theta)*half*CinvA - epsI*eye,  ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! ----------------------------------------------------------

        ! Point to the bottom
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, km1-1, &
             -(one+theta)*half*CinvB - epsI*eye, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Point to the right
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, kp1-1, &
             (one+theta)*half*CinvB - epsI*eye,  ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! -----------------
        ! Right Hand Side:
        ! -----------------

        ! Assemble the the Cinv comonent of the RHS
        call VecSetValuesBlocked(hypRHS, 1, i-1, &
             matmul(Cinv,(/zero, zero, Volume1(i)/)), ADD_VALUES, ierr)
    
        call EChk(ierr, __FILE__, __LINE__)

        ! ============================================
        !             Explicit Smoothing 
        ! ============================================

    
        ! ------------------------------------
        ! Compute the grid angle functions (Equation 6.9)
        ! ------------------------------------
        a_ksi = one
        a_eta = one
     
        ! ------------------------------------
        ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
        ! ------------------------------------
        R_ksi = Sl * d_ksi(i) * a_ksi
        R_eta = Sl * d_eta(i) * a_eta

        ! ------------------------------------
        ! Matrix Norm approximation  (Equation 6.3)
        ! ------------------------------------
        N_ksi = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
             (C0(1 ,1 ,i)**2 + C0(1, 2, i)**2 + C0(1, 3, i)**2))
        N_eta = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
             (C0(2 ,1 ,i)**2 + C0(2, 2, i)**2 + C0(2, 3, i)**2))

        ! ------------------------------------
        ! Final explict dissipation (Equation 6.1 and 6.2)
        ! ------------------------------------
        De = &
             epsE*R_ksi*N_ksi * (X0(:, jm1) - two*X0(:, i) + X0(:, jp1)) + &
             epsE*R_eta*N_eta * (X0(:, km1) - two*X0(:, i) + X0(:, kp1))
        
        if (l == 2) then
           De = zero
        end if

        call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     else

        ! ! Setting the matrix values are a little tricker here now. 

        ! ! (This is the same)
        ! call three_by_three_inverse(C0(:, :, i), Cinv)
        ! CinvA = matmul(Cinv, A0(:, :, i))
        ! CinvB = matmul(Cinv, B0(:, :, i))
        ! MM = nPtr(2, i) 
        ! do m = 0, MM - 1
        !    thetam = 2*pi*m/MM
           
        !    ! For the 'fm' point in ksi
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(3+m, i)-1, &
        !         (one + theta)*CinvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'f0' point in ksi
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !         -(one + theta)*cInvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'fm' point in eta
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(3+m, i)-1, &
        !         (one + theta)*CinvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'f0' point in ksi
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !         -(one + theta)*cInvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! Now we have the second derivative values to do in ksi-ksi

        !    ! ! For the 'fm' point in ksi-ksi 
        !    coef = (two/MM)*(four*cos(thetam)**2 - one)

        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(3+m, i)-1, &
        !         -eye*epsI*coef, ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'f0' point in ksi-ksi
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !         eye*epsI*coef, ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'fm' point in eta-eta
        !    coef = (two/MM)*(four*sin(thetam)**2 - one)

        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(3+m, i)-1, &
        !         -eye*epsI*coef, ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)

        !    ! For the 'f0' point in eta-eta
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !         eye*epsI*coef, ADD_VALUES, ierr)
        !    call EChk(ierr, __FILE__, __LINE__)
        ! end do

        ! ! Finally we need an eye on the diagonal
        ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !      eye, ADD_VALUES, ierr)
        ! call EChk(ierr, __FILE__, __LINE__)

        ! ! -----------------
        ! ! Right Hand Side:
        ! ! -----------------

        ! ! Assemble the the Cinv comonent of the RHS
        ! call VecSetValuesBlocked(hypRHS, 1, i-1, &
        !      matmul(Cinv,(/zero, zero, Volume1(i)/)), ADD_VALUES, ierr)
         
        ! call EChk(ierr, __FILE__, __LINE__)
        
        ! ! ============================================
        ! !             Explicit Smoothing 
        ! ! ============================================
        
        ! ! ------------------------------------
        ! ! Compute the grid angle functions (Equation 6.9)
        ! ! ------------------------------------
        ! a_ksi = one
        ! a_eta = one
           
        ! ! ------------------------------------
        ! ! Combine into 'R_ksi' and 'R_eta' (Equation 6.4)
        ! ! ------------------------------------
        ! R_ksi = Sl * d_ksi(i) * a_ksi
        ! R_eta = Sl * d_eta(i) * a_eta
        
        ! ! ------------------------------------
        ! ! Matrix Norm approximation  (Equation 6.3)
        ! ! ------------------------------------
        ! N_ksi = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
        !      (C0(1 ,1 ,i)**2 + C0(1, 2, i)**2 + C0(1, 3, i)**2))
        ! N_eta = sqrt((A0(1, 1, i)**2 + A0(1, 2, i)**2 + A0(1, 3, i)**2)/ &
        !      (C0(2 ,1 ,i)**2 + C0(2, 2, i)**2 + C0(2, 3, i)**2))
        
        ! ! ------------------------------------
        ! ! Final explict dissipation (Equation 6.1 and 6.2)
        ! ! ------------------------------------
        ! de = zero
        ! r0 = X0(:, i)
        ! MM = nPtr(2, i) 
        ! do m = 0, MM - 1
        !    thetam = 2*pi*m/MM
        !    rm = X0(:, nptr(3 + m, i))

        !    De = De + &
        !         epsE*R_ksi*N_ksi*(two/MM)*(rm - r0)*(four*cos(thetam)**2 - one) + &
        !         epsE*R_eta*N_eta*(two/MM)*(rm - r0)*(four*sin(thetam)**2 - one) 
        ! end do
        
        ! if (l == 2) then
        !    De = zero
        ! end if
        
        ! call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
        ! call EChk(ierr, __FILE__, __LINE__)
        
        !We have to do averaging at the corners
        
        ! ! Get the number of nodes to average:
        nAverage = nPtr(2, i) 

        ! ! Loop over neighbours
        do ii = 1, nAverage
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+ii, i)-1, &
                -one/nAverage*eye, ADD_VALUES, ierr)
        end do

        ! Center Point
        call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
             eye, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Don't set a RHS value for this, it is zero. This is already
        ! taken care of with the vecSet above


     end if
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
  if (mod(l,10) == 0) then
     call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Now solve the system
  call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Take the solution in deltaR and add to X0 to get X1
  do i=1,nx
     call VecGetValues(hypDelta, 3, (/3*(i-1), 3*(i-1)+1, 3*(i-1)+2/), &
          deltaR, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     X1(:, i) = X0(:, i) + deltaR
  end do
 
  contains 

    function dist(p1, p2)

      real(kind=realType) :: p1(3), p2(3)
      real(kind=realType) :: dist
      dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)

    end function dist
  end subroutine assembleAndSolve3D

subroutine create3DPetscVars(nx)

  use hypData
  
  implicit none
  
  ! Input Variables
  integer(kind=intType) :: nx
  
  ! Working Variables
  integer(kind=intType) :: onProc(nx*3), offProc(nx*3), ierr

  ! ----------------------------------------------------------
  !          Linearized Hyperbolic System Variables
  ! ----------------------------------------------------------

  ! Lets to things in the proper way, create the Mat first
  onProc = 35
  offProc = 1

  call MatCreateMPIAIJ(PETSC_COMM_WORLD,  &
       nx*3, nx*3, PETSC_DETERMINE, PETSC_DETERMINE, &
       0, onProc, 0, offProc, hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Make it a bocked matrix of size 3
  call MatSetBlockSize(hypMat, 3, ierr)
  
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

  call KSPGMRESSetRestart(hypKSP, 250, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetTolerances(hypKSP, 1e-10, 1e-16, 1e5, 500, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine create3DPetscVars

subroutine three_by_three_inverse(Jac, Jinv)

  use precision
  implicit none

  real(kind=realType), intent(in) :: Jac(3, 3)
  real(kind=realType), intent(out) :: Jinv(3, 3)
  real(kind=realType) :: invdet, det

  ! Express the inverse of the jacobian explicitly 
  det = Jac(1, 1)*Jac(2, 2)*Jac(3, 3)  &
       +   Jac(1, 2)*Jac(2, 3)*Jac(3, 1)&
       +   Jac(1, 3)*Jac(2, 1)*Jac(3, 2)&
       -   Jac(1, 1)*Jac(2, 3)*Jac(3, 2)&
       -   Jac(1, 2)*Jac(2, 1)*Jac(3, 3)&
       -   Jac(1, 3)*Jac(2, 2)*Jac(3, 1)
  invdet = one/det

  !             [ |a22 a23|   |a12 a13|  |a12 a13|]     */
  !             [ |a32 a33|  -|a32 a33|  |a22 a23|]     */
  !             [                                 ]     */
  !             [ |a21 a23|   |a11 a13|  |a11 a13|]     */
  !    A^(-1) = [-|a31 a33|   |a31 a33| -|a21 a23|] / d */
  !             [                                 ]     */
  !             [ |a21 a22|   |a11 a12|  |a11 a12|]     */
  !             [ |a31 a32|  -|a31 a32|  |a21 a22|]     */

  Jinv(1, 1) =  invdet*(Jac(2, 2)*Jac(3, 3)-Jac(2, 3)*Jac(3, 2))
  Jinv(1, 2) = -invdet*(Jac(1, 2)*Jac(3, 3)-Jac(1, 3)*Jac(3, 2))
  Jinv(1, 3) =  invdet*(Jac(1, 2)*Jac(2, 3)-Jac(1, 3)*Jac(2, 2))

  Jinv(2, 1) = -invdet*(Jac(2, 1)*Jac(3, 3)-Jac(2, 3)*Jac(3, 1))
  Jinv(2, 2) =  invdet*(Jac(1, 1)*Jac(3, 3)-Jac(1, 3)*Jac(3, 1))
  Jinv(2, 3) = -invdet*(Jac(1, 1)*Jac(2, 3)-Jac(1, 3)*Jac(2, 1))

  Jinv(3, 1) =  invdet*(Jac(2, 1)*Jac(3, 2)-Jac(2, 2)*Jac(3, 1))
  Jinv(3, 2) = -invdet*(Jac(1, 1)*Jac(3, 2)-Jac(1, 2)*Jac(3, 1))
  Jinv(3, 3) =  invdet*(Jac(1, 1)*Jac(2, 2)-Jac(1, 2)*Jac(2, 1))

end subroutine three_by_three_inverse

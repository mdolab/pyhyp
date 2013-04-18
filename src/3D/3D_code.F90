subroutine init3d(sizes, nn, nPtr_in, lmax, nglobal)

  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: sizes(nn, 2), nn, lMax
  integer(kind=intType), intent(in) :: nPtr_in(lMax, nGlobal), nGlobal
  ! Working Parameters
  integer(kind=intType) :: i

  ! Set the number of patches:
  nPatch = nn
 
  ! Allocate patch list of patch types
  if (not(allocated(patches))) then
     allocate(patches(nPatch))
  end if
  
  ! Set sizes
  do i=1,nPatch
     patches(i)%il = sizes(i, 1)
     patches(i)%jl = sizes(i, 2)
  end do

  ! Allocate the global node pointer structure
  if (not(allocated(nPtr))) then
     allocate(nptr(lmax, nGlobal))
  end if

  ! And copy in our python data
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
  integer(kind=intType) :: i, j, l, idim, ll
  real(kind=realType), pointer, dimension(:, :) :: X0, X1, Xm1
  real(kind=realType) :: Volume0(nx), Volume1(nx), VBar
  real(kind=realType) :: cmax

  ! First thing we will do is allocate the final grid array, and zero
  if (allocated(pGrid3D)) then
     deallocate(pGrid3D)
  end if
  allocate(pGrid3D(3, nx, NMax))
  pGrid3D = zero
  call writeHeader
  ! Copy Xin into the first slot of pGrid3D
  do i=1,nx
     do idim=1,3
        pGrid3D(idim, i, 1) = Xin(idim, i)
     end do
  end do

  ! Compute the original Volumes, based on initial s0 value
  X0 => pGrid3D(:, :, 1)

  ! Compute the inital 'Radius' value, R0
  call computeR(X0, nx, radius0)

  ! This initializes the volumes for the corners that do not get a
  ! volume
  Volume0 = zero
  Volume1 = zero
  call computeVolumes(X0, Volume0, nx, VBar, cmax)

  ! Create the petsc variables
  call create3DPetscVars(nx)

  ! Set the Xm1 pointer, to first layer. It isn't actually used, just
  ! good practice to point it somewhere as it is passed to a function.
  Xm1 => pGrid3D(:, :, 1) 
  limitStep = .True.
  factorNext = .True.

  ! Determine starting CPU Time
  timeStart = mpi_wtime()

  ! This is the master marching direction loop
  marchLoop: do marchIter=2, Nmax
     
     ! Use 'l' as the working variable
     l = marchIter

     ! Set pointers to the two levels we are currently working with
     X0 => pGrid3D(:, :, l-1)
     X1 => pGrid3D(:, :, l  )
     if (l > 2) then
        Xm1 => pGrid3D(:, :, l-2)
     end if

     ! Compute the length increment in the marching direction. This is
     ! put in a separate function to allow for potential callbacks to
     ! user supplied functions if necessary

     call computeStretch(l)

     ! Compute the nodal volumes on the X0 layer as well as the
     ! average volume VVar
     call computeVolumes(X0, Volume1, nx, Vbar, cmax)
     if (limitStep) then
        if (cmax > three) then
           deltaS = deltas*.5
           call computeVolumes(X0, Volume1, nx, Vbar, cmax)
        end if
     else
        if (cmax > 8.0) then
           deltaS = deltas*.5
           call computeVolumes(X0, Volume1, nx, Vbar, cmax)
        end if
     end if

     ! Now perform the "Volume" smoothing operation.
     call VolumeSmooth(Volume1, nx, VBar, l)

     ! Assemble and solve
     call assembleAndSolve3d(X0, X1, Xm1, Volume1, nx, l)
     
     ! Compute the max radius of the current level, X1
     call computeMinR(X1, nx, Radius)

     ! Shuffle the volumes's backwards
     Volume0 = Volume1

     ! Possibly write header and write iteration info
     if (mod(marchIter, 50) == 0) then
        call writeHeader
     end if

     ! Check the quality of this layer
     call computeQualityLayer(X0, X1, nx)

     call writeIteration

     ! Now check if we have gone far enough:
     if (radius/radius0 > rmin) then
        nLayers = l
        exit marchLoop
     end if
     
     if (l == Nmax) then
        Nlayers = l
     end if

  end do marchLoop

  ! Now Process out the actual desired grid 
  call reparameterize(nx)

  ! Destroy the PETSc variables
  call destroyPetscVars

  Xout = X1

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

subroutine computeVolumes(X, nodalVolume, nx, VBar, cmax)

  use precision
  use hypData
  implicit none

  ! Compute the nodal volumes 

  ! Input Parameters
  real(kind=realType), intent(in) :: X(3, nx)
  integer(kind=intType), intent(in) :: nx
  
  ! Output Parameters
  real(kind=realType), intent(out) :: nodalVolume(nx), VBar
  
  ! Working Variables
  integer(kind=intType) :: i, j, jp1, jm1, kp1, km1, nn, ipatch
  real(kind=realType) :: deltaV, ll(3), ul(3), lr(3), ur(3)
  real(kind=realType) :: v1(3), v2(3), s(3), areadS
  real(kind=realType) :: c1,c2,c3,c4,cmax

  ! Zero nodalVolume and VBar
  nodalVolume = zero
  VBar = zero
  cmax = zero
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

           ! Also compute edge lengths
           c1 = deltaS/dist(ll,lr)
           c2 = deltaS/dist(ul,ur)
           c3 = deltaS/dist(ll,ul)
           c4 = deltaS/dist(lr,ur)

           cmax = max(cmax,c1,c2,c2,c4)

        end do
     end do
  end do
  VBar = VBar / nx

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
  real(kind=realType) :: factor, Vtmp(nx), oneOvrNNeighbour
  integer(kind=intType) :: i, ii, iter

  ! Do a Jacobai volume smooth
  Vtmp = Volume
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
     Vtmp = Volume
  end do
  
  ! Now we do a global volume smooth

  factor = half * (one + ( one - volBlend)**(l-2))
  Volume = factor*Volume + (one - factor)*Vbar

end subroutine volumeSmooth

subroutine assembleAndSolve3D(X0, X1, Xm1, Volume1, nx)
  
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  real(kind=realType), intent(in) :: X0(3, nx), Volume1(nx), Xm1(3, nx)
  integer(kind=intType), intent(in) :: nx
  real(kind=realType), intent(out) :: X1(3, nx)

  ! Working parameters
  real(kind=realType) :: A0(3, 3, nx), B0(3, 3, nx), C0(3, 3, nx)
  integer(kind=intType) :: j, l, idim, ierr
  integer(kind=intType) :: i, jp1, jm1, kp1, km1, nAverage, ii, m, MM
  real(kind=realType) :: eye(3, 3), CInv(3,3), CInvA(3,3), CInvB(3, 3)
  real(kind=realType) :: De(3), Sl, tmp, a_ksi, a_eta
  real(kind=realType) :: exp, Rksi, Reta, N_ksi, N_eta, deltaR(3)
  real(kind=realType) :: r0(3), rm(3), thetam, coef, d_ksi(nx), d_eta(nx),dmax
  real(kind=realType) :: r_ksi(3), r_eta(3), r_zeta(3), v1(3), v2(3), s(3), detC
  real(kind=realType) :: deltaVovrDetC
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
  ! Scaling Function (Equation 6.5) -- Modified
  ! ------------------------------------

  ! Since we don't acutally know how many iterations it will take to
  ! get to the farfield, we really don't want to use an N scaling. Sl
  ! is really quite important, but only near the boundary layer such
  ! that grid orthogonality is maintained. 

  if (l < 50) then
     sl = (2-l)/(50)
  else
     sl = one
  end if


  ! Precompute and store grid metrics and convergence sensors:
  d_ksi = zero
  d_eta = zero
  do i=1,nx
     if (nptr(1, i) == 4) then
        ! Extract pointers to make code easier to read
        jp1 = nPtr(2, i)
        kp1 = nPtr(3, i)
        jm1 = nPtr(4, i)
        km1 = nPtr(5, i)
     
        ! Compute centered difference with respect to ksi and eta. This is the
        ! nominal calculation that we will use in the interior of the domain
        r_ksi = half*(X0(:, jp1) - X0(:, jm1))
        r_eta  = half*(X0(:, kp1) - X0(:, km1))

        ! Compute the grid distribution sensor (eq 6.7 Chen and Steger)
        tmp = (dist(Xm1(:, jp1), Xm1(:, i)) + dist(Xm1(:, jm1), Xm1(:, i))) / &
             (dist(X0(:, jp1), X0(:, i)) + dist(X0(:, jm1), X0(:, i)))
        d_ksi(i) = max(tmp, 0.1)
        
        tmp = (dist(Xm1(:, kp1), Xm1(:, i)) + dist(Xm1(:, km1), Xm1(:, i))) / &
             (dist(X0(:, kp1), X0(:, i)) + dist(X0(:, km1), X0(:, i)))
        d_eta(i) = max(tmp, 0.1)
     else
        ! We have a node with 6 neighbours.
        ! Lets do this brute force:
        r_ksi = zero
        r_eta = zero
        r0 = X0(:, i)
        MM = nPtr(1, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = X0(:, nptr(2 + m, i))
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

  if (maxval(d_ksi) < .999 .and. maxval(d_eta) < .999) then
     limitStep = .False.
  else
     limitStep =.true.
  end if

  ! Get grid divergence sensor from neighbours for extraordinary nodes
  do i=1,nx
     if (nptr(1, i) .ne. 4) then

        nAverage = nPtr(1, i) 
        ! Loop over neighbours
        d_ksi(i) = zero
        do ii = 1, nAverage
           d_ksi(i) = d_ksi(i) + d_ksi(nPtr(1+ii, i)) + d_eta(nPtr(1+ii, i))
        end do

        d_ksi(i) = max(0.1, d_ksi(i) / (two * nAverage))
        d_eta(i) = d_ksi(i)
     end if
  end do
  ! Record the maximum and minimum grid sensors
  gridSensorMax = max(maxval(d_ksi), maxval(d_eta))
  gridSensorMin = min(minval(d_ksi), minval(d_eta))

  ! Next we assemble hypMat which is quite straight forward
  do i=1, nx

     if (nptr(1, i) == 4) then
        ! Extract pointers to make code easier to read

        jp1 = nPtr(2, i)
        kp1 = nPtr(3, i)
        jm1 = nPtr(4, i)
        km1 = nPtr(5, i)

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
        Rksi = Sl * d_ksi(i) * a_ksi
        Reta = Sl * d_eta(i) * a_eta

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
             epsE*Rksi*N_ksi * (X0(:, jm1) - two*X0(:, i) + X0(:, jp1)) + &
             epsE*Reta*N_eta * (X0(:, km1) - two*X0(:, i) + X0(:, kp1))
        
        if (l == 2) then
           De = zero
        end if

        call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     else

        ! Setting the matrix values are a little tricker here now. 

        ! (This is the same)
        call three_by_three_inverse(C0(:, :, i), Cinv)
        CinvA = matmul(Cinv, A0(:, :, i))
        CinvB = matmul(Cinv, B0(:, :, i))
        MM = nPtr(1, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           
           ! For the 'fm' point in ksi
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                (one + theta)*CinvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! For the 'f0' point in ksi
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                -(one + theta)*cInvA*(two/MM)*cos(thetam), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! For the 'fm' point in eta
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(2+m, i)-1, &
                (one + theta)*CinvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! For the 'f0' point in ksi
           call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
                -(one + theta)*cInvB*(two/MM)*sin(thetam), ADD_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Now we have the second derivative values to do in ksi-ksi

           ! ! For the 'fm' point in ksi-ksi 
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
        Rksi = Sl * d_ksi(i) * a_ksi
        Reta = Sl * d_eta(i) * a_eta
        
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
        de = zero
        r0 = X0(:, i)
        MM = nPtr(1, i) 
        do m = 0, MM - 1
           thetam = 2*pi*m/MM
           rm = X0(:, nptr(2 + m, i))

           De = De + &
                epsE*Rksi*N_ksi*(two/MM)*(rm - r0)*(four*cos(thetam)**2 - one) + &
                epsE*Reta*N_eta*(two/MM)*(rm - r0)*(four*sin(thetam)**2 - one) 
        end do
        
        if (l == 2) then
           De = zero
        end if
        
        call VecSetValuesBlocked(hypRHS, 1, i-1, De, ADD_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        !We have to do averaging at the corners
        
        ! ! ! ! Get the number of nodes to average:
        ! nAverage = nPtr(1, i) 

        ! ! ! Loop over neighbours
        ! do ii = 1, nAverage
        !    call MatSetValuesBlocked(hypMat, 1, i-1, 1, nPtr(1+ii, i)-1, &
        !         -one/nAverage*eye, ADD_VALUES, ierr)
        ! end do

       ! ! Center Point
        ! call MatSetValuesBlocked(hypMat, 1, i-1, 1, i-1, &
        !      eye, ADD_VALUES, ierr)
        ! call EChk(ierr, __FILE__, __LINE__)

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
  if (mod(l,preConLag) == 0 .or. factorNext) then
     call KSPSetOperators(hypKSP, hypMat, hypMat, SAME_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     factorNext = .False.
  end if

  ! Now solve the system
  call KSPSolve(hypKSP, hypRHS, hypDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGetIterationNumber(hypKsp, kspIts, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! If we did more iterations that the size of the subspace, flag the
  ! next iteration to refactor the PC.
  if (kspIts > kspsubspacesize) then
     factorNext = .True.
  end if
  
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

  use hypInput
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

  ! Create the KSP Object
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

subroutine writeHeader

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Writes the header to stdout for the output
  !
  use precision
  implicit none
  
  ! Working
  integer(kind=intType) :: i

  ! Write the first line of '-'
  write(*,"(a)",advance="no") "#"
  do i=2,79
     write(*,"(a)",advance="no") "-"
  enddo
  print "(1x)"

  write(*,"(a)",advance="no") "# Grid  |     CPU    |  KSP  | Grid       | Grid       |     Min    |   deltaS   |   min R   | "
  print "(1x)"
  write(*,"(a)",advance="no") "# Level |     Time   |  Its  | Sensor Max | Sensor Min |  Quality   |            |           | "
  print "(1x)"

  ! Write the Last line of '-'
  write(*,"(a)",advance="no") "#"
  do i=2,79
     write(*,"(a)",advance="no") "-"
  enddo
  print "(1x)"

end subroutine writeHeader

subroutine writeIteration
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Writes the information pertaining to the curent iteration to
  !     the screen
  !
  use hypData

  implicit none

  ! Iteration Count
  write(*,"(i8,1x)",advance="no") marchIter
  
  ! CPU time
  write(*,"(e12.5,1x)",advance="no") mpi_wtime() - timeStart

  ! KSP ITerations
  write(*,"(i6,1x)",advance="no") kspIts

  ! maximum Grid convergence sensor
  write(*,"(e12.5,1x)",advance="no") gridSensorMax

  ! maximum Grid convergence sensor
  write(*,"(e12.5,1x)",advance="no") gridSensorMin

  ! minimum grid quality for new layer
  write(*,"(e12.5,1x)",advance="no") minQuality

  ! marching step size for this iteration
  write(*,"(e12.5,1x)",advance="no") deltaS

  ! approximate minimim distance to body
  minR = zero
  write(*,"(e12.5,e12.5,e12.5,1x)",advance="no") radius/radius0
  print "(1x)"

end subroutine writeIteration
  
subroutine computeR(X, nx, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     computeR0 determines the minimum radius of a sphere that will
  !     fully enclose the initial grid.
  !
  use hypData
  implicit none

  ! Input parameters
  real(kind=realType), intent(in) :: X(3, nx)
  integer(kind=intType), intent(in) :: nx

  ! Ouput parameters
  real(kind=realType), intent(out) :: R
  
  ! Working
  integer(kind=intType) :: i
  
  ! First find the average of X
  Xavg = zero
  do i=1,nx
     Xavg = Xavg + X(:, i)
  end do

  Xavg = Xavg/nx

  ! Now find the maximum distance of X from Xavg
  R = zero
  do i=1,nx
     R = max(R, dist(Xavg, X(:, i)))
  end do

contains
  
  function dist(p1, p2)
    
    real(kind=realType) :: p1(3), p2(3)
    real(kind=realType) :: dist
    dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)
    
  end function dist
end subroutine computeR

subroutine computeMinR(X, nx, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     computeMinR determine the cloest point of X to Xavg
  !
  use hypData
  implicit none

  ! Input parameters
  real(kind=realType), intent(in) :: X(3, nx)
  integer(kind=intType), intent(in) :: nx

  ! Ouput parameters
  real(kind=realType), intent(out) :: R
  
  ! Working
  integer(kind=intType) :: i
  
  ! Now find the minimum distance of X from Xavg
  R = huge(R)
  do i=1,nx
     R = min(R, dist(Xavg, X(:, i)))
  end do

contains
  
  function dist(p1, p2)
    
    real(kind=realType) :: p1(3), p2(3)
    real(kind=realType) :: dist
    dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)
    
  end function dist
end subroutine computeMinR

subroutine computeQualityLayer(X0, X1, nx)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Compute the minimum quality measure for the cells defined by
  !     two grid levels
  !
  use hypInput
  use hypData
  implicit none
  
  ! Input parameters
  real(kind=realType), intent(in) :: X0(3, nx), X1(3, nx)
  integer(kind=intType), intent(in) :: nx

  ! Working
  integer(kind=intType) :: i, j, iPatch
  real(kind=realType) :: points(3,8), Q

  ! Since we are dealing with volumes, we need to loop over faces:
  minQuality = one
  do iPatch=1, nPatch
     do j=1,patches(iPatch)%jl-1
        do i=1,patches(iPatch)%il-1
           ! Assemble the 8 points we need for the volume - Note
           ! coodinate coordinate ordering!
           points(:,1) = X0(:,patches(iPatch)%l_index(i  ,j  ))
           points(:,2) = X0(:,patches(iPatch)%l_index(i+1,j  ))
           points(:,3) = X0(:,patches(iPatch)%l_index(i  ,j+1))
           points(:,4) = X0(:,patches(iPatch)%l_index(i+1,j+1))

           points(:,5) = X1(:,patches(iPatch)%l_index(i  ,j  ))
           points(:,6) = X1(:,patches(iPatch)%l_index(i+1,j  ))
           points(:,7) = X1(:,patches(iPatch)%l_index(i  ,j+1))
           points(:,8) = X1(:,patches(iPatch)%l_index(i+1,j+1))
           
           call quality_hexa(points, Q)
            minQuality = min(minQuality, Q)

        end do
     end do
  end do

end subroutine computeQualityLayer



subroutine quality_hexa(points, quality)

  use precision
  implicit none
  
  ! Input/Output
  real(kind=realType), intent(in) :: points(3, 8)
  real(kind=realType), intent(out) :: quality

  ! Working
  integer(kind=intType), parameter :: Ng=2
  real(kind=realType) :: det(Ng*Ng*Ng)
  integer(kind=intType):: i, k, l, m
  real(kind=realType) :: r, s, t
   real(kind=realType) :: largest, d_max, corners(Ng)
  real(kind=realType) :: XX(3), XXd(9), N(8), Na(8), Nb(8), Nc(8), Ua(9)
  real(kind=realType) :: JJac(9), h, pt(3)
    ! Relative Determinant: the ratio of the smallest determinat of the
  ! jacobian matrix divided by the largest determinate of the jacobian
  ! matrix using a 3x3x3 stecil

  corners(1) = -1
  corners(2) = 1

  ! 2x2x2 Stencil
  do k=1, Ng ! Gauss in r
     do l=1, Ng ! Gauss in s
        do m=1, Ng ! Gauss in t
           pt(1) = corners(k)
           pt(2) = corners(l)
           pt(3) = corners(m)

           ! Calculate the Jacobian and shape functions
           call SolidJacobian(XX, XXd, N, Na, Nb, Nc, pt, points)
           call jacobian3d(XXd, JJac, h)

           det((k-1)*Ng*Ng + (l-1)*Ng + m) = h
           
        end do
     end do
  end do

  ! We don't use maxval since we can't cs it
  !largest = maxval(abs(det))
  
  largest = abs(det(1))
  do i=1, Ng*Ng*Ng
     if (abs(det(i)) > largest) then
        largest = abs(det(i))
     end if
  end do

  d_max = abs(largest-det(1))
  do i=1, Ng*Ng*Ng
     if (abs(largest-det(i)) > d_max) then
        d_max = abs(largest-det(i))
     end if
  end do
     
  quality = 1-d_max/largest

end subroutine quality_hexa

subroutine SolidJacobian(X, Xd, N, Na, Nb, Nc, pt, Xpts)

  use precision
  implicit none

  ! Arguments:
  ! X: Point evaluated at parametric location 'pt'
  ! Xd: Derivative evaluated at parametric location 'pt'
  ! N : Shape function at 'pt'
  ! Na: 'r' shape function derivative at 'pt'
  ! Nb: 's' shape function derivative at 'pt'
  ! Nc: 't' shape function derivative at 'pt'
  ! pt: Parametric location to evalulate
  ! Xps: The list of corners defining the element

  ! Input/Output
  real(kind=realType), intent(out) :: X(3), Xd(9)
  real(kind=realType), intent(out) :: N(8), Na(8), Nb(8), Nc(8)
  real(kind=realType), intent(in) :: pt(3), Xpts(3, 8)

  ! Local
  integer(kind=intType) :: i, j, k, counter
  real(kind=realType) :: nna(2), nnb(2), nnc(2), dna(2), dnb(2), dnc(2)

  X = zero
  Xd = zero

  call lagrangeSF(nna, dna, pt(1), 2)
  call lagrangeSF(nnb, dnb, pt(2), 2)
  call lagrangeSF(nnc, dnc, pt(3), 2)

  do k=1, 2
     do j=1, 2
        do i=1, 2
           counter = (k-1)*4 + (j-1)*2 + i

           N(counter) = nna(i)*nnb(j)*nnc(k)
           X(1) = X(1) + Xpts(1, counter)*N(counter)
           X(2) = X(2) + Xpts(2, counter)*N(counter)
           X(3) = X(3) + Xpts(3, counter)*N(counter)

           Na(counter) = dna(i)*nnb(j)*nnc(k)
           Xd(1) = Xd(1) + Xpts(1, counter)*Na(counter)
           Xd(2) = Xd(2) + Xpts(2, counter)*Na(counter)
           Xd(3) = Xd(3) + Xpts(3, counter)*Na(counter)

           Nb(counter) = nna(i)*dnb(j)*nnc(k)
           Xd(4) = Xd(4) + Xpts(1, counter)*Nb(counter)
           Xd(5) = Xd(5) + Xpts(2, counter)*Nb(counter)
           Xd(6) = Xd(6) + Xpts(3, counter)*Nb(counter)

           Nc(counter) = nna(i)*nnb(j)*dnc(k)
           Xd(7) = Xd(7) + Xpts(1, counter)*Nc(counter)
           Xd(8) = Xd(8) + Xpts(2, counter)*Nc(counter)
           Xd(9) = Xd(9) + Xpts(3, counter)*Nc(counter)

        end do
     end do
  end do

end subroutine SolidJacobian

subroutine jacobian3d(Xd, Jinv, h)

  use precision
  implicit none

  ! Input/Output
  real(kind=realType), intent(in)  :: Xd(9)
  real(kind=realType), intent(out) :: Jinv(9), h

  ! Working
  real(kind=realType) :: hinv

  h = (Xd(9)*(Xd(1)*Xd(5) - Xd(4)*Xd(2)) &
       - Xd(8)*(Xd(1)*Xd(6) - Xd(4)*Xd(3)) &
       + Xd(7)*(Xd(2)*Xd(6) - Xd(3)*Xd(5)))

  hinv = one/h 

  Jinv(1) =   (Xd(5)*Xd(9) - Xd(6)*Xd(8))*hinv
  Jinv(2) = - (Xd(2)*Xd(9) - Xd(3)*Xd(8))*hinv
  Jinv(3) =   (Xd(2)*Xd(6) - Xd(3)*Xd(5))*hinv

  Jinv(4) = - (Xd(4)*Xd(9) - Xd(6)*Xd(7))*hinv
  Jinv(5) =   (Xd(1)*Xd(9) - Xd(3)*Xd(7))*hinv
  Jinv(6) = - (Xd(1)*Xd(6) - Xd(3)*Xd(4))*hinv

  Jinv(7) =   (Xd(4)*Xd(8) - Xd(5)*Xd(7))*hinv
  Jinv(8) = - (Xd(1)*Xd(8) - Xd(2)*Xd(7))*hinv
  Jinv(9) =   (Xd(1)*Xd(5) - Xd(2)*Xd(4))*hinv

end subroutine jacobian3d


subroutine lagrangeSF( sf, dsf, a, order)

  use precision
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in):: order
  real(kind=realType), intent(out) :: sf(order), dsf(order)
  real(kind=realType), intent(in)  :: a

  ! Working
  integer(kind=intType) :: i, j, k
  real(kind=realType) :: ki, kj, dn, kk

  sf = one
  dsf = zero

  ! Loop over the shape function control points
  do i=0, order-1
     ki = -one + two*i/(order - one)

     ! Loop over each point again, except for the current control point, 
     ! adding the contribution to the shape function

     do j=0, order-1
        if (i .ne. j) then
           kj = -one + two*j/(order-one)
           sf(i+1) = sf(i+1)*(a-kj)/(ki-kj)


           ! Loop over the whole thing again to determine the contribution to
           ! the derivative of the shape function
           dn = one/(ki - kj)

           do k=0, order-1
              if ( k .ne. i .and. k .ne. j ) then
                 kk = -one + two*k/(order -one)
                 dn = dn*(a - kk)/(ki - kk)
              end if
           end do

           dsf(i+1) = dsf(i+1) + dn

        end if
     end do
  end do

end subroutine lagrangeSF


subroutine reparameterize(nx)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     reparameterize() takes the pseudo grid defined in pgrid3d and
  !     computes the final desired distribution according to N, s0,
  !     and gridRatio. (rMin was already used since that determined
  !     when the marching stopped'. When we are done, the final
  !     desired grid will be in grid3D with the corect dimensions
  !
  use hypData
  use hypInput
  implicit none
  
  ! Input
  integer(kind=intType), intent(in) :: nx

  ! Working
  integer(kind=intType) :: i, l, bin
  real(kind=realType) :: s(Nlayers), sAvg(nLayers), lAvg, R, sp1, factor, sp

  ! The first step is to compute a averaged global parameterization of
  ! the k-direction

  savg = zero
  lAvg = zero
  do i=1,nx
     s(1) = zero
     do l=2,nLayers
        s(l) = s(l-1) + dist(pGrid3d(:, i, l), pGrid3d(:, i, l-1))
     end do

     ! Add the length of this curve:
     lAvg = lAvg + s(nLayers)
     
     ! Normalize
     s = s / s(nLayers)
     ! Add to running total
     savg = savg + s
  end do

  ! Final s average
  savg = savg/nx

  ! Final averge length
  lavg = lAvg/nx

  ! Compute 'r' for the exponential distribution
  sp1 = s0/lavg
  R = -log((N-1)*sp1)/(N-2)

  ! Allocate grid3d

  ! Master loop over final levels:
  allocate(grid3D(3, nx, N))
  grid3D = zero

  ! Special case for l=1 --- Always the boundary:
  grid3d(:, :, 1) = pGrid3d(:, :, 1)

  do l=2,N

     ! Generate the position where we need to interpolate
     sp = sp1*(l-1)*exp(R*(l-2))

     ! Now what we want to find is the "bin" in l where this spacing
     ! falls. What we want is something like l=2.75. This means we
     ! interpolate three quarters of the way between l=2 and l=3. 

     call FindBin(sp, savg, Nlayers, bin)

     ! Bin is index of the lower level, bin+1 is index of the upper level

     ! Now get the factor
     factor = (sp-savg(bin))/(savg(bin+1)-savg(bin))

     ! Now in vector form:
     grid3d(:,:,l) = (one-factor)*pGrid3D(:,:,bin) + factor*pGrid3D(:,:,bin+1)
  end do
  
contains

  function dist(p1, p2)
    
    real(kind=realType) :: p1(3), p2(3)
    real(kind=realType) :: dist
    dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)
    
  end function dist

end subroutine reparameterize

subroutine FindBin(ss, s, nn, bin)

  ! Basic bisection search: search 'ss' in 's' (of length 'nn') and
  ! return 'bin'

  use precision
  implicit none

  ! Input/Output
  real(kind=realType), intent(in) :: ss, s(nn)
  integer(kind=intType), intent(in) :: nn
  integer(kind=intType), intent(out) :: bin

  ! Working
  integer(kind=intType) :: low, high

  ! Explicit check of upper bound:
  if (ss >= s(nn)) then
     bin = nn-1
     return
  end if

  if (ss < s(1)) then
     bin = 1
     return
  end if

  ! Do binary search
  low = 1
  high = nn
  bin = (low+high)/2 ! Integer division
     
  do while(ss < s(bin) .or. ss >= s(bin+1))
     
     if (ss < s(bin)) then
        high = bin
     else
        low = bin
     end if
     
     bin = (low + high)/2
  end do
  
end subroutine FindBin

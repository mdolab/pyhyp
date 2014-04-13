subroutine run3DElliptic
  ! run3DElliptic is the main python interface for generatnig 3D elliptic meshes. 
  !
  !
  ! Notes
  ! -----
  ! The init3d routine must have already been called with the
  ! addtional setup information. 

  use hypData
  use hypInput
  use panel
  implicit none

  ! Working parameters
  real(kind=realType), dimension(:, :), allocatable :: A
  real(kind=realType), dimension(:), allocatable :: rhs
  integer(kind=intType), dimension(:), allocatable :: IPIV
  integer(kind=intType) :: i, j, k, info, ierr
  real(kind=realType) :: phi, V(3)
  real(kind=realType), dimension(:), pointer :: vtmp, phitmp
 
  ! allocate(A(ncells, ncells), rhs(ncells), IPIV(ncells))
  ! print *,'forming', Ncells
  ! do i=1, nCells
  !    do j=1, nCells
  !       !call panelInfluence(i, panels(j)%center, phi, V)
  !       call panelInfluence(i, XSurf(:, j), phi, V)
  !       A(j, i) = phi
  !    end do
  ! end do
  ! rhs(:) = one
  ! print *,'solving....'
  ! call DGESV(ncells, 1, A, ncells, IPIV, rhs, ncells, info)

  ! do i=1,ncells
  !    panels(i)%strength = rhs(i)
  ! end do

  ! print *, 'rhs',rhs
  ! deallocate(A, rhs, ipiv)

  ! ! Lets fill in phi and Vx into some of the metric data for
  ! ! debugging
  ! do k=1,N ! This is the march layer
  !    call VecGetArrayF90(X(k), xxtmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)

  !    call VecGetArrayF90(X_ksi(k), vtmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)

  !    call VecGetArrayF90(Vhist(k), phitmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)
  !    print *,'evaluating on layer..', k
  !    do j=1,nx

  !       call slowEval(xxtmp(3*j-2:3*j), phi, V)
  !       vtmp(3*j-2:3*j) = V
  !       phitmp(3*j) = phi
  !       !phitmp(3*j) = panels(j)%strength
  !    end do

  !    call VecRestoreArrayF90(X(k), xxtmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)

  !    call VecRestoreArrayF90(X_ksi(k), vtmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)

  !    call VecRestoreArrayF90(Vhist(k), phitmp, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)

  ! end do

end subroutine run3DElliptic

! subroutine slowEval(pt, phi, V)
  
!   ! slowEval uses the N^2 algorithm to determine the potential and
!   ! velocity at a given pt. 
!   !  
!   ! Parameters
!   ! ----------
!   ! pt : real, dimension(3)
!   !   The point to use for evaluation
!   !
!   ! Returns
!   ! -------
!   ! phi : real
!   !   The computed potential
!   ! V : real size(3)
!   !   The computed velocity at point pt
!   use hypData
!   use panel
!   implicit none
  
!   ! Input Parameters
!   real(kind=realType), intent(in), dimension(3) :: pt
  
!   ! Output Parameters
!   real(kind=realType), intent(out) :: phi
!   real(kind=realType), dimension(3), intent(out) :: V

!   ! Working parameters
!   integer(kind=intType) :: i
!   real(kind=realType) :: pp, vv(3)

!   ! V(:) = zero
!   ! phi = zero
!   ! do i=1, ncells
!   !    call panelInfluence(i, pt, pp, vv)
!   !    V = V + panels(i)%strength*vv
!   !    phi = phi + panels(i)%strength*pp
!   ! end do

! end subroutine slowEval

subroutine panelInfluence(i, pt, phi, V)

  ! panelInfluence generates the influence of panel 'i' on point
  !  'pt'.
  !  
  ! Parameters
  ! ----------
  ! i : int
  !   The index of the panel to use. 
  ! pt : real size (3)
  !   Three dimensional point to get the influence at
  !
  ! Returns
  ! -------
  ! phi : real
  !   The computed potential
  ! V : real size(3)
  !   The computed velocity at point pt
  !
  ! Notes
  ! -----
  ! This routines uses the gloablly stored surface mesh and the
  ! connectivity array to determine the nodes of panel 'i'

  use hypData
  use panel
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: i
  real(kind=realType), intent(in) :: pt(3)

  ! Output Parameters
  real(kind=realtype) :: phi, V(3)

  ! Working Parameters
  real(kind=realType), dimension(3, NNodesMax+1) :: pts
  real(kind=realType), dimension(nNodesMax+1) :: r
  real(kind=realType), dimension(3) :: c, d, edge
  real(kind=realType) :: dist, dTheta, d12, c12, s12, s12_1
  real(kind=realType) :: r1, r2, s12_2, R12, Q12, J12
  real(kind=realType) ::  tmp
  integer(kind=intType) :: iNode, ii
  type(panelType), pointer :: p

  ! Zero out the panel influence and velocity
  phi = zero
  V(:) = zero

  ! Poiner to panel to make code easier to read
  p => MGP(1)%panels(i)

  ! vector from the panel center to the point of interest:
  d = pt - p%center
 
  ! Cartiesian distance
  dist = sqrt(d(1)**2 + d(2)**2 + d(3)**2)

  ! ! Check if we can do the farfield tolerance:
  ! if (dist < farFieldTol) then
  ! if (.true.) then
  !    ! Do something else:
  !    phi = 2*pi/dist
  !    V(:) = zero
  !    return
  ! end if

  ! Transform the vector of point into the local frame as well as the
  ! 4 corner points
  d = matmul(p%C, d)
 
  do ii=1,p%N
     pts(:, ii) = matmul(p%C, p%x(:, p%N+1-ii) - p%center)
     r(ii) = sqrt((d(1) - pts(1, ii))**2 + (d(2) - pts(2, ii))**2 + d(3)**2)
  end do
  pts(:, p%N+1) = pts(:, 1)
  r(p%N+1) = r(1)
  dTheta = 2*pi
 
  ! Loop over each of the 'N' edges
  do iNode=1, p%N
     edge = pts(:, iNode + 1) - pts(:, iNode)
     d12 = sqrt(edge(1)**2 + edge(2)**2)

     ! Hess equation 4.4.10
     C12 = (pts(1, iNode+1) - pts(1, iNode))/d12
     S12 = (pts(2, iNode+1) - pts(2, iNode))/d12

     ! Hess equation 4.4.12
     s12_1 = (pts(1, iNode)   - d(1))*C12 + (pts(2, iNode  ) - d(2))*S12
     s12_2 = (pts(1, iNode+1) - d(1))*C12 + (pts(2, iNode+1) - d(2))*S12

     ! Hess equation 4.4.13
     R12 = (d(1) - pts(1, iNode))*S12 - (d(2) - pts(2, iNode))*C12
  
     if (R12 < zero) then
        dTheta = zero
     end if

     ! Hess equation 4.4.14
     r1 = r(iNode)
     r2 = r(iNode+1)

     ! Hess equation 4.4.15
     Q12 = log((r2 + s12_2)/(r1 + s12_1))

     ! Hess equation 4.4.16
     J12 = atan2( (R12*abs(d(3))*(r1*s12_2 - r2*s12_1)) , &
          (r1*r2*R12**2 + d(3)**2*s12_1*s12_2))
     ! Hess equation 4.4.17 & 4.4.18
     phi = phi + R12*Q12 + abs(d(3))*J12

     ! Hess equation 4.4.19
     V(1) = V(1) - S12*Q12
     V(2) = V(2) + C12*Q12
     V(3) = V(3) - J12
  end do

  ! dTheta correction to Hess Equation 4.4.18 and 4.4.19
  phi = phi - abs(d(3))*dTheta

  V(3) = dTheta + V(3)
  if (d(3) < zero) then
     V(3) = -V(3)
  end if

  ! Transform the velocity back to the global frame
  V = matmul(p%CT, V)

end subroutine panelInfluence


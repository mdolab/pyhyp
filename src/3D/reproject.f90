subroutine allocateSurfaces(nSurf)

  use bspline
  implicit none

  integer(kind=intType), intent(in) :: nSurf
  allocate(surfs(nsurf))

end subroutine allocateSurfaces

subroutine setSurface(iSurf, ku, kv, tu, tv, coef, nCtlu, nCtlv)
  use bspline
  implicit none

  ! Set the supplied surface data into the specified slot
  integer(kind=intType), intent(in) :: ku, kv, nCtlu, nCtlv, iSurf
  real(kind=realType), intent(in), dimension(nctlu+ku) :: tu
  real(kind=realType), intent(in), dimension(nctlv+kv) :: tv
  real(kind=realType), intent(in), dimension(3, nctlv, nctlu) :: coef

  allocate(surfs(iSurf)%tu(nctlu+ku))
  allocate(surfs(iSurf)%tv(nctlv+kv))
  allocate(surfs(iSurf)%coef(3, nctlv, nctlu))

  ! Now set the data
  surfs(iSurf)%ku = ku
  surfs(iSurf)%kv = kv
  surfs(iSurf)%tu = tu
  surfs(iSurf)%tv = tv
  surfs(iSurf)%coef = coef

end subroutine setSurface



! recursive subroutine getQuad(surf, uLow, uHigh, vLow, vHigh, tol)
!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract: getQuad is the main (recursive) subroutine for
!   !     adaptively "meshing" a b-spline surface patch to a given tolerance. 
!   !
!   !     Description of Arguments
!   !     Input
!   !     surf  - type(surface) B-spline surface type
!   !     uLow  - Real, Lower bound of u defining current parametric box
!   !     uHigh - Real, Upper bound of u
!   !     vLow  - Real, Lower v bound
!   !     vHigh - Real, Upper v bound
!   !     tol   - Real, The target tolerance along the edges in consistent units
!   use bspline
!   implicit none

!   ! Input/Output
!   type(surfType), intent(in) :: surf
!   real(kind=realType), intent(in) :: uLow, uHigh, vLow, vHigh, tol

!   ! Working:
!   real(kind=realType), dimension(3, 4) :: corners, approxMidPts
!   real(kind=realType), dimension(3, 4) :: midPts, d
!   real(kind=realType) :: uMId, vMid, tol2
!   integer(kind=intType) :: i
!   logical :: splitU, splitV

!   ! Step 1: Evaluate the 4 corners:
!   call evalSurface(surf, (/uLow, uHigh, uHigh, uLow/), &
!        (/vLow, vLow, vHigh, vHigh/), 4, corners) 

!   ! Step 2: Evaluate the approximate midPts from the previously computed corners
!   do i=1, 4
!      approxMidPts(:, i) = 0.5_8*(corners(:, i) + corners(:, mod(i, 4)+1))
!   end do

!   ! Step 3: Evaluate the mid-pts of the edges
!   uMid = 0.5_8*(uLow + uHigh)
!   vMid = 0.5_8*(vLow + vHigh)
!   call evalSurface(surf, (/uMid, uHigh, uMid, uLow/), &
!        (/vLow, vMid, vHigh, vMid/), 4, midPts)

!   ! Step 4: Check midPts against approxMidPts for being less than the
!   ! specified tolerance
!   splitU = .False.
!   splitV = .False.
!   tol2 = tol**2
!   d = midPts - approxMidPts
!   if ((d(1,1)**2 + d(2,1)**2 + d(3,1)**2) > tol2 .or. &
!        (d(1,3)**2 + d(2,3)**2 + d(3,3)**2) > tol2) then
!      splitU = .True.
!   end if

!   if ((d(1,2)**2 + d(2,2)**2 + d(3,2)**2) > tol2 .or. &
!        (d(1,4)**2 + d(2,4)**2 + d(3,4)**2) > tol2) then
!      splitV = .True.
!   end if

!   ! Step 5: Determine if quad is OK or we have to recursively call
!   ! ourself with smaller ranges:

!   if (.not. splitU .and. .not. splitV) then 
!      ! This quad is ok...add it to the list:
!      nElements = nElements + 1
!      do i=1,4
!         nodes(:, i + nNodes) = corners(:, i)
!         elements(i, nElements) = nNodes*4 + i
!      end do
!      nNodes = nNodes + 4
!   else
!      ! Otherwise we have to split 4 ways or two ways:
!      if (splitU .and. splitV) then 
!         call getQuad(surf, uLow, uMid, vLow, vMid, tol)
!         call getQuad(surf, uMid, uHigh, vLow, vMid, tol)
!         call getQuad(surf, uLow, uMid , vMid, vHigh, tol)
!         call getQuad(surf, uMid, uHigh, vMid, vHigh, tol)
!      else if(splitU) then
!         call getQuad(surf, uLow, uMid , vLow, vHigh, tol)
!         call getQuad(surf, uMid, uHigh, vLow, vHigh, tol)
!      else if (splitV) then
!         call getQuad(surf, uLow, uHigh, vLow, vMid, tol)
!         call getQuad(surf, uLow, uHigh, vMid, vHigh, tol)
!      end if
!   end if
! end subroutine getQuad

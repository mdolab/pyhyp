subroutine computeStretch(L)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: computeStretch determines the global stretch value,
  !     deltaS, depending on the grid level L. In the future
  !     this function may be more complex or call a user-supplied function.
  !    
  !     Parameters
  !     ----------
  !     L : integer
  !         The marching level to compute stretch for
  !
  !     Returns
  !     -------
  !     deltaS : real
  !         Marching increment set in hypData
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

subroutine calcGridRatio(N, s0, S, ratio)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: calcGridRatio() calculates the exponential
  !     distribution Turns out we need to solve a transendental
  !     equation here. We will do this with a secant iteration.
  !
  !     Parameters
  !     ----------
  !     N : integer
  !         The number of nodes in sequence
  !     s0 : real
  !         The initial grid spacing
  !     S : real
  !         The total integrated length
  !     
  !     Returns
  !     -------
  !     ratio : real
  !         The computed grid ratio that satifies, N, s0, and S.

  use precision

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in) :: s0, S

  ! Output Parameters
  real(kind=realType), intent(out) :: ratio

  ! Working Parameters
  integer(kind=intType) :: i, M
  real(kind=realType) ::  r, a,b, c, f, fa, fb

  ! function 'f' is S - s0*(1-r^n)/(1-r) where S is total length, s0 is
  ! initial ratio and r is the grid ratio. 

  M = N-1

  ! Do a bisection search
  ! Max and min bounds...root must be in here...
  a = one + 1e-8
  b = four

  fa = S - s0*(1-a**M)/(1-a)
  fb = S - s0*(1-b**M)/(1-b)
  do i=1, 100
     c = half*(a + b)
     f = S - s0*(1-c**M)/(1-c)
     if (abs(f) < 1e-6) then ! Converged
        exit
     end if
     
     if (f * fa > 0) then 
        a = c
     else
        b = c
     end if
  end do

  ! Finally set the ratio variable to r
  ratio = c

end subroutine calcGridRatio

subroutine three_by_three_inverse(Jac, Jinv)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: return the inverse of the 3x3 matrix Jac in Jinv
  !    
  !     Parameters
  !     ----------
  !     Jac : array size (3,3)
  !         Matrix to invert
  !
  !     Returns
  !     -------
  !     Jinv : array size (3,3)
  !         Inverse of Jac
  !
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

subroutine cross_prod(a,b,c)

  ! Written by Ney Secco, 2015
  ! This function computes the cross-product between a and b

  use precision

  ! Inputs
  real(kind=realType), dimension(3), intent(in) :: a,b

  ! Outputs
  real(kind=realType), dimension(3), intent(out) :: c

  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)

end subroutine cross_prod

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
  do i=2,140
     write(*,"(a)",advance="no") "-"
  enddo
  print "(1x)"
  write(*,"(a)",advance="no") "# Grid  |     CPU    | Sub  | KSP  |     Sl     |"
  write(*,"(a)",advance="no") " Grid       | Grid       |     Min    |   deltaS   |"
  write(*,"(a)",advance="no") "    cMax    |    min R   |    max     | "
  print "(1x)"
  write(*,"(a)",advance="no") "# Level |     Time   | Iter | Its  |            |"
  write(*,"(a)",advance="no") " Sensor Max | Sensor Min |  Quality   |            |"
  write(*,"(a)",advance="no") "            |            |  KStretch  | "
  print "(1x)"

  ! Write the Last line of '-'
  write(*,"(a)",advance="no") "#"
  do i=2,140
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
  use communication
  use hypData

  implicit none

  ! Working variables
  real(kind=realType) :: gridSensorMaxRed, gridSensorMinRed
  integer(kind=intType) :: ierr

  ! Reduce all mins/maxes
  CALL MPI_REDUCE(gridSensorMax, gridSensorMaxRed, 1, MPI_DOUBLE, MPI_MAX, 0, hyp_comm_world, ierr)
  CALL MPI_REDUCE(gridSensorMin, gridSensorMinRed, 1, MPI_DOUBLE, MPI_MIN, 0, hyp_comm_world, ierr)

  ! Only root proc prints stuff
  if (myid == 0) then

     ! Iteration Count
     write(*,"(i8,1x)",advance="no") marchIter

     ! CPU time
     write(*,"(e12.5,1x)",advance="no") dble(mpi_wtime()) - timeStart

     ! Sub Iteration
     write(*,"(i6,1x)",advance="no") nSubIter

     ! KSP ITerations
     write(*,"(i6,1x)",advance="no") kspIts

     ! Modified Sl factor from Steger and Chan
     write(*,"(e12.5,1x)",advance="no") Sl

     ! maximum Grid convergence sensor
     write(*,"(e12.5,1x)",advance="no") gridSensorMaxRed

     ! minimum Grid convergence sensor
     write(*,"(e12.5,1x)",advance="no") gridSensorMinRed

     ! minimum grid quality for new layer
     write(*,"(e12.5,1x)",advance="no") minQuality ! This was reduced already

     ! marching step size for this iteration
     write(*,"(e12.5,1x)",advance="no") deltaS

     ! marching step size for this iteration
     write(*,"(e12.5,1x)",advance="no") cRatio

     ! approximate minimim distance to body
     minR = zero
     write(*,"(e12.5,1x)",advance="no") scaleDist/radius0

     ! maximum stretch ratio in K-direction
     write(*,"(e12.5,1x)",advance="no") maxKStretch

     ! Jump to next line
     print "(1x)"

  end if

end subroutine writeIteration

subroutine computeR(xVec, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     computeR0 determines the minimum radius of a sphere that will
  !     fully enclose the grid level defined by xVec
  !
  !     Parameters
  !     ----------
  !     xVec : Petsc vector 
  !        Vector of nodes
  !
  !     Returns
  !     -------
  !     R : real
  !        Minimum radius surrounding xVec

  use communication
  use hypData, only : xAvg, xx, nXGlobal, nX
  implicit none

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input parameters
  Vec, intent(in) :: xVec
  
  ! Ouput parameters
  real(kind=realType), intent(out) :: R

  ! Working
  integer(kind=intType) :: i, ierr
  real(kind=realType) :: dist, XavgLocal(3), rLocal

  ! First find the average of X
  XavgLocal = zero

  call VecGetArrayF90(xVec, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do i=1,nX
     XavgLocal = XavgLocal + xx(3*i-2:3*i)
  end do

  ! All reduce across all procs:
  call mpi_AllReduce(XavgLocal, Xavg, 3, MPI_DOUBLE, MPI_SUM, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And Now divide by total number of points nGloabl
  Xavg = Xavg / nXGlobal

  ! Now find the maximum distance of X from Xavg
  rLocal = zero
  do i=1,nx
     rLocal = max(rLocal, dist(Xavg, xx(3*i-2:3*i)))
  end do

  ! All reduce with max
  call mpi_AllReduce(rLocal, R, 1, MPI_DOUBLE, MPI_MAX, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Always remember to restore the array
  call VecRestoreArrayF90(xVec, xx, ierr)

end subroutine computeR

subroutine computeMinR(xVec, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: computeMinR determines the closest point of xVec to
  !     Xavg
  !
  !     Parameters
  !     ----------
  !     xVec : Petsc vector 
  !        Vector of nodes
  !
  !     Returns
  !     -------
  !     R : real
  !        Minimum distance from xVec to xAvg

  use communication
  use hypData, only : nx, xAvg, xx
  implicit none

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif
  ! Input parameters
  Vec, intent(in) :: xVec

  ! Ouput parameters
  real(kind=realType), intent(out) :: R

  ! Working
  integer(kind=intType) :: i, ierr
  real(kind=realType) :: dist, Rlocal

  call VecGetArrayF90(xVec, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  Rlocal = huge(Rlocal)
  ! Now find the minimum distance of X from Xavg
  do i=1,nX
     Rlocal = min(Rlocal, dist(Xavg, xx(3*i-2:3*i)))
  end do

  ! All reduce with min
  call mpi_AllReduce(rLocal, R, 1, MPI_DOUBLE, MPI_MIN, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Always remember to restore the array
  call VecRestoreArrayF90(xVec, xx, ierr)

end subroutine computeMinR

subroutine computeQualityLayer
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Compute the minimum quality measure for the cells defined by
  !     two grid levels
  !
  use communication
  use hypInput
  use hypData
  implicit none

  ! Working
  integer(kind=intType) :: i, j, iPatch, ierr
  real(kind=realType) :: points(3,8), Q
  real(Kind=realType) :: minQuality_local
  ! Since we are dealing with volumes, we need to loop over faces:
  minQuality_local = one

  call VecGhostUpdateBegin(X(marchIter), INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecGhostUpdateEnd(X(marchIter), INSERT_VALUES,SCATTER_FORWARD, ierr)

  call VecGhostGetLocalForm(X(marchIter), XL_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGhostGetLocalForm(X(marchIter-1), XLm1_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(XL_local, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(XLm1_local, xxm1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do i=1,nLocalFace
  
     ! Assemble the 8 points we need for the volume - Note
     ! coodinate coordinate ordering!
     points(:,1) = xxm1(3*conn(1, i)-2: 3*conn(1, i))
     points(:,2) = xxm1(3*conn(2, i)-2: 3*conn(2, i))
     points(:,3) = xxm1(3*conn(4, i)-2: 3*conn(4, i))
     points(:,4) = xxm1(3*conn(3, i)-2: 3*conn(3, i))

     points(:,5) = xx(3*conn(1, i)-2: 3*conn(1, i))
     points(:,6) = xx(3*conn(2, i)-2: 3*conn(2, i))
     points(:,7) = xx(3*conn(4, i)-2: 3*conn(4, i))
     points(:,8) = xx(3*conn(3, i)-2: 3*conn(3, i))
                   
     call quality_hexa(points, Q)
     minQuality_local = min(minQuality_local, Q)
              
  end do

  call VecRestoreArrayF90(XL_local, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(XLm1_local, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGhostRestoreLocalForm(X(marchIter), XL_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGhostRestoreLocalForm(X(marchIter-1), XLm1_local, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Need to reduce the min value to the root proc
  call mpi_Reduce(minQuality_local, minQuality, 1, MPI_DOUBLE, MPI_MIN, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine computeQualityLayer

subroutine quality_hexa(points, quality)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: quality_hexa determine the 2x2x2 quality measure for
  !     a hexa defined by 8 nodes. 
  !    
  !     Parameters
  !     ----------
  !     points : real array size(3, 8)
  !         Nodes defining the hexa in coordinate ording. 
  !
  !     Returns
  !     -------
  !     Quality : real 
  !         The computed quality measure
  use precision
  implicit none

  ! Input/Output
  real(kind=realType), intent(in) :: points(3, 8)
  real(kind=realType), intent(out) :: quality

  ! Working
  integer(kind=intType), parameter :: Ng=2
  real(kind=realType) :: det(Ng*Ng*Ng)
  integer(kind=intType):: i, j, k, l, m, nn, counter
  real(kind=realType) :: r, s, t
  real(kind=realType) :: largest, d_max, corners(Ng)
  real(kind=realType) ::  Xd(9), N(8), Na(8), Nb(8), Nc(8), Ua(9)
  real(kind=realType) :: JJac(9), h, pt(3)
  real(kind=realType) :: shp(2, 2), dshp(2, 2)
  ! Relative Determinant: the ratio of the smallest determinat of the
  ! jacobian matrix divided by the largest determinate of the jacobian
  ! matrix using a 3x3x3 stecil

  corners(1) = -1
  corners(2) = 1
  call lagrangeSF(shp(:, 1), dshp(:, 1), corners(1), 2)
  call lagrangeSF(shp(:, 2), dshp(:, 2), corners(2), 2)

  ! 2x2x2 Stencil
  do nn=1, Ng ! Gauss in r
     do m=1, Ng ! Gauss in s
        do l=1, Ng ! Gauss in t

          Xd = zero
           do k=1, 2
              do j=1, 2
                 do i=1, 2
                    counter = (k-1)*4 + (j-1)*2 + i
                    Na(counter) = dshp(i, l)*shp(j, m)*shp(k, nn)
                    Xd(1) = Xd(1) + points(1, counter)*Na(counter)
                    Xd(2) = Xd(2) + points(2, counter)*Na(counter)
                    Xd(3) = Xd(3) + points(3, counter)*Na(counter)
                    
                    Nb(counter) = shp(i, l)*dshp(j, m)*shp(k, nn)
                    Xd(4) = Xd(4) + points(1, counter)*Nb(counter)
                    Xd(5) = Xd(5) + points(2, counter)*Nb(counter)
                    Xd(6) = Xd(6) + points(3, counter)*Nb(counter)
                    
                    Nc(counter) = shp(i, l)*shp(j, m)*dshp(k, nn)
                    Xd(7) = Xd(7) + points(1, counter)*Nc(counter)
                    Xd(8) = Xd(8) + points(2, counter)*Nc(counter)
                    Xd(9) = Xd(9) + points(3, counter)*Nc(counter)
                    
                 end do
              end do
           end do

           h = (Xd(9)*(Xd(1)*Xd(5) - Xd(4)*Xd(2)) &
                - Xd(8)*(Xd(1)*Xd(6) - Xd(4)*Xd(3)) &
                + Xd(7)*(Xd(2)*Xd(6) - Xd(3)*Xd(5)))

           det((nn-1)*Ng*Ng + (m-1)*Ng + l) = h
        end do
     end do
  end do

  ! We don't use maxval since we can't cs it
  largest = maxval(abs(det))
  d_max = maxval(largest - det)

  quality = 1-d_max/largest

end subroutine quality_hexa

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

function dist(p1, p2)
  use precision

  real(kind=realType) :: p1(3), p2(3)
  real(kind=realType) :: dist
  dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)

end function dist


subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

  ! Given a list of N points (pts) in three space, with possible
  ! duplicates, (to within tol) return a list of the nUnique
  ! uniquePoints of points and a link array of length N, that points
  ! into the unique list
  use precision
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(:, :) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(:,:) :: uniquePts
  integer(kind=intType), intent(out), dimension(:) :: link
  integer(kind=intType), intent(out) :: nUnique

  ! Working Parameters
  real(kind=realType), allocatable, dimension(:) :: dists, tmp
  integer(kind=intType), allocatable, dimension(:) :: ind
  integer(kind=intType) :: i, j, nTmp, link_counter, ii, nSubUnique
  logical cont, cont2
  integer(kind=intType), dimension(:), allocatable :: tmpInd, subLInk
  real(kind=realType), dimension(:, :), allocatable :: subPts, subUniquePts
  real(kind=realType), dimension(3) :: Xavg
  integer(kind=intType) :: maxSubUnique
  interface
     subroutine pointReduceBruteForce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none
       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduceBruteForce

  end interface
  maxSubUnique = 10
  ! Allocate dists, and the ind pointer
  allocate(dists(N), tmp(N), ind(N), tmpInd(maxSubUnique), subLink(maxSubUnique), &
       subPts(3, maxSubUnique), subUniquePts(3, maxSubUnique))

  Xavg(1) = sum(pts(:, 1))/N
  Xavg(2) = sum(pts(:, 2))/N
  Xavg(3) = sum(pts(:, 3))/N

  ! Compute distances of all points from the origin
  do i=1, N
     dists(i) = sqrt((pts(1,i)-Xavg(1))**2 + (pts(2,i)-Xavg(2))**2 + (pts(3, i)-Xavg(3))**2)
     tmp(i) = dists(i)
     ind(i) = i
  end do

  ! Do an argsort on the distances
  call ArgQsort(tmp, N, ind)

  i = 1
  cont = .True.
  link_counter = 0
  nUnique = 0

  do while(cont)
     cont2 = .True.
     j = i
     nTmp = 0
     do while(cont2)
        if (abs(dists(ind(i))-dists(ind(j))) < tol) then
           nTmp = nTmp + 1
           j = j + 1
           if (j == N+1) then ! Overrun check
              cont2 = .False.
           end if
        else
           cont2 = .False.
        end if
     end do

     ! Not enough space...deallocate and reallocate
     if (ntmp > maxSubUnique) then
        deallocate(tmpInd, subLink, subPts, subUniquePts)
        maxSubUnique = nTmp
        allocate(tmpInd(maxSubUnique), subLink(maxSubUnique), &
             subPts(3, maxSubUnique), subUniquePts(3, maxSubUnique))
     end if

     ! Copy the points that have the same distance into subPts. Note
     ! these may NOT be the same, since two points can have the same
     ! distance, but not be co-incident (ie (1,0,0), (-1,0,0))
     do ii=1,nTmp
        tmpInd(ii) = ind(j - nTmp + ii - 1)
        subPts(:, ii) = pts(:, tmpInd(ii))
     end do

     ! Brute Force Search them 
     call pointReduceBruteForce(subPts, nTmp, tol, subUniquePts, subLink, nSubUnique)

     do ii=1,nSubUnique
        nUnique = nUnique + 1
        uniquePts(:, nUnique) = subUniquePts(:,ii)
     end do

     do ii=1,nTmp
        link(tmpInd(ii)) = subLink(ii) + link_counter
     end do

     link_counter = link_counter +  maxval(subLink) 

     i = j - 1 + 1
     if (i == N+1) then
        cont = .False.
     end if
  end do
  deallocate(dists, tmp, ind, tmpInd, subLink, subPts, subUniquePts)

end subroutine pointReduce

subroutine pointReduceBruteForce(pts, N, tol, uniquePts, link, nUnique)

  ! Given a list of N points (pts) in three space, with possible
  ! duplicates, (to within tol) return a list of the nUnqiue
  ! uniquePoints of points and a link array of length N, that points
  ! into the unique list

  use precision 
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(:, :) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(:,:) :: uniquePts
  integer(kind=intType), intent(out), dimension(:) :: link
  integer(kind=intType), intent(out) :: nUnique

  ! Working parameters
  integer(kind=intType) :: i, j
  real(kind=realType) :: dist
  logical :: found_it

  ! First point is *always* unique
  uniquePts(:, 1) = pts(:, 1)
  link(:) = 0
  link(1) = 1
  nUnique = 1

  do i=2,N
     found_it = .False.
     uniqueLoop: do j=1,nUnique
        dist = sqrt((pts(1, i)-uniquePts(1, j))**2 + &
             (pts(2, i) - uniquePts(2, j))**2 + &
             (pts(3, i) - uniquePts(3, j))**2)
        if (dist < tol) then
           link(i) = j
           found_it = .True. 
           exit uniqueLoop
        end if
     end do uniqueLoop

     if (.not. found_it) then
        nUnique = nUnique + 1
        uniquePts(:, nUnique) = pts(:, i)
        link(i) = j 
     end if
  end do

end subroutine pointReduceBruteForce

recursive subroutine ArgQSort(A, nA, ind)

  ! Do an ArgQuickSort. Adapted from
  ! http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#FPr. Modified
  ! such that array 'A' is unchanged and the index 'ind' is initialzed
  ! inside the algorithm

  use precision
  implicit none

  ! DUMMY ARGUMENTS
  integer(kind=intType), intent(in) :: nA
  real(kind=realType), dimension(nA), intent(inout) :: A
  integer(kind=intType), dimension(nA), intent(inout) :: ind

  ! LOCAL VARIABLES
  integer(kind=intType) :: left, right, itemp, i
  real(kind=realType) :: random, pivot, temp
  integer(kind=intType) :: marker

  if (nA > 1) then

     call random_number(random)
     i = int(random*real(nA-1))+1
     pivot = A(i)    ! random pivot (not best performance, but avoids worst-case)
     left = 0
     right = nA + 1

     do while (left < right)
        right = right - 1
        do while (A(right) > pivot)
           right = right - 1
        end do
        left = left + 1
        do while (A(left) < pivot)
           left = left + 1
        end do
        if (left < right) then
           ! Swap value 
           temp = A(left)
           A(left) = A(right)
           A(right) = temp
           ! ! And swap index
           iTemp = ind(left)
           ind(left) = ind(right)
           ind(right) = itemp

        end if
     end do

     if (left == right) then
        marker = left + 1
     else
        marker = left
     end if

     call argQSort(A(:marker-1), marker-1, ind(:marker-1))
     call argQSort(A(marker:), nA-marker+1, ind(marker:))
  end if
end subroutine ArgQSort

subroutine getSurfaceCoordinates(coords, n)

  ! Return the surface coordiantes

  use hypData
  implicit none

  integer(kind=intType), intent(in) :: n
  real(kind=realType), intent(inout), dimension(n) :: coords
  integer(kind=intType) :: ierr

  call VecGetArrayF90(X(n), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Just copy
  coords = xx
  call VecRestoreArrayF90(X(n), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine getSurfaceCoordinates

subroutine setSurfaceCoordinates(coords, n)

  ! Return the surface coordiantes

  use hypData
  implicit none

  integer(kind=intType), intent(in) :: n
  real(kind=realType), intent(in), dimension(n) :: coords
  integer(kind=intType) :: ierr

  call VecGetArrayF90(X(n), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Just copy
  xx = coords

  call VecRestoreArrayF90(X(n), xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setSurfaceCoordinates

subroutine assignNode2NodeConn(nodeConn, node1, node2)

  ! This subroutine assigns connection from node 1 to node 2
  ! in the nodeConn matrix. If this connection is already assigned
  ! then we should flag this node as this indicates flipped normals

  use precision
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: node1, node2

  ! Input/Output variables
  integer(kind=intType), dimension(:,:), intent(inout) :: nodeConn

  ! Working variables
  integer(kind=intType) :: i, numConn

  ! BEGIN EXECUTION

  ! Get the total number of connections allowed
  numConn = size(nodeConn,1)

  do i=1,numConn

     ! Check if node1 to node2 connection already exists
     if (nodeConn(i, node1) .eq. node2) then
        ! Flag this node by setting all connections to -1
        nodeConn(:,node1) = -1
        ! Get out of this loop
        exit

     else if (nodeConn(i, node1) .eq. 0) then
        ! If we reach 0, then we already checked all connections
        ! and we can assign the new one without problems
        nodeConn(i, node1) = node2
        ! Get out of this loop
        exit

     end if

  end do

end subroutine assignNode2NodeConn

subroutine extrapolatePoint(x1,x2,x3)

  ! This subroutine extrapolates the line connecting points x1 and x2
  ! to find x3. x2 will be the middle point between x1 and x3

  use precision
  implicit none

  ! Input variables
  real(kind=realType), intent(in) :: x1(3), x2(3)

  ! Output variables
  real(kind=realType), intent(out) :: x3(3)

  ! BEGIN EXECUTION

  x3 = 2*x2 - x1

end subroutine extrapolatePoint

subroutine cellArea(x1,x2,x3,x4,S)

  ! This subroutine computes the cell area defined by the nodes x1, x2, x3, x4.
  !
  !   ========
  !  |x1    x2|
  !  |        |
  !  |x4    x3|
  !   ========


  use precision
  implicit none

  ! Input variables
  real(kind=realType), intent(in) :: x1(3), x2(3), x3(3), x4(3)

  ! Output variables
  real(kind=realType), intent(out) :: S

  ! Working variables
  real(kind=realType) :: v1(3), v2(3), v3(3)

  ! BEGIN EXECUTION

  ! Compute diagonal vectors
  v1(:) = x3 - x1
  v2(:) = x4 - x2

  ! Cross Product
  call cross_prod(v1,v2,v3)
  
  ! Compute area
  S = half*sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))

end subroutine cellArea

subroutine extrapolateArea(x1,x2,x3,x4,side,areaRatio)

  ! This subroutine computes the ratio of the extrapolated area to the
  ! initial cell area defined by the nodes x1, x2, x3, x4.
  !
  !   ==========================
  !  |        |        |        |
  !  |   5    |   2    |   6    |
  !  |        |        |        |
  !   ==========================
  !  |        |x1    x2|        |
  !  |   1    |   0    |   3    |
  !  |        |x4    x3|        |
  !   ==========================
  !  |        |        |        |
  !  |   8    |   4    |   7    |
  !  |        |        |        |
  !   ==========================
  !
  !  side = 1 computes S(1+0)/S(0)
  !  side = 2 computes S(2+0)/S(0)
  !  side = 3 computes S(3+0)/S(0)
  !  side = 4 computes S(4+0)/S(0)
  !  side = 5 computes S(1+5+2+0)/S(0)
  !  side = 6 computes S(2+6+3+0)/S(0)
  !  side = 7 computes S(3+7+4+0)/S(0)
  !  side = 8 computes S(4+8+1+0)/S(0)

  use precision
  implicit none

  ! Input parameters
  real(kind=realType), intent(in) :: x1(3), x2(3), x3(3), x4(3)
  integer(kind=intType), intent(in) :: side

  ! Outputs parameters
  real, intent(out) :: areaRatio

  ! Working parameters
  real(kind=realType) :: xi(3), xj(3), xk(3), xl(3)
  real(kind=realType) :: xm1(3), xm2(3), xm(3)
  real(kind=realType) :: S0, Stot, deltaS

  ! BEGIN EXECUTION

  ! Compute original cell area
  call cellArea(x1,x2,x3,x4,S0)

  ! Initialize total area
  Stot = S0

  if (side .eq. 1) then
     ! Extrapolate points
     call extrapolatePoint(x2,x1,xi)
     call extrapolatePoint(x3,x4,xj)
     ! Compute area of new cells and add to total
     call cellArea(xi,x1,x4,xj,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 2) then
     ! Extrapolate points
     call extrapolatePoint(x4,x1,xi)
     call extrapolatePoint(x3,x2,xj)
     ! Compute area of new cells and add to total
     call cellArea(xi,xj,x2,x1,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 3) then
     ! Extrapolate points
     call extrapolatePoint(x1,x2,xi)
     call extrapolatePoint(x4,x3,xj)
     ! Compute area of new cells and add to total
     call cellArea(x2,xi,xj,x3,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 4) then
     ! Extrapolate points
     call extrapolatePoint(x1,x4,xi)
     call extrapolatePoint(x2,x3,xj)
     ! Compute area of new cells and add to total
     call cellArea(x4,x3,xj,xi,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 5) then
     ! Extrapolate points
     call extrapolatePoint(x2,x1,xi)
     call extrapolatePoint(x3,x4,xj)
     call extrapolatePoint(x4,x1,xk)
     call extrapolatePoint(x3,x2,xl)
     ! Extapolate the diagonal node taking the average of two extrapolations
     call extrapolatePoint(xj,xi,xm1)
     call extrapolatePoint(xl,xk,xm2)     
     xm = half*(xm1 + xm2)
     ! Compute area of new cells and add to total
     call cellArea(xi,x1,x4,xj,deltaS)
     Stot = Stot + deltaS
     call cellArea(xm,xk,x1,xi,deltaS)
     Stot = Stot + deltaS
     call cellArea(xk,xl,x2,x1,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 6) then
     ! Extrapolate points
     call extrapolatePoint(x4,x1,xi)
     call extrapolatePoint(x3,x2,xj)
     call extrapolatePoint(x1,x2,xk)
     call extrapolatePoint(x4,x3,xl)
     ! Extapolate the diagonal node taking the average of two extrapolations
     call extrapolatePoint(xi,xj,xm1)
     call extrapolatePoint(xl,xk,xm2)     
     xm = half*(xm1 + xm2)
     ! Compute area of new cells and add to total
     call cellArea(xi,xj,x2,x1,deltaS)
     Stot = Stot + deltaS
     call cellArea(xj,xm,xk,x2,deltaS)
     Stot = Stot + deltaS
     call cellArea(x2,xk,xl,x3,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 7) then
     ! Extrapolate points
     call extrapolatePoint(x1,x2,xi)
     call extrapolatePoint(x4,x3,xj)
     call extrapolatePoint(x1,x4,xk)
     call extrapolatePoint(x2,x3,xl)
     ! Extapolate the diagonal node taking the average of two extrapolations
     call extrapolatePoint(xi,xj,xm1)
     call extrapolatePoint(xk,xl,xm2)     
     xm = half*(xm1 + xm2)
     ! Compute area of new cells and add to total
     call cellArea(x2,xi,xj,x3,deltaS)
     Stot = Stot + deltaS
     call cellArea(x3,xj,xm,xl,deltaS)
     Stot = Stot + deltaS
     call cellArea(x4,x3,xl,xk,deltaS)
     Stot = Stot + deltaS
  else if (side .eq. 8) then
     ! Extrapolate points
     call extrapolatePoint(x2,x1,xi)
     call extrapolatePoint(x3,x4,xj)
     call extrapolatePoint(x1,x4,xk)
     call extrapolatePoint(x2,x3,xl)
     ! Extapolate the diagonal node taking the average of two extrapolations
     call extrapolatePoint(xi,xj,xm1)
     call extrapolatePoint(xl,xk,xm2)     
     xm = half*(xm1 + xm2)
     ! Compute area of new cells and add to total
     call cellArea(xi,x1,x4,xj,deltaS)
     Stot = Stot + deltaS
     call cellArea(xj,x4,xk,xm,deltaS)
     Stot = Stot + deltaS
     call cellArea(x4,x3,xl,xk,deltaS)
     Stot = Stot + deltaS
  end if

  ! Finally compute the area ratio
  areaRatio = Stot/S0

end subroutine extrapolateArea

subroutine findKStretch(XL, XLm1, XLm2)
  !***DESCRIPTION
  !
  !     Written by Ney Secco
  !
  !     Abstract: findKStretch computes a distance ratio. Let zeta and eta be
  !     the computational coordinates along the surface, while L represents
  !     the marching direction (the current layer). This functions compute
  !     the following ratio:
  !     Kstretch =  distance between X(zeta, eta, L) and X(zeta, eta, L-1)
  !                --------------------------------------------------------
  !                distance between X(zeta, eta, L-1) and X(zeta, eta, L-2)
  !     This can be used as an indication of how much the grid is stretching in the
  !     marching direction.
  !    
  !     Parameters
  !     ----------
  !     L : integer
  !         The marching level to compute stretch for
  !
  !     Returns
  !     -------
  !     deltaS : real
  !         Marching increment set in hypData
  !
  !     Updates
  !     -------
  !     maxKStretch: real (defined in hypData.F90)
  !         Ratio of distance between layers
  !
  use precision
  use hypInput
  use hypData, only : maxKStretch, nx
  implicit none

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input parameters
  Vec, intent(in) :: XL, XLm1, XLm2

  ! Working
  integer(kind=intType) :: i, ierr
  real(kind=realType) :: rl(3), rlm1(3), rlm2(3), KStretch
  real(kind=realType), pointer, dimension(:) :: xxl, xxlm1, xxlm2
  real(kind=realType), external :: dist

  ! BEGIN EXECUTION

  ! Assign pointers
  call VecGetArrayF90(XL, xxl, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(XLm1, xxlm1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(XLm2, xxlm2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Initialize maximum stretch variable
  maxKStretch = 0.0

  ! Compute strech for every point
  do i=1,nx
     
     ! Get node coordinates on each layer
     rl = xxl(3*i-2:3*i)
     rlm1 = xxlm1(3*i-2:3*i)
     rlm2 = xxlm2(3*i-2:3*i)

     ! Compute stretch ratio (function dist defined in 3D_utilities.F90)
     KStretch = dist(rl, rlm1)/dist(rlm1, rlm2)

     ! Compare with the maximum value
     maxKStretch = max(KStretch, maxKStretch)

  end do

  ! Restore arrays to make petsc happy
  call VecRestoreArrayF90(XL, xxl, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(XLm1, xxlm1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecRestoreArrayF90(XLm2, xxlm2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine findKStretch

subroutine init3d(sizes, nn, nPtr_in, lmax, nNodes)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: init3d performs initialization task for the 3D
  !     code. This function is called from Python.
  !    
  !     Description of Arguments
  !     Input:
  !     sizes - integer array size(nn, 2) : Sizes of each of the patches
  !     nn - integer: The number of patches
  !     nPtr_in - integer array size(lMax, nGlobal) : Nodal pointer. This 
  !               structure is used to store the nodal neighbours in a 
  !               consistent counter clockwise ordering fashion. 
  !     lmax - integer: Maximum number of neighbours (first dimension 
  !            of nPtr_in
  !     nGlobal - N
  !     Xin - array size (nNodes, 3): The (unstructured) list of nodes that 
  !           this processor owns. 
  !     nNodes - integer : The number of nodes on this proc.
  !
  !     Ouput: None

  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: sizes(nn, 2), nn, lMax
  integer(kind=intType), intent(in) :: nPtr_in(lMax, nNodes), nNodes
  ! Working Parameters
  integer(kind=intType) :: i

  ! Set the number of patches:
  nPatch = nn

  ! Allocate patch list of patch types
  if (.not. allocated(patches)) then
     allocate(patches(nPatch))
  end if

  ! Set sizes
  do i=1,nPatch
     patches(i)%il = sizes(i, 1)
     patches(i)%jl = sizes(i, 2)
  end do

  ! Allocate the global node pointer structure
  if (.not. allocated(nPtr)) then
     allocate(nptr(lmax, nNodes))
  end if

  ! And copy in our python data
  do i=1, nNodes
     nPtr(:, i) = nPtr_in(:, i) 
  end do

end subroutine init3d

subroutine setLindex(ind_in, nnx, nny, iPatch)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: sets the local index data for each patch from python.
  !    
  !     Description of Arguments
  !     Input:
  !     ind_in - integer array size(nnx, nny) : Indices to set.
  !     nnx - integer: First dimension of the patch
  !     nny - integer: Second dimension of the patch
  !     iPatch - integer: Index of the patch to set. 
  !
  !     Ouput: None
  use hypInput
  use hypData

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: ind_in(nnx, nny), nnx, nny, iPatch

  ! Working Parameters
  integer(kind=intType) :: i, j

  ! Set the ind_in into the proper patch l_index. Note we will convert
  ! to fortran ordering here by adding a 1

  if (.not. allocated(patches(iPatch + 1)%l_index)) then
     allocate(patches(iPatch + 1)%l_index( &
          patches(iPatch + 1)%il, patches(iPatch + 1)%jl))
  end if

  do j=1, nny
     do i=1, nnx
        patches(iPatch + 1)%l_index(i, j) = ind_in(i, j) + 1
     end do
  end do

end subroutine setlindex

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


subroutine calcGridRatio()
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: calcGridRatio() calculates the exponential
  !     distribution Turns out we need to solve a transendental
  !     equation here. We will do this with a secant iteration.


  use hypInput
  use hypData

  implicit none

  ! Working Parameters
  integer(kind=intType) :: i, M
  real(kind=realType) ::  s, r, a,b, c, f, fa, fb

  ! function 'f' is S - s0*(1-r^n)/(1-r) where S is total length, s0 i
  ! sinitial ratio and r is the grid ratio. 

  S = radius0*rMin
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

  ! Finally set the gridRatio variable to r
  gridRatio = c

end subroutine calcGridRatio

subroutine three_by_three_inverse(Jac, Jinv)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: return the inverse of the 3x3 matrix Jac in Jinv
  !    
  !     Description of Arguments
  !     Input:
  !     Jac - array size (3,3): Matrix to invert
  !
  !     Ouput:
  !     Jinv - array size (3,3) : Inverse of Jac
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
  do i=2,125
     write(*,"(a)",advance="no") "-"
  enddo
  print "(1x)"
  write(*,"(a)",advance="no") "# Grid  |     CPU    | Sub  | KSP  |     Sl     |"
  write(*,"(a)",advance="no") " Grid       | Grid       |     Min    |   deltaS   |    cMax   |   min R   | "
  print "(1x)"
  write(*,"(a)",advance="no") "# Level |     Time   | Iter | Its  |            |"
  write(*,"(a)",advance="no") " Sensor Max | Sensor Min |  Quality   |            |           |           | "
  print "(1x)"

  ! Write the Last line of '-'
  write(*,"(a)",advance="no") "#"
  do i=2,125
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
  write(*,"(e12.5,1x)",advance="no") dble(mpi_wtime()) - timeStart

  ! Sub Iteration
  write(*,"(i6,1x)",advance="no") nSubIter

  ! KSP ITerations
  write(*,"(i6,1x)",advance="no") kspIts

  ! Modified Sl factor from Steger and Chan
  write(*,"(e12.5,1x)",advance="no") Sl

  ! maximum Grid convergence sensor
  write(*,"(e12.5,1x)",advance="no") gridSensorMax

  ! minimum Grid convergence sensor
  write(*,"(e12.5,1x)",advance="no") gridSensorMin

  ! minimum grid quality for new layer
  write(*,"(e12.5,1x)",advance="no") minQuality

  ! marching step size for this iteration
  write(*,"(e12.5,1x)",advance="no") deltaS

  ! marching step size for this iteration
  write(*,"(e12.5,1x)",advance="no") cRatio

  ! approximate minimim distance to body
  minR = zero
  write(*,"(e12.5,e12.5,e12.5,1x)",advance="no") scaleDist/radius0
  print "(1x)"

end subroutine writeIteration

subroutine computeR(xVec, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     computeR0 determines the minimum radius of a sphere that will
  !     fully enclose the grid level defined by xVec
  !
  !     Description of Arguments
  !     Input:
  !     xVec - Petsc vector : Vector of nodes to use in calculation
  !
  !     Ouput: 
  !     R - real : Minimum radius surrounding xVec
  !
  use precision
  use hypData, only : nx, xAvg
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
  real(kind=realType), dimension(:), pointer :: xx
  real(kind=realType) :: dist, xpt(3)
  ! First find the average of X
  Xavg = zero

  call VecGetArrayF90(xVec, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do i=1,nx
     Xavg = Xavg + XX(3*(i-1)+1:3*i)
  end do

  Xavg = Xavg/nx

  ! Now find the maximum distance of X from Xavg
  R = zero
  do i=1,nx
     xpt = xx(3*i-2:3*i)
     R = max(R, dist(Xavg, xpt))
  end do

  ! Always remember to restore the array
  call VecRestoreArrayF90(xVec, xx, ierr)
end subroutine computeR

subroutine computeMinR(xVec, R)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: computeMinR determines the cloest point of xVec to
  !     Xavg
  !
  !     Description of Arguments
  !     Input:
  !     xVec - Petsc vector : Vector of nodes to use in calculation
  !
  !     Ouput: 
  !     R - real : Minimum distance from xVec to Xavg

  use precision
  use hypData, only : nx, xAvg
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
  real(kind=realType), dimension(:), pointer :: xx
  real(kind=realType) :: dist, xpt(3)

  call VecGetArrayF90(xVec, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now find the minimum distance of X from Xavg
  R = huge(R)
  do i=1,nx
     xpt = xx(3*i-2:3*i)
     R = min(R, dist(Xavg, xpt))
  end do

  ! Always remember to restore the array
  call VecGetArrayF90(xVec, xx, ierr)

end subroutine computeMinR

subroutine computeQualityLayer
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

  ! Working
  integer(kind=intType) :: i, j, iPatch, ierr
  real(kind=realType) :: points(3,8), Q

  ! Since we are dealing with volumes, we need to loop over faces:
  minQuality = one
  call VecGetArrayF90(X(marchIter), xxtmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(X(marchIter-1), xxm1tmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do iPatch=1, nPatch
     do j=1,patches(iPatch)%jl-1
        do i=1,patches(iPatch)%il-1
           ! Assemble the 8 points we need for the volume - Note
           ! coodinate coordinate ordering!
           points(:,1) = xxm1tmp(3*patches(iPatch)%l_index(i  ,j  )-2:3*patches(iPatch)%l_index(i  ,j  ))
           points(:,2) = xxm1tmp(3*patches(iPatch)%l_index(i+1,j  )-2:3*patches(iPatch)%l_index(i+1,j  ))
           points(:,3) = xxm1tmp(3*patches(iPatch)%l_index(i  ,j+1)-2:3*patches(iPatch)%l_index(i  ,j+1))
           points(:,4) = xxm1tmp(3*patches(iPatch)%l_index(i+1,j+1)-2:3*patches(iPatch)%l_index(i+1,j+1))

           points(:,5) = xxtmp(3*patches(iPatch)%l_index(i  ,j  )-2:3*patches(iPatch)%l_index(i  ,j  ))
           points(:,6) = xxtmp(3*patches(iPatch)%l_index(i+1,j  )-2:3*patches(iPatch)%l_index(i+1,j  ))
           points(:,7) = xxtmp(3*patches(iPatch)%l_index(i  ,j+1)-2:3*patches(iPatch)%l_index(i  ,j+1))
           points(:,8) = xxtmp(3*patches(iPatch)%l_index(i+1,j+1)-2:3*patches(iPatch)%l_index(i+1,j+1))

           call quality_hexa(points, Q)
           minQuality = min(minQuality, Q)

        end do
     end do
  end do
  call VecRestoreArrayF90(X(marchIter), xxtmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecRestoreArrayF90(X(marchIter-1), xxm1tmp, ierr)
  call EChk(ierr, __FILE__, __LINE__)


end subroutine computeQualityLayer

subroutine quality_hexa(points, quality)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: quality_hexa determine the 2x2x2 quality measure for
  !     a hexa defined by 8 nodes. 
  !    
  !     Description of Arguments
  !     Input:
  !     points - array size(3, 8): Nodes defining the hexa in coordinate ording. 
  !
  !     Ouput: 
  !     Quality - real : Quality measure

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

subroutine volume_hexa(points, volume)
  use precision

  implicit none

  ! Input/Output
  real(kind=realType), intent(in) :: points(3, 8)
  real(kind=realType), intent(out) :: volume

  ! Working
  real(kind=realType) :: center(3), volpymrid, v1, v2, v3, v4, v5, v6
  integer(kind=intType) ::  i, idim

  ! Compute the center of the points
  center = zero
  do i=1, 8
     do idim=1, 3
        center(idim) = center(idim) + points(idim, i)
     end do
  end do
  center = center / eight

  ! Compute the volumes of the 6 sub pyramids. The
  ! arguments of volpym must be such that for a (regular)
  ! right handed hexahedron all volumes are positive.

  v1 = volpymrid(center, points(:, 1), points(:, 2), points(:, 4), points(:, 3))
  v2 = volpymrid(center, points(:, 7), points(:, 8), points(:, 6), points(:, 5))
  v3 = volpymrid(center, points(:, 1), points(:, 3), points(:, 7), points(:, 5))
  v4 = volpymrid(center, points(:, 4), points(:, 2), points(:, 6), points(:, 8)) 
  v5 = volpymrid(center, points(:, 2), points(:, 1), points(:, 5), points(:, 6)) 
  v6 = volpymrid(center, points(:, 3), points(:, 4), points(:, 8), points(:, 7))

  volume = (v1+v2+v3+v4+v5+v6)/six

end subroutine volume_hexa

function volpymrid(p, a, b, c, d)

  use precision
  implicit none
  real(kind=realType) :: p(3), a(3), b(3), c(3), d(3), volpymrid

  ! 6*Volume of a pyrimid -> Counter clockwise ordering
  volpymrid = (p(1) - fourth*(a(1) + b(1)  + c(1) + d(1))) *&
       ((a(2) - c(2))*(b(3) - d(3)) - (a(3) - c(3))*(b(2) - d(2)))   + &
       (p(2) - fourth*(a(2) + b(2)  + c(2) + d(2)))*&
       ((a(3) - c(3))*(b(1) - d(1)) - (a(1) - c(1))*(b(3) - d(3)))   + &
       (p(3) - fourth*(a(3) + b(3)  + c(3) + d(3)))* &
       ((a(1) - c(1))*(b(2) - d(2)) - (a(2) - c(2))*(b(1) - d(1)))

end function volpymrid

function dist(p1, p2)
  use precision

  real(kind=realType) :: p1(3), p2(3)
  real(kind=realType) :: dist
  dist = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)

end function dist


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

subroutine pointReduceWrap(pts, N, tol, uniquePts, link, nUnique)
  
  ! Python wrapped versino of pointReduce,

  use precision
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: N
  real(kind=realType), intent(in), dimension(3, N) :: pts
  real(kind=realType), intent(in) :: tol

  ! Output Parametres
  real(kind=realType), intent(out), dimension(3, N) :: uniquePts
  integer(kind=intType), intent(out), dimension(N) :: link
  integer(kind=intType), intent(out) :: nUnique

  interface 
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none
       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduce
  end interface

  call pointReduce(pts, N, tol, uniquePts, link, nUnique)
  
end subroutine pointReduceWrap

subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

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

  ! Working Parameters
  real(kind=realType), allocatable, dimension(:) :: dists, tmp
  integer(kind=intType), allocatable, dimension(:) :: ind
  integer(kind=intType) :: i, j, nTmp, link_counter, ii
  logical cont, cont2
  integer(kind=intType) :: tmpInd(N), subLink(N), nSubUnique
  real(kind=realType) :: subPts(3, N), subUniquePts(3, N)
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

  ! Allocate dists, and the ind pointer
  allocate(dists(N), tmp(N), ind(N))

  ! Compute distances of all points from the origin
  do i=1, N
     !dists(i) = sqrt((pts(1,i)-pts(1,1))**2 + (pts(2,i)-pts(2,1))**2 + (pts(3, i)-pts(3,1))**2)
     dists(i) = sqrt(pts(1,i)**2 + pts(2,i)**2 + pts(3, i)**2)
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

           tmpInd(nTmp) = ind(j)
           j = j + 1
           if (j == N+1) then ! Overrun check
              cont2 = .False.
           end if
        else
           cont2 = .False.
        end if
     end do

     ! Copy the points that have the same distance into subPts. Note
     ! these may NOT be the same, since two points can have the same
     ! distnace, but not be co-incident (ie (1,0,0), (-1,0,0))
     do ii=1,nTmp
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
  deallocate(dists, tmp, ind)

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

! subroutine reparameterize(nx)
!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     reparameterize() takes the pseudo grid defined in pgrid3d and
!   !     computes the final desired distribution according to N, s0,
!   !     and gridRatio. (rMin was already used since that determined
!   !     when the marching stopped'. When we are done, the final
!   !     desired grid will be in grid3D with the corect dimensions
!   !
!   use hypData
!   use hypInput
!   implicit none

!   ! Input
!   integer(kind=intType), intent(in) :: nx

!   ! Working
!   integer(kind=intType) :: i, l, bin
!   real(kind=realType) :: s(Nlayers), sAvg(nLayers), lAvg, R, sp1, factor, sp

!   ! The first step is to compute a averaged global parameterization of
!   ! the k-direction

!   savg = zero
!   lAvg = zero
!   do i=1,nx
!      s(1) = zero
!      do l=2,nLayers
!         s(l) = s(l-1) + dist(pGrid3d(:, i, l), pGrid3d(:, i, l-1))
!      end do

!      ! Add the length of this curve:
!      lAvg = lAvg + s(nLayers)

!      ! Normalize
!      s = s / s(nLayers)
!      ! Add to running total
!      savg = savg + s
!   end do

!   ! Final s average
!   savg = savg/nx

!   ! Final averge length
!   lavg = lAvg/nx

!   ! Compute 'r' for the exponential distribution
!   sp1 = s0/lavg
!   R = -log((N-1)*sp1)/(N-2)

!   ! Allocate grid3d

!   ! Master loop over final levels:
!   allocate(grid3D(3, nx, N))
!   grid3D = zero

!   ! Special case for l=1 --- Always the boundary:
!   grid3d(:, :, 1) = pGrid3d(:, :, 1)

!   do l=2,N

!      ! Generate the position where we need to interpolate
!      sp = sp1*(l-1)*exp(R*(l-2))

!      ! Now what we want to find is the "bin" in l where this spacing
!      ! falls. What we want is something like l=2.75. This means we
!      ! interpolate three quarters of the way between l=2 and l=3. 

!      call FindBin(sp, savg, Nlayers, bin)

!      ! Bin is index of the lower level, bin+1 is index of the upper level

!      ! Now get the factor
!      factor = (sp-savg(bin))/(savg(bin+1)-savg(bin))

!      ! Now in vector form:
!      grid3d(:,:,l) = (one-factor)*pGrid3D(:,:,bin) + factor*pGrid3D(:,:,bin+1)
!   end do

! end subroutine reparameterize


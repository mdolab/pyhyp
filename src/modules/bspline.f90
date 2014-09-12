module bspline

  use precision 
  type surfType
     ! Holds the data for a B-spline surface
     integer(kind=intType) :: ku, kv, nCtlu, nCtlv
     real(kind=realType), dimension(:), allocatable :: tu, tv
     real(kind=realType), dimension(:, :, :), allocatable :: coef
  end type surfType

  type(surfType), dimension(:), allocatable :: surfs

contains
  
  subroutine evalSurface(surf, u, v, n, pts) 
    ! Simplified routine for evaluating points on the surface. 
    use precision
    implicit none

    ! Input/Output
    type(surfType), intent(in) :: surf
    real(kind=realType), intent(in), dimension(n) :: u, v
    real(kind=realType), intent(out), dimension(3, n) :: pts
    integer(kind=intType), intent(in):: n
    ! Working
    integer(kind=intType) :: ileftu, ileftv, istartu, istartv, ii, i, j, idim
    real(kind=realType) :: basisu(surf%ku)
    real(kind=realType) :: basisv(surf%kv)
       
    do ii=1, n
       ! U
       call findSpan(u(ii), surf%ku, surf%tu, surf%nctlu, ileftu)
       call basis(surf%tu, surf%nctlu, surf%ku, u(ii), ileftu, basisu)
       istartu = ileftu-surf%ku

       ! V
       call findSpan(v(ii), surf%kv, surf%tv, surf%nctlv, ileftv)
       call basis(surf%tv, surf%nctlv, surf%kv, v(ii), ileftv, basisv)
       istartv = ileftv-surf%kv
       pts(:, ii) = 0.0_8
       do i=1, surf%ku
          do j=1, surf%kv
             do idim=1, 3
                pts(idim, ii) = pts(idim, ii) + &
                     basisu(i)*basisv(j)*surf%coef(idim, istartv+j, istartu+i)
             end do
          end do
       end do
    end do

  end subroutine evalSurface

  subroutine findSpan(u, k, t, nctl, ind)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway. Adapted from "The NURBS BooK" Algorithm A2.1
    !
    !     Abstract: Determine the knot span index

    !     Description of Arguments
    !     Input:
    !     u       - Real, parametric location we are looking for
    !     k       - Integer, order of B-spline 
    !     t       - Real, size(nctl + k), knot vector
    !
    !     Ouput: 
    !     ind     - Integer, knot span index

    use precision
    implicit none

    ! Input 
    integer(kind=intType), intent(in) :: k, nctl
    real(kind=realType), intent(in) :: u, t(nctl+k)
    ! Output
    integer(kind=intType), intent(out) :: ind

    ! Working
    integer(kind=intType) :: low, mid, high

    if (u >= t(nctl+1)) then
       ind = nctl
    else if (u < t(k)) then
       ind = k
    else
       low = k
       high = nctl+1

       ! Do a binary search
       mid = (low+high)/2
       do while ( (u < t(mid) .or. u >= t(mid+1)))
          if (u < t(mid)) then
             high = mid
          else
             low = mid
          end if

          mid = (low+high)/2
       end do

       ind = mid
    end if
  end subroutine findSpan

  subroutine basis(t, nctl, k, u, ind, B)

    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: basis evaluates the standard b-spline basis
    !               functions at a location x for knot vector t. Adapted
    !               from The NURBS Book, algoritm A2.2

    !     Description of Arguments
    !     Input
    !     t       - Real, Vector of knots. Size (nctl + k)
    !     nctl    - Integer, number of knots
    !     k       - Integer, order of B-spline 
    !     u       - Real, location for evaluation of basis 
    !     ind     - Integer, position in knot vector such that t(ind) <= x
    !
    !     Ouput 
    !     B       - Real, Vector, The k basis vector values
    use precision
    implicit none

    ! Input
    real(kind=realType), intent(in)  :: t(nctl+k), u
    integer(kind=intType), intent(in)  :: nctl, k, ind

    ! Output
    real(kind=realType), intent(out) :: B(0:k-1)

    ! Working
    real(kind=realType) :: left(0:k-1), right(0:k-1), temp, saved
    integer(kind=intType) :: j, r, p

    ! To be consistent with algorithm in The NURBS Book we will use
    ! 0.0_8-based ordering here
    B(0) = 1.0_8
    p = k-1

    do j=1,p
       left(j)  = u - t(ind+1-j)
       right(j) = t(ind+j) - u
       saved = 0.0_8

       do r=0,j-1
          temp = B(r)/(right(r+1)+left(j-r))
          B(r) = saved+right(r+1)*temp
          saved = left(j-r)*temp
       end do

       B(j) = saved
    end do

  end subroutine basis
end module bspline


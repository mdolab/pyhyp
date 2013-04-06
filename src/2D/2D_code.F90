subroutine run2D(Xin, nx)
  
  use hypInput
  use hypData

  implicit none

#include "include/finclude/petsc.h"

  ! Input Parameters
  real(kind=realType) :: Xin(2, nx)
  integer(kind=intType) :: nx

  ! Working parameters
  integer(kind=intType) :: i, j, idim

  ! First thing we will do is allocate the final grid array:
  if (allocated(grid2D)) then
     deallocate(grid2D)
  end if

  allocate(grid2D(2, nx, N))

  grid2D = zero
  do j=1,N
     do i=1,nx
        grid2D(1, i, j) = dble(i)/(nx-1)
        grid2D(2, i, j) = dble(j)/(N-1)
     end do
  end do

end subroutine run2D

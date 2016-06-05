subroutine setNumberPatches(N)

  ! This subroutine loads a plot3d surface file

  use hypData
  implicit none

  ! Input Arguments
  integer(kind=intType), intent(in) :: N

  allocate(patchIO(N))
  nPatch = N
end subroutine setNumberPatches

subroutine setPatch(ii, xPatch, il, jl)
  use hypData
  implicit none

  ! Input Arguments
  integer(kind=intType), intent(in) :: ii, il, jl
  real(kind=realType), intent(in) :: xPatch(il, jl, 3)

  ! Working
  integer(kind=intType) :: i, j, iDim

  ! Allocate space for the grid coordinates on the patch and read
  allocate(patchIO(ii)%X(3, il, jl))
  patchIO(ii)%il = il
  patchIO(ii)%jl = jl

  do i=1, il
     do j=1, jl
        do idim=1,3
           patchIO(ii)%X(idim, i, j) = XPatch(i, j, idim)
        end do
     end do
  end do

end subroutine setPatch

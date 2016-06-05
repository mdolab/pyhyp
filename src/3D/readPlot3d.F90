subroutine readPlot3d(fileName)

  ! This subroutine loads a plot3d surface file

  use hypData
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: fileName

  ! Working
  integer(kind=intType), dimension(:, :), allocatable :: patchSizes
  integer(kind=intType) :: ii, i, j

  ! Open file and read number of patches
  open(unit=7, form='formatted', file=fileName)
  read(7, *) nPatch

  allocate(patchSizes(3, nPatch), patchIO(nPatch))

  ! Read all block sizes
  read(7,*) (patchSizes(1, i), patchSizes(2, i), patchSizes(3, i), i=1, nPatch)

  ! Make sure that ALL k-index patch sizes are 1
  do ii=1, nPatch
     if (patchSizes(3, ii) /= 1) then
        print *,'Plot3d Read Error: k-dimension of all blocks must be precisely 1'
        stop
     end if
  end do

  ! Now allocate and read all the blocks from the plot3d file
  do  ii=1, nPatch
     patchIO(ii)%il = patchSizes(1, ii)
     patchIO(ii)%jl = patchSizes(2, ii)

     ! Allocate space for the grid coordinates on the patch and read
     allocate(patchIO(ii)%X(3, patchIO(ii)%il, patchIO(ii)%jl))

     read(7, *) (( patchIO(ii)%X(1, i, j), i=1,patchIO(ii)%il), j=1,patchIO(ii)%jl)
     read(7, *) (( patchIO(ii)%X(2, i, j), i=1,patchIO(ii)%il), j=1,patchIO(ii)%jl)
     read(7, *) (( patchIO(ii)%X(3, i, j), i=1,patchIO(ii)%il), j=1,patchIO(ii)%jl)
  enddo
  deallocate(patchSizes)

  ! All the reading is done
  close(7)

end subroutine readPlot3d

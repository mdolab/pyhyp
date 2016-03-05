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

  allocate(patchSizes(3, nPatch), patches(nPatch))

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
     patches(ii)%il = patchSizes(1, ii)
     patches(ii)%jl = patchSizes(2, ii)

     ! Allocate space for the grid coordinates on the patch and read
     allocate(patches(ii)%X(3, patches(ii)%il, patches(ii)%jl))
     allocate(patches(ii)%l_index(patches(ii)%il, patches(ii)%jl))
     allocate(patches(ii)%weights(patches(ii)%il, patches(ii)%jl))
     patches(ii)%weights(:, :) = zero

     read(7, *) (( patches(ii)%X(1, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
     read(7, *) (( patches(ii)%X(2, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
     read(7, *) (( patches(ii)%X(3, i, j), i=1,patches(ii)%il), j=1,patches(ii)%jl)
  enddo
  deallocate(patchSizes)

  ! All the reading is done
  close(7)

end subroutine readPlot3d

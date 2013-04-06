subroutine releaseMemory

  use hypData

  implicit none

  if (allocated(grid2D)) then
     deallocate(grid2D)
  end if


end subroutine releaseMemory

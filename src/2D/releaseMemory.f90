subroutine releaseMemory

  use hypData

  implicit none

  ! Working variables
  integer(kind=intType) :: iPatch

  if (allocated(grid2D)) then
     deallocate(grid2D)
  end if

  ! if (allocated(grid3D)) then
  !    deallocate(grid3D)
  ! end if

  if (allocated(pgrid3D)) then
     deallocate(pgrid3D)
  end if

  if (allocated(patches)) then
     do iPatch=1, nPatch
        if (allocated(patches(iPatch)%l_index)) then
           deallocate(patches(iPatch)%l_index) 
        end if
     end do

     ! Finally deallocate the patches array
     deallocate(patches)
  end if
end subroutine releaseMemory

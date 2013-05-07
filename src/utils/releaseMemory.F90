subroutine releaseMemory

  use hypData

  implicit none

  ! Working variables
  integer(kind=intType) :: iPatch

  if (allocated(grid2D)) then
     deallocate(grid2D)
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


subroutine destroyPetscVars

  use hypData

  implicit none
  
  integer(kind=intType) :: ierr

  ! Destroy hyp system objects
  call KSPDestroy(hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatDestroy(hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(hypDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine destroyPetscVars

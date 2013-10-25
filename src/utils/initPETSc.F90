subroutine initPETSc
  ! Simple wrapper for initializing petsc in case it wasn't done in python  
  use hypData

  implicit none
  integer(kind=intType) :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine initPETSc

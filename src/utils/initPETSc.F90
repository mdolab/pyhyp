subroutine initPETSc(comm)
    ! Simple wrapper for initializing petsc in case it wasn't done in
    ! python. This also get the comm data.
    use communication
    use hypData
    implicit none

    ! Input params
    integer(kind=intType) :: comm

    ! Working variables
    integer(kind=intType) :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    hyp_comm_world = comm
    hyp_comm_self = mpi_comm_self
    call MPI_Comm_size(hyp_comm_world, nProc, ierr)
    call MPI_Comm_rank(hyp_comm_world, myid, ierr)

end subroutine initPETSc

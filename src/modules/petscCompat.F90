module petscCompat
! Include PETSc version macros
#include <petsc/finclude/petsc.h>
    ! Module for PETSc compatibility functions
    ! This module provides compatibility for various subroutines in PETSc
    ! This is because of the changes in the API between different versions of PETSc
    ! This module is used to ensure that the code works with different versions of PETSc
    use petscvec
    use petscmat
    use petscksp
    use precision
    implicit none
    save

contains

    subroutine VecGetArrayCompat(v, array, ierr)
        ! PETSc compatability routine to get array from vector
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecGetArrayF90(v, array, ierr)
#else
        call VecGetArray(v, array, ierr)
#endif

    end subroutine VecGetArrayCompat

    subroutine VecRestoreArrayCompat(v, array, ierr)
        ! PETSc compatability routine to restore array
        implicit none
        Vec :: v
        real(kind=realType), pointer :: array(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_LT(3,20,0)
        call VecRestoreArrayF90(v, array, ierr)
#else
        call VecRestoreArray(v, array, ierr)
#endif

    end subroutine VecRestoreArrayCompat

    subroutine MatCreateDenseCompat(comm, m, n, mGlobal, nGlobal, data, A, ierr)
        ! Create a dense matrix in a PETSc-version-compatible way.
        ! If optional `data` is provided, it will be passed to PETSc,
        ! otherwise the correct null constant will be used.
        implicit none
        integer, intent(in) :: comm
        integer(kind=intType), intent(in) :: m, n, mGlobal, nGlobal
        real(kind=realType), intent(in), optional, dimension(:) :: data
        Mat :: A
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,22,0)
        if (present(data)) then
            call MatCreateDense(comm, m, n, mGlobal, nGlobal, data, A, ierr)
        else
            call MatCreateDense(comm, m, n, mGlobal, nGlobal, PETSC_NULL_SCALAR_ARRAY, A, ierr)
        end if
#else
        if (present(data)) then
            call MatCreateDense(comm, m, n, mGlobal, nGlobal, data, A, ierr)
        else
            call MatCreateDense(comm, m, n, mGlobal, nGlobal, PETSC_NULL_SCALAR, A, ierr)
        end if
#endif

    end subroutine MatCreateDenseCompat

    subroutine PCASMGetSubKSPCompat(pc, nlocal, first, subKSP, ierr)
        implicit none
        PC, intent(in) :: pc
        integer(kind=intType), intent(out) :: nlocal, first
        KSP, pointer, intent(out) :: subKSP(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
        ! PETSc >= 3.23, just call directly
        call PCASMGetSubKSP(pc, nlocal, first, subKSP, ierr)
#else
        ! Older PETSc, first query to get the size, then allocate, then call again
        call PCASMGetSubKSP(pc, nlocal, first, PETSC_NULL_KSP, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (nlocal > 0) then
            allocate (subKSP(nlocal))
            call PCASMGetSubKSP(pc, nlocal, first, subKSP, ierr)
        else
            nullify (subKSP)
        end if
#endif

    end subroutine PCASMGetSubKSPCompat

    subroutine PCASMRestoreSubKSPCompat(pc, nlocal, first, subKSP, ierr)

        implicit none
        PC, intent(in) :: pc
        integer(kind=intType), intent(in) :: nlocal, first
        KSP, pointer, intent(inout) :: subKSP(:)
        integer(kind=intType), intent(out) :: ierr

#if PETSC_VERSION_GE(3,23,0)
#if PETSC_VERSION_GE(3,23,0) && PETSC_VERSION_LT(3,24,0)
        ! PETSc 3.23 specific handling
        ! Note on PETSc 3.23 : It appears that PETSc does not have PCASMRestoreSubKSP, this is added in 3.24.
        ierr = 0
#else
        ! PETSc >= 3.24, just call directly
        call PCASMRestoreSubKSP(pc, nlocal, first, subKSP, ierr)
#endif
#else
        ! Older PETSc, check if pointer is associated, deallocate
        if (associated(subKSP)) deallocate(subKSP)
        ierr = 0
#endif

    end subroutine PCASMRestoreSubKSPCompat

end module petscCompat

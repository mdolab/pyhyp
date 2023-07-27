subroutine releaseMemory

    use communication
    use hypData
    use hypInput
    implicit none

    ! Working variables
    integer(kind=intType) :: iPatch, ierr, i, j

    if (allocated(patches)) then
        do iPatch = 1, nPatch
            deallocate (patches(iPatch)%l_index)
            deallocate (patches(iPatch)%X)
            deallocate (patches(iPatch)%weights)
        end do
        deallocate (patches)
    end if

    if (allocated(fullDeltaS)) then
        deallocate (fullDeltaS)
    end if

    if (allocated(gnPtr)) then
        deallocate (gnPtr)
    end if

    if (allocated(lnPtr)) then
        deallocate (lnPtr)
    end if

    if (allocated(families)) then
        deallocate (families)
    end if

    if (allocated(cptr)) then
        deallocate (cptr)
    end if

    if (allocated(topotype)) then
        deallocate (topotype)
    end if

    if (allocated(bctype)) then
        deallocate (bctype)
    end if

    if (allocated(bcval)) then
        deallocate (bcval)
    end if

    if (allocated(ghost)) then
        deallocate (ghost)
    end if

    if (allocated(conn)) then
        deallocate (conn)
    end if

    if (allocated(epseschedule)) then
        deallocate (epseschedule)
    end if

    if (allocated(epsischedule)) then
        deallocate (epsischedule)
    end if

    if (allocated(thetaschedule)) then
        deallocate (thetaschedule)
    end if

    if (allocated(splayschedule)) then
        deallocate (splayschedule)
    end if

    if (allocated(volSmoothSchedule)) then
        deallocate (volSmoothSchedule)
    end if

    if (allocated(volblendschedule)) then
        deallocate (volblendschedule)
    end if

    if (allocated(growthRatios)) then
        deallocate (growthRatios)
    end if

    if (varsAllocated) then
        ! Destroy hyp system objects
        call KSPDestroy(hypKSP, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatDestroy(hypMat, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(hypRHS, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(hypDelta, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(XL, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(XLm1, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(XLocal, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(XLocalm1, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecScatterDestroy(rootScatter, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Detrory all the PETsc vectors
        do i = 1, N
            call vecDestroy(X(i), ierr)
            call EChk(ierr, __FILE__, __LINE__)

            if (metricsAllocated) then
                do j = 1, nMetric
                    call VecDestroy(metrics(i, j), ierr)
                    call EChk(ierr, __FILE__, __LINE__)
                end do
            end if
        end do
        deallocate (X)

        if (metricsAllocated) then
            deallocate (metrics)
            metricsAllocated = .False.
        end if

        varsAllocated = .False.
    end if

end subroutine releaseMemory

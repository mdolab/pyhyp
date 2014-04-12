subroutine releaseMemory

  use communication
  use hypData
  use hypInput
  implicit none

  ! Working variables
  integer(kind=intType) :: iPatch, ierr, i, j
  
  if (allocated(grid2D)) then
     deallocate(grid2D)
  end if

  if (allocated(patches)) then
     do iPatch=1, nPatch
        deallocate(patches(iPatch)%l_index) 
     end do
     deallocate(patches)
  end if

  if (allocated(gnPtr)) then
     deallocate(gnPtr)
  end if

  if (allocated(lnPtr)) then
     deallocate(lnPtr)
  end if

  if (allocated(Xsurf)) then
     deallocate(Xsurf) 
  end if
 
  if (three_d_vars_allocated) then
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
     do i=1,N
        call vecDestroy(X(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (metricsAllocated) then
           do j=1,nMetric
              call VecDestroy(metrics(i, j), ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
        end if
     end do
     deallocate(X)

     if (metricsAllocated) then
        deallocate(metrics)
        metricsAllocated = .False.
     end if

     !deallocate(volume)
     three_d_vars_allocated = .False.
  end if

  if (two_d_vars_allocated) then
     ! Destroy hyp system objects
     call KSPDestroy(hypKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MatDestroy(hypMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecDestroy(hypRHS, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecDestroy(hypDelta, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

end subroutine releaseMemory

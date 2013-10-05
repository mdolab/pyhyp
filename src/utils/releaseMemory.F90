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
  call destroyPetscVars

end subroutine releaseMemory


subroutine destroyPetscVars
  use hypInput
  use hypData

  implicit none
  
  integer(kind=intType) :: ierr, i

  ! Destroy hyp system objects
  call KSPDestroy(hypKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatDestroy(hypMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(hypRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(hypDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  do i=1,N
     call vecDestroy(X(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     if (writeMetrics) then
        call VecDestroy(X_ksi(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(X_eta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(X_zeta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(X_ksi_ksi(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(X_eta_eta(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call VecDestroy(X_diss(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(Vhist(i), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if
  end do

  deallocate(X)

  if (writeMetrics) then
     deallocate(X_ksi, X_eta, X_zeta, X_ksi_ksi, X_eta_eta, X_diss, Vhist)
  end if

  deallocate(xx, rr, xxm1, xxm2, inds, volume, xxinterp)

end subroutine destroyPetscVars

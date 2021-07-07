subroutine writePlot3d(fileName)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writePLot3d write the current grid to a 3D PLOT3D
  !               file. It does not make any attempt to do grid
  !               connectivities or boundary conditions.
  !
  !     Parameters
  !     ----------
  !     fileName : Character array
  !         The name of the plot3d file
  use communication
  use hypInput
  use hypData
  implicit none

  ! Input Arguments
  character*(*) :: fileName

  ! Working variables
  integer(kind=intType) :: i, j, k, idim, iPatch, idGlobal, ierr

  if (myid == 0) then
     ! Open the output file
     open(unit=7, file=fileName)
     
     ! Write the number of zones
     write(7,5) nPatch
5    format (I4)
     
     ! Write the zone sizes
     do i=1,nPatch
        write(7,6) patches(i)%il, patches(i)%jl, N
     end do
6    format(I5, I5, I5)
  end if

  do iPatch = 1, nPatch
     do idim=1,3
        do k=1, N

           call VecScatterBegin(rootScatter, X(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
           call EChk(ierr,__FILE__,__LINE__)
           call VecScatterEnd(rootScatter, X(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
           call EChk(ierr,__FILE__,__LINE__)

           if (myid == 0) then
              call VecGetArrayF90(xLocal, xx, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           
              do j=1, patches(iPatch)%jl
                 do i=1, patches(iPatch)%il
                    ! We will use l_index to get the correct pointer into grid3D
                    idGlobal = patches(iPatch)%l_index(i, j)
                    write(7,7) xx(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(XLocal, xx, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end if
        end do
     end do
  end do
7 format(g21.14)

  if (myid == 0) then
     ! Finally close the file
     close(7)
  end if
end subroutine writePlot3d

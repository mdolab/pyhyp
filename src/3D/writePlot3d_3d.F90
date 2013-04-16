subroutine writePlot3d_3d(fileName)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writePlot3d_3d writes the current grid to a (volume)
  !               plot3d (ascii) file
  !
  !     Description of Arguments
  !     Input:
  !     fileNmae - Character array: the name of the plot3d file
  !
  !     Ouput: None

  use hypInput
  use hypData
  implicit none

  ! Input Arguments
  character*(*) :: fileName

  ! Working variables
  integer(kind=intType) :: i, j, k, idim, iPatch, idGlobal

  ! Open the output file
  open(unit=7, file=fileName)
  
  ! Write the number of zones
  write(7,5) nPatch
5 format (I4)

  ! Write the zone sizes
  do i=1,nPatch
     write(7,6) patches(i)%il, patches(i)%jl, NDebug
  end do
6 format(I5, I5, I5)

  ! Loop over each patch and write
  do iPatch = 1, nPatch
     do idim=1,3
        do k=1, NDebug
           do j=1, patches(iPatch)%jl
              do i=1, patches(iPatch)%il
                 ! We will use l_index to get the correct pointer into grid3D
                 idGlobal = patches(iPatch)%l_index(i, j)
                 write(7,7), grid3d(iDim, idGlobal, k)
              end do
           end do
        end do
     end do
  end do

7 format(g21.14)

  ! Finally close the file
  close(7)
end subroutine writePlot3d_3d

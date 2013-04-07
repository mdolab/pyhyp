subroutine writePlot3d_2d(fileName)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writePlot3d_2d write the current grid to a (surface)
  !               plot3d (ascii) file
  !
  !     Description of Arguments
  !     Input:
  !     fileNmae - Character array: the name of the plot3d file
  !
  !     Ouput: None

  use hypData

  ! Input Arguments
  character*(*) :: fileName

  ! Working variables
  integer(kind=intType) :: i, j, gridShp(3)

  ! Open the output file
  open(unit=7, file=fileName)
  
  ! Write the number of zones
  write(7,5) 1
5 format (I1)

  ! Write the zone sizes
  gridShp = shape(grid2D)
  write(7,6) gridShp(2)+1, gridShp(3), 1
6 format(I5, I5, I5)

  do idim=1,2
     do j=1,gridShp(3)
        do i=1,gridShp(2)
           write(7,7), grid2D(idim, i, j)
        end do
        write(7,7), grid2D(idim, 1, j)
     end do
  end do

  ! Write out the z-coordinate which is just zero
  do j=1,gridShp(3)
     do i=1,gridShp(2)+1
        write(7,7), zero
     end do
  end do

7 format(g21.14)

  ! Finally close the file
  close(7)
end subroutine writePlot3d_2d

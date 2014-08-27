subroutine writeLayerPlot3d(fileName, layer)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeLayerPlot3D writes a single mesh layer out to a plot3d file. 
  !    
  use communication
  use hypData
  use hypInput
  implicit none

  ! Input Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: layer

  ! Working parameters
  integer(kind=intType) :: i, j, ierr, idim, ipatch, idGlobal

  ! Send everything to the root:
  call VecScatterBegin(rootScatter, X(layer), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(rootScatter, X(layer), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
5    format (I4)
6    format(I5, I5, I5)
7    format(g21.14)

     call VecGetArrayF90(xLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! Open the output file
     open(unit=7, file=fileName)
     
     ! Write the number of zones
     write(7,5) nPatch

     ! Write the zone sizes
     do i=1,nPatch
        write(7,6) patches(i)%il, patches(i)%jl, 1
     end do
     do iPatch = 1, nPatch
        do idim=1,3
           
           do j=1, patches(iPatch)%jl
              do i=1, patches(iPatch)%il
                 ! We will use l_index to get the correct pointer
                 idGlobal = patches(iPatch)%l_index(i, j)
                 write(7,7), xx(3*(idGlobal-1)+iDim)
              end do
           end do
        end do
     end do
     close(7)
     call VecRestoreArrayF90(XLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)

  end if
end subroutine writeLayerPlot3d

subroutine writeLayerFE(fileName, layer, partitions)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeLayerFE writes a single mesh layer out to a tecplot file.
  !    
  use communication
  use hypData
  use hypInput
  implicit none

  ! Input Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: layer
  logical :: partitions

  ! Working parameters
  integer(kind=intType) :: i, j, ierr, idim, ipatch, idGlobal

  ! Send everything to the root:
  call VecScatterBegin(rootScatter, X(layer), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(rootScatter, X(layer), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
5    format (I4)
6    format(I5, I5, I5)
7    format(g21.14)

     call VecGetArrayF90(xLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     ! Open the output file
     open(unit=7, file=fileName)

101  format('ZONE NODES=', I5, ', ELEMENTS=', I5, ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL')
102  format(g20.12, g20.12, g20.12, g20.12)
103  format(I6,I6,I6,I6)

     write(7, 101), nxGlobal, faceTotal
     do i=1,nxGlobal
        write(7, 102) xx(3*i-2), xx(3*i-1), xx(3*i)
     end do
     do i=1,faceTotal
        write(7, 103) FullConn(1, i), FullConn(2, i), FullConn(3, i), FullConn(4, i)
     end do
     close(7)
     
     call VecRestoreArrayF90(XLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if
end subroutine writeLayerFE

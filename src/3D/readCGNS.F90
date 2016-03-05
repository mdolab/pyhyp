subroutine readCGNS(cgnsFile)

  use hypData
  use communication
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgnsFile

  ! Working
  integer(kind=intType) :: cg, ierr, i, j, k
  integer(kind=intType) :: CellDim, PhysDim, nZones, base, iZone, nbases
  integer(kind=intType) :: nbocos, iBC, bocoType, bocoIndex
  integer(kind=intType) :: iMinBC, iMaxBC, jMinBC, jMaxBC
  integer(kind=intType) :: NormalIndex(3), NormalListSize, NormalDataType, ndataset
  integer(kind=intType) :: ptset_type, npnts, pnts(2,2), normalList
  integer(kind=intType) :: zoneType, idataset
  integer(kind=intType) :: dummybocoType, DirichletFlag, NeumannFlag
  integer(kind=intType) :: dims(6)
  real(kind=realType), dimension(:, :), allocatable :: coor
  character*32 :: baseName, zoneName, bocoName, DatasetName

  if (myid == 0) then 
     ! Open and get the number of zones:
     call cg_open_f(trim(cgnsFile), CG_MODE_READ, cg, ierr)
     if (ierr /= CG_OK) call cg_error_exit_f

     call cg_nbases_f(cg, nbases, ierr)
     if (ierr /= CG_OK) call cg_error_exit_f

     if (nbases .gt. 1) then
        print *, ' ** Warning: pyHyp only reads the first base in a cgns file'
     end if

     base = 1_intType

     call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ierr)
     if (ierr /= CG_OK) call cg_error_exit_f

     if (cellDim .ne. 2 .or. PhysDim .ne. 3) then
        print *, 'The Cells must be 2 dimensional'
        stop
     end if

     call cg_nzones_f(cg, base, nZones, ierr);
     if (ierr /= CG_OK) call cg_error_exit_f

     ! Determine if we have all structured zones
     do i=1, nZones
        call cg_zone_type_f(cg, base, i, zoneType, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        if (zoneType == Unstructured) then 
           print *, 'Cannot do unstructured zones!'
           stop
        end if
     end do

     !Allocate patches
     allocate(patches(nZones))

     ! Populate patch with block information
     do iZone=1, nZones

        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        ! Store patch size
        patches(iZone)%il = dims(1)
        patches(iZone)%jl = dims(2)

        ! Allocate coordinates array
        allocate(patches(iZone)%X(3, dims(1), dims(2)))

        ! Allocate extra stuff
        allocate(patches(iZone)%l_index(dims(1), dims(2)))
        allocate(patches(iZone)%weights(dims(1), dims(2)))
        patches(iZone)%weights = 0

        ! Allocate auxiliary variable to load coordinates
        allocate(coor(dims(1), dims(2)))

        ! Load X coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateX", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patches(iZone)%X(1,i,j) = coor(i, j)
           end do
        end do

        ! Load Y coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patches(iZone)%X(2,i,j) = coor(i, j)
           end do
        end do

        ! Load Z coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patches(iZone)%X(3,i,j) = coor(i, j)
           end do
        end do

        !Deallocate coor for next iteration
        deallocate(coor)

     end do

     ! Close file
     call cg_close_f(cg, ierr)
     if (ierr /= CG_OK) call cg_error_exit_f
     nPatch = nZones
  end if

end subroutine readCGNS


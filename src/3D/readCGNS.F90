subroutine readCGNS(cgnsFile)

  use hypData
  use communication
  use cgnsGrid
  implicit none

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
  integer(kind=cgsize_t) :: dims(6)
  real(kind=realType), dimension(:, :), allocatable :: coor
  character(len=32) :: baseName, zoneName, bocoName, DatasetName

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

     !Allocate patchIO
     allocate(patchIO(nZones))

     ! Populate patch with block information
     do iZone=1, nZones

        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        ! Store patch size
        patchIO(iZone)%il = dims(1)
        patchIO(iZone)%jl = dims(2)

        ! Allocate coordinates array
        allocate(patchIO(iZone)%X(3, dims(1), dims(2)))

        ! Allocate auxiliary variable to load coordinates
        allocate(coor(dims(1), dims(2)))

        ! Load X coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateX", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patchIO(iZone)%X(1,i,j) = coor(i, j)
           end do
        end do

        ! Load Y coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patchIO(iZone)%X(2,i,j) = coor(i, j)
           end do
        end do

        ! Load Z coordinates into the auxilary variable (coor)
        call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, &
             (/1,1/), dims(1:2),  coor, ierr)
        if (ierr /= CG_OK) call cg_error_exit_f

        do j=1, dims(2)
           do i=1, dims(1)
              patchIO(iZone)%X(3,i,j) = coor(i, j)
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


subroutine readFamily(cgnsFile, iBlock, family, foundFam)

  use hypData
  use communication
  use cgnsGrid
  implicit none

  ! Input/Output Arguments
  character*(*),intent(in) :: cgnsFile
  integer(kind=intType), intent(in) :: iBlock
  character(len=32),intent(out) :: family
  logical, intent(out) :: foundFam

  ! Working
  integer(kind=intType) :: cg, ierr, i
  integer(kind=intType) :: base, nBoco
  integer(kind=intType) :: nbocos, iBC, bocoType, bocoIndex
  integer(kind=intType) :: NormalIndex(3), NormalDataType, ndataset
  integer(kind=intType) :: ptset_type,  normalList
  integer(cgsize_t) :: pts(2,2), npnts, NormalListSize
  character(len=32) :: zoneName, bocoName

  family = ""
  call cg_open_f(trim(cgnsFile), CG_MODE_READ, cg, ierr)
  if (ierr /= CG_OK) call cg_error_exit_f

  base = 1
  call cg_goto_f(cg, base, ierr, 'end')
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Get number of bocos
  call cg_nbocos_f(cg, base, iBlock, nBoco, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  foundFam = .False.
  do iBC=1, nBoco

     call cg_boco_info_f(cg, base, iBlock, iBC, bocoName, bocoType, &
          ptset_type, npnts, NormalIndex, NormalListSize, NormalDataType, nDataSet, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     if (bocoType == BCWall .or. bocoType == BCWallInviscid .or. &
          bocoType == BCWallViscous .or. bocoType == BCWallViscousHeatFlux .or. &
          bocoType == BCWallViscousIsothermal) then

        call cg_boco_read_f(cg, base, iBlock, iBC, pts, Integer, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Make sure it is a surface condition
        if (pts(1, 2)-pts(1, 1) > 0 .and. pts(2,2) - pts(2, 1) > 0) then
           ! It is a surface BC

           call cg_goto_f(cg, base, ierr, "Zone_t", iBlock, "ZoneBC_t",1, "BC_t", iBC, "end")
           if (ierr == 0) then ! Node exits
              if (.not. foundFam) then

                 ! Get family name
                 call cg_famname_read_f(family, ierr)
                 if (ierr .eq. CG_ERROR) call cg_error_exit_f
                 foundFam = .True.
              else
                 print *, 'Warning: multiple families found for block ', iBlock ,'Only the first is used.'
              end if
           end if
        end if
     end if
  end do

  ! Close file
  call cg_close_f(cg, ierr)
  if (ierr /= CG_OK) call cg_error_exit_f

end subroutine readFamily

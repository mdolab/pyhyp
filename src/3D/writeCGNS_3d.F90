subroutine writeCGNS_3D(fileName)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeCGNS_3d write the current grid to a 3D CGNS
  !               file. It computes grid connectivities such that the
  !               grid and boundary condition information such that
  !               the grid can be used directly in a 3D flow solver
  !               such  as SUmb.
  !
  !     Description of Arguments
  !     Input:
  !     fileName - Character array: the name of the cgns file
  !
  !     Ouput: None

  use hypInput
  use hypData

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*) :: fileName

  ! Working Variables
  integer(kind=intType) :: cg, ierr
  integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
  character*256 :: zoneName
  integer(kind=intType) :: sizes(9), zone, ii, i, j, k, transform(3)
  integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
  integer(kind=intType) :: nPatchToWrite, s, f
  real(kind=realType), dimension(:,:,:), allocatable :: coordArray, solArray
  character*12, dimension(3) :: coorNames, r_ksi, r_eta, r_zeta
  character*12, dimension(3) :: r_ksi_ksi, r_eta_eta, diss
  integer(kind=intType) :: iPatch, idim, idglobal

  ! Set All Names
  coorNames(1) = "CoordinateX"
  coorNames(2) = "CoordinateY"
  coorNames(3) = "CoordinateZ"     
  r_ksi(1) = 'x_ksi_x'
  r_ksi(2) = 'x_ksi_y'
  r_ksi(3) = 'x_ksi_z'

  r_eta(1) = 'x_eta_x'
  r_eta(2) = 'x_eta_y'
  r_eta(3) = 'x_eta_z'

  r_zeta(1) = 'x_zeta_x'
  r_zeta(2) = 'x_zeta_y'
  r_zeta(3) = 'x_zeta_z'

  r_ksi_ksi(1) = 'x_ksi_ksi_x'
  r_ksi_ksi(2) = 'x_ksi_ksi_y'
  r_ksi_ksi(3) = 'x_ksi_ksi_z'

  r_eta_eta(1) = 'x_eta_eta_x'
  r_eta_eta(2) = 'x_eta_eta_y'
  r_eta_eta(3) = 'x_eta_eta_z'

  diss(1) = 'diss_x'
  diss(2) = 'diss_y'
  diss(3) = 'diss_z'

  ! Open the CGNS File:
  call cg_open_f(fileName, MODE_WRITE, cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Write the single base
  cellDim = 3
  physDim = 3
  call cg_base_write_f(cg,"Base#1", Celldim, Physdim, base, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  if (writeMirror) then
     nPatchToWrite = nPatch
  else
     nPatchToWrite = nPatch/2
  end if

  do iPatch=1, nPatchToWrite
     ! Write the single zone
999  FORMAT('domain.',I5.5)
     write(zonename,999) iPatch ! Domains are typically ordered form 1
     sizes(:) = 0_intType
     sizes(1) = patches(iPatch)%il
     sizes(2) = patches(iPatch)%jl
     sizes(3) = N
     sizes(4) = sizes(1) - 1
     sizes(5) = sizes(2) - 1
     sizes(6) = sizes(3) - 1 
     call cg_zone_write_f(cg, base, zonename, sizes, Structured, &
          zone, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! Allocate the coord array
     allocate(coordArray(sizes(1),sizes(2),sizes(3)))
     allocate(solArray(sizes(1), sizes(2), sizes(3)))

     do idim=1,3
        ! Copy values and write:
        do k=1,N
           call VecGetArrayF90(X(k), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           do j=1,patches(iPatch)%jl
              do i=1,patches(iPatch)%il
                 idGlobal = patches(iPatch)%l_index(i,j)
                 coordArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
              end do
           end do
           call VecRestoreArrayF90(X(k), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end do
        call cg_coord_write_f(cg, base, zone, realDouble, &
             coorNames(iDim), coordArray, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
     end do

     if (writeMetrics) then
        call cg_sol_write_f(cg, base, zone, "FlowSolution", Vertex, s, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f


        ! ------------------ X_ksi --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_ksi(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_ksi(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                r_ksi(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ X_eta --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_eta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_eta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                r_eta(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ X_zeta --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_zeta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_zeta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                r_zeta(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ X_ksi_ksi --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_ksi_ksi(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_ksi_ksi(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                r_ksi_ksi(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ X_eta_eta --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_eta_eta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_eta_eta(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                r_eta_eta(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ X_diss --------------------
        do idim=1,3
           do k=1,N
              call VecGetArrayF90(X_diss(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    solArray(i,j,k) = xxtmp(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(X_diss(k), xxtmp, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end do
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                diss(idim), solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end do

        ! ------------------ Volume --------------------
        do k=1,N
           call VecGetArrayF90(Vhist(k), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           do j=1,patches(iPatch)%jl
              do i=1,patches(iPatch)%il
                 idGlobal = patches(iPatch)%l_index(i,j)
                 solArray(i,j,k) = xxtmp(3*idGlobal)
              end do
           end do
           call VecRestoreArrayF90(Vhist(k), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end do
        call cg_field_write_f(cg, base, zone, s, RealDouble, &
             'volume', solArray, F, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
     end if

     ! ! We know we can write the klow wall BC and the k-high farfield
     ! ! ---------- Klow --------------
     ! pnts(:, 1) = (/1, 1, 1/)
     ! pnts(:, 2) = (/sizes(1), sizes(2), 1/)

     ! call cg_boco_write_f(cg, base, zone, "klow", BCWallViscous, PointRange, &
     !      2, pnts, BCout, ierr)
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! call cg_goto_f(cg,base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
     !      "BC_t", BCOut, "end")
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! call cg_famname_write_f("wall", ierr)
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! ! ---------- Khigh --------------
     ! pnts(:, 1) = (/1, 1, sizes(3)/)
     ! pnts(:, 2) = (/sizes(1), sizes(2), sizes(3)/)

     ! call cg_boco_write_f(cg, base, zone, "khigh", BCFarField, PointRange, &
     !      2, pnts, BCout, ierr)
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! call cg_goto_f(cg,base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
     !      "BC_t", BCOut, "end")
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! call cg_famname_write_f("far", ierr)
     ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! Deallocate the coord and solution array for this patch
     deallocate(coordArray, solArray)

  end do ! Patch Loop

  ! Finally close cgns_file
  call cg_close_f(cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeCGNS_3D

subroutine zeroMirrorPlane(fileName, mirrorList, mirrorDim, nMirror)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: zeroMirrorPlane supplies EXACT zero coordinates to
  !               the mirror plane. It uses the inital topology to
  !               determine which edges (which have become faces) to
  !               zero. The mirrorList is and nx2 array containint the
  !               faceID and edgeID (in pySpline/pyGeo ordering) to
  !               zero. We will open the previously written CGNS file,
  !               read, the correct dimension, zero the correct
  !               plane, and the rewrite. 
  !
  !     Description of Arguments
  !     Input:
  !     fileName - Character array: the name of the cgns file to modify
  !     mirrorList - integer array, size(nMirror, 3): List of faces and node indices
  !                  (on surface) to mirror
  !     mirrorDim - Dimension to mirror, 1 for x, 2 for y, 3 for z.
  !     nMirror - integer: The number of faces to modify.
  !
  !     Ouput: None

  use hypInput
  use hypData

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: mirrorList(nMirror, 3)
  integer(kind=intType), intent(in) :: mirrorDim, nMirror

  ! Working Variables
  integer(kind=intType) :: cg, ierr, base, nzones, ii, zone, zonesize(3)
  integer(kind=intType) :: coordID, start(3), zoneType, loadedZone 
  character*32 zonename
  real(kind=realType), dimension(:, :, :), allocatable :: coor

  ! Don't do anything if no mirror nodes
  if (nMirror == 0) then
     return
  end if

  ! Open the CGNS File:
  call cg_open_f(fileName, MODE_MODIFY, cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  base = 1

  ! Goto Base Node
  call cg_goto_f(cg, base, ierr, 'end')
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  call cg_nzones_f(cg, base, nzones, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  loadedZone = -1
  do ii = 1,nMirror 
     ! Extract zone from the mirror list
     zone = mirrorList(ii, 1)

     if (zone .ne. loadedZone) then
        ! We need to switch zones. This can happen for two reasons: It
        ! is either the first zone or we just switched zones:

        if (loadedZone == -1) then ! First Pass:
           ! Read data for this one
           call cg_zone_read_f(cg, base, zone, zonename, zonesize, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           call cg_zone_type_f(cg, base, zone, zonetype, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Allocate enough space for one coordiante
           allocate(coor(zonesize(1),zonesize(2),zonesize(3)))

           ! Read the correct coordinate, based on mirrorDim which was
           ! passed in.
           start(:) = 1
           if (mirrorDim == 1) then
              call cg_coord_read_f(cg,base,zone,'CoordinateX',RealDouble,start,zonesize,coor,ierr)
           else if(mirrorDim == 2) then
              call cg_coord_read_f(cg,base,zone,'CoordinateY',RealDouble,start,zonesize,coor,ierr)
           else
              call cg_coord_read_f(cg,base,zone,'CoordinateZ',RealDouble,start,zonesize,coor,ierr)
           end if
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           loadedZone = zone
        else ! switching zones, write the current one and read the next:
           ! We need to write the data back for current zone:
           if (mirrorDim == 1) then
              call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateX",coor,coordID,ierr)
           else if(mirrorDim == 2) then
              call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateY",coor,coordID,ierr)
           else 
              call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateZ",coor,coordID,ierr)
           end if
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           deallocate(coor)

           ! And read the new zone:
           ! Read data for this one
           call cg_zone_read_f(cg, base, zone, zonename, zonesize, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           call cg_zone_type_f(cg, base, zone, zonetype, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Allocate enough space for one coordiante
           allocate(coor(zonesize(1),zonesize(2),zonesize(3)))

           ! Read the correct coordinate, based on mirrorDim which was
           ! passed in.
           start(:) = 1
           if (mirrorDim == 1) then
              call cg_coord_read_f(cg,base,zone,'CoordinateX',RealDouble,start,zonesize,coor,ierr)
           else if(mirrorDim == 2) then
              call cg_coord_read_f(cg,base,zone,'CoordinateY',RealDouble,start,zonesize,coor,ierr)
           else
              call cg_coord_read_f(cg,base,zone,'CoordinateZ',RealDouble,start,zonesize,coor,ierr)
           end if
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
           loadedZone = zone
        end if
     end if

     ! Now mirror the line based on mirro list
     coor(mirrorList(ii, 2), mirrorList(ii, 3), :) = zero

  end do

  ! Finally write the last zone:

  ! We need to write the data back for current zone:
  if (mirrorDim == 1) then
     call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateX",coor,coordID,ierr)
  else if(mirrorDim == 2) then
     call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateY",coor,coordID,ierr)
  else 
     call cg_coord_write_f(cg,base,loadedZone,RealDouble,"CoordinateZ",coor,coordID,ierr)
  end if
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Deallocate the coord array:
  deallocate(coor)
  call cg_close_f(cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

end subroutine zeroMirrorPlane


! subroutine writeCGNS_3DOrig(fileName)

!   !***DESCRIPTION
!   !
!   !     Written by Gaetan Kenway
!   !
!   !     Abstract: writeCGNS_3d write the current grid to a 3D CGNS
!   !               file. It computes grid connectivities such that the
!   !               grid and boundary condition information such that
!   !               the grid can be used directly in a 3D flow solver
!   !               such  as SUmb.
!   !
!   !     Description of Arguments
!   !     Input:
!   !     fileNmae - Character array: the name of the cgns file
!   !
!   !     Ouput: None

!   use hypInput
!   use hypData

!   implicit none
!   include 'cgnslib_f.h'

!   ! Input Arguments
!   character*(*) :: fileName

!   ! Working Variables
!   integer(kind=intType) :: cg, ierr
!   integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
!   character*256 :: zoneName
!   integer(kind=intType) :: sizes(9), zone, ii, i, j, k, nx, transform(3)
!   integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
!   integer(kind=intType) :: nPatchToWrite
!   real(kind=realType), dimension(:,:,:), allocatable :: coordArray
!   character*12, dimension(3) :: coorNames
!   integer(kind=intType) :: iPatch, idim, idglobal
!   ! Open the CGNS File:
!   call cg_open_f(fileName, MODE_WRITE, cg, ierr)
!   if (ierr .eq. CG_ERROR) call cg_error_exit_f

!   ! Write the single base
!   cellDim = 3
!   physDim = 3
!   call cg_base_write_f(cg,"Base#1", Celldim, Physdim, base, ierr)
!   if (ierr .eq. CG_ERROR) call cg_error_exit_f

!   if (writeMirror) then
!      nPatchToWrite = nPatch
!   else
!      nPatchToWrite = nPatch/2
!   end if

!   do iPatch=1, npatchToWrite
!      ! Write the single zone
! 999  FORMAT('domain.',I5.5)
!      write(zonename,999) iPatch ! Domains are typically ordered form 1
!      sizes(:) = 0_intType
!      sizes(1) = patches(iPatch)%il
!      sizes(2) = patches(iPatch)%jl
!      sizes(3) = NLayers
!      sizes(4) = sizes(1) - 1
!      sizes(5) = sizes(2) - 1
!      sizes(6) = sizes(3) - 1 
!      call cg_zone_write_f(cg, base, zonename, sizes, Structured, &
!           zone, ierr)
!      if (ierr .eq. CG_ERROR) call cg_error_exit_f

!      ! Allocate the coord array
!      allocate(coordArray(sizes(1),sizes(2),sizes(3)))
!      coorNames(1) = "CoordinateX"
!      coorNames(2) = "CoordinateY"
!      coorNames(3) = "CoordinateZ"     

!      do idim=1,3
!         ! Copy values and write:
!         do k=1,Nlayers
!            do j=1,patches(iPatch)%jl
!               do i=1,patches(iPatch)%il
!                  idGlobal = patches(iPatch)%l_index(i,j)
!                  coordArray(i,j,k) = pGrid3d(idim, idGlobal, k)
!               end do
!            end do
!         end do
!         call cg_coord_write_f(cg, base, zone, realDouble, &
!              coorNames(iDim), coordArray, coordID, ierr)
!      end do

!      ! Deallocate the coord array for this patch
!      deallocate(coordArray)
!   end do ! Patch Loop

!   ! Finally close cgns_file
!   call cg_close_f(cg, ierr)
!   if (ierr .eq. CG_ERROR) call cg_error_exit_f

! end subroutine writeCGNS_3DOrig

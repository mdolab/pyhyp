subroutine writeCGNS(fileName)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeCGNS write the current grid to a 3D CGNS
  !               file. It does not make any attempt to do grid
  !               connectivities or boundary conditions. Use
  !               cgns_utils from cgnsUtilities for that. 
  !
  !     Parameters
  !     ----------
  !     fileName : Character array
  !         The name of the cgns file

  use communication
  use hypInput 
  use hypData, only : X, metrics, rootScatter, xLocal, metricsAllocated
  use hypData, only : nPatch, patches, xx, metrics, ix_ksi, ix_eta, ix_zeta, iX_eta_eta, iX_ksi_ksi, iX_diss

  implicit none
  include 'cgnslib_f.h'
#include "include/finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input Arguments
  character*(*) :: fileName

  ! Working Variables
  integer(kind=intType) :: cg, ierr
  integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
  character*256 :: zoneName
  integer(kind=intType) :: sizes(9), zone, ii, i, j, k, cgtransform(3)
  integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
  integer(kind=intType) :: nPatchToWrite, s, f
  real(kind=realType), dimension(:,:,:), allocatable :: coordArray, solArray
  character*12, dimension(3) :: coorNames, r_ksi, r_eta, r_zeta
  character*12, dimension(3) :: r_ksi_ksi, r_eta_eta, diss
  integer(kind=intType) :: iPatch, idim, idglobal

  ! Set All Names
  coorNames = (/"CoordinateX", "CoordinateY", "CoordinateZ"/)
  r_ksi = (/'x_ksi_x', 'x_ksi_y', 'x_ksi_z'/)
  r_eta = (/'x_eta_x', 'x_eta_y', 'x_eta_z'/)
  r_zeta = (/'x_zeta_x', 'x_zeta_y', 'x_zeta_z'/)
  r_ksi_ksi = (/'x_ksi_ksi_x', 'x_ksi_ksi_y', 'x_ksi_ksi_z'/)
  r_eta_eta = (/'x_eta_eta_x', 'x_eta_eta_y', 'x_eta_eta_z'/)
  diss = (/'diss_x', 'diss_y', 'diss_z'/)

  ! Open the CGNS File:
  if (myid == 0) then
     call cg_open_f(fileName, CG_MODE_WRITE, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
     ! Write the single base
     cellDim = 3
     physDim = 3
     call cg_base_write_f(cg,"Base#1", Celldim, Physdim, base, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
  end if
  if (writeMirror) then
     nPatchToWrite = nPatch
  else
     nPatchToWrite = nPatch/2
  end if
 
  do iPatch=1, nPatchToWrite
     if (myid == 0) then
        ! Write the single zone
999     FORMAT('domain.',I5.5)
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
     end if

     do idim=1,3
        ! Copy values and write:
        do k=1,N

           call VecScatterBegin(rootScatter, X(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
           call EChk(ierr,__FILE__,__LINE__)
           call VecScatterEnd(rootScatter, X(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           if (myid == 0) then
              call VecGetArrayF90(XLocal, xx, ierr)
              call EChk(ierr, __FILE__, __LINE__)

              do j=1,patches(iPatch)%jl
                 do i=1,patches(iPatch)%il
                    idGlobal = patches(iPatch)%l_index(i,j)
                    coordArray(i,j,k) = xx(3*(idGlobal-1)+iDim)
                 end do
              end do
              call VecRestoreArrayF90(XLocal, xx, ierr)
              call EChk(ierr, __FILE__, __LINE__)
           end if
        end do

        if (myid == 0) then
           call cg_coord_write_f(cg, base, zone, realDouble, &
                coorNames(iDim), coordArray, coordID, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end if
     end do
     if (writeMetrics .and. metricsAllocated) then
        if (myid == 0) then
           call cg_sol_write_f(cg, base, zone, "FlowSolution", Vertex, s, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end if

        call writeVar(metrics(:, iX_ksi), r_ksi)
        call writeVar(metrics(:, iX_eta), r_eta)
        call writeVar(metrics(:, iX_zeta), r_zeta)
        call writeVar(metrics(:, iX_ksi_ksi), r_ksi_ksi)
        call writeVar(metrics(:, iX_eta_eta), r_eta_eta)
        call writeVar(metrics(:, iX_diss), diss)
        
     !    ! ! Volume is special since there is only 1
     !    ! do k=1,N
     !    !    call VecGetArrayF90(Vhist(k), xxtmp, ierr)
     !    !    call EChk(ierr, __FILE__, __LINE__)

     !    !    do j=1,patches(iPatch)%jl
     !    !       do i=1,patches(iPatch)%il
     !    !          idGlobal = patches(iPatch)%l_index(i,j)
     !    !          solArray(i,j,k) = xxtmp(3*idGlobal)
     !    !       end do
     !    !    end do
     !    !    call VecRestoreArrayF90(Vhist(k), xxtmp, ierr)
     !    !    call EChk(ierr, __FILE__, __LINE__)
     !    ! end do
     !    ! call cg_field_write_f(cg, base, zone, s, RealDouble, &
     !    !      'volume', solArray, F, ierr)
     !    ! if (ierr .eq. CG_ERROR) call cg_error_exit_f
     end if

     ! Deallocate the coord and solution array for this patch
     if (myid == 0) then
        deallocate(coordArray, solArray)
     end if
  end do ! Patch Loop

  if (myid == 0) then
     ! Finally close cgns_file
     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
  end if

  contains
    subroutine writeVar(vectors, names)
      Vec, dimension(:) :: vectors
      character*12, dimension(3) :: names
      
      do idim=1,3
         do k=1,N

            call VecScatterBegin(rootScatter, vectors(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call EChk(ierr,__FILE__,__LINE__)
            call VecScatterEnd(rootScatter, vectors(k), XLocal, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call EChk(ierr,__FILE__,__LINE__)

            if (myid == 0) then
               call VecGetArrayF90(XLocal, xx, ierr)
               call EChk(ierr, __FILE__, __LINE__)

               do j=1,patches(iPatch)%jl
                  do i=1,patches(iPatch)%il
                     idGlobal = patches(iPatch)%l_index(i,j)
                     solArray(i,j,k) = xx(3*(idGlobal-1)+iDim)
                  end do
               end do

               call VecRestoreArrayF90(XLocal, xx, ierr)
               call EChk(ierr, __FILE__, __LINE__)
            end if
         end do
         if (myid == 0) then
            call cg_field_write_f(cg, base, zone, s, RealDouble, &
                 names(idim), solArray, F, ierr)
            if (ierr .eq. CG_ERROR) call cg_error_exit_f
         end if
      end do
    end subroutine writeVar
  end subroutine writeCGNS

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
  !     Parameters
  !     ----------
  !     fileName : Character array
  !        The name of the cgns file to modify
  !     mirrorList : integer array, size(nMirror, 3)
  !        The list of faces and node indices (on surface) to mirror
  !     mirrorDim : int 
  !        Dimension to mirror, 1 for x, 2 for y, 3 for z.
  !     nMirror : integer
  !        The number of faces to modify.
  !
  !     Notes
  !     -----
  !     Since this uses file-i/o it should only be called on the root
  !     processor

  use hypInput
  use hypData
  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*), intent(in) :: fileName
  integer(kind=intType), intent(in) :: mirrorList(nMirror, 3)
  integer(kind=intType), intent(in) :: mirrorDim, nMirror

  ! Working Variables
  integer(kind=intType) :: cg, ierr, base, nzones, ii, zone, zonesize(9)
  integer(kind=intType) :: coordID, start(3), zoneType, loadedZone 
  character*32 zonename
  real(kind=realType), dimension(:, :, :), allocatable :: coor

  ! Don't do anything if no mirror nodes
  if (nMirror == 0) then
     return
  end if

  ! Open the CGNS File:
  call cg_open_f(fileName, CG_MODE_MODIFY, cg, ierr)
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

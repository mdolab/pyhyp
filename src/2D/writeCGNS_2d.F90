subroutine writeCGNS_2d(fileName)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeCGNS_2d write the current grid to a 3D CGNS
  !               file that is exactly 1 cell of unit width. It also
  !               applies default boundary and symmetry conditions
  !               such that the grid can be directly used in a 3D flow
  !               solver such as sumb.
  !
  !     Description of Arguments
  !     Input:
  !     fileNmae - Character array: the name of the cgns file
  !
  !     Ouput: None

  use hypInput
  use hypData, only : grid2D

  implicit none
  include 'cgnslib_f.h'
  
  ! Input Arguments
  character*(*) :: fileName
  
  ! Working Variables
  integer(kind=intType) :: cg, ierr
  integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
  character*256 :: zoneName
  integer(kind=intType) :: sizes(9), zone, ii, i, j, k, nx, transform(3)
  integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
  real(kind=realType), dimension(:), allocatable :: coordArray

  ! Open the CGNS File:
  call cg_open_f(fileName, MODE_WRITE, cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Write the single base
  cellDim = 3
  physDim = 3
  call cg_base_write_f(cg,"Base#1", Celldim, Physdim, base, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Write the single zone
999 FORMAT('domain.',I5.5)
  write(zonename,999) 1 ! Domains are typically ordered form 1
  gridShp = shape(grid2D)
  nx = gridShp(2)
  sizes(:) = 0_intType
  sizes(1) = nx + 1
  sizes(2) = N 
  sizes(3) = 2
  sizes(4) = sizes(1) - 1
  sizes(5) = sizes(2) - 1
  sizes(6) = sizes(3) - 1 
  call cg_zone_write_f(cg, base, zonename, sizes, Structured, &
       zone, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  allocate(coordArray(sizes(1)*sizes(2)*sizes(3)))

  ! Assemble the x-coordinates
  ii = 0
  do k=1, 2
     do j=1, N
        do i=1, nx
           ii = ii + 1
           coordArray(ii) = grid2D(1, i, j)
        end do
        ii = ii + 1
        coordArray(ii) = grid2D(1, 1, j)
     end do
  end do

  call cg_coord_write_f(cg, base, zone, realDouble, &
       'CoordinateX', coordArray, coordID, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Assemble the y-coordinates
  ii = 0
  do k=1, 2
     do j=1, N
        do i=1, nx
           ii = ii + 1
           coordArray(ii) = grid2D(2, i, j)
        end do
        ii = ii + 1
        coordArray(ii) = grid2D(2, 1, j)
     end do
  end do

  call cg_coord_write_f(cg, base, zone, realDouble, &
       'CoordinateY', coordArray, coordID, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  
  ! Assemble the z-coordinates
  ii = 0
  do k=1, 2
     do j=1, N
        do i=1, nx+1
           ii = ii + 1
           if (k == 1) then
              coordArray(ii) = zero
           else
              coordArray(ii) = one
           end if
        end do
     end do
  end do

  call cg_coord_write_f(cg, base, zone, realDouble, &
       'CoordinateZ', coordArray, coordID, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f
  
  ! Deallocate the coordiante array:
  deallocate(coordArray)

  ! That takes care of the coordinates, now we have to deal with
  ! boundary conditions. These are straight forward: The j-low is the
  ! given surface, which we will set as BCWallViscous. The j-high is a
  ! FarField condtion. The k-low and k-high surfaces are mirror planes. 

  ! ---------- Jlow --------------
  pnts(:, 1) = (/1, 1, 1/)
  pnts(:, 2) = (/sizes(1), 1, sizes(3)/)
  
  call cg_boco_write_f(cg, base, zone, "jlow", BCWallViscous, PointRange, &
       2, pnts, BCout, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg,base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
       "BC_t", BCOut, "end")
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f("wall", ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! ---------- Jhigh --------------
  pnts(:, 1) = (/1, sizes(2), 1/)
  pnts(:, 2) = (/sizes(1), sizes(2), sizes(3)/)
  
  call cg_boco_write_f(cg, base, zone, "jhigh", BCFarField, PointRange, &
       2, pnts, BCout, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg, base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
       "BC_t", BCOut, "end")
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f("far", ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! ---------- Klow --------------
  pnts(:, 1) = (/1, 1, 1/)
  pnts(:, 2) = (/sizes(1), sizes(2), 1/)
  
  call cg_boco_write_f(cg, base, zone, "klow", BCSymmetryPlane, PointRange, &
       2, pnts, BCout, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg, base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
       "BC_t", BCOut, "end")
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f("sym", ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! ---------- Khigh --------------
  pnts(:, 1) = (/1, 1, sizes(3)/)
  pnts(:, 2) = (/sizes(1), sizes(2), sizes(3)/)
  
  call cg_boco_write_f(cg, base, zone, "khigh", BCSymmetryPlane, PointRange, &
       2, pnts, BCout, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_goto_f(cg, base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
       "BC_t", BCOut, "end")
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  call cg_famname_write_f("sym", ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! That takes care of the boundary conditions. The last thing we have
  ! to deal with is the grid connectivity. This is a little bit
  ! tricker since we need to get the periodic block to block
  ! connection correct at the period boundary. This takes take of the
  ! "boundary" condition on the iLow and iHigh faces.

  ! ---------- Ilow --------------
  pnts(:, 1) = (/1, 1, 1/)
  pnts(:, 2) = (/1, sizes(2), sizes(3)/)

  pnts_donor(:, 1) = (/sizes(1), 1, 1 /)
  pnts_donor(:, 2) = (/sizes(1), sizes(2), sizes(3) /)
  
  transform(1) = 1
  transform(2) = 2
  transform(3) = 3

  call cg_1to1_write_f(cg, base, zone, "b2b1", "domain.00001", pnts, pnts_donor, &
       transform, nCon, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! ---------- Ihigh --------------
  pnts(:, 1) = (/sizes(1), 1, 1/)
  pnts(:, 2) = (/sizes(1), sizes(2), sizes(3)/)

  pnts_donor(:, 1) = (/1, 1, 1 /)
  pnts_donor(:, 2) = (/1, sizes(2), sizes(3) /)
  
  transform(1) = 1
  transform(2) = 2
  transform(3) = 3

  call cg_1to1_write_f(cg, base, zone, "b2b2", "domain.00001", pnts, pnts_donor, &
       transform, nCon, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Finally close cgns_file
  call cg_close_f(cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeCGNS_2d

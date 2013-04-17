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
  !     fileNmae - Character array: the name of the cgns file
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
  integer(kind=intType) :: sizes(9), zone, ii, i, j, k, nx, transform(3)
  integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
  real(kind=realType), dimension(:,:,:), allocatable :: coordArray
  character*12, dimension(3) :: coorNames
  integer(kind=intType) :: iPatch, idim, idglobal
  ! Open the CGNS File:
  call cg_open_f(fileName, MODE_WRITE, cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! Write the single base
  cellDim = 3
  physDim = 3
  call cg_base_write_f(cg,"Base#1", Celldim, Physdim, base, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

  do iPatch=1, nPatch
     ! Write the single zone
999  FORMAT('domain.',I5.5)
     write(zonename,999) iPatch ! Domains are typically ordered form 1
     sizes(:) = 0_intType
     sizes(1) = patches(iPatch)%il
     sizes(2) = patches(iPatch)%jl
     sizes(3) = NDebug
     sizes(4) = sizes(1) - 1
     sizes(5) = sizes(2) - 1
     sizes(6) = sizes(3) - 1 
     call cg_zone_write_f(cg, base, zonename, sizes, Structured, &
          zone, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
     ! Allocate the coord array
     allocate(coordArray(sizes(1),sizes(2),sizes(3)))
     coorNames(1) = "CoordinateX"
     coorNames(2) = "CoordinateY"
     coorNames(3) = "CoordinateZ"     

     do idim=1,3
        ! Copy values and write:
        do k=1,NDebug
           do j=1,patches(iPatch)%jl
              do i=1,patches(iPatch)%il
                 idGlobal = patches(iPatch)%l_index(i,j)
                 coordArray(i,j,k) = grid3d(idim, idGlobal, k)
              end do
           end do
        end do
        call cg_coord_write_f(cg, base, zone, realDouble, &
             coorNames(iDim), coordArray, coordID, ierr)
     end do
  
     ! Deallocate the coord array for this patch
     deallocate(coordArray)

     ! We know we can write the klow wall BC and the k-high farfield
     ! ---------- Klow --------------
     pnts(:, 1) = (/1, 1, 1/)
     pnts(:, 2) = (/sizes(1), sizes(2), 1/)
  
     call cg_boco_write_f(cg, base, zone, "klow", BCWallViscous, PointRange, &
          2, pnts, BCout, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_goto_f(cg,base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
          "BC_t", BCOut, "end")
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_famname_write_f("wall", ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     ! ---------- Khigh --------------
     pnts(:, 1) = (/1, 1, sizes(3)/)
     pnts(:, 2) = (/sizes(1), sizes(2), sizes(3)/)
  
     call cg_boco_write_f(cg, base, zone, "khigh", BCFarField, PointRange, &
          2, pnts, BCout, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_goto_f(cg,base,ierr,'Zone_t', zone,"ZoneBC_t", 1,&
          "BC_t", BCOut, "end")
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_famname_write_f("far", ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

  end do ! Patch Loop

  ! Finally close cgns_file
  call cg_close_f(cg, ierr)
  if (ierr .eq. CG_ERROR) call cg_error_exit_f

end subroutine writeCGNS_3D

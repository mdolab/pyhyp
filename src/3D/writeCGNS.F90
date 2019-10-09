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
  use hypData, only : X, metrics, rootScatter, xLocal, metricsAllocated, iVHist
  use hypData, only : nPatch, patches, xx, metrics, ix_ksi, ix_eta, ix_zeta, iX_eta_eta, iX_ksi_ksi, iX_diss
  use cgnsGrid

#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

  ! Input Arguments
  character*(*) :: fileName

  ! Working Variables
  integer(kind=intType) :: cg, ierr, iEdge, il, jl
  integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
  character*256 :: zoneName, bcName
  integer(kind=intType) :: zone, ii, i, j, k, cgtransform(3)
  integer(cgsize_t) :: sizes(9), ptRange(3,2)
  integer(kind=intType) ::  fullRange(3,2), iBoco
  integer(kind=intType) :: s, f
  real(kind=realType), dimension(:,:,:), allocatable :: coordArray, solArray
  character*12, dimension(3) :: coorNames, r_ksi, r_eta, r_zeta
  character*12, dimension(3) :: r_ksi_ksi, r_eta_eta, diss
  integer(kind=intType) :: iPatch, idim, idglobal
  real(kind=realType), pointer, dimension(:) :: xxtmp

  iBoco = 0

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

  do iPatch=1, nPatch
     if (myid == 0) then
        il = patches(iPatch)%il
        jl = patches(iPatch)%jl
        fullRange = reshape((/1, 1, 1, il, jl, N/), (/3,2/))

        ! Write the single zone
999     FORMAT('domain.',I5.5)
        write(zonename,999) iPatch ! Domains are typically ordered form 1
        sizes(:) = 0
        sizes(1) = il
        sizes(2) = jl
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

              do j=1, jl
                 do i=1, il
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

        ! Volume is special since there is only 1
        do k=1, N
           call VecGetArrayF90(metrics(k, iVHist), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           do j=1, jl
              do i=1, il
                 idGlobal = patches(iPatch)%l_index(i,j)
                 solArray(i,j,k) = xxtmp(3*idGlobal)
              end do
           end do
           call VecRestoreArrayF90(metrics(k, iVHist), xxtmp, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end do
        if (myid == 0) then
           call cg_field_write_f(cg, base, zone, s, RealDouble, &
                'volume', solArray, F, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f
        end if
     end if

     ! Deallocate the coord and solution array for this patch
     if (myid == 0) then
        deallocate(coordArray, solArray)
     end if

     ! ------- Write the BC for K=1 plane -------------
     ptRange = fullRange
     ptRange(3, :) = 1
     call writeBC(bcWallViscous, ptRange, families(iPatch))

     ! ------- Write the BC for K=N plane -------------
     ptRange(3, :) = N
     if (outerfacetype == outerFaceFarfield) then
        call writeBC(bcFarField, ptRange, "Far")
     else if(outerFaceType == outerFaceOverset) then
        call writeBC(CG_UserDefined, ptRange, "Overset")
     end if

     ! ------- Write Possible BC for other 4 sides ---------------
     do iEdge = 1,4
        ptRange = fullRange
        select case(iEdge)
        case (iLow)
           ptRange(1, :) = 1
        case (iHigh)
           ptRange(1, :) = il
        case (jLow)
           ptRange(2, :) = 1
        case (jHigh)
           ptRange(2, :) = jl
        end select

        if (BCs(iEdge, iPatch) == BCSplay .or. &
           BCs(iEdge, iPatch) == BCXConst .or. &
           BCs(iEdge, iPatch) == BCYConst .or. &
           BCs(iEdge, iPatch) == BCZConst) then
           call writeBC(CG_UserDefined, ptRange, "Overset")
        else if (BCs(iEdge, iPatch) == BCXSymm .or. &
             BCs(iEdge, iPatch) == BCYSymm .or. &
             BCs(iEdge, iPatch) == BCZSymm) then
           call writeBC(BCSymmetryPlane, ptRange, "Sym")
        end if
     end do
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

  subroutine writeBC(bcType, ptRange, famName)
    ! Helper routine for writing a boundary condition
    use hypInput
    implicit none
    integer(kind=intType) :: bcType, BCOut
    integer(kind=cgsize_t) :: ptRange(3,2)
    character*(*) :: famName
    if (myid == 0) then
998    FORMAT('boco.',I5.5)
       iBoco = iBoco + 1
       write (bcName, 998), iBoco

       call cg_boco_write_f(cg, base, iPatch, trim(bcName), &
            BCType, PointRange, 2_cgsize_t, ptRange, BCout, ierr)
       if (ierr == CG_ERROR) call cg_error_exit_f

       call cg_goto_f(cg, base, ierr, 'Zone_t', iPatch, "ZoneBC_t", 1, &
            "BC_t", BCOut, "end")
       if (ierr .eq. CG_ERROR) call cg_error_exit_f

       call cg_famname_write_f(trim(famName), ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f

       !Add an user-defined data node for the overset BC case
       if (bcType == CG_UserDefined) then
          call cg_user_data_write_f('BCOverset', ierr)
          if (ierr .eq. CG_ERROR) call cg_error_exit_f
       end if
    end if
  end subroutine writeBC
end subroutine writeCGNS

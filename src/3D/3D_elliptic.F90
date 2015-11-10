subroutine FormFunction_mf(ctx, stateVec, resVec, ierr)

  use precision
  use communication
  use hypData, only : ellipScatter, localSol, xx
  use panel
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! PETSc Variables
  PetscFortranAddr ctx(*)
  Vec     stateVec, resVec
  integer(kind=intType) :: ierr, iLow, iHigh, nPLocal, j
  real(kind=realType) :: V(3)

  call setStrengths(stateVec)

  ! Now we can pull out a local vector for resVec since we will only
  ! be setting local variables. 
  call VecGetArrayF90(resVec, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecGetOwnershipRange(stateVec, iLow, iHigh, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  nPLocal = iHigh-iLow

  do j=1, nPLocal
     call evalAtPoint(MGP(1)%panels(j + iLow)%center, xx(j), V)
  end do

  call VecRestoreArrayF90(resVec, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! We don't check an error here, so just pass back zero
  ierr = 0

end subroutine FormFunction_mf

subroutine runElliptic
  ! run3DElliptic is the main python interface for generatnig 3D elliptic meshes. 
  !
  !
  ! Notes
  ! -----
  ! The init3d routine must have already been called with the
  ! addtional setup information. 
  use communication
  use hypData
  use hypInput
  use panel
  implicit none
#include "cgnslib_f.h"
  ! Working parameters
  integer(kind=intType) :: i, j, ierr
  real(kind=realType), dimension(100) :: targets
  real(kind=realType) ::  d, endPt(3), fact, startPt(3), dx(3)
  integer(kind=intType) :: nTarget, neval
  integer(kind=intType) :: istart, iend, inds(3)

  integer(kind=intType) ::  cg, nzones, zone, zoneSize(9), s, F, base, k
  character*32 :: zonename
  logical :: successful
  real(kind=realType), dimension(:, :, :), allocatable  :: coorX, coorY, coorZ, valPhi, vX, vY, vZ
  real(kind=realType) :: pt(3), V(3), phi, times(2)
  ! First we need to check if we can use the the source strength file.
  ! If so, it would mean we don't have to solve for the strengths
  call loadSolution(successful)

  if (.not. successful) then
     times(1) = mpi_wtime()
     call generateSolution
     times(2) = mpi_wtime()
     if (myid == 0) then
        print *,'Solution time:', times(2)-times(1)
     end if
  end if
  
  !  if (.not. myid ==0) then
  !    return 
  ! end if
  ! call cg_open_f("L1_with_potential.cgns", CG_MODE_MODIFY, cg, ierr)
  ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! base = 1
  ! ! Goto Base Node
  ! call cg_goto_f(cg, base, ierr, 'end')
  ! if (ierr .eq. CG_ERROR) call cg_error_exit_f
  ! call cg_nzones_f(cg, base, nzones, ierr)
  ! if (ierr .eq. CG_ERROR) call cg_error_exit_f

  ! do zone=1, nzones
  !    call cg_zone_read_f(cg, base, zone, zonename, zonesize, ierr)
  !    print *,'zone:',zone, zonesize(1:3)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f
  !    allocate(coorX(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(coorY(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(coorZ(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(valPhi(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(vX(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(vY(zoneSize(1), zoneSize(2), zoneSize(3)))
  !    allocate(vZ(zoneSize(1), zoneSize(2), zoneSize(3)))

  !    call cg_coord_read_f(cg, base, zone, 'CoordinateX', RealDouble, (/1,1,1/), &
  !         zonesize, coorX, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    call cg_coord_read_f(cg, base, zone, 'CoordinateY', RealDouble, (/1,1,1/), &
  !         zonesize, coorY, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
  !    call cg_coord_read_f(cg, base, zone, 'CoordinateZ', RealDouble, (/1,1,1/), &
  !         zonesize, coorZ, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f
     
  !    do k=1,zonesize(3)!-1
  !       do j=1,zonesize(2)!-1
  !          do i=1,zonesize(1)!-1
  !             ! pt(1) = one/eight * (&
  !             !      coorX(i,j,k)+coorX(i+1,j,k)+coorX(i,j+1,k)+coorX(j+1,j+1,k)+&
  !             !      coorX(i,j,k+1)+coorX(i+1,j,k+1)+coorX(i,j+1,k+1)+coorX(j+1,j+1,k+1))

  !             ! pt(2) = one/eight * (&
  !             !      coorY(i,j,k)+coorY(i+1,j,k)+coorY(i,j+1,k)+coorY(j+1,j+1,k)+&
  !             !      coorY(i,j,k+1)+coorY(i+1,j,k+1)+coorY(i,j+1,k+1)+coorY(j+1,j+1,k+1))

  !             ! pt(3) = one/eight * (&
  !             !      coorZ(i,j,k)+coorZ(i+1,j,k)+coorZ(i,j+1,k)+coorZ(j+1,j+1,k)+&
  !             !      coorZ(i,j,k+1)+coorZ(i+1,j,k+1)+coorZ(i,j+1,k+1)+coorZ(j+1,j+1,k+1))
  !             pt(1) = coorX(i,j,k)
  !             pt(2) = coorY(i,j,k)
  !             pt(3) = coorZ(i,j,k)
  !             call evalAtPoint(pt, phi, V)
  !             valPhi(i,j,k) = phi
  !             Vx(i,j,k) = V(1)
  !             Vy(i,j,k) = V(2)
  !             Vz(i,j,k) = V(3)
  !          end do
  !       end do
  !    end do

  !    ! Now write the field
  !    call cg_sol_write_f(cg, base, zone, "FlowSolution", Vertex, s, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    call cg_field_write_f(cg, base, zone, s, RealDouble, &
  !         "phi", valPhi, F, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    call cg_field_write_f(cg, base, zone, s, RealDouble, &
  !         "Vx", Vx, F, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    call cg_field_write_f(cg, base, zone, s, RealDouble, &
  !         "Vy", Vy, F, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    call cg_field_write_f(cg, base, zone, s, RealDouble, &
  !         "Vz", Vz, F, ierr)
  !    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  !    deallocate(coorX, coorY, coorZ, valPhi, Vx, Vy, Vz)
  ! end do
  ! call cg_close_f(cg, ierr)
  ! if (ierr .eq. CG_ERROR) call cg_error_exit_f
  ! return

  ! Generate some dummy target data
  targets(1:49) = (/&
       0.99910944810, &
       0.99909686620, &
       0.99908468440, &
       0.99907097960, &
       0.99904947960, &
       0.99901576160, &
       0.99896290830, &
       0.99888011560, &
       0.99875053840, &
       0.99854796260, &
       0.99823167790, &
       0.99773861560, &
       0.99697139470, &
       0.99578048330, &
       0.99393837930, &
       0.99109986910, &
       0.98677521000, &
       0.98026199120, &
       0.97065198760, &
       0.95686742290, &
       0.93810733410, &
       0.91401221280, &
       0.88540284820, &
       0.85436210400, &
       0.82390724310, &
       0.79557388480, &
       0.76284935420, &
       0.72584510730, &
       0.68535151250, &
       0.64241604440, &
       0.59802614290, &
       0.55209454580, &
       0.50610900630, &
       0.45923529230, &
       0.41301878710, &
       0.36727009730, &
       0.32269817870, &
       0.28017928560, &
       0.24047220010, &
       0.20295351010, &
       0.16783459550, &
       0.13767111480, &
       0.11036346190, &
       0.08609719161, &
       0.06570453572, &
       0.04894374543, &
       0.03549951477, &
       0.02500139869, &
       0.01712893882  &
       /)
  ! targets(1:23) = (/&
  !      0.99909686620, &
  !      0.99907097960, &
  !      0.99901576160, &
  !      0.99888011560, &
  !      0.99854796260, &
  !      0.99773861560, &
  !      0.99578048330, &
  !      0.99109986910, &
  !      0.98026199120, &
  !      0.95686742290, &
  !      0.91401221280, &
  !      0.82390724310, &
  !      0.76284935420, &
  !      0.68535151250, &
  !      0.59802614290, &
  !      0.50610900630, &
  !      0.36727009730, &
  !      0.28017928560, &
  !      0.20295351010, &
  !      0.13767111480, &
  !      0.08609719161, &
  !      0.04894374543, &
  !      0.02500139869 &
  !      /)


  targets = targets + .00089
!  targets(2) = .95

  ! targets(1) = 1.00
  ! targets(2) = 0.999999
  ! do i=3,N ! N is the number of levels
  !    d = targets(i-2)-targets(i-1)
  !    targets(i) = targets(i-1) - d*1.15
  ! end do
  if (myid == 0) then
     print *, targets(1:N)
  end if

  call VecCreate(hyp_comm_world, ellipDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecSetBlockSize(ellipDelta, 3, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecSetType(ellipDelta, VECMPI, ierr)
  call VecSetSizes(ellipDelta, PETSC_DETERMINE, 3*npGlobal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(ellipDelta, ellipNorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecDuplicate(ellipDelta, ellipLayer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterCreateToAll(ellipDelta, ellipDeltaScatter, ellipDeltaLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecDuplicate(ellipDeltaLocal, ellipNormLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecDuplicate(ellipDeltaLocal, ellipLayerLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Get the range of the start/stop of ellipDelta. This determines the
  ! range of panels we need to loop over in the loop below:
  call VecGetOwnershipRange(ellipDelta, iStart, iEnd, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Eliminate the effect of the block size and convert to fortran
  ! ordering:
  iStart = iStart/3 + 1
  iEnd = iEnd/3
  do j=2, N  ! Nominal march loop
     if (myid == 0) then
        print *,'Percent:',dble(j-1)/N*100
     end if
     
     ! Generate the panel centers of the j-1 layer and store into
     ! ellipLayer. 
     call generateEllipLayer(X(j-1), ellipLayer)

     ! Distributed loop over the owned panels
     neval = 0
     do i=iStart, iEnd

        ! Note that getValues is OK here since we are only getting
        ! values on-processor
        inds = (/3*(i-1), 3*(i-1)+1, 3*(i-1)+2/)
        call VecGetValues(ellipLayer, 3, inds, startPt, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        ! Actual call to advance one point to the target
        call marchPt(startPt, targets(j), endPt, V, neval)

        ! Dump the final values back into the delta and norm arrys
        dx = endPt - startPt
        call vecSetValues(ellipDelta, 3, inds, dx, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call vecSetValues(ellipNorm, 3, inds, V, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
     print *,'finished on proc:', myid, neval, iend-istart

     ! Need to finish assembly of the values we added to ellipDelta
     ! and ellipNorm above (even though we have not added any
     ! off-processor values)
     call VecAssemblyBegin(ellipDelta, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecAssemblyBegin(ellipNorm, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecAssemblyEnd(ellipDelta, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecAssemblyEnd(ellipNorm, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! For simplicitiy here, communication the ellipDelta and
     ! ellipNorm to all procs such that everyone has a full copy of
     ! the values
     call VecScatterBegin(ellipDeltaScatter, ellipDelta, ellipDeltaLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecScatterEnd(ellipDeltaScatter, ellipDelta, ellipDeltaLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecScatterBegin(ellipDeltaScatter, ellipNorm, ellipNormLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecScatterEnd(ellipDeltaScatter, ellipNorm, ellipNormLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecScatterBegin(ellipDeltaScatter, ellipLayer, ellipLayerLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecScatterEnd(ellipDeltaScatter, ellipLayer, ellipLayerLocal, &
          INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr, __FILE__, __LINE__)


     ! Call the reconstruct algorithm to get the nodes from the
     ! centers/normals
     ! if (j == 2) then 
     !    call reconstruct(X(j), X(j-1), ellipLayerLocal, ellipDeltaLocal, ellipNormLocal)
     ! else
     !    call reconstruct2(X(j), X(j-1), ellipLayerLocal, ellipDeltaLocal, ellipNormLocal)
     ! end if
     call reconstruct2(X(j), X(j-1), ellipLayerLocal, ellipDeltaLocal, ellipNormLocal)
     ! Now surface smooth the current level if necessary
     call surfaceSmooth2(X(j), 25,  one-targets(j))
  end do

  ! Destroy all the temporary PETSc variables
  call VecDestroy(ellipDelta, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecDestroy(ellipNorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(ellipDeltaLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecDestroy(ellipNormLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(ellipLayer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(ellipLayerLocal, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterDestroy(ellipDeltaScatter, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine runElliptic

subroutine generateSolution
  ! generateSolution determines the source strength on the panels

  use communication
  use hypData
  use hypInput
  use panel
  implicit none
#include "include/petscversion.h"
  ! Working parameters
  integer(kind=intType) :: i, j, k, info, ierr, iLow, iHigh, nPLocal
  integer(kind=intType), dimension(:), allocatable :: onProc, offProc
  real(kind=realType) :: d, phi, V(3), dist, val
  logical :: assembleMat, successful
  external formFunction_mf

  if ( (evalMode == EVAL_SLOW .or. evalMode == EVAL_EXACT) .and. .not. useMatrixFree) then
     ! We have a dense matrix:
     call MatCreateDense(hyp_comm_world, PETSC_DETERMINE, PETSC_DETERMINE, &
          nPGlobal, nPGlobal, PETSC_NULL_SCALAR, ellipMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#if PETSC_VERSION_MINOR > 5
     call MatCreateVecs(hypMat, hypDelta, hypRHS, ierr)
#else
   call MatGetVecs(hypMat, hypDelta, hypRHS, ierr)
#endif

     call EChk(ierr, __FILE__, __LINE__)
     assembleMat = .True.
  else
     ! Otherwise it will be matrix-free matrix
     call MatCreateMFFD(hyp_comm_world, PETSC_DETERMINE, PETSC_DETERMINE, &
          nPGlobal, nPGlobal, ellipMat, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MatMFFDSetFunction(ellipMat,  formFunction_mf, ctx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#if PETSC_VERSION_MINOR > 5
     call MatCreateVecs(hypMat, hypDelta, hypRHS, ierr)
#else
   call MatGetVecs(hypMat, hypDelta, hypRHS, ierr)
#endif
     call EChk(ierr, __FILE__, __LINE__)

     call MatMFFDSetBase(ellipMat, ellipSol, PETSC_NULL_OBJECT, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     assembleMat = .False.
  end if

  call VecGetOwnershipRange(ellipRHS, iLow, iHigh, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  nPLocal = iHigh-iLow
  if (myid == 0) then
     print *,'Determining PC size...'
  end if

  ! Create the required scatter context
  call VecScatterCreateToAll(ellipSol, ellipScatter, localSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  allocate(onProc(nPLocal), offProc(nPLocal), stat=ierr)
  onProc(:) = 0
  offProc(:) = 0

  do j=1,nPLocal
     do i=1, nPGlobal
        d = dist(MGP(1)%panels(i)%center, MGP(1)%panels(j+iLow)%center)
        if (d < MGP(1)%panels(i)%length) then
           ! Determine if on-proc or off-proc
           if (i>=iLow .and. i < iHigh) then
              onProc(j) = onProc(j) + 1
           else
              offProc(j) = offProc(j) + 1
           end if
        end if
     end do
  end do

  ! We arbitrarily add an extra 10 rows in case PETSc wants to be
  ! dumb. Also we make sure that the number of non-zeros isn't more
  ! than the number of columns. 
  onProc = min(onProc + 10, nplocal)
  offProc = min(offProc + 10, npglobal-nplocal)

  call MatCreateAIJ(hyp_comm_world, nPLocal, nPLocal, PETSC_DETERMINE, PETSC_DETERMINE,&
       0, onProc, 0, offProc, ellipPCMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  deallocate(onProc, offProc)

  call KSPCreate(petsc_comm_world, ellipKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetFromOptions(ellipKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGMRESSetRestart(ellipKSP, kspSubspacesize, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetOperators(ellipKSP, ellipMat, ellipPCMat, ierr)

  call KSPSetTolerances(ellipKSP, kspRelTol, 1e-16, 1e20, kspMaxIts, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPMonitorSet(ellipKSP, kspmonitordefault, PETSC_NULL_OBJECT, &
       PETSC_NULL_FUNCTION, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (myid == 0) then
     print *,'Assembling Matrix ...'
  end if
  do j=1,npLocal
     do i=1, nPGlobal
        ! Only evalaute if assembling matrix
        if (assembleMat) then
           call panelInfluence2(i, MGP(1)%panels(j+iLow)%center, phi, V)
           call MatSetValues(ellipMat, 1, j + iLow - 1, 1, i-1, phi, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if

        d = dist(MGP(1)%panels(i)%center, MGP(1)%panels(j+iLow)%center)

        if (d < MGP(1)%panels(i)%length) then
           ! If we're JUST assembling the PC we need to actually do the evaluation here
           if (.not. assembleMat) then
              call panelInfluence2(i, MGP(1)%panels(j+iLow)%center, phi, V)
           end if

           call MatSetValues(ellipPCMat, 1, j + iLow - 1, 1, i-1, phi, INSERT_VALUES, ierr)
           call EChk(ierr, __FILE__, __LINE__)
        end if
     end do
  end do

  ! Now finish assembly of the Matrix and the preconditioner
  call MatAssemblyBegin(ellipMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyBegin(ellipPCMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(ellipMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(ellipPCMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! This vec set is where we specify that the potential on the surface is exactly 1.0
  call VecSet(ellipRhs, one, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Apparently LEFT preconditioning doesn't work AT ALL for this
  ! problem. So hard-code Right preconditioning. 
  call KSPGetPC(ellipKSP, pc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPsetPCSide(ellipksp, PC_RIGHT, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Solve and use the scatter context to give the solution to everyone. 
  call KSPSolve(ellipKSP, ellipRHS, ellipSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the solution
  call setStrengths(ellipSol)

  ! Save the solution 
  call saveSolution(ellipSol)

  if (.not. assembleMat) then
     call MatMFFDSetBase(ellipMat, ellipSol, PETSC_NULL_OBJECT, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  call MatMult(ellipMat, ellipSol, ellipRHS, ierr)
  call VecSet(ellipSol, -one, ierr)
  call VecAXPY(ellipRHS, one, ellipSol, ierr)
  call VecNorm(ellipRHS, NORM_2, val, ierr)
  if (myid == 0) then
     print *,'norm:',val
  end if

  ! Clean up all the PETSc data:
  call KSPDestroy(ellipKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatDestroy(ellipMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatDestroy(ellipPCMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(ellipSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(ellipRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDestroy(localSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterDestroy(ellipScatter, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine generateSolution

subroutine setStrengths(sigma)

  ! setStrengths takes the gloabl PETSc vector and updates the source
  ! strengths in the multi-grid structure. 
  !  
  ! Parameters
  ! ----------
  ! sigma : PETSc Vec (parallel)
  !    The distributed source strengths. 

  use precision
  use hypData, only : ellipScatter, localSol, xx
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input Parameters
  Vec :: sigma

  ! Working Parameters
  integer(kind=intType) :: ierr

  ! First we have to take the states and scatter them to everyone and
  ! set them. This communication is the only communication for the
  ! entire matrix vector product. 

  call VecScatterBegin(ellipScatter, sigma, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd(ellipScatter, sigma, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Extract a pointer and set the strengths
  call VecGetArrayF90(localSol, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now set the strengths
  call setMGData(xx)

  ! Always remember to restore the array
  call VecRestoreArrayF90(localSol, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine setStrengths

subroutine saveSolution(sigma)
  ! setStrengths takes the gloabl PETSc vector and updates the source
  ! strengths in the multi-grid structure. 
  !  
  ! Parameters
  ! ----------
  ! sigma : PETSc Vec (parallel)
  !    The distributed source strengths. 


  ! First we have to take the states and scatter them to everyone and
  ! set them. This communication is the only communication for the
  ! entire matrix vector product. 

  use precision
  use hypData, only : ellipScatter, localSol, xx, npGlobal
  use hypInput
  use communication
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif
  ! Input Parameters
  Vec :: sigma

  ! Working Parameters
  integer(kind=intType) :: ierr,i

  call VecScatterBegin(ellipScatter, sigma, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd(ellipScatter, sigma, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Only need to write on root proc...it now has all the required information
  if (myid == 0) then
     ! Extract a pointer and set the strengths
     call VecGetArrayF90(localSol, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     open(unit=9,file=trim(sourceStrengthFile), status='replace')
     ! Format is pretty simple: A single interger with the number of
     ! panels, followed by that number of floats. 
5    format (I8)
6    format (g22.14)

     write(9, 5) nPGlobal
     do i=1,nPGlobal
        write(9, 6) xx(i)
     end do
     close(9)
     ! Always remember to restore the array
     call VecRestoreArrayF90(localSol, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if
end subroutine saveSolution

subroutine loadSolution(successful)

  ! loadSolution tries to load in an existing solution file. It then
  ! checks that the file *is* actually a solution. 
  !  
  ! Returns
  ! -------
  ! succesful : bool
  !    Flag signifiying if the load was successful or not. 

  use hypInput
  use hypData
  use panel
  use communication
  implicit none

  ! Ouput Parameters
  logical :: successful

  ! Working parameters
  integer(kind=intType) :: i, ierr
  logical :: failed, file_exists
  real(kind=realType) :: tempStrength(nPGlobal)
  real(kind=realType) :: V(3), phi, err, errLocal
  successful = .False. 

  ! Check 1: Does the file exists?
  if (myid == 0) then
     INQUIRE(FILE=trim(sourceStrengthFile), EXIST=file_exists) 
  end if

  call MPI_bcast(file_exists, 1, MPI_LOGICAL, 0, hyp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (.not. file_exists) then
     if (myid == 0) then
        print *, '-> sourceStrengthFile does not exist.'
     end if
     return
  end if

  ! Check 2: Does it contain the right number of panels?
  failed = .True.
  if (myid == 0) then
     open(unit=9, file=trim(sourceStrengthFile))
     read(9, *) i
     if (i == nPGlobal) then
        do i=1,nPGlobal
           read(9, *) tempStrength(i)
        end do
        failed = .False.
     end if
     close(9)
  end if

  call MPI_bcast(failed, 1, MPI_LOGICAL, 0, hyp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (failed) then
     if (myid == 0) then
        print *, '-> sourceStrengthFile has incorrect size.'
     end if
     return
  end if

  ! Step 3: Broadcast the strengths
  call MPI_bcast(tempStrength, nPGlobal, MPI_DOUBLE, 0, hyp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Stemp 4: Set all the strengths
  call setMGData(tempStrength)

  ! Now we evaluate the potential at the center of each panel (control
  ! point). This should be within kspTolerance of 1.0. (Full loop here
  ! since we don't have the PETSc vector partition yet)
  errLocal = zero
  do i=1, nPGlobal
     if (mod(i-1, nproc) == myid) then
        call evalAtPoint(MGP(1)%panels(i)%center, phi, V)
        errLocal = (phi - one)**2
     end if
  end do
  err = zero
  call mpi_allreduce(errLocal, err, 1, MPI_DOUBLE, MPI_SUM, hyp_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Complete the L2 norm
  err = sqrt(err)

  if (err < 10*kspRelTol) then
     successful = .True. 
  else
     if (myid == 0) then
        print *, '-> sourceStrengthFile has incorrect values.'
     end if
  end if

end subroutine loadSolution

subroutine marchPt(startPt, phiTarget, endPt, V, neval)

  ! marchPt evaluates the tracjectory of a point
  !  
  ! Parameters
  ! ----------
  ! pt : real, dimension(3)
  !   The point to use for evaluation
  !
  ! Returns
  ! -------
  use communication
  use hypData
  use hypInput
  use panel
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in), dimension(3) :: startPt
  real(kind=realType), intent(in) :: phiTarget
  
  ! Output Parameters
  real(kind=realType), intent(out), dimension(3) :: endPt, V
  integer(kind=intTYpe) :: neval

  ! Working Parameters
  real(kind=realType) :: y(3), t0, phi, k1(3), k2(3), V_tmp(3), phi_tmp
  real(kind=realType) :: yprime(3), tol, a(3), b(3), c(3), df(3), phi_prime, V_prime(3)
  real(kind=realType) :: yfinal(3), Vfinal(3), vtmp(3), dyfinal(3), dy(3), dt, fa, fb, ff
  real(kind=realType) :: dy_prime(3), dt0
  logical cont
  integer(kind=intTYpe):: i,j,k, kk 

  y = startPt
  tol = 1e-7
  t0 = zero
  ! Evaluate the starting point:
  call evalAtPoint(y, phi, V)
  neval = neval + 1

  ! Determine the timestep that would get us there:
  dt = (phi - phiTarget)/(sqrt(V(1)**2 + V(2)**2 + V(3)**2))
  dt0 = dt*2

  ! Integrate until we hit phi
  cont = .True.
  k = 0 
  do while(cont)
     k = k + 1

     if (k > 200) then
        print *,'pt march:',startPt
        print *, 'i:',i
        print *,'something screwed...', myid
        print *,'target:',phiTarget
        print *,'phi:',phi
        print *,'y:',y
        print *,'dt0:',dt0
        print *,'dt:',dt
        print *,'V:',V
        stop
     end if

     k1 = dt*V
     call evalAtPoint(y+half*k1, phi_tmp, V_tmp)
     neval = neval + 1

     k2 = dt*V_tmp
     yprime = y + k2
     call evalAtPoint(yprime, phi_prime, V_prime)
     neval = neval + 1

     if (abs(phi_prime - phiTarget) < tol) then
        ! We've just so happened to close enough to tol so we're
        ! done
        yfinal = yprime
        Vfinal = V_prime
        cont = .False.
        exit 
     else if (phi_prime < phiTarget) then
        !print *,'start bisection'
        ! Now just do a bisection search
        fa = phi - phiTarget
        fb = phi_prime - phiTarget
        a = y
        b = yprime

        do kk=1,25
           ! We want to solve f = phi - target = 0.0
           c = half*(a + b)
           call evalAtPoint(c, phi, df)
           neval = neval + 1

           ff = phi - phiTarget

           if (abs(ff) < tol) then
              cont = .False.
              yfinal = c
              Vfinal = df
              exit
           end if
           if (ff*fa > 0) then
              a = c
           else
              b = c
           end if
        end do

        if (cont) then
           if (abs(ff) > tol*100) then
              print *,'Bisection failed...'

           end if
           yfinal = c
           vFinal = df
           cont = .False.
        end if
     else

        ! Step is already ok
        yfinal = yprime
        y = yfinal
        V = V_prime
        phi = phi_prime

        if (mod(k, 5) == 0) then
           !  Adapt time step
           dt = dt * 5
        end if
     end if
  end do

  endPt = yfinal
  V = Vfinal
end subroutine marchPt

subroutine evalAtPoint(pt, phi, V)
  ! evalAtPoint is the main routine for determining the potential and
  ! velocity at a given pt. It uses one of three methods depending on
  ! the value in "evalMode" (in hypInput): 
  !
  ! 1. Exact Evaluation: The routine modifies the farfield tolerance
  ! to ensure that only the exact evaluations are used. (Slow N^2
  ! algorithm)
  !  
  ! 2. Slow Evaluation: Uses farfield approximations, but doesn't
  ! group panels (faster N^2 algorithm)
  !
  ! 3. Fast Evaluation: Uses farfield approximations and panel
  ! groupings (N*long(N) algorithm)
  !
  ! Parameters
  ! ----------
  ! pt : real, dimension(3)
  !   The point to use for evaluation
  !
  ! Return
  ! -------
  ! phi : real
  !   The computed potential
  ! V : real size(3)
  !   The computed velocity at point pt
  use hypInput
  use hypData
  use panel
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in), dimension(3) :: pt

  ! Output Parameters
  real(kind=realType), intent(out) :: phi
  real(kind=realType), dimension(3), intent(out) :: V

  ! Working parameters
  integer(kind=intType) :: i
  real(kind=realType) :: phiTmp, Vtmp(3), tmpSave

  V(:) = zero
  phi = zero

  if (evalMode == EVAL_FAST) then
     call evalAtLevel(levelMax, pt, coarseIndices, nCoarse, phi, V)     

  else if (evalMode == EVAL_SLOW) then
     do i=1, nPGlobal
        call panelInfluence2(i, pt, phiTmp, Vtmp)
        V = V + MGP(1)%panels(i)%strength*Vtmp
        phi = phi + MGP(1)%panels(i)%strength*phiTmp
     end do

  else if (evalMode == EVAL_EXACT) then
     tmpSave = farFieldTol
     farFieldTol = huge(phi)
     do i=1, nPGlobal
        call panelInfluence2(i, pt, phiTmp, Vtmp)
        V = V + MGP(1)%panels(i)%strength*Vtmp
        phi = phi + MGP(1)%panels(i)%strength*phiTmp
     end do
     farFieldTol = tmpSave
  end if

end subroutine evalAtPoint

subroutine setupPanels
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: setupPanels generates all the required data
  !     structures for doing the fast multipole method. This
  !     subroutine should only need to be called once from python.
  !  

  use communication
  use hypInput
  use hypData
  use panel
  implicit none
  
  ! Working 
  integer(kind=intType) :: iPatch, sizes(2), i, j, ii, jj, iPanel, nn
  integer(kind=intType) :: iLevel, ierr, iNode
  real(kind=realType), dimension(3) :: xCen, r , v1, v2, v3
  real(kind=realType) :: lenV3, dist

  ! Before we can do anything we must determine the number of
  ! multigrid levels we will be using. This is required to allocate
  ! the MGP with the correct number of levels.

  if (myid == 0) then
     levelMax = 1000 ! Arbitrary large number 
     do iPatch=1, nPatch
        sizes(1) = patches(iPatch)%il
        sizes(2) = patches(iPatch)%jl
        
        levelLoop: do j=1,levelMax
           if (.not. mod(sizes(1)-1, 2) == 0) then
              levelMax = j 
              exit levelLoop
           end if
           if (.not. mod(sizes(2)-1, 2) == 0) then
              levelMax = j
              exit levelLoop
           end if
           sizes(1) = (sizes(1)-1)/2 + 1
           sizes(2) = (sizes(2)-1)/2 + 1
        end do levelLoop
     end do
     
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") " Multigrid levels:"
     write(*,"(I2,1x)",advance="no") levelMax
     print "(1x)" 
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"   
     npGlobal = faceTotal
  end if

  call MPI_Bcast(npGlobal, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MPI_Bcast(levelMax, 1, MPI_INTEGER, 0, hyp_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterCreateToAll(X(1), allScatter, allGlobalNodes, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterBegin(allScatter, X(1), allGlobalNodes, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(allScatter, X(1), allGlobalNodes, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecScatterDestroy(allScatter, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(allGlobalNodes, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now we can allocate MGP
  allocate(MGP(levelMax))
  ! Now we proceed to fill up the "fine" level, ie. level 1
  allocate(MGP(1)%panels(nPGlobal))
  MGP(1)%N = nPGlobal
  ! Generate the panel mesh for the elliptic method:
  do i=1, nPGlobal
     pp => MGP(1)%panels(i)
     pp%N = 4
     pp%center = zero
     pp%nChildren = 0
     pp%strength = zero
     pp%length = zero
     ! Just take it from conn:
     do j=1, pp%N
        iNode = fullConn(j, i)
        pp%X(:, j) = xx(3*iNode-2:3*iNode)
     end do
     
     ! Compute panel center
     do j=1, pp%N
        pp%center = pp%center + pp%X(:, j) / pp%N
     end do

     ! Compute the "length" of the panel...that is the maximal linear
     ! dimension between nodes ---just do it n^2
     do j=1, pp%N
        do jj=1, pp%N
           pp%length = max(pp%length, dist(pp%X(:, j), pp%X(:, jj)))
        end do
     end do
 
     ! Now compute the normal and area
     v1 = pp%X(:, 3) - pp%X(:, 1)
     v2 = pp%X(:, 4) - pp%X(:, 2)
        
     ! v3 = v1 x v2
     v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
     v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
     v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
        
     ! Sum this much of the area
     lenv3 = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
     v3 = v3 / lenv3
     pp%area = half*lenV3 
     pp%normal = v3

     ! Compute the transformation matrics
     ! Following Hess, take the xaxis to be from p1 to p3 (for quads)
     v1 = pp%X(:, 3) - pp%X(:, 1)

     ! Normalize the first direction
     v1 = v1 / sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
 
     ! Panel normal
     v3 = pp%normal

     ! And we can get the (normalized) axis2 by doing one more cross
     ! product: axis2 = axis3 x axis 1
     v2(1) = v3(2) * v1(3) - v3(3) * v1(2)
     v2(2) = v3(3) * v1(1) - v3(1) * v1(3)
     v2(3) = v3(1) * v1(2) - v3(2) * v1(1)

     ! Finally finish our transformation matrix:
     pp%C(1, :) = v1
     pp%C(2, :) = v2
     pp%C(3, :) = v3
     pp%CT = transpose(pp%C)

     ! We don't actually need the physical coordinates for the panel,
     ! but the coordinates in the local frame, pts. So compute them
     ! here once. 
     pp%X(:, pp%N+1) = pp%X(:,pp%N)
     do j=1, pp%N
        pp%pts(:, j) = matmul(pp%C, pp%X(:, pp%N+1-j) - pp%center)
     end do
     pp%pts(:, pp%N+1) = pp%pts(:, 1)
  end do

  call VecRestoreArrayF90(allGlobalNodes, xx, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(allGlobalNodes, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now we can compute the index grouping information for the coarser
  ! grid levels: Because of the MG approach we can know exactly how
  ! many panels we be on each level:
  if (myid == 0) then
     write(*, "(a,I2,a,I6,a)"), 'Level ',1, ' has ',nPglobal, ' panels'
  end if

  MGLoop: do iLevel = 2, levelMax

     ! First determine the number of grouped panels....This comes
     ! purely from the MG-ptch topology
     ii = 0
     do iPatch=1, nPatch
        if (myid == 0) then
           sizes(1) = patches(iPatch)%il - 1
           sizes(2) = patches(iPatch)%jl - 1
        end if
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! And divide by the the level:
        sizes = sizes / 2**(iLevel - 1)
        
        ! And sum:
        ii = ii + sizes(1)*sizes(2)
     end do
     if (myid == 0) then
        write(*, "(a,I2,a,I6,a)"), 'Level ',iLevel, ' has ',ii, ' panels'
     end if
     ! Now allocate this level
     allocate(MGP(iLevel)%panels(ii))
     MGP(iLevel)%N = ii ! And the size for reference

     iPanel = 0
     ! Loop back over the patches and set algommeration
     do iPatch=1, nPatch
        ! Take CELL Sizes
        if (myid == 0) then
           sizes(1) = patches(iPatch)%il - 1
           sizes(2) = patches(iPatch)%jl - 1
        end if
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, hyp_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! And divide by the the level:
        sizes = sizes / 2**(iLevel - 1)

        do j=1,sizes(2)
           do i=1,sizes(1)
              
              iPanel = iPanel + 1
              pp => MGP(iLevel)%panels(iPanel)
              pp%nChildren = 0
              pp%Children = 0
              pp%nChildren = 4
              do ii=1,4
                 pp%children(ii) = 4*(iPanel-1) + ii
              end do
           end do
        end do
     end do
  end do MGLoop

  ! Last thing: We will allocate the 'coarseIndices' array which is
  ! just an array of 1,2,3...N where N is the number of coarsest grid
  ! panels
  nCoarse = MGP(levelMax)%N
  allocate(coarseIndices(nCoarse))
  do i=1,nCoarse
     coarseIndices(i) = i
  end do
  if (myid==0) then
     open(unit=9,file='test_fortran.dat', status='old')
5    format (a)
6    format (g20.15, g20.15, g20.15)
     write(9, 5) "variables= x y z"
     
     do i=1,nPglobal
        pp => MGP(1)%panels(i)
        if (pp%N == 4) then
           write(9, 5) "zone NODES=4, ELEMENTS=1, DATAPACKING=POINT"
           write(9, 5) "ZONETYPE=FEQUADRILATERAL"
           write(9, 6) pp%X(1, 1), pp%X(2, 1), pp%X(3, 1)
           write(9, 6) pp%X(1, 2), pp%X(2, 2), pp%X(3, 2)
           write(9, 6) pp%X(1, 3), pp%X(2, 3), pp%X(3, 3)
           write(9, 6) pp%X(1, 4), pp%X(2, 4), pp%X(3, 4)
           write(9, 5) "1 2 3 4"
        end if
     end do
     close(9) 
  end if
end subroutine setupPanels

subroutine setMGData(strengths)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: Given a (global) set of strengths, update the
  !     information in each of the MG levels
  !    
  !     Parameters
  !     ----------
  !     strengths : real array size nPGlobal

  use hypData
  use panel
  implicit none

  ! Input Params
  real(kind=realType), intent(in) :: strengths(nPGlobal)

  ! Working Param
  integer(kind=intType) :: iLevel, i, j
  real(kind=realType) :: childArea, childCenter(3), childStrength, lMax


  ! Just copy in level 1
  do i=1,nPGlobal
     MGP(1)%panels(i)%strength = strengths(i)
  end do

  ! Loop over the coarser levels
  do iLevel=2, levelMax
     do i=1,MGP(iLevel)%N
        pp => MGP(iLevel)%panels(i)
        pp%area = zero
        pp%center = zero
        pp%strength = zero
        lMax = zero
        do j=1,pp%nChildren
           ! Extract child data
           childArea = MGP(iLevel-1)%panels(pp%children(j))%area
           childStrength = MGP(iLevel-1)%panels(pp%children(j))%strength
           childCenter = MGP(iLevel-1)%panels(pp%children(j))%center
           lMax = max(lMax, MGP(iLevel-1)%panels(pp%children(j))%length)
           ! Sum contribution to the (coarser) panel
           pp%area = pp%area + childArea
           pp%center = pp%center + childArea*childCenter
           pp%strength = pp%strength + childArea*childStrength
        end do
        
        ! Set the "length" of the grouped panel to be twice the length
        ! of the child
        pp%length = lMax * two

        ! Finally divide by the summed area
        pp%center = pp%center / pp%area
        pp%strength = pp%strength / pp%area
     end do
  end do
end subroutine setMGData

recursive subroutine evalAtLevel(level, pt, indices, nIndices, phi, V)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: This is the core routine of the fast-multipole
  !     method. The purose is to evaluate the influence at all panels
  !     for a given point without actually looping over all
  !     N-panels. This routine will normally be called with the
  !     largest multigrid level. Then the routine will (recursively)
  !     call itself to evalaute any lower levels (all the way down to
  !     level 1) if necessary. 
  !    
  !     Parameters
  !     ----------
  !     level : integer  
  !         The multigrid level to evaluate on 
  !     pt : real size (3)
  !         The spatial coordinate of the point of iterest
  !     indices : integer array size (nIndices)
  !         The list of indices to evaluate on the given level
  !     nIndices : integer
  !         The size of the above array
  !     
  !     Returns
  !     -------
  !     phi : real
  !        The computed potential at point 'pt'
  !     V : real array (3)
  !        The compute velocity (gradient of the potential) at point pt

  use panel
  use hypInput
  use hypData
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, nIndices
  real(kind=realType), intent(in), dimension(3) :: pt
  integer(kind=intType), intent(in), dimension(nIndices) :: indices

  ! Ouptut Parameters
  real(kind=realType), intent(inout) :: phi
  real(kind=realType), dimension(3), intent(inout) :: V

  ! Working parameters
  integer(kind=intType) :: i, iPanel
  real(kind=realType), dimension(3) :: r, Vtmp, center
  real(kind=realType) :: d2, d, phiTmp, d32, area, strength, length
  
  do i=1,nIndices
     iPanel = indices(i)
     center = MGP(level)%panels(iPanel)%center
     area = MGP(level)%panels(iPanel)%area
     strength = MGP(level)%panels(iPanel)%strength
     length = MGP(level)%panels(iPanel)%length

     ! The far-field computation is duplicated here:
     r = pt - center
     d2 = r(1)**2 + r(2)**2 + r(3)**2
     d = sqrt(d2)

     if (d > farFieldTol * length) then 
        phi = phi + strength*area / d
        nFar = nFar + 1
        d32 = d2**(three/two)
        V(:) = V(:) + area*strength * r / d32
     else
        if (level > 1) then
           ! Need to call one level coarsen level
           call evalAtLevel(level-1, pt, MGP(level)%panels(iPanel)%children, &
                MGP(level)%panels(iPanel)%nChildren, phi, V)
        else
           !Finally we call the actual eval function --- only on level 1
           call panelInfluence2(iPanel, pt, phiTmp, VTmp)
           phi = phi + phiTmp*strength
           V = V + Vtmp*strength
        end if
     end if

  end do
end subroutine evalAtLevel

! subroutine panelInfluence(i, pt, phi, V)

!   ! panelInfluence generates the influence of panel 'i' on point
!   !  'pt'.
!   !  
!   ! Parameters
!   ! ----------
!   ! i : int
!   !   The index of the panel to use. 
!   ! pt : real size (3)
!   !   Three dimensional point to get the influence at
!   !
!   ! Returns
!   ! -------
!   ! phi : real
!   !   The computed potential
!   ! V : real size(3)
!   !   The computed velocity at point pt
!   !
!   ! Notes
!   ! -----

!   ! This routines uses the gloablly stored surface mesh and the
!   ! connectivity array to determine the nodes of panel 'i'.
!   use hypInput
!   use hypData
!   use panel
!   implicit none

!   ! Input Parameters
!   integer(kind=intType), intent(in) :: i
!   real(kind=realType), intent(in) :: pt(3)

!   ! Output Parameters
!   real(kind=realtype), intent(out) :: phi, V(3)

!   ! Working Parameters
!   real(kind=realType), dimension(3, NNodesMax+1) :: pts
!   real(kind=realType), dimension(nNodesMax+1) :: r
!   real(kind=realType), dimension(3) :: c, edge, rr
!   real(kind=realType) :: d, dTheta, id12, c12, s12, s12_1, d2, d32
!   real(kind=realType) :: r1, r2, s12_2, R12, Q12, J12
!   real(kind=realType) ::  tmp
!   integer(kind=intType) :: iNode, ii
  
!   V(:) = zero
!   phi = zero
  
!   ! Poiner to panel to make code easier to read
!   pp => MGP(1)%panels(i)

!   ! vector from the panel center to the point of interest:
!   rr = pt - pp%center
 
!   ! Cartiesian distance
!   d2 = rr(1)**2 + rr(2)**2 + rr(3)**2
!   d = sqrt(d2)

!   ! Check if we can do the farfield tolerance:
!   if (d > farFieldTol * pp%length) then
!      phi = pp%area / d
!      d32 = d2**(three/two)
!      V(:) = pp%area * rr /d32
!      nFar = nFar + 1
!      return
!   end if
!   nNear = nNear + 1
!   ! Transform the vector of point into the local frame as well as the
!   ! 4 corner points
!   rr = matmul(pp%C, rr)
 
!   do ii=1,pp%N
!      pts(:, ii) = matmul(pp%C, pp%x(:, pp%N+1-ii) - pp%center)
!      r(ii) = sqrt((rr(1) - pts(1, ii))**2 + (rr(2) - pts(2, ii))**2 + rr(3)**2)
!   end do
!   pts(:, pp%N+1) = pts(:, 1)
!   r(pp%N+1) = r(1)
!   dTheta = 2*pi
 
!   ! Loop over each of the 'N' edges
!   do iNode=1, pp%N
!      edge = pts(:, iNode + 1) - pts(:, iNode)
!      id12 = one/sqrt(edge(1)**2 + edge(2)**2)

!      ! Hess equation 4.4.10
!      C12 = (pts(1, iNode+1) - pts(1, iNode))*id12
!      S12 = (pts(2, iNode+1) - pts(2, iNode))*id12

!      ! Hess equation 4.4.12
!      s12_1 = (pts(1, iNode)   - rr(1))*C12 + (pts(2, iNode  ) - rr(2))*S12
!      s12_2 = (pts(1, iNode+1) - rr(1))*C12 + (pts(2, iNode+1) - rr(2))*S12

!      ! Hess equation 4.4.13
!      R12 = (rr(1) - pts(1, iNode))*S12 - (rr(2) - pts(2, iNode))*C12
  
!      if (R12 < zero) then
!         dTheta = zero
!      end if

!      ! Hess equation 4.4.14
!      r1 = r(iNode)
!      r2 = r(iNode+1)

!      ! Hess equation 4.4.15
!      Q12 = log((r2 + s12_2)/(r1 + s12_1))

!      ! Hess equation 4.4.16
!      J12 = atan2( (R12*abs(rr(3))*(r1*s12_2 - r2*s12_1)) , &
!           (r1*r2*R12**2 + rr(3)**2*s12_1*s12_2))
!      ! Hess equation 4.4.17 & 4.4.18
!      phi = phi + R12*Q12 + abs(rr(3))*J12

!      ! Hess equation 4.4.19
!      V(1) = V(1) - S12*Q12
!      V(2) = V(2) + C12*Q12
!      V(3) = V(3) - J12
!   end do

!   ! dTheta correction to Hess Equation 4.4.18 and 4.4.19
!   phi = phi - abs(rr(3))*dTheta

!   V(3) = dTheta + V(3)
!   if (rr(3) < zero) then
!      V(3) = -V(3)
!   end if

!   ! Transform the velocity back to the global frame
!   V = matmul(pp%CT, V)

! end subroutine panelInfluence

subroutine panelInfluence2(i, pt, phi, V)

  ! panelInfluence generates the influence of panel 'i' on point
  !  'pt'.
  !  
  ! Parameters
  ! ----------
  ! i : int
  !   The index of the panel to use. 
  ! pt : real size (3)
  !   Three dimensional point to get the influence at
  !
  ! Returns
  ! -------
  ! phi : real
  !   The computed potential
  ! V : real size(3)
  !   The computed velocity at point pt
  !
  ! Notes
  ! -----

  ! This routines uses the gloablly stored surface mesh and the
  ! connectivity array to determine the nodes of panel 'i'.
  use hypInput
  use hypData
  use panel
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: i
  real(kind=realType), intent(in) :: pt(3)

  ! Output Parameters
  real(kind=realtype), intent(out) :: phi, V(3)

  ! Working Parameters
  real(kind=realType), dimension(3, NNodesMax+1) :: pts
  real(kind=realType), dimension(nNodesMax+1) :: r
  real(kind=realType), dimension(3) :: edge, rr, rrg, vv
  real(kind=realType) :: d, dTheta, id12, c12, s12, s12_1, d2, d32
  real(kind=realType) :: r1, r2, s12_2, R12, Q12, J12
  real(kind=realType) ::  center(3), area, length, C(3,3)
  integer(kind=intType) :: iNode, ii, nn
  
  ! Extract the values we need (first only for the farfield calc)
  center = MGP(1)%panels(i)%center
  area = MGP(1)%panels(i)%area
  length = MGP(1)%panels(i)%length

  ! vector from the panel center to the point of interest:
  rrg = pt - center
 
  ! Cartiesian distance
  d2 = rrg(1)**2 + rrg(2)**2 + rrg(3)**2
  d = sqrt(d2)

  ! Check if we can do the farfield tolerance:
  if (d > farFieldTol * length) then
     phi = area / d
     d32 = d2**(three/two)
     V(:) = area * rrg /d32
     nFar = nFar + 1
     return
  end if

  ! For the local panel calc, we need to know the transformation
  ! matrix and zero out the potential and velocity in the local frame, vv
  C = MGP(1)%panels(i)%C
  nn = MGP(1)%panels(i)%N
  pts = MGP(1)%panels(i)%pts

  nNear = nNear + 1

  ! Transform the vector of point into the local frame as well as the
  ! 4 corner points
  rr(1) = C(1,1)*rrg(1) + C(1,2)*rrg(2) + C(1,3)*rrg(3)
  rr(2) = C(2,1)*rrg(1) + C(2,2)*rrg(2) + C(2,3)*rrg(3)
  rr(3) = C(3,1)*rrg(1) + C(3,2)*rrg(2) + C(3,3)*rrg(3)
  
  vv = zero
  phi = zero
  
  do ii=1, nn
     r(ii) = sqrt((rr(1) - pts(1, ii))**2 + (rr(2) - pts(2, ii))**2 + rr(3)**2)
  end do
  r(nn+1) = r(1)
  dTheta = 2*pi
 
  ! Loop over each of the 'N' edges
  do iNode=1, nn
     edge = pts(:, iNode + 1) - pts(:, iNode)
     id12 = one/sqrt(edge(1)**2 + edge(2)**2)

     ! Hess equation 4.4.10
     C12 = (pts(1, iNode+1) - pts(1, iNode))*id12
     S12 = (pts(2, iNode+1) - pts(2, iNode))*id12

     ! Hess equation 4.4.12
     s12_1 = (pts(1, iNode)   - rr(1))*C12 + (pts(2, iNode  ) - rr(2))*S12
     s12_2 = (pts(1, iNode+1) - rr(1))*C12 + (pts(2, iNode+1) - rr(2))*S12

     ! Hess equation 4.4.13
     R12 = (rr(1) - pts(1, iNode))*S12 - (rr(2) - pts(2, iNode))*C12
  
     if (R12 < zero) then
        dTheta = zero
     end if

     ! Hess equation 4.4.14
     r1 = r(iNode)
     r2 = r(iNode+1)

     ! Hess equation 4.4.15
     Q12 = log((r2 + s12_2)/(r1 + s12_1))

     ! Hess equation 4.4.16
     J12 = atan2( (R12*abs(rr(3))*(r1*s12_2 - r2*s12_1)) , &
          (r1*r2*R12**2 + rr(3)**2*s12_1*s12_2))
     ! Hess equation 4.4.17 & 4.4.18
     phi = phi + R12*Q12 + abs(rr(3))*J12

     ! Hess equation 4.4.19
     vv(1) = vv(1) - S12*Q12
     vv(2) = vv(2) + C12*Q12
     vv(3) = vv(3) - J12
  end do

  ! dTheta correction to Hess Equation 4.4.18 and 4.4.19
  phi = phi - abs(rr(3))*dTheta

  vv(3) = dTheta + vv(3)
  if (rr(3) < zero) then
     vv(3) = -vv(3)
  end if

  ! Transform the velocity back to the global frame
  V(1) = C(1,1)*vv(1) + C(2,1)*vv(2) + C(3,1)*vv(3)
  V(2) = C(1,2)*vv(1) + C(2,2)*vv(2) + C(3,2)*vv(3)
  V(3) = C(1,3)*vv(1) + C(2,3)*vv(2) + C(3,3)*vv(3)

end subroutine panelInfluence2

subroutine generateEllipLayer(X, layer)

  use precision
  use communication
  use panel
  use hypData, only : fullConn, rootScatter, xLocal, xx, npGlobal
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input/Output
  Vec X, layer

  ! Working
  integer(kind=intType) :: i, j, iNode, ierr, isize
  real(kind=realType), dimension(3) :: xCen
  
  ! Essentially we are just computing the average pt in each panel
  ! defined in fullConn and dumping them back into layer. For
  ! simplicity just do this on root proc for now:
  call VecScatterBegin(rootScatter, X, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(rootScatter, X, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  if (myid == 0) then
     call VecGetArrayF90(XLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     do i=1, npGlobal
        xCen = zero
        do j=1,4
           iNode = fullConn(j, i)
           xCen = xCen + xx(3*iNode-2:3*iNode)
        end do
        xCen = fourth*xCen
        
        call VecSetValuesBlocked(layer, 1, (/i-1/), xCen, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
     call VecRestoreArrayF90(XLocal, xx, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  call VecAssemblyBegin(layer, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(layer, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine generateEllipLayer

subroutine reconstruct(X, Xm1, layer, delta, norm)

  use precision
  use communication
  use panel
  use hypData, only : fullcPtr, rootScatter, xLocal, xxm1, npGlobal,nxGlobal
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input/Output
  Vec X, Xm1, layer, delta, norm

  ! Working
  integer(kind=intType) :: i, j, iCell, nn, ierr, rank
  real(kind=realType) :: ovrNeighbour
  real(kind=realType), dimension(3) :: dx, normal, S, x0, newFacePt
  real(kind=realType), dimension(:), pointer:: deltaPtr, normPtr, yy
  real(kind=realType), dimension(:, :), allocatable :: A
  real(kind=realType), dimension(:), allocatable :: B
  real(kind=realType) :: work(1000)
  integer(kind=intType) :: lwork, iwork(100), info
  ! Again, for simplicity we will just do this on the root proc for
  ! now. Notet that delta and norm are already full vectors. 

  call VecScatterBegin(rootScatter, Xm1, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(rootScatter, Xm1, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
     ! Get pointers to arrays
     call VecGetArrayF90(delta, deltaPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(norm, normPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(layer, yy, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(XLocal, xxm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     print *,'about to do reconstruct'

     do i=1,nXGlobal
        
        ! For now, just average the (cell) delta around each node
        nn = fullcPtr(1, i)
        ovrNeighbour = one/nn
        dx = zero
        do j=1, nn
           iCell = fullcPtr(j+1, i)
           dx = dx + deltaPtr(3*iCell-2:3*iCell)
        end do

        ! This is the average delta
        dx = dx * ovrNeighbour

        ! The (temporary) new node location:
        x0 = xxm1(3*i-2:3*i) + dx
        ! Now do it right
        allocate(A(nn+3, 3), B(nn+3))
        A(:, :) = zero
        B(:) = zero
        ! Add in the normals
        do j=1, nn
           iCell = fullcPtr(j+1, i)
           normal = normPtr(3*iCell-2:3*iCell)
           ! Need to normalize the normal
           normal = normal / sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
           A(j, :) = normal
           
           newFacePt = yy(3*iCell-2:3*iCell) + deltaPtr(3*iCell-2:3*iCell)
           B(j) = dot_product(normal, newFacePt - X0)
        end do

        A(nn+1, 1) = one
        A(nn+2, 2) = one
        A(nn+3, 3) = one
        B(nn+1:nn+3) = zero!xxm1(3*i-2:3*i)

        ! if (i==1449) then
        !    print *, 'shape of A:', i
        !    print *, 'orig x:', xxm1(3*i-2:3*i)
        !    do j=1,nn+3
        !       print *, A(j, :)
        !    end do
        !    do j=1,nn+3
        !       print *, B(j)
        !    END do
        ! end if
        
        ! Now solve the least squares problem
        lwork = 1000
        call dgelsd(nn+3, 3, 1, A, nn+3, B, nn+3, S, -one, rank, work, lwork, iwork, info)

        ! if (i==1449 .or. i == 1) then
        !    print *, ' -- i:', i
        !    print *, 'delta:', B(1:3)
        !    print *, 'orig x:', xxm1(3*i-2:3*i)
        !    print *, 'X0:', X0
        ! end if
        if (nn == 3) then
           B(:) = zero
        end if
        !B(:) = zero
        ! Finally set actual value
        call VecSetValuesBlocked(X, 1, (/i-1/), B(1:3)+X0, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        deallocate(A, B)
     end do

     ! Restore all pointers
     call VecRestoreArrayF90(delta, deltaPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(norm, normPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(layer, yy, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(XLocal, xxm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Finish assembly of our brand new X!
  call VecAssemblyBegin(X, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(X, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine reconstruct

subroutine reconstruct2(X, Xm1, layer, delta, norm)

  use precision
  use communication
  use panel
  use hypData, only : fullcPtr, rootScatter, xLocal, xxm1, npGlobal,nxGlobal
  implicit none
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input/Output
  Vec X, Xm1, layer, delta, norm

  ! Working
  integer(kind=intType) :: i, j, iCell, nn, ierr, rank
  real(kind=realType) :: ovrNeighbour
  real(kind=realType), dimension(3) :: dx, normal, S, x0, newFacePt
  real(kind=realType), dimension(:), pointer:: deltaPtr, normPtr, yy
  real(kind=realType), dimension(:, :), allocatable :: A
  real(kind=realType), dimension(:), allocatable :: B
  real(kind=realType) :: work(1000)
  integer(kind=intType) :: lwork, iwork(100), info
  ! Again, for simplicity we will just do this on the root proc for
  ! now. Notet that delta and norm are already full vectors. 

  call VecScatterBegin(rootScatter, Xm1, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecScatterEnd(rootScatter, Xm1, XLocal, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
     ! Get pointers to arrays
     call VecGetArrayF90(delta, deltaPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(norm, normPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(layer, yy, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecGetArrayF90(XLocal, xxm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     print *,'about to do reconstruct'

     do i=1,nXGlobal
        
        ! For now, just average the (cell) delta around each node
        nn = fullcPtr(1, i)
        ovrNeighbour = one/nn
        dx = zero
        do j=1, nn
           iCell = fullcPtr(j+1, i)
           dx = dx + deltaPtr(3*iCell-2:3*iCell)
        end do

        ! This is the average delta
        dx = dx * ovrNeighbour

        ! The (temporary) new node location:
        x0 = xxm1(3*i-2:3*i) + dx
      
        ! Finally set actual value
        call VecSetValuesBlocked(X, 1, (/i-1/), X0, INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do

     ! Restore all pointers
     call VecRestoreArrayF90(delta, deltaPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(norm, normPtr, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(layer, yy, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecRestoreArrayF90(XLocal, xxm1, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Finish assembly of our brand new X!
  call VecAssemblyBegin(X, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call VecAssemblyEnd(X, ierr)
  call EChk(ierr, __FILE__, __LINE__)

end subroutine reconstruct2

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

  ! Working parameters
  integer(kind=intType) :: i, j, k, info, ierr, iLow, iHigh, iPanel
  integer(kind=intType) :: npLocal
  real(kind=realType) :: phi, timeA, timeB, V(3)
  real(kind=realType) :: d, dist, xcen(3)
  real(kind=realType), dimension(:), pointer :: vtmp, phitmp
  integer(kind=intType), dimension(:), allocatable :: onProc, offProc
  real(kind=realType), dimension(100) :: targets
  real(kind=realType) :: meshLine(3, N)
  integer(kind=intType) :: nTarget

  call MatCreateDense(hyp_comm_world, PETSC_DETERMINE, PETSC_DETERMINE, &
       nPGlobal, nPGlobal, PETSC_NULL_SCALAR, ellipMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatGetVecs(ellipMat, ellipSol, ellipRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecGetOwnershipRange(ellipRHS, iLow, iHigh, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  nPLocal = iHigh-iLow
  if (myid == 0) then
     print *,'Determining PC size...'
  end if
  
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

  onProc = onProc + 10
  offProc = offProc + 10
  print *,'Avg nz:',(sum(onProc) + sum(offProc))/dble(npglobal)
  call MatCreateAIJ(hyp_comm_world, nPLocal, nPLocal, &
       PETSC_DETERMINE, PETSC_DETERMINE,&
       0, onProc, 0, offProc, ellipPCMat, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  deallocate(onProc, offProc)

  ! call MatSetOption(ellipPCMat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  call KSPCreate(petsc_comm_world, ellipKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPSetFromOptions(ellipKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)
     
  call KSPGMRESSetRestart(ellipKSP, kspSubspacesize, ierr)
  call EChk(ierr, __FILE__, __LINE__)
     
  call KSPSetOperators(ellipKSP, ellipMat, ellipPCMat, DIFFERENT_NONZERO_PATTERN, ierr)
  call EChk(ierr, __FILE__, __LINE__)
       
  call KSPSetTolerances(ellipKSP, kspRelTol, 1e-16, 1e5, kspMaxIts, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPMonitorSet(ellipKSP, kspmonitordefault, PETSC_NULL_OBJECT, &
       PETSC_NULL_FUNCTION, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (myid == 0) then
     print *,'Assembling Matrix ...'
  end if
  do j=1,npLocal
       do i=1, nPGlobal
          call panelInfluence(i, MGP(1)%panels(j+iLow)%center, phi, V)
          call MatSetValues(ellipMat, 1, j + iLow - 1, 1, i-1, phi, INSERT_VALUES, ierr)
          call EChk(ierr, __FILE__, __LINE__)
          d = dist(MGP(1)%panels(i)%center, MGP(1)%panels(j+iLow)%center)
          if (d < MGP(1)%panels(i)%length) then
             call MatSetValues(ellipPCMat, 1, j + iLow - 1, 1, i-1, phi, INSERT_VALUES, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if
     end do
  end do
  
  call MatAssemblyBegin(ellipMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyBegin(ellipPCMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call MatAssemblyEnd(ellipMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatAssemblyEnd(ellipPCMat, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSet(ellipRhs, one, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGetPC(ellipKSP, pc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPsetPCSide(ellipksp, PC_RIGHT, ierr)
  call EChk(ierr, __FILE__, __LINE__)
   
  if (myid == 0) then
     print *,'Solving Matrix ...'
  end if

  call KSPSolve(ellipKSP, ellipRHS, ellipSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  call VecScatterCreateToAll(ellipSol, ellipScatter, localSol, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterBegin(ellipScatter, ellipSol, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecScatterEnd(ellipScatter, ellipSol, localSol, &
       INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do this later....
  ! call VecScatterDestroy(ellipScatter, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! call VecDestroy(localSol, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  call VecGetArrayF90(localSol, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now set the strengths
  call setMGData(xx)

  call VecRestoreArrayF90(localSol, xx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Generate some dummy target data
  targets(2) = 0.99
  do i=3,N ! N is the number of levels
     targets(i) = targets(i-1) - 0.01
  end do
  !    targets(98)= 0.025
  !    targets(99)= 0.0225
  !    targets(100)= 0.020

  ! Loop over each panel
  do i=iLow+1, iHigh
     if (myid == 0) then
        print *,'Percent:',dble(i-iLow)/(iHigh-iLow)
     end if
     meshLine(:, :) = zero
     call marchPt(MGP(1)%panels(i)%center, targets, meshLine)
     
     ! Now copy mesh line into each of the 'X' petsc vectors. Since
     ! the X vectors are stored 'k-plane' wise, we have 
     do j=1, N
        call vecSetValuesBlocked(X(j), 1, (/i-1/), meshLine(:, j), INSERT_VALUES, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Now do all the assembly
  do j=1,N
     call VecAssemblyBegin(X(j), ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call VecAssemblyEnd(X(j), ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end do
  ! call MPI_barrier(hyp_comm_world, ierr)



  ! ! Lets fill in phi and Vx into some of the metric data for
  ! ! debugging
  ! timeA = mpi_wtime()
  ! do k=2,N ! This is the march layer
  !    call VecGetArrayF90(X(k), xx, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)
    
  !    ! call VecGetArrayF90(metrics(k, 1), vtmp, ierr)
  !    ! call EChk(ierr, __FILE__, __LINE__)

  !    ! call VecGetArrayF90(metrics(k, 2), phitmp, ierr)
  !    ! call EChk(ierr, __FILE__, __LINE__)
  !     print *,'evaluating on layer..', k

  !    ! do j=1,nx
  !    !    !call exactEval(xx(3*j-2:3*j), phi, V)
  !    !    call slowEval(xx(3*j-2:3*j), phi, V)
  !    !    !call fastEval(xx(3*j-2:3*j), phi, V)
  !    !    vtmp(3*j-2:3*j) = V
  !    !    phitmp(3*j-2) = phi
  !    ! end do


  !    ! call VecRestoreArrayF90(metrics(k, 1), vtmp, ierr)
  !    ! call EChk(ierr, __FILE__, __LINE__)

  !    ! call VecRestoreArrayF90(metrics(k, 2), phitmp, ierr)
  !    ! call EChk(ierr, __FILE__, __LINE__)
  !    call VecSet(metrics(k, 1), zero, ierr)
  !    call VecSet(metrics(k, 2), zero, ierr)
  !    do j=1,nPLocal
  !       !call slowEval(MGP(1)%panels(j+iLow)%center, phi, V)
  !       xcen = fourth*(xx(conn(1, j+ilow)*3-2:conn(1,j+ilow)*3)+ &
  !                      xx(conn(2, j+ilow)*3-2:conn(2,j+ilow)*3)+ &
  !                      xx(conn(3, j+ilow)*3-2:conn(3,j+ilow)*3)+ &
  !                      xx(conn(4, j+ilow)*3-2:conn(4,j+ilow)*3))
  !       call slowEval(Xcen, phi, v)
  !       V = V /four
  !       phi = phi /four
  !       call VecSetValuesBLocked(metrics(k, 1), 1, conn(1, j+iLow)-1, V, ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 1), 1, conn(2, j+iLow)-1, V, ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 1), 1, conn(3, j+iLow)-1, V, ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 1), 1, conn(4, j+iLow)-1, V, ADD_VALUES, ierr)
        
  !       call VecSetValuesBLocked(metrics(k, 2), 1, conn(1, j+iLow)-1, (/phi, 0, 0/), ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 2), 1, conn(2, j+iLow)-1, (/phi, 0, 0/), ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 2), 1, conn(3, j+iLow)-1, (/phi, 0, 0/), ADD_VALUES, ierr)
  !       call VecSetValuesBLocked(metrics(k, 2), 1, conn(4, j+iLow)-1, (/phi, 0, 0/), ADD_VALUES, ierr)
  !    end do
  !    call VecRestoreArrayF90(X(k), xx, ierr)
  !    call EChk(ierr, __FILE__, __LINE__)
     
  !    call VecAssemblyBegin(metrics(k, 1), ierr)
  !    call VecAssemblyEnd(metrics(k, 1), ierr)
  !    call VecAssemblyBegin(metrics(k, 2), ierr)
  !    call VecAssemblyEnd(metrics(k, 2), ierr)

  ! end do

  ! call VecCopy(metrics(2, 1), metrics(1, 1), ierr)
  ! call VecCopy(metrics(2, 2), metrics(1, 2), ierr)

  ! timeB = mpi_wtime()
  ! print *,'time on proc:',timeB-timeA
end subroutine runElliptic

subroutine marchPt(pt, targets, meshLine)

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
  real(kind=realType), intent(in), dimension(3) :: pt
  real(kind=realType), intent(in), dimension(N) :: targets

  ! Output Parameters
  real(kind=realType), intent(out), dimension(3, N) :: meshLine

  ! Working Parameters
  real(kind=realType) :: y0(3), y(3), t0, V(3), phi, k1(3), k2(3), V_tmp(3), phi_tmp
  real(kind=realType) :: yprime(3), tol, a(3), b(3), c(3), df(3), phi_prime, V_prime(3)
  real(kind=realType) :: yfinal(3), Vfinal(3), vtmp(3), dyfinal(3), dy(3), dt, fa, fb, ff
  real(kind=realType) :: dy_prime(3)
  logical cont
  integer(kind=intTYpe):: i,j,k, kk 
  
  y0 = pt
  y =pt
  tol = 1e-6

  meshLine(:, 1) = pt
  do i=2,N
    t0 = zero

    ! Evaluate the starting point:
    call slowEval(y0, phi, V)

    ! Determine the timestep that would get us there:
    dt = (phi - targets(i))/(sqrt(V(1)**2 + V(2)**2 + V(3)**2))
    y = y0

    ! Integrate until we hit phi
    cont = .True.
    k = 0 
    do while(cont)
       k = k + 1

       if (k > 200) then
          print *,'something screwed...'
          print *,y
          stop
       end if

       k1 = dt*V
       call slowEval(y+half*k1, phi_tmp, V_tmp)
       k2 = dt*V_tmp
       yprime = y + k2
       call slowEval(yprime, phi_prime, V_prime)

       if (abs(phi_prime - targets(i)) < tol) then
          yfinal = yprime
          Vfinal = V_prime
          cont = .False.
          exit 
       else if (phi_prime < targets(i)) then

          ! Now just do a bisection search
          fa = phi - targets(i)
          fb = phi_prime - targets(i)
          a = y
          b = yprime

          do kk=1,25
             ! We want to solve f = phi - target = 0.0
             c = half*(a + b)
             call slowEval(c, phi, df)

             ff = phi - targets(i)

             if (abs(ff) < tol) then
                cont = .False.
                yfinal = c
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
                cont = .False.
                yfinal = c
             end if
          end if
       else
          ! Step is already ok
          yfinal = yprime
          phi = phi_prime
          if (mod(k, 5) == 0) then
             !  Adapt time step
             dt = dt * 4
          end if
       end if
    end do
    y = yfinal
    y0 = yfinal
    meshLine(:, i) = yfinal
 end do

end subroutine marchPt

subroutine exactEval(pt, phi, V)
  ! exactEval uses the N^2 algorithm to determine the potential and
  ! velocity at a given pt. It modified the farfield tolerance to
  ! ensure that only the exact evaluations are used. 
  !  
  ! Parameters
  ! ----------
  ! pt : real, dimension(3)
  !   The point to use for evaluation
  !
  ! Returns
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
  tmpSave = farFieldTol
  farFieldTol = huge(phi)
  V(:) = zero
  phi = zero
  do i=1, nPGlobal
     call panelInfluence(i, pt, phiTmp, Vtmp)
     V = V + MGP(1)%panels(i)%strength*Vtmp
     phi = phi + MGP(1)%panels(i)%strength*phiTmp
  end do
  farFieldTol = tmpSave
end subroutine exactEval

subroutine slowEval(pt, phi, V)
  
  ! slowEval uses the N^2 algorithm to determine the potential and
  ! velocity at a given pt. 
  !  
  ! Parameters
  ! ----------
  ! pt : real, dimension(3)
  !   The point to use for evaluation
  !
  ! Returns
  ! -------
  ! phi : real
  !   The computed potential
  ! V : real size(3)
  !   The computed velocity at point pt
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
  real(kind=realType) :: phiTmp, Vtmp(3)

  V(:) = zero
  phi = zero
  do i=1, nPGlobal
     call panelInfluence(i, pt, phiTmp, Vtmp)
     V = V + MGP(1)%panels(i)%strength*Vtmp
     phi = phi + MGP(1)%panels(i)%strength*phiTmp
  end do
end subroutine slowEval

subroutine fastEval(pt, phi, V)
  
  ! slowEval uses the Nlog(N) algorithm to determine the potential and
  ! velocity at a given pt. 
  !  
  ! Parameters
  ! ----------
  ! pt : real, dimension(3)
  !   The point to use for evaluation
  !
  ! Returns
  ! -------
  ! phi : real
  !   The computed potential
  ! V : real size(3)
  !   The computed velocity at point pt
  use hypData
  use panel
  implicit none
  
  ! Input Parameters
  real(kind=realType), intent(in), dimension(3) :: pt
  
  ! Output Parameters
  real(kind=realType), intent(out) :: phi
  real(kind=realType), dimension(3), intent(out) :: V

  ! Working parameters
  integer(kind=intType) :: nIndices

  V(:) = zero
  phi = zero
  call evalAtLevel(levelMax, pt, coarseIndices, nCoarse, phi, V)

end subroutine fastEval

subroutine setupPanels
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: setupPanels generates all the required data
  !     structures for doing the fast multipole method. This
  !     subroutine should only need to be called once from python.
  !  
  !  Parameters
  !  ----------

  !  None. This routine uses the global surface in XSurfGlobal as well
  !  as the patch information and connectivity information, all of
  !  which is set from python.
  use communication
  use hypInput
  use hypData
  use panel
  implicit none
  
  ! Working 
  integer(kind=intType) :: iPatch, sizes(2), i, j, ii, jj, iPanel, nn
  integer(kind=intType) :: iLevel, iInd, jInd, gIndex, iNode
  real(kind=realType), dimension(3) :: xCen, r , v1, v2, v3
  real(kind=realType) :: lenV3, dist
  integer(kind=intType), dimension(:), allocatable :: nodeUsed

  ! Before we can do anything we must determine the number of
  ! multigrid levels we will be using. This is required to allocate
  ! the MGP with the correct number of levels.

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

  if (myid == 0) then
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"  
     write(*,"(a)", advance="no") "Multigrid levels:"
     write(*,"(I2,1x)",advance="no") levelMax
     print "(1x)" 
     write(*,"(a)", advance="no") '#--------------------#'
     print "(1x)"   
   end if

  ! Now we can allocate MGP
  allocate(MGP(levelMax))

  ! Now we proceed to fill up the "fine" level, ie. level 1
  allocate(MGP(1)%panels(nPGlobal))
  MGP(1)%N = nPGlobal

  ! Generate the dual panel mesh for the ellipic method:
  do i=1,nPGlobal
     pp => MGP(1)%panels(i)
     ! "The number of nodes on a panel mesh (dual mesh) is equal to
     ! the number CELLS surrounding a node on the primal mesh. 

     ! Set default values
     pp%N = cPtr(1, i)
     NN = pp%N
     pp%area = zero
     pp%normal = zero
     pp%length = zero
     pp%center = zero
     pp%nChildren = 0
     pp%strength = zero
     pp%X = zero

     ! Now loop over the cells surronding the node, compute their
     ! centers. This forms the nodes for the dual mesh (panel mesh)
     do j=1, NN
        iPanel = cPtr(j+1, i)
        ! Now get the center of this panel:
        xCen = zero
        do ii=1,4
           xCen = xCen + xSurfGlobal(:, conn(ii, iPanel))
        end do
        
        pp%X(:, j) = fourth*xCen
     end do

     ! ! Just take it from conn:
     ! do j=1, NN
     !    iNode = conn(j, i)
     !    pp%X(:, j) = xSurfGlobal(:, iNode)
     ! end do
     
     ! Compute panel center --- this will be close to but not the same
     ! as the node on the primal mesh.
     do j=1, NN
        pp%center = pp%center + pp%X(:, j) / NN
     end do

     ! Compute the "length" of the panel...that is the maximal linear
     ! dimension between nodes ---just do it n^2
     do j=1,NN
        do jj=1,NN
           pp%length = max(pp%length, dist(pp%X(:, j), pp%X(:, jj)))
        end do
     end do
 
     ! Now compute the (average) normal as well as the area:
     do j=1, NN
        ii = j
        jj = j + 1
        if (jj > NN) then
           jj = jj - NN
        end if
      
        v1 = pp%X(:, ii) - pp%center
        v2 = pp%X(:, jj) - pp%center
        
        ! v3 = v1 x v2
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
        
        ! Sum this much of the area
        lenv3 = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
        v3 = v3 / lenv3
        pp%area = pp%area + half*lenV3 
        pp%normal = pp%normal + v3 / NN
     end do

     ! Now "Correct" the dual panels such that they lie
     ! "eps" BELOW the actual nodes, along the normal
     ! r = XSurfGlobal(:, i) - pp%center
     ! do j=1,NN
     !    pp%X(:, j) = pp%X(:, j) + r - eps*pp%normal
     ! end do

     ! Compute the transformation matrics
     if (pp%N == 4) then
        ! Following Hess, take the xaxis to be from p1 to p3 (for quads)
        v1 = pp%X(:, 3) - pp%X(:, 1)
     else
        v1 = pp%X(:, 2) - pp%X(:, 1)
     end if

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

     ! And the tranposer
     pp%CT = transpose(pp%C)
  end do

  ! Now move the "real" node in XSurfGlobal to be eps ABOVE the
  ! actual panel center which will be control point for the solution
  do i=1, nPGlobal 
     pp => MGP(1)%panels(i)
     XSurfGlobal(:, i) = pp%center + eps*pp%normal
  end do

  ! Now we can compute the index grouping information for the coarser
  ! grid levels: Because of the MG approach we can know exactly how
  ! many panels we be on each level:
  if (myid == 0) then
     write(*, "(a,I2,a,I6,a)"), 'Level ',1, ' has ',nPglobal, ' panels'
  end if
  allocate(nodeUsed(nPGlobal))
  MGLoop: do iLevel = 2, levelMax

     nodeUsed(:) = 0

     ! First determine the number of grouped panels....This comes
     ! purely from the MG-ptch topology
     ii = 0
     do iPatch=1, nPatch
        ! Take CELL Sizes
        sizes(1) = patches(iPatch)%il - 1 
        sizes(2) = patches(iPatch)%jl - 1

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
        sizes(1) = patches(iPatch)%il - 1 
        sizes(2) = patches(iPatch)%jl - 1

        ! And divide by the the level:
        sizes = sizes / 2**(iLevel - 1)

        do j=1,sizes(2)
           do i=1,sizes(1)
              
              iPanel = iPanel + 1
              pp => MGP(iLevel)%panels(iPanel)
              pp%nChildren = 0
              pp%Children = 0

              if (iLevel == 2) then
                 ! The sub-loop for level 2 is ALWAYS 3x3. 
                 do jj=-1,1
                    do ii=-1,1
                       iInd = 2*i + ii
                       jInd = 2*j + jj
                       gIndex = patches(iPatch)%l_index(iInd, jInd)
                       if (nodeUsed(gIndex) == 0) then
                          ! Not used:
                          nodeUsed(gIndex) = 1 ! Set as used
                          pp%nChildren = pp%nChildren + 1
                          pp%children(pp%nChildren) = gIndex
                       end if
                    end do
                 end do
              else ! Level 3 or higher. This is a little easier since
                   ! we know exactly the number of groupings on the
                   ! previous level. From here on out, there be
                   ! *precisely* 4 children per panel grouping,
                   ! furthermore the groupings will be
                   ! sequential. Ie. Panel 1 on level three will
                   ! contain panels 1,2,3,4 on level 2. Panel 2 on
                   ! level three will contain 5,6,7,8 etc. 
                 pp%nChildren = 4
                 do ii=1,4
                    pp%children(ii) = 4*(iPanel-1) + ii
                 end do
              end if
           end do
        end do
     end do
  end do MGLoop
  deallocate(nodeUsed)

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
6    format (g20.12, g20.12, g20.12)
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
  do iLevel = 2, levelMax
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
  real(kind=realType), dimension(3) :: r, Vtmp
  real(kind=realType) :: d2, d, phiTmp, d32
  
  do i=1,nIndices
     iPanel = indices(i)
     pp => MGP(level)%panels(iPanel)
     ! The far-field computation is duplicated here:

     r = pt - pp%center
     d2 = r(1)**2 + r(2)**2 + r(3)**2
     d = sqrt(d2)

     if (d > farFieldTol * pp%length) then
        phi = phi + pp%strength*pp%area / d
        nFar = nFar + 1
        d32 = d2**(three/two)
        V(:) = V(:) + pp%area*pp%strength * r / d32
     else
        if (level > 1) then
           ! Need to call one level coarsen level
           call evalAtLevel(level-1, pt, pp%children, pp%nChildren, phi, V)
        else
           !Finally we call the actual eval function --- only on level 1
           call panelInfluence(iPanel, pt, phiTmp, VTmp)
           phi = phi + phiTmp*pp%strength
           V = V + Vtmp*pp%strength
        end if
     end if
  end do
end subroutine evalAtLevel

subroutine panelInfluence(i, pt, phi, V)

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
  real(kind=realType), dimension(3) :: c, edge, rr
  real(kind=realType) :: d, dTheta, d12, c12, s12, s12_1, d2, d32
  real(kind=realType) :: r1, r2, s12_2, R12, Q12, J12
  real(kind=realType) ::  tmp
  integer(kind=intType) :: iNode, ii
  
  V(:) = zero
  phi = zero
  
  ! Poiner to panel to make code easier to read
  pp => MGP(1)%panels(i)

  ! vector from the panel center to the point of interest:
  rr = pt - pp%center
 
  ! Cartiesian distance
  d2 = rr(1)**2 + rr(2)**2 + rr(3)**2
  d = sqrt(d2)

  ! Check if we can do the farfield tolerance:
  if (d > farFieldTol * pp%length) then
     phi = pp%area / d
     d32 = d2**(three/two)
     V(:) = pp%area * rr /d32
     nFar = nFar + 1
     return
  end if
  nNear = nNear + 1
  ! Transform the vector of point into the local frame as well as the
  ! 4 corner points
  rr = matmul(pp%C, rr)
 
  do ii=1,pp%N
     pts(:, ii) = matmul(pp%C, pp%x(:, pp%N+1-ii) - pp%center)
     r(ii) = sqrt((rr(1) - pts(1, ii))**2 + (rr(2) - pts(2, ii))**2 + rr(3)**2)
  end do
  pts(:, pp%N+1) = pts(:, 1)
  r(pp%N+1) = r(1)
  dTheta = 2*pi
 
  ! Loop over each of the 'N' edges
  do iNode=1, pp%N
     edge = pts(:, iNode + 1) - pts(:, iNode)
     d12 = sqrt(edge(1)**2 + edge(2)**2)

     ! Hess equation 4.4.10
     C12 = (pts(1, iNode+1) - pts(1, iNode))/d12
     S12 = (pts(2, iNode+1) - pts(2, iNode))/d12

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
     V(1) = V(1) - S12*Q12
     V(2) = V(2) + C12*Q12
     V(3) = V(3) - J12
  end do

  ! dTheta correction to Hess Equation 4.4.18 and 4.4.19
  phi = phi - abs(rr(3))*dTheta

  V(3) = dTheta + V(3)
  if (rr(3) < zero) then
     V(3) = -V(3)
  end if

  ! Transform the velocity back to the global frame
  V = matmul(pp%CT, V)

end subroutine panelInfluence


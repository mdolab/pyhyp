subroutine surfaceSmooth(layer, nSteps, stepSize)
  use hypData
  implicit none

  ! Input Parameters
  integer (kind=intType), intent(in) :: nSteps, layer
  real(kind=realType), intent(in) :: stepSize

  ! Working Parameters
  integer(kind=intType) :: i, ii, iCell, iter, ierr
  real(kind=realType), dimension(3, nx) :: massCen
  real(kind=realType), dimension(3, nx+nghost) :: normals
  real(kind=realType), dimension(3) :: v1, v2, s, xdot, dx
  real(kind=realType), dimension(3, nLocalFace) :: cellCen
  real(kind=realType), dimension(nLocalFace) :: areas
  real(kind=realType) :: aTot, sNorm, fact

  masterIteration: do iter=1,nSteps

     ! Update the ghost entries before we start
     call VecGhostUpdateBegin(X(layer), INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecGhostUpdateEnd(X(layer), INSERT_VALUES,SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get local form for doing stuff
     call VecGhostGetLocalForm(X(layer), XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecGetArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     cellCen(:, :) = zero
     massCen(:, :) = zero
     normals(:, :) = zero

     do i=1, nLocalFace

        ! Compute face center:
        do ii=1,4
           cellCen(:, i) = cellCen(:, i) + fourth*xx(3*conn(ii, i)-2:3*conn(ii, i))
        end do

        ! Compute the normal of this face
        v1 = xx(3*conn(3, i)-2:3*conn(3, i)) - xx(3*conn(1, i)-2:3*conn(1, i))
        v2 = xx(3*conn(4, i)-2:3*conn(4, i)) - xx(3*conn(2, i)-2:3*conn(2, i))
        
        ! Cross product 
        s(1) = (v1(2)*v2(3) - v1(3)*v2(2))
        s(2) = (v1(3)*v2(1) - v1(1)*v2(3))
        s(3) = (v1(1)*v2(2) - v1(2)*v2(1))

        ! Compute the area
        snorm = sqrt(s(1)**2 + s(2)**2 + s(3)**2)
        
        ! Area:
        areas(i) = half * snorm
        
        ! Now normalize:
        s = s / snorm
     
        ! Scatter to each node
        do ii=1,4
           normals(:, conn(ii, i)) = normals(:, conn(ii, i)) + s
        end do
     end do

     ! Now make the second pass over nodes and re-normalize. Also compute
     ! the mass center for each node
     do i=1, nx
        s = normals(:, i)
        normals(:, i) = s / sqrt(s(1)**2 + s(2)**2 + s(3)**2)
     
        aTot = zero
        ! Loop over each cell around the node
        do ii=1, cPtr(1, i)
           iCell = cPtr(ii+1, i) 
           massCen(:, i) = massCen(:, i) + cellCen(:, iCell)*areas(iCell)
           aTot = aTot + areas(iCell)
        end do
        massCen(:, i) = massCen(:, i) / aTot
     end do
     ! Now perform the explict euler udpate
     do i=1,nx
        !fact = min((dist(i)/ dStar)**1.5, 1.0_8)
        fact = one
        dx = massCen(:, i) - xx(3*i-2:3*i)
        xDot = dx - dot_product(normals(:, i), dx)*normals(:, i)
        xx(3*i-2:3*i) = xx(3*i-2:3*i) + xDot*stepSize*fact
     end do

     ! Restore everything to prepare for the next iteration
     call VecRestoreArrayF90(XL_local, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecGhostRestoreLocalForm(X(layer), XL_local, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
  end do masterIteration
  
end subroutine surfaceSmooth



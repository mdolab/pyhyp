! CGNS reader
! Author: Ney Secco
! Date: 11-2015

module readCGNS

contains

  subroutine readGrid(cgnsFile,patches)

    ! This subroutine opens a cgnsFile and 

    use hypData, only: patchType, intType, realType, maxBCneighbors, SplayBCindex, &
         SymmetryBCindex, ConstXBCindex, ConstYBCindex, ConstZBCindex
    implicit none
    include 'cgnslib_f.h'

    ! Input Arguments
    character*(*),intent(in) :: cgnsFile

    ! Output Arguments
    type(patchType), allocatable, dimension(:), intent(out) :: patches

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

    ! Do the I/O that is common to both types of grids

    print *, ' -> Reading CGNS File: ', cgnsFile

    ! Open and get the number of zones:
    call cg_open_f(trim(cgnsFile), CG_MODE_READ, cg, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    call cg_nbases_f(cg, nbases, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    if (nbases .gt. 1) then
       print *, ' ** Warning: pyHyp only reads the first base in a cgns file'
    end if
  
    base = 1_intType
    
    call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    if (cellDim .ne. 2 .or. PhysDim .ne. 3) then
       print *, 'The Cells must be 2 dimensional'
       stop
    end if
    
    call cg_nzones_f(cg, base, nZones, ierr);
    if (ierr .eq. CG_ERROR) call cg_error_exit_f
    
    print *, '   -> Number of Zones:', nzones
     
    ! Determine if we have all structured zones
    do i=1, nZones
       call cg_zone_type_f(cg, base, i, zoneType, ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
       
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
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
       
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
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
       
       do j=1, dims(2)
          do i=1, dims(1)
             patches(iZone)%X(1,i,j) = coor(i, j)
          end do
       end do
       
       ! Load Y coordinates into the auxilary variable (coor)
       call cg_coord_read_f(cg, base, iZone, "CoordinateY", RealDouble, &
            (/1,1/), dims(1:2),  coor, ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
       
       do j=1, dims(2)
          do i=1, dims(1)
             patches(iZone)%X(2,i,j) = coor(i, j)
          end do
       end do

       ! Load Z coordinates into the auxilary variable (coor)
       call cg_coord_read_f(cg, base, iZone, "CoordinateZ", RealDouble, &
            (/1,1/), dims(1:2),  coor, ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f
       
       do j=1, dims(2)
          do i=1, dims(1)
             patches(iZone)%X(3,i,j) = coor(i, j)
          end do
       end do
       
       !Deallocate coor for next iteration
       deallocate(coor)
       
       ! Reading boundary conditions
       call cg_nbocos_f(cg, base, iZone, nbocos, ierr)
       if (ierr .eq. CG_ERROR) call cg_error_exit_f

       ! Allocate and initialize BC arrays
       allocate(patches(iZone)%BCnodes(dims(1),dims(2),1+2*maxBCneighbors))
       ! dims(1) and dims(2) -> because we need BC info for all i's and j's of the patch (all nodes)
       ! 2*maxBCneighbors+1 -> because we need to store bcType and i and j coordinates of each neighbor
       patches(iZone)%BCnodes = 0

       do iBC = 1,nbocos

          ! Get BC information
          call cg_boco_info_f(cg, base, iZone, iBC, bocoName, bocoType, &
               ptset_type, npnts, NormalIndex, NormalListSize, &
               NormalDataType, ndataset, ierr)
          if (ierr .eq. CG_ERROR) call cg_error_exit_f

          ! Get point range
          call cg_boco_read_f(cg, base, iZone, iBC, pnts, NormalList, ierr)
          if (ierr .eq. CG_ERROR) call cg_error_exit_f

          ! Now we set BC index depending on BC type
          ! The indices are defined in hypData.F90
          ! Corners may have 2 distinct BCs defined in the CGNS file, but we
          ! need to apply only one of them. The criteria is that the higher
          ! indices have higher priority.

          if (bocoType .eq. BCExtrapolate) then ! This is the splay BC
             bocoIndex = SplayBCindex

          else if (bocoType .eq. BCSymmetryPlane) then
             bocoIndex = SymmetryBCindex

          else if (bocoType .eq. BCWall) then
             ! We need to find if the user wants an X, Y or Z Wall
             ! by looking at the Velocity components of the Wall BC

             ! Loop over the datasets to find the wanted names.
             ! The last dataset will replace the other ones.
             do idataset = 1,ndataset

                ! Check the dataset names
                call cg_dataset_read_f(cg, base, iZone, iBC, idataset, &
                     DatasetName, dummybocoType, DirichletFlag, NeumannFlag, ierr)
                if (ierr .eq. CG_ERROR) call cg_error_exit_f

                ! Check if it has Dirichlet data
                if (DirichletFlag .eq. 1) then
                   ! Try to access VelocityX
                   call cg_goto_f(cg, base, ierr, 'Zone_t',iZone, 'ZoneBC',0, &
                        'BC_t',iBC, 'BCDataSet_t',idataset, 'DirichletData',0, &
                        'VelocityX',0, 'end')
                   if (ierr .eq. 0) then ! We found VelocityX so we have WallX BC
                      bocoIndex = ConstXBCindex
                   else
                      ! Try to access VelocityY
                      call cg_goto_f(cg, base, ierr, 'Zone_t',iZone, 'ZoneBC',0, &
                           'BC_t',iBC, 'BCDataSet_t',idataset, 'DirichletData',0, &
                           'VelocityY',0, 'end')
                      if (ierr .eq. 0) then ! We found VelocityY so we have WallY BC
                         bocoIndex = ConstYBCindex
                      else
                         ! Try to access VelocityZ
                         call cg_goto_f(cg, base, ierr, 'Zone_t',iZone, 'ZoneBC',0, &
                              'BC_t',iBC, 'BCDataSet_t',idataset, 'DirichletData',0, &
                              'VelocityZ',0, 'end')
                         if (ierr .eq. 0) then ! We found VelocityZ so we have WallZ BC
                            bocoIndex = ConstZBCindex
                         else
                            print *,'Did not recognize Wall BC'
                            print *,'Check if you defined VelocityX, VelocityY, or VelocityZ'
                         end if
                      end if
                   end if
                end if

             end do

          else
             print *,'BC type not supported'
             stop

          end if

          ! Determine BC interval and orientation
          iMinBC = minval(pnts(1,:))
          iMaxBC = maxval(pnts(1,:))
          jMinBC = minval(pnts(2,:))
          jMaxBC = maxval(pnts(2,:))

          ! Add neighbor nodes to the BCnodes array. The neighbor nodes will depend
          ! on the BC type

          neighbor_case: if (bocoType .eq. BCExtrapolate) then ! We have a splay BC

             ! The neighbors are taken in a counter-clock-wise order starting from
             ! the marching direction shown below
             !
             !         <--
             !     ___________
             ! |  |           |
             ! |  |           |
             ! V  |           | ^
             !    |           | |
             ! j  |___________| |
             ! |
             ! +--i    -->
             !
             ! For instance, for a node in the jlow edge (lower dege), the ordering of
             ! the neighbors is as follows:
             !
             !   --X----2---X--
             !     |    |   |
             !   --3--node--1--
             !
             ! For corner nodes, we store the diagonals. For instance, for the upper-left corner we have
             ! the following ordering:
             !
             ! node--3
             !  |    |
             !  |    |
             !  1----2
             

             if (jMinBC .eq. jMaxBC) then ! Only i changes (j-edge BC)
                j = jMinBC
                if (j .eq. 1) then ! We are at lower j bound
                   
                   ! FIRST NODE (Lower-left corner)
                   i = iMinBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i + 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i + 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j + 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j + 1
                   end if
                   
                   ! MIDDLE NODES
                   do i = iMinBC+1,iMaxBC-1
                      ! Update BC type (we don't need to check because edge nodes
                      ! will only have one BC)
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i + 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j + 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i - 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end do
                   
                   ! LAST NODE (Lower-right corner)
                   i = iMaxBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j + 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i - 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j + 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i - 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end if
                   
                else ! We are at upper j bound
                   
                   ! FIRST NODE  (Upper-left corner)
                   i = iMinBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j - 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i + 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j - 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i + 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end if
                   
                   ! MIDDLE NODES
                   do i = iMinBC+1,iMaxBC-1
                      ! Update BC type (we don't need to check because edge nodes
                      ! will only have one BC)
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i - 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j - 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i + 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end do
                   
                   ! LAST NODE (Upper-right corner)
                   i = iMaxBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i - 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i - 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j - 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j - 1
                   end if
                   
                end if
                
             else if (iMinBC .eq. iMaxBC) then ! Only j changes (i-edge BC)
                i = iMinBC
                if (i .eq. 1) then ! We are at lower i bound
                   
                   ! FIRST NODE (lower-left corner)
                   j = jMinBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i + 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i + 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j + 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j + 1
                   end if
                   
                   ! MIDDLE NODES
                   do j = jMinBC+1,jMaxBC-1
                      ! Update BC type (we don't need to check because edge nodes
                      ! will only have one BC)
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j - 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i + 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j + 1
                   end do
                   
                   ! LAST NODE (Upper-left corner)
                   j = jMaxBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j - 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i + 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j - 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i + 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end if
                   
                else ! We are at upper i bound
                   
                   ! FIRST NODE (lower-right corner)
                   j = jMinBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j + 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i - 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j + 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i - 1
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j
                   end if
                   
                   ! MIDDLE NODES
                   do j = jMinBC+1,jMaxBC-1
                      ! Update BC type (we don't need to check because edge nodes
                      ! will only have one BC)
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j + 1
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i - 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j - 1
                   end do
                   
                   ! LAST NODE (Upper-right corner)
                   j = jMaxBC
                   ! Check the BC type and change it if the new BC has higher
                   ! priority (higher index)
                   if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                      ! Update BC type
                      patches(iZone)%BCnodes(i,j,1) = bocoIndex
                      ! Update neighbors
                      ! i-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,2) = i - 1
                      ! j-index of 1st neighbor
                      patches(iZone)%BCnodes(i,j,3) = j
                      ! i-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,4) = i - 1
                      ! j-index of 2nd neighbor
                      patches(iZone)%BCnodes(i,j,5) = j - 1
                      ! i-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,6) = i
                      ! j-index of 3rd neighbor
                      patches(iZone)%BCnodes(i,j,7) = j - 1
                   end if

                end if

             end if

          else if (bocoType .eq. BCSymmetryPlane) then ! We have a symmetry BC

             ! We need two neighbors in the interior of the mesh, so we follow this
             ! ordering:
             !
             !     X----2---X
             !     |    |   |
             !     X----1---X
             !     |    |   |
             !   --X--node--X--
             !
             ! The third neighbor will be zero.
             ! No special treatment for corners is necessary in this case

             if (jMinBC .eq. jMaxBC) then ! Only i changes (j-edge BC)
                j = jMinBC
                if (j .eq. 1) then ! We are at lower j bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do i = iMinBC,iMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j + 1
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = i
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = j + 2
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                else ! We are at upper j bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do i = iMinBC,iMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j - 1
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = i
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = j - 2
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                end if
                
             else if (iMinBC .eq. iMaxBC) then ! Only j changes (i-edge BC)
                i = iMinBC
                if (i .eq. 1) then ! We are at lower i bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do j = jMinBC,jMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i + 1
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = i + 2
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = j
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                else ! We are at upper i bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do j = jMinBC,jMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i - 1
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = i - 2
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = j
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do

                end if

             end if

          else if (bocoType .eq. BCWall) then ! We have a wall BC

             ! We need one neighbor in the interior of the mesh, so we follow this
             ! ordering:
             !
             !     X----X---X
             !     |    |   |
             !     X----1---X
             !     |    |   |
             !   --X--node--X--
             !
             ! The second and third neighbors will be zero.
             ! No special treatment for corners is necessary in this case

             if (jMinBC .eq. jMaxBC) then ! Only i changes (j-edge BC)
                j = jMinBC
                if (j .eq. 1) then ! We are at lower j bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do i = iMinBC,iMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j + 1
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = 0
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = 0
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                else ! We are at upper j bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do i = iMinBC,iMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j - 1
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = 0
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = 0
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                end if
                
             else if (iMinBC .eq. iMaxBC) then ! Only j changes (i-edge BC)
                i = iMinBC
                if (i .eq. 1) then ! We are at lower i bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do j = jMinBC,jMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i + 1
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = 0
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = 0
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do
                   
                else ! We are at upper i bound
                   
                   ! ALL NODES (no special treatment for corners is needed)
                   do j = jMinBC,jMaxBC
                      ! Check the BC type and change it if the new BC has higher
                      ! priority (higher index)
                      if (patches(iZone)%BCnodes(i,j,1) .lt. bocoIndex) then
                         ! Update BC type
                         patches(iZone)%BCnodes(i,j,1) = bocoIndex
                         ! Update neighbors
                         ! i-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,2) = i - 1
                         ! j-index of 1st neighbor
                         patches(iZone)%BCnodes(i,j,3) = j
                         ! i-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,4) = 0
                         ! j-index of 2nd neighbor
                         patches(iZone)%BCnodes(i,j,5) = 0
                         ! i-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,6) = 0
                         ! j-index of 3rd neighbor
                         patches(iZone)%BCnodes(i,j,7) = 0
                      end if
                   end do

                end if

             end if

          end if neighbor_case
          
       end do

       print *,'    Assigned Zone ',iZone,' with dimensions ',dims(1:2)
    end do

    ! Close file
    call cg_close_f(cg, ierr)
    if (ierr .eq. CG_ERROR) call cg_error_exit_f

  end subroutine readGrid

end module readCGNS


Module kd_tree
  use precision
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric works so far.

  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 
  !

  ! .. Parameters ..
  ! you choose this.
  Integer(kind=IntType), Parameter :: bucket_size = 50
  ! ..
  ! .. Derived Type Declarations ..
  ! Global information about the tree
  ! pointer to the actual
  ! data array
  ! dimensionality and total # of points
  ! permuted index into the
  ! data, so that
  ! indexes[l..u]
  ! of some bucket represent the indexes of the actual
  ! points in that bucket.
  ! root pointer of the tree
  ! an internal tree node
  ! the dimension to cut
  ! where to cut the dimension
  ! indices of points included in this node,
  ! referring back to indices
  ! child pointers
  ! One of these is created for each search.
  ! best squared distance found so far
  ! best index found so far
  ! Query vector
  ! indexes of best found so far
  ! squared distances found so far
  ! how many best distances we are searching
  ! for
  ! how many have been found so far, i.e. il(1:nfound)
  ! are the only valid indexes.
  ! exclude points within
  ! 'correltime'
  ! of 'centeridx'
  Type :: tree_master_record
      real(kind=realType), Dimension (:, :), Pointer :: the_data
      Integer(kind=IntType) :: dim, n
      Integer(kind=IntType), Dimension (:), Pointer :: indexes
      Type (tree_node), Pointer :: root
  End Type tree_master_record

  Type :: tree_node
      Integer(kind=IntType) :: dnum
      real(kind=realType) :: val
      Integer(kind=IntType) :: l, u
      Type (tree_node), Pointer :: left, right
  End Type tree_node

  Type :: tree_search_record
      Private
      real(kind=realType) :: bsd
      Integer(kind=IntType) :: bestind
      real(kind=realType), Dimension (:), Pointer :: qv
      Integer(kind=IntType), Dimension (:), Pointer :: il
      real(kind=realType), Dimension (:), Pointer :: dsl
      Integer(kind=IntType) :: nbst, nfound
      Integer(kind=IntType) :: centeridx, correltime
  End Type tree_search_record
  ! ..
  ! set to true if we're doing linfinity metric
  !
  
  logical, parameter :: vectorizing = .true.
  !
  ! set to .true. to use optimization designed for a vectorizing machine
  ! for example, the Intel ifc compiler with -axK or -axW options to use
  ! the special instructions.
  !
  ! Even still it doesn't always improve performance. :(
  !

Contains

  Subroutine destroy_tree(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    Call destroy_node(tp%root)
    Deallocate (tp%indexes)
    Nullify (tp%indexes)
    Return

  Contains

    Recursive Subroutine destroy_node(np)
      ! .. Structure Arguments ..
      Type (tree_node), Pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      Intrinsic ASSOCIATED
      ! ..
      If (ASSOCIATED(np%left)) Then
         Call destroy_node(np%left)
         Deallocate (np%left)
         Nullify (np%left)
      End If
      If (ASSOCIATED(np%right)) Then
         Call destroy_node(np%right)
         Deallocate (np%right)
         Nullify (np%right)
      End If
      Return
    End Subroutine destroy_node
  End Subroutine destroy_tree

  Function create_tree(input_data) Result (master_record)
    ! create the actual tree structure, given an input array of data.
    ! Arguments
    ! .. Function Return Value ..
    Type (tree_master_record), Pointer :: master_record
    ! ..
    ! .. Array Arguments ..
    real(kind=realType), Target :: input_data(:, :)
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic SIZE
    ! ..
    Allocate (master_record)
    master_record%the_data => input_data
    ! pointer assignment
    master_record%dim = SIZE(input_data, 1)
    master_record%n = SIZE(input_data, 2)

!    Print *, 'Creating KD tree with N = ', master_record%n, &
!     ' and dim = ', master_record%dim
    Call build_tree(master_record)

  Contains

    Subroutine build_tree(tp)
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Local Scalars ..
      Integer(kind=IntType) :: j
      ! ..
      Allocate (tp%indexes(tp%n))
      forall (j=1:tp%n)
         tp%indexes(j) = j
      end forall
!      Do j = 1, tp%n
!         tp%indexes(j) = j
!      End Do
      tp%root => build_tree_for_range(tp, 1, tp%n)
    End Subroutine build_tree

    Recursive Function build_tree_for_range(tp, l, u) Result (res)
      ! .. Function Return Value ..
      Type (tree_node), Pointer :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      Integer(kind=IntType) :: c, m
      ! ..
      if (.false.) then 
         if ((l .lt. 1) .or. (l .gt. tp%n)) then
            stop 'illegal L value in build_tree_for_range'
         end if
         if ((u .lt. 1) .or. (u .gt. tp%n)) then
            stop 'illegal u value in build_tree_for_range'
         end if
         if (u .lt. l) then
            stop 'U is less than L, thats illegal.'
         end if
      endif
      If ((u-l)<=bucket_size) Then
         Allocate (res)
         res%dnum = 0
         res%val = 0.0
         res%l = l
         res%u = u
         Nullify (res%left, res%right)
      Else
         Allocate (res)
         c = most_spread_coordinate(tp, l, u)
         m = (l+u)/2
         Call select_on_coordinate(tp%the_data, tp%indexes, c, m, l, u)
         ! moves indexes around
         res%dnum = c
         res%val = tp%the_data(c, tp%indexes(m))
         res%l = l
         res%u = u
         res%left => build_tree_for_range(tp, l, m)
         res%right => build_tree_for_range(tp, m+1, u)
      End If
    End Function build_tree_for_range

    Subroutine select_on_coordinate(v, ind, c, k, li, ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Structure Arguments ..
!      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: c, k, li, ui
      ! ..
      ! .. Local Scalars ..
      Integer(kind=IntType) :: i, l, m, s, t, u
      ! ..
      ! .. Local Arrays ..
      real(kind=realType) :: v(:, :)
      Integer(kind=IntType) :: ind(:)
      ! ..
      !v => tp%the_data
      !ind => tp%indexes
      l = li
      u = ui
      Do While (l<u)
         t = ind(l)
         m = l
         Do i = l + 1, u
            If (v(c, ind(i))<v(c, t)) Then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            End If
         End Do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         If (m<=k) l = m + 1
         If (m>=k) u = m - 1
      End Do
    End Subroutine select_on_coordinate

    Function most_spread_coordinate(tp, l, u) Result (res)
      ! Of indices in l..u find the axis which has the largest spread, 
      ! and
      ! return its index
      ! .. Function Return Value ..
      Integer(kind=IntType) :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      real(kind=realType) :: bsp, sp
      Integer(kind=IntType) :: i = 0
      ! ..
      res = 0
      bsp = -1.0
      Do i = 1, tp%dim
         sp = spread_in_coordinate(tp, i, l, u)
         If (sp>bsp) Then
            res = i
            bsp = sp
         End If
      End Do
    End Function most_spread_coordinate

    Function spread_in_coordinate(tp, c, l, u) Result (res)
      ! the spread in coordinate 'c', between l and u
      ! for easier local access
      ! ibid
      ! .. Function Return Value ..
      real(kind=realType) :: res
      ! ..
      ! .. Structure Arguments ..
      Type (tree_master_record), Pointer :: tp
      ! ..
      ! .. Scalar Arguments ..
      Integer(kind=IntType), Intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(kind=realType) :: last, lmax, lmin, smax, smin, t
      Integer(kind=IntType) :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(kind=realType), Pointer :: v(:, :)
      Integer(kind=IntType), Pointer :: ind(:)
      ! ..
      v => tp%the_data
      ind => tp%indexes
      smin = v(c, ind(l))
      smax = smin
      
      ! truncate u here or not?????
      ulocal = min(u, l+200) !!!?????

      Do i = l + 2, ulocal, 2
         lmin = v(c, ind(i-1))
         lmax = v(c, ind(i))
         If (lmin>lmax) Then
            t = lmin
            lmin = lmax
            lmax = t
         End If
         If (smin>lmin) smin = lmin
         If (smax<lmax) smax = lmax
      End Do
      If (i==ulocal+1) Then
         last = v(c, ind(ulocal))
         If (smin>last) smin = last
         If (smax<last) smax = last
      End If
      res = smax - smin
      ! write *, "Returning spread in coordinate with res = ", res
    End Function spread_in_coordinate
  End Function create_tree

  ! Search routines:

  ! * n_nearest_to(tp, qv, n, indexes, distances)

  ! Find the 'n' vectors in the tree nearest to 'qv' in euclidean norm
  ! returning their indexes and distances in 'indexes' and 'distances'
  ! arrays already allocated passed to the subroutine.

  
  Subroutine n_nearest_to(tp, qv, n, indexes, distances)
    ! find the 'n' nearest neighbors to 'qv', returning
    ! their indexes in indexes and squared Euclidean distances in
    ! 'distances'
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer(kind=IntType), Intent (In) :: n
    ! ..
    ! .. Array Arguments ..
    real(kind=realType), Target :: distances(n)
    real(kind=realType), Target, Intent (In) :: qv(:)
    Integer(kind=IntType), Target :: indexes(n)
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic HUGE
    ! ..
    ! the largest real number
    sr%bsd = HUGE(1.0)
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr


    Call n_search(tp, psr, tp%root)
    Return
  End Subroutine n_nearest_to

  ! Another main routine you call from your user program
  Subroutine n_nearest_to_around_point(tp, idxin, correltime, n, indexes, &
   distances)
    ! find the 'n' nearest neighbors to point 'idxin', returning
    ! their indexes in indexes and squared Euclidean distances in
    ! 'distances'
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! .. 
    ! .. Scalar Arguments ..
    Integer, Intent (In) :: correltime, idxin, n
    ! .. 
    ! .. Array Arguments ..
    Real, Target :: distances(n)
    Integer, Target :: indexes(n)
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Local Arrays ..
    Real, Allocatable, Target :: qv(:)
    ! ..
    ! .. Intrinsic Functions ..                                                                                   
    Intrinsic HUGE
    ! ..
    Allocate (qv(tp%dim))
    qv = tp%the_data(:, idxin) ! copy the vector
    sr%bsd = HUGE(1.0)       ! the largest real number
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr 

    Call n_search(tp, psr, tp%root)
    Deallocate (qv)
    Return
  End Subroutine n_nearest_to_around_point

  function pt_in_tree(tp, qv)
 ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    real(kind=realType), Target, Intent (In) :: qv(:)
    real(kind=realType), parameter :: tol=1e-6
  
    ! .. Scalar Arguments ..
    Integer(kind=IntType) :: n
    real(kind=realType), target :: distances(1)
    Integer(kind=IntType), target :: indexes(1)
    integer(kind=intType) :: pt_in_tree
    ! ..
    ! .. Local Structures ..
    Type (tree_search_record), Pointer :: psr
    Type (tree_search_record), Target :: sr
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic HUGE
    n = 1
    ! ..
    ! the largest real number
    sr%bsd = HUGE(1.0)
    sr%qv => qv
    sr%nbst = n
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%dsl => distances
    sr%il => indexes
    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    sr%il = -1               ! set to invalid indexes
    psr => sr                ! in C this would be psr = &sr

    Call n_search(tp, psr, tp%root)
    
    if (sqrt(distances(1)) < tol) then
       pt_in_tree = indexes(1)
    else
       pt_in_tree = 0_intType
    end if

  end function pt_in_tree

  Function square_distance(n, iv, qv) Result (res)
    ! distance between v[i, *] and qv[*] 
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(kind=realType) :: res
    ! ..
    ! ..
    ! .. Scalar Arguments ..
    Integer(kind=IntType) :: n
    ! ..
    ! .. Array Arguments ..
    real(kind=realType) :: iv(:), qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:n)-qv(1:n))**2 )
  End Function square_distance

  Recursive Subroutine n_search(tp, sr, node)
    ! This is the innermost core routine of the kd-tree search.
    ! it is thus not programmed in quite the clearest style, but is
    ! designed for speed.  -mbk
    ! .. Structure Arguments ..
    Type (tree_node), Pointer :: node
    Type (tree_search_record), Pointer :: sr
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Local Scalars ..
    real(kind=realType) :: dp, sd, sdp, tmp
    Integer(kind=IntType) :: centeridx, i, ii, j, jmax, k, d, correltime, nbst
    Logical :: not_fully_sized
    ! ..
    ! .. Local Arrays ..
    real(kind=realType), Pointer :: qv(:)
    Integer(kind=IntType), Pointer :: ind(:)
    real(kind=realType), dimension(:, :), pointer :: data
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic ABS, ASSOCIATED
    ! ..
    If ( .Not. (ASSOCIATED(node%left)) .And. ( .Not. ASSOCIATED( &
     node%right))) Then
       ! we are on a terminal node
       ind => tp%indexes     ! save for easy access
       qv => sr%qv
       data => tp%the_data
       centeridx = sr%centeridx
       d = tp%dim
       correltime = sr%correltime
       ! search through terminal bucket.
       nbst = sr%nbst
       mainloop: Do i = node%l, node%u
          ii = ind(i)
          If ( (centeridx<0) .OR. (ABS(ii-centeridx)>=correltime)) Then
             ! 
             ! replace with call to square distance with inline
             ! code, an
             ! a goto.  Tested to be significantly faster.   SPECIFIC
             ! FOR
             ! the EUCLIDEAN METRIC ONLY! BEWARE!

             sd = 0.0
             if (.not. vectorizing)  then
                Do k = 1, d
                   ! comment out??
                   !                      Write(*, *) k
                   !                      Write(*, *) qv(k), ' ii = ', ii
                   !                      Write(*, *) tp%the_data(ii, k), qv(k)
                   ! comment out
                   tmp = data(k, ii) - qv(k)
                   sd = sd + tmp*tmp
                   If (sd>sr%bsd) Then
                      cycle mainloop
                   End If
                End Do
             else
2               sd = square_distance(tp%dim, data(:, ii), qv)
                If (sd > sr%bsd) Then
                   cycle mainloop
                End If
             endif
             ! Note test moved out of loop to improve vectorization
             ! should be semantically identical eitehr way as sr%bsd is
             ! a bound afterwhich it doesn't matter. 

             ! we only consider it if it is better than the 'best' on
             ! the list so far.
             ! if we get here
             ! we know sd is < bsd, and bsd is on the end of the list

             ! special case for nbst = 1 (single nearest neighbor)
             if (nbst .eq. 1) then
                sr%il(1) = ii
                sr%dsl(1) = sd
                sr%bsd = sd
             else
                not_fully_sized = (sr%nfound<nbst)
                If (not_fully_sized) Then
                   jmax = sr%nfound
                   sr%nfound = sr%nfound + 1
                Else
                   jmax = nbst - 1
                End If
                ! add it to the list
                
                ! find the location j where sd >= sr%dsl(j) and sd <
                ! sr%dsl(j+1...jmax+1)
                if (vectorizing) then
                   do j=jmax, 1, -1
                      if (sd>=sr%dsl(j)) Exit
                   end do
                   sr%il(j+2:jmax+1) = sr%il(j+1:jmax)
                   sr%dsl(j+2:jmax+1) = sr%dsl(j+1:jmax)
                   sr%il(j+1) = ii
                   sr%dsl(j+1) = sd
                else
                   Do j = jmax, 1, -1
                      If (sd>=sr%dsl(j)) Exit ! we hit insertion location
                      sr%il(j+1) = sr%il(j)
                      sr%dsl(j+1) = sr%dsl(j)
                   End Do
                   ! if loop falls through j=0 here.
                   sr%il(j+1) = ii
                   sr%dsl(j+1) = sd
                endif
                If ( .Not. not_fully_sized) Then
                   sr%bsd = sr%dsl(nbst)
                End If
             end if
          End If
       End Do mainloop
    Else
       ! we are not on a terminal node

       ! Alrighty, this section is essentially the content of the
       ! the Sproul method for searching a Kd-tree, in other words
       ! the second nested "if" statements in the two halves below
       ! and when they get activated.
       dp = sr%qv(node%dnum) - node%val
       sdp = dp*dp        ! Euclidean
       If (dp<0.0) Then
          Call n_search(tp, sr, node%left)
          If (sdp<sr%bsd) Call n_search(tp, sr, node%right)
          ! if the distance projected to the 'wall boundary' is less
          ! than the radius of the ball we must consider, then perform
          ! search on the 'other side' as well.
       Else
          Call n_search(tp, sr, node%right)
          If (sdp<sr%bsd) Call n_search(tp, sr, node%left)
       End If
    End If
  End Subroutine n_search

  Subroutine n_nearest_to_brute_force(tp, qv, n, indexes, distances)
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    ! .. Structure Arguments ..
    Type (tree_master_record), Pointer :: tp
    ! ..
    ! .. Scalar Arguments ..
    Integer(kind=IntType), Intent (In) :: n
    ! ..
    ! .. Array Arguments ..
    real(kind=realType), intent(out) :: distances(n)
    real(kind=realType), Intent (In) :: qv(:)
    Integer(kind=IntType), intent(out) :: indexes(n)
    ! ..
    ! .. Local Scalars ..
    Integer(kind=IntType) :: i, j, k
    ! ..
    ! .. Local Arrays ..
    real(kind=realType), Allocatable :: all_distances(:)
    ! ..
    ! .. Intrinsic Functions ..
    Intrinsic HUGE
    ! ..
    Allocate (all_distances(tp%n))
    Do i = 1, tp%n
       all_distances(i) = square_distance(tp%dim, tp%the_data(:, i), qv)
    End Do
    ! now find 'n' smallest distances
    distances(1:n) = HUGE(1.0)
    indexes(1:n) = -1 
!    Do i = 1, n
!       distances(i) = HUGE(1.0)
!       indexes(i) = -1
!    End Do
    Do i = 1, tp%n
       If (all_distances(i)<distances(n)) Then
          ! insert it somewhere on the list
          Do j = 1, n
             If (all_distances(i)<distances(j)) Exit
          End Do
          ! now we know 'j'
          if (.not. vectorizing) then
             Do k = n - 1, j, -1
                distances(k+1) = distances(k)
                indexes(k+1) = indexes(k)
             End Do
          else
             distances(j+1:n) = distances(j:n-1)
             indexes(j+1:n) = indexes(j:n)
          end if
          distances(j) = all_distances(i)
          indexes(j) = i
       End If
    End Do
    Deallocate (all_distances)
  End Subroutine n_nearest_to_brute_force

End Module kd_tree

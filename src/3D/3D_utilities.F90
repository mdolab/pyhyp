subroutine computeStretch(L)
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: computeStretch determines the global stretch value,
    !     deltaS, depending on the grid level L. In the future
    !     this function may be more complex or call a user-supplied function.
    !
    !     Parameters
    !     ----------
    !     L : integer
    !         The marching level to compute stretch for
    !
    !     Returns
    !     -------
    !     deltaS : real
    !         Marching increment set in hypData
    !
    use precision
    use hypInput
    use hypData
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: L

    ! Simply return the pre-computed distances
    desiredDeltaS = fullDeltaS(L)

end subroutine computeStretch

subroutine calcGridRatio(N, nStart, nEnd, s0, S, ratio)
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: calcGridRatio() calculates the exponential
    !     distribution Turns out we need to solve a transcendental
    !     equation here. We will do this with a bisection search
    !
    !     Parameters
    !     ----------
    !     N : integer
    !         The number of nodes in sequence
    !     nStart : integer
    !         The number of intervals with constant spacing
    !         's0' at the beginning
    !     nEnd : integer
    !         The number of constant sized intervals at the end.
    !         The actual spacing will be determined.
    !     s0 : real
    !         The initial grid spacing
    !     S : real
    !         The total integrated length
    !
    !     Returns
    !     -------
    !     ratio : real
    !         The computed grid ratio that satifies all the inputs.
    use hypInput, only: fullDeltas, growthRatios
    use precision

    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: N, nStart, nEnd
    real(kind=realType), intent(in) :: s0, S

    ! Output Parameters
    real(kind=realType), intent(out) :: ratio

    ! Working Parameters
    integer(kind=intType) :: i, j
    real(kind=realType) :: r, a, b, c, f, fa, fb, curSize

    ! allocate the array that holds delta S values
    if (.not. allocated(fullDeltaS)) then
        allocate (fullDeltaS(2:N))
    end if

    ! if user provided a growthRatio array, compute delta S values based on this
    if (allocated(growthRatios)) then

        curSize = s0
        do j = 1, N - 1
            curSize = curSize * growthRatios(j)
            fullDeltaS(j + 1) = curSize
        end do

        ratio = 1.05

    else ! compute fullDeltaS based on input parameters

        ! function 'f' is S - s0*(1-r^n)/(1-r) where S is total length, s0 is
        ! initial ratio and r is the grid ratio.

        ! Do a bisection search
        ! Max and min bounds...root must be in here...
        a = one + 1e-8
        b = four

        fa = func(a)
        fb = func(b)
        do i = 1, 100
            c = half * (a + b)
            f = func(c)
            if (abs(f) < 1e-10) then ! Converged
                exit
            end if

            if (f * fa > 0) then
                a = c
            else
                b = c
            end if
        end do

        ! Finally set the ratio variable to r
        ratio = c

        ! And we precompute all stretches:
        curSize = s0
        do j = 1, nStart
            fullDeltaS(j + 1) = curSize
        end do

        ! Next we have N - nStart - nEnd layers of exponential
        do j = 1, N - 1 - nStart - nEnd
            curSize = curSize * ratio
            fullDeltaS(nStart + j + 1) = curSize
        end do

        ! Finally we have the last nEnd constant layers
        curSize = curSize * ratio
        do j = 1, nEnd
            fullDeltaS(N - nEnd + j) = curSize
        end do

    end if

contains
    function func(r)

        ! Evaluate the function we want to solve for:
        real(kind=realType) :: r, func, curSize
        integer(kind=intType) :: j

        ! We will have nStart layers at the beginning.
        func = nStart * s0
        curSize = s0

        ! Next we will have M = N - nStart - nEnd layers of exponential
        ! growth.
        do j = 1, N - 1 - nStart - nEnd
            curSize = curSize * r
            func = func + curSize
        end do

        ! Last stretch
        curSize = curSize * r

        ! Now add the last nEnd layers of constant size
        func = func + nEnd * curSize

        ! Finally the actual functio is S - func
        func = S - func

    end function func

end subroutine calcGridRatio

subroutine three_by_three_inverse(Jac, Jinv)
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: return the inverse of the 3x3 matrix Jac in Jinv
    !
    !     Parameters
    !     ----------
    !     Jac : array size (3,3)
    !         Matrix to invert
    !
    !     Returns
    !     -------
    !     Jinv : array size (3,3)
    !         Inverse of Jac
    !
    use precision
    implicit none

    real(kind=realType), intent(in) :: Jac(3, 3)
    real(kind=realType), intent(out) :: Jinv(3, 3)
    real(kind=realType) :: invdet, det

    ! Express the inverse of the jacobian explicitly
    det = Jac(1, 1) * Jac(2, 2) * Jac(3, 3) &
          + Jac(1, 2) * Jac(2, 3) * Jac(3, 1) &
          + Jac(1, 3) * Jac(2, 1) * Jac(3, 2) &
          - Jac(1, 1) * Jac(2, 3) * Jac(3, 2) &
          - Jac(1, 2) * Jac(2, 1) * Jac(3, 3) &
          - Jac(1, 3) * Jac(2, 2) * Jac(3, 1)
    invdet = one / det

    !             [ |a22 a23|   |a12 a13|  |a12 a13|]     */
    !             [ |a32 a33|  -|a32 a33|  |a22 a23|]     */
    !             [                                 ]     */
    !             [ |a21 a23|   |a11 a13|  |a11 a13|]     */
    !    A^(-1) = [-|a31 a33|   |a31 a33| -|a21 a23|] / d */
    !             [                                 ]     */
    !             [ |a21 a22|   |a11 a12|  |a11 a12|]     */
    !             [ |a31 a32|  -|a31 a32|  |a21 a22|]     */

    Jinv(1, 1) = invdet * (Jac(2, 2) * Jac(3, 3) - Jac(2, 3) * Jac(3, 2))
    Jinv(1, 2) = -invdet * (Jac(1, 2) * Jac(3, 3) - Jac(1, 3) * Jac(3, 2))
    Jinv(1, 3) = invdet * (Jac(1, 2) * Jac(2, 3) - Jac(1, 3) * Jac(2, 2))

    Jinv(2, 1) = -invdet * (Jac(2, 1) * Jac(3, 3) - Jac(2, 3) * Jac(3, 1))
    Jinv(2, 2) = invdet * (Jac(1, 1) * Jac(3, 3) - Jac(1, 3) * Jac(3, 1))
    Jinv(2, 3) = -invdet * (Jac(1, 1) * Jac(2, 3) - Jac(1, 3) * Jac(2, 1))

    Jinv(3, 1) = invdet * (Jac(2, 1) * Jac(3, 2) - Jac(2, 2) * Jac(3, 1))
    Jinv(3, 2) = -invdet * (Jac(1, 1) * Jac(3, 2) - Jac(1, 2) * Jac(3, 1))
    Jinv(3, 3) = invdet * (Jac(1, 1) * Jac(2, 2) - Jac(1, 2) * Jac(2, 1))

end subroutine three_by_three_inverse

subroutine cross_prod(a, b, c)

    ! Written by Ney Secco, 2015
    ! This function computes the cross-product between a and b

    use precision

    ! Inputs
    real(kind=realType), dimension(3), intent(in) :: a, b

    ! Outputs
    real(kind=realType), dimension(3), intent(out) :: c

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

end subroutine cross_prod

subroutine writeHeader
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Writes the header to stdout for the output
    !
    use precision
    implicit none

    ! Working
    integer(kind=intType) :: i

    ! Write the first line of '-'
    write (*, "(a)", advance="no") "#"
    do i = 2, 116
        write (*, "(a)", advance="no") "-"
    end do
    write (*, "(a)", advance="no") "#"
    print "(1x)"
write(*,"(a)",advance="no") "# Grid | CPU  | Sub | KSP | nAvg |  Sl  | Sensor | Sensor | Min     | Min     |  deltaS  |"
    write (*, "(a)", advance="no") " March    | cMax  | Ratio |"
    print "(1x)"
write(*,"(a)",advance="no") "# Lvl  | Time | Its | Its |      |      | Max    | Min    | Quality | Volume  |          |"
    write (*, "(a)", advance="no") " Distance |       | kMax  |"
    print "(1x)"

    ! Write the Last line of '-'
    write (*, "(a)", advance="no") "#"
    do i = 2, 116
        write (*, "(a)", advance="no") "-"
    end do
    write (*, "(a)", advance="no") "#"
    print "(1x)"

end subroutine writeHeader

subroutine writeIteration
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Writes the information pertaining to the curent iteration to
    !     the screen
    !
    use communication
    use hypData

    implicit none

    ! Working variables
    real(kind=realType) :: gridSensorMaxRed, gridSensorMinRed, maxKStretchRed
    integer(kind=intType) :: ierr, nAverageRed

    ! Reduce all mins/maxes
    CALL MPI_REDUCE(gridSensorMax, gridSensorMaxRed, 1, MPI_DOUBLE, MPI_MAX, 0, hyp_comm_world, ierr)
    CALL MPI_REDUCE(gridSensorMin, gridSensorMinRed, 1, MPI_DOUBLE, MPI_MIN, 0, hyp_comm_world, ierr)
    CALL MPI_REDUCE(maxKStretch, maxKStretchRed, 1, MPI_DOUBLE, MPI_MAX, 0, hyp_comm_world, ierr)
    CALL MPI_REDUCE(nAverage, nAverageRed, 1, MPI_INTEGER, MPI_SUM, 0, hyp_comm_world, ierr)

    ! Only root proc prints stuff
    if (myid == 0) then

        ! Iteration Count
        write (*, "(i7,1x)", advance="no") marchIter

        ! CPU time
        write (*, "(f6.1,1x)", advance="no") dble(mpi_wtime()) - timeStart

        ! Sub Iteration
        write (*, "(i5,1x)", advance="no") nSubIter

        ! KSP ITerations
        write (*, "(i5,1x)", advance="no") kspIts

        ! Number of averaged nodes
        write (*, "(i6,1x)", advance="no") nAverageRed

        ! Modified Sl factor from Steger and Chan
        write (*, "(f6.3,1x)", advance="no") Sl

        ! maximum Grid convergence sensor
        write (*, "(f8.5,1x)", advance="no") gridSensorMaxRed

        ! minimum Grid convergence sensor
        write (*, "(f8.5,1x)", advance="no") gridSensorMinRed

        ! minimum grid quality for new layer
        write (*, "(f8.5,1x)", advance="no") minQuality ! This was reduced already

        ! minimum volume for new layer
        write (*, "(e10.3,1x)", advance="no") minVolume ! This was reduced already

        ! marching step size for this iteration
        write (*, "(e10.3,1x)", advance="no") deltaS

        ! March Distance
        write (*, "(e10.3,1x)", advance="no") scaleDist

        ! marching step size for this iteration
        write (*, "(f7.4,1x)", advance="no") cRatio

        ! maximum stretch ratio in K-direction
        write (*, "(f7.4,1x)", advance="no") maxKStretchRed

        ! Jump to next line
        print "(1x)"

    end if
end subroutine writeIteration

subroutine computeQualityLayer
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Compute the minimum quality and volume measure for the cells
    !     defined by two grid levels
    !
    use communication
    use hypInput
    use hypData
    implicit none

    ! Working
    integer(kind=intType) :: i, j, iPatch, ierr
    real(kind=realType) :: points(3, 8), Q, V
    real(Kind=realType) :: minQuality_local, minVolume_local
    ! Since we are dealing with volumes, we need to loop over faces:
    minQuality_local = one
    minVolume_local = huge(one)
    call VecGhostUpdateBegin(X(marchIter), INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGhostUpdateEnd(X(marchIter), INSERT_VALUES, SCATTER_FORWARD, ierr)

    call VecGhostGetLocalForm(X(marchIter), XL_local, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGhostGetLocalForm(X(marchIter - 1), XLm1_local, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(XL_local, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(XLm1_local, xxm1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    do i = 1, nLocalFace

        ! Assemble the 8 points we need for the volume - Note
        ! coodinate coordinate ordering!
        points(:, 1) = xxm1(3 * conn(1, i) - 2:3 * conn(1, i))
        points(:, 2) = xxm1(3 * conn(2, i) - 2:3 * conn(2, i))
        points(:, 3) = xxm1(3 * conn(4, i) - 2:3 * conn(4, i))
        points(:, 4) = xxm1(3 * conn(3, i) - 2:3 * conn(3, i))

        points(:, 5) = xx(3 * conn(1, i) - 2:3 * conn(1, i))
        points(:, 6) = xx(3 * conn(2, i) - 2:3 * conn(2, i))
        points(:, 7) = xx(3 * conn(4, i) - 2:3 * conn(4, i))
        points(:, 8) = xx(3 * conn(3, i) - 2:3 * conn(3, i))

        call quality_hexa(points, Q)
        minQuality_local = min(minQuality_local, Q)

        call volume_hexa(points, V)
        minVolume_local = min(minVolume_local, V)
    end do

    call VecRestoreArrayF90(XL_local, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(XLm1_local, xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGhostRestoreLocalForm(X(marchIter), XL_local, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGhostRestoreLocalForm(X(marchIter - 1), XLm1_local, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Need to reduce the min value to the root proc
    call mpi_Reduce(minQuality_local, minQuality, 1, MPI_DOUBLE, MPI_MIN, 0, hyp_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Need to reduce the min value to the root proc
    call mpi_Reduce(minVolume_local, minVolume, 1, MPI_DOUBLE, MPI_MIN, 0, hyp_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Check if these values with the overall minima
    minQualityOverall = min(minQuality, minQualityOverall)
    minVolumeOverall = min(minVolume, minVolumeOverall)

end subroutine computeQualityLayer

subroutine volume_hexa(points, V)
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: volume_hexa determines the volume for a hexa defined
    !     by 8 nodes.
    !
    !     Parameters
    !     ----------
    !     points : real array size(3, 8)
    !         Nodes defining the hexa in coordinate ordering.
    !
    !     Returns
    !     -------
    !     Volume : real
    !         The computed cell volume
    use precision
    implicit none

    ! Input/Output
    real(kind=realType), intent(in) :: points(3, 8)
    real(kind=realType), intent(out) :: V

    ! Working
    real(kind=realType) :: center(3), volpymrid, v1, v2, v3, v4, v5, v6
    integer(kind=intType) :: i, idim

    ! Compute the center of the points
    center = zero
    do i = 1, 8
        do idim = 1, 3
            center(idim) = center(idim) + points(idim, i)
        end do
    end do
    center = center / eight

    ! Compute the volumes of the 6 sub pyramids. The
    ! arguments of volpym must be such that for a (regular)
    ! right handed hexahedron all volumes are positive.

    v1 = volpymrid(center, points(:, 1), points(:, 2), points(:, 4), points(:, 3))
    v2 = volpymrid(center, points(:, 7), points(:, 8), points(:, 6), points(:, 5))
    v3 = volpymrid(center, points(:, 1), points(:, 3), points(:, 7), points(:, 5))
    v4 = volpymrid(center, points(:, 4), points(:, 2), points(:, 6), points(:, 8))
    v5 = volpymrid(center, points(:, 2), points(:, 1), points(:, 5), points(:, 6))
    v6 = volpymrid(center, points(:, 3), points(:, 4), points(:, 8), points(:, 7))

    V = (v1 + v2 + v3 + v4 + v5 + v6) / six
end subroutine volume_hexa

function volpymrid(p, a, b, c, d)
    use precision
    implicit none
    real(kind=realType) :: p(3), a(3), b(3), c(3), d(3), volpymrid

    ! 6*Volume of a pyrimid -> Counter clockwise ordering
    volpymrid = (p(1) - fourth * (a(1) + b(1) + c(1) + d(1))) * &
                ((a(2) - c(2)) * (b(3) - d(3)) - (a(3) - c(3)) * (b(2) - d(2))) + &
                (p(2) - fourth * (a(2) + b(2) + c(2) + d(2))) * &
                ((a(3) - c(3)) * (b(1) - d(1)) - (a(1) - c(1)) * (b(3) - d(3))) + &
                (p(3) - fourth * (a(3) + b(3) + c(3) + d(3))) * &
                ((a(1) - c(1)) * (b(2) - d(2)) - (a(2) - c(2)) * (b(1) - d(1)))

end function volpymrid

subroutine quality_hexa(points, quality)
    !***DESCRIPTION
    !
    !     Written by Gaetan Kenway
    !
    !     Abstract: quality_hexa determine the 2x2x2 quality measure for
    !     a hexa defined by 8 nodes.
    !
    !     Parameters
    !     ----------
    !     points : real array size(3, 8)
    !         Nodes defining the hexa in coordinate ording.
    !
    !     Returns
    !     -------
    !     Quality : real
    !         The computed quality measure
    use precision
    implicit none

    ! Input/Output
    real(kind=realType), intent(in) :: points(3, 8)
    real(kind=realType), intent(out) :: quality

    ! Working
    integer(kind=intType), parameter :: Ng = 2
    real(kind=realType) :: det(Ng * Ng * Ng)
    integer(kind=intType) :: i, j, k, l, m, nn, counter
    real(kind=realType) :: r, s, t
    real(kind=realType) :: largest, d_max, corners(Ng)
    real(kind=realType) :: Xd(9), N(8), Na(8), Nb(8), Nc(8), Ua(9)
    real(kind=realType) :: JJac(9), h, pt(3)
    real(kind=realType) :: shp(2, 2), dshp(2, 2)
    ! Relative Determinant: the ratio of the smallest determinat of the
    ! jacobian matrix divided by the largest determinate of the jacobian
    ! matrix using a 3x3x3 stecil

    corners(1) = -1
    corners(2) = 1
    call lagrangeSF(shp(:, 1), dshp(:, 1), corners(1), 2)
    call lagrangeSF(shp(:, 2), dshp(:, 2), corners(2), 2)

    ! 2x2x2 Stencil
    do nn = 1, Ng ! Gauss in r
        do m = 1, Ng ! Gauss in s
            do l = 1, Ng ! Gauss in t

                Xd = zero
                do k = 1, 2
                    do j = 1, 2
                        do i = 1, 2
                            counter = (k - 1) * 4 + (j - 1) * 2 + i
                            Na(counter) = dshp(i, l) * shp(j, m) * shp(k, nn)
                            Xd(1) = Xd(1) + points(1, counter) * Na(counter)
                            Xd(2) = Xd(2) + points(2, counter) * Na(counter)
                            Xd(3) = Xd(3) + points(3, counter) * Na(counter)

                            Nb(counter) = shp(i, l) * dshp(j, m) * shp(k, nn)
                            Xd(4) = Xd(4) + points(1, counter) * Nb(counter)
                            Xd(5) = Xd(5) + points(2, counter) * Nb(counter)
                            Xd(6) = Xd(6) + points(3, counter) * Nb(counter)

                            Nc(counter) = shp(i, l) * shp(j, m) * dshp(k, nn)
                            Xd(7) = Xd(7) + points(1, counter) * Nc(counter)
                            Xd(8) = Xd(8) + points(2, counter) * Nc(counter)
                            Xd(9) = Xd(9) + points(3, counter) * Nc(counter)

                        end do
                    end do
                end do

                h = (Xd(9) * (Xd(1) * Xd(5) - Xd(4) * Xd(2)) &
                     - Xd(8) * (Xd(1) * Xd(6) - Xd(4) * Xd(3)) &
                     + Xd(7) * (Xd(2) * Xd(6) - Xd(3) * Xd(5)))

                det((nn - 1) * Ng * Ng + (m - 1) * Ng + l) = h
            end do
        end do
    end do

    ! We don't use maxval since we can't cs it
    largest = maxval(abs(det))
    d_max = maxval(largest - det)

    quality = 1 - d_max / largest

end subroutine quality_hexa

subroutine lagrangeSF(sf, dsf, a, order)

    use precision
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: order
    real(kind=realType), intent(out) :: sf(order), dsf(order)
    real(kind=realType), intent(in) :: a

    ! Working
    integer(kind=intType) :: i, j, k
    real(kind=realType) :: ki, kj, dn, kk

    sf = one
    dsf = zero

    ! Loop over the shape function control points
    do i = 0, order - 1
        ki = -one + two * i / (order - one)

        ! Loop over each point again, except for the current control point,
        ! adding the contribution to the shape function

        do j = 0, order - 1
            if (i .ne. j) then
                kj = -one + two * j / (order - one)
                sf(i + 1) = sf(i + 1) * (a - kj) / (ki - kj)

                ! Loop over the whole thing again to determine the contribution to
                ! the derivative of the shape function
                dn = one / (ki - kj)

                do k = 0, order - 1
                    if (k .ne. i .and. k .ne. j) then
                        kk = -one + two * k / (order - one)
                        dn = dn * (a - kk) / (ki - kk)
                    end if
                end do

                dsf(i + 1) = dsf(i + 1) + dn

            end if
        end do
    end do

end subroutine lagrangeSF

function dist(p1, p2)
    use precision

    real(kind=realType) :: p1(3), p2(3)
    real(kind=realType) :: dist
    dist = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)

end function dist

subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)

    ! Given a list of N points (pts) in three space, with possible
    ! duplicates, (to within tol) return a list of the nUnique
    ! uniquePoints of points and a link array of length N, that points
    ! into the unique list
    use precision
    use kdtree2_module
    use hypdata
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: N
    real(kind=realType), intent(in), dimension(3, N) :: pts
    real(kind=realType), intent(in) :: tol

    ! Output Parametres
    real(kind=realType), intent(out), dimension(3, N) :: uniquePts
    integer(kind=intType), intent(out), dimension(N) :: link
    integer(kind=intType), intent(out) :: nUnique

    ! Working paramters
    type(kdtree2), pointer :: mytree
    real(kind=realType) :: tol2, timeb, timea
    integer(kind=intType) :: nFound, i, j, nAlloc
    type(kdtree2_result), allocatable, dimension(:) :: results

    ! We will use the KD_tree to do most of the heavy lifting here:
    mytree => kdtree2_create(pts, sort=.True.)

    ! KD tree works with the square of the tolerance
    tol2 = tol**2

    ! Unlikely we'll have more than 20 points same, but there is a
    ! safetly check anwyay.
    nalloc = 20
    allocate (results(nalloc))

    link = 0
    nUnique = 0

    ! Loop over all nodes
    do i = 1, N
        if (link(i) == 0) then
            call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)

            ! Expand if necesary and re-run
            if (nfound > nalloc) then
                deallocate (results)
                nalloc = nfound
                allocate (results(nalloc))
                call kdtree2_r_nearest(mytree, pts(:, i), tol2, nFound, nAlloc, results)
            end if

            if (nFound == 1) then
                ! This one is easy, it is already a unique node
                nUnique = nUnique + 1
                link(i) = nUnique
                uniquePts(:, nUnique) = pts(:, i)
            else
                if (link(i) == 0) then
                    ! This node hasn't been assigned yet:
                    nUnique = nUnique + 1
                    uniquePts(:, nUnique) = pts(:, i)

                    do j = 1, nFound
                        link(results(j)%idx) = nUnique
                    end do
                end if
            end if
        end if
    end do

    ! Done with the tree and the result vector
    call kdtree2_destroy(mytree)
    deallocate (results)

end subroutine pointReduce

subroutine getSurfaceCoordinates(coords, n)

    ! Return the surface coordiantes

    use hypData
    implicit none

    integer(kind=intType), intent(in) :: n
    real(kind=realType), intent(inout), dimension(n) :: coords
    integer(kind=intType) :: ierr

    call VecGetArrayF90(X(1), xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Just copy
    coords = xx
    call VecRestoreArrayF90(X(1), xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine getSurfaceCoordinates

subroutine setSurfaceCoordinates(coords, n)

    ! Return the surface coordiantes

    use hypData
    implicit none

    integer(kind=intType), intent(in) :: n
    real(kind=realType), intent(in), dimension(n) :: coords
    integer(kind=intType) :: ierr

    call VecGetArrayF90(X(1), xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Just copy
    xx = coords

    call VecRestoreArrayF90(X(1), xx, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine setSurfaceCoordinates

subroutine assignN2NDirected(nodeConn, nEdge, n1, n2, size1, size2)

    ! This subroutine assigns connection from node 1 to node 2
    ! in the nodeConn matrix. If this connection is already assigned
    ! then we should flag this node as this indicates flipped normals

    use precision
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: n1, n2, size1, size2

    ! Input/Output variables
    integer(kind=intType), dimension(size1, size2), intent(inout) :: nodeConn
    integer(kind=intType), dimension(size2), intent(inout) :: nEdge

    ! Working variables
    integer(kind=intType) :: i
    logical :: foundNode

    foundNode = .False.
    do i = 1, nEdge(n1)
        if (nodeConn(i, n1) == n2) then
            foundNode = .True.
        end if
    end do
    if (foundNode) then
        nodeConn(:, n1) = -1
    else
        nEdge(n1) = nEdge(n1) + 1
        nodeConn(nEdge(n1), n1) = n2
    end if

end subroutine assignN2NDirected

subroutine assignN2N(nodeConn, nEdge, n1, n2, size1, size2)

    ! This subroutine assigns connection from node 1 to node 2
    ! in the nodeConn matrix. If this connection is already assigned
    ! then we should flag this node as this indicates flipped normals

    use precision
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: n1, n2, size1, size2

    ! Input/Output variables
    integer(kind=intType), dimension(size1, size2), intent(inout) :: nodeConn
    integer(kind=intType), dimension(size2), intent(inout) :: nEdge

    ! Working variables
    integer(kind=intType) :: i
    logical :: foundNode

    foundNode = .False.
    do i = 1, nEdge(n1)
        if (nodeConn(i, n1) == n2) then
            foundNode = .True.
        end if
    end do

    if (.not. foundNode) then
        nEdge(n1) = nEdge(n1) + 1
        nodeConn(nEdge(n1), n1) = n2
    end if

end subroutine assignN2N

subroutine findKStretch(XL, XLm1, XLm2)
    !***DESCRIPTION
    !
    !     Written by Ney Secco
    !
    !     Abstract: findKStretch computes a distance ratio. Let zeta and eta be
    !     the computational coordinates along the surface, while L represents
    !     the marching direction (the current layer). This functions compute
    !     the following ratio:
    !     Kstretch =  distance between X(zeta, eta, L) and X(zeta, eta, L-1)
    !                --------------------------------------------------------
    !                distance between X(zeta, eta, L-1) and X(zeta, eta, L-2)
    !     This can be used as an indication of how much the grid is stretching in the
    !     marching direction.
    !
    !     Parameters
    !     ----------
    !     L : integer
    !         The marching level to compute stretch for
    !
    !     Returns
    !     -------
    !     deltaS : real
    !         Marching increment set in hypData
    !
    !     Updates
    !     -------
    !     maxKStretch: real (defined in hypData.F90)
    !         Ratio of distance between layers
    !
    use precision
    use hypInput
    use hypData, only: maxKStretch, nx
#include "petsc/finclude/petsc.h"
    use petsc
    implicit none

    ! Input parameters
    Vec, intent(in) :: XL, XLm1, XLm2

    ! Working
    integer(kind=intType) :: i, ierr
    real(kind=realType) :: rl(3), rlm1(3), rlm2(3), KStretch
    real(kind=realType), pointer, dimension(:) :: xxl, xxlm1, xxlm2
    real(kind=realType), external :: dist

    ! BEGIN EXECUTION

    ! Assign pointers
    call VecGetArrayF90(XL, xxl, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(XLm1, xxlm1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecGetArrayF90(XLm2, xxlm2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Initialize maximum stretch variable
    maxKStretch = 0.0

    ! Compute strech for every point
    do i = 1, nx

        ! Get node coordinates on each layer
        rl = xxl(3 * i - 2:3 * i)
        rlm1 = xxlm1(3 * i - 2:3 * i)
        rlm2 = xxlm2(3 * i - 2:3 * i)

        ! Compute stretch ratio (function dist defined in 3D_utilities.F90)
        KStretch = dist(rl, rlm1) / dist(rlm1, rlm2)

        ! Compare with the maximum value
        maxKStretch = max(KStretch, maxKStretch)

    end do

    ! Restore arrays to make petsc happy
    call VecRestoreArrayF90(XL, xxl, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(XLm1, xxlm1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecRestoreArrayF90(XLm2, xxlm2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

end subroutine findKStretch

subroutine addMissing(nList, n1, n2, n3)

    use precision
    implicit none

    integer(kind=intType), intent(in) :: nList(3), n1, n2
    integer(kind=intType), intent(out) :: n3
    integer(kind=intType) :: i
    ! For a list of 3 integers, nList, and two given integer n1, n2 which
    ! are also in the list, return the other one in the list

    do i = 1, 3
        if (.not. (nList(i) == n1 .or. nList(i) == n2)) then
            n3 = nList(i)
        end if
    end do
end subroutine addMissing

subroutine computeCornerAngle(normal1, normal2, edgeVector, angle)

    ! This subroutine computes the corner angle between two faces. The geometry
    ! is shown below:
    !
    !           B
    !  +--------+--------+
    !  |        |        |
    !  | face1  |  face2 |
    !  |        |        |
    !  |        |<-edge  |
    !  +--------+--------+
    !           A
    !
    ! INPUTS:
    ! normal1 is the normal of face 1 (shouldn't be normalized)
    ! normal2 is the normal of face 2 (shouldn't be normalized)
    ! edgeVector should be a vector defined from A to B.
    !
    ! OUTPUTS:
    ! angle: angle between corners from -pi to pi. Positive for convex corners.
    !
    ! Ney Secco - 2016

    use precision
    implicit none

    ! Input parameters
    real(kind=realType), intent(in) :: normal1(3), normal2(3), edgeVector(3)

    ! Output parameters
    real(kind=realType), intent(out) :: angle

    ! Working variables
    real(kind=realType) :: orthog(3), proj, n1norm(3), n2norm(3), acosArg

    ! EXECUTION

    ! Normalize the normals
    n1norm = normal1 / norm2(normal1)
    n2norm = normal2 / norm2(normal2)

    ! Compute the angle between normals
    acosArg = dot_product(n1norm, n2norm)
    acosArg = min(1.0, max(-1.0, acosArg)) ! This is here to make sure argument is between -1,1
    angle = acos(acosArg)

    ! Find vector orthogonal to normals, from face 1 to face 2
    call cross_prod(n1norm, n2norm, orthog)

    ! Compute projection of the orthogonal vector along the edge vector
    proj = dot_product(orthog, edgeVector)

    ! If the orthogonal vector points in the same direction of the edge vector,
    ! then we have a convex corner. Otherwise, we have a concave corner, and
    ! we need to add a negative sign to the angle
    if (proj < 0) then
        angle = -angle
    end if

end subroutine computeCornerAngle

! subroutine getBCCorner(bcType, p0, p1, p2, p3, p4)
!   ! Apply the boundary condition to compute p3 and p4 in the following
!   ! situation:
!   !
!   !               p2----------+
!   !               |           |
!   !               |           |
!   !               |           |
!   !     p3--------p0----------p1
!   !               |
!   !               |
!   !               |
!   !               p4
!   !
!   ! For this case TWO BCTypes must be given. The first BC determines
!   ! p3 while the second BC determines p4.

!   use hypInput
!   use hypData, only : marchIter
!   implicit none

!   ! Input
!   integer(kind=intType), intent(in), dimension(2) :: bcType
!   real(kind=realType), dimension(3) :: p0, p1, p2

!   ! Output
!   real(kind=realType), dimension(3) :: p3, p4

!   ! Working:
!   real(kind=realType), dimension(3) :: v1, v2, normal, edgeNormal, p0_new
!   real(kind=realType), dimension(3) :: p4_extrap, p4_splay, p3_extrap, p3_splay
!   real(kind=realType) :: blend, dist

!   blend = splayCornerOrthogonality

!   p0_new = p0

!   ! ------------------ Boundary Condition 1 -----------------
!   select case(BCType(1))
!   case (BCSplay)
!      p3_extrap = (2+splay)*p0 - (1+splay)*p1

!      ! Distance between p1 and p0
!      dist = norm2(p1-p0)

!      ! Compute the normal
!      v1 = p1 - p0
!      v1 = v1 / norm2(v1)

!      v2 = p2 - p0
!      v2 = v2 / norm2(v2)

!      call cross_prod(v1, v2, normal)
!      normal = normal / norm2(normal)

!      ! Now cross product for the edge normal
!      call cross_prod(v2, normal, edgeNormal)
!      edgeNormal = edgeNormal / norm2(edgeNormal)

!      p3_splay = p0 - edgeNormal*(1+splay)*dist

!      p3 = (one-blend)*p3_extrap + blend*p3_splay

!   case (BCXSymm)
!      p3 = (/-p1(1), p1(2), p1(3)/)
!      p0_new(1) = zero

!   case (BCYSymm)
!      p3 = (/p1(1), -p1(2), p1(3)/)
!      p0_new(2) = zero

!   case (BCZSymm)
!      p3 = (/p1(1), p1(2), -p1(3)/)
!      p0_new(3) = zero

!   case (BCXConst, BCYConst, BCZConst)
!      p3 = 2*p0 - p1

!      p0_new(2) = zero

!   case (BCAverage)
!      p0_new(2) = zero

!   case default
!      print *,'BC Error'
!      stop
!   end select

!   ! ------------------ Boundary Condition 2 -----------------
!   select case(BCType(2))
!   case (BCSplay)
!      p4_extrap = (2+splay)*p0 - (1+splay)*p2

!      ! Distance between p1 and p0
!      dist = norm2(p2-p0)

!      ! Compute the normal
!      v1 = p1 - p0
!      v1 = v1 / norm2(v1)

!      v2 = p2 - p0
!      v2 = v2 / norm2(v2)

!      call cross_prod(v1, v2, normal)
!      normal = normal / norm2(normal)

!      ! Now cross product for the edge normal
!      call cross_prod(v1, normal, edgeNormal)
!      edgeNormal = edgeNormal / norm2(edgeNormal)

!      p4_splay = p0 + edgeNormal*(1+splay)*dist

!      p4 = (one-blend)*p4_extrap + blend*p4_splay

!   case(BCXSymm)
!      p4 = (/-p2(1), p2(2), p2(3)/)
!      p0_new(1) = zero

!   case (BCYSymm)
!      p4 = (/p2(1), -p2(2), p2(3)/)
!      p0_new(2) = zero

!   case (BCZSymm)
!      p4 = (/p2(1), p2(2), -p2(3)/)
!      p0_new(3) = zero

!   case (BCXConst, BCYConst, BCZConst)
!      !p4 = p4_extrap
!      p4 = 2*p0 - p2
!      p0_new(2) = zero

!   case (BCAverage)
!      p0_new(2) = zero

!   case default
!      print *,'BC Error'
!      stop
!  end select

!  ! Set the updated center coordinate.
!  p0 = p0_new
! end subroutine getBCCorner

! subroutine getBCEdge(bcType, p0, p1, p2, p3, p4)

!   ! Apply the boundary condition to compute p4 in the following
!   ! situation:
!   !
!   !     +---------p2----------+
!   !     |         |           |
!   !     |         |           |
!   !     |         |           |
!   !     p3--------p0----------p1
!   !               |
!   !               |
!   !               |
!   !               p4
!   !
!   use hypInput
!   use hypData, only :marchIter
!   implicit none

!   ! Input
!   integer(kind=intType), intent(in) :: bcType
!   real(kind=realType), dimension(3) :: p0, p1, p2, p3

!   ! Output
!   real(kind=realType), dimension(3) :: p4

!   ! Working:
!   real(kind=realType), dimension(3) :: v1, v2, normal
!   real(kind=realType), dimension(3) :: p4_extrap, p4_splay, edgeNormal
!   real(kind=realType) :: blend, dist

!   blend = splayEdgeOrthogonality

!   select case(bcType)
!   case (BCSplay)

!      ! Regular extrapolation without orthogonality
!      p4_extrap = (2+splay)*p0 - (1+splay)*p2

!      ! Distance between p2 and p0
!      dist = norm2(p2-p0)

!      ! Compute the normal
!      v1 = p1 - p3
!      v1 = v1 / norm2(v1)

!      v2 = p2 - p0
!      v2 = v2/norm2(v2)

!      call cross_prod(v1, v2, normal)
!      normal = normal / norm2(normal)

!      ! Now cross product for the edge normal
!      call cross_prod(v1, normal, edgeNormal)
!      edgeNormal = edgeNormal / norm2(edgeNormal)

!      p4_splay = p0 + edgeNormal*(1+splay)*dist

!      p4 = (one-blend)*p4_extrap + blend*p4_splay

!   case (BCXSymm)
!      !p4 = (/two*BCVal(1) -p2(1), p2(2), p2(3)/)
!      p0(1) = zero
!   case (BCYSymm)
!      !p4 = (/p2(1), two*BCVal(1)-p2(2), p2(3)/)
!      p4 = (/p2(1), zero, p2(3)/)
!      p0(2) = zero
!   case (BCZSymm)
!      !p4 = (/p2(1), p2(2), two*BCVal(1)-p2(3)/)
!      p0(3) = zero
!   case (BCXConst)
!      p4 = two*p0 - p2
!      !p0(1) = BCVal(1)
!   case(BCYConst)
!      p4 = two*p0 - p2
!      !p0(2) = BCVal(1)
!   case(BCZConst)
!      p4 = two*p0 - p2
!      !p0(3) = BCVal(1)
!   case default
!      print *,'BC Error', bcTYpe
!      stop
!  end select

! end subroutine getBCEdge

! subroutine getBCEdgeBlocks(bcType, blk4, blk0, blk2)

!   use hypInput
!   implicit none

!   ! Input
!   integer(kind=intType), intent(in) :: bcType
!   real(kind=realType), dimension(3,3), intent(in) :: blk4

!   ! Output
!   real(kind=realType), dimension(3,3), intent(out) :: blk0, blk2
!   blk2 = zero
!   select case(bcType)
!   case (BCSplay)
!      blk0 = (one + splay)*blk4
!      blk2 = -splay*blk4

!   case (BCXSymm)
!      blk0 = zero
!      blk2 = blk4
!      blk2(:, 1) = -blk2(:, 1)

!   case (BCYSymm)
!      blk0 = zero
!      blk2 = blk4
!      blk2(:, 2) = -blk2(:, 2)
!    case (BCZSymm)
!      blk0 = zero
!      blk2 = blk4
!      blk2(:, 3) = -blk2(:, 3)

!   case (BCXConst, BCYConst, BCZConst)
!      blk0 = zero
!      blk2 = blk4
!   end select

! end subroutine getBCEdgeBlocks

! subroutine getBCCornerBlocks(bcType, blk3, blk4, blk0, blk1, blk2)

!   use hypInput
!   implicit none

!   ! Input
!   integer(kind=intType), intent(in), dimension(2) :: bcType
!   real(kind=realType), dimension(3,3), intent(in) :: blk3, blk4

!   ! Output
!   real(kind=realType), dimension(3,3), intent(out) :: blk0, blk1, blk2

!   ! blk0 needs to be zeroed since we will have two contributions to it
!   blk0 = zero
!   blk1 = zero
!   blk2 = zero
!   ! ------------------ Boundary Condition 1 -----------------
!   select case (bcType(1))
!   case (BCSplay)
!      blk0 = blk0 + (one + splay)*blk3
!      blk1 = -splay*blk3

!   case (BCXSymm)
!      blk1 = blk3
!      blk1(:, 1) = -blk1(:, 1)

!   case (BCYSymm)
!      blk1 = blk3
!      blk1(:, 2) = -blk1(:, 2)

!   case (BCZSymm)
!      blk1 = blk3
!      blk1(:, 3) = -blk1(:, 3)
!   case (BCXConst, BCYConst, BCZConst)
!      blk1 = blk3
!   case (BCAverage)
!      blk1 = blk3
!   end select

!   ! ------------------ Boundary Condition 2 -----------------
!   select case (bcType(2))
!   case (BCSplay)
!      blk0 = blk0 + (one + splay)*blk4
!      blk2 = -splay*blk4
!   case (BCXSymm)
!      blk2 = blk4
!      blk2(:, 1) = -blk2(:, 1)
!   case (BCYSymm)
!      blk2 = blk4
!      blk2(:, 2) = -blk2(:, 2)
!   case (BCZSymm)
!      blk2 = blk4
!      blk2(:, 3) = -blk2(:, 3)
!   case (BCXConst, BCYConst, BCZConst)
!      blk2 = blk4
!   case (BCAverage)
!      blk2 = blk4
!   end select

! end subroutine getBCCornerBlocks

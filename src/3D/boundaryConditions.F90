subroutine getBC(bcType, BCVal, isEdge, blend, p0, p1, p2, p3, p4)

    ! Generic boudnary condition function to apply the boundary
    ! condition to compute p4 in the following edge situation:
    !
    !     +---------p2----------+
    !     |         |           |
    !     |         |           |
    !     |         |           |
    !     p3--------p0----------p1
    !               |
    !               |
    !               |
    !               p4
    !
    !  OR the following corner situation if isEdge is false. Here, p3
    !  is unused
    !
    !               p2----------+
    !               |           |
    !               |           |
    !               |           |
    !     p3--------p0----------p1
    !               |
    !               |
    !               |
    !               p4
    !

    use hypInput
    use hypData, only: marchIter
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: bcType
    real(kind=realType), dimension(3) :: p0, p1, p2, p3, p4
    real(kind=realType), dimension(2) :: bcVal
    logical :: isEdge

    ! Working:
    real(kind=realType), dimension(3) :: v1, v2, normal
    real(kind=realType), dimension(3) :: p4_extrap, p4_splay, edgeNormal
    real(kind=realType) :: blend, dist

    select case (bcType)
    case (BCSplay)

        ! Regular extrapolation without orthogonality
        p4_extrap = (2 + splay(marchIter)) * p0 - (1 + splay(marchIter)) * p2

        ! Distance between p2 and p0
        dist = norm2(p2 - p0)

        ! Compute the normal
        if (isEdge) then
            v1 = p1 - p3 ! Use 2nd order approx for edge
        else
            v1 = p1 - p0 ! Use 1st order approx for corner
        end if
        v1 = v1 / norm2(v1)

        v2 = p2 - p0
        v2 = v2 / norm2(v2)

        call cross_prod(v1, v2, normal)
        normal = normal / norm2(normal)

        ! Now cross product for the edge normal
        call cross_prod(v1, normal, edgeNormal)
        edgeNormal = edgeNormal / norm2(edgeNormal)

        p4_splay = p0 + edgeNormal * (1 + splay(marchIter)) * dist

        p4 = (one - blend) * p4_extrap + blend * p4_splay

    case (BCXSymm)
        p4 = (/two * BCVal(1) - p2(1), p2(2), p2(3)/)
        p0(1) = BCVal(1)

    case (BCYSymm)
        p4 = (/p2(1), two * BCVal(1) - p2(2), p2(3)/)
        p0(2) = BCVal(1)

    case (BCZSymm)
        p4 = (/p2(1), p2(2), two * BCVal(1) - p2(3)/)
        p0(3) = BCVal(1)

    case (BCXConst)
        p4 = two * p0 - p2
        p0(1) = BCVal(1)

    case (BCYConst)
        p4 = two * p0 - p2
        p0(2) = BCVal(1)

    case (BCZConst)
        p4 = two * p0 - p2
        p0(3) = BCVal(1)

    case (BCXYConst)
        p4 = two * p0 - p2
        p0(1) = BCVal(1)
        p0(2) = BCVal(2)

    case (BCYZConst)
        p4 = two * p0 - p2
        p0(2) = BCVal(1)
        p0(3) = BCVal(2)

    case (BCXZConst)
        p4 = two * p0 - p2
        p0(1) = BCVal(1)
        p0(3) = BCVal(2)
    case (bcAverage)
        ! Nothing to do
    case default

        print *, 'BC Error', bcTYpe
        stop
    end select

end subroutine getBC

subroutine getBCBlocks(bcType, blk4, blk0, blk2)

    ! Accumulate the blk4 dependence into blk0 and blk2
    use hypInput
    use hypData, only: marchIter
    implicit none

    ! Input
    integer(kind=intType), intent(in) :: bcType
    real(kind=realType), dimension(3, 3), intent(inout) :: blk4, blk0, blk2

    select case (bcType)
    case (BCSplay)
        blk0 = blk0 + (one + splay(marchIter)) * blk4
        blk2 = blk2 - splay(marchIter) * blk4

    case (BCXSymm)
        blk4(:, 1) = -blk4(:, 1)
        blk2 = blk2 + blk4

    case (BCYSymm)
        blk4(:, 2) = -blk4(:, 2)
        blk2 = blk2 + blk4

    case (BCZSymm)
        blk4(:, 3) = -blk4(:, 3)
        blk2 = blk2 + blk4

    case (BCXConst, BCYConst, BCZConst, &
          BCXYConst, BCXZConst, BCYZConst, BCAverage)
        blk2 = blk2 + blk4
    end select

end subroutine getBCBlocks

subroutine setBCVal(bcType, X, bcVal)

    ! This subroutine computes the 'value' for a BC depending on the BC
    ! Type and the spatial coordinate passed in.
    use precision
    use hypInput
    implicit none

    ! Input/output
    integer(kind=intType), intent(in) :: bcType
    real(kind=realType), dimension(3), intent(in) :: X
    real(kind=realType), intent(out) :: bcVal(2)

    ! Note that not all BCs need a value and are not included here.
    select case (bcType)
    case (BCXSymm, BCXConst)
        bcVal(1) = X(1)
    case (BCYSymm, BCYConst)
        bcVal(1) = X(2)
    case (BCZSymm, BCZConst)
        bcVal(1) = X(3)
    case (BCXYConst)
        bcVal(1) = X(1)
        bcVal(2) = X(2)
    case (BCXZConst)
        bcVal(1) = X(1)
        bcVal(2) = X(3)
    case (BCYZConst)
        bcVal(1) = X(2)
        bcVal(2) = X(3)
    end select

end subroutine setBCVal

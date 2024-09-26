module hypInput
    use precision

    implicit none

    ! This module contains the input information for the hyperbolic grid
    ! generator

    ! Input parameters. See pyHyp.py for what each parameter means

    real(kind=realType) :: ps0, pgridratio
    real(kind=realType) :: slexp
    real(kind=realType) :: kspRelTol, cmax, marchDist
    real(kind=realType) :: nodeTol, symTol
    integer(kind=intType) :: N, nConstantStart, nConstantEnd, nTruncate
    integer(kind=intType) :: kspMaxIts, preConLag, kspSubspaceSize
    integer(kind=intType) :: coarsen

    real(kind=realType), dimension(:), allocatable, target :: fullDeltaS

    integer(kind=intType), dimension(:), allocatable, target :: volSmoothIter
    real(kind=realType), dimension(:), allocatable, target :: volBlend
    real(kind=realType), dimension(:), allocatable, target :: epsE
    real(kind=realType), dimension(:), allocatable, target :: epsI
    real(kind=realType), dimension(:), allocatable, target :: theta
    real(kind=realType), dimension(:), allocatable, target :: splay
    real(kind=realType), dimension(:), allocatable, target :: splayEdgeOrthogonality
    real(kind=realType), dimension(:), allocatable, target :: splayCornerOrthogonality
    real(kind=realType), dimension(:), allocatable, target :: cornerAngle
    real(kind=realType), dimension(:), allocatable, target :: volCoef

    logical :: nonLinear, writeMetrics, unattachedEdgesAreSymmetry
    logical :: noPointReduce


    ! Input boundary condition information
    integer(kind=intType), dimension(:, :), allocatable :: BCs

    ! Family names for wall surfaces
    character(maxCGNSNameLen), dimension(:), allocatable :: families

    ! Elliptic Parameters
    real(kind=realType) :: farFieldTol
    integer(kind=intType) :: evalMode
    integer(kind=intType), parameter :: EVAL_EXACT = 0
    integer(kind=intType), parameter :: EVAL_SLOW = 1
    integer(kind=intType), parameter :: EVAL_FAST = 2
    logical :: useMatrixFree
    character(len=512) :: sourceStrengthFile

    ! Topology types:
    integer(kind=intType), parameter :: topoUnknown = -1
    integer(kind=intType), parameter :: topoInternal = 0
    integer(kind=intType), parameter :: topoCorner = 1
    integer(kind=intType), parameter :: topoEdge = 2
    integer(kind=intType), parameter :: topoLCorner = 3

    ! Boundary condition types:
    integer(kind=intType), parameter :: BCDefault = -1
    integer(kind=intType), parameter :: BCSplay = 0
    integer(kind=intType), parameter :: BCXSymm = 1
    integer(kind=intType), parameter :: BCYSymm = 2
    integer(kind=intType), parameter :: BCZSymm = 3
    integer(kind=intType), parameter :: BCXConst = 4
    integer(kind=intType), parameter :: BCYConst = 5
    integer(kind=intType), parameter :: BCZConst = 6
    integer(kind=intType), parameter :: BCXYConst = 7
    integer(kind=intType), parameter :: BCYZConst = 8
    integer(kind=intType), parameter :: BCXZConst = 9

    integer(kind=intType), parameter :: BCAverage = 10
    ! Boundary condition side:
    integer(kind=intType), parameter :: iLow = 1
    integer(kind=intType), parameter :: iHigh = 2
    integer(kind=intType), parameter :: jLow = 3
    integer(kind=intType), parameter :: jHigh = 4

    ! FileType
    integer(kind=intType), parameter :: cgnsFileType = 1
    integer(kind=intType), parameter :: plot3dFileType = 2
    integer(kind=intType), parameter :: patchInput = 3

    ! Farfield type selection
    integer(kind=intType) :: outerFaceType
    integer(kind=intType), parameter :: outerFaceFarfield = 1
    integer(kind=intType), parameter :: outerFaceOverset = 2

end module hypInput

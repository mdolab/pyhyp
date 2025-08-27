module panel

    ! This module containes the "panel" derived type used for storing
    ! panel information for the ellipic method
    use precision
    implicit none
    save
    integer(kind=intType), parameter :: NNodesMax = 4
    type panelType
        ! Data for actual panel
        integer(kind=intType) :: N
        real(kind=realType), dimension(3, NNodesMax + 1) :: X
        real(kind=realType), dimension(3, NNodesMax + 1) :: pts
        real(kind=realType), dimension(3) :: normal
        real(kind=realType), dimension(3, 3) :: C, CT

        ! Common data
        real(kind=realType) :: area
        real(kind=realType) :: strength
        real(kind=realType), dimension(3) :: center
        real(kind=realType) :: length

        ! Data for grouped panel
        integer(kind=intType), dimension(9) :: children
        integer(kind=intType) :: nChildren

    end type panelType

    type PanelLevel
        ! Number of panels on level
        integer(kind=intType) :: N
        type(panelType), dimension(:), allocatable :: panels
    end type PanelLevel

    ! Array of Panel Levels (MGP = multi-grid-panles)
    type(panelLevel), dimension(:), target, allocatable :: MGP
    integer(kind=intType) :: levelMax
    integer(kind=intType), dimension(:), allocatable :: coarseIndices
    integer(kind=intType) :: nCoarse
    integer(kind=intType) :: nFar, nNear

    ! Generic pointer for a panel
    type(panelType), pointer :: pp
end module panel

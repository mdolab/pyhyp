module panel

  ! This module containes the "panel" derived type used for storing
  ! panel information for the ellipic method
  use precision
  implicit none
  save
  integer(kind=intType), parameter :: NNodesMax = 8
  type panelType
     ! Data for actual panel
     integer(kind=intType) :: N
     real(kind=realType), dimension(3, NNodesMax) :: X
     real(kind=realType), dimension(3) :: normal
     real(kind=realType), dimension(3, 3) :: C
     real(kind=realType), dimension(3, 3) :: CT

     ! Common data
     real(kind=realType) :: area
     real(kind=realType) :: strength
     real(kind=realType), dimension(3) :: center

     ! Data for grouped panel
     integer(kind=intType), dimension(:), allocatable :: children
     integer(kind=intType) :: nChildren

  end type panelType

  type PanelLevel
     ! Number of panels on level
     integer(kind=intType) :: N
     type(panelType), dimension(:), allocatable :: panels
  end type PanelLevel

  ! Array of Panel Levels (MGP = multi-grid-panles)
  type(panelLevel), dimension(:), target, allocatable :: MGP

end module panel


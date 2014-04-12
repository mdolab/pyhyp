module panel

  ! This module containes the "panel" derived type used for storing
  ! panel information for the ellipic method
  use precision
  implicit none
  save
  integer(kind=intType), parameter :: NNodesMax = 8
  type panelType
     integer(kind=intType) :: N
     real(kind=realType), dimension(3, NNodesMax) :: X
     real(kind=realType) :: area
     real(kind=realType) :: strength
     real(kind=realType), dimension(3) :: center
     real(kind=realType), dimension(3) :: normal
     real(kind=realType), dimension(3, 3) :: C
     real(kind=realType), dimension(3, 3) :: CT
  end type panelType

  ! List of the zero-level panels
  type(panelType), dimension(:), target, allocatable :: panels

end module panel

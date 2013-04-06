subroutine writeCGNS_2d(fileName)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeCGNS_2d write the current grid to a 3D CGNS
  !               file that is exactly 1 cell wide. It also applies
  !               default boundary and symmetry conditions such that
  !               the grid can be directly used in a 3D flow solver
  !               such as sumb.
  !
  !     Description of Arguments
  !     Input:
  !     fileNmae - Character array: the name of the cgns file
  !
  !     Ouput: None

  use precision
  use hypInput

  ! Input Arguments
  character*(*) :: fileName

  print *,'This function is not written yet'
  
end subroutine writeCGNS_2d

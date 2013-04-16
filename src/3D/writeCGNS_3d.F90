subroutine writeCGNS_3D(fileName)

  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract: writeCGNS_3d write the current grid to a 3D CGNS
  !               file. It computes grid connectivities such that the
  !               grid and boundary condition information such that
  !               the grid can be used directly in a 3D flow solver
  !               such  as SUmb.
  !
  !     Description of Arguments
  !     Input:
  !     fileNmae - Character array: the name of the cgns file
  !
  !     Ouput: None

  use hypInput
  use hypData

  implicit none
  include 'cgnslib_f.h'
  
  ! Input Arguments
  character*(*) :: fileName
  
  ! Working Variables
  integer(kind=intType) :: cg, ierr
  integer(kind=intType) :: cellDim, physDim, base, coordID, gridShp(3)
  character*256 :: zoneName
  integer(kind=intType) :: sizes(9), zone, ii, i, j, k, nx, transform(3)
  integer(kind=intType) :: pnts(3,2), pnts_donor(3,2), BCOut, nCon
  real(kind=realType), dimension(:), allocatable :: coordArray

  ! Open the CGNS File:
  print *,'Not coded yet!'
  
end subroutine writeCGNS_3D

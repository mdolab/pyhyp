module hypData
  use precision
  implicit None

#include "include/finclude/petsc.h"
  
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! Data for the 2D generator:
  real(kind=realType), dimension(:, :, :), allocatable, target :: grid2D

  ! Data for both generators
  real(kind=realType) :: scaleDist

  ! Petsc Varibles
  Mat A
  Vec rDelta
  Vec RHS
  KSP kspObj

end module hypData

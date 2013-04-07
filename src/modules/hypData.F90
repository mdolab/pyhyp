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

  ! Petsc Varibles for solving linearized hyperbolic system 
  Mat hypMat
  Vec hypDelta
  Vec hypRHS
  KSP hypKSP

 ! ! Petsc Varibles for solving the smoothing problem
 !  Mat smoothMat
 !  Vec smoothDelta
 !  Vec smoothRHS
 !  KSP smoothKSP
 !  logical :: smoothingAssembled

end module hypData

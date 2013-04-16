module hypData
  use precision
  implicit None

#include "include/finclude/petsc.h"
  
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! Data for the 2D generator:
  real(kind=realType), dimension(:, :, :), allocatable, target :: grid2D
  real(kind=realType), dimension(:, :, :), allocatable, target :: grid3D

  ! Data for both generators
  real(kind=realType) :: scaleDist

  ! Petsc Varibles for solving linearized hyperbolic system 
  Mat hypMat
  Vec hypDelta
  Vec hypRHS
  KSP hypKSP

  type patchType

     ! Patch Dimensions
     ! il, jl: Use sumb/pyWarp labeling

     integer(kind=intType) :: il, jl

     ! l_index: Pointer of each node on the patch into the globally
     ! reduced X array. 
     integer(kind=intType), dimension(:, :), allocatable :: l_index
     
  end type patchType

  ! The list of the patchs
  integer(kind=intType) :: nPatch
  type(patchType), dimension(:), allocatable :: patches
  integer(kind=intType), dimension(:, :), allocatable :: nPtr
end module hypData

module hypData
  use precision
  implicit None

#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! Data for the 2D generator:
  real(kind=realType), dimension(:, :, :), allocatable, target :: grid2D

  ! Data for both generators
  real(kind=realType) :: scaleDist
  logical :: factorNext
  integer(kind=intType) :: nx
  real(kind=realType), pointer :: xxm1(:), xx(:), xxp1(:)

  ! Data used for convergence info:
  real(kind=realType) :: timeStart, gridSensorMax, gridSensorMin, minQuality, deltaS, minR
  integer(kind=intType) :: marchIter, kspIts, l_0
  real(kind=realType) :: radius, radius0, Xavg(3), cratio, sl

  integer(kind=intType) :: Nlayers,  smoothIter

  ! Petsc Varibles for solving linearized hyperbolic system 
  Mat hypMat, hypMatFD, hypMatPC
  Vec hypDelta
  Vec hypRHS
  KSP hypKSP
  PC hypPC
  SNES hypSNES
  Vec hypRes
  Vec Volume
  PetscFortranAddr   ctx(1)
  Vec, dimension(:), allocatable :: X
  Vec normalVec
  Vec ovrNNeighbours

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

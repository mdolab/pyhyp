module hypData
  use precision
  implicit None

#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! Data for the 2D generator:
  real(kind=realType), dimension(:, :, :), allocatable, target :: grid2D
  real(kind=realType), dimension(:, :), allocatable :: X0, X1, Xm1

  ! Data for both generators
  real(kind=realType) :: scaleDist
  real(kind=realType) :: gridRatio
  logical :: factorNext
  integer(kind=intType) :: nx
  real(kind=realType), pointer, dimension(:) :: xxm2tmp, xxm1tmp, rrtmp, xxtmp, deltaTmp
  real(kind=realType), dimension(:,:), allocatable :: xxm1, xxm2, xx, rr, pxxm1, xxInterp
  real(kind=realType), dimension(:), allocatable :: volume
  integer(kind=intType), dimension(:), allocatable :: inds
  integer(kind=intType) :: nSubIterPrev
  ! Data used for convergence info:
  real(kind=realType) :: timeStart, gridSensorMax, gridSensorMin, minQuality, deltaS, minR
  integer(kind=intType) :: marchIter, kspIts
  real(kind=realType) :: radius, radius0, Xavg(3), cratio, sl, vBar
  integer(kind=intType) :: Nlayers,  smoothIter
  integer(kind=intType) :: nSubIter, subIter
  real(kind=realType) :: desiredS
  logical :: three_d_vars_allocated = .False.
  logical :: two_d_vars_allocated = .False.

  ! Petsc Varibles for solving linearized hyperbolic system 
  Mat hypMat, hypMatFD, hypMatPC
  Vec hypDelta, hypRHS, hypRes
  KSP hypKSP
  SNES hypSNES
  PetscFortranAddr   ctx(1)
  Vec, dimension(:), allocatable :: X, X_ksi, X_eta, X_zeta, X_ksi_ksi, X_eta_eta, X_diss, Vhist
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

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

  ! -------------------------------------
  !         Surface Patch Data
  ! -------------------------------------
  type patchType
     ! Patch Dimensions
     integer(kind=intType) :: il, jl

     ! l_index: Pointer of each node on the patch into the X array
     integer(kind=intType), dimension(:, :), allocatable :: l_index
  end type patchType
  integer(kind=intType) :: nPatch
  type(patchType), dimension(:), allocatable :: patches

  ! ----------------------------
  !         Grid Data
  ! ----------------------------

  ! The gloabl (total) number of nodes
  integer(kind=intType) :: nXGlobal

  ! The local number of nodes
  integer(kind=intType) :: nX 
  
  ! List of petsc vectors. One for each k-plane
  Vec, dimension(:), allocatable :: X 

  ! PETSc Vectors with ghosting for pseudo marching
  Vec XL, XLm1, XL_local, XLm1_local

  ! Local X vectors for volume/quality calc
  Vec Xlocal, Xlocalm1

  ! Fortran pointers to get data from the above vectors
  real(kind=realType), pointer, dimension(:) :: xxm1, xx

  ! Node pointer array storing neighbours for each node. One is for
  ! local indexing in the ghosted array, the other is the global
  ! ordering for assembly. 
  integer(kind=intType), dimension(:, :), allocatable :: gnPtr
  integer(kind=intType), dimension(:, :), allocatable :: lnPtr
 
  ! The (local) initial surface nodes 
  real(kind=realType), dimension(:, :), allocatable :: xSurf

  ! The indices of the ghost nodes this processor requires
  integer(kind=intType) :: nGhost
  integer(kind=intType), dimension(:), allocatable :: ghost

  ! The "nodal volume" of each local node
  Vec Volume, VolumeLocal
  real(kind=realType), dimension(:), pointer :: Vptr

  ! ---------------------------------------
  !    Marching Iteration/Monitoring Data
  ! ---------------------------------------
  real(kind=realType) :: timeStart
  real(kind=realType) :: scaleDist
  real(kind=realType) :: gridRatio
  real(kind=realType) :: gridSensorMax, gridSensorMin, minQuality, deltaS, minR
  integer(kind=intType) :: marchIter, kspIts
  real(kind=realType) :: radius, radius0, Xavg(3), cratio, sl, vBar
  integer(kind=intType) :: nSubIter, nSubIterPrev
  real(kind=realType) :: desiredS

  ! ---------------------------------------
  !    Miscellaneous Variables
  ! ---------------------------------------
  logical :: three_d_vars_allocated = .False.
  logical :: two_d_vars_allocated = .False.

  ! -------------------------------------
  !         PETSc Linear System Variables
  ! -------------------------------------
  Mat hypMat
  Vec hypDelta, hypRHS, hypRes
  KSP hypKSP
  VecScatter rootScatter
  ! ------------------------------------------------------
  !         Variables for Storing Metrics (for debugging)
  ! ------------------------------------------------------
  Vec, dimension(:, :), allocatable :: metrics
  integer(kind=intType), parameter :: iX_ksi = 1
  integer(kind=intType), parameter :: iX_eta = 2
  integer(kind=intType), parameter :: iX_zeta = 3
  integer(kind=intType), parameter :: iX_ksi_ksi = 4
  integer(kind=intType), parameter :: iX_eta_eta = 5
  integer(kind=intType), parameter :: iX_diss = 6
  integer(kind=intType), parameter :: iVhist = 7
  integer(kind=intType), parameter :: nMetric = 7
  logical :: metricsAllocated = .False.
end module hypData

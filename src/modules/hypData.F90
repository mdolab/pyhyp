module hypData
  use precision
  implicit None

#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! -------------------------------------
  !         Surface Patch Data
  ! -------------------------------------
  type patchType
     ! Patch Dimensions
     integer(kind=intType) :: il, jl

     ! l_index: Pointer of each node on the patch into the X array
     integer(kind=intType), dimension(:, :), allocatable :: l_index

     ! X: Coordinates of the patch
     real(kind=realType), dimension(:, :, :), allocatable :: X

     ! symNodes : list of indices that symmetry nodes and must be
     ! zeroed as such when writing the mesh
     integer(kind=intType) :: nSym
     integer(kind=intType), dimension(:, :), allocatable :: symNodes

     ! Array determining the freezing weights of the nodes. 
     real(kind=realType), dimension(:, :), allocatable :: weights
     
  end type patchType

  integer(kind=intType) :: nPatch
  type(patchType), dimension(:), allocatable :: patches

  ! ----------------------------
  !         Grid Data
  ! ----------------------------

  ! The gloabl (total) number of nodes
  integer(kind=intType) :: nXGlobal

  ! The gloabl (total) number of panels
  integer(kind=intType) :: faceTotal, nPGlobal

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
  integer(kind=intType), dimension(:, :), allocatable :: cPtr, fullcPtr
  integer(kind=intType), dimension(:, :), allocatable :: conn, fullConn
  integer(kind=intTYpe) :: nLocalFace

  ! The indices of the ghost nodes this processor requires
  integer(kind=intType) :: nGhost
  integer(kind=intType), dimension(:), allocatable :: ghost

  ! The "nodal volume" of each local node
  Vec Volume, VolumeLocal
  real(kind=realType), dimension(:), pointer :: Vptr, dptr

  ! ---------------------------------------
  !    Marching Iteration/Monitoring Data
  ! ---------------------------------------
  double precision :: timeStart
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

  ! -------------------------------------
  !         PETSc Linear System Variables
  ! -------------------------------------
  Mat hypMat
  Vec hypDelta, hypRHS, hypRes
  KSP hypKSP
  VecScatter rootScatter
  VecScatter allScatter
  Vec allGlobalNodes

  Mat ellipMat, ellipPCMat
  Vec ellipRHS, ellipSol, localSol
  KSP ellipKSP
  PC pc
  VecScatter ellipScatter, ellipDeltaScatter
  PetscFortranAddr ctx(1)
  Vec ellipDelta,  ellipDeltaLocal
  Vec ellipNorm,   ellipNormLocal
  Vec ellipLayer,  ellipLayerLocal

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

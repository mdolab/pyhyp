module hypData
  use precision
  implicit None
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscksp.h"
#else
#include "include/finclude/petsc.h"
#include "petsc/finclude/petscksp.h"
#endif
  ! This module contains the data and data structures required for
  ! running the hyperbolic grid generator

  ! ----------------------------
  !         Grid Data
  ! ----------------------------

  ! The global (total) number of nodes
  integer(kind=intType) :: nXGlobal

  ! The global (total) number of panels
  integer(kind=intType) :: faceTotal, nPGlobal

  ! The local number of nodes
  integer(kind=intType) :: nx

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
  integer(kind=intType) :: nLocalFace

  ! The indices of the ghost nodes this processor requires
  integer(kind=intType) :: nGhost
  integer(kind=intType), dimension(:), allocatable :: ghost

  ! The "nodal volume" of each local node
  Vec Volume, VolumeLocal
  real(kind=realType), dimension(:), pointer :: Vptr

  ! -------------------------------------
  !         Boundary Condition Data
  ! -------------------------------------
  integer(kind=intType), dimension(:), allocatable :: fullTopoType, topoType
  integer(kind=intType), dimension(:, :), allocatable :: fullBCType, BCType
  real(kind=realType), dimension(:, :, :), allocatable :: fullBCVal, BCVal

  ! -------------------------------------
  !         Surface Patch Data
  ! -------------------------------------
  type patchType
     ! Patch Dimensions
     integer(kind=intType) :: il, jl

     ! l_index: Pointer of each node on the patch into the X array
     integer(kind=intType), dimension(:, :), pointer :: l_index

     ! X: Coordinates of the patch
     real(kind=realType), dimension(:, :, :), pointer :: X

     ! Array determining the freezing weights of the nodes.
     real(kind=realType), dimension(:, :), allocatable :: weights

  end type patchType

  ! This is only needed by the Elliptical generator
  integer(kind=intType) :: nPatch
  type(patchType), dimension(:), allocatable :: patches
  type(patchType), dimension(:), allocatable :: patchIO

  ! ---------------------------------------
  !    Marching Iteration/Monitoring Data
  ! ---------------------------------------
  double precision :: timeStart
  real(kind=realType) :: scaleDist
  real(kind=realType) :: gridRatio
  real(kind=realType) :: gridSensorMax, gridSensorMin, minQuality, deltaS, minR, minVolume
  real(kind=realType) :: minQualityOverall = one, minVolumeOverall = one
  integer(kind=intType) :: marchIter, kspIts, nAverage
  real(kind=realType) ::  Xavg(3), cratio, sl, vBar
  integer(kind=intType) :: nSubIter, nSubIterPrev
  real(kind=realType) :: desiredDeltaS, desiredS, maxKStretch

  ! ---------------------------------------
  !    Miscellaneous Variables
  ! ---------------------------------------
  logical :: varsAllocated = .False.

  ! -------------------------------------
  !         PETSc Linear System Variables
  ! -------------------------------------
  Mat hypMat
  Vec hypDelta, hypRHS, hypRes
  KSP hypKSP
  VecScatter rootScatter
  VecScatter allScatter
  Vec allGlobalNodes

  ! Temp Petsc variables for setting up KSP
  PC globalPC, subPC
  KSP subKSP

  ! Elliptic variables
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

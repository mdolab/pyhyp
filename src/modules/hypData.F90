module hypData
  use precision
  implicit None
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
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
  integer(kind=intType) :: nLocalFace

  ! The indices of the ghost nodes this processor requires
  integer(kind=intType) :: nGhost
  integer(kind=intType), dimension(:), allocatable :: ghost

  ! The "nodal volume" of each local node
  Vec Volume, VolumeLocal
  real(kind=realType), dimension(:), pointer :: Vptr, dptr

  ! -------------------------------------
  !         Boundary Condition Data
  ! -------------------------------------

  ! maximum number of neighbors necessary for boundary conditions
  integer(kind=intType), parameter :: maxBCneighbors = 3

  ! Array that stores the boundary condition type and the neighbor nodes 
  ! the will be used in these boundary conditions
  ! The following structure is used:
  ! BCneighbors(nodeID,:) = (\ bcType, neighbor1, neighbor2, neighbor3 \)
  ! Its dimensions should be (numLocalNodes,maxBCneighbors+1)
  ! This first version stores the local indices of the neighbors
  integer(kind=intType), dimension(:,:), allocatable :: BCneighborsLocal
  ! This other vesion stores the global indices of the neighbors
  integer(kind=intType), dimension(:,:), allocatable :: BCneighborsGlobal

  ! Define the indices of each BC
  ! Higher indices will have higher priority
  integer(kind=intType), parameter :: SplayBCindex = 1
  integer(kind=intType), parameter :: SymmetryBCindex = 2
  integer(kind=intType), parameter :: ConstXBCindex = 3
  integer(kind=intType), parameter :: ConstYBCindex = 4
  integer(kind=intType), parameter :: ConstZBCindex = 5

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

     ! BCnodes
     ! This array will be of size (il,jl,2*maxBCneighbors+1)
     ! This will give i,j coordinates of the neighbor nodes that should be used
     ! when applying a boundary condition.
     ! The following structure is used:
     ! BCnodes(node_i, node_j, :) = (\ bcType,
     !                                 neighbor1_i, neighbor1_j,
     !                                 neighbor2_i, neighbor2_j,
     !                                 neighbor3_i, neighbor3_j /)
     integer(kind=intType), dimension(:,:,:), allocatable :: BCnodes

     ! Array determining the freezing weights of the nodes. 
     real(kind=realType), dimension(:, :), allocatable :: weights
     
  end type patchType

  ! This is only needed by the Elliptical generator
  integer(kind=intType) :: nPatch
  type(patchType), dimension(:), allocatable :: patches

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
  real(kind=realType) :: desiredS, maxKStretch

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

  ! Do we need this?
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

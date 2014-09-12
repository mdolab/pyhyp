module hypInput
  use precision
  
  implicit none

  ! This module contains the input information for the hyperbolic grid
  ! generator

  ! Input parameters. See pyHyp.py for what each parameter means

  real(kind=realType) :: s0, ps0, pgridratio
  real(kind=realType) :: slexp
  real(kind=realType) :: epsE, volCoef, volBlend, epsI, theta
  real(kind=realType) :: kspRelTol,rMin, cmax
  real(kind=realType) :: nodeTol, symTol

  integer(kind=intType) :: N
  integer(kind=intType) :: volSmoothIter, kspMaxIts, preConLag, kspSubspaceSize

  logical ::  nonLinear, writeMetrics

  ! Elliptic Parameters
  real(kind=realType) :: farFieldTol
  integer(kind=intType) :: evalMode
  integer(kind=intType), parameter :: EVAL_EXACT = 0
  integer(kind=intType), parameter :: EVAL_SLOW = 1
  integer(kind=intType), parameter :: EVAL_FAST = 2
  logical :: useMatrixFree
  character(len=512) :: sourceStrengthFile

  ! Grid mirroring parameters
  integer(kind=intType) :: mirrorType
  integer(kind=intType), parameter :: nomirror = 0
  integer(kind=intType), parameter :: xmirror = 1
  integer(kind=intType), parameter :: ymirror = 2
  integer(kind=intType), parameter :: zmirror = 3

end module hypInput

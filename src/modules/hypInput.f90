module hypInput
  use precision
  
  implicit none

  ! This module contains the input information for the hyperbolic grid
  ! generator

  ! Input parameters. See pyHyp.py for what each parameter means
  integer(kind=intType) :: N
  real(kind=realType) :: s0
  real(kind=realType) :: gridRatio
  real(kind=realType) :: slexp
  real(kind=realType) :: epsE, volCoef, volBlend
  integer(kind=intType) :: volSmoothIter, kspMaxIts, preConLag, kspSubspaceSize
  real(kind=realType) :: kspRelTol,rMin
  logical :: writeMirror
end module hypInput

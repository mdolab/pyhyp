module hypInput
  use precision
  
  implicit none

  ! This module contains the input information for the hyperbolic grid
  ! generator

  ! Input parameters. See pyHyp.py for what each parameter means
  integer(kind=intType) :: N, nMax
  real(kind=realType) :: s0, ps0
  real(kind=realType) :: gridRatio, pGridRatio
  real(kind=realType) :: epsE, epsI, theta, volCoef, volBlend
  integer(kind=intType) :: volSmoothIter, kspMaxIts, preConLag, kspSubspaceSize
  real(kind=realType) :: kspRelTol,rMin
  
end module hypInput

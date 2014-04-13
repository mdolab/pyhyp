module hypInput
  use precision
  
  implicit none

  ! This module contains the input information for the hyperbolic grid
  ! generator

  ! Input parameters. See pyHyp.py for what each parameter means
  integer(kind=intType) :: N
  real(kind=realType) :: s0, ps0, pgridratio
  real(kind=realType) :: slexp
  real(kind=realType) :: epsE, volCoef, volBlend, epsI, theta
  integer(kind=intType) :: volSmoothIter, kspMaxIts, preConLag, kspSubspaceSize
  real(kind=realType) :: kspRelTol,rMin, cmax
  logical :: writeMirror, nonLinear, writeMetrics
  real(kind=realType) :: eps
end module hypInput

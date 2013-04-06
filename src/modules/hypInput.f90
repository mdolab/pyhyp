module hypInput
  use precision
  
  implicit none

  ! This module contains the input information for the hyperbolic grid
  ! generator

  ! Parameters that affect both 2D and 3D generation:
  ! N: The number of points in the extrusion direction 
  integer(kind=intType) :: N

  ! s0: Initial offwall spacing
  real(kind=realType) :: s0

  ! Grid Spacing Ratio
  real(kind=realType) :: gridRatio

  ! epsE: The explict smoothing coefficient
  real(kind=realType) :: epsE

  ! epsI: The implicit smoothing coefficient
  real(kind=realType) :: epsI

  ! theta: The barth implicit smoothing coefficient
  real(kind=realType) :: theta

  ! volCoef: The volume smoothing coefficinet
  real(kind=realType) :: volCoef
  
  ! volSmoothIter: The number of point-jacobi volume smoothing
  ! iterations
  integer(kind=intType) :: volSmoothIter

end module hypInput

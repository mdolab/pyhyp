module precision
  save 
  integer,parameter :: realType = selected_real_kind(12)
  integer(kind=4), private :: dummyInt
  integer,parameter :: intType = kind(dummyInt)


  ! Floating point constants.
  real(kind=realType), parameter :: zero  = 0.0_realType
  real(kind=realType), parameter :: one   = 1.0_realType
  real(kind=realType), parameter :: two   = 2.0_realType
  real(kind=realType), parameter :: three = 3.0_realType
  real(kind=realType), parameter :: four  = 4.0_realType
  real(kind=realType), parameter :: five  = 5.0_realType
  real(kind=realType), parameter :: six   = 6.0_realType
  real(kind=realType), parameter :: eight = 8.0_realType
  
  real(kind=realType), parameter :: half   = 0.5_realType
  real(kind=realType), parameter :: third  = one/three
  real(kind=realType), parameter :: fourth = 0.25_realType
  real(kind=realType), parameter :: sixth  = one/six
  real(kind=realType), parameter :: eighth = 0.125_realType
  real(kind=realType), parameter :: pi    = 3.1415926535897931_realType

end module precision

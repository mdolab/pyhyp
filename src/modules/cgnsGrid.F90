module cgnsGrid

  use precision
#ifdef USECGNSMODULE
  use cgns
  implicit none
#else
  implicit none
  include "cgnslib_f.h"
  integer(kind=4), private :: dummyInt
  integer, parameter :: cgsize_t=kind(dummyInt)
#endif

end module cgnsGrid

module communication

  use precision
  implicit none
  save

  ! hyp_comm_world: The communicator of this processor group.
  ! hyp_comm_self : The single processor communicator 
  ! myID:            My processor number in hyp_comm_world.
  ! nProc:           The number of processors in hyp_comm_world.

  integer :: hyp_comm_world, hyp_comm_self, myID, nProc

end module communication

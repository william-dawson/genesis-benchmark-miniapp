! Keep the Fortran MPI communicator and MPI_IN_PLACE conversions in MPI's binding.
subroutine genesis_pme_mesh_sum(mesh, count, communicator)
  use mpi
  use constants_mod, only: dp
  implicit none
  integer, intent(in) :: count, communicator
  real(dp), intent(inout) :: mesh(count)
  integer :: ierr, abort_error

  call MPI_Allreduce(MPI_IN_PLACE, mesh, count, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, communicator, ierr)
  if (ierr /= MPI_SUCCESS) then
    call MPI_Abort(communicator, ierr, abort_error)
    error stop 'GPU PME mesh reduction failed'
  end if
end subroutine

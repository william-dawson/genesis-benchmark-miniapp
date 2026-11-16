! Keep the Fortran MPI communicator and MPI_IN_PLACE conversions in MPI's binding.
subroutine genesis_pme_mesh_sum(mesh, count, communicator)
  use mpi
  use constants_mod, only: wp, mpi_wp_real
  implicit none
  integer, intent(in) :: count, communicator
  real(wp), intent(inout) :: mesh(count)
  integer :: ierr, abort_error

  call MPI_Allreduce(MPI_IN_PLACE, mesh, count, mpi_wp_real, &
                     MPI_SUM, communicator, ierr)
  if (ierr /= MPI_SUCCESS) then
    call MPI_Abort(communicator, ierr, abort_error)
    error stop 'GPU PME mesh reduction failed'
  end if
end subroutine

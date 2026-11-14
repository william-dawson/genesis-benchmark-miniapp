! Use GENESIS's OpenMP runtime for the host layout conversions.
subroutine genesis_pme_host_pack(ncell, stride, atoms, natom, offsets, coord, charge, packed, workers)
  implicit none
  integer, intent(in) :: ncell, stride, atoms, workers
  integer, intent(in) :: natom(ncell), offsets(0:ncell)
  real(8), intent(in) :: coord(3,stride,ncell), charge(stride,ncell)
  real(8), intent(out) :: packed(4*atoms)
  integer :: cell, atom, index

  !$omp parallel do num_threads(workers) schedule(static) private(atom,index)
  do cell = 1, ncell
    do atom = 1, natom(cell)
      index = offsets(cell-1) + atom
      packed(3*index-2:3*index) = coord(:,atom,cell)
      packed(3*atoms+index) = charge(atom,cell)
    end do
  end do
  !$omp end parallel do
end subroutine

subroutine genesis_pme_host_clear(ncell, stride, first, force) bind(c)
  use iso_c_binding, only: c_int, c_double
  implicit none
  integer(c_int), intent(in) :: ncell, stride, first
  real(c_double), intent(inout) :: force(3,stride,ncell)
  integer :: cell

  !$omp parallel do schedule(static)
  do cell = first, ncell
    force(:,:,cell) = 0.0_c_double
  end do
  !$omp end parallel do
end subroutine

subroutine genesis_pme_host_scatter(ncell, stride, atoms, natom, offsets, packed, force, workers)
  implicit none
  integer, intent(in) :: ncell, stride, atoms, workers
  integer, intent(in) :: natom(ncell), offsets(0:ncell)
  real(8), intent(in) :: packed(3,atoms)
  real(8), intent(inout) :: force(3,stride,ncell)
  integer :: cell, atom, index

  !$omp parallel do num_threads(workers) schedule(static) private(atom,index)
  do cell = 1, ncell
    do atom = 1, natom(cell)
      index = offsets(cell-1) + atom
      force(:,atom,cell) = packed(:,index)
    end do
  end do
  !$omp end parallel do
end subroutine

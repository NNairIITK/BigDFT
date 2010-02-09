subroutine initialize_potential()

  use defs
  use bigdft_forces

  implicit none

  integer :: ierr, nat_test
  integer, pointer  :: typa(:) ! Atomic type
  real(8), pointer  :: posa(:)  ! Working positions of the atoms

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  call bigdft_init(nat_test,typa,posa,box, nproc, iproc)
  call geopt_set_verbosity(0)
  if (nat_test /= NATOMS) stop "NAT"
  ! todo: test nat_test and nat
  type(:) = typa(:)
  pos(:) = posa(:)

end subroutine initialize_potential

subroutine finalise_potential()
  use bigdft_forces

  implicit none  
  
  call bigdft_finalise()
end subroutine finalise_potential

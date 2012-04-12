program biiiiip

  use BigDFT_API

  implicit none

  integer :: ierr, iproc, nproc, i_all, i_stat, infocode
  character(len=60) :: radical, posinp
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  real(gp) :: etot, fnoise
  real(gp), dimension(:,:), pointer :: rxyz
  real(gp), dimension(6) :: strten
  real(gp), dimension(:,:), allocatable :: fxyz

  call bigdft_mpi_init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

  ! Read a possible radical format argument.
  call get_command_argument(1, value = radical, status = ierr)
  if (ierr > 0) then
     write(radical, "(A)") "input"
     write(posinp, "(A)") "posinp"
  else
     write(posinp, "(A)") trim(radical)
  end if

  call standard_inputfile_names(inputs, radical, nproc)
  call read_input_variables(iproc, trim(posinp), inputs, atoms, rxyz)

  call init_restart_objects(iproc, inputs%iacceleration, atoms, rst, "main")
  allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz,'fxyz', "main")

  call call_bigdft(nproc, iproc, atoms, rxyz, inputs, etot, fxyz, strten, fnoise, rst, infocode)

  call free_restart_objects(rst, "main")

  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz', "main")
  i_all=-product(shape(fxyz))*kind(fxyz)
  deallocate(fxyz,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz', "main")
  call deallocate_atoms(atoms, "main") 
  call free_input_variables(inputs)

  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)
  
end program biiiiip

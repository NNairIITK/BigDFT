!!****p* BigDFT/BigDFT
!! FUNCTION
!!  Main program to calculate electronic structures
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program BigDFT

  use module_base
  use module_types
  use module_interfaces
  use ab6_symmetry

  !as a general policy, we'll have "implicit none" by assuming the same
  !name convention as "implicit real(kind=8) (a-h,o-z)"

  implicit none
  character(len=*), parameter :: subname='BigDFT'
  integer :: iproc,nproc,iat,j,i_stat,i_all,ierr,infocode
  integer :: ncount_bigdft
  real(gp) :: etot,sumx,sumy,sumz,fnoise
  logical :: exist_list
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  character(len=50), dimension(:), allocatable :: arr_posinp
  character(len=60) :: filename
  ! atomic coordinates, forces
  real(gp), dimension(:,:), allocatable :: fxyz
  real(gp), dimension(:,:), pointer :: rxyz
  integer :: iconfig,nconfig

  ! Start MPI in parallel version
  !in the case of MPIfake libraries the number of processors is automatically adjusted
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! find out which input files will be used
  inquire(file="list_posinp",exist=exist_list)
  if (exist_list) then
     open(54,file="list_posinp")
     read(54,*) nconfig
     if (nconfig > 0) then 
        !allocation not referenced since memocc count not initialised
        allocate(arr_posinp(1:nconfig))

        do iconfig=1,nconfig
           read(54,*) arr_posinp(iconfig)
        enddo
     else
        nconfig=1
        allocate(arr_posinp(1:1))
        arr_posinp(1)='posinp'
     endif
     close(54)
  else
     nconfig=1
     allocate(arr_posinp(1:1))
     arr_posinp(1)='posinp'
  end if

  do iconfig=1,nconfig

     !welcome screen
     if (iproc==0) call print_logo()

     ! Read all input files.
     call read_input_variables(iproc,trim(arr_posinp(iconfig)), &
          & "input.dft", "input.kpt", "input.geopt", "input.perf", inputs, atoms, rxyz)
     if (iproc == 0) then
        call print_general_parameters(inputs,atoms)
     end if

     !initialize memory counting
     !call memocc(0,iproc,'count','start')

     allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)

     call init_restart_objects(iproc,inputs%iacceleration,atoms,rst,subname)

     !if other steps are supposed to be done leave the last_run to minus one
     !otherwise put it to one
     if (inputs%last_run == -1 .and. inputs%ncount_cluster_x <=1) then
        inputs%last_run = 1
     end if
 
     call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,fnoise,rst,infocode)


     if (inputs%ncount_cluster_x > 1) then
        if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
        ! geometry optimization
        open(unit=16,file='geopt.mon',status='unknown',position='append')
        if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'
        call geopt(nproc,iproc,rxyz,atoms,fxyz,etot,rst,inputs,ncount_bigdft)
        close(16)
        filename=trim('final_'//trim(arr_posinp(iconfig)))
        if (iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,' ')
     end if

     !if there is a last run to be performed do it now before stopping
     if (inputs%last_run == -1) then
        inputs%last_run = 1
        call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,fnoise,rst,infocode)
     end if


     if (iproc == 0) then
        sumx=0.d0
        sumy=0.d0
        sumz=0.d0
        write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
        do iat=1,atoms%nat
           write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
                iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
           sumx=sumx+fxyz(1,iat)
           sumy=sumy+fxyz(2,iat)
           sumz=sumz+fxyz(3,iat)
        enddo
!!$        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
!!$           write(*,'(1x,a)')'the sum of the forces is'
!!$           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
!!$           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
!!$           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
!!$        end if
     endif

     call deallocate_atoms(atoms,subname) 

     call free_restart_objects(rst,subname)

     i_all=-product(shape(rxyz))*kind(rxyz)
     deallocate(rxyz,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz',subname)
     i_all=-product(shape(fxyz))*kind(fxyz)
     deallocate(fxyz,stat=i_stat)
     call memocc(i_stat,i_all,'fxyz',subname)

     call free_input_variables(inputs)

     !finalize memory counting
     call memocc(0,0,'count','stop')

  enddo !loop over iconfig

  deallocate(arr_posinp)

  call MPI_FINALIZE(ierr)

end program BigDFT
!!***

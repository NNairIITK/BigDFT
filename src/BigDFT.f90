!!****p* BigDFT/BigDFT
!! FUNCTION
!!  Main program to calculate electronic structures
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA, UNIBAS
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
!  use minimization, only: parameterminimization 

  !implicit real(kind=8) (a-h,o-z)
  !as a general policy, I will put "implicit none" by assuming the same
  !name convention as "implicit real(kind=8) (a-h,o-z)"
  !such that the implicit statement can be commented at will

  implicit none
  character(len=*), parameter :: subname='BigDFT'
  character(len=20) :: units
  logical :: exists
  integer :: iproc,nproc,iat,ityp,j,i_stat,i_all,ierr,infocode
  integer ::  ncount_bigdft
  real(gp) :: etot,sumx,sumy,sumz,tt
  logical :: exist_list
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  character(len=20), dimension(:), allocatable :: atomnames
  character(len=50), dimension(:), allocatable :: arr_posinp
  character(len=60)  :: filename
  ! atomic coordinates, forces
  real(gp), dimension(:,:), allocatable :: fxyz
  real(gp), dimension(:,:), pointer :: rxyz
  integer :: npr,iam,iconfig,nconfig
  integer  :: nfluct
  real(gp) :: fluctsum

 
  !!!!$      interface
  !!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
  !!!!$        end function omp_get_num_threads
  !!!!$      end interface
  !!!!$      interface
  !!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
 !!!! !$        end function omp_get_thread_num
  !!!!!$      end interface

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

     !initialize memory counting
     call memocc(0,iproc,'count','start')

     !welcome screen
     if (iproc==0) call print_logo()

     !read atomic file
     call read_atomic_file(trim(arr_posinp(iconfig)),iproc,atoms,rxyz)

     allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)

     ! read dft input variables
     call dft_input_variables(iproc,'input.dft',inputs,atoms%symObj)
     !call read_input_variables(iproc,'input.dat',inputs)

     ! read k-points input variables (if given)
     call kpt_input_variables(iproc,'input.kpt',inputs,atoms)

     !read geometry optimsation input variables
     !inquire for the file needed for geometry optimisation
     !if not present, perform a simple geometry optimisation
     inquire(file="input.geopt",exist=exists)
     if (exists) then
        call geopt_input_variables(iproc,'input.geopt',inputs)
     else
        call geopt_input_variables_default(inputs)
     end if

     do iat=1,atoms%nat
        if (atoms%ifrztyp(iat) == 0) then
           call random_number(tt)
           rxyz(1,iat)=rxyz(1,iat)+inputs%randdis*tt
           call random_number(tt)
           rxyz(2,iat)=rxyz(2,iat)+inputs%randdis*tt
           call random_number(tt)
           rxyz(3,iat)=rxyz(3,iat)+inputs%randdis*tt
        end if
     enddo

     !atoms inside the box (this can be insertedindise call_bigdft routine
     do iat=1,atoms%nat
        if (atoms%geocode == 'P') then
           rxyz(1,iat)=modulo(rxyz(1,iat),atoms%alat1)
           rxyz(2,iat)=modulo(rxyz(2,iat),atoms%alat2)
           rxyz(3,iat)=modulo(rxyz(3,iat),atoms%alat3)
        else if (atoms%geocode == 'S') then
           rxyz(1,iat)=modulo(rxyz(1,iat),atoms%alat1)
           rxyz(3,iat)=modulo(rxyz(3,iat),atoms%alat3)
        end if
     end do

     call init_restart_objects(atoms,rst,subname)

     call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,rst,infocode)

     if (inputs%ncount_cluster_x > 1) then
        if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
        ! geometry optimization
        open(unit=16,file='geopt.mon',status='unknown')
        if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'
        call geopt(nproc,iproc,rxyz,atoms,fxyz,etot,rst,inputs,ncount_bigdft)
        filename=trim('relaxed_'//trim(arr_posinp(iconfig)))
        if (iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,' ')
     end if


     if (iproc.eq.0) then
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
        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
           write(*,'(1x,a)')'the sum of the forces is'
           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
        end if
     endif


     !deallocations
     ! TODO, move them into input_variables as a free method for atoms
     i_all=-product(shape(atoms%ifrztyp))*kind(atoms%ifrztyp)
     deallocate(atoms%ifrztyp,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%ifrztyp',subname)
     i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
     deallocate(atoms%iatype,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%iatype',subname)
     i_all=-product(shape(atoms%natpol))*kind(atoms%natpol)
     deallocate(atoms%natpol,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%natpol',subname)
     i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
     deallocate(atoms%atomnames,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%atomnames',subname)
     i_all=-product(shape(atoms%amu))*kind(atoms%amu)
     deallocate(atoms%amu,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%amu',subname)
     if (atoms%symObj >= 0) call ab6_symmetry_free(atoms%symObj)

     call free_restart_objects(rst,subname)

     i_all=-product(shape(rxyz))*kind(rxyz)
     deallocate(rxyz,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz',subname)
     i_all=-product(shape(fxyz))*kind(fxyz)
     deallocate(fxyz,stat=i_stat)
     call memocc(i_stat,i_all,'fxyz',subname)


     !finalize memory counting
     call memocc(0,0,'count','stop')

     if (GPUshare .and. GPUconv) call stop_gpu_sharing()

  enddo !loop over iconfig

  deallocate(arr_posinp)

  call free_input_variables(inputs)

  call MPI_FINALIZE(ierr)

end program BigDFT
!!***

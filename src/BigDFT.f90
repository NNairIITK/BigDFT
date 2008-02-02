!!****p* BigDFT/BigDFT
!! NAME
!!   BigDFT
!!
!! FUNCTION
!!  Main program to calculate electronic structures
!!
!! COPYRIGHT
!!    Copyright (C) 2007 CEA, UNIBAS
!!
!! SOURCE
!!
program BigDFT

  use module_types
  use module_interfaces

  !implicit real(kind=8) (a-h,o-z)
  !as a general policy, I will put "implicit none" by assuming the same
  !name convention as "implicit real(kind=8) (a-h,o-z)"
  !such that the implicit statement can be commented at will

  implicit none
  include 'mpif.h'
  ! For parallel MPI execution set parallel=.true., for serial parallel=.false.
  ! this statement wil be changed by using the MPIfake.f90 file
  include 'parameters.h'
  character(len=20) :: units
  integer :: iproc,nproc,n1,n2,n3,iat,ityp,j,i_stat,i_all,ierr,infocode
  integer :: ncount_cluster
  integer :: norb,norbp
  real(kind=8) :: energy,etot,energyold,beta,sumx,sumy,sumz,tt
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(wavefunctions_descriptors) :: wfd
  character(len=20), dimension(:), allocatable :: atomnames
  ! atomic coordinates, forces
  real(kind=8), dimension(:,:), allocatable :: rxyz,fxyz,rxyz_old
 
  real(kind=8), dimension(:), pointer :: eval
  real(kind=8), dimension(:,:), pointer :: psi


  !$      interface
  !$        integer ( kind=4 ) function omp_get_num_threads ( )
  !$        end function omp_get_num_threads
  !$      end interface
  !$      interface
  !$        integer ( kind=4 ) function omp_get_thread_num ( )
  !$        end function omp_get_thread_num
  !$      end interface

  ! Start MPI in parallel version
  if (parallel) then
     call MPI_INIT(ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  else
     nproc=1
     iproc=0
  endif

  !experimental, to see if it works
  if (nproc == 1) parallel=.false.

  !initialize memory counting
  call memocc(0,iproc,'count','start')

  !$omp parallel private(iam)  shared (npr)
  !$       iam=omp_get_thread_num()
  !$       if (iam.eq.0) npr=omp_get_num_threads()
  !$       write(*,*) 'iproc,iam,npr',iproc,iam,npr
  !$omp end parallel

  !welcome screen
  if (iproc==0) call print_logo()

  !read number of atoms
  open(unit=99,file='posinp',status='old')
  read(99,*) atoms%nat,units

  allocate(rxyz_old(3,atoms%nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz_old))*kind(rxyz_old),'rxyz_old','BigDFT')
  allocate(rxyz(3,atoms%nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','BigDFT')
  allocate(fxyz(3,atoms%nat),stat=i_stat)
  call memocc(i_stat,product(shape(fxyz))*kind(fxyz),'fxyz','BigDFT')

  ! read atomic positions
  call read_atomic_positions(iproc,99,units,atoms,rxyz)

  close(99)

  !new way of reading the input variables, use structures
  call read_input_variables(iproc,inputs)
 
  do iat=1,atoms%nat
     if (.not. atoms%lfrztyp(iat)) then
        call random_number(tt)
        rxyz(1,iat)=rxyz(1,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(2,iat)=rxyz(2,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(3,iat)=rxyz(3,iat)+inputs%randdis*tt
     end if
  enddo

  call call_cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
       psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,inputs,infocode)

  if (inputs%ncount_cluster_x > 1) then
     if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
     ! geometry optimization
     !    betax=2.d0   ! Cinchonidine
     !    betax=4.d0  ! Si H_4
     !   betax=7.5d0  ! silicon systems
     !    betax=10.d0  !  Na_Cl clusters
     ncount_cluster=1
     beta=inputs%betax
     energyold=1.d100
     call conjgrad(parallel,nproc,iproc,atoms,rxyz,etot,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,inputs)
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
     write(*,'(1x,a)')'the sum of the forces is'
     write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
     write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
     write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
  endif

  !deallocations
  i_all=-product(shape(atoms%lfrztyp))*kind(atoms%lfrztyp)
  deallocate(atoms%lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp','BigDFT')
  i_all=-product(shape(atoms%nspinat))*kind(atoms%nspinat)
  deallocate(atoms%nspinat,stat=i_stat)
  call memocc(i_stat,i_all,'nspinat','BigDFT')

  call deallocate_wfd(wfd,'BigDFT')

  i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
  deallocate(atoms%atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames','BigDFT')
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi','BigDFT')
  i_all=-product(shape(eval))*kind(eval)
  deallocate(eval,stat=i_stat)
  call memocc(i_stat,i_all,'eval','BigDFT')
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz','BigDFT')
  i_all=-product(shape(rxyz_old))*kind(rxyz_old)
  deallocate(rxyz_old,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz_old','BigDFT')
  i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
  deallocate(atoms%iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype','BigDFT')
  i_all=-product(shape(fxyz))*kind(fxyz)
  deallocate(fxyz,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz','BigDFT')


  !finalize memory counting
  call memocc(0,0,'count','stop')

  if (parallel) call MPI_FINALIZE(ierr)

 end program BigDFT
 !!***

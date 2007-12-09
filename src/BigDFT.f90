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
  logical :: output_wf,output_grid,calc_tail
  character(len=20) :: units
  character(len=80) :: line
  integer :: iproc,nproc,nat,ntypes,n1,n2,n3,iat,ityp,j,i_stat,i_all,ierr,infocode
  integer :: ncount_cluster
  integer :: norb,norbp
  real(kind=8) :: energy,etot,energyold,beta,sumx,sumy,sumz,tt
  !input variables
  type(input_variables) :: inputs
  type(wavefunctions_descriptors) :: wfd

  logical, dimension(:), allocatable :: lfrztyp
  character(len=6), dimension(:), allocatable :: frzsymb
  character(len=20), dimension(:), allocatable :: atomnames
  ! atomic types
  integer, dimension(:), allocatable :: iatype
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
  read(99,*) nat,units

  allocate(rxyz_old(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz_old))*kind(rxyz_old),'rxyz_old','BigDFT')
  allocate(rxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','BigDFT')
  allocate(iatype(nat),stat=i_stat)
  call memocc(i_stat,product(shape(iatype))*kind(iatype),'iatype','BigDFT')
  allocate(fxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(fxyz))*kind(fxyz),'fxyz','BigDFT')
  allocate(atomnames(100),stat=i_stat) 
  call memocc(i_stat,product(shape(atomnames))*kind(atomnames),'atomnames','BigDFT')

  ! read atomic positions
  call read_atomic_positions(iproc,99,units,nat,ntypes,iatype,atomnames,rxyz)

  close(99)

  !new way of reading the input variables, use structures
  call read_input_variables(iproc,inputs)
 
  !this array is useful for frozen atoms
  !we should modify the way in which they must be entrered
  allocate(lfrztyp(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(lfrztyp))*kind(lfrztyp),'lfrztyp','BigDFT')
  lfrztyp(:)=.false.

!!$  !  Read the first line of "input.dat"
!!$  open(unit=1,file='input.dat',status='old')
!!$  !Only the first line for the main routine (the program)
!!$  ! ngeostep Number of steps of geometry optimisation (default 500)
!!$  ! ampl     Amplitude for random displacement away from input file geometry 
!!$  ! (usually equilibrium geom.)
!!$  ! betax    Geometry optimisation
!!$  read(1,'(a80)')line
!!$  close(unit=1)
!!$  
!!$  allocate(lfrztyp(ntypes),stat=i_stat)
!!$  call memocc(i_stat,product(shape(lfrztyp))*kind(lfrztyp),'lfrztyp','BigDFT')
!!$
!!$  
!!$  read(line,*,iostat=ierror)ngeostep,ampl,betax,nfrztyp
!!$  if (ierror /= 0) then
!!$     read(line,*)ngeostep,ampl,betax
!!$     lfrztyp(:)=.false.
!!$  else
!!$     allocate(frzsymb(nfrztyp),stat=i_stat)
!!$     call memocc(i_stat,product(shape(frzsymb))*kind(frzsymb),'frzsymb','BigDFT')
!!$     read(line,*,iostat=ierror)ngeostep,ampl,betax,nfrztyp,(frzsymb(i),i=1,nfrztyp)
!!$     lfrztyp(:)=.false.
!!$     do ityp=1,nfrztyp
!!$        seek_frozen: do jtyp=1,ntypes
!!$           if (trim(frzsymb(ityp))==trim(atomnames(jtyp))) then
!!$              lfrztyp(jtyp)=.true.
!!$              exit seek_frozen
!!$           end if
!!$        end do seek_frozen
!!$     end do
!!$     i_all=-product(shape(frzsymb))*kind(frzsymb)
!!$     deallocate(frzsymb,stat=i_stat)
!!$     call memocc(i_stat,i_all,'frzsymb','BigDFT')
!!$  end if

  do iat=1,nat
     if (.not. lfrztyp(iatype(iat))) then
        call random_number(tt)
        rxyz(1,iat)=rxyz(1,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(2,iat)=rxyz(2,iat)+inputs%randdis*tt
        call random_number(tt)
        rxyz(3,iat)=rxyz(3,iat)+inputs%randdis*tt
     end if
  enddo

  call call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
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
     call conjgrad(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,rxyz,etot,fxyz,&
          psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,inputs)
  end if

  if (iproc.eq.0) then
     sumx=0.d0
     sumy=0.d0
     sumz=0.d0
     write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
     do iat=1,nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
             iat,trim(atomnames(iatype(iat))),(fxyz(j,iat),j=1,3)
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
  !the lfrztyp must be restored in the same way as before
  i_all=-product(shape(lfrztyp))*kind(lfrztyp)
  deallocate(lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp','BigDFT')

  call deallocate_wfd(wfd,'BigDFT')

  i_all=-product(shape(atomnames))*kind(atomnames)
  deallocate(atomnames,stat=i_stat)
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
  i_all=-product(shape(iatype))*kind(iatype)
  deallocate(iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype','BigDFT')
  i_all=-product(shape(fxyz))*kind(fxyz)
  deallocate(fxyz,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz','BigDFT')


  !finalize memory counting
  call memocc(0,0,'count','stop')

  if (parallel) call MPI_FINALIZE(ierr)

 end program BigDFT

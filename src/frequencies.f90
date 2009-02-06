!!****p* BigDFT/frequencies
!! FUNCTION
!!  Calculate vibrational frequencies
!!
!! COPYRIGHT
!!    Copyright (C) 2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program frequencies

  use module_base
  use module_types
  use module_interfaces

  implicit none
  character(len=*), parameter :: subname='BigDFT'
  character(len=20) :: units
  character(len=2) :: cc
  integer :: iproc,nproc,iat,ityp,i,j,i_stat,i_all,ierr,infocode
  real(gp) :: etot,etot_m,etot_p,sumx,sumy,sumz,tt,alat,alpha,dd
  !input variables
  type(atoms_data) :: atoms
  type(input_variables) :: inputs
  type(restart_objects) :: rst
  character(len=20), dimension(:), allocatable :: atomnames
  ! atomic coordinates, forces
  real(gp), dimension(:,:), allocatable :: rxyz,fxyz,rpos,fpos_m,fpos_p
  integer npr,iam
 
  !$      interface
  !$        integer ( kind=4 ) function omp_get_num_threads ( )
  !$        end function omp_get_num_threads
  !$      end interface
  !$      interface
  !$        integer ( kind=4 ) function omp_get_thread_num ( )
  !$        end function omp_get_thread_num
  !$      end interface

  ! Start MPI in parallel version
  !in the case of MPIfake libraries the number of processors is automatically adjusted
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  !initialize memory counting
  call memocc(0,iproc,'count','start')

!**********Commented out by Alexey, 15.11.2008************************************************  
!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$       write(*,*) 'iproc,iam,npr',iproc,iam,npr
!$omp end parallel
!*********************************************************************************************

  !welcome screen
  if (iproc==0) call print_logo()

  !read number of atoms
  open(unit=99,file='posinp',status='old')
  read(99,*) atoms%nat,atoms%units

  allocate(rxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)
  allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz,'fxyz',subname)

  ! read atomic positions
  call read_atomic_positions(iproc,99,atoms,rxyz)

  close(99)

  ! read input variables, use structures
  call read_input_variables(iproc,'input.dat',inputs)
 
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

  ! atoms inside the box (this can be inserted inside call_bigdft routine)
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

  if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode

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

  allocate(rpos(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rpos,'rpos',subname)
  allocate(fpos_m(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fpos_m,'fpos_m',subname)
  allocate(fpos_p(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fpos_p,'fpos_p',subname)

! Move to alpha*hgrid
  alpha=1.d0/real(64,kind(1.d0))
  
  if (iproc ==0 ) then
     write(*,"(1x,a)") '=Frequencies calculation='
     write(10,*) "#step=",alpha*inputs%hgrid
     write(10,*) 0,"--",etot,alpha*inputs%hgrid,fxyz
  end if

  do iat=1,atoms%nat
     if (atoms%lfrztyp(iat)) then
        if (iproc==0) write(*,"(1x,a,i0,a)") '=F:The atom ',iat,' is frozen.'
        cycle
     end if

     do i=1,3
        if (i==1) then
           alat=atoms%alat1
           cc(2:2)='x'
        else if (i==2) then
           alat=atoms%alat2
           cc(2:2)='y'
        else
           alat=atoms%alat3
           cc(2:2)='z'
        end if
        do j=-1,1,2
           if (j==-1) then
              cc(1:1)='-'
           else
              cc(1:1)='+'
           end if
           !We copy atomic positions
           rpos=rxyz
           !!Displacement
           dd=real(j,gp)*alpha*inputs%hgrid
           if (iproc==0) write(*,"(1x,a,i0,a,a)") '=F:Move the atom ',iat,' in the direction ',cc

           if (atoms%geocode == 'P') then
              rpos(i,iat)=modulo(rxyz(i,iat)+dd,alat)
           else if (atoms%geocode == 'S') then
              rpos(i,iat)=modulo(rxyz(i,iat)+dd,alat)
           else
              rpos(i,iat)=rxyz(i,iat)+dd
           end if
           inputs%inputPsiId=1
           inputs%output_grid=0
           inputs%output_wf=.false.
           if (j==-1) then
              call call_bigdft(nproc,iproc,atoms,rpos,inputs,etot_m,fpos_m,rst,infocode)
               if (iproc==0) write(10,*) iat,cc,etot_m-etot,fpos_m-fxyz
           else
              call call_bigdft(nproc,iproc,atoms,rpos,inputs,etot_p,fpos_p,rst,infocode)
               if (iproc==0) write(10,*) iat,cc,etot_p-etot,fpos_p-fxyz
           end if
        end do
     end do

  end do


  !deallocations
  i_all=-product(shape(atoms%lfrztyp))*kind(atoms%lfrztyp)
  deallocate(atoms%lfrztyp,stat=i_stat)
  call memocc(i_stat,i_all,'lfrztyp',subname)
  i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
  deallocate(atoms%iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype',subname)
  i_all=-product(shape(atoms%natpol))*kind(atoms%natpol)
  deallocate(atoms%natpol,stat=i_stat)
  call memocc(i_stat,i_all,'natpol',subname)
  i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
  deallocate(atoms%atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames',subname)

  call free_restart_objects(rst,subname)

  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  i_all=-product(shape(fxyz))*kind(fxyz)
  deallocate(fxyz,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz',subname)

  i_all=-product(shape(rpos))*kind(rpos)
  deallocate(rpos,stat=i_stat)
  call memocc(i_stat,i_all,'rpos',subname)
  i_all=-product(shape(fpos_m))*kind(fpos_m)
  deallocate(fpos_m,stat=i_stat)
  call memocc(i_stat,i_all,'fpos_m',subname)
  i_all=-product(shape(fpos_p))*kind(fpos_p)
  deallocate(fpos_p,stat=i_stat)
  call memocc(i_stat,i_all,'fpos_p',subname)


  !finalize memory counting
  call memocc(0,0,'count','stop')

  if (nproc > 1) call MPI_FINALIZE(ierr)

END PROGRAM frequencies
!!***

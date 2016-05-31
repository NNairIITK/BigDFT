!> @file
!! Test program
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test the overlapgeneral routine
program driver
  use bigdft_run
  use module_base
  use module_types
  use module_interfaces
  use sparsematrix_base, only: deallocate_sparse_matrix, matrices_null, allocate_matrices, deallocate_matrices, &
                               SPARSE_FULL, sparsematrix_malloc, assignment(=)
  use sparsematrix_init, only: sparse_matrix_init_fake, check_symmetry, generate_random_symmetric_sparsity_pattern
  use sparsematrix, only: compress_matrix, uncompress_matrix, extract_taskgroup_inplace, &
                          write_matrix_compressed
  use matrix_operations, only: overlapPowerGeneral, deviation_from_unity_parallel
  use yaml_output
  implicit none

  ! Variables
  integer :: iproc, nproc
!$  integer :: omp_get_num_threads
  integer,parameter :: itype=1
  character(len=1),parameter :: jobz='v', uplo='l'
  integer,parameter :: n=64
  real(kind=8) :: val, max_error, mean_error
  real(kind=8),dimension(:,:),allocatable :: ovrlp, ovrlp2
  integer :: norb, nseg, nvctr, iorb, jorb, iorder, power, blocksize, icheck, imode

  logical :: file_exists, symmetric, perform_check, optional_parameters
  type(orbitals_data) :: orbs
  type(sparse_matrix) :: smat_A, smat_B, smat_test
  type(matrices) :: mat_A
  type(matrices),dimension(1) :: inv_mat_B
  character(len=*),parameter :: filename='inputdata.fake'
  integer :: nconfig, ierr, iseg, iiorb
!!  integer :: lwork
  integer, dimension(4) :: mpi_info
  character(len=60) :: run_id
  integer,parameter :: ncheck=33
  !!integer,dimension(:,:),allocatable :: keyg_tmp
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2

  integer :: ncount1, ncount_rate, ncount_max, ncount2, nn, ios
!! integer :: i, j, start
  real(kind=4) :: tr0, tr1
  real(kind=8) :: time, time2, tt
  real(kind=8),dimension(:),allocatable :: tmparr
!! real(kind=8) :: tmp
  real :: rn
  real(kind=8), external :: ddot, dnrm2
  logical, parameter :: timer_on=.false.        !time the different methods
  logical, parameter :: ortho_check=.false.     !check deviation from orthonormality of input overlap matrix
  logical, parameter :: print_matrices=.true.  !output calculated matrices

  ! Initialize
  call f_lib_initialize()

  call bigdft_init()!mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=bigdft_mpi%iproc!mpi_info(1)
  nproc=bigdft_mpi%nproc!mpi_info(2)


  if (iproc==0) then
      call yaml_comment('Program to check the overlapPowergeneral routine',hfill='/')
      call yaml_map('reading from file',filename)
  end  if

  ! Open file for reading
  inquire(file=filename, exist=file_exists)
  if_file_exists: if (file_exists) then
      ! Read the basis quantities
      optional_parameters=.true.
      open(unit=1, file=filename)
      read(1,*) norb
      read(1,*) nseg
      read(1,*) nvctr
      ! the following lines are optional
      read(1,*,iostat=ios) imode
      if (ios/=0) optional_parameters=.false.
      read(1,*,iostat=ios) iorder
      if (ios/=0) optional_parameters=.false.
      read(1,*,iostat=ios) power
      if (ios/=0) optional_parameters=.false.
      read(1,*,iostat=ios) blocksize
      if (ios/=0) optional_parameters=.false.
      close(unit=1)
  else
      stop 'file does not exist!'
  end if if_file_exists

  if (iproc==0) then
      call yaml_sequence_open('parameters for this test')
      call yaml_map('number of rows and columns',norb)
      call yaml_map('number of segments',nseg)
      call yaml_map('number of non-zero elements',nvctr)
      call yaml_sequence_close
  end if


  ! Fake initialization of the orbitals_data type
  call orbs_init_fake(iproc, nproc, norb, orbs)

  ! Fake initialization of the sparse_matrix type
  call sparse_matrix_init_fake(iproc, nproc, bigdft_mpi%mpi_comm, norb, nseg, nvctr, smat_A)
  call sparse_matrix_init_fake(iproc, nproc, bigdft_mpi%mpi_comm, norb, nseg, nvctr, smat_B)



  symmetric = check_symmetry(smat_A)

  !!if (iproc==0) then
  !!    do iseg=1,smat_A%nseg
  !!        write(*,*) smat_A%keyg(:,iseg)
  !!    end do
  !!    do iseg=1,smat_A%nvctr
  !!        write(*,*) smat_A%orb_from_index(:,iseg)
  !!    end do
  !!end if

  ! Initialize an overlap matrix
  allocate(ovrlp(orbs%norb,orbs%norb))
  if (orbs%norb<30) then
      do iorb=1,orbs%norb
          do jorb=iorb,orbs%norb
              val = 2.d-1*(sin(real((iorb-1)*n+jorb,kind=8)))**2
              if (jorb/=iorb) then
                  ovrlp(jorb,iorb) = val
                  ovrlp(iorb,jorb) = val
              else
                  val = val + 1.d0
                  ovrlp(jorb,iorb) = val
              end if
          end do
      end do
  !DEBUG
  !else if (orbs%norb==984.or..true.) then
  !    start=241
  !    open(100)
  !    do iorb=1,984
  !        do jorb=1,984
  !            read(100,*) i,j,tmp
  !            if (iorb<orbs%norb+start.and.jorb<orbs%norb+start.and.iorb>=start.and.jorb>=start) &
  !                 ovrlp(jorb-start+1,iorb-start+1)=tmp
  !        end do
  !    end do
  !    close(100)
  !END DEBUG
  else
      ! above approach has problems for testing larger matrices
      allocate(ovrlp2(orbs%norb,orbs%norb))
      ! randomly generate vectors
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              call random_number(rn)
              ovrlp2(jorb,iorb)=2.0d0*real(rn,kind=8)-1.0d0
          end do
          tt=dnrm2(orbs%norb, ovrlp2(1,iorb), 1)
          call dscal(orbs%norb, 1/tt, ovrlp2(1,iorb), 1)
      end do

      ! calculate overlap from random vectors
      do iorb=1,orbs%norb
          do jorb=iorb,orbs%norb
             ovrlp(jorb,iorb)=ddot(orbs%norb,ovrlp2(1,iorb),1,ovrlp2(1,jorb),1)
             ovrlp(iorb,jorb)=ovrlp(jorb,iorb)
          end do
      end do
      deallocate(ovrlp2)
  end if

  !!lwork=100*orbs%norb
  !!allocate(work(lwork))
  !!allocate(eval(orbs%norb))
  !!ovrlp2=ovrlp
  !!call dsyev('v', 'l', orbs%norb, ovrlp2, orbs%norb, eval, work, lwork, info)
  !!do iseg=1,orbs%norb
  !!    write(*,*) iseg, eval(iseg)
  !!end do

  mat_A = matrices_null()
  inv_mat_B(1) = matrices_null()

  call allocate_matrices(smat_A, allocate_full=.true., matname='mat_A', mat=mat_A)
  call vcopy(orbs%norb**2, ovrlp(1,1), 1, mat_A%matrix(1,1,1), 1)
  call compress_matrix(iproc, nproc, smat_A, inmat=mat_A%matrix, outmat=mat_A%matrix_compr)
  call allocate_matrices(smat_B, allocate_full=.true., matname='inv_mat_B', mat=inv_mat_B(1))
  ! uncomment for sparse and dense modes to be testing the same matrix
  !call uncompress_matrix(iproc, smat_A)

  if (print_matrices.and.iproc==0) call write_matrix_compressed('initial matrix', smat_A, mat_A)



  ! Check of the overlap manipulation routine

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  !!keyg_tmp=f_malloc((/2,smat_A%nseg/))
  !!do iseg=1,smat_A%nseg
  !!    iorb=smat_A%keyg(1,iseg)
  !!    iiorb=mod(iorb-1,norb)+1
  !!    !!write(*,*) 'iorb, iiorb', iorb, iiorb
  !!    keyg_tmp(1,iseg)=iiorb
  !!    iorb=smat_A%keyg(2,iseg)
  !!    iiorb=mod(iorb-1,norb)+1
  !!    !!write(*,*) 'iorb, iiorb', iorb, iiorb
  !!    keyg_tmp(2,iseg)=iiorb
  !!end do

  if (ortho_check) call deviation_from_unity_parallel(iproc, nproc, bigdft_mpi%mpi_comm, &
                        orbs%norb, orbs%norb, 0, ovrlp, smat_A, max_error, mean_error)
  if (ortho_check.and.iproc==0) call yaml_map('max deviation from unity',max_error)
  if (ortho_check.and.iproc==0) call yaml_map('mean deviation from unity',mean_error)
  if (iproc==0) call yaml_comment('starting the checks',hfill='=')

  if (.not.optional_parameters) then
      ! do all checks
      nn=ncheck
  else
      ! do only the check which was specified
      nn=1
  end if

  do icheck=1,nn
      if (.not.optional_parameters) then
          ! get the default parameters
          call get_parameters()
      end if
      if (iproc==0) then
          call yaml_comment('check:'//yaml_toa(icheck,fmt='(i5)'),hfill='-')
          call yaml_map('check number',icheck)
          call yaml_map('imode',imode)
          call yaml_map('iorder',iorder)
          call yaml_map('power',power)
          call yaml_newline()
      end if
      if (.not.symmetric .and. imode==SPARSE .and. iorder==0) then
          perform_check=.false.
      else
          perform_check=.true.
      end if
      if (iproc==0) call yaml_map('Can perform this test',perform_check)
      if (.not.perform_check) cycle
      if (imode==DENSE) then
          call vcopy(orbs%norb**2, ovrlp(1,1), 1, mat_A%matrix(1,1,1), 1)
          if (timer_on) call cpu_time(tr0)
          if (timer_on) call system_clock(ncount1,ncount_rate,ncount_max)
          call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, iorder, 1, (/power/), blocksize, &
               imode, ovrlp_smat=smat_A, inv_ovrlp_smat=smat_B, ovrlp_mat=mat_A, inv_ovrlp_mat=inv_mat_B, &
               check_accur=.true., max_error=max_error, mean_error=mean_error)
          if (timer_on) call cpu_time(tr1)
          if (timer_on) call system_clock(ncount2,ncount_rate,ncount_max)
          if (timer_on) time=real(tr1-tr0,kind=8)
          if (timer_on) time2=dble(ncount2-ncount1)/dble(ncount_rate)
          call compress_matrix(iproc, nproc, smat_B, inmat=inv_mat_B(1)%matrix, outmat=inv_mat_B(1)%matrix_compr)
      else if (imode==SPARSE) then
          call vcopy(orbs%norb**2, ovrlp(1,1), 1, mat_A%matrix(1,1,1), 1)
          call compress_matrix(iproc, nproc, smat_A, inmat=mat_A%matrix, outmat=mat_A%matrix_compr)
          if (timer_on) call cpu_time(tr0)
          if (timer_on) call system_clock(ncount1,ncount_rate,ncount_max)
          tmparr = sparsematrix_malloc(smat_A,iaction=SPARSE_FULL,id='tmparr')
          call vcopy(smat_A%nvctr*smat_A%nspin, mat_A%matrix_compr(1), 1, tmparr(1), 1)
          call extract_taskgroup_inplace(smat_A, mat_A)
          call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, iorder, 1, (/power/), blocksize, &
               imode, ovrlp_smat=smat_A, inv_ovrlp_smat=smat_B, ovrlp_mat=mat_A, inv_ovrlp_mat=inv_mat_B, &
               check_accur=.true., max_error=max_error, mean_error=mean_error)
          call vcopy(smat_A%nvctr*smat_A%nspin, tmparr(1), 1, mat_A%matrix_compr(1), 1)
          call f_free(tmparr)
               !!foe_nseg=smat_A%nseg, foe_kernel_nsegline=smat_A%nsegline, &
               !!foe_istsegline=smat_A%istsegline, foe_keyg=smat_A%keyg)
           !if (iorder==0) call compress_matrix(iproc, smat_B)
          if (timer_on) call cpu_time(tr1)
          if (timer_on) call system_clock(ncount2,ncount_rate,ncount_max)
          if (timer_on) time=real(tr1-tr0,kind=8)
          if (timer_on) time2=dble(ncount2-ncount1)/dble(ncount_rate)
      end if
      if (print_matrices.and.iproc==0) call write_matrix_compressed('final result', smat_B, inv_mat_B(1))
      if (iproc==0) call yaml_map('Max error of the result',max_error)
      if (iproc==0) call yaml_map('Mean error of the result',mean_error)
      if (timer_on.and.iproc==0) call yaml_map('time taken (cpu)',time)
      if (timer_on.and.iproc==0) call yaml_map('time taken (system)',time2)
  end do

  if (iproc==0) call yaml_comment('checks finished',hfill='=')

  !!call f_free(keyg_tmp)

  call deallocate_orbitals_data(orbs)
  call deallocate_sparse_matrix(smat_A)
  call deallocate_sparse_matrix(smat_B)
  call deallocate_matrices(mat_A)
  call deallocate_matrices(inv_mat_B(1))

  deallocate(ovrlp)

  call bigdft_finalize(ierr)

  call f_lib_finalize()


  !!call mpi_barrier(mpi_comm_world, ierr)
  !!call mpi_finalize(ierr)

  contains

    subroutine get_parameters()
      select case (icheck)
      case (1)
          imode = 2 ; iorder=0 ; power= -2 ; blocksize=-1
      case (2)
          imode = 2 ; iorder=1 ; power= -2 ; blocksize=-1
      case (3)
          imode = 2 ; iorder=6 ; power= -2 ; blocksize=-1
      case (4)
          imode = 2 ; iorder=-1 ; power= -2 ; blocksize=-1
      case (5)
          imode = 2 ; iorder=-6 ; power= -2 ; blocksize=-1
      case (6)
          imode = 2 ; iorder=0 ; power= 1 ; blocksize=-1
      case (7)
          imode = 2 ; iorder=1 ; power= 1 ; blocksize=-1
      case (8)
          imode = 2 ; iorder=6 ; power= 1 ; blocksize=-1
      case (9)
          imode = 2 ; iorder=-1 ; power= 1 ; blocksize=-1
      case (10)
          imode = 2 ; iorder=-6 ; power= 1 ; blocksize=-1
      case (11)
          imode = 2 ; iorder=0 ; power= 2 ; blocksize=-1
      case (12)
          imode = 2 ; iorder=1 ; power= 2 ; blocksize=-1
      case (13)
          imode = 2 ; iorder=6 ; power= 2 ; blocksize=-1
      case (14)
          imode = 2 ; iorder=-1 ; power= 2 ; blocksize=-1
      case (15)
          imode = 2 ; iorder=-6 ; power= 2 ; blocksize=-1
      case (16)
          imode = 1 ; iorder=0 ; power= -2 ; blocksize=-1
      case (17)
          imode = 1 ; iorder=1 ; power= -2 ; blocksize=-1
      case (18)
          imode = 1 ; iorder=6 ; power= -2 ; blocksize=-1
      case (19)
          imode = 1 ; iorder=-1 ; power= -2 ; blocksize=-1
      case (20)
          imode = 1 ; iorder=-6 ; power= -2 ; blocksize=-1
      case (21)
          imode = 1 ; iorder=0 ; power= 1 ; blocksize=-1
      case (22)
          imode = 1 ; iorder=1 ; power= 1 ; blocksize=-1
      case (23)
          imode = 1 ; iorder=6 ; power= 1 ; blocksize=-1
      case (24)
          imode = 1 ; iorder=-1 ; power= 1 ; blocksize=-1
      case (25)
          imode = 1 ; iorder=-6 ; power= 1 ; blocksize=-1
      case (26)
          imode = 1 ; iorder=0 ; power= 2 ; blocksize=-1
      case (27)
          imode = 1 ; iorder=1 ; power= 2 ; blocksize=-1
      case (28)
          imode = 1 ; iorder=6 ; power= 2 ; blocksize=-1
      case (29)
          imode = 1 ; iorder=-1 ; power= 2 ; blocksize=-1
      case (30)
          imode = 1 ; iorder=-6 ; power= 2 ; blocksize=-1
      case (31)
          imode = 1 ; iorder=1025 ; power= -2 ; blocksize=-1
      case (32)
          imode = 1 ; iorder=1025 ; power= 2 ; blocksize=-1
      case (33)
          imode = 1 ; iorder=1025 ; power= 2 ; blocksize=-1
      case default
          stop 'wrong icheck'
      end select
    end subroutine get_parameters
end program driver



!> Fake initialization of the orbitals_data type
subroutine orbs_init_fake(iproc, nproc, norb, orbs)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: iproc, nproc, norb
  type(orbitals_data),intent(out) :: orbs

  ! Nullify the data type
  call nullify_orbitals_data(orbs)

  ! Allocate the arrays which will be needed
  call allocate_arrays()

  ! Initialize the relevant data. First the one from the input.
  orbs%norb = norb

  ! Now intialize the remaning fields if they will be needed.
  orbs%norb_par(:,0) = norb_par_init()
  orbs%norbp = orbs%norb_par(iproc,0)
  orbs%isorb_par = isorb_par_init()
  orbs%isorb = orbs%isorb_par(iproc)

  !!write(*,*) 'iproc, orbs%norb', iproc, orbs%norb
  !!write(*,*) 'iproc, orbs%norbp', iproc, orbs%norbp
  !!write(*,*) 'iproc, orbs%isorb', iproc, orbs%isorb
  !!write(*,*) 'iproc, orbs%norb_par', iproc, orbs%norb_par
  !!write(*,*) 'iproc, orbs%isorb_par', iproc, orbs%isorb_par


  contains

    subroutine allocate_arrays
      implicit none
      orbs%norb_par=f_malloc_ptr((/0.to.nproc-1,0.to.0/),id='orbs%norb_par')
      orbs%isorb_par=f_malloc_ptr(0.to.nproc-1,id='orbs%isorb_par')
    end subroutine allocate_arrays

    function norb_par_init() result(norb_par)
      integer,dimension(0:nproc-1) :: norb_par
      real(kind=8) :: tt
      integer :: ii, jproc
      tt=real(norb,kind=8)/real(nproc,kind=8)
      ii=floor(tt)
      do jproc=0,nproc-1
          norb_par(jproc)=ii
      end do
      ii=norb-nproc*ii
      do jproc=0,ii-1
          norb_par(jproc)=norb_par(jproc)+1
      end do
      if (sum(norb_par)/=orbs%norb) stop 'sum(norb_par)/=orbs%norb'
    end function norb_par_init

    function isorb_par_init() result(isorb_par)
      integer,dimension(0:nproc-1) :: isorb_par
      integer :: jproc
      isorb_par(0)=0
      do jproc=1,nproc-1
          isorb_par(jproc)=isorb_par(jproc-1)+orbs%norb_par(jproc-1,0)
      end do
    end function isorb_par_init
 
end subroutine orbs_init_fake

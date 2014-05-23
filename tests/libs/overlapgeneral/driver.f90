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
  use module_base
  use module_types
  use module_interfaces
  use sparsematrix, only: compress_matrix
  use sparsematrix_base, only: deallocate_sparse_matrix
  use yaml_output
  implicit none

  ! Variables
  integer :: iproc, nproc
!$  integer :: omp_get_num_threads
  integer,parameter :: itype=1
  character(len=1),parameter :: jobz='v', uplo='l'
  integer,parameter :: n=64
  real(kind=8) :: val, error
  real(kind=8),dimension(:,:),allocatable :: ovrlp, ovrlp2
  integer :: norb, nseg, nvctr, iorb, jorb, iorder, power, blocksize, icheck, imode
  logical :: file_exists, symmetric, check_symmetry, perform_check
  type(orbitals_data) :: orbs
  type(sparse_matrix) :: smat_A, smat_B
  character(len=*),parameter :: filename='inputdata.fake'
  integer :: nconfig, ierr, iseg, iiorb
!!  integer :: lwork
  integer, dimension(4) :: mpi_info
  character(len=60) :: run_id
  integer,parameter :: ncheck=18
  integer,dimension(:,:),allocatable :: keyg_tmp
  integer,parameter :: SPARSE=1
  integer,parameter :: DENSE=2



  ! Initialize
  call f_lib_initialize()
  call bigdft_init(mpi_info,nconfig,run_id,ierr)
  !just for backward compatibility
  iproc=mpi_info(1)
  nproc=mpi_info(2)


  if (iproc==0) then
      call yaml_comment('Program to check the overlapPowergeneral routine',hfill='/')
      call yaml_map('reading from file',filename)
  end  if

  ! Open file for reading
  inquire(file=filename, exist=file_exists)
  if_file_exists: if (file_exists) then
      ! Read the basis quantities
      open(unit=1, file=filename)
      read(1,*) norb
      read(1,*) nseg
      read(1,*) nvctr
  else
      stop 'file does not exist!'
  end if if_file_exists

  if (iproc==0) then
      call yaml_open_sequence('parameters for this test')
      call yaml_map('number of rows and columns',norb)
      call yaml_map('number of segments',nseg)
      call yaml_map('number of non-zero elements',nvctr)
      call yaml_close_sequence
  end if


  ! Fake initialization of the orbitals_data type
  call orbs_init_fake(iproc, nproc, norb, orbs)

  ! Fake initialization of the sparse_matrix type
  call sparse_matrix_init_fake(norb, orbs%norbp, orbs%isorb, nseg, nvctr, smat_A)
  call sparse_matrix_init_fake(norb, orbs%norbp, orbs%isorb, nseg, nvctr, smat_B)

  symmetric = check_symmetry(norb, smat_A)

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
  allocate(ovrlp2(orbs%norb,orbs%norb))
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
  !!lwork=100*orbs%norb
  !!allocate(work(lwork))
  !!allocate(eval(orbs%norb))
  !!ovrlp2=ovrlp
  !!call dsyev('v', 'l', orbs%norb, ovrlp2, orbs%norb, eval, work, lwork, info)
  !!do iseg=1,orbs%norb
  !!    write(*,*) iseg, eval(iseg)
  !!end do


  call vcopy(orbs%norb**2, ovrlp(1,1), 1, smat_A%matrix(1,1), 1)
  call compress_matrix(iproc, smat_A)
  if (iproc==0) call write_matrix_compressed('initial matrix', smat_A)


  ! Check of the overlap manipulation routine

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  keyg_tmp=f_malloc((/2,smat_A%nseg/))
  do iseg=1,smat_A%nseg
      iorb=smat_A%keyg(1,iseg)
      iiorb=mod(iorb-1,norb)+1
      !!write(*,*) 'iorb, iiorb', iorb, iiorb
      keyg_tmp(1,iseg)=iiorb
      iorb=smat_A%keyg(2,iseg)
      iiorb=mod(iorb-1,norb)+1
      !!write(*,*) 'iorb, iiorb', iorb, iiorb
      keyg_tmp(2,iseg)=iiorb
  end do

  if (iproc==0) call yaml_comment('starting the checks',hfill='=')

  do icheck=1,ncheck
      !if (icheck==1 .or. icheck==4 .or. icheck==10 .or. icheck==11 .or.  icheck==13 .or. icheck==14 .or. icheck==16 .or. icheck==17) cycle
      call get_parameters()
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
          call vcopy(orbs%norb**2, ovrlp(1,1), 1, smat_A%matrix(1,1), 1)
          call overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, norb, orbs, &
               imode, check_accur=.true., ovrlp=smat_A%matrix, inv_ovrlp=smat_B%matrix, error=error)
          call compress_matrix(iproc, smat_B)
      else if (imode==SPARSE) then
          call vcopy(orbs%norb**2, ovrlp(1,1), 1, smat_A%matrix(1,1), 1)
          call compress_matrix(iproc, smat_A)
          call overlapPowerGeneral(iproc, nproc, iorder, power, blocksize, norb, orbs, &
               imode, check_accur=.true., error=error, &
               ovrlp_smat=smat_A, inv_ovrlp_smat=smat_B)!!, &
               !!foe_nseg=smat_A%nseg, foe_kernel_nsegline=smat_A%nsegline, &
               !!foe_istsegline=smat_A%istsegline, foe_keyg=smat_A%keyg)
           !if (iorder==0) call compress_matrix(iproc, smat_B)
      end if
      if (iproc==0) call write_matrix_compressed('final result', smat_B)
      if (iproc==0) call yaml_map('error of the result',error)
  end do

  if (iproc==0) call yaml_comment('checks finished',hfill='=')

  call f_free(keyg_tmp)

  call deallocate_orbitals_data(orbs, 'driver')
  call deallocate_sparse_matrix(smat_A, 'driver')
  call deallocate_sparse_matrix(smat_B, 'driver')


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
          imode = 2 ; iorder=0 ; power= 1 ; blocksize=-1
      case (5)
          imode = 2 ; iorder=1 ; power= 1 ; blocksize=-1
      case (6)
          imode = 2 ; iorder=6 ; power= 1 ; blocksize=-1
      case (7)
          imode = 2 ; iorder=0 ; power= 2 ; blocksize=-1
      case (8)
          imode = 2 ; iorder=1 ; power= 2 ; blocksize=-1
      case (9)
          imode = 2 ; iorder=6 ; power= 2 ; blocksize=-1
      case (10)
          imode = 1 ; iorder=0 ; power= -2 ; blocksize=-1
      case (11)
          imode = 1 ; iorder=1 ; power= -2 ; blocksize=-1
      case (12)
          imode = 1 ; iorder=6 ; power= -2 ; blocksize=-1
      case (13)
          imode = 1 ; iorder=0 ; power= 1 ; blocksize=-1
      case (14)
          imode = 1 ; iorder=1 ; power= 1 ; blocksize=-1
      case (15)
          imode = 1 ; iorder=6 ; power= 1 ; blocksize=-1
      case (16)
          imode = 1 ; iorder=0 ; power= 2 ; blocksize=-1
      case (17)
          imode = 1 ; iorder=1 ; power= 2 ; blocksize=-1
      case (18)
          imode = 1 ; iorder=6 ; power= 2 ; blocksize=-1
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


!> Fake initialization of the sparse_matrix type
subroutine sparse_matrix_init_fake(norb, norbp, isorb, nseg, nvctr, smat)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null
  use sparsematrix_init, only: init_sparse_matrix_matrix_multiplication
  implicit none

  ! Calling arguments
  integer,intent(in) :: norb, norbp, isorb, nseg, nvctr
  type(sparse_matrix) :: smat

  ! Local variables
  integer,dimension(:),allocatable :: nvctr_per_segment

  ! Some checks whether the arguments are reasonable
  if (nseg > nvctr) stop 'sparse matrix would have more segments than elements'
  if (nseg < norb) stop 'sparse matrix would have less segments than lines'
  if (nvctr > norb**2) stop 'sparse matrix would contain more elements than the dense one'

  ! Nullify the data type
  smat = sparse_matrix_null()

  ! Initialize the relevant data. First the one from the input.
  smat%nvctr = nvctr
  smat%nseg = nseg
  smat%nfvctr = norb

  ! Now some default values
  smat%parallel_compression=0
  smat%store_index=.false.

  ! Allocate the arrays which will be needed
  call allocate_arrays()

  ! Auxiliary array
  nvctr_per_segment = nvctr_per_segment_init()

  ! Now intialize the remaning fields if they will be needed.
  smat%nsegline = nsegline_init()
  smat%istsegline = istsegline_init ()
  smat%keyv = keyv_init()
  smat%keyg = keyg_init()
  call init_orbs_from_index(smat)

  call f_free(nvctr_per_segment)

  ! Initialize the parameters for the spare matrix matrix multiplication
  call init_sparse_matrix_matrix_multiplication(norb, norbp, isorb, smat%nseg, &
       smat%nsegline, smat%istsegline, smat%keyg, smat)

  !!if (iproc==0) then
  !!    do jorb=1,norb
  !!        write(*,*) 'jorb, nsegline, istsegline', jorb, smat%nsegline(jorb), smat%istsegline(jorb) 
  !!    end do
  !!    do jseg=1,smat%nseg
  !!        write(*,*) 'keyv, keyg', smat%keyv(jseg), smat%keyg(:,jseg)
  !!    end do
  !!end if

  contains

    subroutine allocate_arrays
      implicit none
      smat%nsegline=f_malloc_ptr(norb,id='smat%nsegline')
      smat%istsegline=f_malloc_ptr(norb,id='smat%istsegline')
      nvctr_per_segment=f_malloc(nseg,id='nvctr_per_segment')
      smat%keyv=f_malloc_ptr(nseg,id='smat%keyv')
      smat%keyg=f_malloc_ptr((/2,nseg/),id='smat%keyg')
      smat%matrix_compr=f_malloc_ptr(smat%nvctr,id='smat%matrix_compr')
      smat%matrix=f_malloc_ptr((/norb,norb/),id='smat%matrix')
    end subroutine allocate_arrays

    function nsegline_init() result(nsegline)
      integer,dimension(norb) :: nsegline
      real(kind=8) :: tt
      integer :: ii, jorb
      ! Distribute segments evenly among the lines
      tt=real(nseg,kind=8)/real(norb,kind=8)
      ii=floor(tt)
      do jorb=1,norb
          nsegline(jorb)=ii
      end do
      ii=nseg-norb*ii
      do jorb=1,ii
          nsegline(jorb)=nsegline(jorb)+1
      end do
    end function nsegline_init

    function istsegline_init() result(istsegline)
      integer,dimension(norb) :: istsegline
      integer :: jorb
      istsegline(1)=1
      do jorb=2,norb
          istsegline(jorb)=istsegline(jorb-1)+smat%nsegline(jorb-1)
      end do
    end function istsegline_init


    function nvctr_per_segment_init() result(nvctr_per_segment)
      integer,dimension(nseg) :: nvctr_per_segment
      real(kind=8) :: tt
      integer :: ii, jseg
      ! Distribute the elements evenly among the segments
      tt=real(nvctr,kind=8)/real(nseg,kind=8)
      ii=floor(tt)
      do jseg=1,nseg
          nvctr_per_segment(jseg)=ii
      end do
      ii=nvctr-nseg*ii
      do jseg=1,ii
          nvctr_per_segment(jseg)=nvctr_per_segment(jseg)+1
      end do
      if (sum(nvctr_per_segment)/=smat%nvctr) stop 'sum(nvctr_per_segment)/=smat%nvctr'
    end function nvctr_per_segment_init


    function keyv_init() result(keyv)
      integer,dimension(smat%nseg) :: keyv
      integer :: jseg
      keyv(1)=1
      do jseg=2,nseg
          keyv(jseg)=keyv(jseg-1)+nvctr_per_segment(jseg-1)
      end do
    end function keyv_init


    function keyg_init() result(keyg)
      integer,dimension(2,smat%nseg) :: keyg
      integer :: jorb, nempty, jseg, jjseg, ii, j, ist, itot, istart, iend, idiag
      integer :: idist_start, idist_end, ilen
      integer,dimension(:),allocatable :: nempty_arr
      real(kind=8) :: tt
      integer,parameter :: DECREASE=1, INCREASE=2

      itot=1
      do jorb=1,norb
          ! Number of empty elements
          nempty=norb
          do jseg=1,smat%nsegline(jorb)
              jjseg=smat%istsegline(jorb)+jseg-1
              nempty=nempty-nvctr_per_segment(jjseg)
          end do
          if (nempty<0) then
              write(*,*) 'ERROR: nemtpy < 0; reduce number of elements'
              stop
          end if
          ! Number of empty elements between the elements
          allocate(nempty_arr(0:smat%nsegline(jorb)))
          tt=real(nempty,kind=8)/real(smat%nsegline(jorb)+1,kind=8)
          ii=floor(tt)
          do j=0,smat%nsegline(jorb)
              nempty_arr(j)=ii
          end do
          ii=nempty-(smat%nsegline(jorb)+1)*ii
          do j=0,ii-1
              nempty_arr(j)=nempty_arr(j)+1
          end do
          ! Check that the diagonal element is not in an empty region. If so,
          ! shift the elements.
          idiag=(jorb-1)*norb+jorb
          adjust_empty: do
              ist=nempty_arr(0)
              do jseg=1,smat%nsegline(jorb)
                  jjseg=smat%istsegline(jorb)+jseg-1
                  istart=itot+ist
                  iend=istart+nvctr_per_segment(jjseg)-1
                  if (istart<=idiag .and. idiag<=iend) exit adjust_empty
                  ! Determine the distance to the start / end of the segment
                  idist_start=abs(idiag-istart)
                  idist_end=abs(idiag-iend)
                  !!if (j==1 .and. idiag<istart) then
                  !!    ! Diagonal element is before the first segment, 
                  !!    ! so decrease the first empty region
                  !!    iaction=DECREASE
                  !!end if
                  !!if (j==smat%nsegline(jorb) .and. idiag>iend) then
                  !!    ! Diagonal element is after the last segment, 
                  !!    ! so increase the first empty region
                  !!    iaction=INCREASE
                  !!end if
                  ist=ist+nvctr_per_segment(jjseg)
                  ist=ist+nempty_arr(jseg)
              end do
              ! If one arrives here, the diagonal element was in an empty
              ! region. Determine whether it was close to the start or end of a
              ! segment.
              if (istart==iend) then
                  ! Segment has only length one
                  if (istart<idiag) then
                      ! Incrase the first empty region and increase the last one
                      nempty_arr(0)=nempty_arr(0)+1
                      nempty_arr(smat%nsegline(jorb))=nempty_arr(smat%nsegline(jorb))-1
                  else
                      ! Decrase the first empty region and increase the last one
                      nempty_arr(0)=nempty_arr(0)-1
                      nempty_arr(smat%nsegline(jorb))=nempty_arr(smat%nsegline(jorb))+1
                  end if
              else if (idist_start<=idist_end) then
                  ! Closer to the start, so decrase the first empty region and increase the last one
                  nempty_arr(0)=nempty_arr(0)-1
                  nempty_arr(smat%nsegline(jorb))=nempty_arr(smat%nsegline(jorb))+1
              else 
                  ! Closer to the end, so increase the first empty region and decrease the last one
                  nempty_arr(0)=nempty_arr(0)+1
                  nempty_arr(smat%nsegline(jorb))=nempty_arr(smat%nsegline(jorb))-1
              end if
          end do adjust_empty

          ! Now fill the keys
          ist=nempty_arr(0)
          do jseg=1,smat%nsegline(jorb)
              jjseg=smat%istsegline(jorb)+jseg-1
              istart=itot+ist
              iend=istart+nvctr_per_segment(jjseg)-1
              keyg(1,jjseg)=istart
              keyg(2,jjseg)=iend
              ist=ist+nvctr_per_segment(jjseg)
              ist=ist+nempty_arr(jseg)
          end do
          itot=itot+ist
          deallocate(nempty_arr)
      end do

      ! Check that the total number is correct
      itot=0
      do jseg=1,smat%nseg
          ilen=keyg(2,jseg)-keyg(1,jseg)+1
          if (ilen/=nvctr_per_segment(jseg)) stop 'ilen/=nvctr_per_segment(jseg)'
          if (jseg/=smat%nseg) then
              if (ilen/=(smat%keyv(jseg+1)-smat%keyv(jseg))) stop 'ilen/=(smat%keyv(jseg+1)-smat%keyv(jseg))'
          else
              if (ilen/=(smat%nvctr+1-smat%keyv(jseg))) stop 'ilen/=(smat%nvctr+1-smat%keyv(jseg))'
          end if
          itot=itot+ilen
      end do
      if (itot/=smat%nvctr) stop 'itot/=smat%nvctr'
    end function keyg_init


    subroutine init_orbs_from_index(sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(inout) :: sparsemat

      ! local variables
      integer :: ind, iseg, segn, iorb, jorb
      character(len=*),parameter :: subname='init_orbs_from_index'

      sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')

      ind = 0
      do iseg = 1, sparsemat%nseg
         do segn = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            ind=ind+1
            iorb = (segn - 1) / sparsemat%nfvctr + 1
            jorb = segn - (iorb-1)*sparsemat%nfvctr
            sparsemat%orb_from_index(1,ind) = jorb
            sparsemat%orb_from_index(2,ind) = iorb
         end do
      end do

    end subroutine init_orbs_from_index

end subroutine sparse_matrix_init_fake


subroutine write_matrix_compressed(message, smat)
  use yaml_output
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  character(len=*),intent(in) :: message
  type(sparse_matrix),intent(in) :: smat

  ! Local variables
  integer :: iseg, ilen, istart, iend, i, iorb, jorb

  !!call yaml_open_sequence(trim(message))
  !!do iseg=1,smat%nseg
  !!    call yaml_sequence(advance='no')
  !!    ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
  !!    call yaml_open_map(flow=.true.)
  !!    call yaml_map('segment',iseg)
  !!    istart=smat%keyv(iseg)
  !!    iend=smat%keyv(iseg)+ilen
  !!    call yaml_map('values',smat%matrix_compr(istart:iend))
  !!    call yaml_close_map()
  !!    call yaml_newline()
  !!end do
  !!call yaml_close_sequence()

  call yaml_open_sequence(trim(message))
  do iseg=1,smat%nseg
      call yaml_sequence(advance='no')
      ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
      call yaml_open_map(flow=.true.)
      call yaml_map('segment',iseg)
      call yaml_open_sequence('elements')
      istart=smat%keyv(iseg)
      iend=smat%keyv(iseg)+ilen-1
      do i=istart,iend
          call yaml_newline()
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          iorb=smat%orb_from_index(1,i)
          jorb=smat%orb_from_index(2,i)
          call yaml_map('coordinates',(/jorb,iorb/))
          call yaml_map('value',smat%matrix_compr(i))
          call yaml_close_map()
      end do
      call yaml_close_sequence()
      !call yaml_map('values',smat%matrix_compr(istart:iend))
      call yaml_close_map()
      call yaml_newline()
  end do
  call yaml_close_sequence()

end subroutine write_matrix_compressed


function check_symmetry(norb, smat)
  use module_base
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: norb
  type(sparse_matrix),intent(in) :: smat
  logical :: check_symmetry

  ! Local variables
  integer :: i, iorb, jorb
  logical,dimension(:,:),allocatable :: lgrid

  lgrid=f_malloc((/norb,norb/),id='lgrid')
  lgrid=.false.

  do i=1,smat%nvctr
      iorb=smat%orb_from_index(1,i)
      jorb=smat%orb_from_index(2,i)
      lgrid(jorb,iorb)=.true.
  end do

  check_symmetry=.true.
  do iorb=1,norb
      do jorb=1,norb
          if (lgrid(jorb,iorb) .and. .not.lgrid(iorb,jorb)) then
              check_symmetry=.false.
          end if
      end do
  end do

  call f_free(lgrid)

end function check_symmetry


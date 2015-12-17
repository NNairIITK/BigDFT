!> @file
!!  File defining the structures to deal with the sparse matrices
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the basic operations with sparse matrices (initialization)
module sparsematrix_init
  use module_base
  use sparsematrix_base
  implicit none

  private

  !> Public routines
  public :: init_sparse_matrix
  public :: matrixindex_in_compressed
  public :: matrixindex_in_compressed_lowlevel
  public :: init_matrix_taskgroups
  public :: check_local_matrix_extents
  public :: read_ccs_format
  public :: ccs_to_sparsebigdft
  public :: ccs_values_to_bigdft
  public :: bigdft_to_sparsebigdft
  public :: distribute_columns_on_processes_simple
  public :: redistribute
  public :: get_modulo_array
  public :: sparsebigdft_to_ccs
  public :: ccs_to_sparsebigdft_short
  public :: init_matrixindex_in_compressed_fortransposed
  public :: distribute_on_threads


contains



  integer function matrixindex_in_compressed(sparsemat, iorb, jorb, init_, n_)
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: sparsemat
      integer,intent(in) :: iorb, jorb
      !> The optional arguments should only be used for initialization purposes
      !! if one is sure what one is doing. Might be removed later.
      logical,intent(in),optional :: init_
      integer,intent(in),optional :: n_
    
      ! Local variables
      integer :: ispin, iiorb, jjorb
      !integer :: ii
      logical :: lispin, ljspin, init

      if (present(init_)) then
          init = init_
      else
          init = .false.
      end if

      ! Use the built-in function and return, without any check. Can be used for initialization purposes.
      if (init) then
          if (.not.present(n_)) stop 'matrixindex_in_compressed: n_ must be present if init_ is true'
          matrixindex_in_compressed = compressed_index_fn(iorb, jorb, n_, sparsemat)
          return
      end if

      !ii=(jorb-1)*sparsemat%nfvctr+iorb
      !ispin=(ii-1)/sparsemat%nfvctr**2+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)

      ! Determine in which "spin matrix" this entry is located
      lispin = (iorb>sparsemat%nfvctr)
      ljspin = (jorb>sparsemat%nfvctr)
      if (any((/lispin,ljspin/))) then
          if (all((/lispin,ljspin/))) then
              ! both indices belong to the second spin matrix
              ispin=2
          else
              ! there seems to be a mix of the spin matrices
              write(*,*) 'iorb, jorb, nfvctr', iorb, jorb, sparsemat%nfvctr
              stop 'matrixindex_in_compressed: problem in determining spin'
          end if
      else
          ! both indices belong to the first spin matrix
          ispin=1
      end if
      iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
      jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
      
    
      if (sparsemat%store_index) then
          ! Take the value from the array
          matrixindex_in_compressed = sparsemat%matrixindex_in_compressed_arr(iiorb,jjorb)
      else
          ! Recalculate the value
          matrixindex_in_compressed = compressed_index_fn(iiorb, jjorb, sparsemat%nfvctr, sparsemat)
      end if

      ! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
      if (ispin==2) then
          if (matrixindex_in_compressed/=0) then
              matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
          end if
      end if
    
!!$    contains

!!$      ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
!!$      integer function compressed_index_fn(irow, jcol, norb, sparsemat)
!!$        use sparsematrix_base
!!$        implicit none
!!$      
!!$        ! Calling arguments
!!$        integer,intent(in) :: irow, jcol, norb
!!$        type(sparse_matrix),intent(in) :: sparsemat
!!$      
!!$        ! Local variables
!!$        integer(kind=8) :: ii, istart, iend, norb8
!!$        integer :: iseg
!!$      
!!$        norb8 = int(norb,kind=8)
!!$        ii = int((jcol-1),kind=8)*norb8+int(irow,kind=8)
!!$      
!!$        iseg=sparsemat%istsegline(jcol)
!!$        do
!!$            istart = int((sparsemat%keyg(1,2,iseg)-1),kind=8)*norb8 + &
!!$                     int(sparsemat%keyg(1,1,iseg),kind=8)
!!$            if (ii<istart) then
!!$                compressed_index_fn=0
!!$                return
!!$            end if
!!$            iend = int((sparsemat%keyg(2,2,iseg)-1),kind=8)*norb8 + &
!!$                   int(sparsemat%keyg(2,1,iseg),kind=8)
!!$            !if (ii>=istart .and. ii<=iend) then
!!$            if (ii<=iend) then
!!$                ! The matrix element is in sparsemat segment
!!$                 compressed_index_fn = sparsemat%keyv(iseg) + int(ii-istart,kind=4)
!!$                return
!!$            end if
!!$            iseg=iseg+1
!!$            if (iseg>sparsemat%nseg) exit
!!$        end do
!!$      
!!$        ! Not found
!!$        compressed_index_fn=0
!!$      
!!$      end function compressed_index_fn
    end function matrixindex_in_compressed

    ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
    integer function compressed_index_fn(irow, jcol, norb, sparsemat)
      use sparsematrix_base
      implicit none
      ! Calling arguments
      integer,intent(in) :: irow, jcol, norb
      type(sparse_matrix),intent(in) :: sparsemat

      ! Local variables
      integer(kind=8) :: ii, istart, iend, norb8
      integer :: iseg

      norb8 = int(norb,kind=8)
      ii = int((jcol-1),kind=8)*norb8+int(irow,kind=8)

      iseg=sparsemat%istsegline(jcol)
      do
         istart = int((sparsemat%keyg(1,2,iseg)-1),kind=8)*norb8 + &
              int(sparsemat%keyg(1,1,iseg),kind=8)
         if (ii<istart) then
            compressed_index_fn=0
            return
         end if
         iend = int((sparsemat%keyg(2,2,iseg)-1),kind=8)*norb8 + &
              int(sparsemat%keyg(2,1,iseg),kind=8)
         !if (ii>=istart .and. ii<=iend) then
         if (ii<=iend) then
            ! The matrix element is in sparsemat segment
            compressed_index_fn = sparsemat%keyv(iseg) + int(ii-istart,kind=4)
            return
         end if
         iseg=iseg+1
         if (iseg>sparsemat%nseg) exit
      end do

      ! Not found
      compressed_index_fn=0

    end function compressed_index_fn



    !> Does the same as matrixindex_in_compressed, but has different
    ! arguments (at lower level) and is less optimized
    integer function matrixindex_in_compressed_lowlevel(irow, jcol, norb, nseg, keyv, keyg, istsegline) result(micf)
      implicit none

      ! Calling arguments
      integer,intent(in) :: irow, jcol, norb, nseg
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      integer,dimension(norb),intent(in) :: istsegline

      ! Local variables
      integer(kind=8) :: ii, istart, iend, norb8
      integer :: iseg

      norb8=int(norb,kind=8)
      ii = int((jcol-1),kind=8)*norb8+int(irow,kind=8)

      !do iseg=1,nseg
      iseg=istsegline(jcol)
      do
          istart = int((keyg(1,2,iseg)-1),kind=8)*norb8 + &
                   int(keyg(1,1,iseg),kind=8)
          !iend = int((keyg(2,2,iseg)-1),kind=8)*int(norb,kind=8) + &
          !       int(keyg(2,1,iseg),kind=8)
          !if (ii>=istart .and. ii<=iend) then
          if (ii<istart) then
              micf=0
              return
          end if
          !if (ii>=istart) then
             iend = int((keyg(2,2,iseg)-1),kind=8)*norb8 + &
                    int(keyg(2,1,iseg),kind=8)
             if (ii<=iend) then
                ! The matrix element is in this segment
                micf = keyv(iseg) + int(ii-istart,kind=4)
                return
             end if
          !end if
          iseg = iseg + 1
          if (iseg>nseg) exit
      end do

      ! Not found
      micf=0

    end function matrixindex_in_compressed_lowlevel


    !!!integer function matrixindex_in_compressed2(sparsemat, iorb, jorb, init_, n_)
    !!!  use sparsematrix_base, only: sparse_matrix
    !!!  implicit none
    !!!
    !!!  ! Calling arguments
    !!!  type(sparse_matrix),intent(in) :: sparsemat
    !!!  integer,intent(in) :: iorb, jorb
    !!!  !> The optional arguments should only be used for initialization purposes
    !!!  !! if one is sure what one is doing. Might be removed later.
    !!!  logical,intent(in),optional :: init_
    !!!  integer,intent(in),optional :: n_
    !!!
    !!!  ! Local variables
    !!!  integer :: ii, ispin, iiorb, jjorb
    !!!  logical :: lispin, ljspin, init

    !!!  if (present(init_)) then
    !!!      init = init_
    !!!  else
    !!!      init = .false.
    !!!  end if

    !!!  ! Use the built-in function and return, without any check. Can be used for initialization purposes.
    !!!  if (init) then
    !!!      if (.not.present(n_)) stop 'matrixindex_in_compressed2: n_ must be present if init_ is true'
    !!!      matrixindex_in_compressed2 = compressed_index_fn(iorb, jorb, n_, sparsemat)
    !!!      return
    !!!  end if

    !!!  !ii=(jorb-1)*sparsemat%nfvctr+iorb
    !!!  !ispin=(ii-1)/sparsemat%nfvctr**2+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)

    !!!  ! Determine in which "spin matrix" this entry is located
    !!!  lispin = (iorb>sparsemat%nfvctr)
    !!!  ljspin = (jorb>sparsemat%nfvctr)
    !!!  if (any((/lispin,ljspin/))) then
    !!!      if (all((/lispin,ljspin/))) then
    !!!          ! both indices belong to the second spin matrix
    !!!          ispin=2
    !!!      else
    !!!          ! there seems to be a mix up the spin matrices
    !!!          stop 'matrixindex_in_compressed2: problem in determining spin'
    !!!      end if
    !!!  else
    !!!      ! both indices belong to the first spin matrix
    !!!      ispin=1
    !!!  end if
    !!!  iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
    !!!  jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
    !!!
    !!!  if (sparsemat%store_index) then
    !!!      ! Take the value from the array
    !!!      matrixindex_in_compressed2 = sparsemat%matrixindex_in_compressed_arr(iiorb,jjorb)
    !!!  else
    !!!      ! Recalculate the value
    !!!      matrixindex_in_compressed2 = compressed_index_fn(iiorb, jjorb, sparsemat%nfvctr, sparsemat)
    !!!  end if

    !!!  ! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
    !!!  if (ispin==2) then
    !!!      matrixindex_in_compressed2 = matrixindex_in_compressed2 + sparsemat%nvctrp_tg
    !!!  end if
    !!!
    !!!contains

    !!!  ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
    !!!  integer function compressed_index_fn(irow, jcol, norb, sparsemat)
    !!!    implicit none
    !!!  
    !!!    ! Calling arguments
    !!!    integer,intent(in) :: irow, jcol, norb
    !!!    type(sparse_matrix),intent(in) :: sparsemat
    !!!  
    !!!    ! Local variables
    !!!    integer(kind=8) :: ii, istart, iend
    !!!    integer :: iseg
    !!!  
    !!!    ii = int((jcol-1),kind=8)*int(norb,kind=8)+int(irow,kind=8)
    !!!  
    !!!    iseg=sparsemat%istsegline(jcol)
    !!!    do
    !!!        istart = int((sparsemat%keyg(1,2,iseg)-1),kind=8)*int(norb,kind=8) + &
    !!!                 int(sparsemat%keyg(1,1,iseg),kind=8)
    !!!        iend = int((sparsemat%keyg(2,2,iseg)-1),kind=8)*int(norb,kind=8) + &
    !!!               int(sparsemat%keyg(2,1,iseg),kind=8)
    !!!        if (ii>=istart .and. ii<=iend) then
    !!!            ! The matrix element is in sparsemat segment
    !!!             compressed_index_fn = sparsemat%keyv(iseg) + int(ii-istart,kind=4)
    !!!            return
    !!!        end if
    !!!        iseg=iseg+1
    !!!        if (iseg>sparsemat%nseg) exit
    !!!        if (ii<istart) then
    !!!            compressed_index_fn=0
    !!!            return
    !!!        end if
    !!!    end do
    !!!  
    !!!    ! Not found
    !!!    compressed_index_fn=0
    !!!  
    !!!  end function compressed_index_fn
    !!!end function matrixindex_in_compressed2







    !!subroutine init_sparse_matrix_matrix_multiplication(iproc, nproc, norb, norbp, isorb, nseg, &
    !!           nsegline, istsegline, keyv, keyg, sparsemat)
    !!  use yaml_output
    !!  implicit none

    !!  ! Calling arguments
    !!  integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseg
    !!  integer,dimension(norb),intent(in) :: nsegline, istsegline
    !!  integer,dimension(nseg),intent(in) :: keyv
    !!  integer,dimension(2,2,nseg),intent(in) :: keyg
    !!  type(sparse_matrix),intent(inout) :: sparsemat

    !!  integer :: ierr, jproc, iorb, jjproc, iiorb, nseq_min, nseq_max, iseq, ind, ii, iseg, ncount
    !!  integer,dimension(:),allocatable :: nseq_per_line, norb_par_ideal, isorb_par_ideal
    !!  integer,dimension(:,:),allocatable :: istartend_dj, istartend_mm
    !!  integer,dimension(:,:),allocatable :: temparr
    !!  real(kind=8) :: rseq, rseq_ideal, tt, ratio_before, ratio_after
    !!  logical :: printable
    !!  real(kind=8),dimension(:),allocatable :: rseq_per_line

    !!  ! Calculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
    !!  ! the default partitioning of the matrix columns.
    !!  call get_nout(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)
    !!  nseq_per_line = f_malloc0(norb,id='nseq_per_line')
    !!  call determine_sequential_length(norb, norbp, isorb, nseg, &
    !!       nsegline, istsegline, keyg, sparsemat, &
    !!       sparsemat%smmm%nseq, nseq_per_line)
    !!  if (nproc>1) call mpiallred(nseq_per_line(1), norb, mpi_sum, bigdft_mpi%mpi_comm)
    !!  rseq=real(sparsemat%smmm%nseq,kind=8) !real to prevent integer overflow
    !!  if (nproc>1) call mpiallred(rseq, 1, mpi_sum, bigdft_mpi%mpi_comm)

    !!  rseq_per_line = f_malloc(norb,id='rseq_per_line')
    !!  do iorb=1,norb
    !!      rseq_per_line(iorb) = real(nseq_per_line(iorb),kind=8)
    !!  end do


    !!  norb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
    !!  isorb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
    !!  ! Assign the columns of the matrix to the processes such that the load
    !!  ! balancing is optimal
    !!  ! First the default initializations
    !!  !!!!norb_par_ideal(:)=0
    !!  !!!!isorb_par_ideal(:)=norb
    !!  !!!!rseq_ideal = rseq/real(nproc,kind=8)
    !!  !!!!jjproc=0
    !!  !!!!tt=0.d0
    !!  !!!!iiorb=0
    !!  !!!!isorb_par_ideal(0)=0
    !!  !!!!do iorb=1,norb
    !!  !!!!    iiorb=iiorb+1
    !!  !!!!    tt=tt+real(nseq_per_line(iorb),kind=8)
    !!  !!!!    if (tt>=real(jjproc+1,kind=8)*rseq_ideal .and. jjproc/=nproc-1) then
    !!  !!!!        norb_par_ideal(jjproc)=iiorb
    !!  !!!!        isorb_par_ideal(jjproc+1)=iorb
    !!  !!!!        jjproc=jjproc+1
    !!  !!!!        iiorb=0
    !!  !!!!    end if
    !!  !!!!end do
    !!  !!!!norb_par_ideal(jjproc)=iiorb
    !!  rseq_ideal = rseq/real(nproc,kind=8)
    !!  call redistribute(nproc, norb, rseq_per_line, rseq_ideal, norb_par_ideal)
    !!  isorb_par_ideal(0) = 0
    !!  do jproc=1,nproc-1
    !!      isorb_par_ideal(jproc) = isorb_par_ideal(jproc-1) + norb_par_ideal(jproc-1)
    !!  end do


    !!  ! some checks
    !!  if (sum(norb_par_ideal)/=norb) stop 'sum(norb_par_ideal)/=norb'
    !!  if (isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb) stop 'isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb'

    !!  ! Copy the values
    !!  sparsemat%smmm%nfvctrp=norb_par_ideal(iproc)
    !!  sparsemat%smmm%isfvctr=isorb_par_ideal(iproc)


    !!  ! Get the load balancing
    !!  nseq_min = sparsemat%smmm%nseq
    !!  if (nproc>1) call mpiallred(nseq_min, 1, mpi_min, bigdft_mpi%mpi_comm)
    !!  nseq_max = sparsemat%smmm%nseq
    !!  if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
    !!  if (nseq_min>0) then
    !!      ratio_before = real(nseq_max,kind=8)/real(nseq_min,kind=8)
    !!      printable=.true.
    !!  else
    !!      printable=.false.
    !!  end if


    !!  ! Realculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
    !!  ! the optimized partitioning of the matrix columns.
    !!  call get_nout(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)
    !!  call determine_sequential_length(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
    !!       nsegline, istsegline, keyg, sparsemat, &
    !!       sparsemat%smmm%nseq, nseq_per_line)

    !!  ! Get the load balancing
    !!  nseq_min = sparsemat%smmm%nseq
    !!  if (nproc>1) call mpiallred(nseq_min, 1, mpi_min, bigdft_mpi%mpi_comm)
    !!  nseq_max = sparsemat%smmm%nseq
    !!  if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
    !!  ! Not necessary to set the printable flag (if nseq_min was zero before it should be zero here as well)
    !!  if (nseq_min>0) then
    !!      ratio_after = real(nseq_max,kind=8)/real(nseq_min,kind=8)
    !!      if (.not.printable) stop 'this should not happen (sparsematrix)'
    !!  else
    !!      if (printable) stop 'this should not happen (sparsematrix)'
    !!  end if
    !!  if (iproc==0) then
    !!      if (printable) then
    !!          call yaml_map('sparse matmul load balancing naive / optimized',(/ratio_before,ratio_after/),fmt='(f4.2)')
    !!      else
    !!          call yaml_map('sparse matmul load balancing naive / optimized','printing not meaningful (division by zero)')
    !!      end if
    !!  end if
    !!  

    !!  call f_free(nseq_per_line)
    !!  call f_free(rseq_per_line)

    !!  call allocate_sparse_matrix_matrix_multiplication(nproc, norb, nseg, nsegline, istsegline, sparsemat%smmm)


    !!  ! Calculate some auxiliary variables
    !!  temparr = f_malloc0((/0.to.nproc-1,1.to.2/),id='isfvctr_par')
    !!  temparr(iproc,1) = sparsemat%smmm%isfvctr
    !!  temparr(iproc,2) = sparsemat%smmm%nfvctrp
    !!  if (nproc>1) then
    !!      call mpiallred(temparr(0,1), 2*nproc,  mpi_sum, bigdft_mpi%mpi_comm)
    !!  end if
    !!  call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, sparsemat%nseg, sparsemat%nvctr, &
    !!       temparr(0,1), temparr(0,2), sparsemat%istsegline, sparsemat%keyv, &
    !!       sparsemat%smmm%isvctr_mm, sparsemat%smmm%nvctrp_mm, sparsemat%smmm%isvctr_mm_par, sparsemat%smmm%nvctr_mm_par)
    !!  call f_free(temparr)

    !!  sparsemat%smmm%nseg=nseg
    !!  call vcopy(norb, nsegline(1), 1, sparsemat%smmm%nsegline(1), 1)
    !!  call vcopy(norb, istsegline(1), 1, sparsemat%smmm%istsegline(1), 1)
    !!  call init_onedimindices_new(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
    !!       nsegline, istsegline, keyg, &
    !!       sparsemat, sparsemat%smmm%nout, sparsemat%smmm%onedimindices)
    !!  call get_arrays_for_sequential_acces(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
    !!       nsegline, istsegline, keyg, sparsemat, &
    !!       sparsemat%smmm%nseq, sparsemat%smmm%ivectorindex)
    !!  call init_sequential_acces_matrix(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
    !!       nsegline, istsegline, keyg, sparsemat, sparsemat%smmm%nseq, &
    !!       sparsemat%smmm%indices_extract_sequential)

    !!  ! This array gives the starting and ending indices of the submatrix which
    !!  ! is used by a given MPI task
    !!  if (sparsemat%smmm%nseq>0) then
    !!      sparsemat%smmm%istartend_mm(1) = sparsemat%nvctr
    !!      sparsemat%smmm%istartend_mm(2) = 1
    !!      do iseq=1,sparsemat%smmm%nseq
    !!          ind=sparsemat%smmm%indices_extract_sequential(iseq)
    !!          sparsemat%smmm%istartend_mm(1) = min(sparsemat%smmm%istartend_mm(1),ind)
    !!          sparsemat%smmm%istartend_mm(2) = max(sparsemat%smmm%istartend_mm(2),ind)
    !!      end do
    !!  else
    !!      sparsemat%smmm%istartend_mm(1)=sparsemat%nvctr+1
    !!      sparsemat%smmm%istartend_mm(2)=sparsemat%nvctr
    !!  end if

    !!  ! Determine to which segments this corresponds
    !!  do iseg=1,sparsemat%nseg
    !!      if (sparsemat%keyv(iseg)>=sparsemat%smmm%istartend_mm(1)) then
    !!          sparsemat%smmm%istartendseg_mm(1)=iseg
    !!          exit
    !!      end if
    !!  end do
    !!  do iseg=sparsemat%nseg,1,-1
    !!      if (sparsemat%keyv(iseg)<=sparsemat%smmm%istartend_mm(2)) then
    !!          sparsemat%smmm%istartendseg_mm(2)=iseg
    !!          exit
    !!      end if
    !!  end do

    !!  istartend_mm = f_malloc0((/1.to.2,0.to.nproc-1/),id='istartend_mm')
    !!  istartend_mm(1:2,iproc) = sparsemat%smmm%istartend_mm(1:2)
    !!  if (nproc>1) then
    !!      call mpiallred(istartend_mm(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm)
    !!  end if

    !!  ! Partition the entire matrix in disjoint submatrices
    !!  istartend_dj = f_malloc((/1.to.2,0.to.nproc-1/),id='istartend_dj')
    !!  istartend_dj(1,0) = istartend_mm(1,0)
    !!  do jproc=1,nproc-1
    !!      ind = (istartend_mm(2,jproc-1)+istartend_mm(1,jproc))/2
    !!      ! check that this is inside the segment of istartend_mm(:,jproc)
    !!      ind = max(ind,istartend_mm(1,jproc))
    !!      ind = min(ind,istartend_mm(2,jproc))
    !!      ! check that this is not smaller than the beginning of the previous chunk
    !!      ind = max(ind,istartend_dj(1,jproc-1))+1
    !!      ! check that this is not outside the total matrix size (may happen if there are more processes than matrix columns)
    !!      ind = min(ind,sparsemat%nvctr)
    !!      istartend_dj(1,jproc) = ind
    !!      istartend_dj(2,jproc-1) = istartend_dj(1,jproc)-1
    !!  end do
    !!  istartend_dj(2,nproc-1) = istartend_mm(2,nproc-1)
    !!  !if (iproc==0) write(*,'(a,100(2i7,3x))') 'istartend_mm',istartend_mm
    !!  !if (iproc==0) write(*,'(a,100(2i7,3x))') 'istartend_dj',istartend_dj

    !!  ! Some checks
    !!  if (istartend_dj(1,0)/=1) stop 'istartend_dj(1,0)/=1'
    !!  if (istartend_dj(2,nproc-1)/=sparsemat%nvctr) stop 'istartend_dj(2,nproc-1)/=sparsemat%nvctr'
    !!  ii = 0
    !!  do jproc=0,nproc-1
    !!      ncount = istartend_dj(2,jproc)-istartend_dj(1,jproc) + 1
    !!      if (ncount<0) stop 'ncount<0'
    !!      ii = ii + ncount
    !!      if (ii<0) stop 'init_sparse_matrix_matrix_multiplication: ii<0'
    !!      if (jproc>0) then
    !!          if (istartend_dj(1,jproc)/=istartend_dj(2,jproc-1)+1) stop 'istartend_dj(1,jproc)/=istartend_dj(2,jproc-1)'
    !!      end if
    !!  end do
    !!  if (ii/=sparsemat%nvctr) stop 'init_sparse_matrix_matrix_multiplication: ii/=sparsemat%nvctr'

    !!  ! Keep the values of its own task
    !!  sparsemat%smmm%istartend_mm_dj(1) = istartend_dj(1,iproc)
    !!  sparsemat%smmm%istartend_mm_dj(2) = istartend_dj(2,iproc)


    !!  call f_free(norb_par_ideal)
    !!  call f_free(isorb_par_ideal)
    !!  call f_free(istartend_mm)
    !!  call f_free(istartend_dj)
    !!end subroutine init_sparse_matrix_matrix_multiplication



    subroutine init_sparse_matrix_matrix_multiplication_new(iproc, nproc, norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyv, keyg, sparsemat)
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(inout) :: sparsemat

      integer :: ierr, jproc, iorb, jjproc, iiorb, iseq, ind, ii, iseg, ncount
      integer :: iiseg, i, iel, ilen_seg, ist_seg, iend_seg, ispt, iline, icolumn, iseg_start
      integer,dimension(:),allocatable :: nseq_per_line, norb_par_ideal, isorb_par_ideal, nout_par, nseq_per_pt
      integer,dimension(:,:),allocatable :: istartend_dj, istartend_mm
      integer,dimension(:,:),allocatable :: temparr
      real(kind=8) :: rseq, rseq_ideal, ratio_before, ratio_after
      !real(kind=8) :: tt
      logical :: printable
      real(kind=8),dimension(2) :: rseq_max, rseq_average
      real(kind=8),dimension(:),allocatable :: rseq_per_line

      call f_routine(id='init_sparse_matrix_matrix_multiplication_new')


      ! Calculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
      ! the default partitioning of the matrix columns.
      call get_nout(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)


      ! Determine ispt
      ispt = get_offset(iproc, nproc, sparsemat%smmm%nout)
      !!nout_par = f_malloc0(0.to.nproc-1,id='ist_par')
      !!nout_par(iproc) = sparsemat%smmm%nout
      !!call mpiallred(nout_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
      !!ispt = 0
      !!do jproc=0,iproc-1
      !!    ispt = ispt + nout_par(jproc)
      !!end do
      !!call f_free(nout_par)

      nseq_per_line = f_malloc(norb,id='nseq_per_line')
      !!call determine_sequential_length(norb, norbp, isorb, nseg, &
      !!     nsegline, istsegline, keyg, sparsemat, &
      !!     sparsemat%smmm%nseq, nseq_per_line)
      call determine_sequential_length_new2(sparsemat%smmm%nout, ispt, nseg, norb, keyv, keyg, &
           sparsemat, istsegline, sparsemat%smmm%nseq, nseq_per_line)
      !write(*,'(a,i3,3x,200i10)') 'iproc, nseq_per_line', iproc, nseq_per_line
      if (nproc>1) call mpiallred(nseq_per_line(1), norb, mpi_sum, comm=bigdft_mpi%mpi_comm)
      rseq=real(sparsemat%smmm%nseq,kind=8) !real to prevent integer overflow
      if (nproc>1) call mpiallred(rseq, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)


      rseq_per_line = f_malloc(norb,id='rseq_per_line')
      do iorb=1,norb
          rseq_per_line(iorb) = real(nseq_per_line(iorb),kind=8)
      end do


      norb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
      isorb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
      ! Assign the columns of the matrix to the processes such that the load
      ! balancing is optimal
      rseq_ideal = rseq/real(nproc,kind=8)
      call redistribute(nproc, norb, rseq_per_line, rseq_ideal, norb_par_ideal)
      isorb_par_ideal(0) = 0
      do jproc=1,nproc-1
          isorb_par_ideal(jproc) = isorb_par_ideal(jproc-1) + norb_par_ideal(jproc-1)
      end do


      ! some checks
      if (sum(norb_par_ideal)/=norb) stop 'sum(norb_par_ideal)/=norb'
      if (isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb) stop 'isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb'

      ! Copy the values
      sparsemat%smmm%nfvctrp=norb_par_ideal(iproc)
      sparsemat%smmm%isfvctr=isorb_par_ideal(iproc)



      ! Get the load balancing
      rseq_max(1) = real(sparsemat%smmm%nseq,kind=8)
      rseq_average(1) = rseq_max(1)/real(nproc,kind=8)
      !if (nproc>1) call mpiallred(nseq_min, 1, mpi_min, bigdft_mpi%mpi_comm)
      !nseq_max = sparsemat%smmm%nseq
      !if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
      !if (nseq_min>0) then
      !    ratio_before = real(nseq_max,kind=8)/real(nseq_min,kind=8)
      !    printable=.true.
      !else
      !    printable=.false.
      !end if


      ! Realculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
      ! the optimized partitioning of the matrix columns.
      call get_nout(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)
      !call get_nout(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), sparsemat%nseg, &
      !     sparsemat%nsegline, sparsemat%istsegline, sparsemat%keyg, sparsemat%smmm%nout)
      !!call determine_sequential_length(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
      !!     nsegline, istsegline, keyg, sparsemat, &
      !!     sparsemat%smmm%nseq, nseq_per_line)
      !!write(*,*) 'OLD: iproc, nseq', iproc, sparsemat%smmm%nseq

      ! Determine ispt
      ispt = get_offset(iproc, nproc, sparsemat%smmm%nout)
      !!nout_par = f_malloc0(0.to.nproc-1,id='ist_par')
      !!nout_par(iproc) = sparsemat%smmm%nout
      !!call mpiallred(nout_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
      !!ispt = 0
      !!do jproc=0,iproc-1
      !!    ispt = ispt + nout_par(jproc)
      !!end do
      !!nseq_per_pt = f_malloc0(sum(nout_par),id='nseq_per_pt')
      !!call determine_sequential_length_new(sparsemat%smmm%nout, ispt, sparsemat%nseg, sparsemat%keyv, sparsemat%keyg, &
      !!     sparsemat, sum(nout_par), sparsemat%smmm%nseq, nseq_per_pt)
      !!call determine_sequential_length_new(sparsemat%smmm%nout, ispt, nseg, keyv, keyg, &
      !!     sparsemat, sum(nout_par), sparsemat%smmm%nseq, nseq_per_pt)
      !write(*,*) 'norb, sparsemat%nfvctr', norb, sparsemat%nfvctr
      call determine_sequential_length_new2(sparsemat%smmm%nout, ispt, nseg, norb, keyv, keyg, &
           sparsemat, istsegline, sparsemat%smmm%nseq, nseq_per_line)
      !write(*,'(a,i3,3x,200i10)') 'iproc, nseq_per_line', iproc, nseq_per_line
      !!call f_free(nout_par)
      !!write(*,*) 'NEW: iproc, nseq', iproc, sparsemat%smmm%nseq

      ! Get the load balancing
      rseq_max(2) = real(sparsemat%smmm%nseq,kind=8)
      rseq_average(2) = rseq_max(2)/real(nproc,kind=8)
      if (nproc>1) call mpiallred(rseq_max, mpi_max, comm=bigdft_mpi%mpi_comm)
      if (nproc>1) call mpiallred(rseq_average, mpi_sum, comm=bigdft_mpi%mpi_comm)
      !nseq_max = sparsemat%smmm%nseq
      !if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
      ! Not necessary to set the printable flag (if nseq_min was zero before it should be zero here as well)
      !!if (nseq_min>0) then
      !!    ratio_after = real(nseq_max,kind=8)/real(nseq_min,kind=8)
      !!    if (.not.printable) stop 'this should not happen (sparsematrix)'
      !!else
      !!    if (printable) stop 'this should not happen (sparsematrix)'
      !!end if
      ratio_before = rseq_max(1)/rseq_average(1)
      ratio_after = rseq_max(2)/rseq_average(2)
      if (iproc==0) then
          !if (printable) then
              call yaml_map('sparse matmul load balancing naive / optimized',(/ratio_before,ratio_after/),fmt='(f4.2)')
          !!else
          !!    call yaml_map('sparse matmul load balancing naive / optimized','printing not meaningful (division by zero)')
          !!end if
      end if
      

      call f_free(nseq_per_line)
      !!call f_free(nseq_per_pt)
      call f_free(rseq_per_line)


      call allocate_sparse_matrix_matrix_multiplication(nproc, norb, nseg, sparsemat%smmm)
      call vcopy(nseg, keyv(1), 1, sparsemat%smmm%keyv(1), 1)
      call vcopy(4*nseg, keyg(1,1,1), 1, sparsemat%smmm%keyg(1,1,1), 1)
      call vcopy(norb, istsegline(1), 1, sparsemat%smmm%istsegline(1), 1)

      ! Calculate some auxiliary variables
      temparr = f_malloc0((/0.to.nproc-1,1.to.2/),id='isfvctr_par')
      temparr(iproc,1) = sparsemat%smmm%isfvctr
      temparr(iproc,2) = sparsemat%smmm%nfvctrp
      if (nproc>1) then
          call mpiallred(temparr,  mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, sparsemat%nseg, sparsemat%nvctr, &
           temparr(0,1), temparr(0,2), sparsemat%istsegline, sparsemat%keyv, &
           sparsemat%smmm%isvctr_mm, sparsemat%smmm%nvctrp_mm, sparsemat%smmm%isvctr_mm_par, sparsemat%smmm%nvctr_mm_par)

      ! Would be better if this were in the wrapper above...
      sparsemat%smmm%line_and_column_mm = f_malloc_ptr((/2,sparsemat%smmm%nvctrp_mm/),id='smmm%line_and_column_mm')

      call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, nseg, keyv(nseg)+(keyg(2,1,nseg)-keyg(1,1,nseg)), &
           temparr(0,1), temparr(0,2), istsegline, keyv, &
           sparsemat%smmm%isvctr, sparsemat%smmm%nvctrp, sparsemat%smmm%isvctr_par, sparsemat%smmm%nvctr_par)
      call f_free(temparr)

      ! Would be better if this were in the wrapper above...
      sparsemat%smmm%line_and_column = f_malloc_ptr((/2,sparsemat%smmm%nvctrp/),id='smmm%line_and_column')

      ! Init line_and_column
      !!call init_line_and_column()
      call init_line_and_column(sparsemat%smmm%nvctrp_mm, sparsemat%smmm%isvctr_mm, &
           sparsemat%nseg, sparsemat%keyv, sparsemat%keyg, &
           sparsemat%smmm%line_and_column_mm)
      call init_line_and_column(sparsemat%smmm%nvctrp, sparsemat%smmm%isvctr, &
           nseg, keyv, keyg, sparsemat%smmm%line_and_column)
      !!iseg_start = 1
      !!do i=1,sparsemat%smmm%nvctrp_mm
      !!    ii = sparsemat%smmm%isvctr_mm + i
      !!    call get_line_and_column(ii, sparsemat%nseg, sparsemat%keyv, sparsemat%keyg, iseg_start, iline, icolumn)
      !!    sparsemat%smmm%line_and_column_mm(1,i) = iline
      !!    sparsemat%smmm%line_and_column_mm(2,i) = icolumn
      !!end do
      !!iseg_start = 1
      !!do i=1,sparsemat%smmm%nvctrp
      !!    ii = sparsemat%smmm%isvctr + i
      !!    call get_line_and_column(ii, nseg, keyv, keyg, iseg_start, iline, icolumn)
      !!    sparsemat%smmm%line_and_column(1,i) = iline
      !!    sparsemat%smmm%line_and_column(2,i) = icolumn
      !!end do



      ! Get the segments containing the first and last element of a sparse
      ! matrix after a multiplication
      do i=1,2
          if (i==1) then
              iel = sparsemat%smmm%isvctr_mm + 1
          else if (i==2) then
              iel = sparsemat%smmm%isvctr_mm + sparsemat%smmm%nvctrp_mm
          end if
          iiseg = sparsemat%nseg !in case iel is the last element
          do iseg=1,sparsemat%nseg
              ist_seg = sparsemat%keyv(iseg)
              ilen_seg = sparsemat%keyg(2,1,iseg) - sparsemat%keyg(1,1,iseg)
              iend_seg = ist_seg + ilen_seg
              if (iend_seg<iel) cycle
              ! If this point is reached, we are in the correct segment
              iiseg = iseg ; exit
          end do
          if (i==1) then
              sparsemat%smmm%isseg = iiseg
          else if (i==2) then
              sparsemat%smmm%ieseg = iiseg
          end if
      end do

      sparsemat%smmm%nseg=nseg
      call vcopy(norb, nsegline(1), 1, sparsemat%smmm%nsegline(1), 1)
      call vcopy(norb, istsegline(1), 1, sparsemat%smmm%istsegline(1), 1)
      !call init_onedimindices_new(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
      !     nsegline, istsegline, keyg, &
      !     sparsemat, sparsemat%smmm%nout, sparsemat%smmm%onedimindices)
      call init_onedimindices_newnew(sparsemat%smmm%nout, ispt, nseg, &
           keyv, keyg, sparsemat, istsegline, sparsemat%smmm%onedimindices_new)
      !call get_arrays_for_sequential_acces(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
      !     nsegline, istsegline, keyg, sparsemat, &
      !     sparsemat%smmm%nseq, sparsemat%smmm%ivectorindex)
      call get_arrays_for_sequential_acces_new(sparsemat%smmm%nout, ispt, nseg, sparsemat%smmm%nseq, &
           keyv, keyg, sparsemat, istsegline, sparsemat%smmm%ivectorindex_new)
      call determine_consecutive_values(sparsemat%smmm%nout, sparsemat%smmm%nseq, sparsemat%smmm%ivectorindex_new, &
           sparsemat%smmm%onedimindices_new, sparsemat%smmm%nconsecutive_max, sparsemat%smmm%consecutive_lookup)
      !call init_sequential_acces_matrix(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), sparsemat%nseg, &
      !     sparsemat%nsegline, sparsemat%istsegline, sparsemat%keyg, sparsemat, sparsemat%smmm%nseq, &
      !     sparsemat%smmm%indices_extract_sequential)
      call init_sequential_acces_matrix_new(sparsemat%smmm%nout, ispt, nseg, sparsemat%smmm%nseq, keyv, keyg, sparsemat, &
           istsegline, sparsemat%smmm%indices_extract_sequential)

      ! This array gives the starting and ending indices of the submatrix which
      ! is used by a given MPI task
      if (sparsemat%smmm%nseq>0) then
          sparsemat%smmm%istartend_mm(1) = sparsemat%nvctr
          sparsemat%smmm%istartend_mm(2) = 1
          do iseq=1,sparsemat%smmm%nseq
              ind=sparsemat%smmm%indices_extract_sequential(iseq)
              sparsemat%smmm%istartend_mm(1) = min(sparsemat%smmm%istartend_mm(1),ind)
              sparsemat%smmm%istartend_mm(2) = max(sparsemat%smmm%istartend_mm(2),ind)
          end do
      else
          sparsemat%smmm%istartend_mm(1)=sparsemat%nvctr+1
          sparsemat%smmm%istartend_mm(2)=sparsemat%nvctr
      end if

      ! Determine to which segments this corresponds
      sparsemat%smmm%istartendseg_mm(1)=sparsemat%nseg+1
      do iseg=1,sparsemat%nseg
          if (sparsemat%keyv(iseg)+sparsemat%keyg(2,1,iseg)-sparsemat%keyg(1,1,iseg)+1>=sparsemat%smmm%istartend_mm(1)) then
              sparsemat%smmm%istartendseg_mm(1)=iseg
              exit
          end if
      end do
      sparsemat%smmm%istartendseg_mm(2)=0
      do iseg=sparsemat%nseg,1,-1
          if (sparsemat%keyv(iseg)<=sparsemat%smmm%istartend_mm(2)) then
              sparsemat%smmm%istartendseg_mm(2)=iseg
              exit
          end if
      end do

      istartend_mm = f_malloc0((/1.to.2,0.to.nproc-1/),id='istartend_mm')
      istartend_mm(1:2,iproc) = sparsemat%smmm%istartend_mm(1:2)
      if (nproc>1) then
          call mpiallred(istartend_mm, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! Partition the entire matrix in disjoint submatrices
      istartend_dj = f_malloc((/1.to.2,0.to.nproc-1/),id='istartend_dj')
      istartend_dj(1,0) = istartend_mm(1,0)
      do jproc=1,nproc-1
          ind = (istartend_mm(2,jproc-1)+istartend_mm(1,jproc))/2
          ! check that this is inside the segment of istartend_mm(:,jproc)
          ind = max(ind,istartend_mm(1,jproc))
          ind = min(ind,istartend_mm(2,jproc))
          ! check that this is not smaller than the beginning of the previous chunk
          ind = max(ind,istartend_dj(1,jproc-1))+1
          ! check that this is not outside the total matrix size (may happen if there are more processes than matrix columns)
          ind = min(ind,sparsemat%nvctr)
          istartend_dj(1,jproc) = ind
          istartend_dj(2,jproc-1) = istartend_dj(1,jproc)-1
      end do
      istartend_dj(2,nproc-1) = istartend_mm(2,nproc-1)
      !if (iproc==0) write(*,'(a,100(2i7,3x))') 'istartend_mm',istartend_mm
      !if (iproc==0) write(*,'(a,100(2i7,3x))') 'istartend_dj',istartend_dj

      ! Some checks
      if (istartend_dj(1,0)/=1) stop 'istartend_dj(1,0)/=1'
      if (istartend_dj(2,nproc-1)/=sparsemat%nvctr) stop 'istartend_dj(2,nproc-1)/=sparsemat%nvctr'
      ii = 0
      do jproc=0,nproc-1
          ncount = istartend_dj(2,jproc)-istartend_dj(1,jproc) + 1
          if (ncount<0) stop 'ncount<0'
          ii = ii + ncount
          if (ii<0) stop 'init_sparse_matrix_matrix_multiplication: ii<0'
          if (jproc>0) then
              if (istartend_dj(1,jproc)/=istartend_dj(2,jproc-1)+1) stop 'istartend_dj(1,jproc)/=istartend_dj(2,jproc-1)'
          end if
      end do
      if (ii/=sparsemat%nvctr) stop 'init_sparse_matrix_matrix_multiplication: ii/=sparsemat%nvctr'

      ! Keep the values of its own task
      sparsemat%smmm%istartend_mm_dj(1) = istartend_dj(1,iproc)
      sparsemat%smmm%istartend_mm_dj(2) = istartend_dj(2,iproc)

      ! Update the segments...
      !write(*,*) 'sparsemat%smmm%istartend_mm_dj(1)',sparsemat%smmm%istartend_mm_dj(1)
      ii=sparsemat%nseg+1
      do iseg=1,sparsemat%nseg
      !write(*,*) 'sparsemat%smmm%istartend_mm_dj(1)',sparsemat%keyv(iseg), sparsemat%smmm%istartend_mm_dj(1)
          if (sparsemat%keyv(iseg)+sparsemat%keyg(2,1,iseg)-sparsemat%keyg(1,1,iseg)+1>=sparsemat%smmm%istartend_mm_dj(1)) then
              ii=iseg
              exit
          end if
      end do
      if (ii<sparsemat%smmm%istartendseg_mm(1)) sparsemat%smmm%istartendseg_mm(1)=ii
      ii=0
      do iseg=sparsemat%nseg,1,-1
          if (sparsemat%keyv(iseg)<=sparsemat%smmm%istartend_mm_dj(2)) then
              ii=iseg
              exit
          end if
      end do
      if (ii>sparsemat%smmm%istartendseg_mm(2)) sparsemat%smmm%istartendseg_mm(2)=ii


      call f_free(norb_par_ideal)
      call f_free(isorb_par_ideal)
      call f_free(istartend_mm)
      call f_free(istartend_dj)

      call f_release_routine()

    end subroutine init_sparse_matrix_matrix_multiplication_new

    subroutine init_line_and_column(nvctrp, isvctr, nseg, keyv, keyg, line_and_column)
      implicit none

      ! Calling arguments
      integer,intent(in) :: nvctrp, isvctr, nseg
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      integer,dimension(2,nvctrp),intent(out) :: line_and_column

      ! Local variables
      integer :: iseg_start, i, ii, iline, icolumn

      call f_routine(id='init_line_and_column')

      iseg_start = 1
      !$omp parallel default(none) &
      !$omp shared(nvctrp, isvctr, nseg, keyv, keyg, line_and_column) &
      !$omp private(i, ii, iline, icolumn) &
      !$omp firstprivate(iseg_start)
      !$omp do schedule(static)
      do i=1,nvctrp
          ii = isvctr + i
          call get_line_and_column(ii, nseg, keyv, keyg, iseg_start, iline, icolumn)
          line_and_column(1,i) = iline
          line_and_column(2,i) = icolumn
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine init_line_and_column


    !> Calculates the offset of a parallel distribution for each MPI task
    function get_offset(iproc, nproc, n) result(is)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc !< task ID
      integer,intent(in) :: nproc !< total number of tasks
      integer,intent(in) :: n !< size of the distributed quantity on each task
      integer :: is

      ! Local variables
      integer :: jproc
      integer,dimension(1) :: n_, is_
      integer,dimension(:),allocatable :: narr, isarr

      ! Since the wrapper wants arrays
      n_(1) = n
      ! Gather the data on the last process
      narr = f_malloc(0.to.nproc-1,id='narr')
      isarr = f_malloc(0.to.nproc-1,id='n=isarr')
      if (nproc>1) then
          call mpigather(n_, narr, nproc-1)
      else
          narr(0) = n_(1)
      end if
      if (iproc==nproc-1) then
          isarr(0) = 0
          do jproc=1,nproc-1
              isarr(jproc) = isarr(jproc-1) + narr(jproc-1)
          end do
      end if
      if (nproc>1) then
          call mpiscatter(isarr, is_, nproc-1)
      else
          is_(1) = isarr(0)
      end if
      is = is_(1)
      call f_free(narr)
      call f_free(isarr)
    end function get_offset


    subroutine nseg_perline(norb, lut, nseg, nvctr, nsegline)
      implicit none

      ! Calling arguments
      integer,intent(in) :: norb
      logical,dimension(norb),intent(in) :: lut
      integer,intent(inout) :: nseg, nvctr
      integer,intent(out) :: nsegline

      ! Local variables
      integer :: jorb
      logical :: segment_started, newline, overlap

      ! Always start a new segment for each line
      segment_started=.false.
      nsegline=0
      newline=.true.
      do jorb=1,norb
          overlap=lut(jorb)
          if (overlap) then
              if (segment_started) then
                  ! there is no "hole" in between, i.e. we are in the same segment
                  nvctr=nvctr+1
              else
                  ! there was a "hole" in between, i.e. we are in a new segment
                  nseg=nseg+1
                  nsegline=nsegline+1
                  nvctr=nvctr+1
                  newline=.false.
              end if
              segment_started=.true.
          else
              segment_started=.false.
          end if
      end do

    end subroutine nseg_perline


    subroutine keyg_per_line(norb, nseg, iline, istseg, lut, ivctr, keyg)
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: norb, nseg, iline, istseg
      logical,dimension(norb),intent(in) :: lut
      integer,intent(inout) :: ivctr
      integer,dimension(2,2,nseg),intent(out) :: keyg
      
      ! Local variables
      integer :: iseg, jorb, ijorb
      logical :: segment_started, overlap

      ! Always start a new segment for each line
      segment_started=.false.
      !iseg=sparsemat%istsegline(iline)-1
      iseg=istseg-1
      do jorb=1,norb
          overlap=lut(jorb)
          ijorb=(iline-1)*norb+jorb
          if (overlap) then
              if (segment_started) then
                  ! there is no "hole" in between, i.e. we are in the same segment
                  ivctr=ivctr+1
              else
                  ! there was a "hole" in between, i.e. we are in a new segment.
                  iseg=iseg+1
                  ivctr=ivctr+1
                  ! open the current segment
                  keyg(1,1,iseg)=jorb
                  keyg(1,2,iseg)=iline
              end if
              segment_started=.true.
          else
              if (segment_started) then
                  ! close the previous segment
                  keyg(2,1,iseg)=jorb-1
                  keyg(2,2,iseg)=iline
              end if
              segment_started=.false.
          end if
      end do
      ! close the last segment on the line if necessary
      if (segment_started) then
          keyg(2,1,iseg)=norb
          keyg(2,2,iseg)=iline
      end if
    end subroutine keyg_per_line


    !!subroutine keyg_per_line_old(norb, nseg, iline, istseg, lut, ivctr, keyg)
    !!  implicit none
    !!  
    !!  ! Calling arguments
    !!  integer,intent(in) :: norb, nseg, iline, istseg
    !!  logical,dimension(norb),intent(in) :: lut
    !!  integer,intent(inout) :: ivctr
    !!  integer,dimension(2,nseg),intent(out) :: keyg
    !!  
    !!  ! Local variables
    !!  integer :: iseg, jorb, ijorb
    !!  logical :: segment_started, overlap

    !!  ! Always start a new segment for each line
    !!  segment_started=.false.
    !!  !iseg=sparsemat%istsegline(iline)-1
    !!  iseg=istseg-1
    !!  do jorb=1,norb
    !!      overlap=lut(jorb)
    !!      ijorb=(iline-1)*norb+jorb
    !!      if (overlap) then
    !!          if (segment_started) then
    !!              ! there is no "hole" in between, i.e. we are in the same segment
    !!              ivctr=ivctr+1
    !!          else
    !!              ! there was a "hole" in between, i.e. we are in a new segment.
    !!              iseg=iseg+1
    !!              ivctr=ivctr+1
    !!              ! open the current segment
    !!              keyg(1,iseg)=ijorb
    !!          end if
    !!          segment_started=.true.
    !!      else
    !!          if (segment_started) then
    !!              ! close the previous segment
    !!              keyg(2,iseg)=ijorb-1
    !!          end if
    !!          segment_started=.false.
    !!      end if
    !!  end do
    !!  ! close the last segment on the line if necessary
    !!  if (segment_started) then
    !!      keyg(2,iseg)=iline*norb
    !!  end if
    !!end subroutine keyg_per_line_old



    !> Currently assuming square matrices
    subroutine init_sparse_matrix(iproc, nproc, nspin, geocode, norbu, norbup, isorbu, store_index, &
               on_which_atom, nnonzero, nonzero, nnonzero_mult, nonzero_mult, sparsemat, &
               allocate_full, print_info)
      use yaml_output
!      use yaml_strings, only: yaml_toa
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin, norbu, norbup, isorbu, nnonzero, nnonzero_mult
      character(len=1),intent(in) :: geocode
      logical,intent(in) :: store_index
      integer,dimension(norbu),intent(in) :: on_which_atom
      integer,dimension(2,nnonzero),intent(in) :: nonzero
      integer,dimension(2,nnonzero_mult),intent(in) :: nonzero_mult
      type(sparse_matrix), intent(out) :: sparsemat
      logical,intent(in),optional :: allocate_full, print_info
      
      ! Local variables
      integer :: jproc, iorb, jorb, iiorb, iseg
      !integer :: jst_line, jst_seg, segn, ind
      integer :: ist, ivctr
      logical,dimension(:),allocatable :: lut
      integer :: nseg_mult, nvctr_mult, ivctr_mult
      integer,dimension(:),allocatable :: nsegline_mult, istsegline_mult
      integer,dimension(:,:,:),allocatable :: keyg_mult
      integer,dimension(:),allocatable :: keyv_mult
      logical :: allocate_full_, print_info_ !LG: internal variables have the underscore, not the opposite
      integer(kind=8) :: ntot

      real(kind=4) :: tr0, tr1, trt0, trt1
      real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
      logical, parameter :: extra_timing=.false.

      call timing(iproc,'init_matrCompr','ON')
      if (extra_timing) call cpu_time(trt0)
      call f_routine(id='init_sparse_matrix')
      
      allocate_full_=.false.
      print_info_=.true.
      if (present(allocate_full)) allocate_full_=allocate_full
      if (present(print_info)) print_info_=print_info

      lut = f_malloc(norbu,id='lut')
    
      sparsemat=sparse_matrix_null()
    
      sparsemat%nspin=nspin
      sparsemat%geocode = geocode
      sparsemat%nfvctr=norbu
      sparsemat%nfvctrp=norbup
      sparsemat%isfvctr=isorbu
      sparsemat%nfvctr_par=f_malloc0_ptr((/0.to.nproc-1/),id='sparsemat%nfvctr_par')
      sparsemat%isfvctr_par=f_malloc0_ptr((/0.to.nproc-1/),id='sparsemat%isfvctr_par')
      if (extra_timing) call cpu_time(tr0)
      ! Same as isorb_par and norb_par
      do jproc=0,nproc-1
          if (iproc==jproc) then
              sparsemat%isfvctr_par(jproc)=isorbu
              sparsemat%nfvctr_par(jproc)=norbup
          end if
      end do
      if (nproc>1) then
          call mpiallred(sparsemat%isfvctr_par, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nfvctr_par, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      call allocate_sparse_matrix_basic(store_index, norbu, nproc, sparsemat)

      call vcopy(norbu, on_which_atom(1), 1, sparsemat%on_which_atom(1), 1)


      sparsemat%nseg=0
      sparsemat%nvctr=0
      sparsemat%nsegline=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero, nonzero, iiorb, norbu, lut)
          call nseg_perline(norbu, lut, sparsemat%nseg, sparsemat%nvctr, sparsemat%nsegline(iiorb))
      end do

      if (nproc>1) then
          call mpiallred(sparsemat%nvctr, 1, mpi_sum,comm=bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nseg, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nsegline(1), sparsemat%nfvctr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time0=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0)

      ist=1
      do jorb=1,sparsemat%nfvctr
          ! Starting segment for this line
          sparsemat%istsegline(jorb)=ist
          ist=ist+sparsemat%nsegline(jorb)
      end do

    
      if (iproc==0 .and. print_info_) then
          ntot = int(norbu,kind=8)*int(norbu,kind=8)
          call yaml_map('total elements',ntot)
          call yaml_map('non-zero elements',sparsemat%nvctr)
          call yaml_comment('segments: '//sparsemat%nseg)
          call yaml_map('sparsity in %',1.d2*real(ntot-int(sparsemat%nvctr,kind=8),kind=8)/real(ntot,kind=8),fmt='(f5.2)')
      end if
    
      call allocate_sparse_matrix_keys(store_index, sparsemat)
    

      ivctr=0
      sparsemat%keyg=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero, nonzero, iiorb, norbu, lut)
          call keyg_per_line(norbu, sparsemat%nseg, iiorb, sparsemat%istsegline(iiorb), &
               lut, ivctr, sparsemat%keyg)
      end do
    
      ! check whether the number of elements agrees
      if (nproc>1) then
          call mpiallred(ivctr, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (ivctr/=sparsemat%nvctr) then
          write(*,'(a,2i8)') 'ERROR: ivctr/=sparsemat%nvctr', ivctr, sparsemat%nvctr
          stop
      end if
      if (nproc>1) then
          call mpiallred(sparsemat%keyg(1,1,1), 2*2*sparsemat%nseg, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! start of the segments
      sparsemat%keyv(1)=1
      do iseg=2,sparsemat%nseg
          ! A segment is always on one line, therefore no double loop
          sparsemat%keyv(iseg) = sparsemat%keyv(iseg-1) + sparsemat%keyg(2,1,iseg-1) - sparsemat%keyg(1,1,iseg-1) + 1
      end do


      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time1=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0)   

      if (store_index) then
          ! store the indices of the matrices in the sparse format
          sparsemat%store_index=.true.
    
          ! initialize sparsemat%matrixindex_in_compressed
          call f_zero(sparsemat%matrixindex_in_compressed_arr)

          !$omp parallel do default(private) shared(sparsemat,norbu,norbup,isorbu) 
          do iorb=1,norbu
             do jorb=1,norbup
                !sparsemat%matrixindex_in_compressed_arr(iorb,jorb)=compressed_index(iorb,jorb,norbu,sparsemat)
                sparsemat%matrixindex_in_compressed_arr(iorb,jorb+isorbu) = &
                     matrixindex_in_compressed(sparsemat, iorb, jorb+isorbu, .true., norbu)
             end do
          end do
          !$omp end parallel do

          if (nproc>1) then
              call mpiallred(sparsemat%matrixindex_in_compressed_arr(1,1), norbu*norbu, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if

          !!! Initialize sparsemat%orb_from_index
          !!ind = 0
          !!do iseg = 1, sparsemat%nseg
          !!   do segn = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
          !!      ind=ind+1
          !!      iorb = (segn - 1) / sparsemat%nfvctr + 1
          !!      jorb = segn - (iorb-1)*sparsemat%nfvctr
          !!      sparsemat%orb_from_index(1,ind) = jorb
          !!      sparsemat%orb_from_index(2,ind) = iorb
          !!   end do
          !!end do
    
      else
          ! Otherwise always calculate them on-the-fly
          sparsemat%store_index=.false.
      end if
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time2=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0)     


      ! parallelization of matrices, following same idea as norb/norbp/isorb
      !most equal distribution, but want corresponding to norbp for second column
      call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, sparsemat%nseg, sparsemat%nvctr, &
           sparsemat%isfvctr_par, sparsemat%nfvctr_par, sparsemat%istsegline, sparsemat%keyv, &
           sparsemat%isvctr, sparsemat%nvctrp, sparsemat%isvctr_par, sparsemat%nvctr_par)
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time3=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0) 

      !!do jproc=0,nproc-1
      !!    jst_line = sparsemat%isfvctr_par(jproc)+1
      !!    if (sparsemat%nfvctr_par(jproc)==0) then
      !!       sparsemat%isvctr_par(jproc) = sparsemat%nvctr
      !!    else
      !!       jst_seg = sparsemat%istsegline(jst_line)
      !!       sparsemat%isvctr_par(jproc) = sparsemat%keyv(jst_seg)-1
      !!    end if
      !!end do
      !!do jproc=0,nproc-1
      !!   if (jproc==nproc-1) then
      !!      sparsemat%nvctr_par(jproc)=sparsemat%nvctr-sparsemat%isvctr_par(jproc)
      !!   else
      !!      sparsemat%nvctr_par(jproc)=sparsemat%isvctr_par(jproc+1)-sparsemat%isvctr_par(jproc)
      !!   end if
      !!   if (iproc==jproc) sparsemat%isvctr=sparsemat%isvctr_par(jproc)
      !!   if (iproc==jproc) sparsemat%nvctrp=sparsemat%nvctr_par(jproc)
      !!end do

    
      ! 0 - none, 1 - mpiallred, 2 - allgather
      sparsemat%parallel_compression=0


      nsegline_mult = f_malloc0(norbu,id='nsegline_mult')
      istsegline_mult = f_malloc(norbu,id='istsegline_mult')
      nseg_mult=0
      nvctr_mult=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero_mult, nonzero_mult, iiorb, norbu, lut)
          call nseg_perline(norbu, lut, nseg_mult, nvctr_mult, nsegline_mult(iiorb))
      end do
      if (nproc>1) then
          call mpiallred(nvctr_mult, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nseg_mult, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nsegline_mult, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if



      ! Initialize istsegline, which gives the first segment of each line
      istsegline_mult(1)=1
      do iorb=2,norbu
          istsegline_mult(iorb) = istsegline_mult(iorb-1) + nsegline_mult(iorb-1)
      end do

      keyg_mult = f_malloc0((/2,2,nseg_mult/),id='keyg_mult')
      keyv_mult = f_malloc0((/nseg_mult/),id='keyg_mult')

      ivctr_mult=0
      do iorb=1,norbup
         iiorb=isorbu+iorb
         call create_lookup_table(nnonzero_mult, nonzero_mult, iiorb, norbu, lut)
         call keyg_per_line(norbu, nseg_mult, iiorb, istsegline_mult(iiorb), &
              lut, ivctr_mult, keyg_mult)
      end do
      ! check whether the number of elements agrees
      if (nproc>1) then
          call mpiallred(ivctr_mult, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (ivctr_mult/=nvctr_mult) then
          write(*,'(a,2i8)') 'ERROR: ivctr_mult/=nvctr_mult', ivctr_mult, nvctr_mult
          stop
      end if
      if (nproc>1) then
          call mpiallred(keyg_mult, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! start of the segments
      keyv_mult(1)=1
      do iseg=2,nseg_mult
          keyv_mult(iseg) = keyv_mult(iseg-1) + keyg_mult(2,1,iseg-1) - keyg_mult(1,1,iseg-1) + 1
      end do
      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time4=real(tr1-tr0,kind=8)

      if (extra_timing) call cpu_time(tr0) 

      ! Allocate the matrices
      !call allocate_sparse_matrix_matrices(sparsemat, allocate_full_)


      ! Initialize the parameters for the spare matrix matrix multiplication
      call init_sparse_matrix_matrix_multiplication_new(iproc, nproc, norbu, norbup, isorbu, nseg_mult, &
               nsegline_mult, istsegline_mult, keyv_mult, keyg_mult, sparsemat)

      if (extra_timing) call cpu_time(tr1)   
      if (extra_timing) time5=real(tr1-tr0,kind=8)    

      call f_free(nsegline_mult)
      call f_free(istsegline_mult)
      call f_free(keyg_mult)
      call f_free(keyv_mult)
      call f_free(lut)
    
      call f_release_routine()
      call timing(iproc,'init_matrCompr','OF')
      if (extra_timing) call cpu_time(trt1)   
      if (extra_timing) ttime=real(trt1-trt0,kind=8)

      if (extra_timing.and.iproc==0) print*,'imctime',time0,time1,time2,time3,time4,time5,&
           time0+time1+time2+time3+time4+time5,ttime

    end subroutine init_sparse_matrix

    subroutine create_lookup_table(nnonzero, nonzero, iiorb, norbu, lut)
      implicit none
      ! Calling arguments
      integer :: nnonzero, iiorb,norbu
      integer,dimension(2,nnonzero) :: nonzero
      logical, dimension(norbu), intent(inout) :: lut

      ! Local variables
      integer(kind=8) :: ist, iend, ind
      integer :: i, jjorb

      lut = .false.
      ist = int(iiorb-1,kind=8)*int(norbu,kind=8) + int(1,kind=8)
      iend = int(iiorb,kind=8)*int(norbu,kind=8)
      do i=1,nnonzero
         ind = int(nonzero(2,i)-1,kind=8)*int(norbu,kind=8) + int(nonzero(1,i),kind=8)
         if (ind<ist) cycle
         if (ind>iend) exit
         jjorb=nonzero(1,i)
         lut(jjorb)=.true.
      end do
    end subroutine create_lookup_table

    subroutine determine_sequential_length(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, &
               sparsemat, nseq, nseq_per_line)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: sparsemat
      integer,intent(out) :: nseq
      integer,dimension(norb),intent(out) :: nseq_per_line
    
      ! Local variables
      integer :: i,iseg,jorb,iorb,jseg,ii,nseqline
      integer :: isegoffset, istart, iend
    
      nseq=0
      do i = 1,norbp
         ii=isorb+i
         isegoffset=istsegline(ii)-1
         nseqline=0
         do iseg=1,nsegline(ii)
              ! A segment is always on one line, therefore no double loop
              istart=keyg(1,1,isegoffset+iseg)
              iend=keyg(2,1,isegoffset+iseg)
              do iorb=istart,iend
                  !write(*,*) 'old: iproc, iorb, ii', bigdft_mpi%iproc, iorb, ii
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      !write(*,*) 'old: iproc, jseg', bigdft_mpi%iproc, jseg
                      ! A segment is always on one line, therefore no double loop
                      do jorb = sparsemat%keyg(1,1,jseg),sparsemat%keyg(2,1,jseg)
                          nseq=nseq+1
                          nseqline=nseqline+1
                      end do
                  end do
              end do
         end do
         nseq_per_line(ii)=nseqline
      end do 
    
    end subroutine determine_sequential_length
    


    !!subroutine determine_sequential_length_new(npt, ispt, nseg, keyv, keyg, smat, nsize_npp, nseq, nseq_per_pt)
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  integer,intent(in) :: npt, ispt, nseg, nsize_npp
    !!  integer,dimension(nseg),intent(in) :: keyv
    !!  integer,dimension(2,2,nseg),intent(in) :: keyg
    !!  type(sparse_matrix),intent(in) :: smat
    !!  integer,intent(out) :: nseq
    !!  integer,dimension(nsize_npp),intent(out) :: nseq_per_pt
    !!
    !!  ! Local variables
    !!  integer :: ipt, iipt, iline, icolumn, nseq_pt, jseg, jorb, iseg_start


    !!  nseq = 0
    !!  iseg_start = 1
    !!  do ipt=1,npt
    !!      iipt = ispt + ipt
    !!      call get_line_and_column(iipt, nseg, keyv, keyg, iseg_start, iline, icolumn)
    !!      !write(*,'(a,4i8)') 'ipt, iipt, iline, icolumn', ipt, iipt, iline, icolumn
    !!      nseq_pt = 0
    !!      ! Take the column due to the symmetry of the sparsity pattern
    !!      !write(*,*) 'new: iproc, iorb, ii', bigdft_mpi%iproc, icolumn, iline
    !!      do jseg=smat%istsegline(icolumn),smat%istsegline(icolumn)+smat%nsegline(icolumn)-1
    !!          !write(*,*) 'new: iproc, jseg', bigdft_mpi%iproc, jseg
    !!          ! A segment is always on one line, therefore no double loop
    !!          do jorb = smat%keyg(1,1,jseg),smat%keyg(2,1,jseg)
    !!              nseq = nseq + 1
    !!              nseq_pt = nseq_pt + 1
    !!          end do
    !!      end do
    !!      nseq_per_pt(iipt) = nseq_pt
    !!  end do

    !!
    !!  !!nseq=0
    !!  !!do i = 1,norbp
    !!  !!   ii=isorb+i
    !!  !!   isegoffset=istsegline(ii)-1
    !!  !!   nseqline=0
    !!  !!   do iseg=1,nsegline(ii)
    !!  !!        ! A segment is always on one line, therefore no double loop
    !!  !!        istart=keyg(1,1,isegoffset+iseg)
    !!  !!        iend=keyg(2,1,isegoffset+iseg)
    !!  !!        do iorb=istart,iend
    !!  !!            do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
    !!  !!                ! A segment is always on one line, therefore no double loop
    !!  !!                do jorb = sparsemat%keyg(1,1,jseg),sparsemat%keyg(2,1,jseg)
    !!  !!                    nseq=nseq+1
    !!  !!                    nseqline=nseqline+1
    !!  !!                end do
    !!  !!            end do
    !!  !!        end do
    !!  !!   end do
    !!  !!   nseq_per_line(ii)=nseqline
    !!  !!end do 

    !!
    !!end subroutine determine_sequential_length_new



    subroutine determine_sequential_length_new2(npt, ispt, nseg, nline, keyv, keyg, smat, istsegline, nseq, nseq_per_line)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: npt, ispt, nseg, nline
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(smat%nfvctr),intent(in) :: istsegline
      integer,intent(out) :: nseq
      integer,dimension(nline),intent(out) :: nseq_per_line
    
      ! Local variables
      integer :: ipt, iipt, iline, icolumn, nseq_pt, jseg, jorb, ii, iseg_start

      call f_routine(id='determine_sequential_length_new2')

      call f_zero(nseq_per_line)

      ! In the following OMP loop, do a reduction of nseq_per_line to avoid the
      ! need of putting a critical statement around its update.

      nseq = 0
      iseg_start = 1
      !$omp parallel default(none) &
      !$omp shared(npt, ispt, nseg, keyv, keyg, smat, nline, istsegline, nseq, nseq_per_line) &
      !$omp private(ipt, iipt, iline, icolumn, jseg, jorb, ii) &
      !$omp firstprivate(iseg_start)
      !$omp do reduction(+:nseq,nseq_per_line)
      do ipt=1,npt
          iipt = ispt + ipt
          call get_line_and_column(iipt, nseg, keyv, keyg, iseg_start, iline, icolumn)
          ! Take the column due to the symmetry of the sparsity pattern
          do jseg=smat%istsegline(icolumn),smat%istsegline(icolumn)+smat%nsegline(icolumn)-1
              ! A segment is always on one line, therefore no double loop
              do jorb = smat%keyg(1,1,jseg),smat%keyg(2,1,jseg)
                  ! Calculate the index in the large compressed format
                  ii = matrixindex_in_compressed_lowlevel(jorb, iline, nline, nseg, keyv, keyg, istsegline)
                  if (ii>0) then
                      nseq = nseq + 1
                      nseq_per_line(iline) = nseq_per_line(iline) + 1
                  end if
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()
    
    end subroutine determine_sequential_length_new2



    !> Determines the line and column indices on an elements iel for a sparsity
    !! pattern defined by nseg, kev, keyg.
    !! iseg_start is the segment where the search starts and can thus be used to
    !! accelerate the loop (useful if this routine is called several times with
    !! steadily increasing values of iel).
    subroutine get_line_and_column(iel, nseg, keyv, keyg, iseg_start, iline, icolumn)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iel, nseg
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      integer,intent(inout) :: iseg_start
      integer,intent(out) :: iline, icolumn

      ! Local variables
      integer :: iseg, ilen_seg, ist_seg, iend_seg, i, ii
      logical :: found

      found = .false.
      ! Search the segment which contains iel
      search_loop: do iseg=iseg_start,nseg
          ilen_seg = keyg(2,1,iseg) - keyg(1,1,iseg) + 1
          ist_seg = keyv(iseg)
          iend_seg = ist_seg + ilen_seg - 1
          !write(1000+bigdft_mpi%iproc,*) 'iel, iseg, iend_seg', iel, iseg, iend_seg
          if (iend_seg<iel) cycle
          ! If this point is reached, we are in the correct segment
          iline = keyg(1,2,iseg)
          icolumn = keyg(1,1,iseg)
          do i=ist_seg,iend_seg
              !write(1000+bigdft_mpi%iproc,*) 'iline, icolumn', iline, icolumn
              if (i==iel) then
                  ii = iseg
                  found = .true.
                  exit search_loop
              end if
              icolumn = icolumn + 1
          end do
      end do search_loop

      if (.not.found) then
          !write(*,*) 'iseg_start, nseg', iseg_start, nseg
          !do iseg=iseg_start,nseg
          !    write(*,'(a,4i8)') 'iseg, keyv, keyg', iseg, keyv(iseg), keyg(1,1,iseg), keyg(2,1,iseg)
          !end do
          call f_err_throw('get_line_and_column failed to determine the indices, iel='//iel, &
              err_id=BIGDFT_RUNTIME_ERROR)
      end if
      
      iseg_start = ii

    end subroutine get_line_and_column



    subroutine get_nout(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, nout)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,2,nseg),intent(in) :: keyg
      integer,intent(out) :: nout
    
      ! Local variables
      integer :: i, iii, iseg, iorb
      integer :: isegoffset, istart, iend

      call f_routine(id='get_nout')
    
      ! OpenMP for a norbp loop is not ideal, but better than nothing.
      nout=0
      !$omp parallel default(none) &
      !$omp shared(norbp, isorb, istsegline, nsegline, keyg, nout) &
      !$omp private(i, iii, isegoffset, iseg, istart, iend, iorb)
      !$omp do reduction(+:nout)
      do i=1,norbp
         iii=isorb+i
         isegoffset=istsegline(iii)-1
         do iseg=1,nsegline(iii)
              ! A segment is always on one line, therefore no double loop
              istart=keyg(1,1,isegoffset+iseg)
              iend=keyg(2,1,isegoffset+iseg)
              do iorb=istart,iend
                  nout=nout+1
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()
    
    end subroutine get_nout


    subroutine init_onedimindices_new(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, sparsemat, nout, onedimindices)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: sparsemat
      integer,intent(in) :: nout
      integer,dimension(4,nout) :: onedimindices
    
      ! Local variables
      integer :: i, iii, iseg, iorb, ii, jseg, ilen, itot
      integer :: isegoffset, istart, iend
    
    
      ii=0
      itot=1
      do i = 1,norbp
         iii=isorb+i
         isegoffset=istsegline(iii)-1
         do iseg=1,nsegline(iii)
              istart=keyg(1,1,isegoffset+iseg)
              iend=keyg(2,1,isegoffset+iseg)
              ! A segment is always on one line, therefore no double loop
              do iorb=istart,iend
                  ii=ii+1
                  onedimindices(1,ii)=i
                  onedimindices(2,ii)=iorb
                  ilen=0
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      ! A segment is always on one line, therefore no double loop
                      ilen=ilen+sparsemat%keyg(2,1,jseg)-sparsemat%keyg(1,1,jseg)+1
                  end do
                  onedimindices(3,ii)=ilen
                  onedimindices(4,ii)=itot
                  itot=itot+ilen
              end do
          end do
      end do
    
    end subroutine init_onedimindices_new



    subroutine init_onedimindices_newnew(nout, ispt, nseg, keyv, keyg, smat, istsegline, onedimindices)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: nout, ispt, nseg
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(smat%nfvctr),intent(in) :: istsegline
      integer,dimension(4,nout) :: onedimindices
    
      ! Local variables
      integer :: itot, ipt, iipt, iline, icolumn, ilen, jseg, ii, jorb, iseg_start
    
      !!write(*,*) 'iproc, nout, ispt', bigdft_mpi%iproc, nout, ispt
      call f_routine(id='init_onedimindices_newnew')

      ! Handle index 3 separately to enable OpenMP
    
      !itot = 1
      iseg_start = 1
      !$omp parallel default(none) &
      !$omp shared(nout, ispt, nseg, keyv, keyg, onedimindices, smat, istsegline) &
      !$omp firstprivate(iseg_start) &
      !$omp private(ipt, iipt, iline, icolumn, ilen, jseg, jorb, ii)
      !$omp do
      do ipt=1,nout
          iipt = ispt + ipt
          call get_line_and_column(iipt, nseg, keyv, keyg, iseg_start, iline, icolumn)
          onedimindices(1,ipt) = matrixindex_in_compressed_lowlevel(icolumn, iline, smat%nfvctr, nseg, keyv, keyg, istsegline)
          if (onedimindices(1,ipt)>0) then
              onedimindices(1,ipt) = onedimindices(1,ipt) - smat%smmm%isvctr
          else
              stop 'onedimindices(1,ipt)==0'
          end if
          ilen = 0
          ! Take the column due to the symmetry of the sparsity pattern
          do jseg=smat%istsegline(icolumn),smat%istsegline(icolumn)+smat%nsegline(icolumn)-1
              ! A segment is always on one line, therefore no double loop
              do jorb = smat%keyg(1,1,jseg),smat%keyg(2,1,jseg)
                  ! Calculate the index in the large compressed format
                  ii = matrixindex_in_compressed_lowlevel(jorb, iline, smat%nfvctr, nseg, keyv, keyg, istsegline)
                  if (ii>0) then
                      ilen = ilen + 1
                  end if
              end do
          end do
          onedimindices(2,ipt) = ilen
          !!onedimindices(3,ipt) = itot
          !itot = itot + ilen
      end do
      !$omp end do
      !$omp end parallel


      itot = 1
      do ipt=1,nout
          onedimindices(3,ipt) = itot
          itot = itot + onedimindices(2,ipt)
      end do



      call f_release_routine()
    
    end subroutine init_onedimindices_newnew




    subroutine get_arrays_for_sequential_acces(norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat, nseq, &
               ivectorindex)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg, nseq
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: sparsemat
      integer,dimension(nseq),intent(out) :: ivectorindex
    
      ! Local variables
      integer :: i,iseg,jorb,jjorb,iorb,jseg,ii,iii
      integer :: isegoffset, istart, iend
    
    
      ii=1
      do i = 1,norbp
         iii=isorb+i
         isegoffset=istsegline(iii)-1
         do iseg=1,nsegline(iii)
              istart=keyg(1,1,isegoffset+iseg)
              iend=keyg(2,1,isegoffset+iseg)
              ! A segment is always on one line, therefore no double loop
              do iorb=istart,iend
                  !!istindexarr(iorb-istart+1,iseg,i)=ii
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      ! A segment is always on one line, therefore no double loop
                      do jorb = sparsemat%keyg(1,1,jseg),sparsemat%keyg(2,1,jseg)
                          jjorb = jorb
                          ivectorindex(ii)=jjorb
                          ii = ii+1
                      end do
                  end do
              end do
         end do
      end do 
      if (ii/=nseq+1) stop 'ii/=nseq+1'
    
    end subroutine get_arrays_for_sequential_acces



    subroutine get_arrays_for_sequential_acces_new(nout, ispt, nseg, nseq, keyv, keyg, smat, istsegline, ivectorindex)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: nout, ispt, nseg, nseq
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(smat%nfvctr),intent(in) :: istsegline
      integer,dimension(nseq),intent(out) :: ivectorindex
    
      ! Local variables
      integer :: ii, ipt, iipt, iline, icolumn, jseg, jorb, itest, ind, iseg_start
      integer :: ithread, jthread, nthread
      integer,dimension(:),allocatable :: iiarr
      integer,dimension(:,:),pointer :: ise
      integer,dimension(:,:),allocatable :: ivectorindex_work
      !$ integer :: omp_get_thread_num

      call f_routine(id='get_arrays_for_sequential_acces_new')


      ! OpenMP parallelization using a large workarray
      !!nthread = 1
      !!!$ nthread = omp_get_max_threads()

      !!! Determine the number of iterations to be done by each thread
      !!n = f_malloc(0.to.nthread-1,id='n')
      !!ii = nout/nthread
      !!n(0:nthread-1) = ii
      !!ii = nout - nthread*ii
      !!n(0:ii-1) = n(0:ii-1) + 1
      !!! Check
      !!if (sum(n)/=nout) call f_err_throw('sum(n)/=nout',err_name='BIGDFT_RUNTIME_ERROR')

      !!! Determine the first and last iteration for each thread
      !!ise = f_malloc((/1.to.2,0.to.nthread-1/),id='ise')
      !!ise(1,0) = 1
      !!do jthread=1,nthread-1
      !!    ise(1,jthread) = ise(1,jthread-1) + n(jthread-1)
      !!    ise(2,jthread-1) = ise(1,jthread) -1
      !!end do
      !!ise(2,nthread-1) = nout
      !!! Check
      !!ii = 0
      !!do jthread=0,nthread-1
      !!    ii = ii + ise(2,jthread) - ise(1,jthread) + 1
      !!    if (jthread>1) then
      !!        if (ise(1,jthread)/=ise(2,jthread-1)+1) then
      !!            call f_err_throw('ise(1,jthread)/=ise(2,jthread-1)+1',err_name='BIGDFT_RUNTIME_ERROR')
      !!        end if
      !!    end if
      !!end do
      !!if (ii/=nout) call f_err_throw('ii/=nout',err_name='BIGDFT_RUNTIME_ERROR')

      call distribute_on_threads(nout, nthread, ise)
    
      iiarr = f_malloc(0.to.nthread-1,id='iiarr')
      ii = 0
      iseg_start = 1
      ithread = 0
      ivectorindex_work = f_malloc((/1.to.nseq,0.to.nthread-1/),id='ivectorindex_work')
      !$omp parallel &
      !$omp default (none) &
      !$omp shared(ise, ispt, nseg, keyv, keyg, smat, istsegline, iiarr, nthread) &
      !$omp shared(ivectorindex_work, ivectorindex, nseq) &
      !$omp private(ipt, iipt, iline, icolumn, ind, jthread,jseg,jorb) &
      !$omp firstprivate(ii, iseg_start, ithread)
      !$ ithread = omp_get_thread_num()
      do ipt=ise(1,ithread),ise(2,ithread)
          iipt = ispt + ipt
          call get_line_and_column(iipt, nseg, keyv, keyg, iseg_start, iline, icolumn)
          ! Take the column due to the symmetry of the sparsity pattern
          do jseg=smat%istsegline(icolumn),smat%istsegline(icolumn)+smat%nsegline(icolumn)-1
              ! A segment is always on one line, therefore no double loop
              do jorb = smat%keyg(1,1,jseg),smat%keyg(2,1,jseg)
                  ind = matrixindex_in_compressed_lowlevel(jorb, iline, smat%nfvctr, nseg, keyv, keyg, istsegline)
                  if (ind>0) then
                      ii = ii+1
                      ivectorindex_work(ii,ithread) = ind - smat%smmm%isvctr
                      if (ivectorindex_work(ii,ithread)<=0) then
                          stop 'ivectorindex_work(ii,ithread)<=0'
                      end if
                  end if
              end do
          end do
      end do
      iiarr(ithread) = ii
      !$omp barrier
      if (sum(iiarr)/=nseq) stop 'sum(iiarr)/=nseq'

      ii = 1
      do jthread=0,nthread-1
          if (ithread==jthread) then
              if (iiarr(jthread)>0) then
                  call f_memcpy(n=iiarr(jthread), src=ivectorindex_work(1,ithread), dest=ivectorindex(ii))
              end if
          end if
          ii = ii + iiarr(jthread)
      end do
      !$omp end parallel

      call f_free(ivectorindex_work)
      !call f_free(n)
      call f_free_ptr(ise)
      call f_free(iiarr)

      call f_release_routine()

    
    end subroutine get_arrays_for_sequential_acces_new



    subroutine determine_consecutive_values(nout, nseq, ivectorindex, onedimindices_new, &
               nconsecutive_max, consecutive_lookup)
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nout, nseq
      integer,dimension(nseq),intent(in) :: ivectorindex
      integer,dimension(4,nout),intent(inout) :: onedimindices_new
      integer,intent(out) :: nconsecutive_max
      integer,dimension(:,:,:),pointer,intent(out) :: consecutive_lookup

      ! Local variables
      integer :: iout, ilen, ii, iend, nconsecutive, jorb, jjorb, jjorb_prev, iconsec


      nconsecutive_max = 0
      do iout=1,nout
          ilen=onedimindices_new(2,iout)
          ii=onedimindices_new(3,iout)

          iend=ii+ilen-1

          nconsecutive = 1
          do jorb=ii,iend
             jjorb=ivectorindex(jorb)
             if (jorb>ii) then
                 if (jjorb/=jjorb_prev+1) then
                     nconsecutive = nconsecutive + 1
                 end if
             end if
             jjorb_prev = jjorb
          end do
          nconsecutive_max = max(nconsecutive,nconsecutive_max)
          onedimindices_new(4,iout) = nconsecutive
      end do

      consecutive_lookup = f_malloc_ptr((/3,nconsecutive_max,nout/),id='consecutive_lookup')


      do iout=1,nout
          ilen=onedimindices_new(2,iout)
          ii=onedimindices_new(3,iout)

          iend=ii+ilen-1

          nconsecutive = 1
          iconsec = 0
          consecutive_lookup(1,nconsecutive,iout) = ii
          consecutive_lookup(2,nconsecutive,iout) = ivectorindex(ii)
          do jorb=ii,iend
             jjorb=ivectorindex(jorb)
             if (jorb>ii) then
                 if (jjorb/=jjorb_prev+1) then
                     consecutive_lookup(3,nconsecutive,iout) = iconsec
                     nconsecutive = nconsecutive + 1
                     consecutive_lookup(1,nconsecutive,iout) = jorb
                     consecutive_lookup(2,nconsecutive,iout) = jjorb
                     iconsec = 0
                 end if
             end if
             iconsec = iconsec + 1
             jjorb_prev = jjorb
          end do
          consecutive_lookup(3,nconsecutive,iout) = iconsec
          if (nconsecutive>nconsecutive_max) stop 'nconsecutive>nconsecutive_max'
      end do


    end subroutine determine_consecutive_values


    subroutine init_sequential_acces_matrix(norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat, nseq, &
               indices_extract_sequential)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg, nseq
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: sparsemat
      integer,dimension(nseq),intent(out) :: indices_extract_sequential
    
      ! Local variables
      integer :: i,iseg,jorb,jj,iorb,jseg,ii,iii
      integer :: isegoffset, istart, iend
    
    
      ii=1
      do i = 1,norbp
         iii=isorb+i
         isegoffset=istsegline(iii)-1
         do iseg=1,nsegline(iii)
              istart=keyg(1,1,isegoffset+iseg)
              iend=keyg(2,1,isegoffset+iseg)
              ! A segment is always on one line, therefore no double loop
              do iorb=istart,iend
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      ! A segment is always on one line, therefore no double loop
                      jj=1
                      do jorb = sparsemat%keyg(1,1,jseg),sparsemat%keyg(2,1,jseg)
                          indices_extract_sequential(ii)=sparsemat%keyv(jseg)+jj-1
                          jj = jj+1
                          ii = ii+1
                      end do
                  end do
              end do
         end do
      end do 
    
    end subroutine init_sequential_acces_matrix




    subroutine init_sequential_acces_matrix_new(nout, ispt, nseg, nseq, keyv, keyg, smat, istsegline, &
               indices_extract_sequential)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: nout, ispt, nseg, nseq
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(smat%nfvctr),intent(in) :: istsegline
      integer,dimension(nseq),intent(out) :: indices_extract_sequential
    
      ! Local variables
      integer :: ii, ipt, iipt, iline, icolumn, jseg, jj, jorb, ind, iseg_start
      integer :: ithread, jthread, nthread
      integer,dimension(:),allocatable :: iiarr
      integer,dimension(:,:),pointer :: ise
      integer,dimension(:,:),allocatable :: indices_extract_sequential_work
      !$ integer :: omp_get_thread_num

      call f_routine(id='init_sequential_acces_matrix_new')

      ! OpenMP parallelization using a large workarray
      !!nthread = 1
      !!!$ nthread = omp_get_max_threads()

      !!! Determine the number of iterations to be done by each thread
      !!n = f_malloc(0.to.nthread-1,id='n')
      !!ii = nout/nthread
      !!n(0:nthread-1) = ii
      !!ii = nout - nthread*ii
      !!n(0:ii-1) = n(0:ii-1) + 1
      !!! Check
      !!if (sum(n)/=nout) call f_err_throw('sum(n)/=nout',err_name='BIGDFT_RUNTIME_ERROR')

      !!! Determine the first and last iteration for each thread
      !!ise = f_malloc((/1.to.2,0.to.nthread-1/),id='ise')
      !!ise(1,0) = 1
      !!do jthread=1,nthread-1
      !!    ise(1,jthread) = ise(1,jthread-1) + n(jthread-1)
      !!    ise(2,jthread-1) = ise(1,jthread) -1
      !!end do
      !!ise(2,nthread-1) = nout
      !!! Check
      !!ii = 0
      !!do jthread=0,nthread-1
      !!    ii = ii + ise(2,jthread) - ise(1,jthread) + 1
      !!    if (jthread>1) then
      !!        if (ise(1,jthread)/=ise(2,jthread-1)+1) then
      !!            call f_err_throw('ise(1,jthread)/=ise(2,jthread-1)+1',err_name='BIGDFT_RUNTIME_ERROR')
      !!        end if
      !!    end if
      !!end do
      !!if (ii/=nout) call f_err_throw('ii/=nout',err_name='BIGDFT_RUNTIME_ERROR')
      call distribute_on_threads(nout, nthread, ise)
    
      iiarr = f_malloc(0.to.nthread-1,id='iiarr')
      ii = 0
      iseg_start = 1
      ithread = 0
      indices_extract_sequential_work = f_malloc((/1.to.nseq,0.to.nthread-1/),id='indices_extract_sequential_work')
      !$omp parallel &
      !$omp default (none) &
      !$omp shared(ise, ispt, nseg, keyv, keyg, smat, istsegline, iiarr, nthread) &
      !$omp shared(indices_extract_sequential_work, indices_extract_sequential) &
      !$omp private(ipt, iipt, iline, icolumn, ind, jj, jthread,jseg,jorb) &
      !$omp firstprivate(ii, iseg_start, ithread)
      !$ ithread = omp_get_thread_num()
      do ipt=ise(1,ithread),ise(2,ithread)
          iipt = ispt + ipt
          call get_line_and_column(iipt, nseg, keyv, keyg, iseg_start, iline, icolumn)
          ! Take the column due to the symmetry of the sparsity pattern
          do jseg=smat%istsegline(icolumn),smat%istsegline(icolumn)+smat%nsegline(icolumn)-1
              ! A segment is always on one line, therefore no double loop
              jj=1
              do jorb = smat%keyg(1,1,jseg),smat%keyg(2,1,jseg)
                  ! Calculate the index in the large compressed format
                  ind = matrixindex_in_compressed_lowlevel(jorb, iline, smat%nfvctr, nseg, keyv, keyg, istsegline)
                  if (ind>0) then
                      ii = ii + 1
                      indices_extract_sequential_work(ii,ithread)=smat%keyv(jseg)+jj-1
                  end if
                  jj = jj+1
              end do
          end do
      end do
      iiarr(ithread) = ii
      !$omp barrier
      ii = 1
      do jthread=0,nthread-1
          if (ithread==jthread) then
              if (iiarr(jthread)>0) then
                  call f_memcpy(n=iiarr(jthread), &
                       src=indices_extract_sequential_work(1,ithread), &
                       dest=indices_extract_sequential(ii))
              end if
          end if
          ii = ii + iiarr(jthread)
      end do
      !$omp end parallel
      !do ii=1,nseq
      !    write(100,*) 'ii, indices_extract_sequential(ii)', ii, indices_extract_sequential(ii)
      !end do

      call f_free(indices_extract_sequential_work)
      !call f_free(n)
      call f_free_ptr(ise)
      call f_free(iiarr)

      call f_release_routine()

    
    end subroutine init_sequential_acces_matrix_new




    subroutine init_matrix_parallelization(iproc, nproc, nfvctr, nseg, nvctr, &
               isfvctr_par, nfvctr_par, istsegline, keyv, &
               isvctr, nvctrp, isvctr_par, nvctr_par)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nfvctr, nseg, nvctr
      integer,dimension(0:nproc-1),intent(in) :: isfvctr_par, nfvctr_par
      integer,dimension(nfvctr),intent(in) :: istsegline
      integer,dimension(nseg),intent(in) :: keyv
      integer,intent(out) :: isvctr, nvctrp
      integer,dimension(0:nproc-1),intent(out) :: isvctr_par, nvctr_par
      
      ! Local variables
      integer :: jproc, jst_line, jst_seg

      call f_routine(id='init_matrix_parallelization')

      ! parallelization of matrices, following same idea as norb/norbp/isorb
      !most equal distribution, but want corresponding to norbp for second column

      !$omp parallel default(none) &
      !$omp shared(nproc, isfvctr_par, nfvctr_par, nvctr, istsegline, keyv, isvctr_par) &
      !$omp private(jproc, jst_line, jst_seg)
      !$omp do
      do jproc=0,nproc-1
          jst_line = isfvctr_par(jproc)+1
          if (nfvctr_par(jproc)==0) then
             isvctr_par(jproc) = nvctr
          else
             jst_seg = istsegline(jst_line)
             isvctr_par(jproc) = keyv(jst_seg)-1
          end if
      end do
      !$omp end do
      !$omp end parallel

      do jproc=0,nproc-1
         if (jproc==nproc-1) then
            nvctr_par(jproc)=nvctr-isvctr_par(jproc)
         else
            nvctr_par(jproc)=isvctr_par(jproc+1)-isvctr_par(jproc)
         end if
         if (iproc==jproc) isvctr=isvctr_par(jproc)
         if (iproc==jproc) nvctrp=nvctr_par(jproc)
      end do

      call f_release_routine()

    end subroutine init_matrix_parallelization


    subroutine init_matrix_taskgroups(iproc, nproc, nat, parallel_layout, collcom, collcom_sr, smat, iirow, iicol)
      use communications_base, only: comms_linear
      use sparsematrix_base, only: sparse_matrix
      use yaml_output
      implicit none

      ! Caling arguments
      integer,intent(in) :: iproc, nproc, nat
      logical,intent(in) :: parallel_layout
      type(comms_linear),intent(in) :: collcom, collcom_sr
      type(sparse_matrix),intent(inout) :: smat
      integer,dimension(2),intent(in) :: iirow, iicol

      ! Local variables
      integer :: ipt, ii, i0, i0i, iiorb, j, i0j, jjorb, ind, ind_min, ind_max, iseq
      integer :: ntaskgroups, jproc, itaskgroups
      integer :: nfvctrp, isfvctr, isegstart, isegend, jorb, istart, iend, iistg, iietg, itg
      integer, dimension(:,:), allocatable :: iuse_startend, itaskgroups_startend
      integer, dimension(:), allocatable :: tasks_per_taskgroup
      integer :: ntaskgrp_calc, ntaskgrp_use, i, ncount, iitaskgroup, group, ierr, iitaskgroups, newgroup, iseg
      !logical :: go_on
      integer,dimension(:,:),allocatable :: in_taskgroup
      integer :: iproc_start, iproc_end, imin, imax
      logical :: found, found_start, found_end
      !integer :: jstart, kkproc, kproc, jend, lproc, llproc
      !integer :: iprocstart_current, iprocend_current, iprocend_prev, iprocstart_next
      integer :: irow, icol, inc, ist, ind_min1, ind_max1
      integer,dimension(:),pointer :: isvctr_par, nvctr_par
      logical, parameter :: print_full=.false.
      integer,dimension(:),pointer :: moduloarray


      call f_routine(id='init_matrix_taskgroups')
      call timing(iproc,'inittaskgroup','ON')

      ! First determine the minimal and maximal value oft the matrix which is used by each process
      iuse_startend = f_malloc0((/1.to.2,0.to.nproc-1/),id='iuse_startend')


      ! The matrices can be parallelized

      ind_min = smat%nvctr
      ind_max = 0

      ! The operations done in the transposed wavefunction layout
      !call check_transposed_layout()
      call get_modulo_array(smat, moduloarray)
      call find_minmax_transposed(smat%matrixindex_in_compressed_fortransposed,collcom,smat%nfvctr,moduloarray,ind_min,ind_max)
      call find_startendseg_transposed(ind_min,ind_max,smat)


      ! Now check the compress_distributed layout
      call check_compress_distributed_layout(smat,ind_min,ind_max)

      ! Now check the matrix matrix multiplications layout
      call check_matmul_layout(smat%smmm%nseq,smat%smmm%indices_extract_sequential,ind_min,ind_max)
      !!write(*,'(a,3i8)') 'after check_matmul: iproc, ind_min, ind_max', iproc, ind_min, ind_max

      ! Now check the sumrho operations
      !call check_sumrho_layout()
      call check_sumrho_layout(collcom_sr,smat%nfvctr,moduloarray,smat%matrixindex_in_compressed_fortransposed,ind_min,ind_max)
      call f_free_ptr(moduloarray)

      ! Now check the pseudo-exact orthonormalization during the input guess
      !call check_ortho_inguess()
      call check_ortho_inguess(smat,ind_min,ind_max)

      ! Now check the submatrix extraction for the projector charge analysis
      call check_projector_charge_analysis(iproc, nproc, nat, smat, ind_min, ind_max)


      ind_min1 = ind_min
      ind_max1 = ind_max

      !!write(*,'(a,3i8)') 'after init: iproc, ind_min1, ind_max1', iproc, ind_min1, ind_max1

      !@ NEW #####################################################################
      !@ Make sure that the min and max are at least as large as the reference
      do i=1,2
          if (i==1) then
              istart = 1
              iend = smat%nfvctr
              inc = 1
          else
              istart = smat%nfvctr
              iend = 1
              inc = -1
              !!write(*,*) 'iproc, iirow(i)', iproc, iirow(i)
          end if
          search_out: do irow=iirow(i),iend,inc
              if (irow==iirow(i)) then
                  ist = iicol(i)
              else
                  ist = istart
              end if
              do icol=ist,iend,inc
                  ii = matrixindex_in_compressed(smat, icol, irow)
                  if (ii>0) then
                      if (i==1) then
                          ind_min = ii
                      else
                          ind_max = ii
                      end if
                      exit search_out
                  end if
              end do
          end do search_out
      end do
      if (ind_min>ind_min1) then
          write(*,*) 'ind_min, ind_min1', ind_min, ind_min1
          stop 'ind_min>ind_min1'
      end if
      if (ind_max<ind_max1) then
          write(*,*) 'ind_max, ind_max1', ind_max, ind_max1
          stop 'ind_max<ind_max1'
      end if
      !!write(*,'(a,i3,3x,2(2i6,4x))') 'iproc, ind_min, ind_max, ind_min1, ind_max1', iproc,  ind_min, ind_max,  ind_min1, ind_max1
      !@ END NEW #################################################################


      if (.not.parallel_layout) then
          ! The matrices can not be parallelized
          ind_min = 1
          ind_max = smat%nvctr
      end if

      ! Enlarge the values if necessary such that they always start and end with a complete segment
      do iseg=1,smat%nseg
          istart = smat%keyv(iseg)
          iend = smat%keyv(iseg) + smat%keyg(2,1,iseg)-smat%keyg(1,1,iseg)
          if (istart<=ind_min .and. ind_min<=iend) then
              ind_min = istart
          end if
          if (istart<=ind_max .and. ind_max<=iend) then
              ind_max = iend
          end if
      end do

      ! Now the minimal and maximal values are known
      iuse_startend(1,iproc) = ind_min
      iuse_startend(2,iproc) = ind_max
      if (nproc>1) then
          call mpiallred(iuse_startend(1,0), 2*nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! Make sure that the used parts are always "monotonically increasing"
      do jproc=nproc-2,0,-1
          ! The start of part jproc must not be greater than the start of part jproc+1
          iuse_startend(1,jproc) = min(iuse_startend(1,jproc),iuse_startend(1,jproc+1)) 
      end do
      do jproc=1,nproc-1
          ! The end of part jproc must not be smaller than the end of part jproc-1
          iuse_startend(2,jproc) = max(iuse_startend(2,jproc),iuse_startend(2,jproc-1)) 
      end do


      !!smat%istartend_local(1) = ind_min
      !!smat%istartend_local(2) = ind_max
      smat%istartend_local(1) = iuse_startend(1,iproc)
      smat%istartend_local(2) = iuse_startend(2,iproc)

      ! Check to which segments these values belong
      found_start = .false.
      found_end = .false.
      do iseg=1,smat%nseg
          if (smat%keyv(iseg)==smat%istartend_local(1)) then
              smat%istartendseg_local(1) = iseg
              found_start = .true.
          end if
          if (smat%keyv(iseg)+smat%keyg(2,1,iseg)-smat%keyg(1,1,iseg)==smat%istartend_local(2)) then
              smat%istartendseg_local(2) = iseg
              found_end = .true.
          end if
      end do
      if (.not.found_start) stop 'segment corresponding to smat%istartend_local(1) not found!'
      if (.not.found_end) stop 'segment corresponding to smat%istartend_local(2) not found!'


      !!if (iproc==0)  then
      !!    do jproc=0,nproc-1
      !!        call yaml_map('iuse_startend',(/jproc,iuse_startend(1:2,jproc)/))
      !!    end do
      !!end if
 
      !if (iproc==0) write(*,'(a,100(2i7,4x))') 'iuse_startend',iuse_startend
 
!!      ntaskgroups = 1
!!      llproc=0 !the first task of the current taskgroup
!!      ii = 0
!!      do 
!!          jproc = llproc + ii
!!          if (jproc==nproc-1) exit
!!          jstart = iuse_startend(1,jproc) !beginning of part used by task jproc
!!          jend = iuse_startend(2,jproc) !end of part used by task jproc
!!          ii = ii + 1
!!          !!!search the last process whose part ends prior to iend
!!          !!go_on = .true.
!!          !!do lproc=nproc-1,0,-1
!!          !!    if (iuse_startend(1,lproc)<=jend) then
!!          !!        !if (iproc==0) write(*,'(a,3i8)') 'lproc, iuse_startend(1,lproc), iuse_startend(2,llproc)', lproc, iuse_startend(1,lproc), iuse_startend(2,llproc)
!!          !!        if (iuse_startend(1,lproc)<=iuse_startend(2,llproc)) then
!!          !!            go_on = .false.
!!          !!        end if
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          !if (iproc==0) write(*,*) '2: llproc, ii, jproc, go_on', llproc, ii, jproc, go_on
!!          ! Make sure that the beginning of the part used by jproc is larger than
!!          ! the end of the part used by llproc (which is the first task of the current taskgroup)
!!          if (iuse_startend(1,jproc)<=iuse_startend(2,llproc)) then
!!              cycle
!!          end if
!!          ntaskgroups = ntaskgroups + 1
!!          !!! Search the starting point of the next taskgroups, defined as the
!!          !!! largest starting part which is smaller than jend
!!          !!llproc=nproc-1
!!          !!do lproc=nproc-1,0,-1
!!          !!    if (iuse_startend(1,lproc)<=jend) then
!!          !!        llproc = lproc
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          ! jproc is now the start of the new taskgroup
!!          llproc=jproc
!!          !if (iproc==0) write(*,*) 'llproc, ii, jproc, ntaskgroups', llproc, ii, jproc, ntaskgroups
!!          ii = 0
!!          !if (llproc==nproc-1) exit
!!      end do
!!      !if (iproc==0) write(*,*) 'iproc, ntaskgroups', iproc, ntaskgroups

!!      !@NEW ###################
!!      ntaskgroups = 1
!!      iproc_start = 0
!!      iproc_end = 0
!!      do
!!          ! Search the first process whose parts does not overlap any more with
!!          ! the end of the first task of the current taskgroup.
!!          ! This will be the last task of the current taskgroup.
!!          found = .false.
!!          do jproc=iproc_start,nproc-1
!!              !if (iproc==0) write(*,'(a,2i8)') 'iuse_startend(1,jproc), iuse_startend(2,iproc_start)', iuse_startend(1,jproc), iuse_startend(2,iproc_start)
!!              if (iuse_startend(1,jproc)>iuse_startend(2,iproc_start)) then
!!                  iproc_end = jproc
!!                  found = .true.
!!                  exit
!!              end if
!!          end do
!!          if (.not.found) exit
!!          !iproc_end = iproc_start
!!
!!          !!! Search the last process whose part overlaps with the end of the current taskgroup.
!!          !!! This will be the first task of the next taskgroup.
!!          !!found = .false.
!!          !!do jproc=nproc-1,0,-1
!!          !!    !if (iproc==0) write(*,'(a,2i8)') 'iuse_startend(1,jproc), iuse_startend(2,iproc_end)', iuse_startend(1,jproc), iuse_startend(2,iproc_end)
!!          !!    if (iuse_startend(1,jproc)<=iuse_startend(2,iproc_end)) then
!!          !!        ntaskgroups = ntaskgroups + 1
!!          !!        iproc_start = jproc
!!          !!        found = .true.
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          !!if (iproc==0) write(*,*) 'iproc_start, iproc_end', iproc_start, iproc_end
!!          !!if (.not.found) exit
!!          !!if (iproc_start==nproc-1) exit
!!          ! Search the last process whose part overlaps with the start of the current taskgroup.
!!          ! This will be the first task of the next taskgroup.
!!          found = .false.
!!          !do jproc=nproc-1,0,-1
!!          do jproc=0,nproc-1
!!              !if (iproc==0) write(*,'(a,2i8)') 'iuse_startend(1,jproc), iuse_startend(2,iproc_end)', iuse_startend(1,jproc), iuse_startend(2,iproc_end)
!!              if (iuse_startend(1,jproc)>iuse_startend(1,iproc_end)) then
!!                  ntaskgroups = ntaskgroups + 1
!!                  iproc_start = jproc
!!                  found = .true.
!!                  exit
!!              end if
!!          end do
!!          if (.not.found) exit
!!      end do
!!      !@END NEW ###############

!!      !@NEW2 ###############################################
!!      iprocstart_next = 0
!!      iprocstart_current = 0
!!      iprocend_current = 0
!!      ntaskgroups = 1
!!      do
!!          iprocend_prev = iprocend_current
!!          iprocstart_current = iprocstart_next 
!!          !itaskgroups_startend(1,itaskgroups) = iuse_startend(2,iprocstart_current)
!!          ! Search the first process whose part starts later than then end of the part of iprocend_prev. This will be the first task of
!!          ! the next taskgroup
!!          do jproc=0,nproc-1
!!             if (iuse_startend(1,jproc)>iuse_startend(2,iprocend_prev)) then
!!                 iprocstart_next = jproc
!!                 exit
!!             end if
!!          end do
!!          ! Search the first process whose part ends later than then the start of the part of iprocstart_next. This will be the last task of
!!          ! the current taskgroup
!!          do jproc=0,nproc-1
!!             if (iuse_startend(2,jproc)>iuse_startend(1,iprocstart_next)) then
!!                 iprocend_current = jproc
!!                 exit
!!             end if
!!          end do
!!          !itaskgroups_startend(2,itaskgroups) = iuse_startend(2,iprocend_current)
!!          if (iproc==0) write(*,'(a,4i5)') 'iprocend_prev, iprocstart_current, iprocend_current, iprocstart_next', iprocend_prev, iprocstart_current, iprocend_current, iprocstart_next
!!          if (iprocstart_current==nproc-1) exit
!!          ntaskgroups = ntaskgroups + 1
!!      end do
!!      !@END NEW2 ###########################################


        !@NEW3 #############################################
        ntaskgroups = 1
        iproc_start = 0
        iproc_end = 0
        ii = 0
        do

            ! Search the first task whose part starts after the end of the part of the reference task
            found = .false.
            do jproc=0,nproc-1
                if (iuse_startend(1,jproc)>iuse_startend(2,iproc_end)) then
                    iproc_start = jproc
                    found = .true.
                    exit
                end if
            end do

            ! If this search was successful, start a new taskgroup
            if (found) then
                ! Determine the reference task, which is the last task whose part starts before the end of the current taskgroup
                ii = iuse_startend(2,iproc_start-1)
                do jproc=nproc-1,0,-1
                    if (iuse_startend(1,jproc)<=ii) then
                        iproc_end = jproc
                        exit
                    end if
                end do
                ! Increase the number of taskgroups
                ntaskgroups = ntaskgroups + 1
            else
                exit
            end if

        end do
        !@END NEW3 #########################################

      smat%ntaskgroup = ntaskgroups
 
      itaskgroups_startend = f_malloc0((/2,ntaskgroups/),id='itaskgroups_startend')
!!      itaskgroups_startend(1,1) = 1
!!      itaskgroups = 1
!!      llproc=0
!!      ii = 0
!!      do 
!!          jproc = llproc + ii
!!          if (jproc==nproc-1) exit
!!          jstart = iuse_startend(1,jproc) !beginning of part used by task jproc
!!          jend = iuse_startend(2,jproc) !end of part used by task jproc
!!          ii = ii + 1
!!          !!!search the last process whose part ends prior to jend
!!          !!go_on = .true.
!!          !!do lproc=nproc-1,0,-1
!!          !!    if (iuse_startend(1,lproc)<=jend) then
!!          !!        if (iuse_startend(1,lproc)<=iuse_startend(2,llproc)) then
!!          !!            go_on = .false.
!!          !!        end if
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          ! Make sure that the beginning of the part used by jproc is larger than
!!          ! the end of the part used by llproc (which is the first task of the current taskgroup)
!!          if (iuse_startend(1,jproc)<=iuse_startend(2,llproc)) then
!!              cycle
!!          end if
!!          !itaskgroups_startend(2,itaskgroups) = jend
!!          ! The end of the taskgroup is the end of the first task whose end is
!!          ! above the start of the new taskgroup
!!          do lproc=0,nproc-1
!!              if (iuse_startend(2,lproc)>jstart) then
!!                  itaskgroups_startend(2,itaskgroups) = iuse_startend(2,lproc)
!!                  exit
!!              end if
!!          end do
!!          itaskgroups = itaskgroups + 1
!!          !!! Search the starting point of the next taskgroups, defined as the
!!          !!! largest starting part which is smaller than jend
!!          !!llproc=nproc-1
!!          !!do lproc=nproc-1,0,-1
!!          !!    if (iuse_startend(1,lproc)<=jend) then
!!          !!        itaskgroups_startend(1,itaskgroups) = iuse_startend(1,lproc)
!!          !!        llproc = lproc
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          ! jproc is now the start of the new taskgroup
!!          llproc=jproc
!!          itaskgroups_startend(1,itaskgroups) = iuse_startend(1,llproc)
!!          ii = 0
!!          !if (llproc==nproc-1) exit
!!      end do
!!      itaskgroups_startend(2,itaskgroups) = iuse_startend(2,nproc-1)


!!      !@NEW ###################
!!      itaskgroups = 1
!!      iproc_start = 0
!!      iproc_end = 0
!!      itaskgroups_startend(1,1) = 1
!!      do
!!          ! Search the first process whose parts does not overlap any more with
!!          ! the end of the first task of the current taskgroup.
!!          ! This will be the last task of the current taskgroup.
!!          found = .false.
!!          do jproc=iproc_start,nproc-1
!!              if (iuse_startend(1,jproc)>iuse_startend(2,iproc_start)) then
!!                  iproc_end = jproc
!!                  itaskgroups_startend(2,itaskgroups) = iuse_startend(2,jproc)
!!                  found = .true.
!!                  exit
!!              end if
!!          end do
!!          if (.not.found) exit
!!          !!iproc_end = iproc_start
!!          !!itaskgroups_startend(2,itaskgroups) = iuse_startend(2,iproc_end)
!!
!!          ! Search the last process whose part overlaps with the end of the current taskgroup.
!!          ! This will be the first task of the next taskgroup.
!!          found = .false.
!!          !do jproc=nproc-1,0,-1
!!          do jproc=0,nproc-1
!!              if (iuse_startend(1,jproc)>iuse_startend(1,iproc_end)) then
!!                  itaskgroups = itaskgroups + 1
!!                  iproc_start = jproc
!!                  itaskgroups_startend(1,itaskgroups) = iuse_startend(1,jproc)
!!                  found = .true.
!!                  exit
!!              end if
!!          end do
!!          if (.not.found) exit
!!          !!!!if (iproc_start==nproc-1) exit
!!          !!! Search the last process whose part overlaps with the end of the current taskgroup.
!!          !!! This will be the first task of the next taskgroup.
!!          !!found = .false.
!!          !!!do jproc=0,nproc-1
!!          !!do jproc=0,nproc-1
!!          !!    if (iuse_startend(1,jproc)>iuse_startend(1,iproc_end)) then
!!          !!        itaskgroups = itaskgroups + 1
!!          !!        iproc_start = jproc
!!          !!        itaskgroups_startend(1,itaskgroups) = iuse_startend(1,jproc)
!!          !!        found = .true.
!!          !!        exit
!!          !!    end if
!!          !!end do
!!          !!if (.not.found) exit
!!      end do
!!      itaskgroups_startend(2,itaskgroups) = smat%nvctr
!!      !@END NEW ###############

!!      !@NEW2 ###############################################
!!      iprocstart_next = 0
!!      iprocstart_current = 0
!!      iprocend_current = 0
!!      itaskgroups = 1
!!      do
!!          iprocend_prev = iprocend_current
!!          iprocstart_current = iprocstart_next 
!!          itaskgroups_startend(1,itaskgroups) = iuse_startend(1,iprocstart_current)
!!          ! Search the first process whose part starts later than then end of the part of iprocend_prev. This will be the first task of
!!          ! the next taskgroup
!!          do jproc=0,nproc-1
!!             if (iuse_startend(1,jproc)>iuse_startend(2,iprocend_prev)) then
!!                 iprocstart_next = jproc
!!                 exit
!!             end if
!!          end do
!!          ! Search the first process whose part ends later than then the start of the part of iprocstart_next. This will be the last task of
!!          ! the current taskgroup
!!          do jproc=0,nproc-1
!!             if (iuse_startend(2,jproc)>iuse_startend(1,iprocstart_next)) then
!!                 iprocend_current = jproc
!!                 exit
!!             end if
!!          end do
!!          itaskgroups_startend(2,itaskgroups) = iuse_startend(2,iprocend_current)
!!          if (iprocstart_current==nproc-1) exit
!!          itaskgroups = itaskgroups +1
!!      end do
!!      !@END NEW2 ###########################################


        !@NEW3 #############################################
        itaskgroups = 1
        iproc_start = 0
        iproc_end = 0
        itaskgroups_startend(1,1) = 1
        do

            ! Search the first task whose part starts after the end of the part of the reference task
            found = .false.
            do jproc=0,nproc-1
                if (iuse_startend(1,jproc)>iuse_startend(2,iproc_end)) then
                    iproc_start = jproc
                    found = .true.
                    exit
                end if
            end do


            ! If this search was successful, start a new taskgroup
            if (found) then
                ! Store the end of the current taskgroup
                itaskgroups_startend(2,itaskgroups) = iuse_startend(2,iproc_start-1)
                ! Determine the reference task, which is the last task whose part starts before the end of the current taskgroup
                do jproc=nproc-1,0,-1
                    if (iuse_startend(1,jproc)<=itaskgroups_startend(2,itaskgroups)) then
                        iproc_end = jproc
                        exit
                    end if
                end do
                ! Increase the number of taskgroups
                itaskgroups = itaskgroups + 1
                ! Store the beginning of the new taskgroup
                itaskgroups_startend(1,itaskgroups) = iuse_startend(1,iproc_start)
            else
                ! End of the taskgroup if the search was not successful
                itaskgroups_startend(2,itaskgroups) = iuse_startend(2,nproc-1)
                exit
            end if

        end do
        !@END NEW3 #########################################

      !!if (iproc==0)  then
      !!    do jproc=1,smat%ntaskgroup
      !!        call yaml_map('itaskgroups_startend',itaskgroups_startend(1:2,jproc))
      !!    end do
      !!end if
      !call yaml_flash_document()
      call mpi_barrier(bigdft_mpi%mpi_comm,jproc)



      if (itaskgroups/=ntaskgroups) stop 'itaskgroups/=ntaskgroups'
      !if (iproc==0) write(*,'(a,i8,4x,1000(2i7,4x))') 'iproc, itaskgroups_startend', itaskgroups_startend
 
      ! Assign the processes to the taskgroups
      ntaskgrp_calc = 0
      ntaskgrp_use = 0
      do itaskgroups=1,ntaskgroups
          if ( iuse_startend(1,iproc)<=itaskgroups_startend(2,itaskgroups) .and.  &
               iuse_startend(2,iproc)>=itaskgroups_startend(1,itaskgroups) ) then
              !!write(*,'(2(a,i0))') 'USE: task ',iproc,' is in taskgroup ',itaskgroups
               ntaskgrp_use = ntaskgrp_use + 1
          end if
      end do
      if (ntaskgrp_use>2) stop 'ntaskgrp_use>2'

      smat%ntaskgroupp = max(ntaskgrp_calc,ntaskgrp_use)

      smat%taskgroup_startend = f_malloc_ptr((/2,2,smat%ntaskgroup/),id='smat%taskgroup_startend')
      smat%taskgroupid = f_malloc_ptr((/smat%ntaskgroupp/),id='smat%smat%taskgroupid')
      smat%inwhichtaskgroup = f_malloc0_ptr((/1.to.2,0.to.nproc-1/),id='smat%smat%inwhichtaskgroup')


      i = 0
      do itaskgroups=1,smat%ntaskgroup
          i = i + 1
          smat%taskgroup_startend(1,1,i) = itaskgroups_startend(1,itaskgroups)
          smat%taskgroup_startend(2,1,i) = itaskgroups_startend(2,itaskgroups)
      end do
      if (i/=smat%ntaskgroup) then
          write(*,*) 'i, smat%ntaskgroup', i, smat%ntaskgroup
          stop 'i/=smat%ntaskgroup'
      end if



      i = 0
      do itaskgroups=1,smat%ntaskgroup
          if( iuse_startend(1,iproc)<=itaskgroups_startend(2,itaskgroups) .and.  &
               iuse_startend(2,iproc)>=itaskgroups_startend(1,itaskgroups) ) then
               i = i + 1
               smat%taskgroupid(i) = itaskgroups
               smat%inwhichtaskgroup(i,iproc) = itaskgroups
          end if
      end do
      if (i/=smat%ntaskgroupp) then
          write(*,*) 'i, smat%ntaskgroupp', i, smat%ntaskgroupp
          stop 'i/=smat%ntaskgroupp'
      end if

      if (nproc>1) then
          call mpiallred(smat%inwhichtaskgroup(1,0), 2*nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      ! Partition the entire matrix in disjoint submatrices
      smat%taskgroup_startend(1,2,1) = smat%taskgroup_startend(1,1,1)
      do itaskgroups=2,smat%ntaskgroup
          smat%taskgroup_startend(1,2,itaskgroups) = &
              (smat%taskgroup_startend(2,1,itaskgroups-1)+smat%taskgroup_startend(1,1,itaskgroups)) / 2
          smat%taskgroup_startend(2,2,itaskgroups-1) = smat%taskgroup_startend(1,2,itaskgroups)-1
      end do
      smat%taskgroup_startend(2,2,smat%ntaskgroup) = smat%taskgroup_startend(2,1,smat%ntaskgroup)

      !if (iproc==0) write(*,'(a,1000(2i8,4x))') 'iproc, smat%taskgroup_startend(:,2,:)',smat%taskgroup_startend(:,2,:)

      ! Some checks
      ncount = 0
      do itaskgroups=1,smat%ntaskgroup
          ncount = ncount + smat%taskgroup_startend(2,2,itaskgroups)-smat%taskgroup_startend(1,2,itaskgroups)+1
          if (itaskgroups>1) then
              if (smat%taskgroup_startend(1,1,itaskgroups)>smat%taskgroup_startend(2,1,itaskgroups-1)) then
                  stop 'smat%taskgroup_startend(1,1,itaskgroups)>smat%taskgroup_startend(2,1,itaskgroups-1)'
              end if
          end if
      end do
      if (ncount/=smat%nvctr) then
          write(*,*) 'ncount, smat%nvctr', ncount, smat%nvctr
          stop 'ncount/=smat%nvctr'
      end if

      ! Check that the data that task iproc needs is really contained in the
      ! taskgroups to which iproc belongs.
      imin=smat%nvctr
      imax=1
      do itaskgroups=1,smat%ntaskgroupp
          iitaskgroup = smat%taskgroupid(itaskgroups)
          imin = min(imin,smat%taskgroup_startend(1,1,iitaskgroup))
          imax = max(imax,smat%taskgroup_startend(2,1,iitaskgroup))
      end do
      if (iuse_startend(1,iproc)<imin) then
          write(*,*) 'iuse_startend(1,iproc),imin', iuse_startend(1,iproc),imin
          stop 'iuse_startend(1,iproc)<imin'
      end if
      if (iuse_startend(2,iproc)>imax) then
          write(*,*) 'iuse_startend(2,iproc),imax', iuse_startend(2,iproc),imax
          stop 'iuse_startend(2,iproc)>imax'
      end if


      ! Assign the values of nvctrp_tg and iseseg_tg
      ! First and last segment of the matrix
      iistg=smat%taskgroupid(1) !first taskgroup of task iproc
      iietg=smat%taskgroupid(smat%ntaskgroupp) !last taskgroup of task iproc
      found_start = .false.
      found_end = .false.
      do iseg=1,smat%nseg
          if (smat%keyv(iseg)==smat%taskgroup_startend(1,1,iistg)) then
              smat%iseseg_tg(1) = iseg
              smat%isvctrp_tg = smat%keyv(iseg)-1
              found_start = .true.
          end if
          if (smat%keyv(iseg)+smat%keyg(2,1,iseg)-smat%keyg(1,1,iseg)==smat%taskgroup_startend(2,1,iietg)) then
              smat%iseseg_tg(2) = iseg
              found_end = .true.
          end if
      end do
      if (.not.found_start) stop 'first segment of taskgroup matrix not found'
      if (.not.found_end) stop 'last segment of taskgroup matrix not found'
      ! Size of the matrix
      smat%nvctrp_tg = smat%taskgroup_startend(2,1,iietg) - smat%taskgroup_startend(1,1,iistg) + 1




      ! Create the taskgroups
      ! Count the number of tasks per taskgroup
      tasks_per_taskgroup = f_malloc0(smat%ntaskgroup,id='tasks_per_taskgroup')
      do itaskgroups=1,smat%ntaskgroupp
          iitaskgroup = smat%taskgroupid(itaskgroups)
          tasks_per_taskgroup(iitaskgroup) = tasks_per_taskgroup(iitaskgroup) + 1
      end do
      if (nproc>1) then
          call mpiallred(tasks_per_taskgroup, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      !if (iproc==0) write(*,'(a,i7,4x,1000i7)') 'iproc, tasks_per_taskgroup', iproc, tasks_per_taskgroup
      !call mpi_comm_group(bigdft_mpi%mpi_comm, group, ierr)
      group=mpigroup(bigdft_mpi%mpi_comm)

      in_taskgroup = f_malloc0((/0.to.nproc-1,1.to.smat%ntaskgroup/),id='in_taskgroup')
      smat%tgranks = f_malloc_ptr((/0.to.maxval(tasks_per_taskgroup)-1,1.to.smat%ntaskgroup/),id='smat%tgranks')
      smat%nranks = f_malloc_ptr(smat%ntaskgroup,id='smat%nranks')
      !smat%isrank = f_malloc_ptr(smat%ntaskgroup,id='smat%isrank')

      ! number of tasks per taskgroup
      do itg=1,smat%ntaskgroup
          smat%nranks(itg) = tasks_per_taskgroup(itg)
      end do

      do itaskgroups=1,smat%ntaskgroupp
          iitaskgroups = smat%taskgroupid(itaskgroups)
          in_taskgroup(iproc,iitaskgroups) = 1
      end do
      if (nproc>1) then
          call mpiallred(in_taskgroup, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if

      allocate(smat%mpi_groups(smat%ntaskgroup))
      do itaskgroups=1,smat%ntaskgroup
          smat%mpi_groups(itaskgroups) = mpi_environment_null()
      end do
      do itaskgroups=1,smat%ntaskgroup
          ii = 0
          do jproc=0,nproc-1
              if (in_taskgroup(jproc,itaskgroups)>0) then
                  smat%tgranks(ii,itaskgroups) = jproc
                  ii = ii + 1
              end if
          end do
          ! Store the ID of the first task of each taskgroup
          !smat%isrank(itaskgroups) = smat%tgranks(1,itaskgroups)
          if (ii/=tasks_per_taskgroup(itaskgroups)) stop 'ii/=tasks_per_taskgroup(itaskgroups)' !for debugging
          call mpi_env_create_group(itaskgroups,smat%ntaskgroup,bigdft_mpi%mpi_comm,&
               group,ii,smat%tgranks(:,itaskgroups),smat%mpi_groups(itaskgroups))
!!$          call mpi_group_incl(group, ii, smat%tgranks(0,itaskgroups), newgroup, ierr)
!!$          call mpi_comm_create(bigdft_mpi%mpi_comm, newgroup, smat%mpi_groups(itaskgroups)%mpi_comm, ierr)
!!$          if (smat%mpi_groups(itaskgroups)%mpi_comm/=MPI_COMM_NULL) then
!!$              call mpi_comm_size(smat%mpi_groups(itaskgroups)%mpi_comm, smat%mpi_groups(itaskgroups)%nproc, ierr)
!!$              call mpi_comm_rank(smat%mpi_groups(itaskgroups)%mpi_comm, smat%mpi_groups(itaskgroups)%iproc, ierr)
!!$          end if
!!$          smat%mpi_groups(itaskgroups)%igroup = itaskgroups
!!$          smat%mpi_groups(itaskgroups)%ngroup = smat%ntaskgroup
!!$          call mpi_group_free(newgroup, ierr)
      end do
      call mpi_group_free(group, ierr)

      !do itaskgroups=1,smat%ntaskgroup
      !    if (smat%mpi_groups(itaskgroups)%iproc==0) write(*,'(2(a,i0))') 'process ',iproc,' is first in taskgroup ',itaskgroups 
      !end do

      ! Print a summary
      if (iproc==0) then
          call yaml_mapping_open('taskgroup summary')
          call yaml_map('number of taskgroups',smat%ntaskgroup)
          call yaml_sequence_open('taskgroups overview')
          do itaskgroups=1,smat%ntaskgroup
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('number of tasks',tasks_per_taskgroup(itaskgroups))
              !call yaml_map('IDs',smat%tgranks(0:tasks_per_taskgroup(itaskgroups)-1,itaskgroups))
              if (print_full) then
                  call yaml_mapping_open('IDs')
                  do itg=0,tasks_per_taskgroup(itaskgroups)-1
                      call yaml_mapping_open(yaml_toa(smat%tgranks(itg,itaskgroups),fmt='(i0)'))
                      call yaml_map('s',iuse_startend(1,smat%tgranks(itg,itaskgroups)))
                      call yaml_map('e',iuse_startend(2,smat%tgranks(itg,itaskgroups)))
                      call yaml_mapping_close()
                  end do
                  call yaml_mapping_close()
                  call yaml_newline()
              end if
              call yaml_map('start / end',smat%taskgroup_startend(1:2,1,itaskgroups))
              call yaml_map('start / end disjoint',smat%taskgroup_startend(1:2,2,itaskgroups))
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()
      end if


      ! Initialize a "local compress" from the matrix matrix multiplication layout
      !!!smat%smmm%ncl_smmm = 0
      !!!if (smat%smmm%nfvctrp>0) then
      !!!    isegstart=smat%istsegline(smat%smmm%isfvctr+1)
      !!!    isegend=smat%istsegline(smat%smmm%isfvctr+smat%smmm%nfvctrp)+smat%nsegline(smat%smmm%isfvctr+smat%smmm%nfvctrp)-1
      !!!    do iseg=isegstart,isegend
      !!!        ! A segment is always on one line, therefore no double loop
      !!!        do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
      !!!            smat%smmm%ncl_smmm = smat%smmm%ncl_smmm + 1
      !!!        end do
      !!!    end do
      !!!end if
      !!!if (smat%smmm%ncl_smmm/=smat%smmm%nvctrp_mm) then
      !!!    write(*,*) 'smat%smmm%ncl_smmm, smat%smmm%nvctrp_mm', smat%smmm%ncl_smmm, smat%smmm%nvctrp_mm
      !!!    stop
      !!!end if

      do i=1,2
          if (i==1) then
              isvctr_par => smat%smmm%isvctr_mm_par
              nvctr_par => smat%smmm%nvctr_mm_par
          else if (i==2) then
              isvctr_par => smat%isvctr_par
              nvctr_par => smat%nvctr_par
          end if

          !smat%smmm%nccomm_smmm = 0
          ii = 0
          do jproc=0,nproc-1
              !!istart = max(smat%istartend_local(1),smat%smmm%isvctr_mm_par(jproc)+1)
              !!iend = min(smat%istartend_local(2),smat%smmm%isvctr_mm_par(jproc)+smat%smmm%nvctr_mm_par(jproc))
              !!if (istart>iend) cycle
              !!smat%smmm%nccomm_smmm = smat%smmm%nccomm_smmm + 1
              istart = max(smat%istartend_local(1),isvctr_par(jproc)+1)
              iend = min(smat%istartend_local(2),isvctr_par(jproc)+nvctr_par(jproc))
              if (istart>iend) cycle
              ii = ii + 1
          end do

          if (i==1) then
              smat%smmm%nccomm_smmm = ii
              smat%smmm%luccomm_smmm = f_malloc_ptr((/4,smat%smmm%nccomm_smmm/),id='smat%smmm%luccomm_smmm')
          else if (i==2) then
              smat%nccomm = ii
              smat%luccomm = f_malloc_ptr((/4,smat%nccomm/),id='smatluccomm')
          end if

          !!smat%smmm%luccomm_smmm = f_malloc_ptr((/4,smat%smmm%nccomm_smmm/),id='smat%smmm%luccomm_smmm')
          ii = 0
          do jproc=0,nproc-1
              !!istart = max(smat%istartend_local(1),smat%smmm%isvctr_mm_par(jproc)+1)
              !!iend = min(smat%istartend_local(2),smat%smmm%isvctr_mm_par(jproc)+smat%smmm%nvctr_mm_par(jproc))
              !!if (istart>iend) cycle
              !!ii = ii + 1
              !!smat%smmm%luccomm_smmm(1,ii) = jproc !get data from this process
              !!smat%smmm%luccomm_smmm(2,ii) = istart-smat%smmm%isvctr_mm_par(jproc) !starting address on sending process
              !!smat%smmm%luccomm_smmm(3,ii) = istart-smat%isvctrp_tg !starting address on receiving process
              !!smat%smmm%luccomm_smmm(4,ii) = iend-istart+1 !number of elements
              istart = max(smat%istartend_local(1),isvctr_par(jproc)+1)
              iend = min(smat%istartend_local(2),isvctr_par(jproc)+nvctr_par(jproc))
              if (istart>iend) cycle
              ii = ii + 1
              if (i==1) then
                  smat%smmm%luccomm_smmm(1,ii) = jproc !get data from this process
                  smat%smmm%luccomm_smmm(2,ii) = istart-isvctr_par(jproc) !starting address on sending process
                  smat%smmm%luccomm_smmm(3,ii) = istart-smat%isvctrp_tg !starting address on receiving process
                  smat%smmm%luccomm_smmm(4,ii) = iend-istart+1 !number of elements
              else if (i==2) then
                  smat%luccomm(1,ii) = jproc !get data from this process
                  smat%luccomm(2,ii) = istart-isvctr_par(jproc) !starting address on sending process
                  smat%luccomm(3,ii) = istart-smat%isvctrp_tg !starting address on receiving process
                  smat%luccomm(4,ii) = iend-istart+1 !number of elements
              end if
          end do
      end do

      call f_free(in_taskgroup)
      call f_free(iuse_startend)
      call f_free(itaskgroups_startend)
      call f_free(tasks_per_taskgroup)
      call timing(iproc,'inittaskgroup','OF')
      call f_release_routine()


!!$      contains

!!$        subroutine check_transposed_layout()
!!$          logical :: found
!!$          integer :: iiseg1, iiseg2, iorb, jorb
!!$          integer,dimension(:),pointer :: moduloarray
!!$
!!$          call f_routine(id='check_transposed_layout')
!!$
!!$          call get_modulo_array(smat, moduloarray)
!!$
!!$          !$omp parallel &
!!$          !$omp default(none) &
!!$          !$omp shared(collcom, smat, moduloarray, ind_min, ind_max) &
!!$          !$omp private(ipt, ii, i0, i, i0i, iiorb, j, i0j, jjorb, ind, iorb, jorb)
!!$          !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$          do ipt=1,collcom%nptsp_c
!!$              ii=collcom%norb_per_gridpoint_c(ipt)
!!$              i0 = collcom%isptsp_c(ipt)
!!$              do i=1,ii
!!$                  i0i=i0+i
!!$                  iiorb=collcom%indexrecvorbital_c(i0i)
!!$                  iorb=moduloarray(iiorb)
!!$                  do j=1,ii
!!$                      i0j=i0+j
!!$                      jjorb=collcom%indexrecvorbital_c(i0j)
!!$                      jorb=moduloarray(jjorb)
!!$                      !ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
!!$                      ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                      !ind = get_transposed_index(smat,jjorb,iiorb)
!!$                      if (ind==0) write(*,'(a,2i8)') 'coarse iszero: iiorb, jjorb', iiorb, jjorb
!!$                      ind_min = min(ind_min,ind)
!!$                      ind_max = max(ind_max,ind)
!!$                  end do
!!$              end do
!!$          end do
!!$          !$omp end do
!!$          !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$          do ipt=1,collcom%nptsp_f
!!$              ii=collcom%norb_per_gridpoint_f(ipt)
!!$              i0 = collcom%isptsp_f(ipt)
!!$              do i=1,ii
!!$                  i0i=i0+i
!!$                  iiorb=collcom%indexrecvorbital_f(i0i)
!!$                  iorb=moduloarray(iiorb)
!!$                  do j=1,ii
!!$                      i0j=i0+j
!!$                      jjorb=collcom%indexrecvorbital_f(i0j)
!!$                      jorb=moduloarray(jjorb)
!!$                      !ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
!!$                      ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                      !ind = get_transposed_index(smat,jjorb,iiorb)
!!$                      if (ind==0) write(*,'(a,2i8)') 'fine iszero: iiorb, jjorb', iiorb, jjorb
!!$                      ind_min = min(ind_min,ind)
!!$                      ind_max = max(ind_max,ind)
!!$                  end do
!!$              end do
!!$          end do
!!$          !$omp end do
!!$          !$omp end parallel
!!$
!!$          ! Store these values
!!$          smat%istartend_t(1) = ind_min
!!$          smat%istartend_t(2) = ind_max
!!$
!!$          ! Determine to which segments this corresponds
!!$          iiseg1 = smat%nseg
!!$          iiseg2 = 1
!!$          !$omp parallel default(none) shared(smat, iiseg1, iiseg2) private(iseg, found)
!!$          found = .false.
!!$          !$omp do reduction(min: iiseg1)
!!$          do iseg=1,smat%nseg
!!$              ! A segment is always on one line
!!$              if (.not.found) then
!!$                  if (smat%keyv(iseg)+smat%keyg(2,1,iseg)-smat%keyg(1,1,iseg)>=smat%istartend_t(1)) then
!!$                      !smat%istartendseg_t(1)=iseg
!!$                      iiseg1=iseg
!!$                      found = .true.
!!$                  end if
!!$              end if
!!$          end do
!!$          !$omp end do
!!$          found = .false.
!!$          !$omp do reduction(max: iiseg2)
!!$          do iseg=smat%nseg,1,-1
!!$              if (.not.found) then
!!$                  if (smat%keyv(iseg)<=smat%istartend_t(2)) then
!!$                      !smat%istartendseg_t(2)=iseg
!!$                      iiseg2=iseg
!!$                      found = .true.
!!$                  end if
!!$              end if
!!$          end do
!!$          !$omp end do
!!$          !$omp end parallel
!!$          smat%istartendseg_t(1) = iiseg1
!!$          smat%istartendseg_t(2) = iiseg2
!!$
!!$          call f_free_ptr(moduloarray)
!!$
!!$          call f_release_routine()
!!$
!!$
!!$        end subroutine check_transposed_layout


        !function get_transposed_index(jorb,iorb) result(ind)
        !    integer,intent(in) :: jorb, iorb
        !    integer :: ind
        !    integer :: jjorb,iiorb
        !    ! If iorb is smaller than the offset, add a periodic shift
        !    if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
        !        iiorb = iorb + smat%nfvctr
        !    else
        !        iiorb = iorb
        !    end if
        !    if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
        !        jjorb = jorb + smat%nfvctr
        !    else
        !        jjorb = jorb
        !    end if
        !    ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
        !end function get_transposed_index


!!$        subroutine check_compress_distributed_layout()
!!$
!!$          call f_routine(id='check_compress_distributed_layout')
!!$
!!$          do i=1,2
!!$              if (i==1) then
!!$                  nfvctrp = smat%nfvctrp
!!$                  isfvctr = smat%isfvctr
!!$              else if (i==2) then
!!$                  nfvctrp = smat%smmm%nfvctrp
!!$                  isfvctr = smat%smmm%isfvctr
!!$              end if
!!$              if (nfvctrp>0) then
!!$                  isegstart=smat%istsegline(isfvctr+1)
!!$                  isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
!!$                  !$omp parallel default(none) &
!!$                  !$omp shared(isegstart, isegend, smat, ind_min, ind_max) &
!!$                  !$omp private(iseg, ii,jorb)
!!$                  !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$                  do iseg=isegstart,isegend
!!$                      ii=smat%keyv(iseg)-1
!!$                      ! A segment is always on one line, therefore no double loop
!!$                      do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$                          ii=ii+1
!!$                          ind_min = min(ii,ind_min)
!!$                          ind_max = max(ii,ind_max)
!!$                      end do
!!$                  end do
!!$                  !$omp end do
!!$                  !$omp end parallel
!!$              end if
!!$          end do
!!$
!!$          call f_release_routine()
!!$
!!$        end subroutine check_compress_distributed_layout


!!$        subroutine check_matmul_layout()
!!$
!!$          call f_routine(id='check_matmul_layout')
!!$
!!$          !$omp parallel default(none) shared(smat, ind_min, ind_max) private(iseq, ind)
!!$          !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$          do iseq=1,smat%smmm%nseq
!!$              ind=smat%smmm%indices_extract_sequential(iseq)
!!$              ind_min = min(ind_min,ind)
!!$              ind_max = max(ind_max,ind)
!!$          end do
!!$          !$omp end do
!!$          !$omp end parallel
!!$
!!$          call f_release_routine()
!!$
!!$        end subroutine check_matmul_layout

!!$        subroutine check_sumrho_layout()
!!$          integer :: iorb
!!$          integer,dimension(:),pointer :: moduloarray
!!$
!!$          call f_routine(id='check_sumrho_layout')
!!$
!!$          call get_modulo_array(smat, moduloarray)
!!$
!!$          !$omp parallel default(none) &
!!$          !$omp shared(collcom_sr, smat, moduloarray, ind_min, ind_max) private(ipt, ii, i0, i, iiorb, ind, iorb, jorb)
!!$          !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$          do ipt=1,collcom_sr%nptsp_c
!!$              ii=collcom_sr%norb_per_gridpoint_c(ipt)
!!$              i0=collcom_sr%isptsp_c(ipt)
!!$              do i=1,ii
!!$                  iiorb=collcom_sr%indexrecvorbital_c(i0+i)
!!$                  iorb=moduloarray(iiorb)
!!$                  !ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
!!$                  ind=smat%matrixindex_in_compressed_fortransposed(iorb,iorb)
!!$                  !ind=get_transposed_index(smat,iiorb,iiorb)
!!$                  ind_min = min(ind_min,ind)
!!$                  ind_max = max(ind_max,ind)
!!$              end do
!!$          end do
!!$          !$omp end do
!!$          !$omp end parallel
!!$
!!$          call f_free_ptr(moduloarray)
!!$
!!$          call f_release_routine()
!!$
!!$          !contains
!!$
!!$          !  function get_transposed_index(jorb,iorb) res(ind)
!!$          !      integer,intent(in) :: jorb, iorb
!!$          !      integer :: ind
!!$          !      integer :: jjorb,iiorb
!!$          !      ! If iorb is smaller than the offset, add a periodic shift
!!$          !      if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
!!$          !          iiorb = iorb + smat%nfvctr
!!$          !      else
!!$          !          iiorb = iorb
!!$          !      end if
!!$          !      if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
!!$          !          jjorb = jorb + smat%nfvctr
!!$          !      else
!!$          !          jjorb = jorb
!!$          !      end if
!!$          !      ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
!!$          !  end function get_transposed_index
!!$
!!$        end subroutine check_sumrho_layout


      !!  function get_start_of_segment(smat, iiseg) result(ist)

      !!      do iseg=smat%nseg,1,-1
      !!          if (iiseg>=smat%keyv(iseg)) then
      !!              it = smat%keyv(iseg)
      !!              exit
      !!          end if
      !!      end do

      !!  end function get_start_of_segment


!!$      subroutine check_ortho_inguess()
!!$        integer :: iorb, iiorb, isegstart, isegsend, iseg, j, i, jorb, korb, ind, nthread, ithread
!!$        logical,dimension(:,:),allocatable :: in_neighborhood
!!$        !$ integer :: omp_get_max_threads, omp_get_thread_num
!!$
!!$        call f_routine(id='check_ortho_inguess')
!!$
!!$        ! Allocate the array for all threads to avoid that it has to be declared private
!!$        nthread = 1
!!$        !$ nthread = omp_get_max_threads()
!!$        in_neighborhood = f_malloc((/1.to.smat%nfvctr,0.to.nthread-1/),id='in_neighborhood')
!!$        
!!$        ithread = 0
!!$        !$omp parallel default(none) &
!!$        !$omp shared(smat, in_neighborhood, ind_min, ind_max) &
!!$        !$omp private(iorb, iiorb, isegstart, isegend, iseg, j, jorb, korb, ind,i) &
!!$        !$omp firstprivate(ithread)
!!$        !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$        do iorb=1,smat%nfvctrp
!!$            !$ ithread = omp_get_thread_num()
!!$
!!$            iiorb = smat%isfvctr + iorb
!!$            isegstart = smat%istsegline(iiorb)
!!$            isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
!!$            in_neighborhood(:,ithread) = .false.
!!$            do iseg=isegstart,isegend
!!$                ! A segment is always on one line, therefore no double loop
!!$                j = smat%keyg(1,2,iseg)
!!$                do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$                    in_neighborhood(i,ithread) = .true.
!!$                end do
!!$            end do
!!$
!!$            do jorb=1,smat%nfvctr
!!$                if (.not.in_neighborhood(jorb,ithread)) cycle
!!$                do korb=1,smat%nfvctr
!!$                    if (.not.in_neighborhood(korb,ithread)) cycle
!!$                    ind = matrixindex_in_compressed(smat,korb,jorb)
!!$                    if (ind>0) then
!!$                        ind_min = min(ind_min,ind)
!!$                        ind_max = max(ind_max,ind)
!!$                    end if
!!$                end do
!!$            end do
!!$
!!$        end do
!!$        !$omp end do
!!$        !$omp end parallel
!!$
!!$        call f_free(in_neighborhood)
!!$
!!$        !!do iorb=1,smat%nfvctrp
!!$        !!    iiorb = smat%isfvctr + iorb
!!$        !!    isegstart = smat%istsegline(iiorb)
!!$        !!    isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
!!$        !!    do iseg=isegstart,isegend
!!$        !!        ! A segment is always on one line, therefore no double loop
!!$        !!        j = smat%keyg(1,2,iseg)
!!$        !!        do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$        !!            ind = matrixindex_in_compressed(smat,i,j)
!!$        !!            ind_min = min(ind_min,ind)
!!$        !!            ind_max = max(ind_max,ind)
!!$        !!        end do
!!$        !!    end do
!!$        !!end do
!!$
!!$        call f_release_routine()
!!$
!!$      end subroutine check_ortho_inguess

 
    end subroutine init_matrix_taskgroups




    subroutine check_local_matrix_extents(iproc, nproc, nat, collcom, collcom_sr, smat, irow, icol)
          use communications_base, only: comms_linear
          use sparsematrix_base, only: sparse_matrix
          use yaml_output
          implicit none
    
          ! Caling arguments
          integer,intent(in) :: iproc, nproc, nat
          type(comms_linear),intent(in) :: collcom, collcom_sr
          type(sparse_matrix),intent(in) :: smat
          integer,dimension(2),intent(out) :: irow, icol
    
          ! Local variables
          integer :: ind_min, ind_max, i, ii_ref, iorb, jorb, ii, iseg
          logical :: found

          real(kind=4) :: tr0, tr1, trt0, trt1
          real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
          logical, parameter :: extra_timing=.false.
          integer,dimension(:),pointer :: moduloarray
                        
    
          call timing(iproc,'matrix_extents','ON')
          if (extra_timing) call cpu_time(trt0)  

          ind_min = smat%nvctr
          ind_max = 0

          if (extra_timing) call cpu_time(tr0)
          ! The operations done in the transposed wavefunction layout
          !call check_transposed_layout()
          call get_modulo_array(smat, moduloarray)
          call find_minmax_transposed(smat%matrixindex_in_compressed_fortransposed,collcom,smat%nfvctr,moduloarray,ind_min,ind_max)

          !write(*,'(a,2i8)') 'after check_transposed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time0=real(tr1-tr0,kind=8)    


          ! Now check the compress_distributed layout
          !call check_compress_distributed_layout()
          call check_compress_distributed_layout(smat,ind_min,ind_max)

          !write(*,'(a,2i8)') 'after check_compress_distributed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time1=real(tr0-tr1,kind=8)        

          ! Now check the matrix matrix multiplications layout
          call check_matmul_layout(smat%smmm%nseq,smat%smmm%indices_extract_sequential,ind_min,ind_max)
          !write(*,'(a,2i8)') 'after check_matmul_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time2=real(tr1-tr0,kind=8)    
    
          ! Now check the sumrho operations
          !call check_sumrho_layout()
          call check_sumrho_layout(collcom_sr,smat%nfvctr,moduloarray,smat%matrixindex_in_compressed_fortransposed,ind_min,ind_max)

          call f_free_ptr(moduloarray)
          !write(*,'(a,2i8)') 'after check_sumrho_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time3=real(tr0-tr1,kind=8)    
    
          ! Now check the pseudo-exact orthonormalization during the input guess
          !call check_ortho_inguess()
          call check_ortho_inguess(smat,ind_min,ind_max)
          !write(*,'(a,2i8)') 'after check_ortho_inguess: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time4=real(tr1-tr0,kind=8)        

          ! Now check the submatrix extraction for the projector charge analysis
          call check_projector_charge_analysis(iproc, nproc, nat, smat, ind_min, ind_max)

          !!write(*,'(a,3i8)') 'after check_local_matrix_extents: iproc, ind_min, ind_max', iproc, ind_min, ind_max

          ! Get the global indices of ind_min and ind_max
          do i=1,2
              if (i==1) then
                  ii_ref = ind_min
              else
                  ii_ref = ind_max
              end if
              ! Search the indices iorb,jorb corresponding to ii_ref
              found=.false.

              ! not sure if OpenMP is really worth it here
              !$omp parallel default(none) &
              !$omp private(iseg,ii,iorb,jorb) &
              !$omp shared(smat,ii_ref,irow,icol,found,i)
              !$omp do
              outloop: do iseg=1,smat%nseg
                  if (.not. found) then
                     iorb = smat%keyg(1,2,iseg)
                     do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                         ii = matrixindex_in_compressed(smat, jorb, iorb)
                         !if (iproc==0) write(*,'(a,5i9)') 'i, ii_ref, ii, iorb, jorb', i, ii_ref, ii, iorb, jorb
                         if (ii==ii_ref) then
                             irow(i) = jorb
                             icol(i) = iorb
                             !exit outloop
                             !SM: I think one should do this within a critical section since it is shared, just to be sure...
                             !$omp critical
                             found=.true.
                             !$omp end critical
                         end if
                     end do
                  end if
              end do outloop
              !$omp end do
              !$omp end parallel

          end do
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time5=real(tr0-tr1,kind=8)    
          if (extra_timing) call cpu_time(trt1)  
          if (extra_timing) ttime=real(trt1-trt0,kind=8)  

          if (extra_timing.and.iproc==0) print*,'matextent',time0,time1,time2,time3,time4,time5,&
               time0+time1+time2+time3+time4+time5,ttime

          call timing(iproc,'matrix_extents','OF')    
    
!!$          contains
    
!!$            subroutine check_transposed_layout()
!!$              implicit none
!!$              integer :: ipt, ii, i0, i, i0i, iiorb, j, i0j, jjorb, ind, iorb, jorb
!!$              integer,dimension(:),pointer :: moduloarray
!!$              integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
!!$
!!$              call get_modulo_array(smat, moduloarray)
!!$
!!$              !SM: when the pointer within the type is used directly the code crashes with Intel, I have no idea why...
!!$              matrixindex_in_compressed_fortransposed => smat%matrixindex_in_compressed_fortransposed
!!$
!!$              !$omp parallel default(none) &
!!$              !$omp private(ipt,ii,i0,i,i0i,iiorb,iorb,j,i0j,jjorb,jorb,ind) &
!!$              !$omp shared(collcom,moduloarray,smat,ind_min,ind_max,matrixindex_in_compressed_fortransposed)
!!$              !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$              do ipt=1,collcom%nptsp_c
!!$                  ii=collcom%norb_per_gridpoint_c(ipt)
!!$                  i0 = collcom%isptsp_c(ipt)
!!$                  do i=1,ii
!!$                      i0i=i0+i
!!$                      iiorb=collcom%indexrecvorbital_c(i0i)
!!$                      iorb=moduloarray(iiorb)
!!$                      do j=1,ii
!!$                          i0j=i0+j
!!$                          jjorb=collcom%indexrecvorbital_c(i0j)
!!$                          jorb=moduloarray(jjorb)
!!$                          !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                          ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                          ind_min = min(ind_min,ind)
!!$                          ind_max = max(ind_max,ind)
!!$                      end do
!!$                  end do
!!$              end do
!!$              !$omp end do
!!$
!!$              !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$              do ipt=1,collcom%nptsp_f
!!$                  ii=collcom%norb_per_gridpoint_f(ipt)
!!$                  i0 = collcom%isptsp_f(ipt)
!!$                  do i=1,ii
!!$                      i0i=i0+i
!!$                      iiorb=collcom%indexrecvorbital_f(i0i)
!!$                      iorb=moduloarray(iiorb)
!!$                      do j=1,ii
!!$                          i0j=i0+j
!!$                          jjorb=collcom%indexrecvorbital_f(i0j)
!!$                          jorb=moduloarray(jjorb)
!!$                          !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                          ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
!!$                          ind_min = min(ind_min,ind)
!!$                          ind_max = max(ind_max,ind)
!!$                      end do
!!$                  end do
!!$              end do
!!$              !$omp end do
!!$              !$omp end parallel
!!$
!!$
!!$              call f_free_ptr(moduloarray)
!!$    
!!$            end subroutine check_transposed_layout


            !function get_transposed_index(jorb,iorb) result(ind)
            !    integer,intent(in) :: jorb, iorb
            !    integer :: ind
            !    integer :: jjorb,iiorb
            !    ! If iorb is smaller than the offset, add a periodic shift
            !    if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
            !        iiorb = iorb + smat%nfvctr
            !    else
            !        iiorb = iorb
            !    end if
            !    if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
            !        jjorb = jorb + smat%nfvctr
            !    else
            !        jjorb = jorb
            !    end if
            !    ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
            !end function get_transposed_index
    
    
!!$            subroutine check_compress_distributed_layout()
!!$              implicit none
!!$              integer :: i, nfvctrp, isfvctr, isegstart, isegend, iseg, ii, jorb
!!$              do i=1,2
!!$                  if (i==1) then
!!$                      nfvctrp = smat%nfvctrp
!!$                      isfvctr = smat%isfvctr
!!$                  else if (i==2) then
!!$                      nfvctrp = smat%smmm%nfvctrp
!!$                      isfvctr = smat%smmm%isfvctr
!!$                  end if
!!$                  if (nfvctrp>0) then
!!$                      isegstart=smat%istsegline(isfvctr+1)
!!$                      isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
!!$                      !$omp parallel default(none) &
!!$                      !$omp private(iseg,ii,jorb) &
!!$                      !$omp shared(isegstart,isegend,smat,ind_min,ind_max)
!!$                      !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$                      do iseg=isegstart,isegend
!!$                          ii=smat%keyv(iseg)-1
!!$                          ! A segment is always on one line, therefore no double loop
!!$                          do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$                              ii=ii+1
!!$                              ind_min = min(ii,ind_min)
!!$                              ind_max = max(ii,ind_max)
!!$                          end do
!!$                      end do
!!$                      !$omp end do
!!$                      !$omp end parallel
!!$                  end if
!!$              end do
!!$            end subroutine check_compress_distributed_layout
    
    
!!$            subroutine check_matmul_layout()
!!$              implicit none
!!$              integer :: iseq, ind
!!$              do iseq=1,smat%smmm%nseq
!!$                  ind=smat%smmm%indices_extract_sequential(iseq)
!!$                  ind_min = min(ind_min,ind)
!!$                  ind_max = max(ind_max,ind)
!!$              end do
!!$              !!write(*,'(a,3i8)') 'after check_matmul_layout: iproc, ind_min, ind_max', iproc, ind_min, ind_max
!!$            end subroutine check_matmul_layout
    
!!$            subroutine check_sumrho_layout()
!!$              implicit none
!!$              integer :: ipt, ii, i0, i, iiorb, ind, iorb
!!$              integer,dimension(:),pointer :: moduloarray
!!$
!!$              call get_modulo_array(smat, moduloarray)
!!$
!!$              !$omp parallel default(none) &
!!$              !$omp private(ipt,ii,i0,iiorb,iorb,ind) &
!!$              !$omp shared(collcom_sr,moduloarray,smat,ind_min,ind_max)
!!$              !$omp do reduction(min: ind_min) reduction(max: ind_max)
!!$              do ipt=1,collcom_sr%nptsp_c
!!$                  ii=collcom_sr%norb_per_gridpoint_c(ipt)
!!$                  i0=collcom_sr%isptsp_c(ipt)
!!$                  do i=1,ii
!!$                      iiorb=collcom_sr%indexrecvorbital_c(i0+i)
!!$                      iorb=moduloarray(iiorb)
!!$                      !ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
!!$                      ind=smat%matrixindex_in_compressed_fortransposed(iorb,iorb)
!!$                      !ind=get_transposed_index(smat,iiorb,iiorb)
!!$                      ind_min = min(ind_min,ind)
!!$                      ind_max = max(ind_max,ind)
!!$                  end do
!!$              end do
!!$              !$omp end do
!!$              !$omp end parallel
!!$
!!$              call f_free_ptr(moduloarray)
!!$
!!$              !contains
!!$
!!$              !  function get_transposed_index(jorb,iorb) res(ind)
!!$              !      integer,intent(in) :: jorb, iorb
!!$              !      integer :: ind
!!$              !      integer :: jjorb,iiorb
!!$              !      ! If iorb is smaller than the offset, add a periodic shift
!!$              !      if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
!!$              !          iiorb = iorb + smat%nfvctr
!!$              !      else
!!$              !          iiorb = iorb
!!$              !      end if
!!$              !      if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
!!$              !          jjorb = jorb + smat%nfvctr
!!$              !      else
!!$              !          jjorb = jorb
!!$              !      end if
!!$              !      ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
!!$              !  end function get_transposed_index
!!$            end subroutine check_sumrho_layout
    
    
          !!  function get_start_of_segment(smat, iiseg) result(ist)
    
          !!      do iseg=smat%nseg,1,-1
          !!          if (iiseg>=smat%keyv(iseg)) then
          !!              it = smat%keyv(iseg)
          !!              exit
          !!          end if
          !!      end do
    
          !!  end function get_start_of_segment
    
    
!!$          subroutine check_ortho_inguess()
!!$            integer :: iorb, iiorb, isegstart, isegend, iseg, j, i, jorb, korb, ind
!!$            logical,dimension(:),allocatable :: in_neighborhood
!!$    
!!$            in_neighborhood = f_malloc(smat%nfvctr,id='in_neighborhood')
!!$            
!!$            do iorb=1,smat%nfvctrp
!!$    
!!$                iiorb = smat%isfvctr + iorb
!!$                isegstart = smat%istsegline(iiorb)
!!$                isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
!!$                in_neighborhood = .false.
!!$
!!$                !$omp parallel default(none) &
!!$                !$omp private(iseg,j,i) &
!!$                !$omp private(jorb,korb,ind) &
!!$                !$omp shared(isegstart,isegend,smat,in_neighborhood) &
!!$                !$omp shared(ind_min,ind_max)
!!$
!!$                !$omp do
!!$                do iseg=isegstart,isegend
!!$                    ! A segment is always on one line, therefore no double loop
!!$                    j = smat%keyg(1,2,iseg)
!!$                    do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$                        in_neighborhood(i) = .true.
!!$                    end do
!!$                end do
!!$                !$omp end do
!!$    
!!$                !$omp do reduction(min: ind_min) reduction(max: ind_max) schedule(static,10) 
!!$                do jorb=1,smat%nfvctr
!!$                    if (.not.in_neighborhood(jorb)) cycle
!!$                    do korb=1,smat%nfvctr
!!$                        if (.not.in_neighborhood(korb)) cycle
!!$                        ind = matrixindex_in_compressed(smat,korb,jorb)
!!$                        if (ind>0) then
!!$                            ind_min = min(ind_min,ind)
!!$                            ind_max = max(ind_max,ind)
!!$                        end if
!!$                    end do
!!$                end do
!!$                !$omp end do
!!$
!!$                !$omp end parallel
!!$    
!!$            end do
!!$    
!!$            call f_free(in_neighborhood)
!!$    
!!$            !!do iorb=1,smat%nfvctrp
!!$            !!    iiorb = smat%isfvctr + iorb
!!$            !!    isegstart = smat%istsegline(iiorb)
!!$            !!    isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
!!$            !!    do iseg=isegstart,isegend
!!$            !!        ! A segment is always on one line, therefore no double loop
!!$            !!        j = smat%keyg(1,2,iseg)
!!$            !!        do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
!!$            !!            ind = matrixindex_in_compressed(smat,i,j)
!!$            !!            ind_min = min(ind_min,ind)
!!$            !!            ind_max = max(ind_max,ind)
!!$            !!        end do
!!$            !!    end do
!!$            !!end do
!!$    
!!$          end subroutine check_ortho_inguess


          !!!> Copied from projector_for_charge_analysis and extract_matrix
          !!subroutine check_projector_charge_analysis()
          !!  implicit none

          !!  integer :: ii, natp, jj, isat, kat, iatold, kkat, i, iat, j, ind
          !!  logical,dimension(:),allocatable :: neighbor

          !!  ! Parallelization over the number of atoms
          !!  ii = at%astruct%nat/bigdft_mpi%nproc
          !!  natp = ii
          !!  jj = at%astruct%nat - bigdft_mpi%nproc*natp
          !!  if (bigdft_mpi%iproc<jj) then
          !!      natp = natp + 1
          !!  end if
          !!  isat = (bigdft_mpi%iproc)*ii + min(bigdft_mpi%iproc,jj)


          !!  neighbor = f_malloc((/smat%nfvctr,natp/),id='neighbor')
          !!  do kat=1,natp
          !!      ! Determine the "neighbors"
          !!      iatold = 0
          !!      neighbor(:) = .false.
          !!      kkat = kat + isat
          !!      do i=1,smat%nfvctr
          !!           iat = smat%on_which_atom(i)
          !!           ! Only do the following for the first TMB per atom
          !!           if (iat==iatold) cycle
          !!           iatold = iat
          !!           if (iat==kkat) then
          !!               do j=1,smat%nfvctr
          !!                   ind =  matrixindex_in_compressed(smat, j, i)
          !!                   if (ind/=0) then
          !!                      neighbor(j) = .true.
          !!                   end if
          !!               end do
          !!           end if
          !!      end do

          !!      ! Determine the size of the matrix needed
          !!      do i=1,smat%nfvctr
          !!          if (neighbor(i)) then
          !!              do j=1,smat%nfvctr
          !!                  if (neighbor(j)) then
          !!                      ind =  matrixindex_in_compressed(smat, j, i)
          !!                      if (ind>0) then
          !!                          ind_min = min(ind_min,ind)
          !!                          ind_max = max(ind_max,ind)
          !!                      end if
          !!                  end if
          !!              end do
          !!          end if
          !!      end do

          !!  end do

          !!  call f_free(neighbor)

          !!end subroutine check_projector_charge_analysis

    end subroutine check_local_matrix_extents

    subroutine check_matmul_layout(nseq,indices_extract_sequential,ind_min,ind_max)
      implicit none
      integer, intent(in) :: nseq
      integer, dimension(nseq), intent(in) :: indices_extract_sequential
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: iseq,ind
      call f_routine(id='check_matmul_layout')

      !$omp parallel default(none) shared(nseq,indices_extract_sequential, ind_min, ind_max) private(iseq, ind)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do iseq=1,nseq
         ind=indices_extract_sequential(iseq)
         ind_min = min(ind_min,ind)
         ind_max = max(ind_max,ind)
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine check_matmul_layout



    !> Uses the CCS sparsity pattern to create a BigDFT sparse_matrix type
    subroutine ccs_to_sparsebigdft(iproc, nproc, nat, ncol, ncolp, iscol, nnonzero, &
               on_which_atom, row_ind, col_ptr, smat)
      use communications_base, only: comms_linear, comms_linear_null
      implicit none
      integer,intent(in) :: iproc, nproc, nat, ncol, ncolp, iscol, nnonzero
      integer,dimension(ncol),intent(in) :: on_which_atom
      !logical,intent(in) :: store_index
      integer,dimension(nnonzero),intent(in) :: row_ind
      integer,dimension(ncol),intent(in) :: col_ptr
      type(sparse_matrix),intent(out) :: smat

      ! Local variables
      integer :: icol, irow, i, ii
      integer,dimension(:,:),allocatable :: nonzero
      logical,dimension(:,:),allocatable :: mat
      type(comms_linear) :: collcom_dummy

      stop 'must be reworked'

      ! Calculate the values of nonzero and nonzero_mult which are required for
      ! the init_sparse_matrix routine.
      ! For the moment simple and stupid using a workarray of dimension ncol x ncol
      nonzero = f_malloc((/2,nnonzero/),id='nonzero')
      mat = f_malloc((/ncol,ncol/),id='mat')
      mat = .false.
      icol=1
      do i=1,nnonzero
          irow=row_ind(i)
          if (icol<ncol) then
              if (i>=col_ptr(icol+1)) then
                  icol=icol+1
              end if
          end if
          mat(irow,icol) = .true.
      end do
      ii = 0
      do irow=1,ncol
          write(333,*) col_ptr(irow)
          do icol=1,ncol
              if (mat(irow,icol)) then
                  ii = ii + 1
                  nonzero(2,ii) = irow
                  nonzero(1,ii) = icol
              end if
          end do
      end do

      call f_free(mat)

      call init_sparse_matrix(iproc, nproc, 1, 'F', ncol, ncolp, iscol, .false., &
           on_which_atom, nnonzero, nonzero, nnonzero, nonzero, smat)

      collcom_dummy = comms_linear_null()
      ! since no taskgroups are used, the values of iirow and iicol are just set to
      ! the minimum and maximum, respectively.
      call init_matrix_taskgroups(iproc, nproc, nat, .false., collcom_dummy, collcom_dummy, smat, &
           (/1,ncol/), (/1,ncol/))

      call f_free(nonzero)

    end subroutine ccs_to_sparsebigdft


    !> Uses the BigDFT sparsity pattern to create a BigDFT sparse_matrix type
    subroutine bigdft_to_sparsebigdft(iproc, nproc, nat, nspin, geocode, ncol, ncolp, iscol, &
               on_which_atom, nvctr, nseg, keyg, smat)
      use communications_base, only: comms_linear, comms_linear_null
      implicit none
      integer,intent(in) :: iproc, nproc, nat, nspin, ncol, ncolp, iscol, nvctr, nseg
      character(len=1),intent(in) :: geocode
      integer,dimension(ncol),intent(in) :: on_which_atom
      !logical,intent(in) :: store_index
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(out) :: smat

      ! Local variables
      integer :: icol, irow, i, ii, iseg
      !integer :: ncolpx
      integer,dimension(:,:),allocatable :: nonzero
      logical,dimension(:,:),allocatable :: mat
      !real(kind=8) :: tt
      type(comms_linear) :: collcom_dummy


      ! Calculate the values of nonzero and nonzero_mult which are required for
      ! the init_sparse_matrix routine.
      ! For the moment simple and stupid using a workarray of dimension ncol x ncol
      nonzero = f_malloc((/2,nvctr/),id='nonzero')
      mat = f_malloc((/ncol,ncol/),id='mat')
      mat = .false.

      do iseg=1,nseg
          do i=keyg(1,1,iseg),keyg(2,1,iseg)
              mat(keyg(1,2,iseg),i) = .true.
          end do
      end do
      ii = 0
      do irow=1,ncol
          do icol=1,ncol
              if (mat(irow,icol)) then
                  ii = ii + 1
                  nonzero(2,ii) = irow
                  nonzero(1,ii) = icol
              end if
          end do
      end do

      call f_free(mat)

      !!! Determine the number of columns per process
      !!tt = real(ncol,kind=8)/real(nproc,kind=8)
      !!ncolpx = floor(tt)
      !!ii = ncol - nproc*ncolpx
      !!if (iproc<ii) then
      !!    ncolp = ncolpx + 1
      !!else
      !!    ncolp = ncolpx
      !!end if
      !!
      !!! Determine the first column of each process
      !!i = 0
      !!do jproc=0,nproc-1
      !!    if (iproc==jproc) isorb = 1
      !!    if (jproc<ii) then
      !!        i = i + ncolpx + 1
      !!    else
      !!        i = i + ncolpx
      !!    end if
      !!end do

      call init_sparse_matrix(iproc, nproc, nspin, geocode, ncol, ncolp, iscol, .false., &
           on_which_atom, nvctr, nonzero, nvctr, nonzero, smat)

      collcom_dummy = comms_linear_null()
      ! since no taskgroups are used, the values of iirow and iicol are just set to
      ! the minimum and maximum, respectively.
      call init_matrix_taskgroups(iproc, nproc, nat, .false., collcom_dummy, collcom_dummy, smat, &
           (/1,ncol/), (/1,ncol/))

      call f_free(nonzero)

    end subroutine bigdft_to_sparsebigdft



    !> Assign the values of a sparse matrix in CCS format to a sparse matrix in the BigDFT format
    subroutine ccs_values_to_bigdft(ncol, nnonzero, row_ind, col_ptr, smat, val, mat)
      implicit none
      integer,intent(in) :: ncol, nnonzero
      integer,dimension(nnonzero),intent(in) :: row_ind
      integer,dimension(ncol),intent(in) :: col_ptr
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(nnonzero),intent(in) :: val
      type(matrices),intent(inout) :: mat

      ! Local variables
      integer :: icol, irow, i, ii
      logical,dimension(:,:),allocatable :: matg


      ! Calculate the values of nonzero and nonzero_mult which are required for
      ! the init_sparse_matrix routine.
      ! For the moment simple and stupid using a workarray of dimension ncol x ncol
      matg = f_malloc((/ncol,ncol/),id='matg')
      matg = .false.
      icol=1
      do i=1,nnonzero
          irow=row_ind(i)
          if (icol<ncol) then
              if (i>=col_ptr(icol+1)) then
                  icol=icol+1
              end if
          end if
          matg(irow,icol) = .true.
      end do
      ii = 0
      do irow=1,ncol
          do icol=1,ncol
              if (matg(irow,icol)) then
                  ii = ii + 1
                  mat%matrix_compr(ii) = val(ii)
              end if
          end do
      end do

      call f_free(matg)

    end subroutine ccs_values_to_bigdft


    subroutine read_ccs_format(filename, ncol, nnonzero, col_ptr, row_ind, val)
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: ncol, nnonzero
      integer,dimension(:),pointer,intent(out) :: col_ptr, row_ind
      real(kind=8),dimension(:),pointer,intent(out) :: val

      ! Local variables
      integer :: i, ii
      logical :: file_exists
      integer,parameter :: iunit=123

      inquire(file=filename,exist=file_exists)
      if (file_exists) then
          open(unit=iunit,file=filename)
          read(iunit,*) ncol, ii, nnonzero
          col_ptr = f_malloc_ptr(ncol,id='col_ptr')
          row_ind = f_malloc_ptr(nnonzero,id='row_ind')
          val = f_malloc_ptr(nnonzero,id='val')
          read(iunit,*) (col_ptr(i), i=1,ncol)
          read(iunit,*) (row_ind(i), i=1,nnonzero)
          do i=1,nnonzero
              read(iunit,*) val(i)
          end do
      else
          stop 'file not present'
      end if
      close(iunit)
    end subroutine read_ccs_format







    subroutine distribute_columns_on_processes_simple(iproc, nproc, ncol, ncolp, iscol)
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ncol
      integer,intent(out) :: ncolp, iscol
    
      ! Local variables
      integer :: ncolpx, ii, i, jproc
      real(kind=8) :: tt
    
      ! Determine the number of columns per process
      tt = real(ncol,kind=8)/real(nproc,kind=8)
      ncolpx = floor(tt)
      ii = ncol - nproc*ncolpx
      if (iproc<ii) then
          ncolp = ncolpx + 1
      else
          ncolp = ncolpx
      end if
      
      ! Determine the first column of each process
      i = 0
      do jproc=0,nproc-1
          if (iproc==jproc) iscol = i
          if (jproc<ii) then
              i = i + ncolpx + 1
          else
              i = i + ncolpx
          end if
      end do
    end subroutine distribute_columns_on_processes_simple


    !> Given the array workload which indicates the workload on each MPI task for a
    !! given distribution of the orbitals (or a similar quantity), this subroutine
    !! redistributes the orbitals such that the load unbalancing is optimal
    subroutine redistribute(nproc, norb, workload, workload_ideal, norb_par)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: nproc, norb
      real(kind=8),dimension(norb),intent(in) :: workload
      real(kind=8),intent(in) :: workload_ideal
      integer,dimension(0:nproc-1),intent(out) :: norb_par
    
      ! Local variables
      real(kind=8) :: tcount, jcount, wli, ratio, ratio_old, average
      real(kind=8),dimension(:),allocatable :: workload_par
      integer,dimension(:),allocatable :: norb_par_trial
      integer :: jproc, jjorb, jjorbtot, jorb, ii, imin, imax
    
      call f_routine(id='redistribute')
    
      wli = workload_ideal
    
      call f_zero(norb_par)
      if (norb>nproc) then
          workload_par = f_malloc(0.to.nproc-1,id='workload_par')
          norb_par_trial = f_malloc(0.to.nproc-1,id='norbpar_par_trial')
          tcount = 0.d0
          jcount = 0.d0
          jproc = 0
          jjorb = 0
          jjorbtot = 0
          do jorb=1,norb
              if (jproc==nproc-1) exit
              jjorb = jjorb + 1
              if(jorb==norb) exit !just to be sure that no out of bound happens
              tcount = tcount + workload(jorb)
              jcount = jcount + workload(jorb)
              if (abs(tcount-wli*real(jproc+1,kind=8)) <= &
                      abs(tcount+workload(jorb+1)-wli*real(jproc+1,kind=8))) then
              !!if (abs(tcount-workload_ideal*real(jproc+1,kind=8)) <= &
              !!        abs(tcount+workload(jorb+1)-workload_ideal*real(jproc+1,kind=8))) then
              !!if (tcount-workload_ideal*real(jproc+1,kind=8)<0.d0 .and. &
              !!        tcount+workload(jorb+1)-workload_ideal*real(jproc+1,kind=8)>0.d0) then
                  norb_par(jproc) = jjorb
                  workload_par(jproc) = jcount
                  jjorbtot = jjorbtot + jjorb
                  !if (bigdft_mpi%iproc==0) write(*,'(a,2i6,2es14.6)') 'jproc, jjorb, tcount/(jproc+1), wli', jproc, jjorb, tcount/(jproc+1), wli
                  jcount = 0.d0
                  jjorb = 0
                  jproc = jproc + 1
                  wli = get_dynamic_ideal_workload(nproc,jproc, tcount, workload_ideal)
              end if
          end do
          norb_par(nproc-1) = jjorb + (norb - jjorbtot) !take the rest
          workload_par(nproc-1) = sum(workload) - tcount

          if (sum(norb_par)/=norb) then
              call f_err_throw('wrong first partition of the workload; sum of distributed workload is '&
                   &//trim(yaml_toa(sum(norb_par)))//', but should be '//trim(yaml_toa(norb)),&
                   err_name='BIGDFT_RUNTIME_ERROR')
          end if

          !do jproc=0,nproc-1
          !    if (iproc==0) write(*,*) 'jproc, norb_par(jproc)', jproc, norb_par(jproc)
          !end do
          !if (bigdft_mpi%iproc==0) write(*,'(a,2i6,2es14.6)') 'jproc, jjorb, tcount/(jproc+1), workload_ideal', &
          !        jproc, jjorb+(norb-jjorbtot), sum(workload)-tcount, workload_ideal
    
          ! Now take away one element from the maximum and add it to the minimum.
          ! Repeat this as long as the ratio max/average decreases 
          average = sum(workload_par)/real(nproc,kind=8)
          ratio_old = maxval(workload_par)/average
          adjust_loop: do
              imin = minloc(workload_par,1) - 1 !subtract 1 because the array starts a 0
              imax = maxloc(workload_par,1) - 1 !subtract 1 because the array starts a 0
              call vcopy(nproc, norb_par(0), 1, norb_par_trial(0), 1)
              norb_par_trial(imin) = norb_par(imin) + 1
              norb_par_trial(imax) = norb_par(imax) - 1
    
              call f_zero(workload_par)
              ii = 0
              do jproc=0,nproc-1
                  do jorb=1,norb_par_trial(jproc)
                      ii = ii + 1
                      workload_par(jproc) = workload_par(jproc) + workload(ii)
                  end do
              end do
              average = sum(workload_par)/real(nproc,kind=8)
              ratio = maxval(workload_par)/average
              !if (bigdft_mpi%iproc==0) write(*,*) 'ratio, ratio_old', ratio, ratio_old
              if (ratio<ratio_old) then
                  call vcopy(nproc, norb_par_trial(0), 1, norb_par(0), 1)
                  ratio_old = ratio
              else
                  exit adjust_loop
              end if
          end do adjust_loop
    
          call f_free(workload_par)
          call f_free(norb_par_trial)
      else
          ! Equal distribution
          norb_par(0:norb-1) = 1
      end if

      if (sum(norb_par)/=norb) then
          call f_err_throw('wrong second partition of the workload; sum of distributed workload is '&
               &//trim(yaml_toa(sum(norb_par)))//', but should be '//trim(yaml_toa(norb)),&
               err_name='BIGDFT_RUNTIME_ERROR')
      end if
    
      call f_release_routine()
    
    end subroutine redistribute

    ! Get dynamically a new ideal workload
    function get_dynamic_ideal_workload(nproc,jproc, wltot, wli) result(wl)
      implicit none
      integer,intent(in) :: nproc !<number of MPI tasks
      integer,intent(in) :: jproc !<currently handled task
      real(kind=8),intent(in) :: wltot !<total workload assigned so far
      real(kind=8),intent(in) :: wli !< theoretical ideal workload
      real(kind=8) :: wl !<new ideal workload
      real(kind=8) :: wls

      ! Average workload so far
      wls = wltot/real(jproc,kind=8)

      ! The new ideal workload is a weighted sum of the average workload so far
      ! and the theoretical ideal workload
      wl = (nproc-jproc)*wls + jproc*wli
      wl = wl/real(nproc,kind=8) 

    end function get_dynamic_ideal_workload




    !!function get_transposed_index(smat,jorb,iorb) result(ind)
    !!    implicit none
    !!    type(sparse_matrix),intent(in) :: smat
    !!    integer,intent(in) :: jorb, iorb
    !!    integer :: ind
    !!    integer :: jjorb,iiorb,ii,jj
    !!    ! If iorb,jorb is smaller than the offset, add a periodic shift
    !!    ! This is rather slow...
    !!    if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
    !!        iiorb = iorb + smat%nfvctr
    !!    else
    !!        iiorb = iorb
    !!    end if
    !!    if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
    !!        jjorb = jorb + smat%nfvctr
    !!    else
    !!        jjorb = jorb
    !!    end if

    !!    !!!! Hopefully faster
    !!    !!!! ii should be 1 if iorb<smat%offset_matrixindex_in_compressed_fortransposed and 0 otherwise...
    !!    !!!ii = iorb/smat%offset_matrixindex_in_compressed_fortransposed
    !!    !!!! Now ii should be 0 if iorb<smat%offset_matrixindex_in_compressed_fortransposed and >=1 otherwise
    !!    !!!ii = 1 - 1**ii
    !!    !!!! Now ii should be 0 if iorb<smat%offset_matrixindex_in_compressed_fortransposed and non-zero otherwise
    !!    !!!iiorb = iorb + ii*smat%nfvctr

    !!    ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
    !!end function get_transposed_index


    subroutine get_modulo_array(smat, moduloarray)
      implicit none
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(:),pointer :: moduloarray
      ! Local variables
      integer :: i
      moduloarray = f_malloc_ptr(smat%nfvctr,id='moduloarray')
      !$omp parallel default(none) &
      !$omp shared(moduloarray,smat) &
      !$omp private(i)
      !$omp do
      do i=1,smat%nfvctr
          moduloarray(i) = modulo(i-smat%offset_matrixindex_in_compressed_fortransposed,smat%nfvctr)+1
      end do
      !$omp end do
      !$omp end parallel
    end subroutine get_modulo_array


    !> Copied from projector_for_charge_analysis and extract_matrix
    subroutine check_projector_charge_analysis(iproc, nproc, nat, smat, ind_min, ind_max)
      use module_base, only: bigdft_mpi
      use sparsematrix_base, only: sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nat
      type(sparse_matrix),intent(in) :: smat
      integer,intent(inout) :: ind_min, ind_max

      integer :: ii, natp, jj, isat, kat, iatold, kkat, i, iat, j, ind
      logical,dimension(:),allocatable :: neighbor

      ! Parallelization over the number of atoms
      ii = nat/nproc
      natp = ii
      jj = nat - nproc*natp
      if (iproc<jj) then
          natp = natp + 1
      end if
      isat = (iproc)*ii + min(iproc,jj)


      neighbor = f_malloc(smat%nfvctr,id='neighbor')
      do kat=1,natp
          ! Determine the "neighbors"
          iatold = 0
          neighbor(:) = .false.
          kkat = kat + isat
          do i=1,smat%nfvctr
               iat = smat%on_which_atom(i)
               ! Only do the following for the first TMB per atom
               if (iat==iatold) cycle
               iatold = iat
               if (iat==kkat) then
                   do j=1,smat%nfvctr
                       ind =  matrixindex_in_compressed(smat, j, i)
                       if (ind/=0) then
                          neighbor(j) = .true.
                       end if
                   end do
               end if
          end do

          ! Determine the size of the matrix needed
          do i=1,smat%nfvctr
              if (neighbor(i)) then
                  do j=1,smat%nfvctr
                      if (neighbor(j)) then
                          ind =  matrixindex_in_compressed(smat, j, i)
                          if (ind>0) then
                              ind_min = min(ind_min,ind)
                              ind_max = max(ind_max,ind)
                          end if
                      end if
                  end do
              end if
          end do

      end do

      call f_free(neighbor)

    end subroutine check_projector_charge_analysis

    subroutine find_minmax_transposed(matrixindex_in_compressed_fortransposed,collcom,nfvctr,moduloarray,ind_min,ind_max)
      use communications_base, only: comms_linear
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom
      integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, i0i, iiorb, j, i0j, jjorb, ind, iorb, jorb

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,i,i0i,iiorb,iorb,j,i0j,jjorb,jorb,ind) &
      !$omp shared(collcom,moduloarray,ind_min,ind_max,matrixindex_in_compressed_fortransposed,nfvctr)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_c
         ii=collcom%norb_per_gridpoint_c(ipt)
         i0 = collcom%isptsp_c(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_c(i0i)
            iorb=moduloarray(iiorb)
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_c(i0j)
               jorb=moduloarray(jjorb)
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do

      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_f
         ii=collcom%norb_per_gridpoint_f(ipt)
         i0 = collcom%isptsp_f(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_f(i0i)
            iorb=moduloarray(iiorb)
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_f(i0j)
               jorb=moduloarray(jjorb)
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine find_minmax_transposed

    subroutine find_startendseg_transposed(ind_min,ind_max,smat)
      use sparsematrix_base, only: sparse_matrix
      implicit none
      integer, intent(in) :: ind_min,ind_max
      type(sparse_matrix),intent(inout) :: smat
      !local variables    
      logical :: found
      integer :: iiseg1, iiseg2, iorb, jorb,iseg

      ! Store these values
      smat%istartend_t(1) = ind_min
      smat%istartend_t(2) = ind_max

      ! Determine to which segments this corresponds
      iiseg1 = smat%nseg
      iiseg2 = 1
      !$omp parallel default(none) shared(smat, iiseg1, iiseg2) private(iseg, found)
      found = .false.
      !$omp do reduction(min: iiseg1)
      do iseg=1,smat%nseg
         ! A segment is always on one line
         if (.not.found) then
            if (smat%keyv(iseg)+smat%keyg(2,1,iseg)-smat%keyg(1,1,iseg)>=smat%istartend_t(1)) then
               !smat%istartendseg_t(1)=iseg
               iiseg1=iseg
               found = .true.
            end if
         end if
      end do
      !$omp end do
      found = .false.
      !$omp do reduction(max: iiseg2)
      do iseg=smat%nseg,1,-1
         if (.not.found) then
            if (smat%keyv(iseg)<=smat%istartend_t(2)) then
               !smat%istartendseg_t(2)=iseg
               iiseg2=iseg
               found = .true.
            end if
         end if
      end do
      !$omp end do
      !$omp end parallel
      smat%istartendseg_t(1) = iiseg1
      smat%istartendseg_t(2) = iiseg2
    end subroutine find_startendseg_transposed

    subroutine check_compress_distributed_layout(smat,ind_min,ind_max)
      implicit none
      type(sparse_matrix),intent(in) :: smat
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: i,nfvctrp,isfvctr,isegstart,isegend,iseg,jorb,ii
      
      !call f_routine(id='check_compress_distributed_layout')

      do i=1,2
         if (i==1) then
            nfvctrp = smat%nfvctrp
            isfvctr = smat%isfvctr
         else if (i==2) then
            nfvctrp = smat%smmm%nfvctrp
            isfvctr = smat%smmm%isfvctr
         end if
         if (nfvctrp>0) then
            isegstart=smat%istsegline(isfvctr+1)
            isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
            !$omp parallel default(none) &
            !$omp shared(isegstart, isegend, smat, ind_min, ind_max) &
            !$omp private(iseg, ii,jorb)
            !$omp do reduction(min: ind_min) reduction(max: ind_max)
            do iseg=isegstart,isegend
               ii=smat%keyv(iseg)-1
               ! A segment is always on one line, therefore no double loop
               do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                  ii=ii+1
                  ind_min = min(ii,ind_min)
                  ind_max = max(ii,ind_max)
               end do
            end do
            !$omp end do
            !$omp end parallel
         end if
      end do

      !call f_release_routine()

    end subroutine check_compress_distributed_layout

    subroutine check_sumrho_layout(collcom_sr,nfvctr,moduloarray,matrixindex_in_compressed_fortransposed,ind_min,ind_max)
      use communications_base, only: comms_linear
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom_sr
      integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, iiorb, ind, iorb

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,iiorb,iorb,ind,i) &
      !$omp shared(collcom_sr,moduloarray,matrixindex_in_compressed_fortransposed,ind_min,ind_max)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom_sr%nptsp_c
         ii=collcom_sr%norb_per_gridpoint_c(ipt)
         i0=collcom_sr%isptsp_c(ipt)
         do i=1,ii
            iiorb=collcom_sr%indexrecvorbital_c(i0+i)
            iorb=moduloarray(iiorb)
            !ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
            ind=matrixindex_in_compressed_fortransposed(iorb,iorb)
            !ind=get_transposed_index(smat,iiorb,iiorb)
            ind_min = min(ind_min,ind)
            ind_max = max(ind_max,ind)
         end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine check_sumrho_layout

    subroutine check_ortho_inguess(smat,ind_min,ind_max)
      implicit none
      type(sparse_matrix),intent(in) :: smat
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: iorb, iiorb, isegstart, isegend, iseg, j, i, jorb, korb, ind, nthread, ithread
      logical, dimension(:,:), allocatable :: in_neighborhood
      !$ integer :: omp_get_max_threads, omp_get_thread_num

      !call f_routine(id='check_ortho_inguess')

      ! Allocate the array for all threads to avoid that it has to be declared private
      nthread = 1
      !$ nthread = omp_get_max_threads()
      in_neighborhood = f_malloc((/1.to.smat%nfvctr,0.to.nthread-1/),id='in_neighborhood')

      ithread = 0
      !$omp parallel default(none) &
      !$omp shared(smat, in_neighborhood, ind_min, ind_max) &
      !$omp private(iorb, iiorb, isegstart, isegend, iseg, j, jorb, korb, ind,i) &
      !$omp firstprivate(ithread)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do iorb=1,smat%nfvctrp
         !$ ithread = omp_get_thread_num()

         iiorb = smat%isfvctr + iorb
         isegstart = smat%istsegline(iiorb)
         isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
         in_neighborhood(:,ithread) = .false.
         do iseg=isegstart,isegend
            ! A segment is always on one line, therefore no double loop
            j = smat%keyg(1,2,iseg)
            do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
               in_neighborhood(i,ithread) = .true.
            end do
         end do

         do jorb=1,smat%nfvctr
            if (.not.in_neighborhood(jorb,ithread)) cycle
            do korb=1,smat%nfvctr
               if (.not.in_neighborhood(korb,ithread)) cycle
               ind = matrixindex_in_compressed(smat,korb,jorb)
               if (ind>0) then
                  ind_min = min(ind_min,ind)
                  ind_max = max(ind_max,ind)
               end if
            end do
         end do

      end do
      !$omp end do
      !$omp end parallel

      call f_free(in_neighborhood)


      !call f_release_routine()

    end subroutine check_ortho_inguess



    !> Converts the sparse matrix descriptors from BigDFT to those from the CCS format.
    !! It required that each column has at least one non-zero element.
    subroutine sparsebigdft_to_ccs(nfvctr, nvctr, nseg, keyg, row_ind, col_ptr)
      use module_base
      implicit none
      ! Calling arguments
      integer,intent(in) :: nfvctr !< number of columns/rows
      integer,intent(in) :: nvctr !< number of non-zero elements
      integer,intent(in) :: nseg !< number of segments for the BigDFT sparsity pattern
      !integer,dimension(nseg),intent(in) :: keyv !< starting index of each segment within the compressed sparse array
      integer,dimension(2,2,nseg),intent(in) :: keyg !< starting ad ending "coordinates" of each segment within the uncompressed matrix
                                                     !! 1,1,iseg: starting line index of segment iseg
                                                     !! 2,1,iseg: ending line index of segment iseg
                                                     !! 1,2,iseg: starting column index of segment iseg
                                                     !! 2,2,iseg: ending column index of segment iseg
      integer,dimension(nvctr),intent(out) :: row_ind !< row index of each non-zero element
      integer,dimension(nfvctr),intent(out) :: col_ptr !< index of the first element of each column within the compressed sparse array
      ! Local variables
      integer :: ii, icol_old, iseg, icol, i

      call f_routine(id='sparsebigdft_to_ccs')

      ii = 1
      icol_old = 0
      do iseg=1,nseg
          icol = keyg(1,2,iseg) !column index
          if (icol>icol_old) then
              col_ptr(icol) = ii
              icol_old = icol
          end if
          do i=keyg(1,1,iseg),keyg(2,1,iseg)
              !i gives the row index
              row_ind(ii) = i
              ii = ii + 1
          end do
      end do

      call f_release_routine()

    end subroutine sparsebigdft_to_ccs


    !> Converts the sparse matrix descriptors from BigDFT to those from the CCS format.
    !! It required that each column has at least one non-zero element.
    subroutine ccs_to_sparsebigdft_short(nfvctr, nvctr, row_ind, col_ptr, nseg, keyv, keyg)
      use module_base
      implicit none
      ! Calling arguments
      integer,intent(in) :: nfvctr !< number of columns/rows
      integer,intent(in) :: nvctr !< number of non-zero elements
      integer,dimension(nvctr),intent(in) :: row_ind !< row index of each non-zero element
      integer,dimension(nfvctr),intent(in) :: col_ptr !< index of the first element of each column within the compressed sparse array
      integer,intent(out) :: nseg !< number of segments for the BigDFT sparsity pattern
      integer,dimension(:),pointer,intent(out) :: keyv !< starting index of each segment within the compressed sparse array
      integer,dimension(:,:,:),pointer,intent(out) :: keyg !< starting ad ending "coordinates" of each segment within the uncompressed matrix
                                                     !! 1,1,iseg: starting line index of segment iseg
                                                     !! 2,1,iseg: ending line index of segment iseg
                                                     !! 1,2,iseg: starting column index of segment iseg
                                                     !! 2,2,iseg: ending column index of segment iseg
      ! Local variables
      integer :: icol, is, ie, i, ii, ii_old, ivctr, iseg, n


      nseg = 0
      do icol=1,nfvctr !column index
          is = col_ptr(icol)
          if (icol<nfvctr) then
              ie = col_ptr(icol+1) - 1
          else
              ie = nvctr
          end if
          ii_old = -1
          do i=is,ie
              ii = row_ind(i)
              if (ii==ii_old+1) then
                  !in the same segment
              else
                  !new segment in the same column
                  nseg = nseg + 1
              end if
              write(*,*) 'icol, i, ii, nseg', icol, i, ii, nseg
              ii_old = ii
          end do
      end do

      keyv = f_malloc_ptr(nseg,id='keyv')
      keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')

      iseg = 0
      ivctr = 0
      do icol=1,nfvctr !column index
          is = col_ptr(icol)
          if (icol<nfvctr) then
              ie = col_ptr(icol+1) - 1
          else
              ie = nvctr
          end if
          ii_old = -1
          do i=is,ie
              ii = row_ind(i)
              ivctr = ivctr + 1
              if (ii==ii_old+1) then
                  !in the same segment
              else
                  !new segment in the same column
                  !! close previous segment in the same column
                  !if (iseg>=1) then
                  !    keyg(2,1,iseg) = ii_old
                  !    keyg(2,2,iseg) = icol
                  !end if
                  iseg = iseg + 1
                  keyv(iseg) =  ivctr
                  keyg(1,1,iseg) =  ii
                  keyg(1,2,iseg) =  icol
              end if
              ii_old = ii
          end do
      end do

      ! Close the segments
      do iseg=1,nseg
          if (iseg<nseg) then
              ii = keyv(iseg+1)
          else
              ii = nvctr + 1
          end if
          n = ii - keyv(iseg)
          keyg(2,1,iseg) = keyg(1,1,iseg) + n - 1
          keyg(2,2,iseg) = keyg(1,2,iseg)
      end do

    end subroutine ccs_to_sparsebigdft_short


    subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, collcom, collcom_shamop, &
               collcom_sr, sparsemat)
      use module_base
      use communications_base, only: comms_linear
      use sparsematrix_base, only: sparse_matrix
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
      type(sparse_matrix),intent(inout) :: sparsemat
      
      ! Local variables
      integer :: iorb, jorb, istat, imin, imax, nmiddle, imax_old, imin_old, iiorb, jjorb
      integer :: ii, imin_new, imax_new, i, nlen, j
      !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
      !integer :: compressed_index
    !  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
      character(len=*),parameter :: subname='init_sparse_matrix'
    
      call f_routine(id='init_matrixindex_in_compressed_fortransposed')
      call timing(iproc,'init_matrCompr','ON')
    
      ! for the calculation of overlaps and the charge density
      !imin=minval(collcom%indexrecvorbital_c)
      !imin=min(imin,minval(collcom%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
      !imax=maxval(collcom%indexrecvorbital_c)
      !imax=max(imax,maxval(collcom%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))
    
      nmiddle = sparsemat%nfvctr/2 + 1
    
      imin_old = huge(1)
      imax_old = 0
      imin_new = huge(1)
      imax_new = 0
      do i=1,size(collcom%indexrecvorbital_c)
          ii = mod(collcom%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom%indexrecvorbital_f)
          ii = mod(collcom%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_shamop%indexrecvorbital_c)
          ii = mod(collcom_shamop%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_shamop%indexrecvorbital_f)
          ii = mod(collcom_shamop%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_sr%indexrecvorbital_c)
          ii = mod(collcom_sr%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
    
    
      !!write(*,*) 'iproc, imin_old, imax_old', iproc, imin_old, imax_old
      !!write(*,*) 'iproc, imin_new, imax_new', iproc, imin_new, imax_new
    
      !! values regardless of the spin
      !imin=mod(imin-1,sparsemat%nfvctr)+1
      !imax=mod(imax-1,sparsemat%nfvctr)+1
    
    
      ! Determine with which size the array should be allocated
      if (imax_new-imin_new<0) then
          ! everything in either first or second half
          imin = imin_old
          imax = imax_old
          !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
      else
          ! in both half
          if (imax_old-imin_old>imax_new-imin_new) then
              ! wrap around
              imin = imin_new
              imax = imax_new
              !sparsemat%offset_matrixindex_in_compressed_fortransposed = imin_new
          else
              ! no wrap around
              imin = imin_old
              imax = imax_old
              !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
          end if
      end if
    
      !!! Check
      !!if (sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1) then
      !!    stop 'sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1'
      !!end if
    
      nlen = imax - imin + 1
      sparsemat%offset_matrixindex_in_compressed_fortransposed = imin
      !!write(*,*) 'iproc, imin, imax, nlen', iproc, imin, imax, nlen
    
      !!! This is a temporary solution for spin polarized systems
      !!imax=min(imax,orbs%norbu)
    
    
    
      !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
      !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
      !sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
      !    id='sparsemat%matrixindex_in_compressed_fortransposed')
      sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/nlen,nlen/),&
          id='sparsemat%matrixindex_in_compressed_fortransposed')
    
      !$omp parallel do default(private) shared(sparsemat,imin,imax)
      do iorb=imin,imax
          i = iorb - imin + 1
          do jorb=imin,imax
              j = jorb - imin + 1
              !@ii=(jorb-1)*sparsemat%nfvctr+iorb
              !@ispin=(ii-1)/sparsemat%nfvctr+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
              !@iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
              !@jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
              !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iiorb,jjorb,orbs%norbu,sparsemat)
              iiorb = mod(iorb-1,sparsemat%nfvctr)+1
              jjorb = mod(jorb-1,sparsemat%nfvctr)+1
              !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
              sparsemat%matrixindex_in_compressed_fortransposed(i,j)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
              !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
              !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
          end do
      end do
      !$omp end parallel do
    
      !@! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
      !@if (ispin==2) then
      !@    matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
      !@end if
    
      call timing(iproc,'init_matrCompr','OF')
      call f_release_routine()
    
    end subroutine init_matrixindex_in_compressed_fortransposed


    subroutine distribute_on_threads(nout, nthread, ise)
      use module_base
      implicit none

      ! Calling arguments
      integer,intent(in) :: nout
      integer,intent(out) :: nthread
      integer,dimension(:,:),pointer :: ise

      ! Local variables
      integer :: ii, jthread
      integer,dimension(:),allocatable :: n
      !$ integer :: omp_get_max_threads

      call f_routine(id='distribute_on_threads')

      ! OpenMP parallelization using a large workarray
      nthread = 1
      !$ nthread = omp_get_max_threads()

      ! Determine the number of iterations to be done by each thread
      n = f_malloc(0.to.nthread-1,id='n')
      ii = nout/nthread
      n(0:nthread-1) = ii
      ii = nout - nthread*ii
      n(0:ii-1) = n(0:ii-1) + 1
      ! Check
      if (sum(n)/=nout) call f_err_throw('sum(n)/=nout',err_name='BIGDFT_RUNTIME_ERROR')

      ! Determine the first and last iteration for each thread
      ise = f_malloc_ptr((/1.to.2,0.to.nthread-1/),id='ise')
      ise(1,0) = 1
      do jthread=1,nthread-1
          ise(1,jthread) = ise(1,jthread-1) + n(jthread-1)
          ise(2,jthread-1) = ise(1,jthread) -1
      end do
      ise(2,nthread-1) = nout
      ! Check
      ii = 0
      do jthread=0,nthread-1
          ii = ii + ise(2,jthread) - ise(1,jthread) + 1
          if (jthread>1) then
              if (ise(1,jthread)/=ise(2,jthread-1)+1) then
                  call f_err_throw('ise(1,jthread)/=ise(2,jthread-1)+1',err_name='BIGDFT_RUNTIME_ERROR')
              end if
          end if
      end do
      if (ii/=nout) call f_err_throw('ii/=nout',err_name='BIGDFT_RUNTIME_ERROR')

      call f_free(n)

      call f_release_routine()

    end subroutine distribute_on_threads

end module sparsematrix_init

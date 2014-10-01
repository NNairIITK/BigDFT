!> @file
!!  File defining the structures to deal with the sparse matrices
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
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
  !!public :: compressed_index
  public :: matrixindex_in_compressed
  public :: check_kernel_cutoff
  public :: init_matrix_taskgroups

contains

    !!!> Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
    !!function compressed_index(irow, jcol, norb, sparsemat)
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  integer,intent(in) :: irow, jcol, norb
    !!  type(sparse_matrix),intent(in) :: sparsemat
    !!  integer :: compressed_index
    !!
    !!  ! Local variables
    !!  integer :: ii, iseg
    !!
    !!  ii=(jcol-1)*norb+irow
    !!
    !!  iseg=sparsemat%istsegline(jcol)
    !!  do
    !!      if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
    !!          ! The matrix element is in this segment
    !!           compressed_index = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
    !!          return
    !!      end if
    !!      iseg=iseg+1
    !!      if (iseg>sparsemat%nseg) exit
    !!      if (ii<sparsemat%keyg(1,iseg)) then
    !!          compressed_index=0
    !!          return
    !!      end if
    !!  end do
    !!
    !!  ! Not found
    !!  compressed_index=0
    !!
    !!end function compressed_index


    integer function matrixindex_in_compressed(sparsemat, iorb, jorb, init_, n_)
      use sparsematrix_base, only: sparse_matrix
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: sparsemat
      integer,intent(in) :: iorb, jorb
      !> The optional arguments should only be used for initialization purposes
      !! if one is sure what one is doing. Might be removed later.
      logical,intent(in),optional :: init_
      integer,intent(in),optional :: n_
    
      ! Local variables
      integer :: ii, ispin, iiorb, jjorb
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
              ! there seems to be a mix up the spin matrices
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
          matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
      end if
    
    contains

      ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
      integer function compressed_index_fn(irow, jcol, norb, sparsemat)
        implicit none
      
        ! Calling arguments
        integer,intent(in) :: irow, jcol, norb
        type(sparse_matrix),intent(in) :: sparsemat
      
        ! Local variables
        integer :: ii, iseg
      
        ii=(jcol-1)*norb+irow
      
        iseg=sparsemat%istsegline(jcol)
        do
            if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
                ! The matrix element is in sparsemat segment
                 compressed_index_fn = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
                return
            end if
            iseg=iseg+1
            if (iseg>sparsemat%nseg) exit
            if (ii<sparsemat%keyg(1,iseg)) then
                compressed_index_fn=0
                return
            end if
        end do
      
        ! Not found
        compressed_index_fn=0
      
      end function compressed_index_fn
    end function matrixindex_in_compressed




    subroutine check_kernel_cutoff(iproc, orbs, atoms, lzd)
      use module_types
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc
      type(orbitals_data),intent(in) :: orbs
      type(atoms_data),intent(in) :: atoms
      type(local_zone_descriptors),intent(inout) :: lzd
    
      ! Local variables
      integer :: iorb, ilr, iat, iatype
      real(kind=8) :: cutoff_sf, cutoff_kernel
      character(len=20) :: atomname
      logical :: write_data
      logical,dimension(atoms%astruct%ntypes) :: write_atomtype
    
      write_atomtype=.true.
    
      if (iproc==0) then
          call yaml_sequence_open('check of kernel cutoff radius')
      end if
    
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
    
          ! cutoff radius of the support function, including shamop region
          cutoff_sf=lzd%llr(ilr)%locrad+8.d0*lzd%hgrids(1)
    
          ! cutoff of the density kernel
          cutoff_kernel=lzd%llr(ilr)%locrad_kernel
    
          ! check whether the date for this atomtype has already shoudl been written
          iat=orbs%onwhichatom(iorb)
          iatype=atoms%astruct%iatype(iat)
          if (write_atomtype(iatype)) then
              if (iproc==0) then
                  write_data=.true.
              else
                  write_data=.false.
              end if
              write_atomtype(iatype)=.false.
          else
              write_data=.false.
          end if
    
          ! Adjust if necessary
          if (write_data) then
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
              call yaml_map('atom type',atomname)
          end if
          if (cutoff_sf>cutoff_kernel) then
              if (write_data) then
                  call yaml_map('adjustment required',.true.)
                  call yaml_map('new value',cutoff_sf,fmt='(f6.2)')
              end if
              lzd%llr(ilr)%locrad_kernel=cutoff_sf
          else
              if (write_data) then
                  call yaml_map('adjustment required',.false.)
              end if
          end if
          if (write_data) then
              call yaml_mapping_close()
          end if
      end do
    
      if (iproc==0) then
          call yaml_sequence_close
      end if
    
    
    end subroutine check_kernel_cutoff


    subroutine init_sparse_matrix_matrix_multiplication(iproc, nproc, norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat)
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(inout) :: sparsemat

      integer :: ierr, jproc, iorb, jjproc, iiorb, nseq_min, nseq_max, iseq, ind, ii, iseg, ncount
      integer,dimension(:),allocatable :: nseq_per_line, norb_par_ideal, isorb_par_ideal
      integer,dimension(:,:),allocatable :: istartend_dj, istartend_mm
      integer,dimension(:,:),allocatable :: temparr
      real(kind=8) :: rseq, rseq_ideal, tt, ratio_before, ratio_after
      logical :: printable

      ! Calculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
      ! the default partitioning of the matrix columns.
      call get_nout(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)
      nseq_per_line = f_malloc0(norb,id='nseq_per_line')
      call determine_sequential_length(norb, norbp, isorb, nseg, &
           nsegline, istsegline, keyg, sparsemat, &
           sparsemat%smmm%nseq, nseq_per_line)
      if (nproc>1) call mpiallred(nseq_per_line(1), norb, mpi_sum, bigdft_mpi%mpi_comm)
      rseq=real(sparsemat%smmm%nseq,kind=8) !real to prevent integer overflow
      if (nproc>1) call mpiallred(rseq, 1, mpi_sum, bigdft_mpi%mpi_comm)


      norb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
      isorb_par_ideal = f_malloc(0.to.nproc-1,id='norb_par_ideal')
      ! Assign the columns of the matrix to the processes such that the load
      ! balancing is optimal
      ! First the default initializations
      norb_par_ideal(:)=0
      isorb_par_ideal(:)=norb
      rseq_ideal = rseq/real(nproc,kind=8)
      jjproc=0
      tt=0.d0
      iiorb=0
      isorb_par_ideal(0)=0
      do iorb=1,norb
          iiorb=iiorb+1
          tt=tt+real(nseq_per_line(iorb),kind=8)
          if (tt>=real(jjproc+1,kind=8)*rseq_ideal .and. jjproc/=nproc-1) then
              norb_par_ideal(jjproc)=iiorb
              isorb_par_ideal(jjproc+1)=iorb
              jjproc=jjproc+1
              iiorb=0
          end if
      end do
      norb_par_ideal(jjproc)=iiorb


      ! some checks
      if (sum(norb_par_ideal)/=norb) stop 'sum(norb_par_ideal)/=norb'
      if (isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb) stop 'isorb_par_ideal(nproc-1)+norb_par_ideal(nproc-1)/=norb'

      ! Copy the values
      sparsemat%smmm%nfvctrp=norb_par_ideal(iproc)
      sparsemat%smmm%isfvctr=isorb_par_ideal(iproc)

      ! Get the load balancing
      nseq_min = sparsemat%smmm%nseq
      if (nproc>1) call mpiallred(nseq_min, 1, mpi_min, bigdft_mpi%mpi_comm)
      nseq_max = sparsemat%smmm%nseq
      if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
      if (nseq_min>0) then
          ratio_before = real(nseq_max,kind=8)/real(nseq_min,kind=8)
          printable=.true.
      else
          printable=.false.
      end if


      ! Realculate the values of sparsemat%smmm%nout and sparsemat%smmm%nseq with
      ! the optimized partitioning of the matrix columns.
      call get_nout(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, nsegline, istsegline, keyg, sparsemat%smmm%nout)
      call determine_sequential_length(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
           nsegline, istsegline, keyg, sparsemat, &
           sparsemat%smmm%nseq, nseq_per_line)

      ! Get the load balancing
      nseq_min = sparsemat%smmm%nseq
      if (nproc>1) call mpiallred(nseq_min, 1, mpi_min, bigdft_mpi%mpi_comm)
      nseq_max = sparsemat%smmm%nseq
      if (nproc>1) call mpiallred(nseq_max, 1, mpi_max, bigdft_mpi%mpi_comm)
      ! Not necessary to set the pritable flag (if nseq_min was zero before it should be zero here as well)
      if (nseq_min>0) then
          ratio_after = real(nseq_max,kind=8)/real(nseq_min,kind=8)
          if (.not.printable) stop 'this should not happen'
      else
          if (printable) stop 'this should not happen'
      end if
      if (iproc==0) then
          if (printable) then
              call yaml_map('sparse matmul load balancing naive / optimized',(/ratio_before,ratio_after/),fmt='(f4.2)')
          else
              call yaml_map('sparse matmul load balancing naive / optimized','printing not meaningful (division by zero)')
          end if
      end if
      

      call f_free(nseq_per_line)

      call allocate_sparse_matrix_matrix_multiplication(nproc, norb, nseg, nsegline, istsegline, keyg, sparsemat%smmm)


      ! Calculate some auxiliary variables
      temparr = f_malloc0((/0.to.nproc-1,1.to.2/),id='isfvctr_par')
      temparr(iproc,1) = sparsemat%smmm%isfvctr
      temparr(iproc,2) = sparsemat%smmm%nfvctrp
      if (nproc>1) then
          call mpiallred(temparr(0,1), 2*nproc,  mpi_sum, bigdft_mpi%mpi_comm)
      end if
      call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, sparsemat%nseg, sparsemat%nvctr, &
           temparr(0,1), temparr(0,2), sparsemat%istsegline, sparsemat%keyv, &
           sparsemat%smmm%isvctr, sparsemat%smmm%nvctrp, sparsemat%smmm%isvctr_par, sparsemat%smmm%nvctr_par)
      call f_free(temparr)

      sparsemat%smmm%nseg=nseg
      call vcopy(norb, nsegline(1), 1, sparsemat%smmm%nsegline(1), 1)
      call vcopy(norb, istsegline(1), 1, sparsemat%smmm%istsegline(1), 1)
      call vcopy(2*nseg, keyg(1,1), 1, sparsemat%smmm%keyg(1,1), 1)
      call init_onedimindices_new(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
           nsegline, istsegline, keyg, &
           sparsemat, sparsemat%smmm%nout, sparsemat%smmm%onedimindices)
      call get_arrays_for_sequential_acces(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
           nsegline, istsegline, keyg, sparsemat, &
           sparsemat%smmm%nseq, sparsemat%smmm%ivectorindex)
      call init_sequential_acces_matrix(norb, norb_par_ideal(iproc), isorb_par_ideal(iproc), nseg, &
           nsegline, istsegline, keyg, sparsemat, sparsemat%smmm%nseq, &
           sparsemat%smmm%indices_extract_sequential)

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
      do iseg=1,sparsemat%nseg
          if (sparsemat%keyv(iseg)>=sparsemat%smmm%istartend_mm(1)) then
              sparsemat%smmm%istartendseg_mm(1)=iseg
              exit
          end if
      end do
      do iseg=sparsemat%nseg,1,-1
          if (sparsemat%keyv(iseg)<=sparsemat%smmm%istartend_mm(2)) then
              sparsemat%smmm%istartendseg_mm(2)=iseg
              exit
          end if
      end do

      istartend_mm = f_malloc0((/1.to.2,0.to.nproc-1/),id='istartend_mm')
      istartend_mm(1:2,iproc) = sparsemat%smmm%istartend_mm(1:2)
      if (nproc>1) then
          call mpiallred(istartend_mm(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm)
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
          ind = min(ind,sparsemat%nvctr+1) ! +1 since the length shoulb be 0
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


      call f_free(norb_par_ideal)
      call f_free(isorb_par_ideal)
      call f_free(istartend_mm)
      call f_free(istartend_dj)
    end subroutine init_sparse_matrix_matrix_multiplication


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
      integer,dimension(2,nseg),intent(out) :: keyg
      
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
                  keyg(1,iseg)=ijorb
              end if
              segment_started=.true.
          else
              if (segment_started) then
                  ! close the previous segment
                  keyg(2,iseg)=ijorb-1
              end if
              segment_started=.false.
          end if
      end do
      ! close the last segment on the line if necessary
      if (segment_started) then
          keyg(2,iseg)=iline*norb
      end if
    end subroutine keyg_per_line



    !> Currently assuming square matrices
    subroutine init_sparse_matrix(iproc, nproc, nspin, norb, norbp, isorb, norbu, norbup, isorbu, store_index, &
               nnonzero, nonzero, nnonzero_mult, nonzero_mult, sparsemat, &
               allocate_full_, print_info_)
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin, norb, norbp, isorb, norbu, norbup, isorbu, nnonzero, nnonzero_mult
      logical,intent(in) :: store_index
      integer,dimension(nnonzero),intent(in) :: nonzero
      integer,dimension(nnonzero_mult),intent(in) :: nonzero_mult
      type(sparse_matrix), intent(out) :: sparsemat
      logical,intent(in),optional :: allocate_full_, print_info_
      
      ! Local variables
      integer :: jproc, iorb, jorb, iiorb, iseg, segn, ind
      integer :: jst_line, jst_seg
      integer :: ist, ivctr
      logical,dimension(:),allocatable :: lut
      integer :: nseg_mult, nvctr_mult, ivctr_mult
      integer,dimension(:),allocatable :: nsegline_mult, istsegline_mult
      integer,dimension(:,:),allocatable :: keyg_mult
      logical :: allocate_full, print_info

      call timing(iproc,'init_matrCompr','ON')

      call set_value_from_optional()

      lut = f_malloc(norb,id='lut')
    
      sparsemat=sparse_matrix_null()
    
      sparsemat%nspin=nspin
      sparsemat%nfvctr=norbu
      sparsemat%nfvctrp=norbup
      sparsemat%isfvctr=isorbu
      sparsemat%nfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nfvctr_par')
      sparsemat%isfvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isfvctr_par')

      ! Same as isorb_par and norb_par
      call to_zero(nproc, sparsemat%nfvctr_par(0))
      call to_zero(nproc, sparsemat%isfvctr_par(0))
      do jproc=0,nproc-1
          if (iproc==jproc) then
              sparsemat%isfvctr_par(jproc)=isorbu
              sparsemat%nfvctr_par(jproc)=norbup
          end if
      end do
      if (nproc>1) then
          call mpiallred(sparsemat%isfvctr_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nfvctr_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
      end if

      call allocate_sparse_matrix_basic(store_index, norbu, nproc, sparsemat)
    

      sparsemat%nseg=0
      sparsemat%nvctr=0
      sparsemat%nsegline=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero, nonzero, iiorb)
          call nseg_perline(norbu, lut, sparsemat%nseg, sparsemat%nvctr, sparsemat%nsegline(iiorb))
      end do


      if (nproc>1) then
          call mpiallred(sparsemat%nvctr, 1, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nseg, 1, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(sparsemat%nsegline(1), sparsemat%nfvctr, mpi_sum, bigdft_mpi%mpi_comm)
      end if


      ist=1
      do jorb=1,sparsemat%nfvctr
          ! Starting segment for this line
          sparsemat%istsegline(jorb)=ist
          ist=ist+sparsemat%nsegline(jorb)
      end do

    
      if (iproc==0 .and. print_info) then
          call yaml_map('total elements',norbu**2)
          call yaml_map('non-zero elements',sparsemat%nvctr)
          call yaml_map('sparsity in %',1.d2*dble(norbu**2-sparsemat%nvctr)/dble(norbu**2),fmt='(f5.2)')
      end if
    
      call allocate_sparse_matrix_keys(store_index, sparsemat)
    


      ivctr=0
      sparsemat%keyg=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero, nonzero, iiorb)
          call keyg_per_line(norbu, sparsemat%nseg, iiorb, sparsemat%istsegline(iiorb), &
               lut, ivctr, sparsemat%keyg)
      end do
    
      ! check whether the number of elements agrees
      if (nproc>1) then
          call mpiallred(ivctr, 1, mpi_sum, bigdft_mpi%mpi_comm)
      end if
      if (ivctr/=sparsemat%nvctr) then
          write(*,'(a,2i8)') 'ERROR: ivctr/=sparsemat%nvctr', ivctr, sparsemat%nvctr
          stop
      end if
      if (nproc>1) then
          call mpiallred(sparsemat%keyg(1,1), 2*sparsemat%nseg, mpi_sum, bigdft_mpi%mpi_comm)
      end if


      ! start of the segments
      sparsemat%keyv(1)=1
      do iseg=2,sparsemat%nseg
          sparsemat%keyv(iseg) = sparsemat%keyv(iseg-1) + sparsemat%keyg(2,iseg-1) - sparsemat%keyg(1,iseg-1) + 1
      end do

    
    
      if (store_index) then
          ! store the indices of the matrices in the sparse format
          sparsemat%store_index=.true.

    
          ! initialize sparsemat%matrixindex_in_compressed
          !$omp parallel do default(private) shared(sparsemat,norbu) 
          do iorb=1,norbu
             do jorb=1,norbu
                !sparsemat%matrixindex_in_compressed_arr(iorb,jorb)=compressed_index(iorb,jorb,norbu,sparsemat)
                sparsemat%matrixindex_in_compressed_arr(iorb,jorb) = matrixindex_in_compressed(sparsemat, iorb, jorb, .true., norbu)
             end do
          end do
          !$omp end parallel do

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
          ! Otherwise alwyas calculate them on-the-fly
          sparsemat%store_index=.false.
      end if
    

      ! parallelization of matrices, following same idea as norb/norbp/isorb
      !most equal distribution, but want corresponding to norbp for second column
      call init_matrix_parallelization(iproc, nproc, sparsemat%nfvctr, sparsemat%nseg, sparsemat%nvctr, &
           sparsemat%isfvctr_par, sparsemat%nfvctr_par, sparsemat%istsegline, sparsemat%keyv, &
           sparsemat%isvctr, sparsemat%nvctrp, sparsemat%isvctr_par, sparsemat%nvctr_par)
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
      sparsemat%can_use_dense=.false.


      nsegline_mult = f_malloc0(norbu,id='nsegline_mult')
      istsegline_mult = f_malloc(norbu,id='istsegline_mult')
      nseg_mult=0
      nvctr_mult=0
      do iorb=1,norbup
          iiorb=isorbu+iorb
          call create_lookup_table(nnonzero_mult, nonzero_mult, iiorb)
          call nseg_perline(norbu, lut, nseg_mult, nvctr_mult, nsegline_mult(iiorb))
      end do
      if (nproc>1) then
          call mpiallred(nvctr_mult, 1, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(nseg_mult, 1, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(nsegline_mult(1), norbu, mpi_sum, bigdft_mpi%mpi_comm)
      end if



      ! Initialize istsegline, which gives the first segment of each line
      istsegline_mult(1)=1
      do iorb=2,norbu
          istsegline_mult(iorb) = istsegline_mult(iorb-1) + nsegline_mult(iorb-1)
      end do

      keyg_mult = f_malloc0((/2,nseg_mult/),id='keyg_mult')

      ivctr_mult=0
      do iorb=1,norbup
         iiorb=isorbu+iorb
         call create_lookup_table(nnonzero_mult, nonzero_mult, iiorb)
         call keyg_per_line(norbu, nseg_mult, iiorb, istsegline_mult(iiorb), &
              lut, ivctr_mult, keyg_mult)
      end do
      ! check whether the number of elements agrees
      if (nproc>1) then
          call mpiallred(ivctr_mult, 1, mpi_sum, bigdft_mpi%mpi_comm)
      end if
      if (ivctr_mult/=nvctr_mult) then
          write(*,'(a,2i8)') 'ERROR: ivctr_mult/=nvctr_mult', ivctr_mult, nvctr_mult
          stop
      end if
      if (nproc>1) then
          call mpiallred(keyg_mult(1,1), 2*nseg_mult, mpi_sum, bigdft_mpi%mpi_comm)
      end if


      ! Allocate the matrices
      !call allocate_sparse_matrix_matrices(sparsemat, allocate_full)


      ! Initialize the parameters for the spare matrix matrix multiplication
      call init_sparse_matrix_matrix_multiplication(iproc, nproc, norbu, norbup, isorbu, nseg_mult, &
               nsegline_mult, istsegline_mult, keyg_mult, sparsemat)

      call f_free(nsegline_mult)
      call f_free(istsegline_mult)
      call f_free(keyg_mult)
      call f_free(lut)
    
      call timing(iproc,'init_matrCompr','OF')


      contains

        subroutine create_lookup_table(nnonzero, nonzero, iiorb)
          implicit none

          ! Calling arguments
          integer :: nnonzero, iiorb
          integer,dimension(nnonzero) :: nonzero

          ! Local variables
          integer :: ist, iend, i, jjorb

          lut = .false.
          ist=(iiorb-1)*norbu+1
          iend=iiorb*norbu
          do i=1,nnonzero
              if (nonzero(i)<ist) cycle
              if (nonzero(i)>iend) exit
              jjorb=mod(nonzero(i)-1,norbu)+1
              lut(jjorb)=.true.
          end do
        end subroutine create_lookup_table


        subroutine set_value_from_optional()
          if (present(allocate_full_))then
              allocate_full = allocate_full_
          else
              allocate_full = .false.
          end if
          if (present(print_info_))then
              print_info = print_info_
          else
              print_info = .true.
          end if
        end subroutine set_value_from_optional


    
    end subroutine init_sparse_matrix


    subroutine determine_sequential_length(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, &
               sparsemat, nseq, nseq_per_line)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
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
              istart=keyg(1,isegoffset+iseg)
              iend=keyg(2,isegoffset+iseg)
              ! keyg is defined in terms of "global coordinates", so get the
              ! coordinate on a given line by using the mod function
              istart=mod(istart-1,norb)+1
              iend=mod(iend-1,norb)+1
              do iorb=istart,iend
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                          nseq=nseq+1
                          nseqline=nseqline+1
                      end do
                  end do
              end do
         end do
         nseq_per_line(ii)=nseqline
      end do 
    
    end subroutine determine_sequential_length


    subroutine get_nout(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, nout)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
      integer,intent(out) :: nout
    
      ! Local variables
      integer :: i, iii, iseg, iorb
      integer :: isegoffset, istart, iend
    
      nout=0
      do i = 1,norbp
         iii=isorb+i
         isegoffset=istsegline(iii)-1
         do iseg=1,nsegline(iii)
              istart=keyg(1,isegoffset+iseg)
              iend=keyg(2,isegoffset+iseg)
              do iorb=istart,iend
                  nout=nout+1
              end do
          end do
      end do
    
    end subroutine get_nout


    subroutine init_onedimindices_new(norb, norbp, isorb, nseg, nsegline, istsegline, keyg, sparsemat, nout, onedimindices)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
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
              istart=keyg(1,isegoffset+iseg)
              iend=keyg(2,isegoffset+iseg)
              ! keyg is defined in terms of "global coordinates", so get the
              ! coordinate on a given line by using the mod function
              istart=mod(istart-1,norb)+1
              iend=mod(iend-1,norb)+1
              do iorb=istart,iend
                  ii=ii+1
                  onedimindices(1,ii)=i
                  onedimindices(2,ii)=iorb
                  ilen=0
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      ilen=ilen+sparsemat%keyg(2,jseg)-sparsemat%keyg(1,jseg)+1
                  end do
                  onedimindices(3,ii)=ilen
                  onedimindices(4,ii)=itot
                  itot=itot+ilen
              end do
          end do
      end do
    
    end subroutine init_onedimindices_new


    subroutine get_arrays_for_sequential_acces(norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat, nseq, &
               ivectorindex)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg, nseq
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
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
              istart=keyg(1,isegoffset+iseg)
              iend=keyg(2,isegoffset+iseg)
              ! keyg is defined in terms of "global coordinates", so get the
              ! coordinate on a given line by using the mod function
              istart=mod(istart-1,norb)+1
              iend=mod(iend-1,norb)+1
              do iorb=istart,iend
                  !!istindexarr(iorb-istart+1,iseg,i)=ii
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                          jjorb = jorb - (iorb-1)*norb
                          ivectorindex(ii)=jjorb
                          ii = ii+1
                      end do
                  end do
              end do
         end do
      end do 
    
    end subroutine get_arrays_for_sequential_acces


    subroutine init_sequential_acces_matrix(norb, norbp, isorb, nseg, &
               nsegline, istsegline, keyg, sparsemat, nseq, &
               indices_extract_sequential)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norb, norbp, isorb, nseg, nseq
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
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
              istart=keyg(1,isegoffset+iseg)
              iend=keyg(2,isegoffset+iseg)
              ! keyg is defined in terms of "global coordinates", so get the
              ! coordinate on a given line by using the mod function
              istart=mod(istart-1,norb)+1
              iend=mod(iend-1,norb)+1
              do iorb=istart,iend
                  do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                      jj=1
                      do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                          indices_extract_sequential(ii)=sparsemat%keyv(jseg)+jj-1
                          jj = jj+1
                          ii = ii+1
                      end do
                  end do
              end do
         end do
      end do 
    
    end subroutine init_sequential_acces_matrix


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

      ! parallelization of matrices, following same idea as norb/norbp/isorb
      !most equal distribution, but want corresponding to norbp for second column
      do jproc=0,nproc-1
          jst_line = isfvctr_par(jproc)+1
          if (nfvctr_par(jproc)==0) then
             isvctr_par(jproc) = nvctr
          else
             jst_seg = istsegline(jst_line)
             isvctr_par(jproc) = keyv(jst_seg)-1
          end if
      end do
      do jproc=0,nproc-1
         if (jproc==nproc-1) then
            nvctr_par(jproc)=nvctr-isvctr_par(jproc)
         else
            nvctr_par(jproc)=isvctr_par(jproc+1)-isvctr_par(jproc)
         end if
         if (iproc==jproc) isvctr=isvctr_par(jproc)
         if (iproc==jproc) nvctrp=nvctr_par(jproc)
      end do

    end subroutine init_matrix_parallelization


    subroutine init_matrix_taskgroups(iproc, nproc, parallel_layout, collcom, collcom_sr, smat)
      use module_base
      use module_types
      use communications_base, only: comms_linear
      use yaml_output
      implicit none

      ! Caling arguments
      integer,intent(in) :: iproc, nproc
      logical,intent(in) :: parallel_layout
      type(comms_linear),intent(in) :: collcom, collcom_sr
      type(sparse_matrix),intent(inout) :: smat

      ! Local variables
      integer :: ipt, ii, i0, i0i, iiorb, j, i0j, jjorb, ind, ind_min, ind_max, iseq
      integer :: ntaskgroups, jproc, jstart, jend, kkproc, kproc, itaskgroups, lproc, llproc
      integer :: nfvctrp, isfvctr, isegstart, isegend, jorb
      integer,dimension(:,:),allocatable :: iuse_startend, itaskgroups_startend, ranks
      integer,dimension(:),allocatable :: tasks_per_taskgroup
      integer :: ntaskgrp_calc, ntaskgrp_use, i, ncount, iitaskgroup, group, ierr, iitaskgroups, newgroup, iseg
      logical :: go_on
      integer,dimension(:,:),allocatable :: in_taskgroup

      call f_routine(id='init_matrix_taskgroups')

      ! First determine the minimal and maximal value oft the matrix which is used by each process
      iuse_startend = f_malloc0((/1.to.2,0.to.nproc-1/),id='iuse_startend')


      if (parallel_layout) then
          ! The matrices can be parallelized

          ind_min = smat%nvctr
          ind_max = 0

          ! The operations done in the transposed wavefunction layout
          call check_transposed_layout()

          ! Now check the compress_distributed layout
          call check_compress_distributed_layout()

          ! Now check the matrix matrix multiplications layout
          call check_matmul_layout()

          ! Now check the sumrho operations
          call check_sumrho_layout()
      else
          ! The matrices can not be parallelized
          ind_min = 1
          ind_max = smat%nvctr
      end if

      ! Now the minimal and maximla values are known
      iuse_startend(1,iproc) = ind_min
      iuse_startend(2,iproc) = ind_max
      if (nproc>1) then
          call mpiallred(iuse_startend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm)
      end if
 
      !if (iproc==0) write(*,'(a,100(2i7,4x))') 'iuse_startend',iuse_startend
 
      ntaskgroups = 1
      llproc=0
      ii = 0
      do 
          jproc = llproc + ii
          if (jproc==nproc-1) exit
          jstart = iuse_startend(1,jproc) !beginning of part used by task jproc
          jend = iuse_startend(2,jproc) !end of part used by task jproc
          ii = ii + 1
          !search the last process whose calculation stops prior to iend
          go_on = .true.
          do lproc=nproc-1,0,-1
              if (iuse_startend(1,lproc)<=jend) then
                  !if (iproc==0) write(*,'(a,3i8)') 'lproc, iuse_startend(1,lproc), iuse_startend(2,llproc)', lproc, iuse_startend(1,lproc), iuse_startend(2,llproc)
                  if (iuse_startend(1,lproc)<=iuse_startend(2,llproc)) then
                      go_on = .false.
                  end if
                  exit
              end if
          end do
          !if (iproc==0) write(*,*) '2: llproc, ii, jproc, go_on', llproc, ii, jproc, go_on
          if (.not.go_on) cycle
          ntaskgroups = ntaskgroups + 1
          ! Search the starting point of the next taskgroups, defined as the
          ! largest starting part which is smaller than jend
          llproc=nproc-1
          do lproc=nproc-1,0,-1
              if (iuse_startend(1,lproc)<=jend) then
                  llproc = lproc
                  exit
              end if
          end do
          !if (iproc==0) write(*,*) 'llproc, ii, jproc, ntaskgroups', llproc, ii, jproc, ntaskgroups
          ii = 0
          !if (llproc==nproc-1) exit
      end do
      !if (iproc==0) write(*,*) 'iproc, ntaskgroups', iproc, ntaskgroups

      smat%ntaskgroup = ntaskgroups
 
      itaskgroups_startend = f_malloc0((/2,ntaskgroups/),id='itaskgroups_startend')
      itaskgroups_startend(1,1) = 1
      itaskgroups = 1
      llproc=0
      ii = 0
      do 
          jproc = llproc + ii
          if (jproc==nproc-1) exit
          jstart = iuse_startend(1,jproc) !beginning of part used by task jproc
          jend = iuse_startend(2,jproc) !end of part used by task jproc
          ii = ii + 1
          !search the last process whose calculation stops prior to jend
          go_on = .true.
          do lproc=nproc-1,0,-1
              if (iuse_startend(1,lproc)<=jend) then
                  if (iuse_startend(1,lproc)<=iuse_startend(2,llproc)) then
                      go_on = .false.
                  end if
                  exit
              end if
          end do
          if (.not.go_on) cycle
          itaskgroups_startend(2,itaskgroups) = jend
          itaskgroups = itaskgroups + 1
          ! Search the starting point of the next taskgroups, defined as the
          ! largest starting part which is smaller than jend
          llproc=nproc-1
          do lproc=nproc-1,0,-1
              if (iuse_startend(1,lproc)<=jend) then
                  itaskgroups_startend(1,itaskgroups) = iuse_startend(1,lproc)
                  llproc = lproc
                  exit
              end if
          end do
          ii = 0
          !if (llproc==nproc-1) exit
      end do
      itaskgroups_startend(2,itaskgroups) = iuse_startend(2,nproc-1)
      if (itaskgroups/=ntaskgroups) stop 'itaskgroups/=ntaskgroups'
      !if (iproc==0) write(*,'(a,i8,4x,1000(2i7,4x))') 'iproc, itaskgroups_startend', itaskgroups_startend
 
      ! Assign the processes to the taskgroups
      ntaskgrp_calc = 0
      ntaskgrp_use = 0
      do itaskgroups=1,ntaskgroups
          if ( iuse_startend(1,iproc)<=itaskgroups_startend(2,itaskgroups) .and.  &
               iuse_startend(2,iproc)>=itaskgroups_startend(1,itaskgroups) ) then
              !write(*,'(2(a,i0))') 'USE: task ',iproc,' is in taskgroup ',itaskgroups
               ntaskgrp_use = ntaskgrp_use + 1
          end if
      end do
      if (ntaskgrp_use>2) stop 'ntaskgrp_use>2'

      smat%ntaskgroupp = max(ntaskgrp_calc,ntaskgrp_use)

      smat%taskgroup_startend = f_malloc_ptr((/2,2,smat%ntaskgroup/),id='smat%taskgroup_startend')
      smat%inwhichtaskgroup = f_malloc_ptr((/smat%ntaskgroupp/),id='smat%smat%inwhichtaskgroup')


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
               smat%inwhichtaskgroup(i) = itaskgroups
          end if
      end do
      if (i/=smat%ntaskgroupp) then
          write(*,*) 'i, smat%ntaskgroupp', i, smat%ntaskgroupp
          stop 'i/=smat%ntaskgroupp'
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

      !Check
      ncount = 0
      do itaskgroups=1,smat%ntaskgroup
          ncount = ncount + smat%taskgroup_startend(2,2,itaskgroups)-smat%taskgroup_startend(1,2,itaskgroups)+1
      end do
      if (ncount/=smat%nvctr) then
          write(*,*) 'ncount, smat%nvctr', ncount, smat%nvctr
          stop 'ncount/=smat%nvctr'
      end if



      ! Create the taskgroups
      ! Count the number of tasks per taskgroup
      tasks_per_taskgroup = f_malloc0(smat%ntaskgroup,id='tasks_per_taskgroup')
      do itaskgroups=1,smat%ntaskgroupp
          iitaskgroup = smat%inwhichtaskgroup(itaskgroups)
          tasks_per_taskgroup(iitaskgroup) = tasks_per_taskgroup(iitaskgroup) + 1
      end do
      if (nproc>1) then
          call mpiallred(tasks_per_taskgroup(1), smat%ntaskgroup, mpi_sum, bigdft_mpi%mpi_comm)
      end if
      !if (iproc==0) write(*,'(a,i7,4x,1000i7)') 'iproc, tasks_per_taskgroup', iproc, tasks_per_taskgroup
      call mpi_comm_group(bigdft_mpi%mpi_comm, group, ierr)

      in_taskgroup = f_malloc0((/0.to.nproc-1,1.to.smat%ntaskgroup/),id='in_taskgroup')
      ranks = f_malloc((/maxval(tasks_per_taskgroup),smat%ntaskgroup/),id='ranks')
      do itaskgroups=1,smat%ntaskgroupp
          iitaskgroups = smat%inwhichtaskgroup(itaskgroups)
          in_taskgroup(iproc,iitaskgroups) = 1
      end do
      if (nproc>1) then
          call mpiallred(in_taskgroup(0,1), nproc*smat%ntaskgroup, mpi_sum, bigdft_mpi%mpi_comm)
      end if

      allocate(smat%mpi_groups(smat%ntaskgroup))
      do itaskgroups=1,smat%ntaskgroup
          smat%mpi_groups(itaskgroups) = mpi_environment_null()
      end do
      do itaskgroups=1,smat%ntaskgroup
          ii = 0
          do jproc=0,nproc-1
              if (in_taskgroup(jproc,itaskgroups)>0) then
                  ii = ii + 1
                  ranks(ii,itaskgroups) = jproc
              end if
          end do
          if (ii/=tasks_per_taskgroup(itaskgroups)) stop 'ii/=tasks_per_taskgroup(itaskgroups)'
          call mpi_group_incl(group, ii, ranks(1,itaskgroups), newgroup, ierr)
          call mpi_comm_create(bigdft_mpi%mpi_comm, newgroup, smat%mpi_groups(itaskgroups)%mpi_comm, ierr)
          if (smat%mpi_groups(itaskgroups)%mpi_comm/=MPI_COMM_NULL) then
              call mpi_comm_size(smat%mpi_groups(itaskgroups)%mpi_comm, smat%mpi_groups(itaskgroups)%nproc, ierr)
              call mpi_comm_rank(smat%mpi_groups(itaskgroups)%mpi_comm, smat%mpi_groups(itaskgroups)%iproc, ierr)
          end if
          smat%mpi_groups(itaskgroups)%igroup = itaskgroups
          smat%mpi_groups(itaskgroups)%ngroup = smat%ntaskgroup
          call mpi_group_free(newgroup, ierr)
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
              call yaml_map('IDs',ranks(1:tasks_per_taskgroup(itaskgroups),itaskgroups))
              call yaml_newline()
              call yaml_map('start / end',smat%taskgroup_startend(1:2,1,itaskgroups))
              call yaml_map('start / end disjoint',smat%taskgroup_startend(1:2,2,itaskgroups))
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()
      end if

      call f_free(in_taskgroup)
      call f_free(iuse_startend)
      call f_free(itaskgroups_startend)
      call f_free(tasks_per_taskgroup)
      call f_free(ranks)





      call f_release_routine()


      contains

        subroutine check_transposed_layout()
          do ipt=1,collcom%nptsp_c
              ii=collcom%norb_per_gridpoint_c(ipt)
              i0 = collcom%isptsp_c(ipt)
              do i=1,ii
                  i0i=i0+i
                  iiorb=collcom%indexrecvorbital_c(i0i)
                  do j=1,ii
                      i0j=i0+j
                      jjorb=collcom%indexrecvorbital_c(i0j)
                      ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                      ind_min = min(ind_min,ind)
                      ind_max = max(ind_max,ind)
                  end do
              end do
          end do
          do ipt=1,collcom%nptsp_f
              ii=collcom%norb_per_gridpoint_f(ipt)
              i0 = collcom%isptsp_f(ipt)
              do i=1,ii
                  i0i=i0+i
                  iiorb=collcom%indexrecvorbital_f(i0i)
                  do j=1,ii
                      i0j=i0+j
                      jjorb=collcom%indexrecvorbital_f(i0j)
                      ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                      ind_min = min(ind_min,ind)
                      ind_max = max(ind_max,ind)
                  end do
              end do
          end do

          ! Store these values
          smat%istartend_t(1) = ind_min
          smat%istartend_t(2) = ind_max
          ! Determine to which segments this corresponds
          do iseg=1,smat%nseg
              if (smat%keyv(iseg)+smat%keyg(2,iseg)-smat%keyg(1,iseg)>=smat%istartend_t(1)) then
                  smat%istartendseg_t(1)=iseg
                  exit
              end if
          end do
          do iseg=smat%nseg,1,-1
              if (smat%keyv(iseg)<=smat%istartend_t(2)) then
                  smat%istartendseg_t(2)=iseg
                  exit
              end if
          end do
        end subroutine check_transposed_layout


        subroutine check_compress_distributed_layout()
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
                  do iseg=isegstart,isegend
                      ii=smat%keyv(iseg)-1
                      do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                          ii=ii+1
                          ind_min = min(ii,ind_min)
                          ind_max = min(ii,ind_max)
                      end do
                  end do
              end if
          end do
        end subroutine check_compress_distributed_layout


        subroutine check_matmul_layout()
          do iseq=1,smat%smmm%nseq
              ind=smat%smmm%indices_extract_sequential(iseq)
              ind_min = min(ind_min,ind)
              ind_max = max(ind_max,ind)
          end do
        end subroutine check_matmul_layout

        subroutine check_sumrho_layout()
          do ipt=1,collcom_sr%nptsp_c
              ii=collcom_sr%norb_per_gridpoint_c(ipt)
              i0=collcom_sr%isptsp_c(ipt)
              do i=1,ii
                  iiorb=collcom_sr%indexrecvorbital_c(i0+i)
                  ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
                  ind_min = min(ind_min,ind)
                  ind_max = max(ind_max,ind)
              end do
          end do
        end subroutine check_sumrho_layout


      !!  function get_start_of_segment(smat, iiseg) result(ist)

      !!      do iseg=smat%nseg,1,-1
      !!          if (iiseg>=smat%keyv(iseg)) then
      !!              it = smat%keyv(iseg)
      !!              exit
      !!          end if
      !!      end do

      !!  end function get_start_of_segment
 
    end subroutine init_matrix_taskgroups

end module sparsematrix_init

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


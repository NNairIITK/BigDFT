!> @file
!!  File defining the structures to deal with the sparse matrices
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Module to deal with the sparse matrices
module sparsematrix

  use module_base
  use sparsematrix_base

  implicit none

  private



  !> Public routines
  public :: compress_matrix, compress_matrix2
  public :: uncompress_matrix, uncompress_matrix2
  public :: check_matrix_compression
  public :: transform_sparse_matrix, transform_sparse_matrix_local
  public :: compress_matrix_distributed_wrapper
  public :: uncompress_matrix_distributed2
  public :: sequential_acces_matrix_fast, sequential_acces_matrix_fast2
  public :: sparsemm_new
  public :: gather_matrix_from_taskgroups, gather_matrix_from_taskgroups_inplace
  public :: extract_taskgroup_inplace, extract_taskgroup
  public :: write_matrix_compressed
  public :: check_symmetry
  public :: write_sparsematrix
  public :: write_sparsematrix_CCS
  public :: transform_sparsity_pattern
  public :: matrix_matrix_mult_wrapper
  public :: trace_sparse


  interface compress_matrix_distributed_wrapper
      module procedure compress_matrix_distributed_wrapper_1
      module procedure compress_matrix_distributed_wrapper_2
  end interface compress_matrix_distributed_wrapper

  contains

    !> subroutine to compress the matrix to sparse form
    subroutine compress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix),intent(in) :: sparsemat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),target,intent(in) :: inmat
      real(kind=8),dimension(sparsemat%nvctr*sparsemat%nspin),target,intent(out) :: outmat
    
      ! Local variables
      integer :: iseg, j, jj, irow, jcol, jjj, ierr, ishift, ispin
      real(kind=8),dimension(:,:,:),pointer :: inm
      real(kind=8),dimension(:),pointer :: outm
      integer,dimension(2) :: irowcol

      !if (present(outmat)) then
      !    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !        stop 'outmat not allowed for the given options'
      !    end if
          outm => outmat
      !else
      !    outm => sparsemat%matrix_compr
      !end if

      !if (present(inmat)) then
      !    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !        stop 'in not allowed for the given options'
      !    end if
          inm => inmat
      !else
      !    inm => sparsemat%matrix
      !end if
    
      call timing(iproc,'compressd_mcpy','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         do ispin=1,sparsemat%nspin
             ishift=(ispin-1)*sparsemat%nvctr
             !OpenMP broken on Vesta
             !$omp parallel default(none) private(iseg,j,jj,irowcol) &
             !$omp shared(sparsemat,inm,outm,ishift,ispin)
             !$omp do
             do iseg=1,sparsemat%nseg
                 jj=sparsemat%keyv(iseg)
                 ! A segment is always on one line, therefore no double loop
                 do j=sparsemat%keyg(1,1,iseg),sparsemat%keyg(2,1,iseg)
                    !irow = sparsemat%orb_from_index(1,jj)
                    !jcol = sparsemat%orb_from_index(2,jj)
                    !irowcol = orb_from_index(sparsemat, j)
                    !write(*,*) 'iseg, j, jj', iseg, j, jj
                    outm(jj+ishift)=inm(j,sparsemat%keyg(1,2,iseg),ispin)
                    jj=jj+1
                 end do
             end do
             !$omp end do
             !$omp end parallel
         end do
      else if (sparsemat%parallel_compression==1) then
         stop 'this needs to be fixed'
         !!#call to_zero(sparsemat%nvctr, sparsemat%matrix_compr(1))
         !!#!$omp parallel do default(private) shared(sparsemat)
         !!#do jj=1,sparsemat%nvctrp
         !!#   jjj=jj+sparsemat%isvctr
         !!#   irow = sparsemat%orb_from_index(1,jjj)
         !!#   jcol = sparsemat%orb_from_index(2,jjj)
         !!#   sparsemat%matrix_compr(jjj)=sparsemat%matrix(irow,jcol)
         !!#end do
         !!#!$omp end parallel do
         !!#if (bigdft_mpi%nproc > 1) then
         !!#   call mpiallred(sparsemat%matrix_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
         !!#end if
      else
         stop 'this needs to be fixed'
         !!#sparsemat%matrix_comprp=f_malloc_ptr((sparsemat%nvctrp),id='sparsemat%matrix_comprp')
         !!#!$omp parallel do default(private) shared(sparsemat)
         !!#do jj=1,sparsemat%nvctrp
         !!#   jjj=jj+sparsemat%isvctr
         !!#   irow = sparsemat%orb_from_index(1,jjj)
         !!#   jcol = sparsemat%orb_from_index(2,jjj)
         !!#   sparsemat%matrix_comprp(jj)=sparsemat%matrix(irow,jcol)
         !!#end do
         !!#!$omp end parallel do
         !!#if (bigdft_mpi%nproc > 1) &
         !!#   & call mpi_allgatherv(sparsemat%matrix_comprp, sparsemat%nvctrp, &
         !!#   &    mpi_double_precision, sparsemat%matrix_compr, sparsemat%nvctr_par(:), &
         !!#   &    sparsemat%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         !!#call f_free_ptr(sparsemat%matrix_comprp)
      end if
    
      call timing(iproc,'compressd_mcpy','OF')
    
    end subroutine compress_matrix


    subroutine compress_matrix2(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),intent(in) :: inmat
      real(kind=8),dimension(sparsemat%nvctrp_tg*sparsemat%nspin),intent(out) :: outmat

      ! Local variables
      real(kind=8),dimension(:),allocatable :: tmparr

      tmparr = sparsematrix_malloc(sparsemat,iaction=SPARSE_FULL,id='tmparr')
      call compress_matrix(iproc,sparsemat,inmat,tmparr)
      call extract_taskgroup(sparsemat, tmparr, outmat)
      call f_free(tmparr)
    end subroutine compress_matrix2


    !> subroutine to uncompress the matrix from sparse form
    subroutine uncompress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix), intent(inout) :: sparsemat
      real(kind=8),dimension(sparsemat%nvctr*sparsemat%nspin),target,intent(in) :: inmat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),target,intent(inout) :: outmat
      
      ! Local variables
      integer :: iseg, i, ii, irow, jcol, iii, ierr, ishift, ispin
      real(kind=8),dimension(:),pointer :: inm
      real(kind=8),dimension(:,:,:),pointer :: outm
      integer,dimension(2) :: irowcol

      !!if (present(outmat)) then
      !!    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !!        stop 'outmat not allowed for the given options'
      !!    end if
          outm => outmat
      !!else
      !!    outm => sparsemat%matrix
      !!end if

      !!if (present(inmat)) then
      !!    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !!        stop 'inmat not allowed for the given options'
      !!    end if
          inm => inmat
      !!else
      !!    inm => sparsemat%matrix_compr
      !!end if
    
      call timing(iproc,'compressd_mcpy','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         !call to_zero(sparsemat%nfvctr**2*sparsemat%nspin, outm(1,1,1))
         call f_zero(outm)
         do ispin=1,sparsemat%nspin
             ishift=(ispin-1)*sparsemat%nvctr
             !OpenMP broken on Vesta
             !$omp parallel default(none) private(iseg,i,ii,irowcol) shared(sparsemat,inm,outm,ispin,ishift)
             !$omp do
             do iseg=1,sparsemat%nseg
                 ii=sparsemat%keyv(iseg)
                 ! A segment is always on one line, therefore no double loop
                 do i=sparsemat%keyg(1,1,iseg),sparsemat%keyg(2,1,iseg)
                    !irow = sparsemat%orb_from_index(1,ii)
                    !jcol = sparsemat%orb_from_index(2,ii)
                    !irowcol = orb_from_index(sparsemat, i)
                    outm(i,sparsemat%keyg(1,2,iseg),ispin)=inm(ii+ishift)
                    ii=ii+1
                end do
             end do
             !$omp end do
             !$omp end parallel
         end do
      else if (sparsemat%parallel_compression==1) then
         stop 'needs to be fixed'
         !!call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
         !!!$omp parallel do default(private) shared(sparsemat)
         !!do ii=1,sparsemat%nvctrp
         !!   iii=ii+sparsemat%isvctr
         !!   irow = sparsemat%orb_from_index(1,iii)
         !!   jcol = sparsemat%orb_from_index(2,iii)
         !!   sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(iii)
         !!end do
         !!!$omp end parallel do
         !!if (bigdft_mpi%nproc > 1) then
         !!   call mpiallred(sparsemat%matrix(1,1), sparsemat%nfvctr**2,mpi_sum,bigdft_mpi%mpi_comm)
         !!end if
      else
         stop 'needs to be fixed'
         !!sparsemat%matrixp=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
         !!call to_zero(sparsemat%nfvctr*sparsemat%nfvctrp, sparsemat%matrixp(1,1))
         !!!$omp parallel do default(private) shared(sparsemat)
         !!do ii=1,sparsemat%nvctrp
         !!   iii=ii+sparsemat%isvctr
         !!   irow = sparsemat%orb_from_index(1,iii)
         !!   jcol = sparsemat%orb_from_index(2,iii) - sparsemat%isfvctr
         !!   sparsemat%matrixp(irow,jcol)=sparsemat%matrix_compr(iii)
         !!end do
         !!!$omp end parallel do
         !!if (bigdft_mpi%nproc > 1) &
         !!   & call mpi_allgatherv(sparsemat%matrixp, sparsemat%nfvctr*sparsemat%nfvctrp, mpi_double_precision, &
         !!   &   sparsemat%matrix, sparsemat%nfvctr*sparsemat%nfvctr_par(:), sparsemat%nfvctr*sparsemat%isfvctr_par, &
         !!   &   mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         !!call f_free_ptr(sparsemat%matrixp)
      end if
      sparsemat%can_use_dense=.true.  
    
      call timing(iproc,'compressd_mcpy','OF')
    
    end subroutine uncompress_matrix


    subroutine uncompress_matrix2(iproc, nproc, smat, matrix_compr, matrix)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(inout) :: smat
      real(kind=8),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: matrix_compr
      real(kind=8),dimension(smat%nfvctr,smat%nfvctr,smat%nspin),intent(out) :: matrix

      ! Local variables
      real(kind=8),dimension(:),allocatable :: tmparr

      tmparr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='tmparr')
      call gather_matrix_from_taskgroups(iproc, nproc, smat, matrix_compr, tmparr)
      call uncompress_matrix(iproc, smat, inmat=tmparr, outmat=matrix)
      call f_free(tmparr)
    end subroutine uncompress_matrix2


    subroutine check_matrix_compression(iproc, sparsemat, mat)
      use yaml_output
      implicit none
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
      type(matrices),intent(inout) :: mat
      !Local variables
      character(len=*), parameter :: subname='check_matrix_compression'
      real(kind=8), parameter :: tol=1.e-10
      integer :: jorb, irow, icol, iseg, ii
      real(kind=8) :: maxdiff
    
      call f_routine('check_matrix_compression')

      call f_free_ptr(mat%matrix_compr)
      mat%matrix_compr = sparsematrix_malloc_ptr(sparsemat,iaction=SPARSE_FULL,id='mat%matrix_compr')
    
      mat%matrix = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_FULL, id='mat%matrix')
    
      !call to_zero(sparsemat%nfvctr**2*sparsemat%nspin,mat%matrix(1,1,1))
      call f_zero(mat%matrix)
      do iseg = 1, sparsemat%nseg
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            !print *,'jorb, irow,icol',jorb, irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
            !SM: need to fix spin 
            mat%matrix(irow,icol,1) = test_value_matrix(sparsemat%nfvctr, irow, icol)
         end do
      end do
      
      call compress_matrix(iproc, sparsemat, inmat=mat%matrix, outmat=mat%matrix_compr)
      !write(*,*) 'mat%matrix',mat%matrix
      !write(*,*) 'mat%matrix_compr',mat%matrix_compr
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ii=0
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            !write(*,'(a,4i8,2es13.3)') 'jorb, irow, icol, sparsemat%keyv(iseg)+ii, val, ref', jorb, irow, icol, sparsemat%keyv(iseg)+ii, mat%matrix_compr(sparsemat%keyv(iseg)+ii), test_value_matrix(sparsemat%nfvctr, irow, icol)
            maxdiff = max(abs(mat%matrix_compr(sparsemat%keyv(iseg)+ii)&
                 -test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff)
            ii=ii+1
         end do
      end do
    
      if (iproc==0) call yaml_map('Tolerances for this check',tol,fmt='(1pe25.17)')
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('COMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for compress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call uncompress_matrix(iproc, sparsemat, inmat=mat%matrix_compr, outmat=mat%matrix)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            maxdiff = max(abs(mat%matrix(irow,icol,1)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
         end do
      end do
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('UNCOMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for uncompress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call f_free_ptr(mat%matrix)
      !!call f_free_ptr(sparsemat%matrix_compr)

      call f_free_ptr(mat%matrix_compr)
      mat%matrix_compr = sparsematrix_malloc_ptr(sparsemat,iaction=SPARSE_TASKGROUP,id='mat%matrix_compr')

      call f_release_routine()
    
    contains
       !> define a value for the wavefunction which is dependent of the indices
       function test_value_matrix(norb,iorb,jorb)
          use module_base
          implicit none
          integer, intent(in) :: norb,iorb,jorb
          real(kind=8) :: test_value_matrix
    
          test_value_matrix = norb*(iorb-1)+jorb
          !print *,iorb,jorb,test_value_matrix
       END FUNCTION test_value_matrix
    
       subroutine get_indices(ind,irow,icol)
         implicit none
         integer, intent(in) :: ind
         integer, intent(out) :: irow, icol
    
         icol = (ind - 1) / sparsemat%nfvctr + 1
         irow = ind - (icol-1)*sparsemat%nfvctr
         !print *,'irow,icol',irow,icol
       END SUBROUTINE get_indices 
    end subroutine check_matrix_compression


    subroutine transform_sparse_matrix(smat, lmat, smatrix_compr, lmatrix_compr, cmode)
      use module_base
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat, lmat
      real(kind=8),dimension(smat%nspin*smat%nvctr),intent(inout) :: smatrix_compr
      real(kind=8),dimension(lmat%nspin*lmat%nvctr),intent(inout) :: lmatrix_compr
      character(len=14),intent(in) :: cmode
    
      ! Local variables
      integer(kind=8) :: isstart, isend, ilstart, ilend, iostart, ioend
      integer :: imode, icheck, isseg, ilseg
      integer :: ilength, iscostart, ilcostart, i
      integer :: ilsegstart, ispin, isshift, ilshift, isoffset, iloffset
      integer,parameter :: SMALL_TO_LARGE=1
      integer,parameter :: LARGE_TO_SMALL=2
    
      call f_routine(id='transform_sparse_matrix')
    
      ! determine the case:
      ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
      ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
      if (cmode=='small_to_large' .or. cmode=='SMALL_TO_LARGE') then
          imode=SMALL_TO_LARGE
      else if (cmode=='large_to_small' .or. cmode=='LARGE_TO_SMALL') then
          imode=LARGE_TO_SMALL
      else
          stop 'wrong cmode'
      end if
    
      select case (imode)
      case (SMALL_TO_LARGE)
         !call to_zero(lmat%nvctr*lmat%nspin,lmatrix_compr(1))
         call f_zero(lmatrix_compr)
      case (LARGE_TO_SMALL)
         !call to_zero(smat%nvctr*lmat%nspin,smatrix_compr(1))
         call f_zero(smatrix_compr)
      case default
          stop 'wrong imode'
      end select
    
      call timing(bigdft_mpi%iproc,'transform_matr','IR')


      icheck=0
      do ispin=1,smat%nspin

          isshift=(ispin-1)*smat%nvctr
          ilshift=(ispin-1)*lmat%nvctr
    
          ilsegstart=1
          !$omp parallel default(private) &
          !$omp shared(smat, lmat, imode, lmatrix_compr, smatrix_compr, icheck, isshift, ilshift) &
          !$omp firstprivate(ilsegstart)
          !$omp do reduction(+:icheck)
          sloop: do isseg=1,smat%nseg
              isstart = int((smat%keyg(1,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) + int(smat%keyg(1,1,isseg),kind=8)
              isend = int((smat%keyg(2,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) + int(smat%keyg(2,1,isseg),kind=8)
              ! A segment is always on one line, therefore no double loop
              lloop: do ilseg=ilsegstart,lmat%nseg
                  ilstart = int((lmat%keyg(1,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) + int(lmat%keyg(1,1,ilseg),kind=8)
                  ilend = int((lmat%keyg(2,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) + int(lmat%keyg(2,1,ilseg),kind=8)
    
                  ! check whether there is an overlap:
                  ! if not, increase loop counters
                  if (ilstart>isend) then
                      !ilsegstart=ilseg
                      exit lloop
                  end if
                  if (isstart>ilend) then
                      ilsegstart=ilseg
                      cycle lloop
                  end if
                  ! if yes, determine start end end of overlapping segment (in uncompressed form)
                  iostart=max(isstart,ilstart)
                  ioend=min(isend,ilend)
                  ilength=ioend-iostart+1
    
                  ! offset with respect to the starting point of the segment
                  isoffset = int(iostart - &
                             (int((smat%keyg(1,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) &
                               + int(smat%keyg(1,1,isseg),kind=8)),kind=4)
                  iloffset = int(iostart - &
                             (int((lmat%keyg(1,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) &
                               + int(lmat%keyg(1,1,ilseg),kind=8)),kind=4)
    
                  ! determine start end and of the overlapping segment in compressed form
                  iscostart=smat%keyv(isseg)+isoffset
                  ilcostart=lmat%keyv(ilseg)+iloffset
    
                  ! copy the elements
                  select case (imode)
                  case (SMALL_TO_LARGE) 
                      do i=0,ilength-1
                          lmatrix_compr(ilcostart+i+ilshift)=smatrix_compr(iscostart+i+isshift)
                      end do
                  case (LARGE_TO_SMALL) 
                      do i=0,ilength-1
                          smatrix_compr(iscostart+i+isshift)=lmatrix_compr(ilcostart+i+ilshift)
                      end do
                  case default
                      stop 'wrong imode'
                  end select
                  icheck=icheck+ilength
              end do lloop
          end do sloop
          !$omp end do 
          !$omp end parallel

      end do
    
      ! all elements of the small matrix must have been processed, no matter in
      ! which direction the transformation has been executed
      if (icheck/=smat%nvctr*smat%nspin) then
          write(*,'(a,2i8)') 'ERROR: icheck/=smat%nvctr*smat%nspin', icheck, smat%nvctr*smat%nspin
          stop
      end if

      call timing(bigdft_mpi%iproc,'transform_matr','RS')
      call f_release_routine()
    
    end subroutine transform_sparse_matrix



    subroutine transform_sparse_matrix_local(smat, lmat, smatrix_compr, lmatrix_compr, cmode)
      use module_base
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat, lmat
      real(kind=8),dimension(smat%nspin*smat%nvctrp_tg),intent(inout) :: smatrix_compr
      real(kind=8),dimension(lmat%nspin*lmat%nvctrp_tg),intent(inout) :: lmatrix_compr
      character(len=14),intent(in) :: cmode
    
      ! Local variables
      real(kind=8),dimension(:),allocatable :: tmparrs, tmparrl
      integer :: ishift_src, ishift_dst, imode, ispin
      integer,parameter :: SMALL_TO_LARGE=1
      integer,parameter :: LARGE_TO_SMALL=2
    
      call f_routine(id='transform_sparse_matrix_local')


      ! determine the case:
      ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
      ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
      if (cmode=='small_to_large' .or. cmode=='SMALL_TO_LARGE') then
          imode=SMALL_TO_LARGE
      else if (cmode=='large_to_small' .or. cmode=='LARGE_TO_SMALL') then
          imode=LARGE_TO_SMALL
      else
          stop 'wrong cmode'
      end if

    
      select case (imode)
      case (SMALL_TO_LARGE)
          tmparrs = sparsematrix_malloc0(smat,iaction=SPARSE_FULL,id='tmparrs')
          tmparrl = sparsematrix_malloc(lmat,iaction=SPARSE_FULL,id='tmparrl')
          do ispin=1,smat%nspin
              ishift_src = (ispin-1)*smat%nvctrp_tg
              ishift_dst = (ispin-1)*smat%nvctr
              call vcopy(smat%nvctrp_tg, smatrix_compr(ishift_src+1), 1, &
                   tmparrs(ishift_dst+smat%isvctrp_tg+1), 1)
              call transform_sparse_matrix(smat, lmat, tmparrs, tmparrl, cmode)
              call extract_taskgroup(lmat, tmparrl, lmatrix_compr)
          end do
      case (LARGE_TO_SMALL)
          tmparrs = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='tmparrs')
          tmparrl = sparsematrix_malloc0(lmat,iaction=SPARSE_FULL,id='tmparrl')
          do ispin=1,smat%nspin
              ishift_src = (ispin-1)*lmat%nvctrp_tg
              ishift_dst = (ispin-1)*lmat%nvctr
              call vcopy(lmat%nvctrp_tg, lmatrix_compr(ishift_src+1), 1, &
                   tmparrl(ishift_dst+lmat%isvctrp_tg+1), 1)
              call transform_sparse_matrix(smat, lmat, tmparrs, tmparrl, cmode)
              call extract_taskgroup(smat, tmparrs, smatrix_compr)
          end do
      case default
          stop 'wrong imode'
      end select

      call f_free(tmparrs)
      call f_free(tmparrl)

      call f_release_routine()
    
  end subroutine transform_sparse_matrix_local





   subroutine compress_matrix_distributed_wrapper_2(iproc, nproc, smat, layout, matrixp, matrix_compr)
     use module_base
     !!use yaml_output
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(:,:),intent(in) :: matrixp
     real(kind=8),dimension(smat%nvctrp_tg),target,intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr, nvctrp, ierr, isvctr
     integer :: ncount, itg, iitg, ist_send, ist_recv
     integer :: window, sizeof, jproc_send, iorb, jproc, info
     integer,dimension(:),pointer :: isvctr_par, nvctr_par
     integer,dimension(:),allocatable :: request, windows
     real(kind=8),dimension(:),pointer :: matrix_local
     real(kind=8),dimension(:),allocatable :: recvbuf

     call f_routine(id='compress_matrix_distributed_wrapper_2')

     call timing(iproc,'compressd_mcpy','ON')

     ! Check the dimensions of the input array and assign some values
     !if (size(matrixp,1)/=smat%nfvctr) stop 'size(matrixp,1)/=smat%nfvctr'
     if (size(matrixp,1)/=smat%nfvctr) then
         call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp,1),fmt='(i0)'))//&
              &' instead of '//trim(yaml_toa(smat%nfvctr,fmt='(i0)')), &
              err_name='BIGDFT_RUNTIME_ERROR')
     end if
     if (layout==DENSE_PARALLEL) then
         if (size(matrixp,2)/=smat%nfvctrp) stop '(ubound(matrixp,2)/=smat%nfvctrp'
         nfvctrp = smat%nfvctrp
         isfvctr = smat%isfvctr
         nvctrp = smat%nvctrp
         isvctr = smat%isvctr
         isvctr_par => smat%isvctr_par
         nvctr_par => smat%nvctr_par
     else if (layout==DENSE_MATMUL) then
         if (size(matrixp,2)/=smat%smmm%nfvctrp) stop '(ubound(matrixp,2)/=smat%smmm%nfvctrp'
         nfvctrp = smat%smmm%nfvctrp
         isfvctr = smat%smmm%isfvctr
         nvctrp = smat%smmm%nvctrp_mm
         isvctr = smat%smmm%isvctr_mm
         isvctr_par => smat%smmm%isvctr_mm_par
         nvctr_par => smat%smmm%nvctr_mm_par
     else
         call f_err_throw('layout has the value '//trim(yaml_toa(layout,fmt='(i0)'))//&
              &'; allowed are '//trim(yaml_toa(DENSE_PARALLEL,fmt='(i0)'))//&
              &' and '//trim(yaml_toa(DENSE_MATMUL,fmt='(i0)')), &
              err_name='BIGDFT_RUNTIME_ERROR')
     end if



     !@ NEW #####################
     matrix_local = f_malloc_ptr(max(1,nvctrp),id='matrix_local')
     if (layout==DENSE_PARALLEL) then
         ii = 0
         if (nfvctrp>0) then
             isegstart=smat%istsegline(isfvctr+1)
             isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
             do iseg=isegstart,isegend
                 ! A segment is always on one line, therefore no double loop
                 do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                     iorb = smat%keyg(1,2,iseg)
                     ii = ii + 1
                     matrix_local(ii) = matrixp(jorb,iorb-isfvctr)
                 end do
             end do
         end if
         if (ii/=nvctrp) stop 'compress_matrix_distributed: ii/=nvctrp'
     else if (layout==DENSE_MATMUL) then
         ii = 0
         if (nvctrp>0) then
             isegstart=smat%istsegline(isfvctr+1)
             isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
             do iseg=isegstart,isegend
                 ! A segment is always on one line, therefore no double loop
                 do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                     iorb = smat%keyg(1,2,iseg)
                     ii = ii + 1
                     matrix_local(ii) = matrixp(jorb,iorb-isfvctr)
                 end do
             end do
         end if
         if (ii/=nvctrp) stop 'compress_matrix_distributed: ii/=nvctrp'
     else
         stop 'compress_matrix_distributed: wrong data_strategy'
     end if

     call timing(iproc,'compressd_mcpy','OF')

     call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_PARALLEL, matrix_local, matrix_compr)
     call f_free_ptr(matrix_local)
     !@ END NEW #################

     call f_release_routine()

  end subroutine compress_matrix_distributed_wrapper_2



   subroutine compress_matrix_distributed_wrapper_1(iproc, nproc, smat, layout, matrixp, matrix_compr)
     use module_base
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(:),target,intent(inout) :: matrixp
     real(kind=8),dimension(smat%nvctrp_tg),target,intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr, nvctrp, ierr, isvctr
     integer :: ncount, itg, iitg, ist_send, ist_recv, i, iline, icolumn, ind
     integer :: window, sizeof, jproc_send, iorb, jproc, info
     integer,dimension(:),pointer :: isvctr_par, nvctr_par
     integer,dimension(:),allocatable :: request, windows
     real(kind=8),dimension(:),pointer :: matrix_local
     real(kind=8),dimension(:),allocatable :: recvbuf

     call f_routine(id='compress_matrix_distributed_wrapper_1')

     !call timing(iproc,'compressd_mcpy','ON')


     if (layout==SPARSE_MATMUL_SMALL) then
         if (size(matrixp)/=max(smat%smmm%nvctrp_mm,1)) then
             write(*,*) 'CRASH 1'
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(smat%smmm%nvctrp_mm,fmt='(i0)')), &
                  err_name='BIGDFT_RUNTIME_ERROR')
         end if
         matrix_local => matrixp
     else if (layout==SPARSE_MATMUL_LARGE) then
         if (size(matrixp)/=smat%smmm%nvctrp) then
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(smat%smmm%nvctrp,fmt='(i0)')), &
                  err_name='BIGDFT_RUNTIME_ERROR')
         end if
         matrix_local = f_malloc_ptr(max(1,smat%smmm%nvctrp_mm),id='matrix_local')
         call transform_sparsity_pattern(smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
              smat%nseg, smat%keyv, smat%keyg, smat%smmm%line_and_column_mm, &
              smat%smmm%nvctrp, smat%smmm%isvctr, &
              smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
              'large_to_small', matrix_local, matrixp)
         !call f_free_ptr(matrix_local)
     else
             call f_err_throw('layout has the value '//trim(yaml_toa(layout,fmt='(i0)'))//&
                  &'; allowed are '//trim(yaml_toa(SPARSE_MATMUL_SMALL,fmt='(i0)'))//&
                  &' and '//trim(yaml_toa(SPARSE_MATMUL_LARGE,fmt='(i0)')), &
                  err_name='BIGDFT_RUNTIME_ERROR')
     end if
     call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_MATMUL_SMALL, matrix_local, matrix_compr)
     if (layout==SPARSE_MATMUL_LARGE) then
         call f_free_ptr(matrix_local)
     end if

     !!call timing(iproc,'compressd_comm_new','OF')

     call f_release_routine()

  end subroutine compress_matrix_distributed_wrapper_1



   !> Gathers together the matrix parts calculated by other tasks.
   !! Attention: Even if the output array has size smat%nvctrp_tg, only the
   !! relevant part (indicated by smat%istartend_local) is calculated
   subroutine compress_matrix_distributed_core(iproc, nproc, smat, layout, matrixp, matrix_compr)
     use module_base
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout
     type(sparse_matrix),intent(in) :: smat
     !real(kind=8),dimension(smat%smmm%nvctrp_mm),intent(in) :: matrixp
     real(kind=8),dimension(:),intent(in) :: matrixp
     real(kind=8),dimension(smat%nvctrp_tg),target,intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr, nvctrp, ierr, isvctr
     integer :: ncount, itg, iitg, ist_send, ist_recv, i, iline, icolumn, ind
     integer :: window, sizeof, jproc_send, iorb, jproc, info, nccomm
     real(kind=8) :: window_fake
     integer,dimension(:),pointer :: isvctr_par, nvctr_par
     integer,dimension(:),allocatable :: request, windows
     real(kind=8),dimension(:),pointer :: matrix_local
     real(kind=8),dimension(:),allocatable :: recvbuf
     integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
     integer,parameter :: comm_strategy=GET
     integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
     integer,dimension(:,:),pointer :: luccomm

     call f_routine(id='compress_matrix_distributed_core')

     call timing(iproc,'compressd_mcpy','ON')


     if (data_strategy==GLOBAL_MATRIX) then
         stop 'compress_matrix_distributed: option GLOBAL_MATRIX is deprecated'
     else if (data_strategy==SUBMATRIX) then
         !if (layout==DENSE_PARALLEL) then
         if (layout==SPARSE_PARALLEL) then
             !stop 'layout==DENSE_PARALLEL not yet implemented'
             luccomm => smat%luccomm
             nvctrp = smat%nvctrp
             nccomm = smat%nccomm
         !else if (layout==DENSE_MATMUL) then
         else if (layout==SPARSE_MATMUL_SMALL) then
             luccomm => smat%smmm%luccomm_smmm
             nvctrp = smat%smmm%nvctrp_mm
             nccomm = smat%smmm%nccomm_smmm
         end if
         if (size(matrixp)/=max(1,nvctrp)) then
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(nvctrp,fmt='(i0)')), &
                  err_name='BIGDFT_RUNTIME_ERROR')
         end if

             call timing(iproc,'compressd_mcpy','OF')
             call timing(iproc,'compressd_comm','ON')

             if (nproc>1) then
                call f_zero(matrix_compr)

                 ! Create a window for all taskgroups to which iproc belongs (max 2)
                 windows = f_malloc(smat%ntaskgroup)
                 do itg=1,smat%ntaskgroupp
                     iitg = smat%taskgroupid(itg)
                     ! Use a fake window if nvctrp is zero
                     if (nvctrp>0) then
                         windows(iitg) = mpiwindow(nvctrp, matrixp(1), smat%mpi_groups(iitg)%mpi_comm)
                     else
                         windows(iitg) = mpiwindow(1, window_fake, smat%mpi_groups(iitg)%mpi_comm)
                     end if
                 end do
                 do jproc=1,nccomm
                     jproc_send = luccomm(1,jproc)
                     ist_send = luccomm(2,jproc)
                     ist_recv = luccomm(3,jproc)
                     ncount = luccomm(4,jproc)
                     !write(*,'(5(a,i0))') 'task ',iproc,' gets ',ncount,' elements at position ',ist_recv,' from position ',ist_send,' on task ',jproc_send
                     iitg = get_taskgroup_id(iproc,jproc_send)
                     ! Now get the task ID on the taskgroup (subtract the ID of the first task)
                     !jproc_send = jproc_send - smat%isrank(iitg)
                     ii = jproc_send
                     jproc_send = get_rank_on_taskgroup(ii,iitg)
                     !call mpiget(matrix_compr(ist_recv), ncount, jproc_send, int(ist_send-1,kind=mpi_address_kind), window)
                     !write(*,'(3(a,i0))') 'task ',iproc,' gets data from task ',jproc_send,' on window ',iitg
                     call mpiget(matrix_compr(ist_recv), ncount, jproc_send, int(ist_send-1,kind=mpi_address_kind), windows(iitg))
                 end do
             else
                 ist_send = luccomm(2,1)
                 ist_recv = luccomm(3,1)
                 ncount = luccomm(4,1)
                 call vcopy(ncount, matrixp(ist_send), 1, matrix_compr(ist_recv), 1)
             end if

             if (nproc>1) then
                 ! Synchronize the communication
                 do itg=1,smat%ntaskgroupp
                     iitg = smat%taskgroupid(itg)
                     call mpi_fenceandfree(windows(iitg))
                 end do
                 call f_free(windows)
                 !call mpi_fenceandfree(window)
             end if

             !!call f_free_ptr(matrix_local)

         !!end if

         call timing(iproc,'compressd_comm','OF')
     else
         stop 'compress_matrix_distributed: wrong data_strategy'
     end if

     call timing(iproc,'compressd_comm','OF')

     call f_release_routine()


     contains



       !> Get the taskgroup which should be used for the communication, i.e. the
       !! one to which both iproc and jproc belong
       integer function get_taskgroup_id(iproc,jproc)
         implicit none
         integer,intent(in) :: iproc, jproc

         ! Local variables
         integer :: itg, iitg, jtg, jjtg
         logical :: found

         ! A task never belongs to more than 2 taskgroups 
         found = .false.
         iloop: do itg=1,2
             iitg = smat%inwhichtaskgroup(itg,iproc)
             jloop: do jtg=1,2
                 jjtg = smat%inwhichtaskgroup(jtg,jproc)
                 if (iitg==jjtg) then
                     get_taskgroup_id = iitg
                     found = .true.
                     exit iloop
                 end if
             end do jloop
         end do iloop
         if (.not.found) stop 'get_taskgroup_id did not suceed'
       end function get_taskgroup_id


       ! Get the ID of task iiproc on taskgroup iitg
       integer function get_rank_on_taskgroup(iiproc,iitg)
         implicit none
         ! Calling arguments
         integer,intent(in) :: iiproc, iitg
         ! Local variables
         integer :: jproc
         logical :: found

         found = .false.
         do jproc=0,smat%nranks(iitg)-1
             if (smat%tgranks(jproc,iitg) == iiproc) then
                 get_rank_on_taskgroup = jproc
                 found = .true.
                 exit
             end if
         end do
         if (.not.found) stop 'get_rank_on_taskgroup did not suceed'
       end function get_rank_on_taskgroup

  end subroutine compress_matrix_distributed_core



  subroutine uncompress_matrix_distributed2(iproc, smat, layout, matrix_compr, matrixp)
    use module_base
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, layout
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctrp_tg),intent(in) :: matrix_compr
    real(kind=8),dimension(:,:),intent(out) :: matrixp

    ! Local variables
    integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr

      call timing(iproc,'compressd_mcpy','ON')

     ! Check the dimensions of the output array and assign some values
     if (size(matrixp,1)/=smat%nfvctr) stop 'size(matrixp,1)/=smat%nfvctr'
     if (layout==DENSE_PARALLEL) then
         if (size(matrixp,2)/=smat%nfvctrp) stop '(ubound(matrixp,2)/=smat%nfvctrp'
         nfvctrp=smat%nfvctrp
         isfvctr=smat%isfvctr
     else if (layout==DENSE_MATMUL) then
         if (size(matrixp,2)/=smat%smmm%nfvctrp) stop '(ubound(matrixp,2)/=smat%smmm%nfvctrp'
         nfvctrp=smat%smmm%nfvctrp
         isfvctr=smat%smmm%isfvctr
     end if

       if (nfvctrp>0) then

          !call to_zero(smat%nfvctr*nfvctrp,matrixp(1,1))
          call f_zero(matrixp)

           isegstart=smat%istsegline(isfvctr+1)
           isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
           !!isegstart=smat%istsegline(smat%isfvctr_par(iproc)+1)
           !!if (smat%isfvctr_par(iproc)+smat%nfvctrp<smat%nfvctr) then
           !!    isegend=smat%istsegline(smat%isfvctr_par(iproc+1)+1)-1
           !!else
           !!    isegend=smat%nseg
           !!end if
           !$omp parallel do default(private) &
           !$omp shared(isegstart, isegend, smat, matrixp, matrix_compr, isfvctr)
           do iseg=isegstart,isegend
               ii=smat%keyv(iseg)-1
               ! A segment is always on one line, therefore no double loop
               do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                   ii=ii+1
                   iiorb = smat%keyg(1,2,iseg)
                   jjorb = jorb
                   matrixp(jjorb,iiorb-isfvctr) = matrix_compr(ii-smat%isvctrp_tg)
               end do
           end do
           !$omp end parallel do
       end if

      call timing(iproc,'compressd_mcpy','OF')

   end subroutine uncompress_matrix_distributed2



   subroutine sequential_acces_matrix_fast(smat, a, a_seq)
     use module_base
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctr),intent(in) :: a
     real(kind=8),dimension(smat%smmm%nseq),intent(out) :: a_seq
   
     ! Local variables
     integer :: iseq, ii

     call f_routine(id='sequential_acces_matrix_fast')
   
     !$omp parallel do default(none) private(iseq, ii) &
     !$omp shared(smat, a_seq, a)
     do iseq=1,smat%smmm%nseq
         ii=smat%smmm%indices_extract_sequential(iseq)
         a_seq(iseq)=a(ii)
     end do
     !$omp end parallel do

     call f_release_routine()
   
   end subroutine sequential_acces_matrix_fast

   subroutine sequential_acces_matrix_fast2(smat, a, a_seq)
     use module_base
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctrp_tg),intent(in) :: a
     real(kind=8),dimension(smat%smmm%nseq),intent(out) :: a_seq
   
     ! Local variables
     integer :: iseq, ii

     call f_routine(id='sequential_acces_matrix_fast2')
   
     !$omp parallel do default(none) private(iseq, ii) &
     !$omp shared(smat, a_seq, a)
     do iseq=1,smat%smmm%nseq
         ii=smat%smmm%indices_extract_sequential(iseq)
         a_seq(iseq)=a(ii-smat%isvctrp_tg)
     end do
     !$omp end parallel do

     call f_release_routine()
   
   end subroutine sequential_acces_matrix_fast2




   subroutine sparsemm_new(smat, a_seq, b, c)
     use module_base
     use yaml_output
     implicit none
   
     !Calling Arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8), dimension(smat%smmm%nvctrp),intent(in) :: b
     real(kind=8), dimension(smat%smmm%nseq),intent(in) :: a_seq
     real(kind=8), dimension(smat%smmm%nvctrp), intent(out) :: c
   
     !Local variables
     !character(len=*), parameter :: subname='sparsemm'
     integer :: i,jorb,jjorb,m,mp1,ist,iend, icontiguous, j, iline, icolumn, nblock, iblock, ncount
     integer :: iorb, ii, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, iout
     real(kind=8) :: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7, ddot
     integer :: n_dense
     real(kind=8),dimension(:,:),allocatable :: a_dense, b_dense, c_dense
     !real(kind=8),dimension(:),allocatable :: b_dense, c_dense
     integer,parameter :: MATMUL_NEW = 101
     integer,parameter :: MATMUL_OLD = 102
     integer,parameter :: matmul_version = MATMUL_NEW
     logical,parameter :: count_flops = .false.
     real(kind=8) :: ts, te, op, gflops
     real(kind=8),parameter :: flop_per_op = 2.d0 !<number of FLOPS per operations
   
     call f_routine(id='sparsemm')
     if (count_flops) then
         n_dense = nint(sqrt(real(smat%smmm%nseq,kind=8)))
         !n_dense = smat%nfvctr
         a_dense = f_malloc((/n_dense,n_dense/),id='a_dense')
         b_dense = f_malloc((/n_dense,1/),id='b_dense')
         c_dense = f_malloc((/n_dense,1/),id='c_dense')
     end if
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'IR')


     if (matmul_version==MATMUL_NEW) then

         if (count_flops) then
             ! Start time
             ts = mpi_wtime()
         end if
         !$omp parallel default(private) shared(smat, a_seq, b, c)
         !$omp do schedule(guided)
         do iout=1,smat%smmm%nout
             i=smat%smmm%onedimindices_new(1,iout)
             nblock=smat%smmm%onedimindices_new(4,iout)
             tt0=0.d0

             do iblock=1,nblock
                 jorb = smat%smmm%consecutive_lookup(1,iblock,iout)
                 jjorb = smat%smmm%consecutive_lookup(2,iblock,iout)
                 ncount = smat%smmm%consecutive_lookup(3,iblock,iout)
                 !tt0 = tt0 + ddot(ncount, b(jjorb), 1, a_seq(jorb), 1)
                 !avoid calling ddot from OpenMP region on BG/Q as too expensive
                 tt0=tt0+my_dot(ncount,b(jjorb:jjorb+ncount-1),a_seq(jorb:jorb+ncount-1))
             end do

             c(i) = tt0
         end do 
         !$omp end do
         !$omp end parallel

         if (count_flops) then
             ! End time
             te = mpi_wtime()
             ! Count the operations
             op = 0.d0
             do iout=1,smat%smmm%nout
                 nblock=smat%smmm%onedimindices_new(4,iout)
                 do iblock=1,nblock
                     ncount = smat%smmm%consecutive_lookup(3,iblock,iout)
                     op = op + real(ncount,kind=8)
                 end do
             end do
             gflops = 1.d-9*op/(te-ts)*flop_per_op
             call yaml_map('SPARSE: operations',op,fmt='(es9.3)')
             call yaml_map('SPARSE: time',te-ts,fmt='(es9.3)')
             call yaml_map('SPARSE: GFLOPS',gflops)

             ! Compare with dgemm of comparable size
             ts = mpi_wtime()
             call dgemv('n', n_dense, n_dense, 1.d0, a_dense, n_dense, b_dense, 1, 0.d0, c_dense, 1)
             !call dgemm('n', 'n', n_dense, n_dense, n_dense, 1.d0, a_dense, n_dense, &
             !     b_dense, n_dense, 0.d0, c_dense, n_dense)
             te = mpi_wtime()
             op = real(n_dense,kind=8)*real(n_dense,kind=8)
             gflops = 1.d-9*op/(te-ts)*flop_per_op
             call yaml_map('DGEMV: operations',op,fmt='(es9.3)')
             call yaml_map('DGEMV: time',te-ts,fmt='(es9.3)')
             call yaml_map('DGEMV: GFLOPS',gflops)
         end if


     else if (matmul_version==MATMUL_OLD) then

         !$omp parallel default(private) shared(smat, a_seq, b, c)
         !$omp do
         do iout=1,smat%smmm%nout
             i=smat%smmm%onedimindices_new(1,iout)
             ilen=smat%smmm%onedimindices_new(2,iout)
             ii=smat%smmm%onedimindices_new(3,iout)
             tt0=0.d0

             iend=ii+ilen-1

             do jorb=ii,iend
                jjorb=smat%smmm%ivectorindex_new(jorb)
                tt0 = tt0 + b(jjorb)*a_seq(jorb)
             end do
             c(i) = tt0
         end do 
         !$omp end do
         !$omp end parallel


     else

         stop 'wrong value of matmul_version'

     end if

   
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'RS')
     call f_release_routine()

        contains

     pure function my_dot(n,x,y) result(tt)
       implicit none
       integer , intent(in) :: n
       double precision :: tt
       double precision, dimension(n), intent(in) :: x,y
       !local variables
       integer :: i

       tt=0.d0
       do i=1,n
          tt=tt+x(i)*y(i)
       end do
     end function
       
   end subroutine sparsemm_new


   !!function orb_from_index(smat, ival)
   !!  use sparsematrix_base, only: sparse_matrix
   !!  implicit none
   !!  ! Calling arguments
   !!  type(sparse_matrix),intent(in) :: smat
   !!  integer,intent(in) :: ival
   !!  integer,dimension(2) :: orb_from_index

   !!  orb_from_index(2) = (ival-1)/smat%nfvctr + 1
   !!  !orb_from_index(1) = ival - (orb_from_index_fn(2)-1)*smat%nfvctr
   !!  orb_from_index(1) = mod(ival-1,smat%nfvctr) + 1

   !!end function orb_from_index


   subroutine gather_matrix_from_taskgroups(iproc, nproc, smat, mat_tg, mat_global)
     use module_base
     use sparsematrix_base, only: sparse_matrix
     implicit none
   
     ! Calling arguments
     integer,intent(in) :: iproc, nproc
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: mat_tg !< matrix distributed over the taskgroups
     real(kind=8),dimension(smat%nvctr*smat%nspin),intent(out) :: mat_global !< global matrix gathered together
   
     ! Local variables
     integer,dimension(:),allocatable :: recvcounts, recvdspls
     integer :: ncount, ist_send, jproc, ispin, ishift

     call f_routine(id='gather_matrix_from_taskgroups')
   
     if (nproc>1) then
         recvcounts = f_malloc0(0.to.nproc-1,id='recvcounts')
         recvdspls = f_malloc0(0.to.nproc-1,id='recvdspls')
         !call to_zero(nproc, recvcounts(0))
         !call to_zero(nproc, recvdspls(0))
         ncount = smat%smmm%istartend_mm_dj(2) - smat%smmm%istartend_mm_dj(1) + 1
         recvcounts(iproc) = ncount
         call mpiallred(recvcounts(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
         recvdspls(0) = 0
         do jproc=1,nproc-1
             recvdspls(jproc) = recvdspls(jproc-1) + recvcounts(jproc-1)
         end do
         do ispin=1,smat%nspin
             ishift = (ispin-1)*smat%nvctr
             ist_send = smat%smmm%istartend_mm_dj(1) - smat%isvctrp_tg + ishift
             call mpi_get_to_allgatherv_double(mat_tg(ist_send), ncount, &
                  mat_global(ishift+1), recvcounts, recvdspls, bigdft_mpi%mpi_comm)
             !!call mpi_allgatherv(mat_tg(ist_send), ncount, mpi_double_precision, &
             !!                    mat_global(1), recvcounts, recvdspls, mpi_double_precision, &
             !!                    bigdft_mpi%mpi_comm, ierr)
         end do
         call f_free(recvcounts)
         call f_free(recvdspls)
     else
         call vcopy(smat%nvctrp*smat%nspin, mat_tg(1), 1, mat_global(1), 1)
     end if

     call f_release_routine()

   end subroutine gather_matrix_from_taskgroups


   subroutine gather_matrix_from_taskgroups_inplace(iproc, nproc, smat, mat)
     use module_base
     use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), SPARSE_FULL
     implicit none
   
     ! Calling arguments
     integer,intent(in) :: iproc, nproc
     type(sparse_matrix),intent(in) :: smat
     type(matrices),intent(inout) :: mat
   
     ! Local variables
     integer,dimension(:),allocatable :: recvcounts, recvdspls
     integer :: ncount, ist_send, jproc, ispin, ishift
     real(kind=8),dimension(:),allocatable :: mat_global
   
      mat_global = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='mat_global')
     if (nproc>1) then
         recvcounts = f_malloc0(0.to.nproc-1,id='recvcounts')
         recvdspls = f_malloc0(0.to.nproc-1,id='recvdspls')
         !call to_zero(nproc, recvcounts(0))
         !call to_zero(nproc, recvdspls(0))
         ncount = smat%smmm%istartend_mm_dj(2) - smat%smmm%istartend_mm_dj(1) + 1
         recvcounts(iproc) = ncount
         call mpiallred(recvcounts(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
         recvdspls(0) = 0
         do jproc=1,nproc-1
             recvdspls(jproc) = recvdspls(jproc-1) + recvcounts(jproc-1)
         end do
         do ispin=1,smat%nspin
             ishift = (ispin-1)*smat%nvctr
             ist_send = smat%smmm%istartend_mm_dj(1) - smat%isvctrp_tg + ishift
             call mpi_get_to_allgatherv_double(mat%matrix_compr(ist_send), ncount, &
                  mat_global(ishift+1), recvcounts, recvdspls, bigdft_mpi%mpi_comm)
             !!call mpi_allgatherv(mat%matrix_compr(ist_send), ncount, mpi_double_precision, &
             !!                    mat_global(1), recvcounts, recvdspls, mpi_double_precision, &
             !!                    bigdft_mpi%mpi_comm, ierr)
         end do
         call f_free(recvcounts)
         call f_free(recvdspls)
     else
         call vcopy(smat%nvctrp_tg*smat%nspin, mat%matrix_compr(1), 1, mat_global(1), 1)
     end if
     call vcopy(smat%nvctrp*smat%nspin, mat_global(1), 1, mat%matrix_compr(1), 1)
     call f_free(mat_global)

   end subroutine gather_matrix_from_taskgroups_inplace


   subroutine extract_taskgroup_inplace(smat, mat)
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     type(matrices),intent(inout) :: mat

     ! Local variables
     integer :: i, ispin, ishift_tg, ishift_glob

     do ispin=1,smat%nspin
         ishift_tg = (ispin-1)*smat%nvctrp_tg
         ishift_glob = (ispin-1)*smat%nvctr
         do i=1,smat%nvctrp_tg
             mat%matrix_compr(i+ishift_tg) = mat%matrix_compr(i+smat%isvctrp_tg+ishift_glob)
         end do
     end do

   end subroutine extract_taskgroup_inplace


   subroutine extract_taskgroup(smat, mat_glob, mat_tg)
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctr*smat%nspin),intent(in) :: mat_glob
     real(kind=8),dimension(smat%nvctrp_tg*smat%nspin),intent(out) :: mat_tg

     ! Local variables
     integer :: i, ispin, ishift_tg, ishift_glob

     do ispin=1,smat%nspin
         ishift_tg = (ispin-1)*smat%nvctrp_tg
         ishift_glob = (ispin-1)*smat%nvctr
         do i=1,smat%nvctrp_tg
             mat_tg(i+ishift_tg) = mat_glob(i+smat%isvctrp_tg+ishift_glob)
         end do
     end do

   end subroutine extract_taskgroup


    subroutine write_matrix_compressed(message, smat, mat)
      use yaml_output
      implicit none
    
      ! Calling arguments
      character(len=*),intent(in) :: message
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: iseg, i, ii, iorb, jorb
      integer,dimension(2) :: irowcol
    
      !!call yaml_sequence_open(trim(message))
      !!do iseg=1,smat%nseg
      !!    call yaml_sequence(advance='no')
      !!    ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
      !!    call yaml_mapping_open(flow=.true.)
      !!    call yaml_map('segment',iseg)
      !!    istart=smat%keyv(iseg)
      !!    iend=smat%keyv(iseg)+ilen
      !!    call yaml_map('values',smat%matrix_compr(istart:iend))
      !!    call yaml_mapping_close()
      !!    call yaml_newline()
      !!end do
      !!call yaml_sequence_close()
    
      call yaml_sequence_open(trim(message))
      do iseg=1,smat%nseg
          ! A segment is always on one line, therefore no double loop
          call yaml_sequence(advance='no')
          !ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
          call yaml_mapping_open(flow=.true.)
          call yaml_map('segment',iseg)
          call yaml_sequence_open('elements')
          !istart=smat%keyv(iseg)
          !iend=smat%keyv(iseg)+ilen-1
          !do i=istart,iend
          ii=smat%keyv(iseg)
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              !irowcol=orb_from_index(smat,i)
              !iorb=orb_from_index(1,i)
              !jorb=orb_from_index(2,i)
              call yaml_map('coordinates',(/smat%keyg(1,2,iseg),i/))
              call yaml_map('value',mat%matrix_compr(ii))
              call yaml_mapping_close()
              ii=ii+1
          end do
          call yaml_sequence_close()
          !call yaml_map('values',smat%matrix_compr(istart:iend))
          call yaml_mapping_close()
          call yaml_newline()
      end do
      call yaml_sequence_close()
    
    end subroutine write_matrix_compressed



   function check_symmetry(norb, smat)
     use module_base
     implicit none
   
     ! Calling arguments
     integer,intent(in) :: norb
     type(sparse_matrix),intent(in) :: smat
     logical :: check_symmetry
   
     ! Local variables
     integer :: i, iseg, ii, jorb, iorb
     logical,dimension(:,:),allocatable :: lgrid
     integer,dimension(2) :: irowcol
   
     lgrid=f_malloc((/norb,norb/),id='lgrid')
     lgrid=.false.
   
     do iseg=1,smat%nseg
         ii=smat%keyv(iseg)
         ! A segment is always on one line, therefore no double loop
         do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
             !irowcol=orb_from_index(smat,i)
             !!iorb=smat%orb_from_index(1,i)
             !!jorb=smat%orb_from_index(2,i)
             lgrid(smat%keyg(1,2,iseg),i)=.true.
             ii=ii+1
         end do
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


    !> Write a sparse matrix to a file
    subroutine write_sparsematrix(filename, smat, mat)
      use yaml_output
      implicit none
    
      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: iseg, i, ii
      integer,parameter :: iunit=234

      call f_routine(id='write_sparsematrix')

      ! First check that no taskgroups are used. Otherwise this routine does not work
      if (smat%ntaskgroup>1) then
          call f_err_throw('write_sparsematrix has not yet been implememted for matrix taskgroups', &
               err_name='BIGDFT_RUNTIME_ERROR')
      end if
    
      open(unit=iunit,file=filename)

      write(iunit,*) smat%nfvctr, '# number of columns'
      write(iunit,*) smat%nseg, '# number of segments'
      write(iunit,*) smat%nvctr, '# number of non-zero elements'
      do iseg=1,smat%nseg
          if(iseg==1) then
              write(iunit,*) smat%keyv(iseg), '# values of keyv'
          else
              write(iunit,*) smat%keyv(iseg)
          end if
      end do
      do iseg=1,smat%nseg
          if(iseg==1) then
              write(iunit,*) smat%keyg(1:2,1:2,iseg), '# values of keyg'
          else
              write(iunit,*) smat%keyg(1:2,1:2,iseg)
          end if
      end do
    
      do iseg=1,smat%nseg
          ! A segment is always on one line, therefore no double loop
          ii=smat%keyv(iseg)
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
              if (i==1 .and. iseg==1) then
                  write(iunit,*) mat%matrix_compr(ii), '# values of matrix_compr'
              else
                  write(iunit,*) mat%matrix_compr(ii)
              end if
              ii=ii+1
          end do
      end do

      close(unit=iunit)

      call f_release_routine()
    
    end subroutine write_sparsematrix



    !> Write a sparse matrix to a file, using the CCS format
    subroutine write_sparsematrix_CCS(filename, smat, mat)
      use yaml_output
      implicit none
    
      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: iseg, i, j, ii, icol, imat
      integer,dimension(:),allocatable :: col_ptr, row_ind, elements_per_column
      logical,dimension(:,:),allocatable :: matg
      logical :: column_started, first_in_column_set
      real(kind=8),dimension(:),allocatable :: val
      integer,parameter :: iunit=234, iunit2=235
      character(len=10) :: num
      character(len=100) :: frmt

      call f_routine(id='write_sparsematrix_CCS')

      col_ptr = f_malloc(smat%nfvctr,id='col_ptr')
      row_ind = f_malloc(smat%nvctr,id='row_ind')
      val = f_malloc(smat%nvctr,id='val')

      ii = 0
      do icol=1,smat%nfvctr
          imat = 0
          first_in_column_set = .false.
          do iseg=1,smat%nseg
              do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                  imat = imat + 1
                  if (i==icol) then
                      ! We are in column icol
                      ii = ii + 1
                      row_ind(ii) = smat%keyg(1,2,iseg) !row index
                      val(ii) = mat%matrix_compr(imat)
                      if (.not.first_in_column_set) then
                          col_ptr(icol) = ii
                          first_in_column_set = .true.
                      end if
                  end if
              end do
          end do
      end do
      if (ii/=smat%nvctr) stop 'ERROR in write_sparsematrix_CCS: ii/=smat%nvctr'

      !!matg = f_malloc((/smat%nfvctr,smat%nfvctr/),id='matg')
      !!matg = .false.
      !!do iseg=1,smat%nseg
      !!    ! A segment is always on one line, therefore no double loop
      !!    ii=smat%keyv(iseg)
      !!    do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
      !!        matg(smat%keyg(1,2,iseg),i) = .true.
      !!    end do
      !!end do

      !!ii = 0
      !!do i=1,smat%nfvctr
      !!    column_started = .false.
      !!    do j=1,smat%nfvctr
      !!       if(matg(j,i)) then
      !!           ii = ii + 1
      !!           row_ind(ii) = j
      !!           if (.not.column_started) then
      !!               col_ptr(i) = ii
      !!               column_started = .true.
      !!           end if
      !!       end if
      !!    end do
      !!end do
    
      open(unit=iunit,file=trim(filename))
      open(unit=iunit2,file=trim(filename)//'_2')

      write(iunit,*) smat%nfvctr, smat%nvctr, '# number of rows/columns, number of non-zero entries'
      write(iunit2,*) smat%nfvctr, smat%nfvctr, smat%nvctr
      do i=1,smat%nfvctr
          if (i==1) then
              write(iunit,*) col_ptr(i), '# col_ptr'
          else
              write(iunit,*) col_ptr(i)
          end if
      end do
      write(num,'(i0)') smat%nfvctr
      frmt='('//num//'(i0,1x))'
      write(iunit2,trim(frmt)) (col_ptr(i),i=1,smat%nfvctr)

      do i=1,smat%nvctr
          if (i==1) then
              write(iunit,*) row_ind(i), '# row_ind'
          else
              write(iunit,*) row_ind(i)
          end if
      end do
      write(num,'(i0)') smat%nvctr
      frmt='('//num//'(i0,1x))'
      write(iunit2,trim(frmt)) (row_ind(i),i=1,smat%nvctr)
      
      do i=1,smat%nvctr
          if(i==1) then
              !!write(iunit,*) mat%matrix_compr(i), '# values of matrix_compr' 
              !!write(iunit2,*) mat%matrix_compr(i)
              write(iunit,*) val(i), '# values of matrix_compr' 
              write(iunit2,*) val(i)
          else
              !write(iunit,*) mat%matrix_compr(i) 
              !write(iunit2,*) mat%matrix_compr(i) 
              write(iunit,*) val(i) 
              write(iunit2,*) val(i) 
          end if
      end do
      !write(num,'(i0)') smat%nvctr
      !frmt='('//num//'i9)'
      !write(iunit2,trim(frmt)) (mat%matrix_compr(i),i=1,smat%nvctr)

      close(unit=iunit)
      close(unit=iunit2)

      call f_free(col_ptr)
      call f_free(row_ind)
      call f_free(val)
      !!call f_free(matg)

      call f_release_routine()
    
    end subroutine write_sparsematrix_CCS


    !> Transform a matrix from a large parsity pattern *_l to a small sparsity pattern *_s or vice versa.
    !! The small pattern must be contained within the large one.
    subroutine transform_sparsity_pattern(nfvctr, nvctrp_s, isvctr_s, nseg_s, keyv_s, keyg_s, line_and_column_s, &
               nvctrp_l, isvctr_l, nseg_l, keyv_l, keyg_l, istsegline_l, direction, matrix_s, matrix_l)
      use sparsematrix_init, only: matrixindex_in_compressed_lowlevel
      implicit none
      ! Calling arguments
      integer,intent(in) :: nfvctr, nvctrp_s, isvctr_s, nseg_s, nvctrp_l, isvctr_l, nseg_l
      integer,dimension(2,nvctrp_s),intent(in) :: line_and_column_s
      integer,dimension(nseg_s),intent(in) :: keyv_s
      integer,dimension(2,2,nseg_s),intent(in) :: keyg_s
      integer,dimension(nseg_l),intent(in) :: keyv_l
      integer,dimension(2,2,nseg_l),intent(in) :: keyg_l
      integer,dimension(nfvctr),intent(in) :: istsegline_l
      character(len=*),intent(in) :: direction
      real(kind=8),dimension(nvctrp_l),intent(inout) :: matrix_l
      real(kind=8),dimension(nvctrp_s),intent(inout) :: matrix_s
      ! Local variables
      integer :: i, ii, ind, iline, icolumn

      call f_routine(id='transform_sparsity_pattern')
      call timing(bigdft_mpi%iproc, 'transformspars', 'ON')

        if (direction=='large_to_small') then

            ! No need for f_zero since every value will be overwritten.
            !$omp parallel default(none) &
            !$omp shared(nvctrp_s, isvctr_s, isvctr_l, line_and_column_s) &
            !$omp shared(nfvctr, nseg_l, keyv_l, keyg_l, istsegline_l, matrix_s, matrix_l) &
            !$omp private(i, ii, iline, icolumn, ind)
            !$omp do
            do i=1,nvctrp_s
                ii = isvctr_s + i
                !!call get_line_and_column(ii, nseg_s, keyv_s, keyg_s, iline, icolumn)
                iline = line_and_column_s(1,i)
                icolumn = line_and_column_s(2,i)
                ind = matrixindex_in_compressed_lowlevel(icolumn, iline, nfvctr, &
                      nseg_l, keyv_l, keyg_l, istsegline_l)
                ind = ind - isvctr_l
                matrix_s(i) = matrix_l(ind)
            end do
            !$omp end do
            !$omp end parallel

        else if (direction=='small_to_large') then
            call f_zero(matrix_l)
            !$omp parallel default(none) &
            !$omp shared(nvctrp_s, isvctr_s, isvctr_l, line_and_column_s) &
            !$omp shared(nfvctr, nseg_l, keyv_l, keyg_l, istsegline_l, matrix_s, matrix_l) &
            !$omp private(i, ii, iline, icolumn, ind)
            !$omp do
            do i=1,nvctrp_s
                ii = isvctr_s + i
                !call get_line_and_column(ii, nseg_s, keyv_s, keyg_s, iline, icolumn)
                iline = line_and_column_s(1,i)
                icolumn = line_and_column_s(2,i)
                ind = matrixindex_in_compressed_lowlevel(icolumn, iline, nfvctr, &
                      nseg_l, keyv_l, keyg_l, istsegline_l)
                ind = ind - isvctr_l
                matrix_l(ind) = matrix_s(i)
            end do
            !$omp end do
            !$omp end parallel
        else
            stop 'wrong direction'
        end if

      call timing(bigdft_mpi%iproc, 'transformspars', 'OF')
      call f_release_routine()

    end subroutine transform_sparsity_pattern



    !> Calculates c = a*b for matrices a,b,c
    subroutine matrix_matrix_mult_wrapper(iproc, nproc, smat, a, b, c)
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(smat%nvctrp_tg),intent(in) :: a
      real(kind=8),dimension(smat%nvctrp_tg),intent(inout) :: b, c

      ! Local variables
      real(kind=8),dimension(:),allocatable :: b_exp, c_exp, a_seq

      b_exp = f_malloc(smat%smmm%nvctrp, id='b_exp')
      c_exp = f_malloc(smat%smmm%nvctrp, id='c_exp')
      a_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='a_seq')

      call sequential_acces_matrix_fast2(smat, a, a_seq)
      if (smat%smmm%nvctrp_mm>0) then !to avoid out of bounds error...
          call transform_sparsity_pattern(smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
               smat%nseg, smat%keyv, smat%keyg, &
               smat%smmm%line_and_column_mm, &
               smat%smmm%nvctrp, smat%smmm%isvctr, &
               smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
               'small_to_large', b(smat%smmm%isvctr_mm-smat%isvctrp_tg+1), b_exp)
      end if
      call sparsemm_new(smat, a_seq, b_exp, c_exp)
      call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_LARGE, &
           c_exp, c)

      call f_free(b_exp)
      call f_free(c_exp)
      call f_free(a_seq)


    end subroutine matrix_matrix_mult_wrapper


    !< Calculates the trace of the matrix product amat*bmat.
    !< WARNING: It is mandatory that the sparsity pattern of amat be contained
    !< within the sparsity pattern of bmat!
    function trace_sparse(iproc, nproc, orbs, asmat, bsmat, amat, bmat, ispin)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc,  nproc, ispin
      type(orbitals_data),intent(in) :: orbs
      type(sparse_matrix),intent(in) :: asmat, bsmat
      real(kind=8),dimension(asmat%nvctrp_tg),intent(in) :: amat
      real(kind=8),dimension(bsmat%nvctrp_tg),intent(in) :: bmat
    
      ! Local variables
      integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iilarge
      integer :: ierr, iashift, ibshift, iel
      real(kind=8) :: sumn, trace_sparse
    
    
      call f_routine(id='trace_sparse')
    
      iashift = 0!(ispin-1)*asmat%nvctr
      ibshift = 0!(ispin-1)*bsmat%nvctr
    
    
      sumn=0.d0
      !if (asmat%smmm%nfvctrp>0) then
          !$omp parallel default(none) &
          !$omp private(iseg, ii, jorb, iiorb, jjorb, iilarge, iel) &
          !$omp shared(bsmat, asmat, amat, bmat, iashift, ibshift, sumn)
          !$omp do reduction(+:sumn)
          !do iseg=isegstart,isegend
          do iseg=asmat%smmm%isseg,asmat%smmm%ieseg
              iel = asmat%keyv(iseg) - 1
              ii=iashift+asmat%keyv(iseg)-1
              ! A segment is always on one line, therefore no double loop
              do jorb=asmat%keyg(1,1,iseg),asmat%keyg(2,1,iseg)
                  iel = iel + 1
                  if (iel<asmat%smmm%isvctr_mm+1) cycle
                  if (iel>asmat%smmm%isvctr_mm+asmat%smmm%nvctrp_mm) then
                      !write(*,*) 'exit with iel',iel
                      exit
                  end if
                  ii=ii+1
                  iiorb = asmat%keyg(1,2,iseg)
                  jjorb = jorb
                  iilarge = ibshift + matrixindex_in_compressed(bsmat, iiorb, jjorb)
                  !!write(*,'(a,4i8,3es16.8)') 'iproc, ii, iilarge, iend, vals, sumn', &
                  !!    iproc, ii, iilarge, asmat%smmm%isvctr_mm+asmat%smmm%nvctrp_mm, amat(ii-asmat%isvctrp_tg), bmat(iilarge-bsmat%isvctrp_tg), sumn
                  sumn = sumn + amat(ii-asmat%isvctrp_tg)*bmat(iilarge-bsmat%isvctrp_tg)
              end do  
          end do
          !$omp end do
          !$omp end parallel
      !end if
    
      if (nproc > 1) then
          call mpiallred(sumn, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      trace_sparse = sumn
    
      call f_release_routine()
    
    end function trace_sparse

end module sparsematrix

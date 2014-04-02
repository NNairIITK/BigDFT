module sparsematrix
  use module_base
  use sparsematrix_base
  implicit none

  private

  !> Public routines
  public :: compress_matrix
  public :: uncompress_matrix
  public :: check_matrix_compression
  public :: transform_sparse_matrix

  contains

    !> subroutine to compress the matrix to sparse form
    subroutine compress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr),target,intent(in),optional :: inmat
      real(kind=8),dimension(sparsemat%nvctr),target,intent(out),optional :: outmat
    
      ! Local variables
      integer :: jj, irow, jcol, jjj, ierr
      real(kind=8),dimension(:,:),pointer :: inm
      real(kind=8),dimension(:),pointer :: outm

      if (present(outmat)) then
          if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
              stop 'outmat not allowed for the given options'
          end if
          outm => outmat
      else
          outm => sparsemat%matrix_compr
      end if

      if (present(inmat)) then
          if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
              stop 'in not allowed for the given options'
          end if
          inm => inmat
      else
          inm => sparsemat%matrix
      end if
    
      call timing(iproc,'compress_uncom','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         !$omp parallel do default(private) shared(sparsemat,inm,outm)
         do jj=1,sparsemat%nvctr
            irow = sparsemat%orb_from_index(1,jj)
            jcol = sparsemat%orb_from_index(2,jj)
            outm(jj)=inm(irow,jcol)
         end do
         !$omp end parallel do
      else if (sparsemat%parallel_compression==1) then
         call to_zero(sparsemat%nvctr, sparsemat%matrix_compr(1))
         !$omp parallel do default(private) shared(sparsemat)
         do jj=1,sparsemat%nvctrp
            jjj=jj+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,jjj)
            jcol = sparsemat%orb_from_index(2,jjj)
            sparsemat%matrix_compr(jjj)=sparsemat%matrix(irow,jcol)
         end do
         !$omp end parallel do
         call mpiallred(sparsemat%matrix_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
         sparsemat%matrix_comprp=f_malloc_ptr((sparsemat%nvctrp),id='sparsemat%matrix_comprp')
         !$omp parallel do default(private) shared(sparsemat)
         do jj=1,sparsemat%nvctrp
            jjj=jj+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,jjj)
            jcol = sparsemat%orb_from_index(2,jjj)
            sparsemat%matrix_comprp(jj)=sparsemat%matrix(irow,jcol)
         end do
         !$omp end parallel do
         call mpi_allgatherv(sparsemat%matrix_comprp, sparsemat%nvctrp, mpi_double_precision, sparsemat%matrix_compr, &
              sparsemat%nvctr_par(:), sparsemat%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         call f_free_ptr(sparsemat%matrix_comprp)
      end if
    
      call timing(iproc,'compress_uncom','OF')
    
    end subroutine compress_matrix



    !> subroutine to uncompress the matrix from sparse form
    subroutine uncompress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix), intent(inout) :: sparsemat
      real(kind=8),dimension(sparsemat%nvctr),target,intent(out),optional :: inmat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr),target,intent(out),optional :: outmat
      
      ! Local variables
      integer :: ii, irow, jcol, iii, ierr
      real(kind=8),dimension(:),pointer :: inm
      real(kind=8),dimension(:,:),pointer :: outm

      if (present(outmat)) then
          if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
              stop 'outmat not allowed for the given options'
          end if
          outm => outmat
      else
          outm => sparsemat%matrix
      end if

      if (present(inmat)) then
          if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
              stop 'inmat not allowed for the given options'
          end if
          inm => inmat
      else
          inm => sparsemat%matrix_compr
      end if
    
      call timing(iproc,'compress_uncom','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         call to_zero(sparsemat%nfvctr**2, outm(1,1))
         !$omp parallel do default(private) shared(sparsemat,inm,outm)
         do ii=1,sparsemat%nvctr
            irow = sparsemat%orb_from_index(1,ii)
            jcol = sparsemat%orb_from_index(2,ii)
            outm(irow,jcol)=inm(ii)
         end do
         !$omp end parallel do
      else if (sparsemat%parallel_compression==1) then
         call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
         !$omp parallel do default(private) shared(sparsemat)
         do ii=1,sparsemat%nvctrp
            iii=ii+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,iii)
            jcol = sparsemat%orb_from_index(2,iii)
            sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(iii)
         end do
         !$omp end parallel do
         call mpiallred(sparsemat%matrix(1,1), sparsemat%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
         sparsemat%matrixp=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
         call to_zero(sparsemat%nfvctr*sparsemat%nfvctrp, sparsemat%matrixp(1,1))
         !$omp parallel do default(private) shared(sparsemat)
         do ii=1,sparsemat%nvctrp
            iii=ii+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,iii)
            jcol = sparsemat%orb_from_index(2,iii) - sparsemat%isfvctr
            sparsemat%matrixp(irow,jcol)=sparsemat%matrix_compr(iii)
         end do
         !$omp end parallel do
         call mpi_allgatherv(sparsemat%matrixp, sparsemat%nfvctr*sparsemat%nfvctrp, mpi_double_precision, sparsemat%matrix, &
              sparsemat%nfvctr*sparsemat%nfvctr_par(:), sparsemat%nfvctr*sparsemat%isfvctr_par, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         call f_free_ptr(sparsemat%matrixp)
      end if
      sparsemat%can_use_dense=.true.  
    
      call timing(iproc,'compress_uncom','OF')
    
    end subroutine uncompress_matrix



    subroutine check_matrix_compression(iproc,sparsemat)
      use yaml_output
      implicit none
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
      !Local variables
      integer :: i_stat, i_all, jorb, irow, icol, iseg, ii
      character(len=*),parameter :: subname='check_matrix_compression'
      real(kind=8) :: maxdiff
      real(kind=8), parameter :: tol=1.e-10
    
    
      sparsemat%matrix=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctr/),id='sparsemat%matrix')
      sparsemat%matrix_compr=f_malloc_ptr(sparsemat%nvctr,id='sparsemat%matrix_compr')
    
      call to_zero(sparsemat%nfvctr**2,sparsemat%matrix(1,1))
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
            !print *,'irow,icol',irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
            sparsemat%matrix(irow,icol) = test_value_matrix(sparsemat%nfvctr, irow, icol)
         end do
      end do
      
      call compress_matrix(iproc,sparsemat)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ii=0
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
            maxdiff = max(abs(sparsemat%matrix_compr(sparsemat%keyv(iseg)+ii)&
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
    
      call uncompress_matrix(iproc,sparsemat)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
            maxdiff = max(abs(sparsemat%matrix(irow,icol)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
         end do
      end do
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('UNCOMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for uncompress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call f_free_ptr(sparsemat%matrix)
      call f_free_ptr(sparsemat%matrix_compr)
    
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
      type(sparse_matrix),intent(inout) :: smat, lmat
      real(kind=8),dimension(smat%nvctr),intent(inout) :: smatrix_compr
      real(kind=8),dimension(lmat%nvctr),intent(inout) :: lmatrix_compr
      character(len=14),intent(in) :: cmode
    
      ! Local variables
      integer :: imode, icheck, isseg, isstart, isend, ilseg, ilstart, ilend
      integer :: iostart, ioend, ilength, isoffset, iloffset, iscostart, ilcostart, i
      integer,parameter :: SMALL_TO_LARGE=1
      integer,parameter :: LARGE_TO_SMALL=2
    
    
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
          call to_zero(lmat%nvctr,lmatrix_compr(1))
      case (LARGE_TO_SMALL)
          call to_zero(smat%nvctr,smatrix_compr(1))
      case default
          stop 'wrong imode'
      end select
    
    
      icheck=0
      sloop: do isseg=1,smat%nseg
          isstart=smat%keyg(1,isseg)
          isend=smat%keyg(2,isseg)
          lloop: do ilseg=1,lmat%nseg
              ilstart=lmat%keyg(1,ilseg)
              ilend=lmat%keyg(2,ilseg)
    
              !write(*,*) 'isstart, isend, ilstart, ilend', isstart, isend, ilstart, ilend
              ! check whether there is an overlap:
              ! if not, increase loop counters
              if (ilstart>isend) exit lloop
              if (isstart>ilend) cycle lloop
              ! if yes, determine start end end of overlapping segment (in uncompressed form)
              iostart=max(isstart,ilstart)
              ioend=min(isend,ilend)
              !write(*,*) 'iostart, ioend', iostart, ioend
              ilength=ioend-iostart+1
    
              ! offset with respect to the starting point of the segment
              isoffset=iostart-smat%keyg(1,isseg)
              iloffset=iostart-lmat%keyg(1,ilseg)
    
              ! determine start end and of the overlapping segment in compressed form
              iscostart=smat%keyv(isseg)+isoffset
              ilcostart=lmat%keyv(ilseg)+iloffset
    
              ! copy the elements
              select case (imode)
              case (SMALL_TO_LARGE) 
                  do i=0,ilength-1
                      lmatrix_compr(ilcostart+i)=smatrix_compr(iscostart+i)
                  end do
              case (LARGE_TO_SMALL) 
                  do i=0,ilength-1
                      smatrix_compr(iscostart+i)=lmatrix_compr(ilcostart+i)
                  end do
              case default
                  stop 'wrong imode'
              end select
              icheck=icheck+ilength
          end do lloop
      end do sloop
    
      ! all elements of the small matrix must have been processed, no matter in
      ! which direction the transformation has been executed
      if (icheck/=smat%nvctr) then
          write(*,'(a,2i8)') 'ERROR: icheck/=smat%nvctr', icheck, smat%nvctr
          stop
      end if
    
    end subroutine transform_sparse_matrix

end module sparsematrix

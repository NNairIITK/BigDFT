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
  public :: compress_matrix
  public :: uncompress_matrix
  public :: check_matrix_compression
  public :: transform_sparse_matrix
  public :: compress_matrix_distributed
  public :: uncompress_matrix_distributed
  public :: sequential_acces_matrix_fast
  public :: sparsemm

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
         if (bigdft_mpi%nproc > 1) then
            call mpiallred(sparsemat%matrix_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
         end if
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
         if (bigdft_mpi%nproc > 1) then
            call mpiallred(sparsemat%matrix(1,1), sparsemat%nfvctr**2,mpi_sum,bigdft_mpi%mpi_comm)
         end if
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
    
      mat%matrix = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_FULL, id='mat%matrix')
    
      call to_zero(sparsemat%nfvctr**2,mat%matrix(1,1))
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
            !print *,'irow,icol',irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
            mat%matrix(irow,icol) = test_value_matrix(sparsemat%nfvctr, irow, icol)
         end do
      end do
      
      call compress_matrix(iproc, sparsemat, inmat=mat%matrix, outmat=mat%matrix_compr)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ii=0
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
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
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
            maxdiff = max(abs(mat%matrix(irow,icol)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
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
      type(sparse_matrix),intent(inout) :: smat, lmat
      real(kind=8),dimension(smat%nvctr),intent(inout) :: smatrix_compr
      real(kind=8),dimension(lmat%nvctr),intent(inout) :: lmatrix_compr
      character(len=14),intent(in) :: cmode
    
      ! Local variables
      integer :: imode, icheck, isseg, isstart, isend, ilseg, ilstart, ilend
      integer :: iostart, ioend, ilength, isoffset, iloffset, iscostart, ilcostart, i
      integer :: ilsegstart
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
    
      call timing(bigdft_mpi%iproc,'transform_matr','IR')
    
      icheck=0
      ilsegstart=1
      !$omp parallel default(private) &
      !$omp shared(smat, lmat, imode, lmatrix_compr, smatrix_compr, icheck) &
      !$omp firstprivate(ilsegstart)
      !$omp do reduction(+:icheck)
      sloop: do isseg=1,smat%nseg
          isstart=smat%keyg(1,isseg)
          isend=smat%keyg(2,isseg)
          lloop: do ilseg=ilsegstart,lmat%nseg
              ilstart=lmat%keyg(1,ilseg)
              ilend=lmat%keyg(2,ilseg)
    
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
      !$omp end do 
      !$omp end parallel
    
      ! all elements of the small matrix must have been processed, no matter in
      ! which direction the transformation has been executed
      if (icheck/=smat%nvctr) then
          write(*,'(a,2i8)') 'ERROR: icheck/=smat%nvctr', icheck, smat%nvctr
          stop
      end if

      call timing(bigdft_mpi%iproc,'transform_matr','RS')
    
    end subroutine transform_sparse_matrix




   subroutine compress_matrix_distributed(iproc, smat, matrixp, matrix_compr)
     use module_base
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nfvctr,smat%nfvctrp),intent(in) :: matrixp
     real(kind=8),dimension(smat%nvctr),intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, ierr


     call to_zero(smat%nvctr, matrix_compr(1))

     if (smat%nfvctrp>0) then
         isegstart=smat%istsegline(smat%isfvctr_par(iproc)+1)
         if (smat%isfvctr_par(iproc)+smat%nfvctrp<smat%nfvctr) then
             isegend=smat%istsegline(smat%isfvctr_par(iproc+1)+1)-1
         else
             isegend=smat%nseg
         end if
         !$omp parallel default(none) &
         !$omp shared(isegstart, isegend, matrixp, smat, matrix_compr,iproc) &
         !$omp private(iseg, ii, jorb, iiorb, jjorb)
         !$omp do
         do iseg=isegstart,isegend
             ii=smat%keyv(iseg)-1
             do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                 ii=ii+1
                 iiorb = (jorb-1)/smat%nfvctr + 1
                 jjorb = jorb - (iiorb-1)*smat%nfvctr
                 matrix_compr(ii)=matrixp(jjorb,iiorb-smat%isfvctr_par(iproc))
             end do
         end do
         !$omp end do
         !$omp end parallel
     end if

     if (bigdft_mpi%nproc>1) then
         call mpiallred(matrix_compr(1), smat%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
     end if

  end subroutine compress_matrix_distributed


  subroutine uncompress_matrix_distributed(iproc, smat, matrix_compr, matrixp)
    use module_base
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctr),intent(in) :: matrix_compr
    real(kind=8),dimension(smat%nfvctr,smat%nfvctrp),intent(out) :: matrixp

    ! Local variables
    integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


       if (smat%nfvctrp>0) then

           call to_zero(smat%nfvctr*smat%nfvctrp,matrixp(1,1))

           isegstart=smat%istsegline(smat%isfvctr_par(iproc)+1)
           if (smat%isfvctr_par(iproc)+smat%nfvctrp<smat%nfvctr) then
               isegend=smat%istsegline(smat%isfvctr_par(iproc+1)+1)-1
           else
               isegend=smat%nseg
           end if
           !$omp parallel do default(private) &
           !$omp shared(isegstart, isegend, smat, matrixp, matrix_compr, iproc)
           do iseg=isegstart,isegend
               ii=smat%keyv(iseg)-1
               do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                   ii=ii+1
                   iiorb = (jorb-1)/smat%nfvctr + 1
                   jjorb = jorb - (iiorb-1)*smat%nfvctr
                   matrixp(jjorb,iiorb-smat%isfvctr_par(iproc)) = matrix_compr(ii)
               end do
           end do
           !$omp end parallel do
       end if

   end subroutine uncompress_matrix_distributed

   subroutine sequential_acces_matrix_fast(smat, a, a_seq)
     use module_base
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(smat%nvctr),intent(in) :: a
     real(kind=8),dimension(smat%smmm%nseq),intent(out) :: a_seq
   
     ! Local variables
     integer :: iseq, ii
   
     !$omp parallel do default(none) private(iseq, ii) &
     !$omp shared(smat, a_seq, a)
     do iseq=1,smat%smmm%nseq
         ii=smat%smmm%indices_extract_sequential(iseq)
         a_seq(iseq)=a(ii)
     end do
     !$omp end parallel do
   
   end subroutine sequential_acces_matrix_fast


   subroutine sparsemm(smat, a_seq, b, c)
     use module_base
     implicit none
   
     !Calling Arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8), dimension(smat%nfvctr,smat%nfvctrp),intent(in) :: b
     real(kind=8), dimension(smat%smmm%nseq),intent(in) :: a_seq
     real(kind=8), dimension(smat%nfvctr,smat%nfvctrp), intent(out) :: c
   
     !Local variables
     !character(len=*), parameter :: subname='sparsemm'
     integer :: i,jorb,jjorb,m,mp1
     integer :: iorb, ii0, ii2, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, iout
     real(kind=8) :: tt
   
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'IR')

   
     !$omp parallel default(private) shared(smat, a_seq, b, c)
     !$omp do
     do iout=1,smat%smmm%nout
         i=smat%smmm%onedimindices(1,iout)
         iorb=smat%smmm%onedimindices(2,iout)
         ilen=smat%smmm%onedimindices(3,iout)
         ii0=smat%smmm%onedimindices(4,iout)
         ii2=0
         tt=0.d0
   
         m=mod(ilen,7)
         if (m/=0) then
             do jorb=1,m
                jjorb=smat%smmm%ivectorindex(ii0+ii2)
                tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
                ii2=ii2+1
             end do
         end if
         mp1=m+1
         do jorb=mp1,ilen,7
   
            jjorb0=smat%smmm%ivectorindex(ii0+ii2+0)
            tt = tt + b(jjorb0,i)*a_seq(ii0+ii2+0)
   
            jjorb1=smat%smmm%ivectorindex(ii0+ii2+1)
            tt = tt + b(jjorb1,i)*a_seq(ii0+ii2+1)
   
            jjorb2=smat%smmm%ivectorindex(ii0+ii2+2)
            tt = tt + b(jjorb2,i)*a_seq(ii0+ii2+2)
   
            jjorb3=smat%smmm%ivectorindex(ii0+ii2+3)
            tt = tt + b(jjorb3,i)*a_seq(ii0+ii2+3)
   
            jjorb4=smat%smmm%ivectorindex(ii0+ii2+4)
            tt = tt + b(jjorb4,i)*a_seq(ii0+ii2+4)
   
            jjorb5=smat%smmm%ivectorindex(ii0+ii2+5)
            tt = tt + b(jjorb5,i)*a_seq(ii0+ii2+5)
   
            jjorb6=smat%smmm%ivectorindex(ii0+ii2+6)
            tt = tt + b(jjorb6,i)*a_seq(ii0+ii2+6)
   
            ii2=ii2+7
         end do
         c(iorb,i)=tt
     end do 
     !$omp end do
     !$omp end parallel
   
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'RS')
       
   end subroutine sparsemm



end module sparsematrix

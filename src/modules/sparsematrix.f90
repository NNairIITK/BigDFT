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
  public :: orb_from_index

  contains

    !> subroutine to compress the matrix to sparse form
    subroutine compress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
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
             ishift=(ispin-1)*sparsemat%nfvctr**2
             !$omp parallel default(none) private(iseg,j,jj,irowcol) &
             !$omp shared(sparsemat,inm,outm,ishift,ispin)
             !$omp do
             do iseg=1,sparsemat%nseg
                 jj=sparsemat%keyv(iseg)
                 do j=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                    !irow = sparsemat%orb_from_index(1,jj)
                    !jcol = sparsemat%orb_from_index(2,jj)
                    irowcol = orb_from_index(sparsemat, j)
                    !write(*,*) 'iseg, j, jj', iseg, j, jj
                    outm(jj+ishift)=inm(irowcol(1),irowcol(2),ispin)
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



    !> subroutine to uncompress the matrix from sparse form
    subroutine uncompress_matrix(iproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix), intent(inout) :: sparsemat
      real(kind=8),dimension(sparsemat%nvctr*sparsemat%nspin),target,intent(out) :: inmat
      real(kind=8),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),target,intent(out) :: outmat
      
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
         call to_zero(sparsemat%nfvctr**2*sparsemat%nspin, outm(1,1,1))
         do ispin=1,sparsemat%nspin
             ishift=(ispin-1)*sparsemat%nvctr
             !$omp parallel default(none) private(iseg,i,ii,irowcol) shared(sparsemat,inm,outm,ispin,ishift)
             !$omp do
             do iseg=1,sparsemat%nseg
                 ii=sparsemat%keyv(iseg)
                 do i=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                    !irow = sparsemat%orb_from_index(1,ii)
                    !jcol = sparsemat%orb_from_index(2,ii)
                    irowcol = orb_from_index(sparsemat, i)
                    outm(irowcol(1),irowcol(2),ispin)=inm(ii+ishift)
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
    
      call to_zero(sparsemat%nfvctr**2*sparsemat%nspin,mat%matrix(1,1,1))
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
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
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
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
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indices(jorb,irow,icol)
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
      real(kind=8),dimension(smat%nspin*smat%nvctr),intent(inout) :: smatrix_compr
      real(kind=8),dimension(lmat%nspin*lmat%nvctr),intent(inout) :: lmatrix_compr
      character(len=14),intent(in) :: cmode
    
      ! Local variables
      integer :: imode, icheck, isseg, isstart, isend, ilseg, ilstart, ilend
      integer :: iostart, ioend, ilength, isoffset, iloffset, iscostart, ilcostart, i
      integer :: ilsegstart, ispin, isshift, ilshift
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
          call to_zero(lmat%nvctr*lmat%nspin,lmatrix_compr(1))
      case (LARGE_TO_SMALL)
          call to_zero(smat%nvctr*lmat%nspin,smatrix_compr(1))
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
    
    end subroutine transform_sparse_matrix


   subroutine compress_matrix_distributed(iproc, nproc, smat, layout, matrixp, matrix_compr)
     use module_base
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout
     type(sparse_matrix),intent(in) :: smat
     real(kind=8),dimension(:,:),intent(in) :: matrixp
     real(kind=8),dimension(smat%nvctr),target,intent(out) :: matrix_compr

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr, nvctrp, ierr, isvctr
     integer,dimension(:),pointer :: isvctr_par, nvctr_par
     real(kind=8),dimension(:),pointer :: matrix_local
     integer,parameter :: ALLGATHERV=51, GET=52
     integer,parameter :: comm_strategy=GET

     call f_routine(id='compress_matrix_distributed')

     call timing(iproc,'compressd_mcpy','ON')

     ! Check the dimensions of the input array and assign some values
     if (size(matrixp,1)/=smat%nfvctr) stop 'size(matrixp,1)/=smat%nfvctr'
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
         nvctrp = smat%smmM%nvctrp
         isvctr = smat%smmm%isvctr
         isvctr_par => smat%smmm%isvctr_par
         nvctr_par => smat%smmm%nvctr_par
     end if

     !call to_zero(smat%nvctr, matrix_compr(1))
     if (nproc>1) then
         matrix_local = f_malloc0_ptr(nvctrp,id='matrix_local')
     else
         matrix_local => matrix_compr
     end if

     if (nfvctrp>0) then
         isegstart=smat%istsegline(isfvctr+1)
         isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
         !!if (isfvctr+nfvctrp<smat%nfvctr) then
         !!    isegend=smat%istsegline(smat%isfvctr_par(iproc+1)+1)-1
         !!else
         !!    isegend=smat%nseg
         !!end if
         !$omp parallel default(none) &
         !$omp shared(isegstart, isegend, matrixp, smat, matrix_local, isvctr, isfvctr) &
         !$omp private(iseg, ii, jorb, iiorb, jjorb)
         !$omp do
         do iseg=isegstart,isegend
             ii=smat%keyv(iseg)-1
             do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                 ii=ii+1
                 iiorb = (jorb-1)/smat%nfvctr + 1
                 jjorb = jorb - (iiorb-1)*smat%nfvctr
                 matrix_local(ii-isvctr)=matrixp(jjorb,iiorb-isfvctr)
             end do
         end do
         !$omp end do
         !$omp end parallel
     end if

     call timing(iproc,'compressd_mcpy','OF')
     call timing(iproc,'compressd_comm','ON')
     if (bigdft_mpi%nproc>1) then
         if (comm_strategy==ALLGATHERV) then
             !call mpiallred(matrix_compr(1), smat%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
             call mpi_allgatherv(matrix_local(1), nvctrp, mpi_double_precision, &
                  matrix_compr(1), nvctr_par, isvctr_par, mpi_double_precision, &
                  bigdft_mpi%mpi_comm, ierr)
             call f_free_ptr(matrix_local)
         else if (comm_strategy==GET) then
             call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, nvctrp, matrix_local, &
                  nvctr_par, isvctr_par, smat%nvctr, matrix_compr)
         else
             stop 'compress_matrix_distributed: wrong communication strategy'
         end if
         call f_free_ptr(matrix_local)
     end if

     call timing(iproc,'compressd_comm','OF')

     call f_release_routine()

  end subroutine compress_matrix_distributed


  subroutine uncompress_matrix_distributed(iproc, smat, layout, matrix_compr, matrixp)
    use module_base
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, layout
    type(sparse_matrix),intent(in) :: smat
    real(kind=8),dimension(smat%nvctr),intent(in) :: matrix_compr
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

       if (smat%nfvctrp>0) then

           call to_zero(smat%nfvctr*nfvctrp,matrixp(1,1))

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
               do jorb=smat%keyg(1,iseg),smat%keyg(2,iseg)
                   ii=ii+1
                   iiorb = (jorb-1)/smat%nfvctr + 1
                   jjorb = jorb - (iiorb-1)*smat%nfvctr
                   matrixp(jjorb,iiorb-isfvctr) = matrix_compr(ii)
               end do
           end do
           !$omp end parallel do
       end if

      call timing(iproc,'compressd_mcpy','OF')

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
     use yaml_output
     implicit none
   
     !Calling Arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=8), dimension(smat%nfvctr,smat%smmm%nfvctrp),intent(in) :: b
     real(kind=8), dimension(smat%smmm%nseq),intent(in) :: a_seq
     real(kind=8), dimension(smat%nfvctr,smat%smmm%nfvctrp), intent(out) :: c
   
     !Local variables
     !character(len=*), parameter :: subname='sparsemm'
     integer :: i,jorb,jjorb,m,mp1
     integer :: iorb, ii, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, iout
     real(kind=8) :: tt0, tt1, tt2, tt3, tt4, tt5, tt6
   
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'IR')

   
     !$omp parallel default(private) shared(smat, a_seq, b, c)
     !$omp do
     do iout=1,smat%smmm%nout
         i=smat%smmm%onedimindices(1,iout)
         iorb=smat%smmm%onedimindices(2,iout)
         ilen=smat%smmm%onedimindices(3,iout)
         ii=smat%smmm%onedimindices(4,iout)
         tt0=0.d0
         tt1=0.d0
         tt2=0.d0
         tt3=0.d0
         tt4=0.d0
         tt5=0.d0
         tt6=0.d0
   
         m=mod(ilen,7)
         if (m/=0) then
             do jorb=1,m
                jjorb=smat%smmm%ivectorindex(ii)
                tt0 = tt0 + b(jjorb,i)*a_seq(ii)
                ii=ii+1
             end do
         end if
         mp1=m+1
         do jorb=mp1,ilen,7
   
            jjorb0=smat%smmm%ivectorindex(ii+0)
            tt0 = tt0 + b(jjorb0,i)*a_seq(ii+0)
   
            jjorb1=smat%smmm%ivectorindex(ii+1)
            tt1 = tt1 + b(jjorb1,i)*a_seq(ii+1)
   
            jjorb2=smat%smmm%ivectorindex(ii+2)
            tt2 = tt2 + b(jjorb2,i)*a_seq(ii+2)
   
            jjorb3=smat%smmm%ivectorindex(ii+3)
            tt3 = tt3 + b(jjorb3,i)*a_seq(ii+3)
   
            jjorb4=smat%smmm%ivectorindex(ii+4)
            tt4 = tt4 + b(jjorb4,i)*a_seq(ii+4)
   
            jjorb5=smat%smmm%ivectorindex(ii+5)
            tt5 = tt5 + b(jjorb5,i)*a_seq(ii+5)
   
            jjorb6=smat%smmm%ivectorindex(ii+6)
            tt6 = tt6 + b(jjorb6,i)*a_seq(ii+6)
   
            ii=ii+7
         end do
         c(iorb,i) = tt0 + tt1 + tt2 + tt3 + tt4 + tt5 + tt6
     end do 
     !$omp end do
     !$omp end parallel

   
     call timing(bigdft_mpi%iproc, 'sparse_matmul ', 'RS')
       
   end subroutine sparsemm


   function orb_from_index(smat, ival)
     use sparsematrix_base, only: sparse_matrix
     implicit none
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     integer,intent(in) :: ival
     integer,dimension(2) :: orb_from_index

     orb_from_index(2) = (ival-1)/smat%nfvctr + 1
     !orb_from_index(1) = ival - (orb_from_index_fn(2)-1)*smat%nfvctr
     orb_from_index(1) = mod(ival-1,smat%nfvctr) + 1

   end function orb_from_index


end module sparsematrix

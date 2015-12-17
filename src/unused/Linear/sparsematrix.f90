
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
     integer :: jproc_send, iorb, jproc
     !integer :: window
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

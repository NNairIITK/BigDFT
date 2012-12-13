subroutine chebyshev_clean(iproc, nproc, npl, cc, tmb, ham_compr, ovrlp_compr, fermi, penalty_ev)
  use module_base
  use module_types
  use module_interfaces, except_this_one => chebyshev_clean
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl
  real(8),dimension(npl,3),intent(in) :: cc
  type(DFT_wavefunction),intent(in) :: tmb 
  real(kind=8),dimension(tmb%mad%nvctr),intent(in) :: ham_compr, ovrlp_compr
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norbp),intent(out) :: fermi
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norbp,2),intent(out) :: penalty_ev
  ! Local variables
  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb, ierr,i,j, nseq, nmaxsegk, nmaxvalk
  integer :: isegstart, isegend, iseg, ii, jjorb, is, ie, i1, i2, jproc, nout
  character(len=*),parameter :: subname='chebyshev'
  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp, t1_tmp2
  real(kind=8),dimension(:),allocatable :: ham_compr_seq, ovrlp_compr_seq, SHS, SHS_seq, vector
  real(kind=8),dimension(:,:),allocatable :: penalty_ev_seq, matrix
  real(kind=8) :: tt1, tt2, time1,time2 , tt, time_to_zero, time_vcopy, time_sparsemm, time_axpy, time_axbyz, time_copykernel
  integer,dimension(:,:,:),allocatable :: istindexarr
  integer,dimension(:),allocatable :: ivectorindex
  integer,parameter :: one=1, three=3
  integer,parameter :: number_of_matmuls=one
  integer,dimension(:,:),pointer :: onedimindices

  call timing(iproc, 'chebyshev_comp', 'ON')




  norb = tmb%orbs%norb
  norbp = tmb%orbs%norbp
  isorb = tmb%orbs%isorb

  call init_onedimindices(norbp, isorb, tmb%mad, nout, onedimindices)

  call determine_sequential_length(norbp, isorb, norb, tmb%mad, nseq, nmaxsegk, nmaxvalk)
  !write(*,'(a,6i9)') 'iproc, tmb%mad%nvctr, nseq, nmaxvalk, nmaxsegk, nout', iproc, tmb%mad%nvctr, nseq, nmaxvalk, nmaxsegk, nout

  allocate(ham_compr_seq(nseq), stat=istat)
  call memocc(istat, ham_compr_seq, 'ham_compr_seq', subname)
  allocate(ovrlp_compr_seq(nseq), stat=istat)
  call memocc(istat, ovrlp_compr_seq, 'ovrlp_compr_seq', subname)
  allocate(istindexarr(nmaxvalk,nmaxsegk,norbp), stat=istat)
  call memocc(istat, istindexarr, 'istindexarr', subname)
  allocate(ivectorindex(nseq), stat=istat)
  call memocc(istat, ivectorindex, 'ivectorindex', subname)


  if (number_of_matmuls==one) then
      allocate(SHS(tmb%mad%nvctr), stat=istat)
      call memocc(istat, SHS, 'SHS', subname)
      allocate(matrix(tmb%orbs%norb,tmb%orbs%norbp), stat=istat)
      call memocc(istat, matrix, 'matrix', subname)
      allocate(SHS_seq(nseq), stat=istat)
      call memocc(istat, SHS_seq, 'SHS_seq', subname)

      !!call uncompressMatrix(tmb%orbs%norb, tmb%mad, ovrlp_compr, matrix)
      call to_zero(norb*norbp, matrix(1,1))
      if (tmb%orbs%norbp>0) then
          isegstart=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc)+1)
          if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
              isegend=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
          else
              isegend=tmb%mad%nseg
          end if
          do iseg=isegstart,isegend
              ii=tmb%mad%keyv(iseg)-1
              do jorb=tmb%mad%keyg(1,iseg),tmb%mad%keyg(2,iseg)
                  ii=ii+1
                  iiorb = (jorb-1)/tmb%orbs%norb + 1
                  jjorb = jorb - (iiorb-1)*tmb%orbs%norb
                  matrix(jjorb,iiorb-tmb%orbs%isorb)=ovrlp_compr(ii)
              end do
          end do
      end if
  end if

  call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, ham_compr, nseq, nmaxsegk, nmaxvalk, &
       ham_compr_seq, istindexarr, ivectorindex)
  call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, ovrlp_compr, nseq, nmaxsegk, nmaxvalk, &
       ovrlp_compr_seq, istindexarr, ivectorindex)

  allocate(column(norb,norbp), stat=istat)
  call memocc(istat, column, 'column', subname)
  call to_zero(norb*norbp, column(1,1))

  if (number_of_matmuls==one) then

      call sparsemm(nseq, ham_compr_seq, matrix(1,1), column, &
           norb, norbp, ivectorindex, nout, onedimindices)
      call to_zero(norbp*norb, matrix(1,1))
      call sparsemm(nseq, ovrlp_compr_seq, column, matrix(1,1), &
           norb, norbp, ivectorindex, nout, onedimindices)
      call to_zero(tmb%mad%nvctr, SHS(1))
      
      if (tmb%orbs%norbp>0) then
          isegstart=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc)+1)
          if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
              isegend=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
          else
              isegend=tmb%mad%nseg
          end if
          do iseg=isegstart,isegend
              ii=tmb%mad%keyv(iseg)-1
              do jorb=tmb%mad%keyg(1,iseg),tmb%mad%keyg(2,iseg)
                  ii=ii+1
                  iiorb = (jorb-1)/tmb%orbs%norb + 1
                  jjorb = jorb - (iiorb-1)*tmb%orbs%norb
                  SHS(ii)=matrix(jjorb,iiorb-tmb%orbs%isorb)
                  !SHS(ii)=matrix(jjorb,iiorb)
              end do
          end do
      end if


      call mpiallred(SHS(1), tmb%mad%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, SHS, nseq, nmaxsegk, &
           nmaxvalk, SHS_seq, istindexarr, ivectorindex)

  end if
  !!write(700+iproc,*) SHS
  !!write(710+iproc,*) SHS_seq


  !!allocate(column(norb,norbp), stat=istat)
  !!call memocc(istat, column, 'column', subname)
  allocate(column_tmp(norb,norbp), stat=istat)
  call memocc(istat, column_tmp, 'column_tmp', subname)
  allocate(t(norb,norbp), stat=istat)
  call memocc(istat, t, 't', subname)
  allocate(t1(norb,norbp), stat=istat)
  call memocc(istat, t1, 't1', subname)
  allocate(t1_tmp(norb,norbp), stat=istat)
  call memocc(istat, t1_tmp, 't1_tmp', subname)
  allocate(t1_tmp2(norb,norbp), stat=istat)
  call memocc(istat, t1_tmp2, 't1_tmp2', subname)
  allocate(t2(norb,norbp), stat=istat)
  call memocc(istat, t2, 't2', subname)





time_to_zero=0.d0
time_vcopy=0.d0
time_sparsemm=0.d0
time_axpy=0.d0
time_axbyz=0.d0
time_copykernel=0.d0

tt1=mpi_wtime() 
  call to_zero(norb*norbp, column(1,1))
  do iorb=1,norbp
      iiorb=isorb+iorb
      column(iiorb,iorb)=1.d0
  end do



  call to_zero(norb*norbp, t1_tmp2(1,1))

tt2=mpi_wtime() 
time_to_zero=time_to_zero+tt2-tt1

tt1=mpi_wtime() 
  call vcopy(norb*norbp, column(1,1), 1, column_tmp(1,1), 1)
  call vcopy(norb*norbp, column(1,1), 1, t(1,1), 1)
tt2=mpi_wtime() 
time_vcopy=time_vcopy+tt2-tt1

  !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) column
tt1=mpi_wtime() 
  if (number_of_matmuls==three) then
      call sparsemm(nseq, ovrlp_compr_seq, column_tmp, column, &
           norb, norbp, ivectorindex, nout, onedimindices)
      call sparsemm(nseq, ham_compr_seq, column, column_tmp, &
           norb, norbp, ivectorindex, nout, onedimindices)
      call sparsemm(nseq, ovrlp_compr_seq, column_tmp, column, &
           norb, norbp, ivectorindex, nout, onedimindices)
  else if (number_of_matmuls==one) then
      call sparsemm(nseq, SHS_seq, column_tmp, column, &
           norb, norbp, ivectorindex, nout, onedimindices)
  end if
tt2=mpi_wtime() 
time_sparsemm=time_sparsemm+tt2-tt1


tt1=mpi_wtime() 
  call vcopy(norb*norbp, column(1,1), 1, t1(1,1), 1)
  call vcopy(norb*norbp, t1(1,1), 1, t1_tmp(1,1), 1)
tt2=mpi_wtime() 
time_vcopy=time_vcopy+tt2-tt1

  !initialize fermi
tt1=mpi_wtime() 
  call to_zero(norbp*norb, fermi(1,1))
  call to_zero(2*norb*norbp, penalty_ev(1,1,1))
tt2=mpi_wtime() 
time_to_zero=time_to_zero+tt2-tt1

tt1=mpi_wtime() 
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, 0.5d0*cc(1,1), t(1,1), fermi(:,1))
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,2))
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, cc(2,1), t1(1,1), fermi(:,1))
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, cc(2,3), t1(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, -cc(2,3), t1(1,1), penalty_ev(:,1,2))

tt2=mpi_wtime() 
time_axpy=time_axpy+tt2-tt1



  do ipl=3,npl
     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
tt1=mpi_wtime() 
     if (number_of_matmuls==three) then
         call sparsemm(nseq, ovrlp_compr_seq, t1_tmp, t1, &
              norb, norbp, ivectorindex, nout, onedimindices)
         call sparsemm(nseq, ham_compr_seq, t1, t1_tmp2, &
              norb, norbp, ivectorindex, nout, onedimindices)
         call sparsemm(nseq, ovrlp_compr_seq, t1_tmp2, t1, &
              norb, norbp, ivectorindex, nout, onedimindices)
     else if (number_of_matmuls==one) then
         call sparsemm(nseq, SHS_seq, t1_tmp, t1, &
              norb, norbp, ivectorindex, nout, onedimindices)
     end if
tt2=mpi_wtime() 
time_sparsemm=time_sparsemm+tt2-tt1
tt1=mpi_wtime() 
     call axbyz_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
tt2=mpi_wtime() 
time_axbyz=time_axbyz+tt2-tt1
tt1=mpi_wtime() 
     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, cc(ipl,1), t2(1,1), fermi(:,1))
     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, cc(ipl,3), t2(1,1), penalty_ev(:,1,1))

     if (mod(ipl,2)==1) then
         tt=cc(ipl,3)
     else
         tt=-cc(ipl,3)
     end if
     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, nout, onedimindices, tt, t2(1,1), penalty_ev(:,1,2))
tt2=mpi_wtime() 
time_axpy=time_axpy+tt2-tt1

     !update t's
tt1=mpi_wtime() 
     call copy_kernel_vectors(norbp, isorb, norb, tmb%mad, t1_tmp, t)
     call copy_kernel_vectors(norbp, isorb, norb, tmb%mad, t2, t1_tmp)
tt2=mpi_wtime() 
time_copykernel=time_copykernel+tt2-tt1
 end do




 
  call timing(iproc, 'chebyshev_comp', 'OF')
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_to_zero', time_to_zero
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_vcopy', time_vcopy
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_sparsemm', time_sparsemm
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_axpy', time_axpy
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_axbyz', time_axbyz
  !!if(iproc==0) write(*,'(a,es16.7)') 'time_copykernel', time_copykernel


  iall=-product(shape(column))*kind(column)
  deallocate(column, stat=istat)
  call memocc(istat, iall, 'column', subname)
  iall=-product(shape(column_tmp))*kind(column_tmp)
  deallocate(column_tmp, stat=istat)
  call memocc(istat, iall, 'column_tmp', subname)
  iall=-product(shape(t))*kind(t)
  deallocate(t, stat=istat)
  call memocc(istat, iall, 't', subname)
  iall=-product(shape(t1))*kind(t1)
  deallocate(t1, stat=istat)
  call memocc(istat, iall, 't1', subname)
  iall=-product(shape(t1_tmp))*kind(t1_tmp)
  deallocate(t1_tmp, stat=istat)
  call memocc(istat, iall, 't1_tmp', subname)
  iall=-product(shape(t1_tmp2))*kind(t1_tmp)
  deallocate(t1_tmp2, stat=istat)
  call memocc(istat, iall, 't1_tmp2', subname)
  iall=-product(shape(t2))*kind(t2)
  deallocate(t2, stat=istat)
  call memocc(istat, iall, 't2', subname)




  iall=-product(shape(ham_compr_seq))*kind(ham_compr_seq)
  deallocate(ham_compr_seq, stat=istat)
  call memocc(istat, iall, 'ham_compr_seq', subname)
  iall=-product(shape(ovrlp_compr_seq))*kind(ovrlp_compr_seq)
  deallocate(ovrlp_compr_seq, stat=istat)
  call memocc(istat, iall, 'ovrlp_compr_seq', subname)
  iall=-product(shape(istindexarr))*kind(istindexarr)
  deallocate(istindexarr, stat=istat)
  call memocc(istat, iall, 'istindexarr', subname)
  iall=-product(shape(ivectorindex))*kind(ivectorindex)
  deallocate(ivectorindex, stat=istat)
  call memocc(istat, iall, 'ivectorindex', subname)

  iall=-product(shape(onedimindices))*kind(onedimindices)
  deallocate(onedimindices, stat=istat)
  call memocc(istat, iall, 'onedimindices', subname)

  if (number_of_matmuls==one) then
      iall=-product(shape(SHS))*kind(SHS)
      deallocate(SHS, stat=istat)
      call memocc(istat, iall, 'SHS', subname)
      iall=-product(shape(matrix))*kind(matrix)
      deallocate(matrix, stat=istat)
      call memocc(istat, iall, 'matrix', subname)
      iall=-product(shape(SHS_seq))*kind(SHS_seq)
      deallocate(SHS_seq, stat=istat)
      call memocc(istat, iall, 'SHS_seq', subname)
  end if

end subroutine chebyshev_clean






!!subroutine chebyshev(iproc, nproc, npl, cc, tmb, ham_compr, ovrlp_compr, nvctr, orbitalindex, &
!!           sendcounts, recvcounts, senddspls, recvdspls, fermi, penalty_ev)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: iproc, nproc, npl, nvctr
!!  real(8),dimension(npl,3),intent(in) :: cc
!!  type(DFT_wavefunction),intent(in) :: tmb 
!!  real(kind=8),dimension(tmb%mad%nvctr),intent(in) :: ham_compr, ovrlp_compr
!!  integer,dimension(nvctr),intent(in) :: orbitalindex
!!  integer,dimension(0:nproc-1),intent(in) :: sendcounts, recvcounts, senddspls, recvdspls
!!  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norbp),intent(out) :: fermi
!!  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norbp,2),intent(out) :: penalty_ev
!!  ! Local variables
!!  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb, ierr,i,j, nseq, nmaxsegk, nmaxvalk
!!  integer :: isegstart, isegend, iseg, ii, jjorb, is, ie, i1, i2, jproc
!!  character(len=*),parameter :: subname='chebyshev'
!!  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp, t1_tmp2
!!  real(kind=8),dimension(:),allocatable :: ham_compr_seq, ovrlp_compr_seq, SHS, SHS_seq, vector, fermi_v, SHS_seq_v
!!  real(kind=8),dimension(:,:),allocatable :: penalty_ev_seq, matrix, penalty_ev_v
!!  real(kind=8) :: tt1, tt2, time1,time2 , tt, time_to_zero, time_vcopy, time_sparsemm, time_axpy, time_axbyz, time_copykernel
!!  integer,dimension(:,:,:),allocatable :: istindexarr
!!  integer,dimension(:),allocatable :: ivectorindex, ivectorindex_v
!!  real(kind=8),dimension(:),allocatable :: column_v, column_tmp_v, t_v, t1_v, t1_tmp_v, t1_tmp2_v, t2_v
!!  integer :: nseq_v, nmaxsegk_v, nmaxvalk_v
!!  real(kind=8),dimension(:),allocatable :: ham_compr_seq_v, ovrlp_compr_seq_v
!!  integer,dimension(:,:,:),allocatable :: istindexarr_v
!!  integer,parameter :: one=1, three=3
!!  integer,parameter :: number_of_matmuls=one
!!
!!  call timing(iproc, 'chebyshev_comp', 'ON')
!!
!!
!!
!!
!!  norb = tmb%orbs%norb
!!  norbp = tmb%orbs%norbp
!!  isorb = tmb%orbs%isorb
!!
!!  call determine_sequential_length(norbp, isorb, norb, tmb%mad, nseq, nmaxsegk, nmaxvalk)
!!  call determine_sequential_length2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, nseq_v, nmaxsegk_v, nmaxvalk_v)
!!
!!  allocate(ham_compr_seq(nseq), stat=istat)
!!  call memocc(istat, ham_compr_seq, 'ham_compr_seq', subname)
!!  allocate(ovrlp_compr_seq(nseq), stat=istat)
!!  call memocc(istat, ovrlp_compr_seq, 'ovrlp_compr_seq', subname)
!!  allocate(istindexarr(nmaxvalk,nmaxsegk,norbp), stat=istat)
!!  call memocc(istat, istindexarr, 'istindexarr', subname)
!!  allocate(ivectorindex(nseq), stat=istat)
!!  call memocc(istat, ivectorindex, 'ivectorindex', subname)
!!
!!  allocate(ham_compr_seq_v(nseq_v), stat=istat)
!!  call memocc(istat, ham_compr_seq_v, 'ham_compr_seq_v', subname)
!!  allocate(ovrlp_compr_seq_v(nseq_v), stat=istat)
!!  call memocc(istat, ovrlp_compr_seq_v, 'ovrlp_compr_seq_v', subname)
!!  allocate(istindexarr_v(nmaxvalk,nmaxsegk,norb), stat=istat)
!!  call memocc(istat, istindexarr_v, 'istindexarr_v', subname)
!!  allocate(ivectorindex_v(nseq_v), stat=istat)
!!  call memocc(istat, ivectorindex_v, 'ivectorindex_v', subname)
!!
!!  if (number_of_matmuls==one) then
!!      allocate(SHS(tmb%mad%nvctr), stat=istat)
!!      call memocc(istat, SHS, 'SHS', subname)
!!      allocate(matrix(tmb%orbs%norb,tmb%orbs%norbp), stat=istat)
!!      call memocc(istat, matrix, 'matrix', subname)
!!      allocate(SHS_seq(nseq), stat=istat)
!!      call memocc(istat, SHS_seq, 'SHS_seq', subname)
!!      allocate(SHS_seq_v(nseq_v), stat=istat)
!!      call memocc(istat, SHS_seq_v, 'SHS_seq_v', subname)
!!
!!      !!call uncompressMatrix(tmb%orbs%norb, tmb%mad, ovrlp_compr, matrix)
!!      call to_zero(norb*norbp, matrix(1,1))
!!      if (tmb%orbs%norbp>0) then
!!          isegstart=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc)+1)
!!          if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
!!              isegend=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
!!          else
!!              isegend=tmb%mad%nseg
!!          end if
!!          do iseg=isegstart,isegend
!!              ii=tmb%mad%keyv(iseg)-1
!!              do jorb=tmb%mad%keyg(1,iseg),tmb%mad%keyg(2,iseg)
!!                  ii=ii+1
!!                  iiorb = (jorb-1)/tmb%orbs%norb + 1
!!                  jjorb = jorb - (iiorb-1)*tmb%orbs%norb
!!                  matrix(jjorb,iiorb-tmb%orbs%isorb)=ovrlp_compr(ii)
!!              end do
!!          end do
!!      end if
!!  end if
!!
!!  call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, ham_compr, nseq, nmaxsegk, &
!!       nmaxvalk, ham_compr_seq, istindexarr, ivectorindex)
!!  call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, ovrlp_compr, nseq, nmaxsegk, &
!!       nmaxvalk, ovrlp_compr_seq, istindexarr, ivectorindex)
!!
!!  allocate(column(norb,norbp), stat=istat)
!!  call memocc(istat, column, 'column', subname)
!!  call to_zero(norb*norbp, column(1,1))
!!
!!  if (number_of_matmuls==one) then
!!
!!      call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, matrix(1,1), column, &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!      call to_zero(norbp*norb, matrix(1,1))
!!      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column, matrix(1,1), &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!      call to_zero(tmb%mad%nvctr, SHS(1))
!!      
!!      if (tmb%orbs%norbp>0) then
!!          isegstart=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc)+1)
!!          if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
!!              isegend=tmb%mad%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
!!          else
!!              isegend=tmb%mad%nseg
!!          end if
!!          do iseg=isegstart,isegend
!!              ii=tmb%mad%keyv(iseg)-1
!!              do jorb=tmb%mad%keyg(1,iseg),tmb%mad%keyg(2,iseg)
!!                  ii=ii+1
!!                  iiorb = (jorb-1)/tmb%orbs%norb + 1
!!                  jjorb = jorb - (iiorb-1)*tmb%orbs%norb
!!                  SHS(ii)=matrix(jjorb,iiorb-tmb%orbs%isorb)
!!                  !SHS(ii)=matrix(jjorb,iiorb)
!!              end do
!!          end do
!!      end if
!!
!!
!!      call mpiallred(SHS(1), tmb%mad%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!
!!      call enable_sequential_acces_matrix(norbp, isorb, norb, tmb%mad, SHS, nseq, nmaxsegk, &
!!           nmaxvalk, SHS_seq, istindexarr, ivectorindex)
!!      !call enable_sequential_acces_matrix2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, SHS, nseq_v, nmaxsegk_v, nmaxvalk_v, SHS_seq_v, istindexarr_v, ivectorindex_v)
!!
!!  end if
!!  !!write(700+iproc,*) SHS
!!  !!write(710+iproc,*) SHS_seq
!!
!!
!!  !!allocate(column(norb,norbp), stat=istat)
!!  !!call memocc(istat, column, 'column', subname)
!!  allocate(column_tmp(norb,norbp), stat=istat)
!!  call memocc(istat, column_tmp, 'column_tmp', subname)
!!  allocate(t(norb,norbp), stat=istat)
!!  call memocc(istat, t, 't', subname)
!!  allocate(t1(norb,norbp), stat=istat)
!!  call memocc(istat, t1, 't1', subname)
!!  allocate(t1_tmp(norb,norbp), stat=istat)
!!  call memocc(istat, t1_tmp, 't1_tmp', subname)
!!  allocate(t1_tmp2(norb,norbp), stat=istat)
!!  call memocc(istat, t1_tmp2, 't1_tmp2', subname)
!!  allocate(t2(norb,norbp), stat=istat)
!!  call memocc(istat, t2, 't2', subname)
!!
!!
!!
!!  allocate(column_v(nvctr), stat=istat)
!!  call memocc(istat, column_v, 'column_v', subname)
!!  allocate(column_tmp_v(nvctr), stat=istat)
!!  call memocc(istat, column_tmp_v, 'column_tmp_v', subname)
!!  allocate(t_v(nvctr), stat=istat)
!!  call memocc(istat, t_v, 't_v', subname)
!!  allocate(t1_v(nvctr), stat=istat)
!!  call memocc(istat, t1_v, 't1_v', subname)
!!  allocate(t1_tmp_v(nvctr), stat=istat)
!!  call memocc(istat, t1_tmp_v, 't1_tmp_v', subname)
!!  allocate(t1_tmp2_v(nvctr), stat=istat)
!!  call memocc(istat, t1_tmp2_v, 't1_tmp2_v', subname)
!!  allocate(t2_v(nvctr), stat=istat)
!!  call memocc(istat, t2_v, 't2_v', subname)
!!  allocate(fermi_v(nvctr), stat=istat)
!!  call memocc(istat, fermi_v, 'fermi_v', subname)
!!  allocate(penalty_ev_v(nvctr,2), stat=istat)
!!  call memocc(istat, penalty_ev_v, 'penalty_ev_v', subname)
!!
!!
!!
!!time_to_zero=0.d0
!!time_vcopy=0.d0
!!time_sparsemm=0.d0
!!time_axpy=0.d0
!!time_axbyz=0.d0
!!time_copykernel=0.d0
!!
!!tt1=mpi_wtime() 
!!  call to_zero(norb*norbp, column(1,1))
!!  do iorb=1,norbp
!!      iiorb=isorb+iorb
!!      column(iiorb,iorb)=1.d0
!!  end do
!!  call mpi_alltoallv(column, sendcounts, senddspls, mpi_double_precision, column_v, recvcounts, &
!!       recvdspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!  do istat=1,nvctr
!!      write(600+iproc,*) istat, column_v(istat)
!!  end do
!!
!!  !!call to_zero(nvctr, column_v(1))
!!  !!do iorb=1,norbp
!!  !!    iiorb=isorb+iorb
!!  !!    column(iiorb,iorb)=1.d0
!!  !!end do
!!
!!
!!
!!  call to_zero(norb*norbp, t1_tmp2(1,1))
!!  call to_zero(nvctr, t1_tmp2_v(1))
!!
!!tt2=mpi_wtime() 
!!time_to_zero=time_to_zero+tt2-tt1
!!
!!tt1=mpi_wtime() 
!!  call vcopy(norb*norbp, column(1,1), 1, column_tmp(1,1), 1)
!!  call vcopy(nvctr, column_v(1), 1, column_tmp_v(1), 1)
!!  call vcopy(norb*norbp, column(1,1), 1, t(1,1), 1)
!!  call vcopy(nvctr, column_v(1), 1, t_v(1), 1)
!!tt2=mpi_wtime() 
!!time_vcopy=time_vcopy+tt2-tt1
!!
!!  !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) column
!!tt1=mpi_wtime() 
!!  if (number_of_matmuls==three) then
!!      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!      call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column, column_tmp, &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!  else if (number_of_matmuls==one) then
!!      call sparsemm(nseq, SHS_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, &
!!           norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!      !call sparsemm2(nseq_v, SHS_seq_v, nmaxsegk_v, nmaxvalk_v, istindexarr_v, column_tmp_v, column_v, norb, norbp, isorb, nvctr, orbitalindex, tmb%mad, ivectorindex_v)
!!  end if
!!  write(800+iproc,*) column
!!  write(810+iproc,*) column_v
!!tt2=mpi_wtime() 
!!time_sparsemm=time_sparsemm+tt2-tt1
!!
!!
!!tt1=mpi_wtime() 
!!  call vcopy(norb*norbp, column(1,1), 1, t1(1,1), 1)
!!  call vcopy(norb*norbp, t1(1,1), 1, t1_tmp(1,1), 1)
!!  call vcopy(nvctr, column_v(1), 1, t1_v(1), 1)
!!  call vcopy(nvctr, t1_v(1), 1, t1_tmp_v(1), 1)
!!tt2=mpi_wtime() 
!!time_vcopy=time_vcopy+tt2-tt1
!!
!!  !initialize fermi
!!tt1=mpi_wtime() 
!!  call to_zero(norbp*norb, fermi(1,1))
!!  call to_zero(nvctr, fermi_v(1))
!!  !call to_zero(norbp*norb, fermider(1,isorb+1))
!!  call to_zero(2*norb*norbp, penalty_ev(1,1,1))
!!  call to_zero(2*nvctr, penalty_ev_v(1,1))
!!tt2=mpi_wtime() 
!!time_to_zero=time_to_zero+tt2-tt1
!!
!!tt1=mpi_wtime() 
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, 0.5d0*cc(1,1), t(1,1), fermi(:,1))
!!  !call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,2), t(1,1), fermider(:,isorb+1))
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,1))
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,2))
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, cc(2,1), t1(1,1), fermi(:,1))
!!  !call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,2), t1(1,1), fermider(:,isorb+1))
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, cc(2,3), t1(1,1), penalty_ev(:,1,1))
!!  call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, -cc(2,3), t1(1,1), penalty_ev(:,1,2))
!!
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, 0.5d0*cc(1,1), t_v(1), fermi_v(1))
!!  !!!call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,2), t(1,1), fermider(:,isorb+1))
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, 0.5d0*cc(1,3), t_v(1), penalty_ev_v(1,1))
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, 0.5d0*cc(1,3), t_v(1), penalty_ev_v(1,2))
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, cc(2,1), t1_v(1), fermi_v(1))
!!  !!!call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,2), t1(1,1), fermider(:,isorb+1))
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, cc(2,3), t1_v(1), penalty_ev_v(1,1))
!!  !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, -cc(2,3), t1_v(1), penalty_ev_v(1,2))
!!tt2=mpi_wtime() 
!!time_axpy=time_axpy+tt2-tt1
!!
!!  !call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)
!!
!!
!!  do ipl=3,npl
!!     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
!!tt1=mpi_wtime() 
!!     if (number_of_matmuls==three) then
!!         call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp, t1, &
!!              norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!         call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1, t1_tmp2, &
!!              norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!         call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp2, t1, &
!!              norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!     else if (number_of_matmuls==one) then
!!         call sparsemm(nseq, SHS_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp, t1, &
!!              norb, norbp, isorb, tmb%mad, ivectorindex, nout, onedimindices)
!!         !!call sparsemm2(nseq, SHS_seq_v, nmaxsegk_v, nmaxvalk_v, istindexarr_v, t1_tmp_v, t1_v, norb, norbp, isorb, nvctr, orbitalindex, tmb%mad, ivectorindex_v)
!!     end if
!!tt2=mpi_wtime() 
!!time_sparsemm=time_sparsemm+tt2-tt1
!!     !call daxbyz(norb*norbp, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
!!tt1=mpi_wtime() 
!!     call axbyz_kernel_vectors(norbp, isorb, norb, tmb%mad, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
!!     !!call axbyz_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, 2.d0, t1_v(1), -1.d0, t_v(1), t2_v(1))
!!tt2=mpi_wtime() 
!!time_axbyz=time_axbyz+tt2-tt1
!!tt1=mpi_wtime() 
!!     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, cc(ipl,1), t2(1,1), fermi(:,1))
!!     !call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,2), t2(1,1), fermider(:,isorb+1))
!!     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, cc(ipl,3), t2(1,1), penalty_ev(:,1,1))
!!
!!     !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, cc(ipl,1), t2_v(1), fermi_v(1))
!!     !call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,2), t2(1,1), fermider(:,isorb+1))
!!     !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, cc(ipl,3), t2_v(1), penalty_ev_v(1,1))
!!     if (mod(ipl,2)==1) then
!!         tt=cc(ipl,3)
!!     else
!!         tt=-cc(ipl,3)
!!     end if
!!     call axpy_kernel_vectors(norbp, isorb, norb, tmb%mad, tt, t2(1,1), penalty_ev(:,1,2))
!!     !!call axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, tt, t2_v(1), penalty_ev_v(1,2))
!!tt2=mpi_wtime() 
!!time_axpy=time_axpy+tt2-tt1
!!
!!     !update t's
!!tt1=mpi_wtime() 
!!     call copy_kernel_vectors(norbp, isorb, norb, tmb%mad, t1_tmp, t)
!!     call copy_kernel_vectors(norbp, isorb, norb, tmb%mad, t2, t1_tmp)
!!     !!call copy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, t1_tmp_v, t_v)
!!     !!call copy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, tmb%mad, t2_v, t1_tmp_v)
!!tt2=mpi_wtime() 
!!time_copykernel=time_copykernel+tt2-tt1
!! end do
!!
!!
!!
!!  call mpi_alltoallv(fermi_v, recvcounts, recvdspls, mpi_double_precision, column, sendcounts, senddspls, &
!!       mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!  do istat=1,norbp
!!      do iall=1,norb
!!          write(500+iproc,'(2i8,es20.10)') istat, iall, column(iall,istat)
!!          write(520+iproc,'(2i8,es20.10)') istat, iall, fermi(iall,istat)
!!      end do
!!  end do
!!
!!
!! 
!!  call timing(iproc, 'chebyshev_comp', 'OF')
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_to_zero', time_to_zero
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_vcopy', time_vcopy
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_sparsemm', time_sparsemm
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_axpy', time_axpy
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_axbyz', time_axbyz
!!  if(iproc==0) write(*,'(a,es16.7)') 'time_copykernel', time_copykernel
!!
!!
!!  iall=-product(shape(column))*kind(column)
!!  deallocate(column, stat=istat)
!!  call memocc(istat, iall, 'column', subname)
!!  iall=-product(shape(column_tmp))*kind(column_tmp)
!!  deallocate(column_tmp, stat=istat)
!!  call memocc(istat, iall, 'column_tmp', subname)
!!  iall=-product(shape(t))*kind(t)
!!  deallocate(t, stat=istat)
!!  call memocc(istat, iall, 't', subname)
!!  iall=-product(shape(t1))*kind(t1)
!!  deallocate(t1, stat=istat)
!!  call memocc(istat, iall, 't1', subname)
!!  iall=-product(shape(t1_tmp))*kind(t1_tmp)
!!  deallocate(t1_tmp, stat=istat)
!!  call memocc(istat, iall, 't1_tmp', subname)
!!  iall=-product(shape(t1_tmp2))*kind(t1_tmp)
!!  deallocate(t1_tmp2, stat=istat)
!!  call memocc(istat, iall, 't1_tmp2', subname)
!!  iall=-product(shape(t2))*kind(t2)
!!  deallocate(t2, stat=istat)
!!  call memocc(istat, iall, 't2', subname)
!!
!!
!!
!!  iall=-product(shape(column_v))*kind(column_v)
!!  deallocate(column_v, stat=istat)
!!  call memocc(istat, iall, 'column_v', subname)
!!  iall=-product(shape(column_tmp_v))*kind(column_tmp_v)
!!  deallocate(column_tmp_v, stat=istat)
!!  call memocc(istat, iall, 'column_tmp_v', subname)
!!  iall=-product(shape(t_v))*kind(t_v)
!!  deallocate(t_v, stat=istat)
!!  call memocc(istat, iall, 't_v', subname)
!!  iall=-product(shape(t1_v))*kind(t1_v)
!!  deallocate(t1_v, stat=istat)
!!  call memocc(istat, iall, 't1_v', subname)
!!  iall=-product(shape(t1_tmp_v))*kind(t1_tmp_v)
!!  deallocate(t1_tmp_v, stat=istat)
!!  call memocc(istat, iall, 't1_tmp_v', subname)
!!  iall=-product(shape(t1_tmp2_v))*kind(t1_tmp_v)
!!  deallocate(t1_tmp2_v, stat=istat)
!!  call memocc(istat, iall, 't1_tmp2_v', subname)
!!  iall=-product(shape(t2_v))*kind(t2_v)
!!  deallocate(t2_v, stat=istat)
!!  call memocc(istat, iall, 't2_v', subname)
!!  iall=-product(shape(fermi_v))*kind(fermi_v)
!!  deallocate(fermi_v, stat=istat)
!!  call memocc(istat, iall, 'fermi_v', subname)
!!  iall=-product(shape(penalty_ev_v))*kind(penalty_ev_v)
!!  deallocate(penalty_ev_v, stat=istat)
!!  call memocc(istat, iall, 'penalty_ev_v', subname)
!!
!!
!!
!!
!!  iall=-product(shape(ham_compr_seq))*kind(ham_compr_seq)
!!  deallocate(ham_compr_seq, stat=istat)
!!  call memocc(istat, iall, 'ham_compr_seq', subname)
!!  iall=-product(shape(ovrlp_compr_seq))*kind(ovrlp_compr_seq)
!!  deallocate(ovrlp_compr_seq, stat=istat)
!!  call memocc(istat, iall, 'ovrlp_compr_seq', subname)
!!  iall=-product(shape(istindexarr))*kind(istindexarr)
!!  deallocate(istindexarr, stat=istat)
!!  call memocc(istat, iall, 'istindexarr', subname)
!!  iall=-product(shape(ivectorindex))*kind(ivectorindex)
!!  deallocate(ivectorindex, stat=istat)
!!  call memocc(istat, iall, 'ivectorindex', subname)
!!
!!
!!
!!  iall=-product(shape(ham_compr_seq_v))*kind(ham_compr_seq_v)
!!  deallocate(ham_compr_seq_v, stat=istat)
!!  call memocc(istat, iall, 'ham_compr_seq_v', subname)
!!  iall=-product(shape(ovrlp_compr_seq_v))*kind(ovrlp_compr_seq_v)
!!  deallocate(ovrlp_compr_seq_v, stat=istat)
!!  call memocc(istat, iall, 'ovrlp_compr_seq_v', subname)
!!  iall=-product(shape(istindexarr_v))*kind(istindexarr_v)
!!  deallocate(istindexarr_v, stat=istat)
!!  call memocc(istat, iall, 'istindexarr_v', subname)
!!  iall=-product(shape(ivectorindex_v))*kind(ivectorindex_v)
!!  deallocate(ivectorindex_v, stat=istat)
!!  call memocc(istat, iall, 'ivectorindex_v', subname)
!!
!!
!!
!!  if (number_of_matmuls==one) then
!!      iall=-product(shape(SHS))*kind(SHS)
!!      deallocate(SHS, stat=istat)
!!      call memocc(istat, iall, 'SHS', subname)
!!      iall=-product(shape(matrix))*kind(matrix)
!!      deallocate(matrix, stat=istat)
!!      call memocc(istat, iall, 'matrix', subname)
!!      iall=-product(shape(SHS_seq))*kind(SHS_seq)
!!      deallocate(SHS_seq, stat=istat)
!!      call memocc(istat, iall, 'SHS_seq', subname)
!!      iall=-product(shape(SHS_seq_v))*kind(SHS_seq_v)
!!      deallocate(SHS_seq_v, stat=istat)
!!      call memocc(istat, iall, 'SHS_seq_v', subname)
!!  end if
!!
!!end subroutine chebyshev




! Performs z = a*x + b*y
subroutine daxbyz(n, a, x, b, y, z)
  implicit none

  ! Calling arguments
  integer,intent(in) :: n
  real(8),intent(in) :: a, b
  real(kind=8),dimension(n),intent(in) :: x, y
  real(kind=8),dimension(n),intent(out) :: z

  ! Local variables
  integer :: i, m, mp1

  m=mod(n,7)

  if (m/=7) then
      do i=1,m
          z(i) = a*x(i) + b*y(i)
      end do
  end if
  
  mp1=m+1
  do i=mp1,n,7
      z(i+0) = a*x(i+0) + b*y(i+0)
      z(i+1) = a*x(i+1) + b*y(i+1)
      z(i+2) = a*x(i+2) + b*y(i+2)
      z(i+3) = a*x(i+3) + b*y(i+3)
      z(i+4) = a*x(i+4) + b*y(i+4)
      z(i+5) = a*x(i+5) + b*y(i+5)
      z(i+6) = a*x(i+6) + b*y(i+6)
  end do


end subroutine daxbyz




! Performs z = a*x + b*y
subroutine axbyz_kernel_vectors(norbp, isorb, norb, mad, nout, onedimindices, a, x, b, y, z)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nout
  type(matrixDescriptors),intent(in) :: mad
  integer,dimension(4,nout),intent(in) :: onedimindices
  real(8),intent(in) :: a, b
  real(kind=8),dimension(norb,norbp),intent(in) :: x, y
  real(kind=8),dimension(norb,norbp),intent(out) :: z

  ! Local variables
  integer :: i, m, mp1, jorb, iseg, ii, iorb

  !$omp parallel default(private) shared(nout, onedimindices,a, b, x, y, z)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      z(jorb,iorb)=a*x(jorb,iorb)+b*y(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel

  !!do i=1,norbp
  !!    ii=isorb+i
  !!    do iseg=1,mad%kernel_nseg(ii)
  !!        m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
  !!        if (m/=0) then
  !!            do jorb=mad%kernel_segkeyg(1,iseg,ii),mad%kernel_segkeyg(1,iseg,ii)+m-1
  !!                z(jorb,i)=a*x(jorb,i)+b*y(jorb,i)
  !!            end do
  !!        end if
  !!        do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
  !!            z(jorb+0,i)=a*x(jorb+0,i)+b*y(jorb+0,i)
  !!            z(jorb+1,i)=a*x(jorb+1,i)+b*y(jorb+1,i)
  !!            z(jorb+2,i)=a*x(jorb+2,i)+b*y(jorb+2,i)
  !!            z(jorb+3,i)=a*x(jorb+3,i)+b*y(jorb+3,i)
  !!            z(jorb+4,i)=a*x(jorb+4,i)+b*y(jorb+4,i)
  !!            z(jorb+5,i)=a*x(jorb+5,i)+b*y(jorb+5,i)
  !!            z(jorb+6,i)=a*x(jorb+6,i)+b*y(jorb+6,i)
  !!        end do
  !!    end do
  !!end do

end subroutine axbyz_kernel_vectors






subroutine sparsemm(nseq, a_seq, b, c, norb, norbp, ivectorindex, nout, onedimindices)
  use module_base
  use module_types

  implicit none

  !Calling Arguments
  integer, intent(in) :: norb,norbp,nseq
  real(kind=8), dimension(norb,norbp),intent(in) :: b
  real(kind=8), dimension(nseq),intent(in) :: a_seq
  real(kind=8), dimension(norb,norbp), intent(out) :: c
  integer,dimension(nseq),intent(in) :: ivectorindex
  integer,intent(in) :: nout
  integer,dimension(4,nout) :: onedimindices

  !Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,norbthrd,orb_rest,tid,istart,iend, mp1, iii
  integer :: iorb, jseg, ii, ii0, ii2, is, ie, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, iout
  character(len=*),parameter :: subname='sparsemm'
  real(8) :: ncount, t1, t2, ddot, tt, tt2, t3, ncount2
  integer :: ierr, imin, imax, nthreads



      !$omp parallel default(private) shared(ivectorindex, a_seq, b, c, onedimindices, nout)
      !$omp do
      !!dir$ parallel
      do iout=1,nout
          i=onedimindices(1,iout)
          iorb=onedimindices(2,iout)
          ilen=onedimindices(3,iout)
          ii0=onedimindices(4,iout)
          ii2=0
          tt=0.d0

          m=mod(ilen,7)
          if (m/=0) then
              do jorb=1,m
                 jjorb=ivectorindex(ii0+ii2)
                 tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
                 ii2=ii2+1
              end do
          end if
          mp1=m+1
          do jorb=mp1,ilen,7

             jjorb0=ivectorindex(ii0+ii2+0)
             tt = tt + b(jjorb0,i)*a_seq(ii0+ii2+0)

             jjorb1=ivectorindex(ii0+ii2+1)
             tt = tt + b(jjorb1,i)*a_seq(ii0+ii2+1)

             jjorb2=ivectorindex(ii0+ii2+2)
             tt = tt + b(jjorb2,i)*a_seq(ii0+ii2+2)

             jjorb3=ivectorindex(ii0+ii2+3)
             tt = tt + b(jjorb3,i)*a_seq(ii0+ii2+3)

             jjorb4=ivectorindex(ii0+ii2+4)
             tt = tt + b(jjorb4,i)*a_seq(ii0+ii2+4)

             jjorb5=ivectorindex(ii0+ii2+5)
             tt = tt + b(jjorb5,i)*a_seq(ii0+ii2+5)

             jjorb6=ivectorindex(ii0+ii2+6)
             tt = tt + b(jjorb6,i)*a_seq(ii0+ii2+6)

             ii2=ii2+7
          end do
          c(iorb,i)=tt
          !end do
      end do 
      !$omp end do
      !$omp end parallel

    
end subroutine sparsemm



subroutine copy_kernel_vectors(norbp, isorb, norb, mad, a, b)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(norb,norbp),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(out) :: b

  ! Local variables
  integer :: i, m, jorb, mp1, iseg, ii


  do i=1,norbp
      ii=isorb+i
      do iseg=1,mad%kernel_nseg(ii)
          m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
          if (m/=0) then
              do jorb=mad%kernel_segkeyg(1,iseg,ii),mad%kernel_segkeyg(1,iseg,ii)+m-1
                  b(jorb,i)=a(jorb,i)
              end do
          end if
          do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
              b(jorb+0,i)=a(jorb+0,i)
              b(jorb+1,i)=a(jorb+1,i)
              b(jorb+2,i)=a(jorb+2,i)
              b(jorb+3,i)=a(jorb+3,i)
              b(jorb+4,i)=a(jorb+4,i)
              b(jorb+5,i)=a(jorb+5,i)
              b(jorb+6,i)=a(jorb+6,i)
          end do
      end do
  end do


end subroutine copy_kernel_vectors




subroutine axpy_kernel_vectors(norbp, isorb, norb, mad, nout, onedimindices, a, x, y)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nout
  type(matrixDescriptors),intent(in) :: mad
  integer,dimension(4,nout),intent(in) :: onedimindices
  real(kind=8),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(in) :: x
  real(kind=8),dimension(norb,norbp),intent(out) :: y

  ! Local variables
  integer :: i, m, jorb, mp1, iseg, ii, iorb

  !$omp parallel default(private) shared(nout, onedimindices, y, x, a)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      y(jorb,iorb)=y(jorb,iorb)+a*x(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel

  !!do i=1,norbp
  !!    ii=isorb+i
  !!    do iseg=1,mad%kernel_nseg(ii)
  !!        m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
  !!        if (m/=0) then
  !!            do jorb=mad%kernel_segkeyg(1,iseg,ii),mad%kernel_segkeyg(1,iseg,ii)+m-1
  !!                y(jorb,i)=y(jorb,i)+a*x(jorb,i)
  !!            end do
  !!        end if
  !!        do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
  !!            y(jorb+0,i)=y(jorb+0,i)+a*x(jorb+0,i)
  !!            y(jorb+1,i)=y(jorb+1,i)+a*x(jorb+1,i)
  !!            y(jorb+2,i)=y(jorb+2,i)+a*x(jorb+2,i)
  !!            y(jorb+3,i)=y(jorb+3,i)+a*x(jorb+3,i)
  !!            y(jorb+4,i)=y(jorb+4,i)+a*x(jorb+4,i)
  !!            y(jorb+5,i)=y(jorb+5,i)+a*x(jorb+5,i)
  !!            y(jorb+6,i)=y(jorb+6,i)+a*x(jorb+6,i)
  !!        end do
  !!    end do
  !!end do

end subroutine axpy_kernel_vectors




subroutine enable_sequential_acces_matrix(norbp, isorb, norb, mad, a, nseq, nmaxsegk, nmaxvalk, a_seq, &
           istindexarr, ivectorindex)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nseq, nmaxsegk, nmaxvalk
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(mad%nvctr),intent(in) :: a
  real(kind=8),dimension(nseq),intent(out) :: a_seq
  integer,dimension(nmaxvalk,nmaxsegk,norbp),intent(out) :: istindexarr
  integer,dimension(nseq),intent(out) :: ivectorindex

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb, iii


  ii=1
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,mad%kernel_nseg(iii)
          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
              istindexarr(iorb-mad%kernel_segkeyg(1,iseg,iii)+1,iseg,i)=ii
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  jj=1
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      jjorb = jorb - (iorb-1)*norb
                      a_seq(ii)=a(mad%keyv(jseg)+jj-1)
                      ivectorindex(ii)=jjorb
                      jj = jj+1
                      ii = ii+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine enable_sequential_acces_matrix

!!subroutine enable_sequential_acces_vector(norbp, norb, isorb, mad, b, nseq, bseq, indexarr)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: norbp, norb, isorb
!!  type(matrixDescriptors),intent(in) :: mad
!!  real(kind=8),dimension(norb,norbp),intent(in) :: b
!!  integer,intent(out) :: nseq
!!  real(kind=8),dimension(:),pointer,intent(out) :: bseq
!!  integer,dimension(norb,norbp),intent(out) :: indexarr
!!
!!  ! Local variables
!!  integer :: i, iii, iseg, iorb, iseq, istat
!!  character(len=*),parameter :: subname='enable_sequential_acces_vector'
!!
!!  nseq=0
!!  do i = 1,norbp
!!     iii=isorb+i
!!     do iseg=1,mad%kernel_nseg(iii)
!!          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
!!              nseq=nseq+1
!!              indexarr(iorb,i)=iseq
!!          end do
!!     end do
!!  end do 
!!
!!  allocate(bseq(nseq), stat=istat)
!!  call memocc(istat, bseq, 'bseq', subname)
!!
!!  iseq=0
!!  do i = 1,norbp
!!     iii=isorb+i
!!     do iseg=1,mad%kernel_nseg(iii)
!!          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
!!              iseq=iseq+1
!!              bseq(iseq)=b(iorb,i)
!!          end do
!!     end do
!!  end do 
!!
!!end subroutine enable_sequential_acces_vector
!!
!!
!!subroutine compress_vector
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer :: norbp, norb, isorb
!!  type(matrixDescriptors),intent(in) :: mad
!!  real(kind=8),dimension(
!!
!!
!!  iseq=0
!!  do i = 1,norbp
!!     iii=isorb+i
!!     do iseg=1,mad%kernel_nseg(iii)
!!          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
!!              iseq=iseq+1
!!              bseq(iseq)=b(iorb,i)
!!          end do
!!     end do
!!  end do 
!!
!!end subroutine uncompress_vector

!!subroutine enable_sequential_acces_vector(norbp, norb, mad, b, nseq, b_seq)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: norbp, norb, nseq
!!  type(matrixDescriptors),intent(in) :: mad
!!  real(kind=8),dimension(norb,norbp),intent(in) :: b
!!  real(kind=8),dimension(nseq),intent(out) :: b_seq
!!
!!  ! Local variables
!!  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
!!  integer :: iorb, jseg, ii, ist_iorb
!!
!!  stop 'DOES NOT WORK'
!!
!!  ii=1
!!  do i = 1,norbp
!!     do iseg=1,mad%kernel_nseg(i)
!!          do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
!!              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
!!                  jj=0
!!                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
!!                      jjorb = jorb - (iorb-1)*norb
!!                      b_seq(ii)=b(jjorb,i)
!!                      jj = jj+1
!!                      ii = ii+1
!!                  end do
!!              end do
!!          end do
!!     end do
!!  end do 
!!
!!end subroutine enable_sequential_acces_vector



subroutine determine_sequential_length(norbp, isorb, norb, mad, nseq, nmaxsegk, nmaxvalk)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb
  type(matrixDescriptors),intent(in) :: mad
  integer,intent(out) :: nseq, nmaxsegk, nmaxvalk

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb

  nseq=0
  nmaxsegk=0
  nmaxvalk=0
  do i = 1,norbp
     ii=isorb+i
     nmaxsegk=max(nmaxsegk,mad%kernel_nseg(ii))
     do iseg=1,mad%kernel_nseg(ii)
          nmaxvalk=max(nmaxvalk,mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1)
          do iorb=mad%kernel_segkeyg(1,iseg,ii),mad%kernel_segkeyg(2,iseg,ii)
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      nseq=nseq+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine determine_sequential_length




subroutine determine_load_balancing(iproc, nproc, orbs, mad, nvctr, orbitalindex, sendcounts, recvcounts, senddspls, recvdspls)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(matrixDescriptors),intent(in) :: mad
  integer,intent(out) :: nvctr
  integer,dimension(:),pointer,intent(out) :: orbitalindex 
  integer,dimension(0:nproc),intent(out) :: sendcounts, recvcounts, senddspls, recvdspls

  ! Local variables
  integer :: i, iseg, iorb, jseg, ierr, jproc, jjproc, istat, iall, ii, iiorb, jorb, jj, is, ie, i1, i2
  real(kind=8) :: op, op_ideal, opf
  real(kind=8),dimension(:),allocatable :: op_arr
  integer,dimension(:,:),allocatable :: istartend
  character(len=*),parameter :: subname='determine_load_balancing'

  allocate(op_arr(0:nproc-1), stat=istat)
  call memocc(istat, op_arr, 'op_arr', subname)
  allocate(istartend(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend, 'istartend', subname)

  call to_zero(nproc, op_arr(0))

  op=0.d0
  do i=1,orbs%norbp
      iiorb=orbs%isorb+i
      do iseg=1,mad%kernel_nseg(iiorb)
          do iorb=mad%kernel_segkeyg(1,iseg,iiorb),mad%kernel_segkeyg(2,iseg,iiorb)
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  op = op + dble(mad%keyg(2,jseg)-mad%keyg(1,jseg)+1)
              end do
          end do
      end do
  end do
  op_arr(iproc)=op

  call mpiallred(op, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call mpiallred(op_arr(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  ! Ideal load per process
  op_ideal=op/dble(nproc)

  ! "first operation" for each process
  opf=dble(iproc)*op_ideal+1.d0


  ii=0
  jjproc=0
  istartend(1,0)=1
  op=0.d0
  outer_loop: do jproc=0,nproc-1
      do i=1,orbs%norb_par(jproc,0)
          ii=ii+1
          do iseg=1,mad%kernel_nseg(ii)
              do iorb=mad%kernel_segkeyg(1,iseg,ii),mad%kernel_segkeyg(2,iseg,ii)
                  do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                      op = op + dble(mad%keyg(2,jseg)-mad%keyg(1,jseg)+1)
                  end do
                  if (op>=op_ideal) then
                      istartend(2,jjproc)=(ii-1)*orbs%norb+iorb
                      jjproc=jjproc+1
                      if (jjproc==nproc) exit outer_loop
                      istartend(1,jjproc)=istartend(2,jjproc-1)+1
                      op=0.d0
                  end if
              end do
          end do
      end do
  end do outer_loop
  istartend(2,min(jjproc,nproc-1))=orbs%norb**2

  ! Fill remainig processes, if there are any
  do jproc=jjproc+1,nproc-1
      istartend(1,jproc)=0
      istartend(2,jproc)=-1
  end do

  ! Some check
  ii=0
  do jproc=0,nproc-1
      ii = ii + istartend(2,jproc)-istartend(1,jproc)+1
  end do
  if (ii/=orbs%norb**2) stop 'ii/=orbs%norb**2'
  
  nvctr=istartend(2,iproc)-istartend(1,iproc)+1
  allocate(orbitalindex(nvctr), stat=istat)
  call memocc(istat, orbitalindex, 'orbitalindex', subname)

  ii=0
  jj=0
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ii=ii+1
          if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
              jj=jj+1
              orbitalindex(jj)=(iorb-1)*orbs%norb+jorb
          end if    
      end do
  end do



  do jproc=0,nproc-1
      is=orbs%isorb*orbs%norb+1
      ie=is+orbs%norbp*orbs%norb-1
      i1=max(is,istartend(1,jproc))
      i2=min(ie,istartend(2,jproc))
      sendcounts(jproc)=max(i2-i1+1,0)
  end do

  do jproc=0,nproc-1
      is=orbs%isorb_par(jproc)*orbs%norb+1
      ie=is+orbs%norb_par(jproc,0)*orbs%norb-1
      i1=max(is,istartend(1,iproc))
      i2=min(ie,istartend(2,iproc))
      recvcounts(jproc)=max(i2-i1+1,0)
  end do

  senddspls(0)=0
  recvdspls(0)=0
  do jproc=1,nproc-1
      senddspls(jproc)=senddspls(jproc-1)+sendcounts(jproc-1)
      recvdspls(jproc)=recvdspls(jproc-1)+recvcounts(jproc-1)
  end do




  iall=-product(shape(op_arr))*kind(op_arr)
  deallocate(op_arr, stat=istat)
  call memocc(istat, iall, 'op_arr', subname)
  iall=-product(shape(istartend))*kind(istartend)
  deallocate(istartend, stat=istat)
  call memocc(istat, iall, 'istartend', subname)

end subroutine determine_load_balancing





subroutine determine_sequential_length2(norbp, isorb, norb, nvctr, orbitalindex, mad, nseq, nmaxsegk, nmaxvalk)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  type(matrixDescriptors),intent(in) :: mad
  integer,intent(out) :: nseq, nmaxsegk, nmaxvalk

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb, ipt, iminorb, imaxorb, imin, imax, iproc, ierr, is, ie

  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

  nseq=0
  nmaxsegk=0
  nmaxvalk=0
  do i=iminorb,imaxorb
     ii=i
     nmaxsegk=max(nmaxsegk,mad%kernel_nseg(ii))
     do iseg=1,mad%kernel_nseg(ii)
          !!write(*,'(a,6i9)') 'iproc, iseg, imin, mad%kernel_segkeyg(1,iseg,ii), imax, mad%kernel_segkeyg(2,iseg,ii)', &
          !!            iproc, iseg, imin,mad%kernel_segkeyg(1,iseg,ii), imax, mad%kernel_segkeyg(2,iseg,ii) 
          is = (ii-1)*norb+mad%kernel_segkeyg(1,iseg,ii)
          ie = (ii-1)*norb+mad%kernel_segkeyg(2,iseg,ii)
          nmaxvalk=max(nmaxvalk,ie-is+1)
          do iiorb=max(imin,is),min(imax,ie)
              iorb=mod(iiorb-1,norb)+1
              !write(*,'(a,4i9)') 'iproc, iorb, mad%istsegline(iorb), mad%istsegline(iorb)+mad%nsegline(iorb)-1',&
              !            iproc, iorb, mad%istsegline(iorb), mad%istsegline(iorb)+mad%nsegline(iorb)-1
              !!write(*,'(a,5i9)') 'iproc, iiorb, iorb, mad%istsegline(iorb), mad%istsegline(iorb)+mad%nsegline(iorb)-1',&
              !!            iproc, iiorb, iorb, mad%istsegline(iorb), mad%istsegline(iorb)+mad%nsegline(iorb)-1
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      nseq=nseq+1
                  end do
              end do
          end do
     end do
  end do 

  !!write(*,'(a,4i8)') 'iproc, iminorb, imaxorb, nseq', iproc, iminorb, imaxorb, nseq

end subroutine determine_sequential_length2






subroutine enable_sequential_acces_matrix2(norbp, isorb, norb, nvctr, orbitalindex, mad, a, nseq, nmaxsegk, nmaxvalk, a_seq, &
           istindexarr, ivectorindex)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nseq, nmaxsegk, nmaxvalk, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(mad%nvctr),intent(in) :: a
  real(kind=8),dimension(nseq),intent(out) :: a_seq
  integer,dimension(nmaxvalk,nmaxsegk,norb),intent(out) :: istindexarr
  integer,dimension(nseq),intent(out) :: ivectorindex

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb, iii, ipt, imin, imax, iminorb, imaxorb, iproc, ierr, is, ie

  call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  ii=1
  ipt=0
  do i=iminorb,imaxorb
     iii=i
     do iseg=1,mad%kernel_nseg(iii)
          is = (iii-1)*norb+mad%kernel_segkeyg(1,iseg,iii)
          ie = (iii-1)*norb+mad%kernel_segkeyg(2,iseg,iii)
          do iiorb=max(imin,is),min(imax,ie)
              iorb=mod(iiorb-1,norb)+1
              istindexarr(iorb-mad%kernel_segkeyg(1,iseg,iii)+1,iseg,i)=ii
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  jj=1
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      jjorb = jorb - (iorb-1)*norb
                      a_seq(ii)=a(mad%keyv(jseg)+jj-1)
                      !ivectorindex(ii)=jjorb+(iminorb-1)*norb-imin+1 + (i-iminorb)*norb
                      ivectorindex(ii)=jjorb+(i-1)*norb-imin+1
                      write(*,'(a,6i9)') 'iproc, ii, jorb, jjorb, imin, ivectorindex(ii)', &
                                  iproc, ii, jorb, jjorb, imin, ivectorindex(ii)
                      jj = jj+1
                      ii = ii+1
                  end do
              end do
          end do
     end do
  end do 

  !write(*,'(a,5i9)') 'iproc, ii, nseq, minval(ivectorindex), maxval(ivectorindex)', iproc, ii, nseq, minval(ivectorindex), maxval(ivectorindex)


end subroutine enable_sequential_acces_matrix2


subroutine axbyz_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, mad, a, x, b, y, z)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  type(matrixDescriptors),intent(in) :: mad
  real(8),intent(in) :: a, b
  real(kind=8),dimension(nvctr),intent(in) :: x, y
  real(kind=8),dimension(nvctr),intent(out) :: z

  ! Local variables
  integer :: i, m, mp1, jorb, iseg, ii, imin, imax, iminorb, imaxorb, is, ie, ipt


  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  ipt=0
  do i=iminorb,imaxorb
      ii=i
      do iseg=1,mad%kernel_nseg(ii)
          !m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
          !if (m/=0) then
              is = (ii-1)*norb+mad%kernel_segkeyg(1,iseg,ii)
              ie = (ii-1)*norb+mad%kernel_segkeyg(2,iseg,ii)
              do jorb=max(imin,is),min(imax,ie)
                  ipt=ipt+1
                  z(ipt)=a*x(ipt)+b*y(ipt)
              end do
          !end if
          !do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
          !    z(jorb+0,i)=a*x(jorb+0,i)+b*y(jorb+0,i)
          !    z(jorb+1,i)=a*x(jorb+1,i)+b*y(jorb+1,i)
          !    z(jorb+2,i)=a*x(jorb+2,i)+b*y(jorb+2,i)
          !    z(jorb+3,i)=a*x(jorb+3,i)+b*y(jorb+3,i)
          !    z(jorb+4,i)=a*x(jorb+4,i)+b*y(jorb+4,i)
          !    z(jorb+5,i)=a*x(jorb+5,i)+b*y(jorb+5,i)
          !    z(jorb+6,i)=a*x(jorb+6,i)+b*y(jorb+6,i)
          !end do
      end do
  end do

end subroutine axbyz_kernel_vectors2


subroutine axpy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, mad, a, x, y)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),intent(in) :: a
  real(kind=8),dimension(nvctr),intent(in) :: x
  real(kind=8),dimension(nvctr),intent(out) :: y

  ! Local variables
  integer :: i, m, jorb, mp1, iseg, ii, ipt, imin, imax, iminorb, imaxorb, is, ie

  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  ipt=0
  do i=iminorb,imaxorb
      ii=i
      do iseg=1,mad%kernel_nseg(ii)
          !m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
          !if (m/=0) then
              is = (ii-1)*norb+mad%kernel_segkeyg(1,iseg,ii)
              ie = (ii-1)*norb+mad%kernel_segkeyg(2,iseg,ii)
              do jorb=max(imin,is),min(imax,ie)
                  ipt=ipt+1
                  y(ipt)=y(ipt)+a*x(ipt)
              end do
          !end if
          !do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
          !    y(jorb+0,i)=y(jorb+0,i)+a*x(jorb+0,i)
          !    y(jorb+1,i)=y(jorb+1,i)+a*x(jorb+1,i)
          !    y(jorb+2,i)=y(jorb+2,i)+a*x(jorb+2,i)
          !    y(jorb+3,i)=y(jorb+3,i)+a*x(jorb+3,i)
          !    y(jorb+4,i)=y(jorb+4,i)+a*x(jorb+4,i)
          !    y(jorb+5,i)=y(jorb+5,i)+a*x(jorb+5,i)
          !    y(jorb+6,i)=y(jorb+6,i)+a*x(jorb+6,i)
          !end do
      end do
  end do

end subroutine axpy_kernel_vectors2





subroutine copy_kernel_vectors2(norbp, isorb, norb, nvctr, orbitalindex, mad, a, b)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(nvctr),intent(in) :: a
  real(kind=8),dimension(nvctr),intent(out) :: b

  ! Local variables
  integer :: i, m, jorb, mp1, iseg, ii, ipt, imin, imax, iminorb, imaxorb, is, ie

  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  ipt=0
  do i=iminorb,imaxorb
      ii=i
      do iseg=1,mad%kernel_nseg(ii)
          !m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
          !if (m/=0) then
              is = (ii-1)*norb+mad%kernel_segkeyg(1,iseg,ii)
              ie = (ii-1)*norb+mad%kernel_segkeyg(2,iseg,ii)
              do jorb=max(imin,is),min(imax,ie)
                  ipt=ipt+1
                  b(ipt)=a(ipt)
              end do
          !end if
          !do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
          !    b(jorb+0,i)=a(jorb+0,i)
          !    b(jorb+1,i)=a(jorb+1,i)
          !    b(jorb+2,i)=a(jorb+2,i)
          !    b(jorb+3,i)=a(jorb+3,i)
          !    b(jorb+4,i)=a(jorb+4,i)
          !    b(jorb+5,i)=a(jorb+5,i)
          !    b(jorb+6,i)=a(jorb+6,i)
          !end do
      end do
  end do


end subroutine copy_kernel_vectors2






subroutine sparsemm2(nseq, a_seq, nmaxsegk, nmaxvalk, istindexarr, b, c, norb, norbp, isorb, nvctr, orbitalindex, mad, ivectorindex)
  use module_base
  use module_types

  implicit none

  !Calling Arguments
  type(matrixDescriptors),intent(in) :: mad
  integer, intent(in) :: norb,norbp,isorb,nseq,nmaxsegk,nmaxvalk, nvctr
  integer,dimension(nvctr),intent(in) :: orbitalindex
  !real(kind=8), dimension(norb,norbp),intent(in) :: b
  real(kind=8), dimension(nvctr),intent(in) :: b
  real(kind=8), dimension(nseq),intent(in) :: a_seq
  integer,dimension(nmaxvalk,nmaxsegk,norb),intent(in) :: istindexarr
  !real(kind=8), dimension(norb,norbp), intent(out) :: c
  real(kind=8), dimension(nvctr), intent(out) :: c
  integer,dimension(nseq),intent(in) :: ivectorindex

  !Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,norbthrd,orb_rest,tid,istart,iend, mp1, iii
  integer :: iorb, jseg, ii, ii0, ii2, is, ie, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, ipt
  integer,dimension(:), allocatable :: n
  character(len=*),parameter :: subname='sparsemm'
  real(8) :: ncount, t1, t2, ddot, tt, tt2, t3, ncount2
  integer :: iproc, ierr, nthreads, imin, imax, iminorb, imaxorb, korb

  call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

  imin=minval(orbitalindex)
  imax=maxval(orbitalindex)

  iminorb=(imin-1)/norb+1
  imaxorb=(imax-1)/norb+1

  ipt=0
  do i=iminorb,imaxorb
      ii=i
      do iseg=1,mad%kernel_nseg(ii)
          !m=mod(mad%kernel_segkeyg(2,iseg,ii)-mad%kernel_segkeyg(1,iseg,ii)+1,7)
          !if (m/=0) then
              is = (ii-1)*norb+mad%kernel_segkeyg(1,iseg,ii)
              ie = (ii-1)*norb+mad%kernel_segkeyg(2,iseg,ii)
              do jorb=max(imin,is),min(imax,ie)
                  ipt=ipt+1
                  c(ipt)=0.d0
              end do
          !end if
          !do jorb=mad%kernel_segkeyg(1,iseg,ii)+m,mad%kernel_segkeyg(2,iseg,ii),7
          !    c(jorb+0,i)=0.d0
          !    c(jorb+1,i)=0.d0
          !    c(jorb+2,i)=0.d0
          !    c(jorb+3,i)=0.d0
          !    c(jorb+4,i)=0.d0
          !    c(jorb+5,i)=0.d0
          !    c(jorb+6,i)=0.d0
          !end do
      end do
  end do


      !!$omp parallel default(private) shared(norbp, norb, mad, istindexarr, a_seq, b, c) firstprivate (i, iseg)
      tt2=0.d0
      ncount=0.d0
      ipt=0
      do i=iminorb,imaxorb
         iii=i
         do iseg=1,mad%kernel_nseg(iii)
              !!$omp do
              is = (iii-1)*norb+mad%kernel_segkeyg(1,iseg,iii)
              ie = (iii-1)*norb+mad%kernel_segkeyg(2,iseg,iii)
              do iiorb=max(imin,is),min(imax,ie)
                  iorb=mod(iiorb-1,norb)+1
                  ipt=ipt+1
                  ii0=istindexarr(iorb-mad%kernel_segkeyg(1,iseg,iii)+1,iseg,i)
                  ii2=0
                  tt=0.d0
                  ilen=0
                  ii=0
                  do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                      ilen=ilen+mad%keyg(2,jseg)-mad%keyg(1,jseg)+1
                  end do
                  m=mod(ilen,7)
                  if (m/=0) then
                      do jorb=1,m
                  !write(*,'(a,5i9)') 'iproc, ii0, ii2, ivectorindex(ii0+ii2), imin', iproc, ii0, ii2, ivectorindex(ii0+ii2), imin
                         jjorb=ivectorindex(ii0+ii2)
                         tt = tt + b(jjorb)*a_seq(ii0+ii2)
                         ii2=ii2+1
                         ncount=ncount+1.d0
                      end do
                  end if
                  mp1=m+1
                  do jorb=mp1,ilen,7

                  !write(*,'(a,5i9)') 'iproc, ii0, ii2, ivectorindex(ii0+ii2), imin', iproc, ii0, ii2, ivectorindex(ii0+ii2), imin

                     jjorb0=ivectorindex(ii0+ii2+0)
                     tt = tt + b(jjorb0)*a_seq(ii0+ii2+0)

                     jjorb1=ivectorindex(ii0+ii2+1)
                     tt = tt + b(jjorb1)*a_seq(ii0+ii2+1)

                     jjorb2=ivectorindex(ii0+ii2+2)
                     tt = tt + b(jjorb2)*a_seq(ii0+ii2+2)

                     jjorb3=ivectorindex(ii0+ii2+3)
                     tt = tt + b(jjorb3)*a_seq(ii0+ii2+3)

                     jjorb4=ivectorindex(ii0+ii2+4)
                     tt = tt + b(jjorb4)*a_seq(ii0+ii2+4)

                     jjorb5=ivectorindex(ii0+ii2+5)
                     tt = tt + b(jjorb5)*a_seq(ii0+ii2+5)

                     jjorb6=ivectorindex(ii0+ii2+6)
                     tt = tt + b(jjorb6)*a_seq(ii0+ii2+6)

                     ii2=ii2+7
                     ncount=ncount+7.d0
                  end do
                  c(ipt)=tt
              end do
              !!$omp end do
         end do
      end do 
      !!$omp end parallel

    
end subroutine sparsemm2



subroutine init_onedimindices(norbp, isorb, mad, nout, onedimindices)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb
  type(matrixDescriptors),intent(in) :: mad
  integer,intent(out) :: nout
  integer,dimension(:,:),pointer :: onedimindices

  ! Local variables
  integer :: i, iii, iseg, iorb, istat, ii, jseg, ilen, itot
  character(len=*),parameter :: subname='init_onedimindices'


  nout=0
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,mad%kernel_nseg(iii)
          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
              nout=nout+1
          end do
      end do
  end do

  allocate(onedimindices(4,nout), stat=istat)
  call memocc(istat, onedimindices, 'onedimindices', subname)

  ii=0
  itot=1
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,mad%kernel_nseg(iii)
          do iorb=mad%kernel_segkeyg(1,iseg,iii),mad%kernel_segkeyg(2,iseg,iii)
              ii=ii+1
              onedimindices(1,ii)=i
              onedimindices(2,ii)=iorb
              ilen=0
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  ilen=ilen+mad%keyg(2,jseg)-mad%keyg(1,jseg)+1
              end do
              onedimindices(3,ii)=ilen
              onedimindices(4,ii)=itot
              itot=itot+ilen
          end do
      end do
  end do

end subroutine init_onedimindices

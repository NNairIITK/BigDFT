subroutine chebyshev(iproc, nproc, npl, cc, tmb, ham_compr, ovrlp_compr, fermi, penalty_ev)
  use module_base
  use module_types
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
  integer :: isegstart, isegend, iseg, ii, jjorb
  character(len=*),parameter :: subname='chebyshev'
  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp, t1_tmp2
  real(kind=8),dimension(:),allocatable :: ham_compr_seq, ovrlp_compr_seq, SHS, SHS_seq
  real(kind=8),dimension(:,:),allocatable :: penalty_ev_seq, matrix
  real(kind=8) :: tt1, tt2, time1,time2 , tt, time_to_zero, time_vcopy, time_sparsemm, time_axpy, time_axbyz, time_copykernel
  integer,dimension(:,:,:),allocatable :: istindexarr
  integer,dimension(:),allocatable :: ivectorindex
  integer,parameter :: one=1, three=3
  integer,parameter :: number_of_matmuls=one

  call timing(iproc, 'chebyshev_comp', 'ON')


  norb = tmb%orbs%norb
  norbp = tmb%orbs%norbp
  isorb = tmb%orbs%isorb

  call determine_sequential_length(norbp, norb, tmb%mad, nseq, nmaxsegk, nmaxvalk)

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

  call enable_sequential_acces_matrix(norbp, norb, tmb%mad, ham_compr, nseq, nmaxsegk, nmaxvalk, ham_compr_seq, istindexarr, ivectorindex)
  call enable_sequential_acces_matrix(norbp, norb, tmb%mad, ovrlp_compr, nseq, nmaxsegk, nmaxvalk, ovrlp_compr_seq, istindexarr, ivectorindex)

  allocate(column(norb,norbp), stat=istat)
  call memocc(istat, column, 'column', subname)
  call to_zero(norb*norbp, column(1,1))

  if (number_of_matmuls==one) then

      call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, matrix(1,1), column, norb, norbp, tmb%mad, ivectorindex)
      call to_zero(norbp*norb, matrix(1,1))
      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column, matrix(1,1), norb, norbp, tmb%mad, ivectorindex)
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

      call enable_sequential_acces_matrix(norbp, norb, tmb%mad, SHS, nseq, nmaxsegk, nmaxvalk, SHS_seq, istindexarr, ivectorindex)

  end if


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
      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, norb, norbp, tmb%mad, ivectorindex)
      call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column, column_tmp, norb, norbp, tmb%mad, ivectorindex)
      call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, norb, norbp, tmb%mad, ivectorindex)
  else if (number_of_matmuls==one) then
      call sparsemm(nseq, SHS_seq, nmaxsegk, nmaxvalk, istindexarr, column_tmp, column, norb, norbp, tmb%mad, ivectorindex)
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
  !call to_zero(norbp*norb, fermider(1,isorb+1))
  call to_zero(2*norb*norbp, penalty_ev(1,1,1))
tt2=mpi_wtime() 
time_to_zero=time_to_zero+tt2-tt1

tt1=mpi_wtime() 
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,1), t(1,1), fermi(:,1))
  !call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,2), t(1,1), fermider(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,2))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,1), t1(1,1), fermi(:,1))
  !call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,2), t1(1,1), fermider(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,3), t1(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, -cc(2,3), t1(1,1), penalty_ev(:,1,2))
tt2=mpi_wtime() 
time_axpy=time_axpy+tt2-tt1

  !call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)


  do ipl=3,npl
     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
tt1=mpi_wtime() 
     if (number_of_matmuls==three) then
         call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp, t1, norb, norbp, tmb%mad, ivectorindex)
         call sparsemm(nseq, ham_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1, t1_tmp2, norb, norbp, tmb%mad, ivectorindex)
         call sparsemm(nseq, ovrlp_compr_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp2, t1, norb, norbp, tmb%mad, ivectorindex)
     else if (number_of_matmuls==one) then
         call sparsemm(nseq, SHS_seq, nmaxsegk, nmaxvalk, istindexarr, t1_tmp, t1, norb, norbp, tmb%mad, ivectorindex)
     end if
tt2=mpi_wtime() 
time_sparsemm=time_sparsemm+tt2-tt1
     !call daxbyz(norb*norbp, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
tt1=mpi_wtime() 
     call axbyz_kernel_vectors(norbp, norb, tmb%mad, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
tt2=mpi_wtime() 
time_axbyz=time_axbyz+tt2-tt1
tt1=mpi_wtime() 
     call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,1), t2(1,1), fermi(:,1))
     !call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,2), t2(1,1), fermider(:,isorb+1))
     call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,3), t2(1,1), penalty_ev(:,1,1))
     if (mod(ipl,2)==1) then
         tt=cc(ipl,3)
     else
         tt=-cc(ipl,3)
     end if
     call axpy_kernel_vectors(norbp, norb, tmb%mad, tt, t2(1,1), penalty_ev(:,1,2))
tt2=mpi_wtime() 
time_axpy=time_axpy+tt2-tt1

     !update t's
tt1=mpi_wtime() 
     call copy_kernel_vectors(norbp, norb, tmb%mad, t1_tmp, t)
     call copy_kernel_vectors(norbp, norb, tmb%mad, t2, t1_tmp)
tt2=mpi_wtime() 
time_copykernel=time_copykernel+tt2-tt1
 end do
 
  call timing(iproc, 'chebyshev_comp', 'OF')
  if(iproc==0) write(*,'(a,es16.7)') 'time_to_zero', time_to_zero
  if(iproc==0) write(*,'(a,es16.7)') 'time_vcopy', time_vcopy
  if(iproc==0) write(*,'(a,es16.7)') 'time_sparsemm', time_sparsemm
  if(iproc==0) write(*,'(a,es16.7)') 'time_axpy', time_axpy
  if(iproc==0) write(*,'(a,es16.7)') 'time_axbyz', time_axbyz
  if(iproc==0) write(*,'(a,es16.7)') 'time_copykernel', time_copykernel


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

end subroutine chebyshev




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
subroutine axbyz_kernel_vectors(norbp, norb, mad, a, x, b, y, z)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb
  type(matrixDescriptors),intent(in) :: mad
  real(8),intent(in) :: a, b
  real(kind=8),dimension(norb,norbp),intent(in) :: x, y
  real(kind=8),dimension(norb,norbp),intent(out) :: z

  ! Local variables
  integer :: i, m, mp1, jorb, iseg


  do i=1,norbp
      do iseg=1,mad%kernel_nseg(i)
          m=mod(mad%kernel_segkeyg(2,iseg,i)-mad%kernel_segkeyg(1,iseg,i)+1,7)
          if (m/=0) then
              do jorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(1,iseg,i)+m-1
                  z(jorb,i)=a*x(jorb,i)+b*y(jorb,i)
              end do
          end if
          do jorb=mad%kernel_segkeyg(1,iseg,i)+m,mad%kernel_segkeyg(2,iseg,i),7
              z(jorb+0,i)=a*x(jorb+0,i)+b*y(jorb+0,i)
              z(jorb+1,i)=a*x(jorb+1,i)+b*y(jorb+1,i)
              z(jorb+2,i)=a*x(jorb+2,i)+b*y(jorb+2,i)
              z(jorb+3,i)=a*x(jorb+3,i)+b*y(jorb+3,i)
              z(jorb+4,i)=a*x(jorb+4,i)+b*y(jorb+4,i)
              z(jorb+5,i)=a*x(jorb+5,i)+b*y(jorb+5,i)
              z(jorb+6,i)=a*x(jorb+6,i)+b*y(jorb+6,i)
          end do
      end do
  end do

end subroutine axbyz_kernel_vectors






subroutine sparsemm(nseq, a_seq, nmaxsegk, nmaxvalk, istindexarr, b, c, norb, norbp, mad, ivectorindex)

use module_base
use module_types

  implicit none

  !Calling Arguments
  type(matrixDescriptors),intent(in) :: mad
  integer, intent(in) :: norb,norbp,nseq,nmaxsegk,nmaxvalk
  real(kind=8), dimension(norb,norbp),intent(in) :: b
  real(kind=8), dimension(nseq),intent(in) :: a_seq
  integer,dimension(nmaxvalk,nmaxsegk,norbp),intent(in) :: istindexarr
  real(kind=8), dimension(norb,norbp), intent(out) :: c
  integer,dimension(nseq),intent(in) :: ivectorindex

  !Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ii0, ii2, is, ie, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6
  integer,dimension(:), allocatable :: n
  character(len=*),parameter :: subname='sparsemm'
  real(8) :: ncount, t1, t2, ddot, tt, tt2, t3, ncount2
  integer :: iproc, ierr, imin, imax, nthreads
  integer,dimension(2,5000) :: istarr
  integer,dimension(5000) :: istarr2
  !$ integer  :: omp_get_max_threads

  !!call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)

  do i=1,norbp
      do iseg=1,mad%kernel_nseg(i)
          m=mod(mad%kernel_segkeyg(2,iseg,i)-mad%kernel_segkeyg(1,iseg,i)+1,7)
          if (m/=0) then
              do jorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(1,iseg,i)+m-1
                  c(jorb,i)=0.d0
              end do
          end if
          do jorb=mad%kernel_segkeyg(1,iseg,i)+m,mad%kernel_segkeyg(2,iseg,i),7
              c(jorb+0,i)=0.d0
              c(jorb+1,i)=0.d0
              c(jorb+2,i)=0.d0
              c(jorb+3,i)=0.d0
              c(jorb+4,i)=0.d0
              c(jorb+5,i)=0.d0
              c(jorb+6,i)=0.d0
          end do
      end do
  end do


t1=mpi_wtime()


  !!nthreads=1
  !!!$  nthreads = omp_get_max_threads()

  !!if (mod(norbp,nthreads)==0) then

  !!    !$omp parallel default(private) shared(norbp, norb, mad, istindexarr, a_seq, b, c)
  !!    !$omp do
  !!    do i = 1,norbp
  !!       do iseg=1,mad%kernel_nseg(i)
  !!            do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
  !!                ii0=istindexarr(iorb-mad%kernel_segkeyg(1,iseg,i)+1,iseg,i)
  !!                ii2=0
  !!                do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
  !!                    do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
  !!                        jjorb = jorb - (iorb-1)*norb
  !!                        c(iorb,i) = c(iorb,i) + b(jjorb,i)*a_seq(ii0+ii2)
  !!                        ii2=ii2+1
  !!                    end do
  !!                end do
  !!            end do
  !!       end do
  !!    end do 
  !!    !$omp end do
  !!    !$omp end parallel

  !!else

t1=mpi_wtime()

      !!$omp parallel default(private) shared(norbp, norb, mad, istindexarr, a_seq, b, c) firstprivate (i, iseg)
      tt2=0.d0
      ncount=0.d0
      do i = 1,norbp
         do iseg=1,mad%kernel_nseg(i)
              !!$omp do
              do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
                  ii0=istindexarr(iorb-mad%kernel_segkeyg(1,iseg,i)+1,iseg,i)
                  ii2=0
                  tt=0.d0
                  ilen=0
                  ii=0
                  do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                      ilen=ilen+mad%keyg(2,jseg)-mad%keyg(1,jseg)+1
                      !!istarr(1,jseg-mad%istsegline(iorb)+1)=ii
                      !!istarr(2,jseg-mad%istsegline(iorb)+1)=mad%keyg(1,jseg)-(iorb-1)*norb
                      !!ii=ii+mad%keyg(2,jseg)-mad%keyg(1,jseg)+1
                  end do
                  !!ii=1
                  !!$$do jorb=1,ilen
                  !!$$!do jorb=1,500
                  !!$$   !jjorb=istarr(2,ii)+jorb
                  !!$$   jjorb=ivectorindex(ii0+ii2)
                  !!$$   tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
                  !!$$   !tt=tt+dble(jjorb)*dble(ii2) 
                  !!$$   ii2=ii2+1
                  !!$$   ncount=ncount+1.d0
                  !!$$   !if (jorb>istarr(1,ii)) ii=ii+1
                  !!$$end do
                  m=mod(ilen,7)
                  if (m/=0) then
                      do jorb=1,m
                         jjorb=ivectorindex(ii0+ii2)
                         tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
                         ii2=ii2+1
                         ncount=ncount+1.d0
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
                     ncount=ncount+7.d0
                  end do

                  !!do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  !!!do jseg=1,100
                  !!    is=mad%keyg(1,jseg)-(iorb-1)*norb
                  !!    ie=mad%keyg(2,jseg)-(iorb-1)*norb
                  !!    !ilen=mad%keyg(2,jseg)-mad%keyg(1,jseg)+1
                  !!    do jjorb = is,ie
                  !!    !do jorb=1,ilen
                  !!        !jjorb=jorb-(iorb-1)*norb
                  !!    !!do jjorb = 1,10
                  !!        tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
                  !!        !tt=tt+dble(jjorb)*dble(ii2)
                  !!        ii2=ii2+1
                  !!        ncount=ncount+1.d0
                  !!    end do
                  !!end do
                  c(iorb,i)=tt
                  !tt2=tt2+tt
                  !temparr(iorb)=tt
              end do
              !!$omp end do
         end do
      end do 
      !!$omp end parallel
      !c(1,1)=tt2

t2=mpi_wtime()
      
  !!tt2=0.d0
  !!ncount2=0.d0
  !!do iall=1,33
  !!    do istat=1,33
  !!        do i=1,nint(ncount/100000.d0)
  !!            do iseg=1,100
  !!                tt2=tt2+dble(i+istat)*dble(i+iseg-iall)
  !!                tt=tt+1.d0
  !!                ncount2=ncount2+1.d0
  !!            end do
  !!        end do
  !!    end do
  !!end do
  !!tt=tt2

t3=mpi_wtime()

!call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
!if (iproc==0) write(*,'(a,2es12.1,2es14.4,10x,es12.2)') 'ncount, ncount2, times, tt', ncount, ncount2, t2-t1, t3-t2, tt


  !!end if

  !!!!ncount=0.d0
  !!!$omp parallel default(private) shared(mad, c, b, a, norb, norbp) firstprivate (i, iseg)
  !!do i = 1,norbp
  !!   do iseg=1,mad%kernel_nseg(i)
  !!        !$omp do
  !!        do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
  !!            do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
  !!                !!c(iorb,i) = c(iorb,i) + &
  !!                !!  ddot(mad%keyg(2,jseg)-mad%keyg(1,jseg)+1, b(mad%keyg(1,jseg)-(iorb-1)*norb,i), 1, a(mad%keyv(jseg)), 1)
  !!                jj=1
  !!                !!m = mod(mad%keyg(2,jseg)-mad%keyg(1,jseg)+1,7)
  !!                !!if (m/=0) then
  !!                    !do jorb = mad%keyg(1,jseg),mad%keyg(1,jseg)+m-1
  !!                    do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
  !!                        jjorb = jorb - (iorb-1)*norb
  !!                        c(iorb,i) = c(iorb,i) + b(jjorb,i)*a(mad%keyv(jseg)+jj-1)
  !!                        jj = jj+1
  !!                        !!ncount=ncount+1.d0
  !!                    end do
  !!                !!end if
  !!                !!do jorb = mad%keyg(1,jseg)+m,mad%keyg(2,jseg),7
  !!                !!    jjorb = jorb - (iorb-1)*norb
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+0,i)*a(mad%keyv(jseg)+jj+0-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+1,i)*a(mad%keyv(jseg)+jj+1-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+2,i)*a(mad%keyv(jseg)+jj+2-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+3,i)*a(mad%keyv(jseg)+jj+3-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+4,i)*a(mad%keyv(jseg)+jj+4-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+5,i)*a(mad%keyv(jseg)+jj+5-1)
  !!                !!    c(iorb,i) = c(iorb,i) + b(jjorb+6,i)*a(mad%keyv(jseg)+jj+6-1)
  !!                !!    jj = jj+7
  !!                !!    !!ncount=ncount+7.d0
  !!                !!end do
  !!            end do
  !!        end do
  !!        !$omp end do
  !!   end do
  !!end do 
  !!!$omp end parallel


  !!imin=1000000000
  !!imax=-1000000000
  !!do i=1,norbp
  !!    do iseg=1,mad%kernel_nseg(i)
  !!        imin=min(imin,mad%kernel_segkeyg(1,iseg,i))
  !!        imax=max(imax,mad%kernel_segkeyg(2,iseg,i))
  !!    end do
  !!end do

  !!!$omp parallel default(private) shared(imin, imax, mad, norbp, norb, a, b, c)
  !!
  !!!ncount=0.d0
  !!!$omp do
  !!do iorb=imin,imax
  !!    do i=1,norbp
  !!        if (.not.mad%kernel_locreg(i,iorb)) cycle
  !!        do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
  !!            !!jj=1
  !!            !!m = mod(mad%keyg(2,jseg)-mad%keyg(1,jseg)+1,2)
  !!            !!if (m/=0) then
  !!                !!do jorb = mad%keyg(1,jseg),mad%keyg(1,jseg)+m-1
  !!                do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
  !!                    jj=jorb-mad%keyg(1,jseg)
  !!                    jjorb = jorb - (iorb-1)*norb
  !!                    c(iorb,i) = c(iorb,i) + b(jjorb,i)*a(mad%keyv(jseg)+jj)
  !!                    !jj = jj+1
  !!                end do
  !!            !!end if
  !!            !!do jorb = mad%keyg(1,jseg)+m,mad%keyg(2,jseg),2
  !!            !!    jj=jorb-mad%keyg(1,jseg)-m
  !!            !!    jjorb = jorb - (iorb-1)*norb
  !!            !!    c(iorb,i) = c(iorb,i) + b(jjorb+0,i)*a(mad%keyv(jseg)+jj+0)
  !!            !!    c(iorb,i) = c(iorb,i) + b(jjorb+1,i)*a(mad%keyv(jseg)+jj+1)
  !!            !!    !!c(iorb,i) = c(iorb,i) + b(jjorb+2,i)*a(mad%keyv(jseg)+jj+2)
  !!            !!    !!c(iorb,i) = c(iorb,i) + b(jjorb+3,i)*a(mad%keyv(jseg)+jj+3)
  !!            !!    !!c(iorb,i) = c(iorb,i) + b(jjorb+4,i)*a(mad%keyv(jseg)+jj+4-1)
  !!            !!    !!c(iorb,i) = c(iorb,i) + b(jjorb+5,i)*a(mad%keyv(jseg)+jj+5-1)
  !!            !!    !!c(iorb,i) = c(iorb,i) + b(jjorb+6,i)*a(mad%keyv(jseg)+jj+6-1)
  !!            !!    !jj = jj+4
  !!            !!    !ncount=ncount+1.d0
  !!            !!end do
  !!        end do
  !!    end do
  !!end do
  !!!$omp end do
  !!!$omp end parallel





t2=mpi_wtime()


  !!call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
  !!write(10000+iproc,'(a,2es16.6)') 'ncount, time', ncount, t2-t1
  !write(10000+iproc,'(a,3es16.6)') 'ncount, time, tt', ncount, t2-t1, tt
    
end subroutine sparsemm



subroutine copy_kernel_vectors(norbp, norb, mad, a, b)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(norb,norbp),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(out) :: b

  ! Local variables
  integer :: i, m, jorb, mp1, iseg


  do i=1,norbp
      do iseg=1,mad%kernel_nseg(i)
          m=mod(mad%kernel_segkeyg(2,iseg,i)-mad%kernel_segkeyg(1,iseg,i)+1,7)
          if (m/=0) then
              do jorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(1,iseg,i)+m-1
                  b(jorb,i)=a(jorb,i)
              end do
          end if
          do jorb=mad%kernel_segkeyg(1,iseg,i)+m,mad%kernel_segkeyg(2,iseg,i),7
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




subroutine axpy_kernel_vectors(norbp, norb, mad, a, x, y)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(in) :: x
  real(kind=8),dimension(norb,norbp),intent(out) :: y

  ! Local variables
  integer :: i, m, jorb, mp1, iseg

  do i=1,norbp
      do iseg=1,mad%kernel_nseg(i)
          m=mod(mad%kernel_segkeyg(2,iseg,i)-mad%kernel_segkeyg(1,iseg,i)+1,7)
          if (m/=0) then
              do jorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(1,iseg,i)+m-1
                  y(jorb,i)=y(jorb,i)+a*x(jorb,i)
              end do
          end if
          do jorb=mad%kernel_segkeyg(1,iseg,i)+m,mad%kernel_segkeyg(2,iseg,i),7
              y(jorb+0,i)=y(jorb+0,i)+a*x(jorb+0,i)
              y(jorb+1,i)=y(jorb+1,i)+a*x(jorb+1,i)
              y(jorb+2,i)=y(jorb+2,i)+a*x(jorb+2,i)
              y(jorb+3,i)=y(jorb+3,i)+a*x(jorb+3,i)
              y(jorb+4,i)=y(jorb+4,i)+a*x(jorb+4,i)
              y(jorb+5,i)=y(jorb+5,i)+a*x(jorb+5,i)
              y(jorb+6,i)=y(jorb+6,i)+a*x(jorb+6,i)
          end do
      end do
  end do

end subroutine axpy_kernel_vectors




subroutine enable_sequential_acces_matrix(norbp, norb, mad, a, nseq, nmaxsegk, nmaxvalk, a_seq, istindexarr, ivectorindex)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb, nseq, nmaxsegk, nmaxvalk
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(mad%nvctr),intent(in) :: a
  real(kind=8),dimension(nseq),intent(out) :: a_seq
  integer,dimension(nmaxvalk,nmaxsegk,norbp),intent(out) :: istindexarr
  integer,dimension(nseq),intent(out) :: ivectorindex

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb


  ii=1
  do i = 1,norbp
     do iseg=1,mad%kernel_nseg(i)
          do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
              istindexarr(iorb-mad%kernel_segkeyg(1,iseg,i)+1,iseg,i)=ii
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



subroutine enable_sequential_acces_vector(norbp, norb, mad, b, nseq, b_seq)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb, nseq
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(norb,norbp),intent(in) :: b
  real(kind=8),dimension(nseq),intent(out) :: b_seq

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb

  ii=1
  do i = 1,norbp
     do iseg=1,mad%kernel_nseg(i)
          do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  jj=0
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      jjorb = jorb - (iorb-1)*norb
                      b_seq(ii)=b(jjorb,i)
                      jj = jj+1
                      ii = ii+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine enable_sequential_acces_vector



subroutine determine_sequential_length(norbp, norb, mad, nseq, nmaxsegk, nmaxvalk)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb
  type(matrixDescriptors),intent(in) :: mad
  integer,intent(out) :: nseq, nmaxsegk, nmaxvalk

  ! Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer :: iorb, jseg, ii, ist_iorb

  nseq=0
  nmaxsegk=0
  nmaxvalk=0
  do i = 1,norbp
     nmaxsegk=max(nmaxsegk,mad%kernel_nseg(i))
     do iseg=1,mad%kernel_nseg(i)
          nmaxvalk=max(nmaxvalk,mad%kernel_segkeyg(2,iseg,i)-mad%kernel_segkeyg(1,iseg,i)+1)
          do iorb=mad%kernel_segkeyg(1,iseg,i),mad%kernel_segkeyg(2,iseg,i)
              do jseg=mad%istsegline(iorb),mad%istsegline(iorb)+mad%nsegline(iorb)-1
                  do jorb = mad%keyg(1,jseg),mad%keyg(2,jseg)
                      nseq=nseq+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine determine_sequential_length

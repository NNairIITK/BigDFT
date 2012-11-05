subroutine chebyshev(iproc, nproc, npl, cc, tmb, ham, ovrlp, fermi, fermider, penalty_ev)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl
  real(8),dimension(npl,3),intent(in) :: cc
  type(DFT_wavefunction),intent(in) :: tmb 
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in) :: ham, ovrlp
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out) :: fermi, fermider
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb,2),intent(out) :: penalty_ev
  ! Local variables
  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb, ierr,i,j
  character(len=*),parameter :: subname='chebyshev'
  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp, t1_tmp2, ts
  real(kind=8), dimension(tmb%orbs%norb,tmb%orbs%norb) :: ovrlp_tmp,ham_eff
  real(kind=8),dimension(:),allocatable :: ovrlp_compr, ham_compr
  real(kind=8) :: time1,time2 , tt
  norb = tmb%orbs%norb
  norbp = tmb%orbs%norbp
  isorb = tmb%orbs%isorb


  allocate(column(norb,norbp), stat=istat)
  call memocc(istat, column, 'column', subname)
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
  allocate(ts(norb,norbp), stat=istat)
  call memocc(istat, ts, 'ts', subname)
  allocate(ovrlp_compr(tmb%mad%nvctr), stat=istat)
  call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)
  allocate(ham_compr(tmb%mad%nvctr), stat=istat)
  call memocc(istat, ham_compr, 'ham_compr', subname)
 

  call timing(iproc, 'chebyshev_comp', 'ON')

 
  call to_zero(norb*norbp, column(1,1))
  call to_zero(norb*norbp, column_tmp(1,1))
  do iorb=1,norbp
      iiorb=isorb+iorb
      column(iiorb,iorb)=1.d0
  end do

  call vcopy(norb*norbp, column(1,1), 1, column_tmp(1,1), 1)
  !column_tmp = column

  do iorb=1,norb
      do jorb=1,norb
          if (iorb==jorb) then
              ovrlp_tmp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
          else
              ovrlp_tmp(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
          end if
          !!write(3000+iproc,'(2i7,2es18.7)') iorb,jorb,ovrlp_tmp(jorb,iorb), ham(jorb,iorb)
      end do
  end do


  call compress_matrix_for_allreduce(norb, tmb%mad, ovrlp_tmp, ovrlp_compr)
  call compress_matrix_for_allreduce(norb, tmb%mad, ham, ham_compr)
  ! t0
  !t = column
  call vcopy(norb*norbp, column(1,1), 1, t(1,1), 1)

  !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) column
  !call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,column_tmp,norb,0.d0,column,norb)
  !call dsymm('L', 'U', norb, norbp,1.d0,ham,norb,column,norb,0.d0,column_tmp,norb)
  !call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,column_tmp,norb,0.d0,column,norb)
  call sparsemm(ovrlp_compr, column_tmp, column, norb, norbp, tmb%mad)
  call sparsemm(ham_compr, column, column_tmp, norb, norbp, tmb%mad)
  call sparsemm(ovrlp_compr, column_tmp, column, norb, norbp, tmb%mad)


  !t1 = column
  call vcopy(norb*norbp, column(1,1), 1, t1(1,1), 1)
  !t1_tmp = t1
  call vcopy(norb*norbp, t1(1,1), 1, t1_tmp(1,1), 1)

  !!if (iproc==0) write(*,'(a,2es18.7)') 't(1,1), t1(1,1)', t(1,1), t1(1,1)
  !initialize fermi
  call to_zero(norb*norb, fermi(1,1))
  call to_zero(norb*norb, fermider(1,1))
  call to_zero(2*norb*norb, penalty_ev(1,1,1))
  !do iorb = 1,norbp
     !!iiorb = isorb + iorb
     !!fermi(:,isorb+iorb) = cc(1,1)*0.5d0*t(:,iorb) + cc(2,1)*t1(:,iorb)
     !!fermider(:,isorb+iorb) = cc(1,2)*0.5d0*t(:,iorb) + cc(2,2)*t1(:,iorb)
     !!penalty_ev(:,isorb+iorb,1) = cc(1,3)*0.5d0*t(:,iorb) + cc(2,3)*t1(:,iorb)
     !!penalty_ev(:,isorb+iorb,2) = cc(1,3)*0.5d0*t(:,iorb) - cc(2,3)*t1(:,iorb)
     call daxpy(norb*norbp, 0.5d0*cc(1,1), t(1,1), 1, fermi(:,isorb+1), 1)
     call daxpy(norb*norbp, 0.5d0*cc(1,2), t(1,1), 1, fermider(:,isorb+1), 1)
     call daxpy(norb*norbp, 0.5d0*cc(1,3), t(1,1), 1, penalty_ev(:,isorb+1,1), 1)
     call daxpy(norb*norbp, 0.5d0*cc(1,3), t(1,1), 1, penalty_ev(:,isorb+1,2), 1)
     call daxpy(norb*norbp, cc(2,1), t1(1,1), 1, fermi(:,isorb+1), 1)
     call daxpy(norb*norbp, cc(2,2), t1(1,1), 1, fermider(:,isorb+1), 1)
     call daxpy(norb*norbp, cc(2,3), t1(1,1), 1, penalty_ev(:,isorb+1,1), 1)
     call daxpy(norb*norbp, -cc(2,3), t1(1,1), 1, penalty_ev(:,isorb+1,2), 1)
  !end do

  !!time1 = MPI_WTIME() 
  !!call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp,norb,0.d0,t1,norb)
  !!time2 = MPI_WTIME()
  !!write(100,*) time2 - time1  
  !!time1= MPI_WTIME()
  call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)
  !!time2 = MPI_WTIME()
  !!write(200,*) time2 -time1

  !!do i = 1,norb
  !!  do j=1,norbp
  !!     write(315+iproc,*) i,j,t1(i,j)
  !!     write(400+iproc,*) i,j,ts(i,j)
  !!  end do
  !!end do


  do ipl=3,npl
     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
     !!call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp,norb,0.d0,t1,norb)
     !!call dsymm('L', 'U', norb, norbp,1.d0,ham,norb,t1,norb,0.d0,t1_tmp2,norb)
     !!call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp2,norb,0.d0,t1,norb)
     call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)
     call sparsemm(ham_compr, t1, t1_tmp2, norb, norbp, tmb%mad)
     call sparsemm(ovrlp_compr, t1_tmp2, t1, norb, norbp, tmb%mad)
     !calculate t2 = 2 * (3/2 - 1/2 S) H (3/2 - 1/2 S) t1 - t
     !t2 = 2*t1 - t
     call daxbyz(norb*norbp, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
     !call daxpy(norb*norbp, 2.d0, t1(1,1), 1, t2(1,1), 1)
     !call daxpy(norb*norbp, -1.d0, t(1,1), 1, t2(1,1), 1)
     !update fermi-matrix
     !!fermi(:,isorb+1:isorb+norbp)=fermi(:,isorb+1:isorb+norbp) + cc(ipl,1)*t2   
     !!fermider(:,isorb+1:isorb+norbp)=fermider(:,isorb+1:isorb+norbp) + cc(ipl,2)*t2   
     !!penalty_ev(:,isorb+1:isorb+norbp,1)=penalty_ev(:,isorb+1:isorb+norbp,1) + cc(ipl,3)*t2   
     call daxpy(norb*norbp, cc(ipl,1), t2(1,1), 1, fermi(:,isorb+1), 1)
     call daxpy(norb*norbp, cc(ipl,2), t2(1,1), 1, fermider(:,isorb+1), 1)
     call daxpy(norb*norbp, cc(ipl,3), t2(1,1), 1, penalty_ev(:,isorb+1,1), 1)
     if (mod(ipl,2)==1) then
         tt=cc(ipl,3)
     else
         tt=-cc(ipl,3)
     end if
     !penalty_ev(:,isorb+1:isorb+norbp,2)=penalty_ev(:,isorb+1:isorb+norbp,2) + tt*t2   
     call daxpy(norb*norbp, tt, t2(1,1), 1, penalty_ev(:,isorb+1,2), 1)

     !update t's
     !t = t1_tmp
     call vcopy(norb*norbp, t1_tmp(1,1), 1, t(1,1), 1)
     !t1_tmp = t2
     call vcopy(norb*norbp, t2(1,1), 1, t1_tmp(1,1), 1)

     !!do iorb=1,norb
     !!    write(1000+iproc,*) ipl, t2(iorb,isorb+1)
     !!end do
 end do
 
  call timing(iproc, 'chebyshev_comp', 'OF')
 !! call timing(iproc, 'chebyshev_comm', 'ON')

 !!call mpiallred(fermi(1,1), norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
 !!call mpiallred(fermider(1,1), norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
 !!call mpiallred(penalty_ev(1,1,1), 2*norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)

 !! call timing(iproc, 'chebyshev_comm', 'OF')


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
  iall=-product(shape(ts))*kind(ts)
  deallocate(ts, stat=istat)
  call memocc(istat, iall, 'ts', subname)
  iall=-product(shape(ovrlp_compr))*kind(ovrlp_compr)
  deallocate(ovrlp_compr, stat=istat)
  call memocc(istat, iall, 'ovrlp_compr', subname)
  iall=-product(shape(ham_compr))*kind(ham_compr)
  deallocate(ham_compr, stat=istat)
  call memocc(istat, iall, 'ham_compr', subname)

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





subroutine sparsemm(a,b,c,norb,norbp,mad)

use module_base
use module_types


  implicit none

  !Calling Arguments
  type(matrixDescriptors),intent(in) :: mad
  integer, intent(in) :: norb,norbp
  real(kind=8), dimension(norb,norbp),intent(in) :: b
  real(kind=8), dimension(mad%nvctr),intent(in) :: a
  real(kind=8), dimension(norb,norbp), intent(out) :: c

  !Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m
  real(kind=8) :: temp

  !write(*,*) 'mad%nseg',mad%nseg


  call to_zero(norb*norbp,c(1,1))
  do i = 1,norbp
     do iseg = 1,mad%nseg
          jj = 1
          m = mod(mad%keyg(2,iseg)-mad%keyg(1,iseg)+1,4)
          iiorb = mad%keyg(1,iseg)/norb + 1
          !!if (mad%kernel_locreg(iiorb,i)) then
              if(m.ne.0) then
                do jorb = mad%keyg(1,iseg),mad%keyg(1,iseg)+m-1 
                  jjorb = jorb - (iiorb-1)*norb
                  c(iiorb,i) = c(iiorb,i) + b(jjorb,i)*a(mad%keyv(iseg)+jj-1)
                  jj = jj+1
                 end do
              end if

     

              do jorb = mad%keyg(1,iseg)+m, mad%keyg(2,iseg),4
                jjorb = jorb - (iiorb - 1)*norb
                c(iiorb,i) = c(iiorb,i) + b(jjorb,i)*a(mad%keyv(iseg)+jj-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+1,i)*a(mad%keyv(iseg)+jj+1-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+2,i)*a(mad%keyv(iseg)+jj+2-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+3,i)*a(mad%keyv(iseg)+jj+3-1)
                jj = jj + 4
              end do
          !!end if
     end do
  end do 
  
    
end subroutine sparsemm

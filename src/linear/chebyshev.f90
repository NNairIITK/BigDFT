subroutine chebyshev(iproc, nproc, npl, cc, tmb, ham_compr, ovrlp_compr, fermi, fermider, penalty_ev)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl
  real(8),dimension(npl,3),intent(in) :: cc
  type(DFT_wavefunction),intent(in) :: tmb 
  real(kind=8),dimension(tmb%mad%nvctr),intent(in) :: ham_compr, ovrlp_compr
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out) :: fermi, fermider
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norbp,2),intent(out) :: penalty_ev
  ! Local variables
  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb, ierr,i,j
  character(len=*),parameter :: subname='chebyshev'
  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp, t1_tmp2
  real(kind=8), dimension(tmb%orbs%norb,tmb%orbs%norb) :: ovrlp_tmp,ham_eff
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

  call timing(iproc, 'chebyshev_comp', 'ON')

 
  call to_zero(norb*norbp, column(1,1))
  do iorb=1,norbp
      iiorb=isorb+iorb
      column(iiorb,iorb)=1.d0
  end do


  call to_zero(norb*norbp, t1_tmp2(1,1))

  call vcopy(norb*norbp, column(1,1), 1, column_tmp(1,1), 1)
  call vcopy(norb*norbp, column(1,1), 1, t(1,1), 1)

  !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) column
  call sparsemm(ovrlp_compr, column_tmp, column, norb, norbp, tmb%mad)
  call sparsemm(ham_compr, column, column_tmp, norb, norbp, tmb%mad)
  call sparsemm(ovrlp_compr, column_tmp, column, norb, norbp, tmb%mad)


  call vcopy(norb*norbp, column(1,1), 1, t1(1,1), 1)
  call vcopy(norb*norbp, t1(1,1), 1, t1_tmp(1,1), 1)

  !initialize fermi
  call to_zero(norbp*norb, fermi(1,isorb+1))
  call to_zero(norbp*norb, fermider(1,isorb+1))
  call to_zero(2*norb*norbp, penalty_ev(1,1,1))

  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,1), t(1,1), fermi(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,2), t(1,1), fermider(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, 0.5d0*cc(1,3), t(1,1), penalty_ev(:,1,2))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,1), t1(1,1), fermi(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,2), t1(1,1), fermider(:,isorb+1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(2,3), t1(1,1), penalty_ev(:,1,1))
  call axpy_kernel_vectors(norbp, norb, tmb%mad, -cc(2,3), t1(1,1), penalty_ev(:,1,2))

  call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)


  do ipl=3,npl
     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
     call sparsemm(ovrlp_compr, t1_tmp, t1, norb, norbp, tmb%mad)
     call sparsemm(ham_compr, t1, t1_tmp2, norb, norbp, tmb%mad)
     call sparsemm(ovrlp_compr, t1_tmp2, t1, norb, norbp, tmb%mad)
     call daxbyz(norb*norbp, 2.d0, t1(1,1), -1.d0, t(1,1), t2(1,1))
     call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,1), t2(1,1), fermi(:,isorb+1))
     call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,2), t2(1,1), fermider(:,isorb+1))
     call axpy_kernel_vectors(norbp, norb, tmb%mad, cc(ipl,3), t2(1,1), penalty_ev(:,1,1))
     if (mod(ipl,2)==1) then
         tt=cc(ipl,3)
     else
         tt=-cc(ipl,3)
     end if
     call axpy_kernel_vectors(norbp, norb, tmb%mad, tt, t2(1,1), penalty_ev(:,1,2))

     !update t's
     call copy_kernel_vectors(norbp, norb, tmb%mad, t1_tmp, t)
     call copy_kernel_vectors(norbp, norb, tmb%mad, t2, t1_tmp)
 end do
 
  call timing(iproc, 'chebyshev_comp', 'OF')


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
  integer :: i,j,iseg,jorb,iiorb,jjorb,jj,m,istat,iall,nthreads,norbthrd,orb_rest,tid,istart,iend, mp1
  integer,dimension(:), allocatable :: n
  character(len=*),parameter :: subname='sparsemm'
  !$ integer :: OMP_GET_MAX_THREADS, OMP_GET_THREAD_NUM
  nthreads = 1
  !$ nthreads = OMP_GET_MAX_THREADS()

  allocate(n(nthreads),stat=istat)
  call memocc(istat, n, 'n', subname)

  norbthrd = norb/nthreads

  orb_rest = norb - norbthrd*nthreads

  n = norbthrd

  do i = 1,orb_rest
     n(i) = n(i) + 1
  end do

  !!call to_zero(norb*norbp,c(1,1))
  do i=1,norbp
      m = mod(norb,7)
      if (m/=0) then
          do jorb=1,m
              if (mad%kernel_locreg(jorb,i)) then
                  c(jorb,i)=0.d0
              end if
          end do
      end if
      mp1=m+1
      do jorb=mp1,norb,7
          if (mad%kernel_locreg(jorb+0,i)) then
              c(jorb+0,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+1,i)) then
              c(jorb+1,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+2,i)) then
              c(jorb+2,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+3,i)) then
              c(jorb+3,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+4,i)) then
              c(jorb+4,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+5,i)) then
              c(jorb+5,i)=0.d0
          end if
          if (mad%kernel_locreg(jorb+6,i)) then
              c(jorb+6,i)=0.d0
          end if
      end do
  end do

  do i = 1,norbp

     !$OMP parallel default(private) shared(mad,n,norb,a,b,c) firstprivate(i)
     tid = 0
     !$ tid = OMP_GET_THREAD_NUM()
     istart = 1
     iend = 0

     do j = 0,tid-1
        istart = istart + n(j+1)
     end do
     do j = 0,tid
        iend = iend + n(j+1)
     end do

     do iseg = 1,mad%nseg
          jj = 1
          m = mod(mad%keyg(2,iseg)-mad%keyg(1,iseg)+1,7)
          iiorb = mad%keyg(1,iseg)/norb + 1
          if(iiorb < istart .or. iiorb>iend) cycle
          if (mad%kernel_locreg(iiorb,i)) then
              if(m.ne.0) then
                do jorb = mad%keyg(1,iseg),mad%keyg(1,iseg)+m-1 
                  jjorb = jorb - (iiorb-1)*norb
                  c(iiorb,i) = c(iiorb,i) + b(jjorb,i)*a(mad%keyv(iseg)+jj-1)
                  jj = jj+1
                 end do
              end if

              do jorb = mad%keyg(1,iseg)+m, mad%keyg(2,iseg),7
                jjorb = jorb - (iiorb - 1)*norb
                c(iiorb,i) = c(iiorb,i) + b(jjorb,i)*a(mad%keyv(iseg)+jj-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+1,i)*a(mad%keyv(iseg)+jj+1-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+2,i)*a(mad%keyv(iseg)+jj+2-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+3,i)*a(mad%keyv(iseg)+jj+3-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+4,i)*a(mad%keyv(iseg)+jj+4-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+5,i)*a(mad%keyv(iseg)+jj+5-1)
                c(iiorb,i) = c(iiorb,i) + b(jjorb+6,i)*a(mad%keyv(iseg)+jj+6-1)
                jj = jj + 7
              end do
          end if
     end do
     !$OMP end parallel
  end do 


  iall=-product(shape(n))*kind(n)
  deallocate(n, stat=istat)
  call memocc(istat, iall, 'n', subname)
    
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
  integer :: i, m, jorb, mp1

  do i=1,norbp
      m = mod(norb,7)
      if (m/=0) then
          do jorb=1,m
              if (mad%kernel_locreg(jorb,i)) then
                  b(jorb,i)=a(jorb,i)
              end if
          end do
      end if
      mp1=m+1
      do jorb=mp1,norb,7
          if (mad%kernel_locreg(jorb+0,i)) then
              b(jorb+0,i)=a(jorb+0,i)
          end if
          if (mad%kernel_locreg(jorb+1,i)) then
              b(jorb+1,i)=a(jorb+1,i)
          end if
          if (mad%kernel_locreg(jorb+2,i)) then
              b(jorb+2,i)=a(jorb+2,i)
          end if
          if (mad%kernel_locreg(jorb+3,i)) then
              b(jorb+3,i)=a(jorb+3,i)
          end if
          if (mad%kernel_locreg(jorb+4,i)) then
              b(jorb+4,i)=a(jorb+4,i)
          end if
          if (mad%kernel_locreg(jorb+5,i)) then
              b(jorb+5,i)=a(jorb+5,i)
          end if
          if (mad%kernel_locreg(jorb+6,i)) then
              b(jorb+6,i)=a(jorb+6,i)
          end if
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
  integer :: i, m, jorb, mp1

  do i=1,norbp
      m = mod(norb,7)
      if (m/=0) then
          do jorb=1,m
              if (mad%kernel_locreg(jorb,i)) then
                  y(jorb,i)=y(jorb,i)+a*x(jorb,i)
              end if
          end do
      end if
      mp1=m+1
      do jorb=mp1,norb,7
          if (mad%kernel_locreg(jorb+0,i)) then
              y(jorb+0,i)=y(jorb+0,i)+a*x(jorb+0,i)
          end if
          if (mad%kernel_locreg(jorb+1,i)) then
              y(jorb+1,i)=y(jorb+1,i)+a*x(jorb+1,i)
          end if
          if (mad%kernel_locreg(jorb+2,i)) then
              y(jorb+2,i)=y(jorb+2,i)+a*x(jorb+2,i)
          end if
          if (mad%kernel_locreg(jorb+3,i)) then
              y(jorb+3,i)=y(jorb+3,i)+a*x(jorb+3,i)
          end if
          if (mad%kernel_locreg(jorb+4,i)) then
              y(jorb+4,i)=y(jorb+4,i)+a*x(jorb+4,i)
          end if
          if (mad%kernel_locreg(jorb+5,i)) then
              y(jorb+5,i)=y(jorb+5,i)+a*x(jorb+5,i)
          end if
          if (mad%kernel_locreg(jorb+6,i)) then
              y(jorb+6,i)=y(jorb+6,i)+a*x(jorb+6,i)
          end if
      end do
  end do

end subroutine axpy_kernel_vectors

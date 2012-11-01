subroutine chebyshev(iproc, nproc, npl, cc, tmb, ham, ovrlp, fermi)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl
  real(8),dimension(npl),intent(in) :: cc
  type(DFT_wavefunction),intent(in) :: tmb 
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in) :: ham, ovrlp
  real(kind=8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out) :: fermi

  ! Local variables
  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb
  character(len=*),parameter :: subname='chebyshev'
  real(8), dimension(:,:), allocatable :: column,column_tmp, t,t1,t2,t1_tmp,ts
  real(8), dimension(tmb%orbs%norb,tmb%orbs%norb) :: ovrlp_tmp,ham_eff

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
  allocate(t2(norb,norbp), stat=istat)
  call memocc(istat, t2, 't2', subname)
  allocate(ts(norb,norbp), stat=istat)
  call memocc(istat, ts, 'ts', subname)
  
  call to_zero(norb*norbp, column(1,1))
  call to_zero(norb*norbp, column_tmp(1,1))
  do iorb=1,norbp
      iiorb=isorb+iorb
      column(iiorb,iorb)=1.d0
  end do

  column_tmp = column

  do iorb=1,norb
      do jorb=1,norb
          if (iorb==jorb) then
              ovrlp_tmp(jorb,iorb)=1.5d0-.5d0*ovrlp(jorb,iorb)
          else
              ovrlp_tmp(jorb,iorb)=-.5d0*ovrlp(jorb,iorb)
          end if
      end do
  end do


  ! t0
  t = column

  !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) column
  call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,column_tmp,norb,0.d0,column,norb)
  call dsymm('L', 'U', norb, norbp,1.d0,ham,norb,column,norb,0.d0,column_tmp,norb)
  call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,column_tmp,norb,0.d0,column,norb)

  t1 = column
  t1_tmp = t1
  !initialize fermi
  call to_zero(norb*norb, fermi(1,1))
  do iorb = 1,norbp
     iiorb = isorb + iorb
     fermi(:,iorb) = cc(1)*0.5d0*t(:,iorb) + cc(2)*t1(:,iorb)
  end do

   
  call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp,norb,0.d0,t1,norb)
  
  call sparsemm(ovrlp_tmp,t1_tmp,ts,norb,tmb%mad)

  do i=1,norb
    do j=1,norb
    if(iproc==0) write(100,*) i,j,t1(i,j)
    if(iproc==0) write(200,*) i,j,ts(i,j)
    end do
  end do

 

  do ipl=3,npl
     !calculate (3/2 - 1/2 S) H (3/2 - 1/2 S) t
     call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp,norb,0.d0,t1,norb)
     call dsymm('L', 'U', norb, norbp,1.d0,ham,norb,t1,norb,0.d0,t1_tmp,norb)
     call dsymm('L', 'U', norb, norbp,1.d0,ovrlp_tmp,norb,t1_tmp,norb,0.d0,t1,norb)
     !calculate t2 = 2 * (3/2 - 1/2 S) H (3/2 - 1/2 S) t1 - t
     t2 = 2*t1 - t
     !update fermi-matrix
     fermi(:,isorb+1:isorb+norbp+1) = cc(ipl)*t2   
     !update t's
     t = t1
     t1_tmp = t2
 end do


  iall=-product(shape(column))*kind(column)
  deallocate(column, stat=istat)
  call memocc(istat, iall, 'column', subname)
  iall=-product(shape(t))*kind(t)
  deallocate(t, stat=istat)
  call memocc(istat, iall, 't', subname)
  iall=-product(shape(t1))*kind(t1)
  deallocate(t1, stat=istat)
  call memocc(istat, iall, 't1', subname)
  iall=-product(shape(t1_tmp))*kind(t1_tmp)
  deallocate(t1_tmp, stat=istat)
  call memocc(istat, iall, 't1_tmp', subname)
  iall=-product(shape(t2))*kind(t2)
  deallocate(t2, stat=istat)
  call memocc(istat, iall, 't2', subname)
  iall=-product(shape(ts))*kind(ts)
  deallocate(ts, stat=istat)
  call memocc(istat, iall, 'ts', subname)

end subroutine chebyshev

subroutine sparsemm(norb,a,b,c,mad)

use module_base
use module_types


  implicit none

  !Calling Arguments
  integer, intent(in) :: norb
  real(kind=8), dimension(norb,norb),intent(in) :: a,b
  real(kind=8), dimension(norb,norb),intent(out) :: c
  type(matrixDescriptors),intent(in) :: mad

  !Local variables
  integer :: i,j,iseg,jorb,iiorb,jjorb
  real(kind=8) :: temp


  call to_zero(norb*norb,c(1,1))

  do i = 1,norb
    do j = 1,norb

        temp = b(i,j)
  
        do iseg = 1,mad%nseg
          do jorb = mad%keyg(1,iseg), mad%keyg(2,iseg)
            iiorb = (jorb-1)/norb + 1
            jjorb = jorb - (iiorb - 1)*norb
            c(iiorb,j) = c(iiorb,j) + temp*a(j,iiorb) 
          end do
        end do
     end do
  end do 
  
    
end subroutine sparsemm

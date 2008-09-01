! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points phidst and the associated error 
! vectors (preconditione gradients) hpsidst
subroutine diisstp(norb,norbp,nproc,iproc,nspinor,  & 
                   ads,ids,mids,idsx,nvctrp,psit,psidst,hpsidst)
  use module_base
  implicit none
! Arguments
  integer, intent(in) :: norb,norbp,nproc,iproc,nspinor,ids,mids,idsx,nvctrp
  real(wp), dimension(nvctrp,norbp*nproc*nspinor,idsx), intent(in) :: psidst,hpsidst
  real(dp), dimension(idsx+1,idsx+1,3), intent(inout) :: ads
  real(wp), dimension(nvctrp,norbp*nproc*nspinor), intent(out) :: psit
! Local variables
  character(len=*), parameter :: subname='diisstp'
  integer :: i,j,ist,jst,mi,iorb,info,jj,mj,k,i_all,i_stat,ierr
  real(kind=8) :: tt
  real(kind=8), external :: DDOT
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:), allocatable :: rds

  allocate(ipiv(idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,ipiv,'ipiv',subname)
  allocate(rds(idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,rds,'rds',subname)

! set up DIIS matrix (upper triangle)
  if (ids > idsx) then
! shift left up matrix
     do i=1,idsx-1
        do j=1,i
           ads(j,i,1)=ads(j+1,i+1,1)
        end do
     end do
  end if

! calculate new line, use rds as work array for summation
  call razero(idsx+1,rds)
  ist=max(1,ids-idsx+1)
  do i=ist,ids
     mi=mod(i-1,idsx)+1
     do iorb=1,norb*nspinor
        tt=DDOT(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
        rds(i-ist+1)=rds(i-ist+1)+tt
     end do
  end do

  if (nproc > 1) then
     call MPI_ALLREDUCE(rds,ads(1,min(idsx,ids),1),min(ids,idsx),  & 
                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do i=1,min(ids,idsx)
        ads(i,min(idsx,ids),1)=rds(i)
    end do
  endif

! copy to work array, right hand side, boundary elements
  do j=1,min(idsx,ids)
     ads(j,min(idsx,ids)+1,2)=1.d0
     rds(j)=0.d0
     do i=j,min(idsx,ids)
        ads(j,i,2)=ads(j,i,1)
     end do
  end do
  ads(min(idsx,ids)+1,min(idsx,ids)+1,2)=0.0_dp
  rds(min(idsx,ids)+1)=1.0_dp

  !if(iproc==0)  write(6,*) 'DIIS matrix'
  !do i=1,min(idsx,ids)+1
  !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(idsx,ids)+1),rds(i)
  !enddo
  if (ids.gt.1) then
     ! solve linear system:(LAPACK)
     call DSYSV('U',min(idsx,ids)+1,1,ads(1,1,2),idsx+1,  & 
          ipiv,rds,idsx+1,ads(1,1,3),(idsx+1)**2,info)
     
     if (info.ne.0) then
        print*, 'diisstp: DSYSV',info
        stop
     end if
  else
     rds(1)=1.0_dp
  endif
  if (iproc.eq.0) then 
     !write(*,*) 'DIIS weights'
     write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j),j=1,min(idsx,ids)+1)
  endif

! new guess
  do iorb=1,norb*nspinor
     call razero(nvctrp,psit(1,iorb))

     jst=max(1,ids-idsx+1)
     jj=0
     do j=jst,ids
        jj=jj+1
        mj=mod(j-1,idsx)+1
        do k=1,nvctrp
           psit(k,iorb)=psit(k,iorb)+rds(jj)*(psidst(k,iorb,mj)-hpsidst(k,iorb,mj))
        end do
     end do
  end do

  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv',subname)
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds',subname)

END SUBROUTINE diisstp

!!****f* BigDFT/psimix
!! FUNCTION
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine psimix(iproc,nproc,orbs,comms,ads,ids,mids,idsx,energy,energy_old,alpha,&
     hpsit,psidst,hpsidst,psit)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,ids,mids,idsx
  real(gp), intent(in) :: energy,energy_old
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  real(gp), intent(inout) :: alpha
  real(wp), dimension(:), pointer :: psit,hpsit,psidst,hpsidst
  real(wp), dimension(:,:,:), pointer :: ads
  !local variables
  integer :: ikptp,nvctrp,ispsi,ispsidst

  if (idsx > 0) then
     !do not transpose the hpsi wavefunction into the diis array
     !for compatibility with the k-points distribution
     ispsi=1
     ispsidst=1
     do ikptp=1,orbs%nkptsp
        nvctrp=comms%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle
        
     !here we can choose to store the DIIS arrays with single precision
     !psidst=psit
        call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
             psit(ispsi),1,&
             psidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(mids-1)),1)
     !hpsidst=hpsi
        call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
             hpsit(ispsi),1,&
             hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(mids-1)),1)
        ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
        ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*idsx
     end do
     
     !here we should separate between up and down spin orbitals, maybe
     call diisstp(iproc,nproc,orbs,comms,ads,ids,mids,idsx,&
          psit,psidst,hpsidst)

  else
     ! update all wavefunctions with the preconditioned gradient
     if (energy > energy_old) then
        alpha=max(.125_wp,.5_wp*alpha)
        if (alpha == .125_wp) write(*,*) ' WARNING: Convergence problem or limit'
     else
        alpha=min(1.05_wp*alpha,1._wp)
     endif
     if (iproc == 0 .and. verbose > 0) write(*,'(1x,a,1pe11.3)') 'alpha=',alpha

!!     do iorb=1,orbs%norb*orbs%nspinor
!!        call axpy(comms%nvctr_par(iproc),&
!!             -alpha,hpsi(1+comms%nvctr_par(iproc)*(iorb-1)),1,&
!!             psit(1+comms%nvctr_par(iproc)*(iorb-1)),1)
!!     enddo

     call axpy(sum(comms%ncntt(0:nproc-1)),-alpha,hpsit(1),1,psit(1),1)

  endif

END SUBROUTINE psimix
!!***


! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points psidst and the associated error 
! vectors (preconditioned gradients) hpsidst
subroutine diisstp(iproc,nproc,orbs,comms,ads,ids,mids,idsx,psit,psidst,hpsidst)
  use module_base
  use module_types
  implicit none
! Arguments
  integer, intent(in) :: nproc,iproc,ids,mids,idsx
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))*idsx), intent(in) :: psidst,hpsidst
  real(dp), dimension(idsx+1,idsx+1,orbs%nkptsp,3), intent(inout) :: ads
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
! Local variables
  character(len=*), parameter :: subname='diisstp'
  integer :: i,j,ist,jst,mi,iorb,info,jj,mj,k,i_all,i_stat,ierr
  integer :: ikptp,ikpt,ispsi,ispsidst,nvctrp
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:,:), allocatable :: rds

  allocate(ipiv(idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,ipiv,'ipiv',subname)
  allocate(rds(idsx+1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,rds,'rds',subname)

  call razero((idsx+1)*orbs%nkpts,rds)

  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     ! set up DIIS matrix (upper triangle)
     if (ids > idsx) then
        ! shift left up matrix
        do i=1,idsx-1
           do j=1,i
              ads(j,i,ikptp,1)=ads(j+1,i+1,ikptp,1)
           end do
        end do
     end if

     ! calculate new line, use rds as work array for summation
     ist=max(1,ids-idsx+1)
     do i=ist,ids
        mi=mod(i-1,idsx)+1
!!     do iorb=1,norb*nspinor
!!        tt=dot(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
!!        rds(i-ist+1)=rds(i-ist+1)+tt
!!     end do
        !to be corrected for complex wavefunctions
        rds(i-ist+1,ikpt)=dot(nvctrp*orbs%norb*orbs%nspinor,&
             hpsidst(ispsidst+(mids-1)*nvctrp*orbs%norb*orbs%nspinor),1,&
             hpsidst(ispsidst+(mi-1)*nvctrp*orbs%norb*orbs%nspinor),1)
     end do

     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*idsx
  end do

  if (nproc > 1) then

     call MPI_ALLREDUCE(MPI_IN_PLACE,rds,(idsx+1)*orbs%nkpts,  & 
                 mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)

!!     call MPI_ALLREDUCE(rds,ads(1,min(idsx,ids),1),min(ids,idsx),  & 
!!                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  endif

  
  ispsi=1
  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     do i=1,min(ids,idsx)
        ads(i,min(idsx,ids),ikptp,1)=rds(i,ikpt)
     end do

     ! copy to work array, right hand side, boundary elements
     do j=1,min(idsx,ids)
        ads(j,min(idsx,ids)+1,ikptp,2)=1.d0
        rds(j,ikpt)=0.d0
        do i=j,min(idsx,ids)
           ads(j,i,ikptp,2)=ads(j,i,ikptp,1)
        end do
     end do
     ads(min(idsx,ids)+1,min(idsx,ids)+1,ikptp,2)=0.0_dp
     rds(min(idsx,ids)+1,ikpt)=1.0_dp
     
     !if(iproc==0)  write(6,*) 'DIIS matrix'
     !do i=1,min(idsx,ids)+1
     !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(idsx,ids)+1),rds(i)
     !enddo
     if (ids > 1) then
        ! solve linear system:(LAPACK)
        call DSYSV('U',min(idsx,ids)+1,1,ads(1,1,ikptp,2),idsx+1,  & 
             ipiv,rds(1,ikpt),idsx+1,ads(1,1,ikptp,3),(idsx+1)**2,info)
        
        if (info /= 0) then
           print*, 'diisstp: DSYSV',info
        end if
     else
        rds(1,ikpt)=1.0_dp
     endif
     if (iproc == 0 .and. verbose > 0) then 
        !write(*,*) 'DIIS weights'
        !we should print the weights for each k-point
        write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j,ikpt),j=1,min(idsx,ids)+1)
     endif

! new guess
     do iorb=1,orbs%norb
        call razero(nvctrp*orbs%nspinor,psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor))
        
        jst=max(1,ids-idsx+1)
        jj=0
        do j=jst,ids
           jj=jj+1
           mj=mod(j-1,idsx)+1
           do k=1,nvctrp*orbs%nspinor
              psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)=&
                   psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)+&
                   rds(jj,ikpt)*(&
                   psidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp)-&
                   hpsidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp))
           end do
        end do
     end do
     ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*idsx
  end do
  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv',subname)
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds',subname)

END SUBROUTINE diisstp

subroutine mix_rhopot(npoints,alphamix,rhopot_old,rhopot,rpnrm)
  use module_base
  implicit none
  integer, intent(in) :: npoints
  real(gp), intent(in) :: alphamix
  real(dp), dimension(npoints), intent(inout) :: rhopot,rhopot_old
  real(gp), intent(out) :: rpnrm
  !local variables
  integer :: ierr
  
  !vold=>vold-vnew
  call axpy(npoints,-1.0_dp,rhopot(1),1,rhopot_old(1),1)

  !calculate rhopot_norm
  rpnrm=(nrm2(npoints,rhopot_old(1),1))**2
  if (npoints > 0) rpnrm=rpnrm/real(npoints,gp)
  call mpiallred(rpnrm,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  rpnrm=sqrt(rpnrm)
      
  !vnew=vnew+alpha(vold-vnew)
  call axpy(npoints,alphamix,rhopot_old(1),1,rhopot(1),1)
  
  !vold=vnew
  call dcopy(npoints,rhopot(1),1,rhopot_old(1),1)

end subroutine mix_rhopot


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
subroutine psimix(iproc,nproc,orbs,comms,diis,hpsit,psit)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
  !real(wp), dimension(:), pointer :: psit,hpsit
  !local variables
  integer :: ikptp,nvctrp,ispsi,ispsidst,jj,ierr

  if (diis%idsx > 0) then
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
             diis%psidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)),1)
     !hpsidst=hpsi
     !   call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
     !        hpsit(ispsi),1,&
     !        hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(mids-1)),1)

     do jj=0,nvctrp*orbs%norb*orbs%nspinor-1
        diis%hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)+jj)&
             =real(hpsit(ispsi+jj),tp) !diis precision conversion
     end do
        ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
        ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
     end do

     !here we should separate between up and down spin orbitals, maybe
     call diisstp(iproc,nproc,orbs,comms,diis,psit)

  else
     ! update all wavefunctions with the preconditioned gradient
     if (diis%energy > diis%energy_old) then
        diis%alpha=max(.125_wp,.5_wp*diis%alpha)
        if (diis%alpha == .125_wp) write(*,*) ' WARNING: Convergence problem or limit'
     else
        diis%alpha=min(1.05_wp*diis%alpha,1._wp)
     endif
     if (iproc == 0 .and. verbose > 0) write(*,'(1x,a,1pe11.3)') 'alpha=',diis%alpha

!!     do iorb=1,orbs%norb*orbs%nspinor
!!        call axpy(comms%nvctr_par(iproc),&
!!             -alpha,hpsi(1+comms%nvctr_par(iproc)*(iorb-1)),1,&
!!             psit(1+comms%nvctr_par(iproc)*(iorb-1)),1)
!!     enddo

     call axpy(sum(comms%ncntt(0:nproc-1)),-diis%alpha,hpsit(1),1,psit(1),1)

  endif

END SUBROUTINE psimix
!!***


subroutine diis_or_sd(iproc,idsx,nkptsp,diis)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,idsx,nkptsp
  type(diis_objects), intent(inout) :: diis

  !here we should add the DIIS/SD switch
  !add switch between DIIS and SD
  if (diis%energy == diis%energy_min .and. .not. diis%switchSD) diis%idiistol=0
  if (diis%energy > diis%energy_min .and. idsx >0 .and. .not. diis%switchSD) then
     diis%idiistol=diis%idiistol+1
  end if
  if (diis%idiistol > idsx .and. .not. diis%switchSD) then
     !the energy has not decreasing for too much steps, switching to SD for next steps
     if (iproc ==0) write(*,'(1x,a,1pe9.2,a)')&
          'WARNING: The energy value is growing (delta=',diis%energy-diis%energy_min,') switch to SD'
     diis%switchSD=.true.
     !il is not necessary to deallocate these arrays
     !i_all=-product(shape(psidst))*kind(psidst)
     !deallocate(psidst,stat=i_stat)
     !call memocc(i_stat,i_all,'psidst',subname)
     !i_all=-product(shape(hpsidst_sp))*kind(hpsidst_sp)
     !deallocate(hpsidst_sp,stat=i_stat)
     !call memocc(i_stat,i_all,'hpsidst_sp',subname)
     !i_all=-product(shape(ads))*kind(ads)
     !deallocate(ads,stat=i_stat)
     !call memocc(i_stat,i_all,'ads',subname)
     diis%idsx=0
     diis%idiistol=0
  end if

  if ((diis%energy == diis%energy_min) .and. diis%switchSD) then
     diis%idiistol=diis%idiistol+1
  end if
  if (diis%idiistol > idsx .and. diis%switchSD) then
     !if (diis%idiistol > 100*idsx .and. diis%switchSD) then
     !restore the original DIIS
     if (iproc ==0) write(*,'(1x,a,1pe9.2)')&
          'WARNING: The energy value is now decreasing again, coming back to DIIS'
     diis%switchSD=.false.
     diis%idsx=idsx
     diis%ids=0
     diis%idiistol=0

     !no need to reallocate
     !allocate(psidst(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     !call memocc(i_stat,psidst,'psidst',subname)
     !allocate(hpsidst_sp(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     !call memocc(i_stat,hpsidst_sp,'hpsidst_sp',subname)
     !allocate(ads(idsx+1,idsx+1,orbs%nkptsp*3+ndebug),stat=i_stat)
     !call memocc(i_stat,ads,'ads',subname)

     call razero(nkptsp*3*(idsx+1)**2,diis%ads)
  end if

end subroutine diis_or_sd


! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points psidst and the associated error 
! vectors (preconditioned gradients) hpsidst
subroutine diisstp(iproc,nproc,orbs,comms,diis,psit)
  use module_base
  use module_types
  implicit none
! Arguments
  integer, intent(in) :: nproc,iproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
! Local variables
  character(len=*), parameter :: subname='diisstp'
  integer :: i,j,ist,jst,mi,iorb,info,jj,mj,k,i_all,i_stat,ierr
  integer :: ikptp,ikpt,ispsi,ispsidst,nvctrp
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:,:), allocatable :: rds

  allocate(ipiv(diis%idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,ipiv,'ipiv',subname)
  allocate(rds(diis%idsx+1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,rds,'rds',subname)

  call razero((diis%idsx+1)*orbs%nkpts,rds)

  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     ! set up DIIS matrix (upper triangle)
     if (diis%ids > diis%idsx) then
        ! shift left up matrix
        do i=1,diis%idsx-1
           do j=1,i
              diis%ads(j,i,ikptp,1)=diis%ads(j+1,i+1,ikptp,1)
           end do
        end do
     end if

     ! calculate new line, use rds as work array for summation
     ist=max(1,diis%ids-diis%idsx+1)
     do i=ist,diis%ids
        mi=mod(i-1,diis%idsx)+1
!!     do iorb=1,norb*nspinor
!!        tt=dot(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
!!        rds(i-ist+1)=rds(i-ist+1)+tt
!!     end do
        !to be corrected for complex wavefunctions
        rds(i-ist+1,ikpt)=dot(nvctrp*orbs%norb*orbs%nspinor,&
             diis%hpsidst(ispsidst+(diis%mids-1)*nvctrp*orbs%norb*orbs%nspinor),1,&
             diis%hpsidst(ispsidst+(mi-1)*nvctrp*orbs%norb*orbs%nspinor),1)
        !this has to be inserted in module_base
        !call ds_dot(nvctrp*orbs%norb*orbs%nspinor,&
        !     hpsidst_sp,ispsidst+(mids-1)*nvctrp*orbs%norb*orbs%nspinor,1,&
        !     hpsidst_sp,ispsidst+(mi  -1)*nvctrp*orbs%norb*orbs%nspinor,1,&
        !     rds(i-ist+1,ikpt))
     end do
     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
  end do

  if (nproc > 1) then
     call mpiallred(rds(1,1),(diis%idsx+1)*orbs%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call MPI_ALLREDUCE(MPI_IN_PLACE,rds,(diis%idsx+1)*orbs%nkpts,  & 
!                 mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!!     call MPI_ALLREDUCE(rds,ads(1,min(diis%idsx,ids),1),min(ids,diis%idsx),  & 
!!                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  
  ispsi=1
  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     do i=1,min(diis%ids,diis%idsx)
        diis%ads(i,min(diis%idsx,diis%ids),ikptp,1)=rds(i,ikpt)
     end do

     ! copy to work array, right hand side, boundary elements
     do j=1,min(diis%idsx,diis%ids)
        diis%ads(j,min(diis%idsx,diis%ids)+1,ikptp,2)=1.0_wp
        rds(j,ikpt)=0.d0
        do i=j,min(diis%idsx,diis%ids)
           diis%ads(j,i,ikptp,2)=diis%ads(j,i,ikptp,1)
        end do
     end do
     diis%ads(min(diis%idsx,diis%ids)+1,min(diis%idsx,diis%ids)+1,ikptp,2)=0.0_dp
     rds(min(diis%idsx,diis%ids)+1,ikpt)=1.0_dp
     
     !if(iproc==0)  write(6,*) 'DIIS matrix'
     !do i=1,min(diis%idsx,ids)+1
     !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(diis%idsx,ids)+1),rds(i)
     !enddo
     if (diis%ids > 1) then
        ! solve linear system:(LAPACK)
        call DSYSV('U',min(diis%idsx,diis%ids)+1,1,diis%ads(1,1,ikptp,2),diis%idsx+1,  & 
             ipiv,rds(1,ikpt),diis%idsx+1,diis%ads(1,1,ikptp,3),(diis%idsx+1)**2,info)
        
        if (info /= 0) then
           print*, 'diisstp: DSYSV',info
        end if
     else
        rds(1,ikpt)=1.0_dp
     endif

! new guess
     do iorb=1,orbs%norb
        call razero(nvctrp*orbs%nspinor,psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor))
        
        jst=max(1,diis%ids-diis%idsx+1)
        jj=0
        do j=jst,diis%ids
           jj=jj+1
           mj=mod(j-1,diis%idsx)+1
           do k=1,nvctrp*orbs%nspinor
              psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)=&
                   psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)+&
                   rds(jj,ikpt)*(&
                   diis%psidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp)&
                   -real(diis%hpsidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp),wp))
           end do
        end do
     end do
     ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
  end do
  ! Output to screen, depending on policy.
  if (verbose >= 10) then
     call broadcast_kpt_objects(nproc, orbs%nkpts, diis%idsx+1, rds, orbs%ikptproc)
  end if
  if (iproc == 0 .and. verbose > 0) then 
     if (verbose < 10) then
        !we restrict the printing to the first k point only.
        write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j,1),j=1,min(diis%idsx,diis%ids)+1)
     else
        do ikpt = 1, orbs%nkpts
           write(*,'(1x,a,I3.3,a,2x,12(1x,1pe9.2))')'DIIS weights (kpt #', ikpt, &
                & ')', (rds(j,0),j=1,min(diis%idsx,diis%ids)+1)
        end do
     end if
  endif
  
  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv',subname)
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds',subname)

END SUBROUTINE diisstp

! compute a dot product of two single precision vectors 
! returning a double precision result
subroutine ds_dot(ndim,x,x0,dx,y,y0,dy,dot_out)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ndim,x0,dx,y0,dy
  real(tp), dimension(ndim), intent(in) :: x,y
  real(wp), intent(out) :: dot_out
  integer :: jj,ix,iy
  
  dot_out=0.0d0
  ix=x0
  iy=y0
  do jj=1,ndim
    dot_out=dot_out + real(x(ix),wp)*real(y(iy),wp)
    ix=ix+dx
    iy=iy+dy
  end do
end subroutine ds_dot

function s2d_dot(ndim,x,dx,y,dy)
  implicit none
  integer, intent(in) :: ndim,dx,dy
  real(kind=4), dimension(ndim), intent(in) :: x,y
  real(kind=8) :: s2d_dot
  !local variables
  integer :: jj,ix,iy
  real(kind=8) :: dot_out

  !quick return if possible
  if (ndim <= 0) then
     s2d_dot=0.0d0
     return
  end if

  dot_out=0.0d0
  ix=1
  iy=1
  do jj=1,ndim
     dot_out=dot_out + real(x(ix),kind=8)*real(y(iy),kind=8)
     ix=ix+dx
     iy=iy+dy
  end do

  s2d_dot=dot_out

end function s2d_dot


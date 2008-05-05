subroutine transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  integer recvcount,sendcount
  dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
  character(len=*), parameter :: subname='transallwwaves'
  real(kind=8), allocatable :: psiw(:,:,:)
  include 'mpif.h'

  call timing(iproc,'Un-Transall   ','ON')
  allocate(psiw(nvctrp,norbp,nproc+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  sendcount=nvctrp*norbp
  recvcount=nvctrp*norbp

  !  reformatting: psiw(i,iorb,j,jorb) <- psi(ij,iorb,jorb)
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     ij=1
     do j=1,nproc
        do i=1,nvctrp
           if (ij .le. nvctr_c+7*nvctr_f) then
              psiw(i,iorb-iproc*norbp,j)=psi(ij,iorb-iproc*norbp)
           else
              psiw(i,iorb-iproc*norbp,j)=0.d0
           endif
           ij=ij+1
        enddo
     enddo
  enddo

  !  transposition: psit(i,iorb,jorb,j) <- psiw(i,iorb,j,jorb) 
  call MPI_ALLTOALL(psiw,sendcount,MPI_DOUBLE_PRECISION,  &
       psit,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  call timing(iproc,'Un-Transall   ','OF')

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)

END SUBROUTINE transallwaves

subroutine switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psiw,nspinor)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspinor
  real(kind=8), dimension(nvctr_c+7*nvctr_f,nspinor,norbp), intent(in) :: psi
  real(kind=8), dimension(nspinor*nvctrp,norbp,nproc), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij,isp

  if(nspinor==1) then
     isp=1
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        ij=1
        do j=1,nproc
           do i=1,nvctrp
              if (ij <= nvctr_c+7*nvctr_f) then
                 psiw(i,iorb-iproc*norbp,j)=psi(ij,nspinor,iorb-iproc*norbp)
              else
                 psiw(i,iorb-iproc*norbp,j)=0.d0
              endif
              ij=ij+1
           enddo
        enddo
     enddo
  else
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        ij=1
        do j=1,nproc
           do i=1,nvctrp
              if (ij <= nvctr_c+7*nvctr_f) then
                 psiw(2*i-1,iorb-iproc*norbp,j)=psi(ij,1,iorb-iproc*norbp)
                 psiw(2*i,iorb-iproc*norbp,j)=psi(ij,2,iorb-iproc*norbp)
                 psiw(2*i+2*nvctrp-1,iorb-iproc*norbp,j)=psi(ij,3,iorb-iproc*norbp)
                 psiw(2*i+2*nvctrp,iorb-iproc*norbp,j)=psi(ij,4,iorb-iproc*norbp)
              else
                 psiw(2*i-1,iorb-iproc*norbp,j)=0.0d0
                 psiw(2*i,iorb-iproc*norbp,j)=0.0d0
                 psiw(2*i+2*nvctrp-1,iorb-iproc*norbp,j)=0.0d0
                 psiw(2*i+2*nvctrp,iorb-iproc*norbp,j)=0.0d0
              endif
              ij=ij+1
           enddo
        enddo
     enddo
  end if
end subroutine switch_waves

subroutine unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,psi,nspinor)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspinor
  real(kind=8), dimension(nspinor*nvctrp,norbp,nproc), intent(in) :: psiw
  real(kind=8), dimension(nvctr_c+7*nvctr_f,nspinor,norbp), intent(out) :: psi
  !local variables
  integer :: iorb,i,j,ij,isp

  if(nspinor==1) then
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        ij=1
        loop: do j=1,nproc
           do i=1,nvctrp
              psi(ij,nspinor,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
              ij=ij+1
              if (ij > nvctr_c+7*nvctr_f) exit loop
           enddo
        enddo loop
     enddo
  else
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
        ij=1
        loop4: do j=1,nproc
           do i=1,nvctrp
              psi(ij,1,iorb-iproc*norbp)=psiw(2*i-1,iorb-iproc*norbp,j)
              psi(ij,2,iorb-iproc*norbp)=psiw(2*i,iorb-iproc*norbp,j)
              psi(ij,3,iorb-iproc*norbp)=psiw(2*i+2*nvctrp-1,iorb-iproc*norbp,j)
              psi(ij,4,iorb-iproc*norbp)=psiw(2*i+2*nvctrp,iorb-iproc*norbp,j)
              ij=ij+1
              if (ij >  nvctr_c+7*nvctr_f) exit loop4
           enddo
        enddo loop4
     enddo
  end if
  
end subroutine unswitch_waves

subroutine untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  integer recvcount,sendcount
  dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
  character(len=*), parameter :: subname='untransallwaves'
  real(kind=8), allocatable :: psiw(:,:,:)
  include 'mpif.h'

  call timing(iproc,'Un-Transall   ','ON')

  allocate(psiw(nvctrp,norbp,nproc+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  sendcount=nvctrp*norbp
  recvcount=nvctrp*norbp

  ! transposition: psiw(i,iorb,j,jorb) <- psit(i,iorb,jorb,j) 
  call MPI_ALLTOALL(psit,sendcount,MPI_DOUBLE_PRECISION,  &
       psiw,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

  ! reformatting: psi(ij,iorb,jorb) <- psiw(i,iorb,j,jorb)
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     ij=1
     do j=1,nproc
        do i=1,nvctrp
           psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
           ij=ij+1
           if (ij.gt. nvctr_c+7*nvctr_f) goto 333
        enddo
     enddo
333  continue
  enddo

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)
  call timing(iproc,'Un-Transall   ','OF')


END SUBROUTINE untransallwaves

subroutine psitransspi(nvctrp,norb,psi,forward)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp
  logical, intent(in) :: forward
  real(kind=8), dimension(4*nvctrp,norb), intent(inout) :: psi
  real(kind=8), dimension(:,:,:), allocatable :: tpsit
  !local variables
  character(len=*), parameter :: subname='psitransspi'
  integer :: i,iorb,ij,isp,i_all,i_stat

!  call timing(0,'Un-Transall   ','ON')

  allocate(tpsit(nvctrp,4,norb+ndebug),stat=i_stat)
  call memocc(i_stat,tpsit,'tpsit',subname)
!  print '(a,5f10.5,2i5,l1)','p2s',(abs(psi(i,1)),i=1,5),shape(psi),forward
  if(forward) then
!     print *, 'transposing psi forward'
!      print '(a,30f10.5)','tsp',((sum(abs(psi(:,iorb)))),iorb=1,norb)
    do iorb=1,norb
!        ij=1
        do isp=1,4
           do i=1,nvctrp
              tpsit(i,isp,iorb)=psi(i+(isp-1)*nvctrp,iorb)
!              ij=ij+1
           enddo
        enddo
     enddo
!      print '(a,30f10.5)','tsp',((sum(abs(tpsit(:,isp,iorb))),isp=1,4),iorb=1,norb)
    do iorb=1,norb
!        ij=1
        do i=1,nvctrp
           psi(2*i-1,iorb)=tpsit(i,1,iorb)
           psi(2*i,iorb)=tpsit(i,2,iorb)
           psi(2*i+2*nvctrp-1,iorb)=tpsit(i,3,iorb)
           psi(2*i+2*nvctrp,iorb)=tpsit(i,4,iorb)
!           ij=ij+1
        enddo
     enddo
  else
!     print *,'tsp',((sum(abs(psi(:,iorb)))),iorb=1,norb)
!     print *,'transposing psi backward'
     do iorb=1,norb
!        ij=1
        do i=1,nvctrp
           tpsit(i,1,iorb)=psi(2*i-1,iorb)
           tpsit(i,2,iorb)=psi(2*i,iorb)
           tpsit(i,3,iorb)=psi(2*i-1+2*nvctrp,iorb)
           tpsit(i,4,iorb)=psi(2*i+2*nvctrp,iorb)
!           ij=ij+1
        enddo
     enddo
!     print *,'tsp',((sum(abs(tpsit(:,isp,iorb))),isp=1,4),iorb=1,norb)
     do iorb=1,norb
!        ij=1
        do isp=1,4
           do i=1,nvctrp
              psi(i+(isp-1)*nvctrp,iorb)=tpsit(i,isp,iorb)
              ij=ij+1
           enddo
        enddo
     enddo
  end if

  i_all=-product(shape(tpsit))*kind(tpsit)
  deallocate(tpsit,stat=i_stat)
  call memocc(i_stat,i_all,'tpsit',subname)

!   call timing(0,'Un-Transall   ','OF')
!  print '(a,5f10.5,2i5,l1)','p2E',(abs(psi(i,1)),i=1,5),shape(psi),forward


end subroutine psitransspi

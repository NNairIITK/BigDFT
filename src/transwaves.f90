subroutine transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
  implicit real*8 (a-h,o-z)
  logical, parameter :: parallel=.true.
  integer recvcount,sendcount
  dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
  real*8, allocatable :: psiw(:,:,:)
  include 'mpif.h'

  call timing(iproc,'Un-Transall   ','ON')
  allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
  call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','transallwaves')

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
  call memocc(i_stat,i_all,'psiw','transallwaves')

END SUBROUTINE transallwaves

subroutine switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psiw)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(in) :: psi
  real(kind=8), dimension(nvctrp,norbp,nproc), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij

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

end subroutine switch_waves

subroutine unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,psi)
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp
  real(kind=8), dimension(nvctrp,norbp,nproc), intent(in) :: psiw
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
  !local variables
  integer :: iorb,i,j,ij

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     ij=1
     loop: do j=1,nproc
        do i=1,nvctrp
           psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
           ij=ij+1
           if (ij.gt. nvctr_c+7*nvctr_f) exit loop
        enddo
     enddo loop
  enddo

end subroutine unswitch_waves

subroutine untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
  implicit real*8 (a-h,o-z)
  logical, parameter :: parallel=.true.
  integer recvcount,sendcount
  dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
  real*8, allocatable :: psiw(:,:,:)
  include 'mpif.h'

  call timing(iproc,'Un-Transall   ','ON')

  allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
  call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','untransallwaves')

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
  call memocc(i_stat,i_all,'psiw','untransallwaves')
  call timing(iproc,'Un-Transall   ','OF')


END SUBROUTINE untransallwaves

subroutine transpose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
     work,out) !optional
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norbp), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work,out
  !real(wp), dimension(*), intent(out), optional :: out
  !local variables
  include 'mpif.h'
  integer :: mpidatatype,ierr

  call timing(iproc,'Un-TransSwitch','ON')
  if (wp == kind(1.d0)) then
     mpidatatype=MPI_DOUBLE_PRECISION
  else
     mpidatatype=MPI_REAL
  end if

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work) .or. &
          product(shape(work)) < nspinor*nvctrp*norbp*nproc) then
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for transposing in parallel"
        stop
     end if
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,work,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     if (present(out)) then
        call MPI_ALLTOALL(work,nvctrp*norbp*nspinor,mpidatatype,  &
             out,nvctrp*norbp*nspinor,mpidatatype,MPI_COMM_WORLD,ierr)
     else
        call MPI_ALLTOALL(work,nvctrp*norbp*nspinor,mpidatatype,  &
             psi,nvctrp*norbp*nspinor,mpidatatype,MPI_COMM_WORLD,ierr)
     end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     if(nspinor==4) then
        call psitransspi(nvctrp,norb,psi,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

end subroutine transpose

subroutine untranspose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
     work,out) !optional
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
  real(wp), dimension(nspinor*nvctrp,norbp,nproc), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work,out
  !real(wp), dimension(*), intent(out), optional :: out
  !local variables
  include 'mpif.h'
  integer :: mpidatatype,ierr

  call timing(iproc,'Un-TransSwitch','ON')

  if (wp == kind(1.d0)) then
     mpidatatype=MPI_DOUBLE_PRECISION
  else
     mpidatatype=MPI_REAL
  end if

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work) .or. &
          product(shape(work)) < nspinor*nvctrp*norbp*nproc) then
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for untransposing in parallel"
        stop
     end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,mpidatatype,  &
          work,nvctrp*norbp*nspinor,mpidatatype,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     if (present(out)) then
        call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,&
        work,out,nspinor)
     else
        call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,&
        work,psi,nspinor)
     end if
  else
     if(nspinor==4) then
        call psitransspi(nvctrp,norb,psi,.false.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
end subroutine untranspose

subroutine switch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psiw,nspinor)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspinor
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor,norbp), intent(in) :: psi
  real(wp), dimension(nspinor*nvctrp,norbp,nproc), intent(out) :: psiw
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
                 psiw(i,iorb-iproc*norbp,j)=0.0_wp
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
                 psiw(2*i-1,iorb-iproc*norbp,j)=0.0_wp
                 psiw(2*i,iorb-iproc*norbp,j)=0.0_wp
                 psiw(2*i+2*nvctrp-1,iorb-iproc*norbp,j)=0.0_wp
                 psiw(2*i+2*nvctrp,iorb-iproc*norbp,j)=0.0_wp
              endif
              ij=ij+1
           enddo
        enddo
     enddo
  end if
end subroutine switch_waves

subroutine unswitch_waves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psiw,psi,nspinor)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspinor
  real(wp), dimension(nspinor*nvctrp,norbp,nproc), intent(in) :: psiw
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor,norbp), intent(out) :: psi
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

subroutine psitransspi(nvctrp,norb,psi,forward)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp
  logical, intent(in) :: forward
  real(wp), dimension(4*nvctrp,norb), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='psitransspi'
  integer :: i,iorb,ij,isp,i_all,i_stat
  real(wp), dimension(:,:,:), allocatable :: tpsit

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

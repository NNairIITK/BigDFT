!transpose the routine to an array. Does not need interface
subroutine transposeto(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,work,out) 
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norbp), intent(inout) :: psi
  real(wp), dimension(*), intent(inout) :: work
  real(wp), dimension(*), intent(out) :: out
  !local variables
  integer :: ierr

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,work,nspinor)
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(work,nvctrp*norbp*nspinor,mpidtypw,  &
          out,nvctrp*norbp*nspinor,mpidtypw,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     if(nspinor==4) then
        call psitransspi(nspinor,nvctrp,norb,psi,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

end subroutine transposeto


subroutine transpose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norbp), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr

  call timing(iproc,'Un-TransSwitch','ON')

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
     if (present(outadd)) then
        call MPI_ALLTOALL(work,nvctrp*norbp*nspinor,mpidtypw,  &
             outadd,nvctrp*norbp*nspinor,mpidtypw,MPI_COMM_WORLD,ierr)
     else
        call MPI_ALLTOALL(work,nvctrp*norbp*nspinor,mpidtypw,  &
             psi,nvctrp*norbp*nspinor,mpidtypw,MPI_COMM_WORLD,ierr)
     end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     if(nspinor==4) then
        call psitransspi(nspinor,nvctrp,norb,psi,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

end subroutine transpose

subroutine untranspose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
  real(wp), dimension(nspinor*nvctrp,norbp,nproc), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr

  call timing(iproc,'Un-TransSwitch','ON')

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
     call MPI_ALLTOALL(psi,nvctrp*norbp*nspinor,mpidtypw,  &
          work,nvctrp*norbp*nspinor,mpidtypw,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     if (present(outadd)) then
        call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,&
        work,outadd,nspinor)
     else
        call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,&
        work,psi,nspinor)
     end if
  else
     if(nspinor==4) then
        call psitransspi(nspinor,nvctrp,norb,psi,.false.)
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


subroutine psitransspi(nspinor,nvctrp,norb,psi,forward)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp,nspinor
  logical, intent(in) :: forward
  real(wp), dimension(nspinor*nvctrp,norb), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='psitransspi'
  integer :: i,iorb,ij,isp,i_all,i_stat
  real(wp), dimension(:,:,:), allocatable :: tpsit

  allocate(tpsit(nvctrp,nspinor,norb+ndebug),stat=i_stat)
  call memocc(i_stat,tpsit,'tpsit',subname)
  if(forward) then
    do iorb=1,norb
        do isp=1,nspinor
           do i=1,nvctrp
              tpsit(i,isp,iorb)=psi(i+(isp-1)*nvctrp,iorb)
           enddo
        enddo
     enddo
     if (nspinor == 2) then
        do iorb=1,norb
           do i=1,nvctrp
              psi(2*i-1,iorb)=tpsit(i,1,iorb)
              psi(2*i,iorb)=tpsit(i,2,iorb)
           enddo
        enddo
     else if (nspinor == 4) then
        do iorb=1,norb
           do i=1,nvctrp
              psi(2*i-1,iorb)=tpsit(i,1,iorb)
              psi(2*i,iorb)=tpsit(i,2,iorb)
              psi(2*i+2*nvctrp-1,iorb)=tpsit(i,3,iorb)
              psi(2*i+2*nvctrp,iorb)=tpsit(i,4,iorb)
           enddo
        enddo
     end if
  else
     if (nspinor == 2) then
        do iorb=1,norb
           do i=1,nvctrp
              tpsit(i,1,iorb)=psi(2*i-1,iorb)
              tpsit(i,2,iorb)=psi(2*i,iorb)
           enddo
        enddo
     else if (nspinor == 4) then
        do iorb=1,norb
           do i=1,nvctrp
              tpsit(i,1,iorb)=psi(2*i-1,iorb)
              tpsit(i,2,iorb)=psi(2*i,iorb)
              tpsit(i,3,iorb)=psi(2*i-1+2*nvctrp,iorb)
              tpsit(i,4,iorb)=psi(2*i+2*nvctrp,iorb)
           enddo
        enddo
     end if

     do iorb=1,norb
        ij=1
        do isp=1,nspinor
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
end subroutine psitransspi


!transposition of the arrays, variable version (non homogeneous)
subroutine transpose_v(iproc,nproc,orbs,wfd,comms,psi,&
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work)) then 
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for transposing in parallel"
        stop
     end if
     call switch_waves_v(nproc,orbs,&
          wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,psi,work)
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     if (present(outadd)) then
        call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
             outadd,comms%ncntt,comms%ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
     else
        call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
             psi,comms%ncntt,comms%ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
     end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     if(orbs%nspinor /= 1) then
        !for only one processor there is no need to transform this
        call psitransspi(orbs%nspinor,wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp,psi,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

end subroutine transpose_v

subroutine untranspose_v(iproc,nproc,orbs,wfd,comms,psi,&
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work)) then
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for untransposing in parallel"
        stop
     end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALLV(psi,comms%ncntt,comms%ndsplt,mpidtypw,  &
          work,comms%ncntd,comms%ndspld,mpidtypw,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     if (present(outadd)) then
        call unswitch_waves_v(nproc,orbs,&
             wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,work,outadd)
     else
        call unswitch_waves_v(nproc,orbs,&
             wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,work,psi)
     end if
  else
     if(orbs%nspinor /= 1) then
        call psitransspi(orbs%nspinor,wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp,psi,.false.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
end subroutine untranspose_v


subroutine switch_waves_v(nproc,orbs,nvctr,nvctr_par,psi,psiw)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,nvctr
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nproc,orbs%nkptsp), intent(in) :: nvctr_par
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt

  isorb=orbs%isorb+1
  isorbp=0
  ikpt=orbs%iskpts+1
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it=ind+i
                 psiw(it)=psi(ij,1,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikptsp)
           enddo
        enddo
     else if (orbs%nspinor == 2) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 psiw(it1)=psi(ij,1,iorb+isorbp)
                 psiw(it2)=psi(ij,2,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikptsp)
           enddo
        enddo
     else if (orbs%nspinor == 4) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 it3=ind+2*i+2*nvctr_par(j,ikptsp)-1
                 it4=ind+2*i+2*nvctr_par(j,ikptsp)
                 psiw(it1)=psi(ij,1,iorb+isorbp)
                 psiw(it2)=psi(ij,2,iorb+isorbp)
                 psiw(it3)=psi(ij,3,iorb+isorbp)
                 psiw(it4)=psi(ij,4,iorb+isorbp)
                 ij=ij+1
              enddo
              ijproc=ijproc+nvctr_par(j,ikptsp)
           enddo
        enddo
     end if
     !update k-point
     ikpt=ikpt+1
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
end subroutine switch_waves_v

subroutine unswitch_waves_v(nproc,orbs,nvctr,nvctr_par,psiw,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,nvctr
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nproc,orbs%nkptsp), intent(in) :: nvctr_par
  real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(in) :: psiw
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(out) :: psi
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt


  isorb=orbs%isorb+1
  isorbp=0
  ikpt=orbs%iskpts+1
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it=ind+i
                 psi(ij,orbs%nspinor,iorb+isorbp)=psiw(it)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikptsp)
           end do
        end do
     else if (orbs%nspinor == 2) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 psi(ij,1,iorb+isorbp)=psiw(it1)
                 psi(ij,2,iorb+isorbp)=psiw(it2)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikptsp)
           end do
        end do
     else if (orbs%nspinor == 4) then
        do iorb=1,norbp_kpt
           ij=1
           ijproc=0
           do j=1,nproc
              ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
                   ispsi
              do i=1,nvctr_par(j,ikptsp)
                 it1=ind+2*i-1
                 it2=ind+2*i
                 it3=ind+2*i+2*nvctr_par(j,ikptsp)-1
                 it4=ind+2*i+2*nvctr_par(j,ikptsp)
                 psi(ij,1,iorb+isorbp)=psiw(it1)
                 psi(ij,2,iorb+isorbp)=psiw(it2)
                 psi(ij,3,iorb+isorbp)=psiw(it3)
                 psi(ij,4,iorb+isorbp)=psiw(it4)
                 ij=ij+1
              end do
              ijproc=ijproc+nvctr_par(j,ikptsp)
           end do
        end do
     end if
     !update k-point
     ikpt=ikpt+1
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
  
end subroutine unswitch_waves_v

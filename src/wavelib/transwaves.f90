!!****f* BigDFT/psitransspi
!! COPYRIGHT
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!! 
subroutine psitransspi(nvctrp,orbs,psi,forward)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nvctrp
  logical, intent(in) :: forward
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(orbs%nspinor*nvctrp,orbs%norb,orbs%nkpts), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='psitransspi'
  integer :: i,iorb,isp,i_all,i_stat,ikpts
  real(wp), dimension(:,:,:,:), allocatable :: tpsit

  allocate(tpsit(nvctrp,orbs%nspinor,orbs%norb,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,tpsit,'tpsit',subname)
  if(forward) then
     !we can use dcopy here
     do ikpts=1,orbs%nkpts
        do iorb=1,orbs%norb
           do isp=1,orbs%nspinor
              do i=1,nvctrp
                 tpsit(i,isp,iorb,ikpts)=psi(i+(isp-1)*nvctrp,iorb,ikpts)
              enddo
           enddo
        enddo
     end do
     if (orbs%nspinor == 2) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 psi(2*i-1,iorb,ikpts)=tpsit(i,1,iorb,ikpts)
                 psi(2*i,iorb,ikpts)=tpsit(i,2,iorb,ikpts)
              enddo
           enddo
        end do
     else if (orbs%nspinor == 4) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 psi(2*i-1,iorb,ikpts)=tpsit(i,1,iorb,ikpts)
                 psi(2*i,iorb,ikpts)=tpsit(i,2,iorb,ikpts)
                 psi(2*i+2*nvctrp-1,iorb,ikpts)=tpsit(i,3,iorb,ikpts)
                 psi(2*i+2*nvctrp,iorb,ikpts)=tpsit(i,4,iorb,ikpts)
              enddo
           enddo
        end do
     end if
  else
     if (orbs%nspinor == 2) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 tpsit(i,1,iorb,ikpts)=psi(2*i-1,iorb,ikpts)
                 tpsit(i,2,iorb,ikpts)=psi(2*i,iorb,ikpts)
              enddo
           enddo
        end do
     else if (orbs%nspinor == 4) then
        do ikpts=1,orbs%nkpts
           do iorb=1,orbs%norb
              do i=1,nvctrp
                 tpsit(i,1,iorb,ikpts)=psi(2*i-1,iorb,ikpts)
                 tpsit(i,2,iorb,ikpts)=psi(2*i,iorb,ikpts)
                 tpsit(i,3,iorb,ikpts)=psi(2*i-1+2*nvctrp,iorb,ikpts)
                 tpsit(i,4,iorb,ikpts)=psi(2*i+2*nvctrp,iorb,ikpts)
              enddo
           enddo
        end do
     end if

     !here we can use dcopy
     do ikpts=1,orbs%nkpts
        do iorb=1,orbs%norb
           do isp=1,orbs%nspinor
              do i=1,nvctrp
                 psi(i+(isp-1)*nvctrp,iorb,ikpts)=tpsit(i,isp,iorb,ikpts)
              enddo
           enddo
        enddo
     end do
  end if

  i_all=-product(shape(tpsit))*kind(tpsit)
  deallocate(tpsit,stat=i_stat)
  call memocc(i_stat,i_all,'tpsit',subname)
END SUBROUTINE psitransspi
!!***


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
        call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.true.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

END SUBROUTINE transpose_v

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
        !if(iproc == 0) 
             write(*,'(1x,a)')&
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
        call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.false.)
     end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
END SUBROUTINE untranspose_v


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
  ispsi=0
  do ikptsp =1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

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
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
END SUBROUTINE switch_waves_v

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
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt,ierr

  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

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
     !update starting orbitals
     isorb=isorb+norbp_kpt
     isorbp=isorbp+norbp_kpt
     !and starting point for psi
     ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
  
END SUBROUTINE unswitch_waves_v










!!  ****************************************************************************************
!!  ****************************************************************************************
!!  ****************************************************************************************


!transposition of the arrays, variable version (non homogeneous)
!subroutine transpose_vLIN(iproc,nproc,orbs,comms,psi,lr,&
!     work,outadd) !optional
subroutine transpose_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, psi, lr, newComm, &
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, lproc, uproc, norbPerGroup, newComm
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors):: lr
  type(communications_arrays), intent(in) :: comms
  !real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: psi
  !real(8),dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor):: psi
  real(8),dimension(orbs%npsidim):: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr, nproc

  nproc=uproc-lproc+1

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work)) then 
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for transposing in parallel"
        stop
     end if
  
     !call switch_waves_v(nproc,orbs,&
     !     wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,psi,work)
     !call switch_waves_vLIN(iproc,nproc, orbs, comms, lr, psi, work)
     call switch_waves_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, lr, psi, work)

     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !if (present(outadd)) then
     !   call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
     !        outadd,comms%ncntt,comms%ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
     !else
     !   call MPI_ALLTOALLV(work,comms%ncntd,comms%ndspld,mpidtypw, &
     !        psi,comms%ncntt,comms%ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
     !end if
!write(*,'(a,i3,3x,8i9)') 'in transpose_vLIN: iproc, comms%ncnttLIN(0:nproc-1)', iproc, comms%ncnttLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in transpose_vLIN: iproc, comms%ndspltLIN(0:nproc-1)', iproc, comms%ndspltLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in transpose_vLIN: iproc, comms%ncntdLIN(0:nproc-1)', iproc, comms%ncntdLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in transpose_vLIN: iproc, comms%ndspldLIN(0:nproc-1)', iproc, comms%ndspldLIN(0:nproc-1)
call mpi_barrier(newComm, ierr)
     if (present(outadd)) then
        !call MPI_ALLTOALLV(work,comms%ncntdLIN,comms%ndspldLIN,mpidtypw, &
        !     outadd,comms%ncnttLIN,comms%ndspltLIN,mpidtypw,MPI_COMM_WORLD,ierr)
        call mpi_alltoallv(work, comms%ncntdLIN, comms%ndspldLIN, mpidtypw, &
             outadd, comms%ncnttLIN, comms%ndspltLIN, mpidtypw, newComm, ierr)
     else
        !call MPI_ALLTOALLV(work,comms%ncntdLIN,comms%ndspldLIN,mpidtypw, &
        !     psi,comms%ncnttLIN,comms%ndspltLIN,mpidtypw,MPI_COMM_WORLD,ierr)
        call mpi_alltoallv(work, comms%ncntdLIN, comms%ndspldLIN, mpidtypw, &
             psi, comms%ncnttLIN, comms%ndspltLIN, mpidtypw, newComm, ierr)
     end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     write(*,*) 'ERROR: transpose_vLIN not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   !for only one processor there is no need to transform this
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.true.)
     !!end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

END SUBROUTINE transpose_vLIN


     !call untranspose_vLIN(iproc,nproc,orbs,comms,psi,lr,work=psiWorkPointer)
!subroutine untranspose_vLIN(iproc,nproc,orbs,comms,psi,lr,&
!     work,outadd) !optional
subroutine untranspose_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, psi, lr, newComm, &
     work,outadd) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,lproc, uproc, norbPerGroup, newComm
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(locreg_descriptors):: lr
  !real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
  !real(8),dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor):: psi
  real(8),dimension(orbs%npsidim):: psi
  real(wp), dimension(:), pointer, optional :: work
  real(wp), dimension(*), intent(out), optional :: outadd
  !local variables
  integer :: ierr, nproc

  nproc=uproc-lproc+1
  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !control check
     if (.not. present(work) .or. .not. associated(work)) then
        !if(iproc == 0) 
             write(*,'(1x,a)')&
             "ERROR: Unproper work array for untransposing in parallel"
        stop
     end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !call MPI_ALLTOALLV(psi,comms%ncntt,comms%ndsplt,mpidtypw,  &
     !     work,comms%ncntd,comms%ndspld,mpidtypw,MPI_COMM_WORLD,ierr)
!write(*,'(a,i3,3x,8i9)') 'in untranspose_vLIN: iproc, comms%ncnttLIN(0:nproc-1)', iproc, comms%ncnttLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in untranspose_vLIN: iproc, comms%ndspltLIN(0:nproc-1)', iproc, comms%ndspltLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in untranspose_vLIN: iproc, comms%ncntdLIN(0:nproc-1)', iproc, comms%ncntdLIN(0:nproc-1)
!write(*,'(a,i3,3x,8i9)') 'in untranspose_vLIN: iproc, comms%ndspldLIN(0:nproc-1)', iproc, comms%ndspldLIN(0:nproc-1)
call mpi_barrier(newComm, ierr)
     !call MPI_ALLTOALLV(psi,comms%ncnttLIN,comms%ndspltLIN,mpidtypw,  &
     !     work,comms%ncntdLIN,comms%ndspldLIN,mpidtypw,MPI_COMM_WORLD,ierr)
     call mpi_alltoallv(psi, comms%ncnttLIN, comms%ndspltLIN, mpidtypw,  &
          work, comms%ncntdLIN, comms%ndspldLIN, mpidtypw, newComm, ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     !if (present(outadd)) then
     !   call unswitch_waves_v(nproc,orbs,&
     !        wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,work,outadd)
     !else
     !   call unswitch_waves_v(nproc,orbs,&
     !        wfd%nvctr_c+7*wfd%nvctr_f,comms%nvctr_par,work,psi)
     !end if
     if (present(outadd)) then
        call unswitch_waves_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, lr, work, outadd)
     else
        call unswitch_waves_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, lr, work, psi)
     end if
  else
     write(*,*) 'ERROR: untranspose_vLIN not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.false.)
     !!end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
END SUBROUTINE untranspose_vLIN




subroutine switch_waves_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, lr, psi, psiw)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,lproc,uproc, norbPerGroup
  type(communications_arrays), intent(in) :: comms
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors):: lr
  !integer, dimension(nproc,orbs%nkptsp), intent(in) :: nvctr_par
  !real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
  !!! CHECK THIS !!!
  !real(wp), dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor), intent(in) :: psi
  real(wp), dimension(orbs%npsidim), intent(in) :: psi
  !! THIS HAS TO BE CHANGED AS WELL !!
  !real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(out) :: psiw
  !real(wp), dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor), intent(out) :: psiw
  real(wp), dimension(orbs%npsidim), intent(out) :: psiw
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp, nproc, jproc, kproc
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt
integer:: k, ii, idebug
real(8):: tt

  nproc=uproc-lproc+1

  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0


if(orbs%nkptsp>1) then
  write(*,*) 'ERROR: more than 1 k-point!'
  stop
end if
  do ikptsp =1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
     ij=1
      !do jproc=lproc,uproc
        !do iorb=1,orbs%norb_par(jproc)
        do iorb=1,orbs%norbp
           ijproc=0
           do j=lproc,uproc
              ii=0
               do k=1,iorb-1
                   ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
               end do
              ind=ii+ijproc*orbs%nspinor+ispsi
              do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 it=ind+i
                 psiw(it)=psi(ij)
                 ij=ij+1
              enddo
              !ijproc=ijproc+comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
              do k=1,orbs%norbp
                  ijproc=ijproc+comms%nvctr_parLIN(k,iproc,j,ikptsp)
              end do
           enddo
        enddo
      !end do
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
        !!do iorb=1,norbp_kpt
        !!   ij=1
        !!   ijproc=0
        !!   do j=1,nproc
        !!      ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
        !!           ispsi
        !!      do i=1,nvctr_par(j,ikptsp)
        !!         it1=ind+2*i-1
        !!         it2=ind+2*i
        !!         psiw(it1)=psi(ij,1,iorb+isorbp)
        !!         psiw(it2)=psi(ij,2,iorb+isorbp)
        !!         ij=ij+1
        !!      enddo
        !!      ijproc=ijproc+nvctr_par(j,ikptsp)
        !!   enddo
        !!enddo
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
        !!do iorb=1,norbp_kpt
        !!   ij=1
        !!   ijproc=0
        !!   do j=1,nproc
        !!      ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
        !!           ispsi
        !!      do i=1,nvctr_par(j,ikptsp)
        !!         it1=ind+2*i-1
        !!         it2=ind+2*i
        !!         it3=ind+2*i+2*nvctr_par(j,ikptsp)-1
        !!         it4=ind+2*i+2*nvctr_par(j,ikptsp)
        !!         psiw(it1)=psi(ij,1,iorb+isorbp)
        !!         psiw(it2)=psi(ij,2,iorb+isorbp)
        !!         psiw(it3)=psi(ij,3,iorb+isorbp)
        !!         psiw(it4)=psi(ij,4,iorb+isorbp)
        !!         ij=ij+1
        !!      enddo
        !!      ijproc=ijproc+nvctr_par(j,ikptsp)
        !!   enddo
        !!enddo
     end if
     !! NOT YET IMPLEMENTED !!
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
END SUBROUTINE switch_waves_vLIN




!subroutine unswitch_waves_vLIN(nproc,orbs,nvctr,nvctr_par,psiw,psi)
subroutine unswitch_waves_vLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, lr, psiw, psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, lproc, uproc, norbPerGroup
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(locreg_descriptors):: lr
  !real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(in) :: psiw
  !real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(out) :: psi
  !real(wp), dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor), intent(out) :: psi
  real(wp), dimension(orbs%npsidim), intent(out) :: psi
  !! THIS HAS TO BE CHANGED AS WELL !!
  !real(wp), dimension(orbs%nspinor*nvctr*orbs%norbp), intent(out) :: psiw
  !real(wp), dimension((sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_c)+7*sum(lr%wfdLIN(1:orbs%norbp,iproc)%nvctr_f))*orbs%nspinor), intent(in) :: psiw
  real(wp), dimension(orbs%npsidim), intent(in) :: psiw
  !local variables
  integer :: iorb,i,j,ij,ijproc,ind,it,it1,it2,it3,it4,ikptsp, nproc, jproc
  integer :: isorb,isorbp,ispsi,norbp_kpt,ikpt,ierr
integer:: k, ii

  nproc=uproc-lproc+1


  !! CAN I DELETE THIS? !!
  !call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)
     !if (ikpt < orbs%iokpt(1) .or. ikpt > orbs%iokpt(orbs%norbp)) cycle

     !calculate the number of orbitals belonging to k-point ikptstp
     !calculate to which k-point it belongs
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
        !ij=1
       !do jproc=lproc,uproc
        !do iorb=1,norbp_kpt
        !do iorb=1,norbPerGroup
        ij=1
        do iorb=1,orbs%norbp
           !ij=1
           ijproc=0
           !do j=0,nproc-1
           do j=lproc,uproc
              ii=0
              do k=1,iorb-1
                  ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
              end do
              !ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
              !     ispsi
              ind=ii+ijproc*orbs%nspinor+ispsi
              !do i=1,nvctr_par(j,ikptsp)
              do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 it=ind+i
                 !psi(ij,orbs%nspinor,iorb+isorbp)=psiw(it)
                 if(ij>size(psi)) then
                     write(*,'(a,2i12)') 'ERROR: ij>size(psi), ij, size(psi)',ij, size(psi)
                     stop
                 end if
                 if(it>size(psiw)) then
                     write(*,'(a,2i12)') 'ERROR: it>size(psiw), it, size(psiw)',it, size(psiw)
                     stop
                 end if
                 psi(ij)=psiw(it)
                 ij=ij+1
              end do
              !ijproc=ijproc+nvctr_par(j,ikptsp)
              !do k=1,norbp_kpt
              !ijproc=ijproc+comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
              do k=1,orbs%norbp
                  ijproc=ijproc+comms%nvctr_parLIN(k,iproc,j,ikptsp)
              end do
              !ijproc=ijproc+comms%nvctr_parLIN(iorb,j,ikptsp)
           end do
        end do
      !end do
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
        !!do iorb=1,norbp_kpt
        !!   ij=1
        !!   ijproc=0
        !!   do j=1,nproc
        !!      ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
        !!           ispsi
        !!      do i=1,nvctr_par(j,ikptsp)
        !!         it1=ind+2*i-1
        !!         it2=ind+2*i
        !!         psi(ij,1,iorb+isorbp)=psiw(it1)
        !!         psi(ij,2,iorb+isorbp)=psiw(it2)
        !!         ij=ij+1
        !!      end do
        !!      ijproc=ijproc+nvctr_par(j,ikptsp)
        !!   end do
        !!end do
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
        !!do iorb=1,norbp_kpt
        !!   ij=1
        !!   ijproc=0
        !!   do j=1,nproc
        !!      ind=(iorb-1)*orbs%nspinor*nvctr_par(j,ikptsp)+ijproc*orbs%nspinor*norbp_kpt+&
        !!           ispsi
        !!      do i=1,nvctr_par(j,ikptsp)
        !!         it1=ind+2*i-1
        !!         it2=ind+2*i
        !!         it3=ind+2*i+2*nvctr_par(j,ikptsp)-1
        !!         it4=ind+2*i+2*nvctr_par(j,ikptsp)
        !!         psi(ij,1,iorb+isorbp)=psiw(it1)
        !!         psi(ij,2,iorb+isorbp)=psiw(it2)
        !!         psi(ij,3,iorb+isorbp)=psiw(it3)
        !!         psi(ij,4,iorb+isorbp)=psiw(it4)
        !!         ij=ij+1
        !!      end do
        !!      ijproc=ijproc+nvctr_par(j,ikptsp)
        !!   end do
        !!end do
     end if
     !!! NOT YET IMPLEMENTED
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do

  
END SUBROUTINE unswitch_waves_vLIN


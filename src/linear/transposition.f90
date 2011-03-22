subroutine transpose_vLIN(iproc, lproc, uproc, orbs, comms, psi, newComm, &
     work,outadd) !optional
! Purpose:
! ========
!   Transposes the wave function(s) contained in psi. Each wave function may have its
!   own localization region. The transposition is done only among the processes
!   in the MPI communicator newComm.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc              process ID
!     lproc              lowest process ID of the current MPI communicator
!     uproc              highest process ID of the current MPI communicator
!     orbs               type describing the orbitals
!     comms              type containing the communications parameters
!     newComm            the current MPI communicator
!   Input / Output arguments:
!   -------------------------
!     psi                the orbitals to be transposed.
!     work (optional)    work array
!     outadd (optional)  if present, the transposed wave function will be 
!                        assigned to outadd instead of psi
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, lproc, uproc, newComm
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in):: comms
real(8),dimension(orbs%npsidim), intent(in out):: psi
real(wp),dimension(:),pointer,optional:: work
real(wp),dimension(*),intent(out),optional:: outadd

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     ! Control check
     if (.not. present(work) .or. .not. associated(work)) then 
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for transposing in parallel"
        stop
     end if
  
     ! Rearrange the orbitals on the current process such that they can be communicated
     ! more easily.
     call switch_waves_vLIN(iproc, lproc, uproc, orbs, comms, psi, work)

     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     if (present(outadd)) then
        call mpi_alltoallv(work, comms%ncntdLIN, comms%ndspldLIN, mpidtypw, &
             outadd, comms%ncnttLIN, comms%ndspltLIN, mpidtypw, newComm, ierr)
     else
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


subroutine untranspose_vLIN(iproc, lproc, uproc, orbs, comms, psi, newComm, &
     work,outadd) !optional
! Purpose:
! ========
!   Untransposes the wave function(s) contained in psi. Each wave function may have its
!   own localization region. The untransposition is done only among the processes
!   in the MPI communicator newComm.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc              process ID
!     lproc              lowest process ID of the current MPI communicator
!     uproc              highest process ID of the current MPI communicator
!     orbs               type describing the orbitals
!     comms              type containing the communications parameters
!     newComm            the current MPI communicator
!   Input / Output arguments:
!   -------------------------
!     psi                the orbitals to be untransposed.
!     work (optional)    work array
!     outadd (optional)  if present, the untransposed wave function will be 
!                        assigned to outadd instead of psi
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc,lproc, uproc, newComm
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in):: comms
real(8),dimension(orbs%npsidim),intent(in out):: psi
real(wp),dimension(:),pointer,optional :: work
real(wp),dimension(*),intent(out),optional :: outadd

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1
  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     ! Control check
     if (.not. present(work) .or. .not. associated(work)) then
        if(iproc == 0) write(*,'(1x,a)')&
             "ERROR: Unproper work array for untransposing in parallel"
        stop
     end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call mpi_alltoallv(psi, comms%ncnttLIN, comms%ndspltLIN, mpidtypw,  &
          work, comms%ncntdLIN, comms%ndspldLIN, mpidtypw, newComm, ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     if (present(outadd)) then
        call unswitch_waves_vLIN(iproc, lproc, uproc, orbs, comms, work, outadd)
     else
        call unswitch_waves_vLIN(iproc, lproc, uproc, orbs, comms, work, psi)
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




subroutine switch_waves_vLIN(iproc, lproc, uproc, orbs, comms, psi, psiw)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,lproc,uproc
  type(communications_arrays), intent(in) :: comms
  type(orbitals_data), intent(in) :: orbs
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
subroutine unswitch_waves_vLIN(iproc, lproc, uproc, orbs, comms, psiw, psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, lproc, uproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
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


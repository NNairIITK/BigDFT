!> @file 
!!   Transposition in linear
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
 

!> Transposes the wave function(s) contained in psi. Each wave function may have its
!!  own localization region. The transposition is done only among the processes
!!  in the MPI communicator newComm.
!
!! Calling arguments:
!! ==================
!!  Input arguments:
!!  ----------------
!!    @param iproc              process ID
!!    @param lproc              lowest process ID of the current MPI communicator
!!    @param uproc              highest process ID of the current MPI communicator
!!    @param orbs               type describing the orbitals
!!    @param comms              type containing the communications parameters
!!    @param newComm            the current MPI communicator
!!  Input / Output arguments:
!!  -------------------------
!!    @param psi                the orbitals to be transposed.
!!    @param work (optional)    work array
!!    @param outadd (optional)  if present, the transposed wave function will be 
!!                              assigned to outadd instead of psi
subroutine transpose_linear(iproc, lproc, uproc, orbs, collComms, psi, newComm, work)
     !work,outadd) !optional

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc, newComm
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(in out) :: psi
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)) :: work
!real(wp),dimension(*),intent(out),optional :: outadd

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !! Control check
     !if (.not. present(work) .or. .not. associated(work)) then 
     !   if(iproc == 0) write(*,'(1x,a)')&
     !        "ERROR: Unproper work array for transposing in parallel"
     !   stop
     !end if
  
     ! Rearrange the orbitals on the current process such that they can be communicated
     ! more easily.
     call switch_waves_linear(iproc, lproc, uproc, orbs, collComms, psi, work)

     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !if (present(outadd)) then
     !   call mpi_alltoallv(work, collComms%sendcnts, collComms%senddspls, mpidtypw, &
     !        outadd, collComms%recvcnts, collComms%recvdspls, mpidtypw, newComm, ierr)
     !else
        call mpi_alltoallv(work, collComms%sendcnts, collComms%senddspls, mpidtypw, &
             psi, collComms%recvcnts, collComms%recvdspls, mpidtypw, newComm, ierr)
     !end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     write(*,*) 'ERROR: transpose_linear not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   !for only one processor there is no need to transform this
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.true.)
     !!end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')

END SUBROUTINE transpose_linear


!> Untransposes the wave function(s) contained in psi. Each wave function may have its
!!  own localization region. The untransposition is done only among the processes
!!  in the MPI communicator newComm.
!
!! Calling arguments:
!! ==================
!!  Input arguments:
!!  ----------------
!!    @param iproc              process ID
!!    @param lproc              lowest process ID of the current MPI communicator
!!    @param uproc              highest process ID of the current MPI communicator
!!    @param orbs               type describing the orbitals
!!    @param comms              type containing the communications parameters
!!    @param newComm            the current MPI communicator
!!  Input / Output arguments:
!!  -------------------------
!!    @param psi                the orbitals to be untransposed.
!!    @param work (optional)    work array
!!    @param outadd (optional)  if present, the untransposed wave function will be 
!!                              assigned to outadd instead of psi
subroutine untranspose_linear(iproc, lproc, uproc, orbs, collComms, psi, newComm, work)
!     work,outadd) !optional

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc,lproc, uproc, newComm
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in out) :: psi
!real(wp),dimension(:),pointer,optional :: work
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)) :: work
!real(wp),dimension(*),intent(out),optional :: outadd

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1
  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     ! Control check
     !if (.not. present(work) .or. .not. associated(work)) then
     !   if(iproc == 0) write(*,'(1x,a)')&
     !        "ERROR: Unproper work array for untransposing in parallel"
     !   stop
     !end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !call mpi_alltoallv(psi, collComms%sendcnts, collComms%senddspls, mpidtypw,  &
     !     work, collComms%recvcnts, collComms%recvdspls, mpidtypw, newComm, ierr)
     call mpi_alltoallv(psi, collComms%recvcnts, collComms%recvdspls, mpidtypw,  &
          work, collComms%sendcnts, collComms%senddspls, mpidtypw, newComm, ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     !if (present(outadd)) then
     !   call unswitch_waves_linear(iproc, lproc, uproc, orbs, collComms, work, outadd)
     !else
        call unswitch_waves_linear(iproc, lproc, uproc, orbs, collComms, work, psi)
     !end if
  else
     write(*,*) 'ERROR: untranspose_vLIN not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.false.)
     !!end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
END SUBROUTINE untranspose_linear


!> This subroutine rearranges the orbitals such that they can be transposed using
!!  a single call to mpi_alltoallv.
!!  Here is an example how it works:
!!
!!      process 0                process 1            process 2        process 3        process 4
!!  1   9       1  13   |   17  25      17  29   |   33      33   |   41      41   |   49      49
!!  2  10       2  14   |   18  26      18  30   |   34      34   |   42      42   |   50      50
!!  3  11       9   7   |   19  27      25  23   |   35      35   |   43      43   |   51      51
!!  4  12      10   0   |   20  28      26   0   |   36      36   |   44      44   |   52      52
!!  5  13  =>   3   8   |   21  29  =>  19  24   |   37  =>  37   |   45  =>  45   |   53  =>  53
!!  6  14       4   0   |   22  30      20   0   |   38      38   |   46      46   |   54      54
!!  7  15      11  15   |   23  31      27  31   |   39      39   |   47      47   |   55      55
!!  8  16      12   0   |   24  32      28   0   |   40       0   |   48       0   |   56       0 
!!  0   0       5  16   |    0   0      21  32   |    0      40   |    0      48   |    0      56
!!  0   0       6   0   |    0   0      22   0   |    0       0   |    0       0   |    0       0 
!!
!!
!!  After this step, we can transpose it with one call to mpi_alltoallv:
!!
!!      process 0               process 1                process 2                process 3                process 4
!! 1  9 17 25 33 41 49  |   3 11 19 27 35 43 51  |   5 13 21 29 37 45 53  |   7 15 23 31 39 47 55  |   8 16 24 32 40 48 56
!! 2 10 18 26 34 42 50  |   4 12 20 28 36 44 52  |   6 14 22 30 38 46 54  |   0  0  0  0  0  0  0  |   0  0  0  0  0  0  0
!!
!!
!! Calling arguments:
!! ==================
!!  Input arguments:
!!  ----------------
!!    @param iproc      process ID
!!    @param lproc      lowest process ID in the current communicator
!!    @param uproc      highest process ID in the current communicator
!!    @param orbs       type describing the orbitals
!!    @param comms      type containing the communications parameters
!!    @param psi        the orbitals to be rearranged
!!  Output arguments:
!!  -----------------
!!    @param psiw       the rearranged orbitals
subroutine switch_waves_linear(iproc, lproc, uproc, orbs, collComms, psi, psiw)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc
type(collectiveComms), intent(in) :: collComms
type(orbitals_data),intent(in) :: orbs
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in) :: psi
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: psiw

! Local variables
integer :: iorb, i, j, ij, ijproc, ind, it, ikptsp, nproc
integer :: isorb, isorbp, ispsi, norbp_kpt, ikpt, kk, iiorb
integer :: k, ii

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  ! Related to k-points...
  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0


  if(orbs%nkptsp>1) then
    write(*,'(1x,a)') 'ERROR: more than 1 k-point!'
    stop
  end if
  
  ! This k-point loop is fake.
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)

     ! Calculate the number of orbitals belonging to k-point ikptstp.
     ! Calculate to which k-point it belongs.
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
         ij=1
         ! Make a loop over all orbitals belonging to iproc.
         do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ijproc=0
             ! Make a loop over all processes in the communicator.
             do j=lproc,uproc
                 ! Go to the starting index ind in psiw. This consists out of two parts:
                 !   - ii gives the amount of all orbitals up to iorb of process iproc which will
                 !     be passed from process iproc to process j
                 !   - ijproc gives the amount of all orbitals belonging to iproc which will be
                 !     passed to the processes up to j
                 ! ispsi is related to k-points and not used at the moment.
                 ii=0
                 do k=1,iorb-1
                      kk=orbs%isorb+k
                      !ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                      ii=ii+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
                 ind=ii+ijproc+ispsi
                 ! Now copy the values from psi to psiw.
                 !do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 do i=1,collComms%nvctr_par(iiorb,j)
                     it=ind+i
                     psiw(it)=psi(ij)
                     ij=ij+1
                 enddo
                 ! Update the ijproc counter.
                 do k=1,orbs%norbp
                     kk=orbs%isorb+k
                     !ijproc=ijproc+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ijproc=ijproc+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
             enddo
         enddo
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
     end if
     !! NOT YET IMPLEMENTED !!
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
END SUBROUTINE switch_waves_linear


!> This subroutine rearranges the orbitals back. As an input it takes the psiw in the form
!! which is used for the mpi_alltoallv.
!!  The idea is the same as in switch_waves_vLIN, just now the other way around.
!
!! Calling arguments:
!! ==================
!!  Input arguments:
!!  ----------------
!!    @param iproc      process ID
!!    @param lproc      lowest process ID in the current communicator
!!    @param uproc      highest process ID in the current communicator
!!    @param orbs       type describing the orbitals
!!    @param comms      type containing the communications parameters
!!    @param psiw       the orbitals to be rearranged
!!  Output arguments:
!!  -----------------
!!    @param psi        the rearranged orbitals
subroutine unswitch_waves_linear(iproc, lproc, uproc, orbs, collComms, psiw, psi)

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in) :: psiw
real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: psi

! Local variables
integer :: iorb, i, j, ij, ijproc, ind, it, ikptsp, nproc
integer :: isorb, isorbp, ispsi, norbp_kpt, ikpt
integer :: k, ii, kk, iiorb

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  ! This subroutine is not yet implemented for k-points.
  if(orbs%nkptsp>1) then
    write(*,'(1x,a)') 'ERROR: more than 1 k-point!'
    stop
  end if

  ! Related to k-points...
  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0


  ! This k-point loop is fake.
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)

     ! Calculate the number of orbitals belonging to k-point ikptstp.
     ! Calculate to which k-point it belongs.
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
         ij=1
         ! Make a loop over all orbitals belonging to iproc.
         do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ijproc=0
             ! Make a loop over all processes in the communicator.
             do j=lproc,uproc
                 ! Go to the starting index ind in psiw. This consists out of two parts:
                 !   - ii gives the amount of all orbitals up to iorb of process iproc which have been
                 !     be passed from process iproc to process j
                 !   - ijproc gives the amount of all orbitals belonging to iproc which have been
                 !     passed to the processes up to j
                 ! ispsi is related to k-points and not used at the moment.
                 ii=0
                 do k=1,iorb-1
                     kk=orbs%isorb+k
                     !ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ii=ii+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
                 ind=ii+ijproc+ispsi
                 !do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 do i=1,collComms%nvctr_par(iiorb,j)
                     it=ind+i
                     psi(ij)=psiw(it)
                     ij=ij+1
                 end do
                 do k=1,orbs%norbp
                     kk=orbs%isorb+k
                     !ijproc=ijproc+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ijproc=ijproc+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
             end do
         end do
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
     end if
     !!! NOT YET IMPLEMENTED
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do

  
END SUBROUTINE unswitch_waves_linear


!>  Transposes the wave function(s) contained in psi. Each wave function may have its
!!  own localization region. The transposition is done only among the processes
!!  in the MPI communicator newComm.
!!
!! Calling arguments:
!! ==================
!!   Input arguments:
!!   ----------------
!!     @param iproc              process ID
!!     @param lproc              lowest process ID of the current MPI communicator
!!     @param uproc              highest process ID of the current MPI communicator
!!     @param orbs               type describing the orbitals
!!     @param comms              type containing the communications parameters
!!     @param newComm            the current MPI communicator
!!   Input / Output arguments:
!!   -------------------------
!!     @param psi                the orbitals to be transposed.
!!     @param work (optional)    work array
!!     @param outadd (optional)  if present, the transposed wave function will be 
!!                               assigned to outadd instead of psi
subroutine transpose_linear_int(iproc, lproc, uproc, orbs, collComms, psi, newComm, work)
     !work,outadd) !optional

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc, newComm
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(in out) :: psi
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)) :: work
!real(wp),dimension(*),intent(out),optional :: outadd

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     !! Control check
     !if (.not. present(work) .or. .not. associated(work)) then 
     !   if(iproc == 0) write(*,'(1x,a)')&
     !        "ERROR: Unproper work array for transposing in parallel"
     !   stop
     !end if
  
     ! Rearrange the orbitals on the current process such that they can be communicated
     ! more easily.
     call switch_waves_linear_int(iproc, lproc, uproc, orbs, collComms, psi, work)
     !!do ierr=1,orbs%npsidim
     !!    write(50,*) ierr, work(ierr)
     !!end do

     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !if (present(outadd)) then
     !   call mpi_alltoallv(work, collComms%sendcnts, collComms%senddspls, mpidtypw, &
     !        outadd, collComms%recvcnts, collComms%recvdspls, mpidtypw, newComm, ierr)
     !else
        call mpi_alltoallv(work, collComms%sendcnts, collComms%senddspls, mpi_integer, &
             psi, collComms%recvcnts, collComms%recvdspls, mpi_integer, newComm, ierr)
     !end if
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
  else
     write(*,*) 'ERROR: transpose_linear not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   !for only one processor there is no need to transform this
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.true.)
     !!end if
  end if
  !!if(iproc==0) then
  !!    do ierr=1,orbs%npsidim
  !!        write(51,*) ierr, psi(ierr)
  !!    end do
  !!end if

  call timing(iproc,'Un-TransSwitch','OF')

END SUBROUTINE transpose_linear_int


!>  Untransposes the wave function(s) contained in psi. Each wave function may have its
!!  own localization region. The untransposition is done only among the processes
!!  in the MPI communicator newComm.
!!
!! Calling arguments:
!! ==================
!!   Input arguments:
!!   ----------------
!!    @param iproc              process ID
!!    @param lproc              lowest process ID of the current MPI communicator
!!    @param uproc              highest process ID of the current MPI communicator
!!    @param orbs               type describing the orbitals
!!    @param comms              type containing the communications parameters
!!    @param newComm            the current MPI communicator
!!   Input / Output arguments:
!!   -------------------------
!!     @param psi                the orbitals to be untransposed.
!!     @param work (optional)    work array
!!     @param outadd (optional)  if present, the untransposed wave function will be 
!!                               assigned to outadd instead of psi
subroutine untranspose_linear_int(iproc, lproc, uproc, orbs, collComms, psi, newComm, work)
!     work,outadd) !optional

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc,lproc, uproc, newComm
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in out) :: psi
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)) :: work

! Local variables
integer :: ierr, nproc

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1
  call timing(iproc,'Un-TransSwitch','ON')

  if (nproc > 1) then
     ! Control check
     !if (.not. present(work) .or. .not. associated(work)) then
     !   if(iproc == 0) write(*,'(1x,a)')&
     !        "ERROR: Unproper work array for untransposing in parallel"
     !   stop
     !end if
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     !call mpi_alltoallv(psi, collComms%sendcnts, collComms%senddspls, mpidtypw,  &
     !     work, collComms%recvcnts, collComms%recvdspls, mpidtypw, newComm, ierr)
     call mpi_alltoallv(psi, collComms%recvcnts, collComms%recvdspls, mpi_integer,  &
          work, collComms%sendcnts, collComms%senddspls, mpi_integer, newComm, ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     !if (present(outadd)) then
     !   call unswitch_waves_linear(iproc, lproc, uproc, orbs, collComms, work, outadd)
     !else
        call unswitch_waves_linear_int(iproc, lproc, uproc, orbs, collComms, work, psi)
     !end if
  else
     write(*,*) 'ERROR: untranspose_vLIN not yet implemented for nproc==1'
     stop
     !!if(orbs%nspinor /= 1) then
     !!   call psitransspi(wfd%nvctr_c+7*wfd%nvctr_f,orbs,psi,.false.)
     !!end if
  end if

  call timing(iproc,'Un-TransSwitch','OF')
END SUBROUTINE untranspose_linear_int


!> This subroutine rearranges the orbitals such that they can be transposed using
!!  a single call to mpi_alltoallv.
!!  Here is an example how it works:
!!
!!      process 0                process 1            process 2        process 3        process 4
!!  1   9       1  13   |   17  25      17  29   |   33      33   |   41      41   |   49      49
!!  2  10       2  14   |   18  26      18  30   |   34      34   |   42      42   |   50      50
!!  3  11       9   7   |   19  27      25  23   |   35      35   |   43      43   |   51      51
!!  4  12      10   0   |   20  28      26   0   |   36      36   |   44      44   |   52      52
!!  5  13  =>   3   8   |   21  29  =>  19  24   |   37  =>  37   |   45  =>  45   |   53  =>  53
!!  6  14       4   0   |   22  30      20   0   |   38      38   |   46      46   |   54      54
!!  7  15      11  15   |   23  31      27  31   |   39      39   |   47      47   |   55      55
!!  8  16      12   0   |   24  32      28   0   |   40       0   |   48       0   |   56       0 
!!  0   0       5  16   |    0   0      21  32   |    0      40   |    0      48   |    0      56
!!  0   0       6   0   |    0   0      22   0   |    0       0   |    0       0   |    0       0 
!!
!!
!!  After this step, we can transpose it with one call to mpi_alltoallv:
!!
!!      process 0               process 1                process 2                process 3                process 4
!! 1  9 17 25 33 41 49  |   3 11 19 27 35 43 51  |   5 13 21 29 37 45 53  |   7 15 23 31 39 47 55  |   8 16 24 32 40 48 56
!! 2 10 18 26 34 42 50  |   4 12 20 28 36 44 52  |   6 14 22 30 38 46 54  |   0  0  0  0  0  0  0  |   0  0  0  0  0  0  0
!!
!!
!! Calling arguments:
!! ==================
!!  Input arguments:
!!  ----------------
!!    @param iproc      process ID
!!    @param lproc      lowest process ID in the current communicator
!!    @param uproc      highest process ID in the current communicator
!!    @param orbs       type describing the orbitals
!!    @param comms      type containing the communications parameters
!!    @param psi        the orbitals to be rearranged
!!  Output arguments:
!!  -----------------
!!    @param psiw       the rearranged orbitals
subroutine switch_waves_linear_int(iproc, lproc, uproc, orbs, collComms, psi, psiw)

use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc
type(collectiveComms), intent(in) :: collComms
type(orbitals_data),intent(in) :: orbs
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in) :: psi
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: psiw

! Local variables
integer :: iorb, i, j, ij, ijproc, ind, it, ikptsp, nproc
integer :: isorb, isorbp, ispsi, norbp_kpt, ikpt, kk, iiorb
integer :: k, ii

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  ! Related to k-points...
  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0


  if(orbs%nkptsp>1) then
    write(*,'(1x,a)') 'ERROR: more than 1 k-point!'
    stop
  end if
  
  ! This k-point loop is fake.
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)

     ! Calculate the number of orbitals belonging to k-point ikptstp.
     ! Calculate to which k-point it belongs.
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
         ij=1
         ! Make a loop over all orbitals belonging to iproc.
         do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ijproc=0
             ! Make a loop over all processes in the communicator.
             do j=lproc,uproc
                 ! Go to the starting index ind in psiw. This consists out of two parts:
                 !   - ii gives the amount of all orbitals up to iorb of process iproc which will
                 !     be passed from process iproc to process j
                 !   - ijproc gives the amount of all orbitals belonging to iproc which will be
                 !     passed to the processes up to j
                 ! ispsi is related to k-points and not used at the moment.
                 ii=0
                 do k=1,iorb-1
                      kk=orbs%isorb+k
                      !ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                      ii=ii+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
                 ind=ii+ijproc+ispsi
                 ! Now copy the values from psi to psiw.
                 !do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 do i=1,collComms%nvctr_par(iiorb,j)
                     it=ind+i
                     psiw(it)=psi(ij)
                     !if(iiorb==1) then
                     !    write(120,*) it, psiw(it)
                     !else if(iiorb==5) then
                     !    write(520,*) it, psiw(it)
                     !end if
                     ij=ij+1
                 enddo
                 ! Update the ijproc counter.
                 do k=1,orbs%norbp
                     kk=orbs%isorb+k
                     !ijproc=ijproc+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ijproc=ijproc+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
             enddo
         enddo
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
     end if
     !! NOT YET IMPLEMENTED !!
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do
END SUBROUTINE switch_waves_linear_int


!> This subroutine rearranges the orbitals back. As an input it takes the psiw in the form
!! which is used for the mpi_alltoallv.
!! The idea is the same as in switch_waves_vLIN, just now the other way around.
!!
!! Calling arguments:
!! ==================
!!   Input arguments:
!!   ----------------
!!     @param iproc      process ID
!!     @param lproc      lowest process ID in the current communicator
!!     @param uproc      highest process ID in the current communicator
!!     @param orbs       type describing the orbitals
!!     @param comms      type containing the communications parameters
!!     @param psiw       the orbitals to be rearranged
!!   Output arguments:
!!   -----------------
!!     @param psi        the rearranged orbitals
subroutine unswitch_waves_linear_int(iproc, lproc, uproc, orbs, collComms, psiw, psi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, lproc, uproc
type(orbitals_data),intent(in) :: orbs
type(collectiveComms),intent(in) :: collComms
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in) :: psiw
integer,dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: psi

! Local variables
integer :: iorb, i, j, ij, ijproc, ind, it, ikptsp, nproc
integer :: isorb, isorbp, ispsi, norbp_kpt, ikpt
integer :: k, ii, kk, iiorb

  ! Number of processes in the current communicator.
  nproc=uproc-lproc+1

  ! This subroutine is not yet implemented for k-points.
  if(orbs%nkptsp>1) then
    write(*,'(1x,a)') 'ERROR: more than 1 k-point!'
    stop
  end if

  ! Related to k-points...
  isorb=orbs%isorb+1
  isorbp=0
  ispsi=0


  ! This k-point loop is fake.
  do ikptsp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptsp !orbs%ikptsp(ikptsp)

     ! Calculate the number of orbitals belonging to k-point ikptstp.
     ! Calculate to which k-point it belongs.
     norbp_kpt=min(orbs%norb*ikpt,orbs%isorb+orbs%norbp)-isorb+1

     if(orbs%nspinor==1) then
         ij=1
         ! Make a loop over all orbitals belonging to iproc.
         do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ijproc=0
             ! Make a loop over all processes in the communicator.
             do j=lproc,uproc
                 ! Go to the starting index ind in psiw. This consists out of two parts:
                 !   - ii gives the amount of all orbitals up to iorb of process iproc which have been
                 !     be passed from process iproc to process j
                 !   - ijproc gives the amount of all orbitals belonging to iproc which have been
                 !     passed to the processes up to j
                 ! ispsi is related to k-points and not used at the moment.
                 ii=0
                 do k=1,iorb-1
                     kk=orbs%isorb+k
                     !ii=ii+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ii=ii+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
                 ind=ii+ijproc+ispsi
                 !do i=1,comms%nvctr_parLIN(iorb,iproc,j,ikptsp)
                 do i=1,collComms%nvctr_par(iiorb,j)
                     it=ind+i
                     psi(ij)=psiw(it)
                     ij=ij+1
                 end do
                 do k=1,orbs%norbp
                     kk=orbs%isorb+k
                     !ijproc=ijproc+orbs%nspinor*comms%nvctr_parLIN(k,iproc,j,ikptsp)
                     ijproc=ijproc+orbs%nspinor*collComms%nvctr_par(kk,j)
                 end do
             end do
         end do
     else if (orbs%nspinor == 2) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==2!'
        stop
     else if (orbs%nspinor == 4) then
        write(*,*) 'ERROR: not yet implemented for orbs%nspinor==4!'
        stop
     end if
     !!! NOT YET IMPLEMENTED
     !!!update starting orbitals
     !!isorb=isorb+norbp_kpt
     !!isorbp=isorbp+norbp_kpt
     !!!and starting point for psi
     !!ispsi=ispsi+orbs%nspinor*nvctr*norbp_kpt
  end do

END SUBROUTINE unswitch_waves_linear_int

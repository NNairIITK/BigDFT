!> @file
!!  Diagonalisation of the input guess in the gaussian basis set
!! @author
!!    Copyright (C) 2008-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Diagonalise the hamiltonian in a basis set of norbe orbitals and select the first
!!  norb eigenvectors. Works also with the spin-polarisation case and perform also the 
!!  treatment of semicore atoms. 
!!  In the absence of norbe parameters, it simply diagonalize the hamiltonian in the given
!!  orbital basis set.
!!  Works for wavefunctions given in a Gaussian basis set provided bt the structure G
subroutine Gaussian_DiagHam(iproc,nproc,natsc,nspin,orbs,G,mpirequests,&
      &   psigau,hpsigau,orbse,etol,norbsc_arr)
   use module_base
   use module_types
   use module_interfaces
  use yaml_output
   implicit none
   integer, intent(in) :: iproc,nproc,natsc,nspin
   real(gp), intent(in) :: etol
   type(gaussian_basis), intent(in) :: G
   type(orbitals_data), intent(inout) :: orbs
   type(orbitals_data), intent(in) :: orbse
   integer, dimension(nproc-1), intent(in) :: mpirequests
   integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   real(wp), dimension(orbse%nspinor*G%ncoeff,orbse%norbp), intent(in) :: psigau,hpsigau
   !local variables
   character(len=*), parameter :: subname='Gaussian_DiagHam'
   !n(c) real(kind=8), parameter :: eps_mach=1.d-12
   logical :: semicore,minimal
   integer :: i,ndim_hamovr,i_all,i_stat,norbi_max,j
   integer :: natsceff,ndh1,ispin,nspinor !n(c) norbtot,norbsc,npsidim
   real(gp) :: tolerance
   integer, dimension(:,:), allocatable :: norbgrp
   real(wp), dimension(:,:), allocatable :: hamovr

   tolerance=etol

   minimal=.true.!present(orbse)

   semicore=.true.!present(norbsc_arr)

   !define the grouping of the orbitals: for the semicore case, follow the semicore atoms,
   !otherwise use the number of orbitals, separated in the spin-polarised case
   !for the spin polarised case it is supposed that the semicore orbitals are disposed equally
   if (semicore) then
      !if (present(orbsv)) then
      !   norbi_max=max(maxval(norbsc_arr),orbsv%norb)
      !else
      norbi_max=maxval(norbsc_arr)
      !end if

      !calculate the dimension of the overlap matrix
      !take the maximum as the two spin dimensions
      ndim_hamovr=0
      do ispin=1,nspin
         ndh1=0
         !n(c) norbsc=0
         do i=1,natsc+1
            ndh1=ndh1+norbsc_arr(i,ispin)**2
         end do
         ndim_hamovr=max(ndim_hamovr,ndh1)
      end do
      if (natsc > 0) then
         if (nspin == 2) then
            if (sum(norbsc_arr(1:natsc,1)) /= sum(norbsc_arr(1:natsc,2))) then
               call yaml_warning('(Gaussian_DiagHam) The number of semicore orbitals must be the same for both spins')
               !write(*,'(1x,a)')&
               !   &   'ERROR (Gaussian_DiagHam): The number of semicore orbitals must be the same for both spins'
               stop
            end if
         end if
         !n(c) norbsc=sum(norbsc_arr(1:natsc,1))
      else
         !n(c) norbsc=0
      end if

      natsceff=natsc
      allocate(norbgrp(natsceff+1,nspin+ndebug),stat=i_stat)
      call memocc(i_stat,norbgrp,'norbgrp',subname)

      !assign the grouping of the orbitals
      do j=1,nspin
         do i=1,natsceff+1
            norbgrp(i,j)=norbsc_arr(i,j)
         end do
      end do
   else
      !this works also for non spin-polarised since there norbu=norb
      norbi_max=max(orbs%norbu,orbs%norbd) 
      ndim_hamovr=norbi_max**2

      natsceff=0
      allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
      call memocc(i_stat,norbgrp,'norbgrp',subname)

      !n(c) norbsc=0
      norbgrp(1,1)=orbs%norbu
      if (nspin == 2) norbgrp(1,2)=orbs%norbd

   end if

   !assign total orbital number for calculating the overlap matrix and diagonalise the system

   if(minimal) then
      !n(c) norbtot=orbse%norb !beware that norbe is equal both for spin up and down
      !n(c) npsidim=orbse%npsidim
      nspinor=orbse%nspinor
   else
      !n(c) norbtot=orbs%norb
      !n(c) npsidim=orbs%npsidim
      nspinor=orbs%nspinor
   end if

   allocate(hamovr(nspin*ndim_hamovr,2+ndebug),stat=i_stat)
   call memocc(i_stat,hamovr,'hamovr',subname)

   if (iproc.eq.0) call yaml_comment('Overlap Matrix...')
   !if (iproc.eq.0) write(*,'(1x,a)',advance='no') 'Overlap Matrix...'

   call overlap_and_gather(iproc,nproc,mpirequests,G%ncoeff,natsc,nspin,ndim_hamovr,orbse,&
      &   norbsc_arr,psigau,hpsigau,hamovr)

   call solve_eigensystem(norbi_max,&
      &   ndim_hamovr,sum(norbgrp),natsceff,nspin,nspinor,norbgrp,hamovr,orbs%eval)
   !!!
   !!!  !allocate the pointer for virtual orbitals
   !!!  if(present(orbsv) .and. present(psivirt) .and. orbsv%norb > 0) then
   !!!     allocate(psivirt(orbsv%npsidim+ndebug),stat=i_stat)
   !!!     call memocc(i_stat,psivirt,'psivirt',subname)
   !!!  end if
   !!!
   !!!  if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'
   !!!  nvctr=wfd%nvctr_c+7*wfd%nvctr_f
   !!!  if (.not. present(orbsv)) then
   !!!     call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
!!!          natsceff,nspin,orbs%nspinor,ndim_hamovr,norbgrp,hamovr,psi,psit)
   !!!  else
   !!!     call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
      !!!          natsceff,nspin,orbs%nspinor,ndim_hamovr,norbgrp,hamovr,psi,psit,orbsv%norb,psivirt)
!!!  end if
   !!!  
   !!!  !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
   !!!     
   i_all=-product(shape(hamovr))*kind(hamovr)
   deallocate(hamovr,stat=i_stat)
   call memocc(i_stat,i_all,'hamovr',subname)
   !!!  i_all=-product(shape(norbgrp))*kind(norbgrp)
   !!!  deallocate(norbgrp,stat=i_stat)
   !!!  call memocc(i_stat,i_all,'norbgrp',subname)
   !!!
   !!!  if (minimal) then
   !!!     !deallocate the old psi
   !!!     i_all=-product(shape(psi))*kind(psi)
   !!!     deallocate(psi,stat=i_stat)
   !!!     call memocc(i_stat,i_all,'psi',subname)
   !!!  else if (nproc == 1) then
   !!!     !reverse objects for the normal diagonalisation in serial
   !!!     !at this stage hpsi is the eigenvectors and psi is the old wavefunction
   !!!     !this will restore the correct identification
   !!!     nullify(hpsi)
   !!!     hpsi => psi
   !!!!     if(nspinor==4) call psitransspi(nvctrp,norb,psit,.false.) 
   !!!    nullify(psi)
   !!!     psi => psit
   !!!  end if
   !!!
   !!!  !orthogonalise the orbitals in the case of semi-core atoms
   !!!  if (norbsc > 0) then
   !if(nspin==1) then
   !   call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,nspinor) 
   !else
   !!!     call orthon_p(iproc,nproc,orbs%norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,&
         !!!          orbs%nspinor) 
   !!!     if(orbs%norbd > 0) then
   !!!        call orthon_p(iproc,nproc,orbs%norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,&
         !!!             psit(1+nvctrp*orbs%norbu),orbs%nspinor) 
   !   end if
   !!!     end if
   !!!  end if
   !!!
   !!!
   !!!  if (minimal) then
   !!!     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
   !!!     call memocc(i_stat,hpsi,'hpsi',subname)
   !!!!     hpsi=0.0d0
   !!!     if (nproc > 1) then
   !!!        !allocate the direct wavefunction
   !!!        allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
   !!!        call memocc(i_stat,psi,'psi',subname)
   !!!     else
   !!!        psi => psit
   !!!     end if
   !!!  end if
   !!!
   !!!  !this untranspose also the wavefunctions 
   !!!  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,&
        !!!       psit,work=hpsi,outadd=psi(1))
   !!!
   !!!  if (nproc == 1) then
   !!!     nullify(psit)
   !!!  end if

END SUBROUTINE Gaussian_DiagHam


subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,Lzde,comms,&
     psi,hpsi,psit,orthpar,passmat,iscf,Tel,occopt,& !mandatory
     orbse,commse,etol,norbsc_arr) !optional
  use module_base
  use module_types
  use module_interfaces, except_this_one => LDiagHam
  use yaml_output
  use communications_base, only: comms_cubic
  use communications, only: transpose_v, untranspose_v
  implicit none
  integer, intent(in) :: iproc,nproc,natsc,nspin,occopt,iscf
  real(gp), intent(in) :: Tel
  type(local_zone_descriptors) :: Lzd        !< Information about the locregs after LIG
  type(local_zone_descriptors) :: Lzde       !< Information about the locregs for LIG
  type(comms_cubic), intent(in) :: comms
  type(orbitals_data), intent(inout) :: orbs
  type(orthon_data), intent(in):: orthpar 
  real(wp), dimension(*), intent(out) :: passmat !< passage matrix for building the eigenvectors (the size depends of the optional arguments)
  real(wp), dimension(:), pointer :: psi,hpsi,psit
  real(gp), intent(in) :: etol
  type(orbitals_data), intent(inout) :: orbse
  type(comms_cubic), intent(in) :: commse
  integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
  !local variables
  character(len=*), parameter :: subname='LDiagHam'
  !real(kind=8), parameter :: eps_mach=1.d-12
  integer :: ikptp,ikpt,nvctrp,iorb,Gdim!,jproc
  integer :: i,ndim_hamovr,i_all,i_stat,ierr,norbi_max,j,noncoll,ispm,ncplx,idum=0
  integer :: norbtot,natsceff,norbsc,ndh1,ispin,nvctr,npsidim,nspinor,ispsi,ispsie,ispsiv
  real(kind=4) :: tt,builtin_rand
  real(gp) :: tolerance
  integer, dimension(:,:), allocatable :: norbgrp
  real(wp), dimension(:,:,:), allocatable :: hamovr
  real(wp), dimension(:), pointer :: psiw
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake
     
  tolerance=etol

  !assign total orbital number for calculating the overlap matrix and diagonalise the system
  norbtot=orbse%norb !beware that norbe is equal both for spin up and down
  npsidim=max(orbse%npsidim_orbs,orbse%npsidim_comp)
  nspinor=orbse%nspinor

  if (nproc > 1 .or. Lzde%linear) then
     Gdim = (Lzde%Glr%wfd%nvctr_c+7*Lzde%Glr%wfd%nvctr_f)*orbse%norb_par(iproc,0)*orbse%nspinor
     allocate(psiw(max(npsidim,Gdim)+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  else
     nullify(psiw)
  end if

  !transpose all the wavefunctions for having a piece of all the orbitals
  call toglobal_and_transpose(iproc,nproc,orbse,Lzde,commse,psi,psiw)
  call toglobal_and_transpose(iproc,nproc,orbse,Lzde,commse,hpsi,psiw)

  if(nproc > 1.or. Lzde%linear) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  else
     nullify(psiw)
  end if

  !define the grouping of the orbitals: for the semicore case, follow the semicore atoms,
  !otherwise use the number of orbitals, separated in the spin-polarised case
  !for the spin-polarised case it is supposed that the semicore orbitals are disposed equally
  !fon non-collinear spins, the orbitals atoms are doubled
  !calculate the maximum of the dimension for each k-point
  if (orbs%nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  norbi_max=noncoll*maxval(norbsc_arr)

  !calculate the dimension of the overlap matrix
  !take the maximum as the two spin dimensions
  ndim_hamovr=0
  do ispin=1,nspin
     ndh1=0
     norbsc=0
     do i=1,natsc+1
        ndh1=ndh1+(noncoll*norbsc_arr(i,ispin))**2
     end do
     ndim_hamovr=max(ndim_hamovr,ndh1)
  end do
  if (natsc > 0) then
     if (nspin == 2) then
        if (f_err_raise(sum(norbsc_arr(1:natsc,1)) /= sum(norbsc_arr(1:natsc,2)),&
             'ERROR (DiagHam): The number of semicore orbitals must be the same for both spins',&
             err_name='BIGDFT_RUNTIME_ERROR')) return
     end if
     norbsc=noncoll*sum(norbsc_arr(1:natsc,1))
  else
     norbsc=0
  end if
  natsceff=natsc
  
  allocate(norbgrp(natsceff+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbgrp,'norbgrp',subname)

  !assign the grouping of the orbitals
  do j=1,nspin
     do i=1,natsceff+1
        norbgrp(i,j)=noncoll*norbsc_arr(i,j)
     end do
  end do

  !for complex matrices the dimension is doubled
  if (nspinor /=1) then
     ndim_hamovr=2*ndim_hamovr
  end if

  allocate(hamovr(nspin*ndim_hamovr,2,orbse%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,hamovr,'hamovr',subname)

  !initialise hamovr
  call to_zero(nspin*ndim_hamovr*2*orbse%nkpts,hamovr(1,1,1))

  if (iproc == 0 .and. verbose > 1) call yaml_open_map('Input Guess Overlap Matrices',flow=.true.)
  !     'Overlap Matrix...'

  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted
  ispsi=1
  do ikptp=1,orbse%nkptsp
     ikpt=orbse%iskpts+ikptp!orbse%ikptsp(ikptp)
     
     nvctrp=commse%nvctr_par(iproc,ikpt)
     if (nvctrp == 0) cycle
     
     !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)
     call overlap_matrices(norbtot,nvctrp,natsceff,nspin,nspinor,&
          & ndim_hamovr,norbgrp,hamovr(1,1,ikpt),psi(ispsi),hpsi(ispsi))
     
     ispsi=ispsi+nvctrp*norbtot*orbse%nspinor
  end do

  if (iproc == 0 .and. verbose > 1) call yaml_map('Calculated',.true.)

!  if(iproc==0 .and. verbose>1) write(*,'(a)') ' done.'
  !if (iproc == 0) print *,'hamovr,iproc:',iproc,hamovr
  !deallocate hpsi in the case of a minimal basis
  call f_free_ptr(hpsi)

  if (nproc > 1) then
     !reduce the overlap matrix between all the processors
     call mpiallred(hamovr(1,1,1),2*nspin*ndim_hamovr*orbse%nkpts,&
          MPI_SUM,bigdft_mpi%mpi_comm)
  end if

! DEBUG
!  if(iproc == 0) then
!     print *,'size(hamovr)',size(hamovr,1),size(hamovr,2),size(hamovr,3)
!     do i_all=1,size(hamovr,1)
!        print *,'iel, ham, ovr:',i_all,hamovr(i_all,1,:),hamovr(i_all,2,:)
!     end do
!  end if
! END DEBUG

  !in the case of minimal basis allocate now the transposed wavefunction
  !otherwise do it only in parallel
  if(.not. associated(psit)) then
     allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
     call memocc(i_stat,psit,'psit',subname)
  end if

  ! There are two possibilities to generate the input guess
  differentInputGuess: if(.not. orthpar%directDiag) then
     if (iproc == 0 .and. verbose > 1) then
        call yaml_close_map()
        call yaml_newline()
     end if
     if(iproc==0) write(*,'(1x,a)') 'Iterative diagonalization...'

     call inputguessParallel(iproc, nproc, orbs, norbsc_arr, hamovr, &
          psi, psit, orthpar, nspin, nspinor, npsidim, comms, natsc, ndim_hamovr, norbsc)

  else
     
     !if(iproc==0) write(*,'(1x,a)') 'Direct diagonalization...'
     
     call timing(iproc, 'Input_comput', 'ON')

     !it is important that the k-points repartition of the inputguess orbitals
     !coincides with the one of the SCF orbitals
     do ikptp=1,orbse%nkptsp
        ikpt=orbse%iskpts+ikptp!orbs%ikptsp(ikptp)
        call solve_eigensystem(norbi_max,&
             ndim_hamovr,sum(norbgrp),natsceff,nspin,nspinor,norbgrp,hamovr(1,1,ikpt),&
             orbse%eval((ikpt-1)*orbse%norb+1)) !changed from orbs

        !assign the value for the orbital
        call vcopy(orbs%norbu,orbse%eval((ikpt-1)*orbse%norb+1),1,orbs%eval((ikpt-1)*orbs%norb+1),1)
        if (orbs%norbd >0) then
           call vcopy(orbs%norbd,orbse%eval((ikpt-1)*orbse%norb+orbse%norbu+1),1,orbs%eval((ikpt-1)*orbs%norb+orbs%norbu+1),1)
        end if
     end do

     if (iproc == 0 .and. verbose > 1) then
        call yaml_map('Diagonalized',.true.)
        call yaml_close_map()
        call yaml_newline()
     end if

     !broadcast values for k-points 
     call broadcast_kpt_objects(nproc, orbse%nkpts, orbse%norb, &
          & orbse%eval(1), orbse%ikptproc)

     !clean the array of the IG occupation
     call to_zero(orbse%norb*orbse%nkpts,orbse%occup(1))
        !put the actual values on it
     do ikpt=1,orbs%nkpts
        call vcopy(orbs%norbu,orbs%occup((ikpt-1)*orbs%norb+1),1,&
             orbse%occup((ikpt-1)*orbse%norb+1),1)
        if (orbs%norbd > 0) then
           call vcopy(orbs%norbd,orbs%occup((ikpt-1)*orbs%norb+orbs%norbu+1),1,&
                orbse%occup((ikpt-1)*orbse%norb+orbse%norbu+1),1)
        end if
     end do

     !here the value of the IG occupation numbers can be calculated
     if (iscf > SCF_KIND_DIRECT_MINIMIZATION .or. Tel > 0.0_gp) then

        if (iproc==0) call yaml_map('Noise added to input eigenvalues to determine occupation numbers',&
             max(Tel,1.0e-3_gp),fmt='(1pe12.5)')
!        !add a small displacement in the eigenvalues
        do iorb=1,orbse%norb*orbse%nkpts
           tt=builtin_rand(idum)
           orbse%eval(iorb)=orbse%eval(iorb)*(1.0_gp+max(Tel,1.0e-3_gp)*real(tt,gp))
        end do
        
        !correct the occupation numbers wrt fermi level
        call evaltoocc(iproc,nproc,.false.,Tel,orbse,occopt)

     end if

     if (iproc ==0) then 
        nullify(mom_vec_fake)
        call write_eigenvalues_data(tolerance,orbse,mom_vec_fake)
        !call write_ig_eigenvectors(tolerance,orbse,nspin,orbs%norb,orbs%norbu,orbs%norbd)
     end if

     !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'
     nvctr=Lzde%Glr%wfd%nvctr_c+7*Lzde%Glr%wfd%nvctr_f

     ispsi=1
     ispsie=1
     ispsiv=1
     ispm=1
     !number of complex components
     ncplx=1
     if (nspinor > 1) ncplx=2

     !@todo: see to broadcast a smaller array, if possible.
     !@todo: switch to a allgatherv to handle kpoint case.
     ! Broadcast in case of different degenerated eigen vectors.
     if (nproc > 1 .and. orbse%nkpts == 1) then
        !reduce the overlap matrix between all the processors
        call mpi_bcast(hamovr(1,1,1), 2*nspin*ndim_hamovr*orbse%nkpts, mpidtypw, &
             & 0, bigdft_mpi%mpi_comm, ierr)
     else if (nproc > 1) then
        do ikpt=1,orbse%nkpts
           if (iproc /= orbse%ikptproc(ikpt)) then
              hamovr(:,:,ikpt) = 0._gp
           end if
        end do
        call mpiallred(hamovr(1,1,1),2*nspin*ndim_hamovr*orbse%nkpts,&
             MPI_SUM,bigdft_mpi%mpi_comm)
     end if

     do ikptp=1,orbse%nkptsp
        ikpt=orbse%iskpts+ikptp!orbse%ikptsp(ikptp)

        nvctrp=commse%nvctr_par(iproc,ikpt)
        if (nvctrp == 0) cycle

        call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
             natsceff,nspin,nspinor,orbs%nspinor,ndim_hamovr,norbgrp,hamovr(1,1,ikpt),&
             psi(ispsie:),psit(ispsi:),passmat(ispm))

        ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
        ispsie=ispsie+nvctrp*norbtot*orbs%nspinor
        ispm=ispm+ncplx*(orbse%norbu*orbs%norbu+orbse%norbd*orbs%norbd)
     end do
     if (iproc == 0 .and. verbose > 1) call yaml_map('IG wavefunctions defined',.true.)
     !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
     !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

     call timing(iproc, 'Input_comput', 'OF')

  end if differentInputGuess
     
  i_all=-product(shape(hamovr))*kind(hamovr)
  deallocate(hamovr,stat=i_stat)
  call memocc(i_stat,i_all,'hamovr',subname)
  i_all=-product(shape(norbgrp))*kind(norbgrp)
  deallocate(norbgrp,stat=i_stat)
  call memocc(i_stat,i_all,'norbgrp',subname)

  !print * ,'debug2,iproc',iproc,orbsv%norb,orbsv%norbp,orbsv%norbu,orbsv%norbd,orbsv%npsidim
  !deallocate the old psi
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

  !orthogonalise the orbitals in the case of semi-core atoms
  if (norbsc > 0) then
     call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar)
  end if

  hpsi = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug,id='hpsi')
  !     hpsi=0.0d0
  if (nproc > 1) then
     !allocate the direct wavefunction
     allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)
  else
     psi => psit
  end if

  !this untranspose also the wavefunctions 
  call untranspose_v(iproc,nproc,orbs,Lzd%Glr%wfd,comms,&
       psit(1),hpsi(1),out_add=psi(1))

!!$!here the checksum of the wavefunction can be extracted
!!$do jproc=0,bigdft_mpi%nproc-1
!!$   call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$   if (jproc==bigdft_mpi%iproc) then
!!$      ispsi=1
!!$      nvctr=Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!!$      do iorb=1,orbs%norbp
!!$         write(*,'(i4,100(1pe10.2))')iorb+orbs%isorb,sum(psi(ispsi:ispsi+nvctr-1))
!!$         ispsi=ispsi+nvctr
!!$      end do
!!$   end if
!!$end do

  if (nproc == 1) nullify(psit)

  ! reput the good wavefunction dimensions:  
  if(.not. Lzd%linear) call wavefunction_dimension(Lzd,orbs)     

END SUBROUTINE LDiagHam


subroutine overlap_matrices(norbe,nvctrp,natsc,nspin,nspinor,ndim_hamovr,&
      &   norbsc_arr,hamovr,psi,hpsi)
   use module_base
   implicit none
   integer, intent(in) :: norbe,nvctrp,natsc,ndim_hamovr,nspin,nspinor
   integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   real(wp), dimension(nspin*ndim_hamovr,2), intent(out) :: hamovr
   real(wp), dimension(nvctrp*nspinor,norbe), intent(in) :: psi,hpsi
   !local variables
   integer :: iorbst,imatrst,norbi,i,ispin,ncomp,ncplx!,ierr,jproc
   !integer :: iorb,jorb
   !WARNING: here nspin=1 for nspinor=4
   if(nspinor == 1) then
      ncplx=1
      elseif(nspinor == 2) then
      ncplx=2
      ncomp=1
   else if (nspinor == 4) then
      ncplx=2
      ncomp=2
   end if

   !calculate the overlap matrix for each group of the semicore atoms
   !       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
   !       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
   iorbst=1
   imatrst=1
   do ispin=1,nspin !this construct assumes that the semicore is identical for both the spins
      do i=1,natsc+1
         norbi=norbsc_arr(i,ispin)
         if (nspinor ==1) then
            call gemm('T','N',norbi,norbi,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
               &   hpsi(1,iorbst),max(1,nvctrp),&
               &   0.0_wp,hamovr(imatrst,1),norbi)
            !here probably dsyrk can be used
            call gemm('T','N',norbi,norbi,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
               &   psi(1,iorbst),max(1,nvctrp),0.0_wp,hamovr(imatrst,2),norbi)
         else
            call c_gemm('C','N',norbi,norbi,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
               &   max(1,ncomp*nvctrp),hpsi(1,iorbst),max(1,ncomp*nvctrp),&
               &   (0.0_wp,0.0_wp),hamovr(imatrst,1),norbi)
            !here probably zherk can be used
            call c_gemm('C','N',norbi,norbi,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
               &   max(1,ncomp*nvctrp),psi(1,iorbst),max(1,ncomp*nvctrp),&
               &   (0.0_wp,0.0_wp),hamovr(imatrst,2),norbi)
         end if
!!$               open(17)
!!$               open(18)
!!$               print *,'ncplx,nspinor',ncplx,nspinor
!!$               do jorb=1,norbi
!!$                  write(17,'(i0,1x,48(1pe10.1))')jorb,&
!!$                       ((hamovr(2*(iorb-1)+icplx+(jorb-1)*ncplx*norbi,1),icplx=1,ncplx),iorb=1,8)!norbi)
!!$                  write(18,'(i0,1x,48(1pe10.1))')jorb,&
!!$                       ((hamovr(2*(iorb-1)+icplx+(jorb-1)*ncplx*norbi,2),icplx=1,ncplx),iorb=1,8)!norbi)
!!$               end do
!!$         
!!$               close(17)
!!$               close(18)
!!$                 stop
!!$if (i==natsc+1) then
!!$do jproc=0,bigdft_mpi%nproc-1
!!$   call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$   if (jproc==bigdft_mpi%iproc) then
!!$         !print out the passage matrix (valid for one k-point only and ncplx=1)
!!$         print *,'HAMILTONIAN',jproc,ispin 
!!$         do iorb=1,norbi
!!$            write(*,'(i4,100(1pe10.2))')iorb,(hamovr(imatrst+iorb-1+norbi*(jorb-1),1),jorb=1,norbi)
!!$         end do
!!$         print *,'OVERLAP',jproc,ispin 
!!$         do iorb=1,norbi
!!$            write(*,'(i4,100(1pe10.2))')iorb,(hamovr(imatrst+iorb-1+norbi*(jorb-1),2),jorb=1,norbi)
!!$       end do
!!$    end if
!!$ end do
!!$stop 
!!$end if
         iorbst=iorbst+norbi
         imatrst=imatrst+ncplx*norbi**2
      end do
   end do

END SUBROUTINE overlap_matrices


subroutine solve_eigensystem(norbi_max,ndim_hamovr,ndim_eval,&
      &   natsc,nspin,nspinor,&
      &   norbsc_arr,hamovr,eval)
   use module_base
   implicit none
   integer, intent(in) :: norbi_max,ndim_hamovr,natsc,nspin,nspinor,ndim_eval
   integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   real(wp), dimension(nspin*ndim_hamovr,2), intent(inout) :: hamovr
   real(wp), dimension(ndim_eval), intent(out) :: eval
   !local variables
   character(len=*), parameter :: subname='solve_eigensystem'
   !n(c) character(len=25) :: gapstring
   !n(c) character(len=64) :: message
   integer :: iorbst,imatrst,norbi,n_lp,info,i_all,i_stat,i,ncplx !n(c) iorb, ncomp, ndegen
   integer :: norbj,jiorb,jjorb,ihs,ispin,norbij,norbu_ig !,jproc n(c) nwrtmsg
   !n(c) real(wp), dimension(2) :: preval
   real(wp), dimension(:), allocatable :: work_lp,evale,work_rp
   !n(c) real(gp) :: HLIGgap

   !if(iproc==0) write(30100,*) hamovr(1:ndim_hamovr,1)
   !if(iproc==0) write(30110,*) hamovr(1:ndim_hamovr,2)
   !if(iproc==0) write(30101,*) hamovr(ndim_hamovr+1:2*ndim_hamovr,1)
   !if(iproc==0) write(30111,*) hamovr(ndim_hamovr+1:2*ndim_hamovr,2)

   !calculate the total number of up orbitals for the input guess
   norbu_ig=0
   do i=1,natsc+1
      norbu_ig=norbu_ig+norbsc_arr(i,1)
   end do

   !WARNING: here nspin=1 for nspinor=4
   if(nspinor == 1) then
      ncplx=1
      !ncomp=1
      elseif(nspinor == 2) then
      ncplx=2
      !ncomp=1
   else if (nspinor == 4) then
      ncplx=2
      !ncomp=2
   end if

   !find the eigenfunctions for each group
   n_lp=max(10,4*norbi_max)
   allocate(work_lp(ncplx*n_lp+ndebug),stat=i_stat)
   call memocc(i_stat,work_lp,'work_lp',subname)
   allocate(evale(nspin*norbi_max+ndebug),stat=i_stat)
   call memocc(i_stat,evale,'evale',subname)

   if (nspinor /= 1) then
      allocate(work_rp(3*norbi_max+1+ndebug),stat=i_stat)
      call memocc(i_stat,work_rp,'work_rp',subname)
   end if

   !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')'Linear Algebra...'
  
   !n(c) nwrtmsg=0
   !n(c) ndegen=0

   !n(c) preval=0.0_wp
   iorbst=1
   imatrst=1
   do i=1,natsc+1
      norbi=norbsc_arr(i,1)

      if (nspinor == 1) then
         !if(iproc==0) write(*,*) 'imatrst',imatrst
         call sygv(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
            &   norbi,evale(1),work_lp(1),n_lp,info)
         if (info /= 0) write(*,*) 'SYGV ERROR',info,i,natsc+1

         !do the diagonalisation separately in case of spin polarization     
         if (nspin==2) then
           !if(iproc==0) write(*,*) 'imatrst+ndim_hamovr',imatrst+ndim_hamovr
            norbj=norbsc_arr(i,2)
            call sygv(1,'V','U',norbj,hamovr(imatrst+ndim_hamovr,1),&
               &   norbj,hamovr(imatrst+ndim_hamovr,2),norbj,evale(norbi+1),work_lp(1),n_lp,info)
            if (info /= 0) write(*,*) 'SYGV ERROR',info,i,natsc+1
         end if
      else
         call hegv(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
            &   norbi,evale(1),work_lp(1),n_lp,work_rp(1),info)
         if (info /= 0) write(*,*) 'HEGV ERROR',info,i,natsc+1

         !do the diagonalisation separately in case of spin polarization     
         if (nspin==2) then
            norbj=norbsc_arr(i,2)
            call hegv(1,'V','U',norbj,hamovr(imatrst+ndim_hamovr,1),&
               &   norbj,hamovr(imatrst+ndim_hamovr,2),norbj,evale(norbi+1),&
               &   work_lp(1),n_lp,work_rp(1),info)
            if (info /= 0) write(*,*) 'HEGV ERROR',info,i,natsc+1
         end if

      end if

      !check the sign of the eigenvector, control the choice per MPI process
      !this is useful when different MPI processes gave eigenvectors with different phases
      norbij=norbi
      ihs=imatrst
      do ispin=1,nspin
         do jjorb=1,norbij
            !if it is negative change the sign to all the values
            if (hamovr(ihs+(jjorb-1)*norbij*ncplx,1) < 0.0_wp) then
               do jiorb=1,norbij*ncplx
                  hamovr(ihs-1+jiorb+(jjorb-1)*norbij*ncplx,1)=&
                      -hamovr(ihs-1+jiorb+(jjorb-1)*norbij*ncplx,1)
               end do
            end if
         end do
         if (nspin==2) then
            norbij=norbj
            ihs=ihs+ndim_hamovr
         end if
      end do


      !!$     if (iproc == 0) then
!!$     print *,norbi,ncomp,ncplx,imatrst
!!$     !write the matrices on a file
!!$     !open(12)
!!$do jproc=0,bigdft_mpi%nproc-1
!!$call MPI_BARRIER(bigdft_mpi%mpi_comm,i_stat)
!!$if (jproc==bigdft_mpi%iproc) then
!!$  print *,'PASSAGE MATRIX',jproc
!!$        open(33+2*(i-1)+100*bigdft_mpi%iproc)
!!$     do jjorb=1,norbi
        !   do jiorb=1,norbi
        !      write(12,'(1x,2(i0,1x),200(1pe24.17,1x))')jjorb,jiorb,&
        !           hamovr(jjorb+norbi*(jiorb-1),1),hamovr(jjorb+norbi*(jiorb-1),2)
        !   end do
        !end do
        !close(12)
!!$        write(33+2*(i-1)+100*bigdft_mpi%iproc,'(2000(1pe10.2))')&
             !(hamovr(imatrst-1+jiorb+(jjorb-1)*norbi*ncomp*ncplx,1),jiorb=1,8*ncomp*ncplx)
!!$             (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi*ncomp*ncplx,1),jiorb=1,norbi*ncomp*ncplx)
!!$        write(*,'(1x,2(i6),2000(1pe10.2))')jjorb,jiorb,(hamovr(jjorb+norbi*(jiorb-1),1),jiorb=1,norbi)
!!$
!!$     end do
!!$  end if
!!$end do
!!$     close(33+2*(i-1)+100*bigdft_mpi%iproc)
!!$     open(34+2*(i-1)+100*iproc)
!!$     do jjorb=1,8!norbi
!!$        write(34+2*(i-1)+100*iproc,'(2000(1pe10.2))')&
            !!$                (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi*ncomp*ncplx,2),jiorb=1,8*ncomp*ncplx)
!!$        !                (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi*ncomp*ncplx,2),jiorb=1,norbi*ncomp*ncplx)
!!$     end do
!!$     close(34+2*(i-1)+100*iproc+100*iproc)
      !!$
!!$     !end if
      !!$     stop

      !!$     if (iproc ==0) then
      !!$        call write_ig_eigenvectors(etol,norbi,nspin,iorbst,norb,norbu,norbd,evale)
      !!$     end if
      !!$     if (nspin==1) then
      !!$        do iorb=iorbst,min(norbi+iorbst-1,norb)
      !!$           eval(iorb)=evale(iorb-iorbst+1)
      !!$        end do
      !!$     else
      !!$        do iorb=iorbst,min(norbi+iorbst-1,norbu)
      !!$           eval(iorb)=evale(iorb-iorbst+1)
      !!$        end do
      !!$        do iorb=iorbst,min(norbi+iorbst-1,norbd)
      !!$           eval(iorb+norbu)=evale(iorb-iorbst+1+norbi)
      !!$        end do
      !!$     end if

      !here we should copy all the eigenvalues in the eval array
      if (nspin > 2) stop 'ERROR(solve_eigensystem): nspin too high'

      call vcopy(norbi,evale(1),1,eval(iorbst),1)
      if (nspin==2) call vcopy(norbi,evale(norbi+1),1,eval(iorbst+norbu_ig),1)

      iorbst=iorbst+norbi
      imatrst=imatrst+ncplx*norbi**2
   end do

   if (nspinor /= 1) then
      i_all=-product(shape(work_rp))*kind(work_rp)
      deallocate(work_rp,stat=i_stat)
      call memocc(i_stat,i_all,'work_rp',subname)
   end if

   i_all=-product(shape(work_lp))*kind(work_lp)
   deallocate(work_lp,stat=i_stat)
   call memocc(i_stat,i_all,'work_lp',subname)
   i_all=-product(shape(evale))*kind(evale)
   deallocate(evale,stat=i_stat)
   call memocc(i_stat,i_all,'evale',subname)

END SUBROUTINE solve_eigensystem


subroutine build_eigenvectors(norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinore,nspinor,&
      &   ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,passmat)
   use module_base
   implicit none
   !Arguments
   integer, intent(in) :: norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinor,ndim_hamovr,nspinore
   integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   real(wp), dimension(nspin*ndim_hamovr), intent(in) :: hamovr
   real(wp), dimension(nvctrp*nspinore,norbe), intent(in) :: psi
   real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: ppsit
   real(wp), dimension(*), intent(out) :: passmat !< passage matrix between ppsit and psi (the size depends of the complex arguments)
   !Local variables
   !n(c) character(len=*), parameter :: subname='build_eigenvectors'
   !n(c) integer, parameter :: iunit=1978
   integer :: ispin,iorbst,iorbst2,imatrst,norbsc,norbi,norbj
   integer :: ncplx,ncomp,i,ispsiv
   integer:: ispm

!  if(iproc==0) then
!      do j=1,size(hamovr)
!          !write(100001,*) hamovr(j)
!      end do
!  end if

   !WARNING: here nspin=1 for nspinor=4
   if(nspinor == 1) then
      ncplx=1
      ncomp=1
      elseif(nspinor == 2) then
      ncplx=2
      ncomp=1
   else if (nspinor == 4) then
      ncplx=2
      ncomp=2
   end if

   !perform the vector-matrix multiplication for building the input wavefunctions
   ! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
   !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
   !allocate the pointer for virtual orbitals
   iorbst=1
   iorbst2=1
   imatrst=1
   ispsiv=1
   ispm=1
   do ispin=1,nspin
      norbsc=0
      do i=1,natsc
         norbi=norbsc_arr(i,ispin)
         norbsc=norbsc+norbi
         if (nspinor == 1) then
            call gemm('N','N',nvctrp,norbi,norbi,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
               &   hamovr(imatrst),norbi,0.0_wp,ppsit(1,iorbst2),max(1,nvctrp))
         else
            call c_gemm('N','N',ncomp*nvctrp,norbi,norbi,(1.0_wp,0.0_wp),&
               &   psi(1,iorbst),max(1,ncomp*nvctrp),hamovr(imatrst),norbi,&
               &   (0.0_wp,0.0_wp),ppsit(1,iorbst2),max(1,ncomp*nvctrp))
         end if
         !store the values of the passage matrix in passmat array
         call vcopy(ncplx*norbi**2,hamovr(imatrst),1,passmat(imatrst),1)

         iorbst=iorbst+norbi
         iorbst2=iorbst2+norbi
         imatrst=imatrst+ncplx*norbi**2
         ispm=ispm+ncplx*norbi**2
      end do
      norbi=norbsc_arr(natsc+1,ispin)
      if(ispin==1) norbj=norbu-norbsc
      if(ispin==2) norbj=norbd-norbsc
      !        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
      !        norbj=norb-norbsc
      if(norbj>0) then
         if (nspinor == 1) then
            call gemm('N','N',nvctrp,norbj,norbi,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
               &   hamovr(imatrst),norbi,0.0_wp,ppsit(1,iorbst2),max(1,nvctrp))
         else
            call c_gemm('N','N',ncomp*nvctrp,norbj,norbi,(1.0_wp,0.0_wp),&
               &   psi(1,iorbst),max(1,ncomp*nvctrp),hamovr(imatrst),norbi,&
               &   (0.0_wp,0.0_wp),ppsit(1,iorbst2),max(1,ncomp*nvctrp))

         end if
      end if
      !store the values of the passage matrix in passmat array
      !print *,'iproc,BBBB11111,',iproc,ispm,norbi,norbj
      call vcopy(ncplx*norbi*norbj,hamovr(imatrst),1,passmat(ispm),1)
      ispm=ispm+ncplx*norbi*norbj
      !print *,'iproc,BBBB,',iproc,ispm,norbi,norbj

      iorbst=norbi+norbsc+1 !this is equal to norbe+1
      iorbst2=norbu+1
      imatrst=ndim_hamovr+1
   end do

END SUBROUTINE build_eigenvectors


!> Reads magnetic moments from file ('moments') and transforms the
!! atomic orbitals to spinors 
!! @warning Does currently not work for mx<0
subroutine psitospi(iproc,nproc,norbe,norbep, &
      &   nvctr_c,nvctr_f,nat,nspin,spinsgne,otoa,psi)
   use module_base
   use yaml_output
   implicit none
   !Arguments
   integer, intent(in) :: norbe,norbep,iproc,nproc,nat,nspin
   integer, intent(in) :: nvctr_c,nvctr_f
   integer, dimension(norbep), intent(in) :: otoa
   integer, dimension(norbe*nspin), intent(in) :: spinsgne
   real(kind=8), dimension(nvctr_c+7*nvctr_f,4*norbep), intent(out) :: psi
   !local variables
   character(len=*), parameter :: subname='psitospi'
   logical :: myorbital
   integer :: i_all,i_stat,nvctr
   integer :: iorb,jorb,iat,i
   real(kind=8) :: mx,my,mz,mnorm,fac
   real(kind=8), dimension(:,:), allocatable :: mom
   !n(c) integer, dimension(2) :: iorbsc,iorbv

   !initialise the orbital counters
   !n(c) iorbsc(1)=0
   !n(c) iorbv(1)=norbsc
   !used in case of spin-polarisation, ignored otherwise
   !n(c) iorbsc(2)=norbe
   !n(c) iorbv(2)=norbsc+norbe

   !if (iproc ==0) write(*,'(1x,a)',advance='no')'Transforming AIO to spinors...'
   if (iproc ==0) call yaml_map('Transforming AIO to spinors',.true.)

   nvctr=nvctr_c+7*nvctr_f

   allocate(mom(3,nat+ndebug),stat=i_stat)
   call memocc(i_stat,mom,'mom',subname)

   open(unit=1978,file='moments')
   do iat=1,nat
      read(1978,*) mx,my,mz
      mnorm=sqrt(mx**2+my**2+mz**2)
      mom(1,iat)=mx/mnorm
      mom(2,iat)=my/mnorm
      mom(3,iat)=mz/mnorm
   end do
   close(1978)
   fac=0.5d0
   do iorb=norbep*nproc,1,-1
      jorb=iorb-iproc*norbep
      !     print *,'Kolla', shape(psi),4*iorb,shape(spinsgne),iorb
      if (myorbital(iorb,nspin*norbe,iproc,nproc)) then
         mx=mom(1,otoa(iorb))
         my=mom(2,otoa(iorb))
         mz=mom(3,otoa(iorb))
         if(spinsgne(jorb)>0.0d0) then
            do i=1,nvctr
               psi(i,iorb*4-3) = (mz+fac*(my+mx))*psi(i,iorb)
               psi(i,iorb*4-2) = fac*(my-mx)*psi(i,iorb)
               psi(i,iorb*4-1) = (fac*(mx-my))*psi(i,iorb)
               psi(i,iorb*4)   = fac*(my-mx)*psi(i,iorb)
            end do
         else
            do i=1,nvctr
               psi(i,iorb*4-3) = (fac*(mx+my))*psi(i,iorb)
               psi(i,iorb*4-2) = -fac*(my+mx)*psi(i,iorb)
               psi(i,iorb*4-1) = -(mz+fac*(my+mx))*psi(i,iorb)
               psi(i,iorb*4)   = -fac*(my-mx)*psi(i,iorb)
            end do
         end if
      end if
      !     print *,'OtoA',(otoa(iorb),iorb=1,norbe)

   end do
   i_all=-product(shape(mom))*kind(mom)
   deallocate(mom,stat=i_stat)
   call memocc(i_stat,i_all,'mom',subname)

   !if (iproc ==0) write(*,'(1x,a)')'done.'

END SUBROUTINE psitospi


!> Generates an input guess for the wavefunctions. 
!! To do this, the eigenvectors of the Hamiltonian are found by an iterative procedure.
!! This gives a guess for the orbitals in the basis of atomic orbitals. These eigenfunctions are then transformed to the
!! wavelet basis.
subroutine inputguessParallel(iproc, nproc, orbs, norbscArr, hamovr, psi,&
      &   psiGuessWavelet, orthpar, nspin, nspinor, sizePsi, comms, natsc, ndim_hamovr, norbsc)
   use module_base
   use module_types
   use yaml_output
   use communications_base, only: comms_cubic
   implicit none

   ! Calling arguments
   integer, intent(in) ::  iproc           !< Process ID
   integer, intent(in) ::  nproc           !< Number of processes
   integer, intent(in) ::  nspin           !< Nspin==1 -> no spin polarization, Nspin==2 -> spin polarization
   integer, intent(in) ::  nspinor         !< Nspinor==1 -> real wavefunction,  Nspinor==2 -> complex wavefunction
   integer, intent(in) ::  sizePsi         !< Length of the vector psi
   integer, intent(in) ::  natsc           !< Number of semicore atoms
   integer, intent(in) ::  ndim_hamovr     !< First dimension of hamovr
   integer, intent(in) ::  norbsc          !< Number of semicore orbitals
   integer, dimension(natsc+1,nspin), intent(in) :: norbscArr
   type(orbitals_data), intent(inout) :: orbs
   !> Array containing both Hamiltonian and overlap matrix:
   !! hamovr(:,:,1,:) is the Hamiltonian, hamovr(:,:,2,:) is the overlap matrix
   real(kind=8), dimension(ndim_hamovr,nspin,2,orbs%nkpts), intent(inout):: hamovr
   real(kind=8), dimension(sizePsi),intent(in):: psi !< Contains the atomic orbitals
   real(kind=8), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(out):: psiGuessWavelet !< Contains the input guess vectors in wavelet basis
   type(orthon_data),intent(in) :: orthpar
   type(comms_cubic), intent(in):: comms

   ! Local variables
   integer :: i, j, iorb, jorb, ispin, ii, jj, kk, norbtot, norbtotPad, iter, ierr, itermax
   integer :: ist, i_stat, i_all, nprocSub, ishift, jjorb
   integer :: jspin, ist2, ishift2, norb, blocksize, lwork
   integer :: norbi, norbj, iat, imatrst, imatrst2, ihs, jiorb, norbij, ncplx, ncomp
   integer :: istpsi, istpsit, nvctrp, kkSave, iiSave, jproc, ist3, ikptp2, ikpt2
   integer :: blocksizeSmall, nprocSubu, nprocSubd, nspinSub, isp, kp, ikpt, ikptp
   integer :: istpsiS, istpsitS, info, wholeGroup, newGroup, newComm
   integer :: newGroupu, newCommu, newGroupd, newCommd
   integer,dimension(:),allocatable:: sendcounts, sdispls, recvcounts, rdispls, norbpArr,&
      &   norbtotpArr, kArr, kstArr, nkArr, norbArr
   integer,dimension(:),allocatable:: newID, newIDu, newIDd, norbpArrSimul, norbpArrSimulLoc
   real :: ttreal
   real(kind=8):: gradientNorm, cosangle, tol, tt, ddot, dnrm2
   real(kind=8),dimension(2):: gradientMax
   real(kind=8),dimension(:),allocatable:: alphaArr, rayleigh, evale, sceval, work, rwork
   real(kind=8),dimension(:,:),allocatable:: gradient, gradientOld
   real(kind=8),dimension(:,:,:),allocatable:: psiGuessP, overlapPsiGuessP, psiGuess, psiGuessPTrunc, sortArr
   real(kind=8),dimension(:,:,:,:),allocatable:: HamPad, overlapPad
   integer,dimension(:,:,:),allocatable:: kpArr
   logical:: success, warning, simul
   complex(kind=8):: zdotc, zz
  integer :: stat(mpi_status_size)
  character(len=*),parameter :: subname='inputguessParallel'

   ! Start the timing for the input guess.
   call timing(iproc, 'Input_comput', 'ON')

   ! First determine how many processes are necessary to perform the input guess. The number of orbitals that
   ! is to be handled by a process is specified by input%norbpInguess.
   ! Processes 1 to nprocSubu will handle the up orbitals, processes nprocSubu+1 to nprocSubu+nprocSubd will
   ! handle the down orbitals. If npsin==1, set nprocSubu=nprocSubd=nprocSub.
   if(nspin==2) then
      nprocSubu=ceiling(dble((orbs%norbu-norbsc)*orbs%nkpts)/dble(orthpar%norbpInguess))
      nprocSubd=ceiling(dble((orbs%norbd-norbsc)*orbs%nkpts)/dble(orthpar%norbpInguess))
      if(nprocSubu+nprocSubd<=nproc) then
         ! There are enough processes to treat spin up and down simultaneously. This is indicated by setting
         ! sumul=.true. and nspinSub=1, i.e. there will be no spin loop ispin=1,nspinSub.
         if(iproc==0) write(*,'(3x,a)') 'There are enough processors to treat spin up and down simultaneously.'
         simul=.true.
         nspinSub=1
      else
         ! Otherwise simul=.false. and there is a loop ispin=1,nspinSub.
         simul=.false.
         nspinSub=nspin
      end if
   else if(nspin==1) then
      nspinSub=nspin
      nprocSub=ceiling(dble((orbs%norb-norbsc)*orbs%nkpts)/dble(orthpar%norbpInguess))
      nprocSubu=nprocSub
      nprocSubd=nprocSub
      simul=.false.
   end if


   allocate(norbArr(nspin), stat=i_stat)
   call memocc(i_stat, norbArr, 'norbArr', subname)
   do ispin=1,nspin
      if(ispin==1) norbArr(ispin)=(orbs%norbu-norbsc)*orbs%nkpts
      if(ispin==2) norbArr(ispin)=(orbs%norbd-norbsc)*orbs%nkpts
   end do


   ! Now make a loop over spin up and down. If simul is true (i.e. nSpinSub=1), there is actually no loop.
   spinLoop: do isp=1,nspinSub
      ! Determine the value of ispin. If simul is false, then simply ispin=isp.
      if(.not.simul) then
         ispin=isp
      else
         ! Otherwise ispin is 1 for the processes treating the up orbitals and
         ! 2 for the processes treating the down orbitals.
         if(0<=iproc .and. iproc<=nprocSubu-1) then
            ispin=1
         else
            ispin=2
         end if
      end if
      if(iproc==0 .and. .not.simul) then
         if(ispin==1) write(*,'(3x,a)') 'treating up orbitals.'
         if(ispin==2) write(*,'(3x,a)') 'treating down orbitals.'
      end if

      ! Determine the value of norb, again depending wheter the process is treating up or
      ! down orbitals.
      if(ispin==1) then
         norb=(orbs%norbu-norbsc)*orbs%nkpts
      else if(ispin==2) then
         norb=(orbs%norbd-norbsc)*orbs%nkpts
      end if


      ! The variable input%norbpInguess indicates how many orbitals shall be treated by each process.
      ! This requires at least norb/input%norbpInguess processes.
      if(orthpar%norbpInguess>norb) then
         call yaml_warning('You want each process to treat '//trim(yaml_toa(orthpar%norbpInguess))//&
              ' orbitals, whereas there are only '//trim(yaml_toa(norb)))
         !if(iproc==0) write(*,'(5x,a,2(i0,a))') 'WARNING: You want each process to treat ',&
         !     &   orthpar%norbpInguess,' orbitals, whereas there are only ',norb,'.'
        ! if(iproc==0) write(*,'(7x,a,i0,a)') 'The value of orthpar%norbpInguess is adjusted to ',norb, '.'
         !orthpar%norbpInguess=norb
      end if

      ! Determine how many processes are necessary for the input guess. If this number exceeds the total
      ! number of processes, the program does not stop, but adjusts the number of orbitals treated by each
      ! process and prints a warning.
      nprocSub=ceiling(dble(norb)/dble(min(orthpar%norbpInguess,norb)))
      if(nprocSub>nproc) then
         warning=.true.
         nprocSub=nproc
      else
         warning=.false.
      end if

      ! Determine how the orbitals shall be distributed to these processes:
      !   norbpArr(i) (i=0,nproc-1) is the number of orbitals treated by process i if the wavefunction is not transposed.
      if(ispin==1 .or. simul) then
         allocate(norbpArr(0:nproc-1), stat=i_stat)
         call memocc(i_stat, norbpArr, 'norbpArr', subname)
      end if
      norbpArr=0
      tt=dble(norb)/dble(nprocSub)
      ii=floor(tt)
      ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
      if(.not.simul .or. iproc<nprocSubu) then
         ! These are the up orbitals if simul is true. Otherwise theses are the up orbitals 
         ! in the first iteration of spinLoop and the down orbitals in the second iteration.
         norbpArr(0:nprocSub-1)=ii
         kk=norb-nprocSub*ii
         norbpArr(0:kk-1)=ii+1
      else
         ! These are the down orbitals if simul is true. Otherwise this branch is never reached.
         norbpArr(nprocSubu:nprocSubu+nprocSubd-1)=ii
         kk=norb-nprocSub*ii
         norbpArr(nprocSubu:nprocSubu+kk-1)=ii+1
      end if
      ! Save kk and ii for printing them later in.
      kkSave=kk 
      iiSave=ii

      ! If simul is true, then there are two version of norbpArr, one for the up orbitals (and known only by
      ! the processes handling the up orbitals) and one for the down orbitals (and known only by the processes
      ! handling the down orbitals). In the end it is however necessary that all processes know norbpArr of all 
      ! processes. This is achieved by merging the two versions of norbpArr in norbpArrSimul and distributing it
      ! to all processes.
      if(simul) then
         allocate(norbpArrSimul(0:nproc-1), stat=i_stat)
         call memocc(i_stat, norbpArrSimul, 'norbpArrSimul', subname)
         allocate(norbpArrSimulLoc(0:nproc-1), stat=i_stat)
         call memocc(i_stat, norbpArrSimulLoc, 'norbpArrSimulLoc', subname)
         norbpArrSimul=0
         norbpArrSimulLoc=0
         if(iproc<nprocSubu+nprocSubd) norbpArrSimulLoc(iproc)=norbpArr(iproc)
         call mpi_allreduce(norbpArrSimulLoc(0), norbpArrSimul(0), nprocSubu+nprocSubd,&
              mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
         i_all=-product(shape(norbpArrSimulLoc))*kind(norbpArrSimulLoc)
         deallocate(norbpArrSimulLoc, stat=i_stat)
         call memocc(i_stat, i_all, 'norbpArrSimulLoc', subname)
      end if

      ! Determine which orbitals belong to which k-point. kArr(i)=j means that orbital i belongs to k-point j.
      ! Since the k-points are distributed among several processes, the values of kArr are of course distinct for
      ! each process.
      if(.not. simul .and. ispin==2) then
         i_all=-product(shape(kArr))*kind(kArr)
         deallocate(kArr, stat=i_stat)
         call memocc(i_stat, i_all, 'kArr', subname)
      end if
      allocate(kArr(1:norbpArr(iproc)), stat=i_stat)
      call memocc(i_stat, kArr, 'kArr', subname)
      ii=0
      if(.not.simul) then
         jj=nprocSub-1
      else
         jj=nprocSubu+nprocSubd-1
      end if
      ! First make a loop over all processes involved in calculating the input guess.
      do i=0,jj
         ! Now make a loop over the orbitals handled by process i. Together with the outer loop, this will
         ! create a loop over all orbitals of all k-points.
         do iorb=1,norbpArr(i)
            ! kk is now going along all orbitals of all k-points. If iproc==i, then process iproc
            ! is handling this orbital and kk is therefore written to kArr.
            if(ispin==1) kk=floor(real(ii)/real(orbs%norbu-norbsc))+1
            if(ispin==2) kk=floor(real(ii)/real(orbs%norbd-norbsc))+1
            if(iproc==i) kArr(iorb)=kk
            ii=ii+1
         end do
      end do

      ! kp is the number of k-points handled by the process.
      kp=maxval(kArr)-minval(kArr)+1

      ! Determine the starting index of the k-points within kArr (kstArr) and the number of orbitals
      ! belonging to a given k-point and handled by a given process (nkArr):
      !   kstArr(i)=j means that the orbitals belonging to the i-th k-point handled by this process
      !     start at index j in kArr. kstArr(1)=1 is always fulfilled.
      !   nkArr(i)=j means that this process handels j orbitals belonging to the i-th k-point handled
      !     by this process.
      if(.not. simul .and. ispin==2) then
         i_all=-product(shape(kstArr))*kind(kstArr)
         deallocate(kstArr, stat=i_stat)
         call memocc(i_stat, i_all, 'kstArr', subname)
      end if
      allocate(kstArr(kp), stat=i_stat)
      call memocc(i_stat, kstArr, 'kstArr', subname)
      if(.not. simul .and. ispin==2) then
         i_all=-product(shape(nkArr))*kind(nkArr)
         deallocate(nkArr, stat=i_stat)
         call memocc(i_stat, i_all, 'nkArr', subname)
      end if
      allocate(nkArr(kp), stat=i_stat)
      call memocc(i_stat, nkArr, 'nkArr', subname)
      ! The first starting index is of course one.
      kstArr(1)=1
      if(kp>1) then
         ! This means that the process handles more than one k-point.
         ii=1
         jj=1
         ! Make a loop over all orbitals handled by this process.
         do iorb=2,norbpArr(iproc)
            if(kArr(iorb)>kArr(iorb-1)) then
               ! We reached a 'boundary' of two k-points.
               ii=ii+1
               kstArr(ii)=iorb
               nkArr(ii-1)=jj
               jj=0
            end if
            ! Count up the orbitals
            jj=jj+1
         end do
         ! These are the remaining orbitals.
         nkArr(ii)=jj
      else
         ! If the process handles only one k-point, then obviously all orbitals
         ! belong to this k-point.
         nkArr(1)=norbpArr(iproc)
      end if

      ! Write some information on how the input guess is generated.    
      if(.not. warning) then
         if(.not.simul) then
            if(iproc==0) write(*,'(5x,a,i0,a)') 'Generating input guess with ',nprocSub,' processes:'
            if(iproc==0 .and. kkSave/=0) write(*,'(7x,a,5(i0,a))') 'Processes from 0 to ',kkSave-1,' treat ',&
               &   iiSave+1,' orbitals, processes from ',kkSave,' to ',nprocSub-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. kkSave==0) write(*,'(7x,a,2(i0,a))') 'Processes from 0 to ',nprocSub-1,' treat ',&
               &   iiSave,' orbitals.'
            if(iproc==0 .and. nprocSub<nproc) write(*,'(7x,a,2(i0,a))') 'Processes from ',nprocSub,' to ',&
               &   nproc,' will idle.'
         else
            if(iproc==0) write(*,'(3x,a,i0,a)') 'Generating input guess with ',nprocSubu+nprocSubd,' processes:'
            if(iproc==0 .and. kkSave/=0) write(*,'(5x,a,5(i0,a))') 'up orbitals: Processes from 0 to ',kkSave-1,&
               &   ' treat ',iiSave+1,' orbitals, processes from ',&
               &   kkSave,' to ',nprocSub-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. kkSave==0) write(*,'(5x,a,2(i0,a))') 'up orbitals: Processes from 0 to ',&
               &   nprocSub-1,' treat ',iiSave,' orbitals.'

            ! Send some information to root process to ensure that they are printed in the correct order.
            if(iproc==nprocSubu) call mpi_send(kkSave, 1, mpi_integer, 0, 1, bigdft_mpi%mpi_comm, ierr)
            if(iproc==nprocSubu) call mpi_send(iiSave, 1, mpi_integer, 0, 2, bigdft_mpi%mpi_comm, ierr)
            if(iproc==0) call mpi_recv(kkSave, 1, mpi_integer, nprocSubu, 1, bigdft_mpi%mpi_comm, stat, ierr)
            if(iproc==0) call mpi_recv(iiSave, 1, mpi_integer, nprocSubu, 2, bigdft_mpi%mpi_comm, stat, ierr)
            if(iproc==0 .and. kkSave/=0) write(*,'(5x,a,6(i0,a))') 'down orbitals: Processes from ',&
               &   nprocSubu,' to ',nprocSubu+kkSave-1,' treat ',iiSave+1,' orbitals, processes from ',&
               &   nprocSubu+kkSave,' to ',nprocSubu+nprocSubd-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. kkSave==0) write(*,'(5x,a,3(i0,a))') 'down orbitals: Processes from ',&
               &   nprocSubu,' to ',nprocSubu+nprocSubd-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. nprocSub<nproc) write(*,'(5x,a,2(i0,a))') 'Processes from ',nprocSubu+nprocSubd,&
               &   ' to ',nproc,' will idle.'
         end if
      else
         if(.not. simul) then
            if(iproc==0) write(*,'(5x,a,i0,a)') 'Generating input guess with ',nprocSub,' processes:'
            if(iproc==0 .and. kkSave/=0) write(*,'(7x,a,5(i0,a))') 'Processes from 0 to ',kkSave-1,'  treat ',iiSave+1,&
               &   ' orbitals, processes from ',kkSave,' to ',nprocSub-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. kkSave==0) write(*,'(7x,a,2(i0,a))') 'Processes from 0 to ',nprocSub-1,'  treat ',iiSave,' orbitals.'
            if(iproc==0) write(*,'(7x,a)') 'WARNING: this is more than intended!'
         else
            if(iproc==0) write(*,'(3x,a,i0,a)') 'Generating input guess with ',nprocSub,' processes:'
            if(iproc==0 .and. kkSave/=0) write(*,'(5x,a,5(i0,a))') 'Processes from 0 to ',kkSave-1,'  treat ',iiSave+1,&
               &   ' orbitals, processes from ',kkSave,' to ',nprocSub-1,' treat ',iiSave,' orbitals.'
            if(iproc==0 .and. kkSave==0) write(*,'(5x,a,2(i0,a))') 'Processes from 0 to ',nprocSub-1,'  treat ',iiSave,' orbitals.'
            if(iproc==0) write(*,'(5x,a)') 'WARNING: this is more than intended!'
         end if
      end if

      if(nprocSub>norb) then
         if(iproc==0) write(*,'(3x,a)') 'ERROR: the number of MPI processes must not be larger than the number of orbitals.'
         if(iproc==0) write(*,'(3x,a,i0,a)') 'You cannot use more than ',norb,' processes.'
         stop
      end if


      ! Since the input guess is possibly generated with only a part of all MPI processes, we need a new MPI communicator
      ! involving only the active processes. For the case where simul is true, we even need two new MPI communicators
      ! handling up and down orbitals, respectively.
      if(ispin==2 .and. .not.simul) then
         i_all=-product(shape(newID))*kind(newID)
         deallocate(newID, stat=i_stat)
         call memocc(i_stat, i_all, 'newID', subname)
      end if
      if(.not.simul) then
         allocate(newID(0:nprocSub-1), stat=i_stat)
         call memocc(i_stat, newID, 'newID', subname)
         ! Assign the IDs of the active processes to newID.
         do iorb=0,nprocSub-1
            newID(iorb)=iorb
         end do
         ! Create the new communicator newComm.
         call mpi_comm_group(bigdft_mpi%mpi_comm, wholeGroup, ierr)
         call mpi_group_incl(wholeGroup, nprocSub, newID, newGroup, ierr)
         call mpi_comm_create(bigdft_mpi%mpi_comm, newGroup, newComm, ierr)
      else
         allocate(newIDu(0:nprocSubu-1), stat=i_stat)
         call memocc(i_stat, newIDu, 'newIDu', subname)
         allocate(newIDd(0:nprocSubd-1), stat=i_stat)
         call memocc(i_stat, newIDd, 'newIDd', subname)
         ! Assign the IDs of the processes handling the up orbitals to newIDu
         do iorb=0,nprocSubu-1
            newIDu(iorb)=iorb
         end do
         ! Create the new communicator newCommu.
         call mpi_comm_group(bigdft_mpi%mpi_comm, wholeGroup, ierr)
         call mpi_group_incl(wholeGroup, nprocSubu, newIDu, newGroupu, ierr)
         call mpi_comm_create(bigdft_mpi%mpi_comm, newGroupu, newCommu, ierr)
         ! Assign the IDs of the processes handling the down orbitals to newIDd
         do iorb=0,nprocSubd-1
            newIDd(iorb)=nprocSubu+iorb
         end do
         ! Create the new communicator newCommd.
         call mpi_comm_group(bigdft_mpi%mpi_comm, wholeGroup, ierr)
         call mpi_group_incl(wholeGroup, nprocSubd, newIDd, newGroupd, ierr)
         call mpi_comm_create(bigdft_mpi%mpi_comm, newGroupd, newCommd, ierr)
      end if


      ! In order to symplify the transposing/untransposing, the orbitals are padded with zeros such that 
      ! they can be distributed evenly over all processes when being transposed. The new length of the 
      ! orbitals after this padding is then given by norbtotPad.
      norbtot=norbscArr(natsc+1,ispin)
      norbtotPad=norbtot
      do
         if(mod(norbtotPad, nprocSub)==0) exit
         norbtotPad=norbtotPad+1
      end do
      if(.not. simul) then
         if(iproc==0 .and. ispin==1) write(*,'(5x,a,i0,a)') &
            &   'up orbitals: padding the orbitals with ',norbtotPad-norbtot,' zeros.'
      else
         if(iproc==0 .and. ispin==1) write(*,'(3x,a,i0,a)') &
            &   'up orbitals: padding the orbitals with ',norbtotPad-norbtot,' zeros.'
      end if
      if(simul) then
         ! Send some information to root process to ensure that they are printed in the correct order.
         if(iproc==nprocSubu) call mpi_send(norbtotPad, 1, mpi_integer, 0, 1, bigdft_mpi%mpi_comm, ierr)
         if(iproc==0) call mpi_recv(norbtotPad, 1, mpi_integer, nprocSubu, 1, bigdft_mpi%mpi_comm, stat, ierr)
      end if
      if(.not. simul) then
         if((simul .and. iproc==0) .or. (.not. simul .and. ispin==2 .and. iproc==0))&
            &   write(*,'(5x,a,i0,a)') &
            &   'down orbitals: padding the orbitals with ',norbtotPad-norbtot,' zeros.'
      else
         if((simul .and. iproc==0) .or. (.not. simul .and. ispin==2 .and. iproc==0))&
            &   write(*,'(3x,a,i0,a)') &
            &   'down orbitals: padding the orbitals with ',norbtotPad-norbtot,' zeros.'
      end if


      ! Allocate the reamaining arrays.
      ! Allocate them also for the processes which do not treat any orbital; in this case, allocate
      ! them with (norbtotPad,1).
      if(ispin==2 .and. .not.simul) then
         ! Deallocate all arrays to reallocate them with different size
         i_all=-product(shape(gradient))*kind(gradient)
         deallocate(gradient, stat=i_stat)
         call memocc(i_stat, i_all, 'gradient', subname)

         i_all=-product(shape(gradientOld))*kind(gradientOld)
         deallocate(gradientOld, stat=i_stat)
         call memocc(i_stat, i_all, 'gradientOld', subname)

         i_all=-product(shape(overlapPsiGuessP))*kind(overlapPsiGuessP)
         deallocate(overlapPsiGuessP, stat=i_stat)
         call memocc(i_stat, i_all, 'overlapPsiGuessP', subname)

         i_all=-product(shape(overlapPad))*kind(overlapPad)
         deallocate(overlapPad, stat=i_stat)
         call memocc(i_stat, i_all, 'overlapPad', subname)

         i_all=-product(shape(HamPad))*kind(HamPad)
         deallocate(HamPad, stat=i_stat)
         call memocc(i_stat, i_all, 'HamPad', subname)

         i_all=-product(shape(psiGuessP))*kind(psiGuessP)
         deallocate(psiGuessP, stat=i_stat)
         call memocc(i_stat, i_all, 'psiGuessP', subname)

         i_all=-product(shape(norbtotpArr))*kind(norbtotpArr)
         deallocate(norbtotpArr, stat=i_stat)
         call memocc(i_stat, i_all, 'norbtotpArr', subname)

         i_all=-product(shape(alphaArr))*kind(alphaArr)
         deallocate(alphaArr, stat=i_stat)
         call memocc(i_stat, i_all, 'alphaArr', subname)

         i_all=-product(shape(rayleigh))*kind(rayleigh)
         deallocate(rayleigh, stat=i_stat)
         call memocc(i_stat, i_all, 'rayleigh', subname)

         i_all=-product(shape(psiGuess))*kind(psiGuess)
         deallocate(psiGuess, stat=i_stat)
         call memocc(i_stat, i_all, 'psiGuess', subname)
      end if

      allocate(psiGuessP(norbtotPad*nspinor,max(norbpArr(iproc),1),nspin), stat=i_stat)
      call memocc(i_stat, psiGuessP, 'psiGuessP', subname)

      allocate(overlapPad(norbtotPad*nspinor,norbtotPad,nspin,kp), stat=i_stat)
      call memocc(i_stat, overlapPad, 'overlapPad', subname)

      allocate(HamPad(norbtotPad*nspinor,norbtotPad,nspin,kp), stat=i_stat)
      call memocc(i_stat, HamPad, 'HamPad', subname)

      allocate(overlapPsiGuessP(norbtotPad*nspinor,max(norbpArr(iproc),1),nspin), stat=i_stat)
      call memocc(i_stat, overlapPsiGuessP, 'overlapPsiGuessP', subname)

      allocate(gradient(norbtotPad*nspinor,max(norbpArr(iproc),1)), stat=i_stat)
      call memocc(i_stat, gradient, 'gradient', subname)

      allocate(gradientOld(norbtotPad*nspinor,max(norbpArr(iproc),1)), stat=i_stat)
      call memocc(i_stat, gradientOld, 'gradientOldPaddded', subname)

      allocate(psiGuess(norbtot*nspinor,(max(orbs%norbu,orbs%norbd)-norbsc)*orbs%nkpts,nspin), stat=i_stat)
      call memocc(i_stat, psiGuess, 'psiGuess', subname)

      allocate(alphaArr(max(norbpArr(iproc),1)), stat=i_stat)
      call memocc(i_stat, alphaArr, 'alphaArr', subname)

      allocate(rayleigh(max(norbpArr(iproc),1)), stat=i_stat)
      call memocc(i_stat, rayleigh, 'rayleigh', subname)

      if(.not.simul) then
         allocate(norbtotpArr(0:nprocSub-1), stat=i_stat)
         call memocc(i_stat, norbtotpArr, 'norbtotpArr', subname)
      else
         if(0<=iproc .and. iproc<nprocSubu) then
            allocate(norbtotpArr(0:nprocSubu-1), stat=i_stat)
            call memocc(i_stat, norbtotpArr, 'norbtotpArr', subname)
         else 
            allocate(norbtotpArr(nprocSubu:nprocSubu+nprocSubd-1), stat=i_stat)
            call memocc(i_stat, norbtotpArr, 'norbtotpArr', subname)
         end if
      end if

      if(ispin==1 .or. simul) then
         allocate(sortArr(1:max(orbs%norbu,orbs%norbd)-norbsc,orbs%nkpts,nspin), stat=i_stat)
         call memocc(i_stat, sortArr, 'sortArr', subname)
         sortArr=0.d0
      end if


      ! The input guess is possibly generated only with a part of all MPI processes.
      ! The following if construct is therefore only executed by the active processes. The other
      ! 'wait' at the end of this if construct at a MPI barrier.
      processIf: if(iproc<nprocSub .or. (simul .and. iproc<nprocSubu+nprocSubd)) then
         if(simul) then
            ! For the case where simul is true, there are two communicators newCommu and newCommd.
            ! For simplicity only one name (newComm) will be used in the following.
            if(0<=iproc .and. iproc<nprocSubu) then
               newComm=newCommu
            else
               newComm=newCommd
            end if
         end if


         ! NorbtotpArr gives the number of non-zero atomic orbitals that each process has when the wavefunction 
         ! is transposed. The total number (including zeros) would be given by norbtotPad/nprocSub.
         if(.not.simul .or. (0<=iproc .and. iproc<nprocSubu)) then
            ! This is the up case if simul or, if simul is false, the up case in the first iteration of
            ! spinLoop and the down case in the second iteration.
            tt=norbtot/dble(nprocSub)
            ii=floor(tt)
            ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
            norbtotpArr=ii
            kk=norbtot-nprocSub*ii
            norbtotpArr(0:kk-1)=ii+1

            ! Check wheter this distribution is correct
            ii=0
            do iorb=0,nprocSub-1
               ii=ii+norbtotpArr(iorb)
            end do
            if(ii/=norbtot) then
               if(iproc==0) write(*,'(3x,a)') 'ERROR: wrong partition of norbtot'   
               stop
            end if
         else
            ! This is the down case if simul is true; if simul is false, this branch is not reached.
            tt=norbtot/dble(nprocSub)
            ii=floor(tt)
            ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
            norbtotpArr=ii
            kk=norbtot-nprocSub*ii
            norbtotpArr(nprocSubu:nprocSubu+kk-1)=ii+1

            ! Check wheter this distribution is correct
            ii=0
            do iorb=0,nprocSub-1
               ii=ii+norbtotpArr(nprocSubu+iorb)
            end do
            if(ii/=norbtot) then
               if(iproc==0) write(*,'(3x,a)') 'ERROR: wrong partition of norbtot'   
               stop
            end if
         end if


         ! Random initialization of psiGuessP. This array contains the orbitals treated by the current process.
         ! Initialize the random number generator. Make sure that it initialized differently for each iproc and ispin.
         call initRandomSeed(iproc, ispin)
         do iorb=1,norbpArr(iproc)
            do j=1,norbtot*nspinor
               call random_number(ttreal)
               psiGuessP(j,iorb,ispin)=dble(ttreal)
            end do
            ! Pad with zeros.
            do j=norbtot*nspinor+1,norbtotPad*nspinor
               psiGuessP(j,iorb,ispin)=0.d0
            end do
         end do


         ! Pad the Hamiltonian and the overlap matrix with zeroes, for all k-points.
         ! Skip the semi core orbitals.
         ii=0
         do i=1,natsc
            ii=ii+norbscArr(i,ispin)**2*nspinor
         end do
         ! First copy all columns and pad them with zeros.
         do ikpt=1,kp
            kk=ikpt+minval(kArr)-1
            do iorb=1,norbtot
               call vcopy(norbtot*nspinor, hamovr(ii+(iorb-1)*norbtot*nspinor+1,ispin,1,kk), 1,&
                  &   HamPad(1,iorb,ispin,ikpt), 1)
               call vcopy(norbtot*nspinor, hamovr(ii+(iorb-1)*norbtot*nspinor+1,ispin,2,kk), 1,&
                  &   overlapPad(1,iorb,ispin,ikpt), 1)
               do j=norbtot*nspinor+1,norbtotPad*nspinor
                  HamPad(j,iorb,ispin,ikpt)=0.d0
                  overlapPad(j,iorb,ispin,ikpt)=0.d0
               end do
            end do
            ! Now add the columns containing only zeros.
            do iorb=norbtot+1,norbtotPad
               do j=1,norbtotPad*nspinor
                  HamPad(j,iorb,ispin,ikpt)=0.d0
                  overlapPad(j,iorb,ispin,ikpt)=0.d0
               end do
            end do
         end do


         ! Calculate the matrix product overlapPad*psiGuessP=overlapPsiGuessP.
         do ikpt=1,kp
            if(nspinor==1) then
               call dsymm('l', 'u', norbtotPad, nkArr(ikpt), 1.d0, overlapPad(1,1,ispin,ikpt),&
                  &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                  &   norbtotPad, 0.d0, overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
            else
               call zhemm('l', 'l', norbtotPad, nkArr(ikpt), (1.d0,0.d0), overlapPad(1,1,ispin,ikpt),&
                  &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                  &   norbtotPad, (0.d0,0.d0), overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
            end if
         end do
         ! Orthonormalize the orbitals.
         blocksize=-1 ; blocksizeSmall=-1
         if(.not. simul .or. iproc<nprocSubu) then
            call orthonormalizePsi(iproc, nprocSub, norbtotPad, norb, norbpArr(iproc), &
               &   norbpArr(0), norbtotpArr(0), psiGuessP(1,1,ispin), &
               &   overlapPsiGuessP(1,1,ispin), newComm, orthpar, orbs, 0, nspinor, blocksize, blocksizeSmall)
         else
            call orthonormalizePsi(iproc, nprocSub, norbtotPad, norb, norbpArr(iproc), &
               &   norbpArr(nprocSubu), norbtotpArr(nprocSubu), &
               &   psiGuessP(1,1,ispin), overlapPsiGuessP(1,1,ispin), newComm, orthpar, &
               &   orbs, nprocSubu, nspinor, blocksize, blocksizeSmall)
         end if


         ! Improve the input guess for the orbitals. To do this, calculate for each orbital the gradient g which is given as
         !   |g>=H|psi>-e|psi>, where e=<psi|H|psi>/<psi|S|psi> (with psi being an orbital, H the Hamiltonian and S the overlap matrix).
         ! The eigenvectors are then improved by following this gradient using a steepest descent with variabel step size. After each step,
         ! the orbitals are orthonormalized.

         ! Initialize some variables
         itermax=500 
         success=.false.
         alphaArr=5.d-1  ! the initial step size
         tol=orthpar%iguessTol  ! the criterion for exiting the loop

         if(.not. simul) then
            if(iproc==0) write(*,'(5x,a)',advance='no') 'Iteratively determining eigenvectors... '
         else
            if(iproc==0) write(*,'(3x,a)',advance='no') 'Iteratively determining eigenvectors... '
         end if

         ! This is the main loop
         mainLoop: do iter=1,itermax
            do ikpt=1,kp
               ! Calculate the matrix products HamPad*psiGuessP and overlapPad*psiGuessP
               if(nspinor==1) then
                  call dsymm('l', 'u', norbtotPad, nkArr(ikpt), 1.d0, HamPad(1,1,ispin,ikpt),&
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                     &   norbtotPad, 0.d0, gradient(1,kstArr(ikpt)), norbtotPad)
                  call dsymm('l', 'u', norbtotPad, nkArr(ikpt), 1.d0, overlapPad(1,1,ispin,ikpt),&
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                     &   norbtotPad, 0.d0, overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
               else
                  call zhemm('l', 'l', norbtotPad, nkArr(ikpt), (1.d0,0.d0), HamPad(1,1,ispin,ikpt), &
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                     &   norbtotPad, (0.d0,0.d0), gradient(1,kstArr(ikpt)), norbtotPad)
                  call zhemm('l', 'l', norbtotPad, nkArr(ikpt), (1.d0,0.d0), overlapPad(1,1,ispin,ikpt),&
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                     &   norbtotPad, (0.d0,0.d0), overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
               end if
            end do
            ! Calculate the rayleigh quotients rayleigh(i)=<psi|H|psi>/<psi|S|psi>. Since the orbitals are normalized, <psi|S|psi> is equal to one.
            do iorb=1,norbpArr(iproc)
               if(nspinor==1) then
                  rayleigh(iorb)=ddot(norbtotPad, psiGuessP(1,iorb,ispin), 1, gradient(1,iorb), 1)
               else
                  zz=zdotc(norbtotPad, psiGuessP(1,iorb,ispin), 1, gradient(1,iorb), 1)
                  call vcopy(1, zz, 1, rayleigh(iorb), 1)
               end if

               if(iorb==1) then
               end if
            end do

            ! Calculate the gradient |g(i)>=H|psi(i)>-rayleigh(i)*S|psi(i)>
            gradientMax=0.d0
            do iorb=1,norbpArr(iproc)
               call daxpy(norbtotPad*nspinor, -rayleigh(iorb), overlapPsiGuessP(1,iorb,ispin), 1, gradient(1,iorb), 1)

               if(iter>2) then
                  ! Adapt the step size if iter>2.
                  ! First determine the angle between to consecutive gradients to adapt the step size.
                  if(nspinor==1) then
                     cosangle=ddot(norbtotPad, gradient(1,iorb), 1, gradientOld(1,iorb), 1)
                     cosangle=cosangle/(dnrm2(norbtotPad, gradient(1,iorb), 1)*dnrm2(norbtotPad, gradientOld(1,iorb), 1))
                  else
                     zz=zdotc(norbtotPad, gradient(1,iorb), 1, gradientOld(1,iorb), 1)
                     call vcopy(1, zz, 1, cosangle, 1)
                     zz=zdotc(norbtotPad, gradient(1,iorb), 1, gradient(1,iorb), 1)
                     call vcopy(1, zz, 1, tt, 1)
                     tt=sqrt(tt)
                     cosangle=cosangle/tt
                     zz=zdotc(norbtotPad, gradientOld(1,iorb), 1, gradientOld(1,iorb), 1)
                     call vcopy(1, zz, 1, tt, 1)
                     tt=sqrt(tt)
                     cosangle=cosangle/tt
                  end if
                  ! Adapt the gradient depending on the value of this angle.
                  if(cosangle>9.d-1) then
                     alphaArr(iorb)=alphaArr(iorb)*1.1d0
                  else
                     alphaArr(iorb)=alphaArr(iorb)*6.d-1
                  end if
               end if
               ! Copy the current gradient to gradientOld.
               call vcopy(norbtotPad*nspinor, gradient(1,iorb), 1, gradientOld(1,iorb), 1)
               ! Calculate the square of the norm of the gradient.
               if(nspinor==1) then
                  gradientNorm=ddot(norbtotPad,gradient(1,iorb),1,gradient(1,iorb),1)
               else
                  zz=zdotc(norbtotPad,gradient(1,iorb),1,gradient(1,iorb),1)
                  call vcopy(1, zz, 1, gradientNorm, 1)
               end if
               ! Determine the maximal gradient norm among all vectors treated by this processor.
               if(gradientNorm>gradientMax(2)) gradientMax(2)=gradientNorm
            end do
            ! Determine the maximal gradient norm among all orbitals.
            if(nprocSub>1) then
               call timing(iproc, 'Input_comput', 'OF')
               call timing(iproc, 'Input_commun', 'ON')
               call mpi_allreduce(gradientMax(2), gradientMax(1), 1, mpi_double_precision, mpi_max, newComm, ierr)
               call timing(iproc, 'Input_commun', 'OF')
               call timing(iproc, 'Input_comput', 'ON')
            else
               gradientMax(1)=gradientMax(2)
            end if
            ! Exit if the maximal gradient norm among all orbitals is smaller than tol.
            if(gradientMax(1)<tol) success=.true.

            if(success) then
               if(iproc==0) write(*,'(a,i0,a)') 'done in ',iter-1,' iterations.'
               exit mainLoop
            end if

            ! Improve the eigenvectors by following the gradient using steepest descent.
            do iorb=1,norbpArr(iproc)
               call daxpy(norbtotPad*nspinor, -alphaArr(iorb), gradient(1,iorb), 1, psiGuessP(1,iorb,ispin), 1)
            end do


            ! Calculate the matrix product overlapPad*psiGuessP=overlapPsiGuessP
            do ikpt=1,kp
               if(nspinor==1) then
                  call dsymm('l', 'u', norbtotPad, nkArr(ikpt), 1.d0, overlapPad(1,1,ispin,ikpt),&
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin),&
                     &   norbtotPad, 0.d0, overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
               else
                  call zhemm('l', 'l', norbtotPad, nkArr(ikpt), (1.d0,0.d0), overlapPad(1,1,ispin,ikpt),&
                     &   norbtotPad, psiGuessP(1,kstArr(ikpt),ispin), &
                     &   norbtotPad, (0.d0,0.d0), overlapPsiGuessP(1,kstArr(ikpt),ispin), norbtotPad)
               end if
            end do
            ! Orthonormalize the orbitals
            if(.not. simul .or. iproc<nprocSubu) then
               call orthonormalizePsi(iproc, nprocSub, norbtotPad, norb, norbpArr(iproc), norbpArr(0), &
                  &   norbtotpArr(0), psiGuessP(1,1,ispin), &
                  &   overlapPsiGuessP(1,1,ispin), newComm, orthpar, orbs, 0, nspinor, blocksize, blocksizeSmall)
            else
               call orthonormalizePsi(iproc, nprocSub, norbtotPad, norb, norbpArr(iproc), norbpArr(nprocSubu),&
                  &   norbtotpArr(nprocSubu), &
                  &   psiGuessP(1,1,ispin), overlapPsiGuessP(1,1,ispin), newComm,orthpar, orbs, nprocSubu, &
                  &   nspinor, blocksize, blocksizeSmall)
            end if

         end do mainLoop

         ! Write a warning in case no convergence was reached within the allowed number of iterations.
         if(.not. success) then
            if(iproc==0) call yaml_warning('No convergence after' // trim(yaml_toa(itermax)) // &
               & ' iterations. gradientMax=' // trim(yaml_toa(gradientMax(1),fmt='(es9.3)')))
            !if(iproc==0) write(*,'(a,i0,a,es9.3)') 'WARNING: no convergence after ',itermax,&
            !   &   ' iterations. gradientMax=',gradientMax(1)
         end if

      end if processIf

      ! Here the processes that are not involved in the input guess wait for the other processes.
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)


      ! Allocate the arrays needed for distributing the eigenvectors and eigenvalues to all processes.
      if(ispin==1 .or. simul) then
         allocate(sendcounts(0:nproc-1), stat=i_stat)
         call memocc(i_stat, sendcounts, 'sendcounts', subname)

         allocate(recvcounts(0:nproc-1), stat=i_stat)
         call memocc(i_stat, recvcounts, 'recvcounts', subname)

         allocate(sdispls(0:nproc-1), stat=i_stat)
         call memocc(i_stat, sdispls, 'sdispls', subname)

         allocate(rdispls(0:nproc-1), stat=i_stat)
         call memocc(i_stat, rdispls, 'rdispls', subname)
      end if


      ! Send all eigenvalues to all processes. These values were calculated in the last iteration of mainLoop.
      ! First define the arrays needed for mpi_allgatherv.
      if(nproc>1) then
         ii=0
         if(.not.simul) then
            do i=0,nproc-1
               recvcounts(i)=norbpArr(i)
               rdispls(i)=ii
               ii=ii+norbpArr(i)
            end do
         else
            do i=0,nproc-1
               recvcounts(i)=norbpArrSimul(i)
               rdispls(i)=ii
               ii=ii+norbpArrSimul(i)
            end do
         end if
         if(.not.simul) then
            !ist=(ispin-1)*norb+1
            if(ispin==1) ist=1
            if(ispin==2) ist=norbArr(1)+1
            call timing(iproc, 'Input_comput', 'OF')
            call timing(iproc, 'Input_commun', 'ON')
            call mpi_allgatherv(rayleigh(1), norbpArr(iproc), mpi_double_precision, orbs%eval(ist), &
               &   recvcounts, rdispls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
            call timing(iproc, 'Input_commun', 'OF')
            call timing(iproc, 'Input_comput', 'ON')
            ii=ii+norb/orbs%nkpts
         else
            call timing(iproc, 'Input_comput', 'OF')
            call timing(iproc, 'Input_commun', 'ON')
            call mpi_allgatherv(rayleigh(1), norbpArrSimul(iproc), mpi_double_precision, orbs%eval(1), &
               &   recvcounts, rdispls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
            call timing(iproc, 'Input_commun', 'OF')
            call timing(iproc, 'Input_comput', 'ON')
            ii=ii+norb/orbs%nkpts
         end if
      else
         if(simul) stop 'should not happen...'
         ii=0
         do ikpt=1,orbs%nkpts
            call vcopy(norb, rayleigh(1), 1, orbs%eval(1+ii), 1)
            ii=ii+1
         end do
      end if


      ! Send all eigenvectors to all processes. Since these eigenvectors contain some padded zeros, we can first
      ! cut off these zeros.
      allocate(psiGuessPTrunc(norbtot*nspinor,max(norbpArr(iproc),1),nspin), stat=i_stat)
      call memocc(i_stat, psiGuessPTrunc, 'psiGuessPTrunc', subname)
      do iorb=1,norbpArr(iproc)
         call vcopy(norbtot*nspinor, psiGuessP(1,iorb,ispin), 1, psiGuessPTrunc(1,iorb,ispin), 1)
      end do
      ! Define the values necessary for mpi_allgatherv
      if(nproc>1) then
         ii=0
         if(.not.simul) then
            do i=0,nproc-1
               recvcounts(i)=norbtot*nspinor*norbpArr(i)
               rdispls(i)=ii
               ii=ii+recvcounts(i)
            end do
         else
            kk=0
            do i=0,nproc-1
               recvcounts(i)=norbtot*nspinor*norbpArrSimul(i)
               rdispls(i)=ii
               kk=kk+norbpArrSimul(i)
               ii=ii+norbtot*nspinor*norbpArrSimul(i)
               ! Skip 'emtpy' vectors in the case that there are more down orbitals than up orbitals.
               ! This makes sure that the down orbitals are copied to the right place.
               if(orbs%norbd>orbs%norbu .and. mod(kk,orbs%norbu-norbsc)==0) then
                  !if(iproc==0) write(*,'(a,i0)') 'skip at ',kk
                  ii=ii+norbtot*nspinor*(orbs%norbd-orbs%norbu)
               end if
            end do
         end if
         if(.not.simul) then
            call timing(iproc, 'Input_comput', 'OF')
            call timing(iproc, 'Input_commun', 'ON')
            call mpi_allgatherv(psiGuessPTrunc(1,1,ispin), norbtot*nspinor*norbpArr(iproc), mpi_double_precision, &
               &   psiGuess(1,1,ispin), recvcounts, rdispls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
            call timing(iproc, 'Input_commun', 'OF')
            call timing(iproc, 'Input_comput', 'ON')
         else
            call timing(iproc, 'Input_comput', 'OF')
            call timing(iproc, 'Input_commun', 'ON')
            call mpi_allgatherv(psiGuessPTrunc(1,1,ispin), norbtot*nspinor*norbpArrSimul(iproc), mpi_double_precision, &
               &   psiGuess(1,1,1), recvcounts, rdispls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
            call timing(iproc, 'Input_commun', 'OF')
            call timing(iproc, 'Input_comput', 'ON')
         end if
      else
         if(simul) stop 'should not happen...'
         call vcopy(norbtot*nspinor*norb, psiGuessP(1,1,ispin), 1, psiGuess(1,1,ispin), 1)
      end if
      i_all=-product(shape(psiGuessPTrunc))*kind(psiGuessPTrunc)
      deallocate(psiGuessPTrunc, stat=i_stat)
      call memocc(i_stat, i_all, 'psiGuessPTrunc', subname)

      ! Transform the eigenvectors to the wavelet basis.
      ! These are the starting indices of the vectors: istpsi is the starting vector for psi
      ! istpsit that for psiGuessWavelet. For the case where simul is true, we use istpsiS and
      ! istpsitS, respectively. If there are any semicore orbitals, we skip them using ishift
      ! (for psiGuessWavelet) and ishift2 (for psi).
      istpsi=1
      istpsit=1
      istpsitS=1
      istpsiS=1
      ishift=0
      ishift2=0

      !if(.not. simul) then
      !   if(iproc==0) write(*,'(5x,a)',advance='no') 'Transforming to wavelet basis... '
      !else
      !   if(iproc==0) write(*,'(3x,a)',advance='no') 'Transforming to wavelet basis... '
      !end if

      ! First make a loop over the k points handled by this process.
      do ikptp=1,orbs%nkptsp
         ! ikpt is the index of the k point.
         ikpt=orbs%iskpts+ikptp
         ! nvctrp is the length of the vectors.
         nvctrp=comms%nvctr_par(iproc,ikpt)

         ! If simul is false and ispin==2, we are treating the down orbitals. Therefore we have to skip
         ! the up orbitals. We have to do this skipping for each k-point since in memory the spin up and 
         ! down orbitals are consecutive for each k-point.
         ! (i.e. schematically: in memory: (orbUp(k1),orbDown(k1)),(orbUp(k2),orbDown(k2)),...)
         if(ispin==2) istpsi=istpsi+nvctrp*norbtot*orbs%nspinor
         if(ispin==2) istpsit=istpsit+nvctrp*norbArr(1)/orbs%nkpts*orbs%nspinor
         if (nvctrp == 0) cycle

         if(.not.simul) then
            ! Now transform to wavelets.
            ! Skip the semicore orbitals.
            ishift=ishift+norbsc*nvctrp*orbs%nspinor
            ishift2=ishift2+sum(norbscArr(1:natsc,1))*nvctrp*orbs%nspinor
            if(ispin==2) ishift=ishift+norbsc*nvctrp*orbs%nspinor
            if(ispin==2) ishift2=ishift2+sum(norbscArr(1:natsc,2))*nvctrp*orbs%nspinor
            if(nspinor==1) then
               ii=(ikpt-1)*norb/orbs%nkpts+1
               call dgemm('n', 'n', nvctrp, norb, norbtot, 1.d0, psi(ishift2+istpsi), nvctrp, psiGuess(1,ii,ispin), &
                  &   norbtot, 0.d0, psiGuessWavelet(ishift+istpsit), nvctrp)
            else
               ii=(ikpt-1)*norb/orbs%nkpts+1
               call zgemm('n', 'n', nvctrp, norb/orbs%nkpts, norbtot, (1.d0,0.d0), psi(ishift2+istpsi),&
                  &   nvctrp, psiGuess(1,ii,ispin), &
                  &   norbtot, (0.d0,0.d0), psiGuessWavelet(ishift+istpsit), nvctrp)
            end if
            if(ispin==1 .and. nspin==2) ishift=ishift+norbsc*nvctrp*orbs%nspinor
            if(ispin==1 .and. nspin==2) ishift2=ishift2+sum(norbscArr(1:natsc,1))*nvctrp*orbs%nspinor
         else
            ! If simul is true, we can treat spin up and down using a simple loop. Since spin up and down are consecutive in memory,
            ! we can just advance the index without skipping.
            do ispin=1,nspin
               if(ispin==1) then
                  norb=(orbs%norbu-norbsc)*orbs%nkpts
               else
                  norb=(orbs%norbd-norbsc)*orbs%nkpts
               end if
               ishift=ishift+norbsc*nvctrp*orbs%nspinor
               ishift2=ishift2+sum(norbscArr(1:natsc,ispin))*nvctrp*orbs%nspinor
               if(nspinor==1) then
                  ii=(ikpt-1)*norb/orbs%nkpts+1
                  call dgemm('n', 'n', nvctrp, norb, norbtot, 1.d0, psi(ishift2+istpsiS), nvctrp, &
                     &   psiGuess(1,ii,ispin), norbtot, 0.d0, psiGuessWavelet(ishift+istpsitS), nvctrp)
               else
                  ii=(ikpt-1)*(max(orbs%norbu,orbs%norbd)-norbsc)+1
                  call zgemm('n', 'n', nvctrp, norb/orbs%nkpts, norbtot, (1.d0,0.d0), psi(ishift2+istpsiS), nvctrp, &
                     &   psiGuess(1,ii,ispin), norbtot, (0.d0,0.d0), psiGuessWavelet(ishift+istpsitS), nvctrp)
               end if
               istpsitS=istpsitS+norb/orbs%nkpts*nvctrp*nspinor
               istpsiS=istpsiS+norbtot*nvctrp*nspinor
            end do
         end if

         ! This is the second part of the skipping procedure for the case simul=.false.
         ! If ispin==1, we move the indices to the next up orbitals (i.e. skipping the down orbitals);
         ! if ispin==2, we advance the indices to the beginning of the next up orbitals and skip them
         ! in the next iteration of ikptp (this can not be done here since nvctrp may be different for
         ! the next k-point).
         if(ispin==1) istpsi=istpsi+nvctrp*norbtot*orbs%nspinor*nspin
         if(ispin==2) istpsi=istpsi+nvctrp*norbtot*orbs%nspinor
         if(ispin==1) istpsit=istpsit+nvctrp*norbArr(1)/orbs%nkpts*orbs%nspinor
         if(nspin==2 .and. ispin==1) istpsit=istpsit+nvctrp*norbArr(2)/orbs%nkpts*orbs%nspinor
         if(ispin==2) istpsit=istpsit+nvctrp*norb/orbs%nkpts*orbs%nspinor
      end do

      if(iproc==0) call yaml_map('Transforming to wavelet basis',.true.)
      !if(iproc==0) write(*,'(a)') 'done.'

   end do spinLoop

   ! The eigenvalues are now stored in the following way:
   ! (e11)(e21)(e31)(e12)(e22)(e32), where (eij) means the eigenvalues of k-point i and spin j.
   ! For further processing they have to be rearranged:
   ! (e11)(e12)(e21)(e22)(e31)(e32). 
   ! Use alphaArr as temporary array
   i_all=-product(shape(alphaArr))*kind(alphaArr)
   deallocate(alphaArr, stat=i_stat)
   call memocc(i_stat, i_all, 'alphaArr', subname)

   allocate(alphaArr((orbs%norb-norbsc)*orbs%nkpts*nspin), stat=i_stat)
   call memocc(i_stat, alphaArr, 'alphaArr', subname)

   call vcopy((orbs%norb-norbsc)*orbs%nkpts*nspin, orbs%eval(1), 1, alphaArr(1), 1)
   ist=1
   do ikpt=1,orbs%nkpts
      kk=(ikpt-1)*(orbs%norbu-norbsc)+1
      do ispin=1,nspin
         if(ispin==1) ii=orbs%norbu-norbsc
         if(ispin==2) ii=orbs%norbd-norbsc
         call vcopy(ii, alphaArr(kk), 1, orbs%eval(ist), 1)
         kk=kk+ii*orbs%nkpts
         ist=ist+ii
      end do
   end do


   ! Now treat the semicore orbitals, if there are any.
   semicoreIf: if(natsc>0) then
      !if(iproc==0) write(*,'(3x,a)',advance='no') 'Generating input guess for semicore orbitals...'

      if(nspinor == 1) then
         ncplx=1
         ncomp=1
         elseif(nspinor == 2) then
         ncplx=2
         ncomp=1
      else if (nspinor == 4) then
         ncplx=2
         ncomp=2
      end if


      ! Determine the number of semicore orbitals.
      ii=sum(norbscArr(1:natsc,1))
      if(nspin==2) jj=sum(norbscArr(1:natsc,2))

      ! Allocate the arrays which will contain the eigenvalues.
      if(nspin==1) allocate(evale(2*ii*nspin*orbs%nkpts), stat=i_stat)
      if(nspin==2) allocate(evale(2*max(ii,jj)*nspin*orbs%nkpts), stat=i_stat)
      call memocc(i_stat, evale, 'evale', subname)



      ist=1
      ! Diagonalize the semicore Hamiltonian for each k-point handled by this process.
      kptLoop: do ikptp=1,orbs%nkptsp
         ! ikpt is the number of the k-point.
         ikpt=orbs%iskpts+ikptp
         imatrst=1
         imatrst2=ndim_hamovr+1
         ! Make a loop over the number of semicore atoms.
         natscLoop: do iat=1,natsc

            ! norbi is the number of semicore orbitals for this atom.
            norbi=norbscArr(iat,1)
            if(nspinor==1) then

               ! Get the optimal work array size
               allocate(work(1), stat=i_stat)
               call memocc(i_stat, work, 'work', subname)
               call dsygv(1, 'v', 'u', norbi, hamovr(imatrst,1,1,ikpt), norbi, hamovr(imatrst,1,2,ikpt), &
                  &   norbi, evale(ist), work(1), -1, info)
               lwork = int(work(1))
               i_all=-product(shape(work))*kind(work)
               deallocate(work, stat=i_stat)
               call memocc(i_stat, i_all, 'work', subname)
               allocate(work(lwork), stat=i_stat)
               call memocc(i_stat, work, 'work', subname)

               call dsygv(1, 'v', 'u', norbi, hamovr(imatrst,1,1,ikpt), norbi, hamovr(imatrst,1,2,ikpt), &
                  &   norbi, evale(ist), work(1), lwork, info)
               if(info/=0) write(*,'(a,i0)') 'ERROR in dsygv, info=',info
               ist=ist+norbi
               if(nspin==1) then
                  i_all=-product(shape(work))*kind(work)
                  deallocate(work, stat=i_stat)
                  call memocc(i_stat, i_all, 'work', subname)
               end if

               ! Diagonalize the Hamiltonian for the down orbitals if we have spin polarization.
               if (nspin==2) then
                  ! norbj is the number of down orbitals.
                  norbj=norbscArr(iat,2)

                  call dsygv(1, 'v', 'u', norbj, hamovr(imatrst,2,1,ikpt), &
                     &   norbj, hamovr(imatrst,2,2,ikpt), norbj, evale(ist), work(1), lwork, info)


                  if(info/=0) write(*,'(a,i0)') 'ERROR in dsygv, info=',info
                  ist=ist+norbj
                  i_all=-product(shape(work))*kind(work)
                  deallocate(work, stat=i_stat)
                  call memocc(i_stat, i_all, 'work', subname)
               end if

            else

               ! Get the optimal work array size
               allocate(work(1), stat=i_stat)
               call memocc(i_stat, work, 'work', subname)
               allocate(rwork(1), stat=i_stat)
               call memocc(i_stat, rwork, 'rwork', subname)
               call zhegv(1, 'v', 'u', norbi, hamovr(imatrst,1,1,ikpt), norbi, hamovr(imatrst,1,2,ikpt), &
                  &   norbi, evale(ist), work(1), -1, work(1), info)
               lwork = int(work(1))
               i_all=-product(shape(work))*kind(work)
               deallocate(work, stat=i_stat)
               call memocc(i_stat, i_all, 'work', subname)
               i_all=-product(shape(rwork))*kind(rwork)
               deallocate(rwork, stat=i_stat)
               call memocc(i_stat, i_all, 'rwork', subname)
               allocate(work(lwork), stat=i_stat)
               call memocc(i_stat, work, 'work', subname)
               allocate(rwork(lwork), stat=i_stat)
               call memocc(i_stat, rwork, 'rwork', subname)

               call zhegv(1, 'v', 'u', norbi, hamovr(imatrst,1,1,ikpt), norbi, hamovr(imatrst,1,2,ikpt), &
                  &   norbi, evale(ist), work(1), lwork, rwork(1), info)
               if(info/=0) write(*,'(a,i0)') 'ERROR in zhegv, info=',info
               ist=ist+norbi
               if(nspin==1) then
                  i_all=-product(shape(work))*kind(work)
                  deallocate(work, stat=i_stat)
                  call memocc(i_stat, i_all, 'work', subname)
                  i_all=-product(shape(rwork))*kind(rwork)
                  deallocate(rwork, stat=i_stat)
                  call memocc(i_stat, i_all, 'rwork', subname)
               end if

               ! Diagonalize the Hamiltonian for the down orbitals if we have spin polarization.
               if (nspin==2) then
                  norbj=norbscArr(iat,2)
                  call zhegv(1, 'v', 'u', norbj, hamovr(imatrst,2,1,ikpt), &
                     &   norbj, hamovr(imatrst,2,2,ikpt), norbj, evale(ist), &
                     &   work(ist), lwork, rwork(1), info)
                  if(info/=0) write(*,'(a,i0)') 'ERROR in zhegv, info=',info
                  ist=ist+norbj
                  i_all=-product(shape(work))*kind(work)
                  deallocate(work, stat=i_stat)
                  call memocc(i_stat, i_all, 'work', subname)
                  i_all=-product(shape(rwork))*kind(rwork)
                  deallocate(rwork, stat=i_stat)
                  call memocc(i_stat, i_all, 'rwork', subname)
               end if

            end if

            ! Make sure that the eigenvectors of all MPI processes are the same.
            norbij=norbi
            ihs=imatrst
            do jspin=1,nspin
               do jjorb=1,norbij
                  ! If it is negative change the sign to all the values.
                  if (hamovr(ihs+(jjorb-1)*norbij*ncplx,jspin,1,ikpt) < 0.0_wp) then
                     do jiorb=1,norbij*ncplx
                        hamovr(ihs-1+jiorb+(jjorb-1)*norbij*ncplx,jspin,1,ikpt)=&
                           &   -hamovr(ihs-1+jiorb+(jjorb-1)*norbij*ncplx,jspin,1,ikpt)
                     end do
                  end if
               end do
               if (nspin==2) then
                  norbij=norbj
               end if
            end do
            imatrst=imatrst+ncplx*norbi**2
            if(nspin>1) imatrst2=imatrst2+ncplx*norbj**2

         end do natscLoop
      end do kptLoop


      ! Now transform the eigenvectors to the wavelet basis.
      ! ist is the starting index of the atomic orbitals (in wavelet basis), ist2 the starting index of
      ! the eigenvectors in wavelet basis, imatrst is the starting index in atomic orbitals basis.
      ist=1
      ist2=1
      ! Make a loop over all k-points handled by this process.
      do ikptp=1,orbs%nkptsp
         ! ikpt is the number of the k-point.
         ikpt=orbs%iskpts+ikptp
         ! nvctrp is the number of components for this k-point which are handled by this process.
         nvctrp=comms%nvctr_par(iproc,ikpt)
         ! Make a loop over the spins.
         do jspin=1,nspin
            imatrst=1
            ! Make a loop over the number of semicore atoms.
            do iat=1,natsc
               ! norbi is ths number of semicore orbitals for this atom.
               norbi=norbscArr(iat,jspin)
               if(nspinor==1) then
                  call gemm('n', 'n', nvctrp, norbi, norbi, 1.0_wp, psi(ist), max(1,nvctrp), &
                     &   hamovr(imatrst,jspin,1,ikpt), norbi, 0.0_wp, psiGuessWavelet(ist2), max(1,nvctrp))
               else
                  call c_gemm('n', 'n', ncomp*nvctrp, norbi, norbi, (1.0_wp,0.0_wp), &
                     &   psi(ist), max(1,ncomp*nvctrp), hamovr(imatrst,jspin,1,ikpt), norbi, &
                     &   (0.0_wp,0.0_wp), psiGuessWavelet(ist2), max(1,ncomp*nvctrp))
               end if
               ist=ist+norbi*nvctrp*nspinor
               ist2=ist2+norbi*nvctrp*nspinor
               imatrst=imatrst+ncplx*norbi**2
            end do
            ! Skip the non-semicore orbitals.
            if(jspin==1) ist2=ist2+nvctrp*(orbs%norbu-norbsc)*orbs%nspinor
            if(jspin==2) ist2=ist2+nvctrp*(orbs%norbd-norbsc)*orbs%nspinor
            ist=ist+nvctrp*norbscArr(natsc+1,jspin)*orbs%nspinor
         end do
      end do


      ! Now send all eigenvalues of all k-points to all processes.
      ! First find out which process handles which k-points.
      allocate(kpArr(orbs%nkpts,0:nproc-1,2), stat=i_stat)
      call memocc(i_stat, kpArr, 'kpArr', subname)
      allocate(sceval(norbsc*nspin*orbs%nkpts), stat=i_stat)
      call memocc(i_stat, sceval, 'sceval', subname)

      ! If process i handles k-point j, set kpArr(j,i) to 1. Then make a mpi_allreduce
      ! to collect this informations from all processes.
      kpArr=0
      do ikptp=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikptp
         kpArr(ikpt,iproc,2)=1
      end do
      call timing(iproc, 'Input_comput', 'OF')
      call timing(iproc, 'Input_commun', 'ON')
      call mpi_allreduce(kpArr(1,0,2), kparr(1,0,1), orbs%nkpts*nproc, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call timing(iproc, 'Input_commun', 'OF')
      call timing(iproc, 'Input_comput', 'ON')

      ! ist2 is the starting index of the collected orbitals in sceval.
      ist2=1
      ! Make a loop over all k-points.
      do ikpt=1,orbs%nkpts
         ! Find the process with lowest process ID handling this k-point.
         do i=0,nproc-1
            if(kpArr(ikpt,i,1)==1) then
               ! Process jproc is the process we want.
               jproc=i
               exit
            end if
         end do
         ! Send the eigenvalues from process jproc to all processes.
         ist=1
         if(iproc==jproc) then
            ! Determine the starting index of the required k-point among the orbitals
            ! handled by process jproc.
            do ikptp2=1,orbs%nkptsp
               ikpt2=orbs%iskpts+ikptp2
               if(ikpt2==ikpt) exit
               ist=ist+norbsc*nspin
            end do
            ! Copy the orbitals of this k-point to sceval.
            call vcopy(norbsc*nspin, evale(ist), 1, sceval(ist2), 1)
         end if
         ! Send the orbitals to all processes.
         call timing(iproc, 'Input_comput', 'OF')
         call timing(iproc, 'Input_commun', 'ON')
         call mpi_bcast(sceval(ist2), norbsc*nspin, mpi_double_precision, jproc, bigdft_mpi%mpi_comm, ierr)
         call timing(iproc, 'Input_commun', 'OF')
         call timing(iproc, 'Input_comput', 'ON')
         ist2=ist2+norbsc*nspin
      end do

      if(iproc==0) call yaml_map('Generating input guess for semicore orbitals',.true.)
      !if(iproc==0) write(*,'(a)') ' done.'

   end if semicoreIf



   ! Sort the (non-semicore) eigenvalues to print them. This requires to copy them first to sortArr since a given 
   ! eigenvalue has to remain associated with its eigenvector.
   ! This part is not implemented very efficiently, but it should not matter.
   if(.not.simul) then
      sortArr=0.d0
      ist=1
      do ikpt=1,orbs%nkpts
         do ispin=1,nspin
            if(ispin==1) ii=orbs%norbu-norbsc
            if(ispin==2) ii=orbs%norbd-norbsc
            call vcopy(ii, orbs%eval(ist), 1, sortArr(1,ikpt,ispin), 1)
            ist=ist+ii
         end do
      end do
      do ikpt=1,orbs%nkpts
         do ispin=1,nspin
            if(ispin==1) ii=orbs%norbu-norbsc
            if(ispin==2) ii=orbs%norbd-norbsc
            do iorb=1,ii
               do jorb=iorb,ii
                  if(sortArr(jorb,ikpt,ispin)<sortArr(iorb,ikpt,ispin)) then
                     tt=sortArr(iorb,ikpt,ispin)
                     sortArr(iorb,ikpt,ispin)=sortArr(jorb,ikpt,ispin)
                     sortArr(jorb,ikpt,ispin)=tt
                  end if
               end do
            end do
         end do
      end do
   else
      ist=1
      sortArr=0.d0
      do ikpt=1,orbs%nkpts
         do ispin=1,nspin
            if(ispin==1) then
               norb=orbs%norbu-norbsc
            else
               norb=orbs%norbd-norbsc
            end if
            call vcopy(norb, orbs%eval(ist), 1, sortArr(1,ikpt,ispin), 1)
            ist=ist+norb
         end do
      end do
      ist=1
      do ispin=1,nspin
         if(ispin==1) then
            norb=(orbs%norbu-norbsc)*orbs%nkpts
         else
            norb=(orbs%norbd-norbsc)*orbs%nkpts
         end if
         ist=ist+norb
         do ikpt=1,orbs%nkpts
            do iorb=1,norb/orbs%nkpts 
               do jorb=iorb,norb/orbs%nkpts
                  if(sortArr(jorb,ikpt,ispin)<sortArr(iorb,ikpt,ispin)) then
                     tt=sortArr(iorb,ikpt,ispin)
                     sortArr(iorb,ikpt,ispin)=sortArr(jorb,ikpt,ispin)
                     sortArr(jorb,ikpt,ispin)=tt
                  end if
               end do
            end do
         end do
      end do
   end if



   ! Now print out the eigenvalues. The semicore eigenvalues are stored in sceval, the
   ! non-semicore eigenvalues in sortArr.
   if(iproc==0) write(*,'(1x,a)') 'Sorted list of eigenvalues:'
   if(iproc==0) then
      if(nspin==1) then
         do ikpt=1,orbs%nkpts
            ii=0
            if(orbs%nkpts>1) write(*,'(/,3x,a,i0,a)') '------- k-point ',ikpt,' -------'
            do iat=1,natsc
               do j=1,norbscArr(iat,1)
                  ii=ii+1
                  if(iat/=natsc .or. j/=norbscArr(iat,1)) then
                     write(*,'(5x,a,i0,a,es14.5)') 'eval(',ii,')',sceval(ii)
                  else
                     write(*,'(5x,a,i0,a,es14.5,a)') 'eval(',ii,')',sceval(ii),' <--- last semicore orbital'
                  end if
               end do
            end do
            do iorb=1,norb/orbs%nkpts
               write(*,'(5x,a,i0,a,es14.5)') 'eval(',ii+iorb,')',sortArr(iorb,ikpt,1)
            end do
         end do
      else if(nspin==2) then
         do ikpt=1,orbs%nkpts
            ii=0
            if(orbs%nkpts>1) write(*,'(/,20x,a,i0,a)') '------- k-point ',ikpt,' -------'
            do iat=1,natsc
               do j=1,maxval(norbscArr(iat,:))
                  ii=ii+1
                  if(iat/=natsc .or. j/=norbscArr(iat,1)) then
                     write(*,'(5x,a,i0,a,es14.5,8x,a,i0,a,es14.5)') 'eval(',ii,',u)',evale(ii), &
                        &   'eval(',ii,',d)',sceval(norbscArr(iat,1)+ii)
                  else
                     write(*,'(5x,a,i0,a,es14.5,8x,a,i0,a,es14.5,a)') 'eval(',ii,',u)',evale(ii), &
                        &   'eval(',ii,',d)',sceval(norbscArr(iat,1)+ii),' <--- last semicore orbital'
                  end if
               end do
            end do
            do iorb=1,norb/orbs%nkpts
               write(*,'(5x,a,i0,a,es14.5,8x,a,i0,a,es14.5)') 'eval(',ii+iorb,',u)', &
                  &   sortArr(iorb,ikpt,1),'eval(',ii+iorb,',d)',sortArr(iorb,ikpt,2)
            end do
         end do
      end if
   end if


   ! Merge the semicore eigenvalues and non-semicore eigenvalues and store them in orbs%eval.
   ! This is obly necessary if there are some semicore orbitals.
   ! Use sortArr as temporary array.
   if(norbsc>0) then

      i_all=-product(shape(sortArr))*kind(sortArr)
      deallocate(sortArr, stat=i_stat)
      call memocc(i_stat, i_all, 'sortArr', subname)

      allocate(sortArr((orbs%norb-nspin*norbsc)*nspin*orbs%nkpts,1,1), stat=i_stat)
      call memocc(i_stat, sortArr, 'sortArr', subname)

      ! The starting indices:
      !   ist is the starting index of the merged eigenvalues
      !   ist2 is the starting index of the semicore orbitals
      !   ist3 is the starting index of the non-semicore orbitals.
      ist=1
      ist2=1
      ist3=1 
      ! First copy the non-semicore orbitals to the temporary array.
      call vcopy((orbs%norb-nspin*norbsc)*nspin, orbs%eval(1), 1, sortArr(1,1,1), 1)
      ! Make a loop over all k-points.
      do ikpt=1,orbs%nkpts
         ! Make a loop over spin up and down.
         do ispin=1,nspin
            if(ispin==1) then 
               norb=orbs%norbu-norbsc 
            else 
               norb=orbs%norbd-norbsc 
            end if
            ! Copy the semicore orbitals.
            call vcopy(norbsc, sceval(ist2), 1, orbs%eval(ist),1)
            ist=ist+norbsc
            ist2=ist2+norbsc
            ! Copy the non-semicore orbitals.
            call vcopy(norb, sortArr(ist3,1,1), 1, orbs%eval(ist),1)
            ist=ist+norb
            ist3=ist3+norb
         end do
      end do
   end if




   ! Deallocate all arrays.
   i_all=-product(shape(gradient))*kind(gradient)
   deallocate(gradient, stat=i_stat)
   call memocc(i_stat, i_all, 'gradient', subname)

   i_all=-product(shape(gradientOld))*kind(gradientOld)
   deallocate(gradientOld, stat=i_stat)
   call memocc(i_stat, i_all, 'gradientOld', subname)

   i_all=-product(shape(overlapPsiGuessP))*kind(overlapPsiGuessP)
   deallocate(overlapPsiGuessP, stat=i_stat)
   call memocc(i_stat, i_all, 'overlapPsiGuessP', subname)

   i_all=-product(shape(overlapPad))*kind(overlapPad)
   deallocate(overlapPad, stat=i_stat)
   call memocc(i_stat, i_all, 'overlapPad', subname)

   i_all=-product(shape(HamPad))*kind(HamPad)
   deallocate(HamPad, stat=i_stat)
   call memocc(i_stat, i_all, 'HamPad', subname)

   i_all=-product(shape(psiGuessP))*kind(psiGuessP)
   deallocate(psiGuessP, stat=i_stat)
   call memocc(i_stat, i_all, 'psiGuessP', subname)

   i_all=-product(shape(norbpArr))*kind(norbpArr)
   deallocate(norbpArr, stat=i_stat)
   call memocc(i_stat, i_all, 'norbpArr', subname)

   if(simul) then
      i_all=-product(shape(norbpArrSimul))*kind(norbpArrSimul)
      deallocate(norbpArrSimul, stat=i_stat)
      call memocc(i_stat, i_all, 'norbpArrSimul', subname)
   end if

   i_all=-product(shape(norbtotpArr))*kind(norbtotpArr)
   deallocate(norbtotpArr, stat=i_stat)
   call memocc(i_stat, i_all, 'norbtotpArr', subname)

   i_all=-product(shape(sendcounts))*kind(sendcounts)
   deallocate(sendcounts, stat=i_stat)
   call memocc(i_stat, i_all, 'sendcounts', subname)

   i_all=-product(shape(recvcounts))*kind(recvcounts)
   deallocate(recvcounts, stat=i_stat)
   call memocc(i_stat, i_all, 'recvcounts', subname)

   i_all=-product(shape(sdispls))*kind(sdispls)
   deallocate(sdispls, stat=i_stat)
   call memocc(i_stat, i_all, 'sdidpls', subname)

   i_all=-product(shape(rdispls))*kind(rdispls)
   deallocate(rdispls, stat=i_stat)
   call memocc(i_stat, i_all, 'rdidpls', subname)

   i_all=-product(shape(alphaArr))*kind(alphaArr)
   deallocate(alphaArr, stat=i_stat)
   call memocc(i_stat, i_all, 'alphaArr', subname)

   i_all=-product(shape(rayleigh))*kind(rayleigh)
   deallocate(rayleigh, stat=i_stat)
   call memocc(i_stat, i_all, 'rayleigh', subname)

   i_all=-product(shape(psiGuess))*kind(psiGuess)
   deallocate(psiGuess, stat=i_stat)
   call memocc(i_stat, i_all, 'psiGuess', subname)

   i_all=-product(shape(sortArr))*kind(sortArr)
   deallocate(sortArr, stat=i_stat)
   call memocc(i_stat, i_all, 'sortArr', subname)

   i_all=-product(shape(kstArr))*kind(kstArr)
   deallocate(kstArr, stat=i_stat)
   call memocc(i_stat, i_all, 'kstArr', subname)

   i_all=-product(shape(nkArr))*kind(nkArr)
   deallocate(nkArr, stat=i_stat)
   call memocc(i_stat, i_all, 'nkArr', subname)

   i_all=-product(shape(kArr))*kind(kArr)
   deallocate(kArr, stat=i_stat)
   call memocc(i_stat, i_all, 'kArr', subname)

   if(.not.simul) then
      i_all=-product(shape(newID))*kind(newID)
      deallocate(newID, stat=i_stat)
      call memocc(i_stat, i_all, 'newID', subname)
   end if

   if(simul) then
      i_all=-product(shape(newIDu))*kind(newIDu)
      deallocate(newIDu, stat=i_stat)
      call memocc(i_stat, i_all, 'newIDu', subname)

      i_all=-product(shape(newIDd))*kind(newIDd)
      deallocate(newIDd, stat=i_stat)
      call memocc(i_stat, i_all, 'newIDd', subname)
   end if

   i_all=-product(shape(norbArr))*kind(norbArr)
   deallocate(norbArr, stat=i_stat)
   call memocc(i_stat, i_all, 'norbArr', subname)

   if(natsc>0) then
      i_all=-product(shape(kpArr))*kind(kpArr)
      deallocate(kpArr, stat=i_stat)
      call memocc(i_stat, i_all, 'kpArr', subname)

      i_all=-product(shape(sceval))*kind(sceval)
      deallocate(sceval, stat=i_stat)
      call memocc(i_stat, i_all, 'sceval', subname)

      i_all=-product(shape(evale))*kind(evale)
      deallocate(evale, stat=i_stat)
      call memocc(i_stat, i_all, 'evale', subname)
   end if

   if(iproc==0) write(*,'(1x,a)') 'Input guess successfully generated.'

   ! Stop the timing for the input guess.
   call timing(iproc, 'Input_comput', 'OF')

END SUBROUTINE inputguessParallel


!> This subroutine orthonormalizes the orbitals psi in a parallel way. To do so, it first transposes the orbitals to all
!! processors using mpi_alltoallv. The orthonomalization is then done in this data layout using a combination of blockwise Gram-Schmidt
!! and Cholesky orthonomalization. At the end the vectors are again untransposed.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param iproc                  process ID
!!    @param nproc                  total number of processes
!!    @param norbtot                length of the vectors that have to be orthonormalized
!!    @param norb                   total number of vectors that have to be orthonomalized, shared over all processes
!!    @param norbp                  number of vectors held by current process
!!    @param norbpArr               array indicating the number of vectors held by each process (e.g. norbpArr(i) is the number of vectors
!!                                  held by process i)
!!    @param norbtotpArr            array indicating how many orbitals each process treats (without the padded zeros) if the wavefunctions
!!                                  are transposed
!!    @param overlapPsi             the vectors after the application of the overlap matrix
!!    @param newComm                the MPI communicator for the processes handling this orthonormalization (not necessarily all MPI processes)
!!    @param input                  contains some parameters
!!    @param simul                  if simul is true, the up and down orbitals are treated in parallel
!!    @param orbs                   type that contains many parameters concerning the orbitals
!!    @param nprocSt                starting process ID of the the processes in newComm
!!    @param nspinor                real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!    @param blocksize              the block size for the blockwise orthonormalization procesure
!!    @param blocksizeSmall         the block size for the orbitals that did not match into blocksize
!!  Input/Output arguments:
!!    @param psi                    on input: the vectors to be orthonormalized
!!                                  on output: the orthonomalized vectors
subroutine orthonormalizePsi(iproc, nproc, norbtot, norb, norbp, norbpArr,&
      &   norbtotpArr, psi, overlapPsi, newComm, orthpar, &
      &   orbs, nprocSt, nspinor, blocksize, blocksizeSmall)
   use module_base
   use module_types
   implicit none

   ! Calling arguments
   integer, intent(in):: iproc, nproc, norbtot, norb, norbp, newComm, nprocSt, nspinor
   integer, intent(in out):: blocksize
   integer, dimension(nprocSt:nprocSt+nproc-1),intent(in):: norbtotpArr, norbpArr
   real(kind=8), dimension(norbtot*norbp*nspinor),intent(in):: overlapPsi
   real(kind=8), dimension(norbtot*norbp*nspinor),intent(in out):: psi
   type(orthon_data), intent(in):: orthpar
   type(orbitals_data), intent(in) :: orbs

   ! Local variables
   integer:: i, j, iorb, iblock, jblock, ii, jj, ist, jst, iter, iter2, gcd,&
      &   blocksizeSmall, norbtotp, ierr, i_stat, i_all, getBlocksize
   real(kind=8),dimension(:),allocatable:: psiW, overlapPsiW, psiWTrans, overlapPsiWTrans
   integer,dimension(:),allocatable:: sendcounts, recvcounts, sdispls, rdispls
   character(len=*),parameter:: subname='orthonormalizePsi'

   !< This variable is the part of each orbital that will be distributed to each processor and will be used throughout the subroutine.
   norbtotp=norbtot/nproc

   ! Allocate all arrays
   allocate(psiW(norbtot*norbp*nspinor), stat=i_stat)
   call memocc(i_stat, psiW, 'psiW', subname)

   allocate(overlapPsiW(norbtot*norbp*nspinor), stat=i_stat)
   call memocc(i_stat, overlapPsiW, 'overlapPsiW', subname)

   allocate(psiWTrans(norbtotp*norb*nspinor), stat=i_stat)
   call memocc(i_stat, psiWTrans, 'psiWTrans', subname)

   allocate(overlapPsiWTrans(norbtotp*norb*nspinor), stat=i_stat)
   call memocc(i_stat, overlapPsiWTrans, 'overlapPsiWTrans', subname)

   allocate(sendcounts(nprocSt:nprocSt+nproc-1), stat=i_stat)
   call memocc(i_stat, sendcounts, 'sendcounts', subname)

   allocate(recvcounts(nprocSt:nprocSt+nproc-1), stat=i_stat)
   call memocc(i_stat, recvcounts, 'recvcounts', subname)

   allocate(sdispls(nprocSt:nprocSt+nproc-1), stat=i_stat)
   call memocc(i_stat, sdispls, 'sdispls', subname)

   allocate(rdispls(nprocSt:nprocSt+nproc-1), stat=i_stat)
   call memocc(i_stat, rdispls, 'rdispls', subname)


   if(nproc>1) then
      ! Rearrange psi and overlapPsi in memory such that they can be transposed using a single call to mpi_alltoallv.
      ! Here is a small example for illustration:
      ! Assume nproc=5, norbtot=10, norbtotpArr(i)=2 (i=0,1,2), norbtotpArr(i)=1, (i=3,4), norbpArr(i)=2 (i=0,1), norbpArr(i)=1 (i=2,3,4).
      ! Then the two vectors will be modified in the following way; for clarity we also show the same thing for the other processes.
      !       process 0                process 1            process 2        process 3        process 4
      !   1   9       1  13   |   17  25      17  29   |   33      33   |   41      41   |   49      49
      !   2  10       2  14   |   18  26      18  30   |   34      34   |   42      42   |   50      50
      !   3  11       9   7   |   19  27      25  23   |   35      35   |   43      43   |   51      51
      !   4  12      10   0   |   20  28      26   0   |   36      36   |   44      44   |   52      52
      !   5  13  =>   3   8   |   21  29  =>  19  24   |   37  =>  37   |   45  =>  45   |   53  =>  53
      !   6  14       4   0   |   22  30      20   0   |   38      38   |   46      46   |   54      54
      !   7  15      11  15   |   23  31      27  31   |   39      39   |   47      47   |   55      55
      !   8  16      12   0   |   24  32      28   0   |   40       0   |   48       0   |   56       0 
      !   0   0       5  16   |    0   0      21  32   |    0      40   |    0      48   |    0      56
      !   0   0       6   0   |    0   0      22   0   |    0       0   |    0       0   |    0       0 

      jj=1
      ii=0
      do i=nprocSt,nprocSt+nproc-1
         do iorb=0,norbp-1
            ! Copy the non-zero parts
            call vcopy(norbtotpArr(i)*nspinor, psi(ii+iorb*norbtot*nspinor+1), 1, psiW(jj), 1)
            call vcopy(norbtotpArr(i)*nspinor, overlapPsi(ii+iorb*norbtot*nspinor+1), 1, overlapPsiW(jj), 1)
            jj=jj+norbtotpArr(i)*nspinor
            do j=norbtotpArr(i)+1,norbtotp
               ! "Copy" the zeros. This happens only if norbtotpArr(i) < norbtotp
               psiW(jj)=0.d0
               if(nspinor>1) psiW(jj+1)=0.d0
               overlapPsiW(jj)=0.d0
               if(nspinor>1) overlapPsiW(jj+1)=0.d0
               jj=jj+nspinor
            end do
         end do
         ii=ii+norbtotpArr(i)*nspinor
      end do


      ! Define the values used for the call to mpi_alltoallv.
      ! Assume again the above example, then we would have:
      ! for iproc==0,1:   sendcounts(i)=4 (i=0,...,4)  ;  sdidpls(0)=0, sdipls(1)=4, sdipls(2)=8, sdipls(3)=12, sdipls(4)=16
      !                   recvcounts(i)=4 (i=0,1), recvcounts(i)=2 (i=2,3,4) ; rdidpls(0)=0, rdipls(1)=4, rdipls(2)=8, rdipls(3)=10, rdipls(4)=12
      ! for iproc==2,3,4: sendcounts(i)=2 (i=0,...,4)  ;  sdidpls(0)=0, sdipls(1)=2, sdipls(2)=4, sdipls(3)=6, sdipls(4)=8
      !                   recvcounts(i)=4 (i=0,1,2), recvcounts(i)=2 (i=3,4) ; rdidpls(0)=0, rdipls(1)=4, rdipls(2)=8, rdipls(3)=10, rdipls(4)=12
      ii=0 ; jj=0
      do i=nprocSt,nprocSt+nproc-1
         ! A given process sends the same amount if data to each process...
         sendcounts(i)=norbtotp*norbpArr(iproc)*nspinor
         sdispls(i)=ii
         ii=ii+sendcounts(i)
      end do
      do i=nprocSt,nprocSt+nproc-1
         ! ... but receives a different amount from each process.
         recvcounts(i)=norbtotp*norbpArr(i)*nspinor
         rdispls(i)=jj
         jj=jj+recvcounts(i)
      end do

      ! Transpose psiW (to psiWTrans) and overlapPsiW (to overlapPsiWTrans)
      ! Assuming the above example, this would result in the following vectors held by each process:
      !       process 0               process 1                process 2                process 3                process 4
      !  1  9 17 25 33 41 49  |   3 11 19 27 35 43 51  |   5 13 21 29 37 45 53  |   7 15 23 31 39 47 55  |   8 16 24 32 40 48 56
      !  2 10 18 26 34 42 50  |   4 12 20 28 36 44 52  |   6 14 22 30 38 46 54  |   0  0  0  0  0  0  0  |   0  0  0  0  0  0  0

      call timing(iproc, 'Input_comput', 'OF')
      call timing(iproc, 'Input_commun', 'ON')
      call mpi_alltoallv(psiW(1), sendcounts, sdispls, mpi_double_precision, psiWTrans(1), &
         &   recvcounts, rdispls, mpi_double_precision, newComm, ierr)
      call mpi_alltoallv(overlapPsiW(1), sendcounts, sdispls, mpi_double_precision, overlapPsiWTrans(1),&
         &   recvcounts, rdispls, mpi_double_precision, newComm, ierr)
      call timing(iproc, 'Input_commun', 'OF')
      call timing(iproc, 'Input_comput', 'ON')
   else
      call vcopy(norbtot*norbp*nspinor, psi(1), 1, psiWTrans(1), 1)
      call vcopy(norbtot*norbp*nspinor, overlapPsi(1), 1, overlapPsiWTrans(1), 1)
   end if


   ! Now orthonormalize the orbitals.
   ! There are two orthonormalization subroutines: gramschmidtOverlap orthogonalizes a given bunch of vectors to another bunch
   ! of already orthonormal vectors, and the subroutine choleskyOverlap orthonormalizes the given bunch.
   ! First determine how many bunches can be created for the given blocksize.
   if(blocksize==-1) blocksize=getBlocksize(orthpar, norb/orbs%nkpts)
   iter=floor(real(norb)/real(blocksize*orbs%nkpts))
   do iblock=1,iter
      ! ist is the starting vector of the current bunch. 
      ist=blocksize*(iblock-1)+1
      ! Now orthogonalize this bunch to all previous ones.
      do j=1,iblock-1
         ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
         jst=blocksize*(j-1)+1
         call gramschmidtOverlap(iproc, nproc, norbtotp, blocksize, psiWTrans(1),&
            &   overlapPsiwTrans(1), newComm, ist, jst, norb/orbs%nkpts, orbs%nkpts, nspinor)
      end do
      ! Orthonormalize the current bunch of vectors.
      call choleskyOverlap(iproc, nproc, norbtotp, blocksize, psiWTrans(1), overlapPsiWTrans(1), &
         &   newComm, ist, norb/orbs%nkpts, orbs%nkpts, nspinor)
   end do

   ! Orthonormalize the remaining vectors, if there are any.
   remainingIf: if(blocksize*iter/=norb/orbs%nkpts) then
      ! ist is the starting vector of the bunch that still havs to be orthonormalized.
      ist=blocksize*iter+1
      ! We have to find a new block size that matches both the remaining vectors and the already orthonomalized ones. This is done by determining
      ! the greatest common divisor of these two numbers.
      if(blocksizeSmall==-1) blocksizeSmall=gcd(blocksize*iter,norb/orbs%nkpts-ist+1)
      ! Determine how many bunches can be created with this new block size.
      iter2=(norb/orbs%nkpts-ist+1)/blocksizeSmall
      ! Now make a loop over all these blocks
      do iblock=1,iter2
         ! ist is the starting vector of the current bunch.
         ist=iter*blocksize+blocksizeSmall*(iblock-1)+1
         ! Now orthogonalize this bunch to all previous ones.
         do jblock=1,(blocksize*iter)/blocksizeSmall+iblock-1
            ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
            jst=blocksizeSmall*(jblock-1)+1
            call gramschmidtOverlap(iproc, nproc, norbtotp, blocksizeSmall, psiWTrans(1), &
               &   overlapPsiwTrans(1), newComm, ist, jst, norb/orbs%nkpts, orbs%nkpts, nspinor)
         end do
         ! Orthonormalize the current bunch of vectors.
         call choleskyOverlap(iproc, nproc, norbtotp, blocksizeSmall, psiWTrans(1),&
            &   overlapPsiWTrans(1), newComm, ist, norb/orbs%nkpts, orbs%nkpts, nspinor)
      end do
   end if remainingIf


   if(nproc>1) then
      ! Now the orthonormalization is done, so we can untranspose the vectors.
      ! First define the variables for the MPI call.
      ii=0 ; jj=0
      do i=nprocSt,nprocSt+nproc-1
         ! Now a given process sends a different amount if data to each process...
         sendcounts(i)=norbtotp*norbpArr(i)*nspinor
         sdispls(i)=ii
         ii=ii+norbtotp*norbpArr(i)*nspinor
      end do
      do i=nprocSt,nprocSt+nproc-1
         ! ... but receives the same amount from each process.
         recvcounts(i)=norbtotp*norbpArr(iproc)*nspinor
         rdispls(i)=jj
         jj=jj+norbtotp*norbpArr(iproc)*nspinor
      end do
      ! Now untranspose the vectors.
      call timing(iproc, 'Input_comput', 'OF')
      call timing(iproc, 'Input_commun', 'ON')
      call mpi_alltoallv(psiWTrans(1), sendcounts, sdispls, mpi_double_precision, psiW(1),&
         &   recvcounts, rdispls, mpi_double_precision, newComm, ierr)
      call timing(iproc, 'Input_commun', 'OF')
      call timing(iproc, 'Input_comput', 'ON')

      ! Now rearrange back the vectors.
      ii=0
      jj=1
      do i=nprocSt,nprocSt+nproc-1
         do iorb=0,norbp-1
            call vcopy(norbtotp*nspinor, psiW(jj), 1, psi(ii+iorb*norbtot*nspinor+1), 1)
            jj=jj+norbtotp*nspinor
         end do
         ii=ii+norbtotpArr(i)*nspinor
      end do
   else
      call vcopy(norbtot*norbp*nspinor, psiWTrans(1), 1, psi(1), 1)
   end if


   ! Deallocate all arrays
   i_all=-product(shape(psiW))*kind(psiW)
   deallocate(psiW, stat=i_stat)
   call memocc(i_stat, i_all, 'psiW', subname)

   i_all=-product(shape(overlapPsiW))*kind(overlapPsiW)
   deallocate(overlapPsiW, stat=i_stat)
   call memocc(i_stat, i_all, 'overlapPsiW', subname)

   i_all=-product(shape(psiWTrans))*kind(psiWTrans)
   deallocate(psiWTrans, stat=i_stat)
   call memocc(i_stat, i_all, 'psiWTrans', subname)

   i_all=-product(shape(overlapPsiWTrans))*kind(overlapPsiWTrans)
   deallocate(overlapPsiWTrans, stat=i_stat)
   call memocc(i_stat, i_all, 'overlapPsiWTrans', subname)

   i_all=-product(shape(sendcounts))*kind(sendcounts)
   deallocate(sendcounts, stat=i_stat)
   call memocc(i_stat, i_all, 'sendcounts', subname)

   i_all=-product(shape(recvcounts))*kind(recvcounts)
   deallocate(recvcounts, stat=i_stat)
   call memocc(i_stat, i_all, 'recvcounts', subname)

   i_all=-product(shape(sdispls))*kind(sdispls)
   deallocate(sdispls, stat=i_stat)
   call memocc(i_stat, i_all, 'sdidpls', subname)

   i_all=-product(shape(rdispls))*kind(rdispls)
   deallocate(rdispls, stat=i_stat)
   call memocc(i_stat, i_all, 'rdidpls', subname)

END SUBROUTINE orthonormalizePsi


!>  This subroutine orthogonalizes a given bunch of vectors (psi) to another bunch of equal size (psi). These other vectors
!!  are assumed to be orthonomal themselves. The orthonormalization is done in parallel, assuming that each process holds a 
!!  small portion of each vector.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param iproc       process ID
!!    @param nproc       total number of processes
!!    @param norbtot     length of the vectors
!!    @param overlapPsi  the overlap matrix applied to psi. The vectors are orthogonalized with respect to this overlap matrix.
!!    @param newComm     the MPI communicator for the processes handling this orthonormalization (not necessarily all MPI processes)
!!    @param istThis     starting index of the orbitals that shall be orthogonalized
!!    @param istOther    starting index of the orbitals to which the orbitals starting at istThis shall be orthogonalized
!!    @param norb        total number of vectors
!!    @param nkpts       number of k-points
!!    @param nspinor     real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!  Input/Output arguments:
!!    @param psi         the vectors that shall be orthogonalized.
subroutine gramschmidtOverlap(iproc, nproc, norbtot, blocksize, psi, overlapPsi, newComm,&
      &   istThis, istOther, norb, nkpts, nspinor)
   use module_base
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, nproc, norbtot, blocksize, newComm, istThis, istOther, norb, nkpts, nspinor
   real(kind=8),dimension(1:norbtot*nspinor,1:norb,nkpts),intent(in):: overlapPsi
   real(kind=8),dimension(1:norbtot*nspinor,1:norb,nkpts),intent(in out):: psi

   ! Local arguments
   integer:: ikpt, ist, ierr, i_stat, i_all
   real(kind=8),dimension(:,:),allocatable:: A
   real(kind=8),dimension(:,:,:),allocatable:: ovrlp
   character(len=*),parameter:: subname='gramschmidtOverlap'

   ! Allocate the matrix A which will hold some partial results.
   allocate(A(1:norbtot*nspinor,1:blocksize), stat=i_stat)
   call memocc(i_stat, A, 'A', subname)

   ! Allocate the matrix ovrlp which will save the overlap between the orbitals in psi and psi. For the parallel case we add another
   ! dimension: Each process writes its values to ovrlp(:,:,2) and then an mpi_allreduce will sum the contributions from all processes to ovrlp(:,:,1).
   if(nproc>1) then
      ist=2
   else
      ist=1
   end if
   allocate(ovrlp(1:blocksize*blocksize*nspinor,nkpts,ist), stat=i_stat)
   call memocc(i_stat, ovrlp, 'ovrlp', subname)

   ! Now calculate this overlap matrix: ovrlp=<psi|overlap|psi>.
   do ikpt=1,nkpts
      if(nspinor==1) then
         call dgemm('t', 'n', blocksize, blocksize, norbtot, 1.d0, psi(1,istOther,ikpt), norbtot,&
            &   overlapPsi(1,istThis,ikpt), norbtot, 0.d0, ovrlp(1,ikpt,ist), blocksize)
      else
         call zgemm('c', 'n', blocksize, blocksize, norbtot, (1.d0,0.d0), psi(1,istOther,ikpt), norbtot,&
            &   overlapPsi(1,istThis,ikpt), norbtot, (0.d0,0.d0), ovrlp(1,ikpt,ist), blocksize)
      end if
   end do
   if(nproc>1) then
      call timing(iproc, 'Input_comput', 'OF')
      call timing(iproc, 'Input_commun', 'ON')
      call mpi_allreduce (ovrlp(1,1,2), ovrlp(1,1,1), blocksize*blocksize*nspinor*nkpts,&
         &   mpi_double_precision, mpi_sum, newComm, ierr)
      call timing(iproc, 'Input_commun', 'OF')
      call timing(iproc, 'Input_comput', 'ON')
   end if

   do ikpt=1,nkpts
      ! Calculate matrix product psi*ovrlp=A. This will give the components that will be projected out of psi.
      ! We actually calculate -psi*ovrlp=-A, since this is better for further processing with daxpy.
      if(nspinor==1) then
         call dgemm('n', 'n', norbtot, blocksize, blocksize, -1.d0, psi(1,istOther,ikpt), &
            &   norbtot, ovrlp(1,ikpt,1), blocksize, 0.d0, A, norbtot)
      else
         call zgemm('n', 'n', norbtot, blocksize, blocksize, (-1.d0,0.d0), psi(1,istOther,ikpt),&
            &   norbtot, ovrlp(1,ikpt,1), blocksize, (0.d0,0.d0), A, norbtot)
      end if
      ! Now project out: psi=psi-A.
      ! Since we calculated -A, we have to put psi=psi+A and can use daxpy to perform psi=A+psi
      if(nspinor==1) then
         call daxpy(norbtot*blocksize, 1.d0, A, 1, psi(1,istThis,ikpt), 1)
      else
         call zaxpy(norbtot*blocksize, (1.d0,0.d0), A, 1, psi(1,istThis,ikpt), 1)
      end if
   end do

   ! Deallocate all arrays
   i_all=-product(shape(A))*kind(A)
   deallocate(A, stat=i_stat)
   call memocc(i_stat, i_all, 'A', subname)

   i_all=-product(shape(ovrlp))*kind(ovrlp)
   deallocate(ovrlp, stat=i_stat)
   call memocc(i_stat, i_all, 'ovrlp', subname)


END SUBROUTINE gramschmidtOverlap


!>  This subroutine orthonormalizes a given bunch of vectors psi. It first determines the overlap matrix S of the vectors
!!  and then calculates the Cholesky composition S=ovrlp*ovrlp^T. This matrix ovrlp is then inverted to get ovrlp^{-1} and the orthonormal
!!  vectors are finally given by psi=psi*ovrlp^{-1}.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param iproc       process ID
!!    @param nproc       total number of processes
!!    @param norbtot     length of the vectors
!!    @param blocksize   number of vectors
!!    @param overlapPsi  the overlap matrix applied to psi. The vectors are orthogonalized with respect to this overlap matrix.
!!    @param newComm     the communicator that handles the current processes
!!    @param istart      starting index of the orbitals that shall be orthogonalized
!!    @param norb        total number of vectors
!!    @param nkpts       number of k-points
!!    @param nspinor     real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!  Input/Output arguments:
!!    @param psi         the vectors that shall be orthonormalized
subroutine choleskyOverlap(iproc, nproc, norbtot, blocksize, psi, overlapPsi, &
      &   newComm, istart, norb, nkpts, nspinor)
   use module_base
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc,nproc,norbtot,blocksize,newComm, istart, norb, nkpts, nspinor
   real(kind=8),dimension(1:norbtot*nspinor,1:norb,nkpts),intent(in):: overlapPsi
   real(kind=8),dimension(1:norbtot*nspinor,1:norb,nkpts),intent(in out):: psi

   ! Local variables
   integer:: ikpt, ist, info, ierr, i_stat, i_all
   real(kind=8),allocatable,dimension(:,:,:):: ovrlp
   character(len=*),parameter:: subname='choleskyOverlap'


   ! Allocate the matrix ovrlp which will be the overlap between the orbitals in psiThis and psiOther. For the parallel case we add another
   ! dimension: Each process writes its values to S(:,:,2) and then an mpi_allreduce will sum the contributions from all processes to S(:,:,1).
   if(nproc>1) then 
      ist=2 
   else
      ist=1 
   end if
   allocate(ovrlp(1:blocksize*blocksize*nspinor,nkpts,ist), stat=i_stat)
   call memocc(i_stat, ovrlp, 'ovrlp', subname)
   ovrlp=0.d0

   ! Now calculate the overlap matrix ovrlp=<psi|overlap|psi>
   do ikpt=1,nkpts
      if(nspinor==1) then
         call dgemm('t', 'n', blocksize, blocksize, norbtot, 1.d0, psi(1,istart,ikpt), norbtot,&
            &   overlapPsi(1,istart,ikpt), norbtot, 0.d0, ovrlp(1,ikpt,ist), blocksize)
      else
         call zgemm('c', 'n', blocksize, blocksize, norbtot, (1.d0,0.d0), psi(1,istart,ikpt),&
            &   norbtot, overlapPsi(1,istart,ikpt), norbtot, (0.d0,0.d0), ovrlp(1,ikpt,ist), blocksize)
      end if
   end do
   if(nproc>1) then
      call timing(iproc, 'Input_comput', 'OF')
      call timing(iproc, 'Input_commun', 'ON')
      call mpi_allreduce (ovrlp(1,1,2), ovrlp(1,1,1), blocksize*blocksize*nkpts*nspinor,&
         &   mpi_double_precision, mpi_sum, newComm, ierr)
      call timing(iproc, 'Input_commun', 'OF')
      call timing(iproc, 'Input_comput', 'ON')
   end if


   do ikpt=1,nkpts
      ! Make a Cholesky factorization of ovrlp.
      if(nspinor==1) then
         call dpotrf('u', blocksize, ovrlp(1,ikpt,1), blocksize, info)
      else
         call zpotrf('l', blocksize, ovrlp(1,ikpt,1), blocksize, info)
      end if

      ! Invert the Cholesky matrix: ovrlp^{-1}.
      if(nspinor==1) then
         call dtrtri('u', 'n', blocksize, ovrlp(1,ikpt,1), blocksize, info)
      else
         call ztrtri('l', 'n', blocksize, ovrlp(1,ikpt,1), blocksize, info)
      end if

      ! Calculate the matrix product psi*ovrlp^{-1}=psi. This will give the orthonormal orbitals.
      !call dtrmm('r', 'u', 'n', 'n', norbtot, blocksize, 1.d0, ovrlp(1,1,1), blocksize, psi, norbtot)
      if(nspinor==1) then
         call dtrmm('r', 'u', 'n', 'n', norbtot, blocksize, 1.d0, ovrlp(1,ikpt,1),&
            &   blocksize, psi(1,istart,ikpt), norbtot)
      else
         call ztrmm('r', 'l', 'c', 'n', norbtot, blocksize, (1.d0,0.d0), ovrlp(1,ikpt,1),&
            &   blocksize, psi(1,istart,ikpt), norbtot)
      end if
   end do

   ! Deallocate the arrays.
   i_all=-product(shape(ovrlp))*kind(ovrlp)
   deallocate(ovrlp, stat=i_stat)
   call memocc(i_stat, i_all, 'ovrlp', subname)

END SUBROUTINE choleskyOverlap


!>  This function calculates the greatest common divisor of two numbers a and b.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param a    first number
!!    @param b    second number
!!  Output arguments:
!!    @param gcd  greatest common divisor of a and b
function gcd(a, b)
   implicit none

   ! Calling arguments
   integer,intent(in):: a,b
   integer:: gcd

   ! Local variables
   integer:: aSub, bSub, c

   aSub=a ; bSub=b
   ! Make sure that aSub>=bSub
   if(aSub<bSub) then
      c=aSub ; aSub=bSub ; bSub=c
   end if
   do
      c=mod(aSub,bSub)
      if(c==0) then
         gcd=bSub
         exit
      end if
      aSub=bSub
      bSub=c
   end do

END FUNCTION gcd


!>  This function determines a good block size for the Gram Schmidt/Cholesky orthonomalization procedure.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param iproc         process ID
!!    @param input         data type containing many parameters
!!    @param norb          number of vectors that have to be orthonormalized
!! Output arguments:
!!    @param getBlocksize  the good block size
function getBlocksize(orthpar, norb)
   !use module_base
   use module_types
   implicit none

   ! Calling arguments
   integer,intent(in):: norb
   integer:: getBlocksize
   type(orthon_data),intent(in):: orthpar

   ! Local variables
   integer:: remain, gcd, gcdCurr, gcdMax, i

   ! Initialze getBlocksize
   getBlocksize=-1


   if(orthpar%bsLow<orthpar%bsUp) then
      ! choose automatically
      !if(iproc==0) write(*,'(1x,a,2(i0,a))') 'Choose block size automatically between ',orthpar%bsLow,' and ',orthpar%bsUp
      !if(iproc==0) write(*,'(1x,a,i0,a)') '(or between 1 and ',orthpar%bsUp,' if norb<orthpar%bsUp, respectively).'
      getBlocksize=-1
      if(norb<orthpar%bsUp) then
         ! In this case we just take one block
         getBlocksize=norb
      else
         ! Else try to take a divisor of norb which is between orthpar%bsLow and orthpar%bsUp.
         do i=orthpar%bsLow,orthpar%bsUp
            if(mod(norb,i)==0) then
               getBlocksize=i
               exit
            end if
         end do
      end if
      ! If this is not possible, try to take a block size such that the greatest common divisor of the remainding block (i.e. the vectors 
      ! which do not fit into this block size) and the vectors that are already orthonormal is maximal. The block size is again limited
      ! to be bweteen orthpar%bsLow and orthpar%bsUp.
      if(getBlocksize==-1) then
         gcdMax=0
         do i=orthpar%bsLow,orthpar%bsUp
            remain=mod(norb,i)
            gcdCurr=gcd(norb-remain, remain)
            if(gcdCurr>gcdMax) then
               gcdMax=gcdCurr
               getBlocksize=i
            end if
         end do
      end if
   else if(orthpar%bsLow==orthpar%bsUp) then
      ! take the value specified by the user
      !if(iproc==0) write(*,'(1x,a)') 'Take blocksize specified by user.'
      if(orthpar%bsLow<=norb) then
         getBlocksize=orthpar%bsLow
      else
         !if(iproc==0) write(*,'(1x,a)') 'WARNING: specified blocksize is larger than the number of orbitals.'
         !if(iproc==0) write(*,'(1x,a,i0,a)') 'The blocksize is adjusted to ',norb,'.'
         getBlocksize=norb
      end if
   else
      !if(iproc==0) write(*,'(1x,a)') 'ERROR: invalid values of orthpar%bsLow and orthpar%bsUp. Change them in input.perf!'
      stop
   end if

END FUNCTION getBlocksize


!>  This subroutine initializes the random number generator for given iproc and ispin.
!!
!! Calling arguments:
!! =================
!!  Input arguments:
!!    @param iproc   process ID. For this subroutine it is just a number without meaning.
!!    @param ispin   spin up/down. For this subroutine it is just a number without meaning.
subroutine initRandomSeed(iproc, ispin)
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, ispin

   ! Local variables
   integer:: i, n
   integer,dimension(:),allocatable:: seed

   call random_seed(size=n)
   allocate(seed(n))
   i=0
   seed=37*(10*iproc+ispin)*(/(i-1, i=1,n)/)
   call random_seed(put=seed)

   deallocate(seed)
END SUBROUTINE initRandomSeed

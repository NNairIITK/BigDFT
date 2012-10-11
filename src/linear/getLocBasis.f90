!> @file
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine get_coeff(iproc,nproc,scf_mode,lzd,orbs,at,rxyz,denspot,&
    GPU, infoCoeff,ebs,nlpspd,proj,&
    SIC,tmb,fnrm,overlapmatrix,calculate_overlap_matrix,communicate_phi_for_lsumrho,&
    tmblarge, ham, ldiis_coeff)
use module_base
use module_types
use module_interfaces, exceptThisOne => get_coeff, exceptThisOneA => writeonewave
use Poisson_Solver
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc, scf_mode
type(local_zone_descriptors),intent(inout) :: lzd
type(orbitals_data),intent(inout) :: orbs
type(atoms_data),intent(in) :: at
real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(GPU_pointers),intent(inout) :: GPU
integer,intent(out) :: infoCoeff
real(kind=8),intent(out) :: ebs
real(kind=8),intent(inout) :: fnrm
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
type(SIC_data),intent(in) :: SIC
type(DFT_wavefunction),intent(inout) :: tmb
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: overlapmatrix
logical,intent(in):: calculate_overlap_matrix, communicate_phi_for_lsumrho
type(DFT_wavefunction),intent(inout):: tmblarge
!real(8),dimension(:),pointer,intent(inout):: lhphilarge
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in),optional:: ham
type(localizedDIISParameters),intent(inout),optional :: ldiis_coeff

! Local variables 
integer :: istat, iall, iorb, jorb, korb, info, iiorb, ierr
real(kind=8),dimension(:),allocatable :: eval, hpsit_c, hpsit_f
real(kind=8),dimension(:,:),allocatable :: ovrlp
real(kind=8),dimension(:,:,:),allocatable :: matrixElements
real(kind=8) :: tt
type(confpot_data),dimension(:),allocatable :: confdatarrtmp
type(energy_terms) :: energs
character(len=*),parameter :: subname='get_coeff'
!! integer :: ldim,istart,lwork,iiorb,ilr,ind2,ncnt
!! character(len=1) :: num

  ! Allocate the local arrays.  
  allocate(matrixElements(tmb%orbs%norb,tmb%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(eval(tmb%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  ! Calculate the overlap matrix if required.
  if(calculate_overlap_matrix) then
      if(.not.tmb%can_use_transposed) then
          if(.not.associated(tmb%psit_c)) then
              allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
              call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
          end if
          if(.not.associated(tmb%psit_f)) then
              allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
              call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
          end if
          call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psi, tmb%psit_c, tmb%psit_f, lzd)
          tmb%can_use_transposed=.true.
      end if
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, overlapmatrix)
  end if

  ! Post the p2p communications for the density. (must not be done in inputguess)
  !if(tmb%wfnmd%bs%communicate_phi_for_lsumrho) then
  if(communicate_phi_for_lsumrho) then
      call communicate_basis_for_density(iproc, nproc, lzd, tmb%orbs, tmb%psi, tmb%comsr)
  end if

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Calculate the Hamiltonian matrix if it is not already present.
  if(.not.present(ham)) then

      call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))
      call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
           tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)
      call test_p2p_communication(iproc, nproc, tmblarge%comgp)

      allocate(lzd%doHamAppl(lzd%nlr), stat=istat)
      call memocc(istat, lzd%doHamAppl, 'lzd%doHamAppl', subname)
      lzd%doHamAppl=.true.
      allocate(confdatarrtmp(tmb%orbs%norbp))
      call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

      if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%hpsi(1))
      allocate(tmblarge%lzd%doHamAppl(tmblarge%lzd%nlr), stat=istat)
      call memocc(istat, tmblarge%lzd%doHamAppl, 'tmblarge%lzd%doHamAppl', subname)
      tmblarge%lzd%doHamAppl=.true.

      call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
           proj,tmblarge%lzd,nlpspd,tmblarge%psi,tmblarge%hpsi,energs%eproj)
      ! only kinetic as waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call full_local_potential(iproc,nproc,tmblarge%orbs,tmblarge%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmblarge%comgp)
      !call wait_p2p_communication(iproc, nproc, tmblarge%comgp)
      ! only potential
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call timing(iproc,'glsynchham1','ON') !lr408t
      call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,tmblarge%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham1','OF') !lr408t
      deallocate(confdatarrtmp)

      iall=-product(shape(lzd%doHamAppl))*kind(lzd%doHamAppl)
      deallocate(lzd%doHamAppl, stat=istat)
      call memocc(istat, iall, 'lzd%doHamAppl', subname)

      iall=-product(shape(tmblarge%lzd%doHamAppl))*kind(tmblarge%lzd%doHamAppl)
      deallocate(tmblarge%lzd%doHamAppl, stat=istat)
      call memocc(istat, iall, 'tmblarge%lzd%doHamAppl', subname)

      iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
      deallocate(denspot%pot_work, stat=istat)
      call memocc(istat, iall, 'denspot%pot_work', subname)

      if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application done.'


      ! Calculate the matrix elements <phi|H|phi>.
      if(.not.tmblarge%can_use_transposed) then
          if(associated(tmblarge%psit_c)) then
              iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
              deallocate(tmblarge%psit_c, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_c', subname)
          end if
          if(associated(tmblarge%psit_f)) then
              iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
              deallocate(tmblarge%psit_f, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_f', subname)
          end if

          allocate(tmblarge%psit_c(tmblarge%collcom%ndimind_c), stat=istat)
          call memocc(istat, tmblarge%psit_c, 'tmblarge%psit_c', subname)
          allocate(tmblarge%psit_f(7*tmblarge%collcom%ndimind_f), stat=istat)
          call memocc(istat, tmblarge%psit_f, 'tmblarge%psit_f', subname)
          call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
               tmblarge%psi, tmblarge%psit_c, tmblarge%psit_f, tmblarge%lzd)
          tmblarge%can_use_transposed=.true.
      end if

      allocate(hpsit_c(tmblarge%collcom%ndimind_c))
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*tmblarge%collcom%ndimind_f))
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
           tmblarge%hpsi, hpsit_c, hpsit_f, tmblarge%lzd)
      call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%mad, tmblarge%collcom, &
           tmblarge%psit_c, hpsit_c, tmblarge%psit_f, hpsit_f, matrixElements)
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)

  else
      if(iproc==0) write(*,*) 'No Hamiltonian application required.'
      call dcopy(tmb%orbs%norb**2, ham(1,1), 1, matrixElements(1,1,1), 1)
  end if



  ! Diagonalize the Hamiltonian.
  if(scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      call dcopy(tmb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      call dcopy(tmb%orbs%norb**2, overlapmatrix(1,1),1 , ovrlp(1,1), 1)
      if(tmb%wfnmd%bpo%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, tmb%orbs, matrixElements(1,1,2), ovrlp, eval)
      else
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          call dsygv_parallel(iproc, nproc, tmb%wfnmd%bpo%blocksize_pdsyev, tmb%wfnmd%bpo%nproc_pdsyev, &
               bigdft_mpi%mpi_comm, 1, 'v', 'l',tmb%orbs%norb, &
               matrixElements(1,1,2), tmb%orbs%norb, ovrlp, tmb%orbs%norb, eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'

      do iorb=1,orbs%norb
          call dcopy(tmb%orbs%norb, matrixElements(1,iorb,2), 1, tmb%wfnmd%coeff(1,iorb), 1)
      end do
      infoCoeff=0

      ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
      if(iproc==0) then
          write(*,'(1x,a)') '-------------------------------------------------'
          write(*,'(1x,a)') 'some selected eigenvalues:'
          do iorb=max(orbs%norb-8,1),min(orbs%norb+8,tmb%orbs%norb)
              if(iorb==orbs%norb) then
                  write(*,'(3x,a,i0,a,es12.5,a)') 'eval(',iorb,')= ',eval(iorb),'  <-- last occupied orbital'
              else if(iorb==orbs%norb+1) then
                  write(*,'(3x,a,i0,a,es12.5,a)') 'eval(',iorb,')= ',eval(iorb),'  <-- first virtual orbital'
              else
                  write(*,'(3x,a,i0,a,es12.5)') 'eval(',iorb,')= ',eval(iorb)
              end if
          end do
          write(*,'(1x,a)') '-------------------------------------------------'
      end if

      ! keep the eigenvalues for the preconditioning - instead should take h_alpha,alpha for both cases
      call vcopy(tmb%orbs%norb, eval(1), 1, tmb%orbs%eval(1), 1)
      call vcopy(tmb%orbs%norb, eval(1), 1, tmblarge%orbs%eval(1), 1)
      ! instead just use -0.5 everywhere
      !tmb%orbs%eval(:) = -0.5_dp
      !tmblarge%orbs%eval(:) = -0.5_dp
  else
      if(.not.present(ldiis_coeff)) stop 'ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION'
      call optimize_coeffs(iproc, nproc, orbs, matrixElements(1,1,1), overlapmatrix, tmb, ldiis_coeff, fnrm)
  end if

  call calculate_density_kernel(iproc, nproc, .true., tmb%wfnmd%ld_coeff, orbs, tmb%orbs, &
       tmb%wfnmd%coeff, tmb%wfnmd%density_kernel, ovrlp)

  ! Calculate the band structure energy with matrixElements instead of wfnmd%coeff due to the problem mentioned
  ! above (wrong size of wfnmd%coeff)
  ebs=0.d0
  do jorb=1,tmb%orbs%norb
      do korb=1,jorb
          tt = tmb%wfnmd%density_kernel(korb,jorb)*matrixElements(korb,jorb,1)
          if(korb/=jorb) tt=2.d0*tt
          ebs = ebs + tt
      end do
  end do

  ! Calculate the KS eigenvalues - needed for Pulay
  call to_zero(orbs%norb, orbs%eval(1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,tmb%orbs%norb
          do korb=1,tmb%orbs%norb
              orbs%eval(iiorb) = orbs%eval(iiorb) + &
                                 tmb%wfnmd%coeff(jorb,iiorb)*tmb%wfnmd%coeff(korb,iiorb)*matrixElements(jorb,korb,1)
          end do
      end do
  end do
  call mpiallred(orbs%eval(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)


  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,*) orbs%eval(iorb), tmblarge%orbs%eval(iorb)
  !!    end do
  !!end if

  ! If closed shell multiply by two.
  !if(orbs%nspin==1) ebs=2.d0*ebs

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)
  
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)


end subroutine get_coeff



subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,&
    denspot,GPU,trH,fnrm, infoBasisFunctions,nlpspd,proj,ldiis,&
    SIC, tmb, tmblarge, energs_base, ham)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are obtained by adding a
!   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
!   determined by minimizing the trace until the gradient norm is below the convergence criterion.
use module_base
use module_types
use module_interfaces, except_this_one => getLocalizedBasis, except_this_one_A => writeonewave
!  use Poisson_Solver
!use allocModule
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
integer,intent(out) :: infoBasisFunctions
type(atoms_data), intent(in) :: at
type(orbitals_data) :: orbs
real(kind=8),dimension(3,at%nat) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(GPU_pointers), intent(inout) :: GPU
real(kind=8),intent(out) :: trH, fnrm
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
type(localizedDIISParameters),intent(inout) :: ldiis
type(DFT_wavefunction),target,intent(inout) :: tmb
type(SIC_data) :: SIC !<parameters for the SIC methods
type(DFT_wavefunction),target,intent(inout) :: tmblarge
!real(kind=8),dimension(:),pointer,intent(inout) :: lhphilarge2
type(energy_terms),intent(in) :: energs_base
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ham

! Local variables
real(kind=8) :: trHold, fnrmMax, meanAlpha, ediff, noise, alpha_max, delta_energy, delta_energy_prev
integer :: iorb, istat,ierr,it,iall,nsatur, it_tot, ncount, jorb, iiorb
real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp
real(kind=8),dimension(:,:),allocatable :: ovrlp, coeff_old, delta_ham, delta_kernel, lagmat
logical :: energy_increased, overlap_calculated
character(len=*),parameter :: subname='getLocalizedBasis'
real(kind=8),dimension(:),pointer :: lhphi, lhphiold, lphiold, hpsit_c, hpsit_f
type(energy_terms) :: energs
real(8),dimension(2):: reducearr
real(8),save:: trH_old



  ! Allocate all local arrays.
  call allocateLocalArrays()

  ! setting lhphiold to zero for calculate_energy_and_gradient_linear - why is this needed?
  call to_zero(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),lhphiold(1))

  call timing(iproc,'getlocbasinit','ON') !lr408t
  tmb%can_use_transposed=.false.
  if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS


  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
  trHold=1.d100

  nsatur=0
  delta_energy_prev=1.d100
 
  call timing(iproc,'getlocbasinit','OF') !lr408t

  !if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) trH_old=0.d0
  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) trH_old=0.d0
  overlap_calculated=.false.
  it=0
  it_tot=0
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))
  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)
  call test_p2p_communication(iproc, nproc, tmblarge%comgp)


  !iterLoop: do it=1,tmb%wfnmd%bs%nit_basis_optimization
  iterLoop: do
      it=it+1
      it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
      it_tot=it_tot+1


      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif


      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%hpsi(1))
      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, &
           tmb%psi, tmblarge%psi)


      allocate(tmblarge%lzd%doHamAppl(tmblarge%lzd%nlr), stat=istat)
      call memocc(istat, tmblarge%lzd%doHamAppl, 'tmblarge%lzd%doHamAppl', subname)
      tmblarge%lzd%doHamAppl=.true.
      call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
           proj,tmblarge%lzd,nlpspd,tmblarge%psi,tmblarge%hpsi,energs%eproj)
      ! only kinetic because waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,tmblarge%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call full_local_potential(iproc,nproc,tmblarge%orbs,tmblarge%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmblarge%comgp)
      !call wait_p2p_communication(iproc, nproc, tmblarge%comgp)
      ! only potential
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,tmblarge%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call timing(iproc,'glsynchham2','ON') !lr408t
      call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,tmblarge%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF') !lr408t


  iall=-product(shape(tmblarge%lzd%doHamAppl))*kind(tmblarge%lzd%doHamAppl)
  deallocate(tmblarge%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'tmblarge%lzd%doHamAppl', subname)

!!!DEBUG
!!if (iproc==0) then
!!   write(*,'(1x,a,4(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
!!   2*energs%ekin,2*energs%epot,2*energs%eproj,2*energs%ekin+2*energs%epot+2*energs%eproj
!!   !!write(*,'(1x,a,4(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
!!   !!2*av_h_sym_diff1,2*av_h_sym_diff2,2*av_h_sym_diff3,2*av_h_sym_diff1+2*av_h_sym_diff2+2*av_h_sym_diff3
!!endif
!!!END DEBUG


  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') ' Orthoconstraint... '
      end if

      call copy_basis_specifications(tmb%wfnmd%bs, tmblarge%wfnmd%bs, subname)
      call copy_orthon_data(tmb%orthpar, tmblarge%orthpar, subname)

      call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, tmblarge%hpsi, hpsit_c, hpsit_f, tmblarge%lzd)

      ncount=sum(tmblarge%collcom%nrecvcounts_c)
      if(ncount>0) call dcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*sum(tmblarge%collcom%nrecvcounts_f)
      if(ncount>0) call dcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)


      call calculate_energy_and_gradient_linear(iproc, nproc, it, &
           tmb%wfnmd%density_kernel, &
           ldiis, &
           fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, &
           meanAlpha, alpha_max, energy_increased, &
           tmb, lhphi, lhphiold, &
           tmblarge, tmblarge%hpsi, overlap_calculated, ovrlp, lagmat, energs_base, hpsit_c, hpsit_f)

      !! EXERIMENTAL #######################################################
      delta_energy=0.d0
      if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              do jorb=1,tmb%orbs%norb
                  delta_energy = delta_energy - alpha(iorb)*lagmat(jorb,iiorb)*tmb%wfnmd%density_kernel(jorb,iiorb)
              end do
          end do
      else if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              delta_energy = delta_energy - alpha(iorb)*lagmat(iiorb,iiorb)
          end do
      end if
      call mpiallred(delta_energy, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      !!if (iproc==0) write(*,*) 'delta_energy',delta_energy
      !! END EXERIMENTAL ####################################################

      if (energy_increased) then
          tmblarge%can_use_transposed=.false.
          call dcopy(tmb%orbs%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          ! Recalculate the kernel with the old coefficients
          call dcopy(orbs%norb*tmb%orbs%norb, coeff_old(1,1), 1, tmb%wfnmd%coeff(1,1), 1)
          call calculate_density_kernel(iproc, nproc, .true., tmb%wfnmd%ld_coeff, orbs, tmb%orbs, &
               tmb%wfnmd%coeff, tmb%wfnmd%density_kernel, ovrlp)
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
          if(associated(tmblarge%psit_c)) then
              iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
              deallocate(tmblarge%psit_c, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_c', subname)
          end if
          if(associated(tmblarge%psit_f)) then
              iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
              deallocate(tmblarge%psit_f, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_f', subname)
              tmblarge%can_use_transposed=.false.
          end if
          if(iproc==0) write(*,*) 'it_tot',it_tot
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (iproc==0) write(*,'(1x,a,i6,2es15.7,f17.10,2es13.4)') 'iter, fnrm, fnrmMax, ebs, diff, noise level', &
          it, fnrm, fnrmMax, trH, ediff,noise
          if(it_tot<3*tmb%wfnmd%bs%nit_basis_optimization) cycle
      end if 


      ediff=trH-trH_old
      !!!noise=tmb%wfnmd%bs%gnrm_mult*fnrm*tmb%orbs%norb
      !!noise=tmb%wfnmd%bs%gnrm_mult*fnrm*abs(trH)
      !!if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY )then
      !!    if (ediff<0.d0 .and. abs(ediff) < noise) then
      !!        if(iproc==0) write(*,'(a)') 'target function seems to saturate, increase nsatur...'
      !!        nsatur=nsatur+1
      !!    else if (abs(ediff) < noise .and. meanAlpha<.1d0) then
      !!        if(iproc==0) write(*,'(a)') 'target function increases (but smaller than noise) and step size is small.'
      !!        if(iproc==0) write(*,'(a)') 'Consider convergence.'
      !!        nsatur=tmb%wfnmd%bs%nsatur_inner
      !!    else
      !!        nsatur=0
      !!    end if
      !!end if


      ! Write some informations to the screen.
      if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
          write(*,'(1x,a,i6,2es15.7,f17.10,es13.4)') 'iter, fnrm, fnrmMax, trace, diff', &
          it, fnrm, fnrmMax, trH, ediff
      if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
          write(*,'(1x,a,i6,2es15.7,f17.10,es13.4)') 'iter, fnrm, fnrmMax, ebs, diff', &
          it, fnrm, fnrmMax, trH, ediff
      if(it>=tmb%wfnmd%bs%nit_basis_optimization .or. it_tot>=3*tmb%wfnmd%bs%nit_basis_optimization) then
          !!if(ediff<0.d0 .and. ediff>1.d-10*delta_energy_prev) then
          !!    if(iproc==0) write(*,*) 'CONVERGED'
          !!    infoBasisFunctions=it
          !!else if(nsatur>=tmb%wfnmd%bs%nsatur_inner) then
          !!    if(iproc==0) then
          !!        write(*,'(1x,a,i0,a,2es15.7,f15.7)') 'converged in ', it, ' iterations.'
          !!        if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
          !!            write (*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
          !!        if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
          !!            write (*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
          !!    end if
          !!    infoBasisFunctions=it
          if(it>=tmb%wfnmd%bs%nit_basis_optimization .and. .not.energy_increased) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=0
          else if(it_tot>=3*tmb%wfnmd%bs%nit_basis_optimization) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: there seem to be some problems, exiting now...'
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          if (infoBasisFunctions>=0) then
              ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
              ! iteration of get_coeff.
              call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%mad, tmblarge%collcom, &
                   tmblarge%psit_c, hpsit_c_tmp, tmblarge%psit_f, hpsit_f_tmp, ham)
          end if

          exit iterLoop
      end if
      trH_old=trH
      delta_energy_prev=delta_energy


      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lhphi, lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS)
      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmblarge%can_use_transposed) then
          iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
          deallocate(tmblarge%psit_c, stat=istat)
          call memocc(istat, iall, 'tmblarge%psit_c', subname)
          iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
          deallocate(tmblarge%psit_f, stat=istat)
          call memocc(istat, iall, 'tmblarge%psit_f', subname)
          tmblarge%can_use_transposed=.false.
      end if

      ! Copy the coefficients to coeff_ols. The coefficients will be modified in reconstruct_kernel.
      call dcopy(orbs%norb*tmb%orbs%norb, tmb%wfnmd%coeff(1,1), 1, coeff_old(1,1), 1)

      if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
          call reconstruct_kernel(iproc, nproc, 1, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
               orbs, tmb, ovrlp, overlap_calculated, tmb%wfnmd%density_kernel)
      end if
      if(iproc==0) then
          write(*,'(a)') 'done.'
      end if


  end do iterLoop

  ! Deallocate potential
  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)


  ! Keep the values for the next iteration
  reducearr(1)=0.d0
  reducearr(2)=0.d0
  do iorb=1,tmb%orbs%norbp
      reducearr(1)=reducearr(1)+alpha(iorb)
      reducearr(2)=reducearr(2)+alphaDIIS(iorb)
  end do
  call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
  reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)

  ldiis%alphaSD=reducearr(1)
  ldiis%alphaDIIS=reducearr(2)


  ! Deallocate all local arrays.
  call deallocateLocalArrays()

contains




    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
      allocate(alpha(tmb%orbs%norbp), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)

      allocate(alphaDIIS(tmb%orbs%norbp), stat=istat)
      call memocc(istat, alphaDIIS, 'alphaDIIS', subname)

      allocate(fnrmOldArr(tmb%orbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)

      allocate(lhphi(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp)), stat=istat)
      call memocc(istat, lhphi, 'lhphi', subname)
    
      allocate(lhphiold(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp)), stat=istat)
      call memocc(istat, lhphiold, 'lhphiold', subname)

      allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, ovrlp, 'ovrlp', subname)

      allocate(lphiold(size(tmb%psi)), stat=istat)
      call memocc(istat, lphiold, 'lphiold', subname)

      allocate(hpsit_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)

      allocate(hpsit_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)

      allocate(hpsit_c_tmp(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c_tmp, 'hpsit_c_tmp', subname)

      allocate(hpsit_f_tmp(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f_tmp, 'hpsit_f_tmp', subname)

      allocate(coeff_old(tmb%orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, coeff_old, 'coeff_old', subname)

      allocate(delta_ham(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, delta_ham, 'delta_ham', subname)

      allocate(delta_kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, delta_kernel, 'delta_kernel', subname)

      allocate(lagmat(tmblarge%orbs%norb,tmblarge%orbs%norb), stat=istat)
      call memocc(istat, lagmat, 'lagmat', subname)

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(alphaDIIS))*kind(alphaDIIS)
      deallocate(alphaDIIS, stat=istat)
      call memocc(istat, iall, 'alphaDIIS', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)

      iall=-product(shape(lhphi))*kind(lhphi)
      deallocate(lhphi, stat=istat)
      call memocc(istat, iall, 'lhphi', subname)

      iall=-product(shape(lhphiold))*kind(lhphiold)
      deallocate(lhphiold, stat=istat)
      call memocc(istat, iall, 'lhphiold', subname)

      iall=-product(shape(lphiold))*kind(lphiold)
      deallocate(lphiold, stat=istat)
      call memocc(istat, iall, 'lphiold', subname)

      iall=-product(shape(ovrlp))*kind(ovrlp)
      deallocate(ovrlp, stat=istat)
      call memocc(istat, iall, 'ovrlp', subname)

      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)

      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)

      iall=-product(shape(hpsit_c_tmp))*kind(hpsit_c_tmp)
      deallocate(hpsit_c_tmp, stat=istat)
      call memocc(istat, iall, 'hpsit_c_tmp', subname)

      iall=-product(shape(hpsit_f_tmp))*kind(hpsit_f_tmp)
      deallocate(hpsit_f_tmp, stat=istat)
      call memocc(istat, iall, 'hpsit_f_tmp', subname)

      iall=-product(shape(coeff_old))*kind(coeff_old)
      deallocate(coeff_old, stat=istat)
      call memocc(istat, iall, 'coeff_old', subname)

      iall=-product(shape(delta_ham))*kind(delta_ham)
      deallocate(delta_ham, stat=istat)
      call memocc(istat, iall, 'delta_ham', subname)

      iall=-product(shape(delta_kernel))*kind(delta_kernel)
      deallocate(delta_kernel, stat=istat)
      call memocc(istat, iall, 'delta_kernel', subname)

      iall=-product(shape(lagmat))*kind(lagmat)
      deallocate(lagmat, stat=istat)
      call memocc(istat, iall, 'lagmat', subname)


    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine improveOrbitals(iproc, nproc, it, tmb, ldiis, lhphi, alpha)
  use module_base
  use module_types
  use module_interfaces, except_this_one => improveOrbitals
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  type(DFT_wavefunction),intent(inout) :: tmb
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmb%wfnmd%nphi),intent(in) :: lhphi
  real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
  
  ! Local variables
  integer :: istart, iorb, iiorb, ilr, ncount
  
  if (ldiis%isx > 0) then
      ldiis%mis=mod(ldiis%is,ldiis%isx)+1
      ldiis%is=ldiis%is+1
  end if
  
  ! steepest descent
  if(ldiis%isx==0) then
      call timing(iproc,'optimize_SD   ','ON')
      istart=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call daxpy(ncount, -alpha(iorb), lhphi(istart), 1, tmb%psi(istart), 1)
          istart=istart+ncount
      end do
      call timing(iproc,'optimize_SD   ','OF')
  else
      ! DIIS
      if(ldiis%alphaDIIS/=1.d0) then
          call dscal(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp), ldiis%alphaDIIS, lhphi, 1)
      end if
      call optimizeDIIS(iproc, nproc, tmb%orbs, tmb%orbs, tmb%lzd, lhphi, tmb%psi, ldiis, it)
  end if
end subroutine improveOrbitals



subroutine my_geocode_buffers(geocode,nl1,nl2,nl3)
  implicit none
  integer, intent(out) :: nl1,nl2,nl3
  character(len=1), intent(in) :: geocode
  !local variables
  logical :: perx,pery,perz
  integer :: nr1,nr2,nr3

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nl1,nr1)
  call ext_buffers(pery,nl2,nr2)
  call ext_buffers(perz,nl3,nr3)

end subroutine my_geocode_buffers







subroutine diagonalizeHamiltonian2(iproc, orbs, HamSmall, ovrlp, eval)
!
! Purpose:
! ========
!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!   the same result. This is done by requiring that the first entry of each vector
!   is positive.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc     process ID
!     nproc     number of MPI processes
!     orbs      type describing the physical orbitals psi
!   Input / Putput arguments
!     HamSmall  on input: the Hamiltonian
!               on exit: the eigenvectors
!   Output arguments
!     eval      the associated eigenvalues 
!
use module_base
use module_types
implicit none

! Calling arguments
integer, intent(in) :: iproc
type(orbitals_data), intent(inout) :: orbs
real(kind=8),dimension(orbs%norb, orbs%norb),intent(inout) :: HamSmall
real(kind=8),dimension(orbs%norb, orbs%norb),intent(in) :: ovrlp
real(kind=8),dimension(orbs%norb),intent(out) :: eval

! Local variables
integer :: lwork, info, istat, iall
real(kind=8),dimension(:),allocatable :: work
character(len=*),parameter :: subname='diagonalizeHamiltonian'

  call timing(iproc,'diagonal_seq  ','ON')

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=int(work(1))

  ! Deallocate the work array and reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)

  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2






subroutine small_to_large_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, phismall, philarge)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbssmall, orbslarge
  real(kind=8),dimension(orbssmall%npsidim_orbs),intent(in) :: phismall
  real(kind=8),dimension(orbslarge%npsidim_orbs),intent(out) :: philarge
  
  ! Local variables
  integer :: ists, istl, iorb, ilr, ilrlarge, sdim, ldim, nspin
       call timing(iproc,'small2large','ON') ! lr408t 
  ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
  !call to_zero(orbslarge%npsidim_orbs, philarge(1))
  ists=1
  istl=1
  do iorb=1,orbslarge%norbp
      ilr = orbssmall%inWhichLocreg(orbssmall%isorb+iorb)
      ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
      sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      ldim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
      nspin=1 !this must be modified later
      call Lpsi_to_global2(iproc, nproc, sdim, ldim, orbssmall%norb, orbssmall%nspinor, nspin, lzdlarge%llr(ilrlarge), &
           lzdsmall%llr(ilr), phismall(ists), philarge(istl))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
  end do
  if(orbssmall%norbp>0 .and. ists/=orbssmall%npsidim_orbs+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= orbssmall%npsidim_orbs+1=',orbssmall%npsidim_orbs+1
      stop
  end if
  if(orbslarge%norbp>0 .and. istl/=orbslarge%npsidim_orbs+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istk /= orbslarge%npsidim_orbs+1=',orbslarge%npsidim_orbs+1
      stop
  end if
       call timing(iproc,'small2large','OF') ! lr408t 
end subroutine small_to_large_locreg


subroutine large_to_small_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, philarge, phismall)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbssmall, orbslarge
  real(kind=8),dimension(orbslarge%npsidim_orbs),intent(in) :: philarge
  real(kind=8),dimension(orbssmall%npsidim_orbs),intent(out) :: phismall
  
  ! Local variables
  integer :: istl, ists, ilr, ilrlarge, ldim, gdim, iorb
       call timing(iproc,'large2small','ON') ! lr408t   
  ! Transform back to small locreg
  ! No need to this array to zero, since all values will be filled with a value during the copy.
  !!call to_zero(orbssmall%npsidim_orbs, phismall(1))
  ists=1
  istl=1
  do iorb=1,orbssmall%norbp
      ilr = orbssmall%inWhichLocreg(orbssmall%isorb+iorb)
      ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
      ldim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdsmall%llr(ilr), lzdlarge%llr(ilrlarge), &
           philarge(istl:istl+gdim-1), phismall(ists:ists+ldim-1))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
  end do

  if(orbssmall%norbp>0 .and. ists/=orbssmall%npsidim_orbs+1) stop 'ists/=orbssmall%npsidim_orbs+1'
  if(orbslarge%norbp>0 .and. istl/=orbslarge%npsidim_orbs+1) stop 'istl/=orbslarge%npsidim_orbs+1'
       call timing(iproc,'large2small','OF') ! lr408t 
end subroutine large_to_small_locreg









subroutine communicate_basis_for_density(iproc, nproc, lzd, llborbs, lphi, comsr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => communicate_basis_for_density
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: llborbs
  real(kind=8),dimension(llborbs%npsidim_orbs),intent(in) :: lphi
  type(p2pComms),intent(inout) :: comsr
  
  ! Local variables
  integer :: ist, istr, iorb, iiorb, ilr
  type(workarr_sumrho) :: w

  call timing(iproc,'commbasis4dens','ON') !lr408t

  ! Allocate the communication buffers for the calculation of the charge density.
  !call allocateCommunicationbufferSumrho(iproc, comsr, subname)
  ! Transform all orbitals to real space.
  ist=1
  istr=1
  do iorb=1,llborbs%norbp
      iiorb=llborbs%isorb+iorb
      ilr=llborbs%inWhichLocreg(iiorb)
      call initialize_work_arrays_sumrho(lzd%Llr(ilr), w)
      call daub_to_isf(lzd%Llr(ilr), w, lphi(ist), comsr%sendBuf(istr))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
  end do
  if(istr/=comsr%nsendBuf+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=comsr%nsendBuf+1'
      stop
  end if
  
  ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
  ! to point communication, the program will continue immediately. The messages will be gathered
  ! in the subroutine sumrhoForLocalizedBasis2.
  !!call postCommunicationSumrho2(iproc, nproc, comsr, comsr%sendBuf, comsr%recvBuf)
  call post_p2p_communication(iproc, nproc, comsr%nsendbuf, comsr%sendbuf, comsr%nrecvbuf, comsr%recvbuf, comsr)

  call timing(iproc,'commbasis4dens','OF') !lr408t

end subroutine communicate_basis_for_density




subroutine DIISorSD(iproc, nproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it
  real(kind=8),intent(in) :: trH
  type(DFT_wavefunction),intent(inout) :: tmbopt
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  real(kind=8),dimension(tmbopt%wfnmd%nphi),intent(out) :: lphioldopt
  
  ! Local variables
  integer :: idsx, ii, offset, istdest, iorb, iiorb, ilr, ncount, istsource
  
  !
  ! Purpose:
  ! ========
  !   This subroutine decides whether one should use DIIS or variable step size
  !   steepest descent to improve the orbitals. In the beginning we start with DIIS
  !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
  !   steepest descent. If the steepest descent iterations are successful, we switch
  !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
  !   history length is limited to be larger or equal than lin%DIISHistMin.
  !



  ! If we swicthed to SD in the previous iteration, reset this flag.
  if(ldiis%switchSD) ldiis%switchSD=.false.
  !if(iproc==0) write(*,'(a,2es15.6,l5)') 'trH, ldiis%trmin, ldiis%resetDIIS', trH, ldiis%trmin, ldiis%resetDIIS

  ! Now come some checks whether the trace is descreasing or not. This further decides
  ! whether we should use DIIS or SD.

  ! Determine wheter the trace is decreasing (as it should) or increasing.
  ! This is done by comparing the current value with diisLIN%energy_min, which is
  ! the minimal value of the trace so far.
  !if(iproc==0) write(*,*) 'trH, ldiis%trmin', trH, ldiis%trmin
  if(trH<=ldiis%trmin .and. .not.ldiis%resetDIIS) then
      ! Everything ok
      ldiis%trmin=trH
      ldiis%switchSD=.false.
      ldiis%itBest=it
      ldiis%icountSDSatur=ldiis%icountSDSatur+1
      ldiis%icountDIISFailureCons=0
      !if(iproc==0) write(*,*) 'everything ok, copy last psi...'
      call dcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)

      ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
      ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
      if(ldiis%icountSDSatur>=10 .and. ldiis%isx==0 .or. ldiis%immediateSwitchToSD) then
          ldiis%icountSwitch=ldiis%icountSwitch+1
          idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-ldiis%icountSwitch)
          if(idsx>0) then
              if(iproc==0) write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
              ldiis%icountSDSatur=0
              ldiis%icountSwitch=0
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%is=0
              ldiis%switchSD=.false.
              ldiis%trmin=1.d100
              ldiis%trold=1.d100
              alpha=ldiis%alphaSD
              alphaDIIS=ldiis%alphaDIIS
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%immediateSwitchToSD=.false.
          end if
      end if
  else
      ! The trace is growing.
      ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
      ! total failures or after 2 consecutive failures.
      ldiis%icountDIISFailureCons=ldiis%icountDIISFailureCons+1
      ldiis%icountDIISFailureTot=ldiis%icountDIISFailureTot+1
      ldiis%icountSDSatur=0
      if((ldiis%icountDIISFailureCons>=2 .or. ldiis%icountDIISFailureTot>=3 .or. ldiis%resetDIIS) .and. ldiis%isx>0) then
          ! Switch back to SD.
          alpha=ldiis%alphaSD
          if(iproc==0) then
              if(ldiis%icountDIISFailureCons>=2) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
              if(ldiis%icountDIISFailureTot>=3) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  ldiis%icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
              if(ldiis%resetDIIS) write(*,'(1x,a)') 'reset DIIS due to flag'
          end if
          if(ldiis%resetDIIS) then
              ldiis%resetDIIS=.false.
              ldiis%immediateSwitchToSD=.true.
              ldiis%trmin=1.d100
          end if
          ! Otherwise there could be problems due to the orthonormalization (which sligtly increases 
          ! value of the target function)
          ldiis%trmin=1.d100
          ! Try to get back the orbitals of the best iteration. This is possible if
          ! these orbitals are still present in the DIIS history.
          if(it-ldiis%itBest<ldiis%isx) then
             if(iproc==0) then
                 if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                     ldiis%itBest, ' which are the best so far.'
             end if
             ii=modulo(ldiis%mis-(it-ldiis%itBest),ldiis%mis)
             offset=0
             istdest=1
             !if(iproc==0) write(*,*) 'copy DIIS history psi...'
             do iorb=1,tmbopt%orbs%norbp
                 iiorb=tmbopt%orbs%isorb+iorb
                 ilr=tmbopt%orbs%inWhichLocreg(iiorb)
                 ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
                 istsource=offset+ii*ncount+1
                 call dcopy(ncount, ldiis%phiHist(istsource), 1, tmbopt%psi(istdest), 1)
                 call dcopy(ncount, ldiis%phiHist(istsource), 1, lphioldopt(istdest), 1)
                 offset=offset+ldiis%isx*ncount
                 istdest=istdest+ncount
             end do
         else
             !if(iproc==0) write(*,*) 'copy last psi...'
             call dcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)
         end if
         ldiis%isx=0
         ldiis%switchSD=.true.
      end if
      ! to indicate that no orthonormalization is required... (CHECK THIS!)
      if(ldiis%isx==0) ldiis%switchSD=.true. 
  end if

end subroutine DIISorSD






subroutine reconstruct_kernel(iproc, nproc, iorder, blocksize_dsyev, blocksize_pdgemm, orbs, tmb, &
           ovrlp_tmb, overlap_calculated, kernel)
  use module_base
  use module_types
  use module_interfaces, except_this_one => reconstruct_kernel
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, iorder, blocksize_dsyev, blocksize_pdgemm
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ovrlp_tmb
  logical,intent(out):: overlap_calculated
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: kernel

  ! Local variables
  integer:: istat, ierr, iall
  real(8),dimension(:,:),allocatable:: coeff_tmp, ovrlp_tmp, ovrlp_coeff
  character(len=*),parameter:: subname='reconstruct_kernel'

  call timing(iproc,'renormCoefComp','ON')

  allocate(coeff_tmp(tmb%orbs%norb,max(orbs%norbp,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
  allocate(ovrlp_tmp(orbs%norb,max(orbs%norbp,1)), stat=istat)
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  if(iproc==0) then
      write(*,'(a)',advance='no') 'coeff renormalization...'
  end if

  ! Calculate the overlap matrix between the TMBs.
  if(.not.tmb%can_use_transposed) then
      if(associated(tmb%psit_c)) then
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
      end if
      if(associated(tmb%psit_f)) then
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if
      allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
      allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
      call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      tmb%can_use_transposed=.true.
  end if
  call timing(iproc,'renormCoefComp','OF')
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, tmb%collcom, &
       tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, ovrlp_tmb)
  call timing(iproc,'renormCoefComp','ON')
  overlap_calculated=.true.

  !write(*,*) 'iproc, orbs%isorb', iproc, orbs%isorb
  ! Calculate the overlap matrix among the coefficients with resct to ovrlp_tmb.
  if (orbs%norbp>0 )then
      call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0, ovrlp_tmb(1,1), tmb%orbs%norb, &
           tmb%wfnmd%coeff(1,orbs%isorb+1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  end if
  !!do iorb=1,orbs%norbp
  !!    do jorb=1,orbs%norb
  !!        ovrlp_tmp(jorb,iorb)=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb), 1, coeff_tmp(1,iorb), 1)
  !!        !!if(iproc==0) write(210,*) ovrlp_tmp(jorb,iorb)
  !!    end do
  !!end do
  if (orbs%norbp>0 )then
      call dgemm('t', 'n', orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0,  tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
           coeff_tmp(1,1), tmb%orbs%norb, 0.d0, ovrlp_tmp(1,1), orbs%norb)
  end if


  call timing(iproc,'renormCoefComp','OF')
  call timing(iproc,'renormCoefComm','ON')
  ! Gather together the complete matrix
  if (nproc>1) then
     call mpi_allgatherv(ovrlp_tmp(1,1), orbs%norb*orbs%norbp, mpi_double_precision, ovrlp_coeff(1,1), &
          orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(orbs%norb*orbs%norb,ovrlp_tmp(1,1),1,ovrlp_coeff(1,1),1)
  end if
  call timing(iproc,'renormCoefComm','OF')
  call timing(iproc,'renormCoefComp','ON')

  ! Recalculate the kernel.
  call overlapPowerMinusOneHalf(iproc, nproc, bigdft_mpi%mpi_comm, iorder, &
       blocksize_dsyev, blocksize_pdgemm, orbs%norb, ovrlp_coeff)

  ! Build the new linear combinations
  if (orbs%norbp>0 )then
      call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
           ovrlp_coeff(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  end if

  call timing(iproc,'renormCoefComp','OF')
  call timing(iproc,'renormCoefComm','ON')
  if (nproc>1) then
     call mpi_allgatherv(coeff_tmp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%wfnmd%coeff(1,1), &
          tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(tmb%orbs%norb*orbs%norb,coeff_tmp(1,1),1,tmb%wfnmd%coeff(1,1),1)
  end if
  call timing(iproc,'renormCoefComm','OF')

  !call dcopy(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmb%wfnmd%coeff(1,1), 1)

  !!! Check normalization
  !!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp_tmb(1,1), tmb%orbs%norb, &
  !!     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!do iorb=1,orbs%norb
  !!    do jorb=1,orbs%norb
  !!        tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!        tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        if(iproc==0) write(200,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
  !!    end do
  !!end do

  ! Recalculate the kernel
  call calculate_density_kernel(iproc, nproc, .true., tmb%wfnmd%ld_coeff, orbs, tmb%orbs, &
       tmb%wfnmd%coeff, kernel, ovrlp_tmb)



  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
  deallocate(ovrlp_tmp,stat=istat)
  call memocc(istat,iall,'ovrlp_tmp',subname)

  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  deallocate(ovrlp_coeff,stat=istat)
  call memocc(istat,iall,'ovrlp_coeff',subname)


end subroutine reconstruct_kernel

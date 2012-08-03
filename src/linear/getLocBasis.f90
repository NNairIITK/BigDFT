!> @file
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine get_coeff(iproc,nproc,scf_mode,lzd,orbs,at,rxyz,denspot,&
    GPU, infoCoeff,ebs,nlpspd,proj,&
    SIC,tmb,fnrm,overlapmatrix,calculate_overlap_matrix,&
    tmblarge, lhphilarge, ham, ldiis_coeff)
use module_base
use module_types
use module_interfaces, exceptThisOne => get_coeff, exceptThisOneA => writeonewave
use Poisson_Solver
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc, scf_mode
type(local_zone_descriptors),intent(inout) :: lzd
type(orbitals_data),intent(in) :: orbs
type(atoms_data),intent(in) :: at
real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(GPU_pointers),intent(inout) :: GPU
integer,intent(out) :: infoCoeff
real(kind=8),intent(out) :: ebs, fnrm
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
type(SIC_data),intent(in) :: SIC
type(DFT_wavefunction),intent(inout) :: tmb
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: overlapmatrix
logical,intent(in):: calculate_overlap_matrix
type(DFT_wavefunction),intent(inout):: tmblarge
real(8),dimension(:),pointer,intent(inout):: lhphilarge
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in),optional:: ham
type(localizedDIISParameters),intent(inout),optional :: ldiis_coeff

! Local variables 
integer :: istat, iall, iorb, jorb, korb, info, inc, jjorb
real(kind=8),dimension(:),allocatable :: eval, hpsit_c, hpsit_f
real(kind=8),dimension(:,:),allocatable :: ovrlp
real(kind=8),dimension(:,:,:),allocatable :: matrixElements
real(kind=8) :: tt
type(confpot_data),dimension(:),allocatable :: confdatarrtmp
type(energy_terms) :: energs
character(len=*),parameter :: subname='get_coeff'

  ! Allocate the local arrays.  
  allocate(matrixElements(tmb%orbs%norb,tmb%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(eval(tmb%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  if(calculate_overlap_matrix) then
      if(tmb%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
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
      else if(tmb%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call getOverlapMatrix2(iproc, nproc, lzd, tmb%orbs, tmb%comon, tmb%op, tmb%psi, tmb%mad, overlapmatrix)
      else
          stop 'wrong communication_strategy_overlap'
      end if
  end if

  if(tmb%wfnmd%bs%communicate_phi_for_lsumrho) then
      call communicate_basis_for_density(iproc, nproc, lzd, tmb%orbs, tmb%psi, tmb%comsr)
  end if

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'


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


      !if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))

      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

      if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))
      allocate(tmblarge%lzd%doHamAppl(tmblarge%lzd%nlr), stat=istat)
      call memocc(istat, tmblarge%lzd%doHamAppl, 'tmblarge%lzd%doHamAppl', subname)
      tmblarge%lzd%doHamAppl=.true.
      call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
           proj,tmblarge%lzd,nlpspd,tmblarge%psi,lhphilarge,energs%eproj)
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,lhphilarge,&
           energs,SIC,GPU,.false.,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call timing(iproc,'glsynchham1','ON') !lr408t
      call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,lhphilarge,&
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
      if(tmb%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then

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
               lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)
          call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%mad, tmblarge%collcom, &
               tmblarge%psit_c, hpsit_c, tmblarge%psit_f, hpsit_f, matrixElements)
          iall=-product(shape(hpsit_c))*kind(hpsit_c)
          deallocate(hpsit_c, stat=istat)
          call memocc(istat, iall, 'hpsit_c', subname)
          iall=-product(shape(hpsit_f))*kind(hpsit_f)
          deallocate(hpsit_f, stat=istat)
          call memocc(istat, iall, 'hpsit_f', subname)

      else if(tmblarge%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call allocateCommuncationBuffersOrtho(tmblarge%comon, subname)
          call getMatrixElements2(iproc, nproc, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, tmblarge%psi, &
               lhphilarge, tmblarge%mad, matrixElements)
          call deallocateCommuncationBuffersOrtho(tmblarge%comon, subname)
      else
          stop 'wrong communication_strategy_overlap'
      end if

  else
      if(iproc==0) write(*,*) 'No Hamiltonian application required.'
      call dcopy(tmb%orbs%norb**2, ham(1,1), 1, matrixElements(1,1,1), 1)
  end if


  ! Keep the Hamiltonian since it will be overwritten by the diagonalization.
  call dcopy(tmb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
  call dcopy(tmb%orbs%norb**2, overlapmatrix(1,1),1 , ovrlp(1,1), 1)

  ! Diagonalize the Hamiltonian, either iteratively or with lapack.
  ! Make a copy of the matrix elements since dsyev overwrites the matrix and the matrix elements
  ! are still needed later.
  if(scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
      call dcopy(tmb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      do iorb=1,tmb%orbs%norb
      end do
      if(tmb%wfnmd%bpo%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, nproc, tmb%orbs, tmb%op%nsubmax, matrixElements(1,1,2), ovrlp, eval)
      else
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          call dsygv_parallel(iproc, nproc, tmb%wfnmd%bpo%blocksize_pdsyev, tmb%wfnmd%bpo%nproc_pdsyev, &
               mpi_comm_world, 1, 'v', 'l',tmb%orbs%norb, &
               matrixElements(1,1,2), tmb%orbs%norb, ovrlp, tmb%orbs%norb, eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'

      !if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
      !    if(iproc==0) write(*,*) 'copy coeffs'
          do iorb=1,orbs%norb
              call dcopy(tmb%orbs%norb, matrixElements(1,iorb,2), 1, tmb%wfnmd%coeff(1,iorb), 1)
          end do
      !else
      !    if(iproc==0) write(*,*) "don't copy coeffs"
      !end if
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

      ! keep the eigeanvalues for the preconditioning
      call vcopy(tmb%orbs%norb, eval(1), 1, tmb%orbs%eval(1), 1)
  end if

  ! TEST
  if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
      if(.not.present(ldiis_coeff)) stop 'ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION'
      !call dcopy(tmb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      call optimize_coeffs(iproc, nproc, orbs, matrixElements(1,1,1), overlapmatrix, tmb, ldiis_coeff, fnrm)
  end if

  call calculate_density_kernel(iproc, nproc, tmb%wfnmd%ld_coeff, orbs, tmb%orbs, &
       tmb%wfnmd%coeff, tmb%wfnmd%density_kernel, ovrlp)

  ! Calculate the band structure energy with matrixElements instead of wfnmd%coeff sue to the problem mentioned
  ! above (wrong size of wfnmd%coeff)
  ebs=0.d0
  do jorb=1,tmb%orbs%norb
      do korb=1,jorb
          tt = tmb%wfnmd%density_kernel(korb,jorb)*matrixElements(korb,jorb,1)
          if(korb/=jorb) tt=2.d0*tt
          ebs = ebs + tt
      end do
  end do

  ! If closed shell multiply by two.
  if(orbs%nspin==1) ebs=2.d0*ebs


  ! Project the lb coefficients on the smaller subset
  if(tmb%wfnmd%bs%use_derivative_basis) then
      inc=4
      do iorb=1,orbs%norb
          jjorb=1
          do jorb=1,tmb%orbs%norb,inc
              tt=0.d0
              do korb=1,tmb%orbs%norb
                  tt = tt + tmb%wfnmd%coeff(korb,iorb)*overlapmatrix(korb,jorb)
              end do
              tmb%wfnmd%coeff_proj(jjorb,iorb)=tt
              jjorb=jjorb+1
          end do
      end do
  else
      do iorb=1,orbs%norb
          do jorb=1,tmb%orbs%norb
              tmb%wfnmd%coeff_proj(jorb,iorb)=tmb%wfnmd%coeff(jorb,iorb)
          end do
      end do
  end if


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
    SIC, tmb, tmblarge2, lhphilarge2, energs_base, ham)
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
type(DFT_wavefunction),target,intent(inout) :: tmblarge2
real(kind=8),dimension(:),pointer,intent(inout) :: lhphilarge2
type(energy_terms),intent(inout) :: energs_base
real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ham

! Local variables
real(kind=8) :: trHold, fnrmMax, meanAlpha, ediff, noise
integer :: iorb, consecutive_rejections,istat,istart,ierr,it,iall,ilr,jorb,nsatur
real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp
real(kind=8),dimension(:,:),allocatable :: ovrlp
logical :: emergency_exit, overlap_calculated
character(len=*),parameter :: subname='getLocalizedBasis'
real(kind=8),dimension(:),pointer :: lhphi, lhphiold, lphiold, hpsit_c, hpsit_f
type(energy_terms) :: energs
character(len=3):: num
integer :: i,j , k, ncount, ist, iiorb, sdim, ldim
real(8),dimension(:),allocatable:: phiplot
real(8),dimension(2):: reducearr
real(8),save:: trH_old



  ! Allocate all local arrays.
  call allocateLocalArrays()


  call timing(iproc,'getlocbasinit','ON') !lr408t
  tmb%can_use_transposed=.false.
  if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS


  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
  consecutive_rejections=0
  trHold=1.d100

  nsatur=0
 
  call timing(iproc,'getlocbasinit','OF') !lr408t

  if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) trH_old=0.d0
  overlap_calculated=.false.
  iterLoop: do it=1,tmb%wfnmd%bs%nit_basis_optimization


      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif


      ! Orthonormalize the orbitals. If the localization regions are smaller that the global box (which
      ! will be the usual case), the orthogonalization can not be done exactly, but only approximately.
      if(iproc==0) then
          write(*,'(1x,a)',advance='no') 'Orthonormalization...'
      end if

      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if (tmblarge2%orbs%npsidim_orbs > 0) call to_zero(tmblarge2%orbs%npsidim_orbs,lhphilarge2(1))
      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, &
           tmb%psi, tmblarge2%psi)
      if(it==1) then
          call local_potential_dimensions(tmblarge2%lzd,tmblarge2%orbs,denspot%dpbox%ngatherarr(0,1))
          call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
               tmblarge2%comgp%nrecvbuf, tmblarge2%comgp%recvbuf, tmblarge2%comgp)
          call test_p2p_communication(iproc, nproc, tmblarge2%comgp)
      end if

      allocate(tmblarge2%lzd%doHamAppl(tmblarge2%lzd%nlr), stat=istat)
      call memocc(istat, tmblarge2%lzd%doHamAppl, 'tmblarge2%lzd%doHamAppl', subname)
      tmblarge2%lzd%doHamAppl=.true.
      call NonLocalHamiltonianApplication(iproc,at,tmblarge2%orbs,rxyz,&
           proj,tmblarge2%lzd,nlpspd,tmblarge2%psi,lhphilarge2,energs%eproj)
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge2%orbs,&
           tmblarge2%lzd,tmblarge2%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge2%psi,lhphilarge2,&
           energs,SIC,GPU,.false.,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge2%comgp)
      call timing(iproc,'glsynchham2','ON') !lr408t
      call SynchronizeHamiltonianApplication(nproc,tmblarge2%orbs,tmblarge2%lzd,GPU,lhphilarge2,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF') !lr408t

  iall=-product(shape(tmblarge2%lzd%doHamAppl))*kind(tmblarge2%lzd%doHamAppl)
  deallocate(tmblarge2%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'tmblarge2%lzd%doHamAppl', subname)

!DEBUG
if (iproc==0) then
   write(*,'(1x,a,4(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
   2*energs%ekin,2*energs%epot,2*energs%eproj,2*energs%ekin+2*energs%epot+2*energs%eproj
   !!write(*,'(1x,a,4(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
   !!2*av_h_sym_diff1,2*av_h_sym_diff2,2*av_h_sym_diff3,2*av_h_sym_diff1+2*av_h_sym_diff2+2*av_h_sym_diff3
endif
!END DEBUG


  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') ' Orthoconstraint... '
      end if

      call copy_basis_specifications(tmb%wfnmd%bs, tmblarge2%wfnmd%bs, subname)
      call copy_orthon_data(tmb%orthpar, tmblarge2%orthpar, subname)

      call transpose_localized(iproc, nproc, tmblarge2%orbs, tmblarge2%collcom, lhphilarge2, hpsit_c, hpsit_f, tmblarge2%lzd)
      ncount=sum(tmblarge2%collcom%nrecvcounts_c)
      if(ncount>0) call dcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*sum(tmblarge2%collcom%nrecvcounts_f)
      if(ncount>0) call dcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      if(tmblarge2%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
          call calculate_energy_and_gradient_linear(iproc, nproc, it, &
               tmb%wfnmd%density_kernel, &
               ldiis, consecutive_rejections, &
               fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, &
               meanAlpha, emergency_exit, &
               tmb, lhphi, lhphiold, &
               tmblarge2, lhphilarge2, overlap_calculated, ovrlp, energs_base, hpsit_c, hpsit_f)
       else
          call calculate_energy_and_gradient_linear(iproc, nproc, it, &
               tmb%wfnmd%density_kernel, &
               ldiis, consecutive_rejections, &
               fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, &
               meanAlpha, emergency_exit, &
               tmb, lhphi, lhphiold, &
               tmblarge2, lhphilarge2, overlap_calculated, ovrlp, energs_base)
       end if

           if (emergency_exit) then
               tmblarge2%can_use_transposed=.false.
               call dcopy(tmb%orbs%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
           end if 


      !!!! trH is now the total energy (name is misleading, correct this)
      !!!if(orbs%nspin==1) trH=2.d0*trH
      !!!trH=trH-energs_base%eh+energs_base%exc-energs_base%evxc-energs_base%eexctX+energs_base%eion+energs_base%edisp


      !!!plot gradient
      !!allocate(phiplot(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f))
      !!ist=1
      !!do iorb=1,tmbopt%orbs%norbp
      !!    iiorb=tmbopt%orbs%isorb+iorb
      !!    ilr=tmbopt%orbs%inwhichlocreg(iiorb)
      !!    sdim=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
      !!    ldim=tmbopt%lzd%glr%wfd%nvctr_c+7*tmbopt%lzd%glr%wfd%nvctr_f
      !!    call to_zero(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f, phiplot(1))
      !!    call Lpsi_to_global2(iproc, nproc, sdim, ldim, tmbopt%orbs%norb, tmbopt%orbs%nspinor, 1, tmbopt%lzd%glr, &
      !!         tmbopt%lzd%llr(ilr), lhphiopt(ist), phiplot(1))
      !!    !!do istat=1,sdim
      !!    !!    write(300,*) lhphiopt(ist+istat-1)
      !!    !!end do
      !!    !!do istat=1,ldim
      !!    !!    write(400,*) phiplot(istat)
      !!    !!end do
      !!    !!call small_to_large_locreg(iproc, nproc, tmbopt%lzd, tmblarge2%lzd, tmbopt%orbs, tmblarge2%orbs, &
      !!    !!     tmbopt%psi, tmblarge2%psi)
      !!    write(num,'(i3.3)') iiorb
      !!    call plot_wf('gradient'//num,2,at,1.d0,tmbopt%lzd%glr,tmb%lzd%hgrids(1),tmb%lzd%hgrids(2),tmb%lzd%hgrids(3),rxyz,phiplot)
      !!    ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
      !!    ist = ist + ncount
      !!end do
      !!deallocate(phiplot)



  

      ediff=trH-trH_old
      !noise=tmb%wfnmd%bs%gnrm_mult*fnrm*tmb%orbs%norb
      noise=tmb%wfnmd%bs%gnrm_mult*fnrm*abs(trH)
      if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY )then
          if (ediff<0.d0 .and. abs(ediff) < noise) then
              if(iproc==0) write(*,'(a)') 'target function seems to saturate, increase nsatur...'
              nsatur=nsatur+1
          else if (abs(ediff) < noise) then
              if(iproc==0) write(*,'(a)') 'target function increases, but smaller than noise. Consider convergence.'
              nsatur=tmb%wfnmd%bs%nsatur_inner
          else
              nsatur=0
          end if
      end if


      !if(iproc==0) write(*,'(a,2es14.6)') 'fnrm*gnrm_in/gnrm_out, energy diff',fnrm*gnrm_in/gnrm_out, trH-trH_old
      ! Write some informations to the screen.
      if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
          write(*,'(1x,a,i6,2es15.7,f17.10,2es13.4)') 'iter, fnrm, fnrmMax, trace, diff, noise level', &
          it, fnrm, fnrmMax, trH, ediff, noise
      if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
          write(*,'(1x,a,i6,2es15.7,f17.10,2es13.4)') 'iter, fnrm, fnrmMax, ebs, diff, noise level', &
          it, fnrm, fnrmMax, trH, ediff,noise
      !!if((fnrm*gnrm_in/gnrm_out < tmb%wfnmd%bs%conv_crit*tmb%wfnmd%bs%conv_crit_ratio) .or. &
      !!    it>=tmb%wfnmd%bs%nit_basis_optimization .or. emergency_exit) then
      !!    if(fnrm*gnrm_in/gnrm_out < tmb%wfnmd%bs%conv_crit*tmb%wfnmd%bs%conv_crit_ratio) then
      if(it>=tmb%wfnmd%bs%nit_basis_optimization .or. emergency_exit .or. nsatur>=tmb%wfnmd%bs%nsatur_inner) then
          !!if(fnrm*gnrm_in/gnrm_out < tmb%wfnmd%bs%conv_crit*tmb%wfnmd%bs%conv_crit_ratio) then
          if(nsatur>=tmb%wfnmd%bs%nsatur_inner .and. .not.emergency_exit) then
              if(iproc==0) then
                  write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                      write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
                  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                      write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=it
          else if(it>=tmb%wfnmd%bs%nit_basis_optimization .and. .not.emergency_exit) then
          !!if(it>=tmb%wfnmd%bs%nit_basis_optimization) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=0
          else if(emergency_exit) then
              if(iproc==0) then
                  write(*,'(1x,a,i0,a)') 'WARNING: emergency exit after ',it, &
                      ' iterations to keep presumably good TMBs before they deteriorate too much.'
                  write (*,'(1x,a,2es15.7,f12.7)') '>>WRONG OUTPUT<< Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=-1
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          call calculate_overlap_transposed(iproc, nproc, tmblarge2%orbs, tmblarge2%mad, tmblarge2%collcom, &
               tmblarge2%psit_c, hpsit_c_tmp, tmblarge2%psit_f, hpsit_f_tmp, ham)

          exit iterLoop
      end if
      trH_old=trH


      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lhphi, lphiold, alpha, trH, meanAlpha, alphaDIIS)
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmblarge2%can_use_transposed) then
          iall=-product(shape(tmblarge2%psit_c))*kind(tmblarge2%psit_c)
          deallocate(tmblarge2%psit_c, stat=istat)
          call memocc(istat, iall, 'tmblarge2%psit_c', subname)
          iall=-product(shape(tmblarge2%psit_f))*kind(tmblarge2%psit_f)
          deallocate(tmblarge2%psit_f, stat=istat)
          call memocc(istat, iall, 'tmblarge2%psit_f', subname)
          tmblarge2%can_use_transposed=.false.
      end if


      if(iproc==0) WRITE(*,*) 'WARNING: NO RECONSTRUCTION OF KERNEL'
      !!if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      !!    call reconstruct_kernel(iproc, nproc, orbs, tmb, ovrlp, overlap_calculated, tmb%wfnmd%density_kernel)
      !!end if


  end do iterLoop

  ! Keep the values for the next iteration
  reducearr(1)=0.d0
  reducearr(2)=0.d0
  do iorb=1,tmb%orbs%norbp
      reducearr(1)=reducearr(1)+alpha(iorb)
      reducearr(2)=reducearr(2)+alphaDIIS(iorb)
  end do
  call mpiallred(reducearr(1), 2, mpi_sum, mpi_comm_world, ierr)
  reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
  reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)

  ldiis%alphaSD=reducearr(1)
  ldiis%alphaDIIS=reducearr(2)




  ! Deallocate potential
  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)


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

      allocate(hpsit_c(sum(tmblarge2%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)

      allocate(hpsit_f(7*sum(tmblarge2%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)

      allocate(hpsit_c_tmp(sum(tmblarge2%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c_tmp, 'hpsit_c_tmp', subname)

      allocate(hpsit_f_tmp(7*sum(tmblarge2%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f_tmp, 'hpsit_f_tmp', subname)


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
integer :: istart, iorb, iiorb, ilr, ncount, owa, owanext
  
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
          owa=tmb%orbs%onwhichatom(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          if(iiorb<tmb%orbs%norb) then
              owanext=tmb%orbs%onwhichatom(iiorb+1)
          else
              owanext=tmb%lzd%nlr+1
          end if
          !!if(owa==owanext) then
              call daxpy(ncount, -alpha(iorb), lhphi(istart), 1, tmb%psi(istart), 1)
          !!end if
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







subroutine diagonalizeHamiltonian2(iproc, nproc, orbs, nsubmax, HamSmall, ovrlp, eval)
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
integer :: iproc, nproc, nsubmax
type(orbitals_data), intent(inout) :: orbs
real(kind=8),dimension(orbs%norb, orbs%norb),intent(inout) :: HamSmall
real(kind=8),dimension(orbs%norb, orbs%norb),intent(in) :: ovrlp
real(kind=8),dimension(orbs%norb),intent(out) :: eval

! Local variables
integer :: lwork, info, istat, iall, iorb, jorb
real(kind=8),dimension(:),allocatable :: work
character(len=*),parameter :: subname='diagonalizeHamiltonian'

  ! temp change
  real(8),dimension(:),allocatable:: eval1,beta
  real(8),dimension(:,:), allocatable :: vr,vl,ovrlp_copy
  !real(8),dimension(:,:), allocatable :: inv_ovrlp,ks
  integer :: ierr
  real(8) :: temp, tt, ddot

  call timing(iproc,'diagonal_seq  ','ON')

  !! OLD VERSION #####################################################################################################
  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=int(work(1))

  !!! find inverse overlap and premultiply Hamiltonian
!!  allocate(inv_ovrlp(1:orbs%norb,1:orbs%norb))
  !!call dcopy(orbs%norb**2,ovrlp(1,1),1,inv_ovrlp(1,1),1)
  !!! Exact inversion
  !!call dpotrf('l', orbs%norb, inv_ovrlp(1,1), orbs%norb, info)
  !!if(info/=0) then
  !!   write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
  !!   stop
  !!end if
  !!call dpotri('l', orbs%norb, inv_ovrlp(1,1), orbs%norb, info)
  !!if(info/=0) then
  !!   write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
  !!   stop
  !!end if

  !!! fill the upper triangle
  !!do iorb=1,orbs%norb
  !!   do jorb=1,iorb-1
  !!      inv_ovrlp(jorb,iorb)=inv_ovrlp(iorb,jorb)
  !!   end do
  !!end do

  !!allocate(ks(1:orbs%norb,1:orbs%norb))
  !!call dgemm('n','n', orbs%norb,orbs%norb,orbs%norb,1.d0,inv_ovrlp(1,1),orbs%norb,&
  !!     HamSmall(1,1),orbs%norb,0.d0,ks(1,1),orbs%norb)
  !!call dcopy(orbs%norb**2,ks(1,1),1,HamSmall(1,1),1)
  !!deallocate(ks)
  !!deallocate(inv_ovrlp)
  !!!!!!!!!!!
  !!allocate(ham_copy(1:orbs%norb,1:orbs%norb))


  !allocate(ovrlp_copy(1:orbs%norb,1:orbs%norb), stat=istat)
  !call memocc(istat, ovrlp_copy, 'ovrlp_copy', subname)
  !allocate(vl(1:orbs%norb,1:orbs%norb), stat=istat)
  !call memocc(istat, vl, 'vl', subname)
  !allocate(vr(1:orbs%norb,1:orbs%norb), stat=istat)
  !call memocc(istat, vr, 'vr', subname)
  !allocate(eval1(1:orbs%norb), stat=istat)
  !call memocc(istat, eval1, 'eval1', subname)
  !allocate(beta(1:orbs%norb), stat=istat)
  !call memocc(istat, beta, 'beta', subname)

  !call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp_copy(1,1), 1)

  !!$call dggev('v', 'v',orbs%norb,&
  !!$      HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
  !!$      vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
  !!call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
  !!     orbs%norb, WORK, LWORK, ierr )
  !!$lwork=work(1) 

  ! Deallocate the work array and reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
!!  call dcopy(orbs%norb**2, HamSmall(1,1), 1, vl(1,1), 1)
!!  call dcopy(orbs%norb**2, ovrlp(1,1), 1, vr(1,1), 1)
!!  call dcopy(orbs%norb**2, HamSmall(1,1), 1, inv_ovrlp(1,1), 1)
!!$  call dcopy(orbs%norb**2, ovrlp(1,1), 1, ks(1,1), 1)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
!!  do iorb=1,orbs%norb
!!    do jorb=1,orbs%norb
!!      write(200+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!    end do
!!    write(250+iproc,*) iorb,eval(iorb)
!!  end do
!!  call dcopy(orbs%norb**2, vl(1,1), 1, HamSmall(1,1), 1)
!!  call dcopy(orbs%norb**2, vr(1,1), 1, ovrlp(1,1), 1)


  !do iorb=1,orbs%norb
  !  eval(iorb) = eval(iorb) / beta(iorb)
  !end do

!!$$$  lwork=-1
!!$$$  call dggev('v', 'v',orbs%norb,&
!!$$$        HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
!!$$$        vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
!!$$$  lwork=work(1) 
!!$$$  iall=-product(shape(work))*kind(work)
!!$$$  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!$$$  call memocc(istat, iall, 'work', subname)
!!$$$  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!$$$  call memocc(istat, work, 'work', subname)
!!$$$  call dggev('v', 'v',orbs%norb,&
!!$$$        HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
!!$$$        vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
!!$$$
!!$$$        hamsmall=vl
!!$$$  do iorb=1,orbs%norb
!!$$$     do jorb=iorb,orbs%norb
!!$$$        if (eval(jorb)/beta(jorb) < eval(iorb)/beta(iorb)) then
!!$$$           temp = eval(iorb)
!!$$$           temp_vec = HamSmall(:,iorb)
!!$$$           eval(iorb) = eval(jorb)
!!$$$           eval(jorb) = temp
!!$$$           HamSmall(:,iorb) = HamSmall(:,jorb)
!!$$$           HamSmall(:,jorb) = temp_vec
!!$$$           temp=beta(iorb)
!!$$$           beta(iorb)=beta(jorb)
!!$$$           beta(jorb)=temp
!!$$$        end if
!!$$$     end do
!!$$$  end do
!!$$$
!!$$$
!!$$$
!!$$$  call dcopy(orbs%norb**2, ks(1,1), 1, ovrlp(1,1), 1)
!!$$$  do iorb=1,orbs%norb
!!$$$      call dgemv('n', orbs%norb, orbs%norb, 1.d0, ovrlp(1,1), &
!!$$$           orbs%norb, hamsmall(1,iorb), 1, 0.d0, vl(1,iorb), 1)
!!$$$      tt=ddot(orbs%norb, hamsmall(1,iorb),  1, vl(1,iorb), 1)
!!$$$      call dscal(orbs%norb, 1/sqrt(tt), hamsmall(1,iorb), 1)
!!$$$  end do
!!$$$
!!$$$
!!$$$
!!$$$
!!$$$!!  do iorb=1,orbs%norb
!!$$$!!    do jorb=1,orbs%norb
!!$$$!!      write(300+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!$$$!!    end do
!!$$$!!    write(350+iproc,*) iorb,eval(iorb)/beta(iorb)
!!$$$!!  end do
!!$$$!!
!!$$$!!  lwork=-1
!!$$$!!  call dcopy(orbs%norb**2, inv_ovrlp(1,1), 1, HamSmall(1,1), 1)
!!$$$!!  call dcopy(orbs%norb**2, ks(1,1), 1, ovrlp(1,1), 1)
!!$$$!!  call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
!!$$$!!       orbs%norb, WORK, LWORK, ierr )
!!$$$!!  ! Deallocate the work array and reallocate it with the optimal size
!!$$$!!  lwork=work(1) 
!!$$$!!  iall=-product(shape(work))*kind(work)
!!$$$!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!$$$!!  call memocc(istat, iall, 'work', subname)
!!$$$!!  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!$$$!!  call memocc(istat, work, 'work', subname)
!!$$$!!
!!$$$!!  call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
!!$$$!!       orbs%norb, WORK, LWORK, ierr )
!!$$$!!
!!$$$!!  HamSmall=vl
!!$$$!!  do iorb=1,orbs%norb
!!$$$!!     do jorb=iorb,orbs%norb
!!$$$!!        if (eval(jorb) < eval(iorb)) then
!!$$$!!           temp = eval(iorb)
!!$$$!!           temp_vec = HamSmall(:,iorb)
!!$$$!!           eval(iorb) = eval(jorb)
!!$$$!!           eval(jorb) = temp
!!$$$!!           HamSmall(:,iorb) = HamSmall(:,jorb)
!!$$$!!           HamSmall(:,jorb) = temp_vec
!!$$$!!        end if
!!$$$!!     end do
!!$$$!!  end do
!!$$$!!
!!$$$!!  do iorb=1,orbs%norb
!!$$$!!    do jorb=1,orbs%norb
!!$$$!!      write(400+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!$$$!!    end do
!!$$$!!    write(450+iproc,*) iorb,eval(iorb)
!!$$$!!  end do
!!$$$!!
!!$$$!!!  do iorb=1,orbs%norb
!!$$$!!!    write(36,*) vl(:,iorb)
!!$$$!!!    write(37,*) vr(:,iorb)
!!$$$!!!  end do
!!$$$!!!  write(36,*) ''
!!$$$!!!  write(37,*) ''
!!$$$!!!  write(38,*) 'eval',eval
!!$$$!!!  write(38,*) 'eval1',eval1
!!$$$!!!  !write(38,*) 'beta',beta
!!$$$
!!$$$  ! Deallocate the work array.
!!$$$  iall=-product(shape(work))*kind(work)
!!$$$  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!$$$  call memocc(istat, iall, 'work', subname)
!!$$$  
!!$$$  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!$$$  ! the first entry of each vector is positive.
!!$$$  do iorb=1,orbs%norb
!!$$$      if(HamSmall(1,iorb)<0.d0) then
!!$$$          do jorb=1,orbs%norb
!!$$$              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!$$$          end do
!!$$$      end if
!!$$$  end do
!!$$$  !! #################################################################################################################
!!$$$
!!$$$  deallocate(vl)
!!$$$  deallocate(vr)
!!$$$  deallocate(eval1)
!!$$$  deallocate(beta)
!!$$$
!!$$$if (.false.) then
!!$$$  allocate(vl(1:orbs%norb,1:orbs%norb))
!!$$$  allocate(vr(1:orbs%norb,1:orbs%norb))
!!$$$  allocate(eval1(1:orbs%norb))
!!$$$  allocate(eval2(1:orbs%norb))
!!$$$
!!$$$
!!$$$  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!$$$  ! check eigenvalues of overlap matrix
!!$$$      call DGEEV( 'v','v', orbs%norb, ovrlp(1,1), orbs%norb, eval2, eval1, VL, orbs%norb, VR,&
!!$$$                  orbs%norb, WORK, LWORK, ierr )
!!$$$!  write(40,*) 'eval',eval2
!!$$$!  write(40,*) 'eval1',eval1
!!$$$!  write(40,*) 'sum',sum(eval2)
!!$$$
!!$$$!  do iorb=1,orbs%norb
!!$$$!    write(44,*) vl(:,iorb)
!!$$$!    write(45,*) vr(:,iorb)
!!$$$!  end do
!!$$$!  write(44,*) ''
!!$$$!  write(45,*) ''
!!$$$
!!$$$  write(41,*) 'sum olap eigs',sum(eval2)
!!$$$
!!$$$  deallocate(work)
!!$$$  deallocate(vl)
!!$$$  deallocate(vr)
!!$$$  deallocate(eval1)
!!$$$  deallocate(eval2)
!!$$$end if

  !!!! NEW VERSION #####################################################################################################
  !!! Determine the maximal number of non-zero subdiagonals
  !!!!nsubmax=0
  !!!!do iorb=1,orbs%norb
  !!!!    nsub=0
  !!!!    do jorb=orbs%norb,iorb+1,-1
  !!!!        if(Hamsmall(jorb,iorb)/=0.d0) then
  !!!!            nsub=jorb-iorb
  !!!!            exit
  !!!!        end if
  !!!!    end do
  !!!!    if(iproc==0) write(*,*) 'iorb,nsub',iorb,nsub
  !!!!    nsubmax=max(nsub,nsubmax)
  !!!!end do
  !!!!if(iproc==0) write(*,*) 'nsubmax',nsubmax
  !!!!if(iproc==0) then
  !!!!      do iorb=1,orbs%norb
  !!!!           write(*,'(14es10.3)') (hamsmall(iorb,jorb), jorb=1,orbs%norb)
  !!!!      end do
  !!!!end if

  !!! Copy to banded format
  !!allocate(ham_band(nsubmax+1,orbs%norb), stat=istat)
  !!call memocc(istat, ham_band, 'ham_band', subname)
  !!allocate(ovrlp_band(nsubmax+1,orbs%norb), stat=istat)
  !!call memocc(istat, ovrlp_band, 'ovrlp_band', subname)
  !!do iorb=1,orbs%norb
  !!    do jorb=iorb,min(iorb+nsubmax,orbs%norb)
  !!        ham_band(1+jorb-iorb,iorb)=HamSmall(jorb,iorb)
  !!        ovrlp_band(1+jorb-iorb,iorb)=ovrlp(jorb,iorb)
  !!    end do
  !!end do
  !!!!if(iproc==0) then
  !!!!      write(*,*) '+++++++++++++++++++++++++++++'
  !!!!      do iorb=1,nsubmax+1
  !!!!           write(*,'(14es10.3)') (ham_band(iorb,jorb), jorb=1,orbs%norb)
  !!!!      end do
  !!!!end if


  !!!!! Get the optimal work array size
  !!!!lwork=-1 
  !!!!allocate(work(1), stat=istat)
  !!!!call memocc(istat, work, 'work', subname)
  !!!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  !!!!lwork=work(1) 

  !!!!! Deallocate the work array ane reallocate it with the optimal size
  !!!!iall=-product(shape(work))*kind(work)
  !!!!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  !!!!call memocc(istat, iall, 'work', subname)
  !!allocate(work(3*orbs%norb), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  !!call memocc(istat, work, 'work', subname)

  !!! Diagonalize the Hamiltonian
  !!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  !!call dsbgv('v', 'l', orbs%norb, nsubmax, nsubmax, ham_band(1,1), nsubmax+1, ovrlp_band(1,1), nsubmax+1, &
  !!     eval(1), HamSmall(1,1), orbs%norb, work, info)

  !!! Deallocate the work array.
  !!iall=-product(shape(work))*kind(work)
  !!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  !!call memocc(istat, iall, 'work', subname)
  !!
  !!! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
  !!! the first entry of each vector is positive.
  !!do iorb=1,orbs%norb
  !!    if(HamSmall(1,iorb)<0.d0) then
  !!        do jorb=1,orbs%norb
  !!            HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
  !!        end do
  !!    end if
  !!end do


  !!iall=-product(shape(ham_band))*kind(ham_band)
  !!deallocate(ham_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ham_band' 
  !!call memocc(istat, iall, 'ham_band', subname)

  !!iall=-product(shape(ovrlp_band))*kind(ovrlp_band)
  !!deallocate(ovrlp_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ovrlp_band' 
  !!call memocc(istat, iall, 'ovrlp_band', subname)

  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2









!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in) :: keyv_c
  integer,dimension(mseg_f),intent(in) :: keyv_f
  integer,dimension(2,mseg_c),intent(in) :: keyg_c
  integer,dimension(2,mseg_f),intent(in) :: keyg_f
  real(wp),dimension(0:3),intent(in) :: scal
  real(wp),dimension(mvctr_c),intent(in) :: psi_c
  real(wp),dimension(7,mvctr_f),intent(in) :: psi_f
  type(workarrays_quartic_convolutions),intent(out) :: work
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !!!$omp parallel default(private) &
  !!!$omp shared(scal,psig_c,psig_f,x_f1,x_f2,x_f3) &
  !!!$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f)
  !!! coarse part
  !!!$omp do
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
        work%xy_c(i2,i,i3)=psi_c(i-i0+jj)*scal(0)
        work%xz_c(i3,i,i2)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !!!$omp enddo
  !!! fine part
  !!!$omp do
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_f1(i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xx_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xy_f(1,i2,i,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xz_f(1,i3,i,i2)=psi_f(1,i-i0+jj)*scal(1)

        work%xy_f2(i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xx_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xy_f(2,i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xz_f(2,i3,i,i2)=psi_f(2,i-i0+jj)*scal(1)

        work%xx_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xy_f(3,i2,i,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xz_f(3,i3,i,i2)=psi_f(3,i-i0+jj)*scal(2)

        work%xz_f4(i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)
        work%xx_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xy_f(4,i2,i,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xz_f(4,i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)

        work%xx_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xy_f(5,i2,i,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xz_f(5,i3,i,i2)=psi_f(5,i-i0+jj)*scal(2)

        work%xx_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xy_f(6,i2,i,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xz_f(6,i3,i,i2)=psi_f(6,i-i0+jj)*scal(2)

        work%xx_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xy_f(7,i2,i,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xz_f(7,i3,i,i2)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !!!$omp enddo
 !!!$omp end parallel

END SUBROUTINE uncompress_for_quartic_convolutions








subroutine build_new_linear_combinations(iproc, nproc, lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
use module_base
use module_types
implicit none

!Calling arguments
integer,intent(in) :: iproc, nproc
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: orbs
type(overlapParameters),intent(in) :: op
integer,intent(in) :: nrecvbuf
real(kind=8),dimension(nrecvbuf),intent(in) :: recvbuf
real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: omat
logical,intent(in) :: reset
real(kind=8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: lphi

! Local variables
integer :: jst, ilrold, iorb, iiorb, ilr, ncount, jorb, jjorb, ldim, indout, gdim, iorbref
integer :: istart, iend, iseg, start, kseg, kold, kstart, kend 
real(kind=8) :: tt

   call timing(iproc,'build_lincomb ','ON')

      ! Build new lphi
      if(reset) then
          !!lphi=0.d0
          call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1))
      end if

      indout=1
      ilrold=-1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          !if(ilr>ilrold) then
          if(ilr/=ilrold) then
              iorbref=iorb
          end if
          gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)
              !jjorb=op%overlaps(jorb,ilr)
              jst=op%indexInRecvBuf(iorbref,jjorb)
              ldim=op%wfd_overlap(jorb,iorbref)%nvctr_c+7*op%wfd_overlap(jorb,iorbref)%nvctr_f
              tt=omat(jjorb,iiorb)
              !!tt=tt*lzd%cutoffweight(jjorb,iiorb)
              !Coarse part
              kold=1
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_c
                  istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg)
                  iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg)
                  ncount=iend-istart+1
                  inner_loop: do kseg=kold,lzd%llr(ilr)%wfd%nseg_c
                     kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg)
                     kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg)
                     if(kstart <= iend .and. kend >= istart) then 
                        kold = kseg + 1
                        start = lzd%llr(ilr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
                        call daxpy(ncount, tt, recvBuf(jst), 1, lphi(indout+start-1), 1)
                        jst=jst+ncount
                        exit inner_loop
                     end if
                  end do inner_loop
              end do
              ! Fine part
              kold = 1
              jst=op%indexInRecvBuf(iorbref,jjorb)              
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_f
                 istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 start = op%wfd_overlap(jorb,iorbref)%keyvglob(iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 ncount=7*(iend-istart+1)
                 inner_loop2: do kseg=kold,lzd%llr(ilr)%wfd%nseg_f
                    kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    if(kstart <= iend .and. kend >= istart) then 
                       kold = kseg + 1
                       start = lzd%llr(ilr)%wfd%nvctr_c+(lzd%llr(ilr)%wfd%keyvglob(kseg+lzd%llr(ilr)%wfd%nseg_c) +&
                              max(0,istart-kstart)-1)*7
                       call daxpy(ncount, tt, recvBuf(jst+op%wfd_overlap(jorb,iorbref)%nvctr_c), 1,&
                               lphi(indout+start), 1)
                       jst=jst+ncount
                       exit inner_loop2
                    end if
                  end do inner_loop2
              end do
          end do
          indout=indout+gdim
          ilrold=ilr

      end do

   call timing(iproc,'build_lincomb ','OF')
          

end subroutine build_new_linear_combinations









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
  call to_zero(orbslarge%npsidim_orbs, philarge(1))
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
  call to_zero(orbssmall%npsidim_orbs, phismall(1))
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
  real(kind=8),dimension(tmbopt%orbs%norbp),intent(out) :: alpha, alphaDIIS
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





subroutine plot_gradient(iproc, nproc, num, lzd, orbs, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, num
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi

  ! Local variables
  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt, gnrm
  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
  character(len=*),parameter:: subname='flatten_at_boundaries'
  

  istc=1
  istf=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)

      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
      call memocc(istat, psig, 'psig', subname)
      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))

      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           psi(istc), psi(istf), psig)

      ii2=nint(lzd%llr(ilr)%d%n2/2.d0)
      ii3=nint(lzd%llr(ilr)%d%n3/2.d0)
      do i1=0,lzd%llr(ilr)%d%n1
          write(num+iiorb,*) i1, psig(i1,1,ii2,1,ii3,1)
      end do

      !!call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
      !!     0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
      !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
      !!     lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
      !!     lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
      !!     lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
      !!     psig, psi(istc), psi(istf))

      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

      iall=-product(shape(psig))*kind(psig)
      deallocate(psig,stat=istat)
      call memocc(istat,iall,'psig',subname)


  end do

end subroutine plot_gradient





subroutine reconstruct_kernel(iproc, nproc, orbs, tmb, ovrlp_tmb, overlap_calculated, kernel)
  use module_base
  use module_types
  use module_interfaces, except_this_one => reconstruct_kernel
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ovrlp_tmb
  logical,intent(out):: overlap_calculated
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: kernel

  ! Local variables
  integer:: istat, iorb, jorb, ierr, iall
  real(8):: ddot, tt, tt2, tt3
  real(8),dimension(:,:),allocatable:: coeff_tmp, ovrlp_tmp, ovrlp_coeff
  character(len=*),parameter:: subname='reconstruct_kernel'


  allocate(coeff_tmp(tmb%orbs%norb,max(orbs%norbp,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
  allocate(ovrlp_tmp(orbs%norb,max(orbs%norbp,1)), stat=istat)
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  ! Calculate the overlap matrix between the TMBs.
  if(tmb%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
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
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, tmb%collcom, &
           tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, ovrlp_tmb)
  else if (tmb%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
      call getOverlapMatrix2(iproc, nproc, tmb%lzd, tmb%orbs, tmb%comon, tmb%op, tmb%psi, tmb%mad, ovrlp_tmb)
  end if
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


  ! Gather together the complete matrix
  call mpi_allgatherv(ovrlp_tmp(1,1), orbs%norb*orbs%norbp, mpi_double_precision, ovrlp_coeff(1,1), &
       orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, mpi_comm_world, ierr)

  ! WARNING: this is the wrong mad, but it does not matter for iorder=0
  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, tmb%mad, ovrlp_coeff)

  ! Build the new linear combinations
  if (orbs%norbp>0 )then
      call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
           ovrlp_coeff(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  end if

  call mpi_allgatherv(coeff_tmp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%wfnmd%coeff(1,1), &
       tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, mpi_comm_world, ierr)

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
  call calculate_density_kernel(iproc, nproc, tmb%wfnmd%ld_coeff, orbs, tmb%orbs, &
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

               
              

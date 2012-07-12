subroutine get_coeff(iproc,nproc,scf_mode,lzd,orbs,at,rxyz,denspot,&
    GPU, infoCoeff,ebs,nlpspd,proj,&
    SIC,tmbmix,tmb,fnrm,overlapmatrix,calculate_overlap_matrix,&
    tmblarge, lhphilarge, ldiis_coeff)
use module_base
use module_types
use module_interfaces, exceptThisOne => get_coeff, exceptThisOneA => writeonewave
use Poisson_Solver
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, scf_mode
type(local_zone_descriptors),intent(inout):: lzd
type(orbitals_data),intent(in) :: orbs
type(atoms_data),intent(in):: at
real(8),dimension(3,at%nat),intent(in):: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(GPU_pointers),intent(inout):: GPU
integer,intent(out):: infoCoeff
real(8),intent(out):: ebs, fnrm
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(SIC_data),intent(in):: SIC
type(DFT_wavefunction),intent(inout):: tmbmix, tmb
real(8),dimension(tmbmix%orbs%norb,tmbmix%orbs%norb),intent(inout):: overlapmatrix
logical,intent(in):: calculate_overlap_matrix
type(DFT_wavefunction),intent(inout):: tmblarge
real(8),dimension(:),pointer,intent(inout):: lhphilarge
type(localizedDIISParameters),intent(inout),optional:: ldiis_coeff

! Local variables 
integer:: istat, iall, iorb, jorb, korb, info, inc, jjorb!,borb
real(8),dimension(:),allocatable:: eval, lhphi, psit_c, psit_f, hpsit_c, hpsit_f
real(8),dimension(:,:),allocatable:: ovrlp!,test!,ovrlp_cpy
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8):: tt
type(confpot_data),dimension(:),allocatable :: confdatarrtmp
type(energy_terms) :: energs!, energs2
character(len=*),parameter:: subname='get_coeff'
!For debug
integer :: ldim,istart,lwork,iiorb,ilr,ind2,ncnt
character(len=1) :: num
real(8),dimension(:),allocatable :: locrad_tmp
!!type(DFT_wavefunction):: tmblarge
real(8),dimension(:,:),allocatable:: locregCenter

  ! Allocate the local arrays.  
  allocate(matrixElements(tmbmix%orbs%norb,tmbmix%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(eval(tmbmix%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlp(tmbmix%orbs%norb,tmbmix%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  if(calculate_overlap_matrix) then
      if(tmbmix%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
          if(.not.tmbmix%can_use_transposed) then
              if(.not.associated(tmbmix%psit_c)) then
                  allocate(tmbmix%psit_c(sum(tmbmix%collcom%nrecvcounts_c)), stat=istat)
                  call memocc(istat, tmbmix%psit_c, 'tmbmix%psit_c', subname)
              end if
              if(.not.associated(tmbmix%psit_f)) then
                  allocate(tmbmix%psit_f(7*sum(tmbmix%collcom%nrecvcounts_f)), stat=istat)
                  call memocc(istat, tmbmix%psit_f, 'tmbmix%psit_f', subname)
              end if
              call transpose_localized(iproc, nproc, tmbmix%orbs, tmbmix%collcom, tmbmix%psi, tmbmix%psit_c, tmbmix%psit_f, lzd)
              tmbmix%can_use_transposed=.true.
          end if
          call calculate_overlap_transposed(iproc, nproc, tmbmix%orbs, tmbmix%mad, tmbmix%collcom, tmbmix%psit_c, &
               tmbmix%psit_c, tmbmix%psit_f, tmbmix%psit_f, overlapmatrix)
      else if(tmbmix%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
          call getOverlapMatrix2(iproc, nproc, lzd, tmbmix%orbs, tmbmix%comon, tmbmix%op, tmbmix%psi, tmbmix%mad, overlapmatrix)
      else
          stop 'wrong communication_strategy_overlap'
      end if
  end if

  if(tmbmix%wfnmd%bs%communicate_phi_for_lsumrho) then
      call communicate_basis_for_density(iproc, nproc, lzd, tmbmix%orbs, tmbmix%psi, tmbmix%comsr)
  end if

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'


  allocate(lhphi(max(tmbmix%orbs%npsidim_orbs,tmbmix%orbs%npsidim_comp)), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)
  allocate(lzd%doHamAppl(lzd%nlr), stat=istat)
  call memocc(istat, lzd%doHamAppl, 'lzd%doHamAppl', subname)
  lzd%doHamAppl=.true.
  allocate(confdatarrtmp(tmbmix%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmbmix%orbs%norbp)


  if (tmbmix%orbs%npsidim_orbs > 0) call to_zero(tmbmix%orbs%npsidim_orbs,lhphi(1))
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))

  if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))
  if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
  call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmbmix%orbs, tmblarge%orbs, tmbmix%psi, tmblarge%psi)

  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)
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

  if(iproc==0) write(*,'(1x,a)') 'done.'



  ! Calculate the matrix elements <phi|H|phi>.
  if(tmbmix%wfnmd%bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then

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


  ! Symmetrize the Hamiltonian
  call dcopy(tmbmix%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
  call dcopy(tmbmix%orbs%norb**2, overlapmatrix(1,1),1 , ovrlp(1,1), 1)

  ! Diagonalize the Hamiltonian, either iteratively or with lapack.
  ! Make a copy of the matrix elements since dsyev overwrites the matrix and the matrix elements
  ! are still needed later.
  if(scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
      call dcopy(tmbmix%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      do iorb=1,tmbmix%orbs%norb
      end do
      if(tmbmix%wfnmd%bpo%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, nproc, tmbmix%orbs, tmbmix%op%nsubmax, matrixElements(1,1,2), ovrlp, eval)
      else
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          call dsygv_parallel(iproc, nproc, tmbmix%wfnmd%bpo%blocksize_pdsyev, tmbmix%wfnmd%bpo%nproc_pdsyev, &
               mpi_comm_world, 1, 'v', 'l',tmbmix%orbs%norb, &
               matrixElements(1,1,2), tmbmix%orbs%norb, ovrlp, tmbmix%orbs%norb, eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'

      !if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
      !    if(iproc==0) write(*,*) 'copy coeffs'
          do iorb=1,orbs%norb
              call dcopy(tmbmix%orbs%norb, matrixElements(1,iorb,2), 1, tmbmix%wfnmd%coeff(1,iorb), 1)
          end do
      !else
      !    if(iproc==0) write(*,*) "don't copy coeffs"
      !end if
      infoCoeff=0
 
      ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
      if(iproc==0) then
          write(*,'(1x,a)') '-------------------------------------------------'
          write(*,'(1x,a)') 'some selected eigenvalues:'
          do iorb=max(orbs%norb-8,1),min(orbs%norb+8,tmbmix%orbs%norb)
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
      !call dcopy(tmbmix%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      call optimize_coeffs(iproc, nproc, orbs, matrixElements(1,1,1), overlapmatrix, tmbmix, ldiis_coeff, fnrm)
  end if

  call calculate_density_kernel(iproc, nproc, tmbmix%wfnmd%ld_coeff, orbs, tmbmix%orbs, &
       tmbmix%wfnmd%coeff, tmbmix%wfnmd%density_kernel, ovrlp)

  ! Calculate the band structure energy with matrixElements instead of wfnmd%coeff sue to the problem mentioned
  ! above (wrong size of wfnmd%coeff)
  ebs=0.d0
  do jorb=1,tmbmix%orbs%norb
      do korb=1,jorb
          tt = tmbmix%wfnmd%density_kernel(korb,jorb)*matrixElements(korb,jorb,1)
          if(korb/=jorb) tt=2.d0*tt
          ebs = ebs + tt
      end do
  end do

  ! If closed shell multiply by two.
  if(orbs%nspin==1) ebs=2.d0*ebs


  ! Project the lb coefficients on the smaller subset
  if(tmbmix%wfnmd%bs%use_derivative_basis) then
      inc=4
      do iorb=1,orbs%norb
          jjorb=1
          do jorb=1,tmbmix%orbs%norb,inc
              tt=0.d0
              do korb=1,tmbmix%orbs%norb
                  tt = tt + tmbmix%wfnmd%coeff(korb,iorb)*overlapmatrix(korb,jorb)
              end do
              tmb%wfnmd%coeff_proj(jjorb,iorb)=tt
              jjorb=jjorb+1
          end do
      end do
  else
      do iorb=1,orbs%norb
          do jorb=1,tmbmix%orbs%norb
              tmb%wfnmd%coeff_proj(jorb,iorb)=tmbmix%wfnmd%coeff(jorb,iorb)
          end do
      end do
  end if



  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

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
    SIC, tmb, tmblarge2, lhphilarge2)
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
integer,intent(in):: iproc, nproc
integer,intent(out):: infoBasisFunctions
type(atoms_data), intent(in) :: at
type(orbitals_data):: orbs
real(8),dimension(3,at%nat):: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(GPU_pointers), intent(inout) :: GPU
real(8),intent(out):: trH, fnrm
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(localizedDIISParameters),intent(inout):: ldiis
type(DFT_wavefunction),target,intent(inout):: tmb
type(SIC_data) :: SIC !<parameters for the SIC methods
type(DFT_wavefunction),target,intent(inout):: tmblarge2
real(8),dimension(:),pointer,intent(inout):: lhphilarge2

! Local variables
real(8):: trHold, fnrmMax, meanAlpha, gnrm_in, gnrm_out, ediff, noise
integer:: iorb, consecutive_rejections,istat,istart,ierr,it,iall,ilr,jorb,nsatur
real(8),dimension(:),allocatable:: alpha,fnrmOldArr,alphaDIIS
real(8),dimension(:,:),allocatable:: ovrlp
logical:: emergency_exit, overlap_calculated
character(len=*),parameter:: subname='getLocalizedBasis'
real(8),dimension(:),pointer:: lhphi, lhphiold, lphiold
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
      if (tmblarge2%orbs%npsidim_orbs > 0) call to_zero(tmblarge2%orbs%npsidim_orbs,tmblarge2%psi(1))
      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, &
           tmb%psi, tmblarge2%psi)
      if(it==1) then
          call local_potential_dimensions(tmblarge2%lzd,tmblarge2%orbs,denspot%dpbox%ngatherarr(0,1))
          call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
               tmblarge2%comgp%nrecvbuf, tmblarge2%comgp%recvbuf, tmblarge2%comgp)
          call full_local_potential(iproc,nproc,tmblarge2%orbs,tmblarge2%Lzd,2,denspot%dpbox,denspot%rhov,&
               denspot%pot_work,tmblarge2%comgp)
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


      call calculate_energy_and_gradient_linear(iproc, nproc, it, &
           tmb%wfnmd%density_kernel, &
           ldiis, consecutive_rejections, &
           fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, gnrm_in, gnrm_out, &
           meanAlpha, emergency_exit, &
           tmb, lhphi, lhphiold, &
           tmblarge2, lhphilarge2, overlap_calculated, ovrlp)

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
      noise=tmb%wfnmd%bs%gnrm_mult*fnrm*tmb%orbs%norb
      if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY .and. &
          !(trH-trH_old)/fnrm <0.d0 .and. abs((trH-trH_old)/fnrm)<2.0d-2) then
          ediff<0.d0 .and. abs(ediff) < noise) then
          nsatur=nsatur+1
      else
          nsatur=0
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
          if(nsatur>=tmb%wfnmd%bs%nsatur_inner) then
              if(iproc==0) then
                  write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                      write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
                  if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                      write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=it
          else if(it>=tmb%wfnmd%bs%nit_basis_optimization) then
          !!if(it>=tmb%wfnmd%bs%nit_basis_optimization) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          else if(emergency_exit) then
              if(iproc==0) then
                  write(*,'(1x,a,i0,a)') 'WARNING: emergency exit after ',it, &
                      ' iterations to keep presumably good TMBs before they deteriorate too much.'
                  write (*,'(1x,a,2es15.7,f12.7)') '>>WRONG OUTPUT<< Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=-1
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          exit iterLoop
      end if
      trH_old=trH


      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lhphi, lphiold, alpha, trH, meanAlpha, alphaDIIS)

      !if(Iproc==0) WRITE(*,*) 'WARNING: NO RECONSTRUCTION OF KERNEL'
      call reconstruct_kernel(iproc, nproc, orbs, tmb, ovrlp, overlap_calculated, tmb%wfnmd%density_kernel)


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


    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine improveOrbitals(iproc, nproc, it, tmb, ldiis, lhphi, alpha)
use module_base
use module_types
use module_interfaces, except_this_one => improveOrbitals
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, it
type(DFT_wavefunction),intent(inout):: tmb
type(localizedDIISParameters),intent(inout):: ldiis
real(8),dimension(tmb%wfnmd%nphi),intent(in):: lhphi
real(8),dimension(tmb%orbs%norbp),intent(in):: alpha

! Local variables
integer:: istart, iorb, iiorb, ilr, ncount, owa, owanext

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








!!$subroutine transformHam(iproc, nproc, orbs, comms, phi, hphi, HamSmall)
!!$!
!!$! Purpose:
!!$! =======
!!$!   Builds the Hamiltonian in the basis of the localized basis functions phi. To do so, it gets all basis
!!$!   functions |phi_i> and H|phi_i> and then calculates H_{ij}=<phi_i|H|phi_j>. The basis functions phi are
!!$!   provided in the transposed form.
!!$!
!!$! Calling arguments:
!!$! ==================
!!$!   Input arguments:
!!$!   ----------------
!!$!     iproc      process ID
!!$!     nproc      total number of processes
!!$!     orbs       type describing the basis functions psi
!!$!     comms      type containing the communication parameters for the physical orbitals phi
!!$!     phi        basis functions 
!!$!     hphi       the Hamiltonian applied to the basis functions 
!!$!   Output arguments:
!!$!   -----------------
!!$!     HamSmall   Hamiltonian in small basis
!!$!
!!$use module_base
!!$use module_types
!!$implicit none
!!$
!!$! Calling arguments
!!$integer,intent(in):: iproc, nproc
!!$type(orbitals_data), intent(in) :: orbs
!!$type(communications_arrays), intent(in) :: comms
!!$real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb), intent(in) :: phi, hphi
!!$real(8),dimension(orbs%norb,orbs%norb),intent(out):: HamSmall
!!$
!!$! Local variables
!!$integer:: istat, ierr, nvctrp, iall
!!$real(8),dimension(:,:),allocatable:: HamTemp
!!$character(len=*),parameter:: subname='transformHam'
!!$
!!$
!!$
!!$  ! Allocate a temporary array if there are several MPI processes
!!$  if(nproc>1) then
!!$      allocate(HamTemp(orbs%norb,orbs%norb), stat=istat)
!!$      call memocc(istat, HamTemp, 'HamTemp', subname)
!!$  end if
!!$  
!!$  ! nvctrp is the amount of each phi hold by the current process
!!$  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!$  
!!$  ! Build the Hamiltonian. In the parallel case, each process writes its Hamiltonian in HamTemp
!!$  ! and a mpi_allreduce sums up the contribution from all processes.
!!$  if(nproc==1) then
!!$      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
!!$                 hphi(1,1), nvctrp, 0.d0, HamSmall(1,1), orbs%norb)
!!$  else
!!$      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
!!$                 hphi(1,1), nvctrp, 0.d0, HamTemp(1,1), orbs%norb)
!!$  end if
!!$  if(nproc>1) then
!!$      call mpi_allreduce(HamTemp(1,1), HamSmall(1,1), orbs%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!$  end if
!!$  
!!$  if(nproc>1) then
!!$     iall=-product(shape(HamTemp))*kind(HamTemp)
!!$     deallocate(HamTemp,stat=istat)
!!$     call memocc(istat, iall, 'HamTemp', subname)
!!$  end if
!!$
!!$end subroutine transformHam




!!!subroutine diagonalizeHamiltonian(iproc, nproc, orbs, HamSmall, eval)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!!!!   the same result. This is done by requiring that the first entry of each vector
!!!!   is positive.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc     process ID
!!!!     nproc     number of MPI processes
!!!!     orbs      type describing the physical orbitals psi
!!!!   Input / Putput arguments
!!!!     HamSmall  on input: the Hamiltonian
!!!!               on exit: the eigenvectors
!!!!   Output arguments
!!!!     eval      the associated eigenvalues 
!!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(inout) :: orbs
!!!real(8),dimension(orbs%norb, orbs%norb):: HamSmall
!!!real(8),dimension(orbs%norb):: eval
!!!
!!!! Local variables
!!!integer:: lwork, info, istat, iall, i, iorb, jorb
!!!real(8),dimension(:),allocatable:: work
!!!character(len=*),parameter:: subname='diagonalizeHamiltonian'
!!!
!!!  ! Get the optimal work array size
!!!  lwork=-1 
!!!  allocate(work(1), stat=istat)
!!!  call memocc(istat, work, 'work', subname)
!!!  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  lwork=work(1) 
!!!
!!!  ! Deallocate the work array ane reallocate it with the optimal size
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!  call memocc(istat, work, 'work', subname)
!!!
!!!  ! Diagonalize the Hamiltonian
!!!  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!
!!!  ! Deallocate the work array.
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!  
!!!  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!!  ! the first entry of each vector is positive.
!!!  do iorb=1,orbs%norb
!!!      if(HamSmall(1,iorb)<0.d0) then
!!!          do jorb=1,orbs%norb
!!!              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!!          end do
!!!      end if
!!!  end do
!!!
!!!
!!!end subroutine diagonalizeHamiltonian



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
integer:: iproc, nproc, nsubmax
type(orbitals_data), intent(inout) :: orbs
real(8),dimension(orbs%norb, orbs%norb),intent(inout) :: HamSmall
real(8),dimension(orbs%norb, orbs%norb),intent(in) :: ovrlp
real(8),dimension(orbs%norb),intent(out) :: eval

! Local variables
integer:: lwork, info, istat, iall, i, iorb, jorb, nsub
real(8),dimension(:),allocatable:: work
real(8),dimension(:,:),allocatable:: ham_band, ovrlp_band
character(len=*),parameter:: subname='diagonalizeHamiltonian'

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
  lwork=work(1) 

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




!!!subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!!!!
!!!! Purpose:
!!!! =======
!!!!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!!!!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc      process ID
!!!!     nproc      total number of processes
!!!!     orbs       type describing the physical orbitals psi
!!!!     orbsLIN    type describing the basis functions phi
!!!!     comms      type containing the communication parameters for the physical orbitals psi
!!!!     commsLIN   type containing the communication parameters for the basis functions phi
!!!!     phi        the basis functions 
!!!!     HamSmall   the  Hamiltonian matrix
!!!!   Output arguments:
!!!!   -----------------
!!!!     psi        the physical orbitals 
!!!!
!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(in) :: orbs
!!!type(orbitals_data), intent(in) :: orbsLIN
!!!type(communications_arrays), intent(in) :: comms
!!!type(communications_arrays), intent(in) :: commsLIN
!!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
!!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
!!!real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall
!!!
!!!! Local variables
!!!integer:: nvctrp
!!!
!!!
!!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, HamSmall(1,1), &
!!!             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
!!!  
!!!
!!!end subroutine buildWavefunction
!!
!!
!!

!!
!!!subroutine buildWavefunctionModified(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, coeff)
!!!
!!!!
!!!! Purpose:
!!!! =======
!!!!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!!!!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc      process ID
!!!!     nproc      total number of processes
!!!!     orbs       type describing the physical orbitals psi
!!!!     orbsLIN    type describing the basis functions phi
!!!!     comms      type containing the communication parameters for the physical orbitals psi
!!!!     commsLIN   type containing the communication parameters for the basis functions phi
!!!!     phi        the basis functions 
!!!!     coeff      the coefficients for the linear combination
!!!!   Output arguments:
!!!!   -----------------
!!!!     psi        the physical orbitals 
!!!!
!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(in) :: orbs
!!!type(orbitals_data), intent(in) :: orbsLIN
!!!type(communications_arrays), intent(in) :: comms
!!!type(communications_arrays), intent(in) :: commsLIN
!!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
!!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
!!!real(8),dimension(orbsLIN%norb,orbs%norb):: coeff
!!!
!!!! Local variables
!!!integer:: nvctrp
!!!
!!!
!!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, coeff(1,1), &
!!!             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
!!!  
!!!
!!!end subroutine buildWavefunctionModified






subroutine apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, confdatarr, hx, &
           psi, centralLocreg, vpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_orbitaldependent_potential
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, centralLocreg
type(atoms_data),intent(in):: at
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(3,at%nat),intent(in):: rxyz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),intent(in):: hx
!real(8),dimension(lzd%lpsidimtot),intent(in):: psi  !!!! ATENTION, intent should be in !
!real(8),dimension(lzd%lpsidimtot),intent(inout):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, vpsir
type(workarr_precond) :: work
type(workarrays_quartic_convolutions):: work_conv
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_orbitaldependent_potential'



  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), vpsi(1))
  ist_c=1
  ist_f=1
  do iorb=1,orbs%norbp
      iiorb=iorb+orbs%isorb
      ilr = orbs%inwhichlocreg(iiorb)
      if(centralLocreg<0) then
          !icenter=lin%orbs%inWhichLocregp(iorb)
          icenter=orbs%inWhichLocreg(iiorb)
      else
          icenter=centralLocreg
      end if
      !components of the potential
      npot=orbs%nspinor
      if (orbs%nspinor == 2) npot=1
      ist_f=ist_f+lzd%llr(ilr)%wfd%nvctr_c
      call allocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

      call uncompress_for_quartic_convolutions(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
           scal, psi(ist_c), psi(ist_f), &
           work_conv)

      if(confdatarr(iorb)%potorder==4) then
          call ConvolQuartic4(iproc,nproc,lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
               work_conv%y_c, work_conv%y_f)
      else if(confdatarr(iorb)%potorder==6) then
          call ConvolSextic(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
               work_conv%y_c, work_conv%y_f)
      else
          stop 'wronf conf pot'
      end if

      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(Ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
           scal, work_conv%y_c, work_conv%y_f, vpsi(ist_c), vpsi(ist_f))

      call deallocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

      ist_f = ist_f + 7*lzd%llr(ilr)%wfd%nvctr_f
      ist_c = ist_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

  end do



end subroutine apply_orbitaldependent_potential




!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in):: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in):: keyv_c
  integer,dimension(mseg_f),intent(in):: keyv_f
  integer,dimension(2,mseg_c),intent(in):: keyg_c
  integer,dimension(2,mseg_f),intent(in):: keyg_f
  real(wp),dimension(0:3),intent(in):: scal
  real(wp),dimension(mvctr_c),intent(in):: psi_c
  real(wp),dimension(7,mvctr_f),intent(in):: psi_f
  type(workarrays_quartic_convolutions),intent(out):: work
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





function dfactorial(n)
implicit none

! Calling arguments
integer,intent(in):: n
real(8):: dfactorial

! Local variables
integer:: i

  dfactorial=1.d0
  do i=1,n
      dfactorial=dfactorial*dble(i)
  end do

end function dfactorial




subroutine build_new_linear_combinations(iproc, nproc, lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
use module_base
use module_types
implicit none

!Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
integer,intent(in):: nrecvbuf
real(8),dimension(nrecvbuf),intent(in):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(in):: omat
logical,intent(in):: reset
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: lphi

! Local variables
integer:: ist, jst, ilrold, iorb, iiorb, ilr, ncount, jorb, jjorb, i, ldim, ind, indout, gdim, iorbref, m, ii
integer:: istart, iend, iseg, start, ifine, igrid, irecv, kseg, kold, kstart, kend 
real(8):: tt

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



!!!subroutine get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
!!!           confdatarr, hx, psi, potmat)
!!!use module_base
!!!use module_types
!!!use module_interfaces, eccept_this_one => get_potential_matrices
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(atoms_data),intent(in):: at
!!!type(orbitals_data),intent(in):: orbs
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(overlapParameters),intent(inout):: op
!!!type(p2pComms),intent(inout):: comon
!!!type(matrixDescriptors),intent(in):: mad
!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!real(8),intent(in):: hx
!!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!!real(8),dimension(orbs%norb,orbs%norb,at%nat),intent(out):: potmat
!!!
!!!! Local variables
!!!integer:: iorb, ilr, ilrold, istat, iall
!!!real(8),dimension(:,:),allocatable:: ttmat
!!!real(8):: tt1, tt2, tt3, tt4, tt5
!!!real(8),dimension(:),allocatable:: vpsi
!!!character(len=*),parameter:: subname='get_potential_matrices'
!!!
!!!allocate(vpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!call memocc(istat, vpsi, 'vpsi', subname)
!!!
!!!
!!!ilrold=-1
!!!do iorb=1,orbs%norb
!!!    ilr=orbs%inwhichlocreg(iorb)
!!!    if(ilr==ilrold) cycle
!!!    call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!         confdatarr, hx, psi, ilr, vpsi)
!!!
!!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, vpsi, comon%nsendBuf, comon%sendBuf)
!!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, vpsi, comon, tt1, tt2)
!!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!    !deallocate(ttmat)
!!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, psi, vpsi, mad, potmat(1,1,ilr))
!!!    ilrold=ilr
!!!    
!!!end do
!!!
!!!iall=-product(shape(vpsi))*kind(vpsi)
!!!deallocate(vpsi, stat=istat)
!!!call memocc(istat, iall, 'vpsi', subname)
!!!
!!!
!!!
!!!end subroutine get_potential_matrices





subroutine apply_position_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, xpsi, ypsi, zpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_position_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, order
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: xpsi, ypsi, zpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_position_operators'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!interface
!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
!!     hxh, hyh, hzh, dir, &
!!     ibyyzz_r) !optional
!!use module_base
!!implicit none
!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!real(8),intent(in):: hxh, hyh, hzh
!!character(len=1),intent(in):: dir
!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!end subroutine
!!end interface

  ishift=(/0,0,0/)

  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), xpsi(1))
  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), ypsi(1))
  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), zpsi(1))
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

     iiorb=orbs%isorb+iorb
     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psirx,'psirx',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)

     allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psiry,'psiry',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)

     allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psirz,'psirz',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)

     !!do i_stat=1,Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
     !!    write(1000+iproc,'(i9,es18.7,i9)') i_stat, psir(i_stat,1), Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
     !!end do
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*hx
     hyh=.5d0*hy
     hzh=.5d0*hz
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     if(lzd%llr(ilr)%geocode == 'F')then
        call position_operators(lzd%Glr%d%n1i,lzd%Glr%d%n2i,lzd%Glr%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, order, psirx, psiry, psirz, &
                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     else
        call position_operators(lzd%Glr%d%n1i,lzd%Glr%d%n2i,lzd%Glr%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                             psir, order, psirx, psiry, psirz, &
                             confdatarr(iorb)) !optional
     end if

     call isf_to_daub(lzd%llr(ilr), work_sr, psirx, xpsi(1+oidx))
     call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
     call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))

     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     i_all=-product(shape(psirx))*kind(psirx)
     deallocate(psirx,stat=i_stat)
     call memocc(i_stat,i_all,'psirx',subname)

     i_all=-product(shape(psiry))*kind(psiry)
     deallocate(psiry,stat=i_stat)
     call memocc(i_stat,i_all,'psiry',subname)

     i_all=-product(shape(psirz))*kind(psirz)
     deallocate(psirz,stat=i_stat)
     call memocc(i_stat,i_all,'psirz',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine apply_position_operators





subroutine position_operators(Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
     psirx, psiry, psirz, &
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(out) :: psirx, psiry, psirz !< x,y,z operator applied to real-space wfn in lr
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: ii1,ii2,ii3,i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, potx, poty, potz

  write(*,*) 'in position_operators'


  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(psir,psirx,psiry,psirz,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,order,Gn1i,Gn2i,Gn3i)&
  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,potx,poty,potz)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psirx(i1,i2,i3,ispinor)=0.0_wp 
             psiry(i1,i2,i3,ispinor)=0.0_wp 
             psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psirx(i1,i2,i3,ispinor)=0.0_wp 
              psiry(i1,i2,i3,ispinor)=0.0_wp 
              psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     write(*,'(a,3i9)') 'confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3)', &
                confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3)
     do ispinor=1,nspinor
        !$omp do
        do ii3=i3s,i3e
           i3=mod(ii3+confdata%ioffset(3)-1,Gn3i)+1
           do ii2=i2s,i2e
              i2=mod(ii2+confdata%ioffset(2)-1,Gn2i)+1
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,ii2-15,ii3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,ii2-15,ii3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do ii1=i1st,i1et
                 i1=mod(ii1+confdata%ioffset(1)-1,Gn1i)+1
                 psir1=psir(ii1,ii2,ii3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 potx=(confdata%hh(1)*real(i1,wp))**order
                 poty=(confdata%hh(2)*real(i2,wp))**order
                 potz=(confdata%hh(3)*real(i3,wp))**order

                 ttx=potx*psir1
                 tty=poty*psir1
                 ttz=potz*psir1

                 psirx(ii1,ii2,ii3,ispinor)=ttx
                 psiry(ii1,ii2,ii3,ispinor)=tty
                 psirz(ii1,ii2,ii3,ispinor)=ttz
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE position_operators







subroutine commutator(norb, A, B, res)
implicit none

! Calling arguments
integer,intent(in):: norb
real(8),dimension(norb,norb),intent(in):: A, B
real(8),dimension(norb,norb),intent(out):: res

! Local variables
real(8),dimension(norb,norb):: AB, BA

call dgemm('n', 'n', norb, norb, norb, 1.d0, A, norb, B, norb, 0.d0, AB, norb)
call dgemm('n', 'n', norb, norb, norb, 1.d0, B, norb, A, norb, 0.d0, BA, norb)
res=AB-BA

end subroutine commutator


subroutine apply_r_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, vpsi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_r_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, order
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_orbitaldependent_potential'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

  ishift=(/0,0,0/)

  !xpsi=0.d0
  !ypsi=0.d0
  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), vpsi(1))
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)

     iiorb=orbs%isorb+iorb
     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !!allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirx,'psirx',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)

     !!allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psiry,'psiry',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)

     !!allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     !!call memocc(i_stat,psirz,'psirz',subname)
     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !!!write(*,*) 'WARNING DEBUG in r_operator'
     !!!psir=1.d0/sqrt(dble(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i))
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*hx
     hyh=.5d0*hy
     hzh=.5d0*hz
     !icenter=confinementCenter(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     if(lzd%llr(ilr)%geocode == 'F') then
        call r_operator(lzd%Glr%d%n1i, lzd%Glr%d%n2i, lzd%Glr%d%n3i, &
                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                        ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                        psir, order, &
                        confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     else
        call r_operator(lzd%Glr%d%n1i, lzd%Glr%d%n2i, lzd%Glr%d%n3i, &
                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
                        ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
                        psir, order, &
                        confdatarr(iorb)) !optional
     end if

     call isf_to_daub(lzd%llr(ilr), work_sr, psir, vpsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
     !!call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))

     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))

     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)


     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     !!i_all=-product(shape(psirx))*kind(psirx)
     !!deallocate(psirx,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirx',subname)

     !!i_all=-product(shape(psiry))*kind(psiry)
     !!deallocate(psiry,stat=i_stat)
     !!call memocc(i_stat,i_all,'psiry',subname)

     !!i_all=-product(shape(psirz))*kind(psirz)
     !!deallocate(psirz,stat=i_stat)
     !!call memocc(i_stat,i_all,'psirz',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine apply_r_operators








subroutine r_operator(Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ii1,ii2,ii3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(wp):: ttx, tty, ttz, potx, poty, potz

  if(order/=1 .and. order/=2) stop 'wrong order'

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,order,Gn1i,Gn2i,Gn3i)&
  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,potx,poty,potz)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
             !!psiry(i1,i2,i3,ispinor)=0.0_wp 
             !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
      stop 'not yet implemented for nspinor==4!'
     !!!$omp do
     !!do i3=i3s,i3e
     !!   do i2=i2s,i2e
     !!      !thanks to the optional argument the conditional is done at compile time
     !!      if (present(ibyyzz_r)) then
     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!      else
     !!         i1st=i1s
     !!         i1et=i1e
     !!      end if
     !!      !no need of setting up to zero values outside wavefunction bounds
     !!      do i1=i1st,i1et
     !!         !wavefunctions
     !!         psir1=psir(i1,i2,i3,1)
     !!         psir2=psir(i1,i2,i3,2)
     !!         psir3=psir(i1,i2,i3,3)
     !!         psir4=psir(i1,i2,i3,4)
     !!         !potentials + confining term
     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!         !diagonal terms
     !!         tt11=pot1*psir1 !p1
     !!         tt22=pot1*psir2 !p2
     !!         tt33=pot4*psir3 !p3
     !!         tt44=pot4*psir4 !p4
     !!         !Rab*Rb
     !!         tt13=pot2*psir3 !p1
     !!         !Iab*Ib
     !!         tt14=pot3*psir4 !p1
     !!         !Rab*Ib
     !!         tt23=pot2*psir4 !p2
     !!         !Iab*Rb
     !!         tt24=pot3*psir3 !p2
     !!         !Rab*Ra
     !!         tt31=pot2*psir1 !p3
     !!         !Iab*Ia
     !!         tt32=pot3*psir2 !p3
     !!         !Rab*Ia
     !!         tt41=pot2*psir2 !p4
     !!         !Iab*Ra
     !!         tt42=pot3*psir1 !p4

     !!         !value of the potential energy
     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!         !wavefunction update
     !!         !p1=h1p1+h2p3-h3p4
     !!         !p2=h1p2+h2p4+h3p3
     !!         !p3=h2p1+h3p2+h4p3
     !!         !p4=h2p2-h3p1+h4p4
     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!      end do
     !!   end do
     !!end do
     !!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do ii3=i3s,i3e
           i3=mod(ii3+confdata%ioffset(3)-1,Gn3i)+1
           do ii2=i2s,i2e
              i2=mod(ii2+confdata%ioffset(2)-1,Gn2i)+1
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,ii2-15,ii3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,ii2-15,ii3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do ii1=i1st,i1et
                 i1=mod(ii1+confdata%ioffset(1)-1,Gn1i)+1
                 psir1=psir(ii1,ii2,ii3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
                 ttx=(confdata%hh(1)*real(i1,wp))**2
                 tty=(confdata%hh(2)*real(i2,wp))**2
                 ttz=(confdata%hh(3)*real(i3,wp))**2

                 tt = ttx+tty+ttz

                 if(order==1) then
                     tt=sqrt(tt)
                 end if

                 psir(ii1,ii2,ii3,ispinor)=tt*psir1
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel


END SUBROUTINE r_operator




subroutine update_confdatarr(lzd, orbs, locregCenter, confdatarr)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(3,lzd%nlr),intent(in):: locregCenter
  type(confpot_data),dimension(orbs%norbp),intent(inout):: confdatarr
  
  ! Local variables
  integer:: iorb, iiorb, ilr, icenter, nl1, nl2, nl3
  
  ! Update confdatarr...
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr=orbs%inWhichlocreg(iiorb)
     icenter=orbs%inwhichlocreg(iiorb)
     !confdatarr(iorb)%potorder=lin%confpotorder
     !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
     !confdatarr(iorb)%hh(1)=.5_gp*hx
     !confdatarr(iorb)%hh(2)=.5_gp*hy
     !confdatarr(iorb)%hh(3)=.5_gp*hz
     confdatarr(iorb)%rxyzConf(1:3)=locregCenter(1:3,icenter)
     call my_geocode_buffers(lzd%Llr(ilr)%geocode,nl1,nl2,nl3)
     confdatarr(iorb)%ioffset(1)=lzd%llr(ilr)%nsi1-nl1-1
     confdatarr(iorb)%ioffset(2)=lzd%llr(ilr)%nsi2-nl2-1
     confdatarr(iorb)%ioffset(3)=lzd%llr(ilr)%nsi3-nl3-1
  end do

end subroutine update_confdatarr



subroutine small_to_large_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, phismall, philarge)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzdsmall, lzdlarge
  type(orbitals_data),intent(in):: orbssmall, orbslarge
  real(8),dimension(orbssmall%npsidim_orbs),intent(in):: phismall
  real(8),dimension(orbslarge%npsidim_orbs),intent(out):: philarge
  
  ! Local variables
  integer:: ists, istl, iorb, ilr, ilrlarge, sdim, ldim, nspin
  
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

end subroutine small_to_large_locreg


subroutine large_to_small_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, philarge, phismall)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzdsmall, lzdlarge
  type(orbitals_data),intent(in):: orbssmall, orbslarge
  real(8),dimension(orbslarge%npsidim_orbs),intent(in):: philarge
  real(8),dimension(orbssmall%npsidim_orbs),intent(out):: phismall
  
  ! Local variables
  integer:: istl, ists, ilr, ilrlarge, ldim, gdim, iorb
  
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

end subroutine large_to_small_locreg



subroutine check_locregCenters(iproc, lzd, locregCenter, hx, hy, hz)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(local_zone_descriptors),intent(in):: lzd
  real(8),dimension(3,lzd%nlr),intent(in):: locregCenter
  real(8),intent(in):: hx, hy, hz
  
  ! Local variables
  integer:: ilr, ierr
  
  do ilr=1,lzd%nlr
      if( floor(locregCenter(1,ilr)/hx) < 0 .or. ceiling(locregCenter(1,ilr)/hx) > lzd%glr%d%n1 ) then
          if(iproc==0) then
              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                  ' is outside of box in x direction! Box limits=',0,lzd%glr%d%n1,&
                  ', center=',floor(locregCenter(1,ilr)/hx),ceiling(locregCenter(1,ilr)/hx)
          end if
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
      if( floor(locregCenter(2,ilr)/hy) < 0 .or. ceiling(locregCenter(2,ilr)/hy) > lzd%glr%d%n2 ) then
          if(iproc==0) then
              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                  'is outside of box in y direction! Box limits=',0,lzd%glr%d%n2,&
                  ', center=',floor(locregCenter(2,ilr)/hy),ceiling(locregCenter(2,ilr)/hy)
          end if
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
      if( floor(locregCenter(3,ilr)/hz) < 0 .or. ceiling(locregCenter(3,ilr)/hz) > lzd%glr%d%n3 ) then
          if(iproc==0) then
              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
                  'is outside of box in z direction! Box limits=',0,lzd%glr%d%n3,&
                  ', center=',floor(locregCenter(3,ilr)/hz),ceiling(locregCenter(3,ilr)/hz)
          end if
          call mpi_barrier(mpi_comm_world, ierr)
          stop
      end if
  end do

end subroutine check_locregCenters           







subroutine communicate_basis_for_density(iproc, nproc, lzd, llborbs, lphi, comsr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => communicate_basis_for_density
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: llborbs
  real(8),dimension(llborbs%npsidim_orbs),intent(in):: lphi
  type(p2pComms),intent(inout):: comsr
  
  ! Local variables
  integer:: ist, istr, iorb, iiorb, ilr, ierr
  type(workarr_sumrho):: w

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



subroutine update_kernel(norb, Umat, kernel)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  real(8),dimension(norb,norb),intent(in):: Umat
  real(8),dimension(norb,norb),intent(inout):: kernel
  
  ! Local variables
  integer:: iorb, jorb, korb, lorb, istat, iall
  real(8):: tt
  real(8),dimension(:,:),allocatable:: kernelold
  character(len=*),parameter:: subname='update_kernel'
  
  allocate(kernelold(norb,norb), stat=istat)
  call memocc(istat, kernelold, 'kernelold', subname)
  
  call dcopy(norb**2, kernel(1,1), 1, kernelold(1,1), 1)
  do iorb=1,norb
      do jorb=1,norb
          tt=0.d0
          do korb=1,norb
              do lorb=1,norb
                  tt=tt+kernelold(korb,lorb)*Umat(korb,iorb)*Umat(lorb,jorb)
                  !tt=tt+kernelold(korb,lorb)*Umat(iorb,korb)*Umat(jorb,lorb)
              end do
          end do
          kernel(jorb,iorb)=tt
      end do
  end do
  
  iall=-product(shape(kernelold))*kind(kernelold)
  deallocate(kernelold, stat=istat)
  call memocc(istat, iall, 'kernelold', subname)

end subroutine update_kernel



subroutine DIISorSD(iproc, nproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, it
  real(8),intent(in):: trH
  type(DFT_wavefunction),intent(inout):: tmbopt
  type(localizedDIISParameters),intent(inout):: ldiis
  real(8),dimension(tmbopt%orbs%norbp),intent(out):: alpha, alphaDIIS
  real(8),dimension(tmbopt%wfnmd%nphi),intent(out):: lphioldopt
  
  ! Local variables
  integer:: idsx, ii, offset, istdest, iorb, iiorb, ilr, ncount, istsource
  
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



subroutine flatten_at_boundaries(lzd, orbs, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi

  ! Local variables
  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ii1, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt
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

      r0=lzd%llr(ilr)%locrad**2/5.d0
      do i3=0,lzd%llr(ilr)%d%n3
          ii3=lzd%llr(ilr)%ns3+i3
          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
          do i2=0,lzd%llr(ilr)%d%n2
              ii2=lzd%llr(ilr)%ns2+i2
              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
              do i1=0,lzd%llr(ilr)%d%n1
                  ii1=lzd%llr(ilr)%ns1+i1
                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
                  rr=r1+r2+r3
                  tt=exp(-(rr-lzd%llr(ilr)%locrad**2/2.d0)/r0)
                  tt=min(tt,1.d0)
                  psig(i1,1,i2,1,i3,1)=tt*psig(i1,1,i2,1,i3,1)
                  psig(i1,2,i2,1,i3,1)=tt*psig(i1,2,i2,1,i3,1)
                  psig(i1,1,i2,2,i3,1)=tt*psig(i1,1,i2,2,i3,1)
                  psig(i1,2,i2,2,i3,1)=tt*psig(i1,2,i2,2,i3,1)
                  psig(i1,1,i2,1,i3,2)=tt*psig(i1,1,i2,1,i3,2)
                  psig(i1,2,i2,1,i3,2)=tt*psig(i1,2,i2,1,i3,2)
                  psig(i1,1,i2,2,i3,2)=tt*psig(i1,1,i2,2,i3,2)
                  psig(i1,2,i2,2,i3,2)=tt*psig(i1,2,i2,2,i3,2)
              end do
          end do
      end do

      call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
           0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           psig, psi(istc), psi(istf))

      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

      iall=-product(shape(psig))*kind(psig)
      deallocate(psig,stat=istat)
      call memocc(istat,iall,'psig',subname)


  end do

end subroutine flatten_at_boundaries



subroutine get_weighted_gradient(iproc, nproc, lzd, orbs, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi

  ! Local variables
  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii1, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt, gnrm
  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
  character(len=*),parameter:: subname='flatten_at_boundaries'
  

  istc=1
  istf=1
  gnrm=0.d0
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

      r0=lzd%llr(ilr)%locrad**2/3.d0
      do i3=0,lzd%llr(ilr)%d%n3
          ii3=lzd%llr(ilr)%ns3+i3
          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
          do i2=0,lzd%llr(ilr)%d%n2
              ii2=lzd%llr(ilr)%ns2+i2
              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
              do i1=0,lzd%llr(ilr)%d%n1
                  ii1=lzd%llr(ilr)%ns1+i1
                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
                  rr=r1+r2+r3
                  !tt=exp(-(rr-lzd%llr(ilr)%locrad**2/2.d0)/r0)
                  tt=exp(-rr/r0)
                  !tt=min(tt,1.d0)
                  gnrm=gnrm+tt*psig(i1,1,i2,1,i3,1)**2
                  gnrm=gnrm+tt*psig(i1,2,i2,1,i3,1)**2
                  gnrm=gnrm+tt*psig(i1,1,i2,2,i3,1)**2
                  gnrm=gnrm+tt*psig(i1,2,i2,2,i3,1)**2
                  gnrm=gnrm+tt*psig(i1,1,i2,1,i3,2)**2
                  gnrm=gnrm+tt*psig(i1,2,i2,1,i3,2)**2
                  gnrm=gnrm+tt*psig(i1,1,i2,2,i3,2)**2
                  gnrm=gnrm+tt*psig(i1,2,i2,2,i3,2)**2
              end do
          end do
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

  call mpiallred(gnrm, 1, mpi_sum, mpi_comm_world, ierr)
  gnrm=gnrm/dble(orbs%norb)
  if(iproc==0) write(*,*) 'weighted grnm',gnrm

end subroutine get_weighted_gradient




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



subroutine cut_at_boundaries(lzd, orbs, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi

  ! Local variables
  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ii1, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt
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

      r0 = (lzd%llr(ilr)%locrad-8.d0*lzd%hgrids(1))**2
      !r0 = (lzd%llr(ilr)%locrad-2.d0)**2
      do i3=0,lzd%llr(ilr)%d%n3
          ii3=lzd%llr(ilr)%ns3+i3
          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
          do i2=0,lzd%llr(ilr)%d%n2
              ii2=lzd%llr(ilr)%ns2+i2
              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
              do i1=0,lzd%llr(ilr)%d%n1
                  ii1=lzd%llr(ilr)%ns1+i1
                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
                  rr=r1+r2+r3
                  !write(999,'(5es14.3)') r1, r2, r3, rr, r0
                  if(rr>=r0) then
                      psig(i1,1,i2,1,i3,1)=0.d0
                      psig(i1,2,i2,1,i3,1)=0.d0
                      psig(i1,1,i2,2,i3,1)=0.d0
                      psig(i1,2,i2,2,i3,1)=0.d0
                      psig(i1,1,i2,1,i3,2)=0.d0
                      psig(i1,2,i2,1,i3,2)=0.d0
                      psig(i1,1,i2,2,i3,2)=0.d0
                      psig(i1,2,i2,2,i3,2)=0.d0
                  else
                      !write(*,*) 'not zero'
                  end if
              end do
          end do
      end do

      call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
           0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           psig, psi(istc), psi(istf))

      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

      iall=-product(shape(psig))*kind(psig)
      deallocate(psig,stat=istat)
      call memocc(istat,iall,'psig',subname)


  end do

end subroutine cut_at_boundaries






subroutine cut_at_boundaries2(lr, orbs, hx, hy, hz, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  type(orbitals_data),intent(in):: orbs
  real(8),intent(in):: hx, hy, hz
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi

  ! Local variables
  integer:: istc, istf, iorb, i1, i2, i3, istat, iall, ii1, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt
  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
  character(len=*),parameter:: subname='flatten_at_boundaries'
  

  istc=1
  istf=1

  allocate(psig(0:lr%d%n1,2,0:lr%d%n2,2,0:lr%d%n3,2), stat=istat)
  call memocc(istat, psig, 'psig', subname)
  call to_zero(8*(lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), psig(0,1,0,1,0,1))

  istf = istf + lr%wfd%nvctr_c
  call uncompress(lr%d%n1, lr%d%n2, lr%d%n3, &
       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       psi(istc), psi(istf), psig)

  r0 = (lr%locrad-8.d0*hx)**2
  !r0 = (lr%locrad-2.d0)**2
  do i3=0,lr%d%n3
      ii3=lr%ns3+i3
      r3 = (dble(ii3)*hz-lr%locregCenter(3))**2
      do i2=0,lr%d%n2
          ii2=lr%ns2+i2
          r2 = (dble(ii2)*hy-lr%locregCenter(2))**2
          do i1=0,lr%d%n1
              ii1=lr%ns1+i1
              r1 = (dble(ii1)*hx-lr%locregCenter(1))**2
              rr=r1+r2+r3
              !write(999,'(5es14.3)') r1, r2, r3, rr, r0
              if(rr>=r0) then
                  psig(i1,1,i2,1,i3,1)=0.d0
                  psig(i1,2,i2,1,i3,1)=0.d0
                  psig(i1,1,i2,2,i3,1)=0.d0
                  psig(i1,2,i2,2,i3,1)=0.d0
                  psig(i1,1,i2,1,i3,2)=0.d0
                  psig(i1,2,i2,1,i3,2)=0.d0
                  psig(i1,1,i2,2,i3,2)=0.d0
                  psig(i1,2,i2,2,i3,2)=0.d0
              else
                  !write(*,*) 'not zero'
              end if
          end do
      end do
  end do

  call compress(lr%d%n1, lr%d%n2, &
       0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, &
       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  &
       psig, psi(istc), psi(istf))

  istf = istf + 7*lr%wfd%nvctr_f
  istc = istc + lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f

  iall=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=istat)
  call memocc(istat,iall,'psig',subname)



end subroutine cut_at_boundaries2










subroutine flatten_at_boundaries2(lr, orbs, hx, hy, hz, psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  type(orbitals_data),intent(in):: orbs
  real(8),intent(in):: hx, hy, hz
  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi

  ! Local variables
  integer:: istc, istf, iorb, i1, i2, i3, istat, iall, ii1, ii2, ii3
  real(8):: r0, r1, r2, r3, rr, tt
  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
  character(len=*),parameter:: subname='flatten_at_boundaries'
  

  istc=1
  istf=1

  allocate(psig(0:lr%d%n1,2,0:lr%d%n2,2,0:lr%d%n3,2), stat=istat)
  call memocc(istat, psig, 'psig', subname)
  call to_zero(8*(lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), psig(0,1,0,1,0,1))

  istf = istf + lr%wfd%nvctr_c
  call uncompress(lr%d%n1, lr%d%n2, lr%d%n3, &
       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       psi(istc), psi(istf), psig)

  r0=lr%locrad**2/5.d0
  !r0 = (lr%locrad-2.d0)**2
  do i3=0,lr%d%n3
      ii3=lr%ns3+i3
      r3 = (dble(ii3)*hz-lr%locregCenter(3))**2
      do i2=0,lr%d%n2
          ii2=lr%ns2+i2
          r2 = (dble(ii2)*hy-lr%locregCenter(2))**2
          do i1=0,lr%d%n1
              ii1=lr%ns1+i1
              r1 = (dble(ii1)*hx-lr%locregCenter(1))**2
              rr=r1+r2+r3
              tt=exp(-(rr-lr%locrad**2/2.d0)/r0)
              if(tt<1.d0) then
                  psig(i1,1,i2,1,i3,1)=tt*psig(i1,1,i2,1,i3,1)
                  psig(i1,2,i2,1,i3,1)=tt*psig(i1,2,i2,1,i3,1)
                  psig(i1,1,i2,2,i3,1)=tt*psig(i1,1,i2,2,i3,1)
                  psig(i1,2,i2,2,i3,1)=tt*psig(i1,2,i2,2,i3,1)
                  psig(i1,1,i2,1,i3,2)=tt*psig(i1,1,i2,1,i3,2)
                  psig(i1,2,i2,1,i3,2)=tt*psig(i1,2,i2,1,i3,2)
                  psig(i1,1,i2,2,i3,2)=tt*psig(i1,1,i2,2,i3,2)
                  psig(i1,2,i2,2,i3,2)=tt*psig(i1,2,i2,2,i3,2)
              end if
          end do
      end do
  end do

  call compress(lr%d%n1, lr%d%n2, &
       0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, &
       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  &
       psig, psi(istc), psi(istf))

  istf = istf + 7*lr%wfd%nvctr_f
  istc = istc + lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f

  iall=-product(shape(psig))*kind(psig)
  deallocate(psig,stat=istat)
  call memocc(istat,iall,'psig',subname)



end subroutine flatten_at_boundaries2





subroutine get_both_gradients(iproc, nproc, lzd, orbs, psi, gnrm_in, gnrm_out)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
  real(8),intent(out):: gnrm_out, gnrm_in

  ! Local variables
  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii1, ii2, ii3, ipts_out, ipts_in
  real(8):: r0, r1, r2, r3, rr, tt, pts_in, pts_out
  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
  character(len=*),parameter:: subname='flatten_at_boundaries'
  

  istc=1
  istf=1
  gnrm_out=0.d0
  gnrm_in=0.d0
  pts_in=0.d0
  pts_out=0.d0
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

      ipts_in=0
      ipts_out=0
      !r0=(lzd%llr(ilr)%locrad-3.d0)**2
      r0=(lzd%llr(ilr)%locrad-8.d0*lzd%hgrids(1))**2
      do i3=0,lzd%llr(ilr)%d%n3
          ii3=lzd%llr(ilr)%ns3+i3
          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
          do i2=0,lzd%llr(ilr)%d%n2
              ii2=lzd%llr(ilr)%ns2+i2
              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
              do i1=0,lzd%llr(ilr)%d%n1
                  ii1=lzd%llr(ilr)%ns1+i1
                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
                  rr=r1+r2+r3
                  if(rr>r0) then
                      gnrm_out=gnrm_out+psig(i1,1,i2,1,i3,1)**2
                      gnrm_out=gnrm_out+psig(i1,2,i2,1,i3,1)**2
                      gnrm_out=gnrm_out+psig(i1,1,i2,2,i3,1)**2
                      gnrm_out=gnrm_out+psig(i1,2,i2,2,i3,1)**2
                      gnrm_out=gnrm_out+psig(i1,1,i2,1,i3,2)**2
                      gnrm_out=gnrm_out+psig(i1,2,i2,1,i3,2)**2
                      gnrm_out=gnrm_out+psig(i1,1,i2,2,i3,2)**2
                      gnrm_out=gnrm_out+psig(i1,2,i2,2,i3,2)**2
                      if(psig(i1,1,i2,1,i3,1)/=0.d0) then
                          ! point carries scaling function
                          pts_out=pts_out+1.d0
                          ipts_out=ipts_out+1
                      end if
                      if(psig(i1,2,i2,1,i3,1)/=0.d0) then
                          ! point carries wavelets
                          pts_out=pts_out+7.d0
                          ipts_out=ipts_out+7
                      end if
                  else
                      gnrm_in=gnrm_in+psig(i1,1,i2,1,i3,1)**2
                      gnrm_in=gnrm_in+psig(i1,2,i2,1,i3,1)**2
                      gnrm_in=gnrm_in+psig(i1,1,i2,2,i3,1)**2
                      gnrm_in=gnrm_in+psig(i1,2,i2,2,i3,1)**2
                      gnrm_in=gnrm_in+psig(i1,1,i2,1,i3,2)**2
                      gnrm_in=gnrm_in+psig(i1,2,i2,1,i3,2)**2
                      gnrm_in=gnrm_in+psig(i1,1,i2,2,i3,2)**2
                      gnrm_in=gnrm_in+psig(i1,2,i2,2,i3,2)**2
                      if(psig(i1,1,i2,1,i3,1)/=0.d0) then
                          ! point carries scaling function
                          pts_in=pts_in+1.d0
                          ipts_in=ipts_in+1
                      end if
                      if(psig(i1,2,i2,1,i3,1)/=0.d0) then
                          ! point carries wavelets
                          pts_in=pts_in+7.d0
                          ipts_in=ipts_in+7
                      end if
                  end if
              end do
          end do
      end do

      if (ipts_in+ipts_out /= lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f) then
          write(*,'(2(a,i0))') 'ERROR: ',ipts_in+ipts_out,&
                      ' = ipts_in+ipts_out /= lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f = ',&
                      lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          stop
      end if

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

  !if(pts_in>0) gnrm_in=gnrm_in/pts_in
  !if(pts_out>0) gnrm_out=gnrm_out/pts_out
  call mpiallred(gnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(gnrm_in, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(pts_in, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(pts_out, 1, mpi_sum, mpi_comm_world, ierr)
  gnrm_out=sqrt(gnrm_out/dble(orbs%norb))
  gnrm_in=sqrt(gnrm_in/dble(orbs%norb))
  if(iproc==0) write(*,'(a,5es14.4)') 'pts_in, pts_out, gnrm_in, gnrm_out, gnrm_in/gnrm_out', &
      pts_in, pts_out, gnrm_in, gnrm_out, gnrm_in/gnrm_out

end subroutine get_both_gradients





subroutine create_penalty_basis_function(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi)
use module_base
use module_types
use module_interfaces, except_this_one => apply_position_operators
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr, owa, owanext
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
real(8),dimension(:,:),allocatable:: psir
type(workarr_sumrho):: work_sr
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
character(len=*),parameter:: subname='apply_position_operators'
integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!interface
!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
!!     hxh, hyh, hzh, dir, &
!!     ibyyzz_r) !optional
!!use module_base
!!implicit none
!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!real(8),intent(in):: hxh, hyh, hzh
!!character(len=1),intent(in):: dir
!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!end subroutine
!!end interface

  ishift=(/0,0,0/)

  oidx = 0
  do iorb=1,orbs%norbp
     iiorb=orbs%isorb+iorb
     ilr = orbs%inwhichlocreg(iiorb)
     owa = orbs%onwhichatom(iiorb)
     if(iiorb<orbs%norb) then
         owanext=orbs%onwhichatom(iiorb+1)
     else
         owanext=lzd%nlr+1
     end if

     if(owa/=owanext) then
         ! last basis function

  
         !initialise the work arrays
         call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

         ! Wavefunction in real space
         allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
         call memocc(i_stat,psir,'psir',subname)
         call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

         call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)

         !!do i_stat=1,Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
         !!    write(1000+iproc,'(i9,es18.7,i9)') i_stat, psir(i_stat,1), Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
         !!end do
         !apply the potential to the psir wavefunction and calculate potential energy
         hxh=.5d0*hx
         hyh=.5d0*hy
         hzh=.5d0*hz
         !icenter=confinementCenter(iorb)
         !components of the potential
         npot=orbs%nspinor
         if (orbs%nspinor == 2) npot=1

         !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
         !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
         !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
         !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
         write(*,*) 'confdatarr(iorb)%prefac',confdatarr(iorb)%prefac
        call penalty_basis_function(lzd%llr(ilr)%d%n1i,lzd%llr(ilr)%d%n2i,lzd%llr(ilr)%d%n3i,&
             lzd%llr(ilr)%d%n1i,lzd%llr(ilr)%d%n2i,lzd%llr(ilr)%d%n3i,&
             ishift,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,&
             orbs%nspinor,psir(1,1),&
             confdata=confdatarr(iorb),ibyyzz_r=lzd%llr(ilr)%bounds%ibyyzz_r)

         call isf_to_daub(lzd%llr(ilr), work_sr, psir, psi(1+oidx))



         i_all=-product(shape(psir))*kind(psir)
         deallocate(psir,stat=i_stat)
         call memocc(i_stat,i_all,'psir',subname)


         call deallocate_work_arrays_sumrho(work_sr)

     end if

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine create_penalty_basis_function



subroutine penalty_basis_function(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
     confdata,ibyyzz_r) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4


  !write(*,*) 'present(confdata)', present(confdata)
  !write(*,*) 'confdata%prefac, confdata%potorder', confdata%prefac, confdata%potorder
  !write(*,*) 'n1ip*n2ip*n3ip', n1ip*n2ip*n3ip

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))


  !$omp parallel default(none)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)

!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
  !case without bounds


  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do


  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do


  !important part of the array
  if (nspinor==4) then
     stop 'not yet implemented for nspinor==4!'
     !!!!$omp do
     !!!do i3=i3s,i3e
     !!!   do i2=i2s,i2e
     !!!      !thanks to the optional argument the conditional is done at compile time
     !!!      if (present(ibyyzz_r)) then
     !!!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
     !!!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
     !!!      else
     !!!         i1st=i1s
     !!!         i1et=i1e
     !!!      end if
     !!!      !no need of setting up to zero values outside wavefunction bounds
     !!!      do i1=i1st,i1et
     !!!         !wavefunctions
     !!!         psir1=psir(i1,i2,i3,1)
     !!!         psir2=psir(i1,i2,i3,2)
     !!!         psir3=psir(i1,i2,i3,3)
     !!!         psir4=psir(i1,i2,i3,4)
     !!!         !potentials + confining term
     !!!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
     !!!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
     !!!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
     !!!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

     !!!         !diagonal terms
     !!!         tt11=pot1*psir1 !p1
     !!!         tt22=pot1*psir2 !p2
     !!!         tt33=pot4*psir3 !p3
     !!!         tt44=pot4*psir4 !p4
     !!!         !Rab*Rb
     !!!         tt13=pot2*psir3 !p1
     !!!         !Iab*Ib
     !!!         tt14=pot3*psir4 !p1
     !!!         !Rab*Ib
     !!!         tt23=pot2*psir4 !p2
     !!!         !Iab*Rb
     !!!         tt24=pot3*psir3 !p2
     !!!         !Rab*Ra
     !!!         tt31=pot2*psir1 !p3
     !!!         !Iab*Ia
     !!!         tt32=pot3*psir2 !p3
     !!!         !Rab*Ia
     !!!         tt41=pot2*psir2 !p4
     !!!         !Iab*Ra
     !!!         tt42=pot3*psir1 !p4

     !!!         !value of the potential energy
     !!!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
     !!!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

     !!!         !wavefunction update
     !!!         !p1=h1p1+h2p3-h3p4
     !!!         !p2=h1p2+h2p4+h3p3
     !!!         !p3=h2p1+h3p2+h4p3
     !!!         !p4=h2p2-h3p1+h4p4
     !!!         psir(i1,i2,i3,1)=tt11+tt13-tt14
     !!!         psir(i1,i2,i3,2)=tt22+tt23+tt24
     !!!         psir(i1,i2,i3,3)=tt33+tt31+tt32
     !!!         psir(i1,i2,i3,4)=tt44+tt41-tt42
     !!!      end do
     !!!   end do
     !!!end do
     !!!!$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              !write(*,'(a,6i9)') 'i1st, i1et, i2s, i2e, i3s, i3e', i1st, i1et, i2s, i2e, i3s, i3e
              do i1=i1st,i1et
                 pot1=cp(i1,i2,i3)
                 psir(i1,i2,i3,ispinor)=pot1
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  
  !$omp end parallel

contains
  
  !inline the definition of the confining potential
  real(wp) function cp(i1,i2,i3)
    implicit none
    integer, intent(in) :: i1,i2,i3
    !local variables
    real(wp) :: r2
    !to be sure that the conditional is executed at compile time
    if (present(confdata)) then
       r2=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1))**2 +&
            (confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2))**2 +&
            (confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3))**2 
       !if(r2>=81.d0) write(*,'(6i8,3es11.2,es13.4)') i1, i2, i3, confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3), confdata%rxyzConf(1), confdata%rxyzConf(2), confdata%rxyzConf(3), r2 

       cp=confdata%prefac*r2**(confdata%potorder/2)
    else
       cp=0.0_wp
    end if

  end function cp

END SUBROUTINE penalty_basis_function

!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2011 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Calculates the application of the Hamiltonian on the wavefunction. The hamiltonian can be self-consistent or not.
!! In the latter case, the potential should be given in the rhov array of denspot structure. 
!! Otherwise, rhov array is filled by the self-consistent density
subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,iscf,alphamix,ixc,&
     nlpspd,proj,rxyz,linflag,unblock_comms,GPU,wfn,&
     energs,rpnrm,xcstr)
  use module_base
  use module_types
  use module_interfaces, fake_name => psitohpsi
  use Poisson_Solver
  use m_ab6_mixing
  use yaml_output
  implicit none
  logical, intent(in) :: scf
  integer, intent(in) :: iproc,nproc,itrp,iscf,ixc,linflag,itwfn
  character(len=3), intent(in) :: unblock_comms
  real(gp), intent(in) :: alphamix
  type(atoms_data), intent(in) :: atoms
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(DFT_local_fields), intent(inout) :: denspot
  type(energy_terms), intent(inout) :: energs
  type(DFT_wavefunction), intent(inout) :: wfn
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  type(GPU_pointers), intent(inout) :: GPU  
  !real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(gp), intent(out) :: rpnrm
  real(gp), dimension(6), intent(out) :: xcstr
  !real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='psitohpsi'
  logical :: unblock_comms_den,unblock_comms_pot,whilepot,savefields,correcth
  integer :: nthread_max,ithread,nthread,irhotot_add,irho_add,ispin,i_all,i_stat
  !$ integer :: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads

  !in the default case, non local hamiltonian is done after potential creation
  whilepot=.true.

  !flag for saving the local fields (rho,vxc,vh)
  savefields= (iscf==SCF_KIND_GENERALIZED_DIRMIN)
  correcth=.false.
  !do not do that if rho_work is already associated
  if (savefields .and. associated(denspot%rho_work)) then
     !flag for correcting the hamiltonian (either is false or toggles savefields)
     correcth=.true.
     savefields=.false.
  end if

  nthread_max=1
  ithread=0
  nthread=1
  !control if we have the possibility of using OMP to de-synchronize communication
  !$ nthread_max=omp_get_max_threads()
  

  !decide the communication strategy
  unblock_comms_den=mpi_thread_funneled_is_supported .and. unblock_comms=='DEN' .and. &
      nthread_max > 1 .and. scf ! density is done only if scf is present
  if (unblock_comms_den) whilepot=.false. !anticipate the NlHamiltonian if density should be overlapped

  unblock_comms_pot=mpi_thread_funneled_is_supported .and. unblock_comms=='POT' .and. &
      nthread_max > 1 .and. whilepot

  !calculate the self-consistent potential
  if (scf) then
     !print *,'here',savefields,correcth,energs%ekin,energs%epot,dot(wfn%orbs%npsidim_orbs,wfn%psi(1),1,wfn%psi(1),1)
     ! Potential from electronic charge density 
     call sumrho(denspot%dpbox,wfn%orbs,wfn%Lzd,GPU,atoms%sym,denspot%rhod,wfn%psi,denspot%rho_psi)
     !print *,'here',wfn%orbs%occup(:),'there',savefields,correcth,energs%ekin,energs%epot,&
     !     dot(wfn%Lzd%Glr%d%n1i*wfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p*wfn%orbs%nspin,&
     !     denspot%rho_psi(1,1),1,denspot%rho_psi(1,1),1)
     !stop
     !initialize nested approach 
     !this has always to be done for using OMP parallelization in the 
     !projector case
     !if nesting is not supported, a bigdft_nesting routine should be defined
     !$   call OMP_SET_NESTED(.true.) 
     !$   call OMP_SET_MAX_ACTIVE_LEVELS(2)
     !$ if (unblock_comms_den) then
     !$   call OMP_SET_NUM_THREADS(2)
     !$ else
     !$   call OMP_SET_NUM_THREADS(1)
     !$ end if
     !print *,'how many threads ?' ,nthread_max
     !$OMP PARALLEL DEFAULT(shared), PRIVATE(ithread,nthread)
     !$ ithread=omp_get_thread_num()
     !$ nthread=omp_get_num_threads() !this should be 2 if active
     !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
     ! thread 0 does mpi communication 
     if (ithread == 0) then
       !communicate density 
        call communicate_density(denspot%dpbox,wfn%orbs%nspin,denspot%rhod,&
             denspot%rho_psi,denspot%rhov,.false.)
        !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
     end if
     if (ithread > 0 .or. nthread==1 .and. .not. whilepot) then
        ! Only the remaining threads do computations (if active) 
        !$ if (unblock_comms_den) then
        !$ call OMP_SET_NUM_THREADS(nthread_max-1)
        !$ else 
        !$ call OMP_SET_NUM_THREADS(nthread_max)
        !$ end if

        !nonlocal hamiltonian
        !$ if (verbose > 2 .and. iproc==0 .and. unblock_comms_den)&
        !$ & print *,'NonLocalHamiltonian with nthread:, out to:' ,omp_get_max_threads(),nthread_max
        if (wfn%orbs%npsidim_orbs > 0) call to_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
        call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs,rxyz,&
             proj,wfn%Lzd,nlpspd,wfn%psi,wfn%hpsi,energs%eproj)
     end if
     !$OMP END PARALLEL
     !finalize the communication scheme
     !$   call OMP_SET_NESTED(.false.) 
     !$   call OMP_SET_NUM_THREADS(nthread_max)
     ithread=0
     nthread=1

     !here the density can be mixed
     if (iscf > SCF_KIND_DIRECT_MINIMIZATION) then
        if (denspot%mix%kind == AB6_MIXING_DENSITY) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                atoms%alat1*atoms%alat2*atoms%alat3,&
                rpnrm,denspot%dpbox%nscatterarr)
           
           if (iproc == 0 .and. itrp > 1) then
              write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
                   &   'DENSITY iteration,Delta : (Norm 2/Volume)',itrp,rpnrm
              !yaml output
              !write(70,'(1x,a,1pe9.2,a,i5)')'DENSITY variation: &rpnrm',rpnrm,', #itrp: ',itrp
           end if
           ! xc_init_rho should be put in the mixing routines
           denspot%rhov = abs(denspot%rhov) + 1.0d-20
        end if
     end if
     call denspot_set_rhov_status(denspot, ELECTRONIC_DENSITY, itwfn, iproc, nproc)

     !before creating the potential, save the density in the second part 
     !in the case of NK SIC, so that the potential can be created afterwards
     !copy the density contiguously since the GGA is calculated inside the NK routines
     !with the savefield scheme, this can be avoided in the future
     if (wfn%SIC%approach=='NK') then !here the density should be copied somewhere else
        irhotot_add=wfn%Lzd%Glr%d%n1i*wfn%Lzd%Glr%d%n2i*denspot%dpbox%i3xcsh+1
        irho_add=wfn%Lzd%Glr%d%n1i*wfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d*wfn%orbs%nspin+1
        do ispin=1,wfn%orbs%nspin
           call dcopy(denspot%dpbox%ndimpot,&
                denspot%rhov(irhotot_add),1,denspot%rhov(irho_add),1)
           irhotot_add=irhotot_add+denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3d
           irho_add=irho_add+denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p
        end do
     end if

     if(wfn%orbs%nspinor==4) then
        !this wrapper can be inserted inside the XC_potential routine
        call PSolverNC(atoms%geocode,'D',denspot%pkernel%iproc,denspot%pkernel%nproc,&
             denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
             denspot%dpbox%n3d,ixc,&
             denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
             denspot%rhov,denspot%pkernel%kernel,denspot%V_ext,&
             energs%eh,energs%exc,energs%evxc,0.d0,.true.,4)
     else
        call XC_potential(atoms%geocode,'D',denspot%pkernel%iproc,denspot%pkernel%nproc,&
             denspot%pkernel%mpi_comm,&
             denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),ixc,&
             denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
             denspot%rhov,energs%exc,energs%evxc,wfn%orbs%nspin,denspot%rho_C,denspot%V_XC,xcstr)
        call denspot_set_rhov_status(denspot, CHARGE_DENSITY, itwfn, iproc, nproc)
        call H_potential('D',denspot%pkernel,&
             denspot%rhov,denspot%V_ext,energs%eh,0.0_dp,.true.,&
             quiet=denspot%PSquiet) !optional argument
        !this is not true, there is also Vext
        call denspot_set_rhov_status(denspot, HARTREE_POTENTIAL, itwfn, iproc, nproc)

        !sum the two potentials in rhopot array
        !fill the other part, for spin, polarised
        if (wfn%orbs%nspin == 2) then
           call dcopy(denspot%dpbox%ndimpot,denspot%rhov(1),1,&
                denspot%rhov(1+denspot%dpbox%ndimpot),1)
        end if
        !spin up and down together with the XC part
        call axpy(denspot%dpbox%ndimpot*wfn%orbs%nspin,&
             1.0_dp,denspot%V_XC(1,1,1,1),1,&
             denspot%rhov(1),1)

        !here a external potential with spinorial indices can be added
     end if

     !here the potential can be mixed
     if (iscf > SCF_KIND_DIRECT_MINIMIZATION) then
        if (denspot%mix%kind == AB6_MIXING_POTENTIAL) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                atoms%alat1*atoms%alat2*atoms%alat3,&!volume should be used 
                rpnrm,denspot%dpbox%nscatterarr)
           if (iproc == 0 .and. itrp > 1) then
              write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
                   &   'POTENTIAL iteration,Delta P (Norm 2/Volume)',itrp,rpnrm
              !yaml output
              !write(70,'(1x,a,1pe9.2,a,i5)')'POTENTIAL variation: &rpnrm',rpnrm,', #itrp: ',itrp
           end if
        end if
     end if
     call denspot_set_rhov_status(denspot, KS_POTENTIAL, itwfn, iproc, nproc)

     if (savefields) then
        if (associated(denspot%rho_work)) then
           write(*,*)'ERROR: the reference potential should be empty to correct the hamiltonian!'
           stop
        end if

        allocate(denspot%rho_work(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim+ndebug),stat=i_stat)
        call dcopy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,denspot%rhov(1),1,&
             denspot%rho_work(1),1)
     end if

     !if the hamiltonian should have to be just updated, perform the difference between potentials
     if (correcth) then
        if (.not. associated(denspot%rho_work)) then
           write(*,*)'ERROR: need a reference potential to correct the hamiltonian!'
           stop
        end if
        call yaml_newline()
        call yaml_comment('Calculating potential delta')
        !subtract the previous potential from the new one
        call axpy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,&
             -1.0_dp,denspot%rho_work(1),1,denspot%rhov(1),1)

        !deallocation should be deplaced
        i_all=-product(shape(denspot%rho_work))*kind(denspot%rho_work)
        deallocate(denspot%rho_work,stat=i_stat)
        call memocc(i_stat,i_all,'denspot%rho_work',subname)
        nullify(denspot%rho_work)
     end if


  end if

  !non self-consistent case: rhov should be the total potential
  if (denspot%rhov_is/=KS_POTENTIAL) then
     stop 'psitohpsi: KS_potential not available' 
  end if

  !temporary, to be corrected with comms structure
  if (wfn%exctxpar == 'OP2P') energs%eexctX = UNINITIALIZED(1.0_gp)

  !initialize nested approach 
  !this has always to be done for using OMP parallelization in the 
  !projector case
  !if nesting is not supported, a bigdft_nesting routine should be defined
  !$   call OMP_SET_NESTED(.true.) 
  !$   call OMP_SET_MAX_ACTIVE_LEVELS(2)
  !$ if (unblock_comms_pot) then
  !$   call OMP_SET_NUM_THREADS(2)
  !$ else
  !$   call OMP_SET_NUM_THREADS(1)
  !$ end if
  !print *,'how many threads ?' ,nthread_max
  !$OMP PARALLEL DEFAULT(shared), PRIVATE(ithread,nthread)
  !$ ithread=omp_get_thread_num()
  !$ nthread=omp_get_num_threads() !this should be 2 if active
  !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
  ! thread 0 does mpi communication 
  if (ithread == 0) then
     call full_local_potential(iproc,nproc,wfn%orbs,wfn%Lzd,linflag,&
          denspot%dpbox,denspot%rhov,denspot%pot_work)
     !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
  end if
  if (ithread > 0 .or. nthread==1 .and. whilepot) then
     ! Only the remaining threads do computations (if active) 
     !$ if (unblock_comms_pot) then
     !$ call OMP_SET_NUM_THREADS(nthread_max-1)
     !$ else 
     !$ call OMP_SET_NUM_THREADS(nthread_max)
     !$ end if

     !nonlocal hamiltonian
     !$ if (verbose > 2 .and. iproc==0 .and. unblock_comms_pot)&
     !$ & print *,'NonLocalHamiltonian with nthread:, out to:' ,omp_get_max_threads(),nthread_max
     if (wfn%orbs%npsidim_orbs >0) call to_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
     call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs,rxyz,&
          proj,wfn%Lzd,nlpspd,wfn%psi,wfn%hpsi,energs%eproj)
  end if
  !$OMP END PARALLEL
  !finalize the communication scheme
  !$   call OMP_SET_NESTED(.false.) 
  !$   call OMP_SET_NUM_THREADS(nthread_max)
  ithread=0
  nthread=1
 
  !here exctxpar might be passed
  !choose to just add the potential if needed
  call LocalHamiltonianApplication(iproc,nproc,atoms,wfn%orbs,&
       wfn%Lzd,wfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,wfn%psi,wfn%hpsi,&
       energs,wfn%SIC,GPU,correcth,pkernel=denspot%pkernelseq)
  call SynchronizeHamiltonianApplication(nproc,wfn%orbs,wfn%Lzd,GPU,wfn%hpsi,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  ! Emit that hpsi are ready.
  if (wfn%c_obj /= 0) then
     call kswfn_emit_psi(wfn, itwfn, 1, iproc, nproc)
  end if

  !deallocate potential
  call free_full_potential(denspot%dpbox%nproc,linflag,denspot%pot_work,subname)
  !----
  if (iproc==0 .and. verbose > 0) then
     if (correcth) then
        call yaml_map('Hamiltonian Applied',.false.)
        call yaml_comment('Only local potential')
     else
        call yaml_map('Hamiltonian Applied',.true.)
     end if
  end if

end subroutine psitohpsi

subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,rxyz,&
     proj,Lzd,nlpspd,confdatarr,ngatherarr,pot,psi,hpsi,&
     energs,SIC,GPU,&
     pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use module_interfaces, fake_name => FullHamiltonianApplication
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: Lzd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(SIC_data), intent(in) :: SIC
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  real(wp), dimension(lzd%ndimpotisf) :: pot
  type(energy_terms), intent(inout) :: energs
  !real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,evsic
  real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  type(coulomb_operator), intent(in), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc

  !put to zero hpsi array (now important since any of the pieces of the hamiltonian is accumulating)
  if (orbs%npsidim_orbs > 0) call to_zero(orbs%npsidim_orbs,hpsi(1))

  !write(*,*) 'lzd%ndimpotisf', lzd%ndimpotisf
  !do i=1,lzd%ndimpotisf
  !    write(210,*) pot(i)
  !end do

 if (.not. present(pkernel)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
         energs,SIC,GPU,.false.)
 else if (present(pkernel) .and. .not. present(orbsocc)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,.false.,pkernel=pkernel)
 else if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
         energs,SIC,GPU,.false.,pkernel,orbsocc,psirocc)
 else
    stop 'HamiltonianApplication, argument error'
 end if

  call NonLocalHamiltonianApplication(iproc,at,orbs,rxyz,&
       proj,Lzd,nlpspd,psi,hpsi,energs%eproj)

  call SynchronizeHamiltonianApplication(nproc,orbs,Lzd,GPU,hpsi,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)

END SUBROUTINE FullHamiltonianApplication


!> Application of the Local Hamiltonian
subroutine LocalHamiltonianApplication(iproc,nproc,at,orbs,&
     Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
     energs,SIC,GPU,onlypot,pkernel,orbsocc,psirocc)
   use module_base
   use module_types
   use module_xc
   use module_interfaces, except_this_one => LocalHamiltonianApplication
   use yaml_output
   implicit none
   logical, intent(in) :: onlypot !< if true, only the potential operator is applied
   integer, intent(in) :: iproc,nproc
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd 
   type(SIC_data), intent(in) :: SIC
   integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
   real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
   type(confpot_data), dimension(orbs%norbp) :: confdatarr
   !real(wp), dimension(:), pointer :: pot
   real(wp), dimension(*) :: pot
   type(energy_terms), intent(inout) :: energs
   real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout) :: hpsi
   type(GPU_pointers), intent(inout) :: GPU
   type(coulomb_operator), intent(in), optional :: pkernel
   type(orbitals_data), intent(in), optional :: orbsocc
   real(wp), dimension(:), pointer, optional :: psirocc
   !local variables
   character(len=*), parameter :: subname='HamiltonianApplication'
   logical :: exctX,op2p
   integer :: i_stat,n3p,ispot,ipotmethod
   real(gp) :: evsic_tmp
   type(coulomb_operator) :: pkernelSIC
   

   ! local potential and kinetic energy for all orbitals belonging to iproc
   if (iproc==0 .and. verbose > 1) then
      !call yaml_comment('Hamiltonian application, ',advance='no')
      !write(*,'(1x,a)',advance='no')&
      !     'Hamiltonian application...'
   end if

   !initialise exact exchange energy 
   op2p=(energs%eexctX == UNINITIALIZED(1.0_gp))
   energs%eexctX=0.0_gp
   energs%evsic=0.0_gp
   evsic_tmp=0.0_gp

   exctX = xc_exctXfac() /= 0.0_gp

   ispot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1

   !potential method
   !traditional case
   ipotmethod=0
   if (exctX) ipotmethod=1

   !the PZ-SIC correction does not makes sense for virtual orbitals procedure
   !if alphaSIC is zero no SIC correction
   if (SIC%approach == 'PZ' .and. .not. present(orbsocc) .and. SIC%alpha /= 0.0_gp ) ipotmethod=2
   if (SIC%approach == 'NK' .and. SIC%alpha /= 0.0_gp) ipotmethod=3

   !the poisson kernel should be present and associated in the case of SIC
   if ((ipotmethod /= 0) .and. present(pkernel)) then
      if (.not. associated(pkernel%kernel)) then
         if (iproc ==0) write(*,*)&
            &   'ERROR(LocalHamiltonianApplication): Poisson Kernel must be associated in SIC case'
         stop
      end if
   end if

   !associate the poisson kernel pointer in case of SIC
   if (ipotmethod == 2 .or. ipotmethod == 3) then
      pkernelSIC = pkernel
   else
      nullify(pkernelSIC%kernel)
   end if

   !fill the rest of the potential with the exact-exchange terms
   if (ipotmethod==1) then
      n3p=ngatherarr(iproc,1)/(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i)
      !exact exchange for virtual orbitals (needs psirocc)

      !here we have to add the round part
      if (present(psirocc) .and. present(orbsocc)) then
         call exact_exchange_potential_virt(iproc,nproc,at%geocode,orbs%nspin,&
              Lzd%Glr,orbsocc,orbs,ngatherarr(0,1),n3p,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernel,psirocc,psi,pot(ispot))
         energs%eexctX = 0._gp
      else
         !here the condition for the scheme should be chosen
         if (.not. op2p) then
            call exact_exchange_potential(iproc,nproc,at%geocode,orbs%nspin,&
                 Lzd%Glr,orbs,ngatherarr(0,1),n3p,&
                 0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
                 pkernel,psi,pot(ispot),energs%eexctX)
         else
            !the psi should be transformed in real space
            call exact_exchange_potential_round(iproc,nproc,at%geocode,orbs%nspin,Lzd%Glr,orbs,&
                 0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
                 pkernel,psi,pot(ispot),energs%eexctX)

         end if
      end if
      !print *,'iproc,eexctX',iproc,eexctX
   else if (ipotmethod==3) then
      !put fref=1/2 for the moment
      if (present(orbsocc) .and. present(psirocc)) then
         call NK_SIC_potential(Lzd%Glr,orbs,SIC%ixc,SIC%fref,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernelSIC,psi,pot(ispot),evsic_tmp,&
              potandrho=psirocc)
      else
         call NK_SIC_potential(Lzd%Glr,orbs,SIC%ixc,SIC%fref,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernelSIC,psi,pot(ispot),evsic_tmp)
      end if
   end if

   !GPU are supported only for ipotmethod=0
   if ((GPUconv .or. OCLconv) .and. ipotmethod /=0) then
      if (iproc ==0) write(*,*)&
         &   'ERROR(HamiltonianApplication): Accelerated hamiltonian are possible only with ipotmethod==0)'
      stop
   end if

   call timing(iproc,'ApplyLocPotKin','ON') 

   !apply the local hamiltonian for each of the orbitals
   !given to each processor
   !pot=0.d0
   !psi=1.d0
   !switch between GPU/CPU treatment
   !  do i=1,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
   !       call random_number(psi(i))
   !  end do
   if(OCLconv .or. GPUconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
      allocate(GPU%hpsi_ASYNC(max(1,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp)),stat=i_stat)
      call memocc(i_stat,GPU%hpsi_ASYNC,'GPU%hpsi_ASYNC',subname)
!      call to_zero((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp,GPU%hpsi_ASYNC(1))!hpsi(1))
   !else if (OCLconv) then
   !   GPU%hpsi_ASYNC => hpsi
   end if
   if (GPUconv) then
      call local_hamiltonian_GPU(orbs,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
           orbs%nspin,pot,psi,GPU%hpsi_ASYNC,energs%ekin,energs%epot,GPU)
   else if (OCLconv) then
      call local_hamiltonian_OCL(orbs,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
           orbs%nspin,pot,psi,GPU%hpsi_ASYNC,energs%ekin,energs%epot,GPU)
   else

!!$      !temporary allocation
!!$      allocate(fake_pot(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+ndebug),stat=i_stat)
!!$      call memocc(i_stat,fake_pot,'fake_pot',subname)
!!$
!!$      call to_zero(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin,fake_pot(1))

      !local hamiltonian application for different methods
      !print *,'here',ipotmethod,associated(pkernelSIC)
      if (.not. onlypot) then
         call local_hamiltonian(iproc,orbs,Lzd,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
              ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
              SIC%ixc,SIC%alpha,energs%ekin,energs%epot,energs%evsic)
!!$      i_all=-product(shape(fake_pot))*kind(fake_pot)
!!$      deallocate(fake_pot,stat=i_stat)
!!$      call memocc(i_stat,i_all,'fake_pot',subname)
         
      else

         call psi_to_vlocpsi(iproc,orbs,Lzd,&
              ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
              SIC%ixc,SIC%alpha,energs%epot,energs%evsic)
      end if

      !sum the external and the BS double counting terms
      energs%evsic=energs%evsic-SIC%alpha*evsic_tmp
   end if

   if (ipotmethod == 2 .or. ipotmethod==3) then
      nullify(pkernelSIC%kernel)
   end if
   call timing(iproc,'ApplyLocPotKin','OF') 

END SUBROUTINE LocalHamiltonianApplication


!> Routine which calculates the application of nonlocal projectors on the wavefunctions
!! Reduce the wavefunction in case it is needed
subroutine NonLocalHamiltonianApplication(iproc,at,orbs,rxyz,&
     proj,Lzd,nlpspd,psi,hpsi,eproj_sum)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc
   type(atoms_data), intent(in) :: at
   type(orbitals_data),  intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(nonlocal_psp_descriptors), intent(in) :: nlpspd
   real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
   real(gp), dimension(3,at%nat), intent(in) :: rxyz
   real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
   real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
   real(gp), intent(out) :: eproj_sum
   !local variables
   logical :: dosome, overlap
   integer :: ikpt,istart_ck,ispsi_k,isorb,ieorb,nspinor,iorb,iat,nwarnings
   integer :: iproj,ispsi,istart_c,ilr,ilr_skip,mproj

   eproj_sum=0.0_gp

   !quick return if no orbitals on this processor
   if (orbs%norbp == 0) then
      return
   end if

   ! apply all PSP projectors for all orbitals belonging to iproc
   call timing(iproc,'ApplyProj     ','ON')

   nwarnings=0

   !here the localisation region should be changed, temporary only for cubic approach

   !apply the projectors following the strategy (On-the-fly calculation or not)

   !apply the projectors  k-point of the processor
   !starting k-point
   ikpt=orbs%iokpt(1)
   istart_ck=1
   ispsi_k=1
   loop_kpt: do

      call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
      !localisation regions loop
      loop_lr: do ilr=1,Lzd%nlr
         !do something only if at least one of the orbitals lives in the ilr
         dosome=.false.
         do iorb=isorb,ieorb
            !dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
            dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr .and. lzd%doHamAppl(ilr))
            if (dosome) exit
         end do
         if (.not. dosome) cycle loop_lr

         if (DistProjApply) then
            !first create a projector ,then apply it for everyone
            iproj=0
            do iat=1,at%nat
               ! Check if atom has projectors, if not cycle
               call numb_proj(at%iatype(iat),at%ntypes,at%psppar,at%npspcode,mproj) 
               if(mproj == 0) cycle

               !check if the atom projector intersect with the given localisation region
               call check_overlap(Lzd%Llr(ilr), nlpspd%plr(iat), Lzd%Glr, overlap)
               if(.not. overlap) cycle

               ! Now create the projector
               istart_c=1
               call atom_projector(ikpt,iat,0,istart_c,iproj,&
                    nlpspd%nprojel,&
                    Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),rxyz(1,iat),at,orbs,&
                    nlpspd%plr(iat),proj,nwarnings)

               !apply the projector to all the orbitals belonging to the processor
               ispsi=ispsi_k
               do iorb=isorb,ieorb
                  if (orbs%inwhichlocreg(iorb+orbs%isorb) /= ilr) then
                     !increase ispsi to meet orbital index
                     ilr_skip=orbs%inwhichlocreg(iorb+orbs%isorb)
                     ispsi=ispsi+(Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                     cycle
                  end if
                  istart_c=1
                  call apply_atproj_iorb_new(iat,iorb,istart_c,&
                       nlpspd%nprojel,&
                       at,orbs,Lzd%Llr(ilr)%wfd,nlpspd%plr(iat),&
                       proj,psi(ispsi),hpsi(ispsi),eproj_sum)
!                print *,'iorb,iat,eproj',iorb+orbs%isorb,ispsi,iat,eproj_sum
                  ispsi=ispsi+&
                       (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
               end do

            end do

            !if (iproj /= nlpspd%nproj) stop &
            !     'NonLocal HamiltonianApplication: incorrect number of projectors created'     
            !for the moment, localization region method is not implemented with
            !once-and-for-all calculation
         else if (Lzd%nlr == 1) then

            !loop over the interesting my orbitals, and apply all the projectors over all orbitals
            ispsi=ispsi_k
            do iorb=isorb,ieorb
               if (orbs%inwhichlocreg(iorb+orbs%isorb) /= ilr) then
                  !increase ispsi to meet orbital index
                  ilr_skip=orbs%inwhichlocreg(iorb+orbs%isorb)
                  ispsi=ispsi+(Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                  cycle
               end if

               istart_c=istart_ck !TO BE CHANGED IN ONCE-AND-FOR-ALL 
               do iat=1,at%nat
                  ! Check if atom has projectors, if not cycle
                  call numb_proj(at%iatype(iat),at%ntypes,at%psppar,at%npspcode,mproj) 
                  if(mproj == 0) cycle
                  !check if the atom intersect with the given localisation region
                  call check_overlap(Lzd%Llr(ilr), nlpspd%plr(iat), Lzd%Glr, overlap)
                  if(.not. overlap) stop 'ERROR all atoms should be in global'
                  call apply_atproj_iorb_new(iat,iorb,istart_c,nlpspd%nprojel,&
                       at,orbs,Lzd%Llr(ilr)%wfd,nlpspd%plr(iat),&
                       proj,psi(ispsi),hpsi(ispsi),eproj_sum)
                  !print *,'iorb,iat,eproj',iorb+orbs%isorb,iat,eproj_sum
               end do
               ispsi=ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
            end do
            istart_ck=istart_c !TO BE CHANGED IN THIS ONCE-AND-FOR-ALL
         else
           ! COULD CHANGE THIS NOW !!
           stop 'Localization Regions not allowed in once-and-for-all'    
         end if

      end do loop_lr

      !last k-point has been treated
      if (ieorb == orbs%norbp) exit loop_kpt
      
      ikpt=ikpt+1
      ispsi_k=ispsi

   end do loop_kpt

   if (.not. DistProjApply) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
      if (istart_ck-1 /= nlpspd%nprojel) &
         &   stop 'incorrect once-and-for-all psp application'
   end if
   !for the moment it has to be removed. A number of components in orbital distribution should be defined
   !if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'

   !used on the on-the-fly projector creation
   if (nwarnings /= 0 .and. iproc == 0) then
      write(*,'(1x,a,i0,a)')'found ',nwarnings,' warnings.'
      write(*,'(1x,a)')'Some projectors may be too rough.'
      write(*,'(1x,a,f6.3)')&
         &   'Consider the possibility of reducing hgrid for having a more accurate run.'
   end if


   call timing(iproc,'ApplyProj     ','OF')

END SUBROUTINE NonLocalHamiltonianApplication

!> routine which puts a barrier to ensure that both local and nonlocal hamiltonians have been applied
!! in the GPU case puts a barrier to end the overlapped Local and nonlocal applications
subroutine SynchronizeHamiltonianApplication(nproc,orbs,Lzd,GPU,hpsi,ekin_sum,epot_sum,eproj_sum,evsic,eexctX)
   use module_base
   use module_types
   use module_xc
   implicit none
   integer, intent(in) :: nproc
   type(orbitals_data),  intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(GPU_pointers), intent(inout) :: GPU
   real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,evsic,eexctX
   real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
   !local variables
   character(len=*), parameter :: subname='SynchronizeHamiltonianApplication'
   logical :: exctX
   integer :: i_all,i_stat,ierr,iorb,ispsi,ilr
   real(gp), dimension(4) :: wrkallred

   if(OCLconv .or. GPUconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
      if (OCLconv) call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU)
      ispsi=1
      do iorb=1,orbs%norbp
         ilr=orbs%inWhichLocreg(orbs%isorb+iorb)
         call axpy((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor,&
              1.0_wp,GPU%hpsi_ASYNC(ispsi),1,hpsi(ispsi),1)
         ispsi=ispsi+&
             (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
      end do
      i_all=-product(shape(GPU%hpsi_ASYNC))*kind(GPU%hpsi_ASYNC)
      deallocate(GPU%hpsi_ASYNC,stat=i_stat)
      call memocc(i_stat,i_all,'GPU%hpsi_ASYNC',subname)
   endif

   exctX = xc_exctXfac() /= 0.0_gp

   !energies reduction
   if (nproc > 1) then
      wrkallred(1)=ekin_sum 
      wrkallred(2)=epot_sum 
      wrkallred(3)=eproj_sum
      wrkallred(4)=evsic

      call mpiallred(wrkallred(1),4,MPI_SUM,MPI_COMM_WORLD,ierr)

      ekin_sum=wrkallred(1)
      epot_sum=wrkallred(2)
      eproj_sum=wrkallred(3) 
      evsic=wrkallred(4) 
   endif

   !up to this point, the value of the potential energy is 
   !only taking into account the local potential part
   !whereas it should consider also the value coming from the 
   !exact exchange operator (twice the exact exchange energy)
   !this operation should be done only here since the exctX energy is already reduced
   if (exctX) epot_sum=epot_sum+2.0_gp*eexctX

END SUBROUTINE SynchronizeHamiltonianApplication


!> Build the potential in the whole box
!! Control also the generation of an orbital
!! @ param i3rho_add Integer which controls the presence of a density after the potential array
!!                   if different than zero, at the address ndimpot*nspin+i3rho_add starts the spin up component of the density
!!                   the spin down component can be found at the ndimpot*nspin+i3rho_add+ndimpot, contiguously
!!                   the same holds for non-collinear calculations
subroutine full_local_potential(iproc,nproc,orbs,Lzd,iflag,dpbox,potential,pot,comgp)
  !ndimpot,ndimgrid,nspin,&
  !   ndimrhopot,i3rho_add,orbs,&
  !   Lzd,iflag,ngatherarr,potential,pot,comgp)
   use module_base
   use module_types
   use module_xc
   implicit none
   integer, intent(in) :: iproc,nproc,iflag!,nspin,ndimpot,ndimgrid
   !integer, intent(in) :: ndimrhopot,i3rho_add
   type(orbitals_data),intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: Lzd
   type(denspot_distribution), intent(in) :: dpbox
   !integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
   real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), target :: potential !< Distributed potential. Might contain the density for the SIC treatments
   real(wp), dimension(:), pointer :: pot
   !type(p2pCommsGatherPot),intent(inout), optional:: comgp
   type(p2pComms),intent(inout), optional:: comgp
   !local variables
   character(len=*), parameter :: subname='full_local_potential'
   logical :: odp,newvalue !orbital dependent potential
   integer :: npot,ispot,ispotential,ispin,ierr,i_stat,i_all,ii,ilr,iorb,iorb2,nilr
   integer:: istl, ist, size_Lpot, i3s, i3e
   integer,dimension(:),allocatable:: ilrtable
   real(wp), dimension(:), pointer :: pot1
   
   call timing(iproc,'Pot_commun    ','ON')

   odp = (xc_exctXfac() /= 0.0_gp .or. (dpbox%i3rho_add /= 0 .and. orbs%norbp > 0))

   !############################################################################
   ! Build the potential on the whole simulation box
   ! NOTE: in the linear scaling case this should be done for a given localisation
   !       region this routine should then be modified or integrated in HamiltonianApplication
   ! WARNING : orbs%nspin and nspin are not the same !! Check if orbs%nspin should be replaced everywhere
   !#############################################################################
   if (iflag<2) then

      !determine the dimension of the potential array
      if (odp) then
         if (xc_exctXfac() /= 0.0_gp) then
            npot=dpbox%ndimgrid*orbs%nspin+&
                 &   max(max(dpbox%ndimgrid*orbs%norbp,dpbox%ngatherarr(0,1)*orbs%norb),1) !part which refers to exact exchange
         else if (dpbox%i3rho_add /= 0 .and. orbs%norbp > 0) then
            npot=dpbox%ndimgrid*orbs%nspin+&
                 &   dpbox%ndimgrid*max(orbs%norbp,orbs%nspin) !part which refers to SIC correction
         end if
      else
         npot=dpbox%ndimgrid*orbs%nspin
      end if
!      write(*,*) 'dpbox%ndimgrid, orbs%norbp, npot, odp', dpbox%ndimgrid, orbs%norbp, npot, odp
!      write(*,*)'nspin',orbs%nspin,dpbox%i3rho_add,dpbox%ndimpot,dpbox%ndimrhopot,sum(potential)
!      write(*,*)'iproc',iproc,'ngatherarr',dpbox%ngatherarr(:,1),dpbox%ngatherarr(:,2)

      !build the potential on the whole simulation box
      !in the linear scaling case this should be done for a given localisation region
      !this routine should then be modified or integrated in HamiltonianApplication
      if (dpbox%nproc > 1) then

         allocate(pot1(npot+ndebug),stat=i_stat)
         call memocc(i_stat,pot1,'pot1',subname)
         ispot=1
         ispotential=1
         do ispin=1,orbs%nspin
            call MPI_ALLGATHERV(potential(ispotential),dpbox%ndimpot,&
                 &   mpidtypw,pot1(ispot),dpbox%ngatherarr(0,1),&
                 dpbox%ngatherarr(0,2),mpidtypw,dpbox%mpi_comm,ierr)
            ispot=ispot+dpbox%ndimgrid
            ispotential=ispotential+max(1,dpbox%ndimpot)
         end do
         !continue to copy the density after the potential if required
         if (dpbox%i3rho_add >0 .and. orbs%norbp > 0) then
            ispot=ispot+dpbox%i3rho_add-1
            do ispin=1,orbs%nspin
               call MPI_ALLGATHERV(potential(ispotential),dpbox%ndimpot,&
                    &   mpidtypw,pot1(ispot),dpbox%ngatherarr(0,1),&
                    dpbox%ngatherarr(0,2),mpidtypw,dpbox%mpi_comm,ierr)
               ispot=ispot+dpbox%ndimgrid
               ispotential=ispotential+max(1,dpbox%ndimpot)
            end do
         end if
      else
         if (odp) then
            allocate(pot1(npot+ndebug),stat=i_stat)
            call memocc(i_stat,pot1,'pot1',subname)
            call dcopy(dpbox%ndimgrid*orbs%nspin,potential,1,pot1,1)
            if (dpbox%i3rho_add >0 .and. orbs%norbp > 0) then
               ispot=dpbox%ndimgrid*orbs%nspin+1
               call dcopy(dpbox%ndimgrid*orbs%nspin,potential(ispot+dpbox%i3rho_add),1,pot1(ispot),1)
            end if
         else
            pot1 => potential
         end if
      end if
   else
       !!if(.not.comgp%communication_complete) call gatherPotential(iproc, nproc, comgp)
       !!if(.not.comgp%communication_complete) call wait_p2p_communication(iproc, nproc, comgp)
       call wait_p2p_communication(iproc, nproc, comgp)
   end if


   !########################################################################
   ! Determine the dimension of the potential array and orbs%ispot
   !########################################################################
!!$   if(associated(orbs%ispot)) then
!!$      nullify(orbs%ispot)
!!$      !     i_all=-product(shape(orbs%ispot))*kind(orbs%ispot)
!!$      !     deallocate(orbs%ispot,stat=i_stat)
!!$      !     call memocc(i_stat,i_all,'orbs%ispot',subname)
!!$   end if
!!$   allocate(orbs%ispot(orbs%norbp),stat=i_stat)
!!$   call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)

   if(Lzd%nlr > 1) then
      allocate(ilrtable(orbs%norbp),stat=i_stat)
      call memocc(i_stat,ilrtable,'ilrtable',subname)
      !call to_zero(orbs%norbp*2,ilrtable(1,1))
      ilrtable=0
      ii=0
      do iorb=1,orbs%norbp
         newvalue=.true.
         !localization region to which the orbital belongs
         ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
         !spin state of the orbital
         if (orbs%spinsgn(orbs%isorb+iorb) > 0.0_gp) then
            ispin = 1       
         else
            ispin=2
         end if
         !check if the orbitals already visited have the same conditions
         loop_iorb2: do iorb2=1,orbs%norbp
            if(ilrtable(iorb2) == ilr) then
               newvalue=.false.
               exit loop_iorb2
            end if
         end do loop_iorb2
         if (newvalue) then
            ii = ii + 1
            ilrtable(ii)=ilr
         end if
      end do
      !number of inequivalent potential regions
      nilr = ii
   else 
      allocate(ilrtable(1),stat=i_stat)
      call memocc(i_stat,ilrtable,'ilrtable',subname)
      nilr = 1
      ilrtable=1
   end if

!!$   !calculate the dimension of the potential in the gathered form 
!!$   !this part has been deplaced in check_linear_and_create_Lzd routine 
!!$   lzd%ndimpotisf=0
!!$   do iilr=1,nilr
!!$      ilr=ilrtable(iilr,1)
!!$      do iorb=1,orbs%norbp
!!$         !put the starting point
!!$         if (orbs%inWhichLocreg(iorb+orbs%isorb) == ilr) then
!!$            !assignment of ispot array to the value of the starting address of inequivalent
!!$            orbs%ispot(iorb)=lzd%ndimpotisf + 1
!!$            if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
!!$               orbs%ispot(iorb)=lzd%ndimpotisf + &
!!$                    1 + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
!!$            end if
!!$         end if
!!$      end do
!!$      lzd%ndimpotisf = lzd%ndimpotisf + &
!!$           lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*nspin
!!$   end do
!!$   !part which refers to exact exchange
!!$   if (exctX) then
!!$      lzd%ndimpotisf = lzd%ndimpotisf + &
!!$           max(max(ndimgrid*orbs%norbp,ngatherarr(0,1)*orbs%norb),1) 
!!$   end if

   !#################################################################################################################################################
   ! Depending on the scheme, cut out the local pieces of the potential
   !#################################################################################################################################################
   if(iflag==0) then
      !       allocate(pot(lzd%ndimpotisf+ndebug),stat=i_stat)
      !       call dcopy(lzd%ndimpotisf,pot,1,pot,1) 
      pot=>pot1
      !print *,iproc,shape(pot),shape(pot1),'shapes'
      !print *,'potential sum',iproc,sum(pot)
   else if(iflag<2 .and. iflag>0) then

      allocate(pot(lzd%ndimpotisf+ndebug),stat=i_stat)
      call memocc(i_stat,pot,'pot',subname)
      ! Cut potential
      istl=1
      do iorb=1,nilr
         ilr = ilrtable(iorb)
         ! Cut the potential into locreg pieces
         call global_to_local(Lzd%Glr,Lzd%Llr(ilr),orbs%nspin,npot,lzd%ndimpotisf,pot1,pot(istl))
         istl = istl + Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspin
      end do
   else
      allocate(pot(lzd%ndimpotisf+ndebug),stat=i_stat)
      call memocc(i_stat,pot,'pot',subname)

      ist=1
      do iorb=1,nilr
         ilr = ilrtable(iorb)
         !determine the dimension of the potential array (copied from full_local_potential)
         if (xc_exctXfac() /= 0.0_gp) then
            stop 'exctX not yet implemented!'
         else
            size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
         end if

         ! Extract the part of the potential which is needed for the current localization region.
         i3s=lzd%Llr(ilr)%nsi3-comgp%ise3(1,iproc)+2 ! starting index of localized  potential with respect to total potential in comgp%recvBuf
         i3e=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i-comgp%ise3(1,iproc)+1 ! ending index of localized potential with respect to total potential in comgp%recvBuf
         if(i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i) then
            write(*,'(a,i0,3x,i0)') 'ERROR: i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i',i3e-i3s+1, Lzd%Llr(ilr)%d%n3i
            stop
         end if

         call global_to_local_parallel(lzd%Glr, lzd%Llr(ilr), orbs%nspin, comgp%nrecvBuf, size_Lpot,&
              comgp%recvBuf, pot(ist), i3s, i3e)

         ist = ist + size_lpot
      end do
   end if

   i_all=-product(shape(ilrtable))*kind(ilrtable)
   deallocate(ilrtable,stat=i_stat)
   call memocc(i_stat,i_all,'ilrtable',subname)

   ! Deallocate pot.
   if (iflag<2 .and. iflag>0) then
      if (dpbox%nproc > 1) then
         i_all=-product(shape(pot1))*kind(pot1)
         deallocate(pot1,stat=i_stat)
         call memocc(i_stat,i_all,'pot1',subname)
      else
         if (xc_exctXfac() /= 0.0_gp) then
            i_all=-product(shape(pot1))*kind(pot1)
            deallocate(pot1,stat=i_stat)
            call memocc(i_stat,i_all,'pot1',subname)
         else
            nullify(pot1)
         end if
      end if
   end if

   call timing(iproc,'Pot_commun    ','OF') 

END SUBROUTINE full_local_potential


subroutine free_full_potential(nproc,flag,pot,subname)
   use module_base
   use module_xc
   implicit none
   character(len=*), intent(in) :: subname
   integer, intent(in) :: nproc, flag
   real(wp), dimension(:), pointer :: pot
   !local variables
   logical :: odp
   integer :: i_all,i_stat

   odp = xc_exctXfac() /= 0.0_gp
   if (nproc > 1 .or. odp .or. flag > 0 ) then
      i_all=-product(shape(pot))*kind(pot)
      deallocate(pot,stat=i_stat)
      call memocc(i_stat,i_all,'pot',subname)
      nullify(pot)
   else
      nullify(pot)
   end if

END SUBROUTINE free_full_potential

!> Calculate total energies from the energy terms
subroutine total_energies(energs, iter, iproc)
  use module_base
  use module_types
  implicit none
  type(energy_terms), intent(inout) :: energs
  integer, intent(in) :: iter, iproc

  !band structure energy calculated with occupation numbers
  energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
  !this is the Kohn-Sham energy
  energs%eKS=energs%ebs-energs%eh+energs%exc-energs%evxc-&
       energs%eexctX-energs%evsic+energs%eion+energs%edisp!-energs%excrhoc

  if (energs%c_obj /= 0) then
     call timing(iproc,'energs_signals','ON')
     call energs_emit(energs%c_obj, iter, 0) ! 0 is for BIGDFT_E_KS in C.
     call timing(iproc,'energs_signals','OF')
  end if
end subroutine total_energies

!> Extract the energy (the quantity which has to be minimised by the wavefunction)
!! and calculate the corresponding gradient.
!! The energy can be the actual Kohn-Sham energy or the trace of the hamiltonian, 
!! depending of the functional we want to calculate. The gradient wrt the wavefunction
!! is put in hpsi accordingly to the functional
subroutine calculate_energy_and_gradient(iter,iproc,nproc,GPU,ncong,iscf,&
     energs,wfn,gnrm,gnrm_zero)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,iscf,iter
  type(energy_terms), intent(inout) :: energs
  type(GPU_pointers), intent(in) :: GPU
  type(DFT_wavefunction), intent(inout) :: wfn
  real(gp), intent(out) :: gnrm,gnrm_zero
  !local variables
  character(len=*), parameter :: subname='calculate_energy_and_gradient' 
  logical :: lcs
  integer :: ierr,ikpt,iorb,i_all,i_stat,k
  real(gp) :: rzeroorbs,tt,garray(2)
  real(wp), dimension(:,:,:), pointer :: mom_vec


!!$  !calculate the entropy contribution (TO BE VERIFIED for fractional occupation numbers and Fermi-Dirac Smearing)
!!$  eTS=0.0_gp
!!$  do iorb=1,orbs%norbu  ! for closed shell case
!!$     !  if (iproc == 0)  print '("iorb,occup,eval,fermi:  ",i,e10.2,e27.17,e27.17)',iorb,orbs%occup(iorb),orbs%eval(iorb),orbs%efermi
!!$     eTS=eTS+exp(-((orbs%eval(iorb)-orbs%efermi)/in%Tel)**2)
!!$  enddo
!!$  if eTS=eTS*2._gp   ! for closed shell case
!!$  eTS=in%Tel/(2._gp*sqrt(3.1415926535897932_gp))* eTS
!!$  energy=energy-eTS
!!$  if (iproc == 0)  print '(" Free energy (energy-ST) = ",e27.17,"  , ST= ",e27.17," ,energy= " , e27.17)',energy,ST,energy+ST

!!$  !band structure energy calculated with occupation numbers
!!$  energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
!!$  !this is the Kohn-Sham energy
!!$  energs%eKS=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX-energs%evsic+energs%eion+energs%edisp

  !calculate orbital polarisation directions
  if(wfn%orbs%nspinor==4) then
     allocate(mom_vec(4,wfn%orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,wfn%orbs%norb,wfn%orbs%norb_par,&
          wfn%Lzd%Glr%wfd%nvctr_c+7*wfn%Lzd%Glr%wfd%nvctr_f,wfn%orbs%nspinor,wfn%psi,mom_vec)
  else
     nullify(mom_vec)
  end if


  if (iproc==0 .and. verbose > 1) then
!!$     write(*,'(1x,a)',advance='no')&
!!$          &   'done,  orthoconstraint...'
     !call yaml_comment('Orthoconstraint, ',advance='no')
  end if
  

  !transpose the hpsi wavefunction
   call transpose_v2(iproc,nproc,wfn%orbs,wfn%Lzd,wfn%comms,wfn%hpsi,work=wfn%psi)

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     wfn%psit => wfn%psi
     call transpose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%psit)
  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  !takes also into account parallel k-points distribution
  !here the orthogonality with respect to other occupied functions should be 
  !passed as an optional argument
  call orthoconstraint(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%hpsi,energs%trH) !n(m)

  !retranspose the hpsi wavefunction
   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%hpsi,work=wfn%psi)

  !after having calcutated the trace of the hamiltonian, the functional have to be defined
  !new value without the trace, to be added in hpsitopsi
  if (iscf >1) then
     wfn%diis%energy=energs%trH
  else
     wfn%diis%energy=energs%eKS!trH-eh+exc-evxc-eexctX+eion+edisp(not correct for non-integer occnums)
  end if

  if (iproc==0 .and. verbose > 0) call yaml_map('Orthoconstraint',.true.)

  !check that the trace of the hamiltonian is compatible with the 
  !band structure energy 
  !this can be done only if the occupation numbers are all equal
  tt=(energs%ebs-energs%trH)/energs%trH
!print *,'tt,energybs,trH',tt,energybs,trH
  if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
       &   (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
     !write this warning only if the system is closed shell
     call check_closed_shell(wfn%orbs,lcs)
     if (lcs) then
        call yaml_newline()
        call yaml_open_map('Energy inconsistencies')
           call yaml_map('Band Structure Energy',energs%ebs,fmt='(1pe22.14)')
           call yaml_map('Trace of the Hamiltonian',energs%trH,fmt='(1pe22.14)')
           call yaml_map('Relative inconsistency',tt,fmt='(1pe9.2)')
        call yaml_close_map()
!        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
!             &   'ERROR: inconsistency between gradient and energy',tt,energs%ebs,energs%trH
     end if
  endif


  call timing(iproc,'Precondition  ','ON')
  if (iproc==0 .and. verbose > 1) then
     !call yaml_comment('Preconditioning')
     !write(*,'(1x,a)',advance='no')&
     !     &   'done,  preconditioning...'
  end if

  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  !switch between CPU and GPU treatment
  if (GPUconv) then
     call preconditionall_GPU(wfn%orbs,wfn%Lzd%Glr,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),ncong,&
          wfn%hpsi,gnrm,gnrm_zero,GPU)
  else if (OCLconv) then
     call preconditionall_OCL(wfn%orbs,wfn%Lzd%Glr,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),ncong,&
          wfn%hpsi,gnrm,gnrm_zero,GPU)
  else
     !this is the final routine, the confining potential has to be passed to 
     !switch between the global and delocalized preconditioner
     call preconditionall2(iproc,nproc,wfn%orbs,wfn%Lzd,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),&
          ncong,wfn%hpsi,wfn%confdatarr,gnrm,gnrm_zero)
     if(.false.) then
        call preconditionall(wfn%orbs,wfn%Lzd%Glr,&
             wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),&
             ncong,wfn%hpsi,gnrm,gnrm_zero)
     end if
  end if

  call timing(iproc,'Precondition  ','OF')

  !sum over all the partial residues
  if (nproc > 1) then
      garray(1)=gnrm
      garray(2)=gnrm_zero
     call mpiallred(garray(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)
      gnrm     =garray(1)
      gnrm_zero=garray(2)
  endif

  !count the number of orbitals which have zero occupation number
  !weight this with the corresponding k point weight
  rzeroorbs=0.0_gp
  do ikpt=1,wfn%orbs%nkpts
     do iorb=1,wfn%orbs%norb
        if (wfn%orbs%occup(iorb+(ikpt-1)*wfn%orbs%norb) == 0.0_gp) then
           rzeroorbs=rzeroorbs+wfn%orbs%kwgts(ikpt)
        end if
     end do
  end do
  !commented out, the kwgts sum already to one
  !if (orbs%nkpts > 1) nzeroorbs=nint(real(nzeroorbs,gp)/real(orbs%nkpts,gp))

  gnrm=sqrt(gnrm/(real(wfn%orbs%norb,gp)-rzeroorbs))

  if (rzeroorbs /= 0.0_gp) then
     gnrm_zero=sqrt(gnrm_zero/rzeroorbs)
  else
     gnrm_zero=0.0_gp
  end if

  if (iproc==0 .and. verbose > 1) then
     !write(*,'(1x,a)')&
     !     &   'done.'
  end if

  if (iproc==0  .and. verbose > 0) then
     call yaml_map('Preconditioning',.true.)
     call yaml_newline()
  end if

  if (wfn%orbs%nspinor == 4) then
     !only the root process has the correct array
     if(iproc==0 .and. verbose > 0) then
        call yaml_open_sequence('Magnetic polarization per orbital')
        call yaml_newline()
        !write(*,'(1x,a)')&
        !     &   'Magnetic polarization per orbital'
!!$        write(*,'(1x,a)')&
!!$             &   '  iorb    m_x       m_y       m_z'
        do iorb=1,wfn%orbs%norb
           call yaml_sequence(advance='no')
           call yaml_open_map(flow=.true.)
           call yaml_map('iorb',iorb,fmt='(i5)')
           call yaml_map('M',(/(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)/),fmt='(3f10.5)')
           call yaml_close_map()
           if (iorb < wfn%orbs%norb)call yaml_newline()
           !write(*,'(1x,i5,3f10.5)') &
           !     &   iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
        end do
        call yaml_close_sequence()
     end if
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if


  !write the energy information
  if (iproc == 0) then
     call write_energies(iter,iscf,energs,gnrm,gnrm_zero,' ')
  endif

END SUBROUTINE calculate_energy_and_gradient

!> Operations after h|psi> 
!! (transposition, orthonormalisation, inverse transposition)
subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn)
   use module_base
   use module_types
   use module_interfaces, except_this_one_A => hpsitopsi
   use yaml_output
   implicit none
   integer, intent(in) :: iproc,nproc,idsx,iter
   type(DFT_wavefunction), intent(inout) :: wfn
   !local variables
   !n(c) character(len=*), parameter :: subname='hpsitopsi'

   !adjust the save variables for DIIS/SD switch
   if (iter == 1) then
      wfn%diis%ids=0
      wfn%diis%mids=1
      wfn%diis%idiistol=0
   end if
   !update variables at each iteration step
   if (idsx > 0) then
      wfn%diis%mids=mod(wfn%diis%ids,idsx)+1
      wfn%diis%ids=wfn%diis%ids+1
   end if

   wfn%diis%energy_min=min(wfn%diis%energy_min,wfn%diis%energy)

   !transpose the hpsi wavefunction
   call transpose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%hpsi,work=wfn%psi)
   
   !!experimental, orthogonalize the preconditioned gradient wrt wavefunction
   !call orthon_virt_occup(iproc,nproc,orbs,orbs,comms,comms,psit,hpsi,(verbose > 2))

   !apply the minimization method (DIIS or steepest descent)
   call timing(iproc,'Diis          ','ON')

   call psimix(iproc,nproc,sum(wfn%comms%ncntt(0:nproc-1)),wfn%orbs,wfn%comms,wfn%diis,wfn%hpsi,wfn%psit)

   call timing(iproc,'Diis          ','OF')

   if (iproc == 0 .and. verbose > 1) then
      !write(*,'(1x,a)',advance='no')&
      !&   'Orthogonalization...'
      call yaml_map('Orthogonalization Method',wfn%orthpar%methortho,fmt='(i3)')
   end if

   call orthogonalize(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%orthpar)

   !       call checkortho_p(iproc,nproc,norb,nvctrp,psit)

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit,work=wfn%hpsi,outadd=wfn%psi(1))
   if (nproc == 1) then
      nullify(wfn%psit)
   end if

   if (iproc == 0 .and. verbose > 1) then
      !write(*,'(1x,a)')&
      !   &   'done.'
   end if

   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   call diis_or_sd(iproc,idsx,wfn%orbs%nkptsp,wfn%diis)

   !previous value already filled
   wfn%diis%energy_old=wfn%diis%energy

END SUBROUTINE hpsitopsi


!> Choose among the wavefunctions a subset of them
!! Rebuild orbital descriptors for the new space and allocate the psi_as wavefunction
!! By hypothesis the work array is big enough to contain both wavefunctions
!! This routine has to be tested
subroutine select_active_space(iproc,nproc,orbs,comms,mask_array,Glr,orbs_as,comms_as,psi,psi_as)
   use module_base
   use module_types
   use module_interfaces, except_this_one => select_active_space
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(locreg_descriptors), intent(in) :: Glr
   type(communications_arrays), intent(in) :: comms
   logical, dimension(orbs%norb*orbs%nkpts), intent(in) :: mask_array
   real(wp), dimension(orbs%npsidim_comp), intent(in) :: psi
   type(orbitals_data), intent(out) :: orbs_as
   type(communications_arrays), intent(out) :: comms_as
   real(wp), dimension(:), pointer :: psi_as
   !local variables
   character(len=*), parameter :: subname='select_active_space'
   integer :: iorb,ikpt,norbu_as,norbd_as,icnt,ikptp,ispsi,ispsi_as
   integer :: i_stat,nvctrp

   !count the number of orbitals of the active space
   norbu_as=-1
   norbd_as=-1
   do ikpt=1,orbs%nkpts
      icnt=0
      do iorb=1,orbs%norbu
         if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
      end do
      if (norbu_as /= icnt .and. norbu_as /= -1) then
         write(*,*)'ERROR(select_active_space): the mask array should define always the same norbu'
         stop
      end if
      norbu_as=icnt
      icnt=0
      do iorb=orbs%norbu+1,orbs%norbu+orbs%norbd
         if (mask_array(iorb+(ikpt-1)*orbs%norb)) icnt=icnt+1
      end do
      if (norbd_as /= icnt .and. norbd_as /= -1) then
         write(*,*)'ERROR(select_active_space): the mask array should define always the same norbd'
         stop
      end if
      norbd_as=icnt
   end do

   !allocate the descriptors of the active space
   call orbitals_descriptors(iproc,nproc,norbu_as+norbd_as,norbu_as,norbd_as, &
        orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbs_as,&
        .false.,basedist=orbs%norb_par(0:,1))
   !allocate communications arrays for virtual orbitals
   call orbitals_communicators(iproc,nproc,Glr,orbs_as,comms_as,basedist=comms_as%nvctr_par(0:,1))  
   !allocate array of the eigenvalues
   allocate(orbs_as%eval(orbs_as%norb*orbs_as%nkpts+ndebug),stat=i_stat)
   call memocc(i_stat,orbs_as%eval,'orbs_as%eval',subname)

   !fill the orbitals array with the values and the wavefunction in transposed form
   icnt=0
   do iorb=1,orbs%nkpts*orbs%norb
      if (mask_array(iorb)) then
         icnt=icnt+1
         orbs_as%eval(icnt)=orbs%eval(iorb)
      end if
   end do
   if (icnt/=orbs_as%norb*orbs_as%nkpts) stop 'ERROR(select_active_space): icnt/=orbs_as%norb*orbs_as%nkpts'

   allocate(psi_as(orbs_as%npsidim_comp+ndebug),stat=i_stat)
   call memocc(i_stat,psi_as,'psi_as',subname)

   ispsi=1
   do ikptp=1,orbs%nkptsp
      ikpt=orbs%iskpts+ikptp
      nvctrp=comms%nvctr_par(iproc,ikpt) 
      !this should be identical in both the distributions
      if (nvctrp /= comms_as%nvctr_par(iproc,ikpt)) then
         write(*,*)'ERROR(select_active_space): the component distrbution is not identical'
         stop
      end if

      !put all the orbitals which match the active space
      ispsi=1
      ispsi_as=1
      do iorb=1,orbs%norb
         if (mask_array(iorb+(ikpt-1)*orbs%norb)) then
            call dcopy(nvctrp,psi(ispsi),1,psi_as(ispsi_as),1)
            ispsi_as=ispsi_as+nvctrp*orbs_as%nspinor
         end if
         ispsi=ispsi+nvctrp*orbs%nspinor
      end do
   end do

END SUBROUTINE select_active_space


!>   First orthonormalisation
subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit,orthpar)
   use module_base
   use module_types
   use module_interfaces, except_this_one_B => first_orthon
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(communications_arrays), intent(in) :: comms
   type(orthon_data):: orthpar
   real(wp), dimension(:) , pointer :: psi,hpsi,psit
   !local variables
   character(len=*), parameter :: subname='first_orthon'
   integer :: i_stat

   !!!  if(nspin==4) then
   !!!     nspinor=4
   !!!  else
   !!!     nspinor=1
   !!!  end if

   if (nproc > 1) then
      !allocate hpsi array (used also as transposed)
      !allocated in the transposed way such as 
      !it can also be used as the transposed hpsi
      allocate(hpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,hpsi,'hpsi',subname)
      !allocate transposed principal wavefunction
      allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,psit,'psit',subname)
   else
      psit => psi
   end if

   !to be substituted, must pass the wavefunction descriptors to the routine
   call transpose_v(iproc,nproc,orbs,wfd,comms,psi,&
      &   work=hpsi,outadd=psit(1))

   call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar)

   !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

   call untranspose_v(iproc,nproc,orbs,wfd,comms,psit,&
      &   work=hpsi,outadd=psi(1))

   if (nproc == 1) then
      nullify(psit)
      !allocate hpsi array
      allocate(hpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,hpsi,'hpsi',subname)
   end if

END SUBROUTINE first_orthon


!>   Transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
   use module_base
   use module_types
   use module_interfaces, fake_name => last_orthon
   implicit none
   integer, intent(in) :: iproc,nproc,iter
   real(wp), intent(out) :: evsum
   type(DFT_wavefunction), intent(inout) :: wfn
   logical, optional :: opt_keeppsit
   !local variables
   logical :: keeppsit
   character(len=*), parameter :: subname='last_orthon'
   integer :: i_all,i_stat
   real(wp), dimension(:,:,:), pointer :: mom_vec

   if (present(opt_keeppsit)) then
      keeppsit=opt_keeppsit
   else
      keeppsit=.false.
   end if

   call transpose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%hpsi,work=wfn%psi)
   if (nproc==1) then
      wfn%psit => wfn%psi
      call transpose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%psit)
   end if

   call subspace_diagonalisation(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%hpsi,evsum)

   !here we should preserve hpsi and transpose it if we are in ensemble mimimization scheme

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit,work=wfn%hpsi,outadd=wfn%psi(1))
   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   if(.not.  keeppsit) then
      if (nproc > 1  ) then
         i_all=-product(shape(wfn%psit))*kind(wfn%psit)
         deallocate(wfn%psit,stat=i_stat)
         call memocc(i_stat,i_all,'psit',subname)
      else
         nullify(wfn%psit)
      end if

      i_all=-product(shape(wfn%hpsi))*kind(wfn%hpsi)
      deallocate(wfn%hpsi,stat=i_stat)
      call memocc(i_stat,i_all,'hpsi',subname)

   endif

   !call eigensystem_info(iproc,nproc,wfn%Lzd%Glr%wfd%nvctr_c+7*wfn%Lzd%Glr%wfd%nvctr_f,&
   !     wfn%orbs,wfn%psi)

END SUBROUTINE last_orthon

subroutine eigensystem_info(iproc,nproc,tolerance,nvctr,orbs,psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nvctr
  real(gp), intent(in) :: tolerance !< threshold to classify degenerate eigenvalues
  type(orbitals_data), intent(inout) :: orbs
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='eigensystem_info'
  integer :: i_all,i_stat
  real(wp), dimension(:,:,:), pointer :: mom_vec


  !for a non-collinear treatment,
  !we add the calculation of the moments for printing their value
  !close to the corresponding eigenvector
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          nvctr,&
          orbs%nspinor,psi,mom_vec)
  else
     nullify(mom_vec)   
  end if

  ! Send all eigenvalues to all procs.
  call broadcast_kpt_objects(nproc,orbs%nkpts,orbs%norb, &
       orbs%eval(1),orbs%ikptproc)

  !here the new occupation numbers should be recalculated for future needs

  !print the found eigenvalues
  if (iproc == 0) then
     call write_eigenvalues_data(nproc,tolerance,orbs,mom_vec)
  end if

  if (orbs%nspinor ==4) then
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if
end subroutine eigensystem_info


!> Finds the fermi level ef for an error function distribution with a width wf
!! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
subroutine evaltoocc(iproc,nproc,filewrite,wf,orbs,occopt)
   use module_base
   use module_types
   implicit none
   logical, intent(in) :: filewrite
   integer, intent(in) :: iproc, nproc
   integer, intent(in) :: occopt      
   real(gp), intent(in) :: wf
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   integer :: ikpt,iorb,melec,ii,ierr
   real(gp) :: charge, chargef
   real(gp) :: ef,pi,electrons,dlectrons,factor,arg,argu,argd,corr,cutoffu,cutoffd,diff,full,res,resu,resd
   parameter(pi=3.1415926535897932d0)
   real(gp)  ::a, x, xu, xd, f, df, tt  
   real(gp)  ::sqrtpi ; parameter (sqrtpi=sqrt(pi)) 
   !write(*,*)  'ENTER Fermilevel',orbs%norbu,orbs%norbd

   orbs%eTS=0.0_gp  

   a = 0.d0
   select case (occopt)
   case  (SMEARING_DIST_ERF  )
   case  (SMEARING_DIST_FERMI)
   case  (SMEARING_DIST_COLD1) !Marzari's cold smearing  with a=-.5634 (bumb minimization)
      a=-.5634d0
   case  (SMEARING_DIST_COLD2) !Marzari's cold smearing  with a=-.8165 (monotonic tail)
      a=-.8165d0
   case  (SMEARING_DIST_METPX) !Methfessel and Paxton (same as COLD with a=0)
      a=0.d0
   case default
      if(iproc==0) print*, 'unrecognized occopt=', occopt
      stop 
   end select

   if (orbs%norbd==0) then 
      full=2.d0   ! maximum occupation for closed shell  orbital
   else
      full=1.d0   ! maximum occupation for spin polarized orbital
   endif

   if (orbs%nkpts.ne.1 .and. filewrite) stop 'Fermilevel: CANNOT write input.occ with more than one k-point'
   charge=0.0_gp
   do ikpt=1,orbs%nkpts
      !number of zero orbitals for the given k-point
      !overall charge of the system
      do iorb=1,orbs%norb
         charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb) * orbs%kwgts(ikpt)
      end do
   end do
   melec=nint(charge)
   !if (iproc == 0) write(*,*) 'charge',charge,melec

   ! Send all eigenvalues to all procs (presumably not necessary)
   call broadcast_kpt_objects(nproc, orbs%nkpts, orbs%norb, &
      &   orbs%eval(1), orbs%ikptproc)
   
   if (wf > 0.0_gp) then
      ii=0
      if (orbs%efermi == UNINITIALIZED(orbs%efermi)) then
         !last value as a guess
         orbs%efermi = orbs%eval(orbs%norbu)
         ! Take initial value at gamma point.
         do iorb = 1, orbs%norbu
            if (orbs%occup(iorb) < 1.0_gp) then
               orbs%efermi = orbs%eval(iorb)
               exit
            end if
         end do
      end if
      ef=orbs%efermi
      factor=1.d0/(sqrt(pi)*wf)
      !print *,0,ef

      ! electrons is N_electons = sum f_i * Wieght_i
      ! dlectrons is dN_electrons/dEf =dN_electrons/darg * darg/dEf= sum df_i/darg /(-wf) , darg/dEf=-1/wf
      !  f:= occupation # for band i ,  df:=df/darg
      loop_fermi: do
         ii=ii+1
         if (ii > 10000) stop 'error Fermilevel'
         electrons=0.d0
         dlectrons=0.d0
         do ikpt=1,orbs%nkpts
            do iorb=1,orbs%norbd+orbs%norbu
               arg=(orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf
               if (occopt == SMEARING_DIST_ERF) then
                  call derf_ab(res,arg)
                  f =.5d0*(1.d0-res)
                  df=-exp(-arg**2)/sqrtpi 
               else if (occopt == SMEARING_DIST_FERMI) then
                  f =1.d0/(1.d0+exp(arg)) 
                  df=-1.d0/(2.d0+exp(arg)+exp(-arg)) 
               else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &  
                    &  occopt == SMEARING_DIST_METPX ) then
                  x= -arg
                  call derf_ab(res,x)
                  f =.5d0*(1.d0+res +exp(-x**2)*(-a*x**2 + .5d0*a+x)/sqrtpi)
                  df=-exp(-x**2) * (a*x**3 -x**2 -1.5d0*a*x +1.5d0) /sqrtpi   ! df:=df/darg=-df/dx
               else
                  f  = 0.d0
                  df = 0.d0
               end if
               electrons=electrons+ f  * orbs%kwgts(ikpt)  ! electrons := N_e(Ef+corr.)
               dlectrons=dlectrons+ df * orbs%kwgts(ikpt)  ! delectrons:= dN_e/darge ( Well! later we need dN_e/dEf=-1/wf*dN_e/darg
               !if(iproc==0) write(*,*) arg,   f , df
            enddo
         enddo

         dlectrons=dlectrons/(-wf)  ! df/dEf=df/darg * -1/wf
         diff=-real(melec,gp)/full+electrons
         if (abs(diff) < 1.d-12) exit loop_fermi
         if (abs(dlectrons) <= 1d-45) then
            corr=wf
         else
            corr=diff/abs(dlectrons) ! for case of no-monotonic func. abs is needed
         end if
         !if (iproc==0) write(*,'(i5,3(1pe17.8),i4,1pe17.8)') ii,electrons,ef,dlectrons,melec,corr
         if (corr > 1.d0*wf) corr=1.d0*wf
         if (corr < -1.d0*wf) corr=-1.d0*wf
         if (abs(dlectrons) < 1.d-18  .and. electrons > real(melec,gp)/full) corr=3.d0*wf
         if (abs(dlectrons) < 1.d-18  .and. electrons < real(melec,gp)/full) corr=-3.d0*wf
         ef=ef-corr  ! Ef=Ef_guess+corr.
         !call MPI_BARRIER(MPI_COMM_WORLD,ierr) !debug
      end do loop_fermi

      do ikpt=1,orbs%nkpts
         argu=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu)-ef)/wf
         argd=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu+orbs%norbd)-ef)/wf
         if (occopt == SMEARING_DIST_ERF) then
            !error function
            call derf_ab(resu,argu)
            call derf_ab(resd,argd)
            cutoffu=.5d0*(1.d0-resu)
            cutoffd=.5d0*(1.d0-resd)
         else if (occopt == SMEARING_DIST_FERMI) then
            !Fermi function
            cutoffu=1.d0/(1.d0+exp(argu))
            cutoffd=1.d0/(1.d0+exp(argd))
         else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &  
              &  occopt == SMEARING_DIST_METPX ) then
            !Marzari's relation with different a 
            xu=-argu
            xd=-argd
            call derf_ab(resu,xu)
            call derf_ab(resd,xd)
            cutoffu=.5d0*(1.d0+resu +exp(-xu**2)*(-a*xu**2 + .5d0*a+xu)/sqrtpi)
            cutoffd=.5d0*(1.d0+resd +exp(-xd**2)*(-a*xd**2 + .5d0*a+xd)/sqrtpi)
         end if
      enddo
      if (iproc==0) write(*,'(1x,a,1pe21.14,2(1x,e8.1))') 'Fermi level, Fermi distribution cut off at:  ',ef,cutoffu,cutoffd
      orbs%efermi=ef

      !update the occupation number
      do ikpt=1,orbs%nkpts
         do iorb=1,orbs%norbu + orbs%norbd
            arg=(orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf
            if (occopt == SMEARING_DIST_ERF) then
               call derf_ab(res,arg)
               f=.5d0*(1.d0-res)
            else if (occopt == SMEARING_DIST_FERMI) then
               f=1.d0/(1.d0+exp(arg))
            else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &  
                 &  occopt == SMEARING_DIST_METPX ) then
               x=-arg
               call derf_ab(res,x)
               f =.5d0*(1.d0+res +exp(-x**2)*(-a*x**2 + .5d0*a+x)/sqrtpi)
            end if
            orbs%occup((ikpt-1)*orbs%norb+iorb)=full* f 
            !if(iproc==0) print*,  orbs%eval((ikpt-1)*orbs%norb+iorb), orbs%occup((ikpt-1)*orbs%norb+iorb)
         end do
      end do

      !update electronic entropy S; eTS=T_ele*S is the electronic entropy term the negative of which is added to energy: Free energy = energy-T*S 
      orbs%eTS=0.0_gp
      do ikpt=1,orbs%nkpts
         do iorb=1,orbs%norbu + orbs%norbd
            if (occopt == SMEARING_DIST_ERF) then
               !error function
               orbs%eTS=orbs%eTS+full*wf/(2._gp*sqrt(pi))*exp(-((orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf)**2)
            else if (occopt == SMEARING_DIST_FERMI) then
               !Fermi function
               tt=orbs%occup((ikpt-1)*orbs%norb+iorb)
               orbs%eTS=orbs%eTS-full*wf*(tt*log(tt) + (1._gp-tt)*log(1._gp-tt))
            else if (occopt == SMEARING_DIST_COLD1 .or. occopt == SMEARING_DIST_COLD2 .or. &  
                 &  occopt == SMEARING_DIST_METPX ) then
               !cold 
               orbs%eTS=orbs%eTS+0._gp  ! to be completed if needed                                             
            end if
         end do
      end do
      ! Sanity check on sum of occup.
      chargef=0.0_gp
      do ikpt=1,orbs%nkpts
         do iorb=1,orbs%norb
            chargef=chargef+orbs%kwgts(ikpt) * orbs%occup(iorb+(ikpt-1)*orbs%norb)
         end do
      end do
      if (abs(charge - chargef) > 1e-6)  stop 'error occupation update'
   else if(full==1.0_gp) then
      call eFermi_nosmearing(iproc,orbs)
      ! no entropic term when electronc temprature is zero
   end if

   !write on file the results if needed
   if (filewrite) then
      open(unit=11,file='input.occ',status='unknown')
      write(11,*)orbs%norbu,orbs%norbd
      do iorb=1,orbs%norb
         !write(11,'(i5,e19.12)')iorb,orbs%occup((ikpt-1)*orbs%norb+iorb)
         !    write(11,'(i5,e19.12)')iorb,full/(1.d0+exp(arg))  !,orbs%eval((ikpt-1)*orbs%norb+iorb)
         write(11,'(i5,e19.12,f10.6)')iorb,orbs%occup((ikpt-1)*orbs%norb+iorb) &
              &   ,orbs%eval ((ikpt-1)*orbs%norb+iorb)
      end do
      close(unit=11)
   end if

 END SUBROUTINE evaltoocc



subroutine eFermi_nosmearing(iproc,orbs)
   use module_base
   use module_types
   use yaml_output
   implicit none
   integer, intent(in) :: iproc
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   integer :: iu,id,n,nzeroorbs,ikpt,iorb
   real(gp) :: charge
   real(wp) :: eF

   iu=0
   id=0
   eF = 0._wp
   do ikpt=1,orbs%nkpts
      !number of zero orbitals for the given k-point
      nzeroorbs=0
      !overall charge of the system
      charge=0.0_gp
      do iorb=1,orbs%norb
         if (orbs%occup(iorb+(ikpt-1)*orbs%norb) == 0.0_gp) then
            nzeroorbs=nzeroorbs+1
         else
            charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb)
         end if
      end do
      if (nzeroorbs /= 0 .and. orbs%norbd .gt.0) then
         do iorb=1,orbs%norbu-1
            if (orbs%eval((ikpt-1)*orbs%norb+iorb) > orbs%eval((ikpt-1)*orbs%norb+iorb+1)) &
               &   write(*,*) 'wrong ordering of up EVs',iorb,iorb+1
         end do
         do iorb=1,orbs%norbd-1
            if (orbs%eval((ikpt-1)*orbs%norb+iorb+orbs%norbu) > orbs%eval((ikpt-1)*orbs%norb+iorb+1+orbs%norbu))&
               &   write(*,*) 'wrong ordering of dw EVs',iorb+orbs%norbu,iorb+1+orbs%norbu
         enddo

         iu=0
         id=0
         n=0
         do while (real(n,gp) < charge)
            if (orbs%eval((ikpt-1)*orbs%norb+iu+1) <= orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)) then
               iu=iu+1
               eF=orbs%eval((ikpt-1)*orbs%norb+iu+1)
            else
               id=id+1
               eF=orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)
            endif
            n=n+1
         enddo
         if (iproc==0) then
            !write(*,'(1x,a,1pe21.14,a,i4)') 'Suggested Homo energy level',eF,', Spin polarization',iu-id
            call yaml_map('Suggested Fermi Level',ef,fmt='(1pe21.14)')
            call yaml_map('Suggested Spin pol.',iu-id,fmt='(i4)')
         end if
         !write(*,*) 'up,down, up-down',iu,id,iu-id
      end if
   end do
   orbs%efermi=eF
   !assign the values for the occupation numbers
   do iorb=1,iu
      orbs%occup(iorb)=1.0_gp
   end do
   do iorb=iu+1,orbs%norbu
      orbs%occup(iorb)=0.0_gp
   end do
   do iorb=1,id
      orbs%occup(iorb+orbs%norbu)=1.0_gp
   end do
   do iorb=id+1,orbs%norbd
      orbs%occup(iorb+orbs%norbu)=0.0_gp
   end do

END SUBROUTINE eFermi_nosmearing


!>   Calculate magnetic moments
subroutine calc_moments(iproc,nproc,norb,norb_par,nvctr,nspinor,psi,mom_vec)
   use module_base
   implicit none
   integer, intent(in) :: iproc,nproc,norb,nvctr,nspinor
   integer, dimension(0:nproc-1), intent(in) :: norb_par
   real(wp), dimension(nvctr,norb*nspinor), intent(in) :: psi
   real(wp), dimension(4,norb,min(nproc,2)), intent(out) :: mom_vec
   !local variables
   character(len=*), parameter :: subname='calc_moments'
   integer :: i_all,i_stat,ierr,iorb,jproc
   integer :: ndim,oidx
   integer, dimension(:), allocatable :: norb_displ
   real(wp) :: m00,m11,m13,m24,m14,m23
   !real(wp), dimension(:,:,:), allocatable :: mom_vec

   ndim=2
   if (nproc==1) ndim=1

   if(nspinor==4) then

      call razero(4*norb*ndim,mom_vec)

      do iorb=1,norb_par(iproc)
         oidx=(iorb-1)*nspinor+1
         m00=dot(2*nvctr,psi(1,oidx),1,psi(1,oidx),1)
         m11=dot(2*nvctr,psi(1,oidx+2),1,psi(1,oidx+2),1)
         m13=dot(nvctr,psi(1,oidx),1,psi(1,oidx+2),1)
         m24=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+3),1)
         !        m12=dot(nvctr,psi(1,oidx),1,psi(1,oidx+1),1)
         !        m34=dot(nvctr,psi(1,oidx+2),1,psi(1,oidx+3),1)
         m14=dot(nvctr,psi(1,oidx),1,psi(1,oidx+3),1)
         m23=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+2),1)

         mom_vec(1,iorb,ndim)=(m00+m11) !rho
         mom_vec(2,iorb,ndim)=2.0d0*(m13+m24)       !m_x
         !        mom_vec(3,iorb,ndim)=2.0d0*(m12-m34)       !m_y
         mom_vec(3,iorb,ndim)=2.0d0*(m14-m23)       !m_y
         mom_vec(4,iorb,ndim)=(m00-m11)             !m_z
      end do

      if(nproc>1) then
         allocate(norb_displ(0:nproc-1+ndebug),stat=i_stat)
         call memocc(i_stat,norb_displ,'norb_displ',subname)

         norb_displ(0)=0
         do jproc=1,nproc-1
            norb_displ(jproc)=norb_displ(jproc-1)+norb_par(jproc-1)
         end do

         call MPI_GATHERV(mom_vec(1,1,2),4*norb_par(iproc),mpidtypw,&
            &   mom_vec(1,1,1),4*norb_par,4*norb_displ,mpidtypw,&
         0,MPI_COMM_WORLD,ierr)

         i_all=-product(shape(norb_displ))*kind(norb_displ)
         deallocate(norb_displ,stat=i_stat)
         call memocc(i_stat,i_all,'norb_displ',subname)
      end if

   end if

END SUBROUTINE calc_moments


subroutine check_communications(iproc,nproc,orbs,lr,comms)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(locreg_descriptors), intent(in) :: lr
   type(communications_arrays), intent(in) :: comms
   !local variables
   character(len=*), parameter :: subname='check_communications'
   integer :: i,ispinor,iorb,indspin,indorb,jproc,i_stat,i_all,iscomp,idsx,index,ikptsp
   integer :: ikpt,ispsi,nspinor,nvctrp,ierr
   real(wp) :: psival,maxdiff
   real(wp), dimension(:), allocatable :: psi
   real(wp), dimension(:), pointer :: pwork
   real(wp) :: epsilon
   character(len = 25) :: filename
   logical :: abort

   !allocate the "wavefunction" amd fill it, and also the workspace
   allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(pwork(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,pwork,'pwork',subname)


   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
         do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
            !vali=real(i,wp)/512.0_wp  ! *1.d-5
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            psi(i+indspin+indorb)=psival!(valorb+vali)*(-1)**(ispinor-1)
         end do
      end do
   end do

   !transpose the hpsi wavefunction
   call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psi,work=pwork)

   !check the results of the transposed wavefunction
   maxdiff=0.0_wp
   ispsi=0
   do ikptsp=1,orbs%nkptsp
      ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
      !valkpt=real(512*ikpt,wp)
      !calculate the starting point for the component distribution
      iscomp=0
      do jproc=0,iproc-1
         iscomp=iscomp+comms%nvctr_par(jproc,ikpt)
      end do
      nvctrp=comms%nvctr_par(iproc,ikpt)
      nspinor=orbs%nspinor

      do iorb=1,orbs%norb
         !valorb=real(iorb,wp)+valkpt
         indorb=(iorb-1)*nvctrp*nspinor
         do idsx=1,(nspinor-1)/2+1
            do i=1,nvctrp
               !vali=real(i+iscomp,wp)/512.d0  ! *1.d-5
               do ispinor=1,((2+nspinor)/4+1)
                  !psival=(-1)**(ispinor-1)*(valorb+vali)
                  call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                  !this is just to force the IEEE representation of psival
                  !              if (psival .lt. 0.d0) then  
                  !              write(321,*) psival,psival**2
                  !              endif
                  index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                     &   (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                  maxdiff=max(abs(psi(index)-psival),maxdiff)
               end do
            end do
         end do
      end do
      ispsi=ispsi+nvctrp*orbs%norb*nspinor
   end do

   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
      write(*,*)'ERROR: process',iproc,'does not transpose wavefunctions correctly!'
      write(*,*)'       found an error of',maxdiff,'cannot continue.'
      write(*,*)'       data are written in the file transerror.log, exiting...'

      write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
      open(unit=22,file=trim(filename),status='unknown')
      ispsi=0
      do ikptsp=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
         !valkpt=real(512*ikpt,wp)
         !calculate the starting point for the component distribution
         iscomp=0
         do jproc=0,iproc-1
            iscomp=iscomp+comms%nvctr_par(jproc,ikpt)
         end do
         nvctrp=comms%nvctr_par(iproc,ikpt)
         nspinor=orbs%nspinor

         do iorb=1,orbs%norb
            !valorb=real(iorb,wp)+valkpt
            indorb=(iorb-1)*nvctrp*nspinor
            do idsx=1,(nspinor-1)/2+1
               do i=1,nvctrp
                  !vali=real(i+iscomp,wp)/512.d0  !*1.d-5
                  do ispinor=1,((2+nspinor)/4+1)
                     !psival=(-1)**(ispinor-1)*(valorb+vali)
                     call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                     index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                        &   (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                     maxdiff=abs(psi(index)-psival)
                     if (maxdiff > 0.d0) then
                        write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i+iscomp,iorb,ikpt,psival,&
                           &   psi(index),maxdiff
                     end if
                  end do
               end do
            end do
         end do
         ispsi=ispsi+nvctrp*orbs%norb*nspinor
      end do
      close(unit=22)
      abort = .true.
      write(filename, "(A,I0,A)") 'distscheme', iproc, '.log'
      open(unit=22,file=trim(filename),status='unknown')
      call print_distribution_schemes(22,nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      close(unit=22)
   end if

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   if (abort) then
      if (iproc == 0) call print_distribution_schemes(6,nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      call MPI_ABORT(MPI_COMM_WORLD,ierr)
   end if

   !retranspose the hpsi wavefunction
   call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,&
      &   psi,work=pwork)

   maxdiff=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
         do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
            !vali=real(i,wp)/512.d0  !*1.d-5
            !psival=(valorb+vali)*(-1)**(ispinor-1)
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            maxdiff=max(abs(psi(i+indspin+indorb)-psival),maxdiff)
         end do
      end do
   end do

   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
      write(*,*)'ERROR: process',iproc,'does not untranspose wavefunctions correctly!'
      write(*,*)'       found an error of',maxdiff,'cannot continue.'
      write(*,*)'       data are written in the file transerror.log, exiting...'

      write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
      open(unit=22,file=trim(filename),status='unknown')
      maxdiff=0.0_wp
      do iorb=1,orbs%norbp
         ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
         !valkpt=real(512*ikpt,wp)
         !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
         indorb=(iorb-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor
         do ispinor=1,orbs%nspinor
            indspin=(ispinor-1)*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
            do i=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
               !vali=real(i,wp)/512.d0  !*1.d-5
               !psival=(valorb+vali)*(-1)**(ispinor-1)
               call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
               maxdiff=abs(psi(i+indspin+indorb)-psival)
               if (maxdiff > 0.d0) then
                  write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i,iorb,orbs%isorb,psival,&
                     &   psi(ispinor+(i-1)*orbs%nspinor+indorb),maxdiff
               end if
            end do
         end do
      end do
      close(unit=22)
      abort = .true.
   end if

   if (abort) call MPI_ABORT(MPI_COMM_WORLD,ierr)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   i_all=-product(shape(psi))*kind(psi)
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all=-product(shape(pwork))*kind(pwork)
   deallocate(pwork,stat=i_stat)
   call memocc(i_stat,i_all,'pwork',subname)

END SUBROUTINE check_communications


!> define a value for the wavefunction which is dependent of the indices
subroutine test_value(ikpt,iorb,ispinor,icomp,val)
   use module_base
   implicit none
   integer, intent(in) :: ikpt,icomp,iorb,ispinor
   real(wp), intent(out) :: val
   !local variables
   real(wp) :: valkpt,valorb,vali

   ! recognizable pattern, for debugging
   ! valkpt=real(10000*(ikpt-1),wp)!real(512*ikpt,wp)
   ! valorb=real(iorb,wp)+valkpt
   ! vali=real(icomp,wp)*1.e-5_wp  !real(icomp,wp)/512.0_wp  ! *1.d-5
   !
   ! val=(valorb+vali)*(-1)**(ispinor-1)

   valkpt=real(512*ikpt,wp)
   valorb=real(iorb,wp)+valkpt
   vali=real(icomp,wp)/512.0_wp  ! *1.d-5

   val=(valorb+vali)*(-1)**(ispinor-1)

END SUBROUTINE test_value


subroutine broadcast_kpt_objects(nproc, nkpts, ndata, data, ikptproc)
   use module_base
   implicit none
   integer, intent(in) :: nproc, nkpts, ndata
   integer, dimension(nkpts), intent(in) :: ikptproc
   real(gp), dimension(ndata,nkpts), intent(inout) :: data

   integer :: ikpt, ierr

   if (nproc > 1) then
      do ikpt = 1, nkpts
         call MPI_BCAST(data(1,ikpt), ndata,mpidtypg, &
            &   ikptproc(ikpt), MPI_COMM_WORLD, ierr)
         !redundant barrier 
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end do
   end if
END SUBROUTINE broadcast_kpt_objects



!!subroutine minimize_by_orthogonal_transformation(iproc, nproc, orbs, wfd, comms, orthpar, E0, El, stepsize, hpsi, psi, derivative)
!!use module_base
!!use module_types
!!use module_interfaces
!!implicit none
!!
!!  integer, intent(in) :: iproc,nproc
!!  type(orbitals_data), intent(in) :: orbs
!!  type(wavefunctions_descriptors),intent(in):: wfd
!!  type(communications_arrays), intent(in) :: comms
!!  type(orthon_data), intent(in) :: orthpar
!!  !n(c) type(wavefunctions_descriptors), intent(in) :: wfd
!!  !real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(in) :: hpsi
!!  !real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
!!  real(8),intent(in):: E0, El
!!  real(8),intent(inout):: stepsize 
!!  real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(inout) :: hpsi
!!  real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(inout) :: psi
!!  real(8),intent(out):: derivative
!!  !local variables
!!  character(len=*), parameter :: subname='orthoconstraint'
!!  integer :: i_stat,i_all,ierr,iorb !n(c) ise
!!  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor,iiorb,jjorb,jorb,lwork,info,k
!!  real(dp) :: occ !n(c) tt
!!  integer, dimension(:,:), allocatable :: ndimovrlp
!!  real(wp), dimension(:), allocatable :: alag, work, rwork, eval
!!  real(wp),dimension(:,:),allocatable:: gmat
!!  complex(8),dimension(:,:),allocatable:: gmatc, tempmatc, tempmat2c, omatc
!!  complex(8),dimension(:),allocatable:: expDc
!!  real(8):: lstep, dfactorial
!!  complex(8):: ttc
!!  real(8),dimension(:),pointer:: psiwork
!!  real(8),dimension(:,:),allocatable:: omat
!!
!!  write(*,*) 'iproc, orbs%npsidim',iproc,max(orbs%npsidim_orbs,orbs%npsidim_comp)
!!
!!  !separate the orthogonalisation procedure for up and down orbitals 
!!  !and for different k-points
!!  call timing(iproc,'LagrM_comput  ','ON')
!!
!!  !number of components of the overlap matrix for parallel case
!!  !calculate the dimension of the overlap matrix for each k-point
!!  if (orbs%norbd > 0) then
!!     nspin=2
!!  else
!!     nspin=1
!!  end if
!!
!!  !number of components for the overlap matrix in wp-kind real numbers
!!
!!  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
!!  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)
!!
!!  call dimension_ovrlp(nspin,orbs,ndimovrlp)
!!
!!  allocate(alag(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
!!  call memocc(i_stat,alag,'alag',subname)
!!
!!  !put to zero all the k-points which are not needed
!!  call razero(ndimovrlp(nspin,orbs%nkpts),alag)
!!
!!  ! Transpose orbitals
!!  allocate(psiwork(size(psi)), stat=i_stat)
!!  call transpose_v(iproc, nproc, orbs, wfd, comms, psi, work=psiwork)
!!  call transpose_v(iproc, nproc, orbs, wfd, comms, hpsi, work=psiwork)
!!
!!  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
!!  ispsi=1
!!  do ikptp=1,orbs%nkptsp
!!     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!
!!     do ispin=1,nspin
!!
!!        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
!!             nvctrp,norb,norbs,ncomp,nspinor)
!!        if (nvctrp == 0) cycle
!!
!!        if(nspinor==1) then
!!           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
!!                max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
!!                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
!!        else
!!           !this part should be recheck in the case of nspinor == 2
!!           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
!!                max(1,ncomp*nvctrp), &
!!                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
!!                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
!!        end if
!!        ispsi=ispsi+nvctrp*norb*nspinor
!!     end do
!!  end do
!!
!!  if (nproc > 1) then
!!     call timing(iproc,'LagrM_comput  ','OF')
!!     call timing(iproc,'LagrM_commun  ','ON')
!!     call mpiallred(alag(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
!!     call timing(iproc,'LagrM_commun  ','OF')
!!     call timing(iproc,'LagrM_comput  ','ON')
!!  end if
!!
!!  !now each processors knows all the overlap matrices for each k-point
!!  !even if it does not handle it.
!!  !this is somehow redundant but it is one way of reducing the number of communications
!!  !without defining group of processors
!!
!!  i_stat=0
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          i_stat=i_stat+1
!!          if(iproc==0) write(1750,'(a,3i8,es20.8)') 'iorb, jorb, i_stat, alag(i_stat)', iorb, jorb, i_stat, alag(i_stat)
!!      end do
!!  end do 
!!
!!  ! Build the antisymmetric matrix "alag(iorb,jorb)-alag(jorb,iorb)"
!!  allocate(gmat(orbs%norb,orbs%norb), stat=i_stat)
!!  call memocc(i_stat, gmat, 'gmat', subname)
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          iiorb=(iorb-1)*orbs%norb+jorb
!!          jjorb=(jorb-1)*orbs%norb+iorb
!!          gmat(jorb,iorb)=alag(iiorb)-alag(jjorb)
!!      end do
!!  end do
!!
!!
!!
!!
!!  !Build the complex matrix -igmat
!!  allocate(gmatc(orbs%norb,orbs%norb), stat=i_stat)
!!  !call memocc(i_stat, gmatc, 'gmatc', subname) !memocc not working for complex
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
!!      end do
!!  end do 
!!
!!
!!
!!  ! Check whether hermitian
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          if(real(gmatc(jorb,iorb))/=real(gmatc(iorb,jorb)) .or. aimag(gmatc(jorb,iorb))/=-aimag(gmatc(iorb,jorb))) then
!!              write(*,'(a,4es16.7)') 'ERROR: (Gmatc(1,jorb,iorb)), (Gmatc(1,iorb,jorb)), (Gmatc(2,jorb,iorb)), (Gmatc(2,iorb,jorb))', Gmatc(jorb,iorb), Gmatc(iorb,jorb)
!!          end if
!!          if(iproc==0) write(1710,'(a,2i8,2es20.12)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
!!      end do
!!  end do 
!!
!!
!!
!!  ! Diagonalize gmatc
!!  allocate(eval(orbs%norb), stat=i_stat)
!!  call memocc(i_stat, eval, 'eval', subname)
!!  lwork=10*orbs%norb
!!  allocate(work(2*lwork), stat=i_stat) ! factor of 2 since it is assumed to be complex
!!  allocate(rwork(lwork), stat=i_stat)
!!  call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
!!  if(info/=0) stop 'ERROR in zheev'
!!  deallocate(work)
!!  deallocate(rwork)
!!
!!  do iorb=1,orbs%norb
!!      if(iproc==0) write(1720,'(a,i8,es20.8)') 'iorb, eval(iorb)', iorb, eval(iorb)
!!  end do
!!
!!
!!  if(stepsize>0.d0) then
!!      if(iproc==0) write(*,*) 'in first branch'
!!      ! Calculate the derivative
!!      ! Calculate optimal step size
!!      lstep=derivative*stepsize**2/(2.d0*(derivative*stepsize+E0-El))
!!      stepsize=-1.d0
!!      derivative=0.d0
!!  else
!!      if(iproc==0) write(*,*) 'in second branch'
!!      lstep=.1d0/(maxval(abs(eval)))
!!      stepsize=lstep
!!      derivative=0.d0
!!      do iorb=1,orbs%norb
!!          do jorb=1,orbs%norb
!!              derivative=derivative+gmat(jorb,iorb)**2
!!          end do
!!      end do
!!      derivative=-.5d0*derivative
!!  end if
!!  if(iproc==0) write(*,'(a,2es12.4)') '>> LSTEP:, derivative', lstep, derivative
!!
!!
!!  ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!  ! This is also a diagonal matrix, so only calculate the diagonal part.
!!  allocate(expDc(orbs%norb), stat=i_stat)
!!  !call memocc(i_stat, expDc, 'expDc', subname)
!!  do iorb=1,orbs%norb
!!     ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
!!     expDc(iorb)=cmplx(0.d0,0.d0,kind=8)
!!      do k=0,100
!!          expDc(iorb)=expDc(iorb)+ttc**k/dfactorial(k)
!!      end do
!!  end do
!!  do iorb=1,orbs%norb
!!     if(iproc==0) write(1740,'(a,i8,2es20.8)') 'iorb, expDc(iorb)', iorb, expDc(iorb)
!!  end do
!!
!!
!!
!!  ! Calculate the matrix O
!!  allocate(tempmatc(orbs%norb,orbs%norb), stat=i_stat)
!!  allocate(tempmat2c(orbs%norb,orbs%norb), stat=i_stat)
!!  allocate(omatc(orbs%norb,orbs%norb), stat=i_stat)
!!  do iorb=1,orbs%norb
!!     do jorb=1,orbs%norb
!!         if(iorb==jorb) then
!!             tempmat2c(jorb,iorb)=expDc(iorb)
!!         else
!!             tempmat2c(jorb,iorb)=cmplx(0.d0,0.d0,kind=8)
!!         end if
!!     end do
!!  end do
!!  call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmat2c(1,1), orbs%norb, &
!!      gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1), orbs%norb)
!!  do iorb=1,orbs%norb
!!     do jorb=1,orbs%norb
!!         if(iproc==0) write(1730,'(a,2i8,2es20.8)') 'iorb, jorb, tempmatc(jorb,iorb)', iorb, jorb, tempmatc(jorb,iorb)
!!     end do
!!  end do
!!
!!  call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!      tempmatc(1,1), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!
!!  if(iproc==0) then
!!     do iorb=1,orbs%norb
!!         do jorb=1,orbs%norb
!!             write(1700,'(2i8,2es20.8)') iorb,jorb,omatc(jorb,iorb)
!!         end do
!!     end do
!!  end if
!!
!!
!!  ! Build new linear combinations
!!  call memocc(i_stat, psiwork, 'psiwork', subname)
!!
!!  allocate(omat(orbs%norb,orbs%norb), stat=i_stat)
!!  call memocc(i_stat, omat, 'omat', subname)
!!
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          omat(jorb,iorb)=real(omatc(jorb,iorb))
!!          !!if(iorb==jorb) then
!!          !!    omat(jorb,iorb)=1.d0
!!          !!else
!!          !!    omat(jorb,iorb)=0.d0
!!          !!end if
!!      end do
!!  end do
!!
!!  nvctrp=comms%nvctr_par(iproc,0)
!!  call gemm('n', 'n', nvctrp, orbs%norb, orbs%norb, 1.0_wp, psi(1), max(1,nvctrp), &
!!            omat(1,1), orbs%norb, 0.0_wp, psiwork(1), max(1,nvctrp))
!!  call dcopy(size(psi), psiwork(1), 1, psi(1), 1)
!!
!!  ! I think this is not required..
!!  call orthogonalize(iproc,nproc,orbs,comms,psi,orthpar)
!!
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, psi, work=psiwork)
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, hpsi, work=psiwork)
!!
!!
!!
!!
!!
!!
!!
!!  i_all=-product(shape(omat))*kind(omat)
!!  deallocate(omat,stat=i_stat)
!!  call memocc(i_stat,i_all,'omat',subname)
!!
!!  i_all=-product(shape(psiwork))*kind(psiwork)
!!  deallocate(psiwork,stat=i_stat)
!!  call memocc(i_stat,i_all,'psiwork',subname)
!!
!!  i_all=-product(shape(alag))*kind(alag)
!!  deallocate(alag,stat=i_stat)
!!  call memocc(i_stat,i_all,'alag',subname)
!!
!!  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
!!  deallocate(ndimovrlp,stat=i_stat)
!!  call memocc(i_stat,i_all,'ndimovrlp',subname)
!!
!!  i_all=-product(shape(gmat))*kind(gmat)
!!  deallocate(gmat,stat=i_stat)
!!  call memocc(i_stat,i_all,'gmat',subname)
!!
!!  i_all=-product(shape(eval))*kind(eval)
!!  deallocate(eval,stat=i_stat)
!!  call memocc(i_stat,i_all,'eval',subname)
!!
!!  !!memoc not working for complex
!!  !i_all=-product(shape(gmatc))*kind(gmatc)  
!!  deallocate(gmatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'gmatc',subname)
!!
!!  !i_all=-product(shape(tempmatc))*kind(tempmatc)  
!!  deallocate(tempmatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'tempmatc',subname)
!!
!!  !i_all=-product(shape(tempmat2c))*kind(tempmat2c)  
!!  deallocate(tempmat2c,stat=i_stat)
!!  !call memocc(i_stat,i_all,'tempmat2c',subname)
!!
!!  !i_all=-product(shape(omatc))*kind(omatc)  
!!  deallocate(omatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'omatc',subname)
!!
!!  !i_all=-product(shape(expDc))*kind(expDc)  
!!  deallocate(expDc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'expDc',subname)
!!
!!  call timing(iproc,'LagrM_comput  ','OF')
!!
!!
!!
!!end subroutine minimize_by_orthogonal_transformation

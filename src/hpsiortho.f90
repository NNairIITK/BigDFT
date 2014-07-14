!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2013 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the application of the Hamiltonian on the wavefunction. The hamiltonian can be self-consistent or not.
!! In the latter case, the potential should be given in the rhov array of denspot structure. 
!! Otherwise, rhov array is filled by the self-consistent density
subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,iscf,alphamix,&
     nlpsp,rxyz,linflag,unblock_comms,GPU,wfn,&
     energs,rpnrm,xcstr,proj_G,paw)
  use module_base
  use module_types
  use module_interfaces, fake_name => psitohpsi
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use m_ab7_mixing
  use yaml_output
  use psp_projectors, only: PSPCODE_PAW,PSP_APPLY_SKIP
  implicit none
  !Arguments
  logical, intent(in) :: scf  !< If .false. do not calculate the self-consistent potential
  integer, intent(in) :: iproc,nproc,itrp,iscf,linflag,itwfn
  character(len=3), intent(in) :: unblock_comms
  real(gp), intent(in) :: alphamix
  type(atoms_data), intent(in) :: atoms
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(DFT_local_fields), intent(inout) :: denspot
  type(energy_terms), intent(inout) :: energs
  type(DFT_wavefunction), intent(inout) :: wfn
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(GPU_pointers), intent(inout) :: GPU  
  !real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(gp), intent(inout) :: rpnrm
  real(gp), dimension(6), intent(out) :: xcstr
  !real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
  type(gaussian_basis),dimension(atoms%astruct%nat),optional,intent(in)::proj_G
  type(paw_objects),optional,intent(inout)::paw
  !local variables
  character(len=*), parameter :: subname='psitohpsi'
  logical :: unblock_comms_den,unblock_comms_pot,whilepot,savefields
  integer :: nthread_max,ithread,nthread,irhotot_add,irho_add,ispin,i_all,i_stat,correcth
  !integer :: ii,jj
  !$ integer :: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads


  !in the default case, non local hamiltonian is done after potential creation
  whilepot=.true.

  !flag for saving the local fields (rho,vxc,vh)
  savefields= (iscf==SCF_KIND_GENERALIZED_DIRMIN)
  correcth=1
  !do not do that if rho_work is already associated
  if (savefields .and. associated(denspot%rho_work)) then
     !flag for correcting the hamiltonian (either is false or toggles savefields)
     correcth=2
     savefields=.false.
  end if

  nthread_max=1
  ithread=0
  nthread=1
  !control if we have the possibility of using OMP to de-synchronize communication
  !$ nthread_max=omp_get_max_threads()
  

  !decide the communication strategy
  unblock_comms_den=mpi_thread_funneled_is_supported .and. unblock_comms=='DEN' .and. &
      nthread_max > 1 .and. scf .and. .not. GPU%OCLconv! density is done only if scf is present
  if (unblock_comms_den) whilepot=.false. !anticipate the NlHamiltonian if density should be overlapped

  unblock_comms_pot=mpi_thread_funneled_is_supported .and. unblock_comms=='POT' .and. &
      nthread_max > 1 .and. whilepot .and. .not. GPU%OCLconv

  if ((unblock_comms_den .or. unblock_comms_pot) .and. iproc==0) then
     if (unblock_comms_den) call yaml_map('Overlapping communication of','Density')
     if (unblock_comms_pot) call yaml_map('Overlapping communication of','Potential')
     call yaml_map('No. of OMP Threads',nthread_max)
     call yaml_newline()
  end if

  !zero the output wavefunction
  if (wfn%orbs%npsidim_orbs > 0) call to_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))

  !calculate the self-consistent potential
  if (scf) then
     !update the entropic energy
     energs%eTS=wfn%orbs%eTS
     !safe the previous value of the energy
     energs%e_prev=energs%energy
     !print *,'here',savefields,correcth,energs%ekin,energs%epot,dot(wfn%orbs%npsidim_orbs,wfn%psi(1),1,wfn%psi(1),1)
     ! Potential from electronic charge density 
     call sumrho(denspot%dpbox,wfn%orbs,wfn%Lzd,GPU,atoms%astruct%sym,denspot%rhod,denspot%xc,wfn%psi,denspot%rho_psi)
     !print *,'here',wfn%orbs%occup(:),'there',savefields,correcth,energs%ekin,energs%epot,&
     !     dot(wfn%Lzd%Glr%d%n1i*wfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p*wfn%orbs%nspin,&
     !     denspot%rho_psi(1,1),1,denspot%rho_psi(1,1),1)

     !initialize nested approach 
     !this has always to be done for using OMP parallelization in the 
     !projector case
     !if nesting is not supported, a bigdft_nesting routine should not be called
     !$ if (unblock_comms_den) then
     !$ call timing(iproc,'UnBlockDen    ','ON')
     !$ call bigdft_open_nesting(2)
     !$ end if
     !print *,'how many threads ?' ,nthread_max
     !$OMP PARALLEL IF(unblock_comms_den) DEFAULT(shared), PRIVATE(ithread,nthread)
     !$ ithread=omp_get_thread_num()
     !$ nthread=omp_get_num_threads() !this should be 2 if active
     !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
     ! thread 0 does mpi communication 
     if (ithread == 0) then
        !$ if (unblock_comms_den) call OMP_SET_NUM_THREADS(1)
        !communicate density 
        call communicate_density(denspot%dpbox,wfn%orbs%nspin,denspot%rhod,&
             denspot%rho_psi,denspot%rhov,.false.)
        !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
     end if
     !in case of GPU do not overlap density communication and projectors
     if ((ithread > 0 .or. nthread==1) .and. .not. whilepot .and. .not. GPU%OCLconv) then
        ! Only the remaining threads do computations (if active) 
        !$ if (unblock_comms_den) call OMP_SET_NUM_THREADS(nthread_max-1)

        !nonlocal hamiltonian
        !$ if (verbose > 2 .and. iproc==0 .and. unblock_comms_den)&
        !$ & print *,'NonLocalHamiltonian with nthread:, out to:' ,omp_get_max_threads(),nthread_max
        call NL_ham()
     end if
     !$OMP END PARALLEL !if unblock_comms_den
     !$ if (unblock_comms_den) then
     !$ call bigdft_close_nesting(nthread_max)
     !$ call timing(iproc,'UnBlockDen    ','OF')
     !$ end if

     ithread=0
     nthread=1

     !here the density can be mixed
     if (iscf > SCF_KIND_DIRECT_MINIMIZATION) then
        if (denspot%mix%kind == AB7_MIXING_DENSITY) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3),&
                rpnrm,denspot%dpbox%nscatterarr)
           
           if (iproc == 0 .and. itrp > 1) then
              call yaml_newline()
              call yaml_map('itrp',itrp,fmt='(i4)')
              call yaml_map('Mixing on','Density')
              call yaml_map('RhoPot delta per volume unit',rpnrm,fmt='(1pe9.2)',&
                   label='rpnrm'//trim(adjustl(yaml_toa(itrp,fmt='(i4.4)'))))
              call yaml_newline()
              !write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
              !     &   'DENSITY iteration,Delta : (Norm 2/Volume)',itrp,rpnrm
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
           call vcopy(denspot%dpbox%ndimpot,&
                denspot%rhov(irhotot_add),1,denspot%rhov(irho_add),1)
           irhotot_add=irhotot_add+denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3d
           irho_add=irho_add+denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p
        end do
     end if

     if(wfn%orbs%nspinor==4) then
        !this wrapper can be inserted inside the XC_potential routine
        call PSolverNC(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
             denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
             denspot%dpbox%n3d,denspot%xc,&
             denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
             denspot%rhov,denspot%pkernel%kernel,denspot%V_ext,&
             energs%eh,energs%exc,energs%evxc,0.d0,.true.,4)
     else

!!           denspot%rhov denspot%rhov+2.e-7   STEFAN Goedecker
        call XC_potential(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
             denspot%pkernel%mpi_env%mpi_comm,&
             denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),denspot%xc,&
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
           call vcopy(denspot%dpbox%ndimpot,denspot%rhov(1),1,&
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
        if (denspot%mix%kind == AB7_MIXING_POTENTIAL) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3),&!volume should be used 
                rpnrm,denspot%dpbox%nscatterarr)
           if (iproc == 0 .and. itrp > 1) then
              call yaml_newline()
              call yaml_map('itrp',itrp,fmt='(i4)')
              call yaml_map('Mixing on','Potential')
              call yaml_map('RhoPot delta per volume unit',rpnrm,fmt='(1pe9.2)',&
                   label='rpnrm'//trim(adjustl(yaml_toa(itrp,fmt='(i4.4)'))))
              call yaml_newline()
              !write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
              !     &   'POTENTIAL iteration,Delta P (Norm 2/Volume)',itrp,rpnrm
           end if
        end if
     end if
     call denspot_set_rhov_status(denspot, KS_POTENTIAL, itwfn, iproc, nproc)

     if (savefields) then
        if (associated(denspot%rho_work)) then
           write(*,*)'ERROR: the reference potential should be empty to correct the hamiltonian!'
           stop
        end if

        denspot%rho_work = f_malloc_ptr(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim+ndebug,id='denspot%rho_work')
        call vcopy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,denspot%rhov(1),1,&
             denspot%rho_work(1),1)
     end if

     !if the hamiltonian should have to be just updated, perform the difference between potentials
     if (correcth==2) then
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
        call f_free_ptr(denspot%rho_work)
        nullify(denspot%rho_work)
     end if


  end if

  !debug
  !call MPI_BARRIER(MPI_COMM_WORLD,i_stat)
  !write(*,*)'hpsiortho l325,erase me'
  !ifile=iproc+140
  !jj=size(wfn%psi)
  !ii=size(wfn%hpsi)
  !write(ifile,*)'# ',jj,ii
  !do ii=1,jj
  ! write(ifile,'(2f20.8)')wfn%psi(ii),wfn%hpsi(ii)
  !end do
  !ifile=iproc+150
  !jj=size(denspot%rhov)
  !write(ifile,*)'# ',jj
  !do ii=1,jj
  ! write(ifile,'(2f20.8)')denspot%rhov(ii)
  !end do
  !call MPI_BARRIER(MPI_COMM_WORLD,i_stat)
  !end debug



  !non self-consistent case: rhov should be the total potential
  if (denspot%rhov_is /= KS_POTENTIAL) then
     stop 'psitohpsi: KS_potential not available' 
  end if

  !temporary, to be corrected with comms structure
  if (wfn%exctxpar == 'OP2P') energs%eexctX = UNINITIALIZED(1.0_gp)

  !initialize nested approach 
  !this has always to be done for using OMP parallelization in the 
  !projector case
  !if nesting is not supported, bigdft_nesting routine should not be called
  !$ if (unblock_comms_pot) then
  !$ call timing(iproc,'UnBlockPot    ','ON')
  !$ call bigdft_open_nesting(2)
  !$ end if
  !print *,'how many threads ?' ,nthread_max
  !$OMP PARALLEL IF (unblock_comms_pot) DEFAULT(shared), PRIVATE(ithread,nthread)
  !$ ithread=omp_get_thread_num()
  !$ nthread=omp_get_num_threads() !this should be 2 if active
  !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
  ! thread 0 does mpi communication 
  if (ithread == 0) then
     !$ if (unblock_comms_pot) call OMP_SET_NUM_THREADS(1)
     call full_local_potential(iproc,nproc,wfn%orbs,wfn%Lzd,linflag,&
          denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work)
     !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
  end if
  if ((ithread > 0 .or. nthread==1) .and. whilepot .and. .not. GPU%OCLconv) then
     ! Only the remaining threads do computations (if active) 
     !$ if (unblock_comms_pot) call OMP_SET_NUM_THREADS(nthread_max-1)

     !nonlocal hamiltonian
     !$ if (verbose > 2 .and. iproc==0 .and. unblock_comms_pot)&
     !$ & print *,'NonLocalHamiltonian with nthread:, out to:' ,omp_get_max_threads(),nthread_max
     call NL_ham()
  end if
  !$OMP END PARALLEL !if unblock_comms_pot
  !$ if (unblock_comms_pot) then
  !$ call bigdft_close_nesting(nthread_max)
  !$ call timing(iproc,'UnBlockPot    ','OF')
  !$ end if

  ithread=0
  nthread=1
 
  !here exctxpar might be passed
  !choose to just add the potential if needed
  call LocalHamiltonianApplication(iproc,nproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
       wfn%Lzd,wfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,wfn%psi,wfn%hpsi,&
       energs,wfn%SIC,GPU,correcth,denspot%xc,pkernel=denspot%pkernelseq)

  !in the case of OCL GPU the nonlocal hamiltonian can run after the local hamiltonian to overlap GPU-CPU computation
  if (GPU%OCLconv) call NL_ham()

  call SynchronizeHamiltonianApplication(nproc,wfn%orbs%npsidim_orbs,wfn%orbs,wfn%Lzd,&
       & GPU,denspot%xc,wfn%hpsi,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  ! Emit that hpsi are ready.
  if (wfn%c_obj /= 0) then
     call kswfn_emit_psi(wfn, itwfn, 1, iproc, nproc)
  end if

  !deallocate potential
  call free_full_potential(denspot%dpbox%mpi_env%nproc,linflag,denspot%xc,denspot%pot_work,subname)
  !----
  if (iproc==0 .and. verbose > 0) then
     if (correcth==2) then
        call yaml_map('Hamiltonian Applied',.false.)
        call yaml_comment('Only local potential')
     else
        call yaml_map('Hamiltonian Applied',.true.)
     end if
  end if

  contains 

    !> code factorization for the application of the nonlocal hamiltonian
    subroutine NL_ham()
      if(any(atoms%npspcode == PSPCODE_PAW)) then
         call to_zero(wfn%orbs%npsidim_orbs,paw%spsi(1))
         call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,rxyz,&
              wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,energs%eproj,proj_G,paw)
      else
         call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,rxyz,&
              wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,energs%eproj)
      end if
    end subroutine NL_ham

end subroutine psitohpsi


!> Application of the Full Hamiltonian
subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,rxyz,&
     Lzd,nlpsp,confdatarr,ngatherarr,pot,psi,hpsi,&
     energs,SIC,GPU,xc,pkernel,orbsocc,psirocc,proj_G,paw)
  use module_base
  use module_types
  use module_interfaces, fake_name => FullHamiltonianApplication
  use module_xc
  use psp_projectors, only: PSPCODE_PAW,PSP_APPLY_SKIP
  implicit none
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(SIC_data), intent(in) :: SIC
  type(xc_info), intent(in) :: xc
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  !real(wp), dimension(lzd%ndimpotisf) :: pot
  real(wp), dimension(:),pointer :: pot
  type(energy_terms), intent(inout) :: energs
  !real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,evsic
  real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  type(coulomb_operator), intent(in), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !PAW variables:
  type(gaussian_basis),dimension(at%astruct%ntypes),optional,intent(in)::proj_G
  type(paw_objects),optional,intent(inout)::paw

  !put to zero hpsi array (now important since any of the pieces of the hamiltonian is accumulating)
  if (orbs%npsidim_orbs > 0) call to_zero(orbs%npsidim_orbs,hpsi(1))

  !write(*,*) 'lzd%ndimpotisf', lzd%ndimpotisf
  !do i=1,lzd%ndimpotisf
  !    write(210,*) pot(i)
  !end do


 if (.not. present(pkernel)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
         energs,SIC,GPU,1,xc)
 else if (present(pkernel) .and. .not. present(orbsocc)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,1,xc,pkernel=pkernel)
 else if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
    call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
         Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
         energs,SIC,GPU,1,xc,pkernel,orbsocc,psirocc)
 else
    stop 'HamiltonianApplication, argument error'
 end if

 !these two sections have to be inverted to profit of overlapping in GPU accelerated case
 if(any(at%npspcode == PSPCODE_PAW)) then
  call NonLocalHamiltonianApplication(iproc,at,orbs%npsidim_orbs,orbs,rxyz,&
       Lzd,nlpsp,psi,hpsi,energs%eproj,proj_G,paw)
 else
  call NonLocalHamiltonianApplication(iproc,at,orbs%npsidim_orbs,orbs,rxyz,&
       Lzd,nlpsp,psi,hpsi,energs%eproj)
 end if

  call SynchronizeHamiltonianApplication(nproc,orbs%npsidim_orbs,orbs,Lzd,GPU,xc,hpsi,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)



END SUBROUTINE FullHamiltonianApplication


!> Application of the Local Hamiltonian
subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
     Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
     energs,SIC,GPU,PotOrKin,xc,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
   use module_base
   use module_types
   use module_xc
   use module_interfaces, except_this_one => LocalHamiltonianApplication
   use yaml_output
   use communications_base, only: p2pComms
   implicit none
   !logical, intent(in) :: onlypot !< if true, only the potential operator is applied
   integer, intent(in) :: PotOrKin
   integer, intent(in) :: iproc,nproc,npsidim_orbs
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd 
   type(SIC_data), intent(in) :: SIC
   type(xc_info), intent(in) :: xc
   integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
   real(wp), dimension(npsidim_orbs), intent(in) :: psi
   type(confpot_data), dimension(orbs%norbp) :: confdatarr
   real(wp), dimension(:), pointer :: pot
   !real(wp), dimension(*) :: pot
   type(energy_terms), intent(inout) :: energs
   real(wp), target, dimension(max(1,npsidim_orbs)), intent(inout) :: hpsi
   type(GPU_pointers), intent(inout) :: GPU
   type(coulomb_operator), intent(in), optional :: pkernel
   type(orbitals_data), intent(in), optional :: orbsocc
   real(wp), dimension(:), pointer, optional :: psirocc
   type(denspot_distribution),intent(in),optional :: dpbox
   !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
   real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
   type(p2pComms),intent(inout), optional:: comgp
   real(wp), target, dimension(max(1,npsidim_orbs)), intent(inout),optional :: hpsi_noconf
   real(gp),intent(out),optional :: econf
   !local variables
   character(len=*), parameter :: subname='HamiltonianApplication'
   logical :: exctX,op2p
   integer :: i_stat,n3p,ispot,ipotmethod
   real(gp) :: evsic_tmp
   type(coulomb_operator) :: pkernelSIC

   call f_routine(id='LocalHamiltonianApplication')

   ! local potential and kinetic energy for all orbitals belonging to iproc
   !if (iproc==0 .and. verbose > 1) then
      !call yaml_comment('Hamiltonian application, ',advance='no')
      !write(*,'(1x,a)',advance='no')&
      !     'Hamiltonian application...'
   !end if

   !initialise exact exchange energy 
   op2p=(energs%eexctX == UNINITIALIZED(1.0_gp))
   energs%eexctX=0.0_gp
   energs%evsic=0.0_gp
   evsic_tmp=0.0_gp

   exctX = xc_exctXfac(xc) /= 0.0_gp

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
         call exact_exchange_potential_virt(iproc,nproc,at%astruct%geocode,orbs%nspin,&
              Lzd%Glr,orbsocc,orbs,ngatherarr(0,1),n3p,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernel,psirocc,psi,pot(ispot))
         energs%eexctX = 0._gp
      else
         !here the condition for the scheme should be chosen
         if (.not. op2p) then
            call exact_exchange_potential(iproc,nproc,at%astruct%geocode,xc,orbs%nspin,&
                 Lzd%Glr,orbs,ngatherarr(0,1),n3p,&
                 0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
                 pkernel,psi,pot(ispot),energs%eexctX)
         else
            !the psi should be transformed in real space
            call exact_exchange_potential_round(iproc,nproc,xc,orbs%nspin,Lzd%Glr,orbs,&
                0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
                pkernel,psi,pot(ispot),energs%eexctX)

            !call exact_exchange_potential_op2p(iproc,nproc,xc,Lzd%Glr,orbs,pkernel,psi,pot(ispot),energs%eexctX)

         end if
      end if
      !print *,'iproc,eexctX',iproc,eexctX
   else if (ipotmethod==3) then
      !put fref=1/2 for the moment
      if (present(orbsocc) .and. present(psirocc)) then
         call NK_SIC_potential(Lzd%Glr,orbs,xc,SIC%fref,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernelSIC,psi,pot(ispot),evsic_tmp,&
              potandrho=psirocc)
      else
         call NK_SIC_potential(Lzd%Glr,orbs,xc,SIC%fref,&
              0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
              pkernelSIC,psi,pot(ispot),evsic_tmp)
      end if
   end if

   !GPU are supported only for ipotmethod=0
   if ((GPUconv .or. GPU%OCLconv) .and. ipotmethod /=0) then
      if (iproc ==0) write(*,*)&
         &   'ERROR(HamiltonianApplication): Accelerated hamiltonian are possible only with ipotmethod==0)'
      stop
   end if


   !apply the local hamiltonian for each of the orbitals
   !given to each processor
   !pot=0.d0
   !psi=1.d0
   !switch between GPU/CPU treatment
   !  do i=1,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
   !       call random_number(psi(i))
   !  end do
   if(GPU%OCLconv .or. GPUconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
      GPU%hpsi_ASYNC = f_malloc_ptr(max(1, (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),id='GPU%hpsi_ASYNC')
!      call to_zero((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp,GPU%hpsi_ASYNC(1))!hpsi(1))
   !else if (GPU%OCLconv) then
   !   GPU%hpsi_ASYNC => hpsi
   end if
   if (GPUconv) then
      call timing(iproc,'ApplyLocPotKin','ON') 
      call local_hamiltonian_GPU(orbs,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
           orbs%nspin,pot,psi,GPU%hpsi_ASYNC,energs%ekin,energs%epot,GPU)
      call timing(iproc,'ApplyLocPotKin','OF') 
   else if (GPU%OCLconv) then

      !pin potential
      call timing(iproc,'ApplyLocPotKin','ON') 
      call local_hamiltonian_OCL(orbs,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
           orbs%nspin,pot,psi,GPU%hpsi_ASYNC,energs%ekin,energs%epot,GPU)
      call timing(iproc,'ApplyLocPotKin','OF') 
   else

!!$      !temporary allocation
!!$      allocate(fake_pot(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+ndebug),stat=i_stat)
!!$      call memocc(i_stat,fake_pot,'fake_pot',subname)
!!$
!!$      call to_zero(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin,fake_pot(1))

      !local hamiltonian application for different methods
      !print *,'here',ipotmethod,associated(pkernelSIC)
      if (PotOrKin==1) then ! both
         call timing(iproc,'ApplyLocPotKin','ON') 
         if(present(dpbox) .and. present(potential) .and. present(comgp)) then
            call local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
                 ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                 xc,SIC%alpha,energs%ekin,energs%epot,energs%evsic,&
                 dpbox,potential,comgp)
         else
            call local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
                 ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                 xc,SIC%alpha,energs%ekin,energs%epot,energs%evsic)
         end if
         call timing(iproc,'ApplyLocPotKin','OF') 
!!$      i_all=-product(shape(fake_pot))*kind(fake_pot)
!!$      deallocate(fake_pot,stat=i_stat)
!!$      call memocc(i_stat,i_all,'fake_pot',subname)
         
      else if (PotOrKin==2) then !only pot
         call timing(iproc,'ApplyLocPot','ON') 
         if (present(hpsi_noconf)) then
             if (.not.present(econf)) then
                 stop 'ERROR: econf must be present when hpsi_noconf is present'
             end if
             call psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
                  ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                  xc,SIC%alpha,energs%epot,energs%evsic,hpsi_noconf,econf)
         else
             call psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
                  ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                  xc,SIC%alpha,energs%epot,energs%evsic)
         end if
         call timing(iproc,'ApplyLocPot','OF') 
      else if (PotOrKin==3) then !only kin
         call timing(iproc,'ApplyLocKin','ON') 
         call psi_to_kinpsi(iproc,npsidim_orbs,orbs,lzd,psi,hpsi,energs%ekin)
         call timing(iproc,'ApplyLocKin','OF') 
      end if

      !sum the external and the BS double counting terms
      energs%evsic=energs%evsic-SIC%alpha*evsic_tmp
   end if

   if (ipotmethod == 2 .or. ipotmethod==3) then
      nullify(pkernelSIC%kernel)
   end if

   call f_release_routine()

END SUBROUTINE LocalHamiltonianApplication


!> Routine which calculates the application of nonlocal projectors on the wavefunctions
!! Reduce the wavefunction in case it is needed
subroutine NonLocalHamiltonianApplication(iproc,at,npsidim_orbs,orbs,rxyz,&
     Lzd,nl,psi,hpsi,eproj_sum,proj_G,paw)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => NonLocalHamiltonianApplication
  use psp_projectors, only: PSPCODE_PAW,PSP_APPLY_SKIP
  implicit none
  integer, intent(in) :: iproc, npsidim_orbs
  type(atoms_data), intent(in) :: at
  type(orbitals_data),  intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nl
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
  real(gp), intent(out) :: eproj_sum
  !PAW variables:
  type(gaussian_basis),dimension(at%astruct%ntypes),optional,intent(in)::proj_G !projectors in gaussian basis (for PAW)
  type(paw_objects),optional,intent(inout)::paw
  !local variables
  logical :: newmethod
  character(len=*), parameter :: subname='NonLocalHamiltonianApplication' 
  logical :: dosome, overlap
  integer :: ikpt,istart_ck,ispsi_k,isorb,ieorb,nspinor,iorb,iat,nwarnings
  integer :: iproj,ispsi,istart_c,ilr,ilr_skip,mproj,iatype,ispinor
  real(wp) :: hp,eproj
  real(wp), dimension(:), allocatable :: scpr

  newmethod=.true.

  eproj_sum=0.0_gp
  !quick return if no orbitals on this task
  if (orbs%norbp == 0) then
     return
  end if

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')


  call f_routine(id=subname)

  !array of the scalar products
  scpr=f_malloc(orbs%norbp*orbs%nspinor,id='scpr')



  nwarnings=0

  if(any(at%npspcode == PSPCODE_PAW)) then  
     newmethod=.false.
     !initialize to zero in PAW case
     if(.not. present(paw) .or. .not. present(proj_G)) then
        stop 'NonLocalHamiltonianApplication: proj_G or paw are not present'
     end if
     if (nl%nprojel >0) call to_zero(nl%nprojel,nl%proj(1))
  end if
  !here the localisation region should be changed, temporary only for cubic approach

  !apply the projectors following the strategy (On-the-fly calculation or not)

  !apply the projectors  k-point of the processor
  !starting k-point
  ikpt=orbs%iokpt(1)
  istart_ck=1
  ispsi_k=1
  ispsi=1 !to initialize the value in case of no projectors
  loop_kpt: do

     call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
     !localisation regions loop
     loop_lr: do ilr=1,Lzd%nlr
        !do something only if at least one of the orbitals lives in the ilr
        dosome=.false.
        do iorb=isorb,ieorb
           !dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
           !dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr .and. lzd%doHamAppl(ilr))
           dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
           if (dosome) exit
        end do
        if (.not. dosome) cycle loop_lr

        if (nl%on_the_fly) then
           !first create a projector ,then apply it for everyone
           iproj=0
           do iat=1,at%astruct%nat
              iatype=at%astruct%iatype(iat)

              mproj=nl%pspd(iat)%mproj
              !no projector on this atom
              if(mproj == 0) cycle
              !projector not overlapping with the locreg
              if(nl%pspd(iat)%tolr(ilr)%strategy == PSP_APPLY_SKIP) cycle
              !check if the atom projector intersect with the given localisation region
              !this part can be moved at the place of the analysis between psp and lrs
              call check_overlap(Lzd%Llr(ilr), nl%pspd(iat)%plr, Lzd%Glr, overlap)
              if(.not. overlap) cycle

              ! Now create the projector
              istart_c=1
              if(any(at%npspcode == PSPCODE_PAW)) then
                 call atom_projector_paw(ikpt,iat,0,istart_c,iproj,&
                      nl%nprojel,&
                      Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),paw%rpaw(iatype),&
                      rxyz(1,iat),at,orbs,&
                      nl%pspd(iat)%plr,nl%proj,nwarnings,proj_G(iatype))
              else
                 call atom_projector(ikpt,iat,0,istart_c,iproj,&
                      nl%nprojel,&
                      Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),rxyz(1,iat),at,orbs,&
                      nl%pspd(iat)%plr,nl%proj,nwarnings)
              end if

              !apply the projector to all the orbitals belonging to the processor
              !this part can be factorized somewhere else
              if (mproj ==1 .and. all(orbs%kpts(:,ikpt) == 0.0_gp) .and. &
                   newmethod .and. Lzd%nlr == 1) then
                 hp=at%psppar(1,1,iatype) !it is supposed that the only projector is the i=1 case
                 ispsi=ispsi_k
                 istart_c=1
                 
                 call apply_oneproj_operator(nl%pspd(iat)%plr%wfd,nl%proj(istart_c),hp,&
                      (ieorb-isorb+1)*nspinor,Lzd%Llr(ilr)%wfd,psi(ispsi),hpsi(ispsi),scpr)
                 istart_c=istart_c+nl%pspd(iat)%plr%wfd%nvctr_c+7*nl%pspd(iat)%plr%wfd%nvctr_f
!                 call f_malloc_dump_status()
                 ispsi=ispsi+&
                      (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor*(ieorb-isorb+1)
                 do iorb=isorb,ieorb
                    do ispinor=1,nspinor
                       eproj=hp*(scpr((iorb-isorb)*nspinor+ispinor)**2)
                       eproj_sum=eproj_sum+orbs%kwgts(orbs%iokpt(iorb))*&
                            orbs%occup(iorb+orbs%isorb)*eproj
                    end do
                 end do
              else
                 ispsi=ispsi_k
                 do iorb=isorb,ieorb
                    if (orbs%inwhichlocreg(iorb+orbs%isorb) /= ilr) then
                       !increase ispsi to meet orbital index
                       ilr_skip=orbs%inwhichlocreg(iorb+orbs%isorb)
                       ispsi=ispsi+(Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                       cycle
                    end if
                    istart_c=1

                    call nl_psp_application()

                    !                print *,'iorb,iat,eproj',iorb+orbs%isorb,ispsi,iat,eproj_sum
                    ispsi=ispsi+&
                         (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
                 end do
              end if
           end do

           !for the moment, localization region method is not tested with
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
              do iat=1,at%astruct%nat
                 iatype=at%astruct%iatype(iat)
                 ! Check if atom has projectors, if not cycle
                 mproj=nl%pspd(iat)%mproj
                 if(mproj == 0) cycle
                 !projector not overlapping with the locreg
                 if(nl%pspd(iat)%tolr(ilr)%strategy == PSP_APPLY_SKIP) cycle

                 !check if the atom intersect with the given localisation region
                 call check_overlap(Lzd%Llr(ilr), nl%pspd(iat)%plr, Lzd%Glr, overlap)
                 if(.not. overlap) stop 'ERROR all atoms should be in global'

                 call nl_psp_application()

                 !print *,'iorb,iat,eproj',iorb+orbs%isorb,iat,eproj_sum
              end do
              !print *,'TOTALPSI',iorb+orbs%isorb,sum(psi(ispsi:&
              !    ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor-1)),&
              !     dot((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor,&
              !     hpsi(ispsi),1,hpsi(ispsi),1)
              ispsi=ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
           end do
           istart_ck=istart_c !TO BE CHANGED IN THIS ONCE-AND-FOR-ALL
        else
           ! COULD CHANGE THIS NOW !!
           call yaml_warning('Localization Regions not allowed in once-and-for-all')
           stop
        end if

     end do loop_lr

     !last k-point has been treated
     if (ieorb == orbs%norbp) exit loop_kpt

     ikpt=ikpt+1
     ispsi_k=ispsi

  end do loop_kpt

  if (.not. nl%on_the_fly) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
     if (istart_ck-1 /= nl%nprojel) then
        call yaml_warning('Incorrect once-and-for-all psp application')
        stop
     end if
  end if
  !for the moment it has to be removed. A number of components in orbital distribution should be defined
  !if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'

  !used on the on-the-fly projector creation
  if (nwarnings /= 0 .and. iproc == 0) then
     call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
     if (nwarnings /= 0) then
        call yaml_newline()
        call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
     end if
  end if

  call f_free(scpr)
  call f_release_routine()

  call timing(iproc,'ApplyProj     ','OF')
contains

!!$  !> count the number of projectors given a set of psppar
!!$  pure function nproj(psppar)
!!$    use module_base, only: gp
!!$    implicit none
!!$    real(gp), dimension(0:4,0:6), intent(in) :: psppar
!!$    integer :: nproj
!!$    !local variables
!!$    integer :: i,l
!!$    nproj = 0
!!$    !count over all the channels
!!$    do l=1,4
!!$       !loop over all the projectors of the channel
!!$       do i=1,3
!!$          !loop over all the components of the projector
!!$          if (psppar(l,i) /= 0.0_gp) nproj=nproj+2*l-1
!!$       end do
!!$    end do
!!$  end function nproj

  !>code factorization useful for routine restructuring
  subroutine nl_psp_application()
    implicit none
    !local variables
    integer :: ncplx_p,ncplx_w,n_w,nvctr_p
    real(gp), dimension(3,3,4) :: hij
    real(gp) :: eproj
!!$    type(nlpsp_to_wfd) :: tolr 
!!$    integer, dimension(:), allocatable :: nbsegs_cf,keyag_lin_cf
!!$    real(wp), dimension(:,:), allocatable :: psi_pack
!!$    real(wp), dimension(:,:,:), allocatable :: pdpsi,hpdpsi
!!$    real(wp), dimension(:,:,:,:), allocatable :: scpr

    if (newmethod) then

       if (all(orbs%kpts(:,ikpt) == 0.0_gp)) then
          ncplx_p=1
       else
          ncplx_p=2
       end if
       if (nspinor > 1) then !which means 2 or 4
          ncplx_w=2
          n_w=nspinor/2
       else
          ncplx_w=1
          n_w=1
       end if
       
       !extract hij parameters
       call hgh_hij_matrix(at%npspcode(iatype),at%psppar(0,0,iatype),hij)

!!$       !allocate temporary workspace
!!$       keyag_lin_cf=f_malloc(Lzd%Llr(ilr)%wfd%nseg_c+Lzd%Llr(ilr)%wfd%nseg_f,id='keyag_lin_cf')
!!$       nbsegs_cf=f_malloc0(nl%pspd(iat)%plr%wfd%nseg_c+nl%pspd(iat)%plr%wfd%nseg_f,id='nbsegs_cf')  
!!$       call nullify_nlpsp_to_wfd(tolr)
!!$       call init_tolr(tolr,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%plr%wfd,keyag_lin_cf,nbsegs_cf)
!!$       call f_free(keyag_lin_cf,nbsegs_cf)

!!$       !create other workspaces  
!!$       psi_pack=f_malloc0(&
!!$            (/nl%pspd(iat)%plr%wfd%nvctr_c+7*nl%pspd(iat)%plr%wfd%nvctr_f,n_w*ncplx_w/),id='psi_pack')
!!$       scpr=f_malloc((/ncplx_w,n_w,ncplx_p,mproj/),id='scpr')
!!$       pdpsi=f_malloc((/max(ncplx_w,ncplx_p),n_w,mproj/),id='pdpsi')
!!$       hpdpsi=f_malloc((/max(ncplx_w,ncplx_p),n_w,mproj/),id='hpdpsi')
       
       call NL_HGH_application(hij,&
            ncplx_p,mproj,nl%pspd(iat)%plr%wfd,nl%proj(istart_c),&
            ncplx_w,n_w,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%tolr(ilr),nl%wpack,nl%scpr,nl%cproj,nl%hcproj,&
            psi(ispsi),hpsi(ispsi),eproj)

!!$       !free workspaces
!!$       call f_free(psi_pack)
!!$       call f_free(scpr)     
!!$       call f_free(pdpsi)
!!$       call f_free(hpdpsi)
!!$       call deallocate_nlpsp_to_wfd(tolr)

       nvctr_p=nl%pspd(iat)%plr%wfd%nvctr_c+7*nl%pspd(iat)%plr%wfd%nvctr_f
       istart_c=istart_c+nvctr_p*ncplx_p*mproj

       eproj_sum=eproj_sum+&
            orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj
    else
       if(any(at%npspcode == PSPCODE_PAW)) then
          call apply_atproj_iorb_paw(iat,iorb,ispsi,istart_c,&
               nl%nprojel,&
               at,orbs,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%plr,&
               nl%proj,psi(ispsi),hpsi(ispsi),eproj_sum,&
               proj_G(iatype),paw)
       else
          call apply_atproj_iorb_new(iat,iorb,istart_c,&
               nl%nprojel,&
               at,orbs,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%plr,&
               nl%proj,psi(ispsi),hpsi(ispsi),eproj_sum)
       end if
    end if
  end subroutine nl_psp_application

END SUBROUTINE NonLocalHamiltonianApplication

!> routine which puts a barrier to ensure that both local and nonlocal hamiltonians have been applied
!! in the GPU case puts a barrier to end the overlapped Local and nonlocal applications
subroutine SynchronizeHamiltonianApplication(nproc,npsidim_orbs,orbs,Lzd,GPU,xc,&
     & hpsi,ekin_sum,epot_sum,eproj_sum,evsic,eexctX)
   use module_base
   use module_types
   use module_xc
   implicit none
   integer, intent(in) :: nproc,npsidim_orbs
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(GPU_pointers), intent(inout) :: GPU
   type(xc_info), intent(in) :: xc
   real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,evsic,eexctX
   real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
   !local variables
   character(len=*), parameter :: subname='SynchronizeHamiltonianApplication'
   logical :: exctX
   integer :: i_all,i_stat,iorb,ispsi,ilr
   real(gp), dimension(4) :: wrkallred

   call f_routine(id='SynchronizeHamiltonianApplication')

   if(GPU%OCLconv .or. GPUconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
      if (GPU%OCLconv) call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU)
      ispsi=1
      do iorb=1,orbs%norbp
         ilr=orbs%inWhichLocreg(orbs%isorb+iorb)
         call axpy((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor,&
              1.0_wp,GPU%hpsi_ASYNC(ispsi),1,hpsi(ispsi),1)
         ispsi=ispsi+&
             (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
      end do
      call f_free_ptr(GPU%hpsi_ASYNC)
   endif

   exctX = xc_exctXfac(xc) /= 0.0_gp

   !energies reduction
   if (nproc > 1) then
      wrkallred(1)=ekin_sum 
      wrkallred(2)=epot_sum 
      wrkallred(3)=eproj_sum
      wrkallred(4)=evsic

      call mpiallred(wrkallred(1),4,MPI_SUM,bigdft_mpi%mpi_comm)

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

   call f_release_routine()

END SUBROUTINE SynchronizeHamiltonianApplication


!> Build the potential in the whole box
!! Control also the generation of an orbital
subroutine full_local_potential(iproc,nproc,orbs,Lzd,iflag,dpbox,xc,potential,pot,comgp)
  !ndimpot,ndimgrid,nspin,&
  !   ndimrhopot,i3rho_add,orbs,&
  !   Lzd,iflag,ngatherarr,potential,pot,comgp)
   use module_base
   use module_types
   use module_xc
   use communications_base, only: p2pComms
   use communications, only: synchronize_onesided_communication
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc,iflag!,nspin,ndimpot,ndimgrid
   !integer, intent(in) :: ndimrhopot,i3rho_add
   type(orbitals_data),intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: Lzd
   type(denspot_distribution), intent(in) :: dpbox
   type(xc_info), intent(in) :: xc
   !integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
   real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), target :: potential !< Distributed potential. Might contain the density for the SIC treatments
   real(wp), dimension(:), pointer :: pot
   !type(p2pCommsGatherPot),intent(inout), optional:: comgp
   type(p2pComms),intent(inout), optional:: comgp
   !local variables
   character(len=*), parameter :: subname='full_local_potential'
   logical :: odp,newvalue !orbital dependent potential
   integer :: npot,ispot,ispotential,ispin,ierr,i_stat,i_all,ii,ilr,iorb,iorb2,nilr,ni1,ni2
   integer:: istl, ist, size_Lpot, i3s, i3e, i2s, i2e, i1s, i1e
   integer,dimension(:),allocatable:: ilrtable
   real(wp), dimension(:), pointer :: pot1
   
   call timing(iproc,'Pot_commun    ','ON')
   call f_routine(id='full_local_potential')

   odp = (xc_exctXfac(xc) /= 0.0_gp .or. (dpbox%i3rho_add /= 0 .and. orbs%norbp > 0))

   !############################################################################
   ! Build the potential on the whole simulation box
   ! NOTE: in the linear scaling case this should be done for a given localisation
   !       region this routine should then be modified or integrated in HamiltonianApplication
   ! WARNING : orbs%nspin and nspin are not the same !! Check if orbs%nspin should be replaced everywhere
   !#############################################################################
   if (iflag<2) then

      !determine the dimension of the potential array
      if (odp) then
         if (xc_exctXfac(xc) /= 0.0_gp) then
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
      if (dpbox%mpi_env%nproc > 1) then

         pot1 = f_malloc_ptr(npot+ndebug,id='pot1')
         ispot=1
         ispotential=1
         do ispin=1,orbs%nspin
            call MPI_ALLGATHERV(potential(ispotential),dpbox%ndimpot,&
                 &   mpidtypw,pot1(ispot),dpbox%ngatherarr(0,1),&
                 dpbox%ngatherarr(0,2),mpidtypw,dpbox%mpi_env%mpi_comm,ierr)
            ispot=ispot+dpbox%ndimgrid
            ispotential=ispotential+max(1,dpbox%ndimpot)
         end do
         !continue to copy the density after the potential if required
         if (dpbox%i3rho_add >0 .and. orbs%norbp > 0) then
            ispot=ispot+dpbox%i3rho_add-1
            do ispin=1,orbs%nspin
               call MPI_ALLGATHERV(potential(ispotential),dpbox%ndimpot,&
                    &   mpidtypw,pot1(ispot),dpbox%ngatherarr(0,1),&
                    dpbox%ngatherarr(0,2),mpidtypw,dpbox%mpi_env%mpi_comm,ierr)
               ispot=ispot+dpbox%ndimgrid
               ispotential=ispotential+max(1,dpbox%ndimpot)
            end do
         end if
      else
         if (odp) then
            pot1 = f_malloc_ptr(npot+ndebug,id='pot1')
            call vcopy(dpbox%ndimgrid*orbs%nspin,potential(1),1,pot1(1),1)
            if (dpbox%i3rho_add >0 .and. orbs%norbp > 0) then
               ispot=dpbox%ndimgrid*orbs%nspin+1
               call vcopy(dpbox%ndimgrid*orbs%nspin,potential(ispot+dpbox%i3rho_add),1,pot1(ispot),1)
            end if
         else
            pot1 => potential
         end if
      end if
   else
       !!if(.not.comgp%communication_complete) call gatherPotential(iproc, nproc, comgp)
       !!if(.not.comgp%communication_complete) call wait_p2p_communication(iproc, nproc, comgp)
       !!call wait_p2p_communication(iproc, nproc, comgp)
       call synchronize_onesided_communication(iproc, nproc, comgp)
   end if

   call timing(iproc,'Pot_commun    ','OF') 

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

   call timing(iproc,'Pot_after_comm','ON')
   
   if(Lzd%nlr > 1) then
      ilrtable = f_malloc(orbs%norbp,id='ilrtable')
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
      ilrtable = f_malloc(1,id='ilrtable')
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
      !       call vcopy(lzd%ndimpotisf,pot,1,pot,1) 
      ! This is due to the dynamic memory managment. The original version was: pot=>pot1
      !pot = f_malloc_ptr(npot+ndebug,id='pot')
      !pot=pot1
      !call f_free_ptr(pot1)
      pot=>pot1
   else if(iflag>0 .and. iflag<2) then
      pot = f_malloc_ptr(lzd%ndimpotisf+ndebug,id='pot')
      ! Cut potential
      istl=1
      do iorb=1,nilr
         ilr = ilrtable(iorb)
         ! Cut the potential into locreg pieces
         call global_to_local(Lzd%Glr,Lzd%Llr(ilr),orbs%nspin,npot,lzd%ndimpotisf,pot1,pot(istl))
         istl = istl + Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspin
      end do
   else
      if(.not.associated(pot)) then !otherwise this has been done already... Should be improved.
         pot = f_malloc_ptr(lzd%ndimpotisf+ndebug,id='pot')

         ist=1
         do iorb=1,nilr
            ilr = ilrtable(iorb)
            !determine the dimension of the potential array (copied from full_local_potential)
            if (xc_exctXfac(xc) /= 0.0_gp) then
               stop 'exctX not yet implemented!'
            else
               size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
            end if


            ! Extract the part of the potential which is needed for the current localization region.
            i3s=lzd%Llr(ilr)%nsi3-comgp%ise(5,iproc)+2 ! starting index of localized  potential with respect to total potential in comgp%recvBuf
            i3e=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i-comgp%ise(5,iproc)+1 ! ending index of localized potential with respect to total potential in comgp%recvBuf
            i2s=lzd%Llr(ilr)%nsi2-comgp%ise(3,iproc)+2
            i2e=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i-comgp%ise(3,iproc)+1
            i1s=lzd%Llr(ilr)%nsi1-comgp%ise(1,iproc)+2
            i1e=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i-comgp%ise(1,iproc)+1
            ni1=comgp%ise(2,iproc)-comgp%ise(1,iproc)+1
            ni2=comgp%ise(4,iproc)-comgp%ise(3,iproc)+1
            if(i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i) then
               write(*,'(a,i0,3x,i0)') 'ERROR: i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i',i3e-i3s+1, Lzd%Llr(ilr)%d%n3i
               stop
            end if
            call global_to_local_parallel(lzd%Glr, lzd%Llr(ilr), orbs%nspin, comgp%nrecvBuf, size_Lpot,&
                 comgp%recvBuf, pot(ist), i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2)

            ist = ist + size_lpot
         end do
      end if
   end if

   call f_free(ilrtable)

   ! Deallocate pot.
   if (iflag<2 .and. iflag>0) then
      if (dpbox%mpi_env%nproc > 1) then
         call f_free_ptr(pot1)
      else
         if (xc_exctXfac(xc) /= 0.0_gp) then
            call f_free_ptr(pot1)
         else
            nullify(pot1)
         end if
      end if
   end if

   call f_release_routine()
   call timing(iproc,'Pot_after_comm','OF')

   !!call mpi_finalize(ierr)
   !!stop


END SUBROUTINE full_local_potential


subroutine free_full_potential(nproc,flag,xc,pot,subname)
   use module_base
   use module_xc
   implicit none
   character(len=*), intent(in) :: subname
   integer, intent(in) :: nproc, flag
   type(xc_info), intent(in) :: xc
   real(wp), dimension(:), pointer :: pot
   !local variables
   logical :: odp
   integer :: i_all,i_stat

   odp = xc_exctXfac(xc) /= 0.0_gp
   if (nproc > 1 .or. odp .or. flag > 0 ) then
      !call f_free_ptr(pot)
      call f_free_ptr(pot)
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

  ! Gibbs Free Energy
  energs%energy=energs%eKS-energs%eTS+energs%ePV

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
     energs,wfn,gnrm,gnrm_zero,paw)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_energy_and_gradient
  use yaml_output
  use communications, only: transpose_v, untranspose_v
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,iscf,iter
  type(energy_terms), intent(inout) :: energs
  type(GPU_pointers), intent(in) :: GPU
  type(DFT_wavefunction), intent(inout) :: wfn
  real(gp), intent(out) :: gnrm,gnrm_zero
  type(paw_objects),optional,intent(inout)::paw
  !local variables
  character(len=*), parameter :: subname='calculate_energy_and_gradient' 
  logical :: lcs
  integer :: ikpt,iorb,i_all,i_stat,k
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

  !calculate orbital polarisation directions
  if(wfn%orbs%nspinor==4) then
     mom_vec = f_malloc_ptr((/ 4, wfn%orbs%norb, min(nproc, 2) /),id='mom_vec')

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
   call toglobal_and_transpose(iproc,nproc,wfn%orbs,wfn%Lzd,wfn%comms,wfn%hpsi,wfn%psi)
  
  !PAW:
  !transpose the spsi wavefunction
  if(present(paw)) then
     call toglobal_and_transpose(iproc,nproc,wfn%orbs,wfn%Lzd,wfn%comms,paw%spsi,wfn%psi)
  end if

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     wfn%psit => wfn%psi
     ! work array not used for nproc==1, so pass the same address
     call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%psit(1),wfn%psit(1))
  end if

!!$  if (iproc==0 .and. verbose > 0) then
!!$     call yaml_map('Symmetrize Lagr. Multiplier',wfn%SIC%alpha/=0.0_gp)
!!$  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  !takes also into account parallel k-points distribution
  !here the orthogonality with respect to other occupied functions should be 
  !passed as an optional argument
  energs%trH_prev=energs%trH
  if(present(paw)) then
    !PAW: spsi is used.
    call orthoconstraint(iproc,nproc,wfn%orbs,wfn%comms,wfn%SIC%alpha/=0.0_gp,wfn%psit,wfn%hpsi,energs%trH,paw%spsi) !n(m)
  else
    !NC:
    call orthoconstraint(iproc,nproc,wfn%orbs,wfn%comms,wfn%SIC%alpha/=0.0_gp,wfn%psit,wfn%hpsi,energs%trH) !n(m)
  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%hpsi(1),wfn%psi(1))

  if(present(paw)) then
   !retranspose the spsi wavefunction
   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,paw%spsi(1),wfn%psi(1))
  end if

  !deallocate temporary array
  !i_all=-product(shape(work))*kind(work)
  !deallocate(work,stat=i_stat)
  !call memocc(i_stat,i_all,'work',subname)



  !after having calcutated the trace of the hamiltonian, the functional have to be defined
  !new value without the trace, to be added in hpsitopsi
  if (iscf >1) then
     wfn%diis%energy=energs%trH
  else
     wfn%diis%energy=energs%eKS!trH-eh+exc-evxc-eexctX+eion+edisp(not correct for non-integer occnums)
  end if

  if (iproc==0 .and. verbose > 0) then
     call yaml_map('Orthoconstraint',.true.)
  end if


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
        call yaml_mapping_open('Energy inconsistencies')
           call yaml_map('Band Structure Energy',energs%ebs,fmt='(1pe22.14)')
           call yaml_map('Trace of the Hamiltonian',energs%trH,fmt='(1pe22.14)')
           call yaml_map('Relative inconsistency',tt,fmt='(1pe9.2)')
        call yaml_mapping_close()
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
  else if (GPU%OCLconv) then
     call preconditionall_OCL(wfn%orbs,wfn%Lzd%Glr,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),ncong,&
          wfn%hpsi,gnrm,gnrm_zero,GPU)
  else
     !this is the final routine, the confining potential has to be passed to 
     !switch between the global and delocalized preconditioner
     call preconditionall2(iproc,nproc,wfn%orbs,wfn%Lzd,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),&
          ncong,wfn%orbs%npsidim_orbs,wfn%hpsi,wfn%confdatarr,gnrm,gnrm_zero)
  end if

  call timing(iproc,'Precondition  ','OF')

  !sum over all the partial residues
  if (nproc > 1) then
      garray(1)=gnrm
      garray(2)=gnrm_zero
     call mpiallred(garray(1),2,MPI_SUM,bigdft_mpi%mpi_comm)
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

  !if (iproc==0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

  if (iproc==0  .and. verbose > 0) then
     call yaml_map('Preconditioning',.true.)
     call yaml_newline()
  end if

  if (wfn%orbs%nspinor == 4) then
     !only the root process has the correct array
     if(iproc==0 .and. verbose > 0) then
        call yaml_sequence_open('Magnetic polarization per orbital')
        call yaml_newline()
        !write(*,'(1x,a)')&
        !     &   'Magnetic polarization per orbital'
!!$        write(*,'(1x,a)')&
!!$             &   '  iorb    m_x       m_y       m_z'
        do iorb=1,wfn%orbs%norb
           call yaml_sequence(advance='no')
           call yaml_mapping_open(flow=.true.)
           call yaml_map('iorb',iorb,fmt='(i5)')
           call yaml_map('M',(/(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)/),fmt='(3f10.5)')
           call yaml_mapping_close()
           if (iorb < wfn%orbs%norb)call yaml_newline()
           !write(*,'(1x,i5,3f10.5)') &
           !     &   iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
        end do
        call yaml_sequence_close()
     end if
     call f_free_ptr(mom_vec)
  end if


  !write the energy information
  if (iproc == 0) then
     call write_energies(iter,iscf,energs,gnrm,gnrm_zero,' ')
  endif

END SUBROUTINE calculate_energy_and_gradient


!> Operations after h|psi> 
!! (transposition, orthonormalisation, inverse transposition)
subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn,&
   at,nlpsp,paw,rxyz,eproj_sum,proj_G)
   use module_base
   use module_types
   use module_interfaces, except_this_one_A => hpsitopsi
   use yaml_output
   use communications, only: transpose_v, untranspose_v
   use psp_projectors, only: PSPCODE_PAW
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc,idsx,iter
   type(DFT_wavefunction), intent(inout) :: wfn
   type(atoms_data), intent(in) :: at
   type(DFT_PSP_projectors), intent(inout) :: nlpsp 
   type(paw_objects),optional,intent(inout)::paw
   type(gaussian_basis),dimension(at%astruct%ntypes),optional,intent(in)::proj_G !projectors in gaussian basis (for PAW)
   real(gp),optional, intent(out) :: eproj_sum
   real(gp),optional, dimension(3,at%astruct%nat), intent(in) :: rxyz
   !local variables
   !character(len=*), parameter :: subname='hpsitopsi'
   integer :: i_all,i_stat
   !debug
   integer :: jorb,iat
   !end debug

   if(any(at%npspcode == PSPCODE_PAW)) then
     if( (.not. present(paw)) .or. &
&        (.not. present(proj_G)) .or. &
&        (.not. present(eproj_sum)) .or. &
&        (.not. present(rxyz))) then
         write(*,*)'ERROR: hpsitopsi for PAW needs the following optional variables::'
         write(*,*)'       paw, proj_G, eproj, rxyz'
         stop
     end if
   end if
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
   call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,&
        wfn%hpsi(1),wfn%psi(1))
   
   !!experimental, orthogonalize the preconditioned gradient wrt wavefunction
   !call orthon_virt_occup(iproc,nproc,orbs,orbs,comms,comms,psit,hpsi,(verbose > 2))

   !apply the minimization method (DIIS or steepest descent)
   call timing(iproc,'Diis          ','ON')

   call psimix(iproc,nproc,sum(wfn%comms%ncntt(0:nproc-1)),wfn%orbs,wfn%comms,wfn%diis,wfn%hpsi,wfn%psit)
  
   !Update spsi, since psi has changed
   !Pending: make this with the transposed wavefunctions:
   if(any(at%npspcode == PSPCODE_PAW)) then
     !retranspose psit
     call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        &   wfn%psit(1),wfn%hpsi(1),out_add=wfn%psi(1))

     !Calculate  hpsi,spsi and cprj with new psi
     if (wfn%orbs%npsidim_orbs >0) then 
       call to_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
       call to_zero(wfn%orbs%npsidim_orbs,paw%spsi(1))
     end if
     call NonLocalHamiltonianApplication(iproc,at,wfn%orbs%npsidim_orbs,wfn%orbs,rxyz,&
          wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,eproj_sum,proj_G,paw)

!    Transpose spsi:     
     call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,paw%spsi(1),wfn%hpsi(1))
!    Gather cprj:
     call gather_cprj()
   end if

   call timing(iproc,'Diis          ','OF')

   if (iproc == 0 .and. verbose > 1) then
      !write(*,'(1x,a)',advance='no')&
      !&   'Orthogonalization...'
      call yaml_map('Orthogonalization Method',wfn%orthpar%methortho,fmt='(i3)')
   end if

   if(any(at%npspcode == PSPCODE_PAW)) then
     call orthogonalize(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%orthpar,paw)
   else
     call orthogonalize(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%orthpar)
   end if

   !call checkortho_p(iproc,nproc,norb,nvctrp,psit)
   if(any(at%npspcode == PSPCODE_PAW)) call checkortho_paw(iproc,wfn%orbs%norb*wfn%orbs%nspinor,&
     wfn%comms%nvctr_par(iproc,0),wfn%psit,paw%spsi)
   !debug: 
   if(any(at%npspcode == PSPCODE_PAW)) then
     write(*,*)'hpsiortho, l1478 erase me:'
     do iat=1,paw%natom 
       do jorb=1,wfn%orbs%norbu
       write(*,'(a,2(i4,1x),1000f20.12)')' cprj(iat,jorb)%cp(:,:)=',iat,jorb,paw%cprj(iat,jorb)%cp(:,:)
       end do
     end do
   end if

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit(1),wfn%hpsi(1),out_add=wfn%psi(1))

   if (nproc == 1) then
      nullify(wfn%psit)
   end if

   !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   call diis_or_sd(iproc,idsx,wfn%orbs%nkptsp,wfn%diis)

   !previous value already filled
   wfn%diis%energy_old=wfn%diis%energy

   !
   !DEBUG hpsi
   !if(paw%usepaw==1 .and. nproc==1) then
   !   call debug_hpsi(wfn,at) 

   !   !
   !   !Recalculate  hpsi,spsi and cprj with new psi
   !   !PENDING: Here we only should recalculate cprj, hpsi and spsi are already updated.
   !   !         make a routine calculate_cprj based on NonLocalHam..?
   !   write(*,*)'Recalculate hpsi,spsi and cprj with new psi'
   !   if (wfn%orbs%npsidim_orbs >0) call to_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
   !   if (wfn%orbs%npsidim_orbs >0.and. paw%usepaw==1) call to_zero(wfn%orbs%npsidim_orbs,paw%spsi(1))
   !   call NonLocalHamiltonianApplication(iproc,at,wfn%orbs,rxyz,&
   !          proj,wfn%Lzd,nlpspd,wfn%psi,wfn%hpsi,eproj_sum,proj_G,paw)

   !   !check ortho again, to see if psi and spsi are correct
   !   if(paw%usepaw==1) call checkortho_paw(iproc,wfn%orbs%norb*wfn%orbs%nspinor,wfn%comms%nvctr_par(iproc,0),wfn%psi,paw%spsi)

   !end if
   !END DEBUG

contains

   !> In this routine cprj will be communicated to all processors
   subroutine gather_cprj()
     use module_base
     use module_types
     implicit none
   ! Local variables:
     integer::iatom,ilmn,iorb,ierr,ikpts,jproc
     character(len=*), parameter :: subname='gather_cprj'
   ! Tabulated data to be send/received for mpi
     integer,allocatable,dimension(:):: ndsplt
     integer,allocatable,dimension(:):: ncntd 
     integer,allocatable,dimension(:):: ncntt 
   ! auxiliar arrays
     real(dp),allocatable,dimension(:,:,:,:)::raux,raux2  

     if (nproc > 1) then

   !   Allocate temporary arrays
       ndsplt = f_malloc(0.to.nproc-1,id='ndsplt')
       ncntd = f_malloc(0.to.nproc-1,id='ncntd')
       ncntt = f_malloc(0.to.nproc-1,id='ncntt')
       raux = f_malloc((/ 2, paw%lmnmax, paw%natom, wfn%orbs%norbp /),id='raux')
       raux2 = f_malloc((/ 2, paw%lmnmax, paw%natom, wfn%orbs%norb /),id='raux2')

   !   Set tables for mpi operations:
   !   Send buffer:
       do jproc=0,nproc-1
         ncntd(jproc)=0
         do ikpts=1,wfn%orbs%nkpts
           ncntd(jproc)=ncntd(jproc)+&
                2*paw%lmnmax*paw%natom*wfn%orbs%norb_par(iproc,ikpts)*wfn%orbs%nspinor
         end do
       end do
   !   receive buffer:
       do jproc=0,nproc-1
          ncntt(jproc)=0
          do ikpts=1,wfn%orbs%nkpts
             ncntt(jproc)=ncntt(jproc)+&
                  2*paw%lmnmax*paw%natom*wfn%orbs%norb_par(jproc,ikpts)*wfn%orbs%nspinor
          end do
       end do
   !   Displacements table:
   !    ndspld(0)=0
   !    do jproc=1,nproc-1
   !      ndspld(jproc)=ndspld(jproc-1)+ncntd(jproc-1)
   !    end do
       ndsplt(0)=0
       do jproc=1,nproc-1
          ndsplt(jproc)=ndsplt(jproc-1)+ncntt(jproc-1)
       end do

   !   Transfer cprj to raux:
       raux=0.0_dp
   !Missing nspinor
       do iorb=1,wfn%orbs%norbp
         do iatom=1,paw%natom
           do ilmn=1,paw%cprj(iatom,iorb)%nlmn
             raux(:,ilmn,iatom,iorb)=paw%cprj(iatom,iorb)%cp(:,ilmn)
           end do
         end do
       end do
   !
   !    sendcnt=2*lmnmax*natom*orbs%norbp !N. of data to send
   !    recvcnt=2*lmnmax*natom*orbs%norb  !N. of data to receive
   !   
   !    call MPI_ALLGATHER(raux,sendcnt,MPI_DOUBLE_PRECISION,&
   !&     raux2,recvcnt,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHERV(raux,ncntd,mpidtypw,&
   &     raux2,ncntt,ndsplt,mpidtypw,MPI_COMM_WORLD,ierr)
   !
   !   Transfer back, raux2 to cprj:
   !   First set cprj to zero
       do iorb=1,wfn%orbs%norbp
         do iatom=1,paw%natom
           paw%cprj(iatom,iorb)%cp(:,:)=0.0_dp
         end do
       end do
   !
       do iorb=1,wfn%orbs%norb
         do iatom=1,paw%natom
           do ilmn=1,paw%cprj(iatom,iorb)%nlmn
             paw%cprj(iatom,iorb)%cp(:,ilmn)=raux2(:,ilmn,iatom,iorb)
           end do
         end do
       end do
   !   Deallocate arrays:
       call f_free(ndsplt)
       call f_free(ncntd)
       call f_free(ncntt)
       call f_free(raux)
       call f_free(raux2)
     end if


   end subroutine gather_cprj

END SUBROUTINE hpsitopsi



!>   First orthonormalisation
subroutine first_orthon(iproc,nproc,orbs,lzd,comms,psi,hpsi,psit,orthpar,paw)
   use module_base
   use module_types
   use module_interfaces, except_this_one_B => first_orthon
   use communications_base, only: comms_cubic
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: lzd
   type(comms_cubic), intent(in) :: comms
   type(orthon_data):: orthpar
   real(wp), dimension(:) , pointer :: psi,hpsi,psit
   type(paw_objects),optional,intent(inout)::paw
   !local variables
   character(len=*), parameter :: subname='first_orthon'
   integer :: i_stat,usepaw=0

   if(present(paw))usepaw=paw%usepaw
   !!!  if(nspin==4) then
   !!!     nspinor=4
   !!!  else
   !!!     nspinor=1
   !!!  end if

   if (nproc > 1) then
      !allocate hpsi array (used also as transposed)
      !allocated in the transposed way such as 
      !it can also be used as the transposed hpsi
      hpsi =f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug,id='hpsi')
      !allocate transposed principal wavefunction
      psit = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit')
   else
      psit => psi
   end if

   !to be substituted, must pass the wavefunction descriptors to the routine
   if (nproc>1) then
       call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),&
          &   hpsi(1),out_add=psit(1))
   else
       ! work array not nedded for nproc==1, so pass the same address
       call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),&
          &   psi(1),out_add=psit(1))
   end if

   if(usepaw==1) then
     call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar,paw)
   else
     call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar)
   end if

   !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

   if (nproc>1) then
       call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit(1),&
          &   hpsi(1),out_add=psi(1))
   else
       ! work array not nedded for nproc==1, so pass the same address
       call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit(1),&
          &   psit(1),out_add=psi(1))
   end if

   if (nproc == 1) then
      nullify(psit)
      !allocate hpsi array
      hpsi = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug,id='hpsi')
   end if

END SUBROUTINE first_orthon


!>   Transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
   use module_base
   use module_types
   use module_interfaces, fake_name => last_orthon
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc,iter
   real(wp), intent(out) :: evsum
   type(DFT_wavefunction), intent(inout) :: wfn
   logical, optional :: opt_keeppsit
   !local variables
   logical :: keeppsit
   character(len=*), parameter :: subname='last_orthon'
   integer :: i_all,i_stat

   if (present(opt_keeppsit)) then
      keeppsit=opt_keeppsit
   else
      keeppsit=.false.
   end if

   call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%hpsi(1),wfn%psi(1))
   if (nproc==1) then
      wfn%psit => wfn%psi
      ! workarray not used for nporc==1, so pass the sa,e address
      call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%psit(1),wfn%psit(1))
   end if

   call subspace_diagonalisation(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%hpsi,evsum)

   !here we should preserve hpsi and transpose it if we are in ensemble mimimization scheme

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit(1),wfn%hpsi(1),out_add=wfn%psi(1))
   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   if(.not.  keeppsit) then
      if (nproc > 1  ) then
         call f_free_ptr(wfn%psit)
      else
         nullify(wfn%psit)
      end if

      call f_free_ptr(wfn%hpsi)

   endif

   !call eigensystem_info(iproc,nproc,wfn%Lzd%Glr%wfd%nvctr_c+7*wfn%Lzd%Glr%wfd%nvctr_f,&
   !     wfn%orbs,wfn%psi)

END SUBROUTINE last_orthon

subroutine eigensystem_info(iproc,nproc,tolerance,nvctr,orbs,psi)
  use module_base
  use module_types
  use module_interfaces, except_this_one => eigensystem_info
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
     mom_vec = f_malloc_ptr((/ 4, orbs%norb, min(nproc, 2) /),id='mom_vec')

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          nvctr,&
          orbs%nspinor,psi,mom_vec)
  else
     nullify(mom_vec)   
  end if

  ! Send all eigenvalues to all procs.
  call broadcast_kpt_objects(nproc,orbs%nkpts,orbs%norb, &
       orbs%eval,orbs%ikptproc)

  !here the new occupation numbers should be recalculated for future needs

  !print the found eigenvalues
  if (iproc == 0) then
     call write_eigenvalues_data(tolerance,orbs,mom_vec)
  end if

  if (orbs%nspinor ==4) then
     call f_free_ptr(mom_vec)
  end if
end subroutine eigensystem_info


!> Finds the fermi level ef for an error function distribution with a width wf
!! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
subroutine evaltoocc(iproc,nproc,filewrite,wf0,orbs,occopt)
   use module_base
   use module_types
   use dictionaries, only: f_err_throw
   use yaml_output
   use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level
   implicit none
   logical, intent(in) :: filewrite
   integer, intent(in) :: iproc, nproc
   integer, intent(in) :: occopt      
   real(gp), intent(in) :: wf0   ! width of Fermi function, i.e. k*T
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   logical :: exitfermi
   real(gp), parameter :: pi=3.1415926535897932d0
   real(gp), parameter :: sqrtpi=sqrt(pi)
   real(gp), dimension(1,1,1) :: fakepsi
   integer :: ikpt,iorb,melec,ii,info_fermi
   real(gp) :: charge, chargef,wf
   real(gp) :: ef,electrons,dlectrons,factor,arg,argu,argd,corr,cutoffu,cutoffd,diff,full,res,resu,resd
   real(gp) :: a, x, xu, xd, f, df, tt
   !integer :: ierr
   type(fermi_aux) :: ft

   exitfermi=.false.
   !if (iproc.lt.1)  write(1000+iproc,*)  'ENTER Fermilevel',orbs%norbu,orbs%norbd,occopt

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
      if(iproc==0) print *, 'unrecognized occopt=', occopt
      stop 
   end select

   if (orbs%norbd==0) then 
      full=2.d0   ! maximum occupation for closed shell  orbital
   else
      full=1.d0   ! maximum occupation for spin polarized orbital
   endif

   if (orbs%nkpts.ne.1 .and. filewrite) then
      if (iproc == 0) print *,'Fermilevel: CANNOT write input.occ with more than one k-point'
      stop
   end if
   charge=0.0_gp
   do ikpt=1,orbs%nkpts
      !number of zero orbitals for the given k-point
      !overall charge of the system
      do iorb=1,orbs%norb
         charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb) * orbs%kwgts(ikpt)
      end do
   end do
   melec=nint(charge)
   !if (iproc == 0) write(1000+iproc,*) 'charge,wf',charge,melec,wf0
   call init_fermi_level(real(melec,gp)/full, 0.d0, ft, ef_interpol_det=1.d-12, verbosity=1)

   ! Send all eigenvalues to all procs (presumably not necessary)
   call broadcast_kpt_objects(nproc, orbs%nkpts, orbs%norb, &
      &   orbs%eval, orbs%ikptproc)
   
   if (wf0 > 0.0_gp) then
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

      ! electrons is N_electons = sum f_i * Wieght_i
      ! dlectrons is dN_electrons/dEf =dN_electrons/darg * darg/dEf= sum df_i/darg /(-wf) , darg/dEf=-1/wf
      !  f:= occupation # for band i ,  df:=df/darg
      wf=wf0
      loop_fermi: do ii=1,100
   !write(1000+iproc,*) 'iteration',ii,' -------------------------------- '
         factor=1.d0/(sqrt(pi)*wf)
         if (ii == 100 .and. iproc == 0) print *,'WARNING Fermilevel'
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
               !call yaml_map('arg,f,orbs%kwgts(ikpt)',(/arg,f,orbs%kwgts(ikpt)/))
               electrons=electrons+ f  * orbs%kwgts(ikpt)  ! electrons := N_e(Ef+corr.)
               dlectrons=dlectrons+ df * orbs%kwgts(ikpt)  ! delectrons:= dN_e/darg ( Well! later we need dN_e/dEf=-1/wf*dN_e/darg
               !if(iproc==0) write(1000,*) iorb,arg,   f , df,dlectrons
            enddo
         enddo
         !call yaml_map('ef',ef)
         !call yaml_map('electrons',electrons)

         dlectrons=dlectrons/(-wf)  ! df/dEf=df/darg * -1/wf
         diff=-real(melec,gp)/full+electrons
         !if (iproc.lt.1) write(1000+iproc,*) diff,full,melec,real(melec,gp)
!         if (iproc.lt.1) flush(1000+iproc)
         !if (iproc.lt.1) write(1000+iproc,*) diff,1.d-11*sqrt(electrons),wf
         !if (iproc.lt.1) flush(1000+iproc)
         !Exit criterion satiesfied, Nevertheles do one mor update of fermi level
         if (abs(diff) < 1.d-11*sqrt(electrons) .and. wf .eq. wf0 ) exitfermi=.true.     ! Assume noise grows as sqrt(electrons)

          corr=diff/abs(dlectrons) ! for case of no-monotonic func. abs is needed
          if (abs(corr).gt.wf) then   !for such a large correction the linear approximation is not any more valid
           if (corr > 0.d0) corr=1.d0*wf
           if (corr < 0.d0*wf) corr=-1.d0*wf
            if (ii.le.10) wf=2.d0*wf  ! speed up search of approximate Fermi level by using higher Temperature 
          else
            wf=max(wf0,.5d0*wf)
          endif
         ef=ef-corr  ! Ef=Ef_guess+corr.
         !if (iproc.lt.1) write(1000+iproc,'(i5,5(1pe17.8))') ii,electrons,ef,dlectrons,abs(dlectrons),corr
!         if (iproc.lt.1) flush(1000+iproc)
         !call determine_fermi_level(ft, electrons, ef,info_fermi)
         !if (info_fermi /= 0) then
         !   call f_err_throw('Difficulties in guessing the new Fermi energy, info='//trim(yaml_toa(info_fermi)),&
         !        err_name='BIGDFT_RUNTIME_ERROR')
         !end if
         !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr) !debug
         if (exitfermi) exit loop_fermi
      end do loop_fermi

      do ikpt=1,orbs%nkpts
         argu=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu)-ef)/wf0
         argd=(orbs%eval((ikpt-1)*orbs%norb+orbs%norbu+orbs%norbd)-ef)/wf0
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

      if ((cutoffu > 1.d-12 .or. cutoffd > 1.d-12) .and. iproc == 0) then
         call yaml_warning('Occupation numbers do not fill all available levels' // &
              ' lastu=' // trim(yaml_toa(cutoffu,fmt='(1pe8.1)')) // &
              ' lastd=' // trim(yaml_toa(cutoffd,fmt='(1pe8.1)')))
      end if
      !if (iproc.lt.1) write(1000+iproc,'(1x,a,1pe21.14,2(1x,e8.1))') 'Fermi level, Fermi distribution cut off at:  ',ef,cutoffu,cutoffd
!      if (iproc.lt.1) flush(1000+iproc)
      orbs%efermi=ef

      !update the occupation number
      do ikpt=1,orbs%nkpts
         do iorb=1,orbs%norbu + orbs%norbd
            arg=(orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf0
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
               orbs%eTS=orbs%eTS+full*wf0/(2._gp*sqrt(pi))*exp(-((orbs%eval((ikpt-1)*orbs%norb+iorb)-ef)/wf0)**2)
            else if (occopt == SMEARING_DIST_FERMI) then
               !Fermi function
               tt=orbs%occup((ikpt-1)*orbs%norb+iorb)
               orbs%eTS=orbs%eTS-full*wf0*(tt*log(tt) + (1._gp-tt)*log(1._gp-tt))
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
      if (abs(melec - chargef) > 1e-9)  then
      !if (abs(real(melec,gp)- chargef) > 1e-6)  then
         if (orbs%nspinor /= 4) call eigensystem_info(iproc,nproc,1.e-8_gp,0,orbs,fakepsi)
         call f_err_throw('Failed to determine correctly the occupation number, expected='//yaml_toa(charge)// &
              ', found='//yaml_toa(chargef),err_name='BIGDFT_RUNTIME_ERROR')
      end if
      !DEBUG call yaml_map('Electronic charges (expected, found)', (/ charge, chargef /))
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
   use module_types
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

      call to_zero(4*norb*ndim,mom_vec)

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
         norb_displ = f_malloc(0.to.nproc-1,id='norb_displ')

         norb_displ(0)=0
         do jproc=1,nproc-1
            norb_displ(jproc)=norb_displ(jproc-1)+norb_par(jproc-1)
         end do

         call MPI_GATHERV(mom_vec(1,1,2),4*norb_par(iproc),mpidtypw,&
            &   mom_vec(1,1,1),4*norb_par,4*norb_displ,mpidtypw,&
         0,bigdft_mpi%mpi_comm,ierr)

         call f_free(norb_displ)
      end if

   end if

END SUBROUTINE calc_moments


!> Check communications (transpose and untranspose orbitals)
subroutine check_communications(iproc,nproc,orbs,lzd,comms)
   use module_base
   use module_types
   use module_interfaces, except_this_one => check_communications
   use yaml_output
   use communications_base, only: comms_cubic
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: lzd
   type(comms_cubic), intent(in) :: comms
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
   psi = f_malloc(max(orbs%npsidim_orbs, orbs%npsidim_comp),id='psi')
   pwork = f_malloc_ptr(max(orbs%npsidim_orbs, orbs%npsidim_comp),id='pwork')

   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
         do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
            !vali=real(i,wp)/512.0_wp  ! *1.d-5
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            psi(i+indspin+indorb)=psival!(valorb+vali)*(-1)**(ispinor-1)
         end do
      end do
   end do

   !transpose the hpsi wavefunction
   call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),pwork(1))

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
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp) .and. maxdiff < 1.e-6_wp) then
      call yaml_warning('Transposition error of'//&
           trim(yaml_toa(maxdiff,fmt='(1pe15.7)'))//', check whether results are meaningful!')
   else if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
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
      !open(unit=22,file=trim(filename),status='unknown')
      call print_distribution_schemes(nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      !close(unit=22)
   end if

   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   if (abort) then
      if (iproc == 0) call print_distribution_schemes(nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
   end if

   !retranspose the hpsi wavefunction
   call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,&
      &   psi(1),pwork(1))

   maxdiff=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
         do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
            !vali=real(i,wp)/512.d0  !*1.d-5
            !psival=(valorb+vali)*(-1)**(ispinor-1)
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            maxdiff=max(abs(psi(i+indspin+indorb)-psival),maxdiff)
         end do
      end do
   end do

   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp) .and. maxdiff < 1.e-6_wp) then
      call yaml_warning('Inverse transposition error of'//&
           trim(yaml_toa(maxdiff,fmt='(1pe15.7)'))//', check whether results are meaningful!')
   else if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
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
         indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
         do ispinor=1,orbs%nspinor
            indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
            do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
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

   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

   call f_free(psi)
   call f_free_ptr(pwork)

END SUBROUTINE check_communications


!> define a value for the wavefunction which is dependent of the indices
subroutine test_value(ikpt,iorb,ispinor,icomp,val)
   use module_base
   implicit none
   integer, intent(in) :: ikpt,icomp,iorb,ispinor
   real(wp), intent(out) :: val
   !local variables
   integer, parameter :: ilog=6
   real(wp) :: valkpt,valorb,vali

   ! recognizable pattern, for debugging
    valkpt=real(10**ilog*(ikpt-1),wp)!real(512*ikpt,wp)
    valorb=real(iorb,wp)+valkpt
    vali=real(icomp,wp)*10.0_wp**(-ilog)  !real(icomp,wp)/512.0_wp  ! *1.d-5
   !
   ! val=(valorb+vali)*(-1)**(ispinor-1)

   !valkpt=real(512*ikpt,wp)
   !valorb=real(iorb,wp)+valkpt
   !vali=real(icomp,wp)/512.0_wp  ! *1.d-5

   val=(valorb+vali)*(-1)**(ispinor-1)

END SUBROUTINE test_value

!> Determine the components which were not communicated correctly
!! works only with the recognizable pattern of test function
subroutine wrong_components(psival,ikpt,iorb,icomp)
   use module_base
  implicit none
  real(wp), intent(in) :: psival
  integer, intent(out) :: ikpt,iorb,icomp
  integer, parameter :: ilog=6

  icomp=nint((psival-real(floor(psival),wp))*10.0_wp**ilog)
  ikpt=floor(psival)/(10**ilog)
  iorb=floor(psival)-(ikpt-1)*(10**ilog)

end subroutine wrong_components


subroutine debug_hpsi(wfn,at)!,proj_G,paw,nlpspd)
   use module_base
   use module_types
   implicit none
   !Arguments
   !type(nonlocal_psp_descriptors), intent(in) :: nlpspd
   type(DFT_wavefunction), intent(in) :: wfn
   !type(paw_objects),intent(in) :: paw
   type(atoms_data), intent(in) :: at
   !type(gaussian_basis),dimension(at%astruct%ntypes),intent(in)::proj_G !projectors in gaussian basis (for PAW)
   !Local variables   
   integer :: ispsi_a,ispsi_b,ia,ib,iatype,ispinor
   integer :: nvctr_c,nvctr_f,nvctr_tot
   integer :: ncplx,nspinor
   !integer :: i,istart_j,jlmn,j_shell,j_l,j_m
   !integer :: mbseg_c,mbseg_f,mbvctr_c,mbvctr_f,istart_ck
   integer :: ikpt,ispsi_k,iat,ilr,ilr_skip,iorb
   integer :: isorb,ieorb,ispsi
   logical :: dosome
   real(gp),dimension(2) :: scpr
   real(gp) :: scalprod(2,2)
   real(dp) :: ddot
   
   write(*,*)'Calculate <PSI|H|PSI> for debugging'
   
   nvctr_c=wfn%Lzd%glr%wfd%nvctr_c
   nvctr_f=wfn%Lzd%glr%wfd%nvctr_f
   nvctr_tot=nvctr_c+7*nvctr_f
   
   ikpt=wfn%orbs%iokpt(1)
   
   ispsi_k=1
   loop_kpt: do

      call orbs_in_kpt(ikpt,wfn%orbs,isorb,ieorb,nspinor)
      
      loop_lr: do ilr=1,wfn%Lzd%nlr
         !do something only if at least one of the orbitals lives in the ilr
         dosome=.false.
         do iorb=isorb,ieorb
            dosome = (wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb) == ilr)
            if (dosome) exit
         end do
         if (.not. dosome) cycle loop_lr
      
      do iat=1,at%astruct%nat
           iatype=at%astruct%iatype(iat)
              ispsi=ispsi_k
            do iorb=isorb,ieorb
               if (wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb) /= ilr) then
                  !increase ispsi to meet orbital index
                  ilr_skip=wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb)
                  ispsi=ispsi+(wfn%Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                  cycle
               end if
      
               !call plr_segs_and_vctrs(nlpspd%plr(iat),mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
               call ncplx_kpt(wfn%orbs%iokpt(iorb),wfn%orbs,ncplx)
               !loop over all the components of the wavefunction
               do ispinor=1,nspinor,ncplx
                 ispsi_a=ispsi
                 ispsi_b=ispsi
                 do ia=1,ncplx
                    do ib=1,ncplx
                       scalprod(ia,ib)=ddot(nvctr_tot,wfn%psi(ispsi_a),1,wfn%hpsi(ispsi_b),1)
                       ispsi_b=ispsi_b+nvctr_tot
                    end do
                    ispsi_a=ispsi_a+nvctr_tot
                 end do
                 !then define the result
                 if (ncplx == 1) then
                    scpr(1)=scalprod(1,1)
                 else if (ncplx == 2) then
                    scpr(1)=scalprod(1,1)+scalprod(2,2)
                    scpr(2)=scalprod(1,2)-scalprod(2,1)
                 end if
              write(*,'("<psi|H|psi> for kpt=,",i5," iat=",i5,"=>",2(f15.5,1x))')ikpt,iat,scpr
                 write(*,*)'ispsi_a','ispsi_b','nvctr_tot',ispsi_a,ispsi_b,nvctr_tot
               end do !ispinor
            end do !iorb
         end do !iat
      end do loop_lr

      !last k-point has been treated
      if (ieorb == wfn%orbs%norbp) exit loop_kpt
   
      ikpt=ikpt+1
      ispsi_k=ispsi

   end do  loop_kpt !ikpt

   !if (.not. DistProjApply) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
   !   if (istart_ck-1 /= nlpspd%nprojel) &
   !      &   stop 'incorrect once-and-for-all psp application'
   !end if

END SUBROUTINE debug_hpsi


subroutine broadcast_kpt_objects(nproc, nkpts, ndata, data, ikptproc)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: nproc, nkpts, ndata
   integer, dimension(nkpts), intent(in) :: ikptproc
   real(gp), dimension(ndata,nkpts), intent(inout) :: data

   integer :: ikpt, ierr

   if (nproc > 1) then
      do ikpt = 1, nkpts
         call MPI_BCAST(data(1,ikpt), ndata,mpidtypg, &
            &   ikptproc(ikpt), bigdft_mpi%mpi_comm, ierr)
         !redundant barrier 
         call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
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
!!  type(comms_cubic), intent(in) :: comms
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
!!  call to_zero(ndimovrlp(nspin,orbs%nkpts),alag)
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
!!     call mpiallred(alag(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
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
!!  call vcopy(size(psi), psiwork(1), 1, psi(1), 1)
!!
!!  ! I think this is not required..
!!  call orthogonalize(iproc,nproc,orbs,comms,psi,orthpar)
!!
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, psi, work=psiwork)
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, hpsi, work=psiwork)
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


subroutine integral_equation(iproc,nproc,atoms,wfn,ngatherarr,local_potential,GPU,xc,nlpsp,rxyz)
  use module_base
  use module_types
  use module_xc
  use module_interfaces, fake_name => integral_equation
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(in) :: atoms
  type(DFT_wavefunction), intent(in) :: wfn
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(xc_info), intent(in) :: xc
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(:), pointer :: local_potential
  !local variables
  integer :: iorb,nbox,ilr,ist
  real(gp) :: eh_fake,eks
  type(energy_terms) :: energs_tmp
  type(coulomb_operator) :: G_Helmholtz
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: vpsi,vpsir

  call f_routine(id='helmholtz_equation')

  !first, apply the potential operator to the wavefunction
  vpsi=f_malloc0(wfn%orbs%npsidim_orbs,id='vpsi')

  call LocalHamiltonianApplication(iproc,nproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
       wfn%Lzd,wfn%confdatarr,ngatherarr,local_potential,wfn%psi,vpsi,&
       energs_tmp,wfn%SIC,GPU,2,xc) !potential only

  call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,rxyz,&
       wfn%Lzd,nlpsp,wfn%psi,vpsi,energs_tmp%eproj)

  !now vpsi is a wavefunction array in orbitals parallelization scheme which is associated to Vpsi
  !rescale it to match with the Green's function treatment
  call vscal(wfn%orbs%npsidim_orbs,-0.5_gp/pi_param,vpsi(1),1)

  !helmholtz-based preconditioning
  ilr=1 !for the moment only cubic version
  call initialize_work_arrays_sumrho(wfn%Lzd%Llr(ilr),w)
  !box elements size
  nbox=wfn%Lzd%Llr(ilr)%d%n1i*wfn%Lzd%Llr(ilr)%d%n2i*wfn%Lzd%Llr(ilr)%d%n3i

  ! Wavefunction in real space
  vpsir=f_malloc0(nbox,id='vpsir')

  !case for helmholtz-based preconditioning
  ist=0
  do iorb=1,wfn%orbs%norbp*wfn%orbs%nspinor
     ilr = wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb)
     !entering
     eks=wfn%orbs%eval(iorb+wfn%orbs%isorb)
     print *,'iorb,initial',iorb,nrm2(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,wfn%psi(1+ist),1)

!!$     call plot_wf('psi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
!!$          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),rxyz,wfn%psi(1+ist:))

!     call axpy(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,-eks,wfn%psi(1+ist),1,vpsi(1+ist),1)

     call plot_wf('Vpsi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),rxyz,vpsi(1+ist))


     call daub_to_isf(wfn%Lzd%Llr(ilr),w,vpsi(1+ist),vpsir(1))

     !sequential kernel


     G_Helmholtz=pkernel_init(.false.,0,1,0,wfn%Lzd%Llr(ilr)%geocode,&
          (/wfn%Lzd%Llr(ilr)%d%n1i,wfn%Lzd%Llr(ilr)%d%n2i,wfn%Lzd%Llr(ilr)%d%n3i/),&
          0.5_gp*wfn%Lzd%hgrids,16,mu0_screening=sqrt(2.0_gp*abs(eks)))

     call pkernel_set(G_Helmholtz,.true.)

     !apply it to the gradient to smooth it
     call H_potential('D',G_Helmholtz,vpsir(1),vpsir(1),eh_fake,0.d0,.false.)

     !convert the gradient back to the locreg
     call isf_to_daub(wfn%Lzd%Llr(ilr),w,vpsir(1),vpsi(1+ist))
!     call vscal(ncplx*(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f),1.0_gp/(16.0_gp*atan(1.0_gp)),hpsi(1+ist),1)
     print *,'iorb,gradient',iorb,nrm2(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,vpsi(1+ist),1),eh_fake
!!$     call plot_wf('GVpsi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
!!$          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),rxyz,vpsi(1+ist:))


     ist=ist+(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f)

     call pkernel_free(G_Helmholtz,'helmholtz_equation')

  end do
  call deallocate_work_arrays_sumrho(w)

  call f_free(vpsi)
  call f_free(vpsir)

  call f_release_routine()

end subroutine integral_equation

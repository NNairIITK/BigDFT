!> @file
!! Define the module module_interfaces containing all interfaces
!!
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG,DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Modules which contains all interfaces
module module_interfaces

   implicit none

   interface

      subroutine copy_old_wavefunctions(nproc,orbs,psi,&
            &   wfd_old,psi_old)
        use module_defs, only: wp
        use module_types
         implicit none
         integer, intent(in) :: nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd_old
         real(wp), dimension(:), pointer :: psi,psi_old
      END SUBROUTINE copy_old_wavefunctions


      subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor, &
                 nkpt,kpt,wkpt,orbs,linear_partition,basedist,basedistu,basedistd)
         use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: linear_partition !< repartition mode for the linear scaling version
         integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
         integer, intent(in) :: nspinor
         type(orbitals_data), intent(inout) :: orbs
         real(gp), dimension(nkpt), intent(in) :: wkpt
         real(gp), dimension(3,nkpt), intent(in) :: kpt
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedist
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedistu
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedistd
      END SUBROUTINE orbitals_descriptors


      subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,&
            &   crmult,frmult,calculate_bounds,Glr,output_denspot)
        use module_defs, only: gp
        use module_types
         implicit none
         !Arguments
         type(atoms_data), intent(in) :: atoms
         integer, intent(in) :: iproc
         real(gp), intent(in) :: hx,hy,hz,crmult,frmult
         real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
         !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
         logical,intent(in) :: calculate_bounds
         type(locreg_descriptors), intent(inout) :: Glr
         logical, intent(in), optional :: output_denspot
      END SUBROUTINE createWavefunctionsDescriptors

      subroutine createProjectorsArrays(lr,rxyz,at,orbs,&
           cpmult,fpmult,hx,hy,hz,dry_run,nlpsp,&
           init_projectors_completely_)
        !n(c) use module_base
        use module_types
        implicit none
        type(atoms_data), intent(in) :: at
        type(orbitals_data), intent(in) :: orbs
        real(kind=8), intent(in) :: cpmult,fpmult,hx,hy,hz
        type(locreg_descriptors),intent(in) :: lr
        real(kind=8), dimension(3,at%astruct%nat), intent(in) :: rxyz
        !real(kind=8), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
        logical, intent(in) :: dry_run
        type(DFT_PSP_projectors), intent(out) :: nlpsp
        logical,intent(in),optional :: init_projectors_completely_
      END SUBROUTINE createProjectorsArrays


       subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield,&
            & rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,n1,n2,n3,&
            & pot_ion,pkernel,psoffset)
         use module_defs, only: gp,dp
         use module_types
         implicit none
         type(denspot_distribution), intent(in) :: dpbox
         type(atoms_data), intent(in) :: at
         integer, intent(in) :: iproc,nproc,n1,n2,n3,dispersion
         real(gp), dimension(3), intent(in) :: elecfield
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         type(coulomb_operator), intent(inout) :: pkernel
         real(gp), intent(out) :: eion,edisp,psoffset
         real(dp), dimension(6),intent(out) :: ewaldstr
         real(gp), dimension(:,:), pointer :: fion,fdisp
         real(dp), dimension(*), intent(out) :: pot_ion
       END SUBROUTINE IonicEnergyandForces


       subroutine input_wf_empty(iproc, nproc, psi, hpsi, psit, orbs, &
            & band_structure_filename, input_spin, atoms, d, denspot)
         use module_defs, only: wp
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc
         type(orbitals_data), intent(in) :: orbs
         character(len = *), intent(in) :: band_structure_filename
         integer, intent(in) :: input_spin
         type(atoms_data), intent(in) :: atoms
         type(grid_dimensions), intent(in) :: d
         type(DFT_local_fields), intent(inout) :: denspot
         real(wp), dimension(:), pointer :: psi
         real(kind=8), dimension(:), pointer :: hpsi, psit
       END SUBROUTINE input_wf_empty

       subroutine input_wf_random(psi, orbs)
         use module_defs, only: wp
         use module_types
         implicit none
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_random

       subroutine input_wf_cp2k(iproc, nproc, nspin, atoms, rxyz, Lzd, &
            & psi, orbs)
         use module_defs, only: wp,gp
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, nspin
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
         type(local_zone_descriptors), intent(in) :: Lzd
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_cp2k

       subroutine input_wf_memory(iproc, atoms, &
            & rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old, &
            & rxyz, lzd, psi, orbs)
         use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz, rxyz_old
         real(gp), intent(in) :: hx_old, hy_old, hz_old
         type(local_zone_descriptors), intent(in) :: lzd
         type(grid_dimensions), intent(in) :: d_old
         type(wavefunctions_descriptors), intent(in) :: wfd_old
         type(orbitals_data), intent(in) :: orbs
         real(wp), dimension(:), pointer :: psi, psi_old
       END SUBROUTINE input_wf_memory

       subroutine input_wf_disk(iproc, nproc, input_wf_format, d, hx, hy, hz, &
            in, atoms, rxyz, wfd, orbs, psi)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, input_wf_format
         type(grid_dimensions), intent(in) :: d
         real(gp), intent(in) :: hx, hy, hz
         type(input_variables), intent(in) :: in
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
         !real(gp), dimension(3, atoms%astruct%nat), intent(out) :: rxyz_old
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_disk

       subroutine input_wf_diag(iproc,nproc,at,denspot,&
            orbs,nvirt,comms,Lzd,energs,rxyz,&
            nlpsp,ixc,psi,hpsi,psit,G,&
            nspin,GPU,input,onlywf)!,paw)
         ! Input wavefunctions are found by a diagonalization in a minimal basis set
         ! Each processors write its initial wavefunctions into the wavefunction file
         ! The files are then read by readwave
         ! @todo pass GPU to be a local variable of this routine (initialized and freed here)
         use module_defs, only: gp,dp,wp
         use module_types
         use gaussians
         use communications_base, only: comms_cubic
         implicit none
         !Arguments
         integer, intent(in) :: iproc,nproc,ixc
         integer, intent(inout) :: nspin,nvirt
         logical, intent(in) :: onlywf  !if .true. finds only the WaveFunctions and return
         type(atoms_data), intent(in) :: at
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(comms_cubic), intent(in) :: comms
         type(orbitals_data), intent(inout) :: orbs
         type(energy_terms), intent(inout) :: energs
         type(DFT_local_fields), intent(inout) :: denspot
         type(GPU_pointers), intent(in) :: GPU
         type(input_variables), intent(in) :: input
         !type(symmetry_data), intent(in) :: symObj
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G !basis for davidson IG
         real(wp), dimension(:), pointer :: psi,hpsi,psit
         !type(paw_objects),optional,intent(inout)::paw
       end subroutine input_wf_diag

       subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,&
            denspot,denspot0,nlpsp,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,&
            lzd_old,psi_old,rxyz_old,tmb_old,ref_frags,cdft,&
            locregcenters)
         use module_defs, only: gp,wp
         use f_enums, only: f_enumerator
         use module_types
         use module_fragments
         use constrained_dft
         implicit none
         integer, intent(in) :: iproc, nproc, input_wf_format
         type(f_enumerator), intent(in) :: inputpsi
         type(input_variables), intent(in) :: in
         type(GPU_pointers), intent(inout) :: GPU
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), target, intent(in) :: rxyz
         type(DFT_local_fields), intent(inout) :: denspot
         type(DFT_wavefunction), intent(inout) :: KSwfn,tmb,tmb_old !<input wavefunction
         type(energy_terms), intent(inout) :: energs !<energies of the system
         real(gp), dimension(*), intent(out) :: denspot0 !< Initial density / potential, if needed
         real(wp), dimension(:), pointer :: psi_old
         integer, intent(out) :: norbv
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz_old
         type(local_zone_descriptors),intent(in):: lzd_old
         type(system_fragment), dimension(:), pointer :: ref_frags
         type(cdft_data), intent(out) :: cdft
         real(kind=8),dimension(3,atoms%astruct%nat),intent(in),optional :: locregcenters
       END SUBROUTINE input_wf

      subroutine first_orthon(iproc,nproc,orbs,lzd,comms,psi,hpsi,psit,orthpar,paw)
        use module_defs, only: wp
         use module_types
         use communications_base, only: comms_cubic
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors),intent(in) :: lzd
         type(comms_cubic), intent(in) :: comms
         type(orthon_data):: orthpar
         type(paw_objects),optional,intent(inout)::paw
         real(wp), dimension(:) , pointer :: psi,hpsi,psit
      END SUBROUTINE first_orthon

      subroutine density_and_hpot(dpbox,symObj,orbs,Lzd,pkernel,rhodsc,GPU,xc,psi,rho,vh,hstrten)
        use module_defs, only: gp,wp,dp
        use module_types
        use module_atoms, only: symmetry_data
        use module_xc
        implicit none
        type(denspot_distribution), intent(in) :: dpbox
        type(rho_descriptors),intent(inout) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
        type(coulomb_operator), intent(inout) :: pkernel
        type(xc_info), intent(in) :: xc
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        type(GPU_pointers), intent(inout) :: GPU
        real(gp), dimension(6), intent(out) :: hstrten
        real(dp), dimension(:), pointer :: rho,vh
      end subroutine density_and_hpot

      subroutine sumrho(dpbox,orbs,Lzd,GPU,symObj,rhodsc,xc,psi,rho_p,mapping)
        use module_defs, only: wp,gp,dp
        use module_atoms, only: symmetry_data
        use module_types
        use module_xc
        implicit none
        !Arguments
        type(denspot_distribution), intent(in) :: dpbox
        type(rho_descriptors),intent(in) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
        type(xc_info), intent(in) :: xc
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        real(dp), dimension(:,:), pointer :: rho_p
        type(GPU_pointers), intent(inout) :: GPU
        integer,dimension(orbs%norb),intent(in),optional:: mapping
      END SUBROUTINE sumrho

      !starting point for the communication routine of the density
      subroutine communicate_density(dpbox,nspin,rhodsc,rho_p,rho,keep_rhop)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        logical, intent(in) :: keep_rhop !< preserves the total density in the rho_p array
        integer, intent(in) :: nspin
        type(rho_descriptors),intent(in) :: rhodsc
        type(denspot_distribution), intent(in) :: dpbox
        real(dp), dimension(:,:), pointer :: rho_p !< partial density in orbital distribution scheme
        real(dp), dimension(max(dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d,1),nspin), intent(out) :: rho
      END SUBROUTINE communicate_density


       subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
            Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
            energs,SIC,GPU,PotOrKin,xc,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
         use module_defs, only: gp,dp,wp
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: PotOrKin !< if true, only the potential operator is applied
         integer, intent(in) :: iproc,nproc,npsidim_orbs
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(SIC_data), intent(in) :: SIC
         type(xc_info), intent(in) :: xc
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
         type(confpot_data), dimension(orbs%norbp) :: confdatarr
         real(wp), dimension(:), pointer :: pot
         !real(wp), dimension(*) :: pot
         type(energy_terms), intent(inout) :: energs
         real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout) :: hpsi
         type(GPU_pointers), intent(inout) :: GPU
         type(coulomb_operator), intent(in), optional :: pkernel
         type(orbitals_data), intent(in), optional :: orbsocc
         real(wp), dimension(:), pointer, optional :: psirocc
         type(denspot_distribution),intent(in),optional :: dpbox
         real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
         type(p2pComms),intent(inout), optional:: comgp
         real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout),optional :: hpsi_noconf
         real(gp),intent(out),optional :: econf
       end subroutine LocalHamiltonianApplication


      subroutine SynchronizeHamiltonianApplication(nproc,npsidim_orbs,orbs,Lzd,GPU,xc,hpsi,&
           energs,energs_work)
        use module_defs, only: gp,wp
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: nproc,npsidim_orbs
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
        type(energy_terms), intent(inout) :: energs
        !real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        type(work_mpiaccumulate),optional,intent(inout) :: energs_work
      END SUBROUTINE SynchronizeHamiltonianApplication

      subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn,&
           at,nlpsp,eproj_sum)
        use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,idsx,iter
         type(DFT_wavefunction), intent(inout) :: wfn
         type(atoms_data), intent(in) :: at
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         real(gp),optional, intent(out) :: eproj_sum
      END SUBROUTINE hpsitopsi

      subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,iter
        real(wp), intent(out) :: evsum
        type(DFT_wavefunction), intent(inout) :: wfn
        logical, optional :: opt_keeppsit
      END SUBROUTINE last_orthon


      subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
           Glr,nlpsp,ncongt,pot,hgrid,rxyz,crmult,frmult,nspin,&
           psi,output_denspot,ekin_sum,epot_sum,eproj_sum,paw)
        use module_defs, only: gp,wp,dp
         use module_types
         use gaussians, only: gaussian_basis
         implicit none
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(locreg_descriptors), intent(in) :: Glr
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         integer, intent(in) :: iproc,nproc,ncongt,nspin
         logical, intent(in) :: output_denspot
         real(kind=8), dimension(3), intent(in) :: hgrid
         real(kind=8), intent(in) :: crmult,frmult,rbuf
         !real(kind=8), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
         real(kind=8), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
         real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
         real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
         type(paw_objects),optional,intent(inout)::paw
      END SUBROUTINE CalculateTailCorrection

      !added for abinit compatibility
      subroutine reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,&
           n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: n1_old,n2_old,n3_old,n1,n2,n3
         real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz_old,rxyz
         real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
         real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
      END SUBROUTINE reformatonewave
      subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
            &   hx,hy,hz,at,wfd,rxyz_old,rxyz,psi,eval,psifscf)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         logical, intent(in) :: useFormattedInput
         integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(atoms_data), intent(in) :: at
         real(gp), intent(in) :: hx,hy,hz
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), intent(out) :: eval
         real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
         real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
      END SUBROUTINE readonewave
      subroutine writeonewave(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,nat,rxyz,  &
         nseg_c,nvctr_c,keyg_c,keyv_c,  &
         nseg_f,nvctr_f,keyg_f,keyv_f, &
         psi_c,psi_f,eval)
         use module_defs, only: gp,dp,wp
         implicit none
         logical, intent(in) :: useFormattedOutput
         integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f
         real(gp), intent(in) :: hx,hy,hz
         real(wp), intent(in) :: eval
         integer, dimension(nseg_c), intent(in) :: keyv_c
         integer, dimension(nseg_f), intent(in) :: keyv_f
         integer, dimension(2,nseg_c), intent(in) :: keyg_c
         integer, dimension(2,nseg_f), intent(in) :: keyg_f
         real(wp), dimension(nvctr_c), intent(in) :: psi_c
         real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
         real(gp), dimension(3,nat), intent(in) :: rxyz
      END SUBROUTINE writeonewave

      subroutine davidson(iproc,nproc,in,at,&
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           rxyz,rhopot,nlpsp,pkernel,psi,v,dpbox,xc,GPU)
        use module_defs, only: gp,dp,wp
        use module_types
        use communications_base, only: comms_cubic
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc
        integer, intent(in) :: nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        type(local_zone_descriptors), intent(inout) :: Lzd
        type(orbitals_data), intent(inout) :: orbs
        type(comms_cubic), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpbox
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        type(coulomb_operator), intent(inout) :: pkernel
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc)
      end subroutine davidson

      subroutine preconditionall2(iproc,nproc,orbs,Lzd,hx,hy,hz,ncong,npsidim,hpsi,confdatarr,gnrm,gnrm_zero, &
                 linear_precond_convol_workarrays, linear_precond_workarrays)
        use module_defs, only: gp,dp,wp
        use module_types
        use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
        implicit none
        integer, intent(in) :: iproc,nproc,ncong,npsidim
        real(gp), intent(in) :: hx,hy,hz
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        real(dp), intent(out) :: gnrm,gnrm_zero
        real(wp), dimension(npsidim), intent(inout) :: hpsi
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
        type(workarrays_quartic_convolutions),dimension(orbs%norbp),intent(inout),optional :: linear_precond_convol_workarrays !< convolution workarrays for the linear case
        type(workarr_precond),dimension(orbs%norbp),intent(inout),optional :: linear_precond_workarrays !< workarrays for the linear case
      end subroutine preconditionall2

      subroutine partial_density_free(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
            &   hfac,nscatterarr,spinsgn,psir,rho_p,ibyyzz_r) !ex-optional argument
         use module_defs, only: gp,dp,wp
         implicit none
         logical, intent(in) :: rsflag
         integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
         real(gp), intent(in) :: hfac,spinsgn
         integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
         real(wp), dimension(n1i,n2i,n3i,nspinn), intent(in) :: psir
         real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
         integer, dimension(:,:,:), pointer :: ibyyzz_r
      END SUBROUTINE partial_density_free

      subroutine parse_cp2k_files(iproc,basisfile,orbitalfile,nat,ntypes,orbs,iatype,rxyz,&
            &   CP2K,wfn_cp2k)
         use module_base
         use module_types, only: orbitals_data, gaussian_basis
         implicit none
         character(len=*), intent(in) :: basisfile,orbitalfile
         integer, intent(in) :: iproc,nat,ntypes
         type(orbitals_data), intent(in) :: orbs
         integer, dimension(nat), intent(in) :: iatype
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: CP2K
         real(wp), dimension(:,:), pointer :: wfn_cp2k
      END SUBROUTINE parse_cp2k_files

      subroutine read_gaussian_information(orbs,G,coeffs,filename, opt_fillrxyz)
          use module_defs, only: wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         type(orbitals_data), intent(inout) :: orbs
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:), pointer :: coeffs
         logical, optional :: opt_fillrxyz
      END SUBROUTINE read_gaussian_information

      subroutine restart_from_gaussians(iproc,nproc,orbs,Lzd,hx,hy,hz,psi,G,coeffs)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(gaussian_basis), intent(inout) :: G
         real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%norbp), intent(out) :: psi
         real(wp), dimension(:,:), pointer :: coeffs
      END SUBROUTINE restart_from_gaussians

      subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin,&
            orbs,orbse,norbsc_arr,locrad,G,psigau,eks,iversion,mapping,quartic_prefactor)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         integer, intent(inout) :: nvirt
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(gp), intent(out) :: eks
         integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
         real(gp), dimension(at%astruct%nat), intent(out) :: locrad
         type(orbitals_data), intent(out) :: orbse
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:,:), pointer :: psigau
         integer,intent(in) :: iversion !< 1:cubic, 2:linear
         integer,dimension(orbs%norb),intent(in),optional:: mapping
         real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
      END SUBROUTINE inputguess_gaussian_orbitals


     subroutine inputguess_gaussian_orbitals_forLinear(iproc,nproc,norb,at,rxyz,nvirt,nspin,&
          nlr, norbsPerAt, mapping, &
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks,quartic_prefactor)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,nlr,norb
       integer, intent(inout) :: nvirt
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
       integer,dimension(norb),intent(in):: mapping
       integer,dimension(at%astruct%nat),intent(in):: norbsPerAt
       real(gp), intent(out) :: eks
       integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
       real(gp), dimension(at%astruct%nat), intent(out) :: locrad
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
       real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
     END SUBROUTINE inputguess_gaussian_orbitals_forLinear

     subroutine AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,&
          nspin,eks,G,gaucoeff,iorbtolr,mapping,quartic_prefactor)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer, intent(in) :: norbe,iproc
       integer, intent(in) :: norbsc,nspin
       type(atoms_data), intent(in) :: at
       !logical, dimension(4,2,at%natsc), intent(in) :: scorb
       real(gp), dimension(3,at%astruct%nat), intent(in), target :: rxyz
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(gp), intent(out) :: eks
       integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
       real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff !fake interface for passing address
       integer,dimension(orbse%norb), optional, intent(in):: mapping
       real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
      END SUBROUTINE AtomicOrbitals

      subroutine apply_potential(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot,&
            &   ibyyzz_r) !optional
         use module_defs, only: gp,dp,wp
         implicit none
         integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
         real(wp), dimension(-nl1:2*n1+2+nl1,-nl2:2*n2+2+nl2,-nl3:2*n3+2+nl3,nspinor), intent(inout) :: psir
         real(wp), dimension(-nl1:2*n1+2+nl1-4*nbuf,-nl2:2*n2+2+nl2-4*nbuf,-nl3:2*n3+2+nl3-4*nbuf,npot), intent(in) :: pot
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
         real(gp), intent(out) :: epot
      END SUBROUTINE apply_potential

      subroutine read_density(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
            &   nat,rxyz,iatypes, znucl)
        use module_defs, only: gp,dp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(out) :: nspin
         integer, intent(out) ::  n1i,n2i,n3i
         real(gp), intent(out) :: hxh,hyh,hzh
         real(dp), dimension(:,:,:,:), pointer :: rho
         real(gp), dimension(:,:), pointer, optional :: rxyz
         integer, intent(out), optional ::  nat
         integer, dimension(:), pointer, optional :: iatypes, znucl
      END SUBROUTINE read_density

      subroutine read_cube(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
            &   nat,rxyz, iatypes, znucl)
         use module_defs, only: gp,dp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(out) :: nspin
         integer, intent(out) ::  n1i,n2i,n3i
         real(gp), intent(out) :: hxh,hyh,hzh
         real(dp), dimension(:,:,:,:), pointer :: rho
         real(gp), dimension(:,:), pointer   :: rxyz
         integer, intent(out)   ::  nat
         integer, dimension(:), pointer   :: iatypes, znucl
      END SUBROUTINE read_cube

      subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
            &   nat,rxyz, iatypes, znucl)
        use module_defs, only: gp,dp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(out) :: nspin
         integer, intent(out) ::  n1i,n2i,n3i
         real(gp), intent(out) :: hxh,hyh,hzh
         real(dp), dimension(:,:,:,:), pointer :: rho
         real(gp), dimension(:,:), pointer :: rxyz
         integer, intent(out) ::  nat
         integer, dimension(:), pointer :: iatypes, znucl
      END SUBROUTINE read_etsf

      subroutine read_potfile4b2B(filename,n1i,n2i,n3i, rho, alat1, alat2, alat3)
         use module_defs, only: gp,dp,wp
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(out) :: n1i,n2i,n3i
         real(gp) alat1, alat2, alat3, dum, dum1
         ! real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
         real(gp), pointer :: rho(:)
      END SUBROUTINE read_potfile4b2B

      subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc, gaenes, &
            &   iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg )
        use module_defs, only: gp,wp
         use module_types
         implicit none
         logical, intent(in) :: enlargerprb
         integer, intent(in) :: iproc,nspin,ng
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:), pointer :: Gocc
         real(gp), pointer, optional :: gaenes(:)
         integer, pointer, optional :: iorbtolr(:)
         integer, pointer, optional :: iorbto_l(:)
         integer, pointer, optional :: iorbto_m(:)
         integer, pointer, optional :: iorbto_ishell(:)
         integer, pointer, optional :: iorbto_iexpobeg(:)
      END SUBROUTINE gaussian_pswf_basis

      subroutine gaussian_pswf_basis_for_paw(at,rxyz,G,  &
            &   iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels,&
         iorbto_imatrixbeg )
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
         type(gaussian_basis_c), intent(inout) :: G

         integer, pointer :: iorbtolr(:)
         integer, pointer :: iorbto_l(:)
         integer, pointer :: iorbto_paw_nchannels(:)
         integer, pointer :: iorbto_m(:)
         integer, pointer :: iorbto_ishell(:)
         integer, pointer :: iorbto_iexpobeg(:)
         integer, pointer :: iorbto_imatrixbeg(:)

         !local variables
      END SUBROUTINE gaussian_pswf_basis_for_paw


      subroutine local_analysis(iproc,nproc,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
      END SUBROUTINE local_analysis

      subroutine plot_gatom_basis(filename,iat,ngx,G,Gocc,rhocoeff,rhoexpo)
        use module_defs, only: wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(in) :: iat,ngx
         type(gaussian_basis), intent(in) :: G
         real(wp), dimension(:), pointer :: Gocc
         real(wp), dimension((ngx*(ngx+1))/2), intent(out) :: rhoexpo
         real(wp), dimension((ngx*(ngx+1))/2,4), intent(out) :: rhocoeff
      END SUBROUTINE plot_gatom_basis

      subroutine calculate_rhocore(at,d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: i3s,n3d,i3xcsh,n3p
        real(gp), intent(in) :: hxh,hyh,hzh
        type(atoms_data), intent(in) :: at
        type(grid_dimensions), intent(in) :: d
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(wp), dimension(:,:,:,:), pointer :: rhocore
      END SUBROUTINE calculate_rhocore

      subroutine XC_potential(geocode,datacode,iproc,nproc,mpi_comm,n01,n02,n03,xc,hgrids,&
           rho,exc,vxc,nspin,rhocore,potxc,xcstr,dvxcdrho,rhohat)
        use module_defs, only: gp,dp,wp
        use module_xc
        implicit none
        character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
        character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
        integer, intent(in) :: iproc,nproc,n01,n02,n03,nspin,mpi_comm
        type(xc_info), intent(in) :: xc
        real(gp), dimension(3), intent(in) :: hgrids
        real(gp), intent(out) :: exc,vxc
        real(dp), dimension(*), intent(inout) :: rho
        real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
        real(wp), dimension(*), intent(out) :: potxc
        real(dp), dimension(6), intent(out) :: xcstr
        real(dp), dimension(:,:,:,:), target, intent(out), optional :: dvxcdrho
        real(wp), dimension(:,:,:,:), optional :: rhohat
      END SUBROUTINE XC_potential

      subroutine xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
           nxcl,nxcr,xc,hx,hy,hz,rhopot,pot_ion,sumpion,zf,zfionxc,exc,vxc,nproc,nspden)
        use module_defs, only: gp,dp,wp
        use module_xc
        use interfaces_41_xc_lowlevel
        implicit none
        character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
        logical, intent(in) :: sumpion
        integer, intent(in) :: m1,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,nproc,nspden
        integer, intent(in) :: nwbl,nwbr
        real(gp), intent(in) :: hx,hy,hz
        type(xc_info), intent(in) :: xc
        real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rhopot
        real(wp), dimension(*), intent(in) :: pot_ion
        real(dp), dimension(md1,md3,md2/nproc), intent(out) :: zf
        real(wp), dimension(md1,md3,md2/nproc,nspden), intent(out) :: zfionxc
        real(dp), intent(out) :: exc,vxc
      END SUBROUTINE xc_energy

      subroutine direct_minimization(iproc,nproc,in,at,nvirt,rxyz,&
           rhopot,nlpsp,pkernel,dpbox,xc,GPU,KSwfn,VTwfn)
        use module_defs, only: gp,dp,wp
        use module_types
        use module_xc

        implicit none
        integer, intent(in) :: iproc,nproc,nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        type(denspot_distribution), intent(in) :: dpbox
        type(DFT_wavefunction), intent(inout) :: KSwfn,VTwfn
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        type(coulomb_operator), intent(inout) :: pkernel
        real(dp), dimension(*), intent(in), target :: rhopot
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
      end subroutine direct_minimization

      subroutine gaussian_rism_basis(nat,radii,rxyz,G)
        use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: nat
         real(gp), dimension(nat), intent(in) :: radii
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G
      END SUBROUTINE gaussian_rism_basis

      subroutine gaussian_hermite_basis(nhermitemax,nat,radii,rxyz,G)
        use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: nat,nhermitemax
         real(gp), dimension(nat), intent(in) :: radii
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G
      END SUBROUTINE gaussian_hermite_basis

      subroutine write_eigenvalues_data(etol,orbs,mom_vec)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        real(gp), intent(in) :: etol
        type(orbitals_data), intent(in) :: orbs
        real(gp), dimension(:,:,:), intent(in), pointer :: mom_vec
      end subroutine write_eigenvalues_data

      subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         logical, intent(in) :: occorbs
         integer, intent(in) :: iproc,nspin,nvirt,nplot
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
       END SUBROUTINE write_eigen_objects

      subroutine free_full_potential(nproc,flag,xc,pot,subname)
         use module_defs, only: gp,dp,wp
         use module_xc
         implicit none
         character(len=*), intent(in) :: subname
         integer, intent(in) :: nproc,flag
         type(xc_info), intent(in) :: xc
         real(wp), dimension(:), pointer :: pot
      END SUBROUTINE free_full_potential

      subroutine calculate_energy_and_gradient(iter,iproc,nproc,GPU,ncong,iscf,&
           energs,wfn,gnrm,gnrm_zero)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,ncong,iscf,iter
        type(energy_terms), intent(inout) :: energs
        type(GPU_pointers), intent(in) :: GPU
        type(DFT_wavefunction), intent(inout) :: wfn
        real(gp), intent(out) :: gnrm,gnrm_zero
      END SUBROUTINE calculate_energy_and_gradient

      subroutine orthoconstraint(iproc,nproc,orbs,comms,symm,&
            psi,hpsi,scprsum,spsi) !n(c) wfd (arg:5)
        use module_defs, only: gp,dp,wp
        use module_types
        use communications_base, only: comms_cubic
        implicit none
        logical, intent(in) :: symm !< symmetrize the lagrange multiplier after calculation
        integer, intent(in) :: iproc,nproc
        type(orbitals_data), intent(in) :: orbs
        type(comms_cubic), intent(in) :: comms
        !n(c) type(wavefunctions_descriptors), intent(in) :: wfd
        real(wp), dimension(orbs%npsidim_comp), intent(in) :: psi
        real(wp), dimension(orbs%npsidim_comp), intent(inout) :: hpsi
        real(dp), intent(out) :: scprsum
        real(wp), dimension(orbs%npsidim_comp), optional, intent(in) :: spsi
      END SUBROUTINE orthoconstraint


      subroutine constrained_davidson(iproc,nproc,in,at,&
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           hx,hy,hz,rxyz,rhopot,psi,v,dpbox,xc,GPU)
        use module_defs, only: gp,dp,wp
        use module_types
        use communications_base, only: comms_cubic
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc
        integer, intent(in) :: nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        type(comms_cubic), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpbox
        type(xc_info), intent(in) :: xc
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc)
        !v, that is psivirt, is transposed on input and direct on output
      end subroutine constrained_davidson

      subroutine local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,hx,hy,hz,&
           ipotmethod,confdatarr,pot,psi,hpsi,pkernel,xc,alphaSIC,ekin_sum,epot_sum,eSIC_DC,&
           dpbox,potential,comgp)
        use module_defs, only: gp,dp,wp
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,ipotmethod,npsidim_orbs
        real(gp), intent(in) :: hx,hy,hz,alphaSIC
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
        type(xc_info), intent(in) :: xc
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi !this dimension will be modified
        real(wp), dimension(:),pointer :: pot !< the potential, with the dimension compatible with the ipotmethod flag
        !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
        real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        type(coulomb_operator), intent(inout) :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
        type(denspot_distribution),intent(in),optional :: dpbox
        !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
        real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
        type(p2pComms),intent(inout), optional:: comgp
      END SUBROUTINE local_hamiltonian

      subroutine NK_SIC_potential(lr,orbs,xc,fref,hgrids,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
        use module_defs, only: gp,wp,dp
        use module_types
         use module_xc
         implicit none
         real(gp), intent(in) :: fref
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         type(coulomb_operator), intent(inout) :: pkernel
         type(xc_info), intent(in) :: xc
         real(gp), dimension(3), intent(in) :: hgrids
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
         !real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
         real(wp), intent(inout) :: poti
         real(gp), intent(out) :: eSIC_DC
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,2*orbs%nspin), intent(in), optional :: potandrho
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin), intent(out), optional :: wxdsave
      END SUBROUTINE NK_SIC_potential

      subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  &
         wfd,psi,orblist)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,n1,n2,n3, iformat
         real(gp), intent(in) :: hx,hy,hz
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(orbitals_data), intent(inout) :: orbs
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         integer, dimension(orbs%norb), optional :: orblist
         real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
         character(len=*), intent(in) :: filename
      END SUBROUTINE readmywaves

      subroutine writemywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
        use module_types
        use module_defs, only: gp,dp,wp
        use yaml_output
        implicit none
        integer, intent(in) :: iproc,n1,n2,n3,iformat
        real(gp), intent(in) :: hx,hy,hz
        type(atoms_data), intent(in) :: at
        type(orbitals_data), intent(in) :: orbs
        type(wavefunctions_descriptors), intent(in) :: wfd
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
        character(len=*), intent(in) :: filename
      end subroutine writemywaves

      subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out,iiorb)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         logical, intent(in) :: lbin
         integer, intent(in) :: iorb,ispinor,unitfile
         type(orbitals_data), intent(in) :: orbs
         integer, intent(out) :: iorb_out
         integer,intent(in),optional :: iiorb
      END SUBROUTINE open_filename_of_iorb

      subroutine filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out,iiorb)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         logical, intent(in) :: lbin
         integer, intent(in) :: iorb,ispinor
         type(orbitals_data), intent(in) :: orbs
         character(len=*) :: filename_out
         integer, intent(out) :: iorb_out
         integer,intent(in),optional :: iiorb
      END SUBROUTINE filename_of_iorb

      subroutine verify_file_presence(filerad,orbs,iformat,nproc,nforb)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: nproc
        character(len=*), intent(in) :: filerad
        type(orbitals_data), intent(in) :: orbs
        integer, intent(out) :: iformat
        integer, optional, intent(in) :: nforb
      end subroutine verify_file_presence

      subroutine readwavetoisf(lstat, filename, formatted, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        logical, intent(in) :: formatted
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
      END SUBROUTINE readwavetoisf
      subroutine readwavetoisf_etsf(lstat, filename, iorbp, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        integer, intent(in) :: iorbp
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
      END SUBROUTINE readwavetoisf_etsf

      subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: ln
        character(len = ln), intent(in) :: filename
        integer, intent(in) :: iorbp
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
      END SUBROUTINE read_wave_to_isf
      subroutine free_wave_to_isf(psiscf)
        use module_defs, only: gp,dp,wp
        implicit none
        real(wp), dimension(:,:,:,:), pointer :: psiscf
      END SUBROUTINE free_wave_to_isf

      subroutine isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
        use module_defs, only: gp,wp
        use locregs, only: locreg_descriptors
        use locreg_operations, only: workarr_locham
        implicit none
        integer, intent(in) :: nspinor
        real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
        type(locreg_descriptors), intent(in) :: lr
        type(workarr_locham), intent(inout) :: w
        real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psir
        real(gp), intent(out) :: ekin
        real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
        real(wp), dimension(6), optional :: k_strten
      end subroutine isf_to_daub_kinetic

      subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
              fnrm_tmb,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,do_iterative_orthogonalization,sf_per_type,&
          nit_precond,target_function,&
          correction_orthoconstraint,nit_basis,&
          ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
          gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, max_inversion_error, kappa_conv, method_updatekernel,&
          purification_quickreturn, correction_co_contra, &
          precond_convol_workarrays, precond_workarrays, &
          wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi, fnrm, energs_work, frag_calc, &
          cdft, input_frag, ref_frags)
        use module_defs, only: gp,dp,wp
        use module_types
        use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
        use fragment_base, only: fragmentInputParameters
        use module_fragments, only: system_fragment
        use constrained_dft, only: cdft_data
        use communications_base, only: work_transpose
        implicit none

        ! Calling arguments
        integer,intent(in) :: iproc, nproc
        integer,intent(inout) :: order_taylor
        real(kind=8),intent(in) :: max_inversion_error
        integer,intent(out) :: infoBasisFunctions
        type(atoms_data), intent(in) :: at
        type(orbitals_data) :: orbs
        real(kind=8),dimension(3,at%astruct%nat) :: rxyz
        type(DFT_local_fields), intent(inout) :: denspot
        type(GPU_pointers), intent(inout) :: GPU
        real(kind=8),intent(out) :: trH, fnrm_tmb
        real(kind=8),intent(inout) :: trH_old
        type(DFT_PSP_projectors),intent(inout) :: nlpsp
        integer,intent(in) :: scf_mode
        type(localizedDIISParameters),intent(inout) :: ldiis
        type(DFT_wavefunction),target,intent(inout) :: tmb
        type(SIC_data) :: SIC !<parameters for the SIC methods
        type(energy_terms),intent(in) :: energs_base
        logical,intent(in) :: do_iterative_orthogonalization
        integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type
        integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
        real(kind=8),intent(out) :: ratio_deltas
        logical, intent(inout) :: ortho_on
        integer, intent(in) :: extra_states
        integer,intent(in) :: itout
        real(kind=8),intent(in) :: conv_crit, early_stop, gnrm_dynamic, min_gnrm_for_dynamic, kappa_conv
        logical,intent(in) :: experimental_mode, purification_quickreturn
        logical,intent(out) :: can_use_ham
        integer,intent(in) :: method_updatekernel
        logical,intent(in) :: correction_co_contra
        type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
        type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
        type(work_transpose),intent(inout) :: wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi
        type(work_mpiaccumulate),intent(inout) :: fnrm, energs_work
        logical, intent(in) :: frag_calc
        !these must all be present together
        type(cdft_data),intent(inout),optional :: cdft
        type(fragmentInputParameters),optional,intent(in) :: input_frag
        type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
      end subroutine getLocalizedBasis

    subroutine psimix(iproc,nproc,ndim_psi,orbs,comms,diis,hpsit,psit)
      use module_defs, only: gp,dp,wp
      use module_types
      use communications_base, only: comms_cubic
      implicit none
      integer, intent(in) :: iproc,nproc,ndim_psi
      type(orbitals_data), intent(in) :: orbs
      type(comms_cubic), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      real(wp), dimension(ndim_psi), intent(inout) :: psit,hpsit
    end subroutine psimix

    subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
        energs,nlpsp,SIC,tmb,fnrm,calculate_overlap_matrix,invert_overlap_matrix,communicate_phi_for_lsumrho,&
        calculate_ham,extra_states,itout,it_scc,it_cdft,order_taylor,max_inversion_error,purification_quickreturn,&
        calculate_KS_residue,calculate_gap,energs_work,remove_coupling_terms,&
        convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,reorder,cdft, updatekernel)
      use module_defs, only: gp,dp,wp
      use module_types
      use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
      use constrained_dft
      use diis_sd_optimization
      use yaml_output
      use sparsematrix_base, only: sparse_matrix
      implicit none
      integer,intent(in) :: iproc, nproc, scf_mode, itout, it_scc, it_cdft
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(orbitals_data),intent(inout) :: orbs
      type(atoms_data),intent(in) :: at
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers),intent(inout) :: GPU
      integer,intent(out) :: infoCoeff
      type(energy_terms),intent(inout) :: energs
      real(kind=8),intent(inout) :: fnrm
      type(DFT_PSP_projectors),intent(inout) :: nlpsp
      type(SIC_data),intent(in) :: SIC
      type(DFT_wavefunction),intent(inout) :: tmb
      logical,intent(in):: calculate_overlap_matrix, invert_overlap_matrix
      logical,intent(in):: communicate_phi_for_lsumrho, purification_quickreturn
      logical,intent(in) :: calculate_ham, calculate_KS_residue, calculate_gap
      type(work_mpiaccumulate),intent(inout) :: energs_work
      logical,intent(in) :: remove_coupling_terms
      type(DIIS_obj),intent(inout),optional :: ldiis_coeff ! for dmin only
      integer, intent(in), optional :: nitdmin ! for dmin only
      real(kind=gp), intent(in), optional :: convcrit_dmin ! for dmin only
      logical, intent(in), optional :: curvefit_dmin ! for dmin only
      type(cdft_data),intent(inout),optional :: cdft
      integer, intent(in) :: extra_states
      logical, optional, intent(in) :: reorder
      logical, optional, intent(in) :: updatekernel
    end subroutine get_coeff

    subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,rxyz,denspot,rhopotold,nlpsp,GPU,&
           energs,energy,fpulay,infocode,ref_frags,cdft, &
           fdisp, fion)
      use module_defs, only: gp,dp,wp
      use module_types
      use module_fragments
      use constrained_dft
      implicit none
      integer,intent(in):: iproc, nproc
      type(atoms_data),intent(inout):: at
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%astruct%nat),intent(inout):: rxyz
      real(8),dimension(3,at%astruct%nat),intent(out):: fpulay
      type(DFT_local_fields), intent(inout) :: denspot
      real(gp), dimension(*), intent(inout) :: rhopotold
      type(DFT_PSP_projectors),intent(inout):: nlpsp
      type(GPU_pointers),intent(in out):: GPU
      type(energy_terms),intent(inout) :: energs
      real(gp), dimension(:), pointer :: rho,pot
      real(8),intent(out):: energy
      type(DFT_wavefunction),intent(inout),target:: tmb
      type(DFT_wavefunction),intent(inout),target:: KSwfn
      integer,intent(out):: infocode
      type(system_fragment), dimension(:), pointer :: ref_frags
      type(cdft_data), intent(inout) :: cdft
      real(kind=8),dimension(3,at%astruct%nat),intent(in) :: fdisp, fion
    end subroutine linearScaling


   subroutine createDerivativeBasis(n1,n2,n3, &
              nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
              hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
              w_c, w_f, w_f1, w_f2, w_f3, x_c, x_f, y_c, y_f, z_c, z_f)
      use module_defs, only: gp,dp,wp
      implicit none
      integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      real(gp), intent(in) :: hgrid
      integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
      integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
      integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
      real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
      real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
      real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f1
      real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: w_f2
      real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: w_f3
      real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: x_c
      real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f
      real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
      real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
      real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: z_c
      real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: z_f
    end subroutine createDerivativeBasis

    subroutine inputguessConfinement(iproc, nproc, at, input, hx, hy, hz, &
         rxyz, nlpsp, GPU, orbs, kswfn, tmb, denspot, rhopotold, energs,&
         locregcenters)
      ! Input wavefunctions are found by a diagonalization in a minimal basis set
      ! Each processors write its initial wavefunctions into the wavefunction file
      ! The files are then read by readwave
      use module_defs, only: gp,dp,wp
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: iproc,nproc
      real(gp), intent(in) :: hx, hy, hz
      type(atoms_data), intent(inout) :: at
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(GPU_pointers), intent(inout) :: GPU
      type(input_variables),intent(in) :: input
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(orbitals_data),intent(inout) :: orbs
      type(DFT_wavefunction),intent(inout) :: kswfn, tmb
      type(DFT_local_fields), intent(inout) :: denspot
      real(dp), dimension(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin), intent(inout) ::  rhopotold
      type(energy_terms),intent(inout) :: energs
      real(kind=8),dimension(3,at%astruct%nat),intent(in),optional :: locregcenters
    end subroutine inputguessConfinement

    subroutine local_partial_densityLinear(nproc,rsflag,nscatterarr,&
         nrhotot,Lzd,hxh,hyh,hzh,xc,nspin,orbs,mapping,psi,rho)
      use module_defs, only: gp,dp,wp
      use module_types
      use module_xc
      implicit none
      logical, intent(in) :: rsflag
      integer, intent(in) :: nproc
      integer,intent(in):: nrhotot
      integer, intent(in) :: nspin
      real(gp), intent(in) :: hxh,hyh,hzh
      type(xc_info), intent(in) :: xc
      type(local_zone_descriptors), intent(in) :: Lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(orbs%norb),intent(in):: mapping
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
      real(dp),dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1),max(nspin,orbs%nspinor)),intent(out):: rho
    end subroutine local_partial_densityLinear

     subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,Lzde,comms,&
          psi,hpsi,psit,orthpar,passmat,iscf,Tel,occopt,& !mandatory
          orbse,commse,etol,norbsc_arr) !optional
       use module_defs, only: gp,dp,wp
       use module_types
       use communications_base, only: comms_cubic
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
     end subroutine LDiagHam

     subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
       integer,dimension(8),intent(out):: comarr
     end subroutine setCommsParameters

     subroutine optimizeDIIS(iproc, nproc, npsidim, orbs, nspin, lzd, hphi, phi, ldiis, experimental_mode)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nspin
       integer,intent(in):: npsidim
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       real(8),dimension(npsidim),intent(in):: hphi
       real(8),dimension(npsidim),intent(inout):: phi
       type(localizedDIISParameters),intent(inout):: ldiis
       logical,intent(in) :: experimental_mode
     end subroutine optimizeDIIS

     subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, orbs, lzd, comgp, onWhichAtomAll, tag)
       use module_defs, only: gp,dp,wp
       use module_types
       use communications_base, only: p2pComms
       implicit none
       integer,intent(in):: iproc, nproc
       integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       !type(p2pCommsGatherPot),intent(out):: comgp
       type(p2pComms),intent(out):: comgp
       integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
       integer,intent(inout):: tag
     end subroutine initializeCommunicationPotential

    subroutine initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, Glr, locregShape, lborbs)
      use module_defs, only: gp,dp,wp
      use module_atoms, only: atomic_structure
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(local_zone_descriptors),intent(inout):: lzd
      real(8),intent(in):: hx, hy, hz
      type(atomic_structure),intent(in) :: astruct
      type(orbitals_data),intent(in):: orbs
      type(locreg_descriptors),intent(in):: Glr
      character(len=1),intent(in):: locregShape
      type(orbitals_data),optional,intent(in):: lborbs
    end subroutine initLocregs

     subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,&
          Lzd,nlpsp,confdatarr,ngatherarr,Lpot,psi,hpsi,paw,&
          energs,SIC,GPU,xc,pkernel,orbsocc,psirocc)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_xc
       use gaussians, only: gaussian_basis
       implicit none
       integer, intent(in) :: iproc,nproc!,nspin
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors),intent(in) :: Lzd
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       type(SIC_data), intent(in) :: SIC
       type(xc_info), intent(in) :: xc
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
       !real(wp), dimension(lzd%ndimpotisf) :: Lpot
       real(wp), dimension(:),pointer :: Lpot
       type(energy_terms), intent(inout) :: energs
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       type(coulomb_operator), intent(in), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
       !PAW variables:
       type(paw_objects),intent(inout)::paw
     end subroutine FullHamiltonianApplication

       subroutine apply_potential_lr(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,&
            psir,pot,epot,&
            confdata,ibyyzz_r,psir_noconf,econf) !optional
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
         integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
         real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
         real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
         type(confpot_data), intent(in), optional, target :: confdata !< data for the confining potential
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
         real(gp), intent(out) :: epot
         real(wp),dimension(n1i,n2i,n3i,nspinor),intent(inout),optional :: psir_noconf !< real-space wfn in lr where only the potential (without confinement) will be applied
         real(gp), intent(out),optional :: econf
       end subroutine apply_potential_lr

       subroutine psir_to_vpsi(npot,nspinor,lr,pot,vpsir,epot,confdata,vpsir_noconf,econf)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer, intent(in) :: npot,nspinor
         type(locreg_descriptors), intent(in) :: lr !< localization region of the wavefunction
         !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,npot), intent(in) :: pot
         real(wp), intent(in) :: pot
         real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout) :: vpsir
         real(gp), intent(out) :: epot
         type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
         real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout), optional :: vpsir_noconf !< wavefunction with  the potential without confinement applied
         real(gp), intent(out),optional :: econf !< confinement energy
       end subroutine psir_to_vpsi

       subroutine erf_stress(at,rxyz,hxh,hyh,hzh,n1i,n2i,n3i,n3p,iproc,nproc,ngatherarr,rho,tens)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         !passed var
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
         real(gp), intent(in) :: hxh,hyh,hzh
         integer,intent(in) :: n1i,n2i,n3i,n3p,iproc,nproc
         real(kind=8), dimension(n1i*n2i*max(n3p,1)), intent(in), target :: rho
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
         real(dp),dimension(6), intent(out) :: tens
       end subroutine erf_stress

       subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         integer, intent(in) :: linType
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(in) :: atoms
         type(orbitals_data),intent(inout) :: orbs
         real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
       end subroutine check_linear_and_create_Lzd

       subroutine create_LzdLIG(iproc,nproc,nspin,linearmode,hx,hy,hz,Glr,atoms,orbs,rxyz,nl,Lzd)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         real(gp), intent(in):: hx, hy, hz
         type(locreg_descriptors), intent(in) :: Glr
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(in) :: atoms
         type(orbitals_data),intent(inout) :: orbs
         integer, intent(in) :: linearmode
         real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
         type(DFT_PSP_projectors), intent(inout) :: nl
       end subroutine create_LzdLIG

       subroutine export_grids(fname, atoms, rxyz, hx, hy, hz, n1, n2, n3, logrid_c, logrid_f)
         use module_defs, only: gp
         use module_types
         implicit none
         character(len = *), intent(in) :: fname
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
         real(gp), intent(in) :: hx, hy, hz
         integer, intent(in) :: n1, n2, n3
         logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid_c
         logical, dimension(0:n1,0:n2,0:n3), intent(in), optional :: logrid_f
       end subroutine export_grids

       subroutine system_initialization(iproc,nproc,dump,inputpsi,input_wf_format,&
            & dry_run,in,atoms,rxyz,OCLconv,&
            orbs,lnpsidim_orbs,lnpsidim_comp,lorbs,Lzd,Lzd_lin,nlpsp,comms,shift,&
            ref_frags, denspot, locregcenters, inwhichlocreg_old, onwhichatom_old, &
            norb_par_ref, norbu_par_ref, norbd_par_ref,output_grid)
         use module_defs, only: gp,dp,wp
         use f_enums, only: f_enumerator
         use module_types
         use module_fragments
         use communications_base, only: comms_cubic
         implicit none
         integer, intent(in) :: iproc,nproc
         integer, intent(out) :: input_wf_format,lnpsidim_orbs,lnpsidim_comp
         type(f_enumerator), intent(inout) :: inputpsi
         type(input_variables), intent(in) :: in
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
         logical, intent(in) :: OCLconv
         type(orbitals_data), intent(inout) :: orbs,lorbs
         type(local_zone_descriptors), intent(inout) :: Lzd, Lzd_lin
         type(DFT_local_fields), intent(out), optional :: denspot
         type(DFT_PSP_projectors), intent(out) :: nlpsp
         type(comms_cubic), intent(out) :: comms
         real(gp), dimension(3), intent(out) :: shift  !< shift on the initial positions
         !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
         type(system_fragment), dimension(:), pointer :: ref_frags
         real(kind=8),dimension(3,atoms%astruct%nat),intent(inout),optional :: locregcenters
         integer,dimension(:),pointer,optional:: inwhichlocreg_old, onwhichatom_old
         integer,dimension(0:nproc-1),optional:: norb_par_ref, norbu_par_ref, norbd_par_ref !< support function distribution to be used as a reference
         logical, intent(in) :: dry_run, dump
         logical, intent(in), optional :: output_grid
       end subroutine system_initialization

       subroutine input_check_psi_id(inputpsi, input_wf_format, dir_output, orbs, lorbs, iproc, nproc, nfrag, frag_dir, ref_frags)
         use module_types
         use f_enums
         use module_fragments
         implicit none
         integer, intent(out) :: input_wf_format         !< (out) Format of WF
         type(f_enumerator), intent(inout) :: inputpsi   !< (in) indicate how check input psi, (out) give how to build psi
         integer, intent(in) :: iproc                    !< (in)  id proc
         integer, intent(in) :: nproc                    !< (in)  #proc
         integer, intent(in) :: nfrag                    !< number of fragment directories which need checking
         type(system_fragment), dimension(:), pointer :: ref_frags  !< number of orbitals for each fragment
         character(len=100), dimension(nfrag), intent(in) :: frag_dir !< label for fragment subdirectories (blank if not a fragment calculation)
         character(len = *), intent(in) :: dir_output
         type(orbitals_data), intent(in) :: orbs, lorbs
       end subroutine input_check_psi_id

       subroutine extract_potential_for_spectra(iproc,nproc,at,rhod,dpbox,&
            orbs,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
            nlpsp,pkernel,ixc,psi,G,&
            nspin,potshortcut,symObj,GPU,input)
         use module_defs, only: gp,dp,wp
         use module_types
         use communications_base, only: comms_cubic
         implicit none
         !Arguments
         integer, intent(in) :: iproc,nproc,ixc
         integer, intent(inout) :: nspin,nvirt
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(inout) :: at
         type(rho_descriptors),intent(in) :: rhod
         type(denspot_distribution), intent(in) :: dpbox
         type(orbitals_data), intent(inout) :: orbs
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(comms_cubic), intent(in) :: comms
         type(GPU_pointers), intent(inout) :: GPU
         type(input_variables):: input
         type(symmetry_data), intent(in) :: symObj
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
         type(gaussian_basis), intent(out) :: G !basis for davidson IG
         real(wp), dimension(:), pointer :: psi
         real(wp), dimension(:,:,:,:), pointer :: rhocore
         type(coulomb_operator), intent(inout) :: pkernel
         integer, intent(in) ::potshortcut
       end subroutine extract_potential_for_spectra

       subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,iscf,alphamix,&
            nlpsp,linflag,unblock_comms,GPU,wfn,&
            energs,rpnrm,xcstr)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         logical, intent(in) :: scf
         integer, intent(in) :: iproc,nproc,itrp,iscf,linflag,itwfn
         character(len=3), intent(in) :: unblock_comms
         real(gp), intent(in) :: alphamix
         type(atoms_data), intent(in) :: atoms
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         type(DFT_local_fields), intent(inout) :: denspot
         type(energy_terms), intent(inout) :: energs
         type(DFT_wavefunction), intent(inout) :: wfn
         type(GPU_pointers), intent(inout) :: GPU
         real(gp), intent(inout) :: rpnrm
         real(gp), dimension(6), intent(out) :: xcstr
       end subroutine psitohpsi

       subroutine assignToLocreg2(iproc, nproc, norb, norbu, norb_par, natom, nlr, nspin, Localnorb, spinsgn, rxyz, inwhichlocreg)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: nlr,iproc,nproc,nspin,natom,norb,norbu
         integer,dimension(nlr),intent(in):: Localnorb
         real(kind=8),dimension(norb),intent(in):: spinsgn
         integer,dimension(0:nproc-1),intent(in):: norb_par
         real(8),dimension(3,nlr),intent(in):: rxyz
         integer,dimension(:),pointer,intent(out):: inwhichlocreg
       end subroutine assignToLocreg2

       subroutine calc_gradient(geocode,n1,n2,n3,n3grad,deltaleft,deltaright,rhoinp,nspden,hx,hy,hz,&
            gradient,rhocore)
         use module_defs, only: gp,dp,wp
         implicit none
         !Arguments
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(in) :: n1,n2,n3,n3grad,deltaleft,deltaright,nspden
         real(dp), intent(in) :: hx,hy,hz
         real(dp), dimension(n1,n2,n3,nspden), intent(inout) :: rhoinp
         real(dp), dimension(n1,n2,n3grad,2*nspden-1,0:3), intent(out) :: gradient
         real(dp), dimension(:,:,:,:), pointer :: rhocore
       end subroutine calc_gradient


       subroutine destroy_new_locregs(iproc, nproc, tmb)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(DFT_wavefunction),intent(inout):: tmb
       end subroutine destroy_new_locregs

       subroutine define_confinement_data(confdatarr,orbs,rxyz,at,hx,hy,hz,&
                  confpotorder,potentialprefac,Lzd,confinementCenter)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         integer,intent(in):: confpotorder
         real(gp),dimension(at%astruct%ntypes),intent(in):: potentialprefac
         type(local_zone_descriptors), intent(in) :: Lzd
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         integer, dimension(orbs%norb), intent(in) :: confinementCenter
         type(confpot_data), dimension(orbs%norbp), intent(out) :: confdatarr
       end subroutine define_confinement_data

       subroutine update_locreg(iproc, nproc, nlr, locrad, locrad_kernel, locrad_mult, locregCenter, glr_tmp, &
                  useDerivativeBasisFunctions, nscatterarr, hx, hy, hz, astruct, input, &
                  orbs_KS, orbs, lzd, npsidim_orbs, npsidim_comp, lbcomgp, lbcollcom, lfoe, lbcollcom_sr)
         use module_defs, only: gp,dp,wp
         use module_types
         use foe_base, only: foe_data
         use communications_base, only: p2pComms
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         integer,intent(out) :: npsidim_orbs, npsidim_comp
         logical,intent(in):: useDerivativeBasisFunctions
         integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
         real(8),intent(in):: hx, hy, hz
         type(atomic_structure),intent(in) :: astruct
         type(input_variables),intent(in) :: input
         real(8),dimension(nlr),intent(in):: locrad, locrad_kernel, locrad_mult
         type(orbitals_data),intent(in):: orbs_KS, orbs
         real(8),dimension(3,nlr),intent(in):: locregCenter
         type(locreg_descriptors),intent(in):: glr_tmp
         type(local_zone_descriptors),intent(inout):: lzd
         type(p2pComms),intent(inout):: lbcomgp
         type(foe_data),intent(inout),optional :: lfoe
         type(comms_linear),intent(inout):: lbcollcom
         type(comms_linear),intent(inout),optional :: lbcollcom_sr
       end subroutine update_locreg

       subroutine destroy_DFT_wavefunction(wfn)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         type(DFT_wavefunction),intent(inout):: wfn
       end subroutine destroy_DFT_wavefunction

       subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs, &
           norb_par_ref, norbu_par_ref, norbd_par_ref)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspinor
         type(input_variables),intent(in):: input
         type(atomic_structure),intent(in):: astruct
         real(8),dimension(3,astruct%nat),intent(in):: rxyz
         type(orbitals_data),intent(out):: lorbs
         integer,dimension(0:nproc-1),intent(in),optional :: norb_par_ref, norbu_par_ref, norbd_par_ref
       end subroutine init_orbitals_data_for_linear


       subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
                  ldiis, fnrmOldArr, fnrm_old, alpha, trH, trHold, fnrm, alpha_mean, alpha_max, &
                  energy_increased, tmb, lhphiold, overlap_calculated, &
                  energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
                  hpsi_small, experimental_mode, calculate_inverse, correction_co_contra, hpsi_noprecond, &
                  norder_taylor, max_inversion_error, method_updatekernel, precond_convol_workarrays, precond_workarrays,&
                  wt_hphi, wt_philarge, wt_hpsinoprecond, &
                  cdft, input_frag, ref_frags)
         use module_defs, only: gp,dp,wp
         use module_types
         use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
         use communications_base, only: work_transpose
         use sparsematrix_base, only: matrices
         use constrained_dft, only: cdft_data
         use module_fragments, only: system_fragment,fragmentInputParameters
         implicit none
         integer, intent(in) :: iproc, nproc, it, method_updatekernel
         integer,intent(inout) :: norder_taylor
         real(kind=8),intent(in) :: max_inversion_error
         type(DFT_wavefunction), target, intent(inout):: tmb
         type(localizedDIISParameters), intent(inout) :: ldiis
         real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: fnrmOldArr
         real(kind=8),intent(inout) :: fnrm_old
         real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
         real(kind=8), intent(out):: trH, alpha_mean, alpha_max
         type(work_mpiaccumulate), intent(inout):: fnrm
         real(kind=8), intent(in):: trHold
         logical,intent(out) :: energy_increased
         real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
         logical, intent(inout):: overlap_calculated
         type(energy_terms), intent(in) :: energs
         real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c) :: hpsit_c
         real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f) :: hpsit_f
         integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
         logical, intent(in) :: experimental_mode, calculate_inverse, correction_co_contra
         real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
         real(kind=8), dimension(tmb%npsidim_orbs),intent(out) :: hpsi_noprecond
         type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
         type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
         type(work_transpose),intent(inout) :: wt_hphi
         type(work_transpose),intent(inout) :: wt_philarge
         type(work_transpose),intent(out) :: wt_hpsinoprecond
         type(cdft_data),intent(inout),optional :: cdft
         type(fragmentInputParameters), optional, intent(in) :: input_frag
         type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
       end subroutine calculate_energy_and_gradient_linear

       subroutine improveOrbitals(iproc, nproc, tmb, nspin, ldiis, alpha, gradient, experimental_mode)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspin
         type(DFT_wavefunction),intent(inout):: tmb
         type(localizedDIISParameters),intent(inout):: ldiis
         real(8),dimension(tmb%orbs%norbp),intent(in):: alpha
         real(kind=wp),dimension(max(tmb%npsidim_orbs,tmb%npsidim_comp)),intent(inout) :: gradient
         logical,intent(in) :: experimental_mode
       end subroutine improveOrbitals

       subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, at, do_iterative_orthonormalization, sf_per_type, &
                  lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
                  experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in) :: iproc, nproc, it
         integer,intent(inout) :: order_taylor
         real(kind=8),intent(in) :: max_inversion_error
         type(localizedDIISParameters),intent(inout):: ldiis
         type(DFT_wavefunction),target,intent(inout):: tmb
         type(atoms_data),intent(in) :: at
         logical,intent(in) :: do_iterative_orthonormalization
         integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type 
         real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lphiold
         real(8),intent(in):: trH, meanAlpha, alpha_max
         real(8),dimension(tmb%orbs%norbp),intent(inout):: alpha, alphaDIIS
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: hpsi_small
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),optional,intent(out) :: psidiff
         logical, intent(in) :: ortho, experimental_mode
         real(kind=8),intent(out) :: trH_ref
         real(kind=8),dimension(tmb%linmat%l%nvctrp_tg*tmb%linmat%l%nspin),intent(inout) :: kernel_best
         logical,intent(out) :: complete_reset
       end subroutine hpsitopsi_linear

       subroutine DIISorSD(iproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt, trH_ref, kernel_best, complete_reset)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: iproc, it
         real(kind=8),intent(in):: trH
         type(DFT_wavefunction),intent(inout):: tmbopt
         type(localizedDIISParameters),intent(inout):: ldiis
         real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout):: alpha, alphaDIIS
         real(kind=8),dimension(max(tmbopt%npsidim_orbs,tmbopt%npsidim_comp)),intent(out):: lphioldopt
         real(kind=8),intent(out) :: trH_ref
         real(kind=8),dimension(tmbopt%linmat%l%nvctrp_tg*tmbopt%linmat%l%nspin),intent(inout) :: kernel_best
         logical,intent(out) :: complete_reset
       end subroutine DIISorSD

       subroutine psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
            ipotmethod,confdatarr,pot,psi,vpsi,pkernel,xc,alphaSIC,epot_sum,evSIC,vpsi_noconf,econf_sum)
         use module_defs, only: gp,dp,wp
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: iproc,ipotmethod,npsidim_orbs
         real(gp), intent(in) :: alphaSIC
         type(xc_info), intent(in) :: xc
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi !this dimension will be modified
         real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
         real(gp), intent(out) :: epot_sum,evSIC
         real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: vpsi
         type(coulomb_operator), intent(inout) ::  pkernel !< the PSolver kernel which should be associated for the SIC schemes
         real(wp), dimension(orbs%npsidim_orbs), intent(inout),optional :: vpsi_noconf
         real(gp),intent(out),optional :: econf_sum
       end subroutine psi_to_vlocpsi


       subroutine initialize_linear_from_file(iproc,nproc,input_frag,astruct,rxyz,orbs,Lzd,&
              iformat,dir_output,filename,ref_frags,orblist)
         use module_defs, only: gp,dp,wp
         use module_types
         use module_fragments
         implicit none
         integer, intent(in) :: iproc, nproc, iformat
         type(orbitals_data), intent(inout) :: orbs  !< orbs related to the basis functions, inwhichlocreg generated in this routine
         type(atomic_structure), intent(in) :: astruct
         real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
         character(len=*), intent(in) :: filename, dir_output
         type(local_zone_descriptors), intent(inout) :: Lzd !< must already contain Glr and hgrids
         type(fragmentInputParameters), intent(in) :: input_frag
         type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
         integer, dimension(orbs%norb), optional :: orblist
       end subroutine initialize_linear_from_file

        subroutine readmywaves_linear_new(iproc,nproc,dir_output,filename,iformat,at,tmb,rxyz,&
               ref_frags,input_frag,frag_calc,kernel_restart,orblist)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_fragments
          use yaml_output
          implicit none
          integer, intent(in) :: iproc, nproc
          integer, intent(in) :: iformat
          type(atoms_data), intent(in) :: at
          type(DFT_wavefunction), intent(inout) :: tmb
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          !real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
          character(len=*), intent(in) :: dir_output, filename
          type(fragmentInputParameters), intent(in) :: input_frag
          type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
          logical, intent(in) :: frag_calc, kernel_restart
          integer, dimension(tmb%orbs%norb), intent(in), optional :: orblist
        end subroutine readmywaves_linear_new

        subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi)
          use module_defs, only: gp,dp,wp
          implicit none
          integer,intent(in):: npsidim
          real(8),dimension(:),pointer,intent(out):: lphi, lhphi
          character(len=*),intent(in):: subname
        end subroutine allocate_auxiliary_basis_function

        subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi)
          use module_defs, only: gp,dp,wp
          implicit none
          real(8),dimension(:),pointer:: lphi, lhphi
          character(len=*),intent(in):: subname
        end subroutine deallocate_auxiliary_basis_function

        subroutine update_ldiis_arrays(tmb, subname, ldiis)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          type(DFT_wavefunction),intent(in):: tmb
          character(len=*),intent(in):: subname
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine update_ldiis_arrays

        subroutine copy_local_zone_descriptors(lzd_in, lzd_out, subname)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          type(local_zone_descriptors),intent(in):: lzd_in
          type(local_zone_descriptors),intent(inout):: lzd_out
          character(len=*),intent(in):: subname
        end subroutine copy_local_zone_descriptors

        subroutine local_potential_dimensions(iproc,Lzd,orbs,xc,ndimfirstproc)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_xc
          implicit none
          integer, intent(in) :: iproc, ndimfirstproc
          type(local_zone_descriptors), intent(inout) :: Lzd
          type(orbitals_data), intent(inout) :: orbs
          type(xc_info), intent(in) :: xc
        end subroutine local_potential_dimensions

        subroutine initialize_DIIS_coeff(isx, ldiis)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          integer,intent(in):: isx
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine initialize_DIIS_coeff

        subroutine allocate_DIIS_coeff(tmb, ldiis)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          type(DFT_wavefunction),intent(in):: tmb
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine allocate_DIIS_coeff

        subroutine initialize_DFT_local_fields(denspot, ixc, nspden)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
          integer, intent(in) :: ixc, nspden
        end subroutine initialize_DFT_local_fields

        subroutine allocate_diis_objects(idsx,alphadiis,npsidim,nkptsp,nspinor,diis)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          integer, intent(in) :: idsx,npsidim,nkptsp,nspinor !n(m)
          real(gp), intent(in) :: alphadiis
          type(diis_objects), intent(inout) :: diis
        end subroutine allocate_diis_objects

        subroutine check_communications(iproc,nproc,orbs,lzd,comms)
          use module_defs, only: gp,dp,wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: lzd
          type(comms_cubic), intent(in) :: comms
        end subroutine check_communications

        subroutine local_forces(iproc,at,rxyz,hxh,hyh,hzh,&
             n1,n2,n3,n3pi,i3s,n1i,n2i,rho,pot,floc,locstrten,charge)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          !Arguments---------
          type(atoms_data), intent(in) :: at
          integer, intent(in) :: iproc,n1,n2,n3,n3pi,i3s,n1i,n2i
          real(gp), intent(in) :: hxh,hyh,hzh
          real(gp),intent(out) :: charge
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          real(dp), dimension(*), intent(in) :: rho,pot
          real(gp), dimension(3,at%astruct%nat), intent(out) :: floc
          real(gp), dimension(6), intent(out) :: locstrten
        end subroutine local_forces

        subroutine denspot_set_history(denspot, iscf, nspin, &
             & n1i, n2i, & !to be removed arguments when denspot has dimensions
             npulayit)
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
          integer, intent(in) :: iscf, n1i, n2i, nspin
          integer,intent(in),optional :: npulayit
        end subroutine denspot_set_history

        subroutine denspot_free_history(denspot)
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
        end subroutine denspot_free_history

        subroutine kswfn_free_scf_data(KSwfn, freePsit)
          use module_types
          implicit none
          type(DFT_wavefunction), intent(inout) :: KSwfn
          logical, intent(in) :: freePsit
        end subroutine kswfn_free_scf_data

        subroutine evaltoocc(iproc,nproc,filewrite,wf,orbs,occopt)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          logical, intent(in) :: filewrite
          integer, intent(in) :: iproc, nproc
          integer, intent(in) :: occopt
          real(gp), intent(in) :: wf
          type(orbitals_data), intent(inout) :: orbs
        end subroutine evaltoocc

        subroutine cholesky(iproc, nspin,norbIn, psi, &
          orbs, comms, ndim_ovrlp, ovrlp, norbTot, block1, &
          ispinIn, paw)
          !use module_defs, only: gp,dp,wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none

          integer:: iproc,nvctrp,norbIn, nspin, block1, ispinIn
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic):: comms
          real(kind=8),dimension(orbs%npsidim_comp),intent(in out):: psi
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts),1):: ovrlp
          integer,dimension(nspin):: norbTot
          type(paw_objects),optional,intent(inout)::paw
        end subroutine cholesky

        subroutine gsChol(iproc, nproc, psi, orthpar, nspinor,&
          orbs, nspin,ndim_ovrlp,norbArr,comms,paw)
          use module_defs, only: wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc, nproc,nspin
          integer, intent(inout) ::  nspinor
          type(orthon_data), intent(in):: orthpar
          type(orbitals_data):: orbs
          type(comms_cubic), intent(in) :: comms
          integer, dimension(nspin), intent(in) :: norbArr
          integer, dimension(nspin,0:orbs%nkpts), intent(inout) :: ndim_ovrlp
          real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psi
          type(paw_objects),optional,intent(inout)::paw
        end subroutine gsCHol

        subroutine loewdin(iproc, norbIn, block1, ispinIn,&
          orbs, comms, nspin, psit, ovrlp, ndim_ovrlp, norbTot, paw)
          !use module_defs, only: gp,dp,wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer,intent(in):: iproc,norbIn, nspin, block1, ispinIn
          type(orbitals_data),intent(in):: orbs
          type(comms_cubic),intent(in):: comms
          real(kind=8),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in out):: psit
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
          integer,dimension(nspin):: norbTot
          type(paw_objects),optional,intent(inout)::paw
        end subroutine loewdin

        subroutine gramschmidt(iproc, norbIn, psit, ndim_ovrlp, ovrlp, orbs, nspin,&
          nspinor, comms, norbTot, block1, block2, ispinIn,paw)
          use module_defs, only: wp !module_base
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer,intent(in):: iproc, norbIn, nspin, block1, block2, ispinIn
          integer, intent(out) :: nspinor
          type(orbitals_data):: orbs
          type(comms_cubic), intent(in) :: comms
          type(paw_objects),optional,intent(inout)::paw
          real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psit
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
          integer,dimension(nspin):: norbTot
        end subroutine gramschmidt

        subroutine orthogonalize(iproc,nproc,orbs,comms,psi,orthpar,paw)
          use module_defs, only: wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic), intent(in) :: comms
          type(orthon_data), intent(in) :: orthpar
          real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
          type(paw_objects),optional,intent(inout) :: paw
        end subroutine orthogonalize

        subroutine calculate_density_kernel(iproc, nproc, isKernel, orbs, orbs_tmb, &
                   coeff, denskern, denskern_, keep_uncompressed_)
          !use module_defs, only: gp,dp,wp
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in):: iproc, nproc
          logical, intent(in) :: isKernel
          type(orbitals_data),intent(in):: orbs, orbs_tmb
          type(sparse_matrix), intent(in) :: denskern
          real(kind=8),dimension(denskern%nfvctr,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
          type(matrices), intent(out) :: denskern_
          logical,intent(in),optional :: keep_uncompressed_ !< keep the uncompressed kernel in denskern_%matrix (requires that this array is already allocated outside of the routine)
        end subroutine calculate_density_kernel
!!$
!!$        subroutine reconstruct_kernel(iproc, nproc, inversion_method, &
!!$                   blocksize_dsyev, blocksize_pdgemm, orbs, tmb, overlap_calculated)
!!$          !use module_defs, only: gp,dp,wp
!!$          use module_types
!!$          implicit none
!!$          integer,intent(in):: iproc, nproc, blocksize_dsyev, blocksize_pdgemm, inversion_method
!!$          type(orbitals_data),intent(in):: orbs
!!$          type(DFT_wavefunction),intent(inout):: tmb
!!$          logical,intent(inout):: overlap_calculated
!!$        end subroutine reconstruct_kernel


!!$        subroutine reorthonormalize_coeff(iproc, nproc, norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, basis_orbs, &
!!$                   basis_overlap, KS_overlap, basis_overlap_mat, coeff, orbs)
!!$          !use module_defs, only: gp,dp,wp
!!$          use module_types
!!$          use sparsematrix_base, only: sparse_matrix, matrices
!!$          implicit none
!!$          integer, intent(in) :: iproc, nproc, norb
!!$          integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
!!$          type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
!!$          type(sparse_matrix),intent(inout) :: basis_overlap
!!$          type(sparse_matrix),dimension(basis_overlap%nspin),intent(inout) :: KS_overlap
!!$          type(matrices),intent(inout) :: basis_overlap_mat
!!$          real(kind=8),dimension(basis_overlap%nfvctr,norb),intent(inout) :: coeff
!!$          type(orbitals_data), intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
!!$        end subroutine reorthonormalize_coeff

        subroutine copy_old_supportfunctions(iproc,orbs,lzd,phi,lzd_old,phi_old)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          integer,intent(in) :: iproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: lzd
          type(local_zone_descriptors), intent(inout) :: lzd_old
          real(wp), dimension(:), pointer :: phi,phi_old
        end subroutine copy_old_supportfunctions

        subroutine input_memory_linear(iproc, nproc, at, KSwfn, tmb, tmb_old, denspot, input, &
                   rxyz_old, rxyz, denspot0, energs, nlpsp, GPU, ref_frags, cdft)
          use module_defs, only: gp,dp,wp
          use module_types, only: DFT_wavefunction,DFT_local_fields,input_variables,&
               energy_terms,Dft_psp_projectors,GPU_pointers,atoms_data
          use module_fragments, only: system_fragment
          use constrained_dft, only: cdft_data
          implicit none
          integer,intent(in) :: iproc, nproc
          type(atoms_data), intent(inout) :: at
          type(DFT_wavefunction),intent(inout):: KSwfn
          type(DFT_wavefunction),intent(inout):: tmb, tmb_old
          type(DFT_local_fields), intent(inout) :: denspot
          type(input_variables),intent(in):: input
          real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz_old, rxyz
          real(8),dimension(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)),intent(out):: denspot0
          type(energy_terms),intent(inout):: energs
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(GPU_pointers), intent(inout) :: GPU
          type(system_fragment), dimension(:), intent(in) :: ref_frags
          type(cdft_data), intent(inout) :: cdft
        end subroutine input_memory_linear

        subroutine copy_old_coefficients(norb_tmb, nfvctr, coeff, coeff_old)
          use module_defs, only: gp,dp,wp
          implicit none
          integer,intent(in):: norb_tmb, nfvctr
          real(8),dimension(:,:),pointer:: coeff, coeff_old
        end subroutine copy_old_coefficients

        subroutine copy_old_inwhichlocreg(norb_tmb, inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old)
          use module_defs, only: gp,dp,wp
          implicit none
          integer,intent(in):: norb_tmb
          integer,dimension(:),pointer:: inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old
        end subroutine copy_old_inwhichlocreg

        subroutine reformat_supportfunctions(iproc,nproc,at,rxyz_old,rxyz,add_derivatives,tmb,ndim_old,lzd_old,&
               frag_trans,psi_old,input_dir,input_frag,ref_frags,max_shift,phi_array_old)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_fragments
          implicit none
          integer, intent(in) :: iproc,nproc
          integer, intent(in) :: ndim_old
          type(atoms_data), intent(in) :: at
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
          type(DFT_wavefunction), intent(inout) :: tmb
          type(local_zone_descriptors), intent(inout) :: lzd_old
          type(fragment_transformation), dimension(tmb%orbs%norbp), intent(in) :: frag_trans
          real(wp), dimension(:), pointer :: psi_old
          type(phi_array), dimension(tmb%orbs%norbp), optional, intent(in) :: phi_array_old
          logical, intent(in) :: add_derivatives
          character(len=*), intent(in) :: input_dir
          type(fragmentInputParameters), intent(in) :: input_frag
          type(system_fragment), dimension(:), intent(in) :: ref_frags
          real(gp),intent(out) :: max_shift
        end subroutine reformat_supportfunctions

        subroutine reformat_one_supportfunction(llr,llr_old,geocode,hgrids_old,n_old,psigold,&
             hgrids,n,centre_old,centre_new,da,frag_trans,psi,psirold)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_fragments
          implicit none
          integer, dimension(3), intent(in) :: n,n_old
          real(gp), dimension(3), intent(in) :: hgrids,hgrids_old
          !type(wavefunctions_descriptors), intent(in) :: wfd
          type(locreg_descriptors), intent(in) :: llr, llr_old
          character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
          real(gp), dimension(3), intent(inout) :: centre_old,centre_new,da
          type(fragment_transformation), intent(in) :: frag_trans
          real(wp), dimension(0:n_old(1),2,0:n_old(2),2,0:n_old(3),2), intent(in) :: psigold
          real(wp), dimension(llr%wfd%nvctr_c+7*llr%wfd%nvctr_f), intent(out) :: psi
          real(wp), dimension(llr_old%d%n1i,llr_old%d%n2i,llr_old%d%n3i), optional, intent(in) :: psirold
        end subroutine reformat_one_supportfunction

        subroutine input_wf_memory_new(nproc,iproc, atoms, &
                 rxyz_old, hx_old, hy_old, hz_old, psi_old,lzd_old, &
                 rxyz,psi,orbs,lzd)
          use module_defs
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc
          type(atoms_data), intent(in) :: atoms
          real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz, rxyz_old
          real(gp), intent(in) :: hx_old, hy_old, hz_old
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: lzd_old
          type(local_zone_descriptors), intent(in) :: lzd
          real(wp), dimension(:), pointer :: psi, psi_old
        end subroutine input_wf_memory_new

        subroutine integral_equation(iproc,nproc,atoms,wfn,ngatherarr,local_potential,GPU,xc,nlpsp,rxyz,paw)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_xc
          implicit none
          integer, intent(in) :: iproc,nproc
          type(atoms_data), intent(in) :: atoms
          type(DFT_wavefunction), intent(in) :: wfn
          type(GPU_pointers), intent(inout) :: GPU
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(xc_info), intent(in) :: xc
          type(paw_objects), intent(inout) :: paw
          integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
          real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
          real(dp), dimension(:), pointer :: local_potential
        end subroutine integral_equation

        subroutine atoms_new(atoms)
          use module_types
          implicit none
          type(atoms_data), pointer :: atoms
        end subroutine atoms_new

        subroutine inputs_new(in)
          use module_input_keys, only: input_variables
          implicit none
          type(input_variables), pointer :: in
        end subroutine inputs_new

        subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, &
               sd_fit_curve, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, num_extra)
          use module_defs, only: gp,dp,wp
          use module_types
          use diis_sd_optimization
          implicit none
          integer,intent(in):: iproc, nproc, itmax, itout, it_scc, it_cdft
          integer,intent(inout) :: order_taylor
          real(kind=8),intent(in) :: max_inversion_error
          type(orbitals_data),intent(in):: orbs
          type(DFT_wavefunction),intent(inout):: tmb
          type(DIIS_obj), intent(inout) :: ldiis_coeff
          real(kind=gp),intent(in):: fnrm_crit
          real(kind=gp),intent(out):: fnrm
          real(kind=gp), intent(inout) :: energy
          logical, intent(in) :: sd_fit_curve
          real(kind=gp), intent(in) :: factor
          integer, optional, intent(in) :: num_extra
          logical, optional, intent(in) :: reorder
        end subroutine optimize_coeffs

        subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
          use module_types
          implicit none

          ! Calling arguments
          integer, intent(in) :: iproc, nproc, num_extra
          type(dft_wavefunction), intent(inout) :: tmb
          type(orbitals_data), intent(in) :: ksorbs
          real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f
        end subroutine calculate_residue_ks

        subroutine write_energies(iter,iscf,energs,gnrm,gnrm_zero,comment,only_energies)
          use module_defs, only: gp,dp,wp
          use module_types
          use yaml_output
          implicit none
          integer, intent(in) :: iter,iscf
          type(energy_terms), intent(in) :: energs
          real(gp), intent(in) :: gnrm,gnrm_zero
          character(len=*), intent(in) :: comment
          logical,intent(in),optional :: only_energies
        end subroutine write_energies


        subroutine applyprojectorsonthefly(iproc,orbs,at,lr,&
             rxyz,hx,hy,hz,wfd,nlpsp,psi,hpsi,eproj_sum,&
             paw)
          use module_defs, only: gp,dp,wp
          use module_types
          use gaussians, only:gaussian_basis
          implicit none
          integer, intent(in) :: iproc
          real(gp), intent(in) :: hx,hy,hz
          type(atoms_data), intent(in) :: at
          type(orbitals_data), intent(in) :: orbs
          type(wavefunctions_descriptors), intent(in) :: wfd
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(locreg_descriptors),intent(in) :: lr
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
          real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: hpsi
          real(gp), intent(out) :: eproj_sum
          type(paw_objects),optional,intent(inout)::paw
        end subroutine applyprojectorsonthefly


        subroutine allocate_precond_arrays(orbs, lzd, confdatarr, precond_convol_workarrays, precond_workarrays)
          use module_types
          use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
          implicit none
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          type(confpot_data),dimension(orbs%norbp),intent(in) ::  confdatarr
          type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
          type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays
        end subroutine allocate_precond_arrays

        subroutine deallocate_precond_arrays(orbs, lzd, precond_convol_workarrays, precond_workarrays)
          use module_defs, only: gp,dp,wp
          use module_types
          use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
          implicit none
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
          type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays
        end subroutine deallocate_precond_arrays

        subroutine plot_wf(units_provided,orbname,nexpo,at,factor,lr,hx,hy,hz,rxyz,psi, &
                   unit0_, unitx_, unity_, unitz_)
          use module_defs, only: gp,dp,wp
          use locregs, only: locreg_descriptors
          use module_types, only: atoms_data
          implicit none
          logical,intent(in) :: units_provided
          character(len=*) :: orbname
          integer, intent(in) :: nexpo
          real(dp), intent(in) :: factor
          real(gp), intent(in) :: hx,hy,hz
          type(atoms_data), intent(in) :: at
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          type(locreg_descriptors), intent(in) :: lr
          real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
          integer,intent(in),optional :: unit0_, unitx_, unity_, unitz_
        end subroutine plot_wf

        subroutine write_orbital_density(iproc, transform_to_global, iformat, &
                   filename, npsidim, psi, input, orbs, lzd_g, at, rxyz, lzd_l)
          use module_defs, only: gp,dp,wp
          use module_types
          implicit none
          logical,intent(in) :: transform_to_global
          character(len=*),intent(in) :: filename
          integer,intent(in) :: iproc, npsidim, iformat
          real(kind=8),dimension(npsidim),intent(in),target :: psi
          type(input_variables),intent(in) :: input
          type(orbitals_data),intent(in) :: orbs !< orbitals descriptors
          type(local_zone_descriptors),intent(inout) :: lzd_g !< global descriptors
          type(atoms_data),intent(in) :: at
          real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
          type(local_zone_descriptors),intent(in),optional :: lzd_l !< local descriptors
        end subroutine write_orbital_density

  end interface
END MODULE module_interfaces

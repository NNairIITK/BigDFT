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
      subroutine kswfn_optimization_loop(iproc, nproc, o, &
           & alphamix, idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, GPU, xcstr, &
           & in)
        use module_base
        use module_types
        implicit none
        real(dp), dimension(6), intent(out) :: xcstr
        integer, intent(in) :: iproc, nproc, idsx, inputpsi
        real(gp), intent(in) :: alphamix
        type(DFT_optimization_loop), intent(inout) :: o
        type(DFT_wavefunction), intent(inout) :: KSwfn
        type(DFT_local_fields), intent(inout) :: denspot
        type(energy_terms), intent(inout) :: energs
        type(atoms_data), intent(in) :: atoms
        type(GPU_pointers), intent(inout) :: GPU
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        type(input_variables), intent(in) :: in !<todo: Remove me
      END SUBROUTINE kswfn_optimization_loop

     subroutine timing(iproc,category,action)
       implicit none
       integer, intent(in) :: iproc
       character(len=*), intent(in) :: category
       character(len=2), intent(in) :: action
     end subroutine timing

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

      subroutine system_properties(iproc,nproc,in,at,orbs)
        use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(input_variables), intent(in) :: in
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(inout) :: orbs
         !real(gp), dimension(at%astruct%ntypes,3), intent(out) :: radii_cf
      END SUBROUTINE system_properties

      subroutine system_size(atoms,rxyz,crmult,frmult,hx,hy,hz,OCLconv,Glr,shift)
        use module_defs, only: gp
         use module_types
         implicit none
         type(atoms_data), intent(inout) :: atoms
         real(gp), intent(in) :: crmult,frmult
         real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
         !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
         real(gp), intent(inout) :: hx,hy,hz
         logical, intent(in) :: OCLconv
         type(locreg_descriptors), intent(out) :: Glr
         real(gp), dimension(3), intent(out) :: shift
      END SUBROUTINE system_size

      subroutine standard_inputfile_names(inputs, radical)
         use module_types
         implicit none
         type(input_variables), intent(out) :: inputs
         character(len = *), intent(in) :: radical
      END SUBROUTINE standard_inputfile_names

      subroutine read_input_dict_from_files(radical, mpi_env,dict)
        use dictionaries, only: dictionary
        use wrapper_MPI, only: mpi_environment
        implicit none
        character(len = *), intent(in) :: radical
        type(mpi_environment), intent(in) :: mpi_env
        type(dictionary), pointer :: dict
      end subroutine read_input_dict_from_files

      subroutine create_log_file(dict, writing_directory, dir_output, run_name)
        use dictionaries
        implicit none
        type(dictionary), pointer :: dict
        character(len = max_field_length), intent(out) :: writing_directory, dir_output, run_name
      end subroutine create_log_file

      subroutine inputs_from_dict(in, atoms, dict)
        use module_types
        use module_defs
        use dictionaries
        implicit none
        type(input_variables), intent(out) :: in
        type(atoms_data), intent(out) :: atoms
        type(dictionary), pointer :: dict
      end subroutine inputs_from_dict

      subroutine kpt_input_analyse(iproc, in, dict, sym, geocode, alat)
        use module_base, only: gp
        use module_atoms, only: symmetry_data
        use module_types
        use dictionaries
        implicit none
        integer, intent(in) :: iproc
        type(input_variables), intent(inout) :: in
        type(dictionary), pointer, intent(in) :: dict
        type(symmetry_data), intent(in) :: sym
        character(len = 1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
        real(gp), intent(in) :: alat(3)
      end subroutine kpt_input_analyse

      subroutine MemoryEstimator(nproc,idsx,lr,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,mem)
         !n(c) use module_base
         use module_types
         implicit none
         !Arguments
         integer, intent(in) :: nproc,idsx,norb,nspin,nprojel
         integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
         type(locreg_descriptors), intent(in) :: lr
         type(memory_estimation), intent(out) :: mem
      END SUBROUTINE MemoryEstimator

      subroutine check_closed_shell(orbs,lcs)
         !n(c) use module_base
         use module_types
         implicit none
         type(orbitals_data), intent(in) :: orbs
         logical, intent(out) :: lcs
      END SUBROUTINE check_closed_shell

      subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor, &
                 nkpt,kpt,wkpt,orbs,simple,basedist,basedistu,basedistd)
         use module_defs, only: gp
         use module_types
         implicit none
         logical, intent(in) :: simple !< simple calculation of the repartition
         integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
         integer, intent(in) :: nspinor
         type(orbitals_data), intent(inout) :: orbs
         real(gp), dimension(nkpt), intent(in) :: wkpt
         real(gp), dimension(3,nkpt), intent(in) :: kpt
         integer, dimension(0:nproc-1), intent(in), optional :: basedist
         integer, dimension(0:nproc-1), intent(in), optional :: basedistu
         integer, dimension(0:nproc-1), intent(in), optional :: basedistd
      END SUBROUTINE orbitals_descriptors

     subroutine orbitals_descriptors_forLinear(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
       use module_defs, only: gp
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
       integer, intent(in) :: nspinor
       type(orbitals_data), intent(out) :: orbs
       real(gp), dimension(nkpt), intent(in) :: wkpt
       real(gp), dimension(3,nkpt), intent(in) :: kpt
     END SUBROUTINE orbitals_descriptors_forLinear

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

      subroutine dpbox_set(dpbox,Lzd,xc,iproc,nproc,mpi_comm,PS_groupsize,SICapproach,geocode,nspin)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,mpi_comm,PS_groupsize,nspin
        character(len=1), intent(in) :: geocode
        character(len=4), intent(in) :: SICapproach
        type(local_zone_descriptors), intent(in) :: Lzd
        type(xc_info), intent(in) :: xc
        type(denspot_distribution), intent(out) :: dpbox
      end subroutine dpbox_set

      subroutine density_descriptors(iproc,nproc,xc,nspin,crmult,frmult,atoms,dpbox,&
           rho_commun,rxyz,rhodsc)
        use module_defs, only: gp
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        type(xc_info), intent(in) :: xc
        real(gp), intent(in) :: crmult,frmult
        type(atoms_data), intent(in) :: atoms
        type(denspot_distribution), intent(in) :: dpbox
        character(len=3), intent(in) :: rho_commun
        real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
        !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
        type(rho_descriptors), intent(out) :: rhodsc
      end subroutine density_descriptors

      subroutine default_confinement_data(confdatarr,norbp)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: norbp
        type(confpot_data), dimension(norbp), intent(out) :: confdatarr
      end subroutine default_confinement_data

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
         type(coulomb_operator), intent(in) :: pkernel
         real(gp), intent(out) :: eion,edisp,psoffset
         real(dp), dimension(6),intent(out) :: ewaldstr
         real(gp), dimension(:,:), pointer :: fion,fdisp
         real(dp), dimension(*), intent(out) :: pot_ion
       END SUBROUTINE IonicEnergyandForces

       subroutine createIonicPotential(geocode,iproc,nproc,verb,at,rxyz,&
            hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,&
            pot_ion,psoffset)
         use module_defs, only: gp,wp
         use module_types
         implicit none
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
         logical, intent(in) :: verb
         real(gp), intent(in) :: hxh,hyh,hzh,psoffset
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3), intent(in) :: elecfield
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         type(coulomb_operator), intent(in) :: pkernel
         real(wp), dimension(*), intent(inout) :: pot_ion
       END SUBROUTINE createIonicPotential

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
         use module_defs
         use module_types
         implicit none
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_random

       subroutine input_wf_cp2k(iproc, nproc, nspin, atoms, rxyz, Lzd, &
            & psi, orbs)
         use module_defs
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
         use module_defs
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
         use module_base
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
         use module_defs
         use module_types
         use module_fragments
         use constrained_dft
         implicit none
         integer, intent(in) :: iproc, nproc, inputpsi,  input_wf_format
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

       subroutine reformatmywaves(iproc,orbs,at,&
            &   hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
         hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
         use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,n1_old,n2_old,n3_old,n1,n2,n3
         real(gp), intent(in) :: hx_old,hy_old,hz_old,hx,hy,hz
         type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
         real(wp), dimension(wfd_old%nvctr_c+7*wfd_old%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: psi
      END SUBROUTINE reformatmywaves


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
        type(coulomb_operator), intent(in) :: pkernel
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
        use module_base
        use module_types
        implicit none
        logical, intent(in) :: keep_rhop !< preserves the total density in the rho_p array
        integer, intent(in) :: nspin
        type(rho_descriptors),intent(in) :: rhodsc
        type(denspot_distribution), intent(in) :: dpbox
        real(dp), dimension(:,:), pointer :: rho_p !< partial density in orbital distribution scheme
        real(dp), dimension(max(dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d,1),nspin), intent(out) :: rho
      END SUBROUTINE communicate_density

      subroutine rho_segkey(iproc,at,rxyz,crmult,frmult,&
            &   n1i,n2i,n3i,hxh,hyh,hzh,nspin,rho_d,iprint)
        use module_defs, only: gp
        use module_types
         implicit none
         integer,intent(in) :: n1i,n2i,n3i,iproc,nspin
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
         !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
         logical,intent(in) :: iprint
         type(rho_descriptors),intent(inout) :: rho_d
       END SUBROUTINE rho_segkey

       subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
            Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
            energs,SIC,GPU,PotOrKin,xc,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
         use module_base
         use module_types
         use module_xc
         use communications_base, only: p2pComms
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

       subroutine NonLocalHamiltonianApplication(iproc,at,npsidim_orbs,orbs,&
           Lzd,nlpsp,psi,hpsi,eproj_sum,paw)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc, npsidim_orbs
        type(atoms_data), intent(in) :: at
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        real(gp), intent(out) :: eproj_sum
        type(paw_objects),intent(inout)::paw
      END SUBROUTINE NonLocalHamiltonianApplication

      subroutine SynchronizeHamiltonianApplication(nproc,npsidim_orbs,orbs,Lzd,GPU,xc,hpsi,&
           ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: nproc,npsidim_orbs
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
        real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
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
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,iter
        real(wp), intent(out) :: evsum
        type(DFT_wavefunction), intent(inout) :: wfn
        logical, optional :: opt_keeppsit
      END SUBROUTINE last_orthon

      subroutine kswfn_post_treatments(iproc, nproc, KSwfn, tmb, linear, &
           & fxyz, fnoise, fion, fdisp, fpulay, &
           & strten, pressure, ewaldstr, xcstr, &
           & GPU, denspot, atoms, rxyz, nlpsp, &
           & output_denspot, dir_output, gridformat, refill_proj, calculate_dipole)
        use module_base
        use module_types

        implicit none

        type(DFT_wavefunction), intent(in) :: KSwfn
        type(DFT_wavefunction), intent(inout) :: tmb
        type(GPU_pointers), intent(inout) :: GPU
        type(DFT_local_fields), intent(inout) :: denspot
        type(atoms_data), intent(in) :: atoms
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        logical, intent(in) :: refill_proj, calculate_dipole, linear
        integer, intent(in) :: output_denspot, iproc, nproc
        character(len = *), intent(in) :: dir_output
        character(len = *), intent(in) :: gridformat
        real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
        real(gp), dimension(3, atoms%astruct%nat), intent(in) :: fdisp, fion, fpulay
        real(dp), dimension(6), intent(in) :: ewaldstr, xcstr
        real(gp), intent(out) :: fnoise, pressure
        real(gp), dimension(6), intent(out) :: strten
        real(gp), dimension(3, atoms%astruct%nat), intent(out) :: fxyz
      END SUBROUTINE kswfn_post_treatments

      subroutine calculate_forces(iproc,nproc,psolver_groupsize,Glr,atoms,orbs,nlpsp,rxyz,hx,hy,hz,i3s,n3p,nspin,&
           refill_proj,ngatherarr,rho,pot,potxc,nsize_psi,psi,fion,fdisp,fxyz,&
           ewaldstr,hstrten,xcstr,strten,fnoise,pressure,psoffset,imode,tmb,fpulay)
        use module_base
        use module_types
        implicit none
        logical, intent(in) :: refill_proj
        integer, intent(in) :: iproc,nproc,i3s,n3p,nspin,psolver_groupsize,imode,nsize_psi
        real(gp), intent(in) :: hx,hy,hz,psoffset
        type(locreg_descriptors), intent(in) :: Glr
        type(atoms_data), intent(in) :: atoms
        type(orbitals_data), intent(in) :: orbs
        type(DFT_PSP_projectors), intent(inout) :: nlpsp
        integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
        real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: rho,pot,potxc
        real(wp), dimension(nsize_psi), intent(in) :: psi
        real(gp), dimension(6), intent(in) :: ewaldstr,hstrten,xcstr
        real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz,fion,fdisp,fpulay
        real(gp), intent(out) :: fnoise,pressure
        real(gp), dimension(6), intent(out) :: strten
        real(gp), dimension(3,atoms%astruct%nat), intent(out) :: fxyz
        type(DFT_wavefunction),intent(inout) :: tmb
      END SUBROUTINE calculate_forces

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
         real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
         !real(kind=8), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
         real(kind=8), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
         real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
         real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
         type(paw_objects),optional,intent(inout)::paw
      END SUBROUTINE CalculateTailCorrection

      !added for abinit compatilbility
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
         use module_base
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
        use module_base
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
        type(coulomb_operator), intent(in) :: pkernel
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc)
      end subroutine davidson

      subroutine preconditionall(orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
        use module_defs, only: gp,wp,dp
         use module_types
         implicit none
         integer, intent(in) :: ncong
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         real(dp), intent(out) :: gnrm,gnrm_zero
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp,orbs%nspinor), intent(inout) :: hpsi
      END SUBROUTINE preconditionall

      subroutine preconditionall2(iproc,nproc,orbs,Lzd,hx,hy,hz,ncong,npsidim,hpsi,confdatarr,gnrm,gnrm_zero, &
                 linear_precond_convol_workarrays, linear_precond_workarrays)
        use module_base
        use module_types
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
         use module_base
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
        use module_defs, only: gp,wp
         use module_types
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
       use module_base
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
       use module_base
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
         use module_base
         implicit none
         integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
         real(wp), dimension(-nl1:2*n1+2+nl1,-nl2:2*n2+2+nl2,-nl3:2*n3+2+nl3,nspinor), intent(inout) :: psir
         real(wp), dimension(-nl1:2*n1+2+nl1-4*nbuf,-nl2:2*n2+2+nl2-4*nbuf,-nl3:2*n3+2+nl3-4*nbuf,npot), intent(in) :: pot
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
         real(gp), intent(out) :: epot
      END SUBROUTINE apply_potential

      subroutine plot_density(iproc,nproc,filename,at,rxyz,box,nspin,rho)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        type(atoms_data), intent(in) :: at
        type(denspot_distribution), intent(in) :: box
        character(len=*), intent(in) :: filename
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(dp), dimension(max(box%ndimpot,1),nspin), target, intent(in) :: rho
      END SUBROUTINE plot_density

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
         use module_base
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
         use module_base
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
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: i3s,n3d,i3xcsh,n3p
        real(gp), intent(in) :: hxh,hyh,hzh
        type(atoms_data), intent(in) :: at
        type(grid_dimensions), intent(in) :: d
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(wp), dimension(:,:,:,:), pointer :: rhocore
      END SUBROUTINE calculate_rhocore

      subroutine XC_potential(geocode,datacode,iproc,nproc,mpi_comm,n01,n02,n03,xc,hx,hy,hz,&
           rho,exc,vxc,nspin,rhocore,potxc,xcstr,dvxcdrho,rhohat)
        use module_base
        use module_xc
        implicit none
        character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
        character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
        integer, intent(in) :: iproc,nproc,n01,n02,n03,nspin,mpi_comm
        type(xc_info), intent(in) :: xc
        real(gp), intent(in) :: hx,hy,hz
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
        use module_base
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
        use module_base
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
        type(coulomb_operator), intent(in) :: pkernel
        real(dp), dimension(*), intent(in), target :: rhopot
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
      end subroutine direct_minimization

      subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
            &   hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(in) :: iproc,nproc,n3pi,i3s
         real(gp), intent(in) :: hxh,hyh,hzh
         real(gp), dimension(3), intent(in) :: shift
         type(input_variables), intent(in) :: in
         type(grid_dimensions), intent(in) :: grid
         type(coulomb_operator), intent(in) :: pkernel
         real(wp), dimension(*), intent(inout) :: pot_ion
      END SUBROUTINE CounterIonPotential

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
        use module_base
        use module_types
        implicit none
        real(gp), intent(in) :: etol
        type(orbitals_data), intent(in) :: orbs
        real(gp), dimension(:,:,:), intent(in), pointer :: mom_vec
      end subroutine write_eigenvalues_data

      subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt,output_wf_format)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         logical, intent(in) :: occorbs
         integer, intent(in) :: iproc,nspin,nvirt,nplot,output_wf_format
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
       END SUBROUTINE write_eigen_objects

       subroutine full_local_potential(iproc,nproc,orbs,Lzd,iflag,dpbox,xc,potential,pot,comgp)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: iproc,nproc,iflag
         type(orbitals_data),intent(in) :: orbs
         type(local_zone_descriptors),intent(in) :: Lzd
         type(denspot_distribution), intent(in) :: dpbox
         type(xc_info), intent(in) :: xc
         real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), target :: potential
         real(wp), dimension(:), pointer :: pot
         !type(p2pCommsGatherPot),intent(inout), optional:: comgp
         type(p2pComms),intent(inout), optional:: comgp
       END SUBROUTINE full_local_potential

      subroutine free_full_potential(nproc,flag,xc,pot,subname)
         use module_base
         use module_xc
         implicit none
         character(len=*), intent(in) :: subname
         integer, intent(in) :: nproc,flag
         type(xc_info), intent(in) :: xc
         real(wp), dimension(:), pointer :: pot
      END SUBROUTINE free_full_potential

      subroutine calculate_energy_and_gradient(iter,iproc,nproc,GPU,ncong,iscf,&
           energs,wfn,gnrm,gnrm_zero)
        use module_base
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
        use module_base
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
        use module_base
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
        use module_base
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
        type(coulomb_operator), intent(in) :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
        type(denspot_distribution),intent(in),optional :: dpbox
        !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
        real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
        type(p2pComms),intent(inout), optional:: comgp
      END SUBROUTINE local_hamiltonian

      subroutine NK_SIC_potential(lr,orbs,xc,fref,hxh,hyh,hzh,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
        use module_defs, only: gp,wp,dp
        use module_types
         use module_xc
         implicit none
         real(gp), intent(in) :: hxh,hyh,hzh,fref
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         type(coulomb_operator), intent(in) :: pkernel
         type(xc_info), intent(in) :: xc
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
         !real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
         real(wp), intent(inout) :: poti
         real(gp), intent(out) :: eSIC_DC
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,2*orbs%nspin), intent(in), optional :: potandrho
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin), intent(out), optional :: wxdsave
      END SUBROUTINE NK_SIC_potential

      subroutine isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
        use module_defs, only: gp,wp
        use module_types
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

      subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  &
         wfd,psi,orblist)
         use module_base
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
        use module_base
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
         use module_base
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
         use module_base
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
        use module_base
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
        use module_base
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
        use module_base
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
        use module_base
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
        use module_base
        implicit none
        real(wp), dimension(:,:,:,:), pointer :: psiscf
      END SUBROUTINE free_wave_to_isf

      subroutine denspot_communications(iproc,nproc,&
           xc,nspin,geocode,SICapproach,dpbox)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        type(xc_info), intent(in) :: xc
        character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
        character(len=4), intent(in) :: SICapproach
        type(denspot_distribution), intent(inout) :: dpbox
      end subroutine denspot_communications

      subroutine allocateRhoPot(Glr,nspin,atoms,rxyz,denspot)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: nspin
        type(locreg_descriptors), intent(in) :: Glr
        type(atoms_data), intent(in) :: atoms
        real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
        type(DFT_local_fields), intent(inout) :: denspot
      END SUBROUTINE allocateRhoPot

      subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
          fnrm,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,&
          nit_precond,target_function,&
          correction_orthoconstraint,nit_basis,&
          ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
          gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, max_inversion_error, kappa_conv, method_updatekernel,&
          purification_quickreturn, correction_co_contra, cdft, input_frag, ref_frags)
        use module_base
        use module_types
        use module_fragments, only: system_fragment
        use constrained_dft, only: cdft_data
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
        real(kind=8),intent(out) :: trH, fnrm
        real(kind=8),intent(inout) :: trH_old
        type(DFT_PSP_projectors),intent(inout) :: nlpsp
        integer,intent(in) :: scf_mode
        type(localizedDIISParameters),intent(inout) :: ldiis
        type(DFT_wavefunction),target,intent(inout) :: tmb
        type(SIC_data) :: SIC !<parameters for the SIC methods
        type(energy_terms),intent(in) :: energs_base
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
        !these must all be present together
        type(cdft_data),intent(in),optional :: cdft
        type(fragmentInputParameters),optional,intent(in) :: input_frag
        type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
      end subroutine getLocalizedBasis

    subroutine psimix(iproc,nproc,ndim_psi,orbs,comms,diis,hpsit,psit)
      use module_base
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
        energs,nlpsp,SIC,tmb,fnrm,calculate_overlap_matrix,communicate_phi_for_lsumrho,&
        calculate_ham,extra_states,itout,it_scc,it_cdft,order_taylor,max_inversion_error,purification_quickreturn,&
        calculate_KS_residue,calculate_gap,&
        convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,reorder,cdft, updatekernel)
      use module_base
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
      logical,intent(in):: calculate_overlap_matrix, communicate_phi_for_lsumrho, purification_quickreturn
      logical,intent(in) :: calculate_ham, calculate_KS_residue, calculate_gap
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
      use module_base
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
      use module_base
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
      use module_base
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

   subroutine determine_locreg_periodic(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,calculateBounds)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc
      integer, intent(in) :: nlr
      real(gp), intent(in) :: hx,hy,hz
      type(locreg_descriptors), intent(in) :: Glr
      real(gp), dimension(nlr), intent(in) :: locrad
      real(gp), dimension(3,nlr), intent(in) :: cxyz
      type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
      logical,dimension(nlr),intent(in):: calculateBounds
   end subroutine determine_locreg_periodic

    subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr)
      use module_base
      use module_types
      implicit none
      integer,intent(in) :: ilr,nlr
      type(locreg_descriptors),intent(in) :: Glr
      type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr
    end subroutine determine_wfd_periodicity

    subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
     implicit none
     integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
     integer, dimension(nseg), intent(in) :: keyv
     integer, dimension(2,nseg), intent(in) :: keyg
     integer, intent(out) :: nseg_loc,nvctr_loc
     integer, dimension(3),intent(in) :: outofzone
    end subroutine num_segkeys_periodic

    subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
               nseg_loc,nvctr_loc,keygloc,keyglob,keyvloc,keyvglob,outofzone)
      implicit none
      integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
      integer, dimension(nseg), intent(in) :: keyv
      integer, dimension(2,nseg), intent(in) :: keyg
      integer, dimension(3), intent(in) :: outofzone
      integer, dimension(nseg_loc), intent(out) :: keyvloc
      integer, dimension(nseg_loc), intent(out) :: keyvglob
      integer, dimension(2,nseg_loc), intent(out) :: keygloc
      integer, dimension(2,nseg_loc), intent(out) :: keyglob
    end subroutine segkeys_periodic

    subroutine psi_to_locreg2(iproc, ldim, gdim, Llr, Glr, gpsi, lpsi)
      use module_base
      use module_types
      implicit none
      integer,intent(in) :: iproc                  ! process ID
      integer,intent(in) :: ldim          ! dimension of lpsi
      integer,intent(in) :: gdim          ! dimension of gpsi
      type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      real(wp),dimension(gdim),intent(in) :: gpsi       !Wavefunction (compressed format)
      real(wp),dimension(ldim),intent(out) :: lpsi   !Wavefunction in localization region
    end subroutine psi_to_locreg2


    subroutine local_partial_densityLinear(nproc,rsflag,nscatterarr,&
         nrhotot,Lzd,hxh,hyh,hzh,xc,nspin,orbs,mapping,psi,rho)
      use module_base
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

    subroutine global_to_local(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho)
      use module_base
      use module_types
      implicit none
      type(locreg_descriptors),intent(in) :: Llr
      type(locreg_descriptors),intent(in) :: Glr
      integer, intent(in) :: size_rho
      integer, intent(in) :: size_Lrho
      integer, intent(in) :: nspin
      real(wp),dimension(size_rho),intent(in) :: rho
      real(wp),dimension(size_Lrho),intent(out) :: Lrho
     end subroutine global_to_local

     subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,Lzde,comms,&
          psi,hpsi,psit,orthpar,passmat,iscf,Tel,occopt,& !mandatory
          orbse,commse,etol,norbsc_arr) !optional
       use module_base
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

     subroutine updatePotential(nspin,denspot,ehart,eexcu,vexcu)
       use module_base
       use module_types
       implicit none
       ! Calling arguments
       integer, intent(in) :: nspin
       type(DFT_local_fields), intent(inout) :: denspot
       real(8),intent(out):: ehart, eexcu, vexcu
     end subroutine updatePotential

     subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
       integer,dimension(8),intent(out):: comarr
     end subroutine setCommsParameters

     subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, max_inversion_error, npsidim_orbs, &
                orbs, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, lphi, psit_c, psit_f, can_use_transposed, foe_obj)
       use module_base
       use module_types
       use sparsematrix_base, only: sparse_matrix
       use foe_base, only: foe_data
       implicit none
       integer,intent(in) :: iproc,nproc,npsidim_orbs
       integer,intent(inout) :: methTransformOverlap
       real(kind=8),intent(in) :: max_inversion_error
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       type(sparse_matrix),intent(inout):: ovrlp
       type(sparse_matrix),intent(inout):: inv_ovrlp_half
       type(comms_linear),intent(in):: collcom
       type(orthon_data),intent(in):: orthpar
       real(8),dimension(npsidim_orbs), intent(inout) :: lphi
       real(8),dimension(:),pointer:: psit_c, psit_f
       logical,intent(inout):: can_use_transposed
       type(foe_data),intent(in) :: foe_obj
     end subroutine orthonormalizeLocalized

     subroutine optimizeDIIS(iproc, nproc, npsidim, orbs, nspin, lzd, hphi, phi, ldiis, experimental_mode)
       use module_base
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
       use module_base
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

     subroutine initializeRepartitionOrbitals(iproc, nproc, tag, lorbs, llborbs, lzd, comrp)
       use module_base
       use module_types
       use communications_base, only: p2pComms
       implicit none
       integer,intent(in):: iproc, nproc
       integer,intent(inout):: tag
       type(orbitals_data),intent(in):: lorbs, llborbs
       type(local_zone_descriptors),intent(in):: lzd
       !type(p2pCommsRepartition),intent(out):: comrp
       type(p2pComms),intent(out):: comrp
     end subroutine initializeRepartitionOrbitals

     subroutine mixrhopotDIIS(iproc, nproc, n3d, n3p, glr, input, rhopot, rhopotold, mixdiis, alphaMix, ioffset, mixMeth, pnrm, xc)
       use module_base
       use module_types
       use module_xc
       implicit none
       integer,intent(in):: iproc, nproc, n3d, n3p, mixMeth, ioffset
       type(locreg_descriptors),intent(in) :: glr
       type(input_variables),intent(in):: input
       real(8),dimension(max(glr%d%n1i*glr%d%n2i*n3d,1)*input%nspin),intent(in):: rhopotold
       real(8),dimension(max(glr%d%n1i*glr%d%n2i*n3d,1)*input%nspin),intent(out):: rhopot
       type(mixrhopotDIISParameters),intent(inout):: mixdiis
       real(8),intent(in):: alphaMix
       real(8),intent(out):: pnrm
       type(xc_info), intent(in) :: xc
     end subroutine mixrhopotDIIS


     subroutine initializeMixrhopotDIIS(isx, ndimpot, mixdiis)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: isx, ndimpot
       type(mixrhopotDIISParameters),intent(out):: mixdiis
     end subroutine initializeMixrhopotDIIS

     subroutine deallocateMixrhopotDIIS(mixdiis)
       use module_base
       use module_types
       implicit none
       type(mixrhopotDIISParameters),intent(inout):: mixdiis
     end subroutine deallocateMixrhopotDIIS

    subroutine deallocate_local_zone_descriptors(lzd)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(local_zone_descriptors),intent(inout):: lzd
    end subroutine deallocate_local_zone_descriptors

    subroutine deallocate_Lzd_except_Glr(lzd)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(local_zone_descriptors),intent(inout):: lzd
    end subroutine deallocate_Lzd_except_Glr

    subroutine deallocate_orbitals_data(orbs)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(orbitals_data),intent(inout):: orbs
    end subroutine deallocate_orbitals_data

    subroutine deallocate_comms_cubic(comms)
      use module_base
      use module_types
      use communications_base, only: comms_cubic
      !use deallocatePointers
      implicit none
      type(comms_cubic),intent(inout):: comms
    end subroutine deallocate_comms_cubic

    !!subroutine nullify_foe(foe_obj)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(foe_data),intent(out):: foe_obj
    !!end subroutine nullify_foe

    subroutine nullify_sparse_matrix(sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix
      implicit none
      type(sparse_matrix),intent(out):: sparsemat
    end subroutine nullify_sparse_matrix

!    subroutine nullify_comms_linear(collcom)
!      use module_base
!      use module_types
!      implicit none
!      type(comms_linear),intent(inout):: collcom
!    end subroutine nullify_comms_linear

    subroutine nullify_orbitals_data(orbs)
      use module_base
      use module_types
      implicit none
      type(orbitals_data),intent(out):: orbs
    end subroutine nullify_orbitals_data

    subroutine initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, Glr, locregShape, lborbs)
      use module_base
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

    subroutine deallocate_foe(foe_obj, subname)
      use module_base
      use module_types
      use foe_base, only: foe_data
      implicit none
      type(foe_data),intent(inout):: foe_obj
      character(len=*),intent(in):: subname
    end subroutine deallocate_foe

     subroutine initInputguessConfinement(iproc, nproc, at, lzd, orbs, collcom_reference, &
                Glr, input, hx, hy, hz, lin, tmb, rxyz, nscatterarr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc,nproc
       real(gp), intent(in) :: hx, hy, hz
       type(atoms_data),intent(inout) :: at
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs
       type(comms_linear),intent(in):: collcom_reference
       type(locreg_descriptors),intent(in) :: Glr
       type(input_variables), intent(in) ::input
       type(linearInputParameters),intent(in):: lin
       type(DFT_wavefunction),intent(in) :: tmb
       integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       real(gp),dimension(3,at%astruct%nat),intent(in):: rxyz
     end subroutine initInputguessConfinement

      subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, npsidim_orbs, npsidim_comp, orbs, collcom, orthpar, &
                 correction_orthoconstraint, linmat, lphi, lhphi, lagmat, lagmat_, psit_c, psit_f, &
                 hpsit_c, hpsit_f, &
                 can_use_transposed, overlap_calculated, &
                 experimental_mode, calculate_inverse, norder_taylor, max_inversion_error, &
           npsidim_orbs_small, lzd_small, hpsi_noprecond)
        use module_base
        use module_types
        use yaml_output
        implicit none
        integer,intent(in) :: iproc, nproc, npsidim_orbs, npsidim_comp, npsidim_orbs_small
        type(local_zone_descriptors),intent(in) :: lzd, lzd_small
        type(orbitals_Data),intent(inout) :: orbs !temporary inout
        type(comms_linear),intent(in) :: collcom
        type(orthon_data),intent(in) :: orthpar
        integer,intent(in) :: correction_orthoconstraint
        real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(in) :: lphi
        real(kind=8),dimension(max(npsidim_comp,npsidim_orbs)),intent(inout) :: lhphi
        type(sparse_matrix),intent(inout) :: lagmat
        type(matrices),intent(out) :: lagmat_
        real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: hpsit_c
        real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: hpsit_f
        real(kind=8),dimension(:),pointer :: psit_c, psit_f
        logical,intent(inout) :: can_use_transposed, overlap_calculated
        type(linear_matrices),intent(inout) :: linmat ! change to ovrlp and inv_ovrlp, and use inv_ovrlp instead of denskern
        logical,intent(in) :: experimental_mode, calculate_inverse
        integer,intent(inout) :: norder_taylor
        real(kind=8),intent(in) :: max_inversion_error
        real(kind=8),dimension(npsidim_orbs_small),intent(out) :: hpsi_noprecond
      end subroutine orthoconstraintNonorthogonal


      subroutine dsygv_parallel(iproc, nproc, blocksize, nprocMax, comm, itype, jobz, uplo, n, a, lda, b, ldb, w, info)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, blocksize, nprocMax, comm, itype, n, lda, ldb
        integer,intent(out):: info
        character(len=1),intent(in):: jobz, uplo
        real(8),dimension(lda,n),intent(inout):: a
        real(8),dimension(ldb,n),intent(inout):: b
        real(8),dimension(n),intent(out):: w
      end subroutine dsygv_parallel

     subroutine choosePreconditioner2(iproc, nproc, orbs, lr, hx, hy, hz, ncong, hpsi, &
                confpotorder, potentialprefac, iorb, eval_zero)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ncong, iorb, confpotorder
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors), intent(in) :: lr
       type(orbitals_data), intent(in) :: orbs
       real(8),intent(in):: potentialprefac
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
       real(8),intent(in):: eval_zero
     end subroutine choosePreconditioner2


     subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,&
          Lzd,nlpsp,confdatarr,ngatherarr,Lpot,psi,hpsi,paw,&
          energs,SIC,GPU,xc,pkernel,orbsocc,psirocc)
       use module_base
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

       subroutine init_foe(iproc, nproc, input, orbs_KS, foe_obj, reset)
         use module_base
         use module_atoms, only: atomic_structure
         use module_types
         use foe_base, only: foe_data
         implicit none
         integer,intent(in):: iproc, nproc
         type(input_variables),intent(in) :: input
         type(orbitals_data),intent(in):: orbs_KS
         type(foe_data),intent(out):: foe_obj
         logical, intent(in) :: reset
       end subroutine init_foe

      subroutine deallocate_workarrays_quartic_convolutions(work)
        use module_base
        use module_types
        implicit none
        type(workarrays_quartic_convolutions),intent(out):: work
      end subroutine deallocate_workarrays_quartic_convolutions

      subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
                 hgrid, offsetx, offsety, offsetz, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
                 rxyzConf, potentialPrefac, with_kinetic, cprecr, maxdim, &
                 xx_c, xx_f1, xx_f, xy_c, xy_f2, xy_f,  xz_c, xz_f4, xz_f, &
                 aeff0array, beff0array, ceff0array, eeff0array, &
                 aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array, &
                 aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray, &
                 xya_c, xyc_c, xza_c, xzc_c, &
                 yza_c, yzc_c, xya_f, xyb_f, xyc_f, xye_f, &
                 xza_f, xzb_f, xzc_f, xze_f, yza_f, yzb_f, yzc_f, yze_f, &
!                 aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, &
!                 ceff0, ceff1, ceff2, ceff3, eeff0, eeff1, eeff2, eeff3, &
!                 aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2, &
!                 ceff0_2, ceff1_2, ceff2_2, ceff3_2, eeff0_2, eeff1_2, eeff2_2, eeff3_2, &
                 y_c, y_f)
        use module_base
        use module_types
        implicit none
        integer,intent(in) :: iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz, maxdim
        real(gp),intent(in) :: hgrid, potentialPrefac, cprecr
        logical,intent(in) :: with_kinetic
        real(8),dimension(3) :: rxyzConf
        integer,dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
        integer,dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
        integer,dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
        real(wp),dimension(0:n1,0:n2,0:n3),intent(in) :: xx_c
        real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f1
        real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f
        real(wp),dimension(0:n2,0:n1,0:n3),intent(in) :: xy_c
        real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f2
        real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f
        real(wp),dimension(0:n3,0:n1,0:n2),intent(in) :: xz_c
        real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f4
        real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f
        real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0array
        real(wp),dimension(-14:14,0:maxdim),intent(in):: eeff0array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0_2array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0_2array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0_2array
        real(wp),dimension(-14:14,0:maxdim),intent(in):: eeff0_2array
        real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0_2auxarray
        real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0_2auxarray
        real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0_2auxarray
        real(wp),dimension(-17:17,0:maxdim),intent(in):: eeff0_2auxarray
        real(wp),dimension(0:n2,0:n1,0:n3):: xya_c, xyb_c, xyc_c, xye_c
        real(wp),dimension(0:n3,0:n1,0:n2):: xza_c, xzb_c, xzc_c, xze_c, yza_c, yzc_c
        real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xya_f
        real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyb_f
        real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyc_f
        real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xye_f
        real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xza_f
        real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzb_f
        real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzc_f
        real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xze_f
        real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yza_f
        real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzb_f
        real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzc_f
        real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yze_f
!        real(wp),dimension(35):: aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, ceff0, ceff1, ceff2, ceff3
!        real(wp),dimension(29):: eeff0, eeff1, eeff2, eeff3
!        real(wp),dimension(35):: aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2
!        real(wp),dimension(35):: ceff0_2, ceff1_2, ceff2_2, ceff3_2
!        real(wp),dimension(29):: eeff0_2, eeff1_2, eeff2_2, eeff3_2
        real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
        real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
      end subroutine ConvolQuartic4


       subroutine apply_potential_lr(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,&
            psir,pot,epot,&
            confdata,ibyyzz_r,psir_noconf,econf) !optional
         use module_base
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
         use module_base
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
         use module_base
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
         use module_base
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
         use module_base
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
            ref_frags, denspot, locregcenters, inwhichlocreg_old, onwhichatom_old, output_grid)
         use module_base
         use module_types
         use module_fragments
         use communications_base, only: comms_cubic
         implicit none
         integer, intent(in) :: iproc,nproc
         integer, intent(out) :: input_wf_format,lnpsidim_orbs,lnpsidim_comp
         integer, intent(inout) :: inputpsi
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
         logical, intent(in) :: dry_run, dump
         logical, intent(in), optional :: output_grid
       end subroutine system_initialization

       subroutine input_check_psi_id(inputpsi, input_wf_format, dir_output, orbs, lorbs, iproc, nproc, nfrag, frag_dir, ref_frags)
         use module_types
         use module_fragments
         implicit none
         integer, intent(out) :: input_wf_format         !< (out) Format of WF
         integer, intent(inout) :: inputpsi              !< (in) indicate how check input psi, (out) give how to build psi
         integer, intent(in) :: iproc                    !< (in)  id proc
         integer, intent(in) :: nproc                    !< (in)  #proc
         integer, intent(in) :: nfrag                    !< number of fragment directories which need checking
         type(system_fragment), dimension(:), pointer :: ref_frags  !< number of orbitals for each fragment
         character(len=100), dimension(nfrag), intent(in) :: frag_dir !< label for fragment subdirectories (blank if not a fragment calculation)
         character(len = *), intent(in) :: dir_output
         type(orbitals_data), intent(in) :: orbs, lorbs
       end subroutine input_check_psi_id

       subroutine nullify_p2pComms(p2pcomm)
         use module_base
         use module_types
         use communications_base, only: p2pComms
         implicit none
         type(p2pComms),intent(inout):: p2pcomm
       end subroutine nullify_p2pComms

       subroutine extract_potential_for_spectra(iproc,nproc,at,rhod,dpbox,&
            orbs,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
            nlpsp,pkernel,ixc,psi,G,&
            nspin,potshortcut,symObj,GPU,input)
         use module_base
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
         type(coulomb_operator), intent(in) :: pkernel
         integer, intent(in) ::potshortcut
       end subroutine extract_potential_for_spectra

       subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,iscf,alphamix,&
            nlpsp,linflag,unblock_comms,GPU,wfn,&
            energs,rpnrm,xcstr)
         use module_base
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
         use module_base
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
         use module_base
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
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(DFT_wavefunction),intent(inout):: tmb
       end subroutine destroy_new_locregs

       subroutine define_confinement_data(confdatarr,orbs,rxyz,at,hx,hy,hz,&
                  confpotorder,potentialprefac,Lzd,confinementCenter)
         use module_base
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
         use module_base
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
         use module_base
         use module_types
         implicit none
         type(DFT_wavefunction),intent(inout):: wfn
       end subroutine destroy_DFT_wavefunction

       subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspinor
         type(input_variables),intent(in):: input
         type(atomic_structure),intent(in):: astruct
         real(8),dimension(3,astruct%nat),intent(in):: rxyz
         type(orbitals_data),intent(out):: lorbs
       end subroutine init_orbitals_data_for_linear

       subroutine mix_main(iproc, nproc, mix_mode, mixHist, input, glr, alpha_mix, &
                  denspot, mixdiis, rhopotold, pnrm)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, mix_mode, mixHist
         type(input_variables),intent(in):: input
         type(locreg_descriptors),intent(in):: glr
         real(8),intent(in):: alpha_mix
         type(DFT_local_fields),intent(inout):: denspot
         type(mixrhopotDIISParameters),intent(inout):: mixdiis
         real(8),dimension(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout):: rhopotold
         real(8),intent(out):: pnrm
       end subroutine mix_main

       subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
                  ldiis, fnrmOldArr, fnrm_old, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
                  energy_increased, tmb, lhphiold, overlap_calculated, &
                  energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
                  hpsi_small, experimental_mode, calculate_inverse, correction_co_contra, hpsi_noprecond, &
                  norder_taylor, max_inversion_error, method_updatekernel, precond_convol_workarrays, precond_workarrays,&
                  cdft, input_frag, ref_frags)
         use module_base
         use module_types
         use constrained_dft, only: cdft_data
         use module_fragments, only: system_fragment
         implicit none
         integer, intent(in) :: iproc, nproc, it, method_updatekernel
         integer,intent(inout) :: norder_taylor
         real(kind=8),intent(in) :: max_inversion_error
         type(DFT_wavefunction),target,intent(inout):: tmb
         type(localizedDIISParameters),intent(inout) :: ldiis
         real(8),dimension(tmb%orbs%norbp),intent(inout) :: fnrmOldArr
         real(kind=8),intent(inout) :: fnrm_old
         real(8),dimension(tmb%orbs%norbp),intent(inout) :: alpha
         real(8),intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
         real(8),intent(inout):: trHold
         logical,intent(out) :: energy_increased
         real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphiold
         logical,intent(inout):: overlap_calculated
         type(energy_terms),intent(in) :: energs
         real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c) :: hpsit_c
         real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f) :: hpsit_f
         integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
         logical, intent(in) :: experimental_mode, calculate_inverse, correction_co_contra
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(out) :: hpsi_small
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(out) :: hpsi_noprecond
         type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
         type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
         type(cdft_data),intent(in),optional :: cdft
         type(fragmentInputParameters),optional,intent(in) :: input_frag
         type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
       end subroutine calculate_energy_and_gradient_linear

       subroutine improveOrbitals(iproc, nproc, tmb, nspin, ldiis, alpha, gradient, experimental_mode)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspin
         type(DFT_wavefunction),intent(inout):: tmb
         type(localizedDIISParameters),intent(inout):: ldiis
         real(8),dimension(tmb%orbs%norbp),intent(in):: alpha
         real(kind=wp),dimension(max(tmb%npsidim_orbs,tmb%npsidim_comp)),intent(inout) :: gradient
         logical,intent(in) :: experimental_mode
       end subroutine improveOrbitals

       subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
                  lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
                  experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
         use module_base
         use module_types
         implicit none
         integer,intent(in) :: iproc, nproc, it
         integer,intent(inout) :: order_taylor
         real(kind=8),intent(in) :: max_inversion_error
         type(localizedDIISParameters),intent(inout):: ldiis
         type(DFT_wavefunction),target,intent(inout):: tmb
         real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lphiold
         real(8),intent(in):: trH, meanAlpha, alpha_max
         real(8),dimension(tmb%orbs%norbp),intent(inout):: alpha, alphaDIIS
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),intent(inout) :: hpsi_small
         real(kind=8),dimension(tmb%orbs%npsidim_orbs),optional,intent(out) :: psidiff
         logical, intent(in) :: ortho, experimental_mode
         real(kind=8),intent(out) :: trH_ref
         real(kind=8),dimension(tmb%linmat%l%nvctr),intent(out) :: kernel_best
         logical,intent(out) :: complete_reset
       end subroutine hpsitopsi_linear

       subroutine DIISorSD(iproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt, trH_ref, kernel_best, complete_reset)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, it
         real(kind=8),intent(in):: trH
         type(DFT_wavefunction),intent(inout):: tmbopt
         type(localizedDIISParameters),intent(inout):: ldiis
         real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout):: alpha, alphaDIIS
         real(kind=8),dimension(max(tmbopt%npsidim_orbs,tmbopt%npsidim_comp)),intent(out):: lphioldopt
         real(kind=8),intent(out) :: trH_ref
         real(kind=8),dimension(tmbopt%linmat%l%nvctr),intent(out) :: kernel_best
         logical,intent(out) :: complete_reset
       end subroutine DIISorSD

       subroutine psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
            ipotmethod,confdatarr,pot,psi,vpsi,pkernel,xc,alphaSIC,epot_sum,evSIC,vpsi_noconf,econf_sum)
         use module_base
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
         type(coulomb_operator), intent(in) ::  pkernel !< the PSolver kernel which should be associated for the SIC schemes
         real(wp), dimension(orbs%npsidim_orbs), intent(inout),optional :: vpsi_noconf
         real(gp),intent(out),optional :: econf_sum
       end subroutine psi_to_vlocpsi

       subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, at, input, &
                  rxyz, KSwfn, tmb, denspot,nlpsp, ldiis, locreg_increased, lowaccur_converged, locrad)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         real(8),intent(in):: hx, hy, hz
         type(atoms_data),intent(in) :: at
         type(input_variables),intent(in):: input
         real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
         type(DFT_wavefunction),intent(inout):: KSwfn, tmb
         type(DFT_local_fields),intent(inout) :: denspot
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         type(localizedDIISParameters),intent(inout):: ldiis
         logical, intent(out) :: locreg_increased
         logical, intent(in) :: lowaccur_converged
         real(8), dimension(tmb%lzd%nlr), intent(inout) :: locrad
       end subroutine adjust_locregs_and_confinement

       subroutine adjust_DIIS_for_high_accuracy(input, denspot, lowaccur_converged, &
                  ldiis_coeff_hist, ldiis_coeff_changed)
         use module_base
         use module_types
         implicit none
         type(input_variables),intent(in):: input
         type(DFT_local_fields),intent(inout) :: denspot
         logical, intent(in) :: lowaccur_converged
         integer, intent(inout) :: ldiis_coeff_hist
         logical, intent(out) :: ldiis_coeff_changed
       end subroutine adjust_DIIS_for_high_accuracy

       subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, &
                  convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
                  convcrit_dmin, nitdmin, conv_crit_TMB)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: nlr
         type(orbitals_data),intent(in):: lorbs
         type(input_variables),intent(in):: input
         type(atoms_data),intent(in):: at
         integer,dimension(lorbs%norb),intent(in):: onwhichatom
         type(confpot_data),dimension(lorbs%norbp),intent(inout):: confdatarr
         real(kind=8), intent(out) :: convCritMix, alpha_mix, convcrit_dmin, conv_crit_TMB
         logical, intent(in) :: lowaccur_converged
         integer, intent(out) :: nit_scc, mix_hist, nitdmin
         real(kind=8), dimension(nlr), intent(out) :: locrad
         integer, intent(out) :: target_function, nit_basis
       end subroutine set_optimization_variables

       subroutine determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
                  istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
                  weightp_c, weightp_f, nptsp_c, nptsp_f, &
                  norb_per_gridpoint_c, norb_per_gridpoint_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
         real(8),intent(in):: weightp_c, weightp_f
         integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
         integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
       end subroutine determine_num_orbs_per_gridpoint

       subroutine initialize_linear_from_file(iproc,nproc,input_frag,astruct,rxyz,orbs,Lzd,&
              iformat,dir_output,filename,ref_frags,orblist)
         use module_base
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

       subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n_old1, n_old2, n_old3, &
            ns_old1, ns_old2, ns_old3, hgrids_old, lstat, error, onwhichatom, locrad, locregCenter, &
            confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, nat, rxyz_old)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: unitwf
         logical, intent(in) :: formatted
         integer, intent(out) :: iorb_old
         integer, intent(out) :: n_old1, n_old2, n_old3, ns_old1, ns_old2, ns_old3
         real(gp), dimension(3), intent(out) :: hgrids_old
         logical, intent(out) :: lstat
         real(wp), intent(out) :: eval
         real(gp), intent(out) :: locrad
         real(gp), dimension(3), intent(out) :: locregCenter
         character(len =256), intent(out) :: error
         integer, intent(out) :: onwhichatom
         integer, intent(out) :: confPotOrder
         real(gp), intent(out) :: confPotprefac
         ! Optional arguments
         integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
         integer, intent(in), optional :: nat
         real(gp), dimension(:,:), intent(out), optional :: rxyz_old
       end subroutine io_read_descr_linear

        subroutine readmywaves_linear_new(iproc,nproc,dir_output,filename,iformat,at,tmb,rxyz,&
               ref_frags,input_frag,frag_calc,orblist)
          use module_base
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
          logical, intent(in) :: frag_calc
          integer, dimension(tmb%orbs%norb), intent(in), optional :: orblist
        end subroutine readmywaves_linear_new

        subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi)
          use module_base
          implicit none
          integer,intent(in):: npsidim
          real(8),dimension(:),pointer,intent(out):: lphi, lhphi
          character(len=*),intent(in):: subname
        end subroutine allocate_auxiliary_basis_function

        subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi)
          use module_base
          implicit none
          real(8),dimension(:),pointer:: lphi, lhphi
          character(len=*),intent(in):: subname
        end subroutine deallocate_auxiliary_basis_function

        subroutine update_ldiis_arrays(tmb, subname, ldiis)
          use module_base
          use module_types
          implicit none
          type(DFT_wavefunction),intent(in):: tmb
          character(len=*),intent(in):: subname
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine update_ldiis_arrays

        subroutine copy_local_zone_descriptors(lzd_in, lzd_out, subname)
          use module_base
          use module_types
          implicit none
          type(local_zone_descriptors),intent(in):: lzd_in
          type(local_zone_descriptors),intent(inout):: lzd_out
          character(len=*),intent(in):: subname
        end subroutine copy_local_zone_descriptors

        subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, &
            & lstat, error, nat, rxyz_old)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: unitwf
         logical, intent(in) :: formatted
         integer, intent(out) :: norb_old, ntmb_old
         logical, intent(out) :: lstat
         character(len =256), intent(out) :: error
         ! Optional arguments
         integer, intent(in), optional :: nat
         real(gp), dimension(:,:), intent(out), optional :: rxyz_old
        end subroutine io_read_descr_coeff

        subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,ntmb,norb_old,coeff,eval,nat,rxyz_old)
          use module_base
          use module_types
          implicit none
          logical, intent(in) :: useFormattedInput
          integer, intent(in) :: unitwf,iproc,ntmb
          integer, intent(out) :: norb_old
          real(wp), dimension(ntmb,ntmb), intent(out) :: coeff
          real(wp), dimension(ntmb), intent(out) :: eval
          integer, optional, intent(in) :: nat
          real(gp), dimension(:,:), optional, intent(out) :: rxyz_old
        end subroutine read_coeff_minbasis

        subroutine local_potential_dimensions(iproc,Lzd,orbs,xc,ndimfirstproc)
          use module_base
          use module_types
          use module_xc
          implicit none
          integer, intent(in) :: iproc, ndimfirstproc
          type(local_zone_descriptors), intent(inout) :: Lzd
          type(orbitals_data), intent(inout) :: orbs
          type(xc_info), intent(in) :: xc
        end subroutine local_potential_dimensions

        subroutine initialize_DIIS_coeff(isx, ldiis)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: isx
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine initialize_DIIS_coeff

        subroutine allocate_DIIS_coeff(tmb, ldiis)
          use module_base
          use module_types
          implicit none
          type(DFT_wavefunction),intent(in):: tmb
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine allocate_DIIS_coeff

        subroutine initialize_DFT_local_fields(denspot, ixc, nspden)
          use module_base
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
          integer, intent(in) :: ixc, nspden
        end subroutine initialize_DFT_local_fields

        subroutine allocate_diis_objects(idsx,alphadiis,npsidim,nkptsp,nspinor,diis)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: idsx,npsidim,nkptsp,nspinor !n(m)
          real(gp), intent(in) :: alphadiis
          type(diis_objects), intent(inout) :: diis
        end subroutine allocate_diis_objects

        subroutine check_communications(iproc,nproc,orbs,lzd,comms)
          use module_base
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: lzd
          type(comms_cubic), intent(in) :: comms
        end subroutine check_communications

        subroutine nonlocal_forces(lr,hx,hy,hz,at,rxyz,&
             orbs,nlpsp,wfd,psi,fsep,refill,strten)
          use module_base
          use module_types
          implicit none
          !Arguments-------------
          type(atoms_data), intent(in) :: at
          type(wavefunctions_descriptors), intent(in) :: wfd
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          logical, intent(in) :: refill
          real(gp), intent(in) :: hx,hy,hz
          type(locreg_descriptors) :: lr
          type(orbitals_data), intent(in) :: orbs
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor), intent(in) :: psi
          real(gp), dimension(3,at%astruct%nat), intent(inout) :: fsep
          real(gp), dimension(6), intent(out) :: strten
        end subroutine nonlocal_forces

        subroutine local_forces(iproc,at,rxyz,hxh,hyh,hzh,&
             n1,n2,n3,n3pi,i3s,n1i,n2i,rho,pot,floc,locstrten,charge)
          use module_base
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
          use module_base
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
          use module_base
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
          use module_base
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
          use module_base
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
          use module_base
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
          use module_base
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
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in):: iproc, nproc
          logical, intent(in) :: isKernel
          type(orbitals_data),intent(in):: orbs, orbs_tmb
          type(sparse_matrix), intent(inout) :: denskern
          real(kind=8),dimension(denskern%nfvctr,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
          type(matrices), intent(out) :: denskern_
          logical,intent(in),optional :: keep_uncompressed_ !< keep the uncompressed kernel in denskern_%matrix (requires that this array is already allocated outside of the routine)
        end subroutine calculate_density_kernel

        subroutine reconstruct_kernel(iproc, nproc, inversion_method, &
                   blocksize_dsyev, blocksize_pdgemm, orbs, tmb, overlap_calculated)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, blocksize_dsyev, blocksize_pdgemm, inversion_method
          type(orbitals_data),intent(in):: orbs
          type(DFT_wavefunction),intent(inout):: tmb
          logical,intent(inout):: overlap_calculated
        end subroutine reconstruct_kernel

        subroutine reorthonormalize_coeff(iproc, nproc, norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, basis_orbs, &
                   basis_overlap, KS_overlap, basis_overlap_mat, coeff, orbs)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix, matrices
          implicit none
          integer, intent(in) :: iproc, nproc, norb
          integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
          type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
          type(sparse_matrix),intent(inout) :: basis_overlap
          type(sparse_matrix),dimension(basis_overlap%nspin),intent(inout) :: KS_overlap
          type(matrices),intent(inout) :: basis_overlap_mat
          real(kind=8),dimension(basis_overlap%nfvctr,norb),intent(inout) :: coeff
          type(orbitals_data), intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
        end subroutine reorthonormalize_coeff

        subroutine pulay_correction(iproc, nproc, orbs, at, rxyz, nlpsp, SIC, denspot, GPU, tmb, fpulay)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(atoms_data),intent(in):: at
          real(8),dimension(at%astruct%nat),intent(in):: rxyz
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(SIC_data),intent(in):: SIC
          type(DFT_local_fields), intent(inout) :: denspot
          type(GPU_pointers),intent(inout):: GPU
          type(DFT_wavefunction),intent(inout):: tmb
          real(8),dimension(3,at%astruct%nat),intent(out):: fpulay
        end subroutine pulay_correction

        subroutine create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot,nlpsp, input, at, rxyz, lowaccur_converged)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(DFT_Wavefunction),intent(inout):: KSwfn, tmb
          type(DFT_local_fields),intent(in):: denspot
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(input_variables),intent(in):: input
          type(atoms_data),intent(in):: at
          real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
          logical,intent(in):: lowaccur_converged
        end subroutine create_large_tmbs


        subroutine solvePrecondEquation(iproc,nproc,lr,ncplx,ncong,cprecr,&
             hx,hy,hz,kx,ky,kz,x,  rxyzParab, orbs, potentialPrefac, confPotOrder,&
             work_conv, w)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc,ncong,ncplx,confPotOrder
          real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
          type(locreg_descriptors), intent(in) :: lr
          real(wp), intent(inout) :: x
          real(8),dimension(3),intent(in):: rxyzParab
          type(orbitals_data), intent(in):: orbs
          real(8):: potentialPrefac
          type(workarrays_quartic_convolutions),intent(inout):: work_conv !< workarrays for the convolutions
          type(workarr_precond),intent(inout) :: w !< workarrays
        end subroutine solvePrecondEquation

        subroutine init_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work)
          use module_base
          use module_types
          implicit none
          integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
          logical,intent(in):: with_confpot
          type(workarrays_quartic_convolutions),intent(inout):: work
        end subroutine init_local_work_arrays

        subroutine psi_to_kinpsi(iproc,npsidim_orbs,orbs,lzd,psi,hpsi,ekin_sum)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,npsidim_orbs
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: Lzd
          real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
          real(gp), intent(out) :: ekin_sum
          real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        end subroutine psi_to_kinpsi

        subroutine copy_old_supportfunctions(iproc,orbs,lzd,phi,lzd_old,phi_old)
          use module_base
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
          use module_base
          use module_types
          use module_fragments
          use constrained_dft
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

        subroutine copy_old_coefficients(norb_tmb, coeff, coeff_old)
          use module_base
          implicit none
          integer,intent(in):: norb_tmb
          real(8),dimension(:,:),pointer:: coeff, coeff_old
        end subroutine copy_old_coefficients

        subroutine copy_old_inwhichlocreg(norb_tmb, inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old)
          use module_base
          implicit none
          integer,intent(in):: norb_tmb
          integer,dimension(:),pointer:: inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old
        end subroutine copy_old_inwhichlocreg

        subroutine reformat_supportfunctions(iproc,nproc,at,rxyz_old,rxyz,add_derivatives,tmb,ndim_old,lzd_old,&
               frag_trans,psi_old,input_dir,input_frag,ref_frags,phi_array_old)
          use module_base
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
        end subroutine reformat_supportfunctions

        subroutine reformat_one_supportfunction(llr,llr_old,geocode,hgrids_old,n_old,psigold,&
             hgrids,n,centre_old,centre_new,da,frag_trans,psi,psirold)
          use module_base
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

        subroutine get_derivative_supportfunctions(ndim, hgrid, lzd, lorbs, phi, phid)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: ndim
          real(kind=8),intent(in) :: hgrid
          type(local_zone_descriptors),intent(in) :: lzd
          type(orbitals_data),intent(in) :: lorbs
          real(kind=8),dimension(lorbs%npsidim_orbs),intent(in) :: phi !< Basis functions
          real(kind=8),dimension(3*lorbs%npsidim_orbs),intent(inout) :: phid  !< Derivative basis functions
        end subroutine get_derivative_supportfunctions

        subroutine normalize_transposed(iproc, nproc, orbs, nspin, collcom, psit_c, psit_f, norm)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, nspin
          type(orbitals_data),intent(in):: orbs
          type(comms_linear),intent(in):: collcom
          real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
          real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
          real(8),dimension(orbs%norb),intent(out):: norm
        end subroutine normalize_transposed


        subroutine determine_locregSphere_parallel(iproc,nproc,nlr,hx,hy,hz,astruct,orbs,Glr,Llr,calculateBounds)!,outofzone)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc
          integer, intent(in) :: nlr
          real(gp), intent(in) :: hx,hy,hz
          type(atomic_structure),intent(in) :: astruct
          type(orbitals_data),intent(in) :: orbs
          type(locreg_descriptors), intent(in) :: Glr
          type(locreg_descriptors), dimension(nlr), intent(inout) :: Llr
          logical,dimension(nlr),intent(in) :: calculateBounds
        end subroutine determine_locregSphere_parallel

        subroutine communicate_basis_for_density_collective(iproc, nproc, lzd, npsidim, orbs, lphi, collcom_sr)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc, npsidim
          type(local_zone_descriptors),intent(in) :: lzd
          type(orbitals_data),intent(in) :: orbs
          real(kind=8),dimension(npsidim),intent(in) :: lphi
          type(comms_linear),intent(inout) :: collcom_sr
        end subroutine communicate_basis_for_density_collective

        subroutine sumrho_for_TMBs(iproc, nproc, hx, hy, hz, collcom_sr, denskern, denskern_, ndimrho, rho, rho_negative, &
                   print_results)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, ndimrho
          real(kind=8),intent(in) :: hx, hy, hz
          type(comms_linear),intent(in) :: collcom_sr
          type(sparse_matrix),intent(in) :: denskern
          type(matrices),intent(in) :: denskern_
          real(kind=8),dimension(ndimrho),intent(out) :: rho
          logical,intent(out) :: rho_negative
          logical,intent(in),optional :: print_results
        end subroutine sumrho_for_TMBs

        subroutine foe(iproc, nproc, tmprtr, &
                   ebs, itout, it_scc, order_taylor, max_inversion_error, purification_quickreturn, &
                   calculate_minusonehalf, foe_verbosity, &
                   accuracy_level, tmb, foe_obj)
          use module_base
          use module_types
          use foe_base, only: foe_data
          implicit none
          integer,intent(in) :: iproc, nproc, itout, it_scc
          integer,intent(inout) :: order_taylor
          real(kind=8),intent(in) :: max_inversion_error
          real(kind=8),intent(in) :: tmprtr
          real(kind=8),intent(out) :: ebs
          logical,intent(in) :: purification_quickreturn
          logical,intent(in) :: calculate_minusonehalf
          integer :: foe_verbosity
          integer,intent(in) :: accuracy_level
          type(DFT_wavefunction),intent(inout) :: tmb
          type(foe_data),intent(inout) :: foe_obj
        end subroutine foe

        subroutine kswfn_init_comm(wfn, dpbox, iproc, nproc, nspin, imethod_overlap)
          use module_types
          implicit none
          integer, intent(in) :: iproc, nproc, nspin, imethod_overlap
          type(DFT_wavefunction), intent(inout) :: wfn
          type(denspot_distribution), intent(in) :: dpbox
        end subroutine kswfn_init_comm


        subroutine nonlocal_forces_linear(iproc,nproc,npsidim_orbs,lr,hx,hy,hz,at,rxyz,&
             orbs,nlpsp,lzd,phi,denskern,denskern_mat,fsep,refill,strten)
          use module_base
          use module_types
          implicit none
          type(atoms_data), intent(in) :: at
          type(local_zone_descriptors), intent(in) :: lzd
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          logical, intent(in) :: refill
          integer, intent(in) :: iproc, nproc, npsidim_orbs
          real(gp), intent(in) :: hx,hy,hz
          type(locreg_descriptors) :: lr
          type(orbitals_data), intent(in) :: orbs
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          real(wp), dimension(npsidim_orbs), intent(in) :: phi
          type(sparse_matrix),intent(in) :: denskern
          type(matrices),intent(inout) :: denskern_mat
          real(gp), dimension(3,at%astruct%nat), intent(inout) :: fsep
          real(gp), dimension(6), intent(out) :: strten
        end subroutine nonlocal_forces_linear

        subroutine calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
                   psit_c1, psit_c2, psit_f1, psit_f2, smat, ovrlp)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc
          type(orbitals_data),intent(in) :: orbs
          type(comms_linear),intent(in) :: collcom
          real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c1, psit_c2
          real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f1, psit_f2
          type(sparse_matrix),intent(inout) :: smat
          type(matrices),intent(inout) :: ovrlp
        end subroutine calculate_overlap_transposed

        subroutine build_linear_combination_transposed(collcom, sparsemat, mat, psitwork_c, psitwork_f, &
             reset, psit_c, psit_f, iproc)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none

          ! Calling arguments
          type(sparse_matrix),intent(in) :: sparsemat
          type(matrices),intent(in) :: mat
          type(comms_linear),intent(in) :: collcom
          real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
          real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
          logical,intent(in) :: reset
          real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: psit_c
          real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: psit_f
          integer, intent(in) :: iproc
        end subroutine build_linear_combination_transposed

        subroutine axpy_kernel_vectors(norbp, norb, nout, onedimindices, a, x, y)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: norbp, norb, nout
          integer,dimension(4,nout),intent(in) :: onedimindices
          real(kind=8),intent(in) :: a
          real(kind=8),dimension(norb,norbp),intent(in) :: x
          real(kind=8),dimension(norb,norbp),intent(inout) :: y
        end subroutine axpy_kernel_vectors

        subroutine axbyz_kernel_vectors(norbp, norb, nout, onedimindices, a, x, b, y, z)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: norbp, norb, nout
          integer,dimension(4,nout),intent(in) :: onedimindices
          real(8),intent(in) :: a, b
          real(kind=8),dimension(norb,norbp),intent(in) :: x, y
          real(kind=8),dimension(norb,norbp),intent(out) :: z
        end subroutine axbyz_kernel_vectors

        subroutine copy_kernel_vectors(norbp, norb, nout, onedimindices, a, b)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: norbp, norb, nout
          integer,dimension(4,nout),intent(in) :: onedimindices
          real(kind=8),dimension(norb,norbp),intent(in) :: a
          real(kind=8),dimension(norb,norbp),intent(out) :: b
        end subroutine copy_kernel_vectors

        subroutine chebyshev_clean(iproc, nproc, npl, cc, norb, norbp, isorb, kernel, ham_compr, &
                   invovrlp_compr, calculate_SHS, nsize_polynomial, ncalc, fermi, penalty_ev, chebyshev_polynomials, &
                   emergency_stop)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, npl, nsize_polynomial, norb, norbp, isorb, ncalc
          real(8),dimension(npl,3,ncalc),intent(in) :: cc
          type(sparse_matrix), intent(in) :: kernel
          real(kind=8),dimension(kernel%nvctrp_tg),intent(in) :: ham_compr
          real(kind=8),dimension(kernel%nvctr),intent(in) :: invovrlp_compr
          logical,intent(in) :: calculate_SHS
          real(kind=8),dimension(kernel%nfvctr,kernel%smmm%nfvctrp,ncalc),intent(out) :: fermi
          real(kind=8),dimension(kernel%nfvctr,kernel%smmm%nfvctrp,2),intent(out) :: penalty_ev
          real(kind=8),dimension(nsize_polynomial,npl),intent(out) :: chebyshev_polynomials
          logical,intent(out) :: emergency_stop
        end subroutine chebyshev_clean


        subroutine set_variables_for_hybrid(nlr, input, at, orbs, lowaccur_converged, confdatarr, &
                   target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix, &
                   conv_crit_TMB)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: nlr
          type(input_variables),intent(in) :: input
          type(atoms_data),intent(in) :: at
          type(orbitals_data),intent(in) :: orbs
          logical,intent(out) :: lowaccur_converged
          type(confpot_data),dimension(orbs%norbp),intent(inout) :: confdatarr
          integer,intent(out) :: target_function, nit_basis, nit_scc, mix_hist
          real(kind=8),dimension(nlr),intent(out) :: locrad
          real(kind=8),intent(out) :: alpha_mix, convCritMix, conv_crit_TMB
        end subroutine set_variables_for_hybrid

        subroutine locreg_bounds(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds)
          use locregs, only: wavefunctions_descriptors, convolutions_bounds
          implicit none
          integer, intent(in) :: n1,n2,n3
          integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
          type(wavefunctions_descriptors), intent(in) :: wfd
          type(convolutions_bounds), intent(out) :: bounds
        end subroutine locreg_bounds

        subroutine wfd_to_logrids(n1,n2,n3,wfd,logrid_c,logrid_f)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: n1,n2,n3
          type(wavefunctions_descriptors), intent(in) :: wfd
          logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid_c,logrid_f
        end subroutine wfd_to_logrids

        subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
           implicit none
           integer, intent(in) :: n1,n2,n3
           logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid
           integer, dimension(2,0:n2,0:n3), intent(out) :: ibyz
           integer, dimension(2,0:n1,0:n3), intent(out) :: ibxz
           integer, dimension(2,0:n1,0:n2), intent(out) :: ibxy
        end subroutine make_bounds

        subroutine make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
             ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
             ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
          use module_base
          implicit none
          integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
          integer :: i1,i2,i3,i_stat,i_all !n(c) m1,m2,m3
          integer,intent(in):: ibyz_c(2,0:n2,0:n3),ibxy_c(2,0:n1,0:n2)
          integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)
          integer,intent(inout):: ibzzx_c(2,-14:2*n3+16,0:n1)
          integer,intent(out):: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)
          integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
          integer,intent(inout):: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
          integer,intent(out):: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
          integer,intent(out):: ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
          integer,intent(out):: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)
          integer,intent(inout):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
          integer,intent(out):: ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
          integer,intent(out):: ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
          character(len=*), parameter :: subname=' make_all_ib'
          logical,allocatable:: logrid_big(:)
          integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
        end subroutine make_all_ib

        subroutine make_ib_inv(logrid_big,ibxy,ibzzx,ibyyzz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
          implicit none
          integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
          integer,intent(in):: ibxy(2,nfl1:nfu1,nfl2:nfu2)
          integer,intent(inout):: ibzzx(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
          integer,intent(out):: ibyyzz(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
          logical, intent(inout) :: logrid_big(nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array
          integer :: nt
        end subroutine make_ib_inv

        subroutine ib_to_logrid_inv(ib,logrid,nfl,nfu,ndat)
          implicit none
          integer, intent(in) :: ndat,nfl,nfu
          integer, intent(in) :: ib(2,ndat)! input
          logical, intent(out) :: logrid(-14+2*nfl:2*nfu+16,ndat)! output
        end subroutine ib_to_logrid_inv

        subroutine ib_from_logrid_inv(ib,logrid,ml1,mu1,ndat)
          implicit none
          integer, intent(in) :: ml1,mu1,ndat
          integer, intent(out) :: ib(2,ndat)
          logical, intent(in) :: logrid(ndat,ml1:mu1)
        end subroutine ib_from_logrid_inv

        subroutine squares(ib,n2,n3)
          implicit none
          integer,intent(in)::n2,n3
          integer,intent(inout)::ib(2,0:n2,0:n3)
        end subroutine squares

        subroutine make_ib_c(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,n3)
          implicit none
          integer nt,n1,n2,n3
          integer ibyz(2,0:n2,0:n3)! input
          integer ibzxx(2,0:n3,-14:2*n1+16)!output
          integer ibxxyy(2,-14:2*n1+16,-14:2*n2+16)!output
          logical logrid_big(0:n3,-14:2*n1+16,-14:2*n2+16)! work array
        end subroutine make_ib_c

        subroutine ib_to_logrid_rot(ib,logrid,nfl,nfu,ndat)
          implicit none
          integer ndat,nfl,nfu,l,i
          integer ib(2,ndat)! input
          logical logrid(ndat,-14+2*nfl:2*nfu+16)! output
        end subroutine ib_to_logrid_rot

        subroutine ib_from_logrid(ib,logrid,ml1,mu1,ndat)
          implicit none
          integer i,i1
          integer ml1,mu1,ndat
          integer ib(2,ndat)
          logical logrid(ml1:mu1,ndat)
        end subroutine ib_from_logrid

        subroutine make_ib(logrid_big,ibyz,ibzxx,ibxxyy,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
          implicit none
          integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
          integer ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
          integer ibzxx( 2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)!output
          integer ibxxyy(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)!output
          logical logrid_big(           nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array
        end subroutine make_ib

        subroutine squares_1d(ib,nfl2,nfu2,nfl3,nfu3)
          implicit none
          integer,intent(in) :: nfl2,nfu2,nfl3,nfu3
          integer,intent(inout) :: ib(2,nfl2:nfu2,nfl3:nfu3)
        end subroutine squares_1d

        subroutine orthonormalize_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
                   orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
                   lphi, psit_c, psit_f, can_use_transposed)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
          type(orbitals_data),intent(in) :: orbs
          type(atoms_data),intent(in) :: at
          integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
          type(local_zone_descriptors),intent(in) :: lzd
          type(sparse_matrix),intent(inout) :: ovrlp
          type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
          type(comms_linear),intent(in) :: collcom
          type(orthon_data),intent(in) :: orthpar
          real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
          real(kind=8),dimension(:),pointer :: psit_c, psit_f
          logical,intent(inout) :: can_use_transposed
        end subroutine orthonormalize_subset

        subroutine gramschmidt_subset(iproc, nproc, methTransformOverlap, npsidim_orbs, &
                   orbs, at, minorbs_type, maxorbs_type, lzd, ovrlp, inv_ovrlp_half, collcom, orthpar, &
                   lphi, psit_c, psit_f, can_use_transposed)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc,nproc,methTransformOverlap,npsidim_orbs
          type(orbitals_data),intent(in) :: orbs
          type(atoms_data),intent(in) :: at
          integer,dimension(at%astruct%ntypes),intent(in) :: minorbs_type, maxorbs_type
          type(local_zone_descriptors),intent(in) :: lzd
          type(sparse_matrix),intent(inout) :: ovrlp
          type(sparse_matrix),intent(inout) :: inv_ovrlp_half ! technically inv_ovrlp structure, but same pattern
          type(comms_linear),intent(in) :: collcom
          type(orthon_data),intent(in) :: orthpar
          real(kind=8),dimension(npsidim_orbs), intent(inout) :: lphi
          real(kind=8),dimension(:),pointer :: psit_c, psit_f
          logical,intent(inout) :: can_use_transposed
        end subroutine gramschmidt_subset


        subroutine overlapPowerGeneral(iproc, nproc, iorder, ncalc, power, blocksize, imode, &
                   ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, inv_ovrlp_mat, check_accur, &
                   max_error, mean_error, nspinx)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix, SPARSE_FULL, DENSE_PARALLEL, DENSE_FULL, SPARSEMM_SEQ
          use yaml_output
          integer,intent(in) :: iproc, nproc, iorder, blocksize, ncalc
          integer,dimension(ncalc),intent(in) :: power
          integer,intent(in) :: imode
          type(sparse_matrix),intent(inout) :: ovrlp_smat, inv_ovrlp_smat
          type(matrices),intent(inout) :: ovrlp_mat
          type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp_mat
          logical,intent(in) :: check_accur
          real(kind=8),intent(out),optional :: max_error, mean_error
          integer,intent(in),optional :: nspinx !< overwrite the default spin value
        end subroutine overlapPowerGeneral


        !!subroutine overlap_plus_minus_one_half_exact(nproc,norb,blocksize,plusminus,inv_ovrlp_half,smat)
        !!  use module_base
        !!  use module_types
        !!  implicit none
        !!  integer,intent(in) :: nproc,norb,blocksize
        !!  real(kind=8),dimension(:,:),pointer :: inv_ovrlp_half
        !!  logical, intent(in) :: plusminus
        !!  type(sparse_matrix),intent(in) :: smat
        !!end subroutine overlap_plus_minus_one_half_exact

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
          use module_base
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
          use module_types
          implicit none
          type(input_variables), pointer :: in
        end subroutine inputs_new

        subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
                   collcom_sr, sparsemat)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc
          type(orbitals_data),intent(in) :: orbs
          type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
          type(sparse_matrix), intent(inout) :: sparsemat
        end subroutine init_matrixindex_in_compressed_fortransposed

        subroutine compress_polynomial_vector(iproc, nproc, nsize_polynomial, norb, norbp, isorb, &
                   fermi, vector, vector_compressed)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp, isorb
          type(sparse_matrix),intent(in) :: fermi
          real(kind=8),dimension(norb,norbp),intent(in) :: vector
          real(kind=8),dimension(nsize_polynomial),intent(out) :: vector_compressed
        end subroutine compress_polynomial_vector

        subroutine uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                   norb, norbp, isorb, fermi, vector_compressed, vector)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp, isorb
          type(sparse_matrix),intent(in) :: fermi
          real(kind=8),dimension(nsize_polynomial),intent(in) :: vector_compressed
          real(kind=8),dimension(norb,norbp),intent(out) :: vector
        end subroutine uncompress_polynomial_vector

        subroutine check_communication_sumrho(iproc, nproc, orbs, lzd, collcom_sr, denspot, denskern, denskern_, check_sumrho)
          use module_base
          use module_types
          use yaml_output
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc
          type(local_zone_descriptors),intent(in) :: lzd
          type(orbitals_data),intent(in) :: orbs
          type(comms_linear),intent(inout) :: collcom_sr
          type(DFT_local_fields),intent(in) :: denspot
          type(sparse_matrix),intent(inout) :: denskern
          type(matrices),intent(inout) :: denskern_
          integer,intent(in) :: check_sumrho
        end subroutine check_communication_sumrho

        subroutine purify_kernel(iproc, nproc, tmb, overlap_calculated, it_shift, it_opt, order_taylor, &
                   max_inversion_error, purification_quickreturn, ispin)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc
          integer,intent(inout) :: order_taylor
          real(kind=8),intent(in) :: max_inversion_error
          type(DFT_wavefunction),intent(inout):: tmb
          logical,intent(inout):: overlap_calculated
          integer,intent(in) :: it_shift, it_opt
          logical,intent(in) :: purification_quickreturn
          integer,intent(in) :: ispin
        end subroutine purify_kernel

        subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, &
               sd_fit_curve, factor, itout, it_scc, it_cdft, order_taylor, max_inversion_error, reorder, num_extra)
          use module_base
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
          use module_base
          use module_types
          implicit none

          ! Calling arguments
          integer, intent(in) :: iproc, nproc, num_extra
          type(dft_wavefunction), intent(inout) :: tmb
          type(orbitals_data), intent(in) :: ksorbs
          real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f
        end subroutine calculate_residue_ks

        subroutine write_energies(iter,iscf,energs,gnrm,gnrm_zero,comment,only_energies)
          use module_base
          use module_types
          use yaml_output
          implicit none
          integer, intent(in) :: iter,iscf
          type(energy_terms), intent(in) :: energs
          real(gp), intent(in) :: gnrm,gnrm_zero
          character(len=*), intent(in) :: comment
          logical,intent(in),optional :: only_energies
        end subroutine write_energies

        subroutine build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
                 energs, nlpsp, input, order_taylor, &
                 energy, energyDiff, energyold)
          use module_base
          use module_types
          implicit none
          integer:: iproc, nproc
          type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
          type(atoms_data), intent(in) :: at
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          type(DFT_local_fields), intent(inout) :: denspot
          type(GPU_pointers), intent(inout) :: GPU
          type(energy_terms),intent(inout) :: energs
          type(DFT_PSP_projectors), intent(inout) :: nlpsp
          type(input_variables),intent(in) :: input
          integer,intent(inout) :: order_taylor
          real(kind=8),intent(out) :: energy, energyDiff
          real(kind=8), intent(inout) :: energyold
        end subroutine build_ks_orbitals

!TEMPORARY, to be cleaned/removed
subroutine build_ks_orbitals_laura_tmp(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
           energs, nlpsp, input, order_taylor, &
           energy, energyDiff, energyold, npsidim_global, phiwork_global)
  use module_base
  use module_types
  use communications_base, only: comms_cubic
  use communications_init, only: orbitals_communicators
  use communications, only: transpose_v, untranspose_v
  use sparsematrix_base, only: sparse_matrix
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer:: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(input_variables),intent(in) :: input
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(out) :: energy, energyDiff
  real(kind=8), intent(inout) :: energyold
integer, intent(in) :: npsidim_global
real(kind=8),dimension(:),pointer :: phiwork_global
end subroutine build_ks_orbitals_laura_tmp

        subroutine small_to_large_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
               orbs, phismall, philarge, to_global)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
          type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
          type(orbitals_data),intent(in) :: orbs
          real(kind=8),dimension(npsidim_orbs_small),intent(in) :: phismall
          real(kind=8),dimension(npsidim_orbs_large),intent(out) :: philarge
          logical,intent(in),optional :: to_global
        end subroutine small_to_large_locreg

        subroutine get_KS_residue(iproc, nproc, tmb, KSorbs, hpsit_c, hpsit_f, KSres)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc
          type(DFT_wavefunction) :: tmb
          type(orbitals_data),intent(in) :: KSorbs
          real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(in) :: hpsit_c
          real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(in) :: hpsit_f
          real(kind=8),intent(out) :: KSres
        end subroutine get_KS_residue

        subroutine applyprojectorsonthefly(iproc,orbs,at,lr,&
             rxyz,hx,hy,hz,wfd,nlpsp,psi,hpsi,eproj_sum,&
             paw)
          use module_base
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

        subroutine pulay_correction_new(iproc, nproc, tmb, orbs, at, fpulay)
          use module_base
          use module_types
          use yaml_output
          implicit none
          integer,intent(in) :: iproc, nproc
          type(DFT_wavefunction),intent(inout) :: tmb
          type(orbitals_data),intent(in) :: orbs
          type(atoms_data),intent(in) :: at
          real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay
        end subroutine pulay_correction_new

        subroutine toglobal_and_transpose(iproc,nproc,orbs,Lzd,comms,psi,&
             work,outadd) !optional
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: Lzd
          type(comms_cubic), intent(in) :: comms
          real(wp), dimension(:), pointer :: psi
          real(wp), dimension(:), pointer, optional :: work
          real(wp), dimension(*), intent(out), optional :: outadd
        end subroutine

        subroutine increase_FOE_cutoff(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, foe_obj, init)
          use module_base
          use module_types
          use foe_base, only: foe_data
          implicit none
          integer,intent(in) :: iproc, nproc
          type(local_zone_descriptors),intent(in) :: lzd
          type(atomic_structure),intent(in) :: astruct
          type(input_variables),intent(in) :: input
          type(orbitals_data),intent(in) :: orbs_KS, orbs
          type(foe_data),intent(out) :: foe_obj
          logical,intent(in) :: init
        end subroutine increase_FOE_cutoff

        subroutine loewdin_charge_analysis(iproc,tmb,atoms,denspot,&
                   calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap)
          use module_base
          use module_types
          use communications, only: transpose_localized
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc
          type(dft_wavefunction),intent(inout) :: tmb
          type(atoms_data),intent(in) :: atoms
          type(DFT_local_fields), intent(inout) :: denspot
          logical,intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
          integer,intent(in) :: meth_overlap
        end subroutine loewdin_charge_analysis
        
        subroutine astruct_set_n_atoms(astruct, nat)
          use module_base
          use module_atoms, only: atomic_structure
          implicit none
          type(atomic_structure), intent(inout) :: astruct
          integer, intent(in) :: nat
        end subroutine astruct_set_n_atoms

        subroutine astruct_set_n_types(astruct, ntypes)
          use module_base
          use module_atoms, only: atomic_structure
          implicit none
          type(atomic_structure), intent(inout) :: astruct
          integer, intent(in) :: ntypes
        end subroutine astruct_set_n_types

        subroutine allocate_atoms_nat(atoms)
          use module_base
          use module_atoms, only: atoms_data
          use ao_inguess, only : aoig_data_null
          implicit none
          type(atoms_data), intent(inout) :: atoms
        end subroutine allocate_atoms_nat

        subroutine allocate_atoms_ntypes(atoms)
          use module_base
          use module_atoms, only: atoms_data
          implicit none
          type(atoms_data), intent(inout) :: atoms
        end subroutine allocate_atoms_ntypes

        subroutine astruct_set_symmetries(astruct, disableSym, tol, elecfield, nspin)
          use module_base
          use module_atoms, only: atomic_structure,deallocate_symmetry_data
          use defs_basis
          use m_ab6_symmetry
          implicit none
          type(atomic_structure), intent(inout) :: astruct
          logical, intent(in) :: disableSym
          real(gp), intent(in) :: tol
          real(gp), intent(in) :: elecfield(3)
          integer, intent(in) :: nspin
        end  subroutine astruct_set_symmetries

        subroutine sic_input_variables_default(in)
          use module_base
          use module_types
          implicit none
          type(input_variables), intent(inout) :: in
        end subroutine sic_input_variables_default

        subroutine psp_from_data(symbol, nzatom, nelpsp, npspcode, ixc, psppar, exists)
          use module_base
          use module_xc
          implicit none
          character(len = *), intent(in) :: symbol
          integer, intent(inout) :: ixc
          integer, intent(out) :: nzatom, nelpsp, npspcode
          real(gp), intent(out) :: psppar(0:4,0:6)
          logical, intent(out) :: exists
        end subroutine psp_from_data

        subroutine total_energies(energs, iter, iproc)
          use module_base
          use module_types
          implicit none
          type(energy_terms), intent(inout) :: energs
          integer, intent(in) :: iter, iproc
        end subroutine total_energies

        subroutine ext_buffers(periodic,nl,nr)
          implicit none
          logical, intent(in) :: periodic
          integer, intent(out) :: nl,nr
        end subroutine ext_buffers

        subroutine ind_positions(periodic,i,n,j,go)
          implicit none
          logical, intent(in) :: periodic
          integer, intent(in) :: i,n
          logical, intent(out) :: go
          integer, intent(out) :: j
        end subroutine ind_positions

        subroutine eigensystem_info(iproc,nproc,tolerance,nvctr,orbs,psi)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc,nvctr
          real(gp), intent(in) :: tolerance
          type(orbitals_data), intent(inout) :: orbs
          real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
        end subroutine eigensystem_info

        subroutine determine_sparsity_pattern(iproc, nproc, orbs, lzd, nnonzero, nonzero)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          integer,intent(out) :: nnonzero
          integer,dimension(:,:),pointer,intent(out) :: nonzero
        end subroutine determine_sparsity_pattern

        subroutine determine_sparsity_pattern_distance(orbs, lzd, astruct, cutoff, nnonzero, nonzero)
          use module_base
          use module_types
          implicit none
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          type(atomic_structure),intent(in) :: astruct
          real(kind=8),dimension(lzd%nlr),intent(in) :: cutoff
          integer,intent(out) :: nnonzero
          integer,dimension(:,:),pointer,intent(out) :: nonzero
        end subroutine determine_sparsity_pattern_distance

        subroutine init_sparse_matrix_wrapper(iproc, nproc, nspin, orbs, lzd, astruct, store_index, imode, smat)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc, nspin, imode
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          type(atomic_structure),intent(in) :: astruct
          logical,intent(in) :: store_index
          type(sparse_matrix), intent(out) :: smat
        end subroutine init_sparse_matrix_wrapper

        subroutine check_accur_overlap_minus_one_sparse(iproc, nproc, smat, norb, norbp, isorb, nseq, nout, &
                   ivectorindex, onedimindices, amat_seq, bmatp, power, &
                   max_error, mean_error, dmat_seq, cmatp)
          use module_base
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseq, nout, power
          type(sparse_matrix) :: smat
          integer,dimension(nseq),intent(in) :: ivectorindex
          integer,dimension(4,nout) :: onedimindices
          real(kind=8),dimension(nseq),intent(in) :: amat_seq
          real(kind=8),dimension(norb,norbp),intent(in) :: bmatp
          real(kind=8),intent(out) :: max_error, mean_error
          real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
          real(kind=8),dimension(norb,norbp),intent(in),optional :: cmatp
        end subroutine check_accur_overlap_minus_one_sparse

        subroutine calculate_kernel_and_energy(iproc,nproc,denskern,ham,denskern_mat,ham_mat,&
                   energy,coeff,orbs,tmb_orbs,calculate_kernel)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer, intent(in) :: iproc, nproc
          type(sparse_matrix), intent(in) :: ham
          type(sparse_matrix), intent(inout) :: denskern
          type(matrices),intent(in) :: ham_mat
          type(matrices),intent(out) :: denskern_mat
          logical, intent(in) :: calculate_kernel
          real(kind=gp), intent(out) :: energy
          type(orbitals_data), intent(in) :: orbs, tmb_orbs
          real(kind=gp), dimension(denskern%nfvctr,tmb_orbs%norb), intent(in) :: coeff
        end subroutine calculate_kernel_and_energy

        subroutine calc_site_energies_transfer_integrals(iproc,nproc,meth_overlap,input_frag,&
                   ref_frags,orbs,ham,ham_mat,ovrlp,ovrlp_mat,KS_overlap)
          use module_base
          use module_types
          use yaml_output
          use module_fragments
          implicit none
          integer, intent(in) :: iproc, nproc, meth_overlap
          type(fragmentInputParameters), intent(in) :: input_frag
          type(orbitals_data), intent(in) :: orbs
          type(sparse_matrix), intent(inout) :: ham, ovrlp
          type(sparse_matrix),dimension(ham%nspin),intent(inout) :: KS_overlap
          type(matrices), intent(inout) :: ovrlp_mat, ham_mat
          type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
        end subroutine calc_site_energies_transfer_integrals

        subroutine calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ham_mat,ovrlp,ovrlp_mat,&
                   homo_coeffs1,homo_coeffs2,homo_ham,homo_ovrlp)
          use module_defs, only:gp
          use module_types
          use module_fragments
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer, intent(in) :: iproc, nproc, nstates
          type(orbitals_data), intent(in) :: orbs
          type(sparse_matrix), intent(inout) :: ham, ovrlp
          type(matrices),intent(inout) :: ovrlp_mat, ham_mat
          real(kind=gp), dimension(ovrlp%nfvctr,nstates), intent(in) :: homo_coeffs1, homo_coeffs2
          real(kind=gp), dimension(nstates), intent(inout) :: homo_ham, homo_ovrlp
        end subroutine calc_transfer_integral

        subroutine overlap_minus_one_half_serial(iproc, nproc, iorder, power, blocksize, &
                   norb, ovrlp_matrix, inv_ovrlp_matrix, check_accur, &
                   smat, max_error, mean_error)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc, iorder, blocksize, power, norb
          real(kind=8),dimension(norb,norb),intent(in) :: ovrlp_matrix
          real(kind=8),dimension(:,:),pointer,intent(out) :: inv_ovrlp_matrix
          type(sparse_matrix),intent(in) :: smat
          logical,intent(in) :: check_accur
          real(kind=8),intent(out),optional :: max_error, mean_error
        end subroutine overlap_minus_one_half_serial

        subroutine calculate_weight_matrix_lowdin(weight_matrix,weight_matrix_,nfrag_charged,ifrag_charged,tmb,input_frag,&
             ref_frags,calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap)
          use module_defs, only: gp
          use module_types
          use module_fragments
          implicit none
          type(sparse_matrix), intent(inout) :: weight_matrix
           type(matrices), intent(inout) :: weight_matrix_
          type(fragmentInputParameters),intent(in) :: input_frag
          type(dft_wavefunction), intent(inout) :: tmb
          logical, intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
          type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
          integer, intent(in) :: nfrag_charged, meth_overlap
          integer, dimension(2), intent(in) :: ifrag_charged
          !local variables
          integer :: ifrag,iorb,ifrag_ref,isforb,istat,ierr
          real(kind=gp), allocatable, dimension(:,:) :: proj_mat, proj_ovrlp_half, weight_matrixp
        end subroutine calculate_weight_matrix_lowdin

        subroutine calculate_weight_matrix_using_density(iproc,cdft,tmb,at,input,GPU,denspot)
          use module_base
          use module_types
          use constrained_dft, only: cdft_data
          use module_fragments
          implicit none
          integer,intent(in) :: iproc
          type(cdft_data), intent(inout) :: cdft
          type(atoms_data), intent(in) :: at
          type(input_variables),intent(in) :: input
          type(dft_wavefunction), intent(inout) :: tmb
          type(DFT_local_fields), intent(inout) :: denspot
          type(GPU_pointers),intent(inout) :: GPU
        end subroutine calculate_weight_matrix_using_density

        subroutine fragment_coeffs_to_kernel(iproc,input,input_frag_charge,ref_frags,tmb,ksorbs,overlap_calculated,&
          nstates_max,cdft)
          use yaml_output
          use module_base
          use module_types
          use module_fragments
          implicit none
          type(DFT_wavefunction), intent(inout) :: tmb
          type(input_variables), intent(in) :: input
          type(system_fragment), dimension(input%frag%nfrag_ref), intent(inout) :: ref_frags
          type(orbitals_data), intent(inout) :: ksorbs
          logical, intent(inout) :: overlap_calculated
          real(kind=gp), dimension(input%frag%nfrag), intent(in) :: input_frag_charge
          integer, intent(in) :: iproc
          integer, intent(out) :: nstates_max ! number of states in total if we consider all partially occupied fragment states to be fully occupied
          logical, intent(in) :: cdft
        end subroutine fragment_coeffs_to_kernel

        subroutine check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,&
                   smat,max_error,mean_error)
          use module_base
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in) :: iproc, nproc, norb, norbp, isorb, power
          real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
          type(sparse_matrix),intent(in) :: smat
          real(kind=8),intent(out) :: max_error, mean_error
        end subroutine check_accur_overlap_minus_one

        subroutine max_matrix_diff(iproc, norb, mat1, mat2, smat, max_deviation, mean_deviation)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in):: iproc, norb
          real(8),dimension(norb,norb),intent(in):: mat1, mat2
          type(sparse_matrix),intent(in) :: smat
          real(8),intent(out):: max_deviation, mean_deviation
        end subroutine max_matrix_diff

        subroutine max_matrix_diff_parallel(iproc, norb, norbp, isorb, mat1, mat2, &
                   smat, max_deviation, mean_deviation)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in):: iproc, norb, norbp, isorb
          real(8),dimension(norb,norbp),intent(in):: mat1, mat2
          type(sparse_matrix),intent(in) :: smat
          real(8),intent(out):: max_deviation, mean_deviation
        end subroutine max_matrix_diff_parallel

        subroutine deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, ovrlp, &
                   smat, max_deviation, mean_deviation)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer,intent(in):: iproc, nproc, norb, norbp, isorb
          real(8),dimension(norb,norbp),intent(in):: ovrlp
          type(sparse_matrix),intent(in) :: smat
          real(8),intent(out):: max_deviation, mean_deviation
        end subroutine deviation_from_unity_parallel

        subroutine estimate_energy_change(npsidim_orbs, orbs, lzd, nspin, psidiff, hpsi_noprecond, delta_energy)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: npsidim_orbs, nspin
          type(orbitals_data),intent(in) :: orbs
          type(local_zone_descriptors),intent(in) :: lzd
          real(kind=8),dimension(npsidim_orbs),intent(in) :: psidiff, hpsi_noprecond
          real(kind=8),intent(out) :: delta_energy
        end subroutine estimate_energy_change

        subroutine chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                   norb, norbp, isorb, fermi, chebyshev_polynomials, ncalc, cc, kernelp)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
        
          ! Calling arguments
          integer,intent(in) :: iproc, nproc, nsize_polynomial, npl, norb, norbp, isorb, ncalc
          type(sparse_matrix),intent(in) :: fermi
          real(kind=8),dimension(nsize_polynomial,npl),intent(in) :: chebyshev_polynomials
          real(kind=8),dimension(npl,ncalc),intent(in) :: cc
          real(kind=8),dimension(norb,norbp,ncalc),intent(out) :: kernelp
        end subroutine chebyshev_fast

        subroutine init_sparse_matrix_for_KSorbs(iproc, nproc, orbs, input, nextra, smat, smat_extra)
          use module_base
          use module_types
          use sparsematrix_base, only: sparse_matrix
          implicit none
          integer, intent(in) :: iproc, nproc, nextra
          type(orbitals_data), intent(in) :: orbs
          type(input_variables), intent(in) :: input
          type(sparse_matrix),dimension(:),pointer,intent(out) :: smat, smat_extra
        end subroutine init_sparse_matrix_for_KSorbs

        subroutine ice(iproc, nproc, norder_polynomial, ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc, norder_polynomial, ncalc
          type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
          integer,dimension(ncalc) :: ex
          type(matrices),intent(in) :: ovrlp_mat
          type(matrices),dimension(ncalc),intent(out) :: inv_ovrlp
        end subroutine ice
        
        subroutine scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
                   smat1, mat1, i1shift, smat2, mat2, i2shift, &
                   matscal_compr, scale_factor, shift_value)
          use module_base
          use sparsematrix_base, only: sparse_matrix, matrices
          use foe_base, only: foe_data
          integer,intent(in) :: iproc, nproc, ispin, i1shift
          type(foe_data),intent(in) :: foe_obj
          type(sparse_matrix),intent(in) :: smatl, smat1
          type(matrices),intent(in) :: mat1
          type(sparse_matrix),intent(in),optional :: smat2
          type(matrices),intent(in),optional :: mat2
          integer,intent(in),optional :: i2shift
          real(kind=8),dimension(smatl%nvctr),intent(out) :: matscal_compr
          real(kind=8),intent(out) :: scale_factor, shift_value
        end subroutine scale_and_shift_matrix

        subroutine build_gradient(iproc, nproc, tmb, target_function, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f)
          use module_base
          use module_types
          implicit none
          integer,intent(in) :: iproc, nproc, target_function
          type(DFT_wavefunction),intent(inout) :: tmb
          real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(inout) :: hpsit_c
          real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(inout) :: hpsit_f
          real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(out) :: hpsittmp_c !<workarray
          real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(out) :: hpsittmp_f !<workarray
        end subroutine build_gradient

  end interface
END MODULE module_interfaces

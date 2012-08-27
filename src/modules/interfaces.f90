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
!!  Interfaces of:
!!  - call_bigdft
!!  - geopt
!!  - geopt_input_variables
!!  - conjgrad
!!  - copy_old_wavefunctions
!!  - system_size
!!  - MemoryEstimator
!!  - createWavefunctionsDescriptors
!!  - createProjectorsArrays
!!  - createDensPotDescriptors
!!  - createIonicPotential
!!  - import_gaussians
!!  - input_wf_diag
!!  - reformatmywaves
!!  - first_orthon
!!  - sumrho
!!  - LocalHamiltonianApplication
!!  - hpsitopsi
!!  - last_orthon
!!  - local_forces
!!  - orbitals_descriptors
!!  - projectors_derivatives
!!  - nonlocal_forces
!!  - CalculateTailCorrection
!!  - reformatonewave
module module_interfaces

   implicit none

   interface

      subroutine call_bigdft(nproc,iproc,atoms,rxyz,in,energy,fxyz,strten,fnoise,rst,infocode)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(input_variables),intent(inout) :: in
         type(atoms_data), intent(inout) :: atoms
         type(restart_objects), intent(inout) :: rst
         integer, intent(inout) :: infocode
         real(gp), intent(out) :: energy,fnoise
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
         real(gp), dimension(6), intent(out) :: strten
         real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
      END SUBROUTINE call_bigdft

      subroutine geopt(nproc,iproc,pos,at,fxyz,strten,epot,rst,in,ncount_bigdft)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: nproc,iproc
        type(atoms_data), intent(inout) :: at
        type(input_variables), intent(inout) :: in
        type(restart_objects), intent(inout) :: rst
        real(gp), intent(inout) :: epot
        integer, intent(inout) :: ncount_bigdft
        real(gp), dimension(3*at%nat), intent(inout) :: pos
        real(gp), dimension(6), intent(inout) :: strten
        real(gp), dimension(3*at%nat), intent(inout) :: fxyz
      END SUBROUTINE geopt

      subroutine kswfn_optimization_loop(iproc, nproc, o, &
           & alphamix, idsx, inputpsi, KSwfn, denspot, nlpspd, proj, energs, atoms, rxyz, GPU, xcstr, &
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
        type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
        real(kind=8), dimension(:), pointer :: proj
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
        type(input_variables), intent(in) :: in !<todo: Remove me
      END SUBROUTINE kswfn_optimization_loop

     subroutine timing(iproc,category,action)
       implicit none
       integer, intent(in) :: iproc
       character(len=*), intent(in) :: category
       character(len=2), intent(in) :: action
     end subroutine timing

      subroutine copy_old_wavefunctions(nproc,orbs,n1,n2,n3,wfd,psi,&
            &   n1_old,n2_old,n3_old,wfd_old,psi_old)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: nproc,n1,n2,n3
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(inout) :: wfd,wfd_old
         integer, intent(out) :: n1_old,n2_old,n3_old
         real(wp), dimension(:), pointer :: psi,psi_old
      END SUBROUTINE copy_old_wavefunctions

      subroutine system_properties(iproc,nproc,in,at,orbs,radii_cf,nelec)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         integer, intent(out) :: nelec
         type(input_variables), intent(in) :: in
         type(atoms_data), intent(inout) :: at
         type(orbitals_data), intent(inout) :: orbs
         real(gp), dimension(at%ntypes,3), intent(out) :: radii_cf
      END SUBROUTINE system_properties

      subroutine system_size(iproc,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,Glr,shift)
         !n(c) use module_base
         use module_types
         implicit none
         type(atoms_data), intent(inout) :: atoms
         integer, intent(in) :: iproc
         real(gp), intent(in) :: crmult,frmult
         real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
         real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
         real(gp), intent(inout) :: hx,hy,hz
         type(locreg_descriptors), intent(out) :: Glr
         real(gp), dimension(3), intent(out) :: shift
      END SUBROUTINE system_size

      subroutine standard_inputfile_names(inputs, radical, nproc)
         use module_types
         implicit none
         type(input_variables), intent(out) :: inputs
         character(len = *), intent(in) :: radical
         integer, intent(in) :: nproc
      END SUBROUTINE standard_inputfile_names

      subroutine read_input_variables(iproc,posinp,inputs,atoms,rxyz)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: posinp
         integer, intent(in) :: iproc
         type(input_variables), intent(inout) :: inputs
         type(atoms_data), intent(out) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
      END SUBROUTINE read_input_variables

      subroutine read_input_parameters(iproc,inputs,atoms,rxyz)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc
         type(input_variables), intent(inout) :: inputs
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
      END SUBROUTINE read_input_parameters

      subroutine read_atomic_file(file,iproc,at,rxyz,status,comment,energy,fxyz)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: file
         integer, intent(in) :: iproc
         type(atoms_data), intent(inout) :: at
         real(gp), dimension(:,:), pointer :: rxyz
         integer, intent(out), optional :: status
         real(gp), intent(out), optional :: energy
         real(gp), dimension(:,:), pointer, optional :: fxyz
         character(len = 1024), intent(out), optional :: comment
      END SUBROUTINE read_atomic_file

      !> @author
      !! Written by Laurent K Beland 2011 UdeM
      !! For QM/MM implementation of BigDFT-ART
      subroutine initialize_atomic_file(iproc,at,rxyz)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc
         type(atoms_data), intent(inout) :: at
         real(gp), dimension(:,:), pointer :: rxyz
      END SUBROUTINE initialize_atomic_file

      subroutine read_xyz_positions(iproc,ifile,atoms,rxyz,comment_,energy_,fxyz_,getLine)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,ifile
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
         real(gp), intent(out) :: energy_
         real(gp), dimension(:,:), pointer :: fxyz_
         character(len = 1024), intent(out) :: comment_
         interface
            subroutine getline(line,ifile,eof)
               integer, intent(in) :: ifile
               character(len=150), intent(out) :: line
               logical, intent(out) :: eof
            END SUBROUTINE getline
         end interface
      END SUBROUTINE read_xyz_positions

      subroutine read_ascii_positions(iproc,ifile,atoms,rxyz,comment_,energy_,fxyz_,getline)
         ! use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,ifile
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
         real(gp), intent(out) :: energy_
         real(gp), dimension(:,:), pointer :: fxyz_
         character(len = 1024), intent(out) :: comment_
         interface
            subroutine getline(line,ifile,eof)
               integer, intent(in) :: ifile
               character(len=150), intent(out) :: line
               logical, intent(out) :: eof
            END SUBROUTINE getline
         end interface
      END SUBROUTINE read_ascii_positions

      subroutine read_yaml_positions(filename, atoms,rxyz,comment_,energy_,fxyz_)
        use module_base
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        type(atoms_data), intent(inout) :: atoms
        real(gp), dimension(:,:), pointer :: rxyz
        real(gp), intent(out) :: energy_
        real(gp), dimension(:,:), pointer :: fxyz_
        character(len = 1024), intent(out) :: comment_
      END SUBROUTINE read_yaml_positions

      subroutine write_atomic_file(filename,energy,rxyz,atoms,comment,forces)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename,comment
         type(atoms_data), intent(in) :: atoms
         real(gp), intent(in) :: energy
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
         real(gp), dimension(3,atoms%nat), intent(in), optional :: forces
      END SUBROUTINE write_atomic_file

      subroutine MemoryEstimator(nproc,idsx,lr,nat,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,peakmem)
         !n(c) use module_base
         use module_types
         implicit none
         !Arguments
         integer, intent(in) :: nproc,idsx,nat,norb,nspin,nprojel
         integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
         type(locreg_descriptors), intent(in) :: lr
         real(kind=8), intent(out) :: peakmem
      END SUBROUTINE MemoryEstimator

      subroutine check_closed_shell(orbs,lcs)
         !n(c) use module_base
         use module_types
         implicit none
         type(orbitals_data), intent(in) :: orbs
         logical, intent(out) :: lcs
      END SUBROUTINE check_closed_shell

      subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs,simple,basedist)
         !n(c) use module_base
         use module_types
         implicit none
         logical, intent(in) :: simple !< simple calculation of the repartition
         integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
         integer, intent(in) :: nspinor
         type(orbitals_data), intent(out) :: orbs
         real(gp), dimension(nkpt), intent(in) :: wkpt
         real(gp), dimension(3,nkpt), intent(in) :: kpt
         integer, dimension(0:nproc-1), intent(in), optional :: basedist 
      END SUBROUTINE orbitals_descriptors

      subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms,basedist)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(inout) :: orbs
         type(communications_arrays), intent(out) :: comms
         integer, dimension(0:nproc-1,orbs%nkpts), intent(in), optional :: basedist
      END SUBROUTINE orbitals_communicators


     subroutine orbitals_descriptors_forLinear(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
       integer, intent(in) :: nspinor
       type(orbitals_data), intent(out) :: orbs
       real(gp), dimension(nkpt), intent(in) :: wkpt
       real(gp), dimension(3,nkpt), intent(in) :: kpt
     END SUBROUTINE orbitals_descriptors_forLinear

      subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
            &   crmult,frmult,Glr,output_denspot)
         !n(c) use module_base
         use module_types
         implicit none
         !Arguments
         type(atoms_data), intent(in) :: atoms
         integer, intent(in) :: iproc
         real(gp), intent(in) :: hx,hy,hz,crmult,frmult
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
         real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
         type(locreg_descriptors), intent(inout) :: Glr
         logical, intent(in), optional :: output_denspot
      END SUBROUTINE createWavefunctionsDescriptors

     subroutine createProjectorsArrays(iproc,lr,rxyz,at,orbs,&
            &   radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
         !n(c) use module_base
         use module_types
         implicit none
       integer, intent(in) :: iproc
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(kind=8), intent(in) :: cpmult,fpmult,hx,hy,hz
       type(locreg_descriptors),intent(in) :: lr
         real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
         real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
         type(nonlocal_psp_descriptors), intent(out) :: nlpspd
         real(kind=8), dimension(:), pointer :: proj
      END SUBROUTINE createProjectorsArrays

      subroutine density_descriptors(iproc,nproc,nspin,crmult,frmult,atoms,dpbox,&
           rho_commun,rxyz,radii_cf,rhodsc)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        real(gp), intent(in) :: crmult,frmult
        type(atoms_data), intent(in) :: atoms
        type(denspot_distribution), intent(in) :: dpbox
        character(len=3), intent(in) :: rho_commun
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
        real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
        type(rho_descriptors), intent(out) :: rhodsc
      end subroutine density_descriptors
!!$      subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
!!$            &   rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
!!$         n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
!!$         !n(c) use module_base
!!$         use module_types
!!$         implicit none
!!$         !Arguments
!!$         character(len=1), intent(in) :: datacode
!!$         character(len=3), intent(in) :: rho_commun
!!$         integer, intent(in) :: iproc,nproc,ixc,nspin
!!$         real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
!!$         type(atoms_data), intent(in) :: atoms
!!$         type(grid_dimensions), intent(in) :: gdim
!!$         real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
!!$         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!$         integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
!!$         type(rho_descriptors), intent(out) :: rhodsc
!!$         integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
!!$         integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
!!$      END SUBROUTINE createDensPotDescriptors

      subroutine createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs, &
            &   radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
         pcproj_data , Glr)

         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,n1,n2,n3
         real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs

         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
         real(gp), intent(in):: ecut_pc

         type(pcproj_data_type) ::pcproj_data

         type(locreg_descriptors),  intent(in):: Glr

      END SUBROUTINE createPcProjectorsArrays


      subroutine applyPCprojectors(orbs,at,&
            &   hx,hy,hz,Glr,PPD,psi,hpsi, dotest)

         use module_base
         use module_types

         type(orbitals_data), intent(inout) :: orbs
         type(atoms_data) :: at
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: Glr
         type(pcproj_data_type) ::PPD
         real(wp), dimension(:), pointer :: psi, hpsi
         logical, optional :: dotest
      END SUBROUTINE applyPCprojectors


      subroutine applyPAWprojectors(orbs,at,&
            &   hx,hy,hz,Glr,PAWD,psi,hpsi,  paw_matrix, dosuperposition , &
         sup_iatom, sup_l, sup_arraym) !, sup_arraychannel)

         use module_base
         use module_types

         type(orbitals_data), intent(inout) :: orbs
         type(atoms_data) :: at
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: Glr
         type(pawproj_data_type) ::PAWD
         real(wp), dimension(:), pointer :: psi, hpsi, paw_matrix
         logical dosuperposition
         integer, optional :: sup_iatom, sup_l
         real(wp) , dimension(:), pointer, optional :: sup_arraym !, sup_arraychannel

       END SUBROUTINE applyPAWprojectors

       subroutine IonicEnergyandForces(iproc,nproc,at,hxh,hyh,hzh,elecfield,&
            & rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,psoffset,n1,n2,n3,&
            & n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
         use module_base
         use module_types
         implicit none
         type(atoms_data), intent(in) :: at
         integer, intent(in) :: iproc,nproc,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,dispersion
         real(gp), intent(in) :: hxh,hyh,hzh
         real(gp), dimension(3), intent(in) :: elecfield
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         type(coulomb_operator), intent(in) :: pkernel
         real(gp), intent(out) :: eion,edisp,psoffset
         real(dp), dimension(6),intent(out) :: ewaldstr
         real(gp), dimension(:,:), pointer :: fion,fdisp
         real(dp), dimension(*), intent(out) :: pot_ion
       END SUBROUTINE IonicEnergyandForces

       subroutine createIonicPotential(geocode,iproc,nproc,verb,at,rxyz,&
            hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,psoffset)
         use module_base
         use module_types
         implicit none
         character(len=1), intent(in) :: geocode
         integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
         logical, intent(in) :: verb
         real(gp), intent(in) :: hxh,hyh,hzh,psoffset
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3), intent(in) :: elecfield
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         type(coulomb_operator), intent(in) :: pkernel
         real(wp), dimension(*), intent(inout) :: pot_ion
       END SUBROUTINE createIonicPotential

       subroutine input_wf_empty(iproc, nproc, psi, hpsi, psit, orbs, &
            & band_structure_filename, input_spin, atoms, d, denspot)
         use module_defs
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

       subroutine input_wf_random(iproc, nproc, psi, orbs)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_random

       subroutine input_wf_cp2k(iproc, nproc, nspin, atoms, rxyz, Lzd, &
            & hx, hy, hz, psi, orbs)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, nspin
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%nat), intent(in) :: rxyz
         type(local_zone_descriptors), intent(in) :: Lzd
         real(gp), intent(in) :: hx, hy, hz
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_cp2k

       subroutine input_wf_memory(iproc, atoms, &
            & rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old, &
            & rxyz, hx, hy, hz, d, wfd, psi, orbs)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%nat), intent(in) :: rxyz, rxyz_old
         real(gp), intent(in) :: hx, hy, hz, hx_old, hy_old, hz_old
         type(grid_dimensions), intent(in) :: d, d_old
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(wavefunctions_descriptors), intent(inout) :: wfd_old
         type(orbitals_data), intent(in) :: orbs
         real(wp), dimension(:), pointer :: psi, psi_old
       END SUBROUTINE input_wf_memory

       subroutine input_wf_disk(iproc, nproc, input_wf_format, d, hx, hy, hz, &
            & in, atoms, rxyz, rxyz_old, wfd, orbs, psi)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, input_wf_format
         type(grid_dimensions), intent(in) :: d
         real(gp), intent(in) :: hx, hy, hz
         type(input_variables), intent(in) :: in
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%nat), intent(in) :: rxyz
         real(gp), dimension(3, atoms%nat), intent(out) :: rxyz_old
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(orbitals_data), intent(inout) :: orbs
         real(wp), dimension(:), pointer :: psi
       END SUBROUTINE input_wf_disk

       subroutine input_wf_diag(iproc,nproc,at,denspot,&
            orbs,nvirt,comms,Lzd,energs,rxyz,&
            nlpspd,proj,ixc,psi,hpsi,psit,G,&
            nspin,symObj,GPU,input)
         ! Input wavefunctions are found by a diagonalization in a minimal basis set
         ! Each processors write its initial wavefunctions into the wavefunction file
         ! The files are then read by readwave
         ! @todo pass GPU to be a local variable of this routine (initialized and freed here)
         use module_base
         use module_types
         implicit none
         !Arguments
         integer, intent(in) :: iproc,nproc,ixc
         integer, intent(inout) :: nspin,nvirt
         type(atoms_data), intent(in) :: at
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd
         type(local_zone_descriptors), intent(in) :: Lzd
         type(communications_arrays), intent(in) :: comms
         type(orbitals_data), intent(inout) :: orbs
         type(energy_terms), intent(inout) :: energs
         type(DFT_local_fields), intent(inout) :: denspot
         type(GPU_pointers), intent(in) :: GPU
         type(input_variables):: input
         type(symmetry_data), intent(in) :: symObj
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
         type(gaussian_basis), intent(out) :: G !basis for davidson IG
         real(wp), dimension(:), pointer :: psi,hpsi,psit
       end subroutine input_wf_diag

       subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,&
            denspot,denspot0,nlpspd,proj,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,&
            wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,linear_start)
         use module_defs
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, inputpsi,  input_wf_format
         type(input_variables), intent(in) :: in
         type(GPU_pointers), intent(in) :: GPU
         real(gp), intent(in) :: hx_old,hy_old,hz_old
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3, atoms%nat), target, intent(in) :: rxyz
         type(DFT_local_fields), intent(inout) :: denspot
         type(DFT_wavefunction), intent(inout) :: KSwfn,tmb !<input wavefunction
         type(energy_terms), intent(inout) :: energs !<energies of the system
         real(gp), dimension(*), intent(out) :: denspot0 !< Initial density / potential, if needed
         real(wp), dimension(:), pointer :: psi_old
         integer, intent(out) :: norbv
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd
         real(kind=8), dimension(:), pointer :: proj
         type(grid_dimensions), intent(in) :: d_old
         real(gp), dimension(3, atoms%nat), intent(inout) :: rxyz_old
         type(wavefunctions_descriptors), intent(inout) :: wfd_old
         logical, intent(in) :: linear_start
       END SUBROUTINE input_wf

       subroutine reformatmywaves(iproc,orbs,at,&
            &   hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
         hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,n1_old,n2_old,n3_old,n1,n2,n3
         real(gp), intent(in) :: hx_old,hy_old,hz_old,hx,hy,hz
         type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(gp), dimension(3,at%nat), intent(in) :: rxyz,rxyz_old
         real(wp), dimension(wfd_old%nvctr_c+7*wfd_old%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: psi
      END SUBROUTINE reformatmywaves

      subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit,orthpar)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(communications_arrays), intent(in) :: comms
         type(orthon_data):: orthpar
         real(wp), dimension(:) , pointer :: psi,hpsi,psit
      END SUBROUTINE first_orthon

      subroutine density_and_hpot(dpbox,symObj,orbs,Lzd,pkernel,rhodsc,GPU,psi,rho,vh,hstrten)
        use module_base
        use module_types
        implicit none
        type(denspot_distribution), intent(in) :: dpbox
        type(rho_descriptors),intent(inout) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
        type(coulomb_operator), intent(in) :: pkernel
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        type(GPU_pointers), intent(inout) :: GPU
        real(gp), dimension(6), intent(out) :: hstrten
        real(dp), dimension(:), pointer :: rho,vh
      end subroutine density_and_hpot

      subroutine sumrho(dpbox,orbs,Lzd,GPU,symObj,rhodsc,psi,rho_p,mapping)
        use module_base
        use module_types
        implicit none
        !Arguments
        type(denspot_distribution), intent(in) :: dpbox
        type(rho_descriptors),intent(in) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
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

      subroutine rho_segkey(iproc,at,rxyz,crmult,frmult,radii_cf,&
            &   n1i,n2i,n3i,hxh,hyh,hzh,nspin,rho_d,iprint)
         !n(c) use module_base
         use module_types
         implicit none
         integer,intent(in) :: n1i,n2i,n3i,iproc,nspin
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
         real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
         logical,intent(in) :: iprint
         type(rho_descriptors),intent(inout) :: rho_d
       END SUBROUTINE rho_segkey

       subroutine LocalHamiltonianApplication(iproc,nproc,at,orbs,&
            Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
            energs,SIC,GPU,PotOrKin,pkernel,orbsocc,psirocc,dpbox,potential,comgp)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: PotOrKin !< if true, only the potential operator is applied
         integer, intent(in) :: iproc,nproc
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd 
         type(SIC_data), intent(in) :: SIC
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
       end subroutine LocalHamiltonianApplication

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
      END SUBROUTINE NonLocalHamiltonianApplication

      subroutine SynchronizeHamiltonianApplication(nproc,orbs,Lzd,GPU,hpsi,&
           ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: nproc
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd 
        type(GPU_pointers), intent(inout) :: GPU
        real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
      END SUBROUTINE SynchronizeHamiltonianApplication

      subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,idsx,iter
         type(DFT_wavefunction), intent(inout) :: wfn
      END SUBROUTINE hpsitopsi

      subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,wfd,comms,&
            &   psi,hpsi,psit,orthpar,passmat,& !mandatory
         orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,natsc,nspin
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(communications_arrays), target, intent(in) :: comms
         type(orbitals_data), target, intent(inout) :: orbs
         type(orthon_data), intent(in) :: orthpar
         real(wp), dimension(:), pointer :: psi,hpsi,psit
         real(wp), dimension(*), intent(out) :: passmat
         !optional arguments
         real(gp), optional, intent(in) :: etol
         type(orbitals_data), optional, intent(in) :: orbsv
         type(orbitals_data), optional, target, intent(in) :: orbse
         type(communications_arrays), optional, target, intent(in) :: commse
         integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
         real(wp), dimension(:), pointer, optional :: psivirt
      END SUBROUTINE DiagHam

      subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,iter
        real(wp), intent(out) :: evsum
        type(DFT_wavefunction), intent(inout) :: wfn
        logical, optional :: opt_keeppsit
      END SUBROUTINE last_orthon

      subroutine calculate_forces(iproc,nproc,psolver_groupsize,Glr,atoms,orbs,nlpspd,rxyz,hx,hy,hz,proj,i3s,n3p,nspin,&
           refill_proj,ngatherarr,rho,pot,potxc,psi,fion,fdisp,fxyz,&
           ewaldstr,hstrten,xcstr,strten,fnoise,pressure,psoffset)
        use module_base
        use module_types
        implicit none
        logical, intent(in) :: refill_proj
        integer, intent(in) :: iproc,nproc,i3s,n3p,nspin,psolver_groupsize
        real(gp), intent(in) :: hx,hy,hz,psoffset
        type(locreg_descriptors), intent(in) :: Glr
        type(atoms_data), intent(in) :: atoms
        type(locreg_descriptors) :: lr
        type(orbitals_data), intent(in) :: orbs
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd
        integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
        real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
        real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: rho,pot,potxc
        real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
        real(gp), dimension(6), intent(in) :: ewaldstr,hstrten,xcstr
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz,fion,fdisp
        real(gp), intent(out) :: fnoise,pressure
        real(gp), dimension(6), intent(out) :: strten
        real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
      END SUBROUTINE calculate_forces
      
      subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
            &   Glr,nlpspd,ncongt,pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,&
         proj,psi,output_denspot,ekin_sum,epot_sum,eproj_sum)
         !n(c) use module_base
         use module_types
         implicit none
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(locreg_descriptors), intent(in) :: Glr
         type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
         integer, intent(in) :: iproc,nproc,ncongt,nspin
         logical, intent(in) :: output_denspot
         real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
         real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
         real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
         real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
         real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
         real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
         real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
      END SUBROUTINE CalculateTailCorrection

      !added for abinit compatilbility
      subroutine reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,&
           n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: n1_old,n2_old,n3_old,n1,n2,n3
         real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz_old,rxyz
         real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
         real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
      END SUBROUTINE reformatonewave
      subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
            &   hx,hy,hz,at,wfd,rxyz_old,rxyz,psi,eval,psifscf)
         !n(c) use module_base
         use module_types
         implicit none
         logical, intent(in) :: useFormattedInput
         integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(atoms_data), intent(in) :: at
         real(gp), intent(in) :: hx,hy,hz
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), intent(out) :: eval
         real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
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
           rxyz,rhopot,nlpspd,proj,pkernel,psi,v,dpbox,GPU)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc
        integer, intent(in) :: nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd
        type(local_zone_descriptors), intent(inout) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        type(communications_arrays), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpbox
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
        type(coulomb_operator), intent(in) :: pkernel
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
      end subroutine davidson

      subroutine build_eigenvectors(iproc,norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinore,nspinor,&
            &   ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,passmat,nvirte,psivirt)
         use module_base
         implicit none
         !Arguments
         integer, intent(in) :: norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinor,ndim_hamovr,nspinore
         integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
         real(wp), dimension(nspin*ndim_hamovr), intent(in) :: hamovr
         real(wp), dimension(nvctrp,norbe), intent(in) :: psi
         real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: ppsit
         real(wp), dimension(*), intent(out) :: passmat
         integer, dimension(2), intent(in), optional :: nvirte
         real(wp), dimension(*), optional :: psivirt
         integer:: iproc
      END SUBROUTINE build_eigenvectors

      subroutine preconditionall(orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: ncong
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         real(dp), intent(out) :: gnrm,gnrm_zero
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp,orbs%nspinor), intent(inout) :: hpsi
      END SUBROUTINE preconditionall

      subroutine preconditionall2(iproc,nproc,orbs,Lzd,hx,hy,hz,ncong,hpsi,confdatarr,gnrm,gnrm_zero)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,ncong
        real(gp), intent(in) :: hx,hy,hz
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        real(dp), intent(out) :: gnrm,gnrm_zero
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
      end subroutine preconditionall2

      subroutine transpose_v(iproc,nproc,orbs,wfd,comms,psi,&
            &   work,outadd) !optional
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(communications_arrays), intent(in) :: comms
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: psi
         real(wp), dimension(:), pointer, optional :: work
         real(wp), intent(out), optional :: outadd
      END SUBROUTINE transpose_v

     subroutine transpose_v2(iproc,nproc,orbs,Lzd,comms,psi,&
          work,outadd) !optional
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors), intent(in) :: Lzd
       type(communications_arrays), intent(in) :: comms
       real(wp), dimension(:), pointer :: psi
       real(wp), dimension(:), pointer, optional :: work
       real(wp), dimension(*), intent(out), optional :: outadd
     end subroutine

      subroutine untranspose_v(iproc,nproc,orbs,wfd,comms,psi,&
            &   work,outadd) !optional
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(communications_arrays), intent(in) :: comms
         real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
         real(wp), dimension(:), pointer, optional :: work
         real(wp), intent(out), optional :: outadd
      END SUBROUTINE untranspose_v

     subroutine plot_wf(orbname,nexpo,at,factor,lr,hx,hy,hz,rxyz,psi)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*) :: orbname
         integer, intent(in) :: nexpo
       real(dp), intent(in) :: factor
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         type(locreg_descriptors), intent(in) :: lr
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
      END SUBROUTINE plot_wf

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
         !n(c) use module_base
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
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         type(orbitals_data), intent(in) :: orbs
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:), pointer :: coeffs
         logical, optional :: opt_fillrxyz
      END SUBROUTINE read_gaussian_information

      subroutine restart_from_gaussians(iproc,nproc,orbs,Lzd,hx,hy,hz,psi,G,coeffs)
         !n(c) use module_base
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
            &   orbs,orbse,norbsc_arr,locrad,G,psigau,eks)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         integer, intent(inout) :: nvirt
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(gp), intent(out) :: eks
         integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
         real(gp), dimension(at%nat), intent(out) :: locrad
         type(orbitals_data), intent(out) :: orbse
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:,:), pointer :: psigau
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
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       integer,dimension(norb),intent(in):: mapping
       integer,dimension(at%nat),intent(in):: norbsPerAt
       real(gp), intent(out) :: eks
       integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
       real(gp), dimension(at%nat), intent(out) :: locrad
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
       real(gp),dimension(at%ntypes),intent(in),optional:: quartic_prefactor
     END SUBROUTINE inputguess_gaussian_orbitals_forLinear

     subroutine inputguess_gaussian_orbitals_withOnWhichAtom(iproc,nproc,at,rxyz,Glr,nvirt,nspin,&
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks,onWhichAtom)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin
       integer, intent(inout) :: nvirt
       type(atoms_data), intent(inout) :: at
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: Glr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(gp), intent(out) :: eks
       integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
       real(gp), dimension(at%nat), intent(out) :: locrad
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
       integer,dimension(orbse%norb),intent(out):: onWhichAtom
     END SUBROUTINE inputguess_gaussian_orbitals_withOnWhichAtom

     subroutine AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,&
          &   nspin,eks,scorb,G,gaucoeff,iorbtolr,mapping,quartic_prefactor)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: norbe,iproc
       integer, intent(in) :: norbsc,nspin
       type(atoms_data), intent(in) :: at
       logical, dimension(4,2,at%natsc), intent(in) :: scorb
       real(gp), dimension(3,at%nat), intent(in), target :: rxyz
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(gp), intent(out) :: eks
       integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
       real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff !fake interface for passing address
       integer,dimension(orbse%norb), optional, intent(in):: mapping
       real(gp),dimension(at%ntypes),intent(in),optional:: quartic_prefactor
      END SUBROUTINE AtomicOrbitals

      subroutine atomic_occupation_numbers(filename,ityp,nspin,at,nmax,lmax,nelecmax,neleconf,nsccode,mxpl,mxchg)
         use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(in) :: ityp,mxpl,mxchg,nspin,nmax,lmax,nelecmax,nsccode
         type(atoms_data), intent(inout) :: at
         !integer, dimension(nmax,lmax), intent(in) :: neleconf
         real(gp), dimension(nmax,lmax), intent(in) :: neleconf
      END SUBROUTINE atomic_occupation_numbers

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

      subroutine correct_hartree_potential(at,iproc,nproc,n1i,n2i,n3i,n3p,n3pi,n3d,&
            &   i3s,i3xcsh,hxh,hyh,hzh,pkernel,ngatherarr,&
         rhoref,pkernel_ref,pot_ion,rhopot,ixc,nspin,ehart,eexcu,vexcu,PSquiet,correct_offset)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=3), intent(in) :: PSquiet
         logical, intent(in) :: correct_offset
         integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,n3p,n3pi,n3d,nspin,ixc,i3xcsh,i3s
         real(gp), intent(in) :: hxh,hyh,hzh
         type(atoms_data), intent(in) :: at
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
         real(dp), dimension(n1i,n2i,max(n3d,1),nspin), intent(inout) :: rhoref
         real(dp), dimension(n1i,n2i,max(n3pi,1)), intent(inout) :: pot_ion
         real(dp), dimension(n1i,n2i,max(n3d,1),nspin), intent(inout) :: rhopot
         real(gp), intent(out) :: ehart,eexcu,vexcu
         type(coulomb_operator), intent(in) :: pkernel_ref,pkernel
      END SUBROUTINE correct_hartree_potential

      subroutine xabs_lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
           radii_cf,nlpspd,proj,Lzd,dpbox,potential,&
           energs,nspin,GPU,in_iat_absorber,&
           in , PAWD , orbs )
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        real(gp), intent(in) :: hx,hy,hz
        type(atoms_data), intent(in), target :: at
        type(nonlocal_psp_descriptors), intent(in), target :: nlpspd
        type(local_zone_descriptors), intent(inout), target :: Lzd
        type(denspot_distribution), intent(in), target :: dpbox
        real(gp), dimension(3,at%nat), intent(in), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), intent(in), target :: proj
        real(wp), dimension(max(dpbox%ndimpot,1),nspin), target :: potential
        type(energy_terms), intent(inout) :: energs
        type(GPU_pointers), intent(inout) , target :: GPU
        integer, intent(in) :: in_iat_absorber
        type(input_variables),intent(in), target :: in
        type(pawproj_data_type), target ::PAWD
        type(orbitals_data), intent(inout), target :: orbs
      END SUBROUTINE xabs_lanczos

      subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
            &   zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
         aeval,ng,psi,res,chrg,&
            &   Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw, PAWpatch, &
         psipsigrid)
         use module_base, only: gp

         implicit real(gp) (a-h,o-z)
         logical :: noproj, readytoexit
         integer, parameter :: n_int=1000
         dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
            &   hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
         gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
            &   pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
         potgrd(n_int),&
            &   rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
         occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
            &   vh(0:ng,0:ng,4,0:ng,0:ng,4),&
         res(noccmax,lmax+1),xp(0:ng),& 
         psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol),&
            &   psigrid_naked_2(Ngrid,Nsol), projgrid(Ngrid,3), &
         rhogrid(Ngrid), potgrid(Ngrid), psigrid_not_fitted(Ngrid,Nsol),&
            &   psigrid_not_fitted_2(Ngrid,Nsol),&
         vxcgrid(Ngrid), &
            &   Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
         H(Nsol, Nsol), &
            &   H_2(Nsol, Nsol), &
         Hcorrected(Nsol, Nsol), &
            &   Hadd(Nsol, Nsol), Egrid_tmp(Nsol),Egrid_tmp_2(Nsol), Etofit(Nsol), &
         Soverlap(Nsol,Nsol), Tpsigrid(Nsol,Ngrid ),Tpsigrid_dum(Nsol, Ngrid),valuesatp(Nsol), &
            &   PAWpatch(Npaw, Npaw ), Spsitildes(Npaw, Npaw), genS(Nsol,Nsol), genH(Nsol,Nsol) , dumH(Nsol,Nsol)

         real(gp) , optional :: psipsigrid(Ngrid, Nsol)


         real(gp) :: rgrid(Ngrid), ene_m, ene_p, factadd, rcond, fixfact
         real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid), dumgrid3(Ngrid)
         logical dofit
         integer real_start, iocc, iwork(Nsol), INFO, volta, ngrid_box_2
         character(1) EQUED
         integer ipiv(Nsol), Npaw
      END SUBROUTINE gatom_modified

      subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,&
            &   psi, aeval, occup, psp_modifier, &
         Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw,  PAWpatch , psipsigrid )

         use module_base, only: gp, memocc,ndebug
         implicit none
         integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax, Nsol, labs, Ngrid,  Ngrid_box
         real(gp), dimension(0:4,0:6), intent(in) :: psppar
         !! real(gp), dimension(:,:), intent(in) :: psppar
         integer, intent(in) :: psp_modifier, Npaw

         real(gp), dimension(ng+1), intent(out) :: expo

         integer, parameter :: n_int=1000

         real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi, Egrid(Nsol),&
            &   rgrid(Ngrid), psigrid(Ngrid,Nsol  )
         real(gp),   intent(out), optional  :: psipsigrid(Ngrid,Nsol  )
         real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
         real(gp):: PAWpatch(Npaw,Npaw)

         !local variables
      END SUBROUTINE abs_generator_modified

      subroutine xabs_cg(iproc,nproc,at,hx,hy,hz,rxyz,&
           &   radii_cf,nlpspd,proj,Lzd,dpbox,potential,&
           &   energs,nspin,GPU,in_iat_absorber,&
           &   in , rhoXanes, PAWD , PPD, orbs )
        use module_base
        use module_types
        implicit none
        integer  :: iproc,nproc,nspin
        real(gp)  :: hx,hy,hz
        type(atoms_data), target :: at
        type(nonlocal_psp_descriptors), target :: nlpspd
        type(local_zone_descriptors), target :: Lzd
        type(pcproj_data_type), target ::PPD
        type(denspot_distribution), intent(in), target :: dpbox
        real(gp), dimension(3,at%nat), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), target :: proj
        real(wp), dimension(max(dpbox%ndimpot,1),nspin), target :: potential
        real(wp), dimension(max(dpbox%ndimpot,1),nspin), target :: rhoXanes
        type(energy_terms), intent(inout) :: energs
        type(GPU_pointers), intent(inout) , target :: GPU
        integer, intent(in) :: in_iat_absorber
        type(pawproj_data_type), target ::PAWD
        type(input_variables),intent(in), target :: in
        type(orbitals_data), intent(inout), target :: orbs
      end subroutine xabs_cg

      subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
           radii_cf,nlpspd,proj,Lzd,dpbox,potential,&
           energs,nspin,GPU,in_iat_absorber,in, PAWD , orbs  )
        use module_base
        use module_types
        implicit none
        integer  :: iproc,nproc,nspin
        real(gp)  :: hx,hy,hz
        type(atoms_data), target :: at
        type(nonlocal_psp_descriptors), target :: nlpspd
        type(local_zone_descriptors), target :: Lzd
        type(denspot_distribution), intent(in), target :: dpbox
        real(gp), dimension(3,at%nat), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), target :: proj
        real(wp), dimension(max(dpbox%ndimpot,1),nspin), target :: potential
        type(energy_terms), intent(inout) :: energs
        type(GPU_pointers), intent(inout) , target :: GPU
        integer, intent(in) :: in_iat_absorber 
        type(input_variables),intent(in), target :: in
        type(pawproj_data_type), target ::PAWD
        type(orbitals_data), intent(inout), target :: orbs
      end subroutine xabs_chebychev

      subroutine cg_spectra(iproc,nproc,at,hx,hy,hz,rxyz,&
           radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
           energs,nspin,GPU,in_iat_absorber,in , PAWD  )! aggiunger a interface
         !n(c) use module_base
         use module_types
         implicit none
         integer  :: iproc,nproc,ndimpot,nspin
         real(gp)  :: hx,hy,hz
         type(atoms_data), target :: at
         type(nonlocal_psp_descriptors), target :: nlpspd
         type(locreg_descriptors), target :: lr
         integer, dimension(0:nproc-1,2), target :: ngatherarr 
         real(gp), dimension(3,at%nat), target :: rxyz
         real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
         real(wp), dimension(nlpspd%nprojel), target :: proj
         real(wp), dimension(max(ndimpot,1),nspin), target :: potential
         type(energy_terms), intent(inout) :: energs
         type(GPU_pointers), intent(inout) , target :: GPU
         integer, intent(in) :: in_iat_absorber
         type(pawproj_data_type), target ::PAWD

         type(input_variables),intent(in) :: in

      END SUBROUTINE cg_spectra


      subroutine eleconf(nzatom,nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
         implicit none
         integer, intent(in) :: nzatom,nvalelec
         character(len=2), intent(out) :: symbol
         real(kind=8), intent(out) :: rcov,rprb,ehomo,amu
         integer, parameter :: nmax=6,lmax=3
         real(kind=8), intent(out) :: neleconf(nmax,0:lmax)
         integer, intent(out) :: nsccode,mxpl,mxchg
      END SUBROUTINE eleconf

      !     subroutine psimix(iproc,nproc,orbs,comms,ads,ids,mids,idsx,energy,energy_old,alpha,&
      !          hpsit,psidst,hpsidst_sp,psit)
      !       use module_base
      !       use module_types
      !       implicit none
      !       integer, intent(in) :: iproc,nproc,ids,mids,idsx
      !       real(gp), intent(in) :: energy,energy_old
      !       type(orbitals_data), intent(in) :: orbs
      !       type(communications_arrays), intent(in) :: comms
      !       real(gp), intent(inout) :: alpha
      !       real(wp), dimension(:), pointer :: psit,hpsit,psidst
      !       real(sp), dimension(:), pointer :: hpsidst_sp
      !       real(wp), dimension(:,:,:), pointer :: ads
      !     END SUBROUTINE psimix
      !
      subroutine plot_density(iproc,nproc,filename,at,rxyz,box,nspin,rho)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        type(atoms_data), intent(in) :: at
        type(denspot_distribution), intent(in) :: box
        character(len=*), intent(in) :: filename
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(dp), dimension(max(box%ndimpot,1),nspin), target, intent(in) :: rho
      END SUBROUTINE plot_density

      subroutine read_density(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
            &   nat,rxyz,iatypes, znucl)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode
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
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode
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
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         character(len=1), intent(in) :: geocode
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

      !!$     subroutine read_density_cube(filename, n1i,n2i,n3i, nspin, hxh,hyh,hzh, nat, rxyz,  rho)
      !!$       !n(c) use module_base
      !!$       use module_types
      !!$       implicit none
      !!$       character(len=*), intent(in) :: filename
      !!$       integer, intent(out) ::  n1i,n2i,n3i
      !!$       integer, intent(in) :: nspin
      !!$       real(gp), intent(out) :: hxh,hyh,hzh
      !!$       real(gp), pointer :: rxyz(:,:)
      !!$       real(dp), dimension(:), pointer :: rho
      !!$       integer, intent(out) ::  nat
      !!$     END SUBROUTINE read_density_cube

      subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc, gaenes, &
            &   iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg )
         use module_types
         implicit none
         logical, intent(in) :: enlargerprb
         integer, intent(in) :: iproc,nspin,ng
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
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
         real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
         type(gaussian_basis_c), intent(out) :: G

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
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
      END SUBROUTINE local_analysis

      subroutine plot_gatom_basis(filename,iat,ngx,G,Gocc,rhocoeff,rhoexpo)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(in) :: iat,ngx
         type(gaussian_basis), intent(in) :: G
         real(wp), dimension(:), pointer :: Gocc
         real(wp), dimension((ngx*(ngx+1))/2), intent(out) :: rhoexpo
         real(wp), dimension((ngx*(ngx+1))/2,4), intent(out) :: rhocoeff
      END SUBROUTINE plot_gatom_basis

      subroutine calculate_rhocore(iproc,at,d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,i3s,n3d,i3xcsh,n3p
        real(gp), intent(in) :: hxh,hyh,hzh
        type(atoms_data), intent(in) :: at
        type(grid_dimensions), intent(in) :: d
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(:,:,:,:), pointer :: rhocore
      END SUBROUTINE calculate_rhocore

      subroutine XC_potential(geocode,datacode,iproc,nproc,mpi_comm,n01,n02,n03,ixc,hx,hy,hz,&
           rho,exc,vxc,nspin,rhocore,potxc,xcstr,dvxcdrho)
        use module_base
        use module_xc
        implicit none
        character(len=1), intent(in) :: geocode
        character(len=1), intent(in) :: datacode
        integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspin,mpi_comm
        real(gp), intent(in) :: hx,hy,hz
        real(gp), intent(out) :: exc,vxc
        real(dp), dimension(*), intent(inout) :: rho
        real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
        real(wp), dimension(*), intent(out) :: potxc
        real(dp), dimension(6), intent(out) :: xcstr
        real(dp), dimension(:,:,:,:), target, intent(out), optional :: dvxcdrho
      END SUBROUTINE XC_potential

      subroutine xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
           nxcl,nxcr,ixc,hx,hy,hz,rhopot,pot_ion,sumpion,zf,zfionxc,exc,vxc,nproc,nspden)
        use module_base
        use module_xc
        use interfaces_56_xc
        implicit none
        character(len=1), intent(in) :: geocode
        logical, intent(in) :: sumpion
        integer, intent(in) :: m1,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,ixc,nproc,nspden
        integer, intent(in) :: nwbl,nwbr
        real(gp), intent(in) :: hx,hy,hz
        real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rhopot
        real(wp), dimension(*), intent(in) :: pot_ion
        real(dp), dimension(md1,md3,md2/nproc), intent(out) :: zf
        real(wp), dimension(md1,md3,md2/nproc,nspden), intent(out) :: zfionxc
        real(dp), intent(out) :: exc,vxc
      END SUBROUTINE xc_energy

      subroutine direct_minimization(iproc,nproc,in,at,nvirt,rxyz,&
           rhopot,nlpspd,proj,pkernel,dpbox,GPU,KSwfn,VTwfn)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd
        type(denspot_distribution), intent(in) :: dpbox
        type(DFT_wavefunction), intent(inout) :: KSwfn,VTwfn
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
        type(coulomb_operator), intent(in) :: pkernel
        real(dp), dimension(*), intent(in), target :: rhopot
        type(GPU_pointers), intent(inout) :: GPU
      end subroutine direct_minimization

      subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
            &   hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=1), intent(in) :: geocode
         integer, intent(in) :: iproc,nproc,n3pi,i3s
         real(gp), intent(in) :: hxh,hyh,hzh
         real(gp), dimension(3), intent(in) :: shift
         type(input_variables), intent(in) :: in
         type(grid_dimensions), intent(in) :: grid
         type(coulomb_operator), intent(in) :: pkernel
         real(wp), dimension(*), intent(inout) :: pot_ion
      END SUBROUTINE CounterIonPotential

      subroutine gaussian_rism_basis(nat,radii,rxyz,G)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: nat
         real(gp), dimension(nat), intent(in) :: radii
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G
      END SUBROUTINE gaussian_rism_basis

      subroutine gaussian_hermite_basis(nhermitemax,nat,radii,rxyz,G)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: nat,nhermitemax
         real(gp), dimension(nat), intent(in) :: radii
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G  
      END SUBROUTINE gaussian_hermite_basis

      subroutine write_eigenvalues_data(nproc,etol,orbs,mom_vec)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: nproc
        real(gp), intent(in) :: etol
        type(orbitals_data), intent(in) :: orbs
        real(gp), dimension(:,:,:), intent(in), pointer :: mom_vec
      end subroutine write_eigenvalues_data
      
      subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt,output_wf_format)
         !n(c) use module_base
         use module_types
         implicit none
         logical, intent(in) :: occorbs
         integer, intent(in) :: iproc,nspin,nvirt,nplot,output_wf_format
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
       END SUBROUTINE write_eigen_objects

       subroutine full_local_potential(iproc,nproc,orbs,Lzd,iflag,dpbox,potential,pot,comgp)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: iproc,nproc,iflag
         type(orbitals_data),intent(in) :: orbs
         type(local_zone_descriptors),intent(in) :: Lzd
         type(denspot_distribution), intent(in) :: dpbox
         real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), target :: potential
         real(wp), dimension(:), pointer :: pot
         !type(p2pCommsGatherPot),intent(inout), optional:: comgp
         type(p2pComms),intent(inout), optional:: comgp
       END SUBROUTINE full_local_potential

      subroutine free_full_potential(nproc,flag,pot,subname)
         use module_base
         implicit none
         character(len=*), intent(in) :: subname
         integer, intent(in) :: nproc,flag
         real(wp), dimension(:), pointer :: pot
      END SUBROUTINE free_full_potential

      subroutine select_active_space(iproc,nproc,orbs,comms,mask_array,Glr,orbs_as,comms_as,psi,psi_as)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(locreg_descriptors), intent(in) :: Glr
         type(communications_arrays), intent(in) :: comms
         logical, dimension(orbs%norb*orbs%nkpts), intent(in) :: mask_array
         real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(in) :: psi
         type(orbitals_data), intent(out) :: orbs_as
         type(communications_arrays), intent(out) :: comms_as
         real(wp), dimension(:), pointer :: psi_as
      END SUBROUTINE select_active_space

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

      subroutine constrained_davidson(iproc,nproc,in,at,& 
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           hx,hy,hz,rxyz,rhopot,psi,v,dpbox,GPU)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc
        integer, intent(in) :: nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        type(communications_arrays), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpbox
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
        !v, that is psivirt, is transposed on input and direct on output
      end subroutine constrained_davidson

      subroutine local_hamiltonian(iproc,nproc,orbs,Lzd,hx,hy,hz,&
           ipotmethod,confdatarr,pot,psi,hpsi,pkernel,ixc,alphaSIC,ekin_sum,epot_sum,eSIC_DC,&
           dpbox,potential,comgp)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc,ipotmethod,ixc
        real(gp), intent(in) :: hx,hy,hz,alphaSIC
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
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

      subroutine NK_SIC_potential(lr,orbs,ixc,fref,hxh,hyh,hzh,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: ixc
         real(gp), intent(in) :: hxh,hyh,hzh,fref
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         type(coulomb_operator), intent(in) :: pkernel
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
         !real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
         real(wp), intent(inout) :: poti
         real(gp), intent(out) :: eSIC_DC
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,2*orbs%nspin), intent(in), optional :: potandrho 
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin), intent(out), optional :: wxdsave 
      END SUBROUTINE NK_SIC_potential

      subroutine isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
        !use module_base
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
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         integer, dimension(orbs%norb), optional :: orblist
         real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
         character(len=*), intent(in) :: filename
      END SUBROUTINE readmywaves

      
      subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out,iiorb)
         use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         logical, intent(in) :: lbin
         integer, intent(in) :: iorb,ispinor,unitfile
         type(orbitals_data), intent(in) :: orbs
         integer, intent(out) :: iorb_out
         integer,optional :: iiorb   
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
         integer,optional :: iiorb
      END SUBROUTINE filename_of_iorb

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

      subroutine denspot_communications(iproc_world,nproc_world,iproc,nproc,mpi_comm,&
           ixc,nspin,geocode,SICapproach,dpbox)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc_world,iproc,nproc,mpi_comm,ixc,nspin,nproc_world
        character(len=*), intent(in) :: geocode,SICapproach
        type(denspot_distribution), intent(inout) :: dpbox
      end subroutine denspot_communications

      subroutine allocateRhoPot(iproc,Glr,nspin,atoms,rxyz,denspot)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nspin
        type(locreg_descriptors), intent(in) :: Glr
        type(atoms_data), intent(in) :: atoms
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
        type(DFT_local_fields), intent(inout) :: denspot
      END SUBROUTINE allocateRhoPot

      !subroutine SWcalczone(nat,posa,boxl,tmp_force, this_atom,numnei,nei)
      !
      !
      !  !use SWpotential
      !  use defs, only : boundary,maxnei,iproc,MPI_COMM_WORLD
      !  
      !  implicit none
      !  
      !  integer, intent(in)                               :: nat
      !  real(kind=8), intent(in), dimension(3*nat) :: posa
      !  real(kind=8), dimension(3), intent(inout)          :: boxl
      !  integer, intent(in) :: this_atom
      !  real(8), intent(out), dimension(3*nat), target:: tmp_force
      !
      !
      !  integer, dimension(nat),intent(in) :: numnei 
      !  integer, dimension(nat,maxnei),intent(in) :: nei 
      !END SUBROUTINE SWcalczone

!!$    subroutine readmywaves(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
!!$         wfd,psi,orblist)
!!$      use module_base
!!$      use module_types
!!$      implicit none
!!$      integer, intent(in) :: iproc,n1,n2,n3
!!$      real(gp), intent(in) :: hx,hy,hz
!!$      type(wavefunctions_descriptors), intent(in) :: wfd
!!$      type(orbitals_data), intent(inout) :: orbs
!!$      type(atoms_data), intent(in) :: at
!!$      real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!$      integer, dimension(orbs%norb), optional :: orblist
!!$      real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
!!$      real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
!!$      character(len=*), intent(in) :: filename
!!$     end subroutine readmywaves


        subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,fnrm,&
                   infoBasisFunctions,nlpspd,proj,ldiis,SIC,tmb,&
                   tmblarge2, lhphilarge2, energs_base, ham)
          use module_base
          use module_types
          implicit none
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
          type(energy_terms),intent(inout) :: energs_base
          real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ham
        end subroutine getLocalizedBasis


    !!!subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, &
    !!!      input, hx, hy, hz, rxyz, nscatterarr, tag, confdatarr, onwhichatom)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  ! Calling arguments
    !!!  integer,intent(in):: iproc, nproc
    !!!  real(gp),intent(in):: hx, hy, hz
    !!!  type(locreg_descriptors),intent(in):: Glr
    !!!  type(orbitals_data),intent(in):: orbs
    !!!  type(atoms_data),intent(inout):: at
    !!!  type(nonlocal_psp_descriptors),intent(in):: nlpspd
    !!!  type(linearParameters),intent(inout):: lin
    !!!  type(input_variables),intent(in):: input
    !!!  real(8),dimension(3,at%nat),intent(in):: rxyz
    !!!  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
    !!!  integer,intent(inout):: tag
    !!!  type(confpot_data), dimension(:),pointer,intent(out) :: confdatarr
    !!!  integer,dimension(:),pointer:: onwhichatom
    !!!end subroutine allocateAndInitializeLinear



    subroutine transpose_vLIN(iproc, lproc, uproc, orbs, comms, psi, newComm, &
         work,outadd) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc, lproc, uproc, newComm
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)):: psi
      real(wp), dimension(:), pointer, optional :: work
      real(wp), dimension(*), intent(out), optional :: outadd
    end subroutine transpose_vLIN


    subroutine untranspose_vLIN(iproc, lproc, uproc, orbs, comms, psi, newComm, &
         work,outadd) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,lproc, uproc, newComm
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)):: psi
      real(wp), dimension(:), pointer, optional :: work
      real(wp), dimension(*), intent(out), optional :: outadd
    end subroutine untranspose_vLIN


    subroutine inputOrbitals(iproc,nproc,at,&
         orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
         nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
         nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,ixc,symObj
      integer, intent(inout) :: nspin,nvirt
      real(gp), intent(in) :: hx,hy,hz
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(inout) :: orbs
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(locreg_descriptors), intent(in) :: Glr
      type(communications_arrays), intent(in) :: comms
      type(GPU_pointers), intent(inout) :: GPU
      type(input_variables):: input
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
      type(gaussian_basis), intent(out) :: G 
      real(wp), dimension(:), pointer :: hpsi,psit,rhocore
      real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)):: psi
      type(coulomb_operator), intent(in) :: pkernel,pkernelseq
      integer, intent(in) :: potshortcut
      integer, dimension(*), intent(in) :: irrzon
      real(dp), dimension(*), intent(in) :: phnons
    END SUBROUTINE inputOrbitals
    
    subroutine psimix(iproc,nproc,ndim_psi,orbs,comms,diis,hpsit,psit)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,ndim_psi
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      real(wp), dimension(ndim_psi), intent(inout) :: psit,hpsit
    end subroutine psimix
    
    subroutine estimatePerturbedOrbitals(iproc, nproc, at, orbs, lr, input, orbsLIN, commsLIN, rxyz, nspin, &
        nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, rxyzParabola, perturbation)
    use module_base
    use module_types
    implicit none
    
    ! Calling arguments
    integer:: iproc, nproc
    type(atoms_data), intent(in) :: at
    type(orbitals_data):: orbs
    type(locreg_descriptors), intent(in) :: lr
    type(input_variables):: input
    type(orbitals_data):: orbsLIN
    type(communications_arrays):: commsLIN
    real(8),dimension(3,at%nat):: rxyz, rxyzParabola
    integer:: nspin
    type(nonlocal_psp_descriptors), intent(in) :: nlpspd
    real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
    integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
    integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
    real(dp), dimension(*), intent(inout) :: rhopot
    type(GPU_pointers), intent(inout) :: GPU
    type(coulomb_operator), intent(in) :: pkernelseq
    real(8),dimension(max(orbsLIN%npsidim_orbs,orbsLIN%npsidim_comp)):: phi
    real(8),dimension(3,at%nat):: perturbation
    end subroutine estimatePerturbedOrbitals
    
    !!subroutine psimixVariable(iproc,nproc,orbs,comms,diis,diisArr, hpsit,psit, quiet)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer, intent(in) :: iproc,nproc
    !!  type(orbitals_data), intent(in) :: orbs
    !!  type(communications_arrays), intent(in) :: comms
    !!  type(diis_objects), intent(inout) :: diis
    !!  type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
    !!  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
    !!  logical, optional:: quiet ! to avoid that the DIIS weights are written
    !!end subroutine psimixVariable
    
    
    
    !!subroutine diisstpVariable(iproc,nproc,orbs,comms,diis,diisArr,psit,quiet)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!! Arguments
    !!  integer, intent(in) :: nproc,iproc
    !!  type(orbitals_data), intent(in) :: orbs
    !!  type(communications_arrays), intent(in) :: comms
    !!  type(diis_objects), intent(inout) :: diis
    !!  type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
    !!  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
    !!  logical, optional:: quiet ! to avoid that the DIIS weights are written
    !!end subroutine diisstpVariable
    
    
    subroutine apply_potentialConfinement(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, rxyzConfinement, &
         hxh, hyh, hzh, potentialPrefac, confPotOrder, &
         ibyyzz_r) !optional
      use module_base
      implicit none
      integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot, confPotOrder
      real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
      real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
           -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
      integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
      real(gp), intent(out) :: epot
      real(8),dimension(3):: rxyzConfinement
      real(8):: hxh, hyh, hzh, potentialPrefac
    end subroutine apply_potentialConfinement

    !!!subroutine getLinearPsi(iproc,nproc,lzd,orbs,&
    !!!     at,rxyz,denspot,&
    !!!     GPU,&
    !!!     infoBasisFunctions,infoCoeff,itSCC,ebs,nlpspd,proj,&
    !!!     ldiis,orthpar,&
    !!!     blocksize_pdgemm,&
    !!!     comrp,blocksize_pdsyev,nproc_pdsyev,&
    !!!     hx,hy,hz,SIC,locrad,tmb, tmbder, tmbmix)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none

    !!!  ! Calling arguments
    !!!  integer,intent(in):: iproc, nproc, itSCC
    !!!  integer,intent(in):: blocksize_pdgemm
    !!!  integer,intent(in):: blocksize_pdsyev, nproc_pdsyev
    !!!  type(local_zone_descriptors),intent(inout):: lzd
    !!!  type(orbitals_data),intent(in) :: orbs
    !!!  !type(p2pCommsSumrho),intent(inout):: comsr
    !!!  !!type(p2pComms),intent(inout):: comsr
    !!!  !!type(matrixDescriptors),intent(in):: mad, lbmad
    !!!  !!type(overlapParameters),intent(inout):: op, lbop
    !!!  !!type(p2pComms),intent(inout):: comon, lbcomon
    !!!  !type(p2pCommsGatherPot):: comgp, lbcomgp
    !!!  !!!!type(p2pComms):: comgp, lbcomgp
    !!!  type(atoms_data),intent(in):: at
    !!!  real(8),dimension(3,at%nat),intent(in):: rxyz
    !!!  type(DFT_local_fields), intent(inout) :: denspot
    !!!  type(GPU_pointers),intent(inout):: GPU
    !!!  integer,intent(out):: infoBasisFunctions, infoCoeff
    !!!  real(8),intent(out):: ebs
    !!!  real(8),intent(in):: hx, hy, hz
    !!!  !real(8),dimension(llborbs%norb,orbs%norb),intent(in out):: coeff
    !!!  !real(8),dimension(max(llborbs%npsidim_orbs,llborbs%npsidim_comp)),intent(inout):: lphi
    !!!  !real(8),dimension(:),pointer,intent(inout):: lphi
    !!!  type(nonlocal_psp_descriptors),intent(in):: nlpspd
    !!!  real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
    !!!  !real(8),dimension(lorbs%norb,orbs%norb),intent(inout):: coeff_proj
    !!!  type(localizedDIISParameters),intent(inout):: ldiis
    !!!  type(orthon_data),intent(in):: orthpar
    !!!  !!type(confpot_data),dimension(lorbs%norbp),intent(in) :: confdatarr
    !!!  !real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout)::lphiRestart
    !!!  !real(8),dimension(:),pointer,intent(inout)::lphiRestart
    !!!  !type(p2pCommsRepartition),intent(inout):: comrp
    !!!  type(p2pComms),intent(inout):: comrp
    !!!  type(SIC_data),intent(in):: SIC
    !!!  real(8),dimension(lzd%nlr),intent(in):: locrad
    !!!  !!type(wfn_metadata),intent(inout):: wfnmd
    !!!  type(DFT_wavefunction),intent(inout):: tmb, tmbder, tmbmix
    !!!end subroutine getLinearPsi

    subroutine get_coeff(iproc,nproc,scf_mode,lzd,orbs,at,rxyz,denspot,&
               GPU, infoCoeff,ebs,nlpspd,proj,&
               SIC,tmb,fnrm,overlapmatrix,calculate_overlap_matrix,&
               tmblarge, lhphilarge, &
               ham, ldiis_coeff)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, scf_mode
      type(local_zone_descriptors),intent(inout):: lzd
      type(orbitals_data),intent(inout) :: orbs
      type(atoms_data),intent(in):: at
      real(8),dimension(3,at%nat),intent(in):: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers),intent(inout):: GPU
      integer,intent(out):: infoCoeff
      real(8),intent(out):: ebs, fnrm
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      type(SIC_data),intent(in):: SIC
      type(DFT_wavefunction),intent(inout):: tmb
      real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: overlapmatrix
      logical,intent(in):: calculate_overlap_matrix
      type(DFT_wavefunction),intent(inout):: tmblarge
      real(8),dimension(:),pointer,intent(inout):: lhphilarge
      real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in),optional:: ham
      type(localizedDIISParameters),intent(inout),optional:: ldiis_coeff
    end subroutine get_coeff


    !!!subroutine local_hamiltonianConfinement(iproc,orbs,lin,lr,hx,hy,hz,&
    !!!     nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, onWhichAtom, at, centralAtom)
    !!!  use module_base
    !!!  use module_types
    !!!  use libxc_functionals
    !!!  implicit none
    !!!  integer, intent(in) :: iproc,nspin
    !!!  real(gp), intent(in) :: hx,hy,hz
    !!!  type(orbitals_data), intent(in) :: orbs
    !!!  type(linearParameters):: lin
    !!!  type(locreg_descriptors), intent(in) :: lr
    !!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
    !!!  real(wp), dimension(*) :: pot
    !!!  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
    !!!  real(gp), intent(out) :: ekin_sum,epot_sum
    !!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
    !!!integer:: nat
    !!!real(8),dimension(3,nat):: rxyz
    !!!integer,dimension(orbs%norbp),intent(in):: onWhichAtom
    !!!type(atoms_data), intent(in) :: at
    !!!integer,intent(in),optional:: centralAtom
    !!!end subroutine local_hamiltonianConfinement


    !!!subroutine local_hamiltonianConfinementForAllLocregs(iproc,orbs,lin,lr,hx,hy,hz,&
    !!!     nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, at, centralAtom)
    !!!  use module_base
    !!!  use module_types
    !!!  use libxc_functionals
    !!!  implicit none
    !!!  integer, intent(in):: iproc, nspin, nat
    !!!  real(gp), intent(in):: hx,hy,hz
    !!!  type(orbitals_data),intent(in):: orbs
    !!!  type(linearParameters),intent(in):: lin
    !!!  type(locreg_descriptors),intent(in) :: lr
    !!!  type(atoms_data),intent(in):: at
    !!!  real(wp),dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp),intent(in):: psi
    !!!  real(wp),dimension(*):: pot
    !!!  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
    !!!  real(gp),intent(out):: ekin_sum,epot_sum
    !!!  real(wp),dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp,at%nat),intent(out):: hpsi
    !!!  real(8),dimension(3,nat),intent(in):: rxyz
    !!!  integer,intent(in),optional:: centralAtom
    !!!end subroutine local_hamiltonianConfinementForAllLocregs

    
    !!!subroutine deallocateLinear(iproc, lin, lphi, coeff)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc
    !!!  type(linearParameters),intent(inout):: lin
    !!!  real(8),dimension(:),pointer,intent(inout):: lphi
    !!!  real(8),dimension(:,:),pointer,intent(inout):: coeff
    !!!end subroutine deallocateLinear
    

    !!!subroutine initializeLocRegLIN(iproc, nproc, lr, lin, at, input, rxyz, radii_cf)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer:: iproc, nproc
    !!!  type(locreg_descriptors):: lr
    !!!  type(linearParameters):: lin
    !!!  type(atoms_data),intent(in):: at
    !!!  type(input_variables),intent(in):: input
    !!!  real(8),dimension(3,at%nat):: rxyz
    !!!  real(8),dimension(at%ntypes,3):: radii_cf
    !!!  type(communications_arrays):: commsLIN
    !!!end subroutine initializeLocRegLIN
    
    
    
    !!!subroutine orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, newComm, norbPerComm)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer, intent(in) :: iproc, lproc, uproc, newComm, norbPerComm
    !!!  type(linearParameters),intent(in out):: lin
    !!!end subroutine orbitalsCommunicatorsWithGroups
    
    subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,&
           rxyz,fion,fdisp,denspot,rhopotold,nlpspd,proj,GPU,&
           energs,scpot,energy)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(atoms_data),intent(inout):: at
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(inout):: rxyz
      real(8),dimension(3,at%nat),intent(in):: fion, fdisp
      type(DFT_local_fields), intent(inout) :: denspot
      real(gp), dimension(:), intent(inout) :: rhopotold
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      type(GPU_pointers),intent(in out):: GPU
      type(energy_terms),intent(inout) :: energs
      logical,intent(in):: scpot
      real(gp), dimension(:), pointer :: rho,pot
      real(8),intent(out):: energy
      type(DFT_wavefunction),intent(inout),target:: tmb
      type(DFT_wavefunction),intent(inout),target:: KSwfn
    end subroutine linearScaling   

   subroutine createDerivativeBasis(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     w_c, w_f, w_f1, w_f2, w_f3, x_c, x_f, y_c, y_f, z_c, z_f)
     use module_base
     !use filterModule
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

    subroutine readAtomicOrbitals(at,norbe,norbsc,nspin,nspinor,scorb,norbsc_arr,locrad)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: nspin,nspinor
      integer, intent(out) :: norbe,norbsc
      type(atoms_data), intent(in) :: at
      logical, dimension(4,2,at%natsc), intent(out) :: scorb
      integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
      real(gp), dimension(at%nat), intent(out) :: locrad
    end subroutine readAtomicOrbitals


    subroutine readAtomicOrbitals_withOnWhichAtom(at,orbsig,norbe,norbsc,nspin,nspinor,scorb,norbsc_arr,locrad,&
               onWhichAtom)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: nspin,nspinor
      type(orbitals_data),intent(in):: orbsig
      integer, intent(out) :: norbe,norbsc
      type(atoms_data), intent(inout) :: at
      logical, dimension(4,2,at%natsc), intent(out) :: scorb
      integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
      real(gp), dimension(at%nat), intent(out) :: locrad
      integer,dimension(orbsig%norb),intent(out):: onWhichAtom
    end subroutine readAtomicOrbitals_withOnWhichAtom

    subroutine inputguessConfinement(iproc, nproc, inputpsi, at, &
         input, hx, hy, hz, lzd, lorbs, rxyz,denspot, rhopotold,&
         nlpspd, proj, GPU,  &
         lphi,orbs,tmb,energs,overlapmatrix)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: iproc,nproc,inputpsi
      real(gp), intent(in) :: hx, hy, hz
      type(atoms_data), intent(inout) :: at
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(GPU_pointers), intent(inout) :: GPU
      type(DFT_local_fields), intent(inout) :: denspot
      type(input_variables),intent(in) :: input
      type(local_zone_descriptors),intent(inout):: lzd
      type(orbitals_data),intent(in):: lorbs
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
      real(dp),dimension(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout) ::  rhopotold
      real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out):: lphi
      type(orbitals_data),intent(inout):: orbs
      type(DFT_wavefunction),intent(inout):: tmb
      type(energy_terms),intent(out) :: energs
      real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: overlapmatrix
    end subroutine inputguessConfinement


    subroutine initialize_comms_sumrho(iproc,nproc,nscatterarr,lzd,orbs,comsr)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc,nproc
      integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs
      type(p2pComms),intent(out):: comsr
    end subroutine initialize_comms_sumrho


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

    !!$subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)
    !!$ use module_base
    !!$ use module_types
    !!$ implicit none
    !!$ integer, intent(in) :: alr,blr              
    !!$ integer, intent(in) :: nlr                  
    !!$ type(locreg_descriptors),intent(in) :: Glr  
    !!$ integer, intent(out) :: isovrlp           
    !!$ type(locreg_descriptors), dimension(nlr), intent(in) :: Llr       
    !!$end subroutine get_number_of_overlap_region

    !!$subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)
    !!$ use module_base
    !!$ use module_types
    !!$ implicit none
    !!$ integer, intent(in) :: alr,blr           
    !!$ integer, intent(in) :: nlr                
    !!$ type(locreg_descriptors),intent(in) :: Glr 
    !!$ integer, intent(in) :: isovrlp              
    !!$ type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  
    !!$ type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr 
    !!$end subroutine get_overlap_region_periodic

!!$    subroutine nlpspd_to_locreg(input_parameters,iproc,Glr,Llr,rxyz,atoms,orbs,&
!!$       radii_cf,cpmult,fpmult,hx,hy,hz,locregShape,nlpspd,Lnlpspd,projflg)
!!$     use module_base
!!$     use module_types
!!$     implicit none 
!!$     type(input_variables),intent(in) :: input_parameters
!!$     integer,intent(in) :: iproc
!!$     type(locreg_descriptors),intent(in) :: Glr  
!!$     type(locreg_descriptors),intent(in) :: Llr  
!!$     type(atoms_data),intent(in) :: atoms       
!!$     type(orbitals_data),intent(in) :: orbs      
!!$     real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz  
!!$     character(len=1),intent(in):: locregShape
!!$     type(nonlocal_psp_descriptors),intent(in) :: nlpspd  
!!$     type(nonlocal_psp_descriptors),intent(out) :: Lnlpspd   
!!$     integer,dimension(atoms%nat),intent(out) :: projflg
!!$     real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!$     real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
!!$    end subroutine nlpspd_to_locreg

!!$    subroutine apply_local_projectors(iorb,iproc,nspin,atoms,hx,hy,hz,Llr,Lnlpspd,orbs,projflg,psi,rxyz,hpsi,eproj)
!!$     use module_base
!!$     use module_types
!!$     implicit none
!!$     integer,intent(in) :: iorb,nspin,iproc
!!$     real(gp), intent(in) :: hx,hy,hz
!!$     type(atoms_data),intent(in) :: atoms
!!$     type(locreg_descriptors),intent(in) :: Llr
!!$     type(nonlocal_psp_descriptors),intent(in) :: Lnlpspd  ! Local descriptors for the projectors
!!$     type(orbitals_data),intent(in) :: orbs
!!$     real(gp), intent(inout) :: eproj
!!$     integer,dimension(atoms%nat),intent(in) :: projflg
!!$     real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*nspin),intent(in) :: psi  !local wavefunction
!!$     real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*nspin),intent(inout):: hpsi ! local |p><p|Psi>
!!$     real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!$    end subroutine apply_local_projectors


    subroutine psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)
     use module_base
     use module_types
     implicit none
     integer, intent(in) :: nlr    
     integer :: ilr           
     integer :: ldim          
     type(orbitals_data),intent(in) :: orbs      
     type(locreg_descriptors),intent(in) :: Glr  
     type(locreg_descriptors), dimension(nlr), intent(in) :: Olr   
     real(wp),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%norbp*orbs%nspinor),intent(in) :: psi      
     real(wp),dimension(ldim),intent(inout) :: lpsi 
    end subroutine psi_to_locreg


    subroutine psi_to_locreg2(iproc, nproc, ldim, gdim, Llr, Glr, gpsi, lpsi)
      use module_base
      use module_types
      implicit none
      integer,intent(in) :: iproc                  ! process ID
      integer,intent(in) :: nproc                  ! number of processes
      integer,intent(in) :: ldim          ! dimension of lpsi 
      integer,intent(in) :: gdim          ! dimension of gpsi 
      type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      real(wp),dimension(gdim),intent(in) :: gpsi       !Wavefunction (compressed format)
      real(wp),dimension(ldim),intent(out) :: lpsi   !Wavefunction in localization region
    end subroutine psi_to_locreg2



    subroutine partial_density_linear(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
         hfac,nscatterarr,spinsgn,psir,rho_p,&
         ibyyzz_r)
      use module_base
      use module_types
      implicit none
      logical, intent(in) :: rsflag
      integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
      real(gp), intent(in) :: hfac,spinsgn
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
      real(wp), dimension(n1i,n2i,n3i,npsir), intent(in) :: psir
      real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
      integer, dimension(:,:,:),pointer :: ibyyzz_r
    end subroutine partial_density_linear

    subroutine local_partial_densityLinear(nproc,rsflag,nscatterarr,&
         nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,mapping,psi,rho)
      use module_base
      use module_types
      use module_xc
      implicit none
      logical, intent(in) :: rsflag
      integer, intent(in) :: nproc
      integer,intent(inout):: nrhotot
      integer, intent(in) :: nspin
      real(gp), intent(in) :: hxh,hyh,hzh
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

     subroutine LinearDiagHam(iproc,at,etol,Lzd,orbs,nspin,natsc,Lhpsi,Lpsi,psit,orbsv,norbsc_arr)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc                                          
       integer, intent(in) :: nspin                                          
       integer, intent(in) :: natsc                                          
       real(gp),intent(in) :: etol         
       type(atoms_data),intent(in) :: at                                  
       type(local_zone_descriptors) :: Lzd                                  
       type(orbitals_data), intent(in) :: orbs                               
       type(orbitals_data), optional, intent(in) :: orbsv                    
       real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: Lhpsi               
       real(wp),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: Lpsi                
       real(wp),dimension(orbs%npsidim_comp),intent(inout):: psit                 
       integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr 
     end subroutine

     subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,Lzde,comms,&
          psi,hpsi,psit,orthpar,passmat,iscf,Tel,occopt,& !mandatory
          orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,natsc,nspin,occopt,iscf
       real(gp), intent(in) :: Tel
       type(local_zone_descriptors) :: Lzd        !> Information about the locregs after LIG
       type(local_zone_descriptors) :: Lzde       !> Informtation about the locregs for LIG
       type(communications_arrays), target, intent(in) :: comms
       type(orbitals_data), target, intent(inout) :: orbs
       type(orthon_data),intent(in):: orthpar 
       real(wp), dimension(*), intent(out) :: passmat !< passage matrix for building the eigenvectors (the size depends of the optional arguments)
       real(wp), dimension(:), pointer :: psi,hpsi,psit
       !optional arguments
       real(gp), optional, intent(in) :: etol
       type(orbitals_data), optional, intent(in) :: orbsv
       type(orbitals_data), optional, target, intent(inout) :: orbse
       type(communications_arrays), optional, target, intent(in) :: commse
       integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
       real(wp), dimension(:), pointer, optional :: psivirt
     end subroutine LDiagHam
     
     !!!subroutine getDerivativeBasisFunctions(iproc, nproc, hgrid, Glr, lin, nphi, phi, phid)
     !!!  use module_base
     !!!  use module_types
     !!!  implicit none
     !!!  integer,intent(in):: iproc, nproc, nphi
     !!!  real(8),intent(in):: hgrid
     !!!  type(locreg_descriptors),intent(in):: Glr
     !!!  type(linearParameters),intent(in):: lin
     !!!  real(8),dimension(nphi),intent(in):: phi
     !!!  real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)),intent(out):: phid
     !!!end subroutine getDerivativeBasisFunctions

     !!!subroutine orthonormalizeOnlyDerivatives(iproc, nproc, lin, phid)
     !!!  use module_base
     !!!  use module_defs
     !!!  use module_types
     !!!  implicit none
     !!!  integer,intent(in):: iproc, nproc
     !!!  type(linearParameters),intent(in):: lin
     !!!  real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)),intent(inout):: phid
     !!!end subroutine orthonormalizeOnlyDerivatives

!!$     subroutine getMatrixElements(iproc, nproc, Glr, orbs, comms, phi, hphi, matrixElements)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       integer,intent(in):: iproc, nproc
!!$       type(locreg_descriptors),intent(in):: Glr
!!$       type(orbitals_data),intent(in):: orbs
!!$       type(communications_arrays),intent(in):: comms
!!$       real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: phi, hphi
!!$       real(8),dimension(orbs%norb,orbs%norb,2),intent(out):: matrixElements
!!$     end subroutine getMatrixElements

        subroutine sumrhoForLocalizedBasis2(iproc,nproc,lzd,input,hx,hy,hz,orbs,&
             comsr,densKern,nrho,rho,at,nscatterarr)
          use module_base
          use module_types
          implicit none
          
          ! Calling arguments
          integer,intent(in):: iproc, nproc, nrho
          real(gp),intent(in):: hx, hy, hz
          type(local_zone_descriptors),intent(in):: lzd
          type(input_variables),intent(in):: input
          type(orbitals_data),intent(in):: orbs
          !type(p2pCommsSumrho),intent(inout):: comsr
          type(p2pComms),intent(inout):: comsr
          !real(8),dimension(orbs%norb,norb),intent(in):: coeff
          !real(8),dimension(ld_coeff,norb),intent(in):: coeff
          real(8),dimension(orbs%norb,orbs%norb),intent(in):: densKern
          real(8),dimension(nrho),intent(out),target:: rho
          type(atoms_data),intent(in):: at
          integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
        end subroutine sumrhoForLocalizedBasis2

     subroutine updatePotential(ixc,nspin,denspot,ehart,eexcu,vexcu)
     use module_base
     use module_types
     implicit none
     ! Calling arguments
     integer, intent(in) :: ixc,nspin
     type(DFT_local_fields), intent(inout) :: denspot
     real(8),intent(out):: ehart, eexcu, vexcu
   end subroutine updatePotential

     subroutine getIndices(lr, is1, ie1, is2, ie2, is3, ie3)
       use module_base
       use module_types
       implicit none
       type(locreg_descriptors),intent(in):: lr
       integer,intent(out):: is1, ie1, is2, ie2, is3, ie3
     end subroutine getIndices
     
     subroutine countOverlaps(iproc, nproc, orbs, lzd, op, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(out):: op
       type(p2pComms),intent(out):: comon
     end subroutine countOverlaps
     
     subroutine determineOverlaps(iproc, nproc, orbs, lzd, op, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(out):: op
       type(p2pComms),intent(out):: comon
     end subroutine determineOverlaps
     
     subroutine determineOverlapDescriptors(iproc, nproc, orbs, lzd, Glr, onWhichAtom, op)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       type(locreg_descriptors),intent(in):: Glr
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(overlapParameters),intent(inout):: op
     end subroutine determineOverlapDescriptors
     
     subroutine initCommsOrtho(iproc, nproc, nspin, hx, hy, hz, lzd, lzdig, orbs, &
                locregShape, bpo, op, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nspin
       real(8),intent(in):: hx, hy, hz
       type(local_zone_descriptors),intent(in):: lzd, lzdig
       type(orbitals_data),intent(in):: orbs
       character(len=1),intent(in):: locregShape
       type(basis_performance_options),intent(in):: bpo
       type(overlapParameters),intent(out):: op
       type(p2pComms),intent(out):: comon
     end subroutine initCommsOrtho
     
     subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
       integer,dimension(8),intent(out):: comarr
     end subroutine setCommsParameters
     
     
     !!subroutine postCommsOverlap(iproc, nproc, comon)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(p2pComms),intent(inout):: comon
     !!end subroutine postCommsOverlap
     
     
     !!subroutine extractOrbital(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, comon)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, sizePhi
     !!  type(orbitals_data),intent(in):: orbs
     !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
     !!  type(local_zone_descriptors),intent(in):: lzd
     !!  type(overlapParameters),intent(inout):: op
     !!  real(8),dimension(sizePhi),intent(in):: phi
     !!  type(p2pComms),intent(out):: comon
     !!end subroutine extractOrbital
     
     !!subroutine gatherOrbitals(iproc, nproc, comon)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(p2pComms),intent(inout):: comon
     !!end subroutine gatherOrbitals
     
     !!subroutine calculateOverlapMatrix(iproc, nproc, orbs, op, comon, onWhichAtom, lovrlp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(orbitals_data),intent(in):: orbs
     !!  type(overlapParameters),intent(in):: op
     !!  type(p2pComms),intent(inout):: comon
     !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
     !!  real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(out):: lovrlp
     !!end subroutine calculateOverlapMatrix
     
     
     !!subroutine calculateOverlapMatrix2(iproc, nproc, orbs, op, comon, onWhichAtom, mad, ovrlp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(orbitals_data),intent(in):: orbs
     !!  type(overlapParameters),intent(in):: op
     !!  type(p2pCommsOrthonormality),intent(inout):: comon
     !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
     !!  type(matrixDescriptors),intent(in):: mad
     !!  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
     !!end subroutine calculateOverlapMatrix2
     
     !!subroutine transformOverlapMatrix(iproc, nproc, comm, blocksize_dsyev, blocksize_pdgemm, norb, ovrlp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, comm, blocksize_dsyev, blocksize_pdgemm, norb
     !!  real(8),dimension(norb,norb),intent(inout):: ovrlp
     !!end subroutine transformOverlapMatrix
     
     !!subroutine expandOrbital(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(orbitals_data),intent(in):: orbs
     !!  type(input_variables),intent(in):: input
     !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
     !!  type(local_zone_descriptors),intent(in):: lzd
     !!  type(overlapParameters),intent(in):: op
     !!  type(p2pComms),intent(in):: comon
     !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
     !!end subroutine expandOrbital
     
     !!subroutine localGramschmidt(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, comon, lovrlp, lphiovrlp, lphi)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(orbitals_data),intent(in):: orbs, lorbs
     !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
     !!  type(local_zone_descriptors),intent(in):: lzd
     !!  type(overlapParameters),intent(in):: op
     !!  type(p2pComms),intent(in):: comon
     !!  real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(in):: lovrlp
     !!  real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
     !!  real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: lphi
     !!end subroutine localGramschmidt
     
     
     subroutine globalLoewdin(iproc, nproc, orbs, lzd, op, comon, ovrlp, lphiovrlp, lphi)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(in):: comon
       real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
       real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
       real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
     end subroutine globalLoewdin

     subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, &
                orbs, op, comon, lzd, mad, collcom, orthpar, bpo, lphi, psit_c, psit_f, &
                can_use_transposed)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc,nproc,methTransformOverlap,nItOrtho
       type(orbitals_data),intent(in):: orbs
       type(overlapParameters),intent(inout):: op
       type(p2pComms),intent(inout):: comon
       type(local_zone_descriptors),intent(in):: lzd
       type(matrixDescriptors),intent(in):: mad
       type(collective_comms),intent(in):: collcom
       type(orthon_data),intent(in):: orthpar
       type(basis_performance_options),intent(in):: bpo
       real(8),dimension(orbs%npsidim_orbs), intent(inout) :: lphi
       real(8),dimension(:),pointer:: psit_c, psit_f
       logical,intent(out):: can_use_transposed
     end subroutine orthonormalizeLocalized

     subroutine optimizeDIIS(iproc, nproc, orbs, lorbs, lzd, hphi, phi, ldiis, it)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, it
       type(orbitals_data),intent(in):: orbs, lorbs
       type(local_zone_descriptors),intent(in):: lzd
       real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(in):: hphi
       real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: phi
       type(localizedDIISParameters),intent(inout):: ldiis
     end subroutine optimizeDIIS

!!$     subroutine getHamiltonianMatrix(iproc, nproc, lzdig, Glr, input, onWhichAtom, onWhichAtomp, nat, chi, hchi, ham, orbsig)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       integer,intent(in):: iproc, nproc, nat
!!$       type(local_zone_descriptors),intent(in):: lzdig
!!$       type(locreg_descriptors),intent(in):: Glr
!!$       type(input_variables),intent(in):: input
!!$       type(orbitals_data),intent(in):: orbsig
!!$       integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!$       integer,dimension(orbsig%norbp),intent(in):: onWhichAtomp
!!$       !real(8),dimension(orbsig%npsidim),intent(in):: chi
!!$       !real(8),dimension(orbsig%npsidim,nat),intent(in):: hchi
!!$       real(8),dimension(orbsig%npsidim_comp),intent(in):: chi
!!$       real(8),dimension(orbsig%npsidim_comp,nat),intent(in):: hchi
!!$       real(8),dimension(orbsig%norb,orbsig%norb,nat),intent(out):: ham
!!$     end subroutine getHamiltonianMatrix

     subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, orbs, lzd, comgp, onWhichAtomAll, tag)
       use module_base
       use module_types
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
       implicit none
       integer,intent(in):: iproc, nproc
       integer,intent(inout):: tag
       type(orbitals_data),intent(in):: lorbs, llborbs
       type(local_zone_descriptors),intent(in):: lzd
       !type(p2pCommsRepartition),intent(out):: comrp
       type(p2pComms),intent(out):: comrp
     end subroutine initializeRepartitionOrbitals

     !!subroutine postCommsRepartition(iproc, nproc, orbs, comrp, nsendBuf, sendBuf, nrecvBuf, recvBuf)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
     !!  type(orbitals_data),intent(in):: orbs
     !!  !type(p2pCommsRepartition),intent(inout):: comrp
     !!  type(p2pComms),intent(inout):: comrp
     !!  real(8),dimension(nsendBuf),intent(in):: sendBuf
     !!  real(8),dimension(nrecvBuf),intent(out):: recvBuf
     !!end subroutine postCommsRepartition


     !!subroutine gatherDerivativeOrbitals(iproc, nproc, orbs, comrp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(orbitals_data),intent(in):: orbs
     !!  !type(p2pCommsRepartition),intent(inout):: comrp
     !!  type(p2pComms),intent(inout):: comrp
     !!end subroutine gatherDerivativeOrbitals


     subroutine getMatrixElements2(iproc, nproc, lzd, orbs, op_lb, comon_lb, lphi, lhphi, mad, matrixElements)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs
       type(overlapParameters),intent(inout):: op_lb
       type(p2pComms),intent(inout):: comon_lb
       real(8),dimension(orbs%npsidim_orbs),intent(in):: lphi, lhphi
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(orbs%norb,orbs%norb),intent(out):: matrixElements
     end subroutine getMatrixElements2



     subroutine determineLocalizationRegions(iproc, nproc, nlr, norb, at, onWhichAtomALl, &
                locrad, rxyz, lzd, lzdig, hx, hy, hz, mlr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norb
       type(atoms_data),intent(in):: at
       integer,dimension(norb),intent(in):: onWhichAtomAll
       real(8),dimension(at%nat),intent(in):: locrad
       real(8),dimension(3,at%nat),intent(in):: rxyz
       type(local_zone_descriptors),intent(in):: lzd, lzdig
       real(8),intent(in):: hx, hy, hz
       type(matrixLocalizationRegion),dimension(:),pointer,intent(out):: mlr
     end subroutine determineLocalizationRegions


     subroutine extractMatrix(iproc, nproc, norb, norbp, orbstot, onWhichAtomPhi, onWhichMPI, nmat, ham, matmin, hamextract)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nmat, norb, norbp
       type(orbitals_data),intent(in):: orbstot
       integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
       real(8),dimension(orbstot%norb,orbstot%norb,nmat),intent(in):: ham
       type(matrixMinimization),intent(inout):: matmin
       real(8),dimension(:,:,:),pointer,intent(out):: hamextract
     end subroutine extractMatrix


     subroutine orthonormalizeVectors(iproc, nproc, comm, nItOrtho, methTransformOverlap, &
                orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, &
                mad, mlr, vec, opm, comom, collcom, orthpar, bpo)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, comm, nItOrtho, methTransformOverlap
       integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       type(matrixDescriptors),intent(in):: mad
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       real(8),dimension(norbmax,norbp),intent(inout):: vec
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(inout):: comom
       type(collective_comms),intent(in):: collcom
       type(orthon_data),intent(in):: orthpar
       type(basis_performance_options),intent(in):: bpo
     end subroutine orthonormalizeVectors




     subroutine initCommsMatrixOrtho(iproc, nproc, norb, norb_par, isorb_par, onWhichAtomPhi, onWhichMPI, opm, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, norb
       integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
       integer,dimension(nproc),intent(in):: norb_par, isorb_par
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(inout):: comom
     end subroutine initCommsMatrixOrtho



     subroutine determineOverlapRegionMatrix(iproc, nproc, lzd, mlr, orbs, orbstot, onWhichAtom, onWhichAtomPhi, comom, opm)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs, orbstot
       integer,dimension(orbstot%norb),intent(in):: onWhichAtom
       integer,dimension(orbs%norb),intent(in):: onWhichAtomPhi
       type(matrixLocalizationRegion),dimension(lzd%nlr),intent(in):: mlr
       type(p2pComms),intent(out):: comom
       type(overlap_parameters_matrix),intent(out):: opm
     end subroutine determineOverlapRegionMatrix


     !!subroutine postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
     !!use module_base
     !!use module_types
     !!implicit none
     !!integer,intent(in):: iproc, nproc, newComm
     !!type(p2pComms),intent(inout):: comom
     !!end subroutine postCommsVectorOrthonormalization


     !!subroutine gatherVectors(iproc, nproc, newComm, comom)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, newComm
     !!  type(p2pComms),intent(inout):: comom
     !!end subroutine gatherVectors


     subroutine extractToOverlapregion(iproc, nproc, norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, opm, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, norbmax, norb, norbp
       integer,dimension(norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       real(8),dimension(norbmax,norbp),intent(in):: vec
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(inout):: comom
     end subroutine extractToOverlapregion



     subroutine expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, opm, comom, norbmax, noverlaps, vecOvrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, isorb, norbp, norbmax, noverlaps
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(in):: comom
       real(8),dimension(norbmax,noverlaps),intent(out):: vecOvrlp
     end subroutine expandFromOverlapregion

     subroutine calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, opm, comom, mlr,&
                onWhichAtom, vec, vecOvrlp, newComm, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, newComm
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(in):: comom
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       integer,dimension(norb),intent(in):: onWhichAtom
       real(8),dimension(norbmax,norbp),intent(in):: vec
       real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
       real(8),dimension(norb,norb),intent(out):: ovrlp
     end subroutine calculateOverlap


     subroutine orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, opm, comom, &
               mlr, onWhichAtom, vecOvrlp, ovrlp, vec)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(in):: comom
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       integer,dimension(norb),intent(in):: onWhichAtom
       real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
       real(8),dimension(norb,norb),intent(in):: ovrlp
       real(8),dimension(norbmax,norbp),intent(inout):: vec
     end subroutine orthonormalLinearCombinations


     !!!subroutine buildLinearCombinationsLocalized(iproc, nproc, orbsig, orbs, comms, at, Glr, input, norbsPerType, &
     !!!      onWhichAtom, lchi, lphi, rxyz, onWhichAtomPhi, lin, lzdig, ham)
     !!!  use module_base
     !!!  use module_types
     !!!  implicit none
     !!!  integer,intent(in):: iproc, nproc
     !!!  type(orbitals_data),intent(in):: orbsig, orbs
     !!!  type(communications_arrays),intent(in):: comms
     !!!  type(atoms_data),intent(in):: at
     !!!  type(locreg_descriptors),intent(in):: Glr
     !!!  type(input_variables),intent(in):: input
     !!!  type(linearParameters),intent(in):: lin
     !!!  type(local_zone_descriptors),intent(inout):: lzdig
     !!!  integer,dimension(at%ntypes):: norbsPerType
     !!!  integer,dimension(orbsig%norb),intent(in):: onWhichAtom
     !!!  real(8),dimension(max(orbsig%npsidim_comp,orbsig%npsidim_orbs)):: lchi
     !!!  real(8),dimension(max(lin%orbs%npsidim_comp,lin%orbs%npsidim_orbs)):: lphi
     !!!  real(8),dimension(3,at%nat):: rxyz
     !!!  integer,dimension(orbs%norb):: onWhichAtomPhi
     !!!  real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(inout):: ham
     !!!end subroutine buildLinearCombinationsLocalized




     subroutine orthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, &
                orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, mad, vec, grad, &
                opm, comom, trace, collcom, orthpar, bpo)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint
       integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(norbmax,norbp),intent(inout):: vec, grad
       type(overlap_parameters_matrix),intent(in):: opm
       type(p2pComms),intent(inout):: comom
       real(8),intent(out):: trace
       type(collective_comms),intent(in):: collcom
       type(orthon_data),intent(in):: orthpar
       type(basis_performance_options),intent(in):: bpo
     end subroutine orthoconstraintVectors


     !!subroutine cubic_exact_exchange(iproc,nproc,nspin,npsidim,size_potxc,hx,hy,hz,Glr,orbs,&
     !!           ngatherarr,psi,potxc,eexctX,pkernel,orbsocc,psirocc)
     !!  use module_base
     !!  use module_types
     !!  use module_xc
     !!  implicit none
     !!  integer, intent(in) :: iproc,nproc,nspin,npsidim,size_potxc
     !!  real(gp), intent(in) :: hx,hy,hz
     !!  type(locreg_descriptors) :: Glr
     !!  type(orbitals_data) :: orbs
     !!  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
     !!  real(wp), dimension(npsidim), intent(in) :: psi
     !!  real(wp), dimension(size_potxc),intent(out) :: potxc
     !!  real(gp), intent(out) :: eexctX
     !!  real(dp), dimension(*), optional :: pkernel
     !!  type(orbitals_data), intent(in), optional :: orbsocc
     !!  real(wp), dimension(:), pointer, optional :: psirocc
     !!end subroutine


!!$     subroutine getHamiltonianMatrix2(iproc, nproc, lzdig, orbsig, Glr, input, onWhichAtom, onWhichAtomp, nat, lchi, lhchi, ham)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       integer,intent(in):: iproc, nproc, nat
!!$       type(local_zone_descriptors),intent(in):: lzdig
!!$       type(orbitals_data),intent(in):: orbsig
!!$       type(locreg_descriptors),intent(in):: Glr
!!$       type(input_variables),intent(in):: input
!!$       integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!$       integer,dimension(orbsig%norbp),intent(in):: onWhichAtomp
!!$       real(8),dimension(orbsig%npsidim_comp),intent(in):: lchi
!!$       real(8),dimension(orbsig%npsidim_comp,nat),intent(in):: lhchi
!!$       real(8),dimension(orbsig%norb,orbsig%norb,nat),intent(out):: ham
!!$     end subroutine getHamiltonianMatrix2
!!$

     subroutine getDerivativeBasisFunctions(iproc, nproc, hgrid, lzd, lorbs, lborbs, comrp, nphi, phi, phid)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nphi
       real(8),intent(in):: hgrid
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: lorbs, lborbs
       !type(p2pCommsRepartition),intent(inout):: comrp
       type(p2pComms),intent(inout):: comrp
       real(8),dimension(nphi),intent(in):: phi
       real(8),dimension(max(lborbs%npsidim_orbs,lborbs%npsidim_comp)),target,intent(out):: phid
     end subroutine getDerivativeBasisFunctions


     !!subroutine buildLinearCombinations(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, locregShape, &
     !!           tag, comonig, opig, madig, lphi)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(local_zone_descriptors),intent(in):: lzdig, lzd
     !!  type(orbitals_data),intent(in):: orbsig, orbs
     !!  type(input_variables),intent(in):: input
     !!  real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
     !!  real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
     !!  character(len=1),intent(in):: locregShape
     !!  integer,intent(inout):: tag
     !!  type(p2pComms):: comonig
     !!  type(overlapParameters):: opig
     !!  type(matrixDescriptors):: madig
     !!  real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: lphi
     !!end subroutine buildLinearCombinations

     subroutine buildLinearCombinations_new(iproc, nproc, lzdig, lzd, orbsig, orbs, coeff, lchi, &
                collcomig, collcom, lphi)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzdig, lzd
       type(orbitals_data),intent(in):: orbsig, orbs
       real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
       real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
       type(collective_comms),intent(in):: collcomig, collcom
       real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
     end subroutine buildLinearCombinations_new

     !!subroutine postCommunicationsPotential(iproc, nproc, ndimpot, pot, comgp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, ndimpot
     !!  real(8),dimension(ndimpot),intent(in):: pot
     !!  !type(p2pCommsGatherPot),intent(inout):: comgp
     !!  type(p2pComms),intent(inout):: comgp
     !!end subroutine postCommunicationsPotential


     subroutine gatherPotential(iproc, nproc, comgp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       !type(p2pCommsGatherPot),intent(inout):: comgp
       type(p2pComms),intent(inout):: comgp
     end subroutine gatherPotential


     subroutine extractOrbital2(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, sizePhi
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(inout):: op
       real(8),dimension(sizePhi),intent(in):: phi
       type(p2pComms),intent(out):: comon
     end subroutine extractOrbital2


     !!subroutine gatherOrbitals2(iproc, nproc, comon)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc
     !!  type(p2pComms),intent(inout):: comon
     !!end subroutine gatherOrbitals2


     !!!subroutine expandOrbital2(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
     !!!  use module_base
     !!!  use module_types
     !!!  implicit none
     !!!  
     !!!  ! Calling arguments
     !!!  integer,intent(in):: iproc, nproc
     !!!  type(orbitals_data),intent(in):: orbs
     !!!  type(input_variables),intent(in):: input
     !!!  type(local_zone_descriptors),intent(in):: lzd
     !!!  type(overlapParameters),intent(in):: op
     !!!  type(p2pComms),intent(in):: comon
     !!!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
     !!!end subroutine expandOrbital2

     !!!subroutine getOverlapMatrix(iproc, nproc, lin, input, lphi, mad, ovrlp)
     !!!  use module_base
     !!!  use module_types
     !!!  implicit none
     !!!  integer,intent(in):: iproc, nproc
     !!!  type(linearParameters),intent(inout):: lin
     !!!  type(input_variables),intent(in):: input
     !!!  real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
     !!!  type(matrixDescriptors),intent(in):: mad
     !!!  real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(out):: ovrlp
     !!!end subroutine getOverlapMatrix


     subroutine getOverlapMatrix2(iproc, nproc, lzd, orbs, comon_lb, op_lb, lphi, mad, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs
       type(p2pComms),intent(inout):: comon_lb
       type(overlapParameters),intent(inout):: op_lb
       real(8),dimension(orbs%npsidim_orbs),intent(inout):: lphi
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
     end subroutine getOverlapMatrix2



     subroutine mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, alphaMix, mixMeth, pnrm)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, ndimpot, ndimtot, mixMeth
       real(8),dimension(ndimpot),intent(in):: rhopotold
       real(8),dimension(ndimpot),intent(out):: rhopot
       type(mixrhopotDIISParameters),intent(inout):: mixdiis
       real(8),intent(in):: alphaMix
       real(8),intent(out):: pnrm
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


     subroutine allocateCommunicationbufferSumrho(iproc, comsr, subname)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc
       type(p2pComms),intent(inout):: comsr
       character(len=*),intent(in):: subname
     end subroutine allocateCommunicationbufferSumrho


     subroutine deallocateCommunicationbufferSumrho(comsr, subname)
       use module_base
       use module_types
       implicit none
       !type(p2pCommsSumrho),intent(inout):: comsr
       type(p2pComms),intent(inout):: comsr
       character(len=*),intent(in):: subname
     end subroutine deallocateCommunicationbufferSumrho


     subroutine allocateCommuncationBuffersOrtho(comon, subname)
       use module_base
       use module_types
       implicit none
       type(p2pComms),intent(inout):: comon
       character(len=*),intent(in):: subname
     end subroutine allocateCommuncationBuffersOrtho


     subroutine deallocateCommuncationBuffersOrtho(comon, subname)
       use module_base
       use module_types
       implicit none
       type(p2pComms),intent(inout):: comon
       character(len=*),intent(in):: subname
     end subroutine deallocateCommuncationBuffersOrtho


     subroutine allocateCommunicationsBuffersPotential(comgp, subname)
       use module_base
       use module_types
       implicit none
       !type(p2pCommsGatherPot),intent(inout):: comgp
       type(p2pComms),intent(inout):: comgp
       character(len=*),intent(in):: subname
     end subroutine allocateCommunicationsBuffersPotential


     subroutine deallocateCommunicationsBuffersPotential(comgp, subname)
       use module_base
       use module_types
       implicit none
       !type(p2pCommsGatherPot),intent(inout):: comgp
       type(p2pComms),intent(inout):: comgp
       character(len=*),intent(in):: subname
     end subroutine deallocateCommunicationsBuffersPotential


     subroutine copy_locreg_descriptors(glrin, glrout, subname)
       use module_base
       use module_types
       implicit none
       type(locreg_descriptors),intent(in):: glrin
       type(locreg_descriptors),intent(out):: glrout
       character(len=*),intent(in):: subname
     end subroutine copy_locreg_descriptors


     subroutine copy_grid_dimensions(din, dout)
       use module_base
       use module_types
       implicit none
       type(grid_dimensions),intent(in):: din
       type(grid_dimensions),intent(out):: dout
     end subroutine copy_grid_dimensions


     subroutine copy_wavefunctions_descriptors(wfdin, wfdout, subname)
       use module_base
       use module_types
       implicit none
       type(wavefunctions_descriptors),intent(in):: wfdin
       type(wavefunctions_descriptors),intent(out):: wfdout
       character(len=*),intent(in):: subname
     end subroutine copy_wavefunctions_descriptors


     subroutine copy_convolutions_bounds(geocode,boundsin, boundsout, subname)
       use module_base
       use module_types
       implicit none
       character(len=1),intent(in) :: geocode
       type(convolutions_bounds),intent(in):: boundsin
       type(convolutions_bounds),intent(out):: boundsout
       character(len=*),intent(in):: subname
     end subroutine copy_convolutions_bounds


     subroutine copy_kinetic_bounds(geocode,kbin, kbout, subname)
       use module_base
       use module_types
       implicit none
       character(len=1),intent(in) :: geocode
       type(kinetic_bounds),intent(in):: kbin
       type(kinetic_bounds),intent(out):: kbout
       character(len=*),intent(in):: subname
     end subroutine copy_kinetic_bounds


     subroutine copy_shrink_bounds(geocode,sbin, sbout, subname)
       use module_base
       use module_types
       implicit none
       character(len=1),intent(in) :: geocode
       type(shrink_bounds),intent(in):: sbin
       type(shrink_bounds),intent(out):: sbout
       character(len=*),intent(in):: subname
     end subroutine copy_shrink_bounds


     subroutine copy_grow_bounds(geocode,gbin, gbout, subname)
       use module_base
       use module_types
       implicit none
       character(len=1),intent(in) :: geocode
       type(grow_bounds),intent(in):: gbin
       type(grow_bounds),intent(out):: gbout
       character(len=*),intent(in):: subname
     end subroutine copy_grow_bounds


     subroutine copy_nonlocal_psp_descriptors(nlpspin, nlpspout, subname)
       use module_base
       use module_types
       implicit none
       type(nonlocal_psp_descriptors),intent(in):: nlpspin
       type(nonlocal_psp_descriptors),intent(out):: nlpspout
       character(len=*),intent(in):: subname
     end subroutine copy_nonlocal_psp_descriptors


     subroutine copy_orbitals_data(orbsin, orbsout, subname)
       use module_base
       use module_types
       implicit none
       type(orbitals_data),intent(in):: orbsin
       type(orbitals_data),intent(out):: orbsout
       character(len=*),intent(in):: subname
     end subroutine copy_orbitals_data


     !!subroutine deallocate_matrixLocalizationRegion(mlr, subname)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  type(matrixLocalizationRegion),intent(inout):: mlr
     !!  character(len=*),intent(in):: subname
     !!end subroutine deallocate_matrixLocalizationRegion


     !!subroutine deallocate_matrixMinimization(matmin, subname)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  type(matrixMinimization),intent(inout):: matmin
     !!  character(len=*),intent(in):: subname
     !!end subroutine deallocate_matrixMinimization


     !!subroutine deallocate_local_zone_descriptors(lzd, subname)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  type(local_zone_descriptors),intent(inout):: lzd
     !!  character(len=*),intent(in):: subname
     !!end subroutine deallocate_local_zone_descriptors


    subroutine deallocate_local_zone_descriptors(lzd, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(local_zone_descriptors),intent(inout):: lzd
      character(len=*),intent(in):: subname
    end subroutine deallocate_local_zone_descriptors

    subroutine deallocate_Lzd_except_Glr(lzd, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(local_zone_descriptors),intent(inout):: lzd
      character(len=*),intent(in):: subname
    end subroutine deallocate_Lzd_except_Glr


    subroutine deallocate_orbitals_data(orbs, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(orbitals_data),intent(inout):: orbs
      character(len=*),intent(in):: subname
    end subroutine deallocate_orbitals_data

    subroutine deallocate_communications_arrays(comms, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(communications_arrays),intent(inout):: comms
      character(len=*),intent(in):: subname
    end subroutine deallocate_communications_arrays

    subroutine deallocate_locreg_descriptors(lr, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(locreg_descriptors),intent(inout):: lr
      character(len=*),intent(in):: subname
    end subroutine deallocate_locreg_descriptors

    subroutine deallocate_wavefunctions_descriptors(wfd, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(wavefunctions_descriptors),intent(inout):: wfd
      character(len=*),intent(in):: subname
    end subroutine deallocate_wavefunctions_descriptors

    subroutine deallocate_convolutions_bounds(bounds, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(convolutions_bounds),intent(inout):: bounds
      character(len=*),intent(in):: subname
    end subroutine deallocate_convolutions_bounds

    subroutine deallocate_kinetic_bounds(kb, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(kinetic_bounds),intent(inout):: kb
      character(len=*),intent(in):: subname
    end subroutine deallocate_kinetic_bounds

    subroutine deallocate_shrink_bounds(sb, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(shrink_bounds),intent(inout):: sb
      character(len=*),intent(in):: subname
    end subroutine deallocate_shrink_bounds

    subroutine deallocate_grow_bounds(gb, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(grow_bounds),intent(inout):: gb
      character(len=*),intent(in):: subname
    end subroutine deallocate_grow_bounds

    subroutine deallocate_nonlocal_psp_descriptors(nlpspd, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(nonlocal_psp_descriptors),intent(inout):: nlpspd
      character(len=*),intent(in):: subname
    end subroutine deallocate_nonlocal_psp_descriptors

    subroutine deallocate_matrixMinimization(matmin, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(matrixMinimization),intent(inout):: matmin
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixMinimization

    subroutine deallocate_matrixLocalizationRegion(mlr, subname)
      use module_base
      use module_types
      !use deallocatePointers
      implicit none
      type(matrixLocalizationRegion),intent(inout):: mlr
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixLocalizationRegion

    !!!subroutine nullify_linearParameters(lin)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  type(linearParameters),intent(out):: lin
    !!!end subroutine nullify_linearParameters

    !!subroutine nullify_p2pCommsSumrho(comsr)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  !type(p2pCommsSumrho),intent(out):: comsr
    !!  type(p2pComms),intent(out):: comsr
    !!end subroutine nullify_p2pCommsSumrho

    !!subroutine nullify_p2pCommsGatherPot(comgp)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pCommsGatherPot),intent(out):: comgp
    !!end subroutine nullify_p2pCommsGatherPot

    !!subroutine nullify_largeBasis(lb)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(largeBasis),intent(out):: lb
    !!end subroutine nullify_largeBasis

    !!subroutine nullify_p2pCommsRepartition(comrp)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pCommsRepartition),intent(out):: comrp
    !!end subroutine nullify_p2pCommsRepartition

    !!subroutine nullify_p2pCommsOrthonormality(comon)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pCommsOrthonormality),intent(out):: comon
    !!end subroutine nullify_p2pCommsOrthonormality

    subroutine nullify_overlapParameters(op)
      use module_base
      use module_types
      implicit none
      type(overlapParameters),intent(out):: op
    end subroutine nullify_overlapParameters

    !!!subroutine nullify_linearInputGuess(lig)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  type(linearInputGuess),intent(out):: lig
    !!!end subroutine nullify_linearInputGuess

    subroutine nullify_matrixDescriptors(mad)
      use module_base
      use module_types
      implicit none
      type(matrixDescriptors),intent(out):: mad
    end subroutine nullify_matrixDescriptors

    subroutine nullify_local_zone_descriptors(lzd)
      use module_base
      use module_types
      implicit none
      type(local_zone_descriptors),intent(out):: lzd
    end subroutine nullify_local_zone_descriptors
    
    subroutine nullify_orbitals_data(orbs)
      use module_base
      use module_types
      implicit none
      type(orbitals_data),intent(out):: orbs
    end subroutine nullify_orbitals_data
    
    subroutine nullify_communications_arrays(comms)
      use module_base
      use module_types
      implicit none
      type(communications_arrays),intent(out):: comms
    end subroutine nullify_communications_arrays
    
    subroutine nullify_locreg_descriptors(lr)
      use module_base
      use module_types
      implicit none
      type(locreg_descriptors),intent(out):: lr
    end subroutine nullify_locreg_descriptors
    
    subroutine nullify_wavefunctions_descriptors(wfd)
      use module_base
      use module_types
      implicit none
      type(wavefunctions_descriptors),intent(out):: wfd
    end subroutine nullify_wavefunctions_descriptors
    
    subroutine nullify_convolutions_bounds(bounds)
      use module_base
      use module_types
      implicit none
      type(convolutions_bounds),intent(out):: bounds
    end subroutine nullify_convolutions_bounds
    
    subroutine nullify_kinetic_bounds(kb)
      use module_base
      use module_types
      implicit none
      type(kinetic_bounds),intent(out):: kb
    end subroutine nullify_kinetic_bounds
    
    subroutine nullify_shrink_bounds(sb)
      use module_base
      use module_types
      implicit none
      type(shrink_bounds),intent(out):: sb
    end subroutine nullify_shrink_bounds
    
    subroutine nullify_grow_bounds(gb)
      use module_base
      use module_types
      implicit none
      type(grow_bounds),intent(out):: gb
    end subroutine nullify_grow_bounds
    
    subroutine nullify_nonlocal_psp_descriptors(nlpspd)
      use module_base
      use module_types
      implicit none
      type(nonlocal_psp_descriptors),intent(out):: nlpspd
    end subroutine nullify_nonlocal_psp_descriptors

    subroutine nullify_matrixMinimization(matmin)
      use module_base
      use module_types
      implicit none
      type(matrixMinimization),intent(out):: matmin
    end subroutine nullify_matrixMinimization

    subroutine nullify_matrixLocalizationRegion(mlr)
      use module_base
      use module_types
      implicit none
      type(matrixLocalizationRegion),intent(out):: mlr
    end subroutine nullify_matrixLocalizationRegion

    subroutine initLocregs(iproc, nproc, nlr, rxyz, hx, hy, hz, lzd, orbs, Glr, locrad, locregShape, lborbs)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, nlr
      real(8),dimension(3,nlr),intent(in):: rxyz
      real(8),intent(in):: hx, hy, hz
      type(local_zone_descriptors),intent(inout):: lzd
      type(orbitals_data),intent(in):: orbs
      type(locreg_descriptors),intent(in):: Glr
      real(8),dimension(lzd%nlr),intent(in):: locrad
      character(len=1),intent(in):: locregShape
      type(orbitals_data),optional,intent(in):: lborbs
    end subroutine initLocregs

    !!!subroutine deallocate_linearParameters(lin, subname)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  ! Calling arguments
    !!!  type(linearParameters),intent(inout):: lin
    !!!  character(len=*),intent(in):: subname
    !!!end subroutine deallocate_linearParameters

    !!subroutine deallocate_largeBasis(lb, subname)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(largeBasis),intent(inout):: lb
    !!  character(len=*),intent(in):: subname
    !!end subroutine deallocate_largeBasis
    
    !!subroutine deallocate_p2pCommsRepartition(comrp, subname)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pCommsRepartition),intent(inout):: comrp
    !!  character(len=*),intent(in):: subname
    !!end subroutine deallocate_p2pCommsRepartition
    
    !!subroutine deallocate_p2pCommsOrthonormality(comon, subname)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pCommsOrthonormality),intent(inout):: comon
    !!  character(len=*),intent(in):: subname
    !!end subroutine deallocate_p2pCommsOrthonormality
    
    subroutine deallocate_overlapParameters(op, subname)
      use module_base
      use module_types
      implicit none
      type(overlapParameters),intent(inout):: op
      character(len=*),intent(in):: subname
    end subroutine deallocate_overlapParameters

    !!subroutine deallocate_inguessParameters(ip, subname)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(inguessParameters),intent(inout):: ip
    !!  character(len=*),intent(in):: subname
    !!end subroutine deallocate_inguessParameters

    !!subroutine deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  type(p2pComms),intent(inout):: comom
    !!  character(len=*),intent(in):: subname
    !!end subroutine deallocate_p2pCommsOrthonormalityMatrix

    subroutine deallocate_matrixDescriptors(mad, subname)
      use module_base
      use module_types
      implicit none
      type(matrixDescriptors),intent(inout):: mad
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixDescriptors


    !!!subroutine cancelCommunicationPotential(iproc, nproc, comgp)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  !type(p2pCommsGatherPot),intent(inout):: comgp
    !!!  type(p2pComms),intent(inout):: comgp
    !!!end subroutine cancelCommunicationPotential

    !!!subroutine initCommsOrthoVariable(iproc, nproc, lzd, orbs, orbsig, onWhichAtomAll, input, op, comon, tag)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  type(local_zone_descriptors),intent(in):: lzd
    !!!  type(orbitals_data),intent(in):: orbs, orbsig
    !!!  integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
    !!!  type(input_variables),intent(in):: input
    !!!  type(overlapParameters),intent(out):: op
    !!!  type(p2pComms),intent(out):: comon
    !!!  integer,intent(inout):: tag
    !!!end subroutine initCommsOrthoVariable
    !!!
    !!!subroutine countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  type(orbitals_data),intent(in):: orbs, orbsig
    !!!  type(local_zone_descriptors),intent(in):: lzd
    !!!  type(overlapParameters),intent(out):: op
    !!!  type(p2pComms),intent(out):: comon
    !!!end subroutine countOverlapsVariable
    !!!
    !!!subroutine determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  type(orbitals_data),intent(in):: orbs, orbsig
    !!!  type(local_zone_descriptors),intent(in):: lzd
    !!!  type(overlapParameters),intent(out):: op
    !!!  type(p2pComms),intent(out):: comon
    !!!end subroutine determineOverlapsVariable
    !!!
    !!!subroutine determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, Glr, onWhichAtom, op)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  type(orbitals_data),intent(in):: orbs, orbsig
    !!!  type(local_zone_descriptors),intent(in):: lzd
    !!!  type(locreg_descriptors),intent(in):: Glr
    !!!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
    !!!  type(overlapParameters),intent(inout):: op
    !!!end subroutine determineOverlapDescriptorsVariable
    !!!
    !!!subroutine setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  integer,intent(in):: iproc, nproc
    !!!  type(orbitals_data),intent(in):: orbs, orbsig
    !!!  type(local_zone_descriptors),intent(in):: lzd
    !!!  type(overlapParameters),intent(inout):: op
    !!!  type(p2pComms),intent(out):: comon
    !!!  integer,intent(inout):: tag
    !!!end subroutine setCommsOrthoVariable
    
    !!subroutine indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer,intent(in):: iproc, nproc
    !!  type(orbitals_data),intent(in):: orbs
    !!  type(input_variables),intent(in):: input
    !!  type(local_zone_descriptors),intent(in):: lzd
    !!  type(overlapParameters),intent(in):: op
    !!  type(p2pComms),intent(in):: comon
    !!end subroutine indicesForExpansionVariable
    
    !!subroutine indicesForExtractionVariable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, comon)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer,intent(in):: iproc, nproc, sizePhi
    !!  type(orbitals_data),intent(in):: orbs, orbsig
    !!  type(local_zone_descriptors),intent(in):: lzd
    !!  type(overlapParameters),intent(inout):: op
    !!  type(p2pComms),intent(out):: comon
    !!end subroutine indicesForExtractionVariable
    
    !!subroutine extractOrbital2Variable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, phi, comon)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer,intent(in):: iproc, nproc, sizePhi
    !!  type(orbitals_data),intent(in):: orbs, orbsig
    !!  type(local_zone_descriptors),intent(in):: lzd
    !!  type(overlapParameters),intent(inout):: op
    !!  real(8),dimension(sizePhi),intent(in):: phi
    !!  type(p2pComms),intent(out):: comon
    !!end subroutine extractOrbital2Variable
    
    !!subroutine expandOrbital2Variable(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer,intent(in):: iproc, nproc
    !!  type(orbitals_data),intent(in):: orbs
    !!  type(input_variables),intent(in):: input
    !!  type(local_zone_descriptors),intent(in):: lzd
    !!  type(overlapParameters),intent(in):: op
    !!  type(p2pComms),intent(in):: comon
    !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
    !!end subroutine expandOrbital2Variable

    !!subroutine buildLinearCombinationsVariable(iproc, nproc, lzdig, lzd, orbsig, &
    !!           orbs, input, coeff, lchi, tag, lphi)
    !!  use module_base
    !!  use module_types
    !!  implicit none

    !!  ! Calling arguments
    !!  integer,intent(in):: iproc, nproc
    !!  type(local_zone_descriptors),intent(in):: lzdig, lzd
    !!  type(orbitals_data),intent(in):: orbsig, orbs
    !!  type(input_variables),intent(in):: input
    !!  real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
    !!  real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
    !!  integer,intent(inout):: tag
    !!  real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
    !!end subroutine buildLinearCombinationsVariable

    !subroutine index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, indexLpsi)
    !  use module_base
    !  use module_types
    !  implicit none
    !  integer,intent(in):: iproc, nproc
    !  integer :: Gdim          ! dimension of psi 
    !  integer :: Ldim          ! dimension of lpsi
    !  integer :: norb          ! number of orbitals
    !  integer :: nspinor       ! number of spinors
    !  integer :: nspin         ! number of spins 
    !  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
    !  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
    !  integer,dimension(Ldim),intent(out) :: indexLpsi         !Wavefunction in localization region
    !end subroutine index_of_Lpsi_to_global2


     subroutine initInputguessConfinement(iproc, nproc, at, lzd, orbs, collcom_reference, &
                Glr, input, hx, hy, hz, lin, tmbig, tmbgauss, rxyz, nscatterarr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc,nproc
       real(gp), intent(in) :: hx, hy, hz
       type(atoms_data),intent(inout) :: at
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs
       type(collective_comms),intent(in):: collcom_reference
       type(locreg_descriptors),intent(in) :: Glr
       type(input_variables), intent(in) ::input
       type(linearInputParameters),intent(in):: lin
       type(DFT_wavefunction),intent(out):: tmbig, tmbgauss
       integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       real(gp),dimension(3,at%nat),intent(in):: rxyz
     end subroutine initInputguessConfinement


    subroutine orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, methTransformOverlap, nItOrtho, &
               lzd, orbs, comon, op, input, mad, collcom, orthpar, bpo, lchi, can_use_transposed)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs
      type(input_variables),intent(in):: input
      type(p2pComms),intent(inout):: comon
      type(overlapParameters),intent(inout):: op
      type(matrixDescriptors),intent(in):: mad
      type(collective_comms),intent(in):: collcom
      type(orthon_data),intent(in):: orthpar
      type(basis_performance_options),intent(in):: bpo
      real(8),dimension(orbs%npsidim_orbs),intent(inout):: lchi
      logical,intent(inout):: can_use_transposed
    end subroutine orthonormalizeAtomicOrbitalsLocalized2

    subroutine build_input_guess(iproc, nproc, nlocregPerMPI, hx, hy, hz, &
               tmb, tmbig, at, input, lchi, locregCenter, rxyz, ham, lphi)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, nlocregPerMPI
      real(gp), intent(in) :: hx, hy, hz
      type(DFT_wavefunction),intent(in):: tmb, tmbig
      type(atoms_data),intent(in):: at
      type(input_variables),intent(in):: input
      real(8),dimension(tmbig%orbs%npsidim_orbs),intent(in):: lchi
      real(8),dimension(3,tmbig%lzd%nlr),intent(in):: locregCenter
      real(8),dimension(3,at%nat),intent(in):: rxyz
      real(8),dimension(tmbig%orbs%norb,tmbig%orbs%norb,nlocregPerMPI),intent(in):: ham!, ovrlp
      real(8),dimension(tmb%orbs%npsidim_orbs),intent(out):: lphi
    end subroutine build_input_guess


    subroutine extractMatrix3(iproc, nproc, norb, norbp, orbstot, onWhichAtomPhi, onWhichMPI, nmat, ham, matmin, hamextract)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, nmat, norb, norbp
      type(orbitals_data),intent(in):: orbstot
      integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
      real(8),dimension(orbstot%norb,orbstot%norb,nmat),intent(in):: ham
      type(matrixMinimization),intent(inout):: matmin
      real(8),dimension(:,:,:),pointer,intent(out):: hamextract
    end subroutine extractMatrix3

!!$    subroutine getHamiltonianMatrix3(iproc, nproc, nprocTemp, lzdig, orbsig, orbs, norb_parTemp, onWhichMPITemp, &
!!$               Glr, input, onWhichAtom, onWhichAtomp, nat, nlocregPerMPI, lchi, lhchi, ham)
!!$      use module_base
!!$      use module_types
!!$      implicit none
!!$      integer,intent(in):: iproc, nproc, nprocTemp, nat, nlocregPerMPI
!!$      type(local_zone_descriptors),intent(in):: lzdig
!!$      type(orbitals_data),intent(in):: orbsig, orbs
!!$      integer,dimension(0:nprocTemp),intent(in):: norb_parTemp
!!$      integer,dimension(orbs%norb),intent(in):: onWhichMPITemp
!!$      type(locreg_descriptors),intent(in):: Glr
!!$      type(input_variables),intent(in):: input
!!$      integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!$      integer,dimension(orbsig%norbp),intent(in):: onWhichAtomp
!!$      real(8),dimension(orbsig%npsidim_comp),intent(in):: lchi
!!$      real(8),dimension(orbsig%npsidim_comp,nat),intent(in):: lhchi
!!$      real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham
!!$      end subroutine getHamiltonianMatrix3
!!$
!!$      subroutine getHamiltonianMatrix4(iproc, nproc, nprocTemp, lzdig, orbsig, orbs, norb_parTemp, onWhichMPITemp, &
!!$                 Glr, input, onWhichAtom, onWhichAtomp, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad, &
!!$                 memoryForCommunOverlapIG, tag, ham)
!!$        use module_base
!!$        use module_types
!!$        implicit none
!!$        integer,intent(in):: iproc, nproc, nprocTemp, ndim_lhchi, nlocregPerMPI
!!$        type(local_zone_descriptors),intent(in):: lzdig
!!$        type(orbitals_data),intent(in):: orbsig, orbs
!!$        integer,dimension(0:nprocTemp),intent(in):: norb_parTemp
!!$        integer,dimension(orbs%norb),intent(in):: onWhichMPITemp
!!$        type(locreg_descriptors),intent(in):: Glr
!!$        type(input_variables),intent(in):: input
!!$        integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!$        integer,dimension(orbsig%norbp),intent(in):: onWhichAtomp
!!$        real(8),dimension(orbsig%npsidim),intent(in):: lchi
!!$        real(8),dimension(orbsig%npsidim,ndim_lhchi),intent(in):: lhchi
!!$        logical,dimension(lzdig%nlr),intent(in):: skip
!!$        type(matrixDescriptors),intent(in):: mad
!!$        integer,intent(in):: memoryForCommunOverlapIG
!!$        integer,intent(inout):: tag
!!$        !logical,dimension(lin%lig%lzdig%nlr,0:nproc-1),intent(in):: skipGlobal
!!$        real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham
!!$      end subroutine getHamiltonianMatrix4
!!$
!!$
!!$      subroutine getHamiltonianMatrix5(iproc, nproc, nprocTemp, lzdig, orbsig, orbs, norb_parTemp, onWhichMPITemp, &
!!$                 Glr, input, onWhichAtom, onWhichAtomp, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad, &
!!$                 memoryForCommunOverlapIG, tag, ham)
!!$        use module_base
!!$        use module_types
!!$        implicit none
!!$        integer,intent(in):: iproc, nproc, nprocTemp, ndim_lhchi, nlocregPerMPI
!!$        type(local_zone_descriptors),intent(in):: lzdig
!!$        type(orbitals_data),intent(in):: orbsig, orbs
!!$        integer,dimension(0:nprocTemp),intent(in):: norb_parTemp
!!$        integer,dimension(orbs%norb),intent(in):: onWhichMPITemp
!!$        type(locreg_descriptors),intent(in):: Glr
!!$        type(input_variables),intent(in):: input
!!$        integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!$        integer,dimension(orbsig%norbp),intent(in):: onWhichAtomp
!!$        real(8),dimension(orbsig%npsidim),intent(in):: lchi
!!$        real(8),dimension(orbsig%npsidim,ndim_lhchi),intent(in):: lhchi
!!$        logical,dimension(lzdig%nlr),intent(in):: skip
!!$        type(matrixDescriptors),intent(in):: mad
!!$        integer,intent(in):: memoryForCommunOverlapIG
!!$        integer,intent(inout):: tag
!!$        !logical,dimension(lin%lig%lzdig%nlr,0:nproc-1),intent(in):: skipGlobal
!!$        real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham
!!$      end subroutine getHamiltonianMatrix5

      subroutine allocateSendBufferOrtho(comon, subname)
        use module_base
        use module_types
        implicit none
        type(p2pComms),intent(inout):: comon
        character(len=*),intent(in):: subname
      end subroutine allocateSendBufferOrtho
      
      
      subroutine deallocateSendBufferOrtho(comon, subname)
        use module_base
        use module_types
        implicit none
        type(p2pComms),intent(inout):: comon
        character(len=*),intent(in):: subname
      end subroutine deallocateSendBufferOrtho
      
      
      subroutine allocateRecvBufferOrtho(comon, subname)
        use module_base
        use module_types
        implicit none
        type(p2pComms),intent(inout):: comon
        character(len=*),intent(in):: subname
      end subroutine allocateRecvBufferOrtho
      
      
      subroutine deallocateRecvBufferOrtho(comon, subname)
        use module_base
        use module_types
        implicit none
        type(p2pComms),intent(inout):: comon
        character(len=*),intent(in):: subname
      end subroutine deallocateRecvBufferOrtho

      subroutine applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, blocksize_pdgemm, &
                 correction_orthoconstraint, &
                 orbs, lagmat, ovrlp, mad, &
                 ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm, correction_orthoconstraint
        type(orbitals_data),intent(in):: orbs
        real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
        real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
      end subroutine applyOrthoconstraintNonorthogonal2

      !!subroutine gatherOrbitalsOverlapWithComput(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp, expanded)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc
      !!  type(orbitals_data),intent(in):: orbs
      !!  type(input_variables),intent(in):: input
      !!  type(local_zone_descriptors),intent(in):: lzd
      !!  type(overlapParameters),intent(in):: op
      !!  type(p2pComms),intent(inout):: comon
      !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      !!  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
      !!end subroutine gatherOrbitalsOverlapWithComput


      !!subroutine expandOneOrbital(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc, orbsource, orbdest
      !!  type(orbitals_data),intent(in):: orbs
      !!  type(input_variables),intent(in):: input
      !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
      !!  type(local_zone_descriptors),intent(in):: lzd
      !!  type(overlapParameters),intent(in):: op
      !!  type(p2pComms),intent(in):: comon
      !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      !!end subroutine expandOneOrbital

      !!subroutine expandRemainingOrbitals(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, expanded, lphiovrlp)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc
      !!  type(orbitals_data),intent(in):: orbs
      !!  type(input_variables),intent(in):: input
      !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
      !!  type(local_zone_descriptors),intent(in):: lzd
      !!  type(overlapParameters),intent(in):: op
      !!  type(p2pComms),intent(in):: comon
      !!  logical,dimension(orbs%norb,orbs%norbp),intent(in):: expanded
      !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      !!end subroutine expandRemainingOrbitals

      subroutine extractOrbital3(iproc, nproc, orbs, orbsig, sizePhi, lzd, lzdig, op, opig, phi, nsendBuf, sendBuf)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, sizePhi
        type(orbitals_data),intent(in):: orbs, orbsig
        type(local_zone_descriptors),intent(in):: lzd, lzdig
        type(overlapParameters),intent(inout):: op, opig
        real(8),dimension(sizePhi),intent(in):: phi
        integer,intent(in):: nsendBuf
        real(8),dimension(nsendBuf),intent(out):: sendBuf
      end subroutine extractOrbital3

      subroutine calculateOverlapMatrix3(iproc, nproc, orbs, op, nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        real(8),dimension(nsendBuf),intent(in):: sendBuf
        real(8),dimension(nrecvBuf),intent(in):: recvBuf
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
      end subroutine calculateOverlapMatrix3


      subroutine calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, &
                 nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        real(8),dimension(nsendBuf),intent(in):: sendBuf
        real(8),dimension(nrecvBuf),intent(in):: recvBuf
        type(matrixDescriptors),intent(in):: mad
        !logical,dimension(0:nproc-1),intent(in):: skip
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
      end subroutine calculateOverlapMatrix3Partial

      subroutine dgemm_parallel(iproc, nproc, blocksize, comm, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        use module_base
        implicit none
        integer,intent(in):: iproc, nproc, blocksize, comm, m, n, k, lda, ldb, ldc
        character(len=1),intent(in):: transa, transb
        real(8),intent(in):: alpha, beta
        real(8),dimension(lda,k),intent(in):: a
        real(8),dimension(ldb,n),intent(in):: b
        real(8),dimension(ldc,n),intent(out):: c
      end subroutine dgemm_parallel

      subroutine dsymm_parallel(iproc, nproc, blocksize, comm, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        use module_base
        implicit none
        integer,intent(in):: iproc, nproc, blocksize, comm, m, n, lda, ldb, ldc
        character(len=1),intent(in):: side, uplo
        real(8),intent(in):: alpha, beta
        real(8),dimension(lda,m),intent(in):: a
        real(8),dimension(ldb,n),intent(in):: b
        real(8),dimension(ldc,n),intent(out):: c
      end subroutine dsymm_parallel



      subroutine applyOrthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOverlap, blocksize_pdgemm, &
                 comm, norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, ovrlp, &
                 lagmat, comom, mlr, mad, orbs, grad, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOverlap, blocksize_pdgemm
        integer,intent(in):: comm, norb, norbmax, norbp, isorb, nlr, noverlaps
        integer,dimension(norb),intent(in):: onWhichAtom
        real(8),dimension(norb,norb),intent(in):: ovrlp
        real(8),dimension(norb,norb),intent(inout):: lagmat
        type(p2pComms),intent(in):: comom
        type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
        type(matrixDescriptors),intent(in):: mad
        type(orbitals_data),intent(in):: orbs
        real(8),dimension(norbmax,norbp),intent(inout):: grad
        real(8),dimension(norb,norb),intent(out):: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
      end subroutine applyOrthoconstraintVectors


      subroutine dsyev_parallel(iproc, nproc, blocksize, comm, jobz, uplo, n, a, lda, w, info)
        use module_base
        use module_types
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, nproc, blocksize, comm, n, lda, info
        character(len=1),intent(in):: jobz, uplo
        real(8),dimension(lda,n),intent(inout):: a
        real(8),dimension(n),intent(out):: w
      end subroutine dsyev_parallel


      subroutine transformOverlapMatrixParallel(iproc, nproc, norb, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, norb
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine transformOverlapMatrixParallel


      !!subroutine initMatrixCompression(iproc, nproc, orbs, op, mad)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc
      !!  type(orbitals_data),intent(in):: orbs
      !!  type(overlapParameters),intent(in):: op
      !!  type(matrixDescriptors),intent(out):: mad
      !!end subroutine initMatrixCompression

      subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, orbs, op, comon, mad, collcom, orthpar, bpo, bs, &
                 lphi, lhphi, lagmat, ovrlp, psit_c, psit_f, hpsit_c, hpsit_f, can_use_transposed, overlap_calculated)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(local_zone_descriptors),intent(in):: lzd
        type(orbitals_Data),intent(in):: orbs
        type(overlapParameters),intent(inout):: op
        type(p2pComms),intent(inout):: comon
        type(matrixDescriptors),intent(in):: mad
        type(collective_comms),intent(in):: collcom
        type(orthon_data),intent(in):: orthpar
        type(basis_performance_options),intent(in):: bpo
        type(basis_specifications),intent(in):: bs
        real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lphi
        real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lhphi
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: lagmat, ovrlp
        real(8),dimension(:),pointer:: psit_c, psit_f, hpsit_c, hpsit_f
        logical,intent(inout):: can_use_transposed, overlap_calculated
      end subroutine orthoconstraintNonorthogonal

      subroutine dsygv_parallel(iproc, nproc, blocksize, nprocMax, comm, itype, jobz, uplo, n, a, lda, b, ldb, w, info)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, blocksize, nprocMax, comm, itype, n, lda, ldb, info
        character(len=1),intent(in):: jobz, uplo
        real(8),dimension(lda,n),intent(inout):: a
        real(8),dimension(ldb,n),intent(inout):: b
        real(8),dimension(n),intent(out):: w
      end subroutine dsygv_parallel

      !!subroutine getOrbitals(iproc, nproc, comon)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc
      !!  type(p2pComms),intent(inout):: comon
      !!end subroutine getOrbitals

      subroutine initCompressedMatmul(iproc, nproc, norb, mad)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, norb
        type(matrixDescriptors),intent(inout):: mad
      end subroutine initCompressedMatmul

      !!subroutine initCompressedMatmul2(norb, nseg, keyg, nsegmatmul, keygmatmul, keyvmatmul)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: norb, nseg
      !!  integer,dimension(2,nseg),intent(in):: keyg
      !!  integer,intent(out):: nsegmatmul
      !!  integer,dimension(:,:),pointer,intent(out):: keygmatmul
      !!  integer,dimension(:),pointer,intent(out):: keyvmatmul
      !!end subroutine initCompressedMatmul2


      subroutine dgemm_compressed2(iproc, nproc, norb, nsegline, nseglinemax, keygline, nsegmatmul, keygmatmul, a, b, c)
        implicit none
        integer,intent(in):: iproc, nproc, norb, nseglinemax, nsegmatmul
        integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
        integer,dimension(norb):: nsegline
        integer,dimension(2,nseglinemax,norb):: keygline
        real(8),dimension(norb,norb),intent(in):: a, b
        real(8),dimension(norb,norb),intent(out):: c
      end subroutine dgemm_compressed2

      subroutine transformOverlapMatrixTaylorOrder2(iproc, nproc, norb, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, norb
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine transformOverlapMatrixTaylorOrder2

      subroutine overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
                 blocksize_pdgemm, norb, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine overlapPowerMinusOneHalf

      subroutine overlapPowerMinusOne(iproc, nproc, iorder, blocksize, norb, mad, orbs, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, iorder, blocksize, norb
        type(matrixDescriptors),intent(in):: mad
        type(orbitals_data),intent(in):: orbs
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine overlapPowerMinusOne

      subroutine initCompressedMatmul3(iproc, norb, mad)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, norb
        type(matrixDescriptors),intent(inout):: mad
      end subroutine initCompressedMatmul3

    subroutine Linearnonlocal_forces(iproc,nproc,Lzd,nlpspd,hx,hy,hz,at,rxyz,&
      orbs,proj,psi,fsep,refill,linorbs,coeff,phi)
      use module_base
      use module_types
      implicit none
      type(atoms_data), intent(in) :: at
      logical, intent(in) :: refill
      integer, intent(in) :: iproc, nproc
      real(gp), intent(in) :: hx,hy,hz
      type(local_zone_descriptors) :: Lzd
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(orbitals_data), intent(in) :: orbs
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: psi
      real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
      real(gp), dimension(3,at%nat), intent(inout) :: fsep
      type(orbitals_data), intent(in) :: linorbs                         
      real(8),dimension(linorbs%npsidim_orbs),intent(in),optional:: phi          
      real(8),dimension(linorbs%norb,orbs%norb),intent(in),optional:: coeff  
    end subroutine Linearnonlocal_forces

     !!!subroutine local_hamiltonian3(iproc,exctX,orbs,Lzd,hx,hy,hz,&
     !!!     nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,&
     !!!     withConfinement, at, rxyz, istexct, lin, confinementCenter)
     !!!  use module_base
     !!!  use module_types
     !!!  use libxc_functionals
     !!!  implicit none
     !!!  integer, intent(in) :: iproc,nspin, istexct
     !!!  real(gp), intent(in) :: hx,hy,hz
     !!!  logical, intent(in) :: exctX
     !!!  type(orbitals_data), intent(in) :: orbs
     !!!  type(local_zone_descriptors), intent(in) :: Lzd
     !!!  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
     !!!  real(wp), dimension(Lzd%ndimpotisf),target :: Lpot
     !!!  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
     !!!  real(gp), intent(out) :: ekin_sum,epot_sum
     !!!  real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
     !!!  logical,intent(in):: withConfinement
     !!!  type(atoms_data), intent(in) :: at
     !!!  real(gp), dimension(3,at%nat), intent(in) :: rxyz
     !!!  type(linearParameters),intent(in),optional:: lin
     !!!  integer,dimension(orbs%norbp),intent(in),optional:: confinementCenter
     !!!end subroutine local_hamiltonian3


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


     subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,rxyz,&
          proj,Lzd,nlpspd,confdatarr,ngatherarr,Lpot,psi,hpsi,&
          energs,SIC,GPU,&
          pkernel,orbsocc,psirocc)
       use module_base
       use module_types
       use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc!,nspin
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
       !real(wp), dimension(lzd%ndimpotisf) :: Lpot
       real(wp), dimension(:),pointer :: Lpot
       type(energy_terms), intent(inout) :: energs
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       type(coulomb_operator), intent(in), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
     end subroutine FullHamiltonianApplication

     !!subroutine prepare_lnlpspd(iproc, at, input, orbs, rxyz, radii_cf, locregShape, lzd)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc
     !!  type(atoms_data),intent(in):: at
     !!  type(input_variables),intent(in):: input
     !!  type(orbitals_data),intent(in):: orbs
     !!  real(8),dimension(3,at%nat),intent(in):: rxyz
     !!  real(8),dimension(at%ntypes,3),intent(in):: radii_cf
     !!  character(len=1),intent(in):: locregShape
     !!  type(local_zone_descriptors),intent(inout):: lzd
     !!end subroutine prepare_lnlpspd

     !!subroutine free_lnlpspd(orbs, lzd)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  type(orbitals_data),intent(in):: orbs
     !!  type(local_zone_descriptors),intent(inout):: lzd
     !!end subroutine free_lnlpspd


     subroutine transformToGlobal(iproc,nproc,lzd,lorbs,orbs,comms,input,ld_coeff,coeff,lphi,psi,psit)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, ld_coeff
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: lorbs, orbs
       type(communications_arrays):: comms
       type(input_variables),intent(in):: input
       real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
       real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: lphi
       real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out):: psi
       real(8),dimension(:),pointer,intent(inout):: psit
     end subroutine transformToGlobal


     !!subroutine my_iallgatherv(iproc, nproc, sendbuf, sendcount, recvbuf, recvcounts, displs, comm, tagx, requests)
     !!  use module_base
     !!  implicit none

     !!  ! Calling arguments
     !!  integer,intent(in):: iproc, nproc, sendcount, comm
     !!  integer,dimension(0:nproc-1),intent(in):: recvcounts, displs
     !!  real(8),dimension(sendcount),intent(in):: sendbuf
     !!  integer,dimension(2,0:nproc*nproc-1),intent(in):: requests
     !!  integer,intent(in):: tagx
     !!  real(8),dimension(sum(recvcounts)),intent(out):: recvbuf
     !!end subroutine my_iallgatherv


     !!subroutine gatherOrbitalsOverlapWithComput2(iproc, nproc, orbs, input, lzd, op, comon, nsendbuf, sendbuf,&
     !!     nrecvbuf, recvbuf, lphiovrlp, expanded, ovrlp)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf                                                                 
     !!  type(orbitals_data),intent(in):: orbs
     !!  type(input_variables),intent(in):: input
     !!  type(local_zone_descriptors),intent(in):: lzd                                                                         
     !!  type(overlapParameters),intent(in):: op                                                                               
     !!  type(p2pComms),intent(inout):: comon
     !!  real(8),dimension(nsendbuf),intent(in):: sendbuf                                                                      
     !!  real(8),dimension(nrecvbuf),intent(in):: recvbuf                                                                      
     !!  real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
     !!  logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
     !!  real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
     !!end subroutine gatherOrbitalsOverlapWithComput2


     subroutine getStartingIndices(iorb, jorb, op, orbs, ist, jst)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iorb, jorb
       type(overlapParameters),intent(in):: op
       type(orbitals_data),intent(in):: orbs
       integer,intent(out):: ist, jst
     end subroutine getStartingIndices


      subroutine getStartingIndicesGlobal(iiorbx, jjorbx, op, orbs, ist, jst, ncount)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iiorbx, jjorbx
        type(overlapParameters),intent(in):: op
        type(orbitals_data),intent(in):: orbs
        integer,intent(out):: ist, jst, ncount
      end subroutine getStartingIndicesGlobal

       subroutine my_iallgather_collect2(iproc, nproc, sendcount, recvcounts, requests)
         use module_base
         implicit none
         integer,intent(in):: iproc, nproc, sendcount
         integer,dimension(0:nproc-1),intent(in):: recvcounts
         integer,dimension(2,0:nproc-1),intent(inout):: requests
       end subroutine my_iallgather_collect2


       subroutine initMatrixCompression(iproc, nproc, nlr, ndim, orbs, noverlaps, overlaps, mad)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr, ndim
         type(orbitals_data),intent(in):: orbs
         integer,dimension(orbs%norb),intent(in):: noverlaps
         integer,dimension(ndim,orbs%norb),intent(in):: overlaps
         type(matrixDescriptors),intent(out):: mad
       end subroutine initMatrixCompression


       !!subroutine postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer,intent(in):: iproc, nproc, newComm
       !!  type(p2pComms),intent(inout):: comom
       !!end subroutine postCommsVectorOrthonormalizationNew


       !!subroutine gatherVectorsNew(iproc, nproc, comom)
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer,intent(in):: iproc, nproc
       !!  type(p2pComms),intent(inout):: comom
       !!end subroutine gatherVectorsNew


       subroutine compressMatrixPerProcess(iproc, nproc, orbs, mad, mat, size_lmat, lmat)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, size_lmat
         type(orbitals_data),intent(in):: orbs
         type(matrixDescriptors),intent(in):: mad
         real(8),dimension(orbs%norb**2),intent(in):: mat
         real(8),dimension(size_lmat),intent(out):: lmat
       end subroutine compressMatrixPerProcess


       subroutine getCommunArraysMatrixCompression(iproc, nproc, orbs, mad, sendcounts, displs)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(matrixDescriptors),intent(in):: mad
         integer,dimension(0:nproc-1),intent(out):: sendcounts, displs
       end subroutine getCommunArraysMatrixCompression


       subroutine get_hamiltonian_matrices(iproc, nproc, lzd, lzdig, orbsig, orbs, &
                  input, hx, hy, hz, onWhichAtom, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad,&
                  memoryForCommunOverlapIG, locregShape, bpo, ham)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, ndim_lhchi, nlocregPerMPI
         real(gp), intent(in) :: hx, hy, hz
         type(local_zone_descriptors),intent(in):: lzd, lzdig
         type(orbitals_data),intent(in):: orbsig, orbs
         type(input_variables),intent(in):: input
         integer,dimension(orbsig%norb),intent(in):: onWhichAtom
         real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
         real(8),dimension(orbsig%npsidim_orbs,ndim_lhchi),intent(in):: lhchi
         logical,dimension(lzd%nlr),intent(in):: skip
         type(matrixDescriptors),intent(in):: mad
         integer,intent(in):: memoryForCommunOverlapIG
         character(len=1),intent(in):: locregShape
         type(basis_performance_options),intent(inout):: bpo
         real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham
       end subroutine get_hamiltonian_matrices
       
       subroutine dgemm_compressed_parallel(iproc, nproc, norb, nsegline, nseglinemax, keygline, &
                  nsegmatmul, keygmatmul, norb_par, isorb_par, norbp, a, b, c)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norb, norbp, nseglinemax, nsegmatmul
         integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
         integer,dimension(norb):: nsegline
         integer,dimension(2,nseglinemax,norb):: keygline
         integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
         real(8),dimension(norb,norb),intent(in):: a, b
         real(8),dimension(norb,norb),intent(out):: c
       end subroutine dgemm_compressed_parallel


       !!!subroutine getCoefficients_new(iproc, nproc, lin, orbs, hamold, lphi, ovrlp, coeff)
       !!!  use module_base
       !!!  use module_types
       !!!  implicit none
       !!!  integer,intent(in):: iproc, nproc
       !!!  type(linearParameters),intent(inout):: lin
       !!!  type(orbitals_data),intent(in):: orbs
       !!!  real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: hamold
       !!!  real(8),dimension(lin%orbs%npsidim_orbs),intent(in):: lphi
       !!!  real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(inout):: ovrlp
       !!!  real(8),dimension(lin%orbs%norb,orbs%norb),intent(inout):: coeff
       !!!end subroutine getCoefficients_new


       subroutine apply_confinement(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
            rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, &
            ibyyzz_r) !optional
         use module_base
         implicit none
         integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor, confPotOrder, offsetx, offsety, offsetz
         real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
         integer,dimension(2,-14:2*n2+16,-14:2*n3+16),intent(in),optional :: ibyyzz_r
         real(8),dimension(3),intent(in):: rxyzConfinement
         real(8),intent(in):: hxh,hyh,hzh,potentialPrefac
       end subroutine apply_confinement



       subroutine unitary_optimization(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
        newgradient, confdatarr, hx, lphi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nit
         type(local_zone_descriptors),intent(in):: lzd
         type(orbitals_data),intent(in):: orbs
         type(atoms_data),intent(in):: at
         type(overlapParameters),intent(inout):: op
         type(p2pComms),intent(inout):: comon
         type(matrixDescriptors),intent(in):: mad
         real(8),dimension(3,at%nat),intent(in):: rxyz
         real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
         logical,intent(in):: newgradient
         real(8),intent(in):: hx
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),dimension(orbs%npsidim_comp),intent(inout):: lphi
       end subroutine unitary_optimization


      subroutine build_new_linear_combinations(iproc, nproc, lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(local_zone_descriptors),intent(in):: lzd
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        integer,intent(in):: nrecvbuf
        real(8),dimension(nrecvbuf),intent(in):: recvbuf
        real(8),dimension(orbs%norb,orbs%norb),intent(in):: omat
        logical,intent(in):: reset
        real(8),dimension(orbs%npsidim_comp),intent(out):: lphi
      end subroutine build_new_linear_combinations


      !!subroutine indicesForExpansion(iproc, nproc, nspin, orbs, onWhichAtom, lzd, op, comon)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc, nspin
      !!  type(orbitals_data),intent(in):: orbs
      !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
      !!  type(local_zone_descriptors),intent(in):: lzd
      !!  type(overlapParameters),intent(inout):: op
      !!  type(p2pComms),intent(in):: comon
      !!end subroutine indicesForExpansion


      !!subroutine nullify_expansionSegments(expseg)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  type(expansionSegments),intent(out):: expseg
      !!end subroutine nullify_expansionSegments


      !!subroutine indicesForExtraction(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, comon)
      !!  use module_base
      !!  use module_types
      !!  implicit none
      !!  integer,intent(in):: iproc, nproc, sizePhi
      !!  type(orbitals_data),intent(in):: orbs
      !!  integer,dimension(orbs%norb),intent(in):: onWhichAtom
      !!  type(local_zone_descriptors),intent(in):: lzd
      !!  type(overlapParameters),intent(inout):: op
      !!  type(p2pComms),intent(out):: comon
      !!end subroutine indicesForExtraction


      subroutine allocate_workarrays_quartic_convolutions(lr, subname, work)
        use module_base
        use module_types
        implicit none
        type(locreg_descriptors),intent(in):: lr
        character(len=*),intent(in):: subname
        type(workarrays_quartic_convolutions),intent(out):: work
      end subroutine allocate_workarrays_quartic_convolutions


      subroutine deallocate_workarrays_quartic_convolutions(lr, subname, work)
        use module_base
        use module_types
        implicit none
        type(locreg_descriptors),intent(in):: lr
        character(len=*),intent(in):: subname
        type(workarrays_quartic_convolutions),intent(out):: work
      end subroutine deallocate_workarrays_quartic_convolutions

      subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
                 hgrid, offsetx, offsety, offsetz, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
                 rxyzConf, potentialPrefac, with_kinetic, cprecr, maxdim, &
                 xx_c, xx_f1, xx_f, xy_c, xy_f2, xy_f,  xz_c, xz_f4, xz_f, &
                 aeff0array, beff0array, ceff0array, eeff0array, &
                 aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array, &
                 aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray, &
                 xya_c, xyb_c, xyc_c, xye_c, xza_c, xzb_c, xzc_c, xze_c, &
                 yza_c, yzb_c, yzc_c, yze_c, xya_f, xyb_f, xyc_f, xye_f, &
                 xza_f, xzb_f, xzc_f, xze_f, yza_f, yzb_f, yzc_f, yze_f, &
                 aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, &
                 ceff0, ceff1, ceff2, ceff3, eeff0, eeff1, eeff2, eeff3, &
                 aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2, &
                 ceff0_2, ceff1_2, ceff2_2, ceff3_2, eeff0_2, eeff1_2, eeff2_2, eeff3_2, & 
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
        real(wp),dimension(0:n3,0:n1,0:n2):: xza_c, xzb_c, xzc_c, xze_c, yza_c, yzb_c, yzc_c, yze_c
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
        real(wp),dimension(35):: aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, ceff0, ceff1, ceff2, ceff3
        real(wp),dimension(29):: eeff0, eeff1, eeff2, eeff3
        real(wp),dimension(35):: aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2
        real(wp),dimension(35):: ceff0_2, ceff1_2, ceff2_2, ceff3_2
        real(wp),dimension(29):: eeff0_2, eeff1_2, eeff2_2, eeff3_2
        real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
        real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
      end subroutine ConvolQuartic4



       subroutine apply_potential_lr(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,&
            psir,pot,epot,&
            confdata,ibyyzz_r) !optional
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
         integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
         real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
         real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
         type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
         real(gp), intent(out) :: epot
       end subroutine apply_potential_lr

       subroutine psir_to_vpsi(npot,nspinor,lr,pot,vpsir,epot,confdata)
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
       end subroutine psir_to_vpsi

       subroutine erf_stress(at,rxyz,hxh,hyh,hzh,n1i,n2i,n3i,n3p,iproc,nproc,ngatherarr,rho,tens)
         use module_base
         use module_types
         implicit none
         !passed var
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
         real(gp), intent(in) :: hxh,hyh,hzh
         integer,intent(in) :: n1i,n2i,n3i,n3p,iproc,nproc
         real(kind=8), dimension(n1i*n2i*max(n3p,1)), intent(in), target :: rho
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
         real(dp),dimension(6), intent(out) :: tens
       end subroutine erf_stress

       subroutine AtomicOrbitals_forLinear(iproc,at,rxyz,mapping,norbe,orbse,norbsc,&
            &   nspin,eks,scorb,G,gaucoeff,iorbtolr)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: norbe,iproc
         integer, intent(in) :: norbsc,nspin
         type(atoms_data), intent(in) :: at
         logical, dimension(4,2,at%natsc), intent(in) :: scorb
         real(gp), dimension(3,at%nat), intent(in), target :: rxyz
         type(orbitals_data), intent(inout) :: orbse
         integer,dimension(orbse%norb),intent(in):: mapping
         type(gaussian_basis), intent(out) :: G
         real(gp), intent(out) :: eks
         integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
         !real(wp), dimension(norbe,orbse%nspinor,orbse%norbp), intent(out) :: gaucoeff !norbe=G%ncoeff
         real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff
       end subroutine AtomicOrbitals_forLinear


       subroutine apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
                  confdatarr, hx, psi, centralLocreg, vpsi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, centralLocreg
         type(atoms_data),intent(in):: at
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         real(8),dimension(3,at%nat),intent(in):: rxyz
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),intent(in):: hx
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi
       end subroutine apply_orbitaldependent_potential


       subroutine get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
                  confdatarr, hx, psi, potmat)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(atoms_data),intent(in):: at
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         type(overlapParameters),intent(inout):: op
         type(p2pComms),intent(inout):: comon
         type(matrixDescriptors),intent(in):: mad
         real(8),dimension(3,at%nat),intent(in):: rxyz
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),intent(in):: hx
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
         real(8),dimension(orbs%norb,orbs%norb,at%nat),intent(out):: potmat
       end subroutine get_potential_matrices
       
       subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         integer, intent(in) :: linType
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(in) :: atoms
         type(orbitals_data),intent(inout) :: orbs
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
       end subroutine check_linear_and_create_Lzd

       subroutine create_LzdLIG(iproc,nproc,nspin,linearmode,hx,hy,hz,Glr,atoms,orbs,rxyz,Lzd)
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
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
       end subroutine create_LzdLIG

!!       subroutine reinitialize_Lzd_after_LIG(iproc,nproc,input,Lzd,atoms,orbs,rxyz)
!!         use module_base
!!         use module_types
!!         implicit none
!!         integer, intent(in) :: iproc,nproc
!!         type(input_variables), intent(in) :: input
!!         type(local_zone_descriptors), intent(inout) :: Lzd
!!         type(atoms_data), intent(in) :: atoms
!!         type(orbitals_data),intent(inout) :: orbs
!!         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!       end subroutine reinitialize_Lzd_after_LIG

       subroutine system_initialization(iproc,nproc,inputpsi,input_wf_format,in,atoms,rxyz,&
            orbs,lorbs,Lzd,Lzd_lin,denspot,nlpspd,comms,lcomms,shift,proj,radii_cf)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         integer, intent(out) :: inputpsi,input_wf_format
         type(input_variables), intent(in) :: in 
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
         type(orbitals_data), intent(out) :: orbs,lorbs
         type(local_zone_descriptors), intent(out) :: Lzd, Lzd_lin
         type(DFT_local_fields), intent(out) :: denspot
         type(nonlocal_psp_descriptors), intent(out) :: nlpspd
         type(communications_arrays), intent(out) :: comms,lcomms
         real(gp), dimension(3), intent(out) :: shift  !< shift on the initial positions
         real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
         real(wp), dimension(:), pointer :: proj
       end subroutine system_initialization

       subroutine nullify_p2pComms(p2pcomm)
         use module_base
         use module_types
         implicit none
         type(p2pComms),intent(inout):: p2pcomm
       end subroutine nullify_p2pComms

       subroutine extract_potential_for_spectra(iproc,nproc,at,rhod,dpbox,&
            orbs,nvirt,comms,Lzd,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
            nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
            nspin,potshortcut,symObj,GPU,input)
         use module_base
         use module_types
         implicit none
         !Arguments
         integer, intent(in) :: iproc,nproc,ixc
         integer, intent(inout) :: nspin,nvirt
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(inout) :: at
         type(rho_descriptors),intent(in) :: rhod
         type(denspot_distribution), intent(in) :: dpbox
         type(orbitals_data), intent(inout) :: orbs
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(communications_arrays), intent(in) :: comms
         type(GPU_pointers), intent(inout) :: GPU
         type(input_variables):: input
         type(symmetry_data), intent(in) :: symObj
         !integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
         !integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
         real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
         type(gaussian_basis), intent(out) :: G !basis for davidson IG
         real(wp), dimension(:), pointer :: psi,hpsi,psit
         real(wp), dimension(:,:,:,:), pointer :: rhocore
         type(coulomb_operator), intent(in) :: pkernel,pkernelseq
         integer, intent(in) ::potshortcut
       end subroutine extract_potential_for_spectra

       subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,iscf,alphamix,ixc,&
            nlpspd,proj,rxyz,linflag,unblock_comms,GPU,wfn,&
            energs,rpnrm,xcstr)
         use module_base
         use module_types
         use m_ab6_mixing
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
         real(gp), intent(out) :: rpnrm
         real(gp), dimension(6), intent(out) :: xcstr
       end subroutine psitohpsi

       subroutine assignToLocreg2(iproc, nproc, norb, norb_par, natom, nlr, nspin, Localnorb, rxyz, inwhichlocreg)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: nlr,iproc,nproc,nspin,natom,norb
         integer,dimension(nlr),intent(in):: Localnorb
         integer,dimension(0:nproc-1),intent(in):: norb_par
         real(8),dimension(3,nlr),intent(in):: rxyz
         integer,dimension(:),pointer,intent(out):: inwhichlocreg
       end subroutine assignToLocreg2
       
       subroutine calc_gradient(geocode,n1,n2,n3,n3grad,deltaleft,deltaright,rhoinp,nspden,hx,hy,hz,&
            gradient,rhocore)
         use module_base
         implicit none
         !Arguments
         character(len=1), intent(in) :: geocode
         integer, intent(in) :: n1,n2,n3,n3grad,deltaleft,deltaright,nspden
         real(dp), intent(in) :: hx,hy,hz
         real(dp), dimension(n1,n2,n3,nspden), intent(inout) :: rhoinp
         real(dp), dimension(n1,n2,n3grad,2*nspden-1,0:3), intent(out) :: gradient
         real(dp), dimension(:,:,:,:), pointer :: rhocore
       end subroutine calc_gradient

       !!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
       !!!       hxh, hyh, hzh, ioffset, dir, &
       !!!       ibyyzz_r) !optional
       !!!  use module_base
       !!!  implicit none
       !!!  integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
       !!!  integer,dimension(3),intent(in):: ioffset
       !!!  real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
       !!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
       !!!  real(8),intent(in):: hxh, hyh, hzh
       !!!  character(len=1),intent(in):: dir
       !!!end subroutine position_operator

       subroutine apply_position_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, xpsi, ypsi, zpsi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, order
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         real(8),intent(in):: hx, hy, hz
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: xpsi, ypsi, zpsi
       end subroutine apply_position_operators

       !!subroutine MLWF(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
       !!             newgradient, confdatarr, hx, locregCenters, lphi, Umat)
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer,intent(in):: iproc, nproc, nit
       !!  type(local_zone_descriptors),intent(in):: lzd
       !!  type(orbitals_data),intent(in):: orbs
       !!  type(atoms_data),intent(in):: at
       !!  type(overlapParameters),intent(inout):: op
       !!  type(p2pComms),intent(inout):: comon
       !!  type(matrixDescriptors),intent(in):: mad
       !!  real(8),dimension(3,at%nat),intent(in):: rxyz
       !!  real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
       !!  logical,intent(in):: newgradient
       !!  real(8),intent(in):: hx, maxDispl
       !!  type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
       !!  real(8),dimension(3,lzd%nlr),intent(in):: locregCenters
       !!  real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
       !!  real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
       !!end subroutine MLWF




       subroutine MLWFnew(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
                    confdatarr, hx, locregCenters, maxDispl, lphi, Umat, centers)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nit
         type(local_zone_descriptors),intent(in):: lzd
         type(orbitals_data),intent(in):: orbs
         type(atoms_data),intent(in):: at
         type(overlapParameters),intent(inout):: op
         type(p2pComms),intent(inout):: comon
         type(matrixDescriptors),intent(in):: mad
         real(8),dimension(3,at%nat),intent(in):: rxyz
         real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
         !logical,intent(in):: newgradient
         real(8),intent(in):: hx, maxDispl
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),dimension(3,lzd%nlr),intent(in):: locregCenters
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
         real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
         real(8),dimension(3,lzd%nlr),intent(out):: centers
       end subroutine MLWFnew


       !!!subroutine create_new_locregs(iproc, nproc, nlr, hx, hy, hz, lorbs, glr, locregCenter, &
       !!!           locrad, nscatterarr, withder, &
       !!!           inwhichlocreg_reference, ldiis, &
       !!!           lphilarge, lhphilarge, lhphilargeold, lphilargeold,tmb)
       !!!  use module_base
       !!!  use module_types
       !!!  implicit none
       !!!  integer,intent(in):: iproc, nproc, nlr
       !!!  real(8),intent(in):: hx, hy, hz
       !!!  type(orbitals_data),intent(in):: lorbs
       !!!  type(locreg_descriptors),intent(in):: glr
       !!!  real(8),dimension(3,nlr),intent(in):: locregCenter
       !!!  real(8),dimension(nlr):: locrad
       !!!  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       !!!  logical,intent(in):: withder
       !!!  integer,dimension(lorbs%norb),intent(in):: inwhichlocreg_reference
       !!!  type(localizedDIISParameters),intent(inout):: ldiis
       !!!  real(8),dimension(:),pointer,intent(out):: lphilarge, lhphilarge, lhphilargeold, lphilargeold
       !!!  type(DFT_wavefunction),intent(out):: tmb
       !!!end subroutine create_new_locregs

       subroutine destroy_new_locregs(iproc, nproc, tmb)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(DFT_wavefunction),intent(inout):: tmb
       end subroutine destroy_new_locregs

       subroutine get_cutoff_weight(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
            cutoff, weight_in, weight_out, &
            confdata,ibyyzz_r) !optional
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
         integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
         real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
         real(8),intent(in):: cutoff
         real(8),intent(out):: weight_in, weight_out
         type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
       end subroutine get_cutoff_weight


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
       end subroutine position_operators

       subroutine apply_r_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, vpsi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, order
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         real(8),intent(in):: hx, hy, hz
         type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
         real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi
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
       end subroutine r_operator

       !!subroutine apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, centers, vpsi)
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer,intent(in):: iproc, nproc
       !!  type(orbitals_data),intent(in):: orbs
       !!  type(local_zone_descriptors),intent(in):: lzd
       !!  real(8),intent(in):: hx, hy, hz
       !!  type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
       !!  real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
       !!  real(8),dimension(3,lzd%nlr),intent(in):: centers
       !!  real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi
       !!end subroutine apply_rminusmu_operator

       !!subroutine rminusmu_operator(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,mu,&
       !!     confdata,ibyyzz_r) !optional
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
       !!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
       !!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
       !!  real(8),dimension(3):: mu
       !!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
       !!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
       !!end subroutine rminusmu_operator

       subroutine define_confinement_data(confdatarr,orbs,rxyz,at,hx,hy,hz,&
                  confpotorder,potentialprefac,Lzd,confinementCenter)
         use module_base
         use module_types
         implicit none
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         integer,intent(in):: confpotorder
         real(gp),dimension(at%ntypes),intent(in):: potentialprefac
         type(local_zone_descriptors), intent(in) :: Lzd
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         integer, dimension(orbs%norb), intent(in) :: confinementCenter
         type(confpot_data), dimension(orbs%norbp), intent(out) :: confdatarr
       end subroutine define_confinement_data

       subroutine update_locreg(iproc, nproc, nlr, locrad, inwhichlocreg_reference, locregCenter, glr_tmp, &
                  bpo, useDerivativeBasisFunctions, nscatterarr, hx, hy, hz, &
                  orbs_tmp, lzd, llborbs, lbop, lbcomon, lbcomgp, comsr, lbmad, lbcollcom)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         logical,intent(in):: useDerivativeBasisFunctions
         integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
         real(8),intent(in):: hx, hy, hz
         real(8),dimension(nlr),intent(in):: locrad
         type(orbitals_data),intent(in):: orbs_tmp
         integer,dimension(orbs_tmp%norb),intent(in):: inwhichlocreg_reference
         real(8),dimension(3,nlr),intent(in):: locregCenter
         type(locreg_descriptors),intent(in):: glr_tmp
         type(basis_performance_options),intent(in):: bpo
         type(local_zone_descriptors),intent(inout):: lzd
         type(orbitals_data),intent(inout):: llborbs
         type(overlapParameters),intent(inout):: lbop
         type(p2pComms),intent(inout):: lbcomon
         type(p2pComms),intent(inout):: lbcomgp
         type(p2pComms),intent(inout):: comsr
         type(matrixDescriptors),intent(inout):: lbmad
         type(collective_comms),intent(inout):: lbcollcom
       end subroutine update_locreg


       subroutine communicate_basis_for_density(iproc, nproc, lzd, llborbs, lphi, comsr)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(local_zone_descriptors),intent(in):: lzd
         type(orbitals_data),intent(in):: llborbs
         real(8),dimension(llborbs%npsidim_orbs),intent(in):: lphi
         type(p2pComms),intent(inout):: comsr
       end subroutine communicate_basis_for_density

       subroutine create_wfn_metadata(mode, nphi, lnorb, llbnorb, norb, norbp, input, wfnmd)
         use module_base
         use module_types
         implicit none
         character(len=1),intent(in):: mode
         integer,intent(in):: nphi, lnorb, llbnorb, norb, norbp
         type(input_variables),intent(in):: input
         type(wfn_metadata),intent(out):: wfnmd
       end subroutine create_wfn_metadata

       subroutine destroy_wfn_metadata(wfnmd)
         use module_base
         use module_types
         !use deallocatePointers
         implicit none
         type(wfn_metadata),intent(inout):: wfnmd
       end subroutine destroy_wfn_metadata

       subroutine create_DFT_wavefunction(mode, nphi, lnorb, norb, norbp, input, wfn)
         use module_base
         use module_types
         implicit none
         character(len=1),intent(in):: mode
         integer,intent(in):: nphi, lnorb, norb, norbp
         type(input_variables),intent(in):: input
         type(DFT_wavefunction),intent(out):: wfn
       end subroutine create_DFT_wavefunction
       
       subroutine destroy_DFT_wavefunction(wfn)
         use module_base
         use module_types
         implicit none
         type(DFT_wavefunction),intent(inout):: wfn
       end subroutine destroy_DFT_wavefunction

       subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, at, glr, use_derivative_basis, rxyz, &
                  lorbs)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspinor
         type(input_variables),intent(in):: input
         type(atoms_data),intent(in):: at
         type(locreg_descriptors),intent(in):: glr
         logical,intent(in):: use_derivative_basis
         real(8),dimension(3,at%nat),intent(in):: rxyz
         type(orbitals_data),intent(out):: lorbs
       end subroutine init_orbitals_data_for_linear

       subroutine mix_main(iproc, nproc, mixHist, compare_outer_loop, input, glr, alpha_mix, &
                  denspot, mixdiis, rhopotold, rhopotold_out, pnrm, pnrm_out)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, mixHist
         logical,intent(in):: compare_outer_loop
         type(input_variables),intent(in):: input
         type(locreg_descriptors),intent(in):: glr
         real(8),intent(in):: alpha_mix
         type(DFT_local_fields),intent(inout):: denspot
         type(mixrhopotDIISParameters),intent(inout):: mixdiis
         real(8),dimension(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout):: rhopotold, rhopotold_out
         real(8),intent(out):: pnrm, pnrm_out
       end subroutine mix_main

       subroutine redefine_locregs_quantities(iproc, nproc, hx, hy, hz, locrad, transform, lzd, tmb, denspot, &
                  ldiis)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         real(8),intent(in):: hx, hy, hz
         type(local_zone_descriptors),intent(inout):: lzd
         real(8),dimension(lzd%nlr),intent(in):: locrad
         logical,intent(in):: transform
         type(DFT_wavefunction),intent(inout):: tmb
         type(DFT_local_fields),intent(inout):: denspot
         type(localizedDIISParameters),intent(inout),optional:: ldiis
       end subroutine redefine_locregs_quantities

       !!!subroutine enlarge_locreg(iproc, nproc, hx, hy, hz, withder, lzd, locrad, &
       !!!           ldiis, denspot, tmb)
       !!!  use module_base
       !!!  use module_types
       !!!  implicit none
       !!!  integer,intent(in):: iproc, nproc
       !!!  real(8),intent(in):: hx, hy, hz
       !!!  logical,intent(in):: withder
       !!!  type(local_zone_descriptors),intent(inout):: lzd
       !!!  real(8),dimension(lzd%nlr),intent(in):: locrad
       !!!  type(localizedDIISParameters),intent(inout):: ldiis
       !!!  type(DFT_local_fields),intent(inout):: denspot
       !!!  type(DFT_wavefunction),intent(inout):: tmb
       !!!end subroutine enlarge_locreg


       subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, kernel, &
                  ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, &
                  fnrmMax, alpha_mean, alpha_max, energy_increased, tmb, lhphi, lhphiold, &
                  tmblarge, lhphilarge2, overlap_calculated, ovrlp, energs, hpsit_c, hpsit_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, it
         type(DFT_wavefunction),target,intent(inout):: tmblarge, tmb
         real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in):: kernel
         type(localizedDIISParameters),intent(inout):: ldiis
         real(8),dimension(tmb%orbs%norb),intent(inout):: fnrmOldArr
         real(8),dimension(tmb%orbs%norbp),intent(inout):: alpha
         real(8),intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
         real(8),intent(inout):: trHold
         logical,intent(out):: energy_increased
         real(8),dimension(:),target,intent(inout):: lhphilarge2
         real(8),dimension(:),target,intent(inout):: lhphi, lhphiold
         logical,intent(inout):: overlap_calculated
         real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout):: ovrlp
         type(energy_terms),intent(in) :: energs
         real(8),dimension(:),pointer:: hpsit_c, hpsit_f
       end subroutine calculate_energy_and_gradient_linear


       subroutine copy_basis_specifications(bsin, bsout, subname)
         use module_base
         use module_types
         implicit none
         type(basis_specifications),intent(in):: bsin
         type(basis_specifications),intent(out):: bsout
         character(len=*),intent(in):: subname
       end subroutine copy_basis_specifications

       subroutine copy_orthon_data(odin, odout, subname)
         use module_base
         use module_types
         implicit none
         type(orthon_data),intent(in):: odin
         type(orthon_data),intent(out):: odout
         character(len=*),intent(in):: subname
       end subroutine copy_orthon_data

       subroutine improveOrbitals(iproc, nproc, it, tmb, ldiis, lhphi, alpha)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, it
         type(DFT_wavefunction),intent(inout):: tmb
         type(localizedDIISParameters),intent(inout):: ldiis
         real(8),dimension(tmb%wfnmd%nphi),intent(in):: lhphi
         real(8),dimension(tmb%orbs%norbp),intent(in):: alpha
       end subroutine improveOrbitals

       subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
                  lhphi, lphiold, alpha, &
                  trH, meanAlpha, alpha_max, alphaDIIS)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, it
        type(localizedDIISParameters),intent(inout):: ldiis
        type(DFT_wavefunction),target,intent(inout):: tmb
        real(8),dimension(tmb%orbs%npsidim_orbs),intent(inout):: lhphi, lphiold
        real(8),intent(in):: trH, meanAlpha, alpha_max
        real(8),dimension(tmb%orbs%norbp),intent(out):: alpha, alphaDIIS
       end subroutine hpsitopsi_linear
       
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
       end subroutine DIISorSD
 
       subroutine psi_to_vlocpsi(iproc,orbs,Lzd,&
            ipotmethod,confdatarr,pot,psi,vpsi,pkernel,ixc,alphaSIC,epot_sum,evSIC)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,ipotmethod,ixc
         real(gp), intent(in) :: alphaSIC
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi !this dimension will be modified
         real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
         real(gp), intent(out) :: epot_sum,evSIC
         real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: vpsi
         type(coulomb_operator), intent(in) ::  pkernel !< the PSolver kernel which should be associated for the SIC schemes
       end subroutine psi_to_vlocpsi

       subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
                  input, tmb, denspot, ldiis, lscv)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         real(8),intent(in):: hx, hy, hz
         type(input_variables),intent(in):: input
         type(DFT_wavefunction),intent(inout):: tmb
         type(DFT_local_fields),intent(inout) :: denspot
         type(localizedDIISParameters),intent(inout):: ldiis
         type(linear_scaling_control_variables),intent(inout):: lscv
       end subroutine adjust_locregs_and_confinement

       subroutine adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)
         use module_base
         use module_types
         implicit none
         type(input_variables),intent(in):: input
         type(DFT_wavefunction),intent(in):: tmb
         type(DFT_local_fields),intent(inout) :: denspot
         type(localizedDIISParameters),intent(inout):: ldiis
         type(mixrhopotDIISParameters),intent(inout):: mixdiis
         type(linear_scaling_control_variables),intent(inout):: lscv
       end subroutine adjust_DIIS_for_high_accuracy

       subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, wfnmd, lscv)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: nlr
         type(orbitals_data),intent(in):: lorbs
         type(input_variables),intent(in):: input
         type(atoms_data),intent(in):: at
         integer,dimension(lorbs%norb),intent(in):: onwhichatom
         type(confpot_data),dimension(lorbs%norbp),intent(inout):: confdatarr
         type(wfn_metadata),intent(inout):: wfnmd
         type(linear_scaling_control_variables),intent(inout):: lscv
       end subroutine set_optimization_variables

       subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, orbsig, lzd, lzdig, op, comon)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs, orbsig
         type(local_zone_descriptors),intent(in):: lzd, lzdig
         type(overlapParameters),intent(out):: op
         type(p2pComms),intent(inout):: comon
       end subroutine determine_overlap_from_descriptors

       subroutine get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out):: weight_c, weight_f
         real(8),intent(out):: weight_c_tot, weight_f_tot
       end subroutine get_weights

       subroutine init_collective_comms(iproc, nproc, orbs, lzd, collcom, collcom_reference)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         type(collective_comms),intent(out):: collcom
         type(collective_comms),optional,intent(in):: collcom_reference
       end subroutine init_collective_comms

       subroutine deallocate_collective_comms(collcom, subname)
         use module_base
         use module_types
         implicit none
         type(collective_comms),intent(inout):: collcom
         character(len=*),intent(in):: subname
       end subroutine deallocate_collective_comms

       subroutine assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
                  istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
                  weightp_c, weightp_f, nptsp_c, nptsp_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(local_zone_descriptors),intent(in):: lzd
         real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
         real(8),intent(in):: weight_tot_c, weight_tot_f
         integer,dimension(2,0:nproc-1),intent(out):: istartend_c, istartend_f
         integer,intent(out):: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
         real(8),intent(out):: weightp_c, weightp_f
         integer,intent(out):: nptsp_c, nptsp_f
       end subroutine assign_weight_to_process

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

       !!subroutine check_gridpoint(nseg, n1, n2, noffset1, noffset2, noffset3, keyg, itarget1, itarget2, itarget3, found)
       !!  use module_base
       !!  use module_types
       !!  implicit none
       !!  integer,intent(in):: nseg, n1, n2, noffset1, noffset2, noffset3, itarget1, itarget2, itarget3
       !!  integer,dimension(2,nseg),intent(in):: keyg
       !!  logical,intent(out):: found
       !!end  subroutine check_gridpoint

       subroutine get_switch_indices(iproc, nproc, orbs, lzd, ndimpsi_c, ndimpsi_f, istartend_c, istartend_f, &
                  nsendcounts_c, nsenddspls_c, ndimind_c, nrecvcounts_c, nrecvdspls_c, &
                  nsendcounts_f, nsenddspls_f, ndimind_f, nrecvcounts_f, nrecvdspls_f, &
                  index_in_global_c, index_in_global_f, &
                  weightp_c, weightp_f,  isendbuf_c, irecvbuf_c, isendbuf_f, irecvbuf_f, &
                  indexrecvorbital_c, iextract_c, iexpand_c, indexrecvorbital_f, iextract_f, iexpand_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, ndimpsi_c, ndimpsi_f, ndimind_c, ndimind_f
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
         integer,dimension(0:nproc-1),intent(in):: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
         integer,dimension(0:nproc-1),intent(in):: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
         integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: index_in_global_c, index_in_global_f
         real(8),intent(in):: weightp_c, weightp_f
         integer,dimension(ndimpsi_c),intent(out):: isendbuf_c, irecvbuf_c
         integer,dimension(ndimpsi_f),intent(out):: isendbuf_f, irecvbuf_f
         integer,dimension(ndimind_c),intent(out):: indexrecvorbital_c, iextract_c, iexpand_c
         integer,dimension(ndimind_f),intent(out):: indexrecvorbital_f, iextract_f, iexpand_f
       end subroutine get_switch_indices

       subroutine determine_communication_arrays(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
                  index_in_global_c, index_in_global_f, &
                  weightp_c, weightp_f,  nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
                  nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
         integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: index_in_global_c, index_in_global_f
         real(8),intent(in):: weightp_c, weightp_f
         integer,dimension(0:nproc-1),intent(out):: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
         integer,dimension(0:nproc-1),intent(out):: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
       end subroutine determine_communication_arrays

       subroutine assign_weight_to_process2(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
                  npts_par_c, npts_par_f, &
                  istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
                  weightp_c, weightp_f, nptsp_c, nptsp_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(local_zone_descriptors),intent(in):: lzd
         real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
         real(8),intent(in):: weight_tot_c, weight_tot_f
         integer,dimension(0:nproc-1),intent(in):: npts_par_c, npts_par_f
         integer,dimension(2,0:nproc-1),intent(out):: istartend_c, istartend_f
         integer,intent(out):: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
         real(8),intent(out):: weightp_c, weightp_f
         integer,intent(out):: nptsp_c, nptsp_f
       end subroutine assign_weight_to_process2

       subroutine get_gridpoint_start_vectors(iproc, nproc, norbig, ndim, indexrecvbuf, weight, gridpoint_start)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norbig, ndim
         integer,dimension(ndim),intent(in):: indexrecvbuf
         real(8),dimension(norbig),intent(out):: weight
         integer,dimension(norbig),intent(out):: gridpoint_start
       end subroutine get_gridpoint_start_vectors

       subroutine get_switch_indices_vectors(iproc, nproc, nlr, norbig, ndimvec, ndimind, orbs, mlr, &
                  istartend, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, weightp, &
                  isendbuf, irecvbuf, indexrecvorbital, iextract, iexpand)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr, norbig, ndimvec, ndimind
         type(orbitals_data),intent(in):: orbs
         type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
         integer,dimension(2,0:nproc-1),intent(in):: istartend
         integer,dimension(0:nproc-1),intent(in):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
         real(8),intent(in):: weightp
         integer,dimension(ndimvec),intent(out):: isendbuf, irecvbuf
         integer,dimension(ndimind),intent(out):: indexrecvorbital, iextract, iexpand
       end subroutine get_switch_indices_vectors

       subroutine determine_communication_arrays_vectors(iproc, nproc, nlr, orbs, mlr, istartend, weightp, &
                  nsendcounts, nsenddspls, nrecvcounts, nrecvdspls)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         type(orbitals_data),intent(in):: orbs
         type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
         integer,dimension(2,0:nproc-1),intent(in):: istartend
         real(8),intent(in):: weightp
         integer,dimension(0:nproc-1),intent(out):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
       end subroutine determine_communication_arrays_vectors

       subroutine determine_num_orbs_per_gridpoint_vectors(iproc, nproc, norbig, nlr, nptsp, orbs, &
                  istartend, mlr, weightp, norb_per_gridpoint)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norbig, nlr, nptsp
         type(orbitals_data),intent(in):: orbs
         integer,dimension(2,0:nproc-1),intent(in):: istartend
         type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
         real(8),intent(in):: weightp
         integer,dimension(nptsp),intent(out):: norb_per_gridpoint
       end subroutine determine_num_orbs_per_gridpoint_vectors

       subroutine assign_weight_to_process_vectors(iproc, nproc, norbig, weight, weight_tot, istartend, weightp, nptsp)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norbig
         real(8),dimension(norbig),intent(in):: weight
         real(8),intent(in):: weight_tot
         integer,dimension(2,0:nproc-1),intent(out):: istartend
         real(8),intent(out):: weightp
         integer,intent(out):: nptsp
       end subroutine assign_weight_to_process_vectors

       subroutine get_weights_vectors(iproc, nproc, orbs, nlr, mlr, norbig, weight, weight_tot)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norbig, nlr
         type(orbitals_data),intent(in):: orbs
         type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
         real(8),dimension(norbig),intent(out):: weight
         real(8),intent(out):: weight_tot
       end subroutine get_weights_vectors

       subroutine init_collective_comms_vectors(iproc, nproc, nlr, orbs, orbsig, mlr, collcom)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         type(orbitals_data),intent(in):: orbs, orbsig
         type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
         type(collective_comms),intent(out):: collcom
       end subroutine init_collective_comms_vectors

       subroutine transpose_switch_psi(orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
         use module_base
         use module_types
         implicit none
         type(orbitals_Data),intent(in):: orbs
         type(collective_comms),intent(in):: collcom
         real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
         real(8),dimension(collcom%ndimpsi_c),intent(out):: psiwork_c
         real(8),dimension(7*collcom%ndimpsi_f),intent(out):: psiwork_f
         type(local_zone_descriptors),intent(in),optional:: lzd
       end subroutine transpose_switch_psi

       subroutine transpose_communicate_psi(iproc, nproc, collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimpsi_c),intent(in):: psiwork_c
         real(8),dimension(7*collcom%ndimpsi_f),intent(in):: psiwork_f
         real(8),dimension(collcom%ndimind_c),intent(out):: psitwork_c
         real(8),dimension(collcom%ndimind_f),intent(out):: psitwork_f
       end subroutine transpose_communicate_psi

       subroutine transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
         use module_base
         use module_types
         implicit none
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimind_c),intent(in):: psitwork_c
         real(8),dimension(7*collcom%ndimind_f),intent(in):: psitwork_f
         real(8),dimension(collcom%ndimind_c),intent(out):: psit_c
         real(8),dimension(7*collcom%ndimind_f),intent(out):: psit_f
       end subroutine transpose_unswitch_psit

       subroutine transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
         use module_base
         use module_types
         implicit none
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimind_c),intent(in):: psit_c
         real(8),dimension(7*collcom%ndimind_f),intent(in):: psit_f
         real(8),dimension(collcom%ndimind_c),intent(out):: psitwork_c
         real(8),dimension(7*collcom%ndimind_f),intent(out):: psitwork_f
       end subroutine transpose_switch_psit

       subroutine transpose_communicate_psit(iproc, nproc, collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimind_c),intent(in):: psitwork_c
         real(8),dimension(7*collcom%ndimind_f),intent(in):: psitwork_f
         real(8),dimension(collcom%ndimpsi_c),intent(out):: psiwork_c
         real(8),dimension(7*collcom%ndimpsi_f),intent(out):: psiwork_f
       end subroutine transpose_communicate_psit

       subroutine transpose_unswitch_psi(orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
         use module_base
         use module_types
         implicit none
         type(orbitals_data),intent(in):: orbs
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimpsi_c),intent(in):: psiwork_c
         real(8),dimension(7*collcom%ndimpsi_f),intent(in):: psiwork_f
         real(8),dimension(orbs%npsidim_orbs),intent(out):: psi
         type(local_zone_descriptors),intent(in),optional:: lzd
       end subroutine transpose_unswitch_psi

       subroutine transpose_localized(iproc, nproc, orbs, collcom, psi, psit_c, psit_f, lzd)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(collective_comms),intent(in):: collcom
         real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
         real(8),dimension(collcom%ndimind_c),intent(out):: psit_c
         real(8),dimension(7*collcom%ndimind_f),intent(out):: psit_f
         type(local_zone_descriptors),optional,intent(in):: lzd
       end subroutine transpose_localized

       subroutine untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, psi, lzd)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(orbitals_data),intent(in):: orbs
         type(collective_comms),intent(in):: collcom
         real(8),dimension(collcom%ndimind_c),intent(in):: psit_c
         real(8),dimension(7*collcom%ndimind_f),intent(in):: psit_f
         real(8),dimension(orbs%npsidim_orbs),intent(out):: psi
         type(local_zone_descriptors),optional,intent(in):: lzd
       end subroutine untranspose_localized

       subroutine initialize_linear_from_file(iproc,nproc,filename,iformat,Lzd,orbs,at,rxyz,orblist)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc, nproc, iformat
         type(orbitals_data), intent(inout) :: orbs  !< orbs related to the basis functions, inwhichlocreg generated in this routine
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         character(len=*), intent(in) :: filename
         type(local_zone_descriptors), intent(inout) :: Lzd !< must already contain Glr and hgrids
         integer, dimension(orbs%norb), optional :: orblist
       end subroutine initialize_linear_from_file

       subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat, &
       & locrad, locregCenter, confPotOrder, confPotprefac)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: unitwf
         logical, intent(in) :: formatted
         integer, intent(out) :: iorb_old
         integer, intent(out) :: n1_old, n2_old, n3_old
         real(gp), intent(out) :: hx_old, hy_old, hz_old
         logical, intent(out) :: lstat
         real(wp), intent(out) :: eval
         integer, intent(out) :: confPotOrder
         real(gp), intent(out) :: locrad, confPotprefac
         real(gp), dimension(3), intent(out) :: locregCenter
         character(len =256), intent(out) :: error
         ! Optional arguments
         integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
         integer, intent(in), optional :: nat
         real(gp), dimension(:,:), intent(out), optional :: rxyz_old
       end subroutine io_read_descr_linear

       subroutine readmywaves_linear(iproc,filename,iformat,norb,Lzd,orbs,at,rxyz_old,rxyz,  &
           psi,coeff,orblist)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc, iformat,norb
         type(orbitals_data), intent(inout) :: orbs  ! orbs related to the basis functions
         type(local_zone_descriptors), intent(in) :: Lzd
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
         real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psi
         character(len=*), intent(in) :: filename
         real(wp), dimension(norb,orbs%norb), intent(out) :: coeff
         integer, dimension(orbs%norb), optional :: orblist
        end subroutine readmywaves_linear

        subroutine post_p2p_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
          real(8),dimension(nsendbuf),intent(in):: sendbuf
          real(8),dimension(nrecvbuf),intent(out):: recvbuf
          type(p2pComms),intent(inout):: comm
        end subroutine post_p2p_communication

        subroutine wait_p2p_communication(iproc, nproc, comm)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(p2pComms),intent(inout):: comm
        end subroutine wait_p2p_communication

        subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi, lphiold, lhphiold)
          use module_base
          implicit none
          integer,intent(in):: npsidim
          real(8),dimension(:),pointer,intent(out):: lphi, lhphi, lphiold, lhphiold
          character(len=*),intent(in):: subname
        end subroutine allocate_auxiliary_basis_function

        subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi, lphiold, lhphiold)
          use module_base
          implicit none
          real(8),dimension(:),pointer:: lphi, lhphi, lphiold, lhphiold
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
          type(local_zone_descriptors),intent(out):: lzd_out
          character(len=*),intent(in):: subname
        end subroutine copy_local_zone_descriptors

        subroutine update_auxiliary_basis_function(subname, npsidim, lphi, lhphi, lphiold, lhphiold)
          use module_base
          implicit none
          integer,intent(in):: npsidim
          real(8),dimension(:),pointer,intent(out):: lphi, lhphi, lphiold, lhphiold
          character(len=*),intent(in):: subname
        end subroutine update_auxiliary_basis_function

        subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, n1_old, n2_old, n3_old, &
            & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: unitwf
         logical, intent(in) :: formatted
         integer, intent(out) :: norb_old, ntmb_old
         integer, intent(out) :: n1_old, n2_old, n3_old
         real(gp), intent(out) :: hx_old, hy_old, hz_old
         logical, intent(out) :: lstat
         character(len =256), intent(out) :: error
         ! Optional arguments
         integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
         integer, intent(in), optional :: nat
         real(gp), dimension(:,:), intent(out), optional :: rxyz_old
        end subroutine io_read_descr_coeff

        subroutine initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, comgp)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
          type(orbitals_data),intent(in):: orbs
          type(local_zone_descriptors),intent(in):: lzd
          type(p2pComms),intent(out):: comgp
        end subroutine initialize_communication_potential

        subroutine set_comms_ortho(iproc, nproc, orbs, lzd, op, comon)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(local_zone_descriptors),intent(in):: lzd
          type(overlapParameters),intent(inout):: op
          type(p2pComms),intent(out):: comon
        end subroutine set_comms_ortho

        subroutine nullify_overlap_parameters_matrix(opm)
          use module_base
          use module_types
          implicit none
          type(overlap_parameters_matrix),intent(inout):: opm
        end subroutine nullify_overlap_parameters_matrix

        !!subroutine nullify_p2pCommsOrthonormalityMatrix(comom)
        !!  use module_base
        !!  use module_types
        !!  implicit none
        !!  type(p2pComms),intent(inout):: comom
        !!end subroutine nullify_p2pCommsOrthonormalityMatrix

        subroutine deallocate_overlap_parameters_matrix(opm, subname)
          use module_base
          use module_types
          implicit none
          type(overlap_parameters_matrix),intent(inout):: opm
          character(len=*),intent(in):: subname
        end subroutine deallocate_overlap_parameters_matrix

        subroutine local_potential_dimensions(Lzd,orbs,ndimfirstproc)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: ndimfirstproc
          type(local_zone_descriptors), intent(inout) :: Lzd
          type(orbitals_data), intent(inout) :: orbs
        end subroutine local_potential_dimensions

        subroutine optimize_coeffs(iproc, nproc, orbs, ham, ovrlp, tmb, ldiis_coeff, fnrm)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(DFT_wavefunction),intent(inout):: tmb
          real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in):: ham, ovrlp
          type(localizedDIISParameters),intent(inout):: ldiis_coeff
          real(8),intent(out):: fnrm
        end subroutine optimize_coeffs

        subroutine DIIS_coeff(iproc, nproc, orbs, tmb, grad, coeff, ldiis)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(DFT_wavefunction),intent(in):: tmb
          real(8),dimension(tmb%orbs%norb*orbs%norb),intent(in):: grad
          real(8),dimension(tmb%orbs%norb*orbs%norb),intent(inout):: coeff
          type(localizedDIISParameters),intent(inout):: ldiis
        end subroutine DIIS_coeff

        subroutine initialize_DIIS_coeff(isx, ldiis)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: isx
          type(localizedDIISParameters),intent(out):: ldiis
        end subroutine initialize_DIIS_coeff

        subroutine allocate_DIIS_coeff(tmb, orbs, ldiis)
          use module_base
          use module_types
          implicit none
          type(DFT_wavefunction),intent(in):: tmb
          type(orbitals_data),intent(in):: orbs
          type(localizedDIISParameters),intent(out):: ldiis
        end subroutine allocate_DIIS_coeff

        subroutine initialize_DFT_local_fields(denspot)
          use module_base
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
        end subroutine initialize_DFT_local_fields

        subroutine allocate_diis_objects(idsx,alphadiis,npsidim,nkptsp,nspinor,diis,subname) !n(m)
          use module_base
          use module_types
          implicit none
          character(len=*), intent(in) :: subname
          integer, intent(in) :: idsx,npsidim,nkptsp,nspinor !n(m)
          real(gp), intent(in) :: alphadiis
          type(diis_objects), intent(inout) :: diis
        end subroutine allocate_diis_objects

        subroutine check_communications(iproc,nproc,orbs,lr,comms)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(locreg_descriptors), intent(in) :: lr
          type(communications_arrays), intent(in) :: comms
        end subroutine check_communications

        subroutine nonlocal_forces(iproc,lr,hx,hy,hz,at,rxyz,&
             orbs,nlpspd,proj,wfd,psi,fsep,refill,strten)
          use module_base
          use module_types
          implicit none
          !Arguments-------------
          type(atoms_data), intent(in) :: at
          type(wavefunctions_descriptors), intent(in) :: wfd
          type(nonlocal_psp_descriptors), intent(in) :: nlpspd
          logical, intent(in) :: refill
          integer, intent(in) :: iproc
          real(gp), intent(in) :: hx,hy,hz
          type(locreg_descriptors) :: lr
          type(orbitals_data), intent(in) :: orbs
          real(gp), dimension(3,at%nat), intent(in) :: rxyz
          real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor), intent(in) :: psi
          real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
          real(gp), dimension(3,at%nat), intent(inout) :: fsep
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
          real(gp), dimension(3,at%nat), intent(in) :: rxyz
          real(dp), dimension(*), intent(in) :: rho,pot
          real(gp), dimension(3,at%nat), intent(out) :: floc
          real(gp), dimension(6), intent(out) :: locstrten
        end subroutine local_forces

        subroutine atoms_set_symmetries(atoms, rxyz, disableSym, tol, elecfield)
          use module_base
          use module_types
          implicit none
          type(atoms_data), intent(inout) :: atoms
          real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
          logical, intent(in) :: disableSym
          real(gp), intent(in) :: tol
          real(gp), intent(in) :: elecfield(3)
        end subroutine atoms_set_symmetries

        subroutine denspot_set_history(denspot, iscf, nspin, &
             & n1i, n2i) !to be removed arguments when denspot has dimensions
          use module_types
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
          integer, intent(in) :: iscf, n1i, n2i, nspin
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

        subroutine transform_coeffs_to_derivatives(iproc, nproc, orbs, lzd, tmb, tmbder)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(local_zone_descriptors),intent(in):: lzd
          type(DFT_wavefunction),intent(in):: tmb
          type(DFT_wavefunction),intent(inout):: tmbder
        end subroutine transform_coeffs_to_derivatives

        subroutine calculate_density_kernel(iproc, nproc, ld_coeff, orbs, orbs_tmb, coeff, kernel,  ovrlp)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, ld_coeff
          type(orbitals_data),intent(in):: orbs, orbs_tmb
          real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
          real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out):: kernel
          real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),optional,intent(in):: ovrlp
        end subroutine calculate_density_kernel

        subroutine calculate_norm_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(collective_comms),intent(in):: collcom
          real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
          real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
        end subroutine calculate_norm_transposed

        subroutine reconstruct_kernel(iproc, nproc, orbs, tmb, ovrlp_tmb, overlap_calculated, kernel)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(orbitals_data),intent(in):: orbs
          type(DFT_wavefunction),intent(inout):: tmb
          real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ovrlp_tmb
          logical,intent(out):: overlap_calculated
          real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: kernel
        end subroutine reconstruct_kernel

        subroutine determine_num_orbs_per_gridpoint_new(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
                   istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
                   weightp_c, weightp_f, nptsp_c, nptsp_f, weight_c, weight_f, &
                   norb_per_gridpoint_c, norb_per_gridpoint_f)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
          type(orbitals_data),intent(in):: orbs
          type(local_zone_descriptors),intent(in):: lzd
          integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
          real(8),intent(in):: weightp_c, weightp_f
          real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
          integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
          integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
        end subroutine determine_num_orbs_per_gridpoint_new

        subroutine iguess_generator(izatom,ielpsp,zion,psppar,npspcode,ngv,ngc,nlccpar,ng,nl,&
              &   nmax_occ,noccmax,lmax,occup,expo,psiat,enlargerprb,quartic_prefactor)
           use module_base
           implicit none
           logical, intent(in) :: enlargerprb
           integer, intent(in) :: ng,npspcode,nmax_occ,lmax,noccmax,ielpsp,izatom,ngv,ngc
           real(gp), intent(in) :: zion
           integer, dimension(lmax+1), intent(in) :: nl
           !real(gp), dimension(0:4,0:6), intent(in) :: psppar
           real(gp), intent(in) :: psppar
           !real(gp), dimension(0:4,max((ngv*(ngv+1)/2)+(ngc*(ngc+1)/2),1)), intent(in) :: nlccpar
           real(gp),  intent(in) :: nlccpar
           real(gp), dimension(noccmax,lmax+1), intent(in) :: occup
           real(gp), dimension(ng+1), intent(out) :: expo
           real(gp), dimension(ng+1,nmax_occ), intent(out) :: psiat
           real(gp),intent(in),optional:: quartic_prefactor
        end subroutine iguess_generator

        subroutine penalty_basis_function(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
             confdata,ibyyzz_r) !optional
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
          integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
          !real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
          real(wp), intent(inout) :: psir !< real-space wfn in lr
          type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
          integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
        end subroutine penalty_basis_function

        subroutine allocate_convolutions_bounds(ab, subname, bounds)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: ab
          character(len=*),intent(in):: subname
          type(convolutions_bounds),intent(out):: bounds
        end subroutine allocate_convolutions_bounds

        subroutine pulay_correction(iproc, nproc, input, orbs, at, rxyz, nlpspd, proj, SIC, denspot, GPU, tmb, &
                   tmblarge, fpulay)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(input_variables),intent(in):: input
          type(orbitals_data),intent(in):: orbs
          type(atoms_data),intent(in):: at
          real(8),dimension(at%nat),intent(in):: rxyz
          type(nonlocal_psp_descriptors),intent(in):: nlpspd
          real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
          type(SIC_data),intent(in):: SIC
          type(DFT_local_fields), intent(inout) :: denspot
          type(GPU_pointers),intent(inout):: GPU
          type(DFT_wavefunction),intent(in):: tmb
          type(DFT_wavefunction),intent(inout):: tmblarge
          real(8),dimension(3,at%nat),intent(out):: fpulay
        end subroutine pulay_correction

        subroutine create_large_tmbs(iproc, nproc, tmb, eval, denspot, input, at, rxyz, lowaccur_converged, &
                   tmblarge, lhphilarge, lhphilargeold, lphilargeold)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(DFT_Wavefunction),intent(in):: tmb
          real(8),dimension(tmb%orbs%norb),intent(in):: eval
          type(DFT_local_fields),intent(in):: denspot
          type(input_variables),intent(in):: input
          type(atoms_data),intent(in):: at
          real(8),dimension(3,at%nat),intent(in):: rxyz
          logical,intent(in):: lowaccur_converged
          type(DFT_Wavefunction),intent(out):: tmblarge
          real(8),dimension(:),intent(out),pointer:: lhphilarge, lhphilargeold, lphilargeold
        end subroutine create_large_tmbs


        subroutine solvePrecondEquation(iproc,nproc,lr,ncplx,ncong,cprecr,&
             hx,hy,hz,kx,ky,kz,x,  rxyzParab, orbs, potentialPrefac, confPotOrder)
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
        end subroutine solvePrecondEquation

        subroutine derivatives_with_orthoconstraint(iproc, nproc, tmb, tmbder)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc
          type(DFT_wavefunction),intent(in):: tmb
          type(DFT_wavefunction),intent(inout):: tmbder
        end subroutine derivatives_with_orthoconstraint

        subroutine init_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work, subname)
          use module_base
          use module_types
          implicit none
          integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
          logical,intent(in):: with_confpot
          type(workarrays_quartic_convolutions),intent(inout):: work
          character(len=*),intent(in):: subname
        end subroutine init_local_work_arrays

        subroutine calculate_charge_density(iproc, nproc, lzd, hxh, hyh, hzh, orbs, densKern, factor, &
                   nscatterarr, maxval_noverlaps, noverlaps, overlaps, comarr, ise3, nrecvBuf, recvBuf, nrho, rho)
          use module_base
          use module_types
          implicit none
          integer,intent(in):: iproc, nproc, nrho, maxval_noverlaps, nrecvBuf
          real(8),intent(in):: hxh, hyh, hzh, factor
          type(local_zone_descriptors),intent(in) :: lzd
          type(orbitals_data),intent(in) :: orbs
          real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: densKern
          integer, dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
          integer,dimension(0:nproc-1),intent(in):: noverlaps
          integer,dimension(noverlaps(iproc)),intent(in):: overlaps
          integer,dimension(6,maxval_noverlaps,0:nproc-1),intent(in):: comarr
          integer,dimension(noverlaps(iproc),2),intent(in):: ise3
          real(8),dimension(nrecvBuf),intent(in):: recvBuf
          real(kind=8),dimension(nrho),intent(out),target :: rho
        end subroutine calculate_charge_density
        
        subroutine psi_to_kinpsi(iproc,orbs,lzd,psi,hpsi,ekin_sum)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: iproc
          type(orbitals_data), intent(in) :: orbs
          type(local_zone_descriptors), intent(in) :: Lzd
          real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
          real(gp), intent(out) :: ekin_sum
          real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        end subroutine psi_to_kinpsi


   end interface

END MODULE module_interfaces

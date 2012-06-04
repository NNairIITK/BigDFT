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
         type(orbitals_data), intent(out) :: orbs
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

      subroutine standard_inputfile_names(inputs, radical)
         use module_types
         implicit none
         type(input_variables), intent(out) :: inputs
         character(len = *), intent(in) :: radical
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

      subroutine read_atomic_file(file,iproc,at,rxyz,status)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: file
         integer, intent(in) :: iproc
         type(atoms_data), intent(inout) :: at
         real(gp), dimension(:,:), pointer :: rxyz
         integer, intent(out), optional :: status
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

      subroutine read_xyz_positions(iproc,ifile,atoms,rxyz,getLine)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,ifile
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
         interface
            subroutine getline(line,ifile,eof)
               integer, intent(in) :: ifile
               character(len=150), intent(out) :: line
               logical, intent(out) :: eof
            END SUBROUTINE getline
         end interface
      END SUBROUTINE read_xyz_positions

      subroutine read_ascii_positions(iproc,ifile,atoms,rxyz,getline)
         ! use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,ifile
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(:,:), pointer :: rxyz
         interface
            subroutine getline(line,ifile,eof)
               integer, intent(in) :: ifile
               character(len=150), intent(out) :: line
               logical, intent(out) :: eof
            END SUBROUTINE getline
         end interface
      END SUBROUTINE read_ascii_positions

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

      subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs,basedist)
         !n(c) use module_base
         use module_types
         implicit none
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
            &   radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj_G,proj)
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
       type(gaussian_basis),dimension(at%ntypes),intent(in) :: proj_G
         real(kind=8), dimension(:), pointer :: proj
      END SUBROUTINE createProjectorsArrays

      subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
            &   rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
         n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
         !n(c) use module_base
         use module_types
         implicit none
         !Arguments
         character(len=1), intent(in) :: datacode
         character(len=3), intent(in) :: rho_commun
         integer, intent(in) :: iproc,nproc,ixc,nspin
         real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
         type(atoms_data), intent(in) :: atoms
         type(grid_dimensions), intent(in) :: gdim
         real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
         integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
         type(rho_descriptors), intent(out) :: rhodsc
         integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
         integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
      END SUBROUTINE createDensPotDescriptors

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
         real(dp), dimension(*), intent(in) :: pkernel
         real(gp), intent(out) :: eion,edisp,psoffset
         real(dp), dimension(6),intent(out) :: ewaldstr
         real(gp), dimension(:,:), pointer :: fion,fdisp
         real(dp), dimension(*), intent(out) :: pot_ion
       END SUBROUTINE IonicEnergyandForces

       subroutine createIonicPotential(geocode,iproc,nproc,at,rxyz,&
            hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,&
            pot_ion,psoffset,rholoc)
         use module_base
         use module_types
         implicit none
         character(len=1), intent(in) :: geocode
         integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
         real(gp), intent(in) :: hxh,hyh,hzh,psoffset
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3), intent(in) :: elecfield
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(dp), dimension(*), intent(in) :: pkernel
         real(wp), dimension(*), intent(inout) :: pot_ion
         type(rholoc_objects),intent(in)::rholoc
       END SUBROUTINE createIonicPotential

       subroutine input_wf_diag(iproc,nproc,at,denspot,&
            orbs,nvirt,comms,Lzd,hx,hy,hz,rxyz,&
            nlpspd,proj,ixc,psi,hpsi,psit,G,&
            nspin,potshortcut,symObj,GPU,input,proj_G,paw)
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
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(inout) :: at
         type(orbitals_data), intent(inout) :: orbs
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(communications_arrays), intent(in) :: comms
         type(DFT_local_fields), intent(inout) :: denspot
         type(GPU_pointers), intent(inout) :: GPU
         type(input_variables):: input
         type(symmetry_data), intent(in) :: symObj
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
         type(gaussian_basis), intent(out) :: G !basis for davidson IG
         real(wp), dimension(:), pointer :: psi,hpsi,psit
         integer, intent(in) ::potshortcut
         type(gaussian_basis),dimension(at%ntypes),intent(in) :: proj_G
         type(paw_objects),intent(inout)::paw
       end subroutine input_wf_diag

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

      subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit,orthpar,paw)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(communications_arrays), intent(in) :: comms
         type(orthon_data):: orthpar
         type(paw_objects),intent(in)::paw
         real(wp), dimension(:) , pointer :: psi,hpsi,psit
      END SUBROUTINE first_orthon

      subroutine density_and_hpot(iproc,nproc,geocode,symObj,orbs,Lzd,hxh,hyh,hzh,nscatterarr,&
           pkernel,rhodsc,GPU,psi,rho,vh,hstrten)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc
        real(gp), intent(in) :: hxh,hyh,hzh
        type(rho_descriptors),intent(inout) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
        character(len=1), intent(in) :: geocode
        real(dp), dimension(*), intent(in) :: pkernel
        integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        type(GPU_pointers), intent(inout) :: GPU
        real(gp), dimension(6), intent(out) :: hstrten
        real(dp), dimension(:), pointer :: rho,vh
      END SUBROUTINE density_and_hpot

      subroutine sumrho(iproc,nproc,orbs,Lzd,hxh,hyh,hzh,nscatterarr,&
           GPU,symObj,rhodsc,psi,rho_p,mapping)
        use module_base
        use module_types
        implicit none
        !Arguments
        integer, intent(in) :: iproc,nproc
        real(gp), intent(in) :: hxh,hyh,hzh
        type(rho_descriptors),intent(in) :: rhodsc
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(symmetry_data), intent(in) :: symObj
        integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        real(dp), dimension(:,:), pointer :: rho_p
        type(GPU_pointers), intent(inout) :: GPU
        integer,dimension(orbs%norb),intent(in),optional:: mapping
      END SUBROUTINE sumrho

      !starting point for the communication routine of the density
      subroutine communicate_density(iproc,nproc,nspin,hxh,hyh,hzh,Lzd,rhodsc,nscatterarr,rho_p,rho)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        real(gp), intent(in) :: hxh,hyh,hzh
        type(local_zone_descriptors), intent(in) :: Lzd
        type(rho_descriptors),intent(in) :: rhodsc
        integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
        real(dp), dimension(:,:), pointer :: rho_p !< partial density in orbital distribution scheme
        real(dp), dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,1),1),nspin), intent(out) :: rho
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

       subroutine LocalHamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,&
            Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
            ekin_sum,epot_sum,eexctX,eSIC_DC,SIC,GPU,pkernel,orbsocc,psirocc)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd 
         type(SIC_data), intent(in) :: SIC
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
         type(confpot_data), dimension(orbs%norbp) :: confdatarr
         !real(wp), dimension(:), pointer :: pot
         real(wp), dimension(*) :: pot
         real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
         real(gp), intent(inout) :: eexctX !used to activate the OP2P scheme
         real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout) :: hpsi
         type(GPU_pointers), intent(inout) :: GPU
         real(dp), dimension(:), pointer, optional :: pkernel
         type(orbitals_data), intent(in), optional :: orbsocc
         real(wp), dimension(:), pointer, optional :: psirocc
       end subroutine LocalHamiltonianApplication

       subroutine NonLocalHamiltonianApplication(iproc,at,orbs,hx,hy,hz,rxyz,&
           proj,Lzd,nlpspd,psi,hpsi,eproj_sum,proj_G,paw)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc
        real(gp), intent(in) :: hx,hy,hz
        type(atoms_data), intent(in) :: at
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd 
        real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        real(gp), intent(out) :: eproj_sum
        type(gaussian_basis),dimension(at%ntypes),intent(in)::proj_G !projectors in gaussian basis (for PAW)
        type(paw_objects),intent(inout)::paw
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

      subroutine hpsitopsi(iproc,nproc,orbs,comms,iter,diis,idsx,psi,psit,hpsi,orthpar,&
         Lzd,paw,at,hx,hy,hz,rxyz,proj,nlpspd,eproj_sum,proj_G)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,idsx,iter
         type(communications_arrays), intent(in) :: comms
         type(orbitals_data), intent(in) :: orbs
         type(orthon_data), intent(in) :: orthpar
         type(diis_objects), intent(inout) :: diis
         real(wp), dimension(:), pointer :: psi,psit,hpsi
         type(paw_objects),intent(inout)::paw
         type(local_zone_descriptors), intent(in) :: Lzd
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd 
         type(gaussian_basis),dimension(at%ntypes),intent(in)::proj_G !projectors in gaussian basis (for PAW)
         real(gp), intent(out) :: eproj_sum
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
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

      subroutine last_orthon(iproc,nproc,orbs,wfd,&
            &   nspin,comms,psi,hpsi,psit,evsum, keeppsit)
         !n(c) use module_base
         use module_types
         implicit none
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(orbitals_data), intent(in) :: orbs
         type(communications_arrays), intent(in) :: comms
         integer, intent(in) :: iproc,nproc,nspin
         real(wp), intent(out) :: evsum
         real(wp), dimension(:) , pointer :: psi,hpsi,psit
         logical , optional :: keeppsit
      END SUBROUTINE last_orthon

      subroutine calculate_forces(iproc,nproc,Glr,atoms,orbs,nlpspd,rxyz,hx,hy,hz,proj,i3s,n3p,nspin,&
           refill_proj,ngatherarr,rho,pot,potxc,psi,fion,fdisp,fxyz,&
           ewaldstr,hstrten,xcstr,strten,fnoise,pressure,psoffset)
        use module_base
        use module_types
        implicit none
        logical, intent(in) :: refill_proj
        integer, intent(in) :: iproc,nproc,i3s,n3p,nspin
        real(gp), intent(in) :: hx,hy,hz,psoffset
        type(locreg_descriptors), intent(in) :: Glr
        type(atoms_data), intent(in) :: atoms
       type(locreg_descriptors) :: lr
        type(orbitals_data), intent(in) :: orbs
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd
        integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
        real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
        real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: rho,pot,potxc
        real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: psi
        real(gp), dimension(6), intent(in) :: ewaldstr,hstrten,xcstr
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz,fion,fdisp
        real(gp), intent(out) :: fnoise,pressure
        real(gp), dimension(6), intent(out) :: strten
        real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
      END SUBROUTINE calculate_forces
      
      subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
            &   Glr,nlpspd,ncongt,pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,&
         proj,psi,output_denspot,ekin_sum,epot_sum,eproj_sum,proj_G,paw)
         !n(c) use module_base
         use module_types
         implicit none
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         type(locreg_descriptors), intent(in) :: Glr
         type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
         type(gaussian_basis),intent(in),dimension(at%ntypes)::proj_G
         type(paw_objects),intent(inout)::paw
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
           hx,hy,hz,rxyz,rhopot,nlpspd,proj,pkernel,psi,v,dpcom,GPU)
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
        type(denspot_distribution), intent(in) :: dpcom
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
        real(dp), dimension(:), pointer :: pkernel
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
       type(atoms_data), intent(inout) :: at
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
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,nlr,norb
       integer, intent(inout) :: nvirt
       type(atoms_data), intent(inout) :: at
       type(orbitals_data), intent(in) :: orbs
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       integer,dimension(norb),intent(in):: mapping
       integer,dimension(at%nat),intent(in):: norbsPerAt
       real(gp), intent(out) :: eks
       integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
       real(gp), dimension(at%nat), intent(out) :: locrad
       type(orbitals_data), intent(out) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
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
            &   nspin,eks,scorb,G,gaucoeff,iorbtolr)
         !n(c) use module_base
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
         real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff
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
         real(dp), dimension(:), pointer :: pkernel_ref,pkernel
      END SUBROUTINE correct_hartree_potential

      subroutine xabs_lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
           &   radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
           &   ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,&
           &   in , PAWD , orbs )
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nspin
        real(gp), intent(in) :: hx,hy,hz
        type(atoms_data), intent(in), target :: at
        type(nonlocal_psp_descriptors), intent(in), target :: nlpspd
        type(local_zone_descriptors), intent(in), target :: Lzd
        type(denspot_distribution), intent(in), target :: dpcom
        real(gp), dimension(3,at%nat), intent(in), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), intent(in), target :: proj
        real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
        real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum
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
           &   radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
           &   ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,&
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
        type(denspot_distribution), intent(in), target :: dpcom
        real(gp), dimension(3,at%nat), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), target :: proj
        real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
        real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: rhoXanes

        real(gp) :: ekin_sum,epot_sum,eproj_sum
        type(GPU_pointers), intent(inout) , target :: GPU
        integer, intent(in) :: in_iat_absorber
        type(pawproj_data_type), target ::PAWD
        type(input_variables),intent(in), target :: in
        type(orbitals_data), intent(inout), target :: orbs
      end subroutine xabs_cg

      subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
           radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
           ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in, PAWD , orbs  )
        use module_base
        use module_types
        implicit none
        integer  :: iproc,nproc,nspin
        real(gp)  :: hx,hy,hz
        type(atoms_data), target :: at
        type(nonlocal_psp_descriptors), target :: nlpspd
        type(local_zone_descriptors), target :: Lzd
        type(denspot_distribution), intent(in), target :: dpcom
        real(gp), dimension(3,at%nat), target :: rxyz
        real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
        real(wp), dimension(nlpspd%nprojel), target :: proj
        real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
        real(gp) :: ekin_sum,epot_sum,eproj_sum
        type(GPU_pointers), intent(inout) , target :: GPU
        integer, intent(in) :: in_iat_absorber 
        type(input_variables),intent(in), target :: in
        type(pawproj_data_type), target ::PAWD
        type(orbitals_data), intent(inout), target :: orbs
      end subroutine xabs_chebychev

      subroutine cg_spectra(iproc,nproc,at,hx,hy,hz,rxyz,&
            &   radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
         ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in , PAWD  )! aggiunger a interface
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

         real(gp) :: ekin_sum,epot_sum,eproj_sum
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
      subroutine plot_density(filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
            &   hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nspin,nproc
         real(gp), intent(in) :: hxh,hyh,hzh
         type(atoms_data), intent(in) :: at
         integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
         real(gp), dimension(3,at%nat), intent(in) :: rxyz
         real(dp), dimension(max(n1i*n2i*n3p,1),nspin), target, intent(in) :: rho
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


      subroutine local_analysis(iproc,nproc,hx,hy,hz,in,at,rxyz,lr,orbs,orbsv,psi,psivirt)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(input_variables), intent(in) :: in
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

!      subroutine H_potential(geocode,datacode,iproc,nproc,n01,n02,n03,hx,hy,hz,&
!           rhopot,karray,pot_ion,eh,offset,sumpion,&
!           quiet,stress_tensor) !optional argument
!        use module_base
!        implicit none
!        character(len=1), intent(in) :: geocode
!        character(len=1), intent(in) :: datacode
!        logical, intent(in) :: sumpion
!        integer, intent(in) :: iproc,nproc,n01,n02,n03
!        real(gp), intent(in) :: hx,hy,hz
!        real(dp), intent(in) :: offset
!        real(dp), dimension(*), intent(in) :: karray
!        real(gp), intent(out) :: eh
!        real(dp), dimension(*), intent(inout) :: rhopot
!        real(wp), dimension(*), intent(inout) :: pot_ion
!        character(len=3), intent(in), optional :: quiet
!        real(dp), dimension(6), intent(out), optional :: stress_tensor
!      end subroutine H_potential
!
      subroutine XC_potential(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
            &   rho,exc,vxc,nspin,rhocore,potxc,xcstr,dvxcdrho,rhohat)
         use module_base
         implicit none
         character(len=1), intent(in) :: geocode
         character(len=1), intent(in) :: datacode
         integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspin
         real(gp), intent(in) :: hx,hy,hz
         real(gp), intent(out) :: exc,vxc
         real(dp), dimension(*), intent(inout) :: rho
         real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
         real(wp), dimension(:,:,:,:), pointer,optional :: rhohat
         real(wp), dimension(*), intent(out) :: potxc
         real(dp), dimension(:,:,:,:), intent(out), target, optional :: dvxcdrho
         real(dp),dimension(6),intent(out) :: xcstr
      END SUBROUTINE XC_potential
!
!
!subroutine XC_potential_test(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
!     rho,exc,vxc,nspin,rhocore,use_rhocore,potxc,xcstr,use_dvxcdrho,dvxcdrho)
!  use module_base
!  implicit none
!  character(len=1), intent(in) :: geocode
!  character(len=1), intent(in) :: datacode
!  integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspin
!  real(gp), intent(in) :: hx,hy,hz
!  real(gp), intent(out) :: exc,vxc
!  real(dp), dimension(*), intent(inout) :: rho
!  real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
!  real(wp), dimension(*), intent(out) :: potxc
!  real(dp), dimension(6), intent(out) :: xcstr
!  real(dp), dimension(:,:,:,:), intent(out), target, optional :: dvxcdrho
!  logical,intent(in)::use_rhocore,use_dvxcdrho
!end subroutine XC_potential_test

      subroutine direct_minimization(iproc,nproc,in,at,& 
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           hx,hy,hz,rxyz,rhopot,nlpspd,proj, &
           pkernel,psi,psivirt,dpcom,GPU)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(nonlocal_psp_descriptors), intent(in) :: nlpspd
        type(local_zone_descriptors),intent(inout) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        type(communications_arrays), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpcom
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
        real(dp), dimension(:), pointer :: pkernel
        real(dp), dimension(*), intent(in), target :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,psivirt
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
         real(dp), dimension(*), intent(in) :: pkernel
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

      subroutine write_eigenvalues_data(nproc,orbs,mom_vec)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: nproc
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

       subroutine full_local_potential(iproc,nproc,orbs,Lzd,iflag,dpcom,potential,pot,comgp)
         use module_base
         use module_types
         use module_xc
         implicit none
         integer, intent(in) :: iproc,nproc,iflag
         type(orbitals_data),intent(in) :: orbs
         type(local_zone_descriptors),intent(in) :: Lzd
         type(denspot_distribution), intent(in) :: dpcom
         real(wp), dimension(max(dpcom%ndimrhopot,orbs%nspin)), intent(in), target :: potential
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

    subroutine calculate_energy_and_gradient(iter,iproc,nproc,orbs,comms,GPU,Lzd,hx,hy,hz,ncong,iscf,&
            &   ekin,epot,eproj,eSIC_DC,ehart,exc,evxc,eexctX,eion,edisp,psi,psit,hpsi,gnrm,gnrm_zero,energy,&
            &   paw)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,ncong,iscf,iter
         real(gp), intent(in) :: hx,hy,hz,ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp,eSIC_DC
         type(orbitals_data), intent(in) :: orbs
         type(communications_arrays), intent(in) :: comms
      type(local_zone_descriptors), intent(in) :: Lzd
         type(GPU_pointers), intent(in) :: GPU
         real(gp), intent(out) :: gnrm,gnrm_zero,energy
         real(wp), dimension(:), pointer :: psi,psit,hpsi
         type(paw_objects),intent(inout)::paw
      END SUBROUTINE calculate_energy_and_gradient

      !!    subroutine calculate_energy_and_gradient_new(iter,iproc,nproc,orbs,comms,GPU,lr,orthpar,hx,hy,hz,ncong,iscf,&
      !!         energs,psi,psit,hpsi,gnrm,gnrm_zero,energy)
      !!      use module_base
      !!      use module_types
      !!      !use module_interfaces!, except_this_one => calculate_energy_and_gradient_new
      !!      implicit none
      !!      integer, intent(in) :: iproc,nproc,ncong,iscf,iter
      !!      real(gp), intent(in) :: hx,hy,hz
      !!      type(orbitals_data), intent(inout) :: orbs
      !!      type(communications_arrays), intent(in) :: comms
      !!      type(locreg_descriptors), intent(in) :: lr
      !!      type(GPU_pointers), intent(in) :: GPU
      !!      type(orthon_data), intent(in) :: orthpar
      !!      type(energy_terms), intent(inout) :: energs
      !!      real(gp), intent(out) :: gnrm,gnrm_zero,energy
      !!      real(wp), dimension(:), pointer :: psi,psit,hpsi
      !!    END SUBROUTINE calculate_energy_and_gradient_new
      subroutine constrained_davidson(iproc,nproc,in,at,& 
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           hx,hy,hz,rxyz,rhopot,psi,v,dpcom,GPU)
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
        type(denspot_distribution), intent(in) :: dpcom
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
        !v, that is psivirt, is transposed on input and direct on output
      end subroutine constrained_davidson

      subroutine local_hamiltonian(iproc,orbs,Lzd,hx,hy,hz,&
           ipotmethod,confdatarr,pot,psi,hpsi,pkernel,ixc,alphaSIC,ekin_sum,epot_sum,eSIC_DC)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,ipotmethod,ixc
        real(gp), intent(in) :: hx,hy,hz,alphaSIC
        type(orbitals_data), intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
        real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi !this dimension will be modified
        real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
        !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
        real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
        real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
        real(dp), dimension(:), pointer :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
      END SUBROUTINE local_hamiltonian

      subroutine NK_SIC_potential(lr,orbs,ixc,fref,hxh,hyh,hzh,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
         !n(c) use module_base
         use module_types
         implicit none
         integer, intent(in) :: ixc
         real(gp), intent(in) :: hxh,hyh,hzh,fref
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         real(dp), dimension(*), intent(in) :: pkernel
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
         real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
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

      subroutine denspot_communications(iproc,nproc,grid,hxh,hyh,hzh,in,atoms,rxyz,radii_cf,dpcom,rhod)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc, nproc
        type(grid_dimensions), intent(in) :: grid
        real(gp), intent(in) :: hxh, hyh, hzh
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: atoms
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
        real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
        type(denspot_distribution), intent(inout) :: dpcom
        type(rho_descriptors), intent(out) :: rhod
      end subroutine denspot_communications

      subroutine allocateRhoPot(iproc,nproc,Glr,hxh,hyh,hzh,in,atoms,rxyz,radii_cf,denspot)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc, nproc
        type(locreg_descriptors), intent(in) :: Glr
        real(gp), intent(in) :: hxh, hyh, hzh
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: atoms
        real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
        real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
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

      subroutine getLocalizedBasis(iproc,nproc,at,lzd,lorbs,orbs,comon,op,comgp,mad,rxyz,&
           denspot,GPU,lphi,trH,&
           infoBasisFunctions,ovrlp,nlpspd,proj,coeff,ldiis,nit,nItInnerLoop,newgradient,orthpar,&
           confdatarr,methTransformOverlap,blocksize_pdgemm,convCrit,hx,hy,hz,SIC,nItPrecond)
        use module_base
        use module_types

        implicit none

        ! Calling arguments
        integer,intent(in):: iproc, nproc, nit, nItInnerLoop, methTransformOverlap, blocksize_pdgemm
        integer,intent(in):: nItPrecond
        integer,intent(out):: infoBasisFunctions
        type(atoms_data), intent(in) :: at
        type(local_zone_descriptors),intent(inout):: lzd
        type(orbitals_data):: lorbs, orbs
        type(p2pComms):: comon
        type(overlapParameters):: op
        !type(p2pCommsGatherPot):: comgp
        type(p2pComms):: comgp
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(3,at%nat):: rxyz
        type(DFT_local_fields), intent(inout) :: denspot
        !integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
        !integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
        !real(dp), dimension(*), intent(inout) :: rhopot
        type(GPU_pointers), intent(inout) :: GPU
        !real(dp), dimension(:), pointer :: pkernelseq
        real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)):: lphi
        real(8),intent(out):: trH
        real(8),intent(in):: convCrit, hx, hy, hz
        real(8),dimension(lorbs%norb,lorbs%norb),intent(out):: ovrlp
        type(nonlocal_psp_descriptors),intent(in):: nlpspd
        real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
        real(8),dimension(lorbs%norb,orbs%norb),intent(in):: coeff
        type(localizedDIISParameters),intent(inout):: ldiis
        logical,intent(in):: newgradient
        type(orthon_data),intent(in):: orthpar
        type(confpot_data), dimension(lorbs%norbp),intent(in) :: confdatarr
        type(SIC_data) :: SIC !<parameters for the SIC methods
      end subroutine getLocalizedBasis



    subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, &
          input, rxyz, nscatterarr, tag, coeff, lphi)
      use module_base
      use module_types
      implicit none
      ! Calling arguments
      integer,intent(in):: iproc, nproc
      type(locreg_descriptors),intent(in):: Glr
      type(orbitals_data),intent(in):: orbs
      type(atoms_data),intent(inout):: at
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      type(linearParameters),intent(inout):: lin
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(in):: rxyz
      integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,intent(inout):: tag
      real(8),dimension(:,:),pointer,intent(out):: coeff
      real(8),dimension(:),pointer,intent(out):: lphi
    end subroutine allocateAndInitializeLinear



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
      real(dp), dimension(:), pointer :: pkernel,pkernelseq
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
    real(dp), dimension(:), pointer :: pkernelseq
    real(8),dimension(max(orbsLIN%npsidim_orbs,orbsLIN%npsidim_comp)):: phi
    real(8),dimension(3,at%nat):: perturbation
    end subroutine estimatePerturbedOrbitals
    
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

    subroutine getLinearPsi(iproc,nproc,lzd,orbs,lorbs,llborbs,comsr,&
         mad,lbmad,op,lbop,comon,lbcomon,comgp,lbcomgp,at,rxyz,denspot,&
         GPU,updatePhi,&
         infoBasisFunctions,infoCoeff,itSCC,ebs,coeff,lphi,nlpspd,proj,communicate_lphi,coeff_proj,&
         ldiis,nit,nItInnerLoop,newgradient,orthpar,confdatarr,&
         methTransformOverlap,blocksize_pdgemm,convCrit,nItPrecond,&
         useDerivativeBasisFunctions,lphiRestart,comrp,blocksize_pdsyev,nproc_pdsyev,&
         hx,hy,hz,SIC)
      use module_base
      use module_types
      implicit none

      ! Calling arguments
      integer,intent(in):: iproc, nproc, itSCC, nit, nItInnerLoop
      integer,intent(in):: methTransformOverlap, blocksize_pdgemm, nItPrecond
      integer,intent(in):: blocksize_pdsyev, nproc_pdsyev
      type(local_zone_descriptors),intent(inout):: lzd
      type(orbitals_data),intent(in) :: orbs, lorbs, llborbs
      !type(p2pCommsSumrho),intent(inout):: comsr
      type(p2pComms),intent(inout):: comsr
      type(matrixDescriptors),intent(in):: mad, lbmad
      type(overlapParameters),intent(inout):: op, lbop
      type(p2pComms),intent(inout):: comon, lbcomon
      !type(p2pCommsGatherPot):: comgp, lbcomgp
      type(p2pComms):: comgp, lbcomgp
      type(atoms_data),intent(in):: at
      real(8),dimension(3,at%nat),intent(in):: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers),intent(inout):: GPU
      logical,intent(in):: updatePhi, newgradient, useDerivativeBasisFunctions
      integer,intent(out):: infoBasisFunctions, infoCoeff
      real(8),intent(out):: ebs
      real(8),intent(in):: convCrit, hx, hy, hz
      real(8),dimension(llborbs%norb,orbs%norb),intent(in out):: coeff
      real(8),dimension(max(llborbs%npsidim_orbs,llborbs%npsidim_comp)),intent(inout):: lphi
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      logical,intent(in):: communicate_lphi
      real(8),dimension(lorbs%norb,orbs%norb),intent(inout):: coeff_proj
      type(localizedDIISParameters),intent(inout):: ldiis
      type(orthon_data),intent(in):: orthpar
      type(confpot_data),dimension(lorbs%norbp),intent(in) :: confdatarr
      real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout)::lphiRestart
      type(p2pCommsRepartition),intent(inout):: comrp
      type(SIC_data),intent(in):: SIC
    end subroutine getLinearPsi


    subroutine local_hamiltonianConfinement(iproc,orbs,lin,lr,hx,hy,hz,&
         nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, onWhichAtom, at, centralAtom)
      use module_base
      use module_types
      use libxc_functionals
      implicit none
      integer, intent(in) :: iproc,nspin
      real(gp), intent(in) :: hx,hy,hz
      type(orbitals_data), intent(in) :: orbs
      type(linearParameters):: lin
      type(locreg_descriptors), intent(in) :: lr
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
      real(wp), dimension(*) :: pot
      !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
      real(gp), intent(out) :: ekin_sum,epot_sum
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
    integer:: nat
    real(8),dimension(3,nat):: rxyz
    integer,dimension(orbs%norbp),intent(in):: onWhichAtom
    type(atoms_data), intent(in) :: at
    integer,intent(in),optional:: centralAtom
    end subroutine local_hamiltonianConfinement


    subroutine local_hamiltonianConfinementForAllLocregs(iproc,orbs,lin,lr,hx,hy,hz,&
         nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, at, centralAtom)
      use module_base
      use module_types
      use libxc_functionals
      implicit none
      integer, intent(in):: iproc, nspin, nat
      real(gp), intent(in):: hx,hy,hz
      type(orbitals_data),intent(in):: orbs
      type(linearParameters),intent(in):: lin
      type(locreg_descriptors),intent(in) :: lr
      type(atoms_data),intent(in):: at
      real(wp),dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp),intent(in):: psi
      real(wp),dimension(*):: pot
      !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
      real(gp),intent(out):: ekin_sum,epot_sum
      real(wp),dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp,at%nat),intent(out):: hpsi
      real(8),dimension(3,nat),intent(in):: rxyz
      integer,intent(in),optional:: centralAtom
    end subroutine local_hamiltonianConfinementForAllLocregs

    
    subroutine deallocateLinear(iproc, lin, lphi, coeff)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc
      type(linearParameters),intent(inout):: lin
      real(8),dimension(:),pointer,intent(inout):: lphi
      real(8),dimension(:,:),pointer,intent(inout):: coeff
    end subroutine deallocateLinear
    

    subroutine initializeLocRegLIN(iproc, nproc, lr, lin, at, input, rxyz, radii_cf)
      use module_base
      use module_types
      implicit none
      integer:: iproc, nproc
      type(locreg_descriptors):: lr
      type(linearParameters):: lin
      type(atoms_data),intent(in):: at
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat):: rxyz
      real(8),dimension(at%ntypes,3):: radii_cf
      type(communications_arrays):: commsLIN
    end subroutine initializeLocRegLIN
    
    
    
    subroutine orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, newComm, norbPerComm)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc, lproc, uproc, newComm, norbPerComm
      type(linearParameters),intent(in out):: lin
    end subroutine orbitalsCommunicatorsWithGroups
    
    subroutine linearScaling(iproc,nproc,Glr,orbs,comms,at,input,&
         lin,rxyz,fion,fdisp,denspot,nlpspd,proj,GPU,&
         eion,edisp,eexctX,scpot,psi,psit,energy,fxyz)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(locreg_descriptors),intent(in) :: Glr
      type(orbitals_data),intent(inout):: orbs
      type(communications_arrays),intent(in) :: comms
      type(atoms_data),intent(inout):: at
      type(linearParameters),intent(inout):: lin
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(inout):: rxyz
      real(8),dimension(3,at%nat),intent(in):: fion, fdisp
      type(DFT_local_fields), intent(inout) :: denspot
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      type(GPU_pointers),intent(in out):: GPU
      real(gp),intent(in):: eion, edisp, eexctX
      logical,intent(in):: scpot
      !real(8),dimension(orbs),intent(out):: psi
      real(8),dimension(:),pointer,intent(out):: psi, psit
      real(gp), dimension(:), pointer :: rho,pot
      real(8),intent(out):: energy
      real(8),dimension(3,at%nat),intent(out):: fxyz
    end subroutine linearScaling
    
    
    subroutine potentialAndEnergySub(iproc, nproc, n3d, n3p, Glr, orbs, atoms, in, lin, phi, psi, rxyz, rxyzParab, &
        rhopot, nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
        proj, nlpspd, pkernelseq, eion, edisp, eexctX, scpot, coeff, ebsMod, energy)
      use module_base
      use module_types
      implicit none
      integer:: iproc, nproc, n3d, n3p, sizeLphir, sizePhibuffr
      type(locreg_descriptors) :: Glr
      type(orbitals_data):: orbs
      type(atoms_data):: atoms
      type(input_variables):: in
      type(linearParameters):: lin
      real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)):: phi
      real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)):: psi
      real(dp), dimension(lin%as%size_rhopot) :: rhopot
      integer,dimension(0:nproc-1,4) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
      type(GPU_pointers),intent(in out):: GPU
      integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)) :: irrzon
      real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)) :: phnons
      real(dp), dimension(lin%as%size_pkernel):: pkernel
      real(wp), dimension(lin%as%size_pot_ion):: pot_ion
      !real(wp), dimension(lin%as%size_rhocore):: rhocore 
      real(wp), dimension(:),pointer:: rhocore
      real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)):: potxc
      character(len=3):: PSquiet
      type(nonlocal_psp_descriptors),intent(in) :: nlpspd
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(dp),dimension(lin%as%size_pkernelseq),intent(in):: pkernelseq
      real(8),dimension(3,atoms%nat),intent(in):: rxyz
      real(8),dimension(3,atoms%nat),intent(in):: rxyzParab
      real(gp):: eion, edisp, eexctX, energy
      real(8),dimension(lin%lb%orbs%norb,orbs%norb):: coeff
      real(8):: ebsMod
      logical:: scpot
    end subroutine potentialAndEnergySub

    
    

    subroutine calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, atoms, in, comms, lin, nlpspd, proj, &
        ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, phi, coeff, rhopot, fxyz, fnoise, radii_cf)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
      type(locreg_descriptors),intent(in):: Glr
      type(orbitals_data),intent(in):: orbs
      type(atoms_data),intent(in):: atoms
      type(input_variables),intent(in):: in
      type(communications_arrays),intent(in):: comms
      type(linearParameters),intent(in):: lin
      type(nonlocal_psp_descriptors),intent(in) :: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
      integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr   !!! NOT NEEDED
      integer,dimension(0:nproc-1,4),intent(inout) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(GPU_pointers),intent(inout):: GPU
      integer,dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
      real(dp),dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
      real(dp),dimension(lin%as%size_pkernel),intent(in):: pkernel
      real(8),dimension(3,atoms%nat),intent(in):: rxyz, fion, fdisp
      real(8),dimension(3,atoms%nat),intent(out):: fxyz
      real(8),intent(out):: fnoise
      real(8),dimension(max(lin%gorbs%npsidim_orbs,lin%gorbs%npsidim_comp)),intent(inout):: phi
      real(8),dimension(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1)),intent(in):: rhopot
      real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
      real(gp), dimension(atoms%ntypes,3+ndebug), intent(in) :: radii_cf
    end subroutine calculateForcesSub



    subroutine optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nspin
      type(orbitals_data),intent(in):: orbs
      type(linearParameters),intent(in):: lin
      real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: matrixElements
      real(8),dimension(lin%orbs%norb,orbs%norb),intent(inout):: coeff
      integer,intent(out):: infoCoeff
    end subroutine

   !!subroutine determine_locreg_periodic(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,outofzone)
   !!   use module_base
   !!   use module_types
   !!   implicit none
   !!   integer, intent(in) :: iproc,nlr
   !!   real(gp), intent(in) :: hx,hy,hz
   !!   type(locreg_descriptors), intent(in) :: Glr
   !!   real(gp), dimension(nlr), intent(in) :: locrad
   !!   real(gp), dimension(3,nlr), intent(in) :: cxyz
   !!   type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
   !!   integer, dimension(3,nlr),intent(out) :: outofzone
   !!end subroutine

    !!subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr,outofzone)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  integer,intent(in) :: ilr,nlr
    !!  type(locreg_descriptors),intent(in) :: Glr
    !!  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr
    !!  integer,dimension(3,nlr),intent(in) :: outofzone
    !!end subroutine

    !!subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
    !! nseg_loc,nvctr_loc,outofzone)
    !! implicit none
    !! integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
    !! integer, dimension(nseg), intent(in) :: keyv
    !! integer, dimension(2,nseg), intent(in) :: keyg
    !! integer, intent(out) :: nseg_loc,nvctr_loc
    !! integer, dimension(3),intent(in) :: outofzone
    !!end subroutine

    !!subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
    !! nseg_loc,nvctr_loc,keyg_loc,keyv_loc,outofzone)
    !! implicit none
    !! integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
    !! integer, dimension(nseg), intent(in) :: keyv
    !! integer, dimension(2,nseg), intent(in) :: keyg
    !! integer, dimension(3), intent(in) :: outofzone
    !! integer, dimension(nseg_loc), intent(out) :: keyv_loc
    !! integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
    !! end subroutine

    !!subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,outofzone)
    !! use module_base
    !! use module_types
    !! implicit none
    !! integer, intent(in) :: alr,blr
    !! integer, intent(in) :: nlr
    !! type(locreg_descriptors),intent(in) :: Glr
    !! integer, intent(out) :: isovrlp
    !! integer,dimension(3,nlr),intent(in) :: outofzone
    !! type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
    !!end subroutine

    !!subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr,outofzone)
    !! use module_base
    !! use module_types
    !! implicit none
    !! integer, intent(in) :: alr,blr
    !! integer, intent(in) :: nlr
    !! type(locreg_descriptors),intent(in) :: Glr
    !! integer, intent(in) :: isovrlp
    !! type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
    !! type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr
    !! integer,dimension(3,nlr),intent(in) :: outofzone
    !!end subroutine

    !!subroutine nlpspd_to_locreg(input_parameters,iproc,Glr,Llr,rxyz,atoms,orbs,&
    !!   radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,Lnlpspd,projflg)
    !! use module_base
    !! use module_types
    !! implicit none
    !! type(input_variables),intent(in) :: input_parameters
    !! integer,intent(in) :: iproc
    !! type(locreg_descriptors),intent(in) :: Glr
    !! type(locreg_descriptors),intent(in) :: Llr
    !! type(atoms_data),intent(in) :: atoms
    !! type(orbitals_data),intent(in) :: orbs
    !! real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
    !! type(nonlocal_psp_descriptors),intent(in) :: nlpspd
    !! type(nonlocal_psp_descriptors),intent(out) :: Lnlpspd
    !! integer,dimension(atoms%nat),intent(out) :: projflg
    !! real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
    !! real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
    !!end subroutine

    !!subroutine apply_local_projectors(atoms,in,Llr,Lnlpspd,Lproj,orbs,projflg,psi,rxyz,hpsi)
    !! use module_base
    !! use module_types
    !! implicit none
    !! type(atoms_data),intent(in) :: atoms
    !! type(input_variables),intent(in) :: in
    !! type(locreg_descriptors),intent(in) :: Llr
    !! type(nonlocal_psp_descriptors),intent(in) :: Lnlpspd
    !! type(orbitals_data),intent(in) :: orbs
    !! integer,dimension(atoms%nat),intent(in) :: projflg
    !! real(wp),dimension(Lnlpspd%nprojel),intent(out):: Lproj
    !! real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(in) :: psi
    !! real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(out):: hpsi
    !! real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
    !!end subroutine

    !!subroutine psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)
    !! use module_base
    !! use module_types
    !! implicit none
    !! integer, intent(in) :: nlr
    !! integer :: ilr
    !! integer :: ldim
    !! type(orbitals_data),intent(in) :: orbs
    !! type(locreg_descriptors),intent(in) :: Glr
    !! type(locreg_descriptors), dimension(nlr), intent(in) :: Olr
    !! real(wp),dimension(orbs%npsidim),intent(in) :: psi
    !! real(wp),dimension(ldim),intent(inout) :: lpsi
    !!end subroutine





    !!subroutine partial_density_linear(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
    !!     hfac,nscatterarr,spinsgn,psir,rho_p,norb,norbPsi,coeff,&
    !!     ibyyzz_r)
    !!  use module_base
    !!  use module_types
    !!  implicit none
    !!  logical, intent(in) :: rsflag
    !!  integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir, norb,norbPsi
    !!  real(gp), intent(in) :: hfac,spinsgn
    !!  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
    !!  real(wp), dimension(n1i,n2i,n3i,npsir), intent(in) :: psir
    !!  real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
    !!  real(8),dimension(norb,norbPsi),intent(in):: coeff
    !!  integer, dimension(:,:,:),pointer :: ibyyzz_r
    !!end subroutine partial_density_linear

!!$    subroutine local_partial_densityLinear(iproc,nproc,rsflag,nscatterarr,&
!!$         nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,psi,rho)
!!$      use module_base
!!$      use module_types
!!$      use module_interfaces, exceptThisOne => local_partial_densityLinear
!!$      use module_xc
!!$      implicit none
!!$      logical, intent(in) :: rsflag
!!$      integer, intent(in) :: iproc,nproc
!!$      integer,intent(inout):: nrhotot
!!$      integer, intent(in) :: nspin
!!$      real(gp), intent(in) :: hxh,hyh,hzh
!!$      type(local_zone_descriptors), intent(in) :: Lzd
!!$      type(orbitals_data),intent(in) :: orbs
!!$      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!$      real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
!!$      real(dp),dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1),max(nspin,orbs%nspinor)),intent(out):: rho
!!$    end subroutine local_partial_densityLinear


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


subroutine HamiltonianApplicationConfinementForAllLocregs(iproc,nproc,at,orbs,lin,hx,hy,hz,rxyz,&
     nlpspd,proj,lr,ngatherarr,ndimpot,potential,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParabola, onWhichAtom, &
     pkernel,orbsocc,psirocc,centralAtom) ! optional
      use module_base
      use module_types
      use libxc_functionals
      implicit none
      integer, intent(in) :: iproc,nproc,ndimpot,nspin
      real(gp), intent(in) :: hx,hy,hz
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs
      type(linearParameters):: lin
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(locreg_descriptors), intent(in) :: lr
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
      real(wp), dimension(max(ndimpot,1)*nspin), intent(in), target :: potential
      real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
      real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp*at%nat), intent(out) :: hpsi
      type(GPU_pointers), intent(inout) :: GPU
      real(gp), dimension(3,at%nat), intent(in) :: rxyzParabola
      integer,dimension(orbs%norb),intent(in):: onWhichAtom
      real(dp), dimension(*), optional :: pkernel
      type(orbitals_data), intent(in), optional :: orbsocc
      real(wp), dimension(:), pointer, optional :: psirocc
      integer,intent(in),optional:: centralAtom
    end subroutine HamiltonianApplicationConfinementForAllLocregs


    subroutine readAtomicOrbitals(at,norbe,norbsc,nspin,nspinor,scorb,norbsc_arr,locrad)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: nspin,nspinor
      integer, intent(out) :: norbe,norbsc
      type(atoms_data), intent(inout) :: at
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

    subroutine inputguessConfinement(iproc, nproc, at, &
         comms, Glr, input, lin, orbs, rxyz,denspot, rhopotold,&
         nlpspd, proj, GPU,  &
         tag, lphi, ehart, eexcu, vexcu)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: iproc,nproc
      type(atoms_data), intent(inout) :: at
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(locreg_descriptors), intent(in) :: Glr
      type(communications_arrays), intent(in) :: comms
      type(GPU_pointers), intent(inout) :: GPU
      type(DFT_local_fields), intent(inout) :: denspot
      type(input_variables):: input
      type(linearParameters),intent(inout):: lin
      type(orbitals_data),intent(in):: orbs
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
      real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin),intent(inout) ::  rhopotold
      integer,intent(inout):: tag
      real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(out):: lphi
      real(8),intent(out):: ehart, eexcu, vexcu
    end subroutine inputguessConfinement

    !subroutine sumrhoForLocalizedBasis(iproc, nproc, orbs, Glr, input, lin, coeff, phi, nrho, rho, &
    !           at, rxyz, nscatterarr, phibuff)
    !  use module_base
    !  use module_types
    !  use libxc_functionals
    !  implicit none
    !  ! Calling arguments
    !  integer,intent(in):: iproc, nproc, nrho
    !  type(orbitals_data),intent(in):: orbs
    !  type(locreg_descriptors),intent(in):: Glr
    !  type(input_variables),intent(in):: input
    !  type(linearParameters),intent(inout):: lin
    !  real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
    !  real(8),dimension(lin%orbs%npsidim),intent(in):: phi
    !  real(8),dimension(nrho),intent(out),target:: rho
    !  type(atoms_data),intent(in):: at
    !  real(8),dimension(3,at%nat),intent(in):: rxyz
    !  integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
    !  real(8),dimension(lin%comsr%sizePhibuff):: phibuff
    !end subroutine sumrhoForLocalizedBasis


    !!!subroutine initializeCommsSumrho(iproc, nproc, nscatterarr, lin, phibuff)
    !!!  use module_base
    !!!  use module_types
    !!!  implicit none
    !!!  
    !!!  ! Calling arguments
    !!!  integer,intent(in):: iproc, nproc
    !!!  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
    !!!  type(linearParameters),intent(inout):: lin
    !!!  real(8),dimension(:),pointer,intent(out):: phibuff
    !!!end subroutine initializeCommsSumrho


    subroutine initializeCommsSumrho(iproc, nproc, nscatterarr, lzd, orbs, tag, comsr)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc
      integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs
      integer,intent(inout):: tag
      !type(p2pCommsSumrho),intent(out):: comsr
      type(p2pComms),intent(out):: comsr
    end subroutine initializeCommsSumrho


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

    subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)
     use module_base
     use module_types
     implicit none
     integer, intent(in) :: alr,blr              
     integer, intent(in) :: nlr                  
     type(locreg_descriptors),intent(in) :: Glr  
     integer, intent(out) :: isovrlp           
     type(locreg_descriptors), dimension(nlr), intent(in) :: Llr       
    end subroutine get_number_of_overlap_region

    subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)
     use module_base
     use module_types
     implicit none
     integer, intent(in) :: alr,blr           
     integer, intent(in) :: nlr                
     type(locreg_descriptors),intent(in) :: Glr 
     integer, intent(in) :: isovrlp              
     type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  
     type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr 
    end subroutine get_overlap_region_periodic

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

    subroutine local_partial_densityLinear(iproc,nproc,rsflag,nscatterarr,&
         nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,mapping,psi,rho)
      use module_base
      use module_types
      use module_xc
      implicit none
      logical, intent(in) :: rsflag
      integer, intent(in) :: iproc,nproc
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

!!$     subroutine LinearHamiltonianApplication(input,iproc,nproc,at,Lzd,orbs,hx,hy,hz,rxyz,&
!!$        proj,ngatherarr,pot,psi,hpsi,&
!!$        ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf,pkernel,orbsocc,psirocc)
!!$       use module_base
!!$       use module_types
!!$       use libxc_functionals
!!$       implicit none
!!$       integer, intent(in) :: iproc,nproc,nspin
!!$       real(gp), intent(in) :: hx,hy,hz
!!$       type(atoms_data), intent(in) :: at
!!$       type(input_variables), intent(in) :: input
!!$       type(local_zone_descriptors),intent(inout) :: Lzd
!!$       type(orbitals_data),intent(in) :: orbs
!!$       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!$       real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!$       real(wp), dimension(Lzd%Gnlpspd%nprojel), intent(in) :: proj
!!$       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
!!$       real(wp), dimension(:), pointer :: pot
!!$       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
!!$       real(wp), target, dimension(orbs%npsidim_orbs), intent(out) :: hpsi
!!$       type(GPU_pointers), intent(inout) :: GPU
!!$       real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
!!$       real(dp), dimension(*), optional :: pkernel
!!$       type(orbitals_data), intent(in), optional :: orbsocc
!!$       real(wp), dimension(:), pointer, optional :: psirocc
!!$     end subroutine LinearHamiltonianApplication

     
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

     subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,comms,&
          psi,hpsi,psit,orthpar,passmat,& !mandatory
          orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,natsc,nspin
       type(local_zone_descriptors) :: Lzd                                  !> Information about the locregs
       type(communications_arrays), target, intent(in) :: comms
       type(orbitals_data), target, intent(inout) :: orbs
       type(input_variables):: input
       type(orthon_data):: orthpar
       real(wp), dimension(*), intent(out) :: passmat !< passage matrix for building the eigenvectors (the size depends of the optional arguments)
       real(wp), dimension(:), pointer :: psi,hpsi,psit
       !optional arguments
       real(gp), optional, intent(in) :: etol
       type(orbitals_data), optional, intent(in) :: orbsv
       type(orbitals_data), optional, target, intent(in) :: orbse
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

     subroutine orthonormalizeOnlyDerivatives(iproc, nproc, lin, phid)
       use module_base
       use module_defs
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(linearParameters),intent(in):: lin
       real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)),intent(inout):: phid
     end subroutine orthonormalizeOnlyDerivatives

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

     subroutine sumrhoForLocalizedBasis2(iproc, nproc, norb, lzd, input, orbs, comsr, coeff, nrho, rho, at, nscatterarr)
       use module_base
       use module_types
       use libxc_functionals
       implicit none
       integer,intent(in):: iproc, nproc, nrho, norb
       type(local_zone_descriptors),intent(in):: lzd
       type(input_variables),intent(in):: input
       type(orbitals_data),intent(in):: orbs
       !type(p2pCommsSumrho),intent(inout):: comsr
       type(p2pComms),intent(inout):: comsr
       real(8),dimension(orbs%norb,norb),intent(in):: coeff
       real(8),dimension(nrho),intent(out),target:: rho
       type(atoms_data),intent(in):: at
       integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     end subroutine sumrhoForLocalizedBasis2


     subroutine postCommunicationSumrho2(iproc, nproc, comsr, sendBuf, recvBuf)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       !type(p2pCommsSumrho),intent(inout):: comsr
       type(p2pComms),intent(inout):: comsr
       real(8),dimension(comsr%nsendBuf),intent(inout):: sendBuf
       real(8),dimension(comsr%nrecvBuf),intent(out):: recvBuf
     end subroutine postCommunicationSumrho2


     subroutine allocateLinArrays(lin)
       use module_base
       use module_types
       implicit none
       type(linearParameters),intent(inout):: lin
     end subroutine allocateLinArrays


     !!subroutine initLocregs(iproc, nat, rxyz, lin, input, Glr)
     !!  use module_base
     !!  use module_types
     !!  implicit none
     !!  integer,intent(in):: iproc, nat
     !!  real(8),dimension(3,nat),intent(in):: rxyz
     !!  type(linearParameters),intent(inout):: lin
     !!  type(input_variables),intent(in):: input
     !!  type(locreg_descriptors),intent(in):: Glr
     !!end subroutine initLocregs


     subroutine initCoefficients(iproc, orbs, lin, coeff)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc
       type(orbitals_data),intent(in):: orbs
       type(linearParameters),intent(in):: lin
       real(8),dimension(:,:),pointer,intent(out):: coeff
     end subroutine initCoefficients



!!$     subroutine HamiltonianApplicationConfinement2(input,iproc,nproc,at,Lzd,orbs,lin,hx,hy,hz,rxyz,&
!!$          ngatherarr,ndimpot,pot,psi,hpsi,&
!!$          ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf, comgp, onWhichAtomp, withConfinement, energyReductionFlag, &
!!$          doNotCalculate, pkernel,orbsocc,psirocc)
!!$       use module_base
!!$       use module_types
!!$       use libxc_functionals
!!$       implicit none
!!$       integer, intent(in) :: iproc,nproc,nspin,ndimpot
!!$       real(gp), intent(in) :: hx,hy,hz
!!$       type(atoms_data), intent(in) :: at
!!$       type(input_variables), intent(in) :: input
!!$       type(local_zone_descriptors),intent(inout) :: Lzd
!!$       type(orbitals_data),intent(in):: orbs
!!$       type(linearParameters),intent(in):: lin
!!$       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!$       real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!$       real(wp), dimension(orbs%npsidim), intent(in) :: psi
!!$       real(wp), dimension(max(ndimpot,1)*nspin), intent(in) :: pot
!!$       !real(wp), dimension(:), pointer :: pot
!!$       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
!!$       real(wp), target, dimension(orbs%npsidim), intent(out) :: hpsi
!!$       type(GPU_pointers), intent(inout) :: GPU
!!$       real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
!!$       type(p2pCommsGatherPot), intent(in):: comgp
!!$       integer,dimension(orbs%norbp),intent(in):: onWhichAtomp
!!$       logical,intent(in):: withConfinement
!!$       logical,intent(in):: energyReductionFlag
!!$       logical,dimension(lzd%nlr),intent(in),optional:: doNotCalculate
!!$       real(dp), dimension(*), optional :: pkernel
!!$       type(orbitals_data), intent(in), optional :: orbsocc
!!$       real(wp), dimension(:), pointer, optional :: psirocc
!!$     end subroutine HamiltonianApplicationConfinement2


     subroutine local_hamiltonian_LinearConfinement(iproc, nproc, ilr, orbs, lr, norb, hx, hy, hz, &
          nspin, ndimpot, pot, psi, hpsi, ekin_sum, epot_sum, lin, at, rxyz, onWhichAtomp, withConfinement)
       use module_base
       use module_types
       use libxc_functionals
       implicit none
       integer, intent(in) :: iproc, nproc, nspin, ilr, norb, ndimpot
       real(gp), intent(in) :: hx, hy, hz
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: lr
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(in) :: psi
       real(wp), dimension(ndimpot) :: pot
       real(gp), intent(out) :: ekin_sum,epot_sum
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(out) :: hpsi
       type(linearParameters),intent(in):: lin
       type(atoms_data),intent(in):: at
       real(8),dimension(3,at%nat),intent(in):: rxyz
       integer,dimension(orbs%norbp),intent(in):: onWhichAtomp
       logical,intent(in):: withConfinement
     end subroutine local_hamiltonian_LinearConfinement


     subroutine apply_potentialConfinement2(iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, &
            rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, &
            ibyyzz_r) !optional
       use module_base
       implicit none
       integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot, confPotOrder, offsetx, offsety, offsetz
       real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
       real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
            -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
       integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
       real(gp), intent(out) :: epot
       real(8),dimension(3),intent(in):: rxyzConfinement
       real(8),intent(in):: hxh, hyh, hzh, potentialPrefac
     end subroutine apply_potentialConfinement2


     !!subroutine applyprojector(ncplx,l,i,psppar,npspcode,&
     !!     nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
     !!     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj,psi,hpsi,eproj)
     !!  use module_base
     !!  implicit none
     !!  integer, intent(in) :: i,l,npspcode,ncplx
     !!  integer, intent(in) :: nvctr_c,nvctr_f,nseg_c,nseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
     !!  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     !!  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     !!  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keyv_p
     !!  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keyg_p
     !!  real(wp), dimension(*), intent(in) :: proj
     !!  real(gp), dimension(0:4,0:6), intent(in) :: psppar
     !!  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(in) :: psi
     !!  real(gp), intent(inout) :: eproj
     !!  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(inout) :: hpsi
     !!end subroutine applyprojector


     subroutine initializeInguessParameters(iproc, orbs, orbsig, newComm, ip)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc
       type(orbitals_data),intent(in):: orbs, orbsig
       integer,intent(in):: newComm
       type(inguessParameters),intent(inout):: ip
     end subroutine initializeInguessParameters

     subroutine updatePotential(iproc,nproc,geocode,ixc,nspin,hxh,hyh,hzh,Glr,denspot,ehart,eexcu,vexcu)
       use module_base
       use module_types
       implicit none
       integer:: iproc,nproc,ixc,nspin
       real(gp), intent(in) :: hxh,hyh,hzh
       character(len=1), intent(in) :: geocode
       type(locreg_descriptors), intent(in) :: Glr
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
     
     subroutine countOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(overlapParameters),intent(out):: op
       type(p2pComms),intent(out):: comon
     end subroutine countOverlaps
     
     subroutine determineOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(in):: lzd
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
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
     
     subroutine initCommsOrtho(iproc, nproc, lzd, orbs, onWhichAtomAll, input, locregShape, op, comon, tag)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
       type(input_variables),intent(in):: input
       character(len=1),intent(in):: locregShape
       type(overlapParameters),intent(out):: op
       type(p2pComms),intent(out):: comon
       integer,intent(inout):: tag
     end subroutine initCommsOrtho
     
     subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
       integer,dimension(8),intent(out):: comarr
     end subroutine setCommsParameters
     
     
     subroutine postCommsOverlap(iproc, nproc, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(p2pComms),intent(inout):: comon
     end subroutine postCommsOverlap
     
     
     subroutine extractOrbital(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, comon)
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
     end subroutine extractOrbital
     
     subroutine gatherOrbitals(iproc, nproc, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(p2pComms),intent(inout):: comon
     end subroutine gatherOrbitals
     
     subroutine calculateOverlapMatrix(iproc, nproc, orbs, op, comon, onWhichAtom, lovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(inout):: comon
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(out):: lovrlp
     end subroutine calculateOverlapMatrix
     
     
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
     
     subroutine transformOverlapMatrix(iproc, nproc, comm, blocksize_dsyev, blocksize_pdgemm, norb, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, comm, blocksize_dsyev, blocksize_pdgemm, norb
       real(8),dimension(norb,norb),intent(inout):: ovrlp
     end subroutine transformOverlapMatrix
     
     subroutine expandOrbital(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(input_variables),intent(in):: input
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(in):: comon
       real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
     end subroutine expandOrbital
     
     subroutine localGramschmidt(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, comon, lovrlp, lphiovrlp, lphi)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs, lorbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(in):: comon
       real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(in):: lovrlp
       real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
       real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: lphi
     end subroutine localGramschmidt
     
     
     subroutine globalLoewdin(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, comon, ovrlp, lphiovrlp, lphi)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs, lorbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(in):: comon
       real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
       real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
       real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out):: lphi
     end subroutine globalLoewdin


     subroutine orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, &
                blocksize_pdgemm, orbs, op, comon, lzd, mad, lphi, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, blocksize_pdgemm
       type(orbitals_data),intent(in):: orbs
       type(overlapParameters),intent(inout):: op
       type(p2pComms),intent(inout):: comon
       type(local_zone_descriptors),intent(in):: lzd
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
       real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
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

     subroutine initializeRepartitionOrbitals(iproc, nproc, tag, lin)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       integer,intent(inout):: tag
       type(linearParameters),intent(inout):: lin
     end subroutine initializeRepartitionOrbitals

     subroutine postCommsRepartition(iproc, nproc, orbs, comrp, nsendBuf, sendBuf, nrecvBuf, recvBuf)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
       type(orbitals_data),intent(in):: orbs
       type(p2pCommsRepartition),intent(inout):: comrp
       real(8),dimension(nsendBuf),intent(in):: sendBuf
       real(8),dimension(nrecvBuf),intent(out):: recvBuf
     end subroutine postCommsRepartition


     subroutine gatherDerivativeOrbitals(iproc, nproc, orbs, comrp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(p2pCommsRepartition),intent(inout):: comrp
     end subroutine gatherDerivativeOrbitals


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



     subroutine determineLocalizationRegions(iproc, nproc, nlr, norb, at, onWhichAtomALl, locrad, rxyz, lzd, hx, hy, hz, mlr)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norb
       type(atoms_data),intent(in):: at
       integer,dimension(norb),intent(in):: onWhichAtomAll
       real(8),dimension(at%nat),intent(in):: locrad
       real(8),dimension(3,at%nat),intent(in):: rxyz
       type(local_zone_descriptors),intent(in):: lzd
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


     subroutine orthonormalizeVectors(iproc, nproc, comm, nItOrtho, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm, &
                orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mad, mlr, vec, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, comm, nItOrtho, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm
       integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       type(matrixDescriptors),intent(in):: mad
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       real(8),dimension(norbmax,norbp),intent(inout):: vec
       type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
     end subroutine orthonormalizeVectors




     subroutine initCommsMatrixOrtho(iproc, nproc, norb, norb_par, isorb_par, onWhichAtomPhi, onWhichMPI, tag, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, norb
       integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
       integer,dimension(nproc),intent(in):: norb_par, isorb_par
       integer,intent(inout):: tag
       type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
     end subroutine initCommsMatrixOrtho



     subroutine determineOverlapRegionMatrix(iproc, nproc, lzd, mlr, orbs, orbstot, onWhichAtom, onWhichAtomPhi, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzd
       type(orbitals_data),intent(in):: orbs, orbstot
       integer,dimension(orbstot%norb),intent(in):: onWhichAtom
       integer,dimension(orbs%norb),intent(in):: onWhichAtomPhi
       type(matrixLocalizationRegion),dimension(lzd%nlr),intent(in):: mlr
       type(p2pCommsOrthonormalityMatrix),intent(out):: comom
     end subroutine determineOverlapRegionMatrix


     subroutine postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
     use module_base
     use module_types
     implicit none
     integer,intent(in):: iproc, nproc, newComm
     type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
     end subroutine postCommsVectorOrthonormalization


     subroutine gatherVectors(iproc, nproc, newComm, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, newComm
       type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
     end subroutine gatherVectors


     subroutine extractToOverlapregion(iproc, nproc, norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, norbmax, norb, norbp
       integer,dimension(norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       real(8),dimension(norbmax,norbp),intent(in):: vec
       type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
     end subroutine extractToOverlapregion



     subroutine expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, isorb, norbp, norbmax, noverlaps
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(p2pCommsOrthonormalityMatrix),intent(in):: comom
       real(8),dimension(norbmax,noverlaps),intent(out):: vecOvrlp
     end subroutine expandFromOverlapregion

     subroutine calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr,&
                onWhichAtom, vec, vecOvrlp, newComm, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, newComm
       type(p2pCommsOrthonormalityMatrix),intent(in):: comom
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       integer,dimension(norb),intent(in):: onWhichAtom
       real(8),dimension(norbmax,norbp),intent(in):: vec
       real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
       real(8),dimension(norb,norb),intent(out):: ovrlp
     end subroutine calculateOverlap


     subroutine orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom,&
               vecOvrlp, ovrlp, vec)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb
       type(p2pCommsOrthonormalityMatrix),intent(in):: comom
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       integer,dimension(norb),intent(in):: onWhichAtom
       real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
       real(8),dimension(norb,norb),intent(in):: ovrlp
       real(8),dimension(norbmax,norbp),intent(inout):: vec
     end subroutine orthonormalLinearCombinations


     subroutine buildLinearCombinationsLocalized(iproc, nproc, orbsig, orbs, comms, at, Glr, input, norbsPerType, &
           onWhichAtom, lchi, lphi, rxyz, onWhichAtomPhi, lin, lzdig, ham)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbsig, orbs
       type(communications_arrays),intent(in):: comms
       type(atoms_data),intent(in):: at
       type(locreg_descriptors),intent(in):: Glr
       type(input_variables),intent(in):: input
       type(linearParameters),intent(in):: lin
       type(local_zone_descriptors),intent(inout):: lzdig
       integer,dimension(at%ntypes):: norbsPerType
       integer,dimension(orbsig%norb),intent(in):: onWhichAtom
       real(8),dimension(max(orbsig%npsidim_comp,orbsig%npsidim_orbs)):: lchi
       real(8),dimension(max(lin%orbs%npsidim_comp,lin%orbs%npsidim_orbs)):: lphi
       real(8),dimension(3,at%nat):: rxyz
       integer,dimension(orbs%norb):: onWhichAtomPhi
       real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(inout):: ham
     end subroutine buildLinearCombinationsLocalized


     subroutine orthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm, &
                orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, mad, vec, grad, comom, trace)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm
       integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
       type(orbitals_data),intent(in):: orbs
       integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
       integer,dimension(0:nproc-1),intent(in):: isorb_par
       type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(norbmax,norbp),intent(inout):: vec, grad
       type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
       real(8),intent(out):: trace
     end subroutine orthoconstraintVectors


     subroutine cubic_exact_exchange(iproc,nproc,nspin,npsidim,size_potxc,hx,hy,hz,Glr,orbs,&
                ngatherarr,psi,potxc,eexctX,pkernel,orbsocc,psirocc)
       use module_base
       use module_types
       use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,npsidim,size_potxc
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors) :: Glr
       type(orbitals_data) :: orbs
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(wp), dimension(npsidim), intent(in) :: psi
       real(wp), dimension(size_potxc),intent(out) :: potxc
       real(gp), intent(out) :: eexctX
       real(dp), dimension(*), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
     end subroutine


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
     type(p2pCommsRepartition),intent(inout):: comrp
     real(8),dimension(nphi),intent(in):: phi
     real(8),dimension(max(lborbs%npsidim_orbs,lborbs%npsidim_comp)),target,intent(out):: phid
     end subroutine getDerivativeBasisFunctions


     subroutine buildLinearCombinations(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, locregShape, &
                tag, comonig, opig, madig, lphi)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(local_zone_descriptors),intent(in):: lzdig, lzd
       type(orbitals_data),intent(in):: orbsig, orbs
       type(input_variables),intent(in):: input
       real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
       real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
       character(len=1),intent(in):: locregShape
       integer,intent(inout):: tag
       real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: lphi
       type(p2pComms):: comonig
       type(overlapParameters):: opig
       type(matrixDescriptors):: madig
     end subroutine buildLinearCombinations


     subroutine postCommunicationsPotential(iproc, nproc, ndimpot, pot, comgp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, ndimpot
       real(8),dimension(ndimpot),intent(in):: pot
       !type(p2pCommsGatherPot),intent(inout):: comgp
       type(p2pComms),intent(inout):: comgp
     end subroutine postCommunicationsPotential


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


     subroutine gatherOrbitals2(iproc, nproc, comon)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(p2pComms),intent(inout):: comon
     end subroutine gatherOrbitals2


     subroutine expandOrbital2(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
       use module_base
       use module_types
       implicit none
       
       ! Calling arguments
       integer,intent(in):: iproc, nproc
       type(orbitals_data),intent(in):: orbs
       type(input_variables),intent(in):: input
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       type(local_zone_descriptors),intent(in):: lzd
       type(overlapParameters),intent(in):: op
       type(p2pComms),intent(in):: comon
       real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
     end subroutine expandOrbital2

     subroutine getOverlapMatrix(iproc, nproc, lin, input, lphi, mad, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(linearParameters),intent(inout):: lin
       type(input_variables),intent(in):: input
       real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
       type(matrixDescriptors),intent(in):: mad
       real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(out):: ovrlp
     end subroutine getOverlapMatrix


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


     subroutine allocateCommunicationbufferSumrho(iproc, with_auxarray, comsr, subname)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc
       logical,intent(in):: with_auxarray
       !type(p2pCommsSumrho),intent(inout):: comsr
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

    subroutine nullify_linearParameters(lin)
      use module_base
      use module_types
      implicit none
      type(linearParameters),intent(out):: lin
    end subroutine nullify_linearParameters

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

    subroutine nullify_largeBasis(lb)
      use module_base
      use module_types
      implicit none
      type(largeBasis),intent(out):: lb
    end subroutine nullify_largeBasis

    subroutine nullify_p2pCommsRepartition(comrp)
      use module_base
      use module_types
      implicit none
      type(p2pCommsRepartition),intent(out):: comrp
    end subroutine nullify_p2pCommsRepartition

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

    subroutine nullify_linearInputGuess(lig)
      use module_base
      use module_types
      implicit none
      type(linearInputGuess),intent(out):: lig
    end subroutine nullify_linearInputGuess

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

    subroutine initLocregs(iproc, nproc, nat, rxyz, lzd, orbs, input, Glr, locrad, locregShape, lborbs)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, nat
      real(8),dimension(3,nat),intent(in):: rxyz
      type(local_zone_descriptors),intent(inout):: lzd
      type(orbitals_data),intent(inout):: orbs
      type(input_variables),intent(in):: input
      type(locreg_descriptors),intent(in):: Glr
      real(8),dimension(lzd%nlr),intent(in):: locrad
      character(len=1),intent(in):: locregShape
      type(orbitals_data),optional,intent(inout):: lborbs
    end subroutine initLocregs

    subroutine deallocate_linearParameters(lin, subname)
      use module_base
      use module_types
      implicit none
      ! Calling arguments
      type(linearParameters),intent(inout):: lin
      character(len=*),intent(in):: subname
    end subroutine deallocate_linearParameters

    subroutine deallocate_largeBasis(lb, subname)
      use module_base
      use module_types
      implicit none
      type(largeBasis),intent(inout):: lb
      character(len=*),intent(in):: subname
    end subroutine deallocate_largeBasis
    
    subroutine dealloctae_p2pCommsRepartition(comrp, subname)
      use module_base
      use module_types
      implicit none
      type(p2pCommsRepartition),intent(inout):: comrp
      character(len=*),intent(in):: subname
    end subroutine dealloctae_p2pCommsRepartition
    
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

    subroutine deallocate_inguessParameters(ip, subname)
      use module_base
      use module_types
      implicit none
      type(inguessParameters),intent(inout):: ip
      character(len=*),intent(in):: subname
    end subroutine deallocate_inguessParameters

    subroutine deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
      use module_base
      use module_types
      implicit none
      type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
      character(len=*),intent(in):: subname
    end subroutine deallocate_p2pCommsOrthonormalityMatrix

    subroutine deallocate_matrixDescriptors(mad, subname)
      use module_base
      use module_types
      implicit none
      type(matrixDescriptors),intent(inout):: mad
      character(len=*),intent(in):: subname
    end subroutine deallocate_matrixDescriptors


    subroutine cancelCommunicationPotential(iproc, nproc, comgp)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      !type(p2pCommsGatherPot),intent(inout):: comgp
      type(p2pComms),intent(inout):: comgp
    end subroutine cancelCommunicationPotential

    subroutine initCommsOrthoVariable(iproc, nproc, lzd, orbs, orbsig, onWhichAtomAll, input, op, comon, tag)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs, orbsig
      integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
      type(input_variables),intent(in):: input
      type(overlapParameters),intent(out):: op
      type(p2pComms),intent(out):: comon
      integer,intent(inout):: tag
    end subroutine initCommsOrthoVariable
    
    subroutine countOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(out):: op
      type(p2pComms),intent(out):: comon
    end subroutine countOverlapsVariable
    
    subroutine determineOverlapsVariable(iproc, nproc, orbs, orbsig, lzd, op, comon)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(out):: op
      type(p2pComms),intent(out):: comon
    end subroutine determineOverlapsVariable
    
    subroutine determineOverlapDescriptorsVariable(iproc, nproc, orbs, orbsig, lzd, Glr, onWhichAtom, op)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(locreg_descriptors),intent(in):: Glr
      integer,dimension(orbs%norb),intent(in):: onWhichAtom
      type(overlapParameters),intent(inout):: op
    end subroutine determineOverlapDescriptorsVariable
    
    subroutine setCommsOrthoVariable(iproc, nproc, orbs, orbsig, lzd, op, comon, tag)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(inout):: op
      type(p2pComms),intent(out):: comon
      integer,intent(inout):: tag
    end subroutine setCommsOrthoVariable
    
    subroutine indicesForExpansionVariable(iproc, nproc, orbs, input, lzd, op, comon)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs
      type(input_variables),intent(in):: input
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(in):: op
      type(p2pComms),intent(in):: comon
    end subroutine indicesForExpansionVariable
    
    subroutine indicesForExtractionVariable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, comon)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, sizePhi
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(inout):: op
      type(p2pComms),intent(out):: comon
    end subroutine indicesForExtractionVariable
    
    subroutine extractOrbital2Variable(iproc, nproc, orbs, orbsig, sizePhi, lzd, op, phi, comon)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, sizePhi
      type(orbitals_data),intent(in):: orbs, orbsig
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(inout):: op
      real(8),dimension(sizePhi),intent(in):: phi
      type(p2pComms),intent(out):: comon
    end subroutine extractOrbital2Variable
    
    subroutine expandOrbital2Variable(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(orbitals_data),intent(in):: orbs
      type(input_variables),intent(in):: input
      type(local_zone_descriptors),intent(in):: lzd
      type(overlapParameters),intent(in):: op
      type(p2pComms),intent(in):: comon
      real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
    end subroutine expandOrbital2Variable

    subroutine buildLinearCombinationsVariable(iproc, nproc, lzdig, lzd, orbsig, &
               orbs, input, coeff, lchi, tag, lphi)
      use module_base
      use module_types
      implicit none

      ! Calling arguments
      integer,intent(in):: iproc, nproc
      type(local_zone_descriptors),intent(in):: lzdig, lzd
      type(orbitals_data),intent(in):: orbsig, orbs
      type(input_variables),intent(in):: input
      real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
      real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
      integer,intent(inout):: tag
      real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
    end subroutine buildLinearCombinationsVariable

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

    subroutine initInputguessConfinement(iproc, nproc, at, Glr, input, lin, lig, rxyz, nscatterarr, tag)
      use module_base
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: iproc,nproc
      type(atoms_data), intent(inout) :: at
      type(locreg_descriptors), intent(in) :: Glr
      type(input_variables):: input
      type(linearParameters),intent(inout):: lin
      type(linearInputGuess),intent(inout):: lig
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      integer,intent(inout):: tag
    end subroutine initInputguessConfinement

    subroutine orthonormalizeAtomicOrbitalsLocalized(iproc, nproc, lzd, orbs, input, lchi)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs
      type(input_variables),intent(in):: input
      real(8),dimension(orbs%npsidim_comp,orbs%npsidim_orbs),intent(inout):: lchi
    end subroutine orthonormalizeAtomicOrbitalsLocalized


    subroutine orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, &
               blocksize_pdgemm, lzd, orbs, comon, op, input, mad, lchi)
      use module_base
      use module_types
      implicit none
      ! Calling arguments
      integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, blocksize_pdgemm
      type(local_zone_descriptors),intent(in):: lzd
      type(orbitals_data),intent(in):: orbs
      type(input_variables),intent(in):: input
      type(p2pComms),intent(inout):: comon
      type(overlapParameters),intent(inout):: op
      type(matrixDescriptors),intent(in):: mad
      real(8),dimension(orbs%npsidim_comp),intent(inout):: lchi
    end subroutine orthonormalizeAtomicOrbitalsLocalized2

    subroutine buildLinearCombinationsLocalized3(iproc, nproc, orbsig, orbs, comms, at, Glr, input, norbsPerType, &
      onWhichAtom, lchi, lphi, rxyz, onWhichAtomPhi, lin, lzdig, nlocregPerMPI, tag, ham3, comonig, opig, madig)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, nlocregPerMPI
      type(orbitals_data),intent(in):: orbsig, orbs
      type(communications_arrays),intent(in):: comms
      type(atoms_data),intent(in):: at
      type(locreg_descriptors),intent(in):: Glr
      type(input_variables),intent(in):: input
      type(linearParameters),intent(in):: lin
      type(local_zone_descriptors),intent(inout):: lzdig
      integer,dimension(at%ntypes):: norbsPerType
      integer,dimension(orbsig%norb),intent(in):: onWhichAtom
      real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp)):: lchi
      real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)):: lphi
      real(8),dimension(3,at%nat):: rxyz
      integer,dimension(orbs%norb):: onWhichAtomPhi
      integer,intent(inout):: tag
      real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(inout):: ham3
      type(p2pComms):: comonig
      type(overlapParameters):: opig
      type(matrixDescriptors):: madig
    end subroutine buildLinearCombinationsLocalized3

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
                 orbs, lorbs, onWhichAtom, lzd, op, comon, lagmat, ovrlp, lphiovrlp, mad, lhphi)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
        type(orbitals_data),intent(in):: orbs, lorbs
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(in):: op
        type(p2pComms),intent(in):: comon
        real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
        real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
        real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out):: lhphi
      end subroutine applyOrthoconstraintNonorthogonal2

      subroutine gatherOrbitalsOverlapWithComput(iproc, nproc, orbs, input, lzd, op, comon, lphiovrlp, expanded)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(in):: op
        type(p2pComms),intent(inout):: comon
        real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
        logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
      end subroutine gatherOrbitalsOverlapWithComput


      subroutine expandOneOrbital(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, orbsource, orbdest
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(in):: op
        type(p2pComms),intent(in):: comon
        real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      end subroutine expandOneOrbital

      subroutine expandRemainingOrbitals(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, expanded, lphiovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(in):: op
        type(p2pComms),intent(in):: comon
        logical,dimension(orbs%norb,orbs%norbp),intent(in):: expanded
        real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      end subroutine expandRemainingOrbitals

      subroutine extractOrbital3(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, nsendBuf, sendBuf)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, sizePhi
        type(orbitals_data),intent(in):: orbs
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(inout):: op
        real(8),dimension(sizePhi),intent(in):: phi
        integer,intent(in):: nsendBuf
        real(8),dimension(nsendBuf),intent(out):: sendBuf
      end subroutine extractOrbital3

      subroutine calculateOverlapMatrix3(iproc, nproc, orbs, op, onWhichAtom, nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        real(8),dimension(nsendBuf),intent(in):: sendBuf
        real(8),dimension(nrecvBuf),intent(in):: recvBuf
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
      end subroutine calculateOverlapMatrix3


      subroutine calculateOverlapMatrix3Partial(iproc, nproc, orbs, op, onWhichAtom, &
                 nsendBuf, sendBuf, nrecvBuf, recvBuf, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
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
                 comm, norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, vecOvrlp, ovrlp, &
                 lagmat, comom, mlr, mad, orbs, grad)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOverlap, blocksize_pdgemm
        integer,intent(in):: comm, norb, norbmax, norbp, isorb, nlr, noverlaps
        integer,dimension(norb),intent(in):: onWhichAtom
        real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
        real(8),dimension(norb,norb),intent(in):: ovrlp
        real(8),dimension(norb,norb),intent(inout):: lagmat
        type(p2pCommsOrthonormalityMatrix),intent(in):: comom
        type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
        type(matrixDescriptors),intent(in):: mad
        type(orbitals_data),intent(in):: orbs
        real(8),dimension(norbmax,norbp),intent(inout):: grad
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

      subroutine orthoconstraintNonorthogonal(iproc, nproc, lzd, orbs, op, comon, mad, ovrlp, &
                 methTransformOverlap, blocksize_pdgemm, lphi, lhphi, lagmat)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
        type(local_zone_descriptors),intent(in):: lzd
        type(orbitals_Data),intent(in):: orbs
        type(overlapParameters),intent(inout):: op
        type(p2pComms),intent(inout):: comon
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
        real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(in):: lphi
        real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lhphi
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: lagmat
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

      subroutine getOrbitals(iproc, nproc, comon)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(p2pComms),intent(inout):: comon
      end subroutine getOrbitals

      subroutine initCompressedMatmul(iproc, nproc, norb, mad)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, norb
        type(matrixDescriptors),intent(inout):: mad
      end subroutine initCompressedMatmul

      subroutine initCompressedMatmul2(norb, nseg, keyg, nsegmatmul, keygmatmul, keyvmatmul)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: norb, nseg
        integer,dimension(2,nseg),intent(in):: keyg
        integer,intent(out):: nsegmatmul
        integer,dimension(:,:),pointer,intent(out):: keygmatmul
        integer,dimension(:),pointer,intent(out):: keyvmatmul
      end subroutine initCompressedMatmul2


      subroutine dgemm_compressed2(iproc, nproc, norb, nsegline, nseglinemax, keygline, nsegmatmul, keygmatmul, a, b, c)
        implicit none
        integer,intent(in):: iproc, nproc, norb, nseglinemax, nsegmatmul
        integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
        integer,dimension(norb):: nsegline
        !integer,dimension(2,maxval(nsegline),norb):: keygline
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

      subroutine overlapPowerMinusOneHalfTaylor(iproc, nproc, methTransformOrder, norb, mad, ovrlp)
        use module_base
        use module_types
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, nproc, methTransformOrder, norb
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine overlapPowerMinusOneHalfTaylor

      subroutine overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOrder, blocksize_dsyev, &
                 blocksize_pdgemm, norb, mad, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, comm, methTransformOrder, blocksize_dsyev, blocksize_pdgemm, norb
        type(matrixDescriptors),intent(in):: mad
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine overlapPowerMinusOneHalf

      subroutine overlapPowerMinusOne(iproc, nproc, iorder, norb, mad, orbs, ovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, iorder, norb
        type(matrixDescriptors),intent(in):: mad
        type(orbitals_data),intent(in):: orbs
        real(8),dimension(norb,norb),intent(inout):: ovrlp
      end subroutine overlapPowerMinusOne

      subroutine initCompressedMatmul3(norb, mad)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: norb
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

     subroutine local_hamiltonian3(iproc,exctX,orbs,Lzd,hx,hy,hz,&
          nspin,Lpot,psi,hpsi,ekin_sum,epot_sum,&
          withConfinement, at, rxyz, istexct, lin, confinementCenter)
       use module_base
       use module_types
       use libxc_functionals
       implicit none
       integer, intent(in) :: iproc,nspin, istexct
       real(gp), intent(in) :: hx,hy,hz
       logical, intent(in) :: exctX
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors), intent(in) :: Lzd
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       real(wp), dimension(Lzd%ndimpotisf),target :: Lpot
       !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
       real(gp), intent(out) :: ekin_sum,epot_sum
       real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
       logical,intent(in):: withConfinement
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       type(linearParameters),intent(in),optional:: lin
       integer,dimension(orbs%norbp),intent(in),optional:: confinementCenter
     end subroutine local_hamiltonian3

     subroutine HamiltonianApplication3(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
          proj,Lzd,nlpspd,confdatarr,ngatherarr,Lpot,psi,hpsi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,SIC,GPU,&
       !withConfinement,energyReductionFlag, &
          pkernel,orbsocc,psirocc)
       !,lin, confinementCenter)
       use module_base
       use module_types
       use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc!,nspin
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors),intent(in) :: Lzd
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(SIC_data), intent(in) :: SIC
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(Lzd%Lnprojel), intent(in) :: proj
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
       real(wp), dimension(lzd%ndimpotisf) :: Lpot
       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC
       real(wp), target, dimension(orbs%npsidim_orbs), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       !real(dp), dimension(*), optional :: pkernel
       real(dp), dimension(:), pointer, optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
     end subroutine HamiltonianApplication3

     subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
          proj,Lzd,nlpspd,confdatarr,ngatherarr,Lpot,psi,hpsi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,SIC,GPU,&
          proj_G,paw,&
          pkernel,orbsocc,psirocc)
       use module_base
       use module_types
       use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc!,nspin
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors),intent(in) :: Lzd
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(SIC_data), intent(in) :: SIC
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(Lzd%Lnprojel), intent(inout) :: proj
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
       real(wp), dimension(lzd%ndimpotisf) :: Lpot
       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       real(dp), dimension(:), pointer, optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
       !PAW variables:
       type(gaussian_basis),dimension(at%ntypes),intent(in)::proj_G
       type(paw_objects),intent(inout)::paw
     end subroutine FullHamiltonianApplication


!!$     subroutine HamiltonianApplication3(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
!!$          proj,Lzd,ngatherarr,Lpot,psi,hpsi,&
!!$          ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,withConfinement,energyReductionFlag,&
!!$          pkernel,orbsocc,psirocc,lin,confinementCenter)
!!$       use module_base
!!$       use module_types
!!$       use libxc_functionals
!!$       implicit none
!!$       integer, intent(in) :: iproc,nproc,nspin
!!$       real(gp), intent(in) :: hx,hy,hz
!!$       type(atoms_data), intent(in) :: at
!!$       type(orbitals_data), intent(in) :: orbs
!!$       type(local_zone_descriptors),intent(in) :: Lzd
!!$       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!$       real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!$       real(wp), dimension(Lzd%Lnprojel), intent(in) :: proj
!!$       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
!!$       real(wp), dimension(lzd%ndimpotisf) :: Lpot
!!$       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
!!$       real(wp), target, dimension(orbs%npsidim_orbs), intent(out) :: hpsi
!!$       type(GPU_pointers), intent(inout) :: GPU
!!$       logical,intent(in):: withConfinement
!!$       logical,intent(in):: energyReductionFlag
!!$       real(dp), dimension(*), optional :: pkernel
!!$       type(orbitals_data), intent(in), optional :: orbsocc
!!$       real(wp), dimension(:), pointer, optional :: psirocc
!!$       type(linearParameters),intent(in),optional:: lin
!!$       integer,dimension(orbs%norbp),intent(in),optional:: confinementCenter
!!$     end subroutine HamiltonianApplication3

     !!subroutine full_local_potential2(iproc,nproc,ndimpot,ndimgrid,ndimrhopot,nspin,orbs,lzd,ngatherarr,potential,Lpot,flag,comgp)
     !!  use module_base
     !!  use module_types
     !!  use module_xc
     !!  implicit none
     !!  integer, intent(in) :: iproc,nproc,ndimpot,ndimgrid,flag,nspin,ndimrhopot
     !!  type(orbitals_data),intent(inout):: orbs
     !!  type(local_zone_descriptors),intent(inout):: lzd
     !!  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
     !!  real(wp), dimension(max(ndimpot,1)*orbs%nspin), intent(in), target ::potential
     !!  real(wp), dimension(:), pointer, intent(out) :: Lpot
     !!  type(p2pCommsGatherPot),intent(inout), optional:: comgp
     !!end subroutine full_local_potential2

     subroutine prepare_lnlpspd(iproc, at, input, orbs, rxyz, radii_cf, locregShape, lzd)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc
       type(atoms_data),intent(in):: at
       type(input_variables),intent(in):: input
       type(orbitals_data),intent(in):: orbs
       real(8),dimension(3,at%nat),intent(in):: rxyz
       real(8),dimension(at%ntypes,3),intent(in):: radii_cf
       character(len=1),intent(in):: locregShape
       type(local_zone_descriptors),intent(inout):: lzd
     end subroutine prepare_lnlpspd

     subroutine free_lnlpspd(orbs, lzd)
       use module_base
       use module_types
       implicit none
       type(orbitals_data),intent(in):: orbs
       type(local_zone_descriptors),intent(inout):: lzd
     end subroutine free_lnlpspd


     subroutine transformToGlobal(iproc, nproc, lin, orbs, comms, input, coeff, lphi, psi, psit)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc
       type(linearParameters),intent(in):: lin
       type(orbitals_data),intent(in):: orbs
       type(communications_arrays):: comms
       type(input_variables),intent(in):: input
       real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
       real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
       real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: psi, psit
     end subroutine transformToGlobal

     subroutine gatherOrbitalsOverlapWithComput2(iproc, nproc, orbs, input, lzd, op, comon, nsendbuf, sendbuf,&
          nrecvbuf, recvbuf, lphiovrlp, expanded, ovrlp)
       use module_base
       use module_types
       implicit none
       integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf                                                                 
       type(orbitals_data),intent(in):: orbs
       type(input_variables),intent(in):: input
       type(local_zone_descriptors),intent(in):: lzd                                                                         
       type(overlapParameters),intent(in):: op                                                                               
       type(p2pComms),intent(inout):: comon
       real(8),dimension(nsendbuf),intent(in):: sendbuf                                                                      
       real(8),dimension(nrecvbuf),intent(in):: recvbuf                                                                      
       real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
       logical,dimension(orbs%norb,orbs%norbp),intent(out):: expanded
       real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
     end subroutine gatherOrbitalsOverlapWithComput2


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

      subroutine collectAndCalculateOverlap(iproc, nproc, comon, mad, op, orbs, input, &
                 lzd, nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, lphiovrlp, timecommunp2p,&
                 timecommuncoll, timeoverlap, timeexpand, timecompress)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
        type(p2pComms),intent(inout):: comon
        type(matrixDescriptors),intent(in):: mad
        type(overlapParameters),intent(in):: op
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        type(local_zone_descriptors),intent(in):: lzd
        real(8),dimension(nsendbuf),intent(in):: sendbuf
        real(8),dimension(nrecvbuf),intent(inout):: recvbuf
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
        real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
        real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timeexpand, timecompress
      end subroutine collectAndCalculateOverlap


      subroutine postCommsOverlapNew(iproc, nproc, orbs, op, lzd, phi, comon, timecommun, timeextract)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(orbitals_data),intent(in):: orbs
        type(overlapParameters),intent(in):: op
        type(local_zone_descriptors),intent(in):: lzd
        real(8),dimension(orbs%npsidim_orbs),intent(in):: phi
        type(p2pComms),intent(inout):: comon
        real(8),intent(out):: timecommun, timeextract
      end subroutine postCommsOverlapNew

      subroutine expandOneOrbital2(iproc, nproc, orbsource, orbdest, orbs, input, onWhichAtom, lzd, op, &
           nrecvbuf, recvbuf, lphiovrlp)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, orbsource, orbdest, nrecvbuf
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(in):: op
        real(8),dimension(nrecvbuf),intent(in):: recvbuf
        real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp
      end subroutine expandOneOrbital2


      subroutine calculateForcesLinear(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, atoms, in, comms, lin, nlpspd, proj, &
                 ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, rho, psi, fxyz, fnoise)
        use module_base
        use module_types
        implicit none
        
        ! Calling arguments
        integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
        type(locreg_descriptors),intent(in):: Glr
        type(orbitals_data),intent(in):: orbs
        type(atoms_data),intent(in):: atoms
        type(input_variables),intent(in):: in
        type(communications_arrays),intent(in):: comms
        type(linearParameters),intent(inout):: lin
        type(nonlocal_psp_descriptors),intent(in) :: nlpspd
        real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
        integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr   !!! NOT NEEDED
        integer,dimension(0:nproc-1,4),intent(inout) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
        type(GPU_pointers),intent(inout):: GPU
        integer,dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
        real(dp),dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
        real(dp),dimension(lin%as%size_pkernel),intent(in):: pkernel
        real(8),dimension(3,atoms%nat),intent(in):: rxyz, fion, fdisp
        real(8),dimension(3,atoms%nat),intent(out):: fxyz
        real(8),intent(out):: fnoise
        real(8),dimension(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1)),intent(in):: rho
        real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
      end subroutine calculateForcesLinear

      subroutine collectAndCalculateOverlap2(iproc, nproc, comon, mad, op, orbs, input, lzd, &
                 nsendbuf, sendbuf, nrecvbuf, recvbuf, ovrlp, timecommunp2p, timecommuncoll, timeoverlap, timecompress)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
        type(p2pComms),intent(inout):: comon
        type(matrixDescriptors),intent(in):: mad
        type(overlapParameters),intent(in):: op
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        type(local_zone_descriptors),intent(in):: lzd
        real(8),dimension(nsendbuf),intent(in):: sendbuf
        real(8),dimension(nrecvbuf),intent(inout):: recvbuf
        real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp
        real(8),intent(inout):: timecommunp2p, timecommuncoll, timeoverlap, timecompress
      end subroutine collectAndCalculateOverlap2

       subroutine applyOrthoconstraintNonorthogonalCubic(iproc, nproc, methTransformOverlap, blocksize_pdgemm, &
                  orbs, gorbs, comms, lzd, input, &
                  op, ovrlp, mad, lphi, lhphi, trH)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, methTransformOverlap, blocksize_pdgemm
         type(orbitals_data),intent(in):: orbs, gorbs
         type(communications_arrays),intent(in):: comms
         type(local_zone_descriptors),intent(in):: lzd
         type(input_variables),intent(in):: input
         type(overlapParameters),intent(in):: op
         real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
         type(matrixDescriptors),intent(in):: mad
         real(8),dimension(max(orbs%npsidim_comp,orbs%npsidim_orbs)),intent(inout):: lphi, lhphi
         real(8),intent(out):: trH
       end subroutine applyOrthoconstraintNonorthogonalCubic

       subroutine collectnew(iproc, nproc, comon, mad, op, orbs, lzd, &
                  nsendbuf, sendbuf, nrecvbuf, recvbuf, timecommunp2p, &
                  timecommuncoll, timecompress)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
         type(p2pComms),intent(inout):: comon
         type(matrixDescriptors),intent(in):: mad
         type(overlapParameters),intent(in):: op
         type(orbitals_data),intent(in):: orbs
         type(local_zone_descriptors),intent(in):: lzd
         real(8),dimension(nsendbuf),intent(in):: sendbuf
         real(8),dimension(nrecvbuf),intent(inout):: recvbuf
         real(8),intent(inout):: timecommunp2p, timecommuncoll, timecompress
       end subroutine collectnew

       subroutine my_iallgather_collect2(iproc, nproc, sendcount, recvcounts, requests)
         use module_base
         implicit none
         integer,intent(in):: iproc, nproc, sendcount
         integer,dimension(0:nproc-1),intent(in):: recvcounts
         integer,dimension(2,0:nproc-1),intent(inout):: requests
       end subroutine my_iallgather_collect2


       subroutine initMatrixCompression(iproc, nproc, nlr, orbs, noverlaps, overlaps, mad)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         type(orbitals_data),intent(in):: orbs
         !integer,dimension(nlr),intent(in):: noverlaps
         integer,dimension(orbs%norb),intent(in):: noverlaps
         !integer,dimension(maxval(noverlaps(:)),nlr),intent(in):: overlaps
         integer,dimension(maxval(noverlaps(:)),orbs%norb),intent(in):: overlaps
         type(matrixDescriptors),intent(out):: mad
       end subroutine initMatrixCompression


       subroutine postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, newComm
         type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
       end subroutine postCommsVectorOrthonormalizationNew


       subroutine gatherVectorsNew(iproc, nproc, comom)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
       end subroutine gatherVectorsNew


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


       subroutine getHamiltonianMatrix6(iproc, nproc, nprocTemp, lzdig, orbsig, orbs, onWhichMPITemp, &
                  input, onWhichAtom, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad, memoryForCommunOverlapIG, locregShape, &
                  tagout, ham)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nprocTemp, ndim_lhchi, nlocregPerMPI
         type(local_zone_descriptors),intent(in):: lzdig
         type(orbitals_data),intent(in):: orbsig, orbs
         integer,dimension(orbs%norb),intent(in):: onWhichMPITemp
         type(input_variables),intent(in):: input
         integer,dimension(orbsig%norb),intent(in):: onWhichAtom
         real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp)),intent(in):: lchi
         real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp),ndim_lhchi),intent(in):: lhchi
         logical,dimension(lzdig%nlr),intent(in):: skip
         type(matrixDescriptors),intent(in):: mad
         integer,intent(in):: memoryForCommunOverlapIG
         character(len=1),intent(in):: locregShape
         integer,intent(inout):: tagout
         real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham
       end subroutine getHamiltonianMatrix6
       
       subroutine dgemm_compressed_parallel(iproc, nproc, norb, nsegline, nseglinemax, keygline, &
                  nsegmatmul, keygmatmul, norb_par, isorb_par, norbp, a, b, c)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, norb, norbp, nseglinemax, nsegmatmul
         integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
         integer,dimension(norb):: nsegline
         !integer,dimension(2,maxval(nsegline),norb):: keygline
         integer,dimension(2,nseglinemax,norb):: keygline
         integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
         real(8),dimension(norb,norb),intent(in):: a, b
         real(8),dimension(norb,norb),intent(out):: c
       end subroutine dgemm_compressed_parallel


       subroutine getCoefficients_new(iproc, nproc, lin, orbs, hamold, lphi, ovrlp, coeff)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(linearParameters),intent(inout):: lin
         type(orbitals_data),intent(in):: orbs
         real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: hamold
         real(8),dimension(lin%orbs%npsidim_orbs),intent(in):: lphi
         real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(inout):: ovrlp
         real(8),dimension(lin%orbs%norb,orbs%norb),intent(inout):: coeff
       end subroutine getCoefficients_new


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


       subroutine minimize_in_subspace(iproc,nproc,lin,at,input,lpot,GPU,ngatherarr,proj,rxyz,pkernelseq,nlpspd,lphi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc,nproc
         type(linearParameters),intent(inout):: lin
         type(atoms_data),intent(in):: at
         type(input_variables),intent(in):: input
         real(8),dimension(lin%lzd%ndimpotisf),intent(in):: lpot
         type(GPU_pointers),intent(inout):: GPU
         integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
         type(nonlocal_psp_descriptors),intent(in):: nlpspd
         real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
         real(8),dimension(3,at%nat),intent(in):: rxyz
         real(dp),dimension(:),pointer :: pkernelseq
         real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
       end subroutine minimize_in_subspace

       subroutine applyOrthoconstraintlocal(iproc, nproc, lzd, orbs, op, lagmat, lphiovrlp, lhphi)
         use module_base
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc
         type(local_zone_descriptors),intent(in):: lzd
         type(orbitals_data),intent(in):: orbs
         type(overlapParameters),intent(in):: op
         real(8),dimension(orbs%norb,orbs%norb),intent(inout):: lagmat
         real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
         real(8),dimension(orbs%npsidim_comp),intent(out):: lhphi
       end subroutine applyOrthoconstraintlocal


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


      subroutine indicesForExpansion(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(orbitals_data),intent(in):: orbs
        type(input_variables),intent(in):: input
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(inout):: op
        type(p2pComms),intent(in):: comon
      end subroutine indicesForExpansion


      subroutine nullify_expansionSegments(expseg)
        use module_base
        use module_types
        implicit none
        type(expansionSegments),intent(out):: expseg
      end subroutine nullify_expansionSegments


      subroutine indicesForExtraction(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, comon)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, sizePhi
        type(orbitals_data),intent(in):: orbs
        integer,dimension(orbs%norb),intent(in):: onWhichAtom
        type(local_zone_descriptors),intent(in):: lzd
        type(overlapParameters),intent(inout):: op
        type(p2pComms),intent(out):: comon
      end subroutine indicesForExtraction


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


      subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, &
           nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
           hgrid, offsetx, offsety, offsetz, &
           ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
           rxyzConf, potentialPrefac,  withKinetic, cprecr, &
           xx_c, xx_f1, xx_f, &
           xy_c, xy_f2, xy_f, &
           xz_c, xz_f4, xz_f, &
           y_c, y_f)
        use module_base
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, offsetx, offsety, offsetz
        real(gp), intent(in) :: hgrid, potentialPrefac, cprecr
        logical,intent(in):: withKinetic
        real(8),dimension(3):: rxyzConf
        integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
        integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
        integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
        real(wp),dimension(0:n1,0:n2,0:n3),intent(in):: xx_c
        real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f1
        real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f
        real(wp),dimension(0:n2,0:n1,0:n3),intent(in):: xy_c
        real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f2
        real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f
        real(wp),dimension(0:n3,0:n1,0:n2),intent(in):: xz_c
        real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f4
        real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f
        real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
        real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
      end subroutine ConvolQuartic4


      subroutine deallocate_collectiveComms(collComms, subname)
        use module_base
        use module_types
        implicit none
        type(collectiveComms),intent(inout):: collComms
        character(len=*),intent(in):: subname
      end subroutine deallocate_collectiveComms


      subroutine flatten(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
           rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, cut, alpha, &
           ibyyzz_r) !optional
        use module_base
        implicit none
        integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor, confPotOrder, offsetx, offsety, offsetz
        real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
        integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
        real(8),dimension(3),intent(in):: rxyzConfinement
        real(8),intent(in):: hxh, hyh, hzh, potentialPrefac, cut, alpha
      end subroutine flatten


      subroutine sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc
        type(orbitals_data),intent(in):: orbs
        type(locreg_descriptors),intent(in):: Glr
        type(input_variables),intent(in):: input
        type(linearParameters),intent(inout):: lin
        real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
        real(8),dimension(lin%lb%orbs%npsidim_comp),intent(in):: phi
        type(atoms_data),intent(in):: at
        integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      end subroutine sumrholinear_auxiliary


      subroutine sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, nrho, rho, at, nscatterarr)
        use module_base
        use module_types
        implicit none
        integer,intent(in):: iproc, nproc, nrho
        type(orbitals_data),intent(in):: orbs
        type(locreg_descriptors),intent(in):: Glr
        type(input_variables),intent(in):: input
        type(linearParameters),intent(inout):: lin
        real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
        real(8),dimension(nrho),intent(out),target:: rho
        type(atoms_data),intent(in):: at
        integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      end subroutine sumrholinear_withauxiliary



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
         real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,npot), intent(in) :: pot
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
       
       subroutine check_linear_and_create_Lzd(iproc,nproc,input,Lzd,atoms,orbs,rxyz)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(input_variables), intent(in) :: input
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(in) :: atoms
         type(orbitals_data),intent(inout) :: orbs
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
       end subroutine check_linear_and_create_Lzd

       subroutine reinitialize_Lzd_after_LIG(iproc,nproc,input,Lzd,atoms,orbs,rxyz)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(input_variables), intent(in) :: input
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(atoms_data), intent(in) :: atoms
         type(orbitals_data),intent(inout) :: orbs
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
       end subroutine reinitialize_Lzd_after_LIG

       subroutine system_initialization(iproc,nproc,in,atoms,rxyz,&
            inputpsi,input_wf_format,orbs,Lzd,denspot,nlpspd,comms,hgrids,shift,proj,radii_cf)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc 
         type(input_variables), intent(in) :: in 
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
         integer, intent(out) :: inputpsi !<strategy for wavefunction input guess
         integer, intent(out) :: input_wf_format !< format of input wavefunctions
         type(orbitals_data), intent(out) :: orbs
         type(local_zone_descriptors), intent(out) :: Lzd
         type(DFT_local_fields), intent(out) :: denspot
         type(nonlocal_psp_descriptors), intent(out) :: nlpspd
         type(communications_arrays), intent(out) :: comms
         real(gp), dimension(3), intent(out) :: hgrids !< grid spacings of the daubechies grid
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

       subroutine extract_potential_for_spectra(iproc,nproc,at,rhod,dpcom,&
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
         type(denspot_distribution), intent(in) :: dpcom
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
         real(dp), dimension(:), pointer :: pkernel,pkernelseq
         integer, intent(in) ::potshortcut
       end subroutine extract_potential_for_spectra

       subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,hxh,hyh,hzh,itrp,iscf,alphamix,mix,ixc,&
            nlpspd,proj,rxyz,linflag,exctxpar,unblock_comms,hx,hy,hz,Lzd,orbs,SIC,confdatarr,GPU,optscf,psi,&
            ekin_sum,epot_sum,eexctX,eSIC_DC,eproj_sum,ehart,eexcu,vexcu,rpnrm,xcstr,hpsi,proj_G,paw)
         use module_base
         use module_types
         use m_ab6_mixing
         implicit none
         logical, intent(in) :: scf
         integer, intent(in) :: optscf
         integer, intent(in) :: iproc,nproc,itrp,iscf,ixc,linflag
         real(gp), intent(in) :: hx,hy,hz,hxh,hyh,hzh,alphamix
         character(len=3), intent(in) :: unblock_comms
         type(atoms_data), intent(in) :: atoms
         type(nonlocal_psp_descriptors), intent(in) :: nlpspd
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(ab6_mixing_object), intent(inout) :: mix
         type(DFT_local_fields), intent(inout) :: denspot
         type(SIC_data), intent(in) :: SIC
         character(len=*), intent(in) :: exctxpar
         real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
         real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
         type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
         type(GPU_pointers), intent(inout) :: GPU
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
         real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC,ehart,eexcu,vexcu,rpnrm
         real(gp), dimension(6), intent(out) :: xcstr
         real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
         type(gaussian_basis),dimension(atoms%nat),intent(in)::proj_G
         type(paw_objects),intent(inout)::paw
       end subroutine psitohpsi

       subroutine allocate_diis_objects(idsx,alphadiis,npsidim,nkptsp,nspinor,diis,subname) !n(m)
         use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: subname
         integer, intent(in) :: idsx,npsidim,nkptsp,nspinor !n(m)
         real(gp), intent(in) :: alphadiis
         type(diis_objects), intent(inout) :: diis
       end subroutine allocate_diis_objects

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
         real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor), intent(inout) :: psi
         real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
         real(gp), dimension(3,at%nat), intent(inout) :: fsep
         real(gp), dimension(6), intent(out) :: strten
       end subroutine nonlocal_forces

       !> Calculates the local forces acting on the atoms belonging to iproc
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

       subroutine initialize_DFT_local_fields(denspot)
         use module_base
         use module_types
         implicit none
         type(DFT_local_fields), intent(out) :: denspot
       end subroutine initialize_DFT_local_fields

       subroutine local_potential_dimensions(Lzd,orbs,ndimfirstproc)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: ndimfirstproc
         type(local_zone_descriptors), intent(inout) :: Lzd
         type(orbitals_data), intent(inout) :: orbs
       end subroutine local_potential_dimensions

       subroutine check_communications(iproc,nproc,orbs,lr,comms)
         use module_base
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         type(orbitals_data), intent(in) :: orbs
         type(locreg_descriptors), intent(in) :: lr
         type(communications_arrays), intent(in) :: comms
       end subroutine check_communications


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

!> Mix the electronic density or the potential using DIIS
     subroutine mix_rhopot(iproc,nproc,npoints,alphamix,mix,rhopot,istep,&
          & n1,n2,n3,ucvol,rpnrm,nscatterarr)
       use module_base
       use defs_basis, only: AB6_NO_ERROR
       use m_ab6_mixing
       implicit none
       integer, intent(in) :: npoints, istep, n1, n2, n3, nproc, iproc
       real(gp), intent(in) :: alphamix, ucvol
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
       type(ab6_mixing_object), intent(inout) :: mix
       real(dp), dimension(npoints), intent(inout) :: rhopot
       real(gp), intent(out) :: rpnrm
     end subroutine mix_rhopot

   end interface

END MODULE module_interfaces

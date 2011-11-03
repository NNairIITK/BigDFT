!> @file
!! Define the module module_interfaces containing all interfaces
!!
!! @author 
!!    Luigi Genovese, Damien Caliste
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>  Modules which contains all interfaces
!!  Interfaces of:
!!  - call_cluster
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
!!  - HamiltonianApplication
!!  - hpsitopsi
!!  - last_orthon
!!  - local_forces
!!  - projectors_derivatives
!!  - nonlocal_forces
!!  - CalculateTailCorrection
!!  - reformatonewave
module module_interfaces

  implicit none

  interface

     subroutine call_bigdft(nproc,iproc,atoms,rxyz,in,energy,fxyz,fnoise,rst,infocode)
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
       real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
     END SUBROUTINE call_bigdft


     subroutine geopt(nproc,iproc,x,at,f,epot,rst,in,ncount_bigdft)

       !    use module_base
       !    use module_interfaces, except_this_one => geopt
       !    use module_types
       !    use minimization, only:parameterminimization

       !n(c) use module_base
       use module_types
       !    use minimization, only:parameterminimization
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(inout) :: ncount_bigdft
       type(atoms_data), intent(in) :: at
       type(input_variables), intent(in) :: in
       type(restart_objects), intent(inout) :: rst
       real(gp), intent(inout) :: epot
       real(gp), dimension(3*at%nat), intent(inout) :: x
       real(gp), dimension(3*at%nat), intent(out) :: f
     END SUBROUTINE geopt


     subroutine copy_old_wavefunctions(nproc,orbs,n1,n2,n3,wfd,psi,&
          n1_old,n2_old,n3_old,wfd_old,psi_old)
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

     subroutine dft_input_variables(iproc,filename,in)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(out) :: in
     END SUBROUTINE dft_input_variables

     subroutine geopt_input_variables(filename,in)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       type(input_variables), intent(inout) :: in
     END SUBROUTINE geopt_input_variables

     subroutine tddft_input_variables(filename,in)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       type(input_variables), intent(inout) :: in
     END SUBROUTINE tddft_input_variables


     subroutine kpt_input_variables(iproc,filename,in,atoms)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: in
       type(atoms_data), intent(in) :: atoms
     END SUBROUTINE kpt_input_variables

     subroutine perf_input_variables(iproc,filename,inputs)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: inputs
     END SUBROUTINE perf_input_variables

     subroutine read_atomic_file(file,iproc,at,rxyz)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: file
       integer, intent(in) :: iproc
       type(atoms_data), intent(inout) :: at
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_atomic_file

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
          end subroutine getline
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
          end subroutine getline
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
     end subroutine orbitals_descriptors

     subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms,basedist)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       type(locreg_descriptors), intent(in) :: lr
       type(orbitals_data), intent(inout) :: orbs
       type(communications_arrays), intent(out) :: comms
       integer, dimension(0:nproc-1,orbs%nkpts), intent(in), optional :: basedist
     end subroutine orbitals_communicators
     
     subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
          crmult,frmult,Glr,output_grid)
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
       logical, intent(in), optional :: output_grid
     END SUBROUTINE createWavefunctionsDescriptors

     subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
          radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,n1,n2,n3
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       real(kind=8), intent(in) :: cpmult,fpmult,hx,hy,hz
       real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
       real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
       type(nonlocal_psp_descriptors), intent(out) :: nlpspd
       real(kind=8), dimension(:), pointer :: proj
     END SUBROUTINE createProjectorsArrays

     subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
          rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
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
          radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
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

     end subroutine createPcProjectorsArrays

     
     subroutine applyPCprojectors(orbs,at,&
          rxyz,hx,hy,hz,Glr,PPD,psi,hpsi, dotest)
       
       use module_base
       use module_types
       
       type(orbitals_data), intent(inout) :: orbs
       type(atoms_data) :: at
       real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors), intent(in) :: Glr
       type(pcproj_data_type) ::PPD
       real(wp), dimension(:), pointer :: psi, hpsi
       logical, optional :: dotest
     end subroutine applyPCprojectors

     

     subroutine applyPAWprojectors(orbs,at,&
          rxyz,hx,hy,hz,Glr,PAWD,psi,hpsi,  paw_matrix, dosuperposition , &
          sup_iatom, sup_l, sup_arraym) !, sup_arraychannel)
       
       use module_base
       use module_types
       
       type(orbitals_data), intent(inout) :: orbs
       type(atoms_data) :: at
       real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors), intent(in) :: Glr
       type(pawproj_data_type) ::PAWD
       real(wp), dimension(:), pointer :: psi, hpsi, paw_matrix
       logical dosuperposition
       integer, optional :: sup_iatom, sup_l
       real(wp) , dimension(:), pointer, optional :: sup_arraym !, sup_arraychannel

     end subroutine applyPAWprojectors


     subroutine IonicEnergyandForces(iproc,nproc,at,hxh,hyh,hzh,elecfield,rxyz,eion,fion,psoffset,&
          nvacancy,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
       !n(c) use module_base
       use module_types
       implicit none
       type(atoms_data), intent(in) :: at
       integer, intent(in) :: iproc,nproc,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,nvacancy
       real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield
       real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
       real(kind=8), dimension(*), intent(in) :: pkernel
       real(kind=8), intent(out) :: eion,psoffset
       real(kind=8), dimension(3,at%nat), intent(out) :: fion
       real(kind=8), dimension(*), intent(out) :: pot_ion
     END SUBROUTINE IonicEnergyandForces

     subroutine createIonicPotential(geocode,iproc,nproc,at,rxyz,&
          hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,psoffset,nvacancy,&
          correct_offset)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=1), intent(in) :: geocode
       logical, intent(in) :: correct_offset
       integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,nvacancy
       real(gp), intent(in) :: hxh,hyh,hzh,psoffset
       type(atoms_data), intent(in) :: at
       real(gp), intent(in) :: elecfield
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(dp), dimension(*), intent(in) :: pkernel
       real(wp), dimension(*), intent(inout) :: pot_ion
     END SUBROUTINE createIonicPotential

     subroutine input_wf_diag(iproc,nproc,at,rhodsc,&
          orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
          nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ixc,symObj
       integer, intent(inout) :: nspin,nvirt
       real(gp), intent(in) :: hx,hy,hz
       type(rho_descriptors),intent(in) :: rhodsc
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
       real(wp), dimension(:), pointer :: psi,hpsi,psit,rhocore
       real(dp), dimension(:), pointer :: pkernel,pkernelseq
       integer, intent(in) :: potshortcut
       integer, dimension(*), intent(in) :: irrzon
       real(dp), dimension(*), intent(in) :: phnons
     END SUBROUTINE input_wf_diag

     subroutine reformatmywaves(iproc,orbs,at,&
          hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
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

     subroutine sumrho(iproc,nproc,orbs,lr,hxh,hyh,hzh,psi,rho, &
          nscatterarr,nspin,GPU,symObj,irrzon,phnons,rhodsc)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,symObj
       real(gp), intent(in) :: hxh,hyh,hzh
       type(rho_descriptors),intent(in) :: rhodsc
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: lr 
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
       real(dp), dimension(max(lr%d%n1i*lr%d%n2i*nscatterarr(iproc,1),1),nspin), intent(out), target :: rho
       type(GPU_pointers), intent(inout) :: GPU
       integer, dimension(*), intent(in) :: irrzon
       real(dp), dimension(*), intent(in) :: phnons
     END SUBROUTINE sumrho


     subroutine rho_segkey(iproc,at,rxyz,crmult,frmult,radii_cf,&
         n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,nspin,rho_d,iprint)
       !n(c) use module_base
       use module_types
       implicit none
       integer,intent(in) :: n1,n2,n3,n1i,n2i,n3i,iproc,nspin
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
       real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
       logical,intent(in) :: iprint
       type(rho_descriptors),intent(inout) :: rho_d
      end subroutine rho_segkey

      subroutine LocalHamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
           lr,ngatherarr,pot,psi,hpsi,ekin_sum,epot_sum,eexctX,eSIC_DC,SIC,GPU,pkernel,orbsocc,psirocc)
        use module_base
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc
        real(gp), intent(in) :: hx,hy,hz
        type(atoms_data), intent(in) :: at
        type(orbitals_data), intent(in) :: orbs
        type(locreg_descriptors), intent(in) :: lr 
        type(SIC_data), intent(in) :: SIC
        integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
        real(gp), dimension(3,at%nat), intent(in) :: rxyz
        real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
        real(wp), dimension(:), pointer :: pot
        real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
        real(gp), intent(inout) :: eexctX !used to activate the OP2P scheme
        real(wp), target, dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(out) :: hpsi
        type(GPU_pointers), intent(inout) :: GPU
        real(dp), dimension(:), pointer, optional :: pkernel
        type(orbitals_data), intent(in), optional :: orbsocc
        real(wp), dimension(:), pointer, optional :: psirocc
      END SUBROUTINE LocalHamiltonianApplication

     subroutine hpsitopsi(iproc,nproc,orbs,lr,comms,iter,diis,idsx,psi,psit,hpsi,nspin,orthpar)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,idsx,iter,nspin
       type(locreg_descriptors), intent(in) :: lr
       type(communications_arrays), intent(in) :: comms
       type(orbitals_data), intent(in) :: orbs
       type(orthon_data), intent(in) :: orthpar
       type(diis_objects), intent(inout) :: diis
       real(wp), dimension(:), pointer :: psi,psit,hpsi
     END SUBROUTINE hpsitopsi

     subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,wfd,comms,&
          psi,hpsi,psit,orthpar,passmat,& !mandatory
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
          nspin,comms,psi,hpsi,psit,evsum, keeppsit)
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

     subroutine local_forces(iproc,at,rxyz,hxh,hyh,hzh,&
          n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,rho,pot,floc)
       ! Calculates the local forces acting on the atoms belonging to iproc
       use module_types
       implicit none
       !Arguments---------
       type(atoms_data), intent(in) :: at
       integer, intent(in) :: iproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
       real(kind=8), intent(in) :: hxh,hyh,hzh
       real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
       real(kind=8), dimension(*), intent(in) :: rho,pot
       real(kind=8), dimension(3,at%nat), intent(out) :: floc
     END SUBROUTINE local_forces

     subroutine nonlocal_forces(iproc,n1,n2,n3,hx,hy,hz,at,rxyz,&
          orbs,nlpspd,proj,wfd,psi,fsep,refill)
       !n(c) use module_base
       use module_types
       implicit none
       !Arguments-------------
       type(atoms_data), intent(in) :: at
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       logical, intent(in) :: refill
       integer, intent(in) :: iproc,n1,n2,n3
       real(gp), intent(in) :: hx,hy,hz
       type(orbitals_data), intent(in) :: orbs
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
       real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
       real(gp), dimension(3,at%nat), intent(inout) :: fsep
     END SUBROUTINE nonlocal_forces

     subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
          Glr,nlpspd,ncongt,pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,&
          proj,psi,output_grid,ekin_sum,epot_sum,eproj_sum)
       !n(c) use module_base
       use module_types
       implicit none
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: Glr
       type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
       integer, intent(in) :: iproc,nproc,ncongt,nspin
       logical, intent(in) :: output_grid
       real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
       real(kind=8), dimension(at%ntypes,3), intent(in) :: radii_cf
       real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
       real(kind=8), dimension(Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nspin), intent(in) :: pot
       real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
       real(kind=8), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
       real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     END SUBROUTINE CalculateTailCorrection

     !added for abinit compatilbility
     subroutine reformatonewave(iproc,displ,wfd,at,hx_old,hy_old,hz_old,&
          n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,n1_old,n2_old,n3_old,n1,n2,n3
       real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%nat), intent(in) :: rxyz_old,rxyz
       real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
       real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
       real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
     END SUBROUTINE reformatonewave
     subroutine readonewave(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
          & hx,hy,hz,at,wfd,rxyz_old,rxyz,psi,eval,psifscf)
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

     subroutine davidson(iproc,nproc,n1i,n2i,in,at,&
          orbs,orbsv,nvirt,lr,comms,commsv,&
          hx,hy,hz,rxyz,rhopot,nlpspd,proj,pkernel,psi,v,nscatterarr,ngatherarr,GPU)
       !n(c) use module_base
       use module_types
       !n(c) use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc,n1i,n2i
       integer, intent(in) :: nvirt
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: lr 
       type(orbitals_data), intent(in) :: orbs
       type(communications_arrays), intent(in) :: comms, commsv
       real(gp), intent(in) :: hx,hy,hz
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(:), pointer :: pkernel
       real(dp), dimension(*), intent(in) :: rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       real(wp), dimension(:), pointer :: psi,v
     END SUBROUTINE davidson

     subroutine build_eigenvectors(iproc,norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinore,nspinor,&
          ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,passmat,nvirte,psivirt)
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

     subroutine preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ncong
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors), intent(in) :: lr
       type(orbitals_data), intent(in) :: orbs
       real(dp), intent(out) :: gnrm,gnrm_zero
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp,orbs%nspinor), intent(inout) :: hpsi
     END SUBROUTINE preconditionall

     subroutine transpose_v(iproc,nproc,orbs,wfd,comms,psi,&
          work,outadd) !optional
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

     subroutine untranspose_v(iproc,nproc,orbs,wfd,comms,psi,&
          work,outadd) !optional
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

     subroutine plot_wf(orbname,nexpo,at,lr,hx,hy,hz,rxyz,psi,comment)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*) :: comment
       character(len=*) :: orbname
       integer, intent(in) :: nexpo
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       type(locreg_descriptors), intent(in) :: lr
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
     END SUBROUTINE plot_wf

     subroutine partial_density_free(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
          hfac,nscatterarr,spinsgn,psir,rho_p,ibyyzz_r) !ex-optional argument
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
          CP2K,wfn_cp2k)
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

     subroutine restart_from_gaussians(iproc,nproc,orbs,lr,hx,hy,hz,psi,G,coeffs)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       real(gp), intent(in) :: hx,hy,hz
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: lr
       type(gaussian_basis), intent(inout) :: G
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp), intent(out) :: psi
       real(wp), dimension(:,:), pointer :: coeffs
     END SUBROUTINE restart_from_gaussians

     subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin,&
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks)
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

     subroutine AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,&
          nspin,eks,scorb,G,gaucoeff,iorbtolr)
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
          ibyyzz_r) !optional
       use module_base
       implicit none
       integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
       real(wp), dimension(-nl1:2*n1+2+nl1,-nl2:2*n2+2+nl2,-nl3:2*n3+2+nl3,nspinor), intent(inout) :: psir
       real(wp), dimension(-nl1:2*n1+2+nl1-4*nbuf,-nl2:2*n2+2+nl2-4*nbuf,-nl3:2*n3+2+nl3-4*nbuf,npot), intent(in) :: pot
       integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
       real(gp), intent(out) :: epot
     END SUBROUTINE apply_potential

     subroutine correct_hartree_potential(at,iproc,nproc,n1i,n2i,n3i,n3p,n3pi,n3d,&
          i3s,i3xcsh,hxh,hyh,hzh,pkernel,ngatherarr,&
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
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
          ekin_sum,epot_sum,eproj_sum,nspin,GPU, in_iat_absorber, in, PAWD, orbs )
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ndimpot,nspin
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in), target :: at
       type(nonlocal_psp_descriptors), intent(in) , target :: nlpspd
       type(locreg_descriptors), intent(in) , target :: lr 
       integer, dimension(0:nproc-1,2), intent(in) , target :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) , target :: rxyz
       real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
       real(wp), dimension(nlpspd%nprojel), intent(in) , target :: proj
       real(wp), dimension(max(ndimpot,1),nspin), target :: potential
       real(gp), intent(out) :: ekin_sum,epot_sum,eproj_sum
       type(GPU_pointers), intent(inout) , target :: GPU
       integer, intent(in) :: in_iat_absorber
       type(pawproj_data_type), target ::PAWD

       type(input_variables),intent(in), target :: in
       type(orbitals_data), intent(inout), target :: orbs
     END SUBROUTINE xabs_lanczos
     subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
          zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
          aeval,ng,psi,res,chrg,&
          Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw, PAWpatch, &
          psipsigrid)
       use module_base, only: gp
       
       implicit real(gp) (a-h,o-z)
       logical :: noproj, readytoexit
       integer, parameter :: n_int=1000
       dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
            hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
            gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
            pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
            potgrd(n_int),&
            rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
            occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
            vh(0:ng,0:ng,4,0:ng,0:ng,4),&
            res(noccmax,lmax+1),xp(0:ng),& 
            psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol),&
            psigrid_naked_2(Ngrid,Nsol), projgrid(Ngrid,3), &
            rhogrid(Ngrid), potgrid(Ngrid), psigrid_not_fitted(Ngrid,Nsol),&
            psigrid_not_fitted_2(Ngrid,Nsol),&
            vxcgrid(Ngrid), &
            Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
            H(Nsol, Nsol), &
            H_2(Nsol, Nsol), &
            Hcorrected(Nsol, Nsol), &
            Hadd(Nsol, Nsol), Egrid_tmp(Nsol),Egrid_tmp_2(Nsol), Etofit(Nsol), &
            Soverlap(Nsol,Nsol), Tpsigrid(Nsol,Ngrid ),Tpsigrid_dum(Nsol, Ngrid),valuesatp(Nsol), &
            PAWpatch(Npaw, Npaw ), Spsitildes(Npaw, Npaw), genS(Nsol,Nsol), genH(Nsol,Nsol) , dumH(Nsol,Nsol)
       
       real(gp) , optional :: psipsigrid(Ngrid, Nsol)
       
       
       real(gp) :: rgrid(Ngrid), ene_m, ene_p, factadd, rcond, fixfact
       real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid), dumgrid3(Ngrid)
       logical dofit
       integer real_start, iocc, iwork(Nsol), INFO, volta, ngrid_box_2
       character(1) EQUED
       integer ipiv(Nsol), Npaw
     end subroutine gatom_modified

     subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,&
          psi, aeval, occup, psp_modifier, &
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
            rgrid(Ngrid), psigrid(Ngrid,Nsol  )
       real(gp),   intent(out), optional  :: psipsigrid(Ngrid,Nsol  )
       real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
       real(gp):: PAWpatch(Npaw,Npaw)
       
       !local variables
     end subroutine abs_generator_modified

     subroutine xabs_cg(iproc,nproc,at,hx,hy,hz,rxyz,&
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
          ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,&
          in , rhoXanes, PAWD , PPD, orbs )
       use module_base
       use module_types
       ! per togliere il bug 
       
       implicit none
       
       integer  :: iproc,nproc,ndimpot,nspin
       real(gp)  :: hx,hy,hz
       type(atoms_data), target :: at
       type(nonlocal_psp_descriptors), target :: nlpspd
       type(locreg_descriptors), target :: lr
       
       type(pcproj_data_type), target ::PPD
       
       integer, dimension(0:nproc-1,2), target :: ngatherarr 
       real(gp), dimension(3,at%nat), target :: rxyz
       real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
       real(wp), dimension(nlpspd%nprojel), target :: proj
       real(wp), dimension(max(ndimpot,1),nspin), target :: potential
       real(wp), dimension(max(ndimpot,1),nspin), target :: rhoXanes
       
       
       
       real(gp) :: ekin_sum,epot_sum,eproj_sum
       type(GPU_pointers), intent(inout) , target :: GPU
       integer, intent(in) :: in_iat_absorber
       type(pawproj_data_type), target ::PAWD
       type(input_variables),intent(in), target :: in
       type(orbitals_data), intent(inout), target :: orbs
     end subroutine xabs_cg
     
     subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
          ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in, PAWD,orbs   )! aggiunger a interface
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


       type(input_variables),intent(in), target :: in
       type(pawproj_data_type), target ::PAWD
       type(orbitals_data), intent(inout), target :: orbs

     END SUBROUTINE xabs_chebychev

     subroutine cg_spectra(iproc,nproc,at,hx,hy,hz,rxyz,&
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
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
          hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
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
          nat,rxyz,iatypes, znucl)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       character(len=1), intent(in) :: geocode
       integer, intent(out) :: nspin
       integer, intent(out) ::  n1i,n2i,n3i
       real(gp), intent(out) :: hxh,hyh,hzh
       real(dp), dimension(:,:), pointer :: rho
       real(gp), dimension(:,:), pointer, optional :: rxyz
       integer, intent(out), optional ::  nat
       integer, dimension(:), pointer, optional :: iatypes, znucl
     END SUBROUTINE read_density

     subroutine read_cube(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
          nat,rxyz, iatypes, znucl)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       character(len=1), intent(in) :: geocode
       integer, intent(out) :: nspin
       integer, intent(out) ::  n1i,n2i,n3i
       real(gp), intent(out) :: hxh,hyh,hzh
       real(dp), dimension(:,:), pointer :: rho
       real(gp), dimension(:,:), pointer   :: rxyz
       integer, intent(out)   ::  nat
       integer, dimension(:), pointer   :: iatypes, znucl
     END SUBROUTINE read_cube

     subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
          nat,rxyz, iatypes, znucl)
       !n(c) use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       character(len=1), intent(in) :: geocode
       integer, intent(out) :: nspin
       integer, intent(out) ::  n1i,n2i,n3i
       real(gp), intent(out) :: hxh,hyh,hzh
       real(dp), dimension(:,:), pointer :: rho
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
          iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg )
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
     end subroutine gaussian_pswf_basis

     subroutine gaussian_pswf_basis_for_paw(iproc,nspin,at,rxyz,G,  &
          iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels,&
          iorbto_imatrixbeg )
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nspin
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
     end subroutine gaussian_pswf_basis_for_paw


     subroutine local_analysis(iproc,nproc,hx,hy,hz,in,at,rxyz,shift,lr,orbs,orbsv,psi,psivirt)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       real(gp), intent(in) :: hx,hy,hz
       type(input_variables), intent(in) :: in
       type(locreg_descriptors), intent(in) :: lr
       type(orbitals_data), intent(in) :: orbs,orbsv
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3),intent(in) :: shift
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
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,i3s,n3d,i3xcsh,n3p
       real(gp), intent(in) :: hxh,hyh,hzh
       type(atoms_data), intent(in) :: at
       type(grid_dimensions), intent(in) :: d
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(:), pointer :: rhocore
     END SUBROUTINE calculate_rhocore

     subroutine XC_potential(geocode,datacode,iproc,nproc,n01,n02,n03,ixc,hx,hy,hz,&
          rho,exc,vxc,nspin,rhocore,potxc,dvxcdrho)
       use module_base
       implicit none
       character(len=1), intent(in) :: geocode
       character(len=1), intent(in) :: datacode
       integer, intent(in) :: iproc,nproc,n01,n02,n03,ixc,nspin
       real(gp), intent(in) :: hx,hy,hz
       real(gp), intent(out) :: exc,vxc
       real(dp), dimension(*), intent(inout) :: rho
       real(wp), dimension(:), pointer :: rhocore !associated if useful
       real(wp), dimension(*), intent(out) :: potxc
       real(dp), dimension(:,:,:,:), intent(out), target, optional :: dvxcdrho
     END SUBROUTINE XC_potential

     subroutine direct_minimization(iproc,nproc,n1i,n2i,in,at,&
          orbs,orbsv,nvirt,lr,comms,commsv,&
          hx,hy,hz,rxyz,rhopot,nlpspd,proj, &
          pkernel,psi,psivirt,nscatterarr,ngatherarr,GPU)
       !n(c) use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,n1i,n2i,nvirt
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: lr 
       type(orbitals_data), intent(in) :: orbs
       type(communications_arrays), intent(in) :: comms, commsv
       real(gp), intent(in) :: hx,hy,hz
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(:), pointer :: pkernel
       real(dp), dimension(*), intent(in), target :: rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       real(wp), dimension(:), pointer :: psi,psivirt
     END SUBROUTINE direct_minimization

     subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
          hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
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

    subroutine full_local_potential(iproc,nproc,ndimpot,ndimgrid,nspin,ndimrhopot,i3rho_add,norb,norbp,ngatherarr,potential,pot)
      use module_base
      implicit none
      integer, intent(in) :: iproc,nproc,nspin,ndimpot,norb,norbp,ndimgrid
      integer, intent(in) :: ndimrhopot,i3rho_add
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      real(wp), dimension(max(ndimrhopot,1)), intent(in), target :: potential
      real(wp), dimension(:), pointer :: pot
    END SUBROUTINE full_local_potential

    subroutine free_full_potential(nproc,pot,subname)
      use module_base
      implicit none
      character(len=*), intent(in) :: subname
      integer, intent(in) :: nproc
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
      real(wp), dimension(orbs%npsidim), intent(in) :: psi
      type(orbitals_data), intent(out) :: orbs_as
      type(communications_arrays), intent(out) :: comms_as
      real(wp), dimension(:), pointer :: psi_as
    END SUBROUTINE select_active_space

    subroutine calculate_energy_and_gradient(iter,iproc,nproc,orbs,comms,GPU,lr,hx,hy,hz,ncong,iscf,&
         ekin,epot,eproj,eSIC_DC,ehart,exc,evxc,eexctX,eion,edisp,psi,psit,hpsi,gnrm,gnrm_zero,energy)
      !n(c) use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,ncong,iscf,iter
      real(gp), intent(in) :: hx,hy,hz,ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp,eSIC_DC
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(locreg_descriptors), intent(in) :: lr
      type(GPU_pointers), intent(in) :: GPU
      real(gp), intent(out) :: gnrm,gnrm_zero,energy
      real(wp), dimension(:), pointer :: psi,psit,hpsi
    end subroutine calculate_energy_and_gradient

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
!!    end subroutine calculate_energy_and_gradient_new

    subroutine constrained_davidson(iproc,nproc,n1i,n2i,in,at,&
         orbs,orbsv,nvirt,lr,comms,commsv,&
         hx,hy,hz,rxyz,rhopot,nlpspd,proj,pkernel,psi,v,nscatterarr,ngatherarr,GPU)
      !n(c) use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,n1i,n2i
      integer, intent(in) :: nvirt
      type(input_variables), intent(in) :: in
      type(atoms_data), intent(in) :: at
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(locreg_descriptors), intent(in) :: lr 
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms, commsv
      real(gp), intent(in) :: hx,hy,hz
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(dp), dimension(*), intent(in) :: rhopot
      real(dp), dimension(:), pointer :: pkernel
      type(orbitals_data), intent(inout) :: orbsv
      type(GPU_pointers), intent(inout) :: GPU
      real(wp), dimension(:), pointer :: psi,v
    end subroutine constrained_davidson

    subroutine local_hamiltonian(iproc,orbs,lr,hx,hy,hz,&
         ipotmethod,pot,psi,hpsi,pkernel,ixc,alphaSIC,ekin_sum,epot_sum,eSIC_DC)
      !n(c) use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,ipotmethod,ixc
      real(gp), intent(in) :: hx,hy,hz,alphaSIC
      type(orbitals_data), intent(in) :: orbs
      type(locreg_descriptors), intent(in) :: lr
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
      real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
      real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
      real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
      real(dp), dimension(:), pointer :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
    end subroutine local_hamiltonian

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
    end subroutine NK_SIC_potential
    
     subroutine center( vector, vecsize )
     
       !use defs, only : natoms, constr
       !use bigdft_forces, only : in_system
       implicit none
     
       !Arguments
       integer, intent(in) :: vecsize
       real(kind=8), dimension(vecsize), intent(inout), target :: vector
     end subroutine center


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
     !end subroutine SWcalczone

  end interface

end module module_interfaces


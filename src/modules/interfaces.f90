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
!!  - input_occup
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
       use module_base
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

       use module_base
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
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nproc,n1,n2,n3
       type(orbitals_data), intent(in) :: orbs
       type(wavefunctions_descriptors), intent(inout) :: wfd,wfd_old
       integer, intent(out) :: n1_old,n2_old,n3_old
       real(wp), dimension(:), pointer :: psi,psi_old
     END SUBROUTINE copy_old_wavefunctions

     subroutine system_properties(iproc,nproc,in,at,orbs,radii_cf,nelec)
       use module_base
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
       use module_base
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

     subroutine read_input_variables(iproc,posinp,inputs,atoms,rxyz)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: posinp
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: inputs
       type(atoms_data), intent(out) :: atoms
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_input_variables

     subroutine read_input_parameters(iproc,inputs,atoms,rxyz)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: inputs
       type(atoms_data), intent(inout) :: atoms
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_input_parameters

     subroutine dft_input_variables(iproc,filename,in)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(out) :: in
     END SUBROUTINE dft_input_variables

     subroutine geopt_input_variables(filename,in)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       type(input_variables), intent(inout) :: in
     END SUBROUTINE geopt_input_variables

     subroutine tddft_input_variables(filename,in)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       type(input_variables), intent(inout) :: in
     END SUBROUTINE tddft_input_variables


     subroutine kpt_input_variables(iproc,filename,in,atoms)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: in
       type(atoms_data), intent(in) :: atoms
     END SUBROUTINE kpt_input_variables

     subroutine perf_input_variables(iproc,filename,inputs)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc
       type(input_variables), intent(inout) :: inputs
     END SUBROUTINE perf_input_variables

     subroutine read_atomic_file(file,iproc,at,rxyz)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: file
       integer, intent(in) :: iproc
       type(atoms_data), intent(inout) :: at
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_atomic_file

     subroutine read_ascii_positions(iproc,ifile,atoms,rxyz)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,ifile
       type(atoms_data), intent(inout) :: atoms
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_ascii_positions

     subroutine write_atomic_file(filename,energy,rxyz,atoms,comment)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename,comment
       type(atoms_data), intent(in) :: atoms
       real(gp), intent(in) :: energy
       real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
     END SUBROUTINE write_atomic_file

     subroutine MemoryEstimator(nproc,idsx,lr,nat,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,peakmem)
       use module_base
       use module_types
       implicit none
       !Arguments
       integer, intent(in) :: nproc,idsx,nat,norb,nspin,nprojel
       integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
       type(locreg_descriptors), intent(in) :: lr
       real(kind=8), intent(out) :: peakmem
     END SUBROUTINE MemoryEstimator

     subroutine check_closed_shell(orbs,lcs)
       use module_base
       use module_types
       implicit none
       type(orbitals_data), intent(in) :: orbs
       logical, intent(out) :: lcs
     END SUBROUTINE check_closed_shell

     subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
       integer, intent(in) :: nspinor
       type(orbitals_data), intent(out) :: orbs
       real(gp), dimension(nkpt), intent(in) :: wkpt
       real(gp), dimension(3,nkpt), intent(in) :: kpt
     END SUBROUTINE orbitals_descriptors

     subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
          crmult,frmult,Glr,output_grid)
       use module_base
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

     subroutine createProjectorsArrays(iproc,lr,rxyz,at,orbs,&
          radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
       use module_base
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

     subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
          n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)
       use module_base
       implicit none
       character(len=1), intent(in) :: geocode,datacode
       integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,ixc
       integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
       integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
       integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
     END SUBROUTINE createDensPotDescriptors

     subroutine IonicEnergyandForces(iproc,nproc,at,hxh,hyh,hzh,elecfield,rxyz,eion,fion,psoffset,&
          nvacancy,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
       use module_base
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
       use module_base
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

     subroutine input_wf_diag(iproc,nproc,at,&
          orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
          nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ixc,symObj
       integer, intent(inout) :: nspin,nvirt
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(inout) :: at
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
       use module_base
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

     subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit,input)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       type(orbitals_data), intent(in) :: orbs
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(communications_arrays), intent(in) :: comms
       type(input_variables):: input
       real(wp), dimension(:) , pointer :: psi,hpsi,psit
     END SUBROUTINE first_orthon

     subroutine sumrho(iproc,nproc,orbs,lr,ixc,hxh,hyh,hzh,psi,rho,nrho, &
          & nscatterarr,nspin,GPU,symObj,irrzon,phnons)
       use module_base!, only: gp,dp,wp,ndebug,memocc
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nrho,nspin,ixc,symObj
       real(gp), intent(in) :: hxh,hyh,hzh
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: lr 
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
       real(dp), dimension(max(nrho,1),nspin), intent(out), target :: rho
       type(GPU_pointers), intent(inout) :: GPU
       integer, dimension(*), intent(in) :: irrzon
       real(dp), dimension(*), intent(in) :: phnons
     END SUBROUTINE sumrho

     subroutine HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,&
          nlpspd,proj,lr,ngatherarr,pot,psi,hpsi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel,orbsocc,psirocc)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: lr 
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
       real(wp), dimension(:), pointer :: pot
       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
       real(wp), target, dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       real(dp), dimension(*), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
     END SUBROUTINE HamiltonianApplication

     subroutine HamiltonianApplicationConfinement(iproc,nproc,at,orbs,lin,hx,hy,hz,rxyz,&
          nlpspd,proj,lr,ngatherarr,ndimpot,potential,psi,hpsi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParabola, onWhichAtom, &
          pkernel,orbsocc,psirocc,centralAtom)
       use module_base
       use module_types
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
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f*orbs%nspinor*orbs%norbp), intent(in) :: psi
       real(wp), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
       real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f*orbs%nspinor*orbs%norbp), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       real(gp), dimension(3,at%nat), intent(in) :: rxyzParabola
       integer,dimension(orbs%norb),intent(in):: onWhichAtom
       real(dp), dimension(*), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
       integer,intent(in),optional:: centralAtom
     END SUBROUTINE HamiltonianApplicationConfinement

 subroutine hpsitopsi(iproc,nproc,orbs,lr,comms,iter,diis,idsx,psi,psit,hpsi,nspin,input)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,idsx,iter,nspin
       type(locreg_descriptors), intent(in) :: lr
       type(communications_arrays), intent(in) :: comms
       type(orbitals_data), intent(in) :: orbs
       type(input_variables), intent(in) :: input
       type(diis_objects), intent(inout) :: diis
       real(wp), dimension(:), pointer :: psi,psit,hpsi
     END SUBROUTINE hpsitopsi

     subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,wfd,comms,&
          psi,hpsi,psit,input,& !mandatory
          orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,natsc,nspin
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(communications_arrays), target, intent(in) :: comms
       type(orbitals_data), target, intent(inout) :: orbs
       type(input_variables):: input
       real(wp), dimension(:), pointer :: psi,hpsi,psit
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
       use module_base
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

     subroutine nonlocal_forces(iproc,lr,hx,hy,hz,at,rxyz,&
          orbs,nlpspd,proj,wfd,psi,fsep,refill)
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
       real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
       real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
       real(gp), dimension(3,at%nat), intent(inout) :: fsep
     END SUBROUTINE nonlocal_forces

     subroutine CalculateTailCorrection(iproc,nproc,at,rbuf,orbs,&
          Glr,nlpspd,ncongt,pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,&
          proj,psi,output_grid,ekin_sum,epot_sum,eproj_sum)
       use module_base
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
       use module_base
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
       use module_base
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
          hx,hy,hz,rxyz,rhopot,n3p,nlpspd,proj,pkernel,psi,v,ngatherarr,GPU)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,n1i,n2i
       integer, intent(in) :: nvirt,n3p
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: lr 
       type(orbitals_data), intent(in) :: orbs
       type(communications_arrays), intent(in) :: comms, commsv
       real(gp), intent(in) :: hx,hy,hz
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(*), intent(in) :: pkernel,rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       real(wp), dimension(:), pointer :: psi,v
     END SUBROUTINE davidson

     subroutine build_eigenvectors(iproc,norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinore,nspinor,&
          ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,nvirte,psivirt)
       use module_base
       implicit none
       !Arguments
       integer, intent(in) :: norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinor,ndim_hamovr,nspinore
       integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
       real(wp), dimension(nspin*ndim_hamovr), intent(in) :: hamovr
       real(wp), dimension(nvctrp,norbe), intent(in) :: psi
       real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: ppsit
       integer, dimension(2), intent(in), optional :: nvirte
       real(wp), dimension(*), optional :: psivirt
       integer:: iproc
     END SUBROUTINE build_eigenvectors

     subroutine preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
       use module_base
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
       use module_base
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
       use module_base
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
       use module_base
       use module_types
       implicit none
       character(len=10) :: comment
       character(len=11) :: orbname
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
       use module_base
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
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       type(orbitals_data), intent(in) :: orbs
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:), pointer :: coeffs
       logical, optional :: opt_fillrxyz
     END SUBROUTINE read_gaussian_information

     subroutine restart_from_gaussians(iproc,nproc,orbs,lr,hx,hy,hz,psi,G,coeffs)
       use module_base
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

     subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin,&
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks)
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
       type(orbitals_data), intent(out) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
     END SUBROUTINE inputguess_gaussian_orbitals

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
          nspin,eks,scorb,G,gaucoeff,iorbtolr)
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
       real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff
     END SUBROUTINE AtomicOrbitals

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
       use module_base
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
          ekin_sum,epot_sum,eproj_sum,nspin,GPU, in_iat_absorber, in )
       use module_base
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

       type(input_variables),intent(in) :: in
     END SUBROUTINE xabs_lanczos


     subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
          ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in  )! aggiunger a interface
       use module_base
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


       type(input_variables),intent(in) :: in

     END SUBROUTINE xabs_chebychev

     subroutine cg_spectra(iproc,nproc,at,hx,hy,hz,rxyz,&
          radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
          ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in  )! aggiunger a interface
       use module_base
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


       type(input_variables),intent(in) :: in

     END SUBROUTINE cg_spectra

     function GetBottom(  atoms, iproc)
       use module_base
       use module_types
       implicit none
       type(atoms_data), intent(in) :: atoms
       integer iproc
       real(gp) GetBottom
     end function GetBottom

     subroutine eleconf(nzatom,nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
       implicit none
       integer, intent(in) :: nzatom,nvalelec
       character(len=2), intent(out) :: symbol
       real(kind=8), intent(out) :: rcov,rprb,ehomo,amu
       integer, parameter :: nmax=6,lmax=3
       integer, intent(out) :: neleconf(nmax,0:lmax)
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
       use module_base
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
       use module_base
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
       use module_base
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
       use module_base
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
!!$       use module_base
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

     subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc)
       use module_base
       use module_types
       implicit none
       logical, intent(in) :: enlargerprb
       integer, intent(in) :: iproc,nspin,ng
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:), pointer :: Gocc
     END SUBROUTINE gaussian_pswf_basis

     subroutine local_analysis(iproc,nproc,hx,hy,hz,in,at,rxyz,shift,lr,orbs,orbsv,psi,psivirt)
       use module_base
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
       use module_base
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
       real(wp), dimension(*), intent(out), optional :: dvxcdrho
     END SUBROUTINE XC_potential

     subroutine direct_minimization(iproc,nproc,n1i,n2i,in,at,&
          orbs,orbsv,nvirt,lr,comms,commsv,&
          hx,hy,hz,rxyz,rhopot,n3p,nlpspd,proj, &
          pkernel,psi,psivirt,ngatherarr,GPU)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,n1i,n2i,nvirt,n3p
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: lr 
       type(orbitals_data), intent(in) :: orbs
       type(communications_arrays), intent(in) :: comms, commsv
       real(gp), intent(in) :: hx,hy,hz
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(*), intent(in) :: pkernel
       real(dp), dimension(*), intent(in), target :: rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       real(wp), dimension(:), pointer :: psi,psivirt
     END SUBROUTINE direct_minimization

     subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
          hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
       use module_base
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
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nat
       real(gp), dimension(nat), intent(in) :: radii
       real(gp), dimension(3,nat), target, intent(in) :: rxyz
       type(gaussian_basis), intent(out) :: G
     END SUBROUTINE gaussian_rism_basis

     subroutine gaussian_hermite_basis(nhermitemax,nat,radii,rxyz,G)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nat,nhermitemax
       real(gp), dimension(nat), intent(in) :: radii
       real(gp), dimension(3,nat), target, intent(in) :: rxyz
       type(gaussian_basis), intent(out) :: G  
     END SUBROUTINE gaussian_hermite_basis

    subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
      use module_base
      use module_types
      implicit none
      logical, intent(in) :: occorbs
      integer, intent(in) :: iproc,nspin,nvirt,nplot
      real(gp), intent(in) :: hx,hy,hz
      type(atoms_data), intent(in) :: at
      type(locreg_descriptors), intent(in) :: lr
      type(orbitals_data), intent(in) :: orbs,orbsv
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(:), pointer :: psi,psivirt
    END SUBROUTINE write_eigen_objects

    subroutine full_local_potential(iproc,nproc,ndimpot,ndimgrid,nspin,norb,norbp,ngatherarr,potential,pot)
      use module_base
      implicit none
      integer, intent(in) :: iproc,nproc,nspin,ndimpot,norb,norbp,ndimgrid
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      real(wp), dimension(max(ndimpot,1)*nspin), intent(in), target :: potential
      real(wp), dimension(:), pointer :: pot
    END SUBROUTINE full_local_potential

    subroutine free_full_potential(nproc,pot,subname)
      use module_base
      implicit none
      character(len=*), intent(in) :: subname
      integer, intent(in) :: nproc
      real(wp), dimension(:), pointer :: pot
    END SUBROUTINE free_full_potential

    subroutine getLocalizedBasis(iproc, nproc, nlr, Llr, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, &
        proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, trH, rxyzParabola, &
        itSCC, lastAlpha, infoBasisFunctions)
      use module_base
      use module_types
      implicit none
      integer:: iproc, nproc, nlr, idsxMin, idsxMax, infoBasisFunctions, itSCC
      type(locreg_descriptors),dimension(nlr),intent(in):: LLr
      type(atoms_data), intent(in) :: at
      type(orbitals_data):: orbs
      type(locreg_descriptors), intent(in) :: Glr
      type(input_variables):: input
      type(linearParameters):: lin
      real(8),dimension(3,at%nat):: rxyz
      integer:: nspin
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      real(dp), dimension(*), intent(inout) :: rhopot
      type(GPU_pointers), intent(in out) :: GPU
      real(dp), dimension(:), pointer :: pkernelseq
      real(8),dimension(lin%orbs%npsidim):: phi
      real(8):: trH, lastAlpha
      real(8),dimension(3,at%nat):: rxyzParabola
    end subroutine getLocalizedBasis


    subroutine getLocalizedBasisNew(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, &
        proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trH, rxyzParabola, coeff, &
        lastAlpha, infoBasisFunctions)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(atoms_data),intent(in) :: at
      type(orbitals_data),intent(in):: orbs
      type(locreg_descriptors),intent(in) :: Glr
      type(input_variables),intent(in):: input
      type(linearParameters),intent(in):: lin
      real(8),dimension(3,at%nat),intent(in):: rxyz, rxyzParabola
      integer,intent(in):: nspin
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp), dimension(nlpspd%nprojel),intent(in):: proj
      integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcs
      integer, dimension(0:nproc-1,2),intent(in):: ngatherarr  
      real(dp), dimension(*),intent(inout):: rhopot
      type(GPU_pointers),intent(inout):: GPU
      real(dp), dimension(:),pointer:: pkernelseq
      real(8),dimension(lin%orbs%npsidim),intent(inout):: phi
      real(8),dimension(lin%orbs%npsidim),intent(out):: hphi
      real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
      real(8),intent(out):: trH, lastAlpha
      integer,intent(out):: infoBasisFunctions
    end subroutine getLocalizedBasisNew



    subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, lind, phi, phid, &
          input, rxyz, occupForInguess, coeff, coeffd, nlr, Llr, outofzone)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc
      type(locreg_descriptors),intent(in):: Glr
      type(orbitals_data),intent(in):: orbs
      type(atoms_data),intent(in):: at
      type(linearParameters),intent(inout):: lin
      type(linearParameters),intent(inout):: lind
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(in):: rxyz
      real(8),dimension(32,at%nat):: occupForInguess
      real(8),dimension(:),allocatable,intent(out):: phi, phid
      real(8),dimension(:,:),allocatable,intent(out):: coeff, coeffd
      ! new
      integer,intent(out):: nlr
      type(locreg_descriptors),dimension(:),pointer,intent(out):: Llr
      integer,dimension(:,:),pointer,intent(out):: outofzone

    end subroutine allocateAndInitializeLinear



    subroutine transpose_vLIN(iproc, lproc, uproc, orbs, comms, psi, newComm, &
         work,outadd) !optional
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc, lproc, uproc, newComm
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      real(8),dimension(orbs%npsidim):: psi
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
      real(8),dimension(orbs%npsidim):: psi
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
      real(8),dimension(orbs%npsidim):: psi
      real(dp), dimension(:), pointer :: pkernel,pkernelseq
      integer, intent(in) :: potshortcut
      integer, dimension(*), intent(in) :: irrzon
      real(dp), dimension(*), intent(in) :: phnons
    END SUBROUTINE inputOrbitals
    
    
    
    
    subroutine diisstp(iproc,nproc,orbs,comms,diis,psit,quiet)
      use module_base
      use module_types
      implicit none
    ! Arguments
      integer, intent(in) :: nproc,iproc
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
      logical, optional:: quiet ! to avoid that the DIIS weights are written
    end subroutine diisstp
    
    
    subroutine psimix(iproc,nproc,orbs,comms,diis,hpsit,psit, quiet)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
      logical, optional:: quiet ! to avoid that the DIIS weights are written
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
    real(8),dimension(orbsLIN%npsidim):: phi
    real(8),dimension(3,at%nat):: perturbation
    end subroutine estimatePerturbedOrbitals
    
    subroutine psimixVariable(iproc,nproc,orbs,comms,diis,diisArr, hpsit,psit, quiet)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
      real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
      logical, optional:: quiet ! to avoid that the DIIS weights are written
    end subroutine psimixVariable
    
    
    
    subroutine diisstpVariable(iproc,nproc,orbs,comms,diis,diisArr,psit,quiet)
      use module_base
      use module_types
      implicit none
    ! Arguments
      integer, intent(in) :: nproc,iproc
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(diis_objects), intent(inout) :: diis
      type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
      real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
      logical, optional:: quiet ! to avoid that the DIIS weights are written
    end subroutine diisstpVariable
    
    
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


    subroutine getLinearPsi(iproc, nproc, nspin, nlr, Llr, Glr, orbs, comms, at, lin, lind, rxyz, rxyzParab, &
        nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, phid, psi, psit, &
        infoBasisFunctions, infoCoeff, itSCC, n3p, n3pi, n3d, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
        i3s, i3xcsh, fion, fdisp, fxyz, eion, edisp, fnoise, ebsMod, coeff, coeffd)
      use module_base
      use module_types
      !use Poisson_Solver
      implicit none
      integer,intent(in):: iproc, nproc, nspin, nlr, n3p, n3pi, n3d, i3s, i3xcsh, itSCC
      type(locreg_descriptors),dimension(nlr),intent(in):: LLr
      type(locreg_descriptors),intent(in):: Glr
      type(orbitals_data),intent(in) :: orbs
      type(communications_arrays),intent(in) :: comms
      type(atoms_data),intent(in):: at
      type(linearParameters),intent(in):: lin, lind
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(in):: rxyz, fion, fdisp
      real(8),dimension(3,at%nat),intent(inout):: rxyzParab
      integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
      type(GPU_pointers),intent(inout):: GPU
      integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
      real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
      real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
      real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
      !real(wp), dimension(lin%as%size_rhocore):: rhocore 
      real(wp), dimension(:),pointer,intent(in):: rhocore
      real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)),intent(inout):: potxc
      real(dp),dimension(:),pointer,intent(in):: pkernelseq
      real(8),dimension(lin%orbs%npsidim),intent(inout):: phi
      real(8),dimension(lind%orbs%npsidim),intent(inout):: phid
      real(8),dimension(orbs%npsidim),intent(out):: psi, psit
      integer,intent(out):: infoBasisFunctions, infoCoeff
      character(len=3),intent(in):: PSquiet
      real(8),intent(out):: ebsMod
      real(8),dimension(lin%orbs%norb,orbs%norb),intent(in out):: coeff
      real(8),dimension(lind%orbs%norb,orbs%norb),intent(in out):: coeffd
      real(8),dimension(3,at%nat),intent(out):: fxyz
      real(8):: eion, edisp, fnoise
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

    
    
    subroutine deallocateLinear(iproc, lin, lind, phi, coeff, phid, coeffd)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc
      type(linearParameters),intent(inout):: lin, lind
      real(8),dimension(:),allocatable,intent(inout):: phi, phid
      real(8),dimension(:,:),allocatable,intent(inout):: coeff, coeffd

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
    
    subroutine linearScaling(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, comms, at, input, lin, rxyz, &
        fion, fdisp, radii_cf, nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, pkernelseq, irrzon, phnons, &
        pkernel, pot_ion, rhocore, potxc, PSquiet, eion, edisp, eexctX, scpot, psi, psit, energy, fxyz)
      use module_base
      use module_types
      implicit none
      integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
      type(locreg_descriptors),intent(in) :: Glr
      type(orbitals_data),intent(in):: orbs
      type(communications_arrays),intent(in) :: comms
      type(atoms_data),intent(inout):: at
      type(linearParameters),intent(in out):: lin
      type(input_variables),intent(in):: input
      real(8),dimension(3,at%nat),intent(inout):: rxyz
      real(8),dimension(3,at%nat),intent(in):: fion, fdisp
      real(8),dimension(at%ntypes,3),intent(in):: radii_cf
      integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      !integer,dimension(0:nproc-1,2),intent(in):: ngatherarr
      integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
      type(nonlocal_psp_descriptors),intent(in):: nlpspd
      real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
      real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in out):: rhopot
      type(GPU_pointers),intent(in out):: GPU
      real(dp),dimension(:),pointer,intent(in):: pkernelseq
      integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
      real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
      real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
      real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
      !real(wp), dimension(lin%as%size_rhocore):: rhocore 
      real(wp), dimension(:),pointer,intent(in):: rhocore
      real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)),intent(inout):: potxc
      character(len=3),intent(in):: PSquiet
      real(gp),intent(in):: eion, edisp, eexctX
      logical,intent(in):: scpot
      real(8),dimension(orbs%npsidim),intent(out):: psi
      real(8),dimension(:),pointer,intent(out):: psit
      real(8),intent(out):: energy
      real(8),dimension(3,at%nat),intent(out):: fxyz
      !real(8),intent(out):: fnoise
      real(8):: fnoise
    end subroutine linearScaling
    
    
    subroutine potentialAndEnergySub(iproc, nproc, n3d, n3p, nlr, Llr, Glr, orbs, atoms, in, lin, phi, psi, rxyz, rxyzParab, &
        rhopot, nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
        proj, nlpspd, pkernelseq, eion, edisp, eexctX, scpot, coeff, ebsMod, energy)
      use module_base
      use module_types
      implicit none
      integer:: iproc, nproc, n3d, n3p, nlr
      type(locreg_descriptors),dimension(nlr),intent(in):: LLr
      type(locreg_descriptors) :: Glr
      type(orbitals_data):: orbs
      type(atoms_data):: atoms
      type(input_variables):: in
      type(linearParameters):: lin
      real(8),dimension(lin%orbs%npsidim):: phi
      real(8),dimension(orbs%npsidim):: psi
      real(dp), dimension(lin%as%size_rhopot) :: rhopot
      integer,dimension(0:nproc-1,4) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
      type(GPU_pointers),intent(in out):: GPU
      integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)) :: irrzon
      real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)) :: phnons
      real(dp), dimension(lin%as%size_pkernel):: pkernel
      real(wp), dimension(lin%as%size_pot_ion):: pot_ion
      real(wp), dimension(:),pointer:: rhocore
      real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)):: potxc
      character(len=3):: PSquiet
      type(nonlocal_psp_descriptors),intent(in) :: nlpspd
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(dp),dimension(lin%as%size_pkernelseq),intent(in):: pkernelseq
      real(8),dimension(3,atoms%nat),intent(in):: rxyz
      real(8),dimension(3,atoms%nat),intent(in):: rxyzParab
      real(gp):: eion, edisp, eexctX, energy
      real(8),dimension(lin%orbs%norb,orbs%norb):: coeff
      real(8):: ebsMod
      logical:: scpot
    end subroutine potentialAndEnergySub
    
    

    subroutine calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, atoms, in, comms, lin, nlpspd, proj, &
        ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, phi, coeff, fxyz, fnoise)
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
      real(8),dimension(orbs%npsidim),intent(inout):: psi
      real(8),dimension(lin%orbs%npsidim),intent(inout):: phi
      real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
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



    subroutine select_active_space(iproc,nproc,orbs,comms,mask_array,Glr,orbs_as,comms_as,psi,psi_as)
      use module_base
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
         ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp,psi,psit,hpsi,gnrm,gnrm_zero,energy)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,ncong,iscf,iter
      real(gp), intent(in) :: hx,hy,hz,ekin,epot,eproj,ehart,exc,evxc,eexctX,eion,edisp
      type(orbitals_data), intent(in) :: orbs
      type(communications_arrays), intent(in) :: comms
      type(locreg_descriptors), intent(in) :: lr
      type(GPU_pointers), intent(in) :: GPU
      real(gp), intent(out) :: gnrm,gnrm_zero,energy
      real(wp), dimension(:), pointer :: psi,psit,hpsi
    end subroutine calculate_energy_and_gradient



   subroutine determine_locreg_periodic(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr,outofzone)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nlr
      real(gp), intent(in) :: hx,hy,hz
      type(locreg_descriptors), intent(in) :: Glr
      real(gp), dimension(nlr), intent(in) :: locrad
      real(gp), dimension(3,nlr), intent(in) :: cxyz
      type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
      integer, dimension(3,nlr),intent(out) :: outofzone
   end subroutine

    subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr,outofzone)
      use module_base
      use module_types
      implicit none
      integer,intent(in) :: ilr,nlr
      type(locreg_descriptors),intent(in) :: Glr
      type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr
      integer,dimension(3,nlr),intent(in) :: outofzone
    end subroutine

    subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
     implicit none
     integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
     integer, dimension(nseg), intent(in) :: keyv
     integer, dimension(2,nseg), intent(in) :: keyg
     integer, intent(out) :: nseg_loc,nvctr_loc
     integer, dimension(3),intent(in) :: outofzone
    end subroutine

    subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc,outofzone)
     implicit none
     integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
     integer, dimension(nseg), intent(in) :: keyv
     integer, dimension(2,nseg), intent(in) :: keyg
     integer, dimension(3), intent(in) :: outofzone
     integer, dimension(nseg_loc), intent(out) :: keyv_loc
     integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
     end subroutine

    subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr,outofzone)
     use module_base
     use module_types
     implicit none
     integer, intent(in) :: alr,blr
     integer, intent(in) :: nlr
     type(locreg_descriptors),intent(in) :: Glr
     integer, intent(out) :: isovrlp
     integer,dimension(3,nlr),intent(in) :: outofzone
     type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
    end subroutine

    subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr,outofzone)
     use module_base
     use module_types
     implicit none
     integer, intent(in) :: alr,blr
     integer, intent(in) :: nlr
     type(locreg_descriptors),intent(in) :: Glr
     integer, intent(in) :: isovrlp
     type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
     type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr
     integer,dimension(3,nlr),intent(in) :: outofzone
    end subroutine

    subroutine nlpspd_to_locreg(input_parameters,iproc,Glr,Llr,rxyz,atoms,orbs,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,Lnlpspd,projflg)
     use module_base
     use module_types
     implicit none
     type(input_variables),intent(in) :: input_parameters
     integer,intent(in) :: iproc
     type(locreg_descriptors),intent(in) :: Glr
     type(locreg_descriptors),intent(in) :: Llr
     type(atoms_data),intent(in) :: atoms
     type(orbitals_data),intent(in) :: orbs
     real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
     type(nonlocal_psp_descriptors),intent(in) :: nlpspd
     type(nonlocal_psp_descriptors),intent(out) :: Lnlpspd
     integer,dimension(atoms%nat),intent(out) :: projflg
     real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
     real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
    end subroutine

    subroutine apply_local_projectors(atoms,in,Llr,Lnlpspd,Lproj,orbs,projflg,psi,rxyz,hpsi)
     use module_base
     use module_types
     implicit none
     type(atoms_data),intent(in) :: atoms
     type(input_variables),intent(in) :: in
     type(locreg_descriptors),intent(in) :: Llr
     type(nonlocal_psp_descriptors),intent(in) :: Lnlpspd
     type(orbitals_data),intent(in) :: orbs
     integer,dimension(atoms%nat),intent(in) :: projflg
     real(wp),dimension(Lnlpspd%nprojel),intent(out):: Lproj
     real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(in) :: psi
     real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),intent(out):: hpsi
     real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
    end subroutine

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
     real(wp),dimension(orbs%npsidim),intent(in) :: psi
     real(wp),dimension(ldim),intent(inout) :: lpsi
    end subroutine





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

    subroutine local_partial_densityLinear(iproc,nproc,nlr,rsflag,nscatterarr,&
         nrhotot,Glr,Llr,nrho,rho,hxh,hyh,hzh,nspin,orbs,psi)
      use module_base
      use module_types
      use libxc_functionals
      implicit none
      logical, intent(in) :: rsflag
      integer, intent(in) :: iproc,nproc,nlr,nrho
      integer,intent(inout):: nrhotot
      integer, intent(in) :: nspin
      real(dp),dimension(max(nrho,1),nspin):: rho
      real(gp), intent(in) :: hxh,hyh,hzh
      type(orbitals_data), intent(in) :: orbs
      type(locreg_descriptors), intent(in) :: Glr
      type(locreg_descriptors),dimension(nlr),intent(in) :: Llr
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      real(wp), dimension(orbs%npsidim), intent(in) :: psi
    end subroutine local_partial_densityLinear


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
         comms, Glr, input, lin, rxyz, n3p, rhopot, rhocore, pot_ion,&
         nlpspd, proj, pkernel, pkernelseq, &
         nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, &
         phi)
      use module_base
      use module_types
      implicit none
      integer, intent(in) :: iproc,nproc,n3p
      type(atoms_data), intent(inout) :: at
      type(nonlocal_psp_descriptors), intent(in) :: nlpspd
      type(locreg_descriptors), intent(in) :: Glr
      type(communications_arrays), intent(in) :: comms
      type(GPU_pointers), intent(inout) :: GPU
      type(input_variables):: input
      type(linearParameters),intent(inout):: lin
      integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
      real(gp), dimension(3,at%nat), intent(in) :: rxyz
      real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
      real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
      real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
      real(wp), dimension(:), pointer :: rhocore
      real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
      real(dp), dimension(:), pointer :: pkernelseq
      integer, intent(in) ::potshortcut
      integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
      real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
      real(8),dimension(lin%orbs%npsidim),intent(out):: phi
    end subroutine inputguessConfinement


    
  end interface

end module module_interfaces


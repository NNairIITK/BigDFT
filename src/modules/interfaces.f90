!!****m* BigDFT/interfaces
!! FUNCTION
!!  Modules which contains all interfaces
!!
!! DESCRIPTION
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
!!
!! AUTHOR
!!    Luigi Genovese, Damien Caliste
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
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
     end subroutine call_bigdft


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

     subroutine read_input_variables(iproc,posinp, file_dft, file_kpt, file_geopt, file_perf, inputs,atoms,rxyz)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: posinp
       character(len=*), intent(in) :: file_dft, file_geopt, file_kpt, file_perf
       integer, intent(in) :: iproc
       type(input_variables), intent(out) :: inputs
       type(atoms_data), intent(out) :: atoms
       real(gp), dimension(:,:), pointer :: rxyz
     END SUBROUTINE read_input_variables

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

     subroutine MemoryEstimator(geocode,nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,nat,ntypes,&
          iatype,rxyz,radii_cf,crmult,frmult,norb,nkpt,nprojel,atomnames,output_grid,nspin,peakmem)
       use module_base
       implicit none
       !Arguments
       character(len=1), intent(in) :: geocode
       integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin,nprojel,output_grid,nkpt
       integer, dimension(nat), intent(in) :: iatype
       character(len=20), dimension(ntypes), intent(in) :: atomnames
       real(kind=8), intent(in) :: hx,hy,hz,crmult,frmult,alat1,alat2,alat3
       real(kind=8), dimension(3,nat), intent(in) :: rxyz
       real(kind=8), dimension(ntypes,3), intent(in) ::  radii_cf
       real(kind=8), intent(out) :: peakmem
     END SUBROUTINE MemoryEstimator

     subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
          crmult,frmult,Glr)
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
     END SUBROUTINE createWavefunctionsDescriptors

     subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
          radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
       use module_base
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

     subroutine import_gaussians(iproc,nproc,at,orbs,comms,&
          Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,nlpspd,proj,& 
          pkernel,ixc,psi,psit,hpsi,nscatterarr,ngatherarr,nspin,symObj,irrzon,phnons)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ixc,nspin,symObj
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       type(nonlocal_psp_descriptors), intent(in) :: nlpspd
       type(locreg_descriptors), intent(in) :: Glr
       type(communications_arrays), intent(in) :: comms
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(*), intent(in) :: pkernel
       type(orbitals_data), intent(inout) :: orbs
       real(dp), dimension(*), intent(inout) :: rhopot
       real(wp), dimension(*), intent(inout) :: pot_ion
       real(wp), dimension(:), pointer :: psi,psit,hpsi,rhocore
       integer, dimension(*), intent(in) :: irrzon
       real(dp), dimension(*), intent(in) :: phnons
     END SUBROUTINE import_gaussians

     subroutine input_wf_diag(iproc,nproc,at,&
          orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
          nlpspd,proj,pkernel,ixc,psi,hpsi,psit,G,&
          nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU)
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
       integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
       real(dp), dimension(*), intent(in) :: pkernel
       real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
       type(gaussian_basis), intent(out) :: G 
       real(wp), dimension(:), pointer :: psi,hpsi,psit,rhocore
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

     subroutine first_orthon(iproc,nproc,orbs,wfd,comms,psi,hpsi,psit)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc
       type(orbitals_data), intent(in) :: orbs
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(communications_arrays), intent(in) :: comms
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
          nlpspd,proj,lr,ngatherarr,ndimpot,potential,psi,hpsi,&
          ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel,orbsocc,psirocc)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ndimpot,nspin
       real(gp), intent(in) :: hx,hy,hz
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
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
       real(dp), dimension(*), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
     END SUBROUTINE HamiltonianApplication

     subroutine hpsitopsi(iproc,nproc,orbs,hx,hy,hz,lr,comms,&
          ncong,iter,idsx,idsx_actual,ads,energy,energy_old,energy_min,&
          alpha,gnrm,scprsum,psi,psit,hpsi,psidst,hpsidst,nspin,GPU)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ncong,idsx,iter,nspin
       real(gp), intent(in) :: hx,hy,hz,energy,energy_old
       type(locreg_descriptors), intent(in) :: lr
       type(communications_arrays), intent(in) :: comms
       type(orbitals_data), intent(in) :: orbs
       integer, intent(inout) :: idsx_actual
       real(wp), intent(inout) :: alpha
       real(dp), intent(inout) :: gnrm,scprsum
       real(gp), intent(inout) :: energy_min
       real(wp), dimension(:), pointer :: psi,psit,hpsi,psidst,hpsidst
       real(wp), dimension(:,:,:), pointer :: ads
       type(GPU_pointers), intent(inout) :: GPU
     END SUBROUTINE hpsitopsi

     subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,wfd,comms,&
          psi,hpsi,psit,& !mandatory
          orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,natsc,nspin
       type(wavefunctions_descriptors), intent(in) :: wfd
       type(communications_arrays), target, intent(in) :: comms
       type(orbitals_data), target, intent(inout) :: orbs
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

     subroutine nonlocal_forces(iproc,n1,n2,n3,hx,hy,hz,at,rxyz,&
          orbs,nlpspd,proj,wfd,psi,fsep,refill)
       use module_base
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
       use module_base
       use module_types
       implicit none
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(locreg_descriptors), intent(in) :: Glr
       type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
       integer, intent(in) :: iproc,nproc,ncongt,nspin,output_grid
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
          psi_c,psi_f,norb,eval)
       use module_base
       implicit none
       logical, intent(in) :: useFormattedOutput
       integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,norb
       real(gp), intent(in) :: hx,hy,hz
       integer, dimension(nseg_c), intent(in) :: keyv_c
       integer, dimension(nseg_f), intent(in) :: keyv_f
       integer, dimension(2,nseg_c), intent(in) :: keyg_c
       integer, dimension(2,nseg_f), intent(in) :: keyg_f
       real(wp), dimension(norb), intent(in) :: eval
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

     subroutine build_eigenvectors(norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinore,nspinor,&
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
     END SUBROUTINE build_eigenvectors

     subroutine preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ncong
       real(gp), intent(in) :: hx,hy,hz
       type(locreg_descriptors), intent(in) :: lr
       type(orbitals_data), intent(in) :: orbs
       real(dp), intent(out) :: gnrm
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
     end subroutine plot_wf

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
       type(atoms_data), intent(in) :: at
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
       real(wp), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
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

     subroutine psimix(iproc,nproc,orbs,comms,ads,ids,mids,idsx,energy,energy_old,alpha,&
          hpsit,psidst,hpsidst,psit)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,ids,mids,idsx
       real(gp), intent(in) :: energy,energy_old
       type(orbitals_data), intent(in) :: orbs
       type(communications_arrays), intent(in) :: comms
       real(gp), intent(inout) :: alpha
       real(wp), dimension(:), pointer :: psit,hpsit,psidst,hpsidst
       real(wp), dimension(:,:,:), pointer :: ads
     END SUBROUTINE psimix

     subroutine plot_density(geocode,filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
          hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
       use module_base
       use module_types
       implicit none
       character(len=1), intent(in) :: geocode
       character(len=*), intent(in) :: filename
       integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nspin,nproc
       real(gp), intent(in) :: hxh,hyh,hzh
       type(atoms_data), intent(in) :: at
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(gp), dimension(3,at%nat), intent(in) :: rxyz
       real(dp), dimension(n1i*n2i*n3p,nspin), target, intent(in) :: rho
     end subroutine plot_density


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

     subroutine read_cube(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
          nat,rxyz)
       use module_base
       use module_types
       implicit none
       character(len=*), intent(in) :: filename
       character(len=1), intent(in) :: geocode
       integer, intent(in) :: nspin
       integer, intent(out) ::  n1i,n2i,n3i
       real(gp), intent(out) :: hxh,hyh,hzh
       real(dp), dimension(:,:), pointer :: rho
       real(gp), dimension(:,:), pointer, optional :: rxyz
       integer, intent(out), optional ::  nat
     end subroutine read_cube

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
          rho,exc,vxc,nspin,rhocore,potxc)
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
     end subroutine XC_potential

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
       real(dp), dimension(*), intent(in) :: pkernel,rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       real(wp), dimension(:), pointer :: psi,psivirt
     end subroutine direct_minimization

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
     end subroutine CounterIonPotential

     subroutine gaussian_rism_basis(nat,radii,rxyz,G)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nat
       real(gp), dimension(nat), intent(in) :: radii
       real(gp), dimension(3,nat), target, intent(in) :: rxyz
       type(gaussian_basis), intent(out) :: G
     end subroutine gaussian_rism_basis

     subroutine gaussian_hermite_basis(nhermitemax,nat,radii,rxyz,G)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nat,nhermitemax
       real(gp), dimension(nat), intent(in) :: radii
       real(gp), dimension(3,nat), target, intent(in) :: rxyz
       type(gaussian_basis), intent(out) :: G  
     end subroutine gaussian_hermite_basis

  end interface

end module module_interfaces
!!***

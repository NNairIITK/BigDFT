!!****m* BigDFT/interfaces
!! NAME
!!   interfaces
!!
!! FUNCTION
!!  Modules which contains all interfaces
!!
!! DESCRIPTION
!!  Interfaces of:
!!  - call_cluster
!!  - conjgrad
!!  - copy_old_wavefunctions
!!  - read_system_variables
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
!!    Copyright (C) 2007 CEA
!!
!! SOURCE
!!
module module_interfaces

implicit none

interface

   subroutine call_cluster(nproc,iproc,atoms,rxyz,energy,fxyz,&
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
     use module_types
     implicit none
     integer, intent(in) :: iproc,nproc
     type(input_variables),intent(inout) :: in
     type(wavefunctions_descriptors), intent(inout) :: wfd
     type(atoms_data), intent(inout) :: atoms
     integer, intent(inout) :: infocode,n1,n2,n3,norbp,norb
     real(kind=8), intent(out) :: energy
     real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz,rxyz_old
     real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz
     real(kind=8), dimension(:), pointer :: eval,psi
   end subroutine call_cluster

   subroutine conjgrad(nproc,iproc,at,wpos,etot,fxyz, &
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,in)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     integer, intent(in) :: nproc,iproc
     integer, intent(inout) :: n1,n2,n3,ncount_cluster,norbp,norb
     real(kind=8), intent(out) :: etot
     type(input_variables), intent(inout) :: in
     type(wavefunctions_descriptors), intent(inout) :: wfd
     real(kind=8), dimension(3,at%nat), intent(inout) :: wpos,rxyz_old
     real(kind=8), dimension(3,at%nat), intent(out) :: fxyz
     real(kind=8), dimension(:), pointer :: eval,psi
   end subroutine conjgrad

   subroutine copy_old_wavefunctions(iproc,nproc,norb,norbp,nspinor,hx,hy,hz,n1,n2,n3,wfd,psi,&
        hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,wfd_old,psi_old)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors) :: wfd,wfd_old
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nspinor
     real(gp), intent(in) :: hx,hy,hz
     integer, intent(out) :: n1_old,n2_old,n3_old
     real(gp), intent(out) :: hx_old,hy_old,hz_old
     real(wp), dimension(:), pointer :: psi,psi_old
   end subroutine copy_old_wavefunctions

   subroutine read_system_variables(iproc,nproc,in,at,radii_cf,nelec,&
        norb,norbu,norbd,norbp,iunit)
     use module_types
     implicit none
     type(input_variables), intent(in) :: in
     integer, intent(in) :: iproc,nproc
     type(atoms_data), intent(inout) :: at
     integer, intent(out) :: nelec,norb,norbu,norbd,norbp,iunit
     real(kind=8), dimension(at%ntypes,2), intent(out) :: radii_cf
   end subroutine read_system_variables

   subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)
     implicit none
     ! Arguments
     integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
     real(kind=8), intent(out) :: occup(norb),spinsgn(norb)
   end subroutine input_occup

   subroutine system_size(iproc,geocode,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,&
        alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: atoms
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc
     real(kind=8), intent(in) :: crmult,frmult
     real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
     real(kind=8), dimension(atoms%ntypes,2), intent(in) :: radii_cf
     integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
     real(kind=8), intent(inout) :: hx,hy,hz,alat1,alat2,alat3
   end subroutine system_size

   subroutine MemoryEstimator(geocode,nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hx,hy,hz,nat,ntypes,&
        iatype,rxyz,radii_cf,crmult,frmult,norb,nprojel,atomnames,output_grid,nspin,peakmem)
     implicit none
     !Arguments
     character(len=1), intent(in) :: geocode
     logical, intent(in) :: output_grid
     integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin,nprojel
     integer, dimension(nat), intent(in) :: iatype
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     real(kind=8), intent(in) :: hx,hy,hz,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
     real(kind=8), intent(out) :: peakmem
   end subroutine MemoryEstimator

   subroutine createWavefunctionsDescriptors(iproc,nproc,geocode,n1,n2,n3,output_grid,&
        hx,hy,hz,atoms,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
        wfd,nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,nspinor)
     use module_types
     implicit none
     !Arguments
     type(atoms_data), intent(in) :: atoms
     character(len=1), intent(in) :: geocode
     logical, intent(in) :: output_grid
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspinor
     real(kind=8), intent(in) :: hx,hy,hz,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,atoms%nat), intent(in) :: rxyz
     real(kind=8), dimension(atoms%ntypes,2), intent(in) :: radii_cf
     type(wavefunctions_descriptors) , intent(out) :: wfd
     !boundary arrays
     type(convolutions_bounds), intent(out) :: bounds
     integer, intent(out) :: nvctrp
   end subroutine createWavefunctionsDescriptors

   subroutine createProjectorsArrays(geocode,iproc,n1,n2,n3,rxyz,at,&
        radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,n1,n2,n3
     real(kind=8), intent(in) :: cpmult,fpmult,hx,hy,hz
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(at%ntypes,2), intent(in) :: radii_cf
     type(nonlocal_psp_descriptors), intent(out) :: nlpspd
     real(kind=8), dimension(:), pointer :: proj
   end subroutine createProjectorsArrays

   subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1i,n2i,n3i,ixc,&
        n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)
     implicit none
     character(len=1), intent(in) :: geocode,datacode
     integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,ixc
     integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
     integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
     integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
   end subroutine createDensPotDescriptors

   subroutine IonicEnergyandForces(geocode,iproc,nproc,at,hxh,hyh,hzh,alat1,alat2,alat3,rxyz,eion,fion,psoffset,&
        n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi
     real(kind=8), intent(in) :: alat1,alat2,alat3,hxh,hyh,hzh
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), intent(out) :: eion,psoffset
     real(kind=8), dimension(3,at%nat), intent(out) :: fion
     real(kind=8), dimension(*), intent(out) :: pot_ion
   end subroutine IonicEnergyandForces

   subroutine createIonicPotential(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,&
        hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,eion,psoffset)
     implicit none
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
     real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield,psoffset
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), intent(out) :: eion
     real(kind=8), dimension(*), intent(inout) :: pot_ion
   end subroutine createIonicPotential

   subroutine import_gaussians(geocode,iproc,nproc,cpmult,fpmult,radii_cf,at,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
        norb,norbp,occup,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj,& 
        pkernel,ixc,psi,psit,hpsi,eval,nscatterarr,ngatherarr,nspin,spinsgn)
     use module_base
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin
     real(kind=8), intent(in) :: hx,hy,hz,cpmult,fpmult
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: spinsgn,occup
     real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf  
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:), pointer :: psi,psit,hpsi
   end subroutine import_gaussians

   subroutine input_wf_diag(geocode,iproc,nproc,cpmult,fpmult,radii_cf,at,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        norb,norbp,nvirte,nvirtep,nvirt,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,&
        wfd,bounds,nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,eval,&
        nscatterarr,ngatherarr,nspin,spinsgn)
     use module_base
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
     integer, intent(inout) :: nspin,nvirte,nvirtep,nvirt
     real(kind=8), intent(in) :: hx,hy,hz,cpmult,fpmult
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: spinsgn
     real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf  
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:), pointer :: psi,hpsi,psit,psivirt
   end subroutine input_wf_diag

   subroutine reformatmywaves(iproc,norb,norbp,nat,&
        hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
        hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
     integer, intent(in) :: iproc,norb,norbp,nat,n1_old,n2_old,n3_old,n1,n2,n3
     real(gp), intent(in) :: hx_old,hy_old,hz_old,hx,hy,hz
     real(gp), dimension(3,nat), intent(in) :: rxyz,rxyz_old
     real(wp), dimension(wfd_old%nvctr_c+7*wfd_old%nvctr_f,norbp), intent(in) :: psi_old
     real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: psi
   end subroutine reformatmywaves

   subroutine first_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,&
        nspin,psi,hpsi,psit)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctrp,nspin
     real(wp), dimension(:) , pointer :: psi,hpsi,psit
   end subroutine first_orthon

   subroutine sumrho(geocode,iproc,nproc,norb,norbp,ixc,n1,n2,n3,hxh,hyh,hzh,occup,  & 
        wfd,psi,rho,nrho,nscatterarr,nspin,nspinor,spinsgn,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin,nspinor,ixc
     integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hxh,hyh,hzh
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(max(nrho,1),nspinor), intent(out), target :: rho
   end subroutine sumrho

   subroutine HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
        norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds,nlpspd,proj,&
        ngatherarr,ndimpot,potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,nspinor,spinsgn)
     use module_base
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ndimpot
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin,nspinor
     real(kind=8), intent(in) :: hx,hy,hz,cpmult,fpmult
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
     real(gp), dimension(3,at%nat), intent(in) :: rxyz
     real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf  
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: hpsi
   end subroutine HamiltonianApplication

   subroutine hpsitopsi(geocode,iter,iproc,nproc,norb,norbp,occup,hx,hy,hz,n1,n2,n3,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,kbounds,&
        eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
        psi,psit,hpsi,psidst,hpsidst,nspin,nspinor,spinsgn)
     use module_types
     implicit none
     type(kinetic_bounds), intent(in) :: kbounds
     type(wavefunctions_descriptors), intent(in) :: wfd
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin,nspinor
     real(kind=8), intent(in) :: hx,hy,hz,energy,energy_old
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
     real(kind=8), intent(inout) :: alpha
     real(kind=8), intent(inout) :: gnrm,scprsum
     real(kind=8), dimension(:), pointer :: psi,psit,hpsi,psidst,hpsidst
     real(kind=8), dimension(:,:,:), pointer :: ads
   end subroutine hpsitopsi

   subroutine DiagHam(iproc,nproc,natsc,nspin,nspinor,norbu,norbd,norb,norbp,nvctrp,wfd,&
        psi,hpsi,psit,eval,& !mandatory
        norbe,norbep,etol,norbsc_arr,nvirte,nvirtep,psivirt) !optional
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,natsc,nspin,nspinor,norb,norbu,norbd,norbp,nvctrp
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:), pointer :: psi,hpsi,psit
     !optional arguments
     integer, optional, intent(in) :: norbe,norbep,nvirte
     integer, optional, intent(out) :: nvirtep
     real(kind=8), optional, intent(in) :: etol
     integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
     real(wp), dimension(:), pointer, optional :: psivirt
   end subroutine DiagHam

   subroutine last_orthon(iproc,nproc,norbu,norbd,norb,norbp,wfd,nvctrp,&
        nspin,psi,hpsi,psit,occup,evsum,eval)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctrp,nspin
     real(gp), dimension(norb), intent(in) :: occup
     real(wp), intent(out) :: evsum
     real(wp), dimension(norb), intent(out) :: eval
     real(wp), dimension(:) , pointer :: psi,hpsi,psit
   end subroutine last_orthon

   subroutine local_forces(geocode,iproc,nproc,at,rxyz,hxh,hyh,hzh,&
        n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,rho,pot,floc)
     ! Calculates the local forces acting on the atoms belonging to iproc
     use module_types
     implicit none
     !Arguments---------
     type(atoms_data), intent(in) :: at
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
     real(kind=8), intent(in) :: hxh,hyh,hzh
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: rho,pot
     real(kind=8), dimension(3,at%nat), intent(out) :: floc
   end subroutine local_forces

   subroutine projectors_derivatives(geocode,iproc,at,n1,n2,n3,norb,&
        nlpspd,proj,rxyz,radii_cf,cpmult,fpmult,hx,hy,hz,derproj)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     character(len=1), intent(in) :: geocode
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     !Arguments-------------
     integer, intent(in) :: iproc,norb
     integer, intent(in) :: n1,n2,n3
     real(kind=8),intent(in) :: cpmult,fpmult,hx,hy,hz
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(at%ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(nlpspd%nprojel,3), intent(out) :: derproj
   end subroutine projectors_derivatives

   subroutine nonlocal_forces(geocode,iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,at,rxyz,radii_cf,&
        norb,norbp,nspinor,occup,nlpspd,proj,wfd,psi,fsep,refill)
     use module_base
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     character(len=1), intent(in) :: geocode
     logical, intent(in) :: refill
     integer, intent(in) :: iproc,norb,norbp,nspinor,n1,n2,n3
     real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult 
     real(gp), dimension(norb), intent(in) :: occup
     real(gp), dimension(3,at%nat), intent(in) :: rxyz
     real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf
     real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp*nspinor), intent(in) :: psi
     real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
     real(gp), dimension(3,at%nat), intent(inout) :: fsep
   end subroutine nonlocal_forces

   subroutine CalculateTailCorrection(iproc,nproc,at,n1,n2,n3,rbuf,norb,norbp,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
        pot,hgrid,rxyz,radii_cf,crmult,frmult,cpmult,fpmult,nspin,spinsgn,&
        proj,psi,occup,output_grid,ekin_sum,epot_sum,eproj_sum)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
     logical, intent(in) :: output_grid
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ncongt,nspin
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf,cpmult,fpmult
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
     real(kind=8), dimension(at%ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31,nspin), intent(in) :: pot
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
   end subroutine CalculateTailCorrection

   !added for abinit compatilbility
   subroutine reformatonewave(iproc,displ,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,nat,&
        & rxyz_old,psigold,hx,hy,hz,nvctr_c,nvctr_f,n1,n2,n3,rxyz,nseg_c,nseg_f,&
        & keyg,keyv,psifscf,psi)
     use module_base
     implicit none
     integer, intent(in) :: iproc,n1_old,n2_old,n3_old,nat,nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
     real(gp), intent(in) :: hx,hy,hz,displ,hx_old,hy_old,hz_old
     integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
     integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
     real(gp), dimension(3,nat), intent(in) :: rxyz_old,rxyz
     real(wp), dimension(0:n1_old,2,0:n2_old,2,0:n3_old,2), intent(in) :: psigold
     real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(out) :: psifscf
     real(wp), dimension(nvctr_c + 7 * nvctr_f), intent(out) :: psi
   end subroutine reformatonewave

   subroutine davidson(geocode,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,at,&
        cpmult,fpmult,radii_cf,&
        norb,norbu,norbp,nvirte,nvirtep,nvirt,gnrm_cv,nplot,n1,n2,n3,nvctrp,&
        hx,hy,hz,rxyz,rhopot,occup,i3xcsh,n3p,itermax,wfd,bounds,nlpspd,proj,  & 
        pkernel,ixc,psi,v,eval,ncong,nscatterarr,ngatherarr)
     use module_base
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc,n1i,n2i,n3i
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,i3xcsh,nvctrp,norbu
     integer, intent(in) :: nvirte,nvirtep,nvirt,ncong,n3p,itermax,nplot
     real(gp), dimension(norb), intent(in) :: occup
     real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf  
     real(dp), intent(in) :: gnrm_cv !convergence criterion for gradients
     real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(gp), dimension(3,at%nat), intent(in) :: rxyz
     real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
     real(dp), dimension(*), intent(in) :: pkernel,rhopot
     !this is a Fortran 95 standard, should be avoided (it is a pity IMHO)
     !real(kind=8), dimension(:,:,:,:), allocatable :: rhopot 
     real(wp), dimension(norb), intent(in) :: eval
     real(wp), dimension(:), pointer :: psi,v
   end subroutine davidson

   subroutine build_eigenvectors(nproc,norbu,norbd,norbp,norbep,nvctrp,nvctr,natsc,nspin,nspinor,&
        ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,nvirte,psivirt)
     use module_base
     implicit none
     !Arguments
     integer, intent(in) :: nproc,norbu,norbd,norbp,norbep,nvctrp,nvctr,natsc
     integer, intent(in) :: nspin,nspinor,ndim_hamovr
     integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
     real(kind=8), dimension(nspin*ndim_hamovr), intent(in) :: hamovr
     real(kind=8), dimension(nvctrp,norbep*nproc), intent(in) :: psi
     real(kind=8), dimension(nvctrp*nspinor,norbp*nproc), intent(out) :: ppsit
     integer, intent(in), optional :: nvirte
     real(wp), dimension(:), pointer, optional :: psivirt
   end subroutine build_eigenvectors

   subroutine preconditionall(geocode,iproc,nproc,norb,norbp,n1,n2,n3,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        hx,hy,hz,ncong,nspinor,wfd,eval,kb,hpsi,gnrm)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(kinetic_bounds), intent(in) :: kb
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     integer, intent(in) :: nspinor,ncong
     real(gp), intent(in) :: hx,hy,hz
     real(wp), dimension(norb), intent(in) :: eval
     real(dp), intent(out) :: gnrm
     real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp*nspinor), intent(inout) :: hpsi
   end subroutine preconditionall

   subroutine transpose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
        work,outadd) !optional
     use module_base
     use module_types
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
     real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norbp), intent(inout) :: psi
     real(wp), dimension(:), pointer, optional :: work
     real(wp), optional, intent(out) :: outadd !pass only the address to avoid pointer problems 
   end subroutine transpose

   subroutine untranspose(iproc,nproc,norb,norbp,nspinor,wfd,nvctrp,psi,&
        work,outadd) !optional
     use module_base
     use module_types
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,norb,norbp,nspinor,nvctrp
     real(wp), dimension(nspinor*nvctrp,norbp,nproc), intent(inout) :: psi
     real(wp), dimension(:), pointer, optional :: work
     real(wp), optional, intent(out) :: outadd !pass only the address to avoid pointer problems 
   end subroutine untranspose

   subroutine plot_wf(orbname,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,rx,ry,rz,wfd,&
        bounds,psi)
     use module_base
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=10) :: orbname 
     integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(gp), intent(in) :: hgrid,rx,ry,rz
     real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
   end subroutine plot_wf

   subroutine partial_density(rsflag,nproc,n1i,n2i,n3i,nspinor,nspinn,nrhotot,&
        hfac,nscatterarr,spinsgn,psir,rho_p,&
        ibyyzz_r) !optional argument
     use module_base
     implicit none
     logical, intent(in) :: rsflag
     integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinor,nspinn
     real(gp), intent(in) :: hfac,spinsgn
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
     real(wp), dimension(n1i,n2i,n3i,nspinn), intent(in) :: psir
     real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
     integer, dimension(:,:,:), pointer, optional :: ibyyzz_r 
   end subroutine partial_density

end interface

end module module_interfaces
!!***

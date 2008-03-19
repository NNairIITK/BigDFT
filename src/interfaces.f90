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

   subroutine call_cluster(parallel,nproc,iproc,atoms,rxyz,energy,fxyz,&
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
     use module_types
     implicit none
     logical, intent(in) :: parallel
     integer, intent(in) :: iproc,nproc
     type(input_variables),intent(inout) :: in
     type(wavefunctions_descriptors), intent(inout) :: wfd
     type(atoms_data), intent(inout) :: atoms
     integer, intent(inout) :: infocode,n1,n2,n3,norbp,norb
     real(kind=8), intent(out) :: energy
     real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
     real(kind=8), dimension(3,atoms%nat), intent(out) :: fxyz,rxyz_old
     real(kind=8), dimension(:), pointer :: eval
     real(kind=8), dimension(:,:), pointer :: psi
   end subroutine call_cluster

   subroutine conjgrad(parallel,nproc,iproc,at,wpos,etot,gg, &
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,in)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     logical, intent(in) :: parallel
     integer, intent(in) :: nproc,iproc
     integer, intent(inout) :: n1,n2,n3,ncount_cluster,norbp,norb
     real(kind=8), intent(out) :: etot
     type(input_variables), intent(inout) :: in
     type(wavefunctions_descriptors), intent(inout) :: wfd
     real(kind=8), dimension(3,at%nat), intent(inout) :: wpos
     real(kind=8), dimension(3,at%nat), intent(out) :: rxyz_old,gg
     real(kind=8), dimension(:), pointer :: eval
     real(kind=8), dimension(:,:), pointer :: psi
   end subroutine conjgrad

   subroutine copy_old_wavefunctions(iproc,nproc,norb,norbp,hgrid,n1,n2,n3,eval,wfd,psi,&
        hgrid_old,n1_old,n2_old,n3_old,eval_old,wfd_old,psi_old)
     use module_types     
     implicit none
     type(wavefunctions_descriptors) :: wfd,wfd_old
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3
     real(kind=8), intent(in) :: hgrid
     integer, intent(out) :: n1_old,n2_old,n3_old
     real(kind=8), intent(out) :: hgrid_old
     real(kind=8), dimension(:), pointer :: eval,eval_old
     real(kind=8), dimension(:,:), pointer :: psi,psi_old
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
        iatype,rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid,nspin,peakmem)
     implicit none
     !Arguments
     character(len=1), intent(in) :: geocode
     logical, intent(in) :: output_grid
     integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin
     integer, dimension(nat), intent(in) :: iatype
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     real(kind=8), intent(in) :: hx,hy,hz,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
     real(kind=8), intent(out) :: peakmem
   end subroutine MemoryEstimator

   subroutine createWavefunctionsDescriptors(iproc,nproc,geocode,n1,n2,n3,output_grid,&
        hx,hy,hz,atoms,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
        wfd,nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
     use module_types
     implicit none
     !Arguments
     type(atoms_data), intent(in) :: atoms
     character(len=1), intent(in) :: geocode
     logical, intent(in) :: output_grid
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
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

   subroutine createIonicPotential(geocode,iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,&
        hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,eion)
     implicit none
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,nat,ntypes,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
     real(kind=8), intent(in) :: hxh,hyh,hzh,elecfield
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), intent(out) :: eion
     real(kind=8), dimension(*), intent(out) :: pot_ion
   end subroutine createIonicPotential

   subroutine import_gaussians(geocode,iproc,nproc,at,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
        norb,norbp,occup,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj,& 
        pkernel,ixc,psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinsgn)
     use module_types
     use Poisson_Solver
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     character(len=1), intent(in) :: geocode,datacode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin
     real(kind=8), intent(in) :: hx,hy,hz
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: spinsgn,occup
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), intent(out) :: accurex
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
   end subroutine import_gaussians

   subroutine input_wf_diag(geocode,iproc,nproc,at,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        norb,norbp,n1,n2,n3,nvctrp,hx,hy,hz,rxyz,rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
        pkernel,ixc,psi,hpsi,psit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinsgn)
     ! Input wavefunctions are found by a diagonalization in a minimal basis set
     ! Each processors write its initial wavefunctions into the wavefunction file
     ! The files are then read by readwave
     use module_types
     use Poisson_Solver
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: datacode,geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
     integer, intent(in) :: nspin
     real(kind=8), intent(in) :: hx,hy,hz
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: spinsgn
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), intent(out) :: accurex
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:), pointer :: psi,hpsi,psit
   end subroutine input_wf_diag

   subroutine reformatmywaves(iproc,norb,norbp,nat,&
        & hgrid_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
        & hgrid,n1,n2,n3,rxyz,wfd,psi)
     use module_types
     implicit real(kind=8) (a-h,o-z)
     type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
     dimension :: rxyz(3,nat), rxyz_old(3,nat), center(3), center_old(3)
     dimension :: psi_old(wfd_old%nvctr_c + 7 * wfd_old%nvctr_f, norbp), &
          psi(wfd%nvctr_c + 7 * wfd%nvctr_f, norbp)
   end subroutine reformatmywaves

   subroutine first_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
        nspin,psi,hpsi,psit)
     implicit none
     logical, intent(in) :: parallel
     integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
     real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
   end subroutine first_orthon

   subroutine sumrho(geocode,iproc,nproc,norb,norbp,n1,n2,n3,hxh,hyh,hzh,occup,  & 
        wfd,psi,rho,nrho,nscatterarr,nspin,spinsgn,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin
     integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hxh,hyh,hzh
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(max(nrho,1),nspin), intent(out), target :: rho
   end subroutine sumrho

   subroutine HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,&
        norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds,nlpspd,proj,&
        ngatherarr,ndimpot,potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinsgn)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     character(len=1), intent(in) :: geocode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ndimpot
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin
     real(kind=8), intent(in) :: hx,hy,hz
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: occup,spinsgn
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(max(ndimpot,1),nspin), intent(in), target :: potential
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: hpsi
   end subroutine HamiltonianApplication

   subroutine hpsitopsi(iter,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,kbounds,&
        eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
        psi,psit,hpsi,psidst,hpsidst,nspin,spinsgn)
     use module_types
     implicit none
     type(kinetic_bounds), intent(in) :: kbounds
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin
     real(kind=8), intent(in) :: hgrid,energy,energy_old
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
     real(kind=8), intent(inout) :: alpha
     real(kind=8), intent(inout) :: gnrm,scprsum
     real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
     real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
   end subroutine hpsitopsi

   subroutine DiagHam(iproc,nproc,natsc,nspin,norbu,norbd,norb,norbp,nvctrp,wfd,&
        psi,hpsi,psit,eval,& !mandatory
        norbe,norbep,etol,norbsc_arr) !optional
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     integer, intent(in) :: iproc,nproc,natsc,nspin,norb,norbu,norbd,norbp,nvctrp
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:), pointer :: psi,hpsi,psit
     !optional arguments
     integer, optional, intent(in) :: norbe,norbep
     real(kind=8), optional, intent(in) :: etol
     integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   end subroutine DiagHam

   subroutine last_orthon(iproc,nproc,parallel,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,&
        nspin,psi,hpsi,psit,occup,evsum,eval)
     implicit none
     logical, intent(in) :: parallel
     integer, intent(in) :: iproc,nproc,norbu,norbd,norb,norbp,nvctr_c,nvctr_f,nvctrp,nspin
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), intent(out) :: evsum
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:) , pointer :: psi,hpsi,psit
   end subroutine last_orthon

   subroutine local_forces(iproc,nproc,at,rxyz,hgrid,n1,n2,n3,n3pi,i3s,rho,pot,floc)
     use module_types
     implicit none
     !Arguments---------
     type(atoms_data), intent(in) :: at
     integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s
     real(kind=8), intent(in) :: hgrid
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

   subroutine nonlocal_forces(iproc,at,norb,norbp,occup,nlpspd,proj,derproj,wfd,psi,fsep)
     use module_types
     implicit none
     !Arguments-------------
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     integer, intent(in) :: iproc,norb,norbp
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(nlpspd%nprojel,3), intent(in) :: derproj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(3,at%nat), intent(inout) :: fsep
   end subroutine nonlocal_forces

   subroutine CalculateTailCorrection(iproc,nproc,at,n1,n2,n3,rbuf,norb,norbp,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,eval,&
        pot,hgrid,rxyz,radii_cf,crmult,frmult,nspin,spinsgn,&
        proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)
     use module_types
     implicit none
     type(atoms_data), intent(in) :: at
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
     logical, intent(in) :: output_grid,parallel
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,ncongt,nspin
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinsgn
     real(kind=8), dimension(at%ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(3,at%nat), intent(in) :: rxyz
     real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31,nspin), intent(in) :: pot
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
   end subroutine CalculateTailCorrection

   !added for abinit compatilbility
   subroutine reformatonewave(iproc, hgrid_old, n1_old, n2_old, n3_old, &
        & center_old, psigold, hgrid, nvctr_c, nvctr_f, n1, n2, n3, center, nseg_c, nseg_f, &
        & keyg, keyv, psifscf, psi)
     implicit real(kind=8) (a-h,o-z)
     dimension :: center(3), center_old(3)
     dimension :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
     dimension :: psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2), psi(nvctr_c + 7 * nvctr_f)
     dimension :: psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
   end subroutine reformatonewave

end interface

end module module_interfaces
!!***

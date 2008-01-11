module module_interfaces

interface

   subroutine call_cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz,&
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,in,infocode)
     use module_types
     implicit none
     type(input_variables) :: in
     type(wavefunctions_descriptors) :: wfd
     logical, intent(in) :: parallel
     integer, intent(in) :: iproc,nproc,nat,ntypes,norbp,norb
     integer, intent(inout) :: infocode,n1,n2,n3
     real(kind=8), intent(out) :: energy
     character(len=20), dimension(100), intent(in) :: atomnames
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(3,nat), intent(in) :: rxyz_old
     real(kind=8), dimension(3,nat), intent(inout) :: rxyz
     real(kind=8), dimension(3,nat), intent(out) :: fxyz
     real(kind=8), dimension(:), pointer :: eval
     real(kind=8), dimension(:,:), pointer :: psi
   end subroutine call_cluster

   subroutine conjgrad(parallel,nproc,iproc,nat,ntypes,iatype,lfrztyp,atomnames,wpos,etot,gg, &
        psi,wfd,norbp,norb,eval,n1,n2,n3,rxyz_old,ncount_cluster,in)
     use module_types
     implicit none
     logical, intent(in) :: parallel
     integer, intent(in) :: nproc,iproc,nat,ntypes,norbp,norb
     integer, intent(inout) :: n1,n2,n3,ncount_cluster
     real(kind=8), intent(out) :: etot
     type(input_variables), intent(inout) :: in
     type(wavefunctions_descriptors), intent(inout) :: wfd
     character(len=20), dimension(100), intent(in) :: atomnames
     logical, dimension(ntypes), intent(in) :: lfrztyp
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(3,nat), intent(inout) :: wpos
     real(kind=8), dimension(3,nat), intent(out) :: rxyz_old,gg
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

   subroutine read_system_variables(iproc,nproc,nat,ntypes,nspin,ncharge,mpol,ixc,hgrid,atomnames,iatype,&
        psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,norb,norbu,norbd,norbp,iunit)
     implicit none
     integer, intent(in) :: iproc,nproc,nat,ntypes,nspin,ncharge,mpol,ixc
     real(kind=8), intent(in) :: hgrid
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     integer, dimension(nat), intent(in) :: iatype
     integer, intent(out) :: nelec,natsc,norb,norbu,norbd,norbp,iunit
     integer, dimension(ntypes), intent(out) :: npspcode,iasctype,nelpsp,nzatom
     real(kind=8), dimension(ntypes,2), intent(out) :: radii_cf
     real(kind=8), dimension(0:4,0:6,ntypes), intent(out) :: psppar
   end subroutine read_system_variables

   subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinar)
     implicit none
     ! Arguments
     integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
     real(kind=8), intent(out) :: occup(norb),spinar(norb)
   end subroutine input_occup

   subroutine system_size(iproc,nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
        alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)
     ! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
     !and shifts the atoms such that their position is the most symmetric possible
     implicit none
     integer, intent(in) :: iproc,nat,ntypes
     real(kind=8), intent(in) :: hgrid,crmult,frmult
     character(len=20), dimension(ntypes), intent(in) :: atomnames
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(3,nat), intent(inout) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
     real(kind=8), intent(out) :: alat1,alat2,alat3
   end subroutine system_size

   subroutine MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
        rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid,nspin,peakmem)
     implicit none
     !Arguments
     logical, intent(in) :: output_grid
     integer, intent(in) :: nproc,idsx,n1,n2,n3,nat,ntypes,norb,nspin
     integer, dimension(nat), intent(in) :: iatype
     character(len=20), dimension(100), intent(in) :: atomnames
     real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) ::  radii_cf
     real(kind=8), intent(out) :: peakmem
   end subroutine MemoryEstimator

   subroutine createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,&
        hgrid,nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
        wfd,nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
     use module_types
     implicit none
     !Arguments
     integer, intent(in) :: iproc,nproc,idsx,n1,n2,n3,nat,ntypes,norb,norbp
     integer, intent(out) :: nvctrp
     logical, intent(in) :: output_grid
     integer, intent(in) :: iatype(nat)
     real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
     real(kind=8) :: rxyz(3, nat), radii_cf(ntypes, 2)
     character(len=20), intent(in) :: atomnames(100)
     integer,intent(in):: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     type(wavefunctions_descriptors) , intent(out) :: wfd
     !boundary arrays
     type(convolutions_bounds), intent(out) :: bounds
   end subroutine createWavefunctionsDescriptors

   subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,nat,ntypes,iatype,atomnames,&
        & psppar,npspcode,radii_cf,cpmult,fpmult,hgrid,nlpspd,proj)
     use module_types
     implicit none
     type(nonlocal_psp_descriptors), intent(out) :: nlpspd
     character(len=20), dimension(100),intent(in) :: atomnames
     integer, intent(in) :: iproc,n1,n2,n3,nat,ntypes
     real(kind=8), intent(in) :: cpmult,fpmult,hgrid
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: npspcode
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(:), pointer :: proj
   end subroutine createProjectorsArrays

   subroutine createDensPotDescriptors(iproc,nproc,geocode,datacode,n1,n2,n3,ixc,&
        n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr)
     implicit none
     character(len=1), intent(in) :: geocode,datacode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,ixc
     integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
     integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
     integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
   end subroutine createDensPotDescriptors

   subroutine createIonicPotential(iproc,nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,hgrid,&
        elecfield,n1,n2,n3,n3pi,i3s,pkernel,pot_ion,eion)
     implicit none
     integer, intent(in) :: iproc,nproc,nat,ntypes,n1,n2,n3,n3pi,i3s
     real(kind=8), intent(in) :: hgrid,elecfield
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), intent(out) :: eion
     real(kind=8), dimension(*), intent(out) :: pot_ion
   end subroutine createIonicPotential

   subroutine import_gaussians(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
        nat,norb,norbp,occup,n1,n2,n3,nvctrp,hgrid,rxyz, & 
        rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
        atomnames,ntypes,iatype,pkernel,psppar,npspcode,ixc,&
        psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     logical, intent(in) :: parallel
     character(len=20), dimension(100), intent(in) :: atomnames
     character(len=1), intent(in) :: datacode
     integer, intent(in) :: iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
     integer, intent(in) :: nspin
     real(kind=8), dimension(norb), intent(in) :: spinar
     real(kind=8), intent(in) :: hgrid
     real(kind=8), intent(out) :: accurex
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: npspcode
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
   end subroutine import_gaussians

   subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        nat,natsc,norb,norbp,n1,n2,n3,nvctrp,hgrid,rxyz, & 
        rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
        atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,ixc,&
        ppsi,ppsit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     logical, intent(in) :: parallel
     character(len=20), dimension(100), intent(in) :: atomnames
     character(len=1), intent(in) :: datacode
     integer, intent(in) :: iproc,nproc,nat,natsc,ntypes,norb,norbp,n1,n2,n3,ixc
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
     integer, intent(in) :: nspin
     real(kind=8), dimension(norb), intent(in) :: spinar
     real(kind=8), intent(in) :: hgrid
     real(kind=8), intent(out) :: accurex
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: iasctype,npspcode,nzatom,nelpsp
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(*), intent(in) :: pkernel
     real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
     real(kind=8), dimension(norb), intent(out) :: eval
     real(kind=8), dimension(:,:), pointer :: ppsi,ppsit
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
   
   subroutine sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
        wfd,psi,rho,nrho,nscatterarr,nspin,spinar,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(convolutions_bounds), intent(in) :: bounds
     logical, intent(in) ::  parallel
     integer, intent(in) :: iproc,nproc,norb,norbp,nrho,nspin
     integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hgrid
     integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
     real(kind=8), dimension(norb), intent(in) :: occup,spinar
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(max(nrho,1),nspin), intent(out) :: rho
   end subroutine sumrho

   subroutine HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
        psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
        wfd,bounds,nlpspd,proj,ngatherarr,n3p,&
        potential,psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,spinar)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     type(convolutions_bounds), intent(in) :: bounds
     logical, intent(in) :: parallel
     character(len=1), intent(in) :: datacode
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,n3p
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin
     real(kind=8), intent(in) :: hgrid
     integer, dimension(ntypes), intent(in) :: npspcode
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
     real(kind=8), dimension(norb), intent(in) :: occup,spinar
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(*), intent(in) :: potential
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: hpsi
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
   end subroutine HamiltonianApplication

   subroutine hpsitopsi(iter,parallel,iproc,nproc,norb,norbp,occup,hgrid,n1,n2,n3,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,wfd,kbounds,&
        eval,ncong,mids,idsx,ads,energy,energy_old,alpha,gnrm,scprsum,&
        psi,psit,hpsi,psidst,hpsidst,nspin,spinar)
     use module_types               
     implicit none
     logical, intent(in) :: parallel
     integer, intent(in) :: iter,iproc,nproc,n1,n2,n3,norb,norbp,ncong,mids,idsx
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp,nspin
     real(kind=8), intent(in) :: hgrid,energy,energy_old
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinar
     real(kind=8), intent(inout) :: alpha
     real(kind=8), intent(inout) :: gnrm,scprsum
     real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi
     real(kind=8), dimension(:,:,:), pointer :: psidst,hpsidst,ads
     type(kinetic_bounds), intent(in) :: kbounds
     type(wavefunctions_descriptors), intent(in) :: wfd
   end subroutine hpsitopsi

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

   subroutine local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
        n1,n2,n3,n3pi,i3s,rho,pot,floc)
     implicit none
     integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3,n3pi,i3s
     real(kind=8), intent(in) :: hgrid
     character(len=20), dimension(100), intent(in) :: atomnames
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(*), intent(in) :: rho,pot
     integer, dimension(nat), intent(in) :: iatype
     integer, dimension(ntypes), intent(in) :: nelpsp
     real(kind=8), dimension(3,nat), intent(out) :: floc
   end subroutine local_forces

   subroutine projectors_derivatives(iproc,n1,n2,n3,ntypes,nat,norb,iatype,psppar,nlpspd,proj,  &
        rxyz,radii_cf,cpmult,fpmult,hgrid,derproj)
     use module_types
     implicit none
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     integer, intent(in) :: iproc,ntypes,nat,norb
     integer, intent(in) :: n1,n2,n3
     real(kind=8),intent(in) :: cpmult,fpmult,hgrid 
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(nlpspd%nprojel,3), intent(out) :: derproj
   end subroutine projectors_derivatives

   subroutine nonlocal_forces(iproc,ntypes,nat,norb,norbp,iatype,psppar,npspcode,occup,&
        nlpspd,proj,derproj,wfd,psi,fsep)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(in) :: nlpspd
     integer, intent(in) :: iproc,ntypes,nat,norb,norbp
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     integer, dimension(ntypes), intent(in) :: npspcode
     real(kind=8), dimension(norb), intent(in) :: occup
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(nlpspd%nprojel,3), intent(in) :: derproj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
     real(kind=8), dimension(3,nat), intent(inout) :: fsep
   end subroutine nonlocal_forces

   subroutine CalculateTailCorrection(iproc,nproc,n1,n2,n3,rbuf,norb,norbp,nat,ntypes,&
        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,nlpspd,ncongt,psppar,npspcode,eval,&
        pot,hgrid,rxyz,radii_cf,crmult,frmult,iatype,atomnames,nspin,spinar,&
        proj,psi,occup,output_grid,parallel,ekin_sum,epot_sum,eproj_sum)
     use module_types
     implicit none
     type(wavefunctions_descriptors), intent(in) :: wfd
     type(nonlocal_psp_descriptors), intent(inout) :: nlpspd
     logical, intent(in) :: output_grid,parallel
     character(len=20), dimension(100), intent(in) :: atomnames
     integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nat,ntypes,ncongt,nspin
     integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
     real(kind=8), intent(in) :: hgrid,crmult,frmult,rbuf
     real(kind=8), intent(out) :: ekin_sum,epot_sum,eproj_sum
     integer, dimension(ntypes), intent(in) :: npspcode
     integer, dimension(nat), intent(in) :: iatype
     real(kind=8), dimension(norb), intent(in) :: occup,eval,spinar
     real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
     real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
     real(kind=8), dimension(3,nat), intent(in) :: rxyz
     real(kind=8), dimension(2*n1+31,2*n2+31,2*n3+31,nspin), intent(in) :: pot
     real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
     real(kind=8), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
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

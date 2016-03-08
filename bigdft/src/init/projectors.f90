!> @file
!!  Routines to handle projectors
!! @author
!!    Copyright (C) 2010-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Localize the projectors for pseudopotential calculations
subroutine localize_projectors(n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,&
     logrid,at,orbs,nl)
  use module_base
  use module_types, only: atoms_data,orbitals_data
  use gaussians, only: gaussian_basis_iter, gaussian_iter_start, gaussian_iter_next_shell
  use yaml_output
  use psp_projectors_base, only: DFT_PSP_projectors
  use psp_projectors, only: bounds_to_plr_limits, pregion_size
  use public_enums, only: PSPCODE_PAW
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(DFT_PSP_projectors), intent(inout) :: nl
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
  logical, dimension(0:n1,0:n2,0:n3), intent(inout) :: logrid
  !Local variables
  !n(c) logical :: cmplxprojs
  type(gaussian_basis_iter) :: iter
  integer :: istart,ityp,iat,mproj,nl1,nu1,nl2,nu2,nl3,nu3,mvctr,mseg,nprojelat,i,l
  integer :: ikpt,nkptsproj,ikptp,izero
  real(gp) :: maxfullvol,totfullvol,totzerovol,fullvol,maxrad,maxzerovol,rad

  call f_routine(id='localize_projectors')
  
  istart=1
  nl%nproj=0
  nl%nprojel=0


  do iat=1,at%astruct%nat

     mproj = 0
     call gaussian_iter_start(nl%proj_G, iat, iter)
     do
        if (.not. gaussian_iter_next_shell(nl%proj_G, iter)) exit
        mproj = mproj + 2 * iter%l - 1
     end do
     
     !assign the number of projector to the localization region
     nl%pspd(iat)%mproj=mproj

     if (mproj /= 0) then 

        !if some projectors are there at least one locreg interacts with the psp
        nl%pspd(iat)%nlr=1

        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))')&
        !     'projector descriptors for atom with mproj ',iat,mproj
        nl%nproj=nl%nproj+mproj

        ! coarse grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),3),&
             cpmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        !these routines are associated to the handling of a locreg
        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,&
             at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),at%radii_cf(1,3),&
             cpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)

        nl%pspd(iat)%plr%wfd%nseg_c=mseg
        nl%pspd(iat)%plr%wfd%nvctr_c=mvctr

        istart=istart+mvctr*mproj

        nprojelat=mvctr*mproj

        !print *,'iat,mvctr',iat,mvctr,mseg,mproj

        ! fine grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,2,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),&
             at%radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr

        nl%pspd(iat)%plr%wfd%nseg_f=mseg
        nl%pspd(iat)%plr%wfd%nvctr_f=mvctr

        istart=istart+7*mvctr*mproj
        nprojelat=nprojelat+7*mvctr*mproj
        nl%nprojel=max(nl%nprojel,nprojelat)

        !print *,'iat,nprojelat',iat,nprojelat,mvctr,mseg

     else  !(atom has no nonlocal PSP, e.g. H)
        !Pb of inout
        izero=0

        !this is not really needed, see below
        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             izero,izero,izero,izero,izero,izero)

        nl%pspd(iat)%plr%wfd%nseg_c=0
        nl%pspd(iat)%plr%wfd%nvctr_c=0
        nl%pspd(iat)%plr%wfd%nseg_f=0
        nl%pspd(iat)%plr%wfd%nvctr_f=0

        !! the following is necessary to the creation of preconditioning projectors
        !! coarse grid quantities ( when used preconditioners are applied to all atoms
        !! even H if present )
        call pregion_size(at%astruct%geocode,rxyz(1,iat),&
             at%radii_cf(at%astruct%iatype(iat),3),cpmult, &
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        ! fine grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,2,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)
     endif
  enddo

  !control the strategy to be applied following the memory limit
  !if the projectors takes too much memory allocate only one atom at the same time
  !control the memory of the projectors expressed in GB
  if (memorylimit /= 0.e0 .and. .not. DistProjApply .and. &
       real(istart-1,kind=4) > memorylimit*134217728.0e0) then
!!$     if (iproc == 0) then
!!$        write(*,'(44x,a)') '------ On-the-fly projectors application'
!!$     end if
     DistProjApply =.true.
  end if

  !assign the distprojapply value to the structure
  nl%on_the_fly=DistProjApply

! Skip this zerovol test for PAW, since the HGH/GTH parameters 
! are not used.
  if (all(at%npspcode(:) /= PSPCODE_PAW)) then
     !calculate the fraction of the projector array used for allocate zero values
     !control the hardest and the softest gaussian
     totzerovol=0.0_gp
     maxfullvol=0.0_gp
     totfullvol=0.0_gp
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        maxrad=min(maxval(at%psppar(1:4,0,ityp)),cpmult/15.0_gp*at%radii_cf(ityp,3))
        nl%zerovol=0.0_gp
        fullvol=0.0_gp
        do l=1,4
           do i=1,3
              if (at%psppar(l,i,ityp) /= 0.0_gp) then
                 rad=min(at%psppar(l,0,ityp),cpmult/15.0_gp*at%radii_cf(ityp,3))
                 nl%zerovol=nl%zerovol+(maxrad**3-rad**3)
                 fullvol=fullvol+maxrad**3
              end if
           end do
        end do
        if (fullvol >= maxfullvol .and. fullvol > 0.0_gp) then
           maxzerovol=nl%zerovol/fullvol
           maxfullvol=fullvol
        end if
        totzerovol=totzerovol+nl%zerovol
        totfullvol=totfullvol+fullvol
     end do

     !assign the total quantity per atom
     nl%zerovol=0.d0
     if (totfullvol /= 0.0_gp) then
        if (nl%on_the_fly) then
           nl%zerovol=maxzerovol
        else
           nl%zerovol=totzerovol/totfullvol
        end if
     end if
  end if !npspcode == PSPCODE_PAW

  !here is the point in which the projector strategy should be decided
  !DistProjApply shoud never change after this point

  !number of elements of the projectors
  if (.not. nl%on_the_fly) nl%nprojel=istart-1

  !Compute the multiplying coefficient for nprojel in case of imaginary k points.
  !activate the complex projector if there are kpoints
  !TO BE COMMENTED OUT
  !n(c) cmplxprojs= (orbs%kpts(1,1)**2+orbs%kpts(2,1)**2+orbs%kpts(3,1)**2 >0) .or. orbs%nkpts>1
  nkptsproj=1
  if ((.not. nl%on_the_fly) .and. orbs%norbp > 0) then
     nkptsproj = 0
     do ikptp=1,orbs%nkptsp
        ikpt=orbs%iskpts+ikptp
        if (any(orbs%kpts(:,ikpt) /= 0.0_gp) .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = nkptsproj + 2
        else
           nkptsproj = nkptsproj + 1
        end if
     end do
  else if (nl%on_the_fly) then
     do ikptp=1,orbs%nkptsp
        ikpt=orbs%iskpts+ikptp
        if (any(orbs%kpts(:,ikpt) /= 0.0_gp) .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = max(nkptsproj, 2)
        end if
     end do
  end if
  nl%nprojel=nkptsproj*nl%nprojel

  !print *,'iproc,nkptsproj',iproc,nkptsproj,nlpspd%nprojel,orbs%iskpts,orbs%iskpts+orbs%nkptsp

  call f_release_routine()

END SUBROUTINE localize_projectors

!> Localize the projectors for pseudopotential calculations
subroutine localize_projectors_new(n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,&
     logrid,at,nl)
  use module_base
  use module_types, only: atoms_data,orbitals_data
  use gaussians, only: gaussian_basis_iter, gaussian_iter_start, gaussian_iter_next_shell
  use psp_projectors_base, only: DFT_PSP_projectors
  use psp_projectors, only: bounds_to_plr_limits, pregion_size
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(DFT_PSP_projectors), intent(inout) :: nl
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
  logical, dimension(0:n1,0:n2,0:n3), intent(inout) :: logrid
  !Local variables
  !n(c) logical :: cmplxprojs
  type(gaussian_basis_iter) :: iter
  integer :: istart,ityp,iat,mproj,nl1,nu1,nl2,nu2,nl3,nu3,mvctr,mseg,i,l
  integer :: izero

  call f_routine(id='localize_projectors')

  nl%nproj=0

  do iat=1,at%astruct%nat

     mproj = 0
     call gaussian_iter_start(nl%proj_G, iat, iter)
     do
        if (.not. gaussian_iter_next_shell(nl%proj_G, iter)) exit
        mproj = mproj + 2 * iter%l - 1
     end do

     !assign the number of projector to the localization region
     nl%pspd(iat)%mproj=mproj

     if (mproj /= 0) then 

        !if some projectors are there at least one locreg interacts with the psp
        nl%pspd(iat)%nlr=1

        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))')&
        !     'projector descriptors for atom with mproj ',iat,mproj
        nl%nproj=nl%nproj+mproj

        ! coarse grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),3),&
             cpmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        !these routines are associated to the handling of a locreg
        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,&
             at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),at%radii_cf(1,3),&
             cpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)

        nl%pspd(iat)%plr%wfd%nseg_c=mseg
        nl%pspd(iat)%plr%wfd%nvctr_c=mvctr

        !print *,'iat,mvctr',iat,mvctr,mseg,mproj

        ! fine grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,2,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),&
             at%radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr

        nl%pspd(iat)%plr%wfd%nseg_f=mseg
        nl%pspd(iat)%plr%wfd%nvctr_f=mvctr

     else  !(atom has no nonlocal PSP, e.g. H)
        !Pb of inout
        izero=0

        !this is not really needed, see below
        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             izero,izero,izero,izero,izero,izero)

        nl%pspd(iat)%plr%wfd%nseg_c=0
        nl%pspd(iat)%plr%wfd%nvctr_c=0
        nl%pspd(iat)%plr%wfd%nseg_f=0
        nl%pspd(iat)%plr%wfd%nvctr_f=0

        !! the following is necessary to the creation of preconditioning projectors
        !! coarse grid quantities ( when used preconditioners are applied to all atoms
        !! even H if present )
        call pregion_size(at%astruct%geocode,rxyz(1,iat),&
             at%radii_cf(at%astruct%iatype(iat),3),cpmult, &
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,1,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)

        ! fine grid quantities
        call pregion_size(at%astruct%geocode,rxyz(1,iat),at%radii_cf(at%astruct%iatype(iat),2),fpmult,&
             hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

        call bounds_to_plr_limits(.true.,2,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)
     endif
  enddo

  call f_release_routine()

END SUBROUTINE localize_projectors_new


!> Fill the proj array with the PSP projectors or their derivatives, following idir value
subroutine fill_projectors(lr,hgrids,astruct,ob,rxyz,nlpsp,idir)
  use module_base
  use locregs
  use psp_projectors_base
  use psp_projectors
  use yaml_output
  use orbitalbasis
  use module_atoms
  implicit none
  integer, intent(in) :: idir
  real(gp), dimension(3), intent(in) :: hgrids
  type(atomic_structure), intent(in) :: astruct
  type(orbital_basis), intent(in) :: ob
  !type(orbitals_data), intent(in) :: orbs
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(locreg_descriptors),intent(in) :: lr !< Global localisation region
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
  !Local variables
  logical :: overlap
  integer :: istart_c,iat,iproj,nwarnings,ikpt,iskpt,iekpt
  type(ket) :: psi_it
  type(atoms_iterator) :: atit

  call f_routine(id='fill_projectors')

  !if (iproc.eq.0 .and. nlpspd%nproj /=0 .and. idir==0) &
  !write(*,'(1x,a)',advance='no') 'Calculating wavelets expansion of projectors...'
  !warnings related to the projectors norm
  nwarnings=0
  !allocate these vectors up to the maximum size we can get
  istart_c=1

  !initialize the orbital basis object, for psi and hpsi
  !call orbital_basis_associate(ob,orbs=orbs)
  !iterate over the orbital_basis
  psi_it=orbital_basis_iterator(ob)
  loop_kpt: do while(ket_next_kpt(psi_it))
     loop_lr: do while(ket_next_locreg(psi_it,ikpt=psi_it%ikpt))
        iproj=0
        atit = atoms_iter(astruct)
        loop_atoms: do while(atoms_iter_next(atit))
           overlap = projector_has_overlap(atit%iat,psi_it%ilr,psi_it%lr, lr, nlpsp)
           if(.not. overlap) cycle loop_atoms
           !this routine is defined to uniformise the call for on-the-fly application
           call atom_projector(nlpsp,atit%ityp, atit%iat, &
                atit%name, astruct%geocode, idir, lr,&
                hgrids(1),hgrids(2),hgrids(3), &
                psi_it%kpoint(1),psi_it%kpoint(2),psi_it%kpoint(3),&
                istart_c, iproj, nwarnings)
        end do loop_atoms
     end do loop_lr
     !this part should be updated by update_nlpspd according to the overlap
     if (iproj /= nlpsp%nproj) then
        call f_err_throw('Incorrect number of projectors created',&
             err_name='BIGDFT_RUNTIME_ERROR')
     end if
     ! projector part finished
  end do loop_kpt

!call orbital_basis_release(ob)

if (istart_c-1 /= nlpsp%nprojel) then
  call yaml_warning('Incorrect once-and-for-all psp generation')
  stop
end if

if (nwarnings /= 0 .and. bigdft_mpi%iproc == 0 .and. nlpsp%nproj /=0 .and. idir == 0) then
  call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
  if (nwarnings /= 0) then
     call yaml_newline()
     call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
  end if
end if

call f_release_routine()

END SUBROUTINE fill_projectors!_new


!> Fill the proj array with the PSP projectors or their derivatives, following idir value
subroutine fill_projectors_old(lr,hx,hy,hz,at,orbs,rxyz,nlpsp,idir)
use module_base
use module_types
use yaml_output
implicit none
integer, intent(in) :: idir
real(gp), intent(in) :: hx,hy,hz
type(atoms_data), intent(in) :: at
type(orbitals_data), intent(in) :: orbs
type(DFT_PSP_projectors), intent(inout) :: nlpsp
type(locreg_descriptors),intent(in) :: lr
real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!Local variables
!n(c) integer, parameter :: nterm_max=20 !if GTH nterm_max=4
integer :: istart_c,iat,iproj,nwarnings,ikpt,iskpt,iekpt

call f_routine(id='fill_projectors')

!if (iproc.eq.0 .and. nlpspd%nproj /=0 .and. idir==0) &
!write(*,'(1x,a)',advance='no') 'Calculating wavelets expansion of projectors...'
!warnings related to the projectors norm
nwarnings=0
!allocate these vectors up to the maximum size we can get
istart_c=1

!create projectors for any of the k-point hosted by the processor
!starting kpoint
if (orbs%norbp > 0) then
  iskpt=orbs%iokpt(1)
  iekpt=orbs%iokpt(orbs%norbp)
else
  iskpt=1
  iekpt=1
end if

do ikpt=iskpt,iekpt
  iproj=0
  do iat=1,at%astruct%nat
     !this routine is defined to uniformise the call for on-the-fly application
     call atom_projector(nlpsp, at%astruct%iatype(iat), iat, &
          & at%astruct%atomnames(at%astruct%iatype(iat)), &
          & at%astruct%geocode, idir, lr, hx, hy, hz, &
          & orbs%kpts(1,ikpt), orbs%kpts(2,ikpt), orbs%kpts(3,ikpt), &
          & istart_c, iproj, nwarnings)
  enddo
  if (iproj /= nlpsp%nproj) then
     call yaml_warning('Incorrect number of projectors created')
  end if
  ! projector part finished
end do

if (istart_c-1 /= nlpsp%nprojel) then
  call yaml_warning('Incorrect once-and-for-all psp generation')
  stop
end if

if (nwarnings /= 0 .and. bigdft_mpi%iproc == 0 .and. nlpsp%nproj /=0 .and. idir == 0) then
  call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
  if (nwarnings /= 0) then
     call yaml_newline()
     call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
     !write(*,'(1x,a,i0,a)') 'found ',nwarnings,' warnings.'
     !write(*,'(1x,a)') 'Some projectors may be too rough.'
     !write(*,'(1x,a,f6.3)') 'Consider the possibility of modifying hgrid and/or the localisation radii.'
  end if
end if

call f_release_routine()

END SUBROUTINE fill_projectors_old


!> Create projectors from gaussian decomposition.
subroutine atom_projector(nl, ityp, iat, atomname, &
  & geocode, idir, lr, hx, hy, hz, kx, ky, kz, &
  & istart_c, iproj, nwarnings)
use module_base
use locregs
use gaussians, only: gaussian_basis_new, gaussian_basis_iter, &
    & gaussian_iter_start, gaussian_iter_next_shell, gaussian_iter_next_gaussian
use psp_projectors_base, only: DFT_PSP_projectors
use yaml_output, only: yaml_warning
use yaml_strings, only: yaml_toa
implicit none
type(DFT_PSP_projectors), intent(inout) :: nl
integer, intent(in) :: ityp, iat
character(len = 1), intent(in) :: geocode
  type(locreg_descriptors), intent(in) :: lr
  integer, intent(in) :: idir
  character(len = *), intent(in) :: atomname
  real(gp), intent(in) :: hx, hy, hz, kx, ky, kz
  integer, intent(inout) :: istart_c, iproj, nwarnings

  type(gaussian_basis_iter) :: iter, iterM
  integer :: nc, ncplx_k, lmax, np, mbvctr_c, mbvctr_f, mbseg_c, mbseg_f
  real(gp) :: scpr
  real(gp), dimension(nl%proj_G%ncplx) :: coeff, expo
  logical :: use_tmp
  real(wp),allocatable::proj_tmp(:)

  call f_routine(id='atom_projector')

  call plr_segs_and_vctrs(nl%pspd(iat)%plr,mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)

  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx_k=1
  else
     ncplx_k=2
  end if

  ! Start a gaussian iterator.
  call gaussian_iter_start(nl%proj_G, iat, iter)

  ! Maximum number of terms for every projector.
  lmax = 0
  use_tmp = .false.
  iterM = iter
  do
     if (.not. gaussian_iter_next_shell(nl%proj_G, iterM)) exit
     lmax = max(lmax, iterM%l)
     if (iterM%ndoc > 1) use_tmp = .true.
  end do

  if (use_tmp) then
     nc = (mbvctr_c+7*mbvctr_f)*(2*lmax-1)*ncplx_k
     proj_tmp = f_malloc(nc, id = 'proj_tmp')
  end if

  ! Loop on shell.
  do
     if (.not. gaussian_iter_next_shell(nl%proj_G, iter)) exit
     nc = (mbvctr_c+7*mbvctr_f) * (2*iter%l-1) * ncplx_k
     if (istart_c + nc > nl%nprojel+1) stop 'istart_c > nprojel+1'
     ! Loop on contraction, treat the first gaussian separately for performance reasons.
     if (gaussian_iter_next_gaussian(nl%proj_G, iter, coeff, expo)) &
          & call projector(geocode, iat, idir, iter%l, iter%n, coeff, expo, &
          & nl%pspd(iat)%gau_cut, nl%proj_G%rxyz(:, iat), lr%ns1, lr%ns2, lr%ns3, lr%d%n1, lr%d%n2, lr%d%n3, &
          & hx, hy, hz, kx, ky, kz, ncplx_k, nl%proj_G%ncplx, &
          & mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, nl%pspd(iat)%plr%wfd%keyvglob, nl%pspd(iat)%plr%wfd%keyglob, &
          & nl%wpr,nl%proj(istart_c:))
     do
        if (.not. gaussian_iter_next_gaussian(nl%proj_G, iter, coeff, expo)) exit
        call projector(geocode, iat, idir, iter%l, iter%n, coeff, expo, &
             & nl%pspd(iat)%gau_cut, nl%proj_G%rxyz(:, iat), lr%ns1, lr%ns2, lr%ns3, lr%d%n1, lr%d%n2, lr%d%n3, &
             & hx, hy, hz, kx, ky, kz, ncplx_k, nl%proj_G%ncplx, &
             & mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, nl%pspd(iat)%plr%wfd%keyvglob, nl%pspd(iat)%plr%wfd%keyglob, &
             & nl%wpr, proj_tmp)
        call axpy(nc, 1._wp, proj_tmp(1), 1, nl%proj(istart_c), 1)
     end do
     ! Check norm for each proj.
     if (idir == 0 .and. nl%normalized) then
        do np = 1, 2 * iter%l - 1
           !here the norm should be done with the complex components
           call wnrm_wrap(ncplx_k,mbvctr_c,mbvctr_f, &
                & nl%proj(istart_c + (np - 1) * (mbvctr_c+7*mbvctr_f) * ncplx_k),scpr)
           !print '(a,3(i6),1pe14.7,2(i6))','iat,l,m,scpr',iat,l,m,scpr,idir,istart_c
           if (abs(1.d0-scpr) > 1.d-2) then
              if (abs(1.d0-scpr) > 1.d-1) then
                 if (bigdft_mpi%iproc == 0) call yaml_warning( &
                      'Norm of the nonlocal PSP [atom ' // trim(yaml_toa(iat)) // &
                      ' (' // trim(atomname) // ') l=' // trim(yaml_toa(iter%l)) // &
                      ' m=' // trim(yaml_toa(iter%n)) // ' is ' // trim(yaml_toa(scpr)) // &
                      ' while it is supposed to be about 1.0.')
                 !stop commented for the moment
                 !restore the norm of the projector
                 !call wscal_wrap(mbvctr_c,mbvctr_f,1.0_gp/sqrt(scpr),proj(istart_c))
              else
                 nwarnings=nwarnings+1
              end if
           end if
        end do
     end if
     iproj    = iproj + 2*iter%l-1
     istart_c = istart_c + nc
  end do

  if (use_tmp) call f_free(proj_tmp)

  call f_release_routine()

end subroutine atom_projector


subroutine projector(geocode,iat,idir,l,i,factor,gau_a,rpaw,rxyz,&
     ns1,ns2,ns3,n1,n2,n3,hx,hy,hz,kx,ky,kz,ncplx_k,ncplx_g,&
     mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,wpr,proj)
  use module_base
  use module_types
  use psp_projectors_base, only: workarrays_projectors
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: ns1,ns2,ns3,n1,n2,n3
  integer, intent(in) :: iat,idir,l,i,mbvctr_c,mbvctr_f,mseg_c,mseg_f,ncplx_k,ncplx_g
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz,rpaw
  real(gp),dimension(ncplx_g),intent(in)::gau_a,factor
  !integer, dimension(2,3), intent(in) :: nboxp_c,nboxp_f
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p

  real(gp), dimension(3), intent(in) :: rxyz
  type(workarrays_projectors),intent(inout) :: wpr
  real(wp), dimension((mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx_k), intent(out) :: proj
  !Local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: m,iterm
  !integer :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  integer :: istart_c,nterm,idir2
  real(gp) :: fpi,fgamma,rx,ry,rz
  integer, dimension(3) :: nterm_arr
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(3,nterm_max,3) :: lxyz_arr
  real(gp), dimension(ncplx_g,nterm_max) :: factors
  real(gp), dimension(nterm_max,3) :: fac_arr
  real(gp), dimension(ncplx_g) :: gau_c


  !call f_routine(id='projector')

  !this value can also be inserted as a parameter
  if (ncplx_g == 1) then
     !fpi=pi^-1/4 pi^-1/2, pi^-1/4 comes from sqrt(gamma(x)) and pi^-1/2 from Ylm.
     !fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)
     fpi=0.42377720812375763_gp
     ! gau_a is real and given as alpha, need to convert it back as coefficient.
     gau_c(1) = 1._gp / sqrt(2._gp * gau_a(1))
     fgamma=sqrt(2.0_gp)*fpi/(sqrt(gau_c(1))**(2*(l-1)+4*i-1))
  else
     fpi=0.56418958354775628_gp
     if(l==1) then
        fgamma= 0.70710678118654757_gp !1.0/sqrt(2.0)
        !1/sqrt(2.0)*fac_arr(1)=1/sqrt(2.0)* [1/sqrt(2.0)]=1/2
        !1/2*fpi=1/2* [1/sqrt(pi)]=1/2 1/sqrt(pi) Factor for Y_00
     elseif(l==2) then
        fgamma= 0.8660254037844386_gp !sqrt(3)/2.0
     elseif(l==3) then
        fgamma= 1.3693063937629153_gp  !sqrt(3*5)/(2.0*sqrt(2))
     elseif(l==4) then
        fgamma= 2.5617376914898995_gp  !sqrt(7*5*3)/(4.0) 
     else
        write(*,'(1x,a)')'error found!'
        write(*,'(1x,a,i4)')&
             'gamma_factor: l should be between 1 and 3, but l= ',l
        stop
     end if
     fgamma = fgamma * fpi
     gau_c(1) = 1._gp / sqrt(2._gp * gau_a(1))
     gau_c(2) = gau_a(2)
  end if

  rx=rxyz(1) 
  ry=rxyz(2) 
  rz=rxyz(3)

  istart_c=1
  !start of the projectors expansion routine
  do m=1,2*l-1

     if (idir==0) then !normal projector calculation case
        idir2=1
        call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)

        do iterm=1,nterm
           factors(:,iterm)=factor(:)*fgamma*fac_arr(iterm,idir2)
        end do
     else !calculation of projector derivative
        idir2=mod(idir-1,3)+1
        call calc_coeff_derproj(l,i,m,nterm_max,gau_c(1),nterm_arr,lxyz_arr,fac_arr)
        nterm=nterm_arr(idir2)

        do iterm=1,nterm
           factors(:,iterm)=factor(:)*fgamma*fac_arr(iterm,idir2)
           lx(iterm)=lxyz_arr(1,iterm,idir2)
           ly(iterm)=lxyz_arr(2,iterm,idir2)
           lz(iterm)=lxyz_arr(3,iterm,idir2)        

!       nterm=nterm_arr(idir)
!       do iterm=1,nterm
!          factors(iterm)=fgamma*fac_arr(iterm,idir)
!          lx(iterm)=lxyz_arr(1,iterm,idir)
!          ly(iterm)=lxyz_arr(2,iterm,idir)
!          lz(iterm)=lxyz_arr(3,iterm,idir)

! sequence: 11 21 31 12 22 32 13 23 33 

!if (idir > 3) then
!        if (idir < 7) then
!lx(iterm)=lx(iterm)+1
!        else if (idir < 10) then
!ly(iterm)=ly(iterm)+1
!        else 
!lz(iterm)=lz(iterm)+1
!        endif
!endif
!seq : 11 22 33 12 23 13
select case(idir)
case(4,9)
   lx(iterm)=lx(iterm)+1
case(5,7)
   ly(iterm)=ly(iterm)+1
case(6,8)
   lz(iterm)=lz(iterm)+1
end select
!if (idir == 4 .or. idir == 9) lx(iterm)=lx(iterm)+1
!if (idir == 5 .or. idir == 7) ly(iterm)=ly(iterm)+1
!if (idir == 6 .or. idir == 8) lz(iterm)=lz(iterm)+1

        end do
     end if

!!$     write(*,*) geocode,nterm,ns1,ns2,ns3,n1,n2,n3,lx(1:nterm),ly(1:nterm),lz(1:nterm), rpaw
!!$     write(*,*) hx,hy,hz,kx,ky,kz,ncplx_g,ncplx_k
     call crtproj(geocode,nterm,ns1,ns2,ns3,n1,n2,n3,&
          hx,hy,hz,kx,ky,kz,ncplx_g,ncplx_k,&
          gau_c,factors,rx,ry,rz,lx,ly,lz,&
          mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj(istart_c),wpr,rpaw)
 
     !do iterm=1,nterm
     !   if (iproc.eq.0) write(*,'(1x,a,i0,1x,a,1pe10.3,3(1x,i0))') &
     !        'projector: iat,atomname,gau_a,lx,ly,lz ', & 
     !        iat,trim(at%astruct%atomnames(at%astruct%iatype(iat))),gau_a,lx(iterm),ly(iterm),lz(iterm)
     !enddo
     !end testing

  !call f_release_routine()

     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx_k
  enddo

  !call f_release_routine()

END SUBROUTINE projector

!> Returns the compressed form of a Gaussian projector 
!! @f$ x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-rx)^2 + (y-ry)^2 + (z-rz)^2 )) @f$
!! in the array proj.
subroutine crtproj(geocode,nterm,ns1,ns2,ns3,n1,n2,n3, & 
     hx,hy,hz,kx,ky,kz,ncplx_g,ncplx_k,&
     gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
     mvctr_c,mvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj,wpr,gau_cut)
  use module_base
  use module_types
  use psp_projectors_base, only: workarrays_projectors, NCPLX_MAX
  !use gaussians
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: nterm,mvctr_c,mvctr_f,mseg_c,mseg_f
  integer, intent(in) :: ncplx_g,ncplx_k,ns1,ns2,ns3,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz,rx,ry,rz,kx,ky,kz
  real(gp), intent(in) :: gau_cut
  integer, dimension(nterm), intent(in) :: lx,ly,lz
  real(gp), dimension(ncplx_g,nterm), intent(in) :: fac_arr
  real(gp), dimension(ncplx_g),intent(in):: gau_a
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p
  type(workarrays_projectors),intent(inout) :: wpr
  real(wp), dimension((mvctr_c+7*mvctr_f)*ncplx_k), intent(out) :: proj
  !Local variables
  character(len=*), parameter :: subname='crtproj'
  logical :: perx,pery,perz !variables controlling the periodicity in x,y,z
  integer :: iterm,n_gau,ml1,ml2,ml3,mu1,mu2,mu3,i1,i2,i3
  integer :: ncplx_w,n1p1,np,i0jj
  integer :: j1,i0,j0,jj,ii,i,iseg,ind_f,ind_c
  integer :: mvctr1, mvctr2, mvctr_cf, mvctr_cf2, iskip
  !integer :: counter !test
  !real(wp) :: re_cmplx_prod,im_cmplx_prod
  real(gp), dimension(ncplx_g) :: factor, one
  !real(gp) :: err_norm
  !real(wp), allocatable, dimension(:,:,:) :: work
  real(wp) :: wprojyz, wprojyz11, wprojyz12, wprojyz21, wprojyz22
  !Variables for OpenMP
  !!$ integer :: ithread,nthread,ichunk
  !!$ integer :: omp_get_thread_num,omp_get_num_threads

!!  integer :: ncount0,ncount_rate,ncount_max,ncount1,ncount2


  !call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  call f_routine(id='crtproj')

  ! rename region boundaries
  n1p1=n1+1
  np=n1p1*(n2+1)
  mvctr_cf=mvctr_c+7*mvctr_f
  mvctr_cf2=2*mvctr_c+7*mvctr_f

  !wproj is complex for PAW and kpoints.
  ncplx_w=max(ncplx_g,ncplx_k,1)

  ! The workarrays wpr%wprojx, wpr%wprojy, wpr%wprojz are allocated with the
  ! first dimension equala to NCPLX_MAX (which is 2). However the routine gauss_to_daub_k always
  ! assumes the correct value for ncplx_w and thus fills the arrays
  ! contiguously. Therefore in the non-complex case one has to fill the holes in
  ! thw workarrays.
  if (ncplx_w==NCPLX_MAX) then
      iskip = 1
  else
      iskip = 2
  end if

  ! Check the dimensions
  if (size(wpr%wprojx,2)<n1) call f_err_throw('workarray wpr%wprojx too small',err_name='BIGDFT_RUNTIME_ERROR')
  if (size(wpr%wprojy,2)<n2) call f_err_throw('workarray wpr%wprojy too small',err_name='BIGDFT_RUNTIME_ERROR')
  if (size(wpr%wprojz,2)<n3) call f_err_throw('workarray wpr%wprojz too small',err_name='BIGDFT_RUNTIME_ERROR')

  !if(ncplx_wproj==2 .or. nterm>1) proj=0.d0 !initialize to zero in this cases

!  allocate(work(0:nw,2,2+ndebug),stat=i_stat)  !always use complex value
!  call memocc(i_stat,work,'work',subname)

  !check that the number of elements of the projector is coherent
  mvctr1=0
  do iseg=1,mseg_c
     mvctr1=mvctr1+keyg_p(2,iseg)-keyg_p(1,iseg)+1
  end do
  mvctr2=0
  do iseg=mseg_c+1,mseg_c+mseg_f
     mvctr2=mvctr2+keyg_p(2,iseg)-keyg_p(1,iseg)+1
  end do
  
  if (mvctr1 /=  mvctr_c) then
     write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 1): mvctr /= mvctr_c ',mvctr1,mvctr_c
     stop
  end if

  if (mvctr2 /= mvctr_f) then
     write(*,'(1x,a,i0,1x,i0)')' ERROR (crtproj 1): mvctr /= mvctr_f ',mvctr2,mvctr_f
     stop
  end if

  !REALLY SLOW ON VESTA, TEMPORARY CHANGE ONLY
  !!wprojx = f_malloc((/ 1.to.ncplx_w, 0.to.n1, 1.to.2, 1.to.nterm /),id='wprojx')
  !!wprojy = f_malloc((/ 1.to.ncplx_w, 0.to.n2, 1.to.2, 1.to.nterm /),id='wprojy')
  !!wprojz = f_malloc((/ 1.to.ncplx_w, 0.to.n3, 1.to.2, 1.to.nterm /),id='wprojz')

  !!if (size(wpr%wprojx)<size(wprojx)) stop 'size x'
  !!if (size(wpr%wprojy)<size(wprojy)) stop 'size y'
  !!if (size(wpr%wprojz)<size(wprojz)) stop 'size z'
  !!if (size(wpr%wprojx,1)/=size(wprojx,1)) stop 'ncplx not equal'
  !!if (size(wpr%wprojx,2)/=size(wprojx,2)) stop 'n1 not equal'
  !allocate(wprojx(1:ncplx_w,0:n1,1:2,1:nterm))
  !allocate(wprojy(1:ncplx_w,0:n2,1:2,1:nterm))
  !allocate(wprojz(1:ncplx_w,0:n3,1:2,1:nterm))

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')
  one = real(1, gp)
  if (ncplx_g == 2) one(2) = real(0, gp)

  ! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
  !n(c) err_norm=0.0_gp 

!!  call system_clock(ncount0,ncount_rate,ncount_max)

  ! OpenMP commented here as doesn't work on Vesta
  !!$omp parallel default(shared) private(iterm,work,ml1,mu1,ml2,mu2,ml3,mu3) &
  !!$omp private(ithread,ichunk,factor,n_gau)

  !!$omp critical
    !#work = f_malloc((/ 0.to.nw, 1.to.2, 1.to.2 /),id='work')
  !allocate(work(0:nw,1:2,1:2))
  !!$omp end critical

  !!wprojx=0.d0
  !!wpr%wprojx=0.d0

  !!$ ithread=omp_get_thread_num()
  !!$ nthread=omp_get_num_threads() 
  !!$ ichunk=0
  do iterm=1,nterm
     !!$ ichunk=ichunk+1
     !!$ if (mod(ichunk,nthread).eq.ithread) then
     factor(:)=fac_arr(:,iterm)
     n_gau=lx(iterm) 
     !!call gauss_to_daub_k(hx,kx*hx,ncplx_w,ncplx_g,ncplx_k,factor,rx,gau_a,n_gau,ns1,n1,ml1,mu1,&
     !!     wprojx(1,0,1,iterm),wpr%work,nw,perx,gau_cut) 
     call gauss_to_daub_k(hx,kx*hx,ncplx_w,ncplx_g,ncplx_k,factor,rx,gau_a,n_gau,ns1,n1,ml1,mu1,&
          wpr%wproj(1),wpr%work,size(wpr%work, 1),perx,gau_cut) 
     !!$ endif
     call vcopy(ncplx_w*(n1+1), wpr%wproj(1), 1, wpr%wprojx(1,0,1,iterm), iskip)
     call vcopy(ncplx_w*(n1+1), wpr%wproj(ncplx_w*(n1+1)+1), 1, wpr%wprojx(1,0,2,iterm), iskip)

     !!$ ichunk=ichunk+1
     !!$ if (mod(ichunk,nthread).eq.ithread) then
     n_gau=ly(iterm) 
     call gauss_to_daub_k(hy,ky*hy,ncplx_w,ncplx_g,ncplx_k,one,ry,gau_a,n_gau,ns2,n2,ml2,mu2,&
          wpr%wproj(1),wpr%work,size(wpr%work, 1),pery,gau_cut) 
     !!$ endif
     call vcopy(ncplx_w*(n2+1), wpr%wproj(1), 1, wpr%wprojy(1,0,1,iterm), iskip)
     call vcopy(ncplx_w*(n2+1), wpr%wproj(ncplx_w*(n2+1)+1), 1, wpr%wprojy(1,0,2,iterm), iskip)

     !!$ ichunk=ichunk+1
     !!$ if (mod(ichunk,nthread).eq.ithread) then
     n_gau=lz(iterm) 
     call gauss_to_daub_k(hz,kz*hz,ncplx_w,ncplx_g,ncplx_k,one,rz,gau_a,n_gau,ns3,n3,ml3,mu3,&
          wpr%wproj(1),wpr%work,size(wpr%work, 1),perz,gau_cut)
     !!$ endif
     call vcopy(ncplx_w*(n3+1), wpr%wproj(1), 1, wpr%wprojz(1,0,1,iterm), iskip)
     call vcopy(ncplx_w*(n3+1), wpr%wproj(ncplx_w*(n3+1)+1), 1, wpr%wprojz(1,0,2,iterm), iskip)
  end do

  !!$omp critical
    !#call f_free(work) 
  !deallocate(work)
  !!$omp end critical
  !!$omp end parallel

  !write(10000+bigdft_mpi%iproc*10,*) wpr%wprojx
  !!write(10000+bigdft_mpi%iproc*10+1,*) wprojx

!wprojx=wpr%wprojx
!wprojy=wpr%wprojy
!wprojz=wpr%wprojz

  !the filling of the projector should be different if ncplx==1 or 2
  !split such as to avoid intensive call to if statements
!!  call system_clock(ncount1,ncount_rate,ncount_max)
!!  write(20,*) 'TIMING1:', dble(ncount1-ncount0)/dble(ncount_rate)

  if (ncplx_w==1) then
     !$omp parallel default(private) shared(mseg_c,keyv_p,keyg_p,n3,n2) &
     !$omp shared(n1,proj,wpr,mvctr_c) &
     !$omp shared(mvctr_f,mseg_f,nterm,n1p1,np)
     ! coarse part
     !$omp do 
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/np
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        i0jj=jj-i0
        wprojyz=wpr%wprojy(1,i2,1,1)*wpr%wprojz(1,i3,1,1)
        do i=i0,i1
           ind_c=i+i0jj
           proj(ind_c)=&
                wpr%wprojx(1,i,1,1)*wprojyz
  !!write(20000+bigdft_mpi%iproc*10,*) wpr%wprojx(1,i,1,1)
  !!write(20000+bigdft_mpi%iproc*10+1,*) wpr%wpr%wprojx(1,i,1,1)
        enddo
     enddo
     !$omp enddo

     ! First term: fine projector components
     ! fine part (beware of the behaviour with loop of zero size!)
     !$omp do
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/(np)
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        i0jj=7*(jj-i0-1)+mvctr_c
        wprojyz11=wpr%wprojy(1,i2,1,1)*wpr%wprojz(1,i3,1,1)
        wprojyz21=wpr%wprojy(1,i2,2,1)*wpr%wprojz(1,i3,1,1)
        wprojyz12=wpr%wprojy(1,i2,1,1)*wpr%wprojz(1,i3,2,1)
        wprojyz22=wpr%wprojy(1,i2,2,1)*wpr%wprojz(1,i3,2,1)
        do i=i0,i1
           ind_f=7*i+i0jj
           proj(ind_f+1)=wpr%wprojx(1,i,2,1)*wprojyz11
           proj(ind_f+2)=wpr%wprojx(1,i,1,1)*wprojyz21
           proj(ind_f+3)=wpr%wprojx(1,i,2,1)*wprojyz21
           proj(ind_f+4)=wpr%wprojx(1,i,1,1)*wprojyz12
           proj(ind_f+5)=wpr%wprojx(1,i,2,1)*wprojyz12
           proj(ind_f+6)=wpr%wprojx(1,i,1,1)*wprojyz22
           proj(ind_f+7)=wpr%wprojx(1,i,2,1)*wprojyz22
  !!write(30000+bigdft_mpi%iproc*10,*) wpr%wprojx(1,i,2,1)
  !!write(30000+bigdft_mpi%iproc*10+1,*) wpr%wpr%wprojx(1,i,2,1)
        enddo
     enddo
     !$omp enddo

     if (nterm >= 2) then
        ! Other terms: coarse projector components
        ! coarse part
        !$omp do 
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=jj-i0
           do i=i0,i1
              ind_c=i+i0jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+&
                      wpr%wprojx(1,i,1,iterm)*wpr%wprojy(1,i2,1,iterm)*wpr%wprojz(1,i3,1,iterm)
  !!write(40000+bigdft_mpi%iproc*10,*) wpr%wprojx(1,i,1,iterm)
  !!write(40000+bigdft_mpi%iproc*10+1,*) wpr%wpr%wprojx(1,i,1,iterm)
              end do
           end do
        end do
        !$omp enddo

        ! Other terms: fine projector components
        !$omp do
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=7*(jj-i0-1)+mvctr_c
           do i=i0,i1
              ind_f=7*i+i0jj
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+&
                      wpr%wprojx(1,i,2,iterm)*wpr%wprojy(1,i2,1,iterm)*wpr%wprojz(1,i3,1,iterm)
                 proj(ind_f+2)=proj(ind_f+2)+&
                      wpr%wprojx(1,i,1,iterm)*wpr%wprojy(1,i2,2,iterm)*wpr%wprojz(1,i3,1,iterm)
                 proj(ind_f+3)=proj(ind_f+3)+&
                      wpr%wprojx(1,i,2,iterm)*wpr%wprojy(1,i2,2,iterm)*wpr%wprojz(1,i3,1,iterm)
                 proj(ind_f+4)=proj(ind_f+4)+&
                      wpr%wprojx(1,i,1,iterm)*wpr%wprojy(1,i2,1,iterm)*wpr%wprojz(1,i3,2,iterm)
                 proj(ind_f+5)=proj(ind_f+5)+&
                      wpr%wprojx(1,i,2,iterm)*wpr%wprojy(1,i2,1,iterm)*wpr%wprojz(1,i3,2,iterm)
                 proj(ind_f+6)=proj(ind_f+6)+&
                      wpr%wprojx(1,i,1,iterm)*wpr%wprojy(1,i2,2,iterm)*wpr%wprojz(1,i3,2,iterm)
                 proj(ind_f+7)=proj(ind_f+7)+&
                      wpr%wprojx(1,i,2,iterm)*wpr%wprojy(1,i2,2,iterm)*wpr%wprojz(1,i3,2,iterm)
  !!write(50000+bigdft_mpi%iproc*10,*) wpr%wprojx(1,i,2,iterm)
  !!write(50000+bigdft_mpi%iproc*10+1,*) wpr%wpr%wprojx(1,i,2,iterm)
                 !! proj_f(1,i-i0+jj)=proj_f(1,i-i0+jj)+&
                 !!      wpr%wprojx(i,2,iterm)*wpr%wprojy(i2,1,iterm)*wpr%wprojz(i3,1,iterm)
                 !! proj_f(2,i-i0+jj)=proj_f(2,i-i0+jj)+&
                 !!      wpr%wprojx(i,1,iterm)*wpr%wprojy(i2,2,iterm)*wpr%wprojz(i3,1,iterm)
                 !! proj_f(3,i-i0+jj)=proj_f(3,i-i0+jj)+&
                 !!      wpr%wprojx(i,2,iterm)*wpr%wprojy(i2,2,iterm)*wpr%wprojz(i3,1,iterm)
                 !! proj_f(4,i-i0+jj)=proj_f(4,i-i0+jj)+&
                 !!      wpr%wprojx(i,1,iterm)*wpr%wprojy(i2,1,iterm)*wpr%wprojz(i3,2,iterm)
                 !! proj_f(5,i-i0+jj)=proj_f(5,i-i0+jj)+&
                 !!      wpr%wprojx(i,2,iterm)*wpr%wprojy(i2,1,iterm)*wpr%wprojz(i3,2,iterm)
                 !! proj_f(6,i-i0+jj)=proj_f(6,i-i0+jj)+&
                 !!      wpr%wprojx(i,1,iterm)*wpr%wprojy(i2,2,iterm)*wpr%wprojz(i3,2,iterm)
                 !! proj_f(7,i-i0+jj)=proj_f(7,i-i0+jj)+&
                 !!      wpr%wprojx(i,2,iterm)*wpr%wprojy(i2,2,iterm)*wpr%wprojz(i3,2,iterm)
              end do
           end do
        end do
        !$omp enddo
     end if
     !$omp end parallel

  else if (ncplx_w==2) then
     !$omp parallel default(private) shared(mseg_c,keyv_p,keyg_p,n3,n2,ncplx_k) &
     !$omp shared(n1,proj,wpr,mvctr_c) &
     !$omp shared(nterm,mvctr_f,mseg_f,n1p1,np,mvctr_cf,mvctr_cf2)
     !part with real and imaginary part
     !modify the openMP statements such as to benefit from parallelisation
     !Here accumulate only the REAL part,
     !The imaginary part is done below

     ! coarse part
     !$omp do
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/(np)
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        i0jj=jj-i0
        do i=i0,i1
           ind_c=i+i0jj
           proj(ind_c)=&
                re_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,1,1))
        enddo
     enddo
     !$omp end do

     ! First term: fine projector components
     ! fine part
     !$omp do 
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/(np)
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        i0jj=7*(jj-i0-1)+mvctr_c
        do i=i0,i1
           ind_f=7*i+i0jj
           proj(ind_f+1)=re_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+2)=re_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+3)=re_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+4)=re_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+5)=re_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+6)=re_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+7)=re_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,2,1))
        enddo
     enddo
     !$omp end do

     if (nterm >= 2) then
        ! Other terms: coarse projector components
        ! coarse part
        !$omp do 
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=jj-i0
           do i=i0,i1
              ind_c=i+i0jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,1,iterm))
              end do
           end do
        end do
        !$omp enddo

        ! Other terms: fine projector components
        !$omp do 
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=7*(jj-i0-1)+mvctr_c
           do i=i0,i1
              ind_f=7*i+i0jj
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+2)=proj(ind_f+2)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+3)=proj(ind_f+3)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+4)=proj(ind_f+4)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+5)=proj(ind_f+5)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+6)=proj(ind_f+6)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+7)=proj(ind_f+7)+re_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,2,iterm))
              end do
           end do
        end do
        !$omp enddo
     end if

     if(ncplx_k==2) then
        !now the imaginary part, only for complex projectors
        !when ncplx_g==2 and ncplx_k==1 the projectors are real.
        !so we skip this part.


     !$omp do 
     do iseg=1,mseg_c
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/(np)
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        i0jj=jj-i0+mvctr_cf
        do i=i0,i1
           ind_c=i+i0jj
           proj(ind_c)=&
                im_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,1,1))
        enddo
     enddo
     !$omp enddo
  
     ! First term: fine projector components
     ! fine part
     !$omp do 
     do iseg=mseg_c+1,mseg_c+mseg_f
        jj=keyv_p(iseg)
        j0=keyg_p(1,iseg)
        j1=keyg_p(2,iseg)
        ii=j0-1
        i3=ii/(np)
        ii=ii-i3*np
        i2=ii/n1p1
        i0=ii-i2*n1p1
        i1=i0+j1-j0
        !correction for xlf compiler bug
        ind_f=mvctr_cf2+7*(jj-2)
        do i=i0,i1
           ind_f=ind_f+7
           !ind_f=mvctr_c+7*mvctr_f+mvctr_c+7*(i-i0+jj-1)
           proj(ind_f+1)=im_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+2)=im_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+3)=im_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,1,1))
           proj(ind_f+4)=im_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+5)=im_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,1,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+6)=im_cmplx_prod(wpr%wprojx(1:2,i,1,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,2,1))
           proj(ind_f+7)=im_cmplx_prod(wpr%wprojx(1:2,i,2,1),wpr%wprojy(1:2,i2,2,1),wpr%wprojz(1:2,i3,2,1))
        enddo
     enddo
     !$omp enddo

     if (nterm >= 2) then
        ! Other terms: coarse projector components
        ! coarse part
        !$omp do
        do iseg=1,mseg_c
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=jj-i0+mvctr_cf
           do i=i0,i1
              ind_c=i+i0jj
              do iterm=2,nterm
                 proj(ind_c)=proj(ind_c)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,1,iterm))
              end do
           end do
        end do
        !$omp enddo

        ! Other terms: fine projector components
        !$omp do
        do iseg=mseg_c+1,mseg_c+mseg_f
           jj=keyv_p(iseg)
           j0=keyg_p(1,iseg)
           j1=keyg_p(2,iseg)
           ii=j0-1
           i3=ii/(np)
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=7*(jj-i0-1)+mvctr_cf2
           do i=i0,i1
              ind_f=7*i+i0jj
              do iterm=2,nterm
                 proj(ind_f+1)=proj(ind_f+1)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+2)=proj(ind_f+2)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+3)=proj(ind_f+3)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,1,iterm))
                 proj(ind_f+4)=proj(ind_f+4)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+5)=proj(ind_f+5)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,1,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+6)=proj(ind_f+6)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,1,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,2,iterm))
                 proj(ind_f+7)=proj(ind_f+7)+im_cmplx_prod(&
                      wpr%wprojx(1:2,i,2,iterm),wpr%wprojy(1:2,i2,2,iterm),wpr%wprojz(1:2,i3,2,iterm))
              end do
           end do
        end do
        !$omp enddo
     end if  !nterm >= 2
     end if  !ncplx_k==2
     !$omp end parallel

  end if !ncplx_w==2
!!  call system_clock(ncount2,ncount_rate,ncount_max)
!!  write(20,*) 'TIMING2:', dble(ncount2-ncount1)/dble(ncount_rate)

  !!call f_free(wprojx)
  !!call f_free(wprojy)
  !!call f_free(wprojz)
  !deallocate(wprojx)
  !deallocate(wprojy)
  !deallocate(wprojz)

!  i_all=-product(shape(work))*kind(work)
!  deallocate(work,stat=i_stat)
!  call memocc(i_stat,i_all,'work',subname)
  !call finalize_real_space_conversion()

  call f_release_routine()

contains

  !> Real part of the complex product
  pure function re_cmplx_prod(a,b,c)
    use module_base, only: wp
    implicit none
    real(wp), dimension(2), intent(in) :: a,b,c
    real(wp) :: re_cmplx_prod

    re_cmplx_prod=a(1)*(b(1)*c(1)-b(2)*c(2)) &
         -a(2)*(b(1)*c(2)+b(2)*c(1))
         

  END FUNCTION re_cmplx_prod


  !>   Imaginary part of the complex product
  pure function im_cmplx_prod(a,b,c)
    use module_base, only: wp
    implicit none
    real(wp), dimension(2), intent(in) :: a,b,c
    real(wp) :: im_cmplx_prod

    im_cmplx_prod=a(2)*(b(1)*c(1)-b(2)*c(2)) &
         +a(1)*(b(2)*c(1)+b(1)*c(2))

  END FUNCTION im_cmplx_prod

END SUBROUTINE crtproj

subroutine calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  use module_base
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(gp), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.7071067811865475244008444_gp
  else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3651483716701107423046465_gp
     fac_arr(2)=0.3651483716701107423046465_gp
     fac_arr(3)=0.3651483716701107423046465_gp
  else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.09200874124564722903948358_gp
     fac_arr(2)=0.1840174824912944580789672_gp
     fac_arr(3)=0.09200874124564722903948358_gp
     fac_arr(4)=0.1840174824912944580789672_gp
     fac_arr(5)=0.1840174824912944580789672_gp
     fac_arr(6)=0.09200874124564722903948358_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.000000000000000000000000_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3380617018914066310038473_gp
     fac_arr(2)=0.3380617018914066310038473_gp
     fac_arr(3)=0.3380617018914066310038473_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.06795295885835007261827187_gp
     fac_arr(2)=0.1359059177167001452365437_gp
     fac_arr(3)=0.06795295885835007261827187_gp
     fac_arr(4)=0.1359059177167001452365437_gp
     fac_arr(5)=0.1359059177167001452365437_gp
     fac_arr(6)=0.06795295885835007261827187_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.414213562373095048801689_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.7071067811865475244008444_gp
     fac_arr(2)=-0.7071067811865475244008444_gp
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.4082482904638630163662140_gp
     fac_arr(2)=-0.4082482904638630163662140_gp
     fac_arr(3)=0.8164965809277260327324280_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=1
     lx(2)=0 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=1
     lx(2)=1 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=0
     lx(2)=1 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3563483225498991795794046_gp
     fac_arr(2)=0.3563483225498991795794046_gp
     fac_arr(3)=0.3563483225498991795794046_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=4 ; lz(2)=0
     lx(3)=2 ; ly(3)=0 ; lz(3)=2
     lx(4)=0 ; ly(4)=2 ; lz(4)=2
     fac_arr(1)=0.1781741612749495897897023_gp
     fac_arr(2)=-0.1781741612749495897897023_gp
     fac_arr(3)=0.1781741612749495897897023_gp
     fac_arr(4)=-0.1781741612749495897897023_gp
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=-0.1028688999747279401740630_gp
     fac_arr(2)=-0.2057377999494558803481260_gp
     fac_arr(3)=-0.1028688999747279401740630_gp
     fac_arr(4)=0.1028688999747279401740630_gp
     fac_arr(5)=0.1028688999747279401740630_gp
     fac_arr(6)=0.2057377999494558803481260_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=1
     lx(2)=2 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=5 ; lz(3)=1
     lx(4)=2 ; ly(4)=1 ; lz(4)=3
     lx(5)=0 ; ly(5)=3 ; lz(5)=3
     lx(6)=0 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=1
     lx(2)=3 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=4 ; lz(3)=1
     lx(4)=3 ; ly(4)=0 ; lz(4)=3
     lx(5)=1 ; ly(5)=2 ; lz(5)=3
     lx(6)=1 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=0
     lx(2)=3 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=5 ; lz(3)=0
     lx(4)=3 ; ly(4)=1 ; lz(4)=2
     lx(5)=1 ; ly(5)=3 ; lz(5)=2
     lx(6)=1 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.05959868750235655989526993_gp
     fac_arr(2)=0.1191973750047131197905399_gp
     fac_arr(3)=0.05959868750235655989526993_gp
     fac_arr(4)=0.1191973750047131197905399_gp
     fac_arr(5)=0.1191973750047131197905399_gp
     fac_arr(6)=0.05959868750235655989526993_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=4 ; ly(5)=0 ; lz(5)=2
     lx(6)=0 ; ly(6)=4 ; lz(6)=2
     lx(7)=2 ; ly(7)=0 ; lz(7)=4
     lx(8)=0 ; ly(8)=2 ; lz(8)=4
     fac_arr(1)=0.02979934375117827994763496_gp
     fac_arr(2)=0.02979934375117827994763496_gp
     fac_arr(3)=-0.02979934375117827994763496_gp
     fac_arr(4)=-0.02979934375117827994763496_gp
     fac_arr(5)=0.05959868750235655989526993_gp
     fac_arr(6)=-0.05959868750235655989526993_gp
     fac_arr(7)=0.02979934375117827994763496_gp
     fac_arr(8)=-0.02979934375117827994763496_gp
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
     nterm=7
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=2 ; ly(5)=0 ; lz(5)=4
     lx(6)=0 ; ly(6)=2 ; lz(6)=4
     lx(7)=0 ; ly(7)=0 ; lz(7)=6
     fac_arr(1)=-0.01720465913641697233541246_gp
     fac_arr(2)=-0.05161397740925091700623738_gp
     fac_arr(3)=-0.05161397740925091700623738_gp
     fac_arr(4)=-0.01720465913641697233541246_gp
     fac_arr(5)=0.05161397740925091700623738_gp
     fac_arr(6)=0.05161397740925091700623738_gp
     fac_arr(7)=0.03440931827283394467082492_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894_gp
     fac_arr(2)=0.3162277660168379331998894_gp
     fac_arr(3)=-1.264911064067351732799557_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894_gp
     fac_arr(2)=0.3162277660168379331998894_gp
     fac_arr(3)=-1.264911064067351732799557_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.7745966692414833770358531_gp
     fac_arr(2)=0.7745966692414833770358531_gp
     fac_arr(3)=-0.5163977794943222513572354_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.4082482904638630163662140_gp
     fac_arr(2)=-1.224744871391589049098642_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-1.224744871391589049098642_gp
     fac_arr(2)=0.4082482904638630163662140_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.000000000000000000000000_gp
     fac_arr(2)=-1.000000000000000000000000_gp
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.000000000000000000000000_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506_gp
     fac_arr(2)=0.1271283452327456420595701_gp
     fac_arr(3)=0.06356417261637282102978506_gp
     fac_arr(4)=-0.1906925178491184630893552_gp
     fac_arr(5)=-0.1906925178491184630893552_gp
     fac_arr(6)=-0.2542566904654912841191402_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506_gp
     fac_arr(2)=0.1271283452327456420595701_gp
     fac_arr(3)=0.06356417261637282102978506_gp
     fac_arr(4)=-0.1906925178491184630893552_gp
     fac_arr(5)=-0.1906925178491184630893552_gp
     fac_arr(6)=-0.2542566904654912841191402_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.1556997888323045941832351_gp
     fac_arr(2)=0.3113995776646091883664703_gp
     fac_arr(3)=0.1556997888323045941832351_gp
     fac_arr(4)=0.05189992961076819806107838_gp
     fac_arr(5)=0.05189992961076819806107838_gp
     fac_arr(6)=-0.1037998592215363961221568_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
     nterm=5
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     fac_arr(1)=0.08206099398622182182282711_gp
     fac_arr(2)=-0.1641219879724436436456542_gp
     fac_arr(3)=-0.2461829819586654654684813_gp
     fac_arr(4)=0.08206099398622182182282711_gp
     fac_arr(5)=-0.2461829819586654654684813_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
     nterm=5
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     fac_arr(1)=-0.2461829819586654654684813_gp
     fac_arr(2)=-0.1641219879724436436456542_gp
     fac_arr(3)=0.08206099398622182182282711_gp
     fac_arr(4)=-0.2461829819586654654684813_gp
     fac_arr(5)=0.08206099398622182182282711_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=4 ; lz(2)=1
     lx(3)=2 ; ly(3)=0 ; lz(3)=3
     lx(4)=0 ; ly(4)=2 ; lz(4)=3
     fac_arr(1)=0.2010075630518424150978747_gp
     fac_arr(2)=-0.2010075630518424150978747_gp
     fac_arr(3)=0.2010075630518424150978747_gp
     fac_arr(4)=-0.2010075630518424150978747_gp
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=1
     lx(2)=1 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.4020151261036848301957494_gp
     fac_arr(2)=0.4020151261036848301957494_gp
     fac_arr(3)=0.4020151261036848301957494_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
     nterm=10
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     lx(10)=1 ; ly(10)=0 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687_gp
     fac_arr(2)=0.02731154967995675489524106_gp
     fac_arr(3)=0.02731154967995675489524106_gp
     fac_arr(4)=0.009103849893318918298413687_gp
     fac_arr(5)=-0.01820769978663783659682737_gp
     fac_arr(6)=-0.03641539957327567319365475_gp
     fac_arr(7)=-0.01820769978663783659682737_gp
     fac_arr(8)=-0.06372694925323242808889581_gp
     fac_arr(9)=-0.06372694925323242808889581_gp
     fac_arr(10)=-0.03641539957327567319365475_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
     nterm=10
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     lx(10)=0 ; ly(10)=1 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687_gp
     fac_arr(2)=0.02731154967995675489524106_gp
     fac_arr(3)=0.02731154967995675489524106_gp
     fac_arr(4)=0.009103849893318918298413687_gp
     fac_arr(5)=-0.01820769978663783659682737_gp
     fac_arr(6)=-0.03641539957327567319365475_gp
     fac_arr(7)=-0.01820769978663783659682737_gp
     fac_arr(8)=-0.06372694925323242808889581_gp
     fac_arr(9)=-0.06372694925323242808889581_gp
     fac_arr(10)=-0.03641539957327567319365475_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
     nterm=10
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=2 ; ly(6)=2 ; lz(6)=3
     lx(7)=0 ; ly(7)=4 ; lz(7)=3
     lx(8)=2 ; ly(8)=0 ; lz(8)=5
     lx(9)=0 ; ly(9)=2 ; lz(9)=5
     lx(10)=0 ; ly(10)=0 ; lz(10)=7
     fac_arr(1)=0.02229978693352242055222348_gp
     fac_arr(2)=0.06689936080056726165667044_gp
     fac_arr(3)=0.06689936080056726165667044_gp
     fac_arr(4)=0.02229978693352242055222348_gp
     fac_arr(5)=0.02973304924469656073629797_gp
     fac_arr(6)=0.05946609848939312147259594_gp
     fac_arr(7)=0.02973304924469656073629797_gp
     fac_arr(8)=-0.007433262311174140184074493_gp
     fac_arr(9)=-0.007433262311174140184074493_gp
     fac_arr(10)=-0.01486652462234828036814899_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
     nterm=9
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     fac_arr(1)=0.01175301967439877980816756_gp
     fac_arr(2)=-0.01175301967439877980816756_gp
     fac_arr(3)=-0.05876509837199389904083778_gp
     fac_arr(4)=-0.03525905902319633942450267_gp
     fac_arr(5)=0.02350603934879755961633511_gp
     fac_arr(6)=-0.04701207869759511923267022_gp
     fac_arr(7)=-0.07051811804639267884900533_gp
     fac_arr(8)=0.01175301967439877980816756_gp
     fac_arr(9)=-0.03525905902319633942450267_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
     nterm=9
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     fac_arr(1)=-0.03525905902319633942450267_gp
     fac_arr(2)=-0.05876509837199389904083778_gp
     fac_arr(3)=-0.01175301967439877980816756_gp
     fac_arr(4)=0.01175301967439877980816756_gp
     fac_arr(5)=-0.07051811804639267884900533_gp
     fac_arr(6)=-0.04701207869759511923267022_gp
     fac_arr(7)=0.02350603934879755961633511_gp
     fac_arr(8)=-0.03525905902319633942450267_gp
     fac_arr(9)=0.01175301967439877980816756_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=0 ; ly(6)=4 ; lz(6)=3
     lx(7)=2 ; ly(7)=0 ; lz(7)=5
     lx(8)=0 ; ly(8)=2 ; lz(8)=5
     fac_arr(1)=0.02878890113916869875409405_gp
     fac_arr(2)=0.02878890113916869875409405_gp
     fac_arr(3)=-0.02878890113916869875409405_gp
     fac_arr(4)=-0.02878890113916869875409405_gp
     fac_arr(5)=0.05757780227833739750818811_gp
     fac_arr(6)=-0.05757780227833739750818811_gp
     fac_arr(7)=0.02878890113916869875409405_gp
     fac_arr(8)=-0.02878890113916869875409405_gp
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=1
     lx(2)=3 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=5 ; lz(3)=1
     lx(4)=3 ; ly(4)=1 ; lz(4)=3
     lx(5)=1 ; ly(5)=3 ; lz(5)=3
     lx(6)=1 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05757780227833739750818811_gp
     fac_arr(2)=0.1151556045566747950163762_gp
     fac_arr(3)=0.05757780227833739750818811_gp
     fac_arr(4)=0.1151556045566747950163762_gp
     fac_arr(5)=0.1151556045566747950163762_gp
     fac_arr(6)=0.05757780227833739750818811_gp

  else
     stop 'PSP format error'
  endif
  
END SUBROUTINE calc_coeff_proj

subroutine plr_segs_and_vctrs(plr,nseg_c,nseg_f,nvctr_c,nvctr_f)
  use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: plr
  integer, intent(out) :: nseg_c,nseg_f,nvctr_c,nvctr_f
  !local variables
  
  nseg_c=plr%wfd%nseg_c
  nseg_f=plr%wfd%nseg_f
  nvctr_c=plr%wfd%nvctr_c
  nvctr_f=plr%wfd%nvctr_f

end subroutine plr_segs_and_vctrs

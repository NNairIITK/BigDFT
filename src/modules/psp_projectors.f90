!> @file
!! Datatypes and associated methods relativ s to the nonlocal projectors
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Module defining datatypes of the projectors as well as constructirs and destructors
module psp_projectors
  use module_base, only: wp,gp
  use locregs
  implicit none

  !parameters identifying the different strategy for the application of a projector 
  !in a localisation region
  integer, parameter :: PSP_APPLY_SKIP=0 !<the projector is not applied. This might happend when ilr and iat does not interact
  integer, parameter :: PSP_APPLY_MASK=1 !<use mask arrays. The mask array has to be created before.
  integer, parameter :: PSP_APPLY_KEYS=2 !<use keys. No mask nor packing. Equivalend to traditional application
  integer, parameter :: PSP_APPLY_MASK_PACK=3 !<use masking and creates a pack arrays from them. 
                                    !!Most likely this is the common usage for atoms with lots of projectors and localization regions "close" to them
  integer, parameter :: PSP_APPLY_KEYS_PACK=4 !<use keys and pack arrays. Useful especially when there is no memory to create a lot of packing arrays, 
                                    !!for example when lots of lrs interacts with lots of atoms

  !> arrays defining how a given projector and a given wavefunction descriptor should interact
  type, public :: nlpsp_to_wfd
     integer :: strategy !< can be MASK,KEYS,MASK_PACK,KEYS_PACK,SKIP
     integer :: nmseg_c !< number of segments intersecting in the coarse region
     integer :: nmseg_f !< number of segments intersecting in the fine region
     integer, dimension(:,:), pointer :: mask !<mask array of dimesion 3,nmseg_c+nmseg_f for psp applilcation
  end type nlpsp_to_wfd

  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: mproj !< number of projectors for this descriptor
     integer :: nlr !< total no. localization regions potentially interacting with the psp
     type(locreg_descriptors) :: plr !< localization region descriptor of a given projector (null if nlp=0)
     type(nlpsp_to_wfd), dimension(:), pointer :: tolr !<maskings for the locregs, dimension nlr
  end type nonlocal_psp_descriptors

  !> describe the information associated to the non-local part of Pseudopotentials
  type, public :: DFT_PSP_projectors 
     logical :: on_the_fly             !< strategy for projector creation
     integer :: nproj,nprojel,natoms   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     real(wp), dimension(:), pointer :: proj !<storage space of the projectors in wavelet basis
     type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd !<descriptor per projector, of size natom
     !most likely work arrays might go here
     
  end type DFT_PSP_projectors

  contains

    !creators
    pure function nlpsp_to_wfd_null() result(tolr)
      implicit none
      type(nlpsp_to_wfd) :: tolr
      call nullify_nlpsp_to_wfd(tolr)
    end function nlpsp_to_wfd_null
    pure subroutine nullify_nlpsp_to_wfd(tolr)
      implicit none
      type(nlpsp_to_wfd), intent(out) :: tolr
      tolr%strategy=PSP_APPLY_SKIP
      tolr%nmseg_c=0
      tolr%nmseg_f=0
      nullify(tolr%mask)
    end subroutine nullify_nlpsp_to_wfd

    pure function nonlocal_psp_descriptors_null() result(pspd)
      implicit none
      type(nonlocal_psp_descriptors) :: pspd
      call nullify_nonlocal_psp_descriptors(pspd)
    end function nonlocal_psp_descriptors_null
    pure subroutine nullify_nonlocal_psp_descriptors(pspd)
      implicit none
      type(nonlocal_psp_descriptors), intent(out) :: pspd
      pspd%mproj=0
      pspd%nlr=0
      call nullify_locreg_descriptors(pspd%plr)
      nullify(pspd%tolr)
    end subroutine nullify_nonlocal_psp_descriptors

    pure function DFT_PSP_projectors_null() result(nl)
      implicit none
      type(DFT_PSP_projectors) :: nl
      call nullify_DFT_PSP_projectors(nl)
    end function DFT_PSP_projectors_null
    pure subroutine nullify_DFT_PSP_projectors(nl)
      implicit none
      type(DFT_PSP_projectors), intent(out) :: nl
      nl%on_the_fly=.true.
      nl%nproj=0
      nl%nprojel=0
      nl%natoms=0
      nl%zerovol=100.0_gp
      nullify(nl%proj)
      nullify(nl%pspd)
    end subroutine nullify_DFT_PSP_projectors

    !allocators

    !destructors
    subroutine deallocate_nlpsp_to_wfd(tolr)
      implicit none
      type(nlpsp_to_wfd), intent(inout) :: tolr
      call f_free_ptr(tolr%mask)
    end subroutine deallocate_nlpsp_to_wfd

    subroutine deallocate_nonlocal_psp_descriptors(pspd)
      implicit none
      type(nonlocal_psp_descriptors), intent(inout) :: pspd
      !local variables
      integer :: ilr
      if (associated(pspd%tolr)) then
         do ilr=1,pspd%nlr
            call deallocate_nlpsp_to_wfd(pspd%tolr(ilr))
            call nullify_nlpsp_to_wfd(pspd%tolr(ilr))
         end do
         deallocate(pspd%tolr)
         nullify(pspd%tolr)
      end if
      call deallocate_locreg_descriptors(pspd%plr)
    end subroutine deallocate_nonlocal_psp_descriptors

    subroutine deallocate_DFT_PSP_projectors(nl)
      implicit none
      type(DFT_PSP_projectors), intent(inout) :: nl
      !local variables
      integer :: iat
      
      if (associated(nl%pspd)) then
         do iat=1,nl%natoms
            call deallocate_nonlocal_psp_descriptors(nl%pspd(iat))
         end do
         deallocate(nl%pspd)
         nullify(nl%pspd)
      end if
      call f_free_ptr(nl%proj)
    END SUBROUTINE deallocate_DFT_PSP_projectors

    subroutine free_DFT_PSP_projectors(nl)
      implicit none
      type(DFT_PSP_projectors), intent(inout) :: nl
      call deallocate_DFT_PSP_projectors(nl)
      call nullify_DFT_PSP_projectors(nl)
    end subroutine free_DFT_PSP_projectors

    !then routines which are typical of the projector application or creation follow

    !> converts the bound of the lr descriptors in local bounds of the plr locregs
    pure subroutine bounds_to_plr_limits(thatway,icoarse,plr,nl1,nl2,nl3,nu1,nu2,nu3)
      implicit none
      logical, intent(in) :: thatway !< if .true., the plr descriptors has to be filled
                         !! if .false., the nl bounds are filled from the plr
      integer, intent(in) :: icoarse !<controls whether to assign coarse or fine
                                     !!limits. Can be 1 or 2. 
                                     !!The 2 case cannot be doe before 
                                     !!the case with 1 has been filled
      type(locreg_descriptors), intent(inout) :: plr !<projectors locreg
      integer, intent(inout) :: nl1,nl2,nl3,nu1,nu2,nu3 !<lower and upper bounds of locregs

      if (thatway) then
         if (icoarse==1) then !coarse limits (to be done first)
            plr%ns1=nl1     
            plr%ns2=nl2       
            plr%ns3=nl3       

            plr%d%n1=nu1-nl1
            plr%d%n2=nu2-nl2
            plr%d%n3=nu3-nl3
         else if (icoarse == 2) then
            plr%d%nfl1=nl1-plr%ns1
            plr%d%nfl2=nl2-plr%ns2       
            plr%d%nfl3=nl3-plr%ns3

            plr%d%nfu1=nu1-plr%ns1
            plr%d%nfu2=nu2-plr%ns2
            plr%d%nfu3=nu3-plr%ns3
!!$         else
!!$            stop 'WRONG icoarse' 
         end if
      else
         if (icoarse==1) then !coarse limits
            nl1=plr%ns1
            nl2=plr%ns2
            nl3=plr%ns3

            nu1=plr%d%n1+plr%ns1
            nu2=plr%d%n2+plr%ns2
            nu3=plr%d%n3+plr%ns3
         else if (icoarse == 2) then
            nl1=plr%d%nfl1+plr%ns1
            nl2=plr%d%nfl2+plr%ns2
            nl3=plr%d%nfl3+plr%ns3

            nu1=plr%d%nfu1+plr%ns1
            nu2=plr%d%nfu2+plr%ns2
            nu3=plr%d%nfu3+plr%ns3
!!$         else
!!$            stop 'WRONG icoarse, false case' 
         end if
      end if
    end subroutine bounds_to_plr_limits

    !> Finds the size of the smallest subbox that contains a localization region made 
    !! out of atom centered spheres
    subroutine pregion_size(geocode,rxyz,radius,rmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
      !use module_base, only: gp
      implicit none
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n1,n2,n3
      real(gp), intent(in) :: hx,hy,hz,rmult,radius
      real(gp), dimension(3), intent(in) :: rxyz
      integer, intent(out) :: nl1,nu1,nl2,nu2,nl3,nu3
      !Local variables
      double precision, parameter :: eps_mach=1.d-12
      !n(c) real(kind=8) :: onem
      real(gp) :: cxmax,cymax,czmax,cxmin,cymin,czmin,rad

      rad=radius*rmult
      cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
      cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
      czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
      !n(c) onem=1.d0-eps_mach

      nl1=ceiling(real(cxmin/hx,kind=8) - eps_mach)   
      nl2=ceiling(real(cymin/hy,kind=8) - eps_mach)   
      nl3=ceiling(real(czmin/hz,kind=8) - eps_mach)   
      nu1=floor(real(cxmax/hx,kind=8) + eps_mach)  
      nu2=floor(real(cymax/hy,kind=8) + eps_mach)  
      nu3=floor(real(czmax/hz,kind=8) + eps_mach)  

      !for non-free BC the projectors are not allowed to be also outside the box
      if (geocode == 'F') then
         if (nl1 < 0)   stop 'nl1: projector region outside cell'
         if (nl2 < 0)   stop 'nl2: projector region outside cell'
         if (nl3 < 0)   stop 'nl3: projector region outside cell'
         if (nu1 > n1)   stop 'nu1: projector region outside cell'
         if (nu2 > n2)   stop 'nu2: projector region outside cell'
         if (nu3 > n3)   stop 'nu3: projector region outside cell'
      else if (geocode == 'S') then
         !correct the extremes if they run outside the box
         if (nl1 < 0 .or. nu1 > n1) then
            nl1=0
            nu1=n1
         end if
         if (nl2 < 0)   stop 'nl2: projector region outside cell'
         if (nu2 > n2)   stop 'nu2: projector region outside cell'
         if (nl3 < 0 .or. nu3 > n3) then
            nl3=0
            nu3=n3
         end if
      else if (geocode == 'P') then
         !correct the extremes if they run outside the box
         if (nl1 < 0 .or. nu1 > n1) then
            nl1=0
            nu1=n1
         end if
         if (nl2 < 0 .or. nu2 > n2) then
            nl2=0
            nu2=n2
         end if
         if (nl3 < 0 .or. nu3 > n3) then
            nl3=0
            nu3=n3
         end if
      end if

    END SUBROUTINE pregion_size

    !> routine to update the PSP descriptors as soon as the localization regions
    ! are modified
    subroutine update_nlpsp(nl,nlr,lrs,Glr)
      implicit none
      integer, intent(in) :: nlr
      type(locreg_descriptors), intent(in) :: Glr
      type(locreg_descriptors), dimension(nlr), intent(in) :: lrs
      type(DFT_PSP_projectors), intent(inout) :: nl
      !local variables
      integer :: nbseg_dim,nkeyg_dim,iat,ilr
      integer, dimension(:), allocatable :: nbsegs_cf,keyg_lin

      !find allocating dimensions for work arrays
      nbseg_dim=0
      do iat=1,nl%natoms
         nbseg_dim=max(nbseg_dim,&
              nl%pspd(iat)%plr%wfd%nseg_c+nl%pspd(iat)%plr%wfd%nseg_f)
      end do
      nkeyg_dim=0
      do ilr=1,nlr
         nkeyg_dim=max(nkeyg_dim,lrs(ilr)%wfd%nseg_c+lrs(ilr)%wfd%nseg_f)
      end do

      !allocate the work arrays for building tolr array of structures
      nbsegs_cf=f_malloc(nbseg_dim,id='nbsegs_cf')
      keyg_lin=f_malloc(nkeyg_dim,id='keyg_lin')

      !reconstruct the projectors for any of the atoms
      do iat=1,nl%natoms
         
         !free the pre-existing array of structures
         if (associated(nl%pspd(iat)%tolr)) then
            do ilr=1,nl%pspd(iat)%nlr
               call deallocate_nlpsp_to_wfd(nl%pspd(iat)%tolr(ilr))
               call nullify_nlpsp_to_wfd(nl%pspd(iat)%tolr(ilr))
            end do
            deallocate(nl%pspd(iat)%tolr)
            nullify(nl%pspd(iat)%tolr)
         end if
         !then fill it again
         nl%pspd(iat)%nlr=nlr
         call set_nlpsp_to_wfd(Glr,nl%pspd(iat)%plr,&
              keyg_lin,nbsegs_cf,nl%pspd(iat)%tolr,lrs)
      end do

      call f_free(keyg_lin)
      call f_free(nbsegs_cf)

    end subroutine update_nlpsp

    !> initialize the information for matching the localisation region
    !! of each projector to all the localisation regions of the system
    subroutine set_nlpsp_to_wfd(Glr,plr,keyag_lin_cf,nbsegs_cf,tolr,lrs)
      implicit none
      type(locreg_descriptors), intent(in) :: Glr !<global simulation domain
      type(locreg_descriptors), intent(in) :: plr !<locreg of the projector
      !>work array: needed to build the mask array
      integer, dimension(plr%wfd%nseg_c+plr%wfd%nseg_f), intent(inout) :: nbsegs_cf
      !>work array: needed to put the unstrided keygloc of all locregs
      !! the dimension has to be maxval(lrs(:)%nseg_c+lrs(:)%nseg_f)
      integer, dimension(*), intent(inout) :: keyag_lin_cf
      !>structures which have to be filled to prepare projector applications
      type(nlpsp_to_wfd), dimension(:), pointer :: tolr 
      !> descriptors of all the localization regions of the simulation domain
      !! susceptible to interact with the projector
      type(locreg_descriptors), dimension(:), optional, intent(in) :: lrs
      !local variables
      logical :: overlap
      integer :: ilr,nlr

      nlr=1
      overlap=.true.
      if (present(lrs)) nlr=size(lrs)

      if (nlr <=0) return
      !allocate the pointer with the good size
      allocate(tolr(nlr))
      !then for any of the localization regions check the strategy
      !for applying the projectors
      do ilr=1,nlr
         !this will set to PSP_APPLY_SKIP the projector application
         call nullify_nlpsp_to_wfd(tolr(ilr))
         !now control if the projector overlaps with this locreg
         if (present(lrs)) then
            call check_overlap(lrs(ilr),plr,Glr,overlap)
            !if there is overlap, activate the strategy for the application
            if (overlap) then
               !calculate the size of the mask array
               call vcopy(lrs(ilr)%wfd%nseg_c+lrs(ilr)%wfd%nseg_f,&
                    lrs(ilr)%wfd%keyglob(1,1),2,keyag_lin_cf(1),1)
               call to_zero(plr%wfd%nseg_c+plr%wfd%nseg_f,nbsegs_cf(1))
               call mask_sizes(lrs(ilr)%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf,&
                    tolr(ilr)%nmseg_c,tolr(ilr)%nmseg_f)
               !then allocate and fill it
               tolr(ilr)%mask=&
                    f_malloc0_ptr((/3,tolr(ilr)%nmseg_c+tolr(ilr)%nmseg_f/),&
                    id='mask')
               !and filled
               call init_mask(lrs(ilr)%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf,&
                    tolr(ilr)%nmseg_c,tolr(ilr)%nmseg_f,tolr(ilr)%mask)
            end if
         else
            !calculate the size of the mask array
            call vcopy(Glr%wfd%nseg_c+Glr%wfd%nseg_f,&
                 Glr%wfd%keyglob(1,1),2,keyag_lin_cf(1),1)
            call to_zero(plr%wfd%nseg_c+plr%wfd%nseg_f,nbsegs_cf(1))
            call mask_sizes(Glr%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf,&
                 tolr(ilr)%nmseg_c,tolr(ilr)%nmseg_f)
            !then allocate and fill it
            tolr(ilr)%mask=&
                 f_malloc0_ptr((/3,tolr(ilr)%nmseg_c+tolr(ilr)%nmseg_f/),&
                 id='mask')
            !and filled
            call init_mask(Glr%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf,&
                 tolr(ilr)%nmseg_c,tolr(ilr)%nmseg_f,tolr(ilr)%mask)
         end if
            !then the best strategy can be decided according to total number of 
            !common points
            !complete stategy, the packing array is created after first projector
         if (overlap) tolr(ilr)%strategy=PSP_APPLY_MASK_PACK
            !masking is used but packing is not created, 
            !useful when only one projector has to be applied
            !tolr(ilr)%strategy=PSP_APPLY_MASK
            !old scheme, even though mask arrays is created it is never used.
            !most likely this scheme is useful for debugging purposes
            !tolr(ilr)%strategy=PSP_APPLY_KEYS
      end do

    end subroutine set_nlpsp_to_wfd

    !>find the size of the mask array for a given couple plr - llr
    subroutine mask_sizes(wfd_w,wfd_p,keyag_lin_cf,nbsegs_cf,nmseg_c,nmseg_f)
      
      implicit none
      type(wavefunctions_descriptors), intent(in) :: wfd_w,wfd_p
      !> array of the unstrided keyglob starting points of wfd_w, pre-filled
      integer, dimension(wfd_w%nseg_c+wfd_w%nseg_f), intent(in) :: keyag_lin_cf
      !> number of common segments of the wfd_w for each of the segment of wfd_p.
      !! should be initialized to zero at input
      integer, dimension(wfd_p%nseg_c+wfd_p%nseg_f), intent(inout) :: nbsegs_cf
      integer, intent(out) :: nmseg_c,nmseg_f

      call count_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,keyag_lin_cf(1),&
           wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),nbsegs_cf(1))
      !  print *,'no of points',sum(nbsegs_cf),wfd_w%nseg_c,wfd_p%nseg_c
      call integrate_nseg(wfd_p%nseg_c,nbsegs_cf(1),nmseg_c)
      !  print *,'no of points',nmseg_c

      if (wfd_w%nseg_f >0 .and. wfd_p%nseg_f > 0 ) then
         call count_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,keyag_lin_cf(wfd_w%nseg_c+1),&
              wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
              nbsegs_cf(wfd_p%nseg_c+1))
         call integrate_nseg(wfd_p%nseg_f,nbsegs_cf(wfd_p%nseg_c+1),nmseg_f)
      else
         nmseg_f=0
      end if

    contains

      !> count the total number of segments and define the integral array of displacements
      pure subroutine integrate_nseg(mseg,msegs,nseg_tot)
        implicit none
        integer, intent(in) :: mseg
        integer, dimension(mseg), intent(inout) :: msegs
        integer, intent(out) :: nseg_tot
        !local variables
        integer :: iseg,jseg

        nseg_tot=0
        do iseg=1,mseg
           jseg=msegs(iseg)
           msegs(iseg)=nseg_tot
           nseg_tot=nseg_tot+jseg
        end do
      end subroutine integrate_nseg

    end subroutine mask_sizes

    !>fill the mask array which has been previoulsly allocated and cleaned
    subroutine init_mask(wfd_w,wfd_p,keyag_lin_cf,nbsegs_cf,nmseg_c,nmseg_f,mask)
      implicit none
      integer, intent(in) :: nmseg_c,nmseg_f
      type(wavefunctions_descriptors), intent(in) :: wfd_w,wfd_p
      !> array of the unstrided keyglob starting points of wfd_w, pre-filled
      integer, dimension(wfd_w%nseg_c+wfd_w%nseg_f), intent(in) :: keyag_lin_cf
      !> number of common segments of the wfd_w for each of the segment of wfd_p.
      !! should be created by mask_sizes routine
      integer, dimension(wfd_p%nseg_c+wfd_p%nseg_f), intent(in) :: nbsegs_cf
      !>masking array. On output, it indicates for any of the segments 
      !which are common between the wavefunction and the projector
      !the starting positions in the packed arrays of projectors and wavefunction
      !respectively
      integer, dimension(3,nmseg_c+nmseg_f), intent(inout) :: mask

      call fill_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,nmseg_c,&
           nbsegs_cf(1),keyag_lin_cf(1),wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),&
           wfd_w%keyvglob(1),wfd_p%keyvglob(1),mask(1,1))
      if (nmseg_f > 0) then
         call fill_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,nmseg_f,&
              nbsegs_cf(wfd_p%nseg_c+1),keyag_lin_cf(wfd_w%nseg_c+1),&
              wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
              wfd_w%keyvglob(wfd_w%nseg_c+1),wfd_p%keyvglob(wfd_p%nseg_c+1),&
              mask(1,nmseg_c+1))
      end if

    end subroutine init_mask

  end module psp_projectors

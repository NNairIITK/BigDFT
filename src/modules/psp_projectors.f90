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
  integer, parameter :: PSP_APPLY_KEYS=2 !<use keys. No mask nor packing. Evuivalend to traditional application
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


  end module psp_projectors

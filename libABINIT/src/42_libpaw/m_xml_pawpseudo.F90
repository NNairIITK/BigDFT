!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xml_pawpseudo
!! NAME
!! m_xml_pseudo
!!
!! FUNCTION
!! This module reads a pseudopotential file written in XML.
!! A full example of the building up of a data structure using
!! the SAX paradigm.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2013 ABINIT group (FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! 
!! OUTPUT
!!
!! PARENTS
!! 
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"


module m_xml_pawpseudo

 use defs_basis
 use m_profiling
! use m_errors
use interfaces_12_hide_mpi
use interfaces_14_hidewrite
use interfaces_16_hideleave
 use m_xml_pawpseudo_types
!#if defined HAVE_TRIO_FOX
! use fox_sax
!#endif

implicit none

private


!
! It defines the routines that are called from xml_parser in response
! to particular events.
!

public  ::  destroy_paw_setup, nullify_paw_setup
public  ::  copy_paw_setup
public  ::  rdpawpsxml
private ::  paw_rdfromline
integer, public, allocatable :: ipsp2xml(:)
integer, public :: npsp_pawxml
type(paw_setup_t), public, target, allocatable,save :: paw_setup(:)
type(paw_setup_t), public,target,save ::paw_setuploc
!#if defined HAVE_TRIO_FOX
!public  :: paw_begin_element1, paw_end_element1,pawdata_chunk
!private :: ddie_
!logical, private  :: in_valenceStates = .false.,in_data=.false.
!logical, private  :: in_generator =.false.
!integer, private, save  :: ndata
!integer, private  :: ii,ival,igrid,ishpf,lmax,mesh_size
!
!! Pointers to make it easier to manage the data
!type(radialfunc_t), private, pointer  :: rp
!type(state_t), private, pointer   :: valstate (:)
!type(radial_grid_t), private, pointer   :: grids (:)
!type(radialfunc_t),private,pointer :: shpf(:)
!#endif

CONTAINS  !===========================================================
!!***
!#if defined HAVE_TRIO_FOX
!!!****f* m_xml_pawpseudo/paw_begin_element1
!!! NAME
!!! begin_element
!!!
!!! FUNCTION
!!!  Read an XML tag with a given name.
!!!  Fills the present module private data.
!!!
!!! INPUTS
!!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!!  localName = local equivalent of tag name?? Not used.
!!!  name = name of XML tag which has been read in
!!!  attributes = attributes of XML tag
!!!
!!! OUTPUT
!!!  Fills private data in present module.
!!!
!!! PARENTS
!!!
!!! CHILDREN
!!!
!!! SOURCE
!subroutine paw_begin_element1(namespaceURI,localName,name,attributes)
!
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'paw_begin_element1'
!!End of the abilint section
!
!character(len=*),intent(in)   :: namespaceURI,localName,name
!type(dictionary_t),intent(in) :: attributes
!
!character(len=100)  :: value
!integer ::iaewf=0,iproj=0,ipswf=0
!!Just to fool abirules
! value=localName
! value=namespaceURI
!
!select case(name)
!
!      case ("paw_setup")
!        paw_setuploc%tread=.true.
!        igrid=0;ishpf=0
!        ABI_DATATYPE_ALLOCATE(grids,(10))
!        ABI_DATATYPE_ALLOCATE(shpf,(7))
!        value = getValue(attributes,"version")
!!         if (value = "0.6") then
!!            write(std_out,*) "Processing a PSEUDO version 0.6 XML file"
!!          else
!!            write(std_out,*) "Can only work with PSEUDO version 0.6 XML files"
!!            STOP
!!         end if
!        paw_setuploc%version=trim(value)
!
!
!      case ("atom")
!         paw_setuploc%atom%tread=.true.
!         value = getValue(attributes,"symbol")
!         if (value == "" ) call ddie_("Cannot determine atomic symbol")
!         paw_setuploc%atom%symbol = trim(value)
!
!         value = getValue(attributes,"Z") 
!         if (value == "" ) call ddie_("Cannot determine znucl")
!         read(unit=value,fmt=*) paw_setuploc%atom%znucl
!
!         value = getValue(attributes,"core")
!         if (value == "" ) call ddie_("Cannot determine zion")
!         read(unit=value,fmt=*) paw_setuploc%atom%zion
!
!         value = getValue(attributes,"valence")
!         if (value == "" ) call ddie_("Cannot determine zval")
!         read(unit=value,fmt=*) paw_setuploc%atom%zval
!
!
!      case ("xc_functional")
!         paw_setuploc%xc_functional%tread=.true.
!         value = getValue(attributes,"type")
!         if (value == "" ) &
!            call ddie_("Cannot determine xc-functional-type")
!         paw_setuploc%xc_functional%functionaltype = trim(value)
!
!         value = getValue(attributes,"name")
!         if (value == "" ) &
!&            call ddie_("Cannot determine xc-functional-name ")
!         paw_setuploc%xc_functional%name= trim(value)
!
!      case ("generator")
!         paw_setuploc%generator%tread=.true.
!         in_generator =.true.
!         value = getValue(attributes,"type")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%generator%gen = trim(value)
!
!         value = getValue(attributes,"name")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%generator%name = trim(value)
!
!
!      case ("valence_states")
!         paw_setuploc%valence_states%tread=.true.
!         in_valenceStates=.true.
!         ival=0
!         lmax=0
!         ABI_DATATYPE_ALLOCATE(valstate,(50))
!
!      case ("state")
!         ival=ival+1
!
!         value = getValue(attributes,"n")
!         if (value == "" ) then 
!           valstate(ival)%nn=-1    
!         else
!           read(unit=value,fmt=*) valstate(ival)%nn
!         end if
! 
!         value = getValue(attributes,"l")
!         if (value == "" ) call ddie_("Cannot determine l")
!         read(unit=value,fmt=*) valstate(ival)%ll
!         if(valstate(ival)%ll>lmax) lmax=valstate(ival)%ll
!
!         value = getValue(attributes,"f")
!         if (value == "" ) then 
!           valstate(ival)%ff=-1.d0
!         else
!           read(unit=value,fmt=*) valstate(ival)%ff
!         end if
!
!         value = getValue(attributes,"rc")
!         if (value == "" ) call ddie_("Cannot determine rc")
!         read(unit=value,fmt=*) valstate(ival)%rc
!
!         value = getValue(attributes,"e")
!         if (value == "" ) call ddie_("Cannot determine e")
!         read(unit=value,fmt=*) valstate(ival)%ee
!
!         value = getValue(attributes,"id")
!         if (value == "" ) value = "unknown"
!         valstate(ival)%id = trim(value)
!
!      case ("radial_grid")
!         igrid=igrid+1
!         value = getValue(attributes,"eq")
!         if (value == "" ) value = "unknown"
!         grids(igrid)%eq = trim(value)
!           
!         value = getValue(attributes,"a")
!         if (value == "" ) then
!           grids(igrid)%aa=0.d0
!         else
!           read(unit=value,fmt=*) grids(igrid)%aa
!         end if
!
!         value = getValue(attributes,"n")
!         if (value == "" ) then
!           grids(igrid)%nn=0
!         else
!           read(unit=value,fmt=*) grids(igrid)%nn
!         end if
!
!         value = getValue(attributes,"d")
!         if (value == "" ) then
!           grids(igrid)%dd=0.d0
!         else
!           read(unit=value,fmt=*) grids(igrid)%dd
!         end if
!
!         value = getValue(attributes,"b")
!         if (value == "" ) then
!           grids(igrid)%bb=0.d0
!         else
!           read(unit=value,fmt=*) grids(igrid)%bb
!         end if
!
!         value = getValue(attributes,"istart")
!         if (value == "" ) call ddie_("Cannot determine istart")
!         read(unit=value,fmt=*) grids(igrid)%istart
!
!         value = getValue(attributes,"iend")
!         if (value == "" ) call ddie_("Cannot determine iend")
!         read(unit=value,fmt=*) grids(igrid)%iend
!
!         value = getValue(attributes,"id")
!         if (value == "" ) value = "unknown"
!         grids(igrid)%id = trim(value)
!
!end select
!
!select case(name)
!      case ("shape_function")
!         paw_setuploc%shape_function%tread=.true.
!         value = getValue(attributes,"type")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%shape_function%gtype = trim(value)
!
!         value = getValue(attributes,"grid")
!         paw_setuploc%shape_function%grid=trim(value)
!         if (value /= "" ) then
!           paw_setuploc%shape_function%gtype ="num" 
!           do ii=1,igrid
!             if(trim(paw_setuploc%shape_function%grid)==trim(grids(ii)%id)) then
!               mesh_size=grids(ii)%iend-grids(ii)%istart+1
!             end if
!           end do
!           ishpf=ishpf+1
!           ABI_ALLOCATE(shpf(ishpf)%data,(mesh_size))
!           rp=>shpf(ishpf)
!           in_data=.true.
!           ndata = 0
!         end if
!
!         value = getValue(attributes,"rc")
!         if (value == "" ) then
!           if(paw_setuploc%shape_function%gtype /="num") call ddie_("Cannot determine rc")
!         else
!           read(unit=value,fmt=*) paw_setuploc%shape_function%rc
!         end if
!
!         value = getValue(attributes,"lamb")
!         if (value == "" ) then
!           paw_setuploc%shape_function%lamb=0
!         else
!           read(unit=value,fmt=*) paw_setuploc%shape_function%lamb
!         end if
!
!      case ("pseudo_partial_wave")
!         ipswf=ipswf+1
!         paw_setuploc%pseudo_partial_wave(ipswf)%tread=.true.
!         value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%idgrid = trim(value)
!         paw_setuploc%pseudo_partial_wave(ipswf)%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) call ddie_("Cannot determine pseudo_partial_wave state")
!         paw_setuploc%pseudo_partial_wave(ipswf)%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%pseudo_partial_wave(ipswf)%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%pseudo_partial_wave(ipswf)%data,(mesh_size))
!         rp=>paw_setuploc%pseudo_partial_wave(ipswf)
!         if(ipswf==paw_setuploc%valence_states%nval) ipswf=0
!         in_data=.true.
!         ndata = 0
!
!      case ("ae_partial_wave")
!         iaewf=iaewf+1
!         paw_setuploc%ae_partial_wave(iaewf)%tread=.true.
!         value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%ae_partial_wave(iaewf)%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) call ddie_("Cannot determine ae_partial_wave state")
!         paw_setuploc%ae_partial_wave(iaewf)%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%ae_partial_wave(iaewf)%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%ae_partial_wave(iaewf)%data,(mesh_size))
!         rp=>paw_setuploc%ae_partial_wave(iaewf)
!         if(iaewf==paw_setuploc%valence_states%nval) iaewf=0
!         in_data=.true.
!         ndata = 0
!
!      case ("projector_function")
!         iproj=iproj+1
!         paw_setuploc%projector_function(iproj)%tread=.true.
!         value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%projector_function(iproj)%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) call ddie_("Cannot determine projector_function state")
!         paw_setuploc%projector_function(iproj)%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%projector_function(iproj)%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%projector_function(iproj)%data,(mesh_size))
!         rp=>paw_setuploc%projector_function(iproj)
!         if(iproj==paw_setuploc%valence_states%nval) iproj=0
!         in_data=.true.
!         ndata = 0
!
!     case ("ae_core_density")
!         paw_setuploc%ae_core_density%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%ae_core_density%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%ae_core_density%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%ae_core_density%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%ae_core_density%data,(mesh_size))
!         rp=>paw_setuploc%ae_core_density
!         in_data=.true.
!         ndata = 0
!
!     case ("pseudo_core_density")
!         paw_setuploc%pseudo_core_density%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_core_density%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_core_density%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%pseudo_core_density%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%pseudo_core_density%data,(mesh_size))
!         rp=>paw_setuploc%pseudo_core_density
!         in_data=.true.
!         ndata = 0
!
!     case ("pseudo_valence_density")
!         paw_setuploc%pseudo_valence_density%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_valence_density%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_valence_density%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%pseudo_valence_density%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%pseudo_valence_density%data,(mesh_size))
!         rp=>paw_setuploc%pseudo_valence_density
!         in_data=.true.
!         ndata = 0
!
!     case ("zero_potential")
!         paw_setuploc%zero_potential%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%zero_potential%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%zero_potential%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%zero_potential%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%zero_potential%data,(mesh_size))
!         rp=>paw_setuploc%zero_potential
!         in_data=.true.
!         ndata = 0
!
!     case ("ae_core_kinetic_energy_density")
!         paw_setuploc%ae_core_kinetic_energy_density%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%ae_core_kinetic_energy_density%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%ae_core_kinetic_energy_density%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%ae_core_kinetic_energy_density%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%ae_core_kinetic_energy_density%data,(mesh_size))
!         rp=>paw_setuploc%ae_core_kinetic_energy_density
!         in_data=.true.
!         ndata = 0
!
!     case ("pseudo_core_kinetic_energy_density")
!         paw_setuploc%pseudo_core_kinetic_energy_density%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_core_kinetic_energy_density%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%pseudo_core_kinetic_energy_density%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%pseudo_core_kinetic_energy_density%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%pseudo_core_kinetic_energy_density%data,(mesh_size))
!         rp=>paw_setuploc%pseudo_core_kinetic_energy_density
!         in_data=.true.
!         ndata = 0
!
!     case ("kresse_joubert_local_ionic_potential")
!         paw_setuploc%kresse_joubert_local_ionic_potential%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%kresse_joubert_local_ionic_potential%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%kresse_joubert_local_ionic_potential%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%kresse_joubert_local_ionic_potential%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%kresse_joubert_local_ionic_potential%data,(mesh_size))
!         rp=>paw_setuploc%kresse_joubert_local_ionic_potential
!         in_data=.true.
!         ndata = 0
!
!     case ("blochl_local_ionic_potential")
!         paw_setuploc%blochl_local_ionic_potential%tread=.true.
!          value = getValue(attributes,"grid")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%blochl_local_ionic_potential%grid=trim(value)
!
!         value = getValue(attributes,"state")
!         if (value == "" ) value = "unknown"
!         paw_setuploc%blochl_local_ionic_potential%state=trim(value)
!
!         do ii=1,igrid
!           if(trim(paw_setuploc%blochl_local_ionic_potential%grid)==trim(grids(ii)%id)) then
!             mesh_size=grids(ii)%iend-grids(ii)%istart+1
!           end if
!         end do
!
!         ABI_ALLOCATE(paw_setuploc%blochl_local_ionic_potential%data,(mesh_size))
!         rp=>paw_setuploc%blochl_local_ionic_potential
!         in_data=.true.
!         ndata = 0
!
!    case ("kinetic_energy_differences")
!         paw_setuploc%kinetic_energy_differences%tread=.true.
!         mesh_size=paw_setuploc%valence_states%nval*paw_setuploc%valence_states%nval
!         ABI_ALLOCATE(paw_setuploc%kinetic_energy_differences%data,(mesh_size))
!         rp=>paw_setuploc%kinetic_energy_differences
!         in_data=.true.
!         ndata = 0
!
!end select
!
!
!
!end subroutine paw_begin_element1
!!!***
!
!!!****f* m_xml_pawpseudo/paw_end_element1
!!! NAME
!!! end_element
!!!
!!! FUNCTION
!!!  End XML tag effect: switches flags in private data of this module
!!!
!!! INPUTS
!!!  namespaceURI = universal resource indicator for XML namespace??? Not used.
!!!  localName = local equivalent of tag name?? Not used.
!!!  name = name of XML tag which has been read in
!!!
!!! OUTPUT
!!!  side effect: private data flags in present module are turned to .false.
!!!
!!! PARENTS
!!!
!!! CHILDREN
!!!
!!! SOURCE
!subroutine paw_end_element1(namespaceURI,localName,name)
!
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'paw_end_element1'
!!End of the abilint section
!
!character(len=*),intent(in) :: namespaceURI,localName,name
!
!character(len=100)  :: value
!
!!Just to fool abirules 
! value=localName
! value=namespaceURI
! 
!select case(name)
!
!      case ("generator")
!         in_generator = .false.
!
!      case ("valence_states")
!        in_valenceStates = .false.
!        if(ival>50) call ddie_("ival>50")
!        if(ival>0)then
!          ABI_DATATYPE_ALLOCATE(paw_setuploc%valence_states%state,(ival))
!          paw_setuploc%valence_states%state(ival)%tread=.true.
!          paw_setuploc%valence_states%nval=ival
!          do ii=1,ival
!            paw_setuploc%valence_states%state(ii)=valstate(ii)
!          end do
!        end if
!        ABI_DATATYPE_DEALLOCATE(valstate)
!        if(.not.associated(paw_setuploc%ae_partial_wave)) then
!          ABI_DATATYPE_ALLOCATE(paw_setuploc%ae_partial_wave,(paw_setuploc%valence_states%nval))
!        end if
!        if(.not.associated(paw_setuploc%pseudo_partial_wave)) then
!          ABI_DATATYPE_ALLOCATE(paw_setuploc%pseudo_partial_wave,(paw_setuploc%valence_states%nval))
!        end if
!        if(.not.associated(paw_setuploc%projector_function)) then
!          ABI_DATATYPE_ALLOCATE(paw_setuploc%projector_function,(paw_setuploc%valence_states%nval))
!        end if
!
!      case ("paw_setup")
!        if(igrid>10) call ddie_("igrid>10")
!        ABI_DATATYPE_ALLOCATE(paw_setuploc%radial_grid,(igrid))
!        paw_setuploc%radial_grid(igrid)%tread=.true.
!        paw_setuploc%ngrid=igrid
!        do ii=1,igrid
!          paw_setuploc%radial_grid(ii)=grids(ii)
!        end do
!        ABI_DATATYPE_DEALLOCATE(grids)
!        do ii=1,igrid
!          if(trim(paw_setuploc%shape_function%grid)==trim(paw_setuploc%radial_grid(ii)%id)) then
!            mesh_size=paw_setuploc%radial_grid(ii)%iend-paw_setuploc%radial_grid(ii)%istart+1
!          end if
!        end do
!        if(ishpf>10) call ddie_("ishpf>7")
!        ABI_ALLOCATE(paw_setuploc%shape_function%data,(mesh_size,ishpf))
!        do ii=1,ishpf
!          paw_setuploc%shape_function%data(:,ii)=shpf(ii)%data(:)
!          ABI_DEALLOCATE(shpf(ii)%data)
!        end do
!        ABI_DATATYPE_DEALLOCATE(shpf)
!
!      case ("shape_function")
!        in_data=.false.
!
!      case ("pseudo_partial_wave")
!        in_data=.false.
!
!      case ("ae_partial_wave")
!        in_data=.false.
!
!      case ("projector_function")
!        in_data=.false.
!
!      case ("ae_core_density")
!        in_data=.false.
!
!      case ("pseudo_core_density")
!        in_data=.false.
!
!      case ("pseudo_valence_density")
!        in_data=.false.
!
!      case ("zero_potential")
!        in_data=.false.
!
!      case ("ae_core_kinetic_energy_density")
!        in_data=.false.
!
!      case ("pseudo_core_kinetic_energy_density")
!        in_data=.false.
!
!      case ("kresse_joubert_local_ionic_potential")
!        in_data=.false.
!
!      case ("blochl_local_ionic_potential")
!        in_data=.false.
!
!      case ("kinetic_energy_differences")
!        in_data=.false.
!
!
!end select
!
!end subroutine paw_end_element1
!!!***
!
!!!****f* m_xml_pseudo/pawdata_chunk
!!! NAME
!!! pawdata_chunk
!!!
!!! FUNCTION
!!!   Read in data from XML structure, if we are in a valid data field
!!!   for the present XML structures
!!!
!!! INPUTS
!!!   chunk = raw data for chunk of XML data
!!!
!!! OUTPUT
!!!   copied and translated into module data (side effect)
!!!
!!! PARENTS
!!!
!!! CHILDREN
!!!
!!! SOURCE
!subroutine pawdata_chunk(chunk)
!
!use m_xml_converters
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'pawdata_chunk'
!!End of the abilint section
!
!character(len=*),intent(in) :: chunk
!
!
!if (len_trim(chunk) == 0) RETURN     ! skip empty chunk
!
!if (in_data) then
!!
!! Note that we know where we need to put it through the pointer rp...
!!
!
!     call build_data_array(chunk,rp%data,ndata)
!
!end if
!
!end subroutine pawdata_chunk
!!!***
!
!!!****f* m_xml_pawpseudo/ddie_
!!! NAME
!!! die_
!!!
!!! FUNCTION
!!!  If there is an error in reading of XML file, stop with an error message
!!!
!!! INPUTS
!!!   str = error message string, for output before dying.
!!!
!!! OUTPUT
!!!
!!! TODO
!!!  Use a generic abinit routine instead of this local one.
!!!
!!! PARENTS
!!!      m_xml_pseudo
!!!
!!! CHILDREN
!!!
!!! SOURCE
!subroutine ddie_(str)
!
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'ddie_'
!!End of the abilint section
!
!      character(len=*), intent(in), optional   :: str
!      if (present(str)) then
!         write(unit=0,fmt="(a)") trim(str)
!      end if
!      write(unit=0,fmt="(a)") "Stopping Program"
!      stop
!end subroutine ddie_
!!!***
!#endif
!!****f* ABINIT/destroy_paw_setup
!! NAME
!! destroy_paw_setup
!!
!! FUNCTION
!!  Destroy a paw_setup datastructure
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2013 ABINIT group (FJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  paw_setup<paw_setup_type>=Datatype gathering information on XML paw setup.
!!
!! PARENTS
!!      abinit,inpspheads
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_paw_setup(paw_setupin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_paw_setup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(inout) :: paw_setupin

!Local variables-------------------------------
 integer :: ii
! *********************************************************************
 
 paw_setupin%tread=.false.
 paw_setupin%atom%tread=.false.
 paw_setupin%xc_functional%tread=.false.
 paw_setupin%generator%tread=.false.
 paw_setupin%valence_states%tread=.false.
 paw_setupin%shape_function%tread=.false.
 paw_setupin%ae_core_density%tread=.false.
 paw_setupin%pseudo_core_density%tread=.false.
 paw_setupin%pseudo_valence_density%tread=.false.
 paw_setupin%zero_potential%tread=.false.
 paw_setupin%ae_core_kinetic_energy_density%tread=.false.
 paw_setupin%pseudo_core_kinetic_energy_density%tread=.false.
 paw_setupin%kresse_joubert_local_ionic_potential%tread=.false.
 paw_setupin%blochl_local_ionic_potential%tread=.false.
 paw_setupin%kinetic_energy_differences%tread=.false.

 if(associated( paw_setupin%shape_function%data)) then
   ABI_DEALLOCATE(paw_setupin%shape_function%data)
   nullify(paw_setupin%shape_function%data)
 end if
 if(associated( paw_setupin%ae_core_density%data)) then
   ABI_DEALLOCATE(paw_setupin%ae_core_density%data)
   nullify(paw_setupin%ae_core_density%data)
 end if
 if(associated( paw_setupin%pseudo_core_density%data)) then
   ABI_DEALLOCATE(paw_setupin%pseudo_core_density%data)
   nullify(paw_setupin%pseudo_core_density%data)
 end if
 if(associated( paw_setupin%pseudo_valence_density%data)) then
   ABI_DEALLOCATE(paw_setupin%pseudo_valence_density%data)
   nullify(paw_setupin%pseudo_valence_density%data)
 end if
 if(associated( paw_setupin%zero_potential%data)) then
   ABI_DEALLOCATE(paw_setupin%zero_potential%data)
   nullify(paw_setupin%zero_potential%data)
 end if
 if(associated( paw_setupin%ae_core_kinetic_energy_density%data)) then
   ABI_DEALLOCATE(paw_setupin%ae_core_kinetic_energy_density%data)
   nullify(paw_setupin%ae_core_kinetic_energy_density%data)
 end if
 if(associated( paw_setupin%pseudo_core_kinetic_energy_density%data)) then
   ABI_DEALLOCATE(paw_setupin%pseudo_core_kinetic_energy_density%data)
   nullify(paw_setupin%pseudo_core_kinetic_energy_density%data)
 end if
 if(associated( paw_setupin%kresse_joubert_local_ionic_potential%data)) then
   ABI_DEALLOCATE(paw_setupin%kresse_joubert_local_ionic_potential%data)
   nullify(paw_setupin%kresse_joubert_local_ionic_potential%data)
 end if
 if(associated( paw_setupin%blochl_local_ionic_potential%data)) then
   ABI_DEALLOCATE(paw_setupin%blochl_local_ionic_potential%data)
   nullify(paw_setupin%blochl_local_ionic_potential%data)
 end if
 if(associated( paw_setupin%kinetic_energy_differences%data)) then
   ABI_DEALLOCATE(paw_setupin%kinetic_energy_differences%data)
   nullify(paw_setupin%kinetic_energy_differences%data)
 end if
 if (associated( paw_setupin%ae_partial_wave)) then
   do ii=1,paw_setupin%valence_states%nval
     if(associated( paw_setupin%ae_partial_wave(ii)%data)) then
       ABI_DEALLOCATE(paw_setupin%ae_partial_wave(ii)%data)
       nullify(paw_setupin%ae_partial_wave(ii)%data)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(paw_setupin%ae_partial_wave)
   nullify(paw_setupin%ae_partial_wave)
 end if
 if (associated( paw_setupin%pseudo_partial_wave)) then
   do ii=1,paw_setupin%valence_states%nval
     if(associated( paw_setupin%pseudo_partial_wave(ii)%data)) then
       ABI_DEALLOCATE(paw_setupin%pseudo_partial_wave(ii)%data)
       nullify(paw_setupin%pseudo_partial_wave(ii)%data)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(paw_setupin%pseudo_partial_wave)
   nullify(paw_setupin%pseudo_partial_wave)
 end if
 if (associated( paw_setupin%projector_function)) then
   do ii=1,paw_setupin%valence_states%nval
     if(associated( paw_setupin%projector_function(ii)%data)) then
       ABI_DEALLOCATE(paw_setupin%projector_function(ii)%data)
       nullify(paw_setupin%projector_function(ii)%data)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(paw_setupin%projector_function)
   nullify(paw_setupin%projector_function)
 end if
 if(associated(paw_setupin%valence_states%state)) then
   ABI_DATATYPE_DEALLOCATE(paw_setupin%valence_states%state)
   nullify(paw_setupin%valence_states%state)
 end if
 if(associated( paw_setupin%radial_grid)) then
   ABI_DATATYPE_DEALLOCATE(paw_setupin%radial_grid)
   nullify(paw_setupin%radial_grid)
 end if


end subroutine destroy_paw_setup
!!***

!!****f* ABINIT/nullify_paw_setup
!! NAME
!! nullify_paw_setup
!!
!! FUNCTION
!!  nullify a paw_setup datastructure
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2013 ABINIT group (FJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  paw_setup<paw_setup_type>=Datatype gathering information on XML paw setup.
!!
!! PARENTS
!!      inpspheads
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_paw_setup(paw_setupin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_paw_setup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(inout) :: paw_setupin

!Local variables-------------------------------

! *********************************************************************

 nullify(paw_setupin%radial_grid)
 nullify(paw_setupin%valence_states%state)
 nullify(paw_setupin%shape_function%data)
 nullify(paw_setupin%pseudo_partial_wave)
 nullify(paw_setupin%ae_partial_wave)
 nullify(paw_setupin%projector_function)
 nullify(paw_setupin%ae_core_density%data)
 nullify(paw_setupin%pseudo_core_density%data)
 nullify(paw_setupin%pseudo_valence_density%data)
 nullify(paw_setupin%zero_potential%data)
 nullify(paw_setupin%ae_core_kinetic_energy_density%data)
 nullify(paw_setupin%pseudo_core_kinetic_energy_density%data)
 nullify(paw_setupin%kresse_joubert_local_ionic_potential%data)
 nullify(paw_setupin%blochl_local_ionic_potential%data)
 nullify(paw_setupin%kinetic_energy_differences%data)

end subroutine nullify_paw_setup
!!***

!!****f* ABINIT/copy_paw_setup
!! NAME
!! copy_paw_setup
!!
!! FUNCTION
!!  Copy a paw_setup datastructure into another
!!
!! INPUTS
!!  
!!  paw_setupin<paw_setup_type>=input paw_setup datastructure
!!
!! OUTPUT
!!  paw_setupout<paw_setup_type>=output paw_setup datastructure
!!
!! PARENTS
!!      inpspheads
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_paw_setup(paw_setupin,paw_setupout)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_paw_setup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 type(paw_setup_t),intent(in) :: paw_setupin
 type(paw_setup_t),intent(out) :: paw_setupout

!Local variables-------------------------------
 integer :: ii,sz1,sz2
! *********************************************************************

!scalars
 paw_setupout%version=paw_setupin%version
 paw_setupout%tread=paw_setupin%tread
 paw_setupout%ngrid=paw_setupin%ngrid
 paw_setupout%idgrid=paw_setupin%idgrid
 paw_setupout%atom%tread=paw_setupin%atom%tread
 paw_setupout%atom%symbol=paw_setupin%atom%symbol
 paw_setupout%atom%znucl=paw_setupin%atom%znucl
 paw_setupout%atom%zion=paw_setupin%atom%zion
 paw_setupout%atom%zval=paw_setupin%atom%zval
 paw_setupout%xc_functional%tread=paw_setupin%xc_functional%tread
 paw_setupout%xc_functional%functionaltype=paw_setupin%xc_functional%functionaltype
 paw_setupout%xc_functional%name=paw_setupin%xc_functional%name
 paw_setupout%generator%tread=paw_setupin%generator%tread
 paw_setupout%generator%gen=paw_setupin%generator%gen
 paw_setupout%generator%name=paw_setupin%generator%name
 paw_setupout%valence_states%tread=paw_setupin%valence_states%tread
 paw_setupout%valence_states%nval=paw_setupin%valence_states%nval
 paw_setupout%shape_function%tread=paw_setupin%shape_function%tread
 paw_setupout%shape_function%gtype=paw_setupin%shape_function%gtype
 paw_setupout%shape_function%grid=paw_setupin%shape_function%grid
 paw_setupout%shape_function%rc=paw_setupin%shape_function%rc
 paw_setupout%shape_function%lamb=paw_setupin%shape_function%lamb
 paw_setupout%ae_core_density%tread=paw_setupin%ae_core_density%tread
 paw_setupout%ae_core_density%grid=paw_setupin%ae_core_density%grid
 paw_setupout%ae_core_density%state=paw_setupin%ae_core_density%state
 paw_setupout%pseudo_core_density%tread=paw_setupin%pseudo_core_density%tread
 paw_setupout%pseudo_core_density%grid=paw_setupin%pseudo_core_density%grid
 paw_setupout%pseudo_core_density%state=paw_setupin%pseudo_core_density%state
 paw_setupout%pseudo_valence_density%tread=paw_setupin%pseudo_valence_density%tread
 paw_setupout%pseudo_valence_density%grid=paw_setupin%pseudo_valence_density%grid
 paw_setupout%pseudo_valence_density%state=paw_setupin%pseudo_valence_density%state
 paw_setupout%zero_potential%tread=paw_setupin%zero_potential%tread
 paw_setupout%zero_potential%grid=paw_setupin%zero_potential%grid
 paw_setupout%zero_potential%state=paw_setupin%zero_potential%state
 paw_setupout%ae_core_kinetic_energy_density%tread=&
&     paw_setupin%ae_core_kinetic_energy_density%tread
 paw_setupout%ae_core_kinetic_energy_density%grid=&
&     paw_setupin%ae_core_kinetic_energy_density%grid
 paw_setupout%ae_core_kinetic_energy_density%state=&
&     paw_setupin%ae_core_kinetic_energy_density%state
 paw_setupout%pseudo_core_kinetic_energy_density%tread=&
&     paw_setupin%pseudo_core_kinetic_energy_density%tread
 paw_setupout%pseudo_core_kinetic_energy_density%grid=&
&     paw_setupin%pseudo_core_kinetic_energy_density%grid
 paw_setupout%pseudo_core_kinetic_energy_density%state=&
&     paw_setupin%pseudo_core_kinetic_energy_density%state
 paw_setupout%kresse_joubert_local_ionic_potential%tread=&
&    paw_setupin%kresse_joubert_local_ionic_potential%tread
 paw_setupout%kresse_joubert_local_ionic_potential%grid=&
&    paw_setupin%kresse_joubert_local_ionic_potential%grid
 paw_setupout%kresse_joubert_local_ionic_potential%state=&
&    paw_setupin%kresse_joubert_local_ionic_potential%state
 paw_setupout%blochl_local_ionic_potential%tread=&
&    paw_setupin%blochl_local_ionic_potential%tread
 paw_setupout%blochl_local_ionic_potential%grid=&
&    paw_setupin%blochl_local_ionic_potential%grid
 paw_setupout%blochl_local_ionic_potential%state=&
&    paw_setupin%blochl_local_ionic_potential%state
 paw_setupout%kinetic_energy_differences%tread=paw_setupin%kinetic_energy_differences%tread
 paw_setupout%kinetic_energy_differences%grid=paw_setupin%kinetic_energy_differences%grid
 paw_setupout%kinetic_energy_differences%state=paw_setupin%kinetic_energy_differences%state

!pointers
 if (associated(paw_setupin%shape_function%data)) then
   sz1=size(paw_setupin%shape_function%data,1)
   sz2=size(paw_setupin%shape_function%data,2)
   ABI_ALLOCATE(paw_setupout%shape_function%data,(sz1,sz2))
   paw_setupout%shape_function%data=paw_setupin%shape_function%data
 else
   nullify(paw_setupout%shape_function%data)
 end if
 if (associated(paw_setupin%ae_core_density%data)) then
   sz1=size(paw_setupin%ae_core_density%data,1)
   ABI_ALLOCATE(paw_setupout%ae_core_density%data,(sz1))
   paw_setupout%ae_core_density%data=paw_setupin%ae_core_density%data
 else
   nullify(paw_setupout%ae_core_density%data)
 end if
 if (associated(paw_setupin%pseudo_core_density%data)) then
   sz1=size(paw_setupin%pseudo_core_density%data,1)
   ABI_ALLOCATE(paw_setupout%pseudo_core_density%data,(sz1))
   paw_setupout%pseudo_core_density%data=paw_setupin%pseudo_core_density%data
 else
   nullify(paw_setupout%pseudo_core_density%data)
 end if
 if (associated(paw_setupin%pseudo_valence_density%data)) then
   sz1=size(paw_setupin%pseudo_valence_density%data,1)
   ABI_ALLOCATE(paw_setupout%pseudo_valence_density%data,(sz1))
   paw_setupout%pseudo_valence_density%data=paw_setupin%pseudo_valence_density%data
 else
   nullify(paw_setupout%pseudo_valence_density%data)
 end if
 if (associated(paw_setupin%zero_potential%data)) then
   sz1=size(paw_setupin%zero_potential%data,1)
   ABI_ALLOCATE(paw_setupout%zero_potential%data,(sz1))
   paw_setupout%zero_potential%data=paw_setupin%zero_potential%data
 else
   nullify(paw_setupout%zero_potential%data)
 end if
 if (associated(paw_setupin%ae_core_kinetic_energy_density%data)) then
   sz1=size(paw_setupin%ae_core_kinetic_energy_density%data,1)
   ABI_ALLOCATE(paw_setupout%ae_core_kinetic_energy_density%data,(sz1))
   paw_setupout%ae_core_kinetic_energy_density%data=paw_setupin%ae_core_kinetic_energy_density%data
 else
   nullify(paw_setupout%ae_core_kinetic_energy_density%data)
 end if
 if (associated(paw_setupin%pseudo_core_kinetic_energy_density%data)) then
   sz1=size(paw_setupin%pseudo_core_kinetic_energy_density%data,1)
   ABI_ALLOCATE(paw_setupout%pseudo_core_kinetic_energy_density%data,(sz1))
   paw_setupout%pseudo_core_kinetic_energy_density%data=paw_setupin%pseudo_core_kinetic_energy_density%data
 else
   nullify(paw_setupout%pseudo_core_kinetic_energy_density%data)
 end if
 if (associated(paw_setupin%kresse_joubert_local_ionic_potential%data)) then
   sz1=size(paw_setupin%kresse_joubert_local_ionic_potential%data,1)
   ABI_ALLOCATE(paw_setupout%kresse_joubert_local_ionic_potential%data,(sz1))
   paw_setupout%kresse_joubert_local_ionic_potential%data=paw_setupin%kresse_joubert_local_ionic_potential%data
 else
   nullify(paw_setupout%kresse_joubert_local_ionic_potential%data)
 end if
 if (associated(paw_setupin%blochl_local_ionic_potential%data)) then
   sz1=size(paw_setupin%blochl_local_ionic_potential%data,1)
   ABI_ALLOCATE(paw_setupout%blochl_local_ionic_potential%data,(sz1))
   paw_setupout%blochl_local_ionic_potential%data=paw_setupin%blochl_local_ionic_potential%data
 else
   nullify(paw_setupout%blochl_local_ionic_potential%data)
 end if
 if (associated(paw_setupin%kinetic_energy_differences%data)) then
   sz1=size(paw_setupin%kinetic_energy_differences%data,1)
   ABI_ALLOCATE(paw_setupout%kinetic_energy_differences%data,(sz1))
   paw_setupout%kinetic_energy_differences%data=paw_setupin%kinetic_energy_differences%data
 else
   nullify(paw_setupout%kinetic_energy_differences%data)
 end if
 if(associated( paw_setupin%radial_grid)) then
   sz1=size(paw_setupin%radial_grid,1)
   ABI_DATATYPE_ALLOCATE(paw_setupout%radial_grid,(sz1))
   paw_setupout%radial_grid=paw_setupin%radial_grid
 else
   nullify(paw_setupout%radial_grid)
 end if
 if(associated(paw_setupin%valence_states%state)) then
   sz1=size(paw_setupin%valence_states%state,1)
   ABI_DATATYPE_ALLOCATE(paw_setupout%valence_states%state,(sz1))
   paw_setupout%valence_states%state=paw_setupin%valence_states%state
 else
   nullify(paw_setupout%valence_states%state)
 end if

 if (associated( paw_setupin%ae_partial_wave)) then
   sz1=size(paw_setupin%ae_partial_wave,1)
   ABI_DATATYPE_ALLOCATE(paw_setupout%ae_partial_wave,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%ae_partial_wave(ii)%tread=paw_setupin%ae_partial_wave(ii)%tread
     paw_setupout%ae_partial_wave(ii)%grid=paw_setupin%ae_partial_wave(ii)%grid
     paw_setupout%ae_partial_wave(ii)%state=paw_setupin%ae_partial_wave(ii)%state
     if(associated( paw_setupin%ae_partial_wave(ii)%data)) then
       sz1=size(paw_setupin%ae_partial_wave(ii)%data,1)
       ABI_ALLOCATE(paw_setupout%ae_partial_wave(ii)%data,(sz1))
       paw_setupout%ae_partial_wave(ii)%data=paw_setupin%ae_partial_wave(ii)%data
     else
       nullify(paw_setupout%ae_partial_wave(ii)%data)
     end if
   end do
 else
   nullify(paw_setupout%ae_partial_wave)
 end if 
 if (associated( paw_setupin%pseudo_partial_wave)) then
   sz1=size(paw_setupin%pseudo_partial_wave,1)
   ABI_DATATYPE_ALLOCATE(paw_setupout%pseudo_partial_wave,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%pseudo_partial_wave(ii)%tread=paw_setupin%pseudo_partial_wave(ii)%tread
     paw_setupout%pseudo_partial_wave(ii)%grid=paw_setupin%pseudo_partial_wave(ii)%grid
     paw_setupout%pseudo_partial_wave(ii)%state=paw_setupin%pseudo_partial_wave(ii)%state
     if(associated( paw_setupin%pseudo_partial_wave(ii)%data)) then
       sz1=size(paw_setupin%pseudo_partial_wave(ii)%data,1)
       ABI_ALLOCATE(paw_setupout%pseudo_partial_wave(ii)%data,(sz1))
       paw_setupout%pseudo_partial_wave(ii)%data=paw_setupin%pseudo_partial_wave(ii)%data
     else
       nullify(paw_setupout%pseudo_partial_wave(ii)%data)
     end if
   end do
 else
   nullify(paw_setupout%pseudo_partial_wave)
 end if 
  if (associated( paw_setupin%projector_function)) then
   sz1=size(paw_setupin%projector_function,1)
   ABI_DATATYPE_ALLOCATE(paw_setupout%projector_function,(sz1))
   do ii=1,paw_setupin%valence_states%nval
     paw_setupout%projector_function(ii)%tread=paw_setupin%projector_function(ii)%tread
     paw_setupout%projector_function(ii)%grid=paw_setupin%projector_function(ii)%grid
     paw_setupout%projector_function(ii)%state=paw_setupin%projector_function(ii)%state
     if(associated( paw_setupin%projector_function(ii)%data)) then
       sz1=size(paw_setupin%projector_function(ii)%data,1)
       ABI_ALLOCATE(paw_setupout%projector_function(ii)%data,(sz1))
       paw_setupout%projector_function(ii)%data=paw_setupin%projector_function(ii)%data
     else
       nullify(paw_setupout%projector_function(ii)%data)
     end if
   end do
 else
   nullify(paw_setupout%projector_function)
 end if 

end subroutine copy_paw_setup
!!***

!{\src2tex{textfont=tt}}
!!****f* m_xml_pawpseudo/paw_rdfromline
!! NAME
!! paw_rdfromline
!!
!! FUNCTION
!! Read the value of a keyword from a XML line
!!
!! COPYRIGHT
!! Copyright (C) 1998-2013 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  keyword= keyword which value has to be read
!!  line= string from which the data are read (line from a XML)
!!
!! OUTPUT
!!  ierr= error code
!!  output= (string) value of the keyword
!!
!! PARENTS
!!      rdpawpsxml
!!
!! CHILDREN
!!
!! SOURCE

 subroutine paw_rdfromline(keyword,line,output,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_rdfromline'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
  character*(*), intent(in) :: keyword,line
  character*(*), intent(out) ::output
  integer, intent(out) :: ierr

#ifndef HAVE_LIBPAW_ABINIT
!Local variables ---------------------------------------
  character*(fnlen) :: temp
  integer :: pos,pos2
#endif

! *********************************************************************

#if defined HAVE_LIBPAW_ABINIT
 call rdfromline(keyword,line,output,ierr)
#else

 ierr=1;output=""
 pos=index(line,trim(keyword))
 if (pos>0) then
   temp=line(pos+len_trim(keyword):len_trim(line))
   pos=index(temp,char(34))
   if (pos>0) then
     pos2=index(temp(pos+1:len_trim(temp)),char(34))
     if (pos2>0) then
       output=temp(pos+1:pos+pos2-1)
     end if
   end if
 end if
#endif

 end subroutine paw_rdfromline
!!***


!{\src2tex{textfont=tt}}
!!****f* m_xml_pawpseudo/rdpawpsxml
!! NAME
!! rdpawpsxml
!!
!! FUNCTION
!! Read the header of a PAW pseudopotential XML file generated by AtomPAW
!!
!! COPYRIGHT
!! Copyright (C) 1998-2013 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  filename= input file name (atomicdata XML)
!!  funit= input unit number
!!
!! OUTPUT
!!  paw_setup=pseudopotential data structure
!!
!! PARENTS
!!      inpspheads
!!
!! CHILDREN
!!      paw_rdfromline
!!
!! SOURCE

 subroutine rdpawpsxml(filename,paw_setup,funit)

 use m_xml_pawpseudo_types

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdpawpsxml'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
 integer, intent(in) :: funit
 character (len=fnlen),intent(in) :: filename
 type(paw_setup_t),intent(inout) :: paw_setup
!Local variables ---------------------------------------
 integer :: iaewf,ii,ipswf,iproj,ir,igrid,ival,ierr,ishpf,lmax,mesh_size
 logical :: endfile,found
 character(len=500) :: message
 character (len=5000) :: line,readline
 character (len=50000) :: strg
 character (len=30) ::strg1
 real(dp), allocatable :: shpf(:,:)
 type(state_t), pointer   :: valstate (:)
 type(radial_grid_t), pointer   :: grids (:)
! *************************************************************************

!Inits
 

!Open the atomicdata XML file for reading
 open(unit=funit,file=filename,form='formatted',status='old')


!Start a reading loop
 endfile=.false.
 found=.false.
 
 do while ((.not.endfile).and.(.not.found))
   read(funit,'(a)',err=10,end=10) readline
   line=adjustl(readline);goto 20
   10 line="";endfile=.true.
   20 continue

!  --Read VERSION
   if (line(1:10)=='<paw_setup') then
     paw_setup%tread=.true.
     igrid=0;ishpf=0
     ABI_DATATYPE_ALLOCATE(grids,(10))

     call paw_rdfromline(" version",line,strg,ierr)
     paw_setup%version=trim(strg)
     cycle
   end if
!  --Read TITLE, ATOMIC CHARGE AND CORE CHARGE
   if (line(1:12)=='<atom symbol') then

     paw_setup%atom%tread=.true.
     call paw_rdfromline("atom symbol",line,strg,ierr)
     paw_setup%atom%symbol=trim(strg)
     call paw_rdfromline(" Z",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%znucl
     else
       read(unit=strg,fmt=*) paw_setup%atom%znucl
     end if
     call paw_rdfromline(" core",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zion
     else
       read(unit=strg,fmt=*) paw_setup%atom%zion
     end if
     call paw_rdfromline(" valence",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) paw_setup%atom%zval
     else
       read(unit=strg,fmt=*) paw_setup%atom%zval
     end if
     cycle
   end if
!  --Read EXCHANGE-CORRELATION TYPE
   if (line(1:14)=='<xc_functional') then
     paw_setup%xc_functional%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%xc_functional%functionaltype = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%xc_functional%name= trim(strg)
     cycle
   end if
!  --Read GENERATOR
   if (line(1:10)=='<generator') then
     paw_setup%generator%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%generator%gen  = trim(strg)
     call paw_rdfromline(" name",line,strg,ierr)
     paw_setup%generator%name= trim(strg)
     cycle
   end if


!  --Read BASIS SIZE, ORBITALS, RC AND OCCUPATIONS/STATE IDs
   if (line(1:16)=='<valence_states>') then
     paw_setup%valence_states%tread=.true.
     ABI_DATATYPE_ALLOCATE(valstate,(50))
     ival=0
     lmax=0
     do while (line(1:17)/='</valence_states>')
       read(funit,'(a)') readline;line=adjustl(readline)
       if (line(1:6)=='<state') then
         ival=ival+1
         if (ival>50) then
           write(std_out,'(a)') "Error in rdpawps1xml: basis size too large (>50) !"
           close(funit);stop
         end if
         call paw_rdfromline(" n",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%nn=-1    
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%nn
           else
             read(unit=strg,fmt=*) valstate(ival)%nn
           end if
         end if
         call paw_rdfromline(" l",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ll
         else
           read(unit=strg,fmt=*) valstate(ival)%ll
         end if
         if(valstate(ival)%ll>lmax) lmax=valstate(ival)%ll
         call paw_rdfromline(" f",line,strg,ierr)
         if (strg == "" ) then 
           valstate(ival)%ff=-1.d0
         else
           if (len(trim(strg))<=30) then
             strg1=trim(strg)
             read(unit=strg1,fmt=*) valstate(ival)%ff
           else
             read(unit=strg,fmt=*) valstate(ival)%ff
           end if
         end if
         call paw_rdfromline(" rc",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%rc
         else
           read(unit=strg,fmt=*) valstate(ival)%rc
         end if
         call paw_rdfromline(" e",line,strg,ierr)
         if (len(trim(strg))<=30) then
           strg1=trim(strg)
           read(unit=strg1,fmt=*) valstate(ival)%ee
         else
           read(unit=strg,fmt=*) valstate(ival)%ee
         end if
         call paw_rdfromline(" id",line,strg,ierr)
         valstate(ival)%id = trim(strg)
       end if
     end do
     cycle
   end if

!  --Read MESH_STEP AND NUMBER OF POINTS
   if (line(1:12)=='<radial_grid')then
     igrid=igrid+1
     call paw_rdfromline(" eq",line,strg,ierr)
     grids(igrid)%eq = trim(strg)
     call paw_rdfromline(" a",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%aa=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%aa
       else
         read(unit=strg,fmt=*) grids(igrid)%aa
       end if
     end if
     call paw_rdfromline(" n",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%nn=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%nn
       else
         read(unit=strg,fmt=*) grids(igrid)%nn
       end if
     end if
     call paw_rdfromline(" d",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%dd=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%dd
       else
         read(unit=strg,fmt=*) grids(igrid)%dd
       end if
     end if
     call paw_rdfromline(" b",line,strg,ierr)
     if (strg == "" ) then
       grids(igrid)%bb=0.d0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) grids(igrid)%bb
       else
         read(unit=strg,fmt=*) grids(igrid)%bb
       end if
     end if
     call paw_rdfromline("istart",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%istart
     else
       read(unit=strg,fmt=*) grids(igrid)%istart
     end if
     call paw_rdfromline("iend",line,strg,ierr)
     if (len(trim(strg))<=30) then
       strg1=trim(strg)
       read(unit=strg1,fmt=*) grids(igrid)%iend
     else
       read(unit=strg,fmt=*) grids(igrid)%iend
     end if
     call paw_rdfromline(" id",line,strg,ierr)
     grids(igrid)%id = trim(strg)
     if(igrid>10)then
       write(std_out,'(a)')("igrid>10")
       close(funit);stop
     end if
     cycle
   end if

!  --Read SHAPE TYPE
   if (line(1:15)=='<shape_function') then
     paw_setup%shape_function%tread=.true.
     call paw_rdfromline(" type",line,strg,ierr)
     paw_setup%shape_function%gtype = trim(strg)
     call paw_rdfromline(" rc",line,strg,ierr)
     if (strg /= "" ) then
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%rc
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%rc
       end if
     end if
     call paw_rdfromline(" lamb",line,strg,ierr)
     if (strg == "" ) then
       paw_setup%shape_function%lamb=0
     else
       if (len(trim(strg))<=30) then
         strg1=trim(strg)
         read(unit=strg1,fmt=*) paw_setup%shape_function%lamb
       else
         read(unit=strg,fmt=*) paw_setup%shape_function%lamb
       end if
     end if
     found=paw_setup%shape_function%tread
     call paw_rdfromline("grid",line,strg,ierr)
     paw_setup%shape_function%grid=trim(strg)
     if (strg /= "" ) then
       paw_setup%shape_function%gtype ="num" 
       do ii=1,igrid
         if(trim(paw_setup%shape_function%grid)==trim(grids(ii)%id)) then
           mesh_size=grids(ii)%iend-grids(ii)%istart+1
           exit
         end if
       end do
       if(.not.allocated(shpf)) then
         ABI_ALLOCATE(shpf,(mesh_size,7))
       end if
       ishpf=ishpf+1
       read(funit,*) (shpf(ir,ishpf),ir=1,mesh_size)
       call paw_rdfromline(" l",line,strg,ierr)
       if (strg /= "" ) then
         found=.false.
         if(paw_setup%valence_states%tread) then
           if(ishpf==2*lmax+1) found=.true.
         else
           write(message,'(a,a,a)')"the grids and the states must be read before the shapefunction",ch10,&
&           "Action: Modify your XML PAW data file"
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
       end if
     end if
     cycle
   end if

!  End of reading loop
 end do
 
 if(igrid==0.or.ival==0) then
   write(message,'(a,a,a)')"the grids and the states must be read before the shapefunction",ch10,&
&   "Action: Modify your XML PAW data file"
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(ishpf>0)then
   ABI_ALLOCATE(paw_setup%shape_function%data,(mesh_size,ishpf))
   do ii=1,ishpf
     paw_setup%shape_function%data(:,ii)=shpf(:,ii)
   end do
   ABI_DEALLOCATE(shpf)
 end if

 if(ival>0)then
   ABI_DATATYPE_ALLOCATE(paw_setup%valence_states%state,(ival))
   paw_setup%valence_states%state(ival)%tread=.true.
   paw_setup%valence_states%nval=ival
   do ii=1,ival
     paw_setup%valence_states%state(ii)=valstate(ii)
   end do
 end if
 ABI_DATATYPE_DEALLOCATE(valstate)
 if(.not.associated(paw_setup%ae_partial_wave)) then
   ABI_DATATYPE_ALLOCATE(paw_setup%ae_partial_wave,(paw_setup%valence_states%nval))
 end if
 if(.not.associated(paw_setup%pseudo_partial_wave)) then
   ABI_DATATYPE_ALLOCATE(paw_setup%pseudo_partial_wave,(paw_setup%valence_states%nval))
 end if
 if(.not.associated(paw_setup%projector_function)) then
   ABI_DATATYPE_ALLOCATE(paw_setup%projector_function,(paw_setup%valence_states%nval))
 end if

 ABI_DATATYPE_ALLOCATE(paw_setup%radial_grid,(igrid))
 paw_setup%radial_grid(igrid)%tread=.true.
 paw_setup%ngrid=igrid
 do ii=1,igrid
   paw_setup%radial_grid(ii)=grids(ii)
 end do
 ABI_DATATYPE_DEALLOCATE(grids)

!Start a reading loop
 ipswf=0;iaewf=0;iproj=0
 endfile=.false.
 do while (.not.endfile)
   read(funit,'(a)',err=11,end=11) readline
   line=adjustl(readline);goto 21
   11 line="";endfile=.true.
   21 continue

!  --Read core density CORE_DENSITY
   if (line(1:16)=='<ae_core_density') then
     paw_setup%ae_core_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%ae_core_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%ae_core_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%ae_core_density%data,(mesh_size))
     read(funit,*) (paw_setup%ae_core_density%data(ir),ir=1,mesh_size)
     cycle
   end if


!  --Read pseudized core density CORETAIL_DENSITY
   if (line(1:20)=='<pseudo_core_density') then
     paw_setup%pseudo_core_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%pseudo_core_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_core_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%pseudo_core_density%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_core_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read pseudized valence density PSEUDO_VALENCE_DENSITY
   if (line(1:23)=='<pseudo_valence_density') then
     paw_setup%pseudo_valence_density%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%pseudo_valence_density%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_valence_density%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%pseudo_valence_density%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_valence_density%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Vbare potential VLOCFUN
   if (line(1:15)=='<zero_potential') then
     paw_setup%zero_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%zero_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%zero_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%zero_potential%data,(mesh_size))
     read(funit,*) (paw_setup%zero_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Vloc for Abinit potential VLOC_ION
   if (line(1:37)=='<kresse_joubert_local_ionic_potential') then
     paw_setup%kresse_joubert_local_ionic_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%kresse_joubert_local_ionic_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%kresse_joubert_local_ionic_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%kresse_joubert_local_ionic_potential%data,(mesh_size))
     read(funit,*) (paw_setup%kresse_joubert_local_ionic_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

   if (line(1:29)=='<blochl_local_ionic_potential') then
     paw_setup%blochl_local_ionic_potential%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%blochl_local_ionic_potential%grid=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%blochl_local_ionic_potential%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%blochl_local_ionic_potential%data,(mesh_size))
     read(funit,*) (paw_setup%blochl_local_ionic_potential%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read WAVE FUNCTIONS PHI
   if (line(1:16)=='<ae_partial_wave') then
     iaewf=iaewf+1
     paw_setup%ae_partial_wave(iaewf)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%ae_partial_wave(iaewf)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%ae_partial_wave(iaewf)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%ae_partial_wave(iaewf)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%ae_partial_wave(iaewf)%data,(mesh_size))
     read(funit,*) (paw_setup%ae_partial_wave(iaewf)%data(ir),ir=1,mesh_size)
     cycle
   end if


!  --Read PSEUDO WAVE FUNCTIONS TPHI
   if (line(1:20)=='<pseudo_partial_wave') then
     ipswf=ipswf+1
     paw_setup%pseudo_partial_wave(ipswf)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%idgrid = trim(strg)
     paw_setup%pseudo_partial_wave(ipswf)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%pseudo_partial_wave(ipswf)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%pseudo_partial_wave(ipswf)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%pseudo_partial_wave(ipswf)%data,(mesh_size))
     read(funit,*) (paw_setup%pseudo_partial_wave(ipswf)%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read PROJECTORS TPROJ
   if (line(1:19)=='<projector_function') then
     iproj=iproj+1
     paw_setup%projector_function(iproj)%tread=.true.
     call paw_rdfromline(" grid",line,strg,ierr)
     if (strg == "" ) strg = "unknown"
     paw_setup%projector_function(iproj)%grid=trim(strg)
     call paw_rdfromline(" state",line,strg,ierr)
     paw_setup%projector_function(iproj)%state=trim(strg)
     do ii=1,paw_setup%ngrid
       if(trim(paw_setup%projector_function(iproj)%grid)==trim(paw_setup%radial_grid(ii)%id)) then
         mesh_size=paw_setup%radial_grid(ii)%iend-paw_setup%radial_grid(ii)%istart+1
         exit
       end if
     end do
     ABI_ALLOCATE(paw_setup%projector_function(iproj)%data,(mesh_size))
     read(funit,*) (paw_setup%projector_function(iproj)%data(ir),ir=1,mesh_size)
     cycle
   end if

!  --Read Kinetic term KINETIC_ENERGY_MATRIX
   if (line(1:28)=='<kinetic_energy_differences>') then
     paw_setup%kinetic_energy_differences%tread=.true.
     mesh_size=paw_setup%valence_states%nval*paw_setup%valence_states%nval
     ABI_ALLOCATE(paw_setup%kinetic_energy_differences%data,(mesh_size))
     read(funit,*) (paw_setup%kinetic_energy_differences%data(ir),ir=1,mesh_size)
     cycle
   end if


!  End of reading loop
 end do


!Close the XML atomicdata file
 close(funit)

!Test flags: is anything OK ?
 found=paw_setup%atom%tread.and.paw_setup%valence_states%tread.and.&
& paw_setup%xc_functional%tread.and.paw_setup%shape_function%tread
 if (.not.paw_setup%atom%tread) &
& write(std_out,'(a,i2,a)') "Error in rdpawpsxml: ATOM SYMBOL not found !"
 if (.not.paw_setup%valence_states%tread) &
& write(std_out,'(a,i2,a)') "Error in rdpawpsxml: VALENCE STATES not found !"
 if (.not.paw_setup%xc_functional%tread) &
& write(std_out,'(a,i2,a)') "Error in rdpawpsxml: EXCHANGE/CORRELATION not found !"
 if (.not.paw_setup%shape_function%tread) &
& write(std_out,'(a,i2,a)') "Error in rdpawps1xml: SHAPE FUNCTION TYPE not found !"
 if (.not.found) stop


 end subroutine rdpawpsxml
!!***


end module m_xml_pawpseudo
!!***

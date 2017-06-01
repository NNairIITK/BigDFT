!> @file
!! Module implementing the connection algorithm(s)
!!
!! @author 
!!    Copyright (C) 2014 UNIBAS, Bastian Schaefer 
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Data structure for module_connect.
!! Is in own routine, because module_io also needs it,
!! and module_connect uses module_io.
module module_connect_object
    use module_base
    use module_interfaces
    implicit none
    
    private

    public :: connect_object
    public :: allocate_connect_object
    public :: deallocate_connect_object

    !connect_object is used to reduce stack size required by
    !connect subroutine
    type connect_object
        real(gp), allocatable :: saddle(:,:,:)
        real(gp), allocatable :: enersad(:)
        real(gp), allocatable :: fsad(:,:,:)
        real(gp), allocatable :: fpsad(:,:)
        real(gp), allocatable :: rotforce(:,:,:)
        real(gp), allocatable :: minmode(:,:,:)

        real(gp), allocatable :: leftmin(:,:,:)
        real(gp), allocatable :: enerleft(:)
        real(gp), allocatable :: fleft(:,:,:)
        real(gp), allocatable :: fpleft(:,:)

        real(gp), allocatable :: rightmin(:,:,:)
        real(gp), allocatable :: enerright(:)
        real(gp), allocatable :: fright(:,:,:)
        real(gp), allocatable :: fpright(:,:)

        real(gp), allocatable :: rxyz1(:,:)
        real(gp), allocatable :: rxyz2(:,:)
        real(gp), allocatable :: tsgforces(:,:)
        real(gp) :: tsgenergy

        !for non-recursive routine
        real(gp), allocatable :: todorxyz(:,:,:,:)
        real(gp), allocatable :: todofp(:,:,:)
        real(gp), allocatable :: todoenergy(:,:)
        integer :: ntodo
    end type

contains
!=====================================================================
subroutine allocate_connect_object(nat,nid,nsadmax,cobj)
    use dynamic_memory
    implicit none
    !parameters
    integer, intent(in) :: nat, nid, nsadmax
    type(connect_object), intent(inout) :: cobj

    cobj%saddle    = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='saddle')
    cobj%enersad   = f_malloc((/1.to.nsadmax/),id='enersad')
    cobj%fsad      = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fsad')
    cobj%fpsad     = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpsad')
    cobj%rotforce  = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='rotorce')
    cobj%minmode   = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='minmode')
    cobj%leftmin   = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='leftmin')
    cobj%enerleft  = f_malloc((/1.to.nsadmax/),id='enerleft')
    cobj%fleft     = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fleft')
    cobj%fpleft    = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpleft')
    cobj%rightmin  = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='rightmin')
    cobj%enerright = f_malloc((/1.to.nsadmax/),id='enerright')
    cobj%fright    = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fright')
    cobj%fpright   = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpright')
    cobj%rxyz1     = f_malloc((/1.to.3,1.to.nat/),id='rxyz1')
    cobj%rxyz2     = f_malloc((/1.to.3,1.to.nat/),id='rxyz2')
    cobj%tsgforces = f_malloc((/1.to.3,1.to.nat/),id='tsgforces')
    cobj%todorxyz  = f_malloc((/1.to.3,1.to.nat,1.to.2,1.to.nsadmax/),&
                     id='todorxyz')
    cobj%todofp    = f_malloc((/1.to.nid,1.to.2,1 .to. nsadmax/),&
                     id='todofp')
    cobj%todoenergy= f_malloc((/1.to.2,1.to.nsadmax/),id='todoenergy')
    cobj%ntodo     =0
end subroutine allocate_connect_object
!=====================================================================
subroutine deallocate_connect_object(cobj)
    use dynamic_memory
    implicit none
    !parameters
    type(connect_object), intent(inout) :: cobj

    call f_free(cobj%saddle)
    call f_free(cobj%enersad)
    call f_free(cobj%fsad)
    call f_free(cobj%fpsad)
    call f_free(cobj%rotforce)
    call f_free(cobj%minmode)
    call f_free(cobj%leftmin)
    call f_free(cobj%enerleft)
    call f_free(cobj%fleft)
    call f_free(cobj%fpleft)
    call f_free(cobj%rightmin)
    call f_free(cobj%enerright)
    call f_free(cobj%fright)
    call f_free(cobj%fpright)
    call f_free(cobj%rxyz1)
    call f_free(cobj%rxyz2)
    call f_free(cobj%tsgforces)
    call f_free(cobj%todorxyz)
    call f_free(cobj%todofp)
    call f_free(cobj%todoenergy)
end subroutine deallocate_connect_object
end module

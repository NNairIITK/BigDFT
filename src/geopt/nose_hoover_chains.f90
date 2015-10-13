!> @file
!!  Routines for MD
!! @author
!!    Nisanth Nair
!!    Copyright (C) 2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
MODULE Nose_Hoover_Chains
  use module_base
  implicit none

  type, public :: NHC_data
     LOGICAL          :: nhchain
     INTEGER          :: eta_size, nhnc, nmultint, nsuzuki
     REAL (KIND=8)    :: nosefrq, enose
     REAL (KIND=8), DIMENSION(:), pointer :: eta, veta, aeta, noseq, dtys
  end type NHC_data

CONTAINS

  pure subroutine nullify_NHC_data(nhc)
    implicit none
    type(NHC_data), intent(out) :: nhc
    nhc%nhchain=.false.
    nhc%eta_size=0
    nhc%nhnc=0
    nhc%nmultint=0
    nhc%nsuzuki=0
    nhc%nosefrq=0.d0
    nhc%enose=0.d0
    nullify(nhc%eta)
    nullify(nhc%veta)
    nullify(nhc%aeta)
    nullify(nhc%noseq)
    nullify(nhc%dtys)
  end subroutine nullify_NHC_data

  subroutine initialize_NHC_data(nhc,nhchain,nhnc,nsuzuki,nmultint,nosefrq)
    implicit none
    logical, intent(in)  :: nhchain !<Nose Hoover Chain thermostat 
    INTEGER, intent(in)  :: nhnc, nmultint, nsuzuki
    REAL (KIND=8), intent(in) :: nosefrq
    type(NHC_data), intent(out) :: nhc

    !perform checks
    if (nhchain.and.all(nsuzuki /= [3,5,7])) then
       call f_err_throw('Error in the nsuzuki parameter, must be 3,5 or 7, found'//nsuzuki,&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

    call nullify_NHC_data(NHC)
    nhc%nhchain=nhchain
    nhc%nhnc=nhnc
    nhc%nsuzuki=nsuzuki
    nhc%nmultint=nmultint
    nhc%nosefrq=nosefrq

  end subroutine initialize_NHC_data

  !
  !     -----------------------------------------------------
  !>     Function : Nose-Hoover-chains are initialized 
  !!     Note     : Units used for Nose-Hoover chain integration 
  !!                is a.u. 
  !!     -----------------------------------------------------         
  !!
  subroutine nose_init(natoms,ndof,dt,T0ions,amass,vxyz,nhc)
    implicit none 
    integer :: NATOMS,NDOF
    REAL (KIND=8)   :: dt, T0ions, amass(natoms),vxyz(3,natoms) 
    type(NHC_data), intent(inout) :: nhc
    !local variables
    integer         :: i, ii, inhc, iatoms, dof
    logical         :: to_restart
    real (kind=8)   :: gauss, dummy, dum(2) 
    real (kind=8)   :: kt, fkt, akin, SIGMA 
    real (kind=8)   :: wys1, wys2, wys3, wys4

    !
    !
    REAL (KIND=8), PARAMETER  :: au_to_k   = 315774.664550534774D0, & 
         cmi_to_au = 7.26D-7, &
         pi        = 4.d0*ATAN(1.D0),  &
         cubroot   = 1.D0/3.D0

    nhc%eta_size = nhc%nhnc

    !ETA_SIZE=ETA_SIZE

    nhc%noseq=f_malloc0_ptr(nhc%eta_size,id='noseq')
    nhc%eta=f_malloc0_ptr(nhc%eta_size,id='eta')
    nhc%veta=f_malloc0_ptr(nhc%eta_size,id='veta')
    nhc%aeta=f_malloc0_ptr(nhc%eta_size,id='aeta')
    nhc%dtys=f_malloc0_ptr(7,id='dtys')
    !
    dof = ndof

    !   ** To a.u **
    kt  = T0ions / au_to_k
    fkt = kt * REAL(dof,gp)
    
    !   ** Nose freq. in a.u **
    nhc%nosefrq = nhc%nosefrq * cmi_to_au
    
    !   ** Nose masses in a.u **
    nhc%noseq(1) = fkt / nhc%nosefrq**2    
    do inhc = 2, nhc%nhnc
          nhc%noseq(inhc) = kt / nhc%nosefrq**2  
    end do
    
    !  ** Restart Nose Hoover Chain **
    !  FIXME

    ! nhc%veta(1:nhc%nhnc)=0.d0
    !   ** Assign initial nose velocities **
    do inhc = 1, nhc%nhnc, 2
       sigma=dsqrt(KT/nhc%noseq(inhc)) 
       CALL random_number(dum(1))
       CALL random_number(dum(2))
       nhc%veta(inhc) =sigma*dsqrt(-2.D0*dlog(dum(2)))*dcos(2.D0*pi*dum(1))
       if(INHC+1.le.nhc%nhnc)then 
          sigma=DSQRT(kt/nhc%noseq(inhc+1)) 
          nhc%veta(inhc+1)=sigma*dsqrt(-2.d0*dlog(dum(2)))*dsin(2.d0*pi*dum(1))
       end if
    end do
    CALL rescale_velocities(nhc%nhnc,1,nhc%nhnc,nhc%noseq,T0ions,nhc%veta,dum(1))

    nhc%eta(1:nhc%nhnc)=0.d0

    !Initializing factors for Yoshida-Suzuki integration 
    IF(nhc%nsuzuki == 3)THEN
       nhc%dtys(1) = 1.d0 / ( 2.d0 - 2.d0**cubroot )
       nhc%dtys(2) = 1.d0 - 2.d0 * nhc%dtys(1)
       nhc%dtys(3) = nhc%dtys(1)
       nhc%dtys(1) = nhc%dtys(1) * dt  / real(nhc%nmultint,gp) 
       nhc%dtys(2) = nhc%dtys(2) * dt  / real(nhc%nmultint,gp)
       nhc%dtys(3) = nhc%dtys(3) * dt  / real(nhc%nmultint,gp)
       nhc%dtys(4) = 0.d0
       nhc%dtys(5) = 0.d0
       nhc%dtys(6) = 0.d0
       nhc%dtys(7) = 0.d0
    ELSE IF(nhc%nsuzuki == 5)THEN
       nhc%dtys(1) = 1.d0 / ( 4.d0 - 4.d0**cubroot )
       nhc%dtys(3) = 1.d0 - 4.d0 * nhc%dtys(1)
       nhc%dtys(1) = nhc%dtys(1) * dt  / real(nhc%nmultint,gp)
       nhc%dtys(3) = nhc%dtys(3) * dt  / real(nhc%nmultint,gp)
       nhc%dtys(2) = nhc%dtys(1)
       nhc%dtys(4) = nhc%dtys(1)
       nhc%dtys(5) = nhc%dtys(1)
       nhc%dtys(6) = 0.d0
       nhc%dtys(7) = 0.d0
    ELSE IF(nhc%NSUZUKI == 7)THEN
       wys1    =  0.784513610477560D0
       wys2    =  0.235573213359357D0
       wys3    = -0.117767998417887D1
       wys4    =  1.d0 - 2.d0 * (wys1 + wys2 + wys3 )
       nhc%dtys(1) =  wys3 * dt  / real(nhc%nmultint,gp)
       nhc%dtys(2) =  wys2 * dt  / real(nhc%nmultint,gp)       
       nhc%dtys(3) =  wys1 * dt  / real(nhc%nmultint,gp)
       nhc%dtys(4) =  wys4 * dt  / real(nhc%nmultint,gp)
       nhc%dtys(5) =  wys1 * dt  / real(nhc%nmultint,gp)
       nhc%dtys(6) =  wys2 * dt  / real(nhc%nmultint,gp)
       nhc%dtys(7) =  wys3 * dt  / real(nhc%nmultint,gp)
    END IF
    
    RETURN
  end subroutine nose_init
  !
  subroutine nose_evolve(natoms,ndof,T0ions,amass,vxyz,nhc)
    implicit none 
    !   ------------------------------------------------------------------
    !   Function : Integration of Nose Hoover chain equations of motion
    !   ------------------------------------------------------------------
    integer :: natoms,ndof
    real (kind=8) :: T0ions, amass(*), vxyz(3,*)
    type(NHC_data), intent(inout) :: nhc
    !
    integer       ::  i, j, k, m, ii, inhc, iatoms, iii, dof
    real (kind=8) ::  wfun
    real (kind=8) ::  akin, efr, dtys2, dtys4, dtys8, scal, fkt, kt
    !
    real (kind=8), parameter  :: amu_to_au = 1822.888485D0, &
         au_to_k   = 315774.664550534774D0, & 
         cmi_to_au = 7.26D-7

    DOF = ndof
    !
    kt  = T0ions / au_to_k
    fkt = kt * real(dof,gp)

    scal = 1.d0

    !   ** Calculate kinetic energy **
    akin = 0.d0
    do i = 1, natoms
       akin = akin + amass(i)*&
            (vxyz(1,i)**2 +   &
             vxyz(2,i)**2 +   &
             vxyz(3,i)**2 )
    end do


    do k = 1, nhc%nmultint
       DO j = 1, nhc%nsuzuki

          dtys2 = nhc%dtys(j) / 2.d0
          dtys4 = nhc%dtys(j) / 4.d0
          dtys8 = nhc%dtys(j) / 8.d0

          nhc%aeta(nhc%nhnc) = (nhc%noseq(nhc%nhnc-1)*nhc%veta(nhc%nhnc-1)**2 & 
                                        - kt)/ nhc%noseq(nhc%nhnc)
          nhc%veta(nhc%nhnc) = nhc%veta(nhc%nhnc) +  dtys4 * nhc%aeta(nhc%nhnc)

          do ii = 1, nhc%nhnc - 2
             inhc = nhc%nhnc - ii
             efr  = dexp(-dtys8 * nhc%veta(inhc+1))
             nhc%veta(inhc) = nhc%veta(inhc)*efr 
             nhc%aeta(inhc) = (nhc%noseq(inhc-1)*nhc%veta(inhc-1)**2 - kt)/nhc%noseq(inhc)
             nhc%veta(inhc) = nhc%veta(inhc) + dtys4*nhc%aeta(inhc)
             nhc%veta(inhc) = nhc%veta(inhc)*efr
          end do
          efr = dexp(-dtys8*nhc%veta(2))
          nhc%veta(1) = nhc%veta(1)*efr 
          nhc%aeta(1) = (akin-fkt)/nhc%noseq(1)
          nhc%veta(1) = nhc%veta(1) + dtys4*nhc%aeta(1) 
          nhc%veta(1) = nhc%veta(1)*efr

          efr=dexp(-dtys2*nhc%veta(1))
          scal = scal * efr
          akin = akin * scal * scal  

          do inhc = 1, nhc%nhnc
             nhc%eta(inhc) = nhc%eta(inhc) + dtys2 * nhc%veta(inhc)
          end do

          efr = dexp(-dtys8 * nhc%veta(2) )
          nhc%veta(1) = nhc%veta(1)*efr 
          nhc%aeta(1) = (akin - fkt)/nhc%noseq(1) !TODO Martyna's paper has a typo here (index of Q)
          nhc%veta(1) = nhc%veta(1) + dtys4*nhc%aeta(1)
          nhc%veta(1) = nhc%veta(1)*efr 

          do inhc = 2, nhc%nhnc - 1
            efr = dexp(-dtys8*nhc%veta(inhc+1))
            nhc%veta(inhc) = nhc%veta(inhc)*efr
            nhc%aeta(inhc) = (nhc%noseq(inhc-1)*nhc%veta(inhc-1)**2 &
                              - kt)/nhc%noseq(inhc)
            nhc%veta(inhc) = nhc%veta(inhc) + dtys4*nhc%aeta(inhc)
            nhc%veta(inhc) = nhc%veta(inhc)*efr 
          end do

           nhc%aeta(nhc%nhnc) = (nhc%noseq(nhc%nhnc-1)*nhc%veta(nhc%nhnc-1)**2 &
                             - kt)/ nhc%noseq(nhc%nhnc)
           nhc%veta(nhc%nhnc) = nhc%veta(nhc%nhnc) + dtys4 * nhc%aeta(nhc%nhnc)

          end do
       end do
       do i = 1, natoms
          vxyz(1:3,i) = vxyz(1:3,i) * scal
       end do
  end subroutine nose_evolve

  subroutine nose_energy(natoms,ndof,T0ions,nhc)
    implicit none 
    !     Function : To calculate the energy of extended system variables
    integer :: natoms, ndof
    real (kind=8) :: T0ions
    type(NHC_data), intent(inout) :: nhc
    !
    real (kind=8), parameter  :: au_to_k   = 315774.664550534774D0

    !
    integer         :: inhc, iatoms, ii, dof
    real (kind=8)            :: kt, fkt

    dof = ndof
    kt  = T0ions / au_to_k 
    fkt = kt * real(dof,gp)
    !
    nhc%enose = 0.D0
    !
    !     ** Kinetic erergy + potential energy of the thermostat **
    !
    do inhc = 1, nhc%nhnc
       nhc%enose = nhc%enose + 0.5d0 * nhc%noseq(inhc) * nhc%veta(inhc) * nhc%veta(inhc)
       if(inhc == 1) then
          nhc%enose = nhc%enose + fkt * nhc%eta(1)
       else 
          nhc%enose = nhc%enose + kt  * nhc%eta(inhc)
       end if
    end do
    !
  end subroutine nose_energy
  !

end module Nose_Hoover_Chains

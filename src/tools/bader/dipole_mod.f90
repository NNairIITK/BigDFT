!> @file
!! Module for the Bader charge density analysis program
!! @author
!!    Copyright 2009 Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!!    Bader is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!    A copy of the GNU General Public License is available at
!!    http://www.gnu.org/licenses/
!!
!!    Written by Ali Sadeghi 2011
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! @warning
!!    From the forum at http://theory.cm.utexas.edu/forum/viewtopic.php?f=1&t=590&p=1921&hilit=dipole#p1921:
!!    We decided not to keep this up because the bader volumes are not neutral. 
!!    So there is a dipole due to the charge distribution in any one Bader volume, 
!!    but there will also be dipoles due to charge transfer between the volumes. 
!!    It was hard to make sense of the dipoles even for a water molecule.
!!    So the integration is easy to do and could be enabled again, 
!!    but the interpretation is difficult, and we decided not to include it.
!!
!!    A. Sadeghi: Knowing how polarized is an atom in the molecule can be instructive at least qulitively.
!!    In this code I include the contibution from charge transfers to avoid the contradiction mentioned in the forum.


!> Module for calculating dipole moment of individual Bader atoms and whole system
MODULE dipole_mod
  USE kind_mod
  USE matrix_mod
  USE options_mod
  USE ions_mod
  USE charge_mod
  USE io_mod
  USE chgcar_mod
  USE bader_mod
 
  IMPLICIT NONE

! Public, allocatable variables
  TYPE dipole_obj
    REAL(q2),DIMENSION(3) :: tot_elec, tot_core
      ! dipole moment of whole system caused by electronic charge density and ionic cores, respectivly. Since the net charge is zero it is does not dependend on the coordinate system
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: ion_polar, ion_chgtrans  
      ! dipole moment of individual Bader atoms from electronic charge polarization and charge transfer between bader volumes, respectivly. 
      !(see R. W. Bader et al, International Journal of Quantum Chemistry, Vol. 85, 592–607 (2001) )
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ion_netchg  
  END TYPE


  PRIVATE
  PUBLIC :: dipole_obj, dipole_cal, dipole_output


  CONTAINS


!> Calculate the dipole moments of each Bader volume with respect to the ions center
!! as well as of the the whole system (which is origin-independent)
  SUBROUTINE dipole_cal(bdr,ions,chg,dpl,opts)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    TYPE(dipole_obj) :: dpl
    TYPE(options_obj) :: opts

    REAL(q2), DIMENSION(3) ::rxyz,rxyz0, shift
    INTEGER :: n1,n2,n3, atom
    REAL(q2) :: dipoleunits

    WRITE(*,'(/,2x,A)') 'CALCULATING DIPOLE MOMENTS'

    ALLOCATE(dpl%ion_polar(ions%nions,3) , dpl%ion_chgtrans (ions%nions,3), dpl%ion_netchg(ions%nions))

!!   From http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/cubeplugin_8C-source.html
!!   As of VMD version 1.8.3, volumetric data points are 
!!   expected to represent the center of a grid box. cube format 
!!   volumetric data represents the value at the edges of the 
!!   grid boxes, so we need to shift the internal origin by half
!!   a grid box diagonal to have the data at the correct position.

    shift=1._q2 
    dipoleunits=1._q2
    if( opts%in_opt /= opts%in_cube) shift=.5_q2
    if( opts%in_opt /= opts%in_cube) dipoleunits=1._q2/.52917720859_q2
    if( opts%in_opt /= opts%in_cube) print*, & 
      "WARNING, only in cube format, the core charges is presented. Unreliable core contribution!"
    dpl%tot_elec=0_q2
    dpl%tot_core=0_q2
    dpl%ion_polar=0_q2
    dpl%ion_chgtrans=0_q2
    DO n1=1,chg%npts(1)
      DO n2=1,chg%npts(2)
        DO n3=1,chg%npts(3)
          rxyz(1)=(n1-shift(1))*ions%lattice(1,1)/chg%npts(1) + chg%org_car(1)
          rxyz(2)=(n2-shift(2))*ions%lattice(2,2)/chg%npts(2) + chg%org_car(2)
          rxyz(3)=(n3-shift(3))*ions%lattice(3,3)/chg%npts(3) + chg%org_car(3)
          !r123(1)=real (n1,q2)
          !r123(2)=real (n2,q2)
          !r123(3)=real (n3,q2)
          !rxyz(:) = lat2car(chg, r123(:))
          dpl%tot_elec(:)=dpl%tot_elec(:)- chg%rho(n1,n2,n3)*rxyz(:) 
          ! sum up dipole momet of each bader atom with respect to the position of the assigned atom
          ! by defintion dipole_ion = integral rho(r')*(r'-r_ion) d3r' + Q_ion*0
          IF (bdr%volnum(n1,n2,n3) /= bdr%nvols+1) THEN
              atom = bdr%nnion(bdr%volnum(n1,n2,n3))  ! get nearest neighboring atom to this voxel
              rxyz0= (rxyz(:)-ions%r_car(atom,:))
              CALL dpbc_car(ions,rxyz0)
              dpl%ion_polar(atom,:) =  dpl%ion_polar(atom,:) - chg%rho(n1,n2,n3)*rxyz0(:)
          ENDIF
        END DO
      END DO
    END DO
   dpl%tot_elec(:)=dpl%tot_elec(:)/chg%nrho *dipoleunits 
   dpl%ion_polar=dpl%ion_polar/chg%nrho     *dipoleunits 

! total dipole moment of ionic cores associated to whole system with respect to origin
     DO atom = 1,ions%nions  ! by defintion p=sum Q_ion*(r_ion-r_o)
       dpl%tot_core(:)=dpl%tot_core(:) + ions%ion_chg(atom)*ions%r_car(atom,:)
       dpl%ion_netchg(atom) = ions%ion_chg(atom)-bdr%ionchg(atom)
       dpl%ion_chgtrans(atom,:) =  dpl%ion_chgtrans(atom,:) + dpl%ion_netchg(atom)*ions%r_car(atom,:)
     ENDDO

    dpl%tot_elec=     dpl%tot_elec      *dipoleunits  
    dpl%tot_core=     dpl%tot_core      *dipoleunits 
    dpl%ion_polar=    dpl%ion_polar     *dipoleunits 
    dpl%ion_chgtrans= dpl%ion_chgtrans  *dipoleunits 
    RETURN
  END SUBROUTINE dipole_cal


!> Write out a summary of the dipole moment calulations.
!! dipole.dat : Stores the main output to the screen.
!! Use yaml format (TD)
 subroutine dipole_output(bdr,ions,dpl)
   use yaml_strings
    use yaml_output

    implicit none

    !Arguments
    type(bader_obj) :: bdr
    type(ions_obj) :: ions
    type(dipole_obj) :: dpl
    !Local variables
    integer, parameter :: iunit=400
    character(len=*), parameter :: format_line = "(100('-'))"
    character(len=*), parameter :: dfile = 'dipole.yaml'
    real(q2), parameter :: dipoleunits=1._q2  !< atomic units 
    real(q2), dimension(3) :: tmp1, tmp2
    real(q2) :: stmp,stmp2
    integer :: i
 
    !WRITE(*,'(A44,/)') 'WRITING BADER ATOMIC DIPOLES TO dipole.dat'
    !call yaml_map('Writing Bader atomic dipoles',dfile)

    !OPEN(UNIT=iunit,FILE='dipole.dat',STATUS='replace',ACTION='write')
    call yaml_set_stream(unit=iunit,filename=dfile)
    call yaml_sequence_open('Atoms coordinates')
    do i=1,ions%nions
       call yaml_sequence(advance='no')
       call yaml_comment('Atoms' // trim(yaml_toa(i)))
       call yaml_mapping_open(flow=.true.)
       call yaml_map('Coord.',ions%r_car(:,i))
       call yaml_map('Charge core',ions%ion_chg(i))
       call yaml_map('Charge elec.',-bdr%ionchg(i))
       call yaml_map('Charge net',dpl%ion_netchg(i))
       call yaml_mapping_close(unit=iunit)
    end do
    call yaml_sequence_close(unit=iunit)

    !WRITE(iunit,format_line) 
    !WRITE(iunit,'(a)') "Atoms coordinates: " 
    !WRITE(iunit,'(A)')  & 
    !'atom#    coordinates:  X           Y           Z           CHARGE:  core      electronic     net'
    !WRITE(iunit,format_line) 
    !DO i=1,ions%nions
    !  WRITE(iunit,'(I4,10x,3F12.4,12x,SP,3F12.5)') i,ions%r_car(i,:), & 
    !                 ions%ion_chg(i),-bdr%ionchg(i), dpl%ion_netchg(i) 
    !END DO
    !WRITE(iunit,format_line) 
    !WRITE(iunit,'(a)') ''

    !WRITE(iunit,format_line) 
    !WRITE(iunit,'(a)') & 
    ! "Atomic polarization dipole-moments with respect to the corresponding nuclei positions [e.a0]" 
    ! WRITE(iunit,'(3A)')  & 
    ! 'atom#         Intra-atomic:                Px          Py          Pz          |P|'   !(1 D=0.3934 e.a0)
    ! WRITE(iunit,format_line) 

    call yaml_comment('Atomic polarization dipole-moments with respect to the corresponding nuclei positions [e.a0]')
    call yaml_sequence_open('Atomic polarization dipole_moments')
    do i=1,ions%nions
       tmp1= dpl%ion_polar(i,:)*dipoleunits
       stmp = sqrt(dot_product(tmp1,tmp1))
       call yaml_sequence(advance='no')
       call yaml_comment('Atoms' // trim(yaml_toa(i)))
       call yaml_mapping_open(flow=.true.)
       call yaml_map('P',tmp1)
       call yaml_map('Norm P',sqrt(dot_product(tmp1,tmp1)))
       call yaml_mapping_close(unit=iunit)
       !WRITE(iunit,'(I3,33x,3F12.6,1x,F13.6,6x,F12.6,5x,F12.5)') & 
       !&   i , tmp1(:), sqrt(DOT_PRODUCT(tmp1(:),tmp1(:)))
    end do
    call yaml_sequence_close(unit=iunit)

    !WRITE(iunit,format_line) 
    tmp1(1)= sum(dpl%ion_polar(:,1))*dipoleunits
    tmp1(2)= sum(dpl%ion_polar(:,2))*dipoleunits
    tmp1(3)= sum(dpl%ion_polar(:,3))*dipoleunits
     
    tmp2(1)= sum(dpl%ion_chgtrans (:,1))*dipoleunits
    tmp2(2)= sum(dpl%ion_chgtrans (:,2))*dipoleunits
    tmp2(3)= sum(dpl%ion_chgtrans (:,3))*dipoleunits
    stmp2 = sqrt(dot_product(tmp2(:),tmp2(:)))
    stmp  = sqrt(dot_product(tmp1+tmp2,tmp1+tmp2))
!    WRITE(iunit,'(A33,3x,3F12.6,1x,F12.6,"   | ",F12.6, 5x F12.5 )') 'Summation:   ' , & 
!    WRITE(iunit,'(A33,3x,3F12.6,1x,F12.6 )') 'Summation:   ' , & 
!                     tmp1(:),  sqrt(DOT_PRODUCT(tmp1(:),tmp1(:)))   
!    WRITE(iunit,format_line) 

    !WRITE(iunit,'(A33,3x,3F12.6,1x,F13.6 )') 'Charge-transfer contribution:' , & 
    !                 tmp2(:), sqrt(DOT_PRODUCT(tmp2(:),tmp2(:)))   
    !WRITE(iunit,format_line) 
    !WRITE(iunit,'(A33,3x,3F12.6,1x,F13.6 )') 'Total dipole moment:' , & 
    !                 tmp1+tmp2,sqrt(DOT_PRODUCT(tmp1+tmp2,tmp1+tmp2))     

    call yaml_map('Charge-transfer contributions',(/ tmp2(1),tmp2(2),tmp2(3),stmp2 /))
    stmp = sqrt(dot_product(tmp1+tmp2,tmp1+tmp2))
    call yaml_map('Total dipole moment', (/ tmp1(1)+tmp2(1), tmp1(2)+tmp2(2), tmp1(3)+tmp2(3), stmp /)) 

    !WRITE(iunit,format_line) 
    tmp1= dpl%tot_core(1:3)*dipoleunits 
    tmp2= dpl%tot_elec(1:3)*dipoleunits 
    stmp = sqrt(dot_product(tmp1+tmp2,tmp1+tmp2))
    !WRITE(iunit,'(a)') & 
    !"Dipole-moment of the whole system (not decomposed to atoms) with respect to arbitrary origin:" 
    call yaml_comment('Dipole-moment of the whole system (not decomposed to atoms) with respect to arbitrary origin')
 !!   WRITE(iunit,'(A33,3x,3F12.6,1x,F13.6 )') 'Elec. dipole moment:' , & 
 !!                    tmp2     ,sqrt(DOT_PRODUCT(tmp2     ,tmp2     ))     
 !!   WRITE(iunit,'(A33,3x,3F12.6,1x,F13.6 )') 'Cores dipole moment:' , & 
 !!                         tmp1,sqrt(DOT_PRODUCT(     tmp1,     tmp1))     
    !WRITE(iunit,'(A33,3x,3F12.6,1x,F13.6 )') 'Total dipole moment:' , & 
    !                 tmp1+tmp2,sqrt(DOT_PRODUCT(tmp1+tmp2,tmp1+tmp2))     
    !WRITE(iunit,format_line) 

    call yaml_map('Total dipole moment of the whole system', (/ tmp1(1)+tmp2(1), tmp1(2)+tmp2(2), tmp1(3)+tmp2(3), stmp /))
    call yaml_close_stream(unit=iunit)
    !CLOSE(UNIT=iunit)

 END SUBROUTINE dipole_output


END MODULE dipole_mod

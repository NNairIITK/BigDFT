!> @file
!!  Routines to do Van der Waals correction
!! @author
!! Written by Quintin Hill in 2007/8 with assistance from
!! Chris-Kriton Skylaris.\n
!! Forces added July 2008 by Quintin Hill\n
!! Modified for BigDFT in March/April 2009 by Quintin Hill\n
!! Modified for BigDFT (DFT-D2) in April 2012 by vama\n
!! Modified for BigDFT (DFT-D3) in December 2012 by vama\n
!!
!! Copyright (C)  2007-2009  Quintin Hill.
!! This file is distributed under the terms of the
!! GNU General Public License either version 2 of the License, or
!! (at your option) any later version, see ~/COPYING file,
!! http://www.gnu.org/licenses/gpl-2.0.txt (for GPL v2)
!! or http://www.gnu.org/copyleft/gpl.txt (for latest version).
!!
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the Van der Waals empirical correction.
!! This module contains subroutines for calculating a Van der
!! Waals energy correction as a sum of damped London potentials.
module vdwcorrection

  use module_base, only: gp,f_malloc,f_malloc0,f_free,f_malloc_ptr,f_malloc0_ptr,f_free_ptr,assignment(=)

  implicit none

  private


  ! vama: DFT-D3
  integer, parameter :: MAX_ELEM = 94
  integer, parameter :: MAX_CN = 94
  integer, parameter :: ntot_values=254
  REAL(kind=GP), parameter :: K1=16.0_GP
  REAL(kind=GP), parameter :: K2=4.0_GP/3.0_GP
  REAL(kind=GP), parameter :: K3=-4.0_GP
  REAL(kind=GP), parameter :: ANGS2AU=1.889725989_GP
!   jnm6mol2au=(((10*angs2au)**6)/4186)/627.51
  REAL(kind=GP), parameter :: JNM6MOL2AU=17.336967408209503_GP
 
  ! qoh: Structure to store parameters in
  TYPE VDWPARAMETERS
     !qoh: C_6 coefficients
     REAL(kind=GP), DIMENSION(109) :: C6COEFF
     !qoh: R_0 values
     REAL(kind=GP), DIMENSION(109) :: RADZERO
     !qoh: damping coefficient
     REAL(kind=GP), dimension(3) :: DCOEFF
     !qoh: Effective number of electrons
     REAL(kind=GP), DIMENSION(109) :: NEFF
     !vama: Scale facfor for radii values
     REAL(kind=GP):: RADSCALE
     !vama: function parameters, to be transformed as pointers
     REAL(kind=GP):: S6, S8, SR6, SR8, ALPHA
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM) :: QATOM
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM) :: cov_table
!!$     integer, DIMENSION(MAX_ELEM)       :: MAXCN
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM,MAX_ELEM) :: R0AB
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM,MAX_ELEM,MAX_CN,MAX_CN,3) :: C6AB
     !to be allocated only for the Grimme's D3 correction calculation
     integer, dimension(ntot_values) :: ivalues
     real(gp), dimension(ntot_values) :: rvalues
     real(gp), dimension(:), pointer :: coeffs
     real(gp), dimension(:), pointer :: qatom
     real(gp), dimension(:), pointer :: cov_table
     integer, dimension(:), pointer       :: maxcn
     real(gp), dimension(:,:), pointer :: r0ab
     !real(gp), dimension(:,:,:,:,:), pointer :: c6ab
  END TYPE VDWPARAMETERS

  !> Van der Waals corrections.
  integer, parameter, public :: VDW_NONE = 0
  integer, parameter, public :: VDW_DAMP_ELSTNER = 1
  integer, parameter, public :: VDW_DAMP_WU_YANG_1 = 2
  integer, parameter, public :: VDW_DAMP_WU_YANG_2 = 3
  integer, parameter, public :: VDW_DAMP_GRIMME_D2 = 4
  integer, parameter, public :: VDW_DAMP_GRIMME_D3 = 5
  character(len = 34), dimension(6), parameter :: vdw_correction_names = &
       & (/ "none                              ",   &
       &    "Damp from Elstner                 ",   &    
       &    "Damp from Wu & Yang               ",   &
       &    "Damp from Wu & Yang, second method",   &
       &    "Damp from Grimme D2               ",   &
       &    "Damp from Grimme D3               " /)

  public :: vdwcorrection_initializeparams
  public :: vdwcorrection_freeparams
  public :: vdwcorrection_calculate_energy
  public :: vdwcorrection_calculate_forces
  public :: vdwcorrection_warnings

  ! qoh: Share array of Van der Waals parameters accross module

  type(VDWPARAMETERS) :: vdwparams ! Array of parameters
 

contains


  !< This subroutine populates the array of Van der Waals parameters
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  !! Some unoptimised parameters (e.g. for P) were taken from Halgren.
  !! Journal of the American Chemical Society 114(20), 7827–7843, 1992
  subroutine vdwcorrection_initializeparams(ixc, dispersion)
    implicit none

    integer, intent(in) :: ixc, dispersion
    integer             :: ii

    ! qoh: Initialise the vdwparams array
    vdwparams%c6coeff = 0.0_GP
    vdwparams%neff    = 0.0_GP
    ! qoh: Setting to 1 prevents division by zero
    vdwparams%radzero = 1.0_GP
    vdwparams%radscale = 1.0_GP

    ! qoh: Populate the vdwparams array  

    if (dispersion < 4) then

       select case (ixc)
   
       case(11)
   
          pbedamp: select case (dispersion)
   
          case(VDW_DAMP_ELSTNER) pbedamp
   
             vdwparams%dcoeff(1)=3.2607_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0000_GP

             vdwparams%c6coeff(1)=2.9239_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3561_GP
             vdwparams%c6coeff(7)=19.5089_GP
             vdwparams%c6coeff(8)=11.7697_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0368_GP
   
             vdwparams%radzero(1)=2.9635_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.3610_GP
             vdwparams%radzero(7)=3.5136_GP
             vdwparams%radzero(8)=3.7294_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=5.1033_GP
   
          case(VDW_DAMP_WU_YANG_1) pbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0000_GP
   
             vdwparams%c6coeff(1)=2.8450_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3200_GP
             vdwparams%c6coeff(7)=19.4800_GP
             vdwparams%c6coeff(8)=11.7600_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0232_GP
   
             vdwparams%radzero(1)=3.0000_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.8000_GP
             vdwparams%radzero(7)=3.8000_GP
             vdwparams%radzero(8)=3.8000_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=5.5998_GP
   
          case(VDW_DAMP_WU_YANG_2) pbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0025_GP
   
             vdwparams%c6coeff(1)=2.9865_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3784_GP
             vdwparams%c6coeff(7)=19.5223_GP
             vdwparams%c6coeff(8)=11.7733_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0388_GP
   
             vdwparams%radzero(1)=2.8996_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.0001_GP
             vdwparams%radzero(7)=3.2659_GP
             vdwparams%radzero(8)=3.6630_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=4.7948_GP
   
   
          end select pbedamp
   
       case(15)
   
          rpbedamp: select case (dispersion)
   
          case(VDW_DAMP_ELSTNER) rpbedamp
   
             vdwparams%dcoeff(1)=3.4967_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=3.2054_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.4608_GP
             vdwparams%c6coeff(7)=19.5765_GP
             vdwparams%c6coeff(8)=11.7926_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0465_GP
             vdwparams%radzero(1)=2.8062_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.0375_GP
             vdwparams%radzero(7)=3.2484_GP
             vdwparams%radzero(8)=3.6109_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=4.0826_GP
   
   
          case(VDW_DAMP_WU_YANG_1) rpbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.9308_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=3.2290_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.4707_GP
             vdwparams%c6coeff(7)=19.5922_GP
             vdwparams%c6coeff(8)=11.8007_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0388_GP
             vdwparams%radzero(1)=2.9406_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.4751_GP
             vdwparams%radzero(7)=3.6097_GP
             vdwparams%radzero(8)=3.7393_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=4.8007_GP
   
   
   
          case(VDW_DAMP_WU_YANG_2) rpbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0033_GP
             vdwparams%c6coeff(1)=3.0210_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3845_GP
             vdwparams%c6coeff(7)=19.5208_GP
             vdwparams%c6coeff(8)=11.7723_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0437_GP
             vdwparams%radzero(1)=2.8667_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=2.9914_GP
             vdwparams%radzero(7)=3.2931_GP
             vdwparams%radzero(8)=3.6758_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=3.9912_GP
   
          end select rpbedamp
   
       case(14) 
   
          revpbedamp: select case (dispersion)
   
          case(VDW_DAMP_ELSTNER) revpbedamp
             vdwparams%dcoeff(1)=3.4962_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=3.2167_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.4616_GP
             vdwparams%c6coeff(7)=19.5760_GP
             vdwparams%c6coeff(8)=11.7927_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0475_GP
             vdwparams%radzero(1)=2.7985_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.0398_GP
             vdwparams%radzero(7)=3.2560_GP
             vdwparams%radzero(8)=3.6122_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=3.9811_GP
   
   
          case(VDW_DAMP_WU_YANG_1) revpbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.8282_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=3.1160_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.4184_GP
             vdwparams%c6coeff(7)=19.5513_GP
             vdwparams%c6coeff(8)=11.7867_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0389_GP
             vdwparams%radzero(1)=2.9540_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.5630_GP
             vdwparams%radzero(7)=3.6696_GP
             vdwparams%radzero(8)=3.7581_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=4.7980_GP
   
   
          case(VDW_DAMP_WU_YANG_2) revpbedamp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0034_GP
             vdwparams%c6coeff(1)=3.0258_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3854_GP
             vdwparams%c6coeff(7)=19.5210_GP
             vdwparams%c6coeff(8)=11.7724_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0435_GP
             vdwparams%radzero(1)=2.8623_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=2.9923_GP
             vdwparams%radzero(7)=3.2952_GP
             vdwparams%radzero(8)=3.6750_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=3.9910_GP
   
   
          end select revpbedamp
   
       case(200)  ! qoh: FIXME: with proper number for PW91
   
          pw91damp: select case (dispersion)
   
          case(VDW_DAMP_ELSTNER) pw91damp
   
             vdwparams%dcoeff(1)=3.2106_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=2.8701_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3422_GP
             vdwparams%c6coeff(7)=19.5030_GP
             vdwparams%c6coeff(8)=11.7670_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0346_GP
             vdwparams%radzero(1)=3.0013_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.4423_GP
             vdwparams%radzero(7)=3.5445_GP
             vdwparams%radzero(8)=3.7444_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=5.4111_GP
   
          case(VDW_DAMP_WU_YANG_1) pw91damp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.3102_GP
             vdwparams%dcoeff(3)=23.0000_GP
             vdwparams%c6coeff(1)=2.5834_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.2743_GP
             vdwparams%c6coeff(7)=19.4633_GP
             vdwparams%c6coeff(8)=11.7472_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0233_GP
             vdwparams%radzero(1)=3.0759_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.9644_GP
             vdwparams%radzero(7)=3.8390_GP
             vdwparams%radzero(8)=3.8330_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=5.6250_GP
   
   
          case(VDW_DAMP_WU_YANG_2) pw91damp
   
             vdwparams%dcoeff(1)=3.0000_GP
             vdwparams%dcoeff(2)=3.5400_GP
             vdwparams%dcoeff(3)=23.0004_GP
   
             vdwparams%c6coeff(1)=2.9252_GP
             vdwparams%c6coeff(2)=0.0000_GP
             vdwparams%c6coeff(3)=0.0000_GP
             vdwparams%c6coeff(4)=0.0000_GP
             vdwparams%c6coeff(5)=0.0000_GP
             vdwparams%c6coeff(6)=27.3608_GP
             vdwparams%c6coeff(7)=19.5150_GP
             vdwparams%c6coeff(8)=11.7706_GP
             vdwparams%c6coeff(9)=0.0000_GP
             vdwparams%c6coeff(10)=0.0000_GP
             vdwparams%c6coeff(11)=0.0000_GP
             vdwparams%c6coeff(12)=0.0000_GP
             vdwparams%c6coeff(13)=0.0000_GP
             vdwparams%c6coeff(14)=0.0000_GP
             vdwparams%c6coeff(15)=190.5033_GP
             vdwparams%c6coeff(16)=161.0381_GP
   
             vdwparams%radzero(1)=2.9543_GP
             vdwparams%radzero(2)=3.0000_GP
             vdwparams%radzero(3)=3.8000_GP
             vdwparams%radzero(4)=3.8000_GP
             vdwparams%radzero(5)=3.8000_GP
             vdwparams%radzero(6)=3.0785_GP
             vdwparams%radzero(7)=3.3119_GP
             vdwparams%radzero(8)=3.6858_GP
             vdwparams%radzero(9)=3.8000_GP
             vdwparams%radzero(10)=3.8000_GP
             vdwparams%radzero(11)=4.8000_GP
             vdwparams%radzero(12)=4.8000_GP
             vdwparams%radzero(13)=4.8000_GP
             vdwparams%radzero(14)=4.8000_GP
             vdwparams%radzero(15)=4.8000_GP
             vdwparams%radzero(16)=4.9014_GP
   
          end select pw91damp
   
   
       case default  
   
          ! qoh: The damping coefficient for the three damping functions          
          vdwparams%dcoeff(1)=3.0000_GP
          vdwparams%dcoeff(2)=3.5400_GP
          vdwparams%dcoeff(3)=23.0000_GP
   
          ! qoh:  Values from Wu and Yang in hartree/bohr ^ 6
          vdwparams%c6coeff(1)=2.8450_GP
          vdwparams%c6coeff(2)=0.0000_GP
          vdwparams%c6coeff(3)=0.0000_GP
          vdwparams%c6coeff(4)=0.0000_GP
          vdwparams%c6coeff(5)=0.0000_GP
          vdwparams%c6coeff(6)=27.3200_GP
          vdwparams%c6coeff(7)=19.4800_GP
          vdwparams%c6coeff(8)=11.7600_GP
          ! qoh: C6 from Halgren
          vdwparams%c6coeff(16)=161.0388_GP
   
          ! qoh: vdw radii from Elstner in Angstrom 
          vdwparams%radzero(1)=3.0000_GP
          vdwparams%radzero(2)=3.0000_GP
          vdwparams%radzero(3)=3.8000_GP
          vdwparams%radzero(4)=3.8000_GP
          vdwparams%radzero(5)=3.8000_GP
          vdwparams%radzero(6)=3.8000_GP
          vdwparams%radzero(7)=3.8000_GP
          vdwparams%radzero(8)=3.8000_GP
          vdwparams%radzero(16)=4.8000_GP
   
       end select
   
       ! qoh: Unoptimised parameters from Halgren
       vdwparams%c6coeff(9)=6.2413_GP
       vdwparams%c6coeff(15)=190.5033_GP
       vdwparams%c6coeff(17)=103.5612_GP
       vdwparams%c6coeff(35)=201.8972_GP  
   
       vdwparams%radzero(9)=3.09_GP 
       vdwparams%radzero(15)=4.8000_GP
       vdwparams%radzero(17)=4.09_GP     
       vdwparams%radzero(35)=4.33_GP   
   
       ! qoh: Array containing the Neff from Wu and Yang 
       vdwparams%neff(1)=0.5300_GP
       vdwparams%neff(2)=0.0000_GP
       vdwparams%neff(3)=0.0000_GP
       vdwparams%neff(4)=0.0000_GP
       vdwparams%neff(5)=0.0000_GP
       vdwparams%neff(6)=2.0200_GP
       vdwparams%neff(7)=2.5200_GP
       vdwparams%neff(8)=2.6500_GP
       ! qoh: Neff from Halgren
       vdwparams%neff(9)=3.48_GP
       vdwparams%neff(15)=4.5_GP
       vdwparams%neff(16)=4.8_GP
       vdwparams%neff(17)=5.10_GP
       vdwparams%neff(35)=6.00_GP

    elseif ( dispersion == 4 ) then
! vama: added
!! VDW_DAMP_GRIMME_D2

       vdwparams%c6coeff( 1)=0.14_GP
       vdwparams%c6coeff( 2)=0.08_GP
       vdwparams%c6coeff( 3)=1.61_GP
       vdwparams%c6coeff( 4)=1.61_GP
       vdwparams%c6coeff( 5)=3.13_GP
       vdwparams%c6coeff( 6)=1.75_GP
       vdwparams%c6coeff( 7)=1.23_GP
       vdwparams%c6coeff( 8)=0.70_GP
       vdwparams%c6coeff( 9)=0.75_GP
       vdwparams%c6coeff(10)=0.63_GP
       vdwparams%c6coeff(11)=5.71_GP
       vdwparams%c6coeff(12)=5.71_GP
       vdwparams%c6coeff(13)=10.79_GP
       vdwparams%c6coeff(14)=9.23_GP
       vdwparams%c6coeff(15)=7.84_GP
       vdwparams%c6coeff(16)=5.57_GP
       vdwparams%c6coeff(17)=5.07_GP
       vdwparams%c6coeff(18)=4.61_GP
       do ii=19,30
            vdwparams%c6coeff(ii)=10.8_GP
       enddo
       vdwparams%c6coeff(31)=16.99_GP
       vdwparams%c6coeff(32)=17.10_GP
       vdwparams%c6coeff(33)=16.37_GP
       vdwparams%c6coeff(34)=12.64_GP
       vdwparams%c6coeff(35)=12.47_GP
       vdwparams%c6coeff(36)=12.01_GP
       do ii=37,48
            vdwparams%c6coeff(ii)=24.67_GP
       enddo
       vdwparams%c6coeff(49)=37.32_GP
       vdwparams%c6coeff(50)=38.71_GP
       vdwparams%c6coeff(51)=38.44_GP
       vdwparams%c6coeff(52)=31.74_GP
       vdwparams%c6coeff(53)=31.50_GP
       vdwparams%c6coeff(54)=29.99_GP
       vdwparams%c6coeff(55)=315.275_GP
       vdwparams%c6coeff(56)=226.994_GP
       vdwparams%c6coeff(57)=176.252_GP
       vdwparams%c6coeff(58)=140.68_GP
       vdwparams%c6coeff(59)=140.68_GP
       vdwparams%c6coeff(60)=140.68_GP
       vdwparams%c6coeff(61)=140.68_GP
       vdwparams%c6coeff(62)=140.68_GP
       vdwparams%c6coeff(63)=140.68_GP
       vdwparams%c6coeff(64)=140.68_GP
       vdwparams%c6coeff(65)=140.68_GP
       vdwparams%c6coeff(66)=140.68_GP
       vdwparams%c6coeff(67)=140.68_GP
       vdwparams%c6coeff(68)=140.68_GP
       vdwparams%c6coeff(69)=140.68_GP
       vdwparams%c6coeff(70)=140.68_GP
       vdwparams%c6coeff(71)=140.68_GP
       vdwparams%c6coeff(72)=105.112_GP
       vdwparams%c6coeff(73)=81.24_GP
       vdwparams%c6coeff(74)=81.24_GP
       vdwparams%c6coeff(75)=81.24_GP
       vdwparams%c6coeff(76)=81.24_GP
       vdwparams%c6coeff(77)=81.24_GP
       vdwparams%c6coeff(78)=81.24_GP
       vdwparams%c6coeff(79)=81.24_GP
       vdwparams%c6coeff(80)=57.364_GP
       vdwparams%c6coeff(81)=57.254_GP
       vdwparams%c6coeff(82)=63.162_GP
       vdwparams%c6coeff(83)=63.540_GP
       vdwparams%c6coeff(84)=55.283_GP
       vdwparams%c6coeff(85)=57.171_GP
       vdwparams%c6coeff(86)=56.64_GP 

       do ii=1,86
          vdwparams%c6coeff(ii)= vdwparams%c6coeff(ii)*jnm6mol2au
       end do

       vdwparams%radzero(1)=0.91_GP
       vdwparams%radzero(2)=0.92_GP
       vdwparams%radzero(3)=0.75_GP
       vdwparams%radzero(4)=1.28_GP
       vdwparams%radzero(5)=1.35_GP
       vdwparams%radzero(6)=1.32_GP
       vdwparams%radzero(7)=1.27_GP
       vdwparams%radzero(8)=1.22_GP
       vdwparams%radzero(9)=1.17_GP
       vdwparams%radzero(10)=1.13_GP
       vdwparams%radzero(11)=1.04_GP
       vdwparams%radzero(12)=1.24_GP
       vdwparams%radzero(13)=1.49_GP
       vdwparams%radzero(14)=1.56_GP
       vdwparams%radzero(15)=1.55_GP
       vdwparams%radzero(16)=1.53_GP
       vdwparams%radzero(17)=1.49_GP
       vdwparams%radzero(18)=1.45_GP
       vdwparams%radzero(19)=1.35_GP
       vdwparams%radzero(20)=1.34_GP
       vdwparams%radzero(21)=1.42_GP
       vdwparams%radzero(22)=1.42_GP
       vdwparams%radzero(23)=1.42_GP
       vdwparams%radzero(24)=1.42_GP
       vdwparams%radzero(25)=1.42_GP
       vdwparams%radzero(26)=1.42_GP
       vdwparams%radzero(27)=1.42_GP
       vdwparams%radzero(28)=1.42_GP
       vdwparams%radzero(29)=1.42_GP
       vdwparams%radzero(30)=1.42_GP
       vdwparams%radzero(31)=1.50_GP
       vdwparams%radzero(32)=1.57_GP
       vdwparams%radzero(33)=1.60_GP
       vdwparams%radzero(34)=1.61_GP
       vdwparams%radzero(35)=1.59_GP
       vdwparams%radzero(36)=1.57_GP
       vdwparams%radzero(37)=1.48_GP
       vdwparams%radzero(38)=1.46_GP
       vdwparams%radzero(39)=1.49_GP
       vdwparams%radzero(40)=1.49_GP
       vdwparams%radzero(41)=1.49_GP
       vdwparams%radzero(42)=1.49_GP
       vdwparams%radzero(43)=1.49_GP
       vdwparams%radzero(44)=1.49_GP
       vdwparams%radzero(45)=1.49_GP
       vdwparams%radzero(46)=1.49_GP
       vdwparams%radzero(47)=1.49_GP
       vdwparams%radzero(48)=1.49_GP
       vdwparams%radzero(49)=1.52_GP
       vdwparams%radzero(50)=1.64_GP
       vdwparams%radzero(51)=1.71_GP
       vdwparams%radzero(52)=1.72_GP
       vdwparams%radzero(53)=1.72_GP
       vdwparams%radzero(54)=1.71_GP
       vdwparams%radzero(54)=1.71_GP
       vdwparams%radzero(55)=1.638_GP
       vdwparams%radzero(56)=1.602_GP
       vdwparams%radzero(57)=1.564_GP
       do ii=58,71
          vdwparams%radzero(ii)= 1.594_GP
       end do
       vdwparams%radzero(72)=1.625_GP
       vdwparams%radzero(73)=1.611_GP
       vdwparams%radzero(74)=1.611_GP
       vdwparams%radzero(75)=1.611_GP
       vdwparams%radzero(76)=1.611_GP
       vdwparams%radzero(77)=1.611_GP
       vdwparams%radzero(78)=1.611_GP
       vdwparams%radzero(79)=1.611_GP
       vdwparams%radzero(80)=1.598_GP
       vdwparams%radzero(81)=1.805_GP
       vdwparams%radzero(82)=1.767_GP
       vdwparams%radzero(83)=1.725_GP
       vdwparams%radzero(84)=1.823_GP
       vdwparams%radzero(85)=1.810_GP
       vdwparams%radzero(86)=1.749_GP

       vdwparams%dcoeff(1)=0.0000_GP
       vdwparams%dcoeff(2)=0.0000_GP
       vdwparams%dcoeff(3)=0.0000_GP

       vdwparams%alpha=20.0000_GP

       select case (ixc)
!!s6
       case(11)            ! pbe
          vdwparams%s6=0.7500_GP
       case(14)            ! revpbe
          vdwparams%s6=1.2500_GP
       case(-406000)       ! pbe0
          vdwparams%s6=0.6000_GP
       case(-170000)       ! b97-d
          vdwparams%s6=1.2500_GP
       case(-106132)       ! b-p
          vdwparams%s6=1.0500_GP
       case(-416000)       ! b-lyp
          vdwparams%s6=1.2000_GP
       case(-402000)       ! b3-lyp
          vdwparams%s6=1.0500_GP
       case(-202231)       ! tpss
          vdwparams%s6=1.0000_GP
       case default  
          vdwparams%s6=1.0000_GP
       end select

!! radii scale
       vdwparams%radscale=1.1000_GP

!! scale radio
       do ii=1, 86
          vdwparams%radzero(ii)=vdwparams%radzero(ii)*vdwparams%radscale*angs2au
       enddo

    elseif ( dispersion == 5 ) then
      
       !allocate memory space for the parameters
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM) :: QATOM
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM) :: cov_table
!!$     integer, DIMENSION(MAX_ELEM)       :: MAXCN
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM,MAX_ELEM) :: R0AB
!!$     REAL(kind=GP), DIMENSION(MAX_ELEM,MAX_ELEM,MAX_CN,MAX_CN,3) :: C6AB
       vdwparams%qatom=f_malloc_ptr(max_elem,id='qatom')
       vdwparams%cov_table=f_malloc_ptr(max_elem,id='cov_table')
       vdwparams%maxcn=f_malloc_ptr(max_elem,id='maxcn')
       vdwparams%r0ab=f_malloc_ptr((/max_elem,max_elem/),id='r0ab')
       !vdwparams%c6ab=f_malloc0_ptr((/max_elem,max_elem,max_cn,max_cn,3/),id='c6ab')
! vama: added

!! VDW_DAMP_GRIMME_D3
       call init_cov_rad_d3
       call init_r0ab_d3
       call init_c6_params_d3
       call init_qat
       vdwparams%alpha=14.0000_GP

       vdwparams%s6=1.0000_GP
       vdwparams%s8=1.0000_GP
       vdwparams%sr6=1.0000_GP
       vdwparams%sr8=1.0000_GP
!!s6,s8,sr6,sr8
       select case (ixc)
       case(11,-101130)            ! pbe
          vdwparams%sr6 = 1.217_GP
          vdwparams%s8  = 0.722_GP
       case(14)            ! revpbe
          vdwparams%sr6 = 0.923
          vdwparams%s8 = 1.010
       case(-406)       ! pbe0
          vdwparams%sr6 = 1.278_GP
          vdwparams%s8  = 0.928_GP
       case(-170000)       ! b97-d
          vdwparams%sr6 = 0.892_GP
          vdwparams%s8  = 0.909_GP
       case(-106132)       ! b-p
          vdwparams%sr6 = 1.139_GP
          vdwparams%s8  = 1.683_GP
       case(-416)       ! b-lyp
          vdwparams%sr6 = 1.094_GP
          vdwparams%s8  = 1.682_GP
       case(-402)       ! b3-lyp
          vdwparams%sr6 = 1.261_GP
          vdwparams%s8  = 1.703_GP
       case(-202231)       ! tpss
          vdwparams%sr6 = 1.252_GP
          vdwparams%s8  = 1.242_GP
       case default  
          vdwparams%s6=1.0000_GP
          vdwparams%sr8=1.0000_GP
          vdwparams%sr6=1.0000_GP
       end select
    endif

  END SUBROUTINE vdwcorrection_initializeparams

  !> Free the parameters of the module
  !! this routine is useful only for the D3 correction
  subroutine vdwcorrection_freeparams()
    implicit none
    
    if (associated(vdwparams%qatom)) call f_free_ptr(vdwparams%qatom)
    if (associated(vdwparams%cov_table)) call f_free_ptr(vdwparams%cov_table)
    if (associated(vdwparams%maxcn)) call f_free_ptr(vdwparams%maxcn)
    if (associated(vdwparams%r0ab)) call f_free_ptr(vdwparams%r0ab)
    if (associated(vdwparams%coeffs)) call f_free_ptr(vdwparams%coeffs)

  end subroutine vdwcorrection_freeparams


  !< This subroutine calculates the dispersion correction to the total energy.                                                    !
  !! @author
  !! Written by Quintin Hill in 2008, with some modifications in 2008
  !! Modified for BigDFT in March/April 2009 by Quintin Hill.
  subroutine vdwcorrection_calculate_energy(dispersion_energy,rxyz,atoms,dispersion)!,iproc)

    use module_types

    implicit none

    ! Arguments
    type(atoms_data),                 intent(in) :: atoms
    real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
    !integer,                          intent(in) :: iproc
    integer,                          intent(in) :: dispersion
    real(kind=GP),                   intent(out) :: dispersion_energy

    ! Internal variables
    integer       :: atom1, atom2   ! atom counters for loops
    integer       :: nzatom1   ! Atomic number of atom 1
    integer       :: nzatom2   ! Atomic number of atom 2
    real(kind=GP) :: distance  ! Distance between pairs of atoms
    real(kind=GP) :: sqdist    ! Square distance between pairs of atoms
    real(kind=GP) :: c6coeff   ! The c6coefficient of the pair
    real(kind=GP) :: damping   ! The damping for the pair
    real(kind=GP) :: c6d3,cni,cnj,e8,e6,c8
    real(kind=GP) :: dum

    dispersion_energy = 0.0_GP
    dum = 0.0_GP
!!vama!!          write(*,'(1x,a, e12.5)') &
!!vama!!               's6 : ', vdwparams%s6
!!vama!!          write(*,'(1x,a, e12.5)') &
!!vama!!               'sr6 : ', vdwparams%sr6
!!vama!!          write(*,'(1x,a, e12.5)') &
!!vama!!               'sr8 : ', vdwparams%sr8
!!vama!!
    if (dispersion /= VDW_NONE) then 

       ! qoh: Loop over all distinct pairs of atoms

       do atom1=1,atoms%astruct%nat
          do atom2=1,atom1-1

             nzatom1 = atoms%nzatom(atoms%astruct%iatype(atom1))
             nzatom2 = atoms%nzatom(atoms%astruct%iatype(atom2))

             ! qoh: Calculate distance between each pair of atoms
             sqdist=(rxyz(1,atom1) - rxyz(1,atom2))**2 &
                  + (rxyz(2,atom1) - rxyz(2,atom2))**2 &
                  + (rxyz(3,atom1) - rxyz(3,atom2))**2 
             distance = sqrt(sqdist)
             ! qoh: distance**6 = sqdist**3


             if (sqdist < 20000_GP) then
                if (dispersion < 5) then
!! DFT -D and -D2
                   ! qoh : Get damping function
                   damping = vdwcorrection_damping(nzatom1,nzatom2,distance,dispersion)

                   ! qoh: Calculate c6 coefficient
                   c6coeff = vdwcorrection_c6(nzatom1,nzatom2,dispersion)
                endif

                if (dispersion == 4) then
!! DFT-D2
                   dispersion_energy = dispersion_energy &
                    - (c6coeff * damping / sqdist**3)*vdwparams%s6

                else if (dispersion == 5) then
                  e6 = 0.0_GP
                  e8 = 0.0_GP
!! DFT-D3
             
                   cni  = crd_nr(atom1,atoms%astruct%nat,rxyz,atoms)
                   cnj  = crd_nr(atom2,atoms%astruct%nat,rxyz,atoms)
                   c6d3 = c6cn(nzatom1,nzatom2,cni,cnj)
                   c8   = 3.0_GP*c6d3*vdwparams%Qatom(nzatom1)*vdwparams%Qatom(nzatom2)
    
                   damping  = 1.0_GP/(1.0_GP &
                      +6.0_GP*((distance/(vdwparams%r0AB(nzatom1,nzatom2)*&
                      vdwparams%sr6))**(-vdwparams%alpha)))
    
                   e6 = c6d3 * damping / sqdist**3 
    
                   damping = 1.0_GP/(1.0_GP &
                      + 6.0_GP*((distance/(vdwparams%r0AB(nzatom1,nzatom2)*&
                      vdwparams%sr8))**(-vdwparams%alpha-2.0_GP)))

    
                   e8 = c8*damping/ sqdist**4
                   dispersion_energy = dispersion_energy &
                                      - (e6*vdwparams%s6 + e8*vdwparams%s8 )
                else
!! DFT-D

                   dispersion_energy = dispersion_energy &
                      - (c6coeff * damping / sqdist**3)

                endif
             endif
!!vama!!          write(*,'(1x,a, e12.5,1x,a)') &
!!vama!!               'Dispersion Correction Energy: ', dispersion_energy, 'Hartree'
          enddo
       enddo

!!$       if (iproc == 0) then
!!$          write(*,'(1x,a, e12.5,1x,a)') &
!!$               'Dispersion Correction Energy: ', dispersion_energy, 'Hartree'
!!$       end if
    end if

  END SUBROUTINE vdwcorrection_calculate_energy


  !< This subroutine calculates the dispersion correction to the total energy.i
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  !! Modified for BigDFT in March/April 2009 by Quintin Hill.
  SUBROUTINE vdwcorrection_calculate_forces(vdw_forces,rxyz,atoms,dispersion) 

    use module_types

    implicit none

    ! Arguments

    type(atoms_data),                 intent(in)  :: atoms
    real(GP), dimension(3,atoms%astruct%nat), intent(out) :: vdw_forces
    real(GP), dimension(3,atoms%astruct%nat), intent(in)  :: rxyz
    integer,                          intent(in)  :: dispersion

    ! Internal variables

    integer       :: atom1,atom2 ! atom counters for loops
    integer       :: nzatom1     ! Atomic number of atom 1
    integer       :: nzatom2     ! Atomic number of atom 2
    integer       :: nzatom3     ! Atomic number of atom 3
    real(kind=GP) :: distance    ! Distance between pairs of atoms
    real(kind=GP) :: sqdist      ! Square distance between pairs of atoms
    real(kind=GP) :: c6coeff     ! The c6coefficient of the pair
    real(kind=GP) :: damping     ! The damping for the pair
    real(kind=GP) :: dampingdrv  ! The damping derivative for the pair
    real(kind=GP) :: drvcommon   ! The common part of the derivative


    ! vama D3
    integer       :: atom3       ! atom counters for loops
    real(kind=GP) :: dxAj, dyAj, dzAj
    real(kind=GP) :: dxAk, dyAk, dzAk
    real(kind=GP) :: dxjk, dyjk, dzjk
    real(kind=GP) :: r0aj, r0jk, r0ak
    real(kind=GP) :: rAj,rjk,rAk
    real(kind=GP) :: Qfac
    real(kind=GP) :: fac6, fac8
    real(kind=GP) :: fdmp6, fdmp8
    real(kind=GP) :: tmp6, tmp8
    real(kind=GP) :: tmp6a, tmp8a
    real(kind=GP) :: cnA, cnj, c6Aj
    real(kind=GP), DIMENSION(3,atoms%astruct%nat)            :: cnij
    real(kind=GP), DIMENSION(3,atoms%astruct%nat,atoms%astruct%nat)  :: cnijk
    real(kind=GP), DIMENSION(3)                      :: grad_c6

    vdw_forces = 0.0_GP

    if (dispersion /= VDW_NONE) then 

    if (dispersion < 5 ) then 

       ! qoh: Loop over all atoms

       do atom1=1,atoms%astruct%nat

          ! qoh: Loop over all other atoms

          do atom2=1,atoms%astruct%nat

             if ( atom2 .ne. atom1) then

                nzatom1 = atoms%nzatom(atoms%astruct%iatype(atom1))
                nzatom2 = atoms%nzatom(atoms%astruct%iatype(atom2))

                ! qoh: Calculate c6 coefficient
                c6coeff = vdwcorrection_c6(nzatom1,nzatom2,dispersion)

                ! qoh: Calculate distance between each pair of atoms
                sqdist=(rxyz(1,atom2) - rxyz(1,atom1))**2 &
                     + (rxyz(2,atom2) - rxyz(2,atom1))**2 &
                     + (rxyz(3,atom2) - rxyz(3,atom1))**2 
                distance = sqrt(sqdist)

                ! qoh : Get damping function
                damping = vdwcorrection_damping(nzatom1,nzatom2,distance,dispersion)

                dampingdrv = &
                     vdwcorrection_drvdamping(nzatom1,nzatom2,distance,dispersion)

                ! qoh: distance**6 = sqdist**3
                if (dispersion == 4) then
                drvcommon = vdwparams%s6*(c6coeff/sqdist**3)*&
                     (dampingdrv-6.0_GP*damping/sqdist)
                else
                drvcommon = (c6coeff/sqdist**3)*&
                     (dampingdrv-6.0_GP*damping/sqdist)
                end if

                ! cks: this is the *negative* of the derivative of the
                ! cks: dispersion energy w.r.t. atomic coordinates
                vdw_forces(:,atom1) = vdw_forces(:,atom1) +drvcommon*&
                     (rxyz(:,atom1) - rxyz(:,atom2))

             end if

          enddo
       enddo
    else if (dispersion == 5) then

       call crd_nr_der(atoms%astruct%nat,rxyz,cnij,cnijk,atoms)

       do atom1=1,atoms%astruct%nat

          do atom2=1,atoms%astruct%nat

             if (atom2 .eq. atom1) cycle

             nzatom1 = atoms%nzatom(atoms%astruct%iatype(atom1))
             nzatom2 = atoms%nzatom(atoms%astruct%iatype(atom2))

             dxAj = rxyz(1,atom1) - rxyz(1,atom2)
             dyAj = rxyz(2,atom1) - rxyz(2,atom2)
             dzAj = rxyz(3,atom1) - rxyz(3,atom2)
             rAj = dxAj**2+dyAj**2+dzAj**2

             distance = sqrt(rAj)

!         Two center derivatives. Grimme uses screening to reduce 
!         computational work
!         Screening r^2 distance vs threshold of 20000.0
             if (rAj.gt.20000.0_GP) cycle
 
!         Factors
 
             r0aj = vdwparams%r0AB(nzatom1,nzatom2)
             Qfac = vdwparams%Qatom(nzatom1)*vdwparams%Qatom(nzatom2)
             fac6 = (distance/(vdwparams%sr6*r0aj))**(-vdwparams%alpha)
             fac8 = (distance/(vdwparams%sr8*r0aj))**(-(vdwparams%alpha+2.0_GP))
             fdmp6 =1.0_GP/(1.0_GP+6.0_GP*fac6)
             fdmp8 =1.0_GP/(1.0_GP+6.0_GP*fac8)
 
!         Coordination dependent C6_AB value
 
             cnA = crd_nr(atom1,atoms%astruct%nat,rxyz,atoms)
             cnj = crd_nr(atom2,atoms%astruct%nat,rxyz,atoms)
             c6Aj = c6cn(nzatom1,nzatom2,cnA,cnj)
!
!         Get gradient for coordination number dependent C6
!
             call c6_grad(grad_c6,atom1,atom2,atom1,rxyz,atoms%astruct%nat,cnij,cnijk,atoms)
   
             tmp6 = 6.0_GP*fdmp6*vdwparams%s6*c6Aj/(rAj**4.0_GP)
             tmp8 = 6.0_GP*fdmp8*vdwparams%s8*c6Aj*Qfac/(rAj**5.0_GP)

!         dx contribution to A

             tmp6a = tmp6*dxAj
             tmp8a = tmp8*dxAj
             vdw_forces(1,atom1) = vdw_forces(1,atom1) -(&
               +(1.0_GP-fdmp6*fac6*vdwparams%alpha)*tmp6a&
               -fdmp6*vdwparams%s6*grad_c6(1)/(rAj**3.0_GP)&
               +(4.0_GP-3.0_GP*fdmp8*fac8*(vdwparams%alpha+2.0_GP))*tmp8a&
               -3.0_GP*fdmp8*vdwparams%s8*grad_c6(1)*Qfac/(rAj**4.0_GP))

!         dy contribution to A

             tmp6a = tmp6*dyAj
             tmp8a = tmp8*dyAj
             vdw_forces(2,atom1) = vdw_forces(2,atom1) - (&
               +(1.0_GP-fdmp6*fac6*vdwparams%alpha)*tmp6a&
               -fdmp6*vdwparams%s6*grad_c6(2)/(rAj**3.0_GP)&
               +(4.0_GP-3.0_GP*fdmp8*fac8*(vdwparams%alpha+2.0_GP))*tmp8a&
               -3.0_GP*fdmp8*vdwparams%s8*grad_c6(2)*Qfac/(rAj**4.0_GP))
 
!         dz contribution to A
 
             tmp6a = tmp6*dzAj
             tmp8a = tmp8*dzAj
  
             vdw_forces(3,atom1) = vdw_forces(3,atom1) -(&
               +(1.0_GP-fdmp6*fac6*vdwparams%alpha)*tmp6a&
               -fdmp6*vdwparams%s6*grad_c6(3)/(rAj**3.0_GP)&
               +(4.0_GP-3.0_GP*fdmp8*fac8*(vdwparams%alpha+2.0_GP))*tmp8a&
               -3.0_GP*fdmp8*vdwparams%s8*grad_c6(3)*Qfac/(rAj**4.0_GP))
          enddo

!         Three center derivatives. Grimme uses aggressive screening
!         to get this N^3 contribution back to N^2

          do atom2=2,atoms%astruct%nat
             if (atom2 .eq. atom1) cycle
             rAj = sqrt((rxyz(1,atom1) - rxyz(1,atom2))**2 &
                    + (rxyz(2,atom1) - rxyz(2,atom2))**2 &
                    + (rxyz(3,atom1) - rxyz(3,atom2))**2)

             r0aj = vdwparams%r0AB(nzatom1,nzatom2)

!            Screening per Grimme

             if (rAj >= 1600.00_GP*r0aj/vdwparams%r0AB(1,1)) cycle

!            Third center involved
!
             do atom3=1,atom2-1
                 nzatom3 = atoms%nzatom(atoms%astruct%iatype(atom3))
                 if (atom1.eq.atom3) cycle

                 dxAk = rxyz(1,atom1) - rxyz(1,atom3)
                 dyAk = rxyz(2,atom1) - rxyz(2,atom3)
                 dzAk = rxyz(3,atom1) - rxyz(3,atom3)
                 rAk = dxAk**2+dyAk**2+dzAk**2

                 r0ak = vdwparams%r0AB(nzatom1,nzatom3)

                 dxjk = rxyz(1,atom2) - rxyz(1,atom3)
                 dyjk = rxyz(2,atom2) - rxyz(2,atom3)
                 dzjk = rxyz(3,atom2) - rxyz(3,atom3)
                 rjk = dxjk**2+dyjk**2+dzjk**2

                 r0jk = vdwparams%r0AB(nzatom2,nzatom3)
!
!                Screening r^2 distance vs threshold of 1600.0*(radii Ak)
!
                 if ((rAk > 1600.0_GP*r0ak/vdwparams%r0AB(1,1)) .OR.&
                        (rjk > 1600.0_GP*r0jk/vdwparams%r0AB(1,1))) cycle
!
!                Get gradient for coordination number dependent C6 for three centers
!
                 call c6_grad(grad_c6,atom2,atom3,atom1,rxyz,atoms%astruct%nat,cnij,cnijk,atoms)

                 fac6 = (vdwparams%sr6*r0jk/sqrt(rjk))**(vdwparams%alpha)
                 fac8 = (vdwparams%sr8*r0jk/sqrt(rjk))**(vdwparams%alpha+2.0_GP)
                 fdmp6 = 1.0_GP/(1.0_GP+6.0_GP*fac6)
                 fdmp8 = 1.0_GP/(1.0_GP+6.0_GP*fac8)
 
!                dx, dy, and dz contribution to A
 
                 Qfac = vdwparams%Qatom(nzatom2)*vdwparams%Qatom(nzatom3)

                 vdw_forces(1,atom1) = vdw_forces(1,atom1)-(&
                       -fdmp6*vdwparams%s6*grad_c6(1)/(rjk**3.0_GP)&
                       -3.0_GP*fdmp8*vdwparams%s8*grad_c6(1)*Qfac/(rjk**4.0_GP))

                 vdw_forces(2,atom1) = vdw_forces(2,atom1)-(&
                       -fdmp6*vdwparams%s6*grad_c6(2)/(rjk**3.0_GP)&
                       -3.0_GP*fdmp8*vdwparams%s8*grad_c6(2)*Qfac/(rjk**4.0_GP))

                 vdw_forces(3,atom1)=vdw_forces(3,atom1)-(&
                       -fdmp6*vdwparams%s6*grad_c6(3)/(rjk**3.0_GP)&
                       -3.0_GP*fdmp8*vdwparams%s8*grad_c6(3)*Qfac/(rjk**4.0_GP))
             enddo
          enddo
        enddo
    end if
    end if


!        write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom vdw'
!        do atom1=1,atoms%astruct%nat
!           write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
!           atom1,trim(atoms%astruct%atomnames(atoms%astruct%iatype(atom1))),(vdw_forces(atom2,atom1),atom2=1,3)
!        enddo



  END SUBROUTINE vdwcorrection_calculate_forces
  
  !< This subroutine warns about the use of unoptimised or unavailable dispersion parameters.
  !! @author
  ! Written by Quintin Hill on 13/02/2009.
  ! Quintin Hill added functional check on 23/02/2009.
  subroutine vdwcorrection_warnings(atoms,in)

    use module_types, only: input_variables, atoms_data

    implicit none

    type(atoms_data),        intent(in) :: atoms
    type(input_variables),   intent(in) :: in

    integer            :: itype ! Type counter
    ! qoh: Atomic numbers of elements with unoptimised dispersion parameters:
    integer, parameter :: unoptimised(4) = (/9,15,17,35/)
    ! qoh: Atomic numbers of elements with optimised dispersion parameters:
    integer, parameter :: optimised(5) = (/1,6,7,8,16/)
    ! qoh: XC functionals with optimised parameters:
    integer, parameter :: xcfoptimised(4) = (/11,14,15,200/)

    if (in%dispersion /= 0) then 

       ! qoh: Loop over types to check we have parameters
       do itype=1,atoms%astruct%ntypes
          if (any(unoptimised == atoms%nzatom(itype)) .and. &
               any(xcfoptimised == in%ixc)) then 
             write(*,'(a,a2)') 'WARNING: Unoptimised dispersion &
                  &parameters used for ', atoms%astruct%atomnames(itype)
          elseif (.not. any(optimised == atoms%nzatom(itype)) .and. &
               .not. any(unoptimised == atoms%nzatom(itype))) then
             write(*,'(a,a2)') 'WARNING: No dispersion parameters &
                  &available for ', atoms%astruct%atomnames(itype) 
          end if
       end do

       if (.not. any(xcfoptimised == in%ixc)) &
            write(*,'(a,i2)') 'WARNING: No optimised dispersion parameters &
            &available for ixc=', in%ixc
    end if
  END SUBROUTINE vdwcorrection_warnings


  !< This function calculates a heteroatomic C_6 coefficient from
  !! homoatomic C_6 coefficients using the formula given in Elstner's
  !! paper (J. Chem. Phys. 114(12), 5149–5155). 
  !! (Corrected error in the formula appearing in this paper.) 
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  function vdwcorrection_c6(nzatom1,nzatom2,dispersion)

    implicit none

    real(kind=GP):: vdwcorrection_c6
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2
    integer, intent(in) :: dispersion
    real(kind=GP) :: c61 ! c6 coefficient of atom 1
    real(kind=GP) :: c62 ! c6 coefficient of atom 2
    real(kind=GP) :: ne1 ! Effective number of electrons for atom 1
    real(kind=GP) :: ne2 ! Effective number of electrons for atom 2
    real(kind=GP), parameter :: third = 1.0_GP/3.0_GP

    ! qoh: Set up shorthands

    c61 = vdwparams%c6coeff(nzatom1)
    c62 = vdwparams%c6coeff(nzatom2)    
    ne1 = vdwparams%neff(nzatom1)
    ne2 = vdwparams%neff(nzatom2)

    if (dispersion == 4) then
    vdwcorrection_c6 = sqrt(c61*c62)
    else
    if (c61 .lt. epsilon(1.0_GP) .or. c62 .lt.epsilon(1.0_GP) .or. &
         ne1 .lt. epsilon(1.0_GP) .or. ne2 .lt. epsilon(1.0_GP)) then
       vdwcorrection_c6=0.0_GP
    else
       vdwcorrection_c6 = 2.0_GP * (c61**2 * c62**2 * ne1 * ne2)**(third)/&
            (((c61 * ne2**2)**(third)) + ((c62* ne1**2)**(third)))
    end if
    end if

  end function vdwcorrection_c6


  !< This function calculates the damping function specified by in%dispersion:
  !! (1) Damping funtion from Elstner                  (J. Chem. Phys. 114(12), 5149-5155).                         
  !! (2) First damping function from Wu and Yang (I)   (J. Chem. Phys. 116(2), 515-524).
  !! (3) Second damping function from Wu and Yang (II) (J. Chem. Phys. 116(2), 515–524).
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  !! Merged into a single function in July 2008
  function vdwcorrection_damping(nzatom1,nzatom2,separation,dispersion)

    implicit none

    real(kind=GP)             :: vdwcorrection_damping
    integer, intent(in)       :: dispersion
    integer, intent(in)       :: nzatom1 ! Atomic number of atom 1
    integer, intent(in)       :: nzatom2 ! Atomic number of atom 2
    real(kind=GP), intent(in) :: separation
    real(kind=GP)             :: radzero
    real(kind=GP)             :: expo ! Exponent
    integer, parameter        :: mexpo(2) = (/4,2/)
    integer, parameter        :: nexpo(2) = (/7,3/)
    integer                   :: mexp
    integer                   :: nexp    

    radzero = vdwcorrection_radzero(nzatom1,nzatom2,dispersion)

    select case (dispersion)
    case(1,2)

       mexp = mexpo(dispersion)
       nexp = nexpo(dispersion) 

       expo = -vdwparams%dcoeff(dispersion)*(separation/radzero)**nexp
       vdwcorrection_damping = ( 1.0_GP - exp(expo))**mexp

    case(3)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_GP) 
       vdwcorrection_damping = 1.0_GP/( 1.0_GP + exp(expo))

    case(4)

       expo = -vdwparams%alpha*((separation/radzero)-1.0_GP) 
       vdwcorrection_damping = 1.0_GP/( 1.0_GP + exp(expo))
       
    case default
       vdwcorrection_damping = 1.0_GP
    end select

  end function vdwcorrection_damping


  !< This function calculates the derivative with respect to atomic
  !!coordinates of the damping function specified by int%dispersion:  
  !! (1) Damping funtion from Elstner
  !!     (J. Chem. Phys. 114(12), 5149-5155).\n
  !! (2) First damping function from Wu and Yang (I)
  !!     (J. Chem. Phys. 116(2), 515-524).\n
  !! (3) Second damping function from Wu and Yang (II)
  !!     (J. Chem. Phys. 116(2), 515–524).\n
  !! (4) S6 method from S. Grimme 
  !!     (J. Comp. Chem. 27, 1787 (2006).\n
  !! Note: For simplicity the @f$(r_{A,i}-r_{B_i})@f$ (where i is x,y or z)
  !!       term is omitted here and included in the calling routine.
  !! @author
  !! Written by Quintin Hill in July 2008.
  function vdwcorrection_drvdamping(nzatom1,nzatom2,separation,dispersion)

    use module_types, only: input_variables

    implicit none

    real(kind=GP)                     :: vdwcorrection_drvdamping
    integer, intent(in)       :: dispersion
    integer,               intent(in) :: nzatom1 ! Atomic number of atom 1
    integer,               intent(in) :: nzatom2 ! Atomic number of atom 2
    real(kind=GP),         intent(in) :: separation
    real(kind=GP)                     :: radzero
    real(kind=GP)                     :: expo ! Exponent
    integer, parameter                :: mexpo(2) = (/4,2/)
    integer, parameter                :: nexpo(2) = (/7,3/)
    integer                           :: mexp
    integer                           :: nexp    

    radzero = vdwcorrection_radzero(nzatom1,nzatom2,dispersion)

    select case (dispersion)
    case(1,2)

       mexp = mexpo(dispersion)
       nexp = nexpo(dispersion)   

       expo = -vdwparams%dcoeff(dispersion)*(separation/radzero)**nexp

       vdwcorrection_drvdamping = real((mexp*nexp),kind=GP)*exp(expo)* &
            vdwparams%dcoeff(dispersion)*( 1.0_GP - exp(expo))**(mexp-1)*&
            separation**(nexp-2)/radzero**nexp

    case(3)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_GP)

       vdwcorrection_drvdamping = vdwparams%dcoeff(3)*exp(expo)/&
            (separation*radzero*( 1.0_GP + exp(expo))**2)

    case(4)

       expo = -vdwparams%alpha*((separation/radzero)-1.0_GP)

       vdwcorrection_drvdamping = vdwparams%alpha*exp(expo)/&
            (separation*radzero*( 1.0_GP + exp(expo))**2)

    case default
       vdwcorrection_drvdamping = 0.0_GP
    end select

  end function vdwcorrection_drvdamping


  !< Function to calculate the R_0 for an atom pair. Uses expression
  !! found in Elstner's paper (J. Chem. Phys. 114(12), 5149–5155).
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  function vdwcorrection_radzero(nzatom1,nzatom2,dispersion)

    implicit none
    real(kind=GP) :: vdwcorrection_radzero
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2
    integer, intent(in) :: dispersion
    REAL(kind=GP), PARAMETER :: ANGSTROM=1.889726313_GP   

    if (dispersion == 4) then
    vdwcorrection_radzero= (vdwparams%radzero(nzatom1)+vdwparams%radzero(nzatom2))
    else
    vdwcorrection_radzero = ANGSTROM*(vdwparams%radzero(nzatom1)**3 + &
         vdwparams%radzero(nzatom2)**3)/&
         (vdwparams%radzero(nzatom1)**2 + vdwparams%radzero(nzatom2)**2)
    end if

  end function vdwcorrection_radzero

  subroutine init_cov_rad_d3
      implicit none
      integer            :: i
!     Grimme:
!     covalent radii (taken from.or.ivdw.eq.4 Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
!     values for metals decreased by 10 %
      vdwparams%cov_table=(/&
      &0.32_GP, 0.46_GP, 1.20_GP, 0.94_GP, 0.77_GP, 0.75_GP, 0.71_GP, 0.63_GP, 0.64_GP, 0.67_GP,&
      &1.40_GP, 1.25_GP, 1.13_GP, 1.04_GP, 1.10_GP, 1.02_GP, 0.99_GP, 0.96_GP, 1.76_GP, 1.54_GP,&
      &1.33_GP, 1.22_GP, 1.21_GP, 1.10_GP, 1.07_GP, 1.04_GP, 1.00_GP, 0.99_GP, 1.01_GP, 1.09_GP,&
      &1.12_GP, 1.09_GP, 1.15_GP, 1.10_GP, 1.14_GP, 1.17_GP, 1.89_GP, 1.67_GP, 1.47_GP, 1.39_GP,&
      &1.32_GP, 1.24_GP, 1.15_GP, 1.13_GP, 1.13_GP, 1.08_GP, 1.15_GP, 1.23_GP, 1.28_GP, 1.26_GP,&
      &1.26_GP, 1.23_GP, 1.32_GP, 1.31_GP, 2.09_GP, 1.76_GP, 1.62_GP, 1.47_GP, 1.58_GP, 1.57_GP,&
      &1.56_GP, 1.55_GP, 1.51_GP, 1.52_GP, 1.51_GP, 1.50_GP, 1.49_GP, 1.49_GP, 1.48_GP, 1.53_GP,&
      &1.46_GP, 1.37_GP, 1.31_GP, 1.23_GP, 1.18_GP, 1.16_GP, 1.11_GP, 1.12_GP, 1.13_GP, 1.32_GP,&
      &1.30_GP, 1.30_GP, 1.36_GP, 1.31_GP, 1.38_GP, 1.42_GP, 2.01_GP, 1.81_GP, 1.67_GP, 1.58_GP,&
      &1.52_GP, 1.53_GP, 1.54_GP, 1.55_GP /)
      do i=1,MAX_ELEM
         vdwparams%cov_table(i)=vdwparams%cov_table(i)*k2*angs2au
      end do
  END SUBROUTINE init_cov_rad_d3

  subroutine init_r0ab_d3
      implicit none
      !REAL(kind=GP), DIMENSION(4465) :: r0ab_ref_table_grimme
      integer, parameter :: refdim=4465
      real(gp), dimension(:), allocatable :: r0ab_ref_table_grimme
      integer            :: i,j,k
 
      r0ab_ref_table_grimme=f_malloc(refdim,id='r0ab_ref_table_grimme')

       r0ab_ref_table_grimme(   1:  70)=(/&
         2.1823_GP,  1.8547_GP,  1.7347_GP,  2.9086_GP,  2.5732_GP,  3.4956_GP,   2.3550_gp&
      &,  2.5095_GP,  2.9802_GP,  3.0982_GP,  2.5141_GP,  2.3917_GP,  2.9977_GP,  2.9484_gp&
      &,  3.2160_GP,  2.4492_GP,  2.2527_GP,  3.1933_GP,  3.0214_GP,  2.9531_GP,  2.9103_gp&
      &,  2.3667_GP,  2.1328_GP,  2.8784_GP,  2.7660_GP,  2.7776_GP,  2.7063_GP,  2.6225_gp&
      &,  2.1768_GP,  2.0625_GP,  2.6395_GP,  2.6648_GP,  2.6482_GP,  2.5697_GP,  2.4846_gp&
      &,  2.4817_GP,  2.0646_GP,  1.9891_GP,  2.5086_GP,  2.6908_GP,  2.6233_GP,  2.4770_gp&
      &,  2.3885_GP,  2.3511_GP,  2.2996_GP,  1.9892_GP,  1.9251_GP,  2.4190_GP,  2.5473_gp&
      &,  2.4994_GP,  2.4091_GP,  2.3176_GP,  2.2571_GP,  2.1946_GP,  2.1374_GP,  2.9898_gp&
      &,  2.6397_GP,  3.6031_GP,  3.1219_GP,  3.7620_GP,  3.2485_GP,  2.9357_GP,  2.7093_gp&
      &,  2.5781_GP,  2.4839_GP,  3.7082_GP,  2.5129_GP,  2.7321_GP,  3.1052_GP,  3.2962_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(  71: 140)=(/&                                            
         3.1331_GP,  3.2000_GP,  2.9586_GP,  3.0822_GP,  2.8582_GP,  2.7120_GP,   3.2570_gp&
      &,  3.4839_GP,  2.8766_GP,  2.7427_GP,  3.2776_GP,  3.2363_GP,  3.5929_GP,  3.2826_gp&
      &,  3.0911_GP,  2.9369_GP,  2.9030_GP,  2.7789_GP,  3.3921_GP,  3.3970_GP,  4.0106_gp&
      &,  2.8884_GP,  2.6605_GP,  3.7513_GP,  3.1613_GP,  3.3605_GP,  3.3325_GP,  3.0991_gp&
      &,  2.9297_GP,  2.8674_GP,  2.7571_GP,  3.8129_GP,  3.3266_GP,  3.7105_GP,  3.7917_gp&
      &,  2.8304_GP,  2.5538_GP,  3.3932_GP,  3.1193_GP,  3.1866_GP,  3.1245_GP,  3.0465_gp&
      &,  2.8727_GP,  2.7664_GP,  2.6926_GP,  3.4608_GP,  3.2984_GP,  3.5142_GP,  3.5418_gp&
      &,  3.5017_GP,  2.6190_GP,  2.4797_GP,  3.1331_GP,  3.0540_GP,  3.0651_GP,  2.9879_gp&
      &,  2.9054_GP,  2.8805_GP,  2.7330_GP,  2.6331_GP,  3.2096_GP,  3.5668_GP,  3.3684_gp&
      &,  3.3686_GP,  3.3180_GP,  3.3107_GP,  2.4757_GP,  2.4019_GP,  2.9789_GP,  3.1468_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 141: 210)=(/&                                            
         2.9768_GP,  2.8848_GP,  2.7952_GP,  2.7457_GP,  2.6881_GP,  2.5728_GP,   3.0574_gp&
      &,  3.3264_GP,  3.3562_GP,  3.2529_GP,  3.1916_GP,  3.1523_GP,  3.1046_GP,  2.3725_gp&
      &,  2.3289_GP,  2.8760_GP,  2.9804_GP,  2.9093_GP,  2.8040_GP,  2.7071_GP,  2.6386_gp&
      &,  2.5720_GP,  2.5139_GP,  2.9517_GP,  3.1606_GP,  3.2085_GP,  3.1692_GP,  3.0982_gp&
      &,  3.0352_GP,  2.9730_GP,  2.9148_GP,  3.2147_GP,  2.8315_GP,  3.8724_GP,  3.4621_gp&
      &,  3.8823_GP,  3.3760_GP,  3.0746_GP,  2.8817_GP,  2.7552_GP,  2.6605_GP,  3.9740_gp&
      &,  3.6192_GP,  3.6569_GP,  3.9586_GP,  3.6188_GP,  3.3917_GP,  3.2479_GP,  3.1434_gp&
      &,  4.2411_GP,  2.7597_GP,  3.0588_GP,  3.3474_GP,  3.6214_GP,  3.4353_GP,  3.4729_gp&
      &,  3.2487_GP,  3.3200_GP,  3.0914_GP,  2.9403_GP,  3.4972_GP,  3.7993_GP,  3.6773_gp&
      &,  3.8678_GP,  3.5808_GP,  3.8243_GP,  3.5826_GP,  3.4156_GP,  3.8765_GP,  4.1035_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 211: 280)=(/&                                            
         2.7361_GP,  2.9765_GP,  3.2475_GP,  3.5004_GP,  3.4185_GP,  3.4378_GP,   3.2084_gp&
      &,  3.2787_GP,  3.0604_GP,  2.9187_GP,  3.4037_GP,  3.6759_GP,  3.6586_GP,  3.8327_gp&
      &,  3.5372_GP,  3.7665_GP,  3.5310_GP,  3.3700_GP,  3.7788_GP,  3.9804_GP,  3.8903_gp&
      &,  2.6832_GP,  2.9060_GP,  3.2613_GP,  3.4359_GP,  3.3538_GP,  3.3860_GP,  3.1550_gp&
      &,  3.2300_GP,  3.0133_GP,  2.8736_GP,  3.4024_GP,  3.6142_GP,  3.5979_GP,  3.5295_gp&
      &,  3.4834_GP,  3.7140_GP,  3.4782_GP,  3.3170_GP,  3.7434_GP,  3.9623_GP,  3.8181_gp&
      &,  3.7642_GP,  2.6379_GP,  2.8494_GP,  3.1840_GP,  3.4225_GP,  3.2771_GP,  3.3401_gp&
      &,  3.1072_GP,  3.1885_GP,  2.9714_GP,  2.8319_GP,  3.3315_GP,  3.5979_GP,  3.5256_gp&
      &,  3.4980_GP,  3.4376_GP,  3.6714_GP,  3.4346_GP,  3.2723_GP,  3.6859_GP,  3.8985_gp&
      &,  3.7918_GP,  3.7372_GP,  3.7211_GP,  2.9230_GP,  2.6223_GP,  3.4161_GP,  2.8999_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 281: 350)=(/&                                            
         3.0557_GP,  3.3308_GP,  3.0555_GP,  2.8508_GP,  2.7385_GP,  2.6640_GP,   3.5263_gp&
      &,  3.0277_GP,  3.2990_GP,  3.7721_GP,  3.5017_GP,  3.2751_GP,  3.1368_GP,  3.0435_gp&
      &,  3.7873_GP,  3.2858_GP,  3.2140_GP,  3.1727_GP,  3.2178_GP,  3.4414_GP,  2.5490_gp&
      &,  2.7623_GP,  3.0991_GP,  3.3252_GP,  3.1836_GP,  3.2428_GP,  3.0259_GP,  3.1225_gp&
      &,  2.9032_GP,  2.7621_GP,  3.2490_GP,  3.5110_GP,  3.4429_GP,  3.3845_GP,  3.3574_gp&
      &,  3.6045_GP,  3.3658_GP,  3.2013_GP,  3.6110_GP,  3.8241_GP,  3.7090_GP,  3.6496_gp&
      &,  3.6333_GP,  3.0896_GP,  3.5462_GP,  2.4926_GP,  2.7136_GP,  3.0693_GP,  3.2699_gp&
      &,  3.1272_GP,  3.1893_GP,  2.9658_GP,  3.0972_GP,  2.8778_GP,  2.7358_GP,  3.2206_gp&
      &,  3.4566_GP,  3.3896_GP,  3.3257_GP,  3.2946_GP,  3.5693_GP,  3.3312_GP,  3.1670_gp&
      &,  3.5805_GP,  3.7711_GP,  3.6536_GP,  3.5927_GP,  3.5775_GP,  3.0411_GP,  3.4885_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 351: 420)=(/&                                            
         3.4421_GP,  2.4667_GP,  2.6709_GP,  3.0575_GP,  3.2357_GP,  3.0908_GP,   3.1537_gp&
      &,  2.9235_GP,  3.0669_GP,  2.8476_GP,  2.7054_GP,  3.2064_GP,  3.4519_GP,  3.3593_gp&
      &,  3.2921_GP,  3.2577_GP,  3.2161_GP,  3.2982_GP,  3.1339_GP,  3.5606_GP,  3.7582_gp&
      &,  3.6432_GP,  3.5833_GP,  3.5691_GP,  3.0161_GP,  3.4812_GP,  3.4339_GP,  3.4327_gp&
      &,  2.4515_GP,  2.6338_GP,  3.0511_GP,  3.2229_GP,  3.0630_GP,  3.1265_GP,  2.8909_gp&
      &,  3.0253_GP,  2.8184_GP,  2.6764_GP,  3.1968_GP,  3.4114_GP,  3.3492_GP,  3.2691_gp&
      &,  3.2320_GP,  3.1786_GP,  3.2680_GP,  3.1036_GP,  3.5453_GP,  3.7259_GP,  3.6090_gp&
      &,  3.5473_GP,  3.5327_GP,  3.0018_GP,  3.4413_GP,  3.3907_GP,  3.3593_GP,  3.3462_gp&
      &,  2.4413_GP,  2.6006_GP,  3.0540_GP,  3.1987_GP,  3.0490_GP,  3.1058_GP,  2.8643_gp&
      &,  2.9948_GP,  2.7908_GP,  2.6491_GP,  3.1950_GP,  3.3922_GP,  3.3316_GP,  3.2585_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 421: 490)=(/&                                            
         3.2136_GP,  3.1516_GP,  3.2364_GP,  3.0752_GP,  3.5368_GP,  3.7117_GP,  3.5941_gp&
      &,  3.5313_GP,  3.5164_GP,  2.9962_GP,  3.4225_GP,  3.3699_GP,  3.3370_GP,  3.3234_gp&
      &,  3.3008_GP,  2.4318_GP,  2.5729_GP,  3.0416_GP,  3.1639_GP,  3.0196_GP,  3.0843_gp&
      &,  2.8413_GP,  2.7436_GP,  2.7608_GP,  2.6271_GP,  3.1811_GP,  3.3591_GP,  3.3045_gp&
      &,  3.2349_GP,  3.1942_GP,  3.1291_GP,  3.2111_GP,  3.0534_GP,  3.5189_GP,  3.6809_gp&
      &,  3.5635_GP,  3.5001_GP,  3.4854_GP,  2.9857_GP,  3.3897_GP,  3.3363_GP,  3.3027_gp&
      &,  3.2890_GP,  3.2655_GP,  3.2309_GP,  2.8502_GP,  2.6934_GP,  3.2467_GP,  3.1921_gp&
      &,  3.5663_GP,  3.2541_GP,  3.0571_GP,  2.9048_GP,  2.8657_GP,  2.7438_GP,  3.3547_gp&
      &,  3.3510_GP,  3.9837_GP,  3.6871_GP,  3.4862_GP,  3.3389_GP,  3.2413_GP,  3.1708_gp&
      &,  3.6096_GP,  3.6280_GP,  3.6860_GP,  3.5568_GP,  3.4836_GP,  3.2868_GP,  3.3994_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 491: 560)=(/&                                            
         3.3476_GP,  3.3170_GP,  3.2950_GP,  3.2874_GP,  3.2606_GP,  3.9579_GP,  2.9226_gp&
      &,  2.6838_GP,  3.7867_GP,  3.1732_GP,  3.3872_GP,  3.3643_GP,  3.1267_GP,  2.9541_gp&
      &,  2.8505_GP,  2.7781_GP,  3.8475_GP,  3.3336_GP,  3.7359_GP,  3.8266_GP,  3.5733_gp&
      &,  3.3959_GP,  3.2775_GP,  3.1915_GP,  3.9878_GP,  3.8816_GP,  3.5810_GP,  3.5364_gp&
      &,  3.5060_GP,  3.8097_GP,  3.3925_GP,  3.3348_GP,  3.3019_GP,  3.2796_GP,  3.2662_gp&
      &,  3.2464_GP,  3.7136_GP,  3.8619_GP,  2.9140_GP,  2.6271_GP,  3.4771_GP,  3.1774_gp&
      &,  3.2560_GP,  3.1970_GP,  3.1207_GP,  2.9406_GP,  2.8322_GP,  2.7571_GP,  3.5455_gp&
      &,  3.3514_GP,  3.5837_GP,  3.6177_GP,  3.5816_GP,  3.3902_GP,  3.2604_GP,  3.1652_gp&
      &,  3.7037_GP,  3.6283_GP,  3.5858_GP,  3.5330_GP,  3.4884_GP,  3.5789_GP,  3.4094_gp&
      &,  3.3473_GP,  3.3118_GP,  3.2876_GP,  3.2707_GP,  3.2521_GP,  3.5570_GP,  3.6496_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 561: 630)=(/&                                            
         3.6625_GP,  2.7300_GP,  2.5870_GP,  3.2471_GP,  3.1487_GP,  3.1667_GP,  3.0914_gp&
      &,  3.0107_GP,  2.9812_GP,  2.8300_GP,  2.7284_GP,  3.3259_GP,  3.3182_GP,  3.4707_gp&
      &,  3.4748_GP,  3.4279_GP,  3.4182_GP,  3.2547_GP,  3.1353_GP,  3.5116_GP,  3.9432_gp&
      &,  3.8828_GP,  3.8303_GP,  3.7880_GP,  3.3760_GP,  3.7218_GP,  3.3408_GP,  3.3059_gp&
      &,  3.2698_GP,  3.2446_GP,  3.2229_GP,  3.4422_GP,  3.5023_GP,  3.5009_GP,  3.5268_gp&
      &,  2.6026_GP,  2.5355_GP,  3.1129_GP,  3.2863_GP,  3.1029_GP,  3.0108_GP,  2.9227_gp&
      &,  2.8694_GP,  2.8109_GP,  2.6929_GP,  3.1958_GP,  3.4670_GP,  3.4018_GP,  3.3805_gp&
      &,  3.3218_GP,  3.2815_GP,  3.2346_GP,  3.0994_GP,  3.3937_GP,  3.7266_GP,  3.6697_gp&
      &,  3.6164_GP,  3.5730_GP,  3.2522_GP,  3.5051_GP,  3.4686_GP,  3.4355_GP,  3.4084_gp&
      &,  3.3748_GP,  3.3496_GP,  3.3692_GP,  3.4052_GP,  3.3910_GP,  3.3849_GP,  3.3662_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 631: 700)=(/&                                            
         2.5087_GP,  2.4814_GP,  3.0239_GP,  3.1312_GP,  3.0535_GP,  2.9457_GP,  2.8496_gp&
      &,  2.7780_GP,  2.7828_GP,  2.6532_GP,  3.1063_GP,  3.3143_GP,  3.3549_GP,  3.3120_gp&
      &,  3.2421_GP,  3.1787_GP,  3.1176_GP,  3.0613_GP,  3.3082_GP,  3.5755_GP,  3.5222_gp&
      &,  3.4678_GP,  3.4231_GP,  3.1684_GP,  3.3528_GP,  3.3162_GP,  3.2827_GP,  3.2527_gp&
      &,  3.2308_GP,  3.2029_GP,  3.3173_GP,  3.3343_GP,  3.3092_GP,  3.2795_GP,  3.2452_gp&
      &,  3.2096_GP,  3.2893_GP,  2.8991_GP,  4.0388_GP,  3.6100_GP,  3.9388_GP,  3.4475_gp&
      &,  3.1590_GP,  2.9812_GP,  2.8586_GP,  2.7683_GP,  4.1428_GP,  3.7911_GP,  3.8225_gp&
      &,  4.0372_GP,  3.7059_GP,  3.4935_GP,  3.3529_GP,  3.2492_GP,  4.4352_GP,  4.0826_gp&
      &,  3.9733_GP,  3.9254_GP,  3.8646_GP,  3.9315_GP,  3.7837_GP,  3.7465_GP,  3.7211_gp&
      &,  3.7012_GP,  3.6893_GP,  3.6676_GP,  3.7736_GP,  4.0660_GP,  3.7926_GP,  3.6158_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 701: 770)=(/&                                            
         3.5017_GP,  3.4166_GP,  4.6176_GP,  2.8786_GP,  3.1658_GP,  3.5823_GP,  3.7689_gp&
      &,  3.5762_GP,  3.5789_GP,  3.3552_GP,  3.4004_GP,  3.1722_GP,  3.0212_GP,  3.7241_gp&
      &,  3.9604_GP,  3.8500_GP,  3.9844_GP,  3.7035_GP,  3.9161_GP,  3.6751_GP,  3.5075_gp&
      &,  4.1151_GP,  4.2877_GP,  4.1579_GP,  4.1247_GP,  4.0617_GP,  3.4874_GP,  3.9848_gp&
      &,  3.9280_GP,  3.9079_GP,  3.8751_GP,  3.8604_GP,  3.8277_GP,  3.8002_GP,  3.9981_gp&
      &,  3.7544_GP,  4.0371_GP,  3.8225_GP,  3.6718_GP,  4.3092_GP,  4.4764_GP,  2.8997_gp&
      &,  3.0953_GP,  3.4524_GP,  3.6107_GP,  3.6062_GP,  3.5783_GP,  3.3463_GP,  3.3855_gp&
      &,  3.1746_GP,  3.0381_GP,  3.6019_GP,  3.7938_GP,  3.8697_GP,  3.9781_GP,  3.6877_gp&
      &,  3.8736_GP,  3.6451_GP,  3.4890_GP,  3.9858_GP,  4.1179_GP,  4.0430_GP,  3.9563_gp&
      &,  3.9182_GP,  3.4002_GP,  3.8310_GP,  3.7716_GP,  3.7543_GP,  3.7203_GP,  3.7053_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 771: 840)=(/&                                            
         3.6742_GP,  3.8318_GP,  3.7631_GP,  3.7392_GP,  3.9892_GP,  3.7832_GP,  3.6406_gp&
      &,  4.1701_GP,  4.3016_GP,  4.2196_GP,  2.8535_GP,  3.0167_GP,  3.3978_GP,  3.5363_gp&
      &,  3.5393_GP,  3.5301_GP,  3.2960_GP,  3.3352_GP,  3.1287_GP,  2.9967_GP,  3.6659_gp&
      &,  3.7239_GP,  3.8070_GP,  3.7165_GP,  3.6368_GP,  3.8162_GP,  3.5885_GP,  3.4336_gp&
      &,  3.9829_GP,  4.0529_GP,  3.9584_GP,  3.9025_GP,  3.8607_GP,  3.3673_GP,  3.7658_gp&
      &,  3.7035_GP,  3.6866_GP,  3.6504_GP,  3.6339_GP,  3.6024_GP,  3.7708_GP,  3.7283_gp&
      &,  3.6896_GP,  3.9315_GP,  3.7250_GP,  3.5819_GP,  4.1457_GP,  4.2280_GP,  4.1130_gp&
      &,  4.0597_GP,  3.0905_GP,  2.7998_GP,  3.6448_GP,  3.0739_GP,  3.2996_GP,  3.5262_gp&
      &,  3.2559_GP,  3.0518_GP,  2.9394_GP,  2.8658_GP,  3.7514_GP,  3.2295_GP,  3.5643_gp&
      &,  3.7808_GP,  3.6931_GP,  3.4723_GP,  3.3357_GP,  3.2429_GP,  4.0280_GP,  3.5589_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 841: 910)=(/&                                            
         3.4636_GP,  3.4994_GP,  3.4309_GP,  3.6177_GP,  3.2946_GP,  3.2376_GP,  3.2050_gp&
      &,  3.1847_GP,  3.1715_GP,  3.1599_GP,  3.5555_GP,  3.8111_GP,  3.7693_GP,  3.5718_gp&
      &,  3.4498_GP,  3.3662_GP,  4.1608_GP,  3.7417_GP,  3.6536_GP,  3.6154_GP,  3.8596_gp&
      &,  3.0301_GP,  2.7312_GP,  3.5821_GP,  3.0473_GP,  3.2137_GP,  3.4679_GP,  3.1975_gp&
      &,  2.9969_GP,  2.8847_GP,  2.8110_GP,  3.6931_GP,  3.2076_GP,  3.4943_GP,  3.5956_gp&
      &,  3.6379_GP,  3.4190_GP,  3.2808_GP,  3.1860_GP,  3.9850_GP,  3.5105_GP,  3.4330_gp&
      &,  3.3797_GP,  3.4155_GP,  3.6033_GP,  3.2737_GP,  3.2145_GP,  3.1807_GP,  3.1596_gp&
      &,  3.1461_GP,  3.1337_GP,  3.4812_GP,  3.6251_GP,  3.7152_GP,  3.5201_GP,  3.3966_gp&
      &,  3.3107_GP,  4.1128_GP,  3.6899_GP,  3.6082_GP,  3.5604_GP,  3.7834_GP,  3.7543_gp&
      &,  2.9189_GP,  2.6777_GP,  3.4925_GP,  2.9648_GP,  3.1216_GP,  3.2940_GP,  3.0975_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 911: 980)=(/&                                            
         2.9757_GP,  2.8493_GP,  2.7638_GP,  3.6085_GP,  3.1214_GP,  3.4006_GP,   3.4793_gp&
      &,  3.5147_GP,  3.3806_GP,  3.2356_GP,  3.1335_GP,  3.9144_GP,  3.4183_GP,  3.3369_gp&
      &,  3.2803_GP,  3.2679_GP,  3.4871_GP,  3.1714_GP,  3.1521_GP,  3.1101_GP,  3.0843_gp&
      &,  3.0670_GP,  3.0539_GP,  3.3890_GP,  3.5086_GP,  3.5895_GP,  3.4783_GP,  3.3484_gp&
      &,  3.2559_GP,  4.0422_GP,  3.5967_GP,  3.5113_GP,  3.4576_GP,  3.6594_GP,  3.6313_gp&
      &,  3.5690_GP,  2.8578_GP,  2.6334_GP,  3.4673_GP,  2.9245_GP,  3.0732_GP,  3.2435_gp&
      &,  3.0338_GP,  2.9462_GP,  2.8143_GP,  2.7240_GP,  3.5832_GP,  3.0789_GP,  3.3617_gp&
      &,  3.4246_GP,  3.4505_GP,  3.3443_GP,  3.1964_GP,  3.0913_GP,  3.8921_GP,  3.3713_gp&
      &,  3.2873_GP,  3.2281_GP,  3.2165_GP,  3.4386_GP,  3.1164_GP,  3.1220_GP,  3.0761_gp&
      &,  3.0480_GP,  3.0295_GP,  3.0155_GP,  3.3495_GP,  3.4543_GP,  3.5260_GP,  3.4413_gp&
      &/)                                                                              
       r0ab_ref_table_grimme( 981:1050)=(/&                                            
         3.3085_GP,  3.2134_GP,  4.0170_GP,  3.5464_GP,  3.4587_GP,  3.4006_GP,   3.6027_gp&
      &,  3.5730_GP,  3.4945_GP,  3.4623_GP,  2.8240_GP,  2.5960_GP,  3.4635_GP,  2.9032_gp&
      &,  3.0431_GP,  3.2115_GP,  2.9892_GP,  2.9148_GP,  2.7801_GP,  2.6873_GP,  3.5776_gp&
      &,  3.0568_GP,  3.3433_GP,  3.3949_GP,  3.4132_GP,  3.3116_GP,  3.1616_GP,  3.0548_gp&
      &,  3.8859_GP,  3.3719_GP,  3.2917_GP,  3.2345_GP,  3.2274_GP,  3.4171_GP,  3.1293_gp&
      &,  3.0567_GP,  3.0565_GP,  3.0274_GP,  3.0087_GP,  2.9939_GP,  3.3293_GP,  3.4249_gp&
      &,  3.4902_GP,  3.4091_GP,  3.2744_GP,  3.1776_GP,  4.0078_GP,  3.5374_GP,  3.4537_gp&
      &,  3.3956_GP,  3.5747_GP,  3.5430_GP,  3.4522_GP,  3.4160_GP,  3.3975_GP,  2.8004_gp&
      &,  2.5621_GP,  3.4617_GP,  2.9154_GP,  3.0203_GP,  3.1875_GP,  2.9548_GP,  2.8038_gp&
      &,  2.7472_GP,  2.6530_GP,  3.5736_GP,  3.0584_GP,  3.3304_GP,  3.3748_GP,  3.3871_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1051:1120)=(/&                                            
         3.2028_GP,  3.1296_GP,  3.0214_GP,  3.8796_GP,  3.3337_GP,  3.2492_GP,   3.1883_gp&
      &,  3.1802_GP,  3.4050_GP,  3.0756_GP,  3.0478_GP,  3.0322_GP,  3.0323_GP,  3.0163_gp&
      &,  3.0019_GP,  3.3145_GP,  3.4050_GP,  3.4656_GP,  3.3021_GP,  3.2433_GP,  3.1453_gp&
      &,  3.9991_GP,  3.5017_GP,  3.4141_GP,  3.3520_GP,  3.5583_GP,  3.5251_GP,  3.4243_gp&
      &,  3.3851_GP,  3.3662_GP,  3.3525_GP,  2.7846_GP,  2.5324_GP,  3.4652_GP,  2.8759_gp&
      &,  3.0051_GP,  3.1692_GP,  2.9273_GP,  2.7615_GP,  2.7164_GP,  2.6212_GP,  3.5744_gp&
      &,  3.0275_GP,  3.3249_GP,  3.3627_GP,  3.3686_GP,  3.1669_GP,  3.0584_GP,  2.9915_gp&
      &,  3.8773_GP,  3.3099_GP,  3.2231_GP,  3.1600_GP,  3.1520_GP,  3.4023_GP,  3.0426_gp&
      &,  3.0099_GP,  2.9920_GP,  2.9809_GP,  2.9800_GP,  2.9646_GP,  3.3068_GP,  3.3930_gp&
      &,  3.4486_GP,  3.2682_GP,  3.1729_GP,  3.1168_GP,  3.9952_GP,  3.4796_GP,  3.3901_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1121:1190)=(/&                                            
         3.3255_GP,  3.5530_GP,  3.5183_GP,  3.4097_GP,  3.3683_GP,  3.3492_GP,   3.3360_gp&
      &,  3.3308_GP,  2.5424_GP,  2.6601_GP,  3.2555_GP,  3.2807_GP,  3.1384_GP,  3.1737_gp&
      &,  2.9397_GP,  2.8429_GP,  2.8492_GP,  2.7225_GP,  3.3875_GP,  3.4910_GP,  3.4520_gp&
      &,  3.3608_GP,  3.3036_GP,  3.2345_GP,  3.2999_GP,  3.1487_GP,  3.7409_GP,  3.8392_gp&
      &,  3.7148_GP,  3.6439_GP,  3.6182_GP,  3.1753_GP,  3.5210_GP,  3.4639_GP,  3.4265_gp&
      &,  3.4075_GP,  3.3828_GP,  3.3474_GP,  3.4071_GP,  3.3754_GP,  3.3646_GP,  3.3308_gp&
      &,  3.4393_GP,  3.2993_GP,  3.8768_GP,  3.9891_GP,  3.8310_GP,  3.7483_GP,  3.3417_gp&
      &,  3.3019_GP,  3.2250_GP,  3.1832_GP,  3.1578_GP,  3.1564_GP,  3.1224_GP,  3.4620_gp&
      &,  2.9743_GP,  2.8058_GP,  3.4830_GP,  3.3474_GP,  3.6863_GP,  3.3617_GP,  3.1608_gp&
      &,  3.0069_GP,  2.9640_GP,  2.8427_GP,  3.5885_GP,  3.5219_GP,  4.1314_GP,  3.8120_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1191:1260)=(/&                                            
         3.6015_GP,  3.4502_GP,  3.3498_GP,  3.2777_GP,  3.8635_GP,  3.8232_GP,   3.8486_gp&
      &,  3.7215_GP,  3.6487_GP,  3.4724_GP,  3.5627_GP,  3.5087_GP,  3.4757_GP,  3.4517_gp&
      &,  3.4423_GP,  3.4139_GP,  4.1028_GP,  3.8388_GP,  3.6745_GP,  3.5562_GP,  3.4806_gp&
      &,  3.4272_GP,  4.0182_GP,  3.9991_GP,  4.0007_GP,  3.9282_GP,  3.7238_GP,  3.6498_gp&
      &,  3.5605_GP,  3.5211_GP,  3.5009_GP,  3.4859_GP,  3.4785_GP,  3.5621_GP,  4.2623_gp&
      &,  3.0775_GP,  2.8275_GP,  4.0181_GP,  3.3385_GP,  3.5379_GP,  3.5036_GP,  3.2589_gp&
      &,  3.0804_GP,  3.0094_GP,  2.9003_GP,  4.0869_GP,  3.5088_GP,  3.9105_GP,  3.9833_gp&
      &,  3.7176_GP,  3.5323_GP,  3.4102_GP,  3.3227_GP,  4.2702_GP,  4.0888_GP,  3.7560_gp&
      &,  3.7687_GP,  3.6681_GP,  3.6405_GP,  3.5569_GP,  3.4990_GP,  3.4659_GP,  3.4433_gp&
      &,  3.4330_GP,  3.4092_GP,  3.8867_GP,  4.0190_GP,  3.7961_GP,  3.6412_GP,  3.5405_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1261:1330)=(/&                                            
         3.4681_GP,  4.3538_GP,  4.2136_GP,  3.9381_GP,  3.8912_GP,  3.9681_GP,   3.7909_gp&
      &,  3.6774_GP,  3.6262_GP,  3.5999_GP,  3.5823_GP,  3.5727_GP,  3.5419_GP,  4.0245_gp&
      &,  4.1874_GP,  3.0893_GP,  2.7917_GP,  3.7262_GP,  3.3518_GP,  3.4241_GP,  3.5433_gp&
      &,  3.2773_GP,  3.0890_GP,  2.9775_GP,  2.9010_GP,  3.8048_GP,  3.5362_GP,  3.7746_gp&
      &,  3.7911_GP,  3.7511_GP,  3.5495_GP,  3.4149_GP,  3.3177_GP,  4.0129_GP,  3.8370_gp&
      &,  3.7739_GP,  3.7125_GP,  3.7152_GP,  3.7701_GP,  3.5813_GP,  3.5187_GP,  3.4835_gp&
      &,  3.4595_GP,  3.4439_GP,  3.4242_GP,  3.7476_GP,  3.8239_GP,  3.8346_GP,  3.6627_gp&
      &,  3.5479_GP,  3.4639_GP,  4.1026_GP,  3.9733_GP,  3.9292_GP,  3.8667_GP,  3.9513_gp&
      &,  3.8959_GP,  3.7698_GP,  3.7089_GP,  3.6765_GP,  3.6548_GP,  3.6409_GP,  3.5398_gp&
      &,  3.8759_GP,  3.9804_GP,  4.0150_GP,  2.9091_GP,  2.7638_GP,  3.5066_GP,  3.3377_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1331:1400)=(/&                                            
         3.3481_GP,  3.2633_GP,  3.1810_GP,  3.1428_GP,  2.9872_GP,  2.8837_GP,   3.5929_gp&
      &,  3.5183_GP,  3.6729_GP,  3.6596_GP,  3.6082_GP,  3.5927_GP,  3.4224_GP,  3.2997_gp&
      &,  3.8190_GP,  4.1865_GP,  4.1114_GP,  4.0540_GP,  3.6325_GP,  3.5697_GP,  3.5561_gp&
      &,  3.5259_GP,  3.4901_GP,  3.4552_GP,  3.4315_GP,  3.4091_GP,  3.6438_GP,  3.6879_gp&
      &,  3.6832_GP,  3.7043_GP,  3.5557_GP,  3.4466_GP,  3.9203_GP,  4.2919_GP,  4.2196_gp&
      &,  4.1542_GP,  3.7573_GP,  3.7039_GP,  3.6546_GP,  3.6151_GP,  3.5293_GP,  3.4849_gp&
      &,  3.4552_GP,  3.5192_GP,  3.7673_GP,  3.8359_GP,  3.8525_GP,  3.8901_GP,  2.7806_gp&
      &,  2.7209_GP,  3.3812_GP,  3.4958_GP,  3.2913_GP,  3.1888_GP,  3.0990_GP,  3.0394_gp&
      &,  2.9789_GP,  2.8582_GP,  3.4716_GP,  3.6883_GP,  3.6105_GP,  3.5704_GP,  3.5059_gp&
      &,  3.4619_GP,  3.4138_GP,  3.2742_GP,  3.7080_GP,  3.9773_GP,  3.9010_GP,  3.8409_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1401:1470)=(/&                                            
         3.7944_GP,  3.4465_GP,  3.7235_GP,  3.6808_GP,  3.6453_GP,  3.6168_GP,   3.5844_gp&
      &,  3.5576_GP,  3.5772_GP,  3.5959_GP,  3.5768_GP,  3.5678_GP,  3.5486_GP,  3.4228_gp&
      &,  3.8107_GP,  4.0866_GP,  4.0169_GP,  3.9476_GP,  3.6358_GP,  3.5800_GP,  3.5260_gp&
      &,  3.4838_GP,  3.4501_GP,  3.4204_GP,  3.3553_GP,  3.6487_GP,  3.6973_GP,  3.7398_gp&
      &,  3.7405_GP,  3.7459_GP,  3.7380_GP,  2.6848_GP,  2.6740_GP,  3.2925_GP,  3.3386_gp&
      &,  3.2473_GP,  3.1284_GP,  3.0301_GP,  2.9531_GP,  2.9602_GP,  2.8272_GP,  3.3830_gp&
      &,  3.5358_GP,  3.5672_GP,  3.5049_GP,  3.4284_GP,  3.3621_GP,  3.3001_GP,  3.2451_gp&
      &,  3.6209_GP,  3.8299_GP,  3.7543_GP,  3.6920_GP,  3.6436_GP,  3.3598_GP,  3.5701_gp&
      &,  3.5266_GP,  3.4904_GP,  3.4590_GP,  3.4364_GP,  3.4077_GP,  3.5287_GP,  3.5280_gp&
      &,  3.4969_GP,  3.4650_GP,  3.4304_GP,  3.3963_GP,  3.7229_GP,  3.9402_GP,  3.8753_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1471:1540)=(/&                                            
         3.8035_GP,  3.5499_GP,  3.4913_GP,  3.4319_GP,  3.3873_GP,  3.3520_GP,  3.3209_gp&
      &,  3.2948_GP,  3.5052_GP,  3.6465_GP,  3.6696_GP,  3.6577_GP,  3.6388_GP,  3.6142_gp&
      &,  3.5889_GP,  3.3968_GP,  3.0122_GP,  4.2241_GP,  3.7887_GP,  4.0049_GP,  3.5384_gp&
      &,  3.2698_GP,  3.1083_GP,  2.9917_GP,  2.9057_GP,  4.3340_GP,  3.9900_GP,  4.6588_gp&
      &,  4.1278_GP,  3.8125_GP,  3.6189_GP,  3.4851_GP,  3.3859_GP,  4.6531_GP,  4.3134_gp&
      &,  4.2258_GP,  4.1309_GP,  4.0692_GP,  4.0944_GP,  3.9850_GP,  3.9416_GP,  3.9112_gp&
      &,  3.8873_GP,  3.8736_GP,  3.8473_GP,  4.6027_GP,  4.1538_GP,  3.8994_GP,  3.7419_gp&
      &,  3.6356_GP,  3.5548_GP,  4.8353_GP,  4.5413_GP,  4.3891_GP,  4.3416_GP,  4.3243_gp&
      &,  4.2753_GP,  4.2053_GP,  4.1790_GP,  4.1685_GP,  4.1585_GP,  4.1536_GP,  4.0579_gp&
      &,  4.1980_GP,  4.4564_GP,  4.2192_GP,  4.0528_GP,  3.9489_GP,  3.8642_GP,  5.0567_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1541:1610)=(/&                                            
         3.0630_GP,  3.3271_GP,  4.0432_GP,  4.0046_GP,  4.1555_GP,  3.7426_GP,  3.5130_gp&
      &,  3.5174_GP,  3.2884_GP,  3.1378_GP,  4.1894_GP,  4.2321_GP,  4.1725_GP,  4.1833_gp&
      &,  3.8929_GP,  4.0544_GP,  3.8118_GP,  3.6414_GP,  4.6373_GP,  4.6268_GP,  4.4750_gp&
      &,  4.4134_GP,  4.3458_GP,  3.8582_GP,  4.2583_GP,  4.1898_GP,  4.1562_GP,  4.1191_gp&
      &,  4.1069_GP,  4.0639_GP,  4.1257_GP,  4.1974_GP,  3.9532_GP,  4.1794_GP,  3.9660_gp&
      &,  3.8130_GP,  4.8160_GP,  4.8272_GP,  4.6294_GP,  4.5840_GP,  4.0770_GP,  4.0088_gp&
      &,  3.9103_GP,  3.8536_GP,  3.8324_GP,  3.7995_GP,  3.7826_GP,  4.2294_GP,  4.3380_gp&
      &,  4.4352_GP,  4.1933_GP,  4.4580_GP,  4.2554_GP,  4.1072_GP,  5.0454_GP,  5.1814_gp&
      &,  3.0632_GP,  3.2662_GP,  3.6432_GP,  3.8088_GP,  3.7910_GP,  3.7381_GP,  3.5093_gp&
      &,  3.5155_GP,  3.3047_GP,  3.1681_GP,  3.7871_GP,  3.9924_GP,  4.0637_GP,  4.1382_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1611:1680)=(/&                                            
         3.8591_GP,  4.0164_GP,  3.7878_GP,  3.6316_GP,  4.1741_GP,  4.3166_GP,  4.2395_gp&
      &,  4.1831_GP,  4.1107_GP,  3.5857_GP,  4.0270_GP,  3.9676_GP,  3.9463_GP,  3.9150_gp&
      &,  3.9021_GP,  3.8708_GP,  4.0240_GP,  4.1551_GP,  3.9108_GP,  4.1337_GP,  3.9289_gp&
      &,  3.7873_GP,  4.3666_GP,  4.5080_GP,  4.4232_GP,  4.3155_GP,  3.8461_GP,  3.8007_gp&
      &,  3.6991_GP,  3.6447_GP,  3.6308_GP,  3.5959_GP,  3.5749_GP,  4.0359_GP,  4.3124_gp&
      &,  4.3539_GP,  4.1122_GP,  4.3772_GP,  4.1785_GP,  4.0386_GP,  4.7004_GP,  4.8604_gp&
      &,  4.6261_GP,  2.9455_GP,  3.2470_GP,  3.6108_GP,  3.8522_GP,  3.6625_GP,  3.6598_gp&
      &,  3.4411_GP,  3.4660_GP,  3.2415_GP,  3.0944_GP,  3.7514_GP,  4.0397_GP,  3.9231_gp&
      &,  4.0561_GP,  3.7860_GP,  3.9845_GP,  3.7454_GP,  3.5802_GP,  4.1366_GP,  4.3581_gp&
      &,  4.2351_GP,  4.2011_GP,  4.1402_GP,  3.5381_GP,  4.0653_GP,  4.0093_GP,  3.9883_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1681:1750)=(/&                                            
         3.9570_GP,  3.9429_GP,  3.9112_GP,  3.8728_GP,  4.0682_GP,  3.8351_GP,  4.1054_gp&
      &,  3.8928_GP,  3.7445_GP,  4.3415_GP,  4.5497_GP,  4.3833_GP,  4.3122_GP,  3.8051_gp&
      &,  3.7583_GP,  3.6622_GP,  3.6108_GP,  3.5971_GP,  3.5628_GP,  3.5408_GP,  4.0780_gp&
      &,  4.0727_GP,  4.2836_GP,  4.0553_GP,  4.3647_GP,  4.1622_GP,  4.0178_GP,  4.5802_gp&
      &,  4.9125_GP,  4.5861_GP,  4.6201_GP,  2.9244_GP,  3.2241_GP,  3.5848_GP,  3.8293_gp&
      &,  3.6395_GP,  3.6400_GP,  3.4204_GP,  3.4499_GP,  3.2253_GP,  3.0779_GP,  3.7257_gp&
      &,  4.0170_GP,  3.9003_GP,  4.0372_GP,  3.7653_GP,  3.9672_GP,  3.7283_GP,  3.5630_gp&
      &,  4.1092_GP,  4.3347_GP,  4.2117_GP,  4.1793_GP,  4.1179_GP,  3.5139_GP,  4.0426_gp&
      &,  3.9867_GP,  3.9661_GP,  3.9345_GP,  3.9200_GP,  3.8883_GP,  3.8498_GP,  4.0496_gp&
      &,  3.8145_GP,  4.0881_GP,  3.8756_GP,  3.7271_GP,  4.3128_GP,  4.5242_GP,  4.3578_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1751:1820)=(/&                                            
         4.2870_GP,  3.7796_GP,  3.7318_GP,  3.6364_GP,  3.5854_GP,  3.5726_GP,  3.5378_gp&
      &,  3.5155_GP,  4.0527_GP,  4.0478_GP,  4.2630_GP,  4.0322_GP,  4.3449_GP,  4.1421_gp&
      &,  3.9975_GP,  4.5499_GP,  4.8825_GP,  4.5601_GP,  4.5950_GP,  4.5702_GP,  2.9046_gp&
      &,  3.2044_GP,  3.5621_GP,  3.8078_GP,  3.6185_GP,  3.6220_GP,  3.4019_GP,  3.4359_gp&
      &,  3.2110_GP,  3.0635_GP,  3.7037_GP,  3.9958_GP,  3.8792_GP,  4.0194_GP,  3.7460_gp&
      &,  3.9517_GP,  3.7128_GP,  3.5474_GP,  4.0872_GP,  4.3138_GP,  4.1906_GP,  4.1593_gp&
      &,  4.0973_GP,  3.4919_GP,  4.0216_GP,  3.9657_GP,  3.9454_GP,  3.9134_GP,  3.8986_gp&
      &,  3.8669_GP,  3.8289_GP,  4.0323_GP,  3.7954_GP,  4.0725_GP,  3.8598_GP,  3.7113_gp&
      &,  4.2896_GP,  4.5021_GP,  4.3325_GP,  4.2645_GP,  3.7571_GP,  3.7083_GP,  3.6136_gp&
      &,  3.5628_GP,  3.5507_GP,  3.5155_GP,  3.4929_GP,  4.0297_GP,  4.0234_GP,  4.2442_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1821:1890)=(/&                                            
         4.0112_GP,  4.3274_GP,  4.1240_GP,  3.9793_GP,  4.5257_GP,  4.8568_GP,  4.5353_gp&
      &,  4.5733_GP,  4.5485_GP,  4.5271_GP,  2.8878_GP,  3.1890_GP,  3.5412_GP,  3.7908_gp&
      &,  3.5974_GP,  3.6078_GP,  3.3871_GP,  3.4243_GP,  3.1992_GP,  3.0513_GP,  3.6831_gp&
      &,  3.9784_GP,  3.8579_GP,  4.0049_GP,  3.7304_GP,  3.9392_GP,  3.7002_GP,  3.5347_gp&
      &,  4.0657_GP,  4.2955_GP,  4.1705_GP,  4.1424_GP,  4.0800_GP,  3.4717_GP,  4.0043_gp&
      &,  3.9485_GP,  3.9286_GP,  3.8965_GP,  3.8815_GP,  3.8500_GP,  3.8073_GP,  4.0180_gp&
      &,  3.7796_GP,  4.0598_GP,  3.8470_GP,  3.6983_GP,  4.2678_GP,  4.4830_GP,  4.3132_gp&
      &,  4.2444_GP,  3.7370_GP,  3.6876_GP,  3.5935_GP,  3.5428_GP,  3.5314_GP,  3.4958_gp&
      &,  3.4730_GP,  4.0117_GP,  4.0043_GP,  4.2287_GP,  3.9939_GP,  4.3134_GP,  4.1096_gp&
      &,  3.9646_GP,  4.5032_GP,  4.8356_GP,  4.5156_GP,  4.5544_GP,  4.5297_GP,  4.5083_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1891:1960)=(/&                                            
         4.4896_GP,  2.8709_GP,  3.1737_GP,  3.5199_GP,  3.7734_GP,  3.5802_GP,  3.5934_gp&
      &,  3.3724_GP,  3.4128_GP,  3.1877_GP,  3.0396_GP,  3.6624_GP,  3.9608_GP,  3.8397_gp&
      &,  3.9893_GP,  3.7145_GP,  3.9266_GP,  3.6877_GP,  3.5222_GP,  4.0448_GP,  4.2771_gp&
      &,  4.1523_GP,  4.1247_GP,  4.0626_GP,  3.4530_GP,  3.9866_GP,  3.9310_GP,  3.9115_gp&
      &,  3.8792_GP,  3.8641_GP,  3.8326_GP,  3.7892_GP,  4.0025_GP,  3.7636_GP,  4.0471_gp&
      &,  3.8343_GP,  3.6854_GP,  4.2464_GP,  4.4635_GP,  4.2939_GP,  4.2252_GP,  3.7169_gp&
      &,  3.6675_GP,  3.5739_GP,  3.5235_GP,  3.5126_GP,  3.4768_GP,  3.4537_GP,  3.9932_gp&
      &,  3.9854_GP,  4.2123_GP,  3.9765_GP,  4.2992_GP,  4.0951_GP,  3.9500_GP,  4.4811_gp&
      &,  4.8135_GP,  4.4959_GP,  4.5351_GP,  4.5105_GP,  4.4891_GP,  4.4705_GP,  4.4515_gp&
      &,  2.8568_GP,  3.1608_GP,  3.5050_GP,  3.7598_GP,  3.5665_GP,  3.5803_GP,  3.3601_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(1961:2030)=(/&                                            
         3.4031_GP,  3.1779_GP,  3.0296_GP,  3.6479_GP,  3.9471_GP,  3.8262_GP,  3.9773_gp&
      &,  3.7015_GP,  3.9162_GP,  3.6771_GP,  3.5115_GP,  4.0306_GP,  4.2634_GP,  4.1385_gp&
      &,  4.1116_GP,  4.0489_GP,  3.4366_GP,  3.9732_GP,  3.9176_GP,  3.8983_GP,  3.8659_gp&
      &,  3.8507_GP,  3.8191_GP,  3.7757_GP,  3.9907_GP,  3.7506_GP,  4.0365_GP,  3.8235_gp&
      &,  3.6745_GP,  4.2314_GP,  4.4490_GP,  4.2792_GP,  4.2105_GP,  3.7003_GP,  3.6510_gp&
      &,  3.5578_GP,  3.5075_GP,  3.4971_GP,  3.4609_GP,  3.4377_GP,  3.9788_GP,  3.9712_gp&
      &,  4.1997_GP,  3.9624_GP,  4.2877_GP,  4.0831_GP,  3.9378_GP,  4.4655_GP,  4.7974_gp&
      &,  4.4813_GP,  4.5209_GP,  4.4964_GP,  4.4750_GP,  4.4565_GP,  4.4375_GP,  4.4234_gp&
      &,  2.6798_GP,  3.0151_GP,  3.2586_GP,  3.5292_GP,  3.5391_GP,  3.4902_GP,  3.2887_gp&
      &,  3.3322_GP,  3.1228_GP,  2.9888_GP,  3.4012_GP,  3.7145_GP,  3.7830_GP,  3.6665_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2031:2100)=(/&                                            
         3.5898_GP,  3.8077_GP,  3.5810_GP,  3.4265_GP,  3.7726_GP,  4.0307_GP,  3.9763_gp&
      &,  3.8890_GP,  3.8489_GP,  3.2706_GP,  3.7595_GP,  3.6984_GP,  3.6772_GP,  3.6428_gp&
      &,  3.6243_GP,  3.5951_GP,  3.7497_GP,  3.6775_GP,  3.6364_GP,  3.9203_GP,  3.7157_gp&
      &,  3.5746_GP,  3.9494_GP,  4.2076_GP,  4.1563_GP,  4.0508_GP,  3.5329_GP,  3.4780_gp&
      &,  3.3731_GP,  3.3126_GP,  3.2846_GP,  3.2426_GP,  3.2135_GP,  3.7491_GP,  3.9006_gp&
      &,  3.8332_GP,  3.8029_GP,  4.1436_GP,  3.9407_GP,  3.7998_GP,  4.1663_GP,  4.5309_gp&
      &,  4.3481_GP,  4.2911_GP,  4.2671_GP,  4.2415_GP,  4.2230_GP,  4.2047_GP,  4.1908_gp&
      &,  4.1243_GP,  2.5189_GP,  2.9703_GP,  3.3063_GP,  3.6235_GP,  3.4517_GP,  3.3989_gp&
      &,  3.2107_GP,  3.2434_GP,  3.0094_GP,  2.8580_GP,  3.4253_GP,  3.8157_GP,  3.7258_gp&
      &,  3.6132_GP,  3.5297_GP,  3.7566_GP,  3.5095_GP,  3.3368_GP,  3.7890_GP,  4.1298_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2101:2170)=(/&                                            
         4.0190_GP,  3.9573_GP,  3.9237_GP,  3.2677_GP,  3.8480_GP,  3.8157_GP,   3.7656_gp&
      &,  3.7317_GP,  3.7126_GP,  3.6814_GP,  3.6793_GP,  3.6218_GP,  3.5788_GP,  3.8763_gp&
      &,  3.6572_GP,  3.5022_GP,  3.9737_GP,  4.3255_GP,  4.1828_GP,  4.1158_GP,  3.5078_gp&
      &,  3.4595_GP,  3.3600_GP,  3.3088_GP,  3.2575_GP,  3.2164_GP,  3.1856_GP,  3.8522_gp&
      &,  3.8665_GP,  3.8075_GP,  3.7772_GP,  4.1391_GP,  3.9296_GP,  3.7772_GP,  4.2134_gp&
      &,  4.7308_GP,  4.3787_GP,  4.3894_GP,  4.3649_GP,  4.3441_GP,  4.3257_GP,  4.3073_gp&
      &,  4.2941_GP,  4.1252_GP,  4.2427_GP,  3.0481_GP,  2.9584_GP,  3.6919_GP,  3.5990_gp&
      &,  3.8881_GP,  3.4209_GP,  3.1606_GP,  3.1938_GP,  2.9975_GP,  2.8646_GP,  3.8138_gp&
      &,  3.7935_GP,  3.7081_GP,  3.9155_GP,  3.5910_GP,  3.4808_GP,  3.4886_GP,  3.3397_gp&
      &,  4.1336_GP,  4.1122_GP,  3.9888_GP,  3.9543_GP,  3.8917_GP,  3.5894_GP,  3.8131_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2171:2240)=(/&                                            
         3.7635_GP,  3.7419_GP,  3.7071_GP,  3.6880_GP,  3.6574_GP,  3.6546_GP,   3.9375_gp&
      &,  3.6579_GP,  3.5870_GP,  3.6361_GP,  3.5039_GP,  4.3149_GP,  4.2978_GP,  4.1321_gp&
      &,  4.1298_GP,  3.8164_GP,  3.7680_GP,  3.7154_GP,  3.6858_GP,  3.6709_GP,  3.6666_gp&
      &,  3.6517_GP,  3.8174_GP,  3.8608_GP,  4.1805_GP,  3.9102_GP,  3.8394_GP,  3.8968_gp&
      &,  3.7673_GP,  4.5274_GP,  4.6682_GP,  4.3344_GP,  4.3639_GP,  4.3384_GP,  4.3162_gp&
      &,  4.2972_GP,  4.2779_GP,  4.2636_GP,  4.0253_GP,  4.1168_GP,  4.1541_GP,  2.8136_gp&
      &,  3.0951_GP,  3.4635_GP,  3.6875_GP,  3.4987_GP,  3.5183_GP,  3.2937_GP,  3.3580_gp&
      &,  3.1325_GP,  2.9832_GP,  3.6078_GP,  3.8757_GP,  3.7616_GP,  3.9222_GP,  3.6370_gp&
      &,  3.8647_GP,  3.6256_GP,  3.4595_GP,  3.9874_GP,  4.1938_GP,  4.0679_GP,  4.0430_gp&
      &,  3.9781_GP,  3.3886_GP,  3.9008_GP,  3.8463_GP,  3.8288_GP,  3.7950_GP,  3.7790_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2241:2310)=(/&                                            
         3.7472_GP,  3.7117_GP,  3.9371_GP,  3.6873_GP,  3.9846_GP,  3.7709_GP,   3.6210_gp&
      &,  4.1812_GP,  4.3750_GP,  4.2044_GP,  4.1340_GP,  3.6459_GP,  3.5929_GP,  3.5036_gp&
      &,  3.4577_GP,  3.4528_GP,  3.4146_GP,  3.3904_GP,  3.9014_GP,  3.9031_GP,  4.1443_gp&
      &,  3.8961_GP,  4.2295_GP,  4.0227_GP,  3.8763_GP,  4.4086_GP,  4.7097_GP,  4.4064_gp&
      &,  4.4488_GP,  4.4243_GP,  4.4029_GP,  4.3842_GP,  4.3655_GP,  4.3514_GP,  4.1162_gp&
      &,  4.2205_GP,  4.1953_GP,  4.2794_GP,  2.8032_GP,  3.0805_GP,  3.4519_GP,  3.6700_gp&
      &,  3.4827_GP,  3.5050_GP,  3.2799_GP,  3.3482_GP,  3.1233_GP,  2.9747_GP,  3.5971_gp&
      &,  3.8586_GP,  3.7461_GP,  3.9100_GP,  3.6228_GP,  3.8535_GP,  3.6147_GP,  3.4490_gp&
      &,  3.9764_GP,  4.1773_GP,  4.0511_GP,  4.0270_GP,  3.9614_GP,  3.3754_GP,  3.8836_gp&
      &,  3.8291_GP,  3.8121_GP,  3.7780_GP,  3.7619_GP,  3.7300_GP,  3.6965_GP,  3.9253_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2311:2380)=(/&                                            
         3.6734_GP,  3.9733_GP,  3.7597_GP,  3.6099_GP,  4.1683_GP,  4.3572_GP,   4.1862_gp&
      &,  4.1153_GP,  3.6312_GP,  3.5772_GP,  3.4881_GP,  3.4429_GP,  3.4395_GP,  3.4009_gp&
      &,  3.3766_GP,  3.8827_GP,  3.8868_GP,  4.1316_GP,  3.8807_GP,  4.2164_GP,  4.0092_gp&
      &,  3.8627_GP,  4.3936_GP,  4.6871_GP,  4.3882_GP,  4.4316_GP,  4.4073_GP,  4.3858_gp&
      &,  4.3672_GP,  4.3485_GP,  4.3344_GP,  4.0984_GP,  4.2036_GP,  4.1791_GP,  4.2622_gp&
      &,  4.2450_GP,  2.7967_GP,  3.0689_GP,  3.4445_GP,  3.6581_GP,  3.4717_GP,  3.4951_gp&
      &,  3.2694_GP,  3.3397_GP,  3.1147_GP,  2.9661_GP,  3.5898_GP,  3.8468_GP,  3.7358_gp&
      &,  3.9014_GP,  3.6129_GP,  3.8443_GP,  3.6054_GP,  3.4396_GP,  3.9683_GP,  4.1656_gp&
      &,  4.0394_GP,  4.0158_GP,  3.9498_GP,  3.3677_GP,  3.8718_GP,  3.8164_GP,  3.8005_gp&
      &,  3.7662_GP,  3.7500_GP,  3.7181_GP,  3.6863_GP,  3.9170_GP,  3.6637_GP,  3.9641_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2381:2450)=(/&                                            
         3.7503_GP,  3.6004_GP,  4.1590_GP,  4.3448_GP,  4.1739_GP,  4.1029_GP,   3.6224_gp&
      &,  3.5677_GP,  3.4785_GP,  3.4314_GP,  3.4313_GP,  3.3923_GP,  3.3680_GP,  3.8698_gp&
      &,  3.8758_GP,  4.1229_GP,  3.8704_GP,  4.2063_GP,  3.9987_GP,  3.8519_GP,  4.3832_gp&
      &,  4.6728_GP,  4.3759_GP,  4.4195_GP,  4.3952_GP,  4.3737_GP,  4.3551_GP,  4.3364_gp&
      &,  4.3223_GP,  4.0861_GP,  4.1911_GP,  4.1676_GP,  4.2501_GP,  4.2329_GP,  4.2208_gp&
      &,  2.7897_GP,  3.0636_GP,  3.4344_GP,  3.6480_GP,  3.4626_GP,  3.4892_GP,  3.2626_gp&
      &,  3.3344_GP,  3.1088_GP,  2.9597_GP,  3.5804_GP,  3.8359_GP,  3.7251_GP,  3.8940_gp&
      &,  3.6047_GP,  3.8375_GP,  3.5990_GP,  3.4329_GP,  3.9597_GP,  4.1542_GP,  4.0278_gp&
      &,  4.0048_GP,  3.9390_GP,  3.3571_GP,  3.8608_GP,  3.8056_GP,  3.7899_GP,  3.7560_gp&
      &,  3.7400_GP,  3.7081_GP,  3.6758_GP,  3.9095_GP,  3.6552_GP,  3.9572_GP,  3.7436_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2451:2520)=(/&                                            
         3.5933_GP,  4.1508_GP,  4.3337_GP,  4.1624_GP,  4.0916_GP,  3.6126_GP,   3.5582_gp&
      &,  3.4684_GP,  3.4212_GP,  3.4207_GP,  3.3829_GP,  3.3586_GP,  3.8604_GP,  3.8658_gp&
      &,  4.1156_GP,  3.8620_GP,  4.1994_GP,  3.9917_GP,  3.8446_GP,  4.3750_GP,  4.6617_gp&
      &,  4.3644_GP,  4.4083_GP,  4.3840_GP,  4.3625_GP,  4.3439_GP,  4.3253_GP,  4.3112_gp&
      &,  4.0745_GP,  4.1807_GP,  4.1578_GP,  4.2390_GP,  4.2218_GP,  4.2097_GP,  4.1986_gp&
      &,  2.8395_GP,  3.0081_GP,  3.3171_GP,  3.4878_GP,  3.5360_GP,  3.5145_GP,  3.2809_gp&
      &,  3.3307_GP,  3.1260_GP,  2.9940_GP,  3.4741_GP,  3.6675_GP,  3.7832_GP,  3.6787_gp&
      &,  3.6156_GP,  3.8041_GP,  3.5813_GP,  3.4301_GP,  3.8480_GP,  3.9849_GP,  3.9314_gp&
      &,  3.8405_GP,  3.8029_GP,  3.2962_GP,  3.7104_GP,  3.6515_GP,  3.6378_GP,  3.6020_gp&
      &,  3.5849_GP,  3.5550_GP,  3.7494_GP,  3.6893_GP,  3.6666_GP,  3.9170_GP,  3.7150_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2521:2590)=(/&                                            
         3.5760_GP,  4.0268_GP,  4.1596_GP,  4.1107_GP,  3.9995_GP,  3.5574_GP,   3.5103_gp&
      &,  3.4163_GP,  3.3655_GP,  3.3677_GP,  3.3243_GP,  3.2975_GP,  3.7071_GP,  3.9047_gp&
      &,  3.8514_GP,  3.8422_GP,  3.8022_GP,  3.9323_GP,  3.7932_GP,  4.2343_GP,  4.4583_gp&
      &,  4.3115_GP,  4.2457_GP,  4.2213_GP,  4.1945_GP,  4.1756_GP,  4.1569_GP,  4.1424_gp&
      &,  4.0620_GP,  4.0494_GP,  3.9953_GP,  4.0694_GP,  4.0516_GP,  4.0396_GP,  4.0280_gp&
      &,  4.0130_GP,  2.9007_GP,  2.9674_GP,  3.8174_GP,  3.5856_GP,  3.6486_GP,  3.5339_gp&
      &,  3.2832_GP,  3.3154_GP,  3.1144_GP,  2.9866_GP,  3.9618_GP,  3.8430_GP,  3.9980_gp&
      &,  3.8134_GP,  3.6652_GP,  3.7985_GP,  3.5756_GP,  3.4207_GP,  4.4061_GP,  4.2817_gp&
      &,  4.1477_GP,  4.0616_GP,  3.9979_GP,  3.6492_GP,  3.8833_GP,  3.8027_GP,  3.7660_gp&
      &,  3.7183_GP,  3.6954_GP,  3.6525_GP,  3.9669_GP,  3.8371_GP,  3.7325_GP,  3.9160_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2591:2660)=(/&                                            
         3.7156_GP,  3.5714_GP,  4.6036_GP,  4.4620_GP,  4.3092_GP,  4.2122_GP,  3.8478_gp&
      &,  3.7572_GP,  3.6597_GP,  3.5969_GP,  3.5575_GP,  3.5386_GP,  3.5153_GP,  3.7818_gp&
      &,  4.1335_GP,  4.0153_GP,  3.9177_GP,  3.8603_GP,  3.9365_GP,  3.7906_GP,  4.7936_gp&
      &,  4.7410_GP,  4.5461_GP,  4.5662_GP,  4.5340_GP,  4.5059_GP,  4.4832_GP,  4.4604_gp&
      &,  4.4429_GP,  4.2346_GP,  4.4204_GP,  4.3119_GP,  4.3450_GP,  4.3193_GP,  4.3035_gp&
      &,  4.2933_GP,  4.1582_GP,  4.2450_GP,  2.8559_GP,  2.9050_GP,  3.8325_GP,  3.5442_gp&
      &,  3.5077_GP,  3.4905_GP,  3.2396_GP,  3.2720_GP,  3.0726_GP,  2.9467_GP,  3.9644_gp&
      &,  3.8050_GP,  3.8981_GP,  3.7762_GP,  3.6216_GP,  3.7531_GP,  3.5297_GP,  3.3742_gp&
      &,  4.3814_GP,  4.2818_GP,  4.1026_GP,  4.0294_GP,  3.9640_GP,  3.6208_GP,  3.8464_gp&
      &,  3.7648_GP,  3.7281_GP,  3.6790_GP,  3.6542_GP,  3.6117_GP,  3.8650_GP,  3.8010_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2661:2730)=(/&                                            
         3.6894_GP,  3.8713_GP,  3.6699_GP,  3.5244_GP,  4.5151_GP,  4.4517_GP,   4.2538_gp&
      &,  4.1483_GP,  3.8641_GP,  3.7244_GP,  3.6243_GP,  3.5589_GP,  3.5172_GP,  3.4973_gp&
      &,  3.4715_GP,  3.7340_GP,  4.0316_GP,  3.9958_GP,  3.8687_GP,  3.8115_GP,  3.8862_gp&
      &,  3.7379_GP,  4.7091_GP,  4.7156_GP,  4.5199_GP,  4.5542_GP,  4.5230_GP,  4.4959_gp&
      &,  4.4750_GP,  4.4529_GP,  4.4361_GP,  4.1774_GP,  4.3774_GP,  4.2963_GP,  4.3406_gp&
      &,  4.3159_GP,  4.3006_GP,  4.2910_GP,  4.1008_GP,  4.1568_GP,  4.0980_GP,  2.8110_gp&
      &,  2.8520_GP,  3.7480_GP,  3.5105_GP,  3.4346_GP,  3.3461_GP,  3.1971_GP,  3.2326_gp&
      &,  3.0329_GP,  2.9070_GP,  3.8823_GP,  3.7928_GP,  3.8264_GP,  3.7006_GP,  3.5797_gp&
      &,  3.7141_GP,  3.4894_GP,  3.3326_GP,  4.3048_GP,  4.2217_GP,  4.0786_GP,  3.9900_gp&
      &,  3.9357_GP,  3.6331_GP,  3.8333_GP,  3.7317_GP,  3.6957_GP,  3.6460_GP,  3.6197_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2731:2800)=(/&                                            
         3.5779_GP,  3.7909_GP,  3.7257_GP,  3.6476_GP,  3.5729_GP,  3.6304_GP,   3.4834_gp&
      &,  4.4368_GP,  4.3921_GP,  4.2207_GP,  4.1133_GP,  3.8067_GP,  3.7421_GP,  3.6140_gp&
      &,  3.5491_GP,  3.5077_GP,  3.4887_GP,  3.4623_GP,  3.6956_GP,  3.9568_GP,  3.8976_gp&
      &,  3.8240_GP,  3.7684_GP,  3.8451_GP,  3.6949_GP,  4.6318_GP,  4.6559_GP,  4.4533_gp&
      &,  4.4956_GP,  4.4641_GP,  4.4366_GP,  4.4155_GP,  4.3936_GP,  4.3764_GP,  4.1302_gp&
      &,  4.3398_GP,  4.2283_GP,  4.2796_GP,  4.2547_GP,  4.2391_GP,  4.2296_GP,  4.0699_gp&
      &,  4.1083_GP,  4.0319_GP,  3.9855_GP,  2.7676_GP,  2.8078_GP,  3.6725_GP,  3.4804_gp&
      &,  3.3775_GP,  3.2411_GP,  3.1581_GP,  3.1983_GP,  2.9973_GP,  2.8705_GP,  3.8070_gp&
      &,  3.7392_GP,  3.7668_GP,  3.6263_GP,  3.5402_GP,  3.6807_GP,  3.4545_GP,  3.2962_gp&
      &,  4.2283_GP,  4.1698_GP,  4.0240_GP,  3.9341_GP,  3.8711_GP,  3.5489_GP,  3.7798_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2801:2870)=(/&                                            
         3.7000_GP,  3.6654_GP,  3.6154_GP,  3.5882_GP,  3.5472_GP,  3.7289_GP,   3.6510_gp&
      &,  3.6078_GP,  3.5355_GP,  3.5963_GP,  3.4480_GP,  4.3587_GP,  4.3390_GP,  4.1635_gp&
      &,  4.0536_GP,  3.7193_GP,  3.6529_GP,  3.5512_GP,  3.4837_GP,  3.4400_GP,  3.4191_gp&
      &,  3.3891_GP,  3.6622_GP,  3.8934_GP,  3.8235_GP,  3.7823_GP,  3.7292_GP,  3.8106_gp&
      &,  3.6589_GP,  4.5535_GP,  4.6013_GP,  4.3961_GP,  4.4423_GP,  4.4109_GP,  4.3835_gp&
      &,  4.3625_GP,  4.3407_GP,  4.3237_GP,  4.0863_GP,  4.2835_GP,  4.1675_GP,  4.2272_gp&
      &,  4.2025_GP,  4.1869_GP,  4.1774_GP,  4.0126_GP,  4.0460_GP,  3.9815_GP,  3.9340_gp&
      &,  3.8955_GP,  2.6912_GP,  2.7604_GP,  3.6037_GP,  3.4194_GP,  3.3094_GP,  3.1710_gp&
      &,  3.0862_GP,  3.1789_GP,  2.9738_GP,  2.8427_GP,  3.7378_GP,  3.6742_GP,  3.6928_gp&
      &,  3.5512_GP,  3.4614_GP,  3.4087_GP,  3.4201_GP,  3.2607_GP,  4.1527_GP,  4.0977_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2871:2940)=(/&                                            
         3.9523_GP,  3.8628_GP,  3.8002_GP,  3.4759_GP,  3.7102_GP,  3.6466_GP,   3.6106_gp&
      &,  3.5580_GP,  3.5282_GP,  3.4878_GP,  3.6547_GP,  3.5763_GP,  3.5289_GP,  3.5086_gp&
      &,  3.5593_GP,  3.4099_GP,  4.2788_GP,  4.2624_GP,  4.0873_GP,  3.9770_GP,  3.6407_gp&
      &,  3.5743_GP,  3.5178_GP,  3.4753_GP,  3.3931_GP,  3.3694_GP,  3.3339_GP,  3.6002_gp&
      &,  3.8164_GP,  3.7478_GP,  3.7028_GP,  3.6952_GP,  3.7669_GP,  3.6137_GP,  4.4698_gp&
      &,  4.5488_GP,  4.3168_GP,  4.3646_GP,  4.3338_GP,  4.3067_GP,  4.2860_GP,  4.2645_gp&
      &,  4.2478_GP,  4.0067_GP,  4.2349_GP,  4.0958_GP,  4.1543_GP,  4.1302_GP,  4.1141_gp&
      &,  4.1048_GP,  3.9410_GP,  3.9595_GP,  3.8941_GP,  3.8465_GP,  3.8089_GP,  3.7490_gp&
      &,  2.7895_GP,  2.5849_GP,  3.6484_GP,  3.0162_GP,  3.1267_GP,  3.2125_GP,  3.0043_gp&
      &,  2.9572_GP,  2.8197_GP,  2.7261_GP,  3.7701_GP,  3.2446_GP,  3.5239_GP,  3.4696_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(2941:3010)=(/&                                            
         3.4261_GP,  3.3508_GP,  3.1968_GP,  3.0848_GP,  4.1496_GP,  3.6598_GP,   3.5111_gp&
      &,  3.4199_GP,  3.3809_GP,  3.5382_GP,  3.2572_GP,  3.2100_GP,  3.1917_GP,  3.1519_gp&
      &,  3.1198_GP,  3.1005_GP,  3.5071_GP,  3.5086_GP,  3.5073_GP,  3.4509_GP,  3.3120_gp&
      &,  3.2082_GP,  4.2611_GP,  3.8117_GP,  3.6988_GP,  3.5646_GP,  3.6925_GP,  3.6295_gp&
      &,  3.5383_GP,  3.4910_GP,  3.4625_GP,  3.4233_GP,  3.4007_GP,  3.2329_GP,  3.6723_gp&
      &,  3.6845_GP,  3.6876_GP,  3.6197_GP,  3.4799_GP,  3.3737_GP,  4.4341_GP,  4.0525_gp&
      &,  3.9011_GP,  3.8945_GP,  3.8635_GP,  3.8368_GP,  3.8153_GP,  3.7936_GP,  3.7758_gp&
      &,  3.4944_GP,  3.4873_GP,  3.9040_GP,  3.7110_GP,  3.6922_GP,  3.6799_GP,  3.6724_gp&
      &,  3.5622_GP,  3.6081_GP,  3.5426_GP,  3.4922_GP,  3.4498_GP,  3.3984_GP,  3.4456_gp&
      &,  2.7522_GP,  2.5524_GP,  3.5742_GP,  2.9508_GP,  3.0751_GP,  3.0158_GP,  2.9644_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3011:3080)=(/&                                            
         2.8338_GP,  2.7891_GP,  2.6933_GP,  3.6926_GP,  3.1814_GP,  3.4528_GP,   3.4186_gp&
      &,  3.3836_GP,  3.2213_GP,  3.1626_GP,  3.0507_GP,  4.0548_GP,  3.5312_GP,  3.4244_gp&
      &,  3.3409_GP,  3.2810_GP,  3.4782_GP,  3.1905_GP,  3.1494_GP,  3.1221_GP,  3.1128_gp&
      &,  3.0853_GP,  3.0384_GP,  3.4366_GP,  3.4562_GP,  3.4638_GP,  3.3211_GP,  3.2762_gp&
      &,  3.1730_GP,  4.1632_GP,  3.6825_GP,  3.5822_GP,  3.4870_GP,  3.6325_GP,  3.5740_gp&
      &,  3.4733_GP,  3.4247_GP,  3.3969_GP,  3.3764_GP,  3.3525_GP,  3.1984_GP,  3.5989_gp&
      &,  3.6299_GP,  3.6433_GP,  3.4937_GP,  3.4417_GP,  3.3365_GP,  4.3304_GP,  3.9242_gp&
      &,  3.7793_GP,  3.7623_GP,  3.7327_GP,  3.7071_GP,  3.6860_GP,  3.6650_GP,  3.6476_gp&
      &,  3.3849_GP,  3.3534_GP,  3.8216_GP,  3.5870_GP,  3.5695_GP,  3.5584_GP,  3.5508_gp&
      &,  3.4856_GP,  3.5523_GP,  3.4934_GP,  3.4464_GP,  3.4055_GP,  3.3551_GP,  3.3888_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3081:3150)=(/&                                            
         3.3525_GP,  2.7202_GP,  2.5183_GP,  3.4947_GP,  2.8731_GP,  3.0198_GP,   3.1457_gp&
      &,  2.9276_GP,  2.7826_GP,  2.7574_GP,  2.6606_GP,  3.6090_GP,  3.0581_GP,  3.3747_gp&
      &,  3.3677_GP,  3.3450_GP,  3.1651_GP,  3.1259_GP,  3.0147_GP,  3.9498_GP,  3.3857_gp&
      &,  3.2917_GP,  3.2154_GP,  3.1604_GP,  3.4174_GP,  3.0735_GP,  3.0342_GP,  3.0096_gp&
      &,  3.0136_GP,  2.9855_GP,  2.9680_GP,  3.3604_GP,  3.4037_GP,  3.4243_GP,  3.2633_gp&
      &,  3.1810_GP,  3.1351_GP,  4.0557_GP,  3.5368_GP,  3.4526_GP,  3.3699_GP,  3.5707_gp&
      &,  3.5184_GP,  3.4085_GP,  3.3595_GP,  3.3333_GP,  3.3143_GP,  3.3041_GP,  3.1094_gp&
      &,  3.5193_GP,  3.5745_GP,  3.6025_GP,  3.4338_GP,  3.3448_GP,  3.2952_GP,  4.2158_gp&
      &,  3.7802_GP,  3.6431_GP,  3.6129_GP,  3.5853_GP,  3.5610_GP,  3.5406_GP,  3.5204_gp&
      &,  3.5036_GP,  3.2679_GP,  3.2162_GP,  3.7068_GP,  3.4483_GP,  3.4323_GP,  3.4221_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3151:3220)=(/&                                            
         3.4138_GP,  3.3652_GP,  3.4576_GP,  3.4053_GP,  3.3618_GP,  3.3224_GP,   3.2711_gp&
      &,  3.3326_GP,  3.2950_GP,  3.2564_GP,  2.5315_GP,  2.6104_GP,  3.2734_GP,  3.2299_gp&
      &,  3.1090_GP,  2.9942_GP,  2.9159_GP,  2.8324_GP,  2.8350_GP,  2.7216_GP,  3.3994_gp&
      &,  3.4475_GP,  3.4354_GP,  3.3438_GP,  3.2807_GP,  3.2169_GP,  3.2677_GP,  3.1296_gp&
      &,  3.7493_GP,  3.8075_GP,  3.6846_GP,  3.6104_GP,  3.5577_GP,  3.2052_GP,  3.4803_gp&
      &,  3.4236_GP,  3.3845_GP,  3.3640_GP,  3.3365_GP,  3.3010_GP,  3.3938_GP,  3.3624_gp&
      &,  3.3440_GP,  3.3132_GP,  3.4035_GP,  3.2754_GP,  3.8701_GP,  3.9523_GP,  3.8018_gp&
      &,  3.7149_GP,  3.3673_GP,  3.3199_GP,  3.2483_GP,  3.2069_GP,  3.1793_GP,  3.1558_gp&
      &,  3.1395_GP,  3.4097_GP,  3.5410_GP,  3.5228_GP,  3.5116_GP,  3.4921_GP,  3.4781_gp&
      &,  3.4690_GP,  4.0420_GP,  4.1759_GP,  4.0078_GP,  4.0450_GP,  4.0189_GP,  3.9952_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3221:3290)=(/&                                            
         3.9770_GP,  3.9583_GP,  3.9434_GP,  3.7217_GP,  3.8228_GP,  3.7826_GP,   3.8640_gp&
      &,  3.8446_GP,  3.8314_GP,  3.8225_GP,  3.6817_GP,  3.7068_GP,  3.6555_GP,  3.6159_gp&
      &,  3.5831_GP,  3.5257_GP,  3.2133_GP,  3.1689_GP,  3.1196_GP,  3.3599_GP,  2.9852_gp&
      &,  2.7881_GP,  3.5284_GP,  3.3493_GP,  3.6958_GP,  3.3642_GP,  3.1568_GP,  3.0055_gp&
      &,  2.9558_GP,  2.8393_GP,  3.6287_GP,  3.5283_GP,  4.1511_GP,  3.8259_GP,  3.6066_gp&
      &,  3.4527_GP,  3.3480_GP,  3.2713_GP,  3.9037_GP,  3.8361_GP,  3.8579_GP,  3.7311_gp&
      &,  3.6575_GP,  3.5176_GP,  3.5693_GP,  3.5157_GP,  3.4814_GP,  3.4559_GP,  3.4445_gp&
      &,  3.4160_GP,  4.1231_GP,  3.8543_GP,  3.6816_GP,  3.5602_GP,  3.4798_GP,  3.4208_gp&
      &,  4.0542_GP,  4.0139_GP,  4.0165_GP,  3.9412_GP,  3.7698_GP,  3.6915_GP,  3.6043_gp&
      &,  3.5639_GP,  3.5416_GP,  3.5247_GP,  3.5153_GP,  3.5654_GP,  4.2862_GP,  4.0437_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3291:3360)=(/&                                            
         3.8871_GP,  3.7741_GP,  3.6985_GP,  3.6413_GP,  4.2345_GP,  4.3663_GP,   4.3257_gp&
      &,  4.0869_GP,  4.0612_GP,  4.0364_GP,  4.0170_GP,  3.9978_GP,  3.9834_GP,  3.9137_gp&
      &,  3.8825_GP,  3.8758_GP,  3.9143_GP,  3.8976_GP,  3.8864_GP,  3.8768_GP,  3.9190_gp&
      &,  4.1613_GP,  4.0566_GP,  3.9784_GP,  3.9116_GP,  3.8326_GP,  3.7122_GP,  3.6378_gp&
      &,  3.5576_GP,  3.5457_GP,  4.3127_GP,  3.1160_GP,  2.8482_GP,  4.0739_GP,  3.3599_gp&
      &,  3.5698_GP,  3.5366_GP,  3.2854_GP,  3.1039_GP,  2.9953_GP,  2.9192_GP,  4.1432_gp&
      &,  3.5320_GP,  3.9478_GP,  4.0231_GP,  3.7509_GP,  3.5604_GP,  3.4340_GP,  3.3426_gp&
      &,  4.3328_GP,  3.8288_GP,  3.7822_GP,  3.7909_GP,  3.6907_GP,  3.6864_GP,  3.5793_gp&
      &,  3.5221_GP,  3.4883_GP,  3.4649_GP,  3.4514_GP,  3.4301_GP,  3.9256_GP,  4.0596_gp&
      &,  3.8307_GP,  3.6702_GP,  3.5651_GP,  3.4884_GP,  4.4182_GP,  4.2516_GP,  3.9687_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3361:3430)=(/&                                            
         3.9186_GP,  3.9485_GP,  3.8370_GP,  3.7255_GP,  3.6744_GP,  3.6476_GP,   3.6295_gp&
      &,  3.6193_GP,  3.5659_GP,  4.0663_GP,  4.2309_GP,  4.0183_GP,  3.8680_GP,  3.7672_gp&
      &,  3.6923_GP,  4.5240_GP,  4.4834_GP,  4.1570_GP,  4.3204_GP,  4.2993_GP,  4.2804_gp&
      &,  4.2647_GP,  4.2481_GP,  4.2354_GP,  3.8626_GP,  3.8448_GP,  4.2267_GP,  4.1799_gp&
      &,  4.1670_GP,  3.8738_GP,  3.8643_GP,  3.8796_GP,  4.0575_GP,  4.0354_GP,  3.9365_gp&
      &,  3.8611_GP,  3.7847_GP,  3.7388_GP,  3.6826_GP,  3.6251_GP,  3.5492_GP,  4.0889_gp&
      &,  4.2764_GP,  3.1416_GP,  2.8325_GP,  3.7735_GP,  3.3787_GP,  3.4632_GP,  3.5923_gp&
      &,  3.3214_GP,  3.1285_GP,  3.0147_GP,  2.9366_GP,  3.8527_GP,  3.5602_GP,  3.8131_gp&
      &,  3.8349_GP,  3.7995_GP,  3.5919_GP,  3.4539_GP,  3.3540_GP,  4.0654_GP,  3.8603_gp&
      &,  3.7972_GP,  3.7358_GP,  3.7392_GP,  3.8157_GP,  3.6055_GP,  3.5438_GP,  3.5089_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3431:3500)=(/&                                            
         3.4853_GP,  3.4698_GP,  3.4508_GP,  3.7882_GP,  3.8682_GP,  3.8837_GP,   3.7055_gp&
      &,  3.5870_GP,  3.5000_GP,  4.1573_GP,  4.0005_GP,  3.9568_GP,  3.8936_GP,  3.9990_gp&
      &,  3.9433_GP,  3.8172_GP,  3.7566_GP,  3.7246_GP,  3.7033_GP,  3.6900_GP,  3.5697_gp&
      &,  3.9183_GP,  4.0262_GP,  4.0659_GP,  3.8969_GP,  3.7809_GP,  3.6949_GP,  4.2765_gp&
      &,  4.2312_GP,  4.1401_GP,  4.0815_GP,  4.0580_GP,  4.0369_GP,  4.0194_GP,  4.0017_gp&
      &,  3.9874_GP,  3.8312_GP,  3.8120_GP,  3.9454_GP,  3.9210_GP,  3.9055_GP,  3.8951_gp&
      &,  3.8866_GP,  3.8689_GP,  3.9603_GP,  3.9109_GP,  3.9122_GP,  3.8233_GP,  3.7438_gp&
      &,  3.7436_GP,  3.6981_GP,  3.6555_GP,  3.5452_GP,  3.9327_GP,  4.0658_GP,  4.1175_gp&
      &,  2.9664_GP,  2.8209_GP,  3.5547_GP,  3.3796_GP,  3.3985_GP,  3.3164_GP,  3.2364_gp&
      &,  3.1956_GP,  3.0370_GP,  2.9313_GP,  3.6425_GP,  3.5565_GP,  3.7209_GP,  3.7108_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3501:3570)=(/&                                            
         3.6639_GP,  3.6484_GP,  3.4745_GP,  3.3492_GP,  3.8755_GP,  4.2457_GP,   3.7758_gp&
      &,  3.7161_GP,  3.6693_GP,  3.6155_GP,  3.5941_GP,  3.5643_GP,  3.5292_GP,  3.4950_gp&
      &,  3.4720_GP,  3.4503_GP,  3.6936_GP,  3.7392_GP,  3.7388_GP,  3.7602_GP,  3.6078_gp&
      &,  3.4960_GP,  3.9800_GP,  4.3518_GP,  4.2802_GP,  3.8580_GP,  3.8056_GP,  3.7527_gp&
      &,  3.7019_GP,  3.6615_GP,  3.5768_GP,  3.5330_GP,  3.5038_GP,  3.5639_GP,  3.8192_gp&
      &,  3.8883_GP,  3.9092_GP,  3.9478_GP,  3.7995_GP,  3.6896_GP,  4.1165_GP,  4.5232_gp&
      &,  4.4357_GP,  4.4226_GP,  4.4031_GP,  4.3860_GP,  4.3721_GP,  4.3580_GP,  4.3466_gp&
      &,  4.2036_GP,  4.2037_GP,  3.8867_GP,  4.2895_GP,  4.2766_GP,  4.2662_GP,  4.2598_gp&
      &,  3.8408_GP,  3.9169_GP,  3.8681_GP,  3.8250_GP,  3.7855_GP,  3.7501_GP,  3.6753_gp&
      &,  3.5499_GP,  3.4872_GP,  3.5401_GP,  3.8288_GP,  3.9217_GP,  3.9538_GP,  4.0054_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3571:3640)=(/&                                            
         2.8388_GP,  2.7890_GP,  3.4329_GP,  3.5593_GP,  3.3488_GP,  3.2486_GP,   3.1615_gp&
      &,  3.1000_GP,  3.0394_GP,  2.9165_GP,  3.5267_GP,  3.7479_GP,  3.6650_GP,  3.6263_gp&
      &,  3.5658_GP,  3.5224_GP,  3.4762_GP,  3.3342_GP,  3.7738_GP,  4.0333_GP,  3.9568_gp&
      &,  3.8975_GP,  3.8521_GP,  3.4929_GP,  3.7830_GP,  3.7409_GP,  3.7062_GP,  3.6786_gp&
      &,  3.6471_GP,  3.6208_GP,  3.6337_GP,  3.6519_GP,  3.6363_GP,  3.6278_GP,  3.6110_gp&
      &,  3.4825_GP,  3.8795_GP,  4.1448_GP,  4.0736_GP,  4.0045_GP,  3.6843_GP,  3.6291_gp&
      &,  3.5741_GP,  3.5312_GP,  3.4974_GP,  3.4472_GP,  3.4034_GP,  3.7131_GP,  3.7557_gp&
      &,  3.7966_GP,  3.8005_GP,  3.8068_GP,  3.8015_GP,  3.6747_GP,  4.0222_GP,  4.3207_gp&
      &,  4.2347_GP,  4.2191_GP,  4.1990_GP,  4.1811_GP,  4.1666_GP,  4.1521_GP,  4.1401_gp&
      &,  3.9970_GP,  3.9943_GP,  3.9592_GP,  4.0800_GP,  4.0664_GP,  4.0559_GP,  4.0488_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3641:3710)=(/&                                            
         3.9882_GP,  4.0035_GP,  3.9539_GP,  3.9138_GP,  3.8798_GP,  3.8355_GP,   3.5359_gp&
      &,  3.4954_GP,  3.3962_GP,  3.5339_GP,  3.7595_GP,  3.8250_GP,  3.8408_GP,  3.8600_gp&
      &,  3.8644_GP,  2.7412_GP,  2.7489_GP,  3.3374_GP,  3.3950_GP,  3.3076_GP,  3.1910_gp&
      &,  3.0961_GP,  3.0175_GP,  3.0280_GP,  2.8929_GP,  3.4328_GP,  3.5883_GP,  3.6227_gp&
      &,  3.5616_GP,  3.4894_GP,  3.4241_GP,  3.3641_GP,  3.3120_GP,  3.6815_GP,  3.8789_gp&
      &,  3.8031_GP,  3.7413_GP,  3.6939_GP,  3.4010_GP,  3.6225_GP,  3.5797_GP,  3.5443_gp&
      &,  3.5139_GP,  3.4923_GP,  3.4642_GP,  3.5860_GP,  3.5849_GP,  3.5570_GP,  3.5257_gp&
      &,  3.4936_GP,  3.4628_GP,  3.7874_GP,  3.9916_GP,  3.9249_GP,  3.8530_GP,  3.5932_gp&
      &,  3.5355_GP,  3.4757_GP,  3.4306_GP,  3.3953_GP,  3.3646_GP,  3.3390_GP,  3.5637_gp&
      &,  3.7053_GP,  3.7266_GP,  3.7177_GP,  3.6996_GP,  3.6775_GP,  3.6558_GP,  3.9331_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3711:3780)=(/&                                            
         4.1655_GP,  4.0879_GP,  4.0681_GP,  4.0479_GP,  4.0299_GP,  4.0152_GP,   4.0006_gp&
      &,  3.9883_GP,  3.8500_GP,  3.8359_GP,  3.8249_GP,  3.9269_GP,  3.9133_GP,  3.9025_gp&
      &,  3.8948_GP,  3.8422_GP,  3.8509_GP,  3.7990_GP,  3.7570_GP,  3.7219_GP,  3.6762_gp&
      &,  3.4260_GP,  3.3866_GP,  3.3425_GP,  3.5294_GP,  3.7022_GP,  3.7497_GP,  3.7542_gp&
      &,  3.7494_GP,  3.7370_GP,  3.7216_GP,  3.4155_GP,  3.0522_GP,  4.2541_GP,  3.8218_gp&
      &,  4.0438_GP,  3.5875_GP,  3.3286_GP,  3.1682_GP,  3.0566_GP,  2.9746_GP,  4.3627_gp&
      &,  4.0249_GP,  4.6947_GP,  4.1718_GP,  3.8639_GP,  3.6735_GP,  3.5435_GP,  3.4479_gp&
      &,  4.6806_GP,  4.3485_GP,  4.2668_GP,  4.1690_GP,  4.1061_GP,  4.1245_GP,  4.0206_gp&
      &,  3.9765_GP,  3.9458_GP,  3.9217_GP,  3.9075_GP,  3.8813_GP,  3.9947_GP,  4.1989_gp&
      &,  3.9507_GP,  3.7960_GP,  3.6925_GP,  3.6150_GP,  4.8535_GP,  4.5642_GP,  4.4134_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3781:3850)=(/&                                            
         4.3688_GP,  4.3396_GP,  4.2879_GP,  4.2166_GP,  4.1888_GP,  4.1768_GP,   4.1660_gp&
      &,  4.1608_GP,  4.0745_GP,  4.2289_GP,  4.4863_GP,  4.2513_GP,  4.0897_GP,  3.9876_gp&
      &,  3.9061_GP,  5.0690_GP,  5.0446_GP,  4.6186_GP,  4.6078_GP,  4.5780_GP,  4.5538_gp&
      &,  4.5319_GP,  4.5101_GP,  4.4945_GP,  4.1912_GP,  4.2315_GP,  4.5534_GP,  4.4373_gp&
      &,  4.4224_GP,  4.4120_GP,  4.4040_GP,  4.2634_GP,  4.7770_GP,  4.6890_GP,  4.6107_gp&
      &,  4.5331_GP,  4.4496_GP,  4.4082_GP,  4.3095_GP,  4.2023_GP,  4.0501_GP,  4.2595_gp&
      &,  4.5497_GP,  4.3056_GP,  4.1506_GP,  4.0574_GP,  3.9725_GP,  5.0796_GP,  3.0548_gp&
      &,  3.3206_GP,  3.8132_GP,  3.9720_GP,  3.7675_GP,  3.7351_GP,  3.5167_GP,  3.5274_gp&
      &,  3.3085_GP,  3.1653_GP,  3.9500_GP,  4.1730_GP,  4.0613_GP,  4.1493_GP,  3.8823_gp&
      &,  4.0537_GP,  3.8200_GP,  3.6582_GP,  4.3422_GP,  4.5111_GP,  4.3795_GP,  4.3362_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3851:3920)=(/&                                            
         4.2751_GP,  3.7103_GP,  4.1973_GP,  4.1385_GP,  4.1129_GP,  4.0800_GP,   4.0647_gp&
      &,  4.0308_GP,  4.0096_GP,  4.1619_GP,  3.9360_GP,  4.1766_GP,  3.9705_GP,  3.8262_gp&
      &,  4.5348_GP,  4.7025_GP,  4.5268_GP,  4.5076_GP,  3.9562_GP,  3.9065_GP,  3.8119_gp&
      &,  3.7605_GP,  3.7447_GP,  3.7119_GP,  3.6916_GP,  4.1950_GP,  4.2110_GP,  4.3843_gp&
      &,  4.1631_GP,  4.4427_GP,  4.2463_GP,  4.1054_GP,  4.7693_GP,  5.0649_GP,  4.7365_gp&
      &,  4.7761_GP,  4.7498_GP,  4.7272_GP,  4.7076_GP,  4.6877_GP,  4.6730_GP,  4.4274_gp&
      &,  4.5473_GP,  4.5169_GP,  4.5975_GP,  4.5793_GP,  4.5667_GP,  4.5559_GP,  4.3804_gp&
      &,  4.6920_GP,  4.6731_GP,  4.6142_GP,  4.5600_GP,  4.4801_GP,  4.0149_GP,  3.8856_gp&
      &,  3.7407_GP,  4.1545_GP,  4.2253_GP,  4.4229_GP,  4.1923_GP,  4.5022_GP,  4.3059_gp&
      &,  4.1591_GP,  4.7883_GP,  4.9294_GP,  3.3850_GP,  3.4208_GP,  3.7004_GP,  3.8800_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3921:3990)=(/&                                            
         3.9886_GP,  3.9040_GP,  3.6719_GP,  3.6547_GP,  3.4625_GP,  3.3370_GP,   3.8394_gp&
      &,  4.0335_GP,  4.2373_GP,  4.3023_GP,  4.0306_GP,  4.1408_GP,  3.9297_GP,  3.7857_gp&
      &,  4.1907_GP,  4.3230_GP,  4.2664_GP,  4.2173_GP,  4.1482_GP,  3.6823_GP,  4.0711_gp&
      &,  4.0180_GP,  4.0017_GP,  3.9747_GP,  3.9634_GP,  3.9383_GP,  4.1993_GP,  4.3205_gp&
      &,  4.0821_GP,  4.2547_GP,  4.0659_GP,  3.9359_GP,  4.3952_GP,  4.5176_GP,  4.3888_gp&
      &,  4.3607_GP,  3.9583_GP,  3.9280_GP,  3.8390_GP,  3.7971_GP,  3.7955_GP,  3.7674_gp&
      &,  3.7521_GP,  4.1062_GP,  4.3633_GP,  4.2991_GP,  4.2767_GP,  4.4857_GP,  4.3039_gp&
      &,  4.1762_GP,  4.6197_GP,  4.8654_GP,  4.6633_GP,  4.5878_GP,  4.5640_GP,  4.5422_gp&
      &,  4.5231_GP,  4.5042_GP,  4.4901_GP,  4.3282_GP,  4.3978_GP,  4.3483_GP,  4.4202_gp&
      &,  4.4039_GP,  4.3926_GP,  4.3807_GP,  4.2649_GP,  4.6135_GP,  4.5605_GP,  4.5232_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(3991:4060)=(/&                                            
         4.4676_GP,  4.3948_GP,  4.0989_GP,  3.9864_GP,  3.8596_GP,  4.0942_GP,   4.2720_gp&
      &,  4.3270_GP,  4.3022_GP,  4.5410_GP,  4.3576_GP,  4.2235_GP,  4.6545_GP,  4.7447_gp&
      &,  4.7043_GP,  3.0942_GP,  3.2075_GP,  3.5152_GP,  3.6659_GP,  3.8289_GP,  3.7459_gp&
      &,  3.5156_GP,  3.5197_GP,  3.3290_GP,  3.2069_GP,  3.6702_GP,  3.8448_GP,  4.0340_gp&
      &,  3.9509_GP,  3.8585_GP,  3.9894_GP,  3.7787_GP,  3.6365_GP,  4.1425_GP,  4.1618_gp&
      &,  4.0940_GP,  4.0466_GP,  3.9941_GP,  3.5426_GP,  3.8952_GP,  3.8327_GP,  3.8126_gp&
      &,  3.7796_GP,  3.7635_GP,  3.7356_GP,  4.0047_GP,  3.9655_GP,  3.9116_GP,  4.1010_gp&
      &,  3.9102_GP,  3.7800_GP,  4.2964_GP,  4.3330_GP,  4.2622_GP,  4.2254_GP,  3.8195_gp&
      &,  3.7560_GP,  3.6513_GP,  3.5941_GP,  3.5810_GP,  3.5420_GP,  3.5178_GP,  3.8861_gp&
      &,  4.1459_GP,  4.1147_GP,  4.0772_GP,  4.3120_GP,  4.1207_GP,  3.9900_GP,  4.4733_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(4061:4130)=(/&                                            
         4.6157_GP,  4.4580_GP,  4.4194_GP,  4.3954_GP,  4.3739_GP,  4.3531_GP,   4.3343_gp&
      &,  4.3196_GP,  4.2140_GP,  4.2339_GP,  4.1738_GP,  4.2458_GP,  4.2278_GP,  4.2158_gp&
      &,  4.2039_GP,  4.1658_GP,  4.3595_GP,  4.2857_GP,  4.2444_GP,  4.1855_GP,  4.1122_gp&
      &,  3.7839_GP,  3.6879_GP,  3.5816_GP,  3.8633_GP,  4.1585_GP,  4.1402_GP,  4.1036_gp&
      &,  4.3694_GP,  4.1735_GP,  4.0368_GP,  4.5095_GP,  4.5538_GP,  4.5240_GP,  4.4252_gp&
      &,  3.0187_GP,  3.1918_GP,  3.5127_GP,  3.6875_GP,  3.7404_GP,  3.6943_GP,  3.4702_gp&
      &,  3.4888_GP,  3.2914_GP,  3.1643_GP,  3.6669_GP,  3.8724_GP,  3.9940_GP,  4.0816_gp&
      &,  3.8054_GP,  3.9661_GP,  3.7492_GP,  3.6024_GP,  4.0428_GP,  4.1951_GP,  4.1466_gp&
      &,  4.0515_GP,  4.0075_GP,  3.5020_GP,  3.9158_GP,  3.8546_GP,  3.8342_GP,  3.8008_gp&
      &,  3.7845_GP,  3.7549_GP,  3.9602_GP,  3.8872_GP,  3.8564_GP,  4.0793_GP,  3.8835_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(4131:4200)=(/&                                            
         3.7495_GP,  4.2213_GP,  4.3704_GP,  4.3300_GP,  4.2121_GP,  3.7643_GP,   3.7130_gp&
      &,  3.6144_GP,  3.5599_GP,  3.5474_GP,  3.5093_GP,  3.4853_GP,  3.9075_GP,  4.1115_gp&
      &,  4.0473_GP,  4.0318_GP,  4.2999_GP,  4.1050_GP,  3.9710_GP,  4.4320_GP,  4.6706_gp&
      &,  4.5273_GP,  4.4581_GP,  4.4332_GP,  4.4064_GP,  4.3873_GP,  4.3684_GP,  4.3537_gp&
      &,  4.2728_GP,  4.2549_GP,  4.2032_GP,  4.2794_GP,  4.2613_GP,  4.2491_GP,  4.2375_gp&
      &,  4.2322_GP,  4.3665_GP,  4.3061_GP,  4.2714_GP,  4.2155_GP,  4.1416_GP,  3.7660_gp&
      &,  3.6628_GP,  3.5476_GP,  3.8790_GP,  4.1233_GP,  4.0738_GP,  4.0575_GP,  4.3575_gp&
      &,  4.1586_GP,  4.0183_GP,  4.4593_GP,  4.5927_GP,  4.4865_GP,  4.3813_GP,  4.4594_gp&
      &,  2.9875_GP,  3.1674_GP,  3.4971_GP,  3.6715_GP,  3.7114_GP,  3.6692_GP,  3.4446_gp&
      &,  3.4676_GP,  3.2685_GP,  3.1405_GP,  3.6546_GP,  3.8579_GP,  3.9637_GP,  4.0581_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(4201:4270)=(/&                                            
         3.7796_GP,  3.9463_GP,  3.7275_GP,  3.5792_GP,  4.0295_GP,  4.1824_GP,   4.1247_gp&
      &,  4.0357_GP,  3.9926_GP,  3.4827_GP,  3.9007_GP,  3.8392_GP,  3.8191_GP,  3.7851_gp&
      &,  3.7687_GP,  3.7387_GP,  3.9290_GP,  3.8606_GP,  3.8306_GP,  4.0601_GP,  3.8625_gp&
      &,  3.7269_GP,  4.2062_GP,  4.3566_GP,  4.3022_GP,  4.1929_GP,  3.7401_GP,  3.6888_gp&
      &,  3.5900_GP,  3.5350_GP,  3.5226_GP,  3.4838_GP,  3.4594_GP,  3.8888_GP,  4.0813_gp&
      &,  4.0209_GP,  4.0059_GP,  4.2810_GP,  4.0843_GP,  3.9486_GP,  4.4162_GP,  4.6542_gp&
      &,  4.5005_GP,  4.4444_GP,  4.4196_GP,  4.3933_GP,  4.3741_GP,  4.3552_GP,  4.3406_gp&
      &,  4.2484_GP,  4.2413_GP,  4.1907_GP,  4.2656_GP,  4.2474_GP,  4.2352_GP,  4.2236_gp&
      &,  4.2068_GP,  4.3410_GP,  4.2817_GP,  4.2479_GP,  4.1921_GP,  4.1182_GP,  3.7346_gp&
      &,  3.6314_GP,  3.5168_GP,  3.8582_GP,  4.0927_GP,  4.0469_GP,  4.0313_GP,  4.3391_gp&
      &/)                                                                              
       r0ab_ref_table_grimme(4271:4340)=(/&                                            
         4.1381_GP,  3.9962_GP,  4.4429_GP,  4.5787_GP,  4.4731_GP,  4.3588_GP,   4.4270_gp&
      &,  4.3957_GP,  2.9659_GP,  3.1442_GP,  3.4795_GP,  3.6503_GP,  3.6814_GP,  3.6476_gp&
      &,  3.4222_GP,  3.4491_GP,  3.2494_GP,  3.1209_GP,  3.6324_GP,  3.8375_GP,  3.9397_gp&
      &,  3.8311_GP,  3.7581_GP,  3.9274_GP,  3.7085_GP,  3.5598_GP,  4.0080_GP,  4.1641_gp&
      &,  4.1057_GP,  4.0158_GP,  3.9726_GP,  3.4667_GP,  3.8802_GP,  3.8188_GP,  3.7989_gp&
      &,  3.7644_GP,  3.7474_GP,  3.7173_GP,  3.9049_GP,  3.8424_GP,  3.8095_GP,  4.0412_gp&
      &,  3.8436_GP,  3.7077_GP,  4.1837_GP,  4.3366_GP,  4.2816_GP,  4.1686_GP,  3.7293_gp&
      &,  3.6709_GP,  3.5700_GP,  3.5153_GP,  3.5039_GP,  3.4684_GP,  3.4437_GP,  3.8663_gp&
      &,  4.0575_GP,  4.0020_GP,  3.9842_GP,  4.2612_GP,  4.0643_GP,  3.9285_GP,  4.3928_gp&
      &,  4.6308_GP,  4.4799_GP,  4.4244_GP,  4.3996_GP,  4.3737_GP,  4.3547_GP,  4.3358_gp&
      &/)
       r0ab_ref_table_grimme(4341:4410)=(/&
         4.3212_GP,  4.2275_GP,  4.2216_GP,  4.1676_GP,  4.2465_GP,  4.2283_GP,   4.2161_gp&
      &,  4.2045_GP,  4.1841_GP,  4.3135_GP,  4.2562_GP,  4.2226_GP,  4.1667_GP,  4.0932_gp&
      &,  3.7134_GP,  3.6109_GP,  3.4962_GP,  3.8352_GP,  4.0688_GP,  4.0281_GP,  4.0099_gp&
      &,  4.3199_GP,  4.1188_GP,  3.9768_GP,  4.4192_GP,  4.5577_GP,  4.4516_GP,  4.3365_gp&
      &,  4.4058_GP,  4.3745_GP,  4.3539_GP,  2.8763_GP,  3.1294_GP,  3.5598_GP,  3.7465_gp&
      &,  3.5659_GP,  3.5816_GP,  3.3599_GP,  3.4024_GP,  3.1877_GP,  3.0484_GP,  3.7009_gp&
      &,  3.9451_GP,  3.8465_GP,  3.9873_GP,  3.7079_GP,  3.9083_GP,  3.6756_GP,  3.5150_gp&
      &,  4.0829_GP,  4.2780_GP,  4.1511_GP,  4.1260_GP,  4.0571_GP,  3.4865_GP,  3.9744_gp&
      &,  3.9150_GP,  3.8930_GP,  3.8578_GP,  3.8402_GP,  3.8073_GP,  3.7977_GP,  4.0036_gp&
      &,  3.7604_GP,  4.0288_GP,  3.8210_GP,  3.6757_GP,  4.2646_GP,  4.4558_GP,  4.2862_gp&
      &/)
       r0ab_ref_table_grimme(4411:4465)=(/&
         4.2122_GP,  3.7088_GP,  3.6729_GP,  3.5800_GP,  3.5276_GP,  3.5165_GP,   3.4783_gp&
      &,  3.4539_GP,  3.9553_GP,  3.9818_GP,  4.2040_GP,  3.9604_GP,  4.2718_GP,  4.0689_gp&
      &,  3.9253_GP,  4.4869_GP,  4.7792_GP,  4.4918_GP,  4.5342_GP,  4.5090_GP,  4.4868_gp&
      &,  4.4680_GP,  4.4486_GP,  4.4341_GP,  4.2023_GP,  4.3122_GP,  4.2710_GP,  4.3587_gp&
      &,  4.3407_GP,  4.3281_GP,  4.3174_GP,  4.1499_GP,  4.3940_GP,  4.3895_GP,  4.3260_gp&
      &,  4.2725_GP,  4.1961_GP,  3.7361_GP,  3.6193_GP,  3.4916_GP,  3.9115_GP,  3.9914_gp&
      &,  3.9809_GP,  3.9866_GP,  4.3329_GP,  4.1276_GP,  3.9782_GP,  4.5097_GP,  4.6769_gp&
      &,  4.5158_GP,  4.3291_GP,  4.3609_GP,  4.3462_GP,  4.3265_GP,  4.4341_gp&
      &/)
 !
 !     Convert from lower triangle to full square matrix of R0AB
 !
       k=0
       do i=1,MAX_ELEM
          do j=1,i
             k=k+1
             vdwparams%r0AB(i,j)=r0ab_ref_table_grimme(k)*angs2au 
             vdwparams%r0AB(j,i)=r0ab_ref_table_grimme(k)*angs2au 
          enddo
       enddo
       
       call f_free(r0ab_ref_table_grimme)
  END SUBROUTINE init_r0ab_d3


!> initialize coordination number dependent c6ab directly from Grimme's tables
  subroutine init_c6_params_d3
    implicit none
    !local variables
    !REAL(kind=GP), DIMENSION(161925) :: grimme_c6_parameters
    !integer, parameter :: c6dim=161925
    !real(gp), dimension(:), allocatable :: grimme_c6_parameters
    integer            :: nlines, idata,iline,iat,jat,iadr,jadr,il,jl
    !real(gp) :: tti,ttj

    vdwparams%coeffs=f_malloc_ptr((ntot_values*(ntot_values+1))/2,id='coeffs')

!
!     Number of distinct C6 value, atom pair combinations
!
      nlines=32385

      do iline=1,MAX_ELEM
         vdwparams%maxcn(iline)=0
      enddo

    !enter the simplified values for the parameters
    do idata=1,ntot_values
       vdwparams%ivalues(idata)=idata
    end do
    vdwparams%ivalues( 95)=101
    vdwparams%ivalues( 96)=103
    vdwparams%ivalues( 97)=104
    vdwparams%ivalues( 98)=105
    vdwparams%ivalues( 99)=106
    vdwparams%ivalues(100)=107
    vdwparams%ivalues(101)=108
    vdwparams%ivalues(102)=109
    vdwparams%ivalues(103)=111
    vdwparams%ivalues(104)=112
    vdwparams%ivalues(105)=113
    vdwparams%ivalues(106)=114
    vdwparams%ivalues(107)=115
    vdwparams%ivalues(108)=116
    vdwparams%ivalues(109)=117
    vdwparams%ivalues(110)=119
    vdwparams%ivalues(111)=120
    vdwparams%ivalues(112)=121
    vdwparams%ivalues(113)=122
    vdwparams%ivalues(114)=123
    vdwparams%ivalues(115)=124
    vdwparams%ivalues(116)=125
    vdwparams%ivalues(117)=126
    vdwparams%ivalues(118)=127
    vdwparams%ivalues(119)=128
    vdwparams%ivalues(120)=129
    vdwparams%ivalues(121)=130
    vdwparams%ivalues(122)=131
    vdwparams%ivalues(123)=132
    vdwparams%ivalues(124)=133
    vdwparams%ivalues(125)=134
    vdwparams%ivalues(126)=135
    vdwparams%ivalues(127)=137
    vdwparams%ivalues(128)=138
    vdwparams%ivalues(129)=139
    vdwparams%ivalues(130)=140
    vdwparams%ivalues(131)=141
    vdwparams%ivalues(132)=142
    vdwparams%ivalues(133)=143
    vdwparams%ivalues(134)=144
    vdwparams%ivalues(135)=145
    vdwparams%ivalues(136)=146
    vdwparams%ivalues(137)=147
    vdwparams%ivalues(138)=148
    vdwparams%ivalues(139)=149
    vdwparams%ivalues(140)=150
    vdwparams%ivalues(141)=151
    vdwparams%ivalues(142)=152
    vdwparams%ivalues(143)=153
    vdwparams%ivalues(144)=155
    vdwparams%ivalues(145)=156
    vdwparams%ivalues(146)=157
    vdwparams%ivalues(147)=159
    vdwparams%ivalues(148)=160
    vdwparams%ivalues(149)=161
    vdwparams%ivalues(150)=162
    vdwparams%ivalues(151)=163
    vdwparams%ivalues(152)=164
    vdwparams%ivalues(153)=165
    vdwparams%ivalues(154)=166
    vdwparams%ivalues(155)=167
    vdwparams%ivalues(156)=168
    vdwparams%ivalues(157)=169
    vdwparams%ivalues(158)=170
    vdwparams%ivalues(159)=171
    vdwparams%ivalues(160)=172
    vdwparams%ivalues(161)=173
    vdwparams%ivalues(162)=174
    vdwparams%ivalues(163)=175
    vdwparams%ivalues(164)=176
    vdwparams%ivalues(165)=177
    vdwparams%ivalues(166)=178
    vdwparams%ivalues(167)=179
    vdwparams%ivalues(168)=180
    vdwparams%ivalues(169)=181
    vdwparams%ivalues(170)=182
    vdwparams%ivalues(171)=183
    vdwparams%ivalues(172)=184
    vdwparams%ivalues(173)=185
    vdwparams%ivalues(174)=187
    vdwparams%ivalues(175)=188
    vdwparams%ivalues(176)=189
    vdwparams%ivalues(177)=190
    vdwparams%ivalues(178)=191
    vdwparams%ivalues(179)=192
    vdwparams%ivalues(180)=193
    vdwparams%ivalues(181)=194
    vdwparams%ivalues(182)=204
    vdwparams%ivalues(183)=205
    vdwparams%ivalues(184)=206
    vdwparams%ivalues(185)=207
    vdwparams%ivalues(186)=208
    vdwparams%ivalues(187)=212
    vdwparams%ivalues(188)=213
    vdwparams%ivalues(189)=214
    vdwparams%ivalues(190)=215
    vdwparams%ivalues(191)=216
    vdwparams%ivalues(192)=220
    vdwparams%ivalues(193)=221
    vdwparams%ivalues(194)=222
    vdwparams%ivalues(195)=223
    vdwparams%ivalues(196)=224
    vdwparams%ivalues(197)=225
    vdwparams%ivalues(198)=226
    vdwparams%ivalues(199)=227
    vdwparams%ivalues(200)=228
    vdwparams%ivalues(201)=231
    vdwparams%ivalues(202)=232
    vdwparams%ivalues(203)=233
    vdwparams%ivalues(204)=234
    vdwparams%ivalues(205)=238
    vdwparams%ivalues(206)=239
    vdwparams%ivalues(207)=240
    vdwparams%ivalues(208)=241
    vdwparams%ivalues(209)=242
    vdwparams%ivalues(210)=243
    vdwparams%ivalues(211)=244
    vdwparams%ivalues(212)=245
    vdwparams%ivalues(213)=246
    vdwparams%ivalues(214)=249
    vdwparams%ivalues(215)=250
    vdwparams%ivalues(216)=251
    vdwparams%ivalues(217)=252
    vdwparams%ivalues(218)=256
    vdwparams%ivalues(219)=257
    vdwparams%ivalues(220)=272
    vdwparams%ivalues(221)=273
    vdwparams%ivalues(222)=274
    vdwparams%ivalues(223)=275
    vdwparams%ivalues(224)=276
    vdwparams%ivalues(225)=277
    vdwparams%ivalues(226)=278
    vdwparams%ivalues(227)=281
    vdwparams%ivalues(228)=282
    vdwparams%ivalues(229)=283
    vdwparams%ivalues(230)=284
    vdwparams%ivalues(231)=288
    vdwparams%ivalues(232)=305
    vdwparams%ivalues(233)=306
    vdwparams%ivalues(234)=307
    vdwparams%ivalues(235)=313
    vdwparams%ivalues(236)=314
    vdwparams%ivalues(237)=315
    vdwparams%ivalues(238)=327
    vdwparams%ivalues(239)=328
    vdwparams%ivalues(240)=331
    vdwparams%ivalues(241)=332
    vdwparams%ivalues(242)=333
    vdwparams%ivalues(243)=349
    vdwparams%ivalues(244)=350
    vdwparams%ivalues(245)=351
    vdwparams%ivalues(246)=381
    vdwparams%ivalues(247)=382
    vdwparams%ivalues(248)=383
    vdwparams%ivalues(249)=405
    vdwparams%ivalues(250)=406
    vdwparams%ivalues(251)=414
    vdwparams%ivalues(252)=432
    vdwparams%ivalues(253)=450
    vdwparams%ivalues(254)=482

    vdwparams%rvalues=0.0_gp
    vdwparams%rvalues(  1)= 0.91180_gp
    vdwparams%rvalues( 58)= 2.79910_gp
    vdwparams%rvalues( 96)= 0.98650_gp
    vdwparams%rvalues( 97)= 0.98080_gp
    vdwparams%rvalues( 98)= 0.97060_gp
    vdwparams%rvalues( 99)= 0.98680_gp
    vdwparams%rvalues(100)= 0.99440_gp
    vdwparams%rvalues(101)= 0.99250_gp
    vdwparams%rvalues(102)= 0.99820_gp
    vdwparams%rvalues(103)= 0.96840_gp
    vdwparams%rvalues(104)= 0.96280_gp
    vdwparams%rvalues(105)= 0.96480_gp
    vdwparams%rvalues(106)= 0.95070_gp
    vdwparams%rvalues(107)= 0.99470_gp
    vdwparams%rvalues(108)= 0.99480_gp
    vdwparams%rvalues(109)= 0.99720_gp
    vdwparams%rvalues(110)= 0.97670_gp
    vdwparams%rvalues(111)= 0.98310_gp
    vdwparams%rvalues(112)= 1.86270_gp
    vdwparams%rvalues(113)= 1.82990_gp
    vdwparams%rvalues(114)= 1.91380_gp
    vdwparams%rvalues(115)= 1.82690_gp
    vdwparams%rvalues(116)= 1.64060_gp
    vdwparams%rvalues(117)= 1.64830_gp
    vdwparams%rvalues(118)= 1.71490_gp
    vdwparams%rvalues(119)= 1.79370_gp
    vdwparams%rvalues(120)= 0.95760_gp
    vdwparams%rvalues(121)= 1.94190_gp
    vdwparams%rvalues(122)= 0.96010_gp
    vdwparams%rvalues(123)= 0.94340_gp
    vdwparams%rvalues(124)= 0.98890_gp
    vdwparams%rvalues(125)= 0.99010_gp
    vdwparams%rvalues(126)= 0.99740_gp
    vdwparams%rvalues(127)= 0.97380_gp
    vdwparams%rvalues(128)= 0.98010_gp
    vdwparams%rvalues(129)= 1.91530_gp
    vdwparams%rvalues(130)= 1.93550_gp
    vdwparams%rvalues(131)= 1.95450_gp
    vdwparams%rvalues(132)= 1.94200_gp
    vdwparams%rvalues(133)= 1.66820_gp
    vdwparams%rvalues(134)= 1.85840_gp
    vdwparams%rvalues(135)= 1.90030_gp
    vdwparams%rvalues(136)= 1.86300_gp
    vdwparams%rvalues(137)= 0.96790_gp
    vdwparams%rvalues(138)= 1.95390_gp
    vdwparams%rvalues(139)= 0.96330_gp
    vdwparams%rvalues(140)= 0.95140_gp
    vdwparams%rvalues(141)= 0.97490_gp
    vdwparams%rvalues(142)= 0.98110_gp
    vdwparams%rvalues(143)= 0.99680_gp
    vdwparams%rvalues(144)= 0.99090_gp
    vdwparams%rvalues(145)= 0.97970_gp
    vdwparams%rvalues(146)= 1.93730_gp
    vdwparams%rvalues(147)= 2.94250_gp
    vdwparams%rvalues(148)= 2.94550_gp
    vdwparams%rvalues(149)= 2.94130_gp
    vdwparams%rvalues(150)= 2.93000_gp
    vdwparams%rvalues(151)= 1.82860_gp
    vdwparams%rvalues(152)= 2.87320_gp
    vdwparams%rvalues(153)= 2.90860_gp
    vdwparams%rvalues(154)= 2.89650_gp
    vdwparams%rvalues(155)= 2.92420_gp
    vdwparams%rvalues(156)= 2.92820_gp
    vdwparams%rvalues(157)= 2.92460_gp
    vdwparams%rvalues(158)= 2.84820_gp
    vdwparams%rvalues(159)= 2.92190_gp
    vdwparams%rvalues(160)= 1.92540_gp
    vdwparams%rvalues(161)= 1.94590_gp
    vdwparams%rvalues(162)= 1.92920_gp
    vdwparams%rvalues(163)= 1.81040_gp
    vdwparams%rvalues(164)= 1.88580_gp
    vdwparams%rvalues(165)= 1.86480_gp
    vdwparams%rvalues(166)= 1.91880_gp
    vdwparams%rvalues(167)= 0.98460_gp
    vdwparams%rvalues(168)= 1.98960_gp
    vdwparams%rvalues(169)= 0.92670_gp
    vdwparams%rvalues(170)= 0.93830_gp
    vdwparams%rvalues(171)= 0.98200_gp
    vdwparams%rvalues(172)= 0.98150_gp
    vdwparams%rvalues(173)= 0.99540_gp
    vdwparams%rvalues(174)= 0.97050_gp
    vdwparams%rvalues(175)= 0.96620_gp
    vdwparams%rvalues(176)= 2.90700_gp
    vdwparams%rvalues(177)= 2.88440_gp
    vdwparams%rvalues(178)= 2.87380_gp
    vdwparams%rvalues(179)= 2.88780_gp
    vdwparams%rvalues(180)= 2.90950_gp
    vdwparams%rvalues(181)= 1.92090_gp
    vdwparams%rvalues(182)= 1.96970_gp
    vdwparams%rvalues(183)= 1.94410_gp
    vdwparams%rvalues(184)= 1.99850_gp
    vdwparams%rvalues(185)= 2.01430_gp
    vdwparams%rvalues(186)= 1.98870_gp
    vdwparams%rvalues(187)= 1.94960_gp
    vdwparams%rvalues(188)= 1.93110_gp
    vdwparams%rvalues(189)= 1.94350_gp
    vdwparams%rvalues(190)= 2.01020_gp
    vdwparams%rvalues(191)= 1.99030_gp
    vdwparams%rvalues(192)= 1.93490_gp
    vdwparams%rvalues(193)= 2.89990_gp
    vdwparams%rvalues(194)= 3.86750_gp
    vdwparams%rvalues(195)= 2.91100_gp
    vdwparams%rvalues(196)=10.61910_gp
    vdwparams%rvalues(197)= 9.88490_gp
    vdwparams%rvalues(198)= 9.13760_gp
    vdwparams%rvalues(199)= 2.92630_gp
    vdwparams%rvalues(200)= 6.54580_gp
    vdwparams%rvalues(201)= 1.93150_gp
    vdwparams%rvalues(202)= 1.94470_gp
    vdwparams%rvalues(203)= 1.97930_gp
    vdwparams%rvalues(204)= 1.98120_gp
    vdwparams%rvalues(205)= 1.91430_gp
    vdwparams%rvalues(206)= 2.89030_gp
    vdwparams%rvalues(207)= 3.91060_gp
    vdwparams%rvalues(208)= 2.92250_gp
    vdwparams%rvalues(209)=11.05560_gp
    vdwparams%rvalues(210)= 9.54020_gp
    vdwparams%rvalues(211)= 8.88950_gp
    vdwparams%rvalues(212)= 2.96960_gp
    vdwparams%rvalues(213)= 5.70950_gp
    vdwparams%rvalues(214)= 1.93780_gp
    vdwparams%rvalues(215)= 1.95050_gp
    vdwparams%rvalues(216)= 1.95230_gp
    vdwparams%rvalues(217)= 1.96390_gp
    vdwparams%rvalues(218)= 1.84670_gp
    vdwparams%rvalues(219)= 2.91750_gp
    vdwparams%rvalues(220)= 3.88400_gp
    vdwparams%rvalues(221)= 2.89880_gp
    vdwparams%rvalues(222)=10.91530_gp
    vdwparams%rvalues(223)= 9.80540_gp
    vdwparams%rvalues(224)= 9.15270_gp
    vdwparams%rvalues(225)= 2.94240_gp
    vdwparams%rvalues(226)= 6.66690_gp
    vdwparams%rvalues(227)= 1.93020_gp
    vdwparams%rvalues(228)= 1.93560_gp
    vdwparams%rvalues(229)= 1.96550_gp
    vdwparams%rvalues(230)= 1.96390_gp
    vdwparams%rvalues(231)= 1.80750_gp
    vdwparams%rvalues(232)= 2.91280_gp
    vdwparams%rvalues(233)= 2.99870_gp
    vdwparams%rvalues(234)= 2.99030_gp
    vdwparams%rvalues(235)= 2.91460_gp
    vdwparams%rvalues(236)= 2.94070_gp
    vdwparams%rvalues(237)= 2.98590_gp
    vdwparams%rvalues(238)= 7.77850_gp
    vdwparams%rvalues(239)= 6.29180_gp
    vdwparams%rvalues(240)= 2.92330_gp
    vdwparams%rvalues(241)= 2.91860_gp
    vdwparams%rvalues(242)= 2.97090_gp
    vdwparams%rvalues(243)= 2.93530_gp
    vdwparams%rvalues(244)= 2.92590_gp
    vdwparams%rvalues(245)= 2.93150_gp
    vdwparams%rvalues(246)= 2.94200_gp
    vdwparams%rvalues(247)= 2.90810_gp
    vdwparams%rvalues(248)= 2.95000_gp
    vdwparams%rvalues(249)= 4.58560_gp
    vdwparams%rvalues(250)= 3.98440_gp
    vdwparams%rvalues(251)= 3.86770_gp
    vdwparams%rvalues(252)= 3.89720_gp
    vdwparams%rvalues(253)= 3.91230_gp
    vdwparams%rvalues(254)= 3.90980_gp

    vdwparams%coeffs(    1:  100)=(/    3.02670_gp,    2.08350_gp,    1.55830_gp,   38.94480_gp,   22.15080_gp,&
         1163.44540_gp,   24.44150_gp,   14.82460_gp,  494.61900_gp,  257.48630_gp,&
         17.31430_gp,   11.09750_gp,  283.73080_gp,  161.59710_gp,  107.17770_gp,&
         12.14020_gp,    8.18410_gp,  169.90300_gp,  102.95600_gp,   71.27940_gp,&
         49.11300_gp,    8.71710_gp,    6.13800_gp,  108.48540_gp,   68.64580_gp,&
         49.11320_gp,   34.81460_gp,   25.26850_gp,    6.71800_gp,    4.89490_gp,&
         76.96130_gp,   50.12520_gp,   36.72470_gp,   26.59290_gp,   19.65460_gp,&
         15.50590_gp,    5.16160_gp,    3.88250_gp,   55.09330_gp,   36.74530_gp,&
         27.48210_gp,   20.28270_gp,   15.24180_gp,   12.18340_gp,    9.69160_gp,&
         4.01120_gp,    3.10250_gp,   40.47310_gp,   27.48670_gp,   20.90220_gp,&
         15.67400_gp,   11.94790_gp,    9.66060_gp,    7.76910_gp,    6.28960_gp,&
         46.82320_gp,   26.86280_gp, 1367.32720_gp,  587.45630_gp,  338.72120_gp,&
         203.76310_gp,  130.65630_gp,   93.02630_gp,   66.84230_gp,   49.27990_gp,&
         1608.02860_gp,   38.35310_gp,   23.03200_gp,  830.81560_gp,  418.21640_gp,&
         258.13030_gp,  162.60820_gp,  107.61500_gp,   78.22500_gp,   57.16050_gp,&
         42.67710_gp,  985.16970_gp,  683.37580_gp,   36.29090_gp,   22.32240_gp,&
         705.82540_gp,  372.63020_gp,  236.47800_gp,  152.09340_gp,  102.20000_gp,&
         75.07550_gp,   55.34120_gp,   41.59660_gp,  838.96480_gp,  603.46890_gp,&
         540.54060_gp,   29.59470_gp,   18.85000_gp,  495.34490_gp,  279.78630_gp,&
         184.51110_gp,  122.13870_gp,   83.84980_gp,   62.53490_gp,   46.69360_gp/)
    vdwparams%coeffs(  101:  200)=(/   35.45500_gp,  591.04580_gp,  447.64230_gp,  408.96060_gp,  317.85740_gp,&
         23.76040_gp,   15.66890_gp,  350.80300_gp,  208.73310_gp,  142.34810_gp,&
         96.75030_gp,   67.78730_gp,   51.30980_gp,   38.80840_gp,   29.77670_gp,&
         420.00640_gp,  330.78010_gp,  307.29650_gp,  244.35460_gp,  191.68870_gp,&
         20.09480_gp,   13.61080_gp,  273.78670_gp,  167.95130_gp,  117.11210_gp,&
         81.09190_gp,   57.67340_gp,   44.14700_gp,   33.72640_gp,   26.09400_gp,&
         328.59900_gp,  264.66650_gp,  248.50080_gp,  200.53740_gp,  159.48980_gp,&
         134.00660_gp,   16.70520_gp,   11.63020_gp,  210.66260_gp,  132.98080_gp,&
         94.76120_gp,   66.84070_gp,   48.26240_gp,   37.36880_gp,   28.84450_gp,&
         22.51210_gp,  253.51360_gp,  208.49780_gp,  197.75940_gp,  161.86860_gp,&
         130.47250_gp,  110.70060_gp,   92.34600_gp,   13.87000_gp,    9.91300_gp,&
         163.54970_gp,  105.72290_gp,   76.79490_gp,   55.08980_gp,   40.34350_gp,&
         31.57830_gp,   24.61780_gp,   19.37740_gp,  197.34400_gp,  165.10060_gp,&
         157.95950_gp,  130.89270_gp,  106.76980_gp,   91.40140_gp,   76.93830_gp,&
         64.64620_gp,   76.23760_gp,   44.04110_gp, 2387.15740_gp,  972.31970_gp,&
         554.19800_gp,  332.26940_gp,  213.20010_gp,  152.14200_gp,  109.64390_gp,&
         81.08610_gp, 2798.61240_gp, 1642.05870_gp, 1387.93610_gp,  967.62830_gp,&
         684.49680_gp,  535.23840_gp,  413.15120_gp,  322.11550_gp, 4983.50090_gp,&
         65.81800_gp,   39.07010_gp, 1614.47190_gp,  757.90700_gp,  454.84390_gp,&
         281.70350_gp,  184.54980_gp,  133.38470_gp,   97.08240_gp,   72.30720_gp/)
    vdwparams%coeffs(  201:  300)=(/ 1907.70810_gp, 1252.59320_gp, 1088.84200_gp,  790.63470_gp,  575.46270_gp,&
         456.84430_gp,  357.47370_gp,  281.68230_gp, 3240.43930_gp, 2352.68620_gp,&
         54.96690_gp,   32.89950_gp, 1278.11830_gp,  617.33580_gp,  375.04180_gp,&
         234.19500_gp,  154.29470_gp,  111.93590_gp,   81.72420_gp,   61.01540_gp,&
         1512.53380_gp, 1015.54540_gp,  888.60250_gp,  651.20400_gp,  477.27970_gp,&
         380.41070_gp,  298.76860_gp,  236.14170_gp, 2549.94120_gp, 1888.79020_gp,&
         1522.46760_gp,   53.68750_gp,   32.53160_gp, 1192.91280_gp,  587.48980_gp,&
         361.08650_gp,  227.57230_gp,  151.01990_gp,  110.14520_gp,   80.79750_gp,&
         60.55940_gp, 1413.15700_gp,  962.98310_gp,  847.46150_gp,  626.23530_gp,&
         462.32470_gp,  370.26310_gp,  292.18530_gp,  231.92080_gp, 2374.66900_gp,&
         1779.51620_gp, 1438.28410_gp, 1361.91850_gp,   49.48190_gp,   30.18510_gp,&
         1069.04260_gp,  533.34990_gp,  330.10230_gp,  209.14500_gp,  139.35120_gp,&
         101.93250_gp,   74.96710_gp,   56.31080_gp, 1267.31390_gp,  872.25790_gp,&
         770.34320_gp,  572.11290_gp,  424.14700_gp,  340.60200_gp,  269.48630_gp,&
         214.40020_gp, 2124.19850_gp, 1604.91780_gp, 1299.51970_gp, 1232.32350_gp,&
         1116.09840_gp,   39.12210_gp,   24.14630_gp,  844.09360_gp,  418.16350_gp,&
         259.30210_gp,  164.98820_gp,  110.48610_gp,   81.20000_gp,   60.00990_gp,&
         45.28200_gp, 1000.46520_gp,  684.40150_gp,  604.40710_gp,  449.22160_gp,&
         333.83800_gp,  268.80290_gp,  213.36680_gp,  170.34650_gp, 1683.67010_gp,&
         1263.08090_gp, 1021.68710_gp,  968.85650_gp,  877.38550_gp,  690.74250_gp/)
    vdwparams%coeffs(  301:  400)=(/   43.00280_gp,   26.49780_gp,  891.97560_gp,  453.36860_gp,  283.47440_gp,&
         181.00350_gp,  121.32330_gp,   89.13460_gp,   65.81120_gp,   49.59590_gp,&
         1058.51030_gp,  739.04110_gp,  656.07380_gp,  490.81560_gp,  366.12040_gp,&
         295.17370_gp,  234.45330_gp,  187.17280_gp, 1768.03220_gp, 1351.47850_gp,&
         1097.15390_gp, 1042.66250_gp,  945.60720_gp,  743.31790_gp,  802.74840_gp,&
         33.91100_gp,   21.17500_gp,  698.45390_gp,  353.53580_gp,  221.75470_gp,&
         142.34670_gp,   95.97820_gp,   70.89420_gp,   52.63100_gp,   39.86590_gp,&
         828.89360_gp,  576.51210_gp,  512.08230_gp,  383.74450_gp,  287.16060_gp,&
         232.26250_gp,  185.18230_gp,  148.43190_gp, 1388.95220_gp, 1056.51930_gp,&
         857.17300_gp,  814.84020_gp,  739.05420_gp,  581.81140_gp,  627.54030_gp,&
         491.33490_gp,   36.32340_gp,   22.72330_gp,  701.53240_gp,  369.19280_gp,&
         234.86930_gp,  151.87030_gp,  102.75030_gp,   75.99120_gp,   56.43070_gp,&
         42.72700_gp,  834.18260_gp,  598.33020_gp,  535.97060_gp,  405.99070_gp,&
         305.92560_gp,  248.21160_gp,  198.35700_gp,  159.19570_gp, 1382.97710_gp,&
         1081.76900_gp,  882.45030_gp,  841.83350_gp,  765.31780_gp,  601.41780_gp,&
         651.96880_gp,  509.77640_gp,  532.77940_gp,   37.15960_gp,   23.07920_gp,&
         738.85800_gp,  383.91400_gp,  242.41740_gp,  155.85070_gp,  104.98340_gp,&
         77.40100_gp,   57.32480_gp,   43.31290_gp,  877.93540_gp,  623.63050_gp,&
         556.53890_gp,  419.35930_gp,  314.56980_gp,  254.46830_gp,  202.77050_gp,&
         162.32640_gp, 1458.45010_gp, 1132.10740_gp,  921.89900_gp,  878.09390_gp/)
    vdwparams%coeffs(  401:  500)=(/  797.51950_gp,  626.64720_gp,  678.45310_gp,  530.32210_gp,  553.08340_gp,&
         574.74360_gp,   28.59400_gp,   18.02160_gp,  569.45260_gp,  292.79550_gp,&
         185.17410_gp,  119.63500_gp,   81.08550_gp,   60.13350_gp,   44.80980_gp,&
         34.05440_gp,  676.47730_gp,  476.24260_gp,  424.77230_gp,  320.19130_gp,&
         240.79740_gp,  195.40750_gp,  156.31610_gp,  125.67760_gp, 1129.95180_gp,&
         868.36740_gp,  706.07200_gp,  672.41280_gp,  610.57000_gp,  480.69380_gp,&
         519.30780_gp,  406.72680_gp,  423.08770_gp,  439.63650_gp,  337.18080_gp,&
         29.86890_gp,   18.96950_gp,  540.91600_gp,  293.31060_gp,  189.62930_gp,&
         124.10320_gp,   84.73110_gp,   63.07740_gp,   47.11220_gp,   35.84270_gp,&
         644.34080_gp,  472.91430_gp,  427.16100_gp,  327.28100_gp,  248.97430_gp,&
         203.24270_gp,  163.38480_gp,  131.81190_gp, 1062.59280_gp,  846.85270_gp,&
         693.67870_gp,  664.10090_gp,  605.06230_gp,  475.55940_gp,  517.09100_gp,&
         404.55580_gp,  424.87900_gp,  440.09000_gp,  336.64930_gp,  340.52130_gp,&
         35.16970_gp,   22.04580_gp,  656.74410_gp,  351.11410_gp,  225.39330_gp,&
         146.55440_gp,   99.44070_gp,   73.62530_gp,   54.67820_gp,   41.36810_gp,&
         781.35240_gp,  567.16380_gp,  510.58540_gp,  389.26560_gp,  294.79250_gp,&
         239.81830_gp,  192.05780_gp,  154.34230_gp, 1292.81210_gp, 1019.88350_gp,&
         833.74240_gp,  796.86170_gp,  725.25140_gp,  569.76450_gp,  618.83870_gp,&
         483.78150_gp,  507.13160_gp,  525.81620_gp,  401.96130_gp,  405.45650_gp,&
         483.75160_gp,   31.81700_gp,   20.43180_gp,  527.92680_gp,  298.42200_gp/)
    vdwparams%coeffs(  501:  600)=(/  197.34480_gp,  131.09540_gp,   90.33280_gp,   67.59320_gp,   50.64120_gp,&
         38.57460_gp,  630.10410_gp,  477.37440_gp,  436.50600_gp,  339.83280_gp,&
         261.82350_gp,  215.31010_gp,  174.19690_gp,  141.20630_gp, 1032.66510_gp,&
         843.36440_gp,  694.70240_gp,  668.37860_gp,  610.76390_gp,  479.96390_gp,&
         524.18480_gp,  410.21060_gp,  433.85080_gp,  447.98680_gp,  342.42610_gp,&
         349.98850_gp,  415.95300_gp,  363.54740_gp,   27.78840_gp,   18.32840_gp,&
         415.15320_gp,  245.48310_gp,  166.90600_gp,  113.26160_gp,   79.30860_gp,&
         60.02890_gp,   45.41800_gp,   34.86860_gp,  496.93380_gp,  389.49270_gp,&
         361.19130_gp,  286.59090_gp,  224.47740_gp,  186.63480_gp,  152.60350_gp,&
         124.85610_gp,  810.65540_gp,  679.00440_gp,  562.69470_gp,  544.67970_gp,&
         499.49700_gp,  393.19220_gp,  430.91460_gp,  338.00960_gp,  359.71500_gp,&
         370.03590_gp,  283.33230_gp,  292.50300_gp,  346.35080_gp,  307.10100_gp,&
         262.94980_gp,   25.30980_gp,   17.02430_gp,  354.54290_gp,  215.02380_gp,&
         148.80080_gp,  102.43460_gp,   72.54030_gp,   55.36270_gp,   42.19230_gp,&
         32.58430_gp,  425.18640_gp,  339.57420_gp,  317.61820_gp,  255.00990_gp,&
         201.90190_gp,  169.13080_gp,  139.31630_gp,  114.74490_gp,  692.91160_gp,&
         587.96470_gp,  488.87150_gp,  475.03480_gp,  436.57180_gp,  344.29560_gp,&
         377.82760_gp,  297.03960_gp,  317.01740_gp,  325.34190_gp,  249.63850_gp,&
         259.04590_gp,  305.95010_gp,  273.65300_gp,  236.34730_gp,  213.67380_gp,&
         22.48340_gp,   15.45530_gp,  294.78970_gp,  183.37420_gp,  129.24950_gp/)
    vdwparams%coeffs(  601:  700)=(/   90.34240_gp,   64.76320_gp,   49.87950_gp,   38.32250_gp,   29.79580_gp,&
         354.28690_gp,  288.27360_gp,  272.00980_gp,  221.05370_gp,  176.99300_gp,&
         149.45630_gp,  124.09090_gp,  102.94710_gp,  577.23180_gp,  495.97540_gp,&
         413.74680_gp,  403.65760_gp,  371.80710_gp,  293.91850_gp,  322.84330_gp,&
         254.53060_gp,  272.30530_gp,  278.77050_gp,  214.50590_gp,  223.63480_gp,&
         263.34580_gp,  237.62950_gp,  207.07720_gp,  188.36050_gp,  167.12970_gp,&
         19.81820_gp,   13.92110_gp,  244.55540_gp,  155.56160_gp,  111.53710_gp,&
         79.10800_gp,   57.38950_gp,   44.60030_gp,   34.54520_gp,   27.04300_gp,&
         294.55740_gp,  243.59650_gp,  231.69040_gp,  190.39340_gp,  154.05750_gp,&
         131.09280_gp,  109.68390_gp,   91.64300_gp,  480.24720_gp,  416.97770_gp,&
         348.84740_gp,  341.62730_gp,  315.32510_gp,  249.94130_gp,  274.64190_gp,&
         217.20710_gp,  232.75940_gp,  237.74470_gp,  183.53470_gp,  192.04840_gp,&
         225.45370_gp,  205.05980_gp,  180.18180_gp,  164.84940_gp,  147.18300_gp,&
         130.40170_gp,   85.94990_gp,   50.04910_gp, 2647.33310_gp, 1082.70520_gp,&
         620.08780_gp,  373.53840_gp,  240.69400_gp,  172.33980_gp,  124.58790_gp,&
         92.38270_gp, 3104.10160_gp, 1826.46480_gp, 1546.92290_gp, 1082.05630_gp,&
         768.13350_gp,  602.22210_gp,  466.15650_gp,  364.40670_gp, 5530.28060_gp,&
         3599.43490_gp, 2834.09200_gp, 2641.24350_gp, 2363.71750_gp, 1874.03260_gp,&
         1968.76310_gp, 1547.19820_gp, 1541.89970_gp, 1625.10340_gp, 1259.44270_gp,&
         1186.26310_gp, 1442.38070_gp, 1155.15260_gp,  909.42570_gp,  778.92040_gp/)
    vdwparams%coeffs(  701:  800)=(/  650.36060_gp,  542.31600_gp, 6138.77550_gp,   78.39010_gp,   46.66580_gp,&
         1944.22600_gp,  903.92650_gp,  541.57290_gp,  335.45970_gp,  219.96080_gp,&
         159.14600_gp,  115.96470_gp,   86.46270_gp, 2296.19060_gp, 1495.81410_gp,&
         1298.53560_gp,  941.41210_gp,  684.96320_gp,  543.96270_gp,  425.89970_gp,&
         335.86900_gp, 3915.35290_gp, 2819.49350_gp, 2260.42170_gp, 2128.39510_gp,&
         1918.69410_gp, 1511.30430_gp, 1614.71900_gp, 1263.23270_gp, 1290.87180_gp,&
         1351.37020_gp, 1037.72840_gp, 1009.68090_gp, 1216.51380_gp, 1004.48240_gp,&
         808.33700_gp,  700.05320_gp,  590.73260_gp,  496.91200_gp, 4349.25770_gp,&
         3381.36720_gp,   70.05700_gp,   42.41950_gp, 1585.39500_gp,  771.48440_gp,&
         472.43640_gp,  297.21710_gp,  197.06110_gp,  143.66020_gp,  105.34470_gp,&
         78.93160_gp, 1876.86340_gp, 1266.76020_gp, 1112.29800_gp,  819.54180_gp,&
         603.99160_gp,  483.37970_gp,  381.24620_gp,  302.51430_gp, 3166.06250_gp,&
         2350.80180_gp, 1896.78090_gp, 1794.31800_gp, 1622.46890_gp, 1276.43520_gp,&
         1371.44920_gp, 1072.32650_gp, 1105.30160_gp, 1153.57370_gp,  884.16480_gp,&
         870.69350_gp, 1045.55050_gp,  874.78360_gp,  711.75060_gp,  620.33340_gp,&
         526.82990_gp,  445.70130_gp, 3520.94660_gp, 2813.70070_gp, 2365.89250_gp,&
         63.78010_gp,   39.10880_gp, 1363.28420_gp,  681.52830_gp,  423.23020_gp,&
         269.04900_gp,  179.77770_gp,  131.78520_gp,   97.09770_gp,   73.03150_gp,&
         1616.28630_gp, 1113.83390_gp,  985.04940_gp,  733.18070_gp,  544.88820_gp,&
         438.38900_gp,  347.53340_gp,  276.99080_gp, 2711.15420_gp, 2048.52510_gp/)
    vdwparams%coeffs(  801:  900)=(/ 1659.11560_gp, 1574.14010_gp, 1426.06870_gp, 1121.54890_gp, 1208.74510_gp,&
         945.17130_gp,  978.97670_gp, 1019.75890_gp,  781.11330_gp,  774.58710_gp,&
         928.19820_gp,  782.94420_gp,  641.55640_gp,  561.54740_gp,  479.00400_gp,&
         406.87270_gp, 3017.68560_gp, 2449.47230_gp, 2072.61700_gp, 1822.71810_gp,&
         58.67680_gp,   36.34180_gp, 1209.03040_gp,  613.98500_gp,  384.90720_gp,&
         246.52410_gp,  165.69430_gp,  121.98680_gp,   90.22340_gp,   68.07690_gp,&
         1434.66410_gp, 1000.56630_gp,  889.03020_gp,  666.15070_gp,  497.97030_gp,&
         402.19740_gp,  320.06760_gp,  255.97690_gp, 2400.69180_gp, 1830.85230_gp,&
         1486.06200_gp, 1412.68180_gp, 1281.33490_gp, 1007.85710_gp, 1087.98370_gp,&
         851.06340_gp,  883.88900_gp,  919.53540_gp,  704.37890_gp,  701.35490_gp,&
         839.28280_gp,  711.68230_gp,  586.00840_gp,  514.50770_gp,  440.29990_gp,&
         375.13390_gp, 2673.93440_gp, 2188.28640_gp, 1858.63770_gp, 1638.49380_gp,&
         1475.25000_gp,   46.04060_gp,   29.16170_gp,  890.51260_gp,  463.30650_gp,&
         295.35920_gp,  191.95280_gp,  130.60450_gp,   97.06750_gp,   72.41860_gp,&
         55.04870_gp, 1058.38990_gp,  751.66280_gp,  673.11930_gp,  510.23520_gp,&
         385.58950_gp,  313.86570_gp,  251.77230_gp,  202.85720_gp, 1765.89460_gp,&
         1365.48370_gp, 1112.00980_gp, 1060.62550_gp,  963.94340_gp,  758.99550_gp,&
         820.93130_gp,  643.12470_gp,  670.34900_gp,  695.83480_gp,  533.65350_gp,&
         534.52070_gp,  637.92330_gp,  545.78450_gp,  453.41860_gp,  400.49980_gp,&
         344.98260_gp,  295.79940_gp, 1969.47960_gp, 1631.54240_gp, 1394.09760_gp/)
    vdwparams%coeffs(  901: 1000)=(/ 1233.90950_gp, 1114.09240_gp,  845.89720_gp,   51.05270_gp,   32.05880_gp,&
         1011.60990_gp,  521.42490_gp,  330.32570_gp,  213.47990_gp,  144.56970_gp,&
         107.05880_gp,   79.61210_gp,   60.35210_gp, 1201.52340_gp,  847.34970_gp,&
         756.61950_gp,  571.05010_gp,  429.76770_gp,  348.77440_gp,  278.91230_gp,&
         224.07940_gp, 2007.61370_gp, 1543.42580_gp, 1255.34460_gp, 1195.84340_gp,&
         1086.02310_gp,  854.74740_gp,  923.87430_gp,  723.32020_gp,  752.96780_gp,&
         782.24830_gp,  599.61920_gp,  599.31410_gp,  715.98110_gp,  610.53270_gp,&
         505.50910_gp,  445.47780_gp,  382.75910_gp,  327.38270_gp, 2237.96650_gp,&
         1844.37740_gp, 1572.37590_gp, 1389.59790_gp, 1253.34650_gp,  949.64900_gp,&
         1067.01690_gp,   39.57450_gp,   25.46220_gp,  730.55330_gp,  387.26810_gp,&
         249.93820_gp,  164.12790_gp,  112.63780_gp,   84.27250_gp,   63.25990_gp,&
         48.34290_gp,  869.37740_gp,  626.21040_gp,  564.07160_gp,  431.21320_gp,&
         328.40720_gp,  268.78240_gp,  216.80700_gp,  175.58760_gp, 1447.30480_gp,&
         1131.17690_gp,  923.58110_gp,  883.14920_gp,  803.88070_gp,  633.43940_gp,&
         686.17720_gp,  538.14350_gp,  562.48000_gp,  582.90140_gp,  447.42590_gp,&
         450.15430_gp,  536.10930_gp,  461.67950_gp,  385.98830_gp,  342.38140_gp,&
         296.26680_gp,  255.15160_gp, 1615.83870_gp, 1351.19880_gp, 1159.90710_gp,&
         1029.76650_gp,  931.74450_gp,  710.28370_gp,  796.19510_gp,  598.19880_gp,&
         43.16610_gp,   27.57450_gp,  794.20060_gp,  424.41200_gp,  273.71950_gp,&
         179.27080_gp,  122.64260_gp,   91.49590_gp,   68.48880_gp,   52.20720_gp/)
    vdwparams%coeffs( 1001: 1100)=(/  945.42210_gp,  685.66280_gp,  617.92440_gp,  472.37380_gp,  359.25470_gp,&
         293.52140_gp,  236.27500_gp,  190.93520_gp, 1567.49540_gp, 1234.60130_gp,&
         1009.19330_gp,  965.23080_gp,  878.80270_gp,  691.63140_gp,  750.32510_gp,&
         587.75440_gp,  615.43640_gp,  637.76340_gp,  488.75010_gp,  492.63190_gp,&
         586.77320_gp,  505.44520_gp,  422.19710_gp,  374.06320_gp,  323.19410_gp,&
         277.86800_gp, 1749.76150_gp, 1473.54870_gp, 1266.84360_gp, 1125.34030_gp,&
         1018.29150_gp,  775.83530_gp,  869.89620_gp,  653.17750_gp,  713.94270_gp,&
         40.23720_gp,   25.87160_gp,  722.42940_gp,  390.43100_gp,  253.34250_gp,&
         166.71110_gp,  114.47630_gp,   85.64270_gp,   64.27070_gp,   49.09900_gp,&
         860.62090_gp,  629.59120_gp,  569.14100_gp,  436.94630_gp,  333.52190_gp,&
         273.15680_gp,  220.40990_gp,  178.50210_gp, 1424.02710_gp, 1129.66700_gp,&
         924.85560_gp,  885.76570_gp,  807.12580_gp,  635.32570_gp,  689.96470_gp,&
         540.65990_gp,  567.10440_gp,  587.18190_gp,  450.05590_gp,  454.80310_gp,&
         541.15010_gp,  467.70420_gp,  391.85280_gp,  347.83930_gp,  301.13900_gp,&
         259.39640_gp, 1590.44150_gp, 1347.90070_gp, 1161.92890_gp, 1033.87350_gp,&
         936.55280_gp,  714.94950_gp,  801.04420_gp,  602.79560_gp,  658.87530_gp,&
         608.50410_gp,   33.54130_gp,   21.78960_gp,  614.14770_gp,  325.57750_gp,&
         210.68110_gp,  138.87580_gp,   95.71370_gp,   71.89260_gp,   54.19030_gp,&
         41.57800_gp,  731.08320_gp,  526.47760_gp,  474.57130_gp,  363.35850_gp,&
         277.35770_gp,  227.50240_gp,  183.98760_gp,  149.42690_gp, 1218.46990_gp/)
    vdwparams%coeffs( 1101: 1200)=(/  951.49350_gp,  776.89140_gp,  743.20410_gp,  676.65530_gp,  533.68780_gp,&
         577.81270_gp,  453.63270_gp,  473.92440_gp,  490.97830_gp,  377.35310_gp,&
         379.56550_gp,  451.58870_gp,  389.30750_gp,  326.03710_gp,  289.64690_gp,&
         251.11890_gp,  216.73360_gp, 1360.75170_gp, 1136.90990_gp,  976.20250_gp,&
         867.00560_gp,  784.81880_gp,  599.05300_gp,  671.17170_gp,  505.02690_gp,&
         551.08020_gp,  508.77420_gp,  426.74500_gp,   35.68790_gp,   23.15630_gp,&
         621.05360_gp,  340.41220_gp,  222.62730_gp,  147.41380_gp,  101.73770_gp,&
         76.40920_gp,   57.55000_gp,   44.10620_gp,  740.57870_gp,  547.67210_gp,&
         497.03370_gp,  383.67900_gp,  294.24650_gp,  241.76140_gp,  195.70440_gp,&
         158.96500_gp, 1222.55530_gp,  978.42510_gp,  802.60200_gp,  770.03030_gp,&
         702.42130_gp,  553.09190_gp,  601.40530_gp,  471.53290_gp,  495.63130_gp,&
         512.62200_gp,  393.05100_gp,  398.45720_gp,  473.42030_gp,  410.90580_gp,&
         345.60850_gp,  307.55390_gp,  266.97000_gp,  230.54930_gp, 1366.39380_gp,&
         1167.04710_gp, 1009.44180_gp,  900.11540_gp,  816.54950_gp,  624.95060_gp,&
         699.53820_gp,  527.93870_gp,  577.01070_gp,  533.41340_gp,  445.87010_gp,&
         468.19000_gp,   44.04170_gp,   27.95860_gp,  828.96730_gp,  438.33960_gp,&
         281.17370_gp,  183.33680_gp,  124.97980_gp,   92.99680_gp,   69.45580_gp,&
         52.85170_gp,  986.09990_gp,  709.35650_gp,  637.58320_gp,  485.53920_gp,&
         368.02750_gp,  299.98530_gp,  240.91040_gp,  194.26460_gp, 1638.81530_gp,&
         1281.26120_gp, 1045.84010_gp,  999.09710_gp,  908.97530_gp,  715.27990_gp/)
    vdwparams%coeffs( 1201: 1300)=(/  775.27360_gp,  607.09640_gp,  634.72630_gp,  658.25260_gp,  504.38020_gp,&
         507.24520_gp,  604.71830_gp,  519.35010_gp,  432.61670_gp,  382.59680_gp,&
         329.92500_gp,  283.12690_gp, 1828.56960_gp, 1529.68620_gp, 1311.89220_gp,&
         1163.59630_gp, 1051.89170_gp,  799.99940_gp,  897.64510_gp,  672.65160_gp,&
         735.25250_gp,  678.09960_gp,  567.33030_gp,  593.35740_gp,  757.73970_gp,&
         41.53340_gp,   26.78620_gp,  708.45970_gp,  393.79910_gp,  258.71460_gp,&
         171.43670_gp,  118.15090_gp,   88.54080_gp,   66.49960_gp,   50.81100_gp,&
         845.13480_gp,  631.94340_gp,  575.40180_gp,  445.76230_gp,  342.41690_gp,&
         281.32550_gp,  227.56350_gp,  184.58450_gp, 1390.29120_gp, 1123.06160_gp,&
         923.01820_gp,  886.63180_gp,  809.41340_gp,  636.65610_gp,  693.74420_gp,&
         543.36290_gp,  572.81690_gp,  592.04640_gp,  453.20890_gp,  461.19710_gp,&
         548.05490_gp,  477.10340_gp,  401.94150_gp,  357.78710_gp,  310.50930_gp,&
         267.95320_gp, 1554.36980_gp, 1338.56650_gp, 1161.36740_gp, 1037.28120_gp,&
         941.77920_gp,  721.29060_gp,  807.23580_gp,  609.61520_gp,  666.89540_gp,&
         616.75220_gp,  514.44770_gp,  541.54900_gp,  685.61380_gp,  627.56770_gp,&
         37.76810_gp,   24.85830_gp,  586.68880_gp,  340.30300_gp,  229.04210_gp,&
         154.48510_gp,  107.84360_gp,   81.53580_gp,   61.69120_gp,   47.40730_gp,&
         701.69030_gp,  541.98260_gp,  499.72920_gp,  393.68990_gp,  306.69120_gp,&
         254.24180_gp,  207.40430_gp,  169.45330_gp, 1147.67000_gp,  950.65430_gp,&
         785.83230_gp,  758.94060_gp,  695.06830_gp,  547.19810_gp,  598.51680_gp/)
    vdwparams%coeffs( 1301: 1400)=(/  469.45690_gp,  498.04720_gp,  513.05850_gp,  393.02490_gp,  403.86210_gp,&
         478.39930_gp,  421.86160_gp,  359.55640_gp,  322.36080_gp,  281.80500_gp,&
         244.79020_gp, 1286.17500_gp, 1132.17130_gp,  992.40010_gp,  892.08770_gp,&
         813.43450_gp,  627.68490_gp,  700.51630_gp,  533.38040_gp,  583.33460_gp,&
         540.91250_gp,  450.64270_gp,  476.57830_gp,  598.28910_gp,  553.35920_gp,&
         492.93790_gp,   35.48410_gp,   23.69110_gp,  521.25110_gp,  309.54390_gp,&
         211.48980_gp,  144.30320_gp,  101.61240_gp,   77.29700_gp,   58.78810_gp,&
         45.36130_gp,  624.40760_gp,  490.87600_gp,  456.04000_gp,  362.95180_gp,&
         285.27790_gp,  237.90180_gp,  195.18610_gp,  160.26740_gp, 1019.31210_gp,&
         855.18820_gp,  709.13030_gp,  687.12350_gp,  630.49010_gp,  496.91430_gp,&
         544.41560_gp,  427.64390_gp,  455.08360_gp,  467.84930_gp,  358.81670_gp,&
         370.59360_gp,  438.11020_gp,  389.29970_gp,  334.22340_gp,  301.05450_gp,&
         264.45160_gp,  230.74260_gp, 1144.19130_gp, 1018.34130_gp,  897.80430_gp,&
         810.10280_gp,  740.62110_gp,  574.29320_gp,  639.74650_gp,  489.69560_gp,&
         535.21890_gp,  497.09590_gp,  414.15120_gp,  438.88610_gp,  548.12670_gp,&
         509.96980_gp,  457.10570_gp,  425.53550_gp,   32.51710_gp,   22.08130_gp,&
         449.85120_gp,  273.73980_gp,  190.10940_gp,  131.40510_gp,   93.45620_gp,&
         71.60590_gp,   54.79900_gp,   42.49540_gp,  539.85360_gp,  432.17920_gp,&
         404.76310_gp,  325.69120_gp,  258.51210_gp,  217.03370_gp,  179.23060_gp,&
         148.02160_gp,  880.17410_gp,  747.93060_gp,  622.18410_gp,  605.05080_gp/)
    vdwparams%coeffs( 1401: 1500)=(/  556.31660_gp,  439.16600_gp,  481.80770_gp,  379.22120_gp,  404.69290_gp,&
         415.13030_gp,  318.96420_gp,  331.06510_gp,  390.45420_gp,  349.77600_gp,&
         302.67430_gp,  274.06300_gp,  242.05560_gp,  212.28200_gp,  989.89210_gp,&
         890.67230_gp,  790.02570_gp,  715.73180_gp,  656.22760_gp,  511.67720_gp,&
         568.79090_gp,  437.99920_gp,  478.24810_gp,  444.96620_gp,  370.93440_gp,&
         393.76550_gp,  488.97090_gp,  457.71560_gp,  412.98520_gp,  386.09970_gp,&
         351.96670_gp,   29.60550_gp,   20.45940_gp,  386.99690_gp,  240.76890_gp,&
         169.85780_gp,  118.91780_gp,   85.42930_gp,   65.94030_gp,   50.79000_gp,&
         39.59510_gp,  465.28170_gp,  378.61340_gp,  357.29600_gp,  290.49850_gp,&
         232.77860_gp,  196.73740_gp,  163.53680_gp,  135.85920_gp,  758.32530_gp,&
         651.53920_gp,  543.56310_gp,  530.43120_gp,  488.65220_gp,  386.50670_gp,&
         424.41250_gp,  334.82590_gp,  358.09540_gp,  366.55970_gp,  282.28880_gp,&
         294.21380_gp,  346.13450_gp,  312.42650_gp,  272.41100_gp,  247.92750_gp,&
         220.15810_gp,  194.06990_gp,  854.52230_gp,  776.08590_gp,  692.23860_gp,&
         629.52300_gp,  578.78760_gp,  453.79980_gp,  503.38110_gp,  389.95320_gp,&
         425.26840_gp,  396.34680_gp,  330.76460_gp,  351.52850_gp,  434.09320_gp,&
         408.59520_gp,  370.95720_gp,  348.22390_gp,  318.89050_gp,  290.22230_gp,&
         105.00960_gp,   61.31840_gp, 3242.24040_gp, 1316.21890_gp,  754.92430_gp,&
         455.69600_gp,  294.14830_gp,  210.86720_gp,  152.56720_gp,  113.17640_gp,&
         3798.94390_gp, 2220.96520_gp, 1881.71640_gp, 1317.00560_gp,  936.30690_gp/)
    vdwparams%coeffs( 1501: 1600)=(/  735.00610_gp,  569.67010_gp,  445.81710_gp, 6808.39000_gp, 4387.16690_gp,&
         3450.83690_gp, 3215.85440_gp, 2877.52790_gp, 2283.06740_gp, 2396.34450_gp,&
         1884.28240_gp, 1875.72090_gp, 1976.86090_gp, 1533.40630_gp, 1443.00040_gp,&
         1755.51050_gp, 1406.22060_gp, 1108.31450_gp,  950.19750_gp,  794.21840_gp,&
         662.95560_gp, 7560.43610_gp, 5305.16270_gp, 4289.60970_gp, 3674.85560_gp,&
         3256.36330_gp, 2398.75200_gp, 2726.04860_gp, 1968.30660_gp, 2129.92750_gp,&
         1935.83920_gp, 1657.78290_gp, 1663.04470_gp, 2226.48670_gp, 1891.92070_gp,&
         1566.50420_gp, 1394.57090_gp, 1207.56720_gp, 1043.39550_gp, 9330.72940_gp,&
         99.45790_gp,   59.11410_gp, 2551.52600_gp, 1160.44190_gp,  690.68520_gp,&
         426.36670_gp,  279.07420_gp,  201.74110_gp,  146.91740_gp,  109.50000_gp,&
         3009.72330_gp, 1926.35040_gp, 1665.76620_gp, 1201.20210_gp,  871.16470_gp,&
         690.86510_gp,  540.31600_gp,  425.79090_gp, 5171.29780_gp, 3659.16550_gp,&
         2924.43950_gp, 2748.97970_gp, 2475.14270_gp, 1952.11260_gp, 2079.49500_gp,&
         1628.40400_gp, 1656.93840_gp, 1736.39870_gp, 1335.73190_gp, 1292.66300_gp,&
         1559.70120_gp, 1281.91520_gp, 1028.54240_gp,  889.61470_gp,  749.84750_gp,&
         630.24360_gp, 5743.56730_gp, 4394.36070_gp, 3639.47350_gp, 3160.17760_gp,&
         2819.51470_gp, 2098.48130_gp, 2373.91310_gp, 1735.48270_gp, 1890.13710_gp,&
         1727.27380_gp, 1460.58580_gp, 1493.73690_gp, 1964.09130_gp, 1710.88320_gp,&
         1442.51350_gp, 1295.63540_gp, 1131.72210_gp,  985.12660_gp, 7016.49170_gp,&
         5726.98870_gp,   89.10180_gp,   53.79820_gp, 2085.06610_gp,  994.03160_gp/)
    vdwparams%coeffs( 1601: 1700)=(/  604.49150_gp,  378.79280_gp,  250.55770_gp,  182.40750_gp,  133.61230_gp,&
         100.02770_gp, 2465.59650_gp, 1637.15860_gp, 1431.74960_gp, 1049.19340_gp,&
         770.49240_gp,  615.55590_gp,  484.76910_gp,  384.23060_gp, 4185.89680_gp,&
         3059.89410_gp, 2461.76190_gp, 2324.74540_gp, 2099.58530_gp, 1653.44350_gp,&
         1771.76530_gp, 1386.31320_gp, 1423.40230_gp, 1487.14540_gp, 1141.40780_gp,&
         1118.39730_gp, 1344.79680_gp, 1120.00210_gp,  908.36700_gp,  790.48550_gp,&
         670.38430_gp,  566.51010_gp, 4653.74540_gp, 3666.65560_gp, 3069.30370_gp,&
         2682.08790_gp, 2401.93260_gp, 1798.19700_gp, 2029.61630_gp, 1493.88950_gp,&
         1629.91260_gp, 1493.48550_gp, 1257.42970_gp, 1295.93640_gp, 1689.40850_gp,&
         1489.04400_gp, 1268.26220_gp, 1145.55950_gp, 1006.51040_gp,  880.82020_gp,&
         5675.74340_gp, 4754.48260_gp, 3990.61720_gp,   44.51040_gp,   29.30810_gp,&
         704.33450_gp,  403.74270_gp,  270.61870_gp,  182.21640_gp,  127.13550_gp,&
         96.11890_gp,   72.73590_gp,   55.90480_gp,  841.89860_gp,  644.22200_gp,&
         592.40470_gp,  465.26310_gp,  361.80150_gp,  299.74470_gp,  244.44680_gp,&
         199.71380_gp, 1381.01830_gp, 1134.67480_gp,  936.40680_gp,  903.33060_gp,&
         826.70470_gp,  651.14400_gp,  711.14230_gp,  558.01770_gp,  590.73870_gp,&
         608.92420_gp,  466.81100_gp,  478.33500_gp,  566.87830_gp,  498.62840_gp,&
         424.29060_gp,  380.16060_gp,  332.18840_gp,  288.50100_gp, 1547.07800_gp,&
         1352.08030_gp, 1182.01240_gp, 1060.97210_gp,  966.62010_gp,  745.11760_gp,&
         831.83600_gp,  632.66390_gp,  691.56460_gp,  640.96110_gp,  534.59930_gp/)
    vdwparams%coeffs( 1701: 1800)=(/  564.39100_gp,  709.50090_gp,  654.67230_gp,  582.16180_gp,  539.41670_gp,&
         487.03330_gp,  437.28340_gp, 1884.45540_gp, 1724.74150_gp, 1512.25040_gp,&
         688.03530_gp,   88.68770_gp,   53.15230_gp, 2210.31680_gp, 1018.95400_gp,&
         610.76150_gp,  379.10380_gp,  249.26330_gp,  180.83830_gp,  132.15840_gp,&
         98.82200_gp, 2609.04400_gp, 1687.54190_gp, 1464.44460_gp, 1061.52130_gp,&
         773.17210_gp,  614.86520_gp,  482.25090_gp,  381.05030_gp, 4473.77840_gp,&
         3190.20860_gp, 2554.73530_gp, 2405.05610_gp, 2167.63250_gp, 1709.11370_gp,&
         1823.78580_gp, 1428.10390_gp, 1457.03580_gp, 1525.41410_gp, 1172.95490_gp,&
         1139.43520_gp, 1372.86490_gp, 1133.16270_gp,  912.52630_gp,  791.02140_gp,&
         668.33280_gp,  563.01560_gp, 4971.07780_gp, 3828.96260_gp, 3181.44050_gp,&
         2767.98940_gp, 2472.69050_gp, 1844.07150_gp, 2084.59750_gp, 1527.60560_gp,&
         1664.46400_gp, 1522.48060_gp, 1286.04410_gp, 1318.27910_gp, 1728.21620_gp,&
         1511.04590_gp, 1278.34060_gp, 1150.42000_gp, 1007.01980_gp,  878.37980_gp,&
         6073.89350_gp, 4983.31260_gp, 4150.32030_gp, 1527.14610_gp, 4342.23860_gp,&
         85.40800_gp,   51.29520_gp, 2091.50680_gp,  974.61670_gp,  586.18960_gp,&
         364.65190_gp,  240.11390_gp,  174.36780_gp,  127.53080_gp,   95.42020_gp,&
         2470.60200_gp, 1611.78310_gp, 1401.31910_gp, 1018.50980_gp,  743.24070_gp,&
         591.67520_gp,  464.50620_gp,  367.31630_gp, 4216.39980_gp, 3036.36620_gp,&
         2435.08890_gp, 2294.28570_gp, 2068.98090_gp, 1630.60180_gp, 1742.19420_gp,&
         1363.84560_gp, 1394.08220_gp, 1458.73990_gp, 1120.99310_gp, 1091.57220_gp/)
    vdwparams%coeffs( 1801: 1900)=(/ 1314.23060_gp, 1087.25330_gp,  877.00660_gp,  760.89490_gp,  643.43050_gp,&
         542.43830_gp, 4685.08460_gp, 3642.04920_gp, 3032.91830_gp, 2642.05260_gp,&
         2361.69740_gp, 1763.11390_gp, 1992.15340_gp, 1461.64580_gp, 1593.36960_gp,&
         1458.16710_gp, 1230.49690_gp, 1263.34820_gp, 1653.50390_gp, 1448.87110_gp,&
         1227.74550_gp, 1105.81430_gp,  968.79470_gp,  845.67030_gp, 5716.79990_gp,&
         4733.41850_gp, 3951.99560_gp, 1465.96950_gp, 4125.90230_gp, 3924.42110_gp,&
         83.33100_gp,   50.07480_gp, 2031.27400_gp,  949.12550_gp,  571.41130_gp,&
         355.67120_gp,  234.29130_gp,  170.18170_gp,  124.49460_gp,   93.16350_gp,&
         2399.82070_gp, 1568.98850_gp, 1364.84940_gp,  992.75070_gp,  724.81970_gp,&
         577.16960_gp,  453.23130_gp,  358.47210_gp, 4092.05750_gp, 2952.97940_gp,&
         2369.12440_gp, 2232.65030_gp, 2013.71940_gp, 1586.85120_gp, 1696.04410_gp,&
         1327.60640_gp, 1357.73950_gp, 1420.50450_gp, 1091.42140_gp, 1063.48980_gp,&
         1280.17620_gp, 1059.75360_gp,  855.21770_gp,  742.16490_gp,  627.73420_gp,&
         529.30700_gp, 4547.07200_gp, 3541.49300_gp, 2950.93030_gp, 2571.49640_gp,&
         2299.04810_gp, 1716.80380_gp, 1939.62150_gp, 1423.55310_gp, 1552.05220_gp,&
         1420.54690_gp, 1198.41940_gp, 1230.95980_gp, 1610.41420_gp, 1411.95880_gp,&
         1197.01570_gp, 1078.38540_gp,  944.98140_gp,  825.04730_gp, 5547.49000_gp,&
         4601.21710_gp, 3844.03670_gp, 1429.06900_gp, 4011.36350_gp, 3816.11180_gp,&
         3710.93750_gp,   81.41230_gp,   48.94060_gp, 1976.69550_gp,  925.81580_gp,&
         557.83800_gp,  347.39210_gp,  228.90650_gp,  166.30040_gp,  121.67270_gp/)
    vdwparams%coeffs( 1901: 2000)=(/   91.06140_gp, 2335.64740_gp, 1529.91200_gp, 1331.47860_gp,  969.10260_gp,&
         707.86040_gp,  563.78900_gp,  442.80950_gp,  350.28120_gp, 3979.63090_gp,&
         2877.06800_gp, 2308.99620_gp, 2176.41840_gp, 1963.27230_gp, 1546.91950_gp,&
         1653.87730_gp, 1294.49510_gp, 1324.47700_gp, 1385.53130_gp, 1064.37980_gp,&
         1037.74900_gp, 1248.99150_gp, 1034.50240_gp,  835.16220_gp,  724.89870_gp,&
         613.24050_gp,  517.16170_gp, 4422.26960_gp, 3449.99060_gp, 2876.17160_gp,&
         2507.08050_gp, 2241.80850_gp, 1674.43600_gp, 1891.58830_gp, 1388.66680_gp,&
         1514.20160_gp, 1386.06220_gp, 1169.03600_gp, 1201.24730_gp, 1570.96430_gp,&
         1378.08300_gp, 1168.75240_gp, 1053.12550_gp,  923.01940_gp,  806.00020_gp,&
         5394.48460_gp, 4481.04490_gp, 3745.69210_gp, 1395.14710_gp, 3907.19230_gp,&
         3717.54070_gp, 3615.20670_gp, 3522.05080_gp,   79.71380_gp,   47.93650_gp,&
         1928.08100_gp,  905.11010_gp,  545.80150_gp,  340.05800_gp,  224.13830_gp,&
         162.86380_gp,  119.17340_gp,   89.19880_gp, 2278.48980_gp, 1495.17920_gp,&
         1301.84380_gp,  948.12880_gp,  692.83380_gp,  551.93900_gp,  433.58310_gp,&
         343.03070_gp, 3879.47290_gp, 2809.52330_gp, 2255.51570_gp, 2126.41900_gp,&
         1918.42540_gp, 1511.41560_gp, 1616.40230_gp, 1265.06430_gp, 1294.93130_gp,&
         1354.45950_gp, 1040.34900_gp, 1014.89590_gp, 1221.30300_gp, 1012.10530_gp,&
         817.38880_gp,  709.60390_gp,  600.40640_gp,  506.40980_gp, 4311.09610_gp,&
         3368.56460_gp, 2809.68850_gp, 2449.81880_gp, 2190.93870_gp, 1636.79560_gp,&
         1848.90970_gp, 1357.68100_gp, 1480.58780_gp, 1355.44250_gp, 1142.93520_gp/)
    vdwparams%coeffs( 2001: 2100)=(/ 1174.86990_gp, 1535.92540_gp, 1348.02160_gp, 1143.69190_gp, 1030.73800_gp,&
         903.56300_gp,  789.13200_gp, 5258.20680_gp, 4374.08260_gp, 3658.21200_gp,&
         1365.06150_gp, 3814.49210_gp, 3629.82680_gp, 3530.02280_gp, 3439.16010_gp,&
         3358.31220_gp,   64.18800_gp,   39.57860_gp, 1407.73350_gp,  689.92870_gp,&
         426.48970_gp,  270.91830_gp,  181.25810_gp,  133.14380_gp,   98.35730_gp,&
         74.19120_gp, 1667.38950_gp, 1130.81980_gp,  996.79240_gp,  739.02000_gp,&
         548.36330_gp,  441.23030_gp,  350.03690_gp,  279.35150_gp, 2817.63330_gp,&
         2094.53260_gp, 1691.70520_gp, 1602.89270_gp, 1450.71740_gp, 1142.70150_gp,&
         1228.04260_gp,  961.56000_gp,  992.06730_gp, 1034.19320_gp,  793.87630_gp,&
         783.51140_gp,  939.36520_gp,  789.63690_gp,  645.98490_gp,  565.29910_gp,&
         482.31720_gp,  409.98030_gp, 3135.92420_gp, 2507.75890_gp, 2113.20760_gp,&
         1854.48370_gp, 1665.43720_gp, 1253.07300_gp, 1411.67590_gp, 1045.06170_gp,&
         1140.45250_gp, 1047.12790_gp,  880.53130_gp,  911.07960_gp, 1179.97200_gp,&
         1048.09760_gp,  899.51620_gp,  816.30470_gp,  720.98470_gp,  634.21630_gp,&
         3823.84260_gp, 3243.65610_gp, 2740.54780_gp, 1070.95970_gp, 2838.13570_gp,&
         2705.68260_gp, 2632.66570_gp, 2566.06070_gp, 2506.83520_gp, 1891.67190_gp,&
         70.17940_gp,   42.36790_gp, 1802.38520_gp,  808.48660_gp,  482.79770_gp,&
         299.84690_gp,  197.62340_gp,  143.77540_gp,  105.39700_gp,   79.04840_gp,&
         2124.08310_gp, 1343.36510_gp, 1162.15600_gp,  839.17760_gp,  610.83180_gp,&
         486.23190_gp,  381.95840_gp,  302.41490_gp, 3687.63050_gp, 2564.14430_gp/)
    vdwparams%coeffs( 2101: 2200)=(/ 2045.49270_gp, 1922.83020_gp, 1730.98900_gp, 1368.00490_gp, 1454.16140_gp,&
         1140.87130_gp, 1157.90110_gp, 1213.20140_gp,  935.75440_gp,  903.59910_gp,&
         1090.09390_gp,  896.46510_gp,  721.20420_gp,  625.44800_gp,  528.93810_gp,&
         446.20840_gp, 4098.56780_gp, 3083.98990_gp, 2548.59110_gp, 2211.49070_gp,&
         1973.56480_gp, 1470.52650_gp, 1663.13470_gp, 1217.35110_gp, 1323.52230_gp,&
         1209.68470_gp, 1025.68350_gp, 1046.57890_gp, 1375.60670_gp, 1197.55970_gp,&
         1011.29750_gp,  909.91270_gp,  796.69950_gp,  695.43500_gp, 5022.69670_gp,&
         4030.09100_gp, 3335.74350_gp, 1209.62870_gp, 3507.17160_gp, 3325.07970_gp,&
         3231.32410_gp, 3146.17400_gp, 3070.38530_gp, 2276.25320_gp, 2851.66770_gp,&
         67.95800_gp,   41.16220_gp, 1720.73530_gp,  775.77180_gp,  465.24390_gp,&
         289.84550_gp,  191.45500_gp,  139.49370_gp,  102.38110_gp,   76.85570_gp,&
         2027.85680_gp, 1287.24710_gp, 1116.03410_gp,  808.33220_gp,  589.89410_gp,&
         470.29770_gp,  369.98330_gp,  293.29050_gp, 3526.97050_gp, 2451.77340_gp,&
         1957.45810_gp, 1841.58830_gp, 1658.70190_gp, 1310.75380_gp, 1394.48530_gp,&
         1094.00160_gp, 1111.81740_gp, 1164.27550_gp,  897.83730_gp,  868.75770_gp,&
         1047.56400_gp,  863.58600_gp,  696.27600_gp,  604.60480_gp,  511.97160_gp,&
         432.38930_gp, 3921.69620_gp, 2948.56880_gp, 2440.17620_gp, 2119.49080_gp,&
         1892.78580_gp, 1411.73630_gp, 1596.20740_gp, 1169.67850_gp, 1271.89370_gp,&
         1163.03740_gp,  985.56970_gp, 1006.83820_gp, 1321.52570_gp, 1152.79750_gp,&
         975.44350_gp,  878.67530_gp,  770.29080_gp,  673.13940_gp, 4812.33230_gp/)
    vdwparams%coeffs( 2201: 2300)=(/ 3852.29490_gp, 3192.40390_gp, 1166.16330_gp, 3355.87650_gp, 3179.96290_gp,&
         3090.40720_gp, 3009.07780_gp, 2936.70650_gp, 2181.06750_gp, 2729.97390_gp,&
         2617.33100_gp,   72.22030_gp,   43.57110_gp, 1704.96790_gp,  811.36880_gp,&
         491.88750_gp,  307.51980_gp,  203.15540_gp,  147.83550_gp,  108.30760_gp,&
         81.14160_gp, 2016.33920_gp, 1337.47830_gp, 1167.93580_gp,  854.07400_gp,&
         625.93940_gp,  499.45760_gp,  392.93850_gp,  311.24870_gp, 3419.51600_gp,&
         2501.41640_gp, 2012.03520_gp, 1899.20440_gp, 1714.86350_gp, 1350.29320_gp,&
         1446.60200_gp, 1131.77100_gp, 1161.49050_gp, 1213.93480_gp,  931.67930_gp,&
         912.00820_gp, 1096.44500_gp,  911.72230_gp,  738.21290_gp,  641.74520_gp,&
         543.71620_gp,  459.12050_gp, 3800.83040_gp, 2997.02120_gp, 2507.35720_gp,&
         2190.02720_gp, 1960.51940_gp, 1466.80170_gp, 1655.94660_gp, 1218.06510_gp,&
         1329.14260_gp, 1217.65820_gp, 1025.39020_gp, 1056.37180_gp, 1377.89010_gp,&
         1213.10430_gp, 1031.80710_gp,  931.12390_gp,  817.32200_gp,  714.64770_gp,&
         4632.84280_gp, 3885.50930_gp, 3259.91320_gp, 1230.59890_gp, 3391.54110_gp,&
         3229.81680_gp, 3141.60980_gp, 3061.25370_gp, 2989.76600_gp, 2237.90260_gp,&
         2724.36020_gp, 2606.51360_gp, 2664.16680_gp,   70.71540_gp,   42.66580_gp,&
         1664.77000_gp,  793.68770_gp,  481.43850_gp,  301.06920_gp,  198.91850_gp,&
         144.75870_gp,  106.05520_gp,   79.45390_gp, 1968.99280_gp, 1307.98140_gp,&
         1142.56560_gp,  835.89690_gp,  612.78090_gp,  489.01070_gp,  384.74960_gp,&
         304.77480_gp, 3337.22270_gp, 2444.67580_gp, 1966.90500_gp, 1856.87430_gp/)
    vdwparams%coeffs( 2301: 2400)=(/ 1676.81290_gp, 1320.19120_gp, 1414.70480_gp, 1106.72460_gp, 1136.18970_gp,&
         1187.38730_gp,  911.17070_gp,  892.33320_gp, 1072.67680_gp,  892.30170_gp,&
         722.66540_gp,  628.29340_gp,  532.36410_gp,  449.55920_gp, 3709.42970_gp,&
         2928.71640_gp, 2451.17750_gp, 2141.42020_gp, 1917.22030_gp, 1434.61300_gp,&
         1619.51860_gp, 1191.47160_gp, 1300.26370_gp, 1191.29620_gp, 1002.97970_gp,&
         1033.60090_gp, 1347.84970_gp, 1187.10400_gp, 1009.95830_gp,  911.51370_gp,&
         800.19180_gp,  699.72310_gp, 4520.96450_gp, 3796.09240_gp, 3186.23270_gp,&
         1204.41850_gp, 3313.88790_gp, 3156.19920_gp, 3070.08420_gp, 2991.62750_gp,&
         2921.83100_gp, 2187.78030_gp, 2661.14610_gp, 2546.12490_gp, 2603.96530_gp,&
         2545.17130_gp,   69.37910_gp,   41.87270_gp, 1627.09510_gp,  777.46020_gp,&
         471.98600_gp,  295.30140_gp,  195.16300_gp,  142.04780_gp,  104.08000_gp,&
         77.97880_gp, 1924.65930_gp, 1280.79110_gp, 1119.33910_gp,  819.42800_gp,&
         600.97010_gp,  479.69090_gp,  377.48730_gp,  299.06220_gp, 3259.93680_gp,&
         2391.99540_gp, 1925.12630_gp, 1817.78810_gp, 1641.73540_gp, 1292.43840_gp,&
         1385.37120_gp, 1083.69570_gp, 1113.02590_gp, 1163.03690_gp,  892.35160_gp,&
         874.39610_gp, 1050.97050_gp,  874.71320_gp,  708.69500_gp,  616.26580_gp,&
         522.26550_gp,  441.09320_gp, 3623.65110_gp, 2865.26130_gp, 2399.24870_gp,&
         2096.64000_gp, 1877.41740_gp, 1405.13520_gp, 1586.11000_gp, 1167.18780_gp,&
         1273.90260_gp, 1167.27070_gp,  982.52100_gp, 1012.89040_gp, 1320.38870_gp,&
         1163.49880_gp,  990.25880_gp,  893.90800_gp,  784.88410_gp,  686.44520_gp/)
    vdwparams%coeffs( 2401: 2500)=(/ 4415.96940_gp, 3712.88120_gp, 3118.00010_gp, 1180.77760_gp, 3241.72660_gp,&
         3087.84890_gp, 3003.69220_gp, 2927.01180_gp, 2858.79770_gp, 2141.52640_gp,&
         2602.31220_gp, 2489.96430_gp, 2548.18020_gp, 2490.69890_gp, 2437.45390_gp,&
         68.49490_gp,   41.28660_gp, 1612.75810_gp,  769.39830_gp,  466.60990_gp,&
         291.68170_gp,  192.63240_gp,  140.13000_gp,  102.62510_gp,   76.85810_gp,&
         1907.53320_gp, 1267.88660_gp, 1107.52150_gp,  810.18450_gp,  593.79110_gp,&
         473.73900_gp,  372.62670_gp,  295.08360_gp, 3231.71170_gp, 2369.09580_gp,&
         1906.28020_gp, 1799.63770_gp, 1625.14300_gp, 1279.34930_gp, 1371.11980_gp,&
         1072.49590_gp, 1101.21820_gp, 1150.85530_gp,  882.99190_gp,  864.85460_gp,&
         1039.66010_gp,  864.79960_gp,  700.27180_gp,  608.71740_gp,  515.66490_gp,&
         435.35350_gp, 3592.04650_gp, 2837.94230_gp, 2375.46790_gp, 2075.34170_gp,&
         1858.03530_gp, 1390.20590_gp, 1569.44260_gp, 1154.52090_gp, 1260.07640_gp,&
         1154.46390_gp,  971.80930_gp, 1001.62140_gp, 1306.20230_gp, 1150.45710_gp,&
         978.68710_gp,  883.19190_gp,  775.20910_gp,  677.75440_gp, 4377.47660_gp,&
         3677.96060_gp, 3087.51540_gp, 1167.07990_gp, 3210.84890_gp, 3058.22770_gp,&
         2974.82380_gp, 2898.83490_gp, 2831.23440_gp, 2119.97940_gp, 2577.80050_gp,&
         2466.35830_gp, 2523.36300_gp, 2466.41430_gp, 2413.65140_gp, 2390.12270_gp,&
         58.65190_gp,   35.95160_gp, 1294.75830_gp,  634.58810_gp,  391.46520_gp,&
         247.95570_gp,  165.41840_gp,  121.21480_gp,   89.34050_gp,   67.25750_gp,&
         1533.26750_gp, 1040.18250_gp,  916.38240_gp,  678.55410_gp,  502.55170_gp/)
    vdwparams%coeffs( 2501: 2600)=(/  403.67180_gp,  319.62860_gp,  254.59790_gp, 2592.12990_gp, 1925.54530_gp,&
         1555.39390_gp, 1473.40910_gp, 1333.40270_gp, 1049.61420_gp, 1128.53090_gp,&
         883.02170_gp,  911.41500_gp,  950.32780_gp,  728.89260_gp,  719.52410_gp,&
         862.91390_gp,  724.75130_gp,  592.05870_gp,  517.45760_gp,  440.84280_gp,&
         374.14020_gp, 2884.75850_gp, 2304.85580_gp, 1942.13870_gp, 1704.04400_gp,&
         1529.97210_gp, 1150.16040_gp, 1296.27010_gp,  958.69920_gp, 1046.75970_gp,&
         960.91280_gp,  807.44060_gp,  835.83100_gp, 1083.30080_gp,  961.84020_gp,&
         824.69640_gp,  747.75790_gp,  659.70490_gp,  579.59110_gp, 3519.52360_gp,&
         2980.64360_gp, 2518.21920_gp,  981.75310_gp, 2608.91930_gp, 2486.19310_gp,&
         2419.11640_gp, 2357.94390_gp, 2303.55400_gp, 1737.65650_gp, 2091.93800_gp,&
         2005.92430_gp, 2056.53100_gp, 2010.52030_gp, 1968.04740_gp, 1948.31020_gp,&
         1597.47960_gp,   58.85040_gp,   36.57490_gp, 1177.22680_gp,  607.58590_gp,&
         383.50840_gp,  246.68260_gp,  166.24720_gp,  122.59130_gp,   90.77840_gp,&
         68.55050_gp, 1398.12540_gp,  987.52390_gp,  880.79020_gp,  663.33140_gp,&
         497.70150_gp,  402.79470_gp,  321.11450_gp,  257.16510_gp, 2330.05620_gp,&
         1797.10240_gp, 1461.95890_gp, 1391.98480_gp, 1263.86600_gp,  993.62820_gp,&
         1074.72870_gp,  840.45370_gp,  875.43360_gp,  909.84070_gp,  696.42830_gp,&
         696.21350_gp,  832.32490_gp,  708.69710_gp,  585.42500_gp,  514.87840_gp,&
         441.33770_gp,  376.52870_gp, 2596.40330_gp, 2146.40400_gp, 1829.75150_gp,&
         1616.47370_gp, 1457.24920_gp, 1102.56690_gp, 1239.52010_gp,  923.41980_gp/)
    vdwparams%coeffs( 2601: 2700)=(/ 1009.75270_gp,  929.46760_gp,  777.79950_gp,  811.20400_gp, 1042.29930_gp,&
         936.63610_gp,  811.48690_gp,  740.06690_gp,  656.81710_gp,  580.13770_gp,&
         3160.64460_gp, 2760.90940_gp, 2360.91860_gp,  963.47730_gp, 2424.04910_gp,&
         2317.05140_gp, 2256.05840_gp, 2200.29060_gp, 2150.74170_gp, 1640.65670_gp,&
         1930.99490_gp, 1852.95810_gp, 1926.59250_gp, 1884.30480_gp, 1845.50500_gp,&
         1826.21490_gp, 1507.24410_gp, 1441.23940_gp,   54.40150_gp,   34.28340_gp,&
         1027.35770_gp,  544.43550_gp,  348.62850_gp,  226.71590_gp,  154.06960_gp,&
         114.29670_gp,   85.08140_gp,   64.52380_gp, 1222.03350_gp,  880.84690_gp,&
         791.41350_gp,  602.13470_gp,  455.69700_gp,  370.87230_gp,  297.28390_gp,&
         239.22910_gp, 2026.66980_gp, 1589.48830_gp, 1297.78260_gp, 1239.50330_gp,&
         1127.58120_gp,  886.58800_gp,  961.51940_gp,  752.31040_gp,  787.03270_gp,&
         816.34550_gp,  624.84080_gp,  628.67410_gp,  749.97760_gp,  643.73010_gp,&
         535.60150_gp,  473.16850_gp,  407.46910_gp,  349.13050_gp, 2260.84000_gp,&
         1896.97990_gp, 1627.26700_gp, 1443.19730_gp, 1304.32130_gp,  991.16570_gp,&
         1112.45110_gp,  832.81910_gp,  910.84260_gp,  839.83510_gp,  701.92840_gp,&
         734.56590_gp,  938.77020_gp,  849.25790_gp,  740.51010_gp,  677.94300_gp,&
         604.19180_gp,  535.77730_gp, 2751.65090_gp, 2434.42070_gp, 2094.85310_gp,&
         878.10980_gp, 2141.76670_gp, 2049.61250_gp, 1996.28790_gp, 1947.46690_gp,&
         1904.10850_gp, 1462.31310_gp, 1702.67050_gp, 1635.57520_gp, 1708.49120_gp,&
         1671.31100_gp, 1637.32250_gp, 1619.77100_gp, 1342.86970_gp, 1292.58360_gp/)
    vdwparams%coeffs( 2701: 2800)=(/ 1163.82410_gp,   46.48130_gp,   29.82060_gp,  842.41230_gp,  451.98610_gp,&
         292.90150_gp,  192.59840_gp,  132.14100_gp,   98.76660_gp,   74.03090_gp,&
         56.47670_gp, 1002.95080_gp,  729.33520_gp,  658.74710_gp,  505.18900_gp,&
         385.41090_gp,  315.56320_gp,  254.52870_gp,  206.02860_gp, 1664.43210_gp,&
         1311.76410_gp, 1072.80500_gp, 1026.92140_gp,  935.38590_gp,  736.50630_gp,&
         799.17580_gp,  626.33800_gp,  656.22600_gp,  679.63140_gp,  521.08220_gp,&
         525.88940_gp,  626.22030_gp,  540.69840_gp,  452.78370_gp,  401.84290_gp,&
         347.79630_gp,  299.48800_gp, 1858.76830_gp, 1565.94130_gp, 1347.81630_gp,&
         1198.32020_gp, 1085.08780_gp,  827.86250_gp,  927.74070_gp,  697.63460_gp,&
         762.26600_gp,  703.76190_gp,  588.70290_gp,  616.64100_gp,  784.70960_gp,&
         712.94940_gp,  624.90100_gp,  574.15840_gp,  513.83150_gp,  457.58530_gp,&
         2263.90130_gp, 2008.77510_gp, 1733.91290_gp,  740.68630_gp, 1769.41390_gp,&
         1693.89880_gp, 1650.00040_gp, 1609.77980_gp, 1574.06610_gp, 1214.44480_gp,&
         1407.83900_gp, 1353.34330_gp, 1413.28210_gp, 1382.57610_gp, 1354.57790_gp,&
         1339.77110_gp, 1114.30560_gp, 1076.38070_gp,  971.95240_gp,  814.36220_gp,&
         47.48140_gp,   30.54920_gp,  840.30290_gp,  456.86870_gp,  297.69360_gp,&
         196.40730_gp,  135.04060_gp,  101.06690_gp,   75.83370_gp,   57.89700_gp,&
         1001.18330_gp,  735.65530_gp,  666.53080_gp,  513.21450_gp,  392.66040_gp,&
         322.00220_gp,  260.08050_gp,  210.74970_gp, 1656.38290_gp, 1317.29610_gp,&
         1079.31680_gp, 1034.54890_gp,  943.14230_gp,  742.38390_gp,  806.78450_gp/)
    vdwparams%coeffs( 2801: 2900)=(/  632.20880_gp,  663.88710_gp,  687.01430_gp,  526.48590_gp,  532.99620_gp,&
         634.16620_gp,  549.32810_gp,  461.15310_gp,  409.81080_gp,  355.13780_gp,&
         306.13060_gp, 1850.61310_gp, 1571.68630_gp, 1356.85330_gp, 1208.47370_gp,&
         1095.42820_gp,  837.07080_gp,  937.53010_gp,  706.23580_gp,  771.95670_gp,&
         713.18460_gp,  595.99800_gp,  625.41990_gp,  794.23010_gp,  723.67260_gp,&
         635.81610_gp,  584.92700_gp,  524.12750_gp,  467.26420_gp, 2253.50370_gp,&
         2013.47400_gp, 1743.37910_gp,  753.13690_gp, 1775.33570_gp, 1700.57060_gp,&
         1656.77680_gp, 1616.62660_gp, 1580.98350_gp, 1223.42240_gp, 1410.59040_gp,&
         1356.71150_gp, 1420.70810_gp, 1389.99340_gp, 1362.03340_gp, 1346.99510_gp,&
         1122.57060_gp, 1087.69310_gp,  983.75910_gp,  824.96150_gp,  836.33100_gp,&
         43.50270_gp,   28.35460_gp,  731.92340_gp,  407.16460_gp,  268.73040_gp,&
         179.05630_gp,  124.05100_gp,   93.35970_gp,   70.39540_gp,   53.96420_gp,&
         873.35890_gp,  653.07500_gp,  595.56280_gp,  462.68400_gp,  356.70880_gp,&
         294.00210_gp,  238.64150_gp,  194.23140_gp, 1439.66200_gp, 1161.26180_gp,&
         954.45330_gp,  917.46900_gp,  837.84550_gp,  659.82040_gp,  718.50590_gp,&
         563.50590_gp,  593.75330_gp,  613.35130_gp,  470.25530_gp,  478.53750_gp,&
         568.22510_gp,  495.61460_gp,  418.68920_gp,  373.55900_gp,  325.06940_gp,&
         281.30190_gp, 1610.37540_gp, 1384.82540_gp, 1202.10440_gp, 1074.34910_gp,&
         976.08420_gp,  748.94180_gp,  837.52750_gp,  633.78860_gp,  692.67540_gp,&
         640.89680_gp,  535.28690_gp,  563.12050_gp,  711.69730_gp,  652.11350_gp/)
    vdwparams%coeffs( 2901: 3000)=(/  576.10750_gp,  531.78230_gp,  478.25570_gp,  427.87240_gp, 1961.01920_gp,&
         1770.79840_gp, 1541.67470_gp,  681.78680_gp, 1564.30820_gp, 1499.87320_gp,&
         1461.63150_gp, 1426.52810_gp, 1395.37630_gp, 1086.42420_gp, 1241.32050_gp,&
         1195.07830_gp, 1255.67280_gp, 1228.71110_gp, 1204.25230_gp, 1190.65280_gp,&
         996.36440_gp,  970.80130_gp,  881.09380_gp,  740.95590_gp,  752.15710_gp,&
         678.52780_gp,   40.19280_gp,   26.49360_gp,  649.22250_gp,  367.68190_gp,&
         245.23810_gp,  164.75860_gp,  114.88740_gp,   86.87930_gp,   65.79040_gp,&
         50.61590_gp,  775.63370_gp,  587.94870_gp,  538.99920_gp,  421.78390_gp,&
         327.23040_gp,  270.85510_gp,  220.77630_gp,  180.37110_gp, 1275.51980_gp,&
         1039.89930_gp,  856.78860_gp,  825.51250_gp,  754.91930_gp,  594.87280_gp,&
         648.70850_gp,  509.21810_gp,  537.89460_gp,  554.84790_gp,  425.67880_gp,&
         434.87920_gp,  515.50510_gp,  452.12560_gp,  383.92650_gp,  343.68330_gp,&
         300.11780_gp,  260.56720_gp, 1428.24930_gp, 1239.73510_gp, 1080.83400_gp,&
         968.64880_gp,  881.71010_gp,  678.87180_gp,  758.18340_gp,  575.94520_gp,&
         629.29090_gp,  582.96220_gp,  486.81310_gp,  513.03760_gp,  645.86290_gp,&
         594.40490_gp,  527.46810_gp,  488.23190_gp,  440.42470_gp,  395.18890_gp,&
         1739.49920_gp, 1583.13330_gp, 1384.22250_gp,  623.81140_gp, 1400.68130_gp,&
         1343.95140_gp, 1309.94540_gp, 1278.69860_gp, 1250.97670_gp,  978.87040_gp,&
         1110.69690_gp, 1070.19120_gp, 1126.93880_gp, 1102.86140_gp, 1081.08110_gp,&
         1068.64600_gp,  897.29150_gp,  878.05750_gp,  799.16860_gp,  673.70750_gp/)
    vdwparams%coeffs( 3001: 3100)=(/  684.58550_gp,  619.11480_gp,  566.06600_gp,   33.73740_gp,   22.65500_gp,&
         523.29220_gp,  300.16240_gp,  202.54080_gp,  137.57620_gp,   96.86580_gp,&
         73.81830_gp,   56.30530_gp,   43.59300_gp,  625.96250_gp,  478.85350_gp,&
         441.20140_gp,  347.88170_gp,  271.98720_gp,  226.47100_gp,  185.74640_gp,&
         152.65670_gp, 1029.98990_gp,  844.43870_gp,  696.92720_gp,  673.05270_gp,&
         616.30660_gp,  486.52160_gp,  530.65270_gp,  417.43140_gp,  441.40210_gp,&
         454.62950_gp,  349.58030_gp,  358.01260_gp,  423.42920_gp,  373.42980_gp,&
         319.04090_gp,  286.87030_gp,  251.74820_gp,  219.65520_gp, 1154.76380_gp,&
         1007.06230_gp,  880.99930_gp,  791.55440_gp,  721.93400_gp,  558.23890_gp,&
         622.42220_gp,  475.07150_gp,  518.45950_gp,  480.93890_gp,  402.19150_gp,&
         424.03900_gp,  531.44420_gp,  491.01090_gp,  437.82090_gp,  406.61850_gp,&
         368.25920_gp,  331.79040_gp, 1407.17240_gp, 1285.43220_gp, 1127.51570_gp,&
         517.68570_gp, 1138.89940_gp, 1093.24620_gp, 1065.70500_gp, 1040.37420_gp,&
         1017.90320_gp,  800.35340_gp,  903.95940_gp,  871.58450_gp,  917.60410_gp,&
         898.02720_gp,  880.36430_gp,  870.03670_gp,  732.91850_gp,  719.57520_gp,&
         656.82120_gp,  555.62970_gp,  565.02640_gp,  512.42690_gp,  469.66640_gp,&
         391.14480_gp,   31.64310_gp,   21.40570_gp,  482.71200_gp,  278.48740_gp,&
         188.78590_gp,  128.78430_gp,   91.02100_gp,   69.57720_gp,   53.22580_gp,&
         41.31660_gp,  577.75500_gp,  443.84910_gp,  409.78640_gp,  324.09250_gp,&
         254.14960_gp,  212.10390_gp,  174.38080_gp,  143.64860_gp,  950.62180_gp/)
    vdwparams%coeffs( 3101: 3200)=(/  781.60440_gp,  645.57500_gp,  624.06930_gp,  571.77460_gp,  451.67990_gp,&
         492.72680_gp,  387.91960_gp,  410.40680_gp,  422.44780_gp,  325.11980_gp,&
         333.31770_gp,  393.82270_gp,  348.09220_gp,  298.10130_gp,  268.50000_gp,&
         236.07650_gp,  206.37560_gp, 1066.32890_gp,  932.20550_gp,  816.74650_gp,&
         734.61520_gp,  670.55170_gp,  519.41000_gp,  578.74610_gp,  442.59210_gp,&
         482.80890_gp,  448.12160_gp,  374.94860_gp,  395.41490_gp,  494.65360_gp,&
         457.75290_gp,  408.94270_gp,  380.29110_gp,  344.94080_gp,  311.27000_gp,&
         1299.59750_gp, 1189.54050_gp, 1044.89960_gp,  483.48310_gp, 1054.63390_gp,&
         1012.57350_gp,  987.12160_gp,  963.70210_gp,  942.92740_gp,  742.91920_gp,&
         837.29210_gp,  807.54210_gp,  850.29540_gp,  832.17180_gp,  815.83760_gp,&
         806.19030_gp,  680.07310_gp,  668.62000_gp,  611.03860_gp,  517.60830_gp,&
         526.53700_gp,  478.06720_gp,  438.60830_gp,  365.82470_gp,  342.35260_gp,&
         32.64750_gp,   22.08370_gp,  487.40200_gp,  285.06170_gp,  194.16150_gp,&
         132.73170_gp,   93.88360_gp,   71.77600_gp,   54.90230_gp,   42.60850_gp,&
         583.79600_gp,  453.35960_gp,  419.83840_gp,  333.21670_gp,  261.85840_gp,&
         218.71850_gp,  179.91240_gp,  148.23250_gp,  957.36170_gp,  794.65850_gp,&
         657.59170_gp,  636.52100_gp,  583.66570_gp,  460.85090_gp,  503.55400_gp,&
         396.29990_gp,  420.25520_gp,  432.27410_gp,  332.43400_gp,  341.87100_gp,&
         403.70590_gp,  357.84550_gp,  307.03710_gp,  276.76970_gp,  243.49600_gp,&
         212.93560_gp, 1074.41280_gp,  947.20110_gp,  832.42220_gp,  749.98240_gp/)
    vdwparams%coeffs( 3201: 3300)=(/  685.23850_gp,  531.45620_gp,  591.92160_gp,  453.28020_gp,  494.71860_gp,&
         459.42820_gp,  383.95180_gp,  405.66070_gp,  506.64510_gp,  470.08500_gp,&
         420.79320_gp,  391.67400_gp,  355.54970_gp,  321.02510_gp, 1309.29210_gp,&
         1207.03410_gp, 1063.62760_gp,  497.15300_gp, 1071.20930_gp, 1029.08350_gp,&
         1003.38460_gp,  979.72320_gp,  958.73970_gp,  757.56310_gp,  849.22660_gp,&
         819.52610_gp,  865.29630_gp,  846.94850_gp,  830.44190_gp,  820.53560_gp,&
         693.54790_gp,  683.91950_gp,  625.92010_gp,  530.51690_gp,  540.05040_gp,&
         490.85310_gp,  450.68210_gp,  376.00960_gp,  351.93500_gp,  362.07550_gp,&
         45.56590_gp,   29.38320_gp,  841.56200_gp,  445.09900_gp,  287.43460_gp,&
         188.91370_gp,  129.77590_gp,   97.18870_gp,   73.03530_gp,   55.87820_gp,&
         1001.34990_gp,  719.84280_gp,  648.50370_gp,  495.88920_gp,  377.86060_gp,&
         309.40570_gp,  249.71520_gp,  202.36390_gp, 1669.97320_gp, 1301.23250_gp,&
         1062.12570_gp, 1015.66760_gp,  924.50540_gp,  728.75260_gp,  789.17040_gp,&
         619.12140_gp,  646.88450_gp,  670.34360_gp,  514.78030_gp,  517.77130_gp,&
         616.53920_gp,  531.01480_gp,  444.13160_gp,  394.08930_gp,  341.15420_gp,&
         293.94640_gp, 1864.61640_gp, 1554.70360_gp, 1334.14320_gp, 1184.35580_gp,&
         1071.68580_gp,  817.13880_gp,  915.94930_gp,  688.33610_gp,  751.41170_gp,&
         693.49420_gp,  581.27400_gp,  607.45600_gp,  773.88330_gp,  701.33540_gp,&
         613.79650_gp,  563.66330_gp,  504.31850_gp,  449.15730_gp, 2272.40780_gp,&
         1997.76650_gp, 1718.79300_gp,  728.02070_gp, 1758.47310_gp, 1682.12810_gp/)
    vdwparams%coeffs( 3301: 3400)=(/ 1638.21730_gp, 1598.00830_gp, 1562.29390_gp, 1202.52290_gp, 1402.66270_gp,&
         1347.90370_gp, 1401.38320_gp, 1370.74490_gp, 1342.76540_gp, 1328.19070_gp,&
         1103.09320_gp, 1061.99600_gp,  957.83300_gp,  802.62850_gp,  812.49710_gp,&
         729.24020_gp,  662.79090_gp,  546.91720_gp,  509.61460_gp,  521.90910_gp,&
         792.23780_gp,   45.17860_gp,   29.39630_gp,  770.49850_gp,  425.98770_gp,&
         280.17140_gp,  186.22350_gp,  128.82430_gp,   96.88010_gp,   73.02850_gp,&
         55.99090_gp,  919.12570_gp,  684.10110_gp,  622.71340_gp,  482.58810_gp,&
         371.29390_gp,  305.63890_gp,  247.81410_gp,  201.53300_gp, 1516.41930_gp,&
         1218.70090_gp, 1000.86020_gp,  961.37180_gp,  877.56550_gp,  691.08340_gp,&
         752.11560_gp,  589.80850_gp,  620.87440_gp,  641.67640_gp,  491.99760_gp,&
         499.93220_gp,  593.73850_gp,  516.89320_gp,  435.92810_gp,  388.53730_gp,&
         337.76830_gp,  292.04720_gp, 1695.69090_gp, 1453.47450_gp, 1259.86910_gp,&
         1124.96170_gp, 1021.47330_gp,  782.99060_gp,  875.95820_gp,  662.16750_gp,&
         723.70460_gp,  669.38820_gp,  559.27100_gp,  587.93680_gp,  743.87970_gp,&
         680.57520_gp,  600.37140_gp,  553.68260_gp,  497.48720_gp,  444.70070_gp,&
         2064.66870_gp, 1859.40850_gp, 1616.46770_gp,  710.62890_gp, 1641.98850_gp,&
         1573.97630_gp, 1533.74690_gp, 1496.82910_gp, 1464.06240_gp, 1138.09790_gp,&
         1303.43160_gp, 1254.53350_gp, 1317.01770_gp, 1288.69070_gp, 1262.96620_gp,&
         1248.78910_gp, 1043.90340_gp, 1015.48870_gp,  920.83050_gp,  773.85260_gp,&
         785.29450_gp,  707.90090_gp,  645.54700_gp,  534.02070_gp,  498.13010_gp/)
    vdwparams%coeffs( 3401: 3500)=(/  511.31450_gp,  761.96480_gp,  738.81560_gp,   42.24190_gp,   27.93290_gp,&
         661.52780_gp,  381.03070_gp,  256.05370_gp,  172.78310_gp,  120.79930_gp,&
         91.49340_gp,   69.36890_gp,   53.42040_gp,  791.08680_gp,  607.60080_gp,&
         559.43360_gp,  440.13430_gp,  342.77880_gp,  284.29480_gp,  232.12190_gp,&
         189.87560_gp, 1296.62530_gp, 1068.63310_gp,  882.51530_gp,  851.88550_gp,&
         779.93360_gp,  614.45090_gp,  671.30550_gp,  526.93210_gp,  558.17400_gp,&
         575.14870_gp,  441.05510_gp,  452.36780_gp,  535.70130_gp,  471.83760_gp,&
         401.99200_gp,  360.47350_gp,  315.27740_gp,  274.07550_gp, 1452.94180_gp,&
         1273.23700_gp, 1114.39770_gp, 1001.03150_gp,  912.48200_gp,  704.07130_gp,&
         785.75170_gp,  598.27920_gp,  653.93370_gp,  606.31280_gp,  505.75740_gp,&
         534.17140_gp,  670.73540_gp,  619.58270_gp,  551.56050_gp,  511.39760_gp,&
         462.08340_gp,  415.20110_gp, 1769.63900_gp, 1623.47400_gp, 1425.12230_gp,&
         651.71260_gp, 1438.26150_gp, 1380.94510_gp, 1346.26740_gp, 1314.37840_gp,&
         1286.09630_gp, 1010.38600_gp, 1139.01660_gp, 1098.31100_gp, 1159.77740_gp,&
         1135.14240_gp, 1112.91090_gp, 1099.94100_gp,  926.12910_gp,  909.81540_gp,&
         829.82850_gp,  700.43950_gp,  712.41700_gp,  645.37650_gp,  590.85050_gp,&
         490.72340_gp,  458.47260_gp,  471.52310_gp,  688.56550_gp,  672.69040_gp,&
         617.52960_gp,   40.71060_gp,   27.21180_gp,  607.31510_gp,  357.46030_gp,&
         243.30570_gp,  165.73270_gp,  116.66150_gp,   88.77360_gp,   67.56720_gp,&
         52.18810_gp,  727.26310_gp,  567.80720_gp,  526.26000_gp,  417.68440_gp/)
    vdwparams%coeffs( 3501: 3600)=(/  327.71360_gp,  273.09950_gp,  223.98500_gp,  183.92230_gp, 1189.19020_gp,&
         992.21250_gp,  821.78250_gp,  795.53840_gp,  729.56570_gp,  575.17310_gp,&
         629.48100_gp,  494.59140_gp,  525.49790_gp,  540.53270_gp,  414.78690_gp,&
         427.46160_gp,  505.36000_gp,  448.08220_gp,  384.08530_gp,  345.72660_gp,&
         303.54480_gp,  264.79780_gp, 1334.36920_gp, 1181.86040_gp, 1039.86220_gp,&
         937.20360_gp,  856.23980_gp,  663.36080_gp,  739.19570_gp,  565.31410_gp,&
         617.69980_gp,  573.50990_gp,  478.23360_gp,  506.16860_gp,  632.78330_gp,&
         587.60970_gp,  525.87280_gp,  489.16400_gp,  443.53000_gp,  399.83350_gp,&
         1626.04190_gp, 1504.76070_gp, 1327.73970_gp,  620.85480_gp, 1335.43090_gp,&
         1283.23360_gp, 1251.28790_gp, 1221.87750_gp, 1195.80460_gp,  945.15700_gp,&
         1056.98640_gp, 1020.27900_gp, 1079.68890_gp, 1056.88310_gp, 1036.37880_gp,&
         1024.03380_gp,  865.80300_gp,  855.04390_gp,  782.48020_gp,  662.38200_gp,&
         674.49800_gp,  612.77720_gp,  562.31000_gp,  468.25250_gp,  437.92710_gp,&
         450.80740_gp,  650.74290_gp,  638.22510_gp,  588.61880_gp,  562.60110_gp,&
         38.04680_gp,   25.80710_gp,  535.65210_gp,  323.15400_gp,  223.39710_gp,&
         153.98830_gp,  109.36180_gp,   83.74370_gp,   64.07970_gp,   49.70430_gp,&
         642.55990_gp,  511.04720_gp,  477.38160_gp,  382.89140_gp,  303.17170_gp,&
         254.18990_gp,  209.69510_gp,  173.06680_gp, 1048.67540_gp,  886.79380_gp,&
         736.87680_gp,  715.82420_gp,  657.77000_gp,  519.22800_gp,  569.18500_gp,&
         447.94120_gp,  477.40460_gp,  490.02960_gp,  376.54500_gp,  390.05450_gp/)
    vdwparams%coeffs( 3601: 3700)=(/  460.13230_gp,  411.19140_gp,  355.08990_gp,  321.15720_gp,  283.36260_gp,&
         248.31690_gp, 1178.78430_gp, 1056.17810_gp,  934.92800_gp,  845.95890_gp,&
         775.00770_gp,  603.51860_gp,  671.20920_gp,  516.18070_gp,  563.61340_gp,&
         524.17160_gp,  437.16280_gp,  463.63160_gp,  576.48690_gp,  538.57130_gp,&
         485.03660_gp,  452.97260_gp,  412.49290_gp,  373.39600_gp, 1437.52370_gp,&
         1342.75570_gp, 1191.80040_gp,  572.19880_gp, 1194.15520_gp, 1148.48790_gp,&
         1120.16740_gp, 1094.05800_gp, 1070.92190_gp,  852.63060_gp,  945.04500_gp,&
         913.33710_gp,  968.25870_gp,  947.92020_gp,  929.71630_gp,  918.34420_gp,&
         780.31510_gp,  775.23370_gp,  712.28880_gp,  605.24130_gp,  617.10680_gp,&
         562.58800_gp,  517.72850_gp,  432.64880_gp,  405.18020_gp,  417.47370_gp,&
         594.34970_gp,  585.43140_gp,  542.87440_gp,  520.58720_gp,  483.65360_gp,&
         35.27500_gp,   24.28860_gp,  471.29940_gp,  290.40790_gp,  203.66210_gp,&
         141.99460_gp,  101.73400_gp,   78.39930_gp,   60.32080_gp,   46.99700_gp,&
         566.30470_gp,  457.52070_gp,  430.40700_gp,  348.53240_gp,  278.34410_gp,&
         234.75630_gp,  194.77780_gp,  161.57780_gp,  923.48960_gp,  789.46510_gp,&
         657.82090_gp,  641.06890_gp,  590.13650_gp,  466.57860_gp,  512.00690_gp,&
         403.71620_gp,  431.25430_gp,  441.80520_gp,  340.10290_gp,  353.75690_gp,&
         416.40720_gp,  374.73780_gp,  325.84680_gp,  296.05910_gp,  262.46780_gp,&
         231.03920_gp, 1039.87880_gp,  940.37660_gp,  836.82440_gp,  759.86860_gp,&
         697.90750_gp,  546.17620_gp,  606.28810_gp,  468.75270_gp,  511.33700_gp/)
    vdwparams%coeffs( 3701: 3800)=(/  476.29130_gp,  397.51230_gp,  422.13930_gp,  522.25240_gp,  490.45740_gp,&
         444.23510_gp,  416.40380_gp,  380.74630_gp,  346.03400_gp, 1269.15140_gp,&
         1194.25530_gp, 1065.38860_gp,  523.79990_gp, 1064.11130_gp, 1024.16010_gp,&
         999.10040_gp,  975.96710_gp,  955.47530_gp,  765.77540_gp,  842.48080_gp,&
         815.08100_gp,  864.86000_gp,  846.76410_gp,  830.63370_gp,  820.22230_gp,&
         700.09390_gp,  699.08500_gp,  644.67130_gp,  549.83920_gp,  561.21370_gp,&
         513.27610_gp,  473.60870_gp,  397.19060_gp,  372.48250_gp,  384.02410_gp,&
         539.90550_gp,  533.69840_gp,  497.29720_gp,  478.31130_gp,  446.04470_gp,&
         412.82750_gp,  100.74940_gp,   59.65940_gp, 2894.59660_gp, 1218.93480_gp,&
         710.65850_gp,  434.19670_gp,  282.77630_gp,  203.97230_gp,  148.35840_gp,&
         110.51200_gp, 3398.64610_gp, 2044.62960_gp, 1746.54180_gp, 1237.82690_gp,&
         888.82250_gp,  701.97610_gp,  547.27120_gp,  430.44680_gp, 6020.91230_gp,&
         3988.55350_gp, 3153.39330_gp, 2948.32140_gp, 2644.13480_gp, 2095.15630_gp,&
         2209.25160_gp, 1735.97670_gp, 1740.60130_gp, 1830.22860_gp, 1416.87370_gp,&
         1346.49750_gp, 1633.08690_gp, 1321.96690_gp, 1050.90640_gp,  905.47240_gp,&
         760.70860_gp,  637.92100_gp, 6687.95510_gp, 4814.06130_gp, 3923.86300_gp,&
         3377.73480_gp, 3001.20090_gp, 2220.84040_gp, 2518.94570_gp, 1828.60370_gp,&
         1981.98220_gp, 1805.27920_gp, 1540.30200_gp, 1555.09640_gp, 2067.17310_gp,&
         1773.26200_gp, 1480.12550_gp, 1323.63290_gp, 1151.63860_gp,  999.49930_gp,&
         8231.52800_gp, 6339.72350_gp, 5171.56280_gp, 1776.59190_gp, 5496.96010_gp/)
    vdwparams%coeffs( 3801: 3900)=(/ 5188.68260_gp, 5037.68300_gp, 4900.97390_gp, 4779.23850_gp, 3500.46500_gp,&
         4518.01450_gp, 4325.68820_gp, 4221.62160_gp, 4121.10130_gp, 4027.06140_gp,&
         3990.87780_gp, 3218.77040_gp, 2921.81190_gp, 2556.32950_gp, 2107.69910_gp,&
         2103.24540_gp, 1838.38820_gp, 1636.38840_gp, 1327.11810_gp, 1227.09140_gp,&
         1239.49860_gp, 2109.47530_gp, 1933.30250_gp, 1669.99650_gp, 1540.83170_gp,&
         1368.73630_gp, 1213.49590_gp, 7314.73980_gp,   99.41310_gp,   59.67400_gp,&
         2434.06970_gp, 1132.03440_gp,  681.37890_gp,  424.16130_gp,  279.36020_gp,&
         202.82600_gp,  148.25160_gp,  110.81260_gp, 2874.57170_gp, 1871.96000_gp,&
         1627.94930_gp, 1183.71030_gp,  864.34290_gp,  688.37560_gp,  540.58680_gp,&
         427.51470_gp, 4912.85810_gp, 3528.95900_gp, 2829.22730_gp, 2665.59500_gp,&
         2403.68780_gp, 1894.77070_gp, 2023.89190_gp, 1584.58380_gp, 1619.22210_gp,&
         1694.24690_gp, 1302.21360_gp, 1267.79660_gp, 1527.01490_gp, 1263.57300_gp,&
         1019.72510_gp,  885.04330_gp,  748.65020_gp,  631.27570_gp, 5459.40390_gp,&
         4233.97360_gp, 3524.61800_gp, 3070.02060_gp, 2744.26750_gp, 2048.83500_gp,&
         2314.93490_gp, 1698.46720_gp, 1851.10770_gp, 1693.93620_gp, 1429.69500_gp,&
         1467.46240_gp, 1921.02150_gp, 1683.25370_gp, 1426.78310_gp, 1285.46080_gp,&
         1126.53190_gp,  983.63370_gp, 6664.60490_gp, 5505.25820_gp, 4594.30950_gp,&
         1703.71040_gp, 4797.60610_gp, 4562.38650_gp, 4436.25920_gp, 4321.49140_gp,&
         4219.36900_gp, 3144.71840_gp, 3869.34750_gp, 3700.69780_gp, 3753.65340_gp,&
         3667.97460_gp, 3588.43510_gp, 3554.00670_gp, 2889.29830_gp, 2692.12290_gp/)
    vdwparams%coeffs( 3901: 4000)=(/ 2381.33180_gp, 1968.47300_gp, 1976.06720_gp, 1742.85220_gp, 1561.68870_gp,&
         1270.42580_gp, 1176.62220_gp, 1195.75730_gp, 1954.72920_gp, 1828.65150_gp,&
         1604.63490_gp, 1491.40890_gp, 1335.16320_gp, 1190.92860_gp, 6045.71680_gp,&
         5305.43990_gp,   89.40600_gp,   54.46130_gp, 2008.25910_gp,  977.21480_gp,&
         600.05330_gp,  378.66720_gp,  251.80860_gp,  184.03180_gp,  135.27920_gp,&
         101.58570_gp, 2377.46170_gp, 1604.03000_gp, 1409.87110_gp, 1040.56810_gp,&
         768.51320_gp,  616.13190_gp,  486.88350_gp,  387.07200_gp, 4017.47820_gp,&
         2977.14780_gp, 2402.15400_gp, 2273.27390_gp, 2055.95890_gp, 1618.36340_gp,&
         1738.44950_gp, 1360.06630_gp, 1401.76220_gp, 1462.53570_gp, 1121.74660_gp,&
         1104.95230_gp, 1326.36010_gp, 1111.09580_gp,  905.53420_gp,  790.26450_gp,&
         672.16580_gp,  569.54790_gp, 4469.00530_gp, 3564.35100_gp, 2997.77220_gp,&
         2627.04460_gp, 2356.70830_gp, 1769.25380_gp, 1994.87290_gp, 1473.04570_gp,&
         1608.17220_gp, 1475.39240_gp, 1240.22740_gp, 1282.28350_gp, 1665.05020_gp,&
         1475.11210_gp, 1262.09170_gp, 1142.90900_gp, 1006.91410_gp,  883.41400_gp,&
         5447.89290_gp, 4612.02470_gp, 3889.70520_gp, 1503.15690_gp, 4032.39920_gp,&
         3843.35360_gp, 3739.36460_gp, 3644.54690_gp, 3560.22780_gp, 2679.68020_gp,&
         3233.36810_gp, 3096.80780_gp, 3176.94060_gp, 3105.68470_gp, 3039.84610_gp,&
         3009.59990_gp, 2462.37740_gp, 2320.22770_gp, 2064.56700_gp, 1711.55700_gp,&
         1723.22310_gp, 1527.65350_gp, 1374.36400_gp, 1121.39240_gp, 1040.02380_gp,&
         1060.06070_gp, 1694.66050_gp, 1600.86160_gp, 1417.43460_gp, 1323.63600_gp/)
    vdwparams%coeffs( 4001: 4100)=(/ 1191.32720_gp, 1067.50650_gp, 4981.49610_gp, 4466.95840_gp, 3799.65650_gp,&
         82.08360_gp,   51.11030_gp, 1672.47680_gp,  852.52490_gp,  536.08580_gp,&
         344.35760_gp,  232.07230_gp,  171.24180_gp,  126.93930_gp,   95.98450_gp,&
         1985.21030_gp, 1388.37560_gp, 1235.26810_gp,  927.51070_gp,  694.79830_gp,&
         562.06360_gp,  448.04860_gp,  358.92730_gp, 3320.81710_gp, 2537.73080_gp,&
         2060.89010_gp, 1960.28270_gp, 1778.65460_gp, 1399.47780_gp, 1511.08160_gp,&
         1182.50930_gp, 1228.73800_gp, 1277.78960_gp,  979.19840_gp,  975.88560_gp,&
         1167.05830_gp,  991.20150_gp,  817.55520_gp,  718.66790_gp,  615.85050_gp,&
         525.42400_gp, 3699.66370_gp, 3033.12730_gp, 2578.70740_gp, 2274.83250_gp,&
         2049.22990_gp, 1549.17800_gp, 1742.12020_gp, 1296.65510_gp, 1416.84230_gp,&
         1303.60600_gp, 1092.59390_gp, 1137.17450_gp, 1463.15800_gp, 1311.58800_gp,&
         1134.43850_gp, 1033.87690_gp,  917.09200_gp,  809.80180_gp, 4505.56810_gp,&
         3907.12620_gp, 3331.54730_gp, 1347.82230_gp, 3427.78940_gp, 3274.47920_gp,&
         3187.74280_gp, 3108.47690_gp, 3038.03450_gp, 2312.20860_gp, 2735.81130_gp,&
         2624.19220_gp, 2719.13190_gp, 2659.12680_gp, 2604.00090_gp, 2576.97670_gp,&
         2123.68190_gp, 2024.82390_gp, 1813.75640_gp, 1510.13610_gp, 1524.92970_gp,&
         1359.85760_gp, 1229.21640_gp, 1007.42250_gp,  936.09580_gp,  956.80470_gp,&
         1491.56980_gp, 1422.94310_gp, 1272.66010_gp, 1195.14920_gp, 1082.87190_gp,&
         976.14170_gp, 4155.59320_gp, 3804.77500_gp, 3270.42580_gp, 2847.27040_gp,&
         80.72600_gp,   50.00520_gp, 1730.60380_gp,  857.63790_gp,  533.08180_gp/)
    vdwparams%coeffs( 4101: 4200)=(/  339.98630_gp,  228.13330_gp,  167.91870_gp,  124.26420_gp,   93.86630_gp,&
         2051.11510_gp, 1403.09240_gp, 1240.30710_gp,  923.24960_gp,  687.29140_gp,&
         554.12740_gp,  440.44180_gp,  352.07430_gp, 3457.04340_gp, 2589.28880_gp,&
         2094.53890_gp, 1986.92650_gp, 1799.65770_gp, 1417.34210_gp, 1525.10220_gp,&
         1194.17780_gp, 1234.49140_gp, 1285.93000_gp,  986.84670_gp,  976.68650_gp,&
         1169.94940_gp,  986.64000_gp,  809.37540_gp,  709.43800_gp,  606.30260_gp,&
         516.14250_gp, 3848.83100_gp, 3098.78680_gp, 2617.93610_gp, 2300.99430_gp,&
         2068.41110_gp, 1558.74810_gp, 1754.97060_gp, 1301.55980_gp, 1420.72110_gp,&
         1305.33200_gp, 1096.81750_gp, 1136.70420_gp, 1469.07150_gp, 1308.51460_gp,&
         1125.84980_gp, 1023.19420_gp,  905.11410_gp,  797.33130_gp, 4691.69640_gp,&
         4003.70890_gp, 3391.56060_gp, 1339.62480_gp, 3505.93360_gp, 3344.21620_gp,&
         3254.42680_gp, 3172.47830_gp, 3099.62070_gp, 2345.03200_gp, 2808.41010_gp,&
         2691.79850_gp, 2769.09430_gp, 2707.31880_gp, 2650.38670_gp, 2623.45790_gp,&
         2153.85160_gp, 2039.39590_gp, 1820.53310_gp, 1513.38810_gp, 1525.64090_gp,&
         1356.64580_gp, 1223.66760_gp, 1001.47540_gp,  929.99430_gp,  948.96980_gp,&
         1497.60180_gp, 1420.69960_gp, 1264.24230_gp, 1184.16270_gp, 1069.86510_gp,&
         962.17680_gp, 4303.96450_gp, 3886.57160_gp, 3320.07640_gp, 2872.51360_gp,&
         2908.92060_gp,   78.20340_gp,   48.49800_gp, 1671.61560_gp,  829.71320_gp,&
         516.01950_gp,  329.28020_gp,  221.06590_gp,  162.79260_gp,  120.52900_gp,&
         91.08690_gp, 1981.46750_gp, 1357.16280_gp, 1200.02750_gp,  893.64440_gp/)
    vdwparams%coeffs( 4201: 4300)=(/  665.49650_gp,  536.70890_gp,  426.73540_gp,  341.23400_gp, 3337.69610_gp,&
         2503.35760_gp, 2025.47400_gp, 1921.68190_gp, 1740.72940_gp, 1370.92640_gp,&
         1475.36690_gp, 1155.27370_gp, 1194.53890_gp, 1244.20370_gp,  954.83090_gp,&
         945.27750_gp, 1132.11070_gp,  955.06070_gp,  783.70690_gp,  687.08390_gp,&
         587.34500_gp,  500.13680_gp, 3716.06780_gp, 2995.71860_gp, 2531.73730_gp,&
         2225.68070_gp, 2000.94390_gp, 1508.24940_gp, 1697.94840_gp, 1259.61420_gp,&
         1374.96650_gp, 1263.41360_gp, 1061.54570_gp, 1100.33990_gp, 1421.60990_gp,&
         1266.63290_gp, 1090.10530_gp,  990.86960_gp,  876.69020_gp,  772.45190_gp,&
         4529.17270_gp, 3869.76110_gp, 3279.33600_gp, 1297.04610_gp, 3389.03950_gp,&
         3233.09750_gp, 3146.36900_gp, 3067.20530_gp, 2996.82370_gp, 2268.01380_gp,&
         2714.02960_gp, 2601.32260_gp, 2677.56720_gp, 2617.87180_gp, 2562.86450_gp,&
         2536.79150_gp, 2083.03850_gp, 1973.08210_gp, 1761.68440_gp, 1464.65290_gp,&
         1476.63370_gp, 1313.30660_gp, 1184.75300_gp,  969.80220_gp,  900.65620_gp,&
         919.08240_gp, 1449.30610_gp, 1375.27520_gp, 1224.13220_gp, 1146.75090_gp,&
         1036.24960_gp,  932.11210_gp, 4156.51210_gp, 3757.24300_gp, 3210.75600_gp,&
         2778.94320_gp, 2813.62690_gp, 2721.52090_gp,   79.43910_gp,   48.72410_gp,&
         1800.50660_gp,  869.56410_gp,  532.89880_gp,  336.46270_gp,  224.20940_gp,&
         164.29590_gp,  121.16670_gp,   91.30900_gp, 2131.30070_gp, 1429.43080_gp,&
         1254.36040_gp,  924.19700_gp,  682.27500_gp,  547.31890_gp,  433.01050_gp,&
         344.83050_gp, 3608.56020_gp, 2661.11780_gp, 2144.92320_gp, 2028.78460_gp/)
    vdwparams%coeffs( 4301: 4400)=(/ 1834.18000_gp, 1445.06820_gp, 1550.18990_gp, 1213.87080_gp, 1248.79090_gp,&
         1303.33340_gp, 1000.95490_gp,  983.72150_gp, 1180.44990_gp,  987.37510_gp,&
         804.25620_gp,  702.03110_gp,  597.50350_gp,  506.82160_gp, 4013.89730_gp,&
         3187.57340_gp, 2676.64720_gp, 2343.75860_gp, 2101.84170_gp, 1577.87360_gp,&
         1779.03230_gp, 1313.71580_gp, 1433.15630_gp, 1314.65610_gp, 1106.86170_gp,&
         1142.51960_gp, 1484.01500_gp, 1312.66530_gp, 1122.14400_gp, 1016.00250_gp,&
         895.23210_gp,  785.80380_gp, 4893.04360_gp, 4127.93300_gp, 3475.71240_gp,&
         1337.28500_gp, 3607.88250_gp, 3437.85440_gp, 3344.52790_gp, 3259.44310_gp,&
         3183.76270_gp, 2393.55840_gp, 2896.40800_gp, 2773.03510_gp, 2839.64690_gp,&
         2775.74240_gp, 2716.65220_gp, 2689.70330_gp, 2198.66320_gp, 2068.00290_gp,&
         1839.04640_gp, 1525.11150_gp, 1534.79720_gp, 1360.23780_gp, 1223.64480_gp,&
         999.20960_gp,  927.02850_gp,  944.28620_gp, 1511.67140_gp, 1425.82660_gp,&
         1261.24150_gp, 1177.40380_gp, 1059.63770_gp,  949.75540_gp, 4469.77740_gp,&
         3995.26060_gp, 3393.26010_gp, 2917.09920_gp, 2964.66680_gp, 2867.10250_gp,&
         3032.97600_gp,   77.10210_gp,   47.37980_gp, 1727.81400_gp,  839.44780_gp,&
         515.79770_gp,  326.25000_gp,  217.67650_gp,  159.64350_gp,  117.81930_gp,&
         88.83650_gp, 2045.92840_gp, 1378.58760_gp, 1211.44430_gp,  894.32550_gp,&
         661.20920_gp,  530.88280_gp,  420.34980_gp,  334.97590_gp, 3458.68830_gp,&
         2561.27500_gp, 2066.17680_gp, 1955.45050_gp, 1768.56490_gp, 1393.13850_gp,&
         1495.56570_gp, 1171.00130_gp, 1206.01580_gp, 1258.21700_gp,  966.05990_gp/)
    vdwparams%coeffs( 4401: 4500)=(/  950.86030_gp, 1140.51970_gp,  955.50390_gp,  779.29760_gp,  680.73550_gp,&
         579.79560_gp,  492.11280_gp, 3847.73240_gp, 3067.14630_gp, 2579.00220_gp,&
         2260.06390_gp, 2027.74400_gp, 1523.37420_gp, 1717.10410_gp, 1269.06470_gp,&
         1384.73150_gp, 1270.65780_gp, 1069.27910_gp, 1104.74670_gp, 1433.43700_gp,&
         1269.73680_gp, 1086.76710_gp,  984.62910_gp,  868.18680_gp,  762.54070_gp,&
         4689.49090_gp, 3969.44880_gp, 3346.92780_gp, 1294.69970_gp, 3470.81030_gp,&
         3308.27820_gp, 3218.72530_gp, 3137.05760_gp, 3064.42240_gp, 2306.83660_gp,&
         2784.25860_gp, 2666.08160_gp, 2734.29740_gp, 2672.90420_gp, 2616.17240_gp,&
         2590.09570_gp, 2118.98110_gp, 1996.02660_gp, 1776.41870_gp, 1473.75440_gp,&
         1483.67470_gp, 1315.80800_gp, 1184.29550_gp,  967.46430_gp,  897.73860_gp,&
         914.79120_gp, 1460.22470_gp, 1379.02630_gp, 1221.26470_gp, 1140.78270_gp,&
         1027.39480_gp,  921.40470_gp, 4288.58580_gp, 3844.43250_gp, 3269.52880_gp,&
         2814.60330_gp, 2858.18750_gp, 2764.25850_gp, 2921.82230_gp, 2815.23660_gp,&
         4.73790_gp,    3.12870_gp,   68.93910_gp,   41.30780_gp,   28.27670_gp,&
         19.26530_gp,   13.51640_gp,   10.23710_gp,    7.74410_gp,    5.94030_gp,&
         82.56410_gp,   65.37030_gp,   60.85580_gp,   48.51820_gp,   38.13970_gp,&
         31.77130_gp,   26.01700_gp,   21.30550_gp,  134.44410_gp,  113.48570_gp,&
         94.20830_gp,   91.33340_gp,   83.83130_gp,   65.98370_gp,   72.41110_gp,&
         56.79930_gp,   60.57460_gp,   62.25320_gp,   47.64990_gp,   49.34690_gp,&
         58.42070_gp,   51.99050_gp,   44.65120_gp,   40.20100_gp,   35.27400_gp/)
    vdwparams%coeffs( 4501: 4600)=(/   30.72610_gp,  150.93400_gp,  135.06640_gp,  119.29610_gp,  107.73060_gp,&
         98.51830_gp,   76.36610_gp,   85.08000_gp,   65.08630_gp,   71.19850_gp,&
         66.12110_gp,   54.96580_gp,   58.35680_gp,   72.90760_gp,   67.94300_gp,&
         60.94700_gp,   56.74270_gp,   51.46560_gp,   46.37940_gp,  184.02230_gp,&
         171.70330_gp,  152.11060_gp,   71.88300_gp,  152.45020_gp,  146.58510_gp,&
         142.96200_gp,  139.62540_gp,  136.66950_gp,  108.36620_gp,  120.40560_gp,&
         116.31780_gp,  123.52200_gp,  120.93020_gp,  118.60700_gp,  117.18110_gp,&
         99.29840_gp,   98.51360_gp,   90.28580_gp,   76.41500_gp,   77.87930_gp,&
         70.80240_gp,   64.98900_gp,   54.04510_gp,   50.50860_gp,   52.05470_gp,&
         74.88170_gp,   73.66920_gp,   68.10310_gp,   65.16050_gp,   60.33380_gp,&
         55.43450_gp,  174.92120_gp,  170.50820_gp,  151.81390_gp,  137.47240_gp,&
         135.87010_gp,  131.56440_gp,  134.73040_gp,  130.59460_gp,    7.59160_gp,&
         14.31650_gp,    8.77730_gp,  282.21060_gp,  148.13230_gp,   93.67290_gp,&
         60.08600_gp,   40.29440_gp,   29.55750_gp,   21.76050_gp,   16.33880_gp,&
         335.34500_gp,  240.16270_gp,  214.72650_gp,  162.05080_gp,  121.50670_gp,&
         98.12490_gp,   77.98650_gp,   62.22050_gp,  555.15160_gp,  434.15760_gp,&
         354.03170_gp,  337.38960_gp,  306.54820_gp,  240.50330_gp,  260.90190_gp,&
         203.61490_gp,  212.89780_gp,  221.17390_gp,  168.80750_gp,  169.49490_gp,&
         202.69000_gp,  172.94250_gp,  142.84710_gp,  125.47680_gp,  107.33800_gp,&
         91.33330_gp,  618.56810_gp,  517.83930_gp,  442.94480_gp,  391.91790_gp/)
    vdwparams%coeffs( 4601: 4700)=(/  353.50080_gp,  267.37610_gp,  300.65330_gp,  223.88710_gp,  245.25610_gp,&
         225.80420_gp,  188.33740_gp,  197.09400_gp,  253.13700_gp,  228.08140_gp,&
         197.76980_gp,  180.30560_gp,  159.86700_gp,  140.98180_gp,  752.34180_gp,&
         664.57190_gp,  570.42660_gp,  234.52710_gp,  583.97970_gp,  558.69720_gp,&
         544.12320_gp,  530.79090_gp,  518.94980_gp,  396.72230_gp,  463.54420_gp,&
         445.01010_gp,  465.41000_gp,  455.28270_gp,  446.00300_gp,  441.32200_gp,&
         364.74670_gp,  350.06170_gp,  314.24830_gp,  261.38330_gp,  264.36970_gp,&
         236.02580_gp,  213.46440_gp,  174.61320_gp,  162.13420_gp,  166.04890_gp,&
         257.38230_gp,  246.83550_gp,  221.43250_gp,  208.12690_gp,  188.60730_gp,&
         169.89430_gp,  697.64490_gp,  649.03250_gp,  561.35190_gp,  491.07770_gp,&
         493.46510_gp,  477.42540_gp,  499.52020_gp,  482.34660_gp,   24.05730_gp,&
         85.31970_gp,   18.46560_gp,   11.46550_gp,  352.21600_gp,  186.34530_gp,&
         119.08150_gp,   77.09510_gp,   52.08520_gp,   38.41320_gp,   28.40980_gp,&
         21.40670_gp,  418.64880_gp,  301.37770_gp,  270.74610_gp,  205.75440_gp,&
         155.36440_gp,  126.09930_gp,  100.72120_gp,   80.71740_gp,  694.85960_gp,&
         543.58420_gp,  443.71900_gp,  423.62170_gp,  385.27540_gp,  302.60610_gp,&
         328.40120_gp,  256.60430_gp,  268.62840_gp,  278.71160_gp,  212.97970_gp,&
         214.42560_gp,  256.20490_gp,  219.73750_gp,  182.53380_gp,  160.96870_gp,&
         138.27470_gp,  118.12400_gp,  774.96560_gp,  648.64220_gp,  556.15030_gp,&
         493.01000_gp,  445.37360_gp,  337.90350_gp,  379.52540_gp,  283.56950_gp/)
    vdwparams%coeffs( 4701: 4800)=(/  310.37810_gp,  286.03950_gp,  238.69620_gp,  249.99150_gp,  320.09120_gp,&
         289.47130_gp,  252.16490_gp,  230.63090_gp,  205.22370_gp,  181.62600_gp,&
         943.85470_gp,  832.62510_gp,  716.07580_gp,  298.87820_gp,  732.16420_gp,&
         700.44420_gp,  682.19010_gp,  665.48680_gp,  650.65570_gp,  499.18370_gp,&
         582.12240_gp,  559.29170_gp,  583.70670_gp,  571.00390_gp,  559.38980_gp,&
         553.42660_gp,  458.65820_gp,  441.37890_gp,  397.13120_gp,  331.21240_gp,&
         335.22610_gp,  299.95650_gp,  271.81190_gp,  222.91570_gp,  207.18310_gp,&
         212.31320_gp,  326.10680_gp,  313.46750_gp,  282.29350_gp,  266.02000_gp,&
         241.86880_gp,  218.56870_gp,  875.85630_gp,  814.13080_gp,  705.40110_gp,&
         619.06050_gp,  621.33990_gp,  601.14960_gp,  627.63850_gp,  606.21390_gp,&
         30.78660_gp,  107.50830_gp,  135.84450_gp,   14.72370_gp,    9.58360_gp,&
         231.42090_gp,  133.72650_gp,   89.74290_gp,   60.29960_gp,   41.90150_gp,&
         31.53540_gp,   23.73690_gp,   18.14240_gp,  276.56710_gp,  212.99790_gp,&
         196.20140_gp,  154.29010_gp,  119.93530_gp,   99.21590_gp,   80.72720_gp,&
         65.75800_gp,  452.51930_gp,  373.88150_gp,  308.87210_gp,  298.08710_gp,&
         272.87470_gp,  214.63600_gp,  234.80000_gp,  183.97370_gp,  195.17310_gp,&
         201.13510_gp,  153.88240_gp,  158.08480_gp,  187.57970_gp,  165.19300_gp,&
         140.56250_gp,  125.84370_gp,  109.80600_gp,   95.17630_gp,  506.94210_gp,&
         445.26080_gp,  389.88230_gp,  350.19700_gp,  319.11510_gp,  245.85270_gp,&
         274.53340_gp,  208.63620_gp,  228.28220_gp,  211.55590_gp,  176.06030_gp/)
    vdwparams%coeffs( 4801: 4900)=(/  186.21990_gp,  234.21710_gp,  216.42760_gp,  192.55420_gp,  178.39190_gp,&
         160.96790_gp,  144.36970_gp,  617.60580_gp,  567.51360_gp,  498.45340_gp,&
         227.41890_gp,  502.49670_gp,  482.50770_gp,  470.40290_gp,  459.27530_gp,&
         449.41020_gp,  352.90870_gp,  397.44990_gp,  383.29520_gp,  405.32850_gp,&
         396.73770_gp,  388.98900_gp,  384.47010_gp,  323.64500_gp,  318.28720_gp,&
         290.18010_gp,  244.56750_gp,  248.77750_gp,  225.18190_gp,  205.97020_gp,&
         170.65560_gp,  159.26350_gp,  163.88200_gp,  240.00000_gp,  234.60480_gp,&
         215.30270_gp,  205.13340_gp,  188.99910_gp,  172.88540_gp,  582.98150_gp,&
         561.00860_gp,  495.63890_gp,  444.81020_gp,  441.55780_gp,  427.48070_gp,&
         440.27100_gp,  426.32550_gp,   23.84120_gp,   77.66330_gp,   98.99380_gp,&
         75.36860_gp,   11.39320_gp,    7.70650_gp,  160.45900_gp,   96.54310_gp,&
         66.79860_gp,   46.06810_gp,   32.70090_gp,   25.00970_gp,   19.09740_gp,&
         14.77160_gp,  192.37290_gp,  152.61300_gp,  142.60940_gp,  114.44050_gp,&
         90.67670_gp,   76.05520_gp,   62.74650_gp,   51.76840_gp,  314.58000_gp,&
         265.12060_gp,  220.18540_gp,  213.86820_gp,  196.48850_gp,  155.11710_gp,&
         169.97960_gp,  133.77000_gp,  142.51400_gp,  146.27110_gp,  112.39140_gp,&
         116.39870_gp,  137.47230_gp,  122.88050_gp,  106.16860_gp,   96.05990_gp,&
         84.77230_gp,   74.28650_gp,  353.63020_gp,  315.89230_gp,  279.47270_gp,&
         252.82130_gp,  231.59760_gp,  180.33900_gp,  200.54640_gp,  154.19380_gp,&
         168.31150_gp,  156.50060_gp,  130.51440_gp,  138.36720_gp,  172.12290_gp/)
    vdwparams%coeffs( 4901: 5000)=(/  160.79100_gp,  144.84620_gp,  135.31220_gp,  123.25380_gp,  111.59230_gp,&
         431.51720_gp,  401.87740_gp,  356.46100_gp,  170.92620_gp,  357.13990_gp,&
         343.41530_gp,  334.92450_gp,  327.09910_gp,  320.16570_gp,  254.81220_gp,&
         282.85470_gp,  273.33750_gp,  289.39200_gp,  283.30010_gp,  277.84930_gp,&
         274.45000_gp,  233.13940_gp,  231.63340_gp,  212.80410_gp,  180.84620_gp,&
         184.35860_gp,  168.04770_gp,  154.62480_gp,  129.18350_gp,  120.94850_gp,&
         124.59510_gp,  177.48530_gp,  174.76750_gp,  162.07400_gp,  155.45060_gp,&
         144.45320_gp,  133.24320_gp,  410.58600_gp,  399.52250_gp,  356.11270_gp,&
         323.49590_gp,  319.67610_gp,  309.60820_gp,  316.59670_gp,  306.93030_gp,&
         18.05750_gp,   56.34680_gp,   72.31620_gp,   56.52950_gp,   43.24520_gp,&
         8.14170_gp,    5.76010_gp,  101.27010_gp,   63.85800_gp,   45.73190_gp,&
         32.48480_gp,   23.62950_gp,   18.41280_gp,   14.30150_gp,   11.22430_gp,&
         121.96010_gp,  100.14430_gp,   95.09160_gp,   78.04660_gp,   63.16690_gp,&
         53.81070_gp,   45.09520_gp,   37.75220_gp,  199.44050_gp,  172.04550_gp,&
         143.76580_gp,  140.71150_gp,  129.82970_gp,  103.01890_gp,  113.02480_gp,&
         89.48680_gp,   95.70510_gp,   97.77610_gp,   75.59870_gp,   78.91830_gp,&
         92.60250_gp,   84.12760_gp,   73.91420_gp,   67.66360_gp,   60.47280_gp,&
         53.64960_gp,  225.19190_gp,  205.15740_gp,  183.68510_gp,  167.55010_gp,&
         154.43560_gp,  121.80660_gp,  134.78870_gp,  105.09050_gp,  114.34960_gp,&
         106.74380_gp,   89.34190_gp,   94.87700_gp,  116.49590_gp,  110.09510_gp/)
    vdwparams%coeffs( 5001: 5100)=(/  100.51820_gp,   94.75080_gp,   87.20280_gp,   79.77730_gp,  275.25670_gp,&
         260.44810_gp,  233.66190_gp,  118.54700_gp,  232.57830_gp,  224.00400_gp,&
         218.56110_gp,  213.52770_gp,  209.07000_gp,  169.05560_gp,  184.59040_gp,&
         178.80740_gp,  189.45310_gp,  185.49310_gp,  181.98230_gp,  179.62020_gp,&
         154.22920_gp,  154.91020_gp,  143.57630_gp,  123.23020_gp,  125.91210_gp,&
         115.69380_gp,  107.17750_gp,   90.44560_gp,   85.01460_gp,   87.64900_gp,&
         121.08000_gp,  120.12680_gp,  112.64850_gp,  108.81340_gp,  102.05240_gp,&
         94.99280_gp,  264.44300_gp,  260.65290_gp,  234.76520_gp,  216.30740_gp,&
         212.74560_gp,  206.16310_gp,  209.15240_gp,  203.04180_gp,   12.59800_gp,&
         37.49390_gp,   48.48850_gp,   39.04330_gp,   30.54100_gp,   22.12410_gp,&
         6.05750_gp,    4.45930_gp,   67.93120_gp,   44.49680_gp,   32.81360_gp,&
         23.91200_gp,   17.76980_gp,   14.07640_gp,   11.09940_gp,    8.82520_gp,&
         82.18200_gp,   69.36750_gp,   66.74510_gp,   55.81750_gp,   45.99370_gp,&
         39.70910_gp,   33.73130_gp,   28.60000_gp,  134.65080_gp,  118.24110_gp,&
         99.29390_gp,   97.83450_gp,   90.60070_gp,   72.28720_gp,   79.30470_gp,&
         63.18650_gp,   67.71170_gp,   68.90860_gp,   53.64200_gp,   56.28840_gp,&
         65.62140_gp,   60.40410_gp,   53.81530_gp,   49.75520_gp,   44.95100_gp,&
         40.30430_gp,  152.67570_gp,  141.15920_gp,  127.62820_gp,  117.23280_gp,&
         108.63780_gp,   86.68030_gp,   95.49660_gp,   75.39100_gp,   81.76400_gp,&
         76.58710_gp,   64.38800_gp,   68.39430_gp,   83.02780_gp,   79.17040_gp/)
    vdwparams%coeffs( 5101: 5200)=(/   73.07180_gp,   69.39010_gp,   64.41410_gp,   59.44680_gp,  186.87910_gp,&
         178.95490_gp,  162.03710_gp,   86.16650_gp,  160.52330_gp,  154.80670_gp,&
         151.09700_gp,  147.65520_gp,  144.60750_gp,  118.53880_gp,  127.78280_gp,&
         124.01580_gp,  131.30040_gp,  128.56660_gp,  126.16140_gp,  124.44020_gp,&
         107.83510_gp,  109.20420_gp,  101.97340_gp,   88.31560_gp,   90.39180_gp,&
         83.62270_gp,   77.91940_gp,   66.35910_gp,   62.60080_gp,   64.55520_gp,&
         86.95520_gp,   86.73920_gp,   82.05590_gp,   79.71330_gp,   75.32680_gp,&
         70.64800_gp,  181.00550_gp,  180.09880_gp,  163.60420_gp,  152.56500_gp,&
         149.55510_gp,  145.01130_gp,  146.19280_gp,  142.07750_gp,    9.18120_gp,&
         26.27080_gp,   34.19830_gp,   28.23270_gp,   22.51780_gp,   16.67750_gp,&
         12.81610_gp,    4.26720_gp,    3.30770_gp,   42.06270_gp,   28.86150_gp,&
         22.09310_gp,   16.63730_gp,   12.71070_gp,   10.28670_gp,    8.27390_gp,&
         6.69500_gp,   51.24980_gp,   44.71410_gp,   43.73640_gp,   37.44210_gp,&
         31.55770_gp,   27.71370_gp,   23.95050_gp,   20.63940_gp,   84.36880_gp,&
         75.58500_gp,   63.85640_gp,   63.47480_gp,   59.06650_gp,   47.53180_gp,&
         52.07800_gp,   41.89420_gp,   44.94070_gp,   45.51550_gp,   35.81000_gp,&
         37.75590_gp,   43.56850_gp,   40.75000_gp,   36.93890_gp,   34.57920_gp,&
         31.66990_gp,   28.78000_gp,   96.22350_gp,   90.40810_gp,   82.73370_gp,&
         76.67120_gp,   71.54930_gp,   57.97840_gp,   63.50790_gp,   50.97700_gp,&
         55.02400_gp,   51.77400_gp,   43.84800_gp,   46.53400_gp,   55.65890_gp/)
    vdwparams%coeffs( 5201: 5300)=(/   53.62430_gp,   50.14730_gp,   48.05350_gp,   45.08480_gp,   42.06340_gp,&
         117.97180_gp,  114.46760_gp,  104.80950_gp,   59.13790_gp,  103.36160_gp,&
         99.83310_gp,   97.47890_gp,   95.28460_gp,   93.34100_gp,   77.87420_gp,&
         82.71540_gp,   80.46070_gp,   84.94780_gp,   83.18270_gp,   81.64370_gp,&
         80.45750_gp,   70.55220_gp,   72.08870_gp,   67.95900_gp,   59.57580_gp,&
         61.09490_gp,   57.01110_gp,   53.52140_gp,   46.14420_gp,   43.74680_gp,&
         45.12050_gp,   58.91290_gp,   59.10970_gp,   56.50070_gp,   55.26140_gp,&
         52.70330_gp,   49.89270_gp,  115.43210_gp,  116.00620_gp,  106.53030_gp,&
         100.90290_gp,   98.57680_gp,   95.66250_gp,   95.73220_gp,   93.16140_gp,&
         6.30000_gp,   17.17330_gp,   22.54130_gp,   19.20290_gp,   15.68570_gp,&
         11.94700_gp,    9.40480_gp,    7.13410_gp,   20.75670_gp,   12.72870_gp,&
         425.71610_gp,  218.65720_gp,  136.95360_gp,   87.39560_gp,   58.48530_gp,&
         42.89110_gp,   31.60600_gp,   23.77450_gp,  505.43630_gp,  355.91800_gp,&
         316.45880_gp,  237.13260_gp,  176.92050_gp,  142.54300_gp,  113.10250_gp,&
         90.17350_gp,  840.82990_gp,  648.31050_gp,  527.09170_gp,  501.23930_gp,&
         454.80300_gp,  357.15080_gp,  386.34260_gp,  301.74010_gp,  314.17180_gp,&
         326.82710_gp,  249.81970_gp,  249.39250_gp,  298.40070_gp,  253.14330_gp,&
         208.18020_gp,  182.47440_gp,  155.82020_gp,  132.43480_gp,  936.30060_gp,&
         773.98540_gp,  658.79630_gp,  581.24800_gp,  523.40580_gp,  395.00580_gp,&
         444.52690_gp,  330.22730_gp,  361.42870_gp,  332.43920_gp,  277.94280_gp/)
    vdwparams%coeffs( 5301: 5400)=(/  289.84740_gp,  373.37940_gp,  334.70230_gp,  288.98880_gp,  262.88460_gp,&
         232.60530_gp,  204.79860_gp, 1138.92560_gp,  995.36600_gp,  850.07680_gp,&
         343.13820_gp,  873.51340_gp,  834.92630_gp,  812.93500_gp,  792.83360_gp,&
         774.97190_gp,  589.68210_gp,  695.02110_gp,  666.65150_gp,  694.07420_gp,&
         678.84860_gp,  664.85800_gp,  657.99220_gp,  542.08720_gp,  517.42950_gp,&
         463.27240_gp,  384.90910_gp,  388.79710_gp,  346.39970_gp,  312.81320_gp,&
         255.72890_gp,  237.40440_gp,  242.79190_gp,  379.72090_gp,  362.56490_gp,&
         323.93500_gp,  303.83030_gp,  274.73150_gp,  247.06350_gp, 1052.13180_gp,&
         969.74210_gp,  834.86500_gp,  726.93870_gp,  732.65540_gp,  708.79360_gp,&
         743.91900_gp,  717.91000_gp,   34.99520_gp,  125.82450_gp,  158.32680_gp,&
         113.40160_gp,   82.01440_gp,   54.44540_gp,   38.12270_gp,   24.95600_gp,&
         186.10520_gp,   31.99130_gp,   19.46480_gp,  670.03300_gp,  340.26020_gp,&
         212.21800_gp,  134.92270_gp,   89.96460_gp,   65.76090_gp,   48.28650_gp,&
         36.19110_gp,  794.75970_gp,  554.59670_gp,  492.08290_gp,  367.58110_gp,&
         273.52680_gp,  219.94890_gp,  174.13990_gp,  138.51900_gp, 1327.31760_gp,&
         1013.76950_gp,  822.89570_gp,  781.69930_gp,  708.76840_gp,  556.61270_gp,&
         601.44990_gp,  469.64370_gp,  488.17150_gp,  508.15850_gp,  388.42090_gp,&
         386.88840_gp,  463.56910_gp,  392.25330_gp,  321.87320_gp,  281.71730_gp,&
         240.15650_gp,  203.74730_gp, 1477.65560_gp, 1210.95930_gp, 1028.13650_gp,&
         905.77500_gp,  814.92890_gp,  614.02650_gp,  691.44970_gp,  512.67590_gp/)
    vdwparams%coeffs( 5401: 5500)=(/  561.00100_gp,  515.64570_gp,  431.28070_gp,  449.15670_gp,  579.92340_gp,&
         518.65180_gp,  446.92300_gp,  406.07280_gp,  358.80520_gp,  315.45810_gp,&
         1799.01800_gp, 1559.50500_gp, 1328.25100_gp,  530.86500_gp, 1367.27140_gp,&
         1305.90000_gp, 1271.28830_gp, 1239.67490_gp, 1211.58250_gp,  919.63290_gp,&
         1089.82140_gp, 1045.19610_gp, 1084.23130_gp, 1060.33860_gp, 1038.36010_gp,&
         1027.73560_gp,  845.53860_gp,  805.00140_gp,  719.71520_gp,  597.43440_gp,&
         603.09070_gp,  536.62210_gp,  484.07160_gp,  395.24010_gp,  366.70170_gp,&
         374.85520_gp,  589.53950_gp,  561.72980_gp,  500.92230_gp,  469.35230_gp,&
         423.85000_gp,  380.67930_gp, 1657.50700_gp, 1517.34570_gp, 1302.83760_gp,&
         1131.38980_gp, 1141.86680_gp, 1104.50830_gp, 1160.81350_gp, 1119.84660_gp,&
         54.10210_gp,  195.70570_gp,  246.22270_gp,  175.56460_gp,  126.66570_gp,&
         83.77030_gp,   58.43030_gp,   37.98840_gp,  289.60820_gp,  451.27080_gp,&
         33.10480_gp,   20.51580_gp,  627.53090_gp,  334.73350_gp,  213.88480_gp,&
         138.32540_gp,   93.35420_gp,   68.79830_gp,   50.85770_gp,   38.31540_gp,&
         746.34190_gp,  541.02850_gp,  486.23650_gp,  369.62320_gp,  278.92720_gp,&
         226.21560_gp,  180.54300_gp,  144.57980_gp, 1233.22370_gp,  973.04930_gp,&
         795.22060_gp,  759.45320_gp,  690.91560_gp,  542.24590_gp,  589.14750_gp,&
         460.04960_gp,  482.29900_gp,  500.35250_gp,  381.99230_gp,  385.14490_gp,&
         459.95440_gp,  394.66500_gp,  327.72290_gp,  288.84850_gp,  247.96920_gp,&
         211.69490_gp, 1375.25580_gp, 1160.25390_gp,  996.37570_gp,  883.85410_gp/)
    vdwparams%coeffs( 5501: 5600)=(/  798.60940_gp,  605.91350_gp,  680.54270_gp,  508.52890_gp,  557.00430_gp,&
         513.41630_gp,  427.99620_gp,  448.80710_gp,  574.34200_gp,  519.84240_gp,&
         452.84940_gp,  414.05440_gp,  368.27920_gp,  325.76250_gp, 1673.19220_gp,&
         1487.29640_gp, 1281.52290_gp,  536.56260_gp, 1308.60410_gp, 1252.69280_gp,&
         1220.22460_gp, 1190.50030_gp, 1164.10890_gp,  893.95150_gp, 1038.10470_gp,&
         997.38370_gp, 1045.00380_gp, 1022.36740_gp, 1001.67840_gp,  990.98080_gp,&
         821.60050_gp,  791.84560_gp,  712.74220_gp,  594.16160_gp,  601.57800_gp,&
         538.38550_gp,  487.89840_gp,  399.97680_gp,  371.72100_gp,  381.05900_gp,&
         584.69720_gp,  562.68750_gp,  506.86630_gp,  477.58620_gp,  434.09690_gp,&
         392.11690_gp, 1555.88480_gp, 1455.42410_gp, 1263.37240_gp, 1110.06280_gp,&
         1113.07380_gp, 1077.00050_gp, 1123.80080_gp, 1085.68460_gp,   55.26360_gp,&
         193.07070_gp,  243.76120_gp,  177.72540_gp,  129.69220_gp,   86.87280_gp,&
         61.21750_gp,   40.32200_gp,  284.22510_gp,  441.73220_gp,  437.75750_gp,&
         27.49890_gp,   17.63560_gp,  450.95830_gp,  256.58810_gp,  170.18910_gp,&
         113.20460_gp,   78.01900_gp,   58.35570_gp,   43.68550_gp,   33.24160_gp,&
         538.34120_gp,  409.91740_gp,  375.52330_gp,  292.99940_gp,  226.06570_gp,&
         186.00510_gp,  150.52310_gp,  122.00120_gp,  881.26240_gp,  722.51110_gp,&
         595.67250_gp,  573.49880_gp,  524.28080_gp,  411.88790_gp,  450.22340_gp,&
         352.24080_gp,  373.01250_gp,  385.00480_gp,  294.14060_gp,  301.17090_gp,&
         357.94110_gp,  313.39110_gp,  265.07110_gp,  236.33230_gp,  205.29760_gp/)
    vdwparams%coeffs( 5601: 5700)=(/  177.18110_gp,  986.04380_gp,  860.34210_gp,  750.38170_gp,  672.18340_gp,&
         611.31470_gp,  469.10920_gp,  524.64790_gp,  396.98620_gp,  434.72510_gp,&
         402.36620_gp,  334.66730_gp,  353.59760_gp,  446.59370_gp,  410.89600_gp,&
         363.78460_gp,  335.92400_gp,  301.98520_gp,  269.83720_gp, 1200.58030_gp,&
         1097.39800_gp,  960.23320_gp,  429.80740_gp,  970.39570_gp,  931.29170_gp,&
         907.79700_gp,  886.21930_gp,  867.08440_gp,  677.48060_gp,  767.30550_gp,&
         739.39750_gp,  781.37480_gp,  764.76680_gp,  749.74130_gp,  741.20300_gp,&
         621.81690_gp,  609.06080_gp,  553.65420_gp,  465.19420_gp,  472.79230_gp,&
         426.79800_gp,  389.49880_gp,  321.71830_gp,  299.88580_gp,  308.42630_gp,&
         456.58720_gp,  445.00360_gp,  406.70170_gp,  386.46750_gp,  354.86660_gp,&
         323.54110_gp, 1129.89150_gp, 1082.50180_gp,  953.14260_gp,  851.47580_gp,&
         846.79590_gp,  819.68560_gp,  846.55430_gp,  819.36980_gp,   44.90300_gp,&
         148.74780_gp,  189.10900_gp,  142.51230_gp,  106.08060_gp,   72.63730_gp,&
         52.13110_gp,   35.12900_gp,  217.50490_gp,  337.02500_gp,  339.65160_gp,&
         270.26500_gp,   23.09080_gp,   15.26190_gp,  339.66930_gp,  202.22820_gp,&
         138.07740_gp,   93.96630_gp,   65.91270_gp,   49.93740_gp,   37.80350_gp,&
         29.02800_gp,  406.71500_gp,  320.41760_gp,  297.80760_gp,  236.98640_gp,&
         186.06840_gp,  154.92290_gp,  126.83120_gp,  103.86530_gp,  663.11840_gp,&
         557.42550_gp,  462.35110_gp,  447.95900_gp,  411.01210_gp,  323.59140_gp,&
         354.84320_gp,  278.40200_gp,  296.57890_gp,  304.91450_gp,  233.48980_gp/)
    vdwparams%coeffs( 5701: 5800)=(/  241.43800_gp,  285.80720_gp,  253.97440_gp,  217.89200_gp,  196.08060_gp,&
         171.99100_gp,  149.79330_gp,  744.24750_gp,  663.56560_gp,  585.24670_gp,&
         528.08410_gp,  482.70360_gp,  373.94270_gp,  416.70800_gp,  318.59120_gp,&
         348.44490_gp,  323.52760_gp,  269.12200_gp,  285.48020_gp,  356.90690_gp,&
         332.17220_gp,  297.65980_gp,  276.98090_gp,  251.10810_gp,  226.22110_gp,&
         907.29930_gp,  844.01200_gp,  746.59930_gp,  351.15850_gp,  749.14380_gp,&
         720.14500_gp,  702.29590_gp,  685.86200_gp,  671.30090_gp,  531.54240_gp,&
         592.01350_gp,  571.74810_gp,  606.49400_gp,  593.74030_gp,  582.29560_gp,&
         575.32620_gp,  487.06190_gp,  482.46160_gp,  441.86160_gp,  373.86540_gp,&
         380.91200_gp,  346.13260_gp,  317.60850_gp,  264.10910_gp,  246.83080_gp,&
         254.31280_gp,  366.58420_gp,  360.26310_gp,  332.71860_gp,  318.18410_gp,&
         294.47100_gp,  270.46770_gp,  861.38240_gp,  837.51280_gp,  744.76320_gp,&
         673.56880_gp,  666.25820_gp,  645.14050_gp,  661.28050_gp,  640.87610_gp,&
         37.02360_gp,  117.74720_gp,  150.60680_gp,  116.38180_gp,   88.07690_gp,&
         61.42860_gp,   44.77700_gp,   30.76880_gp,  171.44690_gp,  265.03780_gp,&
         270.35800_gp,  219.28630_gp,  180.62380_gp,   19.86670_gp,   13.45980_gp,&
         271.30110_gp,  166.18070_gp,  115.81680_gp,   80.18040_gp,   57.02550_gp,&
         43.65480_gp,   33.35500_gp,   25.81070_gp,  325.59690_gp,  261.94370_gp,&
         245.85680_gp,  198.32660_gp,  157.69710_gp,  132.49220_gp,  109.44860_gp,&
         90.37160_gp,  530.51920_gp,  452.37450_gp,  376.61420_gp,  366.51510_gp/)
    vdwparams%coeffs( 5801: 5900)=(/  337.12500_gp,  266.07620_gp,  292.12570_gp,  229.87760_gp,  245.59880_gp,&
         251.80970_gp,  193.38650_gp,  201.07040_gp,  237.25630_gp,  212.94370_gp,&
         184.54670_gp,  167.22640_gp,  147.76790_gp,  129.61230_gp,  596.87630_gp,&
         538.67240_gp,  478.52040_gp,  433.90360_gp,  398.04190_gp,  310.58610_gp,&
         345.14340_gp,  265.95280_gp,  290.41460_gp,  270.25240_gp,  225.11810_gp,&
         239.17840_gp,  296.82140_gp,  278.27870_gp,  251.43440_gp,  235.25190_gp,&
         214.60020_gp,  194.52370_gp,  728.46530_gp,  684.23530_gp,  609.44200_gp,&
         296.45810_gp,  608.91190_gp,  585.91340_gp,  571.53970_gp,  558.28160_gp,&
         546.53950_gp,  436.77850_gp,  481.58940_gp,  465.77580_gp,  494.52850_gp,&
         484.17910_gp,  474.94430_gp,  469.05670_gp,  399.59520_gp,  398.59150_gp,&
         366.94770_gp,  312.20810_gp,  318.55560_gp,  290.82630_gp,  267.91140_gp,&
         224.01000_gp,  209.79980_gp,  216.32530_gp,  306.15090_gp,  302.35050_gp,&
         281.17730_gp,  270.07790_gp,  251.35500_gp,  232.12580_gp,  695.53530_gp,&
         681.65580_gp,  609.93410_gp,  556.24920_gp,  548.49710_gp,  531.25540_gp,&
         541.87030_gp,  525.57770_gp,   31.41270_gp,   97.08490_gp,  124.75120_gp,&
         98.11390_gp,   75.20390_gp,   53.20860_gp,   39.26830_gp,   27.41170_gp,&
         141.05990_gp,  217.66290_gp,  223.79420_gp,  183.94830_gp,  153.18230_gp,&
         130.99650_gp,   16.52730_gp,   11.50920_gp,  208.77430_gp,  131.64980_gp,&
         93.77430_gp,   66.13540_gp,   47.75400_gp,   36.97830_gp,   28.54680_gp,&
         22.28320_gp,  251.23260_gp,  206.45210_gp,  195.76470_gp,  160.18790_gp/)
    vdwparams%coeffs( 5901: 6000)=(/  129.09560_gp,  109.52700_gp,   91.36670_gp,   76.12560_gp,  409.51050_gp,&
         354.09350_gp,  295.90280_gp,  289.35120_gp,  266.85570_gp,  211.29280_gp,&
         232.14470_gp,  183.36610_gp,  196.37550_gp,  200.75710_gp,  154.77560_gp,&
         161.73370_gp,  190.11390_gp,  172.39330_gp,  150.99990_gp,  137.84430_gp,&
         122.77630_gp,  108.52260_gp,  462.02520_gp,  421.88670_gp,  377.56420_gp,&
         344.13180_gp,  316.91100_gp,  249.26710_gp,  276.14550_gp,  214.63800_gp,&
         233.90310_gp,  218.19030_gp,  182.15470_gp,  193.72730_gp,  238.49880_gp,&
         225.23810_gp,  205.25120_gp,  193.14490_gp,  177.34580_gp,  161.81190_gp,&
         564.60060_gp,  535.26780_gp,  480.12440_gp,  241.92230_gp,  477.71970_gp,&
         460.12500_gp,  448.95170_gp,  438.62430_gp,  429.48100_gp,  346.64610_gp,&
         378.39840_gp,  366.51660_gp,  389.19850_gp,  381.08470_gp,  373.88780_gp,&
         369.07690_gp,  316.53270_gp,  317.92390_gp,  294.29820_gp,  251.96160_gp,&
         257.44470_gp,  236.20650_gp,  218.51340_gp,  183.84350_gp,  172.59550_gp,&
         178.06050_gp,  247.21900_gp,  245.29500_gp,  229.72670_gp,  221.65670_gp,&
         207.50080_gp,  192.73160_gp,  542.29530_gp,  535.47060_gp,  482.17720_gp,&
         443.62620_gp,  436.15260_gp,  422.58040_gp,  428.86190_gp,  416.31150_gp,&
         25.74140_gp,   77.20230_gp,   99.69950_gp,   79.88330_gp,   62.08640_gp,&
         44.62170_gp,   33.38010_gp,   23.70590_gp,  111.98290_gp,  172.41480_gp,&
         178.71270_gp,  148.95590_gp,  125.49300_gp,  108.28930_gp,   90.39850_gp,&
         29.30000_gp,   18.66550_gp,  560.29540_gp,  292.22240_gp,  186.95660_gp/)
    vdwparams%coeffs( 6001: 6100)=(/  121.92920_gp,   83.21530_gp,   61.99480_gp,   46.35200_gp,   35.29690_gp,&
         666.03570_gp,  473.76610_gp,  424.89260_gp,  322.82330_gp,  244.57290_gp,&
         199.46720_gp,  160.33040_gp,  129.42830_gp, 1112.22770_gp,  860.17710_gp,&
         700.72810_gp,  668.75460_gp,  607.99510_gp,  478.97440_gp,  518.05830_gp,&
         406.08760_gp,  423.38090_gp,  439.28450_gp,  337.11170_gp,  337.89860_gp,&
         403.07450_gp,  345.44340_gp,  287.55350_gp,  254.36700_gp,  219.46470_gp,&
         188.48190_gp, 1240.86310_gp, 1027.96710_gp,  879.03460_gp,  778.52380_gp,&
         703.30600_gp,  534.63030_gp,  599.93020_gp,  449.30300_gp,  490.58070_gp,&
         452.24640_gp,  379.08360_gp,  395.51350_gp,  505.68380_gp,  456.46130_gp,&
         397.84160_gp,  364.40900_gp,  325.10680_gp,  288.72750_gp, 1511.88980_gp,&
         1322.26190_gp, 1133.77490_gp,  472.24550_gp, 1162.28350_gp, 1111.26380_gp,&
         1082.08650_gp, 1055.38910_gp, 1031.67160_gp,  790.78540_gp,  927.36720_gp,&
         890.50560_gp,  924.61970_gp,  904.32570_gp,  885.75380_gp,  876.28990_gp,&
         725.64410_gp,  696.18310_gp,  626.34420_gp,  523.68920_gp,  529.62490_gp,&
         474.25390_gp,  430.19640_gp,  354.14110_gp,  329.64710_gp,  337.33350_gp,&
         516.95330_gp,  495.69810_gp,  446.30330_gp,  420.86840_gp,  383.35330_gp,&
         347.34370_gp, 1400.16860_gp, 1291.50730_gp, 1115.92880_gp,  978.19280_gp,&
         983.91310_gp,  952.06270_gp,  995.35270_gp,  961.05430_gp,   48.44990_gp,&
         168.74120_gp,  213.44440_gp,  155.77810_gp,  114.57230_gp,   77.63170_gp,&
         55.40320_gp,   37.19820_gp,  249.22000_gp,  387.33490_gp,  382.64220_gp/)
    vdwparams%coeffs( 6101: 6200)=(/  296.94090_gp,  237.21910_gp,  197.38060_gp,  158.73300_gp,  338.02070_gp,&
         52.98710_gp,   31.87960_gp, 1310.03140_gp,  602.56610_gp,  362.68620_gp,&
         225.96290_gp,  148.97930_gp,  108.27180_gp,   79.22180_gp,   59.27790_gp,&
         1545.70640_gp,  997.24860_gp,  866.95170_gp,  630.04560_gp,  460.26350_gp,&
         366.77110_gp,  288.21820_gp,  228.08920_gp, 2661.08300_gp, 1886.42210_gp,&
         1510.14270_gp, 1422.32440_gp, 1282.16100_gp, 1011.59550_gp, 1079.15320_gp,&
         845.46860_gp,  862.50360_gp,  902.60140_gp,  694.50580_gp,  675.02990_gp,&
         813.49430_gp,  672.69170_gp,  543.01230_gp,  471.48220_gp,  399.02260_gp,&
         336.64830_gp, 2958.04100_gp, 2265.41550_gp, 1882.13380_gp, 1637.96670_gp,&
         1463.88900_gp, 1092.67330_gp, 1234.93160_gp,  905.75150_gp,  986.37270_gp,&
         902.44050_gp,  762.61740_gp,  781.66000_gp, 1024.16140_gp,  896.41550_gp,&
         759.79880_gp,  684.70730_gp,  600.26980_gp,  524.35470_gp, 3619.26600_gp,&
         2951.33620_gp, 2456.76010_gp,  907.40120_gp, 2571.53190_gp, 2442.04950_gp,&
         2374.04450_gp, 2312.21470_gp, 2257.19930_gp, 1680.63380_gp, 2081.61730_gp,&
         1992.87280_gp, 2006.29580_gp, 1960.23910_gp, 1917.46400_gp, 1899.13060_gp,&
         1544.68520_gp, 1435.15030_gp, 1268.75180_gp, 1049.32960_gp, 1052.99220_gp,&
         928.45210_gp,  831.87520_gp,  676.99020_gp,  627.08530_gp,  637.11120_gp,&
         1043.05560_gp,  974.28450_gp,  854.66400_gp,  794.44910_gp,  711.41430_gp,&
         634.79210_gp, 3271.25630_gp, 2840.94470_gp, 2386.76370_gp, 2029.62450_gp,&
         2076.11720_gp, 2006.69150_gp, 2135.22690_gp, 2054.01310_gp,   90.79270_gp/)
    vdwparams%coeffs( 6201: 6300)=(/  345.61200_gp,  434.02840_gp,  298.72640_gp,  212.89570_gp,  139.02190_gp,&
         96.15840_gp,   62.05800_gp,  516.69170_gp,  809.21820_gp,  775.08320_gp,&
         576.29620_gp,  446.01700_gp,  363.19510_gp,  285.49050_gp,  688.99110_gp,&
         1525.18910_gp,   29.99530_gp,   19.38410_gp,  504.82680_gp,  282.44200_gp,&
         186.14190_gp,  123.63790_gp,   85.34040_gp,   64.00720_gp,   48.09360_gp,&
         36.74680_gp,  602.49290_gp,  452.76000_gp,  412.90580_gp,  320.58960_gp,&
         246.73580_gp,  202.96720_gp,  164.36550_gp,  133.44040_gp,  988.58810_gp,&
         803.23800_gp,  660.67800_gp,  635.04570_gp,  579.95840_gp,  456.20230_gp,&
         497.34740_gp,  389.59740_gp,  411.05720_gp,  424.66930_gp,  325.08500_gp,&
         331.23100_gp,  393.53000_gp,  343.16790_gp,  289.55820_gp,  258.00840_gp,&
         224.13750_gp,  193.58680_gp, 1105.46080_gp,  957.22710_gp,  831.66790_gp,&
         743.43180_gp,  675.33140_gp,  517.73400_gp,  579.16230_gp,  437.84320_gp,&
         478.95980_gp,  443.08160_gp,  369.49330_gp,  389.18590_gp,  492.19040_gp,&
         451.13400_gp,  398.31520_gp,  367.39080_gp,  330.04180_gp,  294.87060_gp,&
         1344.52740_gp, 1222.74260_gp, 1065.81400_gp,  471.18540_gp, 1079.96470_gp,&
         1036.20550_gp, 1009.88930_gp,  985.72640_gp,  964.28550_gp,  750.78480_gp,&
         855.14710_gp,  822.84200_gp,  868.11040_gp,  849.53940_gp,  832.69610_gp,&
         823.30890_gp,  688.62420_gp,  671.92580_gp,  609.75150_gp,  512.21050_gp,&
         520.05780_gp,  468.95650_gp,  427.68450_gp,  353.48450_gp,  329.59850_gp,&
         338.54780_gp,  503.63930_gp,  489.25410_gp,  445.92310_gp,  423.20630_gp/)
    vdwparams%coeffs( 6301: 6400)=(/  388.21220_gp,  353.80050_gp, 1263.01910_gp, 1203.98590_gp, 1056.35410_gp,&
         940.62410_gp,  937.65790_gp,  907.69560_gp,  939.86910_gp,  909.28950_gp,&
         48.97530_gp,  163.64460_gp,  207.81070_gp,  155.82280_gp,  115.97260_gp,&
         79.53830_gp,   57.26110_gp,   38.80970_gp,  239.97140_gp,  371.77680_gp,&
         373.18450_gp,  295.61010_gp,  239.37010_gp,  200.76220_gp,  162.68060_gp,&
         327.70920_gp,  640.70760_gp,  324.50680_gp,   29.00090_gp,   18.80840_gp,&
         488.29310_gp,  272.02010_gp,  179.50770_gp,  119.43760_gp,   82.58340_gp,&
         62.03090_gp,   46.67650_gp,   35.71170_gp,  582.52700_gp,  436.10910_gp,&
         397.86180_gp,  309.10810_gp,  238.16440_gp,  196.11280_gp,  158.98920_gp,&
         129.21800_gp,  960.33760_gp,  774.68630_gp,  636.87060_gp,  612.22250_gp,&
         559.11570_gp,  440.07620_gp,  479.49690_gp,  375.82420_gp,  396.27480_gp,&
         409.35550_gp,  313.60260_gp,  319.38340_gp,  379.41090_gp,  330.97100_gp,&
         279.49900_gp,  249.22580_gp,  216.69210_gp,  187.32310_gp, 1074.19730_gp,&
         923.65570_gp,  802.03240_gp,  716.86090_gp,  651.28470_gp,  499.50110_gp,&
         558.72110_gp,  422.57080_gp,  462.03140_gp,  427.45770_gp,  356.72450_gp,&
         375.52610_gp,  474.81120_gp,  435.19480_gp,  384.44870_gp,  354.77930_gp,&
         318.91710_gp,  285.13250_gp, 1308.64180_gp, 1180.94240_gp, 1028.41690_gp,&
         454.82000_gp, 1043.31870_gp, 1000.21550_gp,  974.72270_gp,  951.32570_gp,&
         930.56540_gp,  724.58860_gp,  827.79700_gp,  797.23640_gp,  837.46180_gp,&
         819.49730_gp,  803.20320_gp,  794.14040_gp,  664.77270_gp,  647.87460_gp/)
    vdwparams%coeffs( 6401: 6500)=(/  587.98130_gp,  494.23570_gp,  501.77330_gp,  452.57260_gp,  412.85580_gp,&
         341.45400_gp,  318.46190_gp,  327.07700_gp,  486.21160_gp,  472.15520_gp,&
         430.48050_gp,  408.69510_gp,  375.10290_gp,  342.05750_gp, 1226.39700_gp,&
         1162.36600_gp, 1019.14220_gp,  907.23670_gp,  904.83310_gp,  875.87480_gp,&
         906.96160_gp,  877.36450_gp,   47.27550_gp,  157.66670_gp,  200.36690_gp,&
         150.35200_gp,  112.05230_gp,   76.98340_gp,   55.51600_gp,   37.72220_gp,&
         231.26120_gp,  358.36020_gp,  359.64970_gp,  285.09570_gp,  231.07730_gp,&
         193.98450_gp,  157.36040_gp,  316.25180_gp,  619.43010_gp,  312.92020_gp,&
         301.99880_gp,   28.40670_gp,   18.43590_gp,  479.14020_gp,  266.70500_gp,&
         175.89090_gp,  117.00980_gp,   80.91440_gp,   60.79330_gp,   45.76300_gp,&
         35.02900_gp,  571.64550_gp,  427.71560_gp,  390.04320_gp,  302.89850_gp,&
         233.31760_gp,  192.11130_gp,  155.75070_gp,  126.60270_gp,  941.98370_gp,&
         760.07180_gp,  624.78010_gp,  600.52210_gp,  548.39080_gp,  431.66830_gp,&
         470.25290_gp,  368.61360_gp,  388.57290_gp,  401.43290_gp,  307.57630_gp,&
         313.12860_gp,  371.94270_gp,  324.34350_gp,  273.83520_gp,  244.15540_gp,&
         212.27920_gp,  183.51710_gp, 1053.58340_gp,  906.23640_gp,  786.72550_gp,&
         703.07840_gp,  638.70140_gp,  489.81600_gp,  547.89060_gp,  414.35750_gp,&
         453.02620_gp,  419.11340_gp,  349.82880_gp,  368.18710_gp,  465.56470_gp,&
         426.59330_gp,  376.75440_gp,  347.63740_gp,  312.47140_gp,  279.36120_gp,&
         1283.11930_gp, 1158.66620_gp, 1008.83360_gp,  445.75700_gp, 1023.53950_gp/)
    vdwparams%coeffs( 6501: 6600)=(/  981.34570_gp,  956.33210_gp,  933.37380_gp,  913.00140_gp,  710.72370_gp,&
         811.99050_gp,  781.82720_gp,  821.62740_gp,  803.99890_gp,  788.00640_gp,&
         779.12190_gp,  651.97940_gp,  635.29960_gp,  576.49460_gp,  484.55830_gp,&
         491.91200_gp,  443.64130_gp,  404.68640_gp,  334.71530_gp,  312.18730_gp,&
         320.60120_gp,  476.75960_gp,  462.87410_gp,  421.91810_gp,  400.51720_gp,&
         367.55880_gp,  335.16100_gp, 1202.70240_gp, 1140.34300_gp,  999.65520_gp,&
         889.72590_gp,  887.48800_gp,  859.09970_gp,  889.77190_gp,  860.71780_gp,&
         46.30890_gp,  154.57620_gp,  196.40270_gp,  147.31530_gp,  109.77970_gp,&
         75.43030_gp,   54.41090_gp,   36.99620_gp,  226.78860_gp,  351.39280_gp,&
         352.55750_gp,  279.35630_gp,  226.37670_gp,  190.02990_gp,  154.15760_gp,&
         310.11080_gp,  607.55710_gp,  306.73890_gp,  296.00290_gp,  290.14360_gp,&
         28.05620_gp,   18.17540_gp,  480.59440_gp,  265.40880_gp,  174.37190_gp,&
         115.71820_gp,   79.90240_gp,   59.98290_gp,   45.12900_gp,   34.53500_gp,&
         573.14360_gp,  426.22760_gp,  387.86950_gp,  300.39280_gp,  230.90900_gp,&
         189.90520_gp,  153.80570_gp,  124.92650_gp,  945.87000_gp,  759.36840_gp,&
         623.53320_gp,  598.79810_gp,  546.52440_gp,  430.23530_gp,  468.29610_gp,&
         367.07850_gp,  386.44860_gp,  399.45300_gp,  306.11810_gp,  311.06030_gp,&
         369.63020_gp,  321.64040_gp,  271.07700_gp,  241.45850_gp,  209.73930_gp,&
         181.18190_gp, 1057.57180_gp,  905.61870_gp,  784.74890_gp,  700.54320_gp,&
         635.96390_gp,  487.19110_gp,  545.17240_gp,  411.81870_gp,  450.19430_gp/)
    vdwparams%coeffs( 6601: 6700)=(/  416.32300_gp,  347.67860_gp,  365.55310_gp,  462.83040_gp,  423.32080_gp,&
         373.25910_gp,  344.09820_gp,  309.00630_gp,  276.04060_gp, 1287.96160_gp,&
         1158.66720_gp, 1006.96360_gp,  441.78090_gp, 1022.97730_gp,  980.48290_gp,&
         955.40420_gp,  932.39460_gp,  911.97350_gp,  708.57690_gp,  812.04950_gp,&
         781.62600_gp,  820.30710_gp,  802.66090_gp,  786.63360_gp,  777.82320_gp,&
         650.05720_gp,  632.20990_gp,  573.09450_gp,  481.38140_gp,  488.47190_gp,&
         440.16550_gp,  401.25220_gp,  331.69650_gp,  309.30750_gp,  317.51490_gp,&
         473.84960_gp,  459.37430_gp,  418.10440_gp,  396.57360_gp,  363.60680_gp,&
         331.30140_gp, 1205.44360_gp, 1139.27610_gp,  997.02490_gp,  885.76610_gp,&
         884.42800_gp,  856.10390_gp,  887.70660_gp,  858.53020_gp,   45.81730_gp,&
         153.75390_gp,  195.23170_gp,  145.94770_gp,  108.58200_gp,   74.48900_gp,&
         53.67330_gp,   36.45940_gp,  225.78370_gp,  349.93280_gp,  350.42960_gp,&
         276.95120_gp,  224.03030_gp,  187.85700_gp,  152.23830_gp,  308.39010_gp,&
         607.07000_gp,  304.30940_gp,  293.65980_gp,  287.86570_gp,  285.69050_gp,&
         26.18100_gp,   17.15260_gp,  416.08360_gp,  238.96060_gp,  159.84920_gp,&
         107.33820_gp,   74.69980_gp,   56.36220_gp,   42.57640_gp,   32.68060_gp,&
         497.35820_gp,  381.29680_gp,  350.45220_gp,  274.92570_gp,  213.41010_gp,&
         176.51650_gp,  143.70010_gp,  117.20570_gp,  814.24680_gp,  670.83370_gp,&
         553.79310_gp,  534.12660_gp,  488.79470_gp,  384.69310_gp,  420.42070_gp,&
         329.62730_gp,  349.19410_gp,  360.01720_gp,  275.72060_gp,  282.66670_gp/)
    vdwparams%coeffs( 6701: 6800)=(/  335.06150_gp,  294.51840_gp,  250.28310_gp,  223.98570_gp,  195.45360_gp,&
         169.50970_gp,  911.94330_gp,  799.01250_gp,  698.67600_gp,  627.07900_gp,&
         571.18620_gp,  439.93260_gp,  491.32100_gp,  373.34420_gp,  408.36770_gp,&
         378.42950_gp,  315.34820_gp,  333.15180_gp,  419.07070_gp,  386.61840_gp,&
         343.52050_gp,  318.04850_gp,  286.87200_gp,  257.28030_gp, 1110.35880_gp,&
         1018.66880_gp,  893.50560_gp,  405.87370_gp,  902.04820_gp,  866.04420_gp,&
         844.28560_gp,  824.28300_gp,  806.54320_gp,  632.57430_gp,  713.75520_gp,&
         688.11250_gp,  727.23780_gp,  711.80010_gp,  697.85650_gp,  689.78470_gp,&
         580.12870_gp,  569.40560_gp,  518.79050_gp,  437.20880_gp,  444.59970_gp,&
         402.30970_gp,  367.94540_gp,  305.05010_gp,  284.80180_gp,  292.91770_gp,&
         429.59510_gp,  419.41490_gp,  384.44960_gp,  366.05470_gp,  337.09280_gp,&
         308.29470_gp, 1047.20980_gp, 1006.28590_gp,  888.23000_gp,  796.31050_gp,&
         791.25770_gp,  766.09840_gp,  789.92350_gp,  764.80610_gp,   42.39270_gp,&
         138.74910_gp,  176.68170_gp,  134.21360_gp,  100.64630_gp,   69.62150_gp,&
         50.48580_gp,   34.54810_gp,  202.96760_gp,  314.01190_gp,  317.32870_gp,&
         253.88000_gp,  207.09950_gp,  174.57750_gp,  142.21510_gp,  278.72590_gp,&
         535.81490_gp,  278.19650_gp,  268.46570_gp,  263.11270_gp,  260.77730_gp,&
         239.55540_gp,   24.34520_gp,   16.05600_gp,  379.21170_gp,  219.47750_gp,&
         147.62980_gp,   99.58800_gp,   69.56380_gp,   52.63510_gp,   39.86280_gp,&
         30.66480_gp,  453.54440_gp,  349.71160_gp,  322.27810_gp,  253.76380_gp/)
    vdwparams%coeffs( 6801: 6900)=(/  197.65690_gp,  163.87940_gp,  133.73240_gp,  109.31550_gp,  742.25200_gp,&
         613.92330_gp,  507.34010_gp,  489.89960_gp,  448.62620_gp,  353.27090_gp,&
         386.26020_gp,  303.04800_gp,  321.34480_gp,  331.06140_gp,  253.70380_gp,&
         260.53410_gp,  308.55290_gp,  271.97110_gp,  231.76180_gp,  207.79480_gp,&
         181.68360_gp,  157.86610_gp,  831.79600_gp,  731.23300_gp,  640.66410_gp,&
         575.77520_gp,  524.95880_gp,  405.08170_gp,  452.08310_gp,  344.23220_gp,&
         376.40450_gp,  349.02580_gp,  290.90830_gp,  307.52090_gp,  386.06280_gp,&
         356.91590_gp,  317.86240_gp,  294.73480_gp,  266.29100_gp,  239.21900_gp,&
         1013.04590_gp,  931.86100_gp,  818.89940_gp,  375.45680_gp,  825.77350_gp,&
         793.02390_gp,  773.15750_gp,  754.88610_gp,  738.68360_gp,  580.77730_gp,&
         653.47520_gp,  630.25330_gp,  666.33780_gp,  652.21550_gp,  639.47820_gp,&
         632.01230_gp,  532.44100_gp,  523.61220_gp,  477.73970_gp,  403.18600_gp,&
         410.18490_gp,  371.64480_gp,  340.26720_gp,  282.50560_gp,  263.90350_gp,&
         271.50750_gp,  396.15750_gp,  387.34040_gp,  355.75280_gp,  339.14780_gp,&
         312.79910_gp,  286.50170_gp,  956.82050_gp,  921.49130_gp,  814.80800_gp,&
         732.15470_gp,  726.84330_gp,  703.77870_gp,  724.68840_gp,  701.80630_gp,&
         39.27420_gp,  127.54920_gp,  162.61080_gp,  124.11700_gp,   93.38480_gp,&
         64.84680_gp,   47.18160_gp,   32.43210_gp,  186.46690_gp,  288.33710_gp,&
         292.02470_gp,  234.48430_gp,  191.84290_gp,  162.07500_gp,  132.34800_gp,&
         256.75790_gp,  490.74810_gp,  256.88890_gp,  247.96360_gp,  243.01480_gp/)
    vdwparams%coeffs( 6901: 7000)=(/  240.78340_gp,  221.53400_gp,  204.99250_gp,   23.23100_gp,   15.35450_gp,&
         362.37460_gp,  209.22270_gp,  140.74280_gp,   95.00510_gp,   66.42180_gp,&
         50.30160_gp,   38.13130_gp,   29.35980_gp,  433.37360_gp,  333.47670_gp,&
         307.25400_gp,  241.91590_gp,  188.48520_gp,  156.34230_gp,  127.65130_gp,&
         104.40910_gp,  710.18180_gp,  585.95780_gp,  484.07230_gp,  467.39840_gp,&
         427.99360_gp,  337.14360_gp,  368.47190_gp,  289.19370_gp,  306.49680_gp,&
         315.77130_gp,  242.10300_gp,  248.48550_gp,  294.24080_gp,  259.32240_gp,&
         221.02330_gp,  198.22130_gp,  173.37930_gp,  150.71900_gp,  795.88590_gp,&
         698.07620_gp,  611.34030_gp,  549.32940_gp,  500.83550_gp,  386.52450_gp,&
         431.34780_gp,  328.50150_gp,  359.10570_gp,  332.98960_gp,  277.68230_gp,&
         293.40740_gp,  368.32240_gp,  340.42650_gp,  303.18560_gp,  281.16330_gp,&
         254.08770_gp,  228.32560_gp,  969.54620_gp,  889.93460_gp,  781.63500_gp,&
         358.16770_gp,  788.58090_gp,  757.17110_gp,  738.17480_gp,  720.70520_gp,&
         705.21290_gp,  554.34620_gp,  624.44810_gp,  602.26720_gp,  636.03570_gp,&
         622.53720_gp,  610.36000_gp,  603.23520_gp,  508.16310_gp,  499.47010_gp,&
         455.68660_gp,  384.67230_gp,  391.31100_gp,  354.55820_gp,  324.64910_gp,&
         269.63380_gp,  251.91510_gp,  259.13820_gp,  378.09130_gp,  369.55020_gp,&
         339.39390_gp,  323.57100_gp,  298.47990_gp,  273.45010_gp,  915.16100_gp,&
         879.83390_gp,  777.61650_gp,  698.56190_gp,  693.71430_gp,  671.70070_gp,&
         691.79850_gp,  669.91580_gp,   37.44960_gp,  121.60370_gp,  155.05540_gp/)
    vdwparams%coeffs( 7001: 7100)=(/  118.35300_gp,   89.09920_gp,   61.92560_gp,   45.10020_gp,   31.05160_gp,&
         177.83550_gp,  274.98400_gp,  278.41730_gp,  223.55480_gp,  182.95110_gp,&
         154.62430_gp,  126.33190_gp,  245.01880_gp,  468.68990_gp,  245.00350_gp,&
         236.53250_gp,  231.81800_gp,  229.70220_gp,  211.30250_gp,  195.54300_gp,&
         186.54620_gp,   22.68220_gp,   14.98330_gp,  355.96060_gp,  204.90860_gp,&
         137.62540_gp,   92.81070_gp,   64.85190_gp,   49.09960_gp,   37.21580_gp,&
         28.65540_gp,  425.64550_gp,  326.78640_gp,  300.82710_gp,  236.59520_gp,&
         184.18470_gp,  152.70360_gp,  124.63170_gp,  101.91160_gp,  697.76500_gp,&
         574.75220_gp,  474.62540_gp,  458.11360_gp,  419.40370_gp,  330.38080_gp,&
         360.96920_gp,  283.30110_gp,  300.10550_gp,  309.25430_gp,  237.12000_gp,&
         243.19600_gp,  288.00730_gp,  253.61340_gp,  216.00430_gp,  193.64320_gp,&
         169.31330_gp,  147.14170_gp,  781.84410_gp,  684.77150_gp,  599.26520_gp,&
         538.24800_gp,  490.59730_gp,  378.45940_gp,  422.41440_gp,  321.55060_gp,&
         351.49860_gp,  325.88580_gp,  271.80880_gp,  287.09580_gp,  360.57310_gp,&
         333.03010_gp,  296.40640_gp,  274.77510_gp,  248.22370_gp,  222.98560_gp,&
         952.34010_gp,  873.16220_gp,  766.36770_gp,  350.20540_gp,  773.55490_gp,&
         742.67570_gp,  724.02070_gp,  706.86710_gp,  691.65380_gp,  543.27420_gp,&
         612.63290_gp,  590.76940_gp,  623.69960_gp,  610.45110_gp,  598.49350_gp,&
         591.52550_gp,  498.02530_gp,  489.15250_gp,  446.09160_gp,  376.46650_gp,&
         382.90120_gp,  346.82510_gp,  317.48840_gp,  263.63220_gp,  246.29020_gp/)
    vdwparams%coeffs( 7101: 7200)=(/  253.31490_gp,  370.09330_gp,  361.53570_gp,  331.83930_gp,  316.26530_gp,&
         291.63450_gp,  267.09740_gp,  898.47060_gp,  862.94020_gp,  762.20730_gp,&
         684.24210_gp,  679.74920_gp,  658.17270_gp,  678.18560_gp,  656.68000_gp,&
         36.59010_gp,  119.07250_gp,  151.77980_gp,  115.70000_gp,   87.04320_gp,&
         60.46080_gp,   44.01730_gp,   30.30040_gp,  174.20010_gp,  269.37790_gp,&
         272.54090_gp,  218.60660_gp,  178.77350_gp,  151.02860_gp,  123.34570_gp,&
         239.88490_gp,  459.68120_gp,  239.65630_gp,  231.36310_gp,  226.76010_gp,&
         224.71760_gp,  206.61260_gp,  191.17980_gp,  182.38760_gp,  178.33110_gp,&
         22.15500_gp,   14.51470_gp,  367.43080_gp,  205.93130_gp,  136.37050_gp,&
         91.10250_gp,   63.27090_gp,   47.72390_gp,   36.07490_gp,   27.72820_gp,&
         438.71980_gp,  330.00180_gp,  301.46500_gp,  234.75480_gp,  181.32670_gp,&
         149.64210_gp,  121.63290_gp,   99.13980_gp,  722.00020_gp,  585.37040_gp,&
         481.64720_gp,  463.39240_gp,  423.41930_gp,  333.50140_gp,  363.41710_gp,&
         285.10180_gp,  300.73560_gp,  310.51330_gp,  238.12070_gp,  242.68930_gp,&
         287.85200_gp,  251.53860_gp,  212.83460_gp,  190.07710_gp,  165.58600_gp,&
         143.45250_gp,  807.87860_gp,  697.86230_gp,  606.85850_gp,  542.94990_gp,&
         493.63630_gp,  379.22990_gp,  423.91540_gp,  321.25020_gp,  351.12220_gp,&
         325.04660_gp,  271.46340_gp,  285.81800_gp,  360.66260_gp,  331.00630_gp,&
         292.85890_gp,  270.55260_gp,  243.54750_gp,  218.09290_gp,  983.59160_gp,&
         891.68850_gp,  777.70700_gp,  346.43740_gp,  788.37610_gp,  756.15770_gp/)
    vdwparams%coeffs( 7201: 7300)=(/  736.95090_gp,  719.31180_gp,  703.65920_gp,  548.91190_gp,  625.23960_gp,&
         602.14160_gp,  633.52600_gp,  619.95970_gp,  607.66330_gp,  600.75720_gp,&
         503.38410_gp,  491.22120_gp,  446.29550_gp,  375.60840_gp,  381.45840_gp,&
         344.44290_gp,  314.53250_gp,  260.57850_gp,  243.22040_gp,  249.82310_gp,&
         369.71340_gp,  359.38220_gp,  328.09610_gp,  311.75430_gp,  286.48030_gp,&
         261.59230_gp,  923.44480_gp,  878.38880_gp,  771.37100_gp,  687.98510_gp,&
         685.76140_gp,  663.90850_gp,  686.94090_gp,  664.65810_gp,   35.98710_gp,&
         119.43970_gp,  151.84260_gp,  114.34660_gp,   85.47310_gp,   58.99030_gp,&
         42.74580_gp,   29.28080_gp,  175.25310_gp,  271.27550_gp,  272.63840_gp,&
         216.61540_gp,  175.96280_gp,  148.02440_gp,  120.39270_gp,  240.17910_gp,&
         467.91560_gp,  238.01280_gp,  229.72260_gp,  225.20200_gp,  223.40090_gp,&
         204.44800_gp,  188.94450_gp,  180.27190_gp,  176.33580_gp,  174.99690_gp,&
         21.21390_gp,   14.11850_gp,  317.88260_gp,  187.30120_gp,  127.24880_gp,&
         86.47090_gp,   60.73310_gp,   46.13550_gp,   35.06290_gp,   27.05240_gp,&
         380.65680_gp,  297.54400_gp,  275.62990_gp,  218.52820_gp,  171.18890_gp,&
         142.45910_gp,  116.66370_gp,   95.65860_gp,  621.48410_gp,  519.51820_gp,&
         430.36750_gp,  416.53280_gp,  381.96290_gp,  300.93230_gp,  329.51960_gp,&
         258.72900_gp,  275.03810_gp,  282.96560_gp,  216.95700_gp,  223.66220_gp,&
         264.47240_gp,  234.34070_gp,  200.64080_gp,  180.41890_gp,  158.22130_gp,&
         137.85900_gp,  697.20120_gp,  618.59280_gp,  544.31400_gp,  490.50970_gp/)
    vdwparams%coeffs( 7301: 7400)=(/  448.03280_gp,  346.83870_gp,  386.62560_gp,  295.43070_gp,  322.98310_gp,&
         299.83170_gp,  249.83480_gp,  264.57160_gp,  330.95850_gp,  307.26350_gp,&
         274.78160_gp,  255.42890_gp,  231.40060_gp,  208.40550_gp,  849.34460_gp,&
         787.27710_gp,  694.79590_gp,  324.33650_gp,  698.72430_gp,  671.47320_gp,&
         654.77760_gp,  639.40750_gp,  625.78200_gp,  494.42430_gp,  552.58030_gp,&
         533.39380_gp,  565.08580_gp,  553.16870_gp,  542.45160_gp,  536.00770_gp,&
         453.07620_gp,  447.44320_gp,  409.33690_gp,  346.23600_gp,  352.58880_gp,&
         320.19590_gp,  293.71130_gp,  244.37210_gp,  228.47850_gp,  235.24690_gp,&
         340.06610_gp,  333.55770_gp,  307.48710_gp,  293.76220_gp,  271.63520_gp,&
         249.38090_gp,  804.99020_gp,  780.26750_gp,  692.68280_gp,  625.29330_gp,&
         619.43840_gp,  599.85590_gp,  615.96070_gp,  596.81690_gp,   34.03120_gp,&
         109.01270_gp,  139.24200_gp,  107.20260_gp,   81.06860_gp,   56.61400_gp,&
         41.38620_gp,   28.62910_gp,  159.12280_gp,  245.80540_gp,  250.07300_gp,&
         202.13060_gp,  166.18870_gp,  140.88040_gp,  115.45010_gp,  219.98520_gp,&
         415.56380_gp,  221.24670_gp,  213.61910_gp,  209.34280_gp,  207.29260_gp,&
         191.31250_gp,  177.20340_gp,  169.04610_gp,  165.23610_gp,  162.92360_gp,&
         153.45280_gp,   33.25150_gp,   20.98740_gp,  602.81040_gp,  326.66720_gp,&
         211.27560_gp,  138.14970_gp,   94.13290_gp,   69.90390_gp,   52.04760_gp,&
         39.45870_gp,  717.75310_gp,  526.42360_gp,  475.73180_gp,  364.61100_gp,&
         277.35230_gp,  226.27960_gp,  181.71780_gp,  146.38550_gp, 1184.83270_gp/)
    vdwparams%coeffs( 7401: 7500)=(/  942.58280_gp,  771.97630_gp,  739.02810_gp,  673.28230_gp,  529.00060_gp,&
         575.32290_gp,  449.90400_gp,  472.64010_gp,  489.54990_gp,  374.25370_gp,&
         378.73290_gp,  451.40030_gp,  389.73700_gp,  325.72680_gp,  288.39340_gp,&
         248.82110_gp,  213.48210_gp, 1322.73220_gp, 1123.91670_gp,  969.07770_gp,&
         862.03630_gp,  780.48050_gp,  594.57650_gp,  666.76010_gp,  500.51620_gp,&
         547.83460_gp,  505.66980_gp,  421.74200_gp,  442.86390_gp,  564.15640_gp,&
         513.03150_gp,  449.30060_gp,  412.27930_gp,  368.21830_gp,  327.08380_gp,&
         1609.89510_gp, 1439.36100_gp, 1245.03820_gp,  532.06830_gp, 1268.25810_gp,&
         1214.78670_gp, 1183.49310_gp, 1154.81720_gp, 1129.36230_gp,  871.69660_gp,&
         1006.12730_gp,  967.41930_gp, 1014.74430_gp,  992.83980_gp,  972.87530_gp,&
         962.26810_gp,  800.56070_gp,  774.79490_gp,  699.54820_gp,  585.00690_gp,&
         592.90660_gp,  532.19500_gp,  483.50900_gp,  397.73910_gp,  370.15170_gp,&
         379.70920_gp,  575.62810_gp,  555.83070_gp,  502.98850_gp,  475.31340_gp,&
         433.66250_gp,  393.18500_gp, 1501.45410_gp, 1411.52170_gp, 1229.73070_gp,&
         1085.73940_gp, 1086.53800_gp, 1051.49420_gp, 1094.10310_gp, 1057.50950_gp,&
         54.99190_gp,  188.76220_gp,  238.90680_gp,  176.09760_gp,  129.56700_gp,&
         87.66450_gp,   62.35220_gp,   41.59090_gp,  277.51430_gp,  430.82970_gp,&
         428.95210_gp,  335.52060_gp,  268.93700_gp,  223.84530_gp,  179.86830_gp,&
         375.84950_gp,  751.81440_gp,  368.54170_gp,  355.34980_gp,  348.32930_gp,&
         345.97960_gp,  314.42440_gp,  289.75880_gp,  276.31300_gp,  270.40230_gp/)
    vdwparams%coeffs( 7501: 7600)=(/  269.73440_gp,  248.70030_gp,  421.64970_gp,   29.95270_gp,   19.37610_gp,&
         484.01650_gp,  276.57170_gp,  184.20550_gp,  123.05740_gp,   85.16210_gp,&
         63.92450_gp,   48.02330_gp,   36.66150_gp,  578.08820_gp,  441.52390_gp,&
         405.16110_gp,  316.97000_gp,  245.26830_gp,  202.28630_gp,  164.12660_gp,&
         133.37940_gp,  946.57560_gp,  777.53610_gp,  641.40640_gp,  618.04010_gp,&
         565.26390_gp,  444.41960_gp,  485.76480_gp,  380.38530_gp,  402.91460_gp,&
         415.64390_gp,  317.86240_gp,  325.69280_gp,  386.71290_gp,  339.24040_gp,&
         287.58470_gp,  256.85170_gp,  223.57440_gp,  193.36340_gp, 1059.60440_gp,&
         926.00250_gp,  808.59910_gp,  724.98050_gp,  659.80390_gp,  507.15030_gp,&
         566.82940_gp,  429.69960_gp,  470.30780_gp,  435.52680_gp,  362.50760_gp,&
         383.02040_gp,  482.91150_gp,  444.90480_gp,  394.58220_gp,  364.82430_gp,&
         328.47500_gp,  293.99660_gp, 1290.32490_gp, 1180.97290_gp, 1034.49000_gp,&
         466.18800_gp, 1044.86240_gp, 1002.91380_gp,  977.65230_gp,  954.44310_gp,&
         933.86150_gp,  730.92270_gp,  826.51040_gp,  796.62990_gp,  841.75260_gp,&
         823.86870_gp,  807.70310_gp,  798.43800_gp,  670.60970_gp,  657.55650_gp,&
         598.36780_gp,  503.44310_gp,  511.79750_gp,  462.50420_gp,  422.48760_gp,&
         349.51600_gp,  326.00720_gp,  335.30450_gp,  494.29340_gp,  482.15480_gp,&
         441.28680_gp,  419.74070_gp,  385.94190_gp,  352.37490_gp, 1215.44540_gp,&
         1165.71070_gp, 1027.48590_gp,  919.35170_gp,  913.90230_gp,  884.71530_gp,&
         912.97050_gp,  883.77800_gp,   48.72470_gp,  160.45020_gp,  204.16900_gp/)
    vdwparams%coeffs( 7601: 7700)=(/  154.44090_gp,  115.34880_gp,   79.32780_gp,   57.17520_gp,   38.75830_gp,&
         234.60480_gp,  363.33200_gp,  366.65610_gp,  292.53350_gp,  237.95800_gp,&
         200.05260_gp,  162.42060_gp,  321.15740_gp,  620.67350_gp,  320.13420_gp,&
         308.83040_gp,  302.62720_gp,  299.97060_gp,  275.26070_gp,  254.37350_gp,&
         242.55720_gp,  237.17370_gp,  234.83270_gp,  219.45080_gp,  362.68320_gp,&
         316.85210_gp,   27.35740_gp,   18.08160_gp,  406.34270_gp,  240.75590_gp,&
         163.96550_gp,  111.42570_gp,   78.11380_gp,   59.17640_gp,   44.80790_gp,&
         34.42260_gp,  486.46380_gp,  381.83960_gp,  354.37060_gp,  281.48930_gp,&
         220.71460_gp,  183.64620_gp,  150.27420_gp,  123.03580_gp,  793.63060_gp,&
         665.32210_gp,  551.49990_gp,  534.02660_gp,  489.82100_gp,  385.65770_gp,&
         422.68740_gp,  331.63890_gp,  353.00740_gp,  363.05660_gp,  278.05780_gp,&
         287.17780_gp,  339.96250_gp,  301.68000_gp,  258.52630_gp,  232.50720_gp,&
         203.84000_gp,  177.47160_gp,  890.49310_gp,  792.08160_gp,  697.80370_gp,&
         629.21900_gp,  574.90090_gp,  445.07990_gp,  496.10230_gp,  379.04270_gp,&
         414.54110_gp,  384.81520_gp,  320.22310_gp,  339.48140_gp,  424.70080_gp,&
         394.81650_gp,  353.42450_gp,  328.67560_gp,  297.80640_gp,  268.17000_gp,&
         1085.36770_gp, 1007.79130_gp,  890.47640_gp,  417.03500_gp,  894.28150_gp,&
         859.52480_gp,  838.18370_gp,  818.53810_gp,  801.12860_gp,  633.58090_gp,&
         706.87070_gp,  682.51260_gp,  723.60130_gp,  708.36530_gp,  694.68070_gp,&
         686.40110_gp,  580.61170_gp,  574.40900_gp,  525.73750_gp,  444.64860_gp/)
    vdwparams%coeffs( 7701: 7800)=(/  452.91770_gp,  411.36690_gp,  377.33280_gp,  313.70400_gp,  293.16560_gp,&
         301.97230_gp,  436.15840_gp,  428.26390_gp,  395.14950_gp,  377.68730_gp,&
         349.34070_gp,  320.71970_gp, 1029.53990_gp,  999.43950_gp,  887.91120_gp,&
         802.15530_gp,  793.94730_gp,  768.78210_gp,  788.64020_gp,  764.21100_gp,&
         43.90670_gp,  140.13530_gp,  179.13900_gp,  138.14460_gp,  104.45320_gp,&
         72.80730_gp,   53.06450_gp,   36.47280_gp,  204.19370_gp,  315.66880_gp,&
         321.60430_gp,  260.40440_gp,  214.25140_gp,  181.59050_gp,  148.69430_gp,&
         282.30700_gp,  532.24610_gp,  284.45160_gp,  274.59290_gp,  269.02700_gp,&
         266.29390_gp,  245.97130_gp,  227.81130_gp,  217.26350_gp,  212.32210_gp,&
         209.13220_gp,  197.27710_gp,  319.78710_gp,  282.57170_gp,  254.19240_gp,&
         25.31490_gp,   17.02800_gp,  355.04500_gp,  215.16710_gp,  148.85900_gp,&
         102.46200_gp,   72.55670_gp,   55.37520_gp,   42.20270_gp,   32.59350_gp,&
         425.77330_gp,  339.84270_gp,  317.81330_gp,  255.11530_gp,  201.96030_gp,&
         169.17120_gp,  139.34530_gp,  114.76770_gp,  693.98180_gp,  588.57630_gp,&
         489.33050_gp,  475.44530_gp,  436.92930_gp,  344.58580_gp,  378.11330_gp,&
         297.26970_gp,  317.22330_gp,  325.56670_gp,  249.82100_gp,  259.19160_gp,&
         306.12750_gp,  273.76800_gp,  236.42060_gp,  213.72980_gp,  188.40310_gp,&
         164.88330_gp,  780.09860_gp,  700.80150_gp,  620.89180_gp,  562.00010_gp,&
         514.89450_gp,  400.77220_gp,  445.79160_gp,  342.59700_gp,  374.29000_gp,&
         348.04020_gp,  289.83280_gp,  307.72140_gp,  382.83810_gp,  357.96090_gp/)
    vdwparams%coeffs( 7801: 7900)=(/  322.48120_gp,  301.15060_gp,  274.13770_gp,  247.98670_gp,  951.63050_gp,&
         890.62640_gp,  791.24710_gp,  380.31600_gp,  791.88170_gp,  761.69980_gp,&
         742.94350_gp,  725.65360_gp,  710.33760_gp,  565.79140_gp,  626.17050_gp,&
         605.28280_gp,  642.37950_gp,  628.91080_gp,  616.86660_gp,  609.31460_gp,&
         517.90690_gp,  515.23880_gp,  473.46190_gp,  402.07780_gp,  410.03520_gp,&
         373.74220_gp,  343.83760_gp,  286.99450_gp,  268.61450_gp,  276.87780_gp,&
         394.33850_gp,  388.73570_gp,  360.61850_gp,  345.84750_gp,  321.25280_gp,&
         296.13640_gp,  906.70650_gp,  885.97380_gp,  790.97170_gp,  719.20180_gp,&
         710.02990_gp,  687.65590_gp,  702.68800_gp,  681.35890_gp,   40.21170_gp,&
         125.55640_gp,  161.06370_gp,  125.88880_gp,   96.08720_gp,   67.67960_gp,&
         49.76650_gp,   34.58820_gp,  182.60620_gp,  281.92490_gp,  289.01630_gp,&
         236.42530_gp,  196.13740_gp,  167.26720_gp,  137.87350_gp,  254.53760_gp,&
         471.98210_gp,  258.12950_gp,  249.34370_gp,  244.27270_gp,  241.58030_gp,&
         224.07570_gp,  207.87480_gp,  198.29930_gp,  193.72130_gp,  190.16850_gp,&
         180.48120_gp,  288.54930_gp,  256.95230_gp,  232.57840_gp,  213.78620_gp,&
         22.62300_gp,   15.53870_gp,  297.68350_gp,  184.88320_gp,  130.18810_gp,&
         90.93430_gp,   65.15430_gp,   50.16320_gp,   38.52940_gp,   29.95010_gp,&
         357.72460_gp,  290.72760_gp,  274.18940_gp,  222.68120_gp,  178.19840_gp,&
         150.42030_gp,  124.84950_gp,  103.54670_gp,  582.89140_gp,  500.41520_gp,&
         417.36480_gp,  407.09700_gp,  374.92830_gp,  296.36070_gp,  325.49450_gp/)
    vdwparams%coeffs( 7901: 8000)=(/  256.59390_gp,  274.46200_gp,  281.01540_gp,  216.21400_gp,  225.34570_gp,&
         265.39390_gp,  239.36360_gp,  208.49570_gp,  189.59740_gp,  168.17910_gp,&
         148.06890_gp,  656.65740_gp,  596.02140_gp,  531.34260_gp,  482.98730_gp,&
         443.88480_gp,  347.68000_gp,  385.79860_gp,  298.51740_gp,  325.66500_gp,&
         303.41070_gp,  253.00950_gp,  268.94910_gp,  332.47910_gp,  312.79940_gp,&
         283.77610_gp,  266.23860_gp,  243.62960_gp,  221.53670_gp,  801.86140_gp,&
         756.62710_gp,  676.19300_gp,  334.52620_gp,  674.28760_gp,  649.12820_gp,&
         633.28350_gp,  618.65400_gp,  605.69940_gp,  486.34440_gp,  533.64850_gp,&
         516.49050_gp,  548.46290_gp,  537.00790_gp,  526.81470_gp,  520.16740_gp,&
         444.55160_gp,  444.88510_gp,  410.64300_gp,  350.42340_gp,  357.79320_gp,&
         327.43170_gp,  302.24750_gp,  253.48160_gp,  237.68320_gp,  245.14110_gp,&
         343.74040_gp,  340.24220_gp,  317.48380_gp,  305.60970_gp,  285.22290_gp,&
         264.13340_gp,  767.79060_gp,  755.26170_gp,  677.87670_gp,  620.82600_gp,&
         611.30750_gp,  592.18710_gp,  602.58350_gp,  584.69980_gp,   35.51230_gp,&
         108.20480_gp,  139.36170_gp,  110.58020_gp,   85.32650_gp,   60.83590_gp,&
         45.20190_gp,   31.83130_gp,  157.10050_gp,  242.14760_gp,  249.92480_gp,&
         206.78800_gp,  173.15820_gp,  148.72190_gp,  123.52740_gp,  221.16420_gp,&
         402.53800_gp,  225.77840_gp,  218.27110_gp,  213.82760_gp,  211.28090_gp,&
         196.83950_gp,  182.95500_gp,  174.59100_gp,  170.50040_gp,  166.78160_gp,&
         159.30490_gp,  250.72470_gp,  225.17920_gp,  205.23050_gp,  189.64080_gp/)
    vdwparams%coeffs( 8001: 8100)=(/  169.23750_gp,   35.19630_gp,   22.62310_gp,  657.40760_gp,  345.43780_gp,&
         222.50810_gp,  145.99740_gp,  100.14180_gp,   74.88680_gp,   56.17790_gp,&
         42.89570_gp,  781.85270_gp,  559.13370_gp,  502.99330_gp,  383.90840_gp,&
         292.16040_gp,  239.05720_gp,  192.78770_gp,  156.10100_gp, 1305.88400_gp,&
         1013.07060_gp,  826.10410_gp,  789.41170_gp,  718.21250_gp,  566.18560_gp,&
         612.64440_gp,  480.62360_gp,  501.59400_gp,  519.98260_gp,  399.35340_gp,&
         401.05480_gp,  477.97370_gp,  411.03890_gp,  343.39780_gp,  304.53120_gp,&
         263.45920_gp,  226.85900_gp, 1457.78170_gp, 1210.81080_gp, 1037.44340_gp,&
         920.13580_gp,  832.13680_gp,  633.94710_gp,  710.77790_gp,  533.61350_gp,&
         582.37420_gp,  537.25650_gp,  450.47140_gp,  470.31190_gp,  599.89870_gp,&
         542.87770_gp,  474.56410_gp,  435.55550_gp,  389.46840_gp,  346.67910_gp,&
         1776.88800_gp, 1556.99580_gp, 1337.50100_gp,  563.14800_gp, 1369.52450_gp,&
         1309.70140_gp, 1275.39680_gp, 1243.99500_gp, 1216.10180_gp,  934.60820_gp,&
         1093.09340_gp, 1050.09250_gp, 1090.35000_gp, 1066.44690_gp, 1044.60450_gp,&
         1033.31960_gp,  857.24780_gp,  824.20500_gp,  742.73970_gp,  622.09770_gp,&
         629.47280_gp,  564.55160_gp,  512.79660_gp,  422.90190_gp,  393.92880_gp,&
         403.26450_gp,  614.04690_gp,  589.82750_gp,  532.39070_gp,  502.86730_gp,&
         459.00050_gp,  416.73210_gp, 1647.64560_gp, 1522.36560_gp, 1317.64260_gp,&
         1157.86490_gp, 1163.49760_gp, 1125.90760_gp, 1175.35840_gp, 1135.12090_gp,&
         57.90390_gp,  199.67710_gp,  252.95740_gp,  185.72060_gp,  137.21250_gp/)
    vdwparams%coeffs( 8101: 8200)=(/   93.45840_gp,   67.00450_gp,   45.24890_gp,  294.67950_gp,  457.78310_gp,&
         453.36360_gp,  353.41420_gp,  283.43680_gp,  236.54670_gp,  190.86020_gp,&
         401.03710_gp,  812.34890_gp,  389.90570_gp,  376.39190_gp,  369.06490_gp,&
         366.86860_gp,  332.22150_gp,  306.27200_gp,  292.29920_gp,  286.12410_gp,&
         286.01310_gp,  262.73320_gp,  446.07980_gp,  382.50120_gp,  337.21600_gp,&
         304.72790_gp,  265.46470_gp,  476.25910_gp,   64.40650_gp,   38.83330_gp,&
         1618.24650_gp,  734.79260_gp,  441.12610_gp,  274.70310_gp,  181.19730_gp,&
         131.78460_gp,   96.50800_gp,   72.27030_gp, 1907.95920_gp, 1218.12740_gp,&
         1057.02180_gp,  766.38470_gp,  559.38210_gp,  445.78600_gp,  350.42730_gp,&
         277.47590_gp, 3302.25040_gp, 2314.84200_gp, 1849.73190_gp, 1740.74100_gp,&
         1568.22630_gp, 1238.50360_gp, 1318.80960_gp, 1034.06020_gp, 1052.23360_gp,&
         1101.66660_gp,  848.79520_gp,  822.51670_gp,  991.91230_gp,  818.48870_gp,&
         660.09280_gp,  573.08500_gp,  485.06500_gp,  409.37610_gp, 3670.92610_gp,&
         2782.45770_gp, 2305.60860_gp, 2003.78290_gp, 1789.76380_gp, 1335.07000_gp,&
         1509.29810_gp, 1106.12590_gp, 1203.43010_gp, 1100.53700_gp,  931.60220_gp,&
         952.75330_gp, 1250.11730_gp, 1091.66360_gp,  924.13690_gp,  832.50600_gp,&
         729.71760_gp,  637.46600_gp, 4496.98460_gp, 3631.45400_gp, 3014.07810_gp,&
         1104.43250_gp, 3162.11980_gp, 2999.92830_gp, 2915.79490_gp, 2839.35010_gp,&
         2771.32340_gp, 2059.46920_gp, 2566.87290_gp, 2457.84020_gp, 2460.94700_gp,&
         2404.11110_gp, 2351.28030_gp, 2328.93690_gp, 1892.74470_gp, 1752.95650_gp/)
    vdwparams%coeffs( 8201: 8300)=(/ 1547.98970_gp, 1280.49290_gp, 1284.04340_gp, 1131.29130_gp, 1013.10950_gp,&
         824.66090_gp,  763.88530_gp,  775.50460_gp, 1274.22830_gp, 1187.33590_gp,&
         1040.05230_gp,  966.30350_gp,  865.00840_gp,  771.78770_gp, 4053.38970_gp,&
         3491.08370_gp, 2924.81740_gp, 2481.26270_gp, 2542.55320_gp, 2457.23930_gp,&
         2618.00360_gp, 2517.50530_gp,  110.32750_gp,  421.46100_gp,  529.54940_gp,&
         363.36070_gp,  258.98850_gp,  169.17720_gp,  117.10180_gp,   75.67580_gp,&
         630.83040_gp,  988.78570_gp,  944.75350_gp,  701.00090_gp,  542.12620_gp,&
         441.47190_gp,  347.12550_gp,  841.99830_gp, 1876.94170_gp,  779.94880_gp,&
         754.67530_gp,  740.17630_gp,  739.83150_gp,  651.84050_gp,  596.97670_gp,&
         570.30240_gp,  559.38830_gp,  569.91510_gp,  505.26080_gp,  916.23120_gp,&
         755.05840_gp,  647.02200_gp,  573.71350_gp,  489.34610_gp,  992.79700_gp,&
         2312.66020_gp,   51.23560_gp,   31.98890_gp, 1074.58750_gp,  532.51100_gp,&
         334.13510_gp,  214.77600_gp,  144.92070_gp,  107.04380_gp,   79.40760_gp,&
         60.06280_gp, 1272.70280_gp,  869.37390_gp,  771.86060_gp,  578.05500_gp,&
         433.01710_gp,  350.57750_gp,  279.73550_gp,  224.31260_gp, 2163.70310_gp,&
         1603.63620_gp, 1297.29640_gp, 1232.34330_gp, 1116.95310_gp,  880.60690_gp,&
         947.58340_gp,  742.66410_gp,  768.21070_gp,  799.37100_gp,  614.09830_gp,&
         609.04990_gp,  729.58130_gp,  617.99850_gp,  509.55030_gp,  448.13810_gp,&
         384.29010_gp,  328.13160_gp, 2411.45770_gp, 1920.92950_gp, 1624.36360_gp,&
         1429.35760_gp, 1286.47700_gp,  971.59590_gp, 1093.24240_gp,  812.62860_gp/)
    vdwparams%coeffs( 8301: 8400)=(/  886.27500_gp,  814.84020_gp,  684.97960_gp,  710.21720_gp,  916.15920_gp,&
         818.41300_gp,  707.09330_gp,  644.44940_gp,  571.84320_gp,  505.22620_gp,&
         2950.24750_gp, 2485.63100_gp, 2105.74910_gp,  840.84100_gp, 2178.32290_gp,&
         2074.88070_gp, 2018.92240_gp, 1967.88270_gp, 1922.52260_gp, 1458.19690_gp,&
         1752.11440_gp, 1683.05490_gp, 1717.01230_gp, 1678.57810_gp, 1643.21390_gp,&
         1626.33390_gp, 1339.81490_gp, 1268.94920_gp, 1134.64070_gp,  945.40990_gp,&
         953.56340_gp,  849.38170_gp,  767.30030_gp,  629.16050_gp,  584.63530_gp,&
         596.94450_gp,  935.62990_gp,  888.87560_gp,  793.67010_gp,  745.14710_gp,&
         675.18420_gp,  608.86250_gp, 2698.09700_gp, 2413.42250_gp, 2062.26850_gp,&
         1787.18630_gp, 1809.28610_gp, 1749.73790_gp, 1840.69010_gp, 1774.63910_gp,&
         85.68150_gp,  306.88620_gp,  387.65640_gp,  277.40500_gp,  201.97440_gp,&
         135.18390_gp,   95.44580_gp,   63.18890_gp,  454.99280_gp,  709.69610_gp,&
         693.47320_gp,  530.81800_gp,  419.87180_gp,  346.97100_gp,  276.97730_gp,&
         613.84630_gp, 1293.51590_gp,  586.31350_gp,  566.86700_gp,  555.72150_gp,&
         553.47640_gp,  496.33310_gp,  456.38780_gp,  435.67010_gp,  426.74050_gp,&
         429.37430_gp,  389.58110_gp,  678.25900_gp,  573.21390_gp,  500.00710_gp,&
         448.48720_gp,  387.38760_gp,  726.86130_gp, 1586.54020_gp, 1132.29830_gp,&
         39.66690_gp,   25.70370_gp,  678.48100_gp,  374.01350_gp,  246.01470_gp,&
         163.45480_gp,  112.93390_gp,   84.78760_gp,   63.76800_gp,   48.76000_gp,&
         808.99730_gp,  600.55110_gp,  546.71760_gp,  423.69920_gp,  326.00870_gp/)
    vdwparams%coeffs( 8401: 8500)=(/  268.30770_gp,  217.43220_gp,  176.66900_gp, 1335.98990_gp, 1070.91040_gp,&
         879.01650_gp,  844.15650_gp,  770.41850_gp,  606.70400_gp,  660.10050_gp,&
         517.58120_gp,  544.66430_gp,  562.93440_gp,  431.55440_gp,  438.40800_gp,&
         521.18020_gp,  453.69050_gp,  382.64580_gp,  341.03380_gp,  296.39440_gp,&
         256.15010_gp, 1493.79930_gp, 1277.58130_gp, 1106.66300_gp,  987.82570_gp,&
         896.79520_gp,  687.18480_gp,  768.84280_gp,  580.90070_gp,  634.80240_gp,&
         587.02020_gp,  490.38170_gp,  515.38610_gp,  652.56540_gp,  596.89190_gp,&
         526.54260_gp,  485.61760_gp,  436.30460_gp,  389.93850_gp, 1819.04060_gp,&
         1635.36220_gp, 1420.61920_gp,  623.28410_gp, 1442.96800_gp, 1383.21900_gp,&
         1347.78040_gp, 1315.26700_gp, 1286.41060_gp,  999.36150_gp, 1146.12760_gp,&
         1102.56140_gp, 1156.87290_gp, 1131.94210_gp, 1109.30350_gp, 1096.86810_gp,&
         916.29250_gp,  891.41000_gp,  808.10660_gp,  679.04580_gp,  688.97190_gp,&
         620.88940_gp,  566.04620_gp,  468.02460_gp,  436.43220_gp,  447.96630_gp,&
         668.36720_gp,  647.81650_gp,  589.76480_gp,  559.56430_gp,  513.25740_gp,&
         467.84580_gp, 1702.85920_gp, 1607.81780_gp, 1406.12190_gp, 1249.02500_gp,&
         1247.32120_gp, 1207.33490_gp, 1251.91540_gp, 1210.69810_gp,   64.70680_gp,&
         216.71700_gp,  275.36440_gp,  206.00140_gp,  153.42290_gp,  105.31510_gp,&
         75.89940_gp,   51.51580_gp,  318.16780_gp,  493.30570_gp,  494.02450_gp,&
         390.71920_gp,  316.31690_gp,  265.40890_gp,  215.21010_gp,  435.06490_gp,&
         856.71310_gp,  429.30190_gp,  414.20300_gp,  406.02230_gp,  402.92000_gp/)
    vdwparams%coeffs( 8501: 8600)=(/  367.79370_gp,  339.62780_gp,  323.99210_gp,  316.94170_gp,  314.99080_gp,&
         292.40050_gp,  487.87590_gp,  423.20760_gp,  375.91790_gp,  341.20400_gp,&
         298.56340_gp,  517.69610_gp, 1044.34010_gp,  780.82120_gp,  568.67890_gp,&
         40.03040_gp,   25.95340_gp,  684.10400_gp,  377.96930_gp,  248.42660_gp,&
         164.99800_gp,  113.99890_gp,   85.60240_gp,   64.40100_gp,   49.26350_gp,&
         816.03170_gp,  607.00730_gp,  552.35500_gp,  427.88370_gp,  329.09470_gp,&
         270.80590_gp,  219.44620_gp,  178.31830_gp, 1343.37140_gp, 1081.99300_gp,&
         888.33560_gp,  853.04640_gp,  778.52660_gp,  613.02240_gp,  667.03050_gp,&
         523.00200_gp,  550.40920_gp,  568.90880_gp,  436.10420_gp,  442.98860_gp,&
         526.49640_gp,  458.19200_gp,  386.30910_gp,  344.24250_gp,  299.15350_gp,&
         258.53040_gp, 1501.77890_gp, 1290.48260_gp, 1118.14020_gp,  998.12240_gp,&
         906.08010_gp,  694.29830_gp,  776.73680_gp,  586.89450_gp,  641.41640_gp,&
         593.14140_gp,  495.48000_gp,  520.76290_gp,  659.28650_gp,  602.96920_gp,&
         531.74240_gp,  490.32060_gp,  440.45990_gp,  393.61210_gp, 1826.44610_gp,&
         1650.93340_gp, 1434.89560_gp,  629.49770_gp, 1456.51130_gp, 1397.02980_gp,&
         1361.31780_gp, 1328.54150_gp, 1299.44880_gp, 1009.44190_gp, 1155.39470_gp,&
         1110.69130_gp, 1168.83930_gp, 1143.68570_gp, 1120.84360_gp, 1108.27910_gp,&
         925.22660_gp,  900.68080_gp,  816.49740_gp,  685.98560_gp,  696.00400_gp,&
         627.20650_gp,  571.77930_gp,  472.76430_gp,  440.85970_gp,  452.47610_gp,&
         675.15180_gp,  654.44510_gp,  595.66400_gp,  565.06480_gp,  518.21240_gp/)
    vdwparams%coeffs( 8601: 8700)=(/  472.30580_gp, 1712.70500_gp, 1623.53760_gp, 1420.39830_gp, 1261.98670_gp,&
         1260.01470_gp, 1219.70640_gp, 1264.83780_gp, 1223.26380_gp,   65.31290_gp,&
         218.97060_gp,  278.08810_gp,  207.99430_gp,  154.87640_gp,  106.31410_gp,&
         76.63470_gp,   52.04120_gp,  321.53780_gp,  498.34860_gp,  499.07610_gp,&
         394.54030_gp,  319.30990_gp,  267.88580_gp,  217.20780_gp,  439.52280_gp,&
         864.14110_gp,  433.78650_gp,  418.30680_gp,  410.10330_gp,  406.98940_gp,&
         371.46910_gp,  343.00120_gp,  327.19950_gp,  320.09600_gp,  318.21360_gp,&
         295.27620_gp,  492.83160_gp,  427.36220_gp,  379.50970_gp,  344.41690_gp,&
         301.34630_gp,  522.94490_gp, 1052.86290_gp,  787.49500_gp,  574.45910_gp,&
         580.58240_gp,   37.54780_gp,   24.55300_gp,  628.18170_gp,  349.18970_gp,&
         231.04490_gp,  154.32880_gp,  107.13040_gp,   80.73600_gp,   60.94140_gp,&
         46.74980_gp,  749.49740_gp,  559.88560_gp,  511.08100_gp,  397.66020_gp,&
         307.13840_gp,  253.49850_gp,  206.04590_gp,  167.90140_gp, 1237.96430_gp,&
         996.01570_gp,  818.51610_gp,  787.05650_gp,  718.84690_gp,  566.43000_gp,&
         616.59930_gp,  483.85370_gp,  509.69990_gp,  526.37950_gp,  403.84000_gp,&
         410.98660_gp,  487.98580_gp,  426.07030_gp,  360.44960_gp,  321.94150_gp,&
         280.47020_gp,  242.96840_gp, 1385.12280_gp, 1188.19820_gp, 1031.47450_gp,&
         922.06600_gp,  837.98490_gp,  643.47210_gp,  719.36590_gp,  544.80750_gp,&
         595.15030_gp,  550.75570_gp,  460.23670_gp,  484.02860_gp,  611.39110_gp,&
         560.51040_gp,  495.69480_gp,  457.93020_gp,  412.22750_gp,  369.15310_gp/)
    vdwparams%coeffs( 8701: 8800)=(/ 1687.74900_gp, 1520.13170_gp, 1323.25830_gp,  586.64300_gp, 1342.81160_gp,&
         1287.27630_gp, 1254.40390_gp, 1224.23140_gp, 1197.45640_gp,  932.83320_gp,&
         1066.70770_gp, 1027.10260_gp, 1077.42050_gp, 1054.24840_gp, 1033.23640_gp,&
         1021.53110_gp,  855.24570_gp,  833.45230_gp,  756.75320_gp,  636.90770_gp,&
         646.54950_gp,  583.52570_gp,  532.65790_gp,  441.19610_gp,  411.71360_gp,&
         422.71910_gp,  626.93260_gp,  608.66200_gp,  555.32550_gp,  527.59970_gp,&
         484.78900_gp,  442.66740_gp, 1581.54490_gp, 1496.20770_gp, 1311.20880_gp,&
         1167.63090_gp, 1164.92880_gp, 1127.69270_gp, 1167.63830_gp, 1129.48300_gp,&
         60.97730_gp,  202.51730_gp,  257.60880_gp,  193.75570_gp,  144.87380_gp,&
         99.93810_gp,   72.35200_gp,   49.41660_gp,  297.17130_gp,  460.42340_gp,&
         462.16780_gp,  366.95490_gp,  298.06610_gp,  250.75730_gp,  203.93970_gp,&
         407.59060_gp,  797.47730_gp,  403.15400_gp,  389.17870_gp,  381.48040_gp,&
         378.44740_gp,  346.04230_gp,  319.76370_gp,  305.08870_gp,  298.41290_gp,&
         296.19890_gp,  275.60550_gp,  457.16410_gp,  397.75920_gp,  354.18680_gp,&
         322.09780_gp,  282.49330_gp,  485.41290_gp,  972.12330_gp,  730.14460_gp,&
         534.01080_gp,  539.35320_gp,  502.00730_gp,   41.33850_gp,   26.62170_gp,&
         763.78140_gp,  402.77270_gp,  260.37960_gp,  171.25810_gp,  117.64790_gp,&
         88.06090_gp,   66.10930_gp,   50.50870_gp,  908.40250_gp,  651.25640_gp,&
         586.98620_gp,  449.11910_gp,  342.48830_gp,  280.56130_gp,  226.48330_gp,&
         183.52050_gp, 1518.87550_gp, 1178.33020_gp,  961.36570_gp,  919.30990_gp/)
    vdwparams%coeffs( 8801: 8900)=(/  836.73410_gp,  659.72010_gp,  714.18300_gp,  560.36480_gp,  585.29350_gp,&
         606.47400_gp,  465.81470_gp,  468.45290_gp,  558.15220_gp,  480.89260_gp,&
         402.45400_gp,  357.25110_gp,  309.35300_gp,  266.57770_gp, 1696.12140_gp,&
         1408.38950_gp, 1207.98680_gp, 1072.19350_gp,  970.20700_gp,  739.79150_gp,&
         829.24080_gp,  623.13860_gp,  680.05240_gp,  627.58100_gp,  526.09770_gp,&
         549.63380_gp,  700.42290_gp,  634.78310_gp,  555.77940_gp,  510.56820_gp,&
         456.97070_gp,  407.09700_gp, 2068.85570_gp, 1811.06260_gp, 1557.06740_gp,&
         659.23480_gp, 1593.63310_gp, 1523.87140_gp, 1483.98170_gp, 1447.46650_gp,&
         1415.03580_gp, 1088.98390_gp, 1272.72350_gp, 1223.27390_gp, 1268.91740_gp,&
         1241.11510_gp, 1215.72980_gp, 1202.53140_gp,  998.86140_gp,  961.31920_gp,&
         867.00540_gp,  726.70660_gp,  735.57240_gp,  660.18700_gp,  600.02740_gp,&
         495.11760_gp,  461.30260_gp,  472.42330_gp,  717.22290_gp,  689.64350_gp,&
         623.35390_gp,  589.27250_gp,  538.37470_gp,  489.18620_gp, 1918.40290_gp,&
         1771.56390_gp, 1534.67920_gp, 1350.24000_gp, 1356.05040_gp, 1312.21210_gp,&
         1368.66100_gp, 1321.95940_gp,   67.89240_gp,  232.95390_gp,  295.38160_gp,&
         217.47620_gp,  160.92230_gp,  109.77550_gp,   78.79200_gp,   53.28800_gp,&
         343.55920_gp,  533.70850_gp,  529.30430_gp,  413.59910_gp,  332.27830_gp,&
         277.60020_gp,  224.20950_gp,  468.09770_gp,  945.99430_gp,  455.91090_gp,&
         440.25610_gp,  431.63910_gp,  428.96620_gp,  388.90660_gp,  358.65620_gp,&
         342.29790_gp,  335.02870_gp,  334.56150_gp,  307.88150_gp,  521.14410_gp/)
    vdwparams%coeffs( 8901: 9000)=(/  447.71310_gp,  395.24720_gp,  357.47370_gp,  311.69010_gp,  556.12550_gp,&
         1156.26040_gp,  848.37270_gp,  605.33890_gp,  611.27640_gp,  567.76700_gp,&
         649.71430_gp,   33.18750_gp,   22.06160_gp,  517.33300_gp,  297.51860_gp,&
         200.32990_gp,  135.55170_gp,   95.01900_gp,   72.11730_gp,   54.77810_gp,&
         42.24320_gp,  618.66690_gp,  474.41800_gp,  437.02090_gp,  344.21130_gp,&
         268.54160_gp,  223.08830_gp,  182.47210_gp,  149.52220_gp, 1015.79970_gp,&
         835.26920_gp,  689.60880_gp,  665.81430_gp,  609.61590_gp,  480.63820_gp,&
         524.77480_gp,  412.24680_gp,  436.40540_gp,  449.58060_gp,  345.10430_gp,&
         353.78800_gp,  418.85730_gp,  369.17340_gp,  314.92530_gp,  282.73020_gp,&
         247.62030_gp,  215.56830_gp, 1138.52650_gp,  995.61160_gp,  871.26000_gp,&
         782.71660_gp,  713.64770_gp,  551.12400_gp,  614.81850_gp,  468.56950_gp,&
         511.83150_gp,  474.63690_gp,  396.27420_gp,  418.25790_gp,  524.82360_gp,&
         484.88550_gp,  431.98650_gp,  400.83190_gp,  362.53740_gp,  326.11570_gp,&
         1387.12420_gp, 1270.10020_gp, 1114.59620_gp,  510.58640_gp, 1125.11270_gp,&
         1080.16280_gp, 1052.99740_gp, 1028.01580_gp, 1005.85860_gp,  790.51760_gp,&
         891.87640_gp,  859.97290_gp,  906.92140_gp,  887.62020_gp,  870.20630_gp,&
         860.03450_gp,  724.29040_gp,  711.47790_gp,  649.14710_gp,  548.42120_gp,&
         557.75610_gp,  505.49330_gp,  462.98700_gp,  384.91460_gp,  359.74270_gp,&
         369.89520_gp,  539.29820_gp,  526.77380_gp,  483.81430_gp,  461.39360_gp,&
         425.87210_gp,  390.47110_gp, 1308.69040_gp, 1255.35930_gp, 1108.49340_gp/)
    vdwparams%coeffs( 9001: 9100)=(/  995.48230_gp,  989.13640_gp,  957.77810_gp,  986.73520_gp,  955.42490_gp,&
         53.36480_gp,  172.97080_gp,  220.67460_gp,  168.59440_gp,  127.21250_gp,&
         88.66130_gp,   64.74480_gp,   44.70870_gp,  253.09000_gp,  391.36060_gp,&
         396.06340_gp,  318.17950_gp,  260.70120_gp,  220.64770_gp,  180.59160_gp,&
         349.47530_gp,  668.76460_gp,  349.07020_gp,  337.05530_gp,  330.35670_gp,&
         327.35530_gp,  301.03970_gp,  278.65260_gp,  265.87520_gp,  259.94720_gp,&
         256.94800_gp,  240.91240_gp,  393.35710_gp,  345.38720_gp,  309.60350_gp,&
         282.84600_gp,  249.34040_gp,  417.08310_gp,  814.19750_gp,  621.30580_gp,&
         461.91410_gp,  466.53580_gp,  435.10620_gp,  488.32640_gp,  379.26550_gp,&
         31.10240_gp,   20.77760_gp,  482.43420_gp,  277.41740_gp,  187.13270_gp,&
         126.90850_gp,   89.16590_gp,   67.81140_gp,   51.61170_gp,   39.87620_gp,&
         577.01810_gp,  442.33810_gp,  407.69990_gp,  321.45940_gp,  251.14770_gp,&
         208.90670_gp,  171.11870_gp,  140.42730_gp,  948.51400_gp,  779.05470_gp,&
         643.18930_gp,  621.18070_gp,  568.83480_gp,  448.75270_gp,  489.79330_gp,&
         385.01530_gp,  407.45740_gp,  419.67090_gp,  322.39880_gp,  330.47130_gp,&
         391.05080_gp,  344.91000_gp,  294.54210_gp,  264.67170_gp,  232.05800_gp,&
         202.25670_gp, 1063.34770_gp,  928.82550_gp,  812.93380_gp,  730.49700_gp,&
         666.22590_gp,  514.91280_gp,  574.24560_gp,  438.04260_gp,  478.28480_gp,&
         443.62800_gp,  370.63780_gp,  391.06750_gp,  490.33250_gp,  453.17940_gp,&
         404.03840_gp,  375.13030_gp,  339.56200_gp,  305.72010_gp, 1295.85070_gp/)
    vdwparams%coeffs( 9101: 9200)=(/ 1185.17470_gp, 1040.10660_gp,  477.59450_gp, 1050.07630_gp, 1008.05210_gp,&
         982.68750_gp,  959.35980_gp,  938.66880_gp,  738.13170_gp,  832.96700_gp,&
         803.25110_gp,  846.31030_gp,  828.28260_gp,  812.02140_gp,  802.50260_gp,&
         676.13200_gp,  664.19530_gp,  606.24630_gp,  512.56560_gp,  521.30610_gp,&
         472.68090_gp,  433.12920_gp,  360.41730_gp,  336.97140_gp,  346.45260_gp,&
         504.24040_gp,  492.56890_gp,  452.64510_gp,  431.85880_gp,  398.87680_gp,&
         365.99140_gp, 1222.44850_gp, 1171.53780_gp, 1034.56330_gp,  929.53190_gp,&
         923.66730_gp,  894.41350_gp,  921.28530_gp,  892.06410_gp,   49.90900_gp,&
         161.34990_gp,  205.95300_gp,  157.59490_gp,  119.12650_gp,   83.22270_gp,&
         60.91590_gp,   42.20840_gp,  236.14440_gp,  365.07050_gp,  369.58140_gp,&
         297.23390_gp,  243.84390_gp,  206.62580_gp,  169.35840_gp,  326.59470_gp,&
         624.28500_gp,  326.25030_gp,  315.09820_gp,  308.84820_gp,  306.03400_gp,&
         281.51340_gp,  260.65670_gp,  248.73880_gp,  243.19090_gp,  240.32600_gp,&
         225.43980_gp,  367.31970_gp,  322.78160_gp,  289.59340_gp,  264.78150_gp,&
         233.66080_gp,  389.91830_gp,  760.24370_gp,  580.41450_gp,  431.82880_gp,&
         436.14400_gp,  406.93380_gp,  456.56180_gp,  354.92820_gp,  332.24250_gp,&
         28.96610_gp,   19.46890_gp,  443.91060_gp,  256.23270_gp,  173.43600_gp,&
         118.01370_gp,   83.17100_gp,   63.41250_gp,   48.38240_gp,   37.46470_gp,&
         531.16760_gp,  408.30220_gp,  376.87760_gp,  297.81740_gp,  233.21060_gp,&
         194.33730_gp,  159.49140_gp,  131.13380_gp,  873.29840_gp,  718.49830_gp/)
    vdwparams%coeffs( 9201: 9300)=(/  593.49940_gp,  573.59720_gp,  525.47320_gp,  414.79180_gp,  452.73530_gp,&
         356.13440_gp,  376.99280_gp,  388.12000_gp,  298.38940_gp,  306.06330_gp,&
         361.87800_gp,  319.69460_gp,  273.50100_gp,  246.09210_gp,  216.09450_gp,&
         188.63390_gp,  979.40320_gp,  856.71690_gp,  750.59330_gp,  674.99290_gp,&
         615.97590_gp,  476.70900_gp,  531.37000_gp,  405.93980_gp,  443.06590_gp,&
         411.13610_gp,  343.67070_gp,  362.64320_gp,  454.05760_gp,  420.12680_gp,&
         375.09950_gp,  348.60740_gp,  315.92880_gp,  284.79640_gp, 1193.70400_gp,&
         1093.00440_gp,  960.13460_gp,  443.36290_gp,  968.87580_gp,  930.22500_gp,&
         906.85110_gp,  885.34770_gp,  866.27480_gp,  682.20730_gp,  768.79040_gp,&
         741.51490_gp,  781.20400_gp,  764.57050_gp,  749.57790_gp,  740.73920_gp,&
         624.71570_gp,  614.25310_gp,  561.15170_gp,  474.94930_gp,  483.15690_gp,&
         438.46770_gp,  402.08340_gp,  334.98540_gp,  313.34830_gp,  322.18840_gp,&
         467.35710_gp,  456.86740_gp,  420.32740_gp,  401.33490_gp,  371.07280_gp,&
         340.84580_gp, 1126.95490_gp, 1081.04050_gp,  955.53590_gp,  859.69500_gp,&
         853.93730_gp,  826.94450_gp,  851.19160_gp,  824.29470_gp,   46.34890_gp,&
         149.12090_gp,  190.48260_gp,  146.20200_gp,  110.79460_gp,   77.64790_gp,&
         57.00490_gp,   39.66360_gp,  218.23330_gp,  337.22560_gp,  341.79460_gp,&
         275.49250_gp,  226.46050_gp,  192.21720_gp,  157.85240_gp,  302.46350_gp,&
         576.13740_gp,  302.48980_gp,  292.21350_gp,  286.42780_gp,  283.77850_gp,&
         261.26350_gp,  242.01470_gp,  230.97880_gp,  225.81650_gp,  223.01540_gp/)
    vdwparams%coeffs( 9301: 9400)=(/  209.45460_gp,  340.06750_gp,  299.32800_gp,  268.94270_gp,  246.19370_gp,&
         217.57380_gp,  361.30470_gp,  701.65500_gp,  536.85660_gp,  400.42090_gp,&
         404.43000_gp,  377.55350_gp,  423.12790_gp,  329.68260_gp,  308.70590_gp,&
         286.95190_gp,   27.95180_gp,   18.77310_gp,  432.60780_gp,  248.65730_gp,&
         167.83500_gp,  113.99930_gp,   80.26630_gp,   61.17480_gp,   46.67260_gp,&
         36.14860_gp,  517.56800_gp,  396.61250_gp,  365.52690_gp,  288.28630_gp,&
         225.39450_gp,  187.65970_gp,  153.90320_gp,  126.48250_gp,  851.00880_gp,&
         698.86180_gp,  576.96360_gp,  557.28800_gp,  510.36440_gp,  402.83830_gp,&
         439.50990_gp,  345.69950_gp,  365.69250_gp,  376.62810_gp,  289.55800_gp,&
         296.67280_gp,  350.80100_gp,  309.45430_gp,  264.39570_gp,  237.72170_gp,&
         208.60440_gp,  182.00160_gp,  954.13200_gp,  833.31680_gp,  729.33160_gp,&
         655.43290_gp,  597.85280_gp,  462.33830_gp,  515.49530_gp,  393.50880_gp,&
         429.51430_gp,  398.46620_gp,  333.15830_gp,  351.37120_gp,  440.27330_gp,&
         406.90290_gp,  362.87400_gp,  337.01310_gp,  305.21220_gp,  274.97440_gp,&
         1162.55200_gp, 1063.35440_gp,  933.16800_gp,  428.99800_gp,  942.35470_gp,&
         904.66880_gp,  881.90710_gp,  860.96980_gp,  842.39660_gp,  662.60290_gp,&
         747.72470_gp,  721.01030_gp,  759.50840_gp,  743.32210_gp,  728.72030_gp,&
         720.16540_gp,  606.83130_gp,  595.98520_gp,  544.09850_gp,  460.26700_gp,&
         468.10840_gp,  424.58230_gp,  389.18700_gp,  324.12400_gp,  303.15550_gp,&
         311.63540_gp,  453.04550_gp,  442.50950_gp,  406.70840_gp,  388.10010_gp/)
    vdwparams%coeffs( 9401: 9500)=(/  358.59540_gp,  329.19900_gp, 1096.83690_gp, 1051.15150_gp,  928.33020_gp,&
         834.31790_gp,  829.19230_gp,  802.98090_gp,  827.16130_gp,  800.93320_gp,&
         44.78350_gp,  144.65630_gp,  184.65300_gp,  141.40600_gp,  107.02860_gp,&
         74.93540_gp,   54.98560_gp,   38.25390_gp,  211.84190_gp,  327.35090_gp,&
         331.38320_gp,  266.60080_gp,  218.86170_gp,  185.61970_gp,  152.32700_gp,&
         293.29110_gp,  560.12980_gp,  292.90630_gp,  282.93320_gp,  277.35250_gp,&
         274.84520_gp,  252.82980_gp,  234.15110_gp,  223.48050_gp,  218.50750_gp,&
         215.95810_gp,  202.56850_gp,  329.53870_gp,  289.64110_gp,  259.97060_gp,&
         237.82380_gp,  210.04110_gp,  350.22950_gp,  682.18170_gp,  520.87240_gp,&
         387.73360_gp,  391.66370_gp,  365.52610_gp,  410.06670_gp,  318.96650_gp,&
         298.66780_gp,  277.60230_gp,  268.61300_gp,   27.92810_gp,   18.86070_gp,&
         408.44290_gp,  242.09500_gp,  165.65870_gp,  113.43690_gp,   80.24920_gp,&
         61.32340_gp,   46.86920_gp,   36.34010_gp,  489.51290_gp,  384.17280_gp,&
         356.86540_gp,  284.21390_gp,  223.77830_gp,  187.01630_gp,  153.85920_gp,&
         126.73560_gp,  800.25930_gp,  670.26330_gp,  555.66790_gp,  538.55020_gp,&
         494.22370_gp,  389.99530_gp,  426.85790_gp,  335.76420_gp,  356.91720_gp,&
         366.86680_gp,  281.86750_gp,  290.78950_gp,  343.29430_gp,  305.14680_gp,&
         262.28400_gp,  236.57790_gp,  208.21130_gp,  182.08670_gp,  898.53420_gp,&
         798.45280_gp,  703.79920_gp,  635.13700_gp,  580.84120_gp,  450.97760_gp,&
         502.11680_gp,  384.94360_gp,  420.37460_gp,  390.57700_gp,  325.96620_gp/)
    vdwparams%coeffs( 9501: 9600)=(/  345.05990_gp,  430.36070_gp,  400.34740_gp,  359.05550_gp,  334.49170_gp,&
         303.84100_gp,  274.44070_gp, 1094.99430_gp, 1016.18020_gp,  898.20920_gp,&
         423.89100_gp,  902.62940_gp,  867.59380_gp,  846.06160_gp,  826.22610_gp,&
         808.64130_gp,  640.73930_gp,  714.62060_gp,  690.04230_gp,  730.42860_gp,&
         715.01980_gp,  701.18290_gp,  692.75040_gp,  586.68300_gp,  580.27960_gp,&
         531.79240_gp,  450.93600_gp,  459.35930_gp,  417.90820_gp,  383.96000_gp,&
         320.36400_gp,  299.86260_gp,  308.70590_gp,  443.21180_gp,  435.19640_gp,&
         402.09000_gp,  384.76220_gp,  356.59780_gp,  328.18400_gp, 1039.19200_gp,&
         1008.17920_gp,  896.30090_gp,  811.14820_gp,  803.15200_gp,  777.87630_gp,&
         797.75990_gp,  773.11930_gp,   44.50060_gp,  141.09010_gp,  180.52070_gp,&
         139.87930_gp,  106.43320_gp,   74.89250_gp,   55.13920_gp,   38.49060_gp,&
         205.96470_gp,  317.90020_gp,  324.08300_gp,  263.14380_gp,  217.31990_gp,&
         184.95220_gp,  152.26460_gp,  286.27350_gp,  537.08890_gp,  288.38750_gp,&
         278.58690_gp,  273.03310_gp,  270.28900_gp,  249.82010_gp,  231.61700_gp,&
         221.02480_gp,  216.02110_gp,  212.74680_gp,  200.82240_gp,  323.10800_gp,&
         286.05510_gp,  257.96890_gp,  236.66040_gp,  209.60830_gp,  342.34100_gp,&
         653.27820_gp,  505.67470_gp,  381.34050_gp,  385.12060_gp,  359.90590_gp,&
         401.23820_gp,  315.33650_gp,  295.30360_gp,  274.61540_gp,  265.54500_gp,&
         263.45080_gp,   42.28270_gp,   27.02360_gp,  769.25140_gp,  413.49650_gp,&
         267.47500_gp,  175.45850_gp,  120.12790_gp,   89.65200_gp,   67.12320_gp/)
    vdwparams%coeffs( 9601: 9700)=(/   51.17470_gp,  915.94130_gp,  667.28440_gp,  602.41430_gp,  461.50930_gp,&
         351.51570_gp,  287.39220_gp,  231.45540_gp,  187.09160_gp, 1517.36870_gp,&
         1198.92870_gp,  980.84640_gp,  938.76410_gp,  855.07440_gp,  672.85260_gp,&
         730.52210_gp,  572.17290_gp,  599.82550_gp,  621.33050_gp,  476.01650_gp,&
         480.59720_gp,  572.26580_gp,  493.79840_gp,  413.01260_gp,  366.15190_gp,&
         316.52310_gp,  272.22880_gp, 1694.20910_gp, 1430.66400_gp, 1231.69790_gp,&
         1095.04190_gp,  991.40130_gp,  755.89180_gp,  847.35570_gp,  636.74800_gp,&
         696.13210_gp,  642.64560_gp,  537.20580_gp,  563.03120_gp,  716.79570_gp,&
         651.14030_gp,  570.30150_gp,  523.61310_gp,  468.16800_gp,  416.50530_gp,&
         2062.62730_gp, 1834.25950_gp, 1583.87920_gp,  675.79110_gp, 1615.98820_gp,&
         1547.24920_gp, 1507.22360_gp, 1470.54920_gp, 1437.98440_gp, 1109.22500_gp,&
         1284.56480_gp, 1234.85820_gp, 1291.35130_gp, 1263.35030_gp, 1237.81160_gp,&
         1224.31430_gp, 1018.14070_gp,  983.56470_gp,  887.92510_gp,  743.40350_gp,&
         753.15120_gp,  676.23400_gp,  614.65820_gp,  506.54420_gp,  471.77480_gp,&
         483.63130_gp,  732.58080_gp,  706.44460_gp,  639.13150_gp,  604.10820_gp,&
         551.58440_gp,  500.67860_gp, 1921.00630_gp, 1797.55690_gp, 1563.72650_gp,&
         1379.63610_gp, 1382.29030_gp, 1337.78390_gp, 1393.11010_gp, 1346.26970_gp,&
         69.67540_gp,  239.04190_gp,  302.66780_gp,  223.16290_gp,  164.68120_gp,&
         111.97070_gp,   80.09860_gp,   53.93640_gp,  352.00740_gp,  546.27620_gp,&
         543.21880_gp,  424.82930_gp,  340.93960_gp,  284.33500_gp,  229.12290_gp/)
    vdwparams%coeffs( 9701: 9800)=(/  478.01410_gp,  957.98970_gp,  467.67000_gp,  451.16580_gp,  442.34200_gp,&
         439.47280_gp,  399.10790_gp,  367.95920_gp,  351.03000_gp,  343.56200_gp,&
         342.89310_gp,  315.91660_gp,  534.51270_gp,  459.62850_gp,  405.54530_gp,&
         366.36520_gp,  318.93050_gp,  567.60020_gp, 1168.50360_gp,  862.88330_gp,&
         619.65370_gp,  626.02540_gp,  581.05490_gp,  663.03800_gp,  500.16280_gp,&
         467.37710_gp,  432.99730_gp,  419.68680_gp,  411.15760_gp,  678.99200_gp,&
         39.05000_gp,   25.41050_gp,  642.41900_gp,  362.74340_gp,  240.59340_gp,&
         160.58760_gp,  111.27680_gp,   83.71200_gp,   63.08060_gp,   48.32580_gp,&
         767.09550_gp,  580.45000_gp,  531.07940_gp,  414.14030_gp,  319.93070_gp,&
         263.82770_gp,  214.17230_gp,  174.26390_gp, 1259.59640_gp, 1026.60970_gp,&
         845.54470_gp,  813.89640_gp,  743.92530_gp,  585.40690_gp,  638.75870_gp,&
         500.63040_gp,  528.99650_gp,  546.04730_gp,  418.17300_gp,  427.10250_gp,&
         506.88380_gp,  443.48920_gp,  375.36940_gp,  335.11820_gp,  291.71790_gp,&
         252.44080_gp, 1409.53860_gp, 1223.30730_gp, 1065.43420_gp,  953.93220_gp,&
         867.54750_gp,  666.41970_gp,  744.99930_gp,  564.45890_gp,  617.37910_gp,&
         571.56020_gp,  476.59660_gp,  502.55130_gp,  634.11950_gp,  582.76930_gp,&
         515.94670_gp,  476.68110_gp,  428.99630_gp,  383.94570_gp, 1716.15260_gp,&
         1561.82890_gp, 1364.42900_gp,  610.00420_gp, 1381.16770_gp, 1325.11510_gp,&
         1291.57070_gp, 1260.76080_gp, 1233.42800_gp,  963.22340_gp, 1094.07030_gp,&
         1053.97480_gp, 1111.00050_gp, 1087.29010_gp, 1065.81960_gp, 1053.67730_gp/)
    vdwparams%coeffs( 9801: 9900)=(/  883.59210_gp,  863.80750_gp,  785.17710_gp,  660.55390_gp,  671.10130_gp,&
         606.06460_gp,  553.42010_gp,  458.07960_gp,  427.39480_gp,  439.23750_gp,&
         649.47520_gp,  632.16060_gp,  577.57500_gp,  548.93280_gp,  504.41400_gp,&
         460.43530_gp, 1613.20500_gp, 1539.65490_gp, 1353.97430_gp, 1208.89360_gp,&
         1203.70430_gp, 1165.29870_gp, 1204.56860_gp, 1165.70620_gp,   63.51890_gp,&
         210.38150_gp,  267.50360_gp,  201.67580_gp,  150.60400_gp,  103.70700_gp,&
         74.92890_gp,   51.04920_gp,  308.24700_gp,  477.28480_gp,  480.37350_gp,&
         382.12550_gp,  310.42380_gp,  260.94970_gp,  211.97090_gp,  422.00550_gp,&
         819.99260_gp,  419.14320_gp,  404.43780_gp,  396.40790_gp,  393.12150_gp,&
         360.11140_gp,  332.75950_gp,  317.39780_gp,  310.42370_gp,  307.81940_gp,&
         286.94640_gp,  475.10140_gp,  414.06130_gp,  368.81960_gp,  335.26670_gp,&
         293.82390_gp,  502.49860_gp,  998.19190_gp,  754.11460_gp,  554.41070_gp,&
         559.98500_gp,  521.15630_gp,  587.93550_gp,  452.16820_gp,  422.74380_gp,&
         392.14710_gp,  379.63350_gp,  374.34300_gp,  603.09650_gp,  542.04390_gp,&
         37.11860_gp,   24.50680_gp,  569.50200_gp,  332.06650_gp,  224.24380_gp,&
         151.63340_gp,  106.05380_gp,   80.28990_gp,   60.81740_gp,   46.77750_gp,&
         681.36500_gp,  528.35450_gp,  487.98470_gp,  385.31020_gp,  300.75580_gp,&
         249.64690_gp,  203.91030_gp,  166.77960_gp, 1113.93690_gp,  925.32660_gp,&
         765.43030_gp,  739.77510_gp,  677.80140_gp,  533.72530_gp,  584.00790_gp,&
         458.21350_gp,  486.46450_gp,  500.90110_gp,  383.80100_gp,  394.84340_gp/)
    vdwparams%coeffs( 9901:10000)=(/  467.51190_gp,  412.96610_gp,  352.54270_gp,  316.39840_gp,  276.88560_gp,&
         240.75150_gp, 1248.81090_gp, 1101.96440_gp,  967.16720_gp,  870.13340_gp,&
         793.87710_gp,  613.24700_gp,  684.12920_gp,  521.50710_gp,  570.27900_gp,&
         528.99540_gp,  440.70410_gp,  466.29410_gp,  584.71100_gp,  541.52090_gp,&
         483.05960_gp,  448.33330_gp,  405.44000_gp,  364.50640_gp, 1521.23860_gp,&
         1403.57430_gp, 1235.57230_gp,  570.38460_gp, 1244.36900_gp, 1195.33940_gp,&
         1165.48110_gp, 1138.01160_gp, 1113.65740_gp,  877.23370_gp,  984.36690_gp,&
         949.71600_gp, 1005.00590_gp,  983.75110_gp,  964.60520_gp,  953.26960_gp,&
         804.12310_gp,  792.27200_gp,  723.59370_gp,  611.10290_gp,  621.95650_gp,&
         563.97010_gp,  516.67170_gp,  429.17250_gp,  400.97810_gp,  412.68870_gp,&
         600.16760_gp,  587.60720_gp,  540.48320_gp,  515.68470_gp,  476.05910_gp,&
         436.36860_gp, 1438.77300_gp, 1389.21180_gp, 1230.26140_gp, 1107.39310_gp,&
         1098.31730_gp, 1063.48330_gp, 1093.78360_gp, 1059.45920_gp,   59.77960_gp,&
         193.08010_gp,  246.35800_gp,  188.65770_gp,  142.17370_gp,   98.85830_gp,&
         71.98160_gp,   49.49700_gp,  281.99230_gp,  435.99230_gp,  442.39460_gp,&
         356.16810_gp,  291.92250_gp,  246.88440_gp,  201.78980_gp,  388.78440_gp,&
         739.83800_gp,  389.86450_gp,  376.33360_gp,  368.79080_gp,  365.29400_gp,&
         336.51780_gp,  311.48590_gp,  297.11100_gp,  290.44330_gp,  286.74360_gp,&
         269.42080_gp,  439.26920_gp,  386.42600_gp,  346.56560_gp,  316.51250_gp,&
         278.80760_gp,  463.96320_gp,  899.77960_gp,  690.24350_gp,  515.35410_gp/)
    vdwparams%coeffs(10001:10100)=(/  520.41820_gp,  485.33740_gp,  543.48220_gp,  423.37510_gp,  396.03590_gp,&
         367.74930_gp,  355.70800_gp,  352.21320_gp,  557.62720_gp,  505.19000_gp,&
         473.47030_gp,   35.48610_gp,   23.71000_gp,  519.45180_gp,  308.94260_gp,&
         211.27990_gp,  144.26000_gp,  101.63150_gp,   77.33650_gp,   58.83320_gp,&
         45.40440_gp,  622.30990_gp,  489.78200_gp,  455.24770_gp,  362.55610_gp,&
         285.12390_gp,  237.85710_gp,  195.21310_gp,  160.33310_gp, 1015.77380_gp,&
         852.90360_gp,  707.37690_gp,  685.56870_gp,  629.13900_gp,  495.88010_gp,&
         543.34380_gp,  426.83610_gp,  454.31700_gp,  467.00090_gp,  358.18700_gp,&
         370.06740_gp,  437.44170_gp,  388.89460_gp,  334.02630_gp,  300.96170_gp,&
         264.44460_gp,  230.79470_gp, 1140.33970_gp, 1015.61220_gp,  895.73040_gp,&
         808.42520_gp,  739.20970_gp,  573.37000_gp,  638.64620_gp,  489.01000_gp,&
         534.45110_gp,  496.43130_gp,  413.59130_gp,  438.35340_gp,  547.29130_gp,&
         509.38490_gp,  456.75890_gp,  425.31520_gp,  385.99880_gp,  348.21640_gp,&
         1389.95620_gp, 1292.05050_gp, 1142.80090_gp,  538.97800_gp, 1147.37180_gp,&
         1102.94110_gp, 1075.59890_gp, 1050.41680_gp, 1028.09900_gp,  814.57180_gp,&
         907.49290_gp,  876.40540_gp,  928.81650_gp,  909.26160_gp,  891.71060_gp,&
         881.00370_gp,  746.13080_gp,  738.73540_gp,  676.88650_gp,  573.38870_gp,&
         584.18960_gp,  531.22090_gp,  487.79780_gp,  406.33480_gp,  380.05290_gp,&
         391.45330_gp,  562.88150_gp,  553.06450_gp,  510.99880_gp,  488.88150_gp,&
         452.82060_gp,  416.35580_gp, 1319.62820_gp, 1282.15640_gp, 1140.33810_gp/)
    vdwparams%coeffs(10101:10200)=(/ 1031.96460_gp, 1021.11670_gp,  988.86570_gp, 1013.69010_gp,  982.42950_gp,&
         56.71560_gp,  179.98120_gp,  230.26410_gp,  178.25090_gp,  135.27200_gp,&
         94.77050_gp,   69.43170_gp,   48.10370_gp,  262.36970_gp,  405.25290_gp,&
         413.38560_gp,  335.59140_gp,  276.83720_gp,  235.20590_gp,  193.17040_gp,&
         363.84820_gp,  682.96210_gp,  366.99000_gp,  354.40290_gp,  347.26480_gp,&
         343.70870_gp,  317.77500_gp,  294.50790_gp,  280.94820_gp,  274.55800_gp,&
         270.27830_gp,  255.27130_gp,  411.70180_gp,  364.48820_gp,  328.49090_gp,&
         301.05670_gp,  266.22780_gp,  434.93910_gp,  830.36610_gp,  643.26970_gp,&
         485.08560_gp,  489.77500_gp,  457.47170_gp,  509.87810_gp,  400.53650_gp,&
         374.86560_gp,  348.38020_gp,  336.77710_gp,  334.32410_gp,  522.88630_gp,&
         476.20870_gp,  448.01730_gp,  425.10150_gp,   32.81510_gp,   22.27170_gp,&
         454.39540_gp,  276.43090_gp,  191.92730_gp,  132.62510_gp,   94.29910_gp,&
         72.23580_gp,   55.26890_gp,   42.85120_gp,  545.28550_gp,  436.44520_gp,&
         408.71490_gp,  328.81570_gp,  260.94460_gp,  219.04300_gp,  180.86040_gp,&
         149.34300_gp,  889.00260_gp,  755.34780_gp,  628.33170_gp,  610.99490_gp,&
         561.76430_gp,  443.44010_gp,  486.50250_gp,  382.89000_gp,  408.60640_gp,&
         419.15920_gp,  322.03460_gp,  334.24140_gp,  394.22890_gp,  353.11560_gp,&
         305.52140_gp,  276.61090_gp,  244.27520_gp,  214.20040_gp,  999.78710_gp,&
         899.49300_gp,  797.79010_gp,  722.72400_gp,  662.60680_gp,  516.59320_gp,&
         574.28050_gp,  442.17100_gp,  482.82170_gp,  449.20610_gp,  374.44670_gp/)
    vdwparams%coeffs(10201:10300)=(/  397.49790_gp,  493.66250_gp,  462.07160_gp,  416.87150_gp,  389.70290_gp,&
         355.21750_gp,  321.80290_gp, 1219.62720_gp, 1142.93780_gp, 1016.41550_gp,&
         491.61500_gp, 1016.95880_gp,  978.34710_gp,  954.29670_gp,  932.11660_gp,&
         912.46710_gp,  728.00530_gp,  804.53430_gp,  777.85480_gp,  825.36440_gp,&
         808.06540_gp,  792.60610_gp,  782.84050_gp,  666.14790_gp,  663.19200_gp,&
         610.01430_gp,  518.73530_gp,  529.12210_gp,  482.78010_gp,  444.56400_gp,&
         371.68100_gp,  348.13080_gp,  358.83830_gp,  509.11500_gp,  502.19660_gp,&
         466.41830_gp,  447.66490_gp,  416.30430_gp,  384.23110_gp, 1163.07020_gp,&
         1137.64700_gp, 1016.78030_gp,  925.98120_gp,  913.90750_gp,  885.20210_gp,&
         903.96880_gp,  876.65190_gp,   51.94960_gp,  161.43030_gp,  207.21850_gp,&
         162.49390_gp,  124.39580_gp,   87.98690_gp,   64.97640_gp,   45.46180_gp,&
         234.87760_gp,  362.32580_gp,  371.86150_gp,  304.87150_gp,  253.46770_gp,&
         216.58660_gp,  178.95810_gp,  328.22070_gp,  606.18590_gp,  333.18010_gp,&
         321.94330_gp,  315.43490_gp,  311.93950_gp,  289.57840_gp,  268.79340_gp,&
         256.47210_gp,  250.55360_gp,  245.84350_gp,  233.56320_gp,  371.76550_gp,&
         331.59870_gp,  300.60430_gp,  276.68620_gp,  245.86470_gp,  393.18090_gp,&
         736.90210_gp,  577.38050_gp,  440.44740_gp,  444.64020_gp,  416.12070_gp,&
         461.32230_gp,  365.92820_gp,  342.72790_gp,  318.86290_gp,  308.04750_gp,&
         306.65490_gp,  472.64410_gp,  433.05890_gp,  409.24850_gp,  389.59960_gp,&
         358.49950_gp,   42.30580_gp,   27.46950_gp,  732.91200_gp,  399.40200_gp/)
    vdwparams%coeffs(10301:10400)=(/  262.23870_gp,  174.27180_gp,  120.49820_gp,   90.53030_gp,   68.12390_gp,&
         52.10410_gp,  873.26650_gp,  642.18460_gp,  583.68500_gp,  451.61710_gp,&
         347.40300_gp,  286.04180_gp,  231.94830_gp,  188.59500_gp, 1449.58600_gp,&
         1149.82950_gp,  942.27450_gp,  904.22250_gp,  824.77800_gp,  650.06270_gp,&
         706.14120_gp,  554.09250_gp,  581.84580_gp,  601.58080_gp,  461.70190_gp,&
         467.87310_gp,  556.50160_gp,  483.71170_gp,  407.79580_gp,  363.52860_gp,&
         316.07000_gp,  273.29990_gp, 1620.76040_gp, 1372.85150_gp, 1186.38850_gp,&
         1057.76600_gp,  959.79080_gp,  735.14980_gp,  822.58290_gp,  621.19680_gp,&
         678.27670_gp,  627.00570_gp,  524.50750_gp,  550.25650_gp,  697.37450_gp,&
         636.76120_gp,  561.25330_gp,  517.58070_gp,  465.07030_gp,  415.76900_gp,&
         1975.97380_gp, 1760.06020_gp, 1524.93990_gp,  664.83460_gp, 1552.08450_gp,&
         1486.48970_gp, 1448.14640_gp, 1412.99080_gp, 1381.78630_gp, 1071.57170_gp,&
         1235.57010_gp, 1188.80120_gp, 1241.62670_gp, 1214.71810_gp, 1190.26110_gp,&
         1176.97460_gp,  982.43040_gp,  953.27210_gp,  863.45810_gp,  725.69910_gp,&
         735.87000_gp,  662.79760_gp,  604.05240_gp,  499.59380_gp,  465.88050_gp,&
         477.88590_gp,  714.77590_gp,  691.52340_gp,  628.90510_gp,  596.53870_gp,&
         547.12200_gp,  498.78710_gp, 1844.34050_gp, 1728.41170_gp, 1507.81140_gp,&
         1336.65720_gp, 1336.76730_gp, 1293.81740_gp, 1343.11150_gp, 1298.46750_gp,&
         68.95370_gp,  231.41650_gp,  294.14590_gp,  219.63870_gp,  163.67150_gp,&
         112.42780_gp,   81.08450_gp,   55.07860_gp,  340.06240_gp,  527.56300_gp/)
    vdwparams%coeffs(10401:10500)=(/  527.32820_gp,  416.48610_gp,  337.11350_gp,  282.96690_gp,  229.58370_gp,&
         465.51820_gp,  921.79340_gp,  457.90120_gp,  442.05470_gp,  433.31390_gp,&
         430.10400_gp,  392.08310_gp,  362.04550_gp,  345.44700_gp,  337.94540_gp,&
         336.08210_gp,  311.58000_gp,  520.74900_gp,  451.17170_gp,  400.64390_gp,&
         363.72020_gp,  318.38270_gp,  553.97800_gp, 1124.93450_gp,  837.61350_gp,&
         607.19510_gp,  613.20100_gp,  570.26210_gp,  647.66860_gp,  492.83010_gp,&
         460.81940_gp,  427.32470_gp,  413.77960_gp,  406.54590_gp,  661.84230_gp,&
         591.27640_gp,  549.30440_gp,  517.01090_gp,  469.48180_gp,  649.09240_gp,&
         83.01410_gp,   49.92320_gp, 2180.40150_gp,  961.10490_gp,  572.30020_gp,&
         354.84660_gp,  233.50060_gp,  169.60300_gp,  124.08020_gp,   92.84940_gp,&
         2565.93980_gp, 1599.60010_gp, 1381.43890_gp,  994.89260_gp,  723.25040_gp,&
         575.33170_gp,  451.58250_gp,  357.19350_gp, 4495.83690_gp, 3070.46240_gp,&
         2443.42770_gp, 2294.60620_gp, 2064.03090_gp, 1632.70220_gp, 1732.03750_gp,&
         1359.61560_gp, 1375.96740_gp, 1442.51330_gp, 1113.85410_gp, 1072.06430_gp,&
         1295.45630_gp, 1062.73150_gp,  853.91570_gp,  740.14960_gp,  625.55040_gp,&
         527.35210_gp, 4997.92300_gp, 3697.65690_gp, 3045.16200_gp, 2637.67820_gp,&
         2352.12780_gp, 1750.38290_gp, 1980.95080_gp, 1447.62210_gp, 1572.30230_gp,&
         1436.04000_gp, 1219.49040_gp, 1241.29800_gp, 1635.61580_gp, 1420.00220_gp,&
         1197.37900_gp, 1076.75030_gp,  942.25240_gp,  822.05680_gp, 6141.97840_gp,&
         4845.17760_gp, 3994.22730_gp, 1433.04690_gp, 4213.92250_gp, 3987.36000_gp/)
    vdwparams%coeffs(10501:10600)=(/ 3873.72840_gp, 3770.65230_gp, 3678.91460_gp, 2720.76530_gp, 3441.98520_gp,&
         3298.94220_gp, 3259.86980_gp, 3183.60110_gp, 3112.55360_gp, 3083.50020_gp,&
         2501.92620_gp, 2298.70080_gp, 2023.79910_gp, 1673.26400_gp, 1675.11360_gp,&
         1472.29330_gp, 1316.21170_gp, 1070.69110_gp,  991.39190_gp, 1004.80140_gp,&
         1668.89860_gp, 1546.10300_gp, 1348.74880_gp, 1250.85130_gp, 1117.67020_gp,&
         995.86790_gp, 5497.93930_gp, 4643.31620_gp, 3865.60410_gp, 3259.72570_gp,&
         3353.88260_gp, 3240.29130_gp, 3463.05760_gp, 3327.24780_gp,  142.56800_gp,&
         551.07920_gp,  692.84150_gp,  471.05030_gp,  334.92070_gp,  218.15960_gp,&
         150.70560_gp,   97.19260_gp,  826.95650_gp, 1298.94250_gp, 1233.58170_gp,&
         909.65670_gp,  700.99830_gp,  569.85090_gp,  447.36800_gp, 1104.10680_gp,&
         2505.65370_gp, 1013.02980_gp,  982.36380_gp,  963.27250_gp,  963.65340_gp,&
         845.21190_gp,  773.66980_gp,  739.50900_gp,  725.52030_gp,  740.98720_gp,&
         653.73360_gp, 1194.87730_gp,  979.57420_gp,  836.93640_gp,  741.02200_gp,&
         631.13970_gp, 1301.39990_gp, 3095.75940_gp, 2101.03720_gp, 1360.40000_gp,&
         1369.38250_gp, 1266.08030_gp, 1516.21340_gp, 1056.44220_gp,  986.76430_gp,&
         910.52680_gp,  885.40390_gp,  845.15800_gp, 1525.92780_gp, 1296.64080_gp,&
         1165.36030_gp, 1073.86980_gp,  951.54100_gp, 1469.07050_gp, 4171.78420_gp,&
         64.69060_gp,   40.33290_gp, 1384.12830_gp,  676.73740_gp,  423.08150_gp,&
         271.42480_gp,  182.93910_gp,  135.03220_gp,  100.10900_gp,   75.67910_gp,&
         1637.78780_gp, 1106.72210_gp,  980.43760_gp,  732.12110_gp,  547.46530_gp/)
    vdwparams%coeffs(10601:10700)=(/  442.87990_gp,  353.14670_gp,  283.03510_gp, 2800.95690_gp, 2050.42530_gp,&
         1655.70110_gp, 1571.23630_gp, 1423.11730_gp, 1122.65760_gp, 1206.13930_gp,&
         945.68080_gp,  975.98360_gp, 1016.17450_gp,  781.28810_gp,  772.66500_gp,&
         926.37230_gp,  782.74010_gp,  644.35600_gp,  566.29200_gp,  485.29550_gp,&
         414.16640_gp, 3121.55590_gp, 2458.08450_gp, 2072.79780_gp, 1821.21450_gp,&
         1637.92430_gp, 1235.62260_gp, 1390.99890_gp, 1032.58640_gp, 1125.44920_gp,&
         1034.16280_gp,  870.41020_gp,  900.76730_gp, 1164.04800_gp, 1037.30500_gp,&
         894.70260_gp,  814.81140_gp,  722.49320_gp,  637.95900_gp, 3825.08810_gp,&
         3186.29210_gp, 2690.95910_gp, 1064.59330_gp, 2790.98910_gp, 2655.11710_gp,&
         2582.97220_gp, 2517.22390_gp, 2458.78860_gp, 1860.73700_gp, 2250.96540_gp,&
         2163.45990_gp, 2193.86360_gp, 2144.47040_gp, 2098.97020_gp, 2077.57740_gp,&
         1710.35450_gp, 1614.09070_gp, 1441.33940_gp, 1200.56620_gp, 1210.08180_gp,&
         1076.75830_gp,  971.97370_gp,  796.70010_gp,  740.16750_gp,  755.24870_gp,&
         1189.16620_gp, 1127.05130_gp, 1004.58390_gp,  942.43520_gp,  853.27010_gp,&
         768.99820_gp, 3485.53860_gp, 3089.20710_gp, 2632.25610_gp, 2274.96800_gp,&
         2307.16650_gp, 2230.92670_gp, 2350.22530_gp, 2265.01220_gp,  108.30590_gp,&
         389.90700_gp,  492.58240_gp,  351.11660_gp,  255.35340_gp,  170.69380_gp,&
         120.40040_gp,   79.61510_gp,  578.70070_gp,  903.45630_gp,  880.49000_gp,&
         672.15910_gp,  530.85730_gp,  438.34960_gp,  349.67640_gp,  780.68060_gp,&
         1658.34590_gp,  742.63980_gp,  718.67820_gp,  704.47310_gp,  701.87680_gp/)
    vdwparams%coeffs(10701:10800)=(/  628.23400_gp,  577.52630_gp,  551.42260_gp,  540.16790_gp,  544.06610_gp,&
         492.63420_gp,  860.66320_gp,  725.74740_gp,  632.26020_gp,  566.75130_gp,&
         489.22910_gp,  924.22030_gp, 2036.53620_gp, 1446.23940_gp,  990.10590_gp,&
         997.92220_gp,  925.83590_gp, 1078.79080_gp,  786.58860_gp,  734.89180_gp,&
         679.65720_gp,  659.47100_gp,  639.38320_gp, 1095.45900_gp,  955.22580_gp,&
         873.22870_gp,  813.28280_gp,  729.49490_gp, 1063.32740_gp, 2705.58160_gp,&
         1850.00030_gp,   43.22020_gp,   28.50600_gp,  680.24790_gp,  390.79830_gp,&
         262.31920_gp,  176.83390_gp,  123.49550_gp,   93.43270_gp,   70.74810_gp,&
         54.40610_gp,  813.23530_gp,  623.32230_gp,  573.59020_gp,  450.92530_gp,&
         350.95780_gp,  290.93690_gp,  237.40720_gp,  194.06980_gp, 1333.79580_gp,&
         1097.16590_gp,  905.72410_gp,  874.00660_gp,  800.01520_gp,  630.19610_gp,&
         688.36930_gp,  540.23110_gp,  572.07380_gp,  589.56950_gp,  452.03420_gp,&
         463.41480_gp,  549.06920_gp,  483.31740_gp,  411.55330_gp,  368.92060_gp,&
         322.52720_gp,  280.24380_gp, 1494.40740_gp, 1307.36650_gp, 1143.55240_gp,&
         1026.82260_gp,  935.74670_gp,  721.66730_gp,  805.51110_gp,  612.96820_gp,&
         669.99170_gp,  621.06870_gp,  518.02070_gp,  546.99380_gp,  687.26670_gp,&
         634.51560_gp,  564.57530_gp,  523.32050_gp,  472.70040_gp,  424.59150_gp,&
         1820.42010_gp, 1667.46850_gp, 1462.81420_gp,  667.20180_gp, 1476.73890_gp,&
         1417.68550_gp, 1382.03200_gp, 1349.25240_gp, 1320.18090_gp, 1036.44030_gp,&
         1169.68970_gp, 1127.79990_gp, 1190.29180_gp, 1164.98250_gp, 1142.13700_gp/)
    vdwparams%coeffs(10801:10900)=(/ 1128.85480_gp,  950.04210_gp,  932.85060_gp,  850.51620_gp,  717.66720_gp,&
         729.82300_gp,  660.90480_gp,  604.87530_gp,  502.15100_gp,  469.04220_gp,&
         482.34510_gp,  705.37820_gp,  688.81010_gp,  632.03130_gp,  602.29480_gp,&
         555.31140_gp,  508.53060_gp, 1716.91080_gp, 1647.61210_gp, 1454.38810_gp,&
         1304.89700_gp, 1296.62380_gp, 1255.43590_gp, 1293.90940_gp, 1252.78700_gp,&
         69.73410_gp,  227.05790_gp,  289.44340_gp,  220.51700_gp,  165.88050_gp,&
         115.16050_gp,   83.77650_gp,   57.55990_gp,  332.14950_gp,  513.79560_gp,&
         519.61650_gp,  416.62720_gp,  340.64700_gp,  287.74490_gp,  234.95440_gp,&
         457.43380_gp,  877.55280_gp,  456.70870_gp,  440.87260_gp,  432.08390_gp,&
         428.19360_gp,  393.55560_gp,  364.11700_gp,  347.35660_gp,  339.62320_gp,&
         335.86060_gp,  314.61930_gp,  515.45120_gp,  451.95710_gp,  404.53400_gp,&
         369.06970_gp,  324.78940_gp,  545.59170_gp, 1068.06840_gp,  814.06370_gp,&
         604.12500_gp,  610.13660_gp,  568.71480_gp,  638.73760_gp,  495.20560_gp,&
         463.24120_gp,  430.08630_gp,  416.12810_gp,  411.29410_gp,  654.74970_gp,&
         591.36710_gp,  553.20050_gp,  522.90670_gp,  477.14300_gp,  644.38580_gp,&
         1385.64020_gp, 1030.61630_gp,  647.02500_gp,   42.33430_gp,   27.93240_gp,&
         665.53290_gp,  382.57600_gp,  256.86460_gp,  173.19330_gp,  120.97590_gp,&
         91.54070_gp,   69.32610_gp,   53.32010_gp,  795.68730_gp,  610.16250_gp,&
         561.54980_gp,  441.53690_gp,  343.70190_gp,  284.95290_gp,  232.55130_gp,&
         190.12240_gp, 1304.75630_gp, 1073.82350_gp,  886.52660_gp,  855.53540_gp/)
    vdwparams%coeffs(10901:11000)=(/  783.13860_gp,  616.91000_gp,  673.88620_gp,  528.87800_gp,  560.09170_gp,&
         577.19940_gp,  442.55790_gp,  453.74620_gp,  537.57900_gp,  473.26730_gp,&
         403.04410_gp,  361.32220_gp,  315.91320_gp,  274.52240_gp, 1461.90550_gp,&
         1279.52730_gp, 1119.35270_gp, 1005.17420_gp,  916.06480_gp,  706.55730_gp,&
         788.61510_gp,  600.17770_gp,  656.00920_gp,  608.12890_gp,  507.22690_gp,&
         535.62280_gp,  672.89960_gp,  621.32160_gp,  552.89340_gp,  512.52450_gp,&
         462.98250_gp,  415.89380_gp, 1780.73790_gp, 1631.85590_gp, 1431.77960_gp,&
         653.38750_gp, 1445.26110_gp, 1387.52840_gp, 1352.64490_gp, 1320.57170_gp,&
         1292.12680_gp, 1014.55930_gp, 1144.67700_gp, 1103.67880_gp, 1165.04500_gp,&
         1140.27790_gp, 1117.92340_gp, 1104.91640_gp,  929.96320_gp,  913.26380_gp,&
         832.72290_gp,  702.69660_gp,  714.61940_gp,  647.18140_gp,  592.34880_gp,&
         491.78730_gp,  459.37630_gp,  472.41230_gp,  690.65930_gp,  674.50000_gp,&
         618.95970_gp,  589.86940_gp,  543.89180_gp,  498.10570_gp, 1679.77100_gp,&
         1612.54370_gp, 1423.62140_gp, 1277.47450_gp, 1269.29090_gp, 1228.98090_gp,&
         1266.55060_gp, 1226.31940_gp,   68.29310_gp,  222.28830_gp,  283.37150_gp,&
         215.94380_gp,  162.46650_gp,  112.81340_gp,   82.08490_gp,   56.41170_gp,&
         325.16360_gp,  502.96550_gp,  508.72500_gp,  407.96290_gp,  333.60650_gp,&
         281.82650_gp,  230.14880_gp,  447.86360_gp,  858.83830_gp,  447.22200_gp,&
         431.70970_gp,  423.10550_gp,  419.28980_gp,  385.40430_gp,  356.58500_gp,&
         340.17190_gp,  332.59710_gp,  328.89540_gp,  308.12510_gp,  504.68300_gp/)
    vdwparams%coeffs(11001:11100)=(/  442.57280_gp,  396.17290_gp,  361.46790_gp,  318.12780_gp,  534.19470_gp,&
         1045.25410_gp,  796.88440_gp,  591.56300_gp,  597.46030_gp,  556.90660_gp,&
         625.39660_gp,  484.97040_gp,  453.67340_gp,  421.21350_gp,  407.54200_gp,&
         402.82810_gp,  641.08970_gp,  579.09510_gp,  541.76140_gp,  512.12090_gp,&
         467.33270_gp,  630.96990_gp, 1355.88740_gp, 1008.81760_gp,  633.63330_gp,&
         620.52060_gp,   41.03560_gp,   27.11500_gp,  642.18080_gp,  369.81420_gp,&
         248.60710_gp,  167.79680_gp,  117.30260_gp,   88.81590_gp,   67.29990_gp,&
         51.78600_gp,  767.86760_gp,  589.61410_gp,  542.96830_gp,  427.28580_gp,&
         332.86240_gp,  276.11290_gp,  225.45620_gp,  184.41100_gp, 1259.02800_gp,&
         1037.13400_gp,  856.44340_gp,  826.72570_gp,  756.88430_gp,  596.29710_gp,&
         651.44380_gp,  511.33840_gp,  541.64010_gp,  558.09040_gp,  427.96300_gp,&
         438.95460_gp,  519.95360_gp,  458.03850_gp,  390.31610_gp,  350.05540_gp,&
         306.19590_gp,  266.18930_gp, 1410.85650_gp, 1235.80860_gp, 1081.59630_gp,&
         971.56440_gp,  885.62820_gp,  683.36810_gp,  762.61340_gp,  580.65630_gp,&
         634.62930_gp,  588.39190_gp,  490.78330_gp,  518.33560_gp,  650.88900_gp,&
         601.28620_gp,  535.34260_gp,  496.42060_gp,  448.60200_gp,  403.12330_gp,&
         1718.65650_gp, 1575.93960_gp, 1383.31840_gp,  632.60690_gp, 1395.96480_gp,&
         1340.28480_gp, 1306.61230_gp, 1275.64930_gp, 1248.18960_gp,  980.61140_gp,&
         1105.65770_gp, 1066.15530_gp, 1125.54300_gp, 1101.62510_gp, 1080.04410_gp,&
         1067.45130_gp,  898.77880_gp,  883.03380_gp,  805.41670_gp,  679.86960_gp/)
    vdwparams%coeffs(11101:11200)=(/  691.47700_gp,  626.40460_gp,  573.47180_gp,  476.26530_gp,  444.93320_gp,&
         457.59140_gp,  668.21650_gp,  652.80350_gp,  599.31800_gp,  571.30760_gp,&
         526.95800_gp,  482.75550_gp, 1621.75090_gp, 1557.66030_gp, 1375.72240_gp,&
         1235.13340_gp, 1226.96120_gp, 1188.01360_gp, 1223.95390_gp, 1185.13890_gp,&
         66.14390_gp,  214.91570_gp,  274.04530_gp,  209.06270_gp,  157.40780_gp,&
         109.39440_gp,   79.65660_gp,   54.79470_gp,  314.33160_gp,  486.15820_gp,&
         491.97000_gp,  394.85010_gp,  323.09660_gp,  273.08190_gp,  223.12610_gp,&
         433.20740_gp,  829.64080_gp,  432.82760_gp,  417.83600_gp,  409.50580_gp,&
         405.78350_gp,  373.12230_gp,  345.26760_gp,  329.38120_gp,  322.03750_gp,&
         318.36420_gp,  298.41060_gp,  488.21620_gp,  428.40160_gp,  383.67750_gp,&
         350.19520_gp,  308.33610_gp,  516.80260_gp, 1009.70200_gp,  770.48950_gp,&
         572.52260_gp,  578.22130_gp,  539.06400_gp,  605.08290_gp,  469.61340_gp,&
         439.33610_gp,  407.94200_gp,  394.68140_gp,  390.21360_gp,  620.22960_gp,&
         560.54060_gp,  524.60280_gp,  496.03960_gp,  452.81350_gp,  610.65170_gp,&
         1309.60710_gp,  975.34290_gp,  613.50220_gp,  600.81010_gp,  581.74350_gp,&
         41.11860_gp,   27.09100_gp,  653.71250_gp,  373.35380_gp,  250.06900_gp,&
         168.35070_gp,  117.47400_gp,   88.83320_gp,   67.24150_gp,   51.69770_gp,&
         781.20500_gp,  596.00220_gp,  547.77420_gp,  429.94940_gp,  334.25170_gp,&
         276.91160_gp,  225.83700_gp,  184.53320_gp, 1283.63150_gp, 1050.98480_gp,&
         866.92280_gp,  836.09940_gp,  765.05120_gp,  602.71130_gp,  657.96350_gp/)
    vdwparams%coeffs(11201:11300)=(/  516.37700_gp,  546.33630_gp,  563.22550_gp,  431.90330_gp,  442.25450_gp,&
         524.16760_gp,  460.81790_gp,  392.01390_gp,  351.21820_gp,  306.89590_gp,&
         266.54940_gp, 1437.93000_gp, 1252.63420_gp, 1094.26950_gp,  981.85530_gp,&
         894.39170_gp,  689.30910_gp,  769.59900_gp,  585.20650_gp,  639.57050_gp,&
         592.71510_gp,  494.54430_gp,  521.85710_gp,  656.23190_gp,  605.19030_gp,&
         537.98990_gp,  498.42910_gp,  449.99200_gp,  404.01740_gp, 1752.14100_gp,&
         1598.63960_gp, 1400.50680_gp,  635.91570_gp, 1415.27090_gp, 1358.19040_gp,&
         1323.93290_gp, 1292.44750_gp, 1264.52130_gp,  991.51760_gp, 1121.85490_gp,&
         1081.60010_gp, 1139.68250_gp, 1115.39780_gp, 1093.46020_gp, 1080.79740_gp,&
         908.96190_gp,  891.27890_gp,  812.08250_gp,  684.99030_gp,  696.40000_gp,&
         630.31370_gp,  576.65350_gp,  478.57310_gp,  446.96260_gp,  459.52800_gp,&
         673.46620_gp,  657.02270_gp,  602.34550_gp,  573.74780_gp,  528.72700_gp,&
         483.98270_gp, 1650.13610_gp, 1578.51230_gp, 1391.67280_gp, 1247.07710_gp,&
         1240.04450_gp, 1200.60280_gp, 1238.31870_gp, 1198.77200_gp,   66.40570_gp,&
         216.86950_gp,  276.38280_gp,  210.14480_gp,  157.94060_gp,  109.55230_gp,&
         79.64770_gp,   54.68520_gp,  317.40960_gp,  491.13850_gp,  496.09590_gp,&
         397.17520_gp,  324.42460_gp,  273.88040_gp,  223.50790_gp,  436.88540_gp,&
         841.05190_gp,  435.51720_gp,  420.46710_gp,  412.08740_gp,  408.44600_gp,&
         375.10510_gp,  346.98680_gp,  331.02820_gp,  323.67720_gp,  320.26840_gp,&
         299.70920_gp,  491.91570_gp,  430.81190_gp,  385.30990_gp,  351.36440_gp/)
    vdwparams%coeffs(11301:11400)=(/  309.05750_gp,  520.97030_gp, 1024.01010_gp,  778.56140_gp,  576.24530_gp,&
         581.93050_gp,  542.38030_gp,  609.86550_gp,  471.94410_gp,  441.47690_gp,&
         409.84440_gp,  396.58040_gp,  391.74910_gp,  624.94310_gp,  563.83860_gp,&
         527.09190_gp,  498.01930_gp,  454.22360_gp,  614.77890_gp, 1329.75090_gp,&
         986.04390_gp,  616.65720_gp,  603.88680_gp,  584.67880_gp,  587.78910_gp,&
         39.25640_gp,   25.69130_gp,  653.68460_gp,  365.11120_gp,  241.63430_gp,&
         161.39330_gp,  112.05120_gp,   84.47310_gp,   63.79960_gp,   48.98090_gp,&
         780.26740_gp,  585.23670_gp,  534.37570_gp,  415.91500_gp,  321.21200_gp,&
         265.08690_gp,  215.46030_gp,  175.59080_gp, 1285.69650_gp, 1039.46620_gp,&
         854.82160_gp,  822.18100_gp,  751.08930_gp,  591.66600_gp,  644.44210_gp,&
         505.60640_gp,  533.00400_gp,  550.39520_gp,  422.13760_gp,  429.93040_gp,&
         510.21570_gp,  445.63750_gp,  376.99910_gp,  336.69100_gp,  293.30360_gp,&
         254.08590_gp, 1438.54450_gp, 1239.55040_gp, 1077.08580_gp,  963.27150_gp,&
         875.59310_gp,  672.49140_gp,  751.75670_gp,  569.50600_gp,  622.32140_gp,&
         575.99650_gp,  481.15240_gp,  506.33030_gp,  639.22390_gp,  586.33170_gp,&
         518.59750_gp,  479.06490_gp,  431.22670_gp,  386.14600_gp, 1751.83190_gp,&
         1584.59290_gp, 1380.92390_gp,  613.65370_gp, 1400.38800_gp, 1342.92470_gp,&
         1308.74240_gp, 1277.35680_gp, 1249.50540_gp,  974.06900_gp, 1111.18770_gp,&
         1069.99210_gp, 1124.67710_gp, 1100.55120_gp, 1078.67940_gp, 1066.43870_gp,&
         893.16870_gp,  871.10070_gp,  791.18140_gp,  665.83170_gp,  676.05760_gp/)
    vdwparams%coeffs(11401:11500)=(/  610.29470_gp,  557.18050_gp,  461.53520_gp,  430.73110_gp,  442.32070_gp,&
         655.32760_gp,  636.64300_gp,  580.99830_gp,  551.99550_gp,  507.19420_gp,&
         463.10600_gp, 1643.60580_gp, 1560.42360_gp, 1369.04640_gp, 1220.12770_gp,&
         1216.67380_gp, 1177.86120_gp, 1219.14500_gp, 1179.46900_gp,   63.76290_gp,&
         211.74000_gp,  269.23050_gp,  202.61480_gp,  151.47130_gp,  104.50970_gp,&
         75.69080_gp,   51.75570_gp,  310.68640_gp,  481.12160_gp,  483.24150_gp,&
         383.77640_gp,  311.71580_gp,  262.22190_gp,  213.26140_gp,  425.93700_gp,&
         831.27500_gp,  421.69690_gp,  407.00940_gp,  398.98160_gp,  395.79870_gp,&
         362.04930_gp,  334.56150_gp,  319.19930_gp,  312.22210_gp,  309.89590_gp,&
         288.40040_gp,  478.03960_gp,  416.01900_gp,  370.44180_gp,  336.85390_gp,&
         295.42100_gp,  507.22850_gp, 1012.80780_gp,  762.00330_gp,  558.29510_gp,&
         563.98160_gp,  524.90850_gp,  593.23860_gp,  455.13890_gp,  425.67220_gp,&
         394.96580_gp,  382.42790_gp,  376.61900_gp,  607.61490_gp,  545.18900_gp,&
         507.75680_gp,  478.58010_gp,  435.29630_gp,  595.94990_gp, 1317.60670_gp,&
         965.82370_gp,  594.90700_gp,  582.56410_gp,  563.90150_gp,  567.31710_gp,&
         549.02810_gp,   41.32750_gp,   27.15950_gp,  659.09240_gp,  376.46100_gp,&
         251.84240_gp,  169.31590_gp,  117.99770_gp,   89.13590_gp,   67.40300_gp,&
         51.77560_gp,  787.61870_gp,  601.03840_gp,  552.15450_gp,  433.06620_gp,&
         336.37000_gp,  278.45510_gp,  226.91090_gp,  185.26210_gp, 1292.77330_gp,&
         1059.71920_gp,  874.13140_gp,  842.88400_gp,  771.18400_gp,  607.34870_gp/)
    vdwparams%coeffs(11501:11600)=(/  663.13130_gp,  520.25690_gp,  550.50450_gp,  567.60260_gp,  435.08290_gp,&
         445.49820_gp,  528.13050_gp,  464.06710_gp,  394.50330_gp,  353.25210_gp,&
         308.47630_gp,  267.74640_gp, 1447.93600_gp, 1262.84910_gp, 1103.07420_gp,&
         989.59400_gp,  901.27470_gp,  694.29580_gp,  775.29730_gp,  589.24130_gp,&
         644.13130_gp,  596.86540_gp,  497.84100_gp,  525.41210_gp,  660.98150_gp,&
         609.40140_gp,  541.46000_gp,  501.44410_gp,  452.49160_gp,  406.04790_gp,&
         1763.69170_gp, 1611.37650_gp, 1411.64580_gp,  639.99700_gp, 1426.34050_gp,&
         1368.99450_gp, 1334.48230_gp, 1302.76260_gp, 1274.62830_gp,  999.06170_gp,&
         1130.03100_gp, 1089.26140_gp, 1148.83040_gp, 1124.36570_gp, 1102.26140_gp,&
         1089.52000_gp,  915.93800_gp,  898.12800_gp,  818.11890_gp,  689.77370_gp,&
         701.24130_gp,  634.51700_gp,  580.34520_gp,  481.40120_gp,  449.51810_gp,&
         462.16400_gp,  678.04940_gp,  661.43650_gp,  606.16090_gp,  577.21150_gp,&
         531.69500_gp,  486.48120_gp, 1661.50270_gp, 1591.00200_gp, 1402.61900_gp,&
         1256.53510_gp, 1249.41970_gp, 1209.67210_gp, 1247.87500_gp, 1208.01680_gp,&
         66.82350_gp,  218.61950_gp,  278.50800_gp,  211.54850_gp,  158.83010_gp,&
         110.02650_gp,   79.89490_gp,   54.76000_gp,  319.95720_gp,  495.11470_gp,&
         499.98200_gp,  399.98160_gp,  326.45920_gp,  275.40540_gp,  224.56980_gp,&
         439.97650_gp,  847.44860_gp,  438.55500_gp,  423.30690_gp,  414.87410_gp,&
         411.22300_gp,  377.57470_gp,  349.20950_gp,  333.12400_gp,  325.73280_gp,&
         322.37280_gp,  301.55830_gp,  495.56680_gp,  433.76750_gp,  387.73360_gp/)
    vdwparams%coeffs(11601:11700)=(/  353.39890_gp,  310.65690_gp,  524.53740_gp, 1031.58080_gp,  784.04790_gp,&
         580.16290_gp,  585.94010_gp,  545.92550_gp,  613.98170_gp,  474.85900_gp,&
         444.14170_gp,  412.25320_gp,  398.92850_gp,  394.01610_gp,  629.38940_gp,&
         567.63680_gp,  530.44640_gp,  501.02050_gp,  456.75320_gp,  618.85080_gp,&
         1339.02900_gp,  992.82930_gp,  620.58810_gp,  607.73400_gp,  588.38130_gp,&
         591.53680_gp,  571.04500_gp,  595.36280_gp,   38.91500_gp,   25.66950_gp,&
         612.81660_gp,  352.07530_gp,  236.26590_gp,  159.24220_gp,  111.20430_gp,&
         84.13630_gp,   63.71520_gp,   49.00500_gp,  732.64570_gp,  561.59910_gp,&
         516.72630_gp,  406.15220_gp,  316.05550_gp,  261.97720_gp,  213.76170_gp,&
         174.73800_gp, 1201.33070_gp,  988.50950_gp,  816.03910_gp,  787.43780_gp,&
         720.76790_gp,  567.75360_gp,  620.17110_gp,  486.69710_gp,  515.38390_gp,&
         531.15950_gp,  407.24140_gp,  417.47630_gp,  494.61700_gp,  435.33390_gp,&
         370.64210_gp,  332.21730_gp,  290.41790_gp,  252.33250_gp, 1345.95520_gp,&
         1177.85230_gp, 1030.25770_gp,  925.07300_gp,  843.00140_gp,  650.10950_gp,&
         725.65600_gp,  552.17930_gp,  603.56690_gp,  559.49210_gp,  466.65970_gp,&
         492.76250_gp,  619.13660_gp,  571.57510_gp,  508.51710_gp,  471.31940_gp,&
         425.69420_gp,  382.34290_gp, 1639.39950_gp, 1502.19220_gp, 1317.83870_gp,&
         600.95180_gp, 1330.38100_gp, 1277.22410_gp, 1245.11090_gp, 1215.58530_gp,&
         1189.39910_gp,  933.72280_gp, 1053.65870_gp, 1015.88400_gp, 1072.39780_gp,&
         1049.59990_gp, 1029.01990_gp, 1017.05720_gp,  855.90310_gp,  840.38840_gp/)
    vdwparams%coeffs(11701:11800)=(/  766.18830_gp,  646.47490_gp,  657.42490_gp,  595.32710_gp,  544.84580_gp,&
         452.30750_gp,  422.48760_gp,  434.46390_gp,  635.42770_gp,  620.48770_gp,&
         569.29560_gp,  542.47560_gp,  500.12070_gp,  457.95880_gp, 1546.30420_gp,&
         1484.29320_gp, 1310.25950_gp, 1175.55200_gp, 1168.10930_gp, 1131.01270_gp,&
         1165.72890_gp, 1128.68500_gp,   62.79710_gp,  204.55140_gp,  260.72450_gp,&
         198.60250_gp,  149.38110_gp,  103.70340_gp,   75.44650_gp,   51.84140_gp,&
         299.24680_gp,  462.88270_gp,  468.08680_gp,  375.24590_gp,  306.76780_gp,&
         259.10430_gp,  211.55470_gp,  412.06550_gp,  790.52270_gp,  411.39690_gp,&
         397.11900_gp,  389.20800_gp,  385.71070_gp,  354.49400_gp,  327.97050_gp,&
         312.87350_gp,  305.91120_gp,  302.54470_gp,  283.37650_gp,  464.32020_gp,&
         407.07180_gp,  364.31690_gp,  332.35190_gp,  292.45640_gp,  491.46210_gp,&
         962.09720_gp,  733.26220_gp,  544.16190_gp,  549.59500_gp,  512.26070_gp,&
         575.35150_gp,  446.04100_gp,  417.25200_gp,  387.39090_gp,  374.83220_gp,&
         370.45920_gp,  589.83000_gp,  532.68550_gp,  498.26300_gp,  470.94350_gp,&
         429.69450_gp,  580.38410_gp, 1248.01820_gp,  928.27050_gp,  582.77730_gp,&
         570.71750_gp,  552.58430_gp,  555.42950_gp,  535.87600_gp,  558.97620_gp,&
         524.91950_gp,   39.39670_gp,   25.86580_gp,  638.48500_gp,  361.18220_gp,&
         240.78360_gp,  161.56910_gp,  112.48360_gp,   84.92880_gp,   64.20560_gp,&
         49.31660_gp,  762.54390_gp,  577.48050_gp,  529.42090_gp,  414.17170_gp,&
         321.13630_gp,  265.61240_gp,  216.29400_gp,  176.51170_gp, 1255.34340_gp/)
    vdwparams%coeffs(11801:11900)=(/ 1021.33420_gp,  841.37840_gp,  810.57120_gp,  741.19950_gp,  583.88160_gp,&
         636.84280_gp,  499.70290_gp,  527.94350_gp,  544.62440_gp,  417.63490_gp,&
         426.75460_gp,  506.13150_gp,  443.82260_gp,  376.72490_gp,  337.07580_gp,&
         294.15100_gp,  255.18090_gp, 1405.60430_gp, 1217.60720_gp, 1061.30520_gp,&
         950.99220_gp,  865.53240_gp,  666.08950_gp,  744.08720_gp,  564.89840_gp,&
         617.35350_gp,  571.82020_gp,  477.29120_gp,  503.12050_gp,  633.74930_gp,&
         583.21300_gp,  517.42880_gp,  478.82440_gp,  431.76670_gp,  387.21870_gp,&
         1712.88300_gp, 1555.21200_gp, 1359.38020_gp,  611.83610_gp, 1375.84720_gp,&
         1319.77470_gp, 1286.34330_gp, 1255.63290_gp, 1228.38970_gp,  960.89950_gp,&
         1091.37760_gp, 1051.86090_gp, 1106.47400_gp, 1082.82590_gp, 1061.43290_gp,&
         1049.24400_gp,  881.06130_gp,  861.94390_gp,  784.33470_gp,  660.95780_gp,&
         671.62740_gp,  607.23190_gp,  555.05840_gp,  460.25170_gp,  429.70270_gp,&
         441.59570_gp,  650.09130_gp,  633.14710_gp,  579.41850_gp,  551.35470_gp,&
         507.49380_gp,  464.06150_gp, 1609.97090_gp, 1533.83910_gp, 1349.48460_gp,&
         1206.51020_gp, 1201.10490_gp, 1162.82610_gp, 1201.03270_gp, 1162.36430_gp,&
         63.78160_gp,  209.66850_gp,  266.99220_gp,  202.15950_gp,  151.59630_gp,&
         104.89910_gp,   76.12090_gp,   52.14540_gp,  307.14640_gp,  475.48270_gp,&
         479.19230_gp,  382.42860_gp,  311.66960_gp,  262.71630_gp,  214.07020_gp,&
         422.05350_gp,  817.42610_gp,  419.57740_gp,  405.07580_gp,  397.01800_gp,&
         393.63860_gp,  360.94370_gp,  333.74440_gp,  318.40100_gp,  311.36890_gp/)
    vdwparams%coeffs(11901:12000)=(/  308.44670_gp,  288.04090_gp,  474.67810_gp,  414.68730_gp,  370.23680_gp,&
         337.22390_gp,  296.24340_gp,  503.00950_gp,  995.63320_gp,  753.68720_gp,&
         555.31660_gp,  560.79120_gp,  522.44060_gp,  588.68540_gp,  453.92960_gp,&
         424.57780_gp,  394.05630_gp,  381.38420_gp,  376.31210_gp,  603.08050_gp,&
         542.94040_gp,  506.81590_gp,  478.39610_gp,  435.83490_gp,  592.59730_gp,&
         1294.35010_gp,  954.98290_gp,  593.24060_gp,  580.94220_gp,  562.41030_gp,&
         565.58410_gp,  546.41960_gp,  569.23180_gp,  534.34800_gp,  544.43560_gp,&
         37.02950_gp,   24.48060_gp,  575.74990_gp,  332.86460_gp,  224.08560_gp,&
         151.36070_gp,  105.85820_gp,   80.17040_gp,   60.75980_gp,   46.75950_gp,&
         688.59100_gp,  530.38700_gp,  488.84720_gp,  385.09520_gp,  300.19740_gp,&
         249.09840_gp,  203.45290_gp,  166.44660_gp, 1127.94800_gp,  931.74430_gp,&
         769.83360_gp,  743.40690_gp,  680.76970_gp,  536.27310_gp,  586.13150_gp,&
         460.04150_gp,  487.62020_gp,  502.32060_gp,  385.13270_gp,  395.36540_gp,&
         468.22770_gp,  412.81620_gp,  351.98750_gp,  315.77010_gp,  276.27680_gp,&
         240.22750_gp, 1264.15180_gp, 1110.04420_gp,  972.38640_gp,  873.89940_gp,&
         796.83270_gp,  615.10520_gp,  686.33810_gp,  522.81640_gp,  571.48280_gp,&
         529.93880_gp,  441.89760_gp,  466.94330_gp,  586.03950_gp,  541.79100_gp,&
         482.65820_gp,  447.69600_gp,  404.68130_gp,  363.73600_gp, 1539.84530_gp,&
         1414.99430_gp, 1243.18830_gp,  570.24510_gp, 1253.76910_gp, 1203.97270_gp,&
         1173.78420_gp, 1146.01940_gp, 1121.39780_gp,  881.75490_gp,  992.62670_gp/)
    vdwparams%coeffs(12001:12100)=(/  957.31530_gp, 1011.46770_gp,  990.00700_gp,  970.65280_gp,  959.30580_gp,&
         808.19030_gp,  794.71630_gp,  725.17730_gp,  612.26400_gp,  622.84620_gp,&
         564.42200_gp,  516.85750_gp,  429.31700_gp,  401.10390_gp,  412.59070_gp,&
         601.65160_gp,  588.16150_gp,  540.28410_gp,  515.17620_gp,  475.32260_gp,&
         435.54890_gp, 1454.15440_gp, 1399.20620_gp, 1236.84060_gp, 1111.36080_gp,&
         1103.48210_gp, 1068.47830_gp, 1100.23640_gp, 1065.46250_gp,   59.65860_gp,&
         193.47380_gp,  246.74640_gp,  188.47650_gp,  141.97970_gp,   98.72120_gp,&
         71.91040_gp,   49.47930_gp,  282.86510_gp,  437.42550_gp,  443.00780_gp,&
         355.89820_gp,  291.39070_gp,  246.35930_gp,  201.34780_gp,  389.95180_gp,&
         745.17810_gp,  390.02770_gp,  376.50890_gp,  368.99410_gp,  365.59720_gp,&
         336.35820_gp,  311.27980_gp,  296.94730_gp,  290.31390_gp,  286.89270_gp,&
         269.09650_gp,  439.73990_gp,  386.16340_gp,  346.00850_gp,  315.89290_gp,&
         278.20320_gp,  465.26060_gp,  906.71060_gp,  693.06070_gp,  515.80860_gp,&
         520.93800_gp,  485.72140_gp,  544.79060_gp,  423.34180_gp,  396.04750_gp,&
         367.76460_gp,  355.78930_gp,  351.90760_gp,  558.57850_gp,  505.19270_gp,&
         473.00760_gp,  447.36060_gp,  408.47930_gp,  550.05480_gp, 1175.43320_gp,&
         877.14780_gp,  553.04140_gp,  541.60480_gp,  524.43100_gp,  527.02430_gp,&
         508.12970_gp,  530.36160_gp,  498.13190_gp,  506.88950_gp,  472.79500_gp,&
         36.00840_gp,   23.82890_gp,  557.24540_gp,  322.88400_gp,  217.62790_gp,&
         147.12440_gp,  102.95880_gp,   78.00790_gp,   59.14230_gp,   45.52790_gp/)
    vdwparams%coeffs(12101:12200)=(/  666.55340_gp,  514.28720_gp,  474.30950_gp,  373.95400_gp,  291.71000_gp,&
         242.15880_gp,  197.86420_gp,  161.93000_gp, 1091.50910_gp,  902.84450_gp,&
         746.18100_gp,  720.76670_gp,  660.14630_gp,  520.04920_gp,  568.51110_gp,&
         446.24320_gp,  473.14860_gp,  487.33120_gp,  373.65280_gp,  383.76960_gp,&
         454.41890_gp,  400.89710_gp,  342.01700_gp,  306.93000_gp,  268.63470_gp,&
         233.65510_gp, 1223.46680_gp, 1075.56590_gp,  942.68540_gp,  847.48690_gp,&
         772.91870_gp,  596.87070_gp,  665.89970_gp,  507.45680_gp,  554.68820_gp,&
         514.43460_gp,  428.94240_gp,  453.36080_gp,  568.75080_gp,  526.08060_gp,&
         468.89680_gp,  435.06030_gp,  393.38290_gp,  353.68410_gp, 1490.33290_gp,&
         1370.81200_gp, 1205.00730_gp,  553.93560_gp, 1214.84100_gp, 1166.69030_gp,&
         1137.46450_gp, 1110.58190_gp, 1086.74370_gp,  855.01380_gp,  961.71480_gp,&
         927.60080_gp,  980.34030_gp,  959.55360_gp,  940.81400_gp,  929.79290_gp,&
         783.64610_gp,  770.98740_gp,  703.75210_gp,  594.32800_gp,  604.67390_gp,&
         548.10470_gp,  502.02620_gp,  417.09650_gp,  389.72370_gp,  400.92380_gp,&
         583.98610_gp,  571.11950_gp,  524.86410_gp,  500.60050_gp,  462.01380_gp,&
         423.46810_gp, 1407.99440_gp, 1355.89080_gp, 1199.13830_gp, 1078.08360_gp,&
         1070.15400_gp, 1036.22220_gp, 1066.65600_gp, 1033.00900_gp,   57.97600_gp,&
         187.70500_gp,  239.44310_gp,  183.08850_gp,  138.00410_gp,   96.01930_gp,&
         69.97910_gp,   48.18050_gp,  274.37320_gp,  424.24930_gp,  429.89820_gp,&
         345.64250_gp,  283.15870_gp,  239.49360_gp,  195.81550_gp,  378.42050_gp/)
    vdwparams%coeffs(12201:12300)=(/  722.11110_gp,  378.74080_gp,  365.62310_gp,  358.32130_gp,  354.99450_gp,&
         326.72790_gp,  302.40170_gp,  288.47830_gp,  282.02540_gp,  278.62110_gp,&
         261.47450_gp,  426.84110_gp,  375.06910_gp,  336.21830_gp,  307.04780_gp,&
         270.50300_gp,  451.56840_gp,  878.58370_gp,  672.26710_gp,  500.85680_gp,&
         505.82970_gp,  471.70110_gp,  528.80060_gp,  411.27310_gp,  384.77140_gp,&
         357.31920_gp,  345.66580_gp,  341.98940_gp,  542.19990_gp,  490.64160_gp,&
         459.55210_gp,  434.74240_gp,  397.07280_gp,  534.08180_gp, 1138.73290_gp,&
         850.75470_gp,  537.23990_gp,  526.13310_gp,  509.46310_gp,  511.94390_gp,&
         493.46840_gp,  515.17330_gp,  483.89820_gp,  492.33800_gp,  459.31550_gp,&
         446.23140_gp,   35.74800_gp,   23.62970_gp,  555.25980_gp,  321.23930_gp,&
         216.30990_gp,  146.11730_gp,  102.18880_gp,   77.38740_gp,   58.64650_gp,&
         45.12980_gp,  664.09940_gp,  511.80090_gp,  471.79450_gp,  371.72710_gp,&
         289.80090_gp,  240.47400_gp,  196.40620_gp,  160.67610_gp, 1087.69700_gp,&
         898.86340_gp,  742.73900_gp,  717.28970_gp,  656.88060_gp,  517.43220_gp,&
         565.59590_gp,  443.90510_gp,  470.58200_gp,  484.75110_gp,  371.63760_gp,&
         381.58090_gp,  451.90020_gp,  398.47870_gp,  339.79060_gp,  304.83420_gp,&
         266.70930_gp,  231.90470_gp, 1219.07060_gp, 1070.83510_gp,  938.18730_gp,&
         843.23610_gp,  768.91010_gp,  593.57780_gp,  662.31090_gp,  504.53730_gp,&
         551.52390_gp,  511.44340_gp,  426.43750_gp,  450.65890_gp,  565.56470_gp,&
         522.93560_gp,  465.90590_gp,  432.17320_gp,  390.65890_gp,  351.13440_gp/)
    vdwparams%coeffs(12301:12400)=(/ 1484.96610_gp, 1364.92400_gp, 1199.38740_gp,  550.42640_gp, 1209.46630_gp,&
         1161.44940_gp, 1132.33650_gp, 1105.56040_gp, 1081.81600_gp,  850.75460_gp,&
         957.49480_gp,  923.47850_gp,  975.80850_gp,  955.11040_gp,  936.44570_gp,&
         925.49440_gp,  779.79590_gp,  766.90850_gp,  699.85110_gp,  590.88900_gp,&
         601.12600_gp,  544.76460_gp,  498.87240_gp,  414.37420_gp,  387.14220_gp,&
         398.24650_gp,  580.61670_gp,  567.66880_gp,  521.51170_gp,  497.29730_gp,&
         458.84250_gp,  420.45390_gp, 1402.46130_gp, 1349.79310_gp, 1193.34120_gp,&
         1072.41820_gp, 1064.71990_gp, 1030.94560_gp, 1061.48600_gp, 1027.95650_gp,&
         57.59370_gp,  186.72100_gp,  238.14200_gp,  181.93830_gp,  137.05810_gp,&
         95.29790_gp,   69.41370_gp,   47.75510_gp,  272.96680_gp,  422.11770_gp,&
         427.56550_gp,  343.54750_gp,  281.29780_gp,  237.82830_gp,  194.37320_gp,&
         376.30470_gp,  718.87190_gp,  376.45520_gp,  363.40760_gp,  356.15030_gp,&
         352.86320_gp,  324.67610_gp,  300.47150_gp,  286.63320_gp,  280.22740_gp,&
         276.90380_gp,  259.76130_gp,  424.42030_gp,  372.75950_gp,  334.01900_gp,&
         304.95200_gp,  268.56840_gp,  448.98400_gp,  874.67070_gp,  668.76560_gp,&
         497.83950_gp,  502.78080_gp,  468.80410_gp,  525.74810_gp,  408.62370_gp,&
         382.27390_gp,  354.97280_gp,  343.40860_gp,  339.69190_gp,  539.08680_gp,&
         487.62850_gp,  456.59480_gp,  431.85070_gp,  394.32570_gp,  530.87060_gp,&
         1133.83290_gp,  846.37940_gp,  533.82250_gp,  522.78360_gp,  506.20810_gp,&
         508.70480_gp,  490.43340_gp,  511.92810_gp,  480.82180_gp,  489.26080_gp/)
    vdwparams%coeffs(12401:12500)=(/  456.37280_gp,  443.36390_gp,  440.52280_gp,   37.39860_gp,   24.54350_gp,&
         592.36510_gp,  340.23850_gp,  227.88020_gp,  153.21460_gp,  106.73400_gp,&
         80.58630_gp,   60.90310_gp,   46.75690_gp,  708.06290_gp,  542.77630_gp,&
         499.10710_gp,  391.84510_gp,  304.44600_gp,  251.99750_gp,  205.29690_gp,&
         167.55260_gp, 1160.07670_gp,  955.07870_gp,  788.43250_gp,  760.56390_gp,&
         696.06330_gp,  547.95880_gp,  598.76220_gp,  469.57880_gp,  497.40630_gp,&
         512.74730_gp,  392.80830_gp,  402.72690_gp,  477.37550_gp,  419.82970_gp,&
         357.02170_gp,  319.68520_gp,  279.12740_gp,  242.21470_gp, 1299.44820_gp,&
         1137.75250_gp,  994.98690_gp,  893.15950_gp,  813.67980_gp,  626.96020_gp,&
         700.06760_gp,  532.19550_gp,  581.98660_gp,  539.36520_gp,  449.56060_gp,&
         474.87730_gp,  597.14660_gp,  551.05660_gp,  489.86860_gp,  453.72430_gp,&
         409.44070_gp,  367.38410_gp, 1582.55000_gp, 1450.77940_gp, 1272.58220_gp,&
         578.83400_gp, 1284.66820_gp, 1233.33650_gp, 1202.33710_gp, 1173.83920_gp,&
         1148.56550_gp,  901.12100_gp, 1016.90490_gp,  980.43260_gp, 1035.60000_gp,&
         1013.60380_gp,  993.74150_gp,  982.22660_gp,  826.28510_gp,  811.15710_gp,&
         739.22210_gp,  623.22500_gp,  633.76790_gp,  573.62030_gp,  524.73400_gp,&
         435.20090_gp,  406.35950_gp,  417.92020_gp,  612.37220_gp,  597.91090_gp,&
         548.26140_gp,  522.17940_gp,  481.05700_gp,  440.14200_gp, 1492.46700_gp,&
         1433.26400_gp, 1265.08090_gp, 1134.41750_gp, 1127.20490_gp, 1091.36560_gp,&
         1125.10940_gp, 1089.33700_gp,   60.48710_gp,  197.59460_gp,  251.72230_gp/)
    vdwparams%coeffs(12501:12600)=(/  191.41090_gp,  143.69300_gp,   99.50590_gp,   72.21820_gp,   49.45020_gp,&
         289.01230_gp,  447.16520_gp,  452.00550_gp,  361.91400_gp,  295.46110_gp,&
         249.22870_gp,  203.17400_gp,  397.28560_gp,  763.26450_gp,  396.58450_gp,&
         382.75050_gp,  375.10830_gp,  371.75280_gp,  341.56320_gp,  315.90920_gp,&
         301.32730_gp,  294.62620_gp,  291.46390_gp,  272.84960_gp,  448.03280_gp,&
         392.45010_gp,  350.89060_gp,  319.81290_gp,  281.09900_gp,  473.64740_gp,&
         928.71230_gp,  707.39100_gp,  524.42450_gp,  529.63760_gp,  493.48140_gp,&
         554.47820_gp,  429.41670_gp,  401.59420_gp,  372.73910_gp,  360.66920_gp,&
         356.42060_gp,  568.76740_gp,  513.38040_gp,  479.92240_gp,  453.34590_gp,&
         413.30020_gp,  559.18680_gp, 1204.51190_gp,  895.47410_gp,  561.28730_gp,&
         549.66430_gp,  532.16530_gp,  534.96580_gp,  516.24120_gp,  538.45600_gp,&
         505.56900_gp,  514.73030_gp,  479.73200_gp,  466.00690_gp,  463.06960_gp,&
         487.07530_gp,   34.60120_gp,   22.89260_gp,  533.45810_gp,  309.82450_gp,&
         209.00590_gp,  141.34580_gp,   98.92230_gp,   74.94510_gp,   56.81240_gp,&
         43.72630_gp,  638.16340_gp,  493.28830_gp,  455.19600_gp,  359.11410_gp,&
         280.24060_gp,  232.66870_gp,  190.12190_gp,  155.59160_gp, 1044.47010_gp,&
         865.28290_gp,  715.36520_gp,  691.15730_gp,  633.11630_gp,  498.70530_gp,&
         545.33910_gp,  428.01770_gp,  454.01590_gp,  467.56420_gp,  358.43930_gp,&
         368.35200_gp,  436.14550_gp,  384.97400_gp,  328.54540_gp,  294.88120_gp,&
         258.11290_gp,  224.51100_gp, 1170.84000_gp, 1030.71650_gp,  903.85210_gp/)
    vdwparams%coeffs(12601:12700)=(/  812.81250_gp,  741.41830_gp,  572.66250_gp,  638.84880_gp,  486.94650_gp,&
         532.31850_gp,  493.73180_gp,  411.58220_gp,  435.15940_gp,  545.77490_gp,&
         505.06710_gp,  450.32940_gp,  417.90250_gp,  377.92070_gp,  339.81250_gp,&
         1426.23700_gp, 1313.36180_gp, 1155.12990_gp,  531.93250_gp, 1164.09970_gp,&
         1118.06470_gp, 1090.08720_gp, 1064.35050_gp, 1041.52960_gp,  819.84680_gp,&
         921.33410_gp,  888.74530_gp,  939.68880_gp,  919.78180_gp,  901.84110_gp,&
         891.26080_gp,  751.43080_gp,  739.69160_gp,  675.35320_gp,  570.39690_gp,&
         580.39830_gp,  526.19340_gp,  482.01740_gp,  400.48290_gp,  374.20390_gp,&
         385.00720_gp,  560.37600_gp,  548.25190_gp,  504.02540_gp,  480.80590_gp,&
         443.81570_gp,  406.83040_gp, 1348.01980_gp, 1299.40570_gp, 1149.74840_gp,&
         1034.16690_gp, 1026.25430_gp,  993.71920_gp, 1022.57380_gp,  990.38000_gp,&
         55.70080_gp,  180.12920_gp,  229.80950_gp,  175.85020_gp,  132.57550_gp,&
         92.25120_gp,   67.23080_gp,   46.27720_gp,  263.22190_gp,  406.99260_gp,&
         412.61700_gp,  331.94510_gp,  272.02390_gp,  230.10400_gp,  188.15090_gp,&
         363.07990_gp,  692.00830_gp,  363.63190_gp,  351.03090_gp,  344.01120_gp,&
         340.78980_gp,  313.76160_gp,  290.41450_gp,  277.03380_gp,  270.82780_gp,&
         267.48920_gp,  251.14070_gp,  409.72790_gp,  360.20340_gp,  322.97840_gp,&
         294.99200_gp,  259.90550_gp,  433.29230_gp,  841.85430_gp,  644.80130_gp,&
         480.82360_gp,  485.58420_gp,  452.84900_gp,  507.43830_gp,  394.93050_gp,&
         369.47180_gp,  343.11040_gp,  331.90160_gp,  328.46280_gp,  520.37630_gp/)
    vdwparams%coeffs(12701:12800)=(/  471.10690_gp,  441.37250_gp,  417.60220_gp,  381.46630_gp,  512.66580_gp,&
         1090.84220_gp,  815.91040_gp,  515.90750_gp,  505.24280_gp,  489.24070_gp,&
         491.59590_gp,  473.74820_gp,  494.70080_gp,  464.68280_gp,  472.73450_gp,&
         441.09990_gp,  428.54250_gp,  425.78520_gp,  447.51580_gp,  411.56650_gp,&
         45.34660_gp,   28.99860_gp,  833.44750_gp,  442.56360_gp,  286.28310_gp,&
         187.98820_gp,  128.79820_gp,   96.14240_gp,   71.95970_gp,   54.81670_gp,&
         991.34020_gp,  714.75720_gp,  644.92820_gp,  493.84800_gp,  376.40650_gp,&
         307.98130_gp,  248.21220_gp,  200.73430_gp, 1653.94760_gp, 1289.46500_gp,&
         1053.07890_gp, 1007.37780_gp,  917.14520_gp,  722.38490_gp,  783.08020_gp,&
         613.78620_gp,  642.19290_gp,  665.32650_gp,  510.28620_gp,  514.19020_gp,&
         612.87020_gp,  528.47210_gp,  442.18700_gp,  392.24560_gp,  339.28530_gp,&
         291.96450_gp, 1846.96990_gp, 1540.34240_gp, 1323.03110_gp, 1175.03450_gp,&
         1063.47560_gp,  810.65520_gp,  908.84210_gp,  682.68070_gp,  745.65500_gp,&
         688.14480_gp,  575.95500_gp,  602.64510_gp,  768.01030_gp,  696.81050_gp,&
         610.25270_gp,  560.48120_gp,  501.36210_gp,  446.26410_gp, 2252.85860_gp,&
         1978.99090_gp, 1704.02700_gp,  723.44560_gp, 1742.07290_gp, 1666.16150_gp,&
         1622.69780_gp, 1582.90710_gp, 1547.57500_gp, 1192.02260_gp, 1389.46840_gp,&
         1336.06130_gp, 1388.41730_gp, 1358.10650_gp, 1330.44730_gp, 1315.98590_gp,&
         1093.93510_gp, 1054.31080_gp,  951.18510_gp,  796.76620_gp,  806.79990_gp,&
         724.12960_gp,  658.06160_gp,  542.45470_gp,  505.20370_gp,  517.67100_gp/)
    vdwparams%coeffs(12801:12900)=(/  785.62040_gp,  756.33770_gp,  683.97230_gp,  646.57460_gp,  590.53860_gp,&
         536.25670_gp, 2090.85420_gp, 1937.01660_gp, 1680.42010_gp, 1479.79560_gp,&
         1484.71460_gp, 1436.68170_gp, 1497.33020_gp, 1446.49010_gp,   74.62920_gp,&
         255.92080_gp,  324.43230_gp,  238.97670_gp,  176.54080_gp,  120.10680_gp,&
         85.93700_gp,   57.82720_gp,  376.99540_gp,  585.70530_gp,  581.57180_gp,&
         454.71270_gp,  365.12350_gp,  304.70560_gp,  245.70480_gp,  512.81530_gp,&
         1034.14560_gp,  500.45280_gp,  483.16910_gp,  473.64140_gp,  470.60760_gp,&
         426.99780_gp,  393.70530_gp,  375.65220_gp,  367.64140_gp,  366.93480_gp,&
         337.96290_gp,  572.32270_gp,  491.97380_gp,  434.22870_gp,  392.47750_gp,&
         341.85570_gp,  609.10480_gp, 1263.24590_gp,  929.11650_gp,  664.00370_gp,&
         670.40190_gp,  622.62030_gp,  711.74240_gp,  535.54020_gp,  500.50420_gp,&
         463.68240_gp,  449.30820_gp,  439.97740_gp,  727.21640_gp,  645.40420_gp,&
         596.73110_gp,  559.72270_gp,  506.15680_gp,  710.03790_gp, 1655.12170_gp,&
         1181.11300_gp,  700.92550_gp,  686.27630_gp,  663.95860_gp,  669.17280_gp,&
         650.59990_gp,  673.80220_gp,  631.35970_gp,  645.87040_gp,  597.87350_gp,&
         580.33350_gp,  577.00390_gp,  608.71420_gp,  556.95720_gp,  780.38720_gp,&
         42.59910_gp,   27.62080_gp,  738.80890_gp,  403.22390_gp,  264.58080_gp,&
         175.62680_gp,  121.32320_gp,   91.10250_gp,   68.54250_gp,   52.43580_gp,&
         880.35320_gp,  648.29810_gp,  589.23010_gp,  455.75650_gp,  350.31540_gp,&
         288.21970_gp,  233.53600_gp,  189.76410_gp, 1460.47070_gp, 1159.75040_gp/)
    vdwparams%coeffs(12901:13000)=(/  950.68620_gp,  912.32910_gp,  832.23310_gp,  655.70910_gp,  712.59280_gp,&
         558.95410_gp,  587.24740_gp,  607.19280_gp,  465.80690_gp,  472.25470_gp,&
         561.60790_gp,  488.07900_gp,  411.25650_gp,  366.40670_gp,  318.37520_gp,&
         275.12740_gp, 1632.85260_gp, 1384.34960_gp, 1196.70640_gp, 1067.06980_gp,&
         968.23440_gp,  741.42020_gp,  829.76170_gp,  626.45270_gp,  684.26530_gp,&
         632.55740_gp,  528.92880_gp,  555.16740_gp,  703.65600_gp,  642.57250_gp,&
         566.23950_gp,  521.99830_gp,  468.82570_gp,  418.91390_gp, 1990.65140_gp,&
         1774.24330_gp, 1537.73740_gp,  670.54450_gp, 1565.02270_gp, 1498.89090_gp,&
         1460.27400_gp, 1424.86680_gp, 1393.43930_gp, 1080.73610_gp, 1245.46540_gp,&
         1198.55620_gp, 1252.27980_gp, 1225.17800_gp, 1200.54280_gp, 1187.15470_gp,&
         991.15750_gp,  961.78800_gp,  871.15170_gp,  731.91440_gp,  742.27510_gp,&
         668.51880_gp,  609.21920_gp,  503.70150_gp,  469.68660_gp,  481.88500_gp,&
         720.91450_gp,  697.65140_gp,  634.43450_gp,  601.66450_gp,  551.63790_gp,&
         502.69970_gp, 1858.05080_gp, 1742.45060_gp, 1520.81030_gp, 1348.41470_gp,&
         1348.25880_gp, 1304.94040_gp, 1354.56150_gp, 1309.60990_gp,   69.51720_gp,&
         233.60750_gp,  296.81730_gp,  221.50850_gp,  164.89220_gp,  113.16000_gp,&
         81.56100_gp,   55.39050_gp,  343.28010_gp,  532.47660_gp,  532.28640_gp,&
         420.23260_gp,  339.91220_gp,  285.11970_gp,  231.15710_gp,  469.41050_gp,&
         929.43790_gp,  461.93760_gp,  445.95790_gp,  437.13970_gp,  433.91220_gp,&
         395.61160_gp,  365.27420_gp,  348.51570_gp,  340.95800_gp,  339.10500_gp/)
    vdwparams%coeffs(13001:13100)=(/  314.36170_gp,  525.50400_gp,  455.18330_gp,  404.01890_gp,  366.59830_gp,&
         320.71300_gp,  558.49190_gp, 1133.99630_gp,  844.72340_gp,  612.31920_gp,&
         618.33280_gp,  575.04470_gp,  653.04060_gp,  496.97040_gp,  464.66380_gp,&
         430.87120_gp,  417.26470_gp,  410.02790_gp,  667.80870_gp,  596.61500_gp,&
         554.14990_gp,  521.41070_gp,  473.27510_gp,  654.24880_gp, 1480.60880_gp,&
         1072.21420_gp,  649.90970_gp,  636.38010_gp,  615.87840_gp,  620.06220_gp,&
         601.08450_gp,  624.20110_gp,  585.39210_gp,  597.71440_gp,  554.79820_gp,&
         538.68740_gp,  535.45990_gp,  564.10720_gp,  517.09820_gp,  716.09920_gp,&
         659.75010_gp,   39.35300_gp,   25.89160_gp,  642.22390_gp,  360.75580_gp,&
         240.29820_gp,  161.32870_gp,  112.41380_gp,   84.94280_gp,   64.26040_gp,&
         49.38300_gp,  766.74770_gp,  577.30020_gp,  528.75590_gp,  413.30440_gp,&
         320.49320_gp,  265.21130_gp,  216.10340_gp,  176.47630_gp, 1265.21050_gp,&
         1023.75130_gp,  842.49870_gp,  811.28040_gp,  741.59550_gp,  584.57670_gp,&
         636.89720_gp,  500.04450_gp,  527.55830_gp,  544.33150_gp,  417.76790_gp,&
         426.20470_gp,  505.60850_gp,  442.99990_gp,  375.99840_gp,  336.52400_gp,&
         293.79580_gp,  255.00680_gp, 1416.52560_gp, 1221.14340_gp, 1062.81000_gp,&
         951.67080_gp,  865.88950_gp,  666.29380_gp,  744.29510_gp,  564.98080_gp,&
         617.07470_gp,  571.45890_gp,  477.46470_gp,  502.69480_gp,  633.50610_gp,&
         582.38420_gp,  516.49490_gp,  477.98250_gp,  431.09430_gp,  386.74400_gp,&
         1726.47330_gp, 1561.21120_gp, 1362.42760_gp,  610.99070_gp, 1380.27680_gp/)
    vdwparams%coeffs(13101:13200)=(/ 1323.69430_gp, 1290.02980_gp, 1259.11350_gp, 1231.68390_gp,  962.39750_gp,&
         1096.23410_gp, 1056.07480_gp, 1108.88320_gp, 1085.09740_gp, 1063.56750_gp,&
         1051.38220_gp,  882.06880_gp,  861.87340_gp,  783.90750_gp,  660.76360_gp,&
         671.17640_gp,  606.66990_gp,  554.47220_gp,  459.92100_gp,  429.42760_gp,&
         441.12500_gp,  650.21480_gp,  632.55270_gp,  578.55330_gp,  550.48090_gp,&
         506.71990_gp,  463.45510_gp, 1620.92870_gp, 1538.70120_gp, 1351.57160_gp,&
         1206.99930_gp, 1202.63460_gp, 1164.27370_gp, 1203.40960_gp, 1164.43930_gp,&
         63.65500_gp,  209.43010_gp,  266.76230_gp,  201.80580_gp,  151.43990_gp,&
         104.87960_gp,   76.17470_gp,   52.23350_gp,  306.98940_gp,  475.38840_gp,&
         478.55460_gp,  381.66000_gp,  311.07370_gp,  262.32930_gp,  213.88550_gp,&
         422.24100_gp,  820.19820_gp,  419.03830_gp,  404.59180_gp,  396.56360_gp,&
         393.24210_gp,  360.30430_gp,  333.16370_gp,  317.88470_gp,  310.87480_gp,&
         308.07170_gp,  287.48760_gp,  474.09680_gp,  413.91720_gp,  369.53950_gp,&
         336.67770_gp,  295.88320_gp,  503.28860_gp,  999.60180_gp,  754.69660_gp,&
         554.97370_gp,  560.47340_gp,  522.08970_gp,  588.93300_gp,  453.42630_gp,&
         424.17090_gp,  393.71240_gp,  381.04950_gp,  375.75480_gp,  602.66490_gp,&
         542.10550_gp,  505.90030_gp,  477.55630_gp,  435.15150_gp,  592.56180_gp,&
         1300.77530_gp,  956.59970_gp,  592.41900_gp,  580.13530_gp,  561.63190_gp,&
         564.85240_gp,  545.91020_gp,  568.45330_gp,  533.59460_gp,  543.79510_gp,&
         506.13240_gp,  491.59050_gp,  488.51540_gp,  513.89840_gp,  471.98590_gp/)
    vdwparams%coeffs(13201:13300)=(/  645.82280_gp,  597.47460_gp,  543.40360_gp,   39.39020_gp,   25.66140_gp,&
         677.90210_gp,  371.13480_gp,  243.94890_gp,  162.25790_gp,  112.32650_gp,&
         84.50940_gp,   63.70990_gp,   48.83290_gp,  808.10000_gp,  596.53470_gp,&
         542.54110_gp,  420.13030_gp,  323.33970_gp,  266.32530_gp,  216.07690_gp,&
         175.82140_gp, 1339.38070_gp, 1066.40780_gp,  874.54390_gp,  839.59150_gp,&
         766.07010_gp,  603.76280_gp,  656.18380_gp,  514.91990_gp,  541.09810_gp,&
         559.33980_gp,  429.28690_gp,  435.39440_gp,  517.43050_gp,  450.07850_gp,&
         379.61010_gp,  338.48080_gp,  294.39730_gp,  254.67830_gp, 1497.71740_gp,&
         1272.86620_gp, 1101.14360_gp,  982.33900_gp,  891.65750_gp,  683.33990_gp,&
         764.50150_gp,  577.73560_gp,  630.93880_gp,  583.42570_gp,  488.00010_gp,&
         512.25330_gp,  648.62950_gp,  592.71140_gp,  522.69700_gp,  482.12230_gp,&
         433.31830_gp,  387.49410_gp, 1825.36180_gp, 1630.82130_gp, 1414.54000_gp,&
         618.99030_gp, 1438.99860_gp, 1378.51660_gp, 1343.06280_gp, 1310.54570_gp,&
         1281.68250_gp,  994.95010_gp, 1144.82820_gp, 1101.70610_gp, 1152.09550_gp,&
         1127.18460_gp, 1104.54930_gp, 1092.18650_gp,  912.29470_gp,  885.87570_gp,&
         802.82330_gp,  674.90990_gp,  684.56550_gp,  616.87960_gp,  562.43030_gp,&
         465.38980_gp,  434.11300_gp,  445.39470_gp,  664.87270_gp,  643.74060_gp,&
         585.78780_gp,  555.76230_gp,  509.86150_gp,  464.93840_gp, 1705.33180_gp,&
         1602.28660_gp, 1399.53460_gp, 1242.02050_gp, 1241.49570_gp, 1201.69270_gp,&
         1246.89490_gp, 1205.63530_gp,   64.16430_gp,  215.07820_gp,  273.32730_gp/)
    vdwparams%coeffs(13301:13400)=(/  204.34870_gp,  152.35770_gp,  104.79180_gp,   75.70200_gp,   51.58810_gp,&
         316.07860_gp,  490.05630_gp,  490.22120_gp,  387.46930_gp,  313.76610_gp,&
         263.46580_gp,  213.88020_gp,  432.70980_gp,  854.43090_gp,  426.13390_gp,&
         411.39560_gp,  403.28980_gp,  400.28790_gp,  365.13690_gp,  337.22470_gp,&
         321.77920_gp,  314.79940_gp,  313.00420_gp,  290.33450_gp,  484.28840_gp,&
         419.85020_gp,  372.95910_gp,  338.65860_gp,  296.54870_gp,  514.97470_gp,&
         1042.33410_gp,  777.62370_gp,  564.81820_gp,  570.45940_gp,  530.64150_gp,&
         602.13570_gp,  458.94910_gp,  429.20470_gp,  398.10440_gp,  385.54200_gp,&
         378.94250_gp,  615.75540_gp,  550.49820_gp,  511.59560_gp,  481.59190_gp,&
         437.41870_gp,  603.46930_gp, 1360.11230_gp,  986.79870_gp,  599.98280_gp,&
         587.50600_gp,  568.61140_gp,  572.38360_gp,  554.75530_gp,  576.15700_gp,&
         540.43300_gp,  551.66610_gp,  512.23480_gp,  497.37970_gp,  494.37500_gp,&
         520.66780_gp,  477.43820_gp,  660.08460_gp,  608.52880_gp,  551.46290_gp,&
         561.43380_gp,   35.46460_gp,   23.81250_gp,  530.66030_gp,  310.52240_gp,&
         211.41050_gp,  144.24940_gp,  101.74200_gp,   77.55360_gp,   59.12040_gp,&
         45.72050_gp,  635.30720_gp,  493.55920_gp,  457.21240_gp,  362.83630_gp,&
         284.91700_gp,  237.70450_gp,  195.21900_gp,  160.53010_gp, 1042.06540_gp,&
         864.52710_gp,  715.43320_gp,  692.43390_gp,  634.88380_gp,  500.93670_gp,&
         547.65220_gp,  430.64570_gp,  456.96500_gp,  470.05570_gp,  361.10090_gp,&
         371.62530_gp,  439.32630_gp,  389.41390_gp,  333.96020_gp,  300.83180_gp/)
    vdwparams%coeffs(13401:13500)=(/  264.38450_gp,  230.89140_gp, 1169.37980_gp, 1030.36440_gp,  905.55540_gp,&
         815.80510_gp,  745.26460_gp,  577.58390_gp,  643.48640_gp,  492.30820_gp,&
         537.54740_gp,  499.07810_gp,  416.64440_gp,  440.47370_gp,  550.60050_gp,&
         510.94940_gp,  457.28590_gp,  425.51160_gp,  386.04100_gp,  348.27200_gp,&
         1425.75280_gp, 1313.05620_gp, 1157.08730_gp,  540.14830_gp, 1164.95320_gp,&
         1118.98400_gp, 1091.02750_gp, 1065.29550_gp, 1042.48060_gp,  823.51360_gp,&
         923.39310_gp,  891.29600_gp,  940.85810_gp,  920.91700_gp,  902.98170_gp,&
         892.22740_gp,  754.15480_gp,  743.92130_gp,  680.68180_gp,  576.55350_gp,&
         586.93470_gp,  533.25550_gp,  489.40580_gp,  407.85710_gp,  381.54050_gp,&
         392.59340_gp,  566.72940_gp,  555.34130_gp,  512.07880_gp,  489.50850_gp,&
         453.12570_gp,  416.56370_gp, 1349.17160_gp, 1300.73910_gp, 1152.96640_gp,&
         1040.29350_gp, 1031.52190_gp,  998.93380_gp, 1026.14250_gp,  994.07460_gp,&
         56.65670_gp,  180.83220_gp,  231.24090_gp,  178.34690_gp,  135.37830_gp,&
         94.96120_gp,   69.70800_gp,   48.43830_gp,  264.14790_gp,  408.11310_gp,&
         414.94650_gp,  335.78860_gp,  276.67050_gp,  235.08600_gp,  193.19560_gp,&
         366.55110_gp,  693.17420_gp,  368.02850_gp,  355.52300_gp,  348.41720_gp,&
         345.01650_gp,  318.30430_gp,  294.95420_gp,  281.45410_gp,  275.10220_gp,&
         271.23140_gp,  255.47700_gp,  413.18390_gp,  364.83360_gp,  328.42530_gp,&
         300.94220_gp,  266.17570_gp,  438.09620_gp,  843.71770_gp,  649.42160_gp,&
         486.94320_gp,  491.70440_gp,  459.26500_gp,  513.29930_gp,  401.62410_gp/)
    vdwparams%coeffs(13501:13600)=(/  376.01550_gp,  349.51790_gp,  337.98580_gp,  334.91410_gp,  525.56660_gp,&
         477.35010_gp,  448.45320_gp,  425.27250_gp,  389.63030_gp,  519.48750_gp,&
         1093.54920_gp,  821.79670_gp,  524.02380_gp,  513.21590_gp,  497.08550_gp,&
         499.23160_gp,  480.40280_gp,  502.16590_gp,  471.97170_gp,  479.76830_gp,&
         448.20990_gp,  435.52800_gp,  432.64200_gp,  454.17000_gp,  418.27710_gp,&
         562.87830_gp,  523.75380_gp,  479.21590_gp,  483.86730_gp,  426.19030_gp,&
         33.51170_gp,   22.67790_gp,  491.47210_gp,  289.82980_gp,  198.38070_gp,&
         136.00690_gp,   96.32430_gp,   73.66340_gp,   56.32680_gp,   43.67730_gp,&
         588.83180_gp,  460.10610_gp,  427.26470_gp,  340.27370_gp,  268.10830_gp,&
         224.24930_gp,  184.65170_gp,  152.22030_gp,  965.06260_gp,  804.38060_gp,&
         666.34240_gp,  645.67380_gp,  592.40830_gp,  467.75860_gp,  511.52360_gp,&
         402.59270_gp,  427.50640_gp,  439.43760_gp,  337.88310_gp,  348.20520_gp,&
         411.18360_gp,  365.42010_gp,  314.22920_gp,  283.59750_gp,  249.76170_gp,&
         218.57530_gp, 1083.61290_gp,  958.71100_gp,  844.20740_gp,  761.53710_gp,&
         696.35620_gp,  540.76150_gp,  601.99360_gp,  461.58650_gp,  503.78310_gp,&
         468.03160_gp,  390.92010_gp,  413.43530_gp,  515.70310_gp,  479.48370_gp,&
         430.06080_gp,  400.76390_gp,  364.20700_gp,  329.14390_gp, 1321.12730_gp,&
         1221.15640_gp, 1078.14120_gp,  507.92020_gp, 1084.20870_gp, 1041.82480_gp,&
         1015.88020_gp,  991.98570_gp,  970.80130_gp,  768.77680_gp,  859.38870_gp,&
         829.70230_gp,  876.55810_gp,  858.00850_gp,  841.34750_gp,  831.23230_gp/)
    vdwparams%coeffs(13601:13700)=(/  703.68220_gp,  695.44610_gp,  637.22330_gp,  540.56670_gp,  550.51890_gp,&
         500.82420_gp,  460.15330_gp,  384.10330_gp,  359.55540_gp,  370.03820_gp,&
         531.43810_gp,  521.44650_gp,  481.70840_gp,  461.01200_gp,  427.39750_gp,&
         393.50690_gp, 1252.44720_gp, 1210.99990_gp, 1075.30030_gp,  972.45520_gp,&
         963.45570_gp,  933.11220_gp,  957.33480_gp,  927.63110_gp,   53.32550_gp,&
         168.93310_gp,  216.25690_gp,  167.58690_gp,  127.66650_gp,   89.93570_gp,&
         66.27430_gp,   46.29090_gp,  246.67920_gp,  380.86140_gp,  388.03830_gp,&
         315.10570_gp,  260.39720_gp,  221.77880_gp,  182.73850_gp,  343.34310_gp,&
         645.30690_gp,  345.46620_gp,  333.78450_gp,  327.12730_gp,  323.85110_gp,&
         299.19940_gp,  277.42410_gp,  264.76150_gp,  258.76330_gp,  254.85680_gp,&
         240.52350_gp,  386.98910_gp,  342.59940_gp,  309.08100_gp,  283.69920_gp,&
         251.43160_gp,  410.68500_gp,  785.37910_gp,  606.90870_gp,  457.10620_gp,&
         461.61660_gp,  431.45100_gp,  481.28810_gp,  377.95500_gp,  353.98970_gp,&
         329.21750_gp,  318.31050_gp,  315.69280_gp,  492.62420_gp,  448.36200_gp,&
         421.88850_gp,  400.57160_gp,  367.57650_gp,  487.64860_gp, 1017.11830_gp,&
         767.71820_gp,  492.83570_gp,  482.68880_gp,  467.58270_gp,  469.43970_gp,&
         451.34990_gp,  472.10950_gp,  443.88260_gp,  450.95330_gp,  421.64500_gp,&
         409.75870_gp,  406.99850_gp,  426.96980_gp,  393.53520_gp,  527.55310_gp,&
         491.59800_gp,  450.48000_gp,  454.33550_gp,  401.48890_gp,  378.49390_gp,&
         32.09170_gp,   21.83380_gp,  465.20110_gp,  275.45630_gp,  189.14850_gp/)
    vdwparams%coeffs(13701:13800)=(/  130.07010_gp,   92.37100_gp,   70.79750_gp,   54.25160_gp,   42.14990_gp,&
         557.60480_gp,  437.00990_gp,  406.38720_gp,  324.32370_gp,  256.07430_gp,&
         214.53080_gp,  176.95160_gp,  146.11660_gp,  913.85330_gp,  763.27990_gp,&
         632.64330_gp,  613.44560_gp,  563.06260_gp,  444.82350_gp,  486.47550_gp,&
         383.12130_gp,  406.95490_gp,  418.13480_gp,  321.72230_gp,  331.77490_gp,&
         391.48440_gp,  348.44170_gp,  300.12110_gp,  271.18850_gp,  239.15510_gp,&
         209.57880_gp, 1026.50380_gp,  909.78800_gp,  801.97930_gp,  723.99270_gp,&
         662.40900_gp,  515.04420_gp,  573.09040_gp,  440.03780_gp,  480.10810_gp,&
         446.21570_gp,  372.86390_gp,  394.38550_gp,  491.29160_gp,  457.27960_gp,&
         410.67770_gp,  383.04320_gp,  348.47210_gp,  315.27170_gp, 1251.60790_gp,&
         1158.60400_gp, 1023.94960_gp,  485.00190_gp, 1029.17200_gp,  989.09200_gp,&
         964.50040_gp,  941.84450_gp,  921.75840_gp,  731.00210_gp,  815.92560_gp,&
         787.90030_gp,  832.46930_gp,  814.86390_gp,  799.06290_gp,  789.40200_gp,&
         668.92220_gp,  661.71440_gp,  606.82340_gp,  515.28520_gp,  524.88810_gp,&
         477.89070_gp,  439.38870_gp,  367.16870_gp,  343.85620_gp,  353.90700_gp,&
         506.69330_gp,  497.51750_gp,  460.09860_gp,  440.63800_gp,  408.89250_gp,&
         376.82880_gp, 1187.56080_gp, 1149.64350_gp, 1021.81320_gp,  925.32050_gp,&
         916.38160_gp,  887.57730_gp,  909.98500_gp,  881.86200_gp,   50.93660_gp,&
         160.64720_gp,  205.78590_gp,  159.92840_gp,  122.11080_gp,   86.26550_gp,&
         63.73640_gp,   44.67920_gp,  234.56040_gp,  361.98670_gp,  369.23330_gp/)
    vdwparams%coeffs(13801:13900)=(/  300.45310_gp,  248.74070_gp,  212.16920_gp,  175.11990_gp,  327.11390_gp,&
         612.66150_gp,  329.50500_gp,  318.42170_gp,  312.08150_gp,  308.91510_gp,&
         285.62920_gp,  264.94780_gp,  252.88280_gp,  247.14320_gp,  243.27100_gp,&
         229.84390_gp,  368.60280_gp,  326.82470_gp,  295.24020_gp,  271.28510_gp,&
         240.74110_gp,  391.46850_gp,  745.66370_gp,  577.49450_gp,  436.01580_gp,&
         440.33040_gp,  411.76280_gp,  458.83330_gp,  361.09350_gp,  338.29080_gp,&
         314.73130_gp,  304.28650_gp,  301.92540_gp,  469.50870_gp,  427.83570_gp,&
         402.95210_gp,  382.87770_gp,  351.68470_gp,  465.16200_gp,  965.40540_gp,&
         730.40540_gp,  470.64590_gp,  460.96660_gp,  446.58000_gp,  448.26620_gp,&
         430.78610_gp,  450.75260_gp,  423.90040_gp,  430.51250_gp,  402.72530_gp,&
         391.39780_gp,  388.73440_gp,  407.63660_gp,  375.90100_gp,  502.77570_gp,&
         468.91240_gp,  430.08890_gp,  433.48520_gp,  383.82310_gp,  362.01210_gp,&
         346.36100_gp,   30.67570_gp,   20.94100_gp,  446.19080_gp,  263.22300_gp,&
         180.65820_gp,  124.31330_gp,   88.39150_gp,   67.83720_gp,   52.06100_gp,&
         40.50920_gp,  534.80760_gp,  417.88500_gp,  388.34460_gp,  309.76480_gp,&
         244.60550_gp,  205.02050_gp,  169.22710_gp,  139.86320_gp,  877.57830_gp,&
         730.94660_gp,  605.54530_gp,  587.04740_gp,  538.75790_gp,  425.83430_gp,&
         465.40030_gp,  366.71160_gp,  389.19570_gp,  399.92890_gp,  307.93520_gp,&
         317.23870_gp,  374.21180_gp,  332.90660_gp,  286.73520_gp,  259.15840_gp,&
         228.64850_gp,  200.49290_gp,  985.72350_gp,  871.47550_gp,  767.63160_gp/)
    vdwparams%coeffs(13901:14000)=(/  692.75660_gp,  633.75930_gp,  492.82900_gp,  548.33760_gp,  421.10840_gp,&
         459.28620_gp,  426.86340_gp,  356.97430_gp,  377.30490_gp,  469.98930_gp,&
         437.19460_gp,  392.55510_gp,  366.15260_gp,  333.16950_gp,  301.52720_gp,&
         1201.87930_gp, 1110.27920_gp,  980.44160_gp,  463.71250_gp,  986.14070_gp,&
         947.59800_gp,  923.99680_gp,  902.25470_gp,  882.97590_gp,  699.90580_gp,&
         782.32470_gp,  755.33300_gp,  797.26220_gp,  780.37200_gp,  765.20560_gp,&
         755.96270_gp,  640.35930_gp,  632.92020_gp,  580.31360_gp,  492.90480_gp,&
         502.00340_gp,  457.04880_gp,  420.25190_gp,  351.35070_gp,  329.11380_gp,&
         338.64570_gp,  484.95430_gp,  475.88520_gp,  439.96460_gp,  421.33170_gp,&
         391.00970_gp,  360.42890_gp, 1139.64070_gp, 1101.29800_gp,  978.17100_gp,&
         885.38660_gp,  877.28480_gp,  849.72830_gp,  871.57240_gp,  844.56390_gp,&
         48.64780_gp,  153.52370_gp,  196.65700_gp,  152.78290_gp,  116.73350_gp,&
         82.56980_gp,   61.09710_gp,   42.93390_gp,  224.31640_gp,  346.12620_gp,&
         352.82440_gp,  286.98060_gp,  237.61790_gp,  202.77280_gp,  167.48190_gp,&
         313.03240_gp,  586.98040_gp,  315.00110_gp,  304.44790_gp,  298.41150_gp,&
         295.42220_gp,  273.04550_gp,  253.29860_gp,  241.79700_gp,  236.32330_gp,&
         232.69810_gp,  219.74110_gp,  352.30670_gp,  312.25540_gp,  282.08180_gp,&
         259.25510_gp,  230.16390_gp,  374.64460_gp,  714.61700_gp,  552.72720_gp,&
         416.93630_gp,  421.09580_gp,  393.81950_gp,  439.07200_gp,  345.35640_gp,&
         323.61650_gp,  301.13910_gp,  291.18100_gp,  288.80380_gp,  449.08360_gp/)
    vdwparams%coeffs(14001:14100)=(/  409.02650_gp,  385.17040_gp,  365.99200_gp,  336.23370_gp,  444.89470_gp,&
         925.64640_gp,  699.19040_gp,  449.99540_gp,  440.74490_gp,  426.99740_gp,&
         428.62430_gp,  412.03450_gp,  430.96490_gp,  405.31130_gp,  411.67450_gp,&
         385.04690_gp,  374.21380_gp,  371.66190_gp,  389.68240_gp,  359.37430_gp,&
         480.90340_gp,  448.47840_gp,  411.33730_gp,  414.67050_gp,  367.09540_gp,&
         346.30750_gp,  331.39460_gp,  317.15610_gp,   29.92030_gp,   20.58520_gp,&
         413.10610_gp,  250.30210_gp,  174.02940_gp,  120.76050_gp,   86.33820_gp,&
         66.49370_gp,   51.17020_gp,   39.89740_gp,  496.00030_gp,  395.61740_gp,&
         370.29810_gp,  298.04830_gp,  236.98310_gp,  199.43500_gp,  165.20990_gp,&
         136.93730_gp,  810.50400_gp,  686.33430_gp,  570.63520_gp,  554.96780_gp,&
         510.28000_gp,  403.44500_gp,  441.99080_gp,  348.47560_gp,  371.26550_gp,&
         380.80390_gp,  293.24220_gp,  303.81440_gp,  357.75570_gp,  320.47320_gp,&
         277.62260_gp,  251.75200_gp,  222.82360_gp,  195.92110_gp,  911.74980_gp,&
         817.78560_gp,  724.86420_gp,  656.64200_gp,  602.18630_gp,  470.17490_gp,&
         522.36590_gp,  402.91110_gp,  439.46430_gp,  409.02800_gp,  341.73120_gp,&
         362.19820_gp,  449.15100_gp,  420.20430_gp,  379.28220_gp,  354.82940_gp,&
         323.85100_gp,  293.89350_gp, 1112.03090_gp, 1039.69190_gp,  923.89520_gp,&
         447.56220_gp,  925.41710_gp,  890.16940_gp,  868.25330_gp,  848.03530_gp,&
         830.11710_gp,  662.49390_gp,  733.14800_gp,  708.72400_gp,  750.71540_gp,&
         734.93810_gp,  720.83100_gp,  711.92400_gp,  605.86800_gp,  602.48380_gp/)
    vdwparams%coeffs(14101:14200)=(/  554.37230_gp,  472.11020_gp,  481.47760_gp,  439.62380_gp,  405.14540_gp,&
         339.47380_gp,  318.26990_gp,  327.85490_gp,  464.06660_gp,  457.41320_gp,&
         424.89250_gp,  407.96160_gp,  379.73280_gp,  350.93490_gp, 1059.96930_gp,&
         1034.65730_gp,  924.29380_gp,  841.92530_gp,  831.61820_gp,  805.61970_gp,&
         823.08290_gp,  798.16520_gp,   47.16620_gp,  146.25750_gp,  187.79790_gp,&
         147.52400_gp,  113.36050_gp,   80.64710_gp,   59.93130_gp,   42.31950_gp,&
         213.16680_gp,  328.53430_gp,  336.98220_gp,  276.45860_gp,  230.25850_gp,&
         197.22520_gp,  163.49210_gp,  298.87190_gp,  551.42200_gp,  302.95130_gp,&
         292.86460_gp,  287.01690_gp,  283.89530_gp,  263.47450_gp,  244.69520_gp,&
         233.57610_gp,  228.21280_gp,  224.01780_gp,  212.72240_gp,  337.41410_gp,&
         301.05020_gp,  273.20270_gp,  251.83060_gp,  224.26200_gp,  358.23720_gp,&
         670.72660_gp,  525.11790_gp,  400.71700_gp,  404.64840_gp,  378.97070_gp,&
         420.21610_gp,  333.59820_gp,  312.69220_gp,  291.16220_gp,  281.37530_gp,&
         279.90990_gp,  430.03950_gp,  393.95270_gp,  372.39970_gp,  354.73990_gp,&
         326.81510_gp,  427.30670_gp,  866.55170_gp,  663.55950_gp,  434.45080_gp,&
         425.54240_gp,  412.37370_gp,  413.61800_gp,  396.55650_gp,  415.78790_gp,&
         391.29420_gp,  396.85720_gp,  371.98460_gp,  361.61110_gp,  359.07550_gp,&
         376.09030_gp,  347.34250_gp,  460.39360_gp,  430.75880_gp,  396.39370_gp,&
         398.43040_gp,  355.38470_gp,  335.59480_gp,  321.33060_gp,  307.47100_gp,&
         298.86160_gp,   45.65780_gp,   29.52890_gp,  818.93710_gp,  440.31470_gp/)
    vdwparams%coeffs(14201:14300)=(/  286.25250_gp,  188.89680_gp,  130.07910_gp,   97.55100_gp,   73.37840_gp,&
         56.17460_gp,  975.32750_gp,  710.24040_gp,  642.28990_gp,  493.56290_gp,&
         377.40890_gp,  309.61610_gp,  250.28720_gp,  203.07020_gp, 1619.93980_gp,&
         1276.88420_gp, 1044.62150_gp, 1000.55900_gp,  911.70000_gp,  718.35580_gp,&
         779.38520_gp,  611.31150_gp,  640.53000_gp,  663.11020_gp,  508.88790_gp,&
         513.81200_gp,  611.25030_gp,  528.55090_gp,  443.41250_gp,  394.08000_gp,&
         341.65520_gp,  294.73510_gp, 1809.68600_gp, 1524.55490_gp, 1313.13780_gp,&
         1168.21430_gp, 1058.41510_gp,  808.55780_gp,  905.69280_gp,  682.07110_gp,&
         744.93760_gp,  688.07030_gp,  575.99280_gp,  603.30050_gp,  766.65730_gp,&
         697.24530_gp,  611.99840_gp,  562.88000_gp,  504.38410_gp,  449.80100_gp,&
         2204.63590_gp, 1955.74900_gp, 1689.12850_gp,  725.31910_gp, 1723.50940_gp,&
         1649.97970_gp, 1607.24480_gp, 1568.08170_gp, 1533.30570_gp, 1184.47440_gp,&
         1372.19770_gp, 1319.37090_gp, 1376.84410_gp, 1346.92670_gp, 1319.66080_gp,&
         1305.15970_gp, 1086.52290_gp, 1050.10500_gp,  948.98330_gp,  796.01210_gp,&
         806.54170_gp,  725.04930_gp,  659.78220_gp,  544.89940_gp,  507.92100_gp,&
         520.62730_gp,  784.94960_gp,  757.25660_gp,  686.21680_gp,  649.45010_gp,&
         594.09560_gp,  540.36210_gp, 2053.09940_gp, 1917.21060_gp, 1668.13700_gp,&
         1473.55230_gp, 1476.40210_gp, 1428.94000_gp, 1487.12040_gp, 1437.17680_gp,&
         74.82650_gp,  254.79910_gp,  323.09330_gp,  239.25450_gp,  177.39330_gp,&
         121.33560_gp,   87.29830_gp,   59.25500_gp,  375.28460_gp,  582.20680_gp/)
    vdwparams%coeffs(14301:14400)=(/  579.56000_gp,  454.68720_gp,  366.16070_gp,  306.33030_gp,  247.76880_gp,&
         511.64600_gp, 1022.40250_gp,  500.90730_gp,  483.49150_gp,  474.05210_gp,&
         470.89950_gp,  428.04350_gp,  394.93390_gp,  376.86900_gp,  368.82120_gp,&
         367.76620_gp,  339.40980_gp,  571.25400_gp,  492.37270_gp,  435.51000_gp,&
         394.31080_gp,  344.21330_gp,  608.11700_gp, 1247.83980_gp,  922.81660_gp,&
         664.15540_gp,  670.92560_gp,  623.34490_gp,  710.55610_gp,  537.41150_gp,&
         502.47720_gp,  465.83490_gp,  451.42980_gp,  442.50690_gp,  726.60300_gp,&
         646.40140_gp,  598.61930_gp,  562.15510_gp,  509.17420_gp,  709.76410_gp,&
         1630.75530_gp, 1171.81560_gp,  702.86450_gp,  688.22190_gp,  665.93910_gp,&
         670.82470_gp,  651.73940_gp,  675.37850_gp,  633.15820_gp,  647.13510_gp,&
         599.75110_gp,  582.22680_gp,  578.81040_gp,  610.17640_gp,  558.77490_gp,&
         778.70390_gp,  715.93350_gp,  646.97310_gp,  660.39260_gp,  565.40350_gp,&
         530.44050_gp,  505.86340_gp,  484.03950_gp,  463.95530_gp,  778.61880_gp,&
         42.43400_gp,   27.90730_gp,  688.07560_gp,  389.43570_gp,  259.45980_gp,&
         174.07890_gp,  121.24080_gp,   91.60350_gp,   69.32070_gp,   53.31050_gp,&
         821.97960_gp,  622.84710_gp,  570.78670_gp,  446.34760_gp,  345.97620_gp,&
         286.14180_gp,  233.03850_gp,  190.23650_gp, 1351.46350_gp, 1101.52570_gp,&
         907.51920_gp,  874.24690_gp,  799.42710_gp,  629.78200_gp,  686.87860_gp,&
         539.01900_gp,  569.43390_gp,  587.46310_gp,  450.54860_gp,  460.28610_gp,&
         545.68040_gp,  478.36640_gp,  405.93740_gp,  363.17430_gp,  316.92710_gp/)
    vdwparams%coeffs(14401:14500)=(/  274.97470_gp, 1513.10000_gp, 1313.05530_gp, 1144.54560_gp, 1025.56740_gp,&
         933.37410_gp,  718.32310_gp,  802.42420_gp,  609.24720_gp,  665.83460_gp,&
         616.75520_gp,  514.89300_gp,  542.71640_gp,  683.52170_gp,  628.90890_gp,&
         557.84160_gp,  516.14750_gp,  465.37870_gp,  417.35550_gp, 1842.77110_gp,&
         1676.70030_gp, 1465.77040_gp,  659.62730_gp, 1483.42340_gp, 1423.29130_gp,&
         1387.27420_gp, 1354.18230_gp, 1324.82370_gp, 1036.27450_gp, 1176.16350_gp,&
         1133.25710_gp, 1193.43790_gp, 1167.94680_gp, 1144.88170_gp, 1131.73970_gp,&
         950.08080_gp,  929.48800_gp,  845.77140_gp,  712.70800_gp,  724.21000_gp,&
         654.78150_gp,  598.54060_gp,  496.39040_gp,  463.49930_gp,  476.29930_gp,&
         701.16210_gp,  682.85450_gp,  624.80060_gp,  594.46280_gp,  547.11560_gp,&
         500.27110_gp, 1733.13330_gp, 1653.73570_gp, 1455.24790_gp, 1301.18700_gp,&
         1295.36640_gp, 1254.14220_gp, 1295.52580_gp, 1253.84620_gp,   68.70090_gp,&
         226.05860_gp,  287.75870_gp,  217.82360_gp,  163.33240_gp,  113.06740_gp,&
         82.11150_gp,   56.34450_gp,  331.28460_gp,  512.67930_gp,  516.59790_gp,&
         412.11230_gp,  335.77950_gp,  283.02980_gp,  230.64950_gp,  455.12070_gp,&
         880.99550_gp,  452.42230_gp,  436.71780_gp,  428.07670_gp,  424.46770_gp,&
         389.17610_gp,  359.85440_gp,  343.32520_gp,  335.76340_gp,  332.69860_gp,&
         310.58260_gp,  511.74530_gp,  446.93270_gp,  398.94250_gp,  363.33700_gp,&
         319.18450_gp,  542.38170_gp, 1072.82400_gp,  812.13110_gp,  598.68520_gp,&
         604.72750_gp,  563.26320_gp,  634.69440_gp,  489.44530_gp,  457.83760_gp/)
    vdwparams%coeffs(14501:14600)=(/  424.97290_gp,  411.37260_gp,  405.83020_gp,  650.42200_gp,  585.43790_gp,&
         546.38190_gp,  515.67690_gp,  469.75780_gp,  638.71850_gp, 1393.80820_gp,&
         1028.72910_gp,  639.57980_gp,  626.32970_gp,  606.35080_gp,  609.75710_gp,&
         589.24030_gp,  613.69780_gp,  576.11550_gp,  586.97030_gp,  546.49440_gp,&
         530.80340_gp,  527.48330_gp,  554.93190_gp,  509.65210_gp,  696.21110_gp,&
         644.34370_gp,  586.25070_gp,  594.80320_gp,  517.23940_gp,  486.23500_gp,&
         464.24430_gp,  444.00340_gp,  427.98590_gp,  697.95700_gp,  632.99910_gp,&
         41.59640_gp,   27.60930_gp,  639.13910_gp,  371.40510_gp,  250.81810_gp,&
         169.84210_gp,  119.03470_gp,   90.30520_gp,   68.55930_gp,   52.85010_gp,&
         764.72640_gp,  591.32420_gp,  545.86550_gp,  430.92720_gp,  336.54820_gp,&
         279.61500_gp,  228.67090_gp,  187.30610_gp, 1251.98640_gp, 1037.18980_gp,&
         857.56950_gp,  828.73520_gp,  759.24150_gp,  598.25830_gp,  654.12160_gp,&
         513.59410_gp,  544.75080_gp,  560.93530_gp,  430.21600_gp,  442.12850_gp,&
         523.24780_gp,  462.07160_gp,  394.58520_gp,  354.32910_gp,  310.33860_gp,&
         270.12000_gp, 1403.64110_gp, 1235.57780_gp, 1083.73330_gp,  974.77890_gp,&
         889.33630_gp,  687.26090_gp,  766.55860_gp,  584.63360_gp,  638.98360_gp,&
         592.76900_gp,  494.33450_gp,  522.59600_gp,  655.09190_gp,  606.40770_gp,&
         540.93950_gp,  502.16500_gp,  454.32770_gp,  408.72120_gp, 1709.87640_gp,&
         1574.44450_gp, 1384.98620_gp,  638.93210_gp, 1395.80400_gp, 1340.62350_gp,&
         1307.08100_gp, 1276.22200_gp, 1248.85820_gp,  983.47160_gp, 1105.01920_gp/)
    vdwparams%coeffs(14601:14700)=(/ 1065.97900_gp, 1126.77670_gp, 1102.90190_gp, 1081.38800_gp, 1068.67850_gp,&
         901.28220_gp,  887.29660_gp,  810.34560_gp,  684.71570_gp,  696.76070_gp,&
         631.89030_gp,  579.01680_gp,  481.34860_gp,  449.88030_gp,  462.87080_gp,&
         672.90150_gp,  658.43910_gp,  605.54990_gp,  577.80740_gp,  533.56470_gp,&
         489.31160_gp, 1616.29770_gp, 1557.89290_gp, 1378.79840_gp, 1240.73280_gp,&
         1231.20180_gp, 1192.20190_gp, 1226.63580_gp, 1188.05280_gp,   66.88570_gp,&
         215.99080_gp,  275.62110_gp,  211.10780_gp,  159.30450_gp,  111.00780_gp,&
         81.02380_gp,   55.92160_gp,  315.69470_gp,  487.98840_gp,  494.88070_gp,&
         398.38950_gp,  326.70150_gp,  276.53670_gp,  226.30380_gp,  435.78520_gp,&
         829.78660_gp,  436.55630_gp,  421.48580_gp,  413.07940_gp,  409.21130_gp,&
         376.85240_gp,  348.88250_gp,  332.84050_gp,  325.38910_gp,  321.33350_gp,&
         301.79750_gp,  491.62320_gp,  432.41530_gp,  387.92100_gp,  354.46370_gp,&
         312.48750_gp,  520.15670_gp, 1009.52910_gp,  773.62030_gp,  577.27760_gp,&
         583.00940_gp,  543.83670_gp,  609.23200_gp,  474.49870_gp,  443.99260_gp,&
         412.40560_gp,  398.94980_gp,  394.86670_gp,  624.69100_gp,  565.75460_gp,&
         530.22150_gp,  501.81220_gp,  458.58080_gp,  615.48530_gp, 1308.13720_gp,&
         978.87630_gp,  619.71040_gp,  606.90590_gp,  587.70740_gp,  590.49670_gp,&
         568.98830_gp,  594.17810_gp,  558.18530_gp,  567.79730_gp,  529.88060_gp,&
         514.80880_gp,  511.48010_gp,  537.47970_gp,  494.40600_gp,  668.52110_gp,&
         620.87010_gp,  566.90460_gp,  573.34500_gp,  502.67320_gp,  473.05980_gp/)
    vdwparams%coeffs(14701:14800)=(/  451.95060_gp,  432.15130_gp,  417.77890_gp,  671.00110_gp,  612.23660_gp,&
         594.04570_gp,   40.78400_gp,   27.29270_gp,  603.78280_gp,  356.67730_gp,&
         243.25530_gp,  165.91800_gp,  116.89210_gp,   88.99600_gp,   67.76250_gp,&
         52.35200_gp,  723.18320_gp,  566.18910_gp,  525.32650_gp,  417.51820_gp,&
         327.94200_gp,  273.46760_gp,  224.41640_gp,  184.36010_gp, 1181.98630_gp,&
         988.27090_gp,  818.91600_gp,  793.12700_gp,  727.54910_gp,  573.61570_gp,&
         627.98280_gp,  493.46120_gp,  524.58090_gp,  539.44100_gp,  413.95970_gp,&
         426.96110_gp,  504.66270_gp,  447.93420_gp,  384.30870_gp,  346.11240_gp,&
         304.04010_gp,  265.34650_gp, 1326.57080_gp, 1177.09020_gp, 1036.56600_gp,&
         934.73590_gp,  854.28840_gp,  662.24370_gp,  737.78930_gp,  564.59860_gp,&
         616.91080_gp,  572.89550_gp,  477.65720_gp,  505.75670_gp,  631.86020_gp,&
         587.25260_gp,  525.98230_gp,  489.50040_gp,  444.05350_gp,  400.48270_gp,&
         1616.68890_gp, 1498.30600_gp, 1323.17990_gp,  620.88910_gp, 1330.06570_gp,&
         1278.24590_gp, 1246.47090_gp, 1217.21250_gp, 1191.27650_gp,  942.48920_gp,&
         1052.60250_gp, 1016.22720_gp, 1075.82350_gp, 1053.12210_gp, 1032.72450_gp,&
         1020.38150_gp,  863.29250_gp,  853.31590_gp,  781.30510_gp,  661.65650_gp,&
         673.88930_gp,  612.48620_gp,  562.23340_gp,  468.34590_gp,  438.07020_gp,&
         451.02990_gp,  649.93220_gp,  637.84100_gp,  588.69220_gp,  562.90490_gp,&
         521.11770_gp,  478.99830_gp, 1533.01370_gp, 1485.68370_gp, 1319.57920_gp,&
         1192.55370_gp, 1181.06740_gp, 1143.77230_gp, 1173.67000_gp, 1137.27960_gp/)
    vdwparams%coeffs(14801:14900)=(/   65.21460_gp,  207.72960_gp,  265.61960_gp,  205.16820_gp,  155.61880_gp,&
         109.02870_gp,   79.92300_gp,   55.44500_gp,  303.13050_gp,  468.21060_gp,&
         476.85900_gp,  386.38850_gp,  318.41620_gp,  270.43650_gp,  222.08000_gp,&
         420.21160_gp,  791.41640_gp,  422.99690_gp,  408.51000_gp,  400.32420_gp,&
         396.32850_gp,  366.05400_gp,  339.20620_gp,  323.62470_gp,  316.30020_gp,&
         311.63550_gp,  293.90630_gp,  474.78450_gp,  419.70260_gp,  377.92620_gp,&
         346.23010_gp,  306.09970_gp,  502.20180_gp,  962.52050_gp,  743.49910_gp,&
         559.25120_gp,  564.72600_gp,  527.39750_gp,  588.57680_gp,  461.47450_gp,&
         431.95260_gp,  401.45890_gp,  388.17970_gp,  385.01900_gp,  603.40750_gp,&
         548.77650_gp,  515.84790_gp,  489.23280_gp,  448.19060_gp,  596.16690_gp,&
         1245.57400_gp,  940.22950_gp,  602.35500_gp,  589.93220_gp,  571.39080_gp,&
         573.76380_gp,  551.77910_gp,  577.20750_gp,  542.52310_gp,  551.28420_gp,&
         515.27560_gp,  500.71660_gp,  497.39780_gp,  522.20000_gp,  480.93310_gp,&
         645.84430_gp,  601.27470_gp,  550.39340_gp,  555.42350_gp,  489.76750_gp,&
         461.32380_gp,  440.97230_gp,  421.63540_gp,  408.42210_gp,  648.79600_gp,&
         594.36400_gp,  577.97360_gp,  563.24630_gp,   38.51210_gp,   26.11110_gp,&
         541.78230_gp,  327.04600_gp,  226.12850_gp,  155.86820_gp,  110.68200_gp,&
         84.74010_gp,   64.82870_gp,   50.27410_gp,  649.91410_gp,  517.13720_gp,&
         483.14350_gp,  387.56890_gp,  306.89040_gp,  257.30000_gp,  212.24570_gp,&
         175.15290_gp, 1060.53650_gp,  897.15640_gp,  745.54500_gp,  724.28050_gp/)
    vdwparams%coeffs(14901:15000)=(/  665.56000_gp,  525.34660_gp,  575.94770_gp,  453.23530_gp,  483.11010_gp,&
         495.87210_gp,  380.99950_gp,  394.73600_gp,  465.67730_gp,  416.19670_gp,&
         359.43220_gp,  325.08320_gp,  286.81650_gp,  251.32750_gp, 1192.13680_gp,&
         1068.48880_gp,  945.94700_gp,  855.98480_gp,  784.21710_gp,  610.69650_gp,&
         679.19130_gp,  522.32000_gp,  570.34140_gp,  530.43260_gp,  442.33610_gp,&
         469.16940_gp,  583.36540_gp,  545.06230_gp,  490.91940_gp,  458.47770_gp,&
         417.50700_gp,  377.92660_gp, 1453.84080_gp, 1358.33910_gp, 1205.79060_gp,&
         579.11420_gp, 1208.03290_gp, 1161.85760_gp, 1133.21400_gp, 1106.80650_gp,&
         1083.40690_gp,  862.66010_gp,  955.96420_gp,  923.91610_gp,  979.57880_gp,&
         959.00710_gp,  940.59640_gp,  929.08800_gp,  789.50440_gp,  784.48200_gp,&
         720.81970_gp,  612.48290_gp,  624.50880_gp,  569.34740_gp,  523.95200_gp,&
         437.82610_gp,  410.01810_gp,  422.47690_gp,  601.41130_gp,  592.44910_gp,&
         549.42480_gp,  526.88570_gp,  489.51330_gp,  451.44550_gp, 1384.40680_gp,&
         1350.74260_gp, 1205.35870_gp, 1095.73070_gp, 1082.47810_gp, 1048.46140_gp,&
         1072.03110_gp, 1039.42560_gp,   61.07570_gp,  190.88070_gp,  244.79400_gp,&
         191.30870_gp,  146.21180_gp,  103.28040_gp,   76.21860_gp,   53.30920_gp,&
         278.00940_gp,  428.91920_gp,  439.34570_gp,  359.20390_gp,  298.08060_gp,&
         254.42840_gp,  210.02330_gp,  387.91110_gp,  719.71410_gp,  392.89950_gp,&
         379.62690_gp,  371.98630_gp,  367.97730_gp,  341.17250_gp,  316.58390_gp,&
         302.08580_gp,  295.15410_gp,  289.91570_gp,  274.92600_gp,  438.90000_gp/)
    vdwparams%coeffs(15001:15100)=(/  390.64520_gp,  353.61220_gp,  325.17910_gp,  288.69910_gp,  464.45980_gp,&
         875.07120_gp,  683.19770_gp,  519.43980_gp,  524.44490_gp,  490.61730_gp,&
         544.79300_gp,  431.00110_gp,  403.66930_gp,  375.52180_gp,  362.88510_gp,&
         360.89310_gp,  558.18460_gp,  510.50570_gp,  481.83620_gp,  458.32530_gp,&
         421.36600_gp,  553.69970_gp, 1130.61240_gp,  863.38050_gp,  562.02230_gp,&
         550.46390_gp,  533.32520_gp,  535.11070_gp,  513.28320_gp,  538.11960_gp,&
         506.16170_gp,  513.61300_gp,  481.06970_gp,  467.60130_gp,  464.39200_gp,&
         486.88020_gp,  449.18810_gp,  597.61830_gp,  558.26740_gp,  512.82040_gp,&
         515.97410_gp,  458.59230_gp,  432.54320_gp,  413.80590_gp,  395.69230_gp,&
         384.29790_gp,  601.17960_gp,  553.69590_gp,  540.01160_gp,  527.42780_gp,&
         495.44740_gp,   47.70160_gp,   30.98270_gp,  821.28980_gp,  449.07760_gp,&
         295.29920_gp,  196.40820_gp,  135.86190_gp,  102.09130_gp,   76.82810_gp,&
         58.75990_gp,  978.72720_gp,  721.62860_gp,  656.46580_gp,  508.48770_gp,&
         391.45280_gp,  322.43460_gp,  261.53630_gp,  212.69100_gp, 1623.52410_gp,&
         1290.60070_gp, 1058.12330_gp, 1015.75590_gp,  926.71700_gp,  730.34690_gp,&
         793.66630_gp,  622.73900_gp,  654.32150_gp,  676.36930_gp,  519.02340_gp,&
         526.39910_gp,  626.02660_gp,  544.61780_gp,  459.45110_gp,  409.71550_gp,&
         356.33240_gp,  308.18050_gp, 1815.46720_gp, 1540.73550_gp, 1332.50230_gp,&
         1188.57810_gp, 1078.78130_gp,  826.61910_gp,  924.79920_gp,  698.68860_gp,&
         762.96140_gp,  705.40290_gp,  589.92230_gp,  619.17780_gp,  784.33640_gp/)
    vdwparams%coeffs(15101:15200)=(/  716.71230_gp,  632.13270_gp,  583.14470_gp,  524.15300_gp,  468.71230_gp,&
         2213.35130_gp, 1974.69390_gp, 1712.24690_gp,  748.65870_gp, 1741.75230_gp,&
         1668.37740_gp, 1625.40630_gp, 1586.00200_gp, 1551.02880_gp, 1203.74580_gp,&
         1386.13640_gp, 1333.85270_gp, 1393.99110_gp, 1363.81590_gp, 1336.40290_gp,&
         1321.44620_gp, 1103.60080_gp, 1071.72560_gp,  971.15790_gp,  816.39790_gp,&
         827.99830_gp,  746.02310_gp,  680.07030_gp,  562.55190_gp,  524.62170_gp,&
         538.24540_gp,  803.93350_gp,  778.27560_gp,  708.23270_gp,  671.99770_gp,&
         616.53950_gp,  562.22180_gp, 2067.17660_gp, 1939.93920_gp, 1693.57020_gp,&
         1502.46400_gp, 1501.92520_gp, 1453.68780_gp, 1508.33320_gp, 1458.33420_gp,&
         77.70550_gp,  260.24430_gp,  330.86670_gp,  247.38310_gp,  184.44310_gp,&
         126.75340_gp,   91.43880_gp,   62.12190_gp,  382.27080_gp,  592.97950_gp,&
         593.19050_gp,  468.99060_gp,  379.86120_gp,  318.96080_gp,  258.86520_gp,&
         523.46900_gp, 1034.54240_gp,  515.44090_gp,  497.59520_gp,  487.74150_gp,&
         484.06800_gp,  441.52590_gp,  407.74620_gp,  389.03940_gp,  380.57220_gp,&
         378.31310_gp,  350.99820_gp,  585.93320_gp,  508.06740_gp,  451.40640_gp,&
         409.92700_gp,  358.93240_gp,  623.02850_gp, 1262.32480_gp,  941.31820_gp,&
         683.39070_gp,  690.13390_gp,  641.88840_gp,  728.47950_gp,  554.98020_gp,&
         518.92730_gp,  481.22640_gp,  465.93230_gp,  457.99090_gp,  744.56260_gp,&
         665.67970_gp,  618.72370_gp,  582.51630_gp,  529.12450_gp,  730.44860_gp,&
         1647.86740_gp, 1194.78610_gp,  725.65350_gp,  710.54890_gp,  687.68570_gp/)
    vdwparams%coeffs(15201:15300)=(/  692.26220_gp,  670.81490_gp,  696.84430_gp,  653.57440_gp,  667.19620_gp,&
         619.47600_gp,  601.50620_gp,  597.87770_gp,  629.70880_gp,  577.40830_gp,&
         798.73450_gp,  736.26140_gp,  667.10040_gp,  679.12190_gp,  585.15020_gp,&
         549.32980_gp,  524.01710_gp,  501.14870_gp,  481.50990_gp,  798.48790_gp,&
         719.09940_gp,  693.21150_gp,  671.61670_gp,  623.96020_gp,  822.04820_gp,&
         83.60900_gp,   50.79650_gp, 2063.78360_gp,  940.42740_gp,  567.81070_gp,&
         355.48600_gp,  235.51600_gp,  171.84720_gp,  126.19530_gp,   94.70350_gp,&
         2433.40340_gp, 1557.05010_gp, 1354.45410_gp,  985.80450_gp,  722.44190_gp,&
         577.42680_gp,  455.25930_gp,  361.44930_gp, 4218.95280_gp, 2955.57810_gp,&
         2362.88550_gp, 2225.63170_gp, 2006.05930_gp, 1585.09840_gp, 1688.31640_gp,&
         1324.54090_gp, 1348.78630_gp, 1411.19510_gp, 1087.88780_gp, 1055.84500_gp,&
         1272.66820_gp, 1053.21460_gp,  852.17950_gp,  741.55110_gp,  629.20990_gp,&
         532.28900_gp, 4692.05940_gp, 3553.39920_gp, 2947.79080_gp, 2564.32900_gp,&
         2292.27090_gp, 1712.60020_gp, 1934.99170_gp, 1420.59090_gp, 1544.98010_gp,&
         1413.65040_gp, 1196.83650_gp, 1224.70850_gp, 1604.20500_gp, 1403.79770_gp,&
         1191.53350_gp, 1075.37550_gp,  944.59510_gp,  826.92720_gp, 5752.88770_gp,&
         4638.31340_gp, 3853.20930_gp, 1423.51320_gp, 4040.79230_gp, 3832.79720_gp,&
         3725.32980_gp, 3627.67870_gp, 3540.79340_gp, 2635.98410_gp, 3283.09640_gp,&
         3145.72610_gp, 3144.67020_gp, 3072.03860_gp, 3004.59070_gp, 2975.80180_gp,&
         2422.32570_gp, 2246.05740_gp, 1985.85620_gp, 1644.99660_gp, 1650.24540_gp/)
    vdwparams%coeffs(15301:15400)=(/ 1455.76630_gp, 1305.13720_gp, 1063.88060_gp,  986.00810_gp, 1001.43120_gp,&
         1636.77890_gp, 1527.26550_gp, 1340.82190_gp, 1247.63380_gp, 1119.02600_gp,&
         1000.32620_gp, 5185.16360_gp, 4461.39470_gp, 3741.10850_gp, 3178.81670_gp,&
         3255.39130_gp, 3146.16150_gp, 3348.38180_gp, 3220.23880_gp,  142.55890_gp,&
         539.86610_gp,  679.29310_gp,  468.40690_gp,  335.16710_gp,  219.93730_gp,&
         152.84480_gp,   99.28020_gp,  807.46090_gp, 1265.54660_gp, 1211.43030_gp,&
         902.34430_gp,  700.28070_gp,  571.80250_gp,  450.94310_gp, 1080.60420_gp,&
         2400.40190_gp, 1003.12780_gp,  971.17090_gp,  952.38070_gp,  951.58400_gp,&
         839.82780_gp,  769.65310_gp,  735.32570_gp,  721.11550_gp,  733.54280_gp,&
         652.14040_gp, 1176.38960_gp,  972.40690_gp,  835.49110_gp,  742.34170_gp,&
         634.68070_gp, 1275.12970_gp, 2958.54010_gp, 2035.20980_gp, 1343.48360_gp,&
         1353.91760_gp, 1251.52810_gp, 1485.82930_gp, 1050.03620_gp,  980.72260_gp,&
         905.51290_gp,  880.03940_gp,  843.98030_gp, 1500.56320_gp, 1284.94120_gp,&
         1160.59610_gp, 1072.74230_gp,  953.86310_gp, 1447.52050_gp, 3962.70430_gp,&
         2613.04930_gp, 1376.87620_gp, 1347.48290_gp, 1301.84470_gp, 1319.88440_gp,&
         1303.67580_gp, 1329.33230_gp, 1240.18400_gp, 1282.71440_gp, 1169.18840_gp,&
         1133.06970_gp, 1127.90110_gp, 1196.81680_gp, 1085.79330_gp, 1623.25980_gp,&
         1458.94960_gp, 1287.87830_gp, 1341.15770_gp, 1089.61390_gp, 1014.87600_gp,&
         963.92390_gp,  923.74840_gp,  868.31320_gp, 1603.64020_gp, 1381.83250_gp,&
         1302.20540_gp, 1243.04790_gp, 1132.07280_gp, 1624.53320_gp, 3787.59430_gp/)
    vdwparams%coeffs(15401:15500)=(/   51.10730_gp,   33.63980_gp,  806.02430_gp,  463.02300_gp,  310.58590_gp,&
         209.18580_gp,  145.95370_gp,  110.33330_gp,   83.47610_gp,   64.14430_gp,&
         963.53500_gp,  738.54190_gp,  679.47260_gp,  533.94560_gp,  415.34530_gp,&
         344.13770_gp,  280.65670_gp,  229.28490_gp, 1579.66110_gp, 1299.84280_gp,&
         1073.02330_gp, 1035.32520_gp,  947.61900_gp,  746.30200_gp,  815.29190_gp,&
         639.68180_gp,  677.45550_gp,  698.23130_gp,  535.18630_gp,  548.68120_gp,&
         650.23470_gp,  572.21090_gp,  487.04720_gp,  436.43770_gp,  381.38700_gp,&
         331.23010_gp, 1769.72280_gp, 1548.74580_gp, 1354.58060_gp, 1216.18340_gp,&
         1108.18640_gp,  854.38760_gp,  953.76820_gp,  725.52720_gp,  793.14780_gp,&
         735.16550_gp,  613.02710_gp,  647.39340_gp,  813.66740_gp,  751.10850_gp,&
         668.12780_gp,  619.15920_gp,  559.09410_gp,  502.01500_gp, 2155.59720_gp,&
         1975.19960_gp, 1732.70320_gp,  789.54030_gp, 1749.11230_gp, 1679.21180_gp,&
         1636.98580_gp, 1598.16510_gp, 1563.73640_gp, 1227.35700_gp, 1385.10290_gp,&
         1335.42490_gp, 1409.89020_gp, 1379.92030_gp, 1352.86540_gp, 1337.15100_gp,&
         1125.13500_gp, 1104.75470_gp, 1007.08820_gp,  849.53530_gp,  863.91120_gp,&
         782.18420_gp,  715.74490_gp,  593.97880_gp,  554.73410_gp,  570.48480_gp,&
         834.86440_gp,  815.22210_gp,  747.86620_gp,  712.56000_gp,  656.80530_gp,&
         601.29860_gp, 2033.08230_gp, 1951.57780_gp, 1722.60390_gp, 1545.23830_gp,&
         1535.41510_gp, 1486.61930_gp, 1532.29630_gp, 1483.58560_gp,   82.52730_gp,&
         268.98140_gp,  342.82380_gp,  261.02340_gp,  196.21380_gp,  136.08950_gp/)
    vdwparams%coeffs(15501:15600)=(/   98.90800_gp,   67.85860_gp,  393.43830_gp,  608.67090_gp,  615.47850_gp,&
         493.27920_gp,  403.12400_gp,  340.35860_gp,  277.75450_gp,  541.50520_gp,&
         1039.35420_gp,  540.62540_gp,  521.82910_gp,  511.41880_gp,  506.81950_gp,&
         465.76900_gp,  430.87900_gp,  411.02320_gp,  401.87340_gp,  397.45630_gp,&
         372.25130_gp,  610.37340_gp,  535.02330_gp,  478.72210_gp,  436.61290_gp,&
         384.06800_gp,  645.77530_gp, 1264.87950_gp,  963.84860_gp,  715.06360_gp,&
         722.18160_gp,  673.03520_gp,  756.00910_gp,  585.89660_gp,  548.01980_gp,&
         508.73460_gp,  492.22520_gp,  486.48840_gp,  775.11410_gp,  699.94150_gp,&
         654.63420_gp,  618.66240_gp,  564.35520_gp,  762.65600_gp, 1640.78290_gp,&
         1220.19750_gp,  765.64160_gp,  749.79030_gp,  725.95040_gp,  729.71000_gp,&
         704.02290_gp,  734.40270_gp,  689.61450_gp,  702.03220_gp,  654.40840_gp,&
         635.70180_gp,  631.67050_gp,  664.25610_gp,  610.46730_gp,  829.74450_gp,&
         769.21560_gp,  701.02290_gp,  710.06080_gp,  619.88720_gp,  582.90540_gp,&
         556.59850_gp,  532.13080_gp,  513.69080_gp,  831.88530_gp,  756.84580_gp,&
         733.24790_gp,  712.62070_gp,  664.75080_gp,  858.84400_gp, 1630.42110_gp,&
         906.04480_gp,   58.04200_gp,   37.71710_gp, 1008.87190_gp,  548.42030_gp,&
         360.05530_gp,  239.21200_gp,  165.41170_gp,  124.32800_gp,   93.64080_gp,&
         71.71850_gp, 1201.75620_gp,  881.96280_gp,  801.68460_gp,  620.20400_gp,&
         476.97090_gp,  392.62180_gp,  318.31370_gp,  258.81200_gp, 2000.74170_gp,&
         1579.50490_gp, 1294.15880_gp, 1241.94230_gp, 1132.86670_gp,  892.96930_gp/)
    vdwparams%coeffs(15601:15700)=(/  969.99930_gp,  761.15260_gp,  799.25150_gp,  826.38550_gp,  634.31240_gp,&
         642.80130_gp,  764.33880_gp,  664.30070_gp,  559.96260_gp,  499.07290_gp,&
         433.84090_gp,  375.08740_gp, 2237.27760_gp, 1886.08870_gp, 1629.43760_gp,&
         1452.64480_gp, 1318.14880_gp, 1009.51600_gp, 1129.83590_gp,  853.13360_gp,&
         931.56900_gp,  861.20120_gp,  720.51370_gp,  755.92110_gp,  958.11210_gp,&
         874.79670_gp,  771.06630_gp,  711.00220_gp,  638.78920_gp,  570.99520_gp,&
         2730.68360_gp, 2419.12440_gp, 2094.78410_gp,  913.10340_gp, 2134.16630_gp,&
         2042.69320_gp, 1989.90800_gp, 1941.52700_gp, 1898.58440_gp, 1472.32090_gp,&
         1701.11370_gp, 1637.99750_gp, 1705.71920_gp, 1668.72050_gp, 1635.08310_gp,&
         1616.85300_gp, 1350.55540_gp, 1309.12600_gp, 1185.71000_gp,  996.56680_gp,&
         1010.60390_gp,  910.25670_gp,  829.62710_gp,  686.20440_gp,  639.97490_gp,&
         656.57070_gp,  982.06910_gp,  950.02070_gp,  864.04130_gp,  819.55770_gp,&
         751.62660_gp,  685.16600_gp, 2544.06950_gp, 2374.80280_gp, 2071.39670_gp,&
         1835.92960_gp, 1836.56100_gp, 1777.46360_gp, 1845.32020_gp, 1783.91740_gp,&
         94.64180_gp,  317.81250_gp,  403.99690_gp,  301.52640_gp,  224.58760_gp,&
         154.26790_gp,  111.30120_gp,   75.74610_gp,  467.15670_gp,  724.73020_gp,&
         724.26630_gp,  571.94340_gp,  462.83220_gp,  388.40280_gp,  315.07380_gp,&
         639.24470_gp, 1268.21100_gp,  628.63700_gp,  607.26950_gp,  595.22520_gp,&
         590.85660_gp,  538.66860_gp,  497.44160_gp,  474.68740_gp,  464.39490_gp,&
         461.83290_gp,  428.19990_gp,  715.20850_gp,  619.60740_gp,  550.13410_gp/)
    vdwparams%coeffs(15701:15800)=(/  499.33630_gp,  437.01980_gp,  760.67450_gp, 1548.28770_gp, 1152.22220_gp,&
         833.61520_gp,  841.45560_gp,  783.11230_gp,  889.73240_gp,  676.85180_gp,&
         632.96580_gp,  587.02360_gp,  568.47740_gp,  558.62930_gp,  909.26210_gp,&
         812.33800_gp,  754.64760_gp,  710.21200_gp,  644.84300_gp,  890.97610_gp,&
         2024.93160_gp, 1463.50690_gp,  885.03310_gp,  866.59780_gp,  838.70440_gp,&
         844.47740_gp,  818.50700_gp,  849.98270_gp,  797.15610_gp,  814.06210_gp,&
         755.50250_gp,  733.57320_gp,  729.16890_gp,  768.08740_gp,  704.15860_gp,&
         975.50990_gp,  898.59560_gp,  813.74230_gp,  828.81920_gp,  713.50780_gp,&
         669.76280_gp,  638.93520_gp,  611.16560_gp,  587.08670_gp,  975.11530_gp,&
         877.52310_gp,  845.61490_gp,  819.04220_gp,  760.65040_gp, 1002.64610_gp,&
         1992.66470_gp, 1047.43020_gp, 1224.59330_gp,   52.26450_gp,   34.33290_gp,&
         890.70150_gp,  485.02270_gp,  320.84970_gp,  214.65050_gp,  149.30470_gp,&
         112.73380_gp,   85.26210_gp,   65.53310_gp, 1060.80400_gp,  778.88060_gp,&
         710.34880_gp,  552.19130_gp,  426.80750_gp,  352.64770_gp,  286.99980_gp,&
         234.18470_gp, 1778.94680_gp, 1393.96830_gp, 1142.59350_gp, 1097.97050_gp,&
         1002.24710_gp,  790.96650_gp,  859.10440_gp,  674.98040_gp,  709.00470_gp,&
         732.40820_gp,  562.99640_gp,  571.31500_gp,  678.70160_gp,  591.90790_gp,&
         500.94070_gp,  447.74380_gp,  390.43140_gp,  338.58330_gp, 1991.34470_gp,&
         1665.65470_gp, 1440.77830_gp, 1286.02970_gp, 1168.34620_gp,  896.88390_gp,&
         1003.03230_gp,  759.31940_gp,  828.41280_gp,  766.40860_gp,  641.82000_gp/)
    vdwparams%coeffs(15801:15900)=(/  673.42870_gp,  851.58740_gp,  779.29980_gp,  689.07930_gp,  636.81840_gp,&
         573.61210_gp,  514.07080_gp, 2438.15580_gp, 2138.29770_gp, 1852.69460_gp,&
         815.84440_gp, 1888.99250_gp, 1805.77160_gp, 1759.00050_gp, 1716.14510_gp,&
         1678.11740_gp, 1304.86510_gp, 1510.40780_gp, 1457.67600_gp, 1507.57700_gp,&
         1474.81150_gp, 1445.06280_gp, 1428.77660_gp, 1197.51470_gp, 1160.74980_gp,&
         1052.98590_gp,  886.98680_gp,  899.86720_gp,  811.84020_gp,  741.02110_gp,&
         614.28030_gp,  573.39400_gp,  588.42050_gp,  874.53780_gp,  846.98840_gp,&
         772.32730_gp,  733.86420_gp,  674.59880_gp,  616.36260_gp, 2265.52210_gp,&
         2100.07160_gp, 1833.44060_gp, 1628.05850_gp, 1628.19960_gp, 1575.76490_gp,&
         1633.63030_gp, 1579.41500_gp,   84.72380_gp,  281.45990_gp,  358.59140_gp,&
         269.24380_gp,  201.57990_gp,  139.31210_gp,  101.05790_gp,   69.27070_gp,&
         413.51280_gp,  641.46750_gp,  642.36160_gp,  509.71630_gp,  414.27030_gp,&
         348.84970_gp,  284.07250_gp,  568.38310_gp, 1124.36320_gp,  559.73470_gp,&
         541.64230_gp,  530.75760_gp,  526.67650_gp,  480.98950_gp,  444.58610_gp,&
         424.38350_gp,  415.09660_gp,  412.10000_gp,  383.22820_gp,  635.62870_gp,&
         552.67960_gp,  492.29290_gp,  447.97180_gp,  393.23670_gp,  677.15660_gp,&
         1374.34370_gp, 1026.03970_gp,  742.79230_gp,  748.95660_gp,  698.92350_gp,&
         792.75350_gp,  605.35870_gp,  566.44570_gp,  525.69600_gp,  508.89220_gp,&
         500.75690_gp,  808.85400_gp,  724.62440_gp,  674.73780_gp,  636.19670_gp,&
         579.01430_gp,  794.65520_gp, 1803.19140_gp, 1305.08800_gp,  790.95680_gp/)
    vdwparams%coeffs(15901:16000)=(/  774.48190_gp,  749.70660_gp,  754.73340_gp,  730.34120_gp,  759.25810_gp,&
         712.36590_gp,  727.23060_gp,  675.39950_gp,  655.89880_gp,  651.87600_gp,&
         686.05020_gp,  629.62830_gp,  868.95920_gp,  801.43340_gp,  726.95200_gp,&
         739.34260_gp,  639.31090_gp,  600.57540_gp,  573.28530_gp,  548.47560_gp,&
         527.67310_gp,  868.65890_gp,  783.63670_gp,  756.34500_gp,  733.58480_gp,&
         682.69100_gp,  894.35100_gp, 1771.60490_gp,  935.86320_gp, 1093.49140_gp,&
         980.37580_gp,   46.75590_gp,   31.07460_gp,  731.98050_gp,  419.75640_gp,&
         282.43120_gp,  191.02230_gp,  133.87230_gp,  101.59960_gp,   77.17640_gp,&
         59.52710_gp,  875.18610_gp,  669.58510_gp,  616.54260_gp,  485.33720_gp,&
         378.50250_gp,  314.37000_gp,  257.08750_gp,  210.63660_gp, 1438.79940_gp,&
         1179.82480_gp,  973.72410_gp,  939.93750_gp,  860.49210_gp,  678.50520_gp,&
         740.61370_gp,  581.83780_gp,  615.69900_gp,  634.36640_gp,  487.01500_gp,&
         499.03190_gp,  590.84610_gp,  520.52680_gp,  443.90020_gp,  398.44930_gp,&
         348.91220_gp,  303.70810_gp, 1612.54300_gp, 1406.49780_gp, 1230.11360_gp,&
         1104.76660_gp, 1007.12820_gp,  777.57670_gp,  867.55080_gp,  661.00390_gp,&
         721.98170_gp,  669.45830_gp,  559.04400_gp,  589.89000_gp,  740.43240_gp,&
         683.80180_gp,  609.02530_gp,  565.01830_gp,  510.96120_gp,  459.56760_gp,&
         1965.24280_gp, 1794.90590_gp, 1574.09590_gp,  719.85720_gp, 1589.85950_gp,&
         1525.98220_gp, 1487.54080_gp, 1452.19650_gp, 1420.84750_gp, 1116.11380_gp,&
         1260.97660_gp, 1215.93980_gp, 1280.84100_gp, 1253.54930_gp, 1228.91830_gp/)
    vdwparams%coeffs(16001:16100)=(/ 1214.57730_gp, 1022.68030_gp, 1003.89750_gp,  915.71760_gp,  773.55940_gp,&
         786.64370_gp,  712.80370_gp,  652.78320_gp,  542.66550_gp,  507.16770_gp,&
         521.45390_gp,  760.86200_gp,  742.90150_gp,  682.12270_gp,  650.42490_gp,&
         600.26280_gp,  550.29950_gp, 1852.50590_gp, 1773.47910_gp, 1565.14570_gp,&
         1404.88130_gp, 1396.37900_gp, 1352.07130_gp, 1393.36210_gp, 1349.04950_gp,&
         75.20660_gp,  244.03010_gp,  311.32040_gp,  237.66340_gp,  179.25610_gp,&
         124.89310_gp,   91.18630_gp,   62.98670_gp,  357.15820_gp,  552.32750_gp,&
         558.71280_gp,  448.61130_gp,  367.44990_gp,  310.93340_gp,  254.43970_gp,&
         493.05900_gp,  945.12240_gp,  492.17970_gp,  475.31010_gp,  465.86480_gp,&
         461.66860_gp,  424.43750_gp,  392.86320_gp,  374.86590_gp,  366.51980_gp,&
         362.35940_gp,  339.63850_gp,  554.81070_gp,  486.95430_gp,  436.39410_gp,&
         398.61480_gp,  351.33960_gp,  588.39800_gp, 1150.90480_gp,  877.28230_gp,&
         651.37000_gp,  657.82330_gp,  613.55160_gp,  688.95110_gp,  534.68200_gp,&
         500.38440_gp,  464.79060_gp,  449.70340_gp,  444.50890_gp,  705.53710_gp,&
         637.59540_gp,  596.86560_gp,  564.59640_gp,  515.74020_gp,  695.01880_gp,&
         1494.27710_gp, 1110.91990_gp,  698.16140_gp,  683.72500_gp,  662.06630_gp,&
         665.40380_gp,  641.79300_gp,  669.49760_gp,  628.84090_gp,  640.05220_gp,&
         596.80940_gp,  579.78670_gp,  576.05760_gp,  605.40230_gp,  556.73610_gp,&
         755.52630_gp,  700.91760_gp,  639.34100_gp,  647.27080_gp,  566.14580_gp,&
         532.75200_gp,  508.97860_gp,  486.81770_gp,  470.17820_gp,  758.08150_gp/)
    vdwparams%coeffs(16101:16200)=(/  690.14330_gp,  668.95140_gp,  650.52430_gp,  607.48900_gp,  782.64720_gp,&
         1484.26300_gp,  826.01520_gp,  954.79790_gp,  854.10020_gp,  753.86170_gp,&
         45.12360_gp,   30.09780_gp,  695.01280_gp,  401.90350_gp,  271.45390_gp,&
         184.10950_gp,  129.30040_gp,   98.28120_gp,   74.75840_gp,   57.72900_gp,&
         831.49520_gp,  640.31040_gp,  590.77440_gp,  466.30230_gp,  364.44330_gp,&
         303.11500_gp,  248.21860_gp,  203.61560_gp, 1364.35030_gp, 1125.42920_gp,&
         929.88130_gp,  898.45340_gp,  822.98300_gp,  648.99350_gp,  708.90860_gp,&
         557.05720_gp,  590.15540_gp,  607.71500_gp,  466.59660_gp,  478.90520_gp,&
         566.62400_gp,  500.22560_gp,  427.35080_gp,  384.01570_gp,  336.65520_gp,&
         293.34880_gp, 1529.69780_gp, 1341.32750_gp, 1175.37250_gp, 1056.82350_gp,&
         964.12660_gp,  745.33540_gp,  831.17090_gp,  634.18960_gp,  692.70260_gp,&
         642.61150_gp,  536.50370_gp,  566.57800_gp,  710.09360_gp,  656.92100_gp,&
         586.01680_gp,  544.18120_gp,  492.61320_gp,  443.48890_gp, 1863.88600_gp,&
         1710.41530_gp, 1502.97870_gp,  692.46950_gp, 1516.01600_gp, 1455.72450_gp,&
         1419.19920_gp, 1385.59970_gp, 1355.80150_gp, 1067.19440_gp, 1201.55300_gp,&
         1158.92430_gp, 1222.84540_gp, 1196.86070_gp, 1173.43690_gp, 1159.64570_gp,&
         977.68160_gp,  961.56840_gp,  878.07320_gp,  742.36820_gp,  755.24290_gp,&
         684.98790_gp,  627.77710_gp,  522.30280_gp,  488.30230_gp,  502.21500_gp,&
         730.00580_gp,  713.76740_gp,  656.32120_gp,  626.32590_gp,  578.57790_gp,&
         530.88320_gp, 1760.24620_gp, 1691.72540_gp, 1495.73510_gp, 1345.23450_gp/)
    vdwparams%coeffs(16201:16300)=(/ 1335.79380_gp, 1293.50310_gp, 1331.47370_gp, 1289.43400_gp,   72.42720_gp,&
         233.77270_gp,  298.40760_gp,  228.59970_gp,  172.75920_gp,  120.63570_gp,&
         88.24390_gp,   61.09930_gp,  341.93020_gp,  528.51140_gp,  535.61950_gp,&
         431.17410_gp,  353.82780_gp,  299.79380_gp,  245.65840_gp,  472.72210_gp,&
         901.29920_gp,  472.93970_gp,  456.72220_gp,  447.64490_gp,  443.50180_gp,&
         408.25230_gp,  378.02020_gp,  360.70350_gp,  352.64190_gp,  348.33050_gp,&
         327.02340_gp,  532.35460_gp,  468.17890_gp,  420.16770_gp,  384.16910_gp,&
         338.97840_gp,  564.38650_gp, 1097.12550_gp,  839.43700_gp,  625.73340_gp,&
         631.97490_gp,  589.65940_gp,  660.95050_gp,  514.51210_gp,  481.57200_gp,&
         447.43080_gp,  432.85140_gp,  428.22300_gp,  677.04980_gp,  612.92050_gp,&
         574.43630_gp,  543.80490_gp,  497.21050_gp,  667.50050_gp, 1422.87350_gp,&
         1062.51530_gp,  671.66460_gp,  657.79210_gp,  637.00760_gp,  640.04320_gp,&
         616.87680_gp,  643.93890_gp,  604.97690_gp,  615.45730_gp,  574.28150_gp,&
         557.94580_gp,  554.32010_gp,  582.34520_gp,  535.79050_gp,  724.82070_gp,&
         673.16650_gp,  614.71110_gp,  621.76580_gp,  545.16120_gp,  513.20580_gp,&
         490.42130_gp,  469.07220_gp,  453.40690_gp,  727.71340_gp,  663.65960_gp,&
         643.84950_gp,  626.49490_gp,  585.54410_gp,  751.73280_gp, 1415.39250_gp,&
         794.62310_gp,  916.95930_gp,  820.52080_gp,  725.37300_gp,  698.16440_gp,&
         52.39190_gp,   33.86150_gp,  977.20270_gp,  511.20280_gp,  330.02550_gp,&
         217.09250_gp,  149.29170_gp,  111.90540_gp,   84.16280_gp,   64.43330_gp/)
    vdwparams%coeffs(16301:16400)=(/ 1161.72080_gp,  827.48750_gp,  745.01090_gp,  569.32030_gp,  433.98740_gp,&
         355.59410_gp,  287.20610_gp,  232.91960_gp, 1950.98200_gp, 1501.03940_gp,&
         1223.42610_gp, 1169.41900_gp, 1064.06850_gp,  839.49690_gp,  907.88740_gp,&
         712.75770_gp,  743.44550_gp,  770.54070_gp,  592.37050_gp,  594.75930_gp,&
         708.53630_gp,  609.79120_gp,  510.11170_gp,  452.83130_gp,  392.22240_gp,&
         338.15470_gp, 2178.79960_gp, 1795.03590_gp, 1537.27940_gp, 1363.48500_gp,&
         1233.46170_gp,  940.31770_gp, 1054.18410_gp,  792.00470_gp,  863.90680_gp,&
         797.15370_gp,  668.98050_gp,  698.11240_gp,  890.02560_gp,  805.67160_gp,&
         704.98790_gp,  647.53750_gp,  579.56520_gp,  516.40770_gp, 2660.74750_gp,&
         2310.74450_gp, 1983.10260_gp,  836.46730_gp, 2033.39390_gp, 1942.71580_gp,&
         1891.63270_gp, 1844.89320_gp, 1803.37800_gp, 1386.49190_gp, 1627.14710_gp,&
         1564.80290_gp, 1616.28730_gp, 1580.74840_gp, 1548.27410_gp, 1531.51850_gp,&
         1272.11370_gp, 1221.52520_gp, 1101.10390_gp,  923.10550_gp,  934.07850_gp,&
         838.12920_gp,  761.68140_gp,  628.77650_gp,  585.95420_gp,  599.87110_gp,&
         911.89380_gp,  875.76580_gp,  791.06510_gp,  747.64610_gp,  683.00400_gp,&
         620.65460_gp, 2460.87670_gp, 2258.42010_gp, 1953.63690_gp, 1716.72670_gp,&
         1725.93380_gp, 1670.06200_gp, 1743.23510_gp, 1683.41510_gp,   86.00520_gp,&
         295.68750_gp,  374.96860_gp,  275.66970_gp,  204.03100_gp,  139.32300_gp,&
         100.15200_gp,   67.95140_gp,  436.51000_gp,  678.20670_gp,  671.69720_gp,&
         524.29270_gp,  421.07980_gp,  351.86230_gp,  284.33720_gp,  595.04560_gp/)
    vdwparams%coeffs(16401:16500)=(/ 1207.48960_gp,  578.29290_gp,  558.86520_gp,  547.92340_gp,  544.65600_gp,&
         493.39110_gp,  455.04230_gp,  434.39050_gp,  425.20030_gp,  424.83420_gp,&
         390.59480_gp,  661.38520_gp,  567.67290_gp,  500.98750_gp,  453.12420_gp,&
         395.18880_gp,  706.95870_gp, 1477.10990_gp, 1080.92530_gp,  768.28830_gp,&
         775.51800_gp,  720.86910_gp,  826.04480_gp,  619.78620_gp,  579.63490_gp,&
         537.29800_gp,  520.79250_gp,  509.25590_gp,  842.21840_gp,  746.08910_gp,&
         689.36250_gp,  646.64990_gp,  585.07090_gp,  822.45460_gp, 1941.17440_gp,&
         1375.79290_gp,  810.46870_gp,  793.53070_gp,  767.76280_gp,  774.00610_gp,&
         753.18700_gp,  779.07900_gp,  730.04770_gp,  747.27740_gp,  691.19800_gp,&
         670.89470_gp,  667.02610_gp,  703.43090_gp,  643.73580_gp,  904.52750_gp,&
         829.41860_gp,  747.66610_gp,  764.81820_gp,  651.46000_gp,  610.88630_gp,&
         582.48590_gp,  557.57750_gp,  533.44460_gp,  902.97300_gp,  805.72640_gp,&
         773.05280_gp,  746.71410_gp,  691.11280_gp,  924.97250_gp, 1898.50760_gp,&
         959.15900_gp, 1130.77420_gp, 1008.87290_gp,  874.61490_gp,  838.95690_gp,&
         1051.24680_gp,   12.59310_gp,    8.16280_gp,  194.70360_gp,  113.96770_gp,&
         76.72300_gp,   51.56280_gp,   35.79020_gp,   26.89540_gp,   20.20830_gp,&
         15.41760_gp,  232.80600_gp,  181.16530_gp,  167.29290_gp,  131.89200_gp,&
         102.61760_gp,   84.86910_gp,   69.00690_gp,   56.15170_gp,  379.42250_gp,&
         316.51020_gp,  261.94430_gp,  253.05430_gp,  231.80650_gp,  182.15580_gp,&
         199.64260_gp,  156.28450_gp,  166.21660_gp,  171.20430_gp,  130.80120_gp/)
    vdwparams%coeffs(16501:16600)=(/  134.79330_gp,  159.92730_gp,  141.14890_gp,  120.22240_gp,  107.63710_gp,&
         93.89110_gp,   81.33060_gp,  425.17040_gp,  376.64680_gp,  330.71900_gp,&
         297.47500_gp,  271.26060_gp,  209.09590_gp,  233.46060_gp,  177.51740_gp,&
         194.39840_gp,  180.21780_gp,  149.71980_gp,  158.69370_gp,  199.41110_gp,&
         184.68490_gp,  164.53480_gp,  152.49560_gp,  137.61920_gp,  123.40700_gp,&
         517.87260_gp,  479.36610_gp,  422.27540_gp,  194.17280_gp,  424.78830_gp,&
         408.12030_gp,  397.94790_gp,  388.59250_gp,  380.30110_gp,  299.32330_gp,&
         335.36650_gp,  323.59410_gp,  343.28500_gp,  336.05070_gp,  329.53550_gp,&
         325.68440_gp,  274.59690_gp,  270.80850_gp,  247.15510_gp,  208.29050_gp,&
         212.01830_gp,  192.03170_gp,  175.71480_gp,  145.52190_gp,  135.78630_gp,&
         139.83510_gp,  204.18140_gp,  200.02350_gp,  183.83790_gp,  175.25160_gp,&
         161.52400_gp,  147.75410_gp,  490.03150_gp,  474.51000_gp,  420.36790_gp,&
         378.11690_gp,  374.72790_gp,  362.78690_gp,  373.04730_gp,  361.35530_gp,&
         20.40390_gp,   66.20090_gp,   84.40000_gp,   64.43010_gp,   48.30720_gp,&
         33.32740_gp,   24.05920_gp,   16.31850_gp,   96.51250_gp,  149.38320_gp,&
         151.59310_gp,  121.83330_gp,   99.56560_gp,   83.91810_gp,   68.28060_gp,&
         132.47510_gp,  252.57780_gp,  132.98680_gp,  128.28620_gp,  125.67780_gp,&
         124.46450_gp,  114.64650_gp,  106.02860_gp,  101.08010_gp,   98.80040_gp,&
         97.53400_gp,   91.62270_gp,  150.21480_gp,  131.99340_gp,  118.15200_gp,&
         107.67130_gp,   94.55150_gp,  157.94790_gp,  306.94330_gp,  235.52320_gp/)
    vdwparams%coeffs(16601:16700)=(/  175.65640_gp,  177.33570_gp,  165.20660_gp,  185.02440_gp,  143.88510_gp,&
         134.45760_gp,  124.71510_gp,  120.59750_gp,  119.46230_gp,  190.13130_gp,&
         172.17670_gp,  161.22340_gp,  152.38060_gp,  138.92730_gp,  187.13180_gp,&
         397.21940_gp,  297.90180_gp,  188.28490_gp,  184.38050_gp,  178.50850_gp,&
         179.39000_gp,  172.78660_gp,  180.61150_gp,  169.57170_gp,  172.51880_gp,&
         160.96160_gp,  156.37090_gp,  155.38670_gp,  163.48550_gp,  150.21190_gp,&
         203.49990_gp,  188.78180_gp,  172.12060_gp,  174.12660_gp,  152.26480_gp,&
         143.05840_gp,  136.50000_gp,  130.34490_gp,  125.97100_gp,  203.71920_gp,&
         185.86380_gp,  180.30690_gp,  175.30470_gp,  163.50870_gp,  210.81250_gp,&
         395.75600_gp,  222.89780_gp,  256.94260_gp,  229.39980_gp,  202.81970_gp,&
         195.12320_gp,  234.41020_gp,   55.13640_gp,   12.49520_gp,    8.25860_gp,&
         188.76110_gp,  110.41850_gp,   74.98650_gp,   50.90900_gp,   35.67780_gp,&
         27.02470_gp,   20.45920_gp,   15.71270_gp,  225.78430_gp,  175.39140_gp,&
         162.44730_gp,  128.74590_gp,  100.84320_gp,   83.88470_gp,   68.63310_gp,&
         56.19260_gp,  370.15350_gp,  306.90780_gp,  253.95340_gp,  245.66470_gp,&
         225.17830_gp,  177.40000_gp,  194.13780_gp,  152.38870_gp,  161.87200_gp,&
         166.56140_gp,  127.66860_gp,  131.52520_gp,  155.80800_gp,  137.99560_gp,&
         118.13540_gp,  106.21400_gp,   93.10030_gp,   81.05180_gp,  415.22010_gp,&
         365.65710_gp,  321.26520_gp,  289.28690_gp,  264.12850_gp,  204.30920_gp,&
         227.79590_gp,  173.87670_gp,  190.04450_gp,  176.33810_gp,  146.90170_gp/)
    vdwparams%coeffs(16701:16800)=(/  155.47780_gp,  194.76490_gp,  180.69960_gp,  161.55700_gp,  150.17540_gp,&
         136.02390_gp,  122.46340_gp,  506.50540_gp,  465.97660_gp,  410.51230_gp,&
         190.74630_gp,  413.12090_gp,  396.75810_gp,  386.83630_gp,  377.70920_gp,&
         369.61950_gp,  291.64910_gp,  327.23330_gp,  315.88750_gp,  333.55970_gp,&
         326.49630_gp,  320.14370_gp,  316.35220_gp,  267.24090_gp,  263.68480_gp,&
         241.07990_gp,  203.85670_gp,  207.52370_gp,  188.34080_gp,  172.66610_gp,&
         143.53170_gp,  134.11890_gp,  138.06770_gp,  200.08980_gp,  196.08950_gp,&
         180.68820_gp,  172.61360_gp,  159.58520_gp,  146.47210_gp,  478.96760_gp,&
         461.47150_gp,  408.83210_gp,  368.46880_gp,  365.24790_gp,  353.64010_gp,&
         363.28030_gp,  351.90410_gp,   20.05970_gp,   64.25530_gp,   82.14300_gp,&
         63.16850_gp,   47.74620_gp,   33.26770_gp,   24.24160_gp,   16.65310_gp,&
         93.71590_gp,  144.97830_gp,  147.37660_gp,  119.08890_gp,   97.89470_gp,&
         82.94970_gp,   67.91320_gp,  129.59630_gp,  245.94980_gp,  130.17250_gp,&
         125.70170_gp,  123.15290_gp,  121.93400_gp,  112.46480_gp,  104.14370_gp,&
         99.33430_gp,   97.08160_gp,   95.70070_gp,   90.13510_gp,  146.48900_gp,&
         129.22270_gp,  116.15350_gp,  106.24910_gp,   93.73740_gp,  154.78270_gp,&
         299.29140_gp,  230.13840_gp,  172.18100_gp,  173.79460_gp,  162.19830_gp,&
         181.38850_gp,  141.62010_gp,  132.47760_gp,  123.02270_gp,  118.92360_gp,&
         117.88430_gp,  185.85760_gp,  168.72480_gp,  158.40370_gp,  150.08700_gp,&
         137.30200_gp,  183.65680_gp,  388.02680_gp,  291.27440_gp,  185.01940_gp/)
    vdwparams%coeffs(16801:16900)=(/  181.19130_gp,  175.46970_gp,  176.26590_gp,  169.58240_gp,  177.36130_gp,&
         166.62020_gp,  169.42280_gp,  158.21700_gp,  153.73020_gp,  152.73000_gp,&
         160.46200_gp,  147.66360_gp,  199.19290_gp,  185.14520_gp,  169.18110_gp,&
         170.89180_gp,  150.17550_gp,  141.30850_gp,  134.97080_gp,  128.97340_gp,&
         124.81680_gp,  199.64430_gp,  182.52700_gp,  177.34000_gp,  172.70710_gp,&
         161.53370_gp,  206.89790_gp,  386.45800_gp,  218.94300_gp,  252.14850_gp,&
         225.71000_gp,  199.62770_gp,  192.15040_gp,  230.00030_gp,   53.99190_gp,&
         53.11280_gp,    9.42030_gp,    6.50270_gp,  123.77070_gp,   76.75520_gp,&
         54.08540_gp,   37.84190_gp,   27.17040_gp,   20.95970_gp,   16.13120_gp,&
         12.56290_gp,  148.76180_gp,  120.73860_gp,  113.85600_gp,   92.49310_gp,&
         74.07840_gp,   62.59530_gp,   52.02080_gp,   43.20640_gp,  242.63660_gp,&
         208.02230_gp,  173.46300_gp,  169.20640_gp,  155.83870_gp,  123.25910_gp,&
         135.29940_gp,  106.73260_gp,  114.09130_gp,  116.80720_gp,   89.95150_gp,&
         93.68700_gp,  110.27700_gp,   99.46980_gp,   86.68950_gp,   78.88380_gp,&
         70.03520_gp,   61.72540_gp,  273.37570_gp,  247.82850_gp,  220.88300_gp,&
         200.78170_gp,  184.54770_gp,  144.63550_gp,  160.45200_gp,  124.23890_gp,&
         135.47600_gp,  126.23620_gp,  105.35720_gp,  111.92520_gp,  138.28310_gp,&
         130.07660_gp,  118.03370_gp,  110.77470_gp,  101.42160_gp,   92.28760_gp,&
         333.81300_gp,  314.68370_gp,  281.15040_gp,  139.18110_gp,  280.46490_gp,&
         269.98650_gp,  263.39220_gp,  257.30290_gp,  251.90990_gp,  202.29780_gp/)
    vdwparams%coeffs(16901:17000)=(/  222.09340_gp,  214.93990_gp,  228.08570_gp,  223.31660_gp,  219.07220_gp,&
         216.30440_gp,  184.87050_gp,  184.93620_gp,  170.72850_gp,  145.77740_gp,&
         148.83210_gp,  136.24100_gp,  125.80040_gp,  105.59200_gp,   99.04540_gp,&
         102.12660_gp,  143.07200_gp,  141.57460_gp,  132.11550_gp,  127.19520_gp,&
         118.75500_gp,  110.03140_gp,  319.57380_gp,  314.09480_gp,  281.85090_gp,&
         258.14870_gp,  254.26850_gp,  246.32970_gp,  250.69380_gp,  243.24780_gp,&
         14.76230_gp,   44.93200_gp,   57.88050_gp,   45.96320_gp,   35.52190_gp,&
         25.38270_gp,   18.90340_gp,   13.35290_gp,   65.27390_gp,  100.58170_gp,&
         103.79070_gp,   85.90680_gp,   71.99140_gp,   61.89160_gp,   51.47210_gp,&
         92.02400_gp,  167.41450_gp,   93.89110_gp,   90.78380_gp,   88.94290_gp,&
         87.88890_gp,   81.87310_gp,   76.11320_gp,   72.64430_gp,   70.94440_gp,&
         69.40580_gp,   66.28380_gp,  104.18810_gp,   93.58930_gp,   85.33790_gp,&
         78.90290_gp,   70.47440_gp,  110.48560_gp,  203.57060_gp,  161.11150_gp,&
         124.19200_gp,  125.36250_gp,  117.55240_gp,  129.70790_gp,  103.79560_gp,&
         97.29610_gp,   90.62490_gp,   87.49510_gp,   87.29150_gp,  132.65050_gp,&
         122.20300_gp,  115.97420_gp,  110.77110_gp,  102.34890_gp,  132.46620_gp,&
         262.61810_gp,  203.48290_gp,  135.13800_gp,  132.36910_gp,  128.30130_gp,&
         128.59630_gp,  122.94040_gp,  129.24610_gp,  121.68990_gp,  123.26190_gp,&
         115.76120_gp,  112.55880_gp,  111.74980_gp,  116.93130_gp,  108.14300_gp,&
         142.18370_gp,  133.42480_gp,  123.13620_gp,  123.40630_gp,  110.82750_gp/)
    vdwparams%coeffs(17001:17100)=(/  104.72640_gp,  100.30120_gp,   95.92270_gp,   93.46280_gp,  143.24960_gp,&
         132.82830_gp,  130.04040_gp,  127.39760_gp,  120.19670_gp,  149.32790_gp,&
         264.05250_gp,  159.78300_gp,  181.82240_gp,  163.66940_gp,  146.25060_gp,&
         141.11910_gp,  164.50530_gp,   39.28140_gp,   38.98340_gp,   29.36020_gp,&
         7.66100_gp,    5.45590_gp,   94.00600_gp,   59.59830_gp,   42.83480_gp,&
         30.53050_gp,   22.27940_gp,   17.40930_gp,   13.56130_gp,   10.67400_gp,&
         113.30310_gp,   93.41300_gp,   88.84070_gp,   73.08580_gp,   59.28590_gp,&
         50.59460_gp,   42.48160_gp,   35.63340_gp,  185.21810_gp,  160.27500_gp,&
         134.03140_gp,  131.30150_gp,  121.21190_gp,   96.25470_gp,  105.60800_gp,&
         83.69160_gp,   89.53450_gp,   91.42950_gp,   70.76390_gp,   73.92100_gp,&
         86.61940_gp,   78.82350_gp,   69.37700_gp,   63.59190_gp,   56.91870_gp,&
         50.57460_gp,  209.23780_gp,  191.12680_gp,  171.35740_gp,  156.45400_gp,&
         144.31220_gp,  114.00040_gp,  126.08010_gp,   98.47280_gp,  107.10840_gp,&
         100.03700_gp,   83.78750_gp,   88.98580_gp,  109.08250_gp,  103.21130_gp,&
         94.36540_gp,   89.03590_gp,   82.03740_gp,   75.14270_gp,  255.74990_gp,&
         242.55280_gp,  217.89480_gp,  111.28270_gp,  216.78200_gp,  208.83570_gp,&
         203.77290_gp,  199.08850_gp,  194.93970_gp,  157.92210_gp,  172.08570_gp,&
         166.73560_gp,  176.70260_gp,  173.01200_gp,  169.74290_gp,  167.52470_gp,&
         144.02080_gp,  144.79850_gp,  134.33990_gp,  115.44060_gp,  117.98450_gp,&
         108.51360_gp,  100.61040_gp,   85.01970_gp,   79.96250_gp,   82.45400_gp/)
    vdwparams%coeffs(17101:17200)=(/  113.48820_gp,  112.68290_gp,  105.79110_gp,  102.26460_gp,   96.00710_gp,&
         89.45890_gp,  246.00480_gp,  242.92720_gp,  219.10640_gp,  202.23940_gp,&
         198.81690_gp,  192.68400_gp,  195.32220_gp,  189.64690_gp,   11.82140_gp,&
         35.01940_gp,   45.31640_gp,   36.60630_gp,   28.69380_gp,   20.85010_gp,&
         15.76310_gp,   11.35710_gp,   50.86720_gp,   78.17600_gp,   81.20260_gp,&
         68.05040_gp,   57.66320_gp,   50.02980_gp,   42.03670_gp,   72.68180_gp,&
         129.58410_gp,   74.55870_gp,   72.18280_gp,   70.73470_gp,   69.84710_gp,&
         65.34150_gp,   60.89200_gp,   58.16110_gp,   56.78680_gp,   55.37890_gp,&
         53.20700_gp,   82.04020_gp,   74.36270_gp,   68.34770_gp,   63.60740_gp,&
         57.25780_gp,   87.54470_gp,  157.68360_gp,  126.34700_gp,   98.71650_gp,&
         99.66250_gp,   93.74340_gp,  102.85060_gp,   83.27440_gp,   78.19760_gp,&
         72.99650_gp,   70.45090_gp,   70.44200_gp,  104.89720_gp,   97.28260_gp,&
         92.82620_gp,   89.05820_gp,   82.77130_gp,  105.37260_gp,  203.23960_gp,&
         159.49680_gp,  108.11590_gp,  105.91500_gp,  102.71500_gp,  102.83980_gp,&
         98.06280_gp,  103.26480_gp,   97.35910_gp,   98.44450_gp,   92.69390_gp,&
         90.16270_gp,   89.47740_gp,   93.38080_gp,   86.62020_gp,  112.47170_gp,&
         106.07460_gp,   98.41990_gp,   98.27200_gp,   89.25080_gp,   84.57900_gp,&
         81.16350_gp,   77.71150_gp,   75.94820_gp,  113.75630_gp,  106.14180_gp,&
         104.27680_gp,  102.47400_gp,   97.15790_gp,  118.80190_gp,  205.07880_gp,&
         127.74120_gp,  144.66290_gp,  130.73320_gp,  117.32470_gp,  113.35810_gp/)
    vdwparams%coeffs(17201:17300)=(/  130.59750_gp,   31.24070_gp,   31.22210_gp,   23.89650_gp,   19.67680_gp,&
         5.37170_gp,    4.01790_gp,   58.67030_gp,   38.70270_gp,   28.77130_gp,&
         21.14280_gp,   15.83640_gp,   12.62770_gp,   10.02220_gp,    8.01750_gp,&
         71.09990_gp,   60.29830_gp,   58.19510_gp,   48.90660_gp,   40.51670_gp,&
         35.13770_gp,   29.99180_gp,   25.55140_gp,  116.70690_gp,  102.71260_gp,&
         86.33650_gp,   85.22190_gp,   79.00030_gp,   63.18670_gp,   69.25980_gp,&
         55.33560_gp,   59.26880_gp,   60.25750_gp,   47.05850_gp,   49.38970_gp,&
         57.40160_gp,   53.01100_gp,   47.42020_gp,   43.98190_gp,   39.88210_gp,&
         35.89610_gp,  132.48790_gp,  122.69280_gp,  111.16060_gp,  102.28360_gp,&
         94.92810_gp,   76.02390_gp,   83.64090_gp,   66.30250_gp,   71.80710_gp,&
         67.33540_gp,   56.75110_gp,   60.23340_gp,   72.85710_gp,   69.60460_gp,&
         64.42910_gp,   61.31590_gp,   57.07380_gp,   52.82690_gp,  162.20010_gp,&
         155.54370_gp,  141.09330_gp,   75.99100_gp,  139.74720_gp,  134.80190_gp,&
         131.57850_gp,  128.58460_gp,  125.93290_gp,  103.60070_gp,  111.42770_gp,&
         108.18200_gp,  114.38200_gp,  111.99780_gp,  109.90310_gp,  108.38260_gp,&
         94.14310_gp,   95.45060_gp,   89.31420_gp,   77.58950_gp,   79.43980_gp,&
         73.64270_gp,   68.74810_gp,   58.74970_gp,   55.50200_gp,   57.23280_gp,&
         76.52350_gp,   76.40130_gp,   72.43520_gp,   70.47690_gp,   66.75040_gp,&
         62.75780_gp,  157.36950_gp,  156.73170_gp,  142.66360_gp,  133.47620_gp,&
         130.79990_gp,  126.85500_gp,  127.72950_gp,  124.16310_gp,    8.08480_gp/)
    vdwparams%coeffs(17301:17400)=(/   22.89300_gp,   29.85230_gp,   24.81850_gp,   19.90900_gp,   14.86000_gp,&
         11.50090_gp,    8.53970_gp,   33.26360_gp,   50.86690_gp,   53.43650_gp,&
         45.72730_gp,   39.46170_gp,   34.75100_gp,   29.68220_gp,   48.63800_gp,&
         83.72790_gp,   50.33730_gp,   48.83960_gp,   47.88150_gp,   47.22970_gp,&
         44.48880_gp,   41.62700_gp,   39.81310_gp,   38.85970_gp,   37.70830_gp,&
         36.57640_gp,   54.59640_gp,   50.23430_gp,   46.77530_gp,   43.99310_gp,&
         40.10020_gp,   58.90150_gp,  102.01160_gp,   83.51080_gp,   66.74980_gp,&
         67.41280_gp,   63.74020_gp,   69.28380_gp,   57.19280_gp,   53.86700_gp,&
         50.47090_gp,   48.69120_gp,   48.84970_gp,   70.35110_gp,   65.96350_gp,&
         63.49720_gp,   61.35850_gp,   57.56550_gp,   71.32690_gp,  131.26330_gp,&
         105.32360_gp,   73.90320_gp,   72.41560_gp,   70.28980_gp,   70.25110_gp,&
         66.71660_gp,   70.43310_gp,   66.55590_gp,   67.10600_gp,   63.45270_gp,&
         61.75680_gp,   61.24540_gp,   63.63960_gp,   59.32240_gp,   75.44120_gp,&
         71.75240_gp,   67.16870_gp,   66.66820_gp,   61.66760_gp,   58.71620_gp,&
         56.52850_gp,   54.23520_gp,   53.25560_gp,   76.84240_gp,   72.43530_gp,&
         71.56030_gp,   70.67080_gp,   67.53190_gp,   80.43260_gp,  133.28520_gp,&
         87.21080_gp,   97.99450_gp,   89.14210_gp,   80.57240_gp,   78.01900_gp,&
         88.17680_gp,   21.12820_gp,   21.36480_gp,   16.78550_gp,   14.08070_gp,&
         10.37080_gp,   21.83700_gp,   13.84960_gp,  370.81900_gp,  208.67690_gp,&
         136.99320_gp,   90.33050_gp,   61.83720_gp,   46.03470_gp,   34.33140_gp/)
    vdwparams%coeffs(17401:17500)=(/   26.05330_gp,  442.40190_gp,  334.24930_gp,  304.73370_gp,  236.13880_gp,&
         180.98430_gp,  148.21300_gp,  119.39110_gp,   96.38130_gp,  723.63460_gp,&
         590.87370_gp,  486.46530_gp,  467.45060_gp,  426.88150_gp,  335.01180_gp,&
         366.00240_gp,  285.99690_gp,  302.45180_gp,  312.58620_gp,  238.53700_gp,&
         243.56630_gp,  289.74250_gp,  252.39640_gp,  212.33920_gp,  188.62140_gp,&
         163.22210_gp,  140.36050_gp,  808.85370_gp,  703.39060_gp,  611.69610_gp,&
         546.79350_gp,  496.48070_gp,  379.76260_gp,  425.26560_gp,  320.66780_gp,&
         351.43020_gp,  324.95510_gp,  270.16330_gp,  285.21610_gp,  361.38350_gp,&
         331.28280_gp,  292.02460_gp,  268.85950_gp,  240.89430_gp,  214.54790_gp,&
         983.95840_gp,  897.40290_gp,  783.13020_gp,  345.13820_gp,  792.90630_gp,&
         760.73520_gp,  741.48810_gp,  723.82090_gp,  708.14950_gp,  551.11930_gp,&
         626.46370_gp,  603.29280_gp,  637.81980_gp,  624.24880_gp,  611.93920_gp,&
         605.08430_gp,  506.23650_gp,  494.22900_gp,  448.19920_gp,  375.59750_gp,&
         381.47930_gp,  343.61070_gp,  312.99940_gp,  257.89210_gp,  240.17810_gp,&
         246.90090_gp,  368.75940_gp,  358.55060_gp,  326.50440_gp,  309.51770_gp,&
         283.34810_gp,  257.58240_gp,  924.25730_gp,  883.81110_gp,  776.40740_gp,&
         691.12710_gp,  688.30150_gp,  666.22690_gp,  689.65180_gp,  667.29490_gp,&
         35.91900_gp,  120.78080_gp,  153.15750_gp,  114.41900_gp,   84.61780_gp,&
         57.53980_gp,   41.06050_gp,   27.48410_gp,  176.86830_gp,  274.15900_gp,&
         275.26090_gp,  217.54710_gp,  175.50060_gp,  146.58490_gp,  118.15670_gp/)
    vdwparams%coeffs(17501:17600)=(/  240.17250_gp,  470.30590_gp,  238.15400_gp,  229.57540_gp,  224.98220_gp,&
         223.18550_gp,  204.04460_gp,  188.25250_gp,  179.45530_gp,  175.53430_gp,&
         174.36190_gp,  162.00320_gp,  271.24200_gp,  235.26700_gp,  208.52360_gp,&
         188.70250_gp,  164.43830_gp,  285.42190_gp,  571.94480_gp,  430.39400_gp,&
         314.64020_gp,  317.79230_gp,  295.18160_gp,  333.81660_gp,  255.22890_gp,&
         238.31870_gp,  220.73750_gp,  213.74770_gp,  210.54980_gp,  343.30980_gp,&
         307.54750_gp,  285.72510_gp,  268.54180_gp,  243.20970_gp,  335.24300_gp,&
         742.36730_gp,  545.07520_gp,  334.46520_gp,  327.49910_gp,  316.89690_gp,&
         318.94850_gp,  308.82190_gp,  321.30770_gp,  301.27990_gp,  307.34770_gp,&
         285.59620_gp,  277.31030_gp,  275.68110_gp,  290.74100_gp,  266.29260_gp,&
         367.10790_gp,  338.46150_gp,  306.62140_gp,  311.98680_gp,  268.78310_gp,&
         251.96650_gp,  240.10400_gp,  229.34010_gp,  220.47400_gp,  366.91360_gp,&
         331.31350_gp,  319.55040_gp,  309.34230_gp,  286.80210_gp,  377.42160_gp,&
         735.19780_gp,  396.06860_gp,  460.56730_gp,  409.74230_gp,  359.90110_gp,&
         345.71200_gp,  423.22520_gp,   97.81010_gp,   95.35750_gp,   68.29840_gp,&
         53.88030_gp,   35.96650_gp,  175.56160_gp,   26.35240_gp,   16.70310_gp,&
         456.36130_gp,  253.06490_gp,  165.59380_gp,  109.05410_gp,   74.61170_gp,&
         55.52500_gp,   41.39140_gp,   31.39440_gp,  543.93780_gp,  406.06870_gp,&
         369.38290_gp,  285.48250_gp,  218.53020_gp,  178.89260_gp,  144.06690_gp,&
         116.28170_gp,  894.62280_gp,  721.40890_gp,  592.70500_gp,  568.90770_gp/)
    vdwparams%coeffs(17601:17700)=(/  519.13280_gp,  407.72880_gp,  444.62830_gp,  347.63940_gp,  366.72360_gp,&
         379.22800_gp,  289.68240_gp,  294.90170_gp,  351.11920_gp,  305.16430_gp,&
         256.42310_gp,  227.69360_gp,  196.97430_gp,  169.35470_gp,  999.74000_gp,&
         859.56780_gp,  745.18010_gp,  665.03860_gp,  603.35750_gp,  461.05860_gp,&
         516.47660_gp,  388.99560_gp,  425.98370_gp,  393.67760_gp,  327.74720_gp,&
         345.29970_gp,  438.24740_gp,  400.78940_gp,  352.78490_gp,  324.62780_gp,&
         290.73780_gp,  258.86870_gp, 1217.31410_gp, 1098.74850_gp,  955.54880_gp,&
         417.22080_gp,  969.73660_gp,  929.60420_gp,  905.88600_gp,  884.13220_gp,&
         864.83220_gp,  671.29000_gp,  768.13130_gp,  739.45380_gp,  778.13230_gp,&
         761.46230_gp,  746.32100_gp,  738.02330_gp,  616.47150_gp,  600.03080_gp,&
         543.44620_gp,  455.33130_gp,  462.11940_gp,  415.85260_gp,  378.55240_gp,&
         311.84500_gp,  290.37770_gp,  298.30060_gp,  447.39490_gp,  433.98370_gp,&
         394.57020_gp,  373.81510_gp,  342.02110_gp,  310.81620_gp, 1139.70200_gp,&
         1080.39730_gp,  945.98830_gp,  839.76740_gp,  837.81840_gp,  810.85380_gp,&
         840.69890_gp,  813.10950_gp,   43.36610_gp,  146.44260_gp,  185.73100_gp,&
         138.28380_gp,  102.20950_gp,   69.45010_gp,   49.53190_gp,   33.12650_gp,&
         214.67960_gp,  333.05050_gp,  333.53410_gp,  262.97700_gp,  211.92050_gp,&
         176.93650_gp,  142.58060_gp,  291.60280_gp,  575.47210_gp,  288.08340_gp,&
         277.80740_gp,  272.24960_gp,  270.16230_gp,  246.55880_gp,  227.43270_gp,&
         216.83930_gp,  212.11940_gp,  210.90190_gp,  195.59730_gp,  328.51000_gp/)
    vdwparams%coeffs(17701:17800)=(/  284.37390_gp,  251.81270_gp,  227.79820_gp,  198.44730_gp,  346.49590_gp,&
         700.67460_gp,  524.22670_gp,  381.03090_gp,  384.76780_gp,  357.35760_gp,&
         405.20110_gp,  308.54590_gp,  288.12910_gp,  266.84840_gp,  258.40650_gp,&
         254.22360_gp,  416.00330_gp,  371.87790_gp,  345.13030_gp,  324.23490_gp,&
         293.53370_gp,  406.37240_gp,  911.83210_gp,  664.59040_gp,  404.29100_gp,&
         395.85770_gp,  383.02380_gp,  385.63610_gp,  373.69250_gp,  388.45580_gp,&
         364.15710_gp,  371.74820_gp,  345.11820_gp,  335.07830_gp,  333.12390_gp,&
         351.37350_gp,  321.73170_gp,  445.38150_gp,  410.12060_gp,  371.08910_gp,&
         377.97370_gp,  324.78730_gp,  304.40310_gp,  290.03910_gp,  277.07890_gp,&
         266.11650_gp,  444.75260_gp,  400.65540_gp,  386.05640_gp,  373.56510_gp,&
         346.18570_gp,  457.42940_gp,  900.69660_gp,  478.74300_gp,  558.21910_gp,&
         496.75560_gp,  435.13510_gp,  417.82680_gp,  513.93660_gp,  118.12330_gp,&
         115.24800_gp,   82.43030_gp,   65.01640_gp,   43.37880_gp,  212.19710_gp,&
         256.74970_gp,   25.70210_gp,   16.58060_gp,  412.43230_gp,  236.89910_gp,&
         158.01650_gp,  105.57130_gp,   73.00730_gp,   54.74430_gp,   41.07280_gp,&
         31.31070_gp,  492.65620_gp,  377.83600_gp,  347.11540_gp,  271.88440_gp,&
         210.48400_gp,  173.57980_gp,  140.78240_gp,  114.33480_gp,  805.49240_gp,&
         664.09770_gp,  548.21150_gp,  528.45930_gp,  483.45990_gp,  379.93850_gp,&
         415.61120_gp,  325.31170_gp,  344.94790_gp,  355.76410_gp,  271.89110_gp,&
         278.96960_gp,  331.28490_gp,  290.90920_gp,  246.73790_gp,  220.38060_gp/)
    vdwparams%coeffs(17801:17900)=(/  191.79920_gp,  165.82300_gp,  901.77390_gp,  790.67650_gp,  691.20290_gp,&
         620.07920_gp,  564.49010_gp,  433.96780_gp,  485.00980_gp,  367.73060_gp,&
         402.62390_gp,  372.88980_gp,  310.11860_gp,  327.96220_gp,  413.37540_gp,&
         381.22300_gp,  338.31890_gp,  312.87500_gp,  281.72560_gp,  252.13370_gp,&
         1098.12390_gp, 1007.85030_gp,  883.88150_gp,  399.58530_gp,  891.90420_gp,&
         856.28170_gp,  834.76410_gp,  814.99180_gp,  797.46090_gp,  624.73210_gp,&
         705.00060_gp,  679.65440_gp,  719.03230_gp,  703.78810_gp,  690.01840_gp,&
         682.08270_gp,  573.24330_gp,  562.78660_gp,  512.35040_gp,  431.05180_gp,&
         438.32130_gp,  396.19890_gp,  361.96280_gp,  299.35640_gp,  279.18270_gp,&
         287.24760_gp,  422.96870_gp,  412.96020_gp,  378.20780_gp,  359.84120_gp,&
         330.92490_gp,  302.14620_gp, 1035.38600_gp,  995.37320_gp,  878.25500_gp,&
         786.53170_gp,  781.31330_gp,  756.35270_gp,  779.96360_gp,  755.12290_gp,&
         41.82220_gp,  137.44580_gp,  174.93150_gp,  132.48310_gp,   98.92780_gp,&
         67.98200_gp,   48.93910_gp,   33.10750_gp,  200.80380_gp,  310.98970_gp,&
         314.18090_gp,  250.93760_gp,  204.19980_gp,  171.65370_gp,  139.31250_gp,&
         274.80560_gp,  529.89460_gp,  274.34930_gp,  264.63020_gp,  259.29250_gp,&
         256.96820_gp,  235.96760_gp,  218.06460_gp,  207.90840_gp,  203.27770_gp,&
         201.15730_gp,  188.16130_gp,  310.77040_gp,  271.74380_gp,  242.43590_gp,&
         220.46240_gp,  193.17380_gp,  327.30800_gp,  644.40170_gp,  490.21540_gp,&
         362.56720_gp,  366.09840_gp,  340.73580_gp,  383.17340_gp,  295.96310_gp/)
    vdwparams%coeffs(17901:18000)=(/  276.54270_gp,  256.41330_gp,  248.07500_gp,  245.16140_gp,  393.56170_gp,&
         354.86870_gp,  331.34490_gp,  312.59390_gp,  284.40960_gp,  386.43270_gp,&
         835.48010_gp,  620.51600_gp,  387.38880_gp,  379.34290_gp,  367.19680_gp,&
         369.22530_gp,  356.36770_gp,  371.78190_gp,  348.90330_gp,  355.35340_gp,&
         331.01430_gp,  321.51310_gp,  319.53220_gp,  336.42390_gp,  308.79300_gp,&
         421.24260_gp,  389.87130_gp,  354.62710_gp,  359.55640_gp,  312.68360_gp,&
         293.59280_gp,  280.03950_gp,  267.48530_gp,  257.98760_gp,  421.47770_gp,&
         382.92490_gp,  370.65420_gp,  359.83420_gp,  334.97230_gp,  435.20780_gp,&
         830.00970_gp,  458.62090_gp,  530.65480_gp,  473.29240_gp,  417.25880_gp,&
         401.19370_gp,  485.68000_gp,  113.28200_gp,  110.84470_gp,   80.26190_gp,&
         63.71510_gp,   42.96790_gp,  201.77580_gp,  243.84070_gp,  233.12540_gp,&
         22.74290_gp,   15.02450_gp,  337.46620_gp,  199.95760_gp,  136.24020_gp,&
         92.60800_gp,   64.92090_gp,   49.17290_gp,   37.22050_gp,   28.58050_gp,&
         403.98020_gp,  317.08240_gp,  294.34030_gp,  233.87520_gp,  183.43220_gp,&
         152.64700_gp,  124.91620_gp,  102.27010_gp,  659.23940_gp,  552.48890_gp,&
         457.96080_gp,  443.47230_gp,  406.76720_gp,  320.27340_gp,  351.02240_gp,&
         275.41380_gp,  293.16610_gp,  301.49740_gp,  230.90910_gp,  238.50630_gp,&
         282.38710_gp,  250.64010_gp,  214.83630_gp,  193.24100_gp,  169.43130_gp,&
         147.51930_gp,  739.72690_gp,  657.78770_gp,  579.51050_gp,  522.57210_gp,&
         477.47720_gp,  369.67940_gp,  412.04130_gp,  314.82970_gp,  344.29870_gp/)
    vdwparams%coeffs(18001:18100)=(/  319.60770_gp,  265.95030_gp,  281.94690_gp,  352.72470_gp,  327.94600_gp,&
         293.61680_gp,  273.09200_gp,  247.47320_gp,  222.86600_gp,  901.72930_gp,&
         836.99700_gp,  739.56740_gp,  346.46750_gp,  742.67020_gp,  713.79130_gp,&
         696.06280_gp,  679.74330_gp,  665.28190_gp,  526.18270_gp,  587.09720_gp,&
         566.87340_gp,  600.88550_gp,  588.22970_gp,  576.86430_gp,  569.98510_gp,&
         482.16270_gp,  477.07590_gp,  436.67570_gp,  369.35120_gp,  376.21900_gp,&
         341.71560_gp,  313.44980_gp,  260.59070_gp,  243.52110_gp,  250.84270_gp,&
         362.25560_gp,  355.71630_gp,  328.25370_gp,  313.78020_gp,  290.26290_gp,&
         266.50630_gp,  855.32140_gp,  830.08340_gp,  737.40950_gp,  666.21410_gp,&
         659.37380_gp,  638.46460_gp,  654.89730_gp,  634.60740_gp,   36.49360_gp,&
         116.39410_gp,  148.81940_gp,  114.79680_gp,   86.81280_gp,   60.50850_gp,&
         44.09040_gp,   30.28940_gp,  169.56980_gp,  262.16960_gp,  267.14040_gp,&
         216.36930_gp,  178.06220_gp,  150.93640_gp,  123.60110_gp,  234.49210_gp,&
         442.07020_gp,  236.29650_gp,  228.10870_gp,  223.47840_gp,  221.19860_gp,&
         204.33190_gp,  189.24970_gp,  180.48460_gp,  176.37500_gp,  173.70060_gp,&
         163.88840_gp,  265.64070_gp,  234.77850_gp,  211.23800_gp,  193.29960_gp,&
         170.58580_gp,  280.11830_gp,  537.43520_gp,  415.35780_gp,  312.30500_gp,&
         315.27630_gp,  294.24360_gp,  328.33170_gp,  257.20590_gp,  240.57510_gp,&
         223.41250_gp,  215.94190_gp,  214.29840_gp,  336.82850_gp,  306.37120_gp,&
         287.92660_gp,  272.94170_gp,  249.79920_gp,  332.88760_gp,  695.26390_gp/)
    vdwparams%coeffs(18101:18200)=(/  525.25150_gp,  336.08240_gp,  329.13430_gp,  318.75390_gp,  320.10570_gp,&
         307.71680_gp,  322.11770_gp,  302.66100_gp,  307.57360_gp,  287.45570_gp,&
         279.32320_gp,  277.49490_gp,  291.50550_gp,  268.32640_gp,  360.72920_gp,&
         335.64790_gp,  307.01360_gp,  309.82750_gp,  272.86500_gp,  256.78780_gp,&
         245.28060_gp,  234.33320_gp,  226.97090_gp,  361.72460_gp,  331.38450_gp,&
         322.26070_gp,  313.98590_gp,  293.81590_gp,  375.07170_gp,  694.05520_gp,&
         397.72110_gp,  457.02730_gp,  408.99300_gp,  362.53380_gp,  349.05070_gp,&
         416.13190_gp,   98.18590_gp,   96.53240_gp,   70.92690_gp,   56.79810_gp,&
         38.85780_gp,  173.22270_gp,  209.20610_gp,  201.44650_gp,  175.55340_gp,&
         19.47290_gp,   13.20900_gp,  265.54640_gp,  162.68750_gp,  113.43240_gp,&
         78.57290_gp,   55.91360_gp,   42.82450_gp,   32.73650_gp,   25.34350_gp,&
         318.71200_gp,  256.43150_gp,  240.71650_gp,  194.23140_gp,  154.49350_gp,&
         129.84130_gp,  107.29640_gp,   88.62680_gp,  519.38020_gp,  442.87350_gp,&
         368.71310_gp,  358.85540_gp,  330.09380_gp,  260.56470_gp,  286.05270_gp,&
         225.13570_gp,  240.51750_gp,  246.58640_gp,  189.41220_gp,  196.93310_gp,&
         232.34090_gp,  208.56850_gp,  180.80070_gp,  163.86810_gp,  144.83870_gp,&
         127.07910_gp,  584.37840_gp,  527.38380_gp,  468.52690_gp,  424.87480_gp,&
         389.78890_gp,  304.21090_gp,  338.02960_gp,  260.53350_gp,  284.46820_gp,&
         264.73450_gp,  220.55850_gp,  234.31540_gp,  290.72550_gp,  272.58810_gp,&
         246.33500_gp,  230.51370_gp,  210.31760_gp,  190.68220_gp,  713.22290_gp/)
    vdwparams%coeffs(18201:18300)=(/  669.91060_gp,  596.71860_gp,  290.45550_gp,  596.19860_gp,  573.68400_gp,&
         559.61080_gp,  546.62920_gp,  535.13180_gp,  427.73380_gp,  471.58580_gp,&
         456.10620_gp,  484.21010_gp,  474.07530_gp,  465.03260_gp,  459.26340_gp,&
         391.29390_gp,  390.32880_gp,  359.37970_gp,  305.82610_gp,  312.04680_gp,&
         284.91830_gp,  262.49830_gp,  219.53350_gp,  205.62600_gp,  212.01590_gp,&
         299.91980_gp,  296.20520_gp,  275.49690_gp,  264.64810_gp,  246.33950_gp,&
         227.53380_gp,  681.02410_gp,  667.41950_gp,  597.23200_gp,  544.74490_gp,&
         537.15360_gp,  520.27520_gp,  530.64210_gp,  514.69140_gp,   30.77470_gp,&
         95.05310_gp,  122.15340_gp,   96.11000_gp,   73.70060_gp,   52.17530_gp,&
         38.52760_gp,   26.91630_gp,  138.11660_gp,  213.10360_gp,  219.12850_gp,&
         180.16220_gp,  150.07480_gp,  128.37640_gp,  106.16060_gp,  193.34120_gp,&
         355.61610_gp,  196.66010_gp,  190.02900_gp,  186.15760_gp,  184.02800_gp,&
         171.03060_gp,  158.79330_gp,  151.49830_gp,  147.97530_gp,  145.02500_gp,&
         138.03960_gp,  219.21940_gp,  195.95580_gp,  177.90880_gp,  163.90830_gp,&
         145.77240_gp,  231.72650_gp,  432.27830_gp,  339.80820_gp,  259.99910_gp,&
         262.42950_gp,  245.67190_gp,  271.94340_gp,  216.20660_gp,  202.48030_gp,&
         188.37520_gp,  181.91000_gp,  181.25820_gp,  278.50870_gp,  255.63360_gp,&
         241.88440_gp,  230.47020_gp,  212.26280_gp,  277.21210_gp,  557.99180_gp,&
         429.30150_gp,  281.92380_gp,  276.12630_gp,  267.56270_gp,  268.33820_gp,&
         256.90680_gp,  269.82370_gp,  253.86320_gp,  257.39270_gp,  241.38000_gp/)
    vdwparams%coeffs(18301:18400)=(/  234.65510_gp,  233.02060_gp,  244.17040_gp,  225.45330_gp,  298.46420_gp,&
         279.31410_gp,  257.02350_gp,  258.11720_gp,  230.37700_gp,  217.35740_gp,&
         207.95360_gp,  198.75510_gp,  193.32900_gp,  300.09690_gp,  277.30190_gp,&
         270.95850_gp,  265.00160_gp,  249.34960_gp,  312.47270_gp,  559.92810_gp,&
         333.46380_gp,  380.50360_gp,  341.79510_gp,  304.67440_gp,  293.76970_gp,&
         344.71170_gp,   82.19730_gp,   81.26930_gp,   60.66900_gp,   49.06290_gp,&
         34.10420_gp,  143.55230_gp,  173.27790_gp,  168.13020_gp,  147.87500_gp,&
         125.81090_gp,   30.67830_gp,   19.72470_gp,  525.31960_gp,  291.70210_gp,&
         191.41350_gp,  126.68360_gp,   87.19060_gp,   65.25250_gp,   48.93310_gp,&
         37.32580_gp,  626.55260_gp,  468.15270_gp,  426.05680_gp,  329.81710_gp,&
         253.15190_gp,  207.85310_gp,  168.00650_gp,  136.16380_gp, 1030.53940_gp,&
         832.23790_gp,  683.87860_gp,  656.74740_gp,  599.45210_gp,  471.38320_gp,&
         513.65160_gp,  402.18720_gp,  423.95230_gp,  438.24820_gp,  335.35920_gp,&
         341.19230_gp,  405.64880_gp,  352.93360_gp,  297.14260_gp,  264.37900_gp,&
         229.31470_gp,  197.76450_gp, 1152.01970_gp,  991.91430_gp,  860.32290_gp,&
         768.20000_gp,  697.31480_gp,  533.79780_gp,  597.49150_gp,  450.95970_gp,&
         493.39260_gp,  456.20630_gp,  380.42830_gp,  400.45350_gp,  507.26050_gp,&
         464.11930_gp,  409.01730_gp,  376.80930_gp,  338.04940_gp,  301.62730_gp,&
         1402.12900_gp, 1267.84900_gp, 1103.15560_gp,  483.97150_gp, 1119.48420_gp,&
         1073.39520_gp, 1046.03880_gp, 1020.93510_gp,  998.65850_gp,  776.06950_gp/)
    vdwparams%coeffs(18401:18500)=(/  887.04210_gp,  853.81990_gp,  898.65150_gp,  879.38870_gp,  861.89750_gp,&
         852.25470_gp,  712.26080_gp,  693.44710_gp,  628.55660_gp,  527.45270_gp,&
         535.33140_gp,  482.22940_gp,  439.41270_gp,  362.78720_gp,  338.12750_gp,&
         347.21000_gp,  518.71070_gp,  503.21870_gp,  457.90780_gp,  434.14870_gp,&
         397.75460_gp,  362.06750_gp, 1313.98370_gp, 1247.18710_gp, 1092.56780_gp,&
         970.93620_gp,  968.77080_gp,  937.74220_gp,  972.03520_gp,  940.20330_gp,&
         50.23680_gp,  168.90090_gp,  214.31710_gp,  160.07340_gp,  118.82550_gp,&
         81.25250_gp,   58.34230_gp,   39.40480_gp,  247.82230_gp,  384.13080_gp,&
         384.85550_gp,  303.97250_gp,  245.56540_gp,  205.60060_gp,  166.28780_gp,&
         337.77950_gp,  664.06200_gp,  333.67800_gp,  321.84410_gp,  315.46860_gp,&
         313.05380_gp,  285.80680_gp,  263.79300_gp,  251.58410_gp,  246.11560_gp,&
         244.67160_gp,  227.01180_gp,  379.66500_gp,  329.05540_gp,  291.86010_gp,&
         264.50720_gp,  231.01070_gp,  401.65720_gp,  808.69090_gp,  605.84530_gp,&
         441.49510_gp,  445.98950_gp,  414.48110_gp,  469.57520_gp,  358.37320_gp,&
         334.88270_gp,  310.39090_gp,  300.60940_gp,  295.69410_gp,  481.67290_gp,&
         430.88530_gp,  400.23370_gp,  376.37060_gp,  341.27320_gp,  471.05090_gp,&
         1051.94610_gp,  767.94120_gp,  469.04980_gp,  459.29190_gp,  444.46290_gp,&
         447.38260_gp,  433.52520_gp,  450.52800_gp,  422.51560_gp,  431.17490_gp,&
         400.47540_gp,  388.85070_gp,  386.53820_gp,  407.39760_gp,  373.32130_gp,&
         515.54770_gp,  475.26240_gp,  430.58420_gp,  438.32070_gp,  377.53700_gp/)
    vdwparams%coeffs(18501:18600)=(/  354.20050_gp,  337.73330_gp,  322.85240_gp,  310.20900_gp,  515.64400_gp,&
         464.87620_gp,  448.12230_gp,  433.87110_gp,  402.55590_gp,  530.18890_gp,&
         1039.80070_gp,  555.26920_gp,  646.89110_gp,  576.09270_gp,  505.34480_gp,&
         485.42450_gp,  595.84710_gp,  136.59980_gp,  133.59950_gp,   96.05530_gp,&
         76.13650_gp,   51.24250_gp,  245.09830_gp,  296.55720_gp,  281.98140_gp,&
         242.44540_gp,  201.38950_gp,  343.33430_gp,   29.86150_gp,   19.43040_gp,&
         481.95860_gp,  275.00610_gp,  183.24130_gp,  122.61150_gp,   85.04640_gp,&
         63.98390_gp,   48.18920_gp,   36.88000_gp,  575.74000_gp,  439.20300_gp,&
         402.93570_gp,  315.26550_gp,  244.12340_gp,  201.54660_gp,  163.74820_gp,&
         133.28650_gp,  943.39620_gp,  774.13680_gp,  638.48690_gp,  615.24940_gp,&
         562.71800_gp,  442.67360_gp,  483.59930_gp,  378.93660_gp,  401.12950_gp,&
         413.78410_gp,  316.71020_gp,  324.28910_gp,  384.81720_gp,  337.57610_gp,&
         286.30450_gp,  255.87020_gp,  222.92560_gp,  193.02240_gp, 1056.14440_gp,&
         922.14300_gp,  805.04650_gp,  721.78370_gp,  656.95460_gp,  505.23460_gp,&
         564.55600_gp,  428.26290_gp,  468.53580_gp,  433.94840_gp,  361.50570_gp,&
         381.72980_gp,  481.00660_gp,  443.06050_gp,  393.01150_gp,  363.47770_gp,&
         327.43400_gp,  293.27410_gp, 1285.99570_gp, 1176.26020_gp, 1030.09750_gp,&
         464.46000_gp, 1040.81540_gp,  998.99710_gp,  973.82280_gp,  950.69090_gp,&
         930.17500_gp,  728.11220_gp,  823.68220_gp,  793.85410_gp,  838.37100_gp,&
         820.54230_gp,  804.42380_gp,  795.18560_gp,  667.89240_gp,  654.62000_gp/)
    vdwparams%coeffs(18601:18700)=(/  595.77000_gp,  501.53070_gp,  509.81600_gp,  460.83870_gp,  421.09420_gp,&
         348.66710_gp,  325.33860_gp,  334.53100_gp,  492.68870_gp,  480.44970_gp,&
         439.74760_gp,  418.33640_gp,  384.79250_gp,  351.51420_gp, 1211.23840_gp,&
         1160.98370_gp, 1023.14220_gp,  915.52140_gp,  910.36450_gp,  881.34280_gp,&
         909.64690_gp,  880.54330_gp,   48.49290_gp,  159.57150_gp,  203.06780_gp,&
         153.70500_gp,  114.97100_gp,   79.25860_gp,   57.28050_gp,   38.98990_gp,&
         233.46740_gp,  361.45020_gp,  364.66680_gp,  291.00370_gp,  236.87400_gp,&
         199.33220_gp,  162.05420_gp,  320.00430_gp,  618.16760_gp,  318.80040_gp,&
         307.59130_gp,  301.44200_gp,  298.81900_gp,  274.17030_gp,  253.41680_gp,&
         241.68330_gp,  236.32960_gp,  234.03740_gp,  218.66030_gp,  360.91960_gp,&
         315.33670_gp,  281.33310_gp,  255.97460_gp,  224.52090_gp,  381.21510_gp,&
         752.16410_gp,  570.86180_gp,  421.53880_gp,  425.72990_gp,  396.34700_gp,&
         446.14920_gp,  344.29480_gp,  321.85830_gp,  298.56920_gp,  288.94080_gp,&
         285.28140_gp,  457.81840_gp,  412.39390_gp,  384.90810_gp,  363.14460_gp,&
         330.53650_gp,  449.47740_gp,  975.96650_gp,  722.80680_gp,  450.30760_gp,&
         440.96740_gp,  426.86830_gp,  429.25270_gp,  414.58930_gp,  432.14410_gp,&
         405.60400_gp,  413.17950_gp,  384.77440_gp,  373.72510_gp,  371.41000_gp,&
         390.92160_gp,  358.88950_gp,  489.98500_gp,  453.45540_gp,  412.48260_gp,&
         418.38190_gp,  363.75220_gp,  341.71650_gp,  326.08130_gp,  311.64840_gp,&
         300.46490_gp,  490.70860_gp,  445.39670_gp,  430.91620_gp,  418.29040_gp/)
    vdwparams%coeffs(18701:18800)=(/  389.46080_gp,  506.12520_gp,  968.71930_gp,  533.00150_gp,  617.30540_gp,&
         550.82530_gp,  485.40710_gp,  466.74110_gp,  565.87940_gp,  131.29800_gp,&
         128.67160_gp,   93.36030_gp,   74.32970_gp,   50.39200_gp,  234.00680_gp,&
         282.86280_gp,  270.35480_gp,  233.72820_gp,  195.26710_gp,  327.64580_gp,&
         313.98460_gp,   30.25930_gp,   19.72050_gp,  486.20980_gp,  277.93710_gp,&
         185.40740_gp,  124.18490_gp,   86.21090_gp,   64.90370_gp,   48.91310_gp,&
         37.45540_gp,  580.90670_gp,  443.75570_gp,  407.33310_gp,  318.95500_gp,&
         247.16030_gp,  204.16200_gp,  165.96320_gp,  135.15860_gp,  951.62680_gp,&
         781.76790_gp,  644.93880_gp,  621.62280_gp,  568.63160_gp,  447.37630_gp,&
         488.78880_gp,  383.06060_gp,  405.58090_gp,  418.31130_gp,  320.21980_gp,&
         328.00010_gp,  389.13250_gp,  341.56270_gp,  289.85600_gp,  259.14890_gp,&
         225.88180_gp,  195.66660_gp, 1065.48410_gp,  931.21720_gp,  813.33020_gp,&
         729.42430_gp,  664.04520_gp,  510.89670_gp,  570.79220_gp,  433.19140_gp,&
         473.89730_gp,  438.97520_gp,  365.71240_gp,  386.22410_gp,  486.44970_gp,&
         448.27700_gp,  397.83330_gp,  368.05570_gp,  331.68000_gp,  297.18710_gp,&
         1297.37000_gp, 1187.67890_gp, 1040.55700_gp,  470.13410_gp, 1051.10010_gp,&
         1008.94830_gp,  983.54270_gp,  960.19540_gp,  939.48900_gp,  735.79720_gp,&
         831.77260_gp,  801.71120_gp,  846.85680_gp,  828.85550_gp,  812.58600_gp,&
         803.23470_gp,  674.88900_gp,  661.76650_gp,  602.45920_gp,  507.31550_gp,&
         515.74640_gp,  466.33110_gp,  426.21440_gp,  353.02080_gp,  329.44450_gp/)
    vdwparams%coeffs(18801:18900)=(/  338.77410_gp,  498.37090_gp,  486.15110_gp,  445.15370_gp,  423.59010_gp,&
         389.75590_gp,  356.16530_gp, 1222.42790_gp, 1172.53760_gp, 1033.75200_gp,&
         925.48980_gp,  920.08400_gp,  890.76990_gp,  919.11090_gp,  889.75190_gp,&
         49.09830_gp,  161.30130_gp,  205.31340_gp,  155.56610_gp,  116.44710_gp,&
         80.34640_gp,   58.11270_gp,   39.60150_gp,  235.97180_gp,  365.27360_gp,&
         368.70140_gp,  294.44790_gp,  239.82940_gp,  201.91810_gp,  164.24590_gp,&
         323.61900_gp,  624.30360_gp,  322.57050_gp,  311.24080_gp,  305.02010_gp,&
         302.34740_gp,  277.50050_gp,  256.52910_gp,  244.65710_gp,  239.23230_gp,&
         236.85320_gp,  221.39430_gp,  365.02390_gp,  319.10990_gp,  284.83450_gp,&
         259.25390_gp,  227.49370_gp,  385.58300_gp,  759.59780_gp,  577.02860_gp,&
         426.51570_gp,  430.76130_gp,  401.09190_gp,  451.28660_gp,  348.54900_gp,&
         325.85820_gp,  302.31120_gp,  292.55140_gp,  288.91110_gp,  463.07770_gp,&
         417.33250_gp,  389.65780_gp,  367.72520_gp,  334.81970_gp,  454.77990_gp,&
         985.42930_gp,  730.55840_gp,  455.82340_gp,  446.37170_gp,  432.11200_gp,&
         434.49190_gp,  419.56220_gp,  437.40200_gp,  410.57030_gp,  418.18290_gp,&
         389.50840_gp,  378.33210_gp,  375.97960_gp,  395.67790_gp,  363.31600_gp,&
         495.59850_gp,  458.79570_gp,  417.47930_gp,  423.34010_gp,  368.33170_gp,&
         346.06910_gp,  330.26610_gp,  315.65750_gp,  304.40290_gp,  496.42620_gp,&
         450.79950_gp,  436.25560_gp,  423.55850_gp,  394.48330_gp,  512.10760_gp,&
         978.41480_gp,  539.51560_gp,  624.59130_gp,  557.42270_gp,  491.40260_gp/)
    vdwparams%coeffs(18901:19000)=(/  472.54630_gp,  572.40170_gp,  132.88810_gp,  130.26620_gp,   94.60130_gp,&
         75.36090_gp,   51.14100_gp,  236.72460_gp,  286.13170_gp,  273.58870_gp,&
         236.63740_gp,  197.80370_gp,  331.48510_gp,  317.75890_gp,  321.58900_gp,&
         27.68800_gp,   18.11270_gp,  448.29230_gp,  254.37670_gp,  169.52210_gp,&
         113.61290_gp,   78.97610_gp,   59.54450_gp,   44.94930_gp,   34.47890_gp,&
         535.44830_gp,  406.56420_gp,  372.80190_gp,  291.63180_gp,  225.98490_gp,&
         186.76650_gp,  151.94160_gp,  123.86350_gp,  879.74450_gp,  718.16300_gp,&
         591.86950_gp,  570.22820_gp,  521.46020_gp,  410.56610_gp,  448.06980_gp,&
         351.39390_gp,  371.51490_gp,  383.25570_gp,  293.68160_gp,  300.31300_gp,&
         356.24690_gp,  312.41300_gp,  265.07780_gp,  237.05770_gp,  206.72770_gp,&
         179.19650_gp,  984.95040_gp,  855.90050_gp,  746.43350_gp,  668.96670_gp,&
         608.84420_gp,  468.40590_gp,  523.32230_gp,  397.15580_gp,  434.21760_gp,&
         402.17500_gp,  335.44360_gp,  353.82270_gp,  445.77660_gp,  410.35180_gp,&
         364.01890_gp,  336.77190_gp,  303.54590_gp,  272.07950_gp, 1199.75810_gp,&
         1092.66890_gp,  955.71690_gp,  430.33410_gp,  966.66190_gp,  927.52540_gp,&
         904.07530_gp,  882.53130_gp,  863.42100_gp,  675.51500_gp,  766.07760_gp,&
         738.24890_gp,  777.89910_gp,  761.30400_gp,  746.29430_gp,  737.72690_gp,&
         619.45150_gp,  606.44790_gp,  551.85300_gp,  464.84330_gp,  472.40710_gp,&
         427.06310_gp,  390.30550_gp,  323.44940_gp,  301.91160_gp,  310.33530_gp,&
         457.00230_gp,  445.28620_gp,  407.50360_gp,  387.71530_gp,  356.76420_gp/)
    vdwparams%coeffs(19001:19100)=(/  326.09510_gp, 1128.67500_gp, 1077.93500_gp,  948.92310_gp,  848.62350_gp,&
         844.46350_gp,  817.54620_gp,  844.21690_gp,  817.09710_gp,   44.88420_gp,&
         147.64100_gp,  187.95310_gp,  142.27280_gp,  106.57110_gp,   73.62750_gp,&
         53.34000_gp,   36.45130_gp,  216.18670_gp,  334.67830_gp,  337.41280_gp,&
         269.24130_gp,  219.30330_gp,  184.72410_gp,  150.37570_gp,  296.74310_gp,&
         574.23510_gp,  295.23740_gp,  284.94970_gp,  279.27550_gp,  276.88190_gp,&
         253.93700_gp,  234.76650_gp,  223.94370_gp,  218.99430_gp,  216.92030_gp,&
         202.59800_gp,  334.10780_gp,  291.87130_gp,  260.49700_gp,  237.15920_gp,&
         208.20220_gp,  353.59290_gp,  699.13020_gp,  529.68240_gp,  390.60560_gp,&
         394.48540_gp,  367.37500_gp,  413.82300_gp,  319.15840_gp,  298.45770_gp,&
         276.94720_gp,  268.03530_gp,  264.53140_gp,  424.22720_gp,  381.98180_gp,&
         356.53490_gp,  336.46710_gp,  306.41230_gp,  416.67550_gp,  908.13400_gp,&
         670.92080_gp,  417.23950_gp,  408.58800_gp,  395.54120_gp,  397.76470_gp,&
         384.25760_gp,  400.37930_gp,  375.81950_gp,  382.88660_gp,  356.50840_gp,&
         346.27180_gp,  344.11550_gp,  362.09870_gp,  332.49970_gp,  454.19310_gp,&
         420.31920_gp,  382.36280_gp,  387.88720_gp,  337.26140_gp,  316.93610_gp,&
         302.51590_gp,  289.22160_gp,  278.81840_gp,  455.00920_gp,  412.79290_gp,&
         399.32210_gp,  387.66700_gp,  361.07620_gp,  469.15130_gp,  900.55920_gp,&
         493.79900_gp,  572.34950_gp,  510.99770_gp,  450.00790_gp,  432.70840_gp,&
         525.12640_gp,  121.46010_gp,  119.17710_gp,   86.60600_gp,   69.08090_gp/)
    vdwparams%coeffs(19101:19200)=(/   46.98750_gp,  216.43950_gp,  261.72830_gp,  250.15910_gp,  216.41010_gp,&
         180.97100_gp,  303.39700_gp,  290.73090_gp,  294.24170_gp,  269.33490_gp,&
         21.31340_gp,   14.40510_gp,  314.28070_gp,  184.90790_gp,  126.37360_gp,&
         86.54640_gp,   61.24710_gp,   46.81380_gp,   35.78110_gp,   27.73720_gp,&
         376.49060_gp,  293.67770_gp,  272.50150_gp,  216.80010_gp,  170.67670_gp,&
         142.67950_gp,  117.42580_gp,   96.75930_gp,  617.05660_gp,  513.77680_gp,&
         425.47450_gp,  412.13870_gp,  378.06600_gp,  298.48940_gp,  326.35580_gp,&
         256.82780_gp,  272.63010_gp,  280.29770_gp,  215.50230_gp,  221.96560_gp,&
         262.15300_gp,  232.79920_gp,  200.04860_gp,  180.47050_gp,  158.87040_gp,&
         138.97880_gp,  692.73590_gp,  612.35490_gp,  538.90290_gp,  485.94470_gp,&
         444.23280_gp,  344.81050_gp,  383.91860_gp,  294.22580_gp,  321.13840_gp,&
         298.30140_gp,  249.16130_gp,  263.45110_gp,  328.78360_gp,  305.51080_gp,&
         273.85520_gp,  255.10770_gp,  231.74840_gp,  209.36220_gp,  844.46350_gp,&
         780.08380_gp,  688.34330_gp,  323.46600_gp,  692.46670_gp,  665.35710_gp,&
         648.77270_gp,  633.50050_gp,  619.95990_gp,  490.59980_gp,  548.84660_gp,&
         529.80690_gp,  559.70160_gp,  547.85030_gp,  537.20110_gp,  530.75810_gp,&
         449.08170_gp,  443.57620_gp,  406.28410_gp,  344.53870_gp,  350.83400_gp,&
         319.05950_gp,  293.07050_gp,  244.55970_gp,  228.90260_gp,  235.55590_gp,&
         338.74020_gp,  332.22990_gp,  306.74920_gp,  293.47990_gp,  271.98150_gp,&
         250.33250_gp,  800.27390_gp,  773.36630_gp,  686.35790_gp,  620.32700_gp/)
    vdwparams%coeffs(19201:19300)=(/  614.75560_gp,  595.38630_gp,  611.08390_gp,  592.08300_gp,   33.94270_gp,&
         107.75270_gp,  137.89190_gp,  106.72440_gp,   81.23350_gp,   57.17660_gp,&
         42.10460_gp,   29.39230_gp,  157.38300_gp,  243.00490_gp,  247.43720_gp,&
         200.73440_gp,  165.76280_gp,  141.10920_gp,  116.21000_gp,  218.90440_gp,&
         412.06100_gp,  220.10130_gp,  212.64500_gp,  208.40940_gp,  206.34260_gp,&
         190.55280_gp,  176.66070_gp,  168.59690_gp,  164.78470_gp,  162.35450_gp,&
         153.12950_gp,  246.67960_gp,  218.22070_gp,  196.76310_gp,  180.53650_gp,&
         159.93620_gp,  261.78690_gp,  501.51140_gp,  387.09820_gp,  291.22970_gp,&
         294.11860_gp,  274.84220_gp,  306.75850_gp,  240.66050_gp,  225.38880_gp,&
         209.59760_gp,  202.66840_gp,  200.93810_gp,  314.00260_gp,  285.61070_gp,&
         268.62910_gp,  254.97950_gp,  233.89400_gp,  310.70160_gp,  649.56240_gp,&
         489.69530_gp,  313.84690_gp,  307.38260_gp,  297.75220_gp,  298.95830_gp,&
         287.52360_gp,  300.66880_gp,  282.67060_gp,  287.21540_gp,  268.48670_gp,&
         260.90960_gp,  259.15780_gp,  271.91140_gp,  250.57300_gp,  336.24450_gp,&
         313.21870_gp,  286.91560_gp,  289.46620_gp,  255.57530_gp,  240.90600_gp,&
         230.39700_gp,  220.40320_gp,  213.52200_gp,  338.05130_gp,  309.69410_gp,&
         301.20530_gp,  293.66450_gp,  275.25520_gp,  349.98930_gp,  647.93370_gp,&
         371.21200_gp,  426.72840_gp,  382.54780_gp,  339.23950_gp,  326.76380_gp,&
         389.36500_gp,   91.09850_gp,   89.96100_gp,   66.61310_gp,   53.77260_gp,&
         37.30240_gp,  160.55540_gp,  193.98150_gp,  187.00360_gp,  163.47330_gp/)
    vdwparams%coeffs(19301:19400)=(/  138.29440_gp,  225.68810_gp,  217.64940_gp,  220.41700_gp,  201.86410_gp,&
         153.34350_gp,   18.46940_gp,   12.69850_gp,  256.56890_gp,  154.76870_gp,&
         107.48660_gp,   74.55110_gp,   53.28280_gp,   41.02370_gp,   31.55820_gp,&
         24.59570_gp,  307.95570_gp,  244.75170_gp,  228.91250_gp,  184.09450_gp,&
         146.31750_gp,  123.11700_gp,  101.97460_gp,   84.51070_gp,  503.86210_gp,&
         425.23820_gp,  353.32830_gp,  343.49160_gp,  315.74920_gp,  249.68190_gp,&
         273.39370_gp,  215.57380_gp,  229.50770_gp,  235.44990_gp,  181.34820_gp,&
         187.72140_gp,  221.12170_gp,  197.93940_gp,  171.40770_gp,  155.41540_gp,&
         137.54060_gp,  120.92190_gp,  566.71690_gp,  506.80360_gp,  448.77680_gp,&
         406.32670_gp,  372.52150_gp,  290.74940_gp,  323.05520_gp,  249.07370_gp,&
         271.61880_gp,  252.75950_gp,  211.23900_gp,  223.76440_gp,  277.63500_gp,&
         259.55840_gp,  234.17500_gp,  219.04270_gp,  199.89240_gp,  181.38370_gp,&
         691.30860_gp,  644.66000_gp,  572.26670_gp,  276.39550_gp,  573.57140_gp,&
         551.61300_gp,  537.99770_gp,  525.44060_gp,  514.31120_gp,  410.08630_gp,&
         454.63510_gp,  439.40520_gp,  464.97300_gp,  455.18120_gp,  446.42190_gp,&
         440.91820_gp,  374.99750_gp,  372.59650_gp,  342.69950_gp,  291.80770_gp,&
         297.52890_gp,  271.57970_gp,  250.22110_gp,  209.63300_gp,  196.52040_gp,&
         202.40200_gp,  286.86700_gp,  282.56650_gp,  262.34870_gp,  251.84810_gp,&
         234.38100_gp,  216.58230_gp,  658.35450_gp,  641.23120_gp,  572.24570_gp,&
         520.79400_gp,  514.67310_gp,  498.56630_gp,  509.62770_gp,  494.13580_gp/)
    vdwparams%coeffs(19401:19500)=(/   29.12080_gp,   90.42310_gp,  116.10320_gp,   91.10910_gp,   69.98830_gp,&
         49.76980_gp,   36.96920_gp,   26.09260_gp,  131.82570_gp,  203.21730_gp,&
         208.29020_gp,  170.75150_gp,  142.16740_gp,  121.75410_gp,  100.91430_gp,&
         184.81700_gp,  341.76740_gp,  187.13640_gp,  180.91050_gp,  177.29930_gp,&
         175.38670_gp,  162.68440_gp,  151.07680_gp,  144.21590_gp,  140.90800_gp,&
         138.35720_gp,  131.31000_gp,  208.51100_gp,  185.92460_gp,  168.67800_gp,&
         155.46520_gp,  138.42960_gp,  221.51020_gp,  415.83400_gp,  324.97510_gp,&
         247.59810_gp,  250.02660_gp,  234.13150_gp,  259.80510_gp,  206.00200_gp,&
         193.08940_gp,  179.78240_gp,  173.73810_gp,  172.76980_gp,  265.75750_gp,&
         243.29840_gp,  229.91580_gp,  218.98580_gp,  201.72290_gp,  264.11800_gp,&
         537.58040_gp,  410.76200_gp,  268.29020_gp,  262.78540_gp,  254.64690_gp,&
         255.43700_gp,  244.96120_gp,  256.77460_gp,  241.62910_gp,  245.10920_gp,&
         229.68660_gp,  223.27460_gp,  221.71170_gp,  232.23100_gp,  214.45700_gp,&
         284.61080_gp,  266.19500_gp,  244.86880_gp,  246.20160_gp,  219.42180_gp,&
         207.18460_gp,  198.36510_gp,  189.80990_gp,  184.44130_gp,  286.72590_gp,&
         264.30990_gp,  257.93290_gp,  252.12510_gp,  237.19910_gp,  297.60870_gp,&
         538.30700_gp,  317.22510_gp,  362.80170_gp,  326.06070_gp,  290.35010_gp,&
         279.96200_gp,  329.81320_gp,   77.78230_gp,   77.08120_gp,   57.68930_gp,&
         46.86810_gp,   32.84950_gp,  136.16190_gp,  164.40450_gp,  159.32860_gp,&
         140.14430_gp,  119.34860_gp,  191.63220_gp,  185.55440_gp,  187.98380_gp/)
    vdwparams%coeffs(19501:19600)=(/  172.19800_gp,  131.82880_gp,  113.84630_gp,   18.09770_gp,   12.42770_gp,&
         253.24320_gp,  152.29860_gp,  105.55810_gp,   73.10790_gp,   52.20280_gp,&
         40.17090_gp,   30.89250_gp,   24.07410_gp,  303.91240_gp,  240.99620_gp,&
         225.16630_gp,  180.83460_gp,  143.55850_gp,  120.70600_gp,   99.91190_gp,&
         82.75910_gp,  497.31050_gp,  419.06350_gp,  348.06420_gp,  338.22830_gp,&
         310.83750_gp,  245.76290_gp,  269.04940_gp,  212.11090_gp,  225.73460_gp,&
         231.64220_gp,  178.39240_gp,  184.54030_gp,  217.40400_gp,  194.41570_gp,&
         168.19670_gp,  152.41280_gp,  134.80510_gp,  118.45780_gp,  559.21770_gp,&
         499.43750_gp,  441.92580_gp,  399.92950_gp,  366.53290_gp,  285.89770_gp,&
         317.74320_gp,  244.81840_gp,  267.00480_gp,  248.42030_gp,  207.61810_gp,&
         219.87500_gp,  272.97690_gp,  255.01000_gp,  229.88770_gp,  214.92360_gp,&
         196.02880_gp,  177.79050_gp,  682.05790_gp,  635.38110_gp,  563.62450_gp,&
         271.35450_gp,  565.20430_gp,  543.51170_gp,  530.08260_gp,  517.69890_gp,&
         506.72260_gp,  403.67390_gp,  447.99700_gp,  432.92450_gp,  458.04090_gp,&
         448.38960_gp,  439.75060_gp,  434.34720_gp,  369.18400_gp,  366.53030_gp,&
         336.95440_gp,  286.78360_gp,  292.36350_gp,  266.75590_gp,  245.69670_gp,&
         205.76580_gp,  192.87140_gp,  198.62370_gp,  281.96900_gp,  277.59650_gp,&
         257.55980_gp,  247.14790_gp,  229.89230_gp,  212.33880_gp,  649.16640_gp,&
         631.74380_gp,  563.44000_gp,  512.37500_gp,  506.53490_gp,  490.67580_gp,&
         501.82620_gp,  486.53240_gp,   28.56650_gp,   88.95380_gp,  114.16130_gp/)
    vdwparams%coeffs(19601:19700)=(/   89.43540_gp,   68.62930_gp,   48.75580_gp,   36.19190_gp,   25.53140_gp,&
         129.73320_gp,  200.00410_gp,  204.82820_gp,  167.69170_gp,  139.48020_gp,&
         119.37190_gp,   98.87470_gp,  181.70400_gp,  336.70380_gp,  183.82570_gp,&
         177.70070_gp,  174.15940_gp,  172.30360_gp,  159.74020_gp,  148.31710_gp,&
         141.58150_gp,  138.34240_gp,  135.90190_gp,  128.87600_gp,  204.95230_gp,&
         182.56910_gp,  165.50840_gp,  152.46290_gp,  135.68060_gp,  217.71950_gp,&
         409.67730_gp,  319.68600_gp,  243.20860_gp,  245.60370_gp,  229.93990_gp,&
         255.33300_gp,  202.21300_gp,  189.52850_gp,  176.45130_gp,  170.54030_gp,&
         169.52460_gp,  261.22680_gp,  238.95870_gp,  225.68120_gp,  214.86070_gp,&
         197.82510_gp,  259.41930_gp,  529.71550_gp,  404.10260_gp,  263.38510_gp,&
         257.97950_gp,  249.97990_gp,  250.78300_gp,  240.59220_gp,  252.10750_gp,&
         237.21640_gp,  240.67880_gp,  225.46970_gp,  219.16750_gp,  217.64010_gp,&
         228.00720_gp,  210.50650_gp,  279.71240_gp,  261.49250_gp,  240.43090_gp,&
         241.84390_gp,  215.30790_gp,  203.26950_gp,  194.60130_gp,  186.21520_gp,&
         180.88400_gp,  281.77560_gp,  259.55420_gp,  253.18480_gp,  247.40350_gp,&
         232.65410_gp,  292.30010_gp,  530.18710_gp,  311.43140_gp,  356.39700_gp,&
         320.21770_gp,  285.02300_gp,  274.79730_gp,  324.16540_gp,   76.34940_gp,&
         75.63270_gp,   56.54290_gp,   45.91400_gp,   32.15840_gp,  133.78530_gp,&
         161.53570_gp,  156.44370_gp,  137.50430_gp,  117.01210_gp,  188.26960_gp,&
         182.20620_gp,  184.58560_gp,  169.09000_gp,  129.34440_gp,  111.64910_gp/)
    vdwparams%coeffs(19701:19800)=(/  109.50410_gp,   20.77410_gp,   13.88040_gp,  319.84480_gp,  184.89380_gp,&
         124.90870_gp,   84.75520_gp,   59.56280_gp,   45.30660_gp,   34.49350_gp,&
         26.66270_gp,  382.68210_gp,  294.60910_gp,  271.82260_gp,  214.57260_gp,&
         167.73560_gp,  139.54300_gp,  114.30910_gp,   93.80840_gp,  628.12380_gp,&
         517.87690_gp,  427.88660_gp,  413.44790_gp,  378.73270_gp,  298.71690_gp,&
         326.26030_gp,  256.42380_gp,  271.62880_gp,  279.70520_gp,  214.80930_gp,&
         220.45250_gp,  260.75480_gp,  230.21330_gp,  196.70460_gp,  176.78440_gp,&
         155.01690_gp,  135.11470_gp,  704.27300_gp,  617.24930_gp,  540.87060_gp,&
         486.33250_gp,  443.70070_gp,  343.07100_gp,  382.56460_gp,  291.96350_gp,&
         318.87320_gp,  295.83660_gp,  247.04870_gp,  260.87010_gp,  326.88520_gp,&
         302.41310_gp,  269.80100_gp,  250.56260_gp,  226.85100_gp,  204.26630_gp,&
         858.16140_gp,  787.13930_gp,  691.63680_gp,  318.80260_gp,  697.75630_gp,&
         669.98350_gp,  653.17050_gp,  637.70390_gp,  623.98650_gp,  491.20470_gp,&
         553.14370_gp,  533.53430_gp,  562.78640_gp,  550.82490_gp,  540.04150_gp,&
         533.69230_gp,  449.99220_gp,  442.51450_gp,  404.11740_gp,  341.72070_gp,&
         347.65120_gp,  315.34640_gp,  289.04310_gp,  240.54490_gp,  224.91750_gp,&
         231.32500_gp,  336.11600_gp,  328.63270_gp,  302.20800_gp,  288.41500_gp,&
         266.45950_gp,  244.53080_gp,  810.34850_gp,  778.51930_gp,  688.35160_gp,&
         619.15580_gp,  614.84900_gp,  595.39070_gp,  612.87940_gp,  593.53010_gp,&
         33.33070_gp,  107.55820_gp,  137.29980_gp,  105.20190_gp,   79.52730_gp/)
    vdwparams%coeffs(19801:19900)=(/   55.56850_gp,   40.67990_gp,   28.21220_gp,  157.35780_gp,  243.18090_gp,&
         246.45390_gp,  198.41740_gp,  162.85410_gp,  138.01590_gp,  113.13180_gp,&
         217.59670_gp,  414.83380_gp,  217.69090_gp,  210.24640_gp,  206.07570_gp,&
         204.17500_gp,  187.95640_gp,  174.05360_gp,  166.09110_gp,  162.38310_gp,&
         160.40150_gp,  150.59460_gp,  244.99140_gp,  215.47520_gp,  193.40210_gp,&
         176.85580_gp,  156.08590_gp,  259.80390_gp,  504.99050_gp,  386.37760_gp,&
         288.02430_gp,  290.89880_gp,  271.45590_gp,  304.27490_gp,  236.90460_gp,&
         221.76150_gp,  206.06490_gp,  199.36290_gp,  197.23130_gp,  311.68130_gp,&
         282.17500_gp,  264.47480_gp,  250.38960_gp,  228.96540_gp,  307.22970_gp,&
         654.99420_gp,  489.05800_gp,  309.22990_gp,  302.84470_gp,  293.28050_gp,&
         294.67450_gp,  284.01110_gp,  296.45530_gp,  278.53190_gp,  283.35200_gp,&
         264.40130_gp,  256.88240_gp,  255.21000_gp,  268.09220_gp,  246.67760_gp,&
         333.63020_gp,  309.87820_gp,  283.00050_gp,  286.24380_gp,  251.03410_gp,&
         236.34850_gp,  225.87990_gp,  216.07200_gp,  208.86990_gp,  335.05330_gp,&
         305.58290_gp,  296.47100_gp,  288.49230_gp,  269.66300_gp,  345.99510_gp,&
         651.49430_gp,  365.82470_gp,  422.16210_gp,  377.82790_gp,  334.01640_gp,&
         321.50040_gp,  386.30110_gp,   89.78910_gp,   88.44050_gp,   64.98560_gp,&
         52.23220_gp,   35.98920_gp,  159.09800_gp,  192.27550_gp,  184.62850_gp,&
         160.65780_gp,  135.24560_gp,  223.42730_gp,  214.83990_gp,  217.51780_gp,&
         199.20380_gp,  150.48630_gp,  128.96280_gp,  126.58850_gp,  148.06570_gp/)
    vdwparams%coeffs(19901:20000)=(/   19.46480_gp,   13.09790_gp,  291.04880_gp,  170.59360_gp,  116.09230_gp,&
         79.20560_gp,   55.89540_gp,   42.64620_gp,   32.55530_gp,   25.22100_gp,&
         348.57980_gp,  271.22230_gp,  251.18770_gp,  199.28550_gp,  156.44320_gp,&
         130.51020_gp,  107.19760_gp,   88.18290_gp,  570.86590_gp,  474.82250_gp,&
         393.04870_gp,  380.43770_gp,  348.85100_gp,  275.25510_gp,  300.96400_gp,&
         236.68420_gp,  251.18370_gp,  258.38880_gp,  198.52020_gp,  204.31190_gp,&
         241.36360_gp,  213.91340_gp,  183.41100_gp,  165.19610_gp,  145.18280_gp,&
         126.80940_gp,  640.57820_gp,  565.77900_gp,  497.39960_gp,  448.16000_gp,&
         409.42870_gp,  317.34720_gp,  353.55460_gp,  270.55020_gp,  295.44510_gp,&
         274.33660_gp,  229.05370_gp,  242.18240_gp,  302.63410_gp,  280.84550_gp,&
         251.31640_gp,  233.82360_gp,  212.11680_gp,  191.36210_gp,  780.55200_gp,&
         720.70440_gp,  635.35610_gp,  296.82360_gp,  639.62140_gp,  614.52410_gp,&
         599.19690_gp,  585.08560_gp,  572.57280_gp,  452.39680_gp,  506.69850_gp,&
         489.01870_gp,  516.84540_gp,  505.90420_gp,  496.06180_gp,  490.15300_gp,&
         414.29680_gp,  408.71690_gp,  374.00200_gp,  316.78470_gp,  322.51470_gp,&
         293.05060_gp,  268.98340_gp,  224.21600_gp,  209.78810_gp,  215.87920_gp,&
         311.49490_gp,  305.28260_gp,  281.48900_gp,  269.05740_gp,  249.04140_gp,&
         228.94100_gp,  739.16470_gp,  714.05300_gp,  633.28860_gp,  571.61590_gp,&
         566.72940_gp,  548.85680_gp,  563.81700_gp,  546.23160_gp,   31.09810_gp,&
         99.34540_gp,  126.98120_gp,   97.92680_gp,   74.31490_gp,   52.15320_gp/)
    vdwparams%coeffs(20001:20100)=(/   38.31940_gp,   26.70000_gp,  145.18540_gp,  224.18120_gp,  227.95650_gp,&
         184.41750_gp,  151.91520_gp,  129.07630_gp,  106.09050_gp,  201.37990_gp,&
         380.42010_gp,  202.23920_gp,  195.35160_gp,  191.47150_gp,  189.62030_gp,&
         174.95290_gp,  162.12930_gp,  154.71980_gp,  151.24090_gp,  149.14790_gp,&
         140.45520_gp,  227.00010_gp,  200.40030_gp,  180.36950_gp,  165.25840_gp,&
         146.16810_gp,  240.66290_gp,  462.89020_gp,  356.47370_gp,  267.49730_gp,&
         270.16880_gp,  252.32810_gp,  281.97070_gp,  220.72030_gp,  206.67290_gp,&
         192.14160_gp,  185.84300_gp,  184.14330_gp,  288.88240_gp,  262.36960_gp,&
         246.45550_gp,  233.68760_gp,  214.08530_gp,  285.25530_gp,  599.51160_gp,&
         450.93480_gp,  287.96490_gp,  282.02970_gp,  273.16700_gp,  274.33460_gp,&
         264.03730_gp,  275.94550_gp,  259.37430_gp,  263.63770_gp,  246.31030_gp,&
         239.34130_gp,  237.75270_gp,  249.57470_gp,  229.85230_gp,  309.15100_gp,&
         287.70370_gp,  263.27960_gp,  265.85220_gp,  234.19810_gp,  220.66120_gp,&
         210.98510_gp,  201.83220_gp,  195.39650_gp,  310.77960_gp,  284.33510_gp,&
         276.30910_gp,  269.18960_gp,  252.03630_gp,  321.30350_gp,  597.64750_gp,&
         340.62710_gp,  391.98180_gp,  351.14350_gp,  311.17290_gp,  299.66590_gp,&
         357.95210_gp,   83.60300_gp,   82.44740_gp,   60.86780_gp,   49.05110_gp,&
         33.94270_gp,  147.69120_gp,  178.39840_gp,  171.72690_gp,  149.83280_gp,&
         126.49470_gp,  207.44400_gp,  199.85020_gp,  202.37460_gp,  185.32650_gp,&
         140.47480_gp,  120.62490_gp,  118.38000_gp,  138.02230_gp,  128.78200_gp/)
    vdwparams%coeffs(20101:20200)=(/   26.95620_gp,   17.43920_gp,  441.48880_gp,  250.18890_gp,  166.15250_gp,&
         110.83840_gp,   76.66410_gp,   57.54250_gp,   43.23850_gp,   33.02280_gp,&
         527.06710_gp,  399.92630_gp,  366.33710_gp,  285.97610_gp,  220.97630_gp,&
         182.13800_gp,  147.71690_gp,  120.02340_gp,  865.10660_gp,  706.20410_gp,&
         581.91500_gp,  560.29450_gp,  512.20750_gp,  402.83090_gp,  439.88290_gp,&
         344.53520_gp,  364.43400_gp,  376.11120_gp,  287.76680_gp,  294.31290_gp,&
         349.53890_gp,  306.09340_gp,  259.16500_gp,  231.33530_gp,  201.27060_gp,&
         174.02300_gp,  968.16820_gp,  841.36100_gp,  733.35550_gp,  656.85920_gp,&
         597.47450_gp,  458.89060_gp,  513.04100_gp,  388.60430_gp,  425.20690_gp,&
         393.63700_gp,  327.88520_gp,  346.05680_gp,  436.74070_gp,  401.72510_gp,&
         355.84740_gp,  328.80760_gp,  295.88170_gp,  264.71260_gp, 1179.24920_gp,&
         1073.94580_gp,  938.92860_gp,  420.57330_gp,  949.69780_gp,  911.18510_gp,&
         888.14170_gp,  866.97870_gp,  848.20840_gp,  662.77010_gp,  752.04060_gp,&
         724.68160_gp,  764.14750_gp,  747.86110_gp,  733.12340_gp,  724.75960_gp,&
         608.08820_gp,  595.09760_gp,  541.06220_gp,  455.07210_gp,  462.43470_gp,&
         417.62900_gp,  381.32090_gp,  315.40110_gp,  294.16770_gp,  302.43400_gp,&
         447.06730_gp,  435.47920_gp,  398.09290_gp,  378.43240_gp,  347.75130_gp,&
         317.36560_gp, 1108.83590_gp, 1059.06620_gp,  931.87680_gp,  832.42240_gp,&
         828.34850_gp,  801.86390_gp,  828.36820_gp,  801.70370_gp,   43.88810_gp,&
         145.10430_gp,  184.57440_gp,  139.25460_gp,  103.91800_gp,   71.42340_gp/)
    vdwparams%coeffs(20201:20300)=(/   51.47050_gp,   34.90080_gp,  212.36180_gp,  328.96970_gp,  331.40550_gp,&
         263.87400_gp,  214.39030_gp,  180.13610_gp,  146.18730_gp,  290.57390_gp,&
         564.11470_gp,  288.99240_gp,  278.82930_gp,  273.24400_gp,  270.91680_gp,&
         248.31530_gp,  229.43180_gp,  218.79640_gp,  213.96210_gp,  212.02750_gp,&
         197.84590_gp,  327.67720_gp,  285.80490_gp,  254.63590_gp,  231.43120_gp,&
         202.72250_gp,  345.99210_gp,  686.59850_gp,  519.44370_gp,  382.19520_gp,&
         385.93850_gp,  359.15800_gp,  404.91750_gp,  311.59910_gp,  291.22610_gp,&
         270.05910_gp,  261.36540_gp,  257.90070_gp,  415.46410_gp,  373.70470_gp,&
         348.44250_gp,  328.49050_gp,  298.69480_gp,  407.57240_gp,  891.81060_gp,&
         657.96510_gp,  407.71430_gp,  399.24480_gp,  386.44530_gp,  388.70550_gp,&
         375.62490_gp,  391.36320_gp,  367.23070_gp,  374.25920_gp,  348.30790_gp,&
         338.28140_gp,  336.21030_gp,  354.02460_gp,  324.84810_gp,  444.80230_gp,&
         411.20500_gp,  373.63290_gp,  379.28200_gp,  328.99190_gp,  308.90630_gp,&
         294.67460_gp,  281.59150_gp,  271.30180_gp,  445.09120_gp,  403.39300_gp,&
         389.99300_gp,  378.35750_gp,  351.97820_gp,  458.91490_gp,  884.04110_gp,&
         482.63570_gp,  559.82510_gp,  499.34950_gp,  439.35480_gp,  422.33120_gp,&
         513.59610_gp,  118.96190_gp,  116.48300_gp,   84.26910_gp,   66.95500_gp,&
         45.23360_gp,  212.30440_gp,  256.71210_gp,  245.05850_gp,  211.55710_gp,&
         176.44990_gp,  297.10840_gp,  284.48700_gp,  287.87970_gp,  263.38000_gp,&
         196.77040_gp,  167.56290_gp,  164.55530_gp,  194.38680_gp,  180.72890_gp/)
    vdwparams%coeffs(20301:20400)=(/  257.87270_gp,   28.63590_gp,   18.61860_gp,  454.11250_gp,  261.60340_gp,&
         175.09380_gp,  117.41980_gp,   81.50190_gp,   61.30990_gp,   46.14680_gp,&
         35.28420_gp,  542.66570_gp,  417.03070_gp,  383.63050_gp,  301.13800_gp,&
         233.70440_gp,  193.13150_gp,  157.00310_gp,  127.81130_gp,  887.70750_gp,&
         732.64330_gp,  605.03430_gp,  583.62370_gp,  534.12560_gp,  420.06420_gp,&
         459.43140_gp,  359.91720_gp,  381.65880_gp,  393.45300_gp,  300.98740_gp,&
         308.95170_gp,  366.57420_gp,  322.39740_gp,  273.96550_gp,  245.06900_gp,&
         213.66510_gp,  185.07510_gp,  994.20460_gp,  872.44650_gp,  763.34100_gp,&
         685.27160_gp,  624.20610_gp,  480.55520_gp,  536.77950_gp,  407.63390_gp,&
         446.08760_gp,  413.32440_gp,  344.00650_gp,  363.75360_gp,  457.81190_gp,&
         422.63220_gp,  375.60470_gp,  347.72860_gp,  313.52800_gp,  281.00460_gp,&
         1210.84970_gp, 1112.03760_gp,  976.00660_gp,  443.63880_gp,  984.53560_gp,&
         945.30650_gp,  921.57530_gp,  899.76220_gp,  880.42150_gp,  690.68230_gp,&
         778.59640_gp,  750.73230_gp,  793.95920_gp,  777.12650_gp,  761.93210_gp,&
         753.11720_gp,  633.53080_gp,  622.43820_gp,  567.14200_gp,  477.71970_gp,&
         485.86480_gp,  439.56860_gp,  401.91180_gp,  332.86630_gp,  310.61480_gp,&
         319.58080_gp,  468.94260_gp,  458.11250_gp,  420.04290_gp,  399.96790_gp,&
         368.25120_gp,  336.64170_gp, 1142.39470_gp, 1098.81360_gp,  970.25760_gp,&
         870.02090_gp,  864.01750_gp,  836.47440_gp,  862.05870_gp,  834.68930_gp,&
         46.44090_gp,  151.87590_gp,  193.44570_gp,  146.96060_gp,  110.06510_gp/)
    vdwparams%coeffs(20401:20500)=(/   75.92870_gp,   54.86930_gp,   37.32020_gp,  221.90310_gp,  343.51060_gp,&
         347.38790_gp,  278.06710_gp,  226.76600_gp,  190.99270_gp,  155.36720_gp,&
         304.42220_gp,  585.05150_gp,  304.19020_gp,  293.48990_gp,  287.58510_gp,&
         284.97290_gp,  261.88660_gp,  242.13490_gp,  230.89560_gp,  225.74350_gp,&
         223.25740_gp,  209.07070_gp,  344.02320_gp,  301.31020_gp,  269.23080_gp,&
         245.16040_gp,  215.18160_gp,  362.80200_gp,  711.58620_gp,  542.36550_gp,&
         402.08770_gp,  406.01940_gp,  378.12310_gp,  424.77910_gp,  328.83410_gp,&
         307.37400_gp,  285.13600_gp,  275.84840_gp,  272.72280_gp,  436.05730_gp,&
         393.65950_gp,  367.94630_gp,  347.43580_gp,  316.50110_gp,  428.61700_gp,&
         922.50460_gp,  686.47720_gp,  430.15130_gp,  421.22950_gp,  407.78710_gp,&
         409.95360_gp,  395.49150_gp,  412.71440_gp,  387.42350_gp,  394.45350_gp,&
         367.62060_gp,  357.09490_gp,  354.86500_gp,  373.42540_gp,  342.96160_gp,&
         466.74760_gp,  432.39170_gp,  393.71210_gp,  398.90670_gp,  347.67760_gp,&
         326.65250_gp,  311.70830_gp,  297.81700_gp,  287.42670_gp,  467.37550_gp,&
         425.11580_gp,  411.76820_gp,  399.99560_gp,  372.74130_gp,  482.72480_gp,&
         916.93150_gp,  509.16980_gp,  588.62060_gp,  525.40620_gp,  463.59030_gp,&
         445.86270_gp,  538.55800_gp,  125.62090_gp,  123.09690_gp,   89.44480_gp,&
         71.20170_gp,   48.24760_gp,  223.42490_gp,  269.99430_gp,  258.43380_gp,&
         223.70080_gp,  187.09060_gp,  312.54240_gp,  299.90140_gp,  303.52280_gp,&
         277.60900_gp,  208.03810_gp,  177.49730_gp,  174.26540_gp,  205.21220_gp/)
    vdwparams%coeffs(20501:20600)=(/  190.97590_gp,  271.72720_gp,  286.65450_gp,   26.70420_gp,   17.70510_gp,&
         392.43650_gp,  233.51510_gp,  159.49690_gp,  108.64070_gp,   76.29920_gp,&
         57.87840_gp,   43.87620_gp,   33.73910_gp,  469.95960_gp,  370.06380_gp,&
         343.93430_gp,  273.73340_gp,  215.01030_gp,  179.11600_gp,  146.74050_gp,&
         120.26980_gp,  766.53110_gp,  644.03510_gp,  534.15750_gp,  517.56110_gp,&
         474.88980_gp,  374.00380_gp,  410.01990_gp,  321.80920_gp,  342.72150_gp,&
         352.34030_gp,  269.93300_gp,  279.03750_gp,  330.18810_gp,  293.43220_gp,&
         251.81550_gp,  226.68450_gp,  198.93080_gp,  173.35820_gp,  860.36720_gp,&
         766.74210_gp,  676.20360_gp,  610.17610_gp,  557.78500_gp,  432.24960_gp,&
         481.62290_gp,  368.36950_gp,  402.80170_gp,  374.03530_gp,  311.27990_gp,&
         330.10750_gp,  412.55410_gp,  383.94660_gp,  344.10360_gp,  320.25290_gp,&
         290.42110_gp,  261.73610_gp, 1048.81440_gp,  975.32420_gp,  862.67800_gp,&
         405.98460_gp,  865.79390_gp,  832.26700_gp,  811.63560_gp,  792.63850_gp,&
         775.80510_gp,  614.36540_gp,  684.38050_gp,  660.94210_gp,  700.89270_gp,&
         686.14780_gp,  672.91500_gp,  664.85520_gp,  562.89600_gp,  557.47290_gp,&
         510.61380_gp,  432.17460_gp,  440.31270_gp,  400.18070_gp,  367.27240_gp,&
         305.55600_gp,  285.62940_gp,  294.25470_gp,  423.90230_gp,  416.54960_gp,&
         384.73230_gp,  367.96170_gp,  340.61170_gp,  312.93820_gp,  995.68830_gp,&
         967.79800_gp,  860.60590_gp,  778.42680_gp,  770.07950_gp,  745.69540_gp,&
         764.39730_gp,  740.80980_gp,   42.77890_gp,  135.98320_gp,  173.94150_gp/)
    vdwparams%coeffs(20601:20700)=(/  134.47040_gp,  101.84690_gp,   71.12380_gp,   51.91990_gp,   35.75820_gp,&
         198.06920_gp,  306.12750_gp,  312.25010_gp,  253.31120_gp,  208.73090_gp,&
         177.10850_gp,  145.19610_gp,  274.23250_gp,  515.43230_gp,  276.66340_gp,&
         267.10550_gp,  261.68720_gp,  258.98560_gp,  239.41200_gp,  221.80340_gp,&
         211.54220_gp,  206.71720_gp,  203.48070_gp,  192.16860_gp,  310.71100_gp,&
         274.94970_gp,  247.61740_gp,  226.75190_gp,  200.27860_gp,  327.70320_gp,&
         626.56050_gp,  485.24610_gp,  365.63120_gp,  369.11250_gp,  344.61330_gp,&
         384.16440_gp,  301.49240_gp,  282.04640_gp,  261.98970_gp,  253.21800_gp,&
         251.40820_gp,  394.10730_gp,  358.83980_gp,  337.48280_gp,  320.08780_gp,&
         293.14450_gp,  389.68370_gp,  810.23570_gp,  613.51580_gp,  393.84480_gp,&
         385.70960_gp,  373.56920_gp,  375.09470_gp,  360.42550_gp,  377.42120_gp,&
         354.68510_gp,  360.34200_gp,  336.90980_gp,  327.39560_gp,  325.23710_gp,&
         341.56110_gp,  314.51130_gp,  422.01500_gp,  392.93770_gp,  359.67380_gp,&
         362.77430_gp,  319.99870_gp,  301.24260_gp,  287.80780_gp,  274.98950_gp,&
         266.48330_gp,  423.38910_gp,  388.27410_gp,  377.78300_gp,  368.22760_gp,&
         344.77670_gp,  439.08390_gp,  809.36840_gp,  466.04700_gp,  535.07930_gp,&
         479.04230_gp,  424.95090_gp,  409.22440_gp,  486.94030_gp,  115.01370_gp,&
         113.14120_gp,   83.28650_gp,   66.78180_gp,   45.78990_gp,  202.72440_gp,&
         244.78960_gp,  235.90120_gp,  205.77840_gp,  173.52340_gp,  283.80130_gp,&
         273.76840_gp,  277.19390_gp,  253.51690_gp,  191.75980_gp,  164.52480_gp/)
    vdwparams%coeffs(20701:20800)=(/  161.41780_gp,  188.37020_gp,  175.73990_gp,  247.74210_gp,  262.03540_gp,&
         241.24940_gp,   25.12680_gp,   16.91570_gp,  351.81040_gp,  213.31490_gp,&
         147.65020_gp,  101.67860_gp,   72.03280_gp,   54.99430_gp,   41.92600_gp,&
         32.38920_gp,  421.91950_gp,  336.88690_gp,  315.11470_gp,  253.02890_gp,&
         200.37470_gp,  167.88680_gp,  138.32530_gp,  113.95770_gp,  687.73700_gp,&
         583.40090_gp,  485.05940_gp,  471.34280_gp,  433.18300_gp,  341.66370_gp,&
         374.90260_gp,  294.77780_gp,  314.57060_gp,  322.82380_gp,  247.74600_gp,&
         257.05830_gp,  303.57640_gp,  271.54740_gp,  234.56310_gp,  212.09150_gp,&
         186.99940_gp,  163.69060_gp,  773.12550_gp,  694.65560_gp,  615.53170_gp,&
         557.20710_gp,  510.54640_gp,  397.46440_gp,  442.07860_gp,  339.81570_gp,&
         371.22770_gp,  345.21240_gp,  287.50250_gp,  305.24540_gp,  379.68440_gp,&
         355.06520_gp,  319.93540_gp,  298.81560_gp,  272.05830_gp,  246.14950_gp,&
         943.14870_gp,  882.80690_gp,  784.40000_gp,  377.31490_gp,  784.97740_gp,&
         755.07130_gp,  736.48110_gp,  719.34360_gp,  704.16270_gp,  560.98910_gp,&
         620.74640_gp,  600.05550_gp,  636.81140_gp,  623.45960_gp,  611.52130_gp,&
         604.02840_gp,  513.48430_gp,  510.90340_gp,  469.53520_gp,  398.80630_gp,&
         406.70940_gp,  370.75460_gp,  341.12450_gp,  284.77860_gp,  266.55850_gp,&
         274.75930_gp,  391.14380_gp,  385.62000_gp,  357.78450_gp,  343.16660_gp,&
         318.80970_gp,  293.92920_gp,  898.72340_gp,  878.26560_gp,  784.18190_gp,&
         713.16020_gp,  704.03150_gp,  681.85240_gp,  696.68950_gp,  675.55290_gp/)
    vdwparams%coeffs(20801:20900)=(/   39.89670_gp,  124.48680_gp,  159.70990_gp,  124.88450_gp,   95.35530_gp,&
         67.19360_gp,   49.42910_gp,   34.37290_gp,  181.04870_gp,  279.50350_gp,&
         286.58090_gp,  234.50660_gp,  194.60170_gp,  165.99750_gp,  136.86440_gp,&
         252.44980_gp,  467.88320_gp,  256.04740_gp,  247.34020_gp,  242.31110_gp,&
         239.63570_gp,  222.29610_gp,  206.23670_gp,  196.74030_gp,  192.19700_gp,&
         188.65600_gp,  179.07470_gp,  286.16240_gp,  254.88570_gp,  230.75580_gp,&
         212.14750_gp,  188.22630_gp,  302.25300_gp,  568.74160_gp,  444.73220_gp,&
         338.46100_gp,  341.64920_gp,  319.53400_gp,  354.57680_gp,  280.63740_gp,&
         262.72550_gp,  244.29540_gp,  235.98710_gp,  234.84700_gp,  363.36970_gp,&
         332.58250_gp,  314.02340_gp,  298.72470_gp,  274.58620_gp,  360.80530_gp,&
         734.58650_gp,  562.00060_gp,  366.16290_gp,  358.62210_gp,  347.44280_gp,&
         348.59290_gp,  334.17340_gp,  350.60320_gp,  329.73370_gp,  334.55220_gp,&
         313.41130_gp,  304.63840_gp,  302.55580_gp,  317.27860_gp,  292.67680_gp,&
         389.27450_gp,  363.65250_gp,  334.01770_gp,  335.95050_gp,  298.62210_gp,&
         281.53260_gp,  269.22650_gp,  257.29480_gp,  249.94700_gp,  391.12750_gp,&
         360.46170_gp,  351.69120_gp,  343.55060_gp,  322.70510_gp,  406.64340_gp,&
         735.95700_gp,  433.16650_gp,  495.33320_gp,  444.42440_gp,  395.50150_gp,&
         381.18140_gp,  449.46340_gp,  106.80860_gp,  105.41670_gp,   78.31770_gp,&
         63.15490_gp,   43.70220_gp,  187.15040_gp,  225.92450_gp,  218.68530_gp,&
         191.78480_gp,  162.66590_gp,  262.36130_gp,  253.92720_gp,  257.18410_gp/)
    vdwparams%coeffs(20901:21000)=(/  235.27290_gp,  179.15590_gp,  154.30230_gp,  151.32020_gp,  175.48260_gp,&
         163.98640_gp,  229.56900_gp,  243.20020_gp,  224.98180_gp,  210.52310_gp,&
         36.86060_gp,   23.84010_gp,  628.23990_gp,  348.16900_gp,  228.95120_gp,&
         151.97450_gp,  104.90040_gp,   78.69330_gp,   59.14150_gp,   45.19410_gp,&
         749.30050_gp,  558.78190_gp,  508.80610_gp,  394.34550_gp,  303.25210_gp,&
         249.42670_gp,  201.99970_gp,  164.02750_gp, 1234.93200_gp,  994.46760_gp,&
         816.92750_gp,  784.68860_gp,  716.27780_gp,  563.71140_gp,  613.84180_gp,&
         481.04470_gp,  506.72720_gp,  523.70020_gp,  401.17140_gp,  407.94990_gp,&
         484.88130_gp,  422.18130_gp,  355.93940_gp,  317.09380_gp,  275.44760_gp,&
         237.92190_gp, 1380.84110_gp, 1185.80180_gp, 1028.25360_gp,  918.23780_gp,&
         833.71350_gp,  638.77680_gp,  714.71880_gp,  539.96390_gp,  590.37440_gp,&
         545.97680_gp,  455.72670_gp,  479.37970_gp,  606.81450_gp,  555.32330_gp,&
         489.82280_gp,  451.62950_gp,  405.61580_gp,  362.35560_gp, 1681.32490_gp,&
         1516.53690_gp, 1319.03120_gp,  579.73570_gp, 1338.97890_gp, 1283.65210_gp,&
         1250.87660_gp, 1220.80060_gp, 1194.11040_gp,  928.28660_gp, 1062.15340_gp,&
         1022.37080_gp, 1074.32200_gp, 1051.24370_gp, 1030.29280_gp, 1018.73100_gp,&
         851.60220_gp,  828.98030_gp,  751.66620_gp,  631.36390_gp,  640.74450_gp,&
         577.45940_gp,  526.43720_gp,  435.10780_gp,  405.68950_gp,  416.49740_gp,&
         621.14960_gp,  602.48990_gp,  548.55160_gp,  520.38560_gp,  477.18560_gp,&
         434.81090_gp, 1574.99420_gp, 1491.75350_gp, 1306.25580_gp, 1161.07230_gp/)
    vdwparams%coeffs(21001:21100)=(/ 1158.79960_gp, 1121.70090_gp, 1162.62900_gp, 1124.50570_gp,   60.18920_gp,&
         201.69000_gp,  256.12950_gp,  191.64290_gp,  142.61130_gp,   97.80850_gp,&
         70.42850_gp,   47.74110_gp,  296.00050_gp,  458.80010_gp,  459.73700_gp,&
         363.58720_gp,  294.21340_gp,  246.72930_gp,  199.93530_gp,  404.35370_gp,&
         794.54790_gp,  399.32770_gp,  385.27560_gp,  377.65080_gp,  374.74370_gp,&
         342.18210_gp,  315.93550_gp,  301.36140_gp,  294.80040_gp,  292.97830_gp,&
         271.98180_gp,  453.91820_gp,  393.76600_gp,  349.65980_gp,  317.24980_gp,&
         277.46890_gp,  481.06050_gp,  968.07870_gp,  725.18090_gp,  528.65740_gp,&
         534.02030_gp,  496.51280_gp,  562.43130_gp,  429.53830_gp,  401.50480_gp,&
         372.26120_gp,  360.48910_gp,  354.60960_gp,  576.28960_gp,  515.76060_gp,&
         479.38630_gp,  451.12580_gp,  409.47120_gp,  564.34540_gp, 1260.16190_gp,&
         919.44900_gp,  561.90760_gp,  550.22410_gp,  532.49960_gp,  535.94970_gp,&
         519.22710_gp,  539.62270_gp,  506.14540_gp,  516.46990_gp,  479.78110_gp,&
         465.87230_gp,  463.07400_gp,  487.86280_gp,  447.25110_gp,  617.17070_gp,&
         569.19750_gp,  515.95990_gp,  525.05030_gp,  452.76340_gp,  424.95970_gp,&
         405.31970_gp,  387.54040_gp,  372.47490_gp,  617.41750_gp,  556.80680_gp,&
         536.92140_gp,  520.08520_gp,  482.93520_gp,  635.18450_gp, 1245.23720_gp,&
         665.11670_gp,  774.84830_gp,  690.52460_gp,  605.68150_gp,  581.87050_gp,&
         713.81570_gp,  163.46210_gp,  160.10980_gp,  115.41070_gp,   91.67250_gp,&
         61.91590_gp,  292.91740_gp,  354.55260_gp,  337.37570_gp,  290.47610_gp/)
    vdwparams%coeffs(21101:21200)=(/  241.69490_gp,  410.78570_gp,  392.19670_gp,  396.82160_gp,  363.31250_gp,&
         270.75430_gp,  230.12310_gp,  226.05370_gp,  267.82340_gp,  248.73790_gp,&
         355.56210_gp,  374.09750_gp,  340.06880_gp,  314.69330_gp,  491.74480_gp,&
         36.05430_gp,   23.64620_gp,  567.32020_gp,  326.91860_gp,  219.33470_gp,&
         147.60490_gp,  102.84790_gp,   77.63350_gp,   58.63590_gp,   44.97530_gp,&
         678.17260_gp,  521.14300_gp,  479.68840_gp,  377.06390_gp,  293.23240_gp,&
         242.82560_gp,  197.87960_gp,  161.50560_gp, 1110.54260_gp,  916.08550_gp,&
         756.54430_gp,  730.05840_gp,  668.27290_gp,  526.03240_gp,  575.00600_gp,&
         450.91200_gp,  477.89770_gp,  492.51980_gp,  377.23750_gp,  387.08210_gp,&
         458.90330_gp,  403.96460_gp,  343.79750_gp,  307.97330_gp,  268.98800_gp,&
         233.45810_gp, 1244.15780_gp, 1091.24080_gp,  955.01530_gp,  857.65270_gp,&
         781.53770_gp,  602.42530_gp,  672.55530_gp,  511.47220_gp,  559.34530_gp,&
         518.43740_gp,  431.97050_gp,  456.48750_gp,  573.80210_gp,  529.91860_gp,&
         471.40640_gp,  436.80220_gp,  394.31850_gp,  353.91780_gp, 1515.37550_gp,&
         1391.17640_gp, 1221.21520_gp,  556.96970_gp, 1231.98900_gp, 1182.91410_gp,&
         1153.21660_gp, 1125.91270_gp, 1101.70080_gp,  865.00110_gp,  974.97080_gp,&
         940.11780_gp,  993.51430_gp,  972.43020_gp,  953.40360_gp,  942.32490_gp,&
         793.10590_gp,  779.30040_gp,  710.47900_gp,  599.13420_gp,  609.36150_gp,&
         551.68900_gp,  504.77160_gp,  418.66960_gp,  390.91310_gp,  402.09420_gp,&
         588.45790_gp,  574.89740_gp,  527.48310_gp,  502.56880_gp,  463.16440_gp/)
    vdwparams%coeffs(21201:21300)=(/  423.89390_gp, 1430.02470_gp, 1374.93400_gp, 1214.28600_gp, 1089.59110_gp,&
         1082.19090_gp, 1047.77930_gp, 1079.59980_gp, 1045.35530_gp,   58.27550_gp,&
         189.88940_gp,  242.00800_gp,  184.28470_gp,  138.43530_gp,   95.88970_gp,&
         69.58060_gp,   47.59140_gp,  277.57010_gp,  429.50590_gp,  434.51020_gp,&
         348.31620_gp,  284.58190_gp,  240.14930_gp,  195.82670_gp,  381.77780_gp,&
         732.12430_gp,  381.47190_gp,  368.14810_gp,  360.77080_gp,  357.48370_gp,&
         328.62230_gp,  303.96440_gp,  289.91370_gp,  283.44240_gp,  280.25690_gp,&
         262.57350_gp,  430.79250_gp,  377.69300_gp,  337.90570_gp,  308.09200_gp,&
         270.88100_gp,  455.24600_gp,  890.74400_gp,  679.45910_gp,  504.43250_gp,&
         509.41810_gp,  474.68070_gp,  532.97190_gp,  413.19030_gp,  386.39150_gp,&
         358.61330_gp,  346.93600_gp,  342.99870_gp,  546.66900_gp,  493.81290_gp,&
         461.87660_gp,  436.45200_gp,  398.03820_gp,  537.90930_gp, 1154.94230_gp,&
         860.04650_gp,  540.09920_gp,  528.91450_gp,  512.08670_gp,  514.73080_gp,&
         496.50390_gp,  518.09120_gp,  486.46500_gp,  495.19130_gp,  461.64810_gp,&
         448.45310_gp,  445.61860_gp,  468.67840_gp,  430.67940_gp,  585.21160_gp,&
         542.52680_gp,  494.40990_gp,  500.71270_gp,  437.13590_gp,  410.96170_gp,&
         392.33860_gp,  374.99560_gp,  362.03330_gp,  586.48270_gp,  533.75820_gp,&
         517.18370_gp,  502.62410_gp,  468.78500_gp,  605.78760_gp, 1148.15630_gp,&
         639.19770_gp,  738.61040_gp,  659.77710_gp,  582.47220_gp,  560.32370_gp,&
         675.95320_gp,  157.42770_gp,  154.52800_gp,  112.66800_gp,   89.96260_gp/)
    vdwparams%coeffs(21301:21400)=(/   61.27960_gp,  279.68340_gp,  338.01530_gp,  323.83050_gp,  280.74710_gp,&
         235.27170_gp,  391.83600_gp,  376.16230_gp,  380.74190_gp,  348.37280_gp,&
         261.70420_gp,  223.56370_gp,  219.47140_gp,  257.91900_gp,  240.13260_gp,&
         340.65350_gp,  359.42740_gp,  328.94310_gp,  305.65080_gp,  469.26390_gp,&
         451.08450_gp,   36.54050_gp,   24.08830_gp,  562.15830_gp,  327.21920_gp,&
         220.83610_gp,  149.26580_gp,  104.33800_gp,   78.93390_gp,   59.72990_gp,&
         45.88240_gp,  672.43030_gp,  520.68680_gp,  480.73130_gp,  379.42120_gp,&
         296.08060_gp,  245.72860_gp,  200.66650_gp,  164.07610_gp, 1099.83930_gp,&
         912.52050_gp,  754.61820_gp,  729.16220_gp,  667.96950_gp,  525.95350_gp,&
         575.39200_gp,  451.41300_gp,  479.10890_gp,  493.36920_gp,  377.99640_gp,&
         388.72950_gp,  460.49590_gp,  406.61960_gp,  347.03710_gp,  311.42450_gp,&
         272.49100_gp,  236.88550_gp, 1232.92520_gp, 1086.83940_gp,  953.48210_gp,&
         857.60110_gp,  782.30880_gp,  604.14750_gp,  674.00680_gp,  513.61970_gp,&
         561.61690_gp,  520.87810_gp,  433.92780_gp,  459.01640_gp,  575.80200_gp,&
         533.07320_gp,  475.38210_gp,  441.15130_gp,  398.88950_gp,  358.56750_gp,&
         1502.01190_gp, 1384.60110_gp, 1218.36040_gp,  561.44680_gp, 1227.16180_gp,&
         1178.72220_gp, 1149.25010_gp, 1122.13930_gp, 1098.10310_gp,  864.58960_gp,&
         970.85510_gp,  936.59920_gp,  990.84400_gp,  969.87250_gp,  950.97990_gp,&
         939.81690_gp,  792.50890_gp,  780.62700_gp,  712.79020_gp,  601.88010_gp,&
         612.49070_gp,  555.26090_gp,  508.58990_gp,  422.35350_gp,  394.54210_gp/)
    vdwparams%coeffs(21401:21500)=(/  406.00450_gp,  590.97830_gp,  578.41880_gp,  531.86650_gp,  507.39330_gp,&
         468.32900_gp,  429.22150_gp, 1420.14340_gp, 1370.18490_gp, 1212.75100_gp,&
         1091.05890_gp, 1082.33090_gp, 1047.98200_gp, 1078.06580_gp, 1044.16840_gp,&
         58.86250_gp,  190.22830_gp,  242.72640_gp,  185.77620_gp,  139.98770_gp,&
         97.28810_gp,   70.78650_gp,   48.57690_gp,  277.79360_gp,  429.64620_gp,&
         435.78140_gp,  350.70760_gp,  287.38300_gp,  243.00890_gp,  198.57750_gp,&
         383.02240_gp,  729.62110_gp,  383.86450_gp,  370.51580_gp,  363.07110_gp,&
         359.62670_gp,  331.18730_gp,  306.51240_gp,  292.35140_gp,  285.78230_gp,&
         282.17240_gp,  265.04020_gp,  432.63090_gp,  380.45750_gp,  341.14820_gp,&
         311.53680_gp,  274.38370_gp,  457.07170_gp,  887.48100_gp,  680.26830_gp,&
         507.52160_gp,  512.49950_gp,  477.88700_gp,  535.31530_gp,  416.71520_gp,&
         389.76730_gp,  361.87710_gp,  349.99750_gp,  346.47920_gp,  549.05720_gp,&
         497.24240_gp,  465.92890_gp,  440.83920_gp,  402.64130_gp,  541.12590_gp,&
         1149.71630_gp,  860.75030_gp,  544.52030_gp,  533.25740_gp,  516.35730_gp,&
         518.83450_gp,  499.86370_gp,  522.15120_gp,  490.43320_gp,  498.90510_gp,&
         465.55960_gp,  452.30680_gp,  449.40350_gp,  472.39670_gp,  434.41480_gp,&
         587.78890_gp,  545.72750_gp,  498.09170_gp,  503.76490_gp,  441.35260_gp,&
         415.15170_gp,  396.46860_gp,  378.93340_gp,  366.28550_gp,  589.38990_gp,&
         537.73400_gp,  521.73760_gp,  507.54760_gp,  474.01940_gp,  609.50110_gp,&
         1144.73570_gp,  644.38060_gp,  743.02780_gp,  664.26180_gp,  587.40530_gp/)
    vdwparams%coeffs(21501:21600)=(/  565.28300_gp,  678.85160_gp,  158.74420_gp,  155.97060_gp,  114.13350_gp,&
         91.30720_gp,   62.38930_gp,  281.29500_gp,  339.86220_gp,  326.24310_gp,&
         283.45370_gp,  238.08430_gp,  394.11100_gp,  378.93950_gp,  383.59990_gp,&
         350.96590_gp,  264.33470_gp,  226.16840_gp,  221.98580_gp,  260.21080_gp,&
         242.44220_gp,  343.05200_gp,  362.24070_gp,  332.19680_gp,  309.08380_gp,&
         472.10810_gp,  454.73560_gp,  458.68960_gp,   35.31520_gp,   23.35920_gp,&
         549.61140_gp,  316.71620_gp,  213.37350_gp,  144.25580_gp,  100.94600_gp,&
         76.46580_gp,   57.94670_gp,   44.57760_gp,  657.11510_gp,  504.67530_gp,&
         465.22590_gp,  366.61690_gp,  285.99260_gp,  237.44480_gp,  194.03010_gp,&
         158.79150_gp, 1078.83360_gp,  887.61670_gp,  733.01140_gp,  707.81110_gp,&
         648.11200_gp,  510.73850_gp,  557.95140_gp,  438.06030_gp,  464.06680_gp,&
         478.03720_gp,  366.66480_gp,  376.23990_gp,  445.72790_gp,  393.03660_gp,&
         335.29210_gp,  300.92580_gp,  263.40340_gp,  229.11910_gp, 1209.20810_gp,&
         1057.87480_gp,  926.14390_gp,  832.17430_gp,  758.77890_gp,  585.81260_gp,&
         653.60280_gp,  497.92860_gp,  544.08210_gp,  504.50370_gp,  420.85520_gp,&
         444.49320_gp,  557.93400_gp,  515.73740_gp,  459.57500_gp,  426.42020_gp,&
         385.59230_gp,  346.70630_gp, 1473.79700_gp, 1349.39700_gp, 1184.66100_gp,&
         543.05060_gp, 1195.26540_gp, 1147.47570_gp, 1118.62780_gp, 1092.10250_gp,&
         1068.58050_gp,  840.00980_gp,  947.32140_gp,  913.63480_gp,  963.55250_gp,&
         943.06260_gp,  924.58410_gp,  913.77400_gp,  769.78010_gp,  756.63130_gp/)
    vdwparams%coeffs(21601:21700)=(/  690.39390_gp,  583.08880_gp,  593.08720_gp,  537.46300_gp,  492.19170_gp,&
         408.91560_gp,  382.05000_gp,  392.94050_gp,  573.02740_gp,  559.97550_gp,&
         514.44520_gp,  490.63940_gp,  452.81960_gp,  415.06050_gp, 1390.46370_gp,&
         1333.95560_gp, 1178.21430_gp, 1058.27040_gp, 1051.12310_gp, 1017.73280_gp,&
         1048.10380_gp, 1014.88280_gp,   56.84660_gp,  184.12220_gp,  234.95680_gp,&
         179.52910_gp,  135.34930_gp,   94.16310_gp,   68.60870_gp,   47.19410_gp,&
         269.17780_gp,  416.41830_gp,  421.64610_gp,  338.87850_gp,  277.61880_gp,&
         234.83210_gp,  192.01920_gp,  371.44330_gp,  710.75730_gp,  371.32440_gp,&
         358.52320_gp,  351.34710_gp,  348.10210_gp,  320.22560_gp,  296.37490_gp,&
         282.73960_gp,  276.41130_gp,  273.10650_gp,  256.21710_gp,  418.62090_gp,&
         367.71370_gp,  329.61390_gp,  301.04260_gp,  265.23400_gp,  443.26820_gp,&
         865.24620_gp,  660.73770_gp,  491.30810_gp,  496.11270_gp,  462.65460_gp,&
         519.08450_gp,  403.20320_gp,  377.22490_gp,  350.29100_gp,  338.83350_gp,&
         335.13120_gp,  531.77200_gp,  480.95640_gp,  450.41590_gp,  426.11180_gp,&
         389.21120_gp,  524.15620_gp, 1122.76720_gp,  836.55040_gp,  526.67590_gp,&
         515.77970_gp,  499.43430_gp,  501.91700_gp,  483.85420_gp,  505.06200_gp,&
         474.36130_gp,  482.73800_gp,  450.24130_gp,  437.40790_gp,  434.60070_gp,&
         456.79920_gp,  420.06140_gp,  569.59350_gp,  528.52650_gp,  482.14800_gp,&
         487.93610_gp,  426.97730_gp,  401.68310_gp,  383.65960_gp,  366.80700_gp,&
         354.37800_gp,  571.13440_gp,  520.34420_gp,  504.56970_gp,  490.75210_gp/)
    vdwparams%coeffs(21701:21800)=(/  458.30100_gp,  590.30800_gp, 1116.03770_gp,  623.20270_gp,  719.78040_gp,&
         643.70980_gp,  568.41140_gp,  546.92870_gp,  658.59550_gp,  153.29480_gp,&
         150.76980_gp,  110.36390_gp,   88.39440_gp,   60.52170_gp,  271.78760_gp,&
         328.57270_gp,  315.20730_gp,  273.86490_gp,  230.08710_gp,  381.25920_gp,&
         366.37160_gp,  370.88000_gp,  339.49100_gp,  255.76580_gp,  218.84580_gp,&
         214.80530_gp,  251.78560_gp,  234.55580_gp,  331.67180_gp,  350.08010_gp,&
         320.96670_gp,  298.68080_gp,  456.88220_gp,  439.65480_gp,  443.41820_gp,&
         428.89620_gp,   31.66810_gp,   21.35160_gp,  461.57510_gp,  273.54460_gp,&
         187.46660_gp,  128.49780_gp,   90.90430_gp,   69.41810_gp,   52.98450_gp,&
         41.00350_gp,  553.03570_gp,  433.82780_gp,  403.29140_gp,  321.51630_gp,&
         253.41230_gp,  211.90690_gp,  174.39200_gp,  143.63890_gp,  905.08480_gp,&
         757.08580_gp,  627.55600_gp,  608.30000_gp,  558.23240_gp,  440.55160_gp,&
         482.13620_gp,  379.26880_gp,  403.14840_gp,  414.31560_gp,  318.32390_gp,&
         328.47070_gp,  388.04460_gp,  345.15780_gp,  296.91650_gp,  267.96230_gp,&
         235.92640_gp,  206.36370_gp, 1016.37180_gp,  902.10970_gp,  795.17640_gp,&
         717.66320_gp,  656.37740_gp,  509.74240_gp,  567.44410_gp,  435.08880_gp,&
         475.02730_gp,  441.33190_gp,  368.28980_gp,  389.83000_gp,  486.21160_gp,&
         452.46160_gp,  406.03120_gp,  378.43120_gp,  343.90960_gp,  310.74420_gp,&
         1239.19340_gp, 1148.51060_gp, 1015.11540_gp,  479.42470_gp, 1019.80840_gp,&
         980.14620_gp,  955.78910_gp,  933.35430_gp,  913.46790_gp,  723.93420_gp/)
    vdwparams%coeffs(21801:21900)=(/  807.73280_gp,  779.96820_gp,  825.02430_gp,  807.59920_gp,  791.96020_gp,&
         782.41790_gp,  662.69770_gp,  655.75410_gp,  601.05960_gp,  509.81570_gp,&
         519.31350_gp,  472.49570_gp,  434.13160_gp,  362.21680_gp,  338.98810_gp,&
         348.97520_gp,  500.83330_gp,  491.81440_gp,  454.57810_gp,  435.14050_gp,&
         403.44890_gp,  371.42770_gp, 1175.88730_gp, 1139.55560_gp, 1012.73410_gp,&
         916.52630_gp,  907.43040_gp,  878.83040_gp,  901.04700_gp,  873.18690_gp,&
         50.41730_gp,  159.44040_gp,  204.15190_gp,  158.35630_gp,  120.59460_gp,&
         84.85820_gp,   62.43240_gp,   43.47600_gp,  232.59320_gp,  359.18880_gp,&
         366.31710_gp,  297.74360_gp,  246.10860_gp,  209.55990_gp,  172.57580_gp,&
         323.63810_gp,  607.09190_gp,  326.07870_gp,  314.99830_gp,  308.68010_gp,&
         305.52730_gp,  282.43340_gp,  261.86300_gp,  249.86950_gp,  244.18470_gp,&
         240.36730_gp,  227.04050_gp,  365.28020_gp,  323.62930_gp,  292.04750_gp,&
         268.05260_gp,  237.50190_gp,  387.12360_gp,  738.64010_gp,  571.85300_gp,&
         431.35120_gp,  435.56920_gp,  407.06590_gp,  453.71970_gp,  356.63380_gp,&
         333.93870_gp,  310.50050_gp,  300.14930_gp,  297.84780_gp,  464.56790_gp,&
         423.15700_gp,  398.33180_gp,  378.25790_gp,  347.10000_gp,  460.12630_gp,&
         956.02290_gp,  723.24370_gp,  465.18180_gp,  455.59950_gp,  441.33640_gp,&
         443.05470_gp,  425.78150_gp,  445.61250_gp,  418.95420_gp,  425.55500_gp,&
         397.99830_gp,  386.78540_gp,  384.18230_gp,  403.06860_gp,  371.50160_gp,&
         497.58720_gp,  463.78430_gp,  425.06300_gp,  428.53830_gp,  378.88850_gp/)
    vdwparams%coeffs(21901:22000)=(/  357.11560_gp,  341.49670_gp,  326.57380_gp,  316.56330_gp,  500.05880_gp,&
         458.72390_gp,  446.46550_gp,  435.46460_gp,  408.32830_gp,  518.37750_gp,&
         954.63010_gp,  550.25080_gp,  631.74190_gp,  566.37730_gp,  502.65170_gp,&
         484.21620_gp,  575.61980_gp,  135.24280_gp,  133.50900_gp,   98.89030_gp,&
         79.76720_gp,   55.24840_gp,  238.01690_gp,  287.52250_gp,  277.43120_gp,&
         242.66760_gp,  205.36950_gp,  334.32060_gp,  322.66550_gp,  326.76700_gp,&
         299.14950_gp,  227.29470_gp,  195.44560_gp,  191.73330_gp,  222.95460_gp,&
         208.15700_gp,  291.72410_gp,  308.57810_gp,  284.64110_gp,  265.99860_gp,&
         401.05410_gp,  388.07180_gp,  392.06550_gp,  379.20750_gp,  337.10750_gp,&
         28.26270_gp,   19.35560_gp,  391.42550_gp,  236.84570_gp,  164.57770_gp,&
         114.07740_gp,   81.42170_gp,   62.58860_gp,   48.05400_gp,   37.37230_gp,&
         469.76960_gp,  374.26910_gp,  350.28040_gp,  281.84800_gp,  224.00340_gp,&
         188.40440_gp,  155.94050_gp,  129.11050_gp,  767.93390_gp,  649.50000_gp,&
         539.86620_gp,  524.91730_gp,  482.56220_gp,  381.40060_gp,  417.86140_gp,&
         329.31160_gp,  350.85640_gp,  359.90170_gp,  276.99700_gp,  286.99270_gp,&
         338.26270_gp,  302.93890_gp,  262.34800_gp,  237.82350_gp,  210.38120_gp,&
         184.84770_gp,  863.76550_gp,  773.94630_gp,  685.74890_gp,  621.04380_gp,&
         569.41590_gp,  444.34710_gp,  493.74300_gp,  380.57300_gp,  415.14630_gp,&
         386.29860_gp,  322.58210_gp,  341.92370_gp,  424.31450_gp,  396.89200_gp,&
         358.15190_gp,  335.00590_gp,  305.66670_gp,  277.27710_gp, 1053.78180_gp/)
    vdwparams%coeffs(22001:22100)=(/  984.18610_gp,  874.24490_gp,  422.66280_gp,  875.59270_gp,  842.16510_gp,&
         821.40430_gp,  802.25710_gp,  785.28980_gp,  626.38240_gp,  693.66130_gp,&
         670.51140_gp,  710.07710_gp,  695.14360_gp,  681.79200_gp,  673.37970_gp,&
         572.85580_gp,  569.66440_gp,  524.02130_gp,  446.08380_gp,  454.88600_gp,&
         415.19500_gp,  382.49740_gp,  320.27240_gp,  300.15390_gp,  309.19500_gp,&
         438.22480_gp,  431.87000_gp,  401.07000_gp,  385.04170_gp,  358.31060_gp,&
         331.03010_gp, 1004.07570_gp,  979.25300_gp,  874.28430_gp,  795.90040_gp,&
         786.18960_gp,  761.55790_gp,  778.12490_gp,  754.51690_gp,   44.59630_gp,&
         138.36130_gp,  177.68370_gp,  139.47630_gp,  107.09060_gp,   76.05320_gp,&
         56.39620_gp,   39.66680_gp,  201.54030_gp,  310.80830_gp,  318.73880_gp,&
         261.41050_gp,  217.63470_gp,  186.30930_gp,  154.31180_gp,  282.43840_gp,&
         521.84580_gp,  286.22410_gp,  276.65670_gp,  271.10010_gp,  268.13250_gp,&
         248.78950_gp,  231.00810_gp,  220.48020_gp,  215.40300_gp,  211.42380_gp,&
         200.75620_gp,  318.99680_gp,  284.55050_gp,  258.16170_gp,  237.89530_gp,&
         211.74100_gp,  338.50140_gp,  634.81680_gp,  496.63990_gp,  378.65110_gp,&
         382.32250_gp,  357.96340_gp,  397.03670_gp,  314.91710_gp,  295.10120_gp,&
         274.69180_gp,  265.40330_gp,  264.02010_gp,  406.21170_gp,  372.04100_gp,&
         351.63900_gp,  334.92300_gp,  308.47430_gp,  403.90770_gp,  820.41740_gp,&
         627.70030_gp,  410.25930_gp,  401.83710_gp,  389.38420_gp,  390.58410_gp,&
         374.44730_gp,  392.66670_gp,  369.47670_gp,  374.77260_gp,  351.23330_gp/)
    vdwparams%coeffs(22101:22200)=(/  341.42990_gp,  339.04640_gp,  355.18230_gp,  327.97090_gp,  435.15580_gp,&
         407.01010_gp,  374.38790_gp,  376.34630_gp,  335.43720_gp,  316.64800_gp,&
         303.09840_gp,  289.93290_gp,  281.76770_gp,  438.10480_gp,  404.00650_gp,&
         394.33950_gp,  385.48280_gp,  362.63740_gp,  455.15750_gp,  821.87390_gp,&
         485.13840_gp,  554.58830_gp,  498.32020_gp,  443.80410_gp,  427.90900_gp,&
         503.78650_gp,  119.12060_gp,  117.97580_gp,   88.21580_gp,   71.57230_gp,&
         50.04060_gp,  208.41310_gp,  251.64210_gp,  243.91880_gp,  214.51760_gp,&
         182.61700_gp,  293.11980_gp,  283.87840_gp,  287.58070_gp,  263.34090_gp,&
         201.47070_gp,  173.92920_gp,  170.55570_gp,  197.06800_gp,  184.30860_gp,&
         256.39330_gp,  271.64740_gp,  251.79730_gp,  236.10840_gp,  351.95640_gp,&
         342.03460_gp,  346.02660_gp,  334.71740_gp,  298.85080_gp,  265.86540_gp,&
         26.81420_gp,   18.48640_gp,  364.41020_gp,  222.14890_gp,  155.16280_gp,&
         108.02000_gp,   77.37720_gp,   59.64600_gp,   45.91330_gp,   35.78840_gp,&
         437.64910_gp,  350.61410_gp,  328.94340_gp,  265.58080_gp,  211.73900_gp,&
         178.49350_gp,  148.07630_gp,  122.86250_gp,  715.23760_gp,  607.29210_gp,&
         505.28540_gp,  491.86430_gp,  452.47550_gp,  357.86670_gp,  392.19220_gp,&
         309.33800_gp,  329.81210_gp,  338.08050_gp,  260.41990_gp,  270.17940_gp,&
         318.12240_gp,  285.61020_gp,  247.96070_gp,  225.16700_gp,  199.55490_gp,&
         175.65240_gp,  805.00650_gp,  723.69110_gp,  642.42790_gp,  582.55440_gp,&
         534.62660_gp,  417.98890_gp,  464.12390_gp,  358.48400_gp,  390.89370_gp/)
    vdwparams%coeffs(22201:22300)=(/  363.95130_gp,  304.05500_gp,  322.40840_gp,  399.31090_gp,  374.18700_gp,&
         338.35410_gp,  316.91230_gp,  289.60010_gp,  263.10510_gp,  982.30160_gp,&
         919.91100_gp,  818.62810_gp,  399.23510_gp,  819.04250_gp,  787.98460_gp,&
         768.61590_gp,  750.74320_gp,  734.90660_gp,  587.61560_gp,  648.99040_gp,&
         627.56590_gp,  664.79550_gp,  650.83380_gp,  638.36820_gp,  630.42130_gp,&
         537.19250_gp,  535.11600_gp,  492.90160_gp,  420.19700_gp,  428.65240_gp,&
         391.72610_gp,  361.24810_gp,  302.92430_gp,  284.06330_gp,  292.67440_gp,&
         412.85870_gp,  407.37440_gp,  378.97400_gp,  364.21970_gp,  339.40270_gp,&
         313.98490_gp,  937.38470_gp,  916.24010_gp,  819.40930_gp,  747.58980_gp,&
         737.89240_gp,  714.83440_gp,  729.48730_gp,  707.50950_gp,   42.16270_gp,&
         129.89070_gp,  166.98720_gp,  131.66500_gp,  101.41840_gp,   72.29500_gp,&
         53.78630_gp,   37.99610_gp,  189.13000_gp,  291.49160_gp,  299.52560_gp,&
         246.46800_gp,  205.75330_gp,  176.50780_gp,  146.52980_gp,  265.80330_gp,&
         488.33520_gp,  269.90590_gp,  260.94910_gp,  255.71250_gp,  252.85150_gp,&
         234.92290_gp,  218.25880_gp,  208.33810_gp,  203.52280_gp,  199.56410_gp,&
         189.84580_gp,  300.20530_gp,  268.45380_gp,  244.04670_gp,  225.23300_gp,&
         200.82890_gp,  318.80410_gp,  594.03870_gp,  466.51060_gp,  357.08190_gp,&
         360.54300_gp,  337.81860_gp,  374.03450_gp,  297.66760_gp,  279.03330_gp,&
         259.85840_gp,  251.03670_gp,  249.93600_gp,  382.54570_gp,  351.05910_gp,&
         332.29950_gp,  316.85780_gp,  292.24830_gp,  380.90100_gp,  767.31680_gp/)
    vdwparams%coeffs(22301:22400)=(/  589.47680_gp,  387.57590_gp,  379.63090_gp,  367.91330_gp,  368.93680_gp,&
         353.39960_gp,  370.83500_gp,  349.04790_gp,  353.87000_gp,  331.89420_gp,&
         322.66300_gp,  320.37880_gp,  335.42660_gp,  309.95020_gp,  409.80230_gp,&
         383.80190_gp,  353.52930_gp,  355.00360_gp,  317.37010_gp,  299.78350_gp,&
         287.07700_gp,  274.65620_gp,  267.16680_gp,  412.92460_gp,  381.50770_gp,&
         372.76100_gp,  364.68730_gp,  343.49450_gp,  429.26480_gp,  769.57810_gp,&
         458.25190_gp,  523.05080_gp,  470.40220_gp,  419.48400_gp,  404.59980_gp,&
         474.67630_gp,  112.43910_gp,  111.51080_gp,   83.69460_gp,   68.07200_gp,&
         47.78470_gp,  196.31020_gp,  236.98540_gp,  230.09880_gp,  202.78480_gp,&
         173.02410_gp,  276.28420_gp,  267.90880_gp,  271.43870_gp,  248.60210_gp,&
         190.71700_gp,  164.90150_gp,  161.68180_gp,  186.35590_gp,  174.40700_gp,&
         241.86460_gp,  256.39940_gp,  238.09660_gp,  223.55550_gp,  331.86870_gp,&
         323.01240_gp,  326.94770_gp,  316.29900_gp,  282.88320_gp,  252.00750_gp,&
         239.00710_gp,   26.98030_gp,   18.40260_gp,  390.27660_gp,  230.89960_gp,&
         158.71240_gp,  109.28410_gp,   77.70740_gp,   59.61940_gp,   45.72890_gp,&
         35.55630_gp,  467.80860_gp,  366.32450_gp,  340.74550_gp,  272.09130_gp,&
         215.01450_gp,  180.27000_gp,  148.81480_gp,  122.98220_gp,  767.42050_gp,&
         640.14310_gp,  530.51110_gp,  514.47440_gp,  472.24010_gp,  373.22320_gp,&
         408.04110_gp,  321.48500_gp,  341.37600_gp,  350.71710_gp,  269.98860_gp,&
         278.36200_gp,  328.39780_gp,  292.39200_gp,  252.00000_gp,  227.83110_gp/)
    vdwparams%coeffs(22401:22500)=(/  201.04720_gp,  176.30020_gp,  862.12220_gp,  763.17560_gp,  672.68360_gp,&
         607.30800_gp,  555.72040_gp,  432.27820_gp,  480.90810_gp,  369.43280_gp,&
         402.94980_gp,  374.53940_gp,  313.11430_gp,  331.08150_gp,  412.28890_gp,&
         383.78510_gp,  344.80820_gp,  321.72250_gp,  292.82270_gp,  265.06110_gp,&
         1051.36100_gp,  972.14520_gp,  859.02210_gp,  407.25580_gp,  863.53440_gp,&
         829.85260_gp,  809.20320_gp,  790.17910_gp,  773.31220_gp,  613.39590_gp,&
         684.97150_gp,  661.43730_gp,  698.34480_gp,  683.56120_gp,  670.29400_gp,&
         662.17850_gp,  561.18900_gp,  555.10630_gp,  509.14750_gp,  432.53910_gp,&
         440.58660_gp,  401.22800_gp,  368.98490_gp,  308.49150_gp,  288.95780_gp,&
         297.37620_gp,  425.41930_gp,  417.68300_gp,  386.36610_gp,  370.11470_gp,&
         343.58050_gp,  316.77380_gp,  997.38560_gp,  964.61400_gp,  857.21330_gp,&
         776.36620_gp,  768.96120_gp,  744.79890_gp,  763.57260_gp,  739.96210_gp,&
         42.77140_gp,  134.69310_gp,  172.60210_gp,  134.25050_gp,  102.61550_gp,&
         72.58630_gp,   53.69340_gp,   37.69670_gp,  196.69260_gp,  303.53120_gp,&
         309.63670_gp,  252.11010_gp,  208.87230_gp,  178.28770_gp,  147.27550_gp,&
         274.59280_gp,  514.13850_gp,  276.56690_gp,  267.29850_gp,  261.98070_gp,&
         259.31910_gp,  239.79050_gp,  222.46400_gp,  212.35000_gp,  207.52820_gp,&
         204.25110_gp,  193.02210_gp,  309.23520_gp,  274.29900_gp,  247.91670_gp,&
         227.91320_gp,  202.37580_gp,  328.68890_gp,  625.89490_gp,  484.72570_gp,&
         366.05390_gp,  369.67440_gp,  345.76000_gp,  385.25950_gp,  303.29340_gp/)
    vdwparams%coeffs(22501:22600)=(/  284.18210_gp,  264.43210_gp,  255.64790_gp,  253.66630_gp,  394.03700_gp,&
         359.13800_gp,  338.34720_gp,  321.58990_gp,  295.51790_gp,  390.60480_gp,&
         810.58950_gp,  613.13650_gp,  395.21610_gp,  387.09050_gp,  375.02280_gp,&
         376.42310_gp,  361.71140_gp,  378.48060_gp,  355.95940_gp,  361.49530_gp,&
         338.18990_gp,  328.68340_gp,  326.43780_gp,  342.24650_gp,  315.66390_gp,&
         422.04710_gp,  393.70480_gp,  361.19710_gp,  363.99550_gp,  322.46250_gp,&
         304.19920_gp,  291.08740_gp,  278.53890_gp,  270.11300_gp,  424.70380_gp,&
         389.82000_gp,  379.55200_gp,  370.40460_gp,  347.70470_gp,  440.02170_gp,&
         809.24300_gp,  467.36610_gp,  536.49470_gp,  481.51730_gp,  427.50540_gp,&
         411.94110_gp,  489.14000_gp,  114.55870_gp,  113.35080_gp,   84.32980_gp,&
         68.30290_gp,   47.64220_gp,  201.39880_gp,  243.32100_gp,  235.01080_gp,&
         205.96780_gp,  174.75170_gp,  283.43990_gp,  273.71790_gp,  277.24060_gp,&
         253.98530_gp,  193.59700_gp,  166.75330_gp,  163.58130_gp,  189.73800_gp,&
         177.25080_gp,  247.32620_gp,  261.64130_gp,  241.69560_gp,  226.18930_gp,&
         340.23920_gp,  329.40570_gp,  332.91030_gp,  322.21010_gp,  286.93400_gp,&
         254.77490_gp,  241.34970_gp,  244.65860_gp,   28.06310_gp,   19.02080_gp,&
         412.21710_gp,  242.88250_gp,  166.15330_gp,  113.91260_gp,   80.71180_gp,&
         61.76190_gp,   47.26450_gp,   36.68370_gp,  493.92410_gp,  385.71180_gp,&
         358.02870_gp,  285.01380_gp,  224.52110_gp,  187.80070_gp,  154.66930_gp,&
         127.54830_gp,  809.42550_gp,  674.57400_gp,  558.75930_gp,  541.37410_gp/)
    vdwparams%coeffs(22601:22700)=(/  496.69050_gp,  392.22970_gp,  428.85090_gp,  337.57700_gp,  358.37430_gp,&
         368.40320_gp,  283.32920_gp,  291.87250_gp,  344.56860_gp,  306.11920_gp,&
         263.18210_gp,  237.51920_gp,  209.19710_gp,  183.10950_gp,  908.81520_gp,&
         804.00360_gp,  707.83430_gp,  638.44630_gp,  583.76530_gp,  453.32840_gp,&
         504.66070_gp,  386.97280_gp,  422.32030_gp,  392.35480_gp,  327.80150_gp,&
         346.60550_gp,  432.31830_gp,  401.84330_gp,  360.34200_gp,  335.76050_gp,&
         305.12160_gp,  275.75770_gp, 1107.79880_gp, 1024.09010_gp,  903.99280_gp,&
         425.61430_gp,  909.29190_gp,  873.75470_gp,  851.99520_gp,  831.95410_gp,&
         814.18450_gp,  644.63620_gp,  720.73150_gp,  695.78030_gp,  735.12820_gp,&
         719.56930_gp,  705.59090_gp,  697.11310_gp,  590.04660_gp,  582.95090_gp,&
         534.09830_gp,  453.09060_gp,  461.41030_gp,  419.75000_gp,  385.66670_gp,&
         321.98390_gp,  301.43730_gp,  310.19540_gp,  445.55610_gp,  437.08990_gp,&
         403.69830_gp,  386.31140_gp,  358.11890_gp,  329.72180_gp, 1050.17500_gp,&
         1015.47600_gp,  901.60180_gp,  815.26350_gp,  807.84930_gp,  782.42750_gp,&
         802.90250_gp,  777.97910_gp,   44.64980_gp,  141.56390_gp,  181.18110_gp,&
         140.35990_gp,  106.93200_gp,   75.36360_gp,   55.57440_gp,   38.86970_gp,&
         206.79400_gp,  319.21760_gp,  325.13700_gp,  263.92340_gp,  218.06750_gp,&
         185.73680_gp,  153.07080_gp,  287.82600_gp,  541.07730_gp,  289.49650_gp,&
         279.71050_gp,  274.14900_gp,  271.42500_gp,  250.71870_gp,  232.47500_gp,&
         221.87670_gp,  216.86020_gp,  213.63790_gp,  201.55080_gp,  324.27260_gp/)
    vdwparams%coeffs(22701:22800)=(/  286.98700_gp,  258.86990_gp,  237.60650_gp,  210.59700_gp,  344.26700_gp,&
         658.51960_gp,  508.68390_gp,  383.04250_gp,  386.85740_gp,  361.57810_gp,&
         403.41880_gp,  316.74990_gp,  296.69340_gp,  275.95750_gp,  266.84590_gp,&
         264.59410_gp,  412.94080_gp,  375.73450_gp,  353.48570_gp,  335.59490_gp,&
         307.94100_gp,  408.61880_gp,  852.76800_gp,  643.44430_gp,  412.97590_gp,&
         404.47640_gp,  391.81900_gp,  393.38030_gp,  378.30190_gp,  395.60900_gp,&
         371.96670_gp,  377.90420_gp,  353.32080_gp,  343.35820_gp,  341.04430_gp,&
         357.76850_gp,  329.75260_gp,  442.11340_gp,  411.96690_gp,  377.50500_gp,&
         380.78680_gp,  336.44710_gp,  317.20710_gp,  303.42150_gp,  290.30040_gp,&
         281.28840_gp,  444.68180_gp,  407.53420_gp,  396.43040_gp,  386.55960_gp,&
         362.42110_gp,  460.28590_gp,  850.84330_gp,  488.43100_gp,  561.29360_gp,&
         503.30170_gp,  446.49080_gp,  430.11690_gp,  512.13360_gp,  119.79600_gp,&
         118.34870_gp,   87.73200_gp,   70.89090_gp,   49.26290_gp,  211.08310_gp,&
         254.99480_gp,  245.89890_gp,  215.05640_gp,  182.03880_gp,  296.81460_gp,&
         286.30250_gp,  289.95300_gp,  265.57670_gp,  201.90070_gp,  173.64950_gp,&
         170.37740_gp,  198.10080_gp,  184.95540_gp,  258.78790_gp,  273.62340_gp,&
         252.30300_gp,  235.79430_gp,  356.11480_gp,  344.30530_gp,  347.79920_gp,&
         336.55350_gp,  299.22510_gp,  265.33520_gp,  251.22310_gp,  254.97130_gp,&
         265.87460_gp,   34.89310_gp,   22.90760_gp,  565.26540_gp,  320.28070_gp,&
         213.44840_gp,  143.16020_gp,   99.64350_gp,   75.23510_gp,   56.89220_gp/)
    vdwparams%coeffs(22801:22900)=(/   43.72150_gp,  675.22680_gp,  512.09020_gp,  469.45400_gp,  367.20280_gp,&
         284.61460_gp,  235.33010_gp,  191.58160_gp,  156.31790_gp, 1110.19070_gp,&
         905.11860_gp,  745.82820_gp,  718.55370_gp,  657.10020_gp,  517.54970_gp,&
         564.63620_gp,  442.98550_gp,  468.15800_gp,  482.95770_gp,  370.28130_gp,&
         378.46340_gp,  448.74480_gp,  393.48230_gp,  333.91140_gp,  298.69330_gp,&
         260.59250_gp,  226.02150_gp, 1243.00360_gp, 1078.84470_gp,  940.63630_gp,&
         842.95940_gp,  767.22430_gp,  590.40500_gp,  659.57470_gp,  500.72770_gp,&
         547.33000_gp,  506.98660_gp,  423.10210_gp,  446.11840_gp,  561.90060_gp,&
         517.14790_gp,  458.76760_gp,  424.47180_gp,  382.68090_gp,  343.12810_gp,&
         1514.11660_gp, 1377.51850_gp, 1204.50720_gp,  542.38350_gp, 1218.80630_gp,&
         1169.36560_gp, 1139.78420_gp, 1112.60650_gp, 1088.49650_gp,  851.58160_gp,&
         966.31820_gp,  931.21350_gp,  980.60600_gp,  959.67190_gp,  940.73330_gp,&
         929.93260_gp,  780.85310_gp,  764.11810_gp,  695.34460_gp,  585.88080_gp,&
         595.38900_gp,  538.31140_gp,  492.06190_gp,  407.98400_gp,  380.91340_gp,&
         391.49310_gp,  576.27840_gp,  561.37620_gp,  513.73240_gp,  488.80570_gp,&
         449.85640_gp,  411.28990_gp, 1423.96930_gp, 1358.76840_gp, 1195.96420_gp,&
         1069.52570_gp, 1064.54060_gp, 1030.63770_gp, 1064.44340_gp, 1030.22770_gp,&
         56.51890_gp,  185.91700_gp,  236.67710_gp,  179.17550_gp,  134.30140_gp,&
         92.90860_gp,   67.41930_gp,   46.20820_gp,  272.37350_gp,  421.55730_gp,&
         424.89670_gp,  339.03150_gp,  276.21540_gp,  232.76490_gp,  189.61390_gp/)
    vdwparams%coeffs(22901:23000)=(/  374.05850_gp,  723.98140_gp,  372.01170_gp,  359.10570_gp,  351.98070_gp,&
         348.99300_gp,  320.04310_gp,  295.92090_gp,  282.31300_gp,  276.08710_gp,&
         273.51670_gp,  255.41010_gp,  420.85570_gp,  367.62820_gp,  328.15180_gp,&
         298.82440_gp,  262.44940_gp,  445.76050_gp,  881.57390_gp,  667.64710_gp,&
         492.22910_gp,  497.13820_gp,  463.06700_gp,  521.69780_gp,  402.38250_gp,&
         376.36250_gp,  349.31490_gp,  338.11650_gp,  333.64080_gp,  534.72480_gp,&
         481.41850_gp,  449.34480_gp,  424.08620_gp,  386.28580_gp,  525.07760_gp,&
         1145.39310_gp,  845.71730_gp,  525.89560_gp,  514.99740_gp,  498.56690_gp,&
         501.36870_gp,  484.41020_gp,  504.62370_gp,  473.70770_gp,  482.62190_gp,&
         449.36380_gp,  436.46450_gp,  433.73820_gp,  456.34090_gp,  419.08700_gp,&
         572.39970_gp,  529.74800_gp,  481.97450_gp,  488.96420_gp,  425.23900_gp,&
         399.70400_gp,  381.59490_gp,  364.91310_gp,  351.78910_gp,  573.70310_gp,&
         520.42840_gp,  503.41800_gp,  488.73500_gp,  455.27610_gp,  591.17920_gp,&
         1135.56110_gp,  622.34260_gp,  721.48760_gp,  644.33230_gp,  567.39040_gp,&
         545.60930_gp,  662.24690_gp,  152.92390_gp,  150.12670_gp,  109.20000_gp,&
         87.20760_gp,   59.45110_gp,  272.56770_gp,  329.59200_gp,  315.01770_gp,&
         272.58640_gp,  228.04730_gp,  382.27200_gp,  366.30170_gp,  370.73900_gp,&
         339.43990_gp,  254.57810_gp,  217.24670_gp,  213.33730_gp,  251.22060_gp,&
         233.74770_gp,  331.79130_gp,  349.68000_gp,  319.37190_gp,  296.45520_gp,&
         457.81020_gp,  438.94830_gp,  442.22600_gp,  427.86320_gp,  377.13750_gp/)
    vdwparams%coeffs(23001:23100)=(/  332.10880_gp,  313.58880_gp,  320.40110_gp,  334.98740_gp,  427.91400_gp,&
         37.91220_gp,   24.78330_gp,  611.56120_gp,  348.39490_gp,  232.25130_gp,&
         155.60830_gp,  108.12990_gp,   81.50580_gp,   61.52140_gp,   47.19420_gp,&
         730.65870_gp,  556.62410_gp,  510.60870_gp,  399.57740_gp,  309.59920_gp,&
         255.80130_gp,  208.04020_gp,  169.54570_gp, 1198.29930_gp,  981.82060_gp,&
         809.62470_gp,  780.20170_gp,  713.60600_gp,  561.66940_gp,  613.32300_gp,&
         480.85880_gp,  508.75980_gp,  524.79280_gp,  401.97520_gp,  411.37780_gp,&
         487.87860_gp,  428.01030_gp,  363.15760_gp,  324.71350_gp,  283.10510_gp,&
         245.34090_gp, 1341.61520_gp, 1169.74950_gp, 1020.96270_gp,  915.34550_gp,&
         833.20500_gp,  641.06430_gp,  716.22390_gp,  543.61040_gp,  594.52770_gp,&
         550.71590_gp,  459.12190_gp,  484.57590_gp,  610.33820_gp,  562.12760_gp,&
         498.73950_gp,  461.38680_gp,  415.81730_gp,  372.64800_gp, 1633.63480_gp,&
         1492.46250_gp, 1306.58980_gp,  589.45280_gp, 1320.76830_gp, 1267.60870_gp,&
         1235.63890_gp, 1206.26160_gp, 1180.20390_gp,  923.88260_gp, 1045.82860_gp,&
         1007.90930_gp, 1063.61610_gp, 1040.97360_gp, 1020.49950_gp, 1008.77090_gp,&
         847.32020_gp,  830.14100_gp,  755.59170_gp,  636.37500_gp,  646.85600_gp,&
         584.85060_gp,  534.55500_gp,  442.93230_gp,  413.43530_gp,  425.06050_gp,&
         625.52300_gp,  609.84910_gp,  558.25050_gp,  531.15730_gp,  488.73370_gp,&
         446.66830_gp, 1538.25730_gp, 1472.91170_gp, 1297.80970_gp, 1161.39170_gp,&
         1155.15160_gp, 1118.36270_gp, 1154.43740_gp, 1117.47610_gp,   61.48680_gp/)
    vdwparams%coeffs(23101:23200)=(/  202.20530_gp,  257.35510_gp,  194.88950_gp,  145.92740_gp,  100.78260_gp,&
         72.99230_gp,   49.87780_gp,  296.01550_gp,  458.14250_gp,  462.15210_gp,&
         368.87960_gp,  300.43180_gp,  253.00030_gp,  205.89540_gp,  406.08450_gp,&
         784.42190_gp,  404.40960_gp,  390.26550_gp,  382.49970_gp,  379.20210_gp,&
         347.91510_gp,  321.64840_gp,  306.80410_gp,  300.02320_gp,  297.14060_gp,&
         277.60830_gp,  457.60310_gp,  399.86460_gp,  356.87340_gp,  324.85020_gp,&
         285.12610_gp,  483.84610_gp,  954.64460_gp,  724.29340_gp,  534.82740_gp,&
         540.17920_gp,  503.02360_gp,  566.29020_gp,  437.11840_gp,  408.74340_gp,&
         379.27910_gp,  367.09070_gp,  362.39660_gp,  580.92600_gp,  523.29100_gp,&
         488.48390_gp,  460.96710_gp,  419.74500_gp,  570.29930_gp, 1239.03510_gp,&
         917.12310_gp,  571.51660_gp,  559.66960_gp,  541.79760_gp,  544.80780_gp,&
         526.23920_gp,  548.41680_gp,  514.79430_gp,  524.39740_gp,  488.36000_gp,&
         474.34360_gp,  471.38980_gp,  496.04540_gp,  455.49010_gp,  621.65690_gp,&
         575.41630_gp,  523.55470_gp,  531.03290_gp,  461.90580_gp,  434.06330_gp,&
         414.31110_gp,  396.08640_gp,  381.90670_gp,  622.93810_gp,  565.42470_gp,&
         547.06010_gp,  531.09590_gp,  494.64220_gp,  642.14730_gp, 1229.56280_gp,&
         676.40190_gp,  783.51290_gp,  699.37990_gp,  616.33520_gp,  592.68480_gp,&
         718.54240_gp,  166.42480_gp,  163.21840_gp,  118.59420_gp,   94.56960_gp,&
         64.30220_gp,  296.61470_gp,  358.55000_gp,  342.74070_gp,  296.45720_gp,&
         247.85630_gp,  415.57270_gp,  398.27950_gp,  403.09190_gp,  368.92090_gp/)
    vdwparams%coeffs(23201:23300)=(/  276.46740_gp,  235.83660_gp,  231.59140_gp,  272.86680_gp,  253.87790_gp,&
         360.80170_gp,  380.32910_gp,  347.30490_gp,  322.26350_gp,  497.53160_gp,&
         477.20360_gp,  480.76680_gp,  464.95560_gp,  409.69630_gp,  360.63780_gp,&
         340.45080_gp,  347.81660_gp,  363.73910_gp,  464.96770_gp,  505.41480_gp,&
         36.22770_gp,   24.01410_gp,  546.85470_gp,  321.11700_gp,  217.79050_gp,&
         147.75190_gp,  103.59080_gp,   78.56000_gp,   59.59370_gp,   45.88940_gp,&
         654.57160_gp,  510.28240_gp,  472.33520_gp,  374.05660_gp,  292.71450_gp,&
         243.37850_gp,  199.10820_gp,  163.07830_gp, 1069.42110_gp,  891.82320_gp,&
         738.41210_gp,  714.34970_gp,  654.87180_gp,  515.81520_gp,  564.71100_gp,&
         443.24120_gp,  471.01860_gp,  484.71110_gp,  371.50280_gp,  382.78190_gp,&
         452.96620_gp,  401.00980_gp,  343.04870_gp,  308.28460_gp,  270.15060_gp,&
         235.18820_gp, 1199.46460_gp, 1061.99470_gp,  933.69560_gp,  840.95330_gp,&
         767.84150_gp,  593.96240_gp,  662.26990_gp,  505.60890_gp,  552.81280_gp,&
         513.03450_gp,  427.38710_gp,  452.49980_gp,  566.56420_gp,  525.63100_gp,&
         469.72550_gp,  436.44300_gp,  395.16260_gp,  355.67420_gp, 1461.40830_gp,&
         1352.05270_gp, 1192.22390_gp,  554.49760_gp, 1199.39140_gp, 1152.42600_gp,&
         1123.71950_gp, 1097.29970_gp, 1073.87920_gp,  847.61650_gp,  948.67820_gp,&
         915.60320_gp,  969.49560_gp,  949.02750_gp,  930.61310_gp,  919.59590_gp,&
         776.79900_gp,  766.67010_gp,  700.99760_gp,  592.61750_gp,  603.37400_gp,&
         547.65300_gp,  502.12300_gp,  417.48050_gp,  390.19750_gp,  401.70990_gp/)
    vdwparams%coeffs(23301:23400)=(/  581.91410_gp,  570.46250_gp,  525.53390_gp,  501.88930_gp,  463.84860_gp,&
         425.62050_gp, 1384.04100_gp, 1339.43530_gp, 1188.00750_gp, 1071.36640_gp,&
         1061.69300_gp, 1028.06880_gp, 1056.13960_gp, 1023.20180_gp,   58.19620_gp,&
         186.83530_gp,  238.60370_gp,  183.40000_gp,  138.53490_gp,   96.57240_gp,&
         70.46340_gp,   48.57660_gp,  272.68930_gp,  421.46260_gp,  428.44370_gp,&
         345.92360_gp,  284.14620_gp,  240.67810_gp,  197.03300_gp,  376.67550_gp,&
         713.37770_gp,  378.51240_gp,  365.42570_gp,  358.08830_gp,  354.59770_gp,&
         327.08070_gp,  302.88030_gp,  288.91150_gp,  282.39730_gp,  278.51920_gp,&
         262.16840_gp,  425.84750_gp,  375.44380_gp,  337.27630_gp,  308.39060_gp,&
         272.00760_gp,  449.76690_gp,  867.48280_gp,  667.77030_gp,  500.31680_gp,&
         505.20350_gp,  471.39900_gp,  527.00000_gp,  411.74300_gp,  385.21960_gp,&
         357.80550_gp,  346.02160_gp,  342.93620_gp,  540.66240_gp,  490.72030_gp,&
         460.51500_gp,  436.16550_gp,  398.86600_gp,  533.23930_gp, 1122.90680_gp,&
         844.59250_gp,  537.85080_gp,  526.73890_gp,  510.10490_gp,  512.39170_gp,&
         493.17410_gp,  515.59660_gp,  484.42780_gp,  492.51520_gp,  459.97740_gp,&
         446.93120_gp,  444.02250_gp,  466.51200_gp,  429.27350_gp,  578.61310_gp,&
         537.90480_gp,  491.61910_gp,  496.67380_gp,  436.48690_gp,  410.80070_gp,&
         392.46080_gp,  375.14240_gp,  363.02260_gp,  580.68260_gp,  530.94590_gp,&
         515.73940_gp,  502.11810_gp,  469.48140_gp,  600.69020_gp, 1119.52150_gp,&
         636.42920_gp,  732.57850_gp,  655.41880_gp,  580.44230_gp,  558.78270_gp/)
    vdwparams%coeffs(23401:23500)=(/  668.42250_gp,  156.75410_gp,  154.13350_gp,  113.15480_gp,   90.70430_gp,&
         62.19540_gp,  277.27810_gp,  334.86500_gp,  321.95190_gp,  280.21780_gp,&
         235.81280_gp,  388.44680_gp,  373.99630_gp,  378.64590_gp,  346.45320_gp,&
         261.54130_gp,  224.11260_gp,  219.95380_gp,  257.27540_gp,  239.87350_gp,&
         338.47720_gp,  357.62210_gp,  328.50830_gp,  305.97730_gp,  465.36830_gp,&
         449.02050_gp,  453.15700_gp,  438.04100_gp,  387.88180_gp,  342.76350_gp,&
         324.03570_gp,  329.56900_gp,  344.18570_gp,  436.65200_gp,  474.67200_gp,&
         448.05730_gp,   35.24610_gp,   23.58510_gp,  512.71640_gp,  305.76290_gp,&
         209.45590_gp,  143.19470_gp,  100.97450_gp,   76.88630_gp,   58.52230_gp,&
         45.18360_gp,  614.35230_gp,  484.50380_gp,  450.72670_gp,  359.36400_gp,&
         282.89020_gp,  236.14560_gp,  193.92690_gp,  159.36080_gp, 1002.55590_gp,&
         843.05400_gp,  699.45860_gp,  678.14850_gp,  622.46350_gp,  490.67860_gp,&
         537.74700_gp,  422.50680_gp,  449.86640_gp,  462.32080_gp,  354.64230_gp,&
         366.61620_gp,  433.26650_gp,  385.51090_gp,  331.38440_gp,  298.73310_gp,&
         262.62260_gp,  229.31350_gp, 1125.71340_gp, 1003.86650_gp,  885.95760_gp,&
         799.94710_gp,  731.67410_gp,  567.83410_gp,  632.35020_gp,  484.47400_gp,&
         529.45790_gp,  491.88120_gp,  409.79950_gp,  434.43580_gp,  542.08910_gp,&
         504.87630_gp,  453.02540_gp,  422.01870_gp,  383.18370_gp,  345.82830_gp,&
         1372.24220_gp, 1276.89730_gp, 1130.12760_gp,  534.52440_gp, 1134.16810_gp,&
         1090.35360_gp, 1063.35160_gp, 1038.47920_gp, 1016.43710_gp,  805.96520_gp/)
    vdwparams%coeffs(23501:23600)=(/  897.02580_gp,  866.41260_gp,  918.41900_gp,  899.09510_gp,  881.76000_gp,&
         871.14240_gp,  738.17630_gp,  731.34130_gp,  670.40070_gp,  568.12010_gp,&
         578.90380_gp,  526.61010_gp,  483.71080_gp,  403.07740_gp,  377.05920_gp,&
         388.40990_gp,  557.67540_gp,  548.20960_gp,  506.81290_gp,  485.04920_gp,&
         449.46500_gp,  413.43530_gp, 1303.49170_gp, 1267.56670_gp, 1128.02060_gp,&
         1021.55310_gp, 1010.49480_gp,  978.59690_gp, 1002.71720_gp,  971.86890_gp,&
         56.27690_gp,  178.17450_gp,  228.03400_gp,  176.77650_gp,  134.27420_gp,&
         94.16220_gp,   69.03990_gp,   47.87710_gp,  259.66980_gp,  401.03040_gp,&
         409.36850_gp,  332.69600_gp,  274.67950_gp,  233.51090_gp,  191.89600_gp,&
         360.37990_gp,  675.20960_gp,  363.77550_gp,  351.31790_gp,  344.23720_gp,&
         340.67740_gp,  315.12290_gp,  292.09790_gp,  278.65290_gp,  272.30380_gp,&
         267.95690_gp,  253.25200_gp,  407.86260_gp,  361.39310_gp,  325.90920_gp,&
         298.82550_gp,  264.38680_gp,  430.88900_gp,  820.90550_gp,  636.78340_gp,&
         480.83080_gp,  485.46850_gp,  453.54170_gp,  505.18300_gp,  397.28740_gp,&
         371.84930_gp,  345.61400_gp,  334.07750_gp,  331.75670_gp,  518.04080_gp,&
         472.12430_gp,  444.39750_gp,  421.81780_gp,  386.75500_gp,  512.46920_gp,&
         1061.41850_gp,  805.01180_gp,  518.60750_gp,  507.91390_gp,  491.98270_gp,&
         493.89760_gp,  474.46600_gp,  496.85260_gp,  467.06730_gp,  474.37630_gp,&
         443.71620_gp,  431.21480_gp,  428.33450_gp,  449.58170_gp,  414.22170_gp,&
         554.55640_gp,  516.81140_gp,  473.54520_gp,  477.37190_gp,  421.95370_gp/)
    vdwparams%coeffs(23601:23700)=(/  397.50830_gp,  379.98620_gp,  363.22860_gp,  352.17630_gp,  557.05210_gp,&
         511.33530_gp,  497.76870_gp,  485.42330_gp,  454.93030_gp,  577.42140_gp,&
         1060.74320_gp,  613.56050_gp,  703.96650_gp,  630.75620_gp,  560.00710_gp,&
         539.43990_gp,  640.68280_gp,  151.12750_gp,  148.89980_gp,  110.00870_gp,&
         88.49530_gp,   61.02580_gp,  266.13770_gp,  321.31290_gp,  309.94780_gp,&
         270.80050_gp,  228.81280_gp,  373.02450_gp,  360.07160_gp,  364.62620_gp,&
         333.63080_gp,  253.01950_gp,  217.39930_gp,  213.29230_gp,  248.38210_gp,&
         231.86030_gp,  325.67860_gp,  344.53460_gp,  317.59860_gp,  296.51530_gp,&
         447.15550_gp,  432.84780_gp,  437.27120_gp,  422.65520_gp,  375.37250_gp,&
         332.49710_gp,  314.60820_gp,  319.17350_gp,  333.02450_gp,  420.51550_gp,&
         457.07930_gp,  432.69410_gp,  418.57890_gp,   44.83210_gp,   29.18640_gp,&
         750.32270_gp,  417.86520_gp,  276.34300_gp,  184.34620_gp,  127.74190_gp,&
         96.09500_gp,   72.38520_gp,   55.40920_gp,  895.16260_gp,  669.76420_gp,&
         611.42610_gp,  475.65060_gp,  367.14490_gp,  302.79080_gp,  245.86350_gp,&
         200.10940_gp, 1476.86940_gp, 1190.53060_gp,  978.58840_gp,  940.93240_gp,&
         859.37300_gp,  676.81130_gp,  737.08980_gp,  578.08990_gp,  609.28650_gp,&
         629.24930_gp,  482.41490_gp,  491.20730_gp,  583.52030_gp,  509.45200_gp,&
         430.80210_gp,  384.58730_gp,  334.81150_gp,  289.79810_gp, 1652.28320_gp,&
         1419.94170_gp, 1233.00060_gp, 1102.26380_gp, 1001.66600_gp,  768.83820_gp,&
         859.63920_gp,  650.71280_gp,  711.08920_gp,  657.96470_gp,  549.42990_gp/)
    vdwparams%coeffs(23701:23800)=(/  578.10960_gp,  730.51140_gp,  669.79600_gp,  592.20770_gp,  546.93930_gp,&
         492.13770_gp,  440.46680_gp, 2013.02820_gp, 1816.06750_gp, 1581.46710_gp,&
         700.80410_gp, 1604.09060_gp, 1537.92730_gp, 1498.69410_gp, 1462.68340_gp,&
         1430.73060_gp, 1114.51140_gp, 1273.39160_gp, 1226.12320_gp, 1287.46720_gp,&
         1259.80960_gp, 1234.73440_gp, 1220.75570_gp, 1021.97250_gp,  996.37850_gp,&
         904.60210_gp,  760.98160_gp,  772.54200_gp,  697.08220_gp,  636.15080_gp,&
         526.55380_gp,  491.21220_gp,  504.41110_gp,  748.64790_gp,  727.00690_gp,&
         663.22010_gp,  630.00140_gp,  578.68480_gp,  528.16950_gp, 1887.10610_gp,&
         1787.74060_gp, 1567.07840_gp, 1395.45190_gp, 1391.82920_gp, 1347.30760_gp,&
         1394.84040_gp, 1349.30080_gp,   72.90130_gp,  242.28720_gp,  308.14100_gp,&
         231.65430_gp,  173.02490_gp,  119.14340_gp,   86.07730_gp,   58.58360_gp,&
         355.34390_gp,  550.67610_gp,  552.86340_gp,  438.85970_gp,  356.26970_gp,&
         299.50590_gp,  243.34210_gp,  486.93460_gp,  952.57710_gp,  481.82380_gp,&
         465.02110_gp,  455.78960_gp,  452.13520_gp,  413.43280_gp,  381.95910_gp,&
         364.37600_gp,  356.38920_gp,  353.72270_gp,  329.13430_gp,  546.64020_gp,&
         475.53820_gp,  423.29180_gp,  384.76970_gp,  337.23270_gp,  579.79800_gp,&
         1160.90400_gp,  872.28820_gp,  638.07920_gp,  644.45150_gp,  599.66140_gp,&
         678.12320_gp,  519.57380_gp,  485.80720_gp,  450.60820_gp,  436.21390_gp,&
         429.55510_gp,  694.25600_gp,  622.68260_gp,  579.80000_gp,  546.38720_gp,&
         496.79900_gp,  681.34670_gp, 1511.34300_gp, 1105.96430_gp,  679.35610_gp/)
    vdwparams%coeffs(23801:23900)=(/  665.24220_gp,  643.90080_gp,  647.87420_gp,  626.96840_gp,  652.18100_gp,&
         611.90030_gp,  624.06450_gp,  580.19150_gp,  563.43460_gp,  559.99100_gp,&
         589.59320_gp,  540.94020_gp,  743.94190_gp,  686.99120_gp,  623.58700_gp,&
         633.81630_gp,  548.31610_gp,  514.95440_gp,  491.33470_gp,  469.79990_gp,&
         452.04520_gp,  744.45260_gp,  672.69630_gp,  649.47180_gp,  629.76130_gp,&
         585.66430_gp,  766.95690_gp, 1494.44040_gp,  804.04920_gp,  935.36030_gp,&
         834.46090_gp,  732.60640_gp,  704.02570_gp,  860.65850_gp,  197.58390_gp,&
         193.83810_gp,  140.29040_gp,  111.70570_gp,   75.74170_gp,  353.05310_gp,&
         427.39490_gp,  407.46200_gp,  351.68440_gp,  293.41430_gp,  495.41720_gp,&
         473.69840_gp,  479.34490_gp,  438.91900_gp,  328.04360_gp,  279.29530_gp,&
         274.28730_gp,  324.04980_gp,  301.16290_gp,  429.32210_gp,  452.01960_gp,&
         411.81840_gp,  381.69420_gp,  593.35720_gp,  567.24450_gp,  571.02950_gp,&
         552.67180_gp,  486.02980_gp,  427.18180_gp,  403.02800_gp,  412.53210_gp,&
         431.49070_gp,  553.08270_gp,  600.97790_gp,  563.10510_gp,  541.67760_gp,&
         716.54750_gp,   46.82020_gp,   30.58860_gp,  748.90520_gp,  428.04980_gp,&
         286.02900_gp,  191.93550_gp,  133.44100_gp,  100.55750_gp,   75.82980_gp,&
         58.07950_gp,  894.74480_gp,  683.23710_gp,  627.53480_gp,  491.89290_gp,&
         381.67570_gp,  315.61930_gp,  256.84550_gp,  209.37070_gp, 1467.24050_gp,&
         1204.10910_gp,  993.28050_gp,  957.57000_gp,  876.00500_gp,  689.47940_gp,&
         753.09320_gp,  590.43900_gp,  625.01030_gp,  644.50280_gp,  493.60530_gp/)
    vdwparams%coeffs(23901:24000)=(/  505.58370_gp,  599.83010_gp,  526.86600_gp,  447.55460_gp,  400.47040_gp,&
         349.37100_gp,  302.89600_gp, 1643.08790_gp, 1434.65830_gp, 1253.13170_gp,&
         1124.04570_gp, 1023.50000_gp,  787.90450_gp,  880.03140_gp,  668.28590_gp,&
         730.81350_gp,  677.02710_gp,  564.24900_gp,  595.73220_gp,  749.99110_gp,&
         691.36370_gp,  613.97930_gp,  568.35230_gp,  512.53890_gp,  459.57350_gp,&
         2001.25570_gp, 1830.26840_gp, 1603.54520_gp,  725.68940_gp, 1619.66840_gp,&
         1554.63760_gp, 1515.46110_gp, 1479.45840_gp, 1447.52880_gp, 1134.13260_gp,&
         1282.37790_gp, 1236.06680_gp, 1304.71910_gp, 1276.95670_gp, 1251.87390_gp,&
         1237.43050_gp, 1039.96910_gp, 1019.98090_gp,  928.84540_gp,  782.60680_gp,&
         795.60020_gp,  719.59790_gp,  657.87990_gp,  545.19970_gp,  508.86690_gp,&
         523.22550_gp,  768.79970_gp,  749.98940_gp,  687.06370_gp,  654.05030_gp,&
         602.15770_gp,  550.60100_gp, 1885.58670_gp, 1807.11530_gp, 1593.00870_gp,&
         1426.54000_gp, 1418.23590_gp, 1373.05490_gp, 1416.44570_gp, 1371.18900_gp,&
         75.83940_gp,  248.48670_gp,  316.48150_gp,  240.13610_gp,  180.03260_gp,&
         124.41240_gp,   90.09320_gp,   61.44990_gp,  363.46120_gp,  562.68080_gp,&
         568.14570_gp,  454.21760_gp,  370.38990_gp,  312.15000_gp,  254.18510_gp,&
         499.20670_gp,  962.26620_gp,  497.62570_gp,  480.20280_gp,  470.59030_gp,&
         466.42410_gp,  428.18050_gp,  395.89450_gp,  377.58950_gp,  369.19520_gp,&
         365.39900_gp,  341.73160_gp,  562.76860_gp,  492.36040_gp,  439.83840_gp,&
         400.63200_gp,  351.85020_gp,  594.99260_gp, 1171.12050_gp,  889.95010_gp/)
    vdwparams%coeffs(24001:24100)=(/  658.20980_gp,  664.73460_gp,  619.09740_gp,  696.37520_gp,  538.17200_gp,&
         503.19100_gp,  466.88550_gp,  451.74170_gp,  446.17650_gp,  714.06170_gp,&
         643.81770_gp,  601.43340_gp,  567.86370_gp,  517.38510_gp,  702.11480_gp,&
         1519.79330_gp, 1126.91360_gp,  703.63390_gp,  689.04640_gp,  667.06290_gp,&
         670.68630_gp,  647.47300_gp,  675.11960_gp,  633.75280_gp,  645.44130_gp,&
         601.28460_gp,  584.04680_gp,  580.39470_gp,  610.65540_gp,  560.86530_gp,&
         764.59540_gp,  708.05680_gp,  644.52930_gp,  653.36840_gp,  568.93970_gp,&
         534.65980_gp,  510.29930_gp,  487.73540_gp,  470.44860_gp,  765.80950_gp,&
         695.66310_gp,  673.41930_gp,  654.03890_gp,  609.46360_gp,  790.62920_gp,&
         1508.96090_gp,  832.78590_gp,  963.89140_gp,  860.55570_gp,  758.68030_gp,&
         729.61370_gp,  883.20780_gp,  205.08500_gp,  201.20830_gp,  146.33300_gp,&
         116.68030_gp,   79.28770_gp,  364.91620_gp,  441.20040_gp,  422.12020_gp,&
         365.44980_gp,  305.79960_gp,  511.30630_gp,  490.31810_gp,  496.24460_gp,&
         454.08320_gp,  340.51320_gp,  290.56100_gp,  285.26770_gp,  335.81530_gp,&
         312.48780_gp,  444.15730_gp,  468.38860_gp,  428.09060_gp,  397.44450_gp,&
         612.29940_gp,  587.72220_gp,  592.23960_gp,  572.68150_gp,  504.89560_gp,&
         444.57570_gp,  419.69160_gp,  428.42590_gp,  447.92670_gp,  572.08230_gp,&
         621.93900_gp,  584.52180_gp,  563.11340_gp,  739.91600_gp,  766.02020_gp,&
         35.21950_gp,   23.41060_gp,  525.31240_gp,  309.84670_gp,  210.84170_gp,&
         143.42000_gp,  100.74190_gp,   76.48420_gp,   58.05670_gp,   44.71200_gp/)
    vdwparams%coeffs(24101:24200)=(/  628.93080_gp,  491.88460_gp,  456.02900_gp,  361.94550_gp,  283.82390_gp,&
         236.32710_gp,  193.60200_gp,  158.74630_gp, 1027.45370_gp,  858.73510_gp,&
         711.40130_gp,  688.66520_gp,  631.54370_gp,  497.57540_gp,  544.86260_gp,&
         427.80560_gp,  454.84270_gp,  467.85600_gp,  358.68500_gp,  369.92010_gp,&
         437.68410_gp,  388.10850_gp,  332.55870_gp,  299.19670_gp,  262.48890_gp,&
         228.75780_gp, 1152.81500_gp, 1022.66980_gp,  900.10840_gp,  811.30330_gp,&
         741.16010_gp,  573.92370_gp,  639.63880_gp,  488.87360_gp,  534.38470_gp,&
         496.07390_gp,  413.27390_gp,  437.68030_gp,  547.43800_gp,  508.47990_gp,&
         454.99670_gp,  423.13470_gp,  383.49110_gp,  345.49990_gp, 1404.89710_gp,&
         1301.74880_gp, 1149.09500_gp,  537.12000_gp, 1155.04230_gp, 1109.98290_gp,&
         1082.37510_gp, 1056.96040_gp, 1034.43330_gp,  817.61230_gp,  913.64540_gp,&
         881.98390_gp,  934.09640_gp,  914.38960_gp,  896.67750_gp,  886.00250_gp,&
         749.11050_gp,  740.25210_gp,  677.37500_gp,  573.10850_gp,  583.63200_gp,&
         530.09450_gp,  486.29160_gp,  404.60200_gp,  378.24260_gp,  389.43520_gp,&
         562.59750_gp,  551.95500_gp,  509.03810_gp,  486.47930_gp,  450.00480_gp,&
         413.26650_gp, 1331.70580_gp, 1290.41270_gp, 1145.48090_gp, 1034.23170_gp,&
         1024.36030_gp,  991.94790_gp, 1018.21410_gp,  986.57360_gp,   56.45130_gp,&
         180.35670_gp,  230.52060_gp,  177.69180_gp,  134.51000_gp,   93.95690_gp,&
         68.65550_gp,   47.37020_gp,  263.07010_gp,  406.58200_gp,  413.82490_gp,&
         334.84900_gp,  275.54270_gp,  233.69990_gp,  191.57960_gp,  364.08150_gp/)
    vdwparams%coeffs(24201:24300)=(/  687.23560_gp,  366.30740_gp,  353.66860_gp,  346.54410_gp,  343.08700_gp,&
         316.71570_gp,  293.36460_gp,  279.83390_gp,  273.49300_gp,  269.53150_gp,&
         254.02930_gp,  411.65060_gp,  363.51880_gp,  326.99810_gp,  299.29580_gp,&
         264.27820_gp,  434.94430_gp,  835.72700_gp,  644.79280_gp,  484.25820_gp,&
         488.96320_gp,  456.41490_gp,  509.66540_gp,  398.98040_gp,  373.31410_gp,&
         346.80080_gp,  335.29990_gp,  332.49420_gp,  522.63660_gp,  474.95740_gp,&
         446.16040_gp,  422.89020_gp,  387.08050_gp,  516.25270_gp, 1081.53380_gp,&
         815.48330_gp,  521.03820_gp,  510.28040_gp,  494.20110_gp,  496.32490_gp,&
         477.41750_gp,  499.38720_gp,  469.26950_gp,  476.95560_gp,  445.65950_gp,&
         433.04490_gp,  430.20330_gp,  451.84050_gp,  415.95250_gp,  559.54550_gp,&
         520.58420_gp,  476.17010_gp,  480.71520_gp,  423.22750_gp,  398.43450_gp,&
         380.70180_gp,  363.88290_gp,  352.32230_gp,  561.57040_gp,  514.06470_gp,&
         499.68080_gp,  486.75350_gp,  455.47580_gp,  581.59320_gp, 1079.04440_gp,&
         616.50660_gp,  708.90840_gp,  634.54270_gp,  562.36500_gp,  541.47420_gp,&
         646.30030_gp,  151.87100_gp,  149.45570_gp,  109.95340_gp,   88.23330_gp,&
         60.59270_gp,  268.15940_gp,  323.87400_gp,  311.74680_gp,  271.70710_gp,&
         228.98300_gp,  375.83520_gp,  362.14280_gp,  366.66280_gp,  335.47490_gp,&
         253.63880_gp,  217.51740_gp,  213.43980_gp,  249.27270_gp,  232.49050_gp,&
         327.67280_gp,  346.36060_gp,  318.55200_gp,  296.96280_gp,  450.40960_gp,&
         435.02220_gp,  439.16880_gp,  424.51400_gp,  376.28530_gp,  332.76620_gp/)
    vdwparams%coeffs(24301:24400)=(/  314.65980_gp,  319.73310_gp,  333.78790_gp,  422.75520_gp,  459.55200_gp,&
         434.19590_gp,  419.56530_gp,  545.26090_gp,  566.23920_gp,  420.97660_gp,&
         36.59510_gp,   24.28140_gp,  566.10590_gp,  326.44690_gp,  220.41970_gp,&
         149.32900_gp,  104.67690_gp,   79.39680_gp,   60.23930_gp,   46.38780_gp,&
         676.87130_gp,  519.98460_gp,  479.79410_gp,  378.62690_gp,  295.80530_gp,&
         245.86750_gp,  201.13960_gp,  164.78050_gp, 1112.62640_gp,  914.47030_gp,&
         755.25190_gp,  729.56490_gp,  668.15990_gp,  526.76740_gp,  575.38750_gp,&
         451.95440_gp,  478.78370_gp,  493.06650_gp,  378.38780_gp,  388.37940_gp,&
         459.98760_gp,  406.00970_gp,  346.76900_gp,  311.49230_gp,  272.90380_gp,&
         237.59510_gp, 1247.39440_gp, 1090.10540_gp,  954.67960_gp,  858.10820_gp,&
         782.68490_gp,  604.71660_gp,  674.51160_gp,  514.26690_gp,  561.75990_gp,&
         521.00450_gp,  434.78170_gp,  459.16630_gp,  575.97200_gp,  532.73890_gp,&
         475.16010_gp,  441.17360_gp,  399.23870_gp,  359.25420_gp, 1520.97860_gp,&
         1390.82370_gp, 1221.27130_gp,  561.44200_gp, 1232.13480_gp, 1182.77560_gp,&
         1153.02650_gp, 1125.67190_gp, 1101.41510_gp,  866.43530_gp,  977.22080_gp,&
         942.63090_gp,  993.14900_gp,  972.01510_gp,  952.96410_gp,  941.78690_gp,&
         793.83690_gp,  780.53540_gp,  712.53240_gp,  602.20240_gp,  612.58210_gp,&
         555.38920_gp,  508.82230_gp,  423.01850_gp,  395.32750_gp,  406.61330_gp,&
         591.90940_gp,  578.59300_gp,  531.93140_gp,  507.57870_gp,  468.77360_gp,&
         429.97720_gp, 1434.83370_gp, 1375.14650_gp, 1214.82590_gp, 1091.79910_gp/)
    vdwparams%coeffs(24401:24500)=(/ 1084.32090_gp, 1049.87920_gp, 1080.81280_gp, 1046.58580_gp,   58.80380_gp,&
         189.86110_gp,  242.44290_gp,  185.57490_gp,  140.12490_gp,   97.65680_gp,&
         71.26440_gp,   49.12190_gp,  277.53610_gp,  429.32400_gp,  434.97780_gp,&
         350.08490_gp,  287.16960_gp,  243.16080_gp,  199.05370_gp,  383.51900_gp,&
         732.99270_gp,  383.59320_gp,  370.45440_gp,  363.03280_gp,  359.64520_gp,&
         331.00010_gp,  306.43150_gp,  292.35680_gp,  285.79970_gp,  282.25140_gp,&
         265.01660_gp,  432.12490_gp,  379.97370_gp,  340.92760_gp,  311.61210_gp,&
         274.78910_gp,  457.84580_gp,  892.53360_gp,  682.12820_gp,  507.66120_gp,&
         512.56950_gp,  478.19280_gp,  536.26360_gp,  416.99040_gp,  390.19030_gp,&
         362.40560_gp,  350.51370_gp,  346.80430_gp,  549.09790_gp,  497.00410_gp,&
         465.75860_gp,  440.87250_gp,  402.97870_gp,  541.68560_gp, 1158.62560_gp,&
         863.72770_gp,  544.54800_gp,  533.28550_gp,  516.41630_gp,  518.93410_gp,&
         500.05910_gp,  522.12860_gp,  490.44870_gp,  499.03390_gp,  465.55820_gp,&
         452.30800_gp,  449.38570_gp,  472.21050_gp,  434.37410_gp,  588.33640_gp,&
         546.16030_gp,  498.48100_gp,  504.25710_gp,  441.78790_gp,  415.72770_gp,&
         397.14540_gp,  379.72750_gp,  367.01200_gp,  590.00540_gp,  537.89130_gp,&
         521.81620_gp,  507.73170_gp,  474.44690_gp,  610.06800_gp, 1151.65020_gp,&
         644.31070_gp,  743.92250_gp,  665.68850_gp,  587.86120_gp,  565.70590_gp,&
         680.49130_gp,  158.44380_gp,  155.95490_gp,  114.35140_gp,   91.69350_gp,&
         62.89910_gp,  280.61570_gp,  339.28390_gp,  325.70670_gp,  283.27190_gp/)
    vdwparams%coeffs(24501:24600)=(/  238.25580_gp,  393.78620_gp,  378.62350_gp,  383.30360_gp,  350.90960_gp,&
         264.69190_gp,  226.64920_gp,  222.44630_gp,  260.44220_gp,  242.68500_gp,&
         342.72530_gp,  361.82680_gp,  332.02940_gp,  309.17670_gp,  472.00610_gp,&
         454.50240_gp,  458.50240_gp,  443.54250_gp,  392.44480_gp,  346.62740_gp,&
         327.63810_gp,  333.57090_gp,  348.31840_gp,  442.28570_gp,  480.56320_gp,&
         453.04590_gp,  437.32780_gp,  571.16140_gp,  591.93970_gp,  439.11410_gp,&
         458.76770_gp,   34.43520_gp,   23.23010_gp,  499.48070_gp,  296.73610_gp,&
         203.61200_gp,  139.67130_gp,   98.85410_gp,   75.50880_gp,   57.64380_gp,&
         44.61450_gp,  598.53480_gp,  470.40460_gp,  437.59830_gp,  349.16910_gp,&
         275.38760_gp,  230.36660_gp,  189.64180_gp,  156.23470_gp,  979.19910_gp,&
         820.27660_gp,  680.15970_gp,  659.48500_gp,  605.31050_gp,  477.71160_gp,&
         522.92780_gp,  411.37250_gp,  437.43830_gp,  449.47760_gp,  345.33400_gp,&
         356.54050_gp,  421.15000_gp,  374.85400_gp,  322.64010_gp,  291.26630_gp,&
         256.51700_gp,  224.42600_gp, 1099.74810_gp,  977.35030_gp,  861.99720_gp,&
         778.24330_gp,  711.94490_gp,  553.09640_gp,  615.62580_gp,  472.21430_gp,&
         515.56800_gp,  479.05900_gp,  399.72620_gp,  423.22300_gp,  527.65140_gp,&
         491.29650_gp,  441.10510_gp,  411.23900_gp,  373.83140_gp,  337.86400_gp,&
         1340.91930_gp, 1244.07560_gp, 1100.21380_gp,  520.77890_gp, 1104.86750_gp,&
         1061.99280_gp, 1035.62820_gp, 1011.34170_gp,  989.81510_gp,  784.93070_gp,&
         875.00190_gp,  845.02380_gp,  894.10530_gp,  875.23480_gp,  858.30520_gp/)
    vdwparams%coeffs(24601:24700)=(/  847.94180_gp,  718.50660_gp,  711.39240_gp,  652.27240_gp,  553.38510_gp,&
         563.76710_gp,  513.07710_gp,  471.51500_gp,  393.47970_gp,  368.27210_gp,&
         379.16720_gp,  543.57710_gp,  534.01680_gp,  493.81190_gp,  472.81700_gp,&
         438.50630_gp,  403.79910_gp, 1273.00290_gp, 1234.74440_gp, 1097.90520_gp,&
         994.18690_gp,  984.02510_gp,  953.02030_gp,  976.74380_gp,  946.60700_gp,&
         54.79340_gp,  172.98830_gp,  221.55260_gp,  172.03160_gp,  131.07440_gp,&
         92.27610_gp,   67.91190_gp,   47.30860_gp,  252.28990_gp,  389.57210_gp,&
         397.54110_gp,  323.38900_gp,  267.45480_gp,  227.81180_gp,  187.66470_gp,&
         351.18760_gp,  657.80950_gp,  354.08760_gp,  342.06260_gp,  335.19530_gp,&
         331.74290_gp,  306.79020_gp,  284.47550_gp,  271.44430_gp,  265.25960_gp,&
         261.03120_gp,  246.69600_gp,  396.51060_gp,  351.52430_gp,  317.35850_gp,&
         291.36250_gp,  258.22570_gp,  420.13560_gp,  800.28360_gp,  620.27440_gp,&
         468.37410_gp,  472.94180_gp,  442.04910_gp,  492.45940_gp,  387.41710_gp,&
         362.76960_gp,  337.32470_gp,  326.05860_gp,  323.65690_gp,  504.25970_gp,&
         459.56860_gp,  432.76910_gp,  411.05850_gp,  377.29830_gp,  499.58850_gp,&
         1035.58210_gp,  784.41150_gp,  505.32050_gp,  494.91310_gp,  479.43010_gp,&
         481.26120_gp,  462.37370_gp,  484.03110_gp,  455.10020_gp,  462.20710_gp,&
         432.36450_gp,  420.19360_gp,  417.35840_gp,  437.83510_gp,  403.59860_gp,&
         540.10000_gp,  503.56000_gp,  461.65670_gp,  465.29970_gp,  411.68040_gp,&
         388.05470_gp,  371.09900_gp,  354.87160_gp,  344.08260_gp,  542.82280_gp/)
    vdwparams%coeffs(24701:24800)=(/  498.22260_gp,  485.04830_gp,  473.18900_gp,  443.81290_gp,  562.85700_gp,&
         1034.45190_gp,  597.72660_gp,  685.93330_gp,  615.04850_gp,  546.03460_gp,&
         526.04590_gp,  624.73940_gp,  146.93770_gp,  145.06890_gp,  107.51640_gp,&
         86.74430_gp,   60.10130_gp,  258.46390_gp,  312.19790_gp,  301.36080_gp,&
         263.70240_gp,  223.25700_gp,  363.00510_gp,  350.46740_gp,  354.93010_gp,&
         324.91900_gp,  246.98040_gp,  212.43170_gp,  208.38910_gp,  242.21550_gp,&
         226.17060_gp,  316.84610_gp,  335.20830_gp,  309.32650_gp,  289.13180_gp,&
         435.47330_gp,  421.56130_gp,  425.95070_gp,  411.95670_gp,  366.33030_gp,&
         324.83190_gp,  307.50010_gp,  311.80890_gp,  325.14160_gp,  409.61940_gp,&
         444.99610_gp,  421.45120_gp,  407.93570_gp,  527.80130_gp,  548.42160_gp,&
         408.87060_gp,  426.35470_gp,  398.09730_gp,   31.59280_gp,   21.62100_gp,&
         437.57960_gp,  264.83780_gp,  184.01940_gp,  127.52710_gp,   90.99580_gp,&
         69.92890_gp,   53.67380_gp,   41.73100_gp,  525.14830_gp,  418.47680_gp,&
         391.67000_gp,  315.14880_gp,  250.44750_gp,  210.61910_gp,  174.29820_gp,&
         144.28140_gp,  858.35650_gp,  726.11640_gp,  603.56880_gp,  586.85490_gp,&
         539.50230_gp,  426.36910_gp,  467.16530_gp,  368.13400_gp,  392.25490_gp,&
         402.36920_gp,  309.64490_gp,  320.85180_gp,  378.20130_gp,  338.70950_gp,&
         293.30960_gp,  265.87000_gp,  235.16470_gp,  206.59450_gp,  965.46080_gp,&
         865.21550_gp,  766.64930_gp,  694.31570_gp,  636.59010_gp,  496.73190_gp,&
         551.96750_gp,  425.41540_gp,  464.08940_gp,  431.83300_gp,  360.56290_gp/)
    vdwparams%coeffs(24801:24900)=(/  382.21490_gp,  474.34870_gp,  443.70850_gp,  400.39260_gp,  374.50420_gp,&
         341.68440_gp,  309.92280_gp, 1177.86220_gp, 1100.21490_gp,  977.35910_gp,&
         472.49460_gp,  978.80490_gp,  941.44390_gp,  918.23780_gp,  896.83560_gp,&
         877.87040_gp,  700.22530_gp,  775.37110_gp,  749.50170_gp,  793.80070_gp,&
         777.10890_gp,  762.18580_gp,  752.78270_gp,  640.40610_gp,  636.88250_gp,&
         585.84610_gp,  498.67750_gp,  508.52330_gp,  464.13690_gp,  427.56890_gp,&
         357.97140_gp,  335.46850_gp,  345.58560_gp,  489.85300_gp,  482.77210_gp,&
         448.34220_gp,  430.41860_gp,  400.52020_gp,  370.00240_gp, 1122.33710_gp,&
         1094.71710_gp,  977.40400_gp,  889.77500_gp,  878.87760_gp,  851.33500_gp,&
         869.82820_gp,  843.44130_gp,   49.86160_gp,  154.70980_gp,  198.67640_gp,&
         155.94320_gp,  119.71100_gp,   84.99060_gp,   63.00320_gp,   44.29330_gp,&
         225.33300_gp,  347.51630_gp,  356.39850_gp,  292.29100_gp,  243.32340_gp,&
         208.27540_gp,  172.47660_gp,  315.72770_gp,  583.37200_gp,  319.98810_gp,&
         309.28520_gp,  303.06940_gp,  299.74820_gp,  278.13060_gp,  258.24560_gp,&
         246.47110_gp,  240.79380_gp,  236.33890_gp,  224.42270_gp,  356.65940_gp,&
         318.14490_gp,  288.62680_gp,  265.94970_gp,  236.68530_gp,  378.38750_gp,&
         709.63970_gp,  555.20870_gp,  423.30460_gp,  427.40240_gp,  400.15670_gp,&
         443.82790_gp,  352.02000_gp,  329.85620_gp,  307.03010_gp,  296.64300_gp,&
         295.10960_gp,  454.11220_gp,  415.92040_gp,  393.10920_gp,  374.41140_gp,&
         344.82440_gp,  451.52920_gp,  917.08910_gp,  701.71880_gp,  458.62530_gp/)
    vdwparams%coeffs(24901:25000)=(/  449.20880_gp,  435.28510_gp,  436.62810_gp,  418.57700_gp,  438.96330_gp,&
         413.03190_gp,  418.95270_gp,  392.63770_gp,  381.67800_gp,  379.01540_gp,&
         397.06740_gp,  366.63560_gp,  486.47620_gp,  454.99830_gp,  418.51340_gp,&
         420.70230_gp,  374.94930_gp,  353.92940_gp,  338.77130_gp,  324.04200_gp,&
         314.91790_gp,  489.72910_gp,  451.62020_gp,  440.81790_gp,  430.91310_gp,&
         405.35970_gp,  508.82630_gp,  918.74100_gp,  542.34140_gp,  619.96910_gp,&
         557.04040_gp,  496.09310_gp,  478.31910_gp,  563.13100_gp,  133.19340_gp,&
         131.89570_gp,   98.60190_gp,   79.97900_gp,   55.89460_gp,  233.03610_gp,&
         281.37080_gp,  272.72890_gp,  239.83560_gp,  204.14550_gp,  327.70190_gp,&
         317.37170_gp,  321.50860_gp,  294.39580_gp,  225.19270_gp,  194.39360_gp,&
         190.62290_gp,  220.28000_gp,  206.01270_gp,  286.65620_gp,  303.71450_gp,&
         281.50840_gp,  263.95060_gp,  393.46440_gp,  382.37800_gp,  386.83830_gp,&
         374.17790_gp,  334.05520_gp,  297.16160_gp,  281.66080_gp,  284.75480_gp,&
         296.56750_gp,  371.25760_gp,  403.16950_gp,  383.18380_gp,  371.69830_gp,&
         477.55530_gp,  497.02050_gp,  372.00450_gp,  387.48710_gp,  363.09850_gp,&
         332.14480_gp,   29.95970_gp,   20.73210_gp,  399.46230_gp,  245.60490_gp,&
         172.42600_gp,  120.47520_gp,   86.51500_gp,   66.79830_gp,   51.48390_gp,&
         40.16660_gp,  480.02060_gp,  387.03270_gp,  364.08530_gp,  294.97480_gp,&
         235.85440_gp,  199.18590_gp,  165.51840_gp,  137.52330_gp,  783.98180_gp,&
         668.74580_gp,  557.03780_gp,  542.87970_gp,  499.73800_gp,  395.39650_gp/)
    vdwparams%coeffs(25001:25100)=(/  433.57460_gp,  342.14260_gp,  365.17850_gp,  374.07150_gp,  288.25030_gp,&
         299.57880_gp,  352.52400_gp,  317.31180_gp,  276.13370_gp,  251.11640_gp,&
         222.87460_gp,  196.43010_gp,  882.92900_gp,  796.89500_gp,  708.87790_gp,&
         643.66450_gp,  591.24780_gp,  463.02150_gp,  513.80610_gp,  397.55570_gp,&
         433.40920_gp,  403.74870_gp,  337.29140_gp,  357.90130_gp,  442.52210_gp,&
         415.50540_gp,  376.47850_gp,  353.06060_gp,  323.05850_gp,  293.86160_gp,&
         1077.71580_gp, 1012.46430_gp,  902.81190_gp,  444.10240_gp,  902.05950_gp,&
         868.10970_gp,  846.83980_gp,  827.20360_gp,  809.80750_gp,  649.09210_gp,&
         714.73580_gp,  691.43160_gp,  732.88920_gp,  717.52640_gp,  703.83170_gp,&
         694.99350_gp,  593.20990_gp,  592.12750_gp,  546.13100_gp,  466.12960_gp,&
         475.70880_gp,  435.20440_gp,  401.69710_gp,  337.18850_gp,  316.31490_gp,&
         326.00310_gp,  457.88770_gp,  452.44550_gp,  421.63980_gp,  405.64640_gp,&
         378.47640_gp,  350.52560_gp, 1030.14310_gp, 1009.54850_gp,  904.47790_gp,&
         827.02680_gp,  815.50490_gp,  790.06160_gp,  805.13870_gp,  781.06340_gp,&
         46.97930_gp,  143.71930_gp,  184.97380_gp,  146.47120_gp,  113.11480_gp,&
         80.84140_gp,   60.26260_gp,   42.66360_gp,  209.09190_gp,  322.13980_gp,&
         331.74520_gp,  273.89260_gp,  229.21310_gp,  196.96310_gp,  163.78500_gp,&
         294.55150_gp,  538.11810_gp,  299.79870_gp,  289.89590_gp,  284.06410_gp,&
         280.79860_gp,  261.25710_gp,  242.83990_gp,  231.80890_gp,  226.42260_gp,&
         221.76420_gp,  211.39390_gp,  332.89310_gp,  298.43870_gp,  271.81450_gp/)
    vdwparams%coeffs(25101:25200)=(/  251.18530_gp,  224.28140_gp,  353.52250_gp,  654.52330_gp,  516.11830_gp,&
         396.62190_gp,  400.43580_gp,  375.41960_gp,  414.89900_gp,  331.25900_gp,&
         310.57580_gp,  289.31700_gp,  279.42790_gp,  278.48020_gp,  424.25050_gp,&
         390.14410_gp,  369.84470_gp,  353.02740_gp,  326.00590_gp,  423.07320_gp,&
         844.93410_gp,  651.99750_gp,  431.18430_gp,  422.35310_gp,  409.35960_gp,&
         410.38090_gp,  392.71700_gp,  412.44200_gp,  388.30810_gp,  393.47340_gp,&
         369.31690_gp,  359.07870_gp,  356.50780_gp,  373.07840_gp,  344.95150_gp,&
         454.55470_gp,  426.23560_gp,  393.10340_gp,  394.31360_gp,  353.50380_gp,&
         334.05950_gp,  319.98180_gp,  306.13090_gp,  298.06110_gp,  458.19840_gp,&
         424.17020_gp,  414.89180_gp,  406.23070_gp,  383.04430_gp,  476.84810_gp,&
         848.51510_gp,  509.77850_gp,  580.92120_gp,  522.82120_gp,  466.79740_gp,&
         450.36480_gp,  526.48950_gp,  125.10340_gp,  124.18210_gp,   93.47380_gp,&
         76.13650_gp,   53.56410_gp,  217.93200_gp,  263.04720_gp,  255.82610_gp,&
         225.87090_gp,  193.08290_gp,  306.76340_gp,  297.84850_gp,  301.80320_gp,&
         276.40530_gp,  212.49880_gp,  183.96380_gp,  180.34190_gp,  207.43670_gp,&
         194.24620_gp,  268.82090_gp,  285.15730_gp,  265.25130_gp,  249.32490_gp,&
         368.58100_gp,  359.32890_gp,  363.88390_gp,  352.00850_gp,  315.25780_gp,&
         281.15270_gp,  266.75220_gp,  269.04190_gp,  279.93090_gp,  348.65850_gp,&
         378.51600_gp,  360.77560_gp,  350.56820_gp,  447.86270_gp,  466.73680_gp,&
         350.44320_gp,  364.70460_gp,  342.72390_gp,  314.23420_gp,  297.83380_gp/)
    vdwparams%coeffs(25201:25300)=(/   30.22940_gp,   20.76080_gp,  426.36380_gp,  254.74750_gp,  176.35010_gp,&
         122.10880_gp,   87.19070_gp,   67.09130_gp,   51.58500_gp,   40.18460_gp,&
         511.43160_gp,  403.39630_gp,  376.53390_gp,  302.09560_gp,  239.74340_gp,&
         201.58550_gp,  166.87220_gp,  138.23710_gp,  838.84470_gp,  703.06100_gp,&
         583.40320_gp,  566.63000_gp,  520.55850_gp,  411.71570_gp,  450.35590_gp,&
         355.13960_gp,  377.53770_gp,  387.50230_gp,  298.55250_gp,  308.44330_gp,&
         363.54100_gp,  324.81010_gp,  280.90060_gp,  254.53380_gp,  225.13410_gp,&
         197.84740_gp,  943.12750_gp,  838.25760_gp,  740.70150_gp,  669.84280_gp,&
         613.69150_gp,  478.50440_gp,  531.85130_gp,  409.61050_gp,  446.56770_gp,&
         415.38870_gp,  347.36190_gp,  367.54320_gp,  456.59440_gp,  426.12630_gp,&
         383.93930_gp,  358.89370_gp,  327.31610_gp,  296.86250_gp, 1150.66530_gp,&
         1067.30480_gp,  945.34330_gp,  453.35410_gp,  948.86600_gp,  912.14940_gp,&
         889.52970_gp,  868.67840_gp,  850.19490_gp,  676.50560_gp,  752.86050_gp,&
         727.36460_gp,  768.17290_gp,  751.93870_gp,  737.39820_gp,  728.36400_gp,&
         618.61380_gp,  613.42840_gp,  563.63530_gp,  479.69460_gp,  488.87240_gp,&
         445.89270_gp,  410.58900_gp,  343.85250_gp,  322.28400_gp,  331.78240_gp,&
         471.74830_gp,  463.98210_gp,  430.22560_gp,  412.74680_gp,  383.86920_gp,&
         354.53980_gp, 1093.63530_gp, 1060.46550_gp,  944.39170_gp,  857.75700_gp,&
         848.61880_gp,  822.01600_gp,  841.27950_gp,  815.49360_gp,   47.71070_gp,&
         148.77410_gp,  190.95480_gp,  149.41700_gp,  114.67490_gp,   81.46860_gp/)
    vdwparams%coeffs(25301:25400)=(/   60.47490_gp,   42.63300_gp,  217.06070_gp,  334.78560_gp,  342.47320_gp,&
         280.13350_gp,  232.94100_gp,  199.36180_gp,  165.14130_gp,  304.13400_gp,&
         565.34700_gp,  307.20500_gp,  296.99590_gp,  291.07410_gp,  288.00240_gp,&
         266.80930_gp,  247.70870_gp,  236.46980_gp,  231.06240_gp,  227.06590_gp,&
         215.17070_gp,  342.64070_gp,  304.98940_gp,  276.41200_gp,  254.62020_gp,&
         226.59700_gp,  364.41220_gp,  688.23470_gp,  535.70340_gp,  406.65130_gp,&
         410.63670_gp,  384.41940_gp,  427.30450_gp,  337.87470_gp,  316.68990_gp,&
         294.82600_gp,  284.94240_gp,  283.09650_gp,  436.80210_gp,  399.22460_gp,&
         376.89830_gp,  358.78610_gp,  330.31590_gp,  433.96820_gp,  890.82740_gp,&
         677.44550_gp,  440.03130_gp,  430.99700_gp,  417.62690_gp,  419.01960_gp,&
         402.12590_gp,  421.22030_gp,  396.30640_gp,  402.19400_gp,  376.65190_gp,&
         366.11320_gp,  363.56690_gp,  380.89510_gp,  351.63250_gp,  467.99190_gp,&
         437.30900_gp,  401.91040_gp,  404.42320_gp,  359.70450_gp,  339.57300_gp,&
         325.07870_gp,  311.08440_gp,  302.06450_gp,  471.25110_gp,  433.67490_gp,&
         422.87330_gp,  413.15870_gp,  388.47900_gp,  488.93950_gp,  890.67630_gp,&
         520.29480_gp,  595.99730_gp,  535.52370_gp,  476.21600_gp,  459.07020_gp,&
         542.51670_gp,  127.50670_gp,  126.36580_gp,   94.44220_gp,   76.69360_gp,&
         53.71350_gp,  223.45060_gp,  269.92960_gp,  261.31790_gp,  229.65740_gp,&
         195.42290_gp,  314.65860_gp,  304.40560_gp,  308.37100_gp,  282.52780_gp,&
         216.07310_gp,  186.47130_gp,  182.88010_gp,  211.45410_gp,  197.70250_gp/)
    vdwparams%coeffs(25401:25500)=(/  274.93230_gp,  291.09370_gp,  269.57470_gp,  252.71240_gp,  377.90170_gp,&
         366.67450_gp,  370.83430_gp,  358.92720_gp,  320.31090_gp,  284.89430_gp,&
         270.05380_gp,  273.28420_gp,  284.60940_gp,  356.43230_gp,  386.87350_gp,&
         367.31080_gp,  356.16120_gp,  458.58100_gp,  476.69270_gp,  356.51190_gp,&
         371.70810_gp,  348.12240_gp,  318.41110_gp,  301.21160_gp,  305.52580_gp,&
         31.73680_gp,   21.70640_gp,  448.17370_gp,  268.64910_gp,  185.71950_gp,&
         128.32810_gp,   91.44410_gp,   70.24720_gp,   53.92930_gp,   41.95740_gp,&
         537.64280_gp,  425.31660_gp,  396.91230_gp,  318.22850_gp,  252.21940_gp,&
         211.81210_gp,  175.10170_gp,  144.86250_gp,  879.76340_gp,  740.19810_gp,&
         614.51280_gp,  596.80090_gp,  548.28530_gp,  433.31070_gp,  474.33180_gp,&
         373.75740_gp,  397.65430_gp,  408.19710_gp,  314.18780_gp,  324.82660_gp,&
         382.93500_gp,  342.02100_gp,  295.51180_gp,  267.53780_gp,  236.39130_gp,&
         207.51570_gp,  988.97460_gp,  882.12520_gp,  779.87150_gp,  705.32520_gp,&
         646.12000_gp,  503.47710_gp,  559.75880_gp,  430.81480_gp,  469.96810_gp,&
         437.11170_gp,  365.18920_gp,  386.69930_gp,  480.58480_gp,  448.55300_gp,&
         403.93700_gp,  377.37610_gp,  343.91650_gp,  311.65760_gp, 1206.09070_gp,&
         1122.40060_gp,  994.84110_gp,  476.85250_gp,  997.99760_gp,  959.58830_gp,&
         935.85100_gp,  913.96730_gp,  894.57000_gp,  711.80730_gp,  790.87860_gp,&
         764.13790_gp,  808.48440_gp,  791.44150_gp,  776.17700_gp,  766.68210_gp,&
         651.13390_gp,  645.97680_gp,  593.45200_gp,  504.69500_gp,  514.41860_gp/)
    vdwparams%coeffs(25501:25600)=(/  469.06170_gp,  431.79110_gp,  361.31080_gp,  338.54780_gp,  348.60350_gp,&
         496.10630_gp,  488.12730_gp,  452.49450_gp,  433.95750_gp,  403.35980_gp,&
         372.28610_gp, 1147.18090_gp, 1115.45720_gp,  994.04090_gp,  902.93930_gp,&
         892.95040_gp,  864.95490_gp,  885.11780_gp,  858.05490_gp,   50.18990_gp,&
         156.83270_gp,  201.16870_gp,  157.24660_gp,  120.47330_gp,   85.41130_gp,&
         63.28030_gp,   44.50040_gp,  228.73580_gp,  352.78900_gp,  360.93540_gp,&
         295.00710_gp,  245.03090_gp,  209.47030_gp,  173.28290_gp,  319.91700_gp,&
         594.42400_gp,  323.34130_gp,  312.51320_gp,  306.27330_gp,  303.03740_gp,&
         280.75180_gp,  260.58650_gp,  248.72640_gp,  243.04070_gp,  238.86500_gp,&
         226.30430_gp,  360.87970_gp,  321.06640_gp,  290.76020_gp,  267.62610_gp,&
         237.93600_gp,  383.18340_gp,  723.24770_gp,  563.37460_gp,  427.78570_gp,&
         431.99420_gp,  404.27790_gp,  449.29790_gp,  355.22280_gp,  332.86530_gp,&
         309.80620_gp,  299.43630_gp,  297.52810_gp,  459.73620_gp,  420.13030_gp,&
         396.48380_gp,  377.25130_gp,  347.07680_gp,  456.31660_gp,  935.34640_gp,&
         712.22580_gp,  462.81270_gp,  453.30770_gp,  439.22170_gp,  440.70090_gp,&
         422.96520_gp,  443.07970_gp,  416.83270_gp,  423.02650_gp,  396.14680_gp,&
         385.05430_gp,  382.39220_gp,  400.73810_gp,  369.83970_gp,  492.31890_gp,&
         459.93990_gp,  422.59020_gp,  425.29100_gp,  378.03640_gp,  356.76720_gp,&
         341.46440_gp,  326.70140_gp,  317.19220_gp,  495.66950_gp,  456.14500_gp,&
         444.71840_gp,  434.37700_gp,  408.21270_gp,  514.13610_gp,  935.70360_gp/)
    vdwparams%coeffs(25601:25700)=(/  547.28610_gp,  626.74370_gp,  562.83810_gp,  500.66010_gp,  482.60870_gp,&
         570.30020_gp,  134.24720_gp,  132.88680_gp,   99.14260_gp,   80.38780_gp,&
         56.16470_gp,  235.46090_gp,  284.33280_gp,  275.14210_gp,  241.57090_gp,&
         205.31920_gp,  331.24320_gp,  320.36790_gp,  324.52330_gp,  297.23410_gp,&
         227.02370_gp,  195.79190_gp,  192.03650_gp,  222.28680_gp,  207.79350_gp,&
         289.39310_gp,  306.38780_gp,  283.53080_gp,  265.60990_gp,  397.65110_gp,&
         385.77200_gp,  390.08630_gp,  377.42800_gp,  336.58580_gp,  299.17850_gp,&
         283.51650_gp,  287.00770_gp,  299.01090_gp,  374.93950_gp,  407.08260_gp,&
         386.34370_gp,  374.46930_gp,  482.36480_gp,  501.53880_gp,  374.91100_gp,&
         390.79260_gp,  365.80820_gp,  334.38680_gp,  316.17560_gp,  320.76260_gp,&
         336.87270_gp,   37.64070_gp,   25.08750_gp,  584.87960_gp,  336.72920_gp,&
         226.93930_gp,  153.70200_gp,  107.86020_gp,   81.95470_gp,   62.33170_gp,&
         48.13680_gp,  699.57010_gp,  536.90060_gp,  494.77230_gp,  389.92290_gp,&
         304.37830_gp,  252.98030_gp,  207.04280_gp,  169.77090_gp, 1148.83070_gp,&
         944.94750_gp,  780.31190_gp,  753.56620_gp,  690.06660_gp,  544.17950_gp,&
         594.17250_gp,  466.87620_gp,  494.28940_gp,  509.15080_gp,  390.94350_gp,&
         400.86460_gp,  474.36550_gp,  418.28000_gp,  356.98080_gp,  320.59210_gp,&
         280.90020_gp,  244.65900_gp, 1287.78260_gp, 1126.33910_gp,  985.98710_gp,&
         886.00180_gp,  807.98210_gp,  624.23380_gp,  696.29870_gp,  530.92790_gp,&
         579.90290_gp,  537.85540_gp,  449.15460_gp,  474.10150_gp,  594.60280_gp/)
    vdwparams%coeffs(25701:25800)=(/  549.54150_gp,  489.78390_gp,  454.57430_gp,  411.27480_gp,  370.08590_gp,&
         1569.07570_gp, 1436.76820_gp, 1261.21740_gp,  578.83370_gp, 1273.12100_gp,&
         1222.25930_gp, 1191.53650_gp, 1163.28070_gp, 1138.21890_gp,  894.97400_gp,&
         1009.38160_gp,  973.40780_gp, 1026.34200_gp, 1004.50560_gp,  984.80680_gp,&
         973.27760_gp,  819.98780_gp,  805.59300_gp,  735.21850_gp,  621.34250_gp,&
         631.98300_gp,  572.93310_gp,  524.89900_gp,  436.57950_gp,  408.12070_gp,&
         419.66530_gp,  611.17330_gp,  597.12280_gp,  548.62180_gp,  523.30660_gp,&
         483.15870_gp,  443.13040_gp, 1480.48890_gp, 1420.29560_gp, 1254.62790_gp,&
         1127.25700_gp, 1119.95680_gp, 1084.47850_gp, 1117.05540_gp, 1081.66120_gp,&
         60.47920_gp,  195.81130_gp,  249.84400_gp,  191.03520_gp,  144.23770_gp,&
         100.63880_gp,   73.58410_gp,   50.93070_gp,  286.55560_gp,  442.98480_gp,&
         448.45700_gp,  360.47340_gp,  295.50270_gp,  250.21530_gp,  204.91230_gp,&
         395.86240_gp,  756.78810_gp,  395.55800_gp,  381.99770_gp,  374.41990_gp,&
         371.01600_gp,  341.29280_gp,  315.96540_gp,  301.49960_gp,  294.78150_gp,&
         291.34120_gp,  273.25060_gp,  445.54260_gp,  391.38140_gp,  350.96110_gp,&
         320.72380_gp,  282.84750_gp,  472.50610_gp,  921.36230_gp,  703.55580_gp,&
         523.40250_gp,  528.63200_gp,  493.14640_gp,  553.29240_gp,  430.04960_gp,&
         402.51430_gp,  373.95490_gp,  361.82100_gp,  357.76620_gp,  566.74250_gp,&
         512.55600_gp,  480.04840_gp,  454.24350_gp,  415.11430_gp,  558.34700_gp,&
         1195.48060_gp,  890.67110_gp,  561.41900_gp,  549.82120_gp,  532.42950_gp/)
    vdwparams%coeffs(25801:25900)=(/  535.04270_gp,  515.92290_gp,  538.31350_gp,  505.69480_gp,  514.58640_gp,&
         479.98260_gp,  466.31110_gp,  463.29610_gp,  486.80100_gp,  447.77810_gp,&
         606.70050_gp,  563.15480_gp,  513.97290_gp,  520.14190_gp,  455.49230_gp,&
         428.73640_gp,  409.67910_gp,  391.87930_gp,  378.62560_gp,  609.07890_gp,&
         554.94180_gp,  538.10030_gp,  523.40750_gp,  488.96760_gp,  628.76220_gp,&
         1188.34060_gp,  664.19980_gp,  767.13680_gp,  686.33050_gp,  606.32050_gp,&
         583.50520_gp,  702.42340_gp,  163.03250_gp,  160.51470_gp,  117.75890_gp,&
         94.56190_gp,   65.05490_gp,  289.15800_gp,  349.50920_gp,  335.34010_gp,&
         291.53890_gp,  245.18670_gp,  406.06400_gp,  390.20820_gp,  395.04660_gp,&
         361.78060_gp,  272.98940_gp,  233.78110_gp,  229.48940_gp,  268.71310_gp,&
         250.40430_gp,  353.11120_gp,  372.66030_gp,  341.78700_gp,  318.22560_gp,&
         486.70590_gp,  468.32670_gp,  472.37120_gp,  457.08930_gp,  404.45180_gp,&
         357.27900_gp,  337.77820_gp,  344.11230_gp,  359.34930_gp,  456.22010_gp,&
         495.54180_gp,  466.89650_gp,  450.57090_gp,  588.74480_gp,  609.88230_gp,&
         452.35490_gp,  472.75440_gp,  439.36980_gp,  399.36360_gp,  375.91510_gp,&
         383.38850_gp,  403.05640_gp,  487.74180_gp,   39.79220_gp,   26.44670_gp,&
         614.88410_gp,  355.82050_gp,  240.02620_gp,  162.50360_gp,  113.92620_gp,&
         86.47120_gp,   65.68880_gp,   50.67070_gp,  735.57820_gp,  566.88720_gp,&
         522.87510_gp,  412.41080_gp,  321.96130_gp,  267.49610_gp,  218.79400_gp,&
         179.27010_gp, 1205.86720_gp,  995.78930_gp,  822.87620_gp,  794.94560_gp/)
    vdwparams%coeffs(25901:26000)=(/  728.13060_gp,  573.90570_gp,  627.13920_gp,  492.53590_gp,  522.01020_gp,&
         537.61470_gp,  412.50160_gp,  423.50910_gp,  501.21030_gp,  442.28050_gp,&
         377.53680_gp,  338.99500_gp,  296.91970_gp,  258.48000_gp, 1351.80890_gp,&
         1186.52980_gp, 1039.78740_gp,  934.82320_gp,  852.68980_gp,  658.80460_gp,&
         734.87030_gp,  560.35040_gp,  612.29710_gp,  567.95070_gp,  473.88930_gp,&
         500.66290_gp,  627.79210_gp,  580.70530_gp,  517.76630_gp,  480.56780_gp,&
         434.74640_gp,  391.10730_gp, 1646.90530_gp, 1512.63160_gp, 1329.34370_gp,&
         611.69810_gp, 1340.70870_gp, 1287.44260_gp, 1255.16510_gp, 1225.47470_gp,&
         1199.14430_gp,  943.62380_gp, 1062.05890_gp, 1024.40690_gp, 1081.63820_gp,&
         1058.68020_gp, 1037.98140_gp, 1025.80780_gp,  864.71680_gp,  850.48450_gp,&
         776.45220_gp,  656.06950_gp,  667.47630_gp,  605.20260_gp,  554.49180_gp,&
         461.01730_gp,  430.90160_gp,  443.23910_gp,  644.99260_gp,  630.69680_gp,&
         579.74990_gp,  553.07510_gp,  510.64800_gp,  468.27350_gp, 1555.42130_gp,&
         1496.05600_gp, 1322.94520_gp, 1189.59300_gp, 1181.08950_gp, 1143.67040_gp,&
         1177.31660_gp, 1140.16260_gp,   63.97960_gp,  206.91270_gp,  264.00830_gp,&
         202.01910_gp,  152.44950_gp,  106.26350_gp,   77.60240_gp,   53.61090_gp,&
         302.58820_gp,  467.75760_gp,  473.97540_gp,  381.25380_gp,  312.55060_gp,&
         264.56000_gp,  216.53480_gp,  417.74830_gp,  797.00090_gp,  418.01280_gp,&
         403.62390_gp,  395.59060_gp,  391.93640_gp,  360.75390_gp,  333.96980_gp,&
         318.63900_gp,  311.52110_gp,  307.75700_gp,  288.85340_gp,  470.83440_gp/)
    vdwparams%coeffs(26001:26100)=(/  413.85520_gp,  371.16030_gp,  339.12820_gp,  298.97830_gp,  498.60580_gp,&
         969.92320_gp,  742.05020_gp,  552.89600_gp,  558.39300_gp,  520.87780_gp,&
         583.93470_gp,  454.33520_gp,  425.16750_gp,  394.94400_gp,  382.09250_gp,&
         378.01500_gp,  598.52070_gp,  541.68720_gp,  507.48490_gp,  480.22270_gp,&
         438.81240_gp,  589.60800_gp, 1257.56930_gp,  939.14200_gp,  593.29090_gp,&
         581.03250_gp,  562.64920_gp,  565.36800_gp,  544.95840_gp,  568.86990_gp,&
         534.39680_gp,  543.69520_gp,  507.26200_gp,  492.82300_gp,  489.63950_gp,&
         514.52780_gp,  473.26870_gp,  640.58910_gp,  594.75240_gp,  542.91370_gp,&
         549.25780_gp,  481.24530_gp,  452.91160_gp,  432.72460_gp,  413.82910_gp,&
         399.95160_gp,  642.99800_gp,  586.28560_gp,  568.68040_gp,  553.20800_gp,&
         516.80550_gp,  664.02220_gp, 1251.01720_gp,  701.95930_gp,  810.10680_gp,&
         724.64740_gp,  640.54440_gp,  616.46960_gp,  741.11760_gp,  172.49040_gp,&
         169.70770_gp,  124.43980_gp,   99.82940_gp,   68.56320_gp,  305.83900_gp,&
         369.56280_gp,  354.68400_gp,  308.32680_gp,  259.23020_gp,  429.11330_gp,&
         412.49460_gp,  417.60390_gp,  382.32410_gp,  288.38020_gp,  246.93330_gp,&
         242.39630_gp,  283.87680_gp,  264.54680_gp,  373.31260_gp,  394.06900_gp,&
         361.45320_gp,  336.47870_gp,  514.20470_gp,  495.05130_gp,  499.36110_gp,&
         483.03860_gp,  427.36900_gp,  377.46330_gp,  356.81920_gp,  363.42680_gp,&
         379.57410_gp,  482.04630_gp,  523.74110_gp,  493.59530_gp,  476.34340_gp,&
         621.97270_gp,  644.64210_gp,  478.20960_gp,  499.56010_gp,  464.28480_gp/)
    vdwparams%coeffs(26101:26200)=(/  421.94200_gp,  397.11980_gp,  404.88970_gp,  425.75640_gp,  515.25750_gp,&
         544.45630_gp,   40.76490_gp,   27.17460_gp,  613.18280_gp,  359.84820_gp,&
         244.35680_gp,  166.11890_gp,  116.75030_gp,   88.73970_gp,   67.47500_gp,&
         52.07590_gp,  734.12320_gp,  571.93110_gp,  529.49770_gp,  419.59890_gp,&
         328.72910_gp,  273.65740_gp,  224.21020_gp,  183.93680_gp, 1200.40500_gp,&
         1000.18150_gp,  828.06580_gp,  801.23380_gp,  734.59090_gp,  578.97210_gp,&
         633.56300_gp,  497.63260_gp,  528.56480_gp,  543.85410_gp,  417.19990_gp,&
         429.68510_gp,  508.15670_gp,  450.04640_gp,  385.31250_gp,  346.54960_gp,&
         304.00640_gp,  264.98470_gp, 1346.60510_gp, 1191.29730_gp, 1047.36940_gp,&
         943.46900_gp,  861.62430_gp,  666.99900_gp,  743.48240_gp,  568.09890_gp,&
         620.85810_gp,  576.30050_gp,  480.48740_gp,  508.46720_gp,  636.17130_gp,&
         590.26500_gp,  527.74110_gp,  490.58420_gp,  444.49370_gp,  400.41400_gp,&
         1640.71710_gp, 1516.95260_gp, 1337.53620_gp,  623.10240_gp, 1345.89840_gp,&
         1293.16560_gp, 1260.94140_gp, 1231.28010_gp, 1204.98330_gp,  951.46290_gp,&
         1065.14600_gp, 1028.00010_gp, 1087.81120_gp, 1064.82440_gp, 1044.14470_gp,&
         1031.75600_gp,  871.74390_gp,  860.24520_gp,  786.78670_gp,  665.60880_gp,&
         677.67670_gp,  615.34220_gp,  564.41430_gp,  469.71680_gp,  439.19500_gp,&
         452.07010_gp,  653.91820_gp,  640.98700_gp,  590.69190_gp,  564.28990_gp,&
         521.80500_gp,  479.12360_gp, 1553.83860_gp, 1502.86270_gp, 1332.94660_gp,&
         1202.46500_gp, 1191.81720_gp, 1154.13380_gp, 1185.62950_gp, 1148.65540_gp/)
    vdwparams%coeffs(26201:26300)=(/   65.35360_gp,  209.43720_gp,  267.54910_gp,  205.89850_gp,  155.80320_gp,&
         108.88470_gp,   79.65760_gp,   55.12830_gp,  305.82180_gp,  472.51580_gp,&
         480.37720_gp,  388.13130_gp,  319.14590_gp,  270.63230_gp,  221.88150_gp,&
         423.09090_gp,  800.50250_gp,  425.05990_gp,  410.44180_gp,  402.23240_gp,&
         398.32240_gp,  367.44750_gp,  340.35150_gp,  324.70530_gp,  317.39060_gp,&
         313.02130_gp,  294.68800_gp,  477.79970_gp,  421.44820_gp,  378.86260_gp,&
         346.67280_gp,  306.08650_gp,  505.34960_gp,  973.64740_gp,  749.62580_gp,&
         561.97920_gp,  567.51500_gp,  529.72340_gp,  592.10620_gp,  462.93930_gp,&
         433.25010_gp,  402.55460_gp,  389.32060_gp,  385.81270_gp,  607.15180_gp,&
         551.19600_gp,  517.44360_gp,  490.28150_gp,  448.64430_gp,  599.07190_gp,&
         1260.55370_gp,  948.16680_gp,  604.43760_gp,  591.96230_gp,  573.30320_gp,&
         575.83010_gp,  554.23280_gp,  579.35270_gp,  544.41580_gp,  553.45320_gp,&
         516.96070_gp,  502.31200_gp,  499.01960_gp,  524.12520_gp,  482.44240_gp,&
         649.75910_gp,  604.28050_gp,  552.54180_gp,  558.11920_gp,  490.92230_gp,&
         462.22330_gp,  441.72420_gp,  422.35490_gp,  408.77310_gp,  652.49230_gp,&
         596.73810_gp,  579.73290_gp,  564.55350_gp,  528.11830_gp,  674.82210_gp,&
         1256.72780_gp,  715.13070_gp,  823.07590_gp,  736.72690_gp,  652.61520_gp,&
         628.34540_gp,  751.23490_gp,  175.90820_gp,  173.15440_gp,  127.38320_gp,&
         102.31260_gp,   70.39610_gp,  311.01430_gp,  375.62900_gp,  361.30250_gp,&
         314.74790_gp,  265.18320_gp,  436.13820_gp,  420.00450_gp,  425.25490_gp/)
    vdwparams%coeffs(26301:26400)=(/  389.22030_gp,  294.26590_gp,  252.34820_gp,  247.65980_gp,  289.33940_gp,&
         269.84110_gp,  380.00030_gp,  401.50980_gp,  369.05850_gp,  343.97610_gp,&
         522.67140_gp,  504.41040_gp,  509.13000_gp,  492.30130_gp,  436.28420_gp,&
         385.81740_gp,  364.86800_gp,  370.99700_gp,  387.33000_gp,  490.68810_gp,&
         533.23730_gp,  503.51030_gp,  486.40400_gp,  632.57940_gp,  656.53900_gp,&
         487.97410_gp,  509.23540_gp,  474.03980_gp,  431.28970_gp,  406.28510_gp,&
         413.59970_gp,  434.91300_gp,  525.03950_gp,  554.92310_gp,  566.04970_gp,&
         40.68350_gp,   27.26900_gp,  597.04830_gp,  354.12790_gp,  242.07220_gp,&
         165.37510_gp,  116.63660_gp,   88.86410_gp,   67.69920_gp,   52.32370_gp,&
         715.29180_gp,  561.73090_gp,  521.82740_gp,  415.39360_gp,  326.69480_gp,&
         272.64470_gp,  223.90350_gp,  184.04790_gp, 1168.56360_gp,  979.28960_gp,&
         811.91130_gp,  786.75610_gp,  721.92550_gp,  569.23970_gp,  623.40350_gp,&
         489.93540_gp,  521.13280_gp,  535.72570_gp,  411.14170_gp,  424.43520_gp,&
         501.54410_gp,  445.70000_gp,  382.80030_gp,  344.97610_gp,  303.23440_gp,&
         264.79090_gp, 1311.83790_gp, 1166.32700_gp, 1028.08790_gp,  927.65740_gp,&
         848.16650_gp,  657.96520_gp,  732.82990_gp,  561.23130_gp,  613.20640_gp,&
         569.59290_gp,  474.85530_gp,  502.99460_gp,  627.93100_gp,  584.16010_gp,&
         523.70500_gp,  487.65620_gp,  442.64170_gp,  399.42590_gp, 1598.90910_gp,&
         1484.20300_gp, 1311.98050_gp,  618.10270_gp, 1317.95990_gp, 1266.79560_gp,&
         1235.35570_gp, 1206.40010_gp, 1180.73460_gp,  935.18300_gp, 1042.90690_gp/)
    vdwparams%coeffs(26401:26500)=(/ 1007.06570_gp, 1066.54710_gp, 1044.06540_gp, 1023.87930_gp, 1011.59410_gp,&
         856.51120_gp,  847.44200_gp,  776.39170_gp,  657.82210_gp,  670.12730_gp,&
         609.37170_gp,  559.59870_gp,  466.35070_gp,  436.27520_gp,  449.26040_gp,&
         646.07450_gp,  634.50740_gp,  586.10180_gp,  560.69830_gp,  519.36900_gp,&
         477.63150_gp, 1517.31400_gp, 1472.45100_gp, 1308.95970_gp, 1184.16440_gp,&
         1172.18950_gp, 1135.19890_gp, 1164.11930_gp, 1128.15200_gp,   64.97600_gp,&
         206.31430_gp,  263.93510_gp,  204.26390_gp,  155.10600_gp,  108.79270_gp,&
         79.81880_gp,   55.42640_gp,  300.93900_gp,  464.75070_gp,  473.82140_gp,&
         384.51100_gp,  317.21920_gp,  269.61720_gp,  221.56890_gp,  417.56030_gp,&
         784.39170_gp,  420.82330_gp,  406.43450_gp,  398.27910_gp,  394.24420_gp,&
         364.38140_gp,  337.72700_gp,  322.21420_gp,  314.90220_gp,  310.08700_gp,&
         292.73330_gp,  471.99570_gp,  417.72680_gp,  376.46650_gp,  345.08980_gp,&
         305.27830_gp,  499.17560_gp,  953.88690_gp,  738.24840_gp,  556.34680_gp,&
         561.77040_gp,  524.77420_gp,  585.12270_gp,  459.47870_gp,  430.11170_gp,&
         399.79630_gp,  386.52800_gp,  383.57630_gp,  599.85960_gp,  546.09180_gp,&
         513.67970_gp,  487.40710_gp,  446.76180_gp,  593.03760_gp, 1233.99670_gp,&
         933.45900_gp,  599.68400_gp,  587.32120_gp,  568.88850_gp,  571.17360_gp,&
         549.02970_gp,  574.57480_gp,  540.11010_gp,  548.69770_gp,  513.04540_gp,&
         498.57190_gp,  495.24940_gp,  519.84050_gp,  478.88930_gp,  642.07210_gp,&
         598.09800_gp,  547.79890_gp,  552.52070_gp,  487.85140_gp,  459.60480_gp/)
    vdwparams%coeffs(26501:26600)=(/  439.37670_gp,  420.09830_gp,  407.12160_gp,  645.11470_gp,  591.55410_gp,&
         575.53810_gp,  561.08050_gp,  525.66240_gp,  668.13180_gp, 1232.24480_gp,&
         709.44400_gp,  814.73780_gp,  729.95170_gp,  647.69500_gp,  623.85580_gp,&
         742.28920_gp,  174.55350_gp,  172.02380_gp,  127.05660_gp,  102.26240_gp,&
         70.59260_gp,  307.70750_gp,  371.55110_gp,  358.16390_gp,  312.78040_gp,&
         264.20220_gp,  431.56730_gp,  416.32140_gp,  421.58300_gp,  385.84400_gp,&
         292.55420_gp,  251.31440_gp,  246.58940_gp,  287.27800_gp,  268.12740_gp,&
         376.53350_gp,  398.18770_gp,  366.84470_gp,  342.42490_gp,  517.36980_gp,&
         500.39000_gp,  505.40540_gp,  488.64930_gp,  433.86770_gp,  384.25630_gp,&
         363.58980_gp,  369.07770_gp,  385.10740_gp,  486.43230_gp,  528.60140_gp,&
         500.08910_gp,  483.64170_gp,  626.62040_gp,  651.03510_gp,  484.84790_gp,&
         505.60090_gp,  471.47570_gp,  429.54230_gp,  405.08190_gp,  411.78430_gp,&
         432.90980_gp,  521.15810_gp,  550.85190_gp,  562.29090_gp,  558.97120_gp,&
         49.43590_gp,   32.22050_gp,  822.15000_gp,  458.99560_gp,  304.09170_gp,&
         203.12550_gp,  140.87730_gp,  106.03060_gp,   79.89570_gp,   61.16760_gp,&
         980.96360_gp,  735.28720_gp,  671.85720_gp,  523.30930_gp,  404.36880_gp,&
         333.71710_gp,  271.14080_gp,  220.78740_gp, 1618.32100_gp, 1306.02360_gp,&
         1073.86450_gp, 1032.91430_gp,  943.57380_gp,  743.19420_gp,  809.55010_gp,&
         634.99500_gp,  669.51560_gp,  691.28830_gp,  530.01490_gp,  540.01980_gp,&
         641.42920_gp,  560.53310_gp,  474.41890_gp,  423.76020_gp,  369.11450_gp/)
    vdwparams%coeffs(26601:26700)=(/  319.63860_gp, 1810.84410_gp, 1557.69030_gp, 1353.44870_gp, 1210.43840_gp,&
         1100.28590_gp,  844.97140_gp,  944.57550_gp,  715.40020_gp,  781.73110_gp,&
         723.44940_gp,  604.07200_gp,  635.77580_gp,  802.95170_gp,  736.74350_gp,&
         651.89600_gp,  602.35480_gp,  542.27570_gp,  485.56750_gp, 2206.52940_gp,&
         1992.01980_gp, 1735.70450_gp,  771.35460_gp, 1759.79550_gp, 1687.32820_gp,&
         1644.31720_gp, 1604.83510_gp, 1569.80440_gp, 1223.75470_gp, 1397.02370_gp,&
         1345.35230_gp, 1412.79420_gp, 1382.45900_gp, 1354.96960_gp, 1339.58520_gp,&
         1122.03890_gp, 1094.68760_gp,  994.28480_gp,  836.75230_gp,  849.58700_gp,&
         766.88670_gp,  700.06410_gp,  579.64280_gp,  540.79780_gp,  555.40190_gp,&
         823.07860_gp,  799.69300_gp,  730.00700_gp,  693.72340_gp,  637.52300_gp,&
         582.12360_gp, 2069.34750_gp, 1961.58220_gp, 1720.34150_gp, 1532.97370_gp,&
         1528.49670_gp, 1479.61380_gp, 1531.11370_gp, 1481.22610_gp,   80.30480_gp,&
         266.20380_gp,  338.70130_gp,  255.01410_gp,  190.65070_gp,  131.39800_gp,&
         94.99090_gp,   64.68940_gp,  390.28390_gp,  604.79080_gp,  607.64500_gp,&
         482.92640_gp,  392.40610_gp,  330.09010_gp,  268.35530_gp,  535.22110_gp,&
         1045.26550_gp,  530.05110_gp,  511.59590_gp,  501.42290_gp,  497.34050_gp,&
         455.00680_gp,  420.43750_gp,  401.08240_gp,  392.26860_gp,  389.15420_gp,&
         362.39480_gp,  601.02950_gp,  523.34040_gp,  466.17350_gp,  423.95730_gp,&
         371.77260_gp,  637.44040_gp, 1273.84400_gp,  958.36130_gp,  701.95430_gp,&
         708.92730_gp,  659.79250_gp,  745.63430_gp,  571.94600_gp,  534.79690_gp/)
    vdwparams%coeffs(26701:26800)=(/  496.08770_gp,  480.18350_gp,  473.04820_gp,  763.29290_gp,  685.12580_gp,&
         638.30590_gp,  601.76620_gp,  547.40910_gp,  749.57040_gp, 1658.17710_gp,&
         1215.02910_gp,  747.77800_gp,  732.24530_gp,  708.78020_gp,  713.08040_gp,&
         689.79990_gp,  717.79230_gp,  673.51320_gp,  686.77710_gp,  638.67240_gp,&
         620.24820_gp,  616.43990_gp,  648.92260_gp,  595.50270_gp,  818.02970_gp,&
         755.71890_gp,  686.26340_gp,  697.23760_gp,  603.79550_gp,  567.13440_gp,&
         541.15970_gp,  517.41800_gp,  498.05080_gp,  818.61890_gp,  740.24200_gp,&
         714.98450_gp,  693.50480_gp,  645.21980_gp,  843.79470_gp, 1640.20940_gp,&
         885.01740_gp, 1028.95400_gp,  918.19730_gp,  806.43310_gp,  775.04710_gp,&
         946.27710_gp,  217.52840_gp,  213.46830_gp,  154.65610_gp,  123.19970_gp,&
         83.59120_gp,  388.34390_gp,  470.10710_gp,  448.45450_gp,  387.32870_gp,&
         323.37880_gp,  544.93370_gp,  521.29830_gp,  527.52920_gp,  483.02220_gp,&
         361.26840_gp,  307.72170_gp,  302.18010_gp,  356.73210_gp,  331.60290_gp,&
         472.43170_gp,  497.53330_gp,  453.57900_gp,  420.57420_gp,  652.72580_gp,&
         624.38560_gp,  628.66890_gp,  608.42970_gp,  535.32960_gp,  470.69610_gp,&
         444.13890_gp,  454.38240_gp,  475.18820_gp,  608.62900_gp,  661.34830_gp,&
         620.01070_gp,  596.61110_gp,  788.41700_gp,  814.37400_gp,  600.44820_gp,&
         628.83490_gp,  581.36160_gp,  526.20370_gp,  493.62900_gp,  505.21810_gp,&
         531.38910_gp,  648.08000_gp,  684.68130_gp,  696.50710_gp,  690.09760_gp,&
         867.55840_gp,    9.74720_gp,    6.64340_gp,  129.98750_gp,   80.42760_gp/)
    vdwparams%coeffs(26801:26900)=(/   56.41160_gp,   39.24490_gp,   28.01290_gp,   21.49890_gp,   16.46070_gp,&
         12.75770_gp,  156.11380_gp,  126.54250_gp,  119.15450_gp,   96.53130_gp,&
         77.04290_gp,   64.89180_gp,   53.73430_gp,   44.46130_gp,  254.22370_gp,&
         217.94570_gp,  181.68340_gp,  177.06600_gp,  162.99900_gp,  128.72330_gp,&
         141.40710_gp,  111.35690_gp,  119.10910_gp,  122.01440_gp,   93.76600_gp,&
         97.68450_gp,  115.17060_gp,  103.69470_gp,   90.13730_gp,   81.83810_gp,&
         72.46130_gp,   63.67680_gp,  286.24860_gp,  259.52790_gp,  231.11730_gp,&
         209.90660_gp,  192.77650_gp,  150.74470_gp,  167.37830_gp,  129.27260_gp,&
         141.11090_gp,  131.40190_gp,  109.47410_gp,  116.39220_gp,  144.12590_gp,&
         135.44390_gp,  122.68730_gp,  114.97520_gp,  105.06600_gp,   95.39730_gp,&
         349.48600_gp,  329.47840_gp,  294.16920_gp,  144.62380_gp,  293.44860_gp,&
         282.46410_gp,  275.56100_gp,  269.18990_gp,  263.54840_gp,  211.25650_gp,&
         232.10280_gp,  224.59400_gp,  238.59840_gp,  233.61550_gp,  229.17770_gp,&
         226.30560_gp,  193.18890_gp,  193.17020_gp,  178.12490_gp,  151.79460_gp,&
         154.95640_gp,  141.66610_gp,  130.65360_gp,  109.40400_gp,  102.52010_gp,&
         105.73840_gp,  148.82270_gp,  147.21780_gp,  137.20280_gp,  131.95960_gp,&
         123.01010_gp,  113.77160_gp,  334.35640_gp,  328.67970_gp,  294.72000_gp,&
         269.50210_gp,  265.45030_gp,  257.12540_gp,  261.82830_gp,  254.02630_gp,&
         15.35310_gp,   47.03380_gp,   60.52340_gp,   47.85800_gp,   36.81660_gp,&
         26.14810_gp,   19.35710_gp,   13.55810_gp,   68.27550_gp,  105.30540_gp/)
    vdwparams%coeffs(26901:27000)=(/  108.55170_gp,   89.59480_gp,   74.84970_gp,   64.15730_gp,   53.16380_gp,&
         95.84970_gp,  175.14340_gp,   97.75500_gp,   94.47580_gp,   92.54510_gp,&
         91.45280_gp,   85.13120_gp,   79.08320_gp,   75.45250_gp,   73.68650_gp,&
         72.12160_gp,   68.80770_gp,  108.75080_gp,   97.49530_gp,   88.71080_gp,&
         81.85640_gp,   72.92210_gp,  114.97210_gp,  212.87590_gp,  168.13540_gp,&
         129.24090_gp,  130.43790_gp,  122.19410_gp,  134.97050_gp,  107.71280_gp,&
         100.89580_gp,   93.89970_gp,   90.65080_gp,   90.42770_gp,  138.18550_gp,&
         127.14090_gp,  120.50990_gp,  114.96400_gp,  106.03500_gp,  137.80050_gp,&
         274.59140_gp,  212.35690_gp,  140.39590_gp,  137.51220_gp,  133.26430_gp,&
         133.60610_gp,  127.77250_gp,  134.32640_gp,  126.41840_gp,  128.10090_gp,&
         120.23750_gp,  116.90070_gp,  116.07540_gp,  121.56170_gp,  112.32440_gp,&
         148.12250_gp,  138.81710_gp,  127.92520_gp,  128.30620_gp,  114.89450_gp,&
         108.45800_gp,  103.79770_gp,   99.20500_gp,   96.60020_gp,  148.99930_gp,&
         137.99420_gp,  135.00470_gp,  132.15970_gp,  124.51480_gp,  155.34870_gp,&
         275.95780_gp,  166.04890_gp,  189.11370_gp,  170.01950_gp,  151.77090_gp,&
         146.38890_gp,  171.06450_gp,   40.93580_gp,   40.51920_gp,   30.35360_gp,&
         24.59190_gp,   17.14240_gp,   71.30230_gp,   86.05320_gp,   83.65860_gp,&
         73.73970_gp,   62.87830_gp,  100.06100_gp,   97.16320_gp,   98.43630_gp,&
         90.05740_gp,   68.99580_gp,   59.63130_gp,   58.45150_gp,   67.39320_gp,&
         63.07400_gp,   87.76880_gp,   93.13010_gp,   86.54860_gp,   81.24000_gp/)
    vdwparams%coeffs(27001:27100)=(/  120.13000_gp,  117.15700_gp,  118.62430_gp,  114.63280_gp,  102.48830_gp,&
         91.25180_gp,   86.49860_gp,   87.23760_gp,   90.82840_gp,  113.48260_gp,&
         123.33400_gp,  117.53380_gp,  114.15530_gp,  145.93300_gp,  152.22500_gp,&
         114.17940_gp,  118.73150_gp,  111.42620_gp,  102.00810_gp,   96.57060_gp,&
         97.62380_gp,  102.54510_gp,  122.15410_gp,  129.15240_gp,  132.19470_gp,&
         131.78890_gp,  160.86740_gp,   31.44360_gp,    8.82100_gp,    6.11960_gp,&
         115.94980_gp,   71.68180_gp,   50.53030_gp,   35.41290_gp,   25.47990_gp,&
         19.69430_gp,   15.18830_gp,   11.85130_gp,  139.37530_gp,  112.82100_gp,&
         106.35080_gp,   86.39990_gp,   69.25150_gp,   58.57760_gp,   48.74480_gp,&
         40.54410_gp,  227.66000_gp,  194.67190_gp,  162.26090_gp,  158.27050_gp,&
         145.75800_gp,  115.36970_gp,  126.54160_gp,   99.90170_gp,  106.69170_gp,&
         109.23050_gp,   84.20280_gp,   87.61240_gp,  103.07380_gp,   92.96210_gp,&
         81.05630_gp,   73.80640_gp,   65.58680_gp,   57.86650_gp,  256.51980_gp,&
         232.00320_gp,  206.65890_gp,  187.82030_gp,  172.63950_gp,  135.37160_gp,&
         150.13980_gp,  116.32470_gp,  126.77780_gp,  118.14270_gp,   98.70320_gp,&
         104.76920_gp,  129.38530_gp,  121.65980_gp,  110.40890_gp,  103.65030_gp,&
         94.94890_gp,   86.45810_gp,  313.23470_gp,  294.71680_gp,  263.13810_gp,&
         130.23410_gp,  262.66440_gp,  252.81970_gp,  246.63480_gp,  240.92330_gp,&
         235.86410_gp,  189.38240_gp,  208.16040_gp,  201.42940_gp,  213.51460_gp,&
         209.04160_gp,  205.05940_gp,  202.46690_gp,  173.01760_gp,  172.95930_gp/)
    vdwparams%coeffs(27101:27200)=(/  159.67690_gp,  136.42100_gp,  139.25850_gp,  127.50310_gp,  117.76180_gp,&
         98.92870_gp,   92.82720_gp,   95.68530_gp,  133.97180_gp,  132.50170_gp,&
         123.64240_gp,  119.05310_gp,  111.19300_gp,  103.07880_gp,  299.72370_gp,&
         294.09590_gp,  263.75700_gp,  241.53050_gp,  238.01510_gp,  230.59460_gp,&
         234.75660_gp,  227.76830_gp,   13.79920_gp,   41.97130_gp,   54.07700_gp,&
         42.96460_gp,   33.25400_gp,   23.81360_gp,   17.77500_gp,   12.59700_gp,&
         61.01740_gp,   93.99780_gp,   96.95600_gp,   80.26180_gp,   67.30880_gp,&
         57.92210_gp,   48.23280_gp,   86.14810_gp,  156.77020_gp,   87.81530_gp,&
         84.92530_gp,   83.21110_gp,   82.23310_gp,   76.58460_gp,   71.21100_gp,&
         67.97740_gp,   66.38970_gp,   64.96410_gp,   62.02280_gp,   97.38380_gp,&
         87.47810_gp,   79.79860_gp,   73.82570_gp,   65.99670_gp,  103.45590_gp,&
         190.69760_gp,  150.78940_gp,  116.19730_gp,  117.30510_gp,  110.02550_gp,&
         121.43810_gp,   97.17440_gp,   91.11710_gp,   84.89540_gp,   81.97150_gp,&
         81.75000_gp,  124.10850_gp,  114.30720_gp,  108.48830_gp,  103.64770_gp,&
         95.81390_gp,  123.98030_gp,  246.12560_gp,  190.47850_gp,  126.45720_gp,&
         123.86820_gp,  120.06680_gp,  120.34080_gp,  115.07330_gp,  120.93240_gp,&
         113.87540_gp,  115.34910_gp,  108.32620_gp,  105.33020_gp,  104.56940_gp,&
         109.38790_gp,  101.19050_gp,  133.03630_gp,  124.86130_gp,  115.26060_gp,&
         115.51780_gp,  103.77670_gp,   98.09870_gp,   93.97930_gp,   89.90520_gp,&
         87.59300_gp,  134.10840_gp,  124.32080_gp,  121.70140_gp,  119.24010_gp/)
    vdwparams%coeffs(27201:27300)=(/  112.53730_gp,  139.75100_gp,  247.37300_gp,  149.49980_gp,  170.17100_gp,&
         153.24030_gp,  136.92210_gp,  132.12700_gp,  154.07130_gp,   36.69770_gp,&
         36.46040_gp,   27.50630_gp,   22.42790_gp,   15.80090_gp,   63.79440_gp,&
         77.00870_gp,   74.99670_gp,   66.31940_gp,   56.78260_gp,   89.82850_gp,&
         87.30480_gp,   88.47030_gp,   81.02290_gp,   62.39630_gp,   54.06760_gp,&
         52.99280_gp,   60.85130_gp,   57.00460_gp,   78.78120_gp,   83.61310_gp,&
         77.88740_gp,   73.28140_gp,  107.96880_gp,  105.38580_gp,  106.76290_gp,&
         103.27750_gp,   92.59830_gp,   82.64960_gp,   78.43680_gp,   79.02790_gp,&
         82.19300_gp,  102.18770_gp,  110.93720_gp,  105.85840_gp,  102.93730_gp,&
         131.26760_gp,  136.86860_gp,  102.87240_gp,  107.02160_gp,  100.67270_gp,&
         92.37480_gp,   87.60640_gp,   88.52270_gp,   92.90100_gp,  110.26110_gp,&
         116.48110_gp,  119.21790_gp,  118.92080_gp,  144.70670_gp,   28.41210_gp,&
         25.78090_gp,    6.77460_gp,    4.88360_gp,   81.97700_gp,   52.02370_gp,&
         37.56850_gp,   26.93510_gp,   19.77070_gp,   15.52490_gp,   12.15110_gp,&
         9.60500_gp,   98.88060_gp,   81.53520_gp,   77.65200_gp,   64.05760_gp,&
         52.15090_gp,   44.65290_gp,   37.62950_gp,   31.68040_gp,  162.01600_gp,&
         140.02730_gp,  117.11120_gp,  114.82560_gp,  106.04890_gp,   84.36310_gp,&
         92.46360_gp,   73.41850_gp,   78.46960_gp,   80.08530_gp,   62.12990_gp,&
         64.86330_gp,   75.88000_gp,   69.17100_gp,   61.04250_gp,   56.08120_gp,&
         50.33440_gp,   44.85460_gp,  183.14860_gp,  167.08650_gp,  149.89060_gp/)
    vdwparams%coeffs(27301:27400)=(/  136.95780_gp,  126.43070_gp,  100.10800_gp,  110.61190_gp,   86.61630_gp,&
         94.10070_gp,   87.94300_gp,   73.80570_gp,   78.30160_gp,   95.77270_gp,&
         90.68760_gp,   83.05770_gp,   78.48190_gp,   72.45310_gp,   66.50810_gp,&
         223.90760_gp,  212.12970_gp,  190.63700_gp,   97.98820_gp,  189.70840_gp,&
         182.75600_gp,  178.32410_gp,  174.22140_gp,  170.58700_gp,  138.42980_gp,&
         150.82040_gp,  146.14450_gp,  154.62870_gp,  151.39240_gp,  148.52730_gp,&
         146.57090_gp,  126.14450_gp,  126.85210_gp,  117.82080_gp,  101.45610_gp,&
         103.69670_gp,   95.49190_gp,   88.64130_gp,   75.08830_gp,   70.69100_gp,&
         72.86770_gp,   99.85450_gp,   99.15570_gp,   93.20290_gp,   90.18530_gp,&
         84.79910_gp,   79.15550_gp,  215.46340_gp,  212.55240_gp,  191.79540_gp,&
         177.28870_gp,  174.32070_gp,  168.96800_gp,  171.20630_gp,  166.24200_gp,&
         10.39870_gp,   30.60220_gp,   39.64780_gp,   32.16440_gp,   25.33180_gp,&
         18.51800_gp,   14.07930_gp,   10.22140_gp,   44.49240_gp,   68.31410_gp,&
         71.02310_gp,   59.68870_gp,   50.73980_gp,   44.15810_gp,   37.23810_gp,&
         63.86860_gp,  113.41610_gp,   65.52030_gp,   63.46450_gp,   62.20200_gp,&
         61.42020_gp,   57.49150_gp,   53.61710_gp,   51.23180_gp,   50.02180_gp,&
         48.76300_gp,   46.89110_gp,   71.90790_gp,   65.30410_gp,   60.15190_gp,&
         56.09620_gp,   50.62980_gp,   77.00400_gp,  138.09270_gp,  110.82380_gp,&
         86.80350_gp,   87.65040_gp,   82.52510_gp,   90.46600_gp,   73.42590_gp,&
         69.00050_gp,   64.46510_gp,   62.22090_gp,   62.21310_gp,   92.14030_gp/)
    vdwparams%coeffs(27401:27500)=(/   85.54790_gp,   81.72790_gp,   78.50670_gp,   73.09540_gp,   92.70500_gp,&
         178.05360_gp,  139.91280_gp,   95.21780_gp,   93.28410_gp,   90.48070_gp,&
         90.56820_gp,   86.33860_gp,   90.91000_gp,   85.74750_gp,   86.67400_gp,&
         81.65200_gp,   79.42910_gp,   78.81500_gp,   82.18060_gp,   76.30040_gp,&
         98.80250_gp,   93.30080_gp,   86.69090_gp,   86.49830_gp,   78.77350_gp,&
         74.72690_gp,   71.76240_gp,   68.75380_gp,   67.22930_gp,  100.08700_gp,&
         93.48210_gp,   91.89650_gp,   90.37600_gp,   85.80750_gp,  104.51280_gp,&
         179.71610_gp,  112.46690_gp,  127.28190_gp,  115.17510_gp,  103.44890_gp,&
         99.98780_gp,  114.95380_gp,   27.42250_gp,   27.48420_gp,   21.14880_gp,&
         17.49280_gp,   12.60770_gp,   47.20590_gp,   56.97130_gp,   55.91970_gp,&
         49.98180_gp,   43.31300_gp,   66.87060_gp,   65.34150_gp,   66.26070_gp,&
         60.78340_gp,   47.50330_gp,   41.48620_gp,   40.63760_gp,   46.08410_gp,&
         43.30960_gp,   58.81260_gp,   62.56090_gp,   58.79450_gp,   55.70370_gp,&
         80.58680_gp,   79.15500_gp,   80.37580_gp,   77.86810_gp,   70.42440_gp,&
         63.30810_gp,   60.26410_gp,   60.41490_gp,   62.65320_gp,   76.77980_gp,&
         83.19680_gp,   79.89120_gp,   78.02240_gp,   98.26700_gp,  102.62510_gp,&
         77.73790_gp,   80.80520_gp,   76.58520_gp,   70.73530_gp,   67.42720_gp,&
         67.88850_gp,   71.11200_gp,   83.40950_gp,   88.00520_gp,   90.19900_gp,&
         90.20050_gp,  108.38490_gp,   21.72050_gp,   19.86690_gp,   15.58170_gp,&
         21.05840_gp,   13.69820_gp,  322.79250_gp,  189.66050_gp,  127.96450_gp/)
    vdwparams%coeffs(27501:27600)=(/   86.16050_gp,   59.90740_gp,   45.08760_gp,   33.93340_gp,   25.93350_gp,&
         386.09800_gp,  301.33070_gp,  278.56390_gp,  219.94790_gp,  171.35570_gp,&
         141.85170_gp,  115.45520_gp,   94.04280_gp,  628.95370_gp,  525.83670_gp,&
         435.42050_gp,  420.87110_gp,  385.66110_gp,  303.12660_gp,  332.31560_gp,&
         260.22290_gp,  276.89480_gp,  285.11740_gp,  217.89580_gp,  224.71970_gp,&
         266.45340_gp,  235.43820_gp,  200.75140_gp,  179.86350_gp,  157.01910_gp,&
         136.12380_gp,  704.96090_gp,  625.69960_gp,  549.92420_gp,  494.95070_gp,&
         451.53080_gp,  348.33970_gp,  388.82260_gp,  295.92560_gp,  324.04080_gp,&
         300.49770_gp,  249.67590_gp,  264.72910_gp,  332.34100_gp,  308.08400_gp,&
         274.73330_gp,  254.78050_gp,  230.07980_gp,  206.45810_gp,  858.67990_gp,&
         796.10660_gp,  701.93980_gp,  324.15330_gp,  705.78730_gp,  678.19310_gp,&
         661.31780_gp,  645.79390_gp,  632.03590_gp,  498.02070_gp,  557.18410_gp,&
         537.72660_gp,  570.65080_gp,  558.63790_gp,  547.82520_gp,  541.39780_gp,&
         456.83250_gp,  450.89790_gp,  411.77370_gp,  347.23090_gp,  353.52640_gp,&
         320.38750_gp,  293.31160_gp,  243.07640_gp,  226.88500_gp,  233.69270_gp,&
         340.43160_gp,  333.73140_gp,  306.98840_gp,  292.79590_gp,  270.03130_gp,&
         247.16150_gp,  813.11660_gp,  788.41980_gp,  699.12600_gp,  629.55760_gp,&
         623.64620_gp,  603.79750_gp,  620.51110_gp,  601.13500_gp,   34.07000_gp,&
         110.21360_gp,  140.56320_gp,  107.51580_gp,   80.71110_gp,   55.77930_gp,&
         40.33740_gp,   27.44360_gp,  160.65940_gp,  248.56920_gp,  252.49530_gp/)
    vdwparams%coeffs(27601:27700)=(/  203.22260_gp,  166.26980_gp,  140.26180_gp,  114.24050_gp,  220.72640_gp,&
         419.71290_gp,  221.83390_gp,  214.01930_gp,  209.67420_gp,  207.62990_gp,&
         191.38970_gp,  177.05420_gp,  168.80270_gp,  164.99240_gp,  162.80210_gp,&
         153.07760_gp,  250.34940_gp,  220.23170_gp,  197.30880_gp,  179.91980_gp,&
         158.11820_gp,  263.24390_gp,  509.99270_gp,  392.06600_gp,  292.97650_gp,&
         295.77860_gp,  275.64350_gp,  308.43590_gp,  240.26870_gp,  224.56570_gp,&
         208.34700_gp,  201.46610_gp,  199.66700_gp,  316.98520_gp,  287.33370_gp,&
         269.24000_gp,  254.59630_gp,  232.26170_gp,  312.05400_gp,  659.74760_gp,&
         495.80770_gp,  314.34490_gp,  307.83100_gp,  298.04520_gp,  299.47280_gp,&
         288.33130_gp,  301.48740_gp,  283.10610_gp,  287.95020_gp,  268.76200_gp,&
         261.11000_gp,  259.45460_gp,  272.90680_gp,  250.82830_gp,  339.19420_gp,&
         314.85670_gp,  287.25880_gp,  290.46490_gp,  254.37210_gp,  239.06730_gp,&
         228.15940_gp,  217.89610_gp,  210.69250_gp,  339.74310_gp,  310.27180_gp,&
         301.14890_gp,  292.90320_gp,  273.34560_gp,  351.55860_gp,  657.70110_gp,&
         372.10800_gp,  428.61200_gp,  382.81860_gp,  338.70230_gp,  325.90880_gp,&
         390.82830_gp,   92.01260_gp,   90.14400_gp,   65.69680_gp,   52.31360_gp,&
         35.46020_gp,  163.11000_gp,  196.94110_gp,  189.00340_gp,  163.95650_gp,&
         137.39070_gp,  227.81080_gp,  219.10140_gp,  221.77080_gp,  202.71680_gp,&
         152.23080_gp,  130.08090_gp,  127.68330_gp,  149.99510_gp,  139.71080_gp,&
         198.48130_gp,  209.64300_gp,  192.09160_gp,  178.48350_gp,  272.62520_gp/)
    vdwparams%coeffs(27701:27800)=(/  262.76340_gp,  265.02430_gp,  255.93200_gp,  225.96360_gp,  199.15880_gp,&
         188.04460_gp,  191.49380_gp,  200.21140_gp,  255.28030_gp,  277.79270_gp,&
         261.82240_gp,  252.51890_gp,  329.59150_gp,  342.22420_gp,  253.64330_gp,&
         264.56000_gp,  245.51490_gp,  222.68310_gp,  209.25910_gp,  213.19310_gp,&
         224.44800_gp,  272.30380_gp,  288.09540_gp,  293.86510_gp,  291.67140_gp,&
         362.87330_gp,   68.43300_gp,   61.38290_gp,   45.93990_gp,  153.59450_gp,&
         24.98330_gp,   16.11040_gp,  408.08020_gp,  231.82510_gp,  154.03880_gp,&
         102.72020_gp,   70.97610_gp,   53.20610_gp,   39.91720_gp,   30.43380_gp,&
         487.17610_gp,  370.37780_gp,  339.45020_gp,  265.11890_gp,  204.87700_gp,&
         168.82250_gp,  136.84600_gp,  111.10440_gp,  798.83900_gp,  653.41740_gp,&
         538.57540_gp,  518.63330_gp,  474.15850_gp,  372.78380_gp,  407.24250_gp,&
         318.85890_gp,  337.46240_gp,  348.24770_gp,  266.31200_gp,  272.56060_gp,&
         323.81790_gp,  283.69090_gp,  240.23030_gp,  214.41210_gp,  186.49270_gp,&
         161.17030_gp,  894.00720_gp,  778.35390_gp,  678.75770_gp,  608.08480_gp,&
         553.14470_gp,  424.80720_gp,  474.94180_gp,  359.69290_gp,  393.66130_gp,&
         364.42410_gp,  303.37520_gp,  320.34290_gp,  404.32940_gp,  372.07850_gp,&
         329.66250_gp,  304.62850_gp,  274.10450_gp,  245.17980_gp, 1088.85660_gp,&
         993.27630_gp,  868.86070_gp,  389.56690_gp,  878.32220_gp,  842.82680_gp,&
         821.53240_gp,  801.97510_gp,  784.63100_gp,  613.26270_gp,  695.16170_gp,&
         669.87120_gp,  706.96300_gp,  691.90950_gp,  678.29270_gp,  670.54790_gp/)
    vdwparams%coeffs(27801:27900)=(/  562.67060_gp,  551.04570_gp,  501.07280_gp,  421.36670_gp,  428.22450_gp,&
         386.73130_gp,  353.08580_gp,  291.92880_gp,  272.21920_gp,  279.92390_gp,&
         413.74950_gp,  403.19540_gp,  368.67120_gp,  350.49670_gp,  322.07940_gp,&
         293.90180_gp, 1024.42190_gp,  979.76180_gp,  862.41890_gp,  770.60950_gp,&
         766.55200_gp,  742.02530_gp,  766.30100_gp,  741.67020_gp,   40.69720_gp,&
         134.44900_gp,  171.03980_gp,  129.08980_gp,   96.28720_gp,   66.10520_gp,&
         47.56640_gp,   32.17300_gp,  196.65260_gp,  304.67970_gp,  307.10150_gp,&
         244.62880_gp,  198.76160_gp,  166.95970_gp,  135.42340_gp,  268.97780_gp,&
         521.75350_gp,  267.71200_gp,  258.25540_gp,  253.06610_gp,  250.88270_gp,&
         230.01820_gp,  212.51260_gp,  202.63880_gp,  198.15060_gp,  196.29960_gp,&
         183.25560_gp,  303.59070_gp,  264.89450_gp,  236.02870_gp,  214.49810_gp,&
         187.83780_gp,  320.26650_gp,  634.91250_gp,  480.74960_gp,  354.00180_gp,&
         357.46120_gp,  332.59800_gp,  374.82840_gp,  288.54060_gp,  269.62410_gp,&
         249.98080_gp,  241.89890_gp,  238.76910_gp,  384.67520_gp,  346.14310_gp,&
         322.80720_gp,  304.33880_gp,  276.71610_gp,  377.49200_gp,  824.35520_gp,&
         608.87560_gp,  377.64800_gp,  369.79850_gp,  357.93640_gp,  360.01510_gp,&
         347.80870_gp,  362.50210_gp,  340.13210_gp,  346.61200_gp,  322.61560_gp,&
         313.32960_gp,  311.41380_gp,  327.94410_gp,  300.90150_gp,  411.88700_gp,&
         380.80640_gp,  346.01610_gp,  351.18480_gp,  304.65320_gp,  286.00240_gp,&
         272.77910_gp,  260.60470_gp,  251.11540_gp,  411.99350_gp,  373.52430_gp/)
    vdwparams%coeffs(27901:28000)=(/  361.18830_gp,  350.44170_gp,  326.00770_gp,  425.07230_gp,  817.53510_gp,&
         447.07850_gp,  518.36450_gp,  462.24170_gp,  406.82820_gp,  391.05750_gp,&
         475.25970_gp,  110.31340_gp,  107.96470_gp,   78.05740_gp,   61.95300_gp,&
         41.76870_gp,  196.78410_gp,  237.94340_gp,  227.18820_gp,  196.11960_gp,&
         163.53520_gp,  275.22420_gp,  263.58530_gp,  266.72430_gp,  243.96540_gp,&
         182.18540_gp,  155.10950_gp,  152.31640_gp,  179.96900_gp,  167.31730_gp,&
         238.96240_gp,  251.84820_gp,  229.63620_gp,  212.76790_gp,  329.33860_gp,&
         315.62100_gp,  317.85670_gp,  307.23810_gp,  270.18990_gp,  237.41460_gp,&
         223.92390_gp,  228.93290_gp,  239.56050_gp,  307.25150_gp,  334.20220_gp,&
         313.57940_gp,  301.73570_gp,  397.68030_gp,  411.54360_gp,  303.59840_gp,&
         317.46390_gp,  293.46710_gp,  265.45010_gp,  248.89430_gp,  254.48410_gp,&
         267.90110_gp,  326.92180_gp,  345.69090_gp,  351.94960_gp,  348.76390_gp,&
         437.64130_gp,   81.34780_gp,   72.95300_gp,   54.38570_gp,  184.02910_gp,&
         221.50460_gp,   21.25070_gp,   14.15770_gp,  304.77070_gp,  183.30560_gp,&
         126.01940_gp,   86.24270_gp,   60.76630_gp,   46.19270_gp,   35.07330_gp,&
         26.99900_gp,  365.21100_gp,  289.90900_gp,  270.35350_gp,  216.13030_gp,&
         170.40450_gp,  142.29820_gp,  116.83550_gp,   95.93360_gp,  595.10330_gp,&
         502.96540_gp,  417.74000_gp,  405.34620_gp,  372.23260_gp,  293.26800_gp,&
         321.76790_gp,  252.67400_gp,  269.48070_gp,  276.79850_gp,  212.13400_gp,&
         219.80160_gp,  259.91840_gp,  231.75640_gp,  199.50090_gp,  179.93740_gp/)
    vdwparams%coeffs(28001:28100)=(/  158.21040_gp,  138.10760_gp,  668.44010_gp,  598.74980_gp,  529.42130_gp,&
         478.51960_gp,  437.92800_gp,  340.05920_gp,  378.60580_gp,  290.20770_gp,&
         317.26290_gp,  294.79860_gp,  245.29960_gp,  260.38920_gp,  324.73460_gp,&
         303.00330_gp,  272.28130_gp,  253.82690_gp,  230.58600_gp,  208.15070_gp,&
         815.15180_gp,  761.14330_gp,  674.94890_gp,  321.13780_gp,  676.20310_gp,&
         650.26390_gp,  634.20890_gp,  619.41760_gp,  606.31410_gp,  481.60310_gp,&
         534.41950_gp,  516.38640_gp,  548.08770_gp,  536.58580_gp,  526.28420_gp,&
         519.91030_gp,  441.09130_gp,  438.00800_gp,  401.85550_gp,  340.62380_gp,&
         347.22720_gp,  316.01980_gp,  290.35810_gp,  241.87400_gp,  226.20530_gp,&
         233.13640_gp,  333.97090_gp,  328.79450_gp,  304.37900_gp,  291.51300_gp,&
         270.29090_gp,  248.70320_gp,  775.46700_gp,  756.32330_gp,  674.06280_gp,&
         611.38240_gp,  604.05420_gp,  584.96010_gp,  598.56780_gp,  580.26720_gp,&
         33.91930_gp,  106.84660_gp,  136.86740_gp,  106.39230_gp,   80.84960_gp,&
         56.65080_gp,   41.46000_gp,   28.63620_gp,  155.44950_gp,  240.15820_gp,&
         245.65290_gp,  200.14220_gp,  165.45180_gp,  140.69580_gp,  115.60060_gp,&
         215.84670_gp,  402.82860_gp,  218.43190_gp,  210.92280_gp,  206.62740_gp,&
         204.40810_gp,  189.31020_gp,  175.49190_gp,  167.37650_gp,  163.53010_gp,&
         160.72250_gp,  152.20000_gp,  244.79640_gp,  217.33500_gp,  196.21220_gp,&
         179.98620_gp,  159.26660_gp,  258.15200_gp,  489.59430_gp,  381.14140_gp,&
         288.66370_gp,  291.38130_gp,  272.24210_gp,  302.74930_gp,  238.60140_gp/)
    vdwparams%coeffs(28101:28200)=(/  223.25490_gp,  207.45050_gp,  200.43480_gp,  199.27140_gp,  310.50480_gp,&
         283.48980_gp,  267.14010_gp,  253.72280_gp,  232.74260_gp,  307.65740_gp,&
         632.61590_gp,  481.73950_gp,  311.58240_gp,  305.15280_gp,  295.58710_gp,&
         296.68160_gp,  284.70850_gp,  298.47730_gp,  280.58560_gp,  284.87000_gp,&
         266.61130_gp,  259.11400_gp,  257.37860_gp,  270.13590_gp,  248.93750_gp,&
         332.58120_gp,  310.15250_gp,  284.35150_gp,  286.39060_gp,  253.54990_gp,&
         238.81810_gp,  228.23770_gp,  218.05720_gp,  211.57770_gp,  333.79720_gp,&
         306.90290_gp,  299.03800_gp,  291.78570_gp,  273.60270_gp,  346.71520_gp,&
         632.97580_gp,  368.67540_gp,  422.37050_gp,  378.47190_gp,  336.27890_gp,&
         323.95480_gp,  383.67070_gp,   91.02060_gp,   89.63910_gp,   66.23430_gp,&
         53.20760_gp,   36.58750_gp,  159.95620_gp,  193.11320_gp,  186.50310_gp,&
         163.07640_gp,  137.85400_gp,  223.95660_gp,  216.40380_gp,  219.13910_gp,&
         200.40530_gp,  152.00100_gp,  130.62480_gp,  128.12710_gp,  149.11580_gp,&
         139.22010_gp,  195.76670_gp,  207.23900_gp,  191.22650_gp,  178.59070_gp,&
         268.45260_gp,  260.22540_gp,  262.96750_gp,  254.04600_gp,  225.70300_gp,&
         199.94180_gp,  189.15600_gp,  191.69620_gp,  199.99790_gp,  252.44560_gp,&
         274.52930_gp,  260.15780_gp,  251.79490_gp,  325.33280_gp,  338.53620_gp,&
         252.38610_gp,  262.87220_gp,  245.30820_gp,  223.53440_gp,  210.84220_gp,&
         213.97300_gp,  225.00270_gp,  270.49420_gp,  286.07870_gp,  292.30220_gp,&
         290.76120_gp,  358.40590_gp,   68.80070_gp,   61.94330_gp,   46.86350_gp/)
    vdwparams%coeffs(28201:28300)=(/  152.04520_gp,  181.48400_gp,  151.68860_gp,   18.75000_gp,   12.75230_gp,&
         270.75980_gp,  160.88700_gp,  110.53650_gp,   76.00720_gp,   53.96950_gp,&
         41.36230_gp,   31.69850_gp,   24.63430_gp,  324.62410_gp,  255.15610_gp,&
         237.39610_gp,  189.54880_gp,  149.67090_gp,  125.37250_gp,  103.39310_gp,&
         85.36160_gp,  531.30440_gp,  445.06810_gp,  369.08350_gp,  357.97830_gp,&
         328.64070_gp,  259.57220_gp,  284.01670_gp,  223.63540_gp,  237.69940_gp,&
         244.20530_gp,  187.84400_gp,  193.85700_gp,  228.67850_gp,  203.62780_gp,&
         175.41210_gp,  158.49070_gp,  139.75410_gp,  122.45360_gp,  596.83060_gp,&
         530.35890_gp,  467.86120_gp,  422.52080_gp,  386.64910_gp,  300.67530_gp,&
         334.55830_gp,  256.92950_gp,  280.39530_gp,  260.63370_gp,  217.71120_gp,&
         230.40100_gp,  286.92990_gp,  267.20820_gp,  240.03960_gp,  223.89640_gp,&
         203.68550_gp,  184.26650_gp,  727.57110_gp,  675.09330_gp,  597.11540_gp,&
         283.41080_gp,  599.87830_gp,  576.62050_gp,  562.31220_gp,  549.12780_gp,&
         537.43950_gp,  426.47230_gp,  475.30650_gp,  459.03360_gp,  485.49610_gp,&
         475.24580_gp,  466.04870_gp,  460.40630_gp,  390.29370_gp,  386.33500_gp,&
         354.38040_gp,  300.91280_gp,  306.57700_gp,  279.17480_gp,  256.71260_gp,&
         214.50850_gp,  200.89400_gp,  206.81470_gp,  295.86670_gp,  290.66840_gp,&
         268.89430_gp,  257.54500_gp,  239.00100_gp,  220.25450_gp,  690.85290_gp,&
         670.10980_gp,  596.10010_gp,  540.16650_gp,  534.72260_gp,  517.92390_gp,&
         530.80430_gp,  514.45080_gp,   29.76800_gp,   93.83560_gp,  120.18800_gp/)
    vdwparams%coeffs(28301:28400)=(/   93.45410_gp,   71.33180_gp,   50.38140_gp,   37.21780_gp,   26.10310_gp,&
         136.98410_gp,  211.33880_gp,  215.70660_gp,  175.59630_gp,  145.38010_gp,&
         123.99050_gp,  102.32230_gp,  190.95210_gp,  357.07210_gp,  192.52210_gp,&
         186.03690_gp,  182.33530_gp,  180.47740_gp,  166.94170_gp,  154.86080_gp,&
         147.80570_gp,  144.45180_gp,  142.16340_gp,  134.37060_gp,  215.33740_gp,&
         191.00140_gp,  172.55780_gp,  158.54630_gp,  140.68120_gp,  228.51060_gp,&
         434.45030_gp,  336.91730_gp,  254.66980_gp,  257.19610_gp,  240.51390_gp,&
         267.86660_gp,  210.97970_gp,  197.65180_gp,  183.89000_gp,  177.79340_gp,&
         176.46970_gp,  274.25370_gp,  250.02640_gp,  235.52760_gp,  223.79990_gp,&
         205.56330_gp,  271.60040_gp,  562.14810_gp,  426.01980_gp,  275.02300_gp,&
         269.36720_gp,  260.96100_gp,  261.92990_gp,  251.67010_gp,  263.38860_gp,&
         247.70800_gp,  251.53730_gp,  235.34260_gp,  228.72630_gp,  227.16850_gp,&
         238.21710_gp,  219.67420_gp,  293.56720_gp,  273.85420_gp,  251.23200_gp,&
         253.17100_gp,  224.27020_gp,  211.52860_gp,  202.38660_gp,  193.63610_gp,&
         187.79620_gp,  295.44290_gp,  271.27660_gp,  264.14600_gp,  257.74500_gp,&
         241.87370_gp,  305.97750_gp,  561.60150_gp,  325.25320_gp,  373.17930_gp,&
         334.81340_gp,  297.40820_gp,  286.58400_gp,  340.06000_gp,   79.78450_gp,&
         78.86240_gp,   58.60570_gp,   47.41680_gp,   33.02150_gp,  140.35180_gp,&
         169.49030_gp,  163.67210_gp,  143.35180_gp,  121.52560_gp,  197.32010_gp,&
         190.55180_gp,  193.00120_gp,  176.76920_gp,  134.62870_gp,  115.92230_gp/)
    vdwparams%coeffs(28401:28500)=(/  113.72870_gp,  132.00850_gp,  123.31710_gp,  172.20040_gp,  182.17310_gp,&
         168.21430_gp,  157.34260_gp,  236.76080_gp,  229.25270_gp,  231.67830_gp,&
         224.15900_gp,  199.53250_gp,  177.10190_gp,  167.74480_gp,  170.07730_gp,&
         177.30100_gp,  222.99250_gp,  242.13510_gp,  229.40100_gp,  222.11080_gp,&
         286.97840_gp,  298.14430_gp,  222.47740_gp,  232.03320_gp,  216.83460_gp,&
         197.94790_gp,  186.97210_gp,  189.93130_gp,  199.53230_gp,  239.39990_gp,&
         252.89420_gp,  258.18280_gp,  256.81660_gp,  316.07620_gp,   60.65370_gp,&
         54.90660_gp,   41.91930_gp,  133.38200_gp,  159.40630_gp,  133.39560_gp,&
         118.28630_gp,   20.18530_gp,   13.43050_gp,  320.65850_gp,  182.69550_gp,&
         122.39150_gp,   82.59300_gp,   57.84690_gp,   43.91610_gp,   33.39240_gp,&
         25.79400_gp,  383.37730_gp,  291.92630_gp,  268.12070_gp,  210.41880_gp,&
         163.72990_gp,  135.84700_gp,  111.02480_gp,   90.95530_gp,  630.14120_gp,&
         515.54150_gp,  425.12770_gp,  410.01290_gp,  375.17450_gp,  295.85990_gp,&
         322.68450_gp,  253.53650_gp,  267.94820_gp,  276.23270_gp,  212.14340_gp,&
         216.94470_gp,  256.79600_gp,  225.70800_gp,  192.11030_gp,  172.27180_gp,&
         150.74260_gp,  131.16000_gp,  705.93280_gp,  614.61390_gp,  536.67450_gp,&
         481.50000_gp,  438.65020_gp,  338.35350_gp,  377.63020_gp,  287.46140_gp,&
         313.95770_gp,  291.03090_gp,  243.20090_gp,  256.36430_gp,  322.07710_gp,&
         296.88720_gp,  263.94960_gp,  244.63100_gp,  221.02500_gp,  198.65840_gp,&
         859.72520_gp,  784.52520_gp,  686.99020_gp,  312.11020_gp,  694.69690_gp/)
    vdwparams%coeffs(28501:28600)=(/  666.73430_gp,  649.90940_gp,  634.44150_gp,  620.71840_gp,  486.71800_gp,&
         550.92420_gp,  530.97990_gp,  559.38190_gp,  547.44560_gp,  536.65810_gp,&
         530.43510_gp,  445.99990_gp,  437.00330_gp,  398.22740_gp,  336.18450_gp,&
         341.73640_gp,  309.43260_gp,  283.22570_gp,  235.39680_gp,  219.99670_gp,&
         226.09290_gp,  330.90130_gp,  322.64990_gp,  295.78510_gp,  281.78430_gp,&
         259.80700_gp,  238.01290_gp,  809.81600_gp,  774.54730_gp,  682.71560_gp,&
         611.85740_gp,  608.71940_gp,  589.42580_gp,  608.18860_gp,  588.74820_gp,&
         32.51490_gp,  106.15450_gp,  135.27020_gp,  102.92200_gp,   77.51240_gp,&
         53.96900_gp,   39.41250_gp,   27.27050_gp,  155.57860_gp,  240.53020_gp,&
         242.84190_gp,  194.41840_gp,  158.94530_gp,  134.37370_gp,  109.89000_gp,&
         214.49610_gp,  412.61420_gp,  213.61290_gp,  206.26510_gp,  202.20500_gp,&
         200.45980_gp,  184.04780_gp,  170.31170_gp,  162.52860_gp,  158.94050_gp,&
         157.33960_gp,  147.15970_gp,  241.00160_gp,  211.04400_gp,  188.84500_gp,&
         172.34930_gp,  151.80060_gp,  255.85120_gp,  502.46340_gp,  381.79030_gp,&
         282.70210_gp,  285.59450_gp,  266.25680_gp,  299.44610_gp,  231.82850_gp,&
         216.98060_gp,  201.55360_gp,  195.09170_gp,  192.60990_gp,  306.71200_gp,&
         276.63930_gp,  258.61940_gp,  244.42920_gp,  223.08960_gp,  301.64000_gp,&
         652.39960_gp,  483.48120_gp,  302.68200_gp,  296.42470_gp,  287.01630_gp,&
         288.52060_gp,  278.59570_gp,  290.30390_gp,  272.64920_gp,  277.61700_gp,&
         258.69870_gp,  251.29970_gp,  249.69380_gp,  262.46920_gp,  241.27880_gp/)
    vdwparams%coeffs(28601:28700)=(/  328.23180_gp,  304.26470_gp,  277.31590_gp,  281.03260_gp,  245.28440_gp,&
         230.80460_gp,  220.51410_gp,  210.98620_gp,  203.59070_gp,  329.52070_gp,&
         299.46420_gp,  289.96580_gp,  281.77690_gp,  262.91680_gp,  339.61520_gp,&
         647.57740_gp,  358.09400_gp,  414.47150_gp,  370.53510_gp,  326.89030_gp,&
         314.49130_gp,  380.28110_gp,   87.78520_gp,   86.38950_gp,   63.20920_gp,&
         50.71980_gp,   34.85960_gp,  156.13270_gp,  188.76990_gp,  180.76730_gp,&
         156.85490_gp,  131.67410_gp,  219.36790_gp,  210.45590_gp,  213.05010_gp,&
         195.16290_gp,  146.98850_gp,  125.72050_gp,  123.44170_gp,  144.84240_gp,&
         134.89330_gp,  190.47990_gp,  200.85740_gp,  183.86230_gp,  171.00310_gp,&
         262.89820_gp,  252.42450_gp,  254.45670_gp,  246.30770_gp,  217.64570_gp,&
         192.05910_gp,  181.51960_gp,  185.21540_gp,  193.49000_gp,  246.15700_gp,&
         267.32250_gp,  251.43010_gp,  242.41490_gp,  317.82550_gp,  328.83880_gp,&
         243.49370_gp,  254.69760_gp,  236.40400_gp,  214.67590_gp,  201.91630_gp,&
         206.22720_gp,  216.81400_gp,  262.94210_gp,  277.71100_gp,  282.77180_gp,&
         280.49880_gp,  349.77880_gp,   65.56180_gp,   59.19800_gp,   44.74420_gp,&
         146.61210_gp,  176.31430_gp,  145.41330_gp,  128.86410_gp,  141.87630_gp,&
         22.53900_gp,   14.92960_gp,  332.34640_gp,  197.74540_gp,  134.88560_gp,&
         91.76610_gp,   64.39770_gp,   48.83100_gp,   37.01290_gp,   28.46540_gp,&
         398.01800_gp,  313.46580_gp,  291.17040_gp,  231.54320_gp,  181.69810_gp,&
         151.26110_gp,  123.84420_gp,  101.45720_gp,  648.68790_gp,  545.45480_gp/)
    vdwparams%coeffs(28701:28800)=(/  452.41480_gp,  438.28120_gp,  402.12020_gp,  316.61340_gp,  347.15450_gp,&
         272.39980_gp,  290.12680_gp,  298.31410_gp,  228.48110_gp,  236.16610_gp,&
         279.43460_gp,  248.18660_gp,  212.82990_gp,  191.48490_gp,  167.95000_gp,&
         146.29100_gp,  727.99780_gp,  649.26400_gp,  572.54520_gp,  516.56560_gp,&
         472.14060_gp,  365.73900_gp,  407.59020_gp,  311.62540_gp,  340.83280_gp,&
         316.47210_gp,  263.32580_gp,  279.28860_gp,  349.13300_gp,  324.82290_gp,&
         290.96000_gp,  270.67740_gp,  245.34690_gp,  221.01200_gp,  887.20920_gp,&
         825.70960_gp,  730.32590_gp,  343.26010_gp,  733.03460_gp,  704.67870_gp,&
         687.22230_gp,  671.14820_gp,  656.90430_gp,  520.04840_gp,  579.21270_gp,&
         559.35900_gp,  593.50440_gp,  581.02860_gp,  569.82860_gp,  563.01640_gp,&
         476.57670_gp,  471.86840_gp,  432.10850_gp,  365.58740_gp,  372.47260_gp,&
         338.45240_gp,  310.56410_gp,  258.29950_gp,  241.43850_gp,  248.73340_gp,&
         358.61390_gp,  352.36140_gp,  325.32280_gp,  311.04340_gp,  287.80550_gp,&
         264.31470_gp,  842.29970_gp,  819.26690_gp,  728.58620_gp,  658.87100_gp,&
         651.82590_gp,  631.19280_gp,  647.15530_gp,  627.18910_gp,   36.14470_gp,&
         115.12840_gp,  147.19230_gp,  113.67480_gp,   86.01720_gp,   60.02290_gp,&
         43.79620_gp,   30.15570_gp,  167.72780_gp,  259.20610_gp,  264.29930_gp,&
         214.22700_gp,  176.38160_gp,  149.56740_gp,  122.54290_gp,  231.99890_gp,&
         436.28440_gp,  234.01650_gp,  225.91360_gp,  221.33820_gp,  219.07040_gp,&
         202.47400_gp,  187.55800_gp,  178.87740_gp,  174.80580_gp,  172.12150_gp/)
    vdwparams%coeffs(28801:28900)=(/  162.47350_gp,  262.91340_gp,  232.50850_gp,  209.26950_gp,  191.54240_gp,&
         169.09250_gp,  277.17010_gp,  530.24480_gp,  410.52730_gp,  309.19640_gp,&
         312.15680_gp,  291.39440_gp,  324.90950_gp,  254.88120_gp,  238.43110_gp,&
         221.46490_gp,  214.07940_gp,  212.51480_gp,  333.48150_gp,  303.51790_gp,&
         285.33630_gp,  270.52970_gp,  247.64850_gp,  329.44100_gp,  685.50020_gp,&
         518.98320_gp,  332.98690_gp,  326.10980_gp,  315.83740_gp,  317.14450_gp,&
         304.81740_gp,  319.12880_gp,  299.89300_gp,  304.69800_gp,  284.84920_gp,&
         276.80030_gp,  274.98180_gp,  288.82640_gp,  265.90550_gp,  356.95090_gp,&
         332.27960_gp,  304.07910_gp,  306.77590_gp,  270.45010_gp,  254.57310_gp,&
         243.21020_gp,  232.38620_gp,  225.15610_gp,  358.16930_gp,  328.36650_gp,&
         319.41280_gp,  311.25260_gp,  291.32130_gp,  371.19590_gp,  684.77170_gp,&
         394.03910_gp,  452.47860_gp,  404.99240_gp,  359.26250_gp,  345.95630_gp,&
         411.86560_gp,   97.23460_gp,   95.60390_gp,   70.31840_gp,   56.36120_gp,&
         38.62630_gp,  171.53850_gp,  207.08710_gp,  199.47600_gp,  173.89310_gp,&
         146.53820_gp,  240.08000_gp,  231.51840_gp,  234.40930_gp,  214.38080_gp,&
         162.05470_gp,  138.99280_gp,  136.38180_gp,  159.25840_gp,  148.56420_gp,&
         209.51300_gp,  221.56120_gp,  203.87080_gp,  190.04420_gp,  287.61670_gp,&
         278.11250_gp,  280.82410_gp,  271.31370_gp,  240.51090_gp,  212.69370_gp,&
         201.10400_gp,  204.22770_gp,  213.23450_gp,  270.09880_gp,  293.72900_gp,&
         277.71860_gp,  268.41360_gp,  348.19780_gp,  361.93370_gp,  269.24010_gp/)
    vdwparams%coeffs(28901:29000)=(/  280.63840_gp,  261.36250_gp,  237.78890_gp,  224.00480_gp,  227.73200_gp,&
         239.55770_gp,  288.97680_gp,  305.61160_gp,  311.99930_gp,  310.06550_gp,&
         383.47400_gp,   73.07540_gp,   65.75840_gp,   49.61590_gp,  162.40650_gp,&
         194.17780_gp,  161.55910_gp,  142.16340_gp,  155.48130_gp,  172.31530_gp,&
         25.92550_gp,   17.06140_gp,  392.41090_gp,  230.62870_gp,  156.27830_gp,&
         105.80040_gp,   73.96750_gp,   55.93140_gp,   42.28790_gp,   32.45150_gp,&
         469.54960_gp,  366.33190_gp,  339.10340_gp,  268.43700_gp,  209.85610_gp,&
         174.27160_gp,  142.34130_gp,  116.35890_gp,  766.57140_gp,  639.80000_gp,&
         529.77990_gp,  512.42740_gp,  469.71220_gp,  369.70490_gp,  404.96220_gp,&
         317.59430_gp,  337.69310_gp,  347.54450_gp,  266.09200_gp,  274.33390_gp,&
         324.94560_gp,  287.61030_gp,  245.87190_gp,  220.78460_gp,  193.25910_gp,&
         168.01940_gp,  859.65830_gp,  761.73640_gp,  669.74730_gp,  603.15190_gp,&
         550.60110_gp,  425.58490_gp,  474.66380_gp,  362.03660_gp,  396.01830_gp,&
         367.42870_gp,  305.77550_gp,  323.92920_gp,  405.93030_gp,  376.61080_gp,&
         336.42890_gp,  312.46310_gp,  282.71950_gp,  254.24580_gp, 1047.50770_gp,&
         969.65590_gp,  855.12980_gp,  397.08240_gp,  859.90430_gp,  826.24710_gp,&
         805.66910_gp,  786.73400_gp,  769.95120_gp,  607.50730_gp,  679.77520_gp,&
         656.08890_gp,  695.13050_gp,  680.46640_gp,  667.27580_gp,  659.39100_gp,&
         556.87690_gp,  549.81880_gp,  502.58340_gp,  424.56930_gp,  432.28120_gp,&
         392.18940_gp,  359.41840_gp,  298.49470_gp,  278.84190_gp,  287.12700_gp/)
    vdwparams%coeffs(29001:29100)=(/  416.57100_gp,  408.43950_gp,  376.18150_gp,  359.16750_gp,  331.77440_gp,&
         304.22370_gp,  992.08740_gp,  960.60620_gp,  851.94620_gp,  768.02670_gp,&
         760.90240_gp,  736.74860_gp,  756.79710_gp,  733.18910_gp,   41.73230_gp,&
         134.14080_gp,  171.28910_gp,  131.52340_gp,   99.17630_gp,   68.93230_gp,&
         50.12690_gp,   34.36420_gp,  195.62310_gp,  302.52130_gp,  307.55040_gp,&
         248.19650_gp,  203.68650_gp,  172.32720_gp,  140.84920_gp,  269.84030_gp,&
         511.53740_gp,  271.24750_gp,  261.80740_gp,  256.51490_gp,  253.99050_gp,&
         234.26800_gp,  216.86750_gp,  206.82100_gp,  202.14420_gp,  199.34300_gp,&
         187.64810_gp,  305.45790_gp,  269.22150_gp,  241.71350_gp,  220.85650_gp,&
         194.59280_gp,  322.10960_gp,  621.92630_gp,  478.71320_gp,  358.47720_gp,&
         361.92950_gp,  337.57750_gp,  377.42760_gp,  294.66340_gp,  275.57070_gp,&
         255.84410_gp,  247.37280_gp,  245.21190_gp,  387.34040_gp,  351.52010_gp,&
         329.80470_gp,  312.26190_gp,  285.38200_gp,  382.05620_gp,  804.95990_gp,&
         605.48420_gp,  385.13200_gp,  377.16470_gp,  365.22950_gp,  366.89290_gp,&
         353.08810_gp,  369.24660_gp,  346.85640_gp,  352.67810_gp,  329.34030_gp,&
         319.99040_gp,  317.92410_gp,  334.14390_gp,  307.37050_gp,  414.66360_gp,&
         385.33830_gp,  352.00450_gp,  355.66060_gp,  312.28140_gp,  293.75400_gp,&
         280.52440_gp,  268.02920_gp,  259.33470_gp,  415.73130_gp,  380.05970_gp,&
         369.14810_gp,  359.33580_gp,  335.82630_gp,  430.41160_gp,  802.57380_gp,&
         455.79310_gp,  524.67090_gp,  469.16950_gp,  415.36370_gp,  399.79830_gp/)
    vdwparams%coeffs(29101:29200)=(/  478.45560_gp,  112.47080_gp,  110.46530_gp,   80.90840_gp,   64.69460_gp,&
         44.16110_gp,  198.95420_gp,  240.29350_gp,  230.95630_gp,  200.85430_gp,&
         168.82630_gp,  278.42050_gp,  268.02690_gp,  271.33540_gp,  248.16060_gp,&
         187.02850_gp,  160.12070_gp,  157.14100_gp,  184.02720_gp,  171.52460_gp,&
         242.66150_gp,  256.40090_gp,  235.40190_gp,  219.11610_gp,  333.46480_gp,&
         321.71900_gp,  324.64310_gp,  313.69090_gp,  277.53210_gp,  245.04710_gp,&
         231.55270_gp,  235.53230_gp,  246.04900_gp,  312.61830_gp,  339.98590_gp,&
         320.84600_gp,  309.76360_gp,  403.45600_gp,  418.89570_gp,  310.96060_gp,&
         324.38890_gp,  301.55580_gp,  273.96610_gp,  257.78980_gp,  262.44250_gp,&
         276.11510_gp,  334.01410_gp,  353.22530_gp,  360.36290_gp,  357.87940_gp,&
         444.25160_gp,   84.13980_gp,   75.65140_gp,   56.91340_gp,  187.79330_gp,&
         224.91460_gp,  186.42580_gp,  163.95940_gp,  179.81780_gp,  198.98010_gp,&
         229.94980_gp,   25.81440_gp,   17.19000_gp,  373.09010_gp,  223.60990_gp,&
         153.40550_gp,  104.84620_gp,   73.82460_gp,   56.10570_gp,   42.60130_gp,&
         32.80210_gp,  447.01790_gp,  353.91780_gp,  329.66240_gp,  263.15950_gp,&
         207.24190_gp,  172.94560_gp,  141.92450_gp,  116.49640_gp,  728.61950_gp,&
         614.67050_gp,  510.29120_gp,  494.92760_gp,  454.38220_gp,  357.97570_gp,&
         392.64080_gp,  308.30750_gp,  328.64030_gp,  337.66070_gp,  258.78300_gp,&
         267.91080_gp,  316.82010_gp,  282.18280_gp,  242.67360_gp,  218.75310_gp,&
         192.24210_gp,  167.75050_gp,  818.22470_gp,  731.74620_gp,  646.47740_gp/)
    vdwparams%coeffs(29201:29300)=(/  584.01600_gp,  534.29150_gp,  414.65090_gp,  461.75770_gp,  353.73830_gp,&
         386.72900_gp,  359.28420_gp,  299.01370_gp,  317.28800_gp,  395.91570_gp,&
         369.10130_gp,  331.39460_gp,  308.77390_gp,  280.35770_gp,  252.96780_gp,&
         997.62710_gp,  930.38380_gp,  824.34880_gp,  390.90330_gp,  826.40440_gp,&
         794.61360_gp,  774.97180_gp,  756.87850_gp,  740.84790_gp,  587.91190_gp,&
         653.18150_gp,  631.03160_gp,  669.58290_gp,  655.52070_gp,  642.91690_gp,&
         635.15750_gp,  538.51890_gp,  534.25460_gp,  489.90920_gp,  415.09390_gp,&
         423.06850_gp,  384.89200_gp,  353.53050_gp,  294.42160_gp,  275.32990_gp,&
         283.71810_gp,  407.09580_gp,  400.53490_gp,  370.51780_gp,  354.69940_gp,&
         328.71350_gp,  302.33200_gp,  948.44440_gp,  924.07530_gp,  823.01470_gp,&
         745.85710_gp,  737.24570_gp,  713.93940_gp,  730.98970_gp,  708.57740_gp,&
         41.24340_gp,  130.30340_gp,  166.82950_gp,  129.46150_gp,   98.29300_gp,&
         68.82640_gp,   50.35550_gp,   34.77850_gp,  189.67700_gp,  293.04400_gp,&
         299.46280_gp,  243.64100_gp,  201.21190_gp,  171.00300_gp,  140.42820_gp,&
         263.15470_gp,  492.17040_gp,  266.02940_gp,  256.87630_gp,  251.65870_gp,&
         248.99550_gp,  230.46590_gp,  213.60980_gp,  203.73690_gp,  199.06950_gp,&
         195.75950_gp,  185.20600_gp,  298.30600_gp,  264.55660_gp,  238.66010_gp,&
         218.81490_gp,  193.53150_gp,  314.65180_gp,  598.21020_gp,  464.92540_gp,&
         351.56410_gp,  354.89700_gp,  331.52300_gp,  368.96070_gp,  290.41810_gp,&
         271.73870_gp,  252.49180_gp,  243.99260_gp,  242.46050_gp,  378.45390_gp/)
    vdwparams%coeffs(29301:29400)=(/  345.21760_gp,  325.09990_gp,  308.63590_gp,  282.98000_gp,  374.67690_gp,&
         773.12930_gp,  587.68050_gp,  379.25720_gp,  371.43070_gp,  359.77520_gp,&
         361.14910_gp,  346.73650_gp,  363.34580_gp,  341.54120_gp,  346.82760_gp,&
         324.49690_gp,  315.36050_gp,  313.25730_gp,  328.83440_gp,  302.96350_gp,&
         405.27710_gp,  377.77130_gp,  346.18800_gp,  348.83500_gp,  308.49820_gp,&
         290.54470_gp,  277.66310_gp,  265.30450_gp,  257.31840_gp,  406.78760_gp,&
         373.70420_gp,  363.95410_gp,  355.00400_gp,  332.73420_gp,  422.21560_gp,&
         773.16810_gp,  448.75130_gp,  514.46500_gp,  460.88700_gp,  409.32390_gp,&
         394.28590_gp,  467.65480_gp,  110.74050_gp,  109.03310_gp,   80.49020_gp,&
         64.64460_gp,   44.44170_gp,  194.81980_gp,  235.20410_gp,  226.99490_gp,&
         198.34020_gp,  167.54910_gp,  272.79590_gp,  263.44960_gp,  266.77170_gp,&
         243.98790_gp,  184.93070_gp,  158.85780_gp,  155.83590_gp,  181.50050_gp,&
         169.42270_gp,  238.33640_gp,  252.22610_gp,  232.57440_gp,  217.11610_gp,&
         326.96240_gp,  316.71890_gp,  319.99160_gp,  309.16430_gp,  274.53760_gp,&
         243.11650_gp,  229.98060_gp,  233.20160_gp,  243.33990_gp,  307.38320_gp,&
         334.24830_gp,  316.55690_gp,  306.27350_gp,  396.13570_gp,  412.06100_gp,&
         307.03930_gp,  319.88160_gp,  298.37090_gp,  271.79960_gp,  256.30050_gp,&
         260.23930_gp,  273.66690_gp,  329.26130_gp,  348.20710_gp,  355.68970_gp,&
         353.72780_gp,  436.36490_gp,   83.60360_gp,   75.28080_gp,   56.93980_gp,&
         184.98860_gp,  220.91650_gp,  184.43640_gp,  162.28780_gp,  177.05570_gp/)
    vdwparams%coeffs(29401:29500)=(/  196.51790_gp,  226.80190_gp,  224.28510_gp,   30.31930_gp,   20.34020_gp,&
         444.69780_gp,  264.04390_gp,  180.45680_gp,  123.27050_gp,   86.95690_gp,&
         66.27410_gp,   50.51560_gp,   39.06740_gp,  532.85380_gp,  418.85670_gp,&
         389.06980_gp,  309.67790_gp,  243.51060_gp,  203.20970_gp,  166.88670_gp,&
         137.20100_gp,  869.95980_gp,  729.95940_gp,  605.30470_gp,  586.57520_gp,&
         538.26720_gp,  424.40160_gp,  464.84110_gp,  365.31290_gp,  388.62700_gp,&
         399.51060_gp,  306.59770_gp,  316.53820_gp,  373.94790_gp,  332.29320_gp,&
         285.35990_gp,  257.14220_gp,  226.02240_gp,  197.37620_gp,  976.61940_gp,&
         869.27450_gp,  766.40590_gp,  691.59970_gp,  632.35800_gp,  490.57200_gp,&
         546.39300_gp,  418.48510_gp,  457.27800_gp,  424.77830_gp,  354.12210_gp,&
         375.14860_gp,  468.25210_gp,  435.62520_gp,  390.51070_gp,  363.59740_gp,&
         330.01060_gp,  297.77960_gp, 1190.09270_gp, 1105.93390_gp,  977.85820_gp,&
         460.88570_gp,  982.23880_gp,  944.18200_gp,  920.77280_gp,  899.21070_gp,&
         880.09790_gp,  697.17420_gp,  777.03650_gp,  750.35490_gp,  795.06940_gp,&
         778.32350_gp,  763.28690_gp,  754.12760_gp,  638.57600_gp,  631.85040_gp,&
         578.90360_gp,  490.47740_gp,  499.68040_gp,  454.40140_gp,  417.30580_gp,&
         347.79410_gp,  325.38710_gp,  335.07650_gp,  481.76440_gp,  473.18850_gp,&
         437.07980_gp,  418.10900_gp,  387.26470_gp,  356.12660_gp, 1129.67620_gp,&
         1097.28290_gp,  975.75910_gp,  882.88780_gp,  873.88960_gp,  846.33640_gp,&
         867.87930_gp,  841.09640_gp,   48.42700_gp,  153.82640_gp,  196.74540_gp/)
    vdwparams%coeffs(29501:29600)=(/  152.26150_gp,  115.61300_gp,   81.11100_gp,   59.53440_gp,   41.37460_gp,&
         224.40970_gp,  346.49400_gp,  353.26540_gp,  286.64120_gp,  236.44610_gp,&
         200.95580_gp,  165.14950_gp,  311.31870_gp,  584.48340_gp,  313.79200_gp,&
         303.05660_gp,  296.98540_gp,  293.98400_gp,  271.72730_gp,  251.85120_gp,&
         240.28600_gp,  234.83900_gp,  231.27000_gp,  218.30270_gp,  351.91540_gp,&
         311.42600_gp,  280.63500_gp,  257.22750_gp,  227.54720_gp,  372.15270_gp,&
         710.67900_gp,  550.21640_gp,  414.78140_gp,  418.84820_gp,  391.27220_gp,&
         436.21610_gp,  342.63080_gp,  320.74540_gp,  298.15620_gp,  288.28650_gp,&
         286.07710_gp,  447.32690_gp,  407.22750_gp,  383.02800_gp,  363.40770_gp,&
         333.08010_gp,  442.05370_gp,  919.09280_gp,  695.61840_gp,  447.15510_gp,&
         437.94100_gp,  424.19840_gp,  425.89770_gp,  409.41900_gp,  428.43500_gp,&
         402.74990_gp,  409.14080_gp,  382.56860_gp,  371.77680_gp,  369.29910_gp,&
         387.63420_gp,  357.09720_gp,  478.65480_gp,  445.89750_gp,  408.42780_gp,&
         411.95430_gp,  363.77170_gp,  342.72900_gp,  327.66460_gp,  313.31210_gp,&
         303.63280_gp,  481.06690_gp,  441.15490_gp,  429.19400_gp,  418.38610_gp,&
         391.95070_gp,  498.02340_gp,  917.97160_gp,  528.98700_gp,  607.41450_gp,&
         544.18980_gp,  482.98360_gp,  465.22570_gp,  553.43860_gp,  130.11390_gp,&
         128.21910_gp,   94.71530_gp,   76.25300_gp,   52.66900_gp,  229.44250_gp,&
         276.99400_gp,  266.99740_gp,  233.14220_gp,  196.92250_gp,  321.80650_gp,&
         310.41930_gp,  314.34380_gp,  287.70260_gp,  218.15660_gp,  187.41200_gp/)
    vdwparams%coeffs(29601:29700)=(/  183.89760_gp,  214.24530_gp,  199.97120_gp,  280.72610_gp,  296.85620_gp,&
         273.45920_gp,  255.24180_gp,  385.75940_gp,  373.08750_gp,  376.81730_gp,&
         364.32260_gp,  323.48820_gp,  286.50850_gp,  271.11330_gp,  275.23510_gp,&
         287.20260_gp,  362.75140_gp,  394.18070_gp,  372.88760_gp,  360.59670_gp,&
         467.16080_gp,  485.37390_gp,  361.48780_gp,  376.94960_gp,  351.52420_gp,&
         320.26850_gp,  302.03680_gp,  307.06620_gp,  322.83420_gp,  388.66590_gp,&
         410.79890_gp,  419.30730_gp,  416.80560_gp,  514.46020_gp,   98.22350_gp,&
         88.65720_gp,   67.27040_gp,  217.43180_gp,  259.98580_gp,  216.72230_gp,&
         191.53170_gp,  209.21040_gp,  231.16180_gp,  266.79440_gp,  263.67810_gp,&
         310.83150_gp,   34.55030_gp,   22.88750_gp,  527.18510_gp,  307.88140_gp,&
         208.23790_gp,  141.04070_gp,   98.80250_gp,   74.90290_gp,   56.81510_gp,&
         43.75560_gp,  630.86870_gp,  489.75320_gp,  452.61920_gp,  357.74280_gp,&
         279.53920_gp,  232.24340_gp,  189.88240_gp,  155.46270_gp, 1031.53500_gp,&
         857.45370_gp,  709.44180_gp,  685.88510_gp,  628.54190_gp,  495.09360_gp,&
         541.71840_gp,  425.19040_gp,  451.43600_gp,  464.73810_gp,  356.24130_gp,&
         366.57850_gp,  433.86350_gp,  383.51890_gp,  327.67900_gp,  294.27420_gp,&
         257.71960_gp,  224.26550_gp, 1156.64470_gp, 1021.19920_gp,  896.68220_gp,&
         806.99680_gp,  736.48020_gp,  569.27660_gp,  634.92090_gp,  484.34860_gp,&
         529.53710_gp,  491.30410_gp,  409.43180_gp,  433.19920_gp,  542.84130_gp,&
         502.98910_gp,  448.97690_gp,  416.89450_gp,  377.22550_gp,  339.35050_gp/)
    vdwparams%coeffs(29701:29800)=(/ 1409.03750_gp, 1300.63490_gp, 1145.42360_gp,  530.13690_gp, 1153.37210_gp,&
         1107.99030_gp, 1080.33080_gp, 1054.88030_gp, 1032.31590_gp,  813.71020_gp,&
         912.55420_gp,  880.50980_gp,  931.68500_gp,  911.98360_gp,  894.24260_gp,&
         883.70520_gp,  745.78300_gp,  735.06570_gp,  671.61750_gp,  567.50780_gp,&
         577.64320_gp,  524.00420_gp,  480.23510_gp,  399.15430_gp,  373.02900_gp,&
         383.92690_gp,  557.45160_gp,  545.94690_gp,  502.42860_gp,  479.54670_gp,&
         442.91920_gp,  406.20580_gp, 1333.10000_gp, 1287.64570_gp, 1140.78480_gp,&
         1027.48830_gp, 1018.91580_gp,  986.63370_gp, 1014.44050_gp,  982.65960_gp,&
         55.56430_gp,  179.07060_gp,  228.55840_gp,  175.27450_gp,  132.25620_gp,&
         92.11540_gp,   67.18110_gp,   46.30300_gp,  261.54000_gp,  404.27630_gp,&
         410.41760_gp,  330.75310_gp,  271.34880_gp,  229.67560_gp,  187.90960_gp,&
         360.96750_gp,  685.80030_gp,  362.13940_gp,  349.61060_gp,  342.61200_gp,&
         339.34490_gp,  312.72950_gp,  289.53010_gp,  276.18840_gp,  269.98570_gp,&
         266.47580_gp,  250.50770_gp,  407.73530_gp,  358.95200_gp,  322.14500_gp,&
         294.38040_gp,  259.50110_gp,  430.88110_gp,  834.09790_gp,  640.47400_gp,&
         478.73550_gp,  483.44780_gp,  450.98400_gp,  504.76500_gp,  393.62700_gp,&
         368.27100_gp,  342.04040_gp,  330.83540_gp,  327.65190_gp,  517.80290_gp,&
         469.37510_gp,  440.10850_gp,  416.61070_gp,  380.76150_gp,  510.28570_gp,&
         1080.20490_gp,  810.22120_gp,  514.19410_gp,  503.56860_gp,  487.64360_gp,&
         489.91040_gp,  471.83730_gp,  492.98820_gp,  463.13550_gp,  471.01170_gp/)
    vdwparams%coeffs(29801:29900)=(/  439.69410_gp,  427.20060_gp,  424.43580_gp,  446.01630_gp,  410.29850_gp,&
         554.10220_gp,  514.78350_gp,  470.18380_gp,  475.32640_gp,  417.08260_gp,&
         392.48360_gp,  374.93960_gp,  358.43860_gp,  346.65580_gp,  556.06680_gp,&
         507.82080_gp,  492.95080_gp,  479.71720_gp,  448.28660_gp,  574.77940_gp,&
         1076.07440_gp,  608.43480_gp,  701.09350_gp,  627.07340_gp,  554.92490_gp,&
         534.13870_gp,  640.33600_gp,  149.76800_gp,  147.23850_gp,  107.96360_gp,&
         86.51740_gp,   59.30290_gp,  265.25750_gp,  320.39240_gp,  307.75700_gp,&
         267.62990_gp,  225.03380_gp,  371.71230_gp,  357.61290_gp,  362.04400_gp,&
         331.30930_gp,  249.89390_gp,  214.01180_gp,  210.06230_gp,  245.93880_gp,&
         229.23870_gp,  323.67220_gp,  341.84070_gp,  313.73420_gp,  292.07370_gp,&
         445.29880_gp,  429.22790_gp,  433.06580_gp,  418.69180_gp,  370.51430_gp,&
         327.26370_gp,  309.34100_gp,  314.85220_gp,  328.86980_gp,  417.60930_gp,&
         453.91950_gp,  428.12220_gp,  413.26530_gp,  538.67350_gp,  558.84130_gp,&
         414.80660_gp,  432.99810_gp,  402.55440_gp,  365.84940_gp,  344.33760_gp,&
         350.80960_gp,  368.98910_gp,  446.36650_gp,  471.83450_gp,  481.14100_gp,&
         477.72640_gp,  593.04720_gp,  112.13310_gp,  101.01310_gp,   76.21000_gp,&
         250.14210_gp,  299.82180_gp,  248.37820_gp,  219.15410_gp,  240.46050_gp,&
         265.25500_gp,  306.52410_gp,  302.26660_gp,  356.23590_gp,  409.14530_gp,&
         35.09390_gp,   23.37920_gp,  518.64270_gp,  307.44210_gp,  209.67450_gp,&
         142.83480_gp,  100.45000_gp,   76.34240_gp,   58.01760_gp,   44.74080_gp/)
    vdwparams%coeffs(29901:30000)=(/  621.19430_gp,  487.74180_gp,  452.76050_gp,  359.91590_gp,  282.55640_gp,&
         235.42680_gp,  192.98930_gp,  158.34200_gp, 1013.83230_gp,  850.06180_gp,&
         704.72450_gp,  682.61700_gp,  626.23920_gp,  493.43100_gp,  540.59460_gp,&
         424.51240_gp,  451.68930_gp,  464.46650_gp,  356.10950_gp,  367.66370_gp,&
         434.74170_gp,  385.97850_gp,  331.06180_gp,  298.00720_gp,  261.58890_gp,&
         228.09040_gp, 1137.80990_gp, 1012.15180_gp,  891.91400_gp,  804.49370_gp,&
         735.28120_gp,  569.79780_gp,  634.89910_gp,  485.65920_gp,  530.91130_gp,&
         493.00840_gp,  410.66850_gp,  435.17770_gp,  543.81980_gp,  505.66970_gp,&
         452.91870_gp,  421.41870_gp,  382.13760_gp,  344.44580_gp, 1386.58060_gp,&
         1287.77380_gp, 1138.10530_gp,  534.48010_gp, 1143.26360_gp, 1098.88180_gp,&
         1071.61370_gp, 1046.50460_gp, 1024.25000_gp,  810.59520_gp,  904.07560_gp,&
         872.95200_gp,  925.19070_gp,  905.70520_gp,  888.20360_gp,  877.58750_gp,&
         742.66050_gp,  734.64070_gp,  672.68970_gp,  569.41570_gp,  580.04130_gp,&
         527.13410_gp,  483.80280_gp,  402.72830_gp,  376.58500_gp,  387.84150_gp,&
         559.00180_gp,  548.91370_gp,  506.69310_gp,  484.46800_gp,  448.38730_gp,&
         411.97350_gp, 1315.58980_gp, 1277.30250_gp, 1135.23320_gp, 1026.29430_gp,&
         1015.89640_gp,  983.78580_gp, 1009.12300_gp,  977.91270_gp,   56.19460_gp,&
         179.02740_gp,  228.88940_gp,  176.77430_gp,  133.92400_gp,   93.65620_gp,&
         68.51410_gp,   47.38350_gp,  261.06680_gp,  403.30920_gp,  410.97880_gp,&
         333.03840_gp,  274.32020_gp,  232.80570_gp,  190.97240_gp,  361.49280_gp/)
    vdwparams%coeffs(30001:30100)=(/  680.31840_gp,  364.24560_gp,  351.70890_gp,  344.63160_gp,  341.15480_gp,&
         315.19770_gp,  292.03450_gp,  278.57600_gp,  272.25780_gp,  268.17480_gp,&
         253.01490_gp,  409.01750_gp,  361.62180_gp,  325.54390_gp,  298.10330_gp,&
         263.36490_gp,  431.94800_gp,  827.11060_gp,  639.60320_gp,  481.40710_gp,&
         486.07860_gp,  453.85780_gp,  506.29760_gp,  397.07420_gp,  371.57400_gp,&
         345.25330_gp,  333.80120_gp,  331.21190_gp,  519.38200_gp,  472.51890_gp,&
         444.17700_gp,  421.18710_gp,  385.70840_gp,  513.02910_gp, 1069.79550_gp,&
         808.69310_gp,  518.50620_gp,  507.80630_gp,  491.82970_gp,  493.86870_gp,&
         474.82940_gp,  496.89070_gp,  466.99360_gp,  474.50580_gp,  443.55070_gp,&
         431.01790_gp,  428.17260_gp,  449.62020_gp,  414.01670_gp,  555.86120_gp,&
         517.47490_gp,  473.62850_gp,  477.90620_gp,  421.36930_gp,  396.78090_gp,&
         379.18850_gp,  362.45620_gp,  351.13700_gp,  558.15810_gp,  511.51250_gp,&
         497.47420_gp,  484.77400_gp,  453.83680_gp,  578.00010_gp, 1068.12810_gp,&
         613.48930_gp,  704.81290_gp,  631.07130_gp,  559.73100_gp,  539.03440_gp,&
         642.11740_gp,  151.12150_gp,  148.73950_gp,  109.57100_gp,   87.99580_gp,&
         60.52090_gp,  266.67450_gp,  321.97480_gp,  310.12840_gp,  270.48000_gp,&
         228.11080_gp,  373.64230_gp,  360.26270_gp,  364.78380_gp,  333.75970_gp,&
         252.58170_gp,  216.75670_gp,  212.69690_gp,  248.19240_gp,  231.56140_gp,&
         325.93780_gp,  344.61890_gp,  317.16790_gp,  295.78510_gp,  447.75000_gp,&
         432.82190_gp,  437.05000_gp,  422.43200_gp,  374.66630_gp,  331.51190_gp/)
    vdwparams%coeffs(30101:30200)=(/  313.55000_gp,  318.46020_gp,  332.42910_gp,  420.67140_gp,  457.28600_gp,&
         432.33570_gp,  417.89910_gp,  542.09840_gp,  563.22280_gp,  419.08730_gp,&
         437.00180_gp,  407.13540_gp,  370.60120_gp,  349.25650_gp,  355.16130_gp,&
         373.50160_gp,  450.30970_gp,  476.07850_gp,  485.92250_gp,  482.91440_gp,&
         596.98270_gp,  113.75420_gp,  102.51480_gp,   77.54400_gp,  252.47130_gp,&
         301.96690_gp,  251.32430_gp,  221.65120_gp,  242.37210_gp,  268.09680_gp,&
         309.54100_gp,  305.75380_gp,  360.07570_gp,  413.00240_gp,  417.37900_gp,&
         34.04400_gp,   23.06960_gp,  489.87970_gp,  292.40300_gp,  200.96210_gp,&
         138.05320_gp,   97.88350_gp,   74.90960_gp,   57.31890_gp,   44.47810_gp,&
         587.37390_gp,  463.39590_gp,  431.42990_gp,  344.62380_gp,  272.03200_gp,&
         227.71370_gp,  187.62390_gp,  154.74190_gp,  959.74420_gp,  806.85000_gp,&
         669.51760_gp,  649.52300_gp,  596.39240_gp,  470.75810_gp,  515.51420_gp,&
         405.65730_gp,  431.61180_gp,  443.37980_gp,  340.75230_gp,  352.07290_gp,&
         415.43610_gp,  370.09550_gp,  318.76960_gp,  287.90150_gp,  253.71080_gp,&
         222.13540_gp, 1078.13540_gp,  961.12990_gp,  848.64220_gp,  766.69670_gp,&
         701.69050_gp,  545.56530_gp,  607.11050_gp,  466.12980_gp,  508.94760_gp,&
         473.08430_gp,  394.80510_gp,  418.18980_gp,  520.81600_gp,  485.32990_gp,&
         436.04710_gp,  406.66470_gp,  369.83570_gp,  334.41780_gp, 1314.19250_gp,&
         1222.72090_gp, 1082.57070_gp,  514.66410_gp, 1086.68480_gp, 1044.75740_gp,&
         1018.89630_gp,  995.06410_gp,  973.93920_gp,  773.28240_gp,  860.29910_gp/)
    vdwparams%coeffs(30201:30300)=(/  830.99980_gp,  880.07510_gp,  861.53880_gp,  844.91510_gp,  834.67880_gp,&
         707.86620_gp,  701.35370_gp,  643.46570_gp,  546.18300_gp,  556.58840_gp,&
         506.84330_gp,  466.02790_gp,  389.18560_gp,  364.40110_gp,  375.25010_gp,&
         536.71560_gp,  527.65960_gp,  488.26870_gp,  467.65830_gp,  433.90630_gp,&
         399.73560_gp, 1248.87860_gp, 1214.20110_gp, 1081.06360_gp,  980.16620_gp,&
         969.70180_gp,  939.22300_gp,  962.10980_gp,  932.57400_gp,   54.11210_gp,&
         170.52180_gp,  218.38760_gp,  169.84800_gp,  129.53340_gp,   91.36060_gp,&
         67.38490_gp,   47.13160_gp,  248.74970_gp,  383.82790_gp,  392.01850_gp,&
         319.22290_gp,  264.20730_gp,  225.19300_gp,  185.67430_gp,  346.42790_gp,&
         646.89750_gp,  349.69760_gp,  337.85900_gp,  331.10740_gp,  327.68960_gp,&
         303.24930_gp,  281.26900_gp,  268.41200_gp,  262.30550_gp,  258.06750_gp,&
         244.04160_gp,  391.23290_gp,  347.13060_gp,  313.56550_gp,  287.99690_gp,&
         255.39640_gp,  414.51570_gp,  786.76620_gp,  611.08840_gp,  462.40460_gp,&
         466.95430_gp,  436.60630_gp,  485.96370_gp,  383.00380_gp,  358.72510_gp,&
         333.67710_gp,  322.58500_gp,  320.32730_gp,  497.88730_gp,  454.13030_gp,&
         427.84290_gp,  406.48930_gp,  373.25390_gp,  492.96130_gp, 1017.35010_gp,&
         772.50690_gp,  499.42220_gp,  489.15000_gp,  473.87610_gp,  475.62120_gp,&
         456.85930_gp,  478.32350_gp,  449.82470_gp,  456.72830_gp,  427.39350_gp,&
         415.38260_gp,  412.56270_gp,  432.70380_gp,  398.97500_gp,  532.87050_gp,&
         497.12170_gp,  456.06110_gp,  459.49090_gp,  407.10880_gp,  383.88990_gp/)
    vdwparams%coeffs(30301:30400)=(/  367.23030_gp,  351.26140_gp,  340.73060_gp,  536.10400_gp,  492.51500_gp,&
         479.66590_gp,  468.03600_gp,  439.13860_gp,  555.39950_gp, 1017.00700_gp,&
         590.69510_gp,  677.32330_gp,  607.53220_gp,  539.86600_gp,  520.21880_gp,&
         616.73560_gp,  145.07840_gp,  143.27430_gp,  106.36560_gp,   85.94680_gp,&
         59.72270_gp,  255.19400_gp,  308.08810_gp,  297.53400_gp,  260.49720_gp,&
         220.70250_gp,  358.43680_gp,  346.20560_gp,  350.64040_gp,  321.04220_gp,&
         244.31950_gp,  210.30300_gp,  206.32120_gp,  239.60260_gp,  223.81460_gp,&
         312.90410_gp,  331.07290_gp,  305.65610_gp,  285.80050_gp,  429.95510_gp,&
         416.47270_gp,  420.87950_gp,  407.06770_gp,  362.22600_gp,  321.39710_gp,&
         304.35760_gp,  308.56230_gp,  321.73420_gp,  404.90970_gp,  439.80180_gp,&
         416.70500_gp,  403.41630_gp,  521.07880_gp,  541.60900_gp,  404.13150_gp,&
         421.32440_gp,  393.64640_gp,  359.24180_gp,  339.23690_gp,  344.53850_gp,&
         362.04930_gp,  434.55430_gp,  459.16510_gp,  468.85050_gp,  466.35820_gp,&
         573.92320_gp,  110.15230_gp,   99.61750_gp,   75.93040_gp,  242.51390_gp,&
         289.70710_gp,  242.38610_gp,  214.63710_gp,  233.84620_gp,  258.33420_gp,&
         297.94670_gp,  294.87070_gp,  347.80090_gp,  398.05300_gp,  402.61790_gp,&
         389.60690_gp,   39.20180_gp,   26.15420_gp,  592.30750_gp,  346.59360_gp,&
         235.11450_gp,  159.78570_gp,  112.31150_gp,   85.39020_gp,   64.95480_gp,&
         50.15430_gp,  709.06070_gp,  551.15190_gp,  509.90250_gp,  403.75640_gp,&
         316.18250_gp,  263.18820_gp,  215.63940_gp,  176.93360_gp, 1160.21240_gp/)
    vdwparams%coeffs(30401:30500)=(/  964.82980_gp,  798.48900_gp,  772.40680_gp,  708.04470_gp,  558.13800_gp,&
         610.52950_gp,  479.61190_gp,  509.14700_gp,  523.95450_gp,  402.03690_gp,&
         413.76860_gp,  489.32470_gp,  433.09510_gp,  370.65200_gp,  333.32050_gp,&
         292.38780_gp,  254.87000_gp, 1301.39110_gp, 1149.33690_gp, 1009.83130_gp,&
         909.33820_gp,  830.29650_gp,  642.61650_gp,  716.35260_gp,  547.25930_gp,&
         598.00230_gp,  555.03560_gp,  462.92020_gp,  489.66130_gp,  612.79550_gp,&
         568.24370_gp,  507.83300_gp,  471.98620_gp,  427.58360_gp,  385.15720_gp,&
         1585.58210_gp, 1463.90790_gp, 1289.91280_gp,  599.69740_gp, 1298.63240_gp,&
         1247.61030_gp, 1216.48250_gp, 1187.83300_gp, 1162.43110_gp,  917.33060_gp,&
         1028.05790_gp,  992.08130_gp, 1049.21800_gp, 1027.02310_gp, 1007.04730_gp,&
         995.11980_gp,  840.45040_gp,  828.79200_gp,  757.80230_gp,  641.04200_gp,&
         652.56880_gp,  592.43330_gp,  543.33580_gp,  452.19510_gp,  422.82310_gp,&
         435.14150_gp,  629.94840_gp,  617.18100_gp,  568.51380_gp,  542.99390_gp,&
         502.02460_gp,  460.91780_gp, 1500.82180_gp, 1449.84490_gp, 1285.17190_gp,&
         1158.72200_gp, 1148.89880_gp, 1112.57430_gp, 1143.40300_gp, 1107.66340_gp,&
         62.85300_gp,  201.70300_gp,  257.62390_gp,  198.09520_gp,  149.88210_gp,&
         104.75870_gp,   76.66340_gp,   53.08880_gp,  294.65140_gp,  455.26290_gp,&
         462.54200_gp,  373.45290_gp,  306.96940_gp,  260.28570_gp,  213.40440_gp,&
         407.61910_gp,  772.30000_gp,  409.17160_gp,  395.11390_gp,  387.22670_gp,&
         383.50230_gp,  353.62940_gp,  327.53820_gp,  312.49710_gp,  305.47070_gp/)
    vdwparams%coeffs(30501:30600)=(/  301.36610_gp,  283.55320_gp,  460.02380_gp,  405.53300_gp,  364.44410_gp,&
         333.44240_gp,  294.39140_gp,  486.83630_gp,  939.49220_gp,  722.46140_gp,&
         541.04620_gp,  546.39880_gp,  509.99070_gp,  570.35510_gp,  445.58690_gp,&
         417.03620_gp,  387.50290_gp,  374.79630_gp,  371.28680_gp,  584.73500_gp,&
         530.55010_gp,  497.90060_gp,  471.68750_gp,  431.57390_gp,  576.82000_gp,&
         1216.71620_gp,  913.91750_gp,  581.72760_gp,  569.72170_gp,  551.75870_gp,&
         554.22640_gp,  533.59620_gp,  557.60720_gp,  523.97000_gp,  532.73770_gp,&
         497.51580_gp,  483.40880_gp,  480.24390_gp,  504.41450_gp,  464.26970_gp,&
         625.77080_gp,  581.83450_gp,  531.90280_gp,  537.41510_gp,  472.45060_gp,&
         444.83850_gp,  425.12390_gp,  406.52810_gp,  393.35910_gp,  628.45220_gp,&
         574.44200_gp,  557.91660_gp,  543.22420_gp,  508.08860_gp,  649.72250_gp,&
         1212.50540_gp,  688.24310_gp,  792.52420_gp,  709.36470_gp,  628.16270_gp,&
         604.77390_gp,  723.74190_gp,  169.20360_gp,  166.58210_gp,  122.53210_gp,&
         98.44210_gp,   67.76700_gp,  299.29840_gp,  361.51620_gp,  347.61240_gp,&
         302.75900_gp,  255.04960_gp,  419.85730_gp,  404.20110_gp,  409.25130_gp,&
         374.62180_gp,  283.20230_gp,  242.83500_gp,  238.33220_gp,  278.49510_gp,&
         259.70610_gp,  365.69310_gp,  386.32290_gp,  355.00350_gp,  330.85020_gp,&
         503.18910_gp,  485.40080_gp,  489.89620_gp,  473.77460_gp,  419.81260_gp,&
         371.22640_gp,  351.07190_gp,  357.06620_gp,  372.78680_gp,  472.32770_gp,&
         513.22190_gp,  484.46450_gp,  467.94670_gp,  608.95860_gp,  631.82890_gp/)
    vdwparams%coeffs(30601:30700)=(/  469.49580_gp,  490.06850_gp,  456.12840_gp,  414.97090_gp,  390.89370_gp,&
         398.04290_gp,  418.52820_gp,  505.37700_gp,  534.08270_gp,  544.69760_gp,&
         541.02390_gp,  670.46940_gp,  127.13380_gp,  114.69320_gp,   86.80610_gp,&
         282.66470_gp,  338.66010_gp,  281.13750_gp,  248.47100_gp,  272.23020_gp,&
         300.12720_gp,  346.67580_gp,  342.12760_gp,  403.46440_gp,  462.98340_gp,&
         467.50670_gp,  451.16740_gp,  524.18400_gp,   39.74380_gp,   26.62780_gp,&
         583.21980_gp,  346.23830_gp,  236.62000_gp,  161.59220_gp,  113.93640_gp,&
         86.79140_gp,   66.11310_gp,   51.09610_gp,  698.76860_gp,  549.19660_gp,&
         510.16460_gp,  406.05930_gp,  319.27240_gp,  266.39090_gp,  218.72080_gp,&
         179.75600_gp, 1140.93070_gp,  957.09360_gp,  793.62090_gp,  769.03960_gp,&
         705.68530_gp,  556.35380_gp,  609.39430_gp,  478.86060_gp,  509.44850_gp,&
         523.72160_gp,  401.86070_gp,  414.92100_gp,  490.28010_gp,  435.66600_gp,&
         374.11390_gp,  337.09110_gp,  296.25040_gp,  258.64890_gp, 1280.78710_gp,&
         1139.76420_gp, 1004.83930_gp,  906.72570_gp,  829.02620_gp,  643.06490_gp,&
         716.26640_gp,  548.50410_gp,  599.37630_gp,  556.74870_gp,  464.07170_gp,&
         491.65550_gp,  613.77870_gp,  571.01580_gp,  511.87060_gp,  476.58080_gp,&
         432.52650_gp,  390.23930_gp, 1560.87280_gp, 1450.13560_gp, 1282.13000_gp,&
         604.10230_gp, 1287.81720_gp, 1237.89550_gp, 1207.19560_gp, 1178.91990_gp,&
         1153.85730_gp,  913.95230_gp, 1018.77680_gp,  983.79020_gp, 1042.35250_gp,&
         1020.39580_gp, 1000.68110_gp,  988.67650_gp,  837.14090_gp,  828.36360_gp/)
    vdwparams%coeffs(30701:30800)=(/  758.91210_gp,  642.92980_gp,  654.98530_gp,  595.58850_gp,  546.92440_gp,&
         455.73730_gp,  426.33520_gp,  439.04410_gp,  631.41780_gp,  620.17920_gp,&
         572.84430_gp,  547.97510_gp,  507.52580_gp,  466.68040_gp, 1481.52260_gp,&
         1438.75530_gp, 1279.28850_gp, 1157.41850_gp, 1145.60070_gp, 1109.45820_gp,&
         1137.68210_gp, 1102.56050_gp,   63.49760_gp,  201.70420_gp,  257.99550_gp,&
         199.63690_gp,  151.54880_gp,  106.26960_gp,   77.95280_gp,   54.11970_gp,&
         294.20730_gp,  454.32800_gp,  463.21330_gp,  375.84780_gp,  310.00460_gp,&
         263.43270_gp,  216.44080_gp,  408.08240_gp,  766.38090_gp,  411.33690_gp,&
         397.25380_gp,  389.28370_gp,  385.34090_gp,  356.16420_gp,  330.09770_gp,&
         314.92830_gp,  307.78390_gp,  303.08940_gp,  286.11150_gp,  461.38820_gp,&
         408.30150_gp,  367.91630_gp,  337.20180_gp,  298.24950_gp,  487.81130_gp,&
         931.86330_gp,  721.38820_gp,  543.73140_gp,  549.04170_gp,  512.86200_gp,&
         571.79580_gp,  449.04810_gp,  420.33540_gp,  390.69990_gp,  377.74640_gp,&
         374.86870_gp,  586.34250_gp,  533.78120_gp,  502.06010_gp,  476.33350_gp,&
         436.55340_gp,  579.50930_gp, 1205.23090_gp,  912.05970_gp,  586.09720_gp,&
         574.01590_gp,  555.99730_gp,  558.23140_gp,  536.60520_gp,  561.56960_gp,&
         527.88270_gp,  536.26980_gp,  501.42950_gp,  487.28280_gp,  484.03880_gp,&
         508.09720_gp,  468.04880_gp,  627.49410_gp,  584.51000_gp,  535.34350_gp,&
         539.96900_gp,  476.74360_gp,  449.12500_gp,  429.35080_gp,  410.50840_gp,&
         397.82250_gp,  630.50780_gp,  578.17950_gp,  562.50870_gp,  548.34310_gp/)
    vdwparams%coeffs(30801:30900)=(/  513.67500_gp,  652.89240_gp, 1203.70160_gp,  693.37750_gp,  796.21000_gp,&
         713.28770_gp,  632.98840_gp,  609.69380_gp,  725.37520_gp,  170.61220_gp,&
         168.10260_gp,  124.13050_gp,   99.89020_gp,   68.93890_gp,  300.82770_gp,&
         363.19850_gp,  350.08330_gp,  305.66690_gp,  258.14020_gp,  421.85200_gp,&
         406.92810_gp,  412.06730_gp,  377.11770_gp,  285.88310_gp,  245.56070_gp,&
         240.95060_gp,  280.76020_gp,  262.03990_gp,  368.03530_gp,  389.19220_gp,&
         358.50420_gp,  334.59570_gp,  505.67630_gp,  489.06840_gp,  493.95440_gp,&
         477.54800_gp,  423.96710_gp,  375.45330_gp,  355.25030_gp,  360.64340_gp,&
         376.33490_gp,  475.44030_gp,  516.67410_gp,  488.76370_gp,  472.64780_gp,&
         612.39630_gp,  636.29250_gp,  473.84040_gp,  494.09450_gp,  460.71470_gp,&
         419.70150_gp,  395.77380_gp,  402.34540_gp,  423.01710_gp,  509.34070_gp,&
         538.37690_gp,  549.54870_gp,  546.27450_gp,  674.41620_gp,  128.75890_gp,&
         116.17900_gp,   88.10270_gp,  285.08900_gp,  340.88270_gp,  284.13180_gp,&
         250.96550_gp,  274.14590_gp,  303.03570_gp,  349.76730_gp,  345.67500_gp,&
         407.35630_gp,  466.91700_gp,  471.96350_gp,  455.75440_gp,  528.76490_gp,&
         533.88100_gp,    9.20920_gp,    6.31810_gp,  120.66280_gp,   75.13110_gp,&
         52.96060_gp,   37.00420_gp,   26.50830_gp,   20.39980_gp,   15.65780_gp,&
         12.16080_gp,  145.00200_gp,  118.07480_gp,  111.43770_gp,   90.57400_gp,&
         72.51450_gp,   61.21840_gp,   50.81020_gp,   42.13240_gp,  236.18020_gp,&
         203.06790_gp,  169.41860_gp,  165.29120_gp,  152.25000_gp,  120.32780_gp/)
    vdwparams%coeffs(30901:31000)=(/  132.19840_gp,  104.19900_gp,  111.50630_gp,  114.15060_gp,   87.80520_gp,&
         91.57280_gp,  107.87050_gp,   97.34900_gp,   84.82980_gp,   77.15330_gp,&
         68.44160_gp,   60.25440_gp,  266.10180_gp,  241.85040_gp,  215.72630_gp,&
         196.15430_gp,  180.30420_gp,  141.25220_gp,  156.72560_gp,  121.28760_gp,&
         132.32810_gp,  123.29130_gp,  102.77530_gp,  109.28870_gp,  135.08140_gp,&
         127.15210_gp,  115.40180_gp,  108.29220_gp,   99.11090_gp,   90.12960_gp,&
         324.98720_gp,  306.96860_gp,  274.49320_gp,  136.02870_gp,  273.57430_gp,&
         263.38860_gp,  256.96540_gp,  251.03460_gp,  245.78320_gp,  197.45530_gp,&
         216.47080_gp,  209.53600_gp,  222.58700_gp,  217.94160_gp,  213.81010_gp,&
         211.10760_gp,  180.48560_gp,  180.74340_gp,  166.87410_gp,  142.41400_gp,&
         145.42480_gp,  133.10380_gp,  122.87620_gp,  103.04210_gp,   96.61300_gp,&
         99.65610_gp,  139.64810_gp,  138.28410_gp,  129.08340_gp,  124.28020_gp,&
         116.00970_gp,  107.44210_gp,  311.32220_gp,  306.50670_gp,  275.21710_gp,&
         252.16490_gp,  248.21560_gp,  240.44920_gp,  244.57260_gp,  237.32770_gp,&
         14.45410_gp,   43.97460_gp,   56.65270_gp,   44.98880_gp,   34.72260_gp,&
         24.75250_gp,   18.38310_gp,   12.92910_gp,   63.81230_gp,   98.37240_gp,&
         101.58890_gp,   84.11580_gp,   70.46290_gp,   60.52550_gp,   50.27080_gp,&
         89.85690_gp,  163.35340_gp,   91.79510_gp,   88.73890_gp,   86.92640_gp,&
         85.88060_gp,   80.03750_gp,   74.39300_gp,   70.98710_gp,   69.31940_gp,&
         67.78350_gp,   64.77980_gp,  101.92230_gp,   91.58950_gp,   83.50230_gp/)
    vdwparams%coeffs(31001:31100)=(/   77.17000_gp,   68.87140_gp,  107.86680_gp,  198.56600_gp,  157.35980_gp,&
         121.38420_gp,  122.50750_gp,  114.84620_gp,  126.65710_gp,  101.38250_gp,&
         94.99870_gp,   88.45150_gp,   85.37730_gp,   85.22740_gp,  129.59790_gp,&
         119.45590_gp,  113.38670_gp,  108.28980_gp,  100.02100_gp,  129.44610_gp,&
         256.05470_gp,  198.71970_gp,  132.07060_gp,  129.36140_gp,  125.38070_gp,&
         125.66750_gp,  120.08750_gp,  126.32050_gp,  118.91980_gp,  120.44740_gp,&
         113.13110_gp,  110.00160_gp,  109.21440_gp,  114.30830_gp,  105.69650_gp,&
         138.94310_gp,  130.37550_gp,  120.30200_gp,  120.54080_gp,  108.24460_gp,&
         102.24360_gp,   97.88980_gp,   93.57480_gp,   91.19300_gp,  139.86370_gp,&
         129.75270_gp,  127.06290_gp,  124.48510_gp,  117.42660_gp,  145.93990_gp,&
         257.58460_gp,  156.18080_gp,  177.63750_gp,  159.84520_gp,  142.84510_gp,&
         137.82260_gp,  160.54880_gp,   38.47450_gp,   38.13980_gp,   28.67630_gp,&
         23.29000_gp,   16.29870_gp,   66.86960_gp,   80.70050_gp,   78.58350_gp,&
         69.41150_gp,   59.32390_gp,   93.91940_gp,   91.30720_gp,   92.51540_gp,&
         84.65700_gp,   65.03490_gp,   56.29210_gp,   55.16980_gp,   63.45430_gp,&
         59.42470_gp,   82.44470_gp,   87.52820_gp,   81.48960_gp,   76.59390_gp,&
         112.80930_gp,  110.17370_gp,  111.60800_gp,  107.87110_gp,   96.60300_gp,&
         86.12760_gp,   81.68530_gp,   82.28870_gp,   85.62810_gp,  106.69460_gp,&
         115.92860_gp,  110.62760_gp,  107.54340_gp,  137.12500_gp,  143.10680_gp,&
         107.50530_gp,  111.75770_gp,  105.03500_gp,   96.27590_gp,   91.23250_gp/)
    vdwparams%coeffs(31101:31200)=(/   92.14580_gp,   96.75910_gp,  114.99010_gp,  121.56000_gp,  124.47110_gp,&
         124.15670_gp,  151.17960_gp,   29.68040_gp,   26.85060_gp,   20.58860_gp,&
         64.33320_gp,   76.40380_gp,   64.81270_gp,   57.19870_gp,   61.69740_gp,&
         68.79350_gp,   79.16560_gp,   78.74780_gp,   92.53550_gp,  105.52950_gp,&
         107.12000_gp,  103.85430_gp,  119.70710_gp,  121.29570_gp,   28.03150_gp,&
         7.36620_gp,    5.25670_gp,   87.81830_gp,   56.46620_gp,   40.89620_gp,&
         29.28300_gp,   21.41990_gp,   16.75440_gp,   13.05460_gp,   10.27160_gp,&
         105.92480_gp,   88.26200_gp,   84.30730_gp,   69.72240_gp,   56.78520_gp,&
         48.56810_gp,   40.85160_gp,   34.30500_gp,  172.96850_gp,  150.82030_gp,&
         126.34910_gp,  124.00130_gp,  114.58700_gp,   91.02330_gp,   99.97680_gp,&
         79.26300_gp,   84.95470_gp,   86.66010_gp,   67.08230_gp,   70.28200_gp,&
         82.32480_gp,   75.20760_gp,   66.41470_gp,   60.99130_gp,   54.68240_gp,&
         48.65000_gp,  195.59960_gp,  179.83950_gp,  161.77380_gp,  148.00690_gp,&
         136.70500_gp,  108.23560_gp,  119.59950_gp,   93.62690_gp,  101.81940_gp,&
         95.16120_gp,   79.66570_gp,   84.71390_gp,  103.61980_gp,   98.34310_gp,&
         90.17970_gp,   85.23430_gp,   78.66860_gp,   72.16170_gp,  239.24420_gp,&
         228.05560_gp,  205.54070_gp,  106.30350_gp,  204.00390_gp,  196.61650_gp,&
         191.87330_gp,  187.48180_gp,  183.59400_gp,  149.29180_gp,  161.90910_gp,&
         156.98310_gp,  166.53800_gp,  163.07020_gp,  160.00710_gp,  157.88900_gp,&
         136.08820_gp,  137.29240_gp,  127.62070_gp,  109.83770_gp,  112.32740_gp/)
    vdwparams%coeffs(31201:31300)=(/  103.46160_gp,   96.03330_gp,   81.23630_gp,   76.42820_gp,   78.85000_gp,&
         107.89760_gp,  107.36650_gp,  101.05710_gp,   97.83320_gp,   91.99900_gp,&
         85.84320_gp,  230.75540_gp,  228.82200_gp,  206.94660_gp,  191.64430_gp,&
         188.09650_gp,  182.30090_gp,  184.37780_gp,  179.08370_gp,   11.32990_gp,&
         33.21580_gp,   43.06200_gp,   34.99900_gp,   27.52060_gp,   20.04680_gp,&
         15.17510_gp,   10.93910_gp,   48.16220_gp,   74.00280_gp,   77.13620_gp,&
         64.96720_gp,   55.23750_gp,   48.02210_gp,   40.42090_gp,   69.04340_gp,&
         122.07070_gp,   71.08540_gp,   68.82970_gp,   67.43830_gp,   66.55680_gp,&
         62.39310_gp,   58.17760_gp,   55.56540_gp,   54.23970_gp,   52.79840_gp,&
         50.88560_gp,   78.05050_gp,   71.01430_gp,   65.44140_gp,   61.00410_gp,&
         55.00290_gp,   83.24180_gp,  148.51410_gp,  119.76540_gp,   94.11730_gp,&
         95.00090_gp,   89.42290_gp,   97.84170_gp,   79.57500_gp,   74.72870_gp,&
         69.77300_gp,   67.30700_gp,   67.40280_gp,   99.75220_gp,   92.80210_gp,&
         88.74340_gp,   85.26570_gp,   79.37160_gp,  100.47420_gp,  191.24940_gp,&
         151.14110_gp,  103.29330_gp,  101.19170_gp,   98.14670_gp,   98.22700_gp,&
         93.51930_gp,   98.62130_gp,   93.00700_gp,   93.97710_gp,   88.58220_gp,&
         86.17430_gp,   85.51070_gp,   89.19150_gp,   82.79860_gp,  107.01830_gp,&
         101.10190_gp,   93.96010_gp,   93.66720_gp,   85.39250_gp,   80.95520_gp,&
         77.69970_gp,   74.37630_gp,   72.78320_gp,  108.24170_gp,  101.28820_gp,&
         99.66790_gp,   98.05750_gp,   93.10690_gp,  113.30310_gp,  193.37770_gp/)
    vdwparams%coeffs(31301:31400)=(/  122.04020_gp,  137.87680_gp,  124.70980_gp,  112.09800_gp,  108.34540_gp,&
         124.18240_gp,   29.88190_gp,   29.89160_gp,   22.95170_gp,   18.91720_gp,&
         13.55250_gp,   51.35250_gp,   61.96100_gp,   60.86870_gp,   54.39410_gp,&
         47.09430_gp,   72.55240_gp,   70.96700_gp,   71.95950_gp,   65.94730_gp,&
         51.46100_gp,   44.91790_gp,   43.99020_gp,   49.91480_gp,   46.90710_gp,&
         63.91170_gp,   68.03430_gp,   63.96210_gp,   60.57230_gp,   87.38910_gp,&
         85.96680_gp,   87.31160_gp,   84.50410_gp,   76.38690_gp,   68.62370_gp,&
         65.29050_gp,   65.39490_gp,   67.83630_gp,   83.22700_gp,   90.26830_gp,&
         86.75740_gp,   84.74250_gp,  106.58000_gp,  111.46330_gp,   84.44510_gp,&
         87.68090_gp,   83.08050_gp,   76.68670_gp,   73.06660_gp,   73.48530_gp,&
         77.01240_gp,   90.34670_gp,   95.39730_gp,   97.85170_gp,   97.88170_gp,&
         117.58210_gp,   23.62100_gp,   21.53770_gp,   16.81690_gp,   50.03970_gp,&
         59.15830_gp,   51.00480_gp,   45.38910_gp,   48.41490_gp,   53.96170_gp,&
         61.90380_gp,   61.94380_gp,   72.97030_gp,   82.71700_gp,   84.21770_gp,&
         82.28250_gp,   94.13150_gp,   95.60180_gp,   22.37970_gp,   18.20670_gp,&
         21.00640_gp,   13.85470_gp,  308.79530_gp,  184.25660_gp,  125.77640_gp,&
         85.52230_gp,   59.93190_gp,   45.36880_gp,   34.31870_gp,   26.33580_gp,&
         369.77170_gp,  291.87160_gp,  271.32020_gp,  215.89980_gp,  169.43560_gp,&
         140.99550_gp,  115.35420_gp,   94.40270_gp,  602.21160_gp,  507.27030_gp,&
         420.88950_gp,  407.81600_gp,  374.20690_gp,  294.50370_gp,  323.09480_gp/)
    vdwparams%coeffs(31401:31500)=(/  253.39650_gp,  270.08810_gp,  277.67940_gp,  212.52810_gp,  219.88700_gp,&
         260.30910_gp,  231.33020_gp,  198.40790_gp,  178.47990_gp,  156.47810_gp,&
         136.21020_gp,  675.86350_gp,  603.71340_gp,  532.68760_gp,  480.73180_gp,&
         439.42810_gp,  340.34780_gp,  379.31480_gp,  289.93610_gp,  317.20640_gp,&
         294.52200_gp,  244.86850_gp,  259.88150_gp,  324.93450_gp,  302.48790_gp,&
         271.03870_gp,  252.16040_gp,  228.53700_gp,  205.80950_gp,  823.82010_gp,&
         767.60950_gp,  679.35020_gp,  319.68540_gp,  681.42510_gp,  655.12090_gp,&
         638.90630_gp,  623.97630_gp,  610.74860_gp,  483.68880_gp,  538.20120_gp,&
         519.82010_gp,  551.87870_gp,  540.29010_gp,  529.89220_gp,  523.55130_gp,&
         443.29500_gp,  439.27910_gp,  402.32530_gp,  340.31200_gp,  346.76670_gp,&
         315.08760_gp,  289.09700_gp,  240.31030_gp,  224.56090_gp,  231.41000_gp,&
         333.60810_gp,  327.96730_gp,  302.90300_gp,  289.64470_gp,  268.00240_gp,&
         246.08550_gp,  782.45740_gp,  761.83620_gp,  677.80430_gp,  613.15740_gp,&
         606.32200_gp,  587.10090_gp,  601.67790_gp,  583.14710_gp,   33.71260_gp,&
         107.27210_gp,  137.17830_gp,  105.98260_gp,   80.14270_gp,   55.83940_gp,&
         40.66350_gp,   27.90860_gp,  156.15450_gp,  241.39030_gp,  246.30130_gp,&
         199.75300_gp,  164.46690_gp,  139.40860_gp,  114.13580_gp,  215.87970_gp,&
         405.67650_gp,  217.96280_gp,  210.39040_gp,  206.10620_gp,  203.96330_gp,&
         188.58340_gp,  174.67580_gp,  166.56740_gp,  162.76320_gp,  160.19460_gp,&
         151.31410_gp,  244.94460_gp,  216.72460_gp,  195.08440_gp,  178.53020_gp/)
    vdwparams%coeffs(31501:31600)=(/  157.54230_gp,  257.90200_gp,  492.95740_gp,  382.05370_gp,  287.94440_gp,&
         290.66410_gp,  271.29600_gp,  302.36440_gp,  237.27980_gp,  221.91000_gp,&
         206.06750_gp,  199.15530_gp,  197.79040_gp,  310.41220_gp,  282.66840_gp,&
         265.80530_gp,  252.02720_gp,  230.68720_gp,  306.78390_gp,  637.15160_gp,&
         482.96010_gp,  310.10960_gp,  303.69930_gp,  294.12610_gp,  295.33450_gp,&
         283.74430_gp,  297.20650_gp,  279.27050_gp,  283.72210_gp,  265.27330_gp,&
         257.77900_gp,  256.08900_gp,  269.01820_gp,  247.64990_gp,  332.34580_gp,&
         309.38790_gp,  283.12170_gp,  285.56720_gp,  251.78560_gp,  236.94270_gp,&
         226.31410_gp,  216.17050_gp,  209.48600_gp,  333.26920_gp,  305.67610_gp,&
         297.42190_gp,  289.85610_gp,  271.29160_gp,  345.69720_gp,  636.70530_gp,&
         367.00470_gp,  421.26010_gp,  376.98280_gp,  334.44420_gp,  322.04370_gp,&
         383.13380_gp,   90.69560_gp,   89.11900_gp,   65.48790_gp,   52.41360_gp,&
         35.82470_gp,  159.91150_gp,  193.05340_gp,  185.99960_gp,  162.12890_gp,&
         136.57580_gp,  223.61620_gp,  215.70270_gp,  218.38960_gp,  199.66610_gp,&
         150.83440_gp,  129.33240_gp,  126.89320_gp,  148.22590_gp,  138.26350_gp,&
         195.24610_gp,  206.52060_gp,  190.04830_gp,  177.12870_gp,  267.85250_gp,&
         259.10620_gp,  261.64690_gp,  252.70610_gp,  223.95860_gp,  197.99630_gp,&
         187.16590_gp,  190.02130_gp,  198.42220_gp,  251.47160_gp,  273.56280_gp,&
         258.71400_gp,  250.05790_gp,  324.29550_gp,  337.22250_gp,  250.84450_gp,&
         261.38080_gp,  243.38630_gp,  221.37160_gp,  208.49420_gp,  211.88810_gp/)
    vdwparams%coeffs(31601:31700)=(/  222.92680_gp,  268.96090_gp,  284.51810_gp,  290.53760_gp,  288.76380_gp,&
         357.18170_gp,   68.11060_gp,   61.21690_gp,   46.10310_gp,  151.46030_gp,&
         181.02370_gp,  150.63080_gp,  132.28840_gp,  144.66440_gp,  160.60910_gp,&
         185.47740_gp,  183.19150_gp,  215.24240_gp,  247.05720_gp,  249.75520_gp,&
         240.45180_gp,  279.43950_gp,  282.20780_gp,   64.10830_gp,   50.20530_gp,&
         149.77340_gp,   23.17840_gp,   15.48010_gp,  330.56660_gp,  199.40340_gp,&
         137.26170_gp,   94.03570_gp,   66.32930_gp,   50.47410_gp,   38.36970_gp,&
         29.57400_gp,  396.24310_gp,  315.27200_gp,  294.19700_gp,  235.39440_gp,&
         185.72260_gp,  155.16910_gp,  127.47850_gp,  104.74120_gp,  645.29730_gp,&
         546.47840_gp,  454.07710_gp,  440.76620_gp,  404.85430_gp,  319.00750_gp,&
         350.08820_gp,  274.96140_gp,  293.35770_gp,  301.26690_gp,  230.92580_gp,&
         239.39590_gp,  282.94690_gp,  252.45910_gp,  217.44760_gp,  196.19620_gp,&
         172.58230_gp,  150.72600_gp,  724.93940_gp,  650.48490_gp,  575.57310_gp,&
         520.45920_gp,  476.44890_gp,  370.16770_gp,  412.05940_gp,  316.04420_gp,&
         345.50510_gp,  321.11170_gp,  267.21350_gp,  283.72370_gp,  353.60350_gp,&
         330.13010_gp,  296.81450_gp,  276.77860_gp,  251.52290_gp,  227.13240_gp,&
         883.97990_gp,  826.66440_gp,  733.57140_gp,  350.01970_gp,  734.68830_gp,&
         706.59670_gp,  689.17900_gp,  673.12870_gp,  658.90990_gp,  523.79680_gp,&
         580.56100_gp,  561.04910_gp,  595.75320_gp,  583.26490_gp,  572.08340_gp,&
         565.13830_gp,  479.72770_gp,  476.62680_gp,  437.46730_gp,  370.93970_gp/)
    vdwparams%coeffs(31701:31800)=(/  378.19490_gp,  344.33520_gp,  316.47760_gp,  263.74820_gp,  246.71740_gp,&
         254.30200_gp,  363.74800_gp,  358.27790_gp,  331.83780_gp,  317.89300_gp,&
         294.84760_gp,  271.38550_gp,  841.45020_gp,  821.71730_gp,  732.89850_gp,&
         665.26380_gp,  657.09130_gp,  636.34590_gp,  650.90570_gp,  631.06470_gp,&
         36.96640_gp,  116.25770_gp,  148.94110_gp,  115.91450_gp,   88.15320_gp,&
         61.83990_gp,   45.31410_gp,   31.36030_gp,  169.14050_gp,  261.22720_gp,&
         267.36160_gp,  218.00740_gp,  180.33050_gp,  153.42250_gp,  126.13260_gp,&
         234.97420_gp,  437.67790_gp,  237.97060_gp,  229.80530_gp,  225.13260_gp,&
         222.70380_gp,  206.34830_gp,  191.31890_gp,  182.47980_gp,  178.28580_gp,&
         175.18440_gp,  165.97640_gp,  266.53500_gp,  236.78670_gp,  213.87070_gp,&
         196.24890_gp,  173.73160_gp,  281.07340_gp,  531.87620_gp,  414.62030_gp,&
         314.43980_gp,  317.40770_gp,  296.62500_gp,  329.66960_gp,  260.11760_gp,&
         243.41820_gp,  226.22740_gp,  218.58540_gp,  217.37600_gp,  338.18310_gp,&
         308.94110_gp,  291.23020_gp,  276.66850_gp,  253.87120_gp,  335.05240_gp,&
         686.99570_gp,  523.95750_gp,  339.62030_gp,  332.61730_gp,  322.20350_gp,&
         323.36760_gp,  310.25500_gp,  325.31030_gp,  305.84570_gp,  310.46310_gp,&
         290.63460_gp,  282.47120_gp,  280.57190_gp,  294.43330_gp,  271.37930_gp,&
         362.12170_gp,  337.83560_gp,  309.86550_gp,  311.99870_gp,  276.47840_gp,&
         260.47040_gp,  248.97200_gp,  237.89260_gp,  230.89260_gp,  363.61850_gp,&
         334.53020_gp,  326.04840_gp,  318.19910_gp,  298.45520_gp,  377.59560_gp/)
    vdwparams%coeffs(31801:31900)=(/  687.70170_gp,  401.83270_gp,  460.11490_gp,  412.38990_gp,  366.61110_gp,&
         353.22130_gp,  417.84880_gp,   99.17000_gp,   97.68800_gp,   72.25950_gp,&
         58.09740_gp,   40.01320_gp,  174.23120_gp,  210.29720_gp,  203.17830_gp,&
         177.73730_gp,  150.32780_gp,  243.95900_gp,  235.81030_gp,  238.80070_gp,&
         218.40030_gp,  165.77460_gp,  142.52950_gp,  139.80680_gp,  162.60680_gp,&
         151.84970_gp,  213.28870_gp,  225.81390_gp,  208.45000_gp,  194.73060_gp,&
         292.42890_gp,  283.59490_gp,  286.62000_gp,  276.89930_gp,  246.11820_gp,&
         218.11560_gp,  206.39140_gp,  209.11520_gp,  218.15300_gp,  275.16580_gp,&
         299.21170_gp,  283.64510_gp,  274.57510_gp,  354.40230_gp,  368.86820_gp,&
         275.14660_gp,  286.53550_gp,  267.50210_gp,  243.84710_gp,  230.06870_gp,&
         233.44450_gp,  245.47400_gp,  294.94270_gp,  311.92660_gp,  318.74160_gp,&
         317.09360_gp,  390.42810_gp,   75.03160_gp,   67.58540_gp,   51.18740_gp,&
         165.69080_gp,  197.69460_gp,  165.35710_gp,  145.52910_gp,  158.55360_gp,&
         176.12940_gp,  203.19220_gp,  201.06830_gp,  236.37610_gp,  270.80520_gp,&
         274.04840_gp,  264.42240_gp,  306.57470_gp,  309.87480_gp,   70.69170_gp,&
         55.68570_gp,  164.18110_gp,  180.29120_gp,   31.60200_gp,   21.21540_gp,&
         455.91990_gp,  273.00460_gp,  187.38080_gp,  128.31130_gp,   90.61840_gp,&
         69.09390_gp,   52.66610_gp,   40.71950_gp,  546.50960_gp,  432.37650_gp,&
         402.62520_gp,  321.43630_gp,  253.31660_gp,  211.63180_gp,  173.95170_gp,&
         143.07940_gp,  891.20910_gp,  751.45760_gp,  623.82770_gp,  605.13210_gp/)
    vdwparams%coeffs(31901:32000)=(/  555.62000_gp,  438.06280_gp,  480.22450_gp,  377.40640_gp,  402.04340_gp,&
         413.05650_gp,  316.92740_gp,  327.86560_gp,  387.24130_gp,  344.90650_gp,&
         296.75050_gp,  267.67240_gp,  235.47970_gp,  205.76180_gp, 1000.93190_gp,&
         894.70560_gp,  790.39380_gp,  714.09530_gp,  653.41870_gp,  507.49190_gp,&
         565.00310_gp,  433.25440_gp,  473.45910_gp,  439.98550_gp,  366.59430_gp,&
         388.76210_gp,  484.65320_gp,  451.75780_gp,  405.68490_gp,  378.09630_gp,&
         343.49120_gp,  310.18060_gp, 1220.04590_gp, 1137.61460_gp, 1007.85070_gp,&
         478.59810_gp, 1010.92400_gp,  972.04260_gp,  948.02150_gp,  925.88830_gp,&
         906.27400_gp,  719.42080_gp,  799.40550_gp,  772.27150_gp,  819.10080_gp,&
         801.89050_gp,  786.45970_gp,  776.95490_gp,  658.87080_gp,  653.29340_gp,&
         599.21110_gp,  508.05680_gp,  517.81620_gp,  471.29410_gp,  433.09810_gp,&
         361.12210_gp,  337.90780_gp,  348.12330_gp,  498.76900_gp,  490.62310_gp,&
         453.90390_gp,  434.58510_gp,  402.90830_gp,  370.79470_gp, 1159.90290_gp,&
         1129.87190_gp, 1006.48620_gp,  912.46220_gp,  902.19930_gp,  873.76550_gp,&
         894.81040_gp,  867.39140_gp,   50.39400_gp,  159.13760_gp,  203.72270_gp,&
         158.21370_gp,  120.31900_gp,   84.51380_gp,   62.06650_gp,   43.14320_gp,&
         231.90190_gp,  358.00890_gp,  365.77110_gp,  297.64030_gp,  245.97730_gp,&
         209.27130_gp,  172.13180_gp,  322.12730_gp,  601.86910_gp,  325.50230_gp,&
         314.38040_gp,  308.05200_gp,  304.83850_gp,  282.15030_gp,  261.59540_gp,&
         249.56670_gp,  243.87510_gp,  239.89280_gp,  226.89510_gp,  364.63660_gp/)
    vdwparams%coeffs(32001:32100)=(/  323.40730_gp,  291.86440_gp,  267.75420_gp,  237.05490_gp,  385.25160_gp,&
         731.62570_gp,  568.61590_gp,  430.18180_gp,  434.34220_gp,  405.90280_gp,&
         451.73170_gp,  355.82860_gp,  333.09730_gp,  309.66650_gp,  299.33180_gp,&
         297.36580_gp,  463.29930_gp,  422.59500_gp,  398.00110_gp,  377.92500_gp,&
         346.68560_gp,  458.39400_gp,  945.52740_gp,  718.67700_gp,  464.37650_gp,&
         454.81020_gp,  440.56910_gp,  442.22540_gp,  424.69280_gp,  444.84260_gp,&
         418.24050_gp,  424.68210_gp,  397.37350_gp,  386.19540_gp,  383.60070_gp,&
         402.53610_gp,  370.98070_gp,  495.82350_gp,  462.34390_gp,  423.90050_gp,&
         427.13880_gp,  378.05750_gp,  356.25810_gp,  340.62570_gp,  325.64110_gp,&
         315.86290_gp,  498.32750_gp,  457.83400_gp,  445.87260_gp,  434.94050_gp,&
         407.80020_gp,  516.50850_gp,  945.54140_gp,  549.36590_gp,  629.82470_gp,&
         564.51290_gp,  501.56450_gp,  483.22840_gp,  572.96750_gp,  135.26050_gp,&
         133.32140_gp,   98.65400_gp,   79.45220_gp,   54.90050_gp,  238.05490_gp,&
         287.33370_gp,  277.34590_gp,  242.49400_gp,  205.06870_gp,  333.72060_gp,&
         322.29080_gp,  326.38510_gp,  298.65910_gp,  226.74830_gp,  194.95560_gp,&
         191.27030_gp,  222.52350_gp,  207.78550_gp,  291.43960_gp,  308.37570_gp,&
         284.45060_gp,  265.69090_gp,  400.06500_gp,  387.51810_gp,  391.55320_gp,&
         378.46930_gp,  336.35920_gp,  298.11220_gp,  282.14800_gp,  286.12510_gp,&
         298.49500_gp,  376.50510_gp,  409.20080_gp,  387.57520_gp,  375.04230_gp,&
         484.68920_gp,  504.02900_gp,  375.81480_gp,  391.63870_gp,  365.54820_gp/)
    vdwparams%coeffs(32101:32200)=(/  333.25000_gp,  314.43610_gp,  319.35010_gp,  335.74880_gp,  403.62310_gp,&
         426.69190_gp,  435.76920_gp,  433.37030_gp,  533.85140_gp,  102.32470_gp,&
         92.32600_gp,   70.08000_gp,  226.05040_gp,  269.96790_gp,  225.54180_gp,&
         199.11820_gp,  217.13920_gp,  240.41780_gp,  277.37170_gp,  274.34930_gp,&
         323.15090_gp,  370.16740_gp,  374.39290_gp,  361.63550_gp,  419.24330_gp,&
         423.53200_gp,   96.41790_gp,   76.08100_gp,  223.92830_gp,  245.99220_gp,&
         336.10200_gp,   36.78360_gp,   24.79260_gp,  531.22450_gp,  317.32060_gp,&
         217.79730_gp,  149.29680_gp,  105.60300_gp,   80.64440_gp,   61.57370_gp,&
         47.68440_gp,  636.82150_gp,  502.80810_gp,  468.01840_gp,  373.58140_gp,&
         294.52850_gp,  246.23230_gp,  202.58030_gp,  166.81140_gp, 1039.50150_gp,&
         874.85890_gp,  726.03180_gp,  704.21370_gp,  646.55260_gp,  510.01320_gp,&
         558.78150_gp,  439.38300_gp,  467.73930_gp,  480.56400_gp,  368.99340_gp,&
         381.42390_gp,  450.32390_gp,  401.00460_gp,  345.08980_gp,  311.40410_gp,&
         274.12300_gp,  239.71550_gp, 1167.51810_gp, 1041.87140_gp,  919.96850_gp,&
         831.02350_gp,  760.40220_gp,  590.76540_gp,  657.61940_gp,  504.46970_gp,&
         551.07490_gp,  512.14360_gp,  427.04150_gp,  452.57800_gp,  564.05260_gp,&
         525.57130_gp,  471.95910_gp,  439.93100_gp,  399.79760_gp,  361.19780_gp,&
         1423.05700_gp, 1325.12730_gp, 1173.36690_gp,  556.92770_gp, 1177.54000_gp,&
         1132.14840_gp, 1104.13990_gp, 1078.33210_gp, 1055.45830_gp,  837.68840_gp,&
         931.65810_gp,  899.94220_gp,  953.79600_gp,  933.72910_gp,  915.73150_gp/)
    vdwparams%coeffs(32201:32300)=(/  904.66320_gp,  767.04550_gp,  760.10300_gp,  697.15900_gp,  591.32700_gp,&
         602.61270_gp,  548.53340_gp,  504.15770_gp,  420.63020_gp,  393.69250_gp,&
         405.49310_gp,  580.79990_gp,  571.06960_gp,  528.26530_gp,  505.80110_gp,&
         469.02930_gp,  431.79360_gp, 1352.39730_gp, 1315.84340_gp, 1171.64350_gp,&
         1061.96930_gp, 1050.42960_gp, 1017.36200_gp, 1042.17170_gp, 1010.18370_gp,&
         58.58750_gp,  184.99140_gp,  236.83180_gp,  183.95600_gp,  140.03800_gp,&
         98.52470_gp,   72.48720_gp,   50.52220_gp,  269.73770_gp,  416.33670_gp,&
         425.18440_gp,  345.95860_gp,  286.02000_gp,  243.49690_gp,  200.46870_gp,&
         375.03550_gp,  700.97800_gp,  378.67460_gp,  365.78460_gp,  358.44930_gp,&
         354.74250_gp,  328.25710_gp,  304.38240_gp,  290.42420_gp,  283.81310_gp,&
         279.24290_gp,  264.02220_gp,  424.02640_gp,  376.03190_gp,  339.42410_gp,&
         311.50380_gp,  275.95340_gp,  448.59350_gp,  852.31430_gp,  661.91080_gp,&
         500.57750_gp,  505.46450_gp,  472.45030_gp,  525.94030_gp,  414.23030_gp,&
         387.85700_gp,  360.65710_gp,  348.65520_gp,  346.25000_gp,  539.16270_gp,&
         491.65990_gp,  463.03140_gp,  439.73150_gp,  403.50580_gp,  533.52020_gp,&
         1101.84900_gp,  836.68870_gp,  540.39720_gp,  529.27220_gp,  512.71470_gp,&
         514.64180_gp,  494.34850_gp,  517.63800_gp,  486.72230_gp,  494.23650_gp,&
         462.43040_gp,  449.42300_gp,  446.39320_gp,  468.34080_gp,  431.69210_gp,&
         577.00800_gp,  538.08640_gp,  493.41180_gp,  497.22060_gp,  440.14790_gp,&
         414.87760_gp,  396.75610_gp,  379.39980_gp,  367.97200_gp,  580.18030_gp/)
    vdwparams%coeffs(32301:32385)=(/  532.89470_gp,  518.91070_gp,  506.20180_gp,  474.70440_gp,  601.11980_gp,&
         1101.52770_gp,  639.23770_gp,  733.05960_gp,  657.21210_gp,  583.88650_gp,&
         562.56740_gp,  667.28070_gp,  157.19910_gp,  155.06850_gp,  114.88240_gp,&
         92.65080_gp,   64.17440_gp,  276.67550_gp,  333.98620_gp,  322.38950_gp,&
         281.99150_gp,  238.62230_gp,  388.20850_gp,  374.87440_gp,  379.65070_gp,&
         347.49740_gp,  264.05320_gp,  227.11690_gp,  222.82550_gp,  259.08200_gp,&
         241.94730_gp,  338.91650_gp,  358.57470_gp,  330.82030_gp,  309.11340_gp,&
         465.50130_gp,  450.81160_gp,  455.51320_gp,  440.42440_gp,  391.58340_gp,&
         347.19340_gp,  328.67310_gp,  333.31550_gp,  347.65890_gp,  438.17200_gp,&
         476.08450_gp,  450.91710_gp,  436.38380_gp,  564.01480_gp,  586.33110_gp,&
         437.25330_gp,  455.78380_gp,  425.55250_gp,  388.09680_gp,  366.29140_gp,&
         372.06660_gp,  391.08870_gp,  469.92250_gp,  496.66200_gp,  507.15130_gp,&
         504.36840_gp,  621.20410_gp,  119.07380_gp,  107.55120_gp,   81.78100_gp,&
         262.73830_gp,  313.86840_gp,  262.30480_gp,  231.91300_gp,  252.85690_gp,&
         279.61360_gp,  322.56300_gp,  319.09300_gp,  376.12340_gp,  430.71080_gp,&
         435.59950_gp,  421.08940_gp,  487.97400_gp,  492.91380_gp,  112.22540_gp,&
         88.70300_gp,  260.34680_gp,  286.11790_gp,  391.12500_gp,  455.28540_gp/)
    

    !this line can be used as a test for the c6ab function
    !redo the same loop to verify that the values can be entered differently
    idata=1
    il=0
    jl=0
      do iline=1,nlines
         if (jl==il) then
            jl=1
            il=il+1
         else
            jl=jl+1
         end if
         iat=vdwparams%ivalues(il)
         jat=vdwparams%ivalues(jl)
         !vdwparams%ivalues has to be used now
         call limit(iat,jat,iadr,jadr)
!!$         !check the associations for defining an alternative function
!!$         if (iat+(iadr-1)*100 /= vdwparams%ivalues(find_il(iat,iadr))) stop 'ival'
!!$         if (jat+(jadr-1)*100 /= vdwparams%ivalues(find_il(jat,jadr))) stop 'jval'
         
         vdwparams%maxcn(iat)=max(vdwparams%maxcn(iat),iadr)
         vdwparams%maxcn(jat)=max(vdwparams%maxcn(jat),jadr)
         
!!$         vdwparams%c6ab(iat,jat,iadr,jadr,1)=vdwparams%coeffs(iline)
!!$         if (iline /= find_iline(iat,jat,iadr,jadr) ) stop 'iline'
!!$         vdwparams%c6ab(iat,jat,iadr,jadr,2)=vdwparams%rvalues(il)
!!$         vdwparams%c6ab(iat,jat,iadr,jadr,3)=vdwparams%rvalues(jl)
!!$         vdwparams%c6ab(jat,iat,jadr,iadr,1)=vdwparams%coeffs(iline)
!!$         if (iline /= find_iline(jat,iat,jadr,iadr)) stop 'jline'
!!$         vdwparams%c6ab(jat,iat,jadr,iadr,2)=vdwparams%rvalues(jl)
!!$         vdwparams%c6ab(jat,iat,jadr,iadr,3)=vdwparams%rvalues(il)
!!$
!!$         if (vdwparams%coeffs(iline) /= vdwparams_c6ab(iat,jat,iadr,jadr,1)) stop 'icoeff'
!!$         if (vdwparams_c6ab(iat,jat,iadr,jadr,2)/=vdwparams%rvalues(il)) stop 'ival2'
!!$         if (vdwparams_c6ab(iat,jat,iadr,jadr,3)/=vdwparams%rvalues(jl)) stop 'ival3'
!!$         if (vdwparams_c6ab(jat,iat,jadr,iadr,1)/=vdwparams%coeffs(iline)) stop 'jcoeff'
!!$         if (vdwparams_c6ab(jat,iat,jadr,iadr,2)/=vdwparams%rvalues(jl)) stop 'jval2'
!!$         if (vdwparams_c6ab(jat,iat,jadr,iadr,3)/=vdwparams%rvalues(il)) stop 'jval3'
       idata=(iline*5)+1
    enddo
!!$stop
  END SUBROUTINE init_c6_params_d3

  pure function find_il(iat,iadr) result(il)
    implicit none
    integer, intent(in) :: iat,iadr
    integer :: il
    !local variables
    integer :: it,i

    !default value, useless
    il=1
    it=iat+(iadr-1)*100
    loop_find: do i=1,ntot_values
       if (it==vdwparams%ivalues(i)) then
          il=i
          exit loop_find
       end if
    end do loop_find
    
  end function find_il

  pure function find_iline(iat,jat,iadr,jadr) result(iline)
    implicit none
    integer, intent(in) :: iat,jat,iadr,jadr
    integer :: iline
    !local variables
    integer :: it,jt,imin,imax

    it=find_il(iat,iadr)
    jt=find_il(jat,jadr)

    !take the minimum and the max
    imin=min(it,jt)
    imax=max(it,jt)

    iline=(imax*(imax-1))/2+imin

  end function find_iline

  !> function of the parameters, substitute the array and save memory
  pure function vdwparams_c6ab(iat,jat,iadr,jadr,ind) result(c6ab)
    implicit none
    integer, intent(in) :: iat,jat,iadr,jadr,ind
    real(gp) :: c6ab

    select case(ind)
       case(1)
          c6ab=vdwparams%coeffs(find_iline(iat,jat,iadr,jadr))
       case(2)
          c6ab=vdwparams%rvalues(find_il(iat,iadr))
       case(3)
          c6ab=vdwparams%rvalues(find_il(jat,jadr))
       case default
          c6ab=0.0_gp
    end select
    
  end function vdwparams_c6ab

!
!     Small function to extract coordination number representation number
!     and atomic number from column two and three of the table
!     the convention for this function is that
!     iadr=iat/100+1
!     iat=mod(iat-1,100)+1
!     therefore the input value can be retrieved by iat=iat+(iadr-1)*100
  subroutine limit(iat,jat,iadr,jadr)
      implicit none
      integer             :: iat,jat
      integer, intent(out)            :: iadr,jadr
      iadr=1
      jadr=1
810   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 810
      endif
820   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 820
      endif
  END SUBROUTINE  limit

  subroutine init_qat
      implicit none
      integer            :: i
! Grimme:
! PBE0/def2-QZVP atomic values
      vdwparams%Qatom(1:MAX_ELEM)=(/&
     &  8.0589_GP,  3.4698_GP, 29.0974_GP, 14.8517_GP, 11.8799_GP,  7.8715_GP,  5.5588_GP,&
     &  4.7566_GP,  3.8025_GP,  3.1036_GP, 26.1552_GP, 17.2304_GP, 17.7210_GP, 12.7442_GP,&
     &  9.5361_GP,  8.1652_GP,  6.7463_GP,  5.6004_GP, 29.2012_GP, 22.3934_GP, 19.0598_GP,&
     & 16.8590_GP, 15.4023_GP, 12.5589_GP, 13.4788_GP, 12.2309_GP, 11.2809_GP, 10.5569_GP,&
     & 10.1428_GP,  9.4907_GP, 13.4606_GP, 10.8544_GP,  8.9386_GP,  8.1350_GP,  7.1251_GP,&
     &  6.1971_GP, 30.0162_GP, 24.4103_GP, 20.3537_GP, 17.4780_GP, 13.5528_GP, 11.8451_GP,&
     & 11.0355_GP, 10.1997_GP,  9.5414_GP,  9.0061_GP,  8.6417_GP,  8.9975_GP, 14.0834_GP,&
     & 11.8333_GP, 10.0179_GP,  9.3844_GP,  8.4110_GP,  7.5152_GP, 32.7622_GP, 27.5708_GP,&
     & 23.1671_GP, 21.6003_GP, 20.9615_GP, 20.4562_GP, 20.1010_GP, 19.7475_GP, 19.4828_GP,&
     & 15.6013_GP, 19.2362_GP, 17.4717_GP, 17.8321_GP, 17.4237_GP, 17.1954_GP, 17.1631_GP,&
     & 14.5716_GP, 15.8758_GP, 13.8989_GP, 12.4834_GP, 11.4421_GP, 10.2671_GP,  8.3549_GP,&
     &  7.8496_GP,  7.3278_GP,  7.4820_GP, 13.5124_GP, 11.6554_GP, 10.0959_GP,  9.7340_GP,&
     &  8.8584_GP,  8.0125_GP, 29.8135_GP, 26.3157_GP, 19.1885_GP, 15.8542_GP, 16.1305_GP,&
     & 15.6161_GP, 15.1226_GP, 16.1576_GP /)

      do i=1,MAX_ELEM
         vdwparams%Qatom(i)=sqrt(0.5_GP*vdwparams%Qatom(i)*real(i)**0.5_GP)
      end do
  END SUBROUTINE init_qat

  function c6cn(iat,jat,cni,cnj)
      implicit none
      integer, intent(in) ::iat,jat
      real(kind=GP), intent(in) ::cni,cnj
      real(kind=GP) :: c6cn
      integer       :: i,j
      real(kind=GP) :: top,bottom,dist,c6_ref
      real(kind=GP) :: cna,cnb
 
      c6cn=0.0_GP
      top=0.0_GP
      bottom=0.0_GP
      do i=1,vdwparams%maxcn(iat)
         do j=1,vdwparams%maxcn(jat)
            c6_ref=vdwparams_c6ab(iat,jat,i,j,1)
            if (c6_ref.gt.0.0_GP) then
               cna=vdwparams_c6ab(iat,jat,i,j,2)
               cnb=vdwparams_c6ab(iat,jat,i,j,3)
               dist=(cna-cni)**2+(cnb-cnj)**2
               top=top+exp(k3*dist)*c6_ref
               bottom=bottom+exp(k3*dist)
            endif
         enddo
      enddo
 
      if (bottom.gt.0.0_GP) then
         c6cn=top/bottom
      else
         c6cn=0.0_GP
      endif
  end function c6cn
  function crd_nr(iat,nat,xyz,atoms)
    
      use module_types
 
      implicit none
 
      type(atoms_data),                 intent(in) :: atoms
 
      real(kind=GP)                               :: crd_nr
      integer ,intent(in)                         :: iat,nat
      !!integer ,DIMENSION(nat),intent(in)          :: iz
      real(kind=GP), DIMENSION(3,nat),intent(in)  :: xyz
      real(kind=GP)                               :: dx,dy,dz,r,rcov

      integer                                     :: i
        
      crd_nr=0.0_GP
      do i=1,nat
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            !!rcov=vdwparams%cov_table(iz(i))+vdwparams%cov_table(iz(iat))
            rcov = vdwparams%cov_table(atoms%nzatom(atoms%astruct%iatype(i)))+&
             vdwparams%cov_table(atoms%nzatom(atoms%astruct%iatype(iat)))
            crd_nr=crd_nr+1.0_GP/(1.0_GP+exp(-k1*(rcov/r-1.0_GP)))
         endif
      enddo
  end function crd_nr

  subroutine crd_nr_der(n,xyz,cnij,cnijk,atoms)
 
      use module_types
      implicit none

      type(atoms_data),                 intent(in) :: atoms
 
      integer, intent(in)                         :: n
      real(kind=GP), DIMENSION(3,n),intent(in)    :: xyz
      real(kind=GP), DIMENSION(3,n),intent(out)   :: cnij
      real(kind=GP), DIMENSION(3,n,n),intent(out) :: cnijk
!!      integer, DIMENSION(n),intent(in)            :: iz
      integer                                     :: iat,i
      real(kind=GP)                               :: dx,dy,dz,r,cov_rad,expf,fac2,fac3

      do i=1,n
         cnij(1,i)=0.0_GP
         cnij(2,i)=0.0_GP
         cnij(3,i)=0.0_GP
         do iat=1,n
            if (i.ne.iat) then
               dx=xyz(1,iat)-xyz(1,i)
               dy=xyz(2,iat)-xyz(2,i)
               dz=xyz(3,iat)-xyz(3,i)
               r=sqrt(dx*dx+dy*dy+dz*dz)
               cov_rad=  vdwparams%cov_table(atoms%nzatom(atoms%astruct%iatype(iat)))&
                       + vdwparams%cov_table(atoms%nzatom(atoms%astruct%iatype(i)))
               expf=exp(-k1*((cov_rad/r)-1.0_GP))
               fac2=1.0_GP/(expf+1.0_GP)
               fac3=k1*cov_rad*expf*fac2*fac2/(r*r*r)
               cnij(1,i)=cnij(1,i)-fac3*dx
               cnij(2,i)=cnij(2,i)-fac3*dy
               cnij(3,i)=cnij(3,i)-fac3*dz
               cnijk(1,iat,i)=fac3*dx
               cnijk(2,iat,i)=fac3*dy
               cnijk(3,iat,i)=fac3*dz
            endif
         end do
      end do
  END SUBROUTINE crd_nr_der
  !!subroutine c6_grad(grad,iat,jat,kat,x,z,n,cnij,cnijk,atoms)
  subroutine c6_grad(grad,iat,jat,kat,x,n,cnij,cnijk,atoms)
      use module_types
      implicit none
 
      type(atoms_data),                 intent(in) :: atoms
      integer, intent(in)                          :: n,iat,jat,kat
      real(kind=GP), DIMENSION(3,n),intent(in)     :: x
!!      integer, DIMENSION(n),intent(in)             :: z
      real(kind=GP), DIMENSION(3),intent(out)      :: grad
      real(kind=GP), DIMENSION(3,n),intent(in)    :: cnij
      real(kind=GP), DIMENSION(3,n,n),intent(in)  :: cnijk

      real(kind=GP), DIMENSION(3)                  :: cnik
      real(kind=GP), DIMENSION(3)                  :: cnjk
      integer                                      :: i, j
      real(kind=GP)                                :: cni,cnj
      real(kind=GP)                                :: t1,t2,dt1x,dt1y,dt1z,dt2x,dt2y,dt2z
      real(kind=GP)                                :: tmp1,tmp2,tmp3,tmp4,fac1,fac2
 
      !!   cni=crd_nr(iat,n,x,z)
      cni  = crd_nr(iat,n,x,atoms)
      !!    cnj=crd_nr(jat,n,x,z)
      cnj  = crd_nr(jat,n,x,atoms)
      if (iat.eq.kat) then
         cnik(1)=cnij(1,iat)
         cnik(2)=cnij(2,iat)
         cnik(3)=cnij(3,iat)
         cnjk(1)=cnijk(1,iat,jat)
         cnjk(2)=cnijk(2,iat,jat)
         cnjk(3)=cnijk(3,iat,jat)
      else
         cnik(1)=cnijk(1,kat,iat)
         cnik(2)=cnijk(2,kat,iat)
         cnik(3)=cnijk(3,kat,iat)
         cnjk(1)=cnijk(1,kat,jat)
         cnjk(2)=cnijk(2,kat,jat)
         cnjk(3)=cnijk(3,kat,jat)
      endif
      t1=0.0_GP
      t2=0.0_GP
      dt1x=0.0_GP
      dt1y=0.0_GP
      dt1z=0.0_GP
      dt2x=0.0_GP
      dt2y=0.0_GP
      dt2z=0.0_GP
      do i=1,vdwparams%maxcn(atoms%nzatom(atoms%astruct%iatype(iat)))
        do j=1,vdwparams%maxcn(atoms%nzatom(atoms%astruct%iatype(jat)))
          tmp1=vdwparams_c6AB(atoms%nzatom(atoms%astruct%iatype(iat)),atoms%nzatom(atoms%astruct%iatype(jat)),i,j,3)-cnj
          tmp2=vdwparams_c6AB(atoms%nzatom(atoms%astruct%iatype(iat)),atoms%nzatom(atoms%astruct%iatype(jat)),i,j,2)-cni
          tmp3=exp(k3*(tmp1*tmp1+tmp2*tmp2))
          t1=t1+vdwparams_c6AB(atoms%nzatom(atoms%astruct%iatype(iat)),atoms%nzatom(atoms%astruct%iatype(jat)),i,j,1)*tmp3
          t2=t2+tmp3
          fac1=tmp3*k3*2.0_GP
          fac2=fac1*vdwparams_c6AB(atoms%nzatom(atoms%astruct%iatype(iat)),atoms%nzatom(atoms%astruct%iatype(jat)),i,j,1)
          tmp4=(tmp2*cnik(1)+tmp1*cnjk(1))
          dt1x=dt1x+fac2*tmp4
          dt2x=dt2x+fac1*tmp4
          tmp4=(tmp2*cnik(2)+tmp1*cnjk(2))
          dt1y=dt1y+fac2*tmp4
          dt2y=dt2y+fac1*tmp4
          tmp4=(tmp2*cnik(3)+tmp1*cnjk(3))
          dt1z=dt1z+fac2*tmp4
          dt2z=dt2z+fac1*tmp4
        enddo
      enddo
 
!     Numerical test to avoid division by zero (following Grimme)
 
      if (t2.gt.0.0_GP) then
        grad(1)=(dt1x*t2-dt2x*t1)/(t2**2)
        grad(2)=(dt1y*t2-dt2y*t1)/(t2**2)
        grad(3)=(dt1z*t2-dt2z*t1)/(t2**2)
      else
        grad(1)=0.0_GP
        grad(2)=0.0_GP
        grad(3)=0.0_GP
      endif
  END SUBROUTINE c6_grad




end module vdwcorrection

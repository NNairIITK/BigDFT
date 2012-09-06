!> @file
!!  Routines to do Van der Waals correction
!! @author
!! Written by Quintin Hill in 2007/8 with assistance from
!! Chris-Kriton Skylaris.\n
!! Forces added July 2008 by Quintin Hill\n
!! Modified for BigDFT in March/April 2009 by Quintin Hill\n
!! Modified for BigDFT in April 2012 by Alvaro (vama)\n
!! COPYRIGHT\n
!! Copyright (C)  2007-2009  Quintin Hill.
!! This file is distributed under the terms of the
!! GNU General Public License either version 2 of the License, or
!! (at your option) any later version, see ~/COPYING file,
!! http://www.gnu.org/licenses/gpl-2.0.txt (for GPL v2)
!! or http://www.gnu.org/copyleft/gpl.txt (for latest version).


!> @brief  Van der Waals empirical correction module. 
!! This module contains subroutines for calculating a Van der
!! Waals energy correction as a sum of damped London potentials.
module vdwcorrection

  use module_base, only: GP

  implicit none

  private

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
  END TYPE VDWPARAMETERS

  !> Van der Waals corrections.
  integer, parameter, public :: VDW_NONE = 0
  integer, parameter, public :: VDW_DAMP_ELSTNER = 1
  integer, parameter, public :: VDW_DAMP_WU_YANG_1 = 2
  integer, parameter, public :: VDW_DAMP_WU_YANG_2 = 3
  integer, parameter, public :: VDW_DAMP_GRIMME_D2 = 4
  character(len = 34), dimension(5), parameter :: vdw_correction_names = &
       & (/ "none                              ",   &
       &    "Damp from Elstner                 ",   &    
       &    "Damp from Wu & Yang               ",   &
       &    "Damp from Wu & Yang, second method",   &
       &    "Damp from Grimme D2               " /)

  public :: vdwcorrection_initializeparams
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
    real(kind=GP) :: jnm6mol2au
    real(kind=GP) :: angs2au
 
    angs2au=1.889725989_GP

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
       vdwparams%radzero(55)=315.275_GP
       vdwparams%radzero(56)=226.994_GP
       vdwparams%radzero(57)=176.252_GP
       vdwparams%radzero(58)=140.68_GP
       vdwparams%radzero(59)=140.68_GP
       vdwparams%radzero(60)=140.68_GP
       vdwparams%radzero(61)=140.68_GP
       vdwparams%radzero(62)=140.68_GP
       vdwparams%radzero(63)=140.68_GP
       vdwparams%radzero(64)=140.68_GP
       vdwparams%radzero(65)=140.68_GP
       vdwparams%radzero(66)=140.68_GP
       vdwparams%radzero(67)=140.68_GP
       vdwparams%radzero(68)=140.68_GP
       vdwparams%radzero(69)=140.68_GP
       vdwparams%radzero(70)=140.68_GP
       vdwparams%radzero(71)=140.68_GP
       vdwparams%radzero(72)=105.112_GP
       vdwparams%radzero(73)=81.24_GP
       vdwparams%radzero(74)=81.24_GP
       vdwparams%radzero(75)=81.24_GP
       vdwparams%radzero(76)=81.24_GP
       vdwparams%radzero(77)=81.24_GP
       vdwparams%radzero(78)=81.24_GP
       vdwparams%radzero(79)=81.24_GP
       vdwparams%radzero(80)=57.364_GP
       vdwparams%radzero(81)=57.254_GP
       vdwparams%radzero(82)=63.162_GP
       vdwparams%radzero(83)=63.540_GP
       vdwparams%radzero(84)=55.283_GP
       vdwparams%radzero(85)=57.171_GP
       vdwparams%radzero(86)=56.64_GP 

!!          jnm6mol2au=(((10*angs2au)**6)/4186)/627.51
       jnm6mol2au=17.336967408209503_GP
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
!! alpha
       vdwparams%dcoeff(3)=20.0000_GP
!! radii scale
       vdwparams%radscale=1.1000_GP

       select case (ixc)
   
       case(11)
!pbe
!!s6
          vdwparams%dcoeff(2)=0.7500_GP
       case(14) 
!revpbe
          vdwparams%dcoeff(2)=1.2500_GP
       case(-406000) 
! pbe0
          vdwparams%dcoeff(2)=0.6000_GP
       case(-170000) 
! b97-d
          vdwparams%dcoeff(2)=1.2500_GP
       case(-106132) 
! b-p
          vdwparams%dcoeff(2)=1.0500_GP
       case(-416000) 
! b-lyp
          vdwparams%dcoeff(2)=1.2000_GP
       case(-402000) 
! b3-lyp
          vdwparams%dcoeff(2)=1.0500_GP
       case(-202231) 
! tpss
          vdwparams%dcoeff(2)=1.0000_GP
       case default  
          vdwparams%dcoeff(2)=1.0000_GP
       end select

!! scale radio
       do ii=1, 86
          vdwparams%radzero(ii)=vdwparams%radzero(ii)*vdwparams%radscale*angs2au
       enddo
  endif

  END SUBROUTINE vdwcorrection_initializeparams


  !< This subroutine calculates the dispersion correction to the total energy.                                                    !
  !! @author
  !! Written by Quintin Hill in 2008, with some modifications in 2008
  !! Modified for BigDFT in March/April 2009 by Quintin Hill.
  subroutine vdwcorrection_calculate_energy(dispersion_energy,rxyz,atoms,dispersion,&
       iproc)

    use module_types

    implicit none

    ! Arguments
    type(atoms_data),                 intent(in) :: atoms
    real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
    integer,                          intent(in) :: iproc, dispersion
    real(kind=GP),                   intent(out) :: dispersion_energy

    ! Internal variables
    integer       :: atom1, atom2   ! atom counters for loops
    integer       :: nzatom1   ! Atomic number of atom 1
    integer       :: nzatom2   ! Atomic number of atom 2
    real(kind=GP) :: distance  ! Distance between pairs of atoms
    real(kind=GP) :: sqdist    ! Square distance between pairs of atoms
    real(kind=GP) :: c6coeff   ! The c6coefficient of the pair
    real(kind=GP) :: damping   ! The damping for the pair

    dispersion_energy = 0.0_GP

    if (dispersion /= VDW_NONE) then 


       ! qoh: Loop over all distinct pairs of atoms

       do atom1=1,atoms%nat
          do atom2=1,atom1-1

             nzatom1 = atoms%nzatom(atoms%iatype(atom1))
             nzatom2 = atoms%nzatom(atoms%iatype(atom2))

             ! qoh: Calculate c6 coefficient
             c6coeff = vdwcorrection_c6(nzatom1,nzatom2,dispersion)

             ! qoh: Calculate distance between each pair of atoms
             sqdist=(rxyz(1,atom2) - rxyz(1,atom1))**2 &
                  + (rxyz(2,atom2) - rxyz(2,atom1))**2 &
                  + (rxyz(3,atom2) - rxyz(3,atom1))**2 
             distance = sqrt(sqdist)

             ! qoh : Get damping function
             damping = vdwcorrection_damping(nzatom1,nzatom2,distance,dispersion)

             ! qoh: distance**6 = sqdist**3
             if (dispersion == 4) then
             dispersion_energy = dispersion_energy &
                  - (c6coeff * damping / sqdist**3)*vdwparams%dcoeff(2)
             else
             dispersion_energy = dispersion_energy &
                  - (c6coeff * damping / sqdist**3)
             endif

          enddo
       enddo

       if (iproc == 0) then
          write(*,'(1x,a, e12.5,1x,a)') &
               'Dispersion Correction Energy: ', dispersion_energy, 'Hartree'
       end if
    end if

  END SUBROUTINE vdwcorrection_calculate_energy


  !< This subroutine calculates the dispersion correction to the total energy.i
  !! @author
  !! Written by Quintin Hill in 2007, with some modifications in 2008
  !! Modified for BigDFT in March/April 2009 by Quintin Hill.
  subroutine vdwcorrection_calculate_forces(vdw_forces,rxyz,atoms,dispersion) 

    use module_types

    implicit none

    ! Arguments

    type(atoms_data),                 intent(in)  :: atoms
    real(GP), dimension(3,atoms%nat), intent(out) :: vdw_forces
    real(GP), dimension(3,atoms%nat), intent(in)  :: rxyz
    integer,                          intent(in)  :: dispersion

    ! Internal variables

    integer       :: atom1,atom2 ! atom counters for loops
    integer       :: nzatom1     ! Atomic number of atom 1
    integer       :: nzatom2     ! Atomic number of atom 2
    real(kind=GP) :: distance    ! Distance between pairs of atoms
    real(kind=GP) :: sqdist      ! Square distance between pairs of atoms
    real(kind=GP) :: c6coeff     ! The c6coefficient of the pair
    real(kind=GP) :: damping     ! The damping for the pair
    real(kind=GP) :: dampingdrv  ! The damping derivative for the pair
    real(kind=GP) :: drvcommon   ! The common part of the derivative

    vdw_forces = 0.0_GP

    if (dispersion /= VDW_NONE) then 

       ! qoh: Loop over all atoms

       do atom1=1,atoms%nat

          ! qoh: Loop over all other atoms

          do atom2=1,atoms%nat

             if ( atom2 .ne. atom1) then

                nzatom1 = atoms%nzatom(atoms%iatype(atom1))
                nzatom2 = atoms%nzatom(atoms%iatype(atom2))

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
                drvcommon = vdwparams%dcoeff(2)*(c6coeff/sqdist**3)*&
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
    end if


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
       do itype=1,atoms%ntypes
          if (any(unoptimised == atoms%nzatom(itype)) .and. &
               any(xcfoptimised == in%ixc)) then 
             write(*,'(a,a2)') 'WARNING: Unoptimised dispersion &
                  &parameters used for ', atoms%atomnames(itype)
          elseif (.not. any(optimised == atoms%nzatom(itype)) .and. &
               .not. any(unoptimised == atoms%nzatom(itype))) then
             write(*,'(a,a2)') 'WARNING: No dispersion parameters &
                  &available for ', atoms%atomnames(itype) 
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

    case(3,4)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_GP) 
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

    case(3,4)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_GP)

       vdwcorrection_drvdamping = vdwparams%dcoeff(3)*exp(expo)/&
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

end module vdwcorrection

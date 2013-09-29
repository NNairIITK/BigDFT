!> @file
!! Routine to flush a unit file 
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine to flush a unit file 
subroutine bigdft_utils_flush(unit)
   implicit none
   integer, intent(in) :: unit
   flush(unit=unit)
END SUBROUTINE bigdft_utils_flush

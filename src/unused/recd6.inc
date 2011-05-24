!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-6
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N=3
           real(kind=8), parameter :: W(-N:N) = (/&
           0.d0,0.0858797754503927846d0,1.04723768042233092d0,-0.188678293253531177d0,&
           0.079578122143014546d0,-0.028872131277603392d0,0.0048548465153963165d0&
           /)

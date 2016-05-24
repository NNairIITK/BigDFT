!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-8
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N= 4
           real(kind=8), parameter :: W(-N:N) = (/&
           0.d0,0.030093479272564175d0,0.8829879993841829d0,0.16246932689194357d0,&
          -0.1055654622867131d0,0.036244324645550385d0,-0.006464655184673674d0,&
           0.00013598646352976367d0,0.00009900081361574767d0&
           /)

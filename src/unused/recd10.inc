!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-10
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N=5
           real(kind=8), parameter :: W(-N:N) = (/&
           0.D0,0.008791347634813946D0,0.6384485113914279D0,0.6134449542982949D0,&
           -0.44963335759334067D0,0.31251560169295034D0,-0.1873500095585224D0,&
           0.08696481413462559D0,-0.028350887611356727D0,0.005696590191566821D0,&
           -0.0005275645804628203D0&
           /)

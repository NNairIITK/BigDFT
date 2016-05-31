!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-16
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N=8
           real(kind=8), parameter :: W(-N:N) = (/ 0.D0,&
          -0.00043552487975836995D0,0.13722387938570457D0,1.0097266554807494D0,&
          -0.281508018186976D0,0.4237806086614171D0,-0.7690029996917929D0,&
           1.0844256721218244D0,-1.1929529530682998D0,1.0356323940573589D0,&
          -0.7093626710811805D0,0.3796141707904326D0,-0.15558108556949435D0,&
           0.047191894571362525D0,-0.009987814863618105D0,0.0013173447537444538D0,&
          -0.00008155248178021157D0&
          /)

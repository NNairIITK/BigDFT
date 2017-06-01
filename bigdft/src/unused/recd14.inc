!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-14
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N=7
           real(kind=8), parameter :: W(-N:N) = (/0.d0,&
          -0.00019734150275056012d0,0.24695942557979256d0,1.067322833221752d0,&
          -0.6140032870771961d0,0.6889567761403583d0,-0.8383576343506007d0,&
           0.8619492953989152d0,-0.7042922242789302d0,0.44664455633125044d0,&
          -0.21485364926677875d0,0.07574869628227532d0,-0.018463488836077957d0,&
           0.002781179632359021d0,-0.00019513727462610458d0&
        /)

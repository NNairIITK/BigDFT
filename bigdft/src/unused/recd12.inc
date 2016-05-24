!> @file
!!          THE MAGIC FILTER FOR DAUBECHIES-12
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

           integer, parameter :: N=6
           real(kind=8), parameter :: W(-N:N) = (/0.d0,&
          0.0016835193068583942d0,0.4141977906939418d0,0.9438335333891432d0,&
         -0.6779372263414825d0,0.6222143057725837d0,-0.5531569377014378d0,&
          0.4070136773598989d0,-0.23061394243396194d0,0.09562202100986283d0,&
         -0.02722432135497499d0,0.004750520300849548d0,-0.0003829400012651252d0&
        /)

!> @file
!! @f$VALINTS(N_GAU)=INT(X^(2*N_GAU) EXP(-X^2) dx)@f$ 
!! from @f$-\infty@f$ TO @f$\infty@f$
!! The total squared norm of scfunction expansion coefficients should be equal to it
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

  real(gp), parameter :: valints(0:6) = (/& 
       1.7724538509055160272981674834_gp,0.8862269254527580136490837417_gp,&
       1.3293403881791370204736256126_gp,3.3233509704478425511840640314_gp,&
       11.63172839656744892914422411_gp,52.3427777845535201811490084949_gp,&
       287.8852778150443609963195467219_gp/)

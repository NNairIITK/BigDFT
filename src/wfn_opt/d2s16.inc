!> @file
!!  Second derivative filters for Daubechies 16
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

integer,parameter :: lowfil=-14,lupfil=14
real(kind=8) :: scale
real(kind=8) :: fil(lowfil:lupfil)

scale=-.5d0/hgrid**2

fil(0)=   -3.5536922899131901941296809374d0*scale
fil(1)=    2.2191465938911163898794546405d0*scale
fil(2)=   -0.6156141465570069496314853949d0*scale
fil(3)=    0.2371780582153805636239247476d0*scale
fil(4)=   -0.0822663999742123340987663521d0*scale
fil(5)=    0.02207029188482255523789911295638968409d0*scale
fil(6)=   -0.409765689342633823899327051188315485d-2*scale
fil(7)=    0.45167920287502235349480037639758496d-3*scale
fil(8)=   -0.2398228524507599670405555359023135d-4*scale
fil(9)=    2.0904234952920365957922889447361d-6*scale
fil(10)=  -3.7230763047369275848791496973044d-7*scale
fil(11)=  -1.05857055496741470373494132287d-8*scale
fil(12)=  -5.813879830282540547959250667d-11*scale
fil(13)=   2.70800493626319438269856689037647576d-13*scale
fil(14)=  -6.924474940639200152025730585882d-18*scale

do i=1,14
   fil(-i)=fil(i)
enddo

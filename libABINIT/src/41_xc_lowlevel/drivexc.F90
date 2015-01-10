!{\src2tex{textfont=tt}}
!!****f* ABINIT/drivexc
!! NAME
!! drivexc
!!
!! FUNCTION
!! Driver of XC functionals.
!! Treat spin-polarized as well as non-spin-polarized.
!! Treat local approximations or GGAs.
!! Optionally, deliver the XC kernel, or even the derivative
!! of the XC kernel (the third derivative of the XC energy)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ixc=number of the XC functional
!!   (to be described)
!!  ndvxc= size of dvxc(npts,ndvxc)
!!  ngr2= size of grho2_updn(npts,ngr2)
!!  nvxcdgr= size of vxcgr(npts,nvxcdgr)
!!  npts=number of real space points on which the density
!!   (and its gradients, if needed) is provided
!!  nspden=number of spin-density components (1 or 2)
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  rho_updn(npts,nspden)=the spin-up and spin-down densities
!!    If nspden=1, only the spin-up density must be given.
!!    In the calling routine, the spin-down density must
!!    be equal to the spin-up density,
!!    and both are half the total density.
!!    If nspden=2, the spin-up and spin-down density must be given
!!  Optional inputs :
!!  exexch= choice of local exact exchange. Active if exexch=3
!!  grho2_updn(npts,ngr2)=the square of the gradients of
!!    spin-up, spin-down, and total density
!!    If nspden=1, only the square of the gradient of the spin-up density
!!     must be given. In the calling routine, the square of the gradient
!!     of the spin-down density
!!     must be equal to the square of the gradient of the spin-up density,
!!     and both must be equal to one-quarter of the square of the
!!     gradient of the total density.
!!    If nspden=2, the square of the gradients of
!!     spin-up, spin-down, and total density must be given.
!!     Note that the square of the gradient of the total
!!     density is usually NOT related to the square of the
!!     gradient of the spin-up and spin-down densities,
!!     because the gradients are not usually aligned.
!!     This is not the case when nspden=1.
!!
!! OUTPUT
!!  exc(npts)=exchange-correlation energy density (hartree)
!!  vxc(npts,nspden)= (d($\rho$*exc)/d($\rho_up$)) (hartree)
!!               and  (d($\rho$*exc)/d($\rho_down$)) (hartree)
!!  vxcgr(npts,3)= 1/$|grad \rho_up|$ (d($\rho$*exc)/d($|grad \rho_up|$)) (hartree)
!!                 1/$|grad \rho_dn|$ (d($\rho$*exc)/d($|grad \rho_dn|$)) (hartree)
!!            and  1/$|grad \rho|$ (d($\rho$*exc)/d($|grad \rho|$))       (hartree)
!!     (will be zero if a LDA functional is used)
!!  Optional output :
!!  if(abs(order)>1)
!!   dvxc=partial second derivatives of the xc energy, only if abs(order)>1
!!   (This is a mess, to be rationalized !!)
!!   In case of local energy functional (option=1,-1 or 3):
!!    dvxc(npts,1+nspden)=              (Hartree*bohr^3)
!!     if(nspden=1 .and. order==2): dvxci(:,1)=dvxc/d$\rho$ , dvxc(:,2) empty
!!     if(nspden=1 .and. order==-2): also compute dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$
!!     if(nspden=2): dvxc(:,1)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxc(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   In case of gradient corrected functional (option=2,-2, 4, 5, 6, 7):
!!    dvxc(npts,15)=
!!     dvxc(:,1)= d2Ex/drho_up drho_up
!!     dvxc(:,2)= d2Ex/drho_dn drho_dn
!!     dvxc(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!     dvxc(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!     dvxc(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!     dvxc(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!     dvxc(:,7)= 1/abs(grad(rho_up)) * d/drho_up (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!     dvxc(:,8)= 1/abs(grad(rho_dn)) * d/drho_dn (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!     dvxc(:,9)= d2Ec/drho_up drho_up
!!     dvxc(:,10)=d2Ec/drho_up drho_dn
!!     dvxc(:,11)=d2Ec/drho_dn drho_dn
!!     dvxc(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!     dvxc(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!     dvxc(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!     dvxc(:,15)=1/abs(grad(rho)) * d/drho (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!
!!   if(abs(order)>2)  (only available for LDA and nspden=1)
!!    if nspden=1 d2vxc(npts,1)=second derivative of the XC potential=3rd order derivative of energy
!!    if nspden=2 d2vxc(npts,1), d2vxc(npts,2), d2vxc(npts,3), d2vxc(npts,4) (3rd derivative of energy)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,   & !Mandatory arguments
&                  dvxc,d2vxc,grho2_updn,vxcgr,exexch,lrho_updn,vxclrho,tau_updn,vxctau)    !Optional arguments

 use defs_basis
 use abi_interfaces_lowlevel
 use interfaces_41_xc_lowlevel, except_this_one => drivexc

#if defined HAVE_LIBXC
 use libxc_functionals
#endif

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,ndvxc,ngr2,nd2vxc,npts,nspden,nvxcdgr,order
 integer,intent(in),optional :: exexch
!arrays
 real(dp),intent(in) :: rho_updn(npts,nspden)
 real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
 real(dp),intent(in),optional :: lrho_updn(npts,nspden), tau_updn(npts,nspden)
 real(dp),intent(out) :: exc(npts),vxc(npts,nspden)
 real(dp),intent(out),optional :: d2vxc(npts,nd2vxc),dvxc(npts,ndvxc)
 real(dp),intent(out),optional :: vxcgr(npts,nvxcdgr)
 real(dp),intent(out),optional :: vxclrho(npts,nspden),vxctau(npts,nspden)

!Local variables-------------------------------
!scalars
 integer :: i_all,optpbe
 real(dp),parameter :: rsfac=0.6203504908994000e0_dp
 character(len=500) :: message
!arrays
 real(dp),allocatable :: exci_rpa(:)
 real(dp),allocatable :: rhotot(:),rspts(:),vxci_rpa(:,:),zeta(:)
!no_abirules
!real(dp),allocatable :: d2vxci(:,:),dvxci(:,:),grho2_updn_fake(:,:)

!  *************************************************************************

!Checks the values of order
 if( (order<1 .and. order/=-2) .or. order>4)then
   write(message, '(a,a,a,a,i6,a)' )ch10,&
&   ' drivexc : BUG -',ch10,&
&   '  The only allowed values for order are 1,2, -2, or 3, while it is found to be ',&
&   order,'.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

!Checks the compatibility between the inputs and the presence of the optional arguments

 if(present(dvxc))then
   if(order**2 <= 1 .or. ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27 )then
     write(message, '(8a,i6,a,i6)' )ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  The value of the number of the XC functional ixc',ch10,&
&     '  or the value of order is not compatible with the presence of the array dvxc',ch10,&
&     '  ixc=',ixc,'order=',order
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if
 end if

 if(present(d2vxc))then
   if(order /= 3 .or. (ixc /= 3 .and.  (((ixc > 15) .and. (ixc /=23)) .or. (ixc >= 0 .and. ixc < 7))) )then
     write(message, '(8a,i6,a,i6)' )ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  The value of the number of the XC functional ixc',ch10,&
&     '  or the value of order is not compatible with the presence of the array d2vxc',ch10,&
&     '  ixc=',ixc,'order=',order
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if
 end if

 if(present(vxcgr))then
   if(nvxcdgr == 0 .or. &
&   ((((ixc > 17 .and. ixc /= 23 .and. ixc/=26 .and. ixc/=27) .or. (ixc >= 0 .and. ixc < 7)) .and. nvxcdgr /=3 ))) then
     write(message, '(8a,i6,a,i6)' )ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  The value of the number of the XC functional ixc',ch10,&
&     '  or the value of nvxcdgr is not compatible with the presence of the array vxcgr',ch10,&
&     '  ixc=',ixc,'  nvxcdgr=',nvxcdgr
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if
 end if

 if(present(grho2_updn))then
   if (ngr2/=2*nspden-1 ) then
     write(message, '(4a)' ) ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  ngr2 must be 2*nspden-1 !'
!    call abi_wrtout(std_out,message,'COLL')
!    call abi_leave_new('COLL')
   end if
   if((ixc > 17 .and. ixc /= 23 .and. ixc/=26 .and. ixc/=27) .or. (ixc >= 0 .and. ixc < 11))then
     write(message, '(8a,i6)' )ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  The value of the number of the XC functional ixc',ch10,&
&     '  is not compatible with the presence of the array grho2_updn',ch10,&
&     '  ixc=',ixc
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if
 end if

!If needed, compute rhotot and rs
 if (ixc==1 .or. ixc==2 .or. ixc==3 .or. ixc==4 .or. ixc==5 .or. ixc==6 .or. &
& ixc==21 .or. ixc==22) then
   allocate(rhotot(npts),stat=i_all)
   allocate(rspts(npts),stat=i_all)
   if(nspden==1)then
     rhotot(:)=two*rho_updn(:,1)
   else
     rhotot(:)=rho_updn(:,1)+rho_updn(:,2)
   end if
   call invcb(rhotot,rspts,npts)
   rspts(:)=rsfac*rspts(:)
 end if
!If needed, compute zeta
 if (ixc==1 .or. ixc==21 .or. ixc==22) then
   allocate(zeta(npts),stat=i_all)
   if(nspden==1)then
     zeta(:)=zero
   else
     zeta(:)=two*rho_updn(:,1)/rhotot(:)-one
   end if
 end if

!Default value for vxcgr
 if (present(vxcgr)) vxcgr(:,:)=zero

!!$ !Could be more selective in allocating arrays ...
!!$ allocate(dvxci(npts,15),d2vxci(npts))

 if (ixc==1 .or. ixc==21 .or. ixc==22) then
!  new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option
   if (order**2 <= 1) then
     call xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc)
   else
     if(ndvxc /= nspden + 1 )then
       write(message, '(10a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  the value of the spin polarization',ch10,&
&       '  is not compatible with the ixc',ch10,&
&       '  ixc=',ixc,'nspden=',nspden
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,dvxc)
   end if

 else if (ixc==2) then
!  Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)
   if (order**2 <= 1) then
     call xcpzca(exc,npts,order,rhotot,rspts,vxc(:,1))
   else
     if(ndvxc /= 1 )then
       write(message, '(6a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xcpzca(exc,npts,order,rhotot,rspts,vxc(:,1),dvxc)
   end if

 else if (ixc==3) then
!  Teter fit (4/91) to Ceperley-Alder values (no spin-pol)
   if (order**2 <= 1) then
     call xctetr(exc,npts,order,rhotot,rspts,vxc(:,1))
   else if (order == 2) then
     if(ndvxc /= 1 )then
       write(message, '(6a,i3,a,i3)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xctetr(exc,npts,order,rhotot,rspts,vxc(:,1),dvxc=dvxc)
   else if (order == 3) then
     if(ndvxc /= 1 )then
       write(message, '(6a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xctetr(exc,npts,order,rhotot,rspts,vxc(:,1),d2vxc,dvxc)
   end if

 else if (ixc==4) then
!  Wigner xc (no spin-pol)
   if (order**2 <= 1) then
     call xcwign(exc,npts,order,rspts,vxc(:,1))
   else
     if(ndvxc /= 1 )then
       write(message, '(6a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xcwign(exc,npts,order,rspts,vxc(:,1),dvxc)
   end if

 else if (ixc==5) then
!  Hedin-Lundqvist xc (no spin-pol)
   if (order**2 <= 1) then
     call xchelu(exc,npts,order,rspts,vxc(:,1))
   else
     if(ndvxc /= 1 )then
       write(message, '(6a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xchelu(exc,npts,order,rspts,vxc(:,1),dvxc)
   end if

 else if (ixc==6) then
!  X-alpha (no spin-pol)
   if (order**2 <= 1) then
     call xcxalp(exc,npts,order,rspts,vxc(:,1))
   else
     if(ndvxc /= 1 )then
       write(message, '(6a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     call xcxalp(exc,npts,order,rspts,vxc(:,1),dvxc)
   end if

 else if (((ixc>=7 .and. ixc<=15) .or. (ixc==23)) .and. ixc/=10 .and. ixc/=13) then

!  Perdew-Wang LSD is coded in Perdew-Burke-Ernzerhof GGA, with optpbe=1
   if(ixc==7)optpbe=1
!  x-only part of Perdew-Wang
   if(ixc==8)optpbe=-1
!  Exchange + RPA correlation from Perdew-Wang
   if(ixc==9)optpbe=3
!  Perdew-Burke-Ernzerhof GGA
   if(ixc==11)optpbe=2
!  x-only part of PBE
   if(ixc==12)optpbe=-2
!  revPBE of Zhang and Yang
   if(ixc==14)optpbe=5
!  RPBE of Hammer, Hansen and Norskov
   if(ixc==15)optpbe=6
!  Wu and Cohen
   if(ixc==23)optpbe=7

   if (ixc >= 7 .and. ixc <= 9) then
     if (order**2 <= 1) then
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc)
     else if (order /=3) then
       if(ndvxc /= 1+nspden .or. nvxcdgr /=0)then
         write(message, '(6a,i6,a,i6,a,i6)' )ch10,&
&         ' drivexc : BUG -',ch10,&
&         '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&         '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
         call abi_wrtout(std_out,message,'COLL')
         call abi_leave_new('COLL')
       end if
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,dvxci=dvxc)
     else if (order ==3) then
       if(ndvxc /= 1+nspden .or. nvxcdgr /=0 .or. nd2vxc /=(3*nspden-2))then
         write(message, '(6a,i6,a,i6,a,i6,a,i6)' )ch10,&
&         ' drivexc : BUG -',ch10,&
&         '  Wrong value of ndvxc or nvxcdgr or nd2vxc:',ch10,&
&         '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr,'nd2vxc=',nd2vxc
         call abi_wrtout(std_out,message,'COLL')
         call abi_leave_new('COLL')
       end if
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,d2vxci=d2vxc,dvxci=dvxc)
     end if
   else if ((ixc >= 11 .and. ixc <= 15) .or. (ixc==23)) then
     if (order**2 <= 1) then
         if (present(exexch)) then
          call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,dvxcdgr=vxcgr,exexch=exexch,&
&       grho2_updn=grho2_updn)
         else
          call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,dvxcdgr=vxcgr,&
&       grho2_updn=grho2_updn)
         end if
     else if (order /=3) then
       if(ixc == 12 .and. ndvxc /=8  .or. ixc/=12 .and. ndvxc /= 15 .or. nvxcdgr /= 3)then
         write(message, '(6a,i6,a,i6,a,i6)' )ch10,&
&         ' drivexc : BUG -',ch10,&
&         '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&         '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
         call abi_wrtout(std_out,message,'COLL')
         call abi_leave_new('COLL')
       end if
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,&
       dvxcdgr=vxcgr,dvxci=dvxc,grho2_updn=grho2_updn)
     else if (order ==3) then
       if(ixc == 12 .and. ndvxc /=8  .or. ixc/=12 .and. ndvxc /= 15 .or. nvxcdgr /=3)then
         write(message, '(6a,i6,a,i6,a,i6)' )ch10,&
&         ' drivexc : BUG -',ch10,&
&         '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&         '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
         call abi_wrtout(std_out,message,'COLL')
         call abi_leave_new('COLL')
       end if
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,d2vxci=d2vxc,dvxcdgr=vxcgr,dvxci=dvxc,&
&       grho2_updn=grho2_updn)
     end if
   end if

 else if (ixc==10) then
!  RPA correlation from Perdew-Wang
   if (order**2 <= 1) then
     allocate(exci_rpa(npts),vxci_rpa(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,ngr2,nd2vxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc)
     exc(:)=exc(:)-exci_rpa(:)
!    PMA: second index of vxc is nspden while that of rpa is 2 they can mismatch
     vxc(:,1:min(nspden,2))=vxc(:,1:min(nspden,2))-vxci_rpa(:,1:min(nspden,2))
     deallocate(exci_rpa,vxci_rpa)
   else if (order /=3) then
     if(ndvxc /= 1+nspden .or. nvxcdgr /= 0)then
       write(message,'(6a,i6,a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     allocate(exci_rpa(npts),vxci_rpa(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,ngr2,nd2vxc,dvxci=dvxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,dvxci=dvxc)
     exc(:)=exc(:)-exci_rpa(:)
     vxc(:,:)=vxc(:,:)-vxci_rpa(:,:)
     deallocate(exci_rpa,vxci_rpa)
   else if (order ==3) then
     if(ndvxc /= 1+nspden .or. nvxcdgr /=0)then
       write(message,'(6a,i6,a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     allocate(exci_rpa(npts),vxci_rpa(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,ngr2,nd2vxc,&
&     d2vxci=d2vxc,dvxci=dvxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,&
&     d2vxci=d2vxc,dvxci=dvxc)
     exc(:)=exc(:)-exci_rpa(:)
     vxc(:,:)=vxc(:,:)-vxci_rpa(:,:)
     deallocate(exci_rpa,vxci_rpa)
   end if

 else if(ixc==13) then
!  LDA xc energy like ixc==7, and Leeuwen-Baerends GGA xc potential
   if (order**2 <= 1) then
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxc)
   else if (order /=3) then
     if(ndvxc /= 1+nspden .or. nvxcdgr /= 0)then
       write(message, '(6a,i6,a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,dvxci=dvxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxc)
   else if (order ==3) then
     if(ndvxc /= 1+nspden .or. nvxcdgr /=0)then
       write(message, '(6a,i6,a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  Wrong value of ndvxc or nvxcdgr:',ch10,&
&       '  ixc=',ixc,'ndvxc=',ndvxc,'nvxcdgr=',nvxcdgr
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,d2vxci=d2vxc,dvxci=dvxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxc)
   end if

 else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27) then
   if(nvxcdgr /= 2 )then
     write(message, '(6a,i6,a,i6)' )ch10,&
&     ' drivexc : BUG -',ch10,&
&     '  Wrong value of nvxcdgr:',ch10,&
&     '  ixc=',ixc,'ndvxcdgr=',nvxcdgr
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if
   call xchcth(vxcgr,exc,grho2_updn,ixc,npts,nspden,order,rho_updn,vxc)

 else if( ixc<0 ) then
   if(.false.)write(6,*)lrho_updn ! just to keep lrho_updn as an argument when LIBXC is not called
   if(.false.)write(6,*)tau_updn ! just to keep tau_updn as an argument when LIBXC is not called
#if defined HAVE_LIBXC

!  Check is all the necessary arrays are present and have the correct dimensions
   if (libxc_functionals_isgga() .or. libxc_functionals_ismgga()) then
     if ( (.not. present(grho2_updn)) .or. (.not. present(vxcgr)))  then
       write(message, '(8a,i7,a,i6,a,i6)' )ch10,&
&       ' drivexc : ERROR -',ch10,&
&       '  At least one of the functionals is a GGA or a MGGA,',ch10,&
&       '  but not all the necessary arrays are present.',ch10,&
&       '  ixc=',ixc,'  nvxcdgr=',nvxcdgr,'  ngr2=',ngr2
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if

     if (ngr2==0 .or. nvxcdgr/=3) then
       write(message, '(8a,i7,a,i6,a,i6)' )ch10,&
&       ' drivexc : BUG -',ch10,&
&       '  The value of the number of the XC functional ixc',ch10,&
&       '  is not compatible with the value of nvxcdgr or ngr2',ch10,&
&       '  ixc=',ixc,'  nvxcdgr=',nvxcdgr,'  ngr2=',ngr2
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
   end if

   if (libxc_functionals_ismgga()) then
     if ( (.not. present(lrho_updn)) .or. (.not. present(vxclrho)) .or. &
     (.not. present(tau_updn))  .or. (.not. present(vxctau))          )  then
       write(message, '(8a,i7)' )ch10,&
&       ' drivexc : ERROR -',ch10,&
&       '  At least one of the functionals is a MGGA,',ch10,&
&       '  but not all the necessary arrays are present.',ch10,&
&       '  ixc=',ixc
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
   end if

!  Call LibXC routines
   if (libxc_functionals_ismgga()) then
     call libxc_functionals_getvxc(npts,exc,nspden,rho_updn,vxc,grho2_updn,vxcgr,lrho_updn,vxclrho,tau_updn,vxctau)
   elseif (libxc_functionals_isgga()) then
     call libxc_functionals_getvxc(npts,exc,nspden,rho_updn,vxc,grho2_updn,vxcgr)
   else
     call libxc_functionals_getvxc(npts,exc,nspden,rho_updn,vxc)
   end if

#else
   write(message, '(5a)' )ch10,&
&   ' drivexc : ERROR -',ch10,&
&   '  ABINIT was not compiled with LibXC support',ch10
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
#endif
 end if

!Deallocate arrays
 if(allocated(rhotot))deallocate(rhotot)
 if(allocated(rspts))deallocate(rspts)
 if(allocated(zeta))deallocate(zeta)

end subroutine drivexc
!!***

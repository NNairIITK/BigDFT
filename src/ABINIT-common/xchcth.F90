!{\src2tex{textfont=tt}}
!!****f* ABINIT/xchcth
!! NAME
!! xchcth
!!
!! FUNCTION
!! Treat XC GGA functional of Hamprecht, Cohen, Tozer and Handy,
!! J. Chem. Phys. 109, 6264 (1998).
!!
!! For a series of values of the density and the square of the
!! gradient of the density, return the associated Exc energy,
!! potential, and, in case of response-function, functions needed
!! to build the XC kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2006 ABINIT group (XG,LG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,2*nspden-1)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2)
!!  order=maximal derivative of Exc to be computed
!!   (1 => energy and potential, or 2 => also XC kernel )
!!   Warning : order=2 not yet available
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!
!! OUTPUT
!!
!!  dvxcdgr(npts,3)=partial derivative of the exchange-correlation
!!    energy (exci*$\rho$) with respect to the spin-up (dvxcdgr(:,1)),
!!    spin-down (dvxcdgr(:,2)) square of gradients of the density
!!
!!  exci(npts)=exchange-correlation energy density (hartree)
!!  vxci(npts,nspden)=partial derivative of the exchange-correlation energy (exci*$\rho$)
!!    with respect to the spin-down (vxci(:,1)) and spin-up (vxci(:,2) densities
!! Normalization: Exc=$\int (exc(r)*\rho (r) d^3 r)$ for $\rho$(r)=electron density.
!!
!! TODO
!! Response function not coded yet, but part of it are already present
!!
!! NOTES
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      invcb,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xchcth(dvxcdgr,exci,grho2_updn,npts,nspden,&
& order,rho_updn,vxci)

 use defs_basis

!This section has been created automatically by the script Abilint (TD). Do not modify these by hand.
#ifdef HAVE_FORTRAN_INTERFACES
 use interfaces_01managempi
 use interfaces_03xc, except_this_one => xchcth
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts,nspden,order
!arrays
 real(dp),intent(in) :: grho2_updn(npts,2*nspden-1),rho_updn(npts,nspden)
 real(dp),intent(out) :: dvxcdgr(npts,2),exci(npts),vxci(npts,nspden)

!Local variables-------------------------------
!scalars
 integer,save :: initialized=0
 integer :: debug,ipts,ispden
 real(dp),parameter :: alpha_zeta2=1.0d0-1.0d-6,alpha_zeta=1.0d0-1.0d-6
 real(dp),parameter :: ccab0=0.729974d0,ccab1=0.335287d1,ccab2=-0.11543d2
 real(dp),parameter :: ccab3=0.808564d1,ccab4=-0.447857d1,ccsig0=0.222601d0
 real(dp),parameter :: ccsig1=-0.338622d-1,ccsig2=-0.125170d-1
 real(dp),parameter :: ccsig3=-0.802496d0,ccsig4=0.155396d1,cxsig0=0.109320d1
 real(dp),parameter :: cxsig1=-0.744056d0,cxsig2=0.559920d1,cxsig3=-0.678549d1
 real(dp),parameter :: cxsig4=0.449357d1,fsec_inv=1.0d0/1.709921d0
 real(dp),parameter :: gammacab=0.006d0,gammacsig=0.2d0,gammax=0.004d0
 real(dp),parameter :: rsfac=0.6203504908994000d0
 real(dp),save :: factf_zeta,factfp_zeta,sixpi2_1_3,sixpi2m1_3,sq_rsfac
 real(dp),save :: sq_rsfac_inv,threefourth_divpi,twom1_3
 real(dp) :: coeffss,d2ecrs0_drs2,d2ecrs1_drs2,d2ecrs_drs2,d2ecrs_drsdzeta
 real(dp) :: d2ecrs_dzeta2,d2fzeta4_dzeta2,d2gcrs_drs2,d2macrs_drs2,decrs0_drs
 real(dp) :: decrs1_drho,decrs1_drs,decrs_drs,decrs_dzeta,delta,dfzeta4_dzeta
 real(dp) :: dgcabdss,dgcrs_drs,dgcsigdss,dgxsigdss,divcab,divcsig,divx
 real(dp) :: dmacrs_drs,drhoecab_drhodn,drhoecab_drhoup,drhoecrs1_drhodn
 real(dp) :: drhoecrs1_drhoup,drsdrho,dssdg,dssdndg,dssdndrho,dssdrho,dssupdg
 real(dp) :: dssupdrho,ducabdss,ducsigdss,duxsigdss,ec0_a1,ec0_aa,ec0_b1,ec0_b2
 real(dp) :: ec0_b3,ec0_b4,ec0_den,ec0_log,ec0_q0,ec0_q1,ec0_q1p,ec0_q1pp
 real(dp) :: ec1_a1,ec1_aa,ec1_b1,ec1_b2,ec1_b3,ec1_b4,ec1_den,ec1_log,ec1_q0
 real(dp) :: ec1_q1,ec1_q1p,ec1_q1pp,ecrs,ecrs0,ecrs1,ex_lsd,exc,f_zeta
 real(dp) :: factfpp_zeta,factor,fp_zeta,fpp_zeta,gcab,gcrs,gcsig,grr,gxsig
 real(dp) :: mac_a1,mac_aa,mac_b1,mac_b2,mac_b3,mac_b4,mac_den,mac_log,mac_q0
 real(dp) :: mac_q1,mac_q1p,mac_q1pp,macrs,rho,rho_dn,rho_dnm,rho_dnp,rho_inv
 real(dp) :: rho_up,rho_upm,rho_upp,rhoecab,rhoecrs1_dn,rhoecrs1_up,rhomo6
 real(dp) :: rhomot,rhoo6,rhotmo6,rhotmot,rhoto6,rhotot,rhotot_inv,rs,rsm1_2
 real(dp) :: sqr_rs,ss,ss_dn,ss_up,ssavg,ucab,ucsig,uxsig,vxcadd,zeta,zeta4
 real(dp) :: zeta_mean,zetm_1_3,zetp_1_3
 character(len=500) :: message
!arrays
 real(dp),allocatable :: rho_updnm1_3(:,:),rhoarr(:),rhom1_3(:),zetm(:)
 real(dp),allocatable :: zetmm1_3(:),zetp(:),zetpm1_3(:)

! *************************************************************************

!DEBUG
!write(6,*)' xchcth : enter'
!write(6,*)' nspden=',nspden
!ENDDEBUG

 if (order/=1) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' xchcth : BUG -',ch10,&
&   '  Order must be 1 ; argument was ',order,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if(initialized==0)then
  twom1_3=two**(-third)
  sixpi2_1_3=(six*pi**2)**third
  sixpi2m1_3=one/sixpi2_1_3
  threefourth_divpi=three_quarters*piinv
  factf_zeta= one / ( two**(four/three)-two )
  factfp_zeta= four_thirds * factf_zeta * alpha_zeta2
  sq_rsfac=sqrt(rsfac)
  sq_rsfac_inv=one/sq_rsfac
  initialized=1
 end if

!Parameters for the Perdew-Wang 92 LSD as well as LSD-RPA,
!see Table I of Phys.Rev.B 45,13244 (1992)
 ec0_aa=0.031091d0 ; ec1_aa=0.015545d0 ; mac_aa=0.016887d0
 ec0_a1=0.21370d0  ; ec1_a1=0.20548d0  ; mac_a1=0.11125d0
 ec0_b1=7.5957d0   ; ec1_b1=14.1189d0  ; mac_b1=10.357d0
 ec0_b2=3.5876d0   ; ec1_b2=6.1977d0   ; mac_b2=3.6231d0
 ec0_b3=1.6382d0   ; ec1_b3=3.3662d0   ; mac_b3=0.88026d0
 ec0_b4=0.49294d0  ; ec1_b4=0.62517d0  ; mac_b4=0.49671d0

!DEBUG
!Finite-difference debugging, do not take away
!Note : here work with collinear gradients. Might be generalized ...
!debug=2  ! Choose 1 (rho grads) or 2 (grho grads)
!factor=1.0d0
!zeta_mean=0.98d0
!zeta_mean=zero
!delta=0.000025*factor
!delta=0.0000125*factor
!if(debug/=0)then
! do ipts=1,npts,5
!  rho=ipts*0.01d0*factor
!  rho_up=rho*(1.0d0+zeta_mean)*0.5d0
!  rho_dn=rho*(1.0d0-zeta_mean)*0.5d0
!  rho_upp=rho_up+delta
!  rho_upm=rho_up-delta
!  rho_dnp=rho_dn+delta
!  rho_dnm=rho_dn-delta
!! Here, vary rho
!  if(debug==1)then
!   rho_updn(ipts  ,1)=rho_up ; rho_updn(ipts  ,2)=rho_dn
!   rho_updn(ipts+1,1)=rho_upp; rho_updn(ipts+1,2)=rho_dn
!   rho_updn(ipts+2,1)=rho_upm; rho_updn(ipts+2,2)=rho_dn
!   rho_updn(ipts+3,1)=rho_up ; rho_updn(ipts+3,2)=rho_dnp
!   rho_updn(ipts+4,1)=rho_up ; rho_updn(ipts+4,2)=rho_dnm
!   grho2_updn(ipts:ipts+4,1)=(0.2d0*factor)**2     ! grad2 of spin up density
!   grho2_updn(ipts:ipts+4,2)=(0.2d0*factor)**2     ! grad2 of spin down density
!   grho2_updn(ipts:ipts+4,3)=(0.3d0*factor)**2     ! grad2 of total density
!  else
!!  Here, vary grho (interchange rho and grho)
!   grho2_updn(ipts  ,1)=rho_up**2 ; grho2_updn(ipts  ,2)=rho_dn**2
!   grho2_updn(ipts+1,1)=rho_upp**2; grho2_updn(ipts+1,2)=rho_dn**2
!   grho2_updn(ipts+2,1)=rho_upm**2; grho2_updn(ipts+2,2)=rho_dn**2
!   grho2_updn(ipts+3,1)=rho_up**2 ; grho2_updn(ipts+3,2)=rho_dnp**2
!   grho2_updn(ipts+4,1)=rho_up**2 ; grho2_updn(ipts+4,2)=rho_dnm**2
!   grho2_updn(ipts  ,3)=(ipts*0.01d0*factor)**2
!   grho2_updn(ipts+1,3)=(ipts*0.01d0*factor+delta)**2
!   grho2_updn(ipts+2,3)=(ipts*0.01d0*factor-delta)**2
!   grho2_updn(ipts+3,3)=(ipts*0.01d0*factor+delta)**2   ! identical to ipts+1
!   grho2_updn(ipts+4,3)=(ipts*0.01d0*factor-delta)**2   ! identical to ipts+2
!   rho_updn(ipts:ipts+4,1)=0.2d0*factor*(1.0d0+zeta_mean)*0.5d0    ! spin up density
!   rho_updn(ipts:ipts+4,2)=0.2d0*factor*(1.0d0-zeta_mean)*0.5d0    ! spin down density
!  end if
! end do
!end if
!Usual option :
!nspden=2 ; order=2
!GGA
!nspden=2 ; order=1
!Might take also, although finite difference later is meaningless
!nspden=1 ; order=-2
!ENDDEBUG

 if(order**2 >1)then
  factfpp_zeta= third * factfp_zeta * alpha_zeta2
 end if

 allocate(rhoarr(npts),rhom1_3(npts),rho_updnm1_3(npts,2))
 do ispden=1,nspden
  call invcb(rho_updn(:,ispden),rho_updnm1_3(:,ispden),npts)
 end do
 if(nspden==1)then
  rhoarr(:)=two*rho_updn(:,1)
  rhom1_3(:)=twom1_3*rho_updnm1_3(:,1)
  rho_updnm1_3(:,2)=rho_updnm1_3(:,1)
 else
  rhoarr(:)=rho_updn(:,1)+rho_updn(:,2)
  call invcb(rhoarr,rhom1_3,npts)
  allocate(zetm(npts),zetmm1_3(npts),zetp(npts),zetpm1_3(npts))
  do ipts=1,npts
   rhotmot=rhom1_3(ipts)
   rhotot_inv=rhotmot*rhotmot*rhotmot
   zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
   zetp(ipts)=1.0d0+zeta*alpha_zeta
   zetm(ipts)=1.0d0-zeta*alpha_zeta
  end do
  call invcb(zetp,zetpm1_3,npts)
  call invcb(zetm,zetmm1_3,npts)
 end if
 
 if (nspden==1) then
    
    if(order==-2) then

       do ipts=1,npts

          ! -----------------------------------------------------------------------
          ! First take care of the spin-split part of the functional
          exc=zero
          ispden=1
          rho   =rho_updn(ipts,ispden)
          rhomot=rho_updnm1_3(ipts,ispden)
          rho_inv=rhomot*rhomot*rhomot

          !  Exchange part
          ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
          !  Note that this definition differs from the PBE one
          coeffss=rho_inv*rho_inv*rhomot*rhomot
          ss=grho2_updn(ipts,ispden)*coeffss
          dssdrho=-eight*third*ss*rho_inv
          dssdg=two*coeffss

          divx=one/(one+gammax*ss)
          uxsig=gammax*ss*divx
          duxsigdss=gammax*divx*(one-ss*gammax*divx)

          gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
          dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
               &   *duxsigdss

          exc=exc+ex_lsd*rho*gxsig
          vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
          dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

          !  Spin parallel correlation part
          !  Note that this definition is for fully spin-polarized quantities
          rs=rsfac*rhomot
          rhomo6=sqrt(rhomot)
          sqr_rs=sq_rsfac*rhomo6
          rhoo6=rho*rhomot*rhomot*rhomo6
          rsm1_2=sq_rsfac_inv*rhoo6
          drsdrho=-third*rs*rho_inv

          ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
          ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
          ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
          ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
          !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
          ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
          ecrs1=ec1_q0*ec1_log
          decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
          decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

          !  Store the LSDA corr energy and density derivative
          rhoecrs1_up=rho*ecrs1
          drhoecrs1_drhoup=decrs1_drho
          ss_up=ss
          dssupdrho=dssdrho
          dssupdg=dssdg

          divcsig=one/(one+gammacsig*ss)
          ucsig=gammacsig*ss*divcsig
          ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

          gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
          dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
               &   *ducsigdss

          !DEBUG
          !  gcsig=zero
          !  dgcsigdss=zero
          !ENDDEBUG
          exc=exc+ecrs1*rho*gcsig
          vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
               &                    ecrs1*rho*dgcsigdss*dssdrho
          dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

          rhoecrs1_dn=rhoecrs1_up
          drhoecrs1_drhodn=drhoecrs1_drhoup
          ss_dn=ss_up
          dssdndrho=dssupdrho
          dssdndg=dssupdg
          exc=exc*two

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

          ! From now, the coding of the PW92 functional start. It is identical in xcpbe.f
          rs=rsfac*rhotmot
          sqr_rs=sq_rsfac*rhotmo6
          rsm1_2=sq_rsfac_inv*rhoto6

          ! Formulas A6-A8 of PW92LSD
          ec0_q0=-2.0d0*ec0_aa*(1.0d0+ec0_a1*rs)
          ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
          ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
          ec0_den=1.0d0/(ec0_q1*ec0_q1+ec0_q1)
          ! ec0_log=log( 1.0d0 + 1.0d0 / ec0_q1 )
          ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
          ecrs0=ec0_q0*ec0_log
          decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

          ec0_q1pp=0.5d0*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
          d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
               &               -ec0_q0*ec0_q1pp*ec0_den                        &
               &               +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+1.0d0)



          mac_q0=-2.0d0*mac_aa*(1.0d0+mac_a1*rs)
          mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
          mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
          mac_den=1.0d0/(mac_q1*mac_q1+mac_q1)
          mac_log=-log( mac_q1*mac_q1*mac_den )
          macrs=mac_q0*mac_log
          dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den


          ecrs=ecrs0
          decrs_drs=decrs0_drs
          decrs_dzeta=0.0d0

          d2ecrs_drs2=d2ecrs0_drs2

          d2ecrs_dzeta2=alpha_zeta**2*(-macrs)

          zeta=0.0d0

          ! At this point, the coding of the PW92 functional finishes.

          ! The correlation between different spin from HCTH is now computed
          ! First, the part without gradient correction factor
          rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
          vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
          drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
          drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

          ! Now, the gradient correction factor
          ssavg=half*(ss_up+ss_dn)
          divcab=one/(one+gammacab*ssavg)
          ucab=gammacab*ssavg*divcab
          ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

          gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
          dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
               &   *ducabdss

          exc=exc+rhoecab*gcab

          vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
          dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
          !  If non spin-polarized, treat spin down contribution now, similar to spin up
          !  vxci(ipts,2)=vxci(ipts,1)
          dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          ! Final division by the total density, to give the energy density
          exci(ipts)=exc*rhotot_inv

       end do ! ipts=1,npts

    else if(order**2>1) then

       do ipts=1,npts

          ! -----------------------------------------------------------------------
          ! First take care of the spin-split part of the functional
          exc=zero
          ispden=1
          rho   =rho_updn(ipts,ispden)
          rhomot=rho_updnm1_3(ipts,ispden)
          rho_inv=rhomot*rhomot*rhomot

          !  Exchange part
          ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
          !  Note that this definition differs from the PBE one
          coeffss=rho_inv*rho_inv*rhomot*rhomot
          ss=grho2_updn(ipts,ispden)*coeffss
          dssdrho=-eight*third*ss*rho_inv
          dssdg=two*coeffss

          divx=one/(one+gammax*ss)
          uxsig=gammax*ss*divx
          duxsigdss=gammax*divx*(one-ss*gammax*divx)

          gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
          dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
               &   *duxsigdss

          exc=exc+ex_lsd*rho*gxsig
          vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
          dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

          !  Spin parallel correlation part
          !  Note that this definition is for fully spin-polarized quantities
          rs=rsfac*rhomot
          rhomo6=sqrt(rhomot)
          sqr_rs=sq_rsfac*rhomo6
          rhoo6=rho*rhomot*rhomot*rhomo6
          rsm1_2=sq_rsfac_inv*rhoo6
          drsdrho=-third*rs*rho_inv

          ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
          ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
          ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
          ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
          !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
          ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
          ecrs1=ec1_q0*ec1_log
          decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
          decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

          !  Store the LSDA corr energy and density derivative
          rhoecrs1_up=rho*ecrs1
          drhoecrs1_drhoup=decrs1_drho
          ss_up=ss
          dssupdrho=dssdrho
          dssupdg=dssdg

          divcsig=one/(one+gammacsig*ss)
          ucsig=gammacsig*ss*divcsig
          ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

          gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
          dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
               &   *ducsigdss

          !DEBUG
          !  gcsig=zero
          !  dgcsigdss=zero
          !ENDDEBUG
          exc=exc+ecrs1*rho*gcsig
          vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
               &                    ecrs1*rho*dgcsigdss*dssdrho
          dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

          rhoecrs1_dn=rhoecrs1_up
          drhoecrs1_drhodn=drhoecrs1_drhoup
          ss_dn=ss_up
          dssdndrho=dssupdrho
          dssdndg=dssupdg
          exc=exc*two

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

          ! From now, the coding of the PW92 functional start. It is identical in xcpbe.f
          rs=rsfac*rhotmot
          sqr_rs=sq_rsfac*rhotmo6
          rsm1_2=sq_rsfac_inv*rhoto6
          
          ! Formulas A6-A8 of PW92LSD
          ec0_q0=-2.0d0*ec0_aa*(1.0d0+ec0_a1*rs)
          ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
          ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
          ec0_den=1.0d0/(ec0_q1*ec0_q1+ec0_q1)
          ! ec0_log=log( 1.0d0 + 1.0d0 / ec0_q1 )
          ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
          ecrs0=ec0_q0*ec0_log
          decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
          ec0_q1pp=0.5d0*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
          d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
               &               -ec0_q0*ec0_q1pp*ec0_den                        &
               &               +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+1.0d0)


          ecrs=ecrs0
          decrs_drs=decrs0_drs
          decrs_dzeta=0.0d0
          d2ecrs_drs2=d2ecrs0_drs2
          zeta=0.0d0

          ! At this point, the coding of the PW92 functional finishes.

          ! The correlation between different spin from HCTH is now computed
          ! First, the part without gradient correction factor
          rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
          vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
          drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
          drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

          ! Now, the gradient correction factor
          ssavg=half*(ss_up+ss_dn)
          divcab=one/(one+gammacab*ssavg)
          ucab=gammacab*ssavg*divcab
          ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

          gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
          dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
               &   *ducabdss

          exc=exc+rhoecab*gcab

          vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
          dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
          !  If non spin-polarized, treat spin down contribution now, similar to spin up
          !  vxci(ipts,2)=vxci(ipts,1)
          dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          ! Final division by the total density, to give the energy density
          exci(ipts)=exc*rhotot_inv

       end do ! ipts=1,npts
    else

       do ipts=1,npts

          ! -----------------------------------------------------------------------
          ! First take care of the spin-split part of the functional
          exc=zero
          ispden=1
          rho   =rho_updn(ipts,ispden)
          rhomot=rho_updnm1_3(ipts,ispden)
          rho_inv=rhomot*rhomot*rhomot

          !  Exchange part
          ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
          !  Note that this definition differs from the PBE one
          coeffss=rho_inv*rho_inv*rhomot*rhomot
          ss=grho2_updn(ipts,ispden)*coeffss
          dssdrho=-eight*third*ss*rho_inv
          dssdg=two*coeffss

          divx=one/(one+gammax*ss)
          uxsig=gammax*ss*divx
          duxsigdss=gammax*divx*(one-ss*gammax*divx)

          gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
          dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
               &   *duxsigdss

          exc=exc+ex_lsd*rho*gxsig
          vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
          dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

          !  Spin parallel correlation part
          !  Note that this definition is for fully spin-polarized quantities
          rs=rsfac*rhomot
          rhomo6=sqrt(rhomot)
          sqr_rs=sq_rsfac*rhomo6
          rhoo6=rho*rhomot*rhomot*rhomo6
          rsm1_2=sq_rsfac_inv*rhoo6
          drsdrho=-third*rs*rho_inv

          ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
          ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
          ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
          ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
          !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
          ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
          ecrs1=ec1_q0*ec1_log
          decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
          decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

          !  Store the LSDA corr energy and density derivative
          rhoecrs1_up=rho*ecrs1
          drhoecrs1_drhoup=decrs1_drho
          ss_up=ss
          dssupdrho=dssdrho
          dssupdg=dssdg

          divcsig=one/(one+gammacsig*ss)
          ucsig=gammacsig*ss*divcsig
          ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

          gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
          dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
               &   *ducsigdss

          !DEBUG
          !  gcsig=zero
          !  dgcsigdss=zero
          !ENDDEBUG
          exc=exc+ecrs1*rho*gcsig
          vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
               &                    ecrs1*rho*dgcsigdss*dssdrho
          dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

          rhoecrs1_dn=rhoecrs1_up
          drhoecrs1_drhodn=drhoecrs1_drhoup
          ss_dn=ss_up
          dssdndrho=dssupdrho
          dssdndg=dssupdg
          exc=exc*two

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

          ! From now, the coding of the PW92 functional start. It is identical in xcpbe.f
          rs=rsfac*rhotmot
          sqr_rs=sq_rsfac*rhotmo6
          rsm1_2=sq_rsfac_inv*rhoto6

          ! Formulas A6-A8 of PW92LSD
          ec0_q0=-2.0d0*ec0_aa*(1.0d0+ec0_a1*rs)
          ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
          ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
          ec0_den=1.0d0/(ec0_q1*ec0_q1+ec0_q1)
          ! ec0_log=log( 1.0d0 + 1.0d0 / ec0_q1 )
          ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
          ecrs0=ec0_q0*ec0_log
          decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

          ecrs=ecrs0
          decrs_drs=decrs0_drs
          decrs_dzeta=0.0d0
          zeta=0.0d0

          ! At this point, the coding of the PW92 functional finishes.

          ! The correlation between different spin from HCTH is now computed
          ! First, the part without gradient correction factor
          rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
          vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
          drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
          drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

          ! Now, the gradient correction factor
          ssavg=half*(ss_up+ss_dn)
          divcab=one/(one+gammacab*ssavg)
          ucab=gammacab*ssavg*divcab
          ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

          gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
          dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
               &   *ducabdss

          exc=exc+rhoecab*gcab

          vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
          dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
          !  If non spin-polarized, treat spin down contribution now, similar to spin up
          !  vxci(ipts,2)=vxci(ipts,1)
          dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

          ! Final division by the total density, to give the energy density
          exci(ipts)=exc*rhotot_inv

       end do ! ipts=1,npts
    end if

 else if(nspden==2) then
    
    if(order**2>1) then

       do ipts=1,npts

          ! -----------------------------------------------------------------------
          ! First take care of the spin-split part of the functional
          exc=zero
          do ispden=1,nspden
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             rho_inv=rhomot*rhomot*rhomot

             !  Exchange part
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !  Note that this definition differs from the PBE one
             coeffss=rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             dssdrho=-eight*third*ss*rho_inv
             dssdg=two*coeffss

             divx=one/(one+gammax*ss)
             uxsig=gammax*ss*divx
             duxsigdss=gammax*divx*(one-ss*gammax*divx)

             gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
             dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
                  &   *duxsigdss

             exc=exc+ex_lsd*rho*gxsig
             vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

             !  Spin parallel correlation part
             !  Note that this definition is for fully spin-polarized quantities
             rs=rsfac*rhomot
             rhomo6=sqrt(rhomot)
             sqr_rs=sq_rsfac*rhomo6
             rhoo6=rho*rhomot*rhomot*rhomo6
             rsm1_2=sq_rsfac_inv*rhoo6
             drsdrho=-third*rs*rho_inv

             ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
             ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
             ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
             ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
             !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
             ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
             ecrs1=ec1_q0*ec1_log
             decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
             decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

             !  Store the LSDA corr energy and density derivative
             if(ispden==1)then
                rhoecrs1_up=rho*ecrs1
                drhoecrs1_drhoup=decrs1_drho
                ss_up=ss
                dssupdrho=dssdrho
                dssupdg=dssdg
             else
                rhoecrs1_dn=rho*ecrs1
                drhoecrs1_drhodn=decrs1_drho
                ss_dn=ss
                dssdndrho=dssdrho
                dssdndg=dssdg
             end if

             divcsig=one/(one+gammacsig*ss)
             ucsig=gammacsig*ss*divcsig
             ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

             gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
             dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
                  &   *ducsigdss

             !DEBUG
             !  gcsig=zero
             !  dgcsigdss=zero
             !ENDDEBUG
             exc=exc+ecrs1*rho*gcsig
             vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
                  &                    ecrs1*rho*dgcsigdss*dssdrho
             dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

          end do

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

          ! From now, the coding of the PW92 functional start. It is identical in xcpbe.f
          rs=rsfac*rhotmot
          sqr_rs=sq_rsfac*rhotmo6
          rsm1_2=sq_rsfac_inv*rhoto6

          ! Formulas A6-A8 of PW92LSD
          ec0_q0=-2.0d0*ec0_aa*(1.0d0+ec0_a1*rs)
          ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
          ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
          ec0_den=1.0d0/(ec0_q1*ec0_q1+ec0_q1)
          ! ec0_log=log( 1.0d0 + 1.0d0 / ec0_q1 )
          ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
          ecrs0=ec0_q0*ec0_log
          decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
          ec0_q1pp=0.5d0*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
          d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
               &               -ec0_q0*ec0_q1pp*ec0_den                        &
               &               +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+1.0d0)

          mac_q0=-2.0d0*mac_aa*(1.0d0+mac_a1*rs)
          mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
          mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
          mac_den=1.0d0/(mac_q1*mac_q1+mac_q1)
          mac_log=-log( mac_q1*mac_q1*mac_den )
          macrs=mac_q0*mac_log
          dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

          zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
          ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
          ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
          ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
          ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
          !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
          ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
          ecrs1=ec1_q0*ec1_log
          decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

          !  alpha_zeta is introduced in order to remove singularities for fully
          !  polarized systems.
          zetp_1_3=(1.0d0+zeta*alpha_zeta)*zetpm1_3(ipts)**2
          zetm_1_3=(1.0d0-zeta*alpha_zeta)*zetmm1_3(ipts)**2

          f_zeta=( (1.0d0+zeta*alpha_zeta2)*zetp_1_3 +                      &
               &           (1.0d0-zeta*alpha_zeta2)*zetm_1_3 - 2.0d0 ) * factf_zeta
          fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
          zeta4=zeta**4
          gcrs=ecrs1-ecrs0+macrs*fsec_inv
          !  ecrs=ecrs0+f_zeta*(-macrs*(1.0d0-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
          ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

          dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
          !  decrs_drs=decrs0_drs+f_zeta*&
          !&       (-dmacrs_drs*(1.0d0-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
          decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
          dfzeta4_dzeta=4.0d0*zeta**3*f_zeta+fp_zeta*zeta4
          decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

          ec1_q1pp=0.5d0*ec1_aa*(-ec1_b1*rsm1_2**3+3.d0*ec1_b3*rsm1_2+8.d0*ec1_b4)
          d2ecrs1_drs2= 4.0d0*ec1_aa*ec1_a1*ec1_q1p*ec1_den            &
               &                -ec1_q0*ec1_q1pp*ec1_den                        &
               &                +ec1_q0*ec1_q1p**2*ec1_den**2*(2.d0*ec1_q1+1.0d0)

          mac_q1pp=0.5d0*mac_aa*(-mac_b1*rsm1_2**3+3.d0*mac_b3*rsm1_2+8.d0*mac_b4)
          d2macrs_drs2= 4.0d0*mac_aa*mac_a1*mac_q1p*mac_den            &
               &                -mac_q0*mac_q1pp*mac_den                        &
               &                +mac_q0*mac_q1p**2*mac_den**2*(2.d0*mac_q1+1.0d0)

          d2gcrs_drs2=d2ecrs1_drs2-d2ecrs0_drs2+d2macrs_drs2*fsec_inv
          fpp_zeta=(zetpm1_3(ipts)**2+zetmm1_3(ipts)**2) * factfpp_zeta
          d2fzeta4_dzeta2=12.0d0*zeta**2*f_zeta  &
               &                  + 8.0d0*zeta**3*fp_zeta &
               &                  +       zeta4  *fpp_zeta

          d2ecrs_drs2=d2ecrs0_drs2+&
               &       f_zeta*(zeta4*d2gcrs_drs2-d2macrs_drs2*fsec_inv)
          d2ecrs_drsdzeta=dfzeta4_dzeta*dgcrs_drs-fp_zeta*dmacrs_drs*fsec_inv
          d2ecrs_dzeta2=d2fzeta4_dzeta2*gcrs-fpp_zeta*macrs*fsec_inv


          ! At this point, the coding of the PW92 functional finishes.

          ! The correlation between different spin from HCTH is now computed
          ! First, the part without gradient correction factor
          rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
          vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
          drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
          drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

          ! Now, the gradient correction factor
          ssavg=half*(ss_up+ss_dn)
          divcab=one/(one+gammacab*ssavg)
          ucab=gammacab*ssavg*divcab
          ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

          gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
          dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
               &   *ducabdss

          exc=exc+rhoecab*gcab

          vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
          dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
          vxci(ipts,2)=vxci(ipts,2)+drhoecab_drhodn*gcab+rhoecab*dgcabdss*half*dssdndrho
          dvxcdgr(ipts,2)=dvxcdgr(ipts,2)+rhoecab*dgcabdss*half*dssdndg

          ! Final division by the total density, to give the energy density
          exci(ipts)=exc*rhotot_inv

       end do ! ipts=1,npts

    else

       do ipts=1,npts

          ! -----------------------------------------------------------------------
          ! First take care of the spin-split part of the functional
          exc=zero
          do ispden=1,nspden
             rho   =rho_updn(ipts,ispden)
             rhomot=rho_updnm1_3(ipts,ispden)
             rho_inv=rhomot*rhomot*rhomot

             !  Exchange part
             ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
             !  Note that this definition differs from the PBE one
             coeffss=rho_inv*rho_inv*rhomot*rhomot
             ss=grho2_updn(ipts,ispden)*coeffss
             dssdrho=-eight*third*ss*rho_inv
             dssdg=two*coeffss

             divx=one/(one+gammax*ss)
             uxsig=gammax*ss*divx
             duxsigdss=gammax*divx*(one-ss*gammax*divx)

             gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
             dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
                  &   *duxsigdss

             exc=exc+ex_lsd*rho*gxsig
             vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
             dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

             !  Spin parallel correlation part
             !  Note that this definition is for fully spin-polarized quantities
             rs=rsfac*rhomot
             rhomo6=sqrt(rhomot)
             sqr_rs=sq_rsfac*rhomo6
             rhoo6=rho*rhomot*rhomot*rhomo6
             rsm1_2=sq_rsfac_inv*rhoo6
             drsdrho=-third*rs*rho_inv

             ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
             ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
             ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
             ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
             !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
             ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
             ecrs1=ec1_q0*ec1_log
             decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
             decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

             !  Store the LSDA corr energy and density derivative
             if(ispden==1)then
                rhoecrs1_up=rho*ecrs1
                drhoecrs1_drhoup=decrs1_drho
                ss_up=ss
                dssupdrho=dssdrho
                dssupdg=dssdg
             else
                rhoecrs1_dn=rho*ecrs1
                drhoecrs1_drhodn=decrs1_drho
                ss_dn=ss
                dssdndrho=dssdrho
                dssdndg=dssdg
             end if

             divcsig=one/(one+gammacsig*ss)
             ucsig=gammacsig*ss*divcsig
             ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

             gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
             dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
                  &   *ducsigdss

             !DEBUG
             !  gcsig=zero
             !  dgcsigdss=zero
             !ENDDEBUG
             exc=exc+ecrs1*rho*gcsig
             vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
                  &                    ecrs1*rho*dgcsigdss*dssdrho
             dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

          end do

          ! -----------------------------------------------------------------------------
          ! Then takes care of the LSD correlation part of the functional

          rhotot=rhoarr(ipts)
          rhotmot=rhom1_3(ipts)
          rhotot_inv=rhotmot*rhotmot*rhotmot
          rhotmo6=sqrt(rhotmot)
          rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

          ! From now, the coding of the PW92 functional start. It is identical in xcpbe.f
          rs=rsfac*rhotmot
          sqr_rs=sq_rsfac*rhotmo6
          rsm1_2=sq_rsfac_inv*rhoto6

          ! Formulas A6-A8 of PW92LSD
          ec0_q0=-2.0d0*ec0_aa*(1.0d0+ec0_a1*rs)
          ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
          ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
          ec0_den=1.0d0/(ec0_q1*ec0_q1+ec0_q1)
          ! ec0_log=log( 1.0d0 + 1.0d0 / ec0_q1 )
          ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
          ecrs0=ec0_q0*ec0_log
          decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

          mac_q0=-2.0d0*mac_aa*(1.0d0+mac_a1*rs)
          mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
          mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
          mac_den=1.0d0/(mac_q1*mac_q1+mac_q1)
          mac_log=-log( mac_q1*mac_q1*mac_den )
          macrs=mac_q0*mac_log
          dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

          zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
          ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
          ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
          ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
          ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
          !  ec1_log=log( 1.0d0 + 1.0d0 / ec1_q1 )
          ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
          ecrs1=ec1_q0*ec1_log
          decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

          !  alpha_zeta is introduced in order to remove singularities for fully
          !  polarized systems.
          zetp_1_3=(1.0d0+zeta*alpha_zeta)*zetpm1_3(ipts)**2
          zetm_1_3=(1.0d0-zeta*alpha_zeta)*zetmm1_3(ipts)**2

          f_zeta=( (1.0d0+zeta*alpha_zeta2)*zetp_1_3 +                      &
               &           (1.0d0-zeta*alpha_zeta2)*zetm_1_3 - 2.0d0 ) * factf_zeta
          fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
          zeta4=zeta**4
          gcrs=ecrs1-ecrs0+macrs*fsec_inv
          !  ecrs=ecrs0+f_zeta*(-macrs*(1.0d0-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
          ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

          dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
          !  decrs_drs=decrs0_drs+f_zeta*&
          !&       (-dmacrs_drs*(1.0d0-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
          decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
          dfzeta4_dzeta=4.0d0*zeta**3*f_zeta+fp_zeta*zeta4
          decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

          ! At this point, the coding of the PW92 functional finishes.

          ! The correlation between different spin from HCTH is now computed
          ! First, the part without gradient correction factor
          rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
          vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
          drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
          drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

          ! Now, the gradient correction factor
          ssavg=half*(ss_up+ss_dn)
          divcab=one/(one+gammacab*ssavg)
          ucab=gammacab*ssavg*divcab
          ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

          gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
          dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
               &   *ducabdss

          exc=exc+rhoecab*gcab

          vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
          dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
          vxci(ipts,2)=vxci(ipts,2)+drhoecab_drhodn*gcab+rhoecab*dgcabdss*half*dssdndrho
          dvxcdgr(ipts,2)=dvxcdgr(ipts,2)+rhoecab*dgcabdss*half*dssdndg

          ! Final division by the total density, to give the energy density
          exci(ipts)=exc*rhotot_inv

       end do ! ipts=1,npts

    end if

 else
    ! Disallowed value for nspden
    write(message, '(a,a,a,a,a,a,i12,a)' ) ch10,&
         &   ' xchcth: BUG -',ch10,&
         &   '  Argument nspden must be 1 or 2; ',ch10,&
         &   '  Value provided as argument was ',nspden,'.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')

 end if

!DEBUG
!Finite-difference debugging, do not take away
!Beware that dvxcdgr(:,3) no longer exists
!if(debug/=0)then
! do ipts=1,npts,5

!  rho=rho_updn(ipts,1)+rho_updn(ipts,2)
!  write(6, '(a,i5,a,es16.8)' ) ' Point number',ipts,' with rho=',rho

!! For rho
!  if(debug==1)then
!   write(6, '(3es16.8)' )exci(ipts)*rho,vxci(ipts,1),vxci(ipts,2)
!  else
!!  For grho2
!   write(6, '(4es16.8)' )exci(ipts)*rho,dvxcdgr(ipts,1),&
! &  dvxcdgr(ipts,2),dvxcdgr(ipts,3)
!  end if

!  if(debug==1)then
!!  For rho
!   write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*(rho+delta) - exci(ipts+2)*(rho-delta) )/2.d0/delta,&
!&      ( exci(ipts+3)*(rho+delta) - exci(ipts+4)*(rho-delta) )/2.d0/delta
!   write(6, '(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2.d0/delta,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2.d0/delta,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2.d0/delta
!   write(6, '(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2.d0/delta
!  else
!!  For grho2  (should distinguish exchange and correlation ...)
!   grr=sqrt(grho2_updn(ipts,1)) ! Analysis of exchange
!   grr=sqrt(grho2_updn(ipts,3)) ! Analysis of correlation
!   write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*rho - exci(ipts+2)*rho )/2.d0/delta/grr,&
!&      ( exci(ipts+3)*rho - exci(ipts+4)*rho )/2.d0/delta/grr
!   write(6, '(3es16.8)' )&
! &    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2.d0/delta/grr,&
! &    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2.d0/delta/grr,&
! &    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2.d0/delta/grr
!   write(6, '(4es16.8)' )&
! &    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2.d0/delta/grr,&
! &    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2.d0/delta/grr,&
! &    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2.d0/delta/grr,&
! &    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2.d0/delta/grr
!  end if
! end do
! stop
!end if
!ENDDEBUG

 deallocate(rhoarr,rhom1_3,rho_updnm1_3)
 if(nspden==2)deallocate(zetm,zetmm1_3,zetp,zetpm1_3)

!DEBUG
!write(6,*)' xchcth : exit'
!write(6,*)' nspden=',nspden
!if(order==2)stop
!ENDDEBUG

end subroutine xchcth
!!***

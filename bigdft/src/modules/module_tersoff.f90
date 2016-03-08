!> @file
!!  Tersoff potential
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module module_tersoff
  use module_base
  implicit none

  private
  public :: tersoff
  public :: init_tersoff

     integer, parameter :: natmax=5000
     real(gp), save :: rcut1(2,2)         !cutoff
     real(gp), save :: rcut2(2,2)         !cutoff
     integer, save :: Kinds_tersoff(natmax)   !Kinds, 1 for carbon and 2 for silicon
     logical, save :: only_c
     logical, save :: initialized_tersoff=.false.

 

contains
subroutine init_tersoff(nat,astruct,paramset,paramfile,geocode)
use module_base
use yaml_output
use module_atoms
implicit none
!parameter
integer, intent(in) :: nat
type(atomic_structure), intent(in) :: astruct
character(len=*), intent(in) :: paramset                           
character(len=*), intent(in) :: paramfile                          
character(len=*), intent(in) :: geocode 
!internal
integer:: iat
    call yaml_comment('Initializing Tersoff potential',hfill='-')
    if(nat>natmax)then
        initialized_tersoff=.false.                                       
        call f_err_throw('Tersoff potential: Too many atoms '//&      
             'nat > natmax, natmax: '//yaml_toa(natmax)) 
    endif
    if(trim(geocode)/='P')then                                         
        initialized_tersoff=.false.                                       
        call f_err_throw('Tersoff potential only works with periodic '//&      
             'boundary conditions. Specified boundary conditions are: '//&
             trim(adjustl(geocode)))                                   
    endif

    only_c=.true.
do iat=1,nat
   if(trim(astruct%atomnames(astruct%iatype(iat)))=='C') then
     Kinds_tersoff(iat)=1
   elseif(trim(astruct%atomnames(astruct%iatype(iat)))=='Si') then
     Kinds_tersoff(iat)=2
     only_c=.false.
   else 
     stop "Tersoff only allowed for Silicon and/or Carbon atoms"
   endif
enddo
end subroutine

!        subroutine energyandforces_per(nat,latvec,rxyz,fxyz,stress,pressure,etot,count)
!subroutine tersoff(nat,latvec_bohr,rxyz,fxyz,strten,etot)
subroutine tersoff(nat,cell,rxyzIN,fxyz,strten,etot)
use module_base
        implicit none
        integer, intent(in) :: nat
        integer   :: iat
        real(gp)   :: etot
        real(gp)                                   ::rxyzIN(3,nat),fxyz(3,nat), alat(3)
        real(gp),              dimension(1:2,1:2)  :: R1, R2
        real(gp),              dimension(1:2,1:2)  :: Cr, Ca, alr, ala, X
        real(gp),              dimension(1:2)      :: Pmass

        integer                                   :: Nmol, Npmax, NNmax
        real(gp)                                      :: xbox, ybox, zbox
!        integer,allocatable,  dimension(:)   :: Kinds
!        integer, dimension(nat)   :: Kinds

      real(gp),allocatable,  dimension(:):: XYZRrefdf
      real(gp),allocatable,  dimension(:):: UadUrdf
      real(gp)                                      :: Urtot

      integer i,j
      integer Nptot, NNtot
    real(gp) :: rxyz(3,nat)
!      real(gp) pi
      real(gp) Xi, Yi, Zi
      real(gp) Xij, Yij, Zij, Rij, Rreij
      real(gp) R1ij, R2ij
      real(gp) PL1, PL2

      real(gp) Crij, Caij, alrij, alaij
      real(gp) fij, dfij
      real(gp) xhalf, yhalf, zhalf
      real(gp),              dimension(1:2)      ::Pn, Co_bcd, bcsq, dsq, h
      real(gp)                                   :: Uatot

      real(gp) COSijk, Gi, dGi, fdG, fdGcos
      real(gp) Co_pa, Co_mb1, Co_mb2
      real(gp) Co_cdi, Co_hcosi, Co_dhcosi
      real(gp) Pni, bcsqi, dsqi, hi

      real(gp) Bij, Eij
      real(gp) djEij
      real(gp) dXjEij2, dYjEij2, dZjEij2
      real(gp) Co1_dkEij, Co2_dkEij
      real(gp),allocatable,dimension(:)::dkEij

      integer, ALLOCATABLE, DIMENSION(:,:) :: lsta
      integer, ALLOCATABLE, DIMENSION(:) :: lstb
      integer:: nnbrx,nnbrxt

!Periodic part
      real(gp) :: cell(3)
      real(gp):: latvec_bohr(3,3),latvec(3,3),stress(3,3),strten(6),cut,F(3*nat),R(3*nat),transinv(3,3),trans(3,3),vol
      real(gp):: tmplat(3,3)
      real(gp), ALLOCATABLE, DIMENSION(:,:) :: rel
      integer:: alpha, beta, a
!        if(.not.allocated(Kinds_tersoff)) stop "Allocate Kinds before calling energyandforces" 
!        pi=dacos(-1.0_gp)

        latvec_bohr(1,1)=cell(1)
        latvec_bohr(2,1)=0.0_gp
        latvec_bohr(3,1)=0.0_gp
        latvec_bohr(1,2)=0.0_gp
        latvec_bohr(2,2)=cell(2)
        latvec_bohr(3,2)=0.0_gp
        latvec_bohr(1,3)=0.0_gp
        latvec_bohr(2,3)=0.0_gp
        latvec_bohr(3,3)=cell(3)

        latvec=latvec_bohr*Bohr_Ang
        rxyz=rxyzIN*Bohr_Ang
!        xred=xred0
        call backtocell_cart(nat,latvec,rxyz)
!        call rxyz_int2cart(latvec,xred,rxyz,nat)


        nnbrx=64
        nnbrxt=3*nnbrx/2
        nnmax=nnbrxt*nat
        npmax=nnbrxt*nat
        allocate(XYZRrefdf(1:6*Npmax), UadUrdf(1:3*Npmax), dkEij(1:3*NNmax))
        allocate(lsta(2,nat),lstb(nnbrx*nat),rel(5,nnbrx*nat))
!        do i=1,nat
!          kinds(i)=2
!        enddo
        fxyz=0.0_gp
        call Parameters_per(R1,R2,Cr,Ca,alr,ala,X,Pn,Co_bcd,bcsq,dsq,h,Pmass)
        cut=R2(2,2)-1.e-3_gp
        if(only_c) cut=R2(1,1)-1.e-3_gp
        trans=latvec
        call invertmat(trans,transinv,3)
        call makepairlist(lsta,lstb,rel,nat,rxyz,latvec,cut,nnbrx)
!end of creating pairlist part---------------------------------------------------------

!start calculating forces and energy
      call convert(rxyz,R,nat)
      Nmol=nat
      Urtot=0.0_gp
      Nptot=0

      Uatot=0.0_gp
  
      trans=latvec
      call invertmat(trans,transinv,3)

      DO_I:DO i=1, Nmol
         call subeniat_t(i,Nmol,Npmax,NNmax,X,R,R1,R2,Cr,Ca,alr,ala,XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,rel,pi,stress,transinv)
      END DO DO_I

      Urtot = 0.5_gp*Urtot
!force------------
      F=0.0_gp
      stress=0.0_gp   

      DO_If:DO i=1, Nmol
         call subfiat_t(i,Nmol,Npmax,NNmax,Pn,Co_bcd,bcsq,dsq,h,XYZRrefdf,UadUrdf,F,Uatot,dkEij,lsta,lstb,nnbrx,stress,transinv)
      END DO DO_If
      stress=stress*0.5_gp
      F     = 0.5_gp*F
      Uatot = 0.5_gp*Uatot
      etot=Urtot+Uatot
!------------------------------------------

        call convert(F,fxyz,nat)
        etot=etot/Ha_eV
        fxyz=fxyz/Ha_eV*Bohr_Ang
!Stress
!Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
        do i=1,3
           tmplat(:,i)=latvec(i,:)
        enddo
        call getvol(latvec,vol)
        stress=-matmul(stress,tmplat)/vol
        strten(1) = stress(1,1)
        strten(2) = stress(2,2)
        strten(3) = stress(3,3)
        strten(6) = stress(2,1)
        strten(5) = stress(3,1)
        strten(4) = stress(3,2)
!This is not very clear yet...
        strten=strten/Ha_eV*Bohr_Ang**3
end subroutine

!----------------------------------------------------------------------------!

      Subroutine Parameters_per(R1,R2,Cr,Ca,alr,ala,X,Pn,Co_bcd,bcsq,dsq,h,Pmass)
      implicit none
      real(gp), intent(out), dimension(1:2,1:2) :: R1, R2
      real(gp), intent(out), dimension(1:2,1:2) :: Cr, Ca, alr, ala, X
      real(gp), intent(out), dimension(1:2) :: Pn, Co_bcd, bcsq, dsq, h
      real(gp), intent(out), dimension(1:2) :: Pmass

      real(gp), parameter :: C_Cr = 1.3936e+3_gp
      real(gp), parameter :: C_Ca = 3.4674e+2_gp
      real(gp), parameter :: C_alr= 3.4879_gp
      real(gp), parameter :: C_ala= 2.2119_gp
      real(gp), parameter :: C_b  = 1.5724e-7_gp
      real(gp), parameter :: C_n  = 7.2751e-1_gp
      real(gp), parameter :: C_c  = 3.8049e+4_gp
      real(gp), parameter :: C_d  = 4.3484_gp
      real(gp), parameter :: C_h  =-5.7058e-1_gp
      real(gp), parameter :: C_R1 = 1.8_gp
      real(gp), parameter :: C_R2 = 2.1_gp
      real(gp), parameter :: C_mass = 12.0_gp

      real(gp), parameter :: Si_Cr = 1.8308e+3_gp
      real(gp), parameter :: Si_Ca = 4.7118e+2_gp
      real(gp), parameter :: Si_alr= 2.4799_gp
      real(gp), parameter :: Si_ala= 1.7322_gp
      real(gp), parameter :: Si_b  = 1.1000e-6_gp
      real(gp), parameter :: Si_n  = 7.8734e-1_gp
      real(gp), parameter :: Si_c  = 1.0039e+5_gp
      real(gp), parameter :: Si_d  = 1.6217e+1_gp
      real(gp), parameter :: Si_h  =-5.9825e-1_gp
      real(gp), parameter :: Si_R1 = 2.7_gp
      real(gp), parameter :: Si_R2 = 3.0_gp
      real(gp), parameter :: Si_mass = 28.0855_gp


      Cr(1,1) =  C_Cr
      Cr(2,2) = Si_Cr
      Cr(1,2) = dsqrt( Cr(1,1)*Cr(2,2) )
      Cr(2,1) = Cr(1,2)

      Ca(1,1) =  C_Ca
      Ca(2,2) = Si_Ca
      Ca(1,2) = dsqrt( Ca(1,1)*Ca(2,2) )
      Ca(2,1) = Ca(1,2)

      R1(1,1) =  C_R1
      R1(2,2) = Si_R1
      R1(1,2) = dsqrt( R1(1,1)*R1(2,2) )
      R1(2,1) = R1(1,2)

      R2(1,1) =  C_R2
      R2(2,2) = Si_R2
      R2(1,2) = dsqrt( R2(1,1)*R2(2,2) )
      R2(2,1) = R2(1,2)

      X(1,1)  = 1.0_gp
      X(2,2)  = 1.0_gp
      X(1,2)  = 0.9776_gp
      X(2,1)  = 0.9776_gp

      alr(1,1) =  C_alr
      alr(2,2) = Si_alr
      alr(1,2) = 0.5_gp*( alr(1,1)+alr(2,2) )
      alr(2,1) = alr(1,2)

      ala(1,1) =  C_ala
      ala(2,2) = Si_ala
      ala(1,2) = 0.5_gp*( ala(1,1)+ala(2,2) )
      ala(2,1) = ala(1,2)

      Pn(1) =  C_n
      Pn(2) = Si_n

      Co_bcd(1) =  C_b*( 1.0_gp + C_c*C_c  /(C_d*C_d)   )
      Co_bcd(2) = Si_b*( 1.0_gp + Si_c*Si_c/(Si_d*Si_d) )

      bcsq(1) = C_b  * C_c  * C_c
      bcsq(2) = Si_b * Si_c * Si_c

      dsq(1) =  C_d *  C_d
      dsq(2) = Si_d * Si_d

      h(1) =  C_h
      h(2) = Si_h

      Pmass(1)=C_mass
      Pmass(2)=Si_mass
      rcut1=R1
      rcut2=R2 

      RETURN

      STOP
      END SUBROUTINE
subroutine subeniat_t(i,Nmol,Npmax,NNmax,X,R,R1,R2,Cr,Ca,alr,ala,XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,rel,pi,stress,transinv)
        implicit none
        integer::lsta(2,Nmol),lstb(nnbrx*Nmol)
        real(gp), DIMENSION(:,:) :: rel(5,nnbrx*Nmol)
      integer, intent(in)                       :: Nmol, Npmax, NNmax
!      integer, intent(in),  dimension(1:Nmol)   :: Kinds
      real(gp)   , intent(in),  dimension(1:3*Nmol) :: R
      real(gp)   , intent(in),  dimension(1:2,1:2)  :: R1, R2 
      real(gp), intent(in),     dimension(1:2,1:2)  :: Cr, Ca, alr, ala
      real(gp), intent(in),     dimension(1:2,1:2)  :: X

      real(gp)   , intent(out), dimension(1:6*Npmax):: XYZRrefdf
      real(gp)   , intent(out), dimension(1:3*Npmax):: UadUrdf
      real(gp)   , intent(out):: Urtot

      real(gp) Ur, Ua

      integer i,j,l,nnbrx
      integer Ipt3,Jpt3,Nppt3,Nppt6
      integer Ki,Kj
      integer Nptot, NNtot, nat

      real(gp) pi
      real(gp) Xi, Yi, Zi
      real(gp) Xij, Yij, Zij, Rij, Rreij
      real(gp) R1ij, R2ij
      real(gp) PL1, PL2

      real(gp) Crij, Caij, alrij, alaij
      real(gp) fij, dfij
      real(gp) xhalf, yhalf, zhalf
      real(gp) XRreik, YRreik,  ZRreik, Rreik, fik, dfik
!Periodic part
      real(gp):: stress(3,3),transinv(3,3)
      real(gp):: sil_sjl,si1_sj1,si2_sj2,si3_sj3,tkmsum,tt

!     #######################################
!     # Calculate XYZRrefdf, UadUrdf, Urtot #
!     #######################################


         nat=Nmol
         Ki=Kinds_tersoff(i)
         Ipt3=3*(i-1)
         Xi=R(Ipt3+1)
         Yi=R(Ipt3+2)
         Zi=R(Ipt3+3)
!         NNtot=0


      DO_J:  do l=lsta(1,i),lsta(2,i)

        j=lstb(l)

         Kj=Kinds_tersoff(j)
         R2ij=R2(Ki,Kj)
         Jpt3=3*(j-1)
         Rij=rel(4,l)

!      IF_Rij:IF (Rij <= R2ij) THEN

         Xij=rel(1,l)!*Rij
         Yij=rel(2,l)!*Rij
         Zij=rel(3,l)!*Rij
         Nptot=l

!         NNtot=NNtot+1
 
         Nppt3=3*(Nptot-1)
         Nppt6=6*(Nptot-1)
         Rreij=rel(5,l)

         XYZRrefdf(Nppt6+1)=Xij!*Rreij
         XYZRrefdf(Nppt6+2)=Yij!*Rreij
         XYZRrefdf(Nppt6+3)=Zij!*Rreij
         XYZRrefdf(Nppt6+4)=Rreij


         alrij = alr(Ki,Kj)
         alaij = ala(Ki,Kj)

         Ur= Cr(Ki,Kj)*dexp(-alrij*Rij)
         Ua=-Ca(Ki,Kj)*dexp(-alaij*Rij) * X(Ki,Kj)
         R1ij=R1(Ki,Kj)




      IF (Rij <= R1ij) THEN
          XYZRrefdf(Nppt6+5)=1.0_gp
          XYZRrefdf(Nppt6+6)=0.0_gp
          Urtot=Urtot+Ur
          UadUrdf(Nppt3+1)=Ua
          UadUrdf(Nppt3+2)=-alrij*Ur
          UadUrdf(Nppt3+3)=-alaij*Ua
      ELSEIF (Rij <= R2ij) THEN
          PL1=pi/(R2ij-R1ij)
          PL2=PL1*(Rij-R1ij)
          fij =0.5_gp+0.5_gp    *dcos(PL2)
          dfij=     -0.5_gp*PL1*dsin(PL2)
          XYZRrefdf(Nppt6+5)=fij
          XYZRrefdf(Nppt6+6)=dfij
          Urtot=Urtot+fij*Ur
          UadUrdf(Nppt3+1)=fij*Ua
          UadUrdf(Nppt3+2)=(dfij-alrij*fij)*Ur
          UadUrdf(Nppt3+3)=(dfij-alaij*fij)*Ua
      END IF

!      END IF IF_Rij
      END DO DO_J
end subroutine





subroutine subfiat_t(i,Nmol,Npmax,NNmax,Pn,Co_bcd,bcsq,dsq,h,XYZRrefdf,UadUrdf,F,Uatot,dkEij,lsta,lstb,nnbrx,stress,transinv)
      implicit none
      integer,intent(in)                     ::Nmol, Npmax, NNmax
      real(gp),   intent(in),dimension(1:2) ::Pn, Co_bcd, bcsq, dsq, h

      real(gp),   intent(in),dimension(1:6*Npmax)::XYZRrefdf
      real(gp),   intent(in),dimension(1:3*Npmax)::UadUrdf

      real(gp), intent(out),   dimension(1:3*Nmol)::F
      real(gp), intent(out) :: Uatot
      real(gp), dimension(1:3*NNmax)::dkEij
      integer i
      integer Ki
      integer IJ, IJpt3, IJpt6, IK, IKpt6
      integer Ipb, Ipe
      integer Ipt3, Jpt3, Kpt3
      integer Nkpt3

      real(gp) XRreij, YRreij,  ZRreij, Rreij, fij, dfij 
      real(gp) XRreik, YRreik,  ZRreik, Rreik, fik, dfik

      real(gp) COSijk, Gi, dGi, fdG, fdGcos
      real(gp) Co_pa, Co_mb1, Co_mb2
      real(gp) Co_cdi, Co_hcosi, Co_dhcosi
      real(gp) Pni, bcsqi, dsqi, hi

      real(gp) Bij, Eij
      real(gp) djEij
      real(gp) dXjEij2, dYjEij2, dZjEij2
      real(gp) Co1_dkEij, Co2_dkEij


      real(gp) dFxi, dFyi, dFzi
      real(gp) dFxj, dFyj, dFzj
      real(gp) dFxk, dFyk, dFzk

      real(gp) Ua
      integer::nnbrx
      integer::lsta(2,Nmol),lstb(nnbrx*Nmol)
!Periodic part
      real(gp):: stress(3,3),transinv(3,3)
      real(gp):: sil_sjl,si1_sj1,si2_sj2,si3_sj3,tkmsum,tt
      real(gp):: dkEij_rad(1:4*NNmax),XRrejk,YRrejk,ZRrejk,Rrejk,Rij,Rik
      integer:: Nkpt3_per
         Ipb = lsta(1,i) 
         Ipe = lsta(2,i)




         Ki=Kinds_tersoff(i)
         bcsqi=bcsq(Ki)
         dsqi=dsq(Ki)
         hi  =h(Ki)
         Pni =Pn(Ki)

         Co_cdi=Co_bcd(Ki) 

         dFxi=0.0_gp
         dFyi=0.0_gp
         dFzi=0.0_gp

      DO_J:DO ij=Ipb,Ipe,+1

         IJpt3=3*(ij-1)
         IJpt6=6*(ij-1)

         XRreij = XYZRrefdf(IJpt6+1)
         YRreij = XYZRrefdf(IJpt6+2)
         ZRreij = XYZRrefdf(IJpt6+3)
         Rreij  = XYZRrefdf(IJpt6+4)
         fij    = XYZRrefdf(IJpt6+5)
         dfij   = XYZRrefdf(IJpt6+6)

         Eij     = 0.0_gp
         djEij   = 0.0_gp
         dXjEij2 = 0.0_gp
         dYjEij2 = 0.0_gp
         dZjEij2 = 0.0_gp

         Nkpt3=-3
   !Periodic   
         Nkpt3_per=-4
         Rij=1.0_gp/Rreij

      DO_K:DO ik=Ipb,Ipe,+1

         Nkpt3=Nkpt3+3
   !Periodic
         Nkpt3_per=Nkpt3_per+4

      IKIJ:IF (ik /= ij ) THEN

         IKpt6=6*(ik-1)

         XRreik = XYZRrefdf(IKpt6+1)
         YRreik = XYZRrefdf(IKpt6+2)
         ZRreik = XYZRrefdf(IKpt6+3)
         Rreik  = XYZRrefdf(IKpt6+4)
         fik    = XYZRrefdf(IKpt6+5)
         dfik   = XYZRrefdf(IKpt6+6)

         COSijk = XRreij*XRreik + YRreij*YRreik + ZRreij*ZRreik

         Co_hcosi  = hi - COSijk
         Co_dhcosi = 1.0_gp / ( dsqi + Co_hcosi*Co_hcosi )
         Gi        = -bcsqi*Co_dhcosi
         dGi       = 2.0_gp*Co_hcosi*Co_dhcosi*Gi
         Gi        = Gi + Co_cdi

         Eij = Eij + fik*Gi

         fdG    = fik * dGi
         fdGcos = fdG * COSijk

         djEij    = djEij   + fdGcos

         dXjEij2 = dXjEij2 + fdG*XRreik
         dYjEij2 = dYjEij2 + fdG*YRreik
         dZjEij2 = dZjEij2 + fdG*ZRreik


         Co1_dkEij = -dfik*Gi + fdGcos*Rreik
         Co2_dkEij = -fdG*Rreik

         dkEij(Nkpt3+1) = Co1_dkEij * XRreik + Co2_dkEij * XRreij
         dkEij(Nkpt3+2) = Co1_dkEij * YRreik + Co2_dkEij * YRreij
         dkEij(Nkpt3+3) = Co1_dkEij * ZRreik + Co2_dkEij * ZRreij
   !Periodic
         dkEij_rad(Nkpt3_per+1)=XRreik
         dkEij_rad(Nkpt3_per+2)=YRreik
         dkEij_rad(Nkpt3_per+3)=ZRreik
         dkEij_rad(Nkpt3_per+4)=1.0_gp/Rreik

      ELSE 
         dkEij(Nkpt3+1) = 0.0_gp
         dkEij(Nkpt3+2) = 0.0_gp
         dkEij(Nkpt3+3) = 0.0_gp
   !Periodic
         dkEij_rad(Nkpt3_per+1)=0.0_gp
         dkEij_rad(Nkpt3_per+2)=0.0_gp
         dkEij_rad(Nkpt3_per+3)=0.0_gp
         dkEij_rad(Nkpt3_per+4)=0.0_gp
      END IF IKIJ

      END DO DO_K


      Bij  = 1.0_gp + Eij**Pni
      Ua   = UadUrdf(IJpt3+1) * Bij**(-0.5_gp/Pni)
      Uatot = Uatot + Ua

      Co_pa  = UadUrdf(IJpt3+2) + UadUrdf(IJpt3+3)*Bij**(-0.5_gp/Pni)


      CEij:IF ( Nkpt3 > 0 ) THEN 

         Co_mb1 = Ua * 0.5_gp*Eij**(Pni-1.0_gp) / Bij
         Co_mb2 = Co_mb1*Rreij

            Nkpt3=-3
   !Periodic   
            Nkpt3_per=-4         
         DO ik=Ipb,Ipe,+1

            Nkpt3=Nkpt3+3
   !Periodic
            Nkpt3_per=Nkpt3_per+4         
           
            dFxk = Co_mb1 * dkEij(Nkpt3+1)
            dFyk = Co_mb1 * dkEij(Nkpt3+2)
            dFzk = Co_mb1 * dkEij(Nkpt3+3)

            Kpt3=3*(lstb(ik)-1)
            F(Kpt3+1) = F(Kpt3+1) + dFxk
            F(Kpt3+2) = F(Kpt3+2) + dFyk
            F(Kpt3+3) = F(Kpt3+3) + dFzk

            dFxi = dFxi + dFxk
            dFyi = dFyi + dFyk
            dFzi = dFzi + dFzk

                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*dkEij_rad(Nkpt3_per+1)+&
                           transinv(1,2)*dkEij_rad(Nkpt3_per+2)+&
                           transinv(1,3)*dkEij_rad(Nkpt3_per+3)
                   si2_sj2=transinv(2,1)*dkEij_rad(Nkpt3_per+1)+&
                           transinv(2,2)*dkEij_rad(Nkpt3_per+2)+&
                           transinv(2,3)*dkEij_rad(Nkpt3_per+3)
                   si3_sj3=transinv(3,1)*dkEij_rad(Nkpt3_per+1)+&
                           transinv(3,2)*dkEij_rad(Nkpt3_per+2)+&
                           transinv(3,3)*dkEij_rad(Nkpt3_per+3)
                   stress(1,1)=stress(1,1)-dkEij_rad(Nkpt3_per+4)*dFxk*si1_sj1 
                   stress(1,2)=stress(1,2)-dkEij_rad(Nkpt3_per+4)*dFxk*si2_sj2 
                   stress(1,3)=stress(1,3)-dkEij_rad(Nkpt3_per+4)*dFxk*si3_sj3 
                   stress(2,1)=stress(2,1)-dkEij_rad(Nkpt3_per+4)*dFyk*si1_sj1 
                   stress(2,2)=stress(2,2)-dkEij_rad(Nkpt3_per+4)*dFyk*si2_sj2 
                   stress(2,3)=stress(2,3)-dkEij_rad(Nkpt3_per+4)*dFyk*si3_sj3 
                   stress(3,1)=stress(3,1)-dkEij_rad(Nkpt3_per+4)*dFzk*si1_sj1 
                   stress(3,2)=stress(3,2)-dkEij_rad(Nkpt3_per+4)*dFzk*si2_sj2 
                   stress(3,3)=stress(3,3)-dkEij_rad(Nkpt3_per+4)*dFzk*si3_sj3 

         END DO

         dFxj = Co_pa*XRreij + Co_mb2*( XRreij*djEij - dXjEij2 )
         dFyj = Co_pa*YRreij + Co_mb2*( YRreij*djEij - dYjEij2 )
         dFzj = Co_pa*ZRreij + Co_mb2*( ZRreij*djEij - dZjEij2 )

      ELSE 

         dFxj = Co_pa*XRreij 
         dFyj = Co_pa*YRreij 
         dFzj = Co_pa*ZRreij 

      END IF CEij



      Jpt3=3*(lstb(ij)-1)

      F(Jpt3+1) = F(Jpt3+1) + dFxj
      F(Jpt3+2) = F(Jpt3+2) + dFyj
      F(Jpt3+3) = F(Jpt3+3) + dFzj

      dFxi = dFxi + dFxj
      dFyi = dFyi + dFyj
      dFzi = dFzi + dFzj



                   !Cell gradient part
                   !My own implementation
                   si1_sj1=transinv(1,1)*XRreij+transinv(1,2)*YRreij+transinv(1,3)*ZRreij
                   si2_sj2=transinv(2,1)*XRreij+transinv(2,2)*YRreij+transinv(2,3)*ZRreij
                   si3_sj3=transinv(3,1)*XRreij+transinv(3,2)*YRreij+transinv(3,3)*ZRreij
                   stress(1,1)=stress(1,1)-Rij*dFxj*si1_sj1
                   stress(1,2)=stress(1,2)-Rij*dFxj*si2_sj2
                   stress(1,3)=stress(1,3)-Rij*dFxj*si3_sj3
                   stress(2,1)=stress(2,1)-Rij*dFyj*si1_sj1
                   stress(2,2)=stress(2,2)-Rij*dFyj*si2_sj2
                   stress(2,3)=stress(2,3)-Rij*dFyj*si3_sj3
                   stress(3,1)=stress(3,1)-Rij*dFzj*si1_sj1
                   stress(3,2)=stress(3,2)-Rij*dFzj*si2_sj2
                   stress(3,3)=stress(3,3)-Rij*dFzj*si3_sj3


      END DO DO_J

         Ipt3 = 3*(i-1)
         F(Ipt3+1) = F(Ipt3+1) - dFxi
         F(Ipt3+2) = F(Ipt3+2) - dFyi
         F(Ipt3+3) = F(Ipt3+3) - dFzi
end subroutine


subroutine convert(Fin,Fout,nat)
implicit none
real(gp):: Fin(3*nat),Fout(3*nat)
integer:: nat
Fout=Fin
end subroutine

 subroutine reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
 !This subroutine will calculate neighbors of atom iat within the given periodic image cell and the relative coordinates
 !The relative informations are stored in relinteriat
 !The atomic index jat is stored in lstbiat
 !rxyz holds the positions of iat, rxyzp is one of the periodic images to check
 !d2plane,d2corner,d2egde hold the distance of atom iat to the walls of the cell
 real*8, intent(inout)  :: relinteriat(5,nmax)
 real*8, intent(in)     :: rxyz(3,nat),rxyzp(3,nat),cut
 integer, intent(in)    :: nat,iat,nmax
 integer, intent(inout) :: ninter,lstbiat(nmax)
 real*8                 :: cut2,xrel,yrel,zrel,rr2,tt,tti
 integer                :: jat,Ki,Kj
 Ki=Kinds_tersoff(iat)
! cut2=cut*cut
 do jat=1,nat
 Kj=Kinds_tersoff(jat)
 cut2=(rcut2(Ki,Kj)-1.d-8)*(rcut2(Ki,Kj)-1.d-8)
 xrel= rxyz(1,iat)-rxyzp(1,jat)
 yrel= rxyz(2,iat)-rxyzp(2,jat)
 zrel= rxyz(3,iat)-rxyzp(3,jat)
 rr2=xrel**2 + yrel**2 + zrel**2
! if (rr2.le.cut2 .AND. iat.ne.jat) then
 if (rr2.lt.cut2 .AND. rr2.ne.0.0_gp) then
! write(*,*) xrel,yrel,zrel
  ninter=ninter+1
  if (ninter.gt.nmax) stop 'nmax exceeded'
  lstbiat(ninter)=jat
  tt=dsqrt(rr2)
  tti=1.0_gp/tt
  relinteriat(1,ninter)=xrel*tti
  relinteriat(2,ninter)=yrel*tti
  relinteriat(3,ninter)=zrel*tti
  relinteriat(4,ninter)=tt
  relinteriat(5,ninter)=tti
  
!  write(*,*) 'Inside reliatp,jat found',jat,tt 
 endif
 enddo
  
 
 return
 end subroutine



!********************************************************

 subroutine makepairlist(lsta,lstb,rel,nat,rxyz,latvec,cut,nmax)
 !This subroutine will construct the pairlist lsta and lstb of the atomic positions rxyz in the cell latvec
 !The cutoff is given with cut, all relative coordinates are stored in rel
 implicit none
 integer, intent(in) :: nat
 integer, intent(inout) :: lsta(2,nat),lstb(nat*nmax),nmax
 real*8, intent(inout)  :: rel(5,nat*nmax)
 real*8, intent(in)     :: cut, rxyz(3,nat), latvec(3,3)
 real*8                 :: rxyzexp(3,nat,3,3,3),relinteriat(5,nmax),transvecall(3,3,3,3)
 integer                :: iat,lstbiat(nmax),infocode,ninter,k,l,m
 
 call expand(rxyz,rxyzexp,transvecall,latvec,nat) 
 do iat=1,nat
 ninter=0
    lstbiat=0
    relinteriat=0.0_gp
    call findinter(iat,nat,lstbiat,cut,rxyz,rxyzexp,transvecall,latvec,ninter,nmax,relinteriat,infocode)
    lsta(1,iat)=(iat-1)*nmax+1
    lsta(2,iat)=(iat-1)*nmax+ninter
    lstb(lsta(1,iat):lsta(1,iat)+nmax-1)=lstbiat(:)
    rel(:,lsta(1,iat):lsta(1,iat)+nmax-1)=relinteriat(:,:)
 enddo
 
! do iat=1,nat
!    k=lsta(1,iat)
!    l=lsta(2,iat)
!    do m=k,l
!    write(*,*) iat,lstb(m),m
!    enddo
! enddo
end subroutine


!**************************************************************************************
 subroutine findinter(iat,nat,lstbiat,cut,rxyz,rxyzexp,transvecall,latvec,ninter,nmax,relinteriat,infocode) 
 !This subroutine will find all neighboring atoms of atom iat including those in periodic
 !images within the range of cut. The cell is described by latvec[v1,v2,v3]. rxyz holds the
 !coordinates of all atoms, nat is the number of atoms. The output ist ninter=the number of
 !atoms in the list of neighbors. Relative coordinates are stored in relinter,
 !the index of the interactiong atoms in lstbiat. rxyz holds the pos of all periodic images
 !infocode vill return 1 if an error occured and the number of neighbors exceeds nmax, else 0
 real*8,  intent(in)  :: rxyz(3,nat), latvec(3,3), cut, rxyzexp(3,nat,3,3,3),transvecall(3,3,3,3)
 real*8, intent(inout):: relinteriat(5,nmax)
 integer, intent(in)  :: iat,nat,nmax
 integer,intent(inout):: ninter,infocode,lstbiat(nmax)
 real(gp)              :: a(3),b(3),crossp(3),d2plane(6),d2corner(gp),d2edge(12),d2all(3,3,3),ppoint0(3),tempvec(3)
 real(gp)              :: nvectmp(3),pointtmp(3),ppointtmp(3),latvectmp(3),dist,nvec(3,3),cutplus
 real(gp)              :: tempvece(3,3,3,3,2),tempvecp(3,3,3,3,2)
 integer              :: k,l,m,summe,kt,lt,mt
 real(gp)              :: rxyzp(3,nat)
 !For second nearest images
 real(gp),allocatable  :: d2all_2(:,:,:),rxyzexp_2(:,:,:,:,:),transvecall_2(:,:,:,:)
 !For third nearest images
 real(gp),allocatable  :: d2all_3(:,:,:),rxyzexp_3(:,:,:,:,:),transvecall_3(:,:,:,:)
 !For the fourth nearest images
 real(gp),allocatable  :: d2all_4(:,:,:),rxyzexp_4(:,:,:,:,:),transvecall_4(:,:,:,:)
 real(gp)              :: dx2,dy2,dz2,dd2,cut2,aa,bb,thrd
 integer              :: k1,l1,m1
 real(gp)              :: dproj(6),rotmat(3,3)
 logical              :: nec


 ninter=0 
 ppoint0=(/0.0_gp,0.0_gp,0.0_gp/)
 cutplus=cut+cut*0.01_gp
 d2all=cutplus
 tempvece=0.0_gp

 call nveclatvec(latvec,nvec)
    !Bottom
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall(:,1,1,2),dist)
    d2all(:,:,1)=dist
    !Top
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall(:,1,1,3),dist)
    d2all(:,:,3)=dist

    !Back
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall(:,2,1,1),dist)
    d2all(1,:,:)=dist
    !Front
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall(:,3,1,1),dist)
    d2all(3,:,:)=dist

    !Left
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall(:,1,2,1),dist)
    d2all(:,1,:)=dist
    !Right
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall(:,1,3,1),dist)
    d2all(:,3,:)=dist


 d2all(2,2,2)=0.0_gp !Always calculate in original cell
 !Calculate all information of relinter for atoms within the cell
 ninter=0
 do k=1,3
  do l=1,3
   do m=1,3
      if (d2all(k,l,m).lt.cut) then
      
      rxyzp(:,:)=rxyzexp(:,:,k,l,m)
!      write(*,*) 'calling reliatp' ,iat,k,l,m
!      write(*,*) "ninter1",ninter
      call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
!      write(*,*) "ninter2",ninter
      endif  
   enddo
  enddo
 enddo
 infocode=0
! write(*,*) 'all well'

! return

 !If second layer neccesary
 !Call subroutine to check here...
 call check_rep(latvec,cut,nec)
 if(nec) then
    return
 else
!    write(*,*) "entering layer 2"
    allocate(d2all_2(5,5,5),rxyzexp_2(3,nat,5,5,5),transvecall_2(3,5,5,5))
    call expand_2(rxyz,rxyzexp_2,transvecall_2,latvec,nat)
! call latvec2dproj(dproj,latvec,rotmat,rxyz,nat) 
! open(unit=12,file="secondexp.ascii")
! write(12,*) nat
! write(12,*) dproj(1),dproj(2),dproj(3)
! write(12,*) dproj(4),dproj(5),dproj(6)
! do k=1,5
! do l=1,5
! do m=1,5
! do jat=1,nat
! write(12,*) rxyzexp_2(:,jat,k,l,m),"au"
! enddo
! enddo
! enddo
! enddo
! close(12)
!!! stop

    d2all_2=cutplus
    !Calculate interaction with second nearest cells.
    !Check planes
    !Bottom
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_2(:,1,1,2),dist)
    d2all_2(:,:,1)=dist
    !Top
    call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_2(:,1,1,5),dist)
    d2all_2(:,:,5)=dist

    !Back
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_2(:,2,1,1),dist)
    d2all_2(1,:,:)=dist
    !Front
    call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_2(:,5,1,1),dist)
    d2all_2(5,:,:)=dist

    !Left
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_2(:,1,2,1),dist)
    d2all_2(:,1,:)=dist
    !Right
    call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_2(:,1,5,1),dist)
    d2all_2(:,5,:)=dist

 
!    !Substract corners
!    bb=5.0_gp*(1.0_gp-4.0_gp/3.0_gp)
!    aa=4.0_gp/3.0_gp
!    cut2=cut*cut
!    do k=2,5,3
!       do l=2,5,3
!          do m=2,5,3
!          dx2=rxyz(1,iat)-transvecall_2(1,k,l,m) 
!          dx2=dx2*dx2
!          dy2=rxyz(2,iat)-transvecall_2(2,k,l,m)
!          dy2=dy2*dy2
!          dz2=rxyz(3,iat)-transvecall_2(3,k,l,m)
!          dz2=dz2*dz2
!          dd2=dx2+dy2+dz2
!             if(dd2.gt.cut2) then
!              k1=int(bb+aa*real(k,8))
!              l1=int(bb+aa*real(l,8))
!              m1=int(bb+aa*real(m,8))
!              d2all_2(k1,l1,m1)=cutplus
!             endif   
!          enddo
!       enddo
!    enddo 
   
    do k=1,5
     do l=1,5
      do m=1,5
         if (d2all_2(k,l,m).ge.cut) then
         goto 1010
         else                        
         rxyzp(:,:)=rxyzexp_2(:,:,k,l,m)
   !      write(*,*) 'calling reliatp' ,iat,k,l,m
         call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
         endif
         1010 continue
      enddo
     enddo
    enddo
    infocode=0
    deallocate(d2all_2,rxyzexp_2,transvecall_2)
    nec=.false.
    call check_rep(latvec,1.0_gp/2.0_gp*cut,nec)

    if(nec) then
       return
    else

!    write(*,*) "entering layer 3"
       !Calculate third layer of periodic images if neccesary 
       allocate(d2all_3(7,7,7),rxyzexp_3(3,nat,7,7,7),transvecall_3(3,7,7,7))
       call expand_3(rxyz,rxyzexp_3,transvecall_3,latvec,nat)
! call latvec2dproj(dproj,latvec,rotmat,rxyz,nat) 
! open(unit=12,file="secondexp3.ascii")
! write(12,*) nat
! write(12,*) dproj(1),dproj(2),dproj(3)
! write(12,*) dproj(4),dproj(5),dproj(6)
! do k=1,7
! do l=1,7
! do m=1,7
! do jat=1,nat
! write(12,*) rxyzexp_3(:,jat,k,l,m),"au"
! enddo
! enddo
! enddo
! enddo
! close(12)
! stop


!       write(*,*) "after expand 3"
       d2all_3=cutplus
       !Calculate interaction with third nearest cells.
       !Check planes
       !Bottom
       call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_3(:,1,1,2),dist)
       d2all_3(:,:,1)=dist
       !Top
       call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_3(:,1,1,7),dist)
       d2all_3(:,:,7)=dist
   
       !Back
       call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_3(:,2,1,1),dist)
       d2all_3(1,:,:)=dist
       !Front
       call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_3(:,7,1,1),dist)
       d2all_3(7,:,:)=dist
   
       !Left
       call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_3(:,1,2,1),dist)
       d2all_3(:,1,:)=dist
       !Right
       call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_3(:,1,7,1),dist)
       d2all_3(:,7,:)=dist
   
    
!       !Substract corners
!       bb=7.0_gp*(1.0_gp-6.0_gp/5.0_gp)
!       aa=6.0_gp/5.0_gp
!       cut2=cut*cut
!       do k=2,7,5
!          do l=2,7,5
!             do m=2,7,5
!             dx2=rxyz(1,iat)-transvecall_3(1,k,l,m) 
!             dx2=dx2*dx2
!             dy2=rxyz(2,iat)-transvecall_3(2,k,l,m)
!             dy2=dy2*dy2
!             dz2=rxyz(3,iat)-transvecall_3(3,k,l,m)
!             dz2=dz2*dz2
!             dd2=dx2+dy2+dz2
!                if(dd2.gt.cut2) then
!                 k1=int(bb+aa*real(k,8))
!                 l1=int(bb+aa*real(l,8))
!                 m1=int(bb+aa*real(m,8))
!                 d2all_3(k1,l1,m1)=cutplus
!                endif   
!             enddo
!          enddo
!       enddo 
      
       do k=1,7
        do l=1,7
         do m=1,7
            if (d2all_3(k,l,m).ge.cut) then
            goto 1011
            else                        
            rxyzp(:,:)=rxyzexp_3(:,:,k,l,m)
      !      write(*,*) 'calling reliatp' ,iat,k,l,m
            call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
            endif
            1011 continue
         enddo
        enddo
       enddo
       infocode=0
       deallocate(d2all_3,rxyzexp_3,transvecall_3)
       nec=.false.
       call check_rep(latvec,1.0_gp/3.0_gp*cut,nec)
       if(nec) then
       return
       else
          !Calculate fourth layer of periodic images if neccesary 
          allocate(d2all_4(9,9,9),rxyzexp_4(3,nat,9,9,9),transvecall_4(3,9,9,9))
          call expand_4(rxyz,rxyzexp_4,transvecall_4,latvec,nat)
          d2all_4=cutplus
          !Calculate interaction with fourth nearest cells.
          !Check planes
          !Bottom
          call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_4(:,1,1,2),dist)
          d2all_4(:,:,1)=dist
          !Top
          call dist2plane(rxyz(:,iat),nvec(:,1),transvecall_4(:,1,1,9),dist)
          d2all_4(:,:,9)=dist
   
          !Back
          call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_4(:,2,1,1),dist)
          d2all_4(1,:,:)=dist
          !Front
          call dist2plane(rxyz(:,iat),nvec(:,2),transvecall_4(:,9,1,1),dist)
          d2all_4(9,:,:)=dist
   
          !Left
          call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_4(:,1,2,1),dist)
          d2all_4(:,1,:)=dist
          !Right
          call dist2plane(rxyz(:,iat),nvec(:,3),transvecall_4(:,1,9,1),dist)
          d2all_4(:,9,:)=dist


          do k=1,9
           do l=1,9
            do m=1,9
               if (d2all_4(k,l,m).ge.cut) then
               goto 1012
               else
               rxyzp(:,:)=rxyzexp_4(:,:,k,l,m)
         !      write(*,*) 'calling reliatp' ,iat,k,l,m
               call reliatp(cut,relinteriat,lstbiat,rxyz,rxyzp,nat,iat,ninter,nmax)
               endif
               1012 continue
            enddo
           enddo
          enddo
          infocode=0
          deallocate(d2all_4,rxyzexp_4,transvecall_4)
          nec=.false.
          call check_rep(latvec,1.0_gp/4.0_gp*cut,nec)
          if (nec) then
          return
          else
          write(*,*)  "Cell periodicit exceeded the fourth order"
          stop
          endif
      endif 

    endif


 endif

 end subroutine

 subroutine expand(rxyz,rxyzout,transvecall,latvec,nat)
 !This subroutine will expand the unit cell into 26 periodic cells and store them in rxyzout
 implicit none
 real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
 integer, intent(in) :: nat
 real*8, intent(out) :: rxyzout(3,nat,3,3,3) !26 periodic images plus the main cell
 integer             :: iat,iplane,icorner,iedge,m,k,l
 real*8,intent(inout):: transvecall(3,3,3,3)!,(transvecp(3,6),transvecc(3,8),transvece(3,12)

 do m=-1,1
    do k=-1,1
       do l=-1,1
       transvecall(:,l+2,k+2,m+2)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 do m=1,3
    do k=1,3
       do l=1,3
       do iat=1,nat
       rxyzout(:,iat,m,k,l)=rxyz(:,iat)+transvecall(:,m,k,l)
       enddo
       enddo
    enddo
 enddo
 end subroutine


 
 subroutine expand_2(rxyz,rxyzout_2,transvecall_2,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(gp), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(gp), intent(inout):: rxyzout_2(3,nat,5,5,5)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(gp),intent(inout) :: transvecall_2(3,5,5,5)

 do m=-2,2
    do k=-2,2
       do l=-2,2
       transvecall_2(:,l+3,k+3,m+3)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3) 
       enddo
    enddo
 enddo

 rxyzout_2=0.0_gp
 do m=1,5
    do k=1,5
       do l=1,5
       do iat=1,nat
!       if(m==5 .or. k==5 .or. l==5 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_2(:,iat,m,k,l)=rxyz(:,iat)+transvecall_2(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine

 subroutine expand_3(rxyz,rxyzout_3,transvecall_3,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(gp), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(gp), intent(inout):: rxyzout_3(3,nat,7,7,7)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(gp),intent(inout) :: transvecall_3(3,7,7,7)

 do m=-3,3
    do k=-3,3
       do l=-3,3
       transvecall_3(:,l+4,k+4,m+4)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 rxyzout_3=0.0_gp
 do m=1,7
    do k=1,7
       do l=1,7
       do iat=1,nat
!       if(m==7 .or. k==7 .or. l==7 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_3(:,iat,m,k,l)=rxyz(:,iat)+transvecall_3(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine


 subroutine expand_4(rxyz,rxyzout_4,transvecall_4,latvec,nat)
 !This subroutine will expand the periodic images to the second order, thus creating the periodic cells of 5*5*5 images
 !The 3*3*3 centered cells are nulled bc thay are present in the expansion of the first order (first neighboring cells)
 implicit none
 real(gp), intent(in)   :: latvec(3,3),rxyz(3,nat)
 real(gp), intent(inout):: rxyzout_4(3,nat,9,9,9)
 integer, intent(in)   :: nat
 integer               :: iat,iplane,icorner,iedge,m,k,l
 real(gp),intent(inout) :: transvecall_4(3,9,9,9)

 do m=-4,4
    do k=-4,4
       do l=-4,4
       transvecall_4(:,l+5,k+5,m+5)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 rxyzout_4=0.0_gp
 do m=1,9
    do k=1,9
       do l=1,9
       do iat=1,nat
!       if(m==7 .or. k==7 .or. l==7 .or. m==1 .or. k==1 .or. l==1) then
       rxyzout_4(:,iat,m,k,l)=rxyz(:,iat)+transvecall_4(:,m,k,l)
!       endif
       enddo
       enddo
    enddo
 enddo

 end subroutine

 subroutine check_rep(latvec,cut,nec)
 !This subroutine will check if the cell size is sufficentli large to have periodic boundary conditions
 !with for the given cut. If true, it fits. If false, it does not fit--> new layer neccesary 
 implicit none
 real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist,eps
 integer:: i
 logical:: nec
! eps=1.d-6
 nec=.true.
 call nveclatvec(latvec,nvec)
 point0=(/0.0_gp,0.0_gp,0.0_gp/)
 do i=1,3
 call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist)
! write(*,*) "cut",i,cut, dist
 if (dist.le.cut) nec=.false.
 enddo
 end subroutine




 subroutine backtocell(nat,latvec,rxyz_red)
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 use module_base
 implicit none
 integer:: nat,i,iat,j
 real(gp) :: latvec(3,3), rxyz_red(3,nat), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(gp) :: v(3,3),vol
 logical:: neccesary
!First check if the volume is positive
 call getvol(latvec,vol)
 if(vol.le.0.0_gp) stop "Negative volume during backtocell"
 do iat=1,nat
        rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.0_gp),1.0_gp)
        rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.0_gp),1.0_gp)
        rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.0_gp),1.0_gp)
 enddo
 end subroutine
    subroutine getvol(latvec,vol)
 use module_base
    implicit none
    real(gp):: latvec(3,3),v(3,3),vol
     v=latvec
     vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
          v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
    end subroutine
 subroutine backtocell_cart(nat,latvec,rxyz)
 use module_base
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 implicit none
 integer:: nat,i,iat,j
 real(gp) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(gp) :: v(3,3),vol,rxyz_red(3,nat)
 logical:: neccesary
 call getvol(latvec,vol)
   !First check if the volume is positive
   if(vol.le.0.0_gp) stop "Negative volume during backtocell"
   call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
   do iat=1,nat
       rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.0_gp),1.0_gp)
       rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.0_gp),1.0_gp)
       rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.0_gp),1.0_gp)
   enddo
   call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
end subroutine

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 use module_base
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(gp):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int
 subroutine invertmat(mat,matinv,n)
 use module_base
 implicit none
 real(gp),intent(in) :: mat(n,n)
 integer               :: n
 real(gp),allocatable   :: WORK(:)
 real(gp)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1))
 div=1.0_gp/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
 !General n*n matrix 
 matinv=mat
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=int(WORK(1))
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 endif
 end subroutine
 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(gp):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
 do iat=1,nat
  rxyzcart(:,iat)=matmul(latvec,rxyzint(:,iat))
 enddo
 end subroutine rxyz_int2cart
 subroutine nveclatvec(latvec,nvec)
 use module_base
 !Will calculate the normalized normal vector to the 3 planes of the cell
 implicit none
 real(gp), intent(in) :: latvec(3,3)
 real(gp), intent(out):: nvec(3,3)
 real(gp)             :: a(3),b(3),crossp(3),norm
 integer:: i
 do i=1,3
 a=latvec(:,i)
 b=latvec(:,mod(i,3)+1)
 call cross_product(a,b,crossp)
 norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 nvec(:,i)=crossp(:)/norm
 enddo
 end subroutine
!**********************************************************************************************

 subroutine dist2plane(point,nvec,ppoint,dist)
 !This subroutine will calculate  the distance between a plane and a point in space
 !The point is 'point', the normalized normal vector of the plane is 'nvec', 'ppoint' is an arbitrary point on the plane
 !and the output is the distance 'dist'  
 real(gp), intent(in) :: point(3),nvec(3),ppoint(3)
 real(gp), intent(out):: dist
 integer             :: i
 real(gp)             :: p,nvectmp(3)
 nvectmp(:)=nvec(:)!/sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
 p=DOT_PRODUCT(nvectmp,ppoint)
 p=-p
 dist=DOT_PRODUCT(nvectmp,point)+p
 end subroutine
 subroutine cross_product(a,b,crossp)
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

 subroutine dot_p(a,b,dotp)
 !a very simple implementation of the dot product
 implicit none
 real(8)::a(3),b(3)
 real(8)::dotp(3)
 integer::i
 dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 end subroutine

end module module_tersoff

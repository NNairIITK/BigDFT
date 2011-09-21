!> @file
!!  Data routines for electronic configuration of the atoms
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (TD,LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>   Give electronic configuration of atom
!! SYNOPSIS
!!  Input
!!   @param nzatom    Z number of atom
!!   @param nvalelec  Number of valence electrons
!!  Output
!!   @param symbol    Atomic symbol
!!   @param rcov      Covalent radius
!!   @param rprb      Parabolic radius for the input guess using the subroutine "gatom"
!!   @param ehomo     Highest occupied molecular orbital energy
!!                    See <a>http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html</a>
!!   @param neleconf  Occupation number (electronic configuration of the atom)
!!   @param nsccode   Semicore orbitals, indicated as an integer.
!!                    The integer is the n_s + 4*n_p + 16* n_d + 64* n_f
!!                    where n_l are the number of semicore orbitals for a given angular momentum
!!                    starting from the lower level of course
!!   @param mxpl      Maximum spin polarisation to be placed on the atom
!!   @param mxchg     Maximum charge to be placed on the atom
!!   @param amu       Atomic mass unit (use values coming from ABINIT/11util/atmdata.F90)
subroutine eleconf(nzatom,nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
  implicit none
! Arguments
  integer, intent(in) :: nzatom,nvalelec
  character(len=2), intent(out) :: symbol
  real(kind=8), intent(out) :: rcov,rprb,ehomo,amu
  integer, parameter :: nmax=6,lmax=3
  real(kind=8), intent(out) :: neleconf(nmax,0:lmax)
  integer, intent(out) :: nsccode,mxpl,mxchg
! Local variables
  integer :: n,l,nsum,ipow,lsc,inorbsc,i
  real(kind=8) :: sccode

  neleconf(:,:)=0
  nsccode=0

!Each atomic configuration
select case(nzatom*1000+nvalelec)

case(1*1000+1)
! --------------------           1
! H            1           1     Symbol, Z, Zion
! 0.75 1.21     rad:  covalent, parab. potential
!    n=1    n=2    n=3    n=4    n=5    n=6
!      1      0      0      0      0      0  l=0
!      0      0      0      0      0      0  l=1
!      0      0      0      0      0      0  l=2
!      0      0      0      0      0      0  l=3
symbol = "H"
rcov=0.75d0
rprb=1.21d0
ehomo=-0.233471d0
neleconf(1,0)=1
amu=1.00794d0

case(2*1000+2)
! -----------------------           2
! He           2           2     Symbol, Z, Zion
symbol = "He"
rcov=0.75d0
rprb=1.50d0
ehomo=-0.570425d0
neleconf(1,0)=2
amu=4.002602d0

case(3*1000+1)
! -----------------------           3
! Li           3           1     Symbol, Z, Zion
symbol = "Li"
rcov=3.40d0
rprb=6.80d0
ehomo=-0.10554d0
neleconf(2,0)=1
amu=6.941d0

case(3*1000+3)
! -----------------------           4
! Li           3           3     Symbol, Z, Zion
symbol = "Li"
rcov=1.60d0
rprb=3.61d0
ehomo=-0.10554d0
neleconf(1,0)=2
neleconf(2,0)=1
nsccode=1
amu=6.941d0

case(4*1000+2)
! -----------------------           5
! Be           4           2     Symbol, Z, Zion
symbol = "Be"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.205744d0
neleconf(2,0)=2
amu=9.012182d0

case(4*1000+4)
! -----------------------           6
! Be           4           4     Symbol, Z, Zion
symbol = "Be"
rcov=1.30d0
rprb=3.60d0
ehomo=-0.205744d0
neleconf(1,0)=2
neleconf(2,0)=2
nsccode=1
amu=9.012182d0

case(5*1000+3)
! -----------------------           7
! B            5           3     Symbol, Z, Zion
symbol = "B"
rcov=1.55d0
rprb=3.10d0
ehomo=-0.136603d0
neleconf(2,0)=2
neleconf(2,1)=1
amu=10.811d0

case(6*1000+4)
! -----------------------           8
! C            6           4     Symbol, Z, Zion
symbol = "C"
rcov=1.45d0
rprb=2.90d0
ehomo=-0.199186d0
neleconf(2,0)=2.d0
neleconf(2,1)=2.d0
!neleconf(3,0)=1.d-100
!neleconf(3,1)=1.d-100
!neleconf(3,2)=1.d-100
amu=12.011d0

case(6*1000+6)
! -----------------------           8
! C            6           6     Symbol, Z, Zion
symbol = "C"
rcov=1.45d0
rprb=2.90d0
ehomo=-0.199186d0
neleconf(1,0)=2
neleconf(2,0)=2
neleconf(2,1)=2
amu=12.011d0

case(7*1000+5)
! -----------------------           9
! N            7           5     Symbol, Z, Zion
symbol = "N"
rcov=1.42d0
rprb=2.84d0
ehomo=-0.266297d0
neleconf(2,0)=2
neleconf(2,1)=3
amu=14.00674d0

case(8*1000+6)
! -----------------------          10
! O            8           6     Symbol, Z, Zion
symbol = "O"
rcov=1.38d0
rprb=2.75d0
ehomo=-0.338381d0
neleconf(2,0)=2
neleconf(2,1)=4
amu=15.9994d0

!special case: Test PSP with all electrons
! -----------------------          10b
! O            8           8     Symbol, Z, Zion
case(8*1000+8)
symbol = "O"
rcov=1.38d0
rprb=2.75d0
ehomo=-0.338381d0
neleconf(1,0)=2
neleconf(2,0)=2
neleconf(2,1)=4
amu=15.9994d0

case(9*1000+7)
! -----------------------          11
! F            9           7     Symbol, Z, Zion
symbol = "F"
rcov=1.35d0
rprb=2.72d0
ehomo=-0.415606d0
neleconf(2,0)=2
neleconf(2,1)=5
amu=18.9984032d0

case(10*1000+8)
! -----------------------          12
! Ne          10           8     Symbol, Z, Zion
symbol = "Ne"
rcov=1.35d0
rprb=2.70d0
ehomo=-0.498034d0
neleconf(2,0)=2
neleconf(2,1)=6
amu=20.1797d0

case(11*1000+1)
! -----------------------          13
! Na          11           1     Symbol, Z, Zion
symbol = "Na"
rcov=3.40d0
rprb=6.80d0
ehomo=-0.103415d0
neleconf(3,0)=1
amu=22.989768d0

case(11*1000+9)
! -----------------------          14
! Na          11           9     Symbol, Z, Zion
symbol = "Na"
rcov=1.80d0
rprb=4.36d0
ehomo=-0.103415d0
neleconf(2,0)=2
neleconf(2,1)=6
neleconf(3,0)=1
nsccode=12
amu=22.989768d0

case(12*1000+2)
! -----------------------          16
! Mg          12           2     Symbol, Z, Zion
symbol = "Mg"
rcov=2.65d0
rprb=5.30d0
ehomo=-0.175427d0
neleconf(3,0)=2
amu=24.3050d0

case(12*1000+10)
! -----------------------          15
! Mg          12          10     Symbol, Z, Zion
symbol = "Mg"
rcov=1.20d0
rprb=3.85d0
ehomo=-0.175427d0
neleconf(2,0)=2
neleconf(2,1)=6
neleconf(3,0)=2
nsccode=12
amu=24.3050d0

case(13*1000+3)
! -----------------------          17
! Al          13           3     Symbol, Z, Zion
symbol = "Al"
rcov=2.23d0
rprb=4.45d0
ehomo=-0.102545d0
neleconf(3,0)=2
neleconf(3,1)=1
amu=26.981539d0

case(14*1000+4)
! -----------------------          18
! Si          14           4     Symbol, Z, Zion
symbol = "Si"
rcov=2.09d0
rprb=4.19d0
ehomo=-0.153293d0
neleconf(3,0)=2
neleconf(3,1)=2
amu=28.0855d0

case(14*1000+12)
! -----------------------          18
! Si          14          12     Symbol, Z, Zion
symbol = "Si"
rcov=2.09d0   !< From Si without semi-core
rprb=4.19d0   !< same
ehomo=-0.153293d0
neleconf(2,0)=2
neleconf(2,1)=6
neleconf(3,0)=2
neleconf(3,1)=2
amu=28.0855d0

case(15*1000+5)
! -----------------------          19
! P           15           5     Symbol, Z, Zion
symbol = "P"
rcov=2.00d0
rprb=4.00d0
ehomo=-0.20608d0
neleconf(3,0)=2
neleconf(3,1)=3
amu=30.973762d0

case(16*1000+6)
! -----------------------          20
! S           16           6     Symbol, Z, Zion
symbol = "S"
rcov=1.92d0
rprb=3.85d0
ehomo=-0.261676d0
neleconf(3,0)=2
neleconf(3,1)=4
amu=32.066d0

case(17*1000+7)
! -----------------------          21
! Cl          17           7     Symbol, Z, Zion
symbol = "Cl"
rcov=1.87d0
rprb=3.74d0
ehomo=-0.32038d0
neleconf(3,0)=2
neleconf(3,1)=5
amu=35.4527d0

case(18*1000+8)
! -----------------------          22
! Ar          18           8     Symbol, Z, Zion
symbol = "Ar"
rcov=1.80d0
rprb=3.60d0
ehomo=-0.38233d0
neleconf(3,0)=2
neleconf(3,1)=6
amu=39.948d0

case(19*1000+1)
! -----------------------          23
! K           19           1     Symbol, Z, Zion
symbol = "K"
rcov=4.00d0
rprb=7.00d0
ehomo=-0.088815d0
neleconf(4,0)=1
amu=39.0983d0

case(19*1000+9)
! -----------------------          24
! K           19           9     Symbol, Z, Zion
symbol = "K"
rcov=3.00d0
rprb=5.00d0
ehomo=-0.088815d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(4,0)=1
nsccode=12
amu=39.0983d0

case(20*1000+10)
! -----------------------          25
! Ca          20          10     Symbol, Z, Zion
symbol = "Ca"
rcov=3.00d0
rprb=5.00d0
ehomo=-0.141411d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(4,0)=2
nsccode=12
amu=40.078d0

case(20*1000+2)
! -----------------------          26
! Ca          20           2     Symbol, Z, Zion
symbol = "Ca"
rcov=3.80d0
rprb=7.00d0
ehomo=-0.141411d0
neleconf(4,0)=2
amu=40.078d0

case(21*1000+11)
! -----------------------          27
! Sc          21          11     Symbol, Z, Zion
symbol = "Sc"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.13108d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=1
neleconf(4,0)=2
nsccode=12
amu=44.955910d0

case(21*1000+3)
! -----------------------          28
! Sc          21           3     Symbol, Z, Zion
symbol = "Sc"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.13108d0
neleconf(3,2)=1
neleconf(4,0)=2
amu=44.955910d0

case(22*1000+12)
! -----------------------          29
! Ti          22          12     Symbol, Z, Zion
symbol = "Ti"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.167106d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=2
neleconf(4,0)=2
nsccode=12
amu=47.88d0

case(22*1000+4)
! -----------------------          30
! Ti          22           4     Symbol, Z, Zion
symbol = "Ti"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.167106d0
neleconf(3,2)=2
neleconf(4,0)=2
amu=47.88d0

case(23*1000+13)
! -----------------------          31
! V           23          13     Symbol, Z, Zion
symbol = "V"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.175968d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=3
neleconf(4,0)=2
nsccode=12
amu=50.9415d0

case(23*1000+5)
! -----------------------          32
! V           23           5     Symbol, Z, Zion
symbol = "V"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.175968d0
neleconf(3,2)=3
neleconf(4,0)=2
amu=50.9415d0

case(24*1000+14)
! -----------------------          33
! Cr          24          14     Symbol, Z, Zion
symbol = "Cr"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.118123d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=5
neleconf(4,0)=1
nsccode=12
amu=51.9961d0

case(24*1000+6)
! -----------------------          34
! Cr          24           6     Symbol, Z, Zion
symbol = "Cr"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.118123d0
neleconf(3,2)=5
neleconf(4,0)=1
amu=51.9961d0

case(25*1000+15)
! -----------------------          35
! Mn          25          15     Symbol, Z, Zion
symbol = "Mn"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.191136d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=5
neleconf(4,0)=2
nsccode=12
amu=54.93805d0

case(25*1000+7)
! -----------------------          36
! Mn          25           7     Symbol, Z, Zion
symbol = "Mn"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.191136d0
neleconf(3,2)=5
neleconf(4,0)=2
amu=54.93805d0

case(26*1000+16)
! -----------------------          37
! Fe          26          16     Symbol, Z, Zion
symbol = "Fe"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.197978d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=6
neleconf(4,0)=2
nsccode=12
amu=55.847d0

case(26*1000+8)
! -----------------------          38
! Fe          26           8     Symbol, Z, Zion
symbol = "Fe"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.197978d0
neleconf(3,2)=6
neleconf(4,0)=2
amu=55.847d0

case(27*1000+17)
! -----------------------          39
! Co          27          17     Symbol, Z, Zion
symbol = "Co"
rcov=2.40d0
rprb=4.80d0
ehomo=-0.204497d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=7
neleconf(4,0)=2
nsccode=12
amu=58.93320d0

case(27*1000+9)
! -----------------------          40
! Co          27           9     Symbol, Z, Zion
symbol = "Co"
rcov=2.40d0
rprb=4.80d0
ehomo=-0.204497d0
neleconf(3,2)=7
neleconf(4,0)=2
amu=58.93320d0

case(28*1000+10)
! -----------------------          41
! Ni          28          10     Symbol, Z, Zion
symbol = "Ni"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.210764d0
neleconf(3,2)=8
neleconf(4,0)=2
amu=58.69d0

case(28*1000+18)
! -----------------------          42
! Ni          28          18     Symbol, Z, Zion
symbol = "Ni"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.210764d0
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=8
neleconf(4,0)=2
nsccode=12
amu=58.69d0

case(29*1000+11)
! -----------------------          43
! Cu          29          11     Symbol, Z, Zion
symbol = "Cu"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.172056d0
neleconf(3,2)=10
neleconf(4,0)=1
nsccode=3
amu=63.546d0

case(29*1000+1)
! -----------------------          44
! Cu          29           1     Symbol, Z, Zion
symbol = "Cu"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.172056d0
neleconf(4,0)=1
amu=63.546d0

case(30*1000+12)
! -----------------------          45
! Zn          30          12     Symbol, Z, Zion
symbol = "Zn"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.222725d0
neleconf(3,2)=10
neleconf(4,0)=2
nsccode=3
amu=65.39d0

case(30*1000+2)
! -----------------------          46
! Zn          30           2     Symbol, Z, Zion
symbol = "Zn"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.222725d0
neleconf(4,0)=2
amu=65.39d0

case(31*1000+13)
! -----------------------          47
! Ga          31          13     Symbol, Z, Zion
symbol = "Ga"
rcov=2.10d0
rprb=4.20d0
ehomo=-0.101634d0
neleconf(3,2)=10
neleconf(4,0)=2
neleconf(4,1)=1
nsccode=3
amu=69.723d0

case(31*1000+3)
! -----------------------          48
! Ga          31           3     Symbol, Z, Zion
symbol = "Ga"
rcov=2.40d0
rprb=4.80d0
ehomo=-0.101634d0
neleconf(4,0)=2
neleconf(4,1)=1
amu=69.723d0

case(32*1000+4)
! -----------------------          49
! Ge          32           4     Symbol, Z, Zion
symbol = "Ge"
rcov=2.40d0
rprb=4.80d0
ehomo=-0.149882d0
neleconf(4,0)=2
neleconf(4,1)=2
amu=72.61d0

case(33*1000+5)
! -----------------------          50
! As          33           5     Symbol, Z, Zion
symbol = "As"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.197497d0
neleconf(4,0)=2
neleconf(4,1)=3
amu=74.92159d0

case(34*1000+6)
! -----------------------          51
! Se          34           6     Symbol, Z, Zion
symbol = "Se"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.245806d0
neleconf(4,0)=2
neleconf(4,1)=4
amu=78.96d0

case(35*1000+7)
! -----------------------          52
! Br          35           7     Symbol, Z, Zion
symbol = "Br"
rcov=2.20d0
rprb=4.40d0
ehomo=-0.295334d0
neleconf(4,0)=2
neleconf(4,1)=5
amu=79.904d0

case(36*1000+8)
! -----------------------          53
! Kr          36           8     Symbol, Z, Zion
symbol = "Kr"
rcov=2.20d0
rprb=4.40d0
ehomo=-0.34634d0
neleconf(4,0)=2
neleconf(4,1)=6
amu=83.80d0

case(37*1000+1)
! -----------------------          54
! Rb          37           1     Symbol, Z, Zion
symbol = "Rb"
rcov=4.50d0
rprb=7.00d0
ehomo=-0.085375d0
neleconf(5,0)=1
amu=85.4678d0

case(37*1000+9)
! -----------------------          55
! Rb          37           9     Symbol, Z, Zion
symbol = "Rb"
rcov=3.30d0
rprb=6.60d0
ehomo=-0.085375d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(5,0)=1
nsccode=12
amu=85.4678d0

case(38*1000+10)
! -----------------------          56
! Sr          38          10     Symbol, Z, Zion
symbol = "Sr"
rcov=3.30d0
rprb=6.60d0
ehomo=-0.131793d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(5,0)=2
nsccode=12
amu=87.62d0

case(38*1000+2)
! -----------------------          57
! Sr          38           2     Symbol, Z, Zion
symbol = "Sr"
rcov=4.00d0
rprb=7.00d0
ehomo=-0.131793d0
neleconf(5,0)=2
amu=87.62d0

case(39*1000+11)
! -----------------------          58
! Y           39          11     Symbol, Z, Zion
symbol = "Y"
rcov=3.30d0
rprb=6.60d0
ehomo=-0.108691d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=1
neleconf(5,0)=2
nsccode=12
amu=88.90585d0

case(39*1000+3)
! -----------------------          59
! Y           39           3     Symbol, Z, Zion
symbol = "Y"
rcov=3.50d0
rprb=7.00d0
ehomo=-0.108691d0
neleconf(4,2)=1
neleconf(5,0)=2
amu=88.90585d0

case(40*1000+12)
! -----------------------          60
! Zr          40          12     Symbol, Z, Zion
symbol = "Zr"
rcov=3.00d0
rprb=6.00d0
ehomo=-0.150673d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=2
neleconf(5,0)=2
nsccode=12
amu=91.224d0

case(40*1000+4)
! -----------------------          61
! Zr          40           4     Symbol, Z, Zion
symbol = "Zr"
rcov=3.00d0
rprb=6.00d0
ehomo=-0.150673d0
neleconf(4,2)=2
neleconf(5,0)=2
amu=91.224d0

case(41*1000+13)
! -----------------------          62
! Nb          41          13     Symbol, Z, Zion
symbol = "Nb"
rcov=2.92d0
rprb=5.84d0
ehomo=-0.125252d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=4
neleconf(5,0)=1
nsccode=12
amu=92.90638d0

case(41*1000+5)
! -----------------------          63
! Nb          41           5     Symbol, Z, Zion
symbol = "Nb"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.125252d0
neleconf(4,2)=4
neleconf(5,0)=1
amu=92.90638d0

case(42*1000+14)
! -----------------------          64
! Mo          42          14     Symbol, Z, Zion
symbol = "Mo"
rcov=2.83d0
rprb=5.66d0
ehomo=-0.14788d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=5
neleconf(5,0)=1
nsccode=12
amu=95.94d0

case(42*1000+6)
! -----------------------          65
! Mo          42           6     Symbol, Z, Zion
symbol = "Mo"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.14788d0
neleconf(4,2)=5
neleconf(5,0)=1
amu=95.94d0

case(43*1000+15)
! -----------------------          66
! Tc          43          15     Symbol, Z, Zion
symbol = "Tc"
rcov=2.75d0
rprb=5.50d0
ehomo=-0.183636d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=6
neleconf(5,0)=1
nsccode=12
amu=98.9062d0

case(43*1000+7)
! -----------------------          67
! Tc          43           7     Symbol, Z, Zion
symbol = "Tc"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.183636d0
neleconf(4,2)=6
neleconf(5,0)=1
amu=98.9062d0

case(44*1000+16)
! -----------------------          68
! Ru          44          16     Symbol, Z, Zion
symbol = "Ru"
rcov=2.67d0
rprb=5.34d0
ehomo=-0.152834d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=7
neleconf(5,0)=1
nsccode=12
amu=101.07d0

case(44*1000+8)
! -----------------------          69
! Ru          44           8     Symbol, Z, Zion
symbol = "Ru"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.152834d0
neleconf(4,2)=7
neleconf(5,0)=1
amu=101.07d0

case(45*1000+17)
! -----------------------          70
! Rh          45          17     Symbol, Z, Zion
symbol = "Rh"
rcov=2.58d0
rprb=5.16d0
ehomo=-0.154624d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=8
neleconf(5,0)=1
nsccode=12
amu=102.9055d0

case(45*1000+9)
! -----------------------          71
! Rh          45           9     Symbol, Z, Zion
symbol = "Rh"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.154624d0
neleconf(4,2)=8
neleconf(5,0)=1
amu=102.9055d0

case(46*1000+10)
! -----------------------          72
! Pd          46          10     Symbol, Z, Zion
symbol = "Pd"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.160771d0
neleconf(4,2)=9
neleconf(5,0)=1

case(46*1000+18)
! -----------------------          73
! Pd          46          18     Symbol, Z, Zion
symbol = "Pd"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.154624d0
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=10
nsccode=12
amu=106.42d0

case(47*1000+11)
! -----------------------          74
! Ag          47          11     Symbol, Z, Zion
symbol = "Ag"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.157407d0
neleconf(4,2)=10
neleconf(5,0)=1
nsccode=3
amu=107.8682d0

case(47*1000+1)
! -----------------------          75
! Ag          47           1     Symbol, Z, Zion
symbol = "Ag"
rcov=2.90d0
rprb=5.80d0
ehomo=-0.157407d0
neleconf(5,0)=1
amu=107.8682d0

case(48*1000+12)
! -----------------------          76
! Cd          48          12     Symbol, Z, Zion
symbol = "Cd"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.204228d0
neleconf(4,2)=10
neleconf(5,0)=2
nsccode=3
amu=112.411d0

case(48*1000+2)
! -----------------------          77
! Cd          48           2     Symbol, Z, Zion
symbol = "Cd"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.204228d0
neleconf(5,0)=2
amu=112.411d0

case(49*1000+13)
! -----------------------          78
! In          49          13     Symbol, Z, Zion
symbol = "In"
rcov=2.30d0
rprb=4.60d0
ehomo=-0.101782d0
neleconf(4,2)=10
neleconf(5,0)=2
neleconf(5,1)=1
nsccode=3
amu=114.82d0

case(49*1000+3)
! -----------------------          79
! In          49           3     Symbol, Z, Zion
symbol = "In"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.101782d0
neleconf(5,0)=2
neleconf(5,1)=1
amu=114.82d0

case(50*1000+4)
! -----------------------          80
! Sn          50           4     Symbol, Z, Zion
symbol = "Sn"
rcov=2.66d0
rprb=5.32d0
ehomo=-0.14445d0
neleconf(5,0)=2
neleconf(5,1)=2
amu=118.710d0

case(51*1000+5)
! -----------------------          81
! Sb          51           5     Symbol, Z, Zion
symbol = "Sb"
rcov=2.66d0
rprb=5.32d0
ehomo=-0.185623d0
neleconf(5,0)=2
neleconf(5,1)=3
amu=121.753d0

case(52*1000+6)
! -----------------------          82
! Te          52           6     Symbol, Z, Zion
symbol = "Te"
rcov=2.53d0
rprb=5.06d0
ehomo=-0.226594d0
neleconf(5,0)=2
neleconf(5,1)=4
amu=127.60d0

case(53*1000+7)
! -----------------------          83
! I           53           7     Symbol, Z, Zion
symbol = "I"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.267904d0
neleconf(5,0)=2
neleconf(5,1)=5
amu=126.90447d0

case(54*1000+8)
! -----------------------          84
! Xe          54           8     Symbol, Z, Zion
symbol = "Xe"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.309835d0
neleconf(5,0)=2
neleconf(5,1)=6
amu=131.29d0

case(55*1000+1)
! -----------------------          85
! Cs          55           1     Symbol, Z, Zion
symbol = "Cs"
rcov=4.50d0
rprb=7.00d0
ehomo=-0.078699d0
neleconf(6,0)=1
amu=132.90543d0

case(55*1000+9)
! -----------------------          86
! Cs          55           9     Symbol, Z, Zion
symbol = "Cs"
rcov=3.50d0
rprb=7.00d0
ehomo=-0.078699d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=1
nsccode=12
amu=132.90543d0

case(56*1000+10)
! -----------------------          87
! Ba          56          10     Symbol, Z, Zion
symbol = "Ba"
rcov=3.50d0
rprb=7.00d0
ehomo=-0.118967d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=137.327d0

case(56*1000+2)
! -----------------------          88
! Ba          56           2     Symbol, Z, Zion
symbol = "Ba"
rcov=4.00d0
rprb=7.00d0
ehomo=-0.118967d0
neleconf(6,0)=2
amu=137.327d0

case(57*1000+11)
! -----------------------          89
! La          57          11     Symbol, Z, Zion
symbol = "La"
rcov=3.50d0
rprb=7.00d0
ehomo=-0.132233d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=1
neleconf(6,0)=2
nsccode=12
amu=138.9055d0

case(58*1000+12)
! -----------------------          90
! Ce          58          12     Symbol, Z, Zion
symbol = "Ce"
rcov=3.50d0
rprb=7.00d0
ehomo=-0.133974d0
neleconf(4,3)=2
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=140.115d0

case(59*1000+13)
! -----------------------          91
! Pr          59          13     Symbol, Z, Zion
symbol = "Pr"
rcov=3.44d0
rprb=6.88d0
ehomo=-0.124465d0
neleconf(4,3)=3
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=140.90765d0

case(60*1000+14)
! -----------------------          92
! Nd          60          14     Symbol, Z, Zion
symbol = "Nd"
rcov=3.38d0
rprb=6.77d0
ehomo=-0.125796d0
neleconf(4,3)=4
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=144.24d0

case(61*1000+15)
! -----------------------          93
! Pm          61          15     Symbol, Z, Zion
symbol = "Pm"
rcov=3.33d0
rprb=6.65d0
ehomo=-0.127053d0
neleconf(4,3)=5
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=147.91d0

case(62*1000+16)
! -----------------------          94
! Sm          62          16     Symbol, Z, Zion
symbol = "Sm"
rcov=3.27d0
rprb=6.53d0
ehomo=-0.128259d0
neleconf(4,3)=6
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=150.36d0

case(63*1000+17)
! -----------------------          95
! Eu          63          17     Symbol, Z, Zion
symbol = "Eu"
rcov=3.21d0
rprb=6.42d0
ehomo=-0.129426d0
neleconf(4,3)=7
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=151.965d0

case(64*1000+18)
! -----------------------          96
! Gd          64          18     Symbol, Z, Zion
symbol = "Gd"
rcov=3.15d0
rprb=6.30d0
ehomo=-0.12722d0
neleconf(4,3)=8
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=157.25d0

case(65*1000+19)
! -----------------------          97
! Tb          65          19     Symbol, Z, Zion
symbol = "Tb"
rcov=3.09d0
rprb=6.18d0
ehomo=-0.131677d0
neleconf(4,3)=9
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=158.92534d0

case(66*1000+20)
! -----------------------          98
! Dy          66          20     Symbol, Z, Zion
symbol = "Dy"
rcov=3.03d0
rprb=6.07d0
ehomo=-0.132769d0
neleconf(4,3)=10
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=162.50d0

case(67*1000+21)
! -----------------------          99
! Ho          67          21     Symbol, Z, Zion
symbol = "Ho"
rcov=2.97d0
rprb=5.95d0
ehomo=-0.133845d0
neleconf(4,3)=11
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=164.93032d0

case(68*1000+22)
! -----------------------         100
! Er          68          22     Symbol, Z, Zion
symbol = "Er"
rcov=2.92d0
rprb=5.83d0
ehomo=-0.134905d0
neleconf(4,3)=12
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=167.26d0

case(69*1000+23)
! -----------------------         101
! Tm          69          23     Symbol, Z, Zion
symbol = "Tm"
rcov=2.92d0
rprb=5.83d0
ehomo=-0.135953d0
neleconf(4,3)=13
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=168.93421d0

case(70*1000+24)
! -----------------------         102
! Yb          70          24     Symbol, Z, Zion
symbol = "Yb"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.136989d0
neleconf(4,3)=14
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12
amu=173.04d0

case(71*1000+25)
! -----------------------         103
! Lu          71          25     Symbol, Z, Zion
symbol = "Lu"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.103686d0
neleconf(4,3)=14
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=1
neleconf(6,0)=2
nsccode=12
amu=174.967d0

case(72*1000+12)
! -----------------------         104
! Hf          72          12     Symbol, Z, Zion
symbol = "Hf"
rcov=2.90d0
rprb=5.80d0
ehomo=-0.143805d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=2
neleconf(6,0)=2
nsccode=12
amu=178.49d0

case(73*1000+13)
! -----------------------         105
! Ta          73          13     Symbol, Z, Zion
symbol = "Ta"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.174814d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=3
neleconf(6,0)=2
nsccode=12
amu=180.9479d0

case(73*1000+5)
! -----------------------         106
! Ta          73           5     Symbol, Z, Zion
symbol = "Ta"
rcov=3.10d0
rprb=6.20d0
ehomo=-0.174814d0
neleconf(5,2)=3
neleconf(6,0)=2
amu=180.9479d0

case(74*1000+14)
! -----------------------         107
! W           74          14     Symbol, Z, Zion
symbol = "W"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.181413d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=4
neleconf(6,0)=2
nsccode=12
amu=183.85d0

case(74*1000+6)
! -----------------------         108
! W           74           6     Symbol, Z, Zion
symbol = "W"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.181413d0
neleconf(5,2)=4
neleconf(6,0)=2
amu=183.85d0

case(75*1000+15)
! -----------------------         109
! Re          75          15     Symbol, Z, Zion
symbol = "Re"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.186859d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=5
neleconf(6,0)=2
nsccode=12
amu=186.207d0

case(75*1000+7)
! -----------------------         110
! Re          75           7     Symbol, Z, Zion
symbol = "Re"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.186859d0
neleconf(5,2)=5
neleconf(6,0)=2
amu=186.207d0

case(76*1000+16)
! -----------------------         111
! Os          76          16     Symbol, Z, Zion
symbol = "Os"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.191489d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=6
neleconf(6,0)=2
nsccode=12
amu=190.2d0

case(76*1000+8)
! -----------------------         112
! Os          76           8     Symbol, Z, Zion
symbol = "Os"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.191489d0
neleconf(5,2)=6
neleconf(6,0)=2
amu=190.2d0

case(77*1000+17)
! -----------------------         113
! Ir          77          17     Symbol, Z, Zion
symbol = "Ir"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.195511d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=7
neleconf(6,0)=2
nsccode=12
amu=192.22d0

case(77*1000+9)
! -----------------------         114
! Ir          77           9     Symbol, Z, Zion
symbol = "Ir"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.195511d0
neleconf(5,2)=7
neleconf(6,0)=2
amu=192.22d0

case(78*1000+10)
! -----------------------         115
! Pt          78          10     Symbol, Z, Zion
symbol = "Pt"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.161308d0
neleconf(5,2)=9
neleconf(6,0)=1
amu=195.08d0

case(78*1000+18)
! -----------------------         116
! Pt          78          18     Symbol, Z, Zion
symbol = "Pt"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.161308d0
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=8
neleconf(6,0)=1
neleconf(6,1)=1
nsccode=12
amu=195.08d0

case(79*1000+11)
! -----------------------         117
! Au          79          11     Symbol, Z, Zion
symbol = "Au"
rcov=2.70d0
rprb=5.40d0
ehomo=-0.162334d0
neleconf(5,2)=10
neleconf(6,0)=1
nsccode=3
amu=196.96654d0

case(79*1000+1)
! -----------------------         119
! Au          79           1     Symbol, Z, Zion
symbol = "Au"
rcov=4.00d0
rprb=6.40d0
ehomo=-0.162334d0
neleconf(6,0)=1
amu=196.96654d0

case(80*1000+12)
! -----------------------         120
! Hg          80          12     Symbol, Z, Zion
symbol = "Hg"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.205137d0
neleconf(5,2)=10
neleconf(6,0)=2
nsccode=3
amu=200.59d0

case(80*1000+2)
! -----------------------         121
! Hg          80           2     Symbol, Z, Zion
symbol = "Hg"
rcov=3.20d0
rprb=6.40d0
ehomo=-0.205137d0
neleconf(6,0)=2
amu=200.59d0

case(81*1000+13)
! -----------------------         122
! Tl          81          13     Symbol, Z, Zion
symbol = "Tl"
rcov=2.50d0
rprb=5.00d0
ehomo=-0.101507d0
neleconf(5,2)=10
neleconf(6,0)=2
neleconf(6,1)=1
nsccode=3
amu=204.3833d0

case(81*1000+3)
! -----------------------         123
! Tl          81           3     Symbol, Z, Zion
symbol = "Tl"
rcov=3.20d0
rprb=6.40d0
ehomo=-0.101507d0
neleconf(6,0)=2
neleconf(6,1)=1
amu=204.3833d0

case(82*1000+4)
! -----------------------         124
! Pb          82           4     Symbol, Z, Zion
symbol = "Pb"
rcov=3.30d0
rprb=6.60d0
ehomo=-0.141831d0
neleconf(6,0)=2
neleconf(6,1)=2
amu=207.2d0

case(83*1000+5)
! -----------------------         125
! Bi          83           5     Symbol, Z, Zion
symbol = "Bi"
rcov=2.90d0
rprb=5.80d0
ehomo=-0.180198d0
neleconf(6,0)=2
neleconf(6,1)=3
amu=208.98037d0

case(84*1000+6)
! -----------------------         126
! Po          84           6     Symbol, Z, Zion
symbol = "Po"
rcov=2.80d0
rprb=5.60d0
ehomo=-0.217889d0
neleconf(6,0)=2
neleconf(6,1)=4
amu=209.0d0

case(85*1000+7)
! -----------------------         127
! At          85           7     Symbol, Z, Zion
symbol = "At"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.255453d0
neleconf(6,0)=2
neleconf(6,1)=5
amu=210.0d0

case(86*1000+8)
! -----------------------         128
! Rn          86           8     Symbol, Z, Zion
symbol = "Rn"
rcov=2.60d0
rprb=5.20d0
ehomo=-0.29318d0
neleconf(6,0)=2
neleconf(6,1)=6
amu=222.0d0

case default
    write(*,*) "Electronic configuration ",nzatom,nvalelec," not found!"
    stop
end select

! Test than nvalelec is coherent with neleconf
  nsum = 0
  do l=0,lmax
     do n=1,nmax
        !write(111,*) l,n,neleconf(n,l)
        if ( neleconf(n,l) /= 0 ) nsum = nsum + neleconf(n,l)
     end do
  end do
  if (nsum /= nvalelec) then
     write(*,*) 'BUG: The electronic configuration is not correct'
     stop
  end if
  mxchg=nsum

  !correct the value of the nsccode following the new conventions
  if (nsccode /= 0) then
     sccode=real(nsccode,kind=8)
     nsccode=0
     inorbsc=ceiling(dlog(sccode)/dlog(10.d0))
     if (sccode==1.d0) inorbsc=1
     ipow=inorbsc-1
     do i=1,inorbsc
        lsc=floor(sccode/10.d0**ipow)
        sccode=sccode-real(lsc,kind=8)*10.d0**ipow
        ipow=ipow-1
        nsccode=nsccode+4**(lsc-1)
     end do
  end if

  !calculate the maximum spin polarisation  to be placed on the atom

  mxpl=0
  do l=0,lmax
     do i=1,nmax
        if (neleconf(i,l) /= 0 .and. neleconf(i,l) /= 2*(2*l+1)) then
           mxpl=mxpl+(  (2*l+1) - abs( (2*l+1)- neleconf(i,l)) ) 
        end if
     end do
  end do

END SUBROUTINE eleconf


!>   Give the symbol of element.
subroutine nzsymbol(nzatom, symbol)
  implicit none
! Arguments
  integer, intent(in) :: nzatom
  character(len=2), intent(out) :: symbol
  
 character(len=2), parameter :: symbol_(94)=(/' H','He',        &
      &   'Li','Be',' B',' C',' N',' O',' F','Ne',   &
      &   'Na','Mg','Al','Si',' P',' S','Cl','Ar',   &
      &   ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni',&
      &        'Cu','Zn','Ga','Ge','As','Se','Br','Kr',     &
      &   'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',&
      &        'Ag','Cd','In','Sn','Sb','Te',' I','Xe',     &
      &   'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
      &                       'Tb','Dy','Ho','Er','Tm','Yb',&
      &             'Lu','Hf','Ta',' W','Re','Os','Ir','Pt',&
      &        'Au','Hg','Tl','Pb','Bi','Po','At','Rn',     &
      &   'Fr','Ra','Ac','Th','Pa',' U','Np','Pu'/)

 if (nzatom <= 0 .or. nzatom > 94) then
    stop "Wrong nzatom value"
 end if
 symbol = symbol_(nzatom)
END SUBROUTINE nzsymbol


!>   Correct the electronic configuration for a given atomic charge
subroutine correct_semicore(nmax,lmax,ichg,neleconf,eleconf,nsccode)
  use module_base
  implicit none
  integer, intent(in) :: nmax,lmax,ichg
  real(kind=8) , dimension(nmax,0:lmax), intent(in) :: neleconf
  !integer, dimension(nmax,0:lmax), intent(in) :: neleconf
  real(gp), dimension(nmax,0:lmax), intent(out) :: eleconf
  integer, intent(inout) :: nsccode
  !local variables
  logical :: inocc
  integer :: i,l,nchgres,ichgp,nlsc
  real(gp) :: atchg

  !convert the array in real numbers
  do i=1,nmax
     do l=0,lmax
        eleconf(i,l)=real(neleconf(i,l),gp)
     end do
  end do

  nchgres=ichg !residual charge
  if (ichg >0) then
     !place the charge on the atom starting from the non closed shells
     do i=nmax,1,-1
        do l=lmax,0,-1
           if (neleconf(i,l) /= 2*(2*l+1) .and. neleconf(i,l) /= 0) then
              ichgp=min(nint(neleconf(i,l)),nchgres)
              nchgres=nchgres-ichgp
              eleconf(i,l)=eleconf(i,l)-real(ichgp,gp)
           end if
        end do
     end do
     if (nchgres /= 0) then
        !localise the highest occupied shell and charge it
        do i=nmax,1,-1
           do l=lmax,0,-1
              if (nint(eleconf(i,l)) == 2*(2*l+1)) then
                 ichgp=min(nint(eleconf(i,l)),nchgres)
                 nchgres=nchgres-ichgp
                 eleconf(i,l)=eleconf(i,l)-real(ichgp,gp)
              end if
           end do
        end do
!!        !charge only unoccupied shells 
!!        print *,'Atom ',symbol,': cannot charge occupied shells for the moment'
!!        stop
     end if
  else if (ichg < 0) then
     !place the charge on the atom starting from the non closed shells
     do i=nmax,1,-1
        do l=lmax,0,-1
           if (neleconf(i,l) /= 0) then
              ichgp=min(2*(2*l+1)-nint(neleconf(i,l)),-nchgres)
              nchgres=nchgres+ichgp
              eleconf(i,l)=eleconf(i,l)+real(ichgp,gp)
           end if
        end do
     end do
     if (nchgres /= 0) then
        !localise the highest unoccupied shell and charge it
        inocc=.false.
        do i=1,nmax
           do l=0,lmax
              !once found the first occupied shell search for the first unoccpied
              if (inocc .and. nint(eleconf(i,l)) == 0) then
                 ichgp=min(2*(2*l+1),-nchgres)
                 nchgres=nchgres+ichgp
                 eleconf(i,l)=eleconf(i,l)+real(ichgp,gp)
              end if
              inocc=eleconf(i,l) /= 0.0_gp
           end do
        end do
!!        !charge only occupied shells 
!!        print *,'Atom ',symbol,': cannot charge unoccupied shells for the moment'
!!        stop
     end if
     
  end if

  atchg=0.0_gp
  if (ichg /= 0) then
     !correct the semicore informations for a charged atom
     nsccode=0
     do l=0,lmax
        nlsc=0
        do i=1,nmax
           atchg=atchg+eleconf(i,l)
           if (eleconf(i,l) == real(2*(2*l+1),gp)) then
              nlsc=nlsc+1
              !if (nlsc <= 2) nsccode=nsccode+4**l
           end if
        end do
     end do
     if (atchg==0.0_gp) then
        write(*,*)'ERROR: an Atom must have input charge'
        stop
     end if
  end if

  

!!!  !if the atom has only closed shells we can treat it as semicore atom (commented)
!!!  isccode=nsccode
!!!  do l=lmax,0,-1
!!!     !control whether it is already semicore
!!!     itmp=isccode/((lmax+1)**l)
!!!     isccode=isccode-itmp*((lmax+1)**l)
!!!     !print *,'symbol',symbol,l,itmp,isccode,itmp*(lmax**l)
!!!     do i=1,nmax
!!!        if (neleconf(i,l) == 2*(2*l+1)) then
!!!           if (itmp==1) then
!!!              itmp=0
!!!              cycle
!!!           else
!!!               nsccode=nsccode+4**l !the maximum occupied is noccmax=2
!!!           end if
!!!        end if
!!!     end do
!!!  end do
END SUBROUTINE correct_semicore


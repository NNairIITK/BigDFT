!!****f* BigDFT/eleconf
!! NAME
!!   eleconf
!!
!! FUNCTION
!!   Give electronic configuration of atom
!!
!! SYNOPSIS
!!  Input
!!   nzatom    Z number of atom
!!   nvalelec  Number of valence electrons
!!  Output
!!   symbol    Atomic symbol
!!   rcov      Covalent radius
!!   rprb      Parabolic radius for the input guess using the subroutines "gatom"
!!   neleconf  Occupation number (electronic configuration of the atom)
!!   nsccode    Semicore orbitals, indicated as an integer. Each digit indicates the value(s) 
!!             of the angular momentum of the semicore orbital, increased by one
!!             e.g. if semicore are l=0 and l=2, nsccode=13
!!
!! SOURCE
!!
subroutine eleconf(nzatom,nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode)
  implicit none
! Arguments
  integer, intent(in) :: nzatom,nvalelec
  character(len=2), intent(out) :: symbol
  real(kind=8), intent(out) :: rcov,rprb,ehomo
  integer, parameter :: nmax=6,lmax=3
  integer, intent(out) :: neleconf(nmax,0:lmax)
  integer, intent(out) :: nsccode
! Local variables
  integer :: n,l,nsum

  neleconf(:,:)=0
  nsccode=0

!Eah atomic configuration
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
rcov=0.75
rprb=1.21
ehomo=-0.233471
neleconf(1,0)=1

case(2*1000+2)
! -----------------------           2
! He           2           2     Symbol, Z, Zion
symbol = "He"
rcov=0.75
rprb=1.50
ehomo=-0.570425
neleconf(1,0)=2

case(3*1000+1)
! -----------------------           3
! Li           3           1     Symbol, Z, Zion
symbol = "Li"
rcov=3.40
rprb=6.80
ehomo=-0.10554
neleconf(2,0)=1

case(3*1000+3)
! -----------------------           4
! Li           3           3     Symbol, Z, Zion
symbol = "Li"
rcov=1.60
rprb=3.61
ehomo=-0.10554
neleconf(1,0)=2
neleconf(2,0)=1
nsccode=1

case(4*1000+2)
! -----------------------           5
! Be           4           2     Symbol, Z, Zion
symbol = "Be"
rcov=2.30
rprb=4.60
ehomo=-0.205744
neleconf(2,0)=2

case(4*1000+4)
! -----------------------           6
! Be           4           4     Symbol, Z, Zion
symbol = "Be"
rcov=1.30
rprb=3.60
ehomo=-0.205744
neleconf(1,0)=2
neleconf(2,0)=2
nsccode=1

case(5*1000+3)
! -----------------------           7
! B            5           3     Symbol, Z, Zion
symbol = "B"
rcov=1.55
rprb=3.10
ehomo=-0.136603
neleconf(2,0)=2
neleconf(2,1)=1

case(6*1000+4)
! -----------------------           8
! C            6           4     Symbol, Z, Zion
symbol = "C"
rcov=1.45
rprb=2.90
ehomo=-0.199186
neleconf(2,0)=2
neleconf(2,1)=2

case(7*1000+5)
! -----------------------           9
! N            7           5     Symbol, Z, Zion
symbol = "N"
rcov=1.42
rprb=2.84
ehomo=-0.266297
neleconf(2,0)=2
neleconf(2,1)=3

case(8*1000+6)
! -----------------------          10
! O            8           6     Symbol, Z, Zion
symbol = "O"
rcov=1.38
rprb=2.75
ehomo=-0.338381
neleconf(2,0)=2
neleconf(2,1)=4

case(9*1000+7)
! -----------------------          11
! F            9           7     Symbol, Z, Zion
symbol = "F"
rcov=1.35
rprb=2.72
ehomo=-0.415606
neleconf(2,0)=2
neleconf(2,1)=5

case(10*1000+8)
! -----------------------          12
! Ne          10           8     Symbol, Z, Zion
symbol = "Ne"
rcov=1.35
rprb=2.70
ehomo=-0.498034
neleconf(2,0)=2
neleconf(2,1)=6

case(11*1000+1)
! -----------------------          13
! Na          11           1     Symbol, Z, Zion
symbol = "Na"
rcov=3.40
rprb=6.80
ehomo=-0.103415
neleconf(3,0)=1

case(11*1000+9)
! -----------------------          14
! Na          11           9     Symbol, Z, Zion
symbol = "Na"
rcov=1.80
rprb=4.36
ehomo=-0.103415
neleconf(2,0)=2
neleconf(2,1)=6
neleconf(3,0)=1
nsccode=12

case(12*1000+10)
! -----------------------          15
! Mg          12          10     Symbol, Z, Zion
symbol = "Mg"
rcov=1.20
rprb=3.85
ehomo=-0.175427
neleconf(2,0)=2
neleconf(2,1)=6
neleconf(3,0)=2
nsccode=12

case(12*1000+2)
! -----------------------          16
! Mg          12           2     Symbol, Z, Zion
symbol = "Mg"
rcov=2.65
rprb=5.30
ehomo=-0.175427
neleconf(3,0)=2

case(13*1000+3)
! -----------------------          17
! Al          13           3     Symbol, Z, Zion
symbol = "Al"
rcov=2.23
rprb=4.45
ehomo=-0.102545
neleconf(3,0)=2
neleconf(3,1)=1

case(14*1000+4)
! -----------------------          18
! Si          14           4     Symbol, Z, Zion
symbol = "Si"
rcov=2.09
rprb=4.19
ehomo=-0.153293
neleconf(3,0)=2
neleconf(3,1)=2

case(15*1000+5)
! -----------------------          19
! P           15           5     Symbol, Z, Zion
symbol = "P"
rcov=2.00
rprb=4.00
ehomo=-0.20608
neleconf(3,0)=2
neleconf(3,1)=3

case(16*1000+6)
! -----------------------          20
! S           16           6     Symbol, Z, Zion
symbol = "S"
rcov=1.92
rprb=3.85
ehomo=-0.261676
neleconf(3,0)=2
neleconf(3,1)=4

case(17*1000+7)
! -----------------------          21
! Cl          17           7     Symbol, Z, Zion
symbol = "Cl"
rcov=1.87
rprb=3.74
ehomo=-0.32038
neleconf(3,0)=2
neleconf(3,1)=5

case(18*1000+8)
! -----------------------          22
! Ar          18           8     Symbol, Z, Zion
symbol = "Ar"
rcov=1.80
rprb=3.60
ehomo=-0.38233
neleconf(3,0)=2
neleconf(3,1)=6

case(19*1000+1)
! -----------------------          23
! K           19           1     Symbol, Z, Zion
symbol = "K"
rcov=4.00
rprb=7.00
ehomo=-0.088815
neleconf(4,0)=1

case(19*1000+9)
! -----------------------          24
! K           19           9     Symbol, Z, Zion
symbol = "K"
rcov=3.00
rprb=5.00
ehomo=-0.088815
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(4,0)=1
nsccode=12

case(20*1000+10)
! -----------------------          25
! Ca          20          10     Symbol, Z, Zion
symbol = "Ca"
rcov=3.00
rprb=5.00
ehomo=-0.141411
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(4,0)=2
nsccode=12

case(20*1000+2)
! -----------------------          26
! Ca          20           2     Symbol, Z, Zion
symbol = "Ca"
rcov=3.80
rprb=7.00
ehomo=-0.141411
neleconf(4,0)=2

case(21*1000+11)
! -----------------------          27
! Sc          21          11     Symbol, Z, Zion
symbol = "Sc"
rcov=2.70
rprb=5.40
ehomo=-0.13108
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=1
neleconf(4,0)=2
nsccode=12

case(21*1000+3)
! -----------------------          28
! Sc          21           3     Symbol, Z, Zion
symbol = "Sc"
rcov=2.70
rprb=5.40
ehomo=-0.13108
neleconf(3,2)=1
neleconf(4,0)=2

case(22*1000+12)
! -----------------------          29
! Ti          22          12     Symbol, Z, Zion
symbol = "Ti"
rcov=2.70
rprb=5.40
ehomo=-0.167106
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=2
neleconf(4,0)=2
nsccode=12

case(22*1000+4)
! -----------------------          30
! Ti          22           4     Symbol, Z, Zion
symbol = "Ti"
rcov=2.70
rprb=5.40
ehomo=-0.167106
neleconf(3,2)=2
neleconf(4,0)=2

case(23*1000+13)
! -----------------------          31
! V           23          13     Symbol, Z, Zion
symbol = "V"
rcov=2.60
rprb=5.20
ehomo=-0.175968
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=3
neleconf(4,0)=2
nsccode=12

case(23*1000+5)
! -----------------------          32
! V           23           5     Symbol, Z, Zion
symbol = "V"
rcov=2.60
rprb=5.20
ehomo=-0.175968
neleconf(3,2)=3
neleconf(4,0)=2

case(24*1000+14)
! -----------------------          33
! Cr          24          14     Symbol, Z, Zion
symbol = "Cr"
rcov=2.60
rprb=5.20
ehomo=-0.118123
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=5
neleconf(4,0)=1
nsccode=12

case(24*1000+6)
! -----------------------          34
! Cr          24           6     Symbol, Z, Zion
symbol = "Cr"
rcov=2.60
rprb=5.20
ehomo=-0.118123
neleconf(3,2)=5
neleconf(4,0)=1

case(25*1000+15)
! -----------------------          35
! Mn          25          15     Symbol, Z, Zion
symbol = "Mn"
rcov=2.50
rprb=5.00
ehomo=-0.191136
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=5
neleconf(4,0)=2
nsccode=12

case(25*1000+7)
! -----------------------          36
! Mn          25           7     Symbol, Z, Zion
symbol = "Mn"
rcov=2.50
rprb=5.00
ehomo=-0.191136
neleconf(3,2)=5
neleconf(4,0)=2

case(26*1000+16)
! -----------------------          37
! Fe          26          16     Symbol, Z, Zion
symbol = "Fe"
rcov=2.50
rprb=5.00
ehomo=-0.197978
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=6
neleconf(4,0)=2
nsccode=12

case(26*1000+8)
! -----------------------          38
! Fe          26           8     Symbol, Z, Zion
symbol = "Fe"
rcov=2.50
rprb=5.00
ehomo=-0.197978
neleconf(3,2)=6
neleconf(4,0)=2

case(27*1000+17)
! -----------------------          39
! Co          27          17     Symbol, Z, Zion
symbol = "Co"
rcov=2.40
rprb=4.80
ehomo=-0.204497
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=7
neleconf(4,0)=2
nsccode=12

case(27*1000+9)
! -----------------------          40
! Co          27           9     Symbol, Z, Zion
symbol = "Co"
rcov=2.40
rprb=4.80
ehomo=-0.204497
neleconf(3,2)=7
neleconf(4,0)=2

case(28*1000+10)
! -----------------------          41
! Ni          28          10     Symbol, Z, Zion
symbol = "Ni"
rcov=2.30
rprb=4.60
ehomo=-0.210764
neleconf(3,2)=8
neleconf(4,0)=2

case(28*1000+18)
! -----------------------          42
! Ni          28          18     Symbol, Z, Zion
symbol = "Ni"
rcov=2.30
rprb=4.60
ehomo=-0.210764
neleconf(3,0)=2
neleconf(3,1)=6
neleconf(3,2)=8
neleconf(4,0)=2
nsccode=12

case(29*1000+11)
! -----------------------          43
! Cu          29          11     Symbol, Z, Zion
symbol = "Cu"
rcov=2.30
rprb=4.60
ehomo=-0.172056
neleconf(3,2)=10
neleconf(4,0)=1
nsccode=3

case(29*1000+1)
! -----------------------          44
! Cu          29           1     Symbol, Z, Zion
symbol = "Cu"
rcov=2.80
rprb=5.60
ehomo=-0.172056
neleconf(4,0)=1

case(30*1000+12)
! -----------------------          45
! Zn          30          12     Symbol, Z, Zion
symbol = "Zn"
rcov=2.30
rprb=4.60
ehomo=-0.222725
neleconf(3,2)=10
neleconf(4,0)=2
nsccode=3

case(30*1000+2)
! -----------------------          46
! Zn          30           2     Symbol, Z, Zion
symbol = "Zn"
rcov=2.70
rprb=5.40
ehomo=-0.222725
neleconf(4,0)=2

case(31*1000+13)
! -----------------------          47
! Ga          31          13     Symbol, Z, Zion
symbol = "Ga"
rcov=2.10
rprb=4.20
ehomo=-0.101634
neleconf(3,2)=10
neleconf(4,0)=2
neleconf(4,1)=1
nsccode=3

case(31*1000+3)
! -----------------------          48
! Ga          31           3     Symbol, Z, Zion
symbol = "Ga"
rcov=2.40
rprb=4.80
ehomo=-0.101634
neleconf(4,0)=2
neleconf(4,1)=1

case(32*1000+4)
! -----------------------          49
! Ge          32           4     Symbol, Z, Zion
symbol = "Ge"
rcov=2.40
rprb=4.80
ehomo=-0.149882
neleconf(4,0)=2
neleconf(4,1)=2

case(33*1000+5)
! -----------------------          50
! As          33           5     Symbol, Z, Zion
symbol = "As"
rcov=2.30
rprb=4.60
ehomo=-0.197497
neleconf(4,0)=2
neleconf(4,1)=3

case(34*1000+6)
! -----------------------          51
! Se          34           6     Symbol, Z, Zion
symbol = "Se"
rcov=2.30
rprb=4.60
ehomo=-0.245806
neleconf(4,0)=2
neleconf(4,1)=4

case(35*1000+7)
! -----------------------          52
! Br          35           7     Symbol, Z, Zion
symbol = "Br"
rcov=2.20
rprb=4.40
ehomo=-0.295334
neleconf(4,0)=2
neleconf(4,1)=5

case(36*1000+8)
! -----------------------          53
! Kr          36           8     Symbol, Z, Zion
symbol = "Kr"
rcov=2.20
rprb=4.40
ehomo=-0.34634
neleconf(4,0)=2
neleconf(4,1)=6

case(37*1000+1)
! -----------------------          54
! Rb          37           1     Symbol, Z, Zion
symbol = "Rb"
rcov=4.50
rprb=7.00
ehomo=-0.085375
neleconf(5,0)=1

case(37*1000+9)
! -----------------------          55
! Rb          37           9     Symbol, Z, Zion
symbol = "Rb"
rcov=3.30
rprb=6.60
ehomo=-0.085375
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(5,0)=1
nsccode=12

case(38*1000+10)
! -----------------------          56
! Sr          38          10     Symbol, Z, Zion
symbol = "Sr"
rcov=3.30
rprb=6.60
ehomo=-0.131793
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(5,0)=2
nsccode=12

case(38*1000+2)
! -----------------------          57
! Sr          38           2     Symbol, Z, Zion
symbol = "Sr"
rcov=4.00
rprb=7.00
ehomo=-0.131793
neleconf(5,0)=2

case(39*1000+11)
! -----------------------          58
! Y           39          11     Symbol, Z, Zion
symbol = "Y"
rcov=3.30
rprb=6.60
ehomo=-0.108691
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=1
neleconf(5,0)=2
nsccode=12

case(39*1000+3)
! -----------------------          59
! Y           39           3     Symbol, Z, Zion
symbol = "Y"
rcov=3.50
rprb=7.00
ehomo=-0.108691
neleconf(4,2)=1
neleconf(5,0)=2

case(40*1000+12)
! -----------------------          60
! Zr          40          12     Symbol, Z, Zion
symbol = "Zr"
rcov=3.00
rprb=6.00
ehomo=-0.150673
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=2
neleconf(5,0)=2
nsccode=12

case(40*1000+4)
! -----------------------          61
! Zr          40           4     Symbol, Z, Zion
symbol = "Zr"
rcov=3.00
rprb=6.00
ehomo=-0.150673
neleconf(4,2)=2
neleconf(5,0)=2

case(41*1000+13)
! -----------------------          62
! Nb          41          13     Symbol, Z, Zion
symbol = "Nb"
rcov=2.92
rprb=5.84
ehomo=-0.125252
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=4
neleconf(5,0)=1
nsccode=12

case(41*1000+5)
! -----------------------          63
! Nb          41           5     Symbol, Z, Zion
symbol = "Nb"
rcov=2.70
rprb=5.40
ehomo=-0.125252
neleconf(4,2)=4
neleconf(5,0)=1

case(42*1000+14)
! -----------------------          64
! Mo          42          14     Symbol, Z, Zion
symbol = "Mo"
rcov=2.83
rprb=5.66
ehomo=-0.14788
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=5
neleconf(5,0)=1
nsccode=12

case(42*1000+6)
! -----------------------          65
! Mo          42           6     Symbol, Z, Zion
symbol = "Mo"
rcov=2.60
rprb=5.20
ehomo=-0.14788
neleconf(4,2)=5
neleconf(5,0)=1
nsccode=12

case(43*1000+15)
! -----------------------          66
! Tc          43          15     Symbol, Z, Zion
symbol = "Tc"
rcov=2.75
rprb=5.50
ehomo=-0.183636
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=6
neleconf(5,0)=1
nsccode=12

case(43*1000+7)
! -----------------------          67
! Tc          43           7     Symbol, Z, Zion
symbol = "Tc"
rcov=2.60
rprb=5.20
ehomo=-0.183636
neleconf(4,2)=6
neleconf(5,0)=1

case(44*1000+16)
! -----------------------          68
! Ru          44          16     Symbol, Z, Zion
symbol = "Ru"
rcov=2.67
rprb=5.34
ehomo=-0.152834
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=7
neleconf(5,0)=1
nsccode=12

case(44*1000+8)
! -----------------------          69
! Ru          44           8     Symbol, Z, Zion
symbol = "Ru"
rcov=2.50
rprb=5.00
ehomo=-0.152834
neleconf(4,2)=7
neleconf(5,0)=1

case(45*1000+17)
! -----------------------          70
! Rh          45          17     Symbol, Z, Zion
symbol = "Rh"
rcov=2.58
rprb=5.16
ehomo=-0.154624
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=8
neleconf(5,0)=1
nsccode=12

case(45*1000+9)
! -----------------------          71
! Rh          45           9     Symbol, Z, Zion
symbol = "Rh"
rcov=2.50
rprb=5.00
ehomo=-0.154624
neleconf(4,2)=8
neleconf(5,0)=1

case(46*1000+10)
! -----------------------          72
! Pd          46          10     Symbol, Z, Zion
symbol = "Pd"
rcov=2.50
rprb=5.00
ehomo=-0.160771
neleconf(4,2)=10

case(46*1000+18)
! -----------------------          73
! Pd          46          18     Symbol, Z, Zion
symbol = "Rh"
rcov=2.50
rprb=5.00
ehomo=-0.154624
neleconf(4,0)=2
neleconf(4,1)=6
neleconf(4,2)=10
nsccode=12

case(47*1000+11)
! -----------------------          74
! Ag          47          11     Symbol, Z, Zion
symbol = "Ag"
rcov=2.50
rprb=5.00
ehomo=-0.157407
neleconf(4,2)=10
neleconf(5,0)=1
nsccode=3

case(47*1000+1)
! -----------------------          75
! Ag          47           1     Symbol, Z, Zion
symbol = "Ag"
rcov=2.90
rprb=5.80
ehomo=-0.157407
neleconf(5,0)=1

case(48*1000+12)
! -----------------------          76
! Cd          48          12     Symbol, Z, Zion
symbol = "Cd"
rcov=2.50
rprb=5.00
ehomo=-0.204228
neleconf(4,2)=10
neleconf(5,0)=2
nsccode=3

case(48*1000+2)
! -----------------------          77
! Cd          48           2     Symbol, Z, Zion
symbol = "Cd"
rcov=2.80
rprb=5.60
ehomo=-0.204228
neleconf(5,0)=2

case(49*1000+13)
! -----------------------          78
! In          49          13     Symbol, Z, Zion
symbol = "In"
rcov=2.30
rprb=4.60
ehomo=-0.101782
neleconf(4,2)=10
neleconf(5,0)=2
neleconf(5,1)=1
nsccode=3

case(49*1000+3)
! -----------------------          79
! In          49           3     Symbol, Z, Zion
symbol = "In"
rcov=2.70
rprb=5.40
ehomo=-0.101782
neleconf(5,0)=2
neleconf(5,1)=1

case(50*1000+4)
! -----------------------          80
! Sn          50           4     Symbol, Z, Zion
symbol = "Sn"
rcov=2.66
rprb=5.32
ehomo=-0.14445
neleconf(5,0)=2
neleconf(5,1)=2

case(51*1000+5)
! -----------------------          81
! Sb          51           5     Symbol, Z, Zion
symbol = "Sb"
rcov=2.66
rprb=5.32
ehomo=-0.185623
neleconf(5,0)=2
neleconf(5,1)=3

case(52*1000+6)
! -----------------------          82
! Te          52           6     Symbol, Z, Zion
symbol = "Te"
rcov=2.53
rprb=5.06
ehomo=-0.226594
neleconf(5,0)=2
neleconf(5,1)=4

case(53*1000+7)
! -----------------------          83
! I           53           7     Symbol, Z, Zion
symbol = "I"
rcov=2.50
rprb=5.00
ehomo=-0.267904
neleconf(5,0)=2
neleconf(5,1)=5

case(54*1000+8)
! -----------------------          84
! Xe          54           8     Symbol, Z, Zion
symbol = "Xe"
rcov=2.50
rprb=5.00
ehomo=-0.309835
neleconf(5,0)=2
neleconf(5,1)=6

case(55*1000+1)
! -----------------------          85
! Cs          55           1     Symbol, Z, Zion
symbol = "Cs"
rcov=4.50
rprb=7.00
ehomo=-0.078699
neleconf(6,0)=1

case(55*1000+9)
! -----------------------          86
! Cs          55           9     Symbol, Z, Zion
symbol = "Cs"
rcov=3.50
rprb=7.00
ehomo=-0.078699
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=1
nsccode=12

case(56*1000+10)
! -----------------------          87
! Ba          56          10     Symbol, Z, Zion
symbol = "Ba"
rcov=3.50
rprb=7.00
ehomo=-0.118967
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(56*1000+2)
! -----------------------          88
! Ba          56           2     Symbol, Z, Zion
symbol = "Ba"
rcov=4.00
rprb=7.00
ehomo=-0.118967
neleconf(6,0)=2

case(57*1000+11)
! -----------------------          89
! La          57          11     Symbol, Z, Zion
symbol = "La"
rcov=3.50
rprb=7.00
ehomo=-0.132233
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=1
neleconf(6,0)=2
nsccode=12

case(58*1000+12)
! -----------------------          90
! Ce          58          12     Symbol, Z, Zion
symbol = "Ce"
rcov=3.50
rprb=7.00
ehomo=-0.133974
neleconf(4,3)=2
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(59*1000+13)
! -----------------------          91
! Pr          59          13     Symbol, Z, Zion
symbol = "Pr"
rcov=3.44
rprb=6.88
ehomo=-0.124465
neleconf(4,3)=3
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(60*1000+14)
! -----------------------          92
! Nd          60          14     Symbol, Z, Zion
symbol = "Nd"
rcov=3.38
rprb=6.77
ehomo=-0.125796
neleconf(4,3)=4
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(61*1000+15)
! -----------------------          93
! Pm          61          15     Symbol, Z, Zion
symbol = "Pm"
rcov=3.33
rprb=6.65
ehomo=-0.127053
neleconf(4,3)=5
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(62*1000+16)
! -----------------------          94
! Sm          62          16     Symbol, Z, Zion
symbol = "Sm"
rcov=3.27
rprb=6.53
ehomo=-0.128259
neleconf(4,3)=6
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(63*1000+17)
! -----------------------          95
! Eu          63          17     Symbol, Z, Zion
symbol = "Eu"
rcov=3.21
rprb=6.42
ehomo=-0.129426
neleconf(4,3)=7
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(64*1000+18)
! -----------------------          96
! Gd          64          18     Symbol, Z, Zion
symbol = "Gd"
rcov=3.15
rprb=6.30
ehomo=-0.12722
neleconf(4,3)=8
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(65*1000+19)
! -----------------------          97
! Tb          65          19     Symbol, Z, Zion
symbol = "Tb"
rcov=3.09
rprb=6.18
ehomo=-0.131677
neleconf(4,3)=9
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(66*1000+20)
! -----------------------          98
! Dy          66          20     Symbol, Z, Zion
symbol = "Dy"
rcov=3.03
rprb=6.07
ehomo=-0.132769
neleconf(4,3)=10
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(67*1000+21)
! -----------------------          99
! Ho          67          21     Symbol, Z, Zion
symbol = "Ho"
rcov=2.97
rprb=5.95
ehomo=-0.133845
neleconf(4,3)=11
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(68*1000+22)
! -----------------------         100
! Er          68          22     Symbol, Z, Zion
symbol = "Er"
rcov=2.92
rprb=5.83
ehomo=-0.134905
neleconf(4,3)=12
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(69*1000+23)
! -----------------------         101
! Tm          69          23     Symbol, Z, Zion
symbol = "Tm"
rcov=2.92
rprb=5.83
ehomo=-0.135953
neleconf(4,3)=13
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(70*1000+24)
! -----------------------         102
! Yb          70          24     Symbol, Z, Zion
symbol = "Yb"
rcov=2.80
rprb=5.60
ehomo=-0.136989
neleconf(4,3)=14
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(6,0)=2
nsccode=12

case(71*1000+25)
! -----------------------         103
! Lu          71          25     Symbol, Z, Zion
symbol = "Lu"
rcov=2.80
rprb=5.60
ehomo=-0.103686
neleconf(4,3)=14
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=1
neleconf(6,0)=2
nsccode=12

case(72*1000+12)
! -----------------------         104
! Hf          72          12     Symbol, Z, Zion
symbol = "Hf"
rcov=2.90
rprb=5.80
ehomo=-0.143805
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=2
neleconf(6,0)=2
nsccode=12

case(73*1000+13)
! -----------------------         105
! Ta          73          13     Symbol, Z, Zion
symbol = "Ta"
rcov=2.70
rprb=5.40
ehomo=-0.174814
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=3
neleconf(6,0)=2
nsccode=12

case(73*1000+5)
! -----------------------         106
! Ta          73           5     Symbol, Z, Zion
symbol = "Ta"
rcov=3.10
rprb=6.20
ehomo=-0.174814
neleconf(5,2)=3
neleconf(6,0)=2

case(74*1000+14)
! -----------------------         107
! W           74          14     Symbol, Z, Zion
symbol = "W"
rcov=2.60
rprb=5.20
ehomo=-0.181413
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=4
neleconf(6,0)=2
nsccode=12

case(74*1000+6)
! -----------------------         108
! W           74           6     Symbol, Z, Zion
symbol = "W"
rcov=2.60
rprb=5.20
ehomo=-0.181413
neleconf(5,2)=4
neleconf(6,0)=2

case(75*1000+15)
! -----------------------         109
! Re          75          15     Symbol, Z, Zion
symbol = "Re"
rcov=2.60
rprb=5.20
ehomo=-0.186859
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=5
neleconf(6,0)=2
nsccode=12

case(75*1000+7)
! -----------------------         110
! Re          75           7     Symbol, Z, Zion
symbol = "Re"
rcov=2.60
rprb=5.20
ehomo=-0.186859
neleconf(5,2)=5
neleconf(6,0)=2

case(76*1000+16)
! -----------------------         111
! Os          76          16     Symbol, Z, Zion
symbol = "Os"
rcov=2.50
rprb=5.00
ehomo=-0.191489
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=6
neleconf(6,0)=2
nsccode=12

case(76*1000+8)
! -----------------------         112
! Os          76           8     Symbol, Z, Zion
symbol = "Os"
rcov=2.50
rprb=5.00
ehomo=-0.191489
neleconf(5,2)=6
neleconf(6,0)=2

case(77*1000+17)
! -----------------------         113
! Ir          77          17     Symbol, Z, Zion
symbol = "Ir"
rcov=2.50
rprb=5.00
ehomo=-0.195511
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=7
neleconf(6,0)=2
nsccode=12

case(77*1000+9)
! -----------------------         114
! Ir          77           9     Symbol, Z, Zion
symbol = "Ir"
rcov=2.60
rprb=5.20
ehomo=-0.195511
neleconf(5,2)=7
neleconf(6,0)=2

case(78*1000+10)
! -----------------------         115
! Pt          78          10     Symbol, Z, Zion
symbol = "Pt"
rcov=2.60
rprb=5.20
ehomo=-0.161308
neleconf(5,2)=9
neleconf(6,0)=1

case(78*1000+18)
! -----------------------         116
! Pt          78          18     Symbol, Z, Zion
symbol = "Pt"
rcov=2.60
rprb=5.20
ehomo=-0.161308
neleconf(5,0)=2
neleconf(5,1)=6
neleconf(5,2)=9
neleconf(6,0)=1
nsccode=12

case(79*1000+11)
! -----------------------         117
! Au          79          11     Symbol, Z, Zion
symbol = "Au"
rcov=2.70
rprb=5.40
ehomo=-0.162334
neleconf(5,2)=10
neleconf(6,0)=1
nsccode=3

case(79*1000+1)
! -----------------------         119
! Au          79           1     Symbol, Z, Zion
symbol = "Au"
rcov=4.00
rprb=6.40
ehomo=-0.162334
neleconf(6,0)=1

case(80*1000+12)
! -----------------------         120
! Hg          80          12     Symbol, Z, Zion
symbol = "Hg"
rcov=2.80
rprb=5.60
ehomo=-0.205137
neleconf(5,2)=10
neleconf(6,0)=2
nsccode=3

case(80*1000+2)
! -----------------------         121
! Hg          80           2     Symbol, Z, Zion
symbol = "Hg"
rcov=3.20
rprb=6.40
ehomo=-0.205137
neleconf(6,0)=2

case(81*1000+13)
! -----------------------         122
! Tl          81          13     Symbol, Z, Zion
symbol = "Tl"
rcov=2.50
rprb=5.00
ehomo=-0.101507
neleconf(5,2)=10
neleconf(6,0)=2
neleconf(6,1)=1
nsccode=3

case(81*1000+3)
! -----------------------         123
! Tl          81           3     Symbol, Z, Zion
symbol = "Tl"
rcov=3.20
rprb=6.40
ehomo=-0.101507
neleconf(6,0)=2
neleconf(6,1)=1

case(82*1000+4)
! -----------------------         124
! Pb          82           4     Symbol, Z, Zion
symbol = "Pb"
rcov=3.30
rprb=6.60
ehomo=-0.141831
neleconf(6,0)=2
neleconf(6,1)=2

case(83*1000+5)
! -----------------------         125
! Bi          83           5     Symbol, Z, Zion
symbol = "Bi"
rcov=2.90
rprb=5.80
ehomo=-0.180198
neleconf(6,0)=2
neleconf(6,1)=3

case(84*1000+6)
! -----------------------         126
! Po          84           6     Symbol, Z, Zion
symbol = "Po"
rcov=2.80
rprb=5.60
ehomo=-0.217889
neleconf(6,0)=2
neleconf(6,1)=4

case(85*1000+7)
! -----------------------         127
! At          85           7     Symbol, Z, Zion
symbol = "At"
rcov=2.60
rprb=5.20
ehomo=-0.255453
neleconf(6,0)=2
neleconf(6,1)=5

case(86*1000+8)
! -----------------------         128
! Rn          86           8     Symbol, Z, Zion
symbol = "Rn"
rcov=2.60
rprb=5.20
ehomo=-0.29318
neleconf(6,0)=2
neleconf(6,1)=6

case default
    write(*,*) "Electronic configuration ",nzatom,nvalelec," not found!"
    stop
end select

! Test than nvalelec is coherent with neleconf
  nsum = 0
  do n=1,nmax
     do l=0,lmax
        nsum = nsum + neleconf(n,l)
     end do
  end do
  if (nsum /= nvalelec) then
     write(*,*) 'BUG: The electronic configuration is not correct'
     stop
  end if

end subroutine eleconf
!!***

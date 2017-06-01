!{\src2tex{textfont=tt}}
!!****m* ABINIT/abi_defs_basis
!! NAME
!! abi_defs_basis
!!
!! FUNCTION
!! This module contains definitions for a number of named constants
!! and physical constants
!!
!! COPYRIGHT
!! Copyright (C) 2000-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Of the named constants,
!! by far the most important are those that define the 'kind' types of
!! virtually all the variables used in a (well-written) FORTRAN 90 code
!! the content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module abi_defs_basis

 implicit none

!Keyword 'integer' stands for default integer type
!and may be used whenever integer are presumed to be small

!nb of bytes related to an integer subtype n such as -10^(argument) < n < 10^(argument) (this is standard F90)
 integer, parameter :: i1b=selected_int_kind(2)
 integer, parameter :: i2b=selected_int_kind(4)
 integer, parameter :: i4b=selected_int_kind(9)
 integer, parameter :: i8b=selected_int_kind(18)

!nb of bytes related to default simple-precision real/complex subtypes
!(= 4 for many machine architectures, = 8 for e.g. Cray)
 integer, parameter :: sp=kind(1.0)
 integer, parameter :: spc=kind((1.0,1.0))

!nb of bytes related to default double-precision real/complex subtypes
!(= 8 for many machine architectures)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))

!The default lengths
 integer, parameter :: fnlen=264      ! maximum length of file name variables
 integer, parameter :: strlen=2000000 ! maximum length of input string

!UNIX unit numbers : standard input, standard output, ab_out,
!!                   and a number for temporary access to a file.
 integer, parameter :: std_in=5,ab_in=5
 integer, parameter :: std_out_default=6,ab_out_default=7
 integer, parameter :: std_err=0
 integer, parameter :: dev_null=-1       ! Fake unit number used to skip the printing in abi_wrtout.
 integer, parameter :: ab_xml_out = 50   ! this unit is used to print output into an XML file
 integer, parameter :: tmp_unit=9,tmp_unit2=10

! These vars should be private and only modifiable via an appropriate method (see below)
 integer, public, save :: ab_out = ab_out_default
 integer, public, save :: std_out =  std_out_default
!It should be put to xmpi_world (but it is not possible for the moment - v6.9)
 integer, public, save :: abinit_comm_output = -1 !This default value has to be changed at start !!!

! the maximum length of a record in a file connected for sequential access.
 integer,public,parameter :: ABI_RECL=524288  ! 2**19

!Real constants
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: five=5._dp
 real(dp), parameter :: six=6._dp
 real(dp), parameter :: seven=7._dp
 real(dp), parameter :: eight=8._dp
 real(dp), parameter :: nine=9._dp
 real(dp), parameter :: ten=10._dp

!Fractionary real constants
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: onehalf=1.50_dp
 real(dp), parameter :: third=one/three
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: fifth=0.20_dp
 real(dp), parameter :: sixth=one/six
 real(dp), parameter :: seventh=one/seven
 real(dp), parameter :: eighth=0.125_dp
 real(dp), parameter :: ninth=one/nine
 real(dp), parameter :: two_thirds=two*third
 real(dp), parameter :: four_thirds=four*third
 real(dp), parameter :: five_thirds=five*third
 real(dp), parameter :: three_quarters=0.75_dp
 real(dp), parameter :: three_fifth=three/five

!Real constants related to the golden number
 real(dp), parameter :: gold=1.618033988749894848204586834365638117720309179_dp
 real(dp), parameter :: goldenratio=two-gold

!Real constants derived from pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=two*pi
 real(dp), parameter :: four_pi=four*pi
 real(dp), parameter :: piinv=one/pi

!Real precision
 real(dp), parameter :: greatest_real = huge(one)
 real(dp), parameter :: smallest_real = -greatest_real
 real(dp), parameter :: tol3= 0.001_dp
 real(dp), parameter :: tol4= 0.0001_dp
 real(dp), parameter :: tol5= 0.00001_dp
 real(dp), parameter :: tol6= 0.000001_dp
 real(dp), parameter :: tol7= 0.0000001_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol9= 0.000000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol11=0.00000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol13=0.0000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp
 real(dp), parameter :: tol15=0.000000000000001_dp
 real(dp), parameter :: tol16=0.0000000000000001_dp
 real(dp), parameter :: tol20=0.00000000000000000001_dp

!real constants derived from sqrt(n.)
 real(dp), parameter :: sqrt2=1.4142135623730950488016887242096939_dp
 real(dp), parameter :: half_sqrt2=0.70710678118654752440084436210484697_dp
 real(dp), parameter :: sqrt3=1.7320508075688772935274463415058739_dp
 real(dp), parameter :: half_sqrt3=0.86602540378443864676372317075293693_dp
 real(dp), parameter :: sqrthalf=0.70710678118654752440084436210484697_dp

!Conversion factors of common use, not directly related to physical quantities.
 real(dp), parameter :: b2Mb=one/1024.0_dp**2  ! conversion factor bytes --> Mbytes
 real(dp), parameter :: b2Gb=b2Mb/1024.0_dp    ! conversion factor bytes --> Gbytes

!Real physical constants
!Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
!(from 2006 least squares adjustment)
 real(dp), parameter :: Bohr_Ang=0.52917720859_dp    ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_cmm1=219474.6313705_dp  ! 1 Hartree, in cm^-1
 real(dp), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
 real(dp), parameter :: Ha_meV=Ha_eV*1000_dp ! 1 Hartree, in meV
 real(dp), parameter :: Ha_K=315774.65_dp ! 1Hartree, in Kelvin
 real(dp), parameter :: Ha_THz=6579.683920722_dp ! 1 Hartree, in THz
 real(dp), parameter :: Ha_J=4.35974394d-18    !1 Hartree, in J
 real(dp), parameter :: e_Cb=1.602176487d-19 ! minus the electron charge, in Coulomb
 real(dp), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
 real(dp), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(dp), parameter :: HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0d+21 ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: Avogadro=6.02214179d23 ! per mole
 real(dp), parameter :: Ohmcm=two*pi*Ha_THz*ninth*ten ! 1 Ohm.cm in atomic units
 real(dp), parameter :: eps0=one/(four_pi*0.0000001_dp*299792458.0_dp**2)
 real(dp), parameter :: AmuBohr2_Cm2=e_Cb*1.0d20/(Bohr_Ang*Bohr_Ang)
 real(dp), parameter :: InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
 real(dp), parameter :: Sp_Lt_SI=2.99792458d8 ! speed of light in SI
 real(dp), parameter :: Sp_Lt=Sp_lt_SI/2.1876912633d6 ! speed of light in atomic units
 real(dp), parameter :: Time_Sec=2.418884326505D-17 !  Atomic unit of time, in seconds
 real(dp), parameter :: BField_Tesla=4.254383d-6 ! Tesla in a.u.

!Complex constants
 complex(dpc), parameter :: czero=(0._dp,0._dp)
 complex(dpc), parameter :: cone =(1._dp,0._dp)
 complex(dpc) ,parameter :: j_dpc=(0._dp,1.0_dp)

!Character constants
 character(len=1), parameter :: ch10 = char(10)
 character(len=fnlen),parameter :: ABI_NOFILE="__None__"

!Error codes used by the bindings.
 integer, parameter, public :: AB7_NO_ERROR                 =  0
 integer, parameter, public :: AB7_ERROR_OBJ                =  1
 integer, parameter, public :: AB7_ERROR_ARG                =  2
 integer, parameter, public :: AB7_ERROR_INVARS_ATT         =  3
 integer, parameter, public :: AB7_ERROR_INVARS_ID          =  4
 integer, parameter, public :: AB7_ERROR_INVARS_SIZE        =  5
 integer, parameter, public :: AB7_ERROR_SYM_NOT_PRIMITIVE  =  6
 integer, parameter, public :: AB7_ERROR_SYM_BRAVAIS_XRED   =  7
 integer, parameter, public :: AB7_ERROR_MIXING_ARG         =  8
 integer, parameter, public :: AB7_ERROR_MIXING_CONVERGENCE =  9
 integer, parameter, public :: AB7_ERROR_MIXING_INTERNAL    = 10
 integer, parameter, public :: AB7_ERROR_MIXING_INC_NNSLOOP = 11

!Parameters for LOG files treatment
!This variable tells the code if some lines have to be written in a LOG file
 logical, public, save :: do_write_log   =.true.

CONTAINS
!!***

!!****f* defs_basis/abi_io_redirect
!! NAME
!!  abi_io_redirect
!!
!! FUNCTION
!!  Redirect unit numbers (and|or) change the MPI communicator for the IO (output and log file).
!!  This routine can be used in client code (e.g. bigdft)
!!  that wants to call the abinit routines packed in an external library.
!!
!! INPUTS
!!  new_ab_out=new value for output file unit
!!  new_std_out=new value for standard output unit
!!  new_io_comm=new value for IO MPI communicator
!!
!! PARENTS
!!      abinit,aim,anaddb,band2eps,bsepostproc,conducti,cut3d,driver,fftprof
!!      initmpi_world,ioprof,kss2wfk,lapackprof,m_io_redirect,macroave
!!      memory_eval,mpi_setup,mrgddb,mrggkk,mrgscr,optic,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!
!! SOURCE

 subroutine abi_io_redirect(new_ab_out,new_std_out,new_io_comm,new_leave_comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_io_redirect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: new_std_out,new_ab_out,new_io_comm,new_leave_comm

!************************************************************************

 if (PRESENT(new_ab_out))  ab_out  = new_ab_out
 if (PRESENT(new_std_out)) std_out = new_std_out
 if (PRESENT(new_io_comm)) abinit_comm_output = new_io_comm
 if (.FALSE.) then
   ! Not used anymore
   write(std_out,*)new_leave_comm
 end if

 end subroutine abi_io_redirect
!!***

end module abi_defs_basis
!!***


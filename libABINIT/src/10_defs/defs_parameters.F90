!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_parameters
!! NAME
!! defs_parameters
!!
!! FUNCTION
!! Default parameters.
!! Typical use : parameters that are used in different routines, of which
!! the value is fixed here (still, they are not as fundamental as the
!! definitions in defs_basis).
!! Also, some input variables, of integer type, might be associated
!! with a name for selected values. One might associated here these
!! these particular values with a variable name.
!! Please, make sure that these parameters are easy to trace, in the
!! different routines they are used. Hence, give them a characteristic,
!! informative, name.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (PCasek,FF,XG,YMN)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

module defs_parameters

 use defs_basis

 implicit none

!- Set of parameters for the aim utility -----------------------------------
 real(dp), parameter :: aim_rhocormin=1.d-10  ! the minimal core density
 real(dp), parameter :: aim_epstep=0.5
 real(dp), parameter :: aim_rhomin=1.d-5,aim_dgmin=1.d-9,aim_dmaxcrit=5.d-2
 real(dp), parameter :: aim_dmin=1.d-3,aim_hmax=2.d7,aim_fac0=2.1_dp,aim_facmin=1.d-3
 real(dp), parameter :: aim_hmult=15._dp,aim_tiny=1.d-4,aim_snull=1.d-6
 real(dp), parameter :: aim_deltarmin=1.d-7
!the minimal length of one step following the gradient line
 real(dp), parameter :: aim_fac=1.2_dp,aim_drmin=1.d-5
 real(dp), parameter :: aim_dlimit=1.d-4,aim_dmaxcs=3.d-1
 real(dp), parameter :: aim_dpc0=1.d-2
 integer, parameter :: aim_maxstep=100
 real(dp), parameter :: aim_xymin=1.d-10
 integer, parameter :: aim_npmaxin=17
 real(dp), parameter :: aim_stmax=0.05
 real(dp), parameter :: aim_dmaxc1=1.d-1, aim_dmaxcl=5.d-2

!- Particular values of some input variables -----------------------------------

!- Possible values for ikhxc (choice for the TDDFT kernel) :

 integer, parameter :: ikhxc_NULL = 0 !No kernel.
 integer, parameter :: ikhxc_RPA  = 1 !RPA kernel.
 integer, parameter :: ikhxc_ALDA = 2 !ALDA kernel.
!For the PGG kernel, see M. Petersilka, U. J. Gossmann and E. K. U. Gross,
!Phys. Rev. Lett. 76, 1212 (1996).
 integer, parameter :: ikhxc_PGG  = 3 !PGG kernel.
!For the BPG kernel, see K. Burke, M. Petersilka and E. K. U. Gross,
!in "Recent Advances in Density Functional Methods", Vol. III,
!edited by P. Fantucci and A. Bencini (World Scientific, Singapore, 2002).
 integer, parameter :: ikhxc_BPG  = 4 !BPG kernel.
!For the following two kernels, see J. Dobson and J. Wang, Phys. Rev. B 62, 10038 (2000).
 integer, parameter :: ikhxc_EOK1 = 5 !Energy optimized kernel, scheme I (linear).
 integer, parameter :: ikhxc_EOK2 = 6 !Energy optimized kernel, scheme II (non-linear).

 character(len=4), parameter :: ikhxc_name(0:6) = (/"NULL","RPA ","ALDA","PGG ", &
&                                              "BPG ","EOK1","EOK2"/)

!- Possible values for idyson (choice for the solution of the Dyson equation) :

!Solve the Dyson equation as a linear system (direct inversion).
!See dyson_ls.F90.
 integer, parameter :: idyson_LS = 1
!Solve the Dyson equation as a differential equation with respect to the coupling
!constant. See dyson_de.F90.
 integer, parameter :: idyson_DE = 2
!Solve the Dyson equation as a self-consistent problem.
!See dyson_sc.F90.
 integer, parameter :: idyson_SC = 3

 character(len=25), parameter :: idyson_name(1:3) = (/"linear system          ", &
&                                                "differential equation  ", &
&                                                "self-consistent problem"/)

end module defs_parameters

!!***

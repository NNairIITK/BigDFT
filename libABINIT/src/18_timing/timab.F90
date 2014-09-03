!{\src2tex{textfont=tt}}
!!****f* ABINIT/timab
!! NAME
!!  timab
!!
!! FUNCTION
!!  Timing subroutine.  Calls machine-dependent "timein" which
!!  returns elapsed cpu and wall clock times in sec.
!!
!!  Depending on value of "option" routine will:
!!  (0) zero all accumulators
!!  (1) start with new incremental time slice for accumulator n
!!    also increase by one the counter for this accumulator
!!  (2) stop time slice; add time to accumulator n
!!  (3) not used (use now time_accu)
!!  (4) report time slice for accumlator n (not full time accumlated)
!!  (5) option to suppress timing (nn should be 0) or reenable it (nn /=0)
!!
!!  If, on first entry, subroutine is not being initialized, it
!!  will automatically initialize as well as rezero accumulator n.
!!  However, initialization SHOULD be done explicitly by the user
!!  so that it can be done near the top of his/her main routine.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed); NOT used if option=0
!!  option=see comment above
!!
!! OUTPUT
!!  on option=4 :
!!    tottim(2,nn)=accumulated time for accumulator nn; otherwise
!!     tottim is a dummy variable.
!!    option gives the number of times that the
!!     accumulator has been incremented
!!
!! PARENTS
!!      abinit,acfd_dyson,acfd_intexact,atm2fft,back_wf,bestwfs,bethe_salpeter
!!      calcdensph,cchi0,cgwf,cgwf3,corrmetalwf1,csigme,dielmt,dielmt2,dieltcel
!!      dotprod_g,dotprod_v,dotprod_vn,dotprodm_v,dotprodm_vn,driver,dyfnl3
!!      dyfro3,eltfrhar3,eltfrkin3,eltfrloc3,eltfrnl3,eltfrxc3,energy
!!      entropyrec,etotfor,fermisolverec,fftw3_fourdp,filnam_comm,filterpot
!!      first_rec,forces,forstrnps,forw_wf,fourdp,fourwf,fxphas,getgh1c,getghc
!!      getgsc,getngrec,gran_potrec,green_kernel,gstate,gstateimg,hartre
!!      hartre1,initylmg,inkpts,invars2,inwffil,inwffil3,kpgio,kpgsph,ladielmt
!!      lavnl,leave_test,lobpcgIIwf,lobpcgccIIwf,lobpcgccwf,lobpcgwf,loop3dte
!!      loper3,m_ab6_invars_f90,m_hidecudarec,m_screening,matrixelmt_g
!!      mean_fftr,meanvalue_g,mkcore,mkffnl,mklocl_realspace,mklocl_recipspace
!!      mkresi,mkrho,mkrho3,mkvxc3,mkvxcstr3,newkpt,newocc,newrho,newvtr
!!      newvtr3,nhatgrid,nlenergyrec,nonlinear,nonlop,nstdy3,nstwf3,odamix
!!      opernla_ylm,optics_paw,optics_paw_core,outkss,outscfcv,outwf,pareigocc
!!      partial_dos_fractions_paw,pawdenpot,pawdij,pawinit,pawmknhat,pawmkrhoij
!!      pawpolev,pawxc,pawxc3,pawxcm,pawxcm3,prctfvw1,prctfvw2,precon,precon2
!!      prep_fourwf,prep_getghc,prep_nonlop,projbd,pspheads_comm,pspini
!!      pw_orthon,recursion,recursion_nl,redgr,respfn,rhofermi3,rhohxc,rhotov
!!      rhotov3,rwwf,scfcv,scfcv3,screening,setsym,setvtr,sigma,sqnorm_g
!!      sqnorm_v,sqnormm_v,status,stress,strhar,subdiago,suscep,suscep_dyn
!!      suscep_kxc_dyn,suscep_stat,susk,susk_dyn,susk_dyn_pgg,susk_kxc_dyn
!!      suskmm,suskmm_dyn,suskmm_kxc_dyn,symrhg,symsgcube,tddft,timana
!!      vn_nl_rec,vtorho,vtorho3,vtorhorec,vtorhotf,vtowfk,vtowfk3,wfconv
!!      wfkfermi3,wfsinp,xcden,xcpot,zprecon3
!!
!! CHILDREN
!!      leave_new,papif_flops,papif_perror,timein,wrtout
!!
!! SOURCE
!!
!! TODO use m_iso_c_binding or m_abinit_c_iso_binding

#define C_INT INTEGER 
#define C_FLOAT REAL
#define C_LONG_LONG INTEGER*8

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine timab(nn,option,tottim)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing, except_this_one => timab
!End of the abilint section

 implicit none

#if defined HAVE_PAPI
#include "f90papi.h"
#endif

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nn,option
 !arrays
 real(dp),intent(out) :: tottim(2)

 tottim = (/ 0.0_dp, 0.0_dp /)
end subroutine timab
!!***

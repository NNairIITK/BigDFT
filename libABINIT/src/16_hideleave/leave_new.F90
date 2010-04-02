!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_new
!! NAME
!!  leave_new
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible
!!  parallelization.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, NCJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be
!!     written once only or
!!   'PERS' if the procs are calling the routine with different mesgs
!!     each to be written, or if one proc is calling the routine
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      abi_etsf_init,abinit,acfd_dyson,acfd_intexact,afterscfloop,anaddb
!!      appdig,append_cml,append_cml2,asria9,asrif9,asrprs,berryphase
!!      berryphase_new,besjm,bestwfs,bfactor,bigbx9,bldgrp,blok8
!!      bonds_lgth_angles,bound,bound_new,brdmin,calc_cs,calc_efg,calc_fc
!!      canat9,ccfft,cgwf,cgwf3,chkdilatmx,chkdpr,chkexi,chkgrp,chki8,chkilwf
!!      chkin9,chkinp,chkint_prt,chkneu,chknm8,chkorthsy,chkpawovlp,chkprimit
!!      chkr8,chkrp9,chkvars,clsopn,cmpar8,completeperts,constrf
!!      contract_dp_ge_val,contract_int_ge_val,contract_int_le_val
!!      contract_int_list,cprj_utils,cut3d,cvxclda,datafordmft,ddkten,delocint
!!      der_int,diel9,dielmt2,dieltcel,distrb2,dmft_solve,drivexc,dyson_sc
!!      echo_xc_name,elphon,elpolariz,eltxccore,energy,ep_setupqpt,etot3,ewald3
!!      ewald9,extraprho,extrapwf,fappnd,fermisolverec,fftpac,fftw,fillcell
!!      filterpot,find_getdtset,finddistrproc,findmin,fixsym,forstr,ftgam,ftgkk
!!      ftiaf9,fxphas,gath3,gensymshub,gensymshub4,gensymspgr,get_all_gkq
!!      get_full_kgrid,get_g_tiny,get_gkk_qpt_tr,get_tetra,getattribute,getcut
!!      getfreqsus,getgh1c,getghc,getkgrid,getlambda,getnel,getng,getshell
!!      gstate,gstateimg,gtblk9,handle_err_netcdf,hartre,hartre1,hartrestr
!!      hdr_check,hdr_io,hdr_io_etsf,hdr_io_netcdf,hermit,hessinit,importcml
!!      inarray,ingeo,ingeobld,ini_wf_etsf,init8,init_occ_ent,initang,initberry
!!      initberry3,initmpi_fft,initmpi_grid,initmpi_respfn,initmv,initrhoij
!!      initro,initylmg,inkpts,inpgkk,inprep8,inpspheads,inread,instrng,insy3
!!      int2char,int2char4,int_ang,intagm,integrate_gamma,integrate_gamma_tr
!!      inupper,invars0,invars1,invars1m,invars2,invars9,invcb,inwffil,inwffil3
!!      ioarr,ioddb8_in,iofn1,iofn2,irrzg,isfile,jellium,klocal,kpgsph,kpgstr
!!      kxc_alda,kxc_eok,ladielmt,lavnl,linemin,listkk,lobpcgIIwf,lobpcgccIIwf
!!      loper3,lwf,m_ab6_invars_f90,m_errors,m_green,m_libxc_functionals
!!      m_matlu,m_matrix,m_oper,m_paw_dmft,m_special_funcs,m_wffile
!!      mat_mlms2jmj,mat_slm2ylm,matcginv,matcginv_dpc,mati3inv,matrginv
!!      matrixelmt_g,mean_fftr,meanvalue_g,memana,metcon,metric,metstr,mka2f
!!      mka2fQgrid,mka2f_tr,mkcor3,mkcore,mkdenpos,mkeuler,mkffnl,mkfilename
!!      mkfskgrid,mkkpg,mklocl,mklocl_realspace,mklocl_recipspace
!!      mklocl_wavelets,mknesting,mknormpath,mkph_linwid,mkrho,mkvxc3,mkvxcstr3
!!      mlwfovlp_proj,mlwfovlp_projpaw,mlwfovlp_pw,mlwfovlp_radial
!!      mlwfovlp_setup,mlwfovlp_ylmfac,mlwfovlp_ylmfar,moddiel,moldyn,move
!!      mrgddb,mrggkk,nanal9,nderiv_gen,newkpt,newocc,newrho,newsp,newvtr
!!      newvtr3,nhatgrid,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk
!!      nmsq_pure_gkk_sumfs,nonlinear,nonlop,nonlop_pl,nonlop_ylm,normsq_gkq
!!      nres2vres,occeig,odamix,old_setmesh,opernl2,opernl3,opernl4b,optics_paw
!!      optics_paw_core,out1dm,outelph,outg2f,outphdos,outqmc,outscfcv,outwant
!!      outwf,outxfhist,overlap_g,parsefile,partial_dos_fractions_paw
!!      pawaccrhoij,pawdenpot,pawdensities,pawdij,pawdijhartree
!!      pawfrnhat_recipspace,pawinit,pawlsylm,pawpuxinit,pawxcmpositron
!!      pawxcpositron,pawxcsph,pawxcsphpositron,pawxcsum,ph1d3d,pl_deriv
!!      plm_coeff,plm_d2theta,plm_dphi,plm_dtheta,poslifetime,prcref,prcref_PMA
!!      prctfvw1,prctfvw2,prep_bandfft_tabs,prep_kpgio,print_ierr,print_psps
!!      prteigrs,prtfatbands,prtocc,prtph3,prtrhomxmn,prtspgroup,prttagm,psddb8
!!      psichi_renormalization,psolver_hartree,psolver_kernel,psolver_rhohxc
!!      psp10nl,psp11nl,psp1cc,psp1in,psp1nl,psp2in,psp3nl,psp4cc,psp5in,psp5nl
!!      psp6in,psp7in,psp7nl,psp8cc,psp8in,pspatm,pspini,psxml2ab,q0dy3
!!      ramansus,rchkgsheader,rdddb9,rdm,rdnpw,read_gkk,readeig,relaxpol,respfn
!!      rhofermi3,rhohxc,rhohxcpositron,rhotov3,rotmat,rrho,rsiaf9,rwwan
!!      scalapack,scalewf_nonlop,scfcge,scfcv,scfcv3,scfeig,scfopt,scphon
!!      scprqt,setnoccmmp,setrhoijpbe0,setup1,setup_G_rotation_old
!!      setup_positron,sg_ctrig,sg_fft,sg_fftpad,sg_fftpx,sg_fftrisc
!!      sg_fftrisc_2,sg_fftx,sg_ffty,sg_fftz,sg_fourwf,smallprim,smatrix
!!      smatrix_paw,smatrix_pawinit,smeared_delta,smpbz,sphere,sphereboundary
!!      sphericaldens,spin_current,status,subdiago,suscep_stat,suskmm
!!      suskmm_dyn,suskmm_kxc_dyn,sym_gkk,symatm,symbrav,symdet,symdij,symdm9
!!      symfind,symg,symkchk,symkpt,symlatt,symptgroup,symrelrot,symrhg
!!      symrhoij,symspgr,tddft,testkgrid,tetrahedron,thm9,timab,time_accu
!!      timein,transgrid,uderiv,ujdet,uniformrandom,vlocalstr
!!      vso_realspace_nonlop,vtorho,vtorho3,vtorhorec,vtowfk,vtowfk3,wfconv
!!      wffclose,wffile,wffopen,wffreadnpwrec,wfkfermi3,wfsinp,wght9
!!      wvl_free_type,wvl_init_type_proj,wvl_init_type_wfs,wvl_memory,wvl_mkrho
!!      wvl_newvtr,wvl_nl_gradient,wvl_rwwf,wvl_setboxgeometry,wvl_setngfft
!!      wvl_tail_corrections,wvl_utils,wvl_vtorho,wvl_wfsinp_disk
!!      wvl_wfsinp_reformat,wvl_wfsinp_scratch,xcacfd,xcden,xchcth,xchelu,xcpbe
!!      xcpot,xcpzca,xcspol,xctetr,xcwign,xcxalp,xfpack,xredxcart
!!
!! CHILDREN
!!      dump_config,leave_myproc,wrtout,xbarrier_mpi,xmpi_end,xmpi_me
!!      xmpi_nproc,xmpi_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine leave_new(mode_paral)

  implicit none

  !Arguments ------------------------------------
  character(len=4),intent(in) :: mode_paral

  print *,mode_paral, 'exiting...'
  stop
end subroutine leave_new
!!***

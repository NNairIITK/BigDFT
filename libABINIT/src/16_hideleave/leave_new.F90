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
!!  Copyright (C) 1998-2009 ABINIT group (DCA, XG, GMR, NCJ)
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
!!      abi_etsf_init,abinit,acfd_dyson,acfd_intexact,anaddb,anascr,appdig
!!      append_cml,append_cml2,asria9,asrif9,atm2fft,berryphase,berryphase_new,besjm
!!      bfactor,bigbx9,bldgrp,blok8,bonds_lgth_angles,bound,brdmin,canat9,ccfft
!!      cchi0,cchi0q0,cgwf,cgwf3,chkdilatmx,chkdpr,chkexi,chkgrp
!!      chki8,chkilwf,chkin9,chkinp,chkint,chkneu,chknm8,chkorthsy,chkpawovlp
!!      chkprimit,chkr8,chkrp9,clsopn,cmpar8,completeperts,constrf
!!      contract_dp_ge_val,contract_int_ge_val,contract_int_le_val
!!      contract_int_list,cppm2par,cppm3par,cppm4par,csigme,ctocprj,cvxclda
!!      ddkten,delocint,der_int,diel9,dielmt,dielmt2,dieltcel,distrb2
!!      dotprod_vn,dotprodm_vn,driver,drivergw,drivexc,dyson_sc,eig1fixed
!!      elphon,elpolariz,eltxccore,energy,ewald3,ewald9,fappnd,fermi,fftpac
!!      fftw,filterpot,findk,findmin,findshells,fixsym,forstr,fourdp,fourwf
!!      fresid,ftgam,ftgkk,ftiaf9,fxphas,gensymshub,gensymshub4,gensymspgr
!!      get_all_gkq,get_full_kgrid,get_g_tiny,get_gkk_qpt_tr,get_tetra
!!      getattribute,getcprj,getcut,getfreqsus,getghc,getkgrid,getlambda,getnel
!!      getng,getsc,getshell,gstate,gtblk9,handle_err_netcdf,hartre,hartre1
!!      hartrestr,hdr_check,hdr_init,hdr_io,hdr_io_etsf,hdr_io_netcdf,hermit
!!      hessinit,identk,importcml,inarray,ingeo,ingeobld,ini_wf_etsf,init8
!!      initang,initberry,initberry3,initmpi_fft,initmpi_gs,initmpi_respfn
!!      initmv,inkpts,inprep8,inread,instrng,insy3,int2char,int2char4,intagm
!!      integrate_gamma,integrate_gamma_tr,inupper,invars0,invars1,invars2
!!      invars9,invcb,inwffil,inwffil3,ioarr,ioddb8,iofn1,iofn2,irrzg,isfile
!!      klocal,kpgsph,kpgstr,kxc_alda,kxc_eok,ladielmt,lavnl,linemin,listkk
!!      lobpcgIIwf,lobpcgccIIwf,loper3,lwf,matcginv,mati3inv,matrginv
!!      matrixelmt_g,mean_fftr,meanvalue_g,memana,memerr,metcon,metric,metstr
!!      mka2f,mka2fQgrid,mka2f_tr,mkcor3,mkcore,mkeuler,mkffnl,mkfilename
!!      mkfskgrid,mkkpg,mklocl,mklocl_realspace,mklocl_recipspace
!!      mklocl_wavelets,mknesting,mkph_linwid,mkrho,mkvxc3,mkvxcstr3,moddiel
!!      moldyn,move,mrgddb,mrggkk,mrgscr,nanal9,nderiv_gen,netcdf_data_read
!!      netcdf_data_write,netcdf_def_acfd,netcdf_def_den,netcdf_def_dtset
!!      netcdf_def_wfs,netcdf_dims_read,netcdf_file_close,netcdf_file_create
!!      netcdf_get_acfd,netcdf_get_den,netcdf_get_dims_acfd,netcdf_get_dims_den
!!      netcdf_get_dims_dtset,netcdf_get_dims_wfs,netcdf_get_dtset
!!      netcdf_get_wfs,netcdf_handle_error,netcdf_put_acfd,netcdf_put_den
!!      netcdf_put_dtset,netcdf_put_wfs,newkpt,newocc,newrho,newsp,newvtr
!!      newvtr3,nhatgrid,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk,nonlinear,nonlop
!!      nonlop_pl,nonlop_ylm,normsq_gkq,nstdy3,nstwf3,occeig,operat
!!      opernl2,opernl3,opernl4b,optics_paw,out1dm,outelph,outkss,outqmc
!!      outwant,outwf,overlap_g,pawdenpot,pawdij,pawinit,pawmknhat
!!      pawmkrhoij,pawpuxinit,pawxc,pawxcdenm,pawxcm,ph1d3d,pl_deriv,plm_coeff
!!      plm_d2theta,plm_dphi,plm_dtheta,prcref,prctfvw1,prctfvw2,printbxsf
!!      printxsf,prteigrs,prtocc,prtph3,prtrhomxmn,prtspgroup,prttagm,psddb8
!!      psolver_hartree,psp1cc,psp1in,psp1nl,psp2in,psp3nl,psp4cc,psp5in,psp5nl
!!      psp7in,psp7nl,psp8cc,psp8in,pspatm,pspini,psxml2ab,q0dy3,ramansus
!!      rchkgsheader,rdddb9,rdkss,rdnpw,rdqps,rdscr,read_gkk,readeig,relaxpol
!!      respfn,rhofermi3,rhohxc,rhohxc_coll,rsiaf9,rwwan,rwwf,scalewf_nonlop
!!      scfcge,scfcv,scfcv3,scfeig,scfopt,scprqt,screening,setshells,setup1
!!      setup_hamilt,sg_ctrig,sg_fft,sg_fftpad,sg_fftpx,sg_fftrisc,sg_fftx
!!      sg_ffty,sg_fftz,sg_fourwf,sigma,smallprim,smatrix,smpbz,sphere
!!      sphereboundary,sphericaldens,status,stress,subdiago,surot,suskmm
!!      suskmm_dyn,suskmm_kxc_dyn,sym_gkk,symanal,symatm,symbrav,symdet,symdij
!!      symdm9,symfind,symg,symkchk,symkpt,symrelrot,symrhg,symrhoij,symspgr
!!      tddft,testkgrid,testlda,testscr,thm9,timab,time_accu,timein,transgrid
!!      uderiv,vlocalstr,vtorho,vtorho3,vtowfk,vtowfk3,wfconv
!!      wffclose,wffile,wffopen,wffreadnpwrec,wfkfermi3,wfsinp,wght9,wrscr
!!      wvl_init_type_proj,wvl_init_type_wfs,wvl_mkrho,wvl_rwwf
!!      wvl_setboxgeometry,wvl_vtorho,wvl_wfsinp,xcacfd,xcden,xchcth,xchelu
!!      xcpbe,xcpot,xcpzca,xcspol,xctetr,xcwign,xcxalp,xfpack,xredxcart,ylmc
!!      ylmcd
!!
!! CHILDREN
!!      leave_myproc,mpi_allreduce,mpi_barrier,mpi_comm_rank,mpi_finalize
!!      wrtout
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

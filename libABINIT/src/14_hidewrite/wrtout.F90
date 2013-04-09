!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout
!! NAME
!!  wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
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
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be written once only 
!!   'PERS' if the procs are calling the routine with different messages each to be written, 
!!          or if one proc is calling the routine
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abi_etsf_electrons_put,abi_etsf_geo_put,abi_etsf_init,abinit,acfd_dyson
!!      acfd_intexact,afterscfloop,anaddb,append_cml,append_cml2,asria9,asrif9
!!      asrprs,ass_leg_pol,atm2fft,atomden,berryphase,berryphase_new,besjm
!!      bestwfs,bethe_salpeter,bfactor,bigbx9,bldgrp,blok8,bonds_lgth_angles
!!      bound,bound_new,brdmin,calc_cs,calc_density,calc_efg,calc_exch,calc_fc
!!      calc_rpa_functional,calc_vHxc_braket,calcdensph,canat9,ccfft,cchi0
!!      cchi0q0,cexch_haydock,cgwf,cgwf3,chiscwrt,chkdilatmx,chkdpr,chkexi
!!      chkgrp,chki8,chkilwf,chkin9,chkinp,chkint_prt,chkneu,chknm8,chkorthsy
!!      chkpawovlp,chkph3,chkprimit,chkr8,chkrp9,chkvars,chneu9,clnup1,clnup2
!!      clsopn,cmpar8,completeperts,constrf,contract_dp_ge_val
!!      contract_int_ge_val,contract_int_le_val,contract_int_list,cprj_utils
!!      crho,csigme,cut3d,cvxclda,d3output,datafordmft,ddkten,debug_tools
!!      delocint,denfgr,der_int,diel9,dielmt,dielmt2,dieltcel,diisrelax,distrb2
!!      dmft_solve,dos_hdr_write,driver,drivexc,dsksta,dtsetcopy,dyson_sc
!!      echo_xc_name,eig1fixed,elast9,electrooptic,elphon,elpolariz,eltxccore
!!      energy,entropyrec,ep_setupqpt,etot3,ewald,ewald3,ewald4,ewald9
!!      exc_rmmdiis,exccoupl,excden,exceig,exch,extraprho,extrapwf,fconv,fermi
!!      fermisolverec,fftpac,fftw,fftwfn,fillcell,filterpot,find_getdtset
!!      findmin,first_rec,fixsym,forstr,forstrnps,fsumrule,ftgam,ftgkk,ftiaf9
!!      fxphas,gath3,gensymshub,gensymshub4,gensymspgr,get_all_gkq,get_fs_kpts
!!      get_full_gsphere,get_full_kgrid,get_g_tiny,get_gkk_qpt_tr,get_tetra
!!      getattribute,getcut,getdim_nloc,getfreqsus,getgh1c,getghc,getkgrid
!!      getlambda,getmpw,getnel,getng,getshell,getspinrot,gran_potrec
!!      green_kernel,gstate,gstateimg,gtblk9,gtdyn9,gw2wfk,gw_driver,gw_tools
!!      gwcompleteness,handle_err_netcdf,hartre,hartre1,hartrestr,haydock
!!      hdr_check,hdr_io,hdr_io_etsf,hdr_io_netcdf,hdr_vs_dtset,herald,hermit
!!      hessinit,importcml,inarray,incomprs,ingeo,ingeobld,ini_wf_etsf,init8
!!      init_occ_ent,initang,initberry,initberry3,initmpi_respfn,initmv
!!      initrhoij,initro,initwf,initylmg,inkpts,inpgkk,inprep8,inpspheads
!!      instr9,instrng,insy3,int2char,int2char4,int_ang,intagm,integrate_gamma
!!      integrate_gamma_tr,inupper,invars0,invars1,invars1m,invars2,invars2m
!!      invars9,invcb,inwffil,inwffil3,ioarr,ioddb8_in,ioddb8_out,iofn1,iofn2
!!      ioniondist,irrzg,isfile,jellium,klocal,kpgio,kpgsph,kpgstr
!!      kramerskronig,ks_ddiago,kxc_alda,kxc_eok,ladielmt,lattice,lavnl
!!      leave_new,leave_test,linemin,listkk,lobpcgIIwf,lobpcgccIIwf,lobpcgccwf
!!      lobpcgwf,loop3dte,loper3,lwf,m_ab6_invars_f90,m_abilasi,m_atom
!!      m_bands_sym,m_bs_defs,m_bz_mesh,m_coulombian,m_crystal,m_dyson_solver
!!      m_ebands,m_errors,m_fft_mesh,m_fftw3,m_geometry,m_green,m_gsphere
!!      m_gwdefs,m_hamiltonian,m_hidecudarec,m_initcuda,m_io_kss,m_io_screening
!!      m_libxc_functionals,m_matlu,m_matrix,m_melemts,m_numeric_tools,m_oper
!!      m_paw_dmft,m_paw_pwij,m_paw_slater,m_paw_toolbox,m_phdos,m_ppmodel
!!      m_pretty_rec,m_qparticles,m_radmesh,m_rec,m_screening,m_self
!!      m_sigma_results,m_special_funcs,m_wffile,m_wfs,mat_mlms2jmj,mat_slm2ylm
!!      matcginv,matcginv_dpc,mati3inv,matrginv,matrixelmt_g,mean_fftr
!!      meanvalue_g,memana,memkss,memorf,memory,metcon,metric,metstr,mka2f
!!      mka2fQgrid,mka2f_tr,mkcor3,mkcore,mkdenpos,mkeuler,mkffnl,mkfilename
!!      mkfskgrid,mkifc9,mkkpg,mklocl,mklocl_realspace,mklocl_recipspace
!!      mklocl_wavelets,mknesting,mknormpath,mkph_linwid,mkphbs,mkphdos
!!      mkqptequiv,mkrho,mkrho3,mkvxc3,mkvxcstr3,mlwfovlp,mlwfovlp_proj
!!      mlwfovlp_projpaw,mlwfovlp_pw,mlwfovlp_qp,mlwfovlp_radial
!!      mlwfovlp_seedname,mlwfovlp_setup,mlwfovlp_ylmfac,mlwfovlp_ylmfar
!!      moddiel,moldyn,move,mrgddb,mrggkk,mrgscr,mv_3dte,my_calc_wfwfg,nanal9
!!      nderiv_gen,newfermie1,newkpt,newocc,newrho,newsp,newvtr,newvtr3
!!      nhatgrid,nlenergyrec,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk
!!      nmsq_pure_gkk_sumfs,nonlinear,nonlop,nonlop_pl,nonlop_ylm,normsq_gkq
!!      nres2vres,nselt3,nstdy3,occeig,occred,odamix,old_setmesh,opernl2
!!      opernl3,opernl4b,optics_paw,optics_paw_core,orthonormalize,out1dm
!!      outelph,outg2f,outkss,outphdos,outqmc,outscfcv,outvars,outwant,outwf
!!      outxfhist,overlap_g,pareigocc,parsefile,partial_dos_fractions_paw
!!      paw_mknewh0,paw_qpscgw,pawaccrhoij,pawdenpot,pawdensities,pawdij
!!      pawdijhartree,pawfrnhat_recipspace,pawinit,pawlsylm,pawmkaewf
!!      pawmkrhoij,pawprt,pawpupot,pawpuxinit,pawuenergy,pawuj_det,pawuj_red
!!      pawuj_utils,pawxcmpositron,pawxcpositron,pawxcsph,pawxcsphpositron
!!      pawxcsum,pawxenergy,pawxpot,pclock,ph1d3d,phfrq3,piezo9,pl_deriv
!!      plm_coeff,plm_d2theta,plm_dphi,plm_dtheta,polcart,poslifetime,prcref
!!      prcref_PMA,prctfvw1,prctfvw2,precon,precon2,prep_bandfft_tabs
!!      prep_kpgio,print_ierr,print_ij,print_psps,prmat,projbd,prt_cml2
!!      prteigrs,prtene,prtene3,prtfatbands,prtocc,prtph3,prtrhomxmn,prtspgroup
!!      prttagm,prtxf,prtxvf,psddb8,psichi_renormalization,psolver_hartree
!!      psolver_kernel,psolver_rhohxc,psp10in,psp10nl,psp1cc,psp1in,psp1nl
!!      psp2in,psp2lo,psp3in,psp3nl,psp4cc,psp5in,psp5nl,psp6in,psp7in,psp7nl
!!      psp8cc,psp8in,psp9in,pspatm,pspini,pspnl_hgh_rec,pspnl_operat_rec
!!      psxml2ab,q0dy3,ramansus,randac,rchkgsheader,rdddb9,rdm,rdnpw,read_gkk
!!      readeig,recursion,relaxpol,remove_inversion,respfn,rhofermi3,rhohxc
!!      rhohxcpositron,rhotov3,rotmat,rrho,rsiaf9,rwwan,scalapack
!!      scalewf_nonlop,scfcge,scfcv,scfcv3,scfeig,scfopt,scphon,scprqt
!!      screening,set_gwdistrb,setnoccmmp,setrhoijpbe0,setshells,setsymrhoij
!!      setup1,setup2,setup_G_rotation_old,setup_bethe_salpeter,setup_positron
!!      setup_qmesh,setup_screening,setup_sigma,setvtr,sg_ctrig,sg_fft
!!      sg_fftpad,sg_fftpx,sg_fftrisc,sg_fftrisc_2,sg_fftx,sg_ffty,sg_fftz
!!      sg_fourwf,shellstruct,sigma,smallprim,smatrix,smatrix_paw
!!      smatrix_pawinit,smeared_delta,smpbz,spectra,spectral,sphere
!!      sphereboundary,sphericaldens,spin_current,split_work,split_work2,status
!!      stress,subdiago,sumrule,suscep,suscep_dyn,suscep_kxc_dyn,suscep_stat
!!      suskmm,suskmm_dyn,suskmm_kxc_dyn,sym_gkk,symanal,symatm,symaxes,symbrav
!!      symdet,symdij,symdm9,symfind,symg,symkchk,symkpt,symlatt
!!      symmetrize_afm_chi0,symmultsg,symph3,symplanes,symptgroup,symq3
!!      symrelrot,symrhg,symrhoij,symspgr,tddft,testkgrid,tetrahedron,thm9
!!      thmeig,timab,timana,time_accu,timein,transgrid,uderiv,ujdet,vlocalstr
!!      vso_realspace_nonlop,vtorho,vtorho3,vtorhorec,vtorhotf,vtowfk,vtowfk3
!!      wfconv,wffclose,wffile,wffopen,wffreadnpwrec,wfkfermi3,wfsinp,wght9
!!      write_header_moldynnetcdf,write_moldynvaluenetcdf,wrtloctens
!!      wvl_free_type,wvl_init_type_proj,wvl_init_type_wfs,wvl_memory,wvl_mkrho
!!      wvl_newvtr,wvl_nl_gradient,wvl_rwwf,wvl_setboxgeometry,wvl_setngfft
!!      wvl_tail_corrections,wvl_utils,wvl_vtorho,wvl_wfsinp_disk
!!      wvl_wfsinp_reformat,wvl_wfsinp_scratch,xc_kernel,xcacfd,xcden,xchcth
!!      xchelu,xcpbe,xcpot,xcpzca,xcspol,xctetr,xcwign,xcxalp,xfpack,xredxcart
!!      zorthonormalize,zprecon3
!!
!! CHILDREN
!!      wrtout_myproc,xmpi_me,xmpi_nproc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine wrtout(unit,message,mode_paral)
  use defs_basis
  implicit none
  include 'mpif.h'

  !Arguments ------------------------------------
  integer,intent(in) :: unit
  character(len=4),intent(in) :: mode_paral
  character(len=500),intent(inout) :: message

  integer :: ierr, iproc

  if (unit == ab_out .or. index(message, "ERROR") > 0) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (trim(mode_paral) == "COLL") then
        if (iproc == 0) write(*, "(1x,A)") trim(message)
     else
        write(*, "(I03,2x,A)") iproc, trim(message)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end if
end subroutine wrtout
!!***

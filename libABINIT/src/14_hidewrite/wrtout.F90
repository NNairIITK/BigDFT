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
!!  Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!   The routine uses optional arguments, therefore the interface must be explicit.
!!   Be careful when writing CPP macros that use wrtout since abilint won't see the call
!!   and no interface will be added to the source file.
!!
!! PARENTS
!!      abinit,acfd_dyson,acfd_intexact,afterscfloop,anaddb,append_cml2
!!      append_xyz,atm2fft,atomden,berryphase,berryphase_new,bethe_salpeter
!!      bonds_lgth_angles,bsepostproc,calc_cs,calc_efg,calc_fc
!!      calc_optical_mels,calc_rpa_functional,calc_sigc_me,calc_sigx_me
!!      calc_ucrpa,calc_vhxc_me,calcdensph,cchi0,cchi0q0,cchi0q0_intraband,cgwf
!!      cgwf3,check_completeness,chiscwrt,chkdpr,chkinp,chkint_prt,chkneu
!!      chkpawovlp,clnup1,clnup2,cohsex_me,compute_kgb_indicator,compute_levels
!!      constrf,d3output,datafordmft,debug_tools,defs_scalapack,deloc2xcart
!!      denfgr,dfpt_write_cg,diel9,dielmt,dieltcel,dmft_solve,dos_hdr_write
!!      driver,dsksta,dyson,echo_xc_name,elast9,eliashberg_1d,elphon,elpolariz
!!      entropyrec,ep_fs_weights,ep_setupqpt,evdw_wannier,ewald4
!!      exc_build_block,exc_build_ham,exc_den,exc_diago,exc_iterative_diago
!!      exc_plot,exc_spectra,fconv,fermi_green,fermisolverec,fftprof
!!      find_getdtset,finddistrproc,findmin,findminscf,first_rec,fred2fdeloc
!!      fsumrule,gaus_dos,get_all_gkq,get_fs_bands,get_npert_rbz,get_nv_fs_en
!!      get_nv_fs_temp,get_tau_k,getcgqphase,getcut,getdim_nloc,getfreqsus
!!      getghc,getmpw,getnel,getng,getshell,getspinrot,gran_potrec,green_kernel
!!      gstate,gstateimg,gw_tools,gwcompleteness,haydock,haydock_psherm
!!      hdr_vs_dtset,hermit,hubbard_one,importcml,importxyz,impurity_solve
!!      inarray,ingeo,ingeobld,initberry,initberry3,initmpi_grid,initorbmag
!!      initro,initwf,inkpts,inpgkk,inpspheads,instr9,instrng,insy3,intagm
!!      integrate_gamma,integrate_gamma_alt,integrate_gamma_tr
!!      integrate_gamma_tr_lova,invars1,invars1m,invars2,inwffil,inwffil3,ioarr
!!      ioddb8_out,ioniondist,ioprof,irrzg,isfile,kpgio,kramerskronig,ks_ddiago
!!      kss2wfk,ladielmt,lapackprof,lavnl,ldau_self,leave_new,linemin,lobpcgwf
!!      local_ks_green,loop3dte,loper3,m_abi_etsf,m_abilasi,m_anaddb_dataset
!!      m_atom,m_bands_sym,m_bfgs,m_bs_defs,m_bse_io,m_bz_mesh,m_cgtools,m_chi0
!!      m_commutator_vkbr,m_crystal,m_ddb_blk,m_dynmat,m_dyson_solver,m_ebands
!!      m_eet,m_energy,m_errors,m_eval_lotf,m_exit,m_fft,m_fft_mesh,m_fft_prof
!!      m_fftcore,m_fftw3,m_gamma,m_gaussfit,m_geometry,m_gpu_detect,m_green
!!      m_gsphere,m_hamiltonian,m_header,m_hidecudarec,m_hu,m_ifc,m_initcuda
!!      m_io_gkk,m_io_kss,m_io_screening,m_iterators,m_kxc,m_libxc_functionals
!!      m_lotf,m_matlu,m_matrix,m_melemts,m_mep,m_numeric_tools,m_oper,m_paw_an
!!      m_paw_dmft,m_paw_ij,m_paw_pwij,m_paw_slater,m_pawang,m_pawdij,m_pawfgr
!!      m_pawfgrtab,m_pawio,m_pawpsp,m_pawrad,m_pawrhoij,m_pawtab,m_phdos
!!      m_pimd,m_ppmodel,m_pptools,m_pred_lotf,m_pretty_rec,m_psps,m_ptgroups
!!      m_qparticles,m_rec,m_screen,m_screening,m_self,m_shirley
!!      m_sigma_results,m_sphharm,m_vcoul,m_wffile,m_wfk,m_wfs,m_work_var_lotf
!!      m_xc_vdw,m_xpapi,mag_constr_e,mag_out,mblktyp1,mblktyp5,memana,memorf
!!      memory,metric,mka2f,mka2fQgrid,mka2f_tr,mka2f_tr_lova,mkcore_paw
!!      mkcore_wvl,mkfilename,mkfskgrid,mklocl_recipspace,mklocl_wavelets
!!      mknormpath,mkph_linwid,mkphbs,mkqptequiv,mkrho,mkrho3,mlwfovlp
!!      mlwfovlp_proj,mlwfovlp_projpaw,mlwfovlp_pw,mlwfovlp_qp
!!      mlwfovlp_seedname,mlwfovlp_setup,mover,mpi_setup,mrgddb,mrggkk,mrgscr
!!      multipoles_fftr,mv_3dte,my_calc_wfwfg,new_integrate_gamma
!!      new_integrate_gamma_tr,new_integrate_gamma_tr_lova,newfermie1,newkpt
!!      newocc,newton,nlenergyrec,nonlinear,normsq_gkq,nselt3,nstdy3,nstpaw3
!!      out1dm,outelph,outgkk,outkss,outphbtrap,outphdos,outqmc,outscfcv
!!      outvars,outwant,outwf,paw_mknewh0,paw_qpscgw,pawdenpot,pawdensities
!!      pawmkaewf,pawmkrhoij,pawprt,pawpuxinit,pawuenergy,pawuj_det,pawuj_red
!!      pawuj_utils,pawxenergy,piezo9,pimd_nosehoover_nvt,polcart,posdoppler
!!      poslifetime,precon2,pred_delocint,pred_isokinetic,pred_isothermal
!!      pred_langevin,pred_nose,pred_verlet,predictimg,prep_calc_ucrpa,prt_cml2
!!      prtefield,prteigrs,prtene,prtene3,prtfatbands,prtimg,prtph3,prtrhomxmn
!!      prtspgroup,prtvsound,prtxf,prtxfase,prtxvf,psichi_renormalization
!!      psolver_hartree,psolver_kernel,psolver_rhohxc,psp10in,psp1in,psp2in
!!      psp2lo,psp3in,psp5in,psp6in,psp7wvl2,psp8in,psp9in,pspatm_abinit
!!      pspatm_pspio,pspini,pspnl_hgh_rec,pspnl_operat_rec,psxml2ab
!!      qmc_prep_ctqmc,randac,random_stopping_power,read_gkk,recursion_nl
!!      remove_inversion,respfn,rotate_rho,rotmat,scfcge,scfcv,scfcv3,scfeig
!!      scfopt,scphon,scphon_build_qsym_map,scphon_dynmat_to_freq2
!!      scphon_free_energy,scphon_supercell_vectors_init,scprqt,screening
!!      setnoccmmp,setrhoijpbe0,setsymrhoij,setup1,setup2,setup_bse
!!      setup_bse_interp,setup_positron,setup_screening,setup_sigma,shellstruct
!!      sigma,smpbz,spectral_function,stress,sumrule,suscep,sym_gkk,symanal
!!      symatm,symaxes,symcharac,symkchk,symkpt,symlatt,symmultsg,symph3
!!      symplanes,symq3,symspgr,tddft,testkgrid,tetrahedron,thm9,thmeig,timana
!!      uderiv,ujdet,update_eb_field_vars,vdw_dftd2,vdw_kernelgen,vtorho
!!      vtorhorec,vtorhotf,vtowfk,vtowfk3,wfconv,wfd_mkrho,wfd_pawrhoij
!!      wfkfermi3,wfsinp,wrt_moldyn_netcdf,wrtloctens,wvl_denspot_set
!!      wvl_descr_atoms_set_sym,wvl_hpsitopsi,wvl_initro,wvl_memory,wvl_mkrho
!!      wvl_nl_gradient,wvl_projectors_set,wvl_psitohpsi,wvl_rwwf
!!      wvl_setboxgeometry,wvl_setngfft,wvl_tail_corrections,wvl_wfs_set
!!      wvl_wfsinp_disk,wvl_wfsinp_reformat,wvl_wfsinp_scratch,xcacfd,zprecon3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wrtout(unit,msg,mode_paral,do_flush)

 use defs_basis

 use m_xmpi,      only : xmpi_world, xcomm_rank, xcomm_size
 use m_io_tools,  only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout'
 use interfaces_14_hidewrite, except_this_one => wrtout
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: do_flush

!Local variables-------------------------------
 integer :: comm,me,nproc
 integer,save :: master=0
 logical :: my_flush
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if ((unit == std_out).and.(.not.do_write_log)) RETURN
 if (unit == dev_null) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral
 my_flush = .false.; if (PRESENT(do_flush)) my_flush = do_flush

!Communicator is xmpi_world by default, except for the parallelization over images
 if (abinit_comm_output/=-1) then
   comm=abinit_comm_output
 else
   comm=xmpi_world
 end if

!Determine who I am in COMM_WORLD
 nproc = xcomm_size(comm)
 me    = xcomm_rank(comm)

 if( (my_mode_paral=='COLL') .or. (nproc==1) ) then
   if (me==master) then
     call wrtout_myproc(unit, msg, do_flush=my_flush)
   end if

 else if (my_mode_paral=='PERS') then
   call write_lines(unit,msg)

   ! Flush unit
   if (my_flush) then
     call flush_unit(unit)
   end if

 else if (my_mode_paral=='INIT') then
   master=unit

 else
   write(string,'(7a)')ch10,&
&   'wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)
 end if

end subroutine wrtout
!!***

!!****f* ABINIT/wrtout_myproc
!! NAME
!!  wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc. For parallel or sequential output use wrtout()
!!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!!  Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!  [mpicomm]= Optional argument. If present, no printing is done
!!             Variables iexit, nwarning and ncomment are
!!             summed over the mpicomm communicator
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstateimg,wrtout
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wrtout_myproc(unit,message,mpicomm,do_flush) ! optional argument

 use defs_basis
 use m_profiling

 use m_xmpi,      only : xmpi_sum
 use m_io_tools,  only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout_myproc'
 use interfaces_14_hidewrite, except_this_one => wrtout_myproc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message
 integer,intent(in),optional :: mpicomm
 logical,optional,intent(in) :: do_flush

!Local variables-------------------------------
!scalars
 integer,save :: iexit=0,ncomment=0,nwarning=0
 integer :: ierr
 logical :: print_std_err
!arrays
 integer :: buf(3)

!******************************************************************

!MG: TODO:
! move iexit, ncomment and nwanings to m_errors.
! provide an API to retrieve these values and to MPI sum them.

!When I/O are redirected, it is sometimes necessary to reduce counters (saved) values;
!this can be done by passing mpicomm optional argument to the routine In that case, no printing is done.
 if (present(mpicomm)) then
   buf(1)=iexit;buf(2)=ncomment;buf(3)=nwarning
   call xmpi_sum(buf,mpicomm,ierr)
   iexit=buf(1);ncomment=buf(2);nwarning=buf(3)
   if (iexit/=0) iexit=1
   return
 end if

 print_std_err=(unit==std_out.and.(index(trim(message),'BUG')/=0.or.index(trim(message),'ERROR')/=0))

 call write_lines(unit,message)
 if (print_std_err) then
   call write_lines(std_err,message)
 end if

 if( index(trim(message),'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action : contact ABINIT group.'
   if (print_std_err) write(std_err, '(a)' ) '  Action : contact ABINIT group.'
   write(unit,*)
   if (print_std_err) write(std_err,*)
 end if

 if( index(trim(message),'BUG') /= 0   .or. index(trim(message),'Calculation completed') /= 0 )then
   if(nwarning<10000 .and. ncomment<1000)then
     write(unit, '(a,i5,a,i4,a)' ) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   else
     write(unit, '(a,i6,a,i6,a)' ) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   end if
   if(iexit/=0)then
     write(unit, '(a)' ) ' Note : exit requested by the user.'
   end if
 end if

 if( index(trim(message),'Exit') /= 0 )then
   iexit=1
 end if

!Count the number of warnings and comments. Only take into
!account unit std_out, in order not to duplicate these numbers.
 if( index(trim(message),'WARNING') /= 0 .and. unit==std_out )then
   nwarning=nwarning+1
 end if
 if( index(trim(message),'COMMENT') /= 0 .and. unit==std_out )then
   ncomment=ncomment+1
 end if

 ! Flush unit
 if (present(do_flush)) then
   if (do_flush) then
     call flush_unit(unit)
   end if
 end if

#ifdef DEBUG_MODE
 call flush_unit(unit)
 if (print_std_err) then
   call flush_unit(std_err)
 end if
#endif

end subroutine wrtout_myproc
!!***

!!****f* ABINIT/write_lines
!! NAME
!!  write_lines
!!
!! FUNCTION
!!  This routine receives a string, split the message in lines according to the 
!!  ch10 character and output the text to the specified unit 
!!
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      wrtout
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_lines(unit,message)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_lines'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message

!Local variables-------------------------------
!scalars
 integer :: msg_size,ii,jj,rtnpos

!******************************************************************

 msg_size = len_trim(message)

 if (msg_size == 0) then
   write(unit,*)
   return 
 end if

 ! Here, split the message, according to the char(10) characters (carriage return). 
 ! This technique is portable accross different OS.
 rtnpos = index(message,ch10)

 if (rtnpos == 0) then
   write(unit,"(a)")message(1:msg_size)
   return
 end if 

 ii = 1; jj = rtnpos
 do 
   if (ii == jj) then
     write(unit,*)
   else
     write(unit, '(a)' ) message(ii:jj-1)
   end if
   ii = jj + 1
   if (ii > msg_size) exit
   jj = index(message(ii:msg_size),ch10) 
   if (jj == 0) then 
     ! Will write the last line at the next iteration and exit .
     jj = msg_size + 1
   else
     jj = jj + ii - 1
   end if
   !write(*,*)"ii, jj, msg_size",ii, jj, msg_size
 end do

 ! This is needed to preserve the od behaviour: a ch10 at the 
 ! end of the string was causing an extra newline!
 if (message(msg_size:msg_size) == ch10) write(unit,*)

end subroutine write_lines
!!***

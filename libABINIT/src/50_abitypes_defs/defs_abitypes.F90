!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_abitypes
!! NAME
!! defs_abitypes
!!
!! FUNCTION
!! This module contains definitions of hig-level structured datatypes for the
!! ABINIT package.
!!
!! If you are sure a new high-level structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are
!! well documented, it is a shame ...).
!! Do not forget : you will likely be the major winner if you document
!! properly.
!! Proper documentation of a structured datatype means :
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!
!! List of datatypes :
!! * aim_dataset_type : the "dataset" for aim
!! * anaddb_dataset_type : the "dataset" for anaddb
!! * dataset_type : the "dataset" for the main abinit code
!! * bandfft_kpt_type : the "dataset" for triple band-fft-kpt parallelization
!! * MPI_type : the data related to MPI parallelization
!!
!! COPYRIGHT
!! Copyright (C) 2001-2009 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

module defs_abitypes

 use defs_basis
 use defs_datatypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : atoms_data
#endif

 implicit none

!Structures
!!***

!!****t* defs_abitypes/aim_dataset_type
!! NAME
!! aim_dataset_type
!!
!! FUNCTION
!! The aim_dataset_type structured datatype
!! gathers all the input variables for the aim code
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

 type aim_dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Since all these input variables are described in the aim_help.html
! file, they are not described in length here ...

! Integer
  integer :: crit,denout,dltyp,gpsurf,irho,ivol,lapout,nsa,nsb,nsc
  integer :: ngrid(3)
  integer :: batom  !! Warning : corresponds to the input variable atom
  integer :: foll   !! Warning : corresponds to the input variable follow
  integer :: isurf  !! Warning : corresponds to the input variable surf
  integer :: irsur  !! Warning : corresponds to the input variable rsurf
  integer :: nph    !! Warning : corresponds to the input variable nphi
  integer :: npt    !! Warning : corresponds to the input variable inpt
  integer :: nth    !! Warning : corresponds to the input variable ntheta
  integer :: plden  !! Warning : not documented in help file ?!

! Real
  real(dp) :: atrad,coff1,coff2,dpclim,folstp,lgrad,lgrad2,lstep,lstep2,&
&  maxatd,maxcpd,phimax,phimin
  real(dp) :: foldep(3),scal(3),vpts(3,4)
  real(dp) :: dr0    !! Warning : correspond to the input variable radstp
  real(dp) :: phi0   !! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   !! Warning : correspond to the input variable ratmin
  real(dp) :: th0    !! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax !! Warning : correspond to the input variable thetamax
  real(dp) :: themin !! Warning : correspond to the input variable thetamin

 end type aim_dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/anaddb_dataset_type
!! NAME
!! anaddb_dataset_type
!!
!! FUNCTION
!! The anaddb_dataset_type structured datatype
!! gather all the input variables for the anaddb code.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

 type anaddb_dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Since all these input variables are described in the anaddb_help.html
! file, they are not described in length here ...
! Integer
  integer :: alphon,asr,brav,chneut,dieflag,dipdip,ep_scalprod,eivec,elaflag,elphflag,enunit
  integer :: gkk2write,gkk_rptwrite,gkqwrite
  integer :: iavfrq,ifcana,ifcflag,ifcout,ifltransport,instrflag,natfix,natifc,natom
  integer :: nchan,nfreq,ngrids,nlflag,nph1l,nph2l,nqpath
  integer :: nqshft,nsphere,nstrfix,ntemper,nwchan
  integer :: piezoflag,polflag,prtdos,prtmbm,prtfsurf,prtnest,ramansr
  integer :: relaxat,relaxstr,rfmeth,selectz,symdynmat,telphint,ep_keepbands,thmflag
  integer :: ep_prt_yambo
  integer :: qgrid_type, ep_nqpt
  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: istrfix(6),ng2qpt(3),kptrlatt(3,3)
  integer :: ep_b_min, ep_b_max
  integer :: symgkq

! Real(dp)
  real(dp) :: a2fsmear,dosdeltae,dossmear,dostol,elphsmear,elph_fermie,frmax,frmin
  real(dp) :: temperinc,tempermin,thmtol,mustar,rifcsph
  real(dp) :: q1shft(3,4),q2shft(3),targetpol(3)

! Integer pointers
  integer, pointer :: atifc(:)    ! atifc(natom) WARNING : there is a transformation
                                  ! of this input variable, in chkin9
                                  ! This should be changed ...
  integer, pointer :: iatfix(:)   ! iatfix(natom)

! Real pointers
  real(dp), pointer :: qnrml1(:)  ! qnrml1(nph1l)
  real(dp), pointer :: qnrml2(:)  ! qnrml2(nph2l)
  real(dp), pointer :: qpath(:,:) ! qpath(3,nqpath)
  real(dp), pointer :: qph1l(:,:) ! qph1l(3,nph1l)
  real(dp), pointer :: qph2l(:,:) ! qph2l(3,nph2l)

  real(dp), pointer :: ep_qptlist(:,:) ! qph2l(3,ep_nqpt)

 end type anaddb_dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/dataset_type
!! NAME
!! dataset_type
!!
!! FUNCTION
!! The dataset_type structured datatype gather all the input variables,
!! except those that are labelled NOT INTERNAL.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset.
!! The "evolving" input variables are also stored, with their
!! name appended with _orig, to make clear that this is the original
!! value, decided by the user, and not a possibly modified, intermediate value.
!! The following input variables are NOT INTERNAL, that is, they
!! are input variables used to determine other input variables,
!! after suitable processing, and do not appear anymore afterwards
!! (so, they do not appear as components of a dataset_type variable) :
!! cpuh,cpum(but cpus is present),fband,kptbounds,ndivk,ndism,nobj,
!! objaat,objbat,objaax,objbax,objan,objbn,objarf,objbrf,objaro,objbro
!! objatr,objbtr,vaclst,vacuum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

 type dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Since all these input variables are described in the abinis_help.html
! file, they are not described in length here ...
! Integer
  integer :: accesswff,awtr,bandpp,berryopt,brvltt,ceksph,chkexit,chkprim,&
&  delayperm,dmatpuopt,dmatudiag,enunit,exchn2n3d,fftgw,fft_opt_lob,&
&  frzfermi,getacfd,&
&  getcell,getddk,getden,getkss,getocc,getqps,getscr,getsuscep,getvel,getwfk,&
&  getwfq,getxcart,getxred,get1den,get1wf,gwcalctyp,gwcomp,gwgamma,gw_nqlwl,gw_sigxcore,gwmem,gwpara,iboxcut,&
&  icoulomb,icutcoul,idyson,ieig2rf,iextrapwf,ikhxc,inclvkb,intexact,intxc,ionmov,&
&  iprcch,iprcel,iprctfvw,iprcfc,irdddk,irdden,irdkss,irdqps,irdscr,irdsuscep,irdwfk,irdwfq,ird1wf,&
&  iscf,isecur,istatr,istatshft,ixc,ixcpositron,&
!  jdtset contains the actual number of the dataset
&  jdtset,jellslab,kpara,kptopt,kssform,ldgapp,localrdwf,mband,mffmem,mgfft,mgfftdg,&
&  mkmem,mkqmem,mk1mem,nnos,&
&  mpw,mqgrid,mqgriddg,natom,natpawu,natrd,natsph,natvshift,nbandkss,nbandsus,nbdblock,nbdbuf,&
&  nberry,nconeq,nctime,ndtset,ndyson,&
&  nfft,nfftdg,nfreqim,nfreqre,nfreqsp,nfreqsus,ngeohist,ngroup_rf,nkptgw,nkpt,nline,&
&  nnsclo,nomegasf,nomegasi,nomegasrd,npack,npara,npband,npfft,npkpt,npsp,npspalch,npulayit,&
&  npweps,npwkss,npwsigx,npwwfn,nqpt,nqptdm,nscforder,&
&  nsheps,nshiftk,nshsigx,nshwfn,nspden,nspinor,nsppol,nstep,nsym,ntime,&
&  ntypalch,ntypat,ntyppure,nwfshist,occopt,optcell,optdriver,&
&  optforces,optfreqsus,optnlxccc,optstress,ortalg,&
&  paral_kgb,paral_rf,parareel,&
&  pawcpxocc,pawfatbnd,pawlcutd,pawlmix,pawmixdg,pawnhatxc,pawnphi,pawntheta,pawnzlm,pawoptmix,&
&  pawprtden,pawprtdos,pawprtvol,pawprtwf,pawspnorb,pawstgylm,pawusecp,macro_uj,pawujat,pawxcdev,&
&  positron,ppmodel,prepanl,prepgkk,prtacfd,prtbbb,prtcml,prtcs,&
&  prtden,prtdensph,prtdos,prtdosm,prtefg,prteig,prtelf,prtfc,prtfsurf,prtgeo,prtgkk,prtkden,prtkpt,&
&  prtnabla,prtpmp,prtpot,prtspcur,prtstm,prtvha,prtvhxc,prtvol,prtvxc,&
&  prtwant,prtwf,prtxml,prt1dm,ptgroupma,rdmnb,recgratio,recnpath,recnrec,recptrott,rectesteg,restartxf,&
&  rfasr,rfddk,rfelfd,rfmeth,rfmgfd,rfphon,rfstrs,rfthrd,&
&  rfuser,rf1elfd,rf1phon,rf2elfd,rf2phon,rf3elfd,rf3phon,&
&  signperm,smdelta,spgaxor,spgorig,spgroup,spmeth,suskxcrs,symmorphi,symchi,symsigma,&
&  td_mexcit,tfkinfunc,timopt,tl_nprccg,usedmatpu,useexexch,usepaw,&
&  usepawu,userec,useria,userib,useric,userid,userie,usewvl,useylm,vacnum,&
&  wfoptalg,wvl_nprccg,w90iniprj,w90prtunk,w90nplot,xclevel
! Integer arrays
  integer :: bdberry(4),dsifkpt(3),kptrlatt(3,3),ngfft(18),ngfftdg(18),nloalg(5),&
&  qprtrb(3),rfatpol(2),rfdir(3),rf1atpol(2),rf1dir(3),&
&  rf2atpol(2),rf2dir(3),rf3atpol(2),rf3dir(3),scphon_supercell(3),supercell(3),w90cplot(3)
! Integer pointers
  integer, pointer ::  algalch(:)    ! algalch(ntypalch)
  integer, pointer ::  bdgw(:,:)     ! bdgw(2,nkptgw)
  integer, pointer ::  iatfix(:,:)   ! iatfix(3,natom)
  integer, pointer ::  iatsph(:)     ! iatsph(natsph)
  integer, pointer ::  istwfk(:)     ! istwfk(nkpt)
  integer, pointer ::  kberry(:,:)   ! kberry(3,nberry)
  integer, pointer ::  lexexch(:)    ! lexexch(ntypat)
  integer, pointer ::  lpawu(:)      ! lpawu(ntypat)
  integer, pointer ::  nband(:)      ! nband(nkpt*nsppol)
  integer, pointer ::  so_psp(:)     ! so_psp(npsp)
  integer, pointer ::  symafm(:)     ! symafm(nsym)
  integer, pointer ::  symrel(:,:,:) ! symrel(3,3,nsym)
  integer, pointer ::  typat(:)      ! typat(natom)
  integer, pointer ::  w90lplot(:)   ! w90lplot(w90nplot)

! Real
  real(dp) :: alpha,bmass,boxcutmin,bxctmindg,charge,cpus,dedlnn,diecut,diegap,dielam,&
&  dielng,diemac,diemix,diemixmag,dilatmx,dosdeltae,dtion,&
&  ecut,ecuteps,ecutsigx,ecutsm,ecutwfn,effmass,&
&  eshift,exchmix,fband,fixmom,freqremax,freqspmax,freqsusin,freqsuslo,friction,gwencomp,&
&  kptnrm,kptrlen,mdftemp,mditemp,mdwall,nelect,noseinert,&
&  omegasimax,omegasrdmax,pawecutdg,pawovlp,postoldfe,pawujv,&
&  ppmfrq,qptnrm,recrcut,recefermi,rectolden,rhoqpmix,rcut,&
&  sciss,scphon_temp,slabwsrad,slabzbeg,slabzend,soenergy,spbroad,spnorbscl,stmbias,strfact,strprecon,&
&  td_maxene,tfnewton,tl_radius,toldfe,toldff,tolrff,&
&  tolmxf,tolsym,tolvrs,tolwfr,tphysel,tsmear,userra,userrb,userrc,userrd,&
&  userre,vacwidth,vis,vmass,wvl_hgrid,wvl_crmult,wvl_frmult,wvl_cpmult,wvl_fpmult,&
&  zcut
! Types
  type(wvl_internalVars_type) :: wvl
! Real arrays
  real(dp) :: acell_orig(3),angdeg_orig(3),boxcenter(3),&
&  efield(3),genafm(3),qpt(3),qptn(3),rprim_orig(3,3),&
&  rprimd_orig(3,3),strtarget(6),vcutgeo(3),vprtrb(2)
! Real pointers
  real(dp), pointer :: amu(:)            ! amu(ntypat)
  real(dp), pointer :: atvshift(:,:,:)   ! atvshift(16,nsppol,natom)
  real(dp), pointer :: corecs(:)         ! corecs(ntypat)
  real(dp), pointer :: densty(:,:)       ! densty(ntypat,4)
  real(dp), pointer :: dmatpawu(:,:,:,:) ! dmatpawu(2*lpawu+1,2*lpawu+1,nsppol*nspinor,natpu) where natpu=number of atoms with lpawu/=1
  real(dp), pointer :: gw_qlwl(:,:)      ! gw_qlwl(3,gw_nqlwl)
  real(dp), pointer :: jpawu(:)       ! jpawu(ntypat)
  real(dp), pointer :: kpt(:,:)       ! kpt(3,nkpt)
  real(dp), pointer :: kptgw(:,:)     ! kptgw(3,nkptgw)
  real(dp), pointer :: kptns(:,:)     ! kptns(3,nkpt)
  real(dp), pointer :: mixalch(:,:)   ! mixalch(npspalch,ntypalch)
  real(dp), pointer :: occ_orig(:)    ! occ_orig(mband*nkpt*nsppol)
  real(dp), pointer :: ptcharge(:)    ! ptcharge(ntypat)
  real(dp), pointer :: qmass(:)       ! qmass(nnos)
  real(dp), pointer :: qptdm(:,:)     ! qptdm(3,nqptdm)
  real(dp), pointer :: quadmom(:)     ! quadmom(ntypat)
  real(dp), pointer :: ratsph(:)      ! ratsph(ntypat)
  real(dp), pointer :: shiftk(:,:)    ! shiftk(3,nshiftk)
  real(dp), pointer :: spinat(:,:)    ! spinat(3,natom)
  real(dp), pointer :: tnons(:,:)     ! tnons(3,nsym)
  real(dp), pointer :: upawu(:)       ! upawu(ntypat)
  real(dp), pointer :: vel_orig(:,:)  ! vel_orig(3,natom)
  real(dp), pointer :: wtatcon(:,:,:) ! wtatcon(3,natom,nconeq)
  real(dp), pointer :: wtk(:)         ! wtk(nkpt)
  real(dp), pointer :: xred_orig(:,:) ! xred_orig(3,natom)
  real(dp), pointer :: ziontypat(:)   ! ziontypat(ntypat)
  real(dp), pointer :: znucl(:)       ! znucl(npsp)
 end type dataset_type

!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/bandfft_kpt_type
!! NAME
!! bandfft_kpt_type
!!
!! FUNCTION
!! The bandfft_kpt_type structured datatype gather different information
!! about the triple band-fft-kpt parallelisation :
!! tabs which are distributed over all the three dimensions and stored during
!! the calculation, dimensions of messages exchange during the calculations...
!! i.e.: all the informations which were spread over the entire code before and
!! recomputed at each iline, istep or itime STEP with a large probability to
!! make a mistake.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

 type bandfft_kpt_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.



  integer :: flag1_is_allocated             ! determine if the following data are allocated or not
  integer :: npw_tot                        ! array holding the total number of plane waves for each k point
  integer :: ndatarecv                      ! total number of values received by the processor and sended
                                            ! by the other processors band
  integer, pointer :: kg_k_gather(:,:)      ! planewave coordinates
                                            ! (of the processor + sended by other processors band)
  integer, pointer :: recvcounts(:)         ! number of values received by the  processor from each processor band
  integer, pointer :: sendcounts(:)         ! number of values sended   by the  processor to   each processor band
  integer, pointer :: rdispls   (:)         ! positions of values received by the processor from each processor band
  integer, pointer :: sdispls   (:)         ! postions of values sended by the processor to each processor band
  integer, pointer :: gbound(:,:)           ! sphere boundary info: gbound(2*mgfft+8,2)

  integer, pointer :: ind_kg_mpi_to_seq(:)  ! ind_kg_mpi_to_seq(nkpt)
                                            ! in case of //band and //fft, for each processor,
                                            ! index of kg in the numerotation of the sequentiel mode


  integer :: flag2_is_allocated             ! determine if the following data are allocated or not
  real(dp), pointer :: ffnl_gather(:,:,:,:) ! ffnl tab (of the processor + sended by other processors band)
  real(dp), pointer :: kinpw_gather(:)      ! kinpw tab (of the processor + sended by other processors band)
  real(dp), pointer :: ph3d_gather(:,:,:)   ! ph3d tab (of the processor + sended by other processors band)
  real(dp), pointer :: kpg_k_gather(:,:)    ! kpg_k tab (of the processor + sended by other processors band)


  integer :: flag3_is_allocated             ! determine if the following data are allocated or not
  integer :: istwf_k                        ! input option parameter that describes the storage of wfs
  integer :: idatarecv0                     ! position of the planewave coordinates (0,0,0)
  integer :: ndatarecv_tot                  ! total number of received values by the processor
                                            ! (ndatarecv   + number of received opposited planewave coordinates)
  integer :: ndatasend_sym                  ! number of sended values to the processors fft to create opposited
                                            ! planewave coordinates
  integer, pointer :: kg_k_gather_sym(:,:)  ! planewave coordinates
                                            ! (kg_k_gather + opposited planewave coordinates sended by the processors fft)
  integer, pointer :: rdispls_sym(:)        ! positions of values received by the processor from each processor fft
  integer, pointer :: recvcounts_sym(:)     ! number of values received by the  processor from each processor fft
  integer, pointer :: recvcounts_sym_tot(:) ! number of values received by each processor from the  other processors fft
  integer, pointer :: sdispls_sym(:)        ! postions of values sended by the processor to each processor fft
  integer, pointer :: sendcounts_sym(:)     ! number of values sended   by the  processor to each processor fft
  integer, pointer :: sendcounts_sym_all(:) ! number of values sended   by each processor to the other processors fft
  integer, pointer :: tab_proc(:)           ! positions of opposited planewave coordinates in the list of the processors fft

 end type bandfft_kpt_type

!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/MPI_type
!! NAME
!! MPI_type
!!
!! FUNCTION
!! The MPI_type structured datatype gather different information
!! about the MPI parallelisation : number of processors,
!! the index of my processor, the different groups of processors, etc ...
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

 type MPI_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Integer scalar

!***********************************************************************************

!Set of variables for parallelism, that do NOT depend on input variables.
!These are independent of the dataset, and are initialized at the beginning of
!an ABINIT run. The other should be initialized only inside a dataset.


 !TODO Clean the datatype, removing obsolete entries.
 ! Update initmpi_seq since not all the variables are
 ! correctly initialized in order to have a sequential run.
 !
 ! ************************************************
 ! MG coment:
 !  Since there are modules using a private mpi_enreg we cannot
 !  use the declaration to initialize some variables using the sequential value.
 !  Sincerely I would prefer to initialize here these values as it is clearer and easier to read.
 !  For the moment, every time a new kind of flag is added to this
 !  datatype in order to deal with some new sort of parallelism,
 !  please, remember to modify initmpi_seq accordingly
 !
 ! kpgsph was crashing in the GW part due the execution of the following piece of code
 !
 ! if(mpi_enreg%mode_para=='b') then  !this check happened to be true although no parallelism was used
 !  np_band=mpi_enreg%sizecart(2)     !and sizecart was not allocated
 ! end if

  integer :: world_comm       ! number of the world communicator MPI_COMM_WORLD
  integer :: world_group      ! number of the world group of processor MPI_GROUP_NULL
  integer :: me               ! number of my processor in the group of all processors
  integer :: nproc            ! number of processors

  integer :: paral_compil
   ! paral_compil =0 : no -DMPI flag was activated in the compiling procedure
   ! paral_compil =1 : the -DMPI flag was activated in the compiling procedure
  integer :: paral_compil_mpio
   ! paral_compil_mpio =0 : no -DMPIO flag was activated in the compiling procedure
   ! paral_compil_mpio =1 : the -DMPIO flag was activated in the compiling procedure

!***********************************************************************************

! The other should be initialized only inside a dataset.
  integer :: paral_compil_kpt
  integer :: paral_compil_fft

  integer :: paral_level
   ! level of parallelization at a moment in the code
   ! level = 1 : level parareel
   ! level = 2 : level nkpt
   ! level = 3 : level FFT

  integer :: paralbd
   ! relevant only if paral_compil_kpt=1 . So, in addition to the kpt parallelization :
   ! paralbd=0 : (no //ization on bands)
   ! paralbd=1 : (//ization on bands)
   ! paralbd>1 : (//ization on blocks of bands)

  integer :: gwpara
   ! level of parallelization at a moment in the GW code
   ! level = 0 : no parallelization (seq run)
   ! level = 1 : kpoints
   ! level = 2 : bands
   ! level = 3 : in the future, maybe : mixed (kpoints+bands)

  integer :: me_group         ! number of my processor in my group of kpt
  integer :: nproc_group      ! number of processors in my group of kpt
  integer :: me_fft           ! number of my processor in my group of FFT
  integer :: me_band           ! number of my processor in my group of bands
  integer :: nproc_fft        ! number of processors in my group of FFT
  integer :: master_fft       ! number of master of my fft group (in the world_group)
  integer :: paral_fft        ! set to 1 if the FFT parallelisation is active
  integer :: me_g0            ! if set to 1, means that the current processor is taking care of the G(0 0 0) planewave.
  integer :: num_group_fft    ! number of FFT group of my processor. 0 if my processor is not in a group
  integer :: num_group        ! number of group of my processor. 0 if my processor is not in a group
  integer :: nproc_per_kpt    ! number of processors per kpt

  integer :: fft_master_group
   ! fft_master_group
   ! group of processors of fft_master_comm
   ! exists only when paral_fft = 1

  integer :: fft_master_comm
   ! fft_master_comm
   ! communicator on master processors
   ! (one processor per fft_group or all processors when paral_fft = 0)

  integer :: fft_option_lob
   ! fft_option_lob
   ! option for lob
   ! fft_option_lob=1 : old version of lob
   ! fft_option_lob=2 : new version of lob
   ! exists only when paral_fft = 1

  integer :: has_band_comm
   ! has_band_comm
   ! 1 if mpi_enreg%band_comm(:) is allocated

! Integer arrays

  integer, pointer :: band_comm(:)
   ! band_comm(nproc_per_kpt)
   ! tab of communicators of processors which treat one set ogf bands
   ! exists only when paralbd = 1 and has_band_comm=1

  integer, pointer :: fft_group(:)
   ! fft_group(nkpt*nsppol)
   ! tab of groups of processors which treat ffts
   ! exists only when paral_fft = 1

  integer, pointer :: fft_comm(:)
   ! fft_comm(nkpt*nsppol)
   ! tab of communicators of processors which treat ffts of a kpt
   ! exists only when paral_fft = 1

  integer, pointer :: proc_distrb(:,:,:)
   ! proc_distrb(nkpt,mband,nsppol)
   ! number of the processor that will treat
   ! each band in each k point.

  integer, pointer :: kpt_group(:)
   ! kpt_group(nproc_per_kpt)
   ! tab of groups of processors which treat one nkpt/nsppol
   ! exists only when paralbd > 1

  integer, pointer :: kpt_comm(:)
   ! kpt_comm(nproc_per_kpt)
   ! tab of communicators of processors which treat one nkpt/nsppol
   ! exists only when paralbd > 1

  integer, pointer :: kptdstrb(:,:,:)
   ! kptdstrb(me,ineigh,ikptloc)
   ! tab of processors required for mv_3dte.f and berryphase_new.f

  integer, pointer :: kptdstrbi(:,:,:)
   ! same as kptdstrb, but for k-points in the iBZ
   ! required for MPI // of the finite electric field (see vtorho.f)

  integer, pointer :: nplanes_fft(:)
   ! nplanes_fft(nkpt)
   ! number of planes for my proc me_fft
   ! exists only if mpi_enreg%paral_compil_fft==1

  integer, pointer :: ind_fft_planes(:,:)
   ! ind_fft_planes(nkpt,nplanes_fft)
   ! indice of planes for each kpoint for my proc me_fft
   ! exists only if mpi_enreg%paral_compil_fft==1

  integer           :: flag_ind_kg_mpi_to_seq
   ! flag to activate the building of bandfft_kpt(:)%ind_kg_mpi_to_seq

  type(bandfft_kpt_type), pointer :: bandfft_kpt(:)
   ! bandfft_kpt(mkmem)
   ! Contains all the informations related to the band/FFT parallelism which depends on kpt
   ! exists only if mpi_enreg%paral_kgb==1

  integer, pointer :: tab_kpt_distrib(:)
   ! tab_kpt_distrib(nkpt)
   ! Indicates the correspondence between the ikpt and ikpt_this_proc
   ! exists only if mpi_enreg%paral_kgb==1



! Adds for parallelization over perturbations
  integer :: paral_compil_respfn
   ! paral_compil_respfn =0 : no -DMPI flag was activated in the compiling procedure
   ! paral_compil_respfn =1 : the -DMPI flag was activated in the compiling procedure

  integer :: me_respfn           ! number of my processor in my group of perturbations
  integer :: nproc_respfn        ! number of processors in my group of perturbations
  integer :: my_respfn_group     ! my group for calculating perturbations
  integer :: my_respfn_comm      ! my communicator of my_respfn_group
  integer :: respfn_master_group ! groups for masters of respfn_groups
  integer :: respfn_master_comm  ! communicator for masters of respfn_groups
  integer :: ngroup_respfn       ! number of groups for calculating perturbations
  integer :: spaceComm           ! communicator for calculating responsefunction
                                 ! default is MPI_COMM_WORLD but may be changed in 08seqpar/loper3.F90

  integer, pointer :: respfn_group(:) ! groups for calculating perturbations
  integer, pointer :: respfn_comm(:)  ! communicators for respfn_group

  ! Wavelet paralelisation, use when %paral_compil_fft == 1
  ! Array to store the description of the scaterring in real space of
  ! the potentials and density. It is allocated to (0:nproc-1,4).
  ! The four values are:
  ! - the density size in z direction ( = ngfft(3)) ;
  ! - the potential size in z direction ( <= ngfft(3)) ;
  ! - the position of the first value in the complete array ;
  ! - the shift for the potential in the array.
  integer, pointer :: nscatterarr(:,:)
  ! Array to store the total size (of this proc) of the potentails arrays when
  ! the memory is distributed following nscatterarr.
  integer, pointer :: ngatherarr(:,:)
  ! Store the ionic potential size in z direction.
  integer :: ngfft3_ionic
  ! End wavelet additions

!This is for the bandFFT case
   character :: mode_para
   !If mode_para=='bandFFT', we are in bandFFT mode
   integer :: commcart
   !This is the communicator for the full cartesian array
   integer :: comm_band, comm_fft
   !The communicators over bands and fft respectively
   integer :: me_cart
   !This is the rank of the proc in the full cartesian array
   integer :: dimcart
   !This is the dimension of the cartesian array (2 for 2-dim)
   integer :: nproc_band
   !This is the number of procs on which we distribute bands
   integer, pointer :: sizecart(:)
   !The first dimension is the number of fft processors, the second the number of bands
   integer, pointer :: coords(:)
   !The coordinate of the proc in the cartesian array

!This is for the kpt & bandFFt case
   integer :: commcart_3d      ! 3D communicator
   integer :: comm_kpt         ! communicator of kpt
   integer :: me_kpt           ! number of my processor in my group of kpt
   integer :: nproc_kpt        ! number of procs on which we distribute kpt
   integer :: me_cart_2d       ! This is the rank of the proc in the commcart

! Adds for parareel
  integer :: parareel
   ! parareel = 0 default
   ! parareel = 1 if treats parareel case

! All the following data exist only in the parareel=1 case
  integer :: npara                 ! number of loops on gstate
  integer :: ipara                 ! number of actual internal loop on gstate
  integer :: jpara                 ! number of actual external loop on gstate
  integer :: me_group_para         ! number of my processor in my group of para
  integer :: nproc_group_para      ! number of processors in my group of para
  integer :: num_group_para        ! number of group of my processor. 0 if my processor is not in a group
  integer :: nproc_per_para        ! number of processors per para
  integer :: master_group_para     ! number of the master processor (in the world group) of my group of para

  integer, pointer :: proc_distrb_para(:,:)
   ! proc_distrb_para(npara,nkpt)
   ! exists only when parareel = 1
   ! number of the processor that will treat
   ! each kpt in each para.

  integer, pointer :: kpt_group_para(:)
   ! kpt_group_para(npara)
   ! tab of groups of processors which treat one npara
   ! exists only when parareel = 1

  integer, pointer :: kpt_comm_para(:)
   ! kpt_comm_para(npara)
   ! tab of communicators of processors which treat one npara
   ! exists only when parareel = 1

  integer :: bandpp

! This is for printing the processor distribution in invars1m

  integer, pointer :: keywp(:,:)
   ! key word for parallelisation keywp(5,250)
   ! The first index gives: nproc, nkpt, npband, npfft and bandpp
   ! for each possible choice among 250, for eack dataset.


  integer, pointer :: trialproc(:)
   ! trialproc(2)
   ! gives for each dataset the number of trial processor and if this number is suitable for the calculation

 end type MPI_type

end module defs_abitypes
!!***

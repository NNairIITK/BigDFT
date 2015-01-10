
!!****m* ABINIT/defs_abitypes
!! NAME
!! defs_abitypes
!!
!! FUNCTION
!! This module contains definitions of high-level structured datatypes for the
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
!!  (5) Declare variables on separated lines in order to reduce the occurence of bzr conflicts.
!!
!! List of datatypes :
!! * aim_dataset_type : the "dataset" for aim
!! * anaddb_dataset_type : the "dataset" for anaddb
!! * dataset_type : the "dataset" for the main abinit code
!! * bandfft_kpt_type : the "dataset" for triple band-fft-kpt parallelization
!! * MPI_type : the data related to MPI parallelization
!!
!! COPYRIGHT
!! Copyright (C) 2001-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
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

 type aim_dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.

! Since all these input variables are described in the aim_help.html
! file, they are not described in length here ...

! Integer
  integer :: crit
  integer :: denout
  integer :: dltyp
  integer :: gpsurf
  integer :: irho
  integer :: ivol
  integer :: lapout
  integer :: nsa
  integer :: nsb
  integer :: nsc

  integer :: batom  !! Warning : corresponds to the input variable atom
  integer :: foll   !! Warning : corresponds to the input variable follow
  integer :: isurf  !! Warning : corresponds to the input variable surf
  integer :: irsur  !! Warning : corresponds to the input variable rsurf
  integer :: nph    !! Warning : corresponds to the input variable nphi
  integer :: npt    !! Warning : corresponds to the input variable inpt
  integer :: nth    !! Warning : corresponds to the input variable ntheta
  integer :: plden  !! Warning : not documented in help file ?!

  integer :: ngrid(3)

! Real
  real(dp) :: atrad
  real(dp) :: coff1 
  real(dp) :: coff2
  real(dp) :: dpclim
  real(dp) :: folstp
  real(dp) :: lgrad
  real(dp) :: lgrad2
  real(dp) :: lstep
  real(dp) :: lstep2
  real(dp) :: maxatd 
  real(dp) :: maxcpd
  real(dp) :: phimax 
  real(dp) :: phimin
  
  real(dp) :: dr0    !! Warning : correspond to the input variable radstp
  real(dp) :: phi0   !! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   !! Warning : correspond to the input variable ratmin
  real(dp) :: th0    !! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax !! Warning : correspond to the input variable thetamax
  real(dp) :: themin !! Warning : correspond to the input variable thetamin

  real(dp) :: foldep(3)
  real(dp) :: scal(3)
  real(dp) :: vpts(3,4)

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

 type anaddb_dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.


! Since all these input variables are described in the anaddb_help.html
! file, they are not described in length here ...
! Integer
  integer :: alphon
  integer :: asr
  integer :: brav
  integer :: chneut
  integer :: dieflag
  integer :: dipdip
  integer :: ep_scalprod
  integer :: eivec
  integer :: elaflag
  integer :: elphflag
  integer :: enunit
  integer :: gkk2write 
  integer :: gkk_rptwrite
  integer :: gkqwrite
  integer :: iavfrq
  integer :: ifcana
  integer :: ifcflag
  integer :: ifcout
  integer :: ifltransport
  integer :: instrflag
  integer :: natfix
  integer :: natifc
  integer :: natom
  integer :: nchan
  integer :: nfreq
  integer :: ngrids
  integer :: nlflag
  integer :: nph1l
  integer :: nph2l
  integer :: nqpath
  integer :: nqshft
  integer :: nsphere
  integer :: nstrfix
  integer :: ntemper
  integer :: nwchan
  integer :: piezoflag
  integer :: polflag
  integer :: prtdos
  integer :: prtmbm
  integer :: prtfsurf
  integer :: prtnest
  integer :: ramansr
  integer :: relaxat
  integer :: relaxstr
  integer :: rfmeth
  integer :: selectz
  integer :: symdynmat
  integer :: telphint
  integer :: thmflag
  integer :: qgrid_type
  integer :: ep_b_min
  integer :: ep_b_max
  integer :: ep_keepbands
  integer :: ep_nqpt
  integer :: ep_prt_yambo
  integer :: symgkq

  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: istrfix(6)
  integer :: ng2qpt(3)
  integer :: kptrlatt(3,3)

! Real(dp)
  real(dp) :: a2fsmear
  real(dp) :: dosdeltae
  real(dp) :: dossmear
  real(dp) :: dostol 
  real(dp) :: elphsmear
  real(dp) :: elph_fermie
  real(dp) :: frmax
  real(dp) :: frmin
  real(dp) :: temperinc
  real(dp) :: tempermin
  real(dp) :: thmtol
  real(dp) :: mustar
  real(dp) :: rifcsph

  real(dp) :: q1shft(3,4)
  real(dp) :: q2shft(3)
  real(dp) :: targetpol(3)

! Integer pointers
  integer, pointer :: atifc(:)    
   ! atifc(natom) WARNING : there is a transformation of this input variable, in chkin9
   ! This should be changed ...

  integer, pointer :: iatfix(:)   
  ! iatfix(natom)

! Real pointers
  real(dp), pointer :: qnrml1(:)  
  ! qnrml1(nph1l)

  real(dp), pointer :: qnrml2(:)  
  ! qnrml2(nph2l)

  real(dp), pointer :: qpath(:,:) 
  ! qpath(3,nqpath)

  real(dp), pointer :: qph1l(:,:) 
  ! qph1l(3,nph1l)

  real(dp), pointer :: qph2l(:,:) 
  ! qph2l(3,nph2l)

  real(dp), pointer :: ep_qptlist(:,:) 
  ! qph2l(3,ep_nqpt)

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

type dataset_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.

! Since all these input variables are described in the abinis_help.html
! file, they are not described in length here ...

! Integer
 integer :: accesswff
 integer :: awtr
 integer :: bandpp
 integer :: bdeigrf
 integer :: berryopt
 integer :: brvltt
 integer :: chkexit
 integer :: chkprim
 integer :: delayperm
 integer :: dmatpuopt
 integer :: dmatudiag
 integer :: dmftbandi
 integer :: dmftbandf
 integer :: enunit
 integer :: exchn2n3d
 integer :: fftgw
 integer :: fft_opt_lob
 integer :: frzfermi
 integer :: getcell
 integer :: getddk
 integer :: getden
 integer :: getkss
 integer :: getocc
 integer :: getqps
 integer :: getscr
 integer :: getsuscep
 integer :: getvel
 integer :: getwfk
 integer :: getwfq
 integer :: getxcart
 integer :: getxred
 integer :: get1den
 integer :: get1wf
 integer :: getbseig
 integer :: gw_optimal_level
 integer :: gw_nband_optimal
 integer :: gwcalctyp
 integer :: gwcomp
 integer :: gwgamma
 integer :: gw_nqlwl
 integer :: gw_nstep
 integer :: gw_sigxcore
 integer :: gwmem
 integer :: gwpara
 integer :: gw_sctype
 integer :: iboxcut
 integer :: icoulomb
 integer :: icutcoul
 integer :: idyson
 integer :: ieig2rf
 integer :: iextrapwf
 integer :: ikhxc
 integer :: imgmov
 integer :: inclvkb
 integer :: intexact
 integer :: intxc
 integer :: ionmov
 integer :: iprcch
 integer :: iprcel
 integer :: iprctfvw
 integer :: iprcfc
 integer :: irdddk
 integer :: irdden
 integer :: irdkss
 integer :: irdqps
 integer :: irdscr
 integer :: irdsuscep
 integer :: irdwfk
 integer :: irdwfq
 integer :: ird1wf
 integer :: irdbseig
 integer :: iscf
 integer :: isecur
 integer :: istatr
 integer :: istatshft
 integer :: ixc
 integer :: ixcpositron
 integer :: jdtset !  jdtset contains the actual number of the dataset
 integer :: jellslab
 integer :: kptopt
 integer :: kssform
 integer :: ldgapp
 integer :: localrdwf
 integer :: maxnsym
 integer :: mband
 integer :: mffmem
 integer :: mgfft
 integer :: mgfftdg
 integer :: mkmem
 integer :: mkqmem
 integer :: mk1mem
 integer :: nnos
 integer :: mpw
 integer :: mqgrid
 integer :: mqgriddg
 integer :: natom
 integer :: natpawu
 integer :: natrd
 integer :: natsph
 integer :: natvshift
 integer :: nbandkss
 integer :: nbandsus
 integer :: nbdblock
 integer :: nbdbuf
 integer :: nberry
 integer :: nconeq
 integer :: nctime
 integer :: ndtset
 integer :: ndynimage
 integer :: ndyson
 integer :: nfft
 integer :: nfftdg
 integer :: nfreqim
 integer :: nfreqre
 integer :: nfreqsp
 integer :: nfreqsus
 integer :: ngeohist
 integer :: ngroup_rf
 integer :: nimage
 integer :: nkptgw
 integer :: nkpt
 integer :: nline
 integer :: nnsclo
 integer :: nomegasf
 integer :: nomegasi
 integer :: nomegasrd
 integer :: npband
 integer :: npfft
 integer :: npimage
 integer :: npkpt
 integer :: npsp
 integer :: npspalch
 integer :: npulayit
 integer :: npweps
 integer :: npwkss
 integer :: npwsigx
 integer :: npwwfn
 integer :: nqpt
 integer :: nqptdm
 integer :: nscforder
 integer :: nsheps
 integer :: nshiftk
 integer :: nshsigx
 integer :: nshwfn
 integer :: nspden
 integer :: nspinor
 integer :: nsppol
 integer :: nstep
 integer :: nsym
 integer :: ntime
 integer :: ntimimage
 integer :: ntypalch
 integer :: ntypat
 integer :: ntyppure
 integer :: nwfshist
 integer :: occopt
 integer :: optcell
 integer :: optdriver
 integer :: optforces
 integer :: optfreqsus
 integer :: optnlxccc
 integer :: optstress
 integer :: ortalg
 integer :: paral_kgb
 integer :: paral_rf
 integer :: pawcpxocc
 integer :: pawfatbnd
 integer :: pawlcutd
 integer :: pawlmix
 integer :: pawmixdg
 integer :: pawnhatxc
 integer :: pawnphi
 integer :: pawntheta
 integer :: pawnzlm
 integer :: pawoptmix
 integer :: pawprtden
 integer :: pawprtdos
 integer :: pawprtvol
 integer :: pawprtwf
 integer :: pawspnorb
 integer :: pawstgylm
 integer :: pawusecp
 integer :: macro_uj
 integer :: pawujat
 integer :: pawxcdev
 integer :: positron
 integer :: posnstep
 integer :: ppmodel
 integer :: prepanl
 integer :: prepgkk
 integer :: prtbbb
 integer :: prtcml
 integer :: prtcs
 integer :: prtden
 integer :: prtdensph
 integer :: prtdos
 integer :: prtdosm
 integer :: prtefg
 integer :: prteig
 integer :: prtelf
 integer :: prtfc
 integer :: prtfsurf
 integer :: prtgden
 integer :: prtgeo
 integer :: prtgkk
 integer :: prtkden
 integer :: prtkpt
 integer :: prtlden
 integer :: prtnabla
 integer :: prtpmp
 integer :: prtpot
 integer :: prtspcur
 integer :: prtstm
 integer :: prtvha
 integer :: prtvhxc
 integer :: prtvol
 integer :: prtvxc
 integer :: prtwant
 integer :: prtwf
 integer :: prtxangst
 integer :: prtxcart
 integer :: prtxml
 integer :: prtxred
 integer :: prt1dm
 integer :: ptgroupma
 integer :: rdmnb
 integer :: recgratio
 integer :: recnpath
 integer :: recnrec
 integer :: recptrott
 integer :: rectesteg
 integer :: restartxf
 integer :: rfasr
 integer :: rfddk
 integer :: rfelfd
 integer :: rfmeth
 integer :: rfmgfd
 integer :: rfphon
 integer :: rfstrs
 integer :: rfuser
 integer :: rf1elfd
 integer :: rf1phon
 integer :: rf2elfd
 integer :: rf2phon
 integer :: rf3elfd
 integer :: rf3phon
 integer :: signperm
 integer :: smdelta
 integer :: spgaxor
 integer :: spgorig
 integer :: spgroup
 integer :: spmeth
 integer :: suskxcrs
 integer :: symmorphi
 integer :: symchi
 integer :: symsigma
 integer :: td_mexcit
 integer :: tfkinfunc
 integer :: timopt
 integer :: tl_nprccg
 integer :: usedmatpu
 integer :: usedmft
 integer :: useexexch
 integer :: usekden
 integer :: usepaw
 integer :: usepawu
 integer :: userec
 integer :: useria
 integer :: userib
 integer :: useric
 integer :: userid
 integer :: userie
 integer :: usewvl
 integer :: useylm
 integer :: vacnum
 integer :: wfoptalg
 integer :: wvl_nprccg
 integer :: w90iniprj
 integer :: w90prtunk
 integer :: xclevel

!Integer arrays
 integer :: bdberry(4)
 integer :: kptrlatt(3,3)
 integer :: ngfft(18)
 integer :: ngfftdg(18)
 integer :: nloalg(5)
 integer :: qprtrb(3)
 integer :: rfatpol(2)
 integer :: rfdir(3)
 integer :: rf1atpol(2)
 integer :: rf1dir(3)
 integer :: rf2atpol(2)
 integer :: rf2dir(3)
 integer :: rf3atpol(2)
 integer :: rf3dir(3)
 integer :: scphon_supercell(3)
 integer :: supercell(3)

!Integer pointers
 integer, pointer ::  algalch(:)    ! algalch(ntypalch)
 integer, pointer ::  bdgw(:,:)     ! bdgw(2,nkptgw)
 integer, pointer ::  dynimage(:)   ! dynimage(nimage or mxnimage)
 integer, pointer ::  iatfix(:,:)   ! iatfix(3,natom)
 integer, pointer ::  iatsph(:)     ! iatsph(natsph)
 integer, pointer ::  istwfk(:)     ! istwfk(nkpt)
 integer, pointer ::  kberry(:,:)   ! kberry(3,nberry)
 integer, pointer ::  lexexch(:)    ! lexexch(ntypat)
 integer, pointer ::  lpawu(:)      ! lpawu(ntypat)
 integer, pointer ::  nband(:)      ! nband(nkpt*nsppol)
 integer, pointer ::  normpawu(:)   ! normpawu(ntypat)
 integer, pointer ::  so_psp(:)     ! so_psp(npsp)
 integer, pointer ::  symafm(:)     ! symafm(nsym)
 integer, pointer ::  symrel(:,:,:) ! symrel(3,3,nsym)
 integer, pointer ::  typat(:)      ! typat(natom)

!Real
 real(dp) :: alpha
 real(dp) :: bmass
 real(dp) :: boxcutmin
 real(dp) :: bxctmindg
 real(dp) :: charge
 real(dp) :: cpus
 real(dp) :: diecut
 real(dp) :: diegap
 real(dp) :: dielam
 real(dp) :: dielng
 real(dp) :: diemac
 real(dp) :: diemix
 real(dp) :: diemixmag
 real(dp) :: dilatmx
 real(dp) :: dosdeltae
 real(dp) :: dtion
 real(dp) :: ecut
 real(dp) :: ecuteps
 real(dp) :: ecutsigx
 real(dp) :: ecutsm
 real(dp) :: ecutwfn
 real(dp) :: effmass
 real(dp) :: eshift
 real(dp) :: exchmix
 real(dp) :: fband
 real(dp) :: fixmom
 real(dp) :: freqremax
 real(dp) :: freqspmax
 real(dp) :: freqsusin
 real(dp) :: freqsuslo
 real(dp) :: friction
 real(dp) :: fxcartfactor
 real(dp) :: gwencomp
 real(dp) :: gw_toldfeig
 real(dp) :: kptnrm
 real(dp) :: kptrlen
 real(dp) :: mdftemp
 real(dp) :: mditemp
 real(dp) :: mdwall
 real(dp) :: nelect
 real(dp) :: noseinert
 real(dp) :: omegasimax
 real(dp) :: omegasrdmax
 real(dp) :: pawecutdg
 real(dp) :: pawovlp
 real(dp) :: posocc
 real(dp) :: postoldfe
 real(dp) :: pawujv
 real(dp) :: ppmfrq
 real(dp) :: qptnrm
 real(dp) :: recrcut
 real(dp) :: recefermi
 real(dp) :: rectolden
 real(dp) :: rhoqpmix
 real(dp) :: rcut
 real(dp) :: sciss
 real(dp) :: scphon_temp
 real(dp) :: slabwsrad
 real(dp) :: slabzbeg
 real(dp) :: slabzend
 real(dp) :: soenergy
 real(dp) :: spbroad
 real(dp) :: spnorbscl
 real(dp) :: stmbias
 real(dp) :: strfact
 real(dp) :: strprecon
 real(dp) :: td_maxene
 real(dp) :: tl_radius
 real(dp) :: toldfe
 real(dp) :: toldff
 real(dp) :: tolimg
 real(dp) :: tolmxf
 real(dp) :: tolrff
 real(dp) :: tolsym
 real(dp) :: tolvrs
 real(dp) :: tolwfr
 real(dp) :: tphysel
 real(dp) :: tsmear
 real(dp) :: userra
 real(dp) :: userrb
 real(dp) :: userrc
 real(dp) :: userrd
 real(dp) :: userre
 real(dp) :: vacwidth
 real(dp) :: vis
 real(dp) :: wvl_hgrid
 real(dp) :: wvl_crmult
 real(dp) :: wvl_frmult
 real(dp) :: wvl_cpmult
 real(dp) :: wvl_fpmult
 real(dp) :: zcut

!Real arrays
 real(dp) :: boxcenter(3)
 real(dp) :: efield(3) 
 real(dp) :: genafm(3) 
 real(dp) :: qpt(3) 
 real(dp) :: qptn(3)
 real(dp) :: strtarget(6) 
 real(dp) :: vcutgeo(3)
 real(dp) :: vprtrb(2)

!Real pointers
 real(dp), pointer :: acell_orig(:,:)     ! acell_orig(3,nimage)
 real(dp), pointer :: amu(:)              ! amu(ntypat)
 real(dp), pointer :: atvshift(:,:,:)     ! atvshift(16,nsppol,natom)
 real(dp), pointer :: corecs(:)           ! corecs(ntypat)
 real(dp), pointer :: densty(:,:)         ! densty(ntypat,4)
 real(dp), pointer :: dmatpawu(:,:,:,:)   ! dmatpawu(2*lpawu+1,2*lpawu+1,nsppol*nspinor,natpu) 
                                          !   where natpu=number of atoms with lpawu/=1
 real(dp), pointer :: gw_qlwl(:,:)        ! gw_qlwl(3,gw_nqlwl)
 real(dp), pointer :: jpawu(:)            ! jpawu(ntypat)
 real(dp), pointer :: kpt(:,:)            ! kpt(3,nkpt)
 real(dp), pointer :: kptgw(:,:)          ! kptgw(3,nkptgw)

 real(dp), pointer :: kptns(:,:)          
 ! kptns(3,nkpt)
 ! k-points renormalized and shifted. The ones that should be used inside the code.

 real(dp), pointer :: mixalch(:,:)        ! mixalch(npspalch,ntypalch)
 real(dp), pointer :: occ_orig(:)         ! occ_orig(mband*nkpt*nsppol)
 real(dp), pointer :: ptcharge(:)         ! ptcharge(ntypat)
 real(dp), pointer :: qmass(:)            ! qmass(nnos)
 real(dp), pointer :: qptdm(:,:)          ! qptdm(3,nqptdm)
 real(dp), pointer :: quadmom(:)          ! quadmom(ntypat)
 real(dp), pointer :: ratsph(:)           ! ratsph(ntypat)
 real(dp), pointer :: rprim_orig(:,:,:)   ! rprim_orig(3,3,nimage)
 real(dp), pointer :: rprimd_orig(:,:,:)  ! rprimd_orig(3,3,nimage)
 real(dp), pointer :: shiftk(:,:)         ! shiftk(3,nshiftk)
 real(dp), pointer :: spinat(:,:)         ! spinat(3,natom)
 real(dp), pointer :: tnons(:,:)          ! tnons(3,nsym)
 real(dp), pointer :: upawu(:)            ! upawu(ntypat)
 real(dp), pointer :: vel_orig(:,:,:)     ! vel_orig(3,natom,nimage)
 real(dp), pointer :: wtatcon(:,:,:)      ! wtatcon(3,natom,nconeq)
 real(dp), pointer :: wtk(:)              ! wtk(nkpt)
 real(dp), pointer :: xred_orig(:,:,:)    ! xred_orig(3,natom,nimage)
 real(dp), pointer :: ziontypat(:)        ! ziontypat(ntypat)
 real(dp), pointer :: znucl(:)            ! znucl(npsp)


!BEGIN VARIABLES FOR @Bethe-Salpeter
 integer :: bs_algorithm
 integer :: bs_haydock_niter
 integer :: bs_exchange_term 
 integer :: bs_coulomb_term
 integer :: bs_calctype
 integer :: bs_coupling

 real(dp) :: bs_haydock_tol

 integer :: bs_eh_basis_set(2)

 real(dp) :: bs_eh_cutoff(3)
 real(dp) :: bs_freq_mesh(3)
!END VARIABLES FOR @Bethe-Salpeter.

 ! Types
 type(wvl_internalVars_type) :: wvl

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

 type bandfft_kpt_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: flag1_is_allocated             ! determine if the following data are allocated or not
  integer :: npw_tot                        ! array holding the total number of plane waves for each k point
  integer :: ndatarecv                      ! total number of values received by the processor and sent
                                            ! by the other processors band
  integer, pointer :: kg_k_gather(:,:)      ! planewave coordinates
                                            ! (of the processor + sent by other processors band)
  integer, pointer :: recvcounts(:)         ! number of values received by the  processor from each processor band
  integer, pointer :: sendcounts(:)         ! number of values sent   by the  processor to   each processor band
  integer, pointer :: rdispls   (:)         ! positions of values received by the processor from each processor band
  integer, pointer :: sdispls   (:)         ! postions of values sent by the processor to each processor band
  integer, pointer :: gbound(:,:)           ! sphere boundary info: gbound(2*mgfft+8,2)

  integer, pointer :: ind_kg_mpi_to_seq(:)  ! ind_kg_mpi_to_seq(nkpt)
                                            ! in case of //band and //fft, for each processor,
                                            ! index of kg in the numerotation of the sequentiel mode


  integer :: flag2_is_allocated             ! determine if the following data are allocated or not
  real(dp), pointer :: ffnl_gather(:,:,:,:) ! ffnl tab (of the processor + sent by other processors band)
  real(dp), pointer :: kinpw_gather(:)      ! kinpw tab (of the processor + sent by other processors band)
  real(dp), pointer :: ph3d_gather(:,:,:)   ! ph3d tab (of the processor + sent by other processors band)
  real(dp), pointer :: kpg_k_gather(:,:)    ! kpg_k tab (of the processor + sent by other processors band)


  integer :: flag3_is_allocated             ! determine if the following data are allocated or not
  integer :: istwf_k                        ! input option parameter that describes the storage of wfs
  integer :: idatarecv0                     ! position of the planewave coordinates (0,0,0)
  integer :: ndatarecv_tot                  ! total number of received values by the processor
                                            ! (ndatarecv   + number of received opposited planewave coordinates)
  integer :: ndatasend_sym                  ! number of sent values to the processors fft to create opposited
                                            ! planewave coordinates
  integer, pointer :: kg_k_gather_sym(:,:)  ! planewave coordinates
                                            ! (kg_k_gather + opposited planewave coordinates sent by the processors fft)
  integer, pointer :: rdispls_sym(:)        ! positions of values received by the processor from each processor fft
  integer, pointer :: recvcounts_sym(:)     ! number of values received by the  processor from each processor fft
  integer, pointer :: recvcounts_sym_tot(:) ! number of values received by each processor from the  other processors fft
  integer, pointer :: sdispls_sym(:)        ! postions of values sent by the processor to each processor fft
  integer, pointer :: sendcounts_sym(:)     ! number of values sent   by the  processor to each processor fft
  integer, pointer :: sendcounts_sym_all(:) ! number of values sent   by each processor to the other processors fft
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

 type MPI_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.
! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.

! *****************************************************************************************
! Please make sure that initmpi_seq is changed so that any variable or any flag in MPI_type 
! is initialized with the value used for sequential executions.
! In particular any MPI communicator should be set to MPI_COMM_SELF 
! ************************************************************************************

  ! Set of variables for parallelism, that do NOT depend on input variables.
  ! These are independent of the dataset, and are initialized at the beginning of
  ! an ABINIT run. The other should be initialized only inside a dataset.

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
   ! level = 1 : unused (previously used for parareel)
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

  integer :: me_group         ! number of my processor in my group of kpt
  integer :: nproc_group      ! number of processors in my group of kpt
  integer :: me_fft           ! number of my processor in my group of FFT
  integer :: me_band          ! number of my processor in my group of bands
  integer :: nproc_fft        ! number of processors in my group of FFT
  integer :: master_fft       ! number of master of my fft group (in the world_group)
  integer :: paral_fft        ! set to 1 if the FFT parallelisation is active
  integer :: me_g0            ! if set to 1, means that the current processor is taking care of the G(0 0 0) planewave.
  integer :: num_group_fft    ! number of FFT group of my processor. 0 if my processor is not in a group
  integer :: num_group        ! number of group of my processor. 0 if my processor is not in a group
  integer :: nproc_per_kpt    ! number of processors per kpt

  integer :: fft_master_group
   ! group of processors of fft_master_comm
   ! exists only when paral_fft = 1

  integer :: fft_master_comm
   ! communicator on master processors
   ! (one processor per fft_group or all processors when paral_fft = 0)

  integer :: fft_option_lob
   ! option for lob
   ! fft_option_lob=1 : old version of lob
   ! fft_option_lob=2 : new version of lob
   ! exists only when paral_fft = 1

  integer :: has_band_comm
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

  integer :: flag_ind_kg_mpi_to_seq
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

  !MG: TODO introduce comm_respfn.
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
   ! If mode_para=='bandFFT', we are in bandFFT mode

   integer :: commcart
   ! This is the communicator for the full cartesian array

   integer :: comm_band, comm_fft
   ! The communicators over bands and fft respectively

   integer :: me_cart
   ! This is the rank of the proc in the full cartesian array

   integer :: dimcart
   ! This is the dimension of the cartesian array (2 for 2-dim)

   integer :: nproc_band
   ! This is the number of procs on which we distribute bands

   integer, pointer :: sizecart(:)
   ! The first dimension is the number of fft processors, the second the number of bands

   integer, pointer :: coords(:)
   ! The coordinate of the proc in the cartesian array

!This is for the kpt & bandFFt case
   integer :: commcart_3d      ! 3D communicator
   integer :: comm_kpt         ! communicator of kpt
   integer :: me_kpt           ! number of my processor in my group of kpt
   integer :: nproc_kpt        ! number of procs on which we distribute kpt
   integer :: me_cart_2d       ! This is the rank of the proc in the commcart

!This is for the parallelisation over atoms (PAW)
   integer :: nproc_atom             ! Number of processors on which we distribute atoms
   integer :: natom                  ! Number of atoms treated by current processor
   integer :: comm_atom              ! Communicator over atoms
   integer,pointer :: atom_indx(:)   ! atom_indx(mpi_enreg%natom)= indexes of the atoms treated by current processor

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

!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/datafiles_type
!! NAME
!! datafiles_type
!!
!! FUNCTION
!! The datafiles_type structures datatype gather all the variables related
!! to files, such as filename, and file units.
!! For one dataset, it is initialized in driver.F90, and will not change
!! at all during the treatment of the dataset.
!!
!! SOURCE

 type datafiles_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


  integer :: ireadden
   ! ireadden non-zero  if the den file must be read

  integer :: ireadkden
   ! ireadkden non-zero  if the kden file must be read

  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffk file must be read
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1 file must be read
   !   (if ird1wf non-zero or get1wf non-zero)

  integer :: unchi0  ! unit number for chi0 files
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: unddk   ! unit number for ddk 1WF file
  integer :: unem1ggp! unit number for Epsilon minus one (G,Gp) file
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unkss   ! unit number for KSS file
  integer :: unqps   ! unit number for QPS file
  integer :: unscr   ! unit number for SCR file
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: unpaw   ! unit number for temporary PAW data (for ex. rhoij_nk) (Paw only)
  integer :: unpaw1  ! unit number for temporary PAW first-order cprj1=<c1_k,q|p>(1) data
  integer :: unpawq  ! unit number for temporary PAW cprjq=<c+_k+q|p> at k+qdata
  integer :: unpos   ! unit number for restart molecular dynamics

  character(len=fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input :
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character(len=fnlen) :: filchi0
   ! The name of the KS independent-particle polarizability to be read (see driver.F90) 
   ! if no dataset mode             : abi//'SUS'
   ! if dataset mode, and getsuscep==0 : abi//'_DS'//trim(jdtset)//'SUS'
   ! if dataset mode, and getsuscep/=0 : abo//'_DS'//trim(jgetsuscep)//'SUS'

  character(len=fnlen) :: fildensin
   ! if no dataset mode             : abi//'DEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character(len=fnlen) :: filkdensin
   ! if no dataset mode             : abi//'KDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'KDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'KDEN'

  character(len=fnlen) :: filvhain
   ! if no dataset mode             : abi//'VHA'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'VHA'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'VHA'

  character(len=fnlen) :: filkss
   ! The name of the Kohn-Sham structure file to be read (see driver.F90) 
   ! if no dataset mode             : abi//'KSS'
   ! if dataset mode, and getkss==0 : abi//'_DS'//trim(jdtset)//'KSS'
   ! if dataset mode, and getkss/=0 : abo//'_DS'//trim(jgetkss)//'KSS'

  character(len=fnlen) :: filqps
   ! The name of the Quasi-Particle structure file to be read (see driver.F90) 
   ! if no dataset mode             : abi//'QPS'
   ! if dataset mode, and getqps==0 : abi//'_DS'//trim(jdtset)//'QPS'
   ! if dataset mode, and getqps/=0 : abo//'_DS'//trim(jgetqps)//'QPS'

  character(len=fnlen) :: filscr
   ! The name of the SCReening file (symmetrized inverse dielectric matrix) to be read (see driver.F90) 
   ! if no dataset mode             : abi//'SCR'
   ! if dataset mode, and getscr==0 : abi//'_DS'//trim(jdtset)//'SCR'
   ! if dataset mode, and getscr/=0 : abo//'_DS'//trim(jgetscr)//'SCR'

! character(len=fnlen) :: filpsp(ntypat)
   ! the filenames of the pseudopotential files, from the standard input.

  character(len=fnlen) :: filstat
   ! tmp//'_STATUS'

  character(len=fnlen) :: fnamewffk
   ! the name of the ground-state wavefunction file to be read (see driver.F90)

  character(len=fnlen) :: fnamewffq
   ! the name of the k+q ground-state wavefunction file to be read (see driver.F90)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewff1
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fildens1in   ! to be described by MVeithen

  character(len=fnlen) :: fname_tdwf   

  character(len=fnlen) :: fname_w90

  character(len=fnlen) :: fnameabi_hes
  character(len=fnlen) :: fnameabi_phfrq
  character(len=fnlen) :: fnameabi_phvec

  character(len=fnlen) :: fnametmp_wf1
  character(len=fnlen) :: fnametmp_wf2
  character(len=fnlen) :: fnametmp_1wf1
  character(len=fnlen) :: fnametmp_1wf2
  character(len=fnlen) :: fnametmp_wfgs
  character(len=fnlen) :: fnametmp_wfkq
   ! Set of filemanes formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _WF1, _WF2 ... 
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_kg
  character(len=fnlen) :: fnametmp_kgq
  character(len=fnlen) :: fnametmp_kg1
  character(len=fnlen) :: fnametmp_dum
  character(len=fnlen) :: fnametmp_ylm
  character(len=fnlen) :: fnametmp_ylm1
  character(len=fnlen) :: fnametmp_paw
  character(len=fnlen) :: fnametmp_paw1
  character(len=fnlen) :: fnametmp_pawq
   ! Set of filemanes formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _KG, _DUM, followed
   ! by the index of the processor. 
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_cg
  character(len=fnlen) :: fnametmp_cprj
  character(len=fnlen) :: fnametmp_eig
  character(len=fnlen) :: fnametmp_1wf1_eig
  character(len=fnlen) :: fnametmp_fft
  character(len=fnlen) :: fnametmp_kgs
  character(len=fnlen) :: fnametmp_sustr
  character(len=fnlen) :: fnametmp_tdexcit
  character(len=fnlen) :: fnametmp_tdwf

!@Bethe-Salpeter 
! New files introduced for the Bethe-Salpeter part.

!  character(len=fnlen) :: filresobsham
!   ! if no dataset mode             : abi//'RESONANT_BS_HAM'
!   ! if dataset mode, and get_reso_bsham==0 : abi//'_DS'//trim(jdtset)//'RESONANT_BS_HAM'
!   ! if dataset mode, and get_reso_bsham/=0 : abo//'_DS'//trim(jget_reso_bsham)//'RESONANT_BS_HAM'

!  character(len=fnlen) :: filcoupbsham
!   ! if no dataset mode             : abi//'COUPLING_BS_HAM'
!   ! if dataset mode, and get_coup_bsham==0 : abi//'_DS'//trim(jdtset)//'COUPLING_BS_HAM'
!   ! if dataset mode, and get_coup_bsham/=0 : abo//'_DS'//trim(jget_coup_bsham)//'COUPLING_BS_HAM'
!
  character(len=fnlen) :: filbseig
   ! The name of the file containing the eigenstates and eigenvalues of the Bethe-Salpeter Hamiltonian 
   ! or the iterative subspace generated by the Haydock method. (see driver.F90) 
   ! if no dataset mode             : abi//'BS_EIG'
   ! if dataset mode, and getbseig==0 : abi//'_DS'//trim(jdtset)//'BS_EIG'
   ! if dataset mode, and getbseig/=0 : abo//'_DS'//trim(jget_bseig)//'BS_EIG'

!END @BEthe-Salpeter

!The following filenames do not depend on itimimage, iimage and itime loops.

  character(len=fnlen) :: fnameabo_ae_wfk
  character(len=fnlen) :: fnameabo_ddb
  character(len=fnlen) :: fnameabo_den
  character(len=fnlen) :: fnameabo_dos
  character(len=fnlen) :: fnameabo_eelf
  character(len=fnlen) :: fnameabo_eig
  character(len=fnlen) :: fnameabo_eigi2d
  character(len=fnlen) :: fnameabo_eigr2d
  character(len=fnlen) :: fnameabo_em1
  character(len=fnlen) :: fnameabo_em1_lf
  character(len=fnlen) :: fnameabo_em1_nlf
  character(len=fnlen) :: fnameabo_exc_mdf
  character(len=fnlen) :: fnameabo_gkk
  character(len=fnlen) :: fnameabo_gw
  character(len=fnlen) :: fnameabo_gw_nlf_mdf
  character(len=fnlen) :: fnameabo_kss
  character(len=fnlen) :: fnameabo_moldyn
  character(len=fnlen) :: fnameabo_pot
  character(len=fnlen) :: fnameabo_qps
  character(len=fnlen) :: fnameabo_qps_ks
  character(len=fnlen) :: fnameabo_qp_den
  character(len=fnlen) :: fnameabo_qp_dos
  character(len=fnlen) :: fnameabo_qp_eig
  character(len=fnlen) :: fnameabo_rpa
  character(len=fnlen) :: fnameabo_rpa_nlf_mdf
  character(len=fnlen) :: fnameabo_scr
  character(len=fnlen) :: fnameabo_sgm
  character(len=fnlen) :: fnameabo_sgr
  character(len=fnlen) :: fnameabo_sig
  character(len=fnlen) :: fnameabo_spcur
  character(len=fnlen) :: fnameabo_sus
  character(len=fnlen) :: fnameabo_vso
  character(len=fnlen) :: fnameabo_wan
  character(len=fnlen) :: fnameabo_wfk
  character(len=fnlen) :: fnameabo_wfq
  character(len=fnlen) :: fnameabo_w90
  character(len=fnlen) :: fnameabo_1wf

!The following filenames are initialized only iniside itimimage, iimage and itime loops,
!and are appended with the adequate specifier 'app'.

  character(len=fnlen) :: fnameabo_app
  character(len=fnlen) :: fnameabo_app_atmden_core
  character(len=fnlen) :: fnameabo_app_atmden_full
  character(len=fnlen) :: fnameabo_app_atmden_val
  character(len=fnlen) :: fnameabo_app_bxsf
  character(len=fnlen) :: fnameabo_app_cml_xml
  character(len=fnlen) :: fnameabo_app_den
  character(len=fnlen) :: fnameabo_app_dos
  character(len=fnlen) :: fnameabo_app_elf
  character(len=fnlen) :: fnameabo_app_elf_down
  character(len=fnlen) :: fnameabo_app_elf_up
  character(len=fnlen) :: fnameabo_app_eig
  character(len=fnlen) :: fnameabo_app_fatbands
  character(len=fnlen) :: fnameabo_app_gden1
  character(len=fnlen) :: fnameabo_app_gden2
  character(len=fnlen) :: fnameabo_app_gden3
  character(len=fnlen) :: fnameabo_app_geo
  character(len=fnlen) :: fnameabo_app_kden
  character(len=fnlen) :: fnameabo_app_lden
  character(len=fnlen) :: fnameabo_app_pawden
  character(len=fnlen) :: fnameabo_app_pot
  character(len=fnlen) :: fnameabo_app_opt
  character(len=fnlen) :: fnameabo_app_opt2
  character(len=fnlen) :: fnameabo_app_stm
  character(len=fnlen) :: fnameabo_app_vha
  character(len=fnlen) :: fnameabo_app_vhxc
  character(len=fnlen) :: fnameabo_app_vxc
  character(len=fnlen) :: fnameabo_app_wfk
  character(len=fnlen) :: fnameabo_app_1dm

  character(len=fnlen) :: fnametmp_app_den
  character(len=fnlen) :: fnametmp_app_kden

 end type datafiles_type

end module defs_abitypes
!!***

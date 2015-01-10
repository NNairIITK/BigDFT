!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_datatypes
!! NAME
!! defs_datatypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! ABINIT package.
!!
!! List of datatypes :
!! * pseudopotential_gth_type : part of pseudopotential_type
!! * pseudopotential_type : datas for norm-conserving pseudopotential
!! * pspheader_paw_type : part of pspheader_type
!! * pspheader_type : data for the header of files
!! * dataset_type : the "dataset" for the main abinit code
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

module defs_datatypes

 use defs_basis

 implicit none
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pseudopotential_gth_type
!! NAME
!! pseudopotential_gth_type
!!
!! FUNCTION
!! This structure is a sub-structure of pseudopotential_type used to
!! store parameters from the GTH pseudo-potentials. All arrays have
!! indices running on 1:npsp for each read pseudo-file. The 'set' array
!! is a check array, since several different pseudo can be used in a simulation
!! it set a flag for each npsp if params have been set or not. This is
!! redundant with psps%pspcod in the way that when psps%pspcod(i) is 2,
!! then gth_params%set(i) is .true.. GTH pseudo previous to wavelets introduction
!! doesn't have geometric informations. These have been added on the last line.
!! It is three radius informations, the %hasGeometry flag is there to know
!! which kind of pseudo has been read.
!!
!! SOURCE

 type pseudopotential_gth_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp), pointer :: psppar(:, :, :)
   ! These are {rloc, C(1...4)} coefficients for psppar(0, :, :) indices,
   ! Followed by the h coefficients for psppar(1:2, :, :) indices.
   !  size (0:2, 0:4, npsp)

  real(dp), pointer :: radii_cov(:)
   ! The covalence radii for each pseudo (?) size (npsp)

  real(dp), pointer :: radii_cf(:, :)
   ! Cut-off radii for core part and long-range part.
   ! radii_cf(:, 1) is for the long-range cut-off and
   ! radii_cf(:, 2) is for the core cut-off. size (npsp, 2)

  integer, pointer :: semicore(:)
   ! The semicore code, indicated as an integer.
   ! The integer is the n_s + 4*n_p + 16* n_d + 64* n_f
   ! where n_l are the number of semicore orbitals for a given angular momentum
   ! starting from the lower level of course

  real(dp), pointer :: psp_k_par(:, :, :)
   ! Spin orbit coefficients in HGH/GTH formats: k11p etc... see psp3ini.F90
   !   dimension = num l channels, 3 coeffs, num psp = (1:lmax+1,1:3,npsp)

  logical, pointer :: hasGeometry(:)
   ! Flag for geometric informations in the pseudo. size (npsp)

  logical, pointer :: set(:)
   ! Consistency array, used for checking size (npsp)

 end type pseudopotential_gth_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pseudopotential_type
!! NAME
!! pseudopotential_type
!!
!! FUNCTION
!! This structured datatype contains all the information about one
!! norm-conserving pseudopotential, including the description of the local
!! and non-local parts, the different projectors, the non-linear core
!! correction ...
!!
!! SOURCE

 type pseudopotential_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Integer scalars

  integer :: dimekb
   ! Dimension of Ekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb=lnmax (lnmax: see this file)
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom type
   !                     dimekb=lmnmax*(lmnmax+1)/2 (lmnmax: see this file)

  integer :: lmnmax
   !  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
   !  If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lmnprojso or lnprojso, see pspheader_type

  integer :: lnmax
   !  Max. number of (l,n) components over all type of psps
   !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
   !  so, it is equal to the max of lnprojso, see pspheader_type

  integer :: mproj    ! TO BE SUPPRESSED
   ! Maximum number of non-local projectors over all angular momenta
   !  and type of psps
   ! 0 only if all psps are local

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"

  integer :: mpspso
   ! mpspso is set to 1 if none of the psps is used with a spin-orbit part (that
   !  is, if the user input variable so_psp is not equal
   !  to 1 in at least one case
   ! otherwise, it is set to 2

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1

  integer :: mqgrid_ff
   ! Number of points in the reciprocal space grid on which
   ! the radial functions ffspl are specified

  integer :: mqgrid_vl
   ! Number of points in the reciprocal space grid on which
   ! the radial functions vlspl are specified

  integer :: mtypalch
   ! Maximum number of alchemical pseudo atoms. If non-zero,
   ! the mechanism to generate mixing of pseudopotentials is activated

  integer :: npsp
   ! Number of types of pseudopotentials

  integer :: npspalch
   ! Number of types of pseudopotentials use for alchemical purposes

  integer :: ntypat
   ! Number of types of atoms (might be alchemy wrt pseudopotentials)

  integer :: ntypalch
   ! Number of types of alchemical pseudoatoms

  integer :: ntyppure
   ! Number of types of pure pseudoatoms

  integer :: n1xccc
   ! Number of radial points for the description of the pseudo-core charge
   ! (in the framework of the non-linear XC core correction)

  integer :: optnlxccc
   ! Option for the choice of non-linear XC core correction treatment (see the input variable)

  integer :: positron
   ! Option for the choice of type of GS calculation (electron or positron)

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Logical scalars

  logical :: vlspl_recipSpace
   ! governs if vlspl is compute in reciprocal space or in real
   ! space (when available).

! Integer arrays

  integer, pointer :: algalch(:)   
   ! algalch(ntypalch)
   ! For each type of pseudo atom, the algorithm to mix the pseudopotentials

  integer, pointer :: indlmn(:,:,:)  
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

  integer, pointer :: pspdat(:)  
   ! pspdat(ntypat)
   ! For each type of psp, the date of psp generation, as given by the psp file

  integer, pointer :: pspcod(:)  
   ! pspcod(npsp)
   ! For each type of psp, the format -or code- of psp generation,
   !  as given by the psp file

  integer, pointer :: pspso(:)  
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used

  integer, pointer :: pspxc(:)   
   ! pspxc(ntypat)
   ! For each type of psp, the XC functional that was used to generate it,
   ! as given by the psp file

! Real (real(dp)) arrays

  real(dp), pointer :: ekb(:,:)  
   ! ekb(dimekb,ntypat*(1-usepaw))
   !  ->NORM-CONSERVING PSPS ONLY:
   !    (Real) Kleinman-Bylander energies (hartree)
   !           for number of basis functions (l,n) (lnmax)
   !           and number of atom types (ntypat)
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb=lnmax*(lnmax+1)/2
   !             in the place of dimekb=lmnmax.

  real(dp), pointer :: ffspl(:,:,:,:)  
   ! ffspl(mqgrid_ff,2,lnmax,ntypat)
   ! Gives, on the radial grid, the different non-local projectors,
   ! in both the norm-conserving case, and the PAW case

  real(dp), pointer :: mixalch(:,:)  
   ! mixalch(npspalch,ntypalch)
   ! Mixing coefficients to generate alchemical pseudo atoms

  real(dp), pointer :: qgrid_ff(:)  
   ! qgrid_ff(mqgrid_ff)
   ! The coordinates of all the points of the radial grid for the nl form factors

  real(dp), pointer :: qgrid_vl(:)  
   ! qgrid_vl(mqgrid_vl)
   ! The coordinates of all the points of the radial grid for the local part of psp

  real(dp), pointer :: vlspl(:,:,:)  
   ! vlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the local part of each type of psp.

  real(dp), pointer :: dvlspl(:,:,:)  
   ! dvlspl(mqgrid_vl,2,ntypat)
   ! Gives, on the radial grid, the first derivative of the local
   ! part of each type of psp (computed when the flag 'vlspl_recipSpace' is true).

  real(dp), pointer :: xcccrc(:)  
   ! xcccrc(ntypat)
   ! Gives the maximum radius of the pseudo-core charge, for each type of psp.

  real(dp), pointer :: xccc1d(:,:,:)  
   ! xccc1d(n1xccc*(1-usepaw),6,ntypat)
   ! Norm-conserving psps only
   ! The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
   ! for each type of atom, on the radial grid. The components
   ! xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
   ! pseudo-core charge with respect to the radial distance.

  real(dp), pointer :: zionpsp(:)  
   ! zionpsp(npsp)
   ! For each pseudopotential, the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), pointer :: ziontypat(:)  
   ! ziontypat(ntypat)
   !  For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp), pointer :: znuclpsp(:)  
   ! znuclpsp(npsp)
   ! The atomic number of each pseudopotential

  real(dp), pointer :: znucltypat(:)  
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

! Character arrays

  character(len=fnlen), pointer :: filpsp(:)  
   ! filpsp(ntypat)
   ! The filename of the pseudopotential

  character(len=fnlen), pointer :: title(:)   
   ! title(ntypat)
   ! The content of first line read from the psp file

  type(pseudopotential_gth_type) :: gth_params
   ! Types for pseudo-potentials that are based on parameters. Currently, only
   ! GTH are supported (see pseudopotential_gth_type). To add one, one should
   ! create an initialisation method and a destruction method in 02psp (see
   ! psp2params.F90). These methods are called in driver().

 end type pseudopotential_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_paw_type
!! NAME
!! pspheader_paw_type
!!
!! FUNCTION
!! The pspheader_paw_type structured datatype gather additional information
!! about a PAW pseudopotential file, from its header.
!!
!! SOURCE

 type pspheader_paw_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: basis_size    ! Number of elements of the wf basis ((l,n) quantum numbers)
  integer :: l_size        ! Maximum value of l+1 leading to a non zero Gaunt coefficient
  integer :: lmn_size      ! Number of elements of the paw basis
  integer :: mesh_size     ! Dimension of (main) radial mesh
  integer :: pawver        ! Version number of paw psp format
  integer :: shape_type    ! Type of shape function
  real(dp) :: rpaw         ! Radius for paw spheres
  real(dp) :: rshp         ! Cut-off radius of shape function

 end type pspheader_paw_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pspheader_type
!! NAME
!! pspheader_type
!!
!! FUNCTION
!! The pspheader_type structured datatype gather different information
!! about a pseudopotential file, from its header.
!!
!! SOURCE

 type pspheader_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: nproj(0:3) ! number of scalar projectors for each angular momentum
  integer :: nprojso(3) ! number of spin-orbit projectors for each angular momentum
  integer :: lmax       ! maximum l quantum number (-1 if only local)
                        ! Example : s only       -> lmax=0
                        !           s and p      -> lmax=1
                        !           d only       -> lmax=2
  integer :: pspcod     ! code number of the pseudopotential
  integer :: pspdat     ! date of generation of the pseudopotential
  integer :: pspxc      ! exchange-correlation functional
  integer :: pspso      ! spin-orbit characteristics
  integer :: xccc       ! =0 if no XC core correction, non-zero if XC core correction

  real(dp) :: zionpsp       ! charge of the ion made of core electrons only
  real(dp) :: znuclpsp      ! atomic number of the nuclei

  real(dp) :: GTHradii(0:4) ! Radii values for GTH (and HGH) family potentials

  character(len=fnlen) :: filpsp   ! name of the psp file
  character(len=fnlen) :: title    ! content of first line read from the psp file

  type(pspheader_paw_type) :: pawheader ! only for PAW psps. See above

 end type pspheader_type
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

 end type dataset_type
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

end module defs_datatypes
!!***

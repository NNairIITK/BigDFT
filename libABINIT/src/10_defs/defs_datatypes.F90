!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_datatypes
!! NAME
!! defs_datatypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! ABINIT package.
!! If you want to add one new datatype, please, examine first whether
!! another datatype might meet your need (e.g. adding some records to it).
!! Then, if you are sure your new structured datatype is needed,
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
!! * bandstructure_type : different information about the band structure
!! * bcp_type : a "bonding critical point" for aim
!! * coeff?_type : a small datatype for ?D-arrays with dimensions depending on the type of atom
!! * datafil_type : the data (units,filenames) related to files
!! * efield_type : First-principles calculations in a finite electric field
!! * energies_type : simple datastructure to store parts of total energy.
!! * epsilonm1_results : for GW part of ABINIT, results of screening
!! * gs_hamiltonian_type : datastructure describing an Hamiltonian
!! * hdr_type   : the header of wf, den and pot files
!! * nuclear_type : data (esp. related to different nspden) at each nuclear site
!! * pawang_type : for PAW, ANGular mesh discretization and related data
!! * pawfgr_type : for PAW, Fine rectangular GRid parameters and related data
!! * pawfgrtab_type : for PAW, various arrays giving data related to fine grid for a given atom
!! * pawrad_type : for PAW, RADial mesh discretization and related data
!! * pawtab_type : for PAW, TABulated data initialized at start
!! * paw_an_type : for PAW, various arrays given on ANgular mesh or ANgular moments
!! * paw_ij_type : for PAW, various arrays given on (i,j) (partial waves) channels
!! * pawrhoij_type : for PAW, rhoij quantities and related data
!! * pseudopotential_type : for norm-conserving pseudopotential, all the
!!   information
!! * pspheader_paw_type : for PAW, the header of the atomic file
!! * pspheader_type : for norm-conserving pseudopotentials, the header of
!!   the file
!! * rdm_parameters : contains the parameters used during a RDM calculation
!! * results_gs_type : contains the results of a GS calculation
!! * results_out_type : contains a subset of the results, for internal
!!   tests and timing analysis
!! * scf_history : contains an history of previous SCF cycles (densities...)
!! * sigma_results : for GW part of ABINIT, results of sigma
!! * wvl_internalVars_type : all internal input variables used by wavelets.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! PAW structured datatypes to be described ...
!! * pawang_type : ANGular mesh discretization and related data
!! * pawfgr_type : Fine rectangular GRid parameters and related data
!! * pawfgrtab_type : various arrays giving data related to fine grid for a given atom
!! * pawrad_type :  RADial mesh discretization and related data
!! * pawtab_type : TABulated data initialized at start
!! * paw_an_type : various arrays given on ANgular mesh or
!! * paw_ij_type : various arrays given on (i,j) (partial waves) channels
!! * pawrhoij_type : for PAW, rhoij quantities and related data
!! * pspheader_paw_type: the header of the atomic file
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_datatypes

 use defs_basis
 use defs_parameters

#if defined HAVE_BIGDFT
 use BigDFT_API, only : atoms_data
#endif

 implicit none

!Structures

!!***


!!****t* defs_datatypes/wvl_internalVars_type
!! NAME
!! wvl_internalVars_type
!!
!! FUNCTION
!! This type is a gathering for all internal variables wavelets required. It is
!! included in the datatypes strutcture.
!!
!! NOTES
!! This array should be defined early since it is included in datatype.
!!
!! SOURCE

type wvl_internalVars_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: n(3)
  ! Size of the coarse grid, called n1, n2 and n3 in BigDFT.

  integer :: ni(3)
  ! Size of the fine grid used to expand potentials, densities...
  ! Each dimensions are equal to 2 * n(i) + buffer. The buffer value
  ! depends on the boundary conditions.
  ! They are called n1i, n2i and n3i in BigDFT.

  integer :: ntot
  ! Number of elements in the density/potential grid,
  ! equals to ni(1) * ni(2) * ni(3)

  integer :: fGrid(2, 3)
  ! Localisation information of a cube containing the fine grid.
  ! (1,:) gives the lower point and (2,:) the higher one.

  real(dp) :: h(3)
  ! The hgrid values in each direction, after the application of the
  ! boundary conditions. In free boundary conditions, the three values are equal.

#if defined HAVE_BIGDFT
  type(atoms_data) :: atoms
  ! A copy of the current dtset values.
#endif

end type wvl_internalVars_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/bandstructure_type
!! NAME
!! bandstructure_type
!!
!! FUNCTION
!! It contains different information about the band structure: eigenenergies, residuals, derivative of occupation
!! numbers vs energy in case of metallic occupations and Brillouin zone according to the context: k points,
!! occupation numbers, storage mode of wavefunctions, weights ...
!! For example, the initial Brillouin zone, set up in the dataset, will be treated in the response function part of
!! the code, to give a reduced Brillouin zone different from the original one, due to the breaking of the symmetries
!! related to the existence of a wavevector, or the lack of time-reversal invariance.
!!
!! SOURCE

 type bandstructure_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: bantot                ! total number of bands (sum(nband(:))
  integer :: mband                 ! Max number of bands i.e MAXVAL(nband) (to dimension arrays)
  integer :: nkpt                  ! number of k points
  integer :: nspinor               ! 1 for collinear, 2 for noncollinear.
  integer :: nsppol                ! number of spin-polarizations
  integer :: occopt                ! Occupation option, see input variable.

  !!$integer :: kptopt
  !!$real(dp) :: tolwfr
  !!$real(dp),pointer :: resid(mband*nkpt*nsppol)
  !!$resid(mband*nkpt*nsppol)=residuals (hartree**2)

  real(dp) :: entropy              ! Entropy associated with the smearing (adimensional)
  real(dp) :: fermie               ! Fermi energy
  real(dp) :: nelect               ! Number of electrons.
  real(dp) :: tphysel              ! Physical temperature of electrons.
  real(dp) :: tsmear               ! Temperature of smearing.

  integer,pointer :: istwfk(:)   
  ! istwfk(nkpt)
  ! Storage mode at each k point.

  integer,pointer :: nband(:)  
  ! nband(nkpt*nsppol)
  ! Number of bands at each k point and spin-polarisation.

  integer,pointer :: npwarr(:)   
  ! npwarr(nkpt)
  ! Number of plane waves at each k point.

  real(dp),pointer :: kptns(:,:)  
  ! kptns(3,nkpt)
  ! k-point vectors.

  real(dp),pointer :: eig(:,:,:)   
  ! eig(mband,nkpt,nsppol)
  ! Eigenvalues of each band.

  real(dp),pointer :: occ(:,:,:)   
  ! occ(mband,nkpt,nsppol)
  ! occupation of each band.

  real(dp),pointer :: doccde(:,:,:)   
  ! doccde(mband,nkpt,nsppol)
  ! derivative of the occupation of each band wrt energy (needed for RF).

  real(dp),pointer :: wtk(:)  
  ! wtk(nkpt)
  ! weight of each k point, normalized to one.

 end type bandstructure_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/bcp_type
!! NAME
!! bcp_type
!!
!! FUNCTION
!! a "bonding critical point" for aim
!!
!! SOURCE

 type bcp_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Integer
  integer :: iat     !! number of the bonding atom inside a primitive cell
  integer :: ipos    !! number of the primitive cell of the bonding atom

! Real
  real(dp) :: chg     !! charge at the critical point
  real(dp) :: diff(3) !! three distances : AT-CP,BAT-CP,AT-BAT
  real(dp) :: ev(3)   !! eigenvalues of the Hessian
  real(dp) :: pom(3)  !! position of the bonding atom
  real(dp) :: rr(3)   !! position of the bcp
  real(dp) :: vec(3,3)!! eigenvectors of the Hessian
  real(dp) :: vv(3)   !! position of the bcp relative to the central atom

 end type bcp_type
!!***

!!****t* defs_datatypes/macro_uj_type
!! NAME
!! dtmacro_uj
!!
!! FUNCTION
!! This data type contains the potential shifts and the occupations 
!! for the determination of U and J for the DFT+U calculations. 
!! iuj=1,2: non-selfconsistent calculations. iuj=3,4 selfconsistent calculations. 
!! iuj=2,4  => pawujsh<0 ; iuj=1,3 => pawujsh >0, 
!!
!! SOURCE

 type macro_uj_type

! Integer
  integer :: iuj        !! dataset treated
  integer :: nat        !! number of atoms U (J) is determined on
  integer :: ndtset     !! total number of datasets
  integer :: nspden     !! number of densities treated
  integer :: macro_uj   !! which mode the determination runs in
  integer :: pawujat    !! which atom U (J) is determined on
  integer :: pawprtvol  !! controlling amount of output
  integer :: option     !! controls the determination of U (1 with compensating charge bath, 2 without) 

! Real
  real(dp) :: diemix    !! mixing parameter
  real(dp) :: mdist     !! maximal distance of ions
  real(dp) :: pawujga   !! gamma for inversion of singular matrices

! Integer arrays
  integer , pointer  :: scdim(:)   !! size of supercell
 
! Real arrays
  real(dp) , pointer :: occ(:,:)    !! occupancies after a potential shift: occ(ispden,nat)
  real(dp) , pointer :: rprimd(:,:) !! unit cell for symmetrization
  real(dp) , pointer :: vsh(:,:)    !! potential shifts on atoms, dimensions: nspden,nat
  real(dp) , pointer :: xred(:,:)   !! atomic position for symmetrization 

 end type macro_uj_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/efield_type
!! NAME
!! efield_type
!!
!! FUNCTION
!! First-principles calculations in a finite electric field
!!
!! SOURCE

 type efield_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


! Integer variables
  integer :: has_qijb            ! 2 if qijb computed, 1 if only allocated, zero else
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: maxnstr             ! max number of strings along idir=1,2,3
  integer :: maxnkstr            ! max number of k-points per string
  integer :: mkmem_max           ! max of mkmem
  integer :: nband_occ           ! number of occupied bands
                                 ! this number must be the same for every k

! Integer arrays
  integer :: nstr(3)             ! nstr(idir) = number of strings along idir
  integer :: nkstr(3)            ! nkstr(idir) = number of k-points per string

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1
                                 !                         1 if nsppol = 2

! Real(dp) arrays
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
  real(dp) :: efield_dot(3)      ! reciprocal lattice coordinates of the
                                 ! electric field

! Integer pointers
  integer, pointer :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, pointer :: cgqindex(:,:,:) ! cgqindex(3,6,nkpt*nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cgq and pwnsfacq
                                      ! arrays
                                      ! (see vtorho.f and initberry.f)
  integer, pointer :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, pointer :: idxkstr(:,:,:)  ! idxkstr(maxnkstr,maxnstr,3)
                                      ! idxkstr(ikstr,istr,idir) index (ikpt) of
                                      ! k-point ikstr on string istr along idir
  integer, pointer :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtefield%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, pointer :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list
  integer, pointer :: nneigh(:)          ! nneigh(nkpt)
                                         ! for each k-point, nneigh stores
                                         ! the number of its nearest neighbours
                                         ! that are not related by symmetry
  integer, pointer :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg
  integer, pointer :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, pointer :: sflag(:,:,:,:)  ! sflag(nband_occ,nkpt*nsppol,2,3)
                                      ! sflag = 0 : compute the whole row of
                                      !             smat
                                      ! sflag = 1 : the row is up to date

! Real(dp) pointers
  real(dp), pointer :: fkptns(:,:)       ! fkptns(3,1:dtefield%fnkpt)
                                         ! k-points in FBZ

  real(dp), pointer :: smat(:,:,:,:,:,:)
! smat(2,nband_occ,nband_occ,nkpt*nsppol,2,3)
! Overlap matrix for every k-point. In an electric field calculation,
! smat is updated at every iteration.

! pointer to qijb
 type(qijb_type),pointer :: paw_qijb(:)


 end type efield_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/energies_type
!! NAME
!! energies_type
!!
!! FUNCTION
!! Simple datastructure to gather all part of total energy. Not all
!! attributes may have a value, depending on the scheme used to
!! compute the total energy and several options read from dtset.
!!
!! SOURCE

 type energies_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp) :: e_localpsp
   ! Local psp energy (hartree)

  real(dp) :: e_eigenvalues
   ! Sum of the eigenvalues - Band energy (Hartree)
   ! (valid for double-counting scheme dtset%optene == 1)

  real(dp) :: e_ewald
   ! Ewald energy (hartree), store also the ion/ion energy for free boundary conditions.

  real(dp) :: e_hartree
   ! Hartree part of total energy (hartree units)
  
  real(dp) :: e_corepsp
   ! psp core-core energy

  real(dp) :: e_kinetic
   ! Kinetic energy part of total energy.
   ! (valid for direct scheme, dtset%optene == 0)

  real(dp) :: e_nonlocalpsp
   ! Nonlocal pseudopotential part of total energy.

  real(dp) :: e_entropy
   ! Entropy energy due to the occupation number smearing (if metal)
   ! Value is multiplied by dtset%tsmear, see %entropy for the entropy alone.
   ! (valid for metals, dtset%occopt>=3 .and. dtset%occopt<=7)

  real(dp) :: entropy

  real(dp) :: e_xc
   ! Exchange-correlation energy (hartree)

  real(dp) :: e_vxc
   ! Potential exchange-correlation energy (hartree)

  real(dp) :: e_xcdc
   ! enxcdc=exchange-correlation double-counting energy (hartree)

  real(dp) :: e_paw
   ! PAW spherical part energy

  real(dp) :: e_pawdc
   ! PAW spherical part double-counting energy

  real(dp) :: e_elecfield
   ! Electric enthalpy, by adding both ionic and electronic contributions

  real(dp) :: e_fermie
   !  Fermie energy

  real(dp) :: h0
   !  h0=e_kinetic+e_localpsp+e_nonlocalpsp

  real(dp) :: e_electronpositron
   ! Electron-positron: electron-positron interaction energy

  real(dp) :: edc_electronpositron
   ! Electron-positron: double-counting electron-positron interaction energy

  real(dp) :: e0_electronpositron
   !  Electron-positron: energy only due to unchanged particles
   !                     (if calctype=1, energy due to electrons only)
   !                     (if calctype=2, energy due to positron only)

 end type energies_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/gs_hamiltonian_type
!! NAME
!! gs_hamiltonian_type
!!
!! FUNCTION
!! This datastructure contains the information about one Hamiltonian,
!! needed in the "getghc" routine, that apply the Hamiltonian
!! on a wavefunction.
!! About the non-local part of the Hamiltonian
!! The operator Onl has the following general form:
!! $Onl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!! Operator Onl is -- in the typical case -- the nonlocal potential.
!! - In a classical plane-wave calculation, $Enl^{R}_{lmn,l''m''n''}$ is the
!!   Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!! - In a PAW calculation, $Enl^{R}_{lmn,l''m''n''}$ can either be the nonlocal
!!   contribution to the Hamiltonian $D_{ij}$ or the overlap matrix $S_{ij}$.
!! - The |P_{Rlmn}> are the projector functions.
!!
!! SOURCE

 type gs_hamiltonian_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: dimekb1
   ! First dimension of Ekb (see ekb in this file)
   ! Same as psps%dimekb
   ! ->Norm conserving : Max. number of Kleinman-Bylander energies
   !                     for each atom type
   !                     dimekb1=lnmax
   ! ->PAW : Max. number of Dij coefficients connecting projectors
   !                     for each atom
   !                     dimekb1=lmnmax*(lmnmax+1)/2

  integer :: dimekb2
   ! Second dimension of Ekb (see ekb in this file)
   ! ->Norm conserving psps: dimekb2=ntypat
   ! ->PAW                 : dimekb2=natom

  integer :: istwf_k
   ! option parameter that describes the storage of wfs (k-dependent)

  integer :: lmnmax
   ! Maximum number of different l,m,n components over all types of psps.
   ! same as dtset%lmnmax

  integer :: matblk
   ! dimension of the array ph3d

  integer :: mgfft
   ! maximum size for 1D FFTs (same as dtset%mgfft)

  integer :: mproj  ! TO BE SUPPRESSED LATER
   ! Maximum number of non-local projectors over all angular momenta
   !  and type of psps
   ! 0 only if all psps are local
   ! same as psps%mproj

  integer :: mpsang
   ! Highest angular momentum of non-local projectors over all type of psps.
   ! shifted by 1 : for all local psps, mpsang=0; for largest s, mpsang=1,
   ! for largest p, mpsang=2; for largest d, mpsang=3; for largest f, mpsang=4
   ! This gives also the number of non-local "channels"
   ! same as psps%mpsang

  integer :: mpssoang
   ! Maximum number of channels, including those for treating the spin-orbit coupling
   ! when mpspso=1, mpssoang=mpsang
   ! when mpspso=2, mpssoang=2*mpsang-1
   ! same as psps%mpssoang

  integer :: natom
   ! The number of atoms for this dataset
   ! same as dtset%natom

  integer :: nfft
   ! number of FFT grid points
   ! same as dtset%nfft

  integer :: npw
   ! number of plane waves (k-dependent)

  integer:: nspinor
   ! Number of spinorial components

  integer :: ntypat
   ! Number of types of pseudopotentials
   ! same as dtset%ntypat

  integer :: nvloc
   ! Number of components of vloc
   ! usually, nvloc=1, except in the non-collinear magnetism case, where nvloc=4

  integer :: n4,n5,n6
   ! same as ngfft(4:6)

  integer :: usepaw
   ! if usepaw=0 , use norm-conserving psps part of the code
   ! is usepaw=1 , use paw part of the code

  integer :: useylm
   ! governs the way the nonlocal operator is to be applied:
   !   1=using Ylm, 0=using Legendre polynomials

! Integer arrays

  integer, pointer :: atindx(:)    
   ! atindx(natom)
   ! index table for atoms (see scfcv.f)

  integer, pointer :: atindx1(:)   
   ! atindx1(natom)
   ! index table for atoms, inverse of atindx (see scfcv.f)

  integer, pointer :: gbound(:,:)  
   ! gbound(2*mgfft+8,2)
   ! G sphere boundary

  integer, pointer :: indlmn(:,:,:)  
   ! indlmn(6,lmnmax,ntypat)
   ! For each type of psp,
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

! integer, pointer :: indpw_k(:,:)
   ! indpw_k(4,npw)
   ! array which gives fft box index for given basis sphere
   ! This component was taken away : CPU time problem !

! integer, pointer :: kg_k(:,:)
   ! kg_k(3,npw)
   ! G vector coordinates with respect to reciprocal lattice translations
   ! This component was taken away : CPU time problem !

  integer, pointer :: nattyp(:)  
   ! nattyp(ntypat)
   ! # of atoms of each type

  integer :: ngfft(18)
   ! ngfft(1:3)=integer fft box dimensions
   ! ngfft(4:6)=integer fft box dimensions, might be augmented for CPU speed
   ! ngfft(7)=fftalg
   ! ngfft(8)=fftalg
   ! same as dtset%ngfft

  integer :: nloalg(5)
   ! governs the choice of the algorithm for non-local operator
   ! same as dtset%nloalg

  integer, pointer :: pspso(:)  
   ! pspso(ntypat)
   ! For each type of psp, 1 if no spin-orbit component is taken
   ! into account, 2 if a spin-orbit component is used

  integer, pointer :: typat(:)  
   ! typat(natom)
   ! type of each atom

! Real (real(dp)) scalar

  real(dp) :: ucvol
   ! unit cell volume (Bohr**3)

! Real (real(dp)) arrays

  real(dp), pointer :: ekb(:,:,:)   
   ! ekb(dimekb1,dimekb2,nspinor**2)
   !  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
   !          for number of basis functions (l,n) (lnmax)
   !          and number of atom types (ntypat)
   !          dimekb1=lnmax ; dimekb2=ntypat
   !  ->PAW : (Real, symmetric) Frozen part of Dij coefficients
   !                            to connect projectors
   !          for number of basis functions (l,m,n) (lmnmax)
   !          and number of atom (natom)
   !          dimekb1=lmnmax*(lmnmax+1)/2 ; dimekb2=natom
   ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
   !             lnmax); it would be easy to give it a second
   !             (symmetric) dimension by putting
   !             dimekb1=lnmax*(lnmax+1)/2
   !             in the place of dimekb1=lnmax.
   ! %ekb is spin dependent in the case of PAW calculations.

  real(dp), pointer :: sij(:,:)   
   ! sij(dimekb1,ntypat*usepaw) = overlap matrix for paw calculation

! real(dp), pointer :: ffnl(:,:,:,:)
   ! ffnl(npw,2,lmnmax,ntypat)
   ! nonlocal form factors
   ! This component was taken away : CPU time problem !

  real(dp) :: gmet(3,3)
   ! reciprocal space metric tensor in Bohr**-2

  real(dp) :: gprimd(3,3)
   ! dimensional reciprocal space primitive translations (Bohr^-1)

! real(dp), pointer :: kinpw(:)
   ! kinpw(npw)
   ! (modified) kinetic energy for each plane wave (Hartree)
   ! This component was taken away : CPU time problem !

  real(dp) :: kpoint(3)
   ! dimensionless k point coordinates wrt reciprocal lattice vectors. (k-dependent).

  real(dp), pointer :: phkxred(:,:)  
   ! phkxred(2,natom)
   ! phase factors exp(2 pi k.xred) (k-dependent)

  real(dp), pointer :: ph1d(:,:)   
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim phase arrays for structure factor (see getph.f).

! real(dp), pointer :: ph3d(:,:,:)
   ! ph3d(2,npw,matblk)
   ! 3-dim structure factors, for each atom and plane wave
   ! This component was taken away : CPU time problem !

! real(dp), pointer :: vlocal(:,:,:,:)
   ! vlocal(n4,n5,n6,nvloc)
   ! local potential in real space, on the augmented fft grid
   ! This component was taken away : CPU time problem !

  real(dp),pointer :: xred(:,:)  
   ! xred(3,natom)
   ! reduced coordinates of atoms (dimensionless)

! real(dp),pointer :: ylm(:,:)
   ! ylm(npw,mpsang*mpsang*useylm)
   ! Real spherical harmonics for each G
   ! This component was taken away : CPU time problem !

 end type gs_hamiltonian_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/hdr_type
!! NAME
!! hdr_type
!!
!! FUNCTION
!! It contains all the information needed to write a header for a
!! wf, den or pot file.
!! The structure of the header is explained in the abinis_help.html file.
!! The datatype is considered as an object, to which are attached a whole
!! set of "methods", actually, different subroutines.
!! A few of these subroutines are : hdr_init, hdr_update, hdr_clean,
!! hdr_check, hdr_io, hdr_skip.
!!
!! SOURCE

 type hdr_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
  integer :: date          ! starting date
  integer :: headform      ! format of the header
  integer :: intxc         ! input variable
  integer :: ixc           ! input variable
  integer :: natom         ! input variable
  integer :: nkpt          ! input variable
  integer :: npsp          ! input variable
  integer :: nspden        ! input variable
  integer :: nspinor       ! input variable
  integer :: nsppol        ! input variable
  integer :: nsym          ! input variable
  integer :: ntypat        ! input variable
  integer :: occopt        ! input variable
  integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
  integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
  integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)

! This record is not a part of the hdr_type, although it is present in the
! header of the files. This is because it depends on the kind of file
! that is written, while all other information does not depend on it.
! It was preferred to let it be initialized or defined outside of hdr_type.
! integer :: fform         ! file descriptor (or file format)

  integer :: ngfft(3)      ! input variable
  integer :: nwvlarr(2)    ! nwvlarr(2) array holding the number of wavelets for each resolution.

  ! MG: nullifying the pointers renders the object undefined, maybe intent(out) somewhere?
  integer, pointer :: istwfk(:)   ! 
  ! input variable istwfk(nkpt)

  integer, pointer :: lmn_size(:)   !
  ! lmn_size(npsp) from psps

  integer, pointer :: nband(:)      !
  ! input variable nband(nkpt*nsppol)

  integer, pointer :: npwarr(:)     !
  ! npwarr(nkpt) array holding npw for each k point

  integer, pointer :: pspcod(:)    !
  ! pscod(npsp) from psps

  integer, pointer :: pspdat(:)   !
  ! psdat(npsp) from psps

  integer, pointer :: pspso(:)    !
  ! pspso(npsp) from psps

  integer, pointer :: pspxc(:)    !
  ! pspxc(npsp) from psps

  integer, pointer :: so_psp(:)   !
  ! input variable so_psp(npsp)

  integer, pointer :: symafm(:)   !
  ! input variable symafm(nsym)

  integer, pointer :: symrel(:,:,:)    !
  ! input variable symrel(3,3,nsym)

  integer, pointer :: typat(:)    !
  ! input variable typat(natom)

  real(dp) :: ecut                  ! input variable
  real(dp) :: ecutdg                ! input variable (ecut for NC psps, pawecutdg for paw)
  real(dp) :: ecutsm                ! input variable
  real(dp) :: ecut_eff              ! ecut*dilatmx**2 (dilatmx is an input variable)
  real(dp) :: etot                  ! EVOLVING variable
  real(dp) :: fermie                ! EVOLVING variable
  real(dp) :: residm                ! EVOLVING variable
  real(dp) :: stmbias               ! input variable
  real(dp) :: tphysel               ! input variable
  real(dp) :: tsmear                ! input variable

  real(dp) :: qptn(3)               ! the wavevector, in case of a perturbation
  real(dp) :: rprimd(3,3)           ! EVOLVING variables

  real(dp), pointer :: kptns(:,:)   !
  ! input variable kptns(3,nkpt)

  real(dp), pointer :: occ(:)       !
  ! EVOLVING variable occ(bantot)

  real(dp), pointer :: tnons(:,:)   ! 
  ! input variable tnons(3,nsym)

  real(dp), pointer :: wtk(:)       !
  ! weight of kpoints wtk(nkpt)

  real(dp), pointer :: xred(:,:)    !
  ! EVOLVING variable xred(3,natom)

  real(dp), pointer :: zionpsp(:)   ! 
  ! zionpsp(npsp) from psps

  real(dp), pointer :: znuclpsp(:)  ! 
  ! znuclpsp(npsp) from psps
  ! Note the difference between znucl and znuclpsp !

  real(dp), pointer :: znucltypat(:)  !
  ! znucltypat(ntypat) from alchemy

  character(len=6) :: codvsn              
  ! version of the code

  character(len=132), pointer :: title(:)  !
  ! title(npsp) from psps

  type(pawrhoij_type), pointer :: pawrhoij(:)  !
  ! EVOLVING variable, only for paw

!Should make a list of supplementary infos
! MG: For postprocessing purposes, it is quite useful to
!  have kptrlatt as well as nshiftk and shiftk. also kptopt is useful
!  to know if time reversal can be employed

 end type hdr_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/mttk_type
!! NAME
!! mttk_type
!!
!! FUNCTION
!! For Martyna et al. (TTK) reversible MD integration scheme and related data
!!
!! SOURCE

 type mttk_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Real (double precision) scalars

   real(dp) :: glogv
    !Logarithm of the volume

   real(dp) :: vlogv
    !Derivative of logv

!Real (double precision) arrays

  real(dp) :: gboxg(3,3)
   !Imbalance in pressure (see paper)

  real(dp) :: vboxg(3,3)
   !Velocity of log rprimd (see paper)

  real(dp), pointer :: glogs(:)
   ! glogs(nnos)
   ! Imbalance of kinetic energy

  real(dp), pointer :: vlogs(:)
   ! vlogs(nnos)
   ! Velocities of thermostat variables

  real(dp), pointer :: xlogs(:)
   ! xlogs(nnos)
   ! Positions of thermostat variables

 end type mttk_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/nuclear_type
!! NAME
!! nuclear_type
!!
!! FUNCTION
!! Property results typically at each atom for each nspden. This appears to
!! be necessary because in PAW calculations there can be different nspden values
!! at each nuclear site.
!!
!! SOURCE

 type nuclear_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Real (real(dp)) arrays

  real(dp), pointer :: spden(:)
   ! spden(nspden)
   ! data for each ispden value; note that nspden for a given nucleus will
   ! typically be retrieved from pawrhoij(iatom)%nspden and hence is nucleus
   ! specific

 end type nuclear_type
!!***

!----------------------------------------------------------------------


!!****t* defs_datatypes/pawfgr_type
!! NAME
!! pawfgr_type
!!
!! FUNCTION
!! For PAW, Fine rectangular GRid parameters and related data
!!
!! SOURCE

 type pawfgr_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: mgfft, nfft
   ! Values of mffft and nfft for the fine rectangular grid:
   !   mgfft= max(ngfft(i)) [max. size of 1D FFT grid]
   !   nfft=ngfft1*ngfft2*ngfft3 [number of pts in the FFT box]

  integer :: mgfftc, nfftc
   ! Values of mffft and nfft for the COARSE rectangular grid:
   !   mgfftc= max(ngfftc(i)) [max. size of 1D FFT grid]
   !   nfftc=ngfftc1*ngfftc2*ngfftc3 [number of pts in the FFT box]

  integer :: usefinegrid
   ! Flag: =1 if a double-grid is used to convert spherical data
   !       to Fourier grid. =0 otherwise

!Integer arrays

  integer, pointer :: coatofin(:)
   ! coatofin(nfftc)
   ! Index of the points of the coarse grid on the fine grid

  integer, pointer :: fintocoa(:)
   ! fintocoa(nfft)
   ! Index of the points of the fine grid on the coarse grid
   !  (=0 if the point of the fine grid does not belong to the coarse grid)

  integer :: ngfft(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the fine rectangular grid

  integer :: ngfftc(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the COARSE rectangular grid

!Real (real(dp))

  real(dp) :: gsqcut
   ! Fourier cutoff on G^2 for "large sphere" of radius double
   ! that of the basis sphere corresponding to paw_ecutdg
   ! (concerns the fine rectangular grid)

 end type pawfgr_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pawfgrtab_type
!! NAME
!! pawfgrtab_type
!!
!! FUNCTION
!! For PAW, various arrays giving data related to fine grid for a given atom
!!
!! SOURCE

 type pawfgrtab_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: cplex
   ! cplex=1 if potentials/densities are real, 2 if they are complex

  integer :: l_size
   ! 1+maximum value of l leading to non zero Gaunt coeffs
   ! for the considered atom type

  integer :: gylm_allocated
   ! 1 if gylm() is allocated (and computed)

  integer :: gylmgr_allocated
   ! 1 if gylmgr() is allocated (and computed)

  integer :: gylmgr2_allocated
   ! 1 if gylmgr2() is allocated (and computed)

  integer :: nfgd
   ! Number of Fine rectangular GriD points
   ! in the paw sphere around considered atom

  integer :: nhatfr_allocated
   ! 1 if nhatfr() is allocated (and computed)

  integer :: nspden
   ! Number of spin-density components

  integer :: rfgd_allocated
   ! 1 if rfgd() is allocated (and computed)

  integer :: vlocgr_allocated
   ! 1 if vlocgr() is allocated (and computed)

!Integer arrays

  integer, pointer :: ifftsph(:)
   ! ifftsph(nfgd)
   ! Array giving the FFT index (fine grid) of a point in the paw
   ! sphere around considered atom (ifftsph=ix+n1*(iy-1+n2*(iz-1))

!Real (real(dp)) arrays

  real(dp), pointer :: gylm(:,:)
   ! gylm(nfgd,l_size*l_size)
   ! Gives g_l(r)*Y_lm(r) on the fine rectangular grid
   ! around considered atom

  real(dp), pointer :: gylmgr(:,:,:)
   ! gylmgr(3,nfgd,l_size*l_size)
   ! Gives the gradient of g_l(r)*Y_lm(r) wrt cart. coordinates
   ! on the fine rectangular grid around considered atom

  real(dp), pointer :: gylmgr2(:,:,:)
   ! gylmgr(6,nfgd,l_size*l_size)
   ! Gives the second gradient of g_l(r)*Y_lm(r) wrt cart. coordinates
   ! on the fine rectangular grid around considered atom

  real(dp), pointer :: nhatfr(:,:)
   ! nhatfr(cplex*nfgd,nspden)
   ! Gives the gradient of local potential wrt cart. coordinates
   ! on the fine rectangular grid around considered atom
   ! Only use in response function calculations

  real(dp), pointer :: rfgd(:,:)
   ! r(3,nfgd)
   ! Gives all R vectors (r-r_atom) on the Fine rectangular GriD
   ! around considered atom in Cartesian coordinates.

  real(dp), pointer :: vlocgr(:,:)
   ! vlocgr(3,nfgd)
   ! Gives the gradient of local potential wrt cart. coordinates
   ! on the fine rectangular grid around considered atom
   ! Only use in response function calculations

 end type pawfgrtab_type
!!***

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

!!****t* defs_datatypes/paw_an_type
!! NAME
!! paw_an_type
!!
!! FUNCTION
!! For PAW, various arrays given on ANgular mesh or ANgular moments
!!
!! SOURCE

 type paw_an_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: angl_size
   ! Dimension of paw angular mesh (angl_size=ntheta*nphi)

  integer :: cplex
   ! cplex=1 if potentials/densities are real, 2 if they are complex

  integer :: has_kxc
   ! set to 1 if xc kernels kxc1 and kxct1 are allocated and used
   !        2 if they are already computed

  integer :: has_vhartree
   ! set to 1 if vh1 and vht1 are allocated and used
   !        2 if they are already computed

  integer :: has_vxc
   ! set to 1 if vxc1 and vxct1 are allocated and used
   !        2 if they are already computed

  integer :: has_vxcval
   ! set to 1 if vxc1_val and vxct1_val are allocated and used
   !        2 if they are already computed

  integer :: lm_size
   ! lm_size=(l_size)**2
   ! l is Maximum value of l+1 leading to non zero Gaunt coeffs (l_size=2*l_max+1)

  integer :: mesh_size
   ! Dimension of radial mesh

  integer :: nkxc1
   ! number of independent components of Kxc1 and Kxct1
   ! (usually 3 for LDA, 23 for GGA)

  integer :: nspden
   ! Number of spin-density components

!Logical arrays

  logical, pointer :: lmselect(:)
   ! lmselect(lm_size)
   ! lmselect(ilm)=select the non-zero LM-moments of "one-center" densities/potentials

!Real (real(dp)) arrays

  real(dp), pointer :: kxc1 (:,:,:)
   ! kxc1(cplex*mesh_size,lm_size or angl_size,nkxc1)
   ! Gives xc kernel inside the sphere
   !   (theta,phi) values of kernel if pawxcdev=0
   !   LM-moments of kernel if pawxcdev/=0

  real(dp), pointer :: kxct1 (:,:,:)
   ! kxct1(cplex*mesh_size,lm_size or angl_size,nkxc1)
   ! Gives xc pseudo kernel inside the sphere
   !   (theta,phi) values of kernel if pawxcdev=0
   !   LM-moments of kernel if pawxcdev/=0

  real(dp), pointer :: vh1 (:,:,:)
   ! vh1(cplex*mesh_size,lm_size,nspden)
   ! Gives Hartree potential LM-moments inside the sphere

  real(dp), pointer :: vht1 (:,:,:)
   ! vht1(cplex*mesh_size,lm_size,nspden)
   ! Gives Hartree tilde potential LM-moments inside the sphere

  real(dp), pointer :: vxc1 (:,:,:)
   ! vxc1(cplex*mesh_size,lm_size or angl_size,nspden)
   ! Gives xc potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), pointer :: vxc1_val (:,:,:)
   ! vxc1_val(cplex*mesh_size,lm_size or angl_size,nspden) (Usually real, Mainly used for GW)
   ! Gives xc potential inside the sphere arising from valence only electrons
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), pointer :: vxct1 (:,:,:)
   ! vxct1(cplex*mesh_size,angl_size,nspden)
   ! Gives xc pseudo potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), pointer :: vxct1_val (:,:,:)
   ! vxct1_val(cplex*mesh_size,angl_size,nspden) (Usually real, Mainly used for GW)
   ! Gives xc pseudo potential inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

  real(dp), pointer :: vxc_ex (:,:,:)
   ! vxc_ex(cplex*mesh_size,angl_size,nspden)
   ! Gives xc  potential for local exact exchange inside the sphere
   !   (theta,phi) values of potential if pawxcdev=0
   !   LM-moments of potential if pawxcdev/=0

 end type paw_an_type
!!***

!-------------------------------------------------------------------------

!!****t* defs_datatypes/Paw_an_flags_type
!! NAME
!! Paw_an_flags_type
!!
!! FUNCTION
!! For PAW, the various "has flags" defined in the paw_an_type
!!
!! SOURCE

 type Paw_an_flags_type

  integer :: has_kxc
  integer :: has_vhartree
  integer :: has_vxc
  integer :: has_vxcval

 end type Paw_an_flags_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/paw_ij_type
!! NAME
!! paw_ij_type
!!
!! FUNCTION
!! For PAW, various arrays given on (i,j) (partial waves) channels
!!
!! SOURCE

 type paw_ij_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: cplex
   ! cplex=1 if all on-site PAW quantities are real, 2 if they are complex
   ! cplex=2 is useful for RF calculations

  integer :: cplex_dij
   ! cplex=1 if dij are real, 2 if they are complex

  !!$integer :: has_dijexxcore
   ! 1 if dijexxcore is allocated
   ! 2 if dijexxcore is already computed

  integer :: has_dijfr
   ! 1 if dijfr is allocated
   ! 2 if dijfr is already computed

  integer :: has_dijhartree
   ! 1 if dijhartree is allocated
   ! 2 if dijhartree is already computed

  integer :: has_dijhat
   ! 1 if dijhat is allocated
   ! 2 if dijhat is already computed

  integer :: has_dijso
   ! 1 if dijso is associated and used, 0 otherwise
   ! 2 if dijso is already computed

  integer :: has_dijU
   ! 1 if dijU is associated and used, 0 otherwise
   ! 2 if dijU is already computed

  integer :: has_dijxc
   ! 1 if dijxc is associated and used, 0 otherwise
   ! 2 if dijxc is already computed

  integer :: has_dijxc_val
   ! 1 if dijxc_val is associated and used, 0 otherwise
   ! 2 if dijxc_val is already computed

  integer :: has_dijR
   ! 1 if dijR is associated and used, 0 otherwise
   ! 2 if dijR is already computed

  integer :: has_dijL
   ! 1 if dijL is associated and used, 0 otherwise
   ! 2 if dijL is already computed

  integer :: has_dijLr3
   ! 1 if dijLr3 is associated and used, 0 otherwise
   ! 2 if dijLr3 is already computed

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: ndij
   ! Number of components of dij
   ! Usually ndij=nspden, except for nspinor==2 (where ndij=nspinor**2)

  integer :: nspden
   ! Number of spin-density components (may be different from dtset%nspden if spin-orbit)

  integer :: nsppol
   ! Number of independant spin-components

!Real (real(dp)) arrays

  real(dp), pointer :: dij(:,:)
   ! dij(cplex_dij*lmn2_size,ndij)
   ! Dij term (non-local operator)
   ! May be complex if cplex_dij=2
   !  dij(:,:,1) contains Dij^up-up
   !  dij(:,:,2) contains Dij^dn-dn
   !  dij(:,:,3) contains Dij^up-dn (only if nspinor=2)
   !  dij(:,:,4) contains Dij^dn-up (only if nspinor=2)

  !!$real(dp),pointer :: dijexxcore(:,:)
  ! dijexxcore(cplex_dij*lmn2_size,ndij)
  ! Onsite matrix elements of the Fock operator generated by core electrons

  real(dp), pointer :: dijfr(:,:)
   ! dijhat(cplex_dij*lmn2_size,ndij)
   ! For response function calculation only
   ! RF Frozen part of Dij (depends on q vector but not on 1st-order wave function)
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijhartree(:)
   ! dijhartree(cplex*lmn2_size)
   ! Dij_hartree term
   ! Contains all contributions to Dij from hartree
   ! Warning: Dimensioned by cplex, not cplex_dij
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijhat(:,:)
   ! dijhat(cplex_dij*lmn2_size,ndij)
   ! Dij_hat term (non-local operator) i.e \sum_LM \int_FFT Q_{ij}^{LM} vtrial
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijU(:,:)
   ! dijU(cplex_dij*lmn2_size,ndij)
   ! Onsite matrix elements of the U part of the PAW Hamiltonian.
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijso(:,:)
   ! dijso(cplex_dij*lmn2_size,ndij)
   ! Onsite matrix elements of L.S i.e <phi_i|L.S|phi_j>
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijxc(:,:)
   ! dijxc(cplex_dij*lmn2_size,ndij)
   ! Onsite matrix elements of vxc i.e
   ! <phi_i|vxc[n1+nc]|phi_j> - <tphi_i|vxc(tn1+nhat+tnc]|tphi_j>
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijxc_val(:,:)
   ! dijxc_val(cplex_dij*lmn2_size,ndij)
   ! Onsite matrix elements of valence-only vxc i.e
   ! <phi_i|vxc[n1]|phi_j> - <tphi_i|vxc(tn1+nhat]|tphi_j>
   ! Same storage as Dij (see above)

  real(dp), pointer :: dijR(:,:)
   ! dijR(3,lmn2_size)
   ! Onsite matrix elements of r-R, i.e., <phi_i|r-R|phi_j>
   ! needed for PAW treatment of electric field perturbation
   ! as r-R is a vector, there are three directions

  real(dp), pointer :: dijL(:,:)
   ! dijL(3,lmn2_size)
   ! Onsite matrix elements of L, the orbital angular momentum.
   ! needed for PAW treatment of magnetic field perturbation.
   ! As L is a vector, there are three directions, also these are
   ! always pure imaginary, only the imaginary part is saved here

  real(dp), pointer :: dijLr3(:,:)
   ! dijLr3(3,lmn2_size)
   ! Onsite matrix elements of L/r^3, the orbital angular momentum,
   ! weighted by 1/r^3. Needed for PAW treatment of magnetic field perturbation.
   ! As L is a vector, there are three directions, also these are
   ! always pure imaginary, only the imaginary part is saved here

  real(dp), pointer :: noccmmp(:,:,:,:)
   ! noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,nocc_nspden)
   ! cplex_dij=1 if collinear
   ! cplex_dij=2 if spin orbit is used
   ! cplex_dij=2 is used if non-collinear (for coherence, it is not necessary in this case, however)
   ! gives occupation matrix for lda+u (computed in setnoccmmp)
   ! Stored as: noccmmp(:,:,1)=   n^{up,up}_{m,mp}
   !            noccmmp(:,:,2)=   n^{dn,dn}_{m,mp}
   !            noccmmp(:,:,3)=   n^{up,dn}_{m,mp}
   !            noccmmp(:,:,4)=   n^{dn,up}_{m,mp}
   ! noccmmp(m,mp,:) is computed from rhoij(klmn) with  m=klmntomn(2)>mp=klmntomn(1)

  real(dp), pointer :: nocctot(:)
   ! nocctot(nspden)
   ! gives trace of occupation matrix for lda+u (computed in pawdenpot)
   ! for each value of ispden (1 or 2)

  real(dp), pointer :: vpawx(:,:,:)
   ! vpawx(2*lexexch+1,2*lexexch+1,nspden)
   ! exact exchange potential

 end type paw_ij_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/paw_ij_flags_type
!! NAME
!! paw_ij_flags_type
!!
!! FUNCTION
!! For PAW, the various "has flags" defined in the paw_ij_type
!!
!! SOURCE

 type paw_ij_flags_type

  integer :: has_dijfr
  integer :: has_dijhartree
  integer :: has_dijhat
  integer :: has_dijso
  integer :: has_dijU
  integer :: has_dijxc
  integer :: has_dijxc_val
  integer :: has_dijR
  integer :: has_dijL
  integer :: has_dijLr3

 end type paw_ij_flags_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/pawrhoij_type
!! NAME
!! pawrhoij_type
!!
!! FUNCTION
!! For PAW, rhoij quantities and related data
!!
!! SOURCE

 type pawrhoij_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


!Integer scalars

  integer :: cplex
   ! cplex=1 if rhoij are real, 2 if rhoij are complex

  integer :: itypat
   ! itypat=type of atom corresponding to this rhoij

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
   !           i.e. number of rhoij elements being mixed during SCF cycle
   ! lmnmix_sz=0 if mixing data are not used

  integer :: ngrhoij
   ! First dimension of array grhoij

  integer :: nrhoijsel
   ! nrhoijsel
   ! Number of non-zero value of rhoij
   ! This is the size of rhoijp(:,:) (see below in this datastructure)

  integer :: nspden
   ! Number of spin-density components for rhoij (may be different from nspden for density)

  integer :: nsppol
   ! Number of independant spin-components

  integer :: use_rhoij_
   ! 1 if pawrhoij%rhoij_ is allocated

  integer :: use_rhoijres
   ! 1 if pawrhoij%rhoijres is allocated

!Integer arrays

  integer, pointer :: kpawmix(:)
   ! kpawmix(lmnmix_sz)
   ! Indirect array selecting the elements of rhoij
   ! being mixed during SCF cycle

  integer, pointer :: rhoijselect(:)
   ! rhoijselect(lmn2_size)
   ! Indirect array selecting the non-zero elements of rhoij:
   ! rhoijselect(isel,ispden)=klmn if rhoij(klmn,ispden) is non-zero

!Real (real(dp)) arrays

  real(dp), pointer :: grhoij (:,:,:)
   ! grhoij(ngrhoij,cplex*lmn2_size,nspden)
   ! Gradients of Rho_ij wrt xred, strains, ... (non-packed storage)

  real(dp), pointer :: rhoij_ (:,:)
   ! rhoij_(cplex*lmn2_size,nspden)
   ! Array used to (temporary) store Rho_ij in a non-packed storage mode

  real(dp), pointer :: rhoijp (:,:)
   ! rhoijp(cplex*lmn2_size,nspden)
   ! Augmentation waves occupancies Rho_ij
   ! in PACKED STORAGE (only non-zero elements are stored)

  real(dp), pointer :: rhoijres (:,:)
   ! rhoijres(cplex*lmn2_size,nspden)
   ! Rho_ij residuals during SCF cycle (non-packed storage)

 end type pawrhoij_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/qijb_type
!! NAME
!! qijb_type
!!
!! FUNCTION
!!
!! SOURCE

 type qijb_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: lmn_size
  integer :: lmn2_size

!Real (real(dp)) arrays

  real(dp), pointer :: qijb(:,:,:)

 end type qijb_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/cprj_type
!! NAME
!! cprj_type
!!
!! FUNCTION
!! <p_lmn|Cnk> projected scalars and derivatives
!!             where |p_lmn> are non-local projectors for a given atom
!!                   |Cnk> is a wave function
!! Used only when useylm=1
!!
!! SOURCE

 type cprj_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.


!Integer scalars

  integer :: ncpgr
   ! Number of gradients of cp=<p_lmn|Cnk>

  integer :: nlmn
   ! Number of (l,m,n) non-local projectors

!Real (real(dp)) arrays

  real(dp), pointer :: cp (:,:)
   ! cp(2,nlmn)
   ! <p_lmn|Cnk> projected scalars for a given atom and wave function

  real(dp), pointer :: dcp (:,:,:)
   ! dcp(2,ncpgr,nlmn)
   ! derivatives of <p_lmn|Cnk> projected scalars for a given atom and wave function

 end type cprj_type
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

!!****t* defs_datatypes/rdm_parameters
!! NAME
!! rdm_parameters
!!
!! FUNCTION
!! For the RDM part of ABINIT, the rdm_parameters structured datatype
!! gather different parameters that characterize a RDM calculation.
!!
!!
!! SOURCE

 type rdm_parameters

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: npwvec                      ! Max between npwe and npwwfn, used to pass the dimension of arrays e.g gvec
  integer :: npwwfn                      ! Number of planewaves for wavefunctions
  integer :: npwx                        ! Number of planewaves for the exchange part
  integer :: npwc                        ! Number of planewaves for $\Sigma_c$ and W
  integer :: nbnds                       ! Number of bands kept in the calculation
  integer :: nkibz                       ! Number of k-points in the IBZ
  integer :: nqibz                       ! Number of q-points in the IBZ
  integer :: nkbz                        ! Number of k-points in the BZ
  integer :: nqbz                        ! Number of q-points in the BZ
  integer :: nsym                        ! Number of symmetry operations
                                         ! (operations related through the inversion symmetry are not considered)
  integer :: nsppol                      ! 1 for unpolarized, 2 for spin-polarized calculations
  integer :: time_reversal               ! 2 if time-reversal symmetry is used, 1 otherwise

  integer :: mG0(3)                      
   ! For each reduced direction gives the max G0 component to account for umklapp processes

 end type rdm_parameters
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/results_gs_type
!! NAME
!! results_gs_type
!!
!! FUNCTION
!! This structured datatype contains the results of a GS calculation :
!! energy and its decomposition, forces and their decompositions, stresses
!! and their decompositions
!!
!! SOURCE

 type results_gs_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom
   ! The number of atoms for this dataset

! Real (real(dp)) scalars

! All the energies are in Hartree, obtained "per unit cell".
  type(energies_type) :: energies
!!$  real(dp) :: eei      ! local pseudopotential energy (Hartree)
!!$  real(dp) :: eeig     ! sum of eigenvalue energy (Hartree)
!!$  real(dp) :: eew      ! Ewald energy (Hartree)
!!$  real(dp) :: ehart    ! Hartree part of total energy (Hartree)
!!$  real(dp) :: eii      ! pseudopotential core-core energy
!!$  real(dp) :: ek       ! kinetic energy (Hartree)
!!$  real(dp) :: enefield ! the term of the energy functional that depends
!!$                       ! explicitely on the electric field
!!$                       ! enefield = -ucvol*E*P
!!$  real(dp) :: enl      ! nonlocal pseudopotential energy (Hartree)
  real(dp) :: entropy  ! entropy (Hartree)
!!$  real(dp) :: enxc     ! exchange-correlation energy (Hartree)
!!$  real(dp) :: enxcdc   ! exchange-correlation double-counting energy (Hartree)
!!$  real(dp) :: epaw     ! PAW spherical energy (Hartree)
!!$  real(dp) :: epawdc   ! PAW spherical double-counting energy (Hartree)
  real(dp) :: etotal   ! total energy (Hartree)
                       ! for fixed occupation numbers (occopt==0,1,or 2):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl+PAW_spherical_part
                       ! for varying occupation numbers (occopt>=3):
                       !   etotal=ek+ehart+enxc+eei+eew+eii+enl - tsmear*entropy +PAW_spherical_part
  real(dp) :: fermie   ! Fermi energy (Hartree)
  real(dp) :: residm   ! maximum value for the residual over all bands, all k points,
                       !   and all spins (Hartree or Hartree**2, to be checked !)
  real(dp) :: vxcavg   ! Average of the exchange-correlation energy. The average
                       ! of the local psp pot and the Hartree pot is set to zero (due
                       ! to the usual problem at G=0 for Coulombic system, so vxcavg
                       ! is also the average of the local part of the Hamiltonian

! Real (real(dp)) arrays

  real(dp), pointer :: fcart(:,:)
   ! fcart(3,natom)
   ! Cartesian forces (Hartree/Bohr)
   ! Note : unlike fred, this array has been corrected by enforcing
   ! the translational symmetry, namely that the sum of force
   ! on all atoms is zero.

  real(dp), pointer :: fred(:,:)
   ! fred(3,natom)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp), pointer :: gresid(:,:)
   ! gresid(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the residual
   ! of the potential

  real(dp), pointer :: grewtn(:,:)
   ! grewtn(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the Ewald energy

  real(dp), pointer :: grxc(:,:)
   ! grxc(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the XC energy

  real(dp) :: pel(3)
   ! ucvol times the electronic polarization in reduced coordinates

  real(dp) :: strten(6)
   ! Stress tensor in cartesian coordinates (Hartree/Bohr^3)
   ! 6 unique components of this symmetric 3x3 tensor:
   ! Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).

  real(dp), pointer :: synlgr(:,:)
   ! synlgr(3,natom)
   ! Part of the gradient of the total energy (Hartree) with respect
   ! to change of reduced coordinates, that comes from the non-local energy
   ! The "sy" prefix refer to the fact that this gradient has been
   ! symmetrized.

 end type results_gs_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/results_out_type
!! NAME
!! results_out_type
!!
!! FUNCTION
!! This structured datatype contains a subset of the results of a GS
!! calculation, needed to perform the so-called "internal tests", and
!! to perform the timing analysis
!!
!! SOURCE

 type results_out_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: natom 
   ! The number of atoms for this dataset

! Integer arrays

  integer, pointer :: npwtot(:,:)      
   ! npw(mxnkpt,nimage) Full number of plane waves for each
   ! k point, computed with the "true" rprimd
   ! Not taking into account the decrease due to istwfk
   ! Not taking into account the spread of pws on different procs

! Real (real(dp)) scalars

  real(dp), pointer :: acell(:,:) 
   ! acell(3,nimage)  length of primitive vectors

  real(dp), pointer :: etotal(:)  
   ! etotal(nimage) total energy (Hartree)
   ! All the energies are in Hartree, obtained "per unit cell".

  real(dp), pointer :: fcart(:,:,:) 
   ! fcart(3,natom,nimage) Cartesian forces (Hartree/Bohr)

  real(dp), pointer :: fred(:,:,:)  
   ! fred(3,natom,nimage)
   ! Forces in reduced coordinates (Hartree)
   ! Actually, gradient of the total energy with respect
   ! to change of reduced coordinates

  real(dp), pointer :: occ(:,:)     
   ! occ(mxmband_upper*mxnkpt*mxnsppol,nimage)

  real(dp), pointer :: rprim(:,:,:)  
   ! rprim(3,3,nimage)

  real(dp), pointer :: rprimd(:,:,:) 
   ! rprim(3,3,nimage)

  real(dp), pointer :: strten(:,:)
   ! strten(6,nimage)

  real(dp), pointer :: vel(:,:,:)   
   ! vel(3,natom,nimage)

  real(dp), pointer :: xred(:,:,:)  
   ! xred(3,natom,nimage)

 end type results_out_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/scf_history_type
!! NAME
!! scf_history_type
!!
!! FUNCTION
!! This structured datatype contains various arrays obtained from
!! previous SCF cycles (density, positions...)
!!
!! SOURCE

 type scf_history_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: history_size
   ! Number of previous SCF cycles stored in history

  integer :: icall
   ! Number of call for the routine extraprho

  integer :: natom
   ! Number of atoms in cell

  integer :: nfft
   ! Size of FFT grid (for density)

  integer :: nspden
   ! Number of independant spin components for density

  integer :: usecg
   ! usecg=1 if the extrapolation of the wavefunctions is active

  real(dp) :: alpha
   ! alpha mixing coefficient for the prediction of density and wavefunctions

  real(dp) :: beta
   ! beta mixing coefficient for the prediction of density and wavefunctions

! Integer arrays

  integer,pointer :: hindex(:)
   ! hindex(history_size)
   ! Indexes of SCF cycles in the history
   ! hindex(1) is the newest SCF cycle
   ! hindex(history_size) is the oldest SCF cycle

! Real (real(dp)) arrays

   real(dp) :: rprimd(3,3)

   real(dp),pointer :: cg(:,:,:)
    ! cg(2,mcg,history_size)
    ! wavefunction coefficients needed for each SCF cycle of history

   real(dp),pointer :: deltarhor(:,:,:)
    ! deltarhor(nfft,nspden,history_size)
    ! Diference between electronic density (in real space)
    ! and sum of atomic densities at the end of each SCF cycle of history

   real(dp),pointer :: atmrho_last(:)
    ! atmrho_last(nfft)
    ! Sum of atomic densities at the end of the LAST SCF cycle

   real(dp),pointer :: xreddiff(:,:,:)
    ! xreddiff(3,natom,history_size)
    ! Difference of reduced coordinates of atoms between a
    ! SCF cycle and the previous

! Structured datatypes arrays

  type(pawrhoij_type), pointer :: pawrhoij(:,:)
    ! pawrhoij(natom,history_size)
    ! PAW only: occupancies matrix at the end of each SCF cycle of history

  type(cprj_type),pointer :: cprj(:,:,:)
    !cprj(natom,nspinor*mband*mkmem*nsppol,history_size)


 end type scf_history_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/sigma_results
!! NAME
!! sigma_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_results structured datatype
!! gather the results of a sigma calculation.
!!
!! SOURCE

 type sigma_results

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: b1gw,b2gw      ! min and Max gw band indeces over spin and k-points (used to dimension)
  integer :: gwcalctyp      ! Flag defining the calculation type.
  integer :: nkcalc         ! No. of points calculated
  integer :: nkibz          ! No. of irreducible k-points.
  integer :: nbnds          ! Total number of bands
  integer :: nomega_r       ! No. of real frequencies for the spectral function.
  integer :: nomega_i       ! No. of frequencies along the imaginary axis.
  integer :: nomega4sd      ! No. of real frequencies to evaluate the derivative of $\Sigma(E)$.
  integer :: nsig_ab        ! 1 if nspinor=1,4 for noncollinear case.
  integer :: nsppol         ! No. of spin polarizations.
  integer :: usepawu        ! 1 if we are using LDA+U as starting point (only for PAW)

  real(dp) :: deltae       ! Frequency step for the calculation of d\Sigma/dE
  real(dp) :: maxomega4sd  ! Max frequency around E_ks for d\Sigma/dE.
  real(dp) :: maxomega_r   ! Max frequency for spectral function.
  real(dp) :: scissor_ene  ! Scissor energy value. zero for None.

  integer,pointer :: maxbnd(:)   
  ! maxbnd(nkcalc)
  ! Max band index considered in GW for this k-point.

  integer,pointer :: minbnd(:)   
  ! minbnd(nkcalc)
  ! Min band index considered in GW for this k-point.

  !real(dp),pointer :: ame(:,:,:)
  ! ame(nbnds,nkibz,nomega))
  ! Diagonal matrix elements of the spectral function.
  ! Commented out, it can be calculated from the other quantities

  real(dp),pointer :: degwgap(:,:)   
  ! degwgap(nkibz,nsppol)
  ! Difference btw the QP and the KS optical gap.

  real(dp),pointer :: egwgap(:,:)   
  ! egwgap(nkibz,nsppol))
  ! QP optical gap at each k-point and spin.

  real(dp),pointer :: en_qp_diago(:,:,:)   
  ! en_qp_diago(nbnds,nkibz,nsppol))
  ! QP energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

  real(dp),pointer :: e0(:,:,:)    
  ! e0(nbnds,nkibz,nsppol)
  ! KS eigenvalues for each band, k-point and spin. In case of self-consistent?

  real(dp),pointer :: e0gap(:,:)   
  ! e0gap(nkibz,nsppol),
  ! KS gap at each k-point, for each spin.

  real(dp),pointer :: omega_r(:)   
  ! omega_r(nomega_r)
  ! real frequencies used for the self energy.

  real(dp),pointer :: xkcalc(:,:)  
  ! xkcalc(3,nkcalc)
  ! ! TODO this should be replaced by a table (nkibz)
  ! List of calculated k-points

  real(dp),pointer :: sigxme(:,:,:)  
  ! sigxme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal matrix elements of $\Sigma_x$ i.e $\<nks|\Sigma_x|nks\>$

  real(dp),pointer :: vxcme(:,:,:)  
  ! vxcme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{xc}[n_val]|nks\>$ matrix elements of vxc (valence-only contribution).

  real(dp),pointer :: vUme(:,:,:)   
  ! vUme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{U}|nks\>$ for LDA+U.

  complex(dpc),pointer :: degw(:,:,:)   
  ! degw(b1gw:b2gw,nkibz,nsppol))
  ! Difference between the QP and the KS energies.

  complex(dpc),pointer :: dsigmee0(:,:,:)  
  ! dsigmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

  complex(dpc),pointer :: egw(:,:,:)  
  ! degw(nbnds,nkibz,nsppol))
  ! QP energies, $\epsilon_{nks}^{QP}$.

  complex(dpc),pointer :: eigvec_qp(:,:,:,:)   
  ! eigvec_qp(nbnds,nbnds,nkibz,nsppol))
  ! Expansion of the QP amplitude in the KS basis set.

  complex(dpc),pointer :: hhartree(:,:,:,:)   
  ! hhartree(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab)
  ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$

  complex(dpc),pointer :: sigcme(:,:,:,:)   
  ! sigcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

  complex(dpc),pointer :: sigmee(:,:,:)  
  ! sigmee(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\Sigma_{xc}E_{KS} + (E_{QP}- E_{KS})*dSigma/dE_KS

  complex(dpc),pointer :: sigcmee0(:,:,:)   
  ! sigcmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal mat. elements of $\Sigma_c(E)$ calculated at the KS energy $E_{KS}$

  complex(dpc),pointer :: sigcmesi(:,:,:,:)   
  ! sigcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_c$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),pointer :: sigcme4sd(:,:,:,:)   
  ! sigcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_c around the zeroth order eigenvalue (usually KS).

  complex(dpc),pointer :: sigxcme(:,:,:,:)   
  ! sigxme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{xc}(E)|nks\>$ at each real frequency frequency.

  complex(dpc),pointer :: sigxcmesi(:,:,:,:)   
  ! sigxcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_{xc}$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),pointer :: sigxcme4sd(:,:,:,:)   
  ! sigxcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_xc for frequencies around the zeroth order eigenvalues.

  complex(dpc),pointer :: ze0(:,:,:)   
  ! ze0(b1gw:b2gw,nkibz,nsppol))
  ! renormalization factor. $(1-\dfrac{\partial\Sigma_c} {\partial E_{KS}})^{-1}$

  complex(dpc),pointer :: omega_i(:)  
  ! omegasi(nomega_i)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dpc),pointer :: omega4sd(:,:,:,:)  
  ! omega4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol).
  ! Frequencies used to evaluate the Derivative of Sigma.

 end type sigma_results
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/vardims_type
!! NAME
!!  vardims_type
!!
!! FUNCTION
!!  Stores dimensions of dataset variables.
!!
!! SOURCE

 type vardims_type

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: mband 
  integer :: mproj
  integer :: mpsang
  integer :: mpw
  integer :: ngrid1
  integer :: ngrid2 
  integer :: ngrid3
  integer :: ntypat
  integer :: natom
  integer :: natsph
  integer :: nkpt 
  integer :: nkptgw
  integer :: nshiftk
  integer :: nsppol
  integer :: nberry
  integer :: nsym
  integer :: npsp
  integer :: nconeq
  integer :: ntypalch
  integer :: npspalch
  integer :: nfft
  integer :: nspden
  integer :: wfs_dim1
  integer :: wfs_dim2
  integer :: nfreqsus
  integer :: npw_tiny
  integer :: nqptdm
  integer :: nspinor

 end type vardims_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/ScrHdr_type
!! NAME
!!  ScrHdr_type
!!
!! FUNCTION
!!  For the GW part of ABINIT, this structure defines the header
!!  of the _SCR or _SUSC file.
!!  FIXME:
!!   This is just an hack to cheat the build system.
!!   It should be defined in the appropriated module (m_io_screening.F90)
!!
!! SOURCE


!----------------------------------------------------------------------

 type ScrHdr_type

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  !Other variables that can be added are, for the moment, commented out.
  !Most of them are related to the Abinit implementation  and are not specified in the ETSF specs.

  !Index of the qlwl section?
  !gwcomp, gwencomp  ! Info on the extrapolar algorithm

  integer :: ID           ! Matrix identifier: O if not yet defined, 1 for chi0,
                          ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}
  integer :: ikxc         ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for frequency-dependent TDDFT
  integer :: inclvkb      ! q-->0 treatment, 0 for None, 1-2 for transversal gauge, 3 for longitudinal
  integer :: headform     ! format of the SCR header
  integer :: fform        ! File format:
  integer :: gwcalctyp    ! Calculation type (G0W0, G0W, GW ...)
  integer :: nI,nJ        ! Number of spin components (rows,columns) in chi|eps^-1. (1,1) if collinear.
                          !  The internal representation of the matrix is eps(nI*npwe,nJ*npwe)
  integer :: nqibz        ! Number of q-points in the IBZ.
  integer :: nqlwl        ! Number of points for the treatment of the long wavelength limit.
  integer :: nomega       ! Total number of frequencies.
  integer :: nbnds_used   ! Number of bands used during the screening calculation (only for info)
  integer :: npwe         ! Number of G vectors reported on the file.
  integer :: npwwfn_used  ! Number of G vectors for wavefunctions used during the screening calculation (only for info)
  integer :: spmeth       ! Method used to approximate the delta function in the expression for Im Chi_0
  integer :: test_type    ! 0 for None, 1 for TEST-PARTICLE, 2 for TEST-ELECTRON (only for TDDFT)
  integer :: tordering    ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  real(dp) :: soenergy    ! Scissor Energy, zero if not used
  real(dp) :: spsmear     ! Smearing of the delta in case of spmeth==2
  real(dp) :: zcut        ! Imaginary shift to avoid the poles along the real axis.

  type(Hdr_type) :: Hdr   ! The abinit header.

!arrays
  character(len=80) :: title(2)
  ! Title describing the content of the file.

  integer,pointer :: gvec(:,:)   
  ! gvec(3,npwe)
  ! G vectors in r.l.u.

  real(dp),pointer :: qibz(:,:)  
  ! qibz(3,nqibz)
  ! q-points in r.l.u.

  real(dp),pointer :: qlwl(:,:)  
  ! qlwl(3,nqlwl)
  ! q-points for the long wave-length limit treatment (r.l.u)

  complex(dpc),pointer :: lwing(:,:,:)  
  ! lwing(npwe,nomega,nqlwl)
  ! Lower wings for the different q"s -->0

  complex(dpc),pointer :: omega(:)   
  ! omega(nomega)
  ! All frequencies calculated both along the real and the imaginary axis.

  complex(dpc),pointer :: uwing(:,:,:)   
  ! uwing(npwe,nomega,nqlwl)
  ! Upper wings for the different q"s -->0

 end type ScrHdr_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/epsilonm1_results
!! NAME
!! epsilonm1_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the epsilonm1_results structured datatype
!! gather the results of screening : the inverse dielectric matrix,
!! and the omega matrices .
!!
!! SOURCE

 type Epsilonm1_results

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: ID                          ! Matrix identifier: O if not yet defined, 1 for chi0,
                                         ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}, 5 for W.
  integer :: ikxc                        ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for TDDFT
  integer :: fform                       ! File format: 1002 for SCR|SUSC files.
  integer :: mqmem                       ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.
  integer :: nI,nJ                       ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.
  integer :: nqibz                       ! Number of q-points in the IBZ used.
  integer :: nqlwl                       ! Number of point used for the treatment of the long wave-length limit.
  integer :: nomega                      ! Number of frequencies used.
  integer :: nomega_i                    ! Number of purely imaginary frequencies used.
  integer :: nomega_r                    ! Number of real frequencies used.
  integer :: npwe                        ! Number of G vectors used.
  integer :: test_type                   ! 0 for None, 1 for TEST-PARTICLE, 2 for TEST-ELECTRON (only for TDDFT)
  integer :: Tordering                   ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  character(len=fnlen) :: fname          ! Name of the file from which epsm1 is read.

!arrays
  integer,pointer  :: gvec(:,:)    
  ! gvec(3,npwe)
  ! G-vectors used to describe the two-point function (r.l.u.).

  real(dp),pointer :: qibz(:,:)  
  ! qibz(3,nqibz)
  ! q-points in reduced coordinates

  real(dp),pointer :: qlwl(:,:) 
  ! qlwl(3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  complex(gwpc),pointer :: epsm1(:,:,:,:)  
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(dpc),pointer :: lwing(:,:,:)  
  ! lwing(npwe,nomega,nqlwl)
  ! Lower wings for the different q"s -->0

  complex(dpc),pointer :: omega(:)  
  ! omega(nomega)
  ! Frequencies used both along the real and the imaginary axis.

  complex(dpc),pointer :: uwing(:,:,:)  
  ! uwing(npwe,nomega,nqlwl)
  ! Upper wings for the different q"s -->0

  type(ScrHdr_type) :: Hscr
  ! The header reported in the _SCR of _SUSC file.
  ! This object contains information on the susceptibility or the inverse dielectric matrix
  ! as stored in the external file. These quantities do *NOT* correspond to the quantities
  ! used during the GW calculation since some parameters might differ, actually they might be smaller.
  ! For example, the number of G-vectors used can be smaller than the number of G"s stored on file.

 end type Epsilonm1_results
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/Jpt_gwpc_2D
!! NAME
!! Jpt_gwpc_2D
!!
!! FUNCTION
!!  The Jpt_gwpc_2D data type defines a polymorphic object that extends the concept 
!!  of pointer and allocatable variable. It is mainly used in the GW part 
!!  in order to minimize the memory requirements. The object has two different possible status:
!!   * true pointer storing the address in memory of the pointee.
!!   * array allocated on the heap (same as an allocatable variable).
!!
!! SOURCE

 type :: Jpt_gwpc_2D

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: stat = 1 !!$=JPT_ISPOINTER  =1
  complex(gwpc),pointer :: datum(:,:)   
 end type Jpt_gwpc_2D
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/Jpt_gwpc_3D
!! NAME
!! Jpt_gwpc_3D
!!
!! FUNCTION
!!  The Jpt_gwpc_3D data type defines a polymorphic object that extends the concept 
!!  of pointer and allocatable variable. It is mainly used in the GW part 
!!  in order to minimize the memory requirements. The object has two different possible status:
!!   * true pointer storing the address in memory of the pointee.
!!   * array allocated on the heap (same as an allocatable variable).
!!
!! SOURCE

 type :: Jpt_gwpc_3D

! WARNING : if you modify this datatype, please check there there is no creation/destruction routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: stat = 1 !!$=JPT_ISPOINTER  =1
  complex(gwpc),pointer :: datum(:,:,:)  
 end type Jpt_gwpc_3D
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/PPmodel_type
!! NAME
!! PPmodel_type
!!
!! FUNCTION
!!  For the GW part of ABINIT, this structure datatype gathers all
!!  the information on the Plasmonpole technique used in the calculations
!!
!! SOURCE

 type PPmodel_type

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

 !integers
   integer :: dm2_botsq                          ! =npwc if ppmodel=1,2; =1 if ppmodel=3,4
   integer :: dm_eig                             ! =npwc if ppmodel=3;   =0 if ppmodel=1,2,4
   integer :: dm2_otq                            ! =npwc if ppmodel=1,2; =1 if ppmodel=3,4
   integer :: model                              ! The type of Plasmonpole model
   integer :: mqmem                              ! =nqibz if in-core solutions, =0 for out-of-core
                                                 ! (In the former case the kast dimension in PPm arrays has size 1)
   !integer :: needs_rhog                        ! 1 if the ppmodel requires rho(G).
   integer :: nqibz                              ! Number of q-points in the IBZ
   integer :: npwc                               ! Number of G vectors in $\tilde \epsilon $

   real(dp) :: drude_plsmf                       ! Drude plasma frequency
   logical :: save_memory_devel=.FALSE.          

   ! Use in-place storage of the PPm parameters, should be FALSE for non-developers.
   !!$real(dp) :: zcut
   !!$real(dp),pointer :: qibz(:,:)

   !logical :: has_inversion
   !logical :: has_time_reversal

 !arrays
   type(Jpt_gwpc_3D) :: bigomegatwsq
   ! bigomegatwsq(npwc,dm2_botsq,nqibz)
   ! Plasmon pole parameters $\tilde\Omega^2_{G Gp}(q)$.

   type(Jpt_gwpc_3D) :: omegatw
   ! omegatw(npwc,dm2_otq,nqibz)
   ! Plasmon pole parameters $\tilde\omega_{G Gp}(q)$.

   complex(gwpc),pointer :: eigpot(:,:,:)   
   ! eigpot(dm_eig,dm_eig,nqibz)
   ! Eigvectors of the symmetrized inverse dielectric matrix
 end type
!!***

end module defs_datatypes
!!***

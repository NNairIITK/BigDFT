!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2015 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which contains the Fortran data structures
!! and the routines of allocations and de-allocations
module module_types

  use m_ab7_mixing, only : ab7_mixing_object
  use module_base!, only : gp,wp,dp,tp,uninitialized,mpi_environment,mpi_environment_null,&
  !bigdft_mpi,ndebug,memocc!,vcopy
  use module_xc, only : xc_info
  use gaussians, only: gaussian_basis
  use Poisson_Solver, only: coulomb_operator
  use locregs
  use psp_projectors_base
  use module_atoms, only: atoms_data,symmetry_data,atomic_structure
  use module_dpbox, only: denspot_distribution,dpbox_null
  use communications_base, only: comms_linear, comms_cubic, p2pComms
  use sparsematrix_base, only: matrices, sparse_matrix
  use foe_base, only: foe_data
  use m_pawcprj, only: pawcprj_type
  use m_paw_an, only: paw_an_type
  use m_paw_ij, only: paw_ij_type
  use m_pawfgrtab, only: pawfgrtab_type
  use m_pawrhoij, only: pawrhoij_type
  use module_input_keys, only: SIC_data,orthon_data,input_variables
  use fragment_base, only: fragmentInputParameters

  implicit none

  private

  !> Contains all energy terms
  type, public :: energy_terms
     real(gp) :: eh      !< Hartree energy
     real(gp) :: exc     !< Exchange-correlation energy
     real(gp) :: evxc    !< Energy from the exchange-correlation potential
     real(gp) :: eion    !< Ion-Ion interaction
     real(gp) :: edisp   !< Dispersion force energy
     real(gp) :: ekin    !< Kinetic term
     real(gp) :: epot    !< local potential energy
     real(gp) :: eproj   !< energy of PSP projectors
     real(gp) :: eexctX  !< exact exchange energy
     real(gp) :: eelec   !< electrostatic energy. Replaces the hartree energy for cavities
     real(gp) :: ebs     
     real(gp) :: eKS     
     real(gp) :: trH     !< Trace of Hamiltonian i.e. band structure 
     real(gp) :: evsum   
     real(gp) :: evsic   
     real(gp) :: excrhoc 
     real(gp) :: eTS     
     real(gp) :: ePV     !< Pressure term
     real(gp) :: energy  !< The functional which is minimized
     real(gp) :: e_prev  !< The previous value, to show the delta
     real(gp) :: trH_prev!< The previous value, to show the delta
     !real(gp), dimension(:,:), pointer :: fion,f

     integer(kind = 8) :: c_obj !< Storage of the C wrapper object.
  end type energy_terms


  !> Memory estimation requirements
  type, public :: memory_estimation
     double precision :: submat
     integer :: ncomponents, norb, norbp
     double precision :: oneorb, allpsi_mpi, psistorage
     double precision :: projarr
     double precision :: grid
     double precision :: workarr

     double precision :: kernel, density, psolver, ham

     double precision :: peak
  end type memory_estimation


  !> Used to split between points to be treated in simple or in double precision
  type, public :: rho_descriptors
     character(len=1) :: geocode !< @copydoc poisson_solver::doc::geocode
     integer :: icomm !< method for communicating the density
     integer :: nrhotot !< dimension of the partial density array before communication
     integer :: n_csegs,n_fsegs,dp_size,sp_size
     integer, dimension(:,:), pointer :: spkey,dpkey
     integer, dimension(:), pointer :: cseg_b,fseg_b
  end type rho_descriptors


  !> Structures of basis of gaussian functions of the form exp(-a*r2)cos/sin(b*r2)
  type, public :: gaussian_basis_c
     integer :: nat,ncoeff,nshltot,nexpo
     integer, dimension(:), pointer :: nshell,ndoc,nam
     complex(gp), dimension(:), pointer :: expof,psiat
     real(gp), dimension(:,:), pointer :: rxyz
  end type gaussian_basis_c


  !> All the parameters which are important for describing the orbitals
  !! Add also the objects related to k-points sampling, after symmetries applications
  type, public :: orbitals_data 
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: norbup        !< Total number of orbitals if there were only up orbitals
     integer :: norbdp        !< probably to be deleted...
     integer :: norbu,norbd,nspin,nspinor,isorb,isorbu,isorbd
     integer :: nkpts,nkptsp,iskpts
     real(gp) :: efermi,HLgap,eTS
     integer, dimension(:), pointer :: iokpt,ikptproc,isorb_par,ispot
     integer, dimension(:), pointer :: inwhichlocreg,onwhichatom
     integer, dimension(:,:), pointer :: norb_par, norbu_par, norbd_par
     real(wp), dimension(:), pointer :: eval
     real(gp), dimension(:), pointer :: occup,spinsgn,kwgts
     real(gp), dimension(:,:), pointer :: kpts
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  end type orbitals_data


  !> Contains the pointers to be handled to control GPU information
  !! Given that they are pointers on GPU address, they are C pointers
  !! which take 8 bytes
  !! So they are declared as kind=8 variables either if the GPU works in simple precision
  !! Also other information concerning the GPU runs can be stored in this structure
  type, public :: GPU_pointers
     logical :: OCLconv  !< True if OCL is used for convolutions for this MPI process
     logical :: useDynamic,full_locham
     integer :: id_proc,ndevices
     real(kind=8) :: keys,work1,work2,work3,rhopot,r,d
     real(kind=8) :: rhopot_down, rhopot_up
     real(kind=8) :: work1_i,work2_i,work3_i,d_i
     real(kind=8) :: pinned_in,pinned_out
     real(kind=8), dimension(:), pointer :: psi
     real(kind=8) :: psi_c,psi_f
     real(kind=8) :: psi_c_i,psi_f_i
     real(kind=8) :: psi_c_r,psi_f_r,psi_c_b,psi_f_b,psi_c_d,psi_f_d
     real(kind=8) :: psi_c_r_i,psi_f_r_i,psi_c_b_i,psi_f_b_i,psi_c_d_i,psi_f_d_i
     real(kind=8) :: keyg_c,keyg_f,keyv_c,keyv_f
     real(kind=8) :: keyg_c_host,keyg_f_host,keyv_c_host,keyv_f_host
     real(kind=8) :: context,queue
     !host pointers to be freed
     real(kind=8) :: rhopot_down_host, rhopot_up_host
     real(kind=8), dimension(:,:,:), pointer :: ekinpot_host
     real(kind=8), dimension(:,:), pointer :: psicf_host
     real(kind=8), dimension(:,:), pointer :: hpsicf_host
     real(kind=8), dimension(:), pointer :: bprecond_host

     real(gp), dimension(:,:), pointer :: ekin, epot !< values of the kinetic and potential energies to be passed to local_hamiltonian
     real(wp), dimension(:), pointer :: hpsi_ASYNC !<pointer to the wavefunction allocated in the case of asyncronous local_hamiltonian
  end type GPU_pointers


  !> Contains all the descriptors necessary for splitting the calculation in different locregs 
  type, public :: local_zone_descriptors
     logical :: linear                                      !< if true, use linear part of the code
     integer :: nlr                                         !< Number of localization regions 
     integer :: lintyp                                      !< If 0 cubic, 1 locreg and 2 TMB
     integer :: ndimpotisf                                  !< Total dimension of potential in isf (including exctX)
     real(gp), dimension(3) :: hgrids                       !< Grid spacings of wavelet grid (coarser resolution)
     type(locreg_descriptors) :: Glr                        !< Global region descriptors
     type(locreg_descriptors), dimension(:), pointer :: Llr !< Local region descriptors (dimension = nlr)
  end type local_zone_descriptors


  !> Contains the work arrays needed for expressing wavefunction in real space
  !! with all the BC
  type, public :: workarr_sumrho
     integer :: nw1,nw2,nxc,nxf
     real(wp), dimension(:), pointer :: x_c,x_f,w1,w2
  end type workarr_sumrho


  !> Contains the work arrays needed for hamiltonian application with all the BC
  type, public :: workarr_locham
     integer :: nw1,nw2,nxc,nyc,nxf1,nxf2,nxf3,nxf,nyf
     real(wp), dimension(:), pointer :: w1,w2
     !for the periodic BC case, these arrays substitute 
     !psifscf,psifscfk,psig,ww respectively
     real(wp), dimension(:,:), pointer :: x_c,y_c,x_f1,x_f2,x_f3,x_f,y_f
  end type workarr_locham


  !> Contains the work arrays needed for th preconditioner with all the BC
  !! Take different pointers depending on the boundary conditions
  type, public :: workarr_precond
     integer, dimension(:), pointer :: modul1,modul2,modul3
     real(wp), dimension(:), pointer :: psifscf,ww,x_f1,x_f2,x_f3,kern_k1,kern_k2,kern_k3
     real(wp), dimension(:,:), pointer :: af,bf,cf,ef
     real(wp), dimension(:,:,:), pointer :: xpsig_c,ypsig_c,x_c
     real(wp), dimension(:,:,:,:), pointer :: xpsig_f,ypsig_f,x_f,y_f
     real(wp), dimension(:,:,:,:,:), pointer :: z1,z3 ! work array for FFT
  end type workarr_precond


  !!> Fermi Operator Expansion parameters
  !type, public :: foe_data
  !  integer :: nseg
  !  integer,dimension(:),pointer :: nsegline, istsegline
  !  integer,dimension(:,:),pointer :: keyg
  !  real(kind=8) :: ef                     !< Fermi energy for FOE
  !  real(kind=8) :: evlow, evhigh          !< Eigenvalue bounds for FOE 
  !  real(kind=8) :: bisection_shift        !< Bisection shift to find Fermi energy (FOE)
  !  real(kind=8) :: fscale                 !< Length scale for complementary error function (FOE)
  !  real(kind=8) :: ef_interpol_det        !< FOE: max determinant of cubic interpolation matrix
  !  real(kind=8) :: ef_interpol_chargediff !< FOE: max charge difference for interpolation
  !  real(kind=8) :: charge                 !< Total charge of the system
  !  real(kind=8) :: fscale_lowerbound      !< lower bound for the error function decay length
  !  real(kind=8) :: fscale_upperbound       !< upper bound for the error function decay length
  !  integer :: evbounds_isatur, evboundsshrink_isatur, evbounds_nsatur, evboundsshrink_nsatur !< variables to check whether the eigenvalue bounds might be too big
  !end type foe_data

  
  type,public :: linear_matrices
      type(sparse_matrix) :: s !< small: sparsity pattern given by support function cutoff
      type(sparse_matrix) :: m !< medium: sparsity pattern given by SHAMOP cutoff
      type(sparse_matrix) :: l !< medium: sparsity pattern given by kernel cutoff
      type(sparse_matrix),dimension(:),pointer :: ks !< sparsity pattern for the KS orbitals (i.e. dense); spin up and down
      type(sparse_matrix),dimension(:),pointer :: ks_e !< sparsity pattern for the KS orbitals including extra stated (i.e. dense); spin up and down
      type(matrices) :: ham_, ovrlp_, kernel_
      type(matrices),dimension(3) :: ovrlppowers_
  end type linear_matrices


  type, public :: workarrays_quartic_convolutions
    real(wp), dimension(:,:,:), pointer :: xx_c, xy_c, xz_c
    real(wp), dimension(:,:,:), pointer :: xx_f1
    real(wp), dimension(:,:,:), pointer :: xy_f2
    real(wp), dimension(:,:,:), pointer :: xz_f4
    real(wp), dimension(:,:,:,:), pointer :: xx_f, xy_f, xz_f
    real(wp), dimension(:,:,:), pointer :: y_c
    real(wp), dimension(:,:,:,:), pointer :: y_f
    ! The following arrays are work arrays within the subroutine
    real(wp), dimension(:,:), pointer :: aeff0array, beff0array, ceff0array, eeff0array
    real(wp), dimension(:,:), pointer :: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
    real(wp), dimension(:,:), pointer :: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
    real(wp), dimension(:,:,:), pointer :: xya_c, xyc_c
    real(wp), dimension(:,:,:), pointer :: xza_c, xzc_c
    real(wp), dimension(:,:,:), pointer :: yza_c, yzb_c, yzc_c, yze_c
    real(wp), dimension(:,:,:,:), pointer :: xya_f, xyb_f, xyc_f, xye_f
    real(wp), dimension(:,:,:,:), pointer :: xza_f, xzb_f, xzc_f, xze_f
    real(wp), dimension(:,:,:,:), pointer :: yza_f, yzb_f, yzc_f, yze_f
  end type workarrays_quartic_convolutions


  type,public :: work_mpiaccumulate
    integer :: ncount
    real(wp),dimension(:),pointer :: receivebuf
    real(wp),dimension(:),pointer :: sendbuf
    integer :: window
  end type work_mpiaccumulate


  type, public :: localizedDIISParameters
    integer :: is, isx, mis, DIISHistMax, DIISHistMin
    integer :: icountSDSatur, icountDIISFailureCons, icountSwitch, icountDIISFailureTot, itBest
    real(kind=8), dimension(:), pointer :: phiHist, hphiHist, energy_hist
    real(kind=8) :: alpha_coeff !step size for optimization of coefficients
    real(kind=8), dimension(:,:,:), pointer :: mat
    real(kind=8) :: trmin, trold, alphaSD, alphaDIIS
    logical :: switchSD, immediateSwitchToSD, resetDIIS
  end type localizedDIISParameters


  type, public :: mixrhopotDIISParameters
    integer :: is, isx, mis
    real(kind=8), dimension(:), pointer :: rhopotHist, rhopotresHist
    real(kind=8), dimension(:,:), pointer :: mat
  end type mixrhopotDIISParameters


  !> Contains the arguments needed for the diis procedure
  type, public :: diis_objects
     logical :: switchSD
     integer :: idiistol,mids,ids,idsx
     real(gp) :: energy_min,energy_old,energy,alpha,alpha_max
     real(tp), dimension(:), pointer :: psidst
     real(tp), dimension(:), pointer :: hpsidst
     real(tp), dimension(:,:,:,:,:,:), pointer :: ads
  end type diis_objects


  !> Contains the information needed for the preconditioner
  type, public :: precond_data
    integer :: confPotOrder                                !< The order of the algebraic expression for Confinement potential
    integer :: ncong                                       !< Number of CG iterations for the preconditioning equation
    logical, dimension(:), pointer :: withConfPot          !< Use confinement potentials
    real(kind=8), dimension(:), pointer :: potentialPrefac !< Prefactor for the potentiar : Prefac * f(r) 
  end type precond_data


  !> Information for the confining potential to be used in TMB scheme
  !! The potential is supposed to be defined as prefac*(r-rC)**potorder
  type, public :: confpot_data
     integer :: potorder                !< Order of the confining potential
     integer, dimension(3) :: ioffset   !< Offset for the coordinates of potential lr in global region
     real(gp) :: prefac                 !< Prefactor
     real(gp), dimension(3) :: hh       !< Grid spacings in ISF grid
     real(gp), dimension(3) :: rxyzConf !< Confining potential center in global coordinates
     real(gp) :: damping                !< Damping factor to be used after the restart
  end type confpot_data


  !> Defines the important information needed to reformat a old wavefunctions
  type, public :: old_wavefunction
     type(local_zone_descriptors) :: Lzd       !< Local zone descriptors of the corresponding run
     real(wp), dimension(:), pointer :: psi    !< Wavelets coefficients in compressed form
     real(gp), dimension(:,:), pointer :: rxyz !< Atomic positions of the step
  end type old_wavefunction


  !> Densities and potentials, and related metadata, needed for their creation/application
  !! Not all these quantities are available, some of them may point to the same memory space
  type, public :: DFT_local_fields
     real(dp), dimension(:), pointer :: rhov          !< Generic workspace. What is there is indicated by rhov_is
     type(ab7_mixing_object), pointer :: mix          !< History of rhov, allocated only when using diagonalisation
     !> Local fields which are associated to their name
     !! normally given in parallel distribution
     real(dp), dimension(:,:), pointer :: rho_psi     !< density as given by square of el. WFN
     real(dp), dimension(:,:,:,:), pointer :: rho_C   !< core density
     real(wp), dimension(:,:,:,:), pointer :: V_ext   !< local part of pseudopotientials
     real(wp), dimension(:,:,:,:), pointer :: V_XC    !< eXchange and Correlation potential (local)
     real(wp), dimension(:,:,:,:), pointer :: Vloc_KS !< complete local potential of KS Hamiltonian (might point on rho_psi)
     real(wp), dimension(:,:,:,:), pointer :: f_XC    !< dV_XC[rho]/d_rho
     real(wp), dimension(:,:,:,:), pointer :: rho_ion !< charge density of the ions, to be passed to PSolver
     !temporary arrays
     real(wp), dimension(:), pointer :: rho_work,pot_work !< Full grid arrays
     !metadata
     integer :: rhov_is
     real(gp) :: psoffset                 !< Offset of the Poisson Solver in the case of Periodic BC
     type(rho_descriptors) :: rhod        !< Descriptors of the density for parallel communication
     type(denspot_distribution) :: dpbox  !< Distribution of density and potential box
     type(xc_info) :: xc                  !< Structure about the used xc functionals
     character(len=3) :: PSquiet
     !real(gp), dimension(3) :: hgrids    !< Grid spacings of denspot grid (half of the wvl grid)
     type(coulomb_operator) :: pkernel    !< Kernel of the Poisson Solver used for V_H[rho]
     type(coulomb_operator) :: pkernelseq !< For monoproc PS (useful for exactX, SIC,...)

     integer(kind = 8) :: c_obj = 0       !< Storage of the C wrapper object.
  end type DFT_local_fields

  !> Check if all comms are necessary here
  type, public :: hamiltonian_descriptors
     integer :: npsidim_orbs             !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp             !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd !< Data on the localisation regions, if associated
     type(comms_linear) :: collcom       !< describes collective communication
     type(p2pComms) :: comgp             !< Describing p2p communications for distributing the potential
     real(wp), dimension(:), pointer :: psi,psit_c,psit_f !< these should eventually be eliminated
     logical :: can_use_transposed
  end type hamiltonian_descriptors

  !> Contains the arguments needed for the PAW implementation:
  !> to be better integrated into the other structures.
  type, public :: paw_objects
     logical :: usepaw
     integer :: lmnmax
     integer :: ntypes
     integer :: natom
     type(paw_an_type), dimension(:), pointer :: paw_an
     type(paw_ij_type), dimension(:), pointer :: paw_ij
     type(pawcprj_type), dimension(:,:), pointer :: cprj
     type(pawfgrtab_type), dimension(:), pointer :: pawfgrtab
     type(pawrhoij_type), dimension(:), pointer :: pawrhoij

     real(wp), dimension(:), pointer :: spsi !< Metric operator applied to psi (To be used for PAW)
  end type paw_objects

  !> The wavefunction which have to be considered at the DFT level
  type, public :: DFT_wavefunction
     !coefficients
     real(wp), dimension(:), pointer :: psi,hpsi,psit,psit_c,psit_f !< orbitals, or support functions, in wavelet basis
     real(wp), dimension(:,:), pointer :: gaucoeffs                 !< orbitals in gbd basis
     !basis sets
     type(gaussian_basis) :: gbd         !< gaussian basis description, if associated
     type(local_zone_descriptors) :: Lzd !< data on the localisation regions, if associated
     !restart objects (consider to move them in rst structure)
     type(old_wavefunction), dimension(:), pointer :: oldpsis !< previously calculated wfns
     integer :: istep_history                                 !< present step of wfn history
     !data properties
     logical :: can_use_transposed                           !< true if the transposed quantities are allocated and can be used
     type(orbitals_data) :: orbs                             !< wavefunction specification in terms of orbitals
     type(comms_cubic) :: comms                              !< communication objects for the cubic approach
     type(diis_objects) :: diis                              !< DIIS objects
     type(confpot_data), dimension(:), pointer :: confdatarr !< data for the confinement potential
     type(SIC_data) :: SIC                                   !< control the activation of SIC scheme in the wavefunction
     type(paw_objects) :: paw                                !< PAW objects
     type(orthon_data) :: orthpar                            !< control the application of the orthogonality scheme for cubic DFT wavefunction
     character(len=4) :: exctxpar                            !< Method for exact exchange parallelisation for the wavefunctions, in case
     type(p2pComms) :: comgp                                 !< describing p2p communications for distributing the potential
     type(comms_linear) :: collcom                           !< describes collective communication
     type(comms_linear) :: collcom_sr                        !< describes collective communication for the calculation of the charge density
     integer(kind = 8) :: c_obj                              !< Storage of the C wrapper object. it has to be initialized to zero
     type(foe_data) :: foe_obj                               !< describes the structure of the matrices for the linear method foe
     type(linear_matrices) :: linmat
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(hamiltonian_descriptors) :: ham_descr
     real(kind=8), dimension(:,:), pointer :: coeff          !< Expansion coefficients
     real(kind=8) :: damping_factor_confinement !< damping for the confinement after a restart
  end type DFT_wavefunction

  !> Used to control the optimization of wavefunctions
  type, public :: DFT_optimization_loop
     integer :: iscf !< Kind of optimization scheme.

     integer :: itrpmax !< specify the maximum number of mixing cycle on potential or density
     integer :: nrepmax !< specify the maximum number of restart after re-diagonalization
     integer :: itermax !< specify the maximum number of minimization iterations, self-consistent or not
     integer :: itermin !< specify the minimum number of minimization iterations, self-consistent or not !Bastian

     integer :: itrp    !< actual number of mixing cycle.
     integer :: itrep   !< actual number of re-diagonalisation runs.
     integer :: iter    !< actual number of minimization iterations.

     integer :: infocode !< return value after optimization loop.
                         !! - 0 run successfully succeded
                         !! - 1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
                         !!     forces may be meaningless   
                         !! - 2 (present only for inputPsiId=INPUT_PSI_MEMORY_WVL) gnrm of the first iteration > 1 AND growing in
                         !!     the second iteration OR grnm 1st >2.
                         !!     Input wavefunctions need to be recalculated.
                         !! - 3 (present only for inputPsiId=INPUT_PSI_LCAO) gnrm > 4. SCF error.

     real(gp) :: gnrm   !< actual value of cv criterion of the minimization loop.
     real(gp) :: rpnrm  !< actual value of cv criterion of the mixing loop.

     real(gp) :: gnrm_cv       !< convergence criterion of the minimization loop.
     real(gp) :: rpnrm_cv      !< convergence criterion of the mixing loop.
     real(gp) :: gnrm_startmix !< gnrm value to start mixing after.

     integer(kind = 8) :: c_obj = 0 !< Storage of the C wrapper object.
  end type DFT_optimization_loop


 !> Define generic subroutine

 !> Timing categories
 character(len=*), parameter, private :: tgrp_pot='Potential'
 integer, save, public :: TCAT_EXCHANGECORR=TIMING_UNINITIALIZED
 integer, parameter, private :: ncls_max=6,ncat_bigdft=153   ! define timimg categories and classes
 character(len=14), dimension(ncls_max), parameter, private :: clss = (/ &
      'Communications'    ,  &
      'Convolutions  '    ,  &
      'Linear Algebra'    ,  &
      'Other         '    ,  &
!      'Potential     '    ,  &
      'Initialization'    ,  &
      'Finalization  '    /)
 character(len=14), dimension(3,ncat_bigdft), parameter, private :: cats = reshape((/ &
                               !       Name           Class       Operation Kind
      'ReformatWaves ','Initialization' ,'Small Convol  ' ,  &  !< Reformatting of input waves
      'CrtDescriptors','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of descriptor arrays
      'CrtLocPot     ','Initialization' ,'Miscellaneous ' ,  &  !< Calculation of local potential
      'CrtProjectors ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of projectors
      'CrtPcProjects ','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of preconditioning projectors
      'CrtPawProjects','Initialization' ,'RMA Pattern   ' ,  &  !< Calculation of abscalc-pawprojectors
      'ApplyLocPotKin','Convolutions  ' ,'OpenCL ported ' ,  &  !< Application of PSP, kinetic energy
      'ApplyProj     ','Other         ' ,'RMA pattern   ' ,  &  !< Application of nonlocal PSP
      'Precondition  ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Precondtioning
      'Rho_comput    ','Convolutions  ' ,'OpenCL ported ' ,  &  !< Calculation of charge density (sumrho) computation
      'Rho_commun    ','Communications' ,'AllReduce grid' ,  &  !< Calculation of charge density (sumrho) communication
      'Pot_commun    ','Communications' ,'AllGathrv grid' ,  &  !< Communication of potential
      'Pot_comm start','Communications' ,'MPI_types/_get' ,  &  !< Communication of potential
      'Un-TransSwitch','Other         ' ,'RMA pattern   ' ,  &  !< Transposition of wavefunction, computation
      'Un-TransComm  ','Communications' ,'ALLtoALLV     ' ,  &  !< Transposition of wavefunction, communication
      'GramS_comput  ','Linear Algebra' ,'DPOTRF        ' ,  &  !< Gram Schmidt computation        
      'GramS_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Gram Schmidt communication
      'LagrM_comput  ','Linear Algebra' ,'DGEMM         ' ,  &  !< Lagrange Multipliers computation
      'LagrM_commun  ','Communications' ,'ALLReduce orbs' ,  &  !< Lagrange Multipliers communication
      'Diis          ','Other         ' ,'Other         ' ,  &  
      !       'PSolv_comput  ','Potential     ' ,'3D FFT        ' ,  &  
      !       'PSolv_commun  ','Communications' ,'ALLtoALL      ' ,  &  
      !       'PSolvKernel   ','Initialization' ,'Miscellaneous ' ,  &  
!      'Exchangecorr  ','Potential     ' ,'Miscellaneous ' ,  &  
      'Forces        ','Finalization  ' ,'Miscellaneous ' ,  &  
      'Tail          ','Finalization  ' ,'Miscellaneous ' ,  &
      'Loewdin_comput','Linear Algebra' ,'              ' ,  &
      'Loewdin_commun','Communications' ,'ALLReduce orbs' ,  &
      'Chol_commun   ','Communications' ,'              ' ,  &
      'Chol_comput   ','Linear Algebra' ,'ALLReduce orbs' ,  &
      'GS/Chol_comput','Linear Algebra' ,'              ' ,  &
      'GS/Chol_commun','Communications' ,'ALLReduce orbs' ,  &
      'Input_comput  ','Initialization' ,'Miscellaneous ' ,  &
      'Input_commun  ','Communications' ,'ALLtoALL+Reduc' ,  &
      'Davidson      ','Finalization  ' ,'Complete SCF  ' ,  &
      'check_IG      ','Initialization' ,'Linear Scaling' ,  &
      'constrc_locreg','Initialization' ,'Miscellaneous ' ,  &
      'wavefunction  ','Initialization' ,'Miscellaneous ' ,  &
      'create_nlpspd ','Initialization' ,'RMA pattern   ' ,  &
      'p2pOrtho_post ','Communications' ,'irecv / irsend' ,  &
      'p2pOrtho_wait ','Communications' ,'mpi_waitany   ' ,  &
      'lovrlp_comm   ','Communications' ,'mpi_allgatherv' ,  &
      'lovrlp_comp   ','Linear Algebra' ,'many ddots    ' ,  &
      'lovrlp_compr  ','Other         ' ,'cut out zeros ' ,  &
      'lovrlp_uncompr','Other         ' ,'insert zeros  ' ,  &
      'extract_orbs  ','Other         ' ,'copy to sendb ' ,  &
      'lovrlp^-1/2   ','Linear Algebra' ,'exact or appr ' ,  &
      'lovrlp^-1/2old','Linear Algebra' ,'exact or appr ' ,  &
      'lovrlp^-1/2com','Linear Algebra' ,'exact or appr ' ,  &
      'lovrlp^-1/2par','Linear Algebra' ,'exact or appr ' ,  &
      'build_lincomb ','Linear Algebra' ,'many daxpy    ' ,  &
      'convolQuartic ','Convolutions  ' ,'No OpenCL     ' ,  &
      'p2pSumrho_wait','Communications' ,'mpi_test/wait ' ,  &
      'sumrho_TMB    ','Other         ' ,'port to GPU?  ' ,  &
      'TMB_kernel    ','Linear Algebra' ,'dgemm         ' ,  &
      'diagonal_seq  ','Linear Algebra' ,'dsygv         ' ,  &
      'diagonal_par  ','Linear Algebra' ,'pdsygvx       ' ,  &
      'lovrlp^-1     ','Linear Algebra' ,'exact or appr ' ,  &
      'lagmat_orthoco','Linear Algebra' ,'dgemm seq/par ' ,  &
      'optimize_DIIS ','Other         ' ,'Other         ' ,  &
      'optimize_SD   ','Other         ' ,'Other         ' ,  &
      'mix_linear    ','Other         ' ,'Other         ' ,  &
      'mix_DIIS      ','Other         ' ,'Other         ' ,  &
      'ig_matric_comm','Communications' ,'mpi p2p       ' ,  &
      'wf_signals    ','Communications' ,'Socket transf.' ,  &
      'energs_signals','Communications' ,'Socket transf.' ,  &
      'rhov_signals  ','Communications' ,'Socket transf.' ,  &
      'init_locregs  ','Initialization' ,'Miscellaneous ' ,  &
      'init_commSumro','Initialization' ,'Miscellaneous ' ,  &
      'init_commPot  ','Initialization' ,'Miscellaneous ' ,  &
      'init_commOrtho','Initialization' ,'Miscellaneous ' ,  &
      'init_inguess  ','Initialization' ,'Miscellaneous ' ,  &
      'init_matrCompr','Initialization' ,'Miscellaneous ' ,  &
      'init_collcomm ','Initialization' ,'Miscellaneous ' ,  &
      'init_collco_sr','Initialization' ,'Miscellaneous ' ,  &
      'init_orbs_lin ','Initialization' ,'Miscellaneous ' ,  &
      'init_repart   ','Initialization' ,'Miscellaneous ' ,  &
      'initMatmulComp','Initialization' ,'Miscellaneous ' ,  &
      'Pot_after_comm','Other         ' ,'global_to_loca' ,  & 
      !       'Init to Zero  ','Other         ' ,'Memset        ' ,  &
      'calc_kernel   ','Other         ' ,'Miscellaneous ' ,  &
      'commun_kernel ','Communications' ,'mpi_allgatherv' ,  &
      'getlocbasinit ','Other         ' ,'Miscellaneous ' ,  &
      'updatelocreg1 ','Other         ' ,'Miscellaneous ' ,  &
      'linscalinit   ','Other         ' ,'Miscellaneous ' ,  &
      'commbasis4dens','Communications' ,'Miscellaneous ' ,  &
      'buildgrad_mcpy','Other         ' ,'Miscellaneous ' ,  &
      'buildgrad_comm','Communications' ,'Allgatherv    ' ,  &
      'allocommsumrho','Communications' ,'Miscellaneous ' ,  &
      'ovrlptransComp','Other         ' ,'Miscellaneous ' ,  &
      'ovrlptransComm','Communications' ,'mpi_allreduce ' ,  &
      'lincombtrans  ','Other         ' ,'Miscellaneous ' ,  &
      'glsynchham1   ','Communications' ,'load balancing' ,  &
      'glsynchham2   ','Communications' ,'load balancing' ,  &
      'gauss_proj    ','Other         ' ,'Miscellaneous ' ,  &
      'sumrho_allred ','Communications' ,'mpiallred     ' ,  &
      'deallocprec   ','Other         ' ,'Miscellaneous ' ,  &
      'large2small   ','Other         ' ,'Miscellaneous ' ,  &
      'small2large   ','Other         ' ,'Miscellaneous ' ,  &
      'renormCoefCom1','Linear Algebra' ,'Miscellaneous ' ,  &
      'renormCoefCom2','Linear Algebra' ,'Miscellaneous ' ,  &
      'renormCoefComm','Communications' ,'Miscellaneous ' ,  &
      'waitAllgatKern','Other         ' ,'Miscellaneous ' ,  &
      'UnBlockPot    ','Other         ' ,'Overlap comms ' ,  &
      'UnBlockDen    ','Other         ' ,'Overlap comms ' ,  &
      'global_local  ','Initialization' ,'Unknown       ' ,  &
      'wfd_creation  ','Other         ' ,'Miscellaneous ' ,  & 
      'comm_llr      ','Communications' ,'Miscellaneous ' ,  &
      !       'AllocationProf','Other         ' ,'Allocate arrs ' ,  &
      'dirmin_lagmat1','Linear Algebra' ,'grad calc     ' ,  &
      'dirmin_lagmat2','Linear Algebra' ,'allgatherv    ' ,  &
      'dirmin_dgesv  ','Linear Algebra' ,'dgesv/pdgesv  ' ,  &
      'dirmin_sddiis ','Linear Algebra' ,'Miscellaneous ' ,  &
      'dirmin_allgat ','Linear Algebra' ,'allgatherv    ' ,  &
      'dirmin_sdfit  ','Linear Algebra' ,'allgatherv etc' ,  &
      'chebyshev_comp','Linear Algebra' ,'matmul/matadd ' ,  &
      'chebyshev_comm','Communications' ,'allreduce     ' ,  &
      'chebyshev_coef','Other         ' ,'Miscellaneous ' ,  &
      'FOE_auxiliary ','Other         ' ,'Miscellaneous ' ,  &
      'FOE_init      ','Other         ' ,'Miscellaneous ' ,  &
      'compressd_mcpy','Other         ' ,'Miscellaneous ' ,  &
      'compressd_comm','Communications' ,'Allgatherv    ' ,  &
      'foe_aux_mcpy  ','Other         ' ,'Miscellaneous ' ,  &
      'foe_aux_comm  ','Communications' ,'Allgatherv    ' ,  &
      'norm_trans    ','Other         ' ,'Miscellaneous ' ,  &
      'misc          ','Other         ' ,'Miscellaneous ' ,  &
      'sparse_copy   ','Other         ' ,'Miscellaneous ' ,  &
      'constraineddft','Other         ' ,'Miscellaneous ' ,  &
      'transfer_int  ','Other         ' ,'Miscellaneous ' ,  &
      'Reformatting  ','Initialization' ,'Interpolation ' ,  &
      'restart_wvl   ','Initialization' ,'inguess    rst' ,  &
      'restart_rsp   ','Initialization' ,'inguess    rst' ,  &
      'check_sumrho  ','Initialization' ,'unitary check ' ,  &
      'check_pot     ','Initialization' ,'unitary check ' ,  &
      'ApplyLocPot   ','Convolutions  ' ,'OpenCL ported ' ,  &
      'ApplyLocKin   ','Convolutions  ' ,'OpenCL ported ' ,  &
      'kernel_init   ','Other         ' ,'Fragment calc ' ,  &
      'calc_energy   ','Linear Algebra' ,'allred etc    ' ,  &
      'new_pulay_corr','Other         ' ,'Pulay forces  ' ,  &
      'dev_from_unity','Other         ' ,'Miscellaneous ' ,  &
      'ks_residue    ','Linear Algebra' ,'Miscellaneous ' ,  &
      'weightanalysis','Linear Algebra' ,'Fragment calc ' ,  &
      'tmbrestart    ','Initialization' ,'Miscellaneous ' ,  &
      'readtmbfiles  ','Initialization' ,'Miscellaneous ' ,  &
      'readisffiles  ','Initialization' ,'Miscellaneous ' ,  &
      'purify_kernel ','Linear Algebra' ,'dgemm         ' ,  &
      'potential_dims','Other         ' ,'auxiliary     ' ,  &
      'sparse_matmul ','Linear Algebra' ,'self-made     ' ,  &
      'transform_matr','Other         ' ,'small to large' ,  &
      'calctrace_comp','Other         ' ,'Miscellaneous ' ,  &
      'calctrace_comm','Communications' ,'allreduce     ' ,  &
      'determinespars','Other         ' ,'Miscellaneous ' ,  &
      'inittaskgroup ','Other         ' ,'Miscellaneous ' ,  &
      'write_matrices','Other         ' ,'dump to disk  ' ,  &
      'transformspars','Other         ' ,'Miscellaneous ' ,  &
      'matrix_extents','Other         ' ,'Miscellaneous ' ,  &
      'lin_inputguess','Other         ' ,'Miscellaneous ' ,  &
      'ionic_energy  ','Other         ' ,'Miscellaneous ' ,  &
      'calc_bounds   ','Other         ' ,'Miscellaneous ' /),(/3,ncat_bigdft/))
 integer, dimension(ncat_bigdft), private, save :: cat_ids !< id of the categories to be converted


 public :: gaussian_basis
 public :: nullify_local_zone_descriptors,locreg_descriptors
 public :: wavefunctions_descriptors,atoms_data,DFT_PSP_projectors
 public :: grid_dimensions,p2pComms,comms_linear,sparse_matrix,matrices
 public :: coulomb_operator,symmetry_data,atomic_structure,comms_cubic
 public :: nonlocal_psp_descriptors
 public :: default_lzd,find_category,old_wavefunction_null,old_wavefunction_free
 public :: bigdft_init_errors,bigdft_init_timing_categories
 public :: deallocate_orbs,deallocate_locreg_descriptors,nullify_wfd
 public :: deallocate_paw_objects!,deallocate_wfd,
 public :: old_wavefunction_set
 public :: nullify_locreg_descriptors
 public :: deallocate_rho_descriptors
 public :: nullify_paw_objects
 public :: cprj_to_array,deallocate_gwf_c
 public :: local_zone_descriptors_null
 public :: energy_terms_null, work_mpiaccumulate_null
 public :: allocate_work_mpiaccumulate, deallocate_work_mpiaccumulate
 public :: nullify_orbitals_data
 public :: SIC_data,orthon_data,input_variables
!!$ public :: SIC_data_null,input_psi_names,input_set,output_wf_format_help
!!$ public :: material_acceleration_null,wf_format_names,input_psi_help,basis_params_set_dict,frag_from_dict



contains

  !> Nullify all energy terms
  pure function energy_terms_null() result(en)
    implicit none
    type(energy_terms) :: en
    en%eh      =0.0_gp 
    en%exc     =0.0_gp 
    en%evxc    =0.0_gp 
    en%eion    =0.0_gp 
    en%edisp   =0.0_gp 
    en%ekin    =0.0_gp 
    en%epot    =0.0_gp
    en%eproj   =0.0_gp
    en%eexctX  =0.0_gp
    en%eelec   =0.0_gp
    en%ebs     =0.0_gp
    en%eKS     =0.0_gp
    en%trH     =0.0_gp
    en%evsum   =0.0_gp
    en%evsic   =0.0_gp 
    en%excrhoc =0.0_gp 
    en%eTS     =0.0_gp
    en%ePV     =0.0_gp 
    en%energy  =0.0_gp 
    en%e_prev  =0.0_gp 
    en%trH_prev=0.0_gp 
    en%c_obj   =int(0,kind=8) 
  end function energy_terms_null


  !> Nullify the data structure associated to Self-Interaction Correction (SIC)
  function old_wavefunction_null() result(wfn)
    implicit none
    type(old_wavefunction) :: wfn
    wfn%Lzd=default_lzd()
    nullify(wfn%psi)
    nullify(wfn%rxyz)
  end function old_wavefunction_null


  function default_lzd() result(lzd)
    type(local_zone_descriptors) :: lzd
    lzd%linear=.false.
    lzd%nlr=0
    lzd%lintyp=0
    lzd%ndimpotisf=0
    lzd%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    lzd%Glr=locreg_null()
    nullify(lzd%Llr)
  end function default_lzd
 
  !> Fills the old_wavefunction structure with corresponding data
  !! Deallocate previous workspaces if already existing
  subroutine old_wavefunction_set(wfn,nat,norbp,Lzd,rxyz,psi)
    
    implicit none
    integer, intent(in) :: nat,norbp
    type(local_zone_descriptors), intent(in) :: Lzd
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp), intent(in) :: psi
    type(old_wavefunction), intent(inout) :: wfn
    !local variables
    character(len=*), parameter :: subname='old_wavefunction_set'

    !first, free the workspace if not already done
    call old_wavefunction_free(wfn)
    !then allocate the workspaces and fill them
    wfn%psi = f_malloc_ptr((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp,id='wfn%psi')
    
    if (norbp>0) call vcopy((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norbp,&
         psi(1),1,wfn%psi(1),1)

    wfn%rxyz = f_malloc_ptr((/ 3, nat /),id='wfn%rxyz')
    if (nat>0) call vcopy(3*nat,rxyz(1,1),1,wfn%rxyz(1,1),1)
    call copy_local_zone_descriptors(Lzd,wfn%Lzd,subname)

  end subroutine old_wavefunction_set


  subroutine old_wavefunction_free(wfn)
    implicit none
    type(old_wavefunction), intent(inout) :: wfn
    !local variables

    if (associated(wfn%psi)) then
       call f_free_ptr(wfn%psi)
    end if
    if (associated(wfn%rxyz)) then
       call f_free_ptr(wfn%rxyz)
    end if
    !lzd should be deallocated also (to be checked again)
    call deallocate_local_zone_descriptors(wfn%Lzd)

  end subroutine old_wavefunction_free
   

  !> De-Allocate gaussian_basis type
  subroutine deallocate_gwf_c(G)
    use module_base
    implicit none
    type(gaussian_basis_c) :: G

    !normally positions should be deallocated outside
    call f_free_ptr(G%ndoc)
    call f_free_ptr(G%nam)
    call f_free_ptr(G%nshell)
    call f_free_ptr(G%psiat)
    call f_free_ptr(G%expof)
    call f_free_ptr(G%rxyz)

  END SUBROUTINE 


!!$  subroutine deallocate_abscalc_input(in)
!!$    use module_base
!!$    implicit none
!!$    type(input_variables) :: in
!!$
!!$    call f_free_ptr(in%Gabs_coeffs)
!!$
!!$  END SUBROUTINE deallocate_abscalc_input


  !> De-Allocate orbitals data structure, except eval pointer
  !! which is not allocated in the orbitals_descriptor routine
  subroutine deallocate_orbs(orbs)
    use module_base
    implicit none
    !Arguments
    type(orbitals_data), intent(inout) :: orbs !< Orbital to de-allocate

    call f_free_ptr(orbs%norb_par)
    call f_free_ptr(orbs%norbu_par)
    call f_free_ptr(orbs%norbd_par)

    call f_free_ptr(orbs%occup)
    call f_free_ptr(orbs%spinsgn)
    call f_free_ptr(orbs%kpts)
    call f_free_ptr(orbs%kwgts)

    call f_free_ptr(orbs%iokpt)

    call f_free_ptr(orbs%ikptproc)

    call f_free_ptr(orbs%inwhichlocreg)

    call f_free_ptr(orbs%onwhichatom)

    call f_free_ptr(orbs%isorb_par)
    if (associated(orbs%ispot)) then
       call f_free_ptr(orbs%ispot)
    end if

  END SUBROUTINE deallocate_orbs

  !> Deallocate rho descriptors
  subroutine deallocate_rho_descriptors(rhodsc)
    use module_base
    implicit none
    type(rho_descriptors) :: rhodsc

    if (associated(rhodsc%spkey))then
       call f_free_ptr(rhodsc%spkey)
    end if
    if (associated(rhodsc%dpkey))then
       call f_free_ptr(rhodsc%dpkey)
    end if
    if (associated(rhodsc%cseg_b))then
       call f_free_ptr(rhodsc%cseg_b)
    end if
    if (associated(rhodsc%fseg_b))then
       call f_free_ptr(rhodsc%fseg_b)
    end if

  end subroutine deallocate_rho_descriptors


  subroutine deallocate_Lzd(Lzd)
    use module_base
    !Arguments
    type(local_zone_descriptors) :: Lzd
    !Local variables
    integer :: ilr

!   nullify the bounds of Glr
    if ((Lzd%Glr%geocode == 'P' .and. Lzd%Glr%hybrid_on) .or. Lzd%Glr%geocode == 'F') then
       nullify(Lzd%Glr%bounds%kb%ibyz_f)
       nullify(Lzd%Glr%bounds%kb%ibxz_f)
       nullify(Lzd%Glr%bounds%kb%ibxy_f)
       nullify(Lzd%Glr%bounds%sb%ibxy_ff)
       nullify(Lzd%Glr%bounds%sb%ibzzx_f)
       nullify(Lzd%Glr%bounds%sb%ibyyzz_f)
       nullify(Lzd%Glr%bounds%gb%ibyz_ff)
       nullify(Lzd%Glr%bounds%gb%ibzxx_f)
       nullify(Lzd%Glr%bounds%gb%ibxxyy_f)
    end if
    !the arrays which are needed only for free BC
    if (Lzd%Glr%geocode == 'F') then
       nullify(Lzd%Glr%bounds%kb%ibyz_c)
       nullify(Lzd%Glr%bounds%kb%ibxz_c)
       nullify(Lzd%Glr%bounds%kb%ibxy_c)
       nullify(Lzd%Glr%bounds%sb%ibzzx_c)
       nullify(Lzd%Glr%bounds%sb%ibyyzz_c)
       nullify(Lzd%Glr%bounds%gb%ibzxx_c)
       nullify(Lzd%Glr%bounds%gb%ibxxyy_c)
       nullify(Lzd%Glr%bounds%ibyyzz_r)
    end if

! nullify the wfd of Glr
   nullify(Lzd%Glr%wfd%keyglob)
   nullify(Lzd%Glr%wfd%keygloc)
!   nullify(Lzd%Glr%wfd%keyv)
   nullify(Lzd%Glr%wfd%keyvloc)
   nullify(Lzd%Glr%wfd%keyvglob)

! nullify the Gnlpspd
!   call deallocate_proj_descr(Lzd%Gnlpspd,subname)
!!$   nullify(Lzd%Gnlpspd%nvctr_p)
!!$   nullify(Lzd%Gnlpspd%nseg_p)
!!$   nullify(Lzd%Gnlpspd%keyv_p)
!!$   nullify(Lzd%Gnlpspd%keyg_p)
!!$   nullify(Lzd%Gnlpspd%nboxp_c)
!!$   nullify(Lzd%Gnlpspd%nboxp_f)
 
!Now destroy the Llr
    do ilr = 1, Lzd%nlr
       call deallocate_locreg_descriptors(Lzd%Llr(ilr))
       !call deallocate_lr(Lzd%Llr(ilr))
!       call deallocate_Lnlpspd(Lzd%Lnlpspd(ilr),subname)
    end do
     nullify(Lzd%Llr)
!     nullify(Lzd%Lnlpspd)

  END SUBROUTINE deallocate_Lzd

  !> Nullify a DFT_local_fields structure
  subroutine nullify_DFT_local_fields(denspot)
    implicit none
    type(DFT_local_fields),intent(out) :: denspot

    nullify(denspot%rhov)
    nullify(denspot%mix)
    nullify(denspot%rho_psi)
    nullify(denspot%rho_C)
    nullify(denspot%V_ext)
    nullify(denspot%V_XC)
    nullify(denspot%Vloc_KS)
    nullify(denspot%f_XC)
    nullify(denspot%rho_work)
    nullify(denspot%pot_work)
    call nullify_rho_descriptors(denspot%rhod)
    !call nullify_denspot_distribution(denspot%dpbox)
    denspot%dpbox=dpbox_null()
    call nullify_coulomb_operator(denspot%pkernel)
    call nullify_coulomb_operator(denspot%pkernelseq)
    
  end subroutine nullify_DFT_local_fields


  subroutine nullify_coulomb_operator(coul_op)
    implicit none
    type(coulomb_operator),intent(out) :: coul_op
    nullify(coul_op%kernel)
  end subroutine nullify_coulomb_operator


  subroutine copy_coulomb_operator(coul1,coul2)
    implicit none
    type(coulomb_operator),intent(in) :: coul1
    type(coulomb_operator),intent(inout) :: coul2

    if(associated(coul2%kernel)) then
      call f_free_ptr(coul2%kernel)
    end if
    coul2%kernel   =>coul1%kernel
    coul2%itype_scf= coul1%itype_scf
    coul2%mu       = coul1%mu
    coul2%geocode  = coul1%geocode
    coul2%ndims    = coul1%ndims
    coul2%hgrids   = coul1%hgrids
    coul2%angrad   = coul1%angrad
    coul2%work1_GPU= coul1%work1_GPU
    coul2%work2_GPU= coul1%work2_GPU
    coul2%k_GPU    = coul1%k_GPU
    coul2%plan     = coul1%plan
    coul2%geo      = coul1%geo
    coul2%igpu     = coul1%igpu
    coul2%initCufftPlan=coul1%initCufftPlan
    coul2%keepGPUmemory=coul1%keepGPUmemory
  ! mpi_env:
    coul2%mpi_env%mpi_comm = coul1%mpi_env%mpi_comm
    coul2%mpi_env%iproc    = coul1%mpi_env%iproc
    coul2%mpi_env%nproc    = coul1%mpi_env%nproc
    coul2%mpi_env%igroup   = coul1%mpi_env%igroup
    coul2%mpi_env%ngroup   = coul1%mpi_env%ngroup

  end subroutine copy_coulomb_operator


  subroutine deallocate_coulomb_operator(coul_op)
    implicit none
    type(coulomb_operator),intent(inout) :: coul_op

    if(associated(coul_op%kernel)) then
      call f_free_ptr(coul_op%kernel)
    end if
    call nullify_coulomb_operator(coul_op)
  end subroutine deallocate_coulomb_operator


  subroutine nullify_rho_descriptors(rhod)
    implicit none
    type(rho_descriptors),intent(out) :: rhod

    nullify(rhod%spkey)
    nullify(rhod%dpkey)
    nullify(rhod%cseg_b)
    nullify(rhod%fseg_b)
  end subroutine nullify_rho_descriptors


  subroutine nullify_GPU_pointers(gpup)
    implicit none
    type(GPU_pointers), intent(out) :: gpup

    nullify(gpup%psi)
    nullify(gpup%ekinpot_host)
    nullify(gpup%psicf_host)
    nullify(gpup%hpsicf_host)
    nullify(gpup%bprecond_host)
    nullify(gpup%ekin)
    nullify(gpup%epot)

  end subroutine nullify_GPU_pointers


  !subroutine nullify_wfn_metadata(wfnmd)
  !  implicit none
  !  type(wfn_metadata),intent(inout) :: wfnmd
  !
  !  nullify(wfnmd%coeff)
  !  nullify(wfnmd%coeff_proj)
  !  nullify(wfnmd%alpha_coeff)
  !  nullify(wfnmd%grad_coeff_old)
  !
  !end subroutine nullify_wfn_metadata


  subroutine nullify_diis_objects(diis)
    implicit none
    type(diis_objects),intent(inout) :: diis

    nullify(diis%psidst)
    nullify(diis%hpsidst)
    nullify(diis%ads)

  end subroutine nullify_diis_objects

  pure subroutine nullify_paw_objects(paw)
    implicit none
    type(paw_objects),intent(inout) :: paw
    
    paw%usepaw = .false.
    nullify(paw%spsi)

    nullify(paw%paw_an)
    nullify(paw%paw_ij)
    nullify(paw%cprj)
    nullify(paw%pawfgrtab)
    nullify(paw%pawrhoij)
  end subroutine nullify_paw_objects

  subroutine deallocate_paw_objects(paw)
    use m_paw_an, only: paw_an_destroy
    use m_paw_ij, only: paw_ij_destroy
    use m_pawcprj, only: pawcprj_destroy
    use m_pawfgrtab, only: pawfgrtab_destroy
    use m_pawrhoij, only: pawrhoij_destroy
    implicit none
    type(paw_objects),intent(inout) :: paw
    
    call f_free_ptr(paw%spsi)

    if (associated(paw%paw_an)) then
       call paw_an_destroy(paw%paw_an)
       deallocate(paw%paw_an)
    end if
    if (associated(paw%paw_ij)) then
       call paw_ij_destroy(paw%paw_ij)
       deallocate(paw%paw_ij)
    end if
    if (associated(paw%cprj)) then
       call pawcprj_destroy(paw%cprj)
       deallocate(paw%cprj)
    end if
    if (associated(paw%pawfgrtab)) then
       call pawfgrtab_destroy(paw%pawfgrtab)
       deallocate(paw%pawfgrtab)
    end if
    if (associated(paw%pawrhoij)) then
       call pawrhoij_destroy(paw%pawrhoij)
       deallocate(paw%pawrhoij)
    end if
  end subroutine deallocate_paw_objects

  subroutine cprj_to_array(cprj,array,norb,nspinor,shift,option)
    implicit none
    integer,intent(in) :: option,norb,nspinor,shift
    real(kind=8),intent(inout) :: array(:,:)
    type(pawcprj_type),intent(inout) :: cprj(:)
    !
    integer :: ii,jj,ilmn,iorb
    !
    if(option==1) then
      do iorb=1,norb*nspinor
        ii=0
        do ilmn=1,cprj(iorb+shift)%nlmn
          do jj=1,2
            ii=ii+1
            array(ii,iorb)=cprj(iorb+shift)%cp(jj,ilmn)
          end do
        end do
      end do
    elseif(option==2) then
      do iorb=1,norb*nspinor
        ii=0
        do ilmn=1,cprj(iorb+shift)%nlmn
          do jj=1,2
            ii=ii+1
            cprj(iorb+shift)%cp(jj,ilmn)=array(ii,iorb)
          end do
        end do
      end do
    end if
  end subroutine cprj_to_array


  pure function work_mpiaccumulate_null() result(w)
    implicit none
    type(work_mpiaccumulate) :: w
    call nullify_work_mpiaccumulate(w)
  end function work_mpiaccumulate_null


  pure subroutine nullify_work_mpiaccumulate(w)
    implicit none
    type(work_mpiaccumulate),intent(out) :: w
    w%ncount = 0
    w%window = 0
    nullify(w%receivebuf)
    nullify(w%sendbuf)
  end subroutine nullify_work_mpiaccumulate


  subroutine allocate_work_mpiaccumulate(w)
    implicit none
    type(work_mpiaccumulate),intent(out) :: w
    w%receivebuf = f_malloc_ptr(w%ncount,id='w%receivebuf')
    w%sendbuf = f_malloc_ptr(w%ncount,id='w%sendbuf')
  end subroutine allocate_work_mpiaccumulate


  subroutine deallocate_work_mpiaccumulate(w)
    implicit none
    type(work_mpiaccumulate),intent(out) :: w
    call f_free_ptr(w%receivebuf)
    call f_free_ptr(w%sendbuf)
  end subroutine deallocate_work_mpiaccumulate


  !> create a null Lzd. Note: this is the correct way of defining 
  !! association through prure procedures.
  !! A pure subroutine has to be defined to create a null structure.
  !! this is important when using the nullification inside other
  !! nullification routines since the usage of a pure function is forbidden
  !! otherwise the routine cannot be pure
  pure function local_zone_descriptors_null() result(lzd)
    implicit none
    type(local_zone_descriptors) :: lzd
    call nullify_local_zone_descriptors(lzd)
  end function local_zone_descriptors_null


  pure subroutine nullify_local_zone_descriptors(lzd)
    implicit none
    type(local_zone_descriptors), intent(out) :: lzd

    lzd%linear=.false.
    lzd%nlr=0
    lzd%lintyp=0
    lzd%ndimpotisf=0
    lzd%hgrids=0.0_gp
    call nullify_locreg_descriptors(lzd%glr)
    nullify(lzd%llr) 
  end subroutine nullify_local_zone_descriptors


  !> Define the BigDFT errors
  subroutine bigdft_init_errors()
    use dictionaries
    !use module_input_keys, only: input_keys_errors
    implicit none
    external :: bigdft_severe_abort

    call f_err_define('BIGDFT_RUNTIME_ERROR',&
         'An invalid operation has been done during runtime',&
         BIGDFT_RUNTIME_ERROR,&
         err_action='Check the exact unrolling of runtime operations,'//&
         ' likely something has been initialized/finalized twice')

    call f_err_define('BIGDFT_MPI_ERROR',&
         'An error of MPI library occurred',&
         BIGDFT_MPI_ERROR,&
         err_action='Check if the error is related to MPI library or runtime conditions')

      call f_err_define('BIGDFT_LINALG_ERROR',&
         'An error of linear algebra occurred',&
         BIGDFT_LINALG_ERROR,&
         err_action='Check if the matrix is correct at input, also look at the info value')

      call f_err_define('BIGDFT_INPUT_VARIABLES_ERROR',&
         'An error while parsing the input variables occured',&
         BIGDFT_INPUT_VARIABLES_ERROR,&
         err_action='Check above which input variable has been not correctly parsed, or check their values')

!!$    !define the errors of internal modules
!!$    call input_keys_errors()

      call f_err_define('BIGDFT_INPUT_FILE_ERROR',&
      'The input file does not exist',&
         BIGDFT_INPUT_FILE_ERROR,&
         err_action='Check if the file does exist')

    !define the severe operation via MPI_ABORT
    call f_err_severe_override(bigdft_severe_abort)
  end subroutine bigdft_init_errors


  !> initialize the timing categories for BigDFT runs.
  !! It is of course assumed that f_lib_initialize has already been called
  subroutine bigdft_init_timing_categories()
    use Poisson_Solver, only: PS_initialize_timing_categories
    implicit none
    !local variables
    integer :: icls,icat

    !initialize categories for the Poisson Solver
    call PS_initialize_timing_categories()

    !initialize groups
    call f_timing_category_group(tgrp_pot,'Operations for local potential construction (mainly XC)')

    do icls=2,ncls_max
       call f_timing_category_group(trim(clss(icls)),'Empty description for the moment')
    end do

    !define the timing categories for exchange and correlation
    call f_timing_category('Exchange-Correlation',tgrp_pot,&
         'Operations needed to construct local XC potential',&
         TCAT_EXCHANGECORR)


    !! little by little, these categories should be transformed in the 
    !! new scheme dictated by f_timing API in time_profiling module of f_lib.

    !initialize categories
    do icat=1,ncat_bigdft
       call f_timing_category(trim(cats(1,icat)),trim(cats(2,icat)),trim(cats(3,icat)),&
            cat_ids(icat))
    end do

  end subroutine bigdft_init_timing_categories


  !> Routine to convert timing categories from the old scheme to the API of f_lib
  !! as soon as the timing categories are identified with their ID, this routine should disappear
  subroutine find_category(category,cat_id)
    use yaml_output, only: yaml_warning
    implicit none
    character(len=*), intent(in) :: category
    integer, intent(out) :: cat_id !< id of the found category
    !local variables
    integer :: i
    !controls if the category exists
    cat_id=0
    do i=1,ncat_bigdft
       if (trim(category) == trim(cats(1,i))) then
          cat_id=cat_ids(i)
          exit
       endif
    enddo
    if (cat_id==0) then
       call f_err_throw('Timing routine error,'//&
            ' requested category '//trim(category)//' has not been found',&
            err_id=TIMING_INVALID)
  !!$     if (bigdft_mpi%iproc==0) &
  !!$          call yaml_warning('Requested timing category '//trim(category)//&
  !!$          ' has not been found')
       cat_id=TIMING_UNINITIALIZED
    end if
  end subroutine find_category

  subroutine nullify_orbitals_data(orbs)
    implicit none

    ! Calling arguments
    type(orbitals_data),intent(out):: orbs

    nullify(orbs%norb_par)
    nullify(orbs%norbu_par)
    nullify(orbs%norbd_par)
    nullify(orbs%iokpt)
    nullify(orbs%ikptproc)
    nullify(orbs%inwhichlocreg)
    nullify(orbs%onwhichatom)
    nullify(orbs%isorb_par)
    nullify(orbs%eval)
    nullify(orbs%occup)
    nullify(orbs%spinsgn)
    nullify(orbs%kwgts)
    nullify(orbs%kpts)
    nullify(orbs%ispot)
    orbs%npsidim_orbs=1
    orbs%npsidim_comp=1

  end subroutine nullify_orbitals_data

  !!integer function matrixindex_in_compressed(this, iorb, jorb)
  !!  implicit none
  !!
  !!  ! Calling arguments
  !!  class(sparse_matrix),intent(in) :: this
  !!  integer,intent(in) :: iorb, jorb
  !!
  !!  ! Local variables
  !!  integer :: compressed_index
  !!
  !!  if (this%store_index) then
  !!      ! Take the value from the array
  !!      matrixindex_in_compressed = this%matrixindex_in_compressed_arr(iorb,jorb)
  !!  else
  !!      ! Recalculate the value
  !!      matrixindex_in_compressed = compressed_index_fn(iorb, jorb, this%full_dim1, this)
  !!  end if
  !!
  !!  contains
  !!    ! Function that gives the index of the matrix element (jjorb,iiorb) in the compressed format.
  !!    integer function compressed_index_fn(irow, jcol, norb, sparsemat)
  !!      !use module_base
  !!      !use module_types
  !!      implicit none
  !!    
  !!      ! Calling arguments
  !!      integer,intent(in) :: irow, jcol, norb
  !!      type(sparse_matrix),intent(in) :: sparsemat
  !!    
  !!      ! Local variables
  !!      integer :: ii, iseg
  !!    
  !!      ii=(jcol-1)*norb+irow
  !!    
  !!      iseg=sparsemat%istsegline(jcol)
  !!      do
  !!          if (ii>=sparsemat%keyg(1,iseg) .and. ii<=sparsemat%keyg(2,iseg)) then
  !!              ! The matrix element is in this segment
  !!               compressed_index_fn = sparsemat%keyv(iseg) + ii - sparsemat%keyg(1,iseg)
  !!              return
  !!          end if
  !!          iseg=iseg+1
  !!          if (iseg>sparsemat%nseg) exit
  !!          if (ii<sparsemat%keyg(1,iseg)) then
  !!              compressed_index_fn=0
  !!              return
  !!          end if
  !!      end do
  !!    
  !!      ! Not found
  !!      compressed_index_fn=0
  !!    
  !!    end function compressed_index_fn
  !!
  !!end function matrixindex_in_compressed


end module module_types

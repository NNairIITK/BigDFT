!> @file
!!  Module to store all dictionary keys of the input files.
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define all static strings to store input variables
module module_input_keys
  use dictionaries
  use module_defs, only: gp
  use f_refcnts
  use f_enums
  use public_keys
  use public_enums
  use multipole_base, only: external_potential_descriptors
  use dynamic_memory
  use fragment_base
  implicit none

  private

  !public :: input_keys_init, input_keys_finalize


  type(dictionary), pointer :: parameters=>null()
  type(dictionary), pointer :: parsed_parameters=>null()
  type(dictionary), pointer :: profiles=>null()


 character(len = 12), dimension(0:2), parameter, public :: OUTPUT_DENSPOT_names = &
       (/ "none        ", &
       "density     ", &
       "dens. + pot." /)

  character(len = 4), dimension(0:2), parameter, public :: OUTPUT_DENSPOT_format_names = &
       (/ "text", &
       "ETSF", &
       "cube" /)

  character(len = 11), dimension(5), parameter, public :: smearing_names = &
       (/ "Error func.", &
       "Fermi      ", &
       "Cold (bumb)", &
       "Cold (mono)", &
       "Meth.-Pax. " /) !< Name of the smearing methods 

  !> Type used for the orthogonalisation parameters
  type, public :: orthon_data
     !> directDiag decides which input guess is chosen:
     !!   if .true. -> as usual direct diagonalization of the Hamiltonian with dsyev (suitable for small systems)
     !!   if .false. -> iterative diagonalization (suitable for large systems)
     logical :: directDiag
     !> norbpInguess indicates how many orbitals shall be treated by each process during the input guess
     !! if directDiag=.false.
     integer :: norbpInguess
     !> You have to choose two numbers for the block size, bsLow and bsUp:
     !!   if bsLow<bsUp, then the program will choose an appropriate block size in between these two numbers
     !!   if bsLow==bsUp, then the program will take exactly this blocksize
     integer :: bsLow
     !> Block size up value (see bsLow)
     integer :: bsUp
     !> The variable methOrtho indicates which orthonormalization procedure is used:
     !!   methOrtho==0 -> Gram-Schmidt with Cholesky decomposition
     !!   methOrtho==1 -> combined block wise classical Gram-Schmidt and Cholesky
     !!   methOrtho==2 -> Loewdin
     integer :: methOrtho
     real(gp) :: iguessTol            !< Gives the tolerance to which the input guess will converged (maximal residue of all orbitals).
     integer :: blocksize_pdsyev      !< Size of the block for the Scalapack routine pdsyev (computes eigenval and vectors)
     integer :: blocksize_pdgemm      !< Size of the block for the Scalapack routine pdgemm
     integer :: nproc_pdsyev          !< Number of proc for the Scalapack routine pdsyev (linear version)
  end type orthon_data

  !> Self-Interaction Correction data
  type, public :: SIC_data
     character(len=4) :: approach !< Approach for the Self-Interaction-Correction (PZ, NK)
     integer :: ixc               !< Base for the SIC correction
     real(gp) :: alpha            !< Downscaling coefficient
     real(gp) :: fref             !< Reference occupation (for alphaNK case)
  end type SIC_data

  !> Contains all parameters related to the linear scaling version.
  type, public :: linearInputParameters 
     integer :: DIIS_hist_lowaccur
     integer :: DIIS_hist_highaccur
     integer :: nItPrecond
     integer :: nItSCCWhenOptimizing
     integer :: nItBasis_lowaccuracy
     integer :: nItBasis_highaccuracy
     integer :: mixHist_lowaccuracy
     integer :: mixHist_highaccuracy
     integer :: dmin_hist_lowaccuracy, dmin_hist_highaccuracy
     integer :: blocksize_pdgemm, blocksize_pdsyev
     integer :: correctionOrthoconstraint, nproc_pdsyev, nproc_pdgemm
     integer :: nit_lowaccuracy, nit_highaccuracy, nItdmin_lowaccuracy, nItdmin_highaccuracy
     integer :: nItSCCWhenFixed_lowaccuracy, nItSCCWhenFixed_highaccuracy, nit_extendedIG
     real(kind=8) :: convCrit_lowaccuracy, convCrit_highaccuracy, convCrit_extendedIG
     real(kind=8) :: alphaSD, alphaDIIS, evlow, evhigh, ef_interpol_chargediff
     real(kind=8) :: alpha_mix_lowaccuracy, alpha_mix_highaccuracy, reduce_confinement_factor, ef_interpol_det
     integer :: plotBasisFunctions
     real(kind=8) :: fscale, deltaenergy_multiplier_TMBexit, deltaenergy_multiplier_TMBfix
     real(kind=8) :: lowaccuracy_conv_crit, convCritMix_lowaccuracy, convCritMix_highaccuracy
     real(kind=8) :: highaccuracy_conv_crit, support_functions_converged, alphaSD_coeff
     real(kind=8) :: convCritDmin_lowaccuracy, convCritDmin_highaccuracy
     real(kind=8), dimension(:), pointer :: locrad, locrad_lowaccuracy, locrad_highaccuracy, kernel_cutoff_FOE
     real(kind=8), dimension(:,:), pointer :: locrad_type
     real(kind=8), dimension(:), pointer :: potentialPrefac_lowaccuracy, potentialPrefac_highaccuracy, potentialPrefac_ao
     real(kind=8), dimension(:),pointer :: kernel_cutoff, locrad_kernel, locrad_mult
     real(kind=8) :: early_stop, gnrm_dynamic, min_gnrm_for_dynamic 
     integer, dimension(:), pointer :: norbsPerType
     integer :: kernel_mode, mixing_mode
     integer :: scf_mode, nlevel_accuracy
     logical :: calc_dipole, pulay_correction, iterative_orthogonalization, new_pulay_correction
     logical :: fragment_calculation, calc_transfer_integrals, constrained_dft, curvefit_dmin, diag_end, diag_start
     integer :: extra_states, order_taylor, mixing_after_inputguess
     !> linear scaling: maximal error of the Taylor approximations to calculate the inverse of the overlap matrix
     real(kind=8) :: max_inversion_error
    logical :: calculate_onsite_overlap
     integer :: output_mat_format     !< Output Matrices format
     integer :: output_coeff_format   !< Output Coefficients format
     integer :: output_fragments   !< Output fragments/full system/both
     logical :: charge_multipoles !< Calculate the multipoles expansion coefficients of the charge density
     integer :: kernel_restart_mode !< How to generate the kernel in a restart calculation
  end type linearInputParameters

  !> Structure controlling the nature of the accelerations (Convolutions, Poisson Solver)
  type, public :: material_acceleration
     !> variable for material acceleration
     !! values 0: traditional CPU calculation
     !!        1: CUDA acceleration with CUBLAS
     !!        2: OpenCL acceleration (with CUBLAS one day)
     !!        3: OpenCL accelearation for CPU
     !!        4: OCLACC (OpenCL both ? see input_variables.f90)
     integer :: iacceleration
     integer :: Psolver_igpu            !< Acceleration of the Poisson solver
     character(len=11) :: OCL_platform  !< Name of the OpenCL platform
     character(len=11) :: OCL_devices   !< Name of the OpenCL devices
  end type material_acceleration


  !> Structure of the variables read by input.* files (*.dft, *.geopt...)
  type, public :: input_variables

     !>reference counter
     type(f_reference_counter) :: refcnt
     !> enumerator for the mode of the run
     type(f_enumerator) :: run_mode
     !> Strings of the input files
     character(len=100) :: file_occnum !< Occupation number (input)
     character(len=100) :: file_igpop
     character(len=100) :: file_lin   
     character(len=100) :: file_frag   !< Fragments
     character(len=max_field_length) :: dir_output  !< Strings of the directory which contains all data output files
     !integer :: files                  !< Existing files.

     !> Miscellaneous variables
!     logical :: gaussian_help !to be used as iterator attribute
     integer :: itrpmax
     integer :: iscf
     integer :: norbsempty
     integer :: norbsuempty,norbsdempty
     integer :: occopt
     integer :: last_run
     real(gp) :: frac_fluct,gnrm_sw,alphamix
     real(gp) :: Tel                 !< Electronic temperature for the mixing scheme
     real(gp) :: alphadiis
     real(gp) :: rpnrm_cv
     real(gp) :: gnrm_startmix
     integer :: verbosity            !< Verbosity of the output file
     logical :: multipole_preserving !< Preserve multipole for ionic charge (integrated isf)
     integer :: mp_isf               !< Interpolating scaling function order for multipole preserving

     !> DFT basic parameters.
     integer :: ixc         !< XC functional Id
     real(gp):: qcharge     !< Total charge of the system
     integer :: itermax     !< Maximal number of SCF iterations
     integer :: itermin     !< Minimum number of SCF iterations !Bastian
     integer :: nrepmax
     integer :: ncong       !< Number of conjugate gradient iterations for the preconditioner
     integer :: idsx        !< DIIS history
     integer :: ncongt      !< Number of conjugate garident for the tail treatment
     type(f_enumerator) :: inputpsiid  !< Input PSI choice
     !!   - 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
     !!   - 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
     !!         as definition of the previous system.
     !!   - 2 : read waves from disk
     integer :: nspin       !< Spin components (no spin 1, collinear 2, non collinear 4)
     integer :: mpol        !< Total spin polarisation of the system
     integer :: norbv       !< Number of virtual orbitals to compute after direct minimisation
     !! using the Davidson scheme
     integer :: nvirt
     integer :: nplot
     type(f_enumerator) :: output_denspot        !< 0= No output, 1= density, 2= density+potential
     integer :: dispersion            !< Dispersion term
     type(f_enumerator) :: output_wf!_format      !< Output Wavefunction format
     !integer :: output_denspot_format !< Format for the output density and potential
     real(gp) :: hx,hy,hz   !< Step grid parameter (hgrid)
     real(gp) :: crmult     !< Coarse radius multiplier
     real(gp) :: frmult     !< Fine radius multiplier
     real(gp) :: gnrm_cv    !< Convergence parameters of orbitals
     real(gp) :: rbuf       !< buffer for tail treatment
     real(gp), dimension(3) :: elecfield   !< Electric Field vector
     logical :: disableSym                 !< .true. disable symmetry
     !> boolean to activate the calculation of the stress tensor
     logical :: calculate_strten
     character(len=8) :: set_epsilon !< method for setting the dielectric constant

     !> For absorption calculations
     integer :: iabscalc_type   !< 0 non calc, 1 cheb ,  2 lanc
     !! integer :: iat_absorber, L_absorber, N_absorber, rpower_absorber, Linit_absorber
     integer :: iat_absorber,  L_absorber
     real(gp), dimension(:), pointer :: Gabs_coeffs
     real(gp) :: abscalc_bottomshift
     logical ::  c_absorbtion
     logical ::  abscalc_alterpot
     logical ::  abscalc_eqdiff
     logical ::  abscalc_S_do_cg
     logical ::  abscalc_Sinv_do_cg
     integer ::  potshortcut
     integer ::  nsteps
     character(len=100) :: extraOrbital
     character(len=1000) :: xabs_res_prefix

     !> Frequencies calculations (finite difference)
     real(gp) :: freq_alpha  !< Factor for the finite difference step (step = alpha * hgrid)
     integer :: freq_order   !< Order of the finite difference scheme
     integer :: freq_method  !< Method to calculate the frequencies

     !> kpoints related input variables
     integer :: gen_nkpt                            !< K points to generate the results
     real(gp), dimension(:,:), pointer :: gen_kpt   !< K points coordinates
     real(gp), dimension(:), pointer :: gen_wkpt    !< Weights of k points
     !> Band structure path
     integer :: nkptv
     integer :: ngroups_kptv
     integer, dimension(:), pointer :: nkptsv_group
     real(gp), dimension(:,:), pointer :: kptv
     character(len=100) :: band_structure_filename

     ! Orbitals and input occupation.
     integer :: gen_norb, gen_norbu, gen_norbd
     real(gp), dimension(:), pointer :: gen_occup

     ! Geometry variables from *.geopt
     character(len=10) :: geopt_approach !<id of geopt driver
     integer :: ncount_cluster_x         !< Maximum number of geopt steps 
     integer :: wfn_history              !< Number of previous steps saved for wfn reformatting
     integer :: history                  !< History of DIIS method
     real(gp) :: betax,forcemax,randdis
     integer :: optcell, ionmov, nnos
     real(gp) :: dtion
     real(gp) :: mditemp, mdftemp
     real(gp) :: noseinert, friction, mdwall
     real(gp) :: bmass, vmass, strprecon, strfact
     integer:: sockinet, sockport
     character(len=1032)  :: sockhost
     real(gp), dimension(6) :: strtarget
     real(gp), dimension(:), pointer :: qmass
     real(gp) :: dtinit, dtmax           !< For FIRE
     character(len=10) :: tddft_approach !< TD-DFT variables from *.tddft
     type(SIC_data) :: SIC               !< Parameters for the SIC methods
     !variables for SQNM
     integer  :: nhistx
     logical  :: biomode
     real(gp) :: beta_stretchx
     real(gp) :: maxrise
     real(gp) :: cutoffratio
     real(gp) :: steepthresh
     real(gp) :: trustr

     !Force Field Parameter
     character(len=64) :: mm_paramset
     character(len=64) :: mm_paramfile

     !MD input keywords
     integer :: mdsteps
     integer :: md_printfrq
     real(gp) :: temperature
     real(gp) :: dt
     logical  :: no_translation
     logical  :: nhc
     integer  :: nhnc
     integer  :: nmultint
     integer  :: nsuzuki
     real(gp) :: nosefrq 

     ! Performance variables from input.perf
     logical :: debug      !< Debug option (used by memocc)
     integer :: ncache_fft !< Cache size for FFT
     real(gp) :: projrad   !< Coarse radius of the projectors in units of the maxrad
     real(gp) :: symTol    !< Tolerance for symmetry detection.
     integer :: linear
     logical :: signaling                    !< Expose results on DBus or Inet.
     integer :: signalTimeout                !< Timeout for inet connection.
     character(len = 64) :: domain           !< Domain to get the IP from hostname.
     double precision :: gmainloop           !< Internal C pointer on the signaling structure.
     integer :: inguess_geopt                !< 0= Wavelet input guess, 1 = real space input guess 

     !> Orthogonalisation data
     type(orthon_data) :: orthpar

     !> Acceleration parameters
     type(material_acceleration) :: matacc

     ! Parallelisation parameters
     !> Parallelisation scheme of the exact exchange operator
     !!   BC (Blocking Collective)
     !!   OP2P (Overlap Point-to-Point)
     character(len=4) :: exctxpar
     !> Paradigm for unblocking global communications via OMP_NESTING
     character(len=3) :: unblock_comms
     !> Communication scheme for the density
     !!   DBL traditional scheme with double precision
     !!   MIX mixed single-double precision scheme (requires rho_descriptors)
     character(len=3) :: rho_commun
     !> Number of taskgroups for the poisson solver
     !! works only if the number of MPI processes is a multiple of it
     integer :: PSolver_groupsize
     !> Global MPI group size (will be written in the mpi_environment)
     ! integer :: mpi_groupsize 

     type(external_potential_descriptors) :: ep

     ! Linear scaling parameters
     type(linearInputParameters) :: lin    !< Linear scaling data
     type(fragmentInputParameters) :: frag !< Fragment data
     logical :: store_index                !< (LS) Store indices of the sparse matrices or recalculate them 
     integer :: check_sumrho               !< (LS) Perform a check of sumrho (no check, light check or full check)
     integer :: check_overlap              !< (LS) Perform a check of the overlap calculation
     logical :: experimental_mode          !< (LS) Activate the experimental mode
     integer :: write_orbitals             !< (LS) Write KS orbitals for cubic restart (0: no, 1: wvl, 2: wvl+isf)
     logical :: explicit_locregcenters     !< (LS) Explicitely specify localization centers
     logical :: calculate_KS_residue       !< (LS) Calculate Kohn-Sham residue
     logical :: intermediate_forces        !< (LS) Calculate intermediate forces

     !> linear scaling: exit kappa for extended input guess (experimental mode)
     real(kind=8) :: kappa_conv

     !> linear scaling: number of FOE cycles before the eigenvalue bounds are shrinked
     integer :: evbounds_nsatur

     !> linear scaling: maximal number of unsuccessful eigenvalue bounds shrinkings
     integer :: evboundsshrink_nsatur

     !> linear scaling: how to update the density kernel during the support function optimization (0: purification, 1: FOE)
     integer :: method_updatekernel

     !> linear scaling: quick return in purification
     logical :: purification_quickreturn

     !> linear scaling: dynamic adjustment of the decay length of the FOE error function
     logical :: adjust_FOE_temperature

     !> linear scaling: calculate the HOMO LUMO gap even when FOE is used for the kernel calculation
     logical :: calculate_gap

     !> linear scaling: perform a Loewdin charge analysis at the end of the calculation
     logical :: loewdin_charge_analysis

     !> linear scaling: perform a Loewdin charge analysis of the coefficients for fragment calculations
     logical :: coeff_weight_analysis

     !> linear scaling: perform a check of the matrix compression routines
     logical :: check_matrix_compression

     !> linear scaling: correction covariant / contravariant gradient
     logical :: correction_co_contra

     !> linear scaling: lower bound for the error function decay length
     real(kind=8) :: fscale_lowerbound

     !> linear scaling: upper bound for the error function decay length
     real(kind=8) :: fscale_upperbound

     !> linear scaling: Restart method to be used for the FOE method
     integer :: FOE_restart

     !> linear scaling: method to calculate the overlap matrices (1=old, 2=new)
     integer :: imethod_overlap

     !> linear scaling: enable the matrix taskgroups
     logical :: enable_matrix_taskgroups

     !> linear scaling: radius enlargement for the Hamiltonian application (in grid points)
     integer :: hamapp_radius_incr

     !> linear scaling: enable the addaptive ajustment of the number of kernel iterations
     logical :: adjust_kernel_iterations

     !> linear scaling: perform an analysis of the extent of the support functions (and possibly KS orbitals)
     logical :: wf_extent_analysis

     !> Method for the solution of  generalized poisson Equation
     character(len=4) :: GPS_Method

     !> Use the FOE method to calculate the HOMO-LUMO gap at the end
     logical :: foe_gap

  end type input_variables

  interface input_set
     module procedure input_set_char, input_set_int, input_set_dbl, input_set_bool, &
          & input_set_int_array, input_set_dbl_array, input_set_bool_array, &
          & input_set_dict
  end interface input_set

  public :: SIC_data_null
  public :: free_input_variables,print_dft_parameters,wave_format_from_filename
  public :: inputpsiid_get_policy,inputpsiid_set_policy,set_inputpsiid
  ! Main creation routine
  public :: user_dict_from_files,inputs_from_dict
  public :: input_keys_dump,input_set,input_keys_fill_all,print_general_parameters

contains

  pure function SIC_data_null() result(SIC)
    implicit none
    type(SIC_data) :: SIC

    SIC%approach=repeat(' ',len(SIC%approach))
    SIC%ixc=0
    SIC%alpha=0.0_gp
    SIC%fref=0.0_gp 
  end function SIC_data_null

  function material_acceleration_null() result(ma)
    type(material_acceleration) :: ma
    ma%iacceleration=0
    ma%Psolver_igpu=0
    ma%OCL_platform=repeat(' ',len(ma%OCL_platform))
    ma%OCL_platform=repeat(' ',len(ma%OCL_devices))
  end function material_acceleration_null

!!$  function input_psi_validate(id)
!!$    integer, intent(in) :: id
!!$    logical :: input_psi_validate
!!$
!!$    integer :: i
!!$
!!$    input_psi_validate = .false.
!!$    do i = 1, size(input_psi_values)
!!$       if (id == input_psi_values(i)) then
!!$          input_psi_validate = .true.
!!$          return
!!$       end if
!!$    end do
!!$  end function input_psi_validate
!!$
!!$  function output_wf_format_validate(id)
!!$    integer, intent(in) :: id
!!$    logical :: output_wf_format_validate
!!$
!!$    output_wf_format_validate = (id >= 0 .and. id < size(wf_format_names))
!!$  end function output_wf_format_validate

  subroutine output_denspot_help()
    integer :: i, j

    write(*, "(1x,A)") "Available values of output_denspot are:"
    do i = 0, size(output_denspot_format_names) - 1
       do j = 0, size(output_denspot_names) - 1
          if (j == 0 .and. i == 0) then
             write(*, "(1x,A,I5,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_denspot_names(j)), "."
          else if (j /= 0) then
             write(*, "(1x,A,I5,A,A,A,A,A)") " | ", i * 10 + j, &
                  & " - ", trim(output_denspot_names(j)), &
                  & " in ", trim(output_denspot_format_names(i)), " format."
          end if
       end do
    end do
  end subroutine output_denspot_help

  function output_denspot_validate(id, fid)
    integer, intent(in) :: id, fid
    logical :: output_denspot_validate

    output_denspot_validate = (id >= 0 .and. id < size(output_denspot_names)) .and. &
         & (fid >= 0 .and. fid < size(output_denspot_format_names))
  end function output_denspot_validate

  !> Nullify the linear Input parameters
  subroutine nullifyInputLinparameters(lin)
    implicit none

    ! Calling arguments
    type(linearInputParameters),intent(inout):: lin

    nullify(lin%locrad)
    nullify(lin%potentialPrefac_lowaccuracy)
    nullify(lin%potentialPrefac_highaccuracy)
    nullify(lin%norbsPerType)
    nullify(lin%potentialPrefac_ao)
    nullify(lin%locrad_type)
    nullify(lin%kernel_cutoff_FOE)
    nullify(lin%kernel_cutoff)

  end subroutine nullifyInputLinparameters

  subroutine input_keys_init()
    use yaml_output
    use dynamic_memory
    use yaml_parse
    use f_precisions, only: f_integer
    implicit none
    !local variables
    integer(f_integer) :: params_size
    !integer(kind = 8) :: cbuf_add !< address of c buffer
    character, dimension(:), allocatable :: params

    call f_routine(id='input_keys_init')

    !alternative filling of parameters from hard-coded source file
    !call getstaticinputdef(cbuf_add,params_size)
    call getinputdefsize(params_size)
    !allocate array
    params=f_malloc_str(1,params_size,id='params')
    !fill it and parse dictionary
    !print *,'after', f_loc(params),f_loc(params(1)),'shape',shape(params),params_size
    !print *,'cbuf_add',cbuf_add
    call getinputdef(params)
    !write(*,*)'here definition'
    !write(*,'('//trim(yaml_toa(params_size))//'a)')params

    !call copycbuffer(params,cbuf_add,params_size)
    !print *,'there',params_size
    call yaml_parse_from_char_array(parsed_parameters,params)
    !there is only one document in the input variables specifications
    parameters=>parsed_parameters//0
    profiles => parsed_parameters//1
    call f_free_str(1,params)

    !call yaml_dict_dump(parameters, comment_key = COMMENT)
    
!!$    !in the case the errors have not been initialized before
!!$    call input_keys_errors()

    call f_release_routine()

  END SUBROUTINE input_keys_init


  subroutine input_keys_finalize()
    use dictionaries
    implicit none

    if (associated(parsed_parameters)) then
       call dict_free(parsed_parameters)
       nullify(parameters)
       nullify(profiles)
    else
       call dict_free(parameters)
    end if
  END SUBROUTINE input_keys_finalize

  !> Fill the input_variables and atoms_data structures from the information
  !! contained in the dictionary dict
  !! the dictionary should be completes to fill all the information
  subroutine inputs_from_dict(in, atoms, dict)
    use module_defs, only: DistProjApply,pi_param
    use module_base, only: bigdft_mpi
    use yaml_output
    use dictionaries
    use module_input_dicts
    use dynamic_memory
    use f_utils, only: f_zero
    use f_input_file, only: input_keys_get_profile
    use module_xc
    !  use input_old_text_format, only: dict_from_frag
    use module_atoms!, only: atoms_data,atoms_data_null,atomic_data_set_from_dict,&
                    ! check_atoms_positions,psp_set_from_dict,astruct_set_from_dict
    use yaml_strings, only: f_strcpy
    use m_ab6_symmetry, only: symmetry_get_n_sym
    use interfaces_42_libpaw
    use multipole_base, only: external_potential_descriptors, multipoles_from_dict, lmax
    use public_enums
    use fragment_base
    use f_utils, only: f_get_free_unit
    use wrapper_MPI, only: mpibarrier
    use yaml_strings, only: yaml_toa
    implicit none
    !Arguments
    type(input_variables), intent(out) :: in
    type(atoms_data), intent(out) :: atoms
    type(dictionary), pointer :: dict
    !Local variables
    !type(dictionary), pointer :: profs, dict_frag
    logical :: found, userdef
    integer :: ierr, norb_max, jtype, jxc
    real(gp) :: qelec_up, qelec_down
    character(len = max_field_length) :: msg,filename,run_id,input_id,posinp_id,outdir
    !  type(f_dict) :: dict
    type(dictionary), pointer :: dict_minimal, var, lvl, types

    integer, parameter :: pawlcutd = 10, pawlmix = 10, pawnphi = 13, pawntheta = 12, pawxcdev = 1
    integer, parameter :: xclevel = 1, usepotzero = 0
    integer :: nsym,unt
    real(gp) :: gsqcut_shp, rloc, projr, rlocmin
    real(gp), dimension(2) :: cfrmults
    type(external_potential_descriptors) :: ep
    integer :: impl, l

    !  dict => dict//key

    !  dict = dict//key

    call f_routine(id='inputs_from_dict')


    ! Atoms case.
    !atoms = atoms_data_null()
    call nullify_atoms_data(atoms)

    if (.not. has_key(dict, POSINP)) &
         call f_err_throw("missing posinp",err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    call astruct_set_from_dict(dict // POSINP, atoms%astruct)
    ! Generate the dict of types for later use.
    call astruct_dict_get_types(dict // POSINP, types)

    ! Input variables case.
    call default_input_variables(in)

    !call yaml_map('Dictionary parsed',dict)

    ! extract also the minimal dictionary which is necessary to do this run
    call input_keys_fill_all(dict,dict_minimal)

    ! Add missing pseudo information.
    projr = dict // PERF_VARIABLES // PROJRAD
    cfrmults = dict // DFT_VARIABLES // RMULT
    jxc = dict // DFT_VARIABLES // IXC
    var => dict_iter(types)
    do while(associated(var))
       call psp_dict_fill_all(dict, trim(dict_key(var)), jxc, projr, cfrmults(1), cfrmults(2))
       var => dict_next(var)
    end do

    ! Update interdependant values.
    rlocmin = 999._gp
    var => dict_iter(types)
    do while(associated(var))
       call psp_set_from_dict(dict // ("psppar." // trim(dict_key(var))), rloc = rloc)
       rlocmin = min(rloc, rlocmin)
       var => dict_next(var)
    end do
    select case (trim(input_keys_get_profile(dict // DFT_VARIABLES, HGRIDS, userdef)))
    case ("accurate")
       call set(dict // DFT_VARIABLES // HGRIDS, (/ rlocmin, rlocmin, rlocmin /) * 1.25_gp)
    case ("normal")
       call set(dict // DFT_VARIABLES // HGRIDS, (/ rlocmin, rlocmin, rlocmin /) * 1.75_gp)
    case ("fast")
       call set(dict // DFT_VARIABLES // HGRIDS, (/ rlocmin, rlocmin, rlocmin /) * 2.5_gp)
    end select

    ! Transfer dict values into input_variables structure.
    lvl => dict_iter(dict)
    do while(associated(lvl))
       var => dict_iter(lvl)
       if (.not. associated(var)) then
          ! Scalar case.
          call input_set(in, trim(dict_key(lvl)), lvl)
       else
          do while(associated(var))
             call input_set(in, trim(dict_key(lvl)), var)
             var => dict_next(var)
          end do
       end if
       lvl => dict_next(lvl)
    end do

    ! Generate the dir_output
    !outdir has to be initialized
    call f_zero(outdir)
    call dict_get_run_properties(dict, naming_id = run_id, posinp_id = posinp_id, input_id = input_id, outdir_id = outdir)
    call f_strcpy(dest = in%dir_output, src = trim(outdir) // "data" // trim(run_id))

    call set_cache_size(in%ncache_fft)

    !status of the allocation verbosity and profiling
    !filename of the memory allocation status, if needed
    call f_strcpy(src=trim(outdir) // 'memstatus' // trim(run_id) // '.yaml',&
         dest=filename)
    if (.not. in%debug) then
       if (in%verbosity==3) then
          call f_malloc_set_status(output_level=1, iproc=bigdft_mpi%iproc,logfile_name=filename)
       else
          call f_malloc_set_status(output_level=0, iproc=bigdft_mpi%iproc)
       end if
    else
       call f_malloc_set_status(output_level=2, iproc=bigdft_mpi%iproc,logfile_name=filename)
    end if

    call nullifyInputLinparameters(in%lin)
    call allocateBasicArraysInputLin(in%lin, atoms%astruct%ntypes)

    !First fill all the types by the default, then override by per-type values
    lvl => dict_iter(types)
    do while(associated(lvl))
       jtype = lvl
       var => dict_iter(dict//LIN_BASIS_PARAMS)
       do while(associated(var))
          call basis_params_set_dict(var,in%lin,jtype)
          var => dict_next(var)
       end do
       !then check if the objects exists in separate specifications
       if (has_key(dict//LIN_BASIS_PARAMS,trim(dict_key(lvl)))) then
          var => &
               dict_iter(dict//LIN_BASIS_PARAMS//trim(dict_key(lvl)))
       end if
       do while(associated(var))
          call basis_params_set_dict(var,in%lin,jtype)
          var => dict_next(var)
       end do
       lvl => dict_next(lvl)
    end do

    call allocate_extra_lin_arrays(in%lin,in%nspin,atoms%astruct)

    ! Cross check values of input_variables.
    call input_analyze(in,atoms%astruct)

    ! Shake atoms, if required.
    call astruct_set_displacement(atoms%astruct, in%randdis)
    if (bigdft_mpi%nproc > 1) call mpibarrier(bigdft_mpi%mpi_comm)
    ! Update atoms with symmetry information
    call astruct_set_symmetries(atoms%astruct, in%disableSym, in%symTol, in%elecfield, in%nspin)

    call kpt_input_analyse(bigdft_mpi%iproc, in, dict//KPT_VARIABLES, &
         & atoms%astruct%sym, atoms%astruct%geocode, atoms%astruct%cell_dim)

    ! Update atoms with pseudo information.
    call psp_dict_analyse(dict, atoms)
    call atomic_data_set_from_dict(dict,IG_OCCUPATION, atoms, in%nspin)

    ! Add multipole preserving information
    atoms%multipole_preserving = in%multipole_preserving
    atoms%mp_isf = in%mp_isf

    ! Generate orbital occupation
    call read_n_orbitals(bigdft_mpi%iproc, qelec_up, qelec_down, norb_max, atoms, &
         in%qcharge, in%nspin, in%mpol, in%norbsempty)
    call occupation_set_from_dict(dict, OCCUPATION, &
         in%gen_norbu, in%gen_norbd, in%gen_occup, &
         in%gen_nkpt, in%nspin, in%norbsempty, qelec_up, qelec_down, norb_max)
    in%gen_norb = in%gen_norbu + in%gen_norbd

    ! Complement PAW initialisation.
    if (any(atoms%npspcode == PSPCODE_PAW)) then
       !gsqcut_shp = two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
       gsqcut_shp = 2._gp * 2.2_gp / pi_param ** 2
       call symmetry_get_n_sym(atoms%astruct%sym%symObj, nsym, ierr)
       call pawinit(1, gsqcut_shp, pawlcutd, pawlmix, maxval(atoms%pawtab(:)%lmn_size) + 1, &
            & pawnphi, nsym, pawntheta, atoms%pawang, atoms%pawrad, 0, &
            & atoms%pawtab, pawxcdev, xclevel, usepotzero)
    end if

    if (in%gen_nkpt > 1 .and. (in%inputpsiid .hasattr. 'GAUSSIAN')) then
       call f_err_throw('Gaussian projection is not implemented with k-point support',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

!!$    if(in%inputpsiid == INPUT_PSI_LINEAR_AO .or. &
!!$         in%inputpsiid == INPUT_PSI_MEMORY_LINEAR .or. &
!!$         in%inputpsiid == INPUT_PSI_DISK_LINEAR) &
    if (in%inputpsiid .hasattr. 'LINEAR') DistProjApply=.true.
    if(in%linear /= INPUT_IG_OFF .and. in%linear /= INPUT_IG_LIG) then
       !only on the fly calculation
       DistProjApply=.true.
    end if

    !if other steps are supposed to be done leave the last_run to minus one
    !otherwise put it to one
    if (in%last_run == -1 .and. in%ncount_cluster_x <=1 .or. in%ncount_cluster_x <= 1) then
       in%last_run = 1
    end if

    ! Stop the code if it is trying to run GPU with spin=4
    if (in%nspin == 4 .and. in%matacc%iacceleration /= 0) then
       call f_err_throw('GPU calculation not implemented with non-collinear spin',err_name='BIGDFT_INPUT_VARIABLES_ERROR')
!!$     if (bigdft_mpi%iproc==0) call yaml_warning('GPU calculation not implemented with non-collinear spin')
!!$     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
    end if

    !control atom positions
    call check_atoms_positions(atoms%astruct, (bigdft_mpi%iproc == 0))

!!$  ! Stop code for unproper input variables combination.
!!$  if (in%ncount_cluster_x > 0 .and. .not. in%disableSym .and. atoms%geocode == 'S') then
!!$     if (bigdft_mpi%iproc==0) then
!!$        write(*,'(1x,a)') 'Change "F" into "T" in the last line of "input.dft"'   
!!$        write(*,'(1x,a)') 'Forces are not implemented with symmetry support, disable symmetry please (T)'
!!$     end if
!!$     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
!!$  end if

!!$  if (bigdft_mpi%iproc == 0) then
!!$     profs => input_keys_get_profiles("")
!!$     call yaml_dict_dump(profs)
!!$     call dict_free(profs)
!!$  end if

    !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

    ! Warn for all INPUT_VAR_ILLEGAL errors.
    do while (f_err_pop(err_name='INPUT_VAR_ILLEGAL', add_msg = msg) /= 0)
       call yaml_warning(trim(msg))
    end do
    !check if an error has been found and raise an exception to be handled
    if (f_err_check()) then
       call f_err_throw('Error in reading input variables from dictionary',&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

    ! not sure whether to actually make this an input variable or not so just set to false for now
    in%lin%diag_start=.false.

    !then fill also fragment variables
    in%lin%fragment_calculation=FRAG_VARIABLES .in. dict
    in%lin%calc_transfer_integrals=.false.
    in%lin%constrained_dft=.false.
    if (in%lin%fragment_calculation) then
       in%lin%constrained_dft=CONSTRAINED_DFT .in. dict // FRAG_VARIABLES
       found = TRANSFER_INTEGRALS .in. dict // FRAG_VARIABLES
       if (found) in%lin%calc_transfer_integrals=dict//FRAG_VARIABLES//TRANSFER_INTEGRALS
       call frag_from_dict(dict//FRAG_VARIABLES,in%frag)

!!$     ! again recheck
!!$     call dict_from_frag(in%frag,dict_frag)
!!$     call yaml_map('again',dict_frag)
!!$     call dict_free(dict_frag)
!!$     stop
    else
       call default_fragmentInputParameters(in%frag)
    end if

    ! Process the multipoles for the external potential
    call multipoles_from_dict(dict//DFT_VARIABLES//EXTERNAL_POTENTIAL, in%ep)
    !!do impl=1,in%ep%nmpl
    !!    call yaml_map('rxyz',in%ep%mpl(impl)%rxyz)
    !!    do l=0,lmax
    !!         if(associated(in%ep%mpl(impl)%qlm(l)%q)) then
    !!             call yaml_map(trim(yaml_toa(l)),in%ep%mpl(impl)%qlm(l)%q)
    !!         end if
    !!    end do
    !!end do

    ! No use anymore of the types.
    call dict_free(types)

    if (bigdft_mpi%iproc==0) then
       call input_keys_dump(dict)
    end if

    !check whether a directory name should be associated for the data storage
    call check_for_data_writing_directory(bigdft_mpi%iproc,in)

    if (bigdft_mpi%iproc == 0)  call print_general_parameters(in,atoms,input_id,posinp_id)

    if (associated(dict_minimal) .and. bigdft_mpi%iproc == 0) then
       call dict_get_run_properties(dict, input_id = run_id)
       call f_strcpy(src=trim(run_id)//'_minimal.yaml',dest=filename)
       unt=f_get_free_unit(99971)
       call yaml_set_stream(unit=unt,filename=trim(outdir)//trim(filename),&
            record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
       if (ierr==0) then
          call yaml_comment('Minimal input file',hfill='-',unit=unt)
          call yaml_comment('This file indicates the minimal set of input variables which has to be given '//&
               'to perform the run. The code would produce the same output if this file is used as input.',unit=unt)
          call yaml_dict_dump(dict_minimal,unit=unt)
          call yaml_close_stream(unit=unt)
       else
          call yaml_warning('Failed to create input_minimal.yaml, error code='//trim(yaml_toa(ierr)))
       end if
    end if
    if (associated(dict_minimal)) call dict_free(dict_minimal)

    call f_release_routine()

  end subroutine inputs_from_dict

  !> Check the directory of data (create if not present)
  subroutine check_for_data_writing_directory(iproc,in)
    use yaml_output
    use module_base, only: bigdft_mpi
    use f_utils, only: f_zero,f_mkdir
    use wrapper_MPI, only: mpibcast
    use yaml_strings, only: f_strcpy
    implicit none
    integer, intent(in) :: iproc
    type(input_variables), intent(inout) :: in
    !local variables
    logical :: shouldwrite
    character(len=100) :: dirname

    if (iproc==0) call yaml_comment('|',hfill='-')

    !initialize directory name
    !shouldwrite=.false.

    shouldwrite= &!shouldwrite .or. &
                                !in%output_wf_format /= WF_FORMAT_NONE .or. &    !write wavefunctions
         in%output_wf /= ENUM_EMPTY .or. &    !write wavefunctions
         in%output_denspot /= ENUM_EMPTY .or. & !write output density
         in%ncount_cluster_x > 1 .or. &                  !write posouts or posmds
         !in%inputPsiId == 2 .or. &                       !have wavefunctions to read
         !in%inputPsiId == 12 .or.  &                     !read in gaussian basis
         (in%inputPsiId .hasattr. 'FILE') .or. &
         (in%inputPsiId .hasattr. 'GAUSSIAN') .or. &
         !in%gaussian_help .or. &                         !Mulliken and local density of states
         bigdft_mpi%ngroup > 1   .or. &                  !taskgroups have been inserted
         mod(in%lin%plotBasisFunctions,10) > 0 .or. &    !dumping of basis functions for locreg runs
         !in%inputPsiId == 102 .or. &                     !reading of basis functions
         in%write_orbitals>0 .or. &                      !writing the KS orbitals in the linear case
         mod(in%lin%output_mat_format,10)>0 .or. &       !writing the sparse linear matrices
         mod(in%lin%output_coeff_format,10)>0            !writing the linear KS coefficients

    !here you can check whether the etsf format is compiled

    if (shouldwrite) then
       !!call create_dir_output(iproc, in)
       call f_zero(dirname)
       if (iproc == 0) then
          call f_mkdir(in%dir_output,dirname)
       end if
       call mpibcast(dirname,comm=bigdft_mpi%mpi_comm)
       !in%dir_output=dirname
       call f_strcpy(src=dirname,dest=in%dir_output)
       if (iproc==0) call yaml_map('Data Writing directory',trim(in%dir_output))
    else
       if (iproc==0) call yaml_map('Data Writing directory','./')
       call f_zero(in%dir_output)!=repeat(' ',len(in%dir_output))
    end if

  END SUBROUTINE check_for_data_writing_directory


  !> Fill all the input keys into dict
  subroutine input_keys_fill_all(dict,dict_minimal)
    use dynamic_memory
    use module_defs, only: gp, pi_param
    use f_input_file
    use public_keys
    use yaml_strings, only: operator(.eqv.)
    use yaml_output
    !use yaml_output
    implicit none
    type(dictionary), pointer :: dict,dict_minimal
    !local variables
    type(dictionary), pointer :: as_is,nested,no_check
    character(max_field_length) :: meth, prof
    real(gp) :: dtmax_, betax_
    logical :: user_defined,free,dftvar

    if (f_err_raise(.not. associated(dict),'The input dictionary has to be associated',&
         err_name='BIGDFT_RUNTIME_ERROR')) return

    call f_routine(id='input_keys_fill_all')

    ! Overriding the default for isolated system
    if ((POSINP .in. dict) .and. (DFT_VARIABLES .in. dict) ) then
       free=ASTRUCT_CELL .notin. dict//POSINP
       dftvar=DISABLE_SYM .notin. dict//DFT_VARIABLES
       if (free .and. dftvar) then
          call set(dict // DFT_VARIABLES // DISABLE_SYM,.true.)
       end if
    end if
    nested=>list_new(.item. LIN_BASIS_PARAMS)



    ! Check and complete dictionary.
    call input_keys_init()
! call yaml_map('present status',dict)
    call input_file_complete(parameters,dict,imports=profiles,nocheck=nested)



    !create a shortened dictionary which will be associated to the given run
    !call input_minimal(dict,dict_minimal)
    as_is =>list_new(.item. FRAG_VARIABLES,.item. IG_OCCUPATION, .item. POSINP, .item. OCCUPATION)
    call input_file_minimal(parameters,dict,dict_minimal,nested,as_is)
    call dict_free(nested,as_is)


    ! Additional treatments.
    meth = dict // GEOPT_VARIABLES // GEOPT_METHOD
    if (trim(meth) .eqv. "FIRE") then
       !prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTMAX, user_defined)
       !if (trim(prof) == DEFAULT .and. .not. user_defined) then
       if (input_value_is_default(dict // GEOPT_VARIABLES, DTMAX)) then
          betax_ = dict // GEOPT_VARIABLES // BETAX
          call set(dict // GEOPT_VARIABLES // DTMAX, 0.25 * pi_param * sqrt(betax_), fmt = "(F7.4)")
       end if
       !prof = input_keys_get_source(dict // GEOPT_VARIABLES, DTINIT, user_defined)
       !if (trim(prof) == DEFAULT .and. .not. user_defined) then
       if (input_value_is_default(dict // GEOPT_VARIABLES, DTINIT)) then
          dtmax_ = dict // GEOPT_VARIABLES // DTMAX
          call set(dict // GEOPT_VARIABLES // DTINIT, 0.5 * dtmax_, fmt = "(F7.4)")
       end if
    end if

    call input_keys_finalize()

    call f_release_routine()
  end subroutine input_keys_fill_all

  !> takes the posinp filename from the dictionary. Starting point is dict//POSINP
  subroutine astruct_dict_get_source(dict, source)
    use public_keys, only: POSINP_SOURCE
    use f_utils, only: f_zero
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(inout) :: source !<preserve previous value if present
    !local variables
    type(dictionary), pointer :: dict_tmp

    call f_zero(source)
    dict_tmp=dict .get. ASTRUCT_PROPERTIES
    source=dict_tmp .get. POSINP_SOURCE

!!$    write(source, "(A)") ""
!!$    if (has_key(dict, ASTRUCT_PROPERTIES)) then
!!$       if (has_key(dict // ASTRUCT_PROPERTIES, POSINP_SOURCE)) &
!!$            & source = dict_value(dict // ASTRUCT_PROPERTIES // POSINP_SOURCE)
!!$    end if
  end subroutine astruct_dict_get_source

  !> Dump the dictionary of the input variables.
  !! Should dump only the keys relative to the input variables and
  !! print out warnings for the ignored keys
  subroutine input_keys_dump(dict, userOnly)
    use yaml_strings, only: f_strcpy
    use yaml_output
    use f_input_file, only: input_file_dump
    use dynamic_memory
    use public_keys, only: POSINP,ASTRUCT_PROPERTIES,ASTRUCT_POSITIONS
    
    implicit none
    type(dictionary), pointer :: dict   !< Dictionary to dump
    logical, intent(in), optional :: userOnly

    !local variables
    integer, parameter :: natoms_dump=500
    integer :: i, dlen, skeys,natoms
    character(max_field_length), dimension(:), allocatable :: keys
    character(max_field_length) ::  sourcefile
    logical :: userOnly_
    type(dictionary), pointer :: tmp

    call f_routine(id='input_keys_dump')


    !new mechanism, to see if it works

    userOnly_ = .false.
    if (present(userOnly)) userOnly_ = userOnly

    !create nodump_list
    tmp => list_new(.item. POSINP)
    call input_file_dump(dict, useronly=userOnly_,nodump_list=tmp)

    call dict_free(tmp)
    !then treat separately the posinp list
    call f_strcpy(src='not provided',dest=sourcefile)
    tmp=dict .get. POSINP
    if ( .not. associated(tmp)) return
    !check the number of atoms
    natoms=-1
    if (ASTRUCT_POSITIONS .in. tmp) natoms=dict_len(tmp // ASTRUCT_POSITIONS)
    if (natoms > natoms_dump) then
       call astruct_dict_get_source(tmp, sourcefile)
       call yaml_map(POSINP,sourcefile)
    else
       call yaml_mapping_open(POSINP)
       call input_file_dump(tmp, useronly=userOnly_,msg='Atomic positions')
       call yaml_mapping_close()
    end if
    nullify(tmp)


!!$    call yaml_comment("Input parameters", hfill = "-")
!!$    !TEST (the first dictionary has no key)
!!$    !if (.not. associated(dict%parent)) then
!!$    if (associated(dict%child)) then
!!$       if (dict_len(dict) >= 1) then
!!$          ! List case.
!!$          dlen = dict_len(dict)
!!$          do i = 0,  dlen- 1, 1
!!$             call input_variable_dump(dict // i,userOnly_)
!!$          end do
!!$       else
!!$          ! Dictionary case
!!$          allocate(keys(dict_size(dict)))
!!$          keys = dict_keys(dict)
!!$          skeys = size(keys)
!!$          do i = 1, skeys
!!$             if (POSINP == trim(keys(i))) then
!!$                call f_strcpy(src='not provided',dest=sourcefile)
!!$                tmp=>dict//POSINP
!!$                !check the number of atoms
!!$                natoms=-1
!!$                if (ASTRUCT_POSITIONS .in. tmp) natoms=dict_len(tmp // ASTRUCT_POSITIONS)
!!$                if (natoms > natoms_dump) then
!!$                   call astruct_dict_get_source(tmp, sourcefile)
!!$                   call yaml_map(POSINP,sourcefile)
!!$                else
!!$                   call input_variable_dump(dict // keys(i),userOnly_)
!!$                end if
!!$             else
!!$                call input_variable_dump(dict // keys(i),userOnly_)
!!$             end if
!!$          end do
!!$          deallocate(keys)
!!$       end if
!!$    else
!!$       call yaml_scalar(dict%data%value)
!!$    end if

    call f_release_routine()

  end subroutine input_keys_dump

  subroutine input_set_int(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    integer, intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_int


  subroutine input_set_dbl(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    real(gp), intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_dbl


  subroutine input_set_bool(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    logical, intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_bool

  subroutine input_set_char(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    character(len = *), intent(in) :: val

    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    call set(dict // key(sep + 1:), val)
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_char


  subroutine input_set_int_array(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    integer, dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_int_array


  subroutine input_set_dbl_array(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    real(gp), dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_dbl_array


  subroutine input_set_bool_array(in, key, val)
    implicit none
    type(input_variables), intent(inout) :: in
    character(len = *), intent(in) :: key
    logical, dimension(:), intent(in) :: val

    integer :: i
    type(dictionary), pointer :: dict
    integer :: sep

    sep = index(key, "/")
    if (sep == 0) stop "no level in key"

    call dict_init(dict)
    do i = 0, size(val) - 1, 1
       call set(dict // key(sep + 1:) // i, val(i + 1))
    end do
    call input_set(in, key(1:sep - 1), dict)
    call dict_free(dict)
  END SUBROUTINE input_set_bool_array

  !> set the inputpsiid enumerator, based on the 
  !! profile input variable
  subroutine set_inputpsiid(profile,inputPsiid)
    use yaml_strings, only: yaml_toa
    implicit none
    integer, intent(in) :: profile
    type(f_enumerator), intent(out) :: inputPsiid

    select case(profile)
    case(INPUT_PSI_EMPTY        )
       inputPsiid=f_enumerator('INPUT_PSI_EMPTY',INPUT_PSI_EMPTY,null())
       call f_enum_attr(inputPsiid,ENUM_SCRATCH)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
    case(INPUT_PSI_RANDOM       )
       inputPsiid=f_enumerator('INPUT_PSI_RANDOM',INPUT_PSI_RANDOM,null())
       call f_enum_attr(inputPsiid,ENUM_SCRATCH)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
    case(INPUT_PSI_CP2K         )
       inputPsiid=f_enumerator('INPUT_PSI_CP2K',INPUT_PSI_CP2K,null())
       call f_enum_attr(inputPsiid,ENUM_FILE)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
       call f_enum_attr(inputPsiid,ENUM_GAUSSIAN)
    case(INPUT_PSI_LCAO         )
       inputPsiid=f_enumerator('INPUT_PSI_LCAO',INPUT_PSI_LCAO,null())
       call f_enum_attr(inputPsiid,ENUM_SCRATCH)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
    case(INPUT_PSI_DISK_WVL     )
       inputPsiid=f_enumerator('INPUT_PSI_DISK_WVL',INPUT_PSI_DISK_WVL,null())
       call f_enum_attr(inputPsiid,ENUM_FILE)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
    case(INPUT_PSI_LCAO_GAUSS   )
       inputPsiid=f_enumerator('INPUT_PSI_LCAO_GAUSS',INPUT_PSI_LCAO_GAUSS,null())
       call f_enum_attr(inputPsiid,ENUM_SCRATCH)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
       call f_enum_attr(inputPsiid,ENUM_GAUSSIAN)
    case(INPUT_PSI_DISK_GAUSS   )
       inputPsiid=f_enumerator('INPUT_PSI_DISK_GAUSS',INPUT_PSI_DISK_GAUSS,null())
       call f_enum_attr(inputPsiid,ENUM_FILE)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
       call f_enum_attr(inputPsiid,ENUM_GAUSSIAN)
    case(INPUT_PSI_LINEAR_AO    )
       inputPsiid=f_enumerator('INPUT_PSI_LINEAR_AO',INPUT_PSI_LINEAR_AO,null())
       call f_enum_attr(inputPsiid,ENUM_SCRATCH)
       call f_enum_attr(inputPsiid,ENUM_LINEAR)
    case(INPUT_PSI_DISK_LINEAR  )
       inputPsiid=f_enumerator('INPUT_PSI_DISK_LINEAR',INPUT_PSI_DISK_LINEAR,null())
       call f_enum_attr(inputPsiid,ENUM_FILE)
       call f_enum_attr(inputPsiid,ENUM_LINEAR)
    case(INPUT_PSI_MEMORY_WVL)
       inputPsiid=f_enumerator('INPUT_PSI_MEMORY_WVL',INPUT_PSI_MEMORY_WVL,null())
       call f_enum_attr(inputPsiid,ENUM_MEMORY)
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
    case(INPUT_PSI_MEMORY_GAUSS)
       inputPsiid=f_enumerator('INPUT_PSI_MEMORY_GAUSS',INPUT_PSI_MEMORY_GAUSS,null())
       call f_enum_attr(inputPsiid,ENUM_MEMORY) 
       call f_enum_attr(inputPsiid,ENUM_CUBIC)
       call f_enum_attr(inputPsiid,ENUM_GAUSSIAN) !this attribute should be already in the ENUM_CUBIC
    case(INPUT_PSI_MEMORY_LINEAR)
       inputPsiid=f_enumerator('INPUT_PSI_MEMORY_LINEAR',INPUT_PSI_MEMORY_LINEAR,null())
       call f_enum_attr(inputPsiid,ENUM_MEMORY) 
       call f_enum_attr(inputPsiid,ENUM_LINEAR)
    case default !the other are from memory
       call f_err_throw(err_name="BIGDFT_INPUT_VARIABLES_ERROR",err_msg='The value '//trim(yaml_toa(profile))//&
            ' for the inputpsiid is unknown')
    end select

  end subroutine set_inputpsiid

  !> Get the current run policy
  subroutine inputpsiid_get_policy(inputpsiid, policy)
    implicit none
    ! Calling arguments
    type(f_enumerator), intent(in) :: inputPsiid
    type(f_enumerator), intent(out) :: policy

    if (inputPsiId .hasattr. 'SCRATCH') then
       policy = ENUM_SCRATCH
    else if (inputPsiId .hasattr. 'MEMORY') then
       policy = ENUM_MEMORY
    else if (inputPsiId .hasattr. 'FILE') then
       policy = ENUM_FILE
    end if
  end subroutine inputpsiid_get_policy

  !> Set the current run policy
  !! modify the values of inputpsiid according to
  !! the previously chosen policy
  subroutine inputpsiid_set_policy(policy,inputpsiid)
    implicit none
    ! Calling arguments
    type(f_enumerator), intent(in) :: policy
    type(f_enumerator), intent(inout) :: inputPsiid

    ! Set the new ID for the input guess
    !select case(policy)
    !case(INPUT_POLICY_SCRATCH)
    if (policy == ENUM_SCRATCH) then
       if (inputPsiId .hasattr. 'CUBIC') then
          !if the gaussian help is activated, the scratch mode
          !will become a restart from memory
          if (inputPsiId .hasattr. 'GAUSSIAN') then
             call set_inputpsiid(INPUT_PSI_MEMORY_GAUSS,&
                  inputPsiId)
          else
             call set_inputpsiid(INPUT_PSI_LCAO,inputPsiId)
          end if
       else if (inputPsiId .hasattr. 'LINEAR') then
          call set_inputpsiid(INPUT_PSI_LINEAR_AO,&
               inputPsiId)
       end if
    !case(INPUT_POLICY_MEMORY)
    else if (policy == ENUM_MEMORY) then
       if (inputPsiId .hasattr. 'CUBIC') then
          call set_inputpsiid(INPUT_PSI_MEMORY_WVL,&
               inputPsiId)
       else if (inputPsiId .hasattr. 'LINEAR') then
          call set_inputpsiid(INPUT_PSI_MEMORY_LINEAR,&
               inputPsiId)
       end if
    !case(INPUT_POLICY_DISK)
    else if (policy == ENUM_FILE) then
       if (inputPsiId .hasattr. 'CUBIC') then
          call set_inputpsiid(INPUT_PSI_DISK_WVL,&
               inputPsiId)
       else if (inputPsiId .hasattr. 'LINEAR') then
          call set_inputpsiid(INPUT_PSI_DISK_LINEAR,&
               inputPsiId)
       end if
    end if
  end subroutine inputpsiid_set_policy


  subroutine set_output_wf(profile,output_wf)
    implicit none
    integer, intent(in) :: profile
    type(f_enumerator), intent(out) :: output_wf

    select case(profile)
    case(WF_FORMAT_NONE)
       output_wf=ENUM_EMPTY
    case(WF_FORMAT_PLAIN)
       output_wf=ENUM_TEXT
    case(WF_FORMAT_BINARY)
       output_wf=ENUM_BINARY
    case(WF_FORMAT_ETSF)
       output_wf=ENUM_ETSF
    end select

  end subroutine set_output_wf

  subroutine set_output_denspot(profile,output_denspot)
    implicit none
    integer, intent(in) :: profile
    type(f_enumerator), intent(out) :: output_denspot
    select case(modulo(abs(profile),10))
    case( OUTPUT_DENSPOT_NONE)
       output_denspot=ENUM_EMPTY
    case (OUTPUT_DENSPOT_DENSITY)
       output_denspot=f_enumerator('DENSITY',OUTPUT_DENSPOT_DENSITY,null())
    case (OUTPUT_DENSPOT_DENSPOT)
       output_denspot=f_enumerator('DENSPOT',OUTPUT_DENSPOT_DENSPOT,null())
    end select
    !then set the attribute
    select case(profile/10)
    case( OUTPUT_DENSPOT_FORMAT_TEXT)
       if (profile /= OUTPUT_DENSPOT_NONE)  call f_enum_attr(output_denspot,ENUM_TEXT)
    case( OUTPUT_DENSPOT_FORMAT_ETSF)
       call f_enum_attr(output_denspot,ENUM_ETSF)
    case( OUTPUT_DENSPOT_FORMAT_CUBE)
       call f_enum_attr(output_denspot,ENUM_CUBE)
    end select
  end subroutine set_output_denspot


  !> Set the dictionary from the input variables
  subroutine input_set_dict(in, level, val)
    use module_defs, only: DistProjApply, gp
    use wrapper_linalg, only: GPUblas
    use public_enums
    use dynamic_memory
    use yaml_output, only: yaml_warning
    use yaml_strings, only: operator(.eqv.),is_atoi
    use module_base, only: bigdft_mpi
    implicit none
    type(input_variables), intent(inout) :: in
    type(dictionary), pointer :: val
    character(len = *), intent(in) :: level
    integer, dimension(2) :: dummy_int !<to use as filling for input variables
    real(gp), dimension(3) :: dummy_gp !< to fill the input variables
    logical, dimension(2) :: dummy_log !< to fill the input variables
    character(len=256) :: dummy_char
    character(len = max_field_length) :: str
    integer :: i, ipos

    if (index(dict_key(val), "_attributes") > 0) return

    select case(trim(level))
    case(MODE_VARIABLES)
       select case (trim(dict_key(val)))
       case(METHOD_KEY)
          str=val
          select case(trim(str))
          case('lj')
             in%run_mode=LENNARD_JONES_RUN_MODE
          case('dft')
             in%run_mode=QM_RUN_MODE
          case('lensic')
             in%run_mode=LENOSKY_SI_CLUSTERS_RUN_MODE
          case('lensib')
             in%run_mode=LENOSKY_SI_BULK_RUN_MODE
          case('amber')
             in%run_mode=AMBER_RUN_MODE
          case('morse_bulk')
             in%run_mode=MORSE_BULK_RUN_MODE
          case('morse_slab')
             in%run_mode=MORSE_SLAB_RUN_MODE
          case('tersoff')
             in%run_mode=TERSOFF_RUN_MODE
          case('bmhtf')
             in%run_mode=BMHTF_RUN_MODE
          case('cp2k')
             in%run_mode=CP2K_RUN_MODE
          case('dftbp')
             in%run_mode=DFTBP_RUN_MODE
          case('multi')
             in%run_mode=MULTI_RUN_MODE
          end select
       case(MM_PARAMSET)
          in%mm_paramset=val
       case(MM_PARAMFILE)
          in%mm_paramfile=val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (DFT_VARIABLES)
       ! the DFT variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (HGRIDS)
          !grid spacings (profiles can be used if we already read PSPs)
          dummy_gp(1:3)=val
          in%hx = dummy_gp(1)
          in%hy = dummy_gp(2)
          in%hz = dummy_gp(3)
       case (RMULT)
          !coarse and fine radii around atoms
          dummy_gp(1:2)=val
          in%crmult = dummy_gp(1)
          in%frmult = dummy_gp(2)
       case (IXC)
          in%ixc = val !XC functional (ABINIT XC codes)
       case (NCHARGE) !charge 
          str=val
          !check if the provided value is a integer
          if (is_atoi(str)) then
             ipos=val
             in%qcharge=real(ipos,gp) !exact conversion
          else
             in%qcharge = val 
          end if
       case (ELECFIELD) !electric field
          in%elecfield = val
       case (NSPIN)
          in%nspin = val !spin and polarization
       case (MPOL)
          in%mpol = val
       case (GNRM_CV)
          in%gnrm_cv = val !convergence parameters
       case (ITERMAX)
          in%itermax = val
       case (ITERMIN)
          in%itermin = val
       case (NREPMAX)
          in%nrepmax = val
       case (NCONG)
          in%ncong = val !convergence parameters
       case (IDSX)
          in%idsx = val
       case (DISPERSION)
          in%dispersion = val !dispersion parameter
          ! Now the variables which are to be used only for the last run
       case (INPUTPSIID)
          ipos=val
          call set_inputpsiid(ipos,in%inputPsiId)
          !in%inputPsiId = val
       case (OUTPUT_WF)
          ipos=val
          call set_output_wf(ipos,in%output_wf)
          !in%output_wf = val
       case (OUTPUT_DENSPOT)
          ipos=val
          call set_output_denspot(ipos,in%output_denspot)
          !in%output_denspot = val
       case (RBUF)
          in%rbuf = val ! Tail treatment.
       case (NCONGT)
          in%ncongt = val
       case (NORBV)
          in%norbv = val !davidson treatment
       case (NVIRT)
          in%nvirt = val
       case (NPLOT)
          in%nplot = val
       case (DISABLE_SYM)
          in%disableSym = val ! Line to disable symmetries.
       case (SOLVENT)
          in%set_epsilon= val
!!$          dummy_char = val
!!$          select case(trim(dummy_char))
!!$          case ("vacuum")
!!$             in%set_epsilon =EPSILON_VACUUM
!!$          case("rigid")
!!$             in%set_epsilon =EPSILON_RIGID_CAVITY
!!$          case("sccs")
!!$             in%set_epsilon =EPSILON_SCCS
!!$          end select
       case (EXTERNAL_POTENTIAL)
          ! Do nothing?
       case(CALCULATE_STRTEN)
          in%calculate_strten=val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (PERF_VARIABLES)
       ! the PERF variables -----------------------------------------------------
       select case (trim(dict_key(val)))       
       case (DEBUG)
          in%debug = val
       case (FFTCACHE)
          in%ncache_fft = val
       case (VERBOSITY)
          in%verbosity = val
       case (TOLSYM)
          in%symTol = val
       case (PROJRAD)
          in%projrad = val
       case (EXCTXPAR)
          in%exctxpar = val
       case (INGUESS_GEOPT)
          in%inguess_geopt = val
       case (ACCEL)
          str = dict_value(val)
          if (trim(str) .eqv. "CUDAGPU") then
             in%matacc%iacceleration = 1
          else if (trim(str) .eqv. "OCLGPU") then
             in%matacc%iacceleration = 2
          else if (trim(str) .eqv. "OCLCPU") then
             in%matacc%iacceleration = 3
          else if (trim(str) .eqv. "OCLACC") then
             in%matacc%iacceleration = 4
          else 
             in%matacc%iacceleration = 0
          end if
       case (OCL_PLATFORM)
          !determine desired OCL platform which is used for acceleration
          in%matacc%OCL_platform = val
          ipos=min(len(in%matacc%OCL_platform),len(trim(in%matacc%OCL_platform))+1)
          if (index(in%matacc%OCL_platform,'~') > 0) ipos=1 !restore empty information if not entered by the user
          do i=ipos,len(in%matacc%OCL_platform)
             in%matacc%OCL_platform(i:i)=achar(0)
          end do
       case (OCL_DEVICES)
          in%matacc%OCL_devices = val
          ipos=min(len(in%matacc%OCL_devices),len(trim(in%matacc%OCL_devices))+1)
          if (index(in%matacc%OCL_devices,'~') > 0) ipos=1
          do i=ipos,len(in%matacc%OCL_devices)
             in%matacc%OCL_devices(i:i)=achar(0)
          end do
       case (PSOLVER_ACCEL)
          in%matacc%PSolver_igpu = val
       case (SIGNALING)
          in%signaling = val ! Signaling parameters
       case (SIGNALTIMEOUT)
          in%signalTimeout = val
       case (DOMAIN)
          in%domain = val
       case (BLAS)
          GPUblas = val !!@TODO to relocate
       case (PSP_ONFLY)
          DistProjApply = val
       case (MULTIPOLE_PRESERVING)
          in%multipole_preserving = val
       case (MP_ISF)
          in%mp_isf = val
       case (IG_DIAG)
          in%orthpar%directDiag = val
       case (IG_NORBP)
          in%orthpar%norbpInguess = val
       case (IG_TOL)
          in%orthpar%iguessTol = val
       case (METHORTHO)
          in%orthpar%methOrtho = val
       case (IG_BLOCKS)
          !Block size used for the orthonormalization
          dummy_int(1:2)=val
          in%orthpar%bsLow = dummy_int(1)
          in%orthpar%bsUp  = dummy_int(2)
       case (RHO_COMMUN)
          in%rho_commun = val
       case (PSOLVER_GROUPSIZE)
          in%PSolver_groupsize = val
       case (UNBLOCK_COMMS)
          in%unblock_comms = val
       case (LINEAR)
          !Use Linear scaling methods
          str = dict_value(val)
          if (trim(str) .eqv. "LIG") then
             in%linear = INPUT_IG_LIG
          else if (trim(str) .eqv. "FUL") then
             in%linear = INPUT_IG_FULL
          else if (trim(str) .eqv. "TMO") then
             in%linear = INPUT_IG_TMO
          else
             in%linear = INPUT_IG_OFF
          end if
       case (STORE_INDEX)
          in%store_index = val
       case (PDSYEV_BLOCKSIZE)
          !block size for pdsyev/pdsygv, pdgemm (negative -> sequential)
          in%lin%blocksize_pdsyev = val
       case (PDGEMM_BLOCKSIZE)
          in%lin%blocksize_pdgemm = val
       case (MAXPROC_PDSYEV)
          !max number of process uses for pdsyev/pdsygv, pdgemm
          in%lin%nproc_pdsyev = val
       case (MAXPROC_PDGEMM)
          in%lin%nproc_pdgemm = val
       case (EF_INTERPOL_DET)
          !FOE: if the determinant of the interpolation matrix to find the Fermi energy
          !is smaller than this value, switch from cubic to linear interpolation.
          in%lin%ef_interpol_det = val
       case (EF_INTERPOL_CHARGEDIFF)
          in%lin%ef_interpol_chargediff = val
          !determines whether a mixing step shall be preformed after the input guess !(linear version)
       case (MIXING_AFTER_INPUTGUESS)
          in%lin%mixing_after_inputguess = val
          !determines whether the input guess support functions are orthogonalized iteratively (T) or in the standard way (F)
       case (ITERATIVE_ORTHOGONALIZATION)
          in%lin%iterative_orthogonalization = val
       case (CHECK_SUMRHO)
          in%check_sumrho = val
          !  call input_var("mpi_groupsize",0, "number of MPI processes for BigDFT run (0=nproc)", in%mpi_groupsize)
       case (CHECK_OVERLAP)
          ! perform a check of the overlap calculation
          in%check_overlap = val
       case (EXPERIMENTAL_MODE)
          in%experimental_mode = val
       case (WRITE_ORBITALS)
          ! linear scaling: write KS orbitals for cubic restart
          in%write_orbitals = val
       case (EXPLICIT_LOCREGCENTERS)
          ! linear scaling: explicitely specify localization centers
          in%explicit_locregcenters = val
       case (CALCULATE_KS_RESIDUE)
          ! linear scaling: calculate Kohn-Sham residue
          in%calculate_KS_residue = val
       case (INTERMEDIATE_FORCES)
          ! linear scaling: calculate intermediate forces
          in%intermediate_forces = val
       case (KAPPA_CONV)
          ! linear scaling: exit kappa for extended input guess (experimental mode)
          in%kappa_conv = val
       case (EVBOUNDS_NSATUR)
          ! linear scaling: number of FOE cycles before the eigenvalue bounds are shrinked
          in%evbounds_nsatur = val
       case(EVBOUNDSSHRINK_NSATUR)
          ! linear scaling: maximal number of unsuccessful eigenvalue bounds shrinkings
          in%evboundsshrink_nsatur = val
       case (METHOD_UPDATEKERNEL)
          ! linear scaling: how to update the density kernel during the support function optimization (0: purification, 1: FOE)
          in%method_updatekernel = val
       case (PURIFICATION_QUICKRETURN)
          ! linear scaling: quick return in purification
          in%purification_quickreturn = val
       case (ADJUST_foe_TEMPERATURE)
          ! linear scaling: dynamic adjustment of the decay length of the FOE error function
          in%adjust_FOE_temperature = val
       case (CALCULATE_GAP)
          ! linear scaling: calculate the HOMO LUMO gap even when FOE is used for the kernel calculation
          in%calculate_gap = val
       case (LOEWDIN_CHARGE_ANALYSIS)
          ! linear scaling: perform a Loewdin charge analysis at the end of the calculation
          in%loewdin_charge_analysis = val
       case (COEFF_WEIGHT_ANALYSIS)
          ! linear scaling: perform a Loewdin charge analysis of the coefficients for fragment calculations
          in%coeff_weight_analysis = val
       case (CHECK_MATRIX_COMPRESSION)
          ! linear scaling: perform a check of the matrix compression routines
          in%check_matrix_compression = val
       case (CORRECTION_CO_CONTRA)
          ! linear scaling: correction covariant / contravariant gradient
          in%correction_co_contra = val
       case (FSCALE_LOWERBOUND)
          ! linear scaling: lower bound for the error function decay length
          in%fscale_lowerbound = val
       case (FSCALE_UPPERBOUND)
          ! linear scaling: upper bound for the error function decay length
          in%fscale_upperbound = val
       case (FOE_RESTART)
          ! linear scaling: Restart method to be used for the FOE method
          in%FOE_restart = val
       case (IMETHOD_OVERLAP)
          ! linear scaling: method to calculate the overlap matrices (1=old, 2=new)
          in%imethod_overlap = val
       case (ENABLE_MATRIX_TASKGROUPS) 
          ! linear scaling: enable the matrix taskgroups
          in%enable_matrix_taskgroups = val
       case (HAMAPP_RADIUS_INCR)
          ! linear scaling: radius enlargement for the Hamiltonian application (in grid points)
          in%hamapp_radius_incr = val
       case (ADJUST_KERNEL_ITERATIONS) 
          ! linear scaling: enable the addaptive ajustment of the number of kernel iterations
          in%adjust_kernel_iterations = val
       case(WF_EXTENT_ANALYSIS)
          ! linear scaling: perform an analysis of the extent of the support functions (and possibly KS orbitals)
          in%wf_extent_analysis = val
       case (GPS_METHOD)
          in%GPS_method = val
       case (FOE_GAP)
          ! linear scaling: Use the FOE method to calculate the HOMO-LUMO gap at the end
          in%foe_gap = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case ("geopt")

       ! the GEOPT variables ----------------------------------------------------
       select case (trim(dict_key(val)))
       case (GEOPT_METHOD)
          in%geopt_approach = val !geometry input parameters
       case (NCOUNT_CLUSTER_X)
          in%ncount_cluster_x = val
       case (FRAC_FLUCT)
          in%frac_fluct = val
       case (FORCEMAX)
          in%forcemax = val
       case (RANDDIS)
          in%randdis = val
       case (IONMOV)
          in%ionmov = val
       case (DTION)
          in%dtion = val
       case (MDITEMP)
          in%mditemp = val
       case (MDFTEMP)
          in%mdftemp = val
       case (NOSEINERT)
          in%noseinert = val
       case (FRICTION)
          in%friction = val
       case (MDWALL)
          in%mdwall = val
       case (QMASS)
          in%nnos = dict_len(val)
          call f_free_ptr(in%qmass)
          in%qmass = f_malloc_ptr(in%nnos, id = "in%qmass")
          do i=1,in%nnos-1
             in%qmass(i) = dict_len(val // (i-1))
          end do
       case (BMASS)
          in%bmass = val
       case (VMASS)
          in%vmass = val
       case (BETAX)
          in%betax = val
       case (HISTORY)
          in%history = val
       case (DTINIT)
          in%dtinit = val
       case (DTMAX)
          in%dtmax = val
       case (NHISTX)
          in%nhistx = val
       case (MAXRISE)
          in%maxrise = val
       case (CUTOFFRATIO)
          in%cutoffratio = val
       case (STEEPTHRESH)
          in%steepthresh = val
       case (BIOMODE)
          in%biomode = val
       case (BETA_STRETCHX)
          in%beta_stretchx = val
       case (TRUSTR)
          in%trustr = val
       case (SOCKINET)
          in%sockinet = val
       case (SOCKPORT)
          in%sockport = val
       case (SOCKHOST)
          in%sockhost = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
!NNdbg
    case (MD_VARIABLES)
       select case (trim(dict_key(val)))
       case (MDSTEPS)
          in%mdsteps = val
       case (PRINT_FREQUENCY)
          in%md_printfrq = val
       case (TEMPERATURE)
          in%temperature = val
       case (TIMESTEP)
          in%dt = val
       case (NO_TRANSLATION) !.true. or .false. ?
          in%no_translation = val
       case (THERMOSTAT) !string
         str = dict_value(val) 
         if (trim(str).eqv."nose_hoover_chain") then
           in%nhc=.true.
         else
           in%nhc=.false.
         end if
       case (NOSE_CHAIN_LENGTH) 
         in%nhnc = val
       case (NOSE_MTS_SIZE)
         in%nmultint = val
       case (NOSE_YOSHIDA_FACTOR)
         in%nsuzuki = val
       case (NOSE_FREQUENCY)
         in%nosefrq = val
       case (WAVEFUNCTION_EXTRAPOLATION)
          in%wfn_history = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
!NNdbg
    case (MIX_VARIABLES)
       ! the MIX variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (ISCF)
          in%iscf = val
          !put the startmix if the mixing has to be done
          if (in%iscf >  SCF_KIND_DIRECT_MINIMIZATION) in%gnrm_startmix=1.e300_gp
       case (ITRPMAX)
          in%itrpmax = val
       case (RPNRM_CV)
          in%rpnrm_cv = val
       case (NORBSEMPTY)
          in%norbsempty = val
       case (TEL)
          in%Tel = val
       case (OCCOPT)
          in%occopt = val
       case (ALPHAMIX)
          in%alphamix = val
       case (ALPHADIIS)
          in%alphadiis = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (SIC_VARIABLES)
       ! the SIC variables ------------------------------------------------------
       select case (trim(dict_key(val)))
       case (SIC_APPROACH)
          in%SIC%approach = val
       case (SIC_ALPHA)
          in%SIC%alpha = val
       case (SIC_FREF)
          in%SIC%fref = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (TDDFT_VARIABLES)
       ! the TDDFT variables ----------------------------------------------------
       select case (trim(dict_key(val)))
       case (TDDFT_APPROACH)
          in%tddft_approach = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (LIN_GENERAL)
       ! the variables for the linear version, general section
       select case (trim(dict_key(val)))
       case (HYBRID)
          dummy_log(1) = val
          if (dummy_log(1)) then
             in%lin%nlevel_accuracy = 1
          else
             in%lin%nlevel_accuracy = 2
          end if
       case (NIT)
          dummy_int(1:2) = val
          in%lin%nit_lowaccuracy = dummy_int(1)
          in%lin%nit_highaccuracy = dummy_int(2)
       case (RPNRM_CV)
          dummy_gp(1:2) = val
          in%lin%lowaccuracy_conv_crit = dummy_gp(1)
          in%lin%highaccuracy_conv_crit = dummy_gp(2)
       case (CONF_DAMPING) 
          in%lin%reduce_confinement_factor = val
       case (TAYLOR_ORDER)
          in%lin%order_taylor = val
       case (OUTPUT_WF)
          in%lin%plotBasisFunctions = val
       case (OUTPUT_MAT)
          in%lin%output_mat_format = val
       case (OUTPUT_COEFF)
          in%lin%output_coeff_format = val
       case (OUTPUT_FRAGMENTS)
          in%lin%output_fragments = val
       case (KERNEL_RESTART_MODE)
          in%lin%kernel_restart_mode = val
       case (CALC_DIPOLE)
          in%lin%calc_dipole = val
       case (CALC_PULAY)
          dummy_log(1:2) = val
          in%lin%pulay_correction = dummy_log(1)
          in%lin%new_pulay_correction = dummy_log(2)
       case (SUBSPACE_DIAG)
          in%lin%diag_end = val
       case (EXTRA_STATES)
          in%lin%extra_states = val
       case (MAX_INVERSION_ERROR)
          ! maximal error of the Taylor approximations to calculate the inverse of the overlap matrix
          in%lin%max_inversion_error = val
       case (CALCULATE_ONSITE_OVERLAP)
          in%lin%calculate_onsite_overlap = val
       case (CHARGE_MULTIPOLES)
           in%lin%charge_multipoles = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (LIN_BASIS)
       select case (trim(dict_key(val)))
       case (NIT)
          dummy_int(1:2) = val
          in%lin%nItBasis_lowaccuracy = dummy_int(1)
          in%lin%nItBasis_highaccuracy = dummy_int(2)
       case (IDSX)
          dummy_int(1:2) = val
          in%lin%DIIS_hist_lowaccur = dummy_int(1)
          in%lin%DIIS_hist_highaccur = dummy_int(2)
       case (GNRM_CV)
          dummy_gp(1:2) = val
          in%lin%convCrit_lowaccuracy = dummy_gp(1)
          in%lin%convCrit_highaccuracy = dummy_gp(2)
       case (GNRM_IG)
          in%lin%convCrit_extendedIG = val
       case (NIT_IG)
          in%lin%nit_extendedIG = val
       case (DELTAE_CV)
          in%lin%early_stop = val
       case (GNRM_DYN)
          in%lin%gnrm_dynamic = val
       case (MIN_GNRM_FOR_DYNAMIC)
          in%lin%min_gnrm_for_dynamic = val
       case (ALPHA_DIIS)
          in%lin%alphaDIIS = val
       case (ALPHA_SD)
          in%lin%alphaSD = val
       case (NSTEP_PREC)
          in%lin%nItPrecond = val
       case (fix_basis)
          in%lin%support_functions_converged = val
       case (correction_orthoconstraint)
          in%lin%correctionOrthoconstraint = val
       case DEFAULT
          if (bigdft_mpi%iproc==0) &
               call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
    case (LIN_KERNEL)
       select case (trim(dict_key(val)))
       case (NSTEP)
          dummy_int(1:2) = val
          in%lin%nItdmin_lowaccuracy = dummy_int(1)
          in%lin%nItdmin_highaccuracy = dummy_int(2)
       case (NIT)
          dummy_int(1:2) = val
          in%lin%nItSCCWhenFixed_lowaccuracy = dummy_int(1)
          in%lin%nItSCCWhenFixed_highaccuracy = dummy_int(2)
       case (IDSX_COEFF)
          dummy_int(1:2) = val
          in%lin%dmin_hist_lowaccuracy = dummy_int(1)
          in%lin%dmin_hist_highaccuracy = dummy_int(2)
       case (IDSX)
          dummy_int(1:2) = val
          in%lin%mixHist_lowaccuracy = dummy_int(1)
          in%lin%mixHist_highaccuracy = dummy_int(2)
       case (ALPHAMIX)
          dummy_gp(1:2) = val
          in%lin%alpha_mix_lowaccuracy = dummy_gp(1)
          in%lin%alpha_mix_highaccuracy = dummy_gp(2)
       case (GNRM_CV_COEFF)
          dummy_gp(1:2) = val
          in%lin%convCritdmin_lowaccuracy = dummy_gp(1)
          in%lin%convCritdmin_highaccuracy = dummy_gp(2)
       case (RPNRM_CV)
          dummy_gp(1:2) = val
          in%lin%convCritMix_lowaccuracy = dummy_gp(1)
          in%lin%convCritMix_highaccuracy = dummy_gp(2)
       case (LINEAR_METHOD)
          dummy_char = val
          select case (trim(dummy_char))
          case('DIRMIN')
             in%lin%kernel_mode = KERNELMODE_DIRMIN
          case('DIAG')
             in%lin%kernel_mode = KERNELMODE_DIAG
          case('FOE')
             in%lin%kernel_mode = KERNELMODE_FOE
          end select
       case (MIXING_METHOD)
          dummy_char = val
          select case(trim(dummy_char))
          case('DEN')
             in%lin%mixing_mode = MIXINGMODE_DENS
          case('POT')
             in%lin%mixing_mode = MIXINGMODE_POT
          end select
       case (ALPHA_SD_COEFF)
          in%lin%alphaSD_coeff = val
       case (ALPHA_FIT_COEFF)
          in%lin%curvefit_dmin = val
       case (EVAL_RANGE_FOE)
          dummy_gp(1:2) = val
          in%lin%evlow = dummy_gp(1)
          in%lin%evhigh = dummy_gp(2)
       case (FSCALE_FOE) 
          in%lin%fscale = val
       case DEFAULT
          call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
       end select
       ! the KPT variables ------------------------------------------------------
    case (KPT_VARIABLES)
    case (LIN_BASIS_PARAMS)
    case (OCCUPATION)
    case (IG_OCCUPATION)
    case (FRAG_VARIABLES)
       !case (RUN_NAME_KEY)
    case DEFAULT
    end select
  END SUBROUTINE input_set_dict


  subroutine basis_params_set_dict(dict_basis,lin,jtype)
    implicit none
    integer, intent(in) :: jtype !< local type of which we are filling the values
    type(dictionary), pointer :: dict_basis
    type(linearInputParameters),intent(inout) :: lin
    !local variables
    real(gp), dimension(2) :: dummy_darr
    !-- default parameters for the basis set of linear scaling

    !then update the values of each parameter if present
    select case(trim(dict_key(dict_basis)))
    case(NBASIS)
       lin%norbsPerType(jtype)=dict_basis !npt
    case(AO_CONFINEMENT)
       lin%potentialPrefac_ao(jtype)=dict_basis !ppao
    case(CONFINEMENT)
       dummy_darr=dict_basis
       lin%potentialPrefac_lowaccuracy(jtype)=dummy_darr(1)!ppl
       lin%potentialPrefac_highaccuracy(jtype)=dummy_darr(2)!pph
    case(RLOC)
       dummy_darr=dict_basis
       !locradType(jtype)=dummy_darr(1) !lrl
       lin%locrad_type(jtype,1)=dummy_darr(1) !lrl
       lin%locrad_type(jtype,2)=dummy_darr(2) !lrh
       !locradType_lowaccur(jtype)=dummy_darr(1) !lrl
       !locradType_highaccur(jtype)=dummy_darr(2) !lrh
       !atoms%rloc(jtype,:)=locradType(jtype)
    case(RLOC_KERNEL) 
       lin%kernel_cutoff(jtype)=dict_basis !kco
    case(RLOC_KERNEL_FOE) 
       lin%kernel_cutoff_FOE(jtype)=dict_basis !kco_FOE
    end select

  end subroutine basis_params_set_dict

  !> Creation of the log file (by default log.yaml)
  !>  Free all dynamically allocated memory from the kpt input file.
  subroutine free_kpt_variables(in)
    use dynamic_memory
    implicit none
    type(input_variables), intent(inout) :: in
    character(len=*), parameter :: subname='free_kpt_variables'

    call f_free_ptr(in%gen_kpt)
    call f_free_ptr(in%gen_wkpt)
    call f_free_ptr(in%kptv)
    call f_free_ptr(in%nkptsv_group)
    nullify(in%gen_kpt)
    nullify(in%gen_wkpt)
    nullify(in%kptv)
    nullify(in%nkptsv_group)
  end subroutine free_kpt_variables

  !>  Free all dynamically allocated memory from the geopt input file.
  subroutine free_geopt_variables(in)
    use dynamic_memory
    implicit none
    type(input_variables), intent(inout) :: in
    character(len=*), parameter :: subname='free_geopt_variables'

    call f_free_ptr(in%qmass)
  end subroutine free_geopt_variables
  
  
  !> Set default values for input variables
  subroutine default_input_variables(in)
    use module_base
    use multipole_base, only: external_potential_descriptors_null
    implicit none

    type(input_variables), intent(inout) :: in

    in%refcnt=f_ref_new('inputs')
    in%run_mode=f_enumerator_null()
    in%matacc=material_acceleration_null()

    ! Default values.
    in%dir_output = "data"
    !in%output_wf_format = WF_FORMAT_NONE
    !in%output_denspot_format = output_denspot_FORMAT_CUBE
    call f_zero(in%set_epsilon)
    nullify(in%gen_kpt)
    nullify(in%gen_wkpt)
    nullify(in%kptv)
    nullify(in%nkptsv_group)
    call f_zero(in%calculate_strten)
    in%gen_norb = UNINITIALIZED(0)
    in%gen_norbu = UNINITIALIZED(0)
    in%gen_norbd = UNINITIALIZED(0)
    nullify(in%gen_occup)
    ! Default abscalc variables
    call abscalc_input_variables_default(in)
    ! Default frequencies variables
    call frequencies_input_variables_default(in)
    ! Default values for geopt.
    call geopt_input_variables_default(in) 
    ! Default values for mixing procedure
    call mix_input_variables_default(in) 
    ! Default values for tddft
    call tddft_input_variables_default(in)
    !Default for Self-Interaction Correction variables
    call sic_input_variables_default(in)
    ! Default for signaling
    in%gmainloop = 0.d0
    ! Default for lin.
    nullify(in%lin%potentialPrefac_lowaccuracy)
    nullify(in%lin%potentialPrefac_highaccuracy)
    nullify(in%lin%potentialPrefac_ao)
    nullify(in%lin%norbsPerType)
    nullify(in%lin%locrad)
    nullify(in%lin%locrad_lowaccuracy)
    nullify(in%lin%locrad_highaccuracy)
    nullify(in%lin%locrad_type)
    nullify(in%lin%kernel_cutoff_FOE)
    nullify(in%lin%kernel_cutoff)
    !nullify(in%frag%frag_info)
    nullify(in%frag%label)
    nullify(in%frag%dirname)
    nullify(in%frag%frag_index)
    nullify(in%frag%charge)
    in%ep = external_potential_descriptors_null()
  END SUBROUTINE default_input_variables


  !> Assign default values for mixing variables
  subroutine mix_input_variables_default(in)
    implicit none
    type(input_variables), intent(inout) :: in

    !mixing treatement (hard-coded values)
    in%iscf=0
    in%itrpmax=1
    in%alphamix=0.0_gp
    in%rpnrm_cv=1.e-4_gp
    in%gnrm_startmix=0.0_gp
    in%norbsempty=0
    in%Tel=0.0_gp
    in%occopt=SMEARING_DIST_ERF
    in%alphadiis=2.d0

  END SUBROUTINE mix_input_variables_default


  !> Assign default values for GEOPT variables
  subroutine geopt_input_variables_default(in)
    use module_defs, only: UNINITIALIZED
    implicit none
    type(input_variables), intent(inout) :: in

    !put some fake values for the geometry optimsation case
    in%geopt_approach='SDCG'
    in%ncount_cluster_x=0
    in%frac_fluct=1.0_gp
    in%forcemax=0.0_gp
    in%randdis=0.0_gp
    in%betax=2.0_gp
    in%history = 1
!    in%wfn_history = 1
    in%wfn_history = 0
    in%ionmov = -1
    in%dtion = 0.0_gp
    in%strtarget(:)=0.0_gp
    in%nhistx =0
    in%biomode=.false.
    in%beta_stretchx=0.0_gp
    in%maxrise=0.0_gp
    in%cutoffratio=0.0_gp
    in%steepthresh=0.0_gp
    in%trustr=0.0_gp
    in%mditemp = UNINITIALIZED(in%mditemp)
    in%mdftemp = UNINITIALIZED(in%mdftemp)
    nullify(in%qmass)

  END SUBROUTINE geopt_input_variables_default

  !> Assign default values for MD variables
  subroutine md_input_variables_default(in)
    use module_defs, only: UNINITIALIZED
    implicit none
    type(input_variables), intent(inout) :: in

    in%mdsteps=0
    in%md_printfrq = 1
    in%temperature = 300.d0
    in%dt = 20.d0
    in%no_translation = .false.
    in%nhc=.false.
    in%nhnc = 3
    in%nmultint = 1
    in%nsuzuki  = 7
    in%nosefrq  = 3000.d0
  END SUBROUTINE md_input_variables_default 

  !> Assign default values for self-interaction correction variables
  subroutine sic_input_variables_default(in)
    implicit none
    type(input_variables), intent(inout) :: in

    in%SIC%approach='NONE'
    in%SIC%alpha=0.0_gp
    in%SIC%fref=0.0_gp

  END SUBROUTINE sic_input_variables_default


  !> Assign default values for TDDFT variables
  subroutine tddft_input_variables_default(in)
    implicit none
    type(input_variables), intent(inout) :: in

    in%tddft_approach='NONE'

  END SUBROUTINE tddft_input_variables_default


  !>  Free all dynamically allocated memory from the input variable structure.
  subroutine free_input_variables(in)
    use dynamic_memory, only: f_free_ptr
    use multipole_base, only: deallocate_external_potential_descriptors
    implicit none
    type(input_variables), intent(inout) :: in
    character(len=*), parameter :: subname='free_input_variables'

    !check if freeing is possible
    call f_ref_free(in%refcnt)

    call free_geopt_variables(in)
    call free_kpt_variables(in)
    call f_free_ptr(in%gen_occup)
    call deallocateBasicArraysInput(in%lin)
    call deallocateInputFragArrays(in%frag)
    call deallocate_external_potential_descriptors(in%ep)

    ! Stop the signaling stuff.
    !Destroy C wrappers on Fortran objects,
    ! and stop the GMainLoop.
    if (in%gmainloop /= 0.d0) then
       call bigdft_signals_free(in%gmainloop)
    end if
  END SUBROUTINE free_input_variables


  !> Assign default values for ABSCALC variables
  subroutine abscalc_input_variables_default(in)
    implicit none
    type(input_variables), intent(inout) :: in

    in%c_absorbtion=.false.
    in%potshortcut=0
    in%iat_absorber=0
    in%abscalc_bottomshift=0
    in%abscalc_S_do_cg=.false.
    in%abscalc_Sinv_do_cg=.false.
  END SUBROUTINE abscalc_input_variables_default


  !> Assign default values for frequencies variables
  !!    freq_alpha: frequencies step for finite difference = alpha*hx, alpha*hy, alpha*hz
  !!    freq_order; order of the finite difference (2 or 3 i.e. 2 or 4 points)
  !!    freq_method: 1 - systematic moves of atoms over each direction
  subroutine frequencies_input_variables_default(in)
    implicit none
    type(input_variables), intent(inout) :: in

    in%freq_alpha=1.d0/real(64,kind(1.d0))
    in%freq_order=2
    in%freq_method=1
  END SUBROUTINE frequencies_input_variables_default

  subroutine allocate_extra_lin_arrays(lin,nspin,astruct)
    use module_atoms, only: atomic_structure
    use dynamic_memory
    use yaml_strings, only: yaml_toa
    implicit none
    integer,intent(in) :: nspin
    type(atomic_structure), intent(in) :: astruct
    type(linearInputParameters), intent(inout) :: lin
    !local variables
    character(len=*), parameter :: subname='allocate_extra_lin_arrays'
    integer :: nlr,iat,itype,iiorb,iorb,ispin
    !then perform extra allocations
    nlr=0
    do iat=1,astruct%nat
       itype=astruct%iatype(iat)
       nlr=nlr+nspin*lin%norbsPerType(itype)
    end do
  
    lin%locrad = f_malloc_ptr(nlr,id='lin%locrad')
    lin%locrad_kernel = f_malloc_ptr(nlr,id='lin%locrad_kernel')
    lin%locrad_mult = f_malloc_ptr(nlr,id='lin%locrad_mult')
    lin%locrad_lowaccuracy = f_malloc_ptr(nlr,id='lin%locrad_lowaccuracy')
    lin%locrad_highaccuracy = f_malloc_ptr(nlr,id='lin%locrad_highaccuracy')
  
    ! Assign the localization radius to each atom.
    iiorb=0
    do ispin=1,nspin
        do iat=1,astruct%nat
           itype=astruct%iatype(iat)
           do iorb=1,lin%norbsPerType(itype)
              iiorb=iiorb+1
              lin%locrad(iiorb)=lin%locrad_type(itype,1)
              lin%locrad_kernel(iiorb)=lin%kernel_cutoff(itype)
              lin%locrad_mult(iiorb)=lin%kernel_cutoff_FOE(itype)
              lin%locrad_lowaccuracy(iiorb)=lin%locrad_type(itype,1)
              lin%locrad_highaccuracy(iiorb)=lin%locrad_type(itype,2)
           end do
        end do
    end do
    if (iiorb/=nlr) then
        call f_err_throw('Error in filling the extra_lin_arrays, iiorb/=nlr ('&
        &//trim(yaml_toa(iiorb))//','//trim(yaml_toa(nlr))//')',err_name='BIGDFT_RUNTIME_ERROR')
    end if
  end subroutine allocate_extra_lin_arrays


  subroutine allocateBasicArraysInputLin(lin, ntypes)
    implicit none

    ! Calling arguments
    type(linearInputParameters), intent(inout) :: lin
    integer, intent(in) :: ntypes

    ! Local variables
    character(len=*), parameter :: subname='allocateBasicArrays'

    call f_routine(id='allocateBasicArraysInputLin')

    lin%norbsPerType = f_malloc_ptr(ntypes,id='lin%norbsPerType')
    lin%potentialPrefac_ao = f_malloc_ptr(ntypes,id='lin%potentialPrefac_ao')
    lin%potentialPrefac_lowaccuracy = f_malloc_ptr(ntypes,id='lin%potentialPrefac_lowaccuracy')
    lin%potentialPrefac_highaccuracy = f_malloc_ptr(ntypes,id='lin%potentialPrefac_highaccuracy')
    !added a second dimension to include the low and high accuracy values
    lin%locrad_type = f_malloc_ptr((/ ntypes, 2 /),id='lin%locrad_type')
    lin%kernel_cutoff_FOE = f_malloc_ptr(ntypes,id='lin%kernel_cutoff_FOE')
    lin%kernel_cutoff = f_malloc_ptr(ntypes,id='lin%kernel_cutoff')

    call f_release_routine()

  end subroutine allocateBasicArraysInputLin

  subroutine deallocateBasicArraysInput(lin)
    implicit none

    ! Calling arguments
    type(linearinputParameters), intent(inout) :: lin

    call f_routine(id='deallocateBasicArraysInput')

    call f_free_ptr(lin%potentialPrefac_ao)
    call f_free_ptr(lin%potentialPrefac_lowaccuracy)
    call f_free_ptr(lin%potentialPrefac_highaccuracy)
    call f_free_ptr(lin%norbsPerType)
    call f_free_ptr(lin%locrad)
    call f_free_ptr(lin%locrad_kernel)
    call f_free_ptr(lin%locrad_mult)
    call f_free_ptr(lin%locrad_lowaccuracy)
    call f_free_ptr(lin%locrad_highaccuracy)
    call f_free_ptr(lin%locrad_type)
    call f_free_ptr(lin%kernel_cutoff_FOE)
    call f_free_ptr(lin%kernel_cutoff)

    call f_release_routine()

  end subroutine deallocateBasicArraysInput

  !> Cross check values of input_variables.
  !! and change if necessary
  subroutine input_analyze(in,astruct)
    use module_atoms, only: atomic_structure
    use yaml_strings, only: operator(.eqv.)
    implicit none
    type(input_variables), intent(inout) :: in
    type(atomic_structure), intent(in) :: astruct

    integer :: ierr

    call f_routine(id='input_analyze')

    ! the PERF variables -----------------------------------------------------
    !Check after collecting all values
    if(.not.in%orthpar%directDiag .or. in%orthpar%methOrtho==1) then 
       write(*,'(1x,a)') 'Input Guess: Block size used for the orthonormalization (ig_blocks)'
       if(in%orthpar%bsLow==in%orthpar%bsUp) then
          write(*,'(5x,a,i0)') 'Take block size specified by user: ',in%orthpar%bsLow
       else if(in%orthpar%bsLow<in%orthpar%bsUp) then
          write(*,'(5x,2(a,i0))') 'Choose block size automatically between ',in%orthpar%bsLow,' and ',in%orthpar%bsUp
       else
          call f_err_throw("ERROR: invalid values of inputs%bsLow and inputs%bsUp. Change them in 'inputs.perf'!",&
               err_name='BIGDFT_INPUT_VARIABLES_ERROR')
          !call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
       end if
       write(*,'(5x,a)') 'This values will be adjusted if it is larger than the number of orbitals.'
    end if
    !@todo also the inputguess variable should be checked if BC are nonFree

    ! the DFT variables ------------------------------------------------------
    in%SIC%ixc = in%ixc

    in%idsx = min(in%idsx, in%itermax)

    !project however the wavefunction on gaussians if asking to write them on disk
    ! But not if we use linear scaling version (in%inputPsiId >= 100)
    !in%gaussian_help=(in%inputPsiId >= 10 .and. in%inputPsiId < 100)

!!$    !switch on the gaussian auxiliary treatment 
!!$    !and the zero of the forces
!!$    if (in%inputPsiId == 10) then
!!$       in%inputPsiId = 0
!!$    else if (in%inputPsiId == 13) then !better to insert gaussian_help as a input variable
!!$       in%inputPsiId = 2
!!$    end if

!!$    ! Setup out grid parameters.
!!$    if (in%output_denspot >= 0) then
!!$       in%output_denspot_format = in%output_denspot / 10
!!$    else
!!$       in%output_denspot_format = output_denspot_FORMAT_CUBE
!!$       in%output_denspot = abs(in%output_denspot)
!!$    end if
!!$    in%output_denspot = modulo(in%output_denspot, 10)

    !define whether there should be a last_run after geometry optimization
    !also the mulliken charge population should be inserted
!!$    if ((in%rbuf > 0.0_gp) .or. in%output_wf_format /= WF_FORMAT_NONE .or. &
!!$         in%output_denspot /= output_denspot_NONE .or. in%norbv /= 0) then
    if (in%rbuf > 0.0_gp .or. in%output_wf /= 'NONE' .or. &
         in%output_denspot /= 'NONE' .or. in%norbv /= 0) then
       in%last_run=-1 !last run to be done depending of the external conditions
    else
       in%last_run=0
    end if

    if (astruct%geocode == 'F' .or. astruct%nat == 0) then
       !Disable the symmetry
       in%disableSym = .true.
    end if

    ! the GEOPT variables ----------------------------------------------------
    !target stress tensor
    in%strtarget(:)=0.0_gp

    if (trim(in%geopt_approach) .eqv. "AB6MD") then
       if (in%ionmov /= 13) then
          in%nnos=0
          in%qmass = f_malloc_ptr(in%nnos, id = "in%qmass")
       end if
    end if

    ! Determine the SCF mode
    select case (in%lin%kernel_mode)
    case (KERNELMODE_DIRMIN)
       in%lin%scf_mode = LINEAR_DIRECT_MINIMIZATION
    case (KERNELMODE_DIAG)
       select case (in%lin%mixing_mode)
       case (MIXINGMODE_DENS)
          in%lin%scf_mode = LINEAR_MIXDENS_SIMPLE
       case (MIXINGMODE_POT)
          in%lin%scf_mode = LINEAR_MIXPOT_SIMPLE
       case default
          stop 'wrong value of in%lin%mixing_mode'
       end select
    case (KERNELMODE_FOE)
       in%lin%scf_mode = LINEAR_FOE
    case default
       call f_err_throw('wrong value of in%lin%kernel_mode',&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end select

    ! It is not possible to use both the old and the new Pulay correction at the same time
    if (in%lin%pulay_correction .and. in%lin%new_pulay_correction) then
       call f_err_throw('It is not possible to use both the old and the new Pulay correction at the same time!',&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

    call f_release_routine()
  END SUBROUTINE input_analyze


  !> Analyse the kpt input and calculates k points if needed
  subroutine kpt_input_analyse(iproc, in, dict, sym, geocode, alat)
    use module_base
    use module_atoms, only: symmetry_data
    use defs_basis
    use m_ab6_kpoints
    use yaml_output
    use public_keys
    implicit none
    !Arguments
    integer, intent(in) :: iproc
    type(input_variables), intent(inout) :: in
    type(dictionary), pointer, intent(in) :: dict
    type(symmetry_data), intent(in) :: sym
    character(len = 1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    real(gp), dimension(3), intent(in) :: alat
    !local variables
    logical :: lstat
    character(len=*), parameter :: subname='kpt_input_analyse'
    integer :: ierror,i, nshiftk, ikpt, j, ncount, nseg, iseg_, ngranularity_
    integer, dimension(3) :: ngkpt_
    real(gp), dimension(3) :: alat_
    real(gp), dimension(3,8) :: shiftk_
    real(gp) :: kptrlen_, norm
    character(len = 6) :: method
    real(gp), dimension(:,:), pointer :: gen_kpt   !< K points coordinates
    real(gp), dimension(:), pointer :: gen_wkpt    !< Weights of k points
    ! Set default values.
    in%gen_nkpt=1
    in%nkptv=0
    in%ngroups_kptv=1

    call free_kpt_variables(in)
    nullify(in%kptv, in%nkptsv_group)
    nullify(in%gen_kpt, in%gen_wkpt)

    method = dict // KPT_METHOD
    if (trim(method) .eqv. 'auto') then
       kptrlen_ = dict // KPTRLEN
       if (geocode == 'F') then
          in%gen_nkpt = 1
!!$        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$        in%gen_kpt = 0.
          in%gen_kpt=f_malloc0_ptr([3, in%gen_nkpt],id='gen_kpt')

!!$        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
          in%gen_kpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')

          in%gen_wkpt = 1.
       else
          call kpoints_get_auto_k_grid(sym%symObj, in%gen_nkpt, gen_kpt, gen_wkpt, &
               & kptrlen_, ierror)
          if (ierror /= AB7_NO_ERROR) then
             if (iproc==0) &
                  & call yaml_warning("ERROR: cannot generate automatic k-point grid." // &
                  & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
             stop
          end if
          !assumes that the allocation went through (arrays allocated by abinit routines)
          in%gen_kpt=f_malloc_ptr(src_ptr=gen_kpt,id='gen_kpt')
          in%gen_wkpt=f_malloc_ptr(src_ptr=gen_wkpt,id='gen_wkpt')
          deallocate(gen_kpt,gen_wkpt)
!!$        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
!!$        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
       end if
    else if (trim(method) .eqv. 'mpgrid') then
       !take the points of Monkhorst-pack grid
       ngkpt_(1:3) = dict // NGKPT
       if (geocode == 'S') ngkpt_(2) = 1
       !shift
       nshiftk = dict_len(dict//SHIFTK)
       !read the shifts
       shiftk_=0.0_gp
       do i=1,nshiftk
          shiftk_(1,i) = dict // SHIFTK // (i-1) // 0
          shiftk_(2,i) = dict // SHIFTK // (i-1) // 1
          shiftk_(3,i) = dict // SHIFTK // (i-1) // 2
       end do

       !control whether we are giving k-points to Free BC
       if (geocode == 'F') then
          if (iproc==0 .and. (maxval(ngkpt_) > 1 .or. maxval(abs(shiftk_)) > 0.)) &
               & call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
          in%gen_nkpt = 1
          in%gen_kpt=f_malloc0_ptr([3, in%gen_nkpt],id='gen_kpt')
!!$        allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$        in%gen_kpt = 0.
!!$        allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$        call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
          in%gen_kpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')
          in%gen_wkpt = 1.
       else
          call kpoints_get_mp_k_grid(sym%symObj, in%gen_nkpt, gen_kpt, gen_wkpt, &
               & ngkpt_, nshiftk, shiftk_, ierror)
          if (ierror /= AB7_NO_ERROR) then
             if (iproc==0) &
                  & call yaml_warning("ERROR: cannot generate MP k-point grid." // &
                  & " Error code is " // trim(yaml_toa(ierror,fmt='(i0)')))
             stop
          end if
          !assumes that the allocation went through 
          !(arrays allocated by abinit routines)
          in%gen_kpt=f_malloc_ptr(src_ptr=gen_kpt,id='gen_kpt')
          in%gen_wkpt=f_malloc_ptr(src_ptr=gen_wkpt,id='gen_wkpt')
          deallocate(gen_kpt,gen_wkpt)
!!$        call memocc(0,in%gen_kpt,'in%gen_kpt',subname)
!!$        call memocc(0,in%gen_wkpt,'in%gen_wkpt',subname)
       end if
    else if (trim(method) .eqv. 'manual') then
       in%gen_nkpt = max(1, dict_len(dict//KPT))
       if (geocode == 'F' .and. in%gen_nkpt > 1) then
          if (iproc==0) call yaml_warning('Found input k-points with Free Boundary Conditions, reduce run to Gamma point')
          in%gen_nkpt = 1
       end if
       in%gen_kpt=f_malloc_ptr([3, in%gen_nkpt],id='gen_kpt')
       in%gen_wkpt=f_malloc_ptr(in%gen_nkpt,id='gen_wkpt')

!!$     allocate(in%gen_kpt(3, in%gen_nkpt+ndebug),stat=i_stat)
!!$     call memocc(i_stat,in%gen_kpt,'in%gen_kpt',subname)
!!$     allocate(in%gen_wkpt(in%gen_nkpt+ndebug),stat=i_stat)
!!$     call memocc(i_stat,in%gen_wkpt,'in%gen_wkpt',subname)
       norm=0.0_gp
       do i=1,in%gen_nkpt
          in%gen_kpt(1, i) = dict // KPT // (i-1) // 0
          in%gen_kpt(2, i) = dict // KPT // (i-1) // 1
          in%gen_kpt(3, i) = dict // KPT // (i-1) // 2
          if (geocode == 'S' .and. in%gen_kpt(2,i) /= 0.) then
             in%gen_kpt(2,i) = 0.
             if (iproc==0) call yaml_warning('Surface conditions, suppressing k-points along y.')
          end if
          in%gen_wkpt(i) = dict // WKPT // (i-1)
          if (geocode == 'F') then
             in%gen_kpt = 0.
             in%gen_wkpt = 1.
          end if
          norm=norm+in%gen_wkpt(i)
       end do
       ! We normalise the weights.
       in%gen_wkpt(:)=in%gen_wkpt/norm
    else
       if (iproc==0) &
            & call yaml_warning("ERROR: wrong k-point sampling method (" // &
            & trim(method) // ").")
       stop
    end if

    ! Convert reduced coordinates into BZ coordinates.
    alat_ = alat
    if (geocode /= 'P') alat_(2) = 1.0_gp
    if (geocode == 'F') then
       alat_(1)=1.0_gp
       alat_(3)=1.0_gp
    end if
    do i = 1, in%gen_nkpt, 1
       in%gen_kpt(:, i) = in%gen_kpt(:, i) / alat_(:) * two_pi
    end do

    in%band_structure_filename=''
    lstat = dict // BANDS
    if (lstat) then
       !calculate the number of groups of for the band structure
       in%nkptv=1
       nseg = dict_len(dict // ISEG)
       do i=1,nseg
          iseg_ = dict // ISEG // (i-1)
          in%nkptv=in%nkptv+iseg_
       end do
       ngranularity_ = dict // NGRANULARITY

       in%ngroups_kptv=&
            ceiling(real(in%nkptv,gp)/real(ngranularity_,gp))

       in%nkptsv_group = f_malloc_ptr(in%ngroups_kptv,id='in%nkptsv_group')

       ncount=0
       do i=1,in%ngroups_kptv-1
          !if ngranularity is bigger than nkptv  then ngroups is one
          in%nkptsv_group(i)=ngranularity_
          ncount=ncount+ngranularity_
       end do
       !put the rest in the last group
       in%nkptsv_group(in%ngroups_kptv)=in%nkptv-ncount

       in%kptv = f_malloc_ptr((/ 3, in%nkptv /),id='in%kptv')

       ikpt = 0
       do i=1,nseg
          iseg_ = dict // ISEG // (i-1)
          ikpt=ikpt+iseg_
          in%kptv(1,ikpt) = dict // KPTV // (ikpt - 1) // 0
          in%kptv(2,ikpt) = dict // KPTV // (ikpt - 1) // 1
          in%kptv(3,ikpt) = dict // KPTV // (ikpt - 1) // 2
          !interpolate the values
          do j=ikpt-iseg_+1,ikpt-1
             in%kptv(:,j)=in%kptv(:,ikpt-iseg_) + &
                  (in%kptv(:,ikpt)-in%kptv(:,ikpt-iseg_)) * &
                  real(j-ikpt+iseg_,gp)/real(iseg_, gp)
          end do
       end do

       ! Convert reduced coordinates into BZ coordinates.
       do i = 1, in%nkptv, 1
          in%kptv(:, i) = in%kptv(:, i) / alat_(:) * two_pi
       end do

       if (has_key(dict, BAND_STRUCTURE_FILENAME)) then
          in%band_structure_filename = dict // BAND_STRUCTURE_FILENAME
          !since a file for the local potential is already given, do not perform ground state calculation
          if (iproc==0) then
             write(*,'(1x,a)')'Local Potential read from file, '//trim(in%band_structure_filename)//&
                  ', do not optimise GS wavefunctions'
          end if
          in%nrepmax=0
          in%itermax=0
          in%itrpmax=0
          !in%inputPsiId=-1000 !allocate empty wavefunctions
          call set_inputpsiid(INPUT_PSI_EMPTY,in%inputPsiId)
          call set_output_denspot(OUTPUT_DENSPOT_NONE,in%output_denspot)
          !in%output_denspot=0
       end if
    else
       in%nkptv = 0
       in%kptv = f_malloc_ptr((/ 3, in%nkptv /),id='in%kptv')
    end if

    if (in%nkptv > 0 .and. geocode == 'F' .and. iproc == 0) &
         & call yaml_warning('Defining a k-point path in free boundary conditions.') 

  END SUBROUTINE kpt_input_analyse

  integer function wave_format_from_filename(iproc, filename)
    use yaml_output
    implicit none
    integer, intent(in) :: iproc
    character(len=*), intent(in) :: filename

    integer :: isuffix

    wave_format_from_filename = WF_FORMAT_NONE

    isuffix = index(filename, ".etsf", back = .true.)
    if (isuffix > 0) then
       wave_format_from_filename = WF_FORMAT_ETSF
       if (iproc ==0) call yaml_comment('Reading wavefunctions in ETSF file format.')
       !if (iproc ==0) write(*,*) "Reading wavefunctions in ETSF file format."
    else
       isuffix = index(filename, ".bin", back = .true.)
       if (isuffix > 0) then
          wave_format_from_filename = WF_FORMAT_BINARY
          if (iproc ==0) call yaml_comment('Reading wavefunctions in BigDFT binary file format.')
          !if (iproc ==0) write(*,*) "Reading wavefunctions in BigDFT binary file format."
       else
          wave_format_from_filename = WF_FORMAT_PLAIN
          if (iproc ==0) call yaml_comment('Reading wavefunctions in plain text file format.')
          !if (iproc ==0) write(*,*) "Reading wavefunctions in plain text file format."
       end if
    end if
  end function wave_format_from_filename

  !> Print all general parameters
  subroutine print_general_parameters(in,atoms,input_id,posinp_id)
    use module_atoms, only: atoms_data
    use defs_basis
    use yaml_output
    use yaml_strings, only: operator(.eqv.),yaml_toa
    implicit none
    !Arguments
    type(input_variables), intent(in) :: in
    type(atoms_data), intent(in) :: atoms
    character(len = *), intent(in) :: input_id, posinp_id

    integer :: iat, i
    character(len = 11) :: potden
    character(len = 12) :: dos

    ! Output for atoms
    call yaml_comment('Input Atomic System (file: '//trim(posinp_id)//'.'//trim(atoms%astruct%inputfile_format)//')',hfill='-')

    ! Atomic systems
    call yaml_mapping_open('Atomic System Properties')
    call yaml_map('Number of atomic types', atoms%astruct%ntypes, fmt='(i0)')
    call yaml_map('Number of atoms', atoms%astruct%nat, fmt='(i0)')
    if (atoms%astruct%nat > 0) then
       call yaml_map('Types of atoms',atoms%astruct%atomnames)
       ! Fixed positions
       if (maxval(atoms%astruct%ifrztyp) /= 0) then
          call yaml_sequence_open('Fixed atoms',flow=.true.)
          ! The fixed atom column
          do iat=1,atoms%astruct%nat
             if (atoms%astruct%ifrztyp(iat) /= 0) then
                call yaml_sequence('at.' // trim(yaml_toa(iat,fmt='(i4.4)')) // &
                     & '(' // trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))) // ')' &
                     & // trim(yaml_toa(atoms%astruct%ifrztyp(iat),fmt='(i0)')))
             end if
          end do
          call yaml_sequence_close()
       end if
    end if
    !Boundary Conditions
    select case(atoms%astruct%geocode)
    case('P')
       call yaml_map('Boundary Conditions','Periodic',advance='no')
       call yaml_comment('Code: '//atoms%astruct%geocode)
       call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
            atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
    case('S')
       call yaml_map('Boundary Conditions','Surface',advance='no')
       call yaml_comment('Code: '//atoms%astruct%geocode)
       call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
            atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
    case('W')
       call yaml_map('Boundary Conditions','Wire',advance='no')
       call yaml_comment('Code: '//atoms%astruct%geocode)
       call yaml_map('Box Sizes (AU)',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
            atoms%astruct%cell_dim(3)/),fmt='(1pe12.5)')
    case('F')
       call yaml_map('Boundary Conditions','Free',advance='no')
       call yaml_comment('Code: '//atoms%astruct%geocode)
    end select

    !Symmetries
    call yaml_map('Number of Symmetries',atoms%astruct%sym%nSym)
    call yaml_map('Space group',trim(atoms%astruct%sym%spaceGroup))
    call yaml_mapping_close()

    !Geometry imput Parameters
    if (in%ncount_cluster_x > 0) then
       call yaml_comment('Geometry optimization Input Parameters (file: '//trim(input_id)//'.geopt)',hfill='-')
       call yaml_mapping_open('Geometry Optimization Parameters')
       call yaml_map('Maximum steps',in%ncount_cluster_x)
       call yaml_map('Algorithm', in%geopt_approach)
       call yaml_map('Random atomic displacement', in%randdis, fmt='(1pe7.1)')
       call yaml_map('Fluctuation in forces',in%frac_fluct,fmt='(1pe7.1)')
       call yaml_map('Maximum in forces',in%forcemax,fmt='(1pe7.1)')
       call yaml_map('Steepest descent step',in%betax,fmt='(1pe7.1)')
       if (trim(in%geopt_approach) == "DIIS") then
          call yaml_map('DIIS history', in%history)
       end if
       call yaml_mapping_close()
       if (trim(in%geopt_approach) .eqv. "AB6MD") then
          call yaml_mapping_open('Molecular Dynamics Parameters')
          call yaml_map('ionmov',in%ionmov)
          call yaml_map('dtion', in%dtion,fmt='(0pf7.0)')
          if (in%ionmov > 7) then
             call yaml_map('Start Temperature', in%mditemp, fmt='(f5.0)')
             call yaml_map('Stop Temperature',  in%mdftemp, fmt='(f5.0)')
          end if
          if (in%ionmov == 8) then
             call yaml_map('Nose inertia', in%noseinert,fmt='(f15.5)')
          else if (in%ionmov == 9) then
             call yaml_map('Friction', in%friction,fmt='(f15.5)')
             call yaml_map('MD wall',in%mdwall,fmt='(f15.5)')
          else if (in%ionmov == 13) then
             call yaml_map('nnos', in%nnos,fmt='(f15.5)')
             call yaml_map('qmass',in%qmass,fmt='(f15.5)')
             call yaml_map('bmass',in%bmass,fmt='(f15.5)')
             call yaml_map('vmass',in%vmass,fmt='(f15.5)')
          end if
          call yaml_mapping_close()
       end if
    end if
    !MD input
    if (in%mdsteps > 0) then
       call yaml_comment('Molecular Dynamics Input Parameters',hfill='-')
       call yaml_mapping_open('Molecular Dynamics Parameters')
       call yaml_map('Maximum MD steps',in%mdsteps)
       call yaml_map('Printing Frequency', in%md_printfrq)
       call yaml_map('Initial Temperature (K)', in%temperature, fmt='(1pe7.1)')
       call yaml_map('Time step (a.u.)',in%dt,fmt='(1pe7.1)')
       call yaml_map('Freeze Translation ', in%no_translation)
       call yaml_map('Nose Hoover Chain Thermostat', in%nhc)
       if(in%nhc)then
         call yaml_map('Length of Nose Hoover Chains', in%nhnc)
         call yaml_map('Multiple Time Step for Nose Hoover Chains', in%nmultint)
         call yaml_map('Yoshida-Suzuki factor for Nose Hoover Chains', in%nsuzuki)
         call yaml_map('Frequency of Nose Hoover Chains', in%nosefrq)
       end if
       call yaml_mapping_close()
    end if

    !Output for K points
    if (atoms%astruct%geocode /= 'F') then
       call yaml_comment('K points description (Reduced and Brillouin zone coordinates, Weight)',hfill='-')
       !write(*,'(1x,a)') '--- (file: input.kpt) ----------------------------------------------------- k-points'
       if (in%disableSym .and. in%gen_nkpt > 1) then
          call yaml_warning('symmetries have been disabled, k points are not irreductible.')
          !write(*, "(1x,A)") "WARNING: symmetries have been disabled, k points are not irreductible."
       end if
       call yaml_sequence_open('K points')!,advance='no')
       !call yaml_comment('Reduced coordinates  BZ coordinates  weight',hfill=' ')
       !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
       do i = 1, in%gen_nkpt, 1
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map( 'Rc', &
               & in%gen_kpt(:, i) * atoms%astruct%cell_dim / two_pi,&
               & fmt='(f7.4)')
          call yaml_map( 'Bz', &
               & in%gen_kpt(:, i), &
               & fmt='(f7.4)')
          call yaml_map('Wgt',in%gen_wkpt(i),fmt='(f6.4)')
          call yaml_mapping_close(advance='no')
          call yaml_comment(trim(yaml_toa(i,fmt='(i4.4)')))
          !write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
          !     & in%kpt(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
          !     & in%wkpt(i), i, in%kpt(:, i)
       end do
       call yaml_sequence_close()

       if (in%nkptv > 0) then
          call yaml_sequence_open('K points for band structure calculation')
          !write(*, "(1x,a)")    " K points for band structure calculation"
          !write(*, "(1x,a)")    "       red. coordinates         weight       id        BZ coordinates"
          do i = 1, in%nkptv, 1
             call yaml_sequence(advance='no')
             call yaml_mapping_open(trim(yaml_toa(i,fmt='(i0)')),flow=.true.)
             call yaml_map( 'Red C.', &
                  & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), &
                  & atoms%astruct%cell_dim(3) /) / two_pi,&
                  & fmt='(f9.5)')
             call yaml_map( 'Bz C.', &
                  & in%kptv(:, i), &
                  & fmt='(f9.5)')
             call yaml_map('Weight',1.0d0 / real(size(in%kptv, 2), gp),fmt='(f9.5)')
             call yaml_mapping_close()
             !   write(*, "(1x,3f9.5,2x,f9.5,5x,I4,1x,3f9.5)") &
             !        & in%kptv(:, i) * (/ atoms%astruct%cell_dim(1), atoms%astruct%cell_dim(2), atoms%astruct%cell_dim(3) /) / two_pi, &
             !        & 1.0d0 / real(size(in%kptv, 2), gp), i, in%kptv(:, i)
          end do
          call yaml_sequence_close()
       end if
    end if

    ! Printing for mixing parameters.
    if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) then
       if (in%iscf < 10) then
          write(potden, "(A)") "potential"
       else
          write(potden, "(A)") "density"
       end if
       call yaml_comment('Mixing (file: '//trim(input_id)//'.mix)',hfill='-')
       call yaml_mapping_open('Mixing parameters')
       call yaml_map('Target',trim(potden))
       call yaml_map('Additional bands', in%norbsempty)
       call yaml_map('Mixing Coefficients', in%alphamix,fmt='(0pe10.2)')
       call yaml_map('Scheme',modulo(in%iscf, 10))
       call yaml_map('Electronic temperature',in%Tel,fmt='(1pe12.2)')
       call yaml_map('DIIS',in%alphadiis,fmt='(0pe12.2)')
       call yaml_map('Maximum iterations',in%itrpmax)
       call yaml_map('Occupied scheme',trim(smearing_names(in%occopt)))
       call yaml_map('Rp norm',in%rpnrm_cv,fmt='(1pe12.2)')
       if (in%verbosity > 2) then
          write(dos, "(A)") "dos.gnuplot"
       else
          write(dos, "(A)") "no verb. < 3"
       end if
       call yaml_map('output DOS',trim(dos))
       call yaml_mapping_close()
    end if

  END SUBROUTINE print_general_parameters


  !> Print all dft input parameters
  subroutine print_dft_parameters(in,atoms)
    use module_atoms, only: atoms_data
    use yaml_output
    use module_xc
    use f_enums, only: f_int => int
    use yaml_strings, only: yaml_toa
    implicit none
    type(input_variables), intent(in) :: in
    type(atoms_data), intent(in) :: atoms

    call yaml_comment('Input parameters',hfill='-')

    call yaml_mapping_open('DFT parameters')
    call yaml_mapping_open('eXchange Correlation')
    call yaml_map('XC ID',in%ixc,fmt='(i8)',label='ixc')
    if (in%ixc < 0) then
       call xc_dump(in%ixc, XC_MIXED, in%nspin)
    else ! @todo@ if (in%ixc /= XC_NO_HARTREE) then
       call xc_dump(in%ixc, XC_ABINIT, in%nspin)
    end if
    if (in%nspin>=2) then
       call yaml_map('Polarisation',in%mpol)
    end if
    if (in%nspin==4) then
       call yaml_map('Spin polarization','non collinear')
    else if (in%nspin==2) then
       call yaml_map('Spin polarization','collinear')
    else if (in%nspin==1) then
       call yaml_map('Spin polarization',.false.)
    end if
    call yaml_mapping_close()

    if (in%qcharge /= 0.0_gp) call yaml_map('Net Charge (Ions-Electrons)',in%qcharge,fmt='(f8.5)')!'(i8)')
    if (sqrt(sum(in%elecfield(:)**2)) > 0.0_gp) &
         call yaml_map('External Electric Field (Ha/a0)',in%elecfield(:),fmt='(1pe8.1)')
    call yaml_mapping_close()

    call yaml_mapping_open('Basis set definition')
    call yaml_map('Suggested Grid Spacings (a0)', (/in%hx,in%hy,in%hz/),fmt='(f5.2)')
    call yaml_map('Coarse and Fine Radii Multipliers', (/in%crmult,in%frmult/),fmt='(f4.1)')
    call yaml_mapping_close()


    call yaml_mapping_open('Self-Consistent Cycle Parameters')
    call yaml_mapping_open('Wavefunction')
    call yaml_map('Gradient Norm Threshold',in%gnrm_cv,fmt='(1pe8.1)',label='gnrm_cv')
    call yaml_map('CG Steps for Preconditioner',in%ncong,fmt='(i5)')
    call yaml_map('DIIS History length',in%idsx)
    call yaml_map('Max. Wfn Iterations',in%itermax,label='itermax')
    call yaml_map('Max. Subspace Diagonalizations',in%nrepmax)
!!$    call yaml_map('Input wavefunction policy',  trim(input_psi_names(in%inputPsiId)), advance="no")
    call yaml_map('Input wavefunction policy',  trim(str(in%inputPsiId)), advance="no")
    call yaml_comment(trim(yaml_toa(f_int(in%inputPsiId))))
    !call yaml_comment(trim(yaml_toa(in%inputPsiId)))
    call yaml_map('Output wavefunction policy', trim(str(in%output_wf)), advance="no")
!!$    call yaml_map('Output wavefunction policy', trim(wf_format_names(in%output_wf_format)), advance="no")
!!$    call yaml_comment(trim(yaml_toa(in%output_wf_format)))
    call yaml_comment(trim(yaml_toa(f_int(in%output_wf))))
    call yaml_map('Output grid policy',trim(str(in%output_denspot)),advance='no')
!!$    call yaml_map('Output grid policy',trim(output_denspot_names(in%output_denspot)),advance='no')
!!$    call yaml_comment(trim(yaml_toa(in%output_denspot)))
    call yaml_comment(trim(yaml_toa(f_int(in%output_denspot))))
!!$        call yaml_map('Output grid format',trim(output_denspot_format_names(in%output_denspot_format)),advance='no')
!!$        call yaml_comment(trim(yaml_toa(in%output_denspot_format)))
    if (in%output_denspot .hasattr. 'TEXT') then
       call yaml_map('Output grid format','TEXT',advance='no')
    else if (in%output_denspot .hasattr. 'BINARY') then
       call yaml_map('Output grid format','BINARY',advance='no')
    else if (in%output_denspot .hasattr. 'CUBE') then
       call yaml_map('Output grid format','CUBE',advance='no')
    else if (in%output_denspot .hasattr. 'ETSF') then
       call yaml_map('Output grid format','ETSF',advance='no')
    end if

    call yaml_map('Virtual orbitals',in%nvirt,fmt='(i0)')
    call yaml_map('Number of plotted density orbitals',abs(in%nplot),fmt='(i0)')

    call yaml_mapping_close()

    call yaml_mapping_open('Density/Potential')
    call yaml_map('Max. Iterations',in%itrpmax)
    call yaml_mapping_close()
    call yaml_mapping_close()

    if (atoms%astruct%geocode == 'F') then
       call yaml_mapping_open('Post Optimization Parameters')

       call yaml_mapping_open('Finite-Size Effect estimation')
       call yaml_map('Scheduled',(in%rbuf > 0.0_gp))
       if (in%rbuf > 0.0_gp) then
          call yaml_map('Extension',in%rbuf,fmt='(f4.1)')
          call yaml_map('No. of CG steps',in%ncongt)
       end if
       call yaml_mapping_close()
       call yaml_mapping_close()
    end if

  END SUBROUTINE print_dft_parameters

  !> Read from all input files and build a dictionary
  recursive subroutine user_dict_from_files(dict,radical,posinp_name, mpi_env)
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use wrapper_MPI, only: mpi_environment
    use public_keys, only: POSINP, IG_OCCUPATION, MODE_VARIABLES, SECTIONS, METHOD_KEY
    use yaml_output
    use yaml_strings, only: f_strcpy
    use f_utils, only: f_file_exists
    use module_input_dicts
    !use input_old_text_format
    use module_atoms, only: astruct_file_merge_to_dict,atoms_file_merge_to_dict
    implicit none
    !Arguments
    type(dictionary), pointer :: dict                  !< Contains (out) all the information
    character(len = *), intent(in) :: radical          !< Radical for the input files
    character(len = *), intent(in) :: posinp_name           !< If the dict has no posinp key, use it
    type(mpi_environment), intent(in) :: mpi_env       !< MPI Environment
    !Local variables
    logical :: exists
    type(dictionary), pointer :: at, iter
    character(len = max_field_length) :: str, rad

    !read the input file(s) and transform them into a dictionary
    call read_input_dict_from_files(trim(radical), mpi_env, dict)

    !possible overwrite with a specific posinp file.
    if (len_trim(posinp_name) > 0) then
       call astruct_file_merge_to_dict(dict,POSINP, trim(posinp_name))
    end if
    if (has_key(dict,POSINP)) then
       str = dict_value(dict //POSINP)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          !str contains a file name so add atomic positions from it.
          call astruct_file_merge_to_dict(dict,POSINP, trim(str))
       else
          !The yaml file contains the atomic positions
          !Only add the format
          at => dict //POSINP
          if (.not. has_key(at, ASTRUCT_PROPERTIES)) then
             call set(at // ASTRUCT_PROPERTIES // FORMAT_KEY, FORMAT_YAML)
          else
             at => at // ASTRUCT_PROPERTIES
             if (FORMAT_KEY .notin. at) &
                  call set(at // FORMAT_KEY, FORMAT_YAML)
          end if
       end if
    end if

    ! Add old psppar
    call atoms_file_merge_to_dict(dict)

    call f_strcpy(src = radical, dest = rad)
    if (len_trim(radical) == 0) rad = "input"

    !when the user has not specified the occupation in the input file
    if (.not. has_key(dict,IG_OCCUPATION)) then
       !yaml format should be used even for old method
       call f_file_exists(trim(rad)//".occup",exists)
       if (exists) &
            call merge_input_file_to_dict(dict//IG_OCCUPATION,&
            trim(rad)//".occup",mpi_env)
    else !otherwise the input file always supersedes
       str = dict_value(dict //IG_OCCUPATION)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          !call atomic_data_file_merge_to_dict(dict, ATOMIC_OCC, trim(str))
          call f_file_exists(trim(str),exists)
          if (exists) &
               call merge_input_file_to_dict(dict//IG_OCCUPATION,trim(str),mpi_env)
       end if
    end if

    if (OCCUPATION .notin. dict) then
       ! Add old input.occ
       call occupation_data_file_merge_to_dict(dict,OCCUPATION,trim(rad) // ".occ")
    else
       str = dict_value(dict //OCCUPATION)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          call occupation_data_file_merge_to_dict(dict,OCCUPATION, trim(str))
       end if
    end if

    ! Add section files, if any.
    if (has_key(dict, MODE_VARIABLES)) then
       str = dict_value(dict // MODE_VARIABLES // METHOD_KEY)
       if (trim(str) == 'multi' .and. has_key(dict // MODE_VARIABLES, SECTIONS)) then
          iter => dict_iter(dict // MODE_VARIABLES // SECTIONS)
          do while (associated(iter))
             str = dict_value(dict // dict_value(iter))
             if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
                if (len_trim(str) > 5 .and. str(max(1,len_trim(str)-4):len_trim(str)) == ".yaml") then
                   call user_dict_from_files(dict // dict_value(iter), str(1:len_trim(str)-5), "", mpi_env)
                else
                   call user_dict_from_files(dict // dict_value(iter), str, "", mpi_env)
                end if
             end if
             iter => dict_next(iter)
          end do
       end if
    end if
  end subroutine user_dict_from_files



!!!  !> Write input parameters
!!!  subroutine write_input_parameters(in)!,atoms)
!!!    use module_base
!!!    use module_types
!!!    use yaml_output
!!!    implicit none
!!!    type(input_variables), intent(in) :: in
!!!    !  type(atoms_data), intent(in) :: atoms
!!!    !local variables
!!!    character(len = 11) :: potden
!!!    !start yaml output
!!!    !call yaml_indent_map('Physical System Parameters')
!!!    !write(70,'(a)')repeat(' ',yaml_indent)//'Physical System Parameters:'
!!!    !yaml_indent=yaml_indent+3
!!!    !  write(70,'(a,t55,a)')repeat(' ',yaml_indent)//'Boundary Conditions:',atoms%astruct%geocode
!!!    !  if (atoms%astruct%geocode /= 'F')write(70,'(a,t55,a,3(1x,f5.3,a))')&
!!!    !       repeat(' ',yaml_indent)//'Box Sizes (a0):','[',atoms%astruct%cell_dim(1),',',atoms%astruct%cell_dim(2),',',atoms%astruct%cell_dim(3),' ]'
!!!
!!!    !  yaml_indent=yaml_indent-3
!!!
!!!
!!!    if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) then
!!!       !write(70,'(a)')repeat(' ',yaml_indent)//'Mixing Parameters:'
!!!       !yaml_indent=yaml_indent+3
!!!       if (in%iscf < 10) then
!!!          write(potden, "(A)") "potential"
!!!       else
!!!          write(potden, "(A)") "density"
!!!       end if
!!!       write(70,'(a,t55,a)')'Target:',potden
!!!       !     write(70,'(a,t55,I12)')'Scheme:',modulo(in%iscf, 10)
!!!!!$     write(*,"(1x,A12,A12,1x,A1,1x,A12,I12,1x,A1,1x,A11,F10.2)") &
!!!!!$          & "     Target=", potden,        "|", &
!!!!!$          & " Add. bands=", in%norbsempty, "|", &
!!!!!$          & "    Coeff.=", in%alphamix
!!!!!$     write(*,"(1x,A12,I12,1x,A1,1x,A12,1pe12.2,1x,A1,1x,A11,0pe10.2)") &
!!!!!$          & "     Scheme=", modulo(in%iscf, 10), "|", &
!!!!!$          & "Elec. temp.=", in%Tel,              "|", &
!!!!!$          & "      DIIS=", in%alphadiis
!!!!!$     write(*,"(1x,A12,I12,1x,A1,1x,A12,A12,1x,A1)") &
!!!!!$          & "  Max iter.=", in%itrpmax,    "|", &
!!!!!$          & "Occ. scheme=", smearing_names(in%occopt), "|"
!!!!!$     if (in%verbosity > 2) then
!!!!!$        write(dos, "(A)") "dos.gnuplot"
!!!!!$     else
!!!!!$        write(dos, "(A)") "no verb. < 3"
!!!!!$     end if
!!!!!$     write(*,"(1x,A12,1pe12.2,1x,A1,1x,2A12,1x,A1)") &
!!!!!$          & "   Rp norm.=", in%rpnrm_cv,    "|", " output DOS=", dos, "|"
!!!    end if
!!!    !write(70,'(a)')repeat(' ',yaml_indent)//'Post Optimization Treatments:'
!!!    if (in%rbuf > 0.0_gp) then
!!!       !write(70,'(a)')repeat(' ',yaml_indent)//'Finite-Size Correction Estimation:'
!!!       !write(70,'(a,t55,f4.1)')repeat(' ',yaml_indent)//'Radius (a0):',in%rbuf
!!!       !write(70,'(a,t55,i4)')repeat(' ',yaml_indent)//'CG Steps for the FS Correction:',in%ncongt
!!!    end if
!!!    stop
!!!  end subroutine write_input_parameters
!!!

end module module_input_keys

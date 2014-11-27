!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_global_variables
    use module_base !bigdft base module
    use module_types
    use bigdft_run
    implicit none
    character(len = *), public, parameter :: mhgps_version   = '0.01'
    character(len = *), public, parameter :: inputdir   = 'input'
    character(len = *), public, parameter :: outputdir   = 'input'

    !input parameters for mhgps
    integer          :: mhgps_verbosity     = 3
    character(len=10) :: efmethod            = "LJ" !method/force-field
                                                   !for energies,
                                                   !stresses and
                                                   !forces
    logical          :: external_mini       = .false.
    character(len=5) :: input_dir           = "input"

    !input parameters for sqnm saddle_search
    character(len=20):: operation_mode             = 'connect'
    logical          :: saddle_biomode             = .false.
    logical          :: random_minmode_guess       = .true.
    integer          :: imode                      = 1
    integer          :: saddle_nit_trans           = 5000
    integer          :: saddle_nit_rot             = 5000
    integer          :: saddle_nhistx_trans        = 40
    integer          :: saddle_nhistx_rot          = 4
    real(gp)         :: saddle_alpha0_trans        = 1.e-3_gp
    real(gp)         :: saddle_alpha0_rot          = 1.e-3_gp
    real(gp)         :: saddle_alpha_stretch0      = 4.e-4_gp
    real(gp)         :: saddle_alpha_rot_stretch0  = 2.e-4_gp
    real(gp)         :: saddle_curvgraddiff        = 1.e-3_gp
    real(gp)         :: saddle_rmsdispl0           = 0.04_gp
    real(gp)         :: saddle_trustr              = 0.2_gp
    real(gp)         :: saddle_tolc                = 7.0_gp
    real(gp)         :: saddle_tolf                = 7.0_gp
    logical          :: saddle_tighten             = .true.
!    real(gp)         :: saddle_tightenfac          = -1.0_gp
    real(gp)         :: saddle_maxcurvrise         = 1.e-6_gp
    real(gp)         :: saddle_cutoffratio         = 1.e-4_gp
    real(gp)         :: saddle_minoverlap0         = 0.95_gp
    real(gp)         :: saddle_steepthresh_trans   = 1._gp
    real(gp)         :: saddle_steepthresh_rot     = 1._gp
    integer          :: saddle_recompIfCurvPos     = 5
    real(gp)         :: saddle_fnrmtol             = 1.e-3_gp
    real(gp)         :: saddle_stepoff             = 2.e-2_gp
    real(gp)         :: saddle_scale_stepoff       = 2.e+0_gp
    logical          :: share_rot_history = .false.!not available via
                                                   !input file since
                                                   !sharing tends to
                                                   !introduce a slight
                                                   !inefficiency when
                                                   !compared to no
                                                   !sharing (roughly
                                                   !10% for LJ75)

    !parameters for minimizers implemented in mhgps
    logical :: internal=.true. !use internal or external optimizers?
        !SQNM
        integer  :: mini_nhistx = 15
        integer  :: mini_ncluster_x = 5000
        real(gp) :: mini_frac_fluct = 0.0_gp
        real(gp) :: mini_forcemax = 1.e-4_gp
        real(gp) :: mini_maxrise = 1.e-6_gp
        real(gp) :: mini_betax = 1.0_gp
        real(gp) :: mini_beta_stretchx = 4.e-1_gp
        real(gp) :: mini_cutoffRatio = 1.e-4_gp
        real(gp) :: mini_steepthresh = 1.0_gp
        real(gp) :: mini_trustr = 0.5_gp

    !parameters for inputguess
    !ts guess parameters
    real(gp) :: ts_guess_gammainv=1._gp
    real(gp) :: ts_guess_perpnrmtol=1.e-3_gp
    real(gp) :: ts_guess_trust=0.05_gp
    integer  :: ts_guess_nstepsmax=5
    real(gp) :: lst_interpol_stepfrct=0.1_gp
    real(gp) :: lst_dt_max=5.0_gp
    real(gp) :: lst_fmax_tol=5.e-3_gp

    !variables for connect routine
    integer :: nsadmax=30

    !others
    integer               :: nbond = 1
    integer, allocatable  :: iconnect(:,:) 
    integer, allocatable  :: ixyz_int(:,:)
    real(gp), allocatable :: minmode(:,:)
    integer,parameter     :: usaddle=173
    character(len=60)     :: saddle_filename='saddle.mon'
    logical               :: isForceField=.false.
    real(gp)              :: ef_counter=0.d0
    character(len=8)      :: currDir
    character(len=3), parameter :: prefix='pos'
    character(len=5)      :: isadc
    integer               :: isad
    character(len=5)      :: isadprobc
    integer               :: isadprob=0
    real(gp)              :: en_delta_min, fp_delta_min
    real(gp)              :: en_delta_sad, fp_delta_sad

    !bigdft data types and variables 
    !(these objects must preserve their status in the module)
    type(run_objects), save :: runObj
    type(DFT_global_output), save :: outs
    integer, save :: fdim
    !type(atoms_data), save :: atoms
    type(atomic_structure), pointer, save :: astruct_ptr
    !integer, dimension(4) :: mpi_info
    integer :: infocode
    integer :: inputPsiId=0
    integer :: iproc=0,nproc=1,igroup=0,ngroups=1
    integer :: itermin=0
    real(gp) :: frac_fluct=0.d0


    !following variables might be packed in an own module...
    integer               :: lwork
    real(gp), allocatable :: work(:)
    !variables for rotation
    integer :: nhist_rot,ndim_rot
    real(gp), allocatable :: rxyz_rot(:,:,:)
    real(gp), allocatable :: fxyz_rot(:,:,:)
    real(gp), allocatable :: fxyzraw_rot(:,:,:)
    real(gp), allocatable :: rxyzraw_rot(:,:,:)
    real(gp), allocatable :: fstretch_rot(:,:,:)
    real(gp), allocatable :: eval_rot(:)
    real(gp), allocatable :: res_rot(:)
    real(gp), allocatable :: rrr_rot(:,:,:)
    real(gp), allocatable :: aa_rot(:,:)
    real(gp), allocatable :: ff_rot(:,:,:)
    real(gp), allocatable :: rr_rot(:,:,:)
    real(gp), allocatable :: dd_rot(:,:)
    real(gp), allocatable :: fff_rot(:,:,:)
    real(gp), allocatable :: scpr_rot(:)
    real(gp), allocatable :: wold_rot(:)
    real(gp) :: alpha_rot, alpha_stretch_rot
    !translation
    real(gp), allocatable :: rxyz_trans(:,:,:)
    real(gp), allocatable :: fxyz_trans(:,:,:)
    real(gp), allocatable :: fxyzraw_trans(:,:,:)
    real(gp), allocatable :: rxyzraw_trans(:,:,:)
    real(gp), allocatable :: fstretch_trans(:,:,:)
    real(gp), allocatable :: eval_trans(:)
    real(gp), allocatable :: res_trans(:)
    real(gp), allocatable :: rrr_trans(:,:,:)
    real(gp), allocatable :: aa_trans(:,:)
    real(gp), allocatable :: ff_trans(:,:,:)
    real(gp), allocatable :: rr_trans(:,:,:)
    real(gp), allocatable :: dd_trans(:,:)
    real(gp), allocatable :: fff_trans(:,:,:)
    real(gp), allocatable :: scpr_trans(:)
    real(gp), allocatable :: wold_trans(:)

end module

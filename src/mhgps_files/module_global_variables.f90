!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module module_global_variables
    use module_base, only: gp !bigdft base module
    !use module_types
!    use module_atoms, only: atomic_structure
!    use bigdft_run, only: run_objects, state_properties
    implicit none
    character(len = *), public, parameter :: mhgps_version   = '0.01'
    character(len = *), public, parameter :: inputdir   = 'input'
    character(len = *), public, parameter :: outputdir   = 'input'

    !input parameters for mhgps
    integer, save          :: mhgps_verbosity     = 3
!    character(len=10) :: efmethod            = "LJ" !method/force-field
!                                                   !for energies,
!                                                   !stresses and
                                                   !forces
    logical, save          :: external_mini       = .false.
    character(len=5), save :: input_dir           = "input"

    !input parameters for sqnm saddle_search
    character(len=20), save :: operation_mode             = 'connect'
    logical, save          :: saddle_biomode             = .false.
    logical, save          :: random_minmode_guess       = .true.
    integer, save          :: imode                      = 1
    integer, save          :: saddle_nit_trans           = 5000
    integer, save          :: saddle_nit_rot             = 5000
    integer, save          :: saddle_nhistx_trans        = 40
    integer, save          :: saddle_nhistx_rot          = 4
    real(gp), save         :: saddle_alpha0_trans        = 1.e-3_gp
    real(gp), save         :: saddle_alpha0_rot          = 1.e-3_gp
    real(gp), save         :: saddle_alpha_stretch0      = 4.e-4_gp
    real(gp), save         :: saddle_alpha_rot_stretch0  = 2.e-4_gp
    real(gp), save         :: saddle_curvgraddiff        = 1.e-3_gp
    real(gp), save         :: saddle_rmsdispl0           = 0.04_gp
    real(gp), save         :: saddle_trustr              = 0.2_gp
    real(gp), save         :: saddle_tolc                = 7.0_gp
    real(gp), save         :: saddle_tolf                = 7.0_gp
    logical, save          :: saddle_tighten             = .true.
    real(gp), save         :: saddle_maxcurvrise         = 1.e-6_gp
    real(gp), save         :: saddle_cutoffratio         = 1.e-4_gp
    real(gp), save         :: saddle_minoverlap0         = 0.95_gp
    real(gp), save         :: saddle_steepthresh_trans   = 1._gp
    real(gp), save         :: saddle_steepthresh_rot     = 1._gp
    integer , save         :: saddle_recompIfCurvPos     = 5
    real(gp), save         :: saddle_fnrmtol             = 1.e-3_gp
    real(gp), save         :: saddle_stepoff             = 2.e-2_gp
    real(gp), save         :: saddle_scale_stepoff       = 2.e+0_gp
    logical , save         :: share_rot_history = .false.!not available via
                                                   !input file since
                                                   !sharing tends to
                                                   !introduce a slight
                                                   !inefficiency when
                                                   !compared to no
                                                   !sharing (roughly
                                                   !10% for LJ75)

    !parameters for minimizers implemented in mhgps
    logical, save :: internal=.true. !use internal or external optimizers?
        !SQNM
        integer, save  :: mini_nhistx = 15
        integer, save  :: mini_ncluster_x = 5000
        real(gp), save :: mini_frac_fluct = 0.0_gp
        real(gp), save :: mini_forcemax = 1.e-4_gp
        real(gp), save :: mini_maxrise = 1.e-6_gp
        real(gp), save :: mini_betax = 1.0_gp
        real(gp), save :: mini_beta_stretchx = 4.e-1_gp
        real(gp), save :: mini_cutoffRatio = 1.e-4_gp
        real(gp), save :: mini_steepthresh = 1.0_gp
        real(gp), save :: mini_trustr = 0.5_gp

    !parameters for inputguess
    !ts guess parameters
    real(gp), save :: ts_guess_gammainv=1._gp
    real(gp), save :: ts_guess_perpnrmtol=1.e-3_gp
    real(gp), save :: ts_guess_trust=0.05_gp
    integer, save  :: ts_guess_nstepsmax=5
    real(gp), save :: lst_interpol_stepfrct=0.1_gp
    real(gp), save :: lst_dt_max=5.0_gp
    real(gp), save :: lst_fmax_tol=5.e-3_gp

    !variables for connect routine
    integer, save :: nsadmax=30

    !others
    integer, save               :: nbond = 1
    integer, save, allocatable  :: iconnect(:,:) 
    real(gp), save, allocatable :: minmode(:,:)
    character(len=60), save     :: saddle_filename='saddle.mon'
    logical, save               :: isForceField=.false.
    real(gp), save              :: ef_counter=0.d0
    character(len=8), save      :: currDir
    character(len=3), parameter :: prefix='pos'
    character(len=5), save      :: isadc
    integer, save               :: isad
    integer, save               :: ntodo
    character(len=5), save      :: isadprobc
    integer , save              :: isadprob=0
    real(gp), save              :: en_delta_min, fp_delta_min
    real(gp), save              :: en_delta_sad, fp_delta_sad

    !bigdft data types and variables 
    !(these objects must preserve their status in the module)
!    type(run_objects), save :: runObj
!    type(state_properties), save :: outs
    integer, save :: fdim
    !type(atoms_data), save :: atoms
!    type(atomic_structure), pointer, save :: astruct_ptr
    !integer, dimension(4) :: mpi_info
    integer, save :: infocode
    integer, save :: inputPsiId=0
    integer, save :: iproc=0,nproc=1,igroup=0,ngroups=1
    integer, save :: itermin=0
    real(gp), save :: frac_fluct=0.d0


    !following variables might be packed in an own module...
    integer, save               :: lwork
    real(gp), save, allocatable :: work(:)
    !variables for rotation
    integer, save :: nhist_rot,ndim_rot
    real(gp), save, allocatable :: rxyz_rot(:,:,:)
    real(gp), save, allocatable :: fxyz_rot(:,:,:)
    real(gp), save, allocatable :: fxyzraw_rot(:,:,:)
    real(gp), save, allocatable :: rxyzraw_rot(:,:,:)
    real(gp), save, allocatable :: fstretch_rot(:,:,:)
    real(gp), save, allocatable :: eval_rot(:)
    real(gp), save, allocatable :: res_rot(:)
    real(gp), save, allocatable :: rrr_rot(:,:,:)
    real(gp), save, allocatable :: aa_rot(:,:)
    real(gp), save, allocatable :: ff_rot(:,:,:)
    real(gp), save, allocatable :: rr_rot(:,:,:)
    real(gp), save, allocatable :: dd_rot(:,:)
    real(gp), save, allocatable :: fff_rot(:,:,:)
    real(gp), save, allocatable :: scpr_rot(:)
    real(gp), save, allocatable :: wold_rot(:)
    real(gp), save :: alpha_rot, alpha_stretch_rot
    !translation
    real(gp), save, allocatable :: rxyz_trans(:,:,:)
    real(gp), save, allocatable :: fxyz_trans(:,:,:)
    real(gp), save, allocatable :: fxyzraw_trans(:,:,:)
    real(gp), save, allocatable :: rxyzraw_trans(:,:,:)
    real(gp), save, allocatable :: fstretch_trans(:,:,:)
    real(gp), save, allocatable :: eval_trans(:)
    real(gp), save, allocatable :: res_trans(:)
    real(gp), save, allocatable :: rrr_trans(:,:,:)
    real(gp), save, allocatable :: aa_trans(:,:)
    real(gp), save, allocatable :: ff_trans(:,:,:)
    real(gp), save, allocatable :: rr_trans(:,:,:)
    real(gp), save, allocatable :: dd_trans(:,:)
    real(gp), save, allocatable :: fff_trans(:,:,:)
    real(gp), save, allocatable :: scpr_trans(:)
    real(gp), save, allocatable :: wold_trans(:)

end module

module module_global_variables
    use module_base !bigdft base module
    use module_types
    implicit none
    character(len = *), public, parameter :: mhgps_version   = '0.01'
    character(len=6), public, parameter   :: efmethod_known(2)   = (/"LJ","BIGDFT"/)

    !input parameters for mhgps
    integer          :: mhgps_verbosity     = 3
    character(len=6) :: efmethod            = "LJ" !method/force-field for energies, stresses and forces
    character(len=5) :: input_dir           = "input"

    !input parameters for saddle_search
    logical          :: saddle_biomode             = .false.
    logical          :: saddle_connect             = .false.
    logical          :: random_minmode_guess       = .true.
    integer          :: saddle_imode               = 1
    integer          :: saddle_nit_trans           = 5000
    integer          :: saddle_nit_rot             = 5000
    integer          :: saddle_nhistx_trans        = 40
    integer          :: saddle_nhistx_rot          = 4
    real(gp)         :: saddle_alpha0_trans        = 1.e-3_gp
    real(gp)         :: saddle_alpha0_rot          = 1.e-3_gp
    real(gp)         :: saddle_alpha_stretch0      = 4.e-4_gp
    real(gp)         :: saddle_curvgraddiff        = 1.e-3_gp
    real(gp)         :: saddle_rmsdispl0           = 0.04_gp
    real(gp)         :: saddle_trustr              = 0.2_gp
    real(gp)         :: saddle_tolc                = 7.0_gp
    real(gp)         :: saddle_tolf                = 7.0_gp
    real(gp)         :: saddle_tightenfac          = -1.0_gp
    real(gp)         :: saddle_maxcurvrise         = 1.e-6_gp
    real(gp)         :: saddle_cutoffratio         = 1.e-4_gp
    integer          :: saddle_recompIfCurvPos     = 5
    real(gp)         :: saddle_fnrmtol             = 1.e-3_gp

    !parameters for minimizers implemented in mhgps
    logical :: internal=.true.
        !SBFGS
        integer  :: mini_nhistx
        integer  :: mini_nit
        real(gp) :: mini_maxrise
        real(gp) :: mini_betax
        real(gp) :: mini_cutoffRatio
        real(gp) :: mini_steepthresh
        real(gp) :: mini_trustr

    !other
    integer          :: nbond               = 1
    integer, allocatable :: iconnect(:,:) 
    integer, allocatable :: ixyz_int(:,:)
    real(gp), allocatable :: minmode(:,:)
    integer,parameter :: usaddle=173
    character(len=60) :: saddle_filename='saddle.mon'
    logical :: isForceField=.false.
    real(gp) :: ef_counter=0.d0


    !bigdft data types and variables
    type(run_objects) :: runObj
    type(dictionary), pointer :: user_inputs
    type(DFT_global_output) :: outs
    type(atoms_data) :: atoms
    integer, dimension(4) :: mpi_info
    integer :: infocode
    type(input_variables), target :: inputs_opt
    type(restart_objects) :: rst
    integer :: inputPsiId=0
    integer :: iproc=0,nproc=1,igroup=0,ngroups=1
    integer :: itermin=0
end module

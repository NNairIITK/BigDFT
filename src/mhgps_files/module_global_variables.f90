module module_global_variables
    use module_base !bigdft base module
    use module_types
    implicit none

    type, public :: globals
        !input parameters for mhgps
        integer :: mhgps_verbosity=3
        character(len=6) :: efmethod_known(2)=(/"LJ","BIGDFT"/)
        character(len=6) :: efmethod="LJ" !method/force-field for energies, stresses and forces
        logical :: biomode=.false.
        integer :: nit_trans=5000
        integer :: nit_rot=5000
        integer :: nhistx_trans=40
        integer :: nhistx_rot=4
        real(gp) :: alpha0_trans=1.d-3
        real(gp) :: alpha0_rot=1.d-3
        real(gp) :: alpha_stretch0=4.d-4
        real(gp) :: curvgraddiff=1.d-3
        real(gp) :: rmsdispl0=0.04d0
        real(gp) :: trustr=0.2d0
        real(gp) :: tolc=7.d0
        real(gp) :: tolf=7.d0
        real(gp) :: tightenfac=-1.d0
        real(gp) :: maxcurvrise=1.d-6
        real(gp) :: cutoffratio=1.d-4
        integer  :: recompIfCurvPos=5

        !bigdft data types and variables
        type(run_objects) :: runObj
        type(dictionary), pointer :: user_inputs
        type(DFT_global_output) :: outs
        integer, dimension(4) :: mpi_info
        integer :: infocode
    end type
end module

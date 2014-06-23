module module_global_variables
    use module_base !bigdft base module
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

        !bigdft data types
        type(run_objects) :: runObj
        type(dictionary), pointer :: user_inputs
        type(DFT_global_output) :: outs
    end type
contains

    subroutine init_global_variables(glob)
        implicit none
        !parameter
        type(globals), intent(out) :: glob
        call read_input(glob)
        if(glob%efmethod=='BIGDFT')then
        elseif(glov%efmethod=='LJ')then
!        TODO
        endif
    end subroutine

    subroutine read_input(glob)
        implicit none
        !parameter
        type(globals), intent(out) :: glob
        !internal
        integer, parameter :: u=237
        character(9), parameter :: filename='mhgps.inp'
        logical :: exists

        inquire(file=filename,exist=exists)
        if(.not. exists)then
            call write_input(glob)
            stop 'ERROR: mhgps.inp does not exist.&
                  Wrote default input parameters to mhgps.inp_default.'
        endif
        open(u,file=filename)
            read(u,*)glob%mhgps_verbosity
            read(u,*)glob%efmethod
            read(u,*)glob%biomode
            read(u,*)glob%nit_trans, glob%nit_rot
            read(u,*)glob%nhistx_trans, glob%nhistx_rot
            if(glob%biomode)then
                read(u,*)glob%alpha0_trans, glob%alpha0_rot, glob%alpha_stretch0
            else
                read(u,*)glob%alpha0_trans, glob%alpha0_rot
            endif
            read(u,*)glob%curvgraddiff
            read(u,*)glob%rmsdispl0,glob%trustr
            read(u,*)glob%tolc,glob%tolf
            read(u,*)glob%tightenfac
            read(u,*)glob%maxcurvrise
            read(u,*)glob%cutoffratio
            read(u,*)glob%recompIfCurvPos
        close(u)
    end subroutine
    subroutine write_input(glob)
        implicit none
        !parameter
        type(globals), intent(in) :: glob
        integer, parameter :: u=237
        character(17), parameter :: filename='mhgps.inp_default'
        open(u,file=filename)
            write(u,'(xi0.0,xa)')glob%mhgps_verbosity,' #mhgps_verbosity'
            write(u,'(xa,xa)')trim(adjustl(glob%efmethod)),' #efmethod'
            write(u,'(xL,xa)')glob%biomode,' #biomode'
            write(u,'(xi0,xi0,xa)')glob%nit_trans, glob%nit_rot,'  #nit_trans, not_rot'
            write(u,'(xi0,xi0,xa)')glob%nhistx_trans, glob%nhistx_rot,' #nhistx_trans, nhistx_rot'
            if(glob%biomode)then
                write(u,'(es10.3,2(xes10.3),a)')glob%alpha0_trans, glob%alpha0_rot, glob%alpha_stretch0,' #alpha0_trans, alpha0_rot, alpha_stretch0'
            else
                write(u,'(es10.3,xes10.3,a)')glob%alpha0_trans, glob%alpha0_rot,' #alpha0_trans, alpha0_rot'
            endif
            write(u,'(es10.3,xa)')glob%curvgraddiff,' #curvgraddif'
            write(u,'(es10.3,xes10.3,xa)')glob%rmsdispl0,glob%trustr,' #rmsdispl0, trustr'
            write(u,'(es10.3,xes10.3,xa)')glob%tolc,glob%tolf,' #tolc, tolf'
            write(u,'(es10.3,xa)')glob%tightenfac,' #tightenfac'
            write(u,'(es10.3,xa)')glob%maxcurvrise,' #maxcurvrise'
            write(u,'(es10.3,xa)')glob%cutoffratio,' #cutoffratio'
            write(u,'(xi0,xa)')glob%recompIfCurvPos,' #recompIfCurvPos'
        close(u)

    end subroutine

end module

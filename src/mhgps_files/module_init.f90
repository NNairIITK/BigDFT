module module_init
    use module_base !bigdft base module
    use module_global_variables
    implicit none

contains

    subroutine init_global_variables(glob)
        implicit none
        !parameter
        type(globals), intent(out) :: glob
        !internal
        integer:: nconfig,ierr,run_id
        call print_logo_mhgps()
        call read_input(glob)
        if(glob%efmethod=='BIGDFT')then
            call bigdft_init(glob%mpi_info,nconfig,run_id,ierr)
            if (nconfig < 0) stop 'runs-file not supported for MHGPS executable'
            call print_logo()
        elseif(glob%efmethod=='LJ')then
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

    subroutine print_logo_mhgps()
        use yaml_output
        implicit none

        call yaml_comment('Minima Hopping Guided Path Sampling....',hfill='=')
        
        call yaml_open_map('MHGPS logo')
        call yaml_scalar('      ___           ___           ___           ___           ___     ') 
        call yaml_scalar('     /\__\         /\__\         /\  \         /\  \         /\  \    ')
        call yaml_scalar('    /::|  |       /:/  /        /::\  \       /::\  \       /::\  \   ')
        call yaml_scalar('   /:|:|  |      /:/__/        /:/\:\  \     /:/\:\  \     /:/\ \  \  ')
        call yaml_scalar('  /:/|:|__|__   /::\  \ ___   /:/  \:\  \   /::\~\:\  \   _\:\~\ \  \ ')
        call yaml_scalar(' /:/ |::::\__\ /:/\:\  /\__\ /:/__/_\:\__\ /:/\:\ \:\__\ /\ \:\ \ \__\')
        call yaml_scalar(' \/__/~~/:/  / \/__\:\/:/  / \:\  /\ \/__/ \/__\:\/:/  / \:\ \:\ \/__/')
        call yaml_scalar('       /:/  /       \::/  /   \:\ \:\__\        \::/  /   \:\ \:\__\  ')
        call yaml_scalar('      /:/  /        /:/  /     \:\/:/  /         \/__/     \:\/:/  /  ')
        call yaml_scalar('     /:/  /        /:/  /       \::/  /                     \::/  /   ')
        call yaml_scalar('     \/__/         \/__/         \/__/                       \/__/    ')
        call yaml_scalar('                                                   as post-processing ')
        call yaml_scalar('')
        call yaml_scalar('')
        !call print_logo()
        call yaml_close_map()
        call yaml_map('Reference Paper','The Journal of Chemical Physics 140 (21):214102 (2014)')
    end subroutine print_logo_mhgps


end module

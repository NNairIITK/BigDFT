!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module module_init
    !use module_base !bigdft base module
    use module_global_variables
    implicit none

    private

    public :: read_input
    public :: write_input
    public :: print_input
    public :: print_logo_mhgps
    public :: give_rcov

contains

    subroutine read_input()
        use yaml_output
        use yaml_parse
        use dictionaries
        implicit none
        !parameter
        !internal
        integer, parameter :: u=237
        character(9), parameter :: filename='mhgps.inp'
        logical :: exists
!!! EXAMPLES BY LUIGI FOR READING YAML FILE
!!!        integer :: i
!!!        character(len=32) :: key
!!!        type(dictionary), pointer :: dinp !< the input file
!!!        type(dictionary), pointer :: dtmp !< the parsed input file
!!!        type(dictionary), pointer :: parameters,it !< the needed variables
!!!
!!!
!!!
!!!        call f_file_exists(trim(filename),exists)
!!!        if (exists) then
!!!           call yaml_parse_from_file(dtmp,trim(filename))
!!!           dinp => dtmp .pop. 0
!!!           call dict_free(dtmp)
!!!        else
!!!           call f_err_throw('Input file name "'//trim(filename)//&
!!!                '" does not exists, exiting...')
!!!           return
!!!        end if
!!!        
!!!        !first, dump the values which have been parsed as a dictionary
!!!        call yaml_map('Parsed input file',dinp)
!!!
!!!        !set the wanted parameters
!!!        parameters=>list_new( .item. 'value1',.item. 'value2',.item. 'value3')
!!!
!!!        !check that all the keys of parameters exist (example of sanity check of input file)
!!!        it => dict_iter(parameters)
!!!        do while(associated(it))
!!!           print *,'iterator value',trim(dict_value(it))
!!!           !check for the existence
!!!           print *, trim(dict_value(it)) .in. dinp
!!!           !increase iterator
!!!           it => dict_next(it)
!!!        end do
!!!
!!!        !fortran approach
!!!        do i=0,dict_len(parameters)-1
!!!           key=parameters//i
!!!           call yaml_map('Searching for'//trim(key),trim(key) .in. dinp)
!!!        end do
!!!
!!!        !example of reading of variables
!!!        !retrieve the value of key "value1" and "value2"
!!!        !an put them in fortran variables
!!!        !for example, mhgps_verbosity and operation_mode
!!!        mhgps_verbosity=dinp//'value1'
!!!        operation_mode=dinp//'value2'
!!!        mhgps_verbosity=dinp//'value2'
!!!        !test
!!!        print *,'first case',mhgps_verbosity,operation_mode
!!!        !however, they might be not present in the input file.
!!!        !let us check for the existence of the key
!!!        print *,'value1' .in. dinp,'value2' .in. dinp
!!!        
!!!        stop



        inquire(file=filename,exist=exists)
        if(.not. exists)then
            call write_input()
            call yaml_warning('mhgps.inp does not exist.  Wrote '//&
                              'default input parameters to '//&
                              'mhgps.inp_default.')
             stop
        endif
        open(u,file=filename)
            read(u,*)mhgps_verbosity
            read(u,*)operation_mode, random_minmode_guess
            read(u,*)nsadmax
            read(u,*)external_mini
            read(u,*)en_delta_min,fp_delta_min
            read(u,*)en_delta_sad,fp_delta_sad
            read(u,*)saddle_biomode
            if(saddle_biomode)imode=2
            read(u,*) lst_interpol_stepfrct
            read(u,*) ts_guess_gammainv
            read(u,*) ts_guess_perpnrmtol
            read(u,*) ts_guess_trust
            read(u,*) ts_guess_nstepsmax
            read(u,*) lst_dt_max, lst_fmax_tol
            read(u,*)saddle_nit_trans, saddle_nit_rot
            read(u,*)saddle_nhistx_trans, saddle_nhistx_rot
            read(u,*)saddle_steepthresh_trans,saddle_steepthresh_rot
            read(u,*)saddle_fnrmtol
            if(saddle_biomode)then
                read(u,*)saddle_alpha0_trans, saddle_alpha0_rot,&
                     saddle_alpha_stretch0, saddle_alpha_rot_stretch0
            else
                read(u,*)saddle_alpha0_trans, saddle_alpha0_rot
            endif
            read(u,*)saddle_curvgraddiff
            read(u,*)saddle_rmsdispl0,saddle_trustr
            read(u,*)saddle_tolc,saddle_tolf,saddle_tighten
            read(u,*)saddle_minoverlap0
            read(u,*)saddle_maxcurvrise
            read(u,*)saddle_cutoffratio
            read(u,*)saddle_recompIfCurvPos
            read(u,*)saddle_stepoff, saddle_scale_stepoff
            if(.not.external_mini)then
                read(u,*)mini_nhistx
                read(u,*)mini_ncluster_x
                read(u,*)mini_frac_fluct
                read(u,*)mini_forcemax
                read(u,*)mini_maxrise
                read(u,*)mini_betax
                read(u,*)mini_beta_stretchx
                read(u,*)mini_cutoffRatio
                read(u,*)mini_steepthresh
                read(u,*)mini_trustr
            endif

        close(u)
    end subroutine
    subroutine write_input()
        implicit none
        !parameter
        integer, parameter :: u=237
        character(17), parameter :: filename='mhgps.inp_default'
        open(u,file=filename)
            write(u,'(1x,i0.0,1x,1a)')mhgps_verbosity,' #mhgps_verbosity'
            write(u,'(1x,1L1,1x,1L1,1x,1a)')trim(adjustl(operation_mode)),random_minmode_guess,&
                 ' #mode, random_minmode_guess'
            write(u,'(1x,i0.0,1x,1a)')nsadmax,' #nsadmax'
            write(u,'(1x,1L1,1x,1a)')external_mini,&
                           ' #external minimizer'
            write(u,'(es10.3,1x,es10.3,1x,1a)')en_delta_min,fp_delta_min,' #en_delta_min, fp_delta_min'
            write(u,'(es10.3,1x,es10.3,1x,1a)')en_delta_sad,fp_delta_sad,' #en_delta_sad, fp_delta_sad'
            write(u,'(1x,1L1,1x,1a)')saddle_biomode,' #biomode'
            write(u,'(1x,es10.3,1x,1a)') lst_interpol_stepfrct,&
                 ' #inward interpolation distance as fraction of initial distance'
            write(u,'(1x,es10.3,1x,1a)') ts_guess_gammainv,' #step size for perpedicular optimization in freezing string method'
            write(u,'(1x,es10.3,1x,1a)') ts_guess_perpnrmtol,&
                 ' #convergence criterion perpedicular force in freezing string method'//&
                 ' (disable perpend. optim. by setting this value to a negative number)'
            write(u,'(1x,es10.3,1x,1a)') ts_guess_trust,' #trust radius freezing string method (maximum change of any coordinate'
            write(u,'(1x,i0,1x,1a)') ts_guess_nstepsmax,&
                 ' #maximum number of steps in perpendicular optimization in freezing stringmethod'
            write(u,'(1x,es10.3,1x,es10.3,1x,1a)')lst_dt_max, lst_fmax_tol,&
                 '#max. time step in fire optimizer of lst function, convergence criterion'
            write(u,'(1x,i0,1x,i0,1x,1a)')saddle_nit_trans, saddle_nit_rot,'  #nit_trans, nit_rot'
            write(u,'(1x,i0,1x,i0,1x,1a)')saddle_nhistx_trans, saddle_nhistx_rot,' #nhistx_trans, nhistx_rot'
            write(u,'(es10.3,1x,es10.3,1a)')saddle_steepthresh_trans,saddle_steepthresh_rot,&
                 ' #saddle_steepthresh_trans,saddle_steepthresh_rot'
            write(u,'(es10.3,1x,1a)')saddle_fnrmtol,' #fnrm tolerence convergence criterion for saddle point'
            if(saddle_biomode)then
                write(u,'(es10.3,3(1x,es10.3),1a)')saddle_alpha0_trans, saddle_alpha0_rot,&
                     saddle_alpha_stretch0, saddle_alpha_rot_stretch0,&
                     ' #alpha0_trans, alpha0_rot, alpha_stretch0, alpha_rot_stretch0'
            else
                write(u,'(es10.3,1x,es10.3,1a)')saddle_alpha0_trans, saddle_alpha0_rot,' #alpha0_trans, alpha0_rot'
            endif
            write(u,'(es10.3,1x,1a)')saddle_curvgraddiff,' #curvgraddif'
            write(u,'(es10.3,1x,es10.3,1x,1a)')saddle_rmsdispl0,saddle_trustr,' #rmsdispl0, trustr'
            write(u,'(es10.3,1x,es10.3,1x,1L1,1x,1a)')saddle_tolc,saddle_tolf,saddle_tighten,' #tolc, tolf, tighten'
            write(u,'(es10.3,1x,1a)')saddle_minoverlap0,' #minoverlap0'
            write(u,'(es10.3,1x,1a)')saddle_maxcurvrise,' #maxcurvrise'
            write(u,'(es10.3,1x,1a)')saddle_cutoffratio,' #cutoffratio'
            write(u,'(1x,i0,1x,1a)')saddle_recompIfCurvPos,' #recompIfCurvPos'
            write(u,'(es10.3,1x,es10.3,1x,1a)')saddle_stepoff,saddle_scale_stepoff, ' #stepoff, stepoff_scale'
            write(u,'(1x,i0,1x,1a)')mini_nhistx,'#mini_nhistx'
            write(u,'(1x,i0,1x,1a)')mini_ncluster_x,'#mini_ncluster_x'
            write(u,'(es10.3,1x,1a)')mini_frac_fluct,'#mini_frac_fluct'
            write(u,'(es10.3,1x,1a)')mini_forcemax,'#mini_forcemax'
            write(u,'(es10.3,1x,1a)')mini_maxrise,'#mini_maxrise'
            write(u,'(es10.3,1x,1a)')mini_betax,'#mini_betax'
            write(u,'(es10.3,1x,1a)')mini_beta_stretchx,'#mini_beta_stretchx'
            write(u,'(es10.3,1x,1a)')mini_cutoffRatio,'#mini_cutoffRatio'
            write(u,'(es10.3,1x,1a)')mini_steepthresh,'#mini_steepthresh'
            write(u,'(es10.3,1x,1a)')mini_trustr,'#mini_trustr'
        close(u)
    end subroutine
    subroutine print_input()
        use module_global_variables, only: mhgps_version
        use yaml_output
        call yaml_comment('(MHGPS) Input Parameters',hfill='-')
        call yaml_map('(MHGPS) mhgps_verbosity',mhgps_verbosity)
        call yaml_map('(MHGPS) operation_mode',trim(adjustl(operation_mode)))
        call yaml_map('(MHGPS) random_minmode_guess',random_minmode_guess)
        call yaml_map('(MHGPS) nsadmax',nsadmax)
        call yaml_map('(MHGPS) external minimizer',external_mini)
        call yaml_map('(MHGPS) en_delta_min',en_delta_min)
        call yaml_map('(MHGPS) en_delta_sad',en_delta_sad)
        call yaml_map('(MHGPS) Biomolecule mode',saddle_biomode)
        call yaml_map('(MHGPS) lst_interpol_stepfrct',lst_interpol_stepfrct)
        call yaml_map('(MHGPS) ts_guess_gammainv',ts_guess_gammainv)
        call yaml_map('(MHGPS) ts_guess_perpnrmtol',ts_guess_perpnrmtol)
        call yaml_map('(MHGPS) ts_guess_trust',ts_guess_trust)
        call yaml_map('(MHGPS) ts_guess_nstepsmax',ts_guess_nstepsmax)
        call yaml_map('(MHGPS) lst_dt_max',lst_dt_max)
        call yaml_map('(MHGPS) lst_fmax_tol',lst_fmax_tol)
        call yaml_map('(MHGPS) saddle_nit_trans',saddle_nit_trans)
        call yaml_map('(MHGPS) saddle_nit_rot',saddle_nit_rot)
        call yaml_map('(MHGPS) saddle_nhistx_trans',saddle_nhistx_trans)
        call yaml_map('(MHGPS) saddle_nhistx_rot',saddle_nhistx_rot)
        call yaml_map('(MHGPS) saddle_steepthresh_trans',saddle_steepthresh_trans)
        call yaml_map('(MHGPS) saddle_steepthresh_rot',saddle_steepthresh_rot)
        call yaml_map('(MHGPS) saddle_fnrmtol',saddle_fnrmtol)
        call yaml_map('(MHGPS) saddle_alpha0_trans',saddle_alpha0_trans)
        call yaml_map('(MHGPS) saddle_alpha0_rot',saddle_alpha0_rot)
        if(saddle_biomode)then
            call yaml_map('(MHGPS) saddle_alpha_stretch0',saddle_alpha_stretch0)
            call yaml_map('(MHGPS) saddle_alpha_rot_stretch0',saddle_alpha_rot_stretch0)
        endif
        call yaml_map('(MHGPS) saddle_curvgraddiff',saddle_curvgraddiff)
        call yaml_map('(MHGPS) saddle_rmsdispl0',saddle_rmsdispl0)
        call yaml_map('(MHGPS) saddle_trustr',saddle_trustr)
        call yaml_map('(MHGPS) saddle_tolc',saddle_tolc)
        call yaml_map('(MHGPS) saddle_tolf',saddle_tolf)
        call yaml_map('(MHGPS) saddle_tighten',saddle_tighten)
        call yaml_map('(MHGPS) saddle_maxcurvrise',saddle_maxcurvrise)
        call yaml_map('(MHGPS) saddle_cutoffratio',saddle_cutoffratio)
        call yaml_map('(MHGPS) saddle_recompIfCurvPos',saddle_recompIfCurvPos)
        call yaml_map('(MHGPS) saddle_stepoff',saddle_stepoff)
        call yaml_map('(MHGPS) saddle_scale_stepoff',saddle_scale_stepoff)
        if(.not.external_mini)then
            call yaml_map('(MHGPS) mini_nhistx', mini_nhistx)
            call yaml_map('(MHGPS) mini_ncluster_x', mini_ncluster_x)
            call yaml_map('(MHGPS) mini_frac_fluct', mini_frac_fluct)
            call yaml_map('(MHGPS) mini_forcemax', mini_forcemax)
            call yaml_map('(MHGPS) mini_maxrise', mini_maxrise)
            call yaml_map('(MHGPS) mini_betax', mini_betax)
            call yaml_map('(MHGPS) mini_beta_stretchx', mini_beta_stretchx)
            call yaml_map('(MHGPS) mini_cutoffRatio', mini_cutoffRatio)
            call yaml_map('(MHGPS) mini_steepthresh', mini_steepthresh)
            call yaml_map('(MHGPS) mini_trustr', mini_trustr)
        endif
    end subroutine print_input

    subroutine print_logo_mhgps()
        use module_global_variables, only: mhgps_version
        use yaml_output
        implicit none

        call yaml_comment('(MHGPS) Minima Hopping Guided Path Sampling',hfill='=')
        
        call yaml_mapping_open('(MHGPS) logo')
        call yaml_scalar('(MHGPS)      ___           ___           ___           ___           ___      ') 
        call yaml_scalar('(MHGPS)     /\__\         /\__\         /\  \         /\  \         /\  \     ')
        call yaml_scalar('(MHGPS)    /::|  |       /:/  /        /::\  \       /::\  \       /::\  \    ')
        call yaml_scalar('(MHGPS)   /:|:|  |      /:/__/        /:/\:\  \     /:/\:\  \     /:/\ \  \   ')
        call yaml_scalar('(MHGPS)  /:/|:|__|__   /::\  \ ___   /:/  \:\  \   /::\~\:\  \   _\:\~\ \  \  ')
        call yaml_scalar('(MHGPS) /:/ |::::\__\ /:/\:\  /\__\ /:/__/_\:\__\ /:/\:\ \:\__\ /\ \:\ \ \__\ ')
        call yaml_scalar('(MHGPS) \/__/~~/:/  / \/__\:\/:/  / \:\  /\ \/__/ \/__\:\/:/  / \:\ \:\ \/__/ ')
        call yaml_scalar('(MHGPS)       /:/  /       \::/  /   \:\ \:\__\        \::/  /   \:\ \:\__\   ')
        call yaml_scalar('(MHGPS)      /:/  /        /:/  /     \:\/:/  /         \/__/     \:\/:/  /   ')
        call yaml_scalar('(MHGPS)     /:/  /        /:/  /       \::/  /                     \::/  /    ')
        call yaml_scalar('(MHGPS)     \/__/         \/__/         \/__/                       \/__/     ')
        call yaml_scalar('(MHGPS)                                                   as post-processing  ')
        call yaml_scalar('(MHGPS)')
        call yaml_scalar('(MHGPS)')
        !call print_logo()
        call yaml_mapping_close()
        call yaml_mapping_open('(MHGPS) Reference Paper')
        call yaml_scalar('(MHGPS) The Journal of Chemical Physics 140 (21):214102 (2014)')
        call yaml_mapping_close()
        call yaml_map('(MHGPS) Version Number',mhgps_version)
        call yaml_map('(MHGPS) Timestamp of this run',yaml_date_and_time_toa())
    end subroutine print_logo_mhgps


subroutine give_rcov(astruct,nat,rcov)
  use module_base, only: gp
  use module_types, only: atomic_structure
  use yaml_output
  use module_global_variables, only: iproc
  use dictionaries, only: f_err_throw
  implicit none
  !Arguments
  integer, intent(in) :: nat
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(out) :: rcov(nat)
  !Local variables
  integer :: iat

  do iat=1,nat
     select case(trim(astruct%atomnames(astruct%iatype(iat))))
     case('H') 
        rcov(iat)=0.75d0
!        rcov(iat)=0.75d0*0.529177211d0
     case('LJ') 
        rcov(iat)=0.56d0
     case('He') 
        rcov(iat)=0.75d0
     case('Li') 
        rcov(iat)=3.40d0
     case('Be') 
        rcov(iat)=2.30d0
     case('B' ) 
        rcov(iat)=1.55d0
     case('C' ) 
        rcov(iat)=1.45d0
!        rcov(iat)=1.45d0*0.529177211d0
     case('N' ) 
        rcov(iat)=1.42d0
!        rcov(iat)=1.42d0*0.529177211d0
     case('O' ) 
        rcov(iat)=1.38d0
!        rcov(iat)=1.38d0*0.529177211d0
     case('F' ) 
        rcov(iat)=1.35d0
     case('Ne') 
        rcov(iat)=1.35d0
     case('Na') 
        rcov(iat)=3.40d0
     case('Mg') 
        rcov(iat)=2.65d0
     case('Al') 
        rcov(iat)=2.23d0
     case('Si') 
        rcov(iat)=2.09d0
     case('P' ) 
        rcov(iat)=2.00d0
     case('S' ) 
        rcov(iat)=1.92d0
     case('Cl') 
        rcov(iat)=1.87d0
     case('Ar') 
        rcov(iat)=1.80d0
     case('K' ) 
        rcov(iat)=4.00d0
     case('Ca') 
        rcov(iat)=3.00d0
     case('Sc') 
        rcov(iat)=2.70d0
     case('Ti') 
        rcov(iat)=2.70d0
     case('V' ) 
        rcov(iat)=2.60d0
     case('Cr') 
        rcov(iat)=2.60d0
     case('Mn') 
        rcov(iat)=2.50d0
     case('Fe') 
        rcov(iat)=2.50d0
     case('Co') 
        rcov(iat)=2.40d0
     case('Ni') 
        rcov(iat)=2.30d0
     case('Cu') 
        rcov(iat)=2.30d0
     case('Zn') 
        rcov(iat)=2.30d0
     case('Ga') 
        rcov(iat)=2.10d0
     case('Ge') 
        rcov(iat)=2.40d0
     case('As') 
        rcov(iat)=2.30d0
     case('Se') 
        rcov(iat)=2.30d0
     case('Br') 
        rcov(iat)=2.20d0
     case('Kr') 
        rcov(iat)=2.20d0
     case('Rb') 
        rcov(iat)=4.50d0
     case('Sr') 
        rcov(iat)=3.30d0
     case('Y' ) 
        rcov(iat)=3.30d0
     case('Zr') 
        rcov(iat)=3.00d0
     case('Nb') 
        rcov(iat)=2.92d0
     case('Mo') 
        rcov(iat)=2.83d0
     case('Tc') 
        rcov(iat)=2.75d0
     case('Ru') 
        rcov(iat)=2.67d0
     case('Rh') 
        rcov(iat)=2.58d0
     case('Pd') 
        rcov(iat)=2.50d0
     case('Ag') 
        rcov(iat)=2.50d0
     case('Cd') 
        rcov(iat)=2.50d0
     case('In') 
        rcov(iat)=2.30d0
     case('Sn') 
        rcov(iat)=2.66d0
     case('Sb') 
        rcov(iat)=2.66d0
     case('Te') 
        rcov(iat)=2.53d0
     case('I' ) 
        rcov(iat)=2.50d0
     case('Xe') 
        rcov(iat)=2.50d0
     case('Cs') 
        rcov(iat)=4.50d0
     case('Ba') 
        rcov(iat)=4.00d0
     case('La') 
        rcov(iat)=3.50d0
     case('Ce') 
        rcov(iat)=3.50d0
     case('Pr') 
        rcov(iat)=3.44d0
     case('Nd') 
        rcov(iat)=3.38d0
     case('Pm') 
        rcov(iat)=3.33d0
     case('Sm') 
        rcov(iat)=3.27d0
     case('Eu') 
        rcov(iat)=3.21d0
     case('Gd') 
        rcov(iat)=3.15d0
     case('Td') 
        rcov(iat)=3.09d0
     case('Dy') 
        rcov(iat)=3.03d0
     case('Ho') 
        rcov(iat)=2.97d0
     case('Er') 
        rcov(iat)=2.92d0
     case('Tm') 
        rcov(iat)=2.92d0
     case('Yb') 
        rcov(iat)=2.80d0
     case('Lu') 
        rcov(iat)=2.80d0
     case('Hf') 
        rcov(iat)=2.90d0
     case('Ta') 
        rcov(iat)=2.70d0
     case('W' ) 
        rcov(iat)=2.60d0
     case('Re') 
        rcov(iat)=2.60d0
     case('Os') 
        rcov(iat)=2.50d0
     case('Ir') 
        rcov(iat)=2.50d0
     case('Pt') 
        rcov(iat)=2.60d0
     case('Au') 
        rcov(iat)=2.70d0
     case('Hg') 
        rcov(iat)=2.80d0
     case('Tl') 
        rcov(iat)=2.50d0
     case('Pb') 
        rcov(iat)=3.30d0
     case('Bi') 
        rcov(iat)=2.90d0
     case('Po') 
        rcov(iat)=2.80d0
     case('At') 
        rcov(iat)=2.60d0
     case('Rn') 
        rcov(iat)=2.60d0
     case default
        call f_err_throw('(MH) no covalent radius stored for this atomtype '&
             //trim(astruct%atomnames(astruct%iatype(iat))),&
             err_name='BIGDFT_RUNTIME_ERROR')
     end select
     if (iproc == 0) then
        call yaml_map('(MHGPS) RCOV:'//trim(astruct%atomnames(astruct%iatype(iat))),rcov(iat))
     endif
  enddo
end subroutine give_rcov
end module

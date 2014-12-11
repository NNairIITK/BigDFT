!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module module_userinput
    use module_base !bigdft base module
    use module_global_variables
    implicit none

    private

    public :: userinput

    public :: read_input
    public :: write_input
    public :: print_input
    public :: give_rcov


type userinput
    !input parameters for mhgps
    integer          :: mhgps_verbosity     = 3
    logical          :: external_mini       = .false.

    !input parameters for sqnm saddle_search
    character(len=20) :: operation_mode             = 'connect'
    logical          :: random_minmode_guess       = .true.
    logical          :: saddle_biomode              = .false.
    integer          :: imode                      = 1
    integer          :: saddle_nit_trans           = 5000
    integer          :: saddle_nit_rot             = 5000
    integer          :: saddle_nhistx_trans        = 40
    integer          :: saddle_nhistx_rot          = 4
    real(gp)         :: saddle_alpha0_trans        = 1.e-3_gp
    real(gp)         :: saddle_alpha0_rot          = 1.e-3_gp
    real(gp)         :: saddle_alpha_stretch0      = 4.e-4_gp
    real(gp)         :: saddle_alpha_rot_stretch0  = 2.e-4_gp
    real(gp)         :: saddle_curvforcediff        = 1.e-3_gp
    real(gp)         :: saddle_rmsdispl0           = 0.04_gp
    real(gp)         :: saddle_trustr              = 0.2_gp
    real(gp)         :: saddle_tolc                = 7.0_gp
    real(gp)         :: saddle_tolf                = 7.0_gp
    logical         :: saddle_tighten             = .true.
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
    integer :: ts_guess_nstepsmax=5
    real(gp) :: lst_interpol_stepfrct=0.1_gp
    real(gp) :: lst_dt_max=5.0_gp
    real(gp) :: lst_fmax_tol=5.e-3_gp

    !variables for connect routine
    integer :: nsadmax=30

    !others
    real(gp)              :: en_delta_min, fp_delta_min
    real(gp)              :: en_delta_sad, fp_delta_sad



end type

contains

    subroutine read_input(uinp)
        use yaml_output
        use yaml_parse
        use dictionaries
        implicit none
        !parameter
        type(userinput), intent(inout) :: uinp
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
            call write_input(uinp)
            call yaml_warning('mhgps.inp does not exist.  Wrote '//&
                              'default input parameters to '//&
                              'mhgps.inp_default.')
             stop
        endif
        open(u,file=filename)
            read(u,*)uinp%mhgps_verbosity
            read(u,*)uinp%operation_mode, uinp%random_minmode_guess
            read(u,*)uinp%nsadmax
            read(u,*)uinp%external_mini
            read(u,*)uinp%en_delta_min,uinp%fp_delta_min
            read(u,*)uinp%en_delta_sad,uinp%fp_delta_sad
            read(u,*)uinp%saddle_biomode
            if(uinp%saddle_biomode)uinp%imode=2
            read(u,*) uinp%lst_interpol_stepfrct
            read(u,*) uinp%ts_guess_gammainv
            read(u,*) uinp%ts_guess_perpnrmtol
            read(u,*) uinp%ts_guess_trust
            read(u,*) uinp%ts_guess_nstepsmax
            read(u,*) uinp%lst_dt_max, uinp%lst_fmax_tol
            read(u,*)uinp%saddle_nit_trans, uinp%saddle_nit_rot
            read(u,*)uinp%saddle_nhistx_trans, uinp%saddle_nhistx_rot
            read(u,*)uinp%saddle_steepthresh_trans,uinp%saddle_steepthresh_rot
            read(u,*)uinp%saddle_fnrmtol
            if(uinp%saddle_biomode)then
                read(u,*)uinp%saddle_alpha0_trans, uinp%saddle_alpha0_rot,&
                     uinp%saddle_alpha_stretch0, uinp%saddle_alpha_rot_stretch0
            else
                read(u,*)uinp%saddle_alpha0_trans, uinp%saddle_alpha0_rot
            endif
            read(u,*)uinp%saddle_curvforcediff
            read(u,*)uinp%saddle_rmsdispl0,uinp%saddle_trustr
            read(u,*)uinp%saddle_tolc,uinp%saddle_tolf,uinp%saddle_tighten
            read(u,*)uinp%saddle_minoverlap0
            read(u,*)uinp%saddle_maxcurvrise
            read(u,*)uinp%saddle_cutoffratio
            read(u,*)uinp%saddle_recompIfCurvPos
            read(u,*)uinp%saddle_stepoff, uinp%saddle_scale_stepoff
            if(.not.uinp%external_mini)then
                read(u,*)uinp%mini_nhistx
                read(u,*)uinp%mini_ncluster_x
                read(u,*)uinp%mini_frac_fluct
                read(u,*)uinp%mini_forcemax
                read(u,*)uinp%mini_maxrise
                read(u,*)uinp%mini_betax
                read(u,*)uinp%mini_beta_stretchx
                read(u,*)uinp%mini_cutoffRatio
                read(u,*)uinp%mini_steepthresh
                read(u,*)uinp%mini_trustr
            endif

        close(u)
    end subroutine
    subroutine write_input(uinp)
        implicit none
        !parameter
        type(userinput), intent(in) :: uinp
        !local
        integer, parameter :: u=237
        character(17), parameter :: filename='mhgps.inp_default'
        open(u,file=filename)
            write(u,'(1x,i0.0,1x,1a)')uinp%mhgps_verbosity,' #mhgps_verbosity'
            write(u,'(1x,1L1,1x,1L1,1x,1a)')trim(adjustl(uinp%operation_mode)),uinp%random_minmode_guess,&
                 ' #mode, random_minmode_guess'
            write(u,'(1x,i0.0,1x,1a)')uinp%nsadmax,' #nsadmax'
            write(u,'(1x,1L1,1x,1a)')uinp%external_mini,&
                           ' #external minimizer'
            write(u,'(es10.3,1x,es10.3,1x,1a)')uinp%en_delta_min,uinp%fp_delta_min,' #en_delta_min, fp_delta_min'
            write(u,'(es10.3,1x,es10.3,1x,1a)')uinp%en_delta_sad,uinp%fp_delta_sad,' #en_delta_sad, fp_delta_sad'
            write(u,'(1x,1L1,1x,1a)')uinp%saddle_biomode,' #biomode'
            write(u,'(1x,es10.3,1x,1a)') uinp%lst_interpol_stepfrct,&
                 ' #inward interpolation distance as fraction of initial distance'
            write(u,'(1x,es10.3,1x,1a)') uinp%ts_guess_gammainv,' #step size for perpedicular optimization in freezing string method'
            write(u,'(1x,es10.3,1x,1a)') uinp%ts_guess_perpnrmtol,&
                 ' #convergence criterion perpedicular force in freezing string method'//&
                 ' (disable perpend. optim. by setting this value to a negative number)'
            write(u,'(1x,es10.3,1x,1a)') uinp%ts_guess_trust,' #trust radius freezing string method (maximum change of any coordinate'
            write(u,'(1x,i0,1x,1a)') uinp%ts_guess_nstepsmax,&
                 ' #maximum number of steps in perpendicular optimization in freezing stringmethod'
            write(u,'(1x,es10.3,1x,es10.3,1x,1a)')uinp%lst_dt_max, uinp%lst_fmax_tol,&
                 '#max. time step in fire optimizer of lst function, convergence criterion'
            write(u,'(1x,i0,1x,i0,1x,1a)')uinp%saddle_nit_trans, uinp%saddle_nit_rot,'  #nit_trans, nit_rot'
            write(u,'(1x,i0,1x,i0,1x,1a)')uinp%saddle_nhistx_trans, uinp%saddle_nhistx_rot,' #nhistx_trans, nhistx_rot'
            write(u,'(es10.3,1x,es10.3,1a)')uinp%saddle_steepthresh_trans,uinp%saddle_steepthresh_rot,&
                 ' #saddle_steepthresh_trans,saddle_steepthresh_rot'
            write(u,'(es10.3,1x,1a)')uinp%saddle_fnrmtol,' #fnrm tolerence convergence criterion for saddle point'
            if(uinp%saddle_biomode)then
                write(u,'(es10.3,3(1x,es10.3),1a)')uinp%saddle_alpha0_trans, uinp%saddle_alpha0_rot,&
                     uinp%saddle_alpha_stretch0, uinp%saddle_alpha_rot_stretch0,&
                     ' #alpha0_trans, alpha0_rot, alpha_stretch0, alpha_rot_stretch0'
            else
                write(u,'(es10.3,1x,es10.3,1a)')uinp%saddle_alpha0_trans, uinp%saddle_alpha0_rot,' #alpha0_trans, alpha0_rot'
            endif
            write(u,'(es10.3,1x,1a)')uinp%saddle_curvforcediff,' #curvforcedif'
            write(u,'(es10.3,1x,es10.3,1x,1a)')uinp%saddle_rmsdispl0,uinp%saddle_trustr,' #rmsdispl0, trustr'
            write(u,'(es10.3,1x,es10.3,1x,1L1,1x,1a)')uinp%saddle_tolc,uinp%saddle_tolf,uinp%saddle_tighten,' #tolc, tolf, tighten'
            write(u,'(es10.3,1x,1a)')uinp%saddle_minoverlap0,' #minoverlap0'
            write(u,'(es10.3,1x,1a)')uinp%saddle_maxcurvrise,' #maxcurvrise'
            write(u,'(es10.3,1x,1a)')uinp%saddle_cutoffratio,' #cutoffratio'
            write(u,'(1x,i0,1x,1a)')uinp%saddle_recompIfCurvPos,' #recompIfCurvPos'
            write(u,'(es10.3,1x,es10.3,1x,1a)')uinp%saddle_stepoff,uinp%saddle_scale_stepoff, ' #stepoff, stepoff_scale'
            write(u,'(1x,i0,1x,1a)')uinp%mini_nhistx,'#mini_nhistx'
            write(u,'(1x,i0,1x,1a)')uinp%mini_ncluster_x,'#mini_ncluster_x'
            write(u,'(es10.3,1x,1a)')uinp%mini_frac_fluct,'#mini_frac_fluct'
            write(u,'(es10.3,1x,1a)')uinp%mini_forcemax,'#mini_forcemax'
            write(u,'(es10.3,1x,1a)')uinp%mini_maxrise,'#mini_maxrise'
            write(u,'(es10.3,1x,1a)')uinp%mini_betax,'#mini_betax'
            write(u,'(es10.3,1x,1a)')uinp%mini_beta_stretchx,'#mini_beta_stretchx'
            write(u,'(es10.3,1x,1a)')uinp%mini_cutoffRatio,'#mini_cutoffRatio'
            write(u,'(es10.3,1x,1a)')uinp%mini_steepthresh,'#mini_steepthresh'
            write(u,'(es10.3,1x,1a)')uinp%mini_trustr,'#mini_trustr'
        close(u)
    end subroutine
    subroutine print_input(uinp)
        use module_global_variables, only: mhgps_version
        use yaml_output
        implicit none
        !parameters
        type(userinput), intent(in) :: uinp
        !local
        call yaml_comment('(MHGPS) Input Parameters',hfill='-')
        call yaml_map('(MHGPS) mhgps_verbosity',uinp%mhgps_verbosity)
        call yaml_map('(MHGPS) operation_mode',trim(adjustl(uinp%operation_mode)))
        call yaml_map('(MHGPS) random_minmode_guess',uinp%random_minmode_guess)
        call yaml_map('(MHGPS) nsadmax',uinp%nsadmax)
        call yaml_map('(MHGPS) external minimizer',uinp%external_mini)
        call yaml_map('(MHGPS) en_delta_min',uinp%en_delta_min)
        call yaml_map('(MHGPS) en_delta_sad',uinp%en_delta_sad)
        call yaml_map('(MHGPS) Biomolecule mode',uinp%saddle_biomode)
        call yaml_map('(MHGPS) lst_interpol_stepfrct',uinp%lst_interpol_stepfrct)
        call yaml_map('(MHGPS) ts_guess_gammainv',uinp%ts_guess_gammainv)
        call yaml_map('(MHGPS) ts_guess_perpnrmtol',uinp%ts_guess_perpnrmtol)
        call yaml_map('(MHGPS) ts_guess_trust',uinp%ts_guess_trust)
        call yaml_map('(MHGPS) ts_guess_nstepsmax',uinp%ts_guess_nstepsmax)
        call yaml_map('(MHGPS) lst_dt_max',uinp%lst_dt_max)
        call yaml_map('(MHGPS) lst_fmax_tol',uinp%lst_fmax_tol)
        call yaml_map('(MHGPS) saddle_nit_trans',uinp%saddle_nit_trans)
        call yaml_map('(MHGPS) saddle_nit_rot',uinp%saddle_nit_rot)
        call yaml_map('(MHGPS) saddle_nhistx_trans',uinp%saddle_nhistx_trans)
        call yaml_map('(MHGPS) saddle_nhistx_rot',uinp%saddle_nhistx_rot)
        call yaml_map('(MHGPS) saddle_steepthresh_trans',uinp%saddle_steepthresh_trans)
        call yaml_map('(MHGPS) saddle_steepthresh_rot',uinp%saddle_steepthresh_rot)
        call yaml_map('(MHGPS) saddle_fnrmtol',uinp%saddle_fnrmtol)
        call yaml_map('(MHGPS) saddle_alpha0_trans',uinp%saddle_alpha0_trans)
        call yaml_map('(MHGPS) saddle_alpha0_rot',uinp%saddle_alpha0_rot)
        if(uinp%saddle_biomode)then
            call yaml_map('(MHGPS) saddle_alpha_stretch0',uinp%saddle_alpha_stretch0)
            call yaml_map('(MHGPS) saddle_alpha_rot_stretch0',uinp%saddle_alpha_rot_stretch0)
        endif
        call yaml_map('(MHGPS) saddle_curvgraddiff',uinp%saddle_curvforcediff)
        call yaml_map('(MHGPS) saddle_rmsdispl0',uinp%saddle_rmsdispl0)
        call yaml_map('(MHGPS) saddle_trustr',uinp%saddle_trustr)
        call yaml_map('(MHGPS) saddle_tolc',uinp%saddle_tolc)
        call yaml_map('(MHGPS) saddle_tolf',uinp%saddle_tolf)
        call yaml_map('(MHGPS) saddle_tighten',uinp%saddle_tighten)
        call yaml_map('(MHGPS) saddle_maxcurvrise',uinp%saddle_maxcurvrise)
        call yaml_map('(MHGPS) saddle_cutoffratio',uinp%saddle_cutoffratio)
        call yaml_map('(MHGPS) saddle_recompIfCurvPos',uinp%saddle_recompIfCurvPos)
        call yaml_map('(MHGPS) saddle_stepoff',uinp%saddle_stepoff)
        call yaml_map('(MHGPS) saddle_scale_stepoff',uinp%saddle_scale_stepoff)
        if(.not.uinp%external_mini)then
            call yaml_map('(MHGPS) mini_nhistx', uinp%mini_nhistx)
            call yaml_map('(MHGPS) mini_ncluster_x', uinp%mini_ncluster_x)
            call yaml_map('(MHGPS) mini_frac_fluct', uinp%mini_frac_fluct)
            call yaml_map('(MHGPS) mini_forcemax', uinp%mini_forcemax)
            call yaml_map('(MHGPS) mini_maxrise', uinp%mini_maxrise)
            call yaml_map('(MHGPS) mini_betax', uinp%mini_betax)
            call yaml_map('(MHGPS) mini_beta_stretchx', uinp%mini_beta_stretchx)
            call yaml_map('(MHGPS) mini_cutoffRatio', uinp%mini_cutoffRatio)
            call yaml_map('(MHGPS) mini_steepthresh', uinp%mini_steepthresh)
            call yaml_map('(MHGPS) mini_trustr', uinp%mini_trustr)
        endif
    end subroutine print_input


subroutine give_rcov(astruct,nat,rcov)
  use module_base, only: gp
  use module_types
  use yaml_output
  use module_global_variables, only: iproc
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

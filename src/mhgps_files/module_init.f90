module module_init
    use module_base !bigdft base module
    use module_global_variables
    implicit none

contains

    subroutine read_input()
        use yaml_output
        implicit none
        !parameter
        !internal
        integer, parameter :: u=237
        character(9), parameter :: filename='mhgps.inp'
        logical :: exists

        inquire(file=filename,exist=exists)
        if(.not. exists)then
            call write_input()
            call yaml_warning('mhgps.inp does not exist. &
                Wrote default input parameters to mhgps.inp_default.')
            stop
        endif
        open(u,file=filename)
            read(u,*)mhgps_verbosity
            read(u,*)saddle_connect, random_minmode_guess
            read(u,*)efmethod,external_mini
            read(u,*)saddle_biomode
            if(saddle_biomode)imode=2
            read(u,*)saddle_nit_trans, saddle_nit_rot
            read(u,*)saddle_nhistx_trans, saddle_nhistx_rot
            read(u,*)saddle_fnrmtol
            if(saddle_biomode)then
                read(u,*)saddle_alpha0_trans, saddle_alpha0_rot, saddle_alpha_stretch0, saddle_alpha_rot_stretch0
            else
                read(u,*)saddle_alpha0_trans, saddle_alpha0_rot
            endif
            read(u,*)saddle_curvgraddiff
            read(u,*)saddle_rmsdispl0,saddle_trustr
            read(u,*)saddle_tolc,saddle_tolf
            read(u,*)saddle_tightenfac
            read(u,*)saddle_maxcurvrise
            read(u,*)saddle_cutoffratio
            read(u,*)saddle_recompIfCurvPos
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
            write(u,'(xi0.0,xa)')mhgps_verbosity,' #mhgps_verbosity'
            write(u,'(xL,xL,xa)')saddle_connect,random_minmode_guess, ' #connect, random_minmode_guess'
            write(u,'(xa,xa)')trim(adjustl(efmethod)),' #efmethod'
            write(u,'(xL,xa)')saddle_biomode,' #biomode'
            write(u,'(xi0,xi0,xa)')saddle_nit_trans, saddle_nit_rot,'  #nit_trans, not_rot'
            write(u,'(xi0,xi0,xa)')saddle_nhistx_trans, saddle_nhistx_rot,' #nhistx_trans, nhistx_rot'
            write(u,'(es10.3,xa)')saddle_fnrmtol,' #fnrm tolerence convergence criterion for saddle point'
            if(saddle_biomode)then
                write(u,'(es10.3,3(xes10.3),a)')saddle_alpha0_trans, saddle_alpha0_rot, saddle_alpha_stretch0, saddle_alpha_rot_stretch0,' #alpha0_trans, alpha0_rot, alpha_stretch0, alpha_rot_stretch0'
            else
                write(u,'(es10.3,xes10.3,a)')saddle_alpha0_trans, saddle_alpha0_rot,' #alpha0_trans, alpha0_rot'
            endif
            write(u,'(es10.3,xa)')saddle_curvgraddiff,' #curvgraddif'
            write(u,'(es10.3,xes10.3,xa)')saddle_rmsdispl0,saddle_trustr,' #rmsdispl0, trustr'
            write(u,'(es10.3,xes10.3,xa)')saddle_tolc,saddle_tolf,' #tolc, tolf'
            write(u,'(es10.3,xa)')saddle_tightenfac,' #tightenfac'
            write(u,'(es10.3,xa)')saddle_maxcurvrise,' #maxcurvrise'
            write(u,'(es10.3,xa)')saddle_cutoffratio,' #cutoffratio'
            write(u,'(xi0,xa)')saddle_recompIfCurvPos,' #recompIfCurvPos'
        close(u)
    end subroutine
    subroutine print_input()
        use module_global_variables, only: mhgps_version
        use yaml_output
        call yaml_comment('(MHGPS) Input Parameters',hfill='-')
        call yaml_map('(MHGPS) mhgps_verbosity',mhgps_verbosity)
        call yaml_map('(MHGPS) saddle_connect',saddle_connect)
        call yaml_map('(MHGPS) random_minmode_guess',random_minmode_guess)
        call yaml_map('(MHGPS) Energy and forces method',trim(adjustl(efmethod)))
        call yaml_map('(MHGPS) Biomolecule mode',saddle_biomode)
        call yaml_map('(MHGPS) saddle_nit_trans',saddle_nit_trans)
        call yaml_map('(MHGPS) saddle_nit_rot',saddle_nit_rot)
        call yaml_map('(MHGPS) saddle_nhistx_trans',saddle_nhistx_trans)
        call yaml_map('(MHGPS) saddle_nhistx_rot',saddle_nhistx_rot)
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
        call yaml_map('(MHGPS) saddle_tightenfac',saddle_tightenfac)
        call yaml_map('(MHGPS) saddle_maxcurvrise',saddle_maxcurvrise)
        call yaml_map('(MHGPS) saddle_cutoffratio',saddle_cutoffratio)
        call yaml_map('(MHGPS) saddle_recompIfCurvPos',saddle_recompIfCurvPos)
    end subroutine print_input

    subroutine print_logo_mhgps()
        use module_global_variables, only: mhgps_version
        use yaml_output
        implicit none

        call yaml_comment('(MHGPS) Minima Hopping Guided Path Sampling',hfill='=')
        
        call yaml_mapping_open('(MHGPS) logo')
        call yaml_scalar('(MHGPS)      ___           ___           ___           ___           ___     ') 
        call yaml_scalar('(MHGPS)     /\__\         /\__\         /\  \         /\  \         /\  \    ')
        call yaml_scalar('(MHGPS)    /::|  |       /:/  /        /::\  \       /::\  \       /::\  \   ')
        call yaml_scalar('(MHGPS)   /:|:|  |      /:/__/        /:/\:\  \     /:/\:\  \     /:/\ \  \  ')
        call yaml_scalar('(MHGPS)  /:/|:|__|__   /::\  \ ___   /:/  \:\  \   /::\~\:\  \   _\:\~\ \  \ ')
        call yaml_scalar('(MHGPS) /:/ |::::\__\ /:/\:\  /\__\ /:/__/_\:\__\ /:/\:\ \:\__\ /\ \:\ \ \__\')
        call yaml_scalar('(MHGPS) \/__/~~/:/  / \/__\:\/:/  / \:\  /\ \/__/ \/__\:\/:/  / \:\ \:\ \/__/')
        call yaml_scalar('(MHGPS)       /:/  /       \::/  /   \:\ \:\__\        \::/  /   \:\ \:\__\  ')
        call yaml_scalar('(MHGPS)      /:/  /        /:/  /     \:\/:/  /         \/__/     \:\/:/  /  ')
        call yaml_scalar('(MHGPS)     /:/  /        /:/  /       \::/  /                     \::/  /   ')
        call yaml_scalar('(MHGPS)     \/__/         \/__/         \/__/                       \/__/    ')
        call yaml_scalar('(MHGPS)                                                   as post-processing ')
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


subroutine give_rcov(atoms,nat,rcov)
  use module_base, only: gp
  use module_types
  use yaml_output
  use module_global_variables, only: iproc
  implicit none
  !Arguments
  integer, intent(in) :: nat
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: rcov(nat)
  !Local variables
  integer :: iat

  do iat=1,nat
     if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='H') then
        rcov(iat)=0.75d0
!        rcov(iat)=0.75d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='LJ') then
        rcov(iat)=0.56d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='He') then
        rcov(iat)=0.75d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Li') then
        rcov(iat)=3.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Be') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='B' ) then
        rcov(iat)=1.55d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='C' ) then
        rcov(iat)=1.45d0
!        rcov(iat)=1.45d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='N' ) then
        rcov(iat)=1.42d0
!        rcov(iat)=1.42d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='O' ) then
        rcov(iat)=1.38d0
!        rcov(iat)=1.38d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='F' ) then
        rcov(iat)=1.35d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ne') then
        rcov(iat)=1.35d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Na') then
        rcov(iat)=3.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mg') then
        rcov(iat)=2.65d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Al') then
        rcov(iat)=2.23d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Si') then
        rcov(iat)=2.09d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='P' ) then
        rcov(iat)=2.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='S' ) then
        rcov(iat)=1.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cl') then
        rcov(iat)=1.87d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ar') then
        rcov(iat)=1.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='K' ) then
        rcov(iat)=4.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ca') then
        rcov(iat)=3.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sc') then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ti') then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='V' ) then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cr') then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mn') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Fe') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Co') then
        rcov(iat)=2.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ni') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cu') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Zn') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ga') then
        rcov(iat)=2.10d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ge') then
        rcov(iat)=2.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='As') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Se') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Br') then
        rcov(iat)=2.20d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Kr') then
        rcov(iat)=2.20d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rb') then
        rcov(iat)=4.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sr') then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Y' ) then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Zr') then
        rcov(iat)=3.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Nb') then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mo') then
        rcov(iat)=2.83d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tc') then
        rcov(iat)=2.75d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ru') then
        rcov(iat)=2.67d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rh') then
        rcov(iat)=2.58d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pd') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ag') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cd') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='In') then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sn') then
        rcov(iat)=2.66d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sb') then
        rcov(iat)=2.66d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Te') then
        rcov(iat)=2.53d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='I' ) then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Xe') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cs') then
        rcov(iat)=4.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ba') then
        rcov(iat)=4.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='La') then
        rcov(iat)=3.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ce') then
        rcov(iat)=3.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pr') then
        rcov(iat)=3.44d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Nd') then
        rcov(iat)=3.38d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pm') then
        rcov(iat)=3.33d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sm') then
        rcov(iat)=3.27d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Eu') then
        rcov(iat)=3.21d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Gd') then
        rcov(iat)=3.15d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Td') then
        rcov(iat)=3.09d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Dy') then
        rcov(iat)=3.03d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ho') then
        rcov(iat)=2.97d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Er') then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tm') then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Yb') then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Lu') then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Hf') then
        rcov(iat)=2.90d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ta') then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='W' ) then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Re') then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Os') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ir') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pt') then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Au') then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Hg') then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tl') then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pb') then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Bi') then
        rcov(iat)=2.90d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Po') then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='At') then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rn') then
        rcov(iat)=2.60d0
     else
        call yaml_comment('(MH) no covalent radius stored for this atomtype '&
             //trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))))
        stop
     endif
     if (iproc == 0) then
        call yaml_map('(MHGPS) RCOV:'//trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),rcov(iat))
     endif
  enddo
end subroutine give_rcov
end module

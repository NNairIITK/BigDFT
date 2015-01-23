!> @file
!!    SQNS Saddle Search Method
!! @author Bastian Schaefer
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module saddle for SQNS
module module_saddle
use module_base
implicit none

private

public :: findsad
public :: findsad_work
public :: allocate_finsad_workarrays
public :: deallocate_finsad_workarrays

type findsad_work
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
    real(gp), allocatable :: rxyzold_trans(:,:)
    real(gp), allocatable :: dds_trans(:,:)
    real(gp), allocatable :: dd0_trans(:,:)
    real(gp), allocatable :: delta_trans(:,:)
    real(gp), allocatable :: ftmp_trans(:,:)
    real(gp), allocatable :: minmodeold_trans(:,:)
    real(gp), allocatable :: minmode0_trans(:,:)
    integer               :: lwork_trans
    real(gp), allocatable :: work_trans(:)




    !variables for rotation
    integer              :: nhist_rot,ndim_rot
    real(gp)             :: alpha_rot, alpha_stretch_rot
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
    integer               :: lwork_rot
    real(gp), allocatable :: work_rot(:)
end type

contains
!=====================================================================
subroutine allocate_finsad_workarrays(runObj,uinp,nbond,fsw)
    use module_base
    use dynamic_memory
    use bigdft_run
    use module_userinput
    implicit none
    !parameters
    type(run_objects), intent(in)          :: runObj
    type(findsad_work), intent(out) :: fsw
    type(userinput), intent(in)            :: uinp
    integer, intent(in)            :: nbond
    !internal
    integer :: nat, info
    real(gp) :: wd(1)
    nat = bigdft_nat(runObj)

    fsw%rxyz_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='rxyz_rot')
    fsw%fxyz_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='fxyz_rot')
    fsw%fxyzraw_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='fxyzraw_rot')
    fsw%rxyzraw_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='rxyzraw_rot')
    fsw%fstretch_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='fstretch_rot')
    fsw%eval_rot = f_malloc((/1.to.uinp%saddle_nhistx_rot/),id='eval_rot')
    fsw%res_rot = f_malloc((/1.to.uinp%saddle_nhistx_rot/),id='res_rot')
    fsw%rrr_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='rrr_rot')
    fsw%aa_rot = f_malloc((/1.to.uinp%saddle_nhistx_rot,&
             1.to.uinp%saddle_nhistx_rot/),id='aa_rot')
    fsw%ff_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='ff_rot')
    fsw%rr_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='rr_rot')
    fsw%dd_rot = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='dd_rot')
    fsw%fff_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_rot/),id='fff_rot')
    fsw%scpr_rot = f_malloc((/ 1.to.uinp%saddle_nhistx_rot/),id='scpr_rot')
    fsw%rxyz_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='rxyz_trans')
    fsw%fxyz_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='fxyz_trans')
    fsw%fxyzraw_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='fxyzraw_trans')
    fsw%rxyzraw_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='rxyzraw_trans')
    fsw%fstretch_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='fstretch_trans')
    fsw%eval_trans = f_malloc((/1.to.uinp%saddle_nhistx_trans/),&
                 id='eval_trans')
    fsw%res_trans = f_malloc((/1.to.uinp%saddle_nhistx_trans/),id='res_trans')
    fsw%rrr_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='rrr_trans')
    fsw%aa_trans = f_malloc((/1.to.uinp%saddle_nhistx_trans,&
                1.to.uinp%saddle_nhistx_trans/),id='aa_trans')
    fsw%ff_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='ff_trans')
    fsw%rr_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.uinp%saddle_nhistx_trans/),id='rr_trans')
    fsw%dd_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='dd_trans')
    fsw%fff_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                -1.to.uinp%saddle_nhistx_trans/),id='fff_trans')
    fsw%scpr_trans = f_malloc((/ 1.to.uinp%saddle_nhistx_trans/),&
                 id='scpr_trans')
    fsw%rxyzold_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzold_trans')
    fsw%dds_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='dds_trans')
    fsw%dd0_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='dd0_trans')
    fsw%delta_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='delta_trans')
    fsw%ftmp_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='ftmp_trans')
    fsw%minmodeold_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='minmodeold_trans')
    fsw%minmode0_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='minmode0_trans')
    fsw%wold_trans = f_malloc((/ 1.to.nbond/),id='wold_trans')
    fsw%wold_rot = f_malloc((/ 1.to.nbond/),id='wold_rot')

    call DSYEV('N','L',uinp%saddle_nhistx_trans,fsw%aa_trans,&
         uinp%saddle_nhistx_trans,fsw%eval_trans,wd,-1,info)
    if (info.ne.0) stop 'info query'
    fsw%lwork_trans=nint(wd(1))
    fsw%work_trans = f_malloc((/ 1.to.fsw%lwork_trans/),id='work_trans')

    call DSYEV('N','L',uinp%saddle_nhistx_rot,fsw%aa_rot,&
         uinp%saddle_nhistx_rot,fsw%eval_rot,wd,-1,info)
    if (info.ne.0) stop 'info query'
    fsw%lwork_rot=nint(wd(1))
    fsw%work_rot = f_malloc((/ 1.to.fsw%lwork_rot/),id='work_rot')

end subroutine
!=====================================================================
subroutine deallocate_finsad_workarrays(fsw)
    use dynamic_memory
    implicit none
    !parameters
    type(findsad_work), intent(inout) :: fsw
    call f_free(fsw%rxyz_rot)
    call f_free(fsw%fxyz_rot)
    call f_free(fsw%fxyzraw_rot)
    call f_free(fsw%rxyzraw_rot)
    call f_free(fsw%fstretch_rot)
    call f_free(fsw%eval_rot)
    call f_free(fsw%res_rot)
    call f_free(fsw%rrr_rot)
    call f_free(fsw%aa_rot)
    call f_free(fsw%ff_rot)
    call f_free(fsw%rr_rot)
    call f_free(fsw%dd_rot)
    call f_free(fsw%fff_rot)
    call f_free(fsw%scpr_rot)
    call f_free(fsw%rxyz_trans)
    call f_free(fsw%fxyz_trans)
    call f_free(fsw%fxyzraw_trans)
    call f_free(fsw%rxyzraw_trans)
    call f_free(fsw%fstretch_trans)
    call f_free(fsw%eval_trans)
    call f_free(fsw%res_trans)
    call f_free(fsw%rrr_trans)
    call f_free(fsw%aa_trans)
    call f_free(fsw%ff_trans)
    call f_free(fsw%rr_trans)
    call f_free(fsw%dd_trans)
    call f_free(fsw%fff_trans) 
    call f_free(fsw%scpr_trans)
    call f_free(fsw%rxyzold_trans)
    call f_free(fsw%dds_trans)
    call f_free(fsw%dd0_trans)
    call f_free(fsw%delta_trans)
    call f_free(fsw%ftmp_trans)
    call f_free(fsw%minmodeold_trans)
    call f_free(fsw%minmode0_trans)
    call f_free(fsw%wold_trans)
    call f_free(fsw%wold_rot)
    call f_free(fsw%work_trans)
    call f_free(fsw%work_rot)

end subroutine
!=====================================================================
subroutine findsad(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,iconnect,&
                  wpos,etot,fout,minmode,displ,ener_count,&
                  rotforce,converged)
    !imode=1 for clusters
    !imode=2 for biomolecules
    use module_base
    use module_atoms, only: astruct_dump_to_file
    use bigdft_run
    use yaml_output
    use module_interfaces
    use module_sqn
    use module_userinput
    use module_mhgps_state
    implicit none
    !parameters    
    type(mhgps_state), intent(inout) :: mhgpsst
    type(userinput), intent(in) :: uinp
    type(findsad_work), intent(inout)      :: fsw
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    real(gp), intent(in)      :: rcov(runObj%atoms%astruct%nat)
    real(gp), intent(inout)   :: wpos(3,runObj%atoms%astruct%nat)
    real(gp), intent(out)     :: etot
    real(gp), intent(out)     :: fout(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout)   :: minmode(3,runObj%atoms%astruct%nat)
    real(gp), intent(out)   :: displ
    real(gp), intent(inout)   :: ener_count
    logical, intent(out)      :: converged
    integer, intent(in)       :: nbond
    integer, intent(in)       :: iconnect(2,nbond)
    real(gp), intent(in)      :: rotforce(3,runObj%atoms%astruct%nat)
    !internal
!    real(gp), allocatable, dimension(:,:)   :: rxyzold
!    real(gp), allocatable, dimension(:,:)   :: dds
!    real(gp), allocatable, dimension(:,:)   :: dd0
!    real(gp), allocatable, dimension(:,:)   :: delta
!    real(gp), allocatable, dimension(:,:)   :: ftmp
!    real(gp), allocatable, dimension(:,:)   :: minmodeold
!    real(gp), allocatable, dimension(:,:)   :: minmode0
    logical  :: steep
    real(gp) :: maxd
    real(gp) :: alpha
    real(gp) :: alpha_stretch
    real(gp) :: fnrm
    real(gp) :: fmax
    real(gp) :: etotp
    real(gp) :: etotold
    real(gp) :: detot
    real(gp) :: cosangle
    real(gp) :: scl
    real(gp) :: tt
    real(gp) :: dt
    integer  :: iat
    integer  :: itswitch
    integer  :: nhist
    integer  :: ndim
    integer  :: it
    integer  :: ihist
    integer  :: i
    integer  :: recompute
    integer  :: fc
    integer  :: icheck
    integer  :: icheckmax
    real(gp) :: tol
    real(gp) :: displold
    real(gp) :: curv
    real(gp) :: overlap
    real(gp) :: minoverlap
    real(gp) :: tnatdmy(3,runObj%atoms%astruct%nat)
    real(gp) :: tmp
    logical  :: optCurvConv
    logical  :: tooFar
    logical  :: fixfragmented
    logical  :: subspaceSucc
    logical  :: tighten
    character(len=9)   :: fn9
    character(len=300)  :: comment
    real(gp) :: curvgraddiff_tmp
    !functions
    real(gp) :: ddot,dnrm2

!    if(bigdft_get_geocode(runObj)/='F'.and. .not.&
!            (trim(adjustl(char(runObj%run_mode)))==&
!                        'LENOSKY_SI_CLUSTERS_RUN_MODE'))then
!        stop 'STOP: saddle search only implemented for free BC'
!    endif

    if((uinp%saddle_minoverlap0>=-1.0_gp).and.uinp%saddle_tighten)&
    stop 'STOP: Do not use minoverlap and no tightening in combination'

    if(mhgpsst%iproc==0)then
        call yaml_comment('(MHGPS) Start Saddle Search ....',&
                          hfill='-')
    endif

    fsw%ndim_rot=0
    fsw%nhist_rot=0
    fsw%alpha_rot=uinp%saddle_alpha0_rot
    fsw%alpha_stretch_rot=uinp%saddle_alpha_stretch0
!    flag=.true.
    fc=0
    fixfragmented=.false.
    converged=.false.
    subspaceSucc=.true.
    optCurvConv=.true.
    tol=uinp%saddle_tolc
    displ=0.0_gp
    displold=0.0_gp
    detot=0.0_gp
    curv=1000.0_gp
    icheck=0
!    icheckmax=5
    icheckmax=0
    if(icheckmax==0 .and. uinp%saddle_tighten) icheckmax=1
    tighten=.false.
    alpha_stretch=uinp%saddle_alpha_stretch0

    ! allocate arrays
    fsw%minmode0_trans=minmode
    fsw%minmodeold_trans=minmode
    fsw%wold_trans=0.0_gp
    fsw%fstretch_trans=0.0_gp
    fsw%rxyz_trans(:,:,0)=wpos

    if (bigdft_get_geocode(runObj) == 'F') then
    call fixfrag_posvel(mhgpsst,bigdft_get_geocode(runObj),runObj%atoms%astruct%nat,rcov,&
         fsw%rxyz_trans(1,1,0),tnatdmy,1,fixfragmented)
        if(fixfragmented .and. uinp%mhgps_verbosity >=0.and. &
          mhgpsst%iproc==0) call yaml_comment('fragmentation fixed')
    endif

    runObj%inputs%inputPsiId=0
    call minenergyandforces(mhgpsst,.true.,uinp%imode,runObj,outs,&
         fsw%rxyz_trans(1,1,0),fsw%rxyzraw_trans(1,1,0),&
         fsw%fxyz_trans(1,1,0),fsw%fstretch_trans(1,1,0),&
         fsw%fxyzraw_trans(1,1,0),etot,iconnect,nbond,fsw%wold_trans,&
         uinp%saddle_alpha_stretch0,alpha_stretch)
    fsw%rxyzold_trans=fsw%rxyz_trans(:,:,0)
    ener_count=ener_count+1.0_gp
    if(uinp%imode==2)then
        fsw%rxyz_trans(:,:,0)=fsw%rxyz_trans(:,:,0)+&
                              alpha_stretch*fsw%fstretch_trans(:,:,0)
    endif
    call fnrmandforcemax(fsw%fxyzraw_trans(1,1,0),fnrm,fmax,&
         runObj%atoms%astruct%nat)

    fnrm=sqrt(fnrm)
    etotold=etot
    etotp=etot

!    itswitch=2
    itswitch=-2
    ndim=0
    nhist=0
    alpha=uinp%saddle_alpha0_trans
    if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)then
        write(*,'(a)')&
        '  #(MHGPS) METHOD COUNT  IT  Energy                '//&
        'DIFF      FMAX      FNRM      alpha    ndim dspl         '//&
        'alpha_strtch'
        write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2),1x,'//&
                 'i3.3,1x,es12.5,1x,es9.2)')&
        '   (MHGPS) GEOPT ',nint(ener_count),0,etotp,detot,fmax,&
        fnrm, alpha,ndim,displ,alpha_stretch
    endif
    do it=1,uinp%saddle_nit_trans
        nhist=nhist+1

        if ((.not. subspaceSucc) .or. &
             fnrm .gt. uinp%saddle_steepthresh_trans  .or. &
             it.le.itswitch) then
            ndim=0
            steep=.true.
            !if (it.gt.itswitch) itswitch=it+nhistx_trans
        else
            steep=.false.
            !alpha=alpha0_trans
        endif

        !make space in the history list
        if (nhist.gt.uinp%saddle_nhistx_trans) then
            nhist=uinp%saddle_nhistx_trans
            do ihist=0,nhist-1
                do iat=1,runObj%atoms%astruct%nat
                    do i=1,3
                        fsw%rxyz_trans(i,iat,ihist)=fsw%rxyz_trans(i,iat,ihist+1)
                        fsw%fxyz_trans(i,iat,ihist)=fsw%fxyz_trans(i,iat,ihist+1)
                        fsw%rxyzraw_trans(i,iat,ihist)=fsw%rxyzraw_trans(i,iat,ihist+1)
                        fsw%fxyzraw_trans(i,iat,ihist)=fsw%fxyzraw_trans(i,iat,ihist+1)
                        fsw%fstretch_trans(i,iat,ihist)=fsw%fstretch_trans(i,iat,ihist+1)
                     enddo
                enddo
            enddo
        endif

        !START FINDING LOWEST MODE

        !Walked too far? Then recompute direction of lowest mode!
        tooFar = abs(displ-displold)>uinp%saddle_rmsdispl0*sqrt(dble(3*runObj%atoms%astruct%nat))
 
        !determine if final tightening should be done:
        if(tighten.and.uinp%saddle_tighten)then
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)then
                call yaml_comment('(MHGPS) tightening')
            endif
            tol=uinp%saddle_tolf
            recompute=it
            minoverlap=-2._gp !disable overlap control in opt_curv
        else
            tol=uinp%saddle_tolc
            minoverlap=uinp%saddle_minoverlap0
        endif
        if(tooFar& !recompute lowest mode if walked too far
          .or. it==1& !compute lowest mode at first step
!          .or. (.not. optCurvConv)&
          .or. (curv>=0.0_gp .and. &
               ((mod(it,uinp%saddle_recompIfCurvPos)==0).or. fnrm<uinp%saddle_fnrmtol))&
                                                  !For LJ
                                                  !systems
                                                  !recomputation
                                                  !every
                                                  !nth=recompIfCurvPos
                                                  !step raises
                                                  !stability
          .or.recompute==it)then
            !if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)call yaml_comment(&
            !'(MHGPS) METHOD COUNT  IT  CURVATURE             &
            !DIFF      FMAX      FNRM      alpha    ndim')
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)write(*,'(a)')&
            '  #(MHGPS) METHOD COUNT  IT  CURVATURE             '//&
            'DIFF      GMAX      GNRM      alpha    ndim '//&
            'alpha_strtch overl. displr       displp'
            runObj%inputs%inputPsiId=1
            !determine finite difference
            call clean_minmode(runObj%atoms%astruct%nat,&
                 bigdft_get_geocode(runObj),&
                 fsw%rxyzraw_trans(1,1,nhist-1),minmode(1,1))
            
            dt=0.0_gp
            maxd=-huge(1.0_gp)
            do iat=1,runObj%atoms%astruct%nat
                dt=minmode(1,iat)**2+minmode(2,iat)**2+minmode(3,iat)**2
                maxd=max(maxd,dt)
            enddo
            maxd=sqrt(maxd)
            curvgraddiff_tmp = uinp%saddle_curvforcediff / maxd
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)then
                call yaml_map('  (MHGPS) Finite difference spacing '//&
                     'for curvature computation',curvgraddiff_tmp,&
                     fmt='(1pe21.14)')
            endif

             !inputPsiId=0
            call opt_curv(mhgpsst,it,uinp%imode,fsw,uinp,runObj,outs,uinp%saddle_alpha0_rot,&
                          curvgraddiff_tmp,uinp%saddle_nit_rot,uinp%saddle_nhistx_rot,&
!                          uinp%saddle_curvforcediff,uinp%saddle_nit_rot,uinp%saddle_nhistx_rot,&
                          fsw%rxyzraw_trans(1,1,nhist-1),fsw%fxyzraw_trans(1,1,nhist-1),&
                          minmode(1,1),curv,rotforce(1,1),tol,&
                          ener_count,optCurvConv,iconnect,nbond,&
                          uinp%saddle_alpha_rot_stretch0,uinp%saddle_maxcurvrise,uinp%saddle_cutoffratio,&
                          minoverlap)

            runObj%inputs%inputPsiId=1
            minmode = minmode / dnrm2(3*runObj%atoms%astruct%nat,minmode(1,1),1)
!            if(.not.optCurvConv)then
!                if(mhgpsst%iproc==0)call yaml_warning('(MHGPS) opt_curv '//&
!                                              'failed')
!                converged=.false.
! stop 'opt_curv failed'
!                return
!            endif
            overlap=ddot(3*runObj%atoms%astruct%nat,fsw%minmodeold_trans(1,1),1,minmode(1,1),1)
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
                call yaml_map('  (MHGPS) minmode overlap',overlap)
            if((.not.optCurvConv).and. (overlap <0.85d0))then
                minmode=fsw%minmodeold_trans
            endif
            fsw%minmodeold_trans=minmode
            displold=displ
            recompute=huge(1)
            !if(mhgpsst%iproc==0.and.mhgps_verbosity>=2)call yaml_comment(&
            !'(MHGPS) METHOD COUNT  IT  Energy                &
            !DIFF      FMAX      FNRM      alpha    ndim')
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)write(*,'(a)')&
            '  #(MHGPS) METHOD COUNT  IT  Energy                '//&
            'DIFF      FMAX      FNRM      alpha    ndim dspl   '//&
            '      alpha_strtch'
        endif
        !END FINDING LOWEST MODE
        
        call modify_gradient(runObj%atoms%astruct%nat,ndim,fsw%rrr_trans(1,1,1),fsw%eval_trans(1),&
             fsw%res_trans(1),fsw%fxyz_trans(1,1,nhist-1),alpha,fsw%dd_trans(1,1))
 
        !save a version of dd without minmode direction in fsw%dd0_trans
        !(used for gradient feedback)
        !fsw%dd0_trans=dd-ddot(3*nat,dd(1,1),1,minmode(1,1),1)*minmode
        tmp=-ddot(3*runObj%atoms%astruct%nat,fsw%dd_trans(1,1),1,minmode(1,1),1)
        call vcopy(3*runObj%atoms%astruct%nat,fsw%dd_trans(1,1),1,fsw%dd0_trans(1,1),1) 
        call daxpy(3*runObj%atoms%astruct%nat,tmp, minmode(1,1), 1, fsw%dd0_trans(1,1), 1 )
 
        !invert gradient in minmode direction
        !dd=dd-2.0_gp*ddot(3*nat,dd(1,1),1,minmode(1,1),1)*minmode
        tmp=2.0_gp*tmp
        call daxpy(3*runObj%atoms%astruct%nat,tmp, minmode(1,1), 1, fsw%dd_trans(1,1), 1 )
 
        tt=0.0_gp
        dt=0.0_gp
        maxd=-huge(1.0_gp)
        do iat=1,runObj%atoms%astruct%nat
            dt=fsw%dd_trans(1,iat)**2+fsw%dd_trans(2,iat)**2+fsw%dd_trans(3,iat)**2
            tt=tt+dt
            maxd=max(maxd,dt)
        enddo
        tt=sqrt(tt)
        maxd=sqrt(maxd)
 
        !trust radius approach
        if(maxd>uinp%saddle_trustr .or. (curv>=0.0_gp .and. fnrm<uinp%saddle_fnrmtol))then
!        if(maxd>trustr)then
            if(mhgpsst%iproc==0)call yaml_map('  (MHGPS) resize step ',maxd)
            scl=0.5_gp*uinp%saddle_trustr/maxd
            fsw%dd_trans=fsw%dd_trans*scl
            tt=tt*scl
            maxd=0.5_gp*uinp%saddle_trustr
        endif
        !do the move
        fsw%rxyz_trans(:,:,nhist)=fsw%rxyz_trans(:,:,nhist-1)-fsw%dd_trans(:,:)
        if (bigdft_get_geocode(runObj) == 'F') then
        call fixfrag_posvel(mhgpsst,bigdft_get_geocode(runObj),&
             runObj%atoms%astruct%nat,rcov,fsw%rxyz_trans(1,1,nhist),&
             tnatdmy,1,fixfragmented)
        if(fixfragmented .and. uinp%mhgps_verbosity >=2.and. mhgpsst%iproc==0)&
           call yaml_comment('fragmentation fixed')
        endif
        !displ=displ+tt
 
        fsw%delta_trans=fsw%rxyz_trans(:,:,nhist)-fsw%rxyzold_trans
        displ=displ+dnrm2(3*runObj%atoms%astruct%nat,fsw%delta_trans(1,1),1)
        runObj%inputs%inputPsiId=1
        call minenergyandforces(mhgpsst,.true.,uinp%imode,runObj,outs,&
             fsw%rxyz_trans(1,1,nhist),fsw%rxyzraw_trans(1,1,nhist),fsw%fxyz_trans(1,1,nhist),&
             fsw%fstretch_trans(1,1,nhist),fsw%fxyzraw_trans(1,1,nhist),etotp,iconnect,&
             nbond,fsw%wold_trans,uinp%saddle_alpha_stretch0,alpha_stretch)
        ener_count=ener_count+1.0_gp
        fsw%rxyzold_trans=fsw%rxyz_trans(:,:,nhist)
        detot=etotp-etotold
 
        call fnrmandforcemax(fsw%fxyzraw_trans(1,1,nhist),fnrm,fmax,runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        if (mhgpsst%iproc == 0 .and. uinp%mhgps_verbosity >=4) then
           fc=fc+1
           write(fn9,'(i9.9)') fc
           write(comment,'(a,1pe10.3,5x,1pe10.3)')&
           'ATTENTION! Forces below are no forces but tangents to '//&
           'the guessed reaction path| fnrm, fmax = ',fnrm,fmax
           call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//'_posout_'//&
                fn9,trim(comment),etotp,fsw%rxyz_trans(:,:,nhist),&
                forces=minmode)
        endif

        tmp=-ddot(3*runObj%atoms%astruct%nat,fsw%fxyz_trans(1,1,nhist),1,minmode(1,1),1)
        call vcopy(3*runObj%atoms%astruct%nat,fsw%fxyz_trans(1,1,nhist),1,fsw%ftmp_trans(1,1),1) 
        call daxpy(3*runObj%atoms%astruct%nat,tmp, minmode(1,1), 1, fsw%ftmp_trans(1,1), 1 )
        cosangle=-dot_double(3*runObj%atoms%astruct%nat,fsw%ftmp_trans(1,1),1,fsw%dd0_trans(1,1),1)/&
                 sqrt(dot_double(3*runObj%atoms%astruct%nat,fsw%ftmp_trans(1,1),1,fsw%ftmp_trans(1,1),1)*&
                 dot_double(3*runObj%atoms%astruct%nat,fsw%dd0_trans(1,1),1,fsw%dd0_trans(1,1),1))

        if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
            write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2),'//&
                    '1x,i3.3,1x,es12.5,1x,es9.2)')&
            '   (MHGPS) GEOPT ',nint(ener_count),it,etotp,detot,fmax,&
            fnrm, alpha,ndim,displ,alpha_stretch

        etot=etotp
        etotold=etot
        call convcheck_sad(fnrm,curv,0.0_gp,uinp%saddle_fnrmtol,icheck)
        if(icheck>icheckmax)then
            goto 1000
!        else if(icheck == 1)then
        else if(icheck == icheckmax)then
            tighten=.true.
        else if(icheck == 0)then
            tighten=.false.
        endif

        !now do step in hard directions
        if(uinp%imode==2)then
!           fsw%fstretch_trans(:,:,nhist)=fsw%fstretch_trans(:,:,nhist)-2.0_gp*&
!            ddot(3*nat,fsw%fstretch_trans(1,1,nhist),1,minmode(1,1),1)*minmode
            fsw%dds_trans=alpha_stretch*(fsw%fstretch_trans(:,:,nhist)-&
                2.0_gp*ddot(3*runObj%atoms%astruct%nat,fsw%fstretch_trans(1,1,nhist),1,&
                minmode(1,1),1)*minmode)
            dt=0.0_gp
            maxd=-huge(1.0_gp)
            do iat=1,runObj%atoms%astruct%nat
                !dt=fsw%fstretch_trans(1,iat,nhist)**2+&
                !fsw%fstretch_trans(2,iat,nhist)**2+fsw%fstretch_trans(3,iat,nhist)**2
                dt=fsw%dds_trans(1,iat)**2+fsw%dds_trans(2,iat)**2+fsw%dds_trans(3,iat)**2
                maxd=max(maxd,dt)
            enddo
            maxd=sqrt(maxd)

            !trust radius approach
            if(maxd>uinp%saddle_trustr)then
                if(mhgpsst%iproc==0)write(*,'(a,es10.3,1x,i0,1x,es10.3)')&
                    '(MHGPS) hard direction step too large:maxd,it,alpha_stretch',&
                    maxd,it,alpha_stretch
                scl=0.5_gp*uinp%saddle_trustr/maxd
                fsw%dds_trans=fsw%dds_trans*scl
            endif
            !fsw%rxyz_trans(:,:,nhist)=fsw%rxyz_trans(:,:,nhist)+alpha_stretch*fsw%fstretch_trans(:,:,nhist)
            fsw%rxyz_trans(:,:,nhist)=fsw%rxyz_trans(:,:,nhist)+fsw%dds_trans
!            if (bigdft_get_geocode(runObj) == 'F') then
!            call fixfrag_posvel(nat,rcov,fsw%rxyz_trans(1,1,nhist),tnatdmy,1,fixfragmented)
!                if(fixfragmented .and. uinp%mhgps_verbosity >=2.and. mhgpsst%iproc==0)&
!                  call yaml_comment('fragmentation fixed')
!            endif
        endif

        if (cosangle.gt..20_gp) then
            alpha=alpha*1.10_gp
        else
            alpha=max(alpha*.85_gp,uinp%saddle_alpha0_trans)
        endif

        call getSubSpaceEvecEval('(MHGPS)',mhgpsst%iproc,&
             uinp%mhgps_verbosity,runObj%atoms%astruct%nat,nhist,&
             uinp%saddle_nhistx_trans,ndim,uinp%saddle_cutoffratio,&
             fsw%lwork_trans,fsw%work_trans,fsw%rxyz_trans,&
             fsw%fxyz_trans,fsw%aa_trans,fsw%rr_trans,fsw%ff_trans,&
             fsw%rrr_trans,fsw%fff_trans,fsw%eval_trans,&
             fsw%res_trans,subspaceSucc)

!        fsw%delta_trans=fsw%rxyz_trans(:,:,nhist)-fsw%rxyz_trans(:,:,nhist-1)
!        displ=displ+dnrm2(3*nat,fsw%delta_trans(1,1),1)
  enddo

  if(mhgpsst%iproc==0)call yaml_warning('(MHGPS) No convergence in findsad')
!stop 'no convergence in findsad'
    return

1000 continue
    converged=.true.
    etot=etotp
    if(mhgpsst%iproc==0)call yaml_map( "(MHGPS) convergence at",nint(ener_count))
    
    do iat=1,runObj%atoms%astruct%nat
        do i=1,3
            wpos(i,iat)= fsw%rxyz_trans(i,iat,nhist)
            fout(i,iat)= fsw%fxyzraw_trans(i,iat,nhist)
        enddo
    enddo

end subroutine
!=====================================================================
subroutine opt_curv(mhgpsst,itgeopt,imode,fsw,uinp,runObj,outs,alpha0,curvforcediff,nit,nhistx,rxyz_fix,&
                    fxyz_fix,dxyzin,curv,fout,fnrmtol,ener_count,&
                    converged,iconnect,nbond,alpha_stretch0,&
                    maxcurvrise,cutoffratio,minoverlap)!,mode)
    use module_base
    use yaml_output
    use module_sqn
    use bigdft_run, only: run_objects, state_properties
    use module_userinput
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    integer, intent(in)        :: itgeopt
    integer, intent(in)        :: imode
    type(userinput), intent(in):: uinp
    type(findsad_work), intent(inout)      :: fsw
    integer, intent(in)        :: nbond
    integer, intent(in)        :: iconnect(2,nbond)
    real(gp),intent(in)        :: alpha_stretch0
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    real(gp), intent(in)       :: maxcurvrise,cutoffratio,minoverlap
    integer, intent(in)        :: nit,nhistx
    real(gp), intent(in)       :: alpha0,curvforcediff
    real(gp), dimension(3,runObj%atoms%astruct%nat) :: dxyzin,fout,rxyz_fix,fxyz_fix!,mode
    logical, intent(out)       :: converged
    logical                    :: steep
    integer                    :: iat,l,itswitch
    integer                    :: ihist,it
    real(gp)                   :: ener_count,curv
    real(gp)                   :: fnrmtol,curvold,fnrm,curvp,fmax
    real(gp)                   :: dcurv,tt,cosangle
    real(gp)                   :: overlap
    logical                    :: subspaceSucc
    real(gp), dimension(3,runObj%atoms%astruct%nat) :: dxyzin0
    !internal
    integer, parameter         :: nmaxrise=10
    integer                    :: nrise, itup
    real(gp), dimension(3,runObj%atoms%astruct%nat) :: rxyzOld,delta
    real(gp)                   :: displr,displp
real(gp) :: alpha0int
    !functions
    real(gp)                   :: ddot
    real(gp)                   :: dnrm2

    if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)write(*,'(a,1x,es9.2)')'   (MHGPS) CUOPT minoverlap',minoverlap

    if(.not.uinp%share_rot_history)then
        fsw%alpha_stretch_rot=alpha_stretch0
        fsw%ndim_rot=0
        fsw%nhist_rot=0
        fsw%alpha_rot=alpha0
    endif
    alpha0int=alpha0

    nrise=0
    itup=nint(ener_count)
    converged =.false.
    subspaceSucc=.true.
    displr=0.0_gp
    displp=0.0_gp
    dcurv=0.0_gp
    fsw%wold_rot=0.0_gp
    if(fsw%nhist_rot==0)then
        do iat=1,runObj%atoms%astruct%nat
           do l=1,3
              fsw%rxyz_rot(l,iat,fsw%nhist_rot)=dxyzin(l,iat)
              rxyzOld(l,iat)=dxyzin(l,iat)
              dxyzin0(l,iat)=dxyzin(l,iat)
           enddo
        enddo
    endif

    call mincurvforce(mhgpsst,imode,runObj,outs,curvforcediff,&
         rxyz_fix(1,1),fxyz_fix(1,1),fsw%rxyz_rot(1,1,fsw%nhist_rot),&
         fsw%rxyzraw_rot(1,1,fsw%nhist_rot),fsw%fxyz_rot(1,1,fsw%nhist_rot),fsw%fstretch_rot(1,1,fsw%nhist_rot),&
         fsw%fxyzraw_rot(1,1,fsw%nhist_rot),curv,1,ener_count,iconnect,nbond,fsw%wold_rot,&
         alpha_stretch0,fsw%alpha_stretch_rot)
    if(imode==2)then
        fsw%rxyz_rot(:,:,fsw%nhist_rot)=fsw%rxyz_rot(:,:,fsw%nhist_rot)+&
                        fsw%alpha_stretch_rot*fsw%fstretch_rot(:,:,fsw%nhist_rot)
    endif

    call fnrmandforcemax(fsw%fxyzraw_rot(1,1,fsw%nhist_rot),fnrm,fmax,runObj%atoms%astruct%nat)
    fnrm=sqrt(fnrm)
    curvold=curv
    curvp=curv
    overlap=ddot(3*runObj%atoms%astruct%nat,dxyzin0(1,1),1,fsw%rxyz_rot(1,1,fsw%nhist_rot),1)
 
    if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
     write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2),1x,i3.3,'//&
          '2(1x,es9.2),2(1x,es12.5))')'   (MHGPS) CUOPT ',&
          nint(ener_count),0,curvp,dcurv,fmax,fnrm, fsw%alpha_rot,&
          fsw%ndim_rot,fsw%alpha_stretch_rot,overlap,displr,displp
    itswitch=-2
    minloop: do it=1,nit
        fsw%nhist_rot=fsw%nhist_rot+1
 
        if ((.not. subspaceSucc) .or. fnrm.gt.uinp%saddle_steepthresh_rot  .or. it.le.itswitch ) then
            fsw%ndim_rot=0
            steep=.true.
            if (it.gt.itswitch) itswitch=it+nhistx
        else
            steep=.false.
        endif

        ! make space in the history list
        if (fsw%nhist_rot.gt.nhistx) then
            fsw%nhist_rot=nhistx
            do ihist=0,fsw%nhist_rot-1
                do iat=1,runObj%atoms%astruct%nat
                    do l=1,3
                        fsw%rxyz_rot(l,iat,ihist)=fsw%rxyz_rot(l,iat,ihist+1)
                        fsw%fxyz_rot(l,iat,ihist)=fsw%fxyz_rot(l,iat,ihist+1)
                        fsw%rxyzraw_rot(l,iat,ihist)=fsw%rxyzraw_rot(l,iat,ihist+1)
                        fsw%fxyzraw_rot(l,iat,ihist)=fsw%fxyzraw_rot(l,iat,ihist+1)
                        fsw%fstretch_rot(l,iat,ihist)=fsw%fstretch_rot(l,iat,ihist+1)
                    enddo
                enddo
            enddo
        endif
        500 continue
        call modify_gradient(runObj%atoms%astruct%nat,fsw%ndim_rot,fsw%rrr_rot(1,1,1),fsw%eval_rot(1),fsw%res_rot(1),&
             fsw%fxyz_rot(1,1,fsw%nhist_rot-1),fsw%alpha_rot,fsw%dd_rot(1,1))

        tt=0.0_gp
        do iat=1,runObj%atoms%astruct%nat
            do l=1,3
                tt=tt+fsw%dd_rot(l,iat)**2
            enddo
        enddo
!        displ=displ+sqrt(tt)

        do iat=1,runObj%atoms%astruct%nat
            fsw%rxyz_rot(1,iat,fsw%nhist_rot)=fsw%rxyz_rot(1,iat,fsw%nhist_rot-1)-fsw%dd_rot(1,iat)
            fsw%rxyz_rot(2,iat,fsw%nhist_rot)=fsw%rxyz_rot(2,iat,fsw%nhist_rot-1)-fsw%dd_rot(2,iat)
            fsw%rxyz_rot(3,iat,fsw%nhist_rot)=fsw%rxyz_rot(3,iat,fsw%nhist_rot-1)-fsw%dd_rot(3,iat)
        enddo

        delta=fsw%rxyz_rot(:,:,fsw%nhist_rot)-rxyzOld
        displr=displr+dnrm2(3*runObj%atoms%astruct%nat,delta(1,1),1)
        call mincurvforce(mhgpsst,imode,runObj,outs,curvforcediff,&
             rxyz_fix(1,1),fxyz_fix(1,1),fsw%rxyz_rot(1,1,fsw%nhist_rot),&
             fsw%rxyzraw_rot(1,1,fsw%nhist_rot),fsw%fxyz_rot(1,1,fsw%nhist_rot),fsw%fstretch_rot(1,1,fsw%nhist_rot),&
             fsw%fxyzraw_rot(1,1,fsw%nhist_rot),curvp,1,ener_count,iconnect,nbond,&
             fsw%wold_rot,alpha_stretch0,fsw%alpha_stretch_rot)
        dcurv=curvp-curvold

        call fnrmandforcemax(fsw%fxyzraw_rot(1,1,fsw%nhist_rot),fnrm,fmax,runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        cosangle=-dot_double(3*runObj%atoms%astruct%nat,fsw%fxyz_rot(1,1,fsw%nhist_rot),1,fsw%dd_rot(1,1),1)/&
                  sqrt(dot_double(3*runObj%atoms%astruct%nat,fsw%fxyz_rot(1,1,fsw%nhist_rot),1,&
                  fsw%fxyz_rot(1,1,fsw%nhist_rot),1)*&
                  dot_double(3*runObj%atoms%astruct%nat,fsw%dd_rot(1,1),1,fsw%dd_rot(1,1),1))

        if (dcurv.gt.maxcurvrise .and. fsw%alpha_rot>1.e-1_gp*alpha0int) then 
            itup=nint(ener_count)
            nrise=nrise+1
            if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity>=3)&
                call yaml_comment('INFO: (MHGPS) Curv. raised by'//&
                     ' more than maxcurvrise '//trim(yaml_toa(it))//&
                     ''//trim(yaml_toa(dcurv)))
            overlap=ddot(3*runObj%atoms%astruct%nat,dxyzin0(1,1),1,fsw%rxyz_rot(1,1,fsw%nhist_rot),1)
            if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
                write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2)'//&
                      ',1x,i3.3,2(1x,es9.2),2(1x,es12.5))')&
                      '   (MHGPS) CUOPT ',nint(ener_count),it,curvp,&
                      dcurv,fmax,fnrm, fsw%alpha_rot,fsw%ndim_rot,&
                      fsw%alpha_stretch_rot,overlap,displr,displp
            fsw%alpha_rot=.5_gp*fsw%alpha_rot
            alpha0int=alpha0*0.1_gp
            if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity>=3)&
                call yaml_comment('INFO: (MHGPS) alpha reset (opt. curv): '//&
                     trim(yaml_toa(fsw%alpha_rot)))
            fsw%ndim_rot=0
            if((trim(adjustl(char(runObj%run_mode)))=='QM_RUN_MODE') .and. (runObj%inputs%inputPsiId/=0))then
                if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity>=3)&
                    call yaml_comment('INFO: (MHGPS) Will use LCAO input guess from now on '//&
                    '(until end of current minmode optimization).')
                runObj%inputs%inputPsiId=0
                call mincurvforce(mhgpsst,imode,runObj,outs,&
                     curvforcediff,rxyz_fix(1,1),fxyz_fix(1,1),&
                     fsw%rxyz_rot(1,1,fsw%nhist_rot-1),&
                     fsw%rxyzraw_rot(1,1,fsw%nhist_rot-1),&
                     fsw%fxyz_rot(1,1,fsw%nhist_rot-1),&
                     fsw%fstretch_rot(1,1,fsw%nhist_rot-1),&
                     fsw%fxyzraw_rot(1,1,fsw%nhist_rot-1),curvold,1,&
                     ener_count,iconnect,nbond,fsw%wold_rot,&
                     alpha_stretch0,fsw%alpha_stretch_rot)
                if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
                 write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2)'//&
                      ',1x,i3.3,1x,es9.2,2(1x,es12.5))')&
                      '   (MHGPS)1 CUOPT ',nint(ener_count),it,curvp,&
                      dcurv,fmax,fnrm, fsw%alpha_rot,fsw%ndim_rot,&
                      fsw%alpha_stretch_rot,displr,displp
            endif

            do iat=1,runObj%atoms%astruct%nat
                fsw%rxyz_rot(1,iat,0)=fsw%rxyzraw_rot(1,iat,fsw%nhist_rot-1)
                fsw%rxyz_rot(2,iat,0)=fsw%rxyzraw_rot(2,iat,fsw%nhist_rot-1)
                fsw%rxyz_rot(3,iat,0)=fsw%rxyzraw_rot(3,iat,fsw%nhist_rot-1)
                fsw%rxyzraw_rot(1,iat,0)=fsw%rxyzraw_rot(1,iat,fsw%nhist_rot-1)
                fsw%rxyzraw_rot(2,iat,0)=fsw%rxyzraw_rot(2,iat,fsw%nhist_rot-1)
                fsw%rxyzraw_rot(3,iat,0)=fsw%rxyzraw_rot(3,iat,fsw%nhist_rot-1)
 
                fsw%fxyz_rot(1,iat,0)=fsw%fxyzraw_rot(1,iat,fsw%nhist_rot-1)
                fsw%fxyz_rot(2,iat,0)=fsw%fxyzraw_rot(2,iat,fsw%nhist_rot-1)
                fsw%fxyz_rot(3,iat,0)=fsw%fxyzraw_rot(3,iat,fsw%nhist_rot-1)
                fsw%fxyzraw_rot(1,iat,0)=fsw%fxyzraw_rot(1,iat,fsw%nhist_rot-1)
                fsw%fxyzraw_rot(2,iat,0)=fsw%fxyzraw_rot(2,iat,fsw%nhist_rot-1)
                fsw%fxyzraw_rot(3,iat,0)=fsw%fxyzraw_rot(3,iat,fsw%nhist_rot-1)
            enddo
            fsw%nhist_rot=1
            goto  500
        endif
        if (dcurv.gt.maxcurvrise) then 
            itup=nint(ener_count)
            nrise=nrise+1
        else if(nint(ener_count)-itup>5)then
            nrise=0
        endif
        curv=curvp
        curvold=curv
    
        delta=fsw%rxyz_rot(:,:,fsw%nhist_rot)-rxyzOld
        displp=displp+dnrm2(3*runObj%atoms%astruct%nat,delta(1,1),1)
        rxyzOld=fsw%rxyz_rot(:,:,fsw%nhist_rot)
        overlap=ddot(3*runObj%atoms%astruct%nat,dxyzin0(1,1),1,fsw%rxyz_rot(1,1,fsw%nhist_rot),1)
        if(mhgpsst%iproc==0.and.uinp%mhgps_verbosity>=2)&
            write(*,'(a,1x,i4.4,1x,i4.4,1x,es21.14,4(1x,es9.2),1x,'//&
                 'i3.3,2(1x,es9.2),2(1x,es12.5))')&
                 '   (MHGPS) CUOPT ',nint(ener_count),it,curvp,dcurv,&
                 fmax,fnrm, fsw%alpha_rot,fsw%ndim_rot,&
                 fsw%alpha_stretch_rot,overlap,displr,displp

        do iat=1,runObj%atoms%astruct%nat
            do l=1,3
                dxyzin(l,iat)= fsw%rxyz_rot(l,iat,fsw%nhist_rot) !to be done before stretch modification
            enddo
        enddo

        if(imode==2)then
            fsw%rxyz_rot(:,:,fsw%nhist_rot)=&
                fsw%rxyz_rot(:,:,fsw%nhist_rot)+&
                fsw%alpha_stretch_rot*&
                fsw%fstretch_rot(:,:,fsw%nhist_rot)
        endif

!        if (cosangle.gt..20_gp) then
        if (cosangle.gt..40_gp) then
            fsw%alpha_rot=fsw%alpha_rot*1.10_gp
        else
            fsw%alpha_rot=max(fsw%alpha_rot*.85_gp,alpha0int)
        endif


        call getSubSpaceEvecEval('(MHGPS)',mhgpsst%iproc,&
             uinp%mhgps_verbosity,runObj%atoms%astruct%nat,&
             fsw%nhist_rot,nhistx,fsw%ndim_rot,cutoffratio,&
             fsw%lwork_rot,fsw%work_rot,fsw%rxyz_rot,fsw%fxyz_rot,&
             fsw%aa_rot,fsw%rr_rot,fsw%ff_rot,fsw%rrr_rot,&
             fsw%fff_rot,fsw%eval_rot,fsw%res_rot,subspaceSucc)
        if(.not.subspaceSucc)stop 'subroutine findsad: no success in getSubSpaceEvecEval.'
!        if (fnrm.le.fnrmtol) goto 1000 !has to be in this line for shared history case
        if (fnrm.le.fnrmtol.or.(overlap<minoverlap.and.itgeopt>1)) goto 1000 !has to be in this line for shared history case
        if(nrise>nmaxrise)then
            if(mhgpsst%iproc==0)then
                call yaml_warning('(MHGPS) opt_curv failed because'//&
                     ' nrise > nmaxrise. Wrong finite difference'//&
                     ' or errors in energies and forces.')
                if(trim(adjustl(char(runObj%run_mode)))=='QM_RUN_MODE')then
                    call yaml_warning('(MHGPS) Did QM method converge?')
                endif
            endif
            exit minloop
        endif
    enddo minloop
!    write(*,*) "No convergence in optcurv"
    if(mhgpsst%iproc==0)call yaml_warning('(MHGPS) No convergence in optcurv.')
!stop 'no convergence in optcurv'
    return
    1000 continue
    converged=.true.
    do iat=1,runObj%atoms%astruct%nat
        do l=1,3
!            dxyzin(l,iat)= fsw%rxyz_rot(l,iat,fsw%nhist_rot)
            fout(l,iat)= fsw%fxyzraw_rot(l,iat,fsw%nhist_rot)
        enddo
    enddo
    curv=curvp
end subroutine
!=====================================================================
!> Computes the (curvature along vec) = vec^t H vec / (vec^t*vec)
!! vec mus be normalized
subroutine curvforce(mhgpsst,runObj,outs,diff,rxyz1,fxyz1,vec,curv,rotforce,imethod,ener_count)
    use module_base
    use yaml_output
    use bigdft_run, only: bigdft_get_geocode
    use module_energyandforces
    use bigdft_run, only: run_objects, state_properties
    use module_mhgps_state
    use module_forces
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    integer, intent(in) :: imethod
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    real(gp), intent(in) :: rxyz1(3,runObj%atoms%astruct%nat),fxyz1(3,runObj%atoms%astruct%nat),diff
    real(gp), intent(inout) :: vec(3,runObj%atoms%astruct%nat)
    real(gp), intent(out) :: curv
    real(gp), intent(out) :: rotforce(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout) :: ener_count
    !internal
    integer :: infocode
    real(gp) :: diffinv, etot2,fnoise
    real(gp),allocatable :: rxyz2(:,:), fxyz2(:,:)
    real(gp),allocatable :: drxyz(:,:), dfxyz(:,:)
real(gp) :: maxd(runObj%atoms%astruct%nat)
integer :: i
    !functions
    real(gp), external :: dnrm2,ddot

    
!    diff=1.e-3_gp !lennard jones

    diffinv=1.0_gp/(diff)
    rxyz2 = f_malloc((/ 1.to.3, 1.to.runObj%atoms%astruct%nat/),id='rxyz2')
    fxyz2 = f_malloc((/ 1.to.3, 1.to.runObj%atoms%astruct%nat/),id='fxyz2')
    drxyz = f_malloc((/ 1.to.3, 1.to.runObj%atoms%astruct%nat/),id='drxyz')
    dfxyz = f_malloc((/ 1.to.3, 1.to.runObj%atoms%astruct%nat/),id='dfxyz')
    call elim_moment_fs(runObj%atoms%astruct%nat,vec(1,1))
    if (bigdft_get_geocode(runObj)=='F')&
        call elim_torque_bastian(runObj%atoms%astruct%nat,&
             rxyz1(1,1),vec(1,1))
    

    call clean_forces_base(runObj%atoms,vec) 
    vec = vec / dnrm2(3*runObj%atoms%astruct%nat,vec(1,1),1)
    rxyz2 = rxyz1 + diff * vec

    call mhgpsenergyandforces(mhgpsst,runobj,outs,rxyz2(1,1),&
         fxyz2(1,1),fnoise,etot2,infocode)
    ener_count=ener_count+1.0_gp

    if(imethod==1)then
        dfxyz = fxyz2(:,:)-fxyz1(:,:)
        drxyz = rxyz2(:,:)-rxyz1(:,:)
        drxyz = drxyz * diffinv
        curv  = - diffinv * ddot(3*runObj%atoms%astruct%nat,dfxyz(1,1),1,drxyz(1,1),1)

        rotforce = 2.0_gp*dfxyz*diffinv + 2.0_gp * curv * drxyz
        call elim_moment_fs(runObj%atoms%astruct%nat,rotforce(1,1))
        if (bigdft_get_geocode(runObj)=='F')&
            call elim_torque_bastian(runObj%atoms%astruct%nat,&
                 rxyz1(1,1),rotforce(1,1))
    else
        stop 'unknown method for curvature computation'
    endif
        rotforce = rotforce - ddot(3*runObj%atoms%astruct%nat,rotforce(1,1),1,vec(1,1),1)*vec
    call f_free(rxyz2)
    call f_free(fxyz2)
    call f_free(drxyz)
    call f_free(dfxyz)
    
end subroutine
!=====================================================================
subroutine fixfrag_posvel(mhgpsst,geocode,nat,rcov,pos,vel,option,occured)
!UNTERSCHIED ZUR MH ROUTINE:
!pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.0_gp*bondlength)/mindist)
!ANSTATT
!pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.5_gp*bondlength)/mindist)
!MACHT GROSSEN PERFORMANCE UNTERSCHIED AUS
!

!This subroutine can perform two tasks.
!ATTENTION: it will only work on free BC!!!
!
!option=1
!The atoms in pos are analyzed and, if there is a fragmentation occuring, the
!main fragment will be identified and all neighboring fragments will be moved towards the nearest
!atom of the main fragment. The array pos will then be updated and returned.
!
!option=2
!The fragments are identified and the center of mass of all fragments are computed separately.
!The center of mass of all cluster is also computed.
!Then, the velocities are modified in such a way that the projection of the velocities 
!along the vector pointing towards the center of mass of all fragments are inverted 
!!use module_base
!!use module_types
!!use m_ab6_symmetry
use module_base
use yaml_output
use module_mhgps_state
implicit none
integer, intent(in) :: nat
!type(atoms_data), intent(in) :: at
type(mhgps_state), intent(in) :: mhgpsst
real(gp),dimension(3,nat), INTENT(INOUT) :: pos
real(gp),dimension(3,nat), INTENT(INOUT) :: vel
real(gp),dimension(nat), INTENT(IN) :: rcov
character(len=*), intent(in) :: geocode
integer, INTENT(IN):: option
integer :: nfrag, nfragold
logical :: occured,niter
real(gp)::  dist, mindist, angle, vec(3), cmass(3), velcm(3), bondlength, bfactor,rnrmi,scpr
real(gp):: ekin,vcm1,vcm2,vcm3,ekin0,scale
real(gp), allocatable:: cm_frags(:,:), vel_frags(:,:)
integer::iat, jat, natfragx(1), imin(2),ifrag
integer, allocatable:: fragcount(:)
integer, allocatable:: nat_frags(:)
integer, dimension(nat):: fragarr
logical, allocatable:: invert(:)

!The bondlength (in atomic units) is read from file input.bondcut
! OPTION 1: System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           The two fragment are then brought together such that the minimal distance equals 1.5*bondlength
! OPTION  : System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           the velocities are then inverted
!open(unit=43,file="input.bondcut")
!read(43,*) bondlength
!close(43)

if (option == 1) then 
   bfactor=1.5_gp
else if (option == 2) then 
   bfactor=2.0_gp
else
   stop 'wrong option'
endif



fragarr(:)=0                     !Array, which atom belongs to which fragment
nfrag=0                       !Number of fragments

!Check for correct input
if (option.ne.1 .and. option.ne.2) stop "Wrong option in fixfrag_refvels"

!Calculate number of fragments and fragmentlist of the atoms
loop_nfrag: do
   nfragold=nfrag
   do iat=1,nat                !Check the first atom that isn't part of a cluster yet
      if(fragarr(iat)==0) then
         nfrag=nfrag+1
         fragarr(iat)=nfrag
         exit 
      endif
   enddo
   if (nfragold==nfrag) exit loop_nfrag
7000 niter=.false.
   do iat=1,nat                !Check if all the other atoms are part of the current cluster
      do jat=1,nat
         bondlength=rcov(iat)+rcov(jat)
         if(nfrag==fragarr(iat) .AND. jat.ne.iat .AND. fragarr(jat)==0) then
            dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
            if(dist<(bfactor*bondlength)**2) then
               fragarr(jat)=nfrag
               niter=.true.
            endif
         endif
      enddo
   enddo
   if(niter) then
      goto 7000
   endif
end do loop_nfrag


!   if(mhgpsst%iproc==0) write(*,*) '(MH) nfrag=',nfrag
occured=.false.
if(nfrag.ne.1) then          !"if there is fragmentation..."
   occured=.true.
   if(mhgpsst%iproc==0) then
      call yaml_mapping_open('(MHGPS) FIX')
      call yaml_map('(MHGPS) Number of Fragments counted with option', (/nfrag,option/))
   endif
   if (option==1) then !OPTION=1, FIX FRAGMENTATION
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."

      !Find out which fragment is the main cluster
      fragcount = f_malloc(nfrag,id='fragcount')
      fragcount=0
      do ifrag=1,nfrag
         do iat=1,nat
            if(fragarr(iat)==ifrag) then
               fragcount(ifrag)=fragcount(ifrag)+1
            endif
         enddo
      enddo
      natfragx=maxloc(fragcount(:))
      !if(mhgpsst%iproc==0) call yaml_map('(MH) The main Fragment index is', natfragx(1))

      !Find the minimum distance between the clusters
      do ifrag=1,nfrag
         mindist=1.e100_gp
         if(ifrag.ne.natfragx(1)) then
            do iat=1,nat
               if(fragarr(iat)==ifrag) then
                  do jat=1,nat
                     if(fragarr(jat)==natfragx(1)) then
                        dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
                        if(dist<mindist**2) then
                           mindist=sqrt(dist)
                           imin(1)=jat  !Atom with minimal distance in main fragment
                           imin(2)=iat   !Atom with minimal distance in fragment ifrag
                        endif
                     endif
                  enddo
               endif
            enddo

!            if (mhgpsst%iproc == 0) then
!               write(444,*) nat, 'atomic '
!               write(444,*) 'A fragmented configuration ',imin(1),imin(2)
!               do iat=1,nat
!                  write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
!               enddo
!            endif


            vec(:)=pos(:,imin(1))-pos(:,imin(2))
            bondlength=rcov(imin(1))+rcov(imin(2))
            do iat=1,nat        !Move fragments back towards the main fragment 
               if(fragarr(iat)==ifrag) then
                  pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.0_gp*bondlength)/mindist)
                  fragarr(iat)=natfragx(1)
               endif
            enddo


         endif
      enddo
      call f_free(fragcount)
      if(mhgpsst%iproc==0) then
         call yaml_comment('(MHGPS) FIX: Fragmentation fixed!')
         call yaml_mapping_close()
      end if
!      if (mhgpsst%iproc == 0) then
!         write(444,*) nat, 'atomic '
!         write(444,*) ' fixed configuration '
!         do iat=1,nat
!            write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
!         enddo
!      endif

      !   endif
   elseif(option==2) then !OPTION=2, INVERT VELOCITIES
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."
!      if(mhgpsst%iproc==0) call yaml_map('(MH) FIX: Preparing to invert velocities, option:',option)
      !Compute center of mass of all fragments and the collectiove velocity of each fragment
      cm_frags = f_malloc((/ 3, nfrag /),id='cm_frags')
      vel_frags = f_malloc((/ 3, nfrag /),id='vel_frags')
      nat_frags = f_malloc(nfrag,id='nat_frags')
      invert = f_malloc(nfrag,id='invert')

      cm_frags(:,:)=0.0_gp
      vel_frags(:,:)=0.0_gp
      nat_frags(:)=0         !number of atoms per fragment
      cmass(:)=0.0_gp
      velcm(:)=0.0_gp
      do iat=1,nat
         ifrag=fragarr(iat)
         nat_frags(ifrag)=nat_frags(ifrag)+1
         cm_frags(:,ifrag)=cm_frags(:,ifrag)+pos(:,iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo

      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)/real(nat_frags(ifrag),8)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)/real(nat_frags(ifrag),8)
         cmass(:)=cmass(:)+cm_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
         velcm(:)=velcm(:)+vel_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
      enddo
!      if (mhgpsst%iproc==0) call yaml_map('(MH) CM VELOCITY',sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2))
!      if (velcm(1)**2+velcm(2)**2+velcm(3)**2.gt.1.e-24_gp) then
!         if (mhgpsst%iproc==0) call yaml_comment('(MH) NONZERO CM VELOCITY')
!      endif


      ! now cm_frags contains the unit vector pointing from the center of mass of the entire system to the center of mass of the fragment
      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)-cmass(:)
         rnrmi=1.0_gp/sqrt(cm_frags(1,ifrag)**2+cm_frags(2,ifrag)**2+cm_frags(3,ifrag)**2)
         cm_frags(1,ifrag)=cm_frags(1,ifrag)*rnrmi
         cm_frags(2,ifrag)=cm_frags(2,ifrag)*rnrmi
         cm_frags(3,ifrag)=cm_frags(3,ifrag)*rnrmi
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.0_gp/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
         if (angle.gt.0.0_gp) then
            invert(ifrag)=.true.
         else
            invert(ifrag)=.false.
         endif
!         if (mhgpsst%iproc==0) then
!           write(*,*) '(MH) ifrag, angle ',ifrag, angle,invert(ifrag)
!           call yaml_mapping_open('(MH) Frag. Info',flow=.true.)
!            call yaml_map('ifrag',ifrag)
!            call yaml_map('angle',angle)
!            call yaml_map('ifrag inverted',invert(ifrag))
!           call yaml_mapping_close(advance='yes')
!         endif
      enddo
      !Decompose each atomic velocity into an component parallel and perpendicular to the cm_frags  vector and inter the 
      !paralle part if it point away from the CM

      !Check kinetic energy before inversion
      ekin0=0.0_gp
      vcm1=0.0_gp
      vcm2=0.0_gp
      vcm3=0.0_gp
      do iat=1,nat
         ekin0=ekin0+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
!      if (mhgpsst%iproc==0) then
!          write(*,'(a,e14.7,3(e10.3))') '(MH) EKIN CM before invert',ekin0,vcm1,vcm2,vcm3
!!          call yaml_mapping_open(,flow=.true.)
!          call yaml_map('(MH) EKIN CM before invert',(/ekin0,vcm1,vcm2,vcm3/),fmt='(e10.3)')
!!          call yaml_mapping_close(advance='yes')
!      endif
      !Checkend kinetic energy before inversion

      do iat=1,nat
         ! inversions  by fragment group
         ifrag=fragarr(iat)
         if (invert(ifrag)) then
            scpr=cm_frags(1,ifrag)*vel(1,iat)+cm_frags(2,ifrag)*vel(2,iat)+cm_frags(3,ifrag)*vel(3,iat)
            vel(:,iat)=vel(:,iat)-scpr*cm_frags(:,ifrag)*2.0_gp
         endif
      enddo

      call elim_moment_fs(nat,vel)
      if (trim(adjustl(geocode))=='F')&
          call elim_torque_bastian(nat,pos,vel)

      ! scale velocities to regain initial ekin0
      ekin=0.0_gp
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
      enddo
      scale=sqrt(ekin0/ekin)
      do iat=1,nat
         vel(1,iat)=vel(1,iat)*scale
         vel(2,iat)=vel(2,iat)*scale
         vel(3,iat)=vel(3,iat)*scale
      enddo

      !Check kinetic energy after inversion
      ekin=0.0_gp
      vcm1=0.0_gp
      vcm2=0.0_gp
      vcm3=0.0_gp
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
!      if (mhgpsst%iproc==0) then
!          write(*,'(a,e14.7,3(e10.3))') '(MH) EKIN CM after  invert',ekin,vcm1,vcm2,vcm3
!          !call yaml_mapping_open('(MH) EKIN CM after invert',flow=.true.)
!          call yaml_map('(MH) EKIN CM after invert',(/ekin0,vcm1,vcm2,vcm3/),fmt='(e10.3)')
!          !call yaml_mapping_close(advance='yes')
!      endif
      !Checkend kinetic energy after inversion

      !Check angle  after inversion
      vel_frags(:,:)=0.0_gp
      do iat=1,nat
         ifrag=fragarr(iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo
      do ifrag=1,nfrag
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.0_gp/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
      enddo
      !Checkend kinetic energy after inversion

      call f_free(cm_frags)
      call f_free(vel_frags)
      call f_free(nat_frags)
      call f_free(invert)
   endif
endif
end subroutine fixfrag_posvel
!=====================================================================
subroutine mincurvforce(mhgpsst,imode,runObj,outs,diff,rxyz1,fxyz1,&
           vec,vecraw,rotforce,rotfstretch,rotforceraw,curv,imethod,&
           ec,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch)
    use module_base, only: gp
    use module_sqn
    use bigdft_run, only: run_objects, state_properties
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    integer,  intent(in)     :: imode
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    integer,  intent(in)     :: imethod
    integer,  intent(in)     :: nbond_
    real(gp), intent(in)     :: diff
    real(gp), intent(in)     :: rxyz1(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)     :: fxyz1(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout)  :: vec(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout)  :: vecraw(3,runObj%atoms%astruct%nat)
    real(gp), intent(out)    :: rotforce(3,runObj%atoms%astruct%nat)
    real(gp), intent(out)    :: rotforceraw(3,runObj%atoms%astruct%nat)
    real(gp), intent(out)    :: rotfstretch(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout)  :: wold(nbond_)
    real(gp), intent(out)    :: curv
    real(gp), intent(inout)  :: ec
    real(gp), intent(in)     :: alpha_stretch0
    real(gp), intent(inout)  :: alpha_stretch
    integer,  intent(in)     :: iconnect(2,nbond_)
    !internal
    integer :: infocode
    real(gp) :: rxyz2(3,runObj%atoms%astruct%nat)

     vecraw=vec
!    call mhgpsenergyandforces(nat,rat,fat,epot)
     call curvforce(mhgpsst,runObj,outs,diff,rxyz1,fxyz1,vec,curv,&
          rotforce,imethod,ec)
     rxyz2 =rxyz1+diff*vec
     rotforceraw=rotforce
     rotfstretch=0.0_gp

     if(imode==2)then
         call projectbond(runObj%atoms%astruct%nat,nbond_,rxyz2,rotforce,rotfstretch,&
              iconnect,wold,alpha_stretch0,alpha_stretch)
     endif

end subroutine mincurvforce
!=====================================================================
subroutine minenergyandforces(mhgpsst,eeval,imode,runObj,outs,rat,&
           rxyzraw,fat,fstretch,fxyzraw,epot,iconnect,nbond_,wold,&
           alpha_stretch0,alpha_stretch)
    use module_base, only: gp
    use module_energyandforces
    use module_sqn
    use bigdft_run, only: run_objects, state_properties
    use module_mhgps_state
    implicit none
    !parameter
    type(mhgps_state), intent(inout) :: mhgpsst
    integer, intent(in)           :: imode
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    integer, intent(in)           :: nbond_
    integer, intent(in)           :: iconnect(2,nbond_)
    real(gp),intent(in)           :: rat(3,runObj%atoms%astruct%nat)
    real(gp),intent(out)          :: rxyzraw(3,runObj%atoms%astruct%nat)
    real(gp),intent(out)          :: fxyzraw(3,runObj%atoms%astruct%nat)
    real(gp),intent(inout)          :: fat(3,runObj%atoms%astruct%nat)
    real(gp),intent(out)          :: fstretch(3,runObj%atoms%astruct%nat)
    real(gp), intent(inout)       :: wold(nbond_)
    real(gp), intent(in)          :: alpha_stretch0
    real(gp), intent(inout)       :: alpha_stretch
    real(gp), intent(inout)         :: epot
    logical, intent(in)           :: eeval
    !internal
    integer :: infocode
    real(gp) :: fnoise

    rxyzraw=rat
    if(eeval)call mhgpsenergyandforces(mhgpsst,runObj,outs,rat,fat,&
                  fnoise,epot,infocode)
    fxyzraw=fat
    fstretch=0.0_gp

    if(imode==2)then
        call projectbond(runObj%atoms%astruct%nat,nbond_,rat,fat,fstretch,iconnect,&
             wold,alpha_stretch0,alpha_stretch)
    endif

end subroutine minenergyandforces
!=====================================================================
subroutine convcheck_sad(fmax,curv,fluctfrac_fluct,forcemax,check)
  use module_base
  implicit none
  real(gp), intent(in):: fmax, fluctfrac_fluct,forcemax,curv
  integer, intent(inout)::check

  if ( fmax < max(forcemax,fluctfrac_fluct).and. curv <0.0_gp) then
    check=check+1
  else
    check=0
  endif
end subroutine convcheck_sad
!=====================================================================
subroutine elim_torque_bastian(nat,ratin,vat)
use module_base
!theory:
!(in the following: x represents cross product)
!L = sum(r_i x p_i) = I w
!w = I^-1 L
!v_ortho = w x r
!v_tot = v_ortho + v_radial
!=> vradial = v_tot - w x r
!note: routine must be modified before using for cases
!where masses /= 1
    implicit none
    integer, intent(in) :: nat
    real(gp), intent(in) :: ratin(3,nat)
    real(gp), intent(inout) :: vat(3,nat)
    !internal
    integer :: iat
    real(gp) :: rat(3,nat),cmx,cmy,cmz
    real(gp) :: tt,tinv(3,3),dettinI,w(3),vo(3),angmom(3)
    real(gp) :: tin11,tin22,tin33,tin12,tin13
    real(gp) :: tin23,tin21,tin31,tin32
    real(gp), dimension(:), allocatable :: amass

    amass = f_malloc((/1.to.nat/),id='amass')
    amass(1:nat)=1.0_gp

    !shift cm to origin
    cmx=0.0_gp; cmy=0.0_gp; cmz=0.0_gp
    do iat=1,nat
        cmx = cmx + ratin(1,iat)
        cmy = cmy + ratin(2,iat)
        cmz = cmz + ratin(3,iat)
    enddo
    cmx=cmx/nat; cmy=cmy/nat; cmz=cmz/nat
    do iat=1,nat
        rat(1,iat) = ratin(1,iat)-cmx
        rat(2,iat) = ratin(2,iat)-cmy
        rat(3,iat) = ratin(3,iat)-cmz
    enddo
    
    
    !calculate inertia tensor
    tin11=0.0_gp; tin22=0.0_gp; tin33=0.0_gp
    tin12=0.0_gp; tin13=0.0_gp; tin23=0.0_gp
    tin21=0.0_gp; tin31=0.0_gp; tin32=0.0_gp
    do iat=1,nat
        tt=amass(iat)
        tin11=tin11+tt*(rat(2,iat)*rat(2,iat)+rat(3,iat)*rat(3,iat))
        tin22=tin22+tt*(rat(1,iat)*rat(1,iat)+rat(3,iat)*rat(3,iat))
        tin33=tin33+tt*(rat(1,iat)*rat(1,iat)+rat(2,iat)*rat(2,iat))
        tin12=tin12-tt*(rat(1,iat)*rat(2,iat))
        tin13=tin13-tt*(rat(1,iat)*rat(3,iat))
        tin23=tin23-tt*(rat(2,iat)*rat(3,iat))
        tin21=tin12
        tin31=tin13
        tin32=tin23
    enddo

    !invert inertia tensor
    dettinI = 1.0_gp/(tin11*tin22*tin33+tin12*tin23*tin31&
              +tin13*tin21*tin32-tin13*tin22*tin31-&
              tin12*tin21*tin33-tin11*tin23*tin32)
    tinv(1,1)=dettinI*(tin22*tin33-tin23*tin32)
    tinv(2,2)=dettinI*(tin11*tin33-tin13*tin31)
    tinv(3,3)=dettinI*(tin11*tin22-tin12*tin21)
    tinv(1,2)=dettinI*(tin13*tin32-tin12*tin33)
    tinv(1,3)=dettinI*(tin12*tin23-tin13*tin22)
    tinv(2,3)=dettinI*(tin13*tin21-tin11*tin23)
    tinv(2,1)=dettinI*(tin23*tin31-tin21*tin33)
    tinv(3,1)=dettinI*(tin21*tin32-tin22*tin31)
    tinv(3,2)=dettinI*(tin12*tin31-tin11*tin32)

    !Compute angular momentum
    angmom=0.0_gp
    vo=0.0_gp
    do iat=1,nat
        vo(1)=rat(2,iat)*vat(3,iat)-rat(3,iat)*vat(2,iat)
        vo(2)=rat(3,iat)*vat(1,iat)-rat(1,iat)*vat(3,iat)
        vo(3)=rat(1,iat)*vat(2,iat)-rat(2,iat)*vat(1,iat)
        angmom(1)=angmom(1)+vo(1)
        angmom(2)=angmom(2)+vo(2)
        angmom(3)=angmom(3)+vo(3)
    enddo

    !matrix product w= I^-1 L
    w(1)=tinv(1,1)*angmom(1)+tinv(1,2)*angmom(2)+tinv(1,3)*angmom(3)
    w(2)=tinv(2,1)*angmom(1)+tinv(2,2)*angmom(2)+tinv(2,3)*angmom(3)
    w(3)=tinv(3,1)*angmom(1)+tinv(3,2)*angmom(2)+tinv(3,3)*angmom(3)

    vo=0.0_gp
    do iat=1,nat
        !tangential velocity v_ortho = w x r_i
        vo(1)=w(2)*rat(3,iat)-w(3)*rat(2,iat)
        vo(2)=w(3)*rat(1,iat)-w(1)*rat(3,iat)
        vo(3)=w(1)*rat(2,iat)-w(2)*rat(1,iat)
        !remove tangential velocity from velocity of atom iat
        vat(1,iat) = vat(1,iat) - vo(1)
        vat(2,iat) = vat(2,iat) - vo(2)
        vat(3,iat) = vat(3,iat) - vo(3)
    enddo
    call f_free(amass)
end subroutine
!=====================================================================
subroutine elim_moment_fs(nat,vxyz)
    use module_base
    implicit none
    !parameters
    integer :: nat
    !internal
    integer :: iat
    real(gp) :: vxyz(3,nat),sx,sy,sz
    sx=0.0_gp ; sy=0.0_gp ; sz=0.0_gp
    do iat=1,nat
        sx=sx+vxyz(1,iat)
        sy=sy+vxyz(2,iat)
        sz=sz+vxyz(3,iat)
    enddo
    sx=sx/nat ; sy=sy/nat ; sz=sz/nat
    do iat=1,nat
        vxyz(1,iat)=vxyz(1,iat)-sx
        vxyz(2,iat)=vxyz(2,iat)-sy
        vxyz(3,iat)=vxyz(3,iat)-sz
    enddo
end subroutine
!=====================================================================
subroutine clean_minmode(nat,geocode,rxyz,minmode)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    character(len=*), intent(in) :: geocode
    real(gp), intent(in) :: rxyz(3,nat)
    real(gp), intent(inout) :: minmode(3,nat)
    !local
    !functions
    real(gp) :: dnrm2
    call elim_moment_fs(nat,minmode(1,1))
    if (trim(adjustl(geocode))=='F')call elim_torque_bastian(nat,&
                                          rxyz(1,1),minmode(1,1))
    minmode = minmode / dnrm2(3*nat,minmode(1,1),1)
end subroutine clean_minmode
end module

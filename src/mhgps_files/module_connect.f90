!> @file
!! module implementing the connection algorithm(s)
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_connect
    use module_base
    use module_interfaces
    use module_io
    implicit none
    
    private

    public :: connect_object
    public :: allocate_connect_object
    public :: deallocate_connect_object
    public :: connect_recursively

    type connect_object
        real(gp), allocatable :: saddle(:,:,:)
        real(gp), allocatable :: enersad(:)
        real(gp), allocatable :: fsad(:,:,:)
        real(gp), allocatable :: fpsad(:,:)
        real(gp), allocatable :: rotforce(:,:,:)
        real(gp), allocatable :: minmode(:,:,:)

        real(gp), allocatable :: leftmin(:,:,:)
        real(gp), allocatable :: enerleft(:)
        real(gp), allocatable :: fleft(:,:,:)
        real(gp), allocatable :: fpleft(:,:)

        real(gp), allocatable :: rightmin(:,:,:)
        real(gp), allocatable :: enerright(:)
        real(gp), allocatable :: fright(:,:,:)
        real(gp), allocatable :: fpright(:,:)

        real(gp), allocatable :: rxyz1(:,:)
        real(gp), allocatable :: rxyz2(:,:)
    end type

contains
!=====================================================================
subroutine allocate_connect_object(nat,nid,nsadmax,cobj)
    use dynamic_memory
    implicit none
    !parameters
    integer, intent(in) :: nat, nid, nsadmax
    type(connect_object), intent(inout) :: cobj

    cobj%saddle    = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='saddle')
    cobj%enersad   = f_malloc((/1.to.nsadmax/),id='enersad')
    cobj%fsad      = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fsad')
    cobj%fpsad     = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpsad')
    cobj%rotforce  = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='rotorce')
    cobj%minmode   = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='minmode')
    cobj%leftmin   = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='leftmin')
    cobj%enerleft  = f_malloc((/1.to.nsadmax/),id='enerleft')
    cobj%fleft     = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fleft')
    cobj%fpleft    = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpleft')
    cobj%rightmin  = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='rightmin')
    cobj%enerright = f_malloc((/1.to.nsadmax/),id='enerright')
    cobj%fright    = f_malloc((/1.to.3,1.to.nat,1.to.nsadmax/),&
                             id='fright')
    cobj%fpright   = f_malloc((/1.to.nid,1.to.nsadmax/),id='fpright')
    cobj%rxyz1 = f_malloc((/1.to.3,1.to.nat/),id='rxyz1')
    cobj%rxyz2 = f_malloc((/1.to.3,1.to.nat/),id='rxyz2')
end subroutine
!=====================================================================
subroutine deallocate_connect_object(cobj)
    use dynamic_memory
    implicit none
    !parameters
    type(connect_object), intent(inout) :: cobj

    call f_free(cobj%saddle)
    call f_free(cobj%enersad)
    call f_free(cobj%fsad)
    call f_free(cobj%fpsad)
    call f_free(cobj%rotforce)
    call f_free(cobj%minmode)
    call f_free(cobj%leftmin)
    call f_free(cobj%enerleft)
    call f_free(cobj%fleft)
    call f_free(cobj%fpleft)
    call f_free(cobj%rightmin)
    call f_free(cobj%enerright)
    call f_free(cobj%fright)
    call f_free(cobj%fpright)
    call f_free(cobj%rxyz1)
    call f_free(cobj%rxyz2)
end subroutine
!=====================================================================
recursive subroutine connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz1,rxyz2,ener1,ener2,fp1,fp2,&
                     nsad,cobj,connected)
    !if called from outside recursion, connected has to be set 
    !to .true. and nsad=0
    use module_base
    use module_global_variables,&
       only: atoms,&
             imode,&
             nsadmax,&
             iproc,&
             isad,&
             isadc,&
             currDir,&
             ixyz_int,&
             en_delta_min, fp_delta_min,&
             en_delta_sad, fp_delta_sad
    use module_ls_rmsd
    use module_fingerprints
    use module_minimizers
    use yaml_output
    use module_saddle
    use module_freezingstring
use module_energyandforces
    implicit none
    !parameters
    integer, intent(in)     :: nat
    integer, intent(in)     :: nid
    integer, intent(in)     :: nbond
    real(gp), intent(in)    :: rcov(nat)
    real(gp), intent(inout) :: alat(3)
    integer, intent(in)     :: iconnect(2,nbond)
    real(gp), intent(in)    :: rxyz1(3,nat), rxyz2(3,nat)
    real(gp), intent(in)    :: fp1(nid), fp2(nid)
    real(gp), intent(in)    :: ener1,ener2
    integer, intent(inout)  :: nsad
    type(connect_object), intent(inout) :: cobj
    logical, intent(inout)    :: connected
    !internal
    integer  :: nsad_loc
    real(gp) :: displ=0._gp,ener_count=0._gp
    real(gp) :: fnoise,fnrm,fmax
    logical  :: converged =.false.
    logical  :: lnl, rnr, lnr, rnl 
    character(len=200) :: comment
    real(gp) :: tsgforces(3,nat), tsgenergy

    if(.not.connected)return
    if(nsad>=nsadmax.and.iproc==0)then
        connected=.false.
        return
    endif
    if(iproc==0)then
        call yaml_comment('(MHGPS) nsad:'//&
             trim(adjustl(yaml_toa(nsad)))//'; connect minima with &
             following energies')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener1)))//' and ')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener2))))
    endif


    !check if input structures are distinct 
    if(equal(nid,en_delta_min,fp_delta_min,ener1,ener2,fp1,fp2)&
                                                    .and.iproc==0)then
        call yaml_warning('(MHGPS)  connect: input minima are&
                           identical. Will NOT attempt to find&
                           an intermediate TS. recursion depth: '&
                           //yaml_toa(nsad))
        return
    endif

    call vcopy(3*nat,rxyz1(1,1),1,cobj%rxyz1(1,1), 1)
    call vcopy(3*nat,rxyz2(1,1),1,cobj%rxyz2(1,1), 1)
    !rmsd alignment (optional in mhgps approach)
    call superimpose(nat,cobj%rxyz1,cobj%rxyz2)

    !get input guess for transition state
    nsad=nsad+1
    nsad_loc=nsad
    isad=isad+1
    write(isadc,'(i5.5)')isad

    call get_ts_guess(nat,alat,cobj%rxyz1,cobj%rxyz2,&
          cobj%saddle(1,1,nsad),cobj%minmode(1,1,nsad),tsgenergy,&
          tsgforces(1,1))


    !compute saddle
    ener_count=0.0_gp
    call findsad(nat,alat,rcov,nbond,iconnect,cobj%saddle(1,1,nsad),&
                cobj%enersad(nsad),cobj%fsad(1,1,nsad),&
                cobj%minmode(1,1,nsad),displ,ener_count,&
                cobj%rotforce(1,1,nsad),converged)

    if(.not.converged)then
        nsad=nsad-1!in case we don't want to STOP
        isad=isad-1
        write(isadc,'(i5.5)')isad
        stop 'STOP saddle not converged'
    endif

    call fnrmandforcemax(cobj%fsad(1,1,nsad),fnrm,fmax,nat)
    if (iproc == 0) then
        write(comment,'(a,1pe10.3,5x1pe10.3)')'ATTENTION! Forces &
        below give no forces, but the final minmode| fnrm, fmax = ',&
        fnrm,fmax

        call write_atomic_file(currDir//'/sad'//trim(adjustl(isadc))&
        //'_finalM',cobj%enersad(nsad),cobj%saddle(1,1,nsad),&
        ixyz_int,atoms,comment,forces=cobj%minmode(1,1,nsad))

        write(comment,'(a,1pe10.3,5x1pe10.3)')&
                                            'fnrm, fmax = ',fnrm,fmax
        call write_atomic_file(currDir//'/sad'//trim(adjustl(isadc))&
        //'_finalF',cobj%enersad(nsad),cobj%saddle(1,1,nsad),&
        ixyz_int,atoms,comment,forces=cobj%fsad(1,1,nsad))

        call write_mode(nat,currDir//'/sad'//trim(adjustl(isadc))//&
        '_mode_final',cobj%minmode(1,1,nsad),cobj%rotforce(1,1,nsad))
    endif


    call fingerprint(nat,nid,alat,atoms%astruct%geocode,rcov,&
                    cobj%saddle(1,1,nsad),cobj%fpsad(1,nsad))

    !pushoff and minimize left and right
    call pushoff(nat,cobj%saddle(1,1,nsad),cobj%minmode(1,1,nsad),&
                 cobj%leftmin(1,1,nsad),cobj%rightmin(1,1,nsad))

    if(iproc==0)&
    call yaml_comment('(MHGPS) Relax from left side ',hfill='.')
    ener_count=0.0_gp
    call energyandforces(nat,alat,cobj%leftmin(1,1,nsad),&
    cobj%fleft(1,1,nsad),fnoise,cobj%enerleft(nsad))
    call minimize(imode,nat,alat,nbond,iconnect,&
                        cobj%leftmin(1,1,nsad),cobj%fleft(1,1,nsad),&
                        fnoise,cobj%enerleft(nsad),&
                        ener_count,converged,'L')
    call fnrmandforcemax(cobj%fleft(1,1,nsad),fnrm,fmax,nat)
    write(comment,'(a,1pe10.3,5x1pe10.3)')'fnrm, fmax = ',fnrm,fmax
    if(iproc==0)&
    call write_atomic_file(currDir//'/sad'//trim(adjustl(isadc))&
    //'_minFinalL',cobj%enerleft(nsad),cobj%leftmin(1,1,nsad),&
    ixyz_int,atoms,comment,cobj%fleft(1,1,nsad))

    if(iproc==0)&
    call yaml_comment('(MHGPS) Relax from right side ',hfill='.')
    ener_count=0.0_gp
    call energyandforces(nat,alat,cobj%rightmin(1,1,nsad),&
    cobj%fright(1,1,nsad),fnoise,cobj%enerright(nsad))
    call minimize(imode,nat,alat,nbond,iconnect,&
                        cobj%rightmin(1,1,nsad),cobj%fright(1,1,nsad)&
                       ,fnoise,cobj%enerright(nsad),&
                        ener_count,converged,'R')
    call fnrmandforcemax(cobj%fright(1,1,nsad),fnrm,fmax,nat)
    write(comment,'(a,1pe10.3,5x1pe10.3)')'fnrm, fmax = ',fnrm,fmax
    if(iproc==0)&
    call write_atomic_file(currDir//'/sad'//trim(adjustl(isadc))&
    //'_minFinalR',cobj%enerright(nsad),cobj%rightmin(1,1,nsad),&
    ixyz_int,atoms,comment,cobj%fright(1,1,nsad))

    call fingerprint(nat,nid,alat,atoms%astruct%geocode,rcov,&
                    cobj%leftmin(1,1,nsad),cobj%fpleft(1,nsad))
    call fingerprint(nat,nid,alat,atoms%astruct%geocode,rcov,&
                    cobj%rightmin(1,1,nsad),cobj%fpright(1,nsad))
    !check if relaxed structures are identical to saddle itself
    if(equal(nid,en_delta_sad,fp_delta_sad,cobj%enersad(nsad),&
    cobj%enerright(nsad),cobj%fpsad(1,nsad),cobj%fpright(1,nsad)).or.&
    equal(nid,en_delta_sad,fp_delta_sad,cobj%enersad(nsad),&
    cobj%enerleft(nsad),cobj%fpsad(1,nsad),cobj%fpleft(1,nsad)))then

        connected=.false.
        nsad=nsad-1
        isad=isad-1
        write(isadc,'(i5.5)')isad

        if(iproc==0)&
        call yaml_warning('(MHGPS)  after relaxation from saddle &
                           point the left and/or right minimum are &
                           identical to the saddle point. Stopped &
                           connection attempt. Will proceed with &
                           next connection attempt.')
        return
    endif

    !is minimum, obtained by relaxation from left bar end identical to
    !left input minimum?
    lnl=equal(nid,en_delta_min,fp_delta_min,ener1,&
        cobj%enerleft(nsad),fp1,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !right input minimum?
    rnr=equal(nid,en_delta_min,fp_delta_min,ener2,&
        cobj%enerright(nsad),fp2,cobj%fpright(1,nsad))

    !is minimum obtained by relaxation from left bar end identical to 
    !right input minimum?
    lnr=equal(nid,en_delta_min,fp_delta_min,ener2,&
        cobj%enerleft(nsad),fp2,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !left input minimum?
    rnl=equal(nid,en_delta_min,fp_delta_min,ener1,&
        cobj%enerright(nsad),fp1,cobj%fpright(1,nsad))

    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
        connected=.true.
        return
    endif

    if(lnl .and. (.not. rnr))then
        !connect right input min with right relaxed bar-end
if(iproc==0)write(*,*)'(MHGS) connection check lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
        return
    endif

    if(rnr .and. (.not. lnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',rnr,lnl
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
!write(*,*)rxyz1
!write(*,*)
!write(*,*)cobj%leftmin(:,:,nsad_loc)
!if(sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))<1.d-2)then
!    lnl=equal(nid,en_delta_min,fp_delta_min,ener1,&
!        cobj%enerleft(nsad),fp1,cobj%fpleft(1,nsad))
!!call fingerprint(nat,nid,alat,atoms%astruct%geocode,rcov,&
!!                cobj%leftmin(1,1,nsad_loc),cobj%fpleft(1,nsad_loc))
!!call fingerprint(nat,nid,alat,atoms%astruct%geocode,rcov,&
!!                rxyz1(1,1),cobj%fpright(1,nsad_loc))
!write(*,*)'fprints:'
!write(*,*)cobj%fpleft(:,nsad_loc)
!write(*,*)
!write(*,*)fp1
!
!    stop
!endif
        !connect left relaxed bar end with left input min
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if(lnr .and. (.not. rnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check lnr and not rnl',sqrt(sum((rxyz1-cobj%rightmin(:,:,nsad_loc))**2))
        !connect right relaxed bar end with left input min
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz1,cobj%rightmin(1,1,nsad_loc),&
                     ener1,cobj%enerright(nsad_loc),&
                     fp1,cobj%fpright(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if(.not. lnr .and. rnl)then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnr and rnl',sqrt(sum((rxyz2-cobj%leftmin(:,:,nsad_loc))**2))
        !connect left relaxed bar end with right input min
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz2,cobj%leftmin(1,1,nsad_loc),&
                     ener2,cobj%enerleft(nsad_loc),&
                     fp2,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if((.not. lnl) .and. (.not. rnr))then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
        !connect left input min with left relaxed bar end  and right
        !input min with right relaxed bar end
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
        return
    endif

    !should and must not happen:
    if(iproc==0)&
    call yaml_warning('(MHGPS) Severe error in connect: none of &
                  the checks in connect subroutine were successful! &
                  STOP') 
    stop '(MHGPS) Severe error in connect: none of &
                  the checks in connect subroutine were successful! &
                  STOP'

end subroutine
!=====================================================================
subroutine pushoff(nat,saddle,minmode,left,right)
    use module_base
    use module_misc
    use module_global_variables, only: saddle_stepoff
    implicit none
    !parameters 
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(out) :: left(3,nat)
    real(gp), intent(out) :: right(3,nat)
    !internal
    real(gp)  :: step(3,nat)

    !functions
    real(gp) :: dnrm2

    step = saddle_stepoff*minmode
    left = saddle - step
    right = saddle + step

end subroutine
!=====================================================================


end module

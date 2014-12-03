!> @file
!! module implementing the connection algorithm(s)
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

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
    public :: connect

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
recursive subroutine connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,rxyz1,rxyz2,ener1,ener2,fp1,fp2,&
                     nsad,cobj,connected)
    !if called from outside recursion, connected has to be set 
    !to .true. and nsad=0
  use module_base
  use module_atoms, only: astruct_dump_to_file
    use module_global_variables,&
       only: astruct_ptr,&
             imode,&
             nsadmax,&
             iproc,&
             isad,isadprob,&
             isadc,isadprobc,&
             currDir,&
             en_delta_min, fp_delta_min,&
             en_delta_sad, fp_delta_sad,&
             saddle_scale_stepoff,&
             mhgps_verbosity
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
    integer, intent(inout)  :: nsad,isame
    type(connect_object), intent(inout) :: cobj
    logical, intent(inout)    :: connected
    !internal
    integer  :: nsad_loc,ipush
    real(gp) :: displ=0._gp,ener_count=0._gp
    real(gp) :: fnoise,fnrm,fmax
    logical  :: converged =.false.
    logical  :: lnl, rnr, lnr, rnl 
    character(len=200) :: comment
    real(gp) :: tsgforces(3,nat), tsgenergy,scl

    if(.not.connected)return
    if(nsad>=nsadmax)then
        connected=.false.
        return
    endif
    if(iproc==0)then
        call yaml_comment('(MHGPS) nsad:'//&
             trim(adjustl(yaml_toa(nsad)))//'; connect minima with'//&
             ' following energies')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener1)))&
                           //' and ')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener2))))
    endif


    !check if input structures are distinct 
    if(equal('MM',nid,en_delta_min,fp_delta_min,ener1,ener2,fp1,fp2))then
        if(iproc==0)call yaml_warning('(MHGPS) connect: input '//&
                    'minima are identical. Will NOT attempt to '//&
                    'find an intermediate TS. recursion depth: '//&
                    yaml_toa(nsad))
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
        connected=.false.
        stop 'STOP saddle not converged'
    endif

    call fnrmandforcemax(cobj%fsad(1,1,nsad),fnrm,fmax,nat)
    fnrm=sqrt(fnrm)
    if (iproc == 0) then
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'ATTENTION! Forces '//&
        'below give no forces, but the final minmode| fnrm, fmax = ',&
        fnrm,fmax

        call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_finalM',&
             comment,&
             cobj%enersad(nsad),rxyz=cobj%saddle(:,:,nsad),&
             forces=cobj%minmode(:,:,nsad))

        write(comment,'(a,1pe10.3,5x,1pe10.3)')&
                                            'fnrm, fmax = ',fnrm,fmax

        call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_finalF',&
             comment,&
             cobj%enersad(nsad),rxyz=cobj%saddle(:,:,nsad),&
             forces=cobj%fsad(:,:,nsad))

        call write_mode(nat,currDir//'/sad'//trim(adjustl(isadc))//&
        '_mode_final',cobj%minmode(1,1,nsad),cobj%rotforce(1,1,nsad))
    endif


    call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
                    cobj%saddle(1,1,nsad),cobj%fpsad(1,nsad))

    if(nsad>1)then
        if(equal('SS',nid,en_delta_sad,fp_delta_sad,&
               cobj%enersad(nsad-1),cobj%enersad(nsad),&
               cobj%fpsad(1,nsad-1),cobj%fpsad(1,nsad)))then
            isame=isame+1
            !if we find the same saddle point there times (cosecutively)
            !we are stuck
            if(isame==3)then
                connected=.false.
                nsad=nsad-1
                isad=isad-1
                write(isadc,'(i5.5)')isad
    
                if(iproc==0)&
                call yaml_warning('(MHGPS) found same saddle '//&
                                   'point again. Aborting connection'//&
                                   ' attempt.')
                return
            endif
        else
            isame=0
        endif
    endif

    scl=-1.0_gp
    ipush=1
    loopL: do
        if(iproc==0)&
        call yaml_comment('(MHGPS) Relax from left side ',hfill='.')
    
        call pushoffsingle(nat,cobj%saddle(1,1,nsad),&
        cobj%minmode(1,1,nsad),scl,cobj%leftmin(1,1,nsad))

        ener_count=0.0_gp
        call energyandforces(nat,alat,cobj%leftmin(1,1,nsad),&
        cobj%fleft(1,1,nsad),fnoise,cobj%enerleft(nsad))

        if(iproc==0 .and. mhgps_verbosity >= 3)&
             call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_pushL',&
             comment,&
             cobj%enerleft(nsad),cobj%leftmin(:,:,nsad),&
             cobj%fleft(:,:,nsad))

        call minimize(imode,nat,alat,nbond,iconnect,&
                            cobj%leftmin(1,1,nsad),&
                            cobj%fleft(1,1,nsad),fnoise,&
                            cobj%enerleft(nsad),ener_count,converged,&
                            'L')
        call fnrmandforcemax(cobj%fleft(1,1,nsad),fnrm,fmax,nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(iproc==0)&
             call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_minFinalL',&
             comment,&
             cobj%enerleft(nsad),cobj%leftmin(:,:,nsad),&
             cobj%fleft(:,:,nsad))

        call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
                        cobj%leftmin(1,1,nsad),cobj%fpleft(1,nsad))
        if(.not.equal('MS',nid,en_delta_sad,fp_delta_sad,&
           cobj%enersad(nsad),cobj%enerleft(nsad),cobj%fpsad(1,nsad),&
           cobj%fpleft(1,nsad)))then
           exit loopL 
        elseif(ipush>=3)then
            isadprob=isadprob+1
            write(isadprobc,'(i5.5)')isadprob
            if(iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                'unknown (stepoff converged back to saddle)'
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_finalM',&
                     comment,&
                cobj%enersad(nsad),rxyz=cobj%saddle(:,:,nsad),&
                forces=cobj%minmode(:,:,nsad))
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_Reactant',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_Product',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz2)
            endif
    
            connected=.false.
            nsad=nsad-1
            isad=isad-1
            write(isadc,'(i5.5)')isad
    
            if(iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the left minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            return
        endif
        scl=saddle_scale_stepoff*scl
        if(iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, left side '//&
                       'converged back to saddle. Will retry with '//&
                       'increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopL

    scl=1.0_gp
    ipush=1
    loopR: do
        if(iproc==0)&
        call yaml_comment('(MHGPS) Relax from right side ',hfill='.')
    
        call pushoffsingle(nat,cobj%saddle(1,1,nsad),&
        cobj%minmode(1,1,nsad),scl,cobj%rightmin(1,1,nsad))

        ener_count=0.0_gp
        call energyandforces(nat,alat,cobj%rightmin(1,1,nsad),&
        cobj%fright(1,1,nsad),fnoise,cobj%enerright(nsad))

        if(iproc==0 .and. mhgps_verbosity >= 3)&
             call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_pushR',&
             comment,&
             cobj%enerright(nsad),cobj%rightmin(1,1,nsad),&
             cobj%fright(1,1,nsad))

        call minimize(imode,nat,alat,nbond,iconnect,&
                            cobj%rightmin(1,1,nsad),&
                            cobj%fright(1,1,nsad),fnoise,&
                            cobj%enerright(nsad),ener_count,&
                            converged,'R')
        call fnrmandforcemax(cobj%fright(1,1,nsad),fnrm,fmax,nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(iproc==0)&
             call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_minFinalR',&
             comment,&
             cobj%enerright(nsad),cobj%rightmin(1,1,nsad),&
             cobj%fright(1,1,nsad))
        call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
                        cobj%rightmin(1,1,nsad),cobj%fpright(1,nsad))
        if(.not.equal('MS',nid,en_delta_sad,fp_delta_sad,&
           cobj%enersad(nsad),cobj%enerright(nsad),&
           cobj%fpsad(1,nsad),cobj%fpright(1,nsad)))then
           exit loopR 
        elseif(ipush>=3)then
            isadprob=isadprob+1
            write(isadprobc,'(i5.5)')isadprob
            if(iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                     'unknown (stepoff converged back to saddle)'
                        
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_finalM',&
                     comment,&
                     cobj%enersad(nsad),cobj%saddle(:,:,nsad),&
                     forces=cobj%minmode(:,:,nsad))
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_Reactant',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(astruct_ptr,&
                     currDir//'/sadProb'//trim(adjustl(isadprobc))//'_Product',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz2)

            endif
    
            connected=.false.
            nsad=nsad-1
            isad=isad-1
            write(isadc,'(i5.5)')isad
    
            if(iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the right minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            return
        endif
        scl=saddle_scale_stepoff*scl
        if(iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, right side'//&
                       ' converged back to saddle. Will retry with'//&
                       ' increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopR

    !is minimum, obtained by relaxation from left bar end identical to
    !left input minimum?
    lnl=equal('MM',nid,en_delta_min,fp_delta_min,ener1,&
        cobj%enerleft(nsad),fp1,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !right input minimum?
    rnr=equal('MM',nid,en_delta_min,fp_delta_min,ener2,&
        cobj%enerright(nsad),fp2,cobj%fpright(1,nsad))

    !is minimum obtained by relaxation from left bar end identical to 
    !right input minimum?
    lnr=equal('MM',nid,en_delta_min,fp_delta_min,ener2,&
        cobj%enerleft(nsad),fp2,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !left input minimum?
    rnl=equal('MM',nid,en_delta_min,fp_delta_min,ener1,&
        cobj%enerright(nsad),fp1,cobj%fpright(1,nsad))

    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        connected=.true.
        return
    endif

    if(lnl .and. (.not. rnr))then
        !connect right input min with right relaxed bar-end
if(iproc==0)write(*,*)'(MHGPS) connection check lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
        return
    endif

    if(rnr .and. (.not. lnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',rnr,lnl
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
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
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if(lnr .and. (.not. rnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check lnr and not rnl',sqrt(sum((rxyz1-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect right relaxed bar end with left input min
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%rightmin(1,1,nsad_loc),&
                     ener1,cobj%enerright(nsad_loc),&
                     fp1,cobj%fpright(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if(.not. lnr .and. rnl)then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnr and rnl',sqrt(sum((rxyz2-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect left relaxed bar end with right input min
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,rxyz2,cobj%leftmin(1,1,nsad_loc),&
                     ener2,cobj%enerleft(nsad_loc),&
                     fp2,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        return
    endif

    if((.not. lnl) .and. (.not. rnr))then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect left input min with left relaxed bar end  and right
        !input min with right relaxed bar end
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
        call connect_recursively(nat,nid,alat,rcov,nbond,isame,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
        return
    endif

    !should and must not happen:
    if(iproc==0)&
    call yaml_warning('(MHGPS) Severe error in connect: none of '//&
                      'the checks in connect subroutine were '//&
                       'successful! STOP') 
    call f_err_throw('(MHGPS) Severe error in connect: none of '//&
         'the checks in connect subroutine were '//&
         'successful! STOP')

end subroutine
!=====================================================================
!same as connect_recursively, but in an iterative fashion
subroutine connect(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz1,rxyz2,ener1,ener2,fp1,fp2,&
                     nsad,cobj,connected)
    !if called from outside recursion, connected has to be set 
    !to .true. and nsad=0
    use module_base
    use module_atoms, only: astruct_dump_to_file
    use module_global_variables,&
       only: astruct_ptr,&
             imode,&
             nsadmax,&
             iproc,&
             isad,isadprob,&
             isadc,isadprobc,&
             currDir,&
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
    real(gp) :: todorxyz(3,nat,2,2*nsadmax),todofp(nid,2,2*nsadmax)
    real(gp) :: todoenergy(2,2*nsadmax)
    integer :: ntodo
    real(gp) :: rxyz1cur(3,nat), rxyz2cur(3,nat)
    real(gp) :: fp1cur(nid), fp2cur(nid)
    real(gp) :: ener1cur, ener2cur
integer :: i

    ntodo=1
    todorxyz(:,:,1,ntodo)=rxyz1
    todorxyz(:,:,2,ntodo)=rxyz2
    todofp(:,1,ntodo)=fp1
    todofp(:,2,ntodo)=fp2
    todoenergy(1,ntodo)=ener1
    todoenergy(2,ntodo)=ener2

connectloop: do while(ntodo>=1)
    if(nsad>=nsadmax)then
        connected=.false.
        exit connectloop
    endif

    call vcopy(3*nat,todorxyz(1,1,1,ntodo),1,rxyz1cur(1,1), 1)
    call vcopy(nid,todofp(1,1,ntodo),1,fp1cur(1), 1)
    ener1cur=todoenergy(1,ntodo)
    call vcopy(3*nat,todorxyz(1,1,2,ntodo),1,rxyz2cur(1,1), 1)
    call vcopy(nid,todofp(1,2,ntodo),1,fp2cur(1), 1)
    ener2cur=todoenergy(2,ntodo)
    if(iproc==0)then
        call yaml_comment('(MHGPS) nsad:'//&
             trim(adjustl(yaml_toa(nsad)))//'; connect minima '//&
             'with following energies')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(todoenergy(1,ntodo))))//' and ')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(todoenergy(2,ntodo)))))
    endif


    !check if input structures are distinct 
    if(equal('MM',nid,en_delta_min,fp_delta_min,todoenergy(1,ntodo),todoenergy(2,ntodo),todofp(1,1,ntodo),todofp(1,2,ntodo)))then
        if(iproc==0)call yaml_warning('(MHGPS) connect: input '//&
                    'minima are identical. Will NOT attempt to find '//&
                    'an intermediate TS. recursion depth: '//&
                    yaml_toa(nsad))
        exit connectloop
    endif
    !rmsd alignment (optional in mhgps approach)
    call superimpose(nat,rxyz1cur,rxyz2cur)

    !get input guess for transition state
    nsad=nsad+1
    nsad_loc=nsad
    isad=isad+1
    write(isadc,'(i5.5)')isad

    call get_ts_guess(nat,alat,rxyz1cur,rxyz2cur,&
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
        connected=.false.
        stop 'STOP saddle not converged'
        !exit connectloop
    endif
    ntodo=ntodo-1

    call fnrmandforcemax(cobj%fsad(1,1,nsad),fnrm,fmax,nat)
    fnrm=sqrt(fnrm)
    if (iproc == 0) then
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'ATTENTION! Forces '//&
        'below give no forces, but the final minmode| fnrm, fmax = ',&
        fnrm,fmax

        call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_finalM',&
             comment,&
             cobj%enersad(nsad),cobj%saddle(:,:,nsad),&
             forces=cobj%minmode(:,:,nsad))

        write(comment,'(a,1pe10.3,5x,1pe10.3)')&
                                            'fnrm, fmax = ',fnrm,fmax
        call astruct_dump_to_file(astruct_ptr,&
             currDir//'/sad'//trim(adjustl(isadc))//'_finalF',&
             comment,&
             cobj%enersad(nsad),cobj%saddle(:,:,nsad),&
             forces=cobj%minmode(:,:,nsad))

        call write_mode(nat,currDir//'/sad'//trim(adjustl(isadc))//&
        '_mode_final',cobj%minmode(1,1,nsad),cobj%rotforce(1,1,nsad))
    endif


    call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
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
    fnrm=sqrt(fnrm)
    write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,fmax
    if(iproc==0)&
         call astruct_dump_to_file(astruct_ptr,&
         currDir//'/sad'//trim(adjustl(isadc))//'_minFinalL',&
         comment,&
         cobj%enerleft(nsad),cobj%leftmin(:,:,nsad),&
         forces=cobj%fleft(:,:,nsad))

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
    fnrm=sqrt(fnrm)
    write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,fmax
    if(iproc==0)&
         call astruct_dump_to_file(astruct_ptr,&
         currDir//'/sad'//trim(adjustl(isadc))//'_minFinalR',&
         comment,&
         cobj%enerright(nsad),cobj%rightmin(:,:,nsad),&
         cobj%fright(:,:,nsad))

    call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
                    cobj%leftmin(1,1,nsad),cobj%fpleft(1,nsad))
    call fingerprint(nat,nid,alat,astruct_ptr%geocode,rcov,&
                    cobj%rightmin(1,1,nsad),cobj%fpright(1,nsad))
    !check if relaxed structures are identical to saddle itself
    if(equal('MS',nid,en_delta_sad,fp_delta_sad,cobj%enersad(nsad),&
    cobj%enerright(nsad),cobj%fpsad(1,nsad),cobj%fpright(1,nsad)).or.&
    equal('MS',nid,en_delta_sad,fp_delta_sad,cobj%enersad(nsad),&
    cobj%enerleft(nsad),cobj%fpsad(1,nsad),cobj%fpleft(1,nsad)))then

        isadprob=isadprob+1
        write(isadprobc,'(i5.5)')isadprob
        if(iproc==0)then
            write(comment,'(a)')'Prob: Neighbors '//&
            'unknown (converged back to saddle after stepoff)'
            call astruct_dump_to_file(astruct_ptr,&
                 currDir//'/sadProb'//trim(adjustl(isadprobc))//'_finalM',&
                 comment,&
                 cobj%enersad(nsad),cobj%saddle(:,:,nsad),&
                 forces=cobj%minmode(:,:,nsad))
        endif

        connected=.false.
        nsad=nsad-1
        isad=isad-1
        write(isadc,'(i5.5)')isad

        if(iproc==0)&
        call yaml_warning('(MHGPS)  after relaxation from saddle '//&
                         'point the left and/or right minimum are '//&
                         'identical to the saddle point. Stopped '//&
                         'connection attempt. Will proceed with '//&
                         'next connection attempt.')
        exit connectloop !stop connection
    endif

    !is minimum, obtained by relaxation from left bar end identical to
    !left input minimum?
    lnl=equal('MM',nid,en_delta_min,fp_delta_min,ener1cur,&
        cobj%enerleft(nsad),fp1cur,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !right input minimum?
    rnr=equal('MM',nid,en_delta_min,fp_delta_min,ener2cur,&
        cobj%enerright(nsad),fp2cur,cobj%fpright(1,nsad))

    !is minimum obtained by relaxation from left bar end identical to 
    !right input minimum?
    lnr=equal('MM',nid,en_delta_min,fp_delta_min,ener2cur,&
        cobj%enerleft(nsad),fp2cur,cobj%fpleft(1,nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !left input minimum?
    rnl=equal('MM',nid,en_delta_min,fp_delta_min,ener1cur,&
        cobj%enerright(nsad),fp1cur,cobj%fpright(1,nsad))

    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        connected=.true.
!        return
    elseif(lnl .and. (.not. rnr))then
        !connect right input min with right relaxed bar-end
if(iproc==0)write(*,*)'(MHGPS) connection check lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=cobj%rightmin(:,:,nsad_loc)
        todorxyz(:,:,2,ntodo)=rxyz2cur
        todofp(:,1,ntodo)=cobj%fpright(:,nsad_loc)
        todofp(:,2,ntodo)=fp2cur
        todoenergy(1,ntodo)=cobj%enerright(nsad_loc)
        todoenergy(2,ntodo)=ener2cur
        
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
!                     cobj%enerright(nsad_loc),ener2,&
!                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
!        return
    elseif(rnr .and. (.not. lnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',rnr,lnl
if(iproc==0)write(*,*)'(MHGPS)connection check rnr and not lnl',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
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
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=rxyz1cur
        todorxyz(:,:,2,ntodo)=cobj%leftmin(:,:,nsad_loc)
        todofp(:,1,ntodo)=fp1cur
        todofp(:,2,ntodo)=cobj%fpleft(:,nsad_loc)
        todoenergy(1,ntodo)=ener1cur
        todoenergy(2,ntodo)=cobj%enerleft(nsad_loc)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
!                     ener1,cobj%enerleft(nsad_loc),&
!                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        return
    elseif(lnr .and. (.not. rnl))then
if(iproc==0)write(*,*)'(MHGPS)connection check lnr and not rnl',sqrt(sum((rxyz1-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect right relaxed bar end with left input min
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=rxyz1cur
        todorxyz(:,:,2,ntodo)=cobj%rightmin(:,:,nsad_loc)
        todofp(:,1,ntodo)=fp1cur
        todofp(:,2,ntodo)=cobj%fpright(:,nsad_loc)
        todoenergy(1,ntodo)=ener1cur
        todoenergy(2,ntodo)=cobj%enerright(nsad_loc)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%rightmin(1,1,nsad_loc),&
!                     ener1,cobj%enerright(nsad_loc),&
!                     fp1,cobj%fpright(1,nsad_loc),nsad,cobj,connected)
!        return
    elseif(.not. lnr .and. rnl)then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnr and rnl',sqrt(sum((rxyz2-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect left relaxed bar end with right input min
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=rxyz2cur
        todorxyz(:,:,2,ntodo)=cobj%leftmin(:,:,nsad_loc)
        todofp(:,1,ntodo)=fp2cur
        todofp(:,2,ntodo)=cobj%fpleft(:,nsad_loc)
        todoenergy(1,ntodo)=ener2cur
        todoenergy(2,ntodo)=cobj%enerleft(nsad_loc)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz2,cobj%leftmin(1,1,nsad_loc),&
!                     ener2,cobj%enerleft(nsad_loc),&
!                     fp2,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        return
    elseif((.not. lnl) .and. (.not. rnr))then
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,*)'(MHGPS)connection check not lnl and not rnr',sqrt(sum((rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
if(iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) connection check connected',cobj%enerleft(nsad),cobj%enerright(nsad)
        !connect left input min with left relaxed bar end  and right
        !input min with right relaxed bar end
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=rxyz1cur
        todorxyz(:,:,2,ntodo)=cobj%leftmin(:,:,nsad_loc)
        todofp(:,1,ntodo)=fp1cur
        todofp(:,2,ntodo)=cobj%fpleft(:,nsad_loc)
        todoenergy(1,ntodo)=ener1cur
        todoenergy(2,ntodo)=cobj%enerleft(nsad_loc)
        ntodo=ntodo+1
if(ntodo>nsadmax)stop 'error: ntodo>nsadmax'
        todorxyz(:,:,1,ntodo)=cobj%rightmin(:,:,nsad_loc)
        todorxyz(:,:,2,ntodo)=rxyz2cur
        todofp(:,1,ntodo)=cobj%fpright(:,nsad_loc)
        todofp(:,2,ntodo)=fp2cur
        todoenergy(1,ntodo)=cobj%enerright(nsad_loc)
        todoenergy(2,ntodo)=ener2cur
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
!                     ener1,cobj%enerleft(nsad_loc),&
!                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
!                     cobj%enerright(nsad_loc),ener2,&
!                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
!        return
    else
        !should and must not happen:
        if(iproc==0)&
        call yaml_warning('(MHGPS) Severe error in connect: none '//&
                          'of the checks in connect subroutine '//&
                          'were successful! STOP') 
        call f_err_throw('(MHGPS) Severe error in connect: none '//&
                          'of the checks in connect subroutine '//&
                          'were successful! STOP')
    endif
do i=1,2*nsadmax
write(*,*)'ener',todoenergy(1,i),todoenergy(2,i)
enddo
enddo connectloop

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
subroutine pushoffsingle(nat,saddle,minmode,scl,pushed)
    use module_base, only: gp
    use module_misc
    use module_global_variables, only: saddle_stepoff
    implicit none
    !parameters 
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in) :: scl
    real(gp), intent(out) :: pushed(3,nat)
    !internal
    integer  :: iat 
    real(gp) :: step(3,nat)
    real(gp) :: maxd, tt, dp
    !functions
    real(gp) :: dnrm2

    tt=0.0_gp
    dp=0.0_gp
    maxd=-huge(1.0_gp)
    do iat=1,nat
        dp=minmode(1,iat)**2+minmode(2,iat)**2+minmode(3,iat)**2
        tt=tt+dp
        maxd=max(maxd,dp)
    enddo
    tt=sqrt(tt)
    maxd=sqrt(maxd)

    step = minmode*saddle_stepoff/maxd
    pushed = saddle + scl*step

    !control:
    tt=0.0_gp
    dp=0.0_gp
    maxd=-huge(1.0_gp)
    do iat=1,nat
        dp=step(1,iat)**2+step(2,iat)**2+step(3,iat)**2
        tt=tt+dp
        maxd=max(maxd,dp)
    enddo
    tt=sqrt(tt)
    maxd=sqrt(maxd)
    write(133,*)scl, maxd

! old:
!    step = saddle_stepoff*minmode
!    pushed = saddle + scl*step

end subroutine
!=====================================================================
subroutine pushoff_assym(nat,saddle,minmode,scll,sclr,left,right)
    use module_base
    use module_misc
    use module_global_variables, only: saddle_stepoff
    implicit none
    !parameters 
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in) :: scll,sclr
    real(gp), intent(out) :: left(3,nat)
    real(gp), intent(out) :: right(3,nat)
    !internal
    real(gp)  :: step(3,nat)

    !functions
    real(gp) :: dnrm2

    step = saddle_stepoff*minmode
    left = saddle - scll*step
    right = saddle + sclr*step
end subroutine
!=====================================================================


end module

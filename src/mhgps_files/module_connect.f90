!> @file
!! Module implementing the connection algorithm(s)
!!
!! @author 
!!    Copyright (C) 2014 UNIBAS, Bastian Schaefer 
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module implementing the connection algorithm(s)
module module_connect
    use module_base
    use module_interfaces
    use module_connect_object
    use module_io
    implicit none
    
    private

    public :: connect_recursively
    public :: connect
contains
!=====================================================================
!> This recursive subroutine does not fully support a restart.
!! (subroutine is depracted, use non-recursive version, instead.
!! Before using recursive routine, again, it has to be updated to
!! the full functionality of the non-recursive function
recursive subroutine connect_recursively(mhgpsst,fsw,uinp,runObj,outs,&
                     rcov,nbond,isame,iconnect,rxyz1,rxyz2,ener1,&
                     ener2,fp1,fp2,cobj,connected)
    !if called from outside recursion, connected has to be set 
    !to .true. and nsad=0
    use module_base
    use module_atoms, only: astruct_dump_to_file
    use module_connect_object
    use module_userinput
    use module_mhgps_state
    use module_ls_rmsd
    use module_fingerprints
    use module_minimizers
    use yaml_output
    use module_saddle
    use module_freezingstring
    use module_energyandforces
    use bigdft_run
    implicit none
    !parameters
    type(mhgps_state), intent(inout)      :: mhgpsst
    type(findsad_work), intent(inout)     :: fsw
    type(userinput), intent(in)           :: uinp
    type(run_objects), intent(inout)      :: runObj
    type(state_properties), intent(inout) :: outs
    type(connect_object), intent(inout)   :: cobj
    integer, intent(in)    :: nbond
    real(gp), intent(in)   :: rcov(runObj%atoms%astruct%nat)
    integer, intent(in)    :: iconnect(2,nbond)
    real(gp), intent(in)   :: rxyz1(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)   :: rxyz2(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)   :: fp1(mhgpsst%nid), fp2(mhgpsst%nid)
    real(gp), intent(in)   :: ener1,ener2
    integer, intent(inout) :: isame
    logical, intent(inout) :: connected
    !local
    integer  :: nsad_loc,ipush,infocode
    real(gp) :: displ,ener_count
    real(gp) :: fnoise,fnrm,fmax
    logical  :: converged
    logical  :: lnl, rnr, lnr, rnl 
    character(len=200) :: comment
    real(gp) :: scl

    if(.not.connected)then
        call write_todo(mhgpsst,runObj,outs,rxyz1,rxyz2,ener1,ener2)
        return
    endif
    if(mhgpsst%nsad>=uinp%nsadmax)then
        connected=.false.
        call write_todo(mhgpsst,runObj,outs,rxyz1,rxyz2,ener1,ener2)
        return
    endif

    if(mhgpsst%iproc==0)then
        call yaml_comment('(MHGPS) nsad:'//&
             trim(adjustl(yaml_toa(mhgpsst%nsad)))//'; connect minima with'//&
             ' following energies')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener1)))&
                           //' and ')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener2))))
    endif


    !check if input structures are distinct 
    if(equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,&
                           uinp%fp_delta_min,ener1,ener2,fp1,fp2))then
        if(mhgpsst%iproc==0)call yaml_warning('(MHGPS) connect: '//&
                    'input minima are identical. Will NOT attempt'//&
                    ' to find an intermediate TS. recursion depth: '//&
                    yaml_toa(mhgpsst%nsad))
        return
    endif

    call vcopy(3*runObj%atoms%astruct%nat,rxyz1(1,1),1,cobj%rxyz1(1,1),1)
    call vcopy(3*runObj%atoms%astruct%nat,rxyz2(1,1),1,cobj%rxyz2(1,1),1)
    !rmsd alignment (optional in mhgps approach)
    call superimpose(runObj%atoms%astruct%nat,cobj%rxyz1,cobj%rxyz2)

    !get input guess for transition state
    mhgpsst%nsad=mhgpsst%nsad+1
    nsad_loc=mhgpsst%nsad
    mhgpsst%isad=mhgpsst%isad+1
    write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad

    runObj%inputs%inputPsiId=0
    call get_ts_guess(mhgpsst,uinp,runObj,outs,cobj%rxyz1,cobj%rxyz2,&
          cobj%saddle(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),cobj%tsgenergy,&
          cobj%tsgforces(1,1))


    !compute saddle
    ener_count=0.0_gp
    displ=0.0_gp
    converged=.false.
    call findsad(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,iconnect,&
         cobj%saddle(1,1,mhgpsst%nsad),cobj%enersad(mhgpsst%nsad),&
         cobj%fsad(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),displ,ener_count,&
         cobj%rotforce(1,1,mhgpsst%nsad),converged)

    if(.not.converged)then
        mhgpsst%nsad=mhgpsst%nsad-1!in case we don't want to STOP
        mhgpsst%isad=mhgpsst%isad-1
        write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
        connected=.false.
        call write_todo(mhgpsst,runObj,outs,rxyz1,rxyz2,ener1,ener2)
        call yaml_warning('(MHGPS) Saddle search not converged. '//&
             'Aborting connection attempt.')
!        stop 'STOP saddle not converged'
              return
    endif

    call fnrmandforcemax(cobj%fsad(1,1,mhgpsst%nsad),fnrm,fmax,&
         runObj%atoms%astruct%nat)
    fnrm=sqrt(fnrm)
    if (mhgpsst%iproc == 0) then
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'ATTENTION! Forces '//&
        'below give no forces, but the final minmode| fnrm, fmax = ',&
        fnrm,fmax

        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//&
             '_finalM',comment,cobj%enersad(mhgpsst%nsad),&
             rxyz=cobj%saddle(:,:,mhgpsst%nsad),forces=cobj%minmode(:,:,mhgpsst%nsad))

        write(comment,'(a,1pe10.3,5x,1pe10.3)')&
                                            'fnrm, fmax = ',fnrm,fmax

        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//&
             '_finalF',comment,cobj%enersad(mhgpsst%nsad),&
             rxyz=cobj%saddle(:,:,mhgpsst%nsad),forces=cobj%fsad(:,:,mhgpsst%nsad))

        call write_mode(runObj,outs,mhgpsst%currDir//'/sad'//&
             trim(adjustl(mhgpsst%isadc))//'_mode_final',&
             cobj%minmode(1,1,mhgpsst%nsad),cobj%rotforce(1,1,mhgpsst%nsad))
    endif


    call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,&
         runObj%atoms%astruct%cell_dim,bigdft_get_geocode(runObj),&
         rcov,cobj%saddle(1,1,mhgpsst%nsad),cobj%fpsad(1,mhgpsst%nsad))

    if(mhgpsst%nsad>1)then
        if(equal(mhgpsst%iproc,'(MHGPS)','SS',mhgpsst%nid,uinp%en_delta_sad,&
          uinp%fp_delta_sad,cobj%enersad(mhgpsst%nsad-1),cobj%enersad(mhgpsst%nsad),&
          cobj%fpsad(1,mhgpsst%nsad-1),cobj%fpsad(1,mhgpsst%nsad)))then
            isame=isame+1
            !if we find the same saddle point there times
            !(cosecutively) we are stuck
            if(isame==3)then
                call write_todo(mhgpsst,runObj,outs,rxyz1,rxyz2,&
                     ener1,ener2)
                connected=.false.
                mhgpsst%nsad=mhgpsst%nsad-1
                mhgpsst%isad=mhgpsst%isad-1
                write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
    
                if(mhgpsst%iproc==0)then
                    call yaml_warning('(MHGPS) found same saddle '//&
                                'point again. Aborting connection'//&
                                ' attempt.')
                endif
                isame=0

                return
            endif
        else
            isame=0
        endif
    endif

    scl=-1.0_gp
    ipush=1
    loopL: do
        if(mhgpsst%iproc==0)&
        call yaml_comment('(MHGPS) Relax from left side ',hfill='.')
    
        call pushoffsingle(uinp,runObj%atoms%astruct%nat,&
             cobj%saddle(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),scl,&
             cobj%leftmin(1,1,mhgpsst%nsad))

        ener_count=0.0_gp
        call mhgpsenergyandforces(mhgpsst,runObj,outs,&
             cobj%leftmin(1,1,mhgpsst%nsad),cobj%fleft(1,1,mhgpsst%nsad),fnoise,&
             cobj%enerleft(mhgpsst%nsad),infocode)

        if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity >= 3)&
             call astruct_dump_to_file(&
                  bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                  '/sad'//trim(adjustl(mhgpsst%isadc))//'_pushL',&
                  comment,cobj%enerleft(mhgpsst%nsad),cobj%leftmin(:,:,mhgpsst%nsad),&
                  cobj%fleft(:,:,mhgpsst%nsad))

        call minimize(mhgpsst,uinp,runObj,outs,nbond,iconnect,&
             cobj%leftmin(1,1,mhgpsst%nsad),cobj%fleft(1,1,mhgpsst%nsad),fnoise,&
             cobj%enerleft(mhgpsst%nsad),ener_count,converged,'L')
        call fnrmandforcemax(cobj%fleft(1,1,mhgpsst%nsad),fnrm,fmax,&
             runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(mhgpsst%iproc==0)&
             call astruct_dump_to_file(&
                  bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                  '/sad'//trim(adjustl(mhgpsst%isadc))//'_minFinalL',&
                  comment,cobj%enerleft(mhgpsst%nsad),cobj%leftmin(:,:,mhgpsst%nsad),&
                  cobj%fleft(:,:,mhgpsst%nsad))

        call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,&
             runObj%atoms%astruct%cell_dim,&
             bigdft_get_geocode(runObj),rcov,cobj%leftmin(1,1,mhgpsst%nsad),&
             cobj%fpleft(1,mhgpsst%nsad))
        if(.not.equal(mhgpsst%iproc,'(MHGPS)','MS',mhgpsst%nid,uinp%en_delta_sad,&
           uinp%fp_delta_sad,cobj%enersad(mhgpsst%nsad),cobj%enerleft(mhgpsst%nsad),&
           cobj%fpsad(1,mhgpsst%nsad),cobj%fpleft(1,mhgpsst%nsad)))then
           exit loopL 
        elseif(ipush>=3)then
            mhgpsst%isadprob=mhgpsst%isadprob+1
            write(mhgpsst%isadprobc,'(i5.5)')mhgpsst%isadprob
            if(mhgpsst%iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                'unknown (stepoff converged back to saddle)'
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_finalM',comment,cobj%enersad(mhgpsst%nsad),&
                     rxyz=cobj%saddle(:,:,mhgpsst%nsad),&
                     forces=cobj%minmode(:,:,mhgpsst%nsad))
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_Reactant',comment,0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_Product',comment,0.0_gp,rxyz=cobj%rxyz2)
            endif
    
            connected=.false.
            mhgpsst%nsad=mhgpsst%nsad-1
            mhgpsst%isad=mhgpsst%isad-1
            write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
    
            if(mhgpsst%iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the left minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            return
        endif
        scl=uinp%saddle_scale_stepoff*scl
        if(mhgpsst%iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, left side '//&
                       'converged back to saddle. Will retry with '//&
                       'increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopL

    scl=1.0_gp
    ipush=1
    loopR: do
        if(mhgpsst%iproc==0)&
        call yaml_comment('(MHGPS) Relax from right side ',hfill='.')
    
        call pushoffsingle(uinp,runObj%atoms%astruct%nat,cobj%saddle(1,1,mhgpsst%nsad),&
        cobj%minmode(1,1,mhgpsst%nsad),scl,cobj%rightmin(1,1,mhgpsst%nsad))

        ener_count=0.0_gp
        call mhgpsenergyandforces(mhgpsst,runObj,outs,cobj%rightmin(1,1,mhgpsst%nsad),&
        cobj%fright(1,1,mhgpsst%nsad),fnoise,cobj%enerright(mhgpsst%nsad),infocode)

        if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity >= 3)&
             call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//'_pushR',&
             comment,&
             cobj%enerright(mhgpsst%nsad),cobj%rightmin(1,1,mhgpsst%nsad),&
             cobj%fright(1,1,mhgpsst%nsad))

        call minimize(mhgpsst,uinp,runObj,outs,nbond,iconnect,&
                            cobj%rightmin(1,1,mhgpsst%nsad),&
                            cobj%fright(1,1,mhgpsst%nsad),fnoise,&
                            cobj%enerright(mhgpsst%nsad),ener_count,&
                            converged,'R')
        call fnrmandforcemax(cobj%fright(1,1,mhgpsst%nsad),fnrm,fmax,runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(mhgpsst%iproc==0)&
             call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//'_minFinalR',&
             comment,&
             cobj%enerright(mhgpsst%nsad),cobj%rightmin(1,1,mhgpsst%nsad),&
             cobj%fright(1,1,mhgpsst%nsad))
        call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,runObj%atoms%astruct%cell_dim,bigdft_get_geocode(runObj),rcov,&
                        cobj%rightmin(1,1,mhgpsst%nsad),cobj%fpright(1,mhgpsst%nsad))
        if(.not.equal(mhgpsst%iproc,'(MHGPS)','MS',mhgpsst%nid,uinp%en_delta_sad,uinp%fp_delta_sad,&
           cobj%enersad(mhgpsst%nsad),cobj%enerright(mhgpsst%nsad),&
           cobj%fpsad(1,mhgpsst%nsad),cobj%fpright(1,mhgpsst%nsad)))then
           exit loopR 
        elseif(ipush>=3)then
            mhgpsst%isadprob=mhgpsst%isadprob+1
            write(mhgpsst%isadprobc,'(i5.5)')mhgpsst%isadprob
            if(mhgpsst%iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                     'unknown (stepoff converged back to saddle)'
                        
                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_finalM',&
                     comment,&
                     cobj%enersad(mhgpsst%nsad),cobj%saddle(:,:,mhgpsst%nsad),&
                     forces=cobj%minmode(:,:,mhgpsst%nsad))
                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_Reactant',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_Product',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz2)

            endif
    
            connected=.false.
            mhgpsst%nsad=mhgpsst%nsad-1
            mhgpsst%isad=mhgpsst%isad-1
            write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
    
            if(mhgpsst%iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the right minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            return
        endif
        scl=uinp%saddle_scale_stepoff*scl
        if(mhgpsst%iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, right side'//&
                       ' converged back to saddle. Will retry with'//&
                       ' increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopR

    !is minimum, obtained by relaxation from left bar end identical to
    !left input minimum?
    lnl=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,&
        uinp%en_delta_min,uinp%fp_delta_min,ener1,&
        cobj%enerleft(mhgpsst%nsad),fp1,cobj%fpleft(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !right input minimum?
    rnr=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,&
        uinp%en_delta_min,uinp%fp_delta_min,ener2,&
        cobj%enerright(mhgpsst%nsad),fp2,cobj%fpright(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from left bar end identical to 
    !right input minimum?
    lnr=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,&
        uinp%en_delta_min,uinp%fp_delta_min,ener2,&
        cobj%enerleft(mhgpsst%nsad),fp2,cobj%fpleft(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !left input minimum?
    rnl=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,&
        uinp%en_delta_min,uinp%fp_delta_min,ener1,&
        cobj%enerright(mhgpsst%nsad),fp1,cobj%fpright(1,mhgpsst%nsad))

    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
        if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                           ' connection check connected',&
                            cobj%enerleft(mhgpsst%nsad),&
                            cobj%enerright(mhgpsst%nsad)
        connected=.true.
        return
    endif

    if(lnl .and. (.not. rnr))then
        !connect right input min with right relaxed bar-end
        if(mhgpsst%iproc==0)write(*,*)'(MHGPS) connection check lnl'//&
                           ' and not rnr',sqrt(sum((&
                           rxyz2-cobj%rightmin(:,:,nsad_loc))**2))
        if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')&
                            '(MHGPS) connection check connected',&
                            cobj%enerleft(mhgpsst%nsad),&
                            cobj%enerright(mhgpsst%nsad)
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,cobj,connected)
        return
    endif

    if(rnr .and. (.not. lnl))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check rnr and'//&
                                  ' not lnl',rnr,lnl
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check rnr and'//&
                                  ' not lnl',sqrt(sum((rxyz1-&
                                  cobj%leftmin(:,:,nsad_loc))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                                ' connection check connected',&
                                cobj%enerleft(mhgpsst%nsad),&
                                cobj%enerright(mhgpsst%nsad)
!write(*,*)rxyz1
!write(*,*)
!write(*,*)cobj%leftmin(:,:,nsad_loc)
!if(sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))<1.d-2)then
!    lnl=equal(mhgpsst%iproc,'(MHGPS)',mhgpsst%nid,en_delta_min,uinp%fp_delta_min,ener1,&
!        cobj%enerleft(nsad),fp1,cobj%fpleft(1,nsad))
!!call fingerprint(nat,mhgpsst%nid,alat,atoms%astruct%geocode,rcov,&
!!                cobj%leftmin(1,1,nsad_loc),cobj%fpleft(1,nsad_loc))
!!call fingerprint(nat,mhgpsst%nid,alat,atoms%astruct%geocode,rcov,&
!!                rxyz1(1,1),cobj%fpright(1,nsad_loc))
!write(*,*)'fprints:'
!write(*,*)cobj%fpleft(:,nsad_loc)
!write(*,*)
!write(*,*)fp1
!
!    stop
!endif
        !connect left relaxed bar end with left input min
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),cobj,connected)
        return
    endif

    if(lnr .and. (.not. rnl))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check lnr and'//&
                                  ' not rnl',sqrt(sum((rxyz1-&
                                  cobj%rightmin(:,:,nsad_loc))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                             ' connection check connected',&
                             cobj%enerleft(mhgpsst%nsad),&
                             cobj%enerright(mhgpsst%nsad)
        !connect right relaxed bar end with left input min
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%rightmin(1,1,nsad_loc),&
                     ener1,cobj%enerright(nsad_loc),&
                     fp1,cobj%fpright(1,nsad_loc),cobj,connected)
        return
    endif

    if(.not. lnr .and. rnl)then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnr and'//&
                                  ' rnl',sqrt(sum((rxyz2-&
                                  cobj%leftmin(:,:,nsad_loc))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                        ' connection check connected',&
                        cobj%enerleft(mhgpsst%nsad),cobj%enerright(mhgpsst%nsad)
        !connect left relaxed bar end with right input min
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,rxyz2,cobj%leftmin(1,1,nsad_loc),&
                     ener2,cobj%enerleft(nsad_loc),&
                     fp2,cobj%fpleft(1,nsad_loc),cobj,connected)
        return
    endif

    if((.not. lnl) .and. (.not. rnr))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnl'//&
                       ' and not rnr',sqrt(sum((rxyz1-&
                       cobj%leftmin(:,:,nsad_loc))**2))
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnl'//&
                        ' and not rnr',sqrt(sum((rxyz2-&
                        cobj%rightmin(:,:,nsad_loc))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                        ' connection check connected',&
                        cobj%enerleft(mhgpsst%nsad),&
                        cobj%enerright(mhgpsst%nsad)
        !connect left input min with left relaxed bar end  and right
        !input min with right relaxed bar end
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
                     ener1,cobj%enerleft(nsad_loc),&
                     fp1,cobj%fpleft(1,nsad_loc),cobj,connected)
        call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,isame,&
                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
                     cobj%enerright(nsad_loc),ener2,&
                     cobj%fpright(1,nsad_loc),fp2,cobj,connected)
        return
    endif

    !should and must not happen:
    if(mhgpsst%iproc==0)&
    call yaml_warning('(MHGPS) Severe error in connect: none of '//&
                      'the checks in connect subroutine were '//&
                       'successful! STOP') 
    call f_err_throw('(MHGPS) Severe error in connect: none of '//&
         'the checks in connect subroutine were '//&
         'successful! STOP')

end subroutine
!=====================================================================
subroutine connect(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,&
                     iconnect,rxyz1,rxyz2,ener1,ener2,fp1,fp2,&
                     cobj,connected,premature_exit,nsad)
    use module_base
    use module_atoms, only: astruct_dump_to_file
    use module_connect_object
    use module_io
    use module_userinput
    use module_mhgps_state
    use module_ls_rmsd
    use module_fingerprints
    use module_minimizers
    use yaml_output
    use module_saddle
    use module_freezingstring
    use module_energyandforces
    use bigdft_run
    implicit none
    !parameters
    type(mhgps_state), intent(inout)      :: mhgpsst
    type(findsad_work), intent(inout)      :: fsw
    type(userinput), intent(in) :: uinp
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    integer, intent(in)     :: nbond
    real(gp), intent(in)    :: rcov(runObj%atoms%astruct%nat)
    integer, intent(in)     :: iconnect(2,nbond)
    real(gp), intent(in)    :: rxyz1(3,runObj%atoms%astruct%nat), rxyz2(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)    :: fp1(mhgpsst%nid), fp2(mhgpsst%nid)
    real(gp), intent(in)    :: ener1,ener2
    type(connect_object), intent(inout) :: cobj
    logical, intent(inout)    :: connected
    logical, intent(out)      :: premature_exit
    integer, intent(out)      :: nsad
    !local
    integer  :: infocode
    real(gp) :: displ,ener_count
    real(gp) :: fnoise,fnrm,fmax
    logical  :: converged
    logical  :: lnl, rnr, lnr, rnl 
    character(len=200) :: comment
    integer :: ntodo
    real(gp) :: fp1cur(mhgpsst%nid), fp2cur(mhgpsst%nid)
    real(gp) :: ener1cur, ener2cur
    integer :: isame
    integer :: ipush
    integer :: iloop
    real(gp) :: scl

    connected=.false.
    premature_exit = .false.
!    mhgpsst%nsad=0
    isame=0

    cobj%ntodo=1
    cobj%todorxyz(:,:,1,cobj%ntodo)=rxyz1
    cobj%todorxyz(:,:,2,cobj%ntodo)=rxyz2
    cobj%todofp(:,1,cobj%ntodo)=fp1
    cobj%todofp(:,2,cobj%ntodo)=fp2
    cobj%todoenergy(1,cobj%ntodo)=ener1
    cobj%todoenergy(2,cobj%ntodo)=ener2

    iloop=0
connectloop: do while(cobj%ntodo>=1)
    iloop=iloop+1
    if(mhgpsst%iproc==0)then
        call write_restart(mhgpsst,runObj,cobj)
    endif
    if(mhgpsst%nsad>=uinp%nsadmax)then
        connected=.false.
        exit connectloop
    endif
    !following if must be AFTER check for
    !mhgpsst%nsad>=uinp%nsadmax
    if(iloop>1 .and. uinp%singlestep)then
        premature_exit=.true.
        exit connectloop
    endif

    call vcopy(3*runObj%atoms%astruct%nat,cobj%todorxyz(1,1,1,cobj%ntodo),1,cobj%rxyz1(1,1), 1)
    call vcopy(mhgpsst%nid,cobj%todofp(1,1,cobj%ntodo),1,fp1cur(1), 1)
    ener1cur=cobj%todoenergy(1,cobj%ntodo)
    call vcopy(3*runObj%atoms%astruct%nat,cobj%todorxyz(1,1,2,cobj%ntodo),1,cobj%rxyz2(1,1), 1)
    call vcopy(mhgpsst%nid,cobj%todofp(1,2,cobj%ntodo),1,fp2cur(1), 1)
    ener2cur=cobj%todoenergy(2,cobj%ntodo)
    if(mhgpsst%iproc==0)then
        call yaml_comment('(MHGPS) nsad:'//&
             trim(adjustl(yaml_toa(mhgpsst%nsad)))//'; connect minima with'//&
             ' following energies')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener1cur)))&
                           //' and ')
        call yaml_comment('(MHGPS) '//trim(adjustl(yaml_toa(ener2cur))))
    endif

    !check if input structures are distinct 
    if(equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,&
       uinp%fp_delta_min,cobj%todoenergy(1,cobj%ntodo),cobj%todoenergy(2,cobj%ntodo),&
                             cobj%todofp(1,1,cobj%ntodo),cobj%todofp(1,2,cobj%ntodo)))then
        if(mhgpsst%iproc==0)call yaml_warning('(MHGPS) connect: '//&
                    'minima are identical. Will NOT attempt to find '//&
                    'an intermediate TS.')
!        connected=.true.
        cobj%ntodo=cobj%ntodo-1
        cycle
    endif
    if(bigdft_get_geocode(runObj)=='F')then
    !rmsd alignment (optional in mhgps approach)
    call superimpose(runObj%atoms%astruct%nat,cobj%rxyz1,cobj%rxyz2)

    !check if previously connected
    if(previously_connected(mhgpsst,uinp,runObj,cobj%rxyz1,cobj%rxyz2))then
        if(mhgpsst%iproc==0)call yaml_comment('(MHGPS) connect: '//&
                    'Minima previously connected. Will not connect again.')
!        connected=.true.
        cobj%ntodo=cobj%ntodo-1
        cycle
    endif
    endif

    !get input guess for transition state
    mhgpsst%nsad=mhgpsst%nsad+1
    mhgpsst%isad=mhgpsst%isad+1
    write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad

    runObj%inputs%inputPsiId=0
    call get_ts_guess(mhgpsst,uinp,runObj,outs,cobj%rxyz1,cobj%rxyz2,&
          cobj%saddle(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),cobj%tsgenergy,&
          cobj%tsgforces(1,1))

    !compute saddle
    ener_count=0.0_gp
    displ=0.0_gp
    converged=.false.
    call findsad(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,iconnect,&
         cobj%saddle(1,1,mhgpsst%nsad),cobj%enersad(mhgpsst%nsad),&
         cobj%fsad(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),displ,ener_count,&
         cobj%rotforce(1,1,mhgpsst%nsad),converged)

    if(.not.converged)then
        mhgpsst%nsad=mhgpsst%nsad-1!in case we don't want to STOP
        mhgpsst%isad=mhgpsst%isad-1
        write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
        connected=.false.
        call yaml_warning('(MHGPS) Saddle search not converged. '//&
             'Aborting connection attempt.')
!        stop 'STOP saddle not converged'
        exit connectloop
    endif

    call fnrmandforcemax(cobj%fsad(1,1,mhgpsst%nsad),fnrm,fmax,&
         runObj%atoms%astruct%nat)
    fnrm=sqrt(fnrm)
    if (mhgpsst%iproc == 0) then
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'ATTENTION! Forces '//&
        'below give no forces, but the final minmode| fnrm, fmax = ',&
        fnrm,fmax

        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//&
             '_finalM',comment,cobj%enersad(mhgpsst%nsad),&
             rxyz=cobj%saddle(:,:,mhgpsst%nsad),forces=cobj%minmode(:,:,mhgpsst%nsad))

        write(comment,'(a,1pe10.3,5x,1pe10.3)')&
                                            'fnrm, fmax = ',fnrm,fmax

        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//&
             '_finalF',comment,cobj%enersad(mhgpsst%nsad),&
             rxyz=cobj%saddle(:,:,mhgpsst%nsad),forces=cobj%fsad(:,:,mhgpsst%nsad))

        call write_mode(runObj,outs,mhgpsst%currDir//'/sad'//&
             trim(adjustl(mhgpsst%isadc))//'_mode_final',&
             cobj%minmode(1,1,mhgpsst%nsad),cobj%rotforce(1,1,mhgpsst%nsad))
    endif


    call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,&
         runObj%atoms%astruct%cell_dim,bigdft_get_geocode(runObj),&
         rcov,cobj%saddle(1,1,mhgpsst%nsad),cobj%fpsad(1,mhgpsst%nsad))

    if(mhgpsst%nsad>1 .and. (.not. uinp%singlestep))then
        if(equal(mhgpsst%iproc,'(MHGPS)','SS',mhgpsst%nid,uinp%en_delta_sad,&
          uinp%fp_delta_sad,cobj%enersad(mhgpsst%nsad-1),cobj%enersad(mhgpsst%nsad),&
          cobj%fpsad(1,mhgpsst%nsad-1),cobj%fpsad(1,mhgpsst%nsad)))then
            isame=isame+1
            !if we find the same saddle point there times
            !(cosecutively) we are stuck
            !its three, because:
            !
            !     sad1       relax              sad1
            !               ------>             /  \
            !min1      min2         min1     min3  min4    min2
            !now compute saddle between min4 and min2, if sad1 is
            !found again, it is tried to connect min3 and min2.
            !If now again sad1 is found, min4 and min2 was tried again
            !if we would not stop here.
            if(isame==3)then
                if(mhgpsst%iproc==0)then
                    call yaml_warning('(MHGPS) found same saddle '//&
                                'point again. Aborting connection'//&
                                ' attempt.')
                endif
                call write_todo(mhgpsst,runObj,outs,rxyz1,rxyz2,&
                     ener1,ener2)
                connected=.false.
                mhgpsst%nsad=mhgpsst%nsad-1
                mhgpsst%isad=mhgpsst%isad-1
                write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
                isame=0
                exit connectloop !stop connection
            endif
        else
            isame=0
        endif
    endif

    scl=-1.0_gp
    ipush=1
    loopL: do
        if(mhgpsst%iproc==0)&
        call yaml_comment('(MHGPS) Relax from left side ',hfill='.')

        call pushoffsingle(uinp,runObj%atoms%astruct%nat,&
             cobj%saddle(1,1,mhgpsst%nsad),cobj%minmode(1,1,mhgpsst%nsad),scl,&
             cobj%leftmin(1,1,mhgpsst%nsad))

        ener_count=0.0_gp
        call mhgpsenergyandforces(mhgpsst,runObj,outs,&
             cobj%leftmin(1,1,mhgpsst%nsad),cobj%fleft(1,1,mhgpsst%nsad),fnoise,&
             cobj%enerleft(mhgpsst%nsad),infocode)

        if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity >= 3)&
             call astruct_dump_to_file(&
                  bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                  '/sad'//trim(adjustl(mhgpsst%isadc))//'_pushL',&
                  comment,cobj%enerleft(mhgpsst%nsad),cobj%leftmin(:,:,mhgpsst%nsad),&
                  cobj%fleft(:,:,mhgpsst%nsad))

        call minimize(mhgpsst,uinp,runObj,outs,nbond,iconnect,&
             cobj%leftmin(1,1,mhgpsst%nsad),cobj%fleft(1,1,mhgpsst%nsad),fnoise,&
             cobj%enerleft(mhgpsst%nsad),ener_count,converged,'L')
        call fnrmandforcemax(cobj%fleft(1,1,mhgpsst%nsad),fnrm,fmax,&
             runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(mhgpsst%iproc==0)&
             call astruct_dump_to_file(&
                  bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                  '/sad'//trim(adjustl(mhgpsst%isadc))//'_minFinalL',&
                  comment,cobj%enerleft(mhgpsst%nsad),cobj%leftmin(:,:,mhgpsst%nsad),&
                  cobj%fleft(:,:,mhgpsst%nsad))

        call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,&
             runObj%atoms%astruct%cell_dim,&
             bigdft_get_geocode(runObj),rcov,cobj%leftmin(1,1,mhgpsst%nsad),&
             cobj%fpleft(1,mhgpsst%nsad))
        if(.not.equal(mhgpsst%iproc,'(MHGPS)','MS',mhgpsst%nid,uinp%en_delta_sad,&
           uinp%fp_delta_sad,cobj%enersad(mhgpsst%nsad),cobj%enerleft(mhgpsst%nsad),&
           cobj%fpsad(1,mhgpsst%nsad),cobj%fpleft(1,mhgpsst%nsad)))then
           exit loopL
        elseif(ipush>=3)then
            mhgpsst%isadprob=mhgpsst%isadprob+1
            write(mhgpsst%isadprobc,'(i5.5)')mhgpsst%isadprob
            if(mhgpsst%iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                'unknown (stepoff converged back to saddle)'
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_finalM',comment,cobj%enersad(mhgpsst%nsad),&
                     rxyz=cobj%saddle(:,:,mhgpsst%nsad),&
                     forces=cobj%minmode(:,:,mhgpsst%nsad))
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_Reactant',comment,0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(&
                     bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                     '/sadProb'//trim(adjustl(mhgpsst%isadprobc))//&
                     '_Product',comment,0.0_gp,rxyz=cobj%rxyz2)
            endif

            connected=.false.
            mhgpsst%nsad=mhgpsst%nsad-1
            mhgpsst%isad=mhgpsst%isad-1
            write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad

            if(mhgpsst%iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the left minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            
            exit connectloop !stop connection
        endif
        scl=uinp%saddle_scale_stepoff*scl
        if(mhgpsst%iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, left side '//&
                       'converged back to saddle. Will retry with '//&
                       'increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopL

    scl=1.0_gp
    ipush=1
    loopR: do
        if(mhgpsst%iproc==0)&
        call yaml_comment('(MHGPS) Relax from right side ',hfill='.')

        call pushoffsingle(uinp,runObj%atoms%astruct%nat,cobj%saddle(1,1,mhgpsst%nsad),&
        cobj%minmode(1,1,mhgpsst%nsad),scl,cobj%rightmin(1,1,mhgpsst%nsad))

        ener_count=0.0_gp
        call mhgpsenergyandforces(mhgpsst,runObj,outs,cobj%rightmin(1,1,mhgpsst%nsad),&
        cobj%fright(1,1,mhgpsst%nsad),fnoise,cobj%enerright(mhgpsst%nsad),infocode)

        if(mhgpsst%iproc==0 .and. uinp%mhgps_verbosity >= 3)&
             call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//'_pushR',&
             comment,&
             cobj%enerright(mhgpsst%nsad),cobj%rightmin(1,1,mhgpsst%nsad),&
             cobj%fright(1,1,mhgpsst%nsad))

        call minimize(mhgpsst,uinp,runObj,outs,nbond,iconnect,&
                            cobj%rightmin(1,1,mhgpsst%nsad),&
                            cobj%fright(1,1,mhgpsst%nsad),fnoise,&
                            cobj%enerright(mhgpsst%nsad),ener_count,&
                            converged,'R')
        call fnrmandforcemax(cobj%fright(1,1,mhgpsst%nsad),fnrm,fmax,runObj%atoms%astruct%nat)
        fnrm=sqrt(fnrm)
        write(comment,'(a,1pe10.3,5x,1pe10.3)')'fnrm, fmax = ',fnrm,&
                                              fmax
        if(mhgpsst%iproc==0)&
             call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/sad'//trim(adjustl(mhgpsst%isadc))//'_minFinalR',&
             comment,&
             cobj%enerright(mhgpsst%nsad),cobj%rightmin(1,1,mhgpsst%nsad),&
             cobj%fright(1,1,mhgpsst%nsad))
        call fingerprint(runObj%atoms%astruct%nat,mhgpsst%nid,runObj%atoms%astruct%cell_dim,bigdft_get_geocode(runObj),rcov,&
                        cobj%rightmin(1,1,mhgpsst%nsad),cobj%fpright(1,mhgpsst%nsad))
        if(.not.equal(mhgpsst%iproc,'(MHGPS)','MS',mhgpsst%nid,uinp%en_delta_sad,uinp%fp_delta_sad,&
           cobj%enersad(mhgpsst%nsad),cobj%enerright(mhgpsst%nsad),&
           cobj%fpsad(1,mhgpsst%nsad),cobj%fpright(1,mhgpsst%nsad)))then
           exit loopR
        elseif(ipush>=3)then
            mhgpsst%isadprob=mhgpsst%isadprob+1
            write(mhgpsst%isadprobc,'(i5.5)')mhgpsst%isadprob
            if(mhgpsst%iproc==0)then
                write(comment,'(a)')'Prob: Neighbors '//&
                     'unknown (stepoff converged back to saddle)'

                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_finalM',&
                     comment,&
                     cobj%enersad(mhgpsst%nsad),cobj%saddle(:,:,mhgpsst%nsad),&
                     forces=cobj%minmode(:,:,mhgpsst%nsad))
                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_Reactant',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz1)
                call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                     mhgpsst%currDir//'/sadProb'//trim(adjustl(mhgpsst%isadprobc))//'_Product',&
                     comment,&
                0.0_gp,rxyz=cobj%rxyz2)

            endif

            connected=.false.
            mhgpsst%nsad=mhgpsst%nsad-1
            mhgpsst%isad=mhgpsst%isad-1
            write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad

            if(mhgpsst%iproc==0)&
            call yaml_warning('(MHGPS)  after relaxation from '//&
                              'saddle point the right minimum is '//&
                              'identical to the saddle point. '//&
                              'Stopped connection attempt. Will '//&
                              'proceed with next connection attempt.')
            exit connectloop !stop connection
        endif
        scl=uinp%saddle_scale_stepoff*scl
        if(mhgpsst%iproc==0)&
        call yaml_comment('INFO: (MHGPS) After pushoff, right side'//&
                       ' converged back to saddle. Will retry with'//&
                       ' increased pushoff: '//&
                        yaml_toa(scl))
        ipush=ipush+1
    enddo loopR

    !one more saddle point is done
    cobj%ntodo=cobj%ntodo-1

    !is minimum, obtained by relaxation from left bar end identical to
    !left input minimum?
    lnl=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,uinp%fp_delta_min,ener1cur,&
        cobj%enerleft(mhgpsst%nsad),fp1cur,cobj%fpleft(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !right input minimum?
    rnr=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,uinp%fp_delta_min,ener2cur,&
        cobj%enerright(mhgpsst%nsad),fp2cur,cobj%fpright(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from left bar end identical to 
    !right input minimum?
    lnr=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,uinp%fp_delta_min,ener2cur,&
        cobj%enerleft(mhgpsst%nsad),fp2cur,cobj%fpleft(1,mhgpsst%nsad))

    !is minimum obtained by relaxation from right bar end identical to
    !left input minimum?
    rnl=equal(mhgpsst%iproc,'(MHGPS)','MM',mhgpsst%nid,uinp%en_delta_min,uinp%fp_delta_min,ener1cur,&
        cobj%enerright(mhgpsst%nsad),fp1cur,cobj%fpright(1,mhgpsst%nsad))

    if((lnl .and. rnr) .or. (lnr .and. rnl))then!connection done
        if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')&
                            '(MHGPS) connection check connected',&
                            cobj%enerleft(mhgpsst%nsad),&
                            cobj%enerright(mhgpsst%nsad)
!        connected=.true.
!        return
cycle
    elseif(lnl .and. (.not. rnr))then
        !connect right input min with right relaxed bar-end
        if(mhgpsst%iproc==0)write(*,*)'(MHGPS) connection check'//&
                            ' lnl and not rnr',sqrt(sum((rxyz2-&
                            cobj%rightmin(:,:,mhgpsst%nsad))**2))
        if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')&
                           '(MHGPS) connection check connected',&
                           cobj%enerleft(mhgpsst%nsad),&
                           cobj%enerright(mhgpsst%nsad)
        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rightmin(:,:,mhgpsst%nsad)
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%rxyz2
        cobj%todofp(:,1,cobj%ntodo)=cobj%fpright(:,mhgpsst%nsad)
        cobj%todofp(:,2,cobj%ntodo)=fp2cur
        cobj%todoenergy(1,cobj%ntodo)=cobj%enerright(mhgpsst%nsad)
        cobj%todoenergy(2,cobj%ntodo)=ener2cur
        
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
!                     cobj%enerright(nsad_loc),ener2,&
!                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
!        return
cycle
    elseif(rnr .and. (.not. lnl))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check rnr and'//&
                       ' not lnl',rnr,lnl
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check rnr and'//&
                       ' not lnl',sqrt(sum((cobj%rxyz1-&
                       cobj%leftmin(:,:,mhgpsst%nsad))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                        ' connection check connected',&
                        cobj%enerleft(mhgpsst%nsad),&
                        cobj%enerright(mhgpsst%nsad)
!write(*,*)rxyz1
!write(*,*)
!write(*,*)cobj%leftmin(:,:,nsad_loc)
!if(sqrt(sum((rxyz1-cobj%leftmin(:,:,nsad_loc))**2))<1.d-2)then
!    lnl=equal(mhgpsst%iproc,'(MHGPS)',nid,en_delta_min,fp_delta_min,ener1,&
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
        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rxyz1
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%leftmin(:,:,mhgpsst%nsad)
        cobj%todofp(:,1,cobj%ntodo)=fp1cur
        cobj%todofp(:,2,cobj%ntodo)=cobj%fpleft(:,mhgpsst%nsad)
        cobj%todoenergy(1,cobj%ntodo)=ener1cur
        cobj%todoenergy(2,cobj%ntodo)=cobj%enerleft(mhgpsst%nsad)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
!                     ener1,cobj%enerleft(nsad_loc),&
!                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        return
cycle
    elseif(lnr .and. (.not. rnl))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check lnr and'//&
                        ' not rnl',sqrt(sum((cobj%rxyz1-&
                        cobj%rightmin(:,:,mhgpsst%nsad))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                        ' connection check connected',&
                        cobj%enerleft(mhgpsst%nsad),&
                        cobj%enerright(mhgpsst%nsad)
        !connect right relaxed bar end with left input min
        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rxyz1
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%rightmin(:,:,mhgpsst%nsad)
        cobj%todofp(:,1,cobj%ntodo)=fp1cur
        cobj%todofp(:,2,cobj%ntodo)=cobj%fpright(:,mhgpsst%nsad)
        cobj%todoenergy(1,cobj%ntodo)=ener1cur
        cobj%todoenergy(2,cobj%ntodo)=cobj%enerright(mhgpsst%nsad)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%rightmin(1,1,nsad_loc),&
!                     ener1,cobj%enerright(nsad_loc),&
!                     fp1,cobj%fpright(1,nsad_loc),nsad,cobj,connected)
!        return
cycle
    elseif(.not. lnr .and. rnl)then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnr'//&
                       ' and rnl',sqrt(sum((cobj%rxyz2-&
                       cobj%leftmin(:,:,mhgpsst%nsad))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS) '//&
                       'connection check connected',&
                       cobj%enerleft(mhgpsst%nsad),&
                       cobj%enerright(mhgpsst%nsad)
        !connect left relaxed bar end with right input min
        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rxyz2
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%leftmin(:,:,mhgpsst%nsad)
        cobj%todofp(:,1,cobj%ntodo)=fp2cur
        cobj%todofp(:,2,cobj%ntodo)=cobj%fpleft(:,mhgpsst%nsad)
        cobj%todoenergy(1,cobj%ntodo)=ener2cur
        cobj%todoenergy(2,cobj%ntodo)=cobj%enerleft(mhgpsst%nsad)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz2,cobj%leftmin(1,1,nsad_loc),&
!                     ener2,cobj%enerleft(nsad_loc),&
!                     fp2,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        return
cycle
    elseif((.not. lnl) .and. (.not. rnr))then
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnl'//&
                        ' and not rnr',sqrt(sum((cobj%rxyz1-&
                        cobj%leftmin(:,:,mhgpsst%nsad))**2))
    if(mhgpsst%iproc==0)write(*,*)'(MHGPS)connection check not lnl'//&
                        ' and not rnr',sqrt(sum((cobj%rxyz2-&
                        cobj%rightmin(:,:,mhgpsst%nsad))**2))
    if(mhgpsst%iproc==0)write(*,'(a,es24.17,1x,es24.17)')'(MHGPS)'//&
                        ' connection check connected',&
                        cobj%enerleft(mhgpsst%nsad),&
                        cobj%enerright(mhgpsst%nsad)
        !connect left input min with left relaxed bar end  and right
        !input min with right relaxed bar end
        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rightmin(:,:,mhgpsst%nsad)
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%rxyz2
        cobj%todofp(:,1,cobj%ntodo)=cobj%fpright(:,mhgpsst%nsad)
        cobj%todofp(:,2,cobj%ntodo)=fp2cur
        cobj%todoenergy(1,cobj%ntodo)=cobj%enerright(mhgpsst%nsad)
        cobj%todoenergy(2,cobj%ntodo)=ener2cur

        cobj%ntodo=cobj%ntodo+1
if(cobj%ntodo>uinp%nsadmax)stop 'error: cobj%ntodo>uinp%nsadmax'
        cobj%todorxyz(:,:,1,cobj%ntodo)=cobj%rxyz1
        cobj%todorxyz(:,:,2,cobj%ntodo)=cobj%leftmin(:,:,mhgpsst%nsad)
        cobj%todofp(:,1,cobj%ntodo)=fp1cur
        cobj%todofp(:,2,cobj%ntodo)=cobj%fpleft(:,mhgpsst%nsad)
        cobj%todoenergy(1,cobj%ntodo)=ener1cur
        cobj%todoenergy(2,cobj%ntodo)=cobj%enerleft(mhgpsst%nsad)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,rxyz1,cobj%leftmin(1,1,nsad_loc),&
!                     ener1,cobj%enerleft(nsad_loc),&
!                     fp1,cobj%fpleft(1,nsad_loc),nsad,cobj,connected)
!        call connect_recursively(nat,nid,alat,rcov,nbond,&
!                     iconnect,cobj%rightmin(1,1,nsad_loc),rxyz2,&
!                     cobj%enerright(nsad_loc),ener2,&
!                     cobj%fpright(1,nsad_loc),fp2,nsad,cobj,connected)
!        return
cycle
    else
        !should and must not happen:
        call f_err_throw('(MHGPS) Severe error in connect: none '//&
                          'of the checks in connect subroutine '//&
                          'were successful! STOP')
    endif
enddo connectloop
nsad=mhgpsst%nsad
if(cobj%ntodo<=0)then
    connected=.true.
endif
if(.not. premature_exit)then
if(connected)then
!only write if connection really connected
!(that is, no premature exit)
!if connected, the write_restart inside the connectloop
!has not been callled a last time.
!Therefore, it has to be done here.
if(cobj%ntodo>=1) stop 'bastian'
    if(mhgpsst%iproc==0)then
        call write_restart(mhgpsst,runObj,cobj)
!        call write_restart(mhgpsst,runObj)
    endif
else
!only write if connection really failed
!(that is, no premature exit)
    call write_todoList(uinp,mhgpsst,runObj,cobj)
    if(mhgpsst%iproc==0)then
        call write_restart(mhgpsst,runObj)
    endif
endif
endif


end subroutine
!=====================================================================
function previously_connected(mhgpsst,uinp,runObj,rxyz1,rxyz2)
    !this function compares rxyz1 an rxyz2 with a non-permutational
    !invariant rmsd if the have been tried to be connected previously.
    !a non-permutationally invariant rmsd is used such that
    !transition states between different premuational variants of the 
    !same structures are computed
    use module_base
    use bigdft_run, only: run_objects
    use module_mhgps_state
    use module_userinput
    use module_ls_rmsd
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(userinput), intent(in)  :: uinp
    type(run_objects), intent(in) :: runObj
    real(gp), intent(in) :: rxyz1(3,runObj%atoms%astruct%nat)
    real(gp), intent(in) :: rxyz2(3,runObj%atoms%astruct%nat)
    logical :: previously_connected
    !local
    real(gp), parameter :: rmsdthresh=0.01_gp
    integer :: iatt
    integer :: i
    logical :: match
    real(gp) :: rmsd1, rmsd2, rmsd3
    real(gp),allocatable :: attempted_connections_tmp(:,:,:,:)
!f_malloc((/1.to.3,1.to.runObj%atoms%astruct%nat,2,mhgpsst%nattemptedmax/),id='mhgpsst%attempted_connections')
!  mhgpsst%nattemptedmax = 1000
!  mhgpsst%nattempted    = 0
!  mhgpsst%attempted_connections

    previously_connected = .false.
    outer: do iatt = 1 , mhgpsst%nattempted
        rmsd1=rmsd(runObj%atoms%astruct%nat,&
              mhgpsst%attempted_connections(1,1,1,iatt),rxyz1)
if(mhgpsst%iproc==0)write(*,*)'rmsd ',rmsd1
        if(rmsd1 <= rmsdthresh)then
            rmsd3=rmsd(runObj%atoms%astruct%nat,&
                  mhgpsst%attempted_connections(1,1,2,iatt),rxyz2)
if(mhgpsst%iproc==0)write(*,*)'rmsd ',rmsd3
            if(rmsd3 <= rmsdthresh) then
                previously_connected = .true.
                exit outer
            endif
        endif
        rmsd2=rmsd(runObj%atoms%astruct%nat,&
              mhgpsst%attempted_connections(1,1,1,iatt),rxyz2)
if(mhgpsst%iproc==0)write(*,*)'rmsd ',rmsd2
        if(rmsd2 <= rmsdthresh)then
            rmsd3=rmsd(runObj%atoms%astruct%nat,&
                  mhgpsst%attempted_connections(1,1,2,iatt),rxyz1)
if(mhgpsst%iproc==0)write(*,*)'rmsd ',rmsd3
            if(rmsd3 <= rmsdthresh) then
                previously_connected = .true.
                exit outer
            endif
        endif
    enddo outer

    if(.not. previously_connected)then
        mhgpsst%nattempted = mhgpsst%nattempted + 1
        !check if enough space, if not
        !resize array
        if(mhgpsst%nattempted > mhgpsst%nattemptedmax)then
if(mhgpsst%iproc==0)write(*,*)'prevresize '

            attempted_connections_tmp = f_malloc((/3,&
                                       runObj%atoms%astruct%nat,2,&
                                       mhgpsst%nattemptedmax/),&
                                       id='attempted_connections_tmp')
            attempted_connections_tmp = mhgpsst%attempted_connections
            call f_free(mhgpsst%attempted_connections)
            mhgpsst%nattemptedmax = mhgpsst%nattemptedmax + 1000
            mhgpsst%attempted_connections = f_malloc((/3,&
                                       runObj%atoms%astruct%nat,2,&
                                       mhgpsst%nattemptedmax/),&
                                       id='mhgpsst%attempted_connections')
            do i = 1, mhgpsst%nattemptedmax - 1000
                mhgpsst%attempted_connections(:,:,:,i) = attempted_connections_tmp(:,:,:,i)
            enddo
            call f_free(attempted_connections_tmp)
        endif

        !add pair to list
        mhgpsst%attempted_connections(:,:,1,mhgpsst%nattempted) &
                = rxyz1
        mhgpsst%attempted_connections(:,:,2,mhgpsst%nattempted) &
                = rxyz2
    endif
end function
!=====================================================================
subroutine pushoff(uinp,nat,saddle,minmode,left,right)
    use module_base
    use module_userinput
    use module_misc
    implicit none
    !parameters 
    type(userinput), intent(in) :: uinp
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(out) :: left(3,nat)
    real(gp), intent(out) :: right(3,nat)
    !local
    real(gp)  :: step(3,nat)


    step = uinp%saddle_stepoff*minmode
    left = saddle - step
    right = saddle + step

end subroutine
!=====================================================================
subroutine pushoffsingle(uinp,nat,saddle,minmode,scl,pushed)
    use module_base, only: gp
    use module_misc
    use module_userinput
    implicit none
    !parameters
    type(userinput)     :: uinp
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in) :: scl
    real(gp), intent(out) :: pushed(3,nat)
    !local
    integer  :: iat 
    real(gp) :: step(3,nat)
    real(gp) :: maxd, tt, dp

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

    step = minmode*uinp%saddle_stepoff/maxd
    pushed = saddle + scl*step

!    !control:
!    tt=0.0_gp
!    dp=0.0_gp
!    maxd=-huge(1.0_gp)
!    do iat=1,nat
!        dp=step(1,iat)**2+step(2,iat)**2+step(3,iat)**2
!        tt=tt+dp
!        maxd=max(maxd,dp)
!    enddo
!    tt=sqrt(tt)
!    maxd=sqrt(maxd)
!    write(133,*)scl, maxd

! old:
!    step = saddle_stepoff*minmode
!    pushed = saddle + scl*step

end subroutine
!=====================================================================
subroutine pushoff_assym(uinp,nat,saddle,minmode,scll,sclr,left,right)
    use module_base
    use module_userinput
    use module_misc
    implicit none
    !parameters 
    type(userinput)     :: uinp
    integer, intent(in) :: nat
    real(gp), intent(in) :: saddle(3,nat)
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in) :: scll,sclr
    real(gp), intent(out) :: left(3,nat)
    real(gp), intent(out) :: right(3,nat)
    !local
    real(gp), dimension(3,nat) :: step

    step =uinp% saddle_stepoff*minmode
    left = saddle - scll*step
    right = saddle + sclr*step
end subroutine
!=====================================================================
subroutine write_todo(mhgpsst,runObj,outs,left,right,eleft,eright)
    use module_base, only: gp
    use module_atoms, only: astruct_dump_to_file
    use bigdft_run, only: bigdft_get_astruct_ptr, run_objects,&
                          state_properties
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    real(gp), intent(in)   :: left(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)   :: right(3,runObj%atoms%astruct%nat)
    real(gp), intent(in)   :: eleft
    real(gp), intent(in)   :: eright
    !local
    character(len=1) :: comment=' '
    
    mhgpsst%ntodo=mhgpsst%ntodo+1
    write(mhgpsst%ntodoc,'(i5.5)')mhgpsst%ntodo

    if(mhgpsst%iproc==0)then
        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/todo'//trim(adjustl(mhgpsst%ntodoc))//'_L',comment,&
             energy=eleft,rxyz=left)
        call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
             mhgpsst%currDir//'/todo'//trim(adjustl(mhgpsst%ntodoc))//'_R',comment,&
             energy=eright,rxyz=right)
    endif
end subroutine
!=====================================================================
subroutine write_todoList(uinp,mhgpsst,runObj,cobj)
    use module_base, only: gp
    use module_userinput
    use module_atoms, only: atomic_structure,&
                            deallocate_atomic_structure,&
                            astruct_dump_to_file,&
                            read_atomic_file=>set_astruct_from_file
    use bigdft_run, only: bigdft_get_astruct_ptr, run_objects,&
                          state_properties
    use module_connect_object
    use module_mhgps_state
    implicit none
    !parameters
    type(userinput), intent(in)     :: uinp
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(inout) :: runObj
    type(connect_object), intent(in) :: cobj
    !local
    integer :: itodo,ijob
    character(len=1) :: comment=' '
    type(atomic_structure):: astruct !< Contains all info
   
    do itodo=cobj%ntodo,1,-1 
        mhgpsst%ntodo=mhgpsst%ntodo+1
        write(mhgpsst%ntodoc,'(i5.5)')mhgpsst%ntodo
        if(mhgpsst%iproc==0)then
            call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                 mhgpsst%currDir//'/todo'//trim(adjustl(mhgpsst%ntodoc))&
                 //'_L',comment,rxyz=cobj%todorxyz(1,1,1,itodo))
            call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                 mhgpsst%currDir//'/todo'//trim(adjustl(mhgpsst%ntodoc))&
                 //'_R',comment,rxyz=cobj%todorxyz(1,1,2,itodo))
        endif
    enddo
    if(uinp%singlestep)then
        do ijob=mhgpsst%ijob+1,mhgpsst%njobs
            if(trim(adjustl(mhgpsst%joblist(1,ijob)(10:16)))=='restart')then
                mhgpsst%ntodo=mhgpsst%ntodo+1
                write(mhgpsst%ntodoc,'(i5.5)')mhgpsst%ntodo
                if(mhgpsst%iproc==0)then
                    call read_atomic_file(trim(adjustl(&
                         mhgpsst%joblist(1,ijob))),mhgpsst%iproc,&
                         astruct)
                    call astruct_dump_to_file(astruct,mhgpsst%currDir//&
                         '/todo'//trim(adjustl(mhgpsst%ntodoc))//'_L',&
                         comment)
                    call deallocate_atomic_structure(astruct)
                    call read_atomic_file(trim(adjustl(&
                         mhgpsst%joblist(2,ijob))),mhgpsst%iproc,&
                         astruct)
                    call astruct_dump_to_file(astruct,mhgpsst%currDir//&
                         '/todo'//trim(adjustl(mhgpsst%ntodoc))//'_R',&
                         comment)
                    call deallocate_atomic_structure(astruct)
                endif
            endif
        enddo
    endif
end subroutine
!=====================================================================
end module

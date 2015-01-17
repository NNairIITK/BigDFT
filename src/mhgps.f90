!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program mhgps
    use module_base
    use module_types
    use module_interfaces
    use module_input_dicts
    use yaml_output
    use module_io
    use module_atoms, only: deallocate_atoms_data,&
                            deallocate_atomic_structure,&
                            atomic_structure,&
                            read_atomic_file=>set_astruct_from_file,&
                            astruct_dump_to_file
    use module_connect_object
    use module_mhgps_state
    use module_userinput
    use module_energyandforces, only: mhgpsenergyandforces
    use module_ls_rmsd, only: superimpose
    use module_sqn, only: findbonds !for finding binds
    use module_freezingstring, only: get_ts_guess 
    use module_saddle, only: findsad_work, findsad,&
                             allocate_finsad_workarrays,&
                             deallocate_finsad_workarrays
    use module_connect, only: connect_recursively, connect
    use module_fingerprints, only: fingerprint
    use module_hessian, only: cal_hessian_fd 
    use module_minimizers
    use bigdft_run
    implicit none
    integer                   :: u
    integer                   :: nfree
    integer                   :: iat
    integer                   :: info
    integer                   :: isame
    integer                   :: nbond
    integer                   :: nsad
    integer                   :: infocode
    integer                   :: ifolder
    integer                   :: ifolderstart
    integer                   :: ijob
    integer                   :: ierr
    integer                   :: lwork
    integer, allocatable      :: iconnect(:,:)
    logical                   :: connected
    logical                   :: premature_exit
    character(len=200)        :: filename
    character(len=60)         :: run_id,naming_id
    real(gp), allocatable     :: rcov(:)
    real(8), allocatable      :: work(:)
    real(8)                   :: wd(1)
    character(len=300)        :: comment
    logical                   :: converged
    logical                   :: ltmp
    type(connect_object)      :: cobj
    type(dictionary), pointer :: options
    type(dictionary), pointer :: run
    type(mhgps_state)         :: mhgpsst
    type(run_objects)         :: runObj
    type(state_properties)    :: outs
    type(findsad_work)        :: fsw
    type(userinput)           :: uinp

    !simple atomic datastructre
    real(gp), allocatable :: rxyz(:,:),fxyz(:,:)
    real(gp), allocatable :: fat(:,:)
    real(gp), allocatable :: rxyz2(:,:),fxyz2(:,:)
    real(gp), allocatable :: tsguess(:,:),minmodeguess(:,:)
    real(gp), allocatable :: minmode(:,:)
    real(gp), allocatable :: tsgforces(:,:)
    real(gp)              :: tsgenergy
    real(gp), allocatable :: fp(:),fp2(:)
    real(gp), allocatable :: rotforce(:,:),hess(:,:)
    real(gp), allocatable :: eval(:)
    real(gp)              :: energy, energy2, ec, displ
    real(gp)              :: fnoise, fnrm, fmax
    integer :: idum=0
    real(kind=4) :: builtin_rand

    !functions
    real(gp) :: dnrm2

    nbond=1
    converged = .false.
    premature_exit=.false.

    call f_lib_initialize()

    call bigdft_command_line_options(options)
    call bigdft_init(options)!mpi_info,nconfig,run_id,ierr)
    if (bigdft_nruns(options) > 1) then
        call f_err_throw('runs-file not supported for MHGPS '//&
                         'executable')
    endif
    run => options // 'BigDFT' // 0

    !initalize mhgps internal state
    !(only non-system dependent variables)
    call init_mhgps_state(mhgpsst)
    !read user input file mhgps.inp
    call read_input(uinp)
    !obtain first strucutre (used for initialization of
    !bigdft)
    call get_first_struct_file(mhgpsst,filename)

    if(mhgpsst%iproc==0) call print_logo_mhgps(mhgpsst)

    !reset input and output positions of run
    call bigdft_get_run_properties(run,input_id=run_id,&
         naming_id=naming_id)
    call bigdft_set_run_properties(run,&
         posinp_id=trim(adjustl(filename))//trim(naming_id))

    call run_objects_init(runObj,run)

    !now read state of previous mhgps run (if present)
    call read_restart(mhgpsst,runObj)


    !options and run are not needed
    call dict_free(options)
    nullify(run)

    call init_state_properties(outs, bigdft_nat(runObj))


!    elseif(efmethod=='AMBER')then
!        mhgpsst%iproc=0
!        isForceField=.true.
!        call read_atomic_file(trim(adjustl(filename)),&
!             mhgpsst%iproc,atom_struct)
!        astruct_ptr=>atom_struct
!        fdim=astruct_ptr%nat
!        !alanine stuff ......................START!>
!          l_sat=5
!          allocate(atomnamesdmy(1000))
!          rxyzdmy = f_malloc((/ 1.to.3, 1.to.astruct_ptr%nat/),&
!                            id='rxyzdmy')
!          fxyzdmy = f_malloc((/ 1.to.3, 1.to.astruct_ptr%nat/),&
!                            id='fxyzdmy')
!          fnpdb='ald_new.pdb'
!          nfnpdb=len(trim(fnpdb));
!          call nab_init(astruct_ptr%nat,rxyzdmy,fxyzdmy,&
!                        trim(fnpdb),nfnpdb,l_sat,atomnamesdmy)
!          call f_free(rxyzdmy)
!          call f_free(fxyzdmy)
!          deallocate(atomnamesdmy)
!        !alanine stuff ......................END!>
!
!        call print_logo_mhgps()
!    else
!        call yaml_warning('Following method for evaluation of '//&
!        'energies and forces is unknown: '//trim(adjustl(efmethod)))
!        stop
!    endif
    if(mhgpsst%iproc==0) call print_input(uinp)

    mhgpsst%nid = bigdft_nat(runObj) !s-overlap fingerprints
    
    !allocate arrays
    hess     = f_malloc((/ 1.to.3*bigdft_nat(runObj),&
                1.to.3*bigdft_nat(runObj)/),id='hess')
    eval  = f_malloc((/ 1.to.3*bigdft_nat(runObj)/),id='eval')
    call DSYEV('N','L',3*bigdft_nat(runObj),hess,3*bigdft_nat(runObj),eval,wd,-1,info)
    if (info.ne.0) stop 'info query'
    lwork=nint(wd(1))
    work = f_malloc((/ 1.to.lwork/),id='work')
    tsgforces     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='tsgforces')
    tsguess     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='tsguess')
    minmodeguess  = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='minmodeguess')
    minmode  = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='minmode')
    rxyz     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='rxyz')
    fp       = f_malloc((/ 1.to.mhgpsst%nid/),&
                id='fp')
    fp2      = f_malloc((/ 1.to.mhgpsst%nid/),&
                id='fp2')
    fxyz     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='fxyz')
    rxyz2     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='rxyz2')
    fat       = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='fat')
    fxyz2     = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='fxyz2')
    rcov     = f_malloc((/ 1.to.bigdft_nat(runObj)/),id='rcov')
    iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
    rotforce = f_malloc((/ 1.to.3, 1.to.bigdft_nat(runObj)/),&
                id='rotforce')


    call allocate_connect_object(bigdft_nat(runObj),mhgpsst%nid,uinp%nsadmax,cobj)

    iconnect = 0
    call give_rcov(mhgpsst%iproc,bigdft_get_astruct_ptr(runObj),bigdft_nat(runObj),rcov)
    !if in biomode, determine bonds betweens atoms once and for all
    !(it isassuemed that all conifugrations over which will be
    !iterated have the same bonds)
    if(uinp%saddle_biomode)then
        call findbonds('(MHGPS)',mhgpsst%iproc,uinp%mhgps_verbosity,bigdft_nat(runObj),rcov,&
        bigdft_get_rxyz_ptr(runObj),nbond,iconnect)
    endif
    call allocate_finsad_workarrays(runObj,uinp,nbond,fsw)

    ifolderstart=mhgpsst%ifolder
    outer: do ifolder = ifolderstart,999
        mhgpsst%ifolder=ifolder
        write(mhgpsst%currDir,'(a,i3.3)')trim(adjustl(mhgpsst%dirprefix)),ifolder
        call read_jobs(uinp,mhgpsst)
        if(mhgpsst%njobs==0)cycle
        mhgpsst%ijob=0
        if(mhgpsst%iproc==0)then
           call write_restart(mhgpsst,runObj,writeJobList=.false.)
        endif

        inner: do ijob = 1,mhgpsst%njobs
           if(uinp%singlestep .and. (ifolder-ifolderstart)*mhgpsst%njobs+ijob > 1)then
              premature_exit =.true.
              exit outer
           endif
           mhgpsst%ijob=ijob
           call bigdft_get_rxyz(filename=&
                trim(adjustl(mhgpsst%joblist(1,ijob))),rxyz=rxyz)

           select case(trim(adjustl(uinp%operation_mode)))
           case('guessonly')
              call bigdft_get_rxyz(filename=mhgpsst%joblist(2,ijob),&
                   rxyz=rxyz2)

              mhgpsst%isad=mhgpsst%isad+1
              write(mhgpsst%isadc,'(i5.5)')mhgpsst%isad
              !rmsd alignment (optional in mhgps approach)
              call superimpose(bigdft_nat(runObj),rxyz(1,1),rxyz2(1,1))
              runObj%inputs%inputPsiId=0
              call get_ts_guess(mhgpsst,uinp,runObj,outs,rxyz(1,1),&
                   rxyz2(1,1),tsguess(1,1),minmodeguess(1,1),&
                   tsgenergy,tsgforces(1,1))
              write(comment,'(a)')&
                   'TS guess; forces below give guessed '//&
                   'minimummode.'
              if(mhgpsst%iproc==0)then
              call astruct_dump_to_file(&
                   bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//'/sad'//&
                   trim(adjustl(mhgpsst%isadc))//'_ig_finalM',comment,&
                   tsgenergy,rxyz=tsguess,forces=minmodeguess)

              call astruct_dump_to_file(&
                   bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//'/sad'//&
                   trim(adjustl(mhgpsst%isadc))//'_ig_finalF',comment,&
                   tsgenergy,rxyz=tsguess,forces=tsgforces)
              endif
           case('connect')
              call bigdft_get_rxyz(filename=mhgpsst%joblist(2,ijob),&
                   rxyz=rxyz2)

              !Evalute energies. They are needed in connect
              !for identification
              runObj%inputs%inputPsiId=0
              call mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz,&
                                        fat,fnoise,energy,infocode)
              runObj%inputs%inputPsiId=0
              call mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz2,&
                                        fat,fnoise,energy2,infocode)
              call fingerprint(bigdft_nat(runObj),mhgpsst%nid,&
                              runObj%atoms%astruct%cell_dim,&
                              bigdft_get_geocode(runObj),rcov,&
                              rxyz(1,1),fp(1))
              call fingerprint(bigdft_nat(runObj),mhgpsst%nid,&
                              runObj%atoms%astruct%cell_dim,&
                              bigdft_get_geocode(runObj),rcov,&
                              rxyz2(1,1),fp2(1))
              if(mhgpsst%iproc==0)then
                 call yaml_comment('(MHGPS) Connect '//&
                      trim(adjustl(mhgpsst%joblist(1,ijob)))//' and '//&
                      trim(adjustl(mhgpsst%joblist(2,ijob)))//' ....',&
                      hfill='-')
              endif
              isame=0
              connected=.true.
              if(trim(adjustl(mhgpsst%joblist(1,ijob)(10:16)))/='restart')then
                  mhgpsst%nsad=0
              endif
              call connect(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,&
                   iconnect,rxyz,rxyz2,energy,energy2,fp,fp2,&
                   cobj,connected,premature_exit,nsad)
!              call connect_recursively(mhgpsst,fsw,uinp,runObj,outs,rcov,&
!                   nbond,isame,iconnect,rxyz,rxyz2,energy,energy2,fp,&
!                   fp2,cobj,connected)
              if(connected)then
                ltmp=.true.
                if(ijob+1<=mhgpsst%njobs)then
                    !If connected==.true. then in subroutine connect, NO
                    !new connection jobs for a restart have been added.
                    !This meand, the job in mhgpsst%joblist(1,ijob+1) is identical
                    !to the job that will be read an after a restart.
                    if(trim(adjustl(mhgpsst%joblist(1,ijob+1)(10:16)))=='restart')then
                      ltmp=.false.
                    endif
                endif
                
                if(mhgpsst%iproc==0 .and. ltmp)call yaml_map('(MHGPS) '//&
                      'succesfully connected, intermediate'//&
                      ' transition states',nsad)
              else
                 if(.not.premature_exit)then
                 if(mhgpsst%iproc==0)call yaml_comment('(MHGPS) '//&
                      'Connection not established within '//&
                      trim(adjustl(yaml_toa(nsad)))//&
                      ' transition state computations')
                 endif
              endif
              if(premature_exit)then
                 exit outer
              endif
           case('simple')
              mhgpsst%isad=mhgpsst%isad+1
              write(mhgpsst%isadc,'(i3.3)')mhgpsst%isad
              if(uinp%random_minmode_guess)then
                 do iat=1,bigdft_nat(runObj)
                    minmode(1,iat)=2.0_gp*&
                         (real(builtin_rand(idum),gp)-0.5_gp)
                    minmode(2,iat)=2.0_gp*&
                         (real(builtin_rand(idum),gp)-0.5_gp)
                    minmode(3,iat)=2.0_gp*&
                         (real(builtin_rand(idum),gp)-0.5_gp)
                 enddo
                 if(mhgpsst%iproc==0)then
                 call write_mode(runObj,outs,&
                      trim(adjustl(mhgpsst%joblist(1,ijob)))//'_mode',minmode)
                 endif
              else
                 call read_mode(mhgpsst,bigdft_nat(runObj),trim(adjustl(mhgpsst%joblist(1,ijob)))&
                      //'_mode',minmode)
              endif
              !normalize
              minmode = minmode/dnrm2(3*bigdft_nat(runObj),minmode(1,1),1)
              ec=0.0_gp
              displ=0.0_gp
              call findsad(mhgpsst,fsw,uinp,runObj,outs,rcov,nbond,iconnect,&
                   rxyz(1,1),energy,fxyz(1,1),minmode(1,1),displ,ec,&
                   rotforce(1,1),converged)
              if(.not.converged)then
                 call yaml_warning('Saddle '//yaml_toa(mhgpsst%isad)//&
                      ' not converged')
              endif
              call fnrmandforcemax(fxyz(1,1),fnrm,&
                   fmax,bigdft_nat(runObj))
              fnrm = sqrt(fnrm)
              if (mhgpsst%iproc == 0) then
                 write(comment,'(a,1pe10.3,5x1pe10.3)')&
                      'ATTENTION! Forces below give no forces, '//&
                      'but the final minmode| '//&
                      'fnrm, fmax = ',fnrm,fmax

                 call astruct_dump_to_file(&
                      bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                      '/sad'//trim(adjustl(mhgpsst%isadc))//'_finalM',&
                      comment,energy,rxyz=rxyz,forces=minmode)

                 write(comment,'(a,1pe10.3,5x1pe10.3)')&
                      'fnrm, fmax = ',fnrm,fmax
                 call astruct_dump_to_file(&
                      bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//&
                      '/sad'//trim(adjustl(mhgpsst%isadc))//'_finalF',&
                      comment,energy,rxyz=rxyz,forces=fxyz)

                 call write_mode(runObj,outs,mhgpsst%currDir//'/sad'//&
                      trim(adjustl(mhgpsst%isadc))//'_mode_final',&
                      minmode(1,1),rotforce(1,1))
              endif
           case('minimize')
              mhgpsst%isad=mhgpsst%isad+1
              write(mhgpsst%isadc,'(i3.3)')mhgpsst%isad
              ec=0.0_gp
              runObj%inputs%inputPsiId=0
              call mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz,&
                   fxyz,fnoise,energy,infocode)
              call minimize(mhgpsst,uinp,runObj,outs,nbond,&
                   iconnect,rxyz(1,1),fxyz(1,1),fnoise,energy,ec,&
                   converged,'')
              if(.not.converged)then
                 call yaml_warning('Minimization '//yaml_toa(mhgpsst%isad)&
                      //' not converged')
              endif
              call fnrmandforcemax(fxyz(1,1),fnrm,fmax,bigdft_nat(runObj))
              if (mhgpsst%iproc == 0) then
                 write(comment,'(a,1pe10.3,5x1pe10.3)')&
                      'fnrm, fmax = ',fnrm,fmax

                 call astruct_dump_to_file(&
                      bigdft_get_astruct_ptr(runObj),mhgpsst%currDir//'/min'&
                      //trim(adjustl(mhgpsst%isadc))//'_final',comment,&
                      energy,rxyz=rxyz,forces=fxyz)
              endif
           case('hessian')
              runObj%inputs%inputPsiId=0
              call mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz,&
                   fxyz,fnoise,energy,infocode)
              runObj%inputs%inputPsiId=0
              call cal_hessian_fd(mhgpsst,runObj,outs,rxyz,&
                   hess,nfree)
              if(mhgpsst%iproc==0)then
                 write(*,*)'(hess) HESSIAN:'
                 write(*,*)hess
              endif
!              call DSYEV('V','L',3*bigdft_nat(runObj),hess,3*bigdft_nat(runObj),eval,WORK,LWORK,&
!                   INFO)
              call DSYEV('V','L',nfree,hess,3*bigdft_nat(runObj),eval,WORK,LWORK,&
                   INFO)
              if (info.ne.0) stop 'DSYEV'
              if(mhgpsst%iproc==0)then
                 write(*,*)'(hess) EIGENVECTORS:'
                 write(*,*) hess
                 write(*,'(a,1x,es9.2,1x,es24.17)') '(hess)'//&
                      ' ---   App. eigenvalues in exact --------'//&
                      '--- fnrm:',sqrt(sum(fxyz**2)),energy
                 do iat=1,nfree
                    write(*,*) '(hess) eval ',iat,eval(iat)
                 enddo
              endif
           case default
              call yaml_warning('(MHGPS) operation mode '//&
                   trim(adjustl(uinp%operation_mode))//' unknown STOP')
              stop '(MHGPS) operation mode unknown STOP'
           end select
        enddo inner
        mhgpsst%isad=0
        if(mhgpsst%iproc==0)then
!        call f_delete_file('restart')
        call f_delete_file(trim(adjustl(mhgpsst%currDir))//'/job_list_restart')
        endif
     enddo outer

    if(mhgpsst%iproc==0 .and. (.not. premature_exit))then
       call f_delete_file('restart')
        u=f_get_free_unit()
        open(unit=u,file='finished') 
        write(u,*)'finished'
        close(u)
    endif

    !finalize (dealloctaion etc...)
    call free_run_objects(runObj)
    call deallocate_state_properties(outs)
    call bigdft_finalize(ierr)

    call f_free(work)
    call f_free(eval)
    call f_free(tsguess)
    call f_free(tsgforces)
    call f_free(minmodeguess)
    call f_free(minmode)
    call f_free(fp)
    call f_free(fp2)
    call f_free(rxyz)
    call f_free(fat)
    call f_free(fxyz) 
    call f_free(rxyz2)
    call f_free(fxyz2) 
    call f_free(rcov)
    call f_free(iconnect)
    call f_free(rotforce)
    call f_free(hess)
    call deallocate_connect_object(cobj)
    call deallocate_finsad_workarrays(fsw)

    call finalize_mhgps_state(mhgpsst)

    if(mhgpsst%iproc==0)call yaml_map('(MHGPS) Total calls to energy and '//&
                               'forces',nint(mhgpsst%ef_counter))
    if(mhgpsst%iproc==0)call yaml_map('(MHGPS) Run finished at',&
                               yaml_date_and_time_toa())
    call f_lib_finalize()
end program

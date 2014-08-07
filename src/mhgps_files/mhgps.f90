!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS
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
                            read_atomic_file=>set_astruct_from_file
    use module_global_variables
    use module_init
    use module_energyandforces
    use module_sbfgs !for finding binds
    use module_freezingstring
    use module_saddle
    use module_connect
    use module_fingerprints
    implicit none
    integer :: nstart
    character(len=6) :: filename
    integer :: ifolder, ifile
    logical :: xyzexists,asciiexists
    character(len=60) :: run_id
    integer :: ierr, nconfig
    real(gp), allocatable :: rcov(:)
    character(len=300)  :: comment
    logical :: converged=.false.
    type(connect_object) :: cobj
    integer :: nsad
    logical :: connected

    !simple atomic datastructre
    integer :: nat
    integer :: nid
    real(gp) :: alat(3)
    real(gp), allocatable :: rxyz(:,:),fxyz(:,:)
    real(gp), allocatable :: rxyz2(:,:),fxyz2(:,:)
    real(gp), allocatable :: tsguess(:,:),minmodeguess(:,:)
    real(gp), allocatable :: tsgforces(:,:)
    real(gp)              :: tsgenergy
    real(gp), allocatable :: fp(:),fp2(:)
    real(gp), allocatable :: rotforce(:,:),hess(:,:)
    real(gp), allocatable :: eval(:)
    real(gp) :: energy, energy2, ec, displ
    real(gp) :: fnoise, fnrm, fmax
    integer :: i,j,info
    integer :: idum=0
    real(kind=4) :: builtin_rand

    !functions
    real(gp) :: dnrm2

    !alanine stuff ......................START!>
    real(gp), allocatable :: rxyzdmy(:,:), fxyzdmy(:,:)
    character(len=5), allocatable :: atomnamesdmy(:)
    integer :: l_sat, nfnpdb
    character(len=11) :: fnpdb
    !alanine stuff ......................END!>

!for debugging
real(gp) :: etest
real(gp), allocatable :: fat(:,:)
!debugging end

    ifolder=1
    ifile=1
    ef_counter=0.d0 !from module_global_variables
    isad=0  !from module_global_variables


    call f_lib_initialize()

    !read mhgps.inp
    call read_input()

    !initialize the energy and forces method
    isForceField=.false.
    if(efmethod=='BIGDFT')then
        isForceField=.false.
        call bigdft_init(mpi_info,nconfig,run_id,ierr)
        iproc=mpi_info(1)
        nproc=mpi_info(2)
        igroup=mpi_info(3)
        !number of groups
         ngroups=mpi_info(4)
        !actual value of iproc
        iproc=iproc+igroup*ngroups
        if (nconfig < 0) then
            call yaml_warning('runs-file not supported for MHGPS&
                               executable')
            stop
        endif
        if(iproc==0) call print_logo_mhgps()
        call dict_init(user_inputs)
        write(currDir,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call user_dict_from_files(user_inputs, trim(run_id)//&
             trim(bigdft_run_id_toa()),currDir//'/'//filename//&
             trim(bigdft_run_id_toa()), bigdft_mpi)
        call inputs_from_dict(inputs_opt, atoms, user_inputs)
        call dict_free(user_inputs)
        call init_global_output(outs, atoms%astruct%nat)
        fdim=outs%fdim
        call init_restart_objects(bigdft_mpi%iproc,inputs_opt,atoms,&
                                 rst)
        call run_objects_nullify(runObj)
        call run_objects_associate(runObj, inputs_opt, atoms, rst)
        !set minimum number of wave function optimizations
!        if(runObj%inputs%itermin<5)then
!            itermin=5
!        else
            itermin=runObj%inputs%itermin
!        endif

    elseif(efmethod=='LJ')then
        iproc=0
        isForceField=.true.
        write(currDir,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call deallocate_atomic_structure(atoms%astruct)
        call read_atomic_file(currDir//'/'//filename,iproc,&
                              atoms%astruct)
        fdim=atoms%astruct%nat
        call print_logo_mhgps()
    elseif(efmethod=='AMBER')then
        iproc=0
        isForceField=.true.
        write(currDir,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call deallocate_atomic_structure(atoms%astruct)
        call read_atomic_file(currDir//'/'//filename,iproc,&
                              atoms%astruct)
        fdim=atoms%astruct%nat
        !alanine stuff ......................START!>
          l_sat=5
          allocate(atomnamesdmy(1000))
          rxyzdmy = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                            id='rxyzdmy')
          fxyzdmy = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                            id='fxyzdmy')
          fnpdb='ald_new.pdb'
          nfnpdb=len(trim(fnpdb));
          call nab_init(atoms%astruct%nat,rxyzdmy,fxyzdmy,&
                        trim(fnpdb),nfnpdb,l_sat,atomnamesdmy)
          call f_free(rxyzdmy)
          call f_free(fxyzdmy)
          deallocate(atomnamesdmy)
        !alanine stuff ......................END!>

        call print_logo_mhgps()
    else
        call yaml_warning('Following method for evaluation of &
        energies and forces is unknown: '//trim(adjustl(efmethod)))
        stop
    endif
    if(iproc==0) call print_input()

    nat = atoms%astruct%nat
    nid = nat !s-overlap fingerprints
    alat = atoms%astruct%cell_dim
    
    !allocate more arrays
    lwork=1000+10*nat**2
    work = f_malloc((/1.to.lwork/),id='work')
    eval  = f_malloc((/ 1.to.3*nat/),&
                id='eval')
    tsgforces     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='tsgforces')
    tsguess     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='tsguess')
    minmodeguess  = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='minmodeguess')
    minmode  = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='minmode')
    rxyz     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='rxyz')
    fp       = f_malloc((/ 1.to.nid/),&
                id='fp')
    fp2      = f_malloc((/ 1.to.nid/),&
                id='fp2')
    fxyz     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='fxyz')
    rxyz2     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='rxyz2')
    fxyz2     = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='fxyz2')
    rcov     = f_malloc((/ 1.to.nat/),id='rcov')
    iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
    ixyz_int = f_malloc((/ 1.to.3,1.to.nat/),&
                id='nat')
    rotforce = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='rotforce')
    hess     = f_malloc((/ 1.to.3*nat,&
                1.to.3*nat/),id='hess')
    rxyz_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='rxyz_rot')
    fxyz_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='fxyz_rot') 
    fxyzraw_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='fxyzraw_rot')
    rxyzraw_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='rxyzraw_rot')
    fstretch_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='fstretch_rot')
    eval_rot = f_malloc((/1.to.saddle_nhistx_rot/),id='eval_rot')
    res_rot = f_malloc((/1.to.saddle_nhistx_rot/),id='res_rot')
    rrr_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='rrr_rot')
    aa_rot = f_malloc((/1.to.saddle_nhistx_rot,&
             1.to.saddle_nhistx_rot/),id='aa_rot')
    ff_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='ff_rot')
    rr_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='rr_rot')
    dd_rot = f_malloc((/ 1.to.3, 1.to.nat/),&
                id='dd_rot')
    fff_rot = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_rot/),id='fff_rot')
    scpr_rot = f_malloc((/ 1.to.saddle_nhistx_rot/),id='scpr_rot')
    rxyz_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='rxyz_trans')
    fxyz_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='fxyz_trans') 
    fxyzraw_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='fxyzraw_trans')
    rxyzraw_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='rxyzraw_trans')
    fstretch_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='fstretch_trans')
    eval_trans = f_malloc((/1.to.saddle_nhistx_trans/),&
                 id='eval_trans')
    res_trans = f_malloc((/1.to.saddle_nhistx_trans/),id='res_trans')
    rrr_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='rrr_trans')
    aa_trans = f_malloc((/1.to.saddle_nhistx_trans,&
                1.to.saddle_nhistx_trans/),id='aa_trans')
    ff_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='ff_trans')
    rr_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                0.to.saddle_nhistx_trans/),id='rr_trans')
    dd_trans = f_malloc((/ 1.to.3, 1.to.nat/),id='dd_trans')
    fff_trans = f_malloc((/ 1.to.3, 1.to.nat,&
                -1.to.saddle_nhistx_trans/),id='fff_trans')
    scpr_trans = f_malloc((/ 1.to.saddle_nhistx_trans/),&
                 id='scpr_trans')

    call allocate_connect_object(nat,nid,nsadmax,cobj)

!for debugging:
allocate(fat(3,nat))
    
    iconnect = 0
    ixyz_int = 0
    call give_rcov(atoms,nat,rcov)
    !if in biomode, determine bonds betweens atoms once and for all
    !(it isassuemed that all conifugrations over which will be
    !iterated have the same bonds)
    if(saddle_biomode)then
        call findbonds('(MHGPS)',iproc,mhgps_verbosity,nat,rcov,&
        atoms%astruct%rxyz,nbond,iconnect)
    endif
    wold_trans = f_malloc((/ 1.to.nbond/),id='wold_trans')
    wold_rot = f_malloc((/ 1.to.nbond/),id='wold_rot')

    do ifolder = 1,999
        write(currDir,'(a,i3.3)')'input',ifolder

        if(trim(adjustl(operation_mode))=='simple'.or.&
                       trim(adjustl(operation_mode))=='hessian')then
            nstart=1
        elseif(trim(adjustl(operation_mode))=='connect'.or.&
                       trim(adjustl(operation_mode))=='guessonly')then
            nstart=2
        else
            call yaml_warning('(MHGPS) operation mode unknown STOP')
            stop '(MHGPS) operation mode unknown STOP'
        endif

        do ifile = nstart,999

            !read (first) file
            write(filename,'(a,i3.3)')'pos',ifile
            inquire(file=currDir//'/'//filename//'.xyz',&
                    exist=xyzexists)
            inquire(file=currDir//'/'//filename//'.ascii',&
                    exist=asciiexists)
            if(.not.(xyzexists.or.asciiexists))exit
            call deallocate_atomic_structure(atoms%astruct)
            call read_atomic_file(currDir//'/'//filename,iproc,&
                    atoms%astruct)
            call vcopy(3 * nat,atoms%astruct%rxyz(1,1),1,rxyz(1,1), 1)

            if(trim(adjustl(operation_mode))=='guessonly')then
                !read second file
                write(filename,'(a,i3.3)')'pos',ifile+1
                inquire(file=currDir//'/'//filename//'.xyz',&
                            exist=xyzexists)
                inquire(file=currDir//'/'//filename//'.ascii',&
                            exist=asciiexists)
                if(.not.(xyzexists.or.asciiexists))exit
                call deallocate_atomic_structure(atoms%astruct)
                call read_atomic_file(currDir//'/'//filename,iproc,&
                            atoms%astruct)
                call vcopy(3*nat,atoms%astruct%rxyz(1,1),1,&
                           rxyz2(1,1),1)

                isad=isad+1
                write(isadc,'(i5.5)')isad
                call get_ts_guess(nat,alat,rxyz(1,1),rxyz2(1,1),&
                     tsguess(1,1),minmodeguess(1,1),tsgenergy,&
                     tsgforces(1,1))
                write(comment,'(a)')&
                     'TS guess; forces below give guessed &
                      minimummode.'
                call write_atomic_file(currDir//'/sad'//&
                     trim(adjustl(isadc))//'_ig_finalM',&
                     tsgenergy,tsguess(1,1),ixyz_int,atoms,&
                     comment,forces=minmodeguess(1,1))
                call write_atomic_file(currDir//'/sad'//&
                     trim(adjustl(isadc))//'_ig_finalF',&
                     tsgenergy,tsguess(1,1),ixyz_int,atoms,&
                     comment,forces=tsgforces(1,1))
            else if(trim(adjustl(operation_mode))=='connect')&
                                                                 then
                !read second file
                write(filename,'(a,i3.3)')'pos',ifile+1
                inquire(file=currDir//'/'//filename//'.xyz',&
                            exist=xyzexists)
                inquire(file=currDir//'/'//filename//'.ascii',&
                            exist=asciiexists)
                if(.not.(xyzexists.or.asciiexists))exit
                call deallocate_atomic_structure(atoms%astruct)
                call read_atomic_file(currDir//'/'//filename,iproc,&
                            atoms%astruct)
                call vcopy(3*nat,atoms%astruct%rxyz(1,1),1,&
                           rxyz2(1,1),1)

                !Evalute energyies. They are needed in connect
                !for identification
                call energyandforces(nat,alat,rxyz,fat,fnoise,energy)
                call energyandforces(nat,alat,rxyz2,fat,fnoise,&
                         energy2)
                call fingerprint(nat,nid,alat,atoms%astruct%geocode,&
                         rcov,rxyz(1,1),fp(1))
                call fingerprint(nat,nid,alat,atoms%astruct%geocode,&
                         rcov,rxyz2(1,1),fp2(1))
                nsad=0
                connected=.true.
                call connect_recursively(nat,nid,alat,rcov,nbond,&
                     iconnect,rxyz,rxyz2,energy,energy2,fp,fp2,&
                     nsad,cobj,connected)
                if(connected)then
                    if(iproc==0)call yaml_map('(MHGPS) &
                                succesfully connected',nsad)
                else
                    if(iproc==0)call yaml_map('(MHGPS) &
                                Connection not established',nsad)
                endif
            else if(trim(adjustl(operation_mode))=='simple')then
                isad=isad+1
                write(isadc,'(i3.3)')isad
                if(random_minmode_guess)then
                    do i=1,nat
                        minmode(1,i)=2.0_gp*&
                                (real(builtin_rand(idum),gp)-0.5_gp)
                        minmode(2,i)=2.0_gp*&
                                (real(builtin_rand(idum),gp)-0.5_gp)
                        minmode(3,i)=2.0_gp*&
                                (real(builtin_rand(idum),gp)-0.5_gp)
                    enddo
                    call write_mode(nat,currDir//'/pos'//&
                    trim(adjustl(isadc))//'_mode',minmode)
                else
                    call read_mode(nat,currDir//'/pos'//&
                    trim(adjustl(isadc))//'_mode',minmode)
                endif
                !!call random_seed
                !!call random_number(minmode)
                !!minmode=2.d0*(minmode-0.5d0)
                !normalize
                minmode = minmode/dnrm2(3*nat,minmode(1,1),1)
                ec=0.0_gp
                displ=0.0_gp
                call findsad(nat,alat,rcov,nbond,iconnect,&
                rxyz(1,1),energy,fxyz(1,1),minmode(1,1),displ,ec,&
                rotforce(1,1),converged)
!for debugging:
call energyandforces(nat,alat,rxyz,fat,fnoise,etest)
write(*,*)'energy MAINcheck saddle: ',energy-etest
                    if(.not.converged)then
                        call yaml_warning('Saddle '//yaml_toa(isad)//&
                            ' not converged')
                    endif
                    call fnrmandforcemax(fxyz(1,1),fnrm,&
                                        fmax,nat)
                    if (iproc == 0) then
                        write(comment,'(a,1pe10.3,5x1pe10.3)')&
                       'ATTENTION! Forces below give no forces, &
                        but the final minmode| &
                        fnrm, fmax = ',fnrm,fmax

                        call write_atomic_file(currDir//'/sad'//&
                        trim(adjustl(isadc))//'_finalM',energy,&
                        rxyz(1,1),ixyz_int,atoms,comment,&
                        forces=minmode(1,1))
                        !atoms,comment,forces=fxyz(1,1))

                        write(comment,'(a,1pe10.3,5x1pe10.3)')&
                       'fnrm, fmax = ',fnrm,fmax
                        call write_atomic_file(currDir//'/sad'//&
                        trim(adjustl(isadc))//'_finalF',energy,&
                        rxyz(1,1),ixyz_int,atoms,comment,&
                        forces=fxyz(1,1))
                        
                        call write_mode(nat,currDir//'/sad'//&
                        trim(adjustl(isadc))//'_mode_final',&
                        minmode(1,1),rotforce(1,1))
                    endif
            else if(trim(adjustl(operation_mode))=='hessian')then
                call cal_hessian_fd(iproc,nat,alat,rxyz,hess)
                call DSYEV('V','L',3*nat,hess,3*nat,eval,WORK,LWORK,&
                     INFO)
                if (info.ne.0) stop 'DSYEV'
                    if(iproc==0)then
                        write(*,'(a,1x,es9.2,1x,es24.17)') '(hess) &
                        ---   App. eigenvalues in exact ------------&
                        - fnrm:',sqrt(sum(fxyz**2)),energy
                        do j=1,3*nat
                            write(*,*) '(hess) eval ',j,eval(j)
                        enddo
                    endif
            else
                call yaml_warning('(MHGPS) operation mode unknown &
                                  STOP')
                stop '(MHGPS) operation mode unknown STOP'
            endif
        enddo
    enddo

    !finalize (dealloctaion etc...)
    if(efmethod=='BIGDFT')then
        call free_restart_objects(rst)
        call deallocate_atoms_data(atoms)
        call deallocate_global_output(outs)
        call run_objects_free_container(runObj)
        call free_input_variables(inputs_opt)
        call bigdft_finalize(ierr)
    elseif(efmethod=='LJ'.or.efmethod=='AMBER')then
        call deallocate_atoms_data(atoms)
    endif

    call f_free(work)
    call f_free(eval)
    call f_free(tsguess)
    call f_free(tsgforces)
    call f_free(minmodeguess)
    call f_free(minmode)
    call f_free(fp)
    call f_free(fp2)
    call f_free(rxyz)
    call f_free(fxyz) 
    call f_free(rxyz2)
    call f_free(fxyz2) 
    call f_free(rcov)
    call f_free(iconnect)
    call f_free(ixyz_int)
    call f_free(rotforce)
    call f_free(hess)
    call f_free(rxyz_rot)
    call f_free(fxyz_rot)
    call f_free(fxyzraw_rot)
    call f_free(rxyzraw_rot)
    call f_free(fstretch_rot)
    call f_free(eval_rot)
    call f_free(res_rot)
    call f_free(rrr_rot)
    call f_free(aa_rot)
    call f_free(ff_rot)
    call f_free(rr_rot)
    call f_free(dd_rot)
    call f_free(fff_rot)
    call f_free(scpr_rot)
    call f_free(rxyz_trans)
    call f_free(fxyz_trans)
    call f_free(fxyzraw_trans)
    call f_free(rxyzraw_trans)
    call f_free(fstretch_trans)
    call f_free(eval_trans)
    call f_free(res_trans)
    call f_free(rrr_trans)
    call f_free(aa_trans)
    call f_free(ff_trans)
    call f_free(rr_trans)
    call f_free(dd_trans)
    call f_free(fff_trans)
    call f_free(scpr_trans)
    call f_free(wold_trans)
    call f_free(wold_rot)
    call deallocate_connect_object(cobj)

    


    if(iproc==0)call yaml_map('(MHGPS) Total calls to energy and forces',nint(ef_counter))
    if(iproc==0)call yaml_map('(MHGPS) Run finished at',yaml_date_and_time_toa())
    call f_lib_finalize()
end program

subroutine cal_hessian_fd(iproc,nat,alat,pos,hess)
use module_base
use module_energyandforces
    implicit none
    integer, intent(in):: iproc, nat
    real(8), intent(in) :: pos(3*nat),alat(3)
    real(8), intent(inout) :: hess(3*nat,3*nat)
    !local variables
    integer :: iat
    real(8) :: t1,t2,t3
    !real(8), allocatable, dimension(:,:) :: hess
    real(8), allocatable, dimension(:) :: tpos,grad,eval,work
    real(8) :: h,rlarge,twelfth,twothird,etot,cmx,cmy,cmz,shift,dm,tt,s,fnoise
    integer :: i,j,k,lwork,info

    !allocate(hess(3*nat,3*nat))
    tpos = f_malloc((/1.to.3*nat/),id='tpos')
    grad = f_malloc((/1.to.3*nat/),id='grad')
    eval = f_malloc((/1.to.3*nat/),id='eval')

    lwork=1000*nat
    work = f_malloc((/1.to.lwork/),id='work')

    !h=1.d-1
    !h=7.5d-2
    !h=5.d-2
!    h=1.d-3
!    h=1.d-2
    h=5.d-3
    !h=2.d-2
    rlarge=1.d0*1.d4
    twelfth=-1.d0/(12.d0*h)
    twothird=-2.d0/(3.d0*h)
    if(iproc==0) write(*,*) '(hess) HESSIAN: h',h
    !-------------------------------------------------------
    do i=1,3*nat
        iat=(i-1)/3+1
        do k=1,3*nat
            tpos(k)=pos(k)
            grad(k)=0.d0
        enddo
        !-----------------------------------------
        tpos(i)=tpos(i)-2*h
        call energyandforces(nat,alat,tpos,grad,fnoise,etot)
        do j=1,3*nat
            hess(j,i)=twelfth*grad(j)
        enddo
        !if(iproc==0) write(*,*) 'ALIREZA-6',i,iat
        !-----------------------------------------
        tpos(i)=tpos(i)+h
        call energyandforces(nat,alat,tpos,grad,fnoise,etot)
        do j=1,3*nat
        hess(j,i)=hess(j,i)-twothird*grad(j)
        enddo
        !-----------------------------------------
        tpos(i)=tpos(i)+2*h
        call energyandforces(nat,alat,tpos,grad,fnoise,etot)
        do j=1,3*nat
        hess(j,i)=hess(j,i)+twothird*grad(j)
        enddo
        !-----------------------------------------
        tpos(i)=tpos(i)+h
        call energyandforces(nat,alat,tpos,grad,fnoise,etot)
        do j=1,3*nat
        hess(j,i)=hess(j,i)-twelfth*grad(j)
        !write(*,*) 'HESS ',j,i,hess(j,i)
        enddo
        !-----------------------------------------
    enddo
    !-------------------------------------------------------

    !check symmetry
    dm=0.d0
    do i=1,3*nat
    do j=1,i-1
    s=.5d0*(hess(i,j)+hess(j,i))
    tt=abs(hess(i,j)-hess(j,i))/(1.d0+abs(s))
    dm=max(dm,tt)
    hess(i,j)=s
    hess(j,i)=s
    enddo
    enddo
    if (dm.gt.1.d-1) write(*,*) '(hess) max dev from sym',dm

!    do j=1,3*nat
!    do i=1,3*nat
!    write(*,*) '(hess) hier',nat,hess(i,j)
!    write(499,*) hess(i,j)
!    enddo
!    enddo

    !-------------------------------------------------------
    !project out rotations
    cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
    do i=1,3*nat-2,3
    cmx=cmx+pos(i+0)
    cmy=cmy+pos(i+1)
    cmz=cmz+pos(i+2)
    enddo
    cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat
  
    !x-y plane
    do i=1,3*nat-2,3
    work(i+1)= (pos(i+0)-cmx)
    work(i+0)=-(pos(i+1)-cmy)
    enddo
    call f_free(tpos)
    call f_free(grad)
    call f_free(eval)
    call f_free(work)
end subroutine cal_hessian_fd

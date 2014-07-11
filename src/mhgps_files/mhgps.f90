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
    use module_atoms, only: deallocate_atoms_data,&
                            deallocate_atomic_structure,&
                            atomic_structure,&
                            read_atomic_file=>set_astruct_from_file
    use module_global_variables
    use module_init
    use module_energyandforces
    use module_sbfgs
    use module_saddle
    use module_minimizers
    use module_io
    implicit none
    integer :: bigdft_get_number_of_atoms,bigdft_get_number_of_orbitals
    character(len=8) :: folder
    character(len=6) :: filename
    character(len=6) :: comment='saddle'
    integer :: ifolder, ifile
    logical :: xyzexists,asciiexists

    character(len=60) :: run_id
    integer :: ierr, nconfig

    real(gp) :: displ,ec
    logical :: converged
    real(gp), allocatable :: rcov(:)

integer :: ncount_bigdft
logical :: fail

real(gp) :: curv

!simple atomic datastructre
integer :: nat
real(gp),allocatable :: rxyz(:,:),fxyz(:,:),rotforce(:,:),hess(:,:)
real(gp) :: energy
integer :: i,j,info
integer :: idum=0
real(kind=4) :: tt,builtin_rand
real(gp),allocatable :: eval(:)
        !alanine stuff ......................START!>
        real(gp), allocatable :: rxyzdmy(:,:), fxyzdmy(:,:)
        character(len=5), allocatable :: atomnamesdmy(:)
        integer :: l_sat, nfnpdb
        character(len=11) :: fnpdb
        !alanine stuff ......................END!>

    ifolder=1
    ifile=1
    ef_counter=0.d0



    call f_lib_initialize()

    call read_input()
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
            call yaml_warning('runs-file not supported for MHGPS executable')
            stop
        endif
        if(iproc==0) call print_logo_mhgps()
        call dict_init(user_inputs)
        write(folder,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call user_dict_from_files(user_inputs, trim(run_id)//trim(bigdft_run_id_toa()), &
           & folder//'/'//filename//trim(bigdft_run_id_toa()), bigdft_mpi)
        call inputs_from_dict(inputs_opt, atoms, user_inputs)
        call dict_free(user_inputs)
        call init_global_output(outs, atoms%astruct%nat)
        call init_restart_objects(bigdft_mpi%iproc,inputs_opt,atoms,rst)
        call run_objects_nullify(runObj)
        call run_objects_associate(runObj, inputs_opt, atoms, rst)
!        if(runObj%inputs%itermin<5)then
!            itermin=5
!        else
            itermin=runObj%inputs%itermin
!        endif

    elseif(efmethod=='LJ')then
        iproc=0
        isForceField=.true.
        write(folder,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call deallocate_atomic_structure(atoms%astruct)
        call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
        call init_global_output(outs, atoms%astruct%nat)
        call print_logo_mhgps()
    elseif(efmethod=='AMBER')then
        iproc=0
        isForceField=.true.
        write(folder,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'pos',ifile
        call deallocate_atomic_structure(atoms%astruct)
        call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
        call init_global_output(outs, atoms%astruct%nat)
        !alanine stuff ......................START!>
          l_sat=5
          allocate(rxyzdmy(3,1000),fxyzdmy(3,1000),atomnamesdmy(1000))
          fnpdb='ald_new.pdb'
          nfnpdb=len(trim(fnpdb));
          call nab_init(atoms%astruct%nat,rxyzdmy,fxyzdmy,trim(fnpdb),nfnpdb,l_sat,atomnamesdmy)
          deallocate(rxyzdmy,fxyzdmy,atomnamesdmy)
        !alanine stuff ......................END!>

        call print_logo_mhgps()
    else
        call yaml_warning('Following method for evaluation of energies and forces is unknown: '//trim(adjustl(efmethod)))
        stop
    endif
    if(iproc==0) call print_input()

    !allocate more arrays
    lwork=1000+10*atoms%astruct%nat**2
    work = f_malloc((/1.to.lwork/),id='work')
    minmode  = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                id='minmode')
    rxyz     = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                id='rxyz')
    fxyz     = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                id='fxyz')
    rcov     = f_malloc((/ 1.to.atoms%astruct%nat/),id='rcov')
    iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
    ixyz_int = f_malloc((/ 1.to.3,1.to.atoms%astruct%nat/),&
                id='atoms%astruct%nat')
    rotforce = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                id='rotforce')
    hess     = f_malloc((/ 1.to.3*atoms%astruct%nat,&
                1.to.3*atoms%astruct%nat/),id='hess')
    rxyz_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='rxyz_rot')
    fxyz_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='fxyz_rot') 
    fxyzraw_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='fxyzraw_rot')
    rxyzraw_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='rxyzraw_rot')
    fstretch_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='fstretch_rot')
    eval_rot = f_malloc((/1.to.saddle_nhistx_rot/),id='eval_rot')
    res_rot = f_malloc((/1.to.saddle_nhistx_rot/),id='res_rot')
    rrr_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='rrr_rot')
    aa_rot = f_malloc((/1.to.saddle_nhistx_rot,1.to.saddle_nhistx_rot/),&
                id='aa_rot')
    ff_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='ff_rot')
    rr_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='rr_rot')
    dd_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),&
                id='dd_rot')
    fff_rot = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_rot/),id='fff_rot')
    scpr_rot = f_malloc((/ 1.to.saddle_nhistx_rot/),id='scpr_rot')
    rxyz_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='rxyz_trans')
    fxyz_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='fxyz_trans') 
    fxyzraw_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='fxyzraw_trans')
    rxyzraw_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='rxyzraw_trans')
    fstretch_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='fstretch_trans')
    eval_trans = f_malloc((/1.to.saddle_nhistx_trans/),id='eval_trans')
    res_trans = f_malloc((/1.to.saddle_nhistx_trans/),id='res_trans')
    rrr_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='rrr_trans')
    aa_trans = f_malloc((/1.to.saddle_nhistx_trans,&
                1.to.saddle_nhistx_trans/),id='aa_trans')
    ff_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='ff_trans')
    rr_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                0.to.saddle_nhistx_trans/),id='rr_trans')
    dd_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),id='dd_trans')
    fff_trans = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat,&
                -1.to.saddle_nhistx_trans/),id='fff_trans')
    scpr_trans = f_malloc((/ 1.to.saddle_nhistx_trans/),id='scpr_trans')
    
    iconnect = 0
    ixyz_int = 0
    call give_rcov(atoms,atoms%astruct%nat,rcov)
    !if in biomode, determine bonds betweens atoms once and for all (it is
    !assuemed that all conifugrations over which will be iterated have the same
    !bonds)
    if(saddle_biomode)then
        call findbonds(atoms%astruct%nat,rcov,atoms%astruct%rxyz,nbond,iconnect)
    endif
    wold_trans = f_malloc((/ 1.to.nbond/),id='wold_trans')
    wold_rot = f_malloc((/ 1.to.nbond/),id='wold_rot')
!allocate(rxyz_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot),fxyz_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot),fxyzraw_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot),rxyzraw_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot),fstretch_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot),eval_rot(saddle_nhistx_rot),res_rot(saddle_nhistx_rot),rrr_rot(3,atoms%astruct%nat,0:saddle_nhistx_rot))





!    if (iproc == 0) &
!        call yaml_set_stream(unit=usaddle,filename=trim(saddle_filename),tabbing=0,record_length=100,setdefault=.false.,istat=ierr)
!    if(ierr==0)then
!        call yaml_warning('Error while opening saddle.mon file. STOP.')
!        stop
!    endif
!    if (iproc ==0 ) call yaml_comment('Saddle monitoring file opened, name:'//trim(filename)//', timestamp: '//trim(yaml_date_and_time_toa()),&
!       hfill='-',unit=usaddle)


!   LWORK=3*3*atoms%astruct%nat-1
!   allocate(eval(3*atoms%astruct%nat))

    do ifolder = 1,999
        write(folder,'(a,i3.3)')'input',ifolder
        currDir=folder
        do ifile = 1,999
            write(filename,'(a,i3.3)')'pos',ifile
            inquire(file=folder//'/'//filename//'.xyz',exist=xyzexists)
            inquire(file=folder//'/'//filename//'.ascii',exist=asciiexists)
            currFile=filename
            if(.not.(xyzexists.or.asciiexists))exit
            call deallocate_atomic_structure(atoms%astruct)
            call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
            call vcopy(3 * atoms%astruct%nat,atoms%astruct%rxyz(1,1),1,rxyz(1,1), 1)
            call vcopy(3 * atoms%astruct%nat,outs%fxyz(1,1),1,fxyz(1,1), 1)
!!            call energyandforces(atoms%astruct%nat,atoms%astruct%cell_dim,rxyz,fxyz,energy)
!!call cal_hessian_fd(iproc,atoms%astruct%nat,atoms%astruct%cell_dim,rxyz,hess)
!!        call DSYEV('V','L',3*atoms%astruct%nat,hess,3*atoms%astruct%nat,eval,WORK,LWORK,INFO)
!!        if (info.ne.0) stop 'DSYEV'
!!        write(*,*) '---   App. eigenvalues in exact -------------'
!!        do j=1,10
!!            write(*,*) 'eval ',j,eval(j)
!!        enddo
            rotforce=0.0_gp
            if(random_minmode_guess)then
                do i=1,atoms%astruct%nat
                    minmode(1,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
                    minmode(2,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
                    minmode(3,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
                enddo
                call write_mode(atoms%astruct%nat,currDir//'/'//currFile//'_mode',minmode)
            else
                call read_mode(atoms%astruct%nat,currDir//'/'//currFile//'_mode',minmode)
            endif
            ec=0.0_gp
        
           call findsad(saddle_imode,atoms%astruct%nat,atoms%astruct%cell_dim,rcov,saddle_alpha0_trans,saddle_alpha0_rot,saddle_curvgraddiff,saddle_nit_trans,&
           saddle_nit_rot,saddle_nhistx_trans,saddle_nhistx_rot,saddle_tolc,saddle_tolf,saddle_tightenfac,saddle_rmsdispl0,&
           saddle_trustr,rxyz,energy,fxyz,minmode,saddle_fnrmtol,displ,ec,&
           converged,atoms%astruct%atomnames,nbond,iconnect,saddle_alpha_stretch0,saddle_recompIfCurvPos,saddle_maxcurvrise,&
           saddle_cutoffratio,saddle_alpha_rot_stretch0,rotforce)
           if (iproc == 0) then
               call write_atomic_file(currDir//'/'//currFile//'_final',&
               energy,rxyz(1,1),ixyz_int,&
               atoms,comment,forces=fxyz(1,1))
               call write_mode(atoms%astruct%nat,currDir//'/'//currFile//'_mode_final',minmode,rotforce)
           endif

!call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
!call minimizer_sbfgs(runObj,outs,nproc,iproc,1,ncount_bigdft,fail)
!rxyz=atoms%astruct%rxyz
!fxyz=outs%fxyz
!call curvgrad(atoms%astruct%nat,atoms%astruct%cell_dim,1.d-3,rxyz,fxyz,minmode,curv,rotforce,1,ec)
!rxyz=atoms%astruct%rxyz
!fxyz=outs%fxyz
        enddo
    enddo

!    !compute minmode only:
!    do ifolder = 1,999
!        do ifile = 1,999
!            write(folder,'(a,i3.3)')'input',ifolder
!            write(filename,'(a,i3.3)')'min',ifile
!            inquire(file=folder//'/'//filename//'.xyz',exist=xyzexists)
!            inquire(file=folder//'/'//filename//'.ascii',exist=asciiexists)
!            if(.not.(xyzexists.or.asciiexists))exit
!            call deallocate_atomic_structure(atoms%astruct)
!            call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
!            call vcopy(3 * atoms%astruct%nat,atoms%astruct%rxyz(1,1),1,rxyz(1,1), 1)
!            call vcopy(3 * atoms%astruct%nat,outs%fxyz(1,1),1,fxyz(1,1), 1)
!            call energyandforces(atoms%astruct%nat,atoms%astruct%cell_dim,rxyz,fxyz,energy)
! 
!            do i=1,atoms%astruct%nat
!                minmode(1,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
!                minmode(2,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
!                minmode(3,i)=2.0_gp*(real(builtin_rand(idum),gp)-0.5_gp)
!            enddo
!
!            call opt_curv(saddle_imode,atoms%astruct%nat,atoms%astruct%cell_dim,&
!                 saddle_alpha0_rot,saddle_curvgraddiff,saddle_nit_rot,saddle_nhistx_rot,&
!                 rxyz,fxyz,minmode,curv,gradrot,saddle_tolf,displ,ec,&
!                 .false.,converged,iconnect,nbond,atoms%astruct%atomnames,&
!                 saddle_alpha_stretch0,saddle_maxcurvrise,saddle_cutoffratio)
!       enddo
!    enddo




!    if (iproc==0) call yaml_close_stream(unit=usaddle)



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
        call deallocate_global_output(outs)
    endif

    call f_free(work)
    call f_free(minmode)
    call f_free(rxyz)
    call f_free(fxyz) 
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

    


    if(iproc==0)call yaml_map('(MHGPS) Calls to energy and forces',nint(ef_counter))
    if(iproc==0)call yaml_map('(MHGPS) Run finished at',yaml_date_and_time_toa())
    call f_lib_finalize()
end program

subroutine cal_hessian_fd(iproc,nat,alat,pos,hess)
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
    real(8) :: h,rlarge,twelfth,twothird,etot,cmx,cmy,cmz,shift,dm,tt,s
    integer :: i,j,k,lwork,info

    !allocate(hess(3*nat,3*nat))
    allocate(tpos(3*nat))
    allocate(grad(3*nat))
    allocate(eval(3*nat))

    lwork=1000*nat
    allocate(work(lwork))

    !h=1.d-1
    !h=7.5d-2
    !h=5.d-2
!    h=1.d-3
    h=1.d-2
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
        call energyandforces(nat,alat,tpos,grad,etot)
        do j=1,3*nat
            hess(j,i)=twelfth*grad(j)
        enddo
        !if(iproc==0) write(*,*) 'ALIREZA-6',i,iat
        !-----------------------------------------
        tpos(i)=tpos(i)+h
        call energyandforces(nat,alat,tpos,grad,etot)
        do j=1,3*nat
        hess(j,i)=hess(j,i)-twothird*grad(j)
        enddo
        !-----------------------------------------
        tpos(i)=tpos(i)+2*h
        call energyandforces(nat,alat,tpos,grad,etot)
        do j=1,3*nat
        hess(j,i)=hess(j,i)+twothird*grad(j)
        enddo
        !-----------------------------------------
        tpos(i)=tpos(i)+h
        call energyandforces(nat,alat,tpos,grad,etot)
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
end subroutine cal_hessian_fd

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
    use module_saddle
    use module_minimizers
    implicit none
    integer :: bigdft_get_number_of_atoms,bigdft_get_number_of_orbitals
    character(len=*), parameter :: subname='mhgps'
    integer :: iproc=0,nproc,igroup,ngroups
    character(len=8) :: folder
    character(len=6) :: filename
    integer :: ifolder, ifile
    logical :: xyzexists,asciiexists

    character(len=60) :: run_id
    integer :: ierr, nconfig

    real(gp) :: count,count_sd,displ,ec
    logical :: converged
    real(gp), allocatable :: rcov(:)

integer :: ncount_bigdft
logical :: fail

real(gp) :: curv
real(gp), allocatable, dimension(:,:) :: gradrot

!simple atomic datastructre
integer :: nat
real(gp),allocatable :: rxyz(:,:),fxyz(:,:),rotforce(:,:)
real(gp) :: energy
integer :: i
integer :: idum=0
real(kind=4) :: tt,builtin_rand

    ifolder=1
    ifile=1



    call f_lib_initialize()

    call read_input()
    if(efmethod=='BIGDFT')then
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
        write(filename,'(a,i3.3)')'min',ifile
        call user_dict_from_files(user_inputs, trim(run_id)//trim(bigdft_run_id_toa()), &
           & folder//'/'//filename//trim(bigdft_run_id_toa()), bigdft_mpi)
        call inputs_from_dict(inputs_opt, atoms, user_inputs)
        call dict_free(user_inputs)
        call init_global_output(outs, atoms%astruct%nat)
        call init_restart_objects(bigdft_mpi%iproc,inputs_opt,atoms,rst,subname)
        call run_objects_nullify(runObj)
        call run_objects_associate(runObj, inputs_opt, atoms, rst)

    elseif(efmethod=='LJ')then
        write(folder,'(a,i3.3)')'input',ifolder
        write(filename,'(a,i3.3)')'min',ifile
        call deallocate_atomic_structure(atoms%astruct)
        call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
        call init_global_output(outs, atoms%astruct%nat)
        if(iproc==0) call print_logo_mhgps()
    else
        call yaml_warning('Following method for evaluation of energies and forces is unknown: '//trim(adjustl(efmethod)))
        stop
    endif

    !allocate more arrays
    minmode  = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),id='minmode')
    rxyz     = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),id='rxyz')
    fxyz     = f_malloc((/ 1.to.3, 1.to.atoms%astruct%nat/),id='rxyz')
    rcov     = f_malloc((/ 1.to.atoms%astruct%nat/),id='rcov')
    iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
allocate(gradrot(3,atoms%astruct%nat))
allocate(rotforce(3,atoms%astruct%nat))

    !if in biomode, determine bonds betweens atoms once and for all (it is
    !assuemed that all conifugrations over which will be iterated have the same
    !bonds)
    if(saddle_biomode)then
        call findbonds(atoms%astruct%nat,rcov,atoms%astruct%rxyz,nbond,iconnect)
    endif


    call give_rcov(iproc,atoms,atoms%astruct%nat,rcov)


    do ifolder = 1,999
        do ifile = 1,999
            write(folder,'(a,i3.3)')'input',ifolder
            write(filename,'(a,i3.3)')'min',ifile
            inquire(file=folder//'/'//filename//'.xyz',exist=xyzexists)
            inquire(file=folder//'/'//filename//'.ascii',exist=asciiexists)
            if(.not.(xyzexists.or.asciiexists))exit
            call deallocate_atomic_structure(atoms%astruct)
            call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
!            if(ifolder/=1.or.ifile/=1)then
!                runObj%inputs%inputPsiId=1           
!            else
!                runObj%inputs%inputPsiId=0
!            endif
rxyz=atoms%astruct%rxyz
fxyz=outs%fxyz
            call energyandforces(atoms%astruct%nat,atoms%astruct%cell_dim,rxyz,fxyz,energy)
write(*,*)'BASTIAN',nat
write(*,*)'BASTIAN',rxyz
write(*,*)'BASTIAN',fxyz
 
!            if(iproc==0)write(*,*)'Bastian',outs%energy
!           count=0.0_gp;count_sd=0.0_gp;displ=0.0_gp;ec=0.0_gp
!           call random_seed
!           call random_number(minmode)
do i=1,atoms%astruct%nat
minmode(1,i)=real(builtin_rand(idum),gp)
minmode(2,i)=real(builtin_rand(idum),gp)
minmode(3,i)=real(builtin_rand(idum),gp)
enddo

           minmode=2.d0*(minmode-0.5d0)
!           runObj%inputs%inputPsiId=1
!!           call findsad(saddle_imode,atoms%astruct%nat,atoms%astruct%cell_dim,rcov,saddle_alpha0_trans,saddle_alpha0_rot,saddle_curvgraddiff,saddle_nit_trans,&
!!           saddle_nit_rot,saddle_nhistx_trans,saddle_nhistx_rot,saddle_tolc,saddle_tolf,saddle_tightenfac,saddle_rmsdispl0,&
!!           saddle_trustr,atoms%astruct%rxyz,outs%energy,outs%fxyz,minmode,saddle_fnrmtol,count,count_sd,displ,ec,&
!!           converged,atoms%astruct%atomnames,nbond,iconnect,saddle_alpha_stretch0,saddle_recompIfCurvPos,saddle_maxcurvrise,saddle_cutoffratio)
!call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
!call minimizer_sbfgs(runObj,outs,nproc,iproc,1,ncount_bigdft,fail)
!rxyz=atoms%astruct%rxyz
!fxyz=outs%fxyz
!call curvgrad(atoms%astruct%nat,atoms%astruct%cell_dim,1.d-3,rxyz,fxyz,minmode,curv,rotforce,1,ec)
!rxyz=atoms%astruct%rxyz
!fxyz=outs%fxyz
call opt_curv(saddle_imode,atoms%astruct%nat,atoms%astruct%cell_dim,saddle_alpha0_rot,saddle_curvgraddiff,saddle_nit_rot,saddle_nhistx_rot,rxyz,fxyz,minmode,curv,gradrot,saddle_tolf&
            &,count,count_sd,displ,ec,.false.,converged,iconnect,nbond,atoms%astruct%atomnames,saddle_alpha_stretch0,saddle_maxcurvrise,saddle_cutoffratio)
        enddo
    enddo







    !finalize (dealloctaion etc...)
    if(efmethod=='BIGDFT')then
        call free_restart_objects(rst,subname)
        call deallocate_atoms_data(atoms)
        call deallocate_global_output(outs)
        call run_objects_free_container(runObj)
        call free_input_variables(inputs_opt)
        call bigdft_finalize(ierr)
    elseif(efmethod=='LJ')then
        call deallocate_atoms_data(atoms)
        call deallocate_global_output(outs)
    endif

    call f_free(minmode)
    call f_free(rxyz)
    call f_free(fxyz)
    call f_free(rcov)
    call f_free(iconnect)

    call f_lib_finalize()
end program

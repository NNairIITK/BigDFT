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

real(8) :: count,count_sd,displ,ec
logical :: converged
real(8), allocatable :: rcov(:)

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
    !        TODO
    else
        call yaml_warning('Following method for evaluation of energies and forces is unknown: '//trim(adjustl(efmethod)))
        stop
    endif

    !allocate more
    allocate(iconnect(2,nbond))
    allocate(minmode(3,atoms%astruct%nat))
    allocate(rcov(atoms%astruct%nat))

    call give_rcov(iproc,atoms,atoms%astruct%nat,rcov)



    do ifolder = 1,999
        do ifile = 1,999
            write(folder,'(a,i3.3)')'input',ifolder
            write(filename,'(a,i3.3)')'min',ifile
            inquire(file=folder//'/'//filename//'.xyz',exist=xyzexists)
            inquire(file=folder//'/'//filename//'.ascii',exist=asciiexists)
            if(.not.(xyzexists.or.asciiexists))exit
write(*,*)folder//'/'//filename
            call deallocate_atomic_structure(atoms%astruct)
            call read_atomic_file(folder//'/'//filename,iproc,atoms%astruct)
!            if(ifolder/=1.or.ifile/=1)then
!                runObj%inputs%inputPsiId=1           
!            else
!                runObj%inputs%inputPsiId=0
!            endif
!            call energyandforces(atoms%astruct%nat,atoms%astruct%cell_dim,atoms%astruct%rxyz,outs%fxyz,outs%energy) 
!            if(iproc==0)write(*,*)'Bastian',outs%energy
           count=0.0_gp;count_sd=0.0_gp;displ=0.0_gp;ec=0.0_gp
           call random_seed
           call random_number(minmode)
           minmode=2.d0*(minmode-0.5d0)
           call findsad(saddle_imode,atoms%astruct%nat,atoms%astruct%cell_dim,rcov,saddle_alpha0_trans,saddle_alpha0_rot,saddle_curvgraddiff,saddle_nit_trans,&
           saddle_nit_rot,saddle_nhistx_trans,saddle_nhistx_rot,saddle_tolc,saddle_tolf,saddle_tightenfac,saddle_rmsdispl0,&
           saddle_trustr,atoms%astruct%rxyz,outs%energy,outs%fxyz,minmode,saddle_fnrmtol,count,count_sd,displ,ec,&
           converged,atoms%astruct%atomnames,nbond,iconnect,saddle_alpha_stretch0,saddle_recompIfCurvPos,saddle_maxcurvrise,saddle_cutoffratio)
        enddo
    enddo








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
    call f_lib_finalize()

end program

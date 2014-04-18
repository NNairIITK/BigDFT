!> @file 
!!   Routines to use BigDFT as a blackbox
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
 

!> Routine to use BigDFT as a blackbox
subroutine call_bigdft(runObj,outs,nproc,iproc,infocode)
  use module_base
  use module_types
  use module_interfaces, except_this_one => call_bigdft
  use yaml_output
  use communications_base
  implicit none
  integer, intent(in) :: iproc,nproc
  type(run_objects), intent(inout) :: runObj
  type(DFT_global_output), intent(inout) :: outs
  integer, intent(inout) :: infocode

  !local variables
  character(len=*), parameter :: subname='call_bigdft'
  character(len=40) :: comment
  logical :: exists
  integer :: i_stat,i_all,ierr,inputPsiId_orig,iat,iorb,istep,i,jproc
  real(gp) :: maxdiff

  !temporary interface, not needed anymore since all arguments are structures
!!$  interface
!!$     subroutine cluster(nproc,iproc,atoms,rxyz,energy,energs,fxyz,strten,fnoise,pressure,&
!!$          KSwfn,tmb,&!psi,Lzd,gaucoeffs,gbd,orbs,
!!$          rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
!!$       use module_base
!!$       use module_types
!!$       implicit none
!!$       integer, intent(in) :: nproc,iproc
!!$       integer, intent(out) :: infocode
!!$       real(gp), intent(inout) :: hx_old,hy_old,hz_old
!!$       type(input_variables), intent(in) :: in
!!$       !type(local_zone_descriptors), intent(inout) :: Lzd
!!$       type(atoms_data), intent(inout) :: atoms
!!$       !type(gaussian_basis), intent(inout) :: gbd
!!$       !type(orbitals_data), intent(inout) :: orbs
!!$       type(GPU_pointers), intent(inout) :: GPU
!!$       type(DFT_wavefunction), intent(inout) :: KSwfn,tmb
!!$       type(energy_terms), intent(out), target :: energs
!!$       real(gp), intent(out) :: energy,fnoise,pressure
!!$       real(gp), dimension(3,atoms%astruct%nat), target, intent(inout) :: rxyz
!!$       real(gp), dimension(6), intent(out) :: strten
!!$       real(gp), dimension(3,atoms%astruct%nat), intent(out) :: fxyz
!!$     END SUBROUTINE cluster
!!$  end interface

  !put a barrier for all the processes
  call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  call f_routine(id=subname)

  call check_array_consistency(maxdiff, nproc, runObj%atoms%astruct%rxyz(1,1), &
       & 3 * runObj%atoms%astruct%nat, bigdft_mpi%mpi_comm)
  if (iproc==0 .and. maxdiff > epsilon(1.0_gp)) &
       call yaml_warning('Input positions not identical! '//&
       '(difference:'//trim(yaml_toa(maxdiff))//' )')

  !fill the rxyz array with the positions
  !wrap the atoms in the periodic directions when needed
  do iat=1,runObj%atoms%astruct%nat
     if (runObj%atoms%astruct%geocode == 'P') then
        runObj%rst%rxyz_new(1,iat)=modulo(runObj%atoms%astruct%rxyz(1,iat),runObj%atoms%astruct%cell_dim(1))
        runObj%rst%rxyz_new(2,iat)=modulo(runObj%atoms%astruct%rxyz(2,iat),runObj%atoms%astruct%cell_dim(2))
        runObj%rst%rxyz_new(3,iat)=modulo(runObj%atoms%astruct%rxyz(3,iat),runObj%atoms%astruct%cell_dim(3))
     else if (runObj%atoms%astruct%geocode == 'S') then
        runObj%rst%rxyz_new(1,iat)=modulo(runObj%atoms%astruct%rxyz(1,iat),runObj%atoms%astruct%cell_dim(1))
        runObj%rst%rxyz_new(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
        runObj%rst%rxyz_new(3,iat)=modulo(runObj%atoms%astruct%rxyz(3,iat),runObj%atoms%astruct%cell_dim(3))
     else if (runObj%atoms%astruct%geocode == 'F') then
        runObj%rst%rxyz_new(1,iat)=runObj%atoms%astruct%rxyz(1,iat)
        runObj%rst%rxyz_new(2,iat)=runObj%atoms%astruct%rxyz(2,iat)
        runObj%rst%rxyz_new(3,iat)=runObj%atoms%astruct%rxyz(3,iat)
     end if
  end do

  !assign the verbosity of the output
  !the verbose variables is defined in module_base
  verbose=runObj%inputs%verbosity

  ! Use the restart for the linear scaling version... probably to be modified.
  if(runObj%inputs%inputPsiId==1) then
      if (runObj%rst%version == LINEAR_VERSION) then
          runObj%inputs%inputPsiId=101
          do iorb=1,runObj%rst%tmb%orbs%norb
              if (runObj%inputs%lin%locrad_lowaccuracy(iorb) /=  runObj%inputs%lin%locrad_highaccuracy(iorb))then
                  stop 'ERROR: at the moment the radii for low and high accuracy must be the same &
                        &when using the linear restart!'
              end if
          end do
      end if
  end if
  inputPsiId_orig=runObj%inputs%inputPsiId

  loop_cluster: do
     !allocate history container if it has not been done
     if (runObj%inputs%wfn_history > 1  .and. .not. associated(runObj%rst%KSwfn%oldpsis)) then
        allocate(runObj%rst%KSwfn%oldpsis(0:runObj%inputs%wfn_history+1))
        runObj%rst%KSwfn%istep_history=0
        do istep=0,runObj%inputs%wfn_history+1
          runObj%rst%KSwfn%oldpsis(istep)=old_wavefunction_null() 
        end do
     end if

     if (runObj%inputs%inputPsiId == 0 .and. associated(runObj%rst%KSwfn%psi)) then
        i_all=-product(shape(runObj%rst%KSwfn%psi))*kind(runObj%rst%KSwfn%psi)
        deallocate(runObj%rst%KSwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(runObj%rst%KSwfn%orbs%eval))*kind(runObj%rst%KSwfn%orbs%eval)
        deallocate(runObj%rst%KSwfn%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(runObj%rst%KSwfn%Lzd%Glr%wfd)
     end if
     !experimental, finite difference method for calculating forces on particular quantities
     inquire(file='input.finite_difference_forces',exist=exists)
     if (exists) then
        runObj%inputs%last_run=1 !do the last_run things nonetheless
        runObj%inputs%inputPsiId=0 !the first run always restart from IG
        !experimental_modulebase_var_onlyfion=.true. !put only ionic forces in the forces
     end if
     call cluster(nproc,iproc,runObj%atoms,runObj%rst%rxyz_new, runObj%radii_cf, &
          & outs%energy, outs%energs, outs%fxyz, outs%strten, outs%fnoise, outs%pressure,&
          runObj%rst%KSwfn,runObj%rst%tmb,&!psi,runObj%rst%Lzd,runObj%rst%gaucoeffs,runObj%rst%gbd,runObj%rst%orbs,&
          runObj%rst%rxyz_old,runObj%rst%hx_old,runObj%rst%hy_old,runObj%rst%hz_old,runObj%inputs,runObj%rst%GPU,infocode)
     !save the new atomic positions in the rxyz_old array
     do iat=1,runObj%atoms%astruct%nat
        runObj%rst%rxyz_old(1,iat)=runObj%rst%rxyz_new(1,iat)
        runObj%rst%rxyz_old(2,iat)=runObj%rst%rxyz_new(2,iat)
        runObj%rst%rxyz_old(3,iat)=runObj%rst%rxyz_new(3,iat)
     enddo
     if (exists) then
        call forces_via_finite_differences(iproc,nproc,runObj%atoms,runObj%inputs, &
             & outs%energy,outs%fxyz,outs%fnoise,runObj%rst,infocode)
     end if

     if (runObj%inputs%inputPsiId==1 .and. infocode==2) then
        if (runObj%inputs%gaussian_help) then
           runObj%inputs%inputPsiId=11
        else
           runObj%inputs%inputPsiId=0
        end if
     else if (runObj%inputs%inputPsiId==101 .and. infocode==2) then
         ! problems after restart for linear version
         runObj%inputs%inputPsiId=100
     else if ((runObj%inputs%inputPsiId==1 .or. runObj%inputs%inputPsiId==0) .and. infocode==1) then
        !runObj%inputs%inputPsiId=0 !better to diagonalise than to restart an input guess
        runObj%inputs%inputPsiId=1
        !if (iproc==0) then
        !   call yaml_warning('Self-consistent cycle did not meet convergence criteria')
        !   write(*,*)&
        !        &   ' WARNING: Self-consistent cycle did not meet convergence criteria'
        !end if
        exit loop_cluster
     else if (runObj%inputs%inputPsiId == 0 .and. infocode==3) then
        if (iproc == 0) then
           write( *,'(1x,a)')'Convergence error, cannot proceed.'
           write( *,'(1x,a)')' writing positions in file posfail.xyz then exiting'
           write(comment,'(a)')'UNCONVERGED WF '
           !call wtxyz('posfail',energy,rxyz,atoms,trim(comment))

           call write_atomic_file("posfail",outs%energy,runObj%rst%rxyz_new,runObj%atoms,trim(comment))

        end if

        i_all=-product(shape(runObj%rst%KSwfn%psi))*kind(runObj%rst%KSwfn%psi)
        deallocate(runObj%rst%KSwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        i_all=-product(shape(runObj%rst%KSwfn%orbs%eval))*kind(runObj%rst%KSwfn%orbs%eval)
        deallocate(runObj%rst%KSwfn%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        call deallocate_wfd(runObj%rst%KSwfn%Lzd%Glr%wfd)

        !finalize memory counting (there are still at least positions and the forces allocated)
        call memocc(0,0,'count','stop')

        if (nproc > 1) call MPI_FINALIZE(ierr)

        stop 'unnormal end'
     else
        exit loop_cluster
     end if

  end do loop_cluster

  !preserve the previous value
  runObj%inputs%inputPsiId=inputPsiId_orig

  !put a barrier for all the processes
  call f_release_routine()
  call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)

END SUBROUTINE call_bigdft

!>  Main routine which does self-consistent loop.
!!  Does not parse input file and no geometry optimization.
!!  Does an electronic structure calculation. 
!!  Output is the total energy and the forces 
!!   @warning psi, keyg, keyv and eval should be freed after use outside of the routine.
!!   @param inputPsiId 
!!           - 0 : compute input guess for Psi by subspace diagonalization of atomic orbitals
!!           - 1 : read waves from argument psi, using n1, n2, n3, hgrid and rxyz_old
!!                 as definition of the previous system.
!!           - 2 : read waves from disk
!!   @param infocode -> encloses some information about the status of the run
!!           - 0 run succesfully succeded
!!           - 1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
!!               forces may be meaningless   
!!           - 2 (present only for inputPsiId=1) gnrm of the first iteration > 1 AND growing in
!!               the second iteration OR grnm 1st >2.
!!               Input wavefunctions need to be recalculated. Routine exits.
!!           - 3 (present only for inputPsiId=0) gnrm > 4. SCF error. Routine exits.
subroutine cluster(nproc,iproc,atoms,rxyz,radii_cf,energy,energs,fxyz,strten,fnoise,pressure,&
     KSwfn,tmb,rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
  use module_base
  use module_types
  use module_interfaces
  use module_fragments
  use constrained_dft
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use module_xc
  use communications_init, only: orbitals_communicators
!  use vdwcorrection
  use yaml_output
  use gaussians, only: gaussian_basis
  use psp_projectors
  use sparsematrix_base, only: sparse_matrix_null
  use sparsematrix_init, only: init_sparse_matrix, check_kernel_cutoff
  use sparsematrix, only: check_matrix_compression
  implicit none
  integer, intent(in) :: nproc,iproc
  real(gp), intent(inout) :: hx_old,hy_old,hz_old
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(inout) :: atoms
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_wavefunction), intent(inout) :: KSwfn, tmb
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz_old
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  integer, intent(out) :: infocode
  type(energy_terms), intent(out) :: energs
  real(gp), intent(out) :: energy,fnoise,pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), dimension(3,atoms%astruct%nat), intent(out) :: fxyz
  !local variables
  character(len=*), parameter :: subname='cluster'
  character(len=5) :: gridformat, wfformat
  logical :: refill_proj, calculate_dipole !,potential_from_disk=.false.
  logical :: DoDavidson,DoLastRunThings=.false.
  integer :: nvirt,norbv
  integer :: i, input_wf_format, output_denspot
  integer :: n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i
  integer :: iat,i_all,i_stat,ierr,inputpsi,igroup,ikpt,nproctiming,ifrag
  real :: tcpu0,tcpu1
  real(kind=8) :: tel
  type(grid_dimensions) :: d_old
  type(wavefunctions_descriptors) :: wfd_old
  type(local_zone_descriptors) :: lzd_old
  type(DFT_PSP_projectors) :: nlpsp
  type(DFT_wavefunction) :: VTwfn !< Virtual wavefunction
  type(DFT_wavefunction) :: tmb_old
  !!type(DFT_wavefunction) :: tmb
  type(system_fragment), dimension(:), pointer :: ref_frags
  type(cdft_data) :: cdft
  real(gp), dimension(3) :: shift
  real(dp), dimension(6) :: ewaldstr,xcstr
  real(gp), dimension(:,:), allocatable :: thetaphi,band_structure_eval
  real(gp), dimension(:,:), pointer :: fdisp,fion,fpulay
  ! Charge density/potential,ionic potential, pkernel
  type(DFT_local_fields) :: denspot
  type(DFT_optimization_loop) :: optLoop
  real(gp), dimension(:), allocatable:: denspot0
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(wp), dimension(:), pointer :: psi_old
  type(memory_estimation) :: mem
  ! PSP projectors 
  real(kind=8), dimension(:), pointer :: gbd_occ!,rhocore
  ! Variables for the virtual orbitals and band diagram.
  integer :: nkptv, nvirtu, nvirtd
  real(gp), dimension(:), allocatable :: wkptv

  !Variables for WVL+PAW
  integer:: iatyp
  type(gaussian_basis),dimension(max(atoms%astruct%ntypes,0))::proj_G
  type(rholoc_objects)::rholoc_tmp

  ! testing
  real(kind=8),dimension(:,:),pointer :: locregcenters
  integer :: ilr, nlr, iorb, jorb, ioffset, linear_iscf
  character(len=20) :: comment

  !debug
  !real(kind=8) :: ddot

  call f_routine(id=subname)


  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be 
  !reassigned inside the structure

  if (iproc == 0) then
     !start a new document in the beginning of the output, if the document is closed before
     call yaml_new_document()
     !write( *,'(1x,a,1x,i0)') &
     !     &   '===================== BigDFT Wavefunction Optimization =============== inputPsiId=',&
     !     in%inputPsiId
     call print_dft_parameters(in,atoms)
  end if

  !Time initialization
  if (verbose > 2) then
     nproctiming=-nproc !timing in debug mode
  else
     nproctiming=nproc
  end if
  !call timing(nproctiming,trim(in%dir_output)//'time.yaml','IN')
  call f_timing_reset(filename=trim(in%dir_output)//'time.yaml',master=iproc==0,&
       verbose_mode=verbose>2 .and. nproc>1)
  call cpu_time(tcpu0)
  call system_clock(ncount0,ncount_rate,ncount_max)

  !Nullify for new input guess
  call nullify_local_zone_descriptors(lzd_old)
  !call nullify_wavefunctions_descriptors(wfd_old)
  call nullify_wfd(wfd_old)
  ! We save the variables that defined the previous psi if the restart is active
  inputpsi = in%inputPsiId
  if (in%inputPsiId == INPUT_PSI_MEMORY_WVL) then
     if (associated(KSwfn%psi)) then
        !regenerate grid spacings (this would not be needed if hgrids is in Lzd)
        if (atoms%astruct%geocode == 'P') then
           call correct_grid(atoms%astruct%cell_dim(1),hx_old,KSwfn%Lzd%Glr%d%n1)
           call correct_grid(atoms%astruct%cell_dim(2),hy_old,KSwfn%Lzd%Glr%d%n2)
           call correct_grid(atoms%astruct%cell_dim(3),hz_old,KSwfn%Lzd%Glr%d%n3)
        else if (atoms%astruct%geocode == 'S') then 
           call correct_grid(atoms%astruct%cell_dim(1),hx_old,KSwfn%Lzd%Glr%d%n1)
           call correct_grid(atoms%astruct%cell_dim(3),hz_old,KSwfn%Lzd%Glr%d%n3)
        end if

        call copy_local_zone_descriptors(KSwfn%Lzd, lzd_old, subname)

        !if the history is bigger than two, create the workspace to store the wavefunction
        if (in%wfn_history > 2) then
           call old_wavefunction_set(KSwfn%oldpsis(in%wfn_history+1),&
                atoms%astruct%nat,KSwfn%orbs%norbp*KSwfn%orbs%nspinor,&
                KSwfn%Lzd,rxyz_old,KSwfn%psi)
           !to maintain the same treatment destroy wfd afterwards (to be unified soon)
           !deallocation
           call deallocate_wfd(KSwfn%Lzd%Glr%wfd)
        else
           call copy_old_wavefunctions(nproc,KSwfn%orbs,&
                KSwfn%Lzd%Glr%d%n1,KSwfn%Lzd%Glr%d%n2,KSwfn%Lzd%Glr%d%n3,&
                KSwfn%Lzd%Glr%wfd,KSwfn%psi,d_old%n1,d_old%n2,d_old%n3,wfd_old,psi_old)
        end if
        !already here due to new input guess
        call deallocate_bounds(KSwfn%Lzd%Glr%geocode, KSwfn%Lzd%Glr%hybrid_on, KSwfn%lzd%glr%bounds, subname)
     else
        inputpsi = INPUT_PSI_LCAO
     end if
  else if (in%inputPsiId == INPUT_PSI_MEMORY_GAUSS) then
     if (associated(KSwfn%psi)) then
        !deallocate wavefunction and descriptors for placing the gaussians
        
        call deallocate_wfd(KSwfn%Lzd%Glr%wfd)

        i_all=-product(shape(KSwfn%psi))*kind(KSwfn%psi)
        deallocate(KSwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
     else
        inputpsi = INPUT_PSI_LCAO
     end if
  else if (in%inputPsiId == INPUT_PSI_MEMORY_LINEAR .and. associated(KSwfn%psi)) then
     if (associated(KSwfn%psi) .and. associated(tmb%psi)) then
        call copy_tmbs(iproc, tmb, tmb_old, subname)
        call destroy_DFT_wavefunction(tmb)
        i_all=-product(shape(KSwfn%psi))*kind(KSwfn%psi)
        deallocate(KSwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        call deallocate_wfd(KSwfn%Lzd%Glr%wfd)
     else
        inputpsi = INPUT_PSI_LINEAR_AO
     end if
  end if

  ! Setup all descriptors and allocate what should be.
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. &
      inputpsi == INPUT_PSI_MEMORY_LINEAR .or. &
      inputpsi == INPUT_PSI_DISK_LINEAR) then
      if (in%explicit_locregcenters) then
          locregcenters=f_malloc_ptr((/3,atoms%astruct%nat/),id=' locregcenters')
          open(unit=123, file='locregcenters.xyz')
          read(123,*) nlr
          if (nlr/=atoms%astruct%nat) stop 'ERROR: wrong nlr'
          read(123,*) comment
          do ilr=1,nlr
              read(123,*) comment, locregcenters(1,ilr), locregcenters(2,ilr), locregcenters(3,ilr)
          end do
      else
          locregcenters=f_malloc_ptr((/3,atoms%astruct%nat/),id=' locregcenters')
          locregcenters = rxyz
      end if
  end if

  if(inputpsi == INPUT_PSI_MEMORY_LINEAR) then
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,.false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,shift,radii_cf,ref_frags,denspot,locregcenters,tmb_old%orbs%inwhichlocreg,tmb_old%orbs%onwhichatom)
  else if(inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR) then
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,.false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,shift,radii_cf,ref_frags,denspot,locregcenters)
  else
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,&
         & .false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,shift,radii_cf,ref_frags,denspot)
  end if

  !memory estimation, to be rebuilt in a more modular way
  call MemoryEstimator(nproc,in%idsx,KSwfn%Lzd%Glr,&
       KSwfn%orbs%norb,KSwfn%orbs%nspinor,KSwfn%orbs%nkpts,&
       nlpsp%nprojel,in%nspin,in%itrpmax,in%iscf,mem)
  if (iproc==0 .and. verbose > 0) call print_memory_estimation(mem)

  if (in%lin%fragment_calculation .and. inputpsi == INPUT_PSI_DISK_LINEAR) then
     call output_fragment_rotations(iproc,nproc,atoms%astruct%nat,rxyz,1,trim(in%dir_output),in%frag,ref_frags)
     !call mpi_finalize(i_all)
     !stop
  end if

  ! temporary, really want to just initialize it here rather than copy
  ! but still need to move all cubic references to KSwfn%orbs%npsidim to just KSwfn%npsidim
  KSwfn%npsidim_orbs = KSwfn%orbs%npsidim_orbs
  KSwfn%npsidim_comp = KSwfn%orbs%npsidim_comp

  ! We complete here the definition of DFT_wavefunction structures.
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
      .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     !!call init_p2p_tags(nproc)
     !!tag=0

     call kswfn_init_comm(tmb, in, atoms, denspot%dpbox, iproc, nproc)
     call init_foe(iproc, nproc, tmb%lzd, atoms%astruct, in, KSwfn%orbs, tmb%orbs, tmb%foe_obj, .true.)
     call increase_FOE_cutoff(iproc, nproc, tmb%lzd, atoms%astruct, in, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.true.)

     call create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot,nlpsp,in, atoms, rxyz, .false.)

     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%ham_descr%lzd, atoms%astruct, &
          in%store_index, imode=1, smat=tmb%linmat%ham)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%ham_descr%lzd, atoms%astruct, &
          in%store_index, imode=1, smat=tmb%linmat%m)

     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%m)

     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, atoms%astruct, &
          in%store_index, imode=1, smat=tmb%linmat%ovrlp)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, atoms%astruct, &
          in%store_index, imode=1, smat=tmb%linmat%s)

     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ovrlp)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%s)

     if (in%check_matrix_compression) then
         if (iproc==0) call yaml_open_map('Checking Compression/Uncompression of small sparse matrices')
         call check_matrix_compression(iproc,tmb%linmat%ham)
         call check_matrix_compression(iproc,tmb%linmat%ovrlp)
         if (iproc ==0) call yaml_close_map()
     end if




     ! check the extent of the kernel cutoff (must be at least shamop radius)
     call check_kernel_cutoff(iproc, tmb%orbs, atoms, tmb%lzd)

     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, atoms%astruct, &
          in%store_index, imode=2, smat=tmb%linmat%denskern_large)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, atoms%astruct, &
          in%store_index, imode=2, smat=tmb%linmat%l)

     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%denskern_large)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%l)

     !call nullify_sparse_matrix(tmb%linmat%inv_ovrlp_large)
     tmb%linmat%inv_ovrlp_large=sparse_matrix_null()
     call sparse_copy_pattern(tmb%linmat%denskern_large, tmb%linmat%inv_ovrlp_large, iproc, subname)


     if (in%check_matrix_compression) then
         if (iproc==0) call yaml_open_map('Checking Compression/Uncompression of large sparse matrices')
         call check_matrix_compression(iproc,tmb%linmat%denskern_large)
         if (iproc ==0) call yaml_close_map()
     end if




     if (in%check_sumrho>0) then
         call check_communication_potential(iproc,denspot,tmb)
         call check_communication_sumrho(iproc, nproc, tmb%orbs, tmb%lzd, tmb%collcom_sr, &
              denspot, tmb%linmat%denskern_large, in%check_sumrho)
     end if

     if (in%lin%scf_mode/=LINEAR_FOE .or. in%lin%pulay_correction .or.  in%lin%new_pulay_correction .or. &
         (in%lin%plotBasisFunctions /= WF_FORMAT_NONE) .or. in%lin%diag_end) then
        allocate(tmb%coeff(tmb%orbs%norb,tmb%orbs%norb), stat=i_stat)
        call memocc(i_stat, tmb%coeff, 'tmb%coeff', subname)
     else
        nullify(tmb%coeff)
     end if

     allocate(denspot0(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)), stat=i_stat)
     call memocc(i_stat, denspot0, 'denspot0', subname)
  else
     allocate(denspot0(1+ndebug), stat=i_stat)
     call memocc(i_stat, denspot0, 'denspot0', subname)
  end if

  !the lookup tables for the application of the nonlocal potential can be created from now on

  optLoop%iscf = in%iscf
  optLoop%itrpmax = in%itrpmax
  optLoop%nrepmax = in%nrepmax
  optLoop%itermax = in%itermax
  optLoop%gnrm_cv = in%gnrm_cv
  optLoop%rpnrm_cv = in%rpnrm_cv
  optLoop%gnrm_startmix = in%gnrm_startmix
  optLoop%itrp = 0
  optLoop%itrep = 0
  optLoop%iter = 0
  optLoop%infocode = 0

  call system_signaling(iproc,in%signaling,in%gmainloop,&
       & KSwfn,tmb,energs,denspot,optloop,&
       & atoms%astruct%ntypes,radii_cf,in%crmult,in%frmult)

  !variables substitution for the PSolver part
  n1=KSwfn%Lzd%Glr%d%n1
  n2=KSwfn%Lzd%Glr%d%n2
  n3=KSwfn%Lzd%Glr%d%n3

  !calculate the rhocore contribution to the energy value
  if (associated(denspot%rho_C)) then
     !calculate the XC energy of rhocore, use the rhov array as a temporary variable
     !use Vxc and other quantities as local variables
     call xc_init_rho(denspot%xc, denspot%dpbox%nrhodim,denspot%rhov,1)
     denspot%rhov=1.d-16
     call XC_potential(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
          denspot%pkernel%mpi_env%mpi_comm,&
          denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),denspot%xc,&
          denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
          denspot%rhov,energs%excrhoc,tel,KSwfn%orbs%nspin,denspot%rho_C,denspot%V_XC,xcstr)
     if (iproc==0) call yaml_map('Value for Exc[rhoc]',energs%excrhoc)
     !if (iproc==0) write(*,*)'value for Exc[rhoc]',energs%excrhoc
  end if

  !here calculate the ionic energy and forces accordingly
  call IonicEnergyandForces(iproc,nproc,denspot%dpbox,atoms,in%elecfield,rxyz,&
       energs%eion,fion,in%dispersion,energs%edisp,fdisp,ewaldstr,&
       n1,n2,n3,denspot%V_ext,denspot%pkernel,denspot%psoffset)
  !calculate effective ionic potential, including counter ions if any.
  call nullify_rholoc_objects(rholoc_tmp)
  call createEffectiveIonicPotential(iproc,nproc,(iproc == 0),in,atoms,rxyz,shift,KSwfn%Lzd%Glr,&
       denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
       denspot%dpbox,denspot%pkernel,denspot%V_ext,in%elecfield,denspot%psoffset,&
       rholoc_tmp)
  if (denspot%c_obj /= 0) then
     call denspot_emit_v_ext(denspot, iproc, nproc)
  end if

  norbv=abs(in%norbv)
  if (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
      in%inputPsiId == INPUT_PSI_MEMORY_LINEAR .or. &
      in%inputPsiId == INPUT_PSI_DISK_LINEAR) then
     ! Setup the mixing, if necessary -- NEW
     if (in%lin%mixHist_lowaccuracy /= in%lin%mixHist_highaccuracy) then
         ! This must be fixed later
         stop 'in%lin%mixHist_lowaccuracy /= in%lin%mixHist_highaccuracy'
     end if
     select case (in%lin%scf_mode) 
     case (LINEAR_DIRECT_MINIMIZATION)
         ! still do a density mixing, maybe  to be modified later
         if (in%lin%mixHist_lowaccuracy==0) then
             ! simple mixing
             linear_iscf = 12
         else
             ! Pulay mixing
             linear_iscf = 17
         end if
     case (LINEAR_MIXDENS_SIMPLE) 
         if (in%lin%mixHist_lowaccuracy==0) then
             ! simple mixing
             linear_iscf = 12
         else
             ! Pulay mixing
             linear_iscf = 17
         end if
     case (LINEAR_MIXPOT_SIMPLE) 
         if (in%lin%mixHist_lowaccuracy==0) then
             ! simple mixing
             linear_iscf = 2
         else
             ! Pulay mixing
             linear_iscf = 7
         end if
     case (LINEAR_FOE)
         if (in%lin%mixHist_lowaccuracy==0) then
             ! simple mixing
             linear_iscf = 12
         else
             ! Pulay mixing
             linear_iscf = 17
         end if
     case default
         stop 'ERROR: wrong in%lin%scf_mode'
     end select
     call denspot_set_history(denspot,linear_iscf,in%nspin, &
          KSwfn%Lzd%Glr%d%n1i,KSwfn%Lzd%Glr%d%n2i)
     call input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpsp,KSwfn,tmb,energs,&
          inputpsi,input_wf_format,norbv,lzd_old,wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,tmb_old,ref_frags,cdft,&
          locregcenters)
  else
      call input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpsp,KSwfn,tmb,energs,&
           inputpsi,input_wf_format,norbv,lzd_old,wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,tmb_old,ref_frags,cdft)
      nvirt=in%nvirt
      if(in%nvirt > norbv) then
         nvirt = norbv
      end if
  end if
  
  nvirt=in%nvirt
  if(in%nvirt > norbv) then
     nvirt = norbv
  end if

  if (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
      in%inputPsiId == INPUT_PSI_MEMORY_LINEAR .or. &
      in%inputPsiId == INPUT_PSI_DISK_LINEAR) then
      call f_free_ptr(locregcenters)
  end if

  !call deallocate_wfd(wfd_old,subname)
  ! modified by SM
  call deallocate_wfd(wfd_old)
  call deallocate_local_zone_descriptors(lzd_old,subname)

  !save the new grid spacing into the hgrid_old value
  hx_old=KSwfn%Lzd%hgrids(1)
  hy_old=KSwfn%Lzd%hgrids(2)
  hz_old=KSwfn%Lzd%hgrids(3)

  !end of the initialization part
  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  !start the optimization
  energs%eexctX=0.0_gp
  ! Skip the following part in the linear scaling case.
  skip_if_linear: if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
                     .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
     call kswfn_optimization_loop(iproc, nproc, optLoop, &
     & in%alphamix, in%idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, rxyz, GPU, xcstr, &
     & in)
     infocode = optLoop%infocode

     !if we are in the last_run case, validate the last_run only for the last cycle
     !do the last_run things regardless of infocode
     !nrepmax=0 is needed for the Band Structure calculations
     DoLastRunThings=(in%last_run == 1 .and. optLoop%nrepmax == 0) .or. &
          & (in%last_run == 1 .and. optLoop%itrep >= optLoop%nrepmax)
              !print the energies only if they are meaningful
     energy = energs%energy
     !Davidson is set to false first because used in deallocate_before_exiting
     DoDavidson= .false.

     ! Treat the info code from the optimization routine.
     if (infocode == 2 .or. infocode == 3) then
        call deallocate_before_exiting
        call deallocate_bounds(KSwfn%Lzd%Glr%geocode, KSwfn%Lzd%Glr%hybrid_on, KSwfn%lzd%glr%bounds, subname)
        call f_release_routine()
        return
     end if
  else

     ! Allococation of array for Pulay forces (only needed for linear version)
     allocate(fpulay(3,atoms%astruct%nat),stat=i_stat)
     call memocc(i_stat,fpulay,'fpulay',subname)


     call linearScaling(iproc,nproc,KSwfn,tmb,atoms,in,&
          rxyz,denspot,denspot0,nlpsp,GPU,energs,energy,fpulay,infocode,ref_frags,cdft,&
          fdisp, fion)

     ! Clean denspot parts only needed in the SCF loop -- NEW
     call denspot_free_history(denspot)

     ! maybe not the best place to keep it - think about it!
     if (in%lin%calc_transfer_integrals) then
        if (in%lin%constrained_dft) then
           ! switch excess charge to other fragment, recalculate kernel and density and reset lagrange multiplier
           if (iproc==0) write(*,*) '--------------------------------------------------------------------------------------'
           if (iproc==0) write(*,*) 'Warning: site-energy/transfer integral calculation not yet working for constrained DFT'
           if (iproc==0) write(*,*) '--------------------------------------------------------------------------------------'

           !in_frag_charge=f_malloc_ptr(in%frag%nfrag,id='in_frag_charge')
           !call vcopy(in%frag%nfrag,in%frag%charge(1),1,in_frag_charge(1),1)
           !! assume all other fragments neutral, use total system charge to get correct charge for the other fragment
           !in_frag_charge(cdft%ifrag_charged(1))=in%ncharge - in_frag_charge(cdft%ifrag_charged(2))
           !overlap_calculated=.true.
           !call fragment_coeffs_to_kernel(iproc,in%frag,in_frag_charge,ref_frags,tmb,KSwfn%orbs,overlap_calculated)
           !call f_free_ptr(in_frag_charge)
           !cdft%charge=-cdft%charge

           !call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
           !     KSwfn%orbs, tmb, overlap_calculated)     
           !tmb%can_use_transposed=.false. ! - do we really need to deallocate here?
           !i_all = -product(shape(tmb%psit_c))*kind(tmb%psit_c)                               
           !deallocate(tmb%psit_c,stat=i_stat)                                                 
           !call memocc(i_stat,i_all,'tmb%psit_c',subname)                                     
           !i_all = -product(shape(tmb%psit_f))*kind(tmb%psit_f)                               
           !deallocate(tmb%psit_f,stat=i_stat)                                                 
           !call memocc(i_stat,i_all,'tmb%psit_f',subname)     

           !! Now need to calculate the charge density and the potential related to this inputguess
           !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           !     tmb%orbs, tmb%psi, tmb%collcom_sr)

           !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

           !! Must initialize rhopotold (FOR NOW... use the trivial one)
           !call vcopy(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)*in%nspin, &
           !     denspot%rhov(1), 1, denspot0(1), 1)
           !!!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
           !call updatePotential(in%ixc,in%nspin,denspot,energs%eh,energs%exc,energs%evxc)
           !call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

           !! keep a copy of previous wavefunctions and energies...
           !allocate(psi_constrained(tmb%npsidim_orbs), stat=i_stat)
           !call memocc(i_stat, psi_constrained, 'psi_constrained', subname)
           !call vcopy(tmb%npsidim_orbs,tmb%psi(1),1,psi_constrained(1),1)
           !energy_constrained=energy

           !call linearScaling(iproc,nproc,KSwfn,tmb,atoms,in,&
           !     rxyz,denspot,denspot0,nlpsp,GPU,energs,energy,fpulay,infocode,ref_frags,cdft)

           !! calculate matrix elements here...

           !i_all=-product(shape(psi_constrained))*kind(psi_constrained)
           !deallocate(psi_constrained, stat=i_stat)
           !call memocc(i_stat, i_all, 'psi_constrained', subname)

        else
           if (.not. in%lin%fragment_calculation) stop 'Error, fragment calculation needed for transfer integral calculation'
           !if (input%frag%nfrag==2) call calc_transfer_integrals_old(iproc,nproc,input%frag,ref_frags,tmb%orbs,&
           !     tmb%linmat%ham,tmb%linmat%ovrlp)
           call calc_site_energies_transfer_integrals(iproc,nproc,tmb%orthpar%methTransformOverlap,&
                in%frag,ref_frags,tmb%orbs,tmb%linmat%ham,tmb%linmat%ovrlp)
        end if
     end if

     ! deallocate fragments
     if (inputpsi == INPUT_PSI_DISK_LINEAR) then
        if (in%lin%fragment_calculation) then ! we really need to deallocate
           do ifrag=1,in%frag%nfrag_ref
              call fragment_free(ref_frags(ifrag))
           end do
        else ! we haven't actually allocated anything, so can just nullify - should make this more robust/general
           do ifrag=1,in%frag%nfrag_ref
              ref_frags(ifrag)%astruct_frg%nat=-1
              ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
              call fragment_free(ref_frags(ifrag))
              !ref_frags(ifrag)=fragment_null()
           end do
        end if
        deallocate(ref_frags)
     end if

     !!call finalize_p2p_tags()
  
     !temporary allocation of the density
     !!allocate(denspot%rho_work(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)),stat=i_stat)
     !!call memocc(i_stat,denspot%rho_work,'rho',subname)
     !!call vcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),&
     !!     denspot%rhov(1),1,denspot%rho_work(1),1)

     ! keep only the essential part of the density, without the GGA bufffers
     allocate(denspot%rho_work(denspot%dpbox%ndimpot),stat=i_stat)
     call memocc(i_stat,denspot%rho_work,'rho',subname)
     ioffset=kswfn%lzd%glr%d%n1i*kswfn%lzd%glr%d%n2i*denspot%dpbox%i3xcsh
     call vcopy(denspot%dpbox%ndimpot,denspot%rhov(ioffset+1),1,denspot%rho_work(1),1)

     if (infocode==2) then
        !!! Allocate this array since it will be deallcoated in deallocate_before_exiting
        !!allocate(denspot%V_ext(1,1,1,1),stat=i_stat)
        !!call memocc(i_stat,denspot%V_ext,'denspot%V_ext',subname)
        call deallocate_before_exiting()
        i_all=-product(shape(fpulay))*kind(fpulay)
        deallocate(fpulay,stat=i_stat)
        call memocc(i_stat,i_all,'fpulay',subname)
        call destroy_DFT_wavefunction(tmb)
        i_all=-product(shape(KSwfn%psi))*kind(KSwfn%psi)
        deallocate(KSwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'psi',subname)
        call deallocate_wfd(KSwfn%Lzd%Glr%wfd)
        i_all=-product(shape(denspot%rho_work))*kind(denspot%rho_work)
        deallocate(denspot%rho_work,stat=i_stat)
        call memocc(i_stat,i_all,'denspot%rho',subname)
        i_all=-product(shape(KSwfn%orbs%eval))*kind(KSwfn%orbs%eval)
        deallocate(KSwfn%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'KSwfn%orbs%eval',subname)
        call f_release_routine()
        return
     end if

     !infocode = 0
  end if skip_if_linear


  ! allocate KSwfn%psi here instead for case of linear?!
  !if(inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR .or. &
  !                   inputpsi == INPUT_PSI_LINEAR_LCAO) then
  !   allocate(KSwfn%psi(max(KSwfn%orbs%npsidim_comp,KSwfn%orbs%npsidim_orbs)+ndebug),stat=i_stat)
  !   call memocc(i_stat,KSwfn%psi,'psi',subname)
  !end if

  !last run things has to be done:
  !if it is the last run and the infocode is zero
  !if infocode is not zero but the last run has been done for nrepmax times

  DoLastRunThings= (in%last_run == 1 .and. infocode == 0) .or. DoLastRunThings

  !analyse the possibility to calculate Davidson treatment
  !(nvirt > 0 .and. in%inputPsiId == 0)
  DoDavidson= abs(in%norbv) > 0 .and. DoLastRunThings

  !project the wavefunctions on a gaussian basis and keep in memory
  if (in%gaussian_help) then
     call timing(iproc,'gauss_proj','ON') !lr408t
     if (iproc == 0.and.verbose >1) then
        call yaml_comment('Gaussian Basis Projection',hfill='-')
        !write( *,'(1x,a)') '---------------------------------------------------------- Gaussian Basis Projection'
     end if

     !extract the gaussian basis from the pseudowavefunctions
!!!     if (in%inputPsiId == 11) then
!!!        !extract the gaussian basis from the pseudowavefunctions
!!!        call gaussian_pswf_basis(21,.false.,iproc,atoms,rxyz,gbd)
!!!     else if (in%inputPsiId == 12) then
!!!        !extract the gaussian basis from the pseudopotential
!!!        call gaussian_psp_basis(atoms,rxyz,gbd)
!!!     end if

     !extract the gaussian basis from the pseudowavefunctions
     call gaussian_pswf_basis(21,.false.,iproc,in%nspin,atoms,rxyz,KSwfn%gbd,gbd_occ)

     if (associated(gbd_occ)) then
        i_all=-product(shape(gbd_occ))*kind(gbd_occ)
        deallocate(gbd_occ,stat=i_stat)
        call memocc(i_stat,i_all,'gbd_occ',subname)
        nullify(gbd_occ)
     end if


     if (.not. associated(KSwfn%gaucoeffs)) then
        allocate(KSwfn%gaucoeffs(KSwfn%gbd%ncoeff,KSwfn%orbs%nspinor*KSwfn%orbs%norbp+ndebug),stat=i_stat)
        call memocc(i_stat,KSwfn%gaucoeffs,'gaucoeffs',subname)
     end if

     allocate(thetaphi(2,KSwfn%gbd%nat+ndebug),stat=i_stat)
     call memocc(i_stat,thetaphi,'thetaphi',subname)
     thetaphi=0.0_gp

     call wavelets_to_gaussians(atoms%astruct%geocode,KSwfn%orbs%norbp,KSwfn%orbs%nspinor,&
          n1,n2,n3,KSwfn%gbd,thetaphi,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          KSwfn%Lzd%Glr%wfd,KSwfn%psi,KSwfn%gaucoeffs)

     i_all=-product(shape(thetaphi))*kind(thetaphi)
     deallocate(thetaphi,stat=i_stat)
     call memocc(i_stat,i_all,'thetaphi',subname)
     call timing(iproc,'gauss_proj','OF') !lr408t
  end if

  !  write all the wavefunctions into files
  if (in%output_wf_format /= WF_FORMAT_NONE .and. DoLastRunThings) then
     !add flag for writing waves in the gaussian basis form
     !if (in%gaussian_help) then
     if (in%gaussian_help .and. .not.in%inputPsiId==100 .and. .not.in%inputPsiId==101 ) then

!!!        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
!!!
!!!        call gaussian_orthogonality(iproc,nproc,norb,norbp,gbd,gaucoeffs)
        !write the coefficients and the basis on a file
        if (iproc ==0) call yaml_map('Writing wavefunctions in file','wavefunction.gau')
        !if (iproc ==0) write(*,*)'Writing wavefunctions in wavefunction.gau file'
        call write_gaussian_information(iproc,nproc,KSwfn%orbs,KSwfn%gbd,KSwfn%gaucoeffs,trim(in%dir_output) // 'wavefunctions.gau')

        !build dual coefficients
        call dual_gaussian_coefficients(KSwfn%orbs%norbp*KSwfn%orbs%nspinor,KSwfn%gbd,KSwfn%gaucoeffs)

        !control the accuracy of the expansion
        call check_gaussian_expansion(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,KSwfn%psi,KSwfn%gbd,KSwfn%gaucoeffs)

        call deallocate_gwf(KSwfn%gbd,subname)
        i_all=-product(shape(KSwfn%gaucoeffs))*kind(KSwfn%gaucoeffs)
        deallocate(KSwfn%gaucoeffs,stat=i_stat)
        call memocc(i_stat,i_all,'gaucoeffs',subname)
        nullify(KSwfn%gbd%rxyz)

     else
        call writemywaves(iproc,trim(in%dir_output) // "wavefunction", in%output_wf_format, &
             KSwfn%orbs,n1,n2,n3,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             atoms,rxyz,KSwfn%Lzd%Glr%wfd,KSwfn%psi)
     end if
  end if

  write(gridformat, "(A)") ""
  select case (in%output_denspot_format)
  case (output_denspot_FORMAT_ETSF)
     write(gridformat, "(A)") ".etsf"
  case (output_denspot_FORMAT_CUBE)
     write(gridformat, "(A)") ".cube"
  end select

  !plot the ionic potential, if required by output_denspot
  if (in%output_denspot == output_denspot_DENSPOT .and. DoLastRunThings) then
     if (iproc == 0) call yaml_map('Writing external potential in file', 'external_potential'//gridformat)
     !if (iproc == 0) write(*,*) 'writing external_potential' // gridformat
     call plot_density(iproc,nproc,trim(in%dir_output)//'external_potential' // gridformat,&
          atoms,rxyz,denspot%dpbox,1,denspot%V_ext)
  end if
  if (in%output_denspot == output_denspot_DENSPOT .and. DoLastRunThings) then
     if (iproc == 0) call yaml_map('Writing local potential in file','local_potential'//gridformat)
     !if (iproc == 0) write(*,*) 'writing local_potential' // gridformat
     call plot_density(iproc,nproc,trim(in%dir_output)//'local_potential' // gridformat,&
          atoms,rxyz,denspot%dpbox,in%nspin,denspot%rhov)
  end if

  i_all=-product(shape(denspot%V_ext))*kind(denspot%V_ext)
  deallocate(denspot%V_ext,stat=i_stat)
  call memocc(i_stat,i_all,'denspot%V_ext',subname)
  nullify(denspot%V_ext)

  !variables substitution for the PSolver part
  n1i=KSwfn%Lzd%Glr%d%n1i
  n2i=KSwfn%Lzd%Glr%d%n2i
  n3i=KSwfn%Lzd%Glr%d%n3i

  if (inputpsi /= INPUT_PSI_EMPTY) then
     !------------------------------------------------------------------------
     ! here we start the calculation of the forces
     if (iproc == 0) then
        call yaml_comment('Forces Calculation',hfill='-')
        !write( *,'(1x,a)')'----------------------------------------------------------------- Forces Calculation'
     end if

     !refill projectors for tails, davidson
     refill_proj=((in%rbuf > 0.0_gp) .or. DoDavidson) .and. DoLastRunThings

     if (inputpsi /= INPUT_PSI_LINEAR_AO .and. &
          & inputpsi /= INPUT_PSI_MEMORY_LINEAR .and. &
          & inputpsi /= INPUT_PSI_DISK_LINEAR) then
        allocate(fpulay(3,atoms%astruct%nat+ndebug),stat=i_stat)
        call memocc(i_stat,fpulay,'fpulay',subname)
        if (atoms%astruct%nat > 0) call to_zero(3 * atoms%astruct%nat,fpulay(1, 1))
     end if

     if (DoLastRunThings) then
        if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
             .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
           calculate_dipole=.true.
        else
           calculate_dipole = in%lin%calc_dipole
        end if
        output_denspot = in%output_denspot
     else
        output_denspot = -1
        calculate_dipole = .false.
     end if

     call kswfn_post_treatments(iproc, nproc, KSwfn, tmb, &
          & inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR .or. &
          & inputpsi == INPUT_PSI_MEMORY_LINEAR, fxyz, fnoise, fion, fdisp, fpulay, &
          & strten, pressure, ewaldstr, xcstr, GPU, denspot, atoms, rxyz, nlpsp, &
          & output_denspot, in%dir_output, gridformat, refill_proj, calculate_dipole)

     i_all=-product(shape(fpulay))*kind(fpulay)
     deallocate(fpulay,stat=i_stat)
     call memocc(i_stat,i_all,'fpulay',subname)
  end if

  i_all=-product(shape(fion))*kind(fion)
  deallocate(fion,stat=i_stat)
  call memocc(i_stat,i_all,'fion',subname)
  i_all=-product(shape(fdisp))*kind(fdisp)
  deallocate(fdisp,stat=i_stat)
  call memocc(i_stat,i_all,'fdisp',subname)

  !if (nvirt > 0 .and. in%inputPsiId == 0) then
  if (DoDavidson) then

     !for a band structure calculation allocate the array in which to put the eigenvalues
     if (associated(in%kptv) .and. in%nkptv > 0) then
        allocate(band_structure_eval(KSwfn%orbs%norbu+KSwfn%orbs%norbd+in%nspin*norbv,in%nkptv+ndebug),stat=i_stat)
        call memocc(i_stat,band_structure_eval,'band_structure_eval',subname)
     end if

     !proj_G is dummy here, it is only used for PAW
     do iatyp=1,atoms%astruct%ntypes
        call nullify_gaussian_basis(proj_G(iatyp))
     end do

     !calculate Davidson procedure for all the groups of k-points which are chosen
     ikpt=1
     do igroup=1,in%ngroups_kptv

        ! Set-up number of states and shifting values.
        nvirtu = norbv
        nvirtd = 0
        if (in%nspin==2) nvirtd=nvirtu
        ! Create the orbitals.
        if (associated(in%kptv) .and. in%nkptv > 0) then
           nvirtu = nvirtu + KSwfn%orbs%norbu
           nvirtd = nvirtd + KSwfn%orbs%norbd
           nvirt  = nvirtu+nvirtd

           !number of k-points for this group
           nkptv = in%nkptsv_group(igroup) !size(in%kptv, 2)

           allocate(wkptv(nkptv+ndebug),stat=i_stat)
           call memocc(i_stat,wkptv,'wkptv',subname)
           wkptv(:) = real(1.0, gp) / real(nkptv, gp)
           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                KSwfn%orbs%nspin,KSwfn%orbs%nspinor,nkptv, &
                in%kptv(:,sum(in%nkptsv_group(1:igroup - 1)) + 1:sum(in%nkptsv_group(1:igroup))), &
                wkptv,VTwfn%orbs,.false.)
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,KSwfn%Lzd%Glr,VTwfn%orbs,VTwfn%comms)  

           i_all=-product(shape(wkptv))*kind(wkptv)
           deallocate(wkptv,stat=i_stat)
           call memocc(i_stat,i_all,'wkptv',subname)

           !free projectors
           call free_DFT_PSP_projectors(nlpsp)

           ! Calculate all projectors, or allocate array for on-the-fly calculation
           call timing(iproc,'CrtProjectors ','ON')
           call createProjectorsArrays(KSwfn%Lzd%Glr,rxyz,atoms,VTwfn%orbs,&
                radii_cf,in%frmult,in%frmult,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                .false.,nlpsp,proj_G) 
           call timing(iproc,'CrtProjectors ','OF') 
           if (iproc == 0) call print_nlpsp(nlpsp)

        else
           !the virtual orbitals should be in agreement with the traditional k-points
           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                KSwfn%orbs%nspin,KSwfn%orbs%nspinor,KSwfn%orbs%nkpts,&
                KSwfn%orbs%kpts,KSwfn%orbs%kwgts,VTwfn%orbs,.false.,basedist=KSwfn%orbs%norb_par(0:,1:))
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,KSwfn%Lzd%Glr,VTwfn%orbs,VTwfn%comms,&
                basedist=KSwfn%comms%nvctr_par(0:,1:))  

        end if

        !allocate psivirt pointer (note the orbs dimension)
        allocate(VTwfn%psi(max(VTwfn%orbs%npsidim_comp,VTwfn%orbs%npsidim_orbs)+ndebug),stat=i_stat)
        call memocc(i_stat,VTwfn%psi,'psivirt',subname)
        !to avoid problems with the bindings
        VTwfn%c_obj=0

        !define Local zone descriptors
        VTwfn%Lzd = KSwfn%Lzd
        VTwfn%orthpar=KSwfn%orthpar
        allocate(VTwfn%confdatarr(VTwfn%orbs%norbp))
        call default_confinement_data(VTwfn%confdatarr,VTwfn%orbs%norbp)


        if (in%norbv < 0) then
           call direct_minimization(iproc,nproc,in,atoms,& 
                nvirt,rxyz,denspot%rhov,nlpsp, &
                denspot%pkernelseq,denspot%dpbox,denspot%xc,GPU,KSwfn,VTwfn)

           if(abs(in%nplot)>KSwfn%orbs%norb+nvirt) then
              if(iproc==0) call yaml_warning('More plots requested than orbitals calculated')
           end if
        else if (in%norbv > 0) then
           call davidson(iproc,nproc,in,atoms,& 
                KSwfn%orbs,VTwfn%orbs,in%nvirt,VTwfn%Lzd,&
                KSwfn%comms,VTwfn%comms,&
                rxyz,denspot%rhov,nlpsp, &
                denspot%pkernelseq,KSwfn%psi,VTwfn%psi,denspot%dpbox,denspot%xc,GPU)
!!$           call constrained_davidson(iproc,nproc,in,atoms,&
!!$                orbs,orbsv,in%nvirt,Lzd%Glr,comms,VTwfn%comms,&
!!$                hx,hy,hz,rxyz,denspot%rhov,nlpsp, &
!!$                psi,VTwfn%psi,nscatterarr,ngatherarr,GPU)
           if(abs(in%nplot)>KSwfn%orbs%norb+in%nvirt) then
              if(iproc==0) call yaml_warning('More plots requested than orbitals calculated')
           end if
        end if
        if(in%output_wf_format == 2 .and. abs(in%norbv)>0 ) then
           call dump_eigenfunctions(trim(in%dir_output),in%nplot,atoms,VTwfn%Lzd%hgrids,VTwfn%Lzd%Glr,&
                KSwfn%orbs,VTwfn%orbs,rxyz,KSwfn%psi,VTwfn%psi)
        end if

        deallocate(VTwfn%confdatarr)

        ! Write virtual wavefunctions in ETSF format: WORKS ONLY FOR ONE KPOINT 
        if(in%output_wf_format == 3 .and. abs(in%norbv) > 0) then
           write(wfformat, "(A)") ""
           select case (in%output_wf_format)
           case (WF_FORMAT_ETSF)
              write(wfformat, "(A)") ".etsf"
           case (WF_FORMAT_BINARY)
              write(wfformat, "(A)") ".bin"
           end select

           call  writemywaves(iproc,trim(in%dir_output) // "virtuals" // trim(wfformat),&
                in%output_wf_format, &
                VTwfn%orbs,n1,n2,n3,&
                KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                atoms,rxyz,KSwfn%Lzd%Glr%wfd,VTwfn%psi)
        end if

        ! Write virtual wavefunctions in ETSF format
        if (in%output_wf_format /= WF_FORMAT_NONE  .and. abs(in%norbv) > 0) then
           call  writemywaves(iproc,trim(in%dir_output) // "virtuals", in%output_wf_format, &
                VTwfn%orbs,n1,n2,n3,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                atoms,rxyz,KSwfn%Lzd%Glr%wfd,VTwfn%psi)
        end if

        !start the Casida's treatment 
        if (in%tddft_approach=='TDA') then

           !does it makes sense to use GPU only for a one-shot sumrho?
           if (GPU%OCLconv) then
              call allocate_data_OCL(KSwfn%Lzd%Glr%d%n1,KSwfn%Lzd%Glr%d%n2,KSwfn%Lzd%Glr%d%n3,&
                   atoms%astruct%geocode,&
                   in%nspin,KSwfn%Lzd%Glr%wfd,KSwfn%orbs,GPU)
           end if

           !this could have been calculated before
           ! Potential from electronic charge density
           !WARNING: this is good just because the TDDFT is done with LDA
           call sumrho(denspot%dpbox,KSwfn%orbs,KSwfn%Lzd,&
                GPU,atoms%astruct%sym,denspot%rhod,denspot%xc,KSwfn%psi,denspot%rho_psi)
           call communicate_density(denspot%dpbox,KSwfn%orbs%nspin,&
                denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)
           call denspot_set_rhov_status(denspot, ELECTRONIC_DENSITY, -1,iproc,nproc)

           if (GPU%OCLconv) then
              call free_gpu_OCL(GPU,KSwfn%orbs,in%nspin)
           end if

           !Allocate second Exc derivative
           if (denspot%dpbox%n3p >0) then
              allocate(denspot%f_XC(n1i,n2i,denspot%dpbox%n3p,in%nspin+1+ndebug),stat=i_stat)
              call memocc(i_stat,denspot%f_XC,'f_XC',subname)
           else
              allocate(denspot%f_XC(1,1,1,in%nspin+1+ndebug),stat=i_stat)
              call memocc(i_stat,denspot%f_XC,'denspot%f_XC',subname)
           end if

           call XC_potential(atoms%astruct%geocode,'D',iproc,nproc,bigdft_mpi%mpi_comm,&
                KSwfn%Lzd%Glr%d%n1i,KSwfn%Lzd%Glr%d%n2i,KSwfn%Lzd%Glr%d%n3i,denspot%xc,&
                denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
                denspot%rhov,energs%exc,energs%evxc,in%nspin,denspot%rho_C,denspot%V_XC,xcstr,denspot%f_XC)
           call denspot_set_rhov_status(denspot, CHARGE_DENSITY, -1,iproc,nproc)

           !select the active space if needed

           call tddft_casida(iproc,nproc,atoms,rxyz,&
                denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
                denspot%dpbox%n3p,denspot%dpbox%ngatherarr(0,1),&
                KSwfn%Lzd%Glr,KSwfn%orbs,VTwfn%orbs,denspot%dpbox%i3s+denspot%dpbox%i3xcsh,&
                denspot%f_XC,denspot%pkernelseq,KSwfn%psi,VTwfn%psi)

           i_all=-product(shape(denspot%f_XC))*kind(denspot%f_XC)
           deallocate(denspot%f_XC,stat=i_stat)
           call memocc(i_stat,i_all,'denspot%f_XC',subname)

        end if

        call deallocate_comms(VTwfn%comms,subname)
        call deallocate_orbs(VTwfn%orbs,subname)

        !in the case of band structure calculation, copy the values of the eigenvectors
        !into a new array to write them afterwards
        if (associated(in%kptv) .and. in%nkptv > 0) then
           call vcopy(VTwfn%orbs%norb*nkptv,VTwfn%orbs%eval(1),1,band_structure_eval(1,ikpt),1)
           !increment the value of ikpt
           ikpt=ikpt+in%nkptsv_group(igroup)
        end if

        i_all=-product(shape(VTwfn%orbs%eval))*kind(VTwfn%orbs%eval)
        deallocate(VTwfn%orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'eval',subname)

        !if the local analysis has to be performed the deallocation should not be done
        i_all=-product(shape(VTwfn%psi))*kind(VTwfn%psi)
        deallocate(VTwfn%psi,stat=i_stat)
        call memocc(i_stat,i_all,'VTwfn%psi',subname)

     end do

     if (associated(in%kptv) .and. in%nkptv > 0) then
        !dump the band structure eigenvalue on a file and deallocate it
        if (iproc == 0) then
           open(unit=11,file='band_structure.dat',status='unknown')
           do ikpt=1,in%nkptv
              write(11,'(i5,3(f12.6),10000(1pe12.4))')ikpt,&
                   (in%kptv(i,ikpt),i=1,3),(band_structure_eval(i,ikpt),i=1,VTwfn%orbs%norb)
           end do
           !tentative gnuplot string for the band structure file
           write(11,'(a,9999(a,i6,a))')&
                "#plot 'band_structure.dat' u 1:5 w l t ''",&
                (",'' u 1:",5+i-1," w l t ''" ,i=2,VTwfn%orbs%norb)
           close(unit=11)
        end if
        i_all=-product(shape(band_structure_eval))*kind(band_structure_eval)
        deallocate(band_structure_eval,stat=i_stat)
        call memocc(i_stat,i_all,'band_structure_eval',subname)
     end if

  end if


  !perform here the mulliken charge and density of states
  !localise them on the basis of gatom of a number of atoms
  !if (in%gaussian_help .and. DoLastRunThings) then
  if (in%gaussian_help .and. DoLastRunThings .and.&
&    (.not.inputpsi==INPUT_PSI_LINEAR_AO .and. .not.inputpsi==INPUT_PSI_DISK_LINEAR &
      .and. .not. inputpsi==INPUT_PSI_MEMORY_LINEAR)) then
     !here one must check if psivirt should have been kept allocated
     if (.not. DoDavidson) then
        VTwfn%orbs%norb=0
        VTwfn%orbs%norbp=0
     end if
     call local_analysis(iproc,nproc,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          atoms,rxyz,KSwfn%Lzd%Glr,KSwfn%orbs,VTwfn%orbs,KSwfn%psi,VTwfn%psi)
  else if (DoLastRunThings .and. optLoop%itrpmax /= 1 .and. verbose >= 2) then
     ! Do a full DOS calculation.
     if (iproc == 0) call global_analysis(KSwfn%orbs, in%Tel,in%occopt)
  end if

!!$  i_all=-product(shape(denspot%pkernel))*kind(denspot%pkernel)
!!$  deallocate(denspot%pkernel,stat=i_stat)
!!$  call memocc(i_stat,i_all,'kernel',subname)

  if (((in%exctxpar == 'OP2P' .and. xc_exctXfac(denspot%xc) /= 0.0_gp) &
       .or. in%SIC%alpha /= 0.0_gp) .and. nproc >1) then
     
     !if (loc(denspot%pkernelseq%kernel) /= loc(denspot%pkernel%kernel)) then !this is not standard
     if (.not. associated(denspot%pkernelseq%kernel,target=denspot%pkernel%kernel) .and. &
          associated(denspot%pkernelseq%kernel)) then
        call pkernel_free(denspot%pkernelseq,subname)
     end if
!!$     i_all=-product(shape(denspot%pkernelseq))*kind(denspot%pkernelseq)
!!$     deallocate(denspot%pkernelseq,stat=i_stat)
!!$     call memocc(i_stat,i_all,'kernelseq',subname)
  else if (nproc == 1 .and. (in%exctxpar == 'OP2P' .or. in%SIC%alpha /= 0.0_gp)) then
     nullify(denspot%pkernelseq%kernel)
  end if
  call pkernel_free(denspot%pkernel,subname)


  !------------------------------------------------------------------------
  if ((in%rbuf > 0.0_gp) .and. atoms%astruct%geocode == 'F' .and. DoLastRunThings ) then
     if (in%SIC%alpha /= 0.0_gp) then
        if (iproc==0) call yaml_warning('Tail correction not admitted with SIC corrections for the moment')
        !write(*,*)&
        !     &   'ERROR: Tail correction not admitted with SIC corrections for the moment'
        stop
     end if
     call timing(iproc,'Tail          ','ON')
     !    Calculate energy correction due to finite size effects
     !    ---reformat potential
     allocate(denspot%pot_work(n1i*n2i*n3i*in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%pot_work,'denspot%pot_work',subname)

     if (nproc > 1) then
        call MPI_ALLGATHERV(denspot%rhov,n1i*n2i*denspot%dpbox%n3p,&
             mpidtypd,denspot%pot_work(1),denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2), & 
             mpidtypd,denspot%dpbox%mpi_env%mpi_comm,ierr)
        !print '(a,2f12.6)','RHOup',sum(abs(rhopot(:,:,:,1))),sum(abs(pot(:,:,:,1)))
        if(in%nspin==2) then
           !print '(a,2f12.6)','RHOdw',sum(abs(rhopot(:,:,:,2))),sum(abs(pot(:,:,:,2)))
           call MPI_ALLGATHERV(denspot%rhov(1+n1i*n2i*denspot%dpbox%n3p),n1i*n2i*denspot%dpbox%n3p,&
                mpidtypd,denspot%pot_work(1+n1i*n2i*n3i),&
                denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2), & 
                mpidtypd,denspot%dpbox%mpi_env%mpi_comm,ierr)
        end if
     else
        call vcopy(n1i*n2i*n3i*in%nspin,denspot%rhov(1),1,denspot%pot_work(1),1)
     end if

     call dpbox_free(denspot%dpbox, subname)

     i_all=-product(shape(denspot%rhov))*kind(denspot%rhov)
     deallocate(denspot%rhov,stat=i_stat)
     call memocc(i_stat,i_all,'denspot%rhov',subname)

     i_all=-product(shape(denspot%V_XC))*kind(denspot%V_XC)
     deallocate(denspot%V_XC,stat=i_stat)
     call memocc(i_stat,i_all,'denspot%V_XC',subname)

     !pass hx instead of hgrid since we are only in free BC
     call CalculateTailCorrection(iproc,nproc,atoms,in%rbuf,KSwfn%orbs,&
          KSwfn%Lzd%Glr,nlpsp,in%ncongt,denspot%pot_work,KSwfn%Lzd%hgrids(1),&
          rxyz,radii_cf,in%crmult,in%frmult,in%nspin,&
          KSwfn%psi,(in%output_denspot /= 0),energs%ekin,energs%epot,energs%eproj)

     i_all=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
     deallocate(denspot%pot_work,stat=i_stat)
     call memocc(i_stat,i_all,'denspot%pot_work',subname)

     energs%ebs=energs%ekin+energs%epot+energs%eproj
     energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%evsic+energs%eion+energs%edisp-energs%eTS+energs%ePV

     if (iproc == 0) then
        call yaml_open_map('Corrected Energies', flow=.true.)
        call yaml_map('Ekin', energs%ekin, fmt='(1pe18.11)')
        call yaml_map('Epot', energs%epot, fmt='(1pe18.11)')
        call yaml_map('Eproj',energs%eproj,fmt='(1pe18.11)')
        call yaml_close_map()
        call yaml_map('Total energy with tail correction',energy,fmt='(1pe24.17)')
        call yaml_close_map()
        !write( *,'(1x,a,3(1x,1pe18.11))')&
        !     &   '  Corrected ekin,epot,eproj',energs%ekin,energs%epot,energs%eproj
        !write( *,'(1x,a,1x,1pe24.17)')&
        !     &   'Total energy with tail correction',energy
     endif

     call timing(iproc,'Tail          ','OF')
  else
     !    No tail calculation
     if (nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
     i_all=-product(shape(denspot%rhov))*kind(denspot%rhov)
     deallocate(denspot%rhov,stat=i_stat)
     call memocc(i_stat,i_all,'denspot%rhov',subname)
     i_all=-product(shape(denspot%V_XC))*kind(denspot%V_XC)
     deallocate(denspot%V_XC,stat=i_stat)
     call memocc(i_stat,i_all,'denspot%V_XC',subname)
     call dpbox_free(denspot%dpbox, subname)
  endif
  ! --- End if of tail calculation

  !?!   !Finally, we add the entropic contribution to the energy from non-integer occnums
  !?!   if(orbs%eTS>0_gp) then 
  !?!      energy=energy - orbs%eTS 
  !?! 
  !?!      if (iproc == 0) then
  !?!         write( *,'(1x,a,1(1x,1pe18.11))')&
  !?!              '  Entropic correction due to electronic tempertature',orbs%eTS
  !?!         write( *,'(1x,a,1x,1pe24.17)')&
  !?!              'Free energy (= total energy - T*S)  ',energy
  !?!      endif
  !?!    endif

  call deallocate_before_exiting

  call f_release_routine()
contains

  !> Routine which deallocate the pointers and the arrays before exiting 
  subroutine deallocate_before_exiting
    
  !when this condition is verified we are in the middle of the SCF cycle
    if (infocode /=0 .and. infocode /=1 .and. inputpsi /= INPUT_PSI_EMPTY) then
       i_all=-product(shape(denspot%V_ext))*kind(denspot%V_ext)
       deallocate(denspot%V_ext,stat=i_stat)
       call memocc(i_stat,i_all,'denspot%V_ext',subname)

       if (((in%exctxpar == 'OP2P' .and. xc_exctXfac(denspot%xc) /= 0.0_gp) &
            .or. in%SIC%alpha /= 0.0_gp) .and. nproc >1) then
!          if (loc(denspot%pkernelseq%kernel) /= loc(denspot%pkernel%kernel)) then !not standard
             if (.not. associated(denspot%pkernelseq%kernel,target=denspot%pkernel%kernel) .and. &
                  associated(denspot%pkernelseq%kernel)) then
             call pkernel_free(denspot%pkernelseq,subname)
          end if
       else if (nproc == 1 .and. (in%exctxpar == 'OP2P' .or. in%SIC%alpha /= 0.0_gp)) then
          nullify(denspot%pkernelseq%kernel)
       end if
       call pkernel_free(denspot%pkernel,subname)
!!$       i_all=-product(shape(denspot%pkernel))*kind(denspot%pkernel)
!!$       deallocate(denspot%pkernel,stat=i_stat)
!!$       call memocc(i_stat,i_all,'kernel',subname)

       ! calc_tail false
       i_all=-product(shape(denspot%rhov))*kind(denspot%rhov)
       deallocate(denspot%rhov,stat=i_stat)
       call memocc(i_stat,i_all,'denspot%rhov',subname)
       i_all=-product(shape(denspot%V_XC))*kind(denspot%V_XC)
       deallocate(denspot%V_XC,stat=i_stat)
       call memocc(i_stat,i_all,'denspot%V_XC',subname)

       call dpbox_free(denspot%dpbox, subname)

       i_all=-product(shape(fion))*kind(fion)
       deallocate(fion,stat=i_stat)
       call memocc(i_stat,i_all,'fion',subname)
       i_all=-product(shape(fdisp))*kind(fdisp)
       deallocate(fdisp,stat=i_stat)
       call memocc(i_stat,i_all,'fdisp',subname)
    end if
    call xc_end(denspot%xc)

    !free GPU if it is the case
    if (GPUconv .and. .not.(DoDavidson)) then
       call free_gpu(GPU,KSwfn%orbs%norbp)
    else if (GPU%OCLconv .and. .not.(DoDavidson)) then
       call free_gpu_OCL(GPU,KSwfn%orbs,in%nspin)
    end if

    ! Free all remaining parts of denspot
    call deallocate_rho_descriptors(denspot%rhod,subname)
    if(associated(denspot%rho_C)) then
       i_all=-product(shape(denspot%rho_C))*kind(denspot%rho_C)
       deallocate(denspot%rho_C,stat=i_stat)
       call memocc(i_stat,i_all,'denspot%rho_C',subname)
    end if
    i_all=-product(shape(denspot0))*kind(denspot0)
    deallocate(denspot0, stat=i_stat)
    call memocc(i_stat, i_all, 'denspot0', subname)

    ! Free all remaining parts of KSwfn
!!write(*,*) 'WARNING HERE!!!!!'
    !if(inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
    !                 .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
    if (in%inguess_geopt/=1) then
        call deallocate_bounds(KSwfn%Lzd%Glr%geocode,KSwfn%Lzd%Glr%hybrid_on,&
             KSwfn%Lzd%Glr%bounds,subname)
    end if
    call deallocate_Lzd_except_Glr(KSwfn%Lzd, subname)

!    i_all=-product(shape(KSwfn%Lzd%Glr%projflg))*kind(KSwfn%Lzd%Glr%projflg)
!    deallocate(KSwfn%Lzd%Glr%projflg,stat=i_stat)
!    call memocc(i_stat,i_all,'Glr%projflg',subname)
    call deallocate_comms(KSwfn%comms,subname)
    call deallocate_orbs(KSwfn%orbs,subname)
    if (inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
        .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
       deallocate(KSwfn%confdatarr)
    else
       !deallocate(tmb%confdatarr)
    end if

    ! Free projectors.
    call free_DFT_PSP_projectors(nlpsp)
!!$    call deallocate_proj_descr(nlpsp,subname)
!!$    i_all=-product(shape(proj))*kind(proj)
!!$    deallocate(proj,stat=i_stat)
!!$    call memocc(i_stat,i_all,'proj',subname)

    !end of wavefunction minimisation
    call timing(bigdft_mpi%mpi_comm,'LAST','PR')
    call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm)    
    call cpu_time(tcpu1)
    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    if (iproc == 0) then
       call yaml_comment('Timing for root process',hfill='-')
       call yaml_open_map('Timings for root process')
       call yaml_map('CPU time (s)',tcpu1-tcpu0,fmt='(f12.2)')
       call yaml_map('Elapsed time (s)',tel,fmt='(f12.2)')
       call yaml_close_map()
    end if
!       &
!         &   write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0

    ! Stop signals
    if (in%signaling .and. iproc == 0) then
       call bigdft_signals_rm_denspot(in%gmainloop)
       call bigdft_signals_rm_energs(in%gmainloop)
       call bigdft_signals_rm_wf(in%gmainloop)
       call bigdft_signals_rm_optloop(in%gmainloop)
       call localfields_free_wrapper(denspot%c_obj)
       call energs_free_wrapper(energs%c_obj)
       call optloop_free_wrapper(optLoop%c_obj)
       call wf_free_wrapper(KSwfn%c_obj)
       call wf_free_wrapper(tmb%c_obj)
    end if

!!$    if(inputpsi ==  INPUT_PSI_LINEAR) then
!!$        i_all=-product(shape(atoms%rloc))*kind(atoms%rloc)
!!$        deallocate(atoms%rloc,stat=i_stat)
!!$        call memocc(i_stat,i_all,'atoms%rloc',subname)
!!$    end if

     if (iproc == 0 .and. (in%inputPsiId==1 .or. in%inputPsiId==0) .and. infocode==1) then
        call yaml_warning('Self-consistent cycle did not meet convergence criteria')
     end if
    !release the yaml document
    call yaml_release_document()

  END SUBROUTINE deallocate_before_exiting

END SUBROUTINE cluster


!> Kohn-Sham wavefunction optimization loop
subroutine kswfn_optimization_loop(iproc, nproc, opt, &
     & alphamix, idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, rxyz, GPU, xcstr, &
     & in)
  use module_base
  use module_types
  use module_interfaces, except_this_one => kswfn_optimization_loop
  use yaml_output
  implicit none
  real(dp), dimension(6), intent(out) :: xcstr
  integer, intent(in) :: iproc, nproc, idsx, inputpsi
  real(gp), intent(in) :: alphamix
  type(DFT_optimization_loop), intent(inout) :: opt
  type(DFT_wavefunction), intent(inout) :: KSwfn
  type(DFT_local_fields), intent(inout) :: denspot
  type(energy_terms), intent(inout) :: energs
  type(atoms_data), intent(in) :: atoms
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(input_variables), intent(in) :: in !<todo: Remove me

  character(len = *), parameter :: subname = "kswfn_optimization_loop"
  logical :: endloop, scpot, endlooprp, lcs
  integer :: ndiis_sd_sw, idsx_actual_before, linflag, ierr,iter_for_diis
  real(gp) :: gnrm_zero
  character(len=5) :: final_out
  !temporary variables for PAPI computation
  ! real(kind=4) :: rtime, ptime,  mflops
  ! integer(kind=8) ::flpops

!  !start PAPI counting
!  if (iproc==0) call PAPIF_flops(rtime, ptime, flpops, mflops,ierr)

  ! Setup the mixing, if necessary
  call denspot_set_history(denspot,opt%iscf,in%nspin, &
       KSwfn%Lzd%Glr%d%n1i,KSwfn%Lzd%Glr%d%n2i)

  ! allocate arrays necessary for DIIS convergence acceleration
  call allocate_diis_objects(idsx,in%alphadiis,sum(KSwfn%comms%ncntt(0:nproc-1)),&
       KSwfn%orbs%nkptsp,KSwfn%orbs%nspinor,KSwfn%diis,subname)

  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=KSwfn%diis%idsx

  gnrm_zero=0.0d0
  opt%gnrm=1.d10
  opt%rpnrm=1.d10
  endlooprp=.false.
  energs%e_prev=0.0_gp

  !normal opt%infocode, if everything go through smoothly we should keep this
  opt%infocode=0
  !yaml output
  if (iproc==0) then
     call yaml_comment('Self-Consistent Cycle',hfill='-')
     call yaml_open_sequence('Ground State Optimization')
  end if
  opt%itrp=1
  rhopot_loop: do
  KSwfn%diis%energy_old=1.d100
     if (opt%itrp > opt%itrpmax) exit
     !yaml output 
     if (iproc==0) then
        call yaml_sequence(advance='no')
        call yaml_open_sequence("Hamiltonian Optimization",label=&
             'itrp'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i3.3)'))))

     end if
     !set the opt%infocode to the value it would have in the case of no convergence
     opt%infocode=1
     opt%itrep=1
     iter_for_diis=0 !initialize it here for keeping the history also after a subspace diagonalization
     subd_loop: do
        if (opt%itrep > opt%nrepmax) exit subd_loop
        !yaml output 
        if (iproc==0) then
           call yaml_sequence(advance='no')
           call yaml_open_map("Subspace Optimization",label=&
                'itrep'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i3.3)')))//'-'//&
                trim(adjustl(yaml_toa(opt%itrep,fmt='(i2.2)'))))
        end if

        !yaml output
        if (iproc==0) then
           call yaml_open_sequence("Wavefunctions Iterations")
        end if
        opt%iter=1
        iter_for_diis=0
        wfn_loop: do
           if (opt%iter > opt%itermax) exit wfn_loop

           !control whether the minimisation iterations should end after the hamiltionian application
           endloop= opt%gnrm <= opt%gnrm_cv .or. opt%iter == opt%itermax

           if (iproc == 0) then 
              !yaml output
              if (endloop .and. (opt%itrpmax==1 .or. opt%itrpmax >1 .and. endlooprp)) then
                 call yaml_sequence(label='FINAL'//trim(adjustl(yaml_toa(opt%itrep,fmt='(i3.3)'))),advance='no')
              else if (endloop .and. opt%itrep == opt%nrepmax) then
                 call yaml_sequence(label='final'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))),&
                      advance='no')
              else
                 call yaml_sequence(advance='no')
              end if
              call yaml_open_map(flow=.true.)
              if (verbose > 0) &
                   call yaml_comment('iter:'//yaml_toa(opt%iter,fmt='(i6)'),hfill='-')
           endif

           !control how many times the DIIS has switched into SD
           if (KSwfn%diis%idsx /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

           !let SD runs if the DIIS did not work the second time
           if (ndiis_sd_sw > 1) then
              KSwfn%diis%switchSD=.false.
           end if

           !stop the partial timing counter if necessary
           if (endloop .and. opt%itrpmax==1) call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')
           !logical flag for the self-consistent potential
           scpot=(opt%iscf > SCF_KIND_DIRECT_MINIMIZATION .and. opt%iter==1 .and. opt%itrep==1) .or. & !mixing to be done
                (opt%iscf <= SCF_KIND_DIRECT_MINIMIZATION)!direct minimisation
           !allocate the potential in the full box
           !temporary, should change the use of flag in full_local_potential2
           linflag = 1                                 
           if(in%linear == INPUT_IG_OFF) linflag = 0
           if(in%linear == INPUT_IG_TMO) linflag = 2
           call psitohpsi(iproc,nproc,atoms,scpot,denspot,opt%itrp,opt%iter,opt%iscf,alphamix,&
                nlpsp,rxyz,linflag,in%unblock_comms,GPU,KSwfn,energs,opt%rpnrm,xcstr)

           endlooprp= (opt%itrp > 1 .and. opt%rpnrm <= opt%rpnrm_cv) .or. opt%itrp == opt%itrpmax

           call total_energies(energs, opt%iter, iproc)

           !check for convergence or whether max. numb. of iterations exceeded
           if (endloop) then
              if (opt%gnrm < opt%gnrm_cv) opt%infocode=0
              exit wfn_loop 
           endif

           !evaluate the functional of the wavefunctions and put it into the diis structure
           !the energy values is printed out in this routine
           call calculate_energy_and_gradient(opt%iter,iproc,nproc,GPU,in%ncong,opt%iscf,&
                energs,KSwfn,opt%gnrm,gnrm_zero)

           !control the previous value of idsx_actual
           idsx_actual_before=KSwfn%diis%idsx
           iter_for_diis=iter_for_diis+1
           call hpsitopsi(iproc,nproc,iter_for_diis,idsx,KSwfn,atoms,nlpsp)

           if (inputpsi == INPUT_PSI_LCAO) then
              if ((opt%gnrm > 4.d0 .and. KSwfn%orbs%norbu /= KSwfn%orbs%norbd) .or. &
                   &   (KSwfn%orbs%norbu == KSwfn%orbs%norbd .and. opt%gnrm > 10.d0)) then
                 opt%infocode=3
              end if
           else if (inputpsi == INPUT_PSI_MEMORY_WVL) then
              if (opt%gnrm > 1.d0) then
                 opt%infocode=2
              end if
           end if
           !flush all writings on standart output
           if (iproc==0) then
              !yaml output
              call yaml_close_map() !iteration
              call bigdft_utils_flush(unit=6)
           end if
           ! Emergency exit case
           if (opt%infocode == 2 .or. opt%infocode == 3) then
              if (nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
              !>todo: change this return into a clean out of the routine, so the YAML is clean.
              if (iproc==0) then
                 !call yaml_close_map()
                 call yaml_close_sequence() !wfn iterations
                 call yaml_close_map()
                 call yaml_close_sequence() !itrep
                 if (opt%infocode==2) then
                    call yaml_warning('The norm of the residue is too large, need to recalculate input wavefunctions')
                 else if (opt%infocode ==3) then
                    call yaml_warning('The norm of the residue is too large also with input wavefunctions.')
                 end if
              end if

              exit rhopot_loop
              !return
           end if

           if (opt%c_obj /= 0) then
              call optloop_emit_iter(opt, OPTLOOP_WAVEFUNCTIONS, energs, iproc, nproc)
           end if

           opt%iter = opt%iter + 1
        end do wfn_loop


        if (opt%c_obj /= 0) then
           call optloop_emit_done(opt, OPTLOOP_WAVEFUNCTIONS, energs, iproc, nproc)
        end if

        if (iproc == 0) then
           !if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',opt%iter,' minimization iterations required'
           !write( *,'(1x,a)') &
           !     &   '--------------------------------------------------- End of Wavefunction Optimisation'
           if ((opt%itrpmax >1 .and. endlooprp) .or. opt%itrpmax == 1) then
              write(final_out, "(A5)") "FINAL"
           else
              write(final_out, "(A5)") "final"
           end if
           call write_energies(opt%iter,0,energs,opt%gnrm,gnrm_zero,final_out)
           call yaml_close_map()

           if (opt%itrpmax >1) then
              if ( KSwfn%diis%energy > KSwfn%diis%energy_min) &
                   call yaml_warning('Found an energy value lower than the ' // final_out // &
                   ' energy, delta:' // trim(yaml_toa(KSwfn%diis%energy-KSwfn%diis%energy_min,fmt='(1pe9.2)')))
           else
              !write this warning only if the system is closed shell
              call check_closed_shell(KSwfn%orbs,lcs)
              if (lcs) then
                 if ( energs%eKS > KSwfn%diis%energy_min) &
                      call yaml_warning('Found an energy value lower than the FINAL energy, delta:' // &
                      trim(yaml_toa(energs%eKS-KSwfn%diis%energy_min,fmt='(1pe9.2)')))
              end if
           end if
        end if

        if (iproc==0) then
           call yaml_close_sequence() !wfn iterations
           if (opt%iter == opt%itermax .and. opt%infocode/=0) &
             &  call yaml_comment('No convergence within the allowed number of minimization steps')
             !&   write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'
        end if
        call last_orthon(iproc,nproc,opt%iter,KSwfn,energs%evsum,.true.) !never deallocate psit and hpsi

!!$        !EXPERIMENTAL
!!$        !check if after convergence the integral equation associated with Helmholtz' Green function is satisfied
!!$        !note: valid only for negative-energy eigenstates
!!$        call integral_equation(iproc,nproc,atoms,KSwfn,denspot%dpbox%ngatherarr,denspot%rhov,GPU,proj,nlpspd,rxyz)

        !exit if the opt%infocode is correct
        if (opt%infocode /= 0) then
           if(iproc==0) then
              if (opt%itrp == opt%itrpmax .and. opt%gnrm_cv > 0.0_gp) &
                   call yaml_warning('Wavefunctions not converged after cycle '// trim(yaml_toa(opt%itrep,fmt='(i0)')))
              if (opt%itrep < opt%nrepmax) call yaml_comment('restart after diagonalisation')
              ! write(*,*) ' WARNING: Wavefunctions not converged after cycle',opt%itrep
              ! if (opt%itrep < opt%nrepmax) write(*,*)' restart after diagonalisation'
           end if
           opt%gnrm=1.d10

           if (opt%itrpmax == 1 .and. in%norbsempty > 0) then
              !recalculate orbitals occupation numbers
              call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)
              
              !opt%gnrm =1.d10
              KSwfn%diis%energy_min=1.d10
              !KSwfn%diis%alpha=2.d0
              KSwfn%diis%alpha=in%alphadiis
           end if
        end if

        if (opt%itrpmax ==1) then 
           call eigensystem_info(iproc,nproc,opt%gnrm,&
             KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f,&
             KSwfn%orbs,KSwfn%psi)
           if (opt%infocode /=0) then
              opt%gnrm =1.d10
           end if
        end if

        if (iproc==0) then
           call yaml_close_map()
        end if

        if (opt%infocode ==0) exit subd_loop

        if (opt%c_obj /= 0) then
           call optloop_emit_iter(opt, OPTLOOP_SUBSPACE, energs, iproc, nproc)
        end if
        
        opt%itrep = opt%itrep + 1
     end do subd_loop
     if (opt%c_obj /= 0) then
        call optloop_emit_done(opt, OPTLOOP_SUBSPACE, energs, iproc, nproc)
     end if

     if (iproc==0) then
        call yaml_close_sequence() !itrep
     end if


     if (opt%itrpmax > 1) then

        !recalculate orbitals occupation numbers 
        call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)

        call eigensystem_info(iproc,nproc,opt%gnrm,&
             KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f,&
             KSwfn%orbs,KSwfn%psi)

        !stop the partial timing counter if necessary
        if (endlooprp) then
           call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')
           exit rhopot_loop
        end if

        opt%gnrm =1.d10
        KSwfn%diis%energy_min=1.d10
        ! this line can be commented
        !KSwfn%diis%alpha=in%alphadiis
     end if

     if (iproc == 0) then
        !yaml output
        !summarize the key elements in the opt%itrp element
        if (opt%itrp >1) then
           call yaml_map('RhoPot Delta','*rpnrm'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))))
           call yaml_map('Energies','*final'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))))
!!$           call yaml_comment('End RhoPot Iterations, itrp:'//&
!!$                yaml_toa(opt%itrp,fmt='(i6)'))
        end if
     end if
     if (opt%c_obj /= 0) then
        call optloop_emit_iter(opt, OPTLOOP_HAMILTONIAN, energs, iproc, nproc)
     end if

     opt%itrp = opt%itrp + 1
  end do rhopot_loop

!!$  if (iproc ==0) then
!!$     call PAPIF_flops(rtime, ptime, flpops, mflops,ierr)
!!$
!!$     write (*,90) rtime, ptime, flpops, mflops
!!$
!!$90   format('           Real time (secs) :', f15.3, &
!!$          /'            CPU time (secs) :', f15.3,&
!!$          /'Floating point instructions :', i15,&
!!$          /'                     MFLOPS :', f15.3)
!!$
!!$
!!$  end if


  if (opt%c_obj /= 0) then
     call optloop_emit_done(opt, OPTLOOP_HAMILTONIAN, energs, iproc, nproc)
  end if
  if (iproc==0) call yaml_close_sequence() !opt%itrp
  !recuperate the information coming from the last iteration (useful for post-processing of the document)
  !only if everything got OK
  if (iproc==0 .and. opt%infocode == BIGDFT_SUCCESS) &
       call yaml_map('Last Iteration','*FINAL'//trim(adjustl(yaml_toa(opt%itrep,fmt='(i3.3)'))))

  !!do i_all=1,size(rhopot)
  !!    write(10000+iproc,*) rhopot(i_all)
  !!end do
  !!do i_all=1,size(psi)
  !!    write(11000+iproc,*) psi(i_all)
  !!end do
  !!do i_all=1,size(psi)
  !!    write(12000+iproc,*) psi(i_all)
  !!end do

  !this warning can be deplaced in write_energies
  if (inputpsi /= INPUT_PSI_EMPTY) then
     energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
     !write this warning only if the system is closed shell
     call check_closed_shell(KSwfn%orbs,lcs)  
     if (abs(energs%evsum-energs%ebs) > 1.d-8 .and. iproc==0 .and. lcs) then
        call yaml_newline()
        call yaml_open_map('Energy inconsistencies')
        call yaml_map('Band Structure Energy',energs%ebs,fmt='(1pe22.14)')
        call yaml_map('Sum of Eigenvalues',energs%evsum,fmt='(1pe22.14)')
        if (energs%evsum /= 0.0_gp) call yaml_map('Relative inconsistency',(energs%ebs-energs%evsum)/energs%evsum,fmt='(1pe9.2)')
        call yaml_close_map()
        !write( *,'(1x,a,2(1x,1pe20.13))')&
        !  &   'Difference:evsum,energybs',energs%evsum,energs%ebs
     end if
  end if
  ! Clean KSwfn parts only needed in the SCF loop.
  call kswfn_free_scf_data(KSwfn, (nproc > 1))
  ! Clean denspot parts only needed in the SCF loop.
  call denspot_free_history(denspot)

END SUBROUTINE kswfn_optimization_loop


subroutine kswfn_post_treatments(iproc, nproc, KSwfn, tmb, linear, &
     & fxyz, fnoise, fion, fdisp, fpulay, &
     & strten, pressure, ewaldstr, xcstr, &
     & GPU, denspot, atoms, rxyz, nlpsp, &
     & output_denspot, dir_output, gridformat, refill_proj, calculate_dipole)
  use module_base
  use module_types
  use module_interfaces, except_this_one => kswfn_post_treatments
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  use communications_base, only: deallocate_comms_linear
  use communications, only: synchronize_onesided_communication

  implicit none

  !Arguments
  type(DFT_wavefunction), intent(in) :: KSwfn
  type(DFT_wavefunction), intent(inout) :: tmb
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_local_fields), intent(inout) :: denspot
  type(atoms_data), intent(in) :: atoms
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  logical, intent(in) :: linear, refill_proj, calculate_dipole
  integer, intent(in) :: output_denspot, iproc, nproc
  character(len = *), intent(in) :: dir_output
  character(len = *), intent(in) :: gridformat
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: fdisp, fion, fpulay
  real(dp), dimension(6), intent(in) :: ewaldstr, xcstr
  real(gp), intent(out) :: fnoise, pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), dimension(3, atoms%astruct%nat), intent(out) :: fxyz

  !Local variables
  character(len = *), parameter :: subname = "kswfn_post_treatments"
  integer :: i_all, i_stat, jproc, nsize_psi, imode
  real(dp), dimension(6) :: hstrten
  real(gp) :: ehart_fake


  !manipulate scatter array for avoiding the GGA shift
!!$     call dpbox_repartition(denspot%dpbox%iproc,denspot%dpbox%nproc,atoms%astruct%geocode,'D',1,denspot%dpbox)
  !n3d=n3p
  denspot%dpbox%n3d=denspot%dpbox%n3p
  !i3xcsh=0
  denspot%dpbox%i3s=denspot%dpbox%i3s+denspot%dpbox%i3xcsh
  denspot%dpbox%i3xcsh=0
  do jproc=0,denspot%dpbox%mpi_env%nproc-1
     !n3d=n3p
     denspot%dpbox%nscatterarr(jproc,1)=denspot%dpbox%nscatterarr(jproc,2)
     !i3xcsh=0
     denspot%dpbox%nscatterarr(jproc,4)=0
     !the same for the density
     denspot%dpbox%ngatherarr(:,3)=denspot%dpbox%ngatherarr(:,1)
  end do
  !change communication scheme to LDA case
  !only in the case of no PSolver tasks
  if (denspot%dpbox%mpi_env%nproc < nproc) then
     denspot%rhod%icomm=0
     denspot%rhod%nrhotot=denspot%dpbox%ndims(3)
  else
     denspot%rhod%icomm=1
     denspot%rhod%nrhotot=sum(denspot%dpbox%nscatterarr(:,1))
  end if

  if (linear) then
     if (denspot%dpbox%ndimpot>0) then
        allocate(denspot%pot_work(denspot%dpbox%ndimpot+ndebug),stat=i_stat)
        call memocc(i_stat,denspot%pot_work,'denspot%pot_work',subname)
     else
        allocate(denspot%pot_work(1+ndebug),stat=i_stat)
        call memocc(i_stat,denspot%pot_work,'denspot%pot_work',subname)
     end if
     ! Density already present in denspot%rho_work
     call vcopy(denspot%dpbox%ndimpot,denspot%rho_work(1),1,denspot%pot_work(1),1)
     call H_potential('D',denspot%pkernel,denspot%pot_work,denspot%pot_work,ehart_fake,&
          0.0_dp,.false.,stress_tensor=hstrten)
  else
     call density_and_hpot(denspot%dpbox,atoms%astruct%sym,KSwfn%orbs,KSwfn%Lzd,&
          denspot%pkernel,denspot%rhod, GPU, denspot%xc, &
          & KSwfn%psi,denspot%rho_work,denspot%pot_work,hstrten)
  end if

  !xc stress, diagonal for the moment
  if (atoms%astruct%geocode=='P') then
     if (atoms%astruct%sym%symObj >= 0) call symm_stress((iproc==0),xcstr,atoms%astruct%sym%symObj)
  end if

  if (calculate_dipole) then
     ! calculate dipole moment associated to the charge density
     call calc_dipole(denspot%dpbox,denspot%dpbox%nrhodim,atoms,rxyz,denspot%rho_work,.false.)
  end if
  !plot the density on the cube file
  !to be done either for post-processing or if a restart is to be done with mixing enabled
  if (((output_denspot >= output_denspot_DENSITY))) then
     if (iproc == 0) call yaml_map('Writing electronic density in file','electronic_density'//gridformat)

     call plot_density(iproc,nproc,trim(dir_output)//'electronic_density' // gridformat,&
          atoms,rxyz,denspot%dpbox,denspot%dpbox%nrhodim,denspot%rho_work)

     if (associated(denspot%rho_C)) then
        if (iproc == 0) call yaml_map('Writing core density in file','grid core_density'//gridformat)
        call plot_density(iproc,nproc,trim(dir_output)//'core_density' // gridformat,&
             atoms,rxyz,denspot%dpbox,1,denspot%rho_C(1,1,denspot%dpbox%i3xcsh:,1))
     end if
  end if
  !plot also the electrostatic potential
  if (output_denspot == output_denspot_DENSPOT) then
     if (iproc == 0) call yaml_map('Writing Hartree potential in file','hartree_potential'//gridformat)
     call plot_density(iproc,nproc,trim(dir_output)//'hartree_potential' // gridformat, &
          atoms,rxyz,denspot%dpbox,denspot%dpbox%nrhodim,denspot%pot_work)
  end if

  !     !plot also the electrostatic potential
  !     if (output_denspot == output_denspot_DENSPOT .and. DoLastRunThings) then
  !        if (iproc == 0) write(*,*) 'writing hartree_potential' // gridformat
  !        call plot_density(iproc,nproc,trim(dir_output)//'hartree_potential' // gridformat, &
  !             atoms,rxyz,denspot%dpbox,1,pot)
  !     end if
  !
  call timing(iproc,'Forces        ','ON')

  ! Calculate the forces. Pass the Pulay forces in the linear scaling case.
  if (linear) then
     imode = 1
     nsize_psi=1
     ! This is just to save memory, since calculate_forces will require quite a lot
     call deallocate_comms_linear(tmb%collcom)
     call deallocate_comms_linear(tmb%ham_descr%collcom)
     call deallocate_comms_linear(tmb%collcom_sr)
  else
     imode = 0
     nsize_psi = (KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f)*KSwfn%orbs%nspinor*KSwfn%orbs%norbp
  end if
  call calculate_forces(iproc,nproc,denspot%pkernel%mpi_env%nproc,KSwfn%Lzd%Glr,atoms,KSwfn%orbs,nlpsp,rxyz,&
       KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
       denspot%dpbox%i3s+denspot%dpbox%i3xcsh,denspot%dpbox%n3p,&
       denspot%dpbox%nrhodim,refill_proj,denspot%dpbox%ngatherarr,denspot%rho_work,&
       denspot%pot_work,denspot%V_XC,nsize_psi,KSwfn%psi,fion,fdisp,fxyz,&
       ewaldstr,hstrten,xcstr,strten,fnoise,pressure,denspot%psoffset,imode,tmb,fpulay)

  i_all=-product(shape(denspot%rho_work))*kind(denspot%rho_work)
  deallocate(denspot%rho_work,stat=i_stat)
  call memocc(i_stat,i_all,'denspot%rho',subname)
  i_all=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work,stat=i_stat)
  call memocc(i_stat,i_all,'denspot%pot_work',subname)
  nullify(denspot%rho_work,denspot%pot_work)

  if (linear) then
     ! to eventually be better sorted
     call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
     call deallocate_p2pComms(tmb%ham_descr%comgp, subname)
     call deallocate_local_zone_descriptors(tmb%ham_descr%lzd, subname)
     call deallocate_comms_linear(tmb%ham_descr%collcom)
     call deallocate_auxiliary_basis_function(subname, tmb%ham_descr%psi, tmb%hpsi)

!!!! TEST ##################
     !!fxyz=0.d0
     !!tmb%psi(1:KSwfn%orbs%npsidim_orbs)=KSwfn%psi(1:KSwfn%orbs%npsidim_orbs)
     !!tmb%wfnmd%density_kernel=0.d0
     !!do i_stat=1,KSwfn%orbs%norb
     !!    tmb%wfnmd%density_kernel(i_stat,i_stat)=1.d0
     !!end do
     !!call  nonlocal_forces(iproc,tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
     !! atoms,rxyz,&
     !! KSwfn%orbs,nlpsp,proj,tmb%lzd%glr%wfd,KSwfn%psi,fxyz,refill_proj,strten)
     !!call nonlocal_forces_linear(iproc,nproc,tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),&
     !!     KSwfn%Lzd%hgrids(3),atoms,rxyz,&
     !!     tmb%orbs,nlpsp,proj,tmb%lzd,tmb%psi,tmb%wfnmd%density_kernel,fxyz,refill_proj,strten)
     !!call nonlocal_forces_linear(iproc,nproc,tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),&
     !!     KSwfn%Lzd%hgrids(3),atoms,rxyz,&
     !!     tmb%orbs,nlpsp,proj,tmb%ham_descr%lzd,tmb%ham_descr%psi,tmb%wfnmd%density_kernel,fxyz,refill_proj,strten)
     !!if (nproc > 1) then
     !!   call mpiallred(fxyz(1,1),3*atoms%astruct%nat,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     !!end if
     !!if (iproc==0) then
     !!     do iat=1,atoms%astruct%nat
     !!         write(*,'(a,3es18.8)') 'new forces',fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
     !!     end do 
     !!end if 
!!!! #######################
  end if
  
  !!stop
  call timing(iproc,'Forces        ','OF')
END SUBROUTINE kswfn_post_treatments

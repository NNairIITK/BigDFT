!> @file
!!  Routines related to system properties
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Initialize the objects needed for the computation: basis sets, allocate required space
subroutine system_initialization(iproc,nproc,dump,inputpsi,input_wf_format,dry_run,&
     & in,atoms,rxyz,OCLconv,&
     orbs,lnpsidim_orbs,lnpsidim_comp,lorbs,Lzd,Lzd_lin,nlpsp,comms,shift,radii_cf,&
     ref_frags, denspot, locregcenters, inwhichlocreg_old, onwhichatom_old,output_grid)
  use module_base
  use module_types
  use module_interfaces, fake_name => system_initialization
  use module_xc
  use module_fragments
  use gaussians, only: gaussian_basis, nullify_gaussian_basis
  use vdwcorrection
  use yaml_output
  use module_atoms, only: set_symmetry_data
  use communications_base, only: comms_cubic
  use communications_init, only: orbitals_communicators
  implicit none
  integer, intent(in) :: iproc,nproc 
  logical, intent(in) :: dry_run, dump
  integer, intent(out) :: input_wf_format, lnpsidim_orbs, lnpsidim_comp
  integer, intent(inout) :: inputpsi
  type(input_variables), intent(in) :: in 
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  logical, intent(in) :: OCLconv
  type(orbitals_data), intent(inout) :: orbs, lorbs
  type(local_zone_descriptors), intent(inout) :: Lzd, Lzd_lin
  type(DFT_PSP_projectors), intent(out) :: nlpsp
  type(comms_cubic), intent(out) :: comms
  real(gp), dimension(3), intent(out) :: shift  !< shift on the initial positions
  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  type(system_fragment), dimension(:), pointer :: ref_frags
  real(kind=8),dimension(3,atoms%astruct%nat),intent(inout),optional :: locregcenters
  integer,dimension(:),pointer,optional:: inwhichlocreg_old, onwhichatom_old
  type(DFT_local_fields), intent(out), optional :: denspot
  logical, intent(in), optional :: output_grid
  !local variables
  character(len = *), parameter :: subname = "system_initialization"
  integer :: nB,nKB,nMB,ii,iat,iorb,iatyp,nspin_ig,norbe,norbsc,ifrag,nspinor
  real(gp), dimension(3) :: h_input
  logical:: present_inwhichlocreg_old, present_onwhichatom_old, output_grid_, frag_allocated
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(kind=8), dimension(:), allocatable :: locrad
  !Note proj_G should be filled for PAW:
  type(gaussian_basis),dimension(atoms%astruct%ntypes)::proj_G
  call f_routine(id=subname)
  !nullify dummy variables only used for PAW:
  do iatyp=1,atoms%astruct%ntypes
    call nullify_gaussian_basis(proj_G(iatyp))
  end do

  output_grid_ = .false.
  if (present(output_grid)) output_grid_ = output_grid

  if (iproc == 0 .and. dump) &
       & call print_atomic_variables(atoms, radii_cf, max(in%hx,in%hy,in%hz), &
       & in%ixc, in%dispersion)

  !grid spacings of the zone descriptors (not correct, the set is done by system size)
  Lzd=default_lzd()
  h_input=(/ in%hx, in%hy, in%hz /)
  call lzd_set_hgrids(Lzd,h_input) 

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(atoms,rxyz,radii_cf,in%crmult,in%frmult,&
       Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),OCLconv,Lzd%Glr,shift)
  if (iproc == 0 .and. dump) &
       & call print_atoms_and_grid(Lzd%Glr, atoms, rxyz, shift, &
       & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3))
  if (present(locregcenters)) then
      do iat=1,atoms%astruct%nat
          locregcenters(1:3,iat)=locregcenters(1:3,iat)-shift(1:3)
          if (locregcenters(1,iat)<dble(0)*lzd%hgrids(1) .or. locregcenters(1,iat)>dble(lzd%glr%d%n1)*lzd%hgrids(1) .or. &
              locregcenters(2,iat)<dble(0)*lzd%hgrids(2) .or. locregcenters(2,iat)>dble(lzd%glr%d%n2)*lzd%hgrids(2) .or. &
              locregcenters(3,iat)<dble(0)*lzd%hgrids(3) .or. locregcenters(3,iat)>dble(lzd%glr%d%n3)*lzd%hgrids(3)) then
              stop 'locregcenter outside of global box!'
          end if
      end do
  end if


  if (present(denspot)) then
     call initialize_DFT_local_fields(denspot, in%ixc, in%nspin)

     !here the initialization of dpbox can be set up
     call dpbox_set(denspot%dpbox,Lzd,denspot%xc,iproc,nproc,bigdft_mpi%mpi_comm, &
          & in%PSolver_groupsize, in%SIC%approach, atoms%astruct%geocode, in%nspin)

     ! Create the Poisson solver kernels.
     call system_initKernels(.true.,iproc,nproc,atoms%astruct%geocode,in,denspot)
     call system_createKernels(denspot, (verbose > 1))
  end if

  ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
  call createWavefunctionsDescriptors(iproc,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),atoms,&
       rxyz,radii_cf,in%crmult,in%frmult,Lzd%Glr, output_grid_)
  if (iproc == 0 .and. dump) call print_wfd(Lzd%Glr%wfd)

  ! Create global orbs data structure.
  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  call orbitals_descriptors(iproc, nproc,in%gen_norb,in%gen_norbu,in%gen_norbd,in%nspin,nspinor,&
       in%gen_nkpt,in%gen_kpt,in%gen_wkpt,orbs,.false.)
  orbs%occup(1:orbs%norb*orbs%nkpts) = in%gen_occup
  if (dump .and. iproc==0) call print_orbitals(orbs, atoms%astruct%geocode)
  ! Create linear orbs data structure.
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
      .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     if (present(locregcenters)) then
         call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, locregcenters, lorbs)
     else
         call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, atoms%astruct%rxyz, lorbs)
     end if

     ! There are needed for the restart (at least if the atoms have moved...)
     present_inwhichlocreg_old = present(inwhichlocreg_old)
     present_onwhichatom_old = present(onwhichatom_old)
     if (present_inwhichlocreg_old .and. .not.present_onwhichatom_old &
         .or. present_onwhichatom_old .and. .not.present_inwhichlocreg_old) then
         call yaml_warning('inwhichlocreg_old and onwhichatom_old should be present at the same time')
         stop 
     end if
     if (present_inwhichlocreg_old .and. present_onwhichatom_old) then
         call vcopy(lorbs%norb, onwhichatom_old(1), 1, lorbs%onwhichatom(1), 1)
         !call vcopy(lorbs%norb, inwhichlocreg_old(1), 1, lorbs%inwhichlocreg(1), 1)
         !use onwhichatom to build the new inwhichlocreg (because the old inwhichlocreg can be ordered differently)
         ii = 0
         do iat=1, atoms%astruct%nat
            do iorb=1,lorbs%norb
               if(iat ==  lorbs%onwhichatom(iorb)) then
                  ii = ii + 1
                  lorbs%inwhichlocreg(iorb)= ii
               end if
            end do 
         end do
         !i_all=-product(shape(inwhichlocreg_old))*kind(inwhichlocreg_old)
         !deallocate(inwhichlocreg_old,stat=i_stat)
         !call memocc(i_stat,i_all,'inwhichlocreg_old',subname)
         !i_all=-product(shape(onwhichatom_old))*kind(onwhichatom_old)
         !deallocate(onwhichatom_old,stat=i_stat)
         !call memocc(i_stat,i_all,'onwhichatom_old',subname)
     end if
  end if
  !In the case in which the number of orbitals is not "trivial" check whether they are too many
  if (inputpsi /= INPUT_PSI_RANDOM) then

     ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files)
     norbsc_arr = f_malloc((/ atoms%natsc+1, in%nspin /),id='norbsc_arr')
     locrad = f_malloc(atoms%astruct%nat,id='locrad')

     !calculate the inputguess orbitals
     !spin for inputguess orbitals
     if (in%nspin==4) then
        nspin_ig=1
     else
        nspin_ig=in%nspin
     end if

     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(atoms,norbe,norbsc,nspin_ig,orbs%nspinor,&
          norbsc_arr,locrad)

     if (in%nspin==4) then
        !in that case the number of orbitals doubles
        norbe=2*norbe
     end if

     ! De-allocations
     call f_free(locrad)
     call f_free(norbsc_arr)

     ! Check the maximum number of orbitals
     if (in%nspin==1 .or. in%nspin==4) then
        if (orbs%norb>norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals (',orbs%norb,&
                &   ') must not be greater than the number of orbitals (',norbe,&
                &   ') generated from the input guess.'
           stop
        end if
     else if (in%nspin == 2) then
        if (orbs%norbu > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals up (',orbs%norbu,&
                &   ') must not be greater than the number of orbitals (',norbe,&
                &   ') generated from the input guess.'
           stop
        end if
        if (orbs%norbd > norbe) then
           write(*,'(1x,a,i0,a,i0,a)') 'The number of orbitals down (',orbs%norbd,&
                &   ') must not be greater than the number of orbitals (',norbe,&
                &   ') generated from the input guess.'
           stop
        end if
     end if
  end if

  !allocate communications arrays (allocate it before Projectors because of the definition
  !of iskpts and nkptsp)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)  
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
      .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     if(iproc==0 .and. dump) call print_orbital_distribution(iproc, nproc, lorbs)
  end if

  if (iproc == 0 .and. dump) then
     nB=max(orbs%npsidim_orbs,orbs%npsidim_comp)*8
     nMB=nB/1024/1024
     nKB=(nB-nMB*1024*1024)/1024
     nB=modulo(nB,1024)
     call yaml_map('Wavefunctions memory occupation for root MPI process',&
          trim(yaml_toa(nMB,fmt='(i5)'))//' MB'//trim(yaml_toa(nKB,fmt='(i5)'))//&
          ' KB'//trim(yaml_toa(nB,fmt='(i5)'))//' B')
!!$     write(*,'(1x,a,3(i5,a))') &
!!$       'Wavefunctions memory occupation for root MPI process: ',&
!!$       nMB,' MB ',nKB,' KB ',nB,' B'
  end if
  ! Done orbs

  ! fragment initializations - if not a fragment calculation, set to appropriate dummy values
  if (inputpsi == INPUT_PSI_DISK_LINEAR .or. in%lin%fragment_calculation) then
     allocate(ref_frags(in%frag%nfrag_ref))
     do ifrag=1,in%frag%nfrag_ref
        ref_frags(ifrag)=fragment_null()
     end do
     call init_fragments(in,lorbs,atoms%astruct,ref_frags)
    frag_allocated=.true.
  else
     nullify(ref_frags)
  end if

  call input_check_psi_id(inputpsi, input_wf_format, in%dir_output, &
       orbs, lorbs, iproc, nproc, in%frag%nfrag_ref, in%frag%dirname, ref_frags)

  ! we need to deallocate the fragment arrays we just allocated as not a restart calculation so this is no longer needed
  if (frag_allocated .and. (.not. in%lin%fragment_calculation) .and. inputpsi /= INPUT_PSI_DISK_LINEAR) then
      do ifrag=1,in%frag%nfrag_ref
         ref_frags(ifrag)%astruct_frg%nat=-1
         ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
         call fragment_free(ref_frags(ifrag))
         !ref_frags(ifrag)=fragment_null()
      end do
     deallocate(ref_frags)
  end if

  ! See if linear scaling should be activated and build the correct Lzd 
  call check_linear_and_create_Lzd(iproc,nproc,in%linear,Lzd,atoms,orbs,in%nspin,rxyz)

  lzd_lin=default_lzd()
  call nullify_local_zone_descriptors(lzd_lin)
  lzd_lin%nlr = 0

  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
     .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     call copy_locreg_descriptors(Lzd%Glr, lzd_lin%glr)
     call lzd_set_hgrids(lzd_lin, Lzd%hgrids)
     if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
         !!write(*,*) 'rxyz',rxyz
         !!write(*,*) 'locregcenters',locregcenters
         if (present(locregcenters)) then
            call lzd_init_llr(iproc, nproc, in, atoms%astruct, locregcenters, lorbs, lzd_lin)
        else
            call lzd_init_llr(iproc, nproc, in, atoms%astruct, atoms%astruct%rxyz, lorbs, lzd_lin)
        end if

     else
        call initialize_linear_from_file(iproc,nproc,in%frag,atoms%astruct,rxyz,lorbs,lzd_lin,&
             input_wf_format,in%dir_output,'minBasis',ref_frags)
        !what to do with derivatives?
     end if

     call initLocregs(iproc, nproc, lzd_lin, Lzd_lin%hgrids(1), Lzd_lin%hgrids(2),Lzd_lin%hgrids(3), &
          atoms%astruct, lorbs, Lzd_lin%Glr, 's')
     call update_wavefunctions_size(lzd_lin,lnpsidim_orbs,lnpsidim_comp,lorbs,iproc,nproc)

  end if

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call createProjectorsArrays(Lzd%Glr,rxyz,atoms,orbs,&
       radii_cf,in%frmult,in%frmult,Lzd%hgrids(1),Lzd%hgrids(2),&
       Lzd%hgrids(3),dry_run,nlpsp,proj_G)
  if (iproc == 0 .and. dump) call print_nlpsp(nlpsp)
  !the complicated part of the descriptors has not been filled
  if (dry_run) then
     call f_release_routine()
     return
  end if
  !calculate the partitioning of the orbitals between the different processors
!  print *,'here the localization regions should have been filled already'
!  stop

  if (present(denspot)) then
     !here dpbox can be put as input
     call density_descriptors(iproc,nproc,denspot%xc,in%nspin,in%crmult,in%frmult,atoms,&
          denspot%dpbox,in%rho_commun,rxyz,radii_cf,denspot%rhod)
     !allocate the arrays.
     call allocateRhoPot(Lzd%Glr,in%nspin,atoms,rxyz,denspot)
  end if

  !calculate the irreductible zone for this region, if necessary.
  call set_symmetry_data(atoms%astruct%sym,atoms%astruct%geocode, &
       & Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i, in%nspin)

  ! A message about dispersion forces.
  call vdwcorrection_initializeparams(in%ixc, in%dispersion)

  !check the communication distribution
  if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
     .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
      call check_communications(iproc,nproc,orbs,Lzd,comms)
  else
      ! Do not call check_communication, since the value of orbs%npsidim_orbs is wrong
      if(iproc==0) call yaml_warning('Do not call check_communications in the linear scaling version!')
      !if(iproc==0) write(*,*) 'WARNING: do not call check_communications in the linear scaling version!'
  end if

  !Check if orbitals and electrons
  if (orbs%norb*orbs%nkpts == 0) &
     & call f_err_throw('No electrons in the system! Check your input variables or atomic positions.', &
     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)

  call f_release_routine()
  !---end of system definition routine

END SUBROUTINE system_initialization


subroutine system_initKernels(verb, iproc, nproc, geocode, in, denspot)
  use module_types
  use module_xc
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use module_base
  implicit none
  logical, intent(in) :: verb
  integer, intent(in) :: iproc, nproc
  character, intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  type(input_variables), intent(in) :: in
  type(DFT_local_fields), intent(inout) :: denspot

  integer, parameter :: ndegree_ip = 16

  denspot%pkernel=pkernel_init(verb, iproc,nproc,in%matacc%PSolver_igpu,&
       geocode,denspot%dpbox%ndims,denspot%dpbox%hgrids,ndegree_ip,mpi_env=denspot%dpbox%mpi_env)
  !create the sequential kernel if the exctX parallelisation scheme requires it
  if ((xc_exctXfac(denspot%xc) /= 0.0_gp .and. in%exctxpar=='OP2P' .or. in%SIC%alpha /= 0.0_gp)&
       .and. denspot%dpbox%mpi_env%nproc > 1) then
     !the communicator of this kernel is bigdft_mpi%mpi_comm
     !this might pose problems when using SIC or exact exchange with taskgroups
     denspot%pkernelseq=pkernel_init(iproc==0 .and. verb,0,1,in%matacc%PSolver_igpu,&
          geocode,denspot%dpbox%ndims,denspot%dpbox%hgrids,ndegree_ip)
  else 
     denspot%pkernelseq = denspot%pkernel
  end if

END SUBROUTINE system_initKernels

subroutine system_createKernels(denspot, verb)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  implicit none
  logical, intent(in) :: verb
  type(DFT_local_fields), intent(inout) :: denspot

  call pkernel_set(denspot%pkernel,verb)
    !create the sequential kernel if pkernelseq is not pkernel
  if (denspot%pkernelseq%mpi_env%nproc == 1 .and. denspot%pkernel%mpi_env%nproc /= 1) then
     call pkernel_set(denspot%pkernelseq,.false.)
  else
     denspot%pkernelseq = denspot%pkernel
  end if

END SUBROUTINE system_createKernels


!> Calculate the important objects related to the physical properties of the system
subroutine system_properties(iproc,nproc,in,atoms,orbs,radii_cf)
  use module_base
  use module_types
  use module_interfaces, except_this_one => system_properties
  implicit none
  integer, intent(in) :: iproc,nproc
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(inout) :: atoms
  type(orbitals_data), intent(inout) :: orbs
  real(gp), dimension(atoms%astruct%ntypes,3), intent(out) :: radii_cf
  !local variables
  !n(c) character(len=*), parameter :: subname='system_properties'
  integer :: nspinor

  call read_radii_variables(atoms, radii_cf, in%crmult, in%frmult, in%projrad)
!!$  call read_atomic_variables(atoms, trim(in%file_igpop),in%nspin)
  if (iproc == 0) call print_atomic_variables(atoms, radii_cf, max(in%hx,in%hy,in%hz), in%ixc, in%dispersion)
  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  call orbitals_descriptors(iproc, nproc,in%gen_norb,in%gen_norbu,in%gen_norbd,in%nspin,nspinor,&
       in%gen_nkpt,in%gen_kpt,in%gen_wkpt,orbs,.false.)
  orbs%occup(1:orbs%norb*orbs%nkpts) = in%gen_occup
  if (iproc==0) call print_orbitals(orbs, atoms%astruct%geocode)

  !Check if orbitals and electrons
  if (orbs%norb*orbs%nkpts == 0) &
     & call f_err_throw('No electrons in the system. Check your input variables or atomic positions', &
     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)

END SUBROUTINE system_properties


!> Check for the need of a core density and fill the rhocore array which
!! should be passed at the rhocore pointer
subroutine calculate_rhocore(at,d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: i3s,n3d,i3xcsh,n3p
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  type(grid_dimensions), intent(in) :: d
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(:,:,:,:), pointer :: rhocore
  !local variables
  character(len=*), parameter :: subname='calculate_rhocore'
  integer :: ityp,iat,j3,i1,i2 !,ierr,ind
  real(wp) :: tt
  real(gp) :: rx,ry,rz,rloc,cutoff
  

  !check for the need of a nonlinear core correction
!!$  donlcc=.false.
!!$  chk_nlcc: do ityp=1,at%astruct%ntypes
!!$     filename = 'nlcc.'//at%astruct%atomnames(ityp)
!!$
!!$     inquire(file=filename,exist=exists)
!!$     if (exists) then
!!$        donlcc=.true.
!!$        exit chk_nlcc
!!$     end if
!!$  end do chk_nlcc

  if (at%donlcc) then
     !allocate pointer rhocore
     rhocore = f_malloc_ptr((/ d%n1i , d%n2i , n3d , 10+ndebug /),id='rhocore')
     !initalise it 
     if (n3d > 0) call to_zero(d%n1i*d%n2i*n3d*10,rhocore(1,1,1,1))
     !perform the loop on any of the atoms which have this feature
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
!!$        filename = 'nlcc.'//at%astruct%atomnames(ityp)
!!$        inquire(file=filename,exist=exists)
!!$        if (exists) then
        if (at%nlcc_ngv(ityp)/=UNINITIALIZED(1) .or.&
             at%nlcc_ngc(ityp)/=UNINITIALIZED(1) ) then
           if (bigdft_mpi%iproc == 0) call yaml_map('NLCC, Calculate core density for atom',trim(at%astruct%atomnames(ityp)))
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           cutoff=10.d0*rloc

           call calc_rhocore_iat(bigdft_mpi%iproc,at,ityp,rx,ry,rz,cutoff,hxh,hyh,hzh,&
                d%n1,d%n2,d%n3,d%n1i,d%n2i,d%n3i,i3s,n3d,rhocore)

        end if
     end do

     !calculate total core charge in the grid
     !In general this should be really bad

!!$     do j3=1,n3d
!!$        tt=0.0_wp
!!$        do i2=1,d%n2i
!!$           do i1=1,d%n1i
!!$              !ind=i1+(i2-1)*d%n1i+(j3+i3xcsh-1)*d%n1i*d%n2i
!!$              tt=tt+rhocore(i1,i2,j3,1)
!!$           enddo
!!$        enddo
!!$        write(17+iproc,*)j3+i3s-1,tt
!!$     enddo
!!$call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$stop
     tt=0.0_wp
     do j3=1,n3p
        do i2=1,d%n2i
           do i1=1,d%n1i
              !ind=i1+(i2-1)*d%n1i+(j3+i3xcsh-1)*d%n1i*d%n2i
              tt=tt+rhocore(i1,i2,j3+i3xcsh,1)
           enddo
        enddo
     enddo

     if (bigdft_mpi%nproc > 1) call mpiallred(tt,1,MPI_SUM,bigdft_mpi%mpi_comm)
     tt=tt*hxh*hyh*hzh
     if (bigdft_mpi%iproc == 0) call yaml_map('Total core charge on the grid (To be compared with analytic one)', tt,fmt='(f15.7)')

  else
     !No NLCC needed, nullify the pointer 
     nullify(rhocore)
  end if

END SUBROUTINE calculate_rhocore


subroutine psp_from_file(filename, nzatom, nelpsp, npspcode, &
     & ixcpsp, psppar, donlcc, rcore, qcore, radii_cf, exists, pawpatch)
  use module_base
  use ao_inguess
  implicit none
  
  character(len = *), intent(in) :: filename
  integer, intent(out) :: nzatom, nelpsp, npspcode, ixcpsp
  real(gp), intent(out) :: psppar(0:4,0:6), radii_cf(3), rcore, qcore
  logical, intent(out) :: exists, pawpatch
  logical, intent(inout) ::  donlcc
  !ALEX: Some local variables
  real(gp):: fourpi, sqrt2pi
  character(len=2) :: symbol
  character(len=20) :: skip

  integer :: ierror, ierror1, i, j, nn, nlterms, nprl, l, nzatom_, nelpsp_, npspcode_
  real(dp) :: nelpsp_dp,nzatom_dp
  character(len=100) :: line

  radii_cf = UNINITIALIZED(1._gp)
  pawpatch = .false.
  inquire(file=trim(filename),exist=exists)
  if (.not. exists) return

  ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
  open(unit=11,file=trim(filename),status='old',iostat=ierror)
  !Check the open statement
  if (ierror /= 0) then
     write(*,*) ': Failed to open the file (it must be in ABINIT format!): "',&
          trim(filename),'"'
     stop
  end if
  read(11,*)
  read(11,*) nzatom_dp, nelpsp_dp
  nzatom=int(nzatom_dp); nelpsp=int(nelpsp_dp)
  read(11,*) npspcode, ixcpsp

  psppar(:,:)=0._gp
  if (npspcode == 2) then !GTH case
     read(11,*) (psppar(0,j),j=0,4)
     do i=1,2
        read(11,*) (psppar(i,j),j=0,3-i)
     enddo
  else if (npspcode == 3) then !HGH case
     read(11,*) (psppar(0,j),j=0,4)
     read(11,*) (psppar(1,j),j=0,3)
     do i=2,4
        read(11,*) (psppar(i,j),j=0,3)
        !ALEX: Maybe this can prevent reading errors on CRAY machines?
        read(11,*) skip !k coefficients, not used for the moment (no spin-orbit coupling)
     enddo
  else if (npspcode == 7) then !PAW Pseudos
     ! Need NC psp for input guess.
     call atomic_info(nzatom, nelpsp, symbol = symbol)
     call psp_from_data(symbol, nzatom_, nelpsp_, npspcode_, ixcpsp, &
          & psppar, exists)
     if (.not.exists) stop "Implement here."

     ! PAW format using libPAW.
     write(*,*) 'Reading of PAW atomic-data, under development'
     call psp_from_file_paw()
     
  else if (npspcode == 10) then !HGH-K case
     read(11,*) psppar(0,0),nn,(psppar(0,j),j=1,nn) !local PSP parameters
     read(11,*) nlterms !number of channels of the pseudo
     prjloop: do l=1,nlterms
        read(11,*) psppar(l,0),nprl,psppar(l,1),&
             (psppar(l,j+2),j=2,nprl) !h_ij terms
        do i=2,nprl
           read(11,*) psppar(l,i),(psppar(l,i+j+1),j=i+1,nprl) !h_ij 
        end do
        if (l==1) cycle
        do i=1,nprl
           !ALEX: Maybe this can prevent reading errors on CRAY machines?
           read(11,*)skip !k coefficients, not used
        end do
     end do prjloop
  !ALEX: Add support for reading NLCC from psppar
  else if (npspcode == 12) then !HGH-NLCC: Same as HGH-K + one additional line
     read(11,*) psppar(0,0),nn,(psppar(0,j),j=1,nn) !local PSP parameters
     read(11,*) nlterms !number of channels of the pseudo
     do l=1,nlterms
        read(11,*) psppar(l,0),nprl,psppar(l,1),&
             (psppar(l,j+2),j=2,nprl) !h_ij terms
        do i=2,nprl
           read(11,*) psppar(l,i),(psppar(l,i+j+1),j=i+1,nprl) !h_ij
        end do
        if (l==1) cycle
        do i=1,nprl
           !ALEX: Maybe this can prevent reading errors on CRAY machines?
           read(11,*) skip !k coefficients, not used
        end do
     end do 
     read(11,*) rcore, qcore
     !convert the core charge fraction qcore to the amplitude of the Gaussian
     !multiplied by 4pi. This is the convention used in nlccpar(1,:).
     fourpi=4.0_gp*pi_param!8.0_gp*dacos(0.0_gp)
     sqrt2pi=sqrt(0.5_gp*fourpi)
     qcore=fourpi*qcore*real(nzatom-nelpsp,gp)/&
          (sqrt2pi*rcore)**3
     donlcc=.true.
  else
     !if (iproc == 0) then
     write(*,'(1x,a,a)') trim(filename),&
          'unrecognized pspcode: only GTH, HGH & HGH-K pseudos (ABINIT format)'
     !end if
     stop
  end if

  if (npspcode == 7) return !Skip the rest for PAW

  !old way of calculating the radii, requires modification of the PSP files
  read(11,'(a100)',iostat=ierror) line
  if (ierror /=0) then
     !if (iproc ==0) write(*,*)&
     !     ' WARNING: last line of pseudopotential missing, put an empty line'
     line=''
  end if
  read(line,*,iostat=ierror1) radii_cf(1),radii_cf(2),radii_cf(3)
  if (ierror1 /= 0 ) then
     read(line,*,iostat=ierror) radii_cf(1),radii_cf(2)
     radii_cf(3)=radii_cf(2)
     ! Open64 behaviour, if line is PAWPATCH, then radii_cf(1) = 0.
     if (ierror /= 0) radii_cf = UNINITIALIZED(1._gp)
  end if
  pawpatch = (trim(line) == "PAWPATCH")
  do
     read(11,'(a100)',iostat=ierror) line
     if (ierror /= 0 .or. pawpatch) exit
     pawpatch = (trim(line) == "PAWPATCH")
  end do
  close(11)

contains

   subroutine psp_from_file_paw()
     use module_base
     use m_pawpsp, only: pawpsp_main, pawpsp_read_header, pawpsp_read_header_2
     use defs_basis, only: tol14, fnlen
     use m_pawrad, only: pawrad_type, pawrad_nullify, pawrad_destroy
     use m_pawtab, only: pawtab_type, pawtab_nullify, pawtab_destroy
     implicit none
     integer:: icoulomb,ipsp,ixc,lnmax
     integer:: lloc,l_size,lmax,mmax,pspcod,pspxc
     integer:: pspversion,basis_size,lmn_size
     integer:: mpsang,mqgrid_ff,mqgrid_vl,mqgrid_shp
     integer:: pawxcdev,usewvl,usexcnhat,xclevel
     integer::pspso
     real(dp):: r2well
     real(dp):: xc_denpos,zionpsp,znuclpsp
     real(dp)::epsatm,xcccrc
     character(len=fnlen):: filpsp   ! name of the psp file
     character(len = *), parameter :: subname = "psp_from_file_paw"
     type(pawrad_type):: pawrad
     type(pawtab_type):: pawtab
     integer:: comm_mpi
   !  type(paw_setup_t),optional,intent(in) :: psxml
   !!arrays
    integer:: wvl_ngauss(2)
    real(dp),allocatable:: qgrid_ff(:),qgrid_vl(:)
    real(dp),allocatable:: ffspl(:,:,:)
    real(dp),allocatable:: vlspl(:,:)
    integer:: mesh_size


     !These should be passed as arguments:
     !Defines the number of Gaussian functions for projectors
     !See ABINIT input files documentation
     wvl_ngauss=[10,10]
     icoulomb= 1 !Fake argument, this only indicates that we are inside bigdft..
                 !do not change, even if icoulomb/=1
     ipsp=1      !This is relevant only for XML.
                 !This is not yet working
     xclevel=1 ! xclevel=XC functional level (1=LDA, 2=GGA)
               ! For the moment, it will just work for LDA
     pspso=0 !No spin-orbit for the moment

   ! Read PSP header:
     rewind(11)
     call pawpsp_read_header(lloc,l_size,mmax,pspcod,pspxc,r2well,zionpsp,znuclpsp)
     call pawpsp_read_header_2(pspversion,basis_size,lmn_size)

   ! Problem lnmax are unknown here,
   ! we have to read all of the pseudo files to know it!
   ! We should change the way this is done in ABINIT:
   ! For the moment lnmax=basis_size
   ! The same problem for mpsang
     lnmax=basis_size
     lmax=l_size
   !  do ii=1,psps%npsp
   !   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   !  end do
     mpsang=lmax+1

   ! These are just useful for 
   !reciprocal space approaches (plane-waves):
     mqgrid_shp=0; mqgrid_ff=0; mqgrid_vl=0 
                             
   qgrid_ff = f_malloc(mqgrid_ff,id='qgrid_ff')
   qgrid_vl = f_malloc(mqgrid_vl,id='qgrid_vl')
   ffspl = f_malloc((/ mqgrid_ff, 2, lnmax /),id='ffspl')
   vlspl = f_malloc((/ mqgrid_vl, 2 /),id='vlspl')

   ! Define parameters:
     pawxcdev=1; usewvl=1 ; usexcnhat=0 !default
     xc_denpos=tol14
     filpsp=trim(filename)
     comm_mpi=bigdft_mpi%mpi_comm  
     mesh_size=mmax

     call pawrad_nullify(pawrad)
     call pawtab_nullify(pawtab)

     close(11)

     call pawpsp_main( &
   & pawrad,pawtab,&
   & filpsp,usewvl,icoulomb,ixc,xclevel,pawxcdev,usexcnhat,&
   & qgrid_ff,qgrid_vl,ffspl,vlspl,epsatm,xcccrc,zionpsp,znuclpsp,&
   & wvl_ngauss,comm_mpi=comm_mpi)

   !Print out data to validate this test:
     write(*,'(a)') 'PAW Gaussian projectors:'
     write(*,'("No. of Gaussians:", i4)')pawtab%wvl%pngau
     write(*,'("First five Gaussian complex coefficients:")')
     write(*,'(5("(",f13.7,",",f13.7")"))')pawtab%wvl%parg(:,1:5)
     write(*,'("First five Gaussian complex factors:")')
     write(*,'(5("(",f13.7,",",f13.7")"))')pawtab%wvl%pfac(:,1:5)
   !
     write(*,'(a)') 'GTH parameters (for initial guess):'
     write(*,'("radii_cf= ",3f10.7)')radii_cf(:)
     write(*,'("psppar(0:1,0)= ",2f10.7)')psppar(0:1,0)

   ! Destroy and deallocate objects
     call pawrad_destroy(pawrad)
     call pawtab_destroy(pawtab)

   call f_free(qgrid_ff)
   call f_free(qgrid_vl)
   call f_free(ffspl)
   call f_free(vlspl)

   !PAW is not yet working!
   !Exit here
    stop

   END SUBROUTINE psp_from_file_paw

END SUBROUTINE psp_from_file


subroutine nlcc_dim_from_file(filename, ngv, ngc, dim, read_nlcc)
  use module_base
  implicit none
  
  character(len = *), intent(in) :: filename
  integer, intent(inout) :: dim
  integer, intent(out) :: ngv, ngc
  logical, intent(out) :: read_nlcc

  integer :: ig, j
  real(gp), dimension(0:4) :: fake_nlcc

  inquire(file=filename,exist=read_nlcc)
  if (read_nlcc) then
     !associate the number of gaussians
     open(unit=79,file=filename,status='unknown')
     read(79,*)ngv
     if (ngv==0) then 
        ngv=UNINITIALIZED(1)
     else
        dim=dim+(ngv*(ngv+1)/2)
        do ig=1,(ngv*(ngv+1))/2
           read(79,*) (fake_nlcc(j),j=0,4)!jump the suitable lines (the file is organised with one element per line)
        end do
     end if
     read(79,*)ngc
     if (ngc==0) then
        ngc=UNINITIALIZED(1)
     else
        dim=dim+(ngc*(ngc+1))/2

        !better to read values in a fake array
        do ig=1,(ngc*(ngc+1))/2
           read(79,*) (fake_nlcc(j),j=0,4)!jump the suitable lines (the file is organised with one element per line)
        end do
     end if
     !no need to go further for the moment
     close(unit=79)
  else
     ngv=UNINITIALIZED(1)
     ngc=UNINITIALIZED(1)
  end if
END SUBROUTINE nlcc_dim_from_file


!> Update radii_cf and occupation for each type of atoms (related to pseudopotential)
subroutine read_radii_variables(atoms, radii_cf, crmult, frmult, projrad)
  use module_base
  use ao_inguess, only: atomic_info
  use module_types
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: crmult, frmult, projrad
  real(gp), dimension(atoms%astruct%ntypes,3), intent(out) :: radii_cf
  !Local Variables
  !integer, parameter :: nelecmax=32
  !character(len=2) :: symbol
  integer :: i,ityp!,mxpl,mxchg,nsccode
  real(gp) :: ehomo,maxrad,radfine!,rcov,rprb,amu

  do ityp=1,atoms%astruct%ntypes

     call atomic_info(atoms%nzatom(ityp),atoms%nelpsp(ityp),ehomo=ehomo)
          
     if (any(atoms%radii_cf(ityp, :) == UNINITIALIZED(1.0_gp))) then
        !assigning the radii by calculating physical parameters
        if (radii_cf(ityp,1) == UNINITIALIZED(1.0_gp)) radii_cf(ityp,1)=1._gp/sqrt(abs(2._gp*ehomo))
        radfine=100._gp
        do i=0,4
           if (atoms%psppar(i,0,ityp)/=0._gp) then
              radfine=min(radfine,atoms%psppar(i,0,ityp))
           end if
        end do
        if (radii_cf(ityp,2) == UNINITIALIZED(1.0_gp)) radii_cf(ityp,2)=radfine
        if (radii_cf(ityp,3) == UNINITIALIZED(1.0_gp)) radii_cf(ityp,3)=radfine
     else
        !Everything is already provided
        radii_cf(ityp, :) = atoms%radii_cf(ityp, :)
     end if

     ! Correct radii_cf(:,3) for the projectors.
     maxrad=0.e0_gp ! This line added by Alexey, 03.10.08, to be able to compile with -g -C
     do i=1,4
        !the maximum radii is useful only for projectors
        if (atoms%psppar(i,0,ityp)/=0._gp) then
           maxrad=max(maxrad,atoms%psppar(i,0,ityp))
        end if
     end do
     if (maxrad == 0.0_gp) then
        radii_cf(ityp,3)=0.0_gp
     else
        radii_cf(ityp,3)=max(min(crmult*radii_cf(ityp,1),projrad*maxrad)/frmult,radii_cf(ityp,2))
     end if
  enddo
END SUBROUTINE read_radii_variables


!> Calculate the number of electrons and check the polarisation (mpol)
subroutine read_n_orbitals(iproc, nelec_up, nelec_down, norbe, &
     & atoms, ncharge, nspin, mpol, norbsempty)
  use module_types, only: atoms_data
  use module_base, only: gp, f_err_throw
  use yaml_output, only: yaml_toa , yaml_warning, yaml_comment
  !use ao_inguess, only : count_atomic_shells
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: nelec_up, nelec_down, norbe
  integer, intent(in) :: ncharge, nspin, mpol, norbsempty, iproc
  !Local variables
  integer :: nelec, iat, ityp, ispinsum, ichgsum, ichg, ispol!, nspinor
  !integer, parameter :: nelecmax=32,lmax=4,noccmax=2
  !integer, dimension(lmax) :: nl
  !real(gp), dimension(noccmax,lmax) :: occup

  !calculate number of electrons and orbitals
  ! Number of electrons and number of semicore atoms
  nelec=0
  do iat=1,atoms%astruct%nat
     ityp=atoms%astruct%iatype(iat)
     nelec=nelec+atoms%nelpsp(ityp)
  enddo
  nelec=nelec-ncharge

  if(nelec < 0.0 ) then
    !if(iproc==0) write(*,*)'ERROR: Number of electrons is negative:',nelec,'.'
    !if(iproc==0) write(*,*)'FIX: decrease charge of system.'
    !call mpi_finalize(iat)
    !stop
    call f_err_throw('Number of electrons is negative:' // trim(yaml_toa(nelec)) // &
      & '. FIX: decrease charge of system.', err_name='BIGDFT_RUNTIME_ERROR')
  end if

  ! Number of orbitals
  if (nspin==1) then
     nelec_up=nelec
     nelec_down=0
  else if(nspin==4) then
     nelec_up=nelec
     nelec_down=0
  else 
     if (mod(nelec+mpol,2) /=0) then
          call f_err_throw('The mpol polarization should have the same parity of the number of electrons. ' // &
            & '(mpol=' // trim(yaml_toa(mpol)) // ' and nelec=' // trim(yaml_toa(nelec)) // ')', &
            & err_name='BIGDFT_INPUT_VARIABLES_ERROR')
        !write(*,*)'ERROR: '
        !stop
     end if
     nelec_up=min((nelec+mpol)/2,nelec)
     nelec_down=nelec-nelec_up

     !test if the spin is compatible with the input guess polarisations
     ispinsum=0
     ichgsum=0
     do iat=1,atoms%astruct%nat
        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
        ispinsum=ispinsum+ispol
        ichgsum=ichgsum+ichg
     end do

     if (ispinsum /= nelec_up-nelec_down) then
        !call yaml_warning('Total input polarisation (found ' // trim(yaml_toa(ispinsum)) &
        !     & // ') must be equal to nelec_up-nelec_down.')
        !call yaml_comment('With nelec=' // trim(yaml_toa(nelec)) &
        !     & // ' and mpol=' // trim(yaml_toa(mpol)) // &
        !     & ' nelec_up-nelec_down=' // trim((yaml_toa(nelec_up-nelec_down))))
        !stop
        call f_err_throw('Total polarisation for the input guess (found ' // trim(yaml_toa(ispinsum)) // &
           & ') must be equal to nelec_up-nelec_down ' // &
           & '(nelec=' // trim(yaml_toa(nelec)) // ', mpol=' // trim(yaml_toa(mpol)) // &
           & ', nelec_up-nelec_down=' // trim((yaml_toa(nelec_up-nelec_down))) // &
           & ', nelec_up=' // trim((yaml_toa(nelec_up))) // &
           & ', nelec_down=' // trim((yaml_toa(nelec_down))) // &
           & '). Use the keyword "IGSpin" or add a spin component for the input guess per atom.', &
           & err_name='BIGDFT_INPUT_VARIABLES_ERROR')
     end if

     if (ichgsum /= ncharge .and. ichgsum /= 0) then
        !call yaml_warning('Total input charge (found ' // trim(yaml_toa(ichgsum)) &
        !     & // ') cannot be different than charge.')
        !call yaml_comment('With charge =' // trim(yaml_toa(ncharge)) &
        !     & // ' and input charge=' // trim(yaml_toa(ichgsum)))
        !stop
        call f_err_throw('Total input charge (found ' // trim(yaml_toa(ichgsum)) // &
             & ') cannot be different than charge. With charge =' // trim(yaml_toa(ncharge)) // &
             & ' and input charge=' // trim(yaml_toa(ichgsum)), &
             & err_name='BIGDFT_INPUT_VARIABLES_ERROR')
     end if

     !now warn if there is no input guess spin polarisation
     ispinsum=0
     do iat=1,atoms%astruct%nat
        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
        ispinsum=ispinsum+abs(ispol)
     end do
     if (ispinsum == 0) then
        if (iproc==0 .and. norbsempty == 0) &
             call yaml_warning('Found no input polarisation, add it for a correct input guess')
        !write(*,'(1x,a)')&
        !     'WARNING: Found no input polarisation, add it for a correct input guess'
        !stop
     end if
  end if

  norbe = 0
  !if(nspin==4) then
  !   nspinor=4
  !else
  !   nspinor=1
  !end if
  do iat=1,atoms%astruct%nat
     !ityp=atoms%astruct%iatype(iat)
     !call count_atomic_shells(nspin,atoms%aoig(iat)%aocc,occup,nl)
     norbe=norbe+atoms%aoig(iat)%nao!nl(1)+3*nl(2)+5*nl(3)+7*nl(4)
  end do
end subroutine read_n_orbitals


!> Find the correct position of the nlcc parameters
subroutine nlcc_start_position(ityp,atoms,ngv,ngc,islcc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ityp
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: ngv,ngc,islcc
  !local variables
  integer :: ilcc,jtyp

  ilcc=0
  do jtyp=1,ityp-1
     ngv=atoms%nlcc_ngv(jtyp)
     if (ngv /= UNINITIALIZED(ngv)) ilcc=ilcc+(ngv*(ngv+1)/2)
     ngc=atoms%nlcc_ngc(jtyp)
     if (ngc /= UNINITIALIZED(ngc)) ilcc=ilcc+(ngc*(ngc+1))/2
  end do
  islcc=ilcc

  ngv=atoms%nlcc_ngv(ityp)
  if (ngv==UNINITIALIZED(1)) ngv=0
  ngc=atoms%nlcc_ngc(ityp)
  if (ngc==UNINITIALIZED(1)) ngc=0
END SUBROUTINE nlcc_start_position


!!!!> Define the descriptors of the orbitals from a given norb
!!!!! It uses the cubic strategy for partitioning the orbitals
!!!subroutine orbitals_descriptors_forLinear(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
!!!  integer, intent(in) :: nspinor
!!!  type(orbitals_data), intent(out) :: orbs
!!!  real(gp), dimension(nkpt), intent(in) :: wkpt
!!!  real(gp), dimension(3,nkpt), intent(in) :: kpt
!!!  !local variables
!!!  character(len=*), parameter :: subname='orbitals_descriptors'
!!!  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all,iiorb
!!!  integer :: mpiflag
!!!  logical, dimension(:), allocatable :: GPU_for_orbs
!!!  integer, dimension(:,:), allocatable :: norb_par !(with k-pts)
!!!
!!!
!!!  allocate(orbs%norb_par(0:nproc-1+ndebug,0:nkpt),stat=i_stat)
!!!  call memocc(i_stat,orbs%norb_par,'orbs%norb_par',subname)
!!!
!!!  !assign the value of the k-points
!!!  orbs%nkpts=nkpt
!!!  !allocate vectors related to k-points
!!!  allocate(orbs%kpts(3,orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%kpts,'orbs%kpts',subname)
!!!  allocate(orbs%kwgts(orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%kwgts,'orbs%kwgts',subname)
!!!  orbs%kpts(:,1:nkpt) = kpt(:,:)
!!!  orbs%kwgts(1:nkpt) = wkpt(:)
!!!
!!!  ! Change the wavefunctions to complex if k-points are used (except gamma).
!!!  orbs%nspinor=nspinor
!!!  if (nspinor == 1) then
!!!     if (maxval(abs(orbs%kpts)) > 0._gp) orbs%nspinor=2
!!!     !nspinor=2 !fake, used for testing with gamma
!!!  end if
!!!  orbs%nspin = nspin
!!!
!!!  !initialise the array
!!!  do jproc=0,nproc-1
!!!     orbs%norb_par(jproc,0)=0 !size 0 nproc-1
!!!  end do
!!!
!!!  !create an array which indicate which processor has a GPU associated 
!!!  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
!!!  if (.not. GPUshare) then
!!!     allocate(GPU_for_orbs(0:nproc-1+ndebug),stat=i_stat)
!!!     call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)
!!!     
!!!     if (nproc > 1) then
!!!        call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
!!!             bigdft_mpi%mpi_comm,ierr)
!!!     else
!!!        GPU_for_orbs(0)=GPUconv
!!!     end if
!!!     
!!!     i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
!!!     deallocate(GPU_for_orbs,stat=i_stat)
!!!     call memocc(i_stat,i_all,'GPU_for_orbs',subname)
!!!  end if
!!!
!!!  allocate(norb_par(0:nproc-1,orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,norb_par,'norb_par',subname)
!!!
!!!  !old system for calculating k-point repartition
!!!!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,norb,orbs%norb_par)
!!!!!$
!!!!!$  !check the distribution
!!!!!$  norb_tot=0
!!!!!$  do jproc=0,iproc-1
!!!!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!!!!$  end do
!!!!!$  !reference orbital for process
!!!!!$  orbs%isorb=norb_tot
!!!!!$  do jproc=iproc,nproc-1
!!!!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!!!!$  end do
!!!!!$
!!!!!$  if(norb_tot /= norb*orbs%nkpts) then
!!!!!$     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!!!!$     write(*,*)orbs%norb_par(:),norb*orbs%nkpts
!!!!!$     stop
!!!!!$  end if
!!!!!$
!!!!!$  !calculate the k-points related quantities
!!!!!$  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
!!!!!$  call memocc(i_stat,mykpts,'mykpts',subname)
!!!!!$
!!!!!$  call parallel_repartition_per_kpoints(iproc,nproc,orbs%nkpts,norb,orbs%norb_par,&
!!!!!$       orbs%nkptsp,mykpts,norb_par)
!!!!!$  if (orbs%norb_par(iproc) >0) then
!!!!!$     orbs%iskpts=mykpts(1)-1
!!!!!$  else
!!!!!$     orbs%iskpts=0
!!!!!$  end if
!!!!!$  i_all=-product(shape(mykpts))*kind(mykpts)
!!!!!$  deallocate(mykpts,stat=i_stat)
!!!!!$  call memocc(i_stat,i_all,'mykpts',subname)
!!!
!!!  !new system for k-point repartition
!!!  call kpts_to_procs_via_obj(nproc,orbs%nkpts,norb,norb_par)
!!!  !assign the values for norb_par and check the distribution
!!!  norb_tot=0
!!!  do jproc=0,nproc-1
!!!     if (jproc==iproc) orbs%isorb=norb_tot
!!!     do ikpt=1,orbs%nkpts
!!!        orbs%norb_par(jproc,0)=orbs%norb_par(jproc,0)+norb_par(jproc,ikpt)
!!!     end do
!!!     norb_tot=norb_tot+orbs%norb_par(jproc,0)
!!!  end do
!!!
!!!  if(norb_tot /= norb*orbs%nkpts) then
!!!     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!!     write(*,*)orbs%norb_par(:,0),norb*orbs%nkpts
!!!     stop
!!!  end if
!!!
!!!
!!!  !allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
!!!  !call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
!!!  !orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)
!!!
!!!  !this array will be reconstructed in the orbitals_communicators routine
!!!  i_all=-product(shape(norb_par))*kind(norb_par)
!!!  deallocate(norb_par,stat=i_stat)
!!!  call memocc(i_stat,i_all,'norb_par',subname)
!!!
!!!  !assign the values of the orbitals data
!!!  orbs%norb=norb
!!!  orbs%norbp=orbs%norb_par(iproc,0)
!!!  orbs%norbu=norbu
!!!  orbs%norbd=norbd
!!!
!!!  ! Modify these values
!!!  call repartitionOrbitals2(iproc, nproc, orbs%norb, orbs%norb_par, orbs%norbp, orbs%isorb)
!!!
!!!
!!!  allocate(orbs%iokpt(orbs%norbp+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
!!!
!!!  !assign the k-point to the given orbital, counting one orbital after each other
!!!  jorb=0
!!!  do ikpt=1,orbs%nkpts
!!!     do iorb=1,orbs%norb
!!!        jorb=jorb+1 !this runs over norb*nkpts values
!!!        if (jorb > orbs%isorb .and. jorb <= orbs%isorb+orbs%norbp) then
!!!           orbs%iokpt(jorb-orbs%isorb)=ikpt
!!!        end if
!!!     end do
!!!  end do
!!!
!!!  !allocate occupation number and spinsign
!!!  !fill them in normal way
!!!  allocate(orbs%occup(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%occup,'orbs%occup',subname)
!!!  allocate(orbs%spinsgn(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%spinsgn,'orbs%spinsgn',subname)
!!!  orbs%occup(1:orbs%norb*orbs%nkpts)=1.0_gp 
!!!  do ikpt=1,orbs%nkpts
!!!     do iorb=1,orbs%norbu
!!!        orbs%spinsgn(iorb+(ikpt-1)*orbs%norb)=1.0_gp
!!!     end do
!!!     do iorb=1,orbs%norbd
!!!        orbs%spinsgn(iorb+orbs%norbu+(ikpt-1)*orbs%norb)=-1.0_gp
!!!     end do
!!!  end do
!!!
!!!  !put a default value for the fermi energy
!!!  orbs%efermi = UNINITIALIZED(orbs%efermi)
!!!  !and also for the gap
!!!  orbs%HLgap = UNINITIALIZED(orbs%HLgap)
!!!
!!!  ! allocate inwhichlocreg
!!!
!!!  allocate(orbs%inwhichlocreg(orbs%norb*orbs%nkpts),stat=i_stat)
!!!  call memocc(i_stat,orbs%inwhichlocreg,'orbs%inwhichlocreg',subname)
!!!  ! default for inwhichlocreg
!!!  orbs%inwhichlocreg = 1
!!!
!!!  !nullify(orbs%inwhichlocregP)
!!!
!!!  !allocate the array which assign the k-point to processor in transposed version
!!!  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)
!!!
!!!  !initialize the starting point of the potential for each orbital (to be removed?)
!!!  allocate(orbs%ispot(orbs%norbp),stat=i_stat)
!!!  call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)
!!!
!!!
!!!  ! Define two new arrays:
!!!  ! - orbs%isorb_par is the same as orbs%isorb, but every process also knows
!!!  !   the reference orbital of each other process.
!!!  ! - orbs%onWhichMPI indicates on which MPI process a given orbital
!!!  !   is located.
!!!  allocate(orbs%isorb_par(0:nproc-1), stat=i_stat)
!!!  call memocc(i_stat, orbs%isorb_par, 'orbs%isorb_par', subname)
!!!  allocate(orbs%onWhichMPI(sum(orbs%norb_par(:,0))), stat=i_stat)
!!!  call memocc(i_stat, orbs%onWhichMPI, 'orbs%onWhichMPI', subname)
!!!  iiorb=0
!!!  orbs%isorb_par=0
!!!  do jproc=0,nproc-1
!!!      do iorb=1,orbs%norb_par(jproc,0)
!!!          iiorb=iiorb+1
!!!          orbs%onWhichMPI(iiorb)=jproc
!!!      end do
!!!      if(iproc==jproc) then
!!!          orbs%isorb_par(jproc)=orbs%isorb
!!!      end if
!!!  end do
!!!  call MPI_Initialized(mpiflag,ierr)
!!!  if(mpiflag /= 0 .and. nproc > 1) call mpiallred(orbs%isorb_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!!
!!!END SUBROUTINE orbitals_descriptors_forLinear


!> Routine which assigns to each processor the repartition of nobj*nkpts objects
subroutine kpts_to_procs_via_obj(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nproc !< No. of proc
  integer, intent(in) :: nkpts !< No. K points
  integer, intent(in) :: nobj  !< Object number (i.e. nvctr)
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_par !< iresult of the partition
  !local varaibles
  logical :: intrep
  integer :: jproc,ikpt,iobj,nobjp_max_kpt,nprocs_with_floor,jobj,nobjp
  integer :: jkpt,nproc_per_kpt,nproc_left,kproc,nkpt_per_proc,nkpts_left
  real(gp) :: robjp,rounding_ratio

  !decide the naive number of objects which should go to each processor.
  robjp=real(nobj,gp)*real(nkpts,gp)/real(nproc,gp)
  !print *,'hereweare',robjp,nobj   
  !maximum number of objects which has to go to each processor per k-point
  nobjp_max_kpt=ceiling(modulo(robjp-epsilon(1.0_gp),real(nobj,gp)))


!see the conditions for the integer repartition of k-points
  if (nobjp_max_kpt == nobj .or. (nobjp_max_kpt==1 .and. robjp < 1.0_gp)) then
     intrep=.true.
     rounding_ratio=0.0_gp
     nprocs_with_floor=0
  else
     intrep=.false.
     !the repartition is not obvious, some processors take nobj_max_kpt objects, others take the previous integer.
     !to understand how many, we round the percentage of processors which is given by
     rounding_ratio=(robjp-real(floor(robjp), gp))
     !then this is the number of processors which will take the floor
     nprocs_with_floor=ceiling((1.0_gp-rounding_ratio)*real(nproc,gp))!nproc-(nobj*nkpts-floor(robjp)*nproc)
     !print *,'rounding_ratio,nprocs_with_floor',rounding_ratio,nprocs_with_floor
     if (nprocs_with_floor > nproc) stop 'ERROR: should not happen'
     !if (nprocs_with_floor == nproc) nprocs_with_floor=nproc-1
  end if

  !start separating the objects for the repartition which is suggested by rounding_ratio and nprocs_with_floor
  nobj_par(0:nproc-1,1:nkpts)=0
  !integer repartition
  if (intrep) then 
     !strategy for the repartition
     if (nproc >= nkpts) then
        !decide in how many processors a single k-point can be partitioned
        nproc_per_kpt=max((nproc-1),1)/nkpts !this is the minimum
        !count how many processors are left that way
        !distribute the k-point among these
        nproc_left=nproc-nproc_per_kpt*nkpts
        ikpt=0
        jproc=0
        !print *,'here',nproc_left,nproc_per_kpt
        do kproc=0,nproc_left-1
           ikpt=ikpt+1
           if (ikpt > nkpts) stop 'ERROR: also this should not happen3'
           do iobj=0,nobj-1
              nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)+1
           end do
           jproc=jproc+nproc_per_kpt+1
        end do
        !print *,'debug'
        if ((nproc_per_kpt+1)*nproc_left < nproc) then
           do jproc=(nproc_per_kpt+1)*nproc_left,nproc-1,nproc_per_kpt
              ikpt=ikpt+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen3b'
              do iobj=0,nobj-1
                 nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)+1
              end do
           end do
        end if
        !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     else
        !decide in how many kpoints single processor can be partitioned
        nkpt_per_proc=max((nkpts-1),1)/nproc !this is the minimum
        !count how many k-points are left that way
        !distribute the processors among these
        nkpts_left=nkpts-nkpt_per_proc*nproc
        ikpt=1
        jproc=-1
        !print *,'hello',nkpts_left,nkpts_per_proc
        do jkpt=1,nkpts_left
           jproc=jproc+1
           if (jproc > nproc-1) stop 'ERROR: also this should not happen4'
           do iobj=0,(nobj)*(nkpt_per_proc+1)-1
              nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))+1
           end do
           ikpt=ikpt+nkpt_per_proc+1
        end do
        !print *,'ciao'
        if ((nkpt_per_proc+1)*nkpts_left < nkpts) then
           do ikpt=(nkpt_per_proc+1)*nkpts_left+1,nkpts,nkpt_per_proc
              jproc=jproc+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen4b'
              do iobj=0,(nobj)*(nkpt_per_proc)-1
                 nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))+1
              end do
           end do
        end if
           !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     end if
  else
     !non-integer repartition
     iobj=0
     ikpt=0
     do jproc=0,nproc-2 !leave the last processor at the end
        nobjp=floor(robjp)
        !respect the rounding ratio
        if (nproc-jproc > nprocs_with_floor) nobjp=nobjp+1
        !print *,'jproc,nobjp',jproc,nobjp,nkpts,nobj,nkpts*nobj,iobj,nprocs_with_floor
        do jobj=1,nobjp
           if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
           iobj=iobj+1
           if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen'
           nobj_par(jproc,ikpt)=nobj_par(jproc,ikpt)+1
        end do
     end do
     !in the last processor we put the objects which are lacking
     nobjp=nobj*nkpts-iobj
     do jobj=1,nobjp
        if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
        iobj=iobj+1
        !print *,'finished',jobj,nobjp,iobj,nobj*nkpts,jproc,ikpt
        if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen2'
        nobj_par(nproc-1,ikpt)=nobj_par(nproc-1,ikpt)+1
     end do
  end if
END SUBROUTINE kpts_to_procs_via_obj


subroutine components_kpt_distribution(nproc,nkpts,norb,nvctr,norb_par,nvctr_par)
  use module_base, only: gp, f_err_throw, to_zero,BIGDFT_RUNTIME_ERROR,&
       UNINITIALIZED
  implicit none
  !Arguments
  integer, intent(in) :: nproc,nkpts,nvctr,norb
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nvctr_par
  !local variables
  integer :: ikpt,jsproc,jeproc,kproc,icount,ivctr,jproc,numproc
  real(gp) :: strprc,endprc

  !for any of the k-points find the processors which have such k-point associated
  call to_zero(nproc*nkpts,nvctr_par(0,1))

  !Loop over each k point
  do ikpt=1,nkpts
     jsproc=UNINITIALIZED(1)
     jeproc=UNINITIALIZED(1)
     find_start: do jproc=0,nproc-1
        if(norb_par(jproc,ikpt) > 0) then 
           jsproc=jproc
           exit find_start
        end if
     end do find_start
     if (jsproc == UNINITIALIZED(1)) call f_err_throw('ERROR in kpt assignments',err_id=BIGDFT_RUNTIME_ERROR)
     if(norb_par(jsproc,ikpt) /= norb) then
        strprc=real(norb_par(jsproc,ikpt),gp)/real(norb,gp)     
     else
        strprc=1.0_gp
     end if
     if (ikpt < nkpts) then
        find_end: do jproc=jsproc,nproc-1
           if(norb_par(jproc,ikpt+1) > 0) then
              if (norb_par(jproc,ikpt)==0) then
                 jeproc=jproc-1
              else
                 jeproc=jproc
              end if
              exit find_end
           end if
        end do find_end
        if (jeproc == UNINITIALIZED(1)) call f_err_throw('ERROR in kpt assignments',err_id=BIGDFT_RUNTIME_ERROR)
     else
        jeproc=nproc-1
     end if
     if (jeproc /= jsproc) then
        endprc=real(norb_par(jeproc,ikpt),gp)/real(norb,gp)     
     else
        endprc=0.0_gp
     end if
     !if the number of processors is bigger than the number of orbitals this means 
     !that strprc and endprc are not correctly evaluated
     !evaluate the percentage on the number of components
     if (jeproc-jsproc+1 > norb) then
        strprc=1.0_gp/real(jeproc-jsproc+1,gp)
        endprc=strprc
     end if

     !assign the number of components which corresponds to the same orbital distribution
     numproc=jeproc-jsproc+1
     icount=0

     !print *,'kpoint',ikpt,jsproc,jeproc,strprc,endprc,ceiling(strprc*real(nvctr,gp)),nvctr
     !start filling the first processor
     ivctr=min(ceiling(strprc*real(nvctr,gp)-epsilon(1.0_gp)),nvctr)
     nvctr_par(jsproc,ikpt)=ivctr!min(ceiling(strprc*real(nvctr,gp)),nvctr)
     fill_array: do 
        if (ivctr==nvctr) exit fill_array
        icount=icount+1
        kproc=jsproc+modulo(icount,numproc)
        !put the floor of the components to the first processor
        if (strprc /= 1.0_gp .and. kproc==jsproc .and. &
             nvctr_par(kproc,ikpt)==ceiling(strprc*real(nvctr,gp)-epsilon(1.0_gp))) then
           !do nothing, skip away
        else
           nvctr_par(kproc,ikpt) = nvctr_par(kproc,ikpt)+1
           ivctr=ivctr+1
        end if
     end do fill_array
     !print '(a,i3,i3,i6,2(1pe25.17),i7,20i5)','here',ikpt,jsproc,jeproc,strprc,endprc,sum(nvctr_par(:,ikpt)),nvctr_par(:,ikpt)
     !print '(a,i3,i3,i6,2(1pe25.17),i7,20i5)','there',ikpt,jsproc,jeproc,strprc,endprc,sum(nvctr_par(:,ikpt)),norb_par(:,ikpt)
  end do

END SUBROUTINE components_kpt_distribution


!> Check the distribution of k points over the processors
subroutine check_kpt_distributions(nproc,nkpts,norb,ncomp,norb_par,ncomp_par,info,lub_orbs,lub_comps)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,norb,ncomp
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(in) :: ncomp_par
  integer, intent(inout) :: info
  integer, intent(out) :: lub_orbs,lub_comps
  !local variables
  character(len=*), parameter :: subname='check_kpt_distributions'
  logical :: notcompatible,couldbe
  integer :: ikpt,jproc,norbs,ncomps,kproc,ieproc,isproc,jkpt
  integer, dimension(:,:), allocatable :: load_unbalancing
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  if (info == 0) call print_distribution_schemes(nproc,nkpts,norb_par,ncomp_par)

  do ikpt=1,nkpts
     isproc=UNINITIALIZED(1)
     find_isproc : do kproc=0,nproc-1
        if (ncomp_par(kproc,ikpt) > 0) then
           isproc=kproc
           exit find_isproc
        end if
     end do find_isproc
     if (isproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): isproc cannot be found'
     ieproc=UNINITIALIZED(1)
     find_ieproc : do kproc=nproc-1,0,-1
        if (ncomp_par(kproc,ikpt) > 0) then
           ieproc=kproc
           exit find_ieproc
        end if
     end do find_ieproc
     if (ieproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): ieproc cannot be found'

     norbs=0
     ncomps=0
     do jproc=0,nproc-1
        !count the total number of components
        norbs=norbs+norb_par(jproc,ikpt)
        ncomps=ncomps+ncomp_par(jproc,ikpt)
        notcompatible=(ncomp_par(jproc,ikpt) == 0 .neqv. norb_par(jproc,ikpt) == 0) 
        !check whether there are only 0 orbitals
        if (notcompatible .and. norb_par(jproc,ikpt)==0) then
           !if the processor is the last one then there should not be other k-points on this processors
           couldbe=.false.
           if (jproc == ieproc) then
              couldbe=.true.
              do jkpt=ikpt+1,nkpts
                 couldbe=couldbe .and. (norb_par(jproc,jkpt) ==0 .and. ncomp_par(jproc,jkpt)==0)
              end do
           end if
           if ((isproc < jproc .and. jproc < ieproc) .or. couldbe) notcompatible=.false.
        end if
        if (notcompatible) then     
           if (info == 0) write(*,*)' ERROR: processor ', jproc,' kpt,',ikpt,&
                'have components and orbital distributions not compatible'
           info=1
           return
           !call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
        end if
     end do
     if (norb/=norbs .or. ncomps /= ncomp) then
        if (info == 0) write(*,*)' ERROR: kpt,',ikpt,&
             'has components or orbital distributions not correct'
        info=2
        return
        !call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
     end if
  end do

  load_unbalancing = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='load_unbalancing')

  do jproc=0,nproc-1
     load_unbalancing(jproc,:)=0
     do ikpt=1,nkpts
        load_unbalancing(jproc,1)=load_unbalancing(jproc,1)+norb_par(jproc,ikpt)
        load_unbalancing(jproc,2)=load_unbalancing(jproc,2)+ncomp_par(jproc,ikpt)
     end do
  end do

  !calculate the maximum load_unbalancing
  lub_orbs=0
  lub_comps=0
  do jproc=0,nproc-1
     do kproc=0,nproc-1
        lub_orbs=max(lub_orbs,load_unbalancing(jproc,1)-load_unbalancing(kproc,1))
        lub_comps=max(lub_comps,load_unbalancing(jproc,2)-load_unbalancing(kproc,2))
     end do
  end do

  if (info==0) write(*,*)' Kpoints Distribuitions are compatible, load unbalancings, orbs,comps:',lub_orbs,&
       '/',max(minval(load_unbalancing(:,1)),1),lub_comps,'/',minval(load_unbalancing(:,2))
  info=0
  call f_free(load_unbalancing)


END SUBROUTINE check_kpt_distributions

!>routine which associates to any of the processor a given number of objects
!! depending of the number of processors and k-points
subroutine parallel_repartition_with_kpoints(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nkpts,nobj,nproc
  integer, dimension(0:nproc-1), intent(out) :: nobj_par
  !local variables
  integer :: n_i,n_ip,rs_i,N_a,N_b,N_c,ikpt,jproc,i,ntmp
!!$  real(gp) :: rtmp

  ! Strategy to divide between k points.
  ! There is an nproc length to divide into orbs%nkpts segments.
  ! Segment (ikpt - 1) expand in 0 <= r_i < r_ip <= nproc.
  ! where r_i and r_ip are real values. There are two possibilities:
  !  - We can write r_i <= n_i <= n_ip <= r_ip with n_i and n_ip integers ;
  !  - or r_i <= n_i and n_ip <= r_ip and n_i = n_ip + 1.
  ! For both cases, we can divide nobj into the partition (real values):
  !  - N_a = (n_i - r_i)*nobj*nkpts/nproc (the initial part);
  !  - N_b = max((n_ip - n_i)*nobj*nkpts / nproc, 0) (the naive part, the only one if nkpts is a multiple of nproc);
  !  - N_c = (r_ip - n_ip) * nobj * orbs%nkpts / nproc (the final part);
  ! Before going to integer values, we have r_i = (ikpt - 1) * nproc / orbs%nkpts (the naive division)
  ! and r_ip = (ikpt) * nproc / orbs%nkpts (the segment endpoint)
  ! So N_a and N_b can be simplified and written instead:
  !  - N_a = int(nobj * (n_i * orbs%nkpts - (ikpt - 1) * nproc) / nproc);
  !  - N_c = int(nobj * ((ikpt) * nproc - n_ip * orbs%nkpts) / nproc)
  !  - N_b = nobj - N_a - N_c 
  ! After, if N_a > 0, we put this quantity to proc n_i - 1, if N_c > 0
  ! we put its quantity to proc n_ip ; and finally N_b is distributed
  ! among [n_i;n_ip[ procs.

  nobj_par(:)=0
  do ikpt=1,nkpts
     ! Calculation of n_i and n_ip, rs_i = r_i * orbs%nkpts to avoid rounding.
     rs_i=(ikpt-1)*nproc !integer variable for rounding purposes

     if (mod(rs_i,nkpts) == 0) then
        n_i=rs_i/nkpts 
     else
        n_i=rs_i/nkpts+1
     end if

     rs_i=ikpt*nproc
     n_ip=rs_i/nkpts
!!$     print *,'ikpt,ni,nip',ikpt,n_i,n_ip
     ! Calculation of N_a, N_b and N_c from given n_i and n_ip.
     if (n_ip >= n_i) then
        ntmp = (n_i*nkpts-(ikpt-1)*nproc) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_a = ntmp / nproc
        else
           N_a = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if
!!$        ntmp=n_i*nkpts-(ikpt-1)*nproc
!!$        rtmp=real(nobj,gp)/real(nproc,gp)
!!$        rtmp=rtmp*real(ntmp,gp)
!!$        N_a=floor(rtmp)
!!$        if (iproc == 0) print *,'ikpts,rtmp',ikpt,rtmp
        ntmp = (ikpt*nproc-n_ip*nkpts) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_c = ntmp / nproc
        else
           N_c = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if

!!$        ntmp=ikpt*nproc-n_ip*nkpts
!!$        rtmp=real(nobj,gp)/real(nproc,gp)
!!$        rtmp=rtmp*real(ntmp,gp)
!!$        N_c=ceiling(rtmp)
!!$        if (iproc == 0) print *,'ikpts,rtmp2',ikpt,rtmp,N_a,N_c
        !the corrections above are to avoid the 32 bit integer overflow
        !N_a=nint(real(nobj*(n_i*nkpts-(ikpt-1)*nproc),gp)/real(nproc,gp))
        !N_c=nint(real(nobj*(ikpt*nproc-n_ip*nkpts),gp)/real(nproc,gp))
     else
        N_c=nobj/2
        N_a=nobj-N_c
     end if
     N_b=nobj-N_a-N_c
     if (N_b == -1) then
        N_c = N_c - 1
        N_b = 0
     end if
!!$     write(*,*) ikpt, N_a, N_b, N_c
     if (nkpts > 1 .and. N_b < n_ip - n_i) stop 'ERROR:parallel_repartion_with_kpoints'
     !assign to procs the objects.
     if (N_a>0) nobj_par(n_i-1)=nobj_par(n_i-1)+N_a
     if (N_b>0) then
        do i=0,N_b-1
           jproc=n_i+mod(i,n_ip-n_i)
           nobj_par(jproc)=nobj_par(jproc)+1
        end do
     end if
     if (N_c>0) nobj_par(n_ip)=nobj_par(n_ip)+N_c
  end do
END SUBROUTINE parallel_repartition_with_kpoints


subroutine parallel_repartition_per_kpoints(iproc,nproc,nkpts,nobj,nobj_par,&
     nkptsp,mykpts,nobj_pkpt)
  implicit none
  integer, intent(in) :: iproc,nproc,nkpts,nobj
  integer, dimension(0:nproc-1), intent(in) :: nobj_par
  integer, intent(out) :: nkptsp
  integer, dimension(nkpts), intent(out) :: mykpts
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_pkpt
  !local variables
  integer :: ikpts,jproc,jobj,norb_tot,iorbp

  !initialise the array
  do ikpts=1,nkpts
     do jproc=0,nproc-1
        nobj_pkpt(jproc,ikpts)=0 
     end do
  end do

  !assign the k-point, counting one object after each other
  jobj=1
  ikpts=1
  !print *,'here',nobj_par(:)
  do jproc=0,nproc-1
     do iorbp=1,nobj_par(jproc)
        nobj_pkpt(jproc,ikpts)=nobj_pkpt(jproc,ikpts)+1
        if (mod(jobj,nobj)==0) then
           ikpts=ikpts+1
        end if
        jobj=jobj+1
     end do
  end do
  !some checks
  if (nobj /= 0) then
     !check the distribution
     do ikpts=1,nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',nobj_pkpt(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+nobj_pkpt(jproc,ikpts)
        end do
        if(norb_tot /= nobj) then
           write(*,*)'ERROR: partition of objects incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if

  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  nkptsp=0
  do ikpts=1,nkpts
     if (nobj_pkpt(iproc,ikpts) /= 0) then
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do

END SUBROUTINE parallel_repartition_per_kpoints

subroutine pawpatch_from_file( filename, atoms,ityp, paw_tot_l, &
     paw_tot_q, paw_tot_coefficients, paw_tot_matrices, &
     storeit)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  type(atoms_data), intent(inout) :: atoms
  integer , intent(IN):: ityp 
  integer , intent(INOUT):: paw_tot_l, paw_tot_q, paw_tot_coefficients, paw_tot_matrices
  logical, intent(in) :: storeit

!! local variables  
  character(len=*), parameter :: subname='pawpatch_from_file'
  integer :: npawl, ipawl, paw_l
  integer :: paw_nofgaussians, paw_nofchannels, il, ierror, ig
  real(gp) :: paw_greal, paw_gimag, paw_ccoeff, paw_scoeff, dumpaw
  character(len=100) :: string

  !parameters for abscalc-paw

  if(.not. storeit) then
     !if(ityp == 1) then !this implies that the PSP are all present
     if (.not. associated(atoms%paw_NofL)) then
        atoms%paw_NofL = f_malloc_ptr(atoms%astruct%ntypes+ndebug,id='atoms%paw_NofL')
     end if
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=trim(filename),status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) ': Failed to open the PAWpatch file "',&
             trim(filename),'"'
        stop
     end if
     
     !! search for paw_patch informations
     
     atoms%paw_NofL(ityp)=0
     do while(.true.)
        read(11,'(a)',iostat=ierror, END=110)  string
        if ( trim(string).eq. 'PAWPATCH') then
           exit
        endif
     end do
     !! explain_paw_psp_terms_in_atom_data()
     
     read(11,*) npawl
     
     atoms%paw_NofL(ityp) = npawl
     
     paw_tot_l = paw_tot_l + npawl
     do ipawl=1,npawl
        read(11,*) paw_l
        read(11,*) paw_greal
        read(11,*) paw_nofgaussians
        read(11,*) paw_nofchannels
        read(11,*)  string  !!  follow 300 PAW_Gimag factors
        paw_tot_q = paw_tot_q+paw_nofgaussians
        paw_tot_coefficients = paw_tot_coefficients + paw_nofchannels*paw_nofgaussians*2
        paw_tot_matrices=paw_tot_matrices+paw_nofchannels**2
        
        
        do ig=1, paw_nofgaussians
           read(11,*)  paw_gimag
        enddo
        read(11,*)  string  !!  !!  follow for each of the 7 channels 300 (cos_factor, sin_factor)  pairs
        do il=1, paw_nofchannels
           do ig=1, paw_nofgaussians
              read(11,*)  paw_ccoeff, paw_scoeff
           enddo
        enddo
        read(11,*)  string  !! pawpatch matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
        read(11,*)  string  !! S  matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
        
        read(11,*)  string  !! Sm1  matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
     enddo
110  close(11)

  else
     if(ityp.eq.1) then
        atoms%paw_l   = f_malloc_ptr(paw_tot_l,id='atoms%paw_l  ')
        atoms%paw_nofchannels   = f_malloc_ptr(paw_tot_l,id='atoms%paw_nofchannels  ')
        atoms%paw_nofgaussians   = f_malloc_ptr(paw_tot_l,id='atoms%paw_nofgaussians  ')
        atoms%paw_Greal   = f_malloc_ptr(paw_tot_l,id='atoms%paw_Greal  ')
        atoms%paw_Gimag  = f_malloc_ptr(paw_tot_q   ,id='atoms%paw_Gimag ')
        atoms%paw_Gcoeffs  = f_malloc_ptr(paw_tot_coefficients  ,id='atoms%paw_Gcoeffs ')
        atoms%paw_H_matrices = f_malloc_ptr(paw_tot_matrices,id='atoms%paw_H_matrices')
        atoms%paw_S_matrices  = f_malloc_ptr(paw_tot_matrices  ,id='atoms%paw_S_matrices ')
        atoms%paw_Sm1_matrices  = f_malloc_ptr(paw_tot_matrices  ,id='atoms%paw_Sm1_matrices ')
        
        
        paw_tot_l=0
        paw_tot_q = 0
        paw_tot_coefficients = 0
        paw_tot_matrices= 0
     endif
     
     if( atoms%paw_NofL(ityp).gt.0) then

        open(unit=11,file=trim(filename),status='old')

        do while(.true.)
           read(11,'(a)', END=220)  string
           if ( trim(string).eq. 'PAWPATCH') then
              exit
           endif
        end do
220     continue
        if(trim(string) .ne. 'PAWPATCH') then
           print *, "paw section not found re-reading  file ", filename
           close(11)
           stop
        end if
        
        read(11,*) npawl
        atoms%paw_NofL(ityp) = npawl
        
        
        do ipawl=1,npawl
           !! explain_paw_psp_terms_in_atom_data()
           read(11,*) atoms%paw_l(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_greal(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_nofgaussians(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_nofchannels(paw_tot_l+ipawl  )
           paw_nofchannels = atoms%paw_nofchannels(paw_tot_l+ipawl  )
           paw_nofgaussians = atoms%paw_nofgaussians(paw_tot_l+ipawl  )
           read(11,'(a)')  string  !!  follow  paw_nofchannels PAW_Gimag factors
           
           do ig=1, paw_nofgaussians
              read(11,*)  atoms%paw_Gimag(paw_tot_q + ig )
           enddo
           read(11,'(a)')  string  !!  !!  follow for each of the Npaw channels  nofgaussians (cos_factor, sin_factor)  pairs
           !!print *, string, " reading  " , filename
           
           do il=1, paw_nofchannels
              do ig=1, paw_nofgaussians
                 read(11,*)  atoms%paw_Gcoeffs( paw_tot_coefficients + 2*(il-1)*paw_nofgaussians+2*ig-1) , &
                      atoms%paw_Gcoeffs( paw_tot_coefficients + 2*(il-1)*paw_nofgaussians+2*ig)
              enddo
           enddo
           read(11,'(a)')  string  !! pawpatch matrix
           !!print *, string, " reading  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw 
                 atoms%paw_H_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
              end do
           enddo
           read(11,'(a)')  string  !! S  matrix
           !!print *, string, " reading >>>>>  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw
                 atoms%paw_S_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
              end do
           enddo
           read(11,'(a)')  string  !! Sm1  matrix
           !!print *, string, " reading  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw
                 atoms%paw_Sm1_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
                 !!print *, dumpaw
              end do
           enddo
           paw_tot_q = paw_tot_q+paw_nofgaussians
           paw_tot_coefficients = paw_tot_coefficients + paw_nofchannels*paw_nofgaussians*2
           paw_tot_matrices=paw_tot_matrices+paw_nofchannels**2
        end do
        paw_tot_l = paw_tot_l + npawl
     end if
     close(11)
  endif
end subroutine pawpatch_from_file
 

subroutine system_signaling(iproc, signaling, gmainloop, KSwfn, tmb, energs, denspot, optloop, &
       & ntypes, radii_cf, crmult, frmult)
  use module_defs, only: gp, UNINITIALIZED
  use module_types
  implicit none
  integer, intent(in) :: iproc, ntypes
  logical, intent(in) :: signaling
  double precision, intent(in) :: gmainloop
  type(DFT_wavefunction), intent(inout) :: KSwfn, tmb
  type(DFT_local_fields), intent(inout) :: denspot
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(inout) :: energs
  real(gp), dimension(ntypes,3), intent(in) :: radii_cf
  real(gp), intent(in) :: crmult, frmult

  if (signaling) then
     ! Only iproc 0 has the C wrappers.
     if (iproc == 0) then
        call wf_new_wrapper(KSwfn%c_obj, KSwfn, 0)
        call wf_copy_from_fortran(KSwfn%c_obj, radii_cf, crmult, frmult)
        call wf_new_wrapper(tmb%c_obj, tmb, 1)
        call wf_copy_from_fortran(tmb%c_obj, radii_cf, crmult, frmult)
        call bigdft_signals_add_wf(gmainloop, KSwfn%c_obj, tmb%c_obj)
        call energs_new_wrapper(energs%c_obj, energs)
        call bigdft_signals_add_energs(gmainloop, energs%c_obj)
        call localfields_new_wrapper(denspot%c_obj, denspot)
        call bigdft_signals_add_denspot(gmainloop, denspot%c_obj)
        call optloop_new_wrapper(optLoop%c_obj, optLoop)
        call bigdft_signals_add_optloop(gmainloop, optLoop%c_obj)
     else
        KSwfn%c_obj   = UNINITIALIZED(KSwfn%c_obj)
        tmb%c_obj     = UNINITIALIZED(tmb%c_obj)
        denspot%c_obj = UNINITIALIZED(denspot%c_obj)
        optloop%c_obj = UNINITIALIZED(optloop%c_obj)
     end if
  else
     KSwfn%c_obj  = 0
     tmb%c_obj    = 0
  end if
END SUBROUTINE system_signaling

!> @file
!!  Routines related to system properties
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Initialize the objects needed for the computation: basis sets, allocate required space
subroutine system_initialization(iproc,nproc,inputpsi,input_wf_format,in,atoms,rxyz,&
     orbs,lnpsidim_orbs,lnpsidim_comp,lorbs,Lzd,Lzd_lin,denspot,nlpspd,comms,shift,proj,radii_cf,&
     inwhichlocreg_old, onwhichatom_old)
  use module_base
  use module_types
  use module_interfaces, fake_name => system_initialization
  use module_xc
  use gaussians, only: gaussian_basis
  use vdwcorrection
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc 
  integer, intent(out) :: inputpsi, input_wf_format, lnpsidim_orbs, lnpsidim_comp
  type(input_variables), intent(in) :: in 
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  type(orbitals_data), intent(inout) :: orbs, lorbs
  type(local_zone_descriptors), intent(inout) :: Lzd, Lzd_lin
  type(DFT_local_fields), intent(out) :: denspot
  type(nonlocal_psp_descriptors), intent(out) :: nlpspd
  type(communications_arrays), intent(out) :: comms
  real(gp), dimension(3), intent(out) :: shift  !< shift on the initial positions
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
  real(wp), dimension(:), pointer :: proj
  integer,dimension(:),pointer,optional:: inwhichlocreg_old, onwhichatom_old
  !local variables
  character(len = *), parameter :: subname = "system_initialization"
  integer :: nB,nKB,nMB,ii,iat,iorb,iatyp!,i_stat,i_all
  real(gp) :: peakmem
  real(gp), dimension(3) :: h_input
  logical:: present_inwhichlocreg_old, present_onwhichatom_old
  !Note proj_G should be filled for PAW:
  type(gaussian_basis),dimension(atoms%nat)::proj_G

  !nullify dummy variables only used for PAW:
  do iatyp=1,atoms%ntypes
    call nullify_gaussian_basis(proj_G(iatyp))
  end do


  ! Dump XC functionals (done now in output.f90)
  !if (iproc == 0) call xc_dump()

!!$  if (iproc==0) then
!!$     write( *,'(1x,a)')&
!!$          &   '------------------------------------------------------------------ System Properties'
!!$  end if
  call read_radii_variables(atoms, radii_cf, in%crmult, in%frmult, in%projrad)
  if (iproc == 0) call print_atomic_variables(atoms, radii_cf, max(in%hx,in%hy,in%hz), in%ixc)

  Lzd=default_lzd()

  !grid spacings of the zone descriptors (not correct, the set is done by system size)
  h_input=(/ in%hx, in%hy, in%hz /)
  call lzd_set_hgrids(Lzd,h_input) 

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,&
       Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
       Lzd%Glr,shift)

  ! A message about dispersion forces.
  call vdwcorrection_initializeparams(in%ixc, in%dispersion)
  if (iproc == 0) call vdwcorrection_warnings(atoms, in)

  call initialize_DFT_local_fields(denspot)

  !here the initialization of dpbox can be set up
  call dpbox_set(denspot%dpbox,Lzd,iproc,nproc,bigdft_mpi%mpi_comm,in,atoms%geocode)

  ! Create the Poisson solver kernels.
  call system_initKernels(.true.,iproc,nproc,atoms%geocode,in,denspot)
  call system_createKernels(denspot, (verbose > 1))

  ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
  call createWavefunctionsDescriptors(iproc,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),atoms,&
       rxyz,radii_cf,in%crmult,in%frmult,Lzd%Glr)

  ! Create global orbs data structure.
  call read_orbital_variables(iproc,nproc,(iproc == 0),in,atoms,orbs)
  ! Create linear orbs data structure.
  if (in%inputpsiId == INPUT_PSI_LINEAR_AO .or. in%inputpsiId == INPUT_PSI_DISK_LINEAR &
      .or. in%inputpsiId == INPUT_PSI_MEMORY_LINEAR) then
     call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms, rxyz, lorbs)

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
         do iat=1, atoms%nat
            do iorb=1,lorbs%norb
               if(iat ==  lorbs%onwhichatom(iorb)) then
                  ii = ii + 1
                  lorbs%inwhichlocreg(iorb)= ii
                  !lorbs%onwhichmpi(iorb) = ii-1
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

  !allocate communications arrays (allocate it before Projectors because of the definition
  !of iskpts and nkptsp)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)  
  if (in%inputpsiId == INPUT_PSI_LINEAR_AO .or. in%inputpsiId == INPUT_PSI_DISK_LINEAR &
      .or. in%inputpsiId == INPUT_PSI_MEMORY_LINEAR) then
     if(iproc==0) call print_orbital_distribution(iproc, nproc, lorbs)
  end if

  if (iproc == 0) then
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

  inputpsi = in%inputPsiId

  call input_check_psi_id(inputpsi, input_wf_format, in%dir_output, orbs, lorbs, iproc, nproc)

  ! See if linear scaling should be activated and build the correct Lzd 
  call check_linear_and_create_Lzd(iproc,nproc,in%linear,Lzd,atoms,orbs,in%nspin,rxyz)
  lzd_lin=default_lzd()
  call nullify_local_zone_descriptors(lzd_lin)
  lzd_lin%nlr = 0
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
     .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     call copy_locreg_descriptors(Lzd%Glr, lzd_lin%glr, subname)
     call lzd_set_hgrids(lzd_lin, Lzd%hgrids)
     if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
        call lzd_init_llr(iproc, nproc, in, atoms, rxyz, lorbs, lzd_lin)
     else
        call initialize_linear_from_file(iproc,nproc,trim(in%dir_output)//'minBasis',&
             input_wf_format,lzd_lin,lorbs,atoms,rxyz)
        !what to do with derivatives?
     end if
     call update_wavefunctions_size(lzd_lin,lnpsidim_orbs,lnpsidim_comp,lorbs,iproc,nproc)
  end if

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call createProjectorsArrays(iproc,Lzd%Glr,rxyz,atoms,orbs,&
       radii_cf,in%frmult,in%frmult,Lzd%hgrids(1),Lzd%hgrids(2),&
       Lzd%hgrids(3),nlpspd,proj_G,proj)
  !calculate the partitioning of the orbitals between the different processors
  !memory estimation, to be rebuilt in a more modular way
  if (iproc==0 .and. verbose > 0) then
     call MemoryEstimator(nproc,in%idsx,Lzd%Glr,&
          atoms%nat,orbs%norb,orbs%nspinor,orbs%nkpts,nlpspd%nprojel,&
          in%nspin,in%itrpmax,in%iscf,peakmem)
  end if
  
  !here dpbox can be put as input
  call density_descriptors(iproc,nproc,in%nspin,in%crmult,in%frmult,atoms,&
       denspot%dpbox,in%rho_commun,rxyz,radii_cf,denspot%rhod)

  !allocate the arrays.
  call allocateRhoPot(iproc,Lzd%Glr,in%nspin,atoms,rxyz,denspot)

  !calculate the irreductible zone for this region, if necessary.
  call symmetry_set_irreductible_zone(atoms%sym,atoms%geocode, &
       & Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i, in%nspin)

  !check the communication distribution
  if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
     .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
      call check_communications(iproc,nproc,orbs,Lzd%Glr,comms)
  else
      ! Do not call check_communication, since the value of orbs%npsidim_orbs is wrong
      if(iproc==0) call yaml_warning('Do not call check_communications in the linear scaling version!')
      !if(iproc==0) write(*,*) 'WARNING: do not call check_communications in the linear scaling version!'
  end if

  !---end of system definition routine
END SUBROUTINE system_initialization


subroutine system_initKernels(verb, iproc, nproc, geocode, in, denspot)
  use module_types
  use module_xc
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  implicit none
  logical, intent(in) :: verb
  integer, intent(in) :: iproc, nproc
  character, intent(in) :: geocode
  type(input_variables), intent(in) :: in
  type(DFT_local_fields), intent(inout) :: denspot

  integer, parameter :: ndegree_ip = 16

  denspot%pkernel=pkernel_init(verb, iproc,nproc,in%matacc%PSolver_igpu,&
       geocode,denspot%dpbox%ndims,denspot%dpbox%hgrids,ndegree_ip,mpi_env=denspot%dpbox%mpi_env)

  !create the sequential kernel if the exctX parallelisation scheme requires it
  if ((xc_exctXfac() /= 0.0_gp .and. in%exctxpar=='OP2P' .or. in%SIC%alpha /= 0.0_gp)&
       .and. denspot%dpbox%mpi_env%nproc > 1) then
     !the communicator of this kernel is bigdft_mpi%mpi_comm
     denspot%pkernelseq=pkernel_init(iproc==0 .and. verb,0,1,in%matacc%PSolver_igpu,&
          geocode,denspot%dpbox%ndims,denspot%dpbox%hgrids,ndegree_ip)
  else 
     denspot%pkernelseq = denspot%pkernel
  end if
END SUBROUTINE system_initKernels

subroutine system_createKernels(denspot, verb)
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
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
  !local variables
  !n(c) character(len=*), parameter :: subname='system_properties'

  call read_radii_variables(atoms, radii_cf, in%crmult, in%frmult, in%projrad)
!!$  call read_atomic_variables(atoms, trim(in%file_igpop),in%nspin)
  if (iproc == 0) call print_atomic_variables(atoms, radii_cf, max(in%hx,in%hy,in%hz), in%ixc)
  call read_orbital_variables(iproc,nproc,(iproc == 0),in,atoms,orbs)
END SUBROUTINE system_properties


!> Check for the need of a core density and fill the rhocore array which
!! should be passed at the rhocore pointer
subroutine calculate_rhocore(iproc,at,d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,i3s,n3d,i3xcsh,n3p
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  type(grid_dimensions), intent(in) :: d
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(:,:,:,:), pointer :: rhocore
  !local variables
  character(len=*), parameter :: subname='calculate_rhocore'
  integer :: ityp,iat,i_stat,j3,i1,i2,ierr!,ind
  real(wp) :: tt
  real(gp) :: rx,ry,rz,rloc,cutoff
  

  !check for the need of a nonlinear core correction
!!$  donlcc=.false.
!!$  chk_nlcc: do ityp=1,at%ntypes
!!$     filename = 'nlcc.'//at%atomnames(ityp)
!!$
!!$     inquire(file=filename,exist=exists)
!!$     if (exists) then
!!$        donlcc=.true.
!!$        exit chk_nlcc
!!$     end if
!!$  end do chk_nlcc

  if (at%donlcc) then
     !allocate pointer rhocore
     allocate(rhocore(d%n1i,d%n2i,n3d,10+ndebug),stat=i_stat)
     call memocc(i_stat,rhocore,'rhocore',subname)
     !initalise it 
     call to_zero(d%n1i*d%n2i*n3d*10,rhocore(1,1,1,1))
     !perform the loop on any of the atoms which have this feature
     do iat=1,at%nat
        ityp=at%iatype(iat)
!!$        filename = 'nlcc.'//at%atomnames(ityp)
!!$        inquire(file=filename,exist=exists)
!!$        if (exists) then
        if (at%nlcc_ngv(ityp)/=UNINITIALIZED(1) .or.&
             at%nlcc_ngc(ityp)/=UNINITIALIZED(1) ) then
           if (iproc == 0) call yaml_map('NLCC, Calculate core density for atom:',trim(at%atomnames(ityp)))
           !if (iproc == 0) write(*,'(1x,a)',advance='no') 'NLCC: calculate core density for atom: '// trim(at%atomnames(ityp))//';'
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           cutoff=10.d0*rloc

           call calc_rhocore_iat(iproc,at,ityp,rx,ry,rz,cutoff,hxh,hyh,hzh,&
                d%n1,d%n2,d%n3,d%n1i,d%n2i,d%n3i,i3s,n3d,rhocore)

           !if (iproc == 0) write(*,'(1x,a)')'done.'
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

     call mpiallred(tt,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     tt=tt*hxh*hyh*hzh
     if (iproc == 0) call yaml_map('Total core charge on the grid (To be compared with analytic one)', tt,fmt='(f15.7)')
     !if (iproc == 0) write(*,'(1x,a,f15.7)') 'Total core charge on the grid (To be compared with analytic one): ',tt

  else
     !No NLCC needed, nullify the pointer 
     nullify(rhocore)
  end if

END SUBROUTINE calculate_rhocore


!> Initialization of the atoms_data values especially nullification of the pointers
subroutine init_atomic_values(verb, atoms, ixc)
  use module_base
  use module_types
  implicit none
  
  integer, intent(in) :: ixc
  logical, intent(in) :: verb
  type(atoms_data), intent(inout) :: atoms

  !local variables
  character(len=*), parameter :: subname='init_atomic_values'
  real(gp), dimension(3) :: radii_cf
  integer :: nlcc_dim, ityp, ig, j, ngv, ngc, i_stat,i_all,ierr
  integer :: paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices
  logical :: exists, read_radii,exist_all
  character(len=27) :: filename
  
  ! Read values from pseudo files.
  nlcc_dim=0
  atoms%donlcc=.false.
  paw_tot_l=0
  paw_tot_q=0
  paw_tot_coefficients=0
  paw_tot_matrices=0

  !True if there are atoms
  exist_all=(atoms%ntypes > 0)
  !@todo: eliminate the pawpatch from psppar
  nullify(atoms%paw_NofL)
  do ityp=1,atoms%ntypes
     filename = 'psppar.'//atoms%atomnames(ityp)
     call psp_from_file(filename, atoms%nzatom(ityp), atoms%nelpsp(ityp), &
           & atoms%npspcode(ityp), atoms%ixcpsp(ityp), atoms%psppar(:,:,ityp), &
           & radii_cf, read_radii, exists)
     !To eliminate the runtime warning due to the copy of the array (TD)
     atoms%radii_cf(ityp,:)=radii_cf(:)

     if (exists) then
        !! first time just for dimension ( storeit = . false.)
        call pawpatch_from_file( filename, atoms,ityp,&
             paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .false.)
     end if
     exist_all=exist_all .and. exists

     if (.not. read_radii) atoms%radii_cf(ityp, :) = UNINITIALIZED(1.0_gp)
     if (.not. exists) then
        atoms%ixcpsp(ityp) = ixc
        call psp_from_data(atoms%atomnames(ityp), atoms%nzatom(ityp), &
             & atoms%nelpsp(ityp), atoms%npspcode(ityp), atoms%ixcpsp(ityp), &
             & atoms%psppar(:,:,ityp), exists)
        if (.not. exists) then
           call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
           if (verb) write(*,'(1x,5a)')&
                'ERROR: The pseudopotential parameter file "',trim(filename),&
                '" is lacking, and no registered pseudo found for "', &
                & trim(atoms%atomnames(ityp)), '", exiting...'
           stop
        end if
     end if
     filename ='nlcc.'//atoms%atomnames(ityp)
     call nlcc_dim_from_file(filename, atoms%nlcc_ngv(ityp), &
          atoms%nlcc_ngc(ityp), nlcc_dim, exists)
     atoms%donlcc = (atoms%donlcc .or. exists)
  end do
  !deallocate the paw_array if not all the atoms are present
  if (.not. exist_all .and. associated(atoms%paw_NofL)) then
     i_all=-product(shape(atoms%paw_NofL ))*kind(atoms%paw_NofL )
     deallocate(atoms%paw_NofL,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%paw_NofL',subname)
     nullify(atoms%paw_NofL)
  end if

  if (exist_all) then
     do ityp=1,atoms%ntypes
        filename = 'psppar.'//atoms%atomnames(ityp)
        !! second time allocate and then store
        call pawpatch_from_file( filename, atoms,ityp,&
             paw_tot_l,   paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .true.)
     end do
  else
     nullify(atoms%paw_l,atoms%paw_NofL,atoms%paw_nofchannels)
     nullify(atoms%paw_nofgaussians,atoms%paw_Greal,atoms%paw_Gimag)
     nullify(atoms%paw_Gcoeffs,atoms%paw_H_matrices,atoms%paw_S_matrices,atoms%paw_Sm1_matrices)
  end if
  !process the nlcc parameters if present 
  !(allocation is performed also with zero size)
  allocate(atoms%nlccpar(0:4,max(nlcc_dim,1)+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlccpar,'atoms%nlccpar',subname)
  !start again the file inspection to fill nlcc parameters
  if (atoms%donlcc) then
     nlcc_dim=0
     fill_nlcc: do ityp=1,atoms%ntypes
        filename = 'nlcc.'//atoms%atomnames(ityp)
        inquire(file=filename,exist=exists)
        if (exists) then
           !read the values of the gaussian for valence and core densities
           open(unit=79,file=filename,status='unknown')
           read(79,*)ngv
           do ig=1,(ngv*(ngv+1))/2
              nlcc_dim=nlcc_dim+1
              read(79,*)(atoms%nlccpar(j,nlcc_dim),j=0,4)!rhovxp(ig),(rhovc(ig,j),j=1,4)
           end do
           read(79,*)ngc
           do ig=1,(ngc*(ngc+1))/2
              nlcc_dim=nlcc_dim+1
              read(79,*)(atoms%nlccpar(j,nlcc_dim),j=0,4)!rhocxp(ig),(rhocc(ig,j),j=1,4)
           end do
           close(unit=79)
        end if
     end do fill_nlcc
  end if

END SUBROUTINE init_atomic_values


subroutine psp_from_file(filename, nzatom, nelpsp, npspcode, &
     & ixcpsp, psppar, radii_cf, read_radii, exists)
  use module_base
  implicit none
  
  character(len = *), intent(in) :: filename
  integer, intent(out) :: nzatom, nelpsp, npspcode, ixcpsp
  real(gp), intent(out) :: psppar(0:4,0:6), radii_cf(3)
  logical, intent(out) :: read_radii, exists

  integer :: ierror, ierror1, i, j, nn, nlterms, nprl, l
  character(len=100) :: line

  read_radii = .false.
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
  read(11,*) nzatom, nelpsp
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
        read(11,*) !k coefficients, not used for the moment (no spin-orbit coupling)
     enddo
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
           read(11,*) !k coefficients, not used
        end do
     end do prjloop
  else
     !if (iproc == 0) then
     write(*,'(1x,a,a)') trim(filename),&
          'unrecognized pspcode: only GTH, HGH & HGH-K pseudos (ABINIT format)'
     !end if
     stop
  end if

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
  end if
  close(11)

  read_radii = (ierror == 0)
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
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: crmult, frmult, projrad
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf

  integer, parameter :: nelecmax=32,nmax=6,lmax=4
  character(len=2) :: symbol
  integer :: i,ityp,mxpl,mxchg,nsccode
  real(gp) :: rcov,rprb,ehomo,radfine,amu,maxrad
  real(kind=8), dimension(nmax,0:lmax-1) :: neleconf

  do ityp=1,atoms%ntypes
     !see whether the atom is semicore or not
     !and consider the ground state electronic configuration
     call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
          neleconf,nsccode,mxpl,mxchg,amu)
     if (atoms%radii_cf(ityp, 1) == UNINITIALIZED(1.0_gp)) then
        !assigning the radii by calculating physical parameters
        radii_cf(ityp,1)=1._gp/sqrt(abs(2._gp*ehomo))
        radfine=100._gp
        do i=0,4
           if (atoms%psppar(i,0,ityp)/=0._gp) then
              radfine=min(radfine,atoms%psppar(i,0,ityp))
           end if
        end do
        radii_cf(ityp,2)=radfine
        radii_cf(ityp,3)=radfine
     else
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

subroutine read_orbital_variables(iproc,nproc,verb,in,atoms,orbs)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(in) :: iproc,nproc
  logical, intent(in) :: verb
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(inout) :: orbs
  !local variables
  character(len=*), parameter :: subname='read_orbital_variables'
  integer, parameter :: nelecmax=32,nmax=6,lmax=4,noccmax=2
  logical :: exists
  integer :: iat,iunit,norb,norbu,norbd,nspinor,jpst,norbme,norbyou,jproc,ikpts
  integer :: norbuempty,norbdempty,nelec
  integer :: nt,ntu,ntd,ityp,ierror,ispinsum
  integer :: ispol,ichg,ichgsum,norbe,norbat,nspin
  integer, dimension(lmax) :: nl
  real(gp), dimension(noccmax,lmax) :: occup
  character(len=100) :: radical

  !calculate number of electrons and orbitals
  ! Number of electrons and number of semicore atoms
  nelec=0
  do iat=1,atoms%nat
     ityp=atoms%iatype(iat)
     nelec=nelec+atoms%nelpsp(ityp)
  enddo
  nelec=nelec-in%ncharge

  if(nelec < 0.0 ) then
    if(iproc==0) write(*,*)'ERROR: Number of electrons is negative:',nelec,'.'
    if(iproc==0) write(*,*)'FIX: decrease charge of system.'
    call mpi_finalize(iat)
    stop
  end if

  if (verb) then
     call yaml_comment('Occupation Numbers',hfill='-')
     call yaml_map('Total Number of Electrons',nelec,fmt='(i8)')
     !write(*,'(1x,a,t28,i8)') 'Total Number of Electrons',nelec
  end if

  ! Number of orbitals
  if (in%nspin==1) then
     if (verb) call yaml_map('Spin treatment','Averaged')
     norb=(nelec+1)/2
     norbu=norb
     norbd=0
     if (mod(nelec,2).ne.0 .and. verb) then
        call yaml_warning('Odd number of electrons, no closed shell system')
        !write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
     end if
  else if(in%nspin==4) then
     if (verb) call yaml_map('Spin treatment','Spinorial (non-collinearity possible)')
     !if (verb) write(*,'(1x,a)') 'Spin-polarized non-collinear calculation'
     norb=nelec
     norbu=norb
     norbd=0
  else 
     if (verb) call yaml_map('Spin treatment','Collinear')
     !if (verb) write(*,'(1x,a)') 'Spin-polarized calculation'
     norb=nelec
     if (mod(norb+in%mpol,2) /=0) then
        write(*,*)'ERROR: the mpol polarization should have the same parity of the number of electrons'
        stop
     end if
     norbu=min((norb+in%mpol)/2,norb)
     norbd=norb-norbu

     !test if the spin is compatible with the input guess polarisations
     ispinsum=0
     ichgsum=0
     do iat=1,atoms%nat
        call charge_and_spol(atoms%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+ispol
        ichgsum=ichgsum+ichg
     end do

     if (in%nspin == 2 .and. ispinsum /= norbu-norbd) then
        !if (iproc==0) then 
           call yaml_warning('Total input polarisation (found ' // trim(yaml_toa(ispinsum)) &
                & // ') must be equal to norbu-norbd.')
           call yaml_comment('With norb=' // trim(yaml_toa(norb)) // ' and mpol=' // trim(yaml_toa(in%mpol)) // &
                & ' norbu-norbd=' // trim((yaml_toa(norbu-norbd))))
           !write(*,'(1x,a,i0,a)')&
           !     'ERROR: Total input polarisation (found ',ispinsum,&
           !     ') must be equal to norbu-norbd.'
           !write(*,'(1x,3(a,i0))')&
           !     'With norb=',norb,' and mpol=',in%mpol,' norbu-norbd=',norbu-norbd
        !end if
        stop
     end if

     if (ichgsum /= in%ncharge .and. ichgsum /= 0) then
        !if (iproc==0) then 
           write(*,'(1x,a,i0,a)')&
                'ERROR: Total input charge (found ',ichgsum,&
                ') cannot be different than charge.'
           write(*,'(1x,2(a,i0))')&
                'The charge is=',in%ncharge,' input charge=',ichgsum
        !end if
        stop
     end if


     !now warn if there is no input guess spin polarisation
     ispinsum=0
     do iat=1,atoms%nat
        call charge_and_spol(atoms%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+abs(ispol)
     end do
     if (ispinsum == 0 .and. in%nspin==2) then
        if (iproc==0 .and. in%norbsempty == 0) &
             call yaml_warning('Found no input polarisation, add it for a correct input guess')
        !write(*,'(1x,a)')&
        !     'WARNING: Found no input polarisation, add it for a correct input guess'
        !stop
     end if

  end if


  !initialise the values for the empty orbitals
  norbuempty=0
  norbdempty=0

  ! Test if the file 'input.occ exists
  !this access is performed at each call_bigdft run
  inquire(file=trim(in%file_occnum),exist=exists)
  iunit=0
  if (exists) then
     iunit=25
     open(unit=iunit,file=trim(in%file_occnum),form='formatted',action='read',status='old')
     if (in%nspin==1) then
        !The first line gives the number of orbitals
        read(unit=iunit,fmt=*,iostat=ierror) nt
     else
        !The first line gives the number of orbitals
        read(unit=iunit,fmt=*,iostat=ierror) ntu,ntd
     end if
     if (ierror /=0) then
        !if (iproc==0) 
          write(*,'(1x,a)') &
             'ERROR: reading the number of orbitals in the file "'//trim(in%file_occnum)//'"'
        stop
     end if
     !Check
     if (in%nspin==1) then
        if (nt<norb) then
           !if (iproc==0) 
               write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "'//trim(in%file_occnum)//'" the number of orbitals norb=',nt,&
                ' should be greater or equal than (nelec+1)/2=',norb
           stop
        else
           norb=nt
           norbu=norb
           norbd=0
        end if
     else
        nt=ntu+ntd
        if (nt<norb) then
           !if (iproc==0) 
               write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "'//trim(in%file_occnum)//'" the number of orbitals norb=',nt,&
                ' should be greater or equal than nelec=',norb
           stop
        else
           norb=nt
        end if
        if (ntu<norbu) then
           !if (iproc==0) 
                write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "'//trim(in%file_occnum)//'" the number of orbitals up norbu=',ntu,&
                ' should be greater or equal than min((nelec+mpol)/2,nelec)=',norbu
           stop
        else
           norbu=ntu
        end if
        if (ntd<norbd) then
           !if (iproc==0) 
                  write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "'//trim(in%file_occnum)//'" the number of orbitals down norbd=',ntd,&
                ' should be greater or equal than min((nelec-mpol/2),0)=',norbd
           stop
        else
           norbd=ntd
        end if
     end if
  else if (in%norbsempty > 0) then
     !total number of orbitals
     norbe=0
     if(in%nspin==4) then
        nspin=2
        nspinor=4
     else
        nspin=in%nspin
        nspinor=1
     end if

     do iat=1,atoms%nat
        ityp=atoms%iatype(iat)
        call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,atoms%aocc(1,iat),occup,nl)
        norbat=(nl(1)+3*nl(2)+5*nl(3)+7*nl(4))
        norbe=norbe+norbat
     end do

     !value of empty orbitals up and down, needed to fill occupation numbers
     norbuempty=min(in%norbsempty,norbe-norbu)
     norbdempty=min(in%norbsempty,norbe-norbd)

     if (in%nspin == 4 .or. in%nspin==1) then
        norb=norb+norbuempty
        norbu=norbu+norbuempty
     else if (in%nspin ==2) then
        norbu=norbu+norbuempty
        norbd=norbd+norbdempty
        norb=norbu+norbd
     end if
  end if

  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if

  call orbitals_descriptors(iproc, nproc,norb,norbu,norbd,in%nspin,nspinor,&
       in%gen_nkpt,in%gen_kpt,in%gen_wkpt,orbs,.false.)

  !distribution of wavefunction arrays between processors
  !tuned for the moment only on the cubic distribution
  if (verb .and. nproc > 1) then
     call yaml_open_map('Orbitals Repartition')
     jpst=0
     do jproc=0,nproc-1
        norbme=orbs%norb_par(jproc,0)
        norbyou=orbs%norb_par(min(jproc+1,nproc-1),0)
        if (norbme /= norbyou .or. jproc == nproc-1) then
           call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//trim(yaml_toa(jproc,fmt='(i0)')),norbme,fmt='(i0)')
           !write(*,'(3(a,i0),a)')&
           !     ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
           jpst=jproc+1
        end if
     end do
     !write(*,'(3(a,i0),a)')&
     !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
     call yaml_close_map()
  end if
  
  if (iproc == 0) then
     if (trim(in%run_name) == '') then
        radical = 'input'
     else
        radical = in%run_name
     end if
     if (verb) call yaml_map('Total Number of Orbitals',norb,fmt='(i8)')
     if (verb) then
        if (iunit /= 0) then
           call yaml_map('Occupation numbers coming from', trim(radical) // '.occ')
        else
           call yaml_map('Occupation numbers coming from','System properties')
        end if
     end if
  end if
  !assign to each k-point the same occupation number
  if (verb .and. iproc==0) call yaml_open_sequence('Input Occupation Numbers')
  do ikpts=1,orbs%nkpts
     if (verb .and. iproc == 0 .and. atoms%geocode /= 'F') then
        call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpts,fmt='(i4.4)'))) // ' BZ coord. = ' // &
        & trim(yaml_toa(orbs%kpts(:, ikpts),fmt='(f12.6)')))
     end if
     call occupation_input_variables(verb,iunit,nelec,norb,norbu,norbuempty,norbdempty,in%nspin,&
          orbs%occup(1+(ikpts-1)*orbs%norb),orbs%spinsgn(1+(ikpts-1)*orbs%norb))
  end do
  if (verb .and. iproc == 0) call yaml_close_sequence()
end subroutine read_orbital_variables

subroutine read_atomic_variables(atoms, fileocc, nspin)
  use module_base
  use module_types
  use module_xc
  use m_ab6_symmetry
  implicit none
  character (len=*), intent(in) :: fileocc
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: nspin
  !local variables
  character(len=*), parameter :: subname='read_atomic_variables'
  integer, parameter :: nelecmax=32,nmax=6,lmax=4,noccmax=2
  character(len=2) :: symbol
  integer :: ityp,iat,ierror,mxpl
  integer :: mxchg,nsccode
  real(gp) :: rcov,rprb,ehomo
  real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
  
  do ityp=1,atoms%ntypes
     ! We calculate atoms%aocc and atoms%amu here.
     call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
          neleconf,nsccode,mxpl,mxchg,atoms%amu(ityp))
     call atomic_occupation_numbers(fileocc,ityp,nspin,atoms,nmax,lmax,nelecmax,&
          neleconf,nsccode,mxpl,mxchg)

     !define the localization radius for the Linear input guess
     atoms%rloc(ityp,:) = rcov * 10.0
  end do
  !print *,'iatsctype',atOMS%iasctype(:)
  atoms%natsc = 0
  do iat=1,atoms%nat
     if (atoms%iasctype(iat) /= 0) atoms%natsc=atoms%natsc+1
  enddo

  ! We modify the symmetry object with respect to the spin.
  if (atoms%sym%symObj >= 0) then
     if (nspin == 2) then
        call symmetry_set_collinear_spin(atoms%sym%symObj, atoms%nat, &
             & atoms%natpol, ierror)
!!$     else if (in%nspin == 4) then
!!$        call symmetry_set_spin(atoms%sym%symObj, atoms%nat, &
!!$             & atoms%natpol, ierror)
     end if
  end if
END SUBROUTINE read_atomic_variables


!> Assign some of the physical system variables
!! Performs also some cross-checks with other variables
!! The pointer in atoms structure have to be associated or nullified.
subroutine print_atomic_variables(atoms, radii_cf, hmax, ixc)
  use module_base
  use module_types
  use module_xc
  use yaml_output
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), intent(in) :: hmax
  integer, intent(in) :: ixc
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  !local variables
  character(len=*), parameter :: subname='print_atomic_variables'
  logical :: nonloc
  integer, parameter :: nelecmax=32,nmax=6,lmax=4,noccmax=2
  integer :: i,j,l,ityp,iat,natyp,mproj
  real(gp) :: minrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,3) :: offdiagarr
  character(len=500) :: name_xc1, name_xc2

!!$  write(*,'(1x,a)')&
!!$       ' Atom    N.Electr.  PSP Code  Radii: Coarse     Fine  CoarsePSP    Calculated   File'

!!$  do ityp=1,atoms%ntypes
!!$     !control the hardest gaussian
!!$     minrad=1.e10_gp
!!$     do i=0,4
!!$        if (atoms%psppar(i,0,ityp)/=0._gp) then
!!$           minrad=min(minrad,atoms%psppar(i,0,ityp))
!!$        end if
!!$     end do
!!$     !control whether the grid spacing is too high
!!$     if (hmax > 2.5_gp*minrad) then
!!$        write(*,'(1x,a)')&
!!$             'WARNING: The grid spacing value may be too high to treat correctly the above pseudo.' 
!!$        write(*,'(1x,a,f5.2,a)')&
!!$             '         Results can be meaningless if hgrid is bigger than',2.5_gp*minrad,&
!!$             '. At your own risk!'
!!$     end if
!!$
!!$     if (atoms%radii_cf(ityp, 1) == UNINITIALIZED(1.0_gp)) then
!!$        message='         X              '
!!$     else
!!$        message='                   X ' 
!!$     end if
!!$     write(*,'(1x,a6,8x,i3,5x,i3,10x,3(1x,f8.5),a)')&
!!$          trim(atoms%atomnames(ityp)),atoms%nelpsp(ityp),atoms%npspcode(ityp),&
!!$          radii_cf(ityp,1),radii_cf(ityp,2),radii_cf(ityp,3),message
!!$  end do
  !print *,'iatsctype',atOMS%iasctype(:)

  !If no atoms...
  if (atoms%ntypes == 0) return

  !print the pseudopotential matrices
  do l=1,3
     do i=1,2
        do j=i+1,3
           offdiagarr(i,j-i,l)=0._gp
           if (l==1) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagarr(i,j-i,l)=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagarr(i,j-i,l)=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagarr(i,j-i,l)=-9._gp*sqrt(1._gp/143._gp)
              end if
           end if
        end do
     end do
  end do

!  write(*,'(1x,a)')&
  !       '------------------------------------ Pseudopotential coefficients (Upper Triangular)'
  call yaml_comment('System Properties',hfill='-')
  call yaml_open_sequence('Properties of atoms in the system')
  do ityp=1,atoms%ntypes
     call yaml_sequence(advance='no')
     call yaml_map('Symbol',trim(atoms%atomnames(ityp)),advance='no')
     call yaml_comment('Type No. '//trim(yaml_toa(ityp,fmt='(i2.2)')))
     call yaml_map('No. of Electrons',atoms%nelpsp(ityp))
     natyp=0
     do iat=1,atoms%nat
        if (atoms%iatype(iat) == ityp) natyp=natyp+1
     end do
     call yaml_map('No. of Atoms',natyp)

     call yaml_open_map('Radii of active regions (AU)')!,flow=.true.)
       call yaml_map('Coarse',radii_cf(ityp,1),fmt='(f8.5)')
       call yaml_map('Fine',radii_cf(ityp,2),fmt='(f8.5)')
       call yaml_map('Coarse PSP',radii_cf(ityp,3),fmt='(f8.5)')
       if (atoms%radii_cf(ityp, 1) == UNINITIALIZED(1.0_gp)) then
          call yaml_map('Source','Hard-Coded')
       else
          call yaml_map('Source','PSP File')
       end if
     call yaml_close_map()

     minrad=1.e10_gp
     do i=0,4
        if (atoms%psppar(i,0,ityp)/=0._gp) then
           minrad=min(minrad,atoms%psppar(i,0,ityp))
        end if
     end do
     if (radii_cf(ityp,2) /=0.0_gp) then
        call yaml_map('Grid Spacing threshold (AU)',2.5_gp*minrad,fmt='(f5.2)')
     else
        call yaml_map('Grid Spacing threshold (AU)',1.25_gp*minrad,fmt='(f5.2)')
     end if
     !control whether the grid spacing is too high
     if (hmax > 2.5_gp*minrad) then
        call yaml_warning('Chosen Grids spacings seem too high for this atom. At you own risk!')
!!$        write(*,'(1x,a)')&
!!$             'WARNING: The grid spacing value may be too high to treat correctly the above pseudo.' 
!!$        write(*,'(1x,a,f5.2,a)')&
!!$             '         Results can be meaningless if hgrid is bigger than',2.5_gp*minrad,&
!!$             '. At your own risk!'
     end if

     select case(atoms%npspcode(ityp))
     case(2)
        call yaml_map('Pseudopotential type','GTH')
     case(3)
        call yaml_map('Pseudopotential type','HGH')
     case(10)
        call yaml_map('Pseudopotential type','HGH-K')
     end select
     if (atoms%psppar(0,0,ityp)/=0) then
        call yaml_open_map('Local Pseudo Potential (HGH convention)')
          call yaml_map('Rloc',atoms%psppar(0,0,ityp),fmt='(f9.5)')
          call yaml_map('Coefficients (c1 .. c4)',atoms%psppar(0,1:4,ityp),fmt='(f9.5)')
        call yaml_close_map()
     end if
     !see if nonlocal terms are present
     nonloc=.false.
     verify_nl: do l=1,3
        do i=3,0,-1
           j=i
           if (atoms%psppar(l,i,ityp) /= 0._gp) exit
        end do
        if (j /=0) then
           nonloc=.true.
           exit verify_nl
        end if
     end do verify_nl
     if (nonloc) then
        call yaml_open_sequence('NonLocal PSP Parameters')
        do l=1,3
           do i=3,0,-1
              j=i
              if (atoms%psppar(l,i,ityp) /= 0._gp) exit
           end do
           if (j /=0) then
              call yaml_sequence(advance='no')
              call yaml_map('Channel (l)',l-1)
              call yaml_map('Rloc',atoms%psppar(l,0,ityp),fmt='(f9.5)')
              hij=0._gp
              do i=1,j
                 hij(i,i)=atoms%psppar(l,i,ityp)
              end do
              if (atoms%npspcode(ityp) == 3) then !traditional HGH convention
                 hij(1,2)=offdiagarr(1,1,l)*atoms%psppar(l,2,ityp)
                 hij(1,3)=offdiagarr(1,2,l)*atoms%psppar(l,3,ityp)
                 hij(2,3)=offdiagarr(2,1,l)*atoms%psppar(l,3,ityp)
              else if (atoms%npspcode(ityp) == 10) then !HGH-K convention
                 hij(1,2)=atoms%psppar(l,4,ityp)
                 hij(1,3)=atoms%psppar(l,5,ityp)
                 hij(2,3)=atoms%psppar(l,6,ityp)
              end if
              call yaml_open_sequence('h_ij matrix')
                call yaml_sequence(trim(yaml_toa(hij(1,1:3),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,2),hij(2,2),hij(2,3)/),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,3),hij(2,3),hij(3,3)/),fmt='(f9.5)')))
              call yaml_close_sequence()
           end if
        end do
        call yaml_close_sequence()
     end if
     call numb_proj(ityp,atoms%ntypes,atoms%psppar,atoms%npspcode,mproj)
     call yaml_map('No. of projectors',mproj)

!!$     write(*,'(1x,a)')&
!!$          'Atom Name    rloc      C1        C2        C3        C4  '
!!$     do l=0,4
!!$        if (l==0) then
!!$           do i=4,0,-1
!!$              j=i
!!$              if (atoms%psppar(l,i,ityp) /= 0._gp) exit
!!$           end do
!!$           write(*,'(3x,a6,5(1x,f9.5))')&
!!$                trim(atoms%atomnames(ityp)),(atoms%psppar(l,i,ityp),i=0,j)
!!$        else
!!$           do i=3,0,-1
!!$              j=i
!!$              if (atoms%psppar(l,i,ityp) /= 0._gp) exit
!!$           end do
!!$           if (j /=0) then
!!$              write(*,'(1x,a,i0,a)')&
!!$                   '    l=',l-1,' '//'     rl        h1j       h2j       h3j '
!!$              hij=0._gp
!!$              do i=1,j
!!$                 hij(i,i)=atoms%psppar(l,i,ityp)
!!$              end do
!!$              if (atoms%npspcode(ityp) == 3) then !traditional HGH convention
!!$                 hij(1,2)=offdiagarr(1,1,l)*atoms%psppar(l,2,ityp)
!!$                 hij(1,3)=offdiagarr(1,2,l)*atoms%psppar(l,3,ityp)
!!$                 hij(2,3)=offdiagarr(2,1,l)*atoms%psppar(l,3,ityp)
!!$              else if (atoms%npspcode(ityp) == 10) then !HGH-K convention
!!$                 hij(1,2)=atoms%psppar(l,4,ityp)
!!$                 hij(1,3)=atoms%psppar(l,5,ityp)
!!$                 hij(2,3)=atoms%psppar(l,6,ityp)
!!$              end if
!!$              do i=1,j
!!$                 if (i==1) then
!!$                    write(format,'(a,2(i0,a))')"(9x,(1x,f9.5),",j,"(1x,f9.5))"
!!$                    write(*,format)atoms%psppar(l,0,ityp),(hij(i,k),k=i,j)
!!$                 else
!!$                    write(format,'(a,2(i0,a))')"(19x,",i-1,"(10x),",j-i+1,"(1x,f9.5))"
!!$                    write(*,format)(hij(i,k),k=i,j)
!!$                 end if
!!$
!!$              end do
!!$           end if
!!$        end if
!!$     end do
!!$     !control if the PSP is calculated with the same XC value
     if (atoms%ixcpsp(ityp) < 0) then
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_MIXED)
     else
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_ABINIT)
     end if
     if (ixc < 0) then
        call xc_get_name(name_xc2, ixc, XC_MIXED)
     else
        call xc_get_name(name_xc2, ixc, XC_ABINIT)
     end if
     call yaml_map('PSP XC','"'//trim(name_xc1)//'"')
     if (trim(name_xc1) /= trim(name_xc2)) then
        call yaml_warning('Input ixc parameter corresponds to '//trim(name_xc2)//' XC functional')
!!$        write(*,'(1x,a)')&
!!$             'WARNING: The pseudopotential file psppar."'//trim(atoms%atomnames(ityp))//'"'
!!$        write(*,'(1x,a,i0,a,i0)')&
!!$             '         contains a PSP generated with an XC id=',&
!!$             atoms%ixcpsp(ityp),' while for this run ixc=',ixc
     end if
  end do
  call yaml_close_sequence()
!!!  tt=dble(norb)/dble(nproc)
!!!  norbp=int((1.d0-eps_mach*tt) + tt)
!!!  !if (verb.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp


  ! if linear scaling applied with more then InputGuess, then go read input.lin for radii
  !  if (in%linear /= 'OFF' .and. in%linear /= 'LIG') then
  !     lin%nlr=atoms%nat
  !     call allocateBasicArrays(atoms, lin)
  !     call readLinearParameters(verb, nproc, lin, atoms, atomNames)
  !  end if
END SUBROUTINE print_atomic_variables


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


!> Fix all the atomic occupation numbers of the atoms which has the same type
!! look also at the input polarisation and spin
!! look at the file of the input occupation numbers and, if exists, modify the 
!! occupations accordingly
subroutine atomic_occupation_numbers(filename,ityp,nspin,at,nmax,lmax,nelecmax,neleconf,nsccode,mxpl,mxchg)
  use module_base
  use module_types
  use yaml_output
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: ityp,mxpl,mxchg,nspin,nmax,lmax,nelecmax,nsccode
  type(atoms_data), intent(inout) :: at
  !integer, dimension(nmax,lmax), intent(in) :: neleconf
  real(gp), dimension(nmax,lmax), intent(in) :: neleconf
  !local variables
  integer, parameter :: noccmax=2
  character(len=100) :: string
  logical :: exists,found
  integer :: iat,ichg,ispol,nsp,nspinor,ierror,jat,l,iocc,icoll,noncoll,ispin,nl,inl,m
  real(gp) :: elec
  real(gp), dimension(nmax,lmax) :: eleconf

  !control the spin
  select case(nspin)
     case(1)
        nsp=1
        nspinor=1
        noncoll=1
     case(2)
        nsp=2
        nspinor=1
        noncoll=1
     case(4)
        nsp=1
        nspinor=4
        noncoll=2
     case default
        call yaml_warning('nspin not valid:' // trim(yaml_toa(nspin)))
        !write(*,*)' ERROR: nspin not valid:',nspin
        stop
  end select

  inquire(file=filename,exist=exists)

  !search the corresponding atom
  if (exists) then
     open(unit=91,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        call yaml_warning('Failed to open the existing  file: '// trim(filename))
        !write(*,*)'Failed to open the existing  file: '//filename
        stop
     end if
  end if

  !here we must check of the input guess polarisation
  !control if the values are compatible with the atom configuration
  !do this for all atoms belonging to a given type
  !control the maximum polarisation allowed: consider only non-closed shells   
  do iat=1,at%nat
     !control the atomic input guess occupation number
     !if you want no semicore input guess electrons, uncomment the following line
     !at%iasctype(iat)=0
     if (at%iatype(iat) == ityp) then
        !see whether the input.occup file contains the given atom
        found=.false.
        if (exists) then
           rewind(unit=91)
           parse_inocc: do
              read(91,'(a100)',iostat=ierror)string
              if (ierror /= 0) exit parse_inocc !file ends
              read(string,*,iostat=ierror)jat
              if (ierror /=0) stop 'Error reading line'
              if (jat==iat ) then
                 found=.true.
                 exit parse_inocc
              end if
           end do parse_inocc
        end if
        call charge_and_spol(at%natpol(iat),ichg,ispol)
        if (found) then
           call read_eleconf(string,nsp,nspinor,noccmax,nelecmax,lmax,&
                at%aocc(1,iat),at%iasctype(iat))
        else
           at%iasctype(iat)=nsccode
           if (abs(ispol) > mxpl+abs(ichg)) then
              !if (iproc ==0) 
              write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input polarisation of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxpl,&
                   ', while found ',ispol
              stop 
           end if
           if (abs(ichg) > mxchg) then
              !if (iproc ==0) 
              write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input charge of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxchg,&
                   ', while found ',ichg
              stop
           end if
           !correct the electronic configuration in case there is a charge
           !if (ichg /=0) then
           call correct_semicore(nmax,lmax-1,ichg,&
                neleconf,eleconf,at%iasctype(iat))
           !end if

           call at_occnums(ispol,nsp,nspinor,nmax,lmax,nelecmax,&
                eleconf,at%aocc(1,iat))
        end if

        !check the total number of electrons
        elec=0.0_gp
        iocc=0
        do l=1,lmax
           iocc=iocc+1
           nl=nint(at%aocc(iocc,iat))
           do inl=1,nl
              do ispin=1,nsp
                 do m=1,2*l-1
                    do icoll=1,noncoll !non-trivial only for nspinor=4
                       iocc=iocc+1
                       elec=elec+at%aocc(iocc,iat)
                    end do
                 end do
              end do
           end do
        end do
        if (nint(elec) /= at%nelpsp(ityp) - ichg) then
           write(*,*)'ERROR: the total atomic charge ',elec,&
                ' is different from the PSP charge ',at%nelpsp(ityp),&
                ' plus the charge ',-ichg
           stop
        end if
     end if
  end do

  if (exists) close(unit=91)

END SUBROUTINE atomic_occupation_numbers


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


!> Routine which assign to each processor the repartition of nobj*nkpts objects
subroutine kpts_to_procs_via_obj(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,nobj
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_par
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
     rounding_ratio=(robjp-real(floor(robjp)))
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
        !print *,'qui',nproc_left,nproc_per_kpt
        do kproc=0,nproc_left-1
           ikpt=ikpt+1
           if (ikpt > nkpts) stop 'ERROR: also this should not happen3'
           do iobj=0,nobj-1
              nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)+1
           end do
           jproc=jproc+nproc_per_kpt+1
        end do
        !print *,'ciao'
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
        !print *,'qui',nkpts_left,nkpts_per_proc
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
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,nvctr,norb
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nvctr_par
  !local variables
  integer :: ikpt,jsproc,jeproc,kproc,icount,ivctr,jproc,numproc
  real(gp) :: strprc,endprc

  ! This variable qas not initialized...
  icount=0

  !for any of the k-points find the processors which have such k-point associated
  call to_zero(nproc*nkpts,nvctr_par(0,1))

  do ikpt=1,nkpts
     jsproc=UNINITIALIZED(1)
     jeproc=UNINITIALIZED(1)
     find_start: do jproc=0,nproc-1
        if(norb_par(jproc,ikpt) > 0) then 
           jsproc=jproc
           exit find_start
        end if
     end do find_start
     if (jsproc == UNINITIALIZED(1)) stop 'ERROR in kpt assignments'
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
        if (jeproc == UNINITIALIZED(1)) stop 'ERROR in kpt assignments'
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
     !evaluate the percentace on the number of components
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
           nvctr_par(kproc,ikpt)=&
                nvctr_par(kproc,ikpt)+1
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
  integer :: ikpt,jproc,norbs,ncomps,i_all,i_stat,kproc,ieproc,isproc,jkpt
  integer, dimension(:,:), allocatable :: load_unbalancing
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  if (info == 0) call print_distribution_schemes(6,nproc,nkpts,norb_par,ncomp_par)

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

  allocate(load_unbalancing(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,load_unbalancing,'load_unbalancing',subname)

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
  i_all=-product(shape(load_unbalancing))*kind(load_unbalancing)
  deallocate(load_unbalancing,stat=i_stat)
  call memocc(i_stat,i_all,'load_unbalancing',subname)


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
  integer :: npawl, ipawl, paw_l, i_stat
  integer :: paw_nofgaussians, paw_nofchannels, il, ierror, ig
  real(gp) :: paw_greal, paw_gimag, paw_ccoeff, paw_scoeff, dumpaw
  character(len=100) :: string

  !parameters for abscalc-paw

  if(.not. storeit) then
     !if(ityp == 1) then !this implies that the PSP are all present
     if (.not. associated(atoms%paw_NofL)) then
        allocate(atoms%paw_NofL(atoms%ntypes+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_NofL,'atoms%paw_NofL',subname)
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
        allocate(atoms%paw_l  (paw_tot_l+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_l,'atoms%paw_l',subname)
        
        allocate(atoms%paw_nofchannels  (paw_tot_l+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_nofchannels,'atoms%paw_nofchannels',subname)
        
        allocate(atoms%paw_nofgaussians  (paw_tot_l+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_nofgaussians,'atoms%paw_nofgaussians',subname)
        
        allocate(atoms%paw_Greal  (paw_tot_l+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_Greal,'atoms%paw_Greal',subname)
        
        allocate(atoms%paw_Gimag ( paw_tot_q   +  ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_Gimag,'atoms%paw_Gimag',subname)
        
        allocate(atoms%paw_Gcoeffs ( paw_tot_coefficients  +  ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_Gcoeffs,'atoms%paw_Gcoeffs',subname)
        
        allocate(atoms%paw_H_matrices(paw_tot_matrices+ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_H_matrices,'atoms%paw_H_matrices',subname)
        
        allocate(atoms%paw_S_matrices ( paw_tot_matrices  +  ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_S_matrices,'atoms%paw_S_matrices',subname)
        
        
        allocate(atoms%paw_Sm1_matrices ( paw_tot_matrices  +  ndebug), stat=i_stat)
        call memocc(i_stat,atoms%paw_Sm1_matrices,'atoms%paw_Sm1_matrices',subname)
        
        
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

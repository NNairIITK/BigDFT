!> @file
!! Toy program to use BigDFT API
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Toy program to use BigDFT API
program wvl

  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use BigDFT_API
  use locreg_operations, only: workarr_sumrho,initialize_work_arrays_sumrho,deallocate_work_arrays_sumrho
  use bigdft_run
  use dynamic_memory
  use yaml_output
  use module_input_dicts
  use module_input_keys
  use module_atoms, only: deallocate_atoms_data
  use communications_base, only: deallocate_comms
  use communications_init, only: orbitals_communicators
  use communications, only: transpose_v, untranspose_v
  use rhopotential, only: full_local_potential
  use module_razero
  implicit none

  type(input_variables)             :: inputs
  type(atoms_data)                  :: atoms

  type(local_zone_descriptors)          :: Lzd
  type(orbitals_data)               :: orbs
  type(comms_cubic)       :: comms
  type(workarr_sumrho)              :: wisf
  real(wp), dimension(:), pointer   :: psi, psir

  type(xc_info)                        :: xc
  type(rho_descriptors)                :: rhodsc
  type(denspot_distribution)           :: dpcom
  type(GPU_pointers)                   :: GPU
  integer :: i, j, ierr, iproc, nproc ,nconfig
  real(dp) :: nrm, epot_sum
  real(gp) :: psoffset
!  real(gp), allocatable :: radii_cf(:,:)
  real(gp), dimension(3) :: shift
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(dp), dimension(:), pointer   :: rhor, pot_ion, potential,rho_ion
  real(wp), dimension(:), pointer   :: w
  real(wp), dimension(:,:), pointer :: ovrlp
  real(dp), dimension(:,:), pointer :: rho_p => null() !needs to be nullified
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons
  type(coulomb_operator) :: pkernel
  type(dictionary), pointer :: user_inputs,options
  !temporary variables
  !integer, dimension(4) :: mpi_info
  character(len=60) :: run_id

  call f_lib_initialize()
  nullify(options)
  !-initializes the mpi_environment for each group
  call bigdft_init(options)
  call dict_free(options)
!just for backward compatibility
  iproc=bigdft_mpi%iproc
  nproc=bigdft_mpi%nproc
   call dict_init(user_inputs)
   call user_dict_from_files(user_inputs, 'input', 'posinp', bigdft_mpi)
   call inputs_from_dict(inputs, atoms, user_inputs)
   if (iproc == 0) then
      call print_general_parameters(inputs,atoms,'input','posinp')
   end if
   call dict_free(user_inputs)
   GPU%OCLconv = .false.

!!$  ! Start MPI in parallel version
!!$  call MPI_INIT(ierr)
!!$  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
!!$  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!!$
!!$  call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,0)
!!$
!!$  if (iproc==0) call print_logo()
!!$
!!$  ! Setup names for input and output files.
!!$  call standard_inputfile_names(inputs, "toy",nproc)
!!$  ! Read all input stuff, variables and atomic coordinates and pseudo.
!!$  !the arguments of this routine should be changed
!!$  posinp_name='posinp'
!!$  call read_input_variables(iproc,nproc,posinp_name,inputs, atoms, atoms%astruct%rxyz,1,'input',0)

!  allocate(radii_cf(atoms%astruct%ntypes,3))
  call system_properties(iproc,nproc,inputs,atoms,orbs)!,radii_cf)
  
  call lzd_set_hgrids(Lzd,(/inputs%hx,inputs%hy,inputs%hz/)) 
  call system_size(atoms,atoms%astruct%rxyz,inputs%crmult,inputs%frmult, &
       & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),GPU%OCLconv,Lzd%Glr,shift)
  call print_atoms_and_grid(Lzd%Glr, atoms, atoms%astruct%rxyz, shift, &
       & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3))

  ! Setting up the wavefunction representations (descriptors for the
  !  compressed form...).
  call createWavefunctionsDescriptors(iproc,inputs%hx,inputs%hy,inputs%hz, &
       & atoms,atoms%astruct%rxyz,inputs%crmult,inputs%frmult,.true.,Lzd%Glr)
  call print_wfd(Lzd%Glr%wfd)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)  

  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbs,inputs%nspin,atoms%astruct%rxyz)

  !grid spacings and box of the density
  call dpbox_set(dpcom,Lzd,xc,iproc,nproc,MPI_COMM_WORLD,inputs%PSolver_groupsize, &
       & inputs%SIC%approach,atoms%astruct%geocode, inputs%nspin,inputs%matacc%PSolver_igpu)

  ! Read wavefunctions from disk and store them in psi.
  allocate(orbs%eval(orbs%norb*orbs%nkpts))
  call f_zero(orbs%eval)
  allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)))
  allocate(rxyz_old(3, atoms%astruct%nat))
  call readmywaves(iproc,"data/wavefunction",WF_FORMAT_PLAIN,orbs,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3, &
       & inputs%hx,inputs%hy,inputs%hz,atoms,rxyz_old,atoms%astruct%rxyz,Lzd%Glr%wfd,psi)
  if (nproc>1) then
      call mpiallred(orbs%eval(1),orbs%norb*orbs%nkpts,MPI_SUM)
  end if


  ! Some analysis.
  call yaml_comment("Proc" // trim(yaml_toa(iproc)) // " allocates psi to " // &
                   & trim(yaml_toa(max(orbs%npsidim_orbs,orbs%npsidim_comp))))
  !write(*,*) "Proc", iproc, " allocates psi to",max(orbs%npsidim_orbs,orbs%npsidim_comp)

  call yaml_flush_document()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !-------------------------!
  ! The orbital repartition !
  !-------------------------!
  do i = 1, orbs%norbp, 1
     !write(*,*) "Proc", iproc, " treats orbital", orbs%isorb + i
     call yaml_comment("Proc" // trim(yaml_toa(iproc)) // " treats orbital" // &
                      & trim(yaml_toa(orbs%isorb + i)))
  end do
  
  ! We can do some arithmetic on wavefunctions under compressed form.
  do i = 1, orbs%norbp, 1
     ! Norm calculation.
     call wnrm_wrap(1, Lzd%Glr%wfd%nvctr_c, Lzd%Glr%wfd%nvctr_f, &
          & psi((i - 1) * (Lzd%Glr%wfd%nvctr_c + 7 * Lzd%Glr%wfd%nvctr_f) + 1), nrm)
     !write(*,*) "Proc", iproc, " orbital", orbs%isorb + i, " is of norm ", nrm
     call yaml_comment("Proc" // trim(yaml_toa(iproc)) // " orbital" // &
                      & trim(yaml_toa(orbs%isorb + i)) // " is of norm " // trim(yaml_toa(nrm)))
  end do

  call yaml_flush_document()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !---------------------------!
  ! The component repartition !
  !---------------------------!
  allocate(w(max(orbs%npsidim_orbs,orbs%npsidim_comp)))
  ! Transpose the psi wavefunction
  call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),work_add=w(1))
  !write(*,*) "Proc", iproc, " treats ", comms%nvctr_par(iproc, 0) * orbs%norb, "components of all orbitals."
  call yaml_comment("Proc" // trim(yaml_toa(iproc)) // " treats " // &
                   & trim(yaml_toa(comms%nvctr_par(iproc, 0) * orbs%norb)) // "components of all orbitals.")

  ! Calculate the overlap matrix, thanks to orthogonality of
  ! Daubechies wavelets, one can directly multiply the coefficients.
  !  See getOverlap() function for a complete description.
  allocate(ovrlp(orbs%norb, orbs%norb))
  do j = 1, orbs%norb, 1
     do i = 1, orbs%norb, 1
        ovrlp(i, j) = dot_double(comms%nvctr_par(iproc, 0), &
             & psi((i - 1) * comms%nvctr_par(iproc, 0) + 1), 1, &
             & psi((j - 1) * comms%nvctr_par(iproc, 0) + 1), 1)
     end do
  end do
  ! This double loop can be expressed with BLAS DSYRK function.
  !  call syrk('L','T',orbs%norb,comms%nvctr_par(iproc, 0),1.0_wp,psi(1), &
  !       & max(1,comms%nvctr_par(iproc, 0)),0.0_wp,ovrlp(1,1),orbs%norb)
  if (nproc>1) then
      call mpiallred(ovrlp(1,1),orbs%norb * orbs%norb,MPI_SUM)
  end if
  if (iproc == 0) then
     !uses yaml_output routine to provide example
     call yaml_sequence_open('The overlap matrix is')
          do j = 1, orbs%norb, 1
             call yaml_sequence(trim(yaml_toa(ovrlp(:, j),fmt='(g18.8)')))
          end do
     call yaml_sequence_close()
     !write(*,*) "The overlap matrix is:"
     !do j = 1, orbs%norb, 1
     !   write(*, "(A)", advance = "NO") "("
     !   do i = 1, orbs%norb, 1  
     !      write(*,"(G18.8)", advance = "NO") ovrlp(i, j)
     !   end do
     !   write(*, "(A)") ")"
     !end do
  end if
  deallocate(ovrlp)

  ! Retranspose the psi wavefunction
  call untranspose_v(iproc,nproc,orbs,Lzd%glr%wfd,comms,psi(1),work_add=w(1))
  deallocate(w)

  call yaml_flush_document()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !-------------------------!
  ! Using real-space values !
  !-------------------------!
  ! Wavefunctions can be expressed in interpolating scaling functions, 
  !  in this representation, point coefficients are values on points.
  allocate(psir(Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * Lzd%Glr%d%n3i))
  call razero(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,psir(1))
  ! BigDFT cut rho by slices, while here we keep one array for simplicity.
  allocate(rhor(Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * Lzd%Glr%d%n3i))
  call razero(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,rhor(1))
  call initialize_work_arrays_sumrho(1,[Lzd%Glr],.true.,wisf)
  do i = 1, orbs%norbp, 1
     ! Calculate values of psi_i on each grid points.
     call daub_to_isf(Lzd%Glr,wisf, &
          & psi((i - 1) * (Lzd%Glr%wfd%nvctr_c + 7 * Lzd%Glr%wfd%nvctr_f) + 1),psir)
     ! Compute partial densities coming from occup * psi_i * psi_i.
     do j = 1, Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * Lzd%Glr%d%n3i, 1
        rhor(j) = rhor(j) + orbs%occup(orbs%isorb + i) * psir(j) * psir(j)
     end do
  end do
  if (nproc>1) then
      call mpiallred(rhor(1),Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * Lzd%Glr%d%n3i,MPI_SUM)
  end if
  !if (iproc == 0) write(*,*) "System has", sum(rhor), "electrons."
  if (iproc == 0) call yaml_map("Number of electrons", sum(rhor))
  deallocate(rhor)

  call xc_init(xc, inputs%ixc, XC_ABINIT, inputs%nspin)
  call density_descriptors(iproc,nproc,xc,inputs%nspin,inputs%crmult,inputs%frmult,atoms,&
       dpcom,inputs%rho_commun,atoms%astruct%rxyz,rhodsc)

!!$  ! Equivalent BigDFT routine.
!!$  allocate(nscatterarr(0:nproc-1,4))
!!$  allocate(ngatherarr(0:nproc-1,2))
!!$  call createDensPotDescriptors(iproc,nproc,atoms,Lzd%Glr%d, &
!!$       & inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
!!$       & atoms%astruct%rxyz,inputs%crmult,inputs%frmult,radii_cf,inputs%nspin,'D',inputs%ixc, &
!!$       & inputs%rho_commun,n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
  call local_potential_dimensions(iproc,Lzd,orbs,xc,dpcom%ngatherarr(0,1))

  allocate(rhor(Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * dpcom%n3d))
  allocate(irrzon(1,2,1))
  allocate(phnons(2,1,1))

  !call sumrho(iproc,nproc,orbs,Lzd%Glr,inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
  !     & psi,rhor,nscatterarr,inputs%nspin,GPU,atoms%symObj,irrzon,phnons,rhodsc)
  call sumrho(dpcom,orbs,Lzd,GPU,atoms%astruct%sym,rhodsc,xc,psi,rho_p)
  call communicate_density(dpcom,orbs%nspin,rhodsc,rho_p,rhor,.false.)

  call deallocate_rho_descriptors(rhodsc)

  ! Example of calculation of the energy of the local potential of the pseudos.
  pkernel=pkernel_init(.true.,iproc,nproc,0,&
       atoms%astruct%geocode,(/Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i/),&
       (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/),16)
  call pkernel_set(pkernel,verbose=.false.)
  !call createKernel(iproc,nproc,atoms%astruct%geocode,&
  !     (/Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i/), &
  !     (/inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp/)&
  !     ,16,pkernel,.false.)
  allocate(pot_ion(Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * dpcom%n3p))
  allocate(rho_ion(Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * dpcom%n3p))
  call createIonicPotential(atoms%astruct%geocode,iproc,nproc,(iproc==0),atoms,atoms%astruct%rxyz,&
       & inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
       & inputs%elecfield,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3, &
       & dpcom%n3pi,dpcom%i3s+dpcom%i3xcsh,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i, &
       & pkernel,pot_ion,rho_ion,psoffset)
  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,xc,pot_ion,potential)
!!$  call full_local_potential(iproc,nproc,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*n3p, &
!!$       & Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,inputs%nspin, &
!!$       & Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*n3d,0, &
!!$       & orbs,Lzd,0,ngatherarr,pot_ion,potential)
  epot_sum = 0._dp
  do i = 1, orbs%norbp, 1
     call daub_to_isf(Lzd%Glr,wisf, &
          & psi((i - 1) * (Lzd%Glr%wfd%nvctr_c + 7 * Lzd%Glr%wfd%nvctr_f) + 1),psir)
     do j = 1, Lzd%Glr%d%n1i * Lzd%Glr%d%n2i * Lzd%Glr%d%n3i, 1
        epot_sum = epot_sum + psir(j) * potential(j) * psir(j)
     end do
  end do
  epot_sum = epot_sum * inputs%hx / 2._gp * inputs%hy / 2._gp * inputs%hz / 2._gp
  call free_full_potential(dpcom%mpi_env%nproc,0,xc,potential,"main")
  if (nproc>1) then
      call mpiallred(epot_sum,1,MPI_SUM)
  end if
  
  !if (iproc == 0) write(*,*) "System pseudo energy is", epot_sum, "Ht."
  if (iproc == 0) call yaml_map("System pseudo energy (Ha)", epot_sum)

  deallocate(pot_ion,rho_ion)
  deallocate(psir)
  call deallocate_work_arrays_sumrho(wisf)

  ! Free allocated space.
  deallocate(rxyz_old)
  deallocate(psi)

  call deallocate_comms(comms)
  
  call deallocate_locreg_descriptors(Lzd%Glr)

  call deallocate_Lzd_except_Glr(Lzd)
  !deallocate(Lzd%Glr%projflg)

  call deallocate_orbs(orbs)

  call deallocate_atoms_data(atoms) 
  call xc_end(xc)
  call dpbox_free(dpcom)
  call pkernel_free(pkernel)
  call free_input_variables(inputs)

  !wait all processes before finalisation
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
  call f_lib_finalize()
end program wvl

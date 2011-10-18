program wvl

  use Poisson_Solver
  use BigDFT_API
  
  implicit none

  type(input_variables)             :: inputs
  type(atoms_data)                  :: atoms
  real(gp), dimension(:,:), pointer :: rxyz

  type(locreg_descriptors)          :: Glr
  type(orbitals_data)               :: orbs
  type(communications_arrays)       :: comms
  type(workarr_sumrho)              :: wisf
  real(wp), dimension(:), pointer   :: psi, psir

  type(rho_descriptors)                :: rhodsc
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  type(GPU_pointers)                   :: GPU
  
  integer :: i, j, ierr, iproc, nproc, nelec
  integer :: n3d,n3p,n3pi,i3xcsh,i3s
  real(dp) :: nrm, epot_sum
  real(gp) :: psoffset
  real(gp), allocatable :: radii_cf(:,:)
  real(gp), dimension(3) :: shift
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(dp), dimension(:), pointer   :: rhor, pot_ion, potential
  real(wp), dimension(:), pointer   :: w
  real(wp), dimension(:,:), pointer :: ovrlp
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons
  real(dp), dimension(:), pointer :: pkernel

  ! Start MPI in parallel version
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  if (iproc==0) call print_logo()

  ! Setup names for input and output files.
  call standard_inputfile_names(inputs, "toy")
  ! Read all input stuff, variables and atomic coordinates and pseudo.
  call read_input_variables(iproc,"posinp",inputs, atoms, rxyz)


  ! Setting up the size of the calculation (description of the box and
  !  calculation area).
  allocate(radii_cf(atoms%ntypes,3))
  call system_properties(iproc,nproc,inputs,atoms,orbs,radii_cf,nelec)
  call system_size(iproc,atoms,rxyz,radii_cf,inputs%crmult,inputs%frmult, &
       & inputs%hx,inputs%hy,inputs%hz,Glr,shift)

  ! Setting up the wavefunction representations (descriptors for the
  !  compressed form...).
  call createWavefunctionsDescriptors(iproc,inputs%hx,inputs%hy,inputs%hz, &
       & atoms,rxyz,radii_cf,inputs%crmult,inputs%frmult,Glr)
  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  

  ! Read wavefunctions from disk and store them in psi.
  allocate(orbs%eval(orbs%norb*orbs%nkpts))
  call to_zero(orbs%norb*orbs%nkpts,orbs%eval(1))
  allocate(psi(orbs%npsidim))
  allocate(rxyz_old(3, atoms%nat))
  call readmywaves(iproc,"data/wavefunction", orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3, &
       & inputs%hx,inputs%hy,inputs%hz,atoms,rxyz_old,rxyz,Glr%wfd,psi)
  call mpiallred(orbs%eval(1),orbs%norb*orbs%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)



  ! Some analysis.
  write(*,*) "Proc", iproc, " allocates psi to", orbs%npsidim
  call flush(6)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !-------------------------!
  ! The orbital repartition !
  !-------------------------!
  do i = 1, orbs%norbp, 1
     write(*,*) "Proc", iproc, " treats orbital", orbs%isorb + i
  end do
  
  ! We can do some arithmetic on wavefunctions under compressed form.
  do i = 1, orbs%norbp, 1
     ! Norm calculation.
     call wnrm_wrap(1, Glr%wfd%nvctr_c, Glr%wfd%nvctr_f, &
          & psi((i - 1) * (Glr%wfd%nvctr_c + 7 * Glr%wfd%nvctr_f) + 1), nrm)
     write(*,*) "Proc", iproc, " orbital", orbs%isorb + i, " is of norm ", nrm
  end do

  call flush(6)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !---------------------------!
  ! The component repartition !
  !---------------------------!
  allocate(w(orbs%npsidim))
  ! Transpose the psi wavefunction
  call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,psi, work=w)
  write(*,*) "Proc", iproc, " treats ", comms%nvctr_par(iproc, 0) * orbs%norb, "components of all orbitals."

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
  call mpiallred(ovrlp(1,1),orbs%norb * orbs%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (iproc == 0) then
     write(*,*) "The overlap matrix is:"
     do j = 1, orbs%norb, 1
        write(*, "(A)", advance = "NO") "("
        do i = 1, orbs%norb, 1  
           write(*,"(G18.8)", advance = "NO") ovrlp(i, j)
        end do
        write(*, "(A)") ")"
     end do
  end if
  deallocate(ovrlp)

  ! Retranspose the psi wavefunction
  call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi, work=w)
  deallocate(w)

  call flush(6)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !-------------------------!
  ! Using real-space values !
  !-------------------------!
  ! Wavefunctions can be expressed in interpolating scaling functions, 
  !  in this representation, point coefficients are values on points.
  allocate(psir(Glr%d%n1i * Glr%d%n2i * Glr%d%n3i))
  call razero(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,psir)
  ! BigDFT cut rho by slices, while here we keep one array for simplicity.
  allocate(rhor(Glr%d%n1i * Glr%d%n2i * Glr%d%n3i))
  call razero(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,rhor)
  call initialize_work_arrays_sumrho(Glr,wisf)
  do i = 1, orbs%norbp, 1
     ! Calculate values of psi_i on each grid points.
     call daub_to_isf(Glr,wisf, &
          & psi((i - 1) * (Glr%wfd%nvctr_c + 7 * Glr%wfd%nvctr_f) + 1),psir)
     ! Compute partial densities coming from occup * psi_i * psi_i.
     do j = 1, Glr%d%n1i * Glr%d%n2i * Glr%d%n3i, 1
        rhor(j) = rhor(j) + orbs%occup(orbs%isorb + i) * psir(j) * psir(j)
     end do
  end do
  call mpiallred(rhor(1),Glr%d%n1i * Glr%d%n2i * Glr%d%n3i,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (iproc == 0) write(*,*) "System has", sum(rhor), "electrons."
  deallocate(rhor)
  ! Equivalent BigDFT routine.
  allocate(nscatterarr(0:nproc-1,4))
  allocate(ngatherarr(0:nproc-1,2))
  call createDensPotDescriptors(iproc,nproc,atoms,Glr%d, &
       & inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
       & rxyz,inputs%crmult,inputs%frmult,radii_cf,inputs%nspin,'D',inputs%ixc, &
       & inputs%rho_commun,n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
  allocate(rhor(Glr%d%n1i * Glr%d%n2i * n3d))
  allocate(irrzon(1,2,1))
  allocate(phnons(2,1,1))
  call sumrho(iproc,nproc,orbs,Glr,inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
       & psi,rhor,nscatterarr,inputs%nspin,GPU,atoms%symObj,irrzon,phnons,rhodsc)
  call deallocate_rho_descriptors(rhodsc,"main")

  ! Example of calculation of the energy of the local potential of the pseudos.
  call createKernel(iproc,nproc,atoms%geocode,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i, &
       & inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp,16,pkernel,quiet="yes")
  allocate(pot_ion(Glr%d%n1i * Glr%d%n2i * n3p))
  call createIonicPotential(atoms%geocode,iproc,nproc,atoms,rxyz,&
       & inputs%hx / 2._gp,inputs%hy / 2._gp,inputs%hz / 2._gp, &
       & inputs%elecfield,Glr%d%n1,Glr%d%n2,Glr%d%n3, &
       & n3pi,i3s+i3xcsh,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i, &
       & pkernel,pot_ion,psoffset,0,.false.)
  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p, &
       & Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,inputs%nspin, &
       & Glr%d%n1i*Glr%d%n2i*n3d,0, &
       & orbs%norb,orbs%norbp,ngatherarr,pot_ion,potential)
  epot_sum = 0._dp
  do i = 1, orbs%norbp, 1
     call daub_to_isf(Glr,wisf, &
          & psi((i - 1) * (Glr%wfd%nvctr_c + 7 * Glr%wfd%nvctr_f) + 1),psir)
     do j = 1, Glr%d%n1i * Glr%d%n2i * Glr%d%n3i, 1
        epot_sum = epot_sum + psir(j) * potential(j) * psir(j)
     end do
  end do
  epot_sum = epot_sum * inputs%hx / 2._gp * inputs%hy / 2._gp * inputs%hz / 2._gp
  call free_full_potential(nproc,potential,"main")
  call mpiallred(epot_sum,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (iproc == 0) write(*,*) "System pseudo energy is", epot_sum, "Ht."
  deallocate(pot_ion)
  deallocate(psir)
  call deallocate_work_arrays_sumrho(wisf)

  ! Free allocated space.
  deallocate(rxyz_old)
  deallocate(psi)

  call deallocate_comms(comms,"main")
  call deallocate_wfd(Glr%wfd,"main")

  call deallocate_bounds(Glr%geocode,Glr%hybrid_on,Glr%bounds,"main")
  call deallocate_orbs(orbs,"main")
  call deallocate_atoms_scf(atoms,"main") 

  deallocate(rxyz)
  call deallocate_atoms(atoms,"main") 
  call free_input_variables(inputs)

  !wait all processes before finalisation
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program wvl
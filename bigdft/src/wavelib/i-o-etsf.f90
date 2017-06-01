!> @file
!! Routines to read NetCDF (ETSF) format
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining internal routines for etsf 
module internal_etsf
   implicit none


   contains


   subroutine etsf_error(error)
      use wrapper_MPI, only: MPI_COMM_WORLD
      use etsf_io_low_level

      implicit none

      type(etsf_io_low_error), intent(in) :: error
      integer :: ierr

      call etsf_warning(error)
      call MPI_ABORT(MPI_COMM_WORLD, ierr)
   END SUBROUTINE etsf_error

   subroutine etsf_warning(error)
      use etsf_io_low_level

      implicit none

      type(etsf_io_low_error), intent(in) :: error
      character(len=etsf_io_low_error_len)  :: errmess

      call etsf_io_low_error_to_str(errmess, error)
      write(0,"(A)") trim(errmess)
   END SUBROUTINE etsf_warning

   subroutine etsf_read_astruct(ncid, astruct, lstat, error_data)
     use module_atoms
     use etsf_io_low_level
     use etsf_io
     use etsf_io_file
     use etsf_io_tools
     use dictionaries
     use dynamic_memory
     implicit none
     integer, intent(in) :: ncid
     type(atomic_structure), intent(out) :: astruct
     logical, intent(out) :: lstat
     type(etsf_io_low_error), intent(out) :: error_data
     
     integer :: i
     type(etsf_dims) :: dims
     double precision, dimension(3,3) :: rprimd
     character(len = etsf_charlen), dimension(:), allocatable :: atom_names

     call nullify_atomic_structure(astruct)

     call etsf_io_file_check_crystallographic_data(ncid, lstat, error_data)
     if (.not. lstat) return

     call etsf_io_dims_get(ncid, dims, lstat, error_data)
     if (.not. lstat) return

     call astruct_set_n_atoms(astruct, dims%number_of_atoms)
     call astruct_set_n_types(astruct, dims%number_of_atom_species)
     astruct%units = 'bohr'
     astruct%geocode = "P"
     
     call etsf_io_low_read_var(ncid, "primitive_vectors", &
          &   rprimd, lstat, error_data = error_data)
     if (.not. lstat) return
     astruct%cell_dim(1) = rprimd(1,1)
     astruct%cell_dim(2) = rprimd(2,2)
     astruct%cell_dim(3) = rprimd(3,3)

     call etsf_io_low_read_var(ncid, "atom_species", &
          &   astruct%iatype, lstat, error_data = error_data)
     if (.not. lstat) return

     call etsf_io_low_read_var(ncid, "reduced_atom_positions", &
          &   astruct%rxyz, lstat, error_data = error_data)
     if (.not. lstat) return
     do i = 1, astruct%nat
        astruct%rxyz(:, i) = astruct%rxyz(:, i) * astruct%cell_dim
     end do

     atom_names = f_malloc_str(etsf_charlen, (/ astruct%ntypes /), "atom_names")
     call etsf_io_tools_get_atom_names(ncid, atom_names, lstat, error_data)
     if (lstat) then
        do i = 1, astruct%ntypes
           write(astruct%atomnames(i), "(A)") trim(adjustl(atom_names(i)))
        end do
     end if
     call f_free_str(etsf_charlen, atom_names)
     if (.not. lstat) return

     lstat = .true.
   end subroutine etsf_read_astruct

   subroutine etsf_read_descr(ncid, orbsd, n1_old, n2_old, n3_old, hx_old, hy_old, hz_old, &
         &   lstat, error, nvctr_old, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
      use module_base
      use module_types
      use etsf_io_low_level
      use etsf_io

      implicit none

      integer, intent(in) :: ncid
      type(orbitals_data), intent(out) :: orbsd
      integer, intent(out) :: n1_old, n2_old, n3_old
      real(gp), intent(out) :: hx_old, hy_old, hz_old
      logical, intent(out) :: lstat
      type(etsf_io_low_error), intent(out) :: error
      ! Optional arguments
      integer, pointer, optional :: nvctr_old(:)
      integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
      integer, intent(in), optional :: nat
      real(gp), dimension(:,:), intent(out), optional :: rxyz_old

      character(len = *), parameter :: subname = "etsf_read_descr"
      type(etsf_dims) :: dims
      real(dp) :: rprimd(3,3)
      integer :: i, iat

      call etsf_io_dims_get(ncid, dims, lstat, error)
      if (.not. lstat) return
      ! The number of grid steps.
      n1_old = dims%number_of_grid_points_vector1
      n2_old = dims%number_of_grid_points_vector2
      n3_old = dims%number_of_grid_points_vector3
      ! The hgrid parameters.
      call etsf_io_low_read_var(ncid, "primitive_vectors", &
         &   rprimd, lstat, error_data = error)
      if (.not. lstat) return
      hx_old = rprimd(1,1) / n1_old
      hy_old = rprimd(2,2) / n2_old
      hz_old = rprimd(3,3) / n3_old
      n1_old = n1_old - 1
      n2_old = n2_old - 1
      n3_old = n3_old - 1
      ! We read the eigenvalues & occupations.
      orbsd%eval = f_malloc_ptr(dims%number_of_spins * dims%max_number_of_states * &
         &   dims%number_of_kpoints,id='orbsd%eval')
      !!$    allocate(orbsd%occup(dims%number_of_spins * dims%max_number_of_states * &
      !!$         & dims%number_of_kpoints),stat=i_stat)
      !!$    call memocc(i_stat,orbsd%occup,'orbsd%occup',subname)
      call etsf_io_low_read_var(ncid, "eigenvalues", &
         &   orbsd%eval, lstat, error_data = error)
      if (.not. lstat) return
      !!$    call etsf_io_low_read_var(ncid, "occupations", &
      !!$         & orbsd%occup, lstat, error_data = error)
      !!$    if (.not. lstat) call etsf_error(error)
      ! The orbitals description as on disk.
      orbsd%nspin = dims%number_of_spins
      orbsd%norbu = 0
      orbsd%norbd = 0
      do i = 1, dims%max_number_of_states, 1
         if (orbsd%eval(i) /= UNINITIALIZED(1.d0)) orbsd%norbu = orbsd%norbu + 1
         if (dims%number_of_spins > 1) then
            if (orbsd%eval(i + dims%max_number_of_states * dims%number_of_kpoints) /= &
               &   UNINITIALIZED(1.d0)) orbsd%norbd = orbsd%norbd + 1
         end if
      end do
      orbsd%norb = orbsd%norbu + orbsd%norbd
      orbsd%nspinor = dims%number_of_spinor_components
      orbsd%nkpts = dims%number_of_kpoints
      ! Put back the evals as sorted in BigDFT.
      call sortEvals(orbsd)

      ! Additional information read from the file.
      if (present(nvctr_old) .and. present(nvctr_c_old) .and. present(nvctr_f_old)) then
         ! The number of coarse and fine grid points.
         nvctr_old = f_malloc_ptr(dims%max_number_of_basis_grid_points,id='nvctr_old')
         call etsf_io_low_read_var(ncid, "number_of_coefficients_per_grid_point", &
            &   nvctr_old, lstat, error_data = error)
         if (.not. lstat) return
         nvctr_c_old = dims%max_number_of_basis_grid_points
         nvctr_f_old = 0
         do i = 1, dims%max_number_of_basis_grid_points, 1
            if (nvctr_old(i) == 8) nvctr_f_old = nvctr_f_old + 1
         end do
      end if
      if (present(nat) .And. present(rxyz_old)) then
         ! Sanity checks
         if (dims%number_of_atoms /= nat) call general_error("Mismatch in number of atoms")
         if (size(rxyz_old, 2) /= nat) call general_error("Mismatch in coordinate array size")
         ! The old atomic coordinates.
         call etsf_io_low_read_var(ncid, "reduced_atom_positions", &
            &   rxyz_old, lstat, error_data = error)
         if (.not. lstat) return
         do iat = 1, nat, 1
            rxyz_old(1, iat) = rxyz_old(1, iat) * rprimd(1, 1)
            rxyz_old(2, iat) = rxyz_old(2, iat) * rprimd(2, 2)
            rxyz_old(3, iat) = rxyz_old(3, iat) * rprimd(3, 3)
         end do
      end if

      contains

      subroutine general_error(error)
         character(len = *), intent(in) :: error

         integer :: ierr

         write(0,"(A)") error
         call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
      END SUBROUTINE general_error

      subroutine sortEvals(orbsd)
         type(orbitals_data), intent(inout) :: orbsd

         integer :: i, ik, ikd, isd
         real(wp), dimension(:), allocatable :: eval

         eval = f_malloc(size(orbsd%eval),id='eval')
         ! We transfer the eigenvalues & occupations.
         isd = max(orbsd%norbu, orbsd%norbd) * orbsd%nkpts
         do i = 1, orbsd%nkpts, 1
            ik = (i - 1) * orbsd%norb
            ikd = (i - 1) * max(orbsd%norbu, orbsd%norbd)
            eval(ik + 1:ik + orbsd%norbu) = orbsd%eval(ikd + 1:ikd + orbsd%norbu)
            if (orbsd%nspin > 1) then
               eval(ik + orbsd%norbu + 1:ik + orbsd%norb) = &
                  &   orbsd%eval(isd + ikd + 1:isd + ikd + orbsd%norbd)
            end if
         end do
         orbsd%eval = eval
         call f_free(eval)
      END SUBROUTINE sortEvals
   END SUBROUTINE etsf_read_descr

   subroutine etsf_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr, gcoord, lr)
      use module_defs
      use module_types
      use dynamic_memory

      implicit none

      integer, intent(in) :: n1, n2, n3, nvctr_c
      integer, dimension(nvctr_c), intent(in) :: nvctr
      integer, dimension(3, nvctr_c), intent(in) :: gcoord
      type(locreg_descriptors), intent(out) :: lr

      character(len = *), parameter :: subname = "etsf_gcoordToLocreg"
      logical, dimension(:,:,:), allocatable :: logrid_c, logrid_f
      integer :: i

      lr%geocode = "P"
      lr%hybrid_on = .false.

      lr%ns1 = 0
      lr%ns2 = 0
      lr%ns3 = 0

      lr%d%n1 = n1
      lr%d%n2 = n2
      lr%d%n3 = n3

      lr%d%n1i = 2 * n1 + 2
      lr%d%n2i = 2 * n2 + 2
      lr%d%n3i = 2 * n3 + 2

      logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
      logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')

      lr%d%nfl1 = n1
      lr%d%nfl2 = n2
      lr%d%nfl3 = n3
      lr%d%nfu1 = 0
      lr%d%nfu2 = 0
      lr%d%nfu3 = 0

      logrid_c(:,:,:) = .false.
      logrid_f(:,:,:) = .false.
      do i = 1, nvctr_c, 1
         logrid_c(gcoord(1, i), gcoord(2, i), gcoord(3, i)) = .true.
         if (nvctr(i) == 8) then
            logrid_f(gcoord(1, i), gcoord(2, i), gcoord(3, i)) = .true.
            lr%d%nfl1 = min(lr%d%nfl1, gcoord(1, i))
            lr%d%nfl2 = min(lr%d%nfl2, gcoord(2, i))
            lr%d%nfl3 = min(lr%d%nfl3, gcoord(3, i))
            lr%d%nfu1 = max(lr%d%nfu1, gcoord(1, i))
            lr%d%nfu2 = max(lr%d%nfu2, gcoord(2, i))
            lr%d%nfu3 = max(lr%d%nfu3, gcoord(3, i))
         end if
      end do

      !correct the values of the delimiter if there are no wavelets
      if (lr%d%nfl1 == n1 .and. lr%d%nfu1 == 0) then
         lr%d%nfl1 = n1 / 2
         lr%d%nfu1 = n1 / 2
      end if
      if (lr%d%nfl2 == n2 .and. lr%d%nfu2 == 0) then
         lr%d%nfl2 = n2 / 2
         lr%d%nfu2 = n2 / 2
      end if
      if (lr%d%nfl3 == n3 .and. lr%d%nfu3 == 0) then
         lr%d%nfl3 = n3 / 2
         lr%d%nfu3 = n3 / 2
      end if

      call wfd_from_grids(logrid_c, logrid_f, .true., lr)

      call f_free(logrid_c)
      call f_free(logrid_f)
   END SUBROUTINE etsf_gcoordToLocreg

   subroutine etsf_orbsToStartCount(start, count, iorbp, orbs)
      use module_base
      use module_types

      implicit none

      integer, intent(out) :: start(6), count(6)
      integer, intent(in) :: iorbp
      type(orbitals_data), intent(in) :: orbs
      !  integer,dimension(orbs%norb), intent(in), optional :: orblist

      start(:) = 1
      count(:) = 0
      ! Read/Write one spinor.
      start(3) = modulo(iorbp - 1, orbs%nspinor) + 1
      count(3) = 1
      ! Read/Write one orbital.
      start(4) = modulo(orbs%isorb + (iorbp - 1) / orbs%nspinor, orbs%norb) + 1
      count(4) = 1
      ! Read/Write one kpoint.
      start(5) = (orbs%isorb + (iorbp - 1) / orbs%nspinor) / orbs%norb + 1
      count(5) = 1
      ! Read/Write one spin.
      start(6) = 1
      if (start(4) > orbs%norbu) then
         start(6) = 2
         start(4) = start(4) - orbs%norbu
      end if
      count(6) = 1
   END SUBROUTINE etsf_orbsToStartCount

END MODULE internal_etsf


!> Read a ETSF (NETCDF) file containing wavefunctions.
!! Only import the given iorbp (processor-related number), in a compress form.
subroutine read_psi_compress_etsf(ncid, iorbp, orbs, nvctr, wfd, psi, lstat, error)
   use module_base
   use module_types
   use etsf_io_low_level
   use internal_etsf

   implicit none

   integer, intent(in) :: iorbp, ncid
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(orbitals_data), intent(in) :: orbs
   integer, dimension(wfd%nvctr_c), intent(in) :: nvctr
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
   type(etsf_io_low_error), intent(out) :: error
   logical, intent(out) :: lstat

   integer :: iFine, iCoeff, iGrid, diGrid
   integer :: start(6), count(6)

   ! We read the coefficients.
   call etsf_orbsToStartCount(start, count, iorbp, orbs)

   iFine = wfd%nvctr_c + 1
   iCoeff = 1
   iGrid = 1
   do
      if (iGrid > wfd%nvctr_c) exit
      diGrid = 0
      do
         if (nvctr(iGrid + diGrid) /= 1 .or. &
            &   iGrid + diGrid == wfd%nvctr_c) exit
         diGrid = diGrid + 1
      end do
      ! Read diGrid + 1 coeff.
      start(2) = iCoeff
      count(2) = diGrid + 1
      call etsf_io_low_read_var(ncid, "coefficients_of_wavefunctions", &
         &   psi(iGrid:iGrid + diGrid), lstat, error_data = error, start = start, count = count)
      if (.not. lstat) return
      iCoeff  = iCoeff + diGrid + 1

      if (nvctr(iGrid + diGrid) == 8) then
         ! Read seven coeff.
         start(2) = iCoeff
         count(2) = 7
         call etsf_io_low_read_var(ncid, "coefficients_of_wavefunctions", &
            &   psi(iFine:iFine+6), lstat, error_data = error, start = start, count = count)
         if (.not. lstat) return
         iCoeff = iCoeff + 7
         iFine  = iFine + 7
      end if
      iGrid = iGrid + diGrid + 1
   end do
END SUBROUTINE read_psi_compress_etsf

!>   Read a ETSF (NETCDF) file containing wavefunctions.
!!    Only import the given iorbp (processor-related number), in a grid form.
subroutine read_psi_full_etsf(ncid, iorbp, orbs, n1, n2, n3, &
      &   nvctr_c, nvctr, gcoord, psig, lstat, error)
   use module_base
   use module_types
   use etsf_io_low_level
   use internal_etsf

   implicit none

   integer, intent(in) :: iorbp, n1, n2, n3, nvctr_c, ncid
   type(orbitals_data), intent(in) :: orbs
   real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
   integer, dimension(3,nvctr_c), intent(in) :: gcoord
   integer, dimension(nvctr_c), intent(in) :: nvctr
   type(etsf_io_low_error), intent(out) :: error
   logical, intent(out) :: lstat

   integer :: i, iCoeff
   integer :: start(6), count(6), coord(3)
   real(wp) :: fv(7)

   ! We read the coefficients.
   call etsf_orbsToStartCount(start, count, iorbp, orbs)

   ! We transfer the coefficients in psig.
   iCoeff = 1
   do i = 1, nvctr_c, 1
      coord = gcoord(:, i)
      start(2) = iCoeff
      count(2) = 1
      call etsf_io_low_read_var(ncid, "coefficients_of_wavefunctions", &
         &   psig(coord(1), 1, coord(2), 1, coord(3), 1), &
         & lstat, error_data = error, start = start, count = count)
      if (.not. lstat) return
      iCoeff = iCoeff + 1
      if (nvctr(i) == 8) then
         start(2) = iCoeff
         count(2) = 7
         call etsf_io_low_read_var(ncid, "coefficients_of_wavefunctions", &
            &   fv, lstat, error_data = error, start = start, count = count)
         if (.not. lstat) return
         psig(coord(1), 2, coord(2), 1, coord(3), 1) = fv(1)
         psig(coord(1), 1, coord(2), 2, coord(3), 1) = fv(2)
         psig(coord(1), 2, coord(2), 2, coord(3), 1) = fv(3)
         psig(coord(1), 1, coord(2), 1, coord(3), 2) = fv(4)
         psig(coord(1), 2, coord(2), 1, coord(3), 2) = fv(5)
         psig(coord(1), 1, coord(2), 2, coord(3), 2) = fv(6)
         psig(coord(1), 2, coord(2), 2, coord(3), 2) = fv(7)
         iCoeff = iCoeff + 7
      end if
   end do
END SUBROUTINE read_psi_full_etsf


subroutine read_waves_from_list_etsf(iproc,filename,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz, & 
      &   wfd,psi,norb,nspinor,iorbparr,isorb,eval)
   use module_base
   use module_types
   use etsf_io_low_level
   use etsf_io
   use internal_etsf
   use box

   implicit none

   integer, intent(in) :: iproc,n1,n2,n3,norb,isorb,nspinor
   real(gp), intent(in) :: hx,hy,hz
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norb), intent(out) :: psi
   real(wp), dimension(norb), intent(out) :: eval
   integer, dimension(norb*nspinor), intent(in) :: iorbparr
   character(len=*), intent(in) :: filename
   ! Local variables
   character(len = *), parameter :: subname = "read_waves_from_list_etsf"
   integer, pointer :: nvctr_old(:)
   integer :: n1_old, n2_old, n3_old, nvctr_c_old, nvctr_f_old, ncid, iorb, iat
   integer :: nb1, nb2, nb3, ispinor
   real(gp) :: hx_old, hy_old, hz_old
   real(gp) :: displ
   logical :: perx, pery, perz
   integer, dimension(:,:), allocatable :: gcoord
   real(wp), dimension(:,:,:), allocatable :: psifscf
   real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
   type(orbitals_data) :: orbsd
   type(etsf_io_low_error) :: error
   logical :: lstat
   type(cell) :: mesh

   mesh=cell_new(at%astruct%geocode,[n1,n2,n3],[hx,hy,hz])

   ! We open the ETSF file
   call etsf_io_low_open_read(ncid, filename, lstat, error_data = error)
   if (.not. lstat) call etsf_error(error)

   ! We read the basis set description and the atomic definition.
   call etsf_read_descr(ncid, orbsd, n1_old, n2_old, n3_old, hx_old, hy_old, hz_old, &
      &   lstat, error, nvctr_old, nvctr_c_old, nvctr_f_old, rxyz_old, at%astruct%nat)
   if (.not. lstat) call etsf_error(error)
   orbsd%isorb = isorb

   !conditions for periodicity in the three directions
   perx=(at%astruct%geocode /= 'F')
   pery=(at%astruct%geocode == 'P')
   perz=(at%astruct%geocode /= 'F')

   displ=0.0_gp
   do iat=1,at%astruct%nat
      displ=displ+minimum_distance(mesh,rxyz(:,iat),rxyz_old(:,iat))**2
   enddo
   displ=sqrt(displ)

   if (abs(hx_old - hx) < 1e-6 .and. abs(hy_old - hy) < 1e-6 .and. abs(hz_old - hz) < 1e-6 .and. &
      &   nvctr_c_old == wfd%nvctr_c .and. nvctr_f_old == wfd%nvctr_f .and. & 
      & n1_old == n1 .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then
      if (iproc == 0) write(*,*) 'wavefunctions need NO reformatting'

      do iorb = 1, norb, 1
         do ispinor = 1, nspinor, 1
            call read_psi_compress_etsf(ncid, iorbparr(nspinor * (iorb - 1) + ispinor), &
               &   orbsd, nvctr_old, wfd, psi(1, ispinor, iorb), lstat, error)
            if (.not. lstat) call etsf_error(error)
         end do
      end do
   else
      if (iproc == 0) then
         write(*,*) 'wavefunctions need reformatting'
         if (abs(hx_old - hx) > 1e-6 .and. abs(hy_old - hy) > 1e-6 .and. abs(hz_old - hz) > 1e-6) &
            &   write(*,*) 'because hgrid_old /= hgrid',hx_old,hy_old,hz_old,hx,hy,hz
         if (nvctr_c_old /= wfd%nvctr_c) &
            &   write(*,*) 'because nvctr_c_old /= nvctr_c',nvctr_c_old,wfd%nvctr_c
         if (nvctr_f_old /= wfd%nvctr_f) &
            &   write(*,*) 'because nvctr_f_old /= nvctr_f',nvctr_f_old,wfd%nvctr_f
         if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
            &   write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
         if (displ > 1.d-3 ) &
            &   write(*,*) 'because of large displacement of molecule'
      end if

      ! We read the coordinates of grid points.
      gcoord = f_malloc((/ 3, nvctr_c_old /),id='gcoord')
      call etsf_io_low_read_var(ncid, "coordinates_of_basis_grid_points", &
         &   gcoord, lstat, error_data = error)
      if (.not. lstat) call etsf_error(error)

      !buffers realted to periodicity
      !WARNING: the boundary conditions are not assumed to change between new and old
      call ext_buffers_coarse(perx,nb1)
      call ext_buffers_coarse(pery,nb2)
      call ext_buffers_coarse(perz,nb3)

      psifscf = f_malloc((/ -nb1.to.2*n1+1+nb1, -nb2.to.2*n2+1+nb2, -nb3.to.2*n3+1+nb3 /),id='psifscf')

      psigold = f_malloc0((/ 0.to.n1_old, 1.to.2, 0.to.n2_old, 1.to.2, 0.to.n3_old, 1.to.2 /),id='psigold')
      !call to_zero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)

      do iorb = 1, norb, 1
         do ispinor = 1, nspinor, 1
            call read_psi_full_etsf(ncid, iorbparr(nspinor * (iorb - 1) + ispinor), &
               &   orbsd, n1_old, n2_old, n3_old, nvctr_c_old, nvctr_old, gcoord, psigold, &
               & lstat, error)
            if (.not. lstat) call etsf_error(error)

            call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,&
               &   rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi(1,ispinor,iorb))
         end do
      end do

      call f_free(psigold)
      call f_free(gcoord)
      call f_free(psifscf)
   end if

   call f_free_ptr(nvctr_old)

   ! We transfer the eval values.
   do iorb = 1, norb, 1
      eval(iorb) = orbsd%eval(orbsd%isorb + (iorbparr(nspinor * (iorb - 1) + 1) - 1) / nspinor + 1)
   end do

   call f_free_ptr(orbsd%eval)

   ! We close the file.
   call etsf_io_low_close(ncid, lstat, error)
   if (.not. lstat) call etsf_error(error)
END SUBROUTINE read_waves_from_list_etsf

subroutine readwavedescr_etsf(lstat, filename, norbu, norbd, nkpt, nspinor)
   use module_base
   use module_types
   use etsf_io_low_level
   use etsf_io_file
   use internal_etsf

   implicit none

   character(len = *), intent(in) :: filename
   integer, intent(out) :: norbu, norbd, nkpt, nspinor
   logical, intent(out) :: lstat

   integer :: ncid, n1, n2, n3
   real(gp) :: hx, hy, hz
   type(orbitals_data) :: orbsd
   type(etsf_io_low_error) :: error

   ! We open the ETSF file
   call etsf_io_low_open_read(ncid, filename, lstat, error_data = error)
   if (.not. lstat) then
      call etsf_warning(error)
      return
   end if

   ! Check that we are a valid wavefunction file.
   call etsf_io_file_check_wavefunctions_data(ncid, lstat, error_data = error)
   if (.not. lstat) then
      call etsf_warning(error)
      return
   end if

   ! We read the basis set description and the atomic definition.
   call etsf_read_descr(ncid, orbsd, n1, n2, n3, hx, hy, hz, lstat, error)
   if (.not. lstat) then
      call etsf_warning(error)
      return
   end if

   ! We close the ETSF file.
   call etsf_io_low_close(ncid, lstat, error)
   if (.not. lstat) call etsf_warning(error)

   if (associated(orbsd%eval)) then
      call f_free_ptr(orbsd%eval)
   end if

   norbu   = orbsd%norbu
   norbd   = orbsd%norbd
   nkpt    = orbsd%nkpts
   nspinor = orbsd%nspinor
END SUBROUTINE readwavedescr_etsf

subroutine readwavetoisf_etsf(lstat, filename, iorbp, hx, hy, hz, &
      &   n1, n2, n3, nspinor, psiscf)
   use module_base
   use module_types
   use etsf_io_low_level
   use etsf_io
   use internal_etsf
   use locreg_operations

   implicit none

   character(len = *), intent(in) :: filename
   integer, intent(in) :: iorbp
   integer, intent(out) :: n1, n2, n3, nspinor
   real(gp), intent(out) :: hx, hy, hz
   real(wp), dimension(:,:,:,:), pointer :: psiscf
   logical, intent(out) :: lstat

   integer :: ncid, ispinor, nvctr_c, nvctr_f
   integer, dimension(:,:), allocatable :: gcoord
   integer, dimension(:), pointer :: nvctr
   real(wp), dimension(:), allocatable :: psi
   type(orbitals_data) :: orbsd
   type(locreg_descriptors) :: lr
   type(etsf_io_low_error) :: error
   type(workarr_sumrho) :: w

   call nullify_locreg_descriptors(lr)

   ! We open the ETSF file
   call etsf_io_low_open_read(ncid, filename, lstat, error_data = error)
   if (.not. lstat) then
      call etsf_warning(error)
      return
   end if

   ! We read the basis set description and the atomic definition.
   call etsf_read_descr(ncid, orbsd, n1, n2, n3, hx, hy, hz, lstat, error, &
      &   nvctr, nvctr_c, nvctr_f)
   if (.not. lstat) then
      call etsf_warning(error)
      return
   end if
   nspinor = orbsd%nspinor
   orbsd%isorb = 0

   ! Initial allocations.
   gcoord = f_malloc((/ 3, nvctr_c  /),id='gcoord')
   psi = f_malloc(nvctr_c + 7 * nvctr_f ,id='psi')

   call etsf_io_low_read_var(ncid, "coordinates_of_basis_grid_points", &
      &   gcoord, lstat, error_data = error)
   if (.not. lstat) then
      call etsf_warning(error)
      call deallocate_local()
      return
   end if
   call etsf_gcoordToLocreg(n1, n2, n3, nvctr_c, nvctr, gcoord, lr)

   call f_free(gcoord)

   psiscf = f_malloc_ptr((/ lr%d%n1i, lr%d%n2i, lr%d%n3i, orbsd%nspinor  /),id='psiscf')

   call initialize_work_arrays_sumrho(1, (/ lr /),.true.,w)

   do ispinor = 1, orbsd%nspinor, 1
      call read_psi_compress_etsf(ncid, orbsd%nspinor * (iorbp - 1) + ispinor, &
         &   orbsd, nvctr, lr%wfd, psi, lstat, error)
      if (.not. lstat) then
         call etsf_warning(error)
         call deallocate_local()
         return
      end if
      call daub_to_isf(lr, w, psi, psiscf(1,1,1,ispinor))
   end do

   ! We update the size values to match the allocation of psiscf.
   n1 = lr%d%n1i
   n2 = lr%d%n2i
   n3 = lr%d%n3i
   hx = hx * 0.5d0
   hy = hy * 0.5d0
   hz = hz * 0.5d0

   call deallocate_local()

   contains

   subroutine deallocate_local()
     use locregs, only: deallocate_locreg_descriptors
      character(len = *), parameter :: subname = "read_wave_to_isf_etsf"

      ! We close the ETSF file.
      call etsf_io_low_close(ncid, lstat, error)
      if (.not. lstat) call etsf_warning(error)

      ! Final deallocations.
      call f_free_ptr(nvctr)
      call f_free_ptr(orbsd%eval)
      call f_free(psi)
      call f_free(gcoord)
      call deallocate_locreg_descriptors(lr)
      call deallocate_work_arrays_sumrho(w)
   END SUBROUTINE deallocate_local

END SUBROUTINE readwavetoisf_etsf


!>   Read a ETSF (NETCDF) file containing wavefunctions.
!!    coordinates_of_grid_points is used to store the geometric
!!   position of coefficients of wavelets i, as integer in
!!   dtset%wvl%ni(:) dimensions.
!!   coefficients_of_wavefunctions is used to store the psi values for
!!   each wavelet.
subroutine read_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
      &   wfd,psi)
   use module_base
   use module_types

   implicit none

   integer, intent(in) :: iproc,n1,n2,n3
   real(gp), intent(in) :: hx,hy,hz
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(orbitals_data), intent(inout) :: orbs
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(out) :: psi
   character(len = *), intent(in) :: filename
   ! Local variables
   integer :: i

   i = 0
   call read_waves_from_list_etsf(iproc,filename,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz, & 
      &   wfd,psi,orbs%norbp,orbs%nspinor,(/ (i, i=1, orbs%norbp*orbs%nspinor) /), &
      & orbs%isorb,orbs%eval(orbs%isorb + 1))
END SUBROUTINE read_waves_etsf

subroutine read_one_wave_etsf(iproc,filename,iorbp,isorb,nspinor,n1,n2,n3,&
      &   hx,hy,hz,at,rxyz_old,rxyz,wfd,psi,eval)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iorbp,iproc,n1,n2,n3,nspinor,isorb
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(atoms_data), intent(in) :: at
   real(gp), intent(in) :: hx,hy,hz
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), intent(out) :: eval
   real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor), intent(out) :: psi
   character(len = *), intent(in) :: filename
   real(wp), dimension(1) :: eval_

   if (nspinor == 1) then
      call read_waves_from_list_etsf(iproc,filename,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz, & 
           &   wfd,psi,1,nspinor,(/ iorbp /),isorb,eval_)
   else
      call read_waves_from_list_etsf(iproc,filename,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz, & 
         &   wfd,psi,1,nspinor,(/ 2 * iorbp - 1, 2 * iorbp /),isorb,eval_)
   end if
   eval = eval_(1)
END SUBROUTINE read_one_wave_etsf

subroutine write_psi_compress_etsf(ncid, iorbp, orbs, nvctr, wfd, psi)
   use module_base
   use module_types
   use etsf_io_low_level
   use internal_etsf

   implicit none

   integer, intent(in) :: iorbp, ncid
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(orbitals_data), intent(in) :: orbs
   integer, dimension(wfd%nvctr_c), intent(in) :: nvctr
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi

   integer :: iFine, iCoeff, iGrid, diGrid
   integer :: start(6), count(6)
   type(etsf_io_low_error) :: error
   logical :: lstat

   call etsf_orbsToStartCount(start, count, iorbp, orbs)

   ! iCoeff is the index of the coefficient we are writing in ETSF
   iCoeff  = 1
   ! iFine is the index of the fine part in psi
   iFine = wfd%nvctr_c + 1
   ! iGrid runs on all grid points.
   iGrid = 1
   do
      if (iGrid > wfd%nvctr_c) exit
      diGrid = 0
      do
         if (nvctr(iGrid + diGrid) /= 1 .or. iGrid + diGrid == wfd%nvctr_c) exit
         diGrid = diGrid + 1
      end do
      ! Write diGrid + 1 coeff.
      start(2) = iCoeff
      count(2) = diGrid + 1
      call etsf_io_low_write_var(ncid, "coefficients_of_wavefunctions", &
         &   psi(iGrid:iGrid + diGrid), lstat, error_data = error, &
         & start = start, count = count)
      if (.not. lstat) call etsf_error(error)
      iCoeff  = iCoeff + diGrid + 1

      if (nvctr(iGrid + diGrid) == 8) then
         ! Write seven coeff.
         start(2) = iCoeff
         count(2) = 7
         call etsf_io_low_write_var(ncid, "coefficients_of_wavefunctions", &
            &   psi(iFine:iFine+6), lstat, error_data = error, &
            & start = start, count = count)
         if (.not. lstat) call etsf_error(error)
         iCoeff = iCoeff + 7
         iFine  = iFine  + 7
      end if
      iGrid = iGrid + diGrid + 1
   end do
END SUBROUTINE write_psi_compress_etsf


!>   Write a ETSF file containing wavefunctions.
!!   Write a NetCDF file.
!!    coordinates_of_grid_points is used to store the geometric
!!   position of coefficients of wavelets i, as integer in
!!   (/ n1, n2, n3 /) dimensions.
!!   coefficients_of_wavefunctions is used to store the psi values for
!!   each wavelet.
subroutine write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
   use module_types
   use module_base

   use etsf_io_low_level
   use etsf_io

   use internal_etsf

   implicit none

   integer, intent(in) :: iproc,n1,n2,n3
   real(gp), intent(in) :: hx,hy,hz
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   type(wavefunctions_descriptors), intent(in) :: wfd
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
   character(len = *), intent(in) :: filename

   type(etsf_io_low_error) :: error
   logical :: lstat
   integer :: ncid, ierr, nproc
   integer :: ncount1,ncount_rate,ncount_max, ncount2, i, iorb
   real :: tr0,tr1
   integer, allocatable :: nvctr(:)
   integer, allocatable :: gcoord(:,:)
   real(gp) :: tel
   logical :: sequential
   character(len = *), parameter :: subname = "write_waves_etsf"

   call MPI_COMM_SIZE(bigdft_mpi%mpi_comm,nproc,ierr)

   ! nvctr array will contains the number of coeff per grid point,
   ! required by all processors.
   nvctr = f_malloc(wfd%nvctr_c,id='nvctr')
   gcoord = f_malloc((/ 3, wfd%nvctr_c /),id='gcoord')
   call build_grid(n1, n2, n3, nvctr, gcoord, wfd)

   !!$  sequential = .not. etsf_io_low_check_parallel_io()
   sequential = .true.

   ! Only the master proc create the file.
   if (iproc == 0) then
      call cpu_time(tr0)
      call system_clock(ncount1,ncount_rate,ncount_max)

      call etsf_io_low_open_create(ncid, filename // ".etsf", 1.3, lstat, &
         &   title = "BigDFT wavefunctions", error_data = error, &
         & overwrite = .true., with_etsf_header = .true.)
      if (.not. lstat) call etsf_error(error)

      call etsf_write_global(ncid,orbs, n1,n2,n3,hx,hy,hz,rxyz,at,wfd,gcoord,nvctr)

      ! We close the file.
      call etsf_io_low_close(ncid, lstat, error)
      if (.not. lstat) call etsf_error(error)
   end if

   call f_free(gcoord)

   ! Now that the file is created and writable, we call the writing routines.
   if (sequential) then
      do i = 0, iproc - 1, 1
         call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
      end do
      call etsf_io_low_open_modify(ncid, filename // ".etsf", lstat, error_data = error)
      if (.not. lstat) call etsf_error(error)
      !!$  else
      !!$     call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
      !!$     call etsf_io_low_open_modify(ncid, filename, lstat, error_data = error, &
      !!$          & mpi_comm = bigdft_mpi%mpi_comm, mpi_info = MPI_INFO_NULL)
      !!$     if (.not. lstat) call etsf_error(error)
   end if
   call etsf_io_low_set_write_mode(ncid, lstat, error)
   if (.not. lstat) call etsf_error(error)

   ! We run over a processor independant number of orbitals
   ! to ensure the synchronisation to disk (see later).
   do iorb = 1, (orbs%norb * orbs%nkpts / nproc + 1 ) * orbs%nspinor, 1
      if (iorb <= (orbs%norbp * orbs%nspinor)) then
         call write_psi_compress_etsf(ncid, iorb, orbs, nvctr, wfd, psi(1, iorb))
      end if
   end do

   call f_free(nvctr)

   call etsf_io_low_close(ncid, lstat, error)
   if (.not. lstat) call etsf_error(error)

   ! We wait for all procs to write their waves.
   if (sequential) then
      do i = iproc, nproc - 1, 1
         call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
      end do
   else
      call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   end if

   if (iproc == 0) then
      call cpu_time(tr1)
      call system_clock(ncount2,ncount_rate,ncount_max)
      tel=dble(ncount2-ncount1)/dble(ncount_rate)
      write(*,'(a,l1,a,2(1x,1pe10.3))') '- WRITING WAVES TIME (',sequential,')',tr1-tr0,tel
   end if

   contains

     subroutine etsf_write_global(ncid,orbs, n1,n2,n3,hx,hy,hz,rxyz,at,wfd,gcoord,nvctr)
       use ao_inguess, only: atomic_info
       implicit none
      integer, intent(in) :: ncid, n1, n2, n3
      real(gp), intent(in) :: hx, hy, hz
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(in) :: orbs
      type(wavefunctions_descriptors), intent(in) :: wfd
      real(gp), intent(in) :: rxyz(3,at%astruct%nat)
      integer, target, intent(in) :: nvctr(wfd%nvctr_c)
      integer, target, intent(in) :: gcoord(3,wfd%nvctr_c)

      type(etsf_dims) :: dims
      type(etsf_geometry) :: geo
      type(etsf_basisdata) :: basis
      type(etsf_kpoints) :: kpts
      type(etsf_electrons) :: elec
      integer :: iat, i, ispin, iorb
      double precision, target :: rprimd(3,3)
      double precision, allocatable, target :: xred(:,:)
      double precision, dimension(:), allocatable, target :: znucl
      character(len=etsf_chemlen), allocatable, dimension(:), target :: spnames
      integer, dimension(3,3,1), target :: symId
      real(gp), dimension(3,1), target :: transId
      character(len=etsf_charlen), target :: wvlBasis
      logical :: lstat
      type(etsf_io_low_error) :: error

      ! Unused dims, to be removed later
      dims%number_of_components = etsf_no_dimension
      dims%number_of_atom_species = etsf_no_dimension
      dims%max_number_of_angular_momenta = etsf_no_dimension
      dims%max_number_of_projectors = etsf_no_dimension
      dims%number_of_coefficients_dielectric_function = etsf_no_dimension
      dims%number_of_frequencies_dielectric_function = etsf_no_dimension
      dims%number_of_qpoints_dielectric_function = etsf_no_dimension
      dims%number_of_qpoints_gamma_limit = etsf_no_dimension
      dims%real_or_complex_gw_corrections = etsf_no_dimension
      dims%real_or_complex_wavefunctions = etsf_no_dimension
      dims%real_or_complex_density = etsf_no_dimension
      dims%real_or_complex_potential = etsf_no_dimension

      ! Specific dims of interest.
      dims%number_of_symmetry_operations = 1
      dims%real_or_complex_coefficients = 1
      dims%max_number_of_coefficients = wfd%nvctr_c+7*wfd%nvctr_f
      dims%number_of_spinor_components = orbs%nspinor
      dims%max_number_of_states = max(orbs%norbu, orbs%norbd)
      dims%number_of_kpoints = orbs%nkpts
      if (orbs%norbd > 0) then
         dims%number_of_spins = 2
      else
         dims%number_of_spins = 1
      end if

      dims%max_number_of_basis_grid_points = wfd%nvctr_c
      dims%number_of_localization_regions = 1

      dims%number_of_atoms               = at%astruct%nat
      dims%number_of_atom_species        = at%astruct%ntypes
      !!$    if (at%astruct%geocode == 'P') then
      dims%number_of_grid_points_vector1 = n1 + 1
      dims%number_of_grid_points_vector2 = n2 + 1
      dims%number_of_grid_points_vector3 = n3 + 1
      !!$    else if (at%astruct%geocode == 'S') then
      !!$       dims%number_of_grid_points_vector1 = n1 + 1
      !!$       dims%number_of_grid_points_vector2 = n2
      !!$       dims%number_of_grid_points_vector3 = n3 + 1
      !!$    else
      !!$       dims%number_of_grid_points_vector1 = n1
      !!$       dims%number_of_grid_points_vector2 = n2
      !!$       dims%number_of_grid_points_vector3 = n3
      !!$    end if

      ! We write the dimensions to the file.
      call etsf_io_dims_def(ncid, dims, lstat, error)
      if (.not. lstat) call etsf_error(error)

      ! We set up the required variables.
      call etsf_io_geometry_def(ncid, lstat, error, &
         &   flags = etsf_geometry_primitive_vectors + etsf_geometry_atom_species + &
         & etsf_geometry_red_at_pos + etsf_geometry_atomic_numbers + &
         &   etsf_geometry_chemical_symbols + etsf_geometry_red_sym_matrices + &
         & etsf_geometry_red_sym_trans)
      if (.not. lstat) call etsf_error(error)
      call etsf_io_basisdata_def(ncid, lstat, error, &
         &   flags = etsf_basisdata_basis_set + etsf_basisdata_coord_grid + &
         & etsf_basisdata_n_coeff_grid)
      if (.not. lstat) call etsf_error(error)
      call etsf_io_kpoints_def(ncid, lstat, error, &
         &   flags = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights)
      if (.not. lstat) call etsf_error(error)
      call etsf_io_electrons_def(ncid, lstat, error, &
         &   flags = etsf_electrons_number_of_states + etsf_electrons_eigenvalues + &
         & etsf_electrons_occupations)
      if (.not. lstat) call etsf_error(error)
      call etsf_io_main_def(ncid, lstat, error, flags = etsf_main_wfs_coeff)
      if (.not. lstat) call etsf_error(error)

      ! We write the global informations.
      ! Geometry
      rprimd = reshape((/ (hx * dims%number_of_grid_points_vector1),0.0_gp,0.0_gp, &
         &   0.0_gp,(hy * dims%number_of_grid_points_vector2),0.0_gp, &
         & 0.0_gp,0.0_gp,(hz * dims%number_of_grid_points_vector3) /), (/ 3, 3 /))
      xred = f_malloc((/ 3, at%astruct%nat /),id='xred')
      do iat = 1, at%astruct%nat, 1
         xred(:, iat) = rxyz(:, iat) / &
            &   (/ hx * dims%number_of_grid_points_vector1, &
            &    hy * dims%number_of_grid_points_vector2, &
            &   hz * dims%number_of_grid_points_vector3 /)
      end do
      znucl = f_malloc(at%astruct%ntypes,id='znucl')
      znucl = real(at%nzatom)
      spnames = f_malloc_str(at%astruct%ntypes,id='spnames')
      do iat = 1, at%astruct%ntypes, 1
         !call nzsymbol(at%nzatom(iat), spnames(iat))
         call atomic_info(at%nzatom(iat),at%nelpsp(iat),symbol=spnames(iat))
      end do
      symId = reshape((/1,0,0,0,1,0,0,0,1/), (/3,3,1/))
      transId = reshape((/0.d0,0.d0,0.d0/), (/3,1/))
      geo%chemical_symbols              => spnames
      geo%atom_species                  => at%astruct%iatype
      geo%atomic_numbers                => znucl
      geo%reduced_atom_positions        => xred
      geo%primitive_vectors             => rprimd
      geo%reduced_symmetry_matrices     => symId
      geo%reduced_symmetry_translations => transId
      call etsf_io_geometry_put(ncid, geo, lstat, error)
      if (.not. lstat) call etsf_error(error)
      call f_free(xred)
      call f_free(znucl)
      call f_free_str(etsf_chemlen,spnames)
      ! The eigenvalues & occupation.
      if (dims%number_of_spins == 1) then
         elec%eigenvalues%data1D => orbs%eval
         elec%occupations%data1D => orbs%occup
      else
         elec%eigenvalues%data3D = f_malloc_ptr((/ dims%max_number_of_states , &
             dims%number_of_kpoints , dims%number_of_spins /), &
             id='elec%eigenvalues%data3D')
         elec%eigenvalues%data3D = UNINITIALIZED(1.d0)
         elec%occupations%data3D = f_malloc_ptr((/ dims%max_number_of_states , &
             dims%number_of_kpoints , dims%number_of_spins /), &
             id='elec%occupations%data3D')
         elec%occupations%data3D = UNINITIALIZED(1.d0)
         do i = 1, orbs%norb*orbs%nkpts, 1
            ispin = 1
            iorb = modulo(i - 1, orbs%norb) + 1
            if (iorb > orbs%norbu) then
               ispin = 2
               iorb = iorb - orbs%norbu
            end if
            elec%eigenvalues%data3D(iorb, (i - 1) / orbs%norb + 1, ispin) = orbs%eval(i)
            elec%occupations%data3D(iorb, (i - 1) / orbs%norb + 1, ispin) = orbs%occup(i)
         end do
      end if
      elec%number_of_states%data2D = f_malloc_ptr((/ dims%number_of_kpoints , &
          dims%number_of_spins /),id='elec%number_of_states%data2D')
      do ispin = 1, dims%number_of_spins, 1
         do i = 1, orbs%nkpts, 1
            if (ispin == 1) then
               elec%number_of_states%data2D(i, ispin) = orbs%norbu
            else
               elec%number_of_states%data2D(i, ispin) = orbs%norbd
            end if
         end do
      end do
      call etsf_io_electrons_put(ncid, elec, lstat, error)
      if (.not. lstat) call etsf_error(error)
      if (dims%number_of_spins /= 1) then
         call f_free_ptr(elec%eigenvalues%data3D)
         call f_free_ptr(elec%occupations%data3D)
      end if
      call f_free_ptr(elec%number_of_states%data2D)
      ! Basis set
      write(wvlBasis, "(A)") "daubechies_wavelets"
      basis%coordinates_of_basis_grid_points%data2D => gcoord
      basis%number_of_coefficients_per_grid_point%data1D => nvctr
      basis%basis_set => wvlBasis
      call etsf_io_basisdata_put(ncid, basis, lstat, error)
      if (.not. lstat) call etsf_error(error)
      ! Kpoints
      kpts%reduced_coordinates_of_kpoints => orbs%kpts
      kpts%kpoint_weights => orbs%kwgts
      call etsf_io_kpoints_put(ncid, kpts, lstat, error)
      if (.not. lstat) call etsf_error(error)
   END SUBROUTINE etsf_write_global

   subroutine build_grid(n1,n2,n3,nvctr, gcoord, wfd)
      integer, intent(in) :: n1, n2, n3
      type(wavefunctions_descriptors), intent(in) :: wfd
      integer, intent(out) :: nvctr(wfd%nvctr_c)
      integer, intent(out) :: gcoord(3,wfd%nvctr_c)

      integer :: ii, i0, i1, i2, i3, jj, j0, j1, iGrid, i, iseg
      integer, allocatable :: coeff_map(:,:,:)

      ! Will store the grid index for a given geometric point
      coeff_map = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='coeff_map')
      ! coarse part
      coeff_map = 0
      do iseg = 1, wfd%nseg_c
         jj = wfd%keyvloc(iseg)
         j0 = wfd%keygloc(1, iseg)
         j1 = wfd%keygloc(2, iseg)
         ii = j0 - 1
         i3 = ii / ((n1 + 1) * (n2 + 1))
         ii = ii - i3 * (n1 + 1) * (n2 + 1)
         i2 = ii / (n1 + 1)
         i0 = ii - i2 * (n1 + 1)
         i1 = i0 + j1 - j0
         do i = i0, i1
            iGrid = i - i0 + jj
            coeff_map(i, i2, i3) = iGrid
            gcoord(:, iGrid) = (/ i, i2, i3 /)
            nvctr(iGrid) = 1
         end do
      end do
      ! fine part
      do iseg = 1, wfd%nseg_f
         jj = wfd%keyvloc(wfd%nseg_c + iseg)
         j0 = wfd%keygloc(1, wfd%nseg_c + iseg)
         j1 = wfd%keygloc(2, wfd%nseg_c + iseg)
         ii = j0 - 1
         i3 = ii / ((n1 + 1) * (n2 + 1))
         ii = ii - i3 * (n1 + 1) * (n2 + 1)
         i2 = ii / (n1 + 1)
         i0 = ii - i2 * (n1 + 1)
         i1 = i0 + j1 - j0
         do i = i0, i1
            iGrid = coeff_map(i , i2 , i3)
            nvctr(iGrid) = nvctr(iGrid) + 7
         end do
      end do

      call f_free(coeff_map)
   END SUBROUTINE build_grid
END SUBROUTINE write_waves_etsf


subroutine read_pw_waves(filename, iproc, nproc, at, rxyz, Glr, orbs, psig, rhoij)
  use internal_etsf, only: etsf_read_astruct
  use module_defs, only: gp, wp, pi_param
  use module_types, only: orbitals_data
  use locreg_operations
  use module_atoms
  use locregs
  use etsf_io_low_level
  use etsf_io
  use dynamic_memory
  use f_utils
  use dictionaries
  use yaml_strings, only: yaml_toa
  !use module_interfaces, only: plot_wf
  
  implicit none

  character(len = *), intent(in) :: filename
  integer, intent(in) :: iproc, nproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: Glr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, orbs%nspinor, orbs%norbp), intent(out) :: psig
  real(wp), dimension(:,:,:), pointer, optional :: rhoij

  integer :: ncid, ikpt, isp, iorb, ispinor, inzee, i, j, k, ind, io
  integer :: nd1, nd1b, nd2, nd3, nd3b, n1f, n3f, n1b, n3b, nd1f, nd3f
  type(etsf_io_low_error) :: error
  type(etsf_io_low_var_infos) :: varrhoij
  logical :: lstat
  type(etsf_dims) :: dims
  real(wp) :: ralpha, ialpha, nrm

  integer, dimension(3,3) :: identity
  real(wp), dimension(2) :: cplex
  integer, dimension(:), allocatable :: ncomp
  integer, dimension(:,:), allocatable :: gvect, norbs
  real(wp), dimension(:,:,:), allocatable :: cg, cfft
  real(wp), dimension(:), allocatable :: psi
  type(atomic_structure) :: at_old
  type(workarr_sumrho) :: w

  call nullify_atomic_structure(at_old)
  call f_zero(psig)
  call initialize_work_arrays_sumrho(1, (/ Glr /),.true.,w)
  
  ! We open the ETSF file
  call etsf_io_low_open_read(ncid, filename, lstat, error_data = error)
  if (.not. lstat) call etsf_error(error)

  call etsf_read_astruct(ncid, at_old, lstat, error)
  if (.not. lstat) call etsf_error(error)

  call etsf_io_dims_get(ncid, dims, lstat, error)
  if (.not. lstat) call etsf_error(error)

  ncomp = f_malloc(dims%number_of_kpoints, "ncomp")
  call etsf_io_low_read_var(ncid, "number_of_coefficients", &
       & ncomp, lstat, error_data = error)
  if (.not. lstat) call etsf_error(error)

  norbs = f_malloc((/ dims%number_of_kpoints, dims%number_of_spins /), "norbs")
  call etsf_io_low_read_var(ncid, "number_of_states", &
       & norbs, lstat, error_data = error)
  if (.not. lstat) call etsf_error(error)
  ! Checks on orbitals
  if (f_err_raise(orbs%nkpts /= dims%number_of_kpoints, &
       & "Number of kpoints not matching")) return
  if (f_err_raise(.not. all(sum(norbs, 2) == orbs%norb), &
       & "Number of orbitals not matching")) return

  identity(:,:)=0; identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

  ! We use the Glr box as the FFT box to expand.
  call dimensions_fft(Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, &
       & nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)
  cfft = f_malloc0((/ 2, nd1 * nd2 * nd3, 2 /), "cfft")
  psi = f_malloc((/ Glr%d%n1i * Glr%d%n2i * Glr%d%n3i /), "psi")

  do ikpt = 1, dims%number_of_kpoints
     if (.not. any(orbs%iokpt == ikpt)) cycle

     gvect = f_malloc((/ 3, ncomp(ikpt) /), "gvect")
     call etsf_io_low_read_var(ncid, "reduced_coordinates_of_plane_waves", &
          & gvect, lstat, error_data = error, &
          & start = (/ 1, 1, ikpt /), count = (/ 3, ncomp(ikpt), 1 /))
     if (.not. lstat) call etsf_error(error)

     io = (ikpt - 1) * orbs%norb
     cg = f_malloc((/ 2, ncomp(ikpt), dims%number_of_spinor_components /), "cg")
     do isp = 1, dims%number_of_spins
        do iorb = 1, norbs(ikpt, isp)
           io = io + 1
           if (io <= orbs%isorb .or. io > orbs%isorb + orbs%norbp) cycle

           call etsf_io_low_read_var(ncid, "coefficients_of_wavefunctions", &
                & cg, lstat, error_data = error, &
                & start = (/ 1, 1, 1, iorb, ikpt, isp /), count = (/ 0, ncomp(ikpt), 0, 1, 1, 1 /))
           if (.not. lstat) call etsf_error(error)
           
           do ispinor = 1, dims%number_of_spinor_components
              ! Convert here cg to real space.
              inzee = 1
              call sphere(cg(1, 1, ispinor), 1, ncomp(ikpt), cfft(1,1,inzee), &
                   & Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, nd1, nd2, nd3, gvect, 1, 1, 0, &
                   & (/ 0, 0, 0 /), identity, 1.d0)
              call FFT(Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, nd1, nd2, nd3, cfft, -1, inzee)
              cfft(:,:,inzee) = cfft(:,:,inzee) / &
                   & sqrt(real(Glr%d%n1i * Glr%d%n2i * Glr%d%n3i, wp))
              if (orbs%nspinor == 1 .or. maxval(orbs%kpts(:, ikpt)) == 0.) then
                 ! Remove complex rotation.
                 nrm = 0.d0
                 ralpha = 0.d0
                 ialpha = 0.d0
                 do k = 1, Glr%d%n3i
                    do j = 1, Glr%d%n2i
                       do i = 1, Glr%d%n1i
                          cplex = cfft(:, i + (j - 1) * nd1 + (k - 1) * nd1 * nd2, inzee)
                          nrm = nrm + cplex(1) ** 2 + cplex(2) ** 2
                          ralpha = ralpha + cplex(1) ** 2 - cplex(2) ** 2
                          ialpha = ialpha + 2._wp * cplex(1) * cplex(2)
                       end do
                    end do
                 end do
                 ralpha = ralpha / nrm
                 ialpha = ialpha / nrm
                 if (ralpha >= 0.d0 .and. ialpha >= 0.d0) then
                    ralpha = cos(0.5_wp * acos(ralpha))
                    ialpha = sin(0.5_wp * asin(ialpha))
                 else if (ralpha >= 0.d0 .and. ialpha < 0.d0) then
                    ralpha = cos(0.5_wp * acos(ralpha))
                    ialpha = sin(0.5_wp * asin(ialpha))
                 else if (ralpha < 0.d0 .and. ialpha >= 0.d0) then
                    ralpha = cos(0.5_wp * acos(ralpha))
                    ialpha = sin(0.5_wp * (pi_param - asin(ialpha)))
                 else if (ralpha < 0.d0 .and. ialpha < 0.d0) then
                    ralpha = cos(0.5_wp * acos(ralpha))
                    ialpha = sin(0.5_wp * (pi_param - asin(ialpha)))
                 end if
                 do k = 1, Glr%d%n3i
                    do j = 1, Glr%d%n2i
                       do i = 1, Glr%d%n1i
                          ind = i + (j - 1) * nd1 + (k - 1) * nd1 * nd2
                          psi(i + (j - 1) * Glr%d%n1i + (k - 1) * Glr%d%n1i * Glr%d%n2i) = &
                               & cfft(1, ind, inzee) * ralpha + cfft(2, ind, inzee) * ialpha
                       end do
                    end do
                 end do
                 ! Convert to daub
                 call isf_to_daub(Glr, w, psi(1), psig(1, 1, io - orbs%isorb))
              else
                 ! complex wf case
                 do k = 1, Glr%d%n3i
                    do j = 1, Glr%d%n2i
                       do i = 1, Glr%d%n1i
                          ind = i + (j - 1) * nd1 + (k - 1) * nd1 * nd2
                          psi(i + (j - 1) * Glr%d%n1i + (k - 1) * Glr%d%n1i * Glr%d%n2i) = &
                               & cfft(1, ind, inzee)
                       end do
                    end do
                 end do
                 call isf_to_daub(Glr, w, psi(1), psig(1, 1, io - orbs%isorb))
                 do k = 1, Glr%d%n3i
                    do j = 1, Glr%d%n2i
                       do i = 1, Glr%d%n1i
                          ind = i + (j - 1) * nd1 + (k - 1) * nd1 * nd2
                          psi(i + (j - 1) * Glr%d%n1i + (k - 1) * Glr%d%n1i * Glr%d%n2i) = &
                               & cfft(2, ind, inzee)
                       end do
                    end do
                 end do
                 call isf_to_daub(Glr, w, psi(1), psig(1, 2, io - orbs%isorb))
              end if
!!$              write(*,*) iorb, sum(psi(1,:)**2)
!!$              call plot_wf(.false.,"toto" // trim(yaml_toa(iorb, fmt = "(I2.2)")),&
!!$                   & 1,at,1._wp,Glr,2._gp * at_old%cell_dim(1) / real(Glr%d%n1i, gp),&
!!$                   & 2._gp * at_old%cell_dim(2) / real(Glr%d%n2i, gp),&
!!$                   & 2._gp * at_old%cell_dim(3) / real(Glr%d%n3i, gp),&
!!$                   & at_old%rxyz,psig(1, io))
           end do

        end do
     end do
     call f_free(cg)     
     call f_free(gvect)
  end do
  call f_free(psi)
  call f_free(cfft)

  ! Potentially read rhoij.
  if (present(rhoij)) then
     if(associated(rhoij)) then
        call f_free_ptr(rhoij)
        nullify(rhoij)
     end if
     varrhoij%ncshape = 0
     call etsf_io_low_read_var_infos(ncid, "rhoij", varrhoij, lstat)
     if (lstat .and. varrhoij%ncshape == 4) then
        rhoij = f_malloc_ptr((/ varrhoij%ncdims(1) * varrhoij%ncdims(2), &
             & varrhoij%ncdims(3), varrhoij%ncdims(4) /), id = "rhoij")
        call etsf_io_low_read_var(ncid, "rhoij", rhoij, lstat, error_data = error)
        if (.not. lstat) call etsf_error(error)
     end if
  end if
  
  ! We close the file.
  call etsf_io_low_close(ncid, lstat, error)
  if (.not. lstat) call etsf_error(error)

  call f_free(norbs)
  call f_free(ncomp)

  call deallocate_atomic_structure(at_old)
  call deallocate_work_arrays_sumrho(w)

contains
  
  subroutine etsf_error(error)
    use etsf_io_low_level
    use dictionaries, only: f_err_throw

    implicit none

    type(etsf_io_low_error), intent(in) :: error
    character(len=etsf_io_low_error_len)  :: errmess

    call deallocate_atomic_structure(at_old)
    call deallocate_work_arrays_sumrho(w)

    call etsf_io_low_error_to_str(errmess, error)
    call f_err_throw(errmess, err_name='BIGDFT_RUNTIME_ERROR')
  END SUBROUTINE etsf_error

end subroutine read_pw_waves


!!****f* m_fftcore/sphere
!! NAME
!! sphere
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! INPUTS
!! cg(2,npw*ndat)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! istwf_k=option parameter that describes the storage of wfs
!! iflag=option parameter. Possible values: -1, -2, 1, 2
!! me_g0=1 if this node has G=0.
!! shiftg(3)=The shift in reciprocal space.
!! symm(3,3)=symmetry operation in reciprocal space to be applied.
!! xnorm=Normalization factor.
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!! If istwf_k differs from 1, then special storage modes must be taken
!! into account, for symmetric wavefunctions coming from k=(0 0 0) or other
!! special k points.
!!
!! TODO
!! 1) Order arguments
!! 2) Split the two cases: from and to sphere (merge with cg_box2gpsh and cg_gsph2box?)
!!
!! PARENTS
!!      fourwf,m_fft,m_fftcore,m_fftw3,m_io_kss,wfconv
!!
!! CHILDREN
!!
!! SOURCE

subroutine sphere(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
  use module_defs, only: dp
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iflag,istwf_k,n1,n2,n3,n4,n5,n6,ndat,npw,me_g0
 real(dp),intent(in) :: xnorm
!arrays
 integer,intent(in) :: kg_k(3,npw),shiftg(3),symm(3,3)
 real(dp),intent(inout) :: cfft(2,n4,n5,n6*ndat),cg(2,npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i1inv,i2,i2inv,i3,i3inv,id1,id2,id3,idat,ipw
 integer :: j1,j2,j3,l1,l2,l3,npwmin,use_symmetry,i3dat,i3invdat,i2invdat,ipwdat,i2dat
 character(len=500) :: msg
!arrays
 integer :: identity(3,3)
 integer :: i1inver(n1),i2inver(n2),i3inver(n3)

! *************************************************************************

 !In the case of special k-points, invariant under time-reversal,
 !but not Gamma, initialize the inverse coordinates. !Remember indeed that
 !
 !  u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^* and therefore:
 !  u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.

 if (istwf_k>=2) then
   if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
     i1inver(1)=1
     do i1=2,n1
       i1inver(i1)=n1+2-i1
     end do
   else
     do i1=1,n1
       i1inver(i1)=n1+1-i1
     end do
   end if
   if(istwf_k>=2 .and. istwf_k<=5)then
     i2inver(1)=1
     do i2=2,n2
       i2inver(i2)=n2+2-i2
     end do
   else
     do i2=1,n2
       i2inver(i2)=n2+1-i2
     end do
   end if
   if(istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7)then
     i3inver(1)=1
     do i3=2,n3
       i3inver(i3)=n3+2-i3
     end do
   else
     do i3=1,n3
       i3inver(i3)=n3+1-i3
     end do
   end if
 end if

 if (iflag==1 .or. iflag==2) then
   ! Insert cg into cfft with extra 0 s around outside:
    cfft = 0.d0

   ! Take care of each plane wave, and complete cfft if needed
   if (istwf_k==1) then

     if (iflag==1) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

     if (iflag==2) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

   else if (istwf_k>=2) then

     npwmin=1
     if (istwf_k==2 .and. me_g0==1) then 
       ! If gamma point, then cfft must be completed
       do idat=1,ndat
         cfft(1,1,1,1+n6*(idat-1))=cg(1,1+npw*(idat-1))
         cfft(2,1,1,1+n6*(idat-1))=0.d0
       end do
       npwmin=2
     end if

     if (iflag==1) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
           ! Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i2inv,i3inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i2inv,i3inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

     if (iflag==2) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           ! Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i3inv,i2inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i3inv,i2inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

   end if

 else if (iflag==-1 .or. iflag==-2) then

   use_symmetry=0
   identity(:,:)=0; identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
   if(sum((symm(:,:)-identity(:,:))**2)/=0)use_symmetry=1
   if(sum(shiftg(:)**2)/=0)use_symmetry=1

   ! Extract cg from cfft, ignoring components outside range of cg:
   if (istwf_k==1) then

     if (use_symmetry==0) then
       if (iflag==-1) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,ipwdat,i3dat) IF (ndat>1)
         do idat=1,ndat
           do ipw=1,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
             ipwdat = ipw + (idat-1) * npw
             i3dat = i3 + (idat-1) * n6

             cg(1,ipwdat)=cfft(1,i1,i2,i3dat)*xnorm
             cg(2,ipwdat)=cfft(2,i1,i2,i3dat)*xnorm
           end do
         end do
       else
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,ipwdat,i2dat) IF (ndat>1)
         do idat=1,ndat
           do ipw=1,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ipwdat = ipw + (idat-1) * npw
             i2dat = i2 + (idat-1) * n6

             cg(1,ipwdat)=cfft(1,i1,i3,i2dat)*xnorm
             cg(2,ipwdat)=cfft(2,i1,i3,i2dat)*xnorm
           end do
         end do
       end if
     else
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,j1,j2,j3,l1,l2,l3,ipwdat,i3dat) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           l1=kg_k(1,ipw)+shiftg(1)
           l2=kg_k(2,ipw)+shiftg(2)
           l3=kg_k(3,ipw)+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1; i1=j1+1
           if(j2<0)j2=j2+n2; i2=j2+1
           if(j3<0)j3=j3+n3; i3=j3+1

           ipwdat = ipw + (idat-1) * npw
           i3dat = i3 + (idat-1)*n6

           cg(1,ipwdat)=cfft(1,i1,i2,i3dat)*xnorm
           cg(2,ipwdat)=cfft(2,i1,i2,i3dat)*xnorm
         end do
       end do
     end if

   else if (istwf_k>=2) then

     npwmin=1
     if (istwf_k==2 .and. me_g0==1) then
       ! Extract cg from cfft, in a way that projects on a 
       ! wavefunction with time-reversal symmetry
       do idat=1,ndat
         ipwdat = 1 + (idat-1) * npw
         i3dat = 1 + (idat-1)*n6
         cg(1,ipwdat)=cfft(1,1,1,i3dat)*xnorm
         cg(2,ipwdat)=0.d0
       end do
       npwmin=2
     end if

     if (use_symmetry==0) then

       if (iflag==-1) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipwdat,i3dat,i3invdat) IF (ndat>1)
         do idat=1,ndat
           do ipw=npwmin,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ! Construct the coordinates of -k-G
             i1inv=i1inver(i1); i2inv=i2inver(i2); i3inv=i3inver(i3)

             ipwdat = ipw + (idat-1) * npw
             i3dat = i3 + (idat-1) * n6
             i3invdat = i3inv + (idat-1) * n6

             ! Here the time-reversal symmetry is used to project from cfft
             cg(1,ipwdat)=(cfft(1,i1,i2,i3dat) + cfft(1,i1inv,i2inv,i3invdat))*0.5d0*xnorm
             cg(2,ipwdat)=(cfft(2,i1,i2,i3dat) - cfft(2,i1inv,i2inv,i3invdat))*0.5d0*xnorm
           end do
         end do

       else
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipwdat,i2dat,i2invdat) IF (ndat>1)
         do idat=1,ndat
           do ipw=npwmin,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ! Construct the coordinates of -k-G
             i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

             ipwdat = ipw + (idat-1) * npw
             i2dat = i2 + (idat-1) * n6
             i2invdat = i2inv + (idat-1) * n6

             ! Here the time-reversal symmetry is used to project from cfft
             cg(1,ipwdat)=(cfft(1,i1,i3,i2dat) + cfft(1,i1inv,i3inv,i2invdat))*0.5d0*xnorm
             cg(2,ipwdat)=(cfft(2,i1,i3,i2dat) - cfft(2,i1inv,i3inv,i2invdat))*0.5d0*xnorm
           end do
         end do
       end if

     else ! Use symmetry
       id1=n1/2+2
       id2=n2/2+2
       id3=n3/2+2

!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,j1,j2,j3,l1,l2,l3,ipwdat,i3dat,i3invdat) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw

           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           l1=kg_k(1,ipw)+shiftg(1)
           l2=kg_k(2,ipw)+shiftg(2)
           l3=kg_k(3,ipw)+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1 ; i1=j1+1
           if(j2<0)j2=j2+n2 ; i2=j2+1
           if(j3<0)j3=j3+n3 ; i3=j3+1

           ! Construct the coordinates of -k-G
           l1=i1inv-(i1inv/id1)*n1-1+shiftg(1)
           l2=i2inv-(i2inv/id2)*n2-1+shiftg(2)
           l3=i3inv-(i3inv/id3)*n3-1+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1 ; i1inv=j1+1
           if(j2<0)j2=j2+n2 ; i2inv=j2+1
           if(j3<0)j3=j3+n3 ; i3inv=j3+1

           ipwdat = ipw + (idat-1) * npw
           i3dat = i3 + (idat-1) * n6
           i3invdat = i3inv + (idat-1) * n6

           ! Here the time-reversal symmetry is used to project from cfft
           cg(1,ipwdat)=(cfft(1,i1,i2,i3dat) + cfft(1,i1inv,i2inv,i3invdat))*0.5d0*xnorm
           cg(2,ipwdat)=(cfft(2,i1,i2,i3dat) - cfft(2,i1inv,i2inv,i3invdat))*0.5d0*xnorm
         end do
       end do
     end if

   end if

 else
   write(msg,'(a,i0,a)')'  iflag=',iflag,' not acceptable.'
 end if

end subroutine sphere
!!***

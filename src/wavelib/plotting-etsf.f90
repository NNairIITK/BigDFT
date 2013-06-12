!> @file
!!  Routines to plot in ETSF format
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>   Write a field in the ISF basis in the ETSF format
subroutine write_etsf_density(filename,message,at,rxyz,n1i,n2i,n3i,hxh,hyh,hzh,&
     x, nspin)
  use module_base
  use module_types

  use etsf_io_low_level
  use etsf_io

  implicit none
  character(len=*), intent(in) :: filename,message
  integer, intent(in) :: n1i,n2i,n3i,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(wp), dimension(n1i,n2i,n3i,nspin), target, intent(in) :: x
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  !local variables
  integer :: nl1,nl2,nl3,nbx,nby,nbz,iat, ncid, i3, i2, i_stat, i_all,nc1,nc2,nc3
  double precision, dimension(3, 3), target :: rprim
  double precision, dimension(:,:), allocatable, target :: xred
  double precision, dimension(:), allocatable, target :: znucl
  double precision, dimension(:,:), allocatable :: buffer
  character(len=etsf_chemlen), allocatable, dimension(:), target :: spnames
  logical :: lstat
  character(len = etsf_io_low_error_len) :: error_string
  character(len = *), parameter :: subname = "write_etsf_density"
  
  type(etsf_dims) :: dims
  type(etsf_geometry), target :: geo
  type(etsf_io_low_error) :: error

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%astruct%geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
     nc1=n1i
     nc3=n3i
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
     nc1=n1i-31
     nc3=n3i-31
  end if
  !value of the buffer in the y direction
  if (at%astruct%geocode == 'P') then
     nl2=1
     nby = 1
     nc2=n2i
  else
     nl2=15
     nby = 0
     nc2=n2i-31
  end if

  call etsf_io_low_open_create(ncid, trim(filename) // ".etsf.nc", 3.3, lstat, &
       & error_data = error, title = 'Case for '//trim(message), overwrite = .true.)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if

  ! Unused dims, to be removed later
  dims%max_number_of_angular_momenta = etsf_no_dimension
  dims%max_number_of_basis_grid_points = etsf_no_dimension
  dims%max_number_of_coefficients = etsf_no_dimension
  dims%max_number_of_projectors = etsf_no_dimension
  dims%max_number_of_states = etsf_no_dimension
  dims%number_of_coefficients_dielectric_function = etsf_no_dimension
  dims%number_of_frequencies_dielectric_function = etsf_no_dimension
  dims%number_of_kpoints = etsf_no_dimension
  dims%number_of_localization_regions = etsf_no_dimension
  dims%number_of_qpoints_dielectric_function = etsf_no_dimension
  dims%number_of_qpoints_gamma_limit = etsf_no_dimension
  dims%number_of_spinor_components = etsf_no_dimension
  dims%number_of_spins = etsf_no_dimension
  dims%number_of_symmetry_operations = etsf_no_dimension
  dims%real_or_complex_coefficients = etsf_no_dimension
  dims%real_or_complex_gw_corrections = etsf_no_dimension
  dims%real_or_complex_wavefunctions = etsf_no_dimension

  ! Specific dims of interest.
  dims%number_of_atom_species        = at%astruct%ntypes
  dims%number_of_atoms               = at%astruct%nat
  dims%number_of_grid_points_vector1 = nc1!2*(n1+nbx)
  dims%number_of_grid_points_vector2 = nc2!2*(n2+nby)
  dims%number_of_grid_points_vector3 = nc3!2*(n3+nbz)
  if (nspin == 1) then
     dims%number_of_components       = 1
  else
     dims%number_of_components       = 4
  end if

  ! We write the dimensions to the file.
  call etsf_io_dims_def(ncid, dims, lstat, error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if  

  ! We set up the geometry.
  call etsf_io_geometry_def(ncid, lstat, error, &
       & flags = etsf_geometry_primitive_vectors + etsf_geometry_atom_species + &
       & etsf_geometry_red_at_pos + etsf_geometry_atomic_numbers + etsf_geometry_chemical_symbols)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
  ! We set up the density.
  call etsf_io_main_def(ncid, lstat, error, &
       & flags = etsf_main_density)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if

  ! We fill up the geometry
  rprim = reshape((/ (hxh * nc1),0.0_gp,0.0_gp, &
       & 0.0_gp,(hyh *nc2),0.0_gp, &
       & 0.0_gp,0.0_gp,(hzh *nc3) /), (/ 3, 3 /))
  allocate(xred(3, at%astruct%nat),stat=i_stat)
  call memocc(i_stat,xred,'xred',subname)
  do iat = 1, at%astruct%nat, 1
     xred(:, iat) = rxyz(:, iat) / (/ hxh * nc1, hyh * nc2, hzh * nc3 /)
  end do
  allocate(znucl(at%astruct%ntypes),stat=i_stat)
  call memocc(i_stat,znucl,'znucl',subname)
  znucl = real(at%nzatom)
  allocate(spnames(at%astruct%ntypes),stat=i_stat)
  call memocc(i_stat,spnames,'spnames',subname)
  do iat = 1, at%astruct%ntypes, 1
     call nzsymbol(at%nzatom(iat), spnames(iat))
  end do
  geo%chemical_symbols       => spnames
  geo%primitive_vectors      => rprim
  geo%atom_species           => at%astruct%iatype
  geo%atomic_numbers         => znucl
  geo%reduced_atom_positions => xred
  call etsf_io_geometry_put(ncid, geo, lstat, error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
  i_all=-product(shape(xred))*kind(xred)
  deallocate(xred)
  call memocc(i_stat,i_all,'xred',subname)
  i_all=-product(shape(znucl))*kind(znucl)
  deallocate(znucl)
  call memocc(i_stat,i_all,'znucl',subname)
  i_all=-product(shape(spnames))*kind(spnames)
  deallocate(spnames)
  call memocc(i_stat,i_all,'spnames',subname)

  ! We switch to write mode.
  call etsf_io_low_set_write_mode(ncid, lstat, error_data = error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
  
  ! We fill up the density.
  allocate(buffer(nc1, dims%number_of_components),stat=i_stat)
  call memocc(i_stat,buffer,'buffer',subname)
  do i3=0,nc3 - 1
     do i2=0,nc2 - 1
        buffer(:, 1) = x(nl1:nl1+nc1-1,i2+nl2,i3+nl3, 1)
        if (dims%number_of_components > 1) then
           buffer(:, 2) = x(nl1:nl1+nc1-1,i2+nl2,i3+nl3, 2)
           buffer(:, 3) = buffer(:, 1) + buffer(:, 2)
           buffer(:, 4) = buffer(:, 1) - buffer(:, 2)
        end if
        call etsf_io_low_write_var(ncid, "density", buffer, &
             & lstat, error_data = error, start = (/ 1, 1, i2 + 1, i3 + 1, 1 /), &
             & count = (/ 0, nc1, 1, 1,0 /))
        if (.not. lstat) then
           call etsf_io_low_error_to_str(error_string, error)
           write(0, "(A)") trim(error_string)
           stop
        end if
     end do
  end do
  i_all=-product(shape(buffer))*kind(buffer)
  deallocate(buffer)
  call memocc(i_stat,i_all,'buffer',subname)

  call etsf_io_low_close(ncid, lstat, error_data = error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
END SUBROUTINE write_etsf_density


!>   Read a field in the ISF basis in the ETSF format
subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz, iatypes, znucl)
  use module_base
  use module_types
  use etsf_io
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(out) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:,:,:), pointer :: rho
  real(gp), dimension(:,:), pointer :: rxyz
  integer, intent(out) ::  nat
  integer, dimension(:), pointer :: iatypes, znucl
  !local variables
  character(len=*), parameter :: subname='read_etsf'
  integer :: groupIds, i_stat, i_all, ncid
  integer :: n1t,n2t,n3t,n1,n2,n3, i2, i3, iat
  integer :: nl1,nl2,nl3,nbx,nby,nbz
  double precision, dimension(3, 3), target :: rprim
  double precision, dimension(:), allocatable, target :: znucl_
  logical :: lstat
  character(len = etsf_io_low_error_len) :: error_string

  type(etsf_io_low_error) :: error
  type(etsf_dims) :: dims
  type(etsf_split) :: split
  type(etsf_groups_flags) :: varIds
  type(etsf_groups) :: groups
  type(etsf_geometry), target :: geometry

  call etsf_io_data_contents(trim(filename) // ".etsf.nc", dims, split, &
       & groupIds, varIds, lstat, error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if

  groups%geometry => geometry
  geometry%primitive_vectors => rprim
  nat = dims%number_of_atoms
  allocate(rxyz(3, nat),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)
  geometry%reduced_atom_positions => rxyz
  allocate(iatypes(nat),stat=i_stat)
  call memocc(i_stat,iatypes,'iatypes',subname)
  geometry%atom_species => iatypes
  allocate(znucl_(dims%number_of_atom_species),stat=i_stat)
  call memocc(i_stat,znucl_,'znucl_',subname)
  geometry%atomic_numbers => znucl_
  call etsf_io_data_read(trim(filename) // ".etsf.nc", groups, lstat, error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
  do iat = 1, nat, 1
     rxyz(1, iat) = rxyz(1, iat) * rprim(1, 1)
     rxyz(2, iat) = rxyz(2, iat) * rprim(2, 2)
     rxyz(3, iat) = rxyz(3, iat) * rprim(3, 3)
  end do
  allocate(znucl(dims%number_of_atom_species),stat=i_stat)
  call memocc(i_stat,znucl,'znucl',subname)
  znucl = int(znucl_)
  i_all=-product(shape(znucl_))*kind(znucl_)
  deallocate(znucl_,stat=i_stat)
  call memocc(i_stat,i_all,'znucl_',subname)

  if (dims%number_of_components == 4) then
     nspin = 2
  else
     nspin = 1
  end if

  ! Now read the density 'by hand'.
  if (geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if
  n1t = dims%number_of_grid_points_vector1
  n2t = dims%number_of_grid_points_vector2
  n3t = dims%number_of_grid_points_vector3
  hxh = rprim(1,1) / n1t
  hyh = rprim(2,2) / n2t
  hzh = rprim(3,3) / n3t
  !grid positions
  n1=n1t/2-nbx
  n1i=2*n1+(1-nbx)+2*nl1
  n2=n2t/2-nby
  n2i=2*n2+(1-nby)+2*nl2
  n3=n3t/2-nbz
  n3i=2*n3+(1-nbz)+2*nl3
  if(associated(rho)) then
     i_all=-product(shape(rho))*kind(rho)
     deallocate(rho,stat=i_stat)
     call memocc(i_stat,i_all,'rho',subname)
  end if
  allocate(rho(n1i,n2i,n3i,nspin+ndebug) ,stat=i_stat)
  call memocc(i_stat,rho,'rho',subname)
  
  call etsf_io_low_open_read(ncid, trim(filename) // ".etsf.nc", lstat, error_data = error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if

  do i3=0,2*(n3+nbz) - 1
     do i2=0,2*(n2+nby) - 1
        !i1 = nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
        call etsf_io_low_read_var(ncid, "density", rho(nl1:nl1 + 2*(n1+nbx) - 1,i2+1,i3+1,1:nspin), & !rho(i1:i1 + 2*(n1+nbx) - 1, 1:nspin), &
             & lstat, error_data = error, start = (/ 1, 1, i2 + 1, i3 + 1, 1 /), &
             & count = (/ 1, 2*(n1+nbx), 1, 1, nspin /))
        if (.not. lstat) then
           call etsf_io_low_error_to_str(error_string, error)
           write(0, "(A)") trim(error_string)
           stop
        end if
     end do
  end do
  
  call etsf_io_low_close(ncid, lstat, error_data = error)
  if (.not. lstat) then
     call etsf_io_low_error_to_str(error_string, error)
     write(0, "(A)") trim(error_string)
     stop
  end if
END SUBROUTINE read_etsf

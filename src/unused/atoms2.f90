!> Calculate the symmetries and update
subroutine atoms_set_symmetries(atoms, rxyz, disableSym, tol, elecfield)
  use module_base
  use module_types
  use abi_defs_basis
  use m_ab6_symmetry
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  logical, intent(in) :: disableSym
  real(gp), intent(in) :: tol
  real(gp), intent(in) :: elecfield(3)
  !local variables
  character(len=*), parameter :: subname='atoms_set_symmetries'
  integer :: i_stat, ierr, i_all
  real(gp) :: rprimd(3, 3)
  real(gp), dimension(:,:), allocatable :: xRed

  ! Calculate the symmetries, if needed
  if (atoms%astruct%geocode /= 'F') then
     if (atoms%astruct%sym%symObj < 0) then
        call symmetry_new(atoms%astruct%sym%symObj)
     end if
     ! Adjust tolerance
     if (tol > 0._gp) call symmetry_set_tolerance(atoms%astruct%sym%symObj, tol, ierr)
     ! New values
     rprimd(:,:) = 0
     rprimd(1,1) = atoms%astruct%cell_dim(1)
     rprimd(2,2) = atoms%astruct%cell_dim(2)
     if (atoms%astruct%geocode == 'S') rprimd(2,2) = 1000._gp
     rprimd(3,3) = atoms%astruct%cell_dim(3)
     call symmetry_set_lattice(atoms%astruct%sym%symObj, rprimd, ierr)
     allocate(xRed(3, atoms%astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,xRed,'xRed',subname)
     xRed(1,:) = modulo(rxyz(1, :) / rprimd(1,1), 1._gp)
     xRed(2,:) = modulo(rxyz(2, :) / rprimd(2,2), 1._gp)
     xRed(3,:) = modulo(rxyz(3, :) / rprimd(3,3), 1._gp)
     call symmetry_set_structure(atoms%astruct%sym%symObj, atoms%astruct%nat, atoms%astruct%iatype, xRed, ierr)
     i_all=-product(shape(xRed))*kind(xRed)
     deallocate(xRed,stat=i_stat)
     call memocc(i_stat,i_all,'xRed',subname)
     if (atoms%astruct%geocode == 'S') then
        call symmetry_set_periodicity(atoms%astruct%sym%symObj, &
             & (/ .true., .false., .true. /), ierr)
     else if (atoms%astruct%geocode == 'F') then
        call symmetry_set_periodicity(atoms%astruct%sym%symObj, &
             & (/ .false., .false., .false. /), ierr)
     end if
     !if (all(in%elecfield(:) /= 0)) then
     !     ! I'm not sure what this subroutine does!
     !   call symmetry_set_field(atoms%astruct%sym%symObj, (/ in%elecfield(1) , in%elecfield(2),in%elecfield(3) /), ierr)
     !elseif (in%elecfield(2) /= 0) then
     !   call symmetry_set_field(atoms%astruct%sym%symObj, (/ 0._gp, in%elecfield(2), 0._gp /), ierr)
     if (elecfield(2) /= 0) then
        call symmetry_set_field(atoms%astruct%sym%symObj, (/ 0._gp, elecfield(2), 0._gp /), ierr)
     end if
     if (disableSym) then
        call symmetry_set_n_sym(atoms%astruct%sym%symObj, 1, &
             & reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), (/ 3 ,3, 1 /)), &
             & reshape((/ 0.d0, 0.d0, 0.d0 /), (/ 3, 1/)), (/ 1 /), ierr)
     end if
  else
     call deallocate_symmetry(atoms%astruct%sym, subname)
     atoms%astruct%sym%symObj = -1
  end if
END SUBROUTINE atoms_set_symmetries

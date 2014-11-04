subroutine fnrmandforcemax_old(ff,fnrm,fmax,at)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), intent(in):: ff(3,at%astruct%nat)
  real(gp), intent(out):: fnrm, fmax
  real(gp):: t1,t2,t3
  integer:: iat

  fmax=0._gp
  do iat=1,at%astruct%nat
     call frozen_alpha(at%astruct%ifrztyp(iat),1,ff(1,iat)**2,t1)
     call frozen_alpha(at%astruct%ifrztyp(iat),2,ff(2,iat)**2,t2)
     call frozen_alpha(at%astruct%ifrztyp(iat),3,ff(3,iat)**2,t3)
     fmax=max(fmax,sqrt(t1+t2+t3))
  enddo

  !This is the norm of the forces of non-blocked atoms
  call atomic_dot(at,ff,ff,fnrm)
END SUBROUTINE fnrmandforcemax_old

!> Should we evaluate the translational force also with blocked atoms?
subroutine transforce_forfluct(at,fxyz,sumx,sumy,sumz)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp),intent(in):: fxyz(3,at%astruct%nat)
  real(gp), intent(out) :: sumx,sumy,sumz
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  !atomic_dot with one
  sumx=0._gp 
  sumy=0._gp 
  sumz=0._gp
  do iat=1,at%astruct%nat

     call frozen_alpha(at%astruct%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(at%astruct%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(at%astruct%ifrztyp(iat),3,1.0_gp,alphaz)

     sumx=sumx+alphax*fxyz(1,iat) 
     sumy=sumy+alphay*fxyz(2,iat) 
     sumz=sumz+alphaz*fxyz(3,iat)

  end do
END SUBROUTINE transforce_forfluct


!> Calculate the coefficient for moving atoms following the ifrztyp
subroutine frozen_alpha(ifrztyp,ixyz,alpha,alphai)
  use module_base
  implicit none
  integer, intent(in) :: ifrztyp !< Frozen code of the iat atom 
  integer, intent(in) :: ixyz    !< Direction (1=x,y=2,z=3)
  real(gp), intent(in) :: alpha
  real(gp), intent(out) :: alphai

  if (move_this_coordinate(ifrztyp,ixyz)) then
     alphai=alpha
  else
     alphai=0.0_gp
  end if
 
END SUBROUTINE frozen_alpha

!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: rxyz=txyz+alpha*sxyz
!!   all the shift are inserted into the box if there are periodic directions
!!   if the atom are frozen they are not moved
subroutine atomic_axpy(astruct,txyz,alpha,sxyz,rxyz)
  implicit none
  real(gp), intent(in) :: alpha
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,astruct%nat), intent(inout) :: rxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  do iat=1,astruct%nat
     !adjust the moving of the atoms following the frozen direction
     call frozen_alpha(astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(astruct%ifrztyp(iat),3,alpha,alphaz)

     if (astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),astruct%cell_dim(1))
        rxyz(2,iat)=modulo(txyz(2,iat)+alphay*sxyz(2,iat),astruct%cell_dim(2))
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),astruct%cell_dim(3))
     else if (astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),astruct%cell_dim(1))
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),astruct%cell_dim(3))
     else
        rxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
     end if
  end do

END SUBROUTINE atomic_axpy


!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: fxyz=txyz+alpha*sxyz
!!   update the forces taking into account the frozen atoms
!!   do not apply the modulo operation on forces 
subroutine atomic_axpy_forces(astruct,txyz,alpha,sxyz,fxyz)
  implicit none
  real(gp), intent(in) :: alpha
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,astruct%nat), intent(inout) :: fxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz
  
  do iat=1,astruct%nat
     !adjust the moving of the forces following the frozen direction
     call frozen_alpha(astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(astruct%ifrztyp(iat),3,alpha,alphaz)

     fxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
     fxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
     fxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
  end do
  
END SUBROUTINE atomic_axpy_forces


!>z=alpha*A*x + beta* y
subroutine atomic_gemv(astruct,m,alpha,A,x,beta,y,z)
  implicit none
  integer, intent(in) :: m
  real(gp), intent(in) :: alpha,beta
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: x
  real(gp), dimension(m), intent(in) :: y
  real(gp), dimension(m,3,astruct%nat), intent(in) :: A
  real(gp), dimension(m), intent(out) :: z
  !local variables
  integer :: iat,i,j
  real(gp) :: mv,alphai
  
  do i=1,m
     mv=0.0_gp
     do iat=1,astruct%nat
        do j=1,3
           call frozen_alpha(astruct%ifrztyp(iat),j,A(i,j,iat),alphai)
           mv=mv+alphai*x(j,iat)
        end do
     end do
     z(i)=alpha*mv+beta*y(i)
  end do

END SUBROUTINE atomic_gemv


!> rxyz=txyz+alpha*sxyz
subroutine atomic_coordinate_axpy(astruct,ixyz,iat,t,alphas,r)
  implicit none
  integer, intent(in) :: ixyz,iat
  real(gp), intent(in) :: t,alphas
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(out) :: r
  !local variables
  logical :: periodize
  real(gp) :: alat,alphai

  if (ixyz == 1) then
     alat=astruct%cell_dim(1)
  else if (ixyz == 2) then
     alat=astruct%cell_dim(2)
  else if (ixyz == 3) then
     alat=astruct%cell_dim(3)
  else
     alat = -1
     write(0,*) "Internal error"
     stop
  end if
  
  periodize= astruct%geocode == 'P' .or. &
       (astruct%geocode == 'S' .and. ixyz /= 2)

  call frozen_alpha(astruct%ifrztyp(iat),ixyz,alphas,alphai)

  if (periodize) then
     r=modulo(t+alphai,alat)
  else
     r=t+alphai
  end if

END SUBROUTINE atomic_coordinate_axpy



!> Calculate the symmetries and update
subroutine atoms_set_symmetries(atoms, rxyz, disableSym, tol, elecfield)
  use module_base
  use module_types
  use defs_basis
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

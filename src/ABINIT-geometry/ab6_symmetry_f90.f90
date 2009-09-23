!* * Fortran90 source file *
!*
!* Copyright (c) 2008-2009 ABINIT Group (Damien Caliste)
!* All rights reserved.
!*
!* This file is part of the ABINIT software package. For license information,
!* please see the COPYING file in the top-level directory of the ABINIT source
!* distribution.
!*
!*

module ab6_symmetry

  use defs_basis

  implicit none

  private

  integer, parameter, public :: AB6_MAX_SYMMETRIES = 384

  type, private :: symmetry
     ! The input characteristics
     real(dp) :: tolsym
     real(dp) :: rprimd(3,3), gprimd(3,3), rmet(3,3)
     integer :: nAtoms
     integer, pointer :: typeAt(:)
     real(dp), pointer :: xRed(:,:)

     logical :: withField
     real(dp) :: field(3)

     logical :: withJellium

     integer :: withSpin
     real(dp), pointer :: spinAt(:,:)

     logical :: withSpinOrbit

     ! The output characteristics
     ! The bravais parameters
     integer :: nBravSym
     integer :: bravais(11), bravSym(3, 3, AB6_MAX_SYMMETRIES)
     ! The symmetry matrices
     integer  :: nSym
     integer  :: sym(3, 3, AB6_MAX_SYMMETRIES)
     real(dp) :: transNon(3, AB6_MAX_SYMMETRIES)
     integer  :: symafm(AB6_MAX_SYMMETRIES)
     ! Some additional information
     integer          :: multiplicity
     real(dp)         :: genAfm(3)
     character(len=5) :: pointGroup
     integer          :: spaceGroup, pointGroupMagn
     integer, pointer :: indexingAtoms(:,:,:)
  end type symmetry

  ! We store here a list of symmetry objects to be able to
  ! call several symmetry operations on different objects.
  ! The simplest portable way to do it, is to create
  ! a list of Fortran structure and to use the list index
  ! as an identifier that can be given to the other languages.
  type, private :: symmetry_list
     integer                       :: id
     type(symmetry_list),  pointer :: next
     type(symmetry)                :: data
  end type symmetry_list
  type(symmetry_list), pointer :: my_symmetries
  integer :: n_symmetries = 0

  logical, private, parameter :: AB_DBG = .false.

  public :: ab6_symmetry_new
  public :: ab6_symmetry_free
  public :: ab6_symmetry_set_lattice
  public :: ab6_symmetry_set_structure
  public :: ab6_symmetry_set_spin
  public :: ab6_symmetry_set_spin_orbit
  public :: ab6_symmetry_set_field
  public :: ab6_symmetry_set_jellium

  public :: ab6_symmetry_get_n_atoms
  public :: ab6_symmetry_get_multiplicity
  public :: ab6_symmetry_get_bravais
  public :: ab6_symmetry_get_matrices
  public :: ab6_symmetry_get_group
  public :: ab6_symmetry_get_equivalent_atom
  public :: ab6_symmetry_get_mp_k_grid
  public :: ab6_symmetry_get_auto_k_grid

contains

  subroutine new_item(token)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    type(symmetry_list), pointer :: token

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(0,*) "AB symmetry: create a new token."

    ! Init case, very first call.
    if (n_symmetries == 0) then
       nullify(my_symmetries)
    end if

    ! Normal treatment.
    n_symmetries = n_symmetries + 1

    allocate(token)
    token%id = n_symmetries
    call new_symmetry(token%data)
    token%next => my_symmetries

    my_symmetries => token
    if (AB_DBG) write(0,*) "AB symmetry: creation OK with id ", token%id
  end subroutine new_item

  subroutine free_item(token)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    type(symmetry_list), pointer :: token

    type(symmetry_list), pointer :: tmp

    if (.not. associated(token)) then
       return
    end if

    call free_symmetry(token%data)

    if (AB_DBG) write(0,*) "AB symmetry: free request on token ", token%id
    ! We remove token from the list.
    if (my_symmetries%id == token%id) then
       my_symmetries => token%next
    else
       tmp => my_symmetries
       do
          if (.not.associated(tmp)) then
             return
          end if
          if (associated(tmp%next) .and. tmp%next%id == token%id) then
             exit
          end if
          tmp => tmp%next
       end do
       tmp%next => token%next
    end if
    deallocate(token)
    if (AB_DBG) write(0,*) "AB symmetry: free done"
  end subroutine free_item

  subroutine get_item(token, id)


    type(symmetry_list), pointer :: token
    integer, intent(in) :: id

    type(symmetry_list), pointer :: tmp
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: request list element ", id
    nullify(token)

    tmp => my_symmetries
    do
       if (.not. associated(tmp)) then
          exit
       end if
       if (tmp%id == id) then
          token => tmp
          return
       end if
       tmp => tmp%next
    end do
  end subroutine get_item

  subroutine new_symmetry(sym)


    type(symmetry), intent(out) :: sym

    if (AB_DBG) write(0,*) "AB symmetry: create a new symmetry object."
    nullify(sym%xRed)
    nullify(sym%spinAt)
    nullify(sym%typeAt)
    sym%tolsym   = tol8
    sym%nSym     = -1
    sym%nBravSym = -1
    sym%withField   = .false.
    sym%withJellium = .false.
    sym%withSpin = 0
    sym%withSpinOrbit = .false.
    sym%multiplicity = -1
    nullify(sym%indexingAtoms)
    write(sym%pointGroup, "(A)") "-----"
  end subroutine new_symmetry

  subroutine free_symmetry(sym)


    type(symmetry), intent(inout) :: sym

    if (AB_DBG) write(0,*) "AB symmetry: free a symmetry."

    if (associated(sym%xRed)) deallocate(sym%xRed)
    if (associated(sym%spinAt)) deallocate(sym%spinAt)
    if (associated(sym%typeAt)) deallocate(sym%typeAt)
    if (associated(sym%indexingAtoms)) deallocate(sym%indexingAtoms)
  end subroutine free_symmetry





  subroutine ab6_symmetry_new(id)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(out) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call new symmetry."
    call new_item(token)
    id = token%id
  end subroutine ab6_symmetry_new

  subroutine ab6_symmetry_free(id)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call free symmetry."

    call get_item(token, id)
    if (associated(token)) call free_item(token)
  end subroutine ab6_symmetry_free

  subroutine ab6_symmetry_set_lattice(id, rprimd, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_42_geometry
!End of the abilint section

    integer, intent(in) :: id
    real(dp), intent(in) :: rprimd(3,3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    real(dp) :: ucvol
    real(dp) :: gmet(3,3)

    if (AB_DBG) write(0,*) "AB symmetry: call set lattice."
    if (AB_DBG) write(0, "(A,3F12.6,A)") "  (", rprimd(:,1), ")"
    if (AB_DBG) write(0, "(A,3F12.6,A)") "  (", rprimd(:,2), ")"
    if (AB_DBG) write(0, "(A,3F12.6,A)") "  (", rprimd(:,3), ")"

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    token%data%rprimd = rprimd
    call metric(gmet, token%data%gprimd, -1, token%data%rmet, rprimd, ucvol)

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_lattice

  subroutine ab6_symmetry_set_structure(id, nAtoms, typeAt, xRed, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    integer, intent(in) :: typeAt(nAtoms)
    real(dp), intent(in) :: xRed(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: call set structure."
    if (AB_DBG) write(0, "(A,I3,A)") "  ", nAtoms, " atoms"
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          write(0, "(A,3F12.6,I3)") "  ", xRed(:, i), typeAt(i)
       end do
    end if

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    token%data%nAtoms =  nAtoms
    allocate(token%data%typeAt(nAtoms))
    token%data%typeAt = typeAt
    allocate(token%data%xRed(3, nAtoms))
    token%data%xRed   = xRed

    ! We unset only the symmetries
    token%data%nSym     = -1
    if (associated(token%data%indexingAtoms)) deallocate(token%data%indexingAtoms)
  end subroutine ab6_symmetry_set_structure

  subroutine ab6_symmetry_set_spin(id, nAtoms, spinAt, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    real(dp), intent(in) :: spinAt(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: call set spin."
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          write(0, "(A,3F12.6)") "  ", spinAt(:, i)
       end do
    end if

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if
    if (token%data%nAtoms /= nAtoms) then
       errno = AB6_ERROR_ARG
       return
    end if

    token%data%withSpin = 4
    allocate(token%data%spinAt(3, nAtoms))
    token%data%spinAt = spinAt

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_spin

  subroutine ab6_symmetry_set_collinear_spin(id, nAtoms, spinAt, errno)
    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    integer, intent(in) :: spinAt(nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: call set collinear spin."
    if (AB_DBG) then
       do i = 1, nAtoms, 1
          write(0, "(A,I3)") "  ", spinAt(i)
       end do
    end if

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if
    if (token%data%nAtoms /= nAtoms) then
       errno = AB6_ERROR_ARG
       return
    end if

    token%data%withSpin = 2
    allocate(token%data%spinAt(1, nAtoms))
    token%data%spinAt = real(reshape(spinAt, (/ 1, nAtoms /)), dp)

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_collinear_spin

  subroutine ab6_symmetry_set_spin_orbit(id, withSpinOrbit, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    logical, intent(in) :: withSpinOrbit
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: call set spin orbit."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    token%data%withSpinOrbit = withSpinOrbit

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_spin_orbit

  subroutine ab6_symmetry_set_field(id, field, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    real(dp), intent(in) :: field(3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set field."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    token%data%withField = .true.
    token%data%field = field

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_field

  subroutine ab6_symmetry_set_jellium(id, jellium, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    logical, intent(in) :: jellium
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set jellium."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    token%data%withJellium = jellium

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_jellium





  subroutine ab6_symmetry_get_n_atoms(id, nAtoms, errno)
    !scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nAtoms

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get nAtoms."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    nAtoms = token%data%nAtoms
  end subroutine ab6_symmetry_get_n_atoms

  subroutine compute_bravais(sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_42_geometry
!End of the abilint section

    type(symmetry), intent(inout) :: sym

    integer :: berryopt

    ! We do the computation
    if (sym%withField) then
       berryopt = 4
    else
       berryopt = 0
    end if
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT symbrav."
    call symbrav(berryopt, sym%bravais, AB6_MAX_SYMMETRIES, &
         & sym%nBravSym, sym%bravSym, sym%rmet, sym%rprimd, sym%tolsym)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."
    if (AB_DBG) write(0, "(A,I3)") "  nSymBrav :", sym%nBravSym
    if (AB_DBG) write(0, "(A,I3)") "  holohedry:", sym%bravais(1)
    if (AB_DBG) write(0, "(A,I3)") "  center   :", sym%bravais(2)
  end subroutine compute_bravais

  subroutine ab6_symmetry_get_bravais(id, bravais, holohedry, center, &
       & nBravSym, bravSym, errno)
    !scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nBravSym, holohedry, center
    !arrays
    integer, intent(out) :: bravais(3,3), bravSym(3, 3, AB6_MAX_SYMMETRIES)

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get bravais."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%nBravSym < 0) then
       ! We do the computation
       call compute_bravais(token%data)
    end if
    
    holohedry = token%data%bravais(1)
    center    = token%data%bravais(2)
    bravais   = reshape(token%data%bravais(3:11), (/ 3,3 /))
    nBravSym  = token%data%nBravSym
    bravSym(3, 3, 1:nBravSym) = token%data%bravSym(3, 3, 1:nBravSym)
  end subroutine ab6_symmetry_get_bravais

  subroutine compute_matrices(sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_42_geometry
!End of the abilint section

    type(symmetry), intent(inout) :: sym

    integer :: berryopt, jellslab, noncol, shubnikov, isym, problem
    integer :: nsym_nomagn, isym_nomagn, nspinor, nspden, pawspnorb
    integer, allocatable :: sym_nomagn(:,:,:)
    integer :: identity(3,3)
    real(dp), allocatable :: transNon_nomagn(:,:)
    real(dp), pointer :: spinAt_(:,:)
    character(len=5) :: ptgroupha

    if (sym%nBravSym < 0) then
       ! We do the computation of the Bravais part.
       call compute_bravais(sym)
    end if

    if (sym%withField) then
       berryopt = 4
    else
       berryopt = 0
    end if
    if (sym%withJellium) then
       jellslab = 1
    else
       jellslab = 0
    end if
    if (sym%withSpin == 4) then
       nspden = 4
       noncol = 1
       spinAt_ => sym%spinAt
    else if (sym%withSpin == 2) then
       nspden = 2
       noncol = 0
       spinAt_ => sym%spinAt
    else
       nspden = 1
       noncol = 0
       allocate(spinAt_(3, sym%nAtoms))
       spinAt_ = 0
    end if
    if (sym%withSpinOrbit) then
       nspinor = 2
       pawspnorb = 1
    else
       nspinor = 1
       pawspnorb = 0
    end if

    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT symfind."
    call symfind(berryopt, sym%field, sym%gprimd, jellslab, AB6_MAX_SYMMETRIES, &
         & sym%nAtoms, noncol, sym%nBravSym, nspden, nspinor, sym%nSym, &
         & pawspnorb, sym%bravSym, spinAt_, &
         & sym%symafm, sym%sym, sym%transNon, sym%tolsym, sym%typeAt, sym%xRed)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."
    if (AB_DBG) write(0, "(A,I3)") "  nSym:", sym%nSym

    if (sym%withSpin == 0) then
       deallocate(spinAt_)
    end if

    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT chkprimit."
    call chkprimit(0, sym%multiplicity, sym%nSym, sym%symAfm, sym%sym)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."
    if (AB_DBG) write(0, "(A,I3)") "  multi:", sym%multiplicity

    if (sym%multiplicity == 1) then
       ! The cell is primitive, so that the space group can be
       ! determined. Need to distinguish Fedorov and Shubnikov groups.
       ! Do not distinguish Shubnikov types I and II.
       ! Also identify genafm, in case of Shubnikov type IV
       identity(:,:) = reshape((/1,0,0,0,1,0,0,0,1/), (/3,3/))
       shubnikov = 1
       do isym = 1, sym%nSym, 1
          if(sym%symAfm(isym) == -1)then
             shubnikov = 3
             if(sum(abs(sym%sym(:,:,isym) - identity(:,:))) == 0)then
                shubnikov = 4
                sym%genAfm(:) = sym%transNon(:,isym)
                exit
             end if
          end if
       end do

       if (AB_DBG) write(0, "(A,I3)") "  shubni:", shubnikov
       if(shubnikov == 1 .or. shubnikov == 3)then
          !  Find the point group
          call symanal(sym%bravais, sym%nSym, problem, sym%pointGroup, sym%sym)
          if (AB_DBG) write(0, "(A,I3)") "  problem:", problem
          !  Find the space group
          call symspgr(sym%bravais, sym%nSym, sym%spaceGroup, sym%sym, sym%transNon)
       end if

       if(shubnikov /= 1)then

          !  Determine nonmagnetic symmetry operations
          nsym_nomagn = sym%nSym / 2
          allocate(sym_nomagn(3, 3, nsym_nomagn))
          allocate(transNon_nomagn(3, nsym_nomagn))

          isym_nomagn = 0
          do isym = 1, sym%nSym
             if(sym%symAfm(isym) == 1)then
                isym_nomagn = isym_nomagn + 1
                sym_nomagn(:,:,isym_nomagn)    = sym%sym(:,:,isym)
                transNon_nomagn(:,isym_nomagn) = sym%transNon(:,isym)
             end if
          end do

          if(shubnikov==3)then
             !   Find the point group of the halved symmetry set
             call symanal(sym%bravais, nsym_nomagn, problem, ptgroupha, sym_nomagn)

             !   Deduce the magnetic point group from ptgroup and ptgroupha
             call getptgroupma(sym%pointGroup, ptgroupha, sym%pointGroupMagn)
          else if(shubnikov==4)then
             !   Find the Fedorov space group of the halved symmetry set
             call symspgr(sym%bravais, nsym_nomagn, sym%spaceGroup, &
                  & sym_nomagn, transNon_nomagn)

             !   The magnetic translation generator has already been determined
          end if

          deallocate(sym_nomagn)
          deallocate(transNon_nomagn)
       end if ! Shubnikov groups
    end if

  end subroutine compute_matrices

  subroutine ab6_symmetry_get_matrices(id, nSym, sym, transNon, symAfm, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nSym
    integer, intent(out)  :: sym(3, 3, AB6_MAX_SYMMETRIES)
    integer, intent(out)  :: symAfm(AB6_MAX_SYMMETRIES)
    real(dp), intent(out) :: transNon(3, AB6_MAX_SYMMETRIES)

    type(symmetry_list), pointer :: token
    integer :: berryopt, jellslab

    if (AB_DBG) write(0,*) "AB symmetry: call get matrices."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%nSym < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data)
    end if

    nSym                = token%data%nSym
    sym(:, :, 1:nSym)   = token%data%sym(:, :, 1:nSym)
    symAfm(1:nSym)      = token%data%symAfm(1:nSym)
    transNon(:, 1:nSym) = token%data%transNon(:, 1:nSym)
  end subroutine ab6_symmetry_get_matrices

  subroutine ab6_symmetry_get_multiplicity(id, multiplicity, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: id
    integer, intent(out) :: multiplicity, errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get multiplicity."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%multiplicity < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data)
    end if
    multiplicity = token%data%multiplicity
  end subroutine ab6_symmetry_get_multiplicity

  subroutine ab6_symmetry_get_group(id, pointGroup, spaceGroup, &
       & pointGroupMagn, genAfm, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in)           :: id
    integer, intent(out)          :: errno
    real(dp), intent(out)         :: genAfm(3)
    character(len=5), intent(out) :: pointGroup
    integer, intent(out)          :: spaceGroup, pointGroupMagn

    type(symmetry_list), pointer  :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get group."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%multiplicity < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data)
    end if

    if (token%data%multiplicity /= 1) then
       errno = AB6_ERROR_SYM_NOT_PRIMITIVE
       return
    end if

    write(pointGroup, "(A5)") token%data%pointGroup
    pointGroupMagn = token%data%pointGroupMagn
    spaceGroup     = token%data%spaceGroup
    genAfm         = token%data%genAfm
  end subroutine ab6_symmetry_get_group

  subroutine compute_equivalent_atoms(sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_32_util
! use interfaces_42_geometry
!End of the abilint section

    type(symmetry), intent(inout) :: sym

    integer, allocatable :: symrec(:,:,:)
    integer :: isym

    if (.not. associated(sym%indexingAtoms)) &
         & allocate(sym%indexingAtoms(4, sym%nSym, sym%nAtoms))

    !Get the symmetry matrices in terms of reciprocal basis
    allocate(symrec(3, 3, sym%nSym))
    do isym = 1, sym%nSym, 1
       call mati3inv(sym%sym(:,:,isym), symrec(:,:,isym))
    end do
    
    !Obtain a list of rotated atom labels:
    call symatm(sym%indexingAtoms, sym%nAtoms, sym%nSym, symrec, &
         & sym%transNon, sym%tolsym, sym%typeAt, sym%xRed)

    deallocate(symrec)
  end subroutine compute_equivalent_atoms

  subroutine ab6_symmetry_get_equivalent_atom(id, equiv, iAtom, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in)  :: id
    integer, intent(in)  :: iAtom
    integer, intent(out) :: equiv(4, AB6_MAX_SYMMETRIES)
    integer, intent(out) :: errno

    type(symmetry_list), pointer  :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get equivalent."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (iAtom < 1 .or. iAtom > token%data%nAtoms) then
       errno = AB6_ERROR_ARG
       return
    end if

    if (.not. associated(token%data%indexingAtoms)) then
       ! We do the computation of the matrix part.
       call compute_equivalent_atoms(token%data)
    end if

    equiv(:, 1:token%data%nSym) = token%data%indexingAtoms(:,:,iAtom)
  end subroutine ab6_symmetry_get_equivalent_atom

  subroutine ab6_symmetry_get_mp_k_grid(id, nkpt, kpt, wkpt, &
       & ngkpt, nshiftk, shiftk, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_56_recipspace
!End of the abilint section

    integer, intent(in)  :: id
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(in) :: nshiftk
    real(dp), intent(in) :: shiftk(3, nshiftk)
    integer, intent(out) :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)

    type(symmetry_list), pointer  :: token
    integer :: nshiftk_, nkpt_
    real(dp) :: kptrlen
    integer :: kptrlatt(3,3)
    real(dp) :: shiftk_(3, 8)

    if (AB_DBG) write(0,*) "AB symmetry: call get k grid."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%nSym < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data)
    end if

    ! First, compute the number of kpoints
    kptrlatt(:,:) = 0
    kptrlatt(1,1) = ngkpt(1)
    kptrlatt(2,2) = ngkpt(2)
    kptrlatt(3,3) = ngkpt(3)
    kptrlen = 20.
    nshiftk_ = nshiftk
    shiftk_(:, 1:nshiftk) = shiftk

    call getkgrid((/ 1, 1, 1 /), 6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, 0, nkpt, nshiftk_, token%data%nSym, &
         & token%data%rprimd, shiftk_, token%data%symAfm, token%data%sym, &
         & token%data%transNon, (/ 0, 0, 0 /), wkpt)

    ! Then, we call it again to get the actual values for the k points.
    nshiftk_ = nshiftk
    shiftk_(:, 1:nshiftk) = shiftk
    allocate(kpt(3, nkpt))
    allocate(wkpt(nkpt))
    call getkgrid((/ 1, 1, 1 /), 6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk_, token%data%nSym, &
         & token%data%rprimd, shiftk_, token%data%symAfm, token%data%sym, &
         & token%data%transNon, (/ 0, 0, 0 /), wkpt)
  end subroutine ab6_symmetry_get_mp_k_grid

  subroutine ab6_symmetry_get_auto_k_grid(id, nkpt, kpt, wkpt, &
       & kptrlen, errno)
    integer, intent(in)  :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nkpt
    real(dp), intent(in) :: kptrlen
    real(dp), pointer :: kpt(:,:), wkpt(:)

    type(symmetry_list), pointer  :: token
    integer :: nshiftk_, nkpt_
    integer :: kptrlatt(3,3)
    real(dp) :: shiftk_(3, 8)

    if (AB_DBG) write(0,*) "AB symmetry: call get MP k grid."

    errno = AB6_NO_ERROR
    call get_item(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_OBJ
       return
    end if

    if (token%data%nSym < 0) then
       ! We do the computation of the matrix part.
       call compute_matrices(token%data)
    end if

    !  The parameters of the k lattice are not known, compute
    !  kptrlatt, nshiftk, shiftk.
    call testkgrid(token%data%bravais,6,kptrlatt,kptrlen,&
         & AB6_MAX_SYMMETRIES,nshiftk_,token%data%nSym,0,token%data%rprimd,&
         & shiftk_,token%data%symAfm,token%data%sym,token%data%transNon,(/ 0, 0, 0 /))
    if (AB_DBG) write(0,*) "AB symmetry: testkgrid -> kptrlatt=", kptrlatt

    call getkgrid((/ 1, 1, 1 /), 6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, 0, nkpt, nshiftk_, token%data%nSym, &
         & token%data%rprimd, shiftk_, token%data%symAfm, token%data%sym, &
         & token%data%transNon, (/ 0, 0, 0 /), wkpt)
    if (AB_DBG) write(0,*) "AB symmetry: getkgrid -> nkpt=", nkpt
    allocate(kpt(3, nkpt))
    allocate(wkpt(nkpt))

    ! Then, we call it again to get the actual values for the k points.
    call getkgrid((/ 1, 1, 1 /), 6, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk_, token%data%nSym, &
         & token%data%rprimd, shiftk_, token%data%symAfm, token%data%sym, &
         & token%data%transNon, (/ 0, 0, 0 /), wkpt)
  end subroutine ab6_symmetry_get_auto_k_grid

end module ab6_symmetry

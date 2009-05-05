module ab6_symmetry

  use defs_basis
  use defs_datatypes
  use interfaces_12geometry

  implicit none

  private

  integer, parameter, public :: AB6_MAX_SYMMETRIES = 384

  type, private :: symmetry
     ! The input characteristics
     real(dp) :: rprimd(3,3), gprimd(3,3), rmet(3,3)
     integer :: nAtoms
     integer, pointer :: typeAt(:)
     real(dp), pointer :: xRed(:,:)

     logical :: withField
     real(dp) :: field(3)

     logical :: withJellium

     logical :: withNonColinearSpin
     real(dp), pointer :: spinAt(:,:)

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
     type(symmetry_list),  pointer :: next => null()
     type(symmetry)                :: data
  end type symmetry_list
  type(symmetry_list), pointer :: my_symmetries => null()
  integer :: n_symmetries = 0

  logical, private, parameter :: AB_DBG = .false.

  ! Error codes
  integer, parameter, public :: AB6_NO_ERROR                = 0
  integer, parameter, public :: AB6_ERROR_SYM_OBJ           = 1
  integer, parameter, public :: AB6_ERROR_SYM_ARG           = 2
  integer, parameter, public :: AB6_ERROR_SYM_NOT_PRIMITIVE = 3

  public :: ab6_symmetry_new
  public :: ab6_symmetry_free
  public :: ab6_symmetry_set_lattice
  public :: ab6_symmetry_set_structure
  public :: ab6_symmetry_set_spin
  public :: ab6_symmetry_set_field
  public :: ab6_symmetry_set_jellium

  public :: ab6_symmetry_get_multiplicity
  public :: ab6_symmetry_get_bravais
  public :: ab6_symmetry_get_matrices
  public :: ab6_symmetry_get_group
  public :: ab6_symmetry_get_equivalent_atom
  public :: ab6_symmetry_get_k_grid

contains

  subroutine new_token(token)
    type(symmetry_list), pointer :: token

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(0,*) "AB symmetry: create a new token."
    n_symmetries = n_symmetries + 1

    allocate(token)
    token%id = n_symmetries
    call new_symmetry(token%data)
    token%next => my_symmetries

    my_symmetries => token
    if (AB_DBG) write(0,*) "AB symmetry: creation OK with id ", token%id
  end subroutine new_token

  subroutine free_token(token)
    type(symmetry_list), pointer :: token

    integer :: idtset
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
  end subroutine free_token

  subroutine get_token(token, id)
    type(symmetry_list), pointer :: token
    integer, intent(in) :: id

    type(symmetry_list), pointer :: tmp
    integer :: i

    if (AB_DBG) write(0,*) "AB symmetry: request list element ", id
    nullify(token)
    ! List element are prepended so element id is at (nb - id) position.
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
  end subroutine get_token

  subroutine new_symmetry(sym)
    type(symmetry), intent(out) :: sym

    if (AB_DBG) write(0,*) "AB symmetry: create a new symmetry object."
    nullify(sym%xRed)
    nullify(sym%spinAt)
    nullify(sym%typeAt)
    sym%nSym     = -1
    sym%nBravSym = -1
    sym%withField   = .false.
    sym%withJellium = .false.
    sym%withNonColinearSpin = .false.
    sym%multiplicity = -1
    nullify(sym%indexingAtoms)
  end subroutine new_symmetry

  subroutine free_symmetry(sym)
    type(symmetry), intent(inout) :: sym

    if (AB_DBG) write(0,*) "AB symmetry: free a symmetry."

    if (associated(sym%indexingAtoms)) deallocate(sym%indexingAtoms)
  end subroutine free_symmetry





  subroutine ab6_symmetry_new(id)
    integer, intent(out) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call new symmetry."
    call new_token(token)
    id = token%id
  end subroutine ab6_symmetry_new

  subroutine ab6_symmetry_free(id)
    integer, intent(in) :: id

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call free symmetry."

    call get_token(token, id)
    if (associated(token)) call free_token(token)
  end subroutine ab6_symmetry_free

  subroutine ab6_symmetry_set_lattice(id, rprimd, errno)
    integer, intent(in) :: id
    real(dp), intent(in) :: rprimd(3,3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token
    real(dp) :: ucvol
    real(dp) :: gmet(3,3)

    if (AB_DBG) write(0,*) "AB symmetry: call set lattice."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    token%data%rprimd = rprimd
    call metric(gmet, token%data%gprimd, -1, token%data%rmet, rprimd, ucvol)

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_lattice

  subroutine ab6_symmetry_set_structure(id, nAtoms, typeAt, xRed, errno)
    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    integer, intent(in), target :: typeAt(nAtoms)
    real(dp), intent(in), target :: xRed(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set structure."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    token%data%nAtoms =  nAtoms
    token%data%typeAt => typeAt
    token%data%xRed   => xRed

    ! We unset only the symmetries
    token%data%nSym     = -1
    if (associated(token%data%indexingAtoms)) deallocate(token%data%indexingAtoms)
  end subroutine ab6_symmetry_set_structure

  subroutine ab6_symmetry_set_spin(id, nAtoms, spinAt, errno)
    integer, intent(in) :: id
    integer, intent(in) :: nAtoms
    real(dp), intent(in), target :: spinAt(3,nAtoms)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set spin."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token) .or. token%data%nAtoms /= nAtoms) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    token%data%withNonColinearSpin = .true.
    token%data%spinAt => spinAt

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_spin

  subroutine ab6_symmetry_set_field(id, field, errno)
    integer, intent(in) :: id
    real(dp), intent(in) :: field(3)
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set field."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    token%data%withField = .true.
    token%data%field = field

    ! We unset all the computed symmetries
    token%data%nBravSym = -1
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_field

  subroutine ab6_symmetry_set_jellium(id, errno)
    integer, intent(in) :: id
    integer, intent(out) :: errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call set jellium."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    token%data%withJellium = .true.

    ! We unset only the symmetries
    token%data%nSym     = -1
  end subroutine ab6_symmetry_set_jellium





  subroutine compute_bravais(sym)
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
         & sym%nBravSym, sym%bravSym, sym%rmet, sym%rprimd)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."
  end subroutine compute_bravais

  subroutine ab6_symmetry_get_bravais(id, bravais, holohedry, center, &
       & nBravSym, bravSym, errno)
    !scalars
    integer, intent(in) :: id
    integer, intent(out) :: errno
    integer, intent(out) :: nBravSym, holohedry, center
    !arrays
    integer, intent(out) :: bravais(3,3), bravSym(3, 3, AB6_MAX_SYMMETRIES)

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get bravais."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
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
    type(symmetry), intent(inout) :: sym

    integer :: berryopt, jellslab, noncol, shubnikov, isym, problem
    integer :: nsym_nomagn, isym_nomagn
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
    if (sym%withNonColinearSpin) then
       noncol = 1
       spinAt_ => sym%spinAt
    else
       noncol = 0
       allocate(spinAt_(3, sym%nAtoms))
       spinAt_ = 0
    end if

    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT symfind."
    call symfind(berryopt, sym%field, sym%gprimd, jellslab, AB6_MAX_SYMMETRIES, &
         & sym%nAtoms, noncol, sym%nBravSym, sym%nSym, sym%bravSym, spinAt_, &
         & sym%symafm, sym%sym, sym%transNon, sym%typeAt, sym%xRed)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."

    if (.not. sym%withNonColinearSpin) then
       deallocate(spinAt_)
    end if

    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT chkprimit."
    call chkprimit(0, sym%multiplicity, sym%nSym, sym%symAfm, sym%sym)
    if (AB_DBG) write(0,*) "AB symmetry: call ABINIT OK."

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

       if(shubnikov == 1 .or. shubnikov == 3)then
          !  Find the point group
          call symanal(sym%bravais, sym%nSym, problem, sym%pointGroup, sym%sym)
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
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
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
    integer, intent(in) :: id
    integer, intent(out) :: multiplicity, errno

    type(symmetry_list), pointer :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get multiplicity."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
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
    integer, intent(in)           :: id
    integer, intent(out)          :: errno
    real(dp), intent(out)         :: genAfm(3)
    character(len=5), intent(out) :: pointGroup
    integer, intent(out)          :: spaceGroup, pointGroupMagn

    type(symmetry_list), pointer  :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get group."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
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
         & sym%transNon, sym%typeAt, sym%xRed)

    deallocate(symrec)
  end subroutine compute_equivalent_atoms

  subroutine ab6_symmetry_get_equivalent_atom(id, equiv, iAtom, errno)
    integer, intent(in)  :: id
    integer, intent(in)  :: iAtom
    integer, intent(out) :: equiv(4, AB6_MAX_SYMMETRIES)
    integer, intent(out) :: errno

    type(symmetry_list), pointer  :: token

    if (AB_DBG) write(0,*) "AB symmetry: call get equivalent."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
       return
    end if

    if (iAtom < 1 .or. iAtom > token%data%nAtoms) then
       errno = AB6_ERROR_SYM_ARG
       return
    end if

    if (.not. associated(token%data%indexingAtoms)) then
       ! We do the computation of the matrix part.
       call compute_equivalent_atoms(token%data)
    end if

    equiv(:, 1:token%data%nSym) = token%data%indexingAtoms(:,:,iAtom)
  end subroutine ab6_symmetry_get_equivalent_atom

  subroutine ab6_symmetry_get_k_grid(id, nkpt, kpt, wkpt, &
       & ngkpt, nshiftk, shiftk, errno)
    integer, intent(in)  :: id
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(in) :: nshiftk
    real(dp), intent(in) :: shiftk(3, nshiftk)
    integer, intent(out) :: nkpt
    real(dp), intent(out) :: kpt(3, ngkpt(1) * ngkpt(2) * ngkpt(3))
    real(dp), intent(out) :: wkpt(ngkpt(1) * ngkpt(2) * ngkpt(3))

    type(symmetry_list), pointer  :: token
    integer :: nshiftk_, nkpt_
    real(dp) :: kptrlen
    integer :: kptrlatt(3,3)
    real(dp) :: shiftk_(3, 8)
    real(dp), allocatable :: kpt_(:,:), wkpt_(:)

    if (AB_DBG) write(0,*) "AB symmetry: call get k grid."

    errno = AB6_NO_ERROR
    call get_token(token, id)
    if (.not. associated(token)) then
       errno = AB6_ERROR_SYM_OBJ
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
    allocate(kpt_(3, nkpt))
    allocate(wkpt_(nkpt))
    call getkgrid((/ 1, 1, 1 /), 6, 1, kpt_, 1, kptrlatt, kptrlen, &
         & AB6_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk_, token%data%nSym, &
         & token%data%rprimd, shiftk_, token%data%symAfm, token%data%sym, &
         & token%data%transNon, (/ 0, 0, 0 /), wkpt_)
    kpt(:, 1:nkpt) = kpt_
    wkpt(1:nkpt) = wkpt_
  end subroutine ab6_symmetry_get_k_grid

end module ab6_symmetry

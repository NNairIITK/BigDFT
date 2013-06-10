!> @file
!!  Module to store all dictionary keys of the input files.
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Define all static strings to store input variables
module module_input_keys
  use dictionaries, only: dictionary

  implicit none
  
  private

  character(len = *), parameter, public :: DFT_VARIABLES = "dft"
  character(len = *), parameter, public :: HGRIDS = "hgrids"
  character(len = *), parameter, public :: RMULT = "rmult"
  character(len = *), parameter, public :: IXC = "ixc"
  character(len = *), parameter, public :: NCHARGE = "ncharge"
  character(len = *), parameter, public :: ELECFIELD = "elecfield"
  character(len = *), parameter, public :: NSPIN = "nspin", MPOL = "mpol"
  character(len = *), parameter, public :: GNRM_CV = "gnrm_cv"
  character(len = *), parameter, public :: ITERMAX = "itermax", NREPMAX = "nrepmax"
  character(len = *), parameter, public :: NCONG = "ncong", IDSX = "idsx"
  character(len = *), parameter, public :: DISPERSION = "dispersion"
  character(len = *), parameter, public :: INPUTPSIID = "inputpsiid"
  character(len = *), parameter, public :: OUTPUT_WF = "output_wf"
  character(len = *), parameter, public :: OUTPUT_DENSPOT = "output_denspot"
  character(len = *), parameter, public :: RBUF = "rbuf"
  character(len = *), parameter, public :: NCONGT = "ncongt"
  character(len = *), parameter, public :: NORBV = "norbv", NVIRT = "nvirt"
  character(len = *), parameter, public :: NPLOT = "nplot"
  character(len = *), parameter, public :: DISABLE_SYM = "disablesym"

  !> Error ids for this module.
  integer, parameter, public :: INPUT_VAR_NOT_IN_LIST      = 1
  integer, parameter, public :: INPUT_VAR_NOT_IN_RANGE     = 2
  integer, parameter, public :: INPUT_VAR_CONVERSION_ERROR = 3

  character(len = *), parameter :: RANGE = "__range__", EXCLUSIVE = "__exclusive__"
  character(len = *), parameter :: DEFAULT = "default", COMMENT = "__comment__"

  public :: input_keys_init, input_keys_finalize
  public :: input_keys_set, input_keys_fill

  type(dictionary), pointer :: parameters

contains
  subroutine input_keys_init()
    use yaml_output
    use dictionaries
    implicit none

    call dict_init(parameters)
    call set(parameters // DFT_VARIABLES, get_dft_parameters())

    call yaml_dict_dump(parameters, comment_key = COMMENT)
  END SUBROUTINE input_keys_init
  
  subroutine input_keys_finalize()
    use dictionaries
    implicit none

    call dict_free(parameters)
  END SUBROUTINE input_keys_finalize
  
  function get_dft_parameters()
    use dictionaries
    implicit none
    type(dictionary), pointer :: p, get_dft_parameters

    call dict_init(p)

    ! Settings
    call set(p // HGRIDS, dict_new( (/ &
         & COMMENT   .is. "grid spacing in the three directions (bohr)", &
         & RANGE     .is. list_new((/ .item."0.", .item."2." /)), &
         & DEFAULT   .is. list_new((/ .item."0.45", .item."0.45", .item."0.45" /)), &
         & "fast"    .is. list_new((/ .item."0.55", .item."0.55", .item."0.55" /)), &
         & "accurate".is. list_new((/ .item."0.30", .item."0.30", .item."0.30" /)) /) ))

    call set(p // RMULT, dict_new( (/ &
         & COMMENT   .is. "c(f)rmult*radii_cf(:,1(2))=coarse(fine) atom-basec radius", &
         & RANGE     .is. list_new((/ .item."0.", .item."100." /)), &
         & DEFAULT   .is. list_new((/ .item."5.", .item."8." /)) /)))

    call set(p // IXC, dict_new((/ &
         & COMMENT   .is. "exchange-correlation parameter (LDA=1,PBE=11)", &
         & DEFAULT   .is. "1", &
         & "PBE"     .is. "11", &
         & "Hybrid"  .is. "-406" /)))

    call set(p // NCHARGE, dict_new((/ &
         & COMMENT   .is. "charge of the system", &
         & RANGE     .is. list_new((/ .item."-500.", .item."500." /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // ELECFIELD, dict_new((/ &
         & COMMENT   .is. "electric fielc (Ex,Ey,Ez)", &
         & DEFAULT   .is. list_new((/ .item."0.", .item."0.", .item."0." /)) /)))

    call set(p // NSPIN, dict_new((/ &
         & COMMENT   .is. "spin polarization", &
         & EXCLUSIVE .is. list_new((/ .item."1", .item."2", .item."4" /)), &
         & DEFAULT   .is. "1" /)))

    call set(p // MPOL, dict_new((/ &
         & COMMENT .is. "total magnetic moment", &
         & DEFAULT  . is. "0" /)))

    call set(p // GNRM_CV, dict_new((/ &
         & COMMENT   .is. "convergence criterion gradient", &
         & RANGE     .is. list_new((/ .item."1.e-20", .item."1." /)), &
         & DEFAULT   .is. "1e-4", &
         & "fast"    .is. "1e-3", &
         & "accurate".is. "1e-5"/)))

    call set(p // ITERMAX, dict_new((/ &
         & COMMENT   .is. "max. # of wfn. opt. steps", &
         & RANGE     .is. list_new((/ .item."0", .item."10000" /)), &
         & DEFAULT   .is. "50" /)))

    call set(p // NREPMAX, dict_new((/ &
         & COMMENT   .is. "max. # of re-diag. runs", &
         & RANGE     .is. list_new((/ .item."0", .item."1000" /)), &
         & DEFAULT   .is. "1", &
         & "accurate".is. "10" /)))

    call set(p // NCONG, dict_new((/ &
         & COMMENT   .is. "# of CG it. for preconditioning eq.", &
         & RANGE     .is. list_new((/ .item."0", .item."20" /)), &
         & DEFAULT   .is. "6" /)))

    call set(p // IDSX, dict_new((/ &
         & COMMENT     .is. "wfn. diis history", &
         & RANGE       .is. list_new((/ .item."0", .item."15" /)), &
         & DEFAULT     .is. "6", &
         & "low memory".is. "2" /)))

    call set(p // DISPERSION, dict_new((/ &
         & COMMENT   .is. "dispersion correction potential (values 1,2,3,4,5), 0=none", &
         & RANGE     .is. list_new((/ .item."0", .item."5" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // INPUTPSIID, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. list_new((/ .item."0", .item."1", .item."2", &
         &                            .item."10", .item."11", .item."12", &
         &                            .item."20", .item."21", .item."22" /)) /)))

    call set(p // OUTPUT_WF, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. list_new((/ .item."-2", .item."-1", .item."0", .item."2", &
         &                            .item."10", .item."12", .item."13", .item."100", &
         &                            .item."101", .item."102" /)) /)))

    call set(p // OUTPUT_DENSPOT, dict_new((/ &
         & DEFAULT   .is. "0", &
         & EXCLUSIVE .is. list_new((/ .item."0", .item."1", .item."2", .item."3" /)) /)))

    call set(p // RBUF, dict_new((/ &
         & COMMENT     .is. "length of the tail (AU)", &
         & RANGE       .is. list_new((/ .item."0.", .item."10." /)), &
         & DEFAULT     .is. "0.", &
         & "with tails".is. "5." /)))

    call set(p // NCONGT, dict_new((/ &
         & COMMENT   .is. "# tail CG iterations", &
         & RANGE     .is. list_new((/ .item."0", .item."50" /)), &
         & DEFAULT   .is. "30" /)))

    call set(p // NORBV, dict_new((/ &
         & COMMENT   .is. "Davidson subspace dim.", &
         & RANGE     .is. list_new((/ .item."-9999", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // NVIRT, dict_new((/ &
         & COMMENT   .is. "# of opt. orbs", &
         & RANGE     .is. list_new((/ .item."0", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // NPLOT, dict_new((/ &
         & COMMENT   .is. "# of plottec orbs", &
         & RANGE     .is. list_new((/ .item."0", .item."9999" /)), &
         & DEFAULT   .is. "0" /)))

    call set(p // DISABLE_SYM, dict_new((/ &
         & COMMENT   .is. "disable the symmetry detection", &
         & DEFAULT   .is. "No" /)))

    get_dft_parameters => p
  END FUNCTION get_dft_parameters

  subroutine pouet()
    implicit none
  end subroutine pouet

  !> Read a real or real/real, real:real 
  !! Here the fraction is indicated by the ':' or '/'
  !! The problem is that / is a separator for Fortran
  subroutine read_fraction_string(string,var,ierror)
    use module_base
    implicit none
    !Arguments
    character(len=*), intent(in) :: string
    double precision, intent(out) :: var
    integer, intent(out) :: ierror
    !Local variables
    character(len=max_field_length) :: tmp
    integer :: num,den,pfr,psp

    !First look at the first blank after trim
    tmp(1:max_field_length)=trim(string)
    psp = scan(tmp,' ')

    !see whether there is a fraction in the string
    if(psp==0) psp=len(tmp)
    pfr = scan(tmp(1:psp),':')
    if (pfr == 0) pfr = scan(tmp(1:psp),'/')
    !It is not a fraction
    if (pfr == 0) then
       read(tmp(1:psp),*,iostat=ierror) var
    else 
       read(tmp(1:pfr-1),*,iostat=ierror) num
       read(tmp(pfr+1:psp),*,iostat=ierror) den
       if (ierror == 0) var=real(num,gp)/real(den,gp)
    end if
    !Value by defaut
    if (ierror /= 0) var = huge(1_gp) 
  END SUBROUTINE read_fraction_string

  subroutine input_keys_fill(dict, file, profile)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    character(len = *), intent(in), optional :: profile

    integer :: i
    type(dictionary), pointer :: ref
    character(max_field_length), dimension(:), allocatable :: keys

    ref = parameters // file
    allocate(keys(dict_len(ref)))
    keys = dict_keys(ref)

    if (present(profile)) then
       do i = 1, size(keys), 1
          call input_keys_set(dict // keys(i), file, keys(i), profile)
       end do
    else
       do i = 1, size(keys), 1
          call input_keys_set(dict // keys(i), file, keys(i))
       end do
    end if

!!$  ! Validate inputPsiId value.
!!$  if (f_err_check(INPUT_VAR_NOT_IN_LIST) .and. iproc == 0) then
!!$     call input_psi_help()
!!$     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierror)
!!$  end if

!!$  ! Validate output_wf value.
!!$  if (f_err_check(INPUT_VAR_NOT_IN_LIST) .and. iproc == 0) then
!!$     call output_wf_format_help()
!!$     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierror)
!!$  end if

    deallocate(keys)
  END SUBROUTINE input_keys_fill

  subroutine input_keys_set(dict, file, key, profile)
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file, key
    character(len = *), intent(in), optional :: profile

    integer :: i
    type(dictionary), pointer :: ref, list
    character(len = max_field_length) :: val, item, profile_
    double precision, dimension(2) :: rg
    logical :: found

    profile_(1:max_field_length) = " "
    if (present(profile)) profile_(1:max_field_length) = profile
    if (trim(profile_) == "") profile_(1:max_field_length) = DEFAULT
    
    ! Disable some dict errors to allow default values to be set.
    call f_err_set_callback(pouet)
    val = dict
    call f_err_unset_callback()

    ref => parameters // file // key
    if (.not.f_err_check(DICT_VALUE_ABSENT)) then
       ! There is already a value in dict.
       if (has_key(ref, val)) then
          ! The key was asking for a profile, we copy it.
          call copy(dict, ref // val)
       else
          ! The value is user-defined, we validate it with info from dict from.
          if (has_key(ref, RANGE)) then
             rg(1) = ref // RANGE // 0
             rg(2) = ref // RANGE // 1
             call validate(dict, rg)
          else if (has_key(ref, EXCLUSIVE)) then
             call yaml_dict_dump(ref)
             list => ref // EXCLUSIVE
             found = .false.
             do i = 0, dict_len(list) - 1, 1
                item = list // i
                found = (trim(val) == trim(item))
                if (found) exit
             end do
             if (f_err_raise(.not. found, err_id = INPUT_VAR_NOT_IN_LIST, &
                  & err_msg = trim(dict%data%key) // " = '" // trim(val) // "' is not allowed.")) return
          end if
       end if
    else
       ! There is no value in dict, we take it from dict_from.
       call copy(dict, ref // profile_)
       call f_err_clean()
    end if

  contains

    recursive subroutine validate(dict, rg)
      implicit none
      type(dictionary), pointer :: dict
      double precision, dimension(2), intent(in) :: rg
      
      character(len = max_field_length) :: val
      double precision :: var
      integer :: ierror

      if (.not. associated(dict%child)) then
         val = dict
         call read_fraction_string(val, var, ierror)
         if (f_err_raise(ierror /= 0, err_id = INPUT_VAR_CONVERSION_ERROR, &
              & err_msg = "cannot read a double from '" // trim(val) //"'.")) return

         write(*,*) RANGE, rg, var
         if (f_err_raise(var < rg(1) .or. var > rg(2), err_id = INPUT_VAR_NOT_IN_RANGE, &
              & err_msg = trim(dict%data%key) // " = '" // trim(val) // &
              & "' not in range.")) return
         
         if (associated(dict%next)) call validate(dict%next, rg)
      else
         call validate(dict%child, rg)
      end if
    end subroutine validate

    subroutine copy(dict, ref)
      implicit none
      type(dictionary), pointer :: dict, ref

      character(len = max_field_length) :: val
      integer :: i

      if (dict_len(ref) > 1) then
         do i = 0, dict_len(ref) - 1, 1
            val = ref // i
            call set(dict // i, val)
         end do
      else
         val = ref
         call set(dict, val)
      end if
    end subroutine copy

  END SUBROUTINE input_keys_set

end module module_input_keys

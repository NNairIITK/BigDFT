!>  Modules which contains all interfaces to parse input dictionary.
module module_input_dicts

  implicit none

  private

  !parameters to avoid typos in dictionary keys
  character(len=*), parameter :: ATOMIC_OCC="Atomic occupation"
  character(len=*), parameter :: ASTRUCT_UNITS = 'units' 
  character(len=*), parameter :: ASTRUCT_CELL = 'cell' 
  character(len=*), parameter :: ASTRUCT_POSITIONS = 'positions' 
  character(len=*), parameter :: ASTRUCT_PROPERTIES = 'properties' 
  character(len=*), parameter :: GOUT_ENERGY = 'energy (Ha)' 
  character(len=*), parameter :: GOUT_FORCES = 'forces (Ha/Bohr)' 

  ! update a dictionary from a input file
  public :: merge_input_file_to_dict

  ! Main creation routine
  public :: user_dict_from_files

  ! Dictionary completion
  public :: psp_dict_fill_all, psp_dict_analyse

  ! Dictionary inquire
  public :: astruct_dict_get_source, astruct_dict_get_types

  ! Types from dictionaries
  public :: astruct_set_from_dict
  public :: psp_set_from_dict, nlcc_set_from_dict
  public :: atomic_data_set_from_dict
  public :: occupation_set_from_dict

  ! Types to dictionaries
  public :: psp_data_merge_to_dict
  public :: astruct_merge_to_dict
  public :: global_output_merge_to_dict

  ! Dictionaries from files (old formats).
  public :: psp_file_merge_to_dict, nlcc_file_merge_to_dict
  public :: atoms_file_merge_to_dict
  public :: astruct_file_merge_to_dict
  public :: atomic_data_file_merge_to_dict
  public :: occupation_data_file_merge_to_dict

contains

  !> logical function, true if file is existing (basically a fortran inquire)
  function file_exists(filename) result(exists)
    use module_base, only: f_err_raise
    use yaml_output, only: yaml_toa
    implicit none
    character(len=*), intent(in) :: filename
    logical :: exists
    !local variables
    integer :: ierr
    inquire(file = filename, exist = exists,iostat=ierr)
    if (f_err_raise(ierr /=0,'Error in unit inquiring for filename '//&
         trim(filename)//', ierr='//trim(yaml_toa(ierr)),&
         err_name='BIGDFT_RUNTIME_ERROR')) then
       exists=.false.
       return
    end if
  end function file_exists

  !> Routine to read YAML input files and create input dictionary.
  !! Update the input dictionary with the result of yaml_parse
  subroutine merge_input_file_to_dict(dict, fname, mpi_env)
    use module_base
    !use yaml_output, only :yaml_map
    use dictionaries
    use yaml_parse, only: yaml_parse_from_char_array
    implicit none
    type(dictionary), pointer :: dict !<dictionary of the input files. Should be initialized on entry
    character(len = *), intent(in) :: fname !<name of the file where the dictionaryt has to be read from 
    type(mpi_environment), intent(in) :: mpi_env !<environment of the reading. Used for broadcasting the result
    !local variables
    integer(kind = 8) :: cbuf, cbuf_len
    integer :: ierr
    character(len = max_field_length) :: val
    character, dimension(:), allocatable :: fbuf
    type(dictionary), pointer :: udict
    external :: getFileContent,copyCBuffer,freeCBuffer

    call f_routine(id='merge_input_file_to_dict')
    if (mpi_env%iproc == 0) then
       call getFileContent(cbuf, cbuf_len, fname, len_trim(fname))
       if (mpi_env%nproc > 1) &
            & call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    else
       call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    end if
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')

    if (mpi_env%iproc == 0) then
       call copyCBuffer(fbuf(1), cbuf, cbuf_len)
       call freeCBuffer(cbuf)
       if (mpi_env%nproc > 1) &
            & call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
    else
       call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
    end if

    call f_err_open_try()
    call yaml_parse_from_char_array(udict, fbuf)
    ! Handle with possible partial dictionary.
    call f_free_str(1,fbuf)
    call dict_update(dict, udict // 0)
    call dict_free(udict)
    ierr = 0
    if (f_err_check()) ierr = f_get_last_error(val)
    call f_err_close_try()
    !in the present implementation f_err_check is not cleaned after the close of the try
    if (ierr /= 0) call f_err_throw(err_id = ierr, err_msg = val)
    call f_release_routine()

  end subroutine merge_input_file_to_dict

  subroutine user_dict_from_files(dict,radical,posinp, mpi_env)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use module_defs, only: mpi_environment
    use module_interfaces, only: read_input_dict_from_files
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: radical, posinp
    type(mpi_environment), intent(in) :: mpi_env

    character(len = max_field_length) :: str

    !read the input file(s) and transform them into a dictionary
    call read_input_dict_from_files(trim(radical), mpi_env, dict)

    !consider to move the reading of the atomic position at first place
    if (.not. has_key(dict, "posinp")) then
       !Add old posinp formats
       call astruct_file_merge_to_dict(dict, "posinp", trim(posinp))
    else
       str = dict_value(dict // "posinp")
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          call astruct_file_merge_to_dict(dict, "posinp", trim(str))
       end if
    end if

    ! Add old psppar
    call atoms_file_merge_to_dict(dict)

    !when the user has not specified the occupation in the input file
    if (.not. has_key(dict, ATOMIC_OCC)) then
       ! Add old input.occup
       !call atomic_data_file_merge_to_dict(dict, ATOMIC_OCC, &
       !     & trim(radical) // ".occup")
       !yaml format should be used even for old method
       if (file_exists(trim(radical)//".occup")) &
            call merge_input_file_to_dict(dict//ATOMIC_OCC,trim(radical)//".occup",mpi_env)
    else !otherwise the input file always supersedes
       str = dict_value(dict // ATOMIC_OCC)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          !call atomic_data_file_merge_to_dict(dict, ATOMIC_OCC, trim(str))
          if (file_exists(trim(str))) &
               call merge_input_file_to_dict(dict//ATOMIC_OCC,trim(str),mpi_env)
       end if
    end if

    if (.not. has_key(dict, "occupation")) then
       ! Add old input.occ
       call occupation_data_file_merge_to_dict(dict, "occupation", &
            & trim(radical) // ".occ")
    else
       str = dict_value(dict // "occupation")
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          call occupation_data_file_merge_to_dict(dict, "occupation", trim(str))
       end if
    end if

  end subroutine user_dict_from_files

  subroutine psp_dict_fill_all(dict, atomname, run_ixc)
    use module_defs, only: gp, UNINITIALIZED, bigdft_mpi
    use ao_inguess, only: atomic_info
    use dictionaries
    use dynamic_memory
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: atomname
    integer, intent(in) :: run_ixc

    integer :: ixc, ierr
    character(len=27) :: filename
    logical :: exists
    integer :: nzatom, nelpsp, npspcode
    real(gp) :: psppar(0:4,0:6)
!    integer, parameter :: nmax=6,lmax=4
    !integer, parameter :: nelecmax=32
!    character(len=2) :: symbol
    integer :: i!,mxpl,mxchg,nsccode
    real(gp) :: ehomo,radfine,rad!,amu,rcov,rprb
!    real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
    type(dictionary), pointer :: radii
    real(gp), dimension(3) :: radii_cf
    character(len = max_field_length) :: source

    filename = 'psppar.' // atomname

    radii_cf = UNINITIALIZED(1._gp)
    if (has_key(dict // filename, "Radii of active regions (AU)")) then
       radii => dict // filename // "Radii of active regions (AU)"
       if (has_key(radii, "Coarse")) radii_cf(1) =  radii // "Coarse"
       if (has_key(radii, "Fine")) radii_cf(2) =  radii // "Fine"
       if (has_key(radii, "Coarse PSP")) radii_cf(3) =  radii // "Coarse PSP"
    end if

    exists = has_key(dict // filename, "Local Pseudo Potential (HGH convention)")
    if (.not. exists) then
       ixc = run_ixc
       if (has_key(dict // filename, "Pseudopotential XC")) &
            & ixc = dict // filename // "Pseudopotential XC"
       call psp_from_data(atomname, nzatom, &
            & nelpsp, npspcode, ixc, psppar(:,:), exists)
       call psp_data_merge_to_dict(dict // filename, nzatom, nelpsp, npspcode, ixc, &
            & psppar(0:4,0:6), radii_cf, UNINITIALIZED(1._gp), UNINITIALIZED(1._gp))
       call set(dict // filename // "Source", "Hard-coded")
    else
       nzatom = dict // filename // "Atomic number"
       nelpsp = dict // filename // "No. of Electrons"
    end if

    if (.not. exists) then
       call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
       write(*,'(1x,5a)')&
            'ERROR: The pseudopotential parameter file "',trim(filename),&
            '" is lacking, and no registered pseudo found for "', &
            & trim(atomname), '", exiting...'
       stop
    end if

    write(source, "(A)") "User-defined"
    if (radii_cf(1) == UNINITIALIZED(1.0_gp)) then
       !see whether the atom is semicore or not
       !and consider the ground state electronic configuration
       call atomic_info(nzatom,nelpsp,ehomo=ehomo)
       !call eleconf(nzatom, nelpsp,symbol,rcov,rprb,ehomo,&
       !     neleconf,nsccode,mxpl,mxchg,amu)

       !assigning the radii by calculating physical parameters
       radii_cf(1)=1._gp/sqrt(abs(2._gp*ehomo))
       write(source, "(A)") "Hard-coded"
    end if
    if (radii_cf(2) == UNINITIALIZED(1.0_gp)) then
       radfine = dict // filename // "Local Pseudo Potential (HGH convention)" // "Rloc"
       if (has_key(dict // filename, "NonLocal PSP Parameters")) then
          do i=1, dict_len(dict // filename // "NonLocal PSP Parameters")
             rad = dict // filename // "NonLocal PSP Parameters" // (i - 1) // "Rloc"
             if (rad /= 0._gp) then
                radfine=min(radfine, rad)
             end if
          end do
       end if
       radii_cf(2)=radfine
       write(source, "(A)") "Hard-coded"
    end if
    if (radii_cf(3) == UNINITIALIZED(1.0_gp)) then
       radii_cf(3)=radfine
       write(source, "(A)") "Hard-coded"
    end if
    radii => dict // filename // "Radii of active regions (AU)"
    call set(radii // "Coarse", radii_cf(1))
    call set(radii // "Fine", radii_cf(2))
    call set(radii // "Coarse PSP", radii_cf(3))
    call set(radii // "Source", source)
  end subroutine psp_dict_fill_all

  subroutine psp_dict_analyse(dict, atoms)
    use module_defs, only: gp
    use module_types, only: atoms_data
    use module_atoms, only: allocate_atoms_data
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    type(atoms_data), intent(inout) :: atoms

    integer :: ityp
    character(len = 27) :: filename
    real(gp), dimension(3) :: radii_cf
    logical :: pawpatch, l
    integer :: paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices

    if (.not. associated(atoms%nzatom)) then
       call allocate_atoms_data(atoms)
       !call allocate_atoms_nat(atoms, "psp_dict_analyse")
       !call allocate_atoms_ntypes(atoms, "psp_dict_analyse")
    end if

    pawpatch = .true.
    do ityp=1,atoms%astruct%ntypes
       filename = 'psppar.'//atoms%astruct%atomnames(ityp)
       call psp_set_from_dict(dict // filename, &
            & atoms%nzatom(ityp), atoms%nelpsp(ityp), atoms%npspcode(ityp), &
            & atoms%ixcpsp(ityp), atoms%psppar(:,:,ityp), radii_cf)
       !To eliminate the runtime warning due to the copy of the array (TD)
       atoms%radii_cf(ityp,:)=radii_cf(:)

       l = .false.
       if (has_key(dict // filename, "PAW patch")) l = dict // filename // "PAW patch"
       pawpatch = pawpatch .and. l
    end do
    call nlcc_set_from_dict(dict, atoms)

    if (pawpatch) then
       paw_tot_l=0
       paw_tot_q=0
       paw_tot_coefficients=0
       paw_tot_matrices=0
       do ityp=1,atoms%astruct%ntypes
          filename = 'psppar.'//atoms%astruct%atomnames(ityp)
          call pawpatch_from_file( filename, atoms,ityp,&
               paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .false.)
       end do
       do ityp=1,atoms%astruct%ntypes
          filename = 'psppar.'//atoms%astruct%atomnames(ityp)
          !! second time allocate and then store
          call pawpatch_from_file( filename, atoms,ityp,&
               paw_tot_l, paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .true.)
       end do
    else
       nullify(atoms%paw_l,atoms%paw_NofL,atoms%paw_nofchannels)
       nullify(atoms%paw_nofgaussians,atoms%paw_Greal,atoms%paw_Gimag)
       nullify(atoms%paw_Gcoeffs,atoms%paw_H_matrices,atoms%paw_S_matrices,atoms%paw_Sm1_matrices)
    end if
  end subroutine psp_dict_analyse

  subroutine nlcc_set_from_dict(dict, atoms)
    use module_defs, only: gp
    use module_types, only: atoms_data
    use dictionaries
    use dynamic_memory
    implicit none
    type(dictionary), pointer :: dict
    type(atoms_data), intent(inout) :: atoms

    type(dictionary), pointer :: nloc, coeffs
    integer :: ityp, nlcc_dim, n, i
    character(len=27) :: filename

    nlcc_dim = 0
    do ityp = 1, atoms%astruct%ntypes, 1
       atoms%nlcc_ngc(ityp)=0
       atoms%nlcc_ngv(ityp)=0
       filename = 'psppar.' // trim(atoms%astruct%atomnames(ityp))
       if (.not. has_key(dict, filename)) cycle    
       if (.not. has_key(dict // filename, 'Non Linear Core Correction term')) cycle
       nloc => dict // filename // 'Non Linear Core Correction term'
       if (has_key(nloc, "Valence") .or. has_key(nloc, "Conduction")) then
          n = 0
          if (has_key(nloc, "Valence")) n = dict_len(nloc // "Valence")
          nlcc_dim = nlcc_dim + n
          atoms%nlcc_ngv(ityp) = int((sqrt(real(1 + 8 * n)) - 1) / 2)
          n = 0
          if (has_key(nloc, "Conduction")) n = dict_len(nloc // "Conduction")
          nlcc_dim = nlcc_dim + n
          atoms%nlcc_ngc(ityp) = int((sqrt(real(1 + 8 * n)) - 1) / 2)
       end if
       if (has_key(nloc, "Rcore") .and. has_key(nloc, "Core charge")) then
          nlcc_dim=nlcc_dim+1
          atoms%nlcc_ngc(ityp)=1
          atoms%nlcc_ngv(ityp)=0
       end if
    end do
    atoms%donlcc = (nlcc_dim > 0)
    atoms%nlccpar = f_malloc_ptr((/ 0 .to. 4, 1 .to. max(nlcc_dim,1) /), id = "nlccpar")
    !start again the file inspection to fill nlcc parameters
    if (atoms%donlcc) then
       nlcc_dim=0
       fill_nlcc: do ityp=1,atoms%astruct%ntypes
          filename = 'psppar.' // trim(atoms%astruct%atomnames(ityp))
          !ALEX: These are preferably read from psppar.Xy, as stored in the
          !local variables rcore and qcore
          nloc => dict // filename // 'Non Linear Core Correction term'
          if (has_key(nloc, "Valence") .or. has_key(nloc, "Conduction")) then
             n = 0
             if (has_key(nloc, "Valence")) n = dict_len(nloc // "Valence")
             do i = 1, n, 1
                coeffs => nloc // "Valence" // (i - 1)
                atoms%nlccpar(:, nlcc_dim + i) = coeffs
             end do
             nlcc_dim = nlcc_dim + n
             n = 0
             if (has_key(nloc, "Conduction")) n = dict_len(nloc // "Conduction")
             do i = 1, n, 1
                coeffs => nloc // "Conduction" // (i - 1)
                atoms%nlccpar(0, nlcc_dim + i) = coeffs // 0
                atoms%nlccpar(1, nlcc_dim + i) = coeffs // 1
                atoms%nlccpar(2, nlcc_dim + i) = coeffs // 2
                atoms%nlccpar(3, nlcc_dim + i) = coeffs // 3
                atoms%nlccpar(4, nlcc_dim + i) = coeffs // 4
             end do
             nlcc_dim = nlcc_dim + n
          end if
          if (has_key(nloc, "Rcore") .and. has_key(nloc, "Core charge")) then
             nlcc_dim=nlcc_dim+1
             atoms%nlcc_ngc(ityp)=1
             atoms%nlcc_ngv(ityp)=0
             atoms%nlccpar(0,nlcc_dim)=nloc // "Rcore"
             atoms%nlccpar(1,nlcc_dim)=nloc // "Core charge"
             atoms%nlccpar(2:4,nlcc_dim)=0.0_gp 
          end if
       end do fill_nlcc
    end if
  end subroutine nlcc_set_from_dict

  subroutine psp_set_from_dict(dict, nzatom, nelpsp, npspcode, ixcpsp, &
       & psppar, radii_cf)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(out) :: nzatom, nelpsp, npspcode, ixcpsp
    real(gp), intent(out) :: psppar(0:4,0:6), radii_cf(3)

    type(dictionary), pointer :: loc
    character(len = max_field_length) :: str
    integer :: i, l

    nzatom = -1
    radii_cf(:) = UNINITIALIZED(1._gp)
    psppar(:,:) = 0._gp

    ! We set nzatom at the end as a flag that the psp data are complete.
    if (.not. has_key(dict, "No. of Electrons")) return
    nelpsp = dict // "No. of Electrons"
    if (.not. has_key(dict, "Pseudopotential XC")) return
    ixcpsp = dict // "Pseudopotential XC"
    ! Local terms
    if (.not. has_key(dict, "Local Pseudo Potential (HGH convention)")) return
    loc => dict // "Local Pseudo Potential (HGH convention)"
    if (.not. has_key(loc, "Rloc")) return
    psppar(0,0) = loc // 'Rloc'
    if (.not. has_key(loc, "Coefficients (c1 .. c4)")) return
    psppar(0,1) = loc // 'Coefficients (c1 .. c4)' // 0
    psppar(0,2) = loc // 'Coefficients (c1 .. c4)' // 1
    psppar(0,3) = loc // 'Coefficients (c1 .. c4)' // 2
    psppar(0,4) = loc // 'Coefficients (c1 .. c4)' // 3
    ! Nonlocal terms
    if (has_key(dict, "NonLocal PSP Parameters")) then
       do i = 1, dict_len(dict // "NonLocal PSP Parameters"), 1
          loc => dict // "NonLocal PSP Parameters" // (i - 1)
          if (.not. has_key(loc, "Channel (l)")) return
          l = loc // "Channel (l)"
          l = l + 1
          if (.not. has_key(loc, "Rloc")) return
          psppar(l,0) = loc // 'Rloc'
          if (.not. has_key(loc, "h_ij terms")) return
          psppar(l,1) = loc // 'h_ij terms' // 0
          psppar(l,2) = loc // 'h_ij terms' // 1
          psppar(l,3) = loc // 'h_ij terms' // 2
          psppar(l,4) = loc // 'h_ij terms' // 3
          psppar(l,5) = loc // 'h_ij terms' // 4
          psppar(l,6) = loc // 'h_ij terms' // 5
       end do
    end if
    ! Type
    if (.not. has_key(dict, "Pseudopotential type")) return
    str = dict // "Pseudopotential type"
    select case(trim(str))
    case("GTH")
       npspcode = 2
    case("HGH")
       npspcode = 3
    case("HGH-K")
       npspcode = 10
    case("HGH-K + NLCC")
       npspcode = 12
    case default
       return
    end select
    if (npspcode == 12) then
       if (.not. has_key(dict, 'Non Linear Core Correction term')) return
       loc => dict // 'Non Linear Core Correction term'
       if (.not. has_key(loc, "Rcore")) return
       if (.not. has_key(loc, "Core charge")) return
    end if
    ! Valid pseudo, we set nzatom
    if (.not. has_key(dict, "Atomic number")) return
    nzatom = dict // "Atomic number"

    ! Optional values.
    if (has_key(dict, "Radii of active regions (AU)")) then
       loc => dict // "Radii of active regions (AU)"
       if (has_key(loc, "Coarse")) radii_cf(1) =  loc // "Coarse"
       if (has_key(loc, "Fine")) radii_cf(2) =  loc // "Fine"
       if (has_key(loc, "Coarse PSP")) radii_cf(3) =  loc // "Coarse PSP"
    end if
  end subroutine psp_set_from_dict

  subroutine psp_data_merge_to_dict(dict, nzatom, nelpsp, npspcode, ixcpsp, &
       & psppar, radii_cf, rcore, qcore)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: nzatom, nelpsp, npspcode, ixcpsp
    real(gp), intent(in) :: psppar(0:4,0:6), radii_cf(3), rcore, qcore

    type(dictionary), pointer :: channel, radii
    integer :: l, i

    ! Type
    select case(npspcode)
    case(2)
       call set(dict // "Pseudopotential type", 'GTH')
    case(3)
       call set(dict // "Pseudopotential type", 'HGH')
    case(10)
       call set(dict // "Pseudopotential type", 'HGH-K')
    case(12)
       call set(dict // "Pseudopotential type", 'HGH-K + NLCC')
    end select

    call set(dict // "Atomic number", nzatom)
    call set(dict // "No. of Electrons", nelpsp)
    call set(dict // "Pseudopotential XC", ixcpsp)

    ! Local terms
    if (psppar(0,0)/=0) then
       call dict_init(channel)
       call set(channel // 'Rloc', psppar(0,0))
       do i = 1, 4, 1
          call add(channel // 'Coefficients (c1 .. c4)', psppar(0,i))
       end do
       call set(dict // 'Local Pseudo Potential (HGH convention)', channel)
    end if

    ! nlcc term
    if (npspcode == 12) then
       call set(dict // 'Non Linear Core Correction term', &
            & dict_new( 'Rcore' .is. yaml_toa(rcore), &
            & 'Core charge' .is. yaml_toa(qcore)))
    end if

    ! Nonlocal terms
    do l=1,4
       if (psppar(l,0) /= 0._gp) then
          call dict_init(channel)
          call set(channel // 'Channel (l)', l - 1)
          call set(channel // 'Rloc', psppar(l,0))
          do i = 1, 6, 1
             call add(channel // 'h_ij terms', psppar(l,i))
          end do
          call add(dict // 'NonLocal PSP Parameters', channel)
       end if
    end do

    ! Radii (& carottes)
    if (any(radii_cf /= UNINITIALIZED(1._gp))) then
       call dict_init(radii)
       if (radii_cf(1) /= UNINITIALIZED(1._gp)) call set(radii // "Coarse", radii_cf(1))
       if (radii_cf(2) /= UNINITIALIZED(1._gp)) call set(radii // "Fine", radii_cf(2))
       if (radii_cf(3) /= UNINITIALIZED(1._gp)) call set(radii // "Coarse PSP", radii_cf(3))
       call set(dict // "Radii of active regions (AU)", radii)
    end if
  end subroutine psp_data_merge_to_dict

  subroutine atoms_file_merge_to_dict(dict)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use yaml_output, only: yaml_warning
    implicit none
    type(dictionary), pointer :: dict

    type(dictionary), pointer :: types
    character(len = max_field_length) :: str
    integer :: iat
    character(max_field_length), dimension(:), allocatable :: keys
    character(len=27) :: key
    logical :: exists

    ! Loop on types for atomic data.
    call astruct_dict_get_types(dict // "posinp", types)
    allocate(keys(dict_size(types)))
    keys = dict_keys(types)
    do iat = 1, dict_size(types), 1
       key = 'psppar.' // trim(keys(iat))

       exists = has_key(dict, key)
       if (exists) then
          str = dict_value(dict // key)
          if (trim(str) /= "" .and. trim(str) /= TYPE_LIST .and. trim(str) /= TYPE_DICT) then
             call psp_file_merge_to_dict(dict, key, trim(str))
             if (.not. has_key(dict // key, 'Pseudopotential XC')) then
                call yaml_warning("Pseudopotential file '" // trim(str) // &
                     & "' not found. Fallback to file '" // trim(key) // &
                     & "' or hard-coded pseudopotential.")
             end if
          end if
          exists = has_key(dict // key, 'Pseudopotential XC')
       end if
       if (.not. exists) call psp_file_merge_to_dict(dict, key, key)

       exists = has_key(dict, key)
       if (exists) exists = has_key(dict // key, 'Non Linear Core Correction term')
       if (.not. exists) call nlcc_file_merge_to_dict(dict, key, 'nlcc.' // trim(keys(iat)))
    end do
    deallocate(keys)
    call dict_free(types)
  end subroutine atoms_file_merge_to_dict

  subroutine psp_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    integer :: nzatom, nelpsp, npspcode, ixcpsp
    real(gp) :: psppar(0:4,0:6), radii_cf(3), rcore, qcore
    logical :: exists, donlcc, pawpatch

    !ALEX: if npspcode==12, nlccpar are read from psppar.Xy via rcore and qcore 
    call psp_from_file(filename, nzatom, nelpsp, npspcode, ixcpsp, &
         & psppar, donlcc, rcore, qcore, radii_cf, exists, pawpatch)
    if (.not.exists) return

    call psp_data_merge_to_dict(dict // key, nzatom, nelpsp, npspcode, ixcpsp, &
         & psppar, radii_cf, rcore, qcore)
    call set(dict // key // "PAW patch", pawpatch)
    call set(dict // key // "Source", filename)
  end subroutine psp_file_merge_to_dict

  subroutine nlcc_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    type(dictionary), pointer :: psp, gauss
    logical :: exists
    integer :: i, ig, ngv, ngc
    real(gp), dimension(0:4) :: coeffs

    inquire(file=filename,exist=exists)
    if (.not.exists) return

    psp => dict // key

    !read the values of the gaussian for valence and core densities
    open(unit=79,file=filename,status='unknown')
    read(79,*)ngv
    if (ngv > 0) then
       do ig=1,(ngv*(ngv+1))/2
          call dict_init(gauss)
          read(79,*) coeffs
          do i = 0, 4, 1
             call add(gauss, coeffs(i))
          end do
          call add(psp // 'Non Linear Core Correction term' // "Valence", gauss)
       end do
    end if

    read(79,*)ngc
    if (ngc > 0) then
       do ig=1,(ngc*(ngc+1))/2
          call dict_init(gauss)
          read(79,*) coeffs
          do i = 0, 4, 1
             call add(gauss, coeffs(i))
          end do
          call add(psp // 'Non Linear Core Correction term' // "Conduction", gauss)
       end do
    end if

    close(unit=79)
  end subroutine nlcc_file_merge_to_dict

  !> Convert astruct to dictionary for later dump.
  subroutine astruct_merge_to_dict(dict, astruct, rxyz, comment)
    use module_defs, only: gp, UNINITIALIZED, Bohr_Ang
    use module_atoms, only: atomic_structure
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    type(atomic_structure), intent(in) :: astruct
    real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
    character(len = 1024), intent(in), optional :: comment
    !local variables
    type(dictionary), pointer :: pos, at
    integer :: iat,ichg,ispol
    real(gp) :: factor(3)
    logical :: reduced
    character(len = 4) :: frzstr

    !call dict_init(dict)

    reduced = .false.
    factor=1.0_gp
    Units: select case(trim(astruct%units))
    case('angstroem','angstroemd0')
       call set(dict // ASTRUCT_UNITS, 'angstroem')
       factor=Bohr_Ang
    case('reduced')
       call set(dict // ASTRUCT_UNITS, 'reduced')
       reduced = .true.
    case('atomic','atomicd0','bohr','bohrd0')
       ! Default, store nothing
    end select Units

    !cell information
    BC :select case(astruct%geocode)
    case('S')
       call set(dict // ASTRUCT_CELL // 0, yaml_toa(astruct%cell_dim(1)*factor(1)))
       call set(dict // ASTRUCT_CELL // 1, '.inf')
       call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
       !angdeg to be added
       if (reduced) then
          factor(1) = 1._gp / astruct%cell_dim(1)
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('W')
       call set(dict // ASTRUCT_CELL // 0, '.inf')
       call set(dict // ASTRUCT_CELL // 1, '.inf')
       call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
       if (reduced) then
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('P')
       call set(dict // ASTRUCT_CELL // 0, yaml_toa(astruct%cell_dim(1)*factor(1)))
       call set(dict // ASTRUCT_CELL // 1, yaml_toa(astruct%cell_dim(2)*factor(2)))
       call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
       !angdeg to be added
       if (reduced) then
          factor(1) = 1._gp / astruct%cell_dim(1)
          factor(2) = 1._gp / astruct%cell_dim(2)
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('F')
       ! Default, store nothing and erase key if already exist.
       if (has_key(dict, ASTRUCT_CELL)) call pop(dict, ASTRUCT_CELL)
    end select BC

    if (has_key(dict, ASTRUCT_POSITIONS)) call pop(dict, ASTRUCT_POSITIONS)
    pos => dict // ASTRUCT_POSITIONS
    do iat=1,astruct%nat
       call dict_init(at)
       call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(1,iat) * factor(1))
       call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(2,iat) * factor(2))
       call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(3,iat) * factor(3))
       if (astruct%ifrztyp(iat) /= 0) then
          call frozen_itof(astruct%ifrztyp(iat), frzstr)
          call set(at // "Frozen", adjustl(frzstr))
       end if
       call charge_and_spol(astruct%input_polarization(iat),ichg,ispol)
       if (ichg /= 0) call set(at // "IGChg", ichg)
       if (ispol /= 0) call set(at // "IGSpin", ispol)
       call add(pos, at)
    end do

    if (present(comment)) then
       if (len_trim(comment) > 0) &
            & call add(dict // ASTRUCT_PROPERTIES // "info", comment)
    end if

    if (len_trim(astruct%inputfile_format) > 0) &
         & call set(dict // ASTRUCT_PROPERTIES // "format", astruct%inputfile_format)
  end subroutine astruct_merge_to_dict
  
  subroutine astruct_dict_get_types(dict, types)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict, types

    type(dictionary), pointer :: atoms, at
    character(len = max_field_length) :: str
    integer :: iat, ityp

    call dict_init(types)
    atoms => dict // ASTRUCT_POSITIONS
    ityp = 0
    do iat = 1, dict_len(atoms), 1
       at => dict_iter(atoms // (iat - 1))
       do while(associated(at))
          str = dict_key(at)
          if (dict_len(at) == 3 .and. .not. has_key(types, str)) then
             ityp = ityp + 1
             call set(types // str, ityp)
             nullify(at)
          else
             at => dict_next(at)
          end if
       end do
    end do
  end subroutine astruct_dict_get_types

  subroutine astruct_dict_get_source(dict, source)
    use dictionaries, only: max_field_length, dictionary, has_key, operator(//), dict_value
    implicit none
    type(dictionary), pointer :: dict
    character(len = max_field_length), intent(out) :: source
    
    write(source, "(A)") ""
    if (has_key(dict, ASTRUCT_PROPERTIES)) then
       if (has_key(dict // ASTRUCT_PROPERTIES, "source")) &
            & source = dict_value(dict // ASTRUCT_PROPERTIES // "source")
    end if
  end subroutine astruct_dict_get_source

  subroutine astruct_file_merge_to_dict(dict, key, filename)
    use module_base, only: gp, UNINITIALIZED, bigdft_mpi,f_routine,f_release_routine
    use module_atoms, only: set_astruct_from_file,atomic_structure,&
         nullify_atomic_structure,deallocate_atomic_structure
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key
    
    type(atomic_structure) :: astruct
    integer :: ierr
    call f_routine(id='astruct_file_merge_to_dict')
    ! Read atomic file, old way
    call nullify_atomic_structure(astruct)
    call set_astruct_from_file(filename, bigdft_mpi%iproc, astruct, status = ierr)
    if (ierr == 0) then
       call astruct_merge_to_dict(dict // key, astruct, astruct%rxyz)
       call set(dict // key // ASTRUCT_PROPERTIES // "source", filename)
       call deallocate_atomic_structure(astruct)
    end if
    call f_release_routine()
  end subroutine astruct_file_merge_to_dict

  !> allocate the astruct variable from the dictionary of input data
  !retrieve also other information like the energy and the forces if requested
  !! and presend in the dictionary
  subroutine astruct_set_from_dict(dict, astruct, comment)
    use module_defs, only: gp, Bohr_Ang, UNINITIALIZED
    use module_atoms, only: atomic_structure, nullify_atomic_structure
    use dictionaries
    use dynamic_memory
    implicit none
    type(dictionary), pointer :: dict !< dictionary of the input variables
    !! the keys have to be declared like input_dicts module
    type(atomic_structure), intent(out) :: astruct !<structure created from the file
    !> extra comment retrieved from the file if present
    character(len = 1024), intent(out), optional :: comment

    !local variables
    character(len=*), parameter :: subname='astruct_set_from_dict'
    type(dictionary), pointer :: pos, at, types
    character(len = max_field_length) :: str
    integer :: iat, ityp, units, igspin, igchrg, nsgn, ntyp

    call nullify_atomic_structure(astruct)
    astruct%nat = -1
    if (present(comment)) write(comment, "(A)") " "
    if (.not. has_key(dict, ASTRUCT_POSITIONS)) return

    ! The units
    units = 0
    write(astruct%units, "(A)") "bohr"
    if (has_key(dict, ASTRUCT_UNITS)) astruct%units = dict // ASTRUCT_UNITS
    select case(trim(astruct%units))
    case('atomic','atomicd0','bohr','bohrd0')
       units = 0
    case('angstroem','angstroemd0')
       units = 1
    case('reduced')
       units = 2
    end select
    ! The cell
    astruct%cell_dim = 0.0_gp
    if (.not. has_key(dict, ASTRUCT_CELL)) then
       astruct%geocode = 'F'
    else
       astruct%geocode = 'P'
       ! z
       astruct%cell_dim(3) = dict // ASTRUCT_CELL // 2
       ! y
       str = dict // ASTRUCT_CELL // 1
       if (trim(str) == ".inf") then
          astruct%geocode = 'S'
       else
          astruct%cell_dim(2) = dict // ASTRUCT_CELL // 1
       end if
       ! x
       str = dict // ASTRUCT_CELL // 0
       if (trim(str) == ".inf") then
          astruct%geocode = 'W'
       else
          astruct%cell_dim(1) = dict // ASTRUCT_CELL // 0
       end if
    end if
    if (units == 1) astruct%cell_dim = astruct%cell_dim / Bohr_Ang
    ! The types
    call astruct_dict_get_types(dict, types)
    ntyp = dict_size(types)
    call astruct_set_n_types(astruct, ntyp)
    ! astruct%atomnames = dict_keys(types)
    ityp = 1
    at => dict_iter(types)
    do while (associated(at))
       astruct%atomnames(ityp) = dict_key(at)
       ityp = ityp + 1
       at => dict_next(at)
    end do
    ! The atoms
    pos => dict // ASTRUCT_POSITIONS
    call astruct_set_n_atoms(astruct, dict_len(pos))
    do iat = 1, astruct%nat
       igspin = 0
       igchrg = 0
       nsgn   = 1
       !at => pos // (iat - 1)
       at => dict_iter(pos//(iat-1))!at%child
       do while(associated(at))
          str = dict_key(at)
          if (trim(str) == "Frozen") then
             str = dict_value(at)
             call frozen_ftoi(str(1:4), astruct%ifrztyp(iat))
          else if (trim(str) == "IGSpin") then
             igspin = at
          else if (trim(str) == "IGChg") then
             igchrg = at
             if (igchrg >= 0) then
                nsgn = 1
             else
                nsgn = -1
             end if
          else if (dict_len(at) == 3) then
             astruct%iatype(iat) = types // dict_key(at)
             astruct%rxyz(:, iat) = at
          end if
          at => dict_next(at)
       end do
       astruct%input_polarization(iat) = 1000 * igchrg + nsgn * 100 + igspin
       if (units == 1) then
          astruct%rxyz(1,iat) = astruct%rxyz(1,iat) / Bohr_Ang
          astruct%rxyz(2,iat) = astruct%rxyz(2,iat) / Bohr_Ang
          astruct%rxyz(3,iat) = astruct%rxyz(3,iat) / Bohr_Ang
       endif
       if (units == 2) then !add treatment for reduced coordinates
          if (astruct%cell_dim(1) > 0.) astruct%rxyz(1,iat)=&
               modulo(astruct%rxyz(1,iat),1.0_gp) * astruct%cell_dim(1)
          if (astruct%cell_dim(2) > 0.) astruct%rxyz(2,iat)=&
               modulo(astruct%rxyz(2,iat),1.0_gp) * astruct%cell_dim(2)
          if (astruct%cell_dim(3) > 0.) astruct%rxyz(3,iat)=&
               modulo(astruct%rxyz(3,iat),1.0_gp) * astruct%cell_dim(3)
       else if (astruct%geocode == 'P') then
          astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
          astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),astruct%cell_dim(2))
          astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
       else if (astruct%geocode == 'S') then
          astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
          astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
       else if (astruct%geocode == 'W') then
          astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
       end if
    end do

    if (has_key(dict, ASTRUCT_PROPERTIES)) then
       pos => dict // ASTRUCT_PROPERTIES
       if (has_key(pos, "info") .and. present(comment)) comment = pos // "info"
       if (has_key(pos, "format")) astruct%inputfile_format = pos // "format"
    end if

  end subroutine astruct_set_from_dict

  subroutine aocc_to_dict(dict, nspin, noncoll, nstart, aocc, nelecmax, lmax, nsccode)
    use module_defs, only: gp
    use dictionaries
    implicit none
    integer, intent(in) :: nelecmax, lmax, nsccode, nspin, noncoll, nstart
    type(dictionary), pointer :: dict
    real(gp), dimension(nelecmax), intent(in) :: aocc

    type(dictionary), pointer :: val
    character(len = 4) :: key
    integer :: l, inl, nl, iocc, sccode, nsc, i
    character(len = 1), dimension(4), parameter :: lname = (/ "s", "p", "d", "f" /)

    call dict_init(dict)

    sccode = nsccode
    iocc=0
    do l = 1, lmax
       iocc=iocc+1
       ! Get number of shells for this channel
       nl = int(aocc(iocc))
       ! Get number of semi cores for this channel
       nsc = modulo(sccode, 4)
       sccode = sccode / 4
       if (nl == 0) cycle
       do inl = 1, nl, 1
          if (inl <= nsc) then
             write(key, "(A1,I1,A1,A1)") "(", nstart + inl, lname(l), ")"
          else
             write(key, "(I1, A1)") nstart + inl, lname(l)
          end if
          call dict_init(val)
          do i = 1, nspin * noncoll * (2 * l - 1), 1
             iocc=iocc+1
             call add(val, aocc(iocc))
          end do
          call set(dict // key, val)
       end do
    end do
  end subroutine aocc_to_dict

  subroutine atomic_data_set_from_dict(dict, key, atoms, nspin)
    use module_defs, only: gp
    use ao_inguess, only: ao_ig_charge,atomic_info,aoig_set_from_dict,&
         print_eleconf,aoig_set
    use module_types, only: atoms_data
    use dictionaries
!    use dynamic_memory
    use yaml_output, only: yaml_warning, yaml_toa
    implicit none
    type(dictionary), pointer :: dict
    type(atoms_data), intent(inout) :: atoms
    character(len = *), intent(in) :: key
    integer, intent(in) :: nspin

    integer :: iat, ityp
    real(gp) :: rcov,elec!,rprb,ehomo,elec
    character(len = max_field_length) :: at
    type(dictionary), pointer :: dict_tmp

    do ityp = 1, atoms%astruct%ntypes, 1
       !only amu and rcov are extracted here
       call atomic_info(atoms%nzatom(ityp),atoms%nelpsp(ityp),&
            amu=atoms%amu(ityp),rcov=rcov)
!       atoms%rloc(ityp,:) = rcov * 10.0

       do iat = 1, atoms%astruct%nat, 1
          if (atoms%astruct%iatype(iat) /= ityp) cycle

          !fill the atomic IG configuration from the input_polarization
          atoms%aoig(iat)=aoig_set(atoms%nzatom(ityp),atoms%nelpsp(ityp),&
               atoms%astruct%input_polarization(iat),nspin)

          ! Possible overwrite, if the dictionary has the item
          if (has_key(dict, key)) then
             nullify(dict_tmp)
             at(1:len(at))="Atom "//trim(adjustl(yaml_toa(iat)))
             if (has_key(dict // key,trim(at))) &
                  dict_tmp=>dict//key//trim(at)
             if (has_key(dict // key, trim(atoms%astruct%atomnames(ityp)))) &
                  dict_tmp=>dict // key // trim(atoms%astruct%atomnames(ityp))
             if (associated(dict_tmp)) then
                atoms%aoig(iat)=aoig_set_from_dict(dict_tmp,nspin)
                !check the total number of electrons
                elec=ao_ig_charge(nspin,atoms%aoig(iat)%aocc)
                if (nint(elec) /= atoms%nelpsp(ityp)) then
                   call print_eleconf(nspin,atoms%aoig(iat)%aocc,atoms%aoig(iat)%iasctype)
                   call yaml_warning('The total atomic charge '//trim(yaml_toa(elec))//&
                        ' is different from the PSP charge '//trim(yaml_toa(atoms%nelpsp(ityp))))
                end if
             end if
          end if
       end do

    end do

    !number of atoms with semicore channels
    atoms%natsc = 0
    do iat=1,atoms%astruct%nat
       if (atoms%aoig(iat)%iasctype /= 0) atoms%natsc=atoms%natsc+1
    enddo
  end subroutine atomic_data_set_from_dict

  subroutine atomic_data_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    logical :: exists
    integer :: ierror, jat, nsp, nsccode
    character(len = 1024) :: string
    character(len = max_field_length) :: at
    integer, parameter :: nelecmax = 32, noccmax = 4, lmax = 4
    real(gp), dimension(nelecmax) :: aocc
    type(dictionary), pointer :: val
    
    inquire(file = filename, exist = exists)
    if (.not. exists) return

    open(unit=91,file=filename,status='old',iostat=ierror)
    !Check the open statement
    if (f_err_raise(ierror /= 0,'Failed to open the existing file '// trim(filename),&
         err_name='BIGDFT_RUNTIME_ERROR')) return

    parse_inocc: do
       read(91,'(a1024)',iostat=ierror)string
       if (ierror /= 0) exit parse_inocc !file ends
       read(string,*,iostat=ierror)jat
       if (ierror /=0) stop 'Error reading line'

       write(at, "(A, I0)") "Atom ", jat
       call read_eleconf(string,noccmax,nelecmax,lmax,aocc,nsccode,nsp)
       call aocc_to_dict(val, nsp, 1, 0, aocc, nelecmax, lmax, nsccode)
       call set(dict // key // at, val)
    end do parse_inocc

    close(unit = 91)

  end subroutine atomic_data_file_merge_to_dict

  subroutine occupation_set_from_dict(dict, key, norbu, norbd, occup, &
       & nkpts, nspin, norbsempty, nelec_up, nelec_down, norb_max)
    use module_defs, only: gp
    use dictionaries
    use dynamic_memory
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    real(gp), dimension(:), pointer :: occup
    integer, intent(in) :: nkpts, nspin, norbsempty, nelec_up, nelec_down, norb_max
    integer, intent(out) :: norbu, norbd

    integer :: norb
    integer :: ikpt
    type(dictionary), pointer :: occup_src
    character(len = 12) :: kpt_key

    ! Default case.
    if (nspin == 1) then
       norb  = min((nelec_up + 1) / 2, norb_max)
       norbu = norb
    else
       norb = min(nelec_up + nelec_down, 2 * norb_max)
       if (nspin == 2) then
          norbu = min(nelec_up, norb_max)
       else
          norbu = min(nelec_up, 2 * norb_max)
       end if
    end if
    norbd = norb - norbu
!!$    write(*,*) nelec_up, nelec_down, norbsempty, norb_max
!!$    write(*,*) norbu, norbd, norb
!!$    stop
    ! Modify the default with occupation
    nullify(occup_src)
    if (has_key(dict, key)) then
       occup_src => dict //key
       ! Occupation is provided.
       if (has_key(occup_src, "K point 1")) then
          call count_for_kpt(occup_src // "K point 1")
       else if (nkpts == 1) then
          call count_for_kpt(occup_src)
       end if
       do ikpt = 2, nkpts, 1
          write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
          if (has_key(occup_src, kpt_key)) call count_for_kpt(occup_src // kpt_key)
       end do
    else if (norbsempty > 0) then
       !value of empty orbitals up and down, needed to fill occupation numbers
       if (nspin == 4 .or. nspin == 1) then
          norbu = norbu + min(norbsempty, norb_max - norbu)
       else if (nspin == 2) then
          norbu = norbu + min(norbsempty, norb_max - norbu)
          norbd = norbd + min(norbsempty, norb_max - norbd)
       end if
    end if

    ! Summarize and check.
    norb = norbu + norbd
    if (((nspin == 1 .or. nspin == 2) .and. (norbu > norb_max .or. norbd > norb_max)) &
         & .or. (nspin == 4 .and. (norbu > 2 * norb_max .or. norbd > 0))) then
       call yaml_warning('Total number of orbitals (found ' // trim(yaml_toa(norb)) &
            & // ') exceeds the available input guess orbitals (being ' &
            & // trim(yaml_toa(norb_max)) // ').')
       stop
    end if

    ! Allocate occupation accordingly.
    occup = f_malloc_ptr(norb * nkpts, id = "occup", routine_id = "occupation_set_from_dict")
    ! Setup occupation
    if (nspin==1) then
       do ikpt = 1, nkpts, 1
          call fill_default((ikpt - 1) * norb, 2, nelec_up, norb)
          if (associated(occup_src)) then
             write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
             if (ikpt == 0 .and. .not. has_key(occup_src, kpt_key)) then
                call fill_for_kpt((ikpt - 1) * norb, occup_src)
             else
                call fill_for_kpt((ikpt - 1) * norb, occup_src // kpt_key)
             end if
          end if
       end do
    else
       do ikpt = 0, nkpts - 1, 1
          call fill_default(ikpt * norb, 1, nelec_up, norbu)
          call fill_default(ikpt * norb + norbu, 1, nelec_down, norbd)
          if (associated(occup_src)) then
             write(kpt_key, "(A)") "K point" // trim(yaml_toa(ikpt, fmt = "(I0)"))
             if (ikpt == 0 .and. .not. has_key(occup_src, kpt_key)) then
                call fill_for_kpt((ikpt - 1) * norb, occup_src // "up")
                call fill_for_kpt((ikpt - 1) * norb + norbu, occup_src // "down")
             else
                call fill_for_kpt((ikpt - 1) * norb, occup_src // kpt_key // "up")
                call fill_for_kpt((ikpt - 1) * norb + norbu, occup_src // kpt_key // "down")
             end if
          end if
       end do
    end if

    !Check if sum(occup)=nelec
    if (abs(sum(occup) / nkpts - real(nelec_up + nelec_down,gp))>1.e-6_gp) then
       call yaml_warning('the total number of electrons ' &
            & // trim(yaml_toa(sum(occup) / nkpts,fmt='(f13.6)')) &
            & // ' is not equal to' // trim(yaml_toa(nelec_up + nelec_down)))
       stop
    end if

  contains

    subroutine count_for_kpt(occ)
      implicit none
      type(dictionary), pointer :: occ
      
      if (nspin == 2) then
         if (.not. has_key(occ, "up") .or. &
              & .not. has_key(occ, "down")) stop "missing up or down"
         call count_orbs(norbu, occ // "up")
         call count_orbs(norbd, occ // "down")
      else
         call count_orbs(norbu, occ)
      end if
    end subroutine count_for_kpt

    subroutine count_orbs(n, occ)
      implicit none
      type(dictionary), pointer :: occ
      integer, intent(inout) :: n
      
      type(dictionary), pointer :: it
      character(len = max_field_length) :: key
      integer :: iorb

      it => dict_iter(occ)
      do while(associated(it))
         key = dict_key(it)
         read(key(index(key, " ") + 1:), *) iorb
         n = max(n, iorb)
         it => dict_next(it)
      end do
    end subroutine count_orbs

    subroutine fill_default(isorb, nfill, nelec, norb)
      implicit none
      integer, intent(in) :: isorb, nfill, nelec, norb

      integer :: nt, it, iorb, ne

      nt=0
      ne = (nelec + 1) / nfill
      do iorb=isorb + 1, isorb + min(ne, norb)
         it=min(nfill,nelec-nt)
         occup(iorb)=real(it,gp)
         nt=nt+it
      enddo
      do iorb=isorb+min(ne, norb)+1,isorb+norb
         occup(iorb)=0._gp
      end do
    end subroutine fill_default

    subroutine fill_for_kpt(isorb, occ)
      implicit none
      integer, intent(in) :: isorb
      type(dictionary), pointer :: occ

      type(dictionary), pointer :: it
      character(len = max_field_length) :: key
      integer :: iorb

      it => dict_iter(occ)
      do while(associated(it))
         key = dict_key(it)
         read(key(index(key, " ") + 1:), *) iorb
         occup(isorb + iorb) = it
         it => dict_next(it)
      end do
    end subroutine fill_for_kpt
  end subroutine occupation_set_from_dict

  subroutine occupation_data_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    logical :: exists
    integer :: ierror, ntu, ntd, nt, i, iorb
    character(len = 100) :: line, string
    type(dictionary), pointer :: valu, vald
    
    inquire(file = filename, exist = exists)
    if (.not. exists) return

    open(unit=91,file=filename,status='old',iostat=ierror)
    !Check the open statement
    if (ierror /= 0) then
       call yaml_warning('Failed to open the existing file '// trim(filename))
       stop
    end if

    !The first line gives the number of orbitals
    read(unit=91,fmt='(a100)') line

    read(line,fmt=*,iostat=ierror) ntu, ntd
    if (ierror /= 0) then
       !The first line gives the number of orbitals
       ntd = 0
       read(line,fmt=*,iostat=ierror) ntu
       if (ierror /=0) stop 'ERROR: reading the number of orbitals.'
    end if

    call dict_init(valu)
    if (ntd > 0) call dict_init(vald)

    nt = 1
    do
       read(unit=91,fmt='(a100)',iostat=ierror) line
       if (ierror /= 0) then
          exit
       end if
       !Transform the line in case there are slashes (to ease the parsing)
       do i=1,len(line)
          if (line(i:i) == '/') then
             line(i:i) = ':'
          end if
       end do
       read(line,*,iostat=ierror) iorb,string
       if (ierror /= 0) then
          exit
       end if
       !Transform back the ':' into '/'
       do i=1,len(string)
          if (string(i:i) == ':') then
             string(i:i) = '/'
          end if
       end do
       nt=nt+1

       if (iorb<0 .or. iorb>ntu + ntd) then
          !if (iproc==0) then
          write(*,'(1x,a,i0,a)') 'ERROR in line ',nt,' of the file "[name].occ"'
          write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
          !end if
          stop
       else
          if (iorb <= ntu) then
             call set(valu // ("Orbital" // trim(yaml_toa(iorb, fmt = "(I0)"))), string)
          else
             call set(vald // ("Orbital" // trim(yaml_toa(iorb, fmt = "(I0)"))), string)
          end if
       end if
    end do

    close(unit = 91)

    if (ntd > 0) then
       call set(dict // key // "K point 1" // "up", valu)
       call set(dict // key // "K point 1" // "down", vald)
    else
       call set(dict // key // "K point 1", valu)
    end if

    call set(dict // key // "Source", filename)

  end subroutine occupation_data_file_merge_to_dict

  subroutine global_output_merge_to_dict(dict, outs, astruct)
    use module_defs, only: gp, UNINITIALIZED
    use module_types, only: atomic_structure, DFT_global_output
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    type(DFT_global_output), intent(in) :: outs
    type(atomic_structure), intent(in) :: astruct

    integer :: iat
    type(dictionary), pointer :: pos, fxyz

    if (has_key(dict, GOUT_FORCES)) call pop(dict, GOUT_FORCES)
    pos => dict // GOUT_FORCES
    do iat=1,astruct%nat
       call dict_init(fxyz)
       call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 0, outs%fxyz(1, iat))
       call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 1, outs%fxyz(2, iat))
       call set(fxyz // astruct%atomnames(astruct%iatype(iat)) // 2, outs%fxyz(3, iat))
       call add(pos, fxyz)
    end do

    call set(dict // GOUT_ENERGY, outs%energs%eKS)

  end subroutine global_output_merge_to_dict
end module module_input_dicts

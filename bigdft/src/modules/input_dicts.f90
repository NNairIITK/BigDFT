!> @file
!>  Modules which contains all the interfaces to parse input dictionary.
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Modules which contains all interfaces to parse input dictionary.
module module_input_dicts
  use public_keys
  implicit none

  private

  ! Update a dictionary from a input file
  public :: merge_input_file_to_dict

  ! Main creation routine
  public :: user_dict_from_files

  ! Dictionary completion
  public :: psp_dict_fill_all, psp_dict_analyse

  ! Dictionary inquire
  public :: astruct_dict_get_types

  ! Types from dictionaries
  public :: astruct_set_from_dict
  public :: psp_set_from_dict, nlcc_set_from_dict
  public :: occupation_set_from_dict
  public :: neb_set_from_dict

  ! Types to dictionaries
  public :: psp_data_merge_to_dict

  ! Dictionaries from files (old formats).
  public :: psp_file_merge_to_dict, nlcc_file_merge_to_dict
  public :: atoms_file_merge_to_dict
  public :: astruct_file_merge_to_dict
  public :: occupation_data_file_merge_to_dict

contains

  !> logical function, true if file is existing (basically a fortran inquire)
  function file_exists(filename) result(exists)
    use module_base, only: f_err_raise
    use yaml_output, only: yaml_toa
    implicit none
    !Arguments
    character(len=*), intent(in) :: filename !< Filename
    logical :: exists                        !< Result
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
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict            !< Dictionary of the input files. Should be initialized on entry
    character(len = *), intent(in) :: fname      !< Name of the file where the dictionary has to be read from 
    type(mpi_environment), intent(in) :: mpi_env !< Environment of the reading. Used for broadcasting the result
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
       !if (mpi_env%nproc > 1) &
       !     & call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    !else
       !call mpi_bcast(cbuf_len, 1, MPI_INTEGER8, 0, mpi_env%mpi_comm, ierr)
    end if

    if (mpi_env%nproc > 1) call mpibcast(cbuf_len,comm=mpi_env%mpi_comm)
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')

    if (mpi_env%iproc == 0) then
       call copyCBuffer(fbuf, cbuf, cbuf_len)
       call freeCBuffer(cbuf)
!       if (mpi_env%nproc > 1 .and. cbuf_len > 0) &
!            & call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
!    else
!       if (cbuf_len > 0) call mpi_bcast(fbuf(1), int(cbuf_len), MPI_CHARACTER, 0, mpi_env%mpi_comm, ierr)
    end if

    !this call can be replaced with the size of the character array
    if (mpi_env%nproc > 1) call mpibcast(fbuf,comm=mpi_env%mpi_comm)

    call f_err_open_try()
    call yaml_parse_from_char_array(udict, fbuf)
    call f_free_str(1,fbuf)
    ! Handle with possible partial dictionary.
    if (dict_len(udict) > 0) then
       call dict_update(dict, udict // 0)
    end if
    call dict_free(udict)
    ierr = 0
    if (f_err_check()) ierr = f_get_last_error(val)
    !call f_dump_all_errors()
    call f_err_close_try()

    if (ierr /= 0) call f_err_throw(err_id = ierr, err_msg = val)
    call f_release_routine()

  END SUBROUTINE merge_input_file_to_dict

  !> Read from all input files and build a dictionary
  subroutine user_dict_from_files(dict,radical,posinp_name, mpi_env)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use module_defs, only: mpi_environment
    use module_interfaces, only: read_input_dict_from_files
    use public_keys, only: POSINP,IG_OCCUPATION
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict                  !< Contains (out) all the information
    character(len = *), intent(in) :: radical          !< Radical for the input files
    character(len = *), intent(in) :: posinp_name           !< If the dict has no posinp key, use it
    type(mpi_environment), intent(in) :: mpi_env       !< MPI Environment
    !Local variables
    type(dictionary), pointer :: at
    character(len = max_field_length) :: str

    !read the input file(s) and transform them into a dictionary
    call read_input_dict_from_files(trim(radical), mpi_env, dict)

    !possible overwrite with a specific posinp file.
    call astruct_file_merge_to_dict(dict,POSINP, trim(posinp_name))

    if (has_key(dict,POSINP)) then
       str = dict_value(dict //POSINP)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          !str contains a file name so add atomic positions from it.
          call astruct_file_merge_to_dict(dict,POSINP, trim(str))
       else
          !The yaml file contains the atomic positions
          !Only add the format
          at => dict //POSINP
          if (.not. has_key(at, ASTRUCT_PROPERTIES)) then
             call set(at // ASTRUCT_PROPERTIES // FORMAT_KEY, FORMAT_YAML)
          else
             at => at // ASTRUCT_PROPERTIES
             if (FORMAT_KEY .notin. at) &
                  call set(at // FORMAT_KEY, FORMAT_YAML)
          end if
       end if
    end if

    ! Add old psppar
    call atoms_file_merge_to_dict(dict)

    !when the user has not specified the occupation in the input file
    if (.not. has_key(dict,IG_OCCUPATION)) then
       !yaml format should be used even for old method
       if (file_exists(trim(radical)//".occup")) &
            call merge_input_file_to_dict(dict//IG_OCCUPATION,&
            trim(radical)//".occup",mpi_env)
    else !otherwise the input file always supersedes
       str = dict_value(dict //IG_OCCUPATION)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          !call atomic_data_file_merge_to_dict(dict, ATOMIC_OCC, trim(str))
          if (file_exists(trim(str))) &
               call merge_input_file_to_dict(dict//IG_OCCUPATION,trim(str),mpi_env)
       end if
    end if

    if (OCCUPATION .notin. dict) then
       ! Add old input.occ
       call occupation_data_file_merge_to_dict(dict,OCCUPATION,trim(radical) // ".occ")
    else
       str = dict_value(dict //OCCUPATION)
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          call occupation_data_file_merge_to_dict(dict,OCCUPATION, trim(str))
       end if
    end if

  end subroutine user_dict_from_files


  !> Fill up the dict with all pseudopotential information
  subroutine psp_dict_fill_all(dict, atomname, run_ixc, projrad, crmult, frmult)
    use module_defs, only: gp, UNINITIALIZED
    use ao_inguess, only: atomic_info
    use module_atoms, only : RADII_SOURCE, RADII_SOURCE_HARD_CODED, RADII_SOURCE_FILE
    use dictionaries
    use dynamic_memory
    implicit none
    !Arguments
    type(dictionary), pointer :: dict          !< Input dictionary (inout)
    character(len = *), intent(in) :: atomname !< Atom name
    integer, intent(in) :: run_ixc             !< XC functional
    real(gp), intent(in) :: projrad            !< projector radius
    real(gp), intent(in) :: crmult, frmult     !< radius multipliers
    !Local variables
    integer :: ixc
    !integer :: ierr
    character(len=27) :: filename
    logical :: exists
    integer :: nzatom, nelpsp, npspcode
    real(gp), dimension(0:4,0:6) :: psppar
    integer :: i,nlen
    real(gp) :: ehomo,radfine,rad,maxrad
    type(dictionary), pointer :: radii,dict_psp
    real(gp), dimension(3) :: radii_cf
    character(len = max_field_length) :: source_val

    filename = 'psppar.' // atomname
    dict_psp => dict // filename !inquire for the key?


    exists = has_key(dict_psp, LPSP_KEY)
    if (.not. exists) then
       if (dict_len(dict_psp) > 0) then
          ! Long string case, we parse it.
          call psp_file_merge_to_dict(dict, filename, lstring = dict_psp)
          ! Since it has been overrided.
          dict_psp => dict // filename
          exists = has_key(dict_psp, LPSP_KEY)
          nzatom = dict_psp .get. ATOMIC_NUMBER
          nelpsp = dict_psp .get. ELECTRON_NUMBER
       else
          ixc = run_ixc
          ixc = dict_psp .get. PSPXC_KEY
          call psp_from_data(atomname, nzatom, &
               & nelpsp, npspcode, ixc, psppar(:,:), exists)
          radii_cf(:) = UNINITIALIZED(1._gp)
          call psp_data_merge_to_dict(dict_psp, nzatom, nelpsp, npspcode, ixc, &
               & psppar(0:4,0:6), radii_cf, UNINITIALIZED(1._gp), UNINITIALIZED(1._gp))
          call set(dict_psp // SOURCE_KEY, "Hard-Coded")
       end if
    else
       nzatom = dict_psp // ATOMIC_NUMBER
       nelpsp = dict_psp // ELECTRON_NUMBER
    end if

    if (.not. exists) then
       call f_err_throw('The pseudopotential parameter file "'//&
            trim(filename)//&
            '" is lacking, and no registered pseudo found for "'//&
            trim(atomname),err_name='BIGDFT_INPUT_FILE_ERROR')
       return
    end if

    radii_cf = UNINITIALIZED(1._gp)
    !example with the .get. operator
!    print *,'here',associated(radii)
    nullify(radii)
    radii = dict_psp .get. RADII_KEY
    radii_cf(1) = radii .get. COARSE
    radii_cf(2) = radii .get. FINE
    radii_cf(3) = radii .get. COARSE_PSP

    write(source_val, "(A)") RADII_SOURCE(RADII_SOURCE_FILE)
    if (radii_cf(1) == UNINITIALIZED(1.0_gp)) then
       !see whether the atom is semicore or not
       !and consider the ground state electronic configuration
       call atomic_info(nzatom,nelpsp,ehomo=ehomo)
       !call eleconf(nzatom, nelpsp,symbol,rcov,rprb,ehomo,&
       !     neleconf,nsccode,mxpl,mxchg,amu)

       !assigning the radii by calculating physical parameters
       radii_cf(1)=1._gp/sqrt(abs(2._gp*ehomo))
       write(source_val, "(A)") RADII_SOURCE(RADII_SOURCE_HARD_CODED)
    end if
    if (radii_cf(2) == UNINITIALIZED(1.0_gp)) then
       radfine = dict_psp // LPSP_KEY // "Rloc"
       if (has_key(dict_psp, NLPSP_KEY)) then
          nlen=dict_len(dict_psp // NLPSP_KEY)
          do i=1, nlen
             rad = dict_psp // NLPSP_KEY // (i - 1) // "Rloc"
             if (rad /= 0._gp) then
                radfine=min(radfine, rad)
             end if
          end do
       end if
       radii_cf(2)=radfine
       write(source_val, "(A)") RADII_SOURCE(RADII_SOURCE_HARD_CODED)
    end if
    if (radii_cf(3) == UNINITIALIZED(1.0_gp)) radii_cf(3)=crmult*radii_cf(1)/frmult
    ! Correct radii_cf(3) for the projectors.
    maxrad=0.e0_gp ! This line added by Alexey, 03.10.08, to be able to compile with -g -C
    if (has_key( dict_psp, NLPSP_KEY)) then
       nlen=dict_len(dict_psp // NLPSP_KEY)
       do i=1, nlen
          rad =  dict_psp  // NLPSP_KEY // (i - 1) // "Rloc"
          if (rad /= 0._gp) then
             maxrad=max(maxrad, rad)
          end if
       end do
    end if
    if (maxrad == 0.0_gp) then
       radii_cf(3)=0.0_gp
    else
       radii_cf(3)=max(min(radii_cf(3),projrad*maxrad/frmult),radii_cf(2))
    end if
    radii => dict_psp // RADII_KEY
    call set(radii // COARSE, radii_cf(1))
    call set(radii // FINE, radii_cf(2))
    call set(radii // COARSE_PSP, radii_cf(3))
    call set(radii // SOURCE_KEY, source_val)
    
  end subroutine psp_dict_fill_all

  
  !> Fill up the atoms structure from dict
  subroutine psp_dict_analyse(dict, atoms)
    use module_defs, only: gp
    use module_types, only: atoms_data
    use module_atoms, only: allocate_atoms_data
    use dictionaries
    use m_pawrad, only: pawrad_type, pawrad_nullify
    use m_pawtab, only: pawtab_type, pawtab_nullify
    use psp_projectors, only: PSPCODE_PAW
    use public_keys, only: SOURCE_KEY
    implicit none
    !Arguments
    type(dictionary), pointer :: dict        !< Input dictionary
    type(atoms_data), intent(inout) :: atoms !Atoms structure to fill up
    !Local variables
    integer :: ityp, ityp2
    character(len = 27) :: filename
    real(gp), dimension(3) :: radii_cf
    logical :: pawpatch, l
    integer :: paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices
    character(len = max_field_length) :: fpaw

    if (.not. associated(atoms%nzatom)) then
       call allocate_atoms_data(atoms)
    end if

    pawpatch = .true.
    do ityp=1,atoms%astruct%ntypes
       filename = 'psppar.'//atoms%astruct%atomnames(ityp)
       call psp_set_from_dict(dict // filename, &
            & atoms%nzatom(ityp), atoms%nelpsp(ityp), atoms%npspcode(ityp), &
            & atoms%ixcpsp(ityp), atoms%iradii_source(ityp), atoms%psppar(:,:,ityp), radii_cf)
       !To eliminate the runtime warning due to the copy of the array (TD)
       atoms%radii_cf(ityp,:)=radii_cf(:)

       l = .false.
       if (has_key(dict // filename, "PAW patch")) l = dict // filename // "PAW patch"
       pawpatch = pawpatch .and. l

       ! PAW case.
       if (l .and. atoms%npspcode(ityp) == PSPCODE_PAW) then
          ! Allocate the PAW arrays on the fly.
          if (.not. associated(atoms%pawrad)) then
             allocate(atoms%pawrad(atoms%astruct%ntypes))
             allocate(atoms%pawtab(atoms%astruct%ntypes))
             do ityp2 = 1, atoms%astruct%ntypes
                call pawrad_nullify(atoms%pawrad(ityp2))
                call pawtab_nullify(atoms%pawtab(ityp2))
             end do
          end if
          ! Re-read the pseudo for PAW arrays.
          fpaw = dict // filename // SOURCE_KEY
          !write(*,*) 'Reading of PAW atomic-data, under development', trim(fpaw)
          call paw_from_file(atoms%pawrad(ityp), atoms%pawtab(ityp), trim(fpaw), &
               & atoms%nzatom(ityp), atoms%nelpsp(ityp), atoms%ixcpsp(ityp))
       end if
    end do
    call nlcc_set_from_dict(dict, atoms)

    !For PAW psp
    if (pawpatch.and. any(atoms%npspcode /= PSPCODE_PAW)) then
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


  !> Set the value for atoms_data from the dictionary
  subroutine psp_set_from_dict(dict, &
       & nzatom, nelpsp, npspcode, ixcpsp, iradii_source, psppar, radii_cf)
    use module_defs, only: gp, UNINITIALIZED
    use module_atoms
    use dictionaries
    use psp_projectors, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
    implicit none
    !Arguments
    type(dictionary), pointer :: dict
    integer, intent(out) :: nzatom, nelpsp, npspcode, ixcpsp, iradii_source
    real(gp), dimension(0:4,0:6), intent(out) :: psppar
    !Local variables
    type(dictionary), pointer :: loc
    character(len = max_field_length) :: str
    real(gp), dimension(3) :: radii_cf
    integer :: i, l,nlen

    nzatom = -1
    radii_cf(:) = UNINITIALIZED(1._gp)
    psppar(:,:) = 0._gp

    ! We set nzatom at the end as a flag that the psp data are complete.
    if (.not. has_key(dict, ELECTRON_NUMBER)) return
    nelpsp = dict // ELECTRON_NUMBER
    if (.not. has_key(dict, PSPXC_KEY)) return
    ixcpsp = dict // PSPXC_KEY
    ! Local terms
    if (.not. has_key(dict, LPSP_KEY)) return
    loc => dict // LPSP_KEY
    if (.not. has_key(loc, "Rloc")) return
    psppar(0,0) = loc // 'Rloc'
    if (.not. has_key(loc, "Coefficients (c1 .. c4)")) return
    psppar(0,1:4) = loc // 'Coefficients (c1 .. c4)'
    !psppar(0,1) = loc // 'Coefficients (c1 .. c4)' // 0
    !psppar(0,2) = loc // 'Coefficients (c1 .. c4)' // 1
    !psppar(0,3) = loc // 'Coefficients (c1 .. c4)' // 2
    !psppar(0,4) = loc // 'Coefficients (c1 .. c4)' // 3

    ! Nonlocal terms
    if (has_key(dict, NLPSP_KEY)) then
       nlen=dict_len(dict // NLPSP_KEY)
       do i = 1, nlen
          loc => dict // NLPSP_KEY // (i - 1)
          if (.not. has_key(loc, "Channel (l)")) return
          l = loc // "Channel (l)"
          l = l + 1
          if (.not. has_key(loc, "Rloc")) return
          psppar(l,0) = loc // 'Rloc'
          if (.not. has_key(loc, "h_ij terms")) return
          psppar(l,1:6) = loc // 'h_ij terms'
          !psppar(l,1) = loc // 'h_ij terms' // 0
          !psppar(l,2) = loc // 'h_ij terms' // 1
          !psppar(l,3) = loc // 'h_ij terms' // 2
          !psppar(l,4) = loc // 'h_ij terms' // 3
          !psppar(l,5) = loc // 'h_ij terms' // 4
          !psppar(l,6) = loc // 'h_ij terms' // 5
       end do
    end if
    ! Type
    if (.not. has_key(dict, PSP_TYPE)) return
    str = dict // PSP_TYPE
    select case(trim(str))
    case("GTH")
       npspcode = PSPCODE_GTH
    case("HGH")
       npspcode = PSPCODE_HGH
    case("HGH-K")
       npspcode = PSPCODE_HGH_K
    case("HGH-K + NLCC")
       npspcode = PSPCODE_HGH_K_NLCC
    case("PAW")
       npspcode = PSPCODE_PAW
    case default
       return
    end select
    if (npspcode == PSPCODE_HGH_K_NLCC) then
       if (.not. has_key(dict, 'Non Linear Core Correction term')) return
       loc => dict // 'Non Linear Core Correction term'
       if (.not. has_key(loc, "Rcore")) return
       if (.not. has_key(loc, "Core charge")) return
    end if
    ! Valid pseudo, we set nzatom
    if (.not. has_key(dict, ATOMIC_NUMBER)) return
    nzatom = dict // ATOMIC_NUMBER

    ! Optional values.
    if (has_key(dict, RADII_KEY)) then
       loc => dict // RADII_KEY
       if (has_key(loc, COARSE)) radii_cf(1) =  loc // COARSE
       if (has_key(loc, FINE)) radii_cf(2) =  loc // FINE
       if (has_key(loc, COARSE_PSP)) radii_cf(3) =  loc // COARSE_PSP
       
       if (has_key(loc, SOURCE_KEY)) then
          ! Source of the radii
          str = loc // SOURCE_KEY
          select case(str)
          case(RADII_SOURCE(RADII_SOURCE_HARD_CODED))
             iradii_source = RADII_SOURCE_HARD_CODED
          case(RADII_SOURCE(RADII_SOURCE_FILE))
             iradii_source = RADII_SOURCE_FILE
          case(RADII_SOURCE(RADII_SOURCE_USER))
             iradii_source = RADII_SOURCE_USER
          case default
             !Undefined: we assume this the name of a file
             iradii_source = RADII_SOURCE_UNKNOWN
          end select
       end if

    else
       iradii_source = RADII_SOURCE_HARD_CODED
    end if

  end subroutine psp_set_from_dict


  !> Merge all psp data (coming from a file) in the dictionary
  subroutine psp_data_merge_to_dict(dict, nzatom, nelpsp, npspcode, ixcpsp, &
       & psppar, radii_cf, rcore, qcore)
    use module_defs, only: gp, UNINITIALIZED
    use psp_projectors, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
    use module_atoms, only: RADII_SOURCE_FILE
    use dictionaries
    use yaml_strings
    implicit none
    !Arguments
    type(dictionary), pointer :: dict
    integer, intent(in) :: nzatom, nelpsp, npspcode, ixcpsp
    real(gp), dimension(0:4,0:6), intent(in) :: psppar
    real(gp), dimension(3), intent(in) :: radii_cf
    real(gp), intent(in) :: rcore, qcore
    !Local variables
    type(dictionary), pointer :: channel, radii
    integer :: l, i

    ! Type
    select case(npspcode)
    case(PSPCODE_GTH)
       call set(dict // PSP_TYPE, 'GTH')
    case(PSPCODE_HGH)
       call set(dict // PSP_TYPE, 'HGH')
    case(PSPCODE_HGH_K)
       call set(dict // PSP_TYPE, 'HGH-K')
    case(PSPCODE_HGH_K_NLCC)
       call set(dict // PSP_TYPE, 'HGH-K + NLCC')
    case(PSPCODE_PAW)
       call set(dict // PSP_TYPE, 'PAW')
    end select

    call set(dict // ATOMIC_NUMBER, nzatom)
    call set(dict // ELECTRON_NUMBER, nelpsp)
    call set(dict // PSPXC_KEY, ixcpsp)

    ! Local terms
    if (psppar(0,0)/=0) then
       call dict_init(channel)
       call set(channel // 'Rloc', psppar(0,0))
       do i = 1, 4, 1
          call add(channel // 'Coefficients (c1 .. c4)', psppar(0,i))
       end do
       call set(dict // LPSP_KEY, channel)
    end if

    ! nlcc term
    if (npspcode == PSPCODE_HGH_K_NLCC) then
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
       if (radii_cf(1) /= UNINITIALIZED(1._gp)) call set(radii // COARSE, radii_cf(1))
       if (radii_cf(2) /= UNINITIALIZED(1._gp)) call set(radii // FINE, radii_cf(2))
       if (radii_cf(3) /= UNINITIALIZED(1._gp)) call set(radii // COARSE_PSP, radii_cf(3))
       call set(radii // SOURCE_KEY, RADII_SOURCE_FILE)
       call set(dict // RADII_KEY, radii)
    end if

  end subroutine psp_data_merge_to_dict


  !> Read old psppar file (check if not already in the dictionary) and merge to dict
  subroutine atoms_file_merge_to_dict(dict)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use yaml_output, only: yaml_warning
    use public_keys, only: POSINP,SOURCE_KEY
    implicit none
    type(dictionary), pointer :: dict

    type(dictionary), pointer :: types
    character(len = max_field_length) :: str
    integer :: iat, stypes
    character(len=max_field_length), dimension(:), allocatable :: keys
    character(len=27) :: key
    logical :: exists

    ! Loop on types for atomic data.
    call astruct_dict_get_types(dict //POSINP, types)
    if ( .not. associated(types)) return
    allocate(keys(dict_size(types)))
    keys = dict_keys(types)
    stypes = dict_size(types)
    do iat = 1, stypes, 1
       key = 'psppar.' // trim(keys(iat))

       exists = has_key(dict, key)
       if (exists) then
          if (has_key(dict // key, SOURCE_KEY)) then
             str = dict_value(dict // key // SOURCE_KEY)
          else
             str = dict_value(dict // key)
          end if
          if (trim(str) /= "" .and. trim(str) /= TYPE_DICT) then
             !Read the PSP file and merge to dict
             if (trim(str) /= TYPE_LIST) then
                call psp_file_merge_to_dict(dict, key, filename = trim(str))
             else
                call psp_file_merge_to_dict(dict, key, lstring = dict // key)
             end if
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


  !> Read psp file and merge to dict
  subroutine psp_file_merge_to_dict(dict, key, filename, lstring)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_strings
    use f_utils
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    character(len = *), optional, intent(in) :: filename
    type(dictionary), pointer, optional :: lstring
    !Local variables
    integer :: nzatom, nelpsp, npspcode, ixcpsp
    real(gp) :: psppar(0:4,0:6), radii_cf(3), rcore, qcore
    logical :: exists, donlcc, pawpatch
    type(io_stream) :: ios

    if (present(filename)) then
       inquire(file=trim(filename),exist=exists)
       if (.not. exists) return
       call f_iostream_from_file(ios, filename)
    else if (present(lstring)) then
       call f_iostream_from_lstring(ios, lstring)
    else
       call f_err_throw("Error in psp_file_merge_to_dict, either 'filename' or 'lstring' should be present.", &
            & err_name='BIGDFT_RUNTIME_ERROR')
    end if
    !ALEX: if npspcode==PSPCODE_HGH_K_NLCC, nlccpar are read from psppar.Xy via rcore and qcore 
    call psp_from_stream(ios, nzatom, nelpsp, npspcode, ixcpsp, &
         & psppar, donlcc, rcore, qcore, radii_cf, pawpatch)
    call f_iostream_release(ios)

    if (has_key(dict, key)) call dict_remove(dict, key)
    call psp_data_merge_to_dict(dict // key, nzatom, nelpsp, npspcode, ixcpsp, &
         & psppar, radii_cf, rcore, qcore)
    call set(dict // key // "PAW patch", pawpatch)
    if (present(filename)) then
       call set(dict // key // SOURCE_KEY, filename)
    else
       call set(dict // key // SOURCE_KEY, "In-line")
    end if
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
  
  subroutine astruct_dict_get_types(dict, types)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict, types

    type(dictionary), pointer :: atoms, at
    character(len = max_field_length) :: str
    integer :: ityp

    if (ASTRUCT_POSITIONS .notin. dict) then
       nullify(types)
       return
    end if
    call dict_init(types)
    ityp = 0
    atoms => dict_iter(dict // ASTRUCT_POSITIONS)
    do while(associated(atoms))
       at => dict_iter(atoms)
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
       atoms => dict_next(atoms)
    end do
  end subroutine astruct_dict_get_types


  !> Read Atomic positions and merge into dict
  subroutine astruct_file_merge_to_dict(dict, key, filename)
    use module_base, only: gp, UNINITIALIZED, bigdft_mpi,f_routine,f_release_routine, &
        & BIGDFT_INPUT_FILE_ERROR,f_free_ptr
    use module_atoms, only: set_astruct_from_file,atomic_structure,&
         nullify_atomic_structure,deallocate_atomic_structure,astruct_merge_to_dict
    use public_keys, only: POSINP,RADICAL_NAME
    use dictionaries
    use yaml_strings
    implicit none
    !Arguments
    type(dictionary), pointer :: dict          !< Contains (out) all the information
    character(len = *), intent(in) :: key      !< Key of the dictionary where it should be have the information
    character(len = *), intent(in) :: filename !< Name of the filename where the astruct should be read
    !Local variables
    type(atomic_structure) :: astruct
    !type(DFT_global_output) :: outs
    character(len=max_field_length) :: msg,radical
    integer :: ierr,iat
    real(gp) :: energy
    real(gp), dimension(:,:), pointer :: fxyz
    type(dictionary), pointer :: dict_tmp,pos


    call f_routine(id='astruct_file_merge_to_dict')
    ! Read atomic file, old way
    call nullify_atomic_structure(astruct)
    !call nullify_global_output(outs)
    !Try to read the atomic coordinates from files
    call f_err_open_try()
    nullify(fxyz)
    call set_astruct_from_file(filename, bigdft_mpi%iproc, astruct, &
         energy = energy, fxyz = fxyz)
    !print *,'test2',associated(fxyz)
    !Check if BIGDFT_INPUT_FILE_ERROR
    ierr = f_get_last_error(msg) 
    call f_err_close_try()
    if (ierr == 0) then
       dict_tmp => dict // key
       !No errors: we have all information in astruct and put into dict

       call astruct_merge_to_dict(dict_tmp, astruct, astruct%rxyz)

       call set(dict_tmp // ASTRUCT_PROPERTIES // POSINP_SOURCE, filename)

       if (GOUT_FORCES .in. dict_tmp) call dict_remove(dict_tmp, GOUT_FORCES)
       if (associated(fxyz)) then
          pos => dict_tmp // GOUT_FORCES
          do iat=1,astruct%nat
             call add(pos, dict_new(astruct%atomnames(astruct%iatype(iat)) .is. fxyz(:,iat)))
          end do
       end if

       if (GOUT_ENERGY .in. dict_tmp) call dict_remove(dict_tmp, GOUT_ENERGY)
       if (energy /= UNINITIALIZED(energy)) call set(dict_tmp // GOUT_ENERGY, energy)
       !call global_output_merge_to_dict(dict // key, outs, astruct)
       call deallocate_atomic_structure(astruct)

    else if (ierr == BIGDFT_INPUT_FILE_ERROR) then
       !Found no file: maybe already inside the yaml file ?
       !Check if posinp is in dict
       if ( POSINP .notin.  dict) then
          ! Raise an error
          call f_strcpy(src='input',dest=radical)
          !modify the radical name if it exists
          radical = dict .get. RADICAL_NAME
          msg = "No section 'posinp' for the atomic positions in the file '"//&
               trim(radical) // ".yaml'. " // trim(msg)
          call f_err_throw(err_msg=msg,err_id=ierr)
       end if
    else 
       ! Raise an error
       call f_err_throw(err_msg=msg,err_id=ierr)
    end if
    call f_free_ptr(fxyz)
    !call deallocate_global_output(outs)
    call f_release_routine()

  end subroutine astruct_file_merge_to_dict


  !> Allocate the astruct variable from the dictionary of input data
  !! retrieve also other information like the energy and the forces if requested
  !! and presend in the dictionary
  subroutine astruct_set_from_dict(dict, astruct, comment)
    use module_defs, only: gp, Bohr_Ang, UNINITIALIZED
    use module_atoms, only: atomic_structure, nullify_atomic_structure,frozen_ftoi
    use dictionaries
    use dynamic_memory
    implicit none
    !Arguments
    type(dictionary), pointer :: dict !< dictionary of the input variables
                                      !! the keys have to be declared like input_dicts module
    type(atomic_structure), intent(out) :: astruct          !< Structure created from the file
    character(len = 1024), intent(out), optional :: comment !< Extra comment retrieved from the file if present
    !local variables
    character(len=*), parameter :: subname='astruct_set_from_dict'
    type(dictionary), pointer :: pos, at, atData, types
    character(len = max_field_length) :: str
    integer :: iat, ityp, units, igspin, igchrg, nsgn, ntyp, ierr

    call f_routine(id='astruct_set_from_dict')

    call nullify_atomic_structure(astruct)
    astruct%nat = -1
    if (present(comment)) write(comment, "(A)") " "

    !if (.not. has_key(dict, ASTRUCT_POSITIONS)) return

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
    ntyp = max(dict_size(types),0) !if types is nullified ntyp=-1
    call astruct_set_n_types(astruct, ntyp)
    ! astruct%atomnames = dict_keys(types)
    ityp = 1
    at => dict_iter(types)
    do while (associated(at))
       astruct%atomnames(ityp) = trim(dict_key(at))
       ityp = ityp + 1
       at => dict_next(at)
    end do
    ! The atoms
    if (ASTRUCT_POSITIONS .in. dict) then
       pos => dict // ASTRUCT_POSITIONS
       call astruct_set_n_atoms(astruct, dict_len(pos))
       at => dict_iter(pos)
       iat = 0
       do while(associated(at))
          iat = iat + 1

          igspin = 0
          igchrg = 0
          nsgn   = 1
          !at => pos // (iat - 1)
          atData => dict_iter(at)!at%child
          do while(associated(atData))
             str = dict_key(atData)
             if (trim(str) == "Frozen") then
                str = dict_value(atData)
                call frozen_ftoi(str(1:4), astruct%ifrztyp(iat),ierr)
             else if (trim(str) == "IGSpin") then
                igspin = atData
             else if (trim(str) == "IGChg") then
                igchrg = atData
                if (igchrg >= 0) then
                   nsgn = 1
                else
                   nsgn = -1
                end if
             else if (trim(str) == "int_ref_atoms_1") then
                astruct%ixyz_int(1,iat) = atData
             else if (trim(str) == "int_ref_atoms_2") then
                astruct%ixyz_int(2,iat) = atData
             else if (trim(str) == "int_ref_atoms_3") then
                astruct%ixyz_int(3,iat) = atData
             else if (dict_len(atData) == 3) then
                astruct%iatype(iat) = types // dict_key(atData)
                astruct%rxyz(:, iat) = atData
             end if
             atData => dict_next(atData)
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
          at => dict_next(at)
       end do
    else
       call astruct_set_n_atoms(astruct,0)
    end if

    if (has_key(dict, ASTRUCT_PROPERTIES)) then
       pos => dict // ASTRUCT_PROPERTIES
       if (has_key(pos, "info") .and. present(comment)) comment = pos // "info"
       if (has_key(pos, "format")) astruct%inputfile_format = pos // "format"
    end if

    call dict_free(types)

    call f_release_routine()

  end subroutine astruct_set_from_dict


  !subroutine aocc_to_dict(dict, nspin, noncoll, nstart, aocc, nelecmax, lmax, nsccode)
  !  use module_defs, only: gp
  !  use dictionaries
  !  implicit none
  !  integer, intent(in) :: nelecmax, lmax, nsccode, nspin, noncoll, nstart
  !  type(dictionary), pointer :: dict
  !  real(gp), dimension(nelecmax), intent(in) :: aocc

  !  type(dictionary), pointer :: val
  !  character(len = 4) :: key
  !  integer :: l, inl, nl, iocc, sccode, nsc, i
  !  character(len = 1), dimension(4), parameter :: lname = (/ "s", "p", "d", "f" /)

  !  call dict_init(dict)

  !  sccode = nsccode
  !  iocc=0
  !  do l = 1, lmax
  !     iocc=iocc+1
  !     ! Get number of shells for this channel 
  !     !(to be corrected, the rule is not the same)
  !     nl = int(aocc(iocc))
  !     ! Get number of semi cores for this channel
  !     nsc = modulo(sccode, 4)
  !     sccode = sccode / 4
  !     if (nl == 0) cycle
  !     do inl = 1, nl, 1
  !        if (inl <= nsc) then
  !           write(key, "(A1,I1,A1,A1)") "(", nstart + inl, lname(l), ")"
  !        else
  !           write(key, "(I1, A1)") nstart + inl, lname(l)
  !        end if
  !        call dict_init(val)
  !        do i = 1, nspin * noncoll * (2 * l - 1), 1
  !           iocc=iocc+1
  !           call add(val, aocc(iocc))
  !        end do
  !        call set(dict // key, val)
  !     end do
  !  end do
  !end subroutine aocc_to_dict


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

    call f_routine(id='occupation_set_from_dict')

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

    call f_release_routine()

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
!!$

  subroutine occupation_data_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key

    logical :: exists
    integer :: ierror, ntu, ntd, nt, i, iorb, lline, lstring
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
       lline = len(line)
       do i=1,lline
          if (line(i:i) == '/') then
             line(i:i) = ':'
          end if
       end do
       read(line,*,iostat=ierror) iorb,string
       if (ierror /= 0) then
          exit
       end if
       !Transform back the ':' into '/'
       lstring = len(string)
       do i=1,lstring
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

    call set(dict // key // SOURCE_KEY, filename)

  end subroutine occupation_data_file_merge_to_dict

  subroutine neb_set_from_dict(dict, opt, climbing_, imax, nimg_, &
       & cv, tol, ds_, kmin, kmax, temp_, damp_, meth)
    use module_defs, only: gp
    use dictionaries
    !use module_input_keys
    use public_keys
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(out) :: opt, climbing_
    integer, intent(out) :: imax, nimg_
    real(kind = gp), intent(out) :: cv, tol, ds_, damp_, kmin, kmax, temp_
    character(len = max_field_length), intent(out) :: meth

    if (.not. has_key(dict, GEOPT_VARIABLES)) return
    if (trim(dict_value(dict // GEOPT_VARIABLES // GEOPT_METHOD)) /= "NEB") return

    opt       = dict // GEOPT_VARIABLES // EXTREMA_OPT
    climbing_ = dict // GEOPT_VARIABLES // NEB_CLIMBING
    imax      = dict // GEOPT_VARIABLES // NCOUNT_CLUSTER_X
    nimg_     = dict // GEOPT_VARIABLES // NIMG
    cv        = dict // GEOPT_VARIABLES // FORCEMAX
    tol       = dict // GEOPT_VARIABLES // FIX_TOL
    ds_       = dict // GEOPT_VARIABLES // BETAX
    kmin      = dict // GEOPT_VARIABLES // SPRINGS_K // 0
    kmax      = dict // GEOPT_VARIABLES // SPRINGS_K // 1
    meth      = dict // GEOPT_VARIABLES // NEB_METHOD
    if (has_key(dict // GEOPT_VARIABLES, TEMP)) temp_ = dict // GEOPT_VARIABLES // TEMP
    if (has_key(dict // GEOPT_VARIABLES, NEB_DAMP)) damp_ = dict // GEOPT_VARIABLES // NEB_DAMP
  end subroutine neb_set_from_dict
end module module_input_dicts

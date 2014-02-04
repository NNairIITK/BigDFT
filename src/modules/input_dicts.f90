!> @file
!>  Modules which contains all the interfaces to parse input dictionary.
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>  Modules which contains all interfaces to parse input dictionary.
module module_input_dicts

  implicit none

  private

  ! Main creation routine
  public :: user_dict_from_files

  ! Dictionary completion
  public :: psp_dict_fill_all, psp_dict_analyse

  ! Types from dictionaries
  public :: psp_set_from_dict, nlcc_set_from_dict
  public :: astruct_set_from_dict
  public :: atomic_data_set_from_dict

  ! Types to dictionaries
  public :: psp_data_merge_to_dict
  public :: astruct_merge_to_dict

  ! Dictionaries from files (old formats).
  public :: psp_file_merge_to_dict, nlcc_file_merge_to_dict
  public :: atoms_file_merge_to_dict
  public :: astruct_file_merge_to_dict
  public :: atomic_data_file_merge_to_dict

contains

  subroutine user_dict_from_files(dict,radical,posinp, mpi_env)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use module_defs, only: mpi_environment
    use module_interfaces, only: read_input_dict_from_files
    implicit none
    !Arguments
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: radical, posinp
    type(mpi_environment), intent(in) :: mpi_env
    !Local variables
    character(len = max_field_length) :: str

    nullify(dict)
    !read the input file(s) and transform them into a dictionary
    call read_input_dict_from_files(trim(radical), mpi_env, dict)

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

    if (.not. has_key(dict, "Atomic occupation")) then
       ! Add old input.occup
       call atomic_data_file_merge_to_dict(dict, "Atomic occupation", &
            & trim(radical) // ".occup")
    else
       str = dict_value(dict // "Atomic occupation")
       if (trim(str) /= TYPE_DICT .and. trim(str) /= TYPE_LIST .and. trim(str) /= "") then
          call atomic_data_file_merge_to_dict(dict, "Atomic occupation", trim(str))
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
       call allocate_atoms_nat(atoms, "psp_dict_analyse")
       call allocate_atoms_ntypes(atoms, "psp_dict_analyse")
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
!                atoms%nlccpar(0, nlcc_dim + i) = coeffs // 0
!                atoms%nlccpar(1, nlcc_dim + i) = coeffs // 1
!                atoms%nlccpar(2, nlcc_dim + i) = coeffs // 2
!                atoms%nlccpar(3, nlcc_dim + i) = coeffs // 3
!                atoms%nlccpar(4, nlcc_dim + i) = coeffs // 4
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
       if (radii_cf(1) /= UNINITIALIZED(1._gp)) call set(radii // "Fine", radii_cf(2))
       if (radii_cf(1) /= UNINITIALIZED(1._gp)) call set(radii // "Coarse PSP", radii_cf(3))
       call set(dict // "Radii of active regions (AU)", radii)
    end if
  end subroutine psp_data_merge_to_dict

  subroutine atoms_file_merge_to_dict(dict)
    use dictionaries
    use dictionaries_base, only: TYPE_DICT, TYPE_LIST
    use yaml_output, only: yaml_warning
    implicit none
    type(dictionary), pointer :: dict

    type(dictionary), pointer :: types, at
    character(len = max_field_length) :: str
    integer :: iat
    character(max_field_length), dimension(:), allocatable :: keys
    character(len=27) :: key
    logical :: exists

    ! Loop on types for atomic data.
    call dict_init(types)
    do iat = 0, dict_len(dict // "posinp" // "Positions") - 1, 1
       at => dict_iter(dict // "posinp" // "Positions" // iat)
       do while(associated(at))
          str = dict_key(at)
          if (dict_len(at) == 3 .and. .not. has_key(types, str)) call set(types // str, ".")
          at => dict_next(at)
       end do
    end do
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
  subroutine astruct_merge_to_dict(dict, astruct, rxyz, fxyz, energy, comment)
    use module_types, only: atomic_structure
    use module_defs, only: gp, UNINITIALIZED, Bohr_Ang
    use dictionaries
    use yaml_strings
    implicit none
    type(dictionary), pointer :: dict
    type(atomic_structure), intent(in) :: astruct
    real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
    real(gp), dimension(3, astruct%nat), intent(in), optional :: fxyz
    real(gp), intent(in), optional :: energy
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
       call set(dict // "Units", 'angstroem')
       factor=Bohr_Ang
    case('reduced')
       call set(dict // "Units", 'reduced')
       reduced = .true.
    case('atomic','atomicd0','bohr','bohrd0')
       ! Default, store nothing
    end select Units

    !cell information
    BC :select case(astruct%geocode)
    case('S')
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(1)*factor(1)))
       call add(dict // "Cell", '.inf')
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(3)*factor(3)))
       !angdeg to be added
       if (reduced) then
          factor(1) = 1._gp / astruct%cell_dim(1)
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('W')
       call add(dict // "Cell", '.inf')
       call add(dict // "Cell", '.inf')
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(3)*factor(3)))
       if (reduced) then
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('P')
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(1)*factor(1)))
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(2)*factor(2)))
       call add(dict // "Cell", yaml_toa(astruct%cell_dim(3)*factor(3)))
       !angdeg to be added
       if (reduced) then
          factor(1) = 1._gp / astruct%cell_dim(1)
          factor(2) = 1._gp / astruct%cell_dim(2)
          factor(3) = 1._gp / astruct%cell_dim(3)
       end if
    case('F')
       ! Default, store nothing
    end select BC

    pos => dict // "Positions"
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

    if (present(fxyz)) then
       pos => dict // "Forces (Ha/Bohr)"
       do iat=1,astruct%nat
          call add(pos, yaml_toa(fxyz(:, iat)))
       end do
    end if

    if (present(energy)) then
       if (energy /= UNINITIALIZED(energy)) &
            & call add(dict // "Properties" // "Energy (Ha)", energy)
    end if

    if (present(comment)) then
       if (len_trim(comment) > 0) &
            & call add(dict // "Properties" // "Info", comment)
    end if

    if (len_trim(astruct%inputfile_format) > 0) &
         & call set(dict // "Properties" // "Format", astruct%inputfile_format)
  end subroutine astruct_merge_to_dict


  subroutine astruct_set_from_dict(dict, astruct, comment, energy, fxyz)
    use module_types, only: atomic_structure
    use module_defs, only: gp, Bohr_Ang, UNINITIALIZED
    use dictionaries
    use dynamic_memory
    implicit none
    type(dictionary), pointer :: dict
    type(atomic_structure), intent(out) :: astruct
    real(gp), intent(out), optional :: energy
    real(gp), dimension(:,:), pointer, optional :: fxyz
    character(len = 1024), intent(out), optional :: comment

    !local variables
    character(len=*), parameter :: subname='read_dict_positions'
    type(dictionary), pointer :: pos, at
    character(len = max_field_length) :: str
    integer :: iat, ityp, units, igspin, igchrg, nsgn, ntyp,ierr
    character(len=20), dimension(100) :: atomnames

    call astruct_nullify(astruct)
    astruct%nat = -1
    if (present(energy)) energy = UNINITIALIZED(energy)
    if (present(comment)) write(comment, "(A)") " "
    if(present(fxyz)) nullify(fxyz)
    if (.not. has_key(dict, "Positions")) return

    ! The units
    units = 0
    write(astruct%units, "(A)") "bohr"
    if (has_key(dict, "Units")) astruct%units = dict // "Units"
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
    if (.not. has_key(dict, "Cell")) then
       astruct%geocode = 'F'
    else
       astruct%geocode = 'P'
       ! z
       astruct%cell_dim(3) = dict // "Cell" // 2
       ! y
       str = dict // "Cell" // 1
       if (trim(str) == ".inf") then
          astruct%geocode = 'S'
       else
          astruct%cell_dim(2) = dict // "Cell" // 1
       end if
       ! x
       str = dict // "Cell" // 0
       if (trim(str) == ".inf") then
          astruct%geocode = 'W'
       else
          astruct%cell_dim(1) = dict // "Cell" // 0
       end if
    end if
    if (units == 1) astruct%cell_dim = astruct%cell_dim / Bohr_Ang
    ! The atoms
    if (.not. has_key(dict, "Positions")) return
    pos => dict // "Positions"
    astruct%nat = dict_len(pos)
    call astruct_set_n_atoms(astruct, astruct%nat, subname)
    ntyp = 0
    do iat = 1, astruct%nat
       at => pos // (iat - 1)
       igspin = 0
       igchrg = 0
       nsgn   = 1
       at => at%child
       do while(associated(at))
          str = dict_key(at)
          if (trim(str) == "Frozen") then
             str = dict_value(at)
             call frozen_ftoi(str(1:4), astruct%ifrztyp(iat),ierr)
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
             do ityp=1,ntyp
                if (str(1:20) == atomnames(ityp)) then
                   astruct%iatype(iat)=ityp
                   exit
                endif
             enddo
             if (ityp > ntyp) then
                ntyp=ntyp+1
                if (ntyp > 100) then
                   write(*,*) 'more than 100 atomnames not permitted'
                   astruct%nat = -1
                   return
                end if
                atomnames(ityp)=str(1:20)
                astruct%iatype(iat)=ntyp
             end if
             astruct%rxyz(1, iat) = at // 0
             astruct%rxyz(2, iat) = at // 1
             astruct%rxyz(3, iat) = at // 2
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
          if (astruct%cell_dim(1) > 0.) astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),1.0_gp) * astruct%cell_dim(1)
          if (astruct%cell_dim(2) > 0.) astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),1.0_gp) * astruct%cell_dim(2)
          if (astruct%cell_dim(3) > 0.) astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),1.0_gp) * astruct%cell_dim(3)
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
    if (has_key(dict, "Forces") .and. present(fxyz)) then
       fxyz = f_malloc_ptr((/ 3, astruct%nat /), subname)
       pos => dict // "Forces"
       do iat = 1, astruct%nat
          at => pos // (iat - 1)
          fxyz(1, iat) = at // 0
          fxyz(2, iat) = at // 1
          fxyz(3, iat) = at // 2
       end do
    end if

    call astruct_set_n_types(astruct, ntyp, subname)
    astruct%atomnames(1:ntyp) = atomnames(1:ntyp)

    if (has_key(dict, "Properties")) then
       pos => dict // "Properties"
       if (has_key(pos, "Energy (Ha)") .and. present(energy)) energy = pos // "Energy (Ha)"
       if (has_key(pos, "Info") .and. present(comment)) comment = pos // "Info"
       if (has_key(pos, "Format")) astruct%inputfile_format = pos // "Format"
    end if

  end subroutine astruct_set_from_dict


  subroutine astruct_file_merge_to_dict(dict, key, filename)
    use module_defs, only: gp, UNINITIALIZED, bigdft_mpi
    use module_interfaces, only: read_atomic_file
    use dictionaries
    use yaml_strings
    use module_types, only: atomic_structure
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: filename, key
    
    type(atomic_structure) :: astruct
    integer :: ierr

    ! Read atomic file, old way
    call astruct_nullify(astruct)
    call read_atomic_file(filename, bigdft_mpi%iproc, astruct, status = ierr)
    if (ierr == 0) then
       call astruct_merge_to_dict(dict // key, astruct, astruct%rxyz)
       call set(dict // key // "Properties" // "Source", filename)
       call deallocate_atomic_structure(astruct, "astruct_file_merge_to_dict")
    end if
  end subroutine astruct_file_merge_to_dict

  subroutine aocc_from_dict(dict,nspin,nspinor,nelecmax,lmax,nmax,aocc,nsccode)
    use module_defs, only: gp, UNINITIALIZED
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: nelecmax,lmax,nmax,nspinor,nspin
    integer, intent(out) :: nsccode
    real(gp), dimension(nelecmax), intent(out) :: aocc

    !local variables
    character(len = max_field_length) :: key
    character(max_field_length), dimension(:), allocatable :: keys
    integer :: i, ln
    integer :: m,n,iocc,icoll,inl,noncoll,l,ispin,is,lsc
    integer, dimension(lmax) :: nl,nlsc
    real(gp), dimension(2*(2*lmax-1),nmax,lmax) :: allocc

    nl(:)=0
    nlsc(:)=0
    allocc(:,:,:) = UNINITIALIZED(1._gp)

    !if non-collinear it is like nspin=1 but with the double of orbitals
    if (nspinor == 4) then
       noncoll=2
    else
       noncoll=1
    end if

    allocate(keys(dict_size(dict)))
    keys = dict_keys(dict)
    do i = 1, dict_size(dict), 1
       key = keys(i)
       ln = len_trim(key)
       is = 1
       if (key(1:1) == "(" .and. key(ln:ln) == ")") is = 2
       ! Read the major quantum number
       read(key(is:is), "(I1)") n
       is = is + 1
       ! Read the channel
       select case(key(is:is))
       case('s')
          l=1
       case('p')
          l=2
       case('d')
          l=3
       case('f')
          l=4
       case default
          stop "wrong channel"
       end select
       nl(l) = nl(l) + 1
       if (is == 3) nlsc(l) = nlsc(l) + 1
       if (nlsc(l) > 2) stop 'cannot admit more than two semicore orbitals per channel'

       !read the different atomic occupation numbers
       if (dict_len(dict // key) /= nspin*noncoll*(2*l-1)) then
          write(*,*) "Awaited: ", nspin*noncoll*(2*l-1), nspin, noncoll, l
          write(*,*) "provided", dict_len(dict // key)
          stop 'Not enough aocc'
       end if
       do m = 1, nspin*noncoll*(2*l-1), 1
          allocc(m, n, l) = dict // key // (m - 1)
       end do
    end do
    deallocate(keys)

    !put the values in the aocc array
    aocc(:)=0.0_gp
    iocc=0
    do l=1,lmax
       iocc=iocc+1
       aocc(iocc)=real(nl(l),gp)
       do inl=1,nmax
          if (allocc(1, inl, l) == UNINITIALIZED(1._gp)) cycle
          do ispin=1,nspin
             do m=1,2*l-1
                do icoll=1,noncoll !non-trivial only for nspinor=4
                   iocc=iocc+1
                   aocc(iocc)=allocc(icoll+(m-1)*noncoll+(ispin-1)*(2*l-1)*noncoll,inl,l)
                end do
             end do
          end do
       end do
    end do

    !then calculate the nsccode
    nsccode=0
    do lsc=1,lmax
       nsccode=nsccode+nlsc(lsc) * (4**(lsc-1))
    end do
  end subroutine aocc_from_dict

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
    use ao_inguess
    use module_types, only: atoms_data
    use dictionaries
    use dynamic_memory
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict
    type(atoms_data), intent(inout) :: atoms
    character(len = *), intent(in) :: key
    integer, intent(in) :: nspin

    integer :: iat, ityp, nsccode, mxpl, mxchg, nsp, nspinor
    integer :: ichg, ispol, icoll, iocc, ispin, l, inl, m, nl, noncoll
    integer, parameter :: nelecmax=32,nmax=6,lmax=4
    !character(len=2) :: symbol
    real(gp) :: rcov,elec!,ehomo,rprb
    real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
    character(len = max_field_length) :: at
    real(gp), dimension(nmax,lmax) :: eleconf_

    !control the spin
    select case(nspin)
    case(1)
       nsp=1
       nspinor=1
       noncoll=1
    case(2)
       nsp=2
       nspinor=1
       noncoll=1
    case(4)
       nsp=1
       nspinor=4
       noncoll=2
    case default
       call yaml_warning('nspin not valid. Value=' // trim(yaml_toa(nspin)))
       !write(*,*)' ERROR: nspin not valid:',nspin
       stop
    end select

    do ityp = 1, atoms%astruct%ntypes, 1
       call atomic_info(atoms%nzatom(ityp),atoms%nelpsp(ityp),elconf=neleconf,&
            amu=atoms%amu(ityp),rcov=rcov,nsccode=nsccode,&
            maxpol=mxpl,maxchg=mxchg)
!       call eleconf(atoms%nzatom(ityp), atoms%nelpsp(ityp), symbol,rcov,rprb,ehomo,&
!            neleconf,nsccode,mxpl,mxchg,atoms%amu(ityp))
       !define the localization radius for the Linear input guess
       atoms%rloc(ityp,:) = rcov * 10.0

       do iat = 1, atoms%astruct%nat, 1
          if (atoms%astruct%iatype(iat) /= ityp) cycle

          ! Some checks from input values.
          call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
          if (abs(ispol) > mxpl+abs(ichg)) then
             !if (iproc ==0) 
             write(*,'(1x,a,i0,a,a,2(a,i0))')&
                  'ERROR: Input polarisation of atom No.',iat,&
                  ' (',trim(atoms%astruct%atomnames(ityp)),') must be <=',mxpl,&
                  ', while found ',ispol
             stop 
          end if
          if (abs(ichg) > mxchg) then
             !if (iproc ==0) 
             write(*,'(1x,a,i0,a,a,2(a,i0))')&
                  'ERROR: Input charge of atom No.',iat,&
                  ' (',trim(atoms%astruct%atomnames(ityp)),') must be <=',mxchg,&
                  ', while found ',ichg
             stop
          end if

          ! Fill this atom with default values from eleconf.
          atoms%iasctype(iat)=nsccode
          !correct the electronic configuration in case there is a charge
          !if (ichg /=0) then
          call correct_semicore(nmax,lmax-1,ichg,&
               neleconf,eleconf_,atoms%iasctype(iat))
          !end if
          call at_occnums(ispol,nsp,nspinor,nmax,lmax,nelecmax,&
               eleconf_,atoms%aocc(1:,iat))

          ! Possible overwrite.
          if (has_key(dict, key)) then
             write(at, "(A,I0)") "Atom ", iat
             if (has_key(dict // key, at)) then
                ! Case with an atom specific aocc
                call aocc_from_dict(dict // key // at, &
                     & nsp, nspinor, nelecmax, lmax, nmax, &
                     & atoms%aocc(:,iat), atoms%iasctype(iat))
             else if (has_key(dict // key, trim(atoms%astruct%atomnames(ityp)))) then
                ! Case with a element specific aocc
                call aocc_from_dict(dict // key // trim(atoms%astruct%atomnames(ityp)), &
                     & nsp, nspinor, nelecmax, lmax, nmax, &
                     & atoms%aocc(:,iat), atoms%iasctype(iat))
             end if
          end if

          !check the total number of electrons
          elec=0.0_gp
          iocc=0
          do l=1,lmax
             iocc=iocc+1
             nl=nint(atoms%aocc(iocc,iat))
             do inl=1,nl
                do ispin=1,nsp
                   do m=1,2*l-1
                      do icoll=1,noncoll !non-trivial only for nspinor=4
                         iocc=iocc+1
                         elec=elec+atoms%aocc(iocc,iat)
                      end do
                   end do
                end do
             end do
          end do
          if (nint(elec) /= atoms%nelpsp(ityp) - ichg) then
             call print_eleconf(nsp,nspinor,2,nelecmax,lmax,atoms%aocc(1,iat),atoms%iasctype(iat))
             write(*,*)'ERROR: the total atomic charge ',elec,&
                  ' is different from the PSP charge ',atoms%nelpsp(ityp),&
                  ' plus the charge ',-ichg
             stop
          end if
       end do
    end do

    atoms%natsc = 0
    do iat=1,atoms%astruct%nat
       if (atoms%iasctype(iat) /= 0) atoms%natsc=atoms%natsc+1
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
    if (ierror /= 0) then
       call yaml_warning('Failed to open the existing file '// trim(filename))
       stop
    end if

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

end module module_input_dicts

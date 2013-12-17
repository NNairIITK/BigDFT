!>  Modules which contains all interfaces to parse input dictionary.
module module_input_dicts

  implicit none

  private

  ! Dictionary completion
  public :: psp_dict_fill_all, psp_dict_analyse

  ! Types from dictionaries
  public :: psp_set_from_dict, nlcc_set_from_dict
  public :: astruct_set_from_dict

  ! Dictionaries to types
  public :: psp_data_merge_to_dict
  public :: astruct_merge_to_dict

  ! Dictionaries from files (old formats).
  public :: psp_file_merge_to_dict, nlcc_file_merge_to_dict
  public :: atoms_file_merge_to_dict
  public :: astruct_file_merge_to_dict

contains

  subroutine psp_dict_fill_all(dict, atomname, run_ixc)
    use module_defs, only: gp, UNINITIALIZED, bigdft_mpi
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
    integer, parameter :: nmax=6,lmax=4
    !integer, parameter :: nelecmax=32
    character(len=2) :: symbol
    integer :: i,mxpl,mxchg,nsccode
    real(gp) :: rcov,rprb,ehomo,radfine,amu,rad
    real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
    type(dictionary), pointer :: radii
    real(gp), dimension(3) :: radii_cf

    filename = 'psppar.' // atomname

    radii_cf = UNINITIALIZED(1._gp)
    if (has_key(dict // filename, "Radii of active regions (AU)")) then
       radii => dict // filename // "Radii of active regions (AU)"
       if (has_key(radii, "Coarse")) radii_cf(1) =  radii // "Coarse"
       if (has_key(radii, "Fine")) radii_cf(2) =  radii // "Fine"
       if (has_key(radii, "Coarse PSP")) radii_cf(3) =  radii // "Coarse PSP"
    end if

    exists = has_key(dict // filename, "NonLocal PSP Parameters")
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

    if (any(radii_cf == UNINITIALIZED(1.0_gp))) then
       !see whether the atom is semicore or not
       !and consider the ground state electronic configuration
       call eleconf(nzatom, nelpsp,symbol,rcov,rprb,ehomo,&
            neleconf,nsccode,mxpl,mxchg,amu)

       !assigning the radii by calculating physical parameters
       radii_cf(1)=1._gp/sqrt(abs(2._gp*ehomo))
       radfine = dict // filename // "Local Pseudo Potential (HGH convention)" // "Rloc"
       do i=1, dict_len(dict // filename // "NonLocal PSP Parameters")
          rad = dict // filename // "NonLocal PSP Parameters" // (i - 1) // "Rloc"
          if (rad /= 0._gp) then
             radfine=min(radfine, rad)
          end if
       end do
       radii_cf(2)=radfine
       radii_cf(3)=radfine
       radii => dict // filename // "Radii of active regions (AU)"
       call set(radii // "Coarse", radii_cf(1))
       call set(radii // "Fine", radii_cf(2))
       call set(radii // "Coarse PSP", radii_cf(3))
       call set(radii // "Source", "Hard-coded")
    else
       call set(radii // "Source", "User-defined")
    end if
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

       l = dict // filename // "PAW patch"
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
       filename = 'psppar.' // trim(atoms%astruct%atomnames(ityp))
       if (.not. has_key(dict, filename)) continue
       if (.not. has_key(dict // filename, 'Non Linear Core Correction term')) continue
       nloc => dict // filename // 'Non Linear Core Correction term'
       if (has_key(nloc, "Valence") .and. has_key(nloc, "Conduction")) then
          n = dict_len(nloc // "Valence")
          nlcc_dim = nlcc_dim + n
          atoms%nlcc_ngv(ityp) = int((sqrt(real(1 + 8 * n)) - 1) / 2)
          n = dict_len(nloc // "Conduction")
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
          !ALEX: These are preferably read from psppar.Xy, as stored in the
          !local variables rcore and qcore
          nloc => dict // filename // 'Non Linear Core Correction term'
          if (has_key(nloc, "Valence") .and. has_key(nloc, "Conduction")) then
             n = dict_len(nloc // "Valence")
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
             n = dict_len(nloc // "Conduction")
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
    if (.not. has_key(dict, "NonLocal PSP Parameters")) return
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
       if (any(psppar(l,1:3) /= 0._gp)) then
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
    implicit none
    type(dictionary), pointer :: dict

    type(dictionary), pointer :: types, at
    character(len = max_field_length) :: str
    integer :: iat
    character(max_field_length), dimension(:), allocatable :: keys
    character(len=27) :: filename
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
       filename = 'psppar.' // trim(keys(iat))

       exists = has_key(dict, filename)
       if (exists) exists = has_key(dict // filename, 'Pseudopotential XC')
       if (.not. exists) call psp_file_merge_to_dict(dict, filename, filename)

       exists = has_key(dict, filename)
       if (exists) exists = has_key(dict // filename, 'Non Linear Core Correction term')
       if (.not. exists) call nlcc_file_merge_to_dict(dict, filename, 'nlcc.' // trim(keys(iat)))
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
    call dict_init(gauss)
    do ig=1,(ngv*(ngv+1))/2
       read(79,*) coeffs
       do i = 0, 4, 1
          call add(gauss, coeffs(i))
       end do
    end do
    call set(psp // 'Non Linear Core Correction term' // "Valence", gauss)

    read(79,*)ngc
    call dict_init(gauss)
    do ig=1,(ngc*(ngc+1))/2
       read(79,*) coeffs
       do i = 0, 4, 1
          call add(gauss, coeffs(i))
       end do
    end do
    call set(psp // 'Non Linear Core Correction term' // "Conduction", gauss)

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
          call set(at // "Frozen", frzstr)
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
    integer :: iat, ityp, units, igspin, igchrg, nsgn, ntyp
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
       call set(dict // key // "Properties" // "Format", astruct%inputfile_format)
       call deallocate_atomic_structure(astruct, "astruct_file_merge_to_dict")
    end if
  end subroutine astruct_file_merge_to_dict

end module module_input_dicts

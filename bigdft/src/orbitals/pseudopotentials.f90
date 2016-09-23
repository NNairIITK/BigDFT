!> @file
!!  Define handling of the psp parameters
!! @author
!!    Copyright (C) 2015-2015 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module pseudopotentials

  use module_defs, only: gp
  use m_pawrad, only: pawrad_type
  use m_pawtab, only: pawtab_type
  use dictionaries
  use public_keys
  use public_enums
  
  private 

  integer, parameter :: SKIP=0
  integer, parameter :: DIAGONAL=1
  integer, parameter :: FULL=2

  integer, parameter :: NRLOC_MAX=2 !<maximum number of rlocs for the local potential
  integer, parameter :: NCOEFF_MAX=4 !<maximum number of coefficients
  integer, parameter :: LMAX=3 !<=maximum L
  integer, parameter :: IMAX=3 !<number of proncipla quantum numbers for HGH
  integer, parameter :: NG_CORE_MAX=1 !<max no. of gaussians for core charge density case
  integer, parameter :: NG_VAL_MAX=0 !<max no. of gaussians for valence charge density, in the core region

  type, public :: atomic_proj_coeff
     real(gp) :: hij !<coefficient
     real(gp), dimension(:,:), pointer :: mat !<matrix of nonlocal projectors in the mm' space
  end type atomic_proj_coeff

  type, public :: PSP_data
     integer :: nelpsp
     integer :: npspcode          !< PSP codes (see @link psp_projectors::pspcode_hgh @endlink)
     integer :: ixcpsp            !< PSP ixc code
     integer :: nzatom            !< Atomic number
     integer :: iradii_source     !< Source of the radii coefficients (Hard-Coded, PSP File, ...)
     integer :: ngv !<number of gaussians for core charge
     integer :: ngc !<number of gaussians for valence charge
     real(gp), dimension(NRLOC_MAX) :: rloc !<radii of the local potential
     real(gp), dimension(NCOEFF_MAX) :: c_sr !<coefficients for the local short-range part of the potential
     real(gp), dimension(IMAX,IMAX,0:LMAX) :: h_ij  !< Pseudopotential matrix for HGH case, always symmetric
     real(gp), dimension(NG_CORE_MAX) :: core_c !< core charge case
     real(gp), dimension(NG_VAL_MAX) :: core_v !< valence charge case
  end type PSP_data

  public :: psp_set_from_dict,get_psp,psp_dict_fill_all,f_free_prj_ptr
  public :: apply_hij_coeff,update_psp_dict,psp_from_stream,nullify_atomic_proj_coeff

  contains

    subroutine nullify_atomic_proj_coeff(prj)
      use f_utils, only: f_zero
      implicit none
      type(atomic_proj_coeff), intent(out) :: prj
      call f_zero(prj%hij)
      nullify(prj%mat)
    end subroutine nullify_atomic_proj_coeff
    subroutine free_atomic_proj_coeff(prj)
      use dynamic_memory
      implicit none
      type(atomic_proj_coeff), intent(inout) :: prj
      call f_free_ptr(prj%mat)
      call nullify_atomic_proj_coeff(prj)
    end subroutine free_atomic_proj_coeff

    subroutine f_free_prj_ptr(prj)
      implicit none
      type(atomic_proj_coeff), dimension(:,:,:), pointer :: prj
      !local variables
      integer :: i1,i2,i3
      if (.not. associated(prj)) return
      do i3=lbound(prj,3),ubound(prj,3)
         do i2=lbound(prj,2),ubound(prj,2)
            do i1=lbound(prj,1),ubound(prj,1)
               call free_atomic_proj_coeff(prj(i1,i2,i3))
            end do
         end do
      end do
      deallocate(prj)
      nullify(prj)
    end subroutine f_free_prj_ptr
    

    !> Fill up the dict with all pseudopotential information
    subroutine psp_dict_fill_all(dict, atomname, run_ixc, projrad, crmult, frmult)
      use module_defs, only: gp, UNINITIALIZED
      use ao_inguess, only: atomic_info
      use public_enums, only : RADII_SOURCE, RADII_SOURCE_HARD_CODED, RADII_SOURCE_FILE
      use dynamic_memory
      use dictionaries
      use public_keys, ixc_fake => ixc, projrad_fake => projrad
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

      call f_routine(id='psp_dict_fill_all')


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

      call f_release_routine()

    end subroutine psp_dict_fill_all

    !> Set the value for atoms_data from the dictionary
    subroutine psp_set_from_dict(dict, valid, &
         & nzatom, nelpsp, npspcode, ixcpsp, iradii_source, radii_cf, rloc, lcoeff, psppar)
      use module_defs, only: gp, UNINITIALIZED
      use public_enums
      use dictionaries
      use public_keys, rloc_fake => rloc
      implicit none
      !Arguments
      type(dictionary), pointer :: dict
      logical, intent(out), optional :: valid !< .true. if all required info for a pseudo are present
      integer, intent(out), optional :: nzatom, nelpsp, npspcode, ixcpsp, iradii_source
      real(gp), intent(out), optional :: rloc
      real(gp), dimension(4), intent(out), optional :: lcoeff
      real(gp), dimension(4,0:6), intent(out), optional :: psppar
      real(gp), dimension(3), intent(out), optional :: radii_cf
      !Local variables
      type(dictionary), pointer :: loc
      character(len = max_field_length) :: str
      integer :: l

      ! Default values
      if (present(valid)) valid = .true.

      ! Parameters
      if (present(nzatom)) nzatom = -1
      if (present(nelpsp)) nelpsp = -1
      if (present(ixcpsp)) ixcpsp = -1
      if (has_key(dict, ATOMIC_NUMBER) .and. present(nzatom))   nzatom = dict // ATOMIC_NUMBER
      if (has_key(dict, ELECTRON_NUMBER) .and. present(nelpsp)) nelpsp = dict // ELECTRON_NUMBER
      if (has_key(dict, PSPXC_KEY) .and. present(ixcpsp))       ixcpsp = dict // PSPXC_KEY
      if (present(valid)) valid = valid .and. has_key(dict, ATOMIC_NUMBER) .and. &
           & has_key(dict, ELECTRON_NUMBER) .and. has_key(dict, PSPXC_KEY)

      ! Local terms
      if (present(rloc))   rloc      = 0._gp
      if (present(lcoeff)) lcoeff(:) = 0._gp
      if (has_key(dict, LPSP_KEY)) then
         loc => dict // LPSP_KEY
         if (has_key(loc, "Rloc") .and. present(rloc)) rloc = loc // 'Rloc'
         if (has_key(loc, COEFF_KEY) .and. present(lcoeff)) lcoeff = loc // COEFF_KEY
         ! Validate
         if (present(valid)) valid = valid .and. has_key(loc, "Rloc") .and. &
              & has_key(loc, COEFF_KEY)
      end if

      ! Nonlocal terms
      if (present(psppar))   psppar(:,:) = 0._gp
      if (has_key(dict, NLPSP_KEY) .and. present(psppar)) then
         loc => dict_iter(dict // NLPSP_KEY)
         do while (associated(loc))
            if (has_key(loc, "Channel (l)")) then
               l = loc // "Channel (l)"
               l = l + 1
               if (has_key(loc, "Rloc"))       psppar(l,0)   = loc // 'Rloc'
               if (has_key(loc, "h_ij terms")) psppar(l,1:6) = loc // 'h_ij terms'
               if (present(valid)) valid = valid .and. has_key(loc, "Rloc") .and. &
                    & has_key(loc, "h_ij terms")
            else
               if (present(valid)) valid = .false.
            end if
            loc => dict_next(loc)
         end do
      end if

      ! Type
      if (present(npspcode)) npspcode = UNINITIALIZED(npspcode)
      if (has_key(dict, PSP_TYPE) .and. present(npspcode)) then
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
            if (present(valid) .and. (NLCC_KEY .in. dict)) &
                 valid = valid .and. &
                 ("Rcore" .in. dict // NLCC_KEY) .and.  &
                 ("Core charge" .in. dict // NLCC_KEY)
         case("PAW")
            npspcode = PSPCODE_PAW
         case default
            if (present(valid)) valid = .false.
         end select
      end if

      ! Optional values.
      if (present(iradii_source)) iradii_source = RADII_SOURCE_HARD_CODED
      if (present(radii_cf))      radii_cf(:)   = UNINITIALIZED(1._gp)
      if (has_key(dict, RADII_KEY)) then
         loc => dict // RADII_KEY
         if (has_key(loc, COARSE) .and. present(radii_cf))     radii_cf(1) = loc // COARSE
         if (has_key(loc, FINE) .and. present(radii_cf))       radii_cf(2) = loc // FINE
         if (has_key(loc, COARSE_PSP) .and. present(radii_cf)) radii_cf(3) = loc // COARSE_PSP

         if (has_key(loc, SOURCE_KEY) .and. present(iradii_source)) then
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
      end if

    end subroutine psp_set_from_dict

    subroutine update_psp_dict(dict,type)
      use yaml_output, only: yaml_warning
      use yaml_strings
      implicit none
      character(len=*), intent(in) :: type
      type(dictionary), pointer :: dict
      !local variables
      logical :: exists
      character(len=27) :: key
      character(len = max_field_length) :: str
      type(dictionary), pointer :: tmp
      key = 'psppar.' // trim(type)
      exists=key .in. dict

      if (exists) then
         tmp => dict // key 
         if (SOURCE_KEY .in. tmp) then
            str = dict_value(dict // key // SOURCE_KEY)
         else
            str = dict_value(dict // key)
         end if
         if (len_trim(str) /= 0 .and. trim(str) /= TYPE_DICT) then
            !Read the PSP file and merge to dict
            if (trim(str) /= TYPE_LIST) then
               call psp_file_merge_to_dict(dict, key, filename = trim(str))
               if (PSPXC_KEY .notin. tmp) &
                    call yaml_warning("Pseudopotential file '"+str+ &
                    "' not found. Fallback to file '"+key+&
                    "' or hard-coded pseudopotential.")
            else
               call psp_file_merge_to_dict(dict, key, lstring = dict // key)
            end if
         end if
         exists = PSPXC_KEY .in. dict // key
      end if
      if (.not. exists) call psp_file_merge_to_dict(dict, key, key)
      exists = key .in. dict
      if (exists) exists = NLCC_KEY .in.  dict // key
      if (.not. exists) call nlcc_file_merge_to_dict(dict, key, 'nlcc.' // trim(type))
    end subroutine update_psp_dict


    subroutine get_psp(dict,ityp,ntypes,nzatom,nelpsp,npspcode,ixcpsp,iradii_source,psppar,radii_cf,pawpatch,&
         pawrad,pawtab,epsatm)
      use dynamic_memory
      use m_libpaw_libxc, only: libxc_functionals_init, libxc_functionals_end
      !use libxc_functionals, only: libxc_functionals_init, libxc_functionals_end
      use m_pawtab, only: pawtab_nullify
      implicit none
      type(dictionary), pointer :: dict
      integer, intent(in) :: ntypes,ityp
      logical, intent(inout) :: pawpatch
      integer, intent(out) :: nelpsp
      integer, intent(out) :: npspcode    !< PSP codes (see @link psp_projectors::pspcode_hgh @endlink)
      integer, intent(out) :: ixcpsp            !< PSP ixc code
      integer, intent(out) :: nzatom            !< Atomic number
      integer, intent(out) :: iradii_source     !< Source of the radii coefficients (Hard-Coded, PSP File, ...)
      real(gp), dimension(3) :: radii_cf !< radii of the spheres for the wavelet descriptors
      real(gp), dimension(:), pointer :: epsatm
      real(gp), dimension(0:4,0:6), intent(out) :: psppar !< Pseudopotential parameters (HGH NL section)
      type(pawrad_type), dimension(:), pointer :: pawrad  !< PAW radial objects.
      type(pawtab_type), dimension(:), pointer :: pawtab  !< PAW objects for something.
      !local variables
      logical :: l
      integer :: ityp2
      real(gp) :: rloc
      character(len = max_field_length) :: fpaw
      real(gp), dimension(4) :: lcoeff
      real(gp), dimension(4,0:6) :: psppar_tmp

      call psp_set_from_dict(dict, l, &
           nzatom, nelpsp, npspcode, &
           ixcpsp, iradii_source, radii_cf, rloc, lcoeff, psppar_tmp)
      !To eliminate the runtime warning due to the copy of the array (TD)
      !radii_cf(ityp,:)=radii_cf(:)
      psppar(0,0)=rloc
      psppar(0,1:4)=lcoeff
      psppar(1:4,0:6)=psppar_tmp

      l = .false.
      l= dict .get. "PAW patch"
      pawpatch = pawpatch .and. l

      ! PAW case.
      if (l .and. npspcode == PSPCODE_PAW) then
         ! Allocate the PAW arrays on the fly.
         !TO be removed in the new PSP structure
         if (.not. associated(pawrad)) then
            allocate(pawrad(ntypes))
            allocate(pawtab(ntypes))
            epsatm = f_malloc_ptr(ntypes, id = "epsatm")
            do ityp2 = 1, ntypes
               !call pawrad_nullify(atoms%pawrad(ityp2))
               call pawtab_nullify(pawtab(ityp2))
            end do
         end if
         ! Re-read the pseudo for PAW arrays.
         fpaw = dict // SOURCE_KEY
         call libxc_functionals_init(ixcpsp, 1)
         call paw_from_file(pawrad(ityp), pawtab(ityp), &
              epsatm(ityp), trim(fpaw), &
              & nzatom, nelpsp, ixcpsp)
         radii_cf(3)= pawtab(ityp)%rpaw !/ frmult + 0.01
         call libxc_functionals_end()
      end if

    end subroutine get_psp

    subroutine nlcc_file_merge_to_dict(dict, key, filename)
      use module_defs, only: gp, UNINITIALIZED
      use yaml_strings
      use dictionaries
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
            call add(psp // NLCC_KEY // "Valence", gauss)
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
            call add(psp // NLCC_KEY // "Core", gauss)
         end do
      end if

      close(unit=79)
    end subroutine nlcc_file_merge_to_dict

    !> Read psp file and merge to dict
    subroutine psp_file_merge_to_dict(dict, key, filename, lstring)
      use module_defs, only: gp, UNINITIALIZED
      !      use yaml_strings
      use f_utils
      use yaml_output
      use dictionaries
      use public_keys, only: SOURCE_KEY
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

      if (has_key(dict, key) .and. trim(dict_value(dict // key)) == TYPE_LIST) &
           call dict_remove(dict, key)

      call psp_data_merge_to_dict(dict // key, nzatom, nelpsp, npspcode, ixcpsp, &
           psppar, radii_cf, rcore, qcore)
      call set(dict // key // "PAW patch", pawpatch)
      if (present(filename)) then
         call set(dict // key // SOURCE_KEY, filename)
      else
         call set(dict // key // SOURCE_KEY, "In-line")
      end if
    end subroutine psp_file_merge_to_dict


    !> Merge all psp data (coming from a file) in the dictionary
    subroutine psp_data_merge_to_dict(dict, nzatom, nelpsp, npspcode, ixcpsp, &
         & psppar, radii_cf, rcore, qcore)
      use module_defs, only: gp, UNINITIALIZED
      use public_enums, only: RADII_SOURCE_FILE,PSPCODE_GTH, PSPCODE_HGH, &
           PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
      use yaml_strings
      use dictionaries
      use public_keys
      use yaml_output
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
            call add(channel // COEFF_KEY, psppar(0,i))
         end do
         call set(dict // LPSP_KEY, channel)
      end if

      ! nlcc term
      if (npspcode == PSPCODE_HGH_K_NLCC) then
         call set(dict // NLCC_KEY, &
              dict_new( 'Rcore' .is. yaml_toa(rcore), &
              'Core charge' .is. yaml_toa(qcore)))
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
            call add(dict // NLPSP_KEY, channel)
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
    subroutine psp_from_stream(ios, nzatom, nelpsp, npspcode, &
         & ixcpsp, psppar, donlcc, rcore, qcore, radii_cf, pawpatch)
      use module_defs, only: gp, dp, UNINITIALIZED, pi_param, BIGDFT_INPUT_VARIABLES_ERROR
      use module_xc, only: xc_get_id_from_name
      use yaml_strings, only: yaml_toa
      use public_enums, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_PAW, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC
      use ao_inguess
      use dictionaries
      use f_utils
      use m_pawpsp, only: pawpsp_read_header_2
      implicit none

      type(io_stream), intent(inout) :: ios
      integer, intent(out) :: nzatom, nelpsp, npspcode, ixcpsp
      real(gp), intent(out) :: psppar(0:4,0:6), radii_cf(3), rcore, qcore
      logical, intent(out) :: pawpatch
      logical, intent(inout) ::  donlcc

      !ALEX: Some local variables
      real(gp):: sqrt2pi
      character(len=2) :: symbol

      integer :: ierror, ierror1, i, j, nn, nlterms, nprl, l, nzatom_, nelpsp_, npspcode_
      integer :: lmax,lloc,mmax, ixc_
      integer:: pspversion,basis_size,lmn_size
      real(dp) :: nelpsp_dp,nzatom_dp,r2well
      character(len=max_field_length) :: line
      logical :: exists, eof

      radii_cf = UNINITIALIZED(1._gp)
      pawpatch = .false.
      !inquire(file=trim(filename),exist=exists)
      !if (.not. exists) return

      ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
      !open(unit=11,file=trim(filename),status='old',iostat=ierror)
      !Check the open statement
      !if (ierror /= 0) then
      !   write(*,*) ': Failed to open the file (it must be in ABINIT format!): "',&
      !        trim(filename),'"'
      !   stop
      !end if

      call f_iostream_get_line(ios, line)
      if (line(1:5)/='<?xml') then
         call f_iostream_get_line(ios, line)
         read(line,*) nzatom_dp, nelpsp_dp
         nzatom=int(nzatom_dp); nelpsp=int(nelpsp_dp)
         call f_iostream_get_line(ios, line)
         read(line,*) npspcode, ixcpsp, lmax, lloc, mmax, r2well
      else
         npspcode = PSPCODE_PAW ! XML pseudo, try to get nzatom, nelpsp and ixcpsp by hand
         nzatom = 0
         nelpsp = 0
         ixcpsp = 0
         do

            call f_iostream_get_line(ios, line, eof)
            if (eof .or. (nzatom > 0 .and. nelpsp > 0 .and. ixcpsp > 0)) exit

            if (line(1:6) == "<atom ") then
               i = index(line, " Z")
               i = i + index(line(i:max_field_length), '"')
               j = i + index(line(i:max_field_length), '"') - 2
               read(line(i:j), *) nzatom_dp
               i = index(line, " valence")
               i = i + index(line(i:max_field_length), '"')
               j = i + index(line(i:max_field_length), '"') - 2
               read(line(i:j), *) nelpsp_dp
               nzatom = int(nzatom_dp)
               nelpsp = int(nelpsp_dp)
            end if
            if (line(1:15) == "<xc_functional ") then
               i = index(line, "name")
               i = i + index(line(i:max_field_length), '"')
               j = i + index(line(i:max_field_length), '"') - 2
               call xc_get_id_from_name(ixcpsp, line(i:j))
            end if
         end do
      end if

      psppar(:,:)=0._gp
      if (npspcode == PSPCODE_GTH) then !GTH case
         call f_iostream_get_line(ios, line)
         read(line,*) (psppar(0,j),j=0,4)
         do i=1,2
            call f_iostream_get_line(ios, line)
            read(line,*) (psppar(i,j),j=0,3-i)
         enddo
      else if (npspcode == PSPCODE_HGH) then !HGH case
         call f_iostream_get_line(ios, line)
         read(line,*) (psppar(0,j),j=0,4)
         call f_iostream_get_line(ios, line)
         read(line,*) (psppar(1,j),j=0,3)
         do i=2,4
            call f_iostream_get_line(ios, line)
            read(line,*) (psppar(i,j),j=0,3)
            !ALEX: Maybe this can prevent reading errors on CRAY machines?
            call f_iostream_get_line(ios, line)
            !read(11,*) skip !k coefficients, not used for the moment (no spin-orbit coupling)
         enddo
      else if (npspcode == PSPCODE_PAW) then !PAW Pseudos
         ! Need NC psp for input guess.
         ixc_=ixcpsp
         call atomic_info(nzatom, nelpsp, symbol = symbol)
         call psp_from_data(symbol, nzatom_, nelpsp_, npspcode_, ixc_, &
              psppar, exists)
         if (.not.exists) &
              call f_err_throw('Implement here.',err_name='BIGDFT_RUNTIME_ERROR')
         ! PAW format using libPAW.
         if (ios%iunit /=0) call pawpsp_read_header_2(ios%iunit,pspversion,basis_size,lmn_size)
         ! PAW data will not be saved in the input dictionary,
         ! we keep their reading for later.
         pawpatch = .true.
      else if (npspcode == PSPCODE_HGH_K) then !HGH-K case
         call f_iostream_get_line(ios, line)
         read(line,*) psppar(0,0),nn,(psppar(0,j),j=1,nn) !local PSP parameters
         call f_iostream_get_line(ios, line)
         read(line,*) nlterms !number of channels of the pseudo
         prjloop: do l=1,nlterms
            call f_iostream_get_line(ios, line)
            read(line,*) psppar(l,0),nprl,psppar(l,1),&
                 (psppar(l,j+2),j=2,nprl) !h_ij terms
            do i=2,nprl
               call f_iostream_get_line(ios, line)
               read(line,*) psppar(l,i),(psppar(l,i+j+1),j=i+1,nprl) !h_ij 
            end do
            if (l==1) cycle
            do i=1,nprl
               !ALEX: Maybe this can prevent reading errors on CRAY machines?
               call f_iostream_get_line(ios, line)
               !read(11,*)skip !k coefficients, not used
            end do
         end do prjloop
         !ALEX: Add support for reading NLCC from psppar
      else if (npspcode == PSPCODE_HGH_K_NLCC) then !HGH-NLCC: Same as HGH-K + one additional line
         call f_iostream_get_line(ios, line)
         read(line,*) psppar(0,0),nn,(psppar(0,j),j=1,nn) !local PSP parameters
         call f_iostream_get_line(ios, line)
         read(line,*) nlterms !number of channels of the pseudo
         do l=1,nlterms
            call f_iostream_get_line(ios, line)
            read(line,*) psppar(l,0),nprl,psppar(l,1),&
                 (psppar(l,j+2),j=2,nprl) !h_ij terms
            do i=2,nprl
               call f_iostream_get_line(ios, line)
               read(line,*) psppar(l,i),(psppar(l,i+j+1),j=i+1,nprl) !h_ij
            end do
            if (l==1) cycle
            do i=1,nprl
               call f_iostream_get_line(ios, line)
            end do
         end do
         call f_iostream_get_line(ios, line)
         read(line,*) rcore, qcore
         !convert the core charge fraction qcore to the amplitude of the Gaussian
         !multiplied by 4pi. This is the convention used in nlccpar(1,:).
         !fourpi=4.0_gp*pi_param!8.0_gp*dacos(0.0_gp)
         sqrt2pi=sqrt(0.5_gp*4.0_gp*pi_param)
         qcore=4.0_gp*pi_param*qcore*real(nzatom-nelpsp,gp)/&
              (sqrt2pi*rcore)**3
         donlcc=.true.
      else
         call f_err_throw('PSP code not recognised (' // trim(yaml_toa(npspcode)) // ')', &
              err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      end if

      if (npspcode /= PSPCODE_PAW) then

         !old way of calculating the radii, requires modification of the PSP files
         call f_iostream_get_line(ios, line, eof)
         if (eof) then
            !if (iproc ==0) write(*,*)&
            !     ' WARNING: last line of pseudopotential missing, put an empty line'
            line=''
         end if
         read(line,*,iostat=ierror1) radii_cf(1),radii_cf(2),radii_cf(3)
         if (ierror1 /= 0 ) then
            read(line,*,iostat=ierror) radii_cf(1),radii_cf(2)
            radii_cf(3)=UNINITIALIZED(1._gp)
            ! Open64 behaviour, if line is PAWPATCH, then radii_cf(1) = 0.
            if (ierror /= 0) radii_cf = UNINITIALIZED(1._gp)
         end if
         pawpatch = (trim(line) == "PAWPATCH")
         do
            call f_iostream_get_line(ios, line, eof)
            if (eof .or. pawpatch) exit
            pawpatch = (trim(line) == "PAWPATCH")
         end do
      else
         ! Need NC psp for input guess.
         call atomic_info(nzatom, nelpsp, symbol = symbol)
         ixc_ = ixcpsp ! Because psp_from_data() will change ixc_.
         call psp_from_data(symbol, nzatom_, nelpsp_, npspcode_, ixc_, &
              & psppar, exists)
         ! Fallback to LDA case, anyway, this is only for input guess.
         if (.not.exists) then
            ixc_ = 1
            call psp_from_data(symbol, nzatom_, nelpsp_, npspcode_, ixc_, &
                 psppar, exists)
            if (.not. exists) call f_err_throw('Default PSP data for LDA not found',&
                 err_name='BIGDFT_RUNTIME_ERROR')
         end if
      end if
    END SUBROUTINE psp_from_stream

    subroutine paw_from_file(pawrad, pawtab, epsatm, filename, nzatom, nelpsp, ixc)
      use module_defs, only: dp, gp, pi_param
      use module_base, only: bigdft_mpi
      use module_xc
      use abi_defs_basis, only: tol14, fnlen
      use m_pawpsp, only: pawpsp_main
      use m_pawxmlps, only: paw_setup, rdpawpsxml, ipsp2xml, paw_setup_free
      use m_pawrad, only: pawrad_type !, pawrad_nullify
      use m_pawtab, only: pawtab_type, pawtab_nullify
      use dictionaries, only: max_field_length
      use dynamic_memory
      use f_utils

      implicit none

      type(pawrad_type), intent(out) :: pawrad
      type(pawtab_type), intent(out) :: pawtab
      real(gp), intent(out) :: epsatm
      character(len = *), intent(in) :: filename
      integer, intent(in) :: nzatom, nelpsp, ixc

      integer:: icoulomb,ipsp !, ib, i, ii
      integer:: pawxcdev,usewvl,usexcnhat,xclevel
      integer::pspso
      real(dp) :: xc_denpos
      real(dp) :: xcccrc
      character(len = fnlen) :: filpsp   ! name of the psp file
      character(len = max_field_length) :: line
      type(io_stream) :: ios
      type(xc_info) :: xc
      !!arrays
      integer:: wvl_ngauss(2)
      integer, parameter :: mqgrid_ff = 0, mqgrid_vl = 0
      real(dp):: qgrid_ff(mqgrid_ff),qgrid_vl(mqgrid_vl)
      real(dp):: ffspl(mqgrid_ff,2,1)
      real(dp):: vlspl(mqgrid_vl,2)

      !call pawrad_nullify(pawrad)
      call pawtab_nullify(pawtab)
      !These should be passed as arguments:
      !Defines the number of Gaussian functions for projectors
      !See ABINIT input files documentation
      wvl_ngauss=[10,10]
      icoulomb= 1 !Fake argument, this only indicates that we are inside bigdft..
      !do not change, even if icoulomb/=1
      ipsp=1      !This is relevant only for XML.
      ! For the moment, it will just work for LDA
      pspso=0 !No spin-orbit for the moment

      call xc_init(xc, ixc, XC_MIXED, 1)
      xclevel = 1 ! xclevel=XC functional level (1=LDA, 2=GGA)
      if (xc_isgga(xc)) xclevel = 2
      call xc_end(xc)

      ! Define parameters:
      pawxcdev=1; usewvl=1 ; usexcnhat=0 !default
      xc_denpos=tol14
      filpsp=trim(filename)

      ! Always allocate paw_setup, but fill it only for iproc == 0.
      allocate(paw_setup(1))
      if (bigdft_mpi%iproc == 0) then
         ! Parse the PAW file if in XML format
         call f_iostream_from_file(ios, trim(filename))
         call f_iostream_get_line(ios, line)
         call f_iostream_release(ios)
         if (line(1:5) == "<?xml") then
            call rdpawpsxml(filpsp, paw_setup(1), 789)
            ipsp2xml = f_malloc(1, id = "ipsp2xml")
            ipsp2xml(1) = 1
         end if
      end if

      call pawpsp_main( &
           & pawrad,pawtab,&
           & filpsp,usewvl,icoulomb,ixc,xclevel,pawxcdev,usexcnhat,&
           & qgrid_ff,qgrid_vl,ffspl,vlspl,epsatm,xcccrc,real(nelpsp, dp),real(nzatom, dp),&
           & wvl_ngauss,comm_mpi=bigdft_mpi%mpi_comm,psxml = paw_setup(1))

      call paw_setup_free(paw_setup(1))
      deallocate(paw_setup)
      if (allocated(ipsp2xml)) call f_free(ipsp2xml)

!!$  ii = 0
!!$  do ib = 1, pawtab%basis_size
!!$     do i = 1, pawtab%wvl%pngau(ib)
!!$        ii = ii + 1
!!$        write(80 + ib,*) pawtab%wvl%pfac(:, ii), pawtab%wvl%parg(:, ii)
!!$     end do
!!$     close(80 + ib)
!!$  end do
    END SUBROUTINE paw_from_file
    
    !> routine for applying the coefficients needed HGH-type PSP to the scalar product
    !! among wavefunctions and projectors. The coefficients are real therefore 
    !! there is no need to separate scpr in its real and imaginary part before
    !pure 
    subroutine apply_hij_coeff(prj,n_w,n_p,scpr,hscpr)
      use module_base, only: gp,wp
      use f_blas, only: f_eye,f_gemv
      implicit none
      integer, intent(in) :: n_p,n_w
!!$      real(gp), dimension(3,3,4), intent(in) :: hij
      type(atomic_proj_coeff), dimension(3,3,4), intent(in) :: prj
      real(gp), dimension(n_w,n_p), intent(in) :: scpr
      real(gp), dimension(n_w,n_p), intent(out) :: hscpr
      !local variables
      integer :: i,j,l,m,iproj,iw
      real(gp), dimension(7,3,4) :: cproj,dproj 
      logical, dimension(3,4) :: cont

      hscpr=0.0_gp

      !define the logical array to identify the point from which the block is finished
      do l=1,4
         do i=1,3
            cont(i,l)= prj(i,i,l)%hij /= 0.0_gp
!!$            cont(i,l)=(hij(i,i,l) /= 0.0_gp)
         end do
      end do


      reversed_loop: do iw=1,n_w
         dproj=0.0_gp

         iproj=1
         !fill the complete coefficients
         do l=1,4 !diagonal in l
            do i=1,3
               if (cont(i,l)) then !psppar(l,i) /= 0.0_gp) then
                  do m=1,2*l-1
                     cproj(m,i,l)=scpr(iw,iproj)
                     iproj=iproj+1
                  end do
               else
                  do m=1,2*l-1
                     cproj(m,i,l)=0.0_gp
                  end do
               end if
            end do
         end do

         !applies the hij matrix
         do l=1,4 !diagonal in l
            do i=1,3
               !the projector have always to act symmetrically
               !this will be done by pointing on the same matrices
               do j=1,3
                  if (prj(i,j,l)%hij == 0.0_gp) cycle
                  if (associated(prj(i,j,l)%mat)) then
                     !nondiagonal projector approach
                     call f_gemv(a=prj(i,j,l)%mat,alpha=prj(i,j,l)%hij,beta=1.0_wp,&
                          y=dproj(1,i,l),x=cproj(1,j,l))
                  else
                     !diagonal case
                     call f_gemv(a=f_eye(2*l-1),alpha=prj(i,j,l)%hij,beta=1.0_wp,&
                          y=dproj(1,i,l),x=cproj(1,j,l))
                  end if
!!$                  do m=1,2*l-1 !diagonal in m
!!$                     dproj(m,i,l)=dproj(m,i,l)+&
!!$                          hij(i,j,l)*cproj(m,j,l)
!!$                  end do
               end do
            end do
         end do

         !copy back the coefficient
         iproj=1
         !fill the complete coefficients
         do l=1,4 !diagonal in l
            do i=1,3
               if (cont(i,l)) then !psppar(l,i) /= 0.0_gp) then
                  do m=1,2*l-1
                     hscpr(iw,iproj)=dproj(m,i,l)
                     iproj=iproj+1
                  end do
               end if
            end do
         end do
      end do reversed_loop

    end subroutine apply_hij_coeff

end module pseudopotentials

!> Calculate the core charge described by a sum of spherical harmonics of s-channel with 
!! principal quantum number increased with a given exponent.
!! the principal quantum numbers admitted are from 1 to 4
function charge_from_gaussians(expo,rhoc)
  use module_base, only:gp
  implicit none
  real(gp), intent(in) :: expo
  real(gp), dimension(4), intent(in) :: rhoc
  real(gp) :: charge_from_gaussians

  charge_from_gaussians=rhoc(1)*expo**3+3.0_gp*rhoc(2)*expo**5+&
       15.0_gp*rhoc(3)*expo**7+105.0_gp*rhoc(4)*expo**9

end function charge_from_gaussians


!> Calculate the value of the gaussian described by a sum of spherical harmonics of s-channel with 
!! principal quantum number increased with a given exponent.
!! the principal quantum numbers admitted are from 1 to 4
function spherical_gaussian_value(r2,expo,rhoc,ider)
  use module_base, only: gp,safe_exp
  implicit none
  integer, intent(in) :: ider
  real(gp), intent(in) :: expo,r2
  real(gp), dimension(4), intent(in) :: rhoc
  real(gp) :: spherical_gaussian_value
  !local variables
  real(gp) :: arg
  
  arg=r2/(expo**2)
!added underflow to evaluation of the density to avoid fpe in ABINIT xc routines
  spherical_gaussian_value=&
       (rhoc(1)+r2*rhoc(2)+r2**2*rhoc(3)+r2**3*rhoc(4))*safe_exp(-0.5_gp*arg,underflow=1.d-50)
  if (ider ==1) then !first derivative with respect to r2
     spherical_gaussian_value=-0.5_gp*spherical_gaussian_value/(expo**2)+&
         (rhoc(2)+2.0_gp*r2*rhoc(3)+3.0_gp*r2**2*rhoc(4))*safe_exp(-0.5_gp*arg,underflow=1.d-50)           
     !other derivatives to be implemented
  end if

end function spherical_gaussian_value


!> External routine as the psppar parameters are often passed by address
subroutine hgh_hij_matrix(npspcode,psppar,hij)
  use module_defs, only: gp
  use public_enums, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
  implicit none
  !Arguments
  integer, intent(in) :: npspcode
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(gp), dimension(3,3,4), intent(out) :: hij
  !Local variables
  integer :: l,i,j
  real(gp), dimension(2,2,3) :: offdiagarr

  !enter the coefficients for the off-diagonal terms (HGH case, npspcode=PSPCODE_HGH)
  offdiagarr(1,1,1)=-0.5_gp*sqrt(3._gp/5._gp)
  offdiagarr(2,1,1)=-0.5_gp*sqrt(100._gp/63._gp)
  offdiagarr(1,2,1)=0.5_gp*sqrt(5._gp/21._gp)
  offdiagarr(2,2,1)=0.0_gp !never used
  offdiagarr(1,1,2)=-0.5_gp*sqrt(5._gp/7._gp)  
  offdiagarr(2,1,2)=-7._gp/3._gp*sqrt(1._gp/11._gp)
  offdiagarr(1,2,2)=1._gp/6._gp*sqrt(35._gp/11._gp)
  offdiagarr(2,2,2)=0.0_gp !never used
  offdiagarr(1,1,3)=-0.5_gp*sqrt(7._gp/9._gp)
  offdiagarr(2,1,3)=-9._gp*sqrt(1._gp/143._gp)
  offdiagarr(1,2,3)=0.5_gp*sqrt(63._gp/143._gp)
  offdiagarr(2,2,3)=0.0_gp !never used

  !  call to_zero(3*3*4,hij(1,1,1))
  hij=0.0_gp

  do l=1,4
     !term for all npspcodes
     loop_diag: do i=1,3
        hij(i,i,l)=psppar(l,i) !diagonal term
        if ((npspcode == PSPCODE_HGH .and. l/=4 .and. i/=3) .or. &
             ((npspcode == PSPCODE_HGH_K .or. npspcode == PSPCODE_HGH_K_NLCC) .and. i/=3)) then !HGH(-K) case, offdiagonal terms
           loop_offdiag: do j=i+1,3
              if (psppar(l,j) == 0.0_gp) exit loop_offdiag
              !offdiagonal HGH term
              if (npspcode == PSPCODE_HGH) then !traditional HGH convention
                 hij(i,j,l)=offdiagarr(i,j-i,l)*psppar(l,j)
              else !HGH-K convention
                 hij(i,j,l)=psppar(l,i+j+1)
              end if
              hij(j,i,l)=hij(i,j,l) !symmetrization
           end do loop_offdiag
        end if
     end do loop_diag
  end do
end subroutine hgh_hij_matrix

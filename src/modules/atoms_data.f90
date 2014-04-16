!> @file
!! Routines associated to the generation of data for atoms of the system
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Handling of input guess creation from basis of atomic orbitals
module module_atoms
  use module_defs, only: dp,gp
  use ao_inguess, only: aoig_data
  implicit none
  private

  !> Quantities used for the symmetry operators. To be used in atomic_structure derived type.
  type, public :: symmetry_data
     integer :: symObj    !< The symmetry object from ABINIT
     integer, dimension(:,:,:), pointer :: irrzon
     real(dp), dimension(:,:,:), pointer :: phnons
  end type symmetry_data
  
  !>Structure of the system. This derived type contains the information about the physical properties
  type, public :: atomic_structure
     character(len=1) :: geocode          !< @copydoc poisson_solver::doc::geocode
     character(len=5) :: inputfile_format !< Can be xyz ascii or yaml
     character(len=20) :: units           !< Can be angstroem or bohr 
     integer :: nat                       !< Number of atoms
     integer :: ntypes                    !< Number of atomic species in the structure
     real(gp), dimension(3) :: cell_dim   !< Dimensions of the simulation domain (each one periodic or free according to geocode)
     !pointers
     real(gp), dimension(:,:), pointer :: rxyz !< Atomic positions (always in AU, units variable is considered for I/O only)
     character(len=20), dimension(:), pointer :: atomnames !< Atomic species names
     integer, dimension(:), pointer :: iatype              !< Atomic species id
     integer, dimension(:), pointer :: ifrztyp             !< Freeze atoms while updating structure
     integer, dimension(:), pointer :: input_polarization  !< Used in AO generation for WFN input guess
     type(symmetry_data) :: sym                      !< The symmetry operators
  end type atomic_structure

  !> Data containing the information about the atoms in the system
  type, public :: atoms_data
     type(atomic_structure) :: astruct
     type(aoig_data), dimension(:), pointer :: aoig !< contains the information needed for generating the AO inputguess data for each atom
     integer :: natsc    !< number of atoms with semicore occupations at the input guess
!     integer, dimension(:), pointer :: iasctype
     integer, dimension(:), pointer :: nelpsp
     integer, dimension(:), pointer :: npspcode
     integer, dimension(:), pointer :: ixcpsp
     integer, dimension(:), pointer :: nzatom
     real(gp), dimension(:,:), pointer :: radii_cf  !< user defined radii_cf, overridden in sysprop.f90
     real(gp), dimension(:), pointer :: amu         !< amu(ntypes)  Atomic Mass Unit for each type of atoms
     !real(gp), dimension(:,:), pointer :: rloc !< localization regions for parameters of linear, to be moved somewhere else
     real(gp), dimension(:,:,:), pointer :: psppar  !< pseudopotential parameters (HGH SR section)
     logical :: donlcc                              !< activate non-linear core correction treatment
     integer, dimension(:), pointer :: nlcc_ngv,nlcc_ngc   !<number of valence and core gaussians describing NLCC 
     real(gp), dimension(:,:), pointer :: nlccpar    !< parameters for the non-linear core correction, if present
!     real(gp), dimension(:,:), pointer :: ig_nlccpar !< parameters for the input NLCC

     !! for abscalc with pawpatch
     integer, dimension(:), pointer ::  paw_NofL, paw_l, paw_nofchannels
     integer, dimension(:), pointer ::  paw_nofgaussians
     real(gp), dimension(:), pointer :: paw_Greal, paw_Gimag, paw_Gcoeffs
     real(gp), dimension(:), pointer :: paw_H_matrices, paw_S_matrices, paw_Sm1_matrices
     integer :: iat_absorber 
  end type atoms_data

  public :: atoms_data_null,nullify_atoms_data,deallocate_atoms_data
  public :: atomic_structure_null,nullify_atomic_structure,deallocate_atomic_structure
  public :: deallocate_symmetry_data,set_symmetry_data
  public :: set_astruct_from_file,set_astruct_from_dict
  public :: allocate_atoms_data

  contains

    !> Creators and destructors
    pure function symmetry_data_null() result(sym)
      implicit none
      type(symmetry_data) :: sym
      call nullify_symmetry_data(sym)
    end function symmetry_data_null
    pure subroutine nullify_symmetry_data(sym)
      type(symmetry_data), intent(out) :: sym
      sym%symObj=-1
      nullify(sym%irrzon)
      nullify(sym%phnons)
    end subroutine nullify_symmetry_data

    pure function atomic_structure_null() result(astruct)
      implicit none
      type(atomic_structure) :: astruct
      call nullify_atomic_structure(astruct)
    end function atomic_structure_null
    pure subroutine nullify_atomic_structure(astruct)
      implicit none
      type(atomic_structure), intent(out) :: astruct
      astruct%geocode='X'
      astruct%inputfile_format=repeat(' ',len(astruct%inputfile_format))
      astruct%units=repeat(' ',len(astruct%units))
      astruct%nat=-1
      astruct%ntypes=-1
      astruct%cell_dim=0.0_gp
      nullify(astruct%input_polarization)
      nullify(astruct%ifrztyp)
      nullify(astruct%atomnames)
      nullify(astruct%iatype)
      nullify(astruct%rxyz)
      call nullify_symmetry_data(astruct%sym)
    end subroutine nullify_atomic_structure

    pure function atoms_data_null() result(at)
      type(atoms_data) :: at
      call nullify_atoms_data(at)
    end function atoms_data_null
    pure subroutine nullify_atoms_data(at)
      implicit none
      type(atoms_data), intent(out) :: at
      call nullify_atomic_structure(at%astruct)
      nullify(at%aoig)
      at%natsc=0
      at%donlcc=.false.
      at%iat_absorber=-1
      !     nullify(at%iasctype)
      nullify(at%nelpsp)
      nullify(at%npspcode)
      nullify(at%ixcpsp)
      nullify(at%nzatom)
      nullify(at%radii_cf)
      nullify(at%amu)
      !     nullify(at%aocc)
      !nullify(at%rloc)
      nullify(at%psppar)
      nullify(at%nlcc_ngv)
      nullify(at%nlcc_ngc)
      nullify(at%nlccpar)
      nullify(at%paw_NofL)
      nullify(at%paw_l)
      nullify(at%paw_nofchannels)
      nullify(at%paw_nofgaussians)
      nullify(at%paw_Greal)
      nullify(at%paw_Gimag)
      nullify(at%paw_Gcoeffs)
      nullify(at%paw_H_matrices)
      nullify(at%paw_S_matrices)
      nullify(at%paw_Sm1_matrices)
    end subroutine nullify_atoms_data

    !> destructor of symmetry data operations
    subroutine deallocate_symmetry_data(sym)
      use dynamic_memory, only: f_free_ptr
      use m_ab6_symmetry, only: symmetry_free
      implicit none
      type(symmetry_data), intent(inout) :: sym
      !character(len = *), intent(in) :: subname
!      integer :: i_stat, i_all

      if (sym%symObj >= 0) then
         call symmetry_free(sym%symObj)
      end if
      call f_free_ptr(sym%irrzon)
      call f_free_ptr(sym%phnons)

!!$      if (associated(sym%irrzon)) then
!!$         i_all=-product(shape(sym%irrzon))*kind(sym%irrzon)
!!$         deallocate(sym%irrzon,stat=i_stat)
!!$         call memocc(i_stat,i_all,'irrzon',subname)
!!$         nullify(sym%irrzon)
!!$      end if
!!$
!!$      if (associated(sym%phnons)) then
!!$         i_all=-product(shape(sym%phnons))*kind(sym%phnons)
!!$         deallocate(sym%phnons,stat=i_stat)
!!$         call memocc(i_stat,i_all,'phnons',subname)
!!$         nullify(sym%phnons)
!!$      end if
    end subroutine deallocate_symmetry_data


    !> Deallocate the structure atoms_data.
    subroutine deallocate_atomic_structure(astruct)!,subname) 
      use dynamic_memory, only: f_free_ptr
      use module_base, only: memocc
      implicit none
      !character(len=*), intent(in) :: subname
      type(atomic_structure), intent(inout) :: astruct
      !local variables
      character(len=*), parameter :: subname='dellocate_atomic_structure' !remove
      integer :: i_stat, i_all


      ! Deallocations for the geometry part.
      if (astruct%nat > 0) then
         i_all=-product(shape(astruct%ifrztyp))*kind(astruct%ifrztyp)
         deallocate(astruct%ifrztyp,stat=i_stat)
         call memocc(i_stat,i_all,'astruct%ifrztyp',subname)
         i_all=-product(shape(astruct%iatype))*kind(astruct%iatype)
         deallocate(astruct%iatype,stat=i_stat)
         call memocc(i_stat,i_all,'astruct%iatype',subname)
         i_all=-product(shape(astruct%input_polarization))*kind(astruct%input_polarization)
         deallocate(astruct%input_polarization,stat=i_stat)
         call memocc(i_stat,i_all,'astruct%input_polarization',subname)
         i_all=-product(shape(astruct%rxyz))*kind(astruct%rxyz)
         deallocate(astruct%rxyz,stat=i_stat)
         call memocc(i_stat,i_all,'astruct%rxyz',subname)
      end if
      if (astruct%ntypes > 0) then
         i_all=-product(shape(astruct%atomnames))*kind(astruct%atomnames)
         deallocate(astruct%atomnames,stat=i_stat)
         call memocc(i_stat,i_all,'astruct%atomnames',subname)
      end if
      ! Free additional stuff.
      call deallocate_symmetry_data(astruct%sym)

    END SUBROUTINE deallocate_atomic_structure

    !> Deallocate the structure atoms_data.
    subroutine deallocate_atoms_data(atoms) 
      use module_base
      use dynamic_memory
      implicit none
      type(atoms_data), intent(inout) :: atoms
      !local variables
      character(len=*), parameter :: subname='dellocate_atoms_data' !remove
      integer :: i_stat, i_all

      ! Deallocate atomic structure
      call deallocate_atomic_structure(atoms%astruct) 

      ! Deallocations related to pseudos.
      if (associated(atoms%nzatom)) then
         i_all=-product(shape(atoms%nzatom))*kind(atoms%nzatom)
         deallocate(atoms%nzatom,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%nzatom',subname)
         i_all=-product(shape(atoms%psppar))*kind(atoms%psppar)
         deallocate(atoms%psppar,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%psppar',subname)
         i_all=-product(shape(atoms%nelpsp))*kind(atoms%nelpsp)
         deallocate(atoms%nelpsp,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%nelpsp',subname)
         i_all=-product(shape(atoms%ixcpsp))*kind(atoms%ixcpsp)
         deallocate(atoms%ixcpsp,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%ixcpsp',subname)
         i_all=-product(shape(atoms%npspcode))*kind(atoms%npspcode)
         deallocate(atoms%npspcode,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%npspcode',subname)
         i_all=-product(shape(atoms%nlcc_ngv))*kind(atoms%nlcc_ngv)
         deallocate(atoms%nlcc_ngv,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%nlcc_ngv',subname)
         i_all=-product(shape(atoms%nlcc_ngc))*kind(atoms%nlcc_ngc)
         deallocate(atoms%nlcc_ngc,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%nlcc_ngc',subname)
         i_all=-product(shape(atoms%radii_cf))*kind(atoms%radii_cf)
         deallocate(atoms%radii_cf,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%radii_cf',subname)
         ! Parameters for Linear input guess
         !i_all=-product(shape(atoms%rloc))*kind(atoms%rloc)
         !deallocate(atoms%rloc,stat=i_stat)
         !call memocc(i_stat,i_all,'atoms%rloc',subname)
         i_all=-product(shape(atoms%amu))*kind(atoms%amu)
         deallocate(atoms%amu,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%amu',subname)
      end if
      if (associated(atoms%aoig)) then
         !no flib calls for derived types for the moment
         deallocate(atoms%aoig)
         nullify(atoms%aoig)
      end if
!!$      if (associated(atoms%iasctype)) then
!!$         i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
!!$         deallocate(atoms%iasctype,stat=i_stat)
!!$         call memocc(i_stat,i_all,'atoms%iasctype',subname)
!!$         i_all=-product(shape(atoms%aocc))*kind(atoms%aocc)
!!$         deallocate(atoms%aocc,stat=i_stat)
!!$         call memocc(i_stat,i_all,'atoms%aocc',subname)
!!$      end if
!      if (associated(atoms%nlccpar)) then
         call f_free_ptr(atoms%nlccpar)
!      end if

      !  Free data for pawpatch
      if(associated(atoms%paw_l)) then
         i_all=-product(shape(atoms%paw_l ))*kind(atoms%paw_l )
         deallocate(atoms%paw_l,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_l',subname)
      end if
      if(associated(atoms%paw_NofL)) then
         i_all=-product(shape(  atoms%paw_NofL ))*kind(atoms%paw_NofL )
         deallocate(atoms%paw_NofL,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_NofL',subname)
      end if
      if(associated(atoms%paw_nofchannels)) then
         i_all=-product(shape(  atoms%paw_nofchannels ))*kind(atoms%paw_nofchannels )
         deallocate(atoms%paw_nofchannels,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_nofchannels',subname)
      end if
      if(associated(atoms%paw_nofgaussians)) then
         i_all=-product(shape(  atoms%paw_nofgaussians ))*kind(atoms%paw_nofgaussians )
         deallocate(atoms%paw_nofgaussians,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_nofgaussians',subname)
      end if
      if(associated(atoms%paw_Greal)) then
         i_all=-product(shape(  atoms%paw_Greal ))*kind(atoms%paw_Greal )
         deallocate(atoms%paw_Greal,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_Greal',subname)
      end if
      if(associated(atoms%paw_Gimag)) then
         i_all=-product(shape(  atoms%paw_Gimag ))*kind(atoms%paw_Gimag )
         deallocate(atoms%paw_Gimag,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_Gimag',subname)
      end if
      if(associated(atoms%paw_Gcoeffs)) then
         i_all=-product(shape(  atoms%paw_Gcoeffs ))*kind(atoms%paw_Gcoeffs )
         deallocate(atoms%paw_Gcoeffs,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_Gcoeffs',subname)
      end if
      if(associated(atoms%paw_H_matrices)) then
         i_all=-product(shape(  atoms%paw_H_matrices ))*kind(atoms%paw_H_matrices )
         deallocate(atoms%paw_H_matrices,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_H_matrices',subname)
      end if
      if(associated(atoms%paw_S_matrices)) then
         i_all=-product(shape(  atoms%paw_S_matrices ))*kind(atoms%paw_S_matrices )
         deallocate(atoms%paw_S_matrices,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_S_matrices',subname)
      end if
      if(associated(atoms%paw_Sm1_matrices)) then
         i_all=-product(shape(  atoms%paw_Sm1_matrices ))*kind(atoms%paw_Sm1_matrices )
         deallocate(atoms%paw_Sm1_matrices,stat=i_stat)
         call memocc(i_stat,i_all,'atoms%paw_Sm1_matrices',subname)
      end if

    END SUBROUTINE deallocate_atoms_data

    !> set irreductible Brillouin zone
    subroutine set_symmetry_data(sym, geocode, n1i, n2i, n3i, nspin)
      use module_base
      use m_ab6_kpoints
      use m_ab6_symmetry
      implicit none
      type(symmetry_data), intent(inout) :: sym
      integer, intent(in) :: n1i, n2i, n3i, nspin
      character, intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode

      character(len = *), parameter :: subname = "symmetry_set_irreductible_zone"
      integer :: i_stat, nsym, i_third
!!$      integer :: i_all
      integer, dimension(:,:,:), allocatable :: irrzon
      real(dp), dimension(:,:,:), allocatable :: phnons

      call f_free_ptr(sym%irrzon)
      call f_free_ptr(sym%phnons)

!!$      if (associated(sym%irrzon)) then
!!$         i_all=-product(shape(sym%irrzon))*kind(sym%irrzon)
!!$         deallocate(sym%irrzon,stat=i_stat)
!!$         call memocc(i_stat,i_all,'irrzon',subname)
!!$         nullify(sym%irrzon)
!!$      end if
!!$
!!$      if (associated(sym%phnons)) then
!!$         i_all=-product(shape(sym%phnons))*kind(sym%phnons)
!!$         deallocate(sym%phnons,stat=i_stat)
!!$         call memocc(i_stat,i_all,'phnons',subname)
!!$         nullify(sym%phnons)
!!$      end if

      if (sym%symObj >= 0) then
         call symmetry_get_n_sym(sym%symObj, nsym, i_stat)
         if (nsym > 1) then
            ! Current third dimension is set to 1 always
            ! since nspin == nsppol always in BigDFT
            i_third = 1
            if (geocode == "S") i_third = n2i
            sym%irrzon=f_malloc_ptr((/n1i*(n2i - i_third + 1)*n3i,2,i_third/),id='sym%irrzon')
            sym%phnons=f_malloc_ptr((/2,n1i*(n2i - i_third + 1)*n3i,i_third/),id='sym%phnons')
!!$            allocate(sym%irrzon(n1i*(n2i - i_third + 1)*n3i,2,i_third+ndebug),stat=i_stat)
!!$            call memocc(i_stat,sym%irrzon,'irrzon',subname)
!!$            allocate(sym%phnons(2,n1i*(n2i - i_third + 1)*n3i,i_third+ndebug),stat=i_stat)
!!$            call memocc(i_stat,sym%phnons,'phnons',subname)
            if (geocode /= "S") then
               call kpoints_get_irreductible_zone(sym%irrzon, sym%phnons, &
                    &   n1i, n2i, n3i, nspin, nspin, sym%symObj, i_stat)
            else
               irrzon=f_malloc((/n1i*n3i,2,1/),id='irrzon')
               phnons=f_malloc((/2,n1i*n3i,1/),id='phnons')
!!$               allocate(irrzon(n1i*n3i,2,1+ndebug),stat=i_stat)
!!$               call memocc(i_stat,irrzon,'irrzon',subname)
!!$               allocate(phnons(2,n1i*n3i,1+ndebug),stat=i_stat)
!!$               call memocc(i_stat,phnons,'phnons',subname)
               do i_third = 1, n2i, 1
                  call kpoints_get_irreductible_zone(irrzon, phnons, n1i, 1, n3i, &
                       & nspin, nspin, sym%symObj, i_stat)
                  sym%irrzon(:,:,i_third:i_third) = irrzon
                  call vcopy(2*n1i*n3i, phnons(1,1,1), 1, sym%phnons(1,1,i_third), 1)
               end do
               call f_free(irrzon)
               call f_free(phnons)
!!$               i_all=-product(shape(irrzon))*kind(irrzon)
!!$               deallocate(irrzon,stat=i_stat)
!!$               call memocc(i_stat,i_all,'irrzon',subname)
!!$               i_all=-product(shape(phnons))*kind(phnons)
!!$               deallocate(phnons,stat=i_stat)
!!$               call memocc(i_stat,i_all,'phnons',subname)
            end if
         end if
      end if

      if (.not. associated(sym%irrzon)) then
         ! Allocate anyway to small size otherwise the bounds check does not pass.
         sym%irrzon=f_malloc_ptr((/1,2,1/),id='sym%irrzon')
         sym%phnons=f_malloc_ptr((/2,1,1/),id='sym%phnons')
!!$         allocate(sym%irrzon(1,2,1+ndebug),stat=i_stat)
!!$         call memocc(i_stat,sym%irrzon,'irrzon',subname)
!!$         allocate(sym%phnons(2,1,1+ndebug),stat=i_stat)
!!$         call memocc(i_stat,sym%phnons,'phnons',subname)
      end if
    END SUBROUTINE set_symmetry_data


    ! allocations, and setters


    !> Read atomic file
    subroutine set_astruct_from_file(file,iproc,astruct,status,comment,energy,fxyz)
      use module_base
      use m_ab6_symmetry
      !use position_files
      implicit none
      character(len=*), intent(in) :: file
      integer, intent(in) :: iproc
      type(atomic_structure), intent(inout) :: astruct
      integer, intent(out), optional :: status
      real(gp), intent(out), optional :: energy
      real(gp), dimension(:,:), pointer, optional :: fxyz
      character(len = *), intent(out), optional :: comment
      !Local variables
      character(len=*), parameter :: subname='read_atomic_file'
      integer :: l, extract, i_all, i_stat
      logical :: file_exists, archive
      character(len = 128) :: filename
      character(len = 15) :: arFile
      character(len = 6) :: ext
      real(gp) :: energy_
      real(gp), dimension(:,:), pointer :: fxyz_
      character(len = 1024) :: comment_
      external :: openNextCompress, check_atoms_positions

      file_exists = .false.
      archive = .false.
      if (present(status)) status = 0
      nullify(fxyz_)

      ! Extract from archive
      if (index(file, "posout_") == 1 .or. index(file, "posmd_") == 1) then
         write(arFile, "(A)") "posout.tar.bz2"
         if (index(file, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
         inquire(FILE = trim(arFile), EXIST = file_exists)
         if (file_exists) then
!!$     call extractNextCompress(trim(arFile), len(trim(arFile)), &
!!$          & trim(file), len(trim(file)), extract, ext)
            call openNextCompress(trim(arFile), len(trim(arFile)), &
                 & trim(file), len(trim(file)), extract, ext)
            if (extract == 0) then
               write(*,*) "Can't find '", file, "' in archive."
               if (present(status)) then
                  status = 1
                  return
               else
                  stop
               end if
            end if
            archive = .true.
            write(filename, "(A)") file//'.'//trim(ext)
            write(astruct%inputfile_format, "(A)") trim(ext)
         end if
      end if

      ! Test posinp.xyz
      if (.not. file_exists) then
         inquire(FILE = file//'.xyz', EXIST = file_exists)
         if (file_exists) then
            write(filename, "(A)") file//'.xyz'!"posinp.xyz"
            write(astruct%inputfile_format, "(A)") "xyz"
            open(unit=99,file=trim(filename),status='old')
         end if
      end if
      ! Test posinp.ascii
      if (.not. file_exists) then
         inquire(FILE = file//'.ascii', EXIST = file_exists)
         if (file_exists) then
            write(filename, "(A)") file//'.ascii'!"posinp.ascii"
            write(astruct%inputfile_format, "(A)") "ascii"
            open(unit=99,file=trim(filename),status='old')
         end if
      end if
      ! Test posinp.yaml
      if (.not. file_exists) then
         inquire(FILE = file//'.yaml', EXIST = file_exists)
         if (file_exists) then
            write(filename, "(A)") file//'.yaml'!"posinp.yaml
            write(astruct%inputfile_format, "(A)") "yaml"
            ! Pb if toto.yaml because means that there is no key posinp!!
         end if
      end if
      ! Test the name directly
      if (.not. file_exists) then
         inquire(FILE = file, EXIST = file_exists)
         if (file_exists) then
            write(filename, "(A)") file
            l = len(file)
            if (file(l-3:l) == ".xyz") then
               write(astruct%inputfile_format, "(A)") "xyz"
            else if (file(l-5:l) == ".ascii") then
               write(astruct%inputfile_format, "(A)") "ascii"
            else if (file(l-4:l) == ".yaml") then
               write(astruct%inputfile_format, "(A)") "yaml"
            else
               write(*,*) "Atomic input file '" // trim(file) // "', format not recognised."
               write(*,*) " File should be *.yaml, *.ascii or *.xyz."
               if (present(status)) then
                  status = 1
                  return
               else
                  stop
               end if
            end if
            if (trim(astruct%inputfile_format) /= "yaml") then
               open(unit=99,file=trim(filename),status='old')
            end if
         end if
      end if

      if (.not. file_exists) then
         if (present(status)) then
            status = 1
            return
         else
            write(*,*) "Atomic input file not found."
            write(*,*) " Files looked for were '"//file//".yaml', '"//file//".ascii', '"//file//".xyz' and '"//file//"'."
            stop 
         end if
      end if

      if (astruct%inputfile_format == "xyz") then
         !read atomic positions
         if (.not.archive) then
            call read_xyz_positions(iproc,99,astruct,comment_,energy_,fxyz_,directGetLine)
         else
            call read_xyz_positions(iproc,99,astruct,comment_,energy_,fxyz_,archiveGetLine)
         end if
      else if (astruct%inputfile_format == "ascii") then
         i_stat = iproc
         if (present(status)) i_stat = 1
         !read atomic positions
         if (.not.archive) then
            call read_ascii_positions(i_stat,99,astruct,comment_,energy_,fxyz_,directGetLine)
         else
            call read_ascii_positions(i_stat,99,astruct,comment_,energy_,fxyz_,archiveGetLine)
         end if
      else if (astruct%inputfile_format == "yaml" .and. index(file,'posinp') /= 0) then
         ! Pb if toto.yaml because means that there is already no key posinp in the file toto.yaml!!
         write(*,*) "Atomic input file in YAML not yet supported, call 'set_astruct_from_dict()' instead."
         stop
      end if

      !Check the number of atoms
      if (astruct%nat < 0) then
         if (present(status)) then
            status = 1
            return
         else
            write(*,'(1x,3a,i0,a)') "In the file '",trim(filename),&
                 &  "', the number of atoms (",astruct%nat,") < 0 (should be >= 0)."
            stop 
         end if
      end if

      !control atom positions
      call check_atoms_positions(astruct,(iproc == 0))

      ! We delay the calculation of the symmetries.
      !this should be already in the atoms_null routine
      astruct%sym=symmetry_data_null()
      !   astruct%sym%symObj = -1
      !   nullify(astruct%sym%irrzon)
      !   nullify(astruct%sym%phnons)

      ! close open file.
      if (.not.archive .and. trim(astruct%inputfile_format) /= "yaml") then
         close(99)
!!$  else
!!$     call unlinkExtract(trim(filename), len(trim(filename)))
      end if

      ! We transfer optionals.
      if (present(energy)) then
         energy = energy_
      end if
      if (present(comment)) then
         write(comment, "(A)") comment_
      end if
      if (present(fxyz)) then
         fxyz => fxyz_
      else if (associated(fxyz_)) then
         i_all=-product(shape(fxyz_))*kind(fxyz_)
         deallocate(fxyz_,stat=i_stat)
         call memocc(i_stat,i_all,'fxyz_',subname)
      end if

    END SUBROUTINE set_astruct_from_file


    !> allocate the astruct variable from the dictionary of input data
    !retrieve also other information like the energy and the forces if requested
    !! and presend in the dictionary
    subroutine set_astruct_from_dict(dict, astruct, comment, energy, fxyz)
      use module_defs, only: gp, Bohr_Ang, UNINITIALIZED
      use dictionaries
      use dynamic_memory
      implicit none
      type(dictionary), pointer :: dict !< dictionary of the input variables
                                        !! the keys have to be declared like input_dicts module
      type(atomic_structure), intent(out) :: astruct !<structure created from the file
      !> energy, potentially written in the dictionary, atomic units
      real(gp), intent(out), optional :: energy
      !> forces: the pointer is undefined (or nullified) on entry, otherwise
      !! an allocation error might be produced
      real(gp), dimension(:,:), pointer, optional :: fxyz
      !> extra comment retrieved from the file if present
      character(len = 1024), intent(out), optional :: comment

      !local variables
      character(len=*), parameter :: subname='astruct_set_from_dict'
      type(dictionary), pointer :: pos, at
      character(len = max_field_length) :: str
      integer :: iat, ierr, ityp, units, igspin, igchrg, nsgn, ntyp
      character(len=20), dimension(100) :: atomnames

      call nullify_atomic_structure(astruct)
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
      call astruct_set_n_atoms(astruct, astruct%nat)
      ntyp = 0
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

      call astruct_set_n_types(astruct, ntyp)
      astruct%atomnames(1:ntyp) = atomnames(1:ntyp)

      if (has_key(dict, "Properties")) then
         pos => dict // "Properties"
         if (has_key(pos, "Energy (Ha)") .and. present(energy)) energy = pos // "Energy (Ha)"
         if (has_key(pos, "Info") .and. present(comment)) comment = pos // "Info"
         if (has_key(pos, "Format")) astruct%inputfile_format = pos // "Format"
      end if

    end subroutine set_astruct_from_dict

    include 'astruct-inc.f90'

    !> terminate the allocation of the memory in the pointers of atoms
    subroutine allocate_atoms_data(atoms)
      implicit none
      type(atoms_data), intent(inout) :: atoms
      external :: allocate_atoms_nat,allocate_atoms_ntypes
      
      call allocate_atoms_nat(atoms)
      call allocate_atoms_ntypes(atoms)
    end subroutine allocate_atoms_data

end module module_atoms

!> Allocation of the arrays inside the structure atoms_data, considering the part which is associated to astruct%nat
!! this routine is external to the module as it has to be called from C
subroutine astruct_set_n_atoms(astruct, nat)
  use module_base
  use module_atoms, only: atomic_structure
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: nat
  !local variables
  character(len=*), parameter :: subname='astruct_set_n_atoms' !<remove
  integer :: i_stat

  astruct%nat = nat

  ! Allocate geometry related stuff.
  allocate(astruct%iatype(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%iatype,'astruct%iatype',subname)
  allocate(astruct%ifrztyp(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%ifrztyp,'astruct%ifrztyp',subname)
  allocate(astruct%input_polarization(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%input_polarization,'astruct%input_polarization',subname)
  allocate(astruct%rxyz(3, astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%rxyz,'astruct%rxyz',subname)

  !this array is useful for frozen atoms, no atom is frozen by default
  astruct%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  astruct%input_polarization(:)=100

  if (astruct%nat > 0) call to_zero(3 * astruct%nat, astruct%rxyz(1,1))
END SUBROUTINE astruct_set_n_atoms

!> allocation of the memoey space associated to the number of types astruct%ntypes
subroutine astruct_set_n_types(astruct, ntypes)
  use module_base
  use module_atoms, only: atomic_structure
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: ntypes
  !character(len = *), intent(in) :: subname
  !local variables
  character(len=*), parameter :: subname='astruct_set_n_types' !<remove
  integer :: i, i_stat

  astruct%ntypes = ntypes

  ! Allocate geometry related stuff.
  allocate(astruct%atomnames(astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%atomnames,'astruct%atomnames',subname)

  do i = 1, astruct%ntypes, 1
     write(astruct%atomnames(i), "(A)") " "
  end do
END SUBROUTINE astruct_set_n_types

!> initialize the astruct variable from a file
subroutine astruct_set_from_file(lstat, astruct, filename)
  use module_base
  use module_atoms, only: atomic_structure,read_atomic_file=>set_astruct_from_file
  implicit none
  logical, intent(out) :: lstat
  type(atomic_structure), intent(inout) :: astruct
  character(len = *), intent(in) :: filename

  integer :: status

  call read_atomic_file(filename, 0, astruct, status)
  lstat = (status == 0)
END SUBROUTINE astruct_set_from_file

!> Calculate the symmetries and update
subroutine astruct_set_symmetries(astruct, disableSym, tol, elecfield, nspin)
  use module_base
  use module_atoms, only: atomic_structure,deallocate_symmetry_data
  use defs_basis
  use m_ab6_symmetry
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  logical, intent(in) :: disableSym
  real(gp), intent(in) :: tol
  real(gp), intent(in) :: elecfield(3)
  integer, intent(in) :: nspin
  !local variables
  character(len=*), parameter :: subname='astruct_set_symmetries'
  integer :: i_stat, ierr, i_all
  real(gp) :: rprimd(3, 3)
  real(gp), dimension(:,:), allocatable :: xRed

  ! Calculate the symmetries, if needed
  if (astruct%geocode /= 'F') then
     if (astruct%sym%symObj < 0) then
        call symmetry_new(astruct%sym%symObj)
     end if
     ! Adjust tolerance
     if (tol > 0._gp) call symmetry_set_tolerance(astruct%sym%symObj, tol, ierr)
     ! New values
     rprimd(:,:) = 0
     rprimd(1,1) = astruct%cell_dim(1)
     rprimd(2,2) = astruct%cell_dim(2)
     if (astruct%geocode == 'S') rprimd(2,2) = 1000._gp
     rprimd(3,3) = astruct%cell_dim(3)
     call symmetry_set_lattice(astruct%sym%symObj, rprimd, ierr)
     allocate(xRed(3, astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,xRed,'xRed',subname)
     xRed(1,:) = modulo(astruct%rxyz(1, :) / rprimd(1,1), 1._gp)
     xRed(2,:) = modulo(astruct%rxyz(2, :) / rprimd(2,2), 1._gp)
     xRed(3,:) = modulo(astruct%rxyz(3, :) / rprimd(3,3), 1._gp)
     call symmetry_set_structure(astruct%sym%symObj, astruct%nat, astruct%iatype, xRed, ierr)
     i_all=-product(shape(xRed))*kind(xRed)
     deallocate(xRed,stat=i_stat)
     call memocc(i_stat,i_all,'xRed',subname)
     if (astruct%geocode == 'S') then
        call symmetry_set_periodicity(astruct%sym%symObj, &
             & (/ .true., .false., .true. /), ierr)
     else if (astruct%geocode == 'F') then
        call symmetry_set_periodicity(astruct%sym%symObj, &
             & (/ .false., .false., .false. /), ierr)
     end if
     !if (all(in%elecfield(:) /= 0)) then
     !     ! I'm not sure what this subroutine does!
     !   call symmetry_set_field(astruct%sym%symObj, (/ in%elecfield(1) , in%elecfield(2),in%elecfield(3) /), ierr)
     !elseif (in%elecfield(2) /= 0) then
     !   call symmetry_set_field(astruct%sym%symObj, (/ 0._gp, in%elecfield(2), 0._gp /), ierr)
     if (elecfield(2) /= 0) then
        call symmetry_set_field(astruct%sym%symObj, (/ 0._gp, elecfield(2), 0._gp /), ierr)
     end if
     if (nspin == 2) then
        call symmetry_set_collinear_spin(astruct%sym%symObj, astruct%nat, &
             & astruct%input_polarization, ierr)
!!$     else if (in%nspin == 4) then
!!$        call symmetry_set_spin(atoms%astruct%sym%symObj, atoms%astruct%nat, &
!!$             & atoms%astruct%input_polarization, ierror)
     end if
     if (disableSym) then
        call symmetry_set_n_sym(astruct%sym%symObj, 1, &
             & reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), (/ 3 ,3, 1 /)), &
             & reshape((/ 0.d0, 0.d0, 0.d0 /), (/ 3, 1/)), (/ 1 /), ierr)
     end if
  else
     call deallocate_symmetry_data(astruct%sym)
     astruct%sym%symObj = -1
  end if
END SUBROUTINE astruct_set_symmetries

!> Allocation of the arrays inside the structure atoms_data
subroutine allocate_atoms_nat(atoms)
  use module_base
  use module_atoms, only: atoms_data
  use ao_inguess, only : aoig_data_null
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer :: iat

  allocate(atoms%aoig(atoms%astruct%nat))
  do iat=1,atoms%astruct%nat
     atoms%aoig(iat)=aoig_data_null()
  end do

!!$  ! Allocate geometry related stuff.
!!$  ! semicores useful only for the input guess
!!$  allocate(atoms%iasctype(atoms%astruct%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,atoms%iasctype,'atoms%iasctype',subname)
!!$
!!$  allocate(atoms%aocc(nelecmax,atoms%astruct%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,atoms%aocc,'atoms%aocc',subname)
END SUBROUTINE allocate_atoms_nat


subroutine allocate_atoms_ntypes(atoms)
  use module_base
  use module_atoms, only: atoms_data
  implicit none
  type(atoms_data), intent(inout) :: atoms
  !local variables
  character(len = *), parameter :: subname='allocate_atoms_ntypes'
  integer :: i_stat

  ! Allocate pseudo related stuff.
  ! store PSP parameters, modified to accept both GTH and HGHs pseudopotential types
  allocate(atoms%amu(atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%amu,'atoms%amu',subname)
  allocate(atoms%psppar(0:4,0:6,atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%psppar,'atoms%psppar',subname)
  allocate(atoms%nelpsp(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nelpsp,'atoms%nelpsp',subname)
  allocate(atoms%npspcode(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%npspcode,'atoms%npspcode',subname)
  allocate(atoms%nzatom(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nzatom,'atoms%nzatom',subname)
  allocate(atoms%ixcpsp(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%ixcpsp,'atoms%ixcpsp',subname)
  allocate(atoms%radii_cf(atoms%astruct%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%radii_cf,'atoms%radii_cf',subname)
  ! parameters for NLCC
  allocate(atoms%nlcc_ngv(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngv,'atoms%nlcc_ngv',subname)
  allocate(atoms%nlcc_ngc(atoms%astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngc,'atoms%nlcc_ngc',subname)
  ! Parameters for Linear input guess
  !allocate(atoms%rloc(atoms%astruct%ntypes,3),stat=i_stat)
  !call memocc(i_stat,atoms%rloc,'atoms%rloc',subname)
END SUBROUTINE allocate_atoms_ntypes


! Init and free routines
!> Allocate a new atoms_data type, for bindings.
subroutine atoms_new(atoms)
  use module_atoms, only: atoms_data,nullify_atoms_data
  implicit none
  type(atoms_data), pointer :: atoms

  type(atoms_data), pointer :: intern
  
  allocate(intern)
  call nullify_atoms_data(intern)
  atoms => intern
END SUBROUTINE atoms_new

!> Free an allocated atoms_data type.
subroutine atoms_free(atoms)
  use module_atoms, only: atoms_data,deallocate_atoms_data
  implicit none
  type(atoms_data), pointer :: atoms
  
  call deallocate_atoms_data(atoms)
  deallocate(atoms)
END SUBROUTINE atoms_free

!> @file
!! Routines for handling the structure atoms_data 
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Deallocate the structure atoms_data.
subroutine deallocate_atoms(atoms,subname) 
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: subname
  type(atoms_data), intent(inout) :: atoms
  !local variables
  integer :: i_stat, i_all

  ! Deallocate atomic structure
  call deallocate_atomic_structure(atoms%astruct,subname) 

  ! Deallocations for the geometry part.
  if (atoms%astruct%nat > 0) then
     i_all=-product(shape(atoms%amu))*kind(atoms%amu)
     deallocate(atoms%amu,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%amu',subname)
  end if
  if (atoms%astruct%ntypes > 0) then
     ! Parameters for Linear input guess
     i_all=-product(shape(atoms%rloc))*kind(atoms%rloc)
     deallocate(atoms%rloc,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%rloc',subname)
  end if

  ! Deallocations related to pseudos.
  if (atoms%astruct%ntypes > 0) then
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
     i_all=-product(shape(atoms%iasctype))*kind(atoms%iasctype)
     deallocate(atoms%iasctype,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%iasctype',subname)
     i_all=-product(shape(atoms%aocc))*kind(atoms%aocc)
     deallocate(atoms%aocc,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%aocc',subname)
  end if
  if (associated(atoms%nlccpar)) then
     i_all=-product(shape(atoms%nlccpar))*kind(atoms%nlccpar)
     deallocate(atoms%nlccpar,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%nlccpar',subname)
  end if

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

END SUBROUTINE deallocate_atoms

!> Deallocate the structure atoms_data.
subroutine deallocate_atomic_structure(astruct,subname) 
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: subname
  type(atomic_structure), intent(inout) :: astruct
  !local variables
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
     !i_all=-product(shape(astruct%rxyz))*kind(astruct%rxyz)
     !deallocate(astruct%rxyz,stat=i_stat)
     !call memocc(i_stat,i_all,'astruct%rxyz',subname)
  end if
  if (astruct%ntypes > 0) then
     i_all=-product(shape(astruct%atomnames))*kind(astruct%atomnames)
     deallocate(astruct%atomnames,stat=i_stat)
     call memocc(i_stat,i_all,'astruct%atomnames',subname)
  end if
  ! Free additional stuff.
  call deallocate_symmetry(astruct%sym, subname)

END SUBROUTINE deallocate_atomic_structure

!> Allocation of the arrays inside the structure atoms_data
subroutine allocate_atoms_nat(atoms, subname)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat
  integer, parameter :: nelecmax=32

  ! Allocate geometry related stuff.
  allocate(atoms%amu(atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%amu,'atoms%amu',subname)
  ! semicores useful only for the input guess
  allocate(atoms%iasctype(atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%iasctype,'atoms%iasctype',subname)

  allocate(atoms%aocc(nelecmax,atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%aocc,'atoms%aocc',subname)
END SUBROUTINE allocate_atoms_nat

!> Allocation of the arrays inside the structure atoms_data
subroutine allocate_astruct_nat(astruct, nat, subname)
  use module_base
  use module_types
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: nat
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat
  integer, parameter :: nelecmax=32

  astruct%nat = nat

  ! Allocate geometry related stuff.
  allocate(astruct%iatype(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%iatype,'astruct%iatype',subname)
  allocate(astruct%ifrztyp(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%ifrztyp,'astruct%ifrztyp',subname)
  allocate(astruct%input_polarization(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%input_polarization,'astruct%input_polarization',subname)

  !this array is useful for frozen atoms, no atom is frozen by default
  astruct%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  astruct%input_polarization(:)=100

END SUBROUTINE allocate_astruct_nat

subroutine allocate_atoms_ntypes(atoms, subname)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat

  ! Allocate pseudo related stuff.
  ! store PSP parameters, modified to accept both GTH and HGHs pseudopotential types
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
  allocate(atoms%rloc(atoms%astruct%ntypes,3),stat=i_stat)
  call memocc(i_stat,atoms%rloc,'atoms%rloc',subname)
END SUBROUTINE allocate_atoms_ntypes

subroutine allocate_astruct_ntypes(astruct, ntypes, subname)
  use module_base
  use module_types
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: ntypes
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat

  astruct%ntypes = ntypes

  ! Allocate geometry related stuff.
  allocate(astruct%atomnames(astruct%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%atomnames,'astruct%atomnames',subname)

END SUBROUTINE allocate_astruct_ntypes


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

!> Add a displacement of atomic positions and put in the box
!! @param atom    atoms_data structure
!! @param rxyz    atomic positions
!! @param randdis random displacement
subroutine atoms_set_displacement(atoms, rxyz, randdis)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  real(gp), intent(in) :: randdis

  integer :: iat
  real(gp) :: tt
  
  !Shake atoms if required.
  if (randdis > 0.d0) then
     do iat=1,atoms%astruct%nat
        if (atoms%astruct%ifrztyp(iat) == 0) then
           call random_number(tt)
           rxyz(1,iat)=rxyz(1,iat)+randdis*tt
           call random_number(tt)
           rxyz(2,iat)=rxyz(2,iat)+randdis*tt
           call random_number(tt)
           rxyz(3,iat)=rxyz(3,iat)+randdis*tt
        end if
     enddo
  end if

  !atoms inside the box.
  do iat=1,atoms%astruct%nat
     if (atoms%astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=modulo(rxyz(2,iat),atoms%astruct%cell_dim(2))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     end if
  end do
END SUBROUTINE atoms_set_displacement


!> Read atomic positions
subroutine read_xyz_positions(iproc,ifile,astruct,comment,energy,fxyz,getLine)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(out) :: energy
  real(gp), dimension(:,:), pointer :: fxyz
  character(len = 1024), intent(out) :: comment
  !Routine as argument
  interface
     subroutine getline(line,ifile,eof)
       integer, intent(in) :: ifile
       character(len=150), intent(out) :: line
       logical, intent(out) :: eof
     END SUBROUTINE getline
  end interface
  !local variables
  character(len=*), parameter :: subname='read_atomic_positions'
  character(len=20) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl, eof
  integer :: iat,ityp,ntyp,i,ierrsfx,i_stat
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  character(len=20), dimension(100) :: atomnames

  call getLine(line, ifile, eof)
  if (eof) then
     write(*,*) "Error: unexpected end of file."
     stop
  end if
  read(line,*, iostat = ierrsfx) iat,astruct%units,energy,comment
  if (ierrsfx /= 0) then
     read(line,*, iostat = ierrsfx) iat,astruct%units,energy
     write(comment, "(A)") ""
     if (ierrsfx /= 0) then
        read(line,*, iostat = ierrsfx) iat,astruct%units
        energy = UNINITIALIZED(energy)
        if (ierrsfx /= 0) then
           read(line,*, iostat = ierrsfx) iat
           write(astruct%units, "(A)") "bohr"
        end if
     end if
  else
     i = index(line, trim(comment))
     write(comment, "(A)") line(i:)
  end if

  allocate(astruct%rxyz(3,iat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%rxyz,'astruct%rxyz',subname)
  call allocate_astruct_nat(astruct, iat, subname)

  !controls if the positions are provided with machine precision
  if (astruct%units == 'angstroemd0' .or. astruct%units== 'atomicd0' .or. &
       astruct%units== 'bohrd0' .or. astruct%units=='reduced') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !read from positions of .xyz format, but accepts also the old .ascii format
  call getLine(line, ifile, eof)
  if (eof) then
     write(*,*) "Error: unexpected end of file."
     stop
  end if

!!!  !old format, still here for backward compatibility
!!!  !admits only simple precision calculation
!!!  read(line,*,iostat=ierror) rx,ry,rz,tatonam

!!!  !in case of old format, put geocode to F and alat to 0.
!!!  if (ierror == 0) then
!!!     astruct%geocode='F'
!!!     alat1d0=0.0_gp
!!!     alat2d0=0.0_gp
!!!     alat3d0=0.0_gp
!!!  else
  if (lpsdbl) then
     read(line,*,iostat=ierrsfx) tatonam,alat1d0,alat2d0,alat3d0
  else
     read(line,*,iostat=ierrsfx) tatonam,alat1,alat2,alat3
  end if
  if (ierrsfx == 0) then
     if (trim(tatonam)=='periodic') then
        astruct%geocode='P'
     else if (trim(tatonam)=='surface') then 
        astruct%geocode='S'
        astruct%cell_dim(2)=0.0_gp
     else !otherwise free bc
        astruct%geocode='F'
        astruct%cell_dim(1)=0.0_gp
        astruct%cell_dim(2)=0.0_gp
        astruct%cell_dim(3)=0.0_gp
     end if
     if (.not. lpsdbl) then
        alat1d0=real(alat1,gp)
        alat2d0=real(alat2,gp)
        alat3d0=real(alat3,gp)
     end if
  else
     astruct%geocode='F'
     alat1d0=0.0_gp
     alat2d0=0.0_gp
     alat3d0=0.0_gp
  end if
!!!  end if

  !reduced coordinates are possible only with periodic units
  if (astruct%units == 'reduced' .and. astruct%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

  !convert the values of the cell sizes in bohr
  if (astruct%units=='angstroem' .or. astruct%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     astruct%cell_dim(1)=alat1d0/Bohr_Ang
     astruct%cell_dim(2)=alat2d0/Bohr_Ang
     astruct%cell_dim(3)=alat3d0/Bohr_Ang
  else if  (astruct%units=='atomic' .or. astruct%units=='bohr'  .or.&
       astruct%units== 'atomicd0' .or. astruct%units== 'bohrd0') then
     astruct%cell_dim(1)=alat1d0
     astruct%cell_dim(2)=alat2d0
     astruct%cell_dim(3)=alat3d0
  else if (astruct%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     astruct%cell_dim(1)=alat1d0
     astruct%cell_dim(2)=alat2d0
     astruct%cell_dim(3)=alat3d0
  else
     write(*,*) 'length units in input file unrecognized'
     write(*,*) 'recognized units are angstroem or atomic = bohr'
     stop 
  endif

  ntyp=0
  do iat=1,astruct%nat
     !xyz input file, allow extra information
     call getLine(line, ifile, eof)
     if (eof) then
        write(*,*) "Error: unexpected end of file."
        stop
     end if
     if (lpsdbl) then
        read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
     else
        read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
     end if
     !print *,'extra',iat,extra
     call find_extra_info(line,extra)
     !print *,'then',iat,extra
     call parse_extra_info(iat,extra,astruct)

     tatonam=trim(symbol)
!!!     end if
     if (lpsdbl) then
        astruct%rxyz(1,iat)=rxd0
        astruct%rxyz(2,iat)=ryd0
        astruct%rxyz(3,iat)=rzd0
     else
        astruct%rxyz(1,iat)=real(rx,gp)
        astruct%rxyz(2,iat)=real(ry,gp)
        astruct%rxyz(3,iat)=real(rz,gp)
     end if

     

     if (astruct%units == 'reduced') then !add treatment for reduced coordinates
        astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),1.0_gp)
        if (astruct%geocode == 'P') astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),1.0_gp)
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),1.0_gp)
     else if (astruct%geocode == 'P') then
        astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),alat1d0)
        astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),alat2d0)
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),alat3d0)
     else if (astruct%geocode == 'S') then
        astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),alat1d0)
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),alat3d0)
     end if
 
     do ityp=1,ntyp
        if (tatonam == atomnames(ityp)) then
           astruct%iatype(iat)=ityp
           goto 200
        endif
     enddo
     ntyp=ntyp+1
     if (ntyp > 100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     astruct%iatype(iat)=ntyp
200  continue

     if (astruct%units=='angstroem' .or. astruct%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           astruct%rxyz(i,iat)=astruct%rxyz(i,iat)/Bohr_Ang
        enddo
     else if (astruct%units == 'reduced') then 
        astruct%rxyz(1,iat)=astruct%rxyz(1,iat)*astruct%cell_dim(1)
        if (astruct%geocode == 'P') astruct%rxyz(2,iat)=astruct%rxyz(2,iat)*astruct%cell_dim(2)
        astruct%rxyz(3,iat)=astruct%rxyz(3,iat)*astruct%cell_dim(3)
     endif
  enddo
  ! Try forces
  call getLine(line, ifile, eof)
  if ((.not. eof) .and. (adjustl(trim(line)) == "forces")) then
     allocate(fxyz(3,iat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)
     do iat=1,astruct%nat
        !xyz input file, allow extra information
        call getLine(line, ifile, eof)
        if (eof) then
           write(*,*) "Error: unexpected end of file."
           stop
        end if
        read(line,*,iostat=ierrsfx) symbol,fxyz(:,iat)
     end do
  end if
  !now that ntypes is determined allocate atoms%astruct%atomnames and copy the values
  call allocate_astruct_ntypes(astruct, ntyp, subname)

  astruct%atomnames(1:astruct%ntypes)=atomnames(1:astruct%ntypes)
END SUBROUTINE read_xyz_positions

!> Read atomic positions of ascii files.
subroutine read_ascii_positions(iproc,ifile,astruct,comment,energy,fxyz,getline)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(out) :: energy
  real(gp), dimension(:,:), pointer :: fxyz
  character(len = 1024), intent(out) :: comment
  interface
     subroutine getline(line,ifile,eof)
       integer, intent(in) :: ifile
       character(len=150), intent(out) :: line
       logical, intent(out) :: eof
     END SUBROUTINE getline
  end interface
  !local variables
  character(len=*), parameter :: subname='read_ascii_positions'
  character(len=20) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl, reduced, eof, forces
  integer :: iat,ntyp,ityp,i,i_stat,nlines,istart,istop,count
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3,alat4,alat5,alat6
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0,alat4d0,alat5d0,alat6d0
  character(len=20), dimension(100) :: atomnames
  ! Store the file.
  character(len = 150), dimension(5000) :: lines

  ! First pass to store the file in a string buffer.
  nlines = 1
  do
     call getline(lines(nlines), ifile, eof)
     if (eof) then
        exit
     end if
     nlines = nlines + 1
     if (nlines > 5000) then
        if (iproc==0) write(*,*) 'Atomic input file too long (> 5000 lines).'
        astruct%nat = -1
        return
     end if
  end do
  nlines = nlines - 1

  if (nlines < 4) then
     if (iproc==0) write(*,*) 'Error in ASCII file format, file has less than 4 lines.'
     astruct%nat = -1
     return
  end if

  ! Try to determine the number atoms and the keywords.
  write(astruct%units, "(A)") "bohr"
  if (lines(1)(1:1) == "#" .or. lines(1)(1:1) == "!") then
     write(comment, "(A)") adjustl(lines(1)(1:))
  else
     write(comment, "(A)") lines(1)
  end if
  reduced = .false.
  forces = .false.
  astruct%geocode = 'P'
  iat     = 0
  do i = 4, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (line(1:1) /= '#' .and. line(1:1) /= '!' .and. len(trim(line)) /= 0) then
        iat = iat + 1
     else if (line(1:8) == "#keyword" .or. line(1:8) == "!keyword") then
        if (index(line, 'bohr') > 0)        write(astruct%units, "(A)") "bohr"
        if (index(line, 'bohrd0') > 0)      write(astruct%units, "(A)") "bohrd0"
        if (index(line, 'atomic') > 0)      write(astruct%units, "(A)") "atomicd0"
        if (index(line, 'angstroem') > 0)   write(astruct%units, "(A)") "angstroem"
        if (index(line, 'angstroemd0') > 0) write(astruct%units, "(A)") "angstroemd0"
        if (index(line, 'reduced') > 0)     reduced = .true.
        if (index(line, 'periodic') > 0) astruct%geocode = 'P'
        if (index(line, 'surface') > 0)  astruct%geocode = 'S'
        if (index(line, 'freeBC') > 0)   astruct%geocode = 'F'
     else if (line(1:9) == "#metaData" .or. line(1:9) == "!metaData") then
        if (index(line, 'totalEnergy') > 0) then
           read(line(index(line, 'totalEnergy') + 12:), *, iostat = i_stat) energy
           if (i_stat /= 0) then
              energy = UNINITIALIZED(energy)
           end if
        end if
        if (index(line, 'forces') > 0) forces = .true.
     end if
  end do

  allocate(astruct%rxyz(3,iat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%rxyz,'astruct%rxyz',subname)
  call allocate_astruct_nat(astruct, iat, subname)

  !controls if the positions are provided within machine precision
  if (index(astruct%units, 'd0') > 0 .or. reduced) then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  ! Read the box definition
  astruct%cell_dim(1) = 0.0_gp
  astruct%cell_dim(2) = 0.0_gp
  astruct%cell_dim(3) = 0.0_gp
  if (lpsdbl) then
     read(lines(2),*) alat1d0,alat2d0,alat3d0
     read(lines(3),*) alat4d0,alat5d0,alat6d0
     if (alat2d0 /= 0.d0 .or. alat4d0 /= 0.d0 .or. alat5d0 /= 0.d0) then
        !if (iproc==0) 
        write(*,*) 'Only orthorombic boxes are possible.'
        astruct%nat = -1
        return
     end if
     astruct%cell_dim(1) = real(alat1d0,gp)
     astruct%cell_dim(2) = real(alat3d0,gp)
     astruct%cell_dim(3) = real(alat6d0,gp)
  else
     read(lines(2),*) alat1,alat2,alat3
     read(lines(3),*) alat4,alat5,alat6
     if (alat2 /= 0. .or. alat4 /= 0. .or. alat5 /= 0.) then
        !if (iproc==0) 
           write(*,*) 'Only orthorombic boxes are possible.'
        !if (iproc==0) 
           write(*,*) ' but alat2, alat4 and alat5 = ', alat2, alat4, alat5
        astruct%nat = -1
        return
     end if
     astruct%cell_dim(1) = real(alat1,gp)
     astruct%cell_dim(2) = real(alat3,gp)
     astruct%cell_dim(3) = real(alat6,gp)
  end if
  
  !Convert the values of the cell sizes in bohr
  if (astruct%units=='angstroem' .or. astruct%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     astruct%cell_dim(1) = astruct%cell_dim(1) / Bohr_Ang
     astruct%cell_dim(2) = astruct%cell_dim(2) / Bohr_Ang
     astruct%cell_dim(3) = astruct%cell_dim(3) / Bohr_Ang
  endif

  ntyp=0
  iat = 1
  do i = 4, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (line(1:1) /= '#' .and. line(1:1) /= '!' .and. len(trim(line)) /= 0) then
        write(extra, "(A)") "nothing"
        if (lpsdbl) then
           read(line,*, iostat = i_stat) rxd0,ryd0,rzd0,symbol,extra
           if (i_stat /= 0) read(line,*) rxd0,ryd0,rzd0,symbol
        else
           read(line,*, iostat = i_stat) rx,ry,rz,symbol,extra
           if (i_stat /= 0) read(line,*) rx,ry,rz,symbol
        end if
        call find_extra_info(line,extra)
        call parse_extra_info(iat,extra,astruct)

        tatonam=trim(symbol)

        if (lpsdbl) then
           astruct%rxyz(1,iat)=rxd0
           astruct%rxyz(2,iat)=ryd0
           astruct%rxyz(3,iat)=rzd0
        else
           astruct%rxyz(1,iat)=real(rx,gp)
           astruct%rxyz(2,iat)=real(ry,gp)
           astruct%rxyz(3,iat)=real(rz,gp)
        end if
        if (astruct%units=='angstroem' .or. astruct%units=='angstroemd0') then
           ! if Angstroem convert to Bohr
           astruct%rxyz(1,iat)=astruct%rxyz(1,iat) / Bohr_Ang
           astruct%rxyz(2,iat)=astruct%rxyz(2,iat) / Bohr_Ang
           astruct%rxyz(3,iat)=astruct%rxyz(3,iat) / Bohr_Ang
        end if

        if (reduced) then !add treatment for reduced coordinates
           astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),1.0_gp)*astruct%cell_dim(1)
           astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),1.0_gp)*astruct%cell_dim(2)
           astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),1.0_gp)*astruct%cell_dim(3)
        else if (astruct%geocode == 'P') then
           astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
           astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),astruct%cell_dim(2))
           astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        else if (astruct%geocode == 'S') then
           astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
           astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        end if

        do ityp=1,ntyp
           if (tatonam == atomnames(ityp)) then
              astruct%iatype(iat)=ityp
              goto 200
           endif
        enddo
        ntyp=ntyp+1
        if (ntyp > 100) then
           write(*,*) 'more than 100 atomnames not permitted'
           astruct%nat = -1
           return
        end if
        atomnames(ityp)=tatonam
        astruct%iatype(iat)=ntyp
200     continue

        iat = iat + 1
     end if
  enddo

  if (astruct%geocode == 'S') then
     astruct%cell_dim(2) = 0.0_gp
  else if (astruct%geocode == 'F') then
     astruct%cell_dim(1) = 0.0_gp
     astruct%cell_dim(2) = 0.0_gp
     astruct%cell_dim(3) = 0.0_gp
  end if

  if (reduced) then
     write(astruct%units, "(A)") "reduced"
  end if

  if (forces) then
     allocate(fxyz(3,astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)

     count = 0
     forces = .false.
     do i = 4, nlines, 1
        write(line, "(a150)") adjustl(lines(i))
        if ((line(1:9) == "#metaData" .or. line(1:9) == "!metaData") .and. index(line, 'forces') > 0) then
           forces = .true.
        end if
        if (forces) then
           istart = index(line, "[") + 1
           if (istart == 1) istart = index(line, "#") + 1
           do
              istop = index(line(istart:), ";") + istart - 2
              if (istop == istart - 2) exit
              read(line(istart:istop), *) fxyz(modulo(count, 3) + 1, count / 3 + 1)
              count = count + 1
              istart = istop + 2
           end do
           if (count > astruct%nat * 3) exit
        end if
     end do
  end if

  !now that ntypes is determined copy the values
  call allocate_astruct_ntypes(astruct, ntyp, subname)
  astruct%atomnames(1:astruct%ntypes)=atomnames(1:astruct%ntypes)
END SUBROUTINE read_ascii_positions

subroutine read_yaml_positions(filename, astruct, comment, energy, fxyz)
  use module_base
  use module_types
  implicit none
  character(len = *), intent(in) :: filename
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(out) :: energy
  real(gp), dimension(:,:), pointer :: fxyz
  character(len = 1024), intent(out) :: comment

  !local variables
  character(len=*), parameter :: subname='read_yaml_positions'
  integer(kind = 8) :: lst
  integer :: bc, units, i_stat, iat, i, i_all, nsgn, eunits, conv
  double precision :: acell(3), angdeg(3), gnrm, fnrm, maxval
  integer, allocatable :: igspin(:), igchrg(:)

  call f90_posinp_yaml_parse(lst, filename, len(filename))

  call f90_posinp_yaml_get_cell(lst, 0, bc, units, acell, angdeg)
  if (bc == 3) then
     astruct%geocode = 'P'
  else if (bc == 0) then
     astruct%geocode = 'F'
  else if (bc == 2) then
     astruct%geocode = 'S'
  else if (bc == 1) then
     astruct%geocode = 'W'
  end if
  if (units == 1) then
     write(astruct%units, "(A)") "angstroem"
  else if (units == 0) then
     write(astruct%units, "(A)") "bohr"
  end if
  astruct%cell_dim(1) = acell(1)
  astruct%cell_dim(2) = acell(2)
  astruct%cell_dim(3) = acell(3)
  !Convert the values of the cell sizes in bohr
  if (astruct%units=='angstroem') then
     ! if Angstroem convert to Bohr
     astruct%cell_dim(1) = astruct%cell_dim(1) / Bohr_Ang
     astruct%cell_dim(2) = astruct%cell_dim(2) / Bohr_Ang
     astruct%cell_dim(3) = astruct%cell_dim(3) / Bohr_Ang
  endif
  if (angdeg(1) /= 90. .or. angdeg(2) /= 90. .or. angdeg(3) /= 90.) then
     write(*,*) 'Only orthorombic boxes are possible.'
     write(*,*) ' but angdeg(1), angdeg(2) and angdeg(3) = ', angdeg
     astruct%nat = -1
     return
  end if
  if (astruct%geocode == 'S') then
     astruct%cell_dim(2) = 0.0_gp
  else if (astruct%geocode == 'W') then
     astruct%cell_dim(1) = 0.0_gp
     astruct%cell_dim(2) = 0.0_gp
  else if (astruct%geocode == 'F') then
     astruct%cell_dim(1) = 0.0_gp
     astruct%cell_dim(2) = 0.0_gp
     astruct%cell_dim(3) = 0.0_gp
  end if

  call f90_posinp_yaml_get_dims(lst, 0, astruct%nat, astruct%ntypes)
  if (astruct%nat == 0) then
     astruct%nat = -1
     return
  end if

  allocate(astruct%rxyz(3,astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,astruct%rxyz,'astruct%rxyz',subname)
  call allocate_astruct_nat(astruct, astruct%nat, subname)
  allocate(igspin(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,igspin,'igspin',subname)
  allocate(igchrg(astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,igchrg,'igchrg',subname)

  call f90_posinp_yaml_get_atoms(lst, 0, units, astruct%rxyz, astruct%iatype, astruct%ifrztyp, igspin, igchrg)
  if (units == 2) then
     write(astruct%units, "(A)") "reduced"
  end if
  do iat = 1, astruct%nat, 1
     if (units == 1) then
        astruct%rxyz(1,iat)=astruct%rxyz(1,iat) / Bohr_Ang
        astruct%rxyz(2,iat)=astruct%rxyz(2,iat) / Bohr_Ang
        astruct%rxyz(3,iat)=astruct%rxyz(3,iat) / Bohr_Ang
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
     if (igchrg(iat) >= 0) then
        nsgn = 1
     else
        nsgn = -1
     end if
     astruct%input_polarization(iat) = 1000 * igchrg(iat) + nsgn * 100 + igspin(iat)
  end do

  call allocate_astruct_ntypes(astruct, astruct%ntypes, subname)
  do i = 1, astruct%ntypes, 1
     call f90_posinp_yaml_get_atomname(lst, 0, i - 1, astruct%atomnames(i))
  end do

  call f90_posinp_yaml_get_comment(lst, 0, comment, 1024)
  call f90_posinp_yaml_get_properties(lst, 0, eunits, energy, gnrm, conv)
  call f90_posinp_yaml_has_forces(lst, 0, conv)
  if (conv /= 0) then
     allocate(fxyz(3,astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fxyz,'fxyz',subname)
     call f90_posinp_yaml_get_forces(lst, 0, eunits, fnrm, maxval, fxyz)
  end if

  call f90_posinp_yaml_free_list(lst)

  i_all=-product(shape(igspin))*kind(igspin)
  deallocate(igspin,stat=i_stat)
  call memocc(i_stat,i_all,'igspin',subname)

  i_all=-product(shape(igchrg))*kind(igchrg)
  deallocate(igchrg,stat=i_stat)
  call memocc(i_stat,i_all,'igchrg',subname)
END SUBROUTINE read_yaml_positions

!> Find extra information
subroutine find_extra_info(line,extra)
  implicit none
  character(len=150), intent(in) :: line
  character(len=50), intent(out) :: extra
  !local variables
  logical :: space
  integer :: i,nspace
  i=1
  space=.true.
  nspace=-1
  !print *,'line',line
  find_space : do
     !toggle the space value for each time
     if ((line(i:i) == ' ' .or. line(i:i) == char(9)) .neqv. space) then
        nspace=nspace+1
        space=.not. space
     end if
     !print *,line(i:i),nspace
     if (nspace==8) then
        extra=line(i:min(150,i+49))
        exit find_space
     end if
     if (i==150) then
        !print *,'AAA',extra
        extra='nothing'
        exit find_space
     end if
     i=i+1
  end do find_space
END SUBROUTINE find_extra_info


!> Parse extra information
subroutine parse_extra_info(iat,extra,astruct)
  use module_types
  implicit none
  !Arguments
  integer, intent(in) :: iat
  character(len=50), intent(in) :: extra
  type(atomic_structure), intent(inout) :: astruct
  !Local variables
  character(len=4) :: suffix
  logical :: go
  integer :: ierr,ierr1,ierr2,nspol,nchrg,nsgn
  !case with all the information
  !print *,iat,'ex'//trim(extra)//'ex'
  read(extra,*,iostat=ierr)nspol,nchrg,suffix
  if (extra == 'nothing') then !case with empty information
     nspol=0
     nchrg=0
     suffix='    '
  else if (ierr /= 0) then !case with partial information
     read(extra,*,iostat=ierr1)nspol,suffix
     if (ierr1 /=0) then
        call valid_frzchain(trim(extra),go)
        if (go) then
           suffix=trim(extra)
           nspol=0
           nchrg=0
        else
           read(extra,*,iostat=ierr2)nspol
           if (ierr2 /=0) then
              call error
           end if
           suffix='    '
           nchrg=0
        end if
     else
        nchrg=0
        call valid_frzchain(trim(suffix),go)
        if (.not. go) then
           read(suffix,*,iostat=ierr2)nchrg
           if (ierr2 /= 0) then
              call error
           else
              suffix='    '
           end if
        else

        end if
     end if
  end if

  !now assign the array, following the rule
  if(nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if
  astruct%input_polarization(iat)=1000*nchrg+nsgn*100+nspol

  !print *,'natpol atomic',iat,astruct%input_polarization(iat),suffix

  !convert the suffix into ifrztyp
  call frozen_ftoi(suffix,astruct%ifrztyp(iat))

!!!  if (trim(suffix) == 'f') then
!!!     !the atom is considered as blocked
!!!     astruct%ifrztyp(iat)=1
!!!  end if

contains

 subroutine error
   !if (iproc == 0) then
      print *,extra
      write(*,'(1x,a,i0,a)')&
           'ERROR in input file for atom number ',iat,&
           ': after 4th column you can put the input polarisation(s) or the frzchain: f,fxz,fy'
   !end if
   stop
 END SUBROUTINE error
  
END SUBROUTINE parse_extra_info

subroutine valid_frzchain(frzchain,go)
  implicit none
  character(len=*), intent(in) :: frzchain
  logical, intent(out) :: go

  go= trim(frzchain) == 'f' .or. &
       trim(frzchain) == 'fy' .or. &
       trim(frzchain) == 'fxz'
  
END SUBROUTINE valid_frzchain


subroutine frozen_ftoi(frzchain,ifrztyp)
  implicit none
  character(len=4), intent(in) :: frzchain
  integer, intent(out) :: ifrztyp

  if (trim(frzchain)=='') then
     ifrztyp = 0
  else if (trim(frzchain)=='f') then
     ifrztyp = 1
  else if (trim(frzchain)=='fy') then
     ifrztyp = 2
  else if (trim(frzchain)=='fxz') then
     ifrztyp = 3
  end if
        
END SUBROUTINE frozen_ftoi


!> Check the position of atoms
subroutine check_atoms_positions(iproc,astruct)
  use module_base
  use module_types
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(atomic_structure), intent(in) :: astruct
  !local variables
  integer, parameter :: iunit=9
  logical :: dowrite
  integer :: iat,nateq,jat,j

  nateq=0
  do iat=1,astruct%nat
     do jat=iat+1,astruct%nat
        if ((astruct%rxyz(1,iat)-astruct%rxyz(1,jat))**2+(astruct%rxyz(2,iat)-astruct%rxyz(2,jat))**2+&
             (astruct%rxyz(3,iat)-astruct%rxyz(3,jat))**2 ==0.0_gp) then
           nateq=nateq+1
           call yaml_warning('ERROR: atoms' // trim(yaml_toa(iat)) // ' (' // &
                & trim(astruct%atomnames(astruct%iatype(iat))) // ') and ' // &
                & trim(yaml_toa(jat)) // ' (' // trim(astruct%atomnames(astruct%iatype(jat))) // &
                & ') have the same positions')
           !write(*,'(1x,a,2(i0,a,a6,a))')'ERROR: astruct ',iat,&
           !     ' (',trim(astruct%atomnames(astruct%iatype(iat))),') and ',&
           !     jat,' (',trim(astruct%atomnames(astruct%iatype(jat))),&
           !     ') have the same positions'
        end if
     end do
  end do
  if (nateq /= 0) then
     if (iproc == 0) then
        call yaml_warning('Control your posinp file, cannot proceed')
        write(*,'(1x,a)',advance='no') 'Writing tentative alternative positions in the file posinp_alt...'
        !write(*,'(1x,a)')'Control your posinp file, cannot proceed'
        !write(*,'(1x,a)',advance='no') 'Writing tentative alternative positions in the file posinp_alt...'
        open(unit=iunit,file='posinp_alt')
        write(iunit,'(1x,a)')' ??? atomicd0'
        write(iunit,*)
        do iat=1,astruct%nat
           dowrite=.true.
           do jat=iat+1,astruct%nat
              if ((astruct%rxyz(1,iat)-astruct%rxyz(1,jat))**2+(astruct%rxyz(2,iat)-astruct%rxyz(2,jat))**2+&
                   (astruct%rxyz(3,iat)-astruct%rxyz(3,jat))**2 ==0.0_gp) then
                 dowrite=.false.
              end if
           end do
           if (dowrite) & 
                write(iunit,'(a2,4x,3(1x,1pe21.14))')trim(astruct%atomnames(astruct%iatype(iat))),&
                (astruct%rxyz(j,iat),j=1,3)
        end do
        close(unit=iunit)
        call yaml_map('Writing tentative alternative positions in the file posinp_alt',.true.)
        call yaml_warning('Replace ??? in the file heading with the actual atoms number')               
        !write(*,'(1x,a)')' done.'
        !write(*,'(1x,a)')' Replace ??? in the file heading with the actual atoms number'               
     end if
     stop 'check_atoms_positions'
  end if
END SUBROUTINE check_atoms_positions


!>Write xyz atomic file.
subroutine wtxyz(iunit,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  !local variables
  character(len=20) :: symbol
  character(len=10) :: name
  character(len=11) :: units
  character(len=50) :: extra
  integer :: iat,j
  real(gp) :: xmax,ymax,zmax,factor

  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%astruct%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%astruct%units) == 'angstroem' .or. trim(atoms%astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if

  if (energy /= 0. .and. energy /= UNINITIALIZED(energy)) then
     write(iunit,'(i6,2x,a,2x,1pe24.17,2x,a)') atoms%astruct%nat,trim(units),energy,trim(comment)
  else
     write(iunit,'(i6,2x,a,2x,a)') atoms%astruct%nat,trim(units),trim(comment)
  end if

  if (atoms%astruct%geocode == 'P') then
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  else if (atoms%astruct%geocode == 'S') then
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  else
     write(iunit,*)'free'
  end if
  do iat=1,atoms%astruct%nat
     name=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:min(len(name),5))
     end if

     call write_extra_info(extra,atoms%astruct%input_polarization(iat),atoms%astruct%ifrztyp(iat))

     write(iunit,'(a5,1x,3(1x,1pe24.17),2x,a50)')symbol,(rxyz(j,iat)*factor,j=1,3),extra
  enddo

END SUBROUTINE wtxyz

!> Add the forces in the position file for the xyz system
subroutine wtxyz_forces(iunit,fxyz,at)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=20) :: symbol
  character(len=10) :: name

  write(iunit,*)'forces'
  
  do iat=1,at%astruct%nat
     name=trim(at%astruct%atomnames(at%astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:min(len(name),5))
     end if

     write(iunit,'(a5,1x,3(1x,1pe24.17))')symbol,(fxyz(j,iat),j=1,3)
  end do
  
end subroutine wtxyz_forces

!>Write ascii file (atomic position). 
subroutine wtascii(iunit,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  !local variables
  character(len=2) :: symbol
  character(len=50) :: extra
  character(len=10) :: name
  integer :: iat,j
  real(gp) :: xmax,ymax,zmax,factor(3)

  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%astruct%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%astruct%units) == 'angstroem' .or. trim(atoms%astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
  else
     factor=1.0_gp
  end if

  write(iunit, "(A,A)") "# BigDFT file - ", trim(comment)
  write(iunit, "(3e24.17)") atoms%astruct%cell_dim(1)*factor(1), 0.d0, atoms%astruct%cell_dim(2)*factor(2)
  write(iunit, "(3e24.17)") 0.d0,                  0.d0, atoms%astruct%cell_dim(3)*factor(3)

  write(iunit, "(A,A)") "#keyword: ", trim(atoms%astruct%units)
  if (trim(atoms%astruct%units) == "reduced") write(iunit, "(A,A)") "#keyword: bohr"
  if (atoms%astruct%geocode == 'P') write(iunit, "(A)") "#keyword: periodic"
  if (atoms%astruct%geocode == 'S') write(iunit, "(A)") "#keyword: surface"
  if (atoms%astruct%geocode == 'F') write(iunit, "(A)") "#keyword: freeBC"
  if (energy /= 0.d0 .and. energy /= UNINITIALIZED(energy)) then
     write(iunit, "(A,e24.17,A)") "#metaData: totalEnergy= ", energy, " Ht"
  end if

  if (trim(atoms%astruct%units) == "reduced") then
     if (atoms%astruct%geocode == 'P' .or. atoms%astruct%geocode == 'S') factor(1) = 1._gp / atoms%astruct%cell_dim(1)
     if (atoms%astruct%geocode == 'P') factor(2) = 1._gp / atoms%astruct%cell_dim(2)
     if (atoms%astruct%geocode == 'P' .or. atoms%astruct%geocode == 'S') factor(3) = 1._gp / atoms%astruct%cell_dim(3)
  end if

  do iat=1,atoms%astruct%nat
     name=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     call write_extra_info(extra,atoms%astruct%input_polarization(iat),atoms%astruct%ifrztyp(iat))     

     write(iunit,'(3(1x,1pe24.17),2x,a2,2x,a50)') (rxyz(j,iat)*factor(j),j=1,3),symbol,extra
  end do


END SUBROUTINE wtascii

subroutine wtascii_forces(iunit,fxyz,at)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=1) :: endline
  
  if (at%astruct%nat ==0) return
  !write the first position
  iat=1
  if (at%astruct%nat==iat) then
     endline=']'
  else
     endline=char(92)
  end if

  write(iunit, "(A,3(1pe25.17,A),a)") "#metaData: forces=[",(fxyz(j,iat), ";",j=1,3),' '//endline
  !then the rest until the second-last
  do iat=2,at%astruct%nat
     if (at%astruct%nat==iat) then
        endline=']'
     else
        endline=char(92)
     end if
     write(iunit, "(A,3(1pe25.17,A),a)") "# ",(fxyz(j,iat), ";",j=1,3),' '//endline
  end do
end subroutine wtascii_forces

!>Write the extra info necessary for the output file
subroutine write_extra_info(extra,natpol,ifrztyp)
  use module_base
  implicit none 
  integer, intent(in) :: natpol,ifrztyp
  character(len=50), intent(out) :: extra
  !local variables
  character(len=4) :: frzchain
  integer :: ispol,ichg

  call charge_and_spol(natpol,ichg,ispol)

  call frozen_itof(ifrztyp,frzchain)
  
  !takes into account the blocked atoms and the input polarisation
  if (ispol == 0 .and. ichg == 0 ) then
     write(extra,'(2x,a4)')frzchain
  else if (ispol /= 0 .and. ichg == 0) then
     write(extra,'(i7,2x,a4)')ispol,frzchain
  else if (ichg /= 0) then
     write(extra,'(2(i7),2x,a4)')ispol,ichg,frzchain
  else
     write(extra,'(2x,a4)') ''
  end if
  
END SUBROUTINE write_extra_info


subroutine frozen_itof(ifrztyp,frzchain)
  implicit none
  integer, intent(in) :: ifrztyp
  character(len=4), intent(out) :: frzchain

  if (ifrztyp == 0) then
     frzchain='    '
  else if (ifrztyp == 1) then
     frzchain='   f'
  else if (ifrztyp == 2) then
     frzchain='  fy'
  else if (ifrztyp == 3) then
     frzchain=' fxz'
  end if
        
END SUBROUTINE frozen_itof

!>Write yaml atomic file.
subroutine wtyaml(iunit,energy,rxyz,atoms,comment,wrtforces,forces)
  use module_base
  use module_types
  use yaml_output
  implicit none
  logical, intent(in) :: wrtforces
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz,forces
  !local variables
  logical :: reduced
  integer :: iunit_def,iostat,ierr,iat,ichg,ispol
  real(gp) :: factor
  real(gp) :: xred(3)
  
  iostat=-1 !no changement of the default stream
  !associate iunit with a yaml_stream (do not crash if already associated)
  !first get the default stream
  call yaml_get_default_stream(iunit_def)
  if (iunit_def /= iunit) then
     call yaml_set_stream(unit=iunit,tabbing=0,record_length=100,istat=iostat)
     if (iostat /=0) then
        call yaml_set_default_stream(iunit,ierr)
     end if
     !if the stream was not already present just set back the default to iunit_def
  end if
  !start the writing of the file
  call yaml_new_document(unit=iunit)
  ! Possible comment.
  if (len(trim(comment)) > 0) then
     call yaml_map("Comment", trim(comment(1:min(85, len(trim(comment))))))
  end if
  !cell information
  call yaml_open_map('Cell')
  factor=1.0_gp
  Cell_Units: select case(trim(atoms%astruct%units))
  case('angstroem','angstroemd0')
     call yaml_map('Units','angstroem')
     factor=Bohr_Ang
  case('atomic','atomicd0','bohr','bohrd0','reduced')
     call yaml_map('Units','bohr')
     factor=1.0_gp
  end select Cell_Units
  BC :select case(atoms%astruct%geocode)
  case('F')
     call yaml_map('BC','free')
  case('S')
     call yaml_map('BC','surface')
     call yaml_open_sequence('acell',flow=.true.)
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(1)*factor)) !x
       call yaml_sequence('.inf')             !y
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(3)*factor)) !z
     call yaml_close_sequence()
     !angdeg to be added
  case('W')
     call yaml_map('BC','wire')
     call yaml_open_sequence('acell',flow=.true.)
       call yaml_sequence('.inf')             !x
       call yaml_sequence('.inf')             !y
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(3)*factor)) !z
     call yaml_close_sequence()
  case('P')
     call yaml_map('BC','periodic')
     call yaml_map('acell',(/atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,&
          atoms%astruct%cell_dim(3)*factor/))
     !angdeg to be added
  end select BC
  call yaml_close_map() !cell

  call yaml_open_map('Positions')
  reduced=.false.
  Pos_Units: select case(trim(atoms%astruct%units))
  case('angstroem','angstroemd0')
     call yaml_map('Units','angstroem')
  case('atomic','atomicd0','bohr','bohrd0')
     call yaml_map('Units','bohr')
  case('reduced')
     call yaml_map('Units','reduced')
     reduced=.true.
  end select Pos_Units
  call yaml_open_sequence('Values')
  do iat=1,atoms%astruct%nat
     call yaml_sequence(advance='no')
     if (extra_info(iat)) call yaml_open_map(flow=.true.)
     xred(1:3)=rxyz(1:3,iat)
     if (reduced) then
        if (atoms%astruct%geocode == 'P' .or. atoms%astruct%geocode =='S') xred(1)=rxyz(1,iat)/atoms%astruct%cell_dim(1)
        if (atoms%astruct%geocode == 'P') xred(2)=rxyz(2,iat)/atoms%astruct%cell_dim(2)
        if (atoms%astruct%geocode /='F') xred(3)=rxyz(3,iat)/atoms%astruct%cell_dim(3)
     end if
     call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),xred,fmt='(g25.17)')
     if (extra_info(iat)) then
        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
        if (ispol /=0) call yaml_map('IGSpin',ispol)
        if (ichg /=0) call yaml_map('IGChg',ichg)
        select case(atoms%astruct%ifrztyp(iat))
        case(1)
           call yaml_map('Frozen',.true.)
        case(2)
           call yaml_map('Frozen','fy')
        case(3)
           call yaml_map('Frozen','fxz')
        end select
        call yaml_close_map()
     end if
  end do
  call yaml_close_sequence() !values
  call yaml_close_map() !positions
  call yaml_open_map('Properties')
  call yaml_map('Timestamp',yaml_date_and_time_toa())
  if (energy /= 0. .and. energy /= UNINITIALIZED(energy)) then
     call yaml_map("Energy (Ha)", energy)
  end if
  call yaml_close_map() !properties
  if (wrtforces) then
     call yaml_open_map('Forces')
     call yaml_map('Units','Ha/Bohr')
     call yaml_open_sequence('Values')
     do iat=1,atoms%astruct%nat
        call yaml_sequence(advance='no')
        call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),forces(:,iat),fmt='(g25.17)')
     end do
     call yaml_close_sequence() !values
     call yaml_close_map() !forces
  end if

  !restore the default stream
  if (iostat==0) then
     call yaml_set_default_stream(iunit_def,ierr)
  end if

contains

  function extra_info(iat)
    implicit none
    integer, intent(in) :: iat
    logical extra_info
    extra_info=atoms%astruct%input_polarization(iat) /=100 .or. atoms%astruct%ifrztyp(iat)/=0
  end function extra_info

end subroutine wtyaml


!>Calculate the charge and the spin polarisation to be placed on a given atom
!!   RULE: natpol = c*1000 + sgn(c)*100 + s: charged and polarised atom (charge c, polarisation s)
subroutine charge_and_spol(natpol,nchrg,nspol)
  implicit none
  integer, intent(in) :: natpol
  integer, intent(out) :: nchrg,nspol
  !local variables
  integer :: nsgn

  nchrg=natpol/1000
  if (nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if

  nspol=natpol-1000*nchrg-nsgn*100

END SUBROUTINE charge_and_spol



! Init routine for bindings
!> Allocate a new atoms_data type, for bindings.
subroutine atoms_new(atoms, sym)
  use module_types
  implicit none
  type(atoms_data), pointer :: atoms
  type(symmetry_data), pointer :: sym
  
  allocate(atoms)
  atoms%astruct%geocode = "F"
  atoms%astruct%units = "bohr"
  atoms%astruct%inputfile_format = "none"
  atoms%astruct%nat = -1
  atoms%astruct%ntypes = -1
  atoms%astruct%sym%symObj = -1
  nullify(atoms%astruct%sym%irrzon)
  nullify(atoms%astruct%sym%phnons)
  sym => atoms%astruct%sym
  nullify(atoms%nlccpar)
  nullify(atoms%paw_l)
  nullify(atoms%paw_NofL)
  nullify(atoms%paw_nofchannels)
  nullify(atoms%paw_nofgaussians)
  nullify(atoms%paw_Greal)
  nullify(atoms%paw_Gimag)
  nullify(atoms%paw_Gcoeffs)
  nullify(atoms%paw_H_matrices)
  nullify(atoms%paw_S_matrices)
  nullify(atoms%paw_Sm1_matrices)
END SUBROUTINE atoms_new
subroutine atoms_set_from_file(lstat, atoms, rxyz, filename, ln)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   logical, intent(out) :: lstat
   type(atoms_data), intent(inout) :: atoms
   integer, intent(in) :: ln
   character, intent(in) :: filename(ln)
   real(gp), dimension(:,:), pointer :: rxyz

   integer :: status, i
   character(len = 1024) :: filename_
   character(len=*), parameter :: subname='atoms_set_from_file'

   write(filename_, "(A)") " "
   do i = 1, ln
      write(filename_(i:i), "(A1)") filename(i)
   end do

   lstat = .true.
   call read_atomic_file(trim(filename_), 0, atoms%astruct, status)
   rxyz=>atoms%astruct%rxyz
   call allocate_atoms_nat(atoms, subname)
   call allocate_atoms_ntypes(atoms, subname)
   lstat = (status == 0)
 END SUBROUTINE atoms_set_from_file
!> Deallocate a new atoms_data type, for bindings.
subroutine atoms_empty(atoms)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms

  call deallocate_atoms(atoms, "atoms_free")
END SUBROUTINE atoms_empty
subroutine atoms_free(atoms)
  use module_types
  implicit none
  type(atoms_data), pointer :: atoms
  
  call deallocate_atoms(atoms, "atoms_free")
  deallocate(atoms)
END SUBROUTINE atoms_free

! Set routines for bindings
subroutine atoms_read_variables(atoms, nspin, occup, ln)
  use module_types
  use memory_profiling
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: nspin, ln
  character, intent(in) :: occup(ln)

  integer :: i
  character(len = 1024) :: filename_

  write(filename_, "(A)") " "
  do i = 1, ln
     write(filename_(i:i), "(A1)") occup(i)
  end do

  call read_atomic_variables(atoms, trim(filename_), nspin)
END SUBROUTINE atoms_read_variables
subroutine atoms_set_n_atoms(atoms, rxyz, nat)
  use module_types
  use memory_profiling
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz
  integer, intent(in) :: nat

  integer :: i, i_stat

  call allocate_astruct_nat(atoms%astruct, nat, "atoms_set_n_atoms")
  atoms%astruct%iatype = (/ (i, i=1,nat) /)
  allocate(rxyz(3, atoms%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',"atoms_set_n_atoms")
  call allocate_atoms_nat(atoms, "atoms_set_n_atoms")
END SUBROUTINE atoms_set_n_atoms
subroutine atoms_set_n_types(atoms, ntypes)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ntypes

  call allocate_astruct_ntypes(atoms%astruct, ntypes, "atoms_set_n_types")
  call allocate_atoms_ntypes(atoms, "atoms_set_n_types")
END SUBROUTINE atoms_set_n_types
subroutine atoms_set_name(atoms, ityp, name)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ityp
  character, intent(in) :: name(20)

  write(atoms%astruct%atomnames(ityp), "(20A1)") name
END SUBROUTINE atoms_set_name
subroutine atoms_sync(atoms, alat1, alat2, alat3, geocode, format, units)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  real(gp), intent(in) :: alat1, alat2, alat3
  character, intent(in) :: geocode(1)
  character, intent(in) :: format(5)
  character, intent(in) :: units(20)

  atoms%astruct%cell_dim(1) = alat1
  atoms%astruct%cell_dim(2) = alat2
  atoms%astruct%cell_dim(3) = alat3
  atoms%astruct%geocode = geocode(1)
  write(atoms%astruct%inputfile_format, "(5A1)") format
  write(atoms%astruct%units, "(20A1)") units
END SUBROUTINE atoms_sync

! Accessors for bindings.
subroutine atoms_copy_nat(atoms, nat)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: nat

  nat = atoms%astruct%nat
END SUBROUTINE atoms_copy_nat
subroutine atoms_copy_ntypes(atoms, ntypes)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: ntypes

  ntypes = atoms%astruct%ntypes
END SUBROUTINE atoms_copy_ntypes
subroutine atoms_get_iatype(atoms, iatype)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: iatype

  iatype => atoms%astruct%iatype
END SUBROUTINE atoms_get_iatype
subroutine atoms_get_iasctype(atoms, iasctype)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: iasctype

  iasctype => atoms%iasctype
END SUBROUTINE atoms_get_iasctype
subroutine atoms_get_natpol(atoms, natpol)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: natpol

  natpol => atoms%astruct%input_polarization
END SUBROUTINE atoms_get_natpol
subroutine atoms_get_ifrztyp(atoms, ifrztyp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: ifrztyp

  ifrztyp => atoms%astruct%ifrztyp
END SUBROUTINE atoms_get_ifrztyp
subroutine atoms_get_nelpsp(atoms, nelpsp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nelpsp

  nelpsp => atoms%nelpsp
END SUBROUTINE atoms_get_nelpsp
subroutine atoms_get_npspcode(atoms, npspcode)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: npspcode

  npspcode => atoms%npspcode
END SUBROUTINE atoms_get_npspcode
subroutine atoms_get_nzatom(atoms, nzatom)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nzatom

  nzatom => atoms%nzatom
END SUBROUTINE atoms_get_nzatom
subroutine atoms_get_nlcc_ngv(atoms, nlcc_ngv)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nlcc_ngv

  nlcc_ngv => atoms%nlcc_ngv
END SUBROUTINE atoms_get_nlcc_ngv
subroutine atoms_get_nlcc_ngc(atoms, nlcc_ngc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: nlcc_ngc

  nlcc_ngc => atoms%nlcc_ngc
END SUBROUTINE atoms_get_nlcc_ngc
subroutine atoms_get_ixcpsp(atoms, ixcpsp)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: ixcpsp

  ixcpsp => atoms%ixcpsp
END SUBROUTINE atoms_get_ixcpsp
subroutine atoms_get_amu(atoms, amu)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:), pointer :: amu

  amu => atoms%amu
END SUBROUTINE atoms_get_amu


subroutine atoms_get_aocc(atoms, aocc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: aocc

  aocc => atoms%aocc
END SUBROUTINE atoms_get_aocc


!> get radii_cf values
subroutine atoms_get_radii_cf(atoms, radii_cf)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: radii_cf

  radii_cf => atoms%radii_cf
END SUBROUTINE atoms_get_radii_cf


subroutine atoms_get_psppar(atoms, psppar)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:,:), pointer :: psppar

  psppar => atoms%psppar
END SUBROUTINE atoms_get_psppar


subroutine atoms_get_nlccpar(atoms, nlccpar)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: nlccpar

  nlccpar => atoms%nlccpar
END SUBROUTINE atoms_get_nlccpar


subroutine atoms_get_ig_nlccpar(atoms, ig_nlccpar)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: ig_nlccpar

  ig_nlccpar => atoms%ig_nlccpar
END SUBROUTINE atoms_get_ig_nlccpar


subroutine atoms_copy_geometry_data(atoms, geocode, format, units)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  character(len = 1), intent(out) :: geocode
  character(len = 5), intent(out) :: format
  character(len = 20), intent(out) :: units

  write(geocode, "(A1)") atoms%astruct%geocode
  write(format,  "(A5)") atoms%astruct%inputfile_format
  write(units,  "(A20)") atoms%astruct%units
END SUBROUTINE atoms_copy_geometry_data


subroutine atoms_copy_psp_data(atoms, natsc, donlcc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: natsc
  logical, intent(out) :: donlcc

  natsc = atoms%natsc
  donlcc = atoms%donlcc
END SUBROUTINE atoms_copy_psp_data


subroutine atoms_copy_name(atoms, ityp, name, ln)
  use module_types
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: ityp
  character(len=1), dimension(20), intent(out) :: name
!  character(len=*), intent(out) :: name
  integer, intent(out) :: ln
  !Local variables 
  integer :: i,lname

  lname = len(name)
  ln=min(len(trim(atoms%astruct%atomnames(ityp))),20)
  !print *,'lnt2',lnt
  do i = 1, ln, 1
     name(i:i) = atoms%astruct%atomnames(ityp)(i:i)
  end do
  do i = ln + 1, lname, 1
     name(i:i) = ' '
  end do
END SUBROUTINE atoms_copy_name


subroutine atoms_copy_alat(atoms, alat1, alat2, alat3)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: alat1, alat2, alat3

  alat1 = atoms%astruct%cell_dim(1)
  alat2 = atoms%astruct%cell_dim(2)
  alat3 = atoms%astruct%cell_dim(3)
END SUBROUTINE atoms_copy_alat
subroutine atoms_write(atoms, filename, filelen, rxyz, forces, energy, comment, ln)
  use module_types
  implicit none
  integer, intent(in) :: ln, filelen
  character, intent(in) :: comment(ln)
  character, intent(in) :: filename(filelen)
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(:,:), pointer :: forces

  integer :: iunit, i
  character(len = 1024) :: comment_, filename_

  write(filename_, "(A)") " "
  do i = 1, filelen
     write(filename_(i:i), "(A1)") filename(i)
  end do
  if (trim(filename_) == "stdout") then
     iunit = 6
  else
     open(unit=9,file=trim(filename_)//'.'//trim(atoms%astruct%inputfile_format))
     iunit = 9
  end if
  write(comment_, "(A)") " "
  do i = 1, ln
     write(comment_(i:i), "(A1)") comment(i)
  end do

  if (trim(atoms%astruct%inputfile_format) == "xyz") then
     call wtxyz(iunit,energy,rxyz,atoms,comment_)
     if (associated(forces)) call wtxyz_forces(iunit,forces,atoms)
  else if (trim(atoms%astruct%inputfile_format) == "ascii") then
     call wtascii(iunit,energy,rxyz,atoms,comment_)
     if (associated(forces)) call wtascii_forces(iunit,forces,atoms)
  else
     write(*,*) "Error, unknown file format."
     stop
  end if

  if (trim(filename_) /= "stdout") then
     close(unit=9)
  end if
END SUBROUTINE atoms_write

subroutine symmetry_set_irreductible_zone(sym, geocode, n1i, n2i, n3i, nspin)
  use module_base
  use module_types
  use m_ab6_kpoints
  use m_ab6_symmetry
  implicit none
  type(symmetry_data), intent(inout) :: sym
  integer, intent(in) :: n1i, n2i, n3i, nspin
  character, intent(in) :: geocode

  character(len = *), parameter :: subname = "symmetry_set_irreductible_zone"
  integer :: i_stat, nsym, i_all, i_third
  integer, dimension(:,:,:), allocatable :: irrzon
  real(dp), dimension(:,:,:), allocatable :: phnons

  if (associated(sym%irrzon)) then
     i_all=-product(shape(sym%irrzon))*kind(sym%irrzon)
     deallocate(sym%irrzon,stat=i_stat)
     call memocc(i_stat,i_all,'irrzon',subname)
     nullify(sym%irrzon)
  end if

  if (associated(sym%phnons)) then
     i_all=-product(shape(sym%phnons))*kind(sym%phnons)
     deallocate(sym%phnons,stat=i_stat)
     call memocc(i_stat,i_all,'phnons',subname)
     nullify(sym%phnons)
  end if

  if (sym%symObj >= 0) then
     call symmetry_get_n_sym(sym%symObj, nsym, i_stat)
     if (nsym > 1) then
        ! Current third dimension is set to 1 always
        ! since nspin == nsppol always in BigDFT
        i_third = 1
        if (geocode == "S") i_third = n2i
        allocate(sym%irrzon(n1i*(n2i - i_third + 1)*n3i,2,i_third+ndebug),stat=i_stat)
        call memocc(i_stat,sym%irrzon,'irrzon',subname)
        allocate(sym%phnons(2,n1i*(n2i - i_third + 1)*n3i,i_third+ndebug),stat=i_stat)
        call memocc(i_stat,sym%phnons,'phnons',subname)
        if (geocode /= "S") then
           call kpoints_get_irreductible_zone(sym%irrzon, sym%phnons, &
                &   n1i, n2i, n3i, nspin, nspin, sym%symObj, i_stat)
        else
           allocate(irrzon(n1i*n3i,2,1+ndebug),stat=i_stat)
           call memocc(i_stat,irrzon,'irrzon',subname)
           allocate(phnons(2,n1i*n3i,1+ndebug),stat=i_stat)
           call memocc(i_stat,phnons,'phnons',subname)
           do i_third = 1, n2i, 1
              call kpoints_get_irreductible_zone(irrzon, phnons, n1i, 1, n3i, &
                   & nspin, nspin, sym%symObj, i_stat)
              sym%irrzon(:,:,i_third:i_third) = irrzon
              call dcopy(2*n1i*n3i, phnons, 1, sym%phnons(1,1,i_third), 1)
           end do
           i_all=-product(shape(irrzon))*kind(irrzon)
           deallocate(irrzon,stat=i_stat)
           call memocc(i_stat,i_all,'irrzon',subname)
           i_all=-product(shape(phnons))*kind(phnons)
           deallocate(phnons,stat=i_stat)
           call memocc(i_stat,i_all,'phnons',subname)
        end if
     end if
  end if

  if (.not. associated(sym%irrzon)) then
     ! Allocate anyway to small size otherwise the bounds check does not pass.
     allocate(sym%irrzon(1,2,1+ndebug),stat=i_stat)
     call memocc(i_stat,sym%irrzon,'irrzon',subname)
     allocate(sym%phnons(2,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,sym%phnons,'phnons',subname)
  end if
END SUBROUTINE symmetry_set_irreductible_zone

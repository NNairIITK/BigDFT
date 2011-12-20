!> Deallocate the structure atoms_data.
subroutine deallocate_atoms(atoms,subname) 
  use module_base
  use module_types
  use m_ab6_symmetry
  implicit none
  character(len=*), intent(in) :: subname
  type(atoms_data), intent(inout) :: atoms
  !local variables
  integer :: i_stat, i_all

  ! Deallocations for the geometry part.
  if (atoms%nat > 0) then
     i_all=-product(shape(atoms%ifrztyp))*kind(atoms%ifrztyp)
     deallocate(atoms%ifrztyp,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%ifrztyp',subname)
     i_all=-product(shape(atoms%iatype))*kind(atoms%iatype)
     deallocate(atoms%iatype,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%iatype',subname)
     i_all=-product(shape(atoms%natpol))*kind(atoms%natpol)
     deallocate(atoms%natpol,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%natpol',subname)
     i_all=-product(shape(atoms%amu))*kind(atoms%amu)
     deallocate(atoms%amu,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%amu',subname)
  end if
  if (atoms%ntypes > 0) then
     i_all=-product(shape(atoms%atomnames))*kind(atoms%atomnames)
     deallocate(atoms%atomnames,stat=i_stat)
     call memocc(i_stat,i_all,'atoms%atomnames',subname)
  end if
  if (atoms%symObj >= 0) then
     call symmetry_free(atoms%symObj)
  end if

  ! Deallocations related to pseudos.
  if (atoms%ntypes > 0) then
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

subroutine allocate_atoms_nat(atoms, nat, subname)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: nat
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat

  atoms%nat = nat

  ! Allocate geometry related stuff.
  allocate(atoms%iatype(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%iatype,'atoms%iatype',subname)
  allocate(atoms%ifrztyp(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%ifrztyp,'atoms%ifrztyp',subname)
  allocate(atoms%natpol(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%natpol,'atoms%natpol',subname)
  allocate(atoms%amu(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%amu,'atoms%amu',subname)

  !this array is useful for frozen atoms, no atom is frozen by default
  atoms%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  atoms%natpol(:)=100
END SUBROUTINE allocate_atoms_nat

subroutine allocate_atoms_ntypes(atoms, ntypes, subname)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ntypes
  character(len = *), intent(in) :: subname
  !local variables
  integer :: i_stat
  integer, parameter :: nelecmax=32

  atoms%ntypes = ntypes

  ! Allocate geometry related stuff.
  allocate(atoms%atomnames(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%atomnames,'atoms%atomnames',subname)

  ! Allocate pseudo related stuff.
  ! store PSP parameters, modified to accept both GTH and HGHs pseudopotential types
  allocate(atoms%psppar(0:4,0:6,atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%psppar,'atoms%psppar',subname)
  allocate(atoms%nelpsp(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nelpsp,'atoms%nelpsp',subname)
  allocate(atoms%npspcode(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%npspcode,'atoms%npspcode',subname)
  allocate(atoms%nzatom(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nzatom,'atoms%nzatom',subname)
  allocate(atoms%ixcpsp(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%ixcpsp,'atoms%ixcpsp',subname)
  allocate(atoms%radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%radii_cf,'atoms%radii_cf',subname)
  ! parameters for NLCC
  allocate(atoms%nlcc_ngv(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngv,'atoms%nlcc_ngv',subname)
  allocate(atoms%nlcc_ngc(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngc,'atoms%nlcc_ngc',subname)
  ! semicores useful only for the input guess
  allocate(atoms%iasctype(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%iasctype,'atoms%iasctype',subname)
  allocate(atoms%aocc(nelecmax,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%aocc,'atoms%aocc',subname)
END SUBROUTINE allocate_atoms_ntypes





!> Read atomic positions
subroutine read_xyz_positions(iproc,ifile,atoms,rxyz,getLine)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz
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
  read(line,*, iostat = ierrsfx) iat,atoms%units
  if (ierrsfx /= 0) then
     read(line,*, iostat = ierrsfx) iat
     write(atoms%units, "(A)") "bohr"
  end if

  allocate(rxyz(3,iat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)
  call allocate_atoms_nat(atoms, iat, subname)

  !controls if the positions are provided with machine precision
  if (atoms%units == 'angstroemd0' .or. atoms%units== 'atomicd0' .or. &
       atoms%units== 'bohrd0' .or. atoms%units=='reduced') then
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
!!!     atoms%geocode='F'
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
        atoms%geocode='P'
     else if (trim(tatonam)=='surface') then 
        atoms%geocode='S'
        atoms%alat2=0.0_gp
     else !otherwise free bc
        atoms%geocode='F'
        atoms%alat1=0.0_gp
        atoms%alat2=0.0_gp
        atoms%alat3=0.0_gp
     end if
     if (.not. lpsdbl) then
        alat1d0=real(alat1,gp)
        alat2d0=real(alat2,gp)
        alat3d0=real(alat3,gp)
     end if
  else
     atoms%geocode='F'
     alat1d0=0.0_gp
     alat2d0=0.0_gp
     alat3d0=0.0_gp
  end if
!!!  end if

  !reduced coordinates are possible only with periodic units
  if (atoms%units == 'reduced' .and. atoms%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

  !convert the values of the cell sizes in bohr
  if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     atoms%alat1=alat1d0/bohr2ang
     atoms%alat2=alat2d0/bohr2ang
     atoms%alat3=alat3d0/bohr2ang
  else if  (atoms%units=='atomic' .or. atoms%units=='bohr'  .or.&
       atoms%units== 'atomicd0' .or. atoms%units== 'bohrd0') then
     atoms%alat1=alat1d0
     atoms%alat2=alat2d0
     atoms%alat3=alat3d0
  else if (atoms%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     atoms%alat1=alat1d0
     atoms%alat2=alat2d0
     atoms%alat3=alat3d0
  else
     write(*,*) 'length units in input file unrecognized'
     write(*,*) 'recognized units are angstroem or atomic = bohr'
     stop 
  endif

  ntyp=0
  do iat=1,atoms%nat
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
     call parse_extra_info(iat,extra,atoms)

     tatonam=trim(symbol)
!!!     end if
     if (lpsdbl) then
        rxyz(1,iat)=rxd0
        rxyz(2,iat)=ryd0
        rxyz(3,iat)=rzd0
     else
        rxyz(1,iat)=real(rx,gp)
        rxyz(2,iat)=real(ry,gp)
        rxyz(3,iat)=real(rz,gp)
     end if

     

     if (atoms%units == 'reduced') then !add treatment for reduced coordinates
        rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
        if (atoms%geocode == 'P') rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
        rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
     else if (atoms%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(2,iat)=modulo(rxyz(2,iat),alat2d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
     else if (atoms%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
     end if
 
     do ityp=1,ntyp
        if (tatonam == atomnames(ityp)) then
           atoms%iatype(iat)=ityp
           goto 200
        endif
     enddo
     ntyp=ntyp+1
     if (ntyp > 100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     atoms%iatype(iat)=ntyp
200  continue

     if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/bohr2ang
        enddo
     else if (atoms%units == 'reduced') then 
        rxyz(1,iat)=rxyz(1,iat)*atoms%alat1
        if (atoms%geocode == 'P') rxyz(2,iat)=rxyz(2,iat)*atoms%alat2
        rxyz(3,iat)=rxyz(3,iat)*atoms%alat3
     endif
  enddo

  !now that ntypes is determined allocate atoms%atomnames and copy the values
  call allocate_atoms_ntypes(atoms, ntyp, subname)
  atoms%atomnames(1:atoms%ntypes)=atomnames(1:atoms%ntypes)
END SUBROUTINE read_xyz_positions

!> Read atomic positions of ascii files.
subroutine read_ascii_positions(iproc,ifile,atoms,rxyz,getline)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz
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
  logical :: lpsdbl, reduced, eof
  integer :: iat,ntyp,ityp,i,i_stat,j,nlines
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
        stop 
     end if
  end do
  nlines = nlines - 1

  if (nlines < 4) then
     if (iproc==0) write(*,*) 'Error in ASCII file format, file has less than 4 lines.'
     stop 
  end if

  ! Try to determine the number atoms and the keywords.
  write(atoms%units, "(A)") "bohr"
  reduced = .false.
  atoms%geocode = 'P'
  iat     = 0
  do i = 4, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (line(1:1) /= '#' .and. line(1:1) /= '!' .and. len(trim(line)) /= 0) then
        iat = iat + 1
     else if (line(1:8) == "#keyword" .or. line(1:8) == "!keyword") then
        if (index(line, 'bohr') > 0)        write(atoms%units, "(A)") "bohr"
        if (index(line, 'bohrd0') > 0)      write(atoms%units, "(A)") "bohrd0"
        if (index(line, 'atomic') > 0)      write(atoms%units, "(A)") "atomicd0"
        if (index(line, 'angstroem') > 0)   write(atoms%units, "(A)") "angstroem"
        if (index(line, 'angstroemd0') > 0) write(atoms%units, "(A)") "angstroemd0"
        if (index(line, 'reduced') > 0)     reduced = .true.
        if (index(line, 'periodic') > 0) atoms%geocode = 'P'
        if (index(line, 'surface') > 0)  atoms%geocode = 'S'
        if (index(line, 'freeBC') > 0)   atoms%geocode = 'F'
     end if
  end do

  allocate(rxyz(3,iat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz,'rxyz',subname)
  call allocate_atoms_nat(atoms, iat, subname)

  !controls if the positions are provided within machine precision
  if (index(atoms%units, 'd0') > 0 .or. reduced) then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  ! Read the box definition
  atoms%alat1 = 0.0_gp
  atoms%alat2 = 0.0_gp
  atoms%alat3 = 0.0_gp
  if (lpsdbl) then
     read(lines(2),*) alat1d0,alat2d0,alat3d0
     read(lines(3),*) alat4d0,alat5d0,alat6d0
     if (alat2d0 /= 0.d0 .or. alat4d0 /= 0.d0 .or. alat5d0 /= 0.d0) then
        !if (iproc==0) 
        write(*,*) 'Only orthorombic boxes are possible.'
        stop 
     end if
     atoms%alat1 = real(alat1d0,gp)
     atoms%alat2 = real(alat3d0,gp)
     atoms%alat3 = real(alat6d0,gp)
  else
     read(lines(2),*) alat1,alat2,alat3
     read(lines(3),*) alat4,alat5,alat6
     if (alat2 /= 0. .or. alat4 /= 0. .or. alat5 /= 0.) then
        !if (iproc==0) 
           write(*,*) 'Only orthorombic boxes are possible.'
        !if (iproc==0) 
           write(*,*) ' but alat2, alat4 and alat5 = ', alat2, alat4, alat5
        stop 
     end if
     atoms%alat1 = real(alat1,gp)
     atoms%alat2 = real(alat3,gp)
     atoms%alat3 = real(alat6,gp)
  end if
  
  !Convert the values of the cell sizes in bohr
  if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     atoms%alat1 = atoms%alat1 / bohr2ang
     atoms%alat2 = atoms%alat2 / bohr2ang
     atoms%alat3 = atoms%alat3 / bohr2ang
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
        call parse_extra_info(iat,extra,atoms)

        tatonam=trim(symbol)

        if (lpsdbl) then
           rxyz(1,iat)=rxd0
           rxyz(2,iat)=ryd0
           rxyz(3,iat)=rzd0
        else
           rxyz(1,iat)=real(rx,gp)
           rxyz(2,iat)=real(ry,gp)
           rxyz(3,iat)=real(rz,gp)
        end if

        if (reduced) then !add treatment for reduced coordinates
           rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
           rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
           rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
        else if (atoms%geocode == 'P') then
           rxyz(1,iat)=modulo(rxyz(1,iat),atoms%alat1)
           rxyz(2,iat)=modulo(rxyz(2,iat),atoms%alat2)
           rxyz(3,iat)=modulo(rxyz(3,iat),atoms%alat3)
        else if (atoms%geocode == 'S') then
           rxyz(1,iat)=modulo(rxyz(1,iat),atoms%alat1)
           rxyz(3,iat)=modulo(rxyz(3,iat),atoms%alat3)
        end if

        do ityp=1,ntyp
           if (tatonam == atomnames(ityp)) then
              atoms%iatype(iat)=ityp
              goto 200
           endif
        enddo
        ntyp=ntyp+1
        if (ntyp > 100) stop 'more than 100 atomnames not permitted'
        atomnames(ityp)=tatonam
        atoms%iatype(iat)=ntyp
200     continue

        if (reduced) then
           rxyz(1,iat)=rxyz(1,iat)*atoms%alat1
           rxyz(2,iat)=rxyz(2,iat)*atoms%alat2
           rxyz(3,iat)=rxyz(3,iat)*atoms%alat3
        else if (atoms%units=='angstroem' .or. atoms%units=='angstroemd0') then
           ! if Angstroem convert to Bohr
           do j=1,3 
              rxyz(j,iat)=rxyz(j,iat) / bohr2ang
           enddo
        endif
        iat = iat + 1
     end if
  enddo

  if (atoms%geocode == 'S') then
     atoms%alat2 = 0.0_gp
  else if (atoms%geocode == 'F') then
     atoms%alat1 = 0.0_gp
     atoms%alat2 = 0.0_gp
     atoms%alat3 = 0.0_gp
  end if

  !now that ntypes is determined copy the values
  call allocate_atoms_ntypes(atoms, ntyp, subname)
  atoms%atomnames(1:atoms%ntypes)=atomnames(1:atoms%ntypes)
END SUBROUTINE read_ascii_positions

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
subroutine parse_extra_info(iat,extra,atoms)
  use module_types
  implicit none
  !Arguments
  integer, intent(in) :: iat
  character(len=50), intent(in) :: extra
  type(atoms_data), intent(inout) :: atoms
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
  atoms%natpol(iat)=1000*nchrg+nsgn*100+nspol

  !print *,'natpol atomic',iat,atoms%natpol(iat),suffix

  !convert the suffix into ifrztyp
  call frozen_ftoi(suffix,atoms%ifrztyp(iat))

!!!  if (trim(suffix) == 'f') then
!!!     !the atom is considered as blocked
!!!     atoms%ifrztyp(iat)=1
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

!>Write xyz atomic file.
subroutine wtxyz(iunit,energy,rxyz,atoms,comment)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
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

  do iat=1,atoms%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%units) == 'angstroem' .or. trim(atoms%units) == 'angstroemd0') then
     factor=bohr2ang
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if

  write(iunit,'(i6,2x,a,2x,1pe24.17,2x,a)') atoms%nat,trim(units),energy,comment

  if (atoms%geocode == 'P') then
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          atoms%alat1*factor,atoms%alat2*factor,atoms%alat3*factor
  else if (atoms%geocode == 'S') then
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          atoms%alat1*factor,atoms%alat2*factor,atoms%alat3*factor
  else
     write(iunit,*)'free'
  end if
  do iat=1,atoms%nat
     name=trim(atoms%atomnames(atoms%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:min(len(name),5))
     end if

     call write_extra_info(extra,atoms%natpol(iat),atoms%ifrztyp(iat))

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
  real(gp), dimension(3,at%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=20) :: symbol
  character(len=10) :: name

  write(iunit,*)'forces'
  
  do iat=1,at%nat
     name=trim(at%atomnames(at%iatype(iat)))
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
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !local variables
  character(len=2) :: symbol
  character(len=50) :: extra
  character(len=10) :: name
  integer :: iat,j
  real(gp) :: xmax,ymax,zmax,factor

  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,atoms%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(atoms%units) == 'angstroem' .or. trim(atoms%units) == 'angstroemd0') then
     factor=bohr2ang
  else
     factor=1.0_gp
  end if

  write(iunit, "(A,A)") "# BigDFT file - ", trim(comment)
  write(iunit, "(3e24.17)") atoms%alat1*factor, 0.d0, atoms%alat2*factor
  write(iunit, "(3e24.17)") 0.d0,               0.d0, atoms%alat3*factor

  write(iunit, "(A,A)") "#keyword: ", trim(atoms%units)
  if (atoms%geocode == 'P') write(iunit, "(A)") "#keyword: periodic"
  if (atoms%geocode == 'S') write(iunit, "(A)") "#keyword: surface"
  if (atoms%geocode == 'F') write(iunit, "(A)") "#keyword: freeBC"
  if (energy /= 0.d0) then
     write(iunit, "(A,e24.17,A)") "#metaData: totalEnergy=", energy, "Ht"
  end if

  do iat=1,atoms%nat
     name=trim(atoms%atomnames(atoms%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     call write_extra_info(extra,atoms%natpol(iat),atoms%ifrztyp(iat))     

     write(iunit,'(3(1x,1pe24.17),2x,a2,2x,a50)') (rxyz(j,iat)*factor,j=1,3),symbol,extra
  end do


END SUBROUTINE wtascii

subroutine wtascii_forces(iunit,fxyz,at)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=1) :: endline
  
  if (at%nat ==0) return
  !write the first position
  iat=1
  if (at%nat==iat) then
     endline=']'
  else
     endline=char(92)
  end if

  write(iunit, "(A,3(1pe25.17,A),a)") "#metaData: forces=[",(fxyz(j,iat), ";",j=1,3),' '//endline
  !then the rest until the second-last
  do iat=2,at%nat
     if (at%nat==iat) then
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
subroutine atoms_new(atoms)
  use module_types
  implicit none
  type(atoms_data), pointer :: atoms
  
  allocate(atoms)
  atoms%geocode = "F"
  atoms%units = "bohr"
  atoms%format = "none"
  atoms%nat = -1
  atoms%ntypes = -1
  atoms%symObj = -1
END SUBROUTINE atoms_new
subroutine atoms_new_from_file(lstat, atoms, rxyz, filename, ln)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   logical, intent(out) :: lstat
   type(atoms_data), pointer :: atoms
   integer, intent(in) :: ln
   character(len = ln), intent(in) :: filename
   real(gp), dimension(:,:), pointer :: rxyz

   integer :: status

   lstat = .true.
   allocate(atoms)
   call read_atomic_file(filename, 0, atoms, rxyz, status)
   lstat = (status == 0)
END SUBROUTINE atoms_new_from_file
!> Deallocate a new atoms_data type, for bindings.
subroutine atoms_free(atoms)
  use module_types
  implicit none
  type(atoms_data), pointer :: atoms
  
  call deallocate_atoms(atoms, "atoms_free")
  deallocate(atoms)
END SUBROUTINE atoms_free

! Set routines for bindings
subroutine atoms_set_n_atoms(atoms, nat)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: nat

  integer :: i

  call allocate_atoms_nat(atoms, nat, "atoms_set_n_atoms")
  atoms%iatype = (/ (i, i=1,nat) /)
END SUBROUTINE atoms_set_n_atoms
subroutine atoms_set_n_types(atoms, ntypes)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ntypes

  call allocate_atoms_ntypes(atoms, ntypes, "atoms_set_n_types")
END SUBROUTINE atoms_set_n_types

! Accessors for bindings.
subroutine atoms_get_nat(atoms, nat)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: nat

  nat = atoms%nat
END SUBROUTINE atoms_get_nat
subroutine atoms_get_ntypes(atoms, ntypes)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: ntypes

  ntypes = atoms%ntypes
END SUBROUTINE atoms_get_ntypes
subroutine atoms_get_iatype(atoms, iatype)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: iatype

  iatype => atoms%iatype
END SUBROUTINE atoms_get_iatype
subroutine atoms_get_geocode(atoms, geocode)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  character(len = 1), intent(out) :: geocode

  geocode = atoms%geocode
END SUBROUTINE atoms_get_geocode
subroutine atoms_get_name(atoms, ityp, name, ln)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: ityp
  character(len = 20), intent(out) :: name
  integer, intent(out) :: ln

  name = atoms%atomnames(ityp)
  ln = len(trim(name))
END SUBROUTINE atoms_get_name
subroutine atoms_get_alat(atoms, alat1, alat2, alat3)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: alat1, alat2, alat3

  alat1 = atoms%alat1
  alat2 = atoms%alat2
  alat3 = atoms%alat3
END SUBROUTINE atoms_get_alat

!> @file
!!  Private Routines for the setting and creation of the astruct structure
!!  included in the module module_atoms
!! @author
!!    Copyright (C) 2007-2013 BigDFT group (TD,LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Read atomic positions from xyz file and create astruct structure from it
subroutine read_xyz_positions(ifile,filename,astruct,comment,energy,fxyz,getLine)
  use module_defs, only: gp,UNINITIALIZED,Bohr_Ang, BIGDFT_INPUT_VARIABLES_ERROR
  use dictionaries, only: f_err_raise, f_err_throw
  use module_base, only: ndebug,memocc
  use dynamic_memory
  implicit none
  !Arguments
  integer, intent(in) :: ifile
  character(len=*), intent(in) :: filename
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
  character(len=*), parameter :: subname='read_atomic_positions'
  character(len=20) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl, eof
  integer :: iat,ityp,ntyp,i,ierrsfx
  ! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
  ! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  character(len=20), dimension(100) :: atomnames

  call getLine(line, ifile, eof)
  if (f_err_raise(eof,"Unexpected end of file '"//trim(filename)//"'.",err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return
  !if (eof) then
  !   write(*,*) "Error: unexpected end of file."
  !   stop
  !end if
  energy = UNINITIALIZED(energy)
  read(line,*, iostat = ierrsfx) iat,astruct%units,energy,comment
  if (ierrsfx /= 0) then
     read(line,*, iostat = ierrsfx) iat,astruct%units,energy
     write(comment, "(A)") ""
     if (ierrsfx /= 0) then
        read(line,*, iostat = ierrsfx) iat,astruct%units
        if (ierrsfx /= 0) then
           read(line,*, iostat = ierrsfx) iat
           write(astruct%units, "(A)") "bohr"
        end if
     end if
  else
     i = index(line, trim(comment))
     write(comment, "(A)") line(i:)
  end if

  call astruct_set_n_atoms(astruct, iat)

  !controls if the positions are provided with machine precision
  if (astruct%units == 'angstroemd0' .or. astruct%units== 'atomicd0' .or. &
       astruct%units== 'bohrd0' .or. astruct%units=='reduced') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !read from positions of .xyz format, but accepts also the old .ascii format
  call getLine(line, ifile, eof)
  if (f_err_raise(eof,"Unexpected end of file '"//trim(filename)//"'.",err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return
  !if (eof) then
     !write(*,*) "Error: unexpected end of file."
     !stop
  !end if

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
  if (f_err_raise( (astruct%units == 'reduced' .and. astruct%geocode == 'F'), &
     & 'Reduced coordinates are not allowed with isolated BC', &
       err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return
  !if (astruct%units == 'reduced' .and. astruct%geocode == 'F') then
  !   if (iproc==0) write(*,'(1x,a)')&
  !        'ERROR: Reduced coordinates are not allowed with isolated BC'
  !end if

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
     call f_err_throw('Length units in input file unrecognized.' // &
          'Recognized units are angstroem or atomic = bohr',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
     return
     !write(*,*) 'length units in input file unrecognized'
     !write(*,*) 'recognized units are angstroem or atomic = bohr'
     !stop 
  endif

  ntyp=0
  do iat=1,astruct%nat
     !xyz input file, allow extra information
     call getLine(line, ifile, eof)
     if (f_err_raise(eof,"Unexpected end of file '"//trim(filename)//"'.",err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return

     !!if (lpsdbl) then
     !!   read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
     !!else
     !!   read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
     !!end if
     call check_line_integrity()
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
     if (f_err_raise(ntyp > 100, "More than 100 atomnames not permitted", err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return
     !if (ntyp > 100) stop 'more than 100 atomnames not permitted'
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
     fxyz = f_malloc_ptr((/ 3, iat /),id='fxyz')
     do iat=1,astruct%nat
        !xyz input file, allow extra information
        call getLine(line, ifile, eof)
        if (f_err_raise(eof,"Unexpected end of file '"//trim(filename)//"'.",err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return
        !if (eof) then
        !   write(*,*) "Error: unexpected end of file."
        !   stop
        !end if
        read(line,*,iostat=ierrsfx) symbol,fxyz(:,iat)
     end do
  end if
  !now that ntypes is determined allocate atoms%astruct%atomnames and copy the values
  call astruct_set_n_types(astruct, ntyp)

  astruct%atomnames(1:astruct%ntypes)=atomnames(1:astruct%ntypes)


contains

  !> stop the code and warns if the status of the line is not good
  subroutine check_line_integrity()
    use yaml_output, only: yaml_toa
    use dictionaries, only: f_err_raise
    implicit none


    if (lpsdbl) then
       read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0
    else
       read(line,*,iostat=ierrsfx)symbol,rx,ry,rz
    end if

    if (f_err_raise(ierrsfx/=0,'The line'//trim(yaml_toa(iat+2))//&
         ' of the atomic position is not valid, check if it is in DOS format!',&
         err_name='BIGDFT_LINALG_ERROR')) return

  end subroutine check_line_integrity

END SUBROUTINE read_xyz_positions


!> Read atomic positions of ascii files.
subroutine read_ascii_positions(ifile,filename,astruct,comment,energy,fxyz,getline)
  use module_base
  use dynamic_memory
  use yaml_output
  implicit none
  integer, intent(in) :: ifile
  character(len=*), intent(in) :: filename
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
     if (f_err_raise(nlines > 5000,  "Atomic input file '"//trim(filename)//"' too long (> 5000 lines).", &
        & err_id=BIGDFT_INPUT_VARIABLES_ERROR)) then
        astruct%nat = -1
        return
     end if
  end do
  nlines = nlines - 1

  if (f_err_raise(nlines < 4, "Error in ASCII file format, file '" // trim(filename) // "' has less than 4 lines.", &
     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)) then
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

  call astruct_set_n_atoms(astruct, iat)

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
     if (f_err_raise( (alat2d0 /= 0.d0 .or. alat4d0 /= 0.d0 .or. alat5d0 /= 0.d0), &
        & "File '" // trim(filename) // "': Only orthorombic boxes are possible.", err_id=BIGDFT_INPUT_VARIABLES_ERROR)) then
        astruct%nat = -1
        return
     end if
     astruct%cell_dim(1) = real(alat1d0,gp)
     astruct%cell_dim(2) = real(alat3d0,gp)
     astruct%cell_dim(3) = real(alat6d0,gp)
  else
     read(lines(2),*) alat1,alat2,alat3
     read(lines(3),*) alat4,alat5,alat6
     if (f_err_raise( (alat2 /= 0. .or. alat4 /= 0. .or. alat5 /= 0.), &
        & "File '" // trim(filename) // "': Only orthorombic boxes are possible but alat2, alat4 and alat5 = " // &
        & trim(yaml_toa( (/ alat2, alat4, alat5 /) )), err_id=BIGDFT_INPUT_VARIABLES_ERROR)) then
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
        if (f_err_raise(ntyp>100, "File '"//trim(filename)//"' more than 100 atomnames not permitted", &
           & err_id=BIGDFT_INPUT_VARIABLES_ERROR)) then
        !if (ntyp > 100) then
        !   write(*,*) 'more than 100 atomnames not permitted'
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
     fxyz = f_malloc_ptr((/ 3, astruct%nat /),id='fxyz')

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
  call astruct_set_n_types(astruct, ntyp)
  astruct%atomnames(1:astruct%ntypes)=atomnames(1:astruct%ntypes)
END SUBROUTINE read_ascii_positions


!> Local routines with explicit interface 
subroutine directGetLine(line, ifile, eof)
  !Arguments
  integer, intent(in) :: ifile
  character(len=150), intent(out) :: line
  logical, intent(out) :: eof
  !Local variables
  integer :: i_stat

  eof = .false.
  read(ifile,'(a150)', iostat = i_stat) line
  if (i_stat /= 0) eof = .true.
END SUBROUTINE directGetLine


subroutine archiveGetLine(line, ifile, eof)
  !Arguments
  integer, intent(in) :: ifile
  character(len=150), intent(out) :: line
  logical, intent(out) :: eof
  !Local variables
  integer :: i_stat
  external :: extractNextLine
  !The argument ifile is not used but it is used as argument routine
  !eof = .false.
  eof = (ifile /= ifile)
  call extractNextLine(line, i_stat)
  if (i_stat /= 0) eof = .true.
END SUBROUTINE archiveGetLine

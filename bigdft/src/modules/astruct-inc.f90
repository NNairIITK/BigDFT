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
  use dictionaries, only: f_err_raise, f_err_throw, max_field_length
  use dynamic_memory
  use yaml_strings, only: yaml_toa
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
  character(len=120) :: extra
  character(len=150) :: line
  logical :: lpsdbl, eof
  integer :: iat,ityp,ntyp,i,ierrsfx,nspol,nchrg
  ! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
  ! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  character(len=20), dimension(100) :: atomnames
  character(len = max_field_length) :: errmess

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
     call find_extra_info(line,extra,8)
     !print *,'then',iat,extra
     call parse_extra_info(astruct%attributes(iat),extra,errmess)
     if (len_trim(errmess) > 0) then
        call f_err_throw('At atom ' // trim(yaml_toa(iat)) // ': ' // trim(errmess),&
             & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
     else
        call astruct_at_from_dict(astruct%attributes(iat)%impl, &
             & ifrztyp = astruct%ifrztyp(iat), igspin = nspol, igchrg = nchrg)
        !now assign the array, following the rule
        astruct%input_polarization(iat)=1000*nchrg+sign(1, nchrg)*100+nspol
     end if

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
  character(len=120) :: extra
  character(len=150) :: line
  logical :: lpsdbl, reduced, eof, forces
  integer :: iat,ntyp,ityp,i,i_stat,nlines,istart,istop,count,nspol,nchrg
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3,alat4,alat5,alat6
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0,alat4d0,alat5d0,alat6d0
  character(len=20), dimension(100) :: atomnames
  character(max_field_length) :: errmess
  ! Store the file.
  character(len = 150), dimension(5000) :: lines

  energy = UNINITIALIZED(energy)
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
        call find_extra_info(line,extra,8)
        call parse_extra_info(astruct%attributes(iat),extra,errmess)
        if (len_trim(errmess) > 0) then
           call f_err_throw('At atom ' // trim(yaml_toa(iat)) // ': ' // trim(errmess),&
                & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
        else
           call astruct_at_from_dict(astruct%attributes(iat)%impl, &
                & ifrztyp = astruct%ifrztyp(iat), igspin = nspol, igchrg = nchrg)
           !now assign the array, following the rule
           astruct%input_polarization(iat)=1000*nchrg+sign(1, nchrg)*100+nspol
        end if

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

!> Read atomic positions from int file and create astruct structure from it
subroutine read_int_positions(iproc,ifile,astruct,comment,energy,fxyz,getLine)
  use module_defs, only: gp,UNINITIALIZED,Bohr_Ang, Radian_Degree,BIGDFT_INPUT_VARIABLES_ERROR
  use dictionaries, only: f_err_raise, max_field_length, f_err_throw
  use dynamic_memory
  use yaml_strings, only: yaml_toa
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
  character(len=*), parameter :: subname='read_atomic_positions'
  character(len=20) :: symbol
  character(len=20) :: tatonam
  character(len=120) :: extra
  character(len=150) :: line
  logical :: lpsdbl, eof
  integer :: iat,ityp,ntyp,i,ierrsfx,nchrg, nspol
  ! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
  ! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  integer :: na, nb, nc
  character(len=20), dimension(100) :: atomnames
  character(len = max_field_length) :: errmess

  call getLine(line, ifile, eof)
  if (eof) then
     write(*,*) "Error: unexpected end of file."
     stop
  end if
  energy = UNINITIALIZED(energy)
  read(line,*, iostat = ierrsfx) iat,astruct%units,astruct%angle,energy,comment
  if (ierrsfx /= 0) then
     read(line,*, iostat = ierrsfx) iat,astruct%units,astruct%angle,energy
     write(comment, "(A)") ""
     if (ierrsfx /= 0) then
        read(line,*, iostat = ierrsfx) iat,astruct%units,astruct%angle
        if (ierrsfx /= 0) then
           read(line,*, iostat = ierrsfx) iat
           write(astruct%units, "(A)") "bohr"
           write(astruct%angle, "(A)") "radian"
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
     !int input file, allow extra information
     call getLine(line, ifile, eof)
     if (f_err_raise(eof,"Unexpected end of file.",err_name='BIGDFT_RUNTIME_ERROR')) return


     !!if (lpsdbl) then
     !!   read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
     !!else
     !!   read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
     !!end if
     call check_line_integrity()
     !print *,'extra',iat,extra
     call find_extra_info(line,extra,14)
     !print *,'then',iat,extra
     call parse_extra_info(astruct%attributes(iat),extra,errmess)
     if (len_trim(errmess) > 0) then
        call f_err_throw('At atom ' // trim(yaml_toa(iat)) // ': ' // trim(errmess),&
             & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
     else
        call astruct_at_from_dict(astruct%attributes(iat)%impl, &
             & ifrztyp = astruct%ifrztyp(iat), igspin = nspol, igchrg = nchrg)
        !now assign the array, following the rule
        astruct%input_polarization(iat)=1000*nchrg+sign(1, nchrg)*100+nspol
     end if

     tatonam=trim(symbol)
!!!     end if
     if (lpsdbl) then
        astruct%rxyz_int(1,iat)=rxd0
        astruct%rxyz_int(2,iat)=ryd0
        astruct%rxyz_int(3,iat)=rzd0
        astruct%ixyz_int(1,iat)=na
        astruct%ixyz_int(2,iat)=nb
        astruct%ixyz_int(3,iat)=nc
     else
        astruct%rxyz_int(1,iat)=real(rx,gp)
        astruct%rxyz_int(2,iat)=real(ry,gp)
        astruct%rxyz_int(3,iat)=real(rz,gp)
        astruct%ixyz_int(1,iat)=na
        astruct%ixyz_int(2,iat)=nb
        astruct%ixyz_int(3,iat)=nc
     end if
     if (astruct%units=='angstroem' .or. astruct%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        astruct%rxyz_int(1,iat)=astruct%rxyz_int(1,iat) / Bohr_Ang
     end if
     if (astruct%angle=='degree') then
        ! if degree to radian
        astruct%rxyz_int(2:3,iat)=astruct%rxyz_int(2:3,iat) / Radian_Degree
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
     fxyz = f_malloc_ptr((/ 3, iat /),id='fxyz')
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
  call astruct_set_n_types(astruct, ntyp)

  astruct%atomnames(1:astruct%ntypes)=atomnames(1:astruct%ntypes)


contains

  !> stop the code and warns if the status of the line is not good
  subroutine check_line_integrity()
    use yaml_output, only: yaml_toa
    use dictionaries, only: f_err_raise
    implicit none


    if (lpsdbl) then
       read(line,*,iostat=ierrsfx)symbol,na,rxd0,nb,ryd0,nc,rzd0
    else
       read(line,*,iostat=ierrsfx)symbol,na,rx,nb,ry,nc,rz
    end if

    if (f_err_raise(ierrsfx/=0,'The line'//trim(yaml_toa(iat+2))//&
         ' of the atomic position is not valid, check if it is in DOS format!',&
         err_name='BIGDFT_LINALG_ERROR')) return

  end subroutine check_line_integrity

END SUBROUTINE read_int_positions


!> Local routines with explicit interface (used as arguments)
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


!> Used as arguments
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


!> put the atomic positions inside the box
!! accoding to geocode value and cell if present
subroutine rxyz_inside_box(astruct,rxyz)
  use module_defs, only: gp
  use dynamic_memory, only: f_memcpy
  implicit none
  !> description of the atomic structure
  type(atomic_structure), intent(inout) :: astruct
  !> position to work on alternatively.
  !! if present, the positions given by astruct are only considered as input
  real(gp), dimension(3,astruct%nat), intent(out), target, optional :: rxyz
  !local variables
  integer :: iat
  real(gp), dimension(:,:), pointer :: rxyz_

  rxyz_ => astruct%rxyz
  if (present(rxyz)) rxyz_ => rxyz
  
  !atoms inside the box.
  select case(astruct%geocode)
  case('P')
     do iat=1,astruct%nat
        rxyz_(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
        rxyz_(2,iat)=modulo(astruct%rxyz(2,iat),astruct%cell_dim(2))
        rxyz_(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
     end do
  case('S')
     if (present(rxyz)) then
        do iat=1,astruct%nat
           rxyz_(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
           rxyz_(2,iat)=       astruct%rxyz(2,iat)
           rxyz_(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        end do
     else
        do iat=1,astruct%nat
           rxyz_(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
           rxyz_(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        end do
     end if
  case('W')
     if (present(rxyz)) then
        do iat=1,astruct%nat
           rxyz_(1,iat)=       astruct%rxyz(1,iat)
           rxyz_(2,iat)=       astruct%rxyz(2,iat)
           rxyz_(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        end do
     else
        do iat=1,astruct%nat
           rxyz_(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
        end do
     end if
  case('F')
     if (present(rxyz)) call f_memcpy(src=astruct%rxyz,dest=rxyz_)
     !Do nothing!
  end select

end subroutine rxyz_inside_box

!> Find extra information
subroutine find_extra_info(line,extra,nspace)
  implicit none
  character(len=*), intent(in) :: line
  character(len=120), intent(out) :: extra
  integer,intent(in) :: nspace
  !local variables
  logical :: space
  integer :: i,ispace

  i=1
  space=.true.
  ispace=-1
  !print *,'line',line
  find_space : do
     !toggle the space value for each time
     if ((line(i:i) == ' ' .or. line(i:i) == char(9)) .neqv. space) then
        ispace=ispace+1
        space=.not. space
     end if
     !print *,line(i:i),ispace
     if (ispace==nspace) then
        extra=line(i:min(len(line),i+len(extra)-1))
        exit find_space
     end if
     if (i==len(line)) then
        !print *,'AAA',extra
        extra='nothing'
        exit find_space
     end if
     i=i+1
  end do find_space
END SUBROUTINE find_extra_info


!> Parse extra information
subroutine parse_extra_info(att, extra, errmess)
  use yaml_parse
  use dictionaries
  implicit none
  !Arguments
  type(f_tree), intent(out) :: att
  character(len=*), intent(in) :: extra
  character(len=max_field_length), intent(out) :: errmess
  !Local variables
  character(len=4) :: suffix
  logical :: go
  integer :: ierr,ierr1,ierr2,nspol,nchrg
  type(dictionary), pointer :: dict
  !case with all the information
  !print *,iat,'ex'//trim(extra)//'ex'

  write(errmess, "(A)") " "
  nullify(att%impl)
  if (index(extra, ":") > 0) then
     ! YAML case.
     call yaml_parse_from_string(dict, extra)
     if (dict_len(dict) > 0) att%impl => dict .pop. 0
     call dict_free(dict)
  else
     ! Old case.
     read(extra,*,iostat=ierr) nspol,nchrg,suffix
     if (extra == 'nothing') then !case with empty information
        nspol=0
        nchrg=0
        suffix='    '
     else if (ierr /= 0) then !case with partial information
        read(extra,*,iostat=ierr1) nspol,suffix
        if (ierr1 == 0) then
           !Format nspol frzchain
           nchrg=0
           call valid_frzchain(trim(suffix),go)
           if (.not. go) then
              suffix='    '
              read(suffix,*,iostat=ierr2) nchrg
              if (ierr2 /= 0) then
                 nchrg = 0
                 call error
              end if
           end if
        else
           !Format frzchain
           call valid_frzchain(trim(extra),go)
           if (go) then
              suffix=trim(extra)
              nspol=0
              nchrg=0
           else
              read(extra,*,iostat=ierr2) nspol
              if (ierr2 /=0) call error
              suffix='    '
              nchrg=0
              nspol=0
           end if
        end if
     end if

     ! convert everything into a dict.
     call dict_init(att%impl)
     if (nspol /= 0) call set(att%impl // ASTRUCT_ATT_IGSPIN, nspol)
     if (nchrg /= 0) call set(att%impl // ASTRUCT_ATT_IGCHRG, nchrg)
     if (len_trim(suffix) > 0) call set(att%impl // ASTRUCT_ATT_FROZEN, suffix)
     
     if (dict_size(att%impl) == 0) then
        call dict_free(att%impl)
        nullify(att%impl)
     end if
  end if

contains

  subroutine error
    implicit none
    write(errmess, '(a)') 'wrong additional data, read was "' // trim(extra) // '".'
  END SUBROUTINE error

END SUBROUTINE parse_extra_info


!> Check the validity of the chain for frozen atom option
subroutine valid_frzchain(frzchain,go)
  implicit none
  character(len=*), intent(in) :: frzchain
  logical, intent(out) :: go

  ! x: fix x direction
  ! y: fix z direction
  ! z: fix z direction
  ! b: fix bond length (internal coordinates)
  ! p: fix angle phi (internal coordinates)
  ! t: fix angle theta (internal coordinates)

  go = frzchain == 'f'    .or. &
       frzchain == 'fx'   .or. &
       frzchain == 'fy'   .or. &
       frzchain == 'fz'   .or. &
       frzchain == 'fxy'  .or. &
       frzchain == 'fxz'  .or. &
       frzchain == 'fyz'  .or. &
       frzchain == 'fxyz' .or. &
       frzchain == 'f'    .or. &
       frzchain == 'fb'   .or. &
       frzchain == 'fp'   .or. &
       frzchain == 'ft'   .or. &
       frzchain == 'fbp'  .or. &
       frzchain == 'fbt'  .or. &
       frzchain == 'fyt'  .or. &
       frzchain == 'fbpt'
  if (.not.go .and. len_trim(frzchain) >= 3) then
     go = (frzchain(1:1) == 'f' .and. verify(frzchain(2:), '0123456789') == 0) .or. &
          (frzchain(1:2) == 'fb' .and. verify(frzchain(3:), '12') == 0)
  end if

END SUBROUTINE valid_frzchain


!> Define the frozen type for the given atom
!! f: all atoms are frozen
!! fx: x direction frozen
!! fy: y direction frozen
!! fz: z direction frozen
!! fxz, fxy, fyz,fxyz
!! Move the atom also in a plane given by f111 (Miller indices)
!! Frozen also atoms per block given by fb#if_of_the_block (modified by FL)
!! This function is related to move_this_coordinate
subroutine frozen_ftoi(frzchain,ifrztyp,ierr)
  implicit none
  character(len=4), intent(in) :: frzchain !< Chain to be read
  integer, intent(out) :: ifrztyp          !< Integer coding the frozen type
  integer, intent(out) :: ierr             !< Error code

  ierr = 0
  select case(frzchain)
  case('')
     ifrztyp = 0
  case('f','fxyz')
     ifrztyp = 111
  case('fx')
     ifrztyp = 100
  case('fy')
     ifrztyp = 010
  case('fz')
     ifrztyp = 001
  case('fxz')
     ifrztyp = 101
  case('fxy')
     ifrztyp = 110
  case('fyz')
     ifrztyp = 011
  case('fbpt')
     ifrztyp = 222
  case('fb')
     ifrztyp = 200
  case('fp')
     ifrztyp = 020
  case('ft')
     ifrztyp = 002
  case('fbt')
     ifrztyp = 202
  case('fbp')
     ifrztyp = 220
  case('fpt')
     ifrztyp = 022
  case default
     !Check if we freeze the displacement of the atom only in a plane given by the Miller indices
     if (frzchain(1:1) == 'f' .and. verify(frzchain(2:), '0123456789') == 0) then
        read(frzchain(2:4), *) ifrztyp
        ! Check if 1 <= ifrztyp <= 999
        if (ifrztyp < 1 .or. ifrztyp > 999) ierr = 2
        ! f001 will give 9001 value.
        ifrztyp = 9000 + ifrztyp
     else if (frzchain(1:2) == 'fb' .and. verify(frzchain(3:), '12 ') == 0) then !space nedded since frzchain is a 4 character string
        ! (FL) atom possibly frozen in moving blocks
        read(frzchain(3:), *) ifrztyp
        ! Two blocks are possible
        if (ifrztyp < 1 .or. ifrztyp > 2) ierr = 2
        ! fb1 will give 1001 value.
        ifrztyp = 1000 + ifrztyp
     else
        !The frozen type is not correct!
        ierr = 1
     end if
  end select

END SUBROUTINE frozen_ftoi


!> Convert ifrztyp into the chain format
subroutine frozen_itof(ifrztyp,frzchain)
  use yaml_output, only: yaml_toa
  implicit none
  integer, intent(in) :: ifrztyp
  character(len=4), intent(out) :: frzchain

  select case(ifrztyp)
  case(0)
     frzchain = '    '
  case(111)
     frzchain = 'fxyz'
  case(100)
     frzchain = 'fx  '
  case(010)
     frzchain = 'fy  '
  case(001)
     frzchain = 'fz  '
  case(101)
     frzchain = 'fxz '
  case(110)
     frzchain = 'fxy '
  case(011)
     frzchain = 'fyz '
  case(222)
     frzchain = 'fbpt'
  case(200)
     frzchain = 'fb  '
  case(020)
     frzchain = 'fp  '
  case(002)
     frzchain = 'ft  '
  case(202)
     frzchain = 'fbt '
  case(220)
     frzchain = 'fbp '
  case(022)
     frzchain = 'fpt '
  case(1001)
     frzchain = 'fb1 '
  case(1002)
     frzchain = 'fb2 '
  case(9000:9999)
     frzchain ='f'//adjustl(yaml_toa(ifrztyp))
  case default
     print *,'Bug in frozen_itof'
     stop
  end select

END SUBROUTINE frozen_itof

!> The function which controls all the moving positions
!! This function is related to frozen_ftoi
pure function move_this_coordinate(ifrztyp,ixyz)
  implicit none
  integer, intent(in) :: ifrztyp !< Type of frozen atom
  integer, intent(in) :: ixyz    !w coordinates (x=1, y=2; z=3)
  logical :: move_this_coordinate

  move_this_coordinate = &
       ifrztyp == 0 .or.                     & !Not frozen at all!
       (ifrztyp == 100 .and. ixyz /= 1) .or. & !fx
       (ifrztyp == 010 .and. ixyz /= 2) .or. & !fy
       (ifrztyp == 001 .and. ixyz /= 3) .or. & !fz
       (ifrztyp == 110 .and. ixyz == 3) .or. & !fxy
       (ifrztyp == 101 .and. ixyz == 2) .or. & !fxz
       (ifrztyp == 011 .and. ixyz == 1)        !fyz
  !print *,"MOVE",ifrztyp,ixyz,move_this_coordinate

END FUNCTION move_this_coordinate


!>Write the extra info necessary for the output file
subroutine write_extra_info(extra,natpol,ifrztyp)
  use ao_inguess, only: charge_and_spol
  implicit none 
  integer, intent(in) :: natpol,ifrztyp
  character(len=120), intent(out) :: extra
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


!> Write xyz atomic file.
subroutine wtxyz(iunit,energy,rxyz,astruct,comment)
  use module_defs, only: Bohr_Ang,UNINITIALIZED
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(in) :: energy
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz

  !local variables
  character(len=20) :: symbol
  character(len=10) :: name
  character(len=11) :: units
  character(len=120) :: extra
  integer :: iat,j
  real(gp) :: xmax,ymax,zmax,factor


  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,astruct%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(astruct%units) == 'angstroem' .or. trim(astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if

  if (energy /= 0.0_gp .and. energy /= UNINITIALIZED(energy)) then
     write(iunit,'(i6,2x,a,2x,1pe24.17,1x,a,2x,a)') astruct%nat,trim(units),energy,'(Ha)',trim(comment)
  else
     write(iunit,'(i6,2x,a,2x,a)') astruct%nat,trim(units),trim(comment)
  end if

  select case(astruct%geocode)
  case('P')
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('S')
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('W')
     write(iunit,'(a,3(1x,1pe24.17))')'wire',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('F')
     write(iunit,*)'free'
  end select

  do iat=1,astruct%nat
     name=trim(astruct%atomnames(astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:min(len(name),5))
     end if

     call write_extra_info(extra,astruct%input_polarization(iat),astruct%ifrztyp(iat))

     write(iunit,'(a5,1x,3(1x,1pe24.17),2x,a)')symbol,(rxyz(j,iat)*factor,j=1,3),trim(extra)
  enddo

END SUBROUTINE wtxyz


!> Add the forces in the position file for the xyz system
subroutine wtxyz_forces(iunit,fxyz,astruct)
  implicit none
  integer, intent(in) :: iunit
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=20) :: symbol
  character(len=10) :: name

  ! Please don't change the keyword here.
  ! It is for this stupid XYZ file format, and this keyword is
  ! needed for force recognition in V_Sim for instance.
  write(iunit,*)'forces'

  do iat=1,astruct%nat
     name=trim(astruct%atomnames(astruct%iatype(iat)))
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


!> Write ascii file (atomic position). 
subroutine wtascii(iunit,energy,rxyz,astruct,comment)
  use module_defs, only: Bohr_Ang,UNINITIALIZED
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(in) :: energy
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
  !local variables
  character(len=2) :: symbol
  character(len=120) :: extra
  character(len=10) :: name
  integer :: iat,j
  real(gp) :: xmax,ymax,zmax,factor(3)

  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,astruct%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(astruct%units) == 'angstroem' .or. trim(astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
  else
     factor=1.0_gp
  end if

  write(iunit, "(A,A)") "# BigDFT file - ", trim(comment)
  write(iunit, "(3e24.17)") astruct%cell_dim(1)*factor(1), 0.d0, astruct%cell_dim(2)*factor(2)
  write(iunit, "(3e24.17)") 0.d0,                                0.d0, astruct%cell_dim(3)*factor(3)

  write(iunit, "(A,A)") "#keyword: ", trim(astruct%units)
  if (trim(astruct%units) == "reduced") write(iunit, "(A,A)") "#keyword: bohr"
  select case(astruct%geocode)
  case('P')
     write(iunit, "(A)") "#keyword: periodic"
  case('S')
     write(iunit, "(A)") "#keyword: surface"
  case('W')
     write(iunit, "(A)") "#keyword: wire"
  case('F')
     write(iunit, "(A)") "#keyword: freeBC"
  end select

  if (energy /= 0.d0 .and. energy /= UNINITIALIZED(energy)) then
     write(iunit, "(A,e24.17,A)") "#metaData: totalEnergy= ", energy, " Ht"
  end if

  if (trim(astruct%units) == "reduced") then
     select case(astruct%geocode)
     case('P')
        factor(1) = 1._gp / astruct%cell_dim(1)
        factor(2) = 1._gp / astruct%cell_dim(2)
        factor(3) = 1._gp / astruct%cell_dim(3)
     case('S')
        factor(1) = 1._gp / astruct%cell_dim(1)
        factor(3) = 1._gp / astruct%cell_dim(3)
     case('W')
        factor(3) = 1._gp / astruct%cell_dim(3)
     end select
  end if

  do iat=1,astruct%nat
     name=trim(astruct%atomnames(astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:2)
     end if

     call write_extra_info(extra,astruct%input_polarization(iat),astruct%ifrztyp(iat))     

     write(iunit,'(3(1x,1pe24.17),2x,a2,2x,a)') (rxyz(j,iat)*factor(j),j=1,3),symbol,trim(extra)
  end do

END SUBROUTINE wtascii


subroutine wtascii_forces(iunit,fxyz,astruct)
  implicit none
  integer, intent(in) :: iunit
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: fxyz
  !local variables
  integer :: iat,j
  character(len=1) :: endline

  if (astruct%nat ==0) return
  !write the first position
  iat=1
  if (astruct%nat==iat) then
     endline=']'
  else
     endline=char(92)
  end if

  write(iunit, "(A,3(1pe25.17,A),a)") "#metaData: forces (Ha/Bohr) =[",(fxyz(j,iat), ";",j=1,3),' '//endline
  !then the rest until the second-last
  do iat=2,astruct%nat
     if (astruct%nat==iat) then
        endline=']'
     else
        endline=char(92)
     end if
     write(iunit, "(A,3(1pe25.17,A),a)") "# ",(fxyz(j,iat), ";",j=1,3),' '//endline
  end do
end subroutine wtascii_forces

!> Write int atomic file.
subroutine wtint(iunit,energy,rxyz,astruct,comment,na,nb,nc)
  use module_defs, only: Bohr_Ang,UNINITIALIZED,Radian_Degree
  use module_base, only: f_err_throw
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(in) :: energy
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
  integer,dimension(astruct%nat),intent(in) :: na, nb, nc

  !local variables
  character(len=20) :: symbol
  character(len=10) :: name
  character(len=11) :: units, angle
  character(len=120) :: extra
  integer :: iat
  real(gp) :: xmax,ymax,zmax,factor,factor_angle


  xmax=0.0_gp
  ymax=0.0_gp
  zmax=0.0_gp

  do iat=1,astruct%nat
     xmax=max(rxyz(1,iat),xmax)
     ymax=max(rxyz(2,iat),ymax)
     zmax=max(rxyz(3,iat),zmax)
  enddo
  if (trim(astruct%units) == 'angstroem' .or. trim(astruct%units) == 'angstroemd0') then
     factor=Bohr_Ang
     units='angstroemd0'
  else
     factor=1.0_gp
     units='atomicd0'
  end if
  if (trim(astruct%angle) == 'degree') then
     factor_angle=Radian_Degree
     angle='degree'
  else
     factor_angle=1.0_gp
     angle='radian'
  end if
  write(*,*) '(trim(astruct%angle)), angle',(trim(astruct%angle)), angle

  if (energy /= 0.0_gp .and. energy /= UNINITIALIZED(energy)) then
     write(iunit,'(i6,2x,a,2x,a,2x,1pe24.17,2x,a)') astruct%nat,trim(units),&
          trim(angle),energy,trim(comment)
  else
     write(iunit,'(i6,2x,a,2x,a,2x,a)') astruct%nat,trim(units),trim(angle),trim(comment)
  end if

  select case(astruct%geocode)
  case('P')
     call f_err_throw("Internal coordinates not implemented for periodic BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('S')
     call f_err_throw("Internal coordinates not implemented for surface BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('W')
     call f_err_throw("Internal coordinates not implemented for wire BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'wire',&
          astruct%cell_dim(1)*factor,astruct%cell_dim(2)*factor,astruct%cell_dim(3)*factor
  case('F')
     write(iunit,*)'free'
  end select

  do iat=1,astruct%nat
     name=trim(astruct%atomnames(astruct%iatype(iat)))
     if (name(3:3)=='_') then
        symbol=name(1:2)
     else if (name(2:2)=='_') then
        symbol=name(1:1)
     else
        symbol=name(1:min(len(name),5))
     end if

     call write_extra_info(extra,astruct%input_polarization(iat),astruct%ifrztyp(iat))

     write(iunit,'(a5,1x,3(1x,i6,2x,1pe24.17),2x,a)')symbol,na(iat),rxyz(1,iat)*factor,nb(iat),rxyz(2,iat)*factor_angle,&
          nc(iat),rxyz(3,iat)*factor_angle,trim(extra)
  enddo

END SUBROUTINE wtint

!> Check the position of atoms, verify no atoms have the same coordinates
subroutine check_atoms_positions(astruct, simplify)
  use module_defs, only: gp
  use yaml_output
  implicit none
  !Arguments
  logical, intent(in) :: simplify
  type(atomic_structure), intent(in) :: astruct
  !local variables
  integer, parameter :: iunit=9
  logical :: dowrite
  integer :: iat,nateq,jat,j

  nateq=0
  do iat=1,astruct%nat
     do jat=iat+1,astruct%nat
        if ((astruct%rxyz(1,iat)-astruct%rxyz(1,jat))**2+&
             (astruct%rxyz(2,iat)-astruct%rxyz(2,jat))**2+&
             (astruct%rxyz(3,iat)-astruct%rxyz(3,jat))**2 < 1.e-10_gp) then
           nateq=nateq+1
           call yaml_warning('ERROR: atoms' // trim(yaml_toa(iat)) // ' (' // &
                & trim(astruct%atomnames(astruct%iatype(iat))) // ') and ' // &
                & trim(yaml_toa(jat)) // ' (' // trim(astruct%atomnames(astruct%iatype(jat))) // &
                & ') have the same positions')
        end if
     end do
  end do

  if (nateq /= 0) then
     if (simplify) then
        call yaml_warning('Control your posinp file, cannot proceed')
        write(*,'(1x,a)',advance='no') 'Writing tentative alternative positions in the file posinp_alt...'
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
     end if
     stop 'check_atoms_positions'
  end if
END SUBROUTINE check_atoms_positions

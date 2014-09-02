!> @file
!! Routines for handling the structure atoms_data 
!! @author
!!    Copyright (C) 2011-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!!$subroutine astruct_nullify(astruct)
!!$  use module_types
!!$  implicit none
!!$  type(atomic_structure), intent(out) :: astruct
!!$
!!$  astruct%geocode = "F"
!!$  astruct%units = "bohr"
!!$  astruct%inputfile_format = "none"
!!$  astruct%nat = -1
!!$  astruct%ntypes = -1
!!$  astruct%sym%symObj = -1
!!$  nullify(astruct%sym%irrzon)
!!$  nullify(astruct%sym%phnons)
!!$end subroutine astruct_nullify
!!$
!!$!> Nullify a new atoms_data type.
!!$subroutine atoms_nullify(atoms)
!!$  use module_types
!!$  implicit none
!!$  type(atoms_data), intent(out) :: atoms
!!$
!!$  call astruct_nullify(atoms%astruct)
!!$
!!$  ! Arrays related to ntypes.
!!$  nullify(atoms%psppar)
!!$  nullify(atoms%nelpsp)
!!$  nullify(atoms%npspcode)
!!$  nullify(atoms%nzatom)
!!$  nullify(atoms%ixcpsp)
!!$  nullify(atoms%radii_cf)
!!$  ! parameters for NLCC
!!$  nullify(atoms%nlccpar)
!!$  nullify(atoms%nlcc_ngv)
!!$  nullify(atoms%nlcc_ngc)
!!$  ! Parameters for Linear input guess
!!$  nullify(atoms%rloc)
!!$
!!$  ! Arrays related to nat.
!!$  nullify(atoms%iasctype)
!!$  nullify(atoms%aocc)
!!$  nullify(atoms%amu)
!!$
!!$  nullify(atoms%paw_l)
!!$  nullify(atoms%paw_NofL)
!!$  nullify(atoms%paw_nofchannels)
!!$  nullify(atoms%paw_nofgaussians)
!!$  nullify(atoms%paw_Greal)
!!$  nullify(atoms%paw_Gimag)
!!$  nullify(atoms%paw_Gcoeffs)
!!$  nullify(atoms%paw_H_matrices)
!!$  nullify(atoms%paw_S_matrices)
!!$  nullify(atoms%paw_Sm1_matrices)
!!$end subroutine atoms_nullify

!> Add a displacement of atomic positions and put in the box
subroutine astruct_set_displacement(astruct, randdis)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(in) :: randdis !< random displacement

  integer :: iat
  real(gp) :: tt
  
  !Shake atoms if required.
  if (randdis > 0.d0) then
     do iat=1,astruct%nat
        if (astruct%ifrztyp(iat) == 0) then
           call random_number(tt)
           astruct%rxyz(1,iat)=astruct%rxyz(1,iat)+randdis*tt
           call random_number(tt)
           astruct%rxyz(2,iat)=astruct%rxyz(2,iat)+randdis*tt
           call random_number(tt)
           astruct%rxyz(3,iat)=astruct%rxyz(3,iat)+randdis*tt
        end if
     enddo
  end if

  !atoms inside the box.
  select case(astruct%geocode)
  case('P')
     do iat=1,astruct%nat
        astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
        astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),astruct%cell_dim(2))
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
     end do
  case('S')
     do iat=1,astruct%nat
        astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
     end do
  case('W')
     do iat=1,astruct%nat
        astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
     end do
  case('F')
     !Do nothing!
  end select
END SUBROUTINE astruct_set_displacement


!> For bindings only, use input_dicts module for Fortran usage.
subroutine astruct_merge_to_dict_binding(dict, astruct)
  use module_input_dicts, only: wrapper => astruct_merge_to_dict
  use module_types, only: atomic_structure
  use dictionaries, only: dictionary
  implicit none
  type(dictionary), pointer :: dict
  type(atomic_structure), intent(in) :: astruct

  call wrapper(dict, astruct,astruct%rxyz)
END SUBROUTINE astruct_merge_to_dict_binding


!> Find extra information
subroutine find_extra_info(line,extra,nspace)
  implicit none
  character(len=150), intent(in) :: line
  character(len=50), intent(out) :: extra
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
           read(suffix,*,iostat=ierr2) nchrg
           if (ierr2 /= 0) then
              call error
           else
              suffix='    '
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
           if (ierr2 /=0) then
              call error
           end if
           suffix='    '
           nchrg=0
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
  call frozen_ftoi(suffix,astruct%ifrztyp(iat),ierr)
  if (ierr /= 0) call error

!!!  if (trim(suffix) == 'f') then
!!!     !the atom is considered as blocked
!!!     astruct%ifrztyp(iat)=1
!!!  end if

contains

  subroutine error
    implicit none
    print *,extra
    write(*,'(1x,a,i0,a)')&
         'ERROR in input file for atom number ',iat,&
         ': after 4th column you can put the input polarisation(s) or the frozen chain (f,fxz,fy or f111, fb1, ...)'
    stop
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
function move_this_coordinate(ifrztyp,ixyz)
  use module_base
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


!> Write xyz atomic file.
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

  select case(atoms%astruct%geocode)
  case('P')
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('S')
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('W')
     write(iunit,'(a,3(1x,1pe24.17))')'wire',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('F')
     write(iunit,*)'free'
  end select

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


!> Write ascii file (atomic position). 
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
  write(iunit, "(3e24.17)") 0.d0,                                0.d0, atoms%astruct%cell_dim(3)*factor(3)

  write(iunit, "(A,A)") "#keyword: ", trim(atoms%astruct%units)
  if (trim(atoms%astruct%units) == "reduced") write(iunit, "(A,A)") "#keyword: bohr"
  select case(atoms%astruct%geocode)
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

  if (trim(atoms%astruct%units) == "reduced") then
     select case(atoms%astruct%geocode)
     case('P')
        factor(1) = 1._gp / atoms%astruct%cell_dim(1)
        factor(2) = 1._gp / atoms%astruct%cell_dim(2)
        factor(3) = 1._gp / atoms%astruct%cell_dim(3)
     case('S')
        factor(1) = 1._gp / atoms%astruct%cell_dim(1)
        factor(3) = 1._gp / atoms%astruct%cell_dim(3)
     case('W')
        factor(3) = 1._gp / atoms%astruct%cell_dim(3)
     end select
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



!> Write int atomic file.
subroutine wtint(iunit,energy,rxyz,atoms,comment,na,nb,nc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  integer,dimension(atoms%astruct%nat),intent(in) :: na, nb, nc

  !local variables
  character(len=20) :: symbol
  character(len=10) :: name
  character(len=11) :: units, angle
  character(len=50) :: extra
  integer :: iat
  real(gp) :: xmax,ymax,zmax,factor,factor_angle


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
  if (trim(atoms%astruct%angle) == 'degree') then
     factor_angle=Radian_Degree
     angle='degree'
  else
     factor_angle=1.0_gp
     angle='radian'
  end if
  write(*,*) '(trim(atoms%astruct%angle)), angle',(trim(atoms%astruct%angle)), angle

  if (energy /= 0. .and. energy /= UNINITIALIZED(energy)) then
     write(iunit,'(i6,2x,a,2x,a,2x,1pe24.17,2x,a)') atoms%astruct%nat,trim(units),&
          trim(angle),energy,trim(comment)
  else
     write(iunit,'(i6,2x,a,2x,a,2x,a)') atoms%astruct%nat,trim(units),trim(angle),trim(comment)
  end if

  select case(atoms%astruct%geocode)
  case('P')
     call f_err_throw("Internal coordinates not implemented for periodic BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'periodic',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('S')
     call f_err_throw("Internal coordinates not implemented for surface BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'surface',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('W')
     call f_err_throw("Internal coordinates not implemented for wire BC", err_name='BIGDFT_RUNTIME_ERROR')
     write(iunit,'(a,3(1x,1pe24.17))')'wire',&
          atoms%astruct%cell_dim(1)*factor,atoms%astruct%cell_dim(2)*factor,atoms%astruct%cell_dim(3)*factor
  case('F')
     write(iunit,*)'free'
  end select

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

     write(iunit,'(a5,1x,3(1x,i6,2x,1pe24.17),2x,a50)')symbol,na(iat),rxyz(1,iat)*factor,nb(iat),rxyz(2,iat)*factor_angle,&
                                                     nc(iat),rxyz(3,iat)*factor_angle,extra
  enddo

END SUBROUTINE wtint


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



!> Write atomic file in yaml format
subroutine wtyaml(iunit,energy,rxyz,atoms,wrtforces,forces, &
     & wrtlog, shift, hgrids)
  use module_base, only: f_err_throw
  use module_defs, only: Bohr_Ang, gp, UNINITIALIZED
  use module_types, only: atoms_data
  use yaml_output
  implicit none
  !Arguments
  logical, intent(in) :: wrtforces !< True if write the atomic forces
  logical, intent(in) :: wrtlog
  integer, intent(in) :: iunit
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz,forces
  real(gp), dimension(3), intent(in) :: shift, hgrids
  !local variables
  logical :: reduced, perx, pery, perz
  character(len=4) :: frzchain
  real(gp), dimension(3) :: xred
  real(gp) :: factor
  integer :: iat,ichg,ispol
  
  reduced=.false.
  Units: select case(trim(atoms%astruct%units))
  case('angstroem','angstroemd0')
     call yaml_map('Units','angstroem', unit = iunit)
     factor=Bohr_Ang
  case('reduced')
     if (.not. wrtlog) then
        call yaml_map('Units','reduced', unit = iunit)
        reduced=.true.
     end if
     factor = 1.0_gp
  case('atomic','atomicd0','bohr','bohrd0')
     ! Default
     factor=1.0_gp
     !call yaml_map('Units','bohr')
  case default
     call f_err_throw('Writing the atomic file. Error, unknown units ("'// trim(atoms%astruct%units)//'")', & 
          & err_name='BIGDFT_RUNTIME_ERROR')
  end select Units

  !cell information
  perx = .false.
  pery = .false.
  perz = .false.
  BC :select case(atoms%astruct%geocode)
  case('S')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(1)*factor), unit = iunit) !x
       call yaml_sequence('.inf', unit = iunit)             !y
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .false.
     perz = .true.
  case('W')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
       call yaml_sequence('.inf', unit = iunit)             !x
       call yaml_sequence('.inf', unit = iunit)             !y
       call yaml_sequence(yaml_toa(atoms%astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     perx = .false.
     pery = .false.
     perz = .true.
  case('P')
     call yaml_map('Cell',(/atoms%astruct%cell_dim(1)*factor, &
          & atoms%astruct%cell_dim(2)*factor, atoms%astruct%cell_dim(3)*factor/), unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .true.
     perz = .true.
  case('F')
     ! Default
     !call yaml_map('BC','free')
  end select BC

  !Write atomic positions
  call yaml_sequence_open('Positions', unit = iunit)
  do iat=1,atoms%astruct%nat
     call yaml_sequence(advance='no', unit = iunit)
     if (extra_info(iat)) then
        call yaml_mapping_open(flow=.true., unit = iunit)
     end if
     xred=rxyz(:,iat)
     if (reduced) then
        if (perx) xred(1)=rxyz(1,iat)/atoms%astruct%cell_dim(1)
        if (pery) xred(2)=rxyz(2,iat)/atoms%astruct%cell_dim(2)
        if (perz) xred(3)=rxyz(3,iat)/atoms%astruct%cell_dim(3)
     else
        !Multiply by the factor to have the right units
        xred = xred*factor
     end if
     if (wrtlog) then
        call print_one_atom(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
             xred,hgrids,iat)
!!$        call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
!!$             & xred,fmt="(g18.10)", unit = iunit, advance = "no")
!!$        xred(1:3) = rxyz(1:3,iat) / hgrids
!!$        write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, iat
!!$        call yaml_comment(gu, unit = iunit)
     else
        call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
             & xred,fmt="(1pg25.17)", unit = iunit)
     end if
     if (extra_info(iat)) then
        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
        if (ispol /=0) call yaml_map('IGSpin',ispol, unit = iunit)
        if (ichg /=0) call yaml_map('IGChg',ichg, unit = iunit)
        if (atoms%astruct%ifrztyp(iat) /= 0) then
           call frozen_itof(atoms%astruct%ifrztyp(iat),frzchain)
           call yaml_map('Frozen',frzchain, unit = iunit)
        end if
        call yaml_mapping_close(unit = iunit)
     end if
  end do
  call yaml_sequence_close(unit = iunit) !positions

  !Write atomic forces
  if (wrtforces) then
     call yaml_sequence_open('Forces (Ha/Bohr)', unit = iunit)
     do iat=1,atoms%astruct%nat
        call yaml_sequence(advance='no', unit = iunit)
        call yaml_map(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),forces(:,iat),fmt='(1pg25.17)', unit = iunit)
     end do
     call yaml_sequence_close(unit = iunit) !values
  end if
  if (wrtlog) then
     call yaml_map('Rigid Shift Applied (AU)',(/-shift(1),-shift(2),-shift(3)/),fmt='(1pg12.5)')
  else
     call yaml_mapping_open('Properties', unit = iunit)
     call yaml_map('Timestamp',yaml_date_and_time_toa(), unit = iunit)
     if (energy /= 0. .and. energy /= UNINITIALIZED(energy)) then
        call yaml_map("Energy (Ha)", energy, unit = iunit)
     end if
     call yaml_mapping_close(unit = iunit) !properties
  end if

contains

  function extra_info(iat)
    implicit none
    integer, intent(in) :: iat
    logical extra_info
    extra_info=atoms%astruct%input_polarization(iat) /=100 .or. atoms%astruct%ifrztyp(iat)/=0
  end function extra_info

  subroutine print_one_atom(atomname,rxyz,hgrids,id)
    implicit none
    integer, intent(in) :: id
    character(len=*), intent(in) :: atomname
    double precision, dimension(3), intent(in) :: rxyz,hgrids
    !local variables
    character(len=*), parameter :: fmtat='(1pg18.10)',fmtg='(F6.2)',fmti='(i4.4)'
    integer :: i

    call yaml_sequence_open(atomname,flow=.true.)
    do i=1,3
       call yaml_sequence(yaml_toa(rxyz(i),fmt=fmtat))
    end do
    call yaml_sequence_close(advance='no')
    call yaml_comment(trim(yaml_toa(rxyz/hgrids/factor,fmt=fmtg))//trim(yaml_toa(id,fmt=fmti))) !we can also put tabbing=

  end subroutine print_one_atom

END SUBROUTINE wtyaml


!> Calculate the charge and the spin polarisation to be placed on a given atom
!! RULE: natpol = c*1000 + sgn(c)*100 + s: charged and polarised atom (charge c, polarisation s)
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


subroutine atoms_write(atoms, filename, forces, energy, comment)
  use module_defs, only: gp
  use module_types
  use module_interfaces, only: write_atomic_file
  implicit none
  character(len = *), intent(in) :: comment
  character(len = *), intent(in) :: filename
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(:,:), pointer :: forces

  if (associated(forces)) then
     call write_atomic_file(filename,energy,atoms%astruct%rxyz,atoms%astruct%ixyz_int,atoms,comment,forces=forces)
  else
     call write_atomic_file(filename,energy,atoms%astruct%rxyz,atoms%astruct%ixyz_int,atoms,comment)
  end if
END SUBROUTINE atoms_write


!> Deallocate a new atoms_data type, for bindings.
subroutine atoms_empty(atoms)
  use module_atoms, only: atoms_data, deallocate_atoms_data
  implicit none
  type(atoms_data), intent(inout) :: atoms

  call deallocate_atoms_data(atoms)
END SUBROUTINE atoms_empty


subroutine atoms_set_name(atoms, ityp, name)
  use module_types
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: ityp
  character(len=1), dimension(20), intent(in) :: name

  write(atoms%astruct%atomnames(ityp), "(20A1)") name
END SUBROUTINE atoms_set_name


subroutine astruct_set_geometry(astruct, alat, geocode, format, units)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(in) :: alat(3)
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character, intent(in) :: format(5)
  character, intent(in) :: units(20)

  astruct%cell_dim(:) = alat(:)
  astruct%geocode = geocode
  write(astruct%inputfile_format, "(5A1)") format
  write(astruct%units, "(20A1)") units
END SUBROUTINE astruct_set_geometry

!> Accessors for bindings.
subroutine atoms_get(atoms, astruct, symObj)
  use module_types
  implicit none
  type(atoms_data), intent(in), target :: atoms
  type(atomic_structure), pointer :: astruct
  type(symmetry_data), pointer :: symObj

  astruct => atoms%astruct
  symObj => atoms%astruct%sym
END SUBROUTINE atoms_get


subroutine astruct_copy_nat(astruct, nat)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  integer, intent(out) :: nat

  nat = astruct%nat
END SUBROUTINE astruct_copy_nat


subroutine astruct_copy_ntypes(astruct, ntypes)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  integer, intent(out) :: ntypes

  ntypes = astruct%ntypes
END SUBROUTINE astruct_copy_ntypes


subroutine atoms_get_iatype(atoms, iatype)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, dimension(:), pointer :: iatype

  iatype => atoms%astruct%iatype
END SUBROUTINE atoms_get_iatype


!!$subroutine atoms_get_iasctype(atoms, iasctype)
!!$  use module_types
!!$  implicit none
!!$  type(atoms_data), intent(in) :: atoms
!!$  integer, dimension(:), pointer :: iasctype
!!$  !local variables
!!$  integer
!!$
!!$  iasctype => atoms%iasctype
!!$END SUBROUTINE atoms_get_iasctype


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


subroutine atoms_get_rxyz(atoms, rxyz)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz

  rxyz => atoms%astruct%rxyz
END SUBROUTINE atoms_get_rxyz


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
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:), pointer :: amu

  amu => atoms%amu
END SUBROUTINE atoms_get_amu


!!$subroutine atoms_get_aocc(atoms, aocc)
!!$  use module_types
!!$  implicit none
!!$  type(atoms_data), intent(in) :: atoms
!!$  real(gp), dimension(:,:), pointer :: aocc
!!$
!!$  aocc => atoms%aocc
!!$END SUBROUTINE atoms_get_aocc


!> get radii_cf values
subroutine atoms_get_radii_cf(atoms, radii_cf)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: radii_cf

  radii_cf => atoms%radii_cf
END SUBROUTINE atoms_get_radii_cf


subroutine atoms_get_psppar(atoms, psppar)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:,:), pointer :: psppar

  psppar => atoms%psppar
END SUBROUTINE atoms_get_psppar

subroutine atoms_get_nlccpar(atoms, nlccpar)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(:,:), pointer :: nlccpar

  nlccpar => atoms%nlccpar
END SUBROUTINE atoms_get_nlccpar

!subroutine atoms_get_ig_nlccpar(atoms, ig_nlccpar)
!  use module_types
!  implicit none
!  type(atoms_data), intent(in) :: atoms
!  real(gp), dimension(:,:), pointer :: ig_nlccpar
!
!  ig_nlccpar => atoms%ig_nlccpar
!END SUBROUTINE atoms_get_ig_nlccpar


subroutine astruct_copy_geometry_data(astruct, geocode, format, units)
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  character(len = 1), intent(out) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len = 5), intent(out) :: format
  character(len = 20), intent(out) :: units

  write(geocode, "(A1)") astruct%geocode
  write(format,  "(A5)") astruct%inputfile_format
  write(units,  "(A20)") astruct%units
END SUBROUTINE astruct_copy_geometry_data


subroutine atoms_copy_psp_data(atoms, natsc, donlcc)
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: natsc
  logical, intent(out) :: donlcc

  natsc = atoms%natsc
  donlcc = atoms%donlcc
END SUBROUTINE atoms_copy_psp_data


subroutine astruct_copy_name(astruct, ityp, name, ln)
  use module_types
  implicit none
  !Arguments
  type(atomic_structure), intent(in) :: astruct
  integer, intent(in) :: ityp
  character(len=1), dimension(20), intent(out) :: name
!  character(len=*), intent(out) :: name
  integer, intent(out) :: ln
  !Local variables 
  integer :: i,lname

  if (astruct%ntypes > 0) then
     lname = len(name)
     ln=min(len(trim(astruct%atomnames(ityp))),20)
     !print *,'lnt2',lnt
     do i = 1, ln, 1
     !name(i:i) = astruct%atomnames(ityp)(i:i)
     write(name(i),'(a1)') astruct%atomnames(ityp)(i:i)
     end do
     do i = ln + 1, lname, 1
        name(i) = ' '
     end do
  end if
END SUBROUTINE astruct_copy_name


subroutine astruct_copy_alat(astruct, alat)
  use module_defs, only: gp
  use module_types
  implicit none
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(out) :: alat(3)

  alat(1) = astruct%cell_dim(1)
  alat(2) = astruct%cell_dim(2)
  alat(3) = astruct%cell_dim(3)
END SUBROUTINE astruct_copy_alat


!!$!> Module used for the input positions lines variables
!!$module position_files
!!$   implicit none
!!$   contains
!!$   subroutine directGetLine(line, ifile, eof)
!!$      !Arguments
!!$      integer, intent(in) :: ifile
!!$      character(len=150), intent(out) :: line
!!$      logical, intent(out) :: eof
!!$      !Local variables
!!$      integer :: i_stat
!!$
!!$      eof = .false.
!!$      read(ifile,'(a150)', iostat = i_stat) line
!!$      if (i_stat /= 0) eof = .true.
!!$   END SUBROUTINE directGetLine
!!$
!!$   subroutine archiveGetLine(line, ifile, eof)
!!$      !Arguments
!!$      integer, intent(in) :: ifile
!!$      character(len=150), intent(out) :: line
!!$      logical, intent(out) :: eof
!!$      !Local variables
!!$      integer :: i_stat
!!$      !The argument ifile is not used but it is used as argument routine
!!$      !eof = .false.
!!$      eof = (ifile /= ifile)
!!$      call extractNextLine(line, i_stat)
!!$      if (i_stat /= 0) eof = .true.
!!$   END SUBROUTINE archiveGetLine
!!$end module position_files


!> Write an atomic file
!! Yaml output included
subroutine write_atomic_file(filename,energy,rxyz,ixyz,atoms,comment,forces)
  use module_base
  use module_types
  use module_input_dicts
  use yaml_output
  use internal_coordinates, only : xyzint
  implicit none
  character(len=*), intent(in) :: filename,comment
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: energy
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  integer,dimension(3,atoms%astruct%nat),intent(in) :: ixyz
  real(gp), dimension(3,atoms%astruct%nat), intent(in), optional :: forces
  !local variables
  character(len = 15) :: arFile
  integer :: iunit
  character(len = 1024) :: fname
  real(gp), dimension(3), parameter :: dummy = (/ 0._gp, 0._gp, 0._gp /)
  type(dictionary), pointer :: dict
  real(gp),dimension(:,:),allocatable :: rxyz_int
  real(kind=8),parameter :: degree=1.d0


  if (trim(filename) == "stdout") then
     iunit = 6
  else
     iunit = 9
     write(fname,"(A)") trim(filename)//'.'//trim(atoms%astruct%inputfile_format)
     if (atoms%astruct%inputfile_format == 'yaml') then
        call yaml_set_stream(unit = iunit, filename = trim(fname), &
             & record_length = 92, setdefault = .false., tabbing = 0)
     else
        open(unit = iunit, file = trim(fname))
     end if
  end if

  select case(atoms%astruct%inputfile_format)
     case('xyz')
        call wtxyz(iunit,energy,rxyz,atoms,comment)
        if (present(forces)) call wtxyz_forces(9,forces,atoms)
     case('ascii')
        call wtascii(iunit,energy,rxyz,atoms,comment)
        if (present(forces)) call wtascii_forces(9,forces,atoms)
     case ('int')
       !if (.not.present(na) .or. .not.present(nb) .or. .not.present(nc)) then
       !    call f_err_throw('na, nb, nc must be present to write a file in internal coordinates', &
       !         err_name='BIGDFT_RUNTIME_ERROR')
       !end if
        rxyz_int = f_malloc((/3,atoms%astruct%nat/),id='rxyz_int')
        !call wtint(iunit,energy,rxyz,atoms,comment,ixyz(1,:),ixyz(2,:),ixyz(3,:))
        call xyzint(atoms%astruct%rxyz, atoms%astruct%nat, &
             atoms%astruct%ixyz_int(1,:), atoms%astruct%ixyz_int(2,:),atoms%astruct%ixyz_int(3,:), &
             degree, rxyz_int)
        call wtint(iunit,energy,rxyz_int,atoms,comment,ixyz(1,:),ixyz(2,:),ixyz(3,:))
        call f_free(rxyz_int)
     case('yaml')
        call yaml_new_document(unit = iunit)
        call dict_init(dict)
        call astruct_merge_to_dict(dict, atoms%astruct, rxyz, comment)
        call yaml_dict_dump(dict, unit = iunit)
        call dict_free(dict)
!!$     if (len_trim(comment) > 0) call yaml_comment(comment, unit = iunit)
!!$     if (present(forces)) then
!!$        call wtyaml(iunit,energy,rxyz,atoms,.true.,forces, .false., dummy, dummy)
!!$     else
!!$        call wtyaml(iunit,energy,rxyz,atoms,.false.,rxyz, .false., dummy, dummy)
!!$     end if
     case default
        call f_err_throw('Writing the atomic file. Error, unknown file format ("'//&
             trim(atoms%astruct%inputfile_format)//'")', &
             err_name='BIGDFT_RUNTIME_ERROR')
   end select

  if (iunit /= 6) then
     if (atoms%astruct%inputfile_format == 'yaml') then
        call yaml_close_stream(unit = iunit)
     else
        close(unit = iunit)
     end if
     ! Add to archive
     if (index(filename, "posout_") == 1 .or. index(filename, "posmd_") == 1) then
        write(arFile, "(A)") "posout.tar.bz2"
        if (index(filename, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
        call addToCompress(trim(arFile), len(trim(arFile)), trim(fname), len(trim(fname)))
     end if
  end if
END SUBROUTINE write_atomic_file


!> Calculate the coefficient for moving atoms following the ifrztyp
subroutine frozen_alpha(ifrztyp,ixyz,alpha,alphai)
  use module_base
  implicit none
  integer, intent(in) :: ifrztyp !< Frozen code of the iat atom 
  integer, intent(in) :: ixyz    !< Direction (1=x,y=2,z=3)
  real(gp), intent(in) :: alpha
  real(gp), intent(out) :: alphai
  !Local variables
  logical :: move_this_coordinate

  if (move_this_coordinate(ifrztyp,ixyz)) then
     alphai=alpha
  else
     alphai=0.0_gp
  end if
 
END SUBROUTINE frozen_alpha

!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: rxyz=txyz+alpha*sxyz
!!   all the shift are inserted into the box if there are periodic directions
!!   if the atom are frozen they are not moved
subroutine atomic_axpy(atoms,txyz,alpha,sxyz,rxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz

  do iat=1,atoms%astruct%nat
     !adjust the moving of the atoms following the frozen direction
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,alpha,alphaz)

     if (atoms%astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=modulo(txyz(2,iat)+alphay*sxyz(2,iat),atoms%astruct%cell_dim(2))
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(txyz(1,iat)+alphax*sxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=modulo(txyz(3,iat)+alphaz*sxyz(3,iat),atoms%astruct%cell_dim(3))
     else
        rxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
        rxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
        rxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
     end if
  end do

END SUBROUTINE atomic_axpy


!>Routine for moving atomic positions, takes into account the 
!!   frozen atoms and the size of the cell
!!   synopsis: fxyz=txyz+alpha*sxyz
!!   update the forces taking into account the frozen atoms
!!   do not apply the modulo operation on forces 
subroutine atomic_axpy_forces(atoms,txyz,alpha,sxyz,fxyz)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: alpha
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: txyz,sxyz
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: fxyz
  !local variables
  integer :: iat
  real(gp) :: alphax,alphay,alphaz
  
  do iat=1,atoms%astruct%nat
     !adjust the moving of the forces following the frozen direction
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,alpha,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,alpha,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,alpha,alphaz)

     fxyz(1,iat)=txyz(1,iat)+alphax*sxyz(1,iat)
     fxyz(2,iat)=txyz(2,iat)+alphay*sxyz(2,iat)
     fxyz(3,iat)=txyz(3,iat)+alphaz*sxyz(3,iat)
  end do
  
END SUBROUTINE atomic_axpy_forces


!>Calculate the scalar product between atomic positions by considering
!!   only non-blocked atoms
subroutine atomic_dot(atoms,x,y,scpr)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: x,y
  real(gp), intent(out) :: scpr
  !local variables
  integer :: iat
  real(gp) :: scpr1,scpr2,scpr3
  real(gp) :: alphax,alphay,alphaz

  scpr=0.0_gp

  do iat=1,atoms%astruct%nat
     call frozen_alpha(atoms%astruct%ifrztyp(iat),1,1.0_gp,alphax)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),2,1.0_gp,alphay)
     call frozen_alpha(atoms%astruct%ifrztyp(iat),3,1.0_gp,alphaz)
     scpr1=alphax*x(1,iat)*y(1,iat)
     scpr2=alphay*x(2,iat)*y(2,iat)
     scpr3=alphaz*x(3,iat)*y(3,iat)
     scpr=scpr+scpr1+scpr2+scpr3
  end do
  
END SUBROUTINE atomic_dot


!>z=alpha*A*x + beta* y
subroutine atomic_gemv(atoms,m,alpha,A,x,beta,y,z)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: m
  real(gp), intent(in) :: alpha,beta
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: x
  real(gp), dimension(m), intent(in) :: y
  real(gp), dimension(m,3,atoms%astruct%nat), intent(in) :: A
  real(gp), dimension(m), intent(out) :: z
  !local variables
  integer :: iat,i,j
  real(gp) :: mv,alphai
  
  do i=1,m
     mv=0.0_gp
     do iat=1,atoms%astruct%nat
        do j=1,3
           call frozen_alpha(atoms%astruct%ifrztyp(iat),j,A(i,j,iat),alphai)
           mv=mv+alphai*x(j,iat)
        end do
     end do
     z(i)=alpha*mv+beta*y(i)
  end do

END SUBROUTINE atomic_gemv


!> rxyz=txyz+alpha*sxyz
subroutine atomic_coordinate_axpy(atoms,ixyz,iat,t,alphas,r)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ixyz,iat
  real(gp), intent(in) :: t,alphas
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: r
  !local variables
  logical :: periodize
  real(gp) :: alat,alphai

  if (ixyz == 1) then
     alat=atoms%astruct%cell_dim(1)
  else if (ixyz == 2) then
     alat=atoms%astruct%cell_dim(2)
  else if (ixyz == 3) then
     alat=atoms%astruct%cell_dim(3)
  else
     alat = -1
     write(0,*) "Internal error"
     stop
  end if
  
  periodize= atoms%astruct%geocode == 'P' .or. &
       (atoms%astruct%geocode == 'S' .and. ixyz /= 2)

  call frozen_alpha(atoms%astruct%ifrztyp(iat),ixyz,alphas,alphai)

  if (periodize) then
     r=modulo(t+alphai,alat)
  else
     r=t+alphai
  end if

END SUBROUTINE atomic_coordinate_axpy


!> This routine does the same operations as read_atomic_file but uses inputs from memory
!! as input positions instead of inputs from file
!! Useful for QM/MM implementation of BigDFT-ART
!! @author Written by Laurent K Beland 2011 UdeM
subroutine initialize_atomic_file(iproc,atoms,rxyz)
  use module_base
  use module_types
  use module_interfaces, except_this_one => initialize_atomic_file
  use m_ab6_symmetry
  use yaml_output
  implicit none
  integer, intent(in) :: iproc
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(:,:), pointer :: rxyz
  !local variables
  character(len=*), parameter :: subname='initialize_atomic_file'
  integer :: iat,i,ierr

  atoms%amu = f_malloc_ptr(atoms%astruct%nat+ndebug)

  if (atoms%astruct%geocode=='S') then 
        atoms%astruct%cell_dim(2)=0.0_gp
  else if (atoms%astruct%geocode=='F') then !otherwise free bc    
        atoms%astruct%cell_dim(1)=0.0_gp
        atoms%astruct%cell_dim(2)=0.0_gp
        atoms%astruct%cell_dim(3)=0.0_gp
  else
        atoms%astruct%cell_dim(1)=0.0_gp
        atoms%astruct%cell_dim(2)=0.0_gp
        atoms%astruct%cell_dim(3)=0.0_gp
  end if

  !reduced coordinates are possible only with periodic units
  if (atoms%astruct%units == 'reduced' .and. atoms%astruct%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

   !convert the values of the cell sizes in bohr
  if (atoms%astruct%units=='angstroem' .or. atoms%astruct%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     atoms%astruct%cell_dim(1)=atoms%astruct%cell_dim(1)/Bohr_Ang
     atoms%astruct%cell_dim(2)=atoms%astruct%cell_dim(2)/Bohr_Ang
     atoms%astruct%cell_dim(3)=atoms%astruct%cell_dim(3)/Bohr_Ang
  else if (atoms%astruct%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     atoms%astruct%cell_dim(1)=real(atoms%astruct%cell_dim(1),gp)
     atoms%astruct%cell_dim(2)=real(atoms%astruct%cell_dim(2),gp)
     atoms%astruct%cell_dim(3)=real(atoms%astruct%cell_dim(3),gp)
  else
     call yaml_warning('Length units in input file unrecognized')
     call yaml_warning('recognized units are angstroem or atomic = bohr')
     call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
  endif
  
  do iat=1,atoms%astruct%nat
     !xyz input file, allow extra information
     
     if (atoms%astruct%units == 'reduced') then !add treatment for reduced coordinates
        rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
        if (atoms%astruct%geocode == 'P') rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
        rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
     else if (atoms%astruct%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(2,iat)=modulo(rxyz(2,iat),atoms%astruct%cell_dim(2))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     else if (atoms%astruct%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),atoms%astruct%cell_dim(1))
        rxyz(3,iat)=modulo(rxyz(3,iat),atoms%astruct%cell_dim(3))
     end if
 
     if (atoms%astruct%units=='angstroem' .or. atoms%astruct%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/Bohr_Ang
        enddo
     else if (atoms%astruct%units == 'reduced') then 
        rxyz(1,iat)=rxyz(1,iat)*atoms%astruct%cell_dim(1)
        if (atoms%astruct%geocode == 'P') rxyz(2,iat)=rxyz(2,iat)*atoms%astruct%cell_dim(2)
        rxyz(3,iat)=rxyz(3,iat)*atoms%astruct%cell_dim(3)
     endif
  enddo

  !control atom positions
  call check_atoms_positions(atoms%astruct,(iproc == 0))

  ! We delay the calculation of the symmetries.
  atoms%astruct%sym%symObj = -1
  nullify(atoms%astruct%sym%irrzon)
  nullify(atoms%astruct%sym%phnons)

END SUBROUTINE initialize_atomic_file


!> Check the position of atoms, verify no atoms have the same coordinates
subroutine check_atoms_positions(astruct, simplify)
  use module_defs, only: gp
  use module_atoms, only: atomic_structure
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

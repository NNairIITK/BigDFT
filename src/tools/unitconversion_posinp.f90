!> @file
!!  Routines to handle posinp files
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Rotates the molecular structure in the input file posinp.xyz and 
!!  writes the result in the file rot_posinp.xyz
!!  Comparing the BigDFT energies of the original and rotated configuration 
!!  can help to estimate the accuracy the the chosen parameter set( hgrid, crmult etc).
PROGRAM unitconversion_posinp

   implicit none

   integer, parameter :: natx=2000
   character(len=5) :: atomname(natx)
   character(len=12) :: unitsold, unitsnew
   character(len=20) :: extra(natx) 
   character(len=100) :: line,line2
   real(kind=8), dimension(3,natx) :: pos
   real(kind=8) :: scale
   integer :: iat,nat,ierror

   write(*,*) 'reading atomic positions from file posinp'
   open(unit=9,file='posinp.xyz',status='unknown')
   read(9,*) nat,unitsold
   if (nat.gt.natx) stop 'increase natx'
   read(9,*) line2
   do iat=1,nat
      read(9,'(a100)')line
      read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat) ,extra(iat)
      if (ierror .ne. 0) then
         read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat) 
         extra(iat)='  '
      end if
   enddo
   close(unit=9)

   write(*,*) 'old units:',unitsold
   if (trim(unitsold)=='angstroem') then
      unitsnew='bohr'
      scale=1.889725989d0
   else if (trim(unitsold)=='angstroemd0') then
      unitsnew='bohrd0'
      scale=1.889725989d0
   else if (trim(unitsold)=='bohr') then
      unitsnew='angstroem'
      scale=1.d0/1.889725989d0
   else if (trim(unitsold)=='bohrd0') then
      unitsnew='angstroemd0'
      scale=1.d0/1.889725989d0
   else if (trim(unitsold)=='atomic') then
      unitsnew='angstroem'
      scale=1.d0/1.889725989d0
   else if (trim(unitsold)=='atomicd0') then
      unitsnew='angstroemd0'
      scale=1.d0/1.889725989d0
   else
      stop 'units not recognized'
   endif
   write(*,*) 'new units:',unitsnew

   write(*,*) 'writing atomic positions to file convert_posinp'
   open(unit=9,file='convert_posinp.xyz',status='unknown')
   write(9,*) nat, unitsnew
   write(9,'(a100)') line2
   do iat=1,nat
      write(9,'(a5,3x,3(1x,e24.17),4x,a)') atomname(iat),pos(1,iat)*scale,pos(2,iat)*scale,pos(3,iat)*scale,extra(iat)
   end do
   close(unit=9)

END PROGRAM unitconversion_posinp

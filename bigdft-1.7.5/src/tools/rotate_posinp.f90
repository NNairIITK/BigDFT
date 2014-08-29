!> @file
!!  Program to rotate atomic coordinates 
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
PROGRAM rotate_posinp

   implicit none
   integer, parameter :: natx=2000
   character(len=5) :: atomname(natx)
   character(len=12) :: units
   character(len=8) :: boundary
   character(len=20) :: extra(natx) 
   character(len=120) :: line
!   character(len=120) :: line2
   real(kind=8), dimension(3,natx) :: pos
   real(kind=8), dimension(3) :: pos_s
   real(kind=8) :: scale,t1,t2,t3,phi_1,phi_2,phi_3,alat1,alat2,alat3
   real(kind=8), parameter :: PI=3.141592654d0
   integer :: nat,iat,ierror

   write(*,*) 'reading atomic positions from file posinp'
   open(unit=9,file='posinp.xyz',status='unknown')
   read(9,*) nat,units
   write(*,*) nat,units
   if (nat.gt.natx) stop 'increase natx'
!   read(9,*) line2
!   write(*,*) line2
   read(9,*) boundary,alat1,alat2,alat3
   write(*,*) boundary,alat1,alat2,alat3

   ! put center of mass at origin
   pos_s(1)=0.d0
   pos_s(2)=0.d0
   pos_s(3)=0.d0
   do iat=1,nat
      write(*,*)  '      ---------------------------------- '
      read(9,'(a100)')line
      write(*,*) iat, line

      read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat) ,extra(iat)
      if (ierror .ne. 0) then
         read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat) 
         extra(iat)='  '
      end if
      write(*,*) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat)

      pos_s(1)=pos_s(1)+pos(1,iat)
      pos_s(2)=pos_s(2)+pos(2,iat)
      pos_s(3)=pos_s(3)+pos(3,iat)
   end do
   close(unit=9)
   pos_s(1)=pos_s(1)/nat
   pos_s(2)=pos_s(2)/nat
   pos_s(3)=pos_s(3)/nat  
   do iat=1,nat
      pos(1,iat)=pos(1,iat)-pos_s(1)
      pos(2,iat)=pos(2,iat)-pos_s(2)        
      pos(3,iat)=pos(3,iat)-pos_s(3)
   end do

   write(*,*) 'rotations in degrees (0<= Phi <=360):'
   write(*,*)
   write(*,*) 'around z axis / in xy-plane:'
   read(*,*) phi_1
   phi_1=2.d0*PI*phi_1/360.d0
   write(*,*) 'around y axis / in xz-plane:'
   read(*,*) phi_2
   phi_2=2.d0*PI*phi_2/360.d0
   write(*,*) 'around x axis / in yz-plane:'
   read(*,*) phi_3
   phi_3=2.d0*PI*phi_3/360.d0

   do iat=1,nat
      t1=cos(phi_1)*pos(1,iat)+sin(phi_1)*pos(2,iat)
      t2=-sin(phi_1)*pos(1,iat)+cos(phi_1)*pos(2,iat)
      pos(1,iat)=t1
      pos(2,iat)=t2

      t1=cos(phi_2)*pos(1,iat)+sin(phi_2)*pos(3,iat)
      t3=-sin(phi_2)*pos(1,iat)+cos(phi_2)*pos(3,iat)
      pos(1,iat)=t1
      pos(3,iat)=t3

      t2=cos(phi_3)*pos(2,iat)+sin(phi_3)*pos(3,iat)
      t3=-sin(phi_3)*pos(2,iat)+cos(phi_3)*pos(3,iat)
      pos(2,iat)=t2
      pos(3,iat)=t3
   end do

   ! shift back center of mass to original position
   do iat=1,nat
      pos(1,iat)=pos(1,iat)+pos_s(1)
      pos(2,iat)=pos(2,iat)+pos_s(2)
      pos(3,iat)=pos(3,iat)+pos_s(3)
   end do

   write(*,*) 'scaling factor=?'
   read(*,*) scale

   write(*,*) 'writing atomic positions to file rot_posinp'
   open(unit=9,file='rotate_posinp.xyz',status='unknown')
   write(9,*) nat, units
!   write(9,'(a100)') line2
   write(9,'(a,3(2x,e21.14))') boundary,scale*alat1,scale*alat2,scale*alat3
   do iat=1,nat
      write(9,'(a5,3x,3(1x,e17.10),4x,a)') atomname(iat),   & 
      pos(1,iat)*scale,pos(2,iat)*scale,pos(3,iat)*scale,extra(iat)
   end do
   close(unit=9)

END PROGRAM rotate_posinp

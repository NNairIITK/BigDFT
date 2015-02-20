!> @file
!!  Program transforming the coordinates of a molecule
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Rotates the molecular structure such that its principal axis coincide with the x,y,z axis
PROGRAM p_axis_posinp

   implicit none
   integer, parameter :: natx=2000
   character(len=5) :: atomname(natx)
   character(len=12) :: units
   character(len=20) :: extra(natx) 
   character(len=100) :: line,line2
   real(kind=8), dimension(3,natx) :: pos
   real(kind=8) :: haratio
   integer :: nat,iat,ierror

   write(*,*) 'reading atomic positions from file posinp.xyz'
   open(unit=9,file='posinp.xyz',status='unknown')
   read(9,*) nat,units
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

   call ha_trans(nat,pos,haratio)
   write(*,*) 'ratio between largest/smallest EV ',haratio

   write(*,*) 'writing atomic positions to file axis_posinp.xyz'
   open(unit=9,file='axis_posinp.xyz',status='unknown')
   write(9,*) nat, units,  haratio
   write(9,'(a100)') line2
   do iat=1,nat
      write(9,'(a5,3x,3(1x,e24.17),4x,a)') atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat),extra(iat)
   end do
   close(unit=9)

END PROGRAM p_axis_posinp


subroutine ha_trans(nat,pos,haratio)
   implicit none
   !Arguments
   integer, intent(in) :: nat
   real(kind=8), dimension(3,nat), intent(out) :: pos
   real(kind=8), intent(out) :: haratio
   !Local variables
   integer, parameter :: lwork=100
   real(kind=8), dimension(lwork) :: work
   real(kind=8), dimension(3,1000) :: pos_n
   real(kind=8), dimension(3,3) :: theta
   real(kind=8), dimension(3) :: pos_s,theta_e
   integer :: iat,i,j,info

   if (nat.gt.1000) stop 'ha_trans'

   ! positions relative to center of mass
   pos_s(1)=0.d0
   pos_s(2)=0.d0
   pos_s(3)=0.d0
   do iat=1,nat
      pos_s(1)=pos_s(1)+pos(1,iat)
      pos_s(2)=pos_s(2)+pos(2,iat)
      pos_s(3)=pos_s(3)+pos(3,iat)
   enddo
   pos_s(1)=pos_s(1)/nat
   pos_s(2)=pos_s(2)/nat
   pos_s(3)=pos_s(3)/nat  

   do iat=1,nat
      pos(1,iat)=pos(1,iat)-pos_s(1)
      pos(2,iat)=pos(2,iat)-pos_s(2)        
      pos(3,iat)=pos(3,iat)-pos_s(3)
   enddo

   ! Calculate inertia tensor theta
   do 10,j=1,3
      do 10,i=1,3
         10 theta(i,j)=0.d0
         do iat=1,nat
            theta(1,1)=theta(1,1) + pos(2,iat)* pos(2,iat) + &  
            pos(3,iat)* pos(3,iat)
            theta(2,2)=theta(2,2) + pos(1,iat)* pos(1,iat) + &  
            pos(3,iat)* pos(3,iat)
            theta(3,3)=theta(3,3) + pos(1,iat)* pos(1,iat) + &   
            pos(2,iat)* pos(2,iat)

            theta(1,2)=theta(1,2) - pos(1,iat)*pos(2,iat)
            theta(1,3)=theta(1,3) - pos(1,iat)*pos(3,iat)
            theta(2,3)=theta(2,3) - pos(2,iat)*pos(3,iat)
            theta(2,1)=theta(1,2)
            theta(3,1)=theta(1,3)
            theta(3,2)=theta(2,3)
         enddo
         ! diagonalize theta
         call DSYEV('V','U',3,theta(1,1),3,theta_e(1),work(1),lwork,info)        
         haratio=theta_e(3)/theta_e(1)

         do iat=1,nat
            pos_n(1,iat) = theta(1,1)*pos(1,iat)+ &
               &   theta(2,1)*pos(2,iat)+ &
               &   theta(3,1)*pos(3,iat)
            pos_n(2,iat) = theta(1,2)*pos(1,iat)+ &
               &   theta(2,2)*pos(2,iat)+ &
               &   theta(3,2)*pos(3,iat)
            pos_n(3,iat) = theta(1,3)*pos(1,iat)+ &
               &   theta(2,3)*pos(2,iat)+ &
               &   theta(3,3)*pos(3,iat)
         enddo

         do iat=1,nat
            pos(1,iat)=pos_n(1,iat)
            pos(2,iat)=pos_n(2,iat)
            pos(3,iat)=pos_n(3,iat)
         enddo

      END SUBROUTINE ha_trans

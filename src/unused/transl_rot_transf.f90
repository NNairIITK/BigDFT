!> @file
!!  Unused program
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


        PROGRAM TR_ROT_TRANSF

        implicit real(kind=8) (a-h,o-z)
        character(len=20) atomname(1000)    
        dimension pos(3,1000),pos_n(3,1000),pos_s(3),alat(3)
        dimension theta(3,3)
        parameter(PI=3.141592654d0)       

        write(*,*) 'reading atomic positions from file totoin'
        open(unit=9,file='totoin',status='unknown')
        read(9,*) nat
!        write(*,*) nat
        read(9,*) alat(1),tt,alat(2)
!        write(*,*) alat(1),tt,alat(2)
        read(9,*) tt,tt,alat(3)
!        write(*,*) tt,tt,alat(3)
        pos_s(1:3)=0.d0
        do iat=1,nat
             read(9,*) pos(1,iat),pos(2,iat),pos(3,iat),atomname(iat)
!             write(*,*) pos(1,iat),pos(2,iat),pos(3,iat),atomname(iat)
             pos_s(1)=pos_s(1)+pos(1,iat)
             pos_s(2)=pos_s(2)+pos(2,iat)
             pos_s(3)=pos_s(3)+pos(3,iat)
        enddo
        close(9)
        pos_s(1)=pos_s(1)/nat
        pos_s(2)=pos_s(2)/nat
        pos_s(3)=pos_s(3)/nat  
        do iat=1,nat
             pos(1,iat)=pos(1,iat)-pos_s(1)
             pos(2,iat)=pos(2,iat)-pos_s(2)        
             pos(3,iat)=pos(3,iat)-pos_s(3)
        enddo
        write(*,*) 'Drehungen in Grad (0<= Phi <=360):'
        write(*,*)
        write(*,*) 'Um die x-Achse / in der yz-Ebene:'
        read(*,*) phi_1
        phi_1=2*PI*phi_1/360
        write(*,*) 'Um die y-Achse / in der xz-Ebene:'
        read(*,*) phi_2
        phi_2=2*PI*phi_2/360
        write(*,*) 'Um die z-Achse / in der xy-Ebene:'
        read(*,*) phi_3
        phi_3=2*PI*phi_3/360

        write(*,*) 'Translationen:'
        write(*,*)
        write(*,*) 'Entlang der x-Achse:'
        read(*,*) tr_x
        write(*,*) 'Entlang der y-Achse:'
        read(*,*) tr_y
        write(*,*) 'Entlang der z-Achse:'
        read(*,*) tr_z        


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

             
        enddo

        do iat=1,nat
             pos(1,iat)=pos(1,iat)+pos_s(1)+tr_x
             pos(2,iat)=pos(2,iat)+pos_s(2)+tr_y        
             pos(3,iat)=pos(3,iat)+pos_s(3)+tr_z
        enddo
      
        write(*,*) 'writeing atomic positions to file totout'
        open(unit=9,file='totout',status='unknown')
        write(9,*) nat
        write(9,'(3(e17.10))') alat(1),0.d0,alat(2)
        write(9,'(3(e17.10))') 0.d0,0.d0,alat(3)
        do iat=1,nat
        write(9,'(3(1x,e17.10),3x,a20)') pos(1,iat),pos(2,iat),pos(3,iat),atomname(iat)
        enddo
        close(9)
        end



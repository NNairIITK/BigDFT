!> @file
!!   
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine to_cf(x,xc,xf,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ibyz_c,ibyz_f)
    implicit none

    integer n1,n2,n3
    integer nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

    real(kind=8) x(0:n1,2,0:n2,2,0:n3,2)
    real(kind=8) xc(0:n1,0:n2,0:n3)
    real(kind=8) xf(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

    integer ibyz_c(2,0:n2,0:n3)
    integer ibyz_f(2,0:n2,0:n3)

    integer i1,i2,i3

!    xc=0.d0
!    xf=0.d0

    do i2=0,n2
        do i3=0,n3
            do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
                xc(i1,i2,i3)=x(i1,1,i2,1,i3,1)
            enddo
        enddo
    enddo

    do i2=nfl2,nfu2
        do i3=nfl3,nfu3
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                xf(1,i1,i2,i3)=x(i1,2,i2,1,i3,1)
                     
                xf(2,i1,i2,i3)=x(i1,1,i2,2,i3,1)
                     
                xf(3,i1,i2,i3)=x(i1,2,i2,2,i3,1)
                     
                xf(4,i1,i2,i3)=x(i1,1,i2,1,i3,2)
                     
                xf(5,i1,i2,i3)=x(i1,2,i2,1,i3,2)
                     
                xf(6,i1,i2,i3)=x(i1,1,i2,2,i3,2)
                     
                xf(7,i1,i2,i3)=x(i1,2,i2,2,i3,2)
            enddo
        enddo
    enddo

END SUBROUTINE to_cf

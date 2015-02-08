!> @file
!!  Common routines of convolutions
!! @author
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> A periodic synthesis (backward) wavelet transformation
!! the input array x is not overwritten
subroutine synthese_per_old(nd1,nd2,nd3,x,y,ww)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  syn_rot_per_old(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per_old(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per_old(nd3,nt,ww,y)

END SUBROUTINE synthese_per_old


!>   A periodic synthesis (backward) wavelet transformation
!!   the input array x is not overwritten
subroutine synthese_per_old_self(nd1,nd2,nd3,x,y,ww)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  syn_rot_per_old(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per_old(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per_old(nd3,nt,ww,x)

END SUBROUTINE synthese_per_old_self


!>   An analysis (forward) periodic wavelet transformation
!!   the input array y is not overwritten
subroutine analyse_per_old(nd1,nd2,nd3,y,x,ww)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: y
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  ana_rot_per_old(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per_old(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per_old(nd3,nt,ww,x)

END SUBROUTINE analyse_per_old


!>   An analysis (forward) periodic wavelet transformation
!!   the input array y is not overwritten
subroutine analyse_per_old_self(nd1,nd2,nd3,y,x,ww)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(in) :: y
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: ww
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,i1
  nt=(nd2+1)*(nd3+1)
  call  ana_rot_per_old(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per_old(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per_old(nd3,nt,ww,y)

END SUBROUTINE analyse_per_old_self


!> Synthesis repeated in periodic conditions
subroutine syn_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='syn_repeated_per'
  integer :: nn1,nn2,nn3,i_trans,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  if (num_trans >= 1)  then

     yy = f_malloc((nd1+1)*(nd2+1)*(nd3+1),id='yy')
     xx = f_malloc((nd1+1)*(nd2+1)*(nd3+1),id='xx')

  endif

  if (num_trans >= 2) then

     nn1=(nd1+1)/2-1
     nn2=(nd2+1)/2-1
     nn3=(nd3+1)/2-1

     ww = f_malloc((nn1+1)*(nn2+1)*(nn3+1),id='ww')

     do i_trans=1,num_trans-1

        n1=2*(n1+1)-1
        n2=2*(n2+1)-1
        n3=2*(n3+1)-1

        if (n1.gt.nd1) stop 'n1 beyond borders'
        if (n2.gt.nd2) stop 'n2 beyond borders'
        if (n3.gt.nd3) stop 'n3 beyond borders'

        !$omp parallel do default(private) shared (x,xx,n1,n2,n3)
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 i=1+i1+i2*(n1+1)+i3*(n1+1)*(n2+1) 
                 xx(i)=x(i1,i2,i3)
              enddo
           enddo
        enddo
        !$omp end parallel do

        call synthese_per_old(n1,n2,n3,xx,yy,ww)

        !$omp parallel do default(private) shared (x,yy,n1,n2,n3)
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 i=1+i1+i2*(n1+1)+i3*(n1+1)*(n2+1) 
                 x(i1,i2,i3)=yy(i)
              enddo
           enddo
        enddo
        !$omp end parallel do
     enddo

     call f_free(ww)

  endif

  if (num_trans >= 1) then

     n1=2*(n1+1)-1
     n2=2*(n2+1)-1
     n3=2*(n3+1)-1

     call synthese_per_old_self(n1,n2,n3,x,xx,yy)

     call f_free(xx)
     call f_free(yy)

  endif

END SUBROUTINE syn_repeated_per


!> Periodic Analysis repeated
subroutine ana_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='ana_repeated_per'
  integer :: i_trans,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  n1=nd1
  n2=nd2
  n3=nd3

  !write(21,*) 'ana_repeated_per'

  if (num_trans.ge.1)  then

     yy = f_malloc((nd1+1)*(nd2+1)*(nd3+1),id='yy')
     xx = f_malloc((nd1+1)*(nd2+1)*(nd3+1),id='xx')

     call analyse_per_old_self(n1,n2,n3,x,yy,xx)

     n1=(n1+1)/2-1
     n2=(n2+1)/2-1
     n3=(n3+1)/2-1

  endif

  if (num_trans.ge.2) then

     ww = f_malloc((n1+1)*(n2+1)*(n3+1),id='ww')

     do i_trans=2,num_trans

        !$omp parallel do default(private) shared (x,xx,n1,n2,n3)
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 i=1+i1+i2*(n1+1)+i3*(n1+1)*(n2+1)
                 xx(i)=x(i1,i2,i3)
                 i=i+1
              enddo
           enddo
        enddo
        !$omp end parallel do    
        call analyse_per_old(n1,n2,n3,xx,yy,ww)

        !$omp parallel do default(private) shared (x,yy,n1,n2,n3)
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 i=1+i1+i2*(n1+1)+i3*(n1+1)*(n2+1)
                 x(i1,i2,i3)=yy(i)
              enddo
           enddo
        enddo
        !$omp end parallel do

        n1=(n1+1)/2-1
        n2=(n2+1)/2-1
        n3=(n3+1)/2-1

     enddo

     call f_free(ww)

  endif

  if (num_trans.ge.1) then 
     call f_free(xx)
     call f_free(yy)
  endif

END SUBROUTINE ana_repeated_per

!> @file
!!  Optimzed convolution routines
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine ana_rot_per(n,ndat,x,y)
  use module_base
  implicit none

!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: ci,di
  real(wp), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp /
  data cg  / -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp  /

  integer mod_arr(-7:2*n+8)   
  real(wp) :: ci1,ci2,ci3,ci4,ci5,ci6,ci7,ci8
  real(wp) :: di1,di2,di3,di4,di5,di6,di7,di8

  !write(*,*) 'ana_rot_per executed'

  call fill_mod_arr(mod_arr,-7,2*n+8,2*n+2)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,cg,ch,ndat,n,mod_arr)
!$omp do 

  do j=0,ndat/8-1
     do i=0,n
        ci1=0.e0_wp
        ci2=0.e0_wp
        ci3=0.e0_wp
        ci4=0.e0_wp
        ci5=0.e0_wp
        ci6=0.e0_wp
        ci7=0.e0_wp
        ci8=0.e0_wp

        di1=0.e0_wp
        di2=0.e0_wp
        di3=0.e0_wp
        di4=0.e0_wp
        di5=0.e0_wp
        di6=0.e0_wp
        di7=0.e0_wp
        di8=0.e0_wp

        do l=-7,8
           k= mod_arr(l+2*i)

           ci1=ci1+ch(l)*x(k,j*8+1)
           ci2=ci2+ch(l)*x(k,j*8+2)
           ci3=ci3+ch(l)*x(k,j*8+3)
           ci4=ci4+ch(l)*x(k,j*8+4)
           ci5=ci5+ch(l)*x(k,j*8+5)
           ci6=ci6+ch(l)*x(k,j*8+6)
           ci7=ci7+ch(l)*x(k,j*8+7)
           ci8=ci8+ch(l)*x(k,j*8+8)

           di1=di1+cg(l)*x(k,j*8+1)
           di2=di2+cg(l)*x(k,j*8+2)
           di3=di3+cg(l)*x(k,j*8+3)
           di4=di4+cg(l)*x(k,j*8+4)
           di5=di5+cg(l)*x(k,j*8+5)
           di6=di6+cg(l)*x(k,j*8+6)
           di7=di7+cg(l)*x(k,j*8+7)
           di8=di8+cg(l)*x(k,j*8+8)
        end do
        y(j*8+1,    i)=ci1
        y(j*8+2,    i)=ci2
        y(j*8+3,    i)=ci3
        y(j*8+4,    i)=ci4
        y(j*8+5,    i)=ci5
        y(j*8+6,    i)=ci6
        y(j*8+7,    i)=ci7
        y(j*8+8,    i)=ci8

        y(j*8+1,n+1+i)=di1
        y(j*8+2,n+1+i)=di2
        y(j*8+3,n+1+i)=di3
        y(j*8+4,n+1+i)=di4
        y(j*8+5,n+1+i)=di5
        y(j*8+6,n+1+i)=di6
        y(j*8+7,n+1+i)=di7
        y(j*8+8,n+1+i)=di8
     end do
  end do

  !$omp end do
  
  !$omp do  

  do j=(ndat/8)*8+1,ndat
     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           k= mod_arr(l+2*i)
           ci=ci+ch(l)*x(k    ,j)
           di=di+cg(l)*x(k    ,j)
        end do
        y(j,i)=ci
        y(j,n+1+i)=di
     end do
  end do
  !$omp end do

!$omp end parallel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on ana rot per',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE ana_rot_per



subroutine syn_rot_per(n,ndat,x,y)
  use module_base
  implicit none
!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: so,se
  real(wp), dimension(-8:9) :: ch,cg
  !       Daubechy S16
  data ch  /  0.e0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.e0_wp /
  data cg  / 0.e0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.e0_wp /

  integer mod_arr(-4:n+4)
  real(wp) :: so1,so2,so3,so4,so5,so6,so7,so8
  real(wp) :: se1,se2,se3,se4,se5,se6,se7,se8

!   real(gp)::tel
!   integer::ncount1,ncount2,ncount_rate,ncount_max
!   integer::mflop
!   
!   ! (n+1): the range of i (average)
!   ! 8: filter length for l (average)
!   ! 4: number of flops in one line
!   ! 2: even and odd 
!   mflop=ndat*(n+1)*8*4*2
!  call system_clock(ncount1,ncount_rate,ncount_max)




  call fill_mod_arr(mod_arr,-4,n+4,n+1)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,cg,ch,ndat,n,mod_arr)
!$omp do 
  do j=0,ndat/8-1
     do i=0,n
      l=4 
        k=mod_arr(i-l)
        se1=ch(8)*x(k,j*8+1)+cg(8)*x(n+1+k,j*8+1)
        se2=ch(8)*x(k,j*8+2)+cg(8)*x(n+1+k,j*8+2)
        se3=ch(8)*x(k,j*8+3)+cg(8)*x(n+1+k,j*8+3)
        se4=ch(8)*x(k,j*8+4)+cg(8)*x(n+1+k,j*8+4)
        se5=ch(8)*x(k,j*8+5)+cg(8)*x(n+1+k,j*8+5)
        se6=ch(8)*x(k,j*8+6)+cg(8)*x(n+1+k,j*8+6)
        se7=ch(8)*x(k,j*8+7)+cg(8)*x(n+1+k,j*8+7)
        se8=ch(8)*x(k,j*8+8)+cg(8)*x(n+1+k,j*8+8)

      l=-4 
        k=mod_arr(i-l)
        so1=ch(-7)*x(k,j*8+1)+cg(-7)*x(n+1+k,j*8+1)
        so2=ch(-7)*x(k,j*8+2)+cg(-7)*x(n+1+k,j*8+2)
        so3=ch(-7)*x(k,j*8+3)+cg(-7)*x(n+1+k,j*8+3)
        so4=ch(-7)*x(k,j*8+4)+cg(-7)*x(n+1+k,j*8+4)
        so5=ch(-7)*x(k,j*8+5)+cg(-7)*x(n+1+k,j*8+5)
        so6=ch(-7)*x(k,j*8+6)+cg(-7)*x(n+1+k,j*8+6)
        so7=ch(-7)*x(k,j*8+7)+cg(-7)*x(n+1+k,j*8+7)
        so8=ch(-7)*x(k,j*8+8)+cg(-7)*x(n+1+k,j*8+8)

        do l=-3,3
           k=mod_arr(i-l)

           se1=se1+ch(2*l  )*x(k,j*8+1)+cg(2*l  )*x(n+1+k,j*8+1)
           se2=se2+ch(2*l  )*x(k,j*8+2)+cg(2*l  )*x(n+1+k,j*8+2)
           se3=se3+ch(2*l  )*x(k,j*8+3)+cg(2*l  )*x(n+1+k,j*8+3)
           se4=se4+ch(2*l  )*x(k,j*8+4)+cg(2*l  )*x(n+1+k,j*8+4)
           se5=se5+ch(2*l  )*x(k,j*8+5)+cg(2*l  )*x(n+1+k,j*8+5)
           se6=se6+ch(2*l  )*x(k,j*8+6)+cg(2*l  )*x(n+1+k,j*8+6)
           se7=se7+ch(2*l  )*x(k,j*8+7)+cg(2*l  )*x(n+1+k,j*8+7)
           se8=se8+ch(2*l  )*x(k,j*8+8)+cg(2*l  )*x(n+1+k,j*8+8)

           so1=so1+ch(2*l+1)*x(k,j*8+1)+cg(2*l+1)*x(n+1+k,j*8+1)
           so2=so2+ch(2*l+1)*x(k,j*8+2)+cg(2*l+1)*x(n+1+k,j*8+2)
           so3=so3+ch(2*l+1)*x(k,j*8+3)+cg(2*l+1)*x(n+1+k,j*8+3)
           so4=so4+ch(2*l+1)*x(k,j*8+4)+cg(2*l+1)*x(n+1+k,j*8+4)
           so5=so5+ch(2*l+1)*x(k,j*8+5)+cg(2*l+1)*x(n+1+k,j*8+5)
           so6=so6+ch(2*l+1)*x(k,j*8+6)+cg(2*l+1)*x(n+1+k,j*8+6)
           so7=so7+ch(2*l+1)*x(k,j*8+7)+cg(2*l+1)*x(n+1+k,j*8+7)
           so8=so8+ch(2*l+1)*x(k,j*8+8)+cg(2*l+1)*x(n+1+k,j*8+8)
        end do
        y(j*8+1,2*i  )=se1
        y(j*8+2,2*i  )=se2
        y(j*8+3,2*i  )=se3
        y(j*8+4,2*i  )=se4
        y(j*8+5,2*i  )=se5
        y(j*8+6,2*i  )=se6
        y(j*8+7,2*i  )=se7
        y(j*8+8,2*i  )=se8

        y(j*8+1,2*i+1)=so1
        y(j*8+2,2*i+1)=so2
        y(j*8+3,2*i+1)=so3
        y(j*8+4,2*i+1)=so4
        y(j*8+5,2*i+1)=so5
        y(j*8+6,2*i+1)=so6
        y(j*8+7,2*i+1)=so7
        y(j*8+8,2*i+1)=so8
     end do

  end do

!$omp  end do

!$omp  do  
  do j=(ndat/8)*8+1,ndat
     do i=0,n
      l=4 
        k=mod_arr(i-l)
        se=ch( 8)*x(  k,j)+cg( 8)*x(n+1+k  ,j)
      l=-4 
        k=mod_arr(i-l)
        so=ch(-7)*x(  k,j)+cg(-7)*x(n+1+k  ,j)
        do l=-3,3
           k=mod_arr(i-l)
           se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
           so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
        end do
        y(j,2*i  )=se
        y(j,2*i+1)=so
     end do

  end do

!$omp  end do
!$omp end parallel
!  call system_clock(ncount2,ncount_rate,ncount_max)
!  tel=dble(ncount2-ncount1)/dble(ncount_rate)
!  write(97,'(a40,1x,e10.3,1x,f6.1)') 'syn_rot_per:',tel,1.d-6*mflop/tel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on syn rot per',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE syn_rot_per


subroutine syn_rot_per_temp(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n-1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: so,se
  !       Daubechy S16
  real(wp), dimension(-6:9), parameter :: ch=(/&
                                    0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp,-0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, 0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, -0.051945838107881800736_wp, &
       0.36444189483617893676_wp, 0.77718575169962802862_wp, &
       0.48135965125905339159_wp, -0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, 0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, -0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp/)
  real(wp), dimension(-6:9), parameter :: cg=(/&
                                    -0.0033824159510050025955_wp, & 
       0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       -0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       -0.77718575169962802862_wp,0.36444189483617893676_wp, &
       0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       -0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       -0.0018899503327676891843_wp  /)
  
  do j=1,ndat

!!$     so=0.0_wp
!!$     se=0.0_wp
!!$     do l=-3,4
!!$        k=modulo(n-1+l,n)
!!$       ! se=se+ch(2*l)*x(k,j)+cg(2*l)*x(n+k,j)
!!$       ! so=so+ch(2*l+1)*x(k,j)+cg(2*l+1)*x(n+k,j)
!!$
!!$        se=se+ch(2*l)*x(k,j)+ch(l*-2+3)*x(n+k,j)
!!$        so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)
!!$
!!$     end do
!!$     y(j,2*n-1)=so
!!$     y(j,0)=se
!!$
!!$     do i=0,n-2
!!$        so=0.0_wp
!!$        se=0.0_wp
!!$        do l=-3,4
!!$           k=modulo(i+l,n)
!!$           !se=se+ch(2*l)*x(k,j)+cg(2*l)*x(n+k,j)
!!$           !so=so+ch(2*l+1)*x(k,j)+cg(2*l+1)*x(n+k,j)
!!$
!!$        se=se+ch(2*l)*x(k,j)+ch(-2*l+3)*x(n+k,j)
!!$        so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)
!!$
!!$        end do
!!$        y(j,2*i+1)=so
!!$        y(j,2*i+2)=se
!!$     end do

     so=0.0_wp
     se=0.0_wp
     do l=-3,4
        k=modulo(n-1+l,n)
        se=se+ch(2*l)*x(k,j)+ch(-l*2+3)*x(n+k,j)
        so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)

     end do
     y(j,2*n-1)=so
     y(j,0)=se

     do i=0,2
        so=0.0_wp
        se=0.0_wp
        do l=-3,4
           k=modulo(i+l,n)
           se=se+ch(2*l)*x(k,j)+ch(-2*l+3)*x(n+k,j)
           so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)

        end do
        y(j,2*i+1)=so
        y(j,2*i+2)=se
     end do
     do i=3,n-5
        so=0.0_wp
        se=0.0_wp
        do l=-3,4
           k=i+l
           se=se+ch(2*l)*x(k,j)+ch(-2*l+3)*x(n+k,j)
           so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)
        end do
        y(j,2*i+1)=so
        y(j,2*i+2)=se
     end do
     do i=n-5,n-2
        so=0.0_wp
        se=0.0_wp
        do l=-3,4
           k=modulo(i+l,n)
           se=se+ch(2*l)*x(k,j)+ch(-2*l+3)*x(n+k,j)
           so=so+ch(2*l+1)*x(k,j)-ch(-2*l+2)*x(n+k,j)

        end do
        y(j,2*i+1)=so
        y(j,2*i+2)=se
     end do
  end do

END SUBROUTINE syn_rot_per_temp

subroutine syn_rot_per_simple(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: so,se
  real(wp), dimension(-8:9) :: ch,cg
  !       Daubechy S16
  data ch  /  0.e0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.e0_wp /
  data cg  / 0.e0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.e0_wp /

  do j=1,ndat

     do i=0,n
        se=0.e0_wp
        so=0.e0_wp
        do l=-4,4
           k=modulo(i-l,n+1)
           se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
           so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

  enddo

END SUBROUTINE syn_rot_per_simple


subroutine convrot_n_per(n1,ndat,x,y)
  use module_base
  implicit none

!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test


  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) :: fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  integer :: mod_arr(lowfil:n1+lupfil)   
  integer :: i,j,l,k
  real(wp) :: fill,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8

  call fill_mod_arr(mod_arr,lowfil,n1+lupfil,n1+1)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,fil,ndat,n1,mod_arr)
!$omp do 
  do j=0,ndat/8-1

     do i=0,n1

        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
!        tt9 =0.e0_wp
!        tt10=0.e0_wp
!        tt11=0.e0_wp
!        tt12=0.e0_wp

        do l=lowfil,lupfil
           k=mod_arr(i+l)
           fill=fil(l)

           tt1=tt1+x(  k,j*8+1)*fill
           tt2=tt2+x(  k,j*8+2)*fill
           tt3=tt3+x(  k,j*8+3)*fill
           tt4=tt4+x(  k,j*8+4)*fill
           tt5=tt5+x(  k,j*8+5)*fill
           tt6=tt6+x(  k,j*8+6)*fill
           tt7=tt7+x(  k,j*8+7)*fill
           tt8=tt8+x(  k,j*8+8)*fill

!           tt9 =tt9 +x(  k,j*12+9 )*fill
!           tt10=tt10+x(  k,j*12+10)*fill
!           tt11=tt11+x(  k,j*12+11)*fill
!           tt12=tt12+x(  k,j*12+12)*fill
        end do
        y(j*8+1,i)=tt1
        y(j*8+2,i)=tt2
        y(j*8+3,i)=tt3
        y(j*8+4,i)=tt4
        y(j*8+5,i)=tt5
        y(j*8+6,i)=tt6
        y(j*8+7,i)=tt7
        y(j*8+8,i)=tt8

!        y(j*12+9 ,i)=tt9 
!        y(j*12+10,i)=tt10
!        y(j*12+11,i)=tt11
!        y(j*12+12,i)=tt12

     end do
  end do

  !$omp end do

  !$omp do

  do j=(ndat/8)*8+1,ndat
     do i=0,n1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr(i+l)   
           tt=tt+x(  k,j)*fil(l)
        end do
        y(j,i)=tt

     end do
  end do

  !$omp end do
  !$omp end parallel


!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on convrot n per',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE convrot_n_per


subroutine convrot_t_per(n1,ndat,x,y)
  use module_base
  implicit none

!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

  integer,parameter::lowfil=-7,lupfil=8
  !  dimension x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) :: fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228D-6,&
       -0.5185986881173432922848639136911487D-4,&
       0.49443227688689919192282259476750972D-3,&
       -0.344128144493493857280881509686821861D-2,&
       0.1337263414854794752733423467013220997D-1,&
       -0.2103025160930381434955489412839065067D-1,&
       -0.604895289196983516002834636D-1,&
       0.9940415697834003993178616713D0,&
       0.612625895831207982195380597D-1,&
       0.2373821463724942397566389712597274535D-1,&
       -0.942047030201080385922711540948195075D-2,&
       0.174723713672993903449447812749852942D-2,&
       -0.30158038132690463167163703826169879D-3,&
       0.8762984476210559564689161894116397D-4,&
       -0.1290557201342060969516786758559028D-4,&
       8.4334247333529341094733325815816D-7 /

  integer :: i,j,l,k
  integer :: mod_arr(lowfil:n1+lupfil)   
  real(wp) :: fill,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,tt

  call fill_mod_arr(mod_arr,lowfil,n1+lupfil,n1+1)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,fil,ndat,n1,mod_arr)
!$omp do 
  do j=0,ndat/12-1

     do i=0,n1

        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           k=mod_arr(i+l)   
           fill=fil(l)

           tt1=tt1+x(  k,j*12+1)*fill
           tt2=tt2+x(  k,j*12+2)*fill
           tt3=tt3+x(  k,j*12+3)*fill
           tt4=tt4+x(  k,j*12+4)*fill
           tt5=tt5+x(  k,j*12+5)*fill
           tt6=tt6+x(  k,j*12+6)*fill
           tt7=tt7+x(  k,j*12+7)*fill
           tt8=tt8+x(  k,j*12+8)*fill

           tt9 =tt9 +x(  k,j*12+9 )*fill
           tt10=tt10+x(  k,j*12+10)*fill
           tt11=tt11+x(  k,j*12+11)*fill
           tt12=tt12+x(  k,j*12+12)*fill
        end do
        y(j*12+1,i)=tt1
        y(j*12+2,i)=tt2
        y(j*12+3,i)=tt3
        y(j*12+4,i)=tt4
        y(j*12+5,i)=tt5
        y(j*12+6,i)=tt6
        y(j*12+7,i)=tt7
        y(j*12+8,i)=tt8

        y(j*12+9 ,i)=tt9 
        y(j*12+10,i)=tt10
        y(j*12+11,i)=tt11
        y(j*12+12,i)=tt12

     end do
  end do

  !$omp end do

  !$omp do
  do j=(ndat/12)*12+1,ndat
     do i=0,n1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr(i+l)   
           tt=tt+x(  k,j)*fil(l)
        end do
        y(j,i)=tt

     end do
  end do
  !$omp end do

  !$omp end parallel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE convrot_t_per


subroutine convrot_t_per_test(n1,ndat,x,y)
  use module_base
  implicit none

!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

  integer,parameter::lowfil=-7,lupfil=8
  !  dimension x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) :: fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228D-6,&
       -0.5185986881173432922848639136911487D-4,&
       0.49443227688689919192282259476750972D-3,&
       -0.344128144493493857280881509686821861D-2,&
       0.1337263414854794752733423467013220997D-1,&
       -0.2103025160930381434955489412839065067D-1,&
       -0.604895289196983516002834636D-1,&
       0.9940415697834003993178616713D0,&
       0.612625895831207982195380597D-1,&
       0.2373821463724942397566389712597274535D-1,&
       -0.942047030201080385922711540948195075D-2,&
       0.174723713672993903449447812749852942D-2,&
       -0.30158038132690463167163703826169879D-3,&
       0.8762984476210559564689161894116397D-4,&
       -0.1290557201342060969516786758559028D-4,&
       8.4334247333529341094733325815816D-7 /

  integer :: i,j,l,k
  integer :: mod_arr(lowfil:n1+lupfil)   
  real(wp) :: fill,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,tt

  call fill_mod_arr(mod_arr,lowfil,n1+lupfil,n1+1)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,fil,ndat,n1,mod_arr)
!$omp do 
  do j=0,ndat/12-1

     do i=0,n1

        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           k=mod_arr(i+l)   
           fill=fil(l)

           tt1=tt1+x(  k,j*12+1)*fill
           tt2=tt2+x(  k,j*12+2)*fill
           tt3=tt3+x(  k,j*12+3)*fill
           tt4=tt4+x(  k,j*12+4)*fill
           tt5=tt5+x(  k,j*12+5)*fill
           tt6=tt6+x(  k,j*12+6)*fill
           tt7=tt7+x(  k,j*12+7)*fill
           tt8=tt8+x(  k,j*12+8)*fill

           tt9 =tt9 +x(  k,j*12+9 )*fill
           tt10=tt10+x(  k,j*12+10)*fill
           tt11=tt11+x(  k,j*12+11)*fill
           tt12=tt12+x(  k,j*12+12)*fill
        end do
        y(j*12+1,i)=tt1
        y(j*12+2,i)=tt2
        y(j*12+3,i)=tt3
        y(j*12+4,i)=tt4
        y(j*12+5,i)=tt5
        y(j*12+6,i)=tt6
        y(j*12+7,i)=tt7
        y(j*12+8,i)=tt8

        y(j*12+9 ,i)=tt9 
        y(j*12+10,i)=tt10
        y(j*12+11,i)=tt11
        y(j*12+12,i)=tt12

     end do
  end do

  !$omp end do

  !$omp do
  do j=(ndat/12)*12+1,ndat
     do i=0,n1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=mod_arr(i+l)   
           tt=tt+x(  k,j)*fil(l)
        end do
        y(j,i)=tt

     end do
  end do
  !$omp end do

  !$omp end parallel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on comb',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE convrot_t_per_test


!> Applies the kinetic energy operator onto x to get y. Works for periodic BC
!! This routines is used by OpenCL/conv_check.f90
subroutine convolut_kinetic_per_c(n1,n2,n3,hgrid,x,y,c)
  use module_base
  implicit none

  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y

!  stop 'convolut_kinetic_per_c should never be called'

  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer, dimension(lowfil:n1+lupfil) :: mod_arr1   
  integer, dimension(lowfil:n2+lupfil) :: mod_arr2   
  integer, dimension(lowfil:n3+lupfil) :: mod_arr3   
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  integer :: k

  !$omp parallel default(private) shared(x,y,n1,n2,n3,c,hgrid,fil,mod_arr1,mod_arr2,mod_arr3)
  call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do k=1,14
     fil(-k,:)=fil(k,:)
  end do

  call conv_kin_x(mod_arr1,fil,c,x,y,n1,(n2+1)*(n3+1))   
  call conv_kin_y(mod_arr2,fil,x,y,n1,n2,n3)
  call conv_kin_z(mod_arr3,fil,x,y,n3,(n1+1)*(n2+1))
  !$omp end parallel


  contains


  subroutine conv_kin_x(mod_arr,fil0,c0,x0,y0,n,ndat)
    implicit none
    integer, intent(in) :: n,ndat
    integer, dimension(lowfil:n+lupfil) :: mod_arr   
    real(wp), dimension(lowfil:lupfil,3) :: fil0
    real(gp), intent(in) :: c0
    real(wp), intent(in) :: x0(0:n,ndat)
    real(wp), intent(out) :: y0(0:n,ndat)
    integer :: i1,i,l,j
    real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12

    !$omp do 
    do i=0,ndat/12-1
       do i1=0,n1
          tt1 =x0(i1,i*12+1)*c0
          tt2 =x0(i1,i*12+2)*c0
          tt3 =x0(i1,i*12+3)*c0
          tt4 =x0(i1,i*12+4)*c0
          tt5 =x0(i1,i*12+5)*c0
          tt6 =x0(i1,i*12+6)*c0
          tt7 =x0(i1,i*12+7)*c0
          tt8 =x0(i1,i*12+8)*c0
          tt9 =x0(i1,i*12+9 )*c0
          tt10=x0(i1,i*12+10)*c0
          tt11=x0(i1,i*12+11)*c0
          tt12=x0(i1,i*12+12)*c0

          do l=lowfil,lupfil
             j=mod_arr(i1+l)

             tt1=tt1+x0(j,i*12+1)*fil0(l,1)
             tt2=tt2+x0(j,i*12+2)*fil0(l,1)
             tt3=tt3+x0(j,i*12+3)*fil0(l,1)
             tt4=tt4+x0(j,i*12+4)*fil0(l,1)
             tt5=tt5+x0(j,i*12+5)*fil0(l,1)
             tt6=tt6+x0(j,i*12+6)*fil0(l,1)
             tt7=tt7+x0(j,i*12+7)*fil0(l,1)
             tt8=tt8+x0(j,i*12+8)*fil0(l,1)
             tt9 =tt9 +x0(j,i*12+9 )*fil0(l,1)
             tt10=tt10+x0(j,i*12+10)*fil0(l,1)
             tt11=tt11+x0(j,i*12+11)*fil0(l,1)
             tt12=tt12+x0(j,i*12+12)*fil0(l,1)
          end do
          y0(i1,i*12+1)=tt1
          y0(i1,i*12+2)=tt2
          y0(i1,i*12+3)=tt3
          y0(i1,i*12+4)=tt4
          y0(i1,i*12+5)=tt5
          y0(i1,i*12+6)=tt6
          y0(i1,i*12+7)=tt7
          y0(i1,i*12+8)=tt8
          y0(i1,i*12+9 )=tt9 
          y0(i1,i*12+10)=tt10
          y0(i1,i*12+11)=tt11
          y0(i1,i*12+12)=tt12
       end do
    end do
    !$omp end do

    !$omp do 
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=x0(i1,i)*c0
          do l=lowfil,lupfil
             j=mod_arr(i1+l)
             tt=tt+x0(j   ,i)*fil0(l,1)
          end do
          y0(i1,i)=tt
       end do
    end do
    !$omp end do
  END SUBROUTINE conv_kin_x


  subroutine conv_kin_y(mod_arr,fil0,x0,y0,n1,n2,n3)
    implicit none
    !Arguments
    integer, intent(in) :: n1,n2,n3
    integer, dimension(lowfil:n2+lupfil), intent(in) :: mod_arr
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
    real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x0
    real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y0
    integer :: i1,i2,i3,l,j
    real(wp) :: tt,tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7

    !$omp do
    do i3=0,n3/8-1
        do i1=0,n1
           do i2=0,n2
              tt0=0.e0_wp
              tt1=0.e0_wp
              tt2=0.e0_wp
              tt3=0.e0_wp
              tt4=0.e0_wp
              tt5=0.e0_wp
              tt6=0.e0_wp
              tt7=0.e0_wp
    
              do l=lowfil,lupfil
                 j=mod_arr(i2+l)
    
                 tt0=tt0+x0(i1,j,i3*8+0)*fil0(l,2)
                 tt1=tt1+x0(i1,j,i3*8+1)*fil0(l,2)
                 tt2=tt2+x0(i1,j,i3*8+2)*fil0(l,2)
                 tt3=tt3+x0(i1,j,i3*8+3)*fil0(l,2)
                 tt4=tt4+x0(i1,j,i3*8+4)*fil0(l,2)
                 tt5=tt5+x0(i1,j,i3*8+5)*fil0(l,2)
                 tt6=tt6+x0(i1,j,i3*8+6)*fil0(l,2)
                 tt7=tt7+x0(i1,j,i3*8+7)*fil0(l,2)
              end do
              y0(i1,i2,i3*8+0)=y0(i1,i2,i3*8+0)+tt0
              y0(i1,i2,i3*8+1)=y0(i1,i2,i3*8+1)+tt1
              y0(i1,i2,i3*8+2)=y0(i1,i2,i3*8+2)+tt2
              y0(i1,i2,i3*8+3)=y0(i1,i2,i3*8+3)+tt3
              y0(i1,i2,i3*8+4)=y0(i1,i2,i3*8+4)+tt4
              y0(i1,i2,i3*8+5)=y0(i1,i2,i3*8+5)+tt5
              y0(i1,i2,i3*8+6)=y0(i1,i2,i3*8+6)+tt6
              y0(i1,i2,i3*8+7)=y0(i1,i2,i3*8+7)+tt7
           end do
        end do
     end do
     !$omp end do

     !$omp do 
     do i3=(n3/8)*8,n3
        do i1=0,n1
           do i2=0,n2
              tt=0.e0_wp
              do l=lowfil,lupfil
                 j=mod_arr(i2+l)
                 tt=tt+x0(i1,j   ,i3)*fil0(l,2)
              end do
              y0(i1,i2,i3)=y0(i1,i2,i3)+tt
           end do
        end do
     end do
     !$omp end do
  END SUBROUTINE conv_kin_y


  subroutine conv_kin_z(mod_arr,fil0,x0,y0,n,ndat)
    implicit none
    integer, intent(in) :: n,ndat
    integer, dimension(lowfil:n+lupfil), intent(in) :: mod_arr   
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
    real(wp),intent(in):: x0(ndat,0:n)
    real(wp),intent(inout)::y0(ndat,0:n)
    integer :: i3,i,l,j
    real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12

    !$omp do 
    do i=0,ndat/12-1
       do i3=0,n3
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
          tt5=0.e0_wp
          tt6=0.e0_wp
          tt7=0.e0_wp
          tt8=0.e0_wp
          tt9 =0.e0_wp
          tt10=0.e0_wp
          tt11=0.e0_wp
          tt12=0.e0_wp

          do l=lowfil,lupfil
             j=mod_arr(i3+l)

             tt1=tt1+x0(i*12+1,j)*fil0(l,3)
             tt2=tt2+x0(i*12+2,j)*fil0(l,3)
             tt3=tt3+x0(i*12+3,j)*fil0(l,3)
             tt4=tt4+x0(i*12+4,j)*fil0(l,3)
             tt5=tt5+x0(i*12+5,j)*fil0(l,3)
             tt6=tt6+x0(i*12+6,j)*fil0(l,3)
             tt7=tt7+x0(i*12+7,j)*fil0(l,3)
             tt8=tt8+x0(i*12+8,j)*fil0(l,3)
             tt9 =tt9 +x0(i*12+9 ,j)*fil0(l,3)
             tt10=tt10+x0(i*12+10,j)*fil0(l,3)
             tt11=tt11+x0(i*12+11,j)*fil0(l,3)
             tt12=tt12+x0(i*12+12,j)*fil0(l,3)
          end do

          y0(i*12+1,i3)=y0(i*12+1,i3)+tt1
          y0(i*12+2,i3)=y0(i*12+2,i3)+tt2
          y0(i*12+3,i3)=y0(i*12+3,i3)+tt3
          y0(i*12+4,i3)=y0(i*12+4,i3)+tt4
          y0(i*12+5,i3)=y0(i*12+5,i3)+tt5
          y0(i*12+6,i3)=y0(i*12+6,i3)+tt6
          y0(i*12+7,i3)=y0(i*12+7,i3)+tt7
          y0(i*12+8,i3)=y0(i*12+8,i3)+tt8
          y0(i*12+9 ,i3)=y0(i*12+9 ,i3)+tt9 
          y0(i*12+10,i3)=y0(i*12+10,i3)+tt10
          y0(i*12+11,i3)=y0(i*12+11,i3)+tt11
          y0(i*12+12,i3)=y0(i*12+12,i3)+tt12
       end do
    end do
    !$omp end do

    !$omp do 
    do i=(ndat/12)*12+1,ndat
       do i3=0,n3
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr(i3+l)
             tt=tt+x0(i,j)*fil0(l,3)
          end do
          y0(i,i3)=y0(i,i3)+tt
       end do
    end do
    !$omp end do
  END SUBROUTINE conv_kin_z


END SUBROUTINE convolut_kinetic_per_c


subroutine convolut_kinetic_per_T(n1,n2,n3,hgrid,x,y,kstrten)
  !   applies the kinetic energy operator onto x to get y. Works for periodic BC
  !   y:=y-1/2Delta x
  use module_base
  implicit none

  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !real(wp),intent(out)::ekin_out
  real(wp), dimension(6), intent(out) :: kstrten
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: k!,i !for non OMP case
  integer, dimension(lowfil:n1+lupfil) :: mod_arr1
  integer, dimension(lowfil:n2+lupfil) :: mod_arr2
  integer, dimension(lowfil:n3+lupfil) :: mod_arr3
  real(wp) :: ekin1,ekin2,ekin3
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  !real(wp), dimension(8,3) :: ekin_array
!$ integer :: ithread=0
!$ integer :: omp_get_thread_num

  !ekin_out=0._wp
  kstrten=0.0_gp

  !do i=1,8
  !ekin_array(i,1)=10._wp
  !ekin_array(i,2)=10._wp
  !ekin_array(i,3)=10._wp
  !end do
 !$omp parallel default(private) shared(x,y,n1,n2,n3,hgrid,kstrten,fil,mod_arr1,mod_arr2,mod_arr3)
!$  ithread = omp_get_thread_num()
  call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do k=1,14
     fil(-k,:)=fil(k,:)
  end do

  ekin2=0.0_wp
  ekin3=0.0_wp
  ekin1=0.0_wp
  call conv_kin_x(x,y,n1,n2,n3,ekin1,fil,mod_arr1)
  call conv_kin_y(x,y,n1,n2,n3,ekin2,fil,mod_arr2)
  call conv_kin_z(x,y,n1,n2,n3,ekin3,fil,mod_arr3)
  !ekin_array(ithread+1,1)=ekin1
  !ekin_array(ithread+1,2)=ekin2
  !ekin_array(ithread+1,3)=ekin3
  !$omp critical 

!yk
kstrten(1)=kstrten(1)+ekin1
kstrten(2)=kstrten(2)+ekin2
kstrten(3)=kstrten(3)+ekin3

!     ekin_out=ekin_out+ekin1+ekin2+ekin3
  !$omp end critical
  !$omp end parallel


!dee
!!$!  open(unit=97,file='check_ekin3',status='unknown')
!!$    write(197,*) '-------------------------------------------------------------------'
!!$do i=1,8
!!$    write(197,'(3(1X,e24.17))') ekin_array(i,1),ekin_array(i,2),ekin_array(i,3)
!!$end do
!  close(97)

END SUBROUTINE convolut_kinetic_per_T


subroutine conv_kin_x(x,y,n1,n2,n3,ekin,fil,mod_arr1)
  use module_base
  implicit none
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mod_arr1(lowfil:n1+lupfil)
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  integer :: ndat
  real(wp),intent(in) :: x(0:n1,(n2+1)*(n3+1))
  real(wp),intent(inout) :: y(0:n1,(n2+1)*(n3+1))
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
  integer :: i,l,j,i1

  ndat=(n2+1)*(n3+1)
  !$omp do
  do i=0,ndat/12-1
     do i1=0,n1
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           j=mod_arr1(i1+l)

           tt1=tt1+x(j,i*12+1)*fil(l,1)
           tt2=tt2+x(j,i*12+2)*fil(l,1)
           tt3=tt3+x(j,i*12+3)*fil(l,1)
           tt4=tt4+x(j,i*12+4)*fil(l,1)
           tt5=tt5+x(j,i*12+5)*fil(l,1)
           tt6=tt6+x(j,i*12+6)*fil(l,1)
           tt7=tt7+x(j,i*12+7)*fil(l,1)
           tt8=tt8+x(j,i*12+8)*fil(l,1)
           tt9 =tt9 +x(j,i*12+9 )*fil(l,1)
           tt10=tt10+x(j,i*12+10)*fil(l,1)
           tt11=tt11+x(j,i*12+11)*fil(l,1)
           tt12=tt12+x(j,i*12+12)*fil(l,1)
        end do
        y(i1,i*12+1)=y(i1,i*12+1)+tt1;    ekin=ekin+tt1*x(i1,i*12+1)
        y(i1,i*12+2)=y(i1,i*12+2)+tt2;    ekin=ekin+tt2*x(i1,i*12+2)
        y(i1,i*12+3)=y(i1,i*12+3)+tt3;    ekin=ekin+tt3*x(i1,i*12+3)
        y(i1,i*12+4)=y(i1,i*12+4)+tt4;    ekin=ekin+tt4*x(i1,i*12+4)
        y(i1,i*12+5)=y(i1,i*12+5)+tt5;    ekin=ekin+tt5*x(i1,i*12+5)
        y(i1,i*12+6)=y(i1,i*12+6)+tt6;    ekin=ekin+tt6*x(i1,i*12+6)
        y(i1,i*12+7)=y(i1,i*12+7)+tt7;    ekin=ekin+tt7*x(i1,i*12+7)
        y(i1,i*12+8)=y(i1,i*12+8)+tt8;    ekin=ekin+tt8*x(i1,i*12+8)
        y(i1,i*12+9 )=y(i1,i*12+9 )+tt9 ;    ekin=ekin+tt9 *x(i1,i*12+9 )
        y(i1,i*12+10)=y(i1,i*12+10)+tt10;    ekin=ekin+tt10*x(i1,i*12+10)
        y(i1,i*12+11)=y(i1,i*12+11)+tt11;    ekin=ekin+tt11*x(i1,i*12+11)
        y(i1,i*12+12)=y(i1,i*12+12)+tt12;    ekin=ekin+tt12*x(i1,i*12+12)
       end do
    end do
   !$omp end do

   !$omp do
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr1(i1+l)
             tt=tt+x(j   ,i)*fil(l,1)
          end do
          y(i1,i)=y(i1,i)+tt ; ekin=ekin+tt*x(i1,i)
       end do
    end do
   !$omp end do
END SUBROUTINE conv_kin_x

subroutine conv_kin_y(x,y,n1,n2,n3,ekin,fil,mod_arr2)
  use module_base
  implicit none
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mod_arr2(lowfil:n2+lupfil)
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  real(wp),intent(in) :: x(0:n1,0:n2,0:n3)
  real(wp),intent(inout) :: y(0:n1,0:n2,0:n3)
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7
  integer :: l,j,i1,i2,i3

  !$omp do
  do i3=0,n3/8-1
     do i1=0,n1
        do i2=0,n2
           tt0=0.e0_wp
           tt1=0.e0_wp
           tt2=0.e0_wp
           tt3=0.e0_wp
           tt4=0.e0_wp
           tt5=0.e0_wp
           tt6=0.e0_wp
           tt7=0.e0_wp

           do l=lowfil,lupfil
              j=mod_arr2(i2+l)

              tt0=tt0+x(i1,j,i3*8+0)*fil(l,2)
              tt1=tt1+x(i1,j,i3*8+1)*fil(l,2)
              tt2=tt2+x(i1,j,i3*8+2)*fil(l,2)
              tt3=tt3+x(i1,j,i3*8+3)*fil(l,2)
              tt4=tt4+x(i1,j,i3*8+4)*fil(l,2)
              tt5=tt5+x(i1,j,i3*8+5)*fil(l,2)
              tt6=tt6+x(i1,j,i3*8+6)*fil(l,2)
              tt7=tt7+x(i1,j,i3*8+7)*fil(l,2)
           end do
           y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;    ekin=ekin+tt0*x(i1,i2,i3*8+0)
           y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;    ekin=ekin+tt1*x(i1,i2,i3*8+1)
           y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;    ekin=ekin+tt2*x(i1,i2,i3*8+2)
           y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;    ekin=ekin+tt3*x(i1,i2,i3*8+3)
           y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;    ekin=ekin+tt4*x(i1,i2,i3*8+4)
           y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;    ekin=ekin+tt5*x(i1,i2,i3*8+5)
           y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;    ekin=ekin+tt6*x(i1,i2,i3*8+6)
           y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;    ekin=ekin+tt7*x(i1,i2,i3*8+7)
        end do
     end do
  end do
   
  !$omp end do

  !$omp do
  do i3=(n3/8)*8,n3
     do i1=0,n1
        do i2=0,n2
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=mod_arr2(i2+l)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           end do
           y(i1,i2,i3)=y(i1,i2,i3)+tt
           ekin=ekin+tt*x(i1,i2,i3)
        end do
     end do
  end do
  !$omp end do
END SUBROUTINE conv_kin_y


subroutine conv_kin_z(x,y,n1,n2,n3,ekin,fil,mod_arr3)
  use module_base
  implicit none
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  integer, intent(in) :: mod_arr3(lowfil:n3+lupfil)
  real(wp),intent(in) :: x((n2+1)*(n1+1),0:n3)
  real(wp),intent(inout) :: y((n2+1)*(n1+1),0:n3)
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
  integer :: ndat
  integer ::i,l,j,i3

  ndat=(n2+1)*(n1+1)
  !$omp do
  do i=0,ndat/12-1
     do i3=0,n3
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           j=mod_arr3(i3+l)

           !print *,'j,mod_arr',j,'aa',mod_arr3(:)

           tt1=tt1+x(i*12+1,j)*fil(l,3)
           tt2=tt2+x(i*12+2,j)*fil(l,3)
           tt3=tt3+x(i*12+3,j)*fil(l,3)
           tt4=tt4+x(i*12+4,j)*fil(l,3)
           tt5=tt5+x(i*12+5,j)*fil(l,3)
           tt6=tt6+x(i*12+6,j)*fil(l,3)
           tt7=tt7+x(i*12+7,j)*fil(l,3)
           tt8=tt8+x(i*12+8,j)*fil(l,3)
           tt9 =tt9 +x(i*12+9 ,j)*fil(l,3)
           tt10=tt10+x(i*12+10,j)*fil(l,3)
           tt11=tt11+x(i*12+11,j)*fil(l,3)
           tt12=tt12+x(i*12+12,j)*fil(l,3)
        end do

        y(i*12+1,i3)=y(i*12+1,i3)+tt1;    ekin=ekin+tt1*x(i*12+1,i3)
        y(i*12+2,i3)=y(i*12+2,i3)+tt2;    ekin=ekin+tt2*x(i*12+2,i3)
        y(i*12+3,i3)=y(i*12+3,i3)+tt3;    ekin=ekin+tt3*x(i*12+3,i3)
        y(i*12+4,i3)=y(i*12+4,i3)+tt4;    ekin=ekin+tt4*x(i*12+4,i3)
        y(i*12+5,i3)=y(i*12+5,i3)+tt5;    ekin=ekin+tt5*x(i*12+5,i3)
        y(i*12+6,i3)=y(i*12+6,i3)+tt6;    ekin=ekin+tt6*x(i*12+6,i3)
        y(i*12+7,i3)=y(i*12+7,i3)+tt7;    ekin=ekin+tt7*x(i*12+7,i3)
        y(i*12+8,i3)=y(i*12+8,i3)+tt8;    ekin=ekin+tt8*x(i*12+8,i3)
        y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 ;    ekin=ekin+tt9*x(i*12+9 ,i3)
        y(i*12+10,i3)=y(i*12+10,i3)+tt10;    ekin=ekin+tt10*x(i*12+10,i3)
        y(i*12+11,i3)=y(i*12+11,i3)+tt11;    ekin=ekin+tt11*x(i*12+11,i3)
        y(i*12+12,i3)=y(i*12+12,i3)+tt12;    ekin=ekin+tt12*x(i*12+12,i3)
     end do
  end do
  !$omp end do

  !$omp do
  do i=(ndat/12)*12+1,ndat
     do i3=0,n3
        tt=0.e0_wp
        do l=lowfil,lupfil
           j=mod_arr3(i3+l)
           tt=tt+x(i,j)*fil(l,3)
        end do
        y(i,i3)=y(i,i3)+tt; ekin=ekin+tt*x(i,i3)
     end do
  end do
  !$omp end do
END SUBROUTINE conv_kin_z


subroutine fill_mod_arr(arr,nleft,nright,n)
  implicit none
  integer,intent(in) :: nleft,nright,n
  integer,intent(out) :: arr(nleft:nright)
  integer :: i
  
  if (nleft >= -n) then
     do i=nleft,-1
        arr(i)=n+i
     end do
  else
     do i=nleft,-1
        arr(i)=modulo(i,n)
     end do
  endif
  
  do i=max(0,nleft),min(n-1,nright)
     arr(i)=i
  end do
  
  if (nright < 2*n) then
     do i=n,nright
        arr(i)=i-n
     end do
  else
     do i=n,nright
        arr(i)=modulo(i,n)
     end do
  endif
END SUBROUTINE fill_mod_arr


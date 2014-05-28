!> @file
!!  Optimized convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Forward wavelet transform, analysis, periodic
subroutine ana_rot_per_old(right,nt,c,cd_1)
  
  use module_base
  implicit none
  integer, intent(in) :: right,nt
  real(wp), dimension(0:right,nt), intent(in) :: c
  real(wp), dimension(nt,0:right), intent(out) :: cd_1
  !local variables
  character(len=*), parameter :: subname='ana_rot_per_old'
  integer, parameter :: m=8
  integer :: i_all,i_stat,lenc,len_2,mod_left,mod_right,i,it,i2,it0,j,ji2,il2
!$  integer :: ithread,omp_get_thread_num
  real(wp) :: ci_0,ci_1,ci_2, ci_3,ci_4,ci_5,ci_6,ci_7,ci_8,ci_9,ci_10,ci_11,ci,cgj,chj
  real(wp) :: di_0,di_1,di_2, di_3,di_4,di_5,di_6,di_7,di_8,di_9,di_10,di_11,di
  integer, dimension(:), allocatable :: mod_my
  real(wp) ch(-8:9) ,cg(-8:9)
  !       daubechy s16
  data ch  /  0.0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.0_wp /
  data cg  / 0.0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.0_wp /

!   write(*,*) 'ana_rot_per_old EXECUTED',right,nt,c(0,1),c(right,nt)

  lenc=right+1
  len_2=lenc/2

  mod_left=1-m
  mod_right=2*len_2-2+m

  mod_my = f_malloc(mod_left.to.mod_right,id='mod_my')

  do i=mod_left,mod_right
     mod_my(i)=modulo(i,lenc)
  enddo

  !       nflop=nt*len_2*2*m*4
  !       call system_clock(ncount1,ncount_rate,ncount_max)

!$omp parallel default (private) shared(nt,len_2,ch,cg,mod_my,c,cd_1)
!$  ithread = omp_get_thread_num()

!$omp do schedule(static,1)
  do it=1,nt-11,12
     do i=0,len_2-1
        i2=2*i

        ci_0 =0.0_wp
        ci_1 =0.0_wp
        ci_2 =0.0_wp
        ci_3 =0.0_wp
        ci_4 =0.0_wp
        ci_5 =0.0_wp
        ci_6 =0.0_wp
        ci_7 =0.0_wp
        ci_8 =0.0_wp
        ci_9 =0.0_wp
        ci_10=0.0_wp
        ci_11=0.0_wp

        di_0 =0.0_wp
        di_1 =0.0_wp
        di_2 =0.0_wp
        di_3 =0.0_wp
        di_4 =0.0_wp
        di_5 =0.0_wp
        di_6 =0.0_wp
        di_7 =0.0_wp
        di_8 =0.0_wp
        di_9 =0.0_wp
        di_10=0.0_wp
        di_11=0.0_wp

        do j=1-m,m
           ji2=mod_my(j+i2)

           chj=ch(j)
           cgj=cg(j) 

           ci_0 =ci_0 +chj*c(ji2,it+0 )
           ci_1 =ci_1 +chj*c(ji2,it+1 )
           ci_2 =ci_2 +chj*c(ji2,it+2 )
           ci_3 =ci_3 +chj*c(ji2,it+3 )

           ci_4 =ci_4 +chj*c(ji2,it+4 )
           ci_5 =ci_5 +chj*c(ji2,it+5 )
           ci_6 =ci_6 +chj*c(ji2,it+6 )
           ci_7 =ci_7 +chj*c(ji2,it+7 )
           ci_8 =ci_8  +chj*c(ji2,it+8  )
           ci_9 =ci_9  +chj*c(ji2,it+9  )
           ci_10=ci_10 +chj*c(ji2,it+10 )
           ci_11=ci_11 +chj*c(ji2,it+11 )

           di_0 =di_0 +cgj*c(ji2,it+0 )
           di_1 =di_1 +cgj*c(ji2,it+1 )
           di_2 =di_2 +cgj*c(ji2,it+2 )
           di_3 =di_3 +cgj*c(ji2,it+3 )

           di_4 =di_4 +cgj*c(ji2,it+4 )
           di_5 =di_5 +cgj*c(ji2,it+5 )
           di_6 =di_6 +cgj*c(ji2,it+6 )
           di_7 =di_7 +cgj*c(ji2,it+7 )
           di_8 =di_8 +cgj*c(ji2,it+8 )
           di_9 =di_9 +cgj*c(ji2,it+9 )
           di_10=di_10+cgj*c(ji2,it+10)
           di_11=di_11+cgj*c(ji2,it+11)

        enddo

        cd_1(it+0,i)=ci_0
        cd_1(it+1,i)=ci_1
        cd_1(it+2,i)=ci_2
        cd_1(it+3,i)=ci_3
        cd_1(it+4,i)=ci_4
        cd_1(it+5,i)=ci_5
        cd_1(it+6,i)=ci_6
        cd_1(it+7,i)=ci_7
        cd_1(it+8 ,i)=ci_8 
        cd_1(it+9 ,i)=ci_9 
        cd_1(it+10,i)=ci_10
        cd_1(it+11,i)=ci_11

        il2=len_2+i

        cd_1(it+0,il2)=di_0
        cd_1(it+1,il2)=di_1
        cd_1(it+2,il2)=di_2
        cd_1(it+3,il2)=di_3
        cd_1(it+4,il2)=di_4
        cd_1(it+5,il2)=di_5
        cd_1(it+6,il2)=di_6
        cd_1(it+7,il2)=di_7
        cd_1(it+8 ,il2)=di_8 
        cd_1(it+9 ,il2)=di_9 
        cd_1(it+10,il2)=di_10
        cd_1(it+11,il2)=di_11

     enddo
  enddo
!$omp enddo
  it0=it
!$ if (ithread.eq.0) then
!$  it0=12*int(nt/12)+1
  do it=it0,nt
     do i=0,len_2-1
        i2=2*i
        ci=0.0_wp
        di=0.0_wp
        do j=1-m,m
           ji2=mod_my(j+i2)
           ci=ci+ch(j)*c(ji2,it)
           di=di+cg(j)*c(ji2,it)
        enddo
        cd_1(it,i)=ci
        cd_1(it,len_2+i)=di
     enddo
  enddo
!$ endif
!$omp end parallel
  !        call system_clock(ncount2,ncount_rate,ncount_max)
  !        tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !        write(95,'(a40,1x,e11.4,1x,f10.1,1x,i9)') 'ana_rot_per_old',tel,1.d-6*nflop/tel,nflop
  call f_free(mod_my)

!write(*,*) 'ana_rot_per_old finished'

END SUBROUTINE ana_rot_per_old


!> Backward wavelet transform, synthesis, periodic
subroutine syn_rot_per_old(right1,nt,cd,c1)

  use module_base
  implicit none
  integer, intent(in) :: right1,nt
  real(wp), dimension(0:right1,nt), intent(in) :: cd
  real(wp), dimension(nt,0:right1), intent(out) :: c1
  !local variables
  character(len=*), parameter :: subname='syn_rot_per_old'
  !n(c) integer, parameter :: m=8
  integer, parameter :: m_2=4
  integer :: i_all,i_stat,len_2,mod_left,mod_right,i,it,it0,i2,i_j,j,i21,i_j2,j2,j21
!$  integer :: ithread,omp_get_thread_num
  real(wp) :: ci2_0,ci2_1,ci2_2, ci2_3,ci2_4,ci2_5,ci2_6,ci2_7,ci2_8,ci2_9,ci2_10,ci2_11,ci2
  real(wp) :: ci21_0,ci21_1,ci21_2, ci21_3,ci21_4,ci21_5,ci21_6,ci21_7,ci21_8,ci21_9,ci21_10
  real(wp) :: ci21_11,ci21,cgj2,chj2,cgj21,chj21
  integer, dimension(:), allocatable :: mod_my
  real(wp) ch(-8:9) ,cg(-8:9)
  !       daubechy s16
  data ch  /  0.0_wp , -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp , 0.0_wp /
  data cg  / 0.0_wp , -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp , 0.0_wp /

!  m_2=m/2
  len_2=(right1+1)/2

  mod_left=-m_2
  mod_right=len_2-1+m_2

  mod_my = f_malloc(mod_left.to.mod_right,id='mod_my')


  do i=mod_left,mod_right
     mod_my(i)=modulo(i,len_2)
  enddo

  !       nflop=nt*len_2*(2*m_2+1)*8
  !       call system_clock(ncount1,ncount_rate,ncount_max)

!$omp parallel default (private) shared(nt,len_2,ch,cg,mod_my,cd,c1)
!$  ithread = omp_get_thread_num()


!$omp do schedule(static,1)
  do it=1,nt-11,12
     do i=0,len_2-1

        ci2_0 =0.0_wp
        ci2_1 =0.0_wp
        ci2_2 =0.0_wp
        ci2_3 =0.0_wp
        ci2_4 =0.0_wp
        ci2_5 =0.0_wp
        ci2_6 =0.0_wp
        ci2_7 =0.0_wp
        ci2_8  =0.0_wp
        ci2_9  =0.0_wp
        ci2_10 =0.0_wp
        ci2_11 =0.0_wp

        ci21_0=0.0_wp
        ci21_1=0.0_wp
        ci21_2=0.0_wp
        ci21_3=0.0_wp
        ci21_4=0.0_wp
        ci21_5=0.0_wp
        ci21_6=0.0_wp
        ci21_7=0.0_wp
        ci21_8 =0.0_wp
        ci21_9 =0.0_wp
        ci21_10=0.0_wp
        ci21_11=0.0_wp

        do j=-m_2,m_2
           i_j=mod_my(i-j)

           i_j2=i_j+len_2

           j2=2*j
           j21=j2+1

           chj2=ch(j2)
           cgj2=cg(j2)

           chj21=ch(j21)
           cgj21=cg(j21)

           ci2_0  = ci2_0  + chj2*cd(i_j,it+0) + cgj2*cd(i_j2,it+0)
           ci2_1  = ci2_1  + chj2*cd(i_j,it+1) + cgj2*cd(i_j2,it+1)
           ci2_2  = ci2_2  + chj2*cd(i_j,it+2) + cgj2*cd(i_j2,it+2)
           ci2_3  = ci2_3  + chj2*cd(i_j,it+3) + cgj2*cd(i_j2,it+3)
           ci2_4  = ci2_4  + chj2*cd(i_j,it+4) + cgj2*cd(i_j2,it+4)
           ci2_5  = ci2_5  + chj2*cd(i_j,it+5) + cgj2*cd(i_j2,it+5)
           ci2_6  = ci2_6  + chj2*cd(i_j,it+6) + cgj2*cd(i_j2,it+6)
           ci2_7  = ci2_7  + chj2*cd(i_j,it+7) + cgj2*cd(i_j2,it+7)
           ci2_8   = ci2_8   + chj2*cd(i_j,it+8 ) + cgj2*cd(i_j2,it+8 )
           ci2_9   = ci2_9   + chj2*cd(i_j,it+9 ) + cgj2*cd(i_j2,it+9 )
           ci2_10  = ci2_10  + chj2*cd(i_j,it+10) + cgj2*cd(i_j2,it+10)
           ci2_11  = ci2_11  + chj2*cd(i_j,it+11) + cgj2*cd(i_j2,it+11)

           ci21_0 = ci21_0 + chj21*cd(i_j,it+0) + cgj21*cd(i_j2,it+0)
           ci21_1 = ci21_1 + chj21*cd(i_j,it+1) + cgj21*cd(i_j2,it+1)
           ci21_2 = ci21_2 + chj21*cd(i_j,it+2) + cgj21*cd(i_j2,it+2)
           ci21_3 = ci21_3 + chj21*cd(i_j,it+3) + cgj21*cd(i_j2,it+3)
           ci21_4 = ci21_4 + chj21*cd(i_j,it+4) + cgj21*cd(i_j2,it+4)
           ci21_5 = ci21_5 + chj21*cd(i_j,it+5) + cgj21*cd(i_j2,it+5)
           ci21_6 = ci21_6 + chj21*cd(i_j,it+6) + cgj21*cd(i_j2,it+6)
           ci21_7 = ci21_7 + chj21*cd(i_j,it+7) + cgj21*cd(i_j2,it+7)
           ci21_8  = ci21_8  + chj21*cd(i_j,it+8 ) + cgj21*cd(i_j2,it+8 )
           ci21_9  = ci21_9  + chj21*cd(i_j,it+9 ) + cgj21*cd(i_j2,it+9 )
           ci21_10 = ci21_10 + chj21*cd(i_j,it+10) + cgj21*cd(i_j2,it+10)
           ci21_11 = ci21_11 + chj21*cd(i_j,it+11) + cgj21*cd(i_j2,it+11)

        enddo

        i2=2*i
        i21=i2+1

        c1(it+0,i2 ) = ci2_0
        c1(it+1,i2 ) = ci2_1
        c1(it+2,i2 ) = ci2_2
        c1(it+3,i2 ) = ci2_3
        c1(it+4,i2 ) = ci2_4
        c1(it+5,i2 ) = ci2_5
        c1(it+6,i2 ) = ci2_6
        c1(it+7,i2 ) = ci2_7
        c1(it+8 ,i2 ) = ci2_8 
        c1(it+9 ,i2 ) = ci2_9 
        c1(it+10,i2 ) = ci2_10
        c1(it+11,i2 ) = ci2_11

        c1(it+0,i21) = ci21_0 
        c1(it+1,i21) = ci21_1 
        c1(it+2,i21) = ci21_2 
        c1(it+3,i21) = ci21_3 
        c1(it+4,i21) = ci21_4 
        c1(it+5,i21) = ci21_5 
        c1(it+6,i21) = ci21_6 
        c1(it+7,i21) = ci21_7 
        c1(it+8 ,i21) = ci21_8  
        c1(it+9 ,i21) = ci21_9  
        c1(it+10,i21) = ci21_10 
        c1(it+11,i21) = ci21_11 

     enddo
  enddo
!$omp enddo
  it0=it
!$ if (ithread.eq.0) then
!$  it0=12*int(nt/12)+1
  do it=it0,nt
     do i=0,len_2-1
        ci2 =0.0_wp
        ci21=0.0_wp
        do j=-m_2,m_2
           i_j=mod_my(i-j)
           ci2  = ci2  + ch(2*j  )*cd(i_j,it) + cg(2*j  )*cd(i_j+len_2,it)
           ci21 = ci21 + ch(2*j+1)*cd(i_j,it) + cg(2*j+1)*cd(i_j+len_2,it)
        enddo
        c1(it,2*i  ) = ci2
        c1(it,2*i+1) = ci21
     enddo
  enddo
!$ endif
!$omp end parallel
  !        call system_clock(ncount2,ncount_rate,ncount_max)
  !        tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !        write(95,'(a40,1x,e11.4,1x,f10.1,1x,i9)') 'syn_rot_per_old',tel,1.d-6*nflop/tel,nflop
  call f_free(mod_my)

END SUBROUTINE syn_rot_per_old

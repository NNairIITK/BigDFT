!> @file
!!  Optimized convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine ana_rot_shrink(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(-7:2*n+8,ndat), intent(in) :: x
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

   real(wp) :: ci1,ci2,ci3,ci4,ci5,ci6,ci7,ci8
   real(wp) :: di1,di2,di3,di4,di5,di6,di7,di8
!$omp parallel default (private) shared(ndat,n,ch,x,y,cg)
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
         k= l+2*i
         
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
        enddo
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
     enddo
  enddo
!$omp enddo


!$omp do  
  do j=(ndat/8)*8+1,ndat
     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
         k= l+2*i
            ci=ci+ch(l)*x(k    ,j)
            di=di+cg(l)*x(k    ,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo
  enddo
!$omp enddo
!$omp end parallel
END SUBROUTINE ana_rot_shrink


subroutine syn_rot_grow(n,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-7:2*n+8), intent(out) :: y
  !local variables
  integer :: i,j,k,l,lmax,lmin
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

    real(wp) :: so1,so2,so3,so4,so5,so6,so7,so8
    real(wp) :: se1,se2,se3,se4,se5,se6,se7,se8
   
!   real(gp)::tel
!   integer::ncount1,ncount2,ncount_rate,ncount_max
!   integer::mflop
   
!   ! (n+3): the range of i (average)
!   ! 8: filter length for l (average)
!   ! 4: number of flops in one line
!   ! 2: even and odd 
!   mflop=ndat*(n+3)*8*4*2
!   

!dee
!  open(unit=97,file='time_check',status='unknown')
!  call system_clock(ncount1,ncount_rate,ncount_max)

!$omp parallel default (private) shared(ndat,ch,x,cg,n,y)
!$omp do


  do j=0,ndat/8-1
     
     y(j*8+1,-7)=ch(-7)*x(0,j*8+1)+cg(-7)*x(n+1,j*8+1)
     y(j*8+2,-7)=ch(-7)*x(0,j*8+2)+cg(-7)*x(n+1,j*8+2)
     y(j*8+3,-7)=ch(-7)*x(0,j*8+3)+cg(-7)*x(n+1,j*8+3)
     y(j*8+4,-7)=ch(-7)*x(0,j*8+4)+cg(-7)*x(n+1,j*8+4)
     y(j*8+5,-7)=ch(-7)*x(0,j*8+5)+cg(-7)*x(n+1,j*8+5)
     y(j*8+6,-7)=ch(-7)*x(0,j*8+6)+cg(-7)*x(n+1,j*8+6)
     y(j*8+7,-7)=ch(-7)*x(0,j*8+7)+cg(-7)*x(n+1,j*8+7)
     y(j*8+8,-7)=ch(-7)*x(0,j*8+8)+cg(-7)*x(n+1,j*8+8)
    
     do i=-3,n+3
      if (i-4.ge.0) then
         k=i-4
         se1=ch(8)*x(k,j*8+1)+cg(8)*x(n+1+k,j*8+1)
         se2=ch(8)*x(k,j*8+2)+cg(8)*x(n+1+k,j*8+2)
         se3=ch(8)*x(k,j*8+3)+cg(8)*x(n+1+k,j*8+3)
         se4=ch(8)*x(k,j*8+4)+cg(8)*x(n+1+k,j*8+4)
         se5=ch(8)*x(k,j*8+5)+cg(8)*x(n+1+k,j*8+5)
         se6=ch(8)*x(k,j*8+6)+cg(8)*x(n+1+k,j*8+6)
         se7=ch(8)*x(k,j*8+7)+cg(8)*x(n+1+k,j*8+7)         
         se8=ch(8)*x(k,j*8+8)+cg(8)*x(n+1+k,j*8+8)         
         lmax=3 ! if i-4 >=0, then min(i,3)=3
      else
           se1=0.e0_wp
           se2=0.e0_wp
           se3=0.e0_wp
           se4=0.e0_wp
           se5=0.e0_wp
           se6=0.e0_wp
           se7=0.e0_wp
           se8=0.e0_wp
         lmax=i! if i-4<0 then min(i,3)=i 
      endif
      
      if (i+4.le.n) then
         ! l=-4 ;2*l+1=-7; k=i+4=<n
         k=i+4
            so1=ch(-7)*x(k,j*8+1)+cg(-7)*x(n+1+k,j*8+1)
            so2=ch(-7)*x(k,j*8+2)+cg(-7)*x(n+1+k,j*8+2)
            so3=ch(-7)*x(k,j*8+3)+cg(-7)*x(n+1+k,j*8+3)
            so4=ch(-7)*x(k,j*8+4)+cg(-7)*x(n+1+k,j*8+4)
            so5=ch(-7)*x(k,j*8+5)+cg(-7)*x(n+1+k,j*8+5)
            so6=ch(-7)*x(k,j*8+6)+cg(-7)*x(n+1+k,j*8+6)
            so7=ch(-7)*x(k,j*8+7)+cg(-7)*x(n+1+k,j*8+7)
            so8=ch(-7)*x(k,j*8+8)+cg(-7)*x(n+1+k,j*8+8)
         lmin=-3 ! if i+4=<n then max(i-n,-3)=-3
      else
           so1=0.e0_wp
           so2=0.e0_wp
           so3=0.e0_wp
           so4=0.e0_wp
           so5=0.e0_wp
           so6=0.e0_wp
           so7=0.e0_wp
           so8=0.e0_wp
         lmin=i-n! if i+4>n then max(i-n,-3)=i-n
      endif
      
       !do l=max(i-n,-3),min(i,3)
       do l=lmin,lmax
         k=i-l
         
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
        enddo

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
     enddo

     y(j*8+1,2*n+8)=ch(8)*x(n,j*8+1)+cg(8)*x(2*n+1,j*8+1)
     y(j*8+2,2*n+8)=ch(8)*x(n,j*8+2)+cg(8)*x(2*n+1,j*8+2)
     y(j*8+3,2*n+8)=ch(8)*x(n,j*8+3)+cg(8)*x(2*n+1,j*8+3)
     y(j*8+4,2*n+8)=ch(8)*x(n,j*8+4)+cg(8)*x(2*n+1,j*8+4)
     y(j*8+5,2*n+8)=ch(8)*x(n,j*8+5)+cg(8)*x(2*n+1,j*8+5)
     y(j*8+6,2*n+8)=ch(8)*x(n,j*8+6)+cg(8)*x(2*n+1,j*8+6)
     y(j*8+7,2*n+8)=ch(8)*x(n,j*8+7)+cg(8)*x(2*n+1,j*8+7)
     y(j*8+8,2*n+8)=ch(8)*x(n,j*8+8)+cg(8)*x(2*n+1,j*8+8)
    
  enddo
!$omp enddo

!$omp do

  do j=(ndat/8)*8+1,ndat

     y(j,-7)=ch(-7)*x(0,j)+cg(-7)*x(n+1,j)
    
     do i=-3,n+3

        se=0.e0_wp
        so=0.e0_wp
       do l=max(i-n,-4),min(i,4)
         k=i-l
            se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
            so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

     y(j,2*n+8)=ch(8)*x(n,j)+cg(8)*x(2*n+1,j)
  enddo
!$omp enddo
!$omp end parallel

!  call system_clock(ncount2,ncount_rate,ncount_max)
!  tel=dble(ncount2-ncount1)/dble(ncount_rate)
!  write(97,'(a40,1x,e10.3,1x,f6.1)') 'syn_rot_grow:',tel
!  close(97)
END SUBROUTINE syn_rot_grow



subroutine convrot_grow(n1,ndat,x,y)
  use module_base
  implicit none
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-lupfil:n1-lowfil), intent(out) :: y
  !local variables
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
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

   integer i,j,l,k
   real(wp) fill,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
!$omp parallel default (private) shared(fil,ndat,x,y,n1)
!$omp do


  do j=0,ndat/12-1
     
     do i=-lupfil,n1-lowfil
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
      
        do l=max(-i,lowfil),min(lupfil,n1-i)
         k=i+l   
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
        enddo
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

     enddo
  enddo
!$omp enddo

!$omp do
  do j=(ndat/12)*12+1,ndat
     do i=-lupfil,n1-lowfil
        tt=0.e0_wp
        do l=max(-i,lowfil),min(lupfil,n1-i)
         k=i+l   
            tt=tt+x(  k,j)*fil(l)
        enddo
      y(j,i)=tt

     enddo
  enddo
!$omp enddo
!$omp end parallel

END SUBROUTINE convrot_grow



subroutine convrot_shrink(n1,ndat,x,y)
  use module_base
  implicit none
  
  integer,parameter::lowfil=-7,lupfil=8
  integer, intent(in) :: n1,ndat
  real(wp), dimension(lowfil:n1+lupfil,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
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

   integer i,j,l,k
   real(wp) fill,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,tt
   

!$omp parallel default (private) shared(x,y,n1,ndat,fil)
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
         k=i+l   
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
        enddo
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

     enddo
  enddo

!$omp enddo

!$omp do

  do j=(ndat/12)*12+1,ndat
     do i=0,n1

        tt=0.e0_wp
        do l=lowfil,lupfil
         k=i+l   
            tt=tt+x(  k,j)*fil(l)
        enddo
      y(j,i)=tt

     enddo
  enddo
!$omp enddo
!$omp end parallel

  return
END SUBROUTINE convrot_shrink


subroutine convolut_kinetic_slab_c(n1,n2,n3,hgrid,x,y,c)
!   applies the kinetic energy operator onto x to get y. Works for surface BC
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
   integer mod_arr1(lowfil:n1+lupfil)   
   integer mod_arr3(lowfil:n3+lupfil)   

   call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
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

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  
  call conv_kin_x(x,y,(n2+1)*(n3+1))   
  call conv_kin_y
  call conv_kin_z(x,y,(n1+1)*(n2+1))
  
contains
   
   subroutine conv_kin_y
      implicit none
      real(wp) tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7
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
              
                  do l=max(lowfil,-i2),min(lupfil,n2-i2)
                  j=i2+l
                
                  tt0=tt0+x(i1,j,i3*8+0)*fil(l,2)
                  tt1=tt1+x(i1,j,i3*8+1)*fil(l,2)
                  tt2=tt2+x(i1,j,i3*8+2)*fil(l,2)
                  tt3=tt3+x(i1,j,i3*8+3)*fil(l,2)
                  tt4=tt4+x(i1,j,i3*8+4)*fil(l,2)
                  tt5=tt5+x(i1,j,i3*8+5)*fil(l,2)
                  tt6=tt6+x(i1,j,i3*8+6)*fil(l,2)
                  tt7=tt7+x(i1,j,i3*8+7)*fil(l,2)
               enddo
               y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0
               y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1
               y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2
               y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3
               y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4
               y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5
               y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6
               y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7
            enddo
         enddo
      enddo
!$omp enddo

!$omp do
        do i3=(n3/8)*8,n3
         do i1=0,n1
            do i2=0,n2
               tt=0.e0_wp
                  do l=max(lowfil,-i2),min(lupfil,n2-i2)
                  j=i2+l
                  tt=tt+x(i1,j   ,i3)*fil(l,2)
               enddo
               y(i1,i2,i3)=y(i1,i2,i3)+tt
            enddo
         enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_y

   
   subroutine conv_kin_x(x,y,ndat)
      implicit none
      integer,intent(in)::ndat
      real(wp),intent(in):: x(0:n1,ndat)
      real(wp),intent(out)::y(0:n1,ndat)
      real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
!$omp do
        do i=0,ndat/12-1
           do i1=0,n1
              tt1 =x(i1,i*12+1)*c
              tt2 =x(i1,i*12+2)*c
              tt3 =x(i1,i*12+3)*c
              tt4 =x(i1,i*12+4)*c
              tt5 =x(i1,i*12+5)*c
              tt6 =x(i1,i*12+6)*c
              tt7 =x(i1,i*12+7)*c
              tt8 =x(i1,i*12+8)*c
              tt9 =x(i1,i*12+9 )*c
              tt10=x(i1,i*12+10)*c
              tt11=x(i1,i*12+11)*c
              tt12=x(i1,i*12+12)*c
            
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
              enddo
              y(i1,i*12+1)=tt1
              y(i1,i*12+2)=tt2
              y(i1,i*12+3)=tt3
              y(i1,i*12+4)=tt4
              y(i1,i*12+5)=tt5
              y(i1,i*12+6)=tt6
              y(i1,i*12+7)=tt7
              y(i1,i*12+8)=tt8
              y(i1,i*12+9 )=tt9 
              y(i1,i*12+10)=tt10
              y(i1,i*12+11)=tt11
              y(i1,i*12+12)=tt12
           enddo
      enddo
!$omp enddo
!$omp do
        do i=(ndat/12)*12+1,ndat
           do i1=0,n1
              tt=x(i1,i)*c
              do l=lowfil,lupfil
                 j=mod_arr1(i1+l)
                 tt=tt+x(j   ,i)*fil(l,1)
              enddo
              y(i1,i)=tt
           enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_x
   
   subroutine conv_kin_z(x,y,ndat)
      implicit none
      integer,intent(in)::ndat
      real(wp),intent(in):: x(ndat,0:n1)
      real(wp),intent(inout)::y(ndat,0:n1)
      real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
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
              enddo
            
              y(i*12+1,i3)=y(i*12+1,i3)+tt1
              y(i*12+2,i3)=y(i*12+2,i3)+tt2
              y(i*12+3,i3)=y(i*12+3,i3)+tt3
              y(i*12+4,i3)=y(i*12+4,i3)+tt4
              y(i*12+5,i3)=y(i*12+5,i3)+tt5
              y(i*12+6,i3)=y(i*12+6,i3)+tt6
              y(i*12+7,i3)=y(i*12+7,i3)+tt7
              y(i*12+8,i3)=y(i*12+8,i3)+tt8
              y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 
              y(i*12+10,i3)=y(i*12+10,i3)+tt10
              y(i*12+11,i3)=y(i*12+11,i3)+tt11
              y(i*12+12,i3)=y(i*12+12,i3)+tt12
           enddo
      enddo
!$omp enddo
!$omp do
        do i=(ndat/12)*12+1,ndat
           do i3=0,n3
              tt=0.e0_wp
              do l=lowfil,lupfil
                 j=mod_arr3(i3+l)
                 tt=tt+x(i,j)*fil(l,3)
              enddo
              y(i,i3)=y(i,i3)+tt
           enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_z
  
END SUBROUTINE convolut_kinetic_slab_c



subroutine convolut_kinetic_slab_T(n1,n2,n3,hgrid,x,y,ekin)
!   applies the kinetic energy operator onto x to get y. Works for surface BC
!   y:=y-1/2Delta x
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp),intent(out)::ekin
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
   integer mod_arr1(lowfil:n1+lupfil)   
   integer mod_arr3(lowfil:n3+lupfil)   
      integer :: ncount0,ncount1,ncount_max,ncount_rate

   call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
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

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  ekin=0.0_wp

  !call system_clock(ncount0,ncount_rate,ncount_max)

  !call conv_kin_x(x,y,(n2+1)*(n3+1))   
  call conv_kin_x_new(n1,x,y,(n2+1)*(n3+1),lowfil,lupfil,fil,ekin,mod_arr1)

  !call conv_kin_y
  call conv_kin_y_new(n1,n2,n3,x,y,lowfil,lupfil,fil,ekin)

  !call conv_kin_z(x,y,(n1+1)*(n2+1))
  call conv_kin_z_new(n3,x,y,(n1+1)*(n2+1),lowfil,lupfil,fil,mod_arr3,ekin)
  
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:convolut_kinetic_slab_T',real(ncount1-ncount0)/real(ncount_rate)
contains
   
   subroutine conv_kin_x(x,y,ndat)
      implicit none
      integer,intent(in)::ndat
      real(wp),intent(in):: x(0:n1,ndat)
      real(wp),intent(inout)::y(0:n1,ndat)
      real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
!$omp do
  !omp parallel do default(private)&
  !omp shared(ndat,n1,lowfil,lupfil,mod_arr1,x,y,fil)&
  !omp reduction(+:ekin)
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
              enddo
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
           enddo
      enddo
  !omp end parallel do
!$omp enddo

!$omp do
        do i=(ndat/12)*12+1,ndat
           do i1=0,n1
              tt=0.e0_wp
              do l=lowfil,lupfil
                 j=mod_arr1(i1+l)
                 tt=tt+x(j   ,i)*fil(l,1)
              enddo
              y(i1,i)=y(i1,i)+tt ; ekin=ekin+tt*x(i1,i)
           enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_x
   
   subroutine conv_kin_y
      implicit none
      real(wp) tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7

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
              
                 do l=max(lowfil,-i2),min(lupfil,n2-i2)
                  j=i2+l
                
                  tt0=tt0+x(i1,j,i3*8+0)*fil(l,2)
                  tt1=tt1+x(i1,j,i3*8+1)*fil(l,2)
                  tt2=tt2+x(i1,j,i3*8+2)*fil(l,2)
                  tt3=tt3+x(i1,j,i3*8+3)*fil(l,2)
                  tt4=tt4+x(i1,j,i3*8+4)*fil(l,2)
                  tt5=tt5+x(i1,j,i3*8+5)*fil(l,2)
                  tt6=tt6+x(i1,j,i3*8+6)*fil(l,2)
                  tt7=tt7+x(i1,j,i3*8+7)*fil(l,2)
               enddo
               y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;    ekin=ekin+tt0*x(i1,i2,i3*8+0)
               y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;    ekin=ekin+tt1*x(i1,i2,i3*8+1)
               y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;    ekin=ekin+tt2*x(i1,i2,i3*8+2)
               y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;    ekin=ekin+tt3*x(i1,i2,i3*8+3)
               y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;    ekin=ekin+tt4*x(i1,i2,i3*8+4)
               y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;    ekin=ekin+tt5*x(i1,i2,i3*8+5)
               y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;    ekin=ekin+tt6*x(i1,i2,i3*8+6)
               y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;    ekin=ekin+tt7*x(i1,i2,i3*8+7)
            enddo                                 
         enddo
      enddo
!$omp enddo

!$omp do
        do i3=(n3/8)*8,n3
         do i1=0,n1
            do i2=0,n2
               tt=0.e0_wp
                 do l=max(lowfil,-i2),min(lupfil,n2-i2)
                  j=i2+l
                  tt=tt+x(i1,j   ,i3)*fil(l,2)
               enddo
               y(i1,i2,i3)=y(i1,i2,i3)+tt;   ekin=ekin+tt*x(i1,i2,i3)
            enddo
         enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_y

   subroutine conv_kin_z(x,y,ndat)
      implicit none
      integer,intent(in)::ndat
      real(wp),intent(in):: x(ndat,0:n3)
      real(wp),intent(inout)::y(ndat,0:n3)
      real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
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
              enddo
            
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
           enddo
      enddo
!$omp enddo
!$omp do
        do i=(ndat/12)*12+1,ndat
           do i3=0,n3
              tt=0.e0_wp
              do l=lowfil,lupfil
                 j=mod_arr3(i3+l)
                 tt=tt+x(i,j)*fil(l,3)
              enddo
              y(i,i3)=y(i,i3)+tt; ekin=ekin+tt*x(i,i3)
           enddo
      enddo
!$omp enddo
   END SUBROUTINE conv_kin_z
  
END SUBROUTINE convolut_kinetic_slab_T


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! New conv routine taken out from convolut_kinetic_slab_T. These routines 
! significantly speedup the version included in convolut_kinetic_slab_T.
subroutine conv_kin_x_new(n1,x,y,ndat,lowfil,lupfil,fil,ekin,mod_arr1)
  use module_base
  implicit none
  integer,intent(in)::n1,ndat,lowfil,lupfil
  real(wp), intent(in), dimension(lowfil:lupfil,3) :: fil
  real(wp),intent(in):: x(0:n1,ndat)
  real(wp),intent(inout)::y(0:n1,ndat),ekin
  integer,intent(in) :: mod_arr1(lowfil:n1+lupfil)   
  real(wp) tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,ekin_tmp
  integer :: i,i1,l,j,tt

  ekin_tmp=0.0_wp
  !$omp parallel do default(private)&
  !$omp shared(ndat,n1,lowfil,lupfil,mod_arr1,x,y,fil)&
  !$omp reduction(+:ekin_tmp)
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
      enddo
      y(i1,i*12+1)=y(i1,i*12+1)+tt1;    ekin_tmp=ekin_tmp+tt1*x(i1,i*12+1)
      y(i1,i*12+2)=y(i1,i*12+2)+tt2;    ekin_tmp=ekin_tmp+tt2*x(i1,i*12+2)
      y(i1,i*12+3)=y(i1,i*12+3)+tt3;    ekin_tmp=ekin_tmp+tt3*x(i1,i*12+3)
      y(i1,i*12+4)=y(i1,i*12+4)+tt4;    ekin_tmp=ekin_tmp+tt4*x(i1,i*12+4)
      y(i1,i*12+5)=y(i1,i*12+5)+tt5;    ekin_tmp=ekin_tmp+tt5*x(i1,i*12+5)
      y(i1,i*12+6)=y(i1,i*12+6)+tt6;    ekin_tmp=ekin_tmp+tt6*x(i1,i*12+6)
      y(i1,i*12+7)=y(i1,i*12+7)+tt7;    ekin_tmp=ekin_tmp+tt7*x(i1,i*12+7)
      y(i1,i*12+8)=y(i1,i*12+8)+tt8;    ekin_tmp=ekin_tmp+tt8*x(i1,i*12+8)
      y(i1,i*12+9 )=y(i1,i*12+9 )+tt9 ; ekin_tmp=ekin_tmp+tt9 *x(i1,i*12+9 )
      y(i1,i*12+10)=y(i1,i*12+10)+tt10; ekin_tmp=ekin_tmp+tt10*x(i1,i*12+10)
      y(i1,i*12+11)=y(i1,i*12+11)+tt11; ekin_tmp=ekin_tmp+tt11*x(i1,i*12+11)
      y(i1,i*12+12)=y(i1,i*12+12)+tt12; ekin_tmp=ekin_tmp+tt12*x(i1,i*12+12)
    enddo
  enddo
  !$omp end parallel do

  ekin=ekin+ekin_tmp

  ekin_tmp=0.0_wp
  !$omp parallel do default(private) &
  !$omp shared(ndat,n1,lowfil,lupfil,mod_arr1,x,y,fil) &
  !$omp reduction(+:ekin_tmp)
  do i=(ndat/12)*12+1,ndat
    do i1=0,n1
      tt=0.e0_wp
      do l=lowfil,lupfil
        j=mod_arr1(i1+l)
        tt=tt+x(j   ,i)*fil(l,1)
      enddo
      y(i1,i)=y(i1,i)+tt ; ekin_tmp=ekin_tmp+tt*x(i1,i)
    enddo
  enddo
  !$omp end parallel do

  ekin=ekin+ekin_tmp
END SUBROUTINE conv_kin_x_new

subroutine conv_kin_y_new(n1,n2,n3,x,y,lowfil,lupfil,fil,ekin)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,lowfil,lupfil
  real(wp), intent(in), dimension(lowfil:lupfil,3) :: fil
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  real(wp), intent(inout) :: ekin
  integer :: i1,i2,i3,j,l
  real(wp) :: tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7,ekin_tmp,tt

  ekin_tmp=0.0_wp
  !$omp parallel default(private) &
  !$omp shared(n1,n2,n3,lowfil,lupfil,fil,x,y) &
  !$omp reduction(+:ekin_tmp) 
  !$omp do schedule(static,1)
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
           do l=max(lowfil,-i2),min(lupfil,n2-i2)
            j=i2+l
            tt0=tt0+x(i1,j,i3*8+0)*fil(l,2)
            tt1=tt1+x(i1,j,i3*8+1)*fil(l,2)
            tt2=tt2+x(i1,j,i3*8+2)*fil(l,2)
            tt3=tt3+x(i1,j,i3*8+3)*fil(l,2)
            tt4=tt4+x(i1,j,i3*8+4)*fil(l,2)
            tt5=tt5+x(i1,j,i3*8+5)*fil(l,2)
            tt6=tt6+x(i1,j,i3*8+6)*fil(l,2)
            tt7=tt7+x(i1,j,i3*8+7)*fil(l,2)
         enddo
         y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;    ekin_tmp=ekin_tmp+tt0*x(i1,i2,i3*8+0)
         y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;    ekin_tmp=ekin_tmp+tt1*x(i1,i2,i3*8+1)
         y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;    ekin_tmp=ekin_tmp+tt2*x(i1,i2,i3*8+2)
         y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;    ekin_tmp=ekin_tmp+tt3*x(i1,i2,i3*8+3)
         y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;    ekin_tmp=ekin_tmp+tt4*x(i1,i2,i3*8+4)
         y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;    ekin_tmp=ekin_tmp+tt5*x(i1,i2,i3*8+5)
         y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;    ekin_tmp=ekin_tmp+tt6*x(i1,i2,i3*8+6)
         y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;    ekin_tmp=ekin_tmp+tt7*x(i1,i2,i3*8+7)
      enddo                                 
    enddo
  enddo
  !$omp enddo
  !$omp end parallel

  ekin=ekin+ekin_tmp

  ekin_tmp=0.0_wp

  !$omp parallel default(private) &
  !$omp shared(n1,n2,n3,lowfil,lupfil,fil,x,y)&
  !$omp reduction(+:ekin_tmp)
  !$omp do schedule(static,1)
  do i3=(n3/8)*8,n3
    do i1=0,n1
       do i2=0,n2
          tt=0.e0_wp
            do l=max(lowfil,-i2),min(lupfil,n2-i2)
             j=i2+l
             tt=tt+x(i1,j,i3)*fil(l,2)
            enddo
          y(i1,i2,i3)=y(i1,i2,i3)+tt
          ekin_tmp=ekin_tmp+tt*x(i1,i2,i3)
        enddo
     enddo
  enddo
  !$omp enddo
  !$omp end parallel
  ekin=ekin+ekin_tmp

END SUBROUTINE conv_kin_y_new


subroutine conv_kin_z_new(n3,x,y,ndat,lowfil,lupfil,fil,mod_arr3,ekin)
  use module_base
  implicit none
  integer,intent(in)::ndat,lowfil,lupfil,n3
  real(wp), intent(in), dimension(lowfil:lupfil,3) :: fil
  real(wp),intent(in):: x(ndat,0:n3)
  integer,intent(in) :: mod_arr3(lowfil:n3+lupfil)
  real(wp),intent(inout)::y(ndat,0:n3),ekin
  real(wp) ::  tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,ekin_tmp,tt
  integer :: i,i3,l,j

  ekin_tmp=0.0_wp

  !$omp parallel default(private) &
  !$omp shared(ndat,n3,lowfil,lupfil,fil,mod_arr3,x,y)&
  !$omp reduction(+:ekin_tmp)
  !$omp do schedule(static,1)
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
       enddo
       y(i*12+1,i3)=y(i*12+1,i3)+tt1;       ekin_tmp=ekin_tmp+tt1*x(i*12+1,i3)
       y(i*12+2,i3)=y(i*12+2,i3)+tt2;       ekin_tmp=ekin_tmp+tt2*x(i*12+2,i3)
       y(i*12+3,i3)=y(i*12+3,i3)+tt3;       ekin_tmp=ekin_tmp+tt3*x(i*12+3,i3)
       y(i*12+4,i3)=y(i*12+4,i3)+tt4;       ekin_tmp=ekin_tmp+tt4*x(i*12+4,i3)
       y(i*12+5,i3)=y(i*12+5,i3)+tt5;       ekin_tmp=ekin_tmp+tt5*x(i*12+5,i3)
       y(i*12+6,i3)=y(i*12+6,i3)+tt6;       ekin_tmp=ekin_tmp+tt6*x(i*12+6,i3)
       y(i*12+7,i3)=y(i*12+7,i3)+tt7;       ekin_tmp=ekin_tmp+tt7*x(i*12+7,i3)
       y(i*12+8,i3)=y(i*12+8,i3)+tt8;       ekin_tmp=ekin_tmp+tt8*x(i*12+8,i3)
       y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 ;    ekin_tmp=ekin_tmp+tt9*x(i*12+9 ,i3)
       y(i*12+10,i3)=y(i*12+10,i3)+tt10;    ekin_tmp=ekin_tmp+tt10*x(i*12+10,i3)
       y(i*12+11,i3)=y(i*12+11,i3)+tt11;    ekin_tmp=ekin_tmp+tt11*x(i*12+11,i3)
       y(i*12+12,i3)=y(i*12+12,i3)+tt12;    ekin_tmp=ekin_tmp+tt12*x(i*12+12,i3)
    enddo
  enddo
  !$omp enddo
  !$omp end parallel

  ekin=ekin+ekin_tmp
  ekin_tmp=0.0_wp
  
  !$omp parallel default(private) &
  !$omp shared(ndat,n3,lowfil,lupfil,fil,x,y,mod_arr3)&
  !$omp reduction(+:ekin_tmp)
  !$omp do schedule(static,1)
  do i=(ndat/12)*12+1,ndat
    do i3=0,n3
       tt=0.e0_wp
       do l=lowfil,lupfil
          j=mod_arr3(i3+l)
          tt=tt+x(i,j)*fil(l,3)
       enddo
       y(i,i3)=y(i,i3)+tt; 
       ekin_tmp=ekin_tmp+tt*x(i,i3)
    enddo
  enddo
  !$omp enddo
  !$omp end parallel
  ekin=ekin+ekin_tmp

END SUBROUTINE conv_kin_z_new

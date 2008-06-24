subroutine ana_rot_per_old(right,nt,c,cd_1)
  !
  !      forward wavelet transform, analysis, periodic
  !
  use module_base
  implicit none
  integer, intent(in) :: right,nt
  real(wp), dimension(0:right,nt), intent(in) :: c
  real(wp), dimension(nt,0:right), intent(out) :: cd_1
  !local variables
  integer, parameter :: m=8
  integer :: lenc,len_2,it,i,i2,ji2
  real(wp) :: ci,di
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

  lenc=right+1
  len_2=lenc/2

  do it=1,nt
     !      *nt       
     do i=0,len_2-1
        !        *len_2
        i2=2*i
        ci=0.0_wp
        di=0.0_wp
        do j=1-m,m
           !          *2*m (i.e.,16)
           ji2=modulo(j+i2,lenc)
           ci=ci+ch(j)*c(ji2,it)
           di=di+cg(j)*c(ji2,it)
           !            *4: do not count modulo             
        enddo
        cd_1(it,i)=ci
        cd_1(it,len_2+i)=di
     enddo
  enddo
  !      ana_rot_per_old: nt*len_2*2*m*4 flops

end subroutine ana_rot_per_old

subroutine syn_rot_per_old(right1,nt,cd,c1)
  !
  !     backward wavelet transform, synthesis, periodic
  !
  use module_base
  implicit none
  integer, intent(in) :: right1,nt
  real(wp), dimension(0:right1,nt), intent(in) :: cd
  real(wp), dimension(nt,0:right1), intent(out) :: c1
  !local variables
  integer, parameter :: m=8
  integer :: m_2,len_2,i,j,it
  real(wp) :: ci2,ci21
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

  m_2=m/2
  len_2=(right1+1)/2

  do it=1,nt
     !       *nt
     do i=0,len_2-1
        !         *len_2
        ci2 =0.0_wp
        ci21=0.0_wp
        do j=-m_2,m_2
           !           *(2*m_2+1)
           i_j=modulo(i-j,len_2)
           ci2  = ci2  + ch(2*j  )*cd(i_j,it) + cg(2*j  )*cd(i_j+len_2,it)
           ci21 = ci21 + ch(2*j+1)*cd(i_j,it) + cg(2*j+1)*cd(i_j+len_2,it)
           !             *8: do not count modulo
        enddo
        c1(it,2*i  ) = ci2
        c1(it,2*i+1) = ci21
     enddo
  enddo
  !       syn_rot_per_old:  nt*len_2*(2*m_2+1)*8 flops

end subroutine syn_rot_per_old


subroutine analyse_shrink(n1,n2,n3,ww,y,x)
  ! A analysis wavelet transformation where the size of the data is forced to shrink
  ! The input array y is overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: x
  !local variables
  integer :: nt

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+16)*(2*n3+16)
  call  ana_rot_shrink(n1,nt,y,ww)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+16)*(2*n1+2)
  call  ana_rot_shrink(n2,nt,ww,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink(n3,nt,y,x)

end subroutine analyse_shrink

subroutine synthese_grow(n1,n2,n3,ww,x,y)
  ! A synthesis wavelet transformation where the size of the data is allowed to grow
  ! The input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: x
  real(wp), dimension(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8), intent(inout) :: ww
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow(n1,nt,x,y)
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+16)
  call  syn_rot_grow(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+16)*(2*n2+16)
  call  syn_rot_grow(n3,nt,ww,y)

end subroutine synthese_grow

subroutine analyse_shrink_per(n1,n2,n3,ww,y,x)
  ! A analysis wavelet transformation where the size of the data is forced to shrink
  ! The input array y is NOT overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:2*n2+1,0:2*n3+1,0:2*n1+1), intent(inout) :: ww
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(inout) :: y
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: x
  !local variables
  integer :: nt

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  ana_rot_shrink_per(n1,nt,y,x)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_shrink_per(n2,nt,x,ww)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_shrink_per(n3,nt,ww,x)

end subroutine analyse_shrink_per

subroutine synthese_grow_per(n1,n2,n3,ww,x,y)
  ! A synthesis wavelet transformation where the size of the data is allowed to grow
  ! The input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: x
  real(wp), dimension(0:2*n2+1,0:2*n3+1,0:2*n1+1), intent(inout) :: ww
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_grow_per(n1,nt,x,y) 
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+2)
  call  syn_rot_grow_per(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call  syn_rot_grow_per(n3,nt,ww,y)

end subroutine synthese_grow_per


subroutine ana_rot_shrink(n,ndat,x,y)
  use module_base
  integer, intent(in) :: n,ndat
  real(wp), dimension(-7:2*n+8,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j
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

  do j=1,ndat
     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           ci=ci+ch(l)*x(l+2*i,j)
           di=di+cg(l)*x(l+2*i,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo
  enddo

end subroutine ana_rot_shrink

subroutine ana_rot_shrink_per(n,ndat,x,y)
  use module_base
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j
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

  do j=1,ndat

     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           k=modulo(l+2*i,2*n+2)
!           ci=ci+ch(l)*x(l+2*i,j)
!           di=di+cg(l)*x(l+2*i,j)
            ci=ci+ch(l)*x(k    ,j)
            di=di+cg(l)*x(k    ,j)
        enddo
        y(j,i)=ci
        y(j,n+1+i)=di
     enddo

  enddo
end subroutine ana_rot_shrink_per

subroutine syn_rot_grow(n,ndat,x,y)
  use module_base
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-7:2*n+8), intent(out) :: y
  !local variables
  integer :: i,j
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

     i=-4
     so=0.e0_wp
     do l=max(i-n,-4),min(i,4)
        so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
     enddo
     y(j,2*i+1)=so

     do i=-3,n+3
        se=0.e0_wp
        so=0.e0_wp
        do l=max(i-n,-4),min(i,4)
           se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
           so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

     i=n+4
     se=0.e0_wp
     do l=max(i-n,-4),min(i,4)
        se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
     enddo
     y(j,2*i  )=se

  enddo

end subroutine syn_rot_grow


subroutine syn_rot_grow_per(n,ndat,x,y)
  use module_base
  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j
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
!       do l=max(i-n,-4),min(i,4)
        do l=-4,4
           k=modulo(i-l,n+1)
           !se=se+ch(2*l  )*x(i-l,j)+cg(2*l  )*x(n+1+i-l,j)
           !so=so+ch(2*l+1)*x(i-l,j)+cg(2*l+1)*x(n+1+i-l,j)
           se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
           so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
        enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
     enddo

  enddo

end subroutine syn_rot_grow_per

! Applies the magic filter matrix ( no transposition) ; data set grows
! The input array x is not overwritten
! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_per(n1,n2,n3,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_n_per'
  integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat,i_stat,i_all
  real(wp), dimension(:,:,:), allocatable :: ww
  real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  data filCUDA / &
       8.4334247333529341094733325815816e-7_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       0.612625895831207982195380597e-1_4,&
       0.9940415697834003993178616713_4,&
       -0.604895289196983516002834636e-1_4, &
       -0.2103025160930381434955489412839065067e-1_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       2.72734492911979659657715313017228e-6_4 /


  if (.not. GPUcomputing) then !traditional CPU computation
     allocate(ww(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     !  (i1,i2*i3) -> (i2*i3,I1)
     ndat=(n2+1)*(n3+1)
     call convrot_grow_per(n1,ndat,x,y)
     !  (i2,i3*I1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_grow_per(n2,ndat,y,ww)
     !  (i3,I1*I2) -> (iI*I2,I3)
     ndat=(n1+1)*(n2+1)
     call convrot_grow_per(n3,ndat,ww,y)

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)
  else
     allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wx,'wx',subname)

     allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wy,'wy',subname)

     !copy input data
     wx=real(x,kind=4)

     !call cuda C interface
     call intertamponcGPU(n1,n2,n3,wx,wy,filCUDA,lowfil,lupfil)

     !restore output data
     y=real(wy,kind=wp)

     i_all=-product(shape(wx))*kind(wx)
     deallocate(wx,stat=i_stat)
     call memocc(i_stat,i_all,'wx',subname)

     i_all=-product(shape(wy))*kind(wy)
     deallocate(wy,stat=i_stat)
     call memocc(i_stat,i_all,'wy',subname)
  end if
end subroutine convolut_magic_n_per

! Applies the magic filter matrix transposed ; data set shrinks
! The input array x is overwritten
! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_per(n1,n2,n3,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_t_per'
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: ndat,i_stat,i_all
  real(wp), dimension(:,:,:), allocatable :: ww
  real(kind=4), dimension(:,:,:), allocatable :: wx,wy !these are used for copy in GPU case
  real(kind=4) filCUDA(lowfil:lupfil) !array of filters to be passed to CUDA interface
  data filCUDA / &
       2.72734492911979659657715313017228e-6_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.2103025160930381434955489412839065067e-1_4,&
       -0.604895289196983516002834636e-1_4,&
       0.9940415697834003993178616713_4,&
       0.612625895831207982195380597e-1_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       8.4334247333529341094733325815816e-7_4 /

  
  if (.not. GPUcomputing) then
     allocate(ww(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     !  (I1,I2*I3) -> (I2*I3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_shrink_per(n1,ndat,x,ww)
     !  (I2,I3*i1) -> (I3*i1,i2)
     ndat=(n3+1)*(n1+1)
     call convrot_shrink_per(n2,ndat,ww,x)
     !  (I3,i1*i2) -> (i1*i2,i3)
     ndat=(n1+1)*(n2+1)
     call convrot_shrink_per(n3,ndat,x,y)

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)
  else
     allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wx,'wx',subname)

     allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wy,'wy',subname)

     !copy input data
     wx=real(x,kind=4)

     !call cuda C interface
     call intertamponcGPU(n1,n2,n3,wx,wy,filCUDA,lowfil,lupfil)

     !restore output data
     y=real(wy,kind=wp)

     i_all=-product(shape(wx))*kind(wx)
     deallocate(wx,stat=i_stat)
     call memocc(i_stat,i_all,'wx',subname)

     i_all=-product(shape(wy))*kind(wy)
     deallocate(wy,stat=i_stat)
     call memocc(i_stat,i_all,'wy',subname)
  end if

end subroutine convolut_magic_t_per

! Simple non-optimized version of the major convolution routines
subroutine convrot_grow(n1,ndat,x,y)
  use module_base
  implicit none
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-lupfil:n1-lowfil), intent(out) :: y
  !local variables
  integer :: i,j,l
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

  do j=1,ndat
     do i=-lupfil,n1-lowfil
        tt=0.e0_wp
        do l=max(-i,lowfil),min(lupfil,n1-i)
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

end subroutine convrot_grow

! Simple non-optimized version of the major convolution routines

subroutine convrot_grow_per(n1,ndat,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
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

  do j=1,ndat
!     do i=-lupfil,n1-lowfil
     do i=0,n1
        tt=0.e0_wp
!        do l=max(-i,lowfil),min(lupfil,n1-i)
        do l=lowfil,lupfil
           k=modulo(i+l,n1+1)	
           !tt=tt+x(i+l,j)*fil(l)
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

end subroutine convrot_grow_per

subroutine convrot_shrink(n1,ndat,x,y)
  use module_base
  implicit none
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(lowfil:n1+lupfil,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /

  do j=1,ndat
     do i=0,n1
        tt=0.e0_wp
        do l=lowfil,lupfil
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

end subroutine convrot_shrink

subroutine convrot_shrink_per(n1,ndat,x,y)
  use module_base
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: i,j,k,l
  real(wp) :: tt
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /


  do j=1,ndat
     do i=0,n1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1+1)	
           !tt=tt+x(i+l,j)*fil(l)
           tt=tt+x(k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

end subroutine convrot_shrink_per


subroutine convolut_kinetic_per(n1,n2,n3,hgrid,x,y)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
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
  
  do i3=0,n3
     ! (1/2) d^2/dx^2
     do i2=0,n2
        do i1=0,n1
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i1+l,n1+1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1
        do i2=0,n2
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2+1)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=0,n3
           tt=0.0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3+1)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
  enddo
  
end subroutine convolut_kinetic_per

subroutine comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w2,w1,xc,xf,y,ibyz_c,ibzxx_c,ibxxyy_c,&
     ibyz_f,ibzxx_f,ibxxyy_f,ibyyzz_r)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx_c
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) :: ibxxyy_c
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz_f
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx_f
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy_f
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in):: ibyyzz_r
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: xf
  real(wp), dimension(max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       (2*n1+31)*(n2+1)*(n3+1))), intent(inout) :: w1 !work
  real(wp), dimension(max((n3+1)*(2*n1+31)*(2*n2+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))), intent(inout) :: w2 ! work
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(out) :: y

  call comb_grow_c(n1,n2,n3,w1,w2,xc,y,ibyz_c,ibzxx_c,ibxxyy_c,ibyyzz_r)

  call comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       w1,w2,xf,y,ibyz_f,ibzxx_f,ibxxyy_f)                

end subroutine comb_grow_all


subroutine comb_grow_c(n1,n2,n3,w1,w2,x,y,ibyz,ibzxx,ibxxyy,ibyyzz_r)
  ! In 3d,            
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  !  the size of the data is allowed to grow
  ! The input array x is not overwritten
  ! However, the output array y contains nonphysical values
  ! outside of the localization region
  ! that remain from the first comb_grow
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) ::  ibxxyy
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in):: ibyyzz_r
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n2,0:n3,-14:2*n1+16), intent(inout) :: w1
  real(wp), dimension(0:n3,-14:2*n1+16,-14:2*n2+16), intent(inout) :: w2
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(out) :: y

  ! i1,i2,i3 -> i2,i3,I1
  call  comb_rot_grow_loc_square_1(n1,n2,n3,x,w1,ibyz,ibzxx,.true.) 

  ! i2,i3,I1 -> i3,I1,I2
  call  comb_rot_grow_loc_square_1(n2,n3,2*n1+30,w1,w2,ibzxx,ibxxyy,.true.) 

  ! i3,I1,I2  -> I1,I2,I3
  call  comb_rot_grow_loc_square_1(n3,2*n1+30,2*n2+30,w2,y,ibxxyy,ibyyzz_r,.false.) 

END SUBROUTINE comb_grow_c

subroutine comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w1,w2,x,y,ibyz,ibzxx,ibxxyy)
  ! In 3d,            
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  !  the size of the data is allowed to grow
  ! The input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
  real(wp), dimension(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: w1
  real(wp), dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(in) :: w2
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: y
  !local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,w1,ibyz,ibzxx) 

  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*m1+31)
  call comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,ibzxx,ibxxyy) 

  ! i3,I1,I2  -> I1,I2,I3: add the result to y
  call comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w2,y,ibxxyy)

end subroutine comb_grow_tree

subroutine comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
  ! In one dimesnion,    
  ! with optimised cycles
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  !  the size of the data is allowed to grow
  use module_base
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ib
  integer, dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: ib2 
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
  real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(out) :: y
  !local variables
  !integer ncount0,ncount1,ncount_rate,ncount_max,nflop
  !real(kind=8) tel
  integer :: l2,l3,i,t,l1
  real(wp) y2i__11,y2i__21,y2i1_11,y2i1_21
  real(wp) y2i__12,y2i__22,y2i1_12,y2i1_22

  include 'v_17.inc'


  !    open(unit=20,file='tree.flop')

  !    nflop=0
  !    do l2=nfl2,nfu2
  !        do l3=nfl3,nfu3
  !            if (ib(2,l2,l3).ge.ib(1,l2,l3)) nflop=nflop+(ib(2,l2,l3)-ib(1,l2,l3)+1)*31*2*7
  !        enddo
  !    enddo

  !    call system_clock(ncount0,ncount_rate,ncount_max)

  !   y=0.d0
  do l1=-14+2*nfl1,2*nfu1+16
     do l3=nfl3,nfu3
        y(:,:,ib2(1,l3,l1):ib2(2,l3,l1),l3,l1)=0.d0
     enddo
  enddo

  do l3=nfl3,nfu3
     do l2=nfl2,nfu2

        if (ib(1,l2,l3).le.ib(2,l2,l3)) then

           y(1,1,l2,l3,2*ib(2,l2,l3)+16)=    fil2(16,2)*x(1,ib(2,l2,l3),l2,l3)
           y(2,1,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(2,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(3,ib(2,l2,l3),l2,l3)

           y(1,2,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(4,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(5,ib(2,l2,l3),l2,l3)
           y(2,2,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(6,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(7,ib(2,l2,l3),l2,l3)

           do i=ib(1,l2,l3)-7,ib(2,l2,l3)+7 
              y2i__11=0.d0
              y2i__21=0.d0
              y2i__12=0.d0
              y2i__22=0.d0

              y2i1_11=0.d0
              y2i1_21=0.d0
              y2i1_12=0.d0
              y2i1_22=0.d0
              do t=max(i-8,ib(1,l2,l3)),min(i+7,ib(2,l2,l3))
                 y2i__11=y2i__11                               +fil2(2*(i-t)  ,2)*x(1,t,l2,l3)
                 y2i__21=y2i__21+fil2(2*(i-t)  ,1)*x(2,t,l2,l3)+fil2(2*(i-t)  ,2)*x(3,t,l2,l3)
                 y2i__12=y2i__12+fil2(2*(i-t)  ,1)*x(4,t,l2,l3)+fil2(2*(i-t)  ,2)*x(5,t,l2,l3)
                 y2i__22=y2i__22+fil2(2*(i-t)  ,1)*x(6,t,l2,l3)+fil2(2*(i-t)  ,2)*x(7,t,l2,l3)

                 y2i1_11=y2i1_11                               +fil2(2*(i-t)+1,2)*x(1,t,l2,l3)
                 y2i1_21=y2i1_21+fil2(2*(i-t)+1,1)*x(2,t,l2,l3)+fil2(2*(i-t)+1,2)*x(3,t,l2,l3)
                 y2i1_12=y2i1_12+fil2(2*(i-t)+1,1)*x(4,t,l2,l3)+fil2(2*(i-t)+1,2)*x(5,t,l2,l3)
                 y2i1_22=y2i1_22+fil2(2*(i-t)+1,1)*x(6,t,l2,l3)+fil2(2*(i-t)+1,2)*x(7,t,l2,l3)
              enddo
              y(1,1,l2,l3,2*i  )=y2i__11
              y(2,1,l2,l3,2*i  )=y2i__21
              y(1,2,l2,l3,2*i  )=y2i__12
              y(2,2,l2,l3,2*i  )=y2i__22

              y(1,1,l2,l3,2*i+1)=y2i1_11
              y(2,1,l2,l3,2*i+1)=y2i1_21
              y(1,2,l2,l3,2*i+1)=y2i1_12
              y(2,2,l2,l3,2*i+1)=y2i1_22
           enddo
        endif

     enddo
  enddo

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !
  !    write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_1


subroutine comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
  ! In one dimesnion,    
  ! with optimised cycles
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  !  the size of the data is allowed to grow
  use module_base
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: ib 
  integer, dimension(2,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(in) :: ib2
  real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: x
  real(wp), dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(out) :: y
  !integer ncount0,ncount1,ncount_rate,ncount_max,nflop
  !real(kind=8) :: tel
  integer l1,l3,i,t,l2
  real(wp) y2i__1,y2i__2,y2i1_1,y2i1_2

  include 'v_17.inc'

  !    open(unit=20,file='tree.flop')
  !    nflop=0
  !    do l3=nfl3,nfu3
  !        do l1=-14+2*nfl1,2*nfu1+16
  !            if (ib(2,l3,l1).ge.ib(1,l3,l1)) nflop=nflop+(ib(2,l3,l1)-ib(1,l3,l1)+1)
  !        enddo
  !    enddo
  !    nflop=nflop*31*2*4

  !    call system_clock(ncount0,ncount_rate,ncount_max)

  !     y=0.d0
  do l2=-14+2*nfl2,2*nfu2+16
     do l1=-14+2*nfl1,2*nfu1+16
        y(:,ib2(1,l1,l2):ib2(2,l1,l2),l1,l2)=0.d0
     enddo
  enddo

  do l1=-14+2*nfl1,2*nfu1+16
     do l3=nfl3,nfu3

        if (ib(1,l3,l1).le.ib(2,l3,l1)) then
           y(1,l3,l1,2*ib(2,l3,l1)+16)=&
                fil2(16,1)*x(1,1,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,1,ib(2,l3,l1),l3,l1)
           y(2,l3,l1,2*ib(2,l3,l1)+16)=&
                fil2(16,1)*x(1,2,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,2,ib(2,l3,l1),l3,l1)

           do i=ib(1,l3,l1)-7,ib(2,l3,l1)+7 
              y2i__1=0.d0
              y2i__2=0.d0
              y2i1_1=0.d0
              y2i1_2=0.d0
              do t=max(i-8,ib(1,l3,l1)),min(i+7,ib(2,l3,l1))
                 y2i__1=y2i__1+fil2(2*(i-t)  ,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,1,t,l3,l1)
                 y2i__2=y2i__2+fil2(2*(i-t)  ,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,2,t,l3,l1)
                 y2i1_1=y2i1_1+fil2(2*(i-t)+1,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,1,t,l3,l1)
                 y2i1_2=y2i1_2+fil2(2*(i-t)+1,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,2,t,l3,l1)
              enddo
              y(1,l3,l1,2*i  )=y2i__1
              y(2,l3,l1,2*i  )=y2i__2
              y(1,l3,l1,2*i+1)=y2i1_1
              y(2,l3,l1,2*i+1)=y2i1_2
           enddo
        endif

     enddo
  enddo

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !    write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_2

subroutine comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,&
     ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibzzx_f,ibyyzz_f,xc,xf)
  ! In 3d,            
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  ! The input array y is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c
  integer, dimension(2,-14:2*n3+16,0:n1), intent(in) :: ibzzx_c
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_c
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy_f
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx_f
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz_f
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: y
  real(wp), dimension(max(2*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1),&
       (2*n2+31)*(2*n3+31)*(n1+1))), intent(inout) :: w1
  real(wp), dimension(max(4*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1)*(nfu2-nfl2+1),&
       (2*n3+31)*(n1+1)*(n2+1))), intent(inout) :: w2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xf

  !    perform the combined transform    
  call comb_shrink_loc_c(0,n1,0,n2,0,n3,w1,w2,y,xc,1,1,1,&
       ibxy_c,ibzzx_c,ibyyzz_c) ! for scfunctions
  !    for wavelets:

  call comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,xf,&
       ibxy_f,ibzzx_f,ibyyzz_f)

end subroutine comb_shrink


subroutine comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,&
     ibxy,ibzzx,ibyyzz)
  ! In 3d,            
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The output is only the l1,l2,l3 wavelet component
  ! The size of the data is forced to shrink
  ! The input array y is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: y ! input
  real(wp), dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(inout) :: w1
  real(wp), dimension(4,-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2), intent(inout) :: w2
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x
  !local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! I1,I2,I3 -> I2,I3,i1
  call comb_rot_shrink_loc_1(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,y,w1,ibyyzz)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*m3+31)*(m1+1)
  call comb_rot_shrink_loc_2(nt,w1,w2,nfl2,nfu2,ibzzx)

  ! I3,i1,i2 -> i1,i2,i3
  nt=(m1+1)*(m2+1)
  call comb_rot_shrink_loc_3(nt,w2,x,nfl3,nfu3,ibxy)

end subroutine comb_shrink_loc_f

subroutine comb_shrink_loc_c(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,l1,l2,l3,&
     ibxy,ibzzx,ibyyzz)
  ! In 3d,            
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The output is only the l1,l2,l3 wavelet component
  ! The size of the data is forced to shrink
  ! The input array y is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,l1,l2,l3
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy  
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx 
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz
  real(wp), dimension(-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16),&
       intent(in) :: y!input
  real(wp), dimension(-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(inout) :: w1
  real(wp), dimension(-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2), intent(inout) :: w2
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x!output
  !local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*m2+31)*(2*m3+31)
  call comb_rot_shrink_loc(nt,y,w1,l1,nfl1,nfu1,ibyyzz)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*m3+31)*(m1+1)
  call comb_rot_shrink_loc(nt,w1,w2,l2,nfl2,nfu2,ibzzx)

  ! I3,i1,i2 -> i1,i2,i3
  nt=(m1+1)*(m2+1)
  call comb_rot_shrink_loc(nt,w2,x,l3,nfl3,nfu3,ibxy)

end subroutine comb_shrink_loc_c

subroutine comb_rot_shrink_loc_3(ndat,x,y,nfl,nfu,ib)
  ! In one dimension,    
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  use module_base
  implicit none
  integer, parameter :: lowfil2=-14,lupfil2=16
  integer, intent(in) :: ndat,nfl,nfu
  integer, dimension(2,ndat), intent(in) :: ib
  real(wp), dimension(2,2,lowfil2+2*nfl:2*nfu+lupfil2,ndat), intent(in) :: x
  real(wp), dimension(7,ndat,nfl:nfu), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: ci112,ci121,ci122,ci211,ci212,ci221,ci222
  include 'v.inc'

  !nflop=0
  !open(unit=20,file='long.flop')
  !call system_clock(ncount0,ncount_rate,ncount_max)

  do j=1,ndat
     do i=ib(1,j),ib(2,j)
        ci112=0._wp
        ci121=0._wp
        ci122=0._wp
        ci211=0._wp
        ci212=0._wp
        ci221=0._wp
        ci222=0._wp

        !nflop=nflop+(lupfil2-lowfil2+1)*2*7
        do l=lowfil2+2*i,lupfil2+2*i
           ci112=ci112+fil2(l-2*i,2)*x(1,1,l,j)
           ci121=ci121+fil2(l-2*i,1)*x(1,2,l,j)
           ci122=ci122+fil2(l-2*i,2)*x(1,2,l,j)
           ci211=ci211+fil2(l-2*i,1)*x(2,1,l,j)
           ci212=ci212+fil2(l-2*i,2)*x(2,1,l,j)
           ci221=ci221+fil2(l-2*i,1)*x(2,2,l,j)
           ci222=ci222+fil2(l-2*i,2)*x(2,2,l,j)
        enddo
        y(1,j,i)=ci211
        y(2,j,i)=ci121    
        y(3,j,i)=ci221    
        y(4,j,i)=ci112    
        y(5,j,i)=ci212    
        y(6,j,i)=ci122    
        y(7,j,i)=ci222    
     enddo
  enddo

  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)
  ! write(20,*) tel, 1.d-6*nflop/tel

end subroutine comb_rot_shrink_loc_3

subroutine syn_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='syn_repeated_per'
  integer :: nn1,nn2,nn3,i_trans,i_all,i_stat,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  if (num_trans >= 1)  then

     allocate(yy((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,yy,'yy',subname)

     allocate(xx((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,xx,'xx',subname)

  endif

  if (num_trans >= 2) then

     nn1=(nd1+1)/2-1
     nn2=(nd2+1)/2-1
     nn3=(nd3+1)/2-1

     allocate(ww((nn1+1)*(nn2+1)*(nn3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     do i_trans=1,num_trans-1

        n1=2*(n1+1)-1
        n2=2*(n2+1)-1
        n3=2*(n3+1)-1

        if (n1.gt.nd1) stop 'n1 beyond borders'
        if (n2.gt.nd2) stop 'n2 beyond borders'
        if (n3.gt.nd3) stop 'n3 beyond borders'

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 xx(i)=x(i1,i2,i3)
                 i=i+1
              enddo
           enddo
        enddo

        call synthese_per(n1,n2,n3,xx,yy,ww)

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 x(i1,i2,i3)=yy(i)
                 i=i+1
              enddo
           enddo
        enddo

     enddo

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)

  endif

  if (num_trans >= 1) then

     n1=2*(n1+1)-1
     n2=2*(n2+1)-1
     n3=2*(n3+1)-1

     call synthese_per_self(n1,n2,n3,x,xx,yy)

     i_all=-product(shape(xx))*kind(xx)
     deallocate(xx,stat=i_stat)
     call memocc(i_stat,i_all,'xx',subname)

     i_all=-product(shape(yy))*kind(yy)
     deallocate(yy,stat=i_stat)
     call memocc(i_stat,i_all,'yy',subname)

  endif

end subroutine syn_repeated_per



subroutine ana_repeated_per(nd1,nd2,nd3,x,num_trans,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  character(len=*), parameter :: subname='ana_repeated_per'
  integer :: nn1,nn2,nn3,i_trans,i_all,i_stat,i1,i2,i3,i
  real(wp), dimension(:), allocatable :: xx,yy,ww

  n1=nd1
  n2=nd2
  n3=nd3

  if (num_trans.ge.1)  then

     allocate(yy((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,yy,'yy',subname)

     allocate(xx((nd1+1)*(nd2+1)*(nd3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,xx,'xx',subname)

     call analyse_per_self(n1,n2,n3,x,yy,xx)

     n1=(n1+1)/2-1
     n2=(n2+1)/2-1
     n3=(n3+1)/2-1

  endif

  if (num_trans.ge.2) then

     allocate(ww((n1+1)*(n2+1)*(n3+1)+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     do i_trans=2,num_trans

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 xx(i)=x(i1,i2,i3)
                 i=i+1
              enddo
           enddo
        enddo

        call analyse_per(n1,n2,n3,xx,yy,ww)

        i=1
        do i3=0,n3
           do i2=0,n2
              do i1=0,n1
                 x(i1,i2,i3)=yy(i)
                 i=i+1
              enddo
           enddo
        enddo

        n1=(n1+1)/2-1
        n2=(n2+1)/2-1
        n3=(n3+1)/2-1

     enddo

     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww',subname)

  endif

  if (num_trans.ge.1) then 
     i_all=-product(shape(xx))*kind(xx)
     deallocate(xx,stat=i_stat)
     call memocc(i_stat,i_all,'xx',subname)

     i_all=-product(shape(yy))*kind(yy)
     deallocate(yy,stat=i_stat)
     call memocc(i_stat,i_all,'yy',subname)
  endif

end subroutine ana_repeated_per

subroutine synthese_per(nd1,nd2,nd3,x,y,ww)
  ! a periodic synthesis (backward) wavelet transformation
  ! the input array x is not overwritten
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
  call  syn_rot_per(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per(nd3,nt,ww,y)

end subroutine synthese_per

subroutine synthese_per_self(nd1,nd2,nd3,x,y,ww)
  ! a periodic synthesis (backward) wavelet transformation
  ! the input array x is not overwritten
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
  call  syn_rot_per(nd1,nt,x,y)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  syn_rot_per(nd2,nt,y,ww)
  ! i3,i1,i2  -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  syn_rot_per(nd3,nt,ww,x)

end subroutine synthese_per_self


subroutine analyse_per(nd1,nd2,nd3,y,x,ww)
  ! an analysis (forward) periodic wavelet transformation
  ! the input array y is not overwritten
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
  call  ana_rot_per(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per(nd3,nt,ww,x)

end subroutine analyse_per

subroutine analyse_per_self(nd1,nd2,nd3,y,x,ww)
  ! an analysis (forward) periodic wavelet transformation
  ! the input array y is not overwritten
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
  call  ana_rot_per(nd1,nt,y,x)
  ! i2,i3,i1 -> i3,i1,i2
  nt=(nd3+1)*(nd1+1)
  call  ana_rot_per(nd2,nt,x,ww)
  ! i3,i1,i2 -> i1,i2,i3
  nt=(nd1+1)*(nd2+1)
  call  ana_rot_per(nd3,nt,ww,y)

end subroutine analyse_per_self

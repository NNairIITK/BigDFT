subroutine convolut_magic_n_per_cuda(n1,n2,n3,x,y)
  ! Applies the magic filter matrix ( no transposition) ; data set grows
  ! The input array x is not overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_n_per'
  integer :: ndat,i_stat,i_all
  real(kind=4), dimension(:,:,:), allocatable :: wx,wy
  integer, parameter :: lowfil=-8,lupfil=7
  real(kind=4) fil(lowfil:lupfil)

  !CUDA modif
  
  !1. define right data 
 
  DATA fil / &
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

  allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,wx,'wx',subname)

  allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,wy,'wy',subname)


  wx=real(x,kind=4)

  !2.call cuda C interface

  call intertamponcGPU(n1,n2,n3,wx,wy,fil,lowfil,lupfil)

  y=real(wy,kind=wp)

 i_all=-product(shape(wx))*kind(wx)
 deallocate(wx,stat=i_stat)
 call memocc(i_stat,i_all,'wx',subname)

 i_all=-product(shape(wy))*kind(wy)
 deallocate(wy,stat=i_stat)
 call memocc(i_stat,i_all,'wy',subname)

end subroutine convolut_magic_n_per_cuda


!!$
!!$subroutine convolut_magic_n_per(n1,n2,n3,x,y)
!!$  ! Applies the magic filter matrix ( no transposition) ; data set grows
!!$  ! The input array x is not overwritten
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: n1,n2,n3
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
!!$  !local variables
!!$  character(len=*), parameter :: subname='convolut_magic_n_per'
!!$  integer :: ndat,i_stat,i_all
!!$  real(4), dimension(:,:,:), allocatable :: ww,wx,wy
!!$
!!$  allocate(ww(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,ww,'ww',subname)
!!$
!!$  allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,wx,'wx',subname)
!!$
!!$  allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,wy,'wy',subname)
!!$
!!$  wx=real(x,kind=4)
!!$
!!$  !  (i1,i2*i3) -> (i2*i3,I1)
!!$  ndat=(n2+1)*(n3+1)
!!$  call convrot_grow_per_sp(n1,ndat,wx,wy)
!!$  !  (i2,i3*I1) -> (i3*i1,I2)
!!$  ndat=(n3+1)*(n1+1)
!!$  call convrot_grow_per_sp(n2,ndat,wy,ww)
!!$  !  (i3,I1*I2) -> (iI*I2,I3)
!!$  ndat=(n1+1)*(n2+1)
!!$  call convrot_grow_per_sp(n3,ndat,ww,wy)
!!$
!!$  y=real(wy,kind=wp)
!!$
!!$  i_all=-product(shape(ww))*kind(ww)
!!$  deallocate(ww,stat=i_stat)
!!$  call memocc(i_stat,i_all,'ww',subname)
!!$
!!$  i_all=-product(shape(wx))*kind(wx)
!!$  deallocate(wx,stat=i_stat)
!!$  call memocc(i_stat,i_all,'wx',subname)
!!$
!!$  i_all=-product(shape(wy))*kind(wy)
!!$  deallocate(wy,stat=i_stat)
!!$  call memocc(i_stat,i_all,'wy',subname)
!!$
!!$end subroutine convolut_magic_n_per
!!$
!!$
!!$subroutine convolut_magic_t_per(n1,n2,n3,x,y)
!!$  ! Applies the magic filter matrix transposed ; data set shrinks
!!$  ! The input array x is overwritten
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: n1,n2,n3
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: x
!!$  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
!!$  !local variables
!!$  character(len=*), parameter :: subname='convolut_magic_t_per'
!!$  integer :: ndat,i_stat,i_all
!!$  real(4), dimension(:,:,:), allocatable :: ww,wx,wy
!!$
!!$  allocate(ww(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,ww,'ww',subname)
!!$
!!$  allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,wx,'wx',subname)
!!$
!!$  allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!!$  call memocc(i_stat,wy,'wy',subname)
!!$
!!$  wx=real(x,kind=4)
!!$  
!!$  !  (I1,I2*I3) -> (I2*I3,i1)
!!$  ndat=(n2+1)*(n3+1)
!!$  call convrot_shrink_per_sp(n1,ndat,wx,ww)
!!$  !  (I2,I3*i1) -> (I3*i1,i2)
!!$  ndat=(n3+1)*(n1+1)
!!$  call convrot_shrink_per_sp(n2,ndat,ww,wx)
!!$  !  (I3,i1*i2) -> (i1*i2,i3)
!!$  ndat=(n1+1)*(n2+1)
!!$  call convrot_shrink_per_sp(n3,ndat,wx,wy)
!!$  
!!$  y=real(wy,kind=wp)
!!$
!!$  i_all=-product(shape(ww))*kind(ww)
!!$  deallocate(ww,stat=i_stat)
!!$  call memocc(i_stat,i_all,'ww',subname)
!!$
!!$  i_all=-product(shape(wx))*kind(wx)
!!$  deallocate(wx,stat=i_stat)
!!$  call memocc(i_stat,i_all,'wx',subname)
!!$
!!$  i_all=-product(shape(wy))*kind(wy)
!!$  deallocate(wy,stat=i_stat)
!!$  call memocc(i_stat,i_all,'wy',subname)
!!$
!!$end subroutine convolut_magic_t_per

subroutine convolut_magic_t_per_cuda(n1,n2,n3,x,y)
  ! Applies the magic filter matrix transposed ; data set shrinks
  ! The input array x is overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-7,lupfil=8
  character(len=*), parameter :: subname='convolut_magic_t_per'
  integer :: ndat,i_stat,i_all
  real(kind=4), dimension(:,:,:), allocatable :: wx,wy
  real(kind=4) fil(lowfil:lupfil)
  DATA fil / &
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

  allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,wx,'wx',subname)

  allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,wy,'wy',subname)


  wx=real(x,kind=4)

  !2.call cuda C interface

  call intertamponcGPU(n1,n2,n3,wx,wy,fil,lowfil,lupfil)

  y=real(wy,kind=wp)
    
 i_all=-product(shape(wx))*kind(wx)
 deallocate(wx,stat=i_stat)
 call memocc(i_stat,i_all,'wx',subname)

 i_all=-product(shape(wy))*kind(wy)
 deallocate(wy,stat=i_stat)
 call memocc(i_stat,i_all,'wy',subname)

end subroutine convolut_magic_t_per_cuda




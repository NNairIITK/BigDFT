
subroutine analyse_per(n1,n2,n3,ww,y,x)
  ! Analysis wavelet transformation in periodic BC
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
  call  ana_rot_per(n1,nt,y,x)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_per(n2,nt,x,ww)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,ww,x)

end subroutine analyse_per

subroutine analyse_per_self(n1,n2,n3,y,x)
  ! Analysis wavelet transformation  in periodic BC
  ! The input array y is overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(inout) :: y
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: x
  !local variables
  integer :: nt

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+2)*(2*n3+2)
  call  ana_rot_per(n1,nt,y,x)
  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(2*n1+2)
  call  ana_rot_per(n2,nt,x,y)
  ! I3,i1,i2 -> i1,i2,i3
  nt=(2*n1+2)*(2*n2+2)
  call  ana_rot_per(n3,nt,y,x)

end subroutine analyse_per_self


subroutine synthese_per(n1,n2,n3,ww,x,y)
  ! A synthesis wavelet transformation  in periodic BC
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
  call  syn_rot_per(n1,nt,x,y) 
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+2)
  call  syn_rot_per(n2,nt,y,ww)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call  syn_rot_per(n3,nt,ww,y)

end subroutine synthese_per


subroutine synthese_per_self(n1,n2,n3,x,y)
  ! A synthesis wavelet transformation in periodic BC
  ! The input array x is  overwritten
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: x
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(inout) :: y
  !local variables
  integer :: nt

  ! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+2)*(2*n3+2)
  call  syn_rot_per(n1,nt,x,y) 
  ! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+2)*(2*n1+2)
  call  syn_rot_per(n2,nt,y,x)
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call  syn_rot_per(n3,nt,x,y)

end subroutine synthese_per_self


! Applies the magic filter matrix in periodic BC ( no transposition)
! The input array x is not overwritten
! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_per(n1,n2,n3,x,y,ww)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_n_per'
  integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat,i_stat,i_all
  real(wp), dimension(0:n1,0:n2,0:n3):: ww ! work array
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


  if (.not. GPUconv) then !traditional CPU computation

     !  (i1,i2*i3) -> (i2*i3,I1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*I1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_n_per(n2,ndat,y,ww)
     !  (i3,I1*I2) -> (iI*I2,I3)
     ndat=(n1+1)*(n2+1)
     call convrot_n_per(n3,ndat,ww,y)

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


! Applies the magic filter matrix in periodic BC ( no transposition)
! The input array x is overwritten in the usual case
! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_n_per_self(n1,n2,n3,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_n_per'
  integer, parameter :: lowfil=-8,lupfil=7 !for GPU computation
  integer :: ndat,i_stat,i_all
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


  if (.not. GPUconv) then !traditional CPU computation

     !  (i1,i2*i3) -> (i2*i3,I1)
     ndat=(n2+1)*(n3+1)
     call convrot_n_per(n1,ndat,x,y)
     !  (i2,i3*I1) -> (i3*i1,I2)
     ndat=(n3+1)*(n1+1)
     call convrot_n_per(n2,ndat,y,x)
     !  (i3,I1*I2) -> (iI*I2,I3)
     ndat=(n1+1)*(n2+1)
     call convrot_n_per(n3,ndat,x,y)

  else
     allocate(wx(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wx,'wx',subname)

     allocate(wy(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,wy,'wy',subname)

     !copy input data
     !here we should put a loop to avoid the ndebug issue
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
end subroutine convolut_magic_n_per_self


! Applies the magic filter matrix transposed in periodic BC 
! The input array x is overwritten
! this routine is modified to accept the GPU convolution if it is the case
subroutine convolut_magic_t_per_self(n1,n2,n3,x,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y
  !local variables
  character(len=*), parameter :: subname='convolut_magic_t_per'
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: ndat,i_stat,i_all
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

  
  if (.not. GPUconv) then

     !  (I1,I2*I3) -> (I2*I3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_t_per(n1,ndat,x,y)
     !  (I2,I3*i1) -> (I3*i1,i2)
     ndat=(n3+1)*(n1+1)
     call convrot_t_per(n2,ndat,y,x)
     !  (I3,i1*i2) -> (i1*i2,i3)
     ndat=(n1+1)*(n2+1)
     call convrot_t_per(n3,ndat,x,y)

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

end subroutine convolut_magic_t_per_self





! Applies the magic filter matrix transposed  in periodic BC
! The input array x is not overwritten
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

  
  if (.not. GPUconv) then
     allocate(ww(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,ww,'ww',subname)

     !  (I1,I2*I3) -> (I2*I3,i1)
     ndat=(n2+1)*(n3+1)
     call convrot_t_per(n1,ndat,x,ww)
     !  (I2,I3*i1) -> (I3*i1,i2)
     ndat=(n3+1)*(n1+1)
     call convrot_t_per(n2,ndat,ww,x)
     !  (I3,i1*i2) -> (i1*i2,i3)
     ndat=(n1+1)*(n2+1)
     call convrot_t_per(n3,ndat,x,y)

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


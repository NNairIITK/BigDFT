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

!prepare the keys to be passed in GPU global or constant memory (presumably texture)
subroutine adjust_keys_for_gpu(nseg_c,nseg_f,keyv_c,keyg_c,keyv_f,keyg_f,nvctr_c,keys_GPU)
  use module_base
  implicit none
  integer, intent(in) :: nseg_c,nseg_f,nvctr_c
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f  
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(kind=8), intent(out) :: keys_GPU
  !local variables
  character(len=*), parameter :: subname='adjust_keys_for_gpu'
  !maximum number of elements to be passed into a block
  integer :: iseg,i_stat,i_all,segment_elements,nseg_tot,j,nseggpu,nblocks,nseg_block
  integer :: elems_block,nblocks_max,ncuts
  integer, dimension(:,:), allocatable :: keys
  
  call readkeysinfo(elems_block,nblocks_max)
  
  nseg_tot=0
  !assign the keys values, coarse part
  do iseg=1,nseg_c
     segment_elements=keyg_c(2,iseg)-keyg_c(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts
           nseg_tot=nseg_tot+1
        end do
     end if
  end do
  !assign the keys values, fine part
  do iseg=1,nseg_f
     segment_elements=keyg_f(2,iseg)-keyg_f(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts
           nseg_tot=nseg_tot+1
        end do
     end if
  end do

  !number of segments treated by each block;
  nseg_block =(nseg_tot+1)/nblocks_max +1

  nblocks=min(nblocks_max,nseg_tot+1)

  !print *,'nseg,nblocks',nseg_block,nblocks
  call copynseginfo(nblocks,nseg_block)

  nseggpu=nseg_block*nblocks

  !allocate the array to be copied
  allocate(keys(4,nseggpu+ndebug),stat=i_stat)
  call memocc(i_stat,keys,'keys',subname)
  
  nseg_tot=0
  !assign the keys values, coarse part
  do iseg=1,nseg_c
     segment_elements=keyg_c(2,iseg)-keyg_c(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1;
        keys(1,nseg_tot)=1 !which means coarse segment
        keys(2,nseg_tot)=segment_elements
        keys(3,nseg_tot)=keyv_c(iseg)
        keys(4,nseg_tot)=keyg_c(1,iseg)
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts-1
           nseg_tot=nseg_tot+1
           keys(1,nseg_tot)=1
           keys(2,nseg_tot)=elems_block
           keys(3,nseg_tot)=keyv_c(iseg)+(j-1)*elems_block
           keys(4,nseg_tot)=keyg_c(1,iseg)+(j-1)*elems_block
        end do
        j=ncuts
        nseg_tot=nseg_tot+1
        keys(1,nseg_tot)=1
        keys(2,nseg_tot)=segment_elements-(j-1)*elems_block
        keys(3,nseg_tot)=keyv_c(iseg)+(j-1)*elems_block
        keys(4,nseg_tot)=keyg_c(1,iseg)+(j-1)*elems_block
     end if
  end do
  !assign the keys values, fine part
  do iseg=1,nseg_f
     segment_elements=keyg_f(2,iseg)-keyg_f(1,iseg)+1
     if (segment_elements <= elems_block) then
        nseg_tot=nseg_tot+1;
        keys(1,nseg_tot)=nvctr_c !which means fine segment
        keys(2,nseg_tot)=segment_elements
        keys(3,nseg_tot)=keyv_f(iseg)
        keys(4,nseg_tot)=keyg_f(1,iseg)
     else
        ncuts=(segment_elements-1)/elems_block+1
        do j=1,ncuts-1
           nseg_tot=nseg_tot+1
           keys(1,nseg_tot)=nvctr_c
           keys(2,nseg_tot)=elems_block
           keys(3,nseg_tot)=keyv_f(iseg)+(j-1)*elems_block
           keys(4,nseg_tot)=keyg_f(1,iseg)+(j-1)*elems_block
        end do
        j=ncuts
        nseg_tot=nseg_tot+1
        keys(1,nseg_tot)=nvctr_c
        keys(2,nseg_tot)=segment_elements-(j-1)*elems_block
        keys(3,nseg_tot)=keyv_f(iseg)+(j-1)*elems_block
        keys(4,nseg_tot)=keyg_f(1,iseg)+(j-1)*elems_block
     end if
  end do

  !fill the last part of the segments
  do iseg=nseg_tot+1,nblocks*nseg_block
     keys(1,iseg)=0
     keys(2,iseg)=0
     keys(3,iseg)=0
     keys(4,iseg)=0
  end do

  !allocate the gpu pointer and copy the values
  call GPU_int_allocate(4*nseggpu,keys_GPU,i_stat)

  call GPU_int_send(4*nseggpu,keys,keys_GPU,i_stat)


  i_all=-product(shape(keys))
  deallocate(keys,stat=i_stat)
  call memocc(i_stat,i_all,'keys',subname)


end subroutine adjust_keys_for_gpu

subroutine prepare_gpu_for_locham(n1,n2,n3,wfd,orbs,GPU)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: i_stat

  !after this call, all memory operations are in double precision, 
  !call set_gpu_simple() in order to have simple memory operations
  call set_gpu_double() 

  call adjust_keys_for_gpu(wfd%nseg_c,wfd%nseg_f,wfd%keyv(1),wfd%keyg(1,1),&
       wfd%keyv(wfd%nseg_c+1),wfd%keyg(1,wfd%nseg_c+1),wfd%nvctr_c,GPU%keys)

  !create parameters
  call creategpuparameters(n1,n2,n3)

  !allocate space on the card
  call GPU_allocate((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp,GPU%psi,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work1,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work2,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work3,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%work4,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),GPU%pot,i_stat)
  call GPU_allocate((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp,GPU%hpsi,i_stat)
  
end subroutine prepare_gpu_for_locham



subroutine free_gpu(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: i_stat
  

  call GPU_deallocate(GPU%work1,i_stat)
  call GPU_deallocate(GPU%work2,i_stat)
  call GPU_deallocate(GPU%work3,i_stat)
  call GPU_deallocate(GPU%work4,i_stat)
  call GPU_deallocate(GPU%keys,i_stat)
  call GPU_deallocate(GPU%pot,i_stat)
  call GPU_deallocate(GPU%psi,i_stat)
  call GPU_deallocate(GPU%hpsi,i_stat)

end subroutine free_gpu

subroutine gpu_locham(n1,n2,n3,hx,hy,hz,orbs,GPU,ekin_sum,epot_sum)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  real(gp), intent(out) :: ekin_sum,epot_sum
  type(orbitals_data), intent(in) :: orbs
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: i_stat,iorb
  real(gp) :: ekin,epot

  do iorb=1,orbs%norbp

     call gpufulllocham(n1,n2,n3,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
          GPU%psi,GPU%hpsi,GPU%pot,GPU%keys,&
          GPU%work4,GPU%work3,GPU%work1,GPU%work2,epot,ekin)

     ekin_sum=ekin_sum+orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%occup(iorb+orbs%isorb)*epot

  end do

end subroutine gpu_locham


